#include "saft.hpp"
#include <algorithm>
#include <bitset>
#include <charconv>
#include <cmath>
#include <cstring>
#include <fstream>
#include <mutex>
#include <numeric>
#include <thread>

namespace saft {
const QStringList GroupNames = {"left", "weld", "right"};
const QStringList MuteNames = {"Rayleigh", "L surface", "S surface",
                               "ID LL",    "ID LS/SL",  "ID SS"};
const QStringList PulseNames = {
    "Rayleigh",    "L surface",
    "S surface",   "Backwall LL",
    "Backwall LS", "Backwall SS",
    "Backwall SL", "Repeated SSSS (backwall-top-backwall)"};
const QStringList ReflectedScreenNames = {"Backwall LL", "Backwall LS/SL",
                                          "Backwall SS", "Repeated SSSS"};
const QStringList &shortlistPaths() {
  static const QStringList ids{"1-1_LL",   "1-1_LS",    "1-1_SL",   "1-1_SS",
                               "1-2_SSS",  "2-1_SSS",   "1-3_SSSS", "2-2_SSSS",
                               "3-1_SSSS", "2-3_SSSSS", "3-2_SSSSS"};
  return ids;
}
static void require(bool ok, const char *why) {
  if (!ok)
    throw std::runtime_error(why);
}
static void require(bool ok, const QString &why) {
  if (!ok)
    throw std::runtime_error(why.toStdString());
}
QString sha(const QByteArray &b) {
  return QCryptographicHash::hash(b, QCryptographicHash::Sha256).toHex();
}
QString shaFile(const QString &path) {
  QFile f(path);
  require(f.open(QIODevice::ReadOnly), "Cannot read " + path);
  QCryptographicHash h(QCryptographicHash::Sha256);
  h.addData(&f);
  return h.result().toHex();
}
QJsonObject readJson(const QString &path) {
  QFile f(path);
  require(f.open(QIODevice::ReadOnly), "Cannot read " + path);
  QJsonParseError e;
  auto d = QJsonDocument::fromJson(f.readAll(), &e);
  require(e.error == QJsonParseError::NoError && d.isObject(),
          "Invalid JSON: " + path);
  return d.object();
}
void writeJson(const QString &path, const QJsonObject &o) {
  QSaveFile f(path);
  require(f.open(QIODevice::WriteOnly), "Cannot write " + path);
  f.write(QJsonDocument(o).toJson());
  require(f.commit(), "Cannot save " + path);
}
void parallel(int count, int threads, const std::function<void(int)> &fn,
              std::atomic_bool *cancel) {
  std::atomic_int next{0};
  std::atomic_bool failed{false};
  std::exception_ptr failure;
  std::mutex lock;
  auto work = [&] {
    try {
      for (;;) {
        if (failed || (cancel && cancel->load()))
          break;
        int i = next.fetch_add(1);
        if (i >= count)
          break;
        fn(i);
      }
    } catch (...) {
      std::lock_guard<std::mutex> g(lock);
      if (!failure)
        failure = std::current_exception();
      failed = true;
    }
  };
  std::vector<std::thread> pool;
  int n = std::clamp(threads, 1, std::max(1, count));
  for (int i = 1; i < n; ++i)
    pool.emplace_back(work);
  work();
  for (auto &t : pool)
    t.join();
  if (failure)
    std::rethrow_exception(failure);
  if (cancel && cancel->load())
    throw Cancelled();
}
Config::Config(int count) {
  receivers.resize(count);
  std::iota(receivers.begin(), receivers.end(), 0);
}
static Vec axis(double a, double b, double step) {
  Vec v;
  int n = int(std::floor((b - a) / step + 1e-6)) + 1;
  if (n < 1 || n > 2000000)
    return v;
  for (int i = 0; i < n; ++i)
    v.push_back(a + i * step);
  return v;
}
Vec Config::xs() const { return axis(x0, x1, spacing); }
Vec Config::zs() const { return axis(z0, z1, spacing); }
void Config::validate(const Material &m) const {
  for (double v :
       {low, high, x0, x1, z0, z1, spacing, gate, wall, offset, muteWidth,
        cfPower, fraction, minimumEffective, screenTailDb, screenMargin})
    require(std::isfinite(v), "Numeric settings must be finite");
  require(low > 0 && high > low && high < 50,
          "Use 0 < low frequency < high frequency < 50 MHz");
  require(x0 >= -m.active * 1000 && x1 <= m.active * 1000 && x1 > x0 &&
              z0 > 0 && z1 < m.h * 1000 && z1 > z0,
          "Image bounds must lie inside the active plate and have increasing "
          "limits");
  require(spacing >= .01 && spacing <= 5, "Pixel spacing must be 0.01 to 5 mm");
  require(!xs().empty() && !zs().empty() && xs().size() * zs().size() <= 250000,
          "Grid exceeds 250,000 pixels; increase spacing or reduce extent");
  require(!sources.empty() && sources.size() <= 5,
          "Select at least one source");
  std::set<int> ss(sources.begin(), sources.end());
  require(ss.size() == sources.size() && *ss.begin() >= 0 && *ss.rbegin() < 5,
          "Invalid or duplicate source ID");
  require(!receivers.empty() &&
              std::all_of(receivers.begin(), receivers.end(),
                          [&](int i) { return i >= 0 && i < m.receiverCount; }),
          QString("Receiver IDs must be R0 to R%1").arg(m.receiverCount - 1));
  require(std::any_of(groups.begin(), groups.end(), [](bool b) { return b; }),
          "Select at least one receiver group");
  for (const auto &[source, ids] : sourceReceivers) {
    require(source >= 0 && source < 5, "Invalid per-source receiver key");
    require(
        std::set<int>(ids.begin(), ids.end()).size() == ids.size() &&
            std::all_of(ids.begin(), ids.end(),
                        [&](int r) { return r >= 0 && r < m.receiverCount; }),
        "Invalid or duplicate per-source receiver ID");
  }
  require(gate >= 0 && gate <= 1 && wall >= 0 && wall <= 1 && offset >= 0 &&
              offset <= 120,
          "Invalid response gate or source exclusion");
  require(minimumSources >= 1 && minimumSources <= 5 && minimumReceivers >= 1 &&
              minimumReceivers <= m.receiverCount && minimumEffective >= 0 &&
              minimumEffective <= m.receiverCount && fraction >= 0 &&
              fraction <= 1,
          "Invalid coverage limits");
  require(muteWidth >= .01 && muteWidth <= 10 && cfPower >= 0 && cfPower <= 8 &&
              threads >= 1 && threads <= 64,
          "Invalid suppression or worker setting");
  require(screenTailDb >= 6 && screenTailDb <= 60 && screenMargin >= 0 &&
              screenMargin <= 3,
          "Invalid pulse-screen tail threshold or timing margin");
  require(
      QStringList{"vy", "vx"}.contains(component) &&
          QStringList{"all", "left", "right", "even", "odd"}.contains(
              aperture) &&
          QStringList{"rms", "coherent"}.contains(combination) &&
          QStringList{"minimum", "fixed"}.contains(policy) &&
          QStringList{"weighted", "raw"}.contains(processing) &&
          QStringList{"case", "paired"}.contains(repairs) &&
          QStringList{"off", "L", "S", "R", "LS", "LR", "SR", "LSR"}.contains(
              screening),
      "Unknown reconstruction option");
  require(!cases.empty(), "Choose an image signal");
  for (const auto &c : cases)
    require(QStringList{"lof", "reference", "difference"}.contains(c),
            "Unknown image signal");
  if (freeze)
    require(std::isfinite(anchor.x) && std::isfinite(anchor.z) &&
                std::abs(anchor.x) <= m.active * 1000 && anchor.z > 0 &&
                anchor.z < m.h * 1000,
            "Set an inspection point inside the plate before freezing its "
            "aperture");
}
QJsonObject Config::json() const {
  QJsonObject o{{"path", path},         {"component", component},
                {"aperture", aperture}, {"combination", combination},
                {"policy", policy},     {"processing", processing},
                {"repairs", repairs},   {"pulse_screening", screening}};
  for (auto p : std::initializer_list<std::pair<QString, double>>{
           {"low_MHz", low},
           {"high_MHz", high},
           {"x_min_mm", x0},
           {"x_max_mm", x1},
           {"depth_min_mm", z0},
           {"depth_max_mm", z1},
           {"spacing_mm", spacing},
           {"angular_gate", gate},
           {"wall_gate", wall},
           {"min_offset_mm", offset},
           {"mute_halfwidth_us", muteWidth},
           {"screen_tail_db", screenTailDb},
           {"screen_margin_us", screenMargin},
           {"coherence_power", cfPower},
           {"receiver_fraction", fraction},
           {"minimum_effective_receivers", minimumEffective},
           {"anchor_x_mm", anchor.x},
           {"anchor_depth_mm", anchor.z}})
    o[p.first] = p.second;
  o["minimum_sources"] = minimumSources;
  o["minimum_receivers"] = minimumReceivers;
  o["threads"] = threads;
  o["clean"] = clean;
  o["coherence"] = coherence;
  o["freeze"] = freeze;
  QJsonArray a;
  for (int s : sources)
    a.append(s);
  o["sources"] = a;
  a = {};
  for (int r : receivers)
    a.append(r);
  o["receivers"] = a;
  if (!sourceReceivers.empty()) {
    QJsonObject overrides;
    for (const auto &[source, ids] : sourceReceivers) {
      QJsonArray values;
      for (int r : ids)
        values.append(r);
      overrides[QString::number(source)] = values;
    }
    o["source_receiver_ids"] = overrides;
  }
  a = {};
  for (bool g : groups)
    a.append(g);
  o["groups"] = a;
  a = {};
  for (bool g : mutes)
    a.append(g);
  o["mutes"] = a;
  a = {};
  for (bool on : reflectedScreen)
    a.append(on);
  o["pulse_reflections"] = a;
  o["cases"] = QJsonArray::fromStringList(cases);
  return o;
}
Config Config::fromJson(QJsonObject o, int receiverCount) {
  Config c(receiverCount);
  for (auto p : std::initializer_list<std::pair<QString, QString *>>{
           {"path", &c.path},
           {"component", &c.component},
           {"aperture", &c.aperture},
           {"combination", &c.combination},
           {"policy", &c.policy},
           {"processing", &c.processing},
           {"repairs", &c.repairs},
           {"pulse_screening", &c.screening}})
    if (o.contains(p.first))
      *p.second = o[p.first].toString();
  for (auto p : std::initializer_list<std::pair<QString, double *>>{
           {"low_MHz", &c.low},
           {"high_MHz", &c.high},
           {"x_min_mm", &c.x0},
           {"x_max_mm", &c.x1},
           {"depth_min_mm", &c.z0},
           {"depth_max_mm", &c.z1},
           {"spacing_mm", &c.spacing},
           {"angular_gate", &c.gate},
           {"wall_gate", &c.wall},
           {"min_offset_mm", &c.offset},
           {"mute_halfwidth_us", &c.muteWidth},
           {"screen_tail_db", &c.screenTailDb},
           {"screen_margin_us", &c.screenMargin},
           {"coherence_power", &c.cfPower},
           {"receiver_fraction", &c.fraction},
           {"minimum_effective_receivers", &c.minimumEffective},
           {"anchor_x_mm", &c.anchor.x},
           {"anchor_depth_mm", &c.anchor.z}})
    if (o.contains(p.first))
      *p.second = o[p.first].toDouble(std::nan(""));
  for (auto p : std::initializer_list<std::pair<QString, int *>>{
           {"minimum_sources", &c.minimumSources},
           {"minimum_receivers", &c.minimumReceivers},
           {"threads", &c.threads}})
    if (o.contains(p.first))
      *p.second = o[p.first].toInt(-1);
  c.clean = o["clean"].toBool();
  if (o.contains("source_receiver_ids")) {
    require(o["source_receiver_ids"].isObject(),
            "Expected per-source receiver lists");
    auto values = o["source_receiver_ids"].toObject();
    for (auto it = values.begin(); it != values.end(); ++it) {
      bool ok = false;
      int s = it.key().toInt(&ok);
      require(ok && s >= 0 && s < 5 && it.value().isArray(),
              "Invalid per-source receiver entry");
      auto &ids = c.sourceReceivers[s];
      for (auto id : it.value().toArray()) {
        require(id.isDouble() && id.toDouble() == id.toInt(-1),
                "Receiver IDs must be integers");
        ids.push_back(id.toInt(-1));
      }
    }
  }
  c.coherence = o["coherence"].toBool();
  c.freeze = o["freeze"].toBool();
  if (o.contains("cases")) {
    c.cases.clear();
    for (auto v : o["cases"].toArray())
      c.cases.append(v.toString());
  }
  for (auto p : std::initializer_list<std::pair<QString, std::vector<int> *>>{
           {"sources", &c.sources}, {"receivers", &c.receivers}})
    if (o.contains(p.first)) {
      p.second->clear();
      for (auto v : o[p.first].toArray())
        p.second->push_back(v.toInt(-1));
    }
  if (o.contains("groups")) {
    auto a = o["groups"].toArray();
    require(a.size() == 3, "Expected three receiver groups");
    for (int i = 0; i < 3; ++i)
      c.groups[i] = a[i].toBool();
  }
  if (o.contains("mutes")) {
    auto a = o["mutes"].toArray();
    require(a.size() == 6, "Expected six arrival families");
    for (int i = 0; i < 6; ++i)
      c.mutes[i] = a[i].toBool();
  }
  if (o.contains("pulse_reflections")) {
    auto a = o["pulse_reflections"].toArray();
    require(a.size() == 4,
            "Expected four reflected pulse selections: LL, LS/SL, SS, SSSS");
    for (int i = 0; i < 4; ++i) {
      require(a[i].isBool(),
              "Reflected pulse selections must be true or false");
      c.reflectedScreen[i] = a[i].toBool();
    }
  }
  return c;
}
std::vector<int> receiverIds(const QString &text, int count) {
  if (text.trimmed().isEmpty())
    return Config(count).receivers;
  std::set<int> found;
  QRegularExpression re("^\\s*[Rr]?(\\d+)\\s*(?:-\\s*[Rr]?(\\d+))?\\s*$");
  for (auto part : text.split(',')) {
    auto m = re.match(part);
    require(m.hasMatch(), "Receiver IDs: use 0-156, R200-R240");
    int a = m.captured(1).toInt(),
        b = m.captured(2).isEmpty() ? a : m.captured(2).toInt();
    require(a >= 0 && a <= b && b < count,
            QString("Receiver IDs must lie between R0 and R%1").arg(count - 1));
    for (int i = a; i <= b; ++i)
      found.insert(i);
  }
  return {found.begin(), found.end()};
}

static std::array<Complex, 2> polar(double p, Complex q, char mode,
                                    const Material &m) {
  return mode == 'L' ? std::array<Complex, 2>{m.cl * p, m.cl * q}
                     : std::array<Complex, 2>{-m.cs * q, m.cs * p};
}
static std::array<Complex, 2> traction(double p, Complex q, char mode,
                                       const Material &m) {
  auto e = polar(p, q, mode, m);
  double mu = m.cs * m.cs, lam = m.cl * m.cl - 2 * mu;
  return {mu * (q * e[0] + p * e[1]),
          lam * (p * e[0] + q * e[1]) + 2 * mu * q * e[1]};
}
static std::array<Complex, 2> solve(std::array<Complex, 2> a,
                                    std::array<Complex, 2> b,
                                    std::array<Complex, 2> rhs) {
  auto det = a[0] * b[1] - b[0] * a[1];
  return {(rhs[0] * b[1] - b[0] * rhs[1]) / det,
          (a[0] * rhs[1] - rhs[0] * a[1]) / det};
}
Reflection reflect(double angle, char mode, const Material &m) {
  double p = std::sin(angle) / (mode == 'L' ? m.cl : m.cs);
  Complex ql = std::conj(std::sqrt(Complex(1 / (m.cl * m.cl) - p * p, 0))),
          qs = std::conj(std::sqrt(Complex(1 / (m.cs * m.cs) - p * p, 0)));
  auto a = traction(p, ql, 'L', m), b = traction(p, qs, 'S', m),
       inc = traction(p, mode == 'L' ? -ql : -qs, mode, m);
  auto r = solve(a, b, {-inc[0], -inc[1]});
  auto surf = polar(p, mode == 'L' ? -ql : -qs, mode, m),
       el = polar(p, ql, 'L', m), es = polar(p, qs, 'S', m);
  for (int i = 0; i < 2; ++i)
    surf[i] += r[0] * el[i] + r[1] * es[i];
  return {r, surf};
}
Ray ray(Point a, Point b, const QString &modes, const Material &m,
        bool outgoing, bool points, double top) {
  Ray r;
  std::array<double, 5> heights{};
  heights[0] = a.z;
  int n = modes.size();
  require(n >= 1 && n <= 4, "Invalid ray branch");
  if (n == 2)
    heights[1] = m.h;
  if (n == 3) {
    heights[1] = outgoing ? top : m.h;
    heights[2] = outgoing ? m.h : top;
  }
  if (n == 4) {
    heights[1] = heights[3] = m.h;
    heights[2] = top;
  }
  heights[n] = b.z;
  double hs[4], cs[4], total = 0, maxc = 0;
  bool same = true;
  for (int i = 0; i < n; ++i) {
    hs[i] = std::abs(heights[i + 1] - heights[i]);
    require(hs[i] > 0, "Nonpositive ray leg");
    cs[i] = modes[i] == 'L' ? m.cl : m.cs;
    total += hs[i];
    maxc = std::max(maxc, cs[i]);
    if (i && modes[i] != modes[0])
      same = false;
  }
  double dx = b.x - a.x, p = 0;
  if (same)
    p = dx / (std::hypot(dx, total) * cs[0]);
  else {
    double lo = 0, hi = (1 - 1e-13) / maxc;
    for (int i = 0; i < 44; ++i) {
      double mid = (lo + hi) / 2, advance = 0;
      for (int j = 0; j < n; ++j) {
        double pc = mid * cs[j];
        advance += hs[j] * pc / std::sqrt(1 - pc * pc);
      }
      if (advance < std::abs(dx))
        lo = mid;
      else
        hi = mid;
    }
    p = (dx > 0 ? 1 : dx < 0 ? -1 : 0) * (lo + hi) / 2;
  }
  if (points)
    r.points.push_back(a);
  double offset = 0;
  for (int i = 0; i < n; ++i) {
    double pc = p * cs[i], cosine = std::sqrt(std::max(0., 1 - pc * pc));
    offset += hs[i] * pc / cosine;
    r.time += hs[i] / (cs[i] * cosine);
    if (i < n - 1) {
      auto f =
          reflect(std::asin(std::clamp(pc, -1., 1.)), modes[i].toLatin1(), m);
      double ratio = std::sqrt(std::max(0., 1 - std::pow(p * cs[i + 1], 2))) /
                     std::max(cosine, 1e-30);
      Complex q = f.amplitudes[modes[i + 1] == 'L' ? 0 : 1] *
                  std::sqrt(cs[i + 1] / cs[i] * ratio);
      if (heights[i + 1] == m.h && modes[i] != modes[i + 1])
        q = -q;
      r.wall *= q;
    }
    if (points)
      r.points.push_back(i == n - 1 ? b : Point{a.x + offset, heights[i + 1]});
  }
  r.first = std::asin(std::clamp(p * cs[0], -1., 1.));
  r.last = std::asin(std::clamp(p * cs[n - 1], -1., 1.));
  return r;
}
void fft(CVec &a, bool inverse) {
  int n = int(a.size());
  for (int i = 1, j = 0; i < n; ++i) {
    int bit = n >> 1;
    for (; j & bit; bit >>= 1)
      j ^= bit;
    j ^= bit;
    if (i < j)
      std::swap(a[i], a[j]);
  }
  for (int len = 2; len <= n; len <<= 1) {
    Complex wl = std::polar(1., (inverse ? 2 : -2) * Pi / len);
    for (int i = 0; i < n; i += len) {
      Complex w = 1;
      for (int j = 0; j < len / 2; ++j) {
        auto u = a[i + j], v = a[i + j + len / 2] * w;
        a[i + j] = u + v;
        a[i + j + len / 2] = u - v;
        w *= wl;
      }
    }
  }
  if (inverse)
    for (auto &z : a)
      z /= n;
}
CVec kernel(double low, double high) {
  CVec h(2 * Half + 1);
  double width = (high - low) * 1e6, centre = (high + low) * .5e6,
         beta = .1102 * (55 - 8.7), normal = 0;
  for (int n = -Half; n <= Half; ++n) {
    double t = width * n * Dt, sinc = n ? std::sin(Pi * t) / (Pi * t) : 1,
           window = std::cyl_bessel_i(
                        0, beta * std::sqrt(std::max(
                                      0., 1 - std::pow(double(n) / Half, 2)))) /
                    std::cyl_bessel_i(0, beta);
    double value = 2 * width * Dt * sinc * window;
    h[n + Half] = value * std::polar(1., 2 * Pi * centre * n * Dt);
    normal += value;
  }
  for (auto &z : h)
    z *= 2 / normal;
  return h;
}

Engine::Engine(QString directory, QString dataset)
    : root(QDir(directory).absolutePath()) {
  cache = root + "/build/saft_cache";
  require(!dataset.isEmpty(),
          "Choose an external SAFT dataset folder containing dataset.json");
  {
    data = QDir(dataset).absolutePath();
    const auto description = readJson(data + "/dataset.json");
    require(description["schema"].toInt() == 1 &&
                description["kind"] == "SAFT acquisition dataset",
            "Choose a SAFT dataset folder containing dataset.json");
    bundled = true;
    manifestHash = shaFile(data + "/dataset.json");
    datasetName = description["name"].toString();
    datasetKind = description["model_kind"].toString();
    require(QStringList{"full", "partial"}.contains(datasetKind),
            "Unsupported specimen geometry");
    audit = description["audit"].toObject();
    arrivalAudit = {
        {"baseline_input_fingerprint", description["input_fingerprint"]},
        {"material", description["material"]},
        {"model_inputs", QJsonArray{}},
        {"model_evidence", description["model_evidence"]}};
    geometry = description["geometry"].toObject();
    for (const auto value :
         description["surface_profile_x_depth_m"].toArray()) {
      const auto xy = value.toArray();
      require(xy.size() == 2, "Invalid exposed-surface coordinates");
      surfaceProfile.push_back({xy[0].toDouble(), xy[1].toDouble()});
    }
    require(surfaceProfile.size() >= 2 &&
                std::is_sorted(surfaceProfile.begin(), surfaceProfile.end(),
                               [](Point a, Point b) { return a.x < b.x; }),
            "Invalid surface profile");
  }
  fingerprint = audit["input_fingerprint"].toString();
  if (bundled)
    cache += "/" + fingerprint;
  require(QDir().mkpath(cache), "Cannot create dataset cache");
  require(audit["status"] == "passed" &&
              arrivalAudit["baseline_input_fingerprint"] == fingerprint,
          "Saved audit baseline mismatch");
  for (auto v : audit["inputs"].toArray()) {
    auto o = v.toObject();
    inputHashes[o["path"].toString()] = o["sha256"].toString();
  }
  for (QString name :
       {"metadata/setup_manifest.json", "metadata/multisource.json"})
    if (inputHashes.count(name))
      require(shaFile(data + "/" + name) == inputHashes[name],
              "Acquisition geometry changed; rerun the audit");
  auto m = arrivalAudit["material"].toObject(),
       sp = audit["source_parameters"].toObject();
  material = {m["cl_m_s"].toDouble(),
              m["cs_m_s"].toDouble(),
              m["cr_m_s"].toDouble(),
              m["thickness_m"].toDouble(),
              sp["line_half_width_m"].toDouble(),
              sp["pulse_fwhm_s"].toDouble(),
              .075 - sp["absorber_width_m"].toDouble()};
  require(material.cl > material.cs && material.cs > 0 && material.h > 0,
          "Invalid saved material");
  auto ss = audit["sources"].toArray();
  require(ss.size() == 5, "Expected five source positions");
  for (int i = 0; i < 5; ++i) {
    auto o = ss[i].toObject();
    sources[i] = {o["id"].toString(), o["id"].toString().replace('_', ' '),
                  o["x_mm"].toDouble() / 1000, o["source_trigger_s"].toDouble(),
                  material.h - o["xy_m"].toArray()[1].toDouble()};
  }
  for (auto v : audit["receivers"].toArray()) {
    auto o = v.toObject();
    receivers.push_back({o["receiver_id"].toString().mid(1).toInt(),
                         int(GroupNames.indexOf(o["group"].toString())),
                         o["group_index"].toInt(), o["x_m"].toDouble(),
                         material.h - o["y_m"].toDouble()});
  }
  std::sort(receivers.begin(), receivers.end(),
            [](auto a, auto b) { return a.id < b.id; });
  require(!receivers.empty() && receivers.size() <= Receivers,
          "Unsupported receiver count (maximum 337)");
  material.receiverCount = receiverCount();
  for (int i = 0; i < receiverCount(); ++i)
    require(receivers[i].id == i && receivers[i].group >= 0,
            "Invalid receiver mapping");
  auto tips = audit["defect_benchmark_only"].toObject()["tips_xy_mm"].toArray();
  for (int i = 0; i < 2; ++i)
    truth[i] = {tips[i].toArray()[0].toDouble(),
                material.h * 1000 - tips[i].toArray()[1].toDouble()};
  if (!bundled)
    geometry = readJson(data + "/metadata/setup_manifest.json")["polygons_xy_m"]
                   .toObject();
  for (auto rec : audit["records"].toArray()) {
    auto o = rec.toObject();
    for (auto group : o["groups"].toArray()) {
      auto g = group.toObject();
      auto &ranges =
          repairRanges[o["source"].toString() + "/" + o["case"].toString() +
                       "/" + g["group"].toString()];
      for (auto gap : g["repair_intervals"].toArray())
        ranges.emplace_back(gap.toObject()["start_us"].toDouble() * 1e-6 - Dt,
                            gap.toObject()["end_us"].toDouble() * 1e-6 + Dt);
    }
  }
  for (auto pair : std::vector<std::pair<int, int>>{
           {1, 1}, {2, 1}, {1, 2}, {2, 2}, {1, 3}, {3, 1}, {2, 3}, {3, 2}}) {
    int n = pair.first + pair.second;
    for (int j = 0; j < (1 << n); ++j) {
      QString mode;
      for (int k = n - 1; k >= 0; --k)
        mode += (j & (1 << k)) ? 'S' : 'L';
      QString
          a = mode.left(pair.first),
          b = mode.mid(pair.first),
          id = QString("%1-%2_%3").arg(pair.first).arg(pair.second).arg(mode);
      paths.push_back({id, a, b, a + " | " + b, {}, {}, {}});
    }
  }
  addGroovePaths();
  for (int s = 0; s < 5; ++s)
    for (int r = 0; r < receiverCount(); ++r) {
      double dx = std::abs(receivers[r].x - sources[s].x), t = sources[s].t0;
      auto ls = ray({sources[s].x, 0}, {receivers[r].x, 0}, "LS", material);
      nuisance[s][r] = {
          (t + dx / material.cr) * 1e6,
          (t + dx / material.cl) * 1e6,
          (t + dx / material.cs) * 1e6,
          (t + std::hypot(dx, 2 * material.h) / material.cl) * 1e6,
          (t + ls.time) * 1e6,
          (t + std::hypot(dx, 2 * material.h) / material.cs) * 1e6};
      nuisanceSL[s][r] = nuisance[s][r][4];
      if (datasetKind == "partial") {
        const Point a{sources[s].x, sources[s].z},
            b{receivers[r].x, receivers[r].z};
        const double missing = std::numeric_limits<double>::quiet_NaN();
        nuisanceSL[s][r] = (t + trace(a, b, "SL").time) * 1e6;
        const bool visible = visibleInSteel(a, b);
        const double direct = std::hypot(a.x - b.x, a.z - b.z);
        nuisance[s][r] = {visible && std::abs(a.z - b.z) < 1e-12
                              ? (t + dx / material.cr) * 1e6
                              : missing,
                          visible ? (t + direct / material.cl) * 1e6 : missing,
                          visible ? (t + direct / material.cs) * 1e6 : missing,
                          (t + trace(a, b, "LL").time) * 1e6,
                          (t + trace(a, b, "LS").time) * 1e6,
                          (t + trace(a, b, "SS").time) * 1e6};
      }
      for (const auto &echo : shearReverberations(
               {sources[s].x, sources[s].z}, {receivers[r].x, receivers[r].z}))
        nuisanceSSSS[s][r].push_back((t + echo.time) * 1e6);
    }
}
const Path &Engine::path(const QString &id) const {
  auto it = std::find_if(paths.begin(), paths.end(),
                         [&](auto &p) { return p.id == id; });
  require(it != paths.end(), "Unknown wave path: " + id);
  return *it;
}
std::shared_ptr<Vec> Engine::raw(int si, QString signal, QString component) {
  require(si >= 0 && si < 5 &&
              QStringList{"lof", "reference", "difference"}.contains(signal) &&
              QStringList{"vx", "vy"}.contains(component),
          "Invalid recording request");
  QString key = sources[si].id + "/" + signal + "/" + component;
  if (rawCache.count(key))
    return rawCache[key];
  auto result = std::make_shared<Vec>(receiverCount() * Samples);
  if (signal == "difference") {
    auto a = raw(si, "lof", component), b = raw(si, "reference", component);
    for (size_t i = 0; i < result->size(); ++i)
      (*result)[i] = (*a)[i] - (*b)[i];
    rawCache[key] = result;
    return result;
  }
  for (int group = 0; group < 3; ++group) {
    if (cancel)
      throw Cancelled();
    QString relative = "sources/" + sources[si].id + "/" + signal +
                       "/receivers/" + GroupNames[group] + "_" + component +
                       ".csv";
    require(inputHashes.count(relative),
            "Recording absent from audit: " + relative);
    QString filename = data + "/" + relative;
    QFileInfo info(filename);
    require(info.isFile(), "Missing recording: " + relative);
    std::vector<int> ids;
    for (auto &r : receivers)
      if (r.group == group)
        ids.push_back(r.id);
    qsizetype bytes = qsizetype(ids.size()) * Samples * sizeof(double);
    QString stem = cache + "/" + sha(relative.toUtf8());
    QByteArray binary;
    bool loaded = false;
    if (QFileInfo::exists(stem + ".json") && QFileInfo::exists(stem + ".bin")) {
      try {
        auto meta = readJson(stem + ".json");
        if (meta["source_sha256"] == inputHashes[relative] &&
            meta["bytes"].toInteger() == info.size() &&
            meta["modified_ms"].toInteger() ==
                info.lastModified().toMSecsSinceEpoch()) {
          QFile f(stem + ".bin");
          if (f.open(QIODevice::ReadOnly)) {
            binary = f.readAll();
            loaded = binary.size() == bytes &&
                     sha(binary) == meta["binary_sha256"].toString();
          }
        }
      } catch (...) {
        loaded = false;
      }
    }
    if (loaded)
      ++cacheReads;
    else {
      QFile f(filename);
      require(f.open(QIODevice::ReadOnly), "Cannot read " + filename);
      auto csv = f.readAll();
      require(sha(csv) == inputHashes[relative],
              "Recording checksum changed: " + relative +
                  ". Re-audit before imaging.");
      const char *ptr = csv.constData();
      const char *end = ptr + csv.size();
      const char *newline = std::find(ptr, end, '\n');
      require(newline != end, "CSV header missing");
      QByteArray header(ptr, newline - ptr);
      header = header.trimmed();
      if (header.startsWith("\xef\xbb\xbf"))
        header = header.mid(3);
      QByteArray expected = "time_s";
      for (int id : ids)
        expected += ",R" + QByteArray::number(id);
      require(header == expected, "CSV receiver order changed: " + relative);
      ptr = newline + 1;
      std::vector<double> parsed(ids.size() * Samples);
      for (int t = 0; t < Samples; ++t) {
        for (int c = -1; c < int(ids.size()); ++c) {
          while (ptr < end &&
                 (*ptr == '\r' || *ptr == '\n' || *ptr == ' ' || *ptr == '\t'))
            ++ptr;
          double value = 0;
          auto conv = std::from_chars(ptr, end, value);
          require(conv.ec == std::errc() && std::isfinite(value),
                  "Invalid CSV sample: " + relative);
          ptr = conv.ptr;
          if (c < 0)
            require(std::abs(value - t * Dt) <= 1e-15, "CSV time grid changed");
          else
            parsed[size_t(c) * Samples + t] = value;
          if (c + 1 < int(ids.size())) {
            require(ptr < end && *ptr == ',', "CSV column count changed");
            ++ptr;
          }
        }
        if (t + 1 < Samples)
          require(ptr < end && (*ptr == '\r' || *ptr == '\n'),
                  "CSV row length changed");
      }
      while (ptr < end && std::isspace(static_cast<unsigned char>(*ptr)))
        ++ptr;
      require(ptr == end, "Extra CSV data");
      binary = QByteArray(reinterpret_cast<const char *>(parsed.data()), bytes);
      QSaveFile out(stem + ".bin");
      require(out.open(QIODevice::WriteOnly), "Cannot create recording cache");
      out.write(binary);
      require(out.commit(), "Cannot commit recording cache");
      writeJson(stem + ".json",
                {{"source_sha256", inputHashes[relative]},
                 {"bytes", info.size()},
                 {"modified_ms", info.lastModified().toMSecsSinceEpoch()},
                 {"binary_sha256", sha(binary)}});
      ++csvParses;
    }
    for (size_t i = 0; i < ids.size(); ++i)
      std::memcpy(result->data() + size_t(ids[i]) * Samples,
                  binary.constData() + i * Samples * sizeof(double),
                  Samples * sizeof(double));
    ++loads[relative];
  }
  rawCache[key] = result;
  return result;
}
std::shared_ptr<CVec> Engine::filtered(int si, QString signal,
                                       const Config &c) {
  if (signalProvider)
    return signalProvider(si, signal, c);
  QString key = QString("%1/%2/%3/%4/%5")
                    .arg(si)
                    .arg(signal, c.component)
                    .arg(c.low, 0, 'g', 17)
                    .arg(c.high, 0, 'g', 17);
  if (filteredCache.count(key))
    return filteredCache[key];
  // Bound retained filtered data to two complete paired five-source
  // acquisitions.
  if (filteredCache.size() >= 20)
    filteredCache.clear();
  auto values = raw(si, signal, c.component);
  auto h = kernel(c.low, c.high);
  h.resize(8192);
  fft(h);
  auto out = std::make_shared<CVec>(receiverCount() * Samples);
  parallel(
      receiverCount(), c.threads,
      [&](int r) {
        CVec a(8192);
        for (int t = 0; t < Samples; ++t)
          a[t] = (*values)[size_t(r) * Samples + t];
        fft(a);
        for (int t = 0; t < 8192; ++t)
          a[t] *= h[t];
        fft(a, true);
        std::copy_n(a.begin() + Half, Samples,
                    out->begin() + size_t(r) * Samples);
      },
      &cancel);
  filteredCache[key] = out;
  ++filterBuilds;
  return out;
}
std::shared_ptr<Curves> Engine::curves(const Config &c) {
  QString key =
      QString::number(c.low, 'g', 17) + "/" + QString::number(c.high, 'g', 17);
  if (curveCache.count(key))
    return curveCache[key];
  auto out = std::make_shared<Curves>();
  auto h = kernel(c.low, c.high);
  h.resize(16384);
  fft(h);
  Vec omega(8191), response(8191);
  for (int k = 1; k < 8192; ++k) {
    omega[k - 1] = 2 * Pi * k / (16384 * Dt);
    response[k - 1] = std::norm(h[k]);
  }
  double sigma = material.pulse / (2 * std::sqrt(2 * std::log(2.)));
  for (int mode = 0; mode < 2; ++mode) {
    out->launch[mode].resize(1601);
    for (int d = 0; d < 2; ++d)
      out->receive[mode][d].resize(1601);
    char letter = mode ? 'S' : 'L';
    double speed = mode ? material.cs : material.cl;
    parallel(
        1601, c.threads,
        [&](int i) {
          double theta = i * .05 * Pi / 180, p = std::sin(theta) / speed;
          Complex ql = std::conj(std::sqrt(
                      Complex(1 / (material.cl * material.cl) - p * p, 0))),
                  qs = std::conj(std::sqrt(
                      Complex(1 / (material.cs * material.cs) - p * p, 0)));
          auto amp = solve(traction(p, ql, 'L', material),
                           traction(p, qs, 'S', material), {1., 0.});
          double integral = 0, a = p * p * material.width * material.width / 4 +
                                   sigma * sigma;
          for (int k = 0; k < 8191; ++k)
            integral += omega[k] * response[k] *
                        std::exp(-omega[k] * omega[k] * a) *
                        (k == 0 || k == 8190 ? .5 : 1);
          integral *= omega[1] - omega[0];
          out->launch[mode][i] =
              std::cos(theta) * std::abs(amp[mode] * p) * std::sqrt(integral);
          auto r = reflect(theta, letter, material);
          for (int d = 0; d < 2; ++d)
            out->receive[mode][d][i] = std::abs(r.surface[d]);
        },
        &cancel);
    double peak =
        *std::max_element(out->launch[mode].begin(), out->launch[mode].end());
    for (auto &v : out->launch[mode])
      v /= peak;
    for (int d = 0; d < 2; ++d) {
      double max = *std::max_element(out->receive[mode][d].begin(),
                                     out->receive[mode][d].end());
      out->norm[mode][d] = max;
      for (auto &v : out->receive[mode][d])
        v /= max;
    }
  }
  curveCache[key] = out;
  return out;
}
static double gain(double angle, const Vec &curve) {
  double deg = std::abs(angle) * 180 / Pi;
  if (deg > 80)
    return 0;
  double f = deg / .05;
  int i = std::min(1599, int(f));
  return curve[i] + (f - i) * (curve[i + 1] - curve[i]);
}
bool Engine::supportClean(double time, int si, QString signal, int group,
                          const Config &c) const {
  int j = int(std::floor(time / Dt));
  double a = (j - signalHalf) * Dt, b = (j + 1 + signalHalf) * Dt;
  if (!signalSampleSupported(j))
    return false;
  if (!c.clean)
    return true;
  QStringList cases = (signal == "difference" || c.repairs == "paired")
                          ? QStringList{"lof", "reference"}
                          : QStringList{signal};
  for (auto signalCase : cases) {
    auto it = repairRanges.find(sources[si].id + "/" + signalCase + "/" +
                                GroupNames[group]);
    if (it != repairRanges.end())
      for (auto gap : it->second)
        if (a <= gap.second && b >= gap.first)
          return false;
  }
  return true;
}
double Engine::keep(double time, int si, int ri, const Config &c,
                    const Screening *screen) const {
  if (screen && screen->rejects(time, si, ri))
    return 0;
  double g = 1;
  for (int k = 0; k < (datasetKind == "partial" ? 7 : 6); ++k) {
    const double centre = k == 6 ? nuisanceSL[si][ri] : nuisance[si][ri][k];
    if (c.mutes[k == 6 ? 4 : k] && std::isfinite(centre)) {
      double d = std::abs(time * 1e6 - centre) / c.muteWidth,
             ramp = std::clamp((d - .75) / .25, 0., 1.);
      g = std::min(g, .5 - .5 * std::cos(Pi * ramp));
    }
  }
  return g;
}
static bool manual(const Config &c, Point p, const Receiver &r,
                   const std::array<bool, Receivers> &ids) {
  return ids[r.id] && c.groups[r.group] &&
         (c.aperture == "all" || (c.aperture == "left" && r.x < p.x) ||
          (c.aperture == "right" && r.x >= p.x) ||
          (c.aperture == "even" && r.id % 2 == 0) ||
          (c.aperture == "odd" && r.id % 2 == 1));
}
void Engine::prepareGeometry(const Config &c,
                             const std::vector<Point> &points) {
  // Keep dense grids bounded: compute each pixel's geometry inside its worker
  // instead of retaining several gigabytes of receiver rays.
  if (points.size() * receiverCount() * sizeof(Ray) > 256 * 1024 * 1024) {
    geo = Geometry{};
    ++geometryBuilds;
    return;
  }
  QString key = c.path + "/" +
                sha(QByteArray(reinterpret_cast<const char *>(points.data()),
                               points.size() * sizeof(Point)));
  if (geo.key == key)
    return;
  Geometry next;
  next.key = key;
  next.receiving.resize(points.size() * receiverCount());
  for (auto &v : next.incoming)
    v.resize(points.size());
  auto route = path(c.path);
  parallel(
      int(points.size()), c.threads,
      [&](int i) {
        for (int s = 0; s < 5; ++s)
          next.incoming[s][i] =
              trace({sources[s].x, sources[s].z}, points[i], route.before, false, false, route.beforeVia);
        for (int r = 0; r < receiverCount(); ++r)
          next.receiving[size_t(i) * receiverCount() + r] = trace(
              points[i], {receivers[r].x, receivers[r].z}, route.after, true, false, route.afterVia);
      },
      &cancel);
  geo = std::move(next);
  ++geometryBuilds;
}
std::array<Support, 5> Engine::support(const Config &c, Point point,
                                       QString signal) {
  auto curve = curves(c);
  auto route = path(c.path);
  int d = c.component == "vy" ? 1 : 0, mode = route.after.back() == 'S';
  std::array<bool, Receivers> ids{};
  for (int r : c.receivers)
    ids[r] = true;
  Point p{point.x / 1000, point.z / 1000};
  std::array<Support, 5> result;
  for (int s = 0; s < 5; ++s) {
    auto &out = result[s];
    auto inc = trace({sources[s].x, sources[s].z}, p, route.before, false, false, route.beforeVia);
    if (!inc.valid)
      continue;
    out.gain = gain(inc.first, curve->launch[route.before.front() == 'S']);
    out.angle = inc.first * 180 / Pi;
    double eligible = 0, w = 0, w2 = 0;
    for (int r = 0; r < receiverCount(); ++r) {
      auto &rr = receivers[r];
      auto rec = trace(p, {rr.x, rr.z}, route.after, true, false, route.afterVia);
      if (!rec.valid)
        continue;
      auto surface = reflect(rec.last, route.after.back().toLatin1(), material)
                         .surface[d] /
                     curve->norm[mode][d];
      auto wall = inc.wall * rec.wall;
      bool allowed = out.gain >= c.gate && std::abs(wall) >= c.wall &&
                     std::abs(surface) >= c.gate &&
                     std::abs(rec.last) * 180 / Pi <= 80 &&
                     sourceDistance(s, r) >= c.offset / 1000 - 1e-14 &&
                     std::abs(rr.x) < material.active &&
                     std::abs(p.x) < material.active && manual(c, p, rr, ids) &&
                     (!c.sourceReceivers.count(s) ||
                      std::find(c.sourceReceivers.at(s).begin(),
                                c.sourceReceivers.at(s).end(),
                                r) != c.sourceReceivers.at(s).end());
      if (!allowed)
        continue;
      ++eligible;
      double time = sources[s].t0 + inc.time + rec.time;
      if (!supportClean(time, s, signal, rr.group, c))
        continue;
      double weight = c.processing == "weighted" ? std::abs(wall * surface) : 1;
      ++out.count;
      w += weight;
      w2 += weight * weight;
      out.ids.push_back(r);
    }
    out.effective = w2 ? w * w / w2 : 0;
    out.valid = out.count >= std::max(double(c.minimumReceivers),
                                      std::ceil(c.fraction * eligible)) &&
                out.effective >= c.minimumEffective - 1e-10;
  }
  return result;
}
QJsonObject Result::timing() const {
  return {{"total_s", seconds},
          {"raw_load_s", loadSeconds},
          {"filter_and_curves_s", filterSeconds},
          {"geometry_s", geometrySeconds},
          {"sum_s", sumSeconds}};
}
std::shared_ptr<Result>
Engine::reconstruct(const Config &c, std::function<void(QString)> progress,
                    const std::vector<Point> *points) {
  c.validate(material);
  path(c.path);
  QElapsedTimer total, phase;
  total.start();
  phase.start();
  auto tell = [&](QString t) {
    if (progress)
      progress(t);
    if (cancel)
      throw Cancelled();
  };
  auto out = std::make_shared<Result>();
  out->config = c;
  out->x = c.xs();
  out->z = c.zs();
  if (points)
    out->points = *points;
  else
    for (double z : out->z)
      for (double x : out->x)
        out->points.push_back({x / 1000, z / 1000});
  size_t n = out->points.size();
  for (auto signal : c.cases)
    for (int s : c.sources) {
      tell("Loading " + sources[s].name + " / " + signal);
      raw(s, signal, c.component);
    }
  out->loadSeconds = phase.nsecsElapsed() / 1e9;
  phase.restart();
  auto curve = curves(c);
  std::map<QString, std::vector<std::shared_ptr<CVec>>> signals;
  for (auto signal : c.cases)
    for (int s : c.sources) {
      tell("Filtering " + sources[s].name + " / " + signal);
      signals[signal].push_back(filtered(s, signal, c));
    }
  for (auto signal : c.cases)
    out->screening[signal] = screeningFor(c, signal);
  out->filterSeconds = phase.nsecsElapsed() / 1e9;
  phase.restart();
  tell("Computing paths across the full image");
  prepareGeometry(c, out->points);
  out->geometrySeconds = phase.nsecsElapsed() / 1e9;
  phase.restart();
  std::array<bool, Receivers> ids{};
  for (int r : c.receivers)
    ids[r] = true;
  std::array<std::array<bool, Receivers>, 5> frozen{};
  if (c.freeze) {
    auto p = support(c, c.anchor, c.cases.front());
    for (int s = 0; s < 5; ++s)
      for (int r : p[s].ids)
        frozen[s][r] = true;
  }
  auto route = path(c.path);
  int column = c.component == "vy" ? 1 : 0, mode = route.after.back() == 'S';
  bool weighted = c.processing == "weighted";
  for (auto signal : c.cases) {
    auto &planes = out->planes[signal];
    for (auto &p : planes) {
      p.amplitude.resize(n, std::nan(""));
      p.count.resize(n);
      p.cf.resize(n);
      p.retained.resize(n, std::nan(""));
      p.screened.resize(n, std::nan(""));
      p.receiverUsage.resize(5 * receiverCount());
      p.stats.resize(n * c.sources.size());
      p.effective.resize(n * c.sources.size());
      p.valid.resize(n * c.sources.size());
      p.means.resize(n * c.sources.size());
    }
  }
  // Resolve Qt strings, recording buffers and repair histories once. The
  // parallel pixel loop contains no QString ref-count traffic or map lookups.
  struct CaseContext {
    std::array<Plane, 3> *planes;
    const Screening *screen;
    std::vector<const CVec *> signals;
    std::array<std::array<std::vector<std::pair<double, double>>, 3>, 5> gaps;
  };
  std::vector<CaseContext> contexts;
  for (const auto &signal : c.cases) {
    CaseContext ctx;
    ctx.planes = &out->planes.at(signal);
    ctx.screen = out->screening.at(signal).get();
    for (const auto &data : signals.at(signal))
      ctx.signals.push_back(data.get());
    if (c.clean) {
      QStringList repairCases = signal == "difference" || c.repairs == "paired"
                                    ? QStringList{"lof", "reference"}
                                    : QStringList{signal};
      for (int source : c.sources)
        for (int group = 0; group < 3; ++group)
          for (const auto &repairCase : repairCases) {
            auto found =
                repairRanges.find(sources[source].id + "/" + repairCase + "/" +
                                  GroupNames[group]);
            if (found != repairRanges.end())
              ctx.gaps[source][group].insert(ctx.gaps[source][group].end(),
                                             found->second.begin(),
                                             found->second.end());
          }
    }
    contexts.push_back(std::move(ctx));
  }
  std::array<std::array<bool, Receivers>, 5> offsetAllowed{};
  std::array<std::array<bool, Receivers>, 5> sourceIdsAllowed{};
  std::array<bool, Receivers> manuallyEnabled{};
  for (int r = 0; r < receiverCount(); ++r) {
    manuallyEnabled[r] = ids[r] && c.groups[receivers[r].group];
    for (int source : c.sources) {
      auto found = c.sourceReceivers.find(source);
      sourceIdsAllowed[source][r] =
          found == c.sourceReceivers.end() ||
          std::find(found->second.begin(), found->second.end(), r) !=
              found->second.end();
      offsetAllowed[source][r] =
          sourceDistance(source, r) >= c.offset / 1000 - 1e-14 &&
          std::abs(receivers[r].x) < material.active;
    }
  }
  const int sector =
      QStringList{"all", "left", "right", "even", "odd"}.indexOf(c.aperture);
  tell("Summing admitted measurements");
  std::mutex usageMutex;
  using UsedIds = std::array<std::array<std::bitset<Receivers>, 5>, 3>;
  parallel(
      int((n + 31) / 32), c.threads,
      [&](int block) {
        std::vector<std::array<std::vector<quint64>, 3>> usage(c.cases.size());
        for (auto &entry : usage)
          for (auto &variant : entry)
            variant.resize(5 * receiverCount());
        std::vector<UsedIds> used(c.cases.size());
        for (size_t i = size_t(block) * 32;
             i < std::min(n, size_t(block + 1) * 32); ++i) {
          if (cancel)
            throw Cancelled();
          for (auto &entry : used)
            for (auto &variant : entry)
              for (auto &source : variant)
                source.reset();
          Point p = out->points[i];
          if (datasetKind == "partial" && !containsSteel(p))
            continue;
          std::array<Ray, Receivers> localReceiving;
          std::array<Ray, 5> localIncoming;
          bool cached = !geo.receiving.empty();
          if (!cached) {
            for (int r = 0; r < receiverCount(); ++r)
              localReceiving[r] =
                  trace(p, {receivers[r].x, receivers[r].z}, route.after, true, false, route.afterVia);
            for (int s : c.sources)
              localIncoming[s] =
                  trace({sources[s].x, sources[s].z}, p, route.before, false, false, route.beforeVia);
          }
          std::array<Complex, Receivers> surface;
          for (int r = 0; r < receiverCount(); ++r) {
            auto &rec = cached ? geo.receiving[size_t(i) * receiverCount() + r]
                               : localReceiving[r];
            if (!rec.valid) {
              surface[r] = 0;
              continue;
            }
            surface[r] =
                (column ? -1. : 1.) *
                reflect(rec.last, route.after.back().toLatin1(), material)
                    .surface[column] /
                curve->norm[mode][column];
          }
          for (size_t oi = 0; oi < c.sources.size(); ++oi) {
            int s = c.sources[oi];
            auto &inc = cached ? geo.incoming[s][i] : localIncoming[s];
            if (!inc.valid)
              continue;
            double launch =
                gain(inc.first, curve->launch[route.before.front() == 'S']);
            if (launch < c.gate)
              continue;
            for (int r = 0; r < receiverCount(); ++r) {
              auto &rr = receivers[r];
              auto &rec = cached
                              ? geo.receiving[size_t(i) * receiverCount() + r]
                              : localReceiving[r];
              if (!rec.valid)
                continue;
              auto wall = inc.wall * rec.wall;
              if (std::abs(wall) < c.wall || std::abs(surface[r]) < c.gate ||
                  std::abs(rec.last) * 180 / Pi > 80 || !offsetAllowed[s][r] ||
                  std::abs(p.x) >= material.active)
                continue;
              bool selected =
                  c.freeze ? frozen[s][r]
                           : manuallyEnabled[r] && sourceIdsAllowed[s][r] &&
                                 (sector == 0 || (sector == 1 && rr.x < p.x) ||
                                  (sector == 2 && rr.x >= p.x) ||
                                  (sector == 3 && r % 2 == 0) ||
                                  (sector == 4 && r % 2 == 1));
              double time = sources[s].t0 + inc.time + rec.time;
              double fraction = time / Dt;
              int j = int(std::floor(fraction));
              fraction -= j;
              bool completeFilter = signalSampleSupported(j);
              double supportFirst = (j - signalHalf) * Dt,
                     supportLast = (j + 1 + signalHalf) * Dt;
              j = std::clamp(j, 0, Samples - 2);
              double g = keep(time, s, r, c);
              Complex response = weighted ? wall * surface[r] : Complex(1.);
              double weight = std::abs(response);
              int caseIndex = 0;
              for (const auto &ctx : contexts) {
                bool rejected = ctx.screen && ctx.screen->rejects(time, s, r);
                auto &planes = *ctx.planes;
                bool usable = completeFilter;
                if (usable && c.clean)
                  for (const auto &gap : ctx.gaps[s][rr.group])
                    if (supportFirst <= gap.second &&
                        supportLast >= gap.first) {
                      usable = false;
                      break;
                    }
                Complex value = 0;
                if (usable) {
                  const auto &a = *ctx.signals[oi];
                  value = (1 - fraction) * a[size_t(r) * Samples + j] +
                          fraction * a[size_t(r) * Samples + j + 1];
                }
                for (int v = 0; v < 3; ++v) {
                  if (v && !selected)
                    continue;
                  auto &st = planes[v].stats[oi * n + i];
                  ++st.eligible;
                  if (usable) {
                    ++st.candidates;
                    if (v == 2 && rejected)
                      ++st.screened;
                  }
                  if (!usable || (v == 2 && (g <= 1e-12 || rejected)))
                    continue;
                  Complex val = value * (v == 2 ? g : 1.);
                  ++st.count;
                  used[caseIndex][v][s].set(r);
                  st.weight += weight;
                  st.weight2 += weight * weight;
                  st.sum += val * std::conj(response);
                  st.power += weight * std::norm(val);
                }
                ++caseIndex;
              }
            }
          }
          int caseIndex = 0;
          for (const auto &ctx : contexts) {
            auto &planes = *ctx.planes;
            for (int v = 0; v < 3; ++v) {
              auto &plane = planes[v];
              double sumPower = 0, cfMean = 0, candidates = 0, retained = 0,
                     screened = 0;
              Complex sum = 0;
              int count = 0;
              for (size_t oi = 0; oi < c.sources.size(); ++oi) {
                size_t k = oi * n + i;
                auto &st = plane.stats[k];
                candidates += st.candidates;
                retained += st.count;
                screened += st.screened;
                double eff =
                           st.weight2 ? st.weight * st.weight / st.weight2 : 0,
                       cf = st.weight * st.power > 0
                                ? std::clamp(std::norm(st.sum) /
                                                 (st.weight * st.power),
                                             0., 1.)
                                : 0;
                bool valid =
                    st.count >= std::max(double(c.minimumReceivers),
                                         std::ceil(c.fraction * st.eligible)) &&
                    eff >= c.minimumEffective - 1e-10;
                Complex mean = st.weight ? st.sum / st.weight : Complex(0);
                if (v == 2 && c.coherence)
                  mean *= std::pow(cf, c.cfPower);
                plane.effective[k] = eff;
                plane.valid[k] = valid;
                plane.means[k] = mean;
                if (valid) {
                  ++count;
                  sumPower += std::norm(mean);
                  sum += mean;
                  cfMean += cf;
                }
              }
              if (candidates > 0) {
                plane.retained[i] = 100 * retained / candidates;
                plane.screened[i] = 100 * screened / candidates;
              }
              plane.count[i] = count;
              plane.cf[i] = count ? cfMean / count : 0;
              int required = c.policy == "fixed"     ? int(c.sources.size())
                             : c.sources.size() == 1 ? 1
                                                     : c.minimumSources;
              if (count >= required) {
                for (size_t oi = 0; oi < c.sources.size(); ++oi)
                  if (plane.valid[oi * n + i]) {
                    int source = c.sources[oi];
                    for (int r = 0; r < receiverCount(); ++r)
                      if (used[caseIndex][v][source].test(r))
                        ++usage[caseIndex][v][source * receiverCount() + r];
                  }
                plane.amplitude[i] = c.combination == "rms"
                                         ? std::sqrt(sumPower / count)
                                         : std::abs(sum) / count;
              }
            }
            ++caseIndex;
          }
        }
        std::lock_guard<std::mutex> lock(usageMutex);
        for (int ci = 0; ci < c.cases.size(); ++ci)
          for (int v = 0; v < 3; ++v)
            for (int k = 0; k < 5 * receiverCount(); ++k)
              (*contexts[ci].planes)[v].receiverUsage[k] += usage[ci][v][k];
      },
      &cancel);
  for (auto &entry : out->planes)
    for (auto &plane : entry.second)
      plane.supportedPixels =
          std::count_if(plane.amplitude.begin(), plane.amplitude.end(),
                        [](double v) { return std::isfinite(v); });
  out->sumSeconds = phase.nsecsElapsed() / 1e9;
  out->seconds = total.nsecsElapsed() / 1e9;
  return out;
}
std::shared_ptr<Inspection> Engine::inspect(const Config &c, int s, int rid,
                                            QString signal, QString competitor,
                                            bool hasPoint, Point point) {
  c.validate(material);
  require(rid >= 0 && rid < receiverCount(), "Invalid receiver");
  QElapsedTimer t;
  t.start();
  auto out = std::make_shared<Inspection>();
  out->config = c;
  out->source = s;
  out->receiver = rid;
  out->signal = signal;
  out->competitor = competitor;
  out->hasPoint = hasPoint;
  out->point = point;
  out->raw = *raw(s, signal, c.component);
  out->filtered = *filtered(s, signal, c);
  auto screenConfig = c;
  if (std::find(screenConfig.sources.begin(), screenConfig.sources.end(), s) ==
      screenConfig.sources.end())
    screenConfig.sources.push_back(s);
  out->screening = screeningFor(screenConfig, signal);
  out->times.resize(receiverCount());
  out->competing.resize(receiverCount());
  for (auto &a : out->arrivals)
    a.resize(receiverCount());
  out->orderedIds.resize(receiverCount());
  std::iota(out->orderedIds.begin(), out->orderedIds.end(), 0);
  std::sort(out->orderedIds.begin(), out->orderedIds.end(),
            [&](int a, int b) { return receivers[a].x < receivers[b].x; });
  if (hasPoint) {
    require(point.z > 0 && point.z < material.h * 1000,
            "Inspection depth must lie inside the plate");
    out->support = support(c, point, signal);
  }
  for (int r = 0; r < receiverCount(); ++r) {
    for (int k = 0; k < 6; ++k)
      out->arrivals[k][r] = nuisance[s][r][k];
    out->arrivals[6][r] = nuisanceSL[s][r];
    if (hasPoint) {
      for (int v = 0; v < 2; ++v) {
        auto route = path(v ? competitor : c.path);
        Point p{point.x / 1000, point.z / 1000};
        auto a = trace({sources[s].x, sources[s].z}, p, route.before, false,
                       r == rid, route.beforeVia),
             b = trace(p, {receivers[r].x, receivers[r].z}, route.after, true,
                       r == rid, route.afterVia);
        (v ? out->competing : out->times)[r] =
            (sources[s].t0 + a.time + b.time) * 1e6;
        if (r == rid && a.valid && b.valid) {
          auto &list = v ? out->competingPath : out->path;
          list = a.points;
          list.insert(list.end(), b.points.begin() + 1, b.points.end());
          for (auto &p : list) {
            p.x *= 1000;
            p.z *= 1000;
          }
        }
      }
    }
  }
  out->seconds = t.nsecsElapsed() / 1e9;
  return out;
}
QJsonObject Engine::provenance() const {
  QJsonObject hashes;
  for (auto &[p, h] : inputHashes)
    hashes[p] = h;
  QJsonObject used;
  for (auto &[p, n] : loads)
    used[p] = n;
  auto result = QJsonObject{
      {"dataset_name", datasetName},
      {"dataset_kind", datasetKind},
      {"dataset_folder", data},
      {"dataset_manifest_sha256", manifestHash},
      {"receiver_count", receiverCount()},
      {"sources", audit["sources"]},
      {"receivers", audit["receivers"]},
      {"geometry", geometry},
      {"depth_origin", "Depth = plate thickness - saved model y; zero is the "
                       "original plate OD"},
      {"ray_model",
       datasetKind == "partial"
           ? "Homogeneous bulk rays inside steel; flat backwall and exposed "
             "OD/fill-top reflections only. Earliest valid flat upper branch. "
             "No sloped groove reflection or corner diffraction. Rayleigh "
             "guides only along one unobstructed flat face."
           : "Homogeneous bulk rays with flat OD/ID reflections"},
      {"model_evidence", arrivalAudit["model_evidence"]},
      {"input_fingerprint", fingerprint},
      {"input_checksums", hashes},
      {"recordings_accessed", used},
      {"material", arrivalAudit["material"]},
      {"implementation", "C++17 / Qt 6.8.3"},
      {"cache_policy", "Raw cache verified by payload SHA256 and expected "
                       "audited source SHA256; reuse requires unchanged source "
                       "size and mtime. Full input recheck is explicit."},
      {"saved_model_audit", arrivalAudit["model_inputs"]}};
  QString manifest =
      QCoreApplication::applicationDirPath() + "/implementation.json";
  if (QFileInfo::exists(manifest))
    result["implementation_sha256"] = readJson(manifest);
  return result;
}
// Small, dependency-free NPZ writer. Stored ZIP entries keep the exact float64
// and complex128 arrays and remain directly readable with numpy.load.
static quint32 crc32(const QByteArray &bytes) {
  static const auto table = [] {
    std::array<quint32, 256> a{};
    for (quint32 i = 0; i < 256; ++i) {
      quint32 c = i;
      for (int j = 0; j < 8; ++j)
        c = (c >> 1) ^ ((c & 1) ? 0xedb88320U : 0);
      a[i] = c;
    }
    return a;
  }();
  quint32 c = 0xffffffffU;
  for (unsigned char b : bytes)
    c = table[(c ^ b) & 255] ^ (c >> 8);
  return c ^ 0xffffffffU;
}
static void exportNpz(const QString &filename, const Result &result) {
  QSaveFile file(filename);
  require(file.open(QIODevice::WriteOnly), "Cannot export numerical arrays");
  QDataStream stream(&file);
  stream.setByteOrder(QDataStream::LittleEndian);
  struct Entry {
    QByteArray name;
    quint32 crc, size, offset;
  };
  std::vector<Entry> entries;
  auto write = [&](QByteArray name, QByteArray dtype, std::vector<size_t> shape,
                   const void *memory, size_t bytes) {
    QByteArray dims;
    for (size_t n : shape)
      dims += QByteArray::number(qulonglong(n)) + ",";
    QByteArray header = "{'descr': '" + dtype +
                        "', 'fortran_order': False, 'shape': (" + dims + ")}";
    header += QByteArray((16 - (10 + header.size() + 1) % 16) % 16, ' ');
    header += '\n';
    QByteArray data("\x93NUMPY\x01\x00", 8);
    data += char(header.size() & 255);
    data += char((header.size() >> 8) & 255);
    data += header;
    data.append(reinterpret_cast<const char *>(memory), qsizetype(bytes));
    name += ".npy";
    require(quint64(file.pos()) + quint64(data.size()) < 0xffffffffULL,
            "NPZ exceeds supported ZIP size");
    Entry e{name, crc32(data), quint32(data.size()), quint32(file.pos())};
    entries.push_back(e);
    stream << quint32(0x04034b50) << quint16(20) << quint16(0) << quint16(0)
           << quint16(0) << quint16(0) << e.crc << e.size << e.size
           << quint16(name.size()) << quint16(0);
    file.write(name);
    file.write(data);
  };
  write("x_mm", "<f8", {result.x.size()}, result.x.data(),
        result.x.size() * sizeof(double));
  write("depth_mm", "<f8", {result.z.size()}, result.z.data(),
        result.z.size() * sizeof(double));
  Vec points;
  for (auto p : result.points) {
    points.push_back(p.x);
    points.push_back(p.z);
  }
  write("points_m", "<f8", {result.points.size(), 2}, points.data(),
        points.size() * sizeof(double));
  for (auto signal : result.config.cases) {
    auto &planes = result.planes.at(signal);
    size_t n = result.points.size(), sources = result.config.sources.size();
    for (auto key :
         {"maps", "counts", "mean_coherence", "effective", "coherence",
          "receiver_retention_pct", "pulse_rejection_pct"}) {
      Vec values;
      for (auto &p : planes) {
        if (QString(key) == "coherence") {
          for (auto &s : p.stats)
            values.push_back(
                s.weight * s.power > 0
                    ? std::clamp(std::norm(s.sum) / (s.weight * s.power), 0.,
                                 1.)
                    : 0);
        } else {
          auto &a = QString(key) == "maps"                     ? p.amplitude
                    : QString(key) == "counts"                 ? p.count
                    : QString(key) == "mean_coherence"         ? p.cf
                    : QString(key) == "receiver_retention_pct" ? p.retained
                    : QString(key) == "pulse_rejection_pct"    ? p.screened
                                                               : p.effective;
          values.insert(values.end(), a.begin(), a.end());
        }
      }
      write(signal.toUtf8() + "_" + key, "<f8",
            QString(key) == "effective" || QString(key) == "coherence"
                ? std::vector<size_t>{3, sources, n}
                : std::vector<size_t>{3, n},
            values.data(), values.size() * sizeof(double));
    }
    std::vector<quint64> usage;
    for (const auto &p : planes)
      usage.insert(usage.end(), p.receiverUsage.begin(), p.receiverUsage.end());
    write(signal.toUtf8() + "_receiver_usage", "<u8",
          {3, 5, planes[0].receiverUsage.size() / 5}, usage.data(),
          usage.size() * sizeof(quint64));
    CVec means;
    std::vector<unsigned char> valid;
    for (auto &p : planes) {
      means.insert(means.end(), p.means.begin(), p.means.end());
      valid.insert(valid.end(), p.valid.begin(), p.valid.end());
    }
    write(signal.toUtf8() + "_means", "<c16", {3, sources, n}, means.data(),
          means.size() * sizeof(Complex));
    write(signal.toUtf8() + "_valid", "|u1", {3, sources, n}, valid.data(),
          valid.size());
  }
  quint32 begin = quint32(file.pos());
  for (auto e : entries) {
    stream << quint32(0x02014b50) << quint16(20) << quint16(20) << quint16(0)
           << quint16(0) << quint16(0) << quint16(0) << e.crc << e.size
           << e.size << quint16(e.name.size()) << quint16(0) << quint16(0)
           << quint16(0) << quint16(0) << quint32(0) << e.offset;
    file.write(e.name);
  }
  quint32 size = quint32(file.pos()) - begin;
  stream << quint32(0x06054b50) << quint16(0) << quint16(0)
         << quint16(entries.size()) << quint16(entries.size()) << size << begin
         << quint16(0);
  require(stream.status() == QDataStream::Ok && file.commit(),
          "Numerical export failed");
}
void Engine::exportResult(const QString &folder, const Result &r,
                          QJsonObject display) {
  QDir().mkpath(folder);
  QSaveFile file(folder + "/images.csv");
  require(file.open(QIODevice::WriteOnly), "Cannot export images");
  QTextStream f(&file);
  f.setRealNumberPrecision(17);
  f << "x_mm,depth_mm";
  QStringList variants{"wider", "selected", "processed"};
  for (auto signal : r.config.cases)
    for (int v = 0; v < 3; ++v)
      f << "," << signal << "_" << variants[v] << "_m_s," << signal << "_"
        << variants[v] << "_sources," << signal << "_" << variants[v]
        << "_mean_cf," << signal << "_" << variants[v]
        << "_receiver_retention_pct," << signal << "_" << variants[v]
        << "_pulse_rejection_pct";
  f << "\n";
  for (size_t i = 0; i < r.points.size(); ++i) {
    f << r.points[i].x * 1000 << "," << r.points[i].z * 1000;
    for (auto signal : r.config.cases)
      for (auto &p : r.planes.at(signal))
        f << "," << p.amplitude[i] << "," << p.count[i] << "," << p.cf[i] << ","
          << p.retained[i] << "," << p.screened[i];
    f << "\n";
  }
  f.flush();
  require(file.commit(), "Image export failed");
  exportNpz(folder + "/images.npz", r);
  QSaveFile usageFile(folder + "/receiver_usage.csv");
  require(usageFile.open(QIODevice::WriteOnly), "Cannot export receiver usage");
  QTextStream usageStream(&usageFile);
  usageStream << "signal,variant,source,source_x_mm,receiver_id,receiver_x_mm,"
                 "used_pixels,supported_image_pixels,source_depth_mm,receiver_"
                 "depth_mm\n";
  for (auto signal : r.config.cases)
    for (int v = 0; v < 3; ++v)
      for (int s = 0; s < 5; ++s)
        for (int id = 0; id < receiverCount(); ++id) {
          const auto &plane = r.planes.at(signal)[v];
          usageStream << signal << "," << variants[v] << "," << sources[s].id
                      << "," << sources[s].x * 1000 << ",R" << id << ","
                      << receivers[id].x * 1000 << ","
                      << plane.receiverUsage[s * receiverCount() + id] << ","
                      << plane.supportedPixels << "," << sources[s].z * 1000
                      << "," << receivers[id].z * 1000 << "\n";
        }
  usageStream.flush();
  require(usageFile.commit(), "Receiver usage export failed");
  QJsonObject windows;
  for (const auto &entry : r.screening)
    if (entry.second)
      windows[entry.first] = entry.second->json();
  writeJson(folder + "/pulse_windows.json", windows);
  auto meta = provenance();
  meta["schema"] = 3;
  meta["selection"] = r.config.json();
  meta["display"] = display;
  meta["timing"] = r.timing();
  meta["image_csv_sha256"] = shaFile(folder + "/images.csv");
  meta["image_npz_sha256"] = shaFile(folder + "/images.npz");
  meta["created_utc"] = QDateTime::currentDateTimeUtc().toString(Qt::ISODate);
  meta["interpretation"] =
      "Suppression is optional. Coherence-weighted scores are not calibrated "
      "amplitudes. Known LOF is display-only.";
  const auto &route = path(r.config.path);
  meta["path_model"] = QJsonObject{
      {"name", route.name}, {"family", route.family.isEmpty() ? "conventional" : route.family},
      {"before_modes", route.before}, {"after_modes", route.after},
      {"before_via", route.beforeVia}, {"after_via", route.afterVia},
      {"corner_coefficient", route.family == "corner" ? "unit timing approximation; amplitude and phase not calibrated" : "not applicable"}};
  writeJson(folder + "/selection.json", meta);
  writeJson(folder + "/settings.json", {{"schema", 4},
                                        {"input_fingerprint", fingerprint},
                                        {"dataset_name", datasetName},
                                        {"selection", r.config.json()},
                                        {"display", display}});
}
void Engine::exportInspection(const QString &folder, const Inspection &e) {
  QDir().mkpath(folder);
  QFile file(folder + "/trace.csv");
  require(file.open(QIODevice::WriteOnly), "Cannot export trace");
  QTextStream f(&file);
  f.setRealNumberPrecision(17);
  f << "time_us,raw_m_s,filtered_real_m_s,filtered_imag_m_s,arrival_keep\n";
  for (int t = 0; t < Samples; ++t) {
    size_t k = size_t(e.receiver) * Samples + t;
    f << t * .01 << "," << e.raw[k] << "," << e.filtered[k].real() << ","
      << e.filtered[k].imag() << ","
      << keep(t * Dt, e.source, e.receiver, e.config, e.screening.get())
      << "\n";
  }
  writeJson(folder + "/trace_settings.json",
            {{"input_fingerprint", fingerprint},
             {"dataset_name", datasetName},
             {"selection", e.config.json()},
             {"source", e.source},
             {"receiver", e.receiver},
             {"signal", e.signal},
             {"has_point", e.hasPoint},
             {"point_x_mm", e.point.x},
             {"point_depth_mm", e.point.z},
             {"competing_path", e.competitor}});
}
QJsonObject Engine::verifyInputs(std::function<void(QString)> progress) {
  int count = 0;
  if (bundled) {
    require(shaFile(data + "/dataset.json") == manifestHash,
            "Dataset geometry changed; reopen the folder");
    ++count;
  }
  for (auto &[p, h] : inputHashes) {
    if (cancel)
      throw Cancelled();
    if (progress)
      progress("Verifying " + p);
    require(shaFile(data + "/" + p) == h, "Input checksum changed: " + p);
    ++count;
  }
  for (auto v : arrivalAudit["model_inputs"].toArray()) {
    auto o = v.toObject();
    QString p = o["path"].toString();
    if (cancel)
      throw Cancelled();
    if (progress)
      progress("Verifying " + p);
    require(shaFile(data + "/" + p) == o["sha256"].toString(),
            "Saved model changed: " + p);
    ++count;
  }
  return {
      {"status", "passed"},
      {"checksums", count},
      {"input_fingerprint", fingerprint},
      {"created_utc", QDateTime::currentDateTimeUtc().toString(Qt::ISODate)}};
}
void Engine::setSignalProvider(SignalProvider provider, int half, bool rest) {
  require(half >= 0 && half < Samples / 2, "Invalid signal support half-width");
  signalProvider = std::move(provider);
  signalHalf = signalProvider ? half : Half;
  initialRest = signalProvider && rest;
  filteredCache.clear();
  screeningCache.clear();
}
void Engine::clearMemory() {
  rawCache.clear();
  filteredCache.clear();
  curveCache.clear();
  screeningCache.clear();
  geo = Geometry{};
}
} // namespace saft
