#include "saft.hpp"
#include <algorithm>
#include <cmath>

namespace saft {
// Times are seconds internally. The search identifies a packet near a known
// direct-arrival guide; it is not a wave-mode decomposition or defect detector.
PulseWindow measurePulse(const CVec &trace, int receiver, int family,
                         double predicted, double reach, const Config &c, int samples, double dt, double start) {
  PulseWindow w;
  w.family = family;
  w.predicted = w.peakTime = predicted;
  const double search = .25e-6, margin = c.screenMargin * 1e-6;
  w.first = predicted - reach - margin;
  w.last = predicted + reach + margin;
  if (!std::isfinite(predicted) || predicted+search < start || predicted-search > start+(samples-1)*dt) {
    w.limited=true;
    return w;
  }
  int lo = std::max(0, int(std::ceil((predicted - search - start) / dt)));
  int hi = std::min(samples - 1, int(std::floor((predicted + search - start) / dt)));
  if (lo > hi) {
    w.limited = true;
    return w;
  }
  auto envelope = [&](int j) {
    return std::abs(trace[size_t(receiver) * samples + j]);
  };
  int seed = lo;
  for (int j = lo; j <= hi; ++j)
    if (envelope(j) > envelope(seed))
      seed = j;
  w.peak = envelope(seed);
  w.peakTime = start + seed * dt;
  // An exactly zero packet cannot contaminate a recording.
  if (w.peak == 0) {
    w.first = w.last = predicted;
    w.measured = true;
    return w;
  }
  w.measured = true;
  w.limited = seed == lo || seed == hi;
  const double threshold = w.peak * std::pow(10., -c.screenTailDb / 20.);
  const int maximum = int(std::ceil(reach / dt)), quietsamples = std::max(1, int(std::ceil(.15e-6/dt))); // 0.15 us
  auto edge = [&](int direction) {
    int lastAbove = seed, quiet = 0;
    for (int d = 1; d <= maximum; ++d) {
      int j = seed + direction * d;
      if (j < 0 || j >= samples) {
        w.limited = true;
        return seed + direction * maximum;
      }
      if (envelope(j) >= threshold) {
        lastAbove = j;
        quiet = 0;
      } else if (++quiet >= quietsamples) {
        return lastAbove;
      }
    }
    w.limited = true;
    return seed + direction * maximum;
  };
  w.first = std::min(predicted, start + edge(-1) * dt) - margin;
  w.last = std::max(predicted, start + edge(1) * dt) + margin;
  return w;
}

bool Screening::rejects(double time, int source, int receiver) const {
  for (const auto &w : windows[source][receiver])
    if ((!w.measured || w.peak > 0) && time >= w.first && time <= w.last)
      return true;
  return false;
}

QJsonObject Screening::json() const {
  QJsonArray rows;
  for (int s = 0; s < int(windows.size()); ++s)
    for (int r = 0; r < int(windows[s].size()); ++r)
      for (const auto &w : windows[s][r])
        rows.append(QJsonObject{{"source", s},
                                {"receiver", r},
                                {"family", PulseNames[w.family]},
                                {"predicted_us", w.predicted * 1e6},
                                {"peak_time_us", w.peakTime * 1e6},
                                {"peak_m_s", w.peak},
                                {"start_us", w.first * 1e6},
                                {"end_us", w.last * 1e6},
                                {"measured", w.measured},
                                {"search_or_extent_limited", w.limited}});
  return {{"basis", basis},
          {"windows", rows},
          {"search_halfwidth_us", .25},
          {"quiet_duration_us", .15},
          {"interpretation",
           "Each case uses its own pulse measurement; a difference uses its "
           "declared minuend. Limited windows do not certify uncontaminated "
           "samples outside the guard."}};
}

std::shared_ptr<Screening> Engine::screeningFor(const Config &c,
                                                QString signal) {
  if (c.screening == "off" &&
      std::none_of(c.reflectedScreen.begin(), c.reflectedScreen.end(),
                   [](bool on) { return on; }))
    return {};
  QString basis = signal == "difference" ? primaryCase : signal;
  QJsonArray selected;
  for (int s : c.sources)
    selected.append(s);
  QJsonArray reflected;
  for (bool on : c.reflectedScreen)
    reflected.append(on);
  QString key =
      QString::fromUtf8(QJsonDocument(QJsonObject{{"basis", basis},
                                                  {"component", c.component},
                                                  {"low", c.low},
                                                  {"high", c.high},
                                                  {"mode", c.screening},
                                                  {"reflections", reflected},
                                                  {"tail", c.screenTailDb},
                                                  {"margin", c.screenMargin},
                                                  {"sources", selected}})
                            .toJson(QJsonDocument::Compact));
  auto found = screeningCache.find(key);
  if (found != screeningCache.end())
    return found->second;
  auto result = std::make_shared<Screening>(sourceCount(),receiverCount());
  result->basis = basis;
  // Choose nuisance families, independently of the desired imaging route.
  // Family order follows MuteNames: Rayleigh, direct L, direct S.
  const std::array<bool, 3> families{c.screening.contains('R'),
                                     c.screening.contains('L'),
                                     c.screening.contains('S')};
  for (int s : c.sources) {
    auto trace = filtered(s, basis, c);
    for (int r = 0; r < receiverCount(); ++r) {
      if (cancel)
        throw Cancelled();
      for (int family = 0; family < 3; ++family) {
        if (!families[family] || !std::isfinite(nuisance[s][r][family]))
          continue;
        double speed = family == 0   ? material.cr
                       : family == 1 ? material.cl
                                     : material.cs;
        // Bound the packet search by the FIR support plus source footprint
        // travel time and three saved pulse widths. This is a finite search
        // limit, not a claim that all clutter ends at this time.
        double reach = signalHalfSamples() * dt + material.width / speed + 3 * material.pulse;
        result->windows[s][r].push_back(measurePulse(
            *trace, r, family, nuisance[s][r][family] * 1e-6, reach, c, samples, dt, startTime));
      }
      auto append = [&](int family, double predictedUs, double speed) {
        if (!std::isfinite(predictedUs))
          return;
        const double reach =
            signalHalfSamples() * dt + material.width / speed + 3 * material.pulse;
        result->windows[s][r].push_back(
            measurePulse(*trace, r, family, predictedUs * 1e-6, reach, c, samples, dt, startTime));
      };
      if (c.reflectedScreen[0])
        append(3, nuisance[s][r][3], material.cl);
      if (c.reflectedScreen[1]) {
        append(4, nuisance[s][r][4], material.cs);
        if (!std::isfinite(nuisance[s][r][4]) ||
            std::abs(nuisanceSL[s][r] - nuisance[s][r][4]) > 1e-8)
          append(6, nuisanceSL[s][r], material.cs);
      }
      if (c.reflectedScreen[2])
        append(5, nuisance[s][r][5], material.cs);
      if (c.reflectedScreen[3])
        for (double centre : nuisanceSSSS[s][r])
          append(7, centre, material.cs);
    }
  }
  // Keep bounded: current and recent parameter combinations only.
  if (screeningCache.size() >= 12)
    screeningCache.clear();
  screeningCache[key] = result;
  return result;
}

QJsonObject screeningSelfCheck() {
  auto check = [](bool condition, const char *message) {
    if (!condition)
      throw std::runtime_error(message);
  };
  CVec trace(Samples);
  const double centre = 20e-6, sigma = .15e-6;
  for (int j = 0; j < Samples; ++j) {
    double t = j * Dt;
    trace[j] = std::exp(-.5 * std::pow((t - centre) / sigma, 2)) *
               std::exp(Complex(0, 2 * Pi * 2e6 * t));
  }
  Config c;
  c.screenMargin = 0;
  auto w = measurePulse(trace, 0, 1, centre - .1e-6, 3.5e-6, c);
  double radius = sigma * std::sqrt(2 * std::log(10.) * c.screenTailDb / 20.);
  check(std::abs(w.first - (centre - radius)) <= Dt &&
            std::abs(w.last - (centre + radius)) <= Dt,
        "Gaussian packet edges disagree with analytic threshold");
  c.screenMargin = .12;
  auto margin = measurePulse(trace, 0, 1, centre - .1e-6, 3.5e-6, c);
  check(std::abs(margin.first - w.first + .12e-6) < 1e-12 &&
            std::abs(margin.last - w.last - .12e-6) < 1e-12,
        "Timing margin is not applied on both sides");
  for (auto &v : trace)
    v *= 1e-9;
  auto scaled = measurePulse(trace, 0, 1, centre - .1e-6, 3.5e-6, c);
  check(std::abs(scaled.first - margin.first) < 1e-12 &&
            std::abs(scaled.last - margin.last) < 1e-12,
        "Relative packet widths changed with amplitude scale");
  c.screenTailDb = 40;
  auto lower = measurePulse(trace, 0, 1, centre - .1e-6, 3.5e-6, c);
  check(lower.first <= scaled.first && lower.last >= scaled.last,
        "Lower tail threshold shortened a Gaussian window");
  auto limited = measurePulse(trace, 0, 1, centre, .05e-6, c);
  check(limited.limited, "Finite search limit was not reported");
  auto outside = measurePulse(trace, 0, 1, 80e-6, 3.5e-6, c);
  check(!outside.measured && outside.limited,
        "Out-of-record guide was not flagged");
  Screening screen;
  screen.windows[0][0].push_back(w);
  check(screen.rejects(centre, 0, 0) && !screen.rejects(centre + 2e-6, 0, 0),
        "Screen interval membership failed");
  std::fill(trace.begin(), trace.end(), Complex(0));
  screen.windows[0][0] = {measurePulse(trace, 0, 1, centre, 3.5e-6, c)};
  check(!screen.rejects(centre, 0, 0),
        "Zero trace was treated as contamination");
  return {{"status", "passed"}, {"checks", 8}};
}

} // namespace saft
