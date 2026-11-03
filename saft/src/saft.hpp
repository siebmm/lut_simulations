#pragma once
#include <QtCore>
#include <array>
#include <atomic>
#include <complex>
#include <functional>
#include <map>
#include <memory>
#include <set>
#include <stdexcept>
#include <vector>

namespace saft {
using Complex = std::complex<double>;
using Vec = std::vector<double>;
using CVec = std::vector<Complex>;
constexpr int Samples = 5001, Receivers = 337, Half = 328;
constexpr double Dt = 1e-8, Pi = 3.14159265358979323846;
struct Cancelled : std::runtime_error {
  Cancelled() : std::runtime_error("Cancelled") {}
};
struct Point {
  double x = 0, z = 0;
};
struct Material {
  double cl = 0, cs = 0, cr = 0, h = 0, width = 0, pulse = 0, active = .06;
  int receiverCount = Receivers;
  int sourceCount = 5, groupCount = 3;
  double dt = Dt, xMin = -.06, xMax = .06, zMin = 0, zMax = .0254;
  bool general = false;
  QStringList signals = {"lof", "reference", "difference"};
};
struct Receiver {
  int id = 0, group = 0, groupIndex = 0;
  double x = 0, z = 0;
  Point normal{0, 1};
  QString label;
};
struct Source {
  QString id, name;
  double x = 0, t0 = 0, z = 0;
  Point normal{0, 1};
};
struct Path {
  QString id, before, after, name;
  // Empty via strings preserve the original planar branch implementation.
  QString beforeVia, afterVia, family;
};
struct Config {
  QString path = "1-1_LL", component = "vy", aperture = "all",
          combination = "rms", policy = "minimum", processing = "weighted",
          repairs = "case", screening = "off";
  QStringList cases = {"lof"};
  std::vector<int> sources, receivers;
  // Missing source key inherits the shared list; an explicit empty list
  // admits no receivers for that source. Always intersect with shared IDs.
  std::map<int, std::vector<int>> sourceReceivers;
  std::vector<bool> groups;
  std::array<bool, 6> mutes = {false, false, false, false, false, false};
  std::array<bool, 4> reflectedScreen{}; // LL, LS/SL, SS, repeated SSSS
  double low = .5, high = 5, x0 = -15, x1 = 15, z0 = .25, z1 = 25.25,
         spacing = .5, gate = .2, wall = .05, offset = 2, muteWidth = .6,
         cfPower = 1, screenTailDb = 30, screenMargin = .1;
  int minimumSources = 2, minimumReceivers = 8, threads = 8;
  double fraction = .25, minimumEffective = 8;
  bool clean = false, coherence = false, freeze = false;
  Point anchor{2.5, 20};
  explicit Config(int receiverCount = Receivers, int sourceCount = 5, int groupCount = 3);
  void validate(const Material &) const;
  QJsonObject json() const;
  static Config fromJson(QJsonObject, int receiverCount = Receivers, int sourceCount = 5, int groupCount = 3);
  Vec xs() const;
  Vec zs() const;
};
struct Ray {
  double time = 0, first = 0, last = 0;
  Complex wall = 1;
  std::vector<Point> points;
  bool valid = true;
  Point departure, arrival;
};
struct Reflection {
  std::array<Complex, 2> amplitudes, surface;
};
Reflection reflect(double angle, char mode, const Material &);
Ray ray(Point a, Point b, const QString &modes, const Material &,
        bool outgoing = false, bool points = false, double top = 0);
CVec kernel(double low, double high, double dt = Dt, int half = Half);
void fft(CVec &, bool inverse = false);
void parallel(int count, int threads, const std::function<void(int)> &,
              std::atomic_bool *cancel = nullptr);
QString sha(const QByteArray &);
QString shaFile(const QString &);
QJsonObject readJson(const QString &);
void writeJson(const QString &, const QJsonObject &);
std::vector<int> receiverIds(const QString &, int count = Receivers);
extern const QStringList GroupNames, MuteNames;
extern const QStringList PulseNames, ReflectedScreenNames;
const QStringList &shortlistPaths();

struct Stats {
  double count = 0, eligible = 0, weight = 0, weight2 = 0, power = 0,
         candidates = 0, screened = 0;
  Complex sum = 0;
};
struct Plane {
  Vec amplitude, count, cf, retained, screened;
  std::vector<quint64> receiverUsage;
  quint64 supportedPixels = 0;
  std::vector<Stats> stats;
  Vec effective;
  std::vector<unsigned char> valid;
  CVec means;
};
struct PulseWindow {
  int family = 0;
  double predicted = 0, peakTime = 0, peak = 0, first = 0, last = 0;
  bool limited = false, measured = false;
};
struct Screening {
  QString basis;
  std::vector<std::vector<std::vector<PulseWindow>>> windows;
  explicit Screening(int sources = 5, int receivers = Receivers)
      : windows(sources, std::vector<std::vector<PulseWindow>>(receivers)) {}
  bool rejects(double time, int source, int receiver) const;
  QJsonObject json() const;
};
PulseWindow measurePulse(const CVec &, int receiver, int family,
                         double predicted, double reach, const Config &,
                         int samples = Samples, double dt = Dt, double start = 0);
QJsonObject screeningSelfCheck();
struct Result {
  int sourceCount = 5;
  Config config;
  Vec x, z;
  std::vector<Point> points;
  std::map<QString, std::array<Plane, 3>> planes;
  std::map<QString, std::shared_ptr<Screening>> screening;
  double seconds = 0, loadSeconds = 0, filterSeconds = 0, geometrySeconds = 0,
         sumSeconds = 0;
  QJsonObject timing() const;
};
struct Support {
  double gain = 0, angle = 0, count = 0, effective = 0;
  bool valid = false;
  std::vector<int> ids;
};
struct Inspection {
  Config config;
  int source = 0, receiver = 0;
  QString signal = "lof", competitor = "1-2_SSS";
  bool hasPoint = false;
  Point point;
  Vec raw;
  CVec filtered;
  Vec times, competing;
  std::array<Vec, 7> arrivals;
  std::vector<int> orderedIds;
  std::vector<Point> path, competingPath;
  std::vector<Support> support;
  std::shared_ptr<Screening> screening;
  double seconds = 0;
};
struct Curves {
  std::array<Vec, 2> launch;
  std::array<std::array<Vec, 2>, 2> receive;
  std::array<std::array<double, 2>, 2> norm;
};
class Engine {
public:
  explicit Engine(QString root, QString dataset = "");
  QString root, data, cache, fingerprint, datasetName = "Full weld",
                                          datasetKind = "full", manifestHash;
  bool bundled = false;
  bool general = false, isotropic = false;
  int samples = Samples, filterHalf = Half;
  double dt = Dt, startTime = 0;
  QStringList groupNames = GroupNames, signalCases = {"lof", "reference"};
  QString primaryCase = "lof", referenceCase = "reference";
  struct Boundary { QString id; Point a, b, normal; };
  std::vector<std::vector<Point>> polygons;
  std::vector<Boundary> boundaries;
  QString backwallBoundary, surfaceBoundary;
  QJsonObject description, acquisition;
  struct Recording { QString file; int source; QString signal, component; };
  std::vector<Recording> recordings;
  std::vector<Point> surfaceProfile;
  int receiverCount() const { return int(receivers.size()); }
  int sourceCount() const { return int(sources.size()); }
  Config config(QJsonObject) const;
  void loadGeneral(const QJsonObject &);
  std::shared_ptr<Vec> rawGeneral(int, const QString &, const QString &);
  bool containsPolygon(Point) const;
  bool visibleInPolygon(Point, Point) const;
  Ray polygonRay(Point, Point, const QString &, const QString &, bool) const;
  double sourceAngle(const Ray &, int) const;
  double receiverAngle(const Ray &, int) const;
  Complex receiverResponse(const Ray &, int, const Config &, const Curves &, char) const;
  bool containsSteel(Point) const;
  bool visibleInSteel(Point, Point) const;
  Ray trace(Point, Point, const QString &, bool outgoing = false,
            bool points = false, const QString &via = {}) const;
  Ray grooveRay(Point, Point, const QString &, const QString &via,
                bool outgoing = false, bool points = false) const;
  void addGroovePaths();
  QJsonObject grooveCatalogue() const;
  std::vector<Ray> shearReverberations(Point source, Point receiver) const;
  double sourceDistance(int source, int receiver) const;
  QJsonObject audit, arrivalAudit, geometry;
  Material material;
  std::vector<Receiver> receivers;
  std::vector<Source> sources;
  std::vector<Path> paths;
  std::array<Point, 2> truth;
  std::map<QString, QString> inputHashes;
  std::map<QString, int> loads;
  int csvParses = 0, cacheReads = 0, geometryBuilds = 0, filterBuilds = 0;
  std::atomic_bool cancel{false};
  std::shared_ptr<Result>
  reconstruct(const Config &, std::function<void(QString)> progress = {},
              const std::vector<Point> *points = nullptr);
  std::shared_ptr<Inspection> inspect(const Config &, int source, int receiver,
                                      QString signal, QString competitor,
                                      bool hasPoint, Point);
  std::shared_ptr<Vec> raw(int si, QString signal, QString component);
  std::shared_ptr<CVec> filtered(int si, QString signal, const Config &);
  std::shared_ptr<Curves> curves(const Config &);
  std::vector<Support> support(const Config &, Point, QString signal);
  double keep(double time, int source, int receiver, const Config &,
              const Screening * = nullptr) const;
  std::shared_ptr<Screening> screeningFor(const Config &, QString signal);
  const Path &path(const QString &) const;
  void exportResult(const QString &, const Result &, QJsonObject display = {});
  void exportInspection(const QString &, const Inspection &);
  QJsonObject verifyInputs(std::function<void(QString)> progress);
  QJsonObject provenance() const;
  void clearMemory();

  // Research injection point: a supplied analytic (or explicitly real) signal
  // bank can be evaluated with exactly the production geometry and sums.
  // The desktop never installs a provider; its original filtering is unchanged.
  using SignalProvider = std::function<std::shared_ptr<CVec>(int, QString, const Config &)>;
  void setSignalProvider(SignalProvider, int supportHalfSamples = Half,
                         bool knownInitialRest = false);
  int signalHalfSamples() const { return signalHalf; }
  bool signalSampleSupported(int j) const {
    return j >= 0 && j + 1 < samples &&
           (initialRest || j - signalHalf >= 0) && j + 1 + signalHalf < samples;
  }

private:
  SignalProvider signalProvider;
  int signalHalf = Half;
  bool initialRest = false;
  std::map<QString, std::shared_ptr<Vec>> rawCache;
  std::map<QString, std::shared_ptr<CVec>> filteredCache;
  std::map<QString, std::shared_ptr<Curves>> curveCache;
  std::map<QString, std::shared_ptr<Screening>> screeningCache;
  std::map<QString, std::vector<std::pair<double, double>>> repairRanges;
  std::vector<std::vector<std::array<double, 6>>> nuisance;
  std::vector<Vec> nuisanceSL;
  std::vector<std::vector<Vec>> nuisanceSSSS;
  struct Geometry {
    QString key;
    std::vector<Ray> receiving;
    std::vector<std::vector<Ray>> incoming;
  } geo;
  void prepareGeometry(const Config &, const std::vector<Point> &);
  bool supportClean(double time, int si, QString signal, int group,
                    const Config &) const;
};
} // namespace saft
