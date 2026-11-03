#include "saft.hpp"
#include <algorithm>
#include <cmath>
#include <limits>

namespace saft {
namespace {
constexpr double tolerance = 1e-10;
double surfaceDepth(const std::vector<Point> &profile, double x) {
  double result = std::numeric_limits<double>::infinity();
  for (const auto &p : profile)
    if (std::abs(p.x - x) <= 1e-12)
      result = std::min(result, p.z);
  if (std::isfinite(result))
    return result; // exposed vertical-wall boundary
  for (size_t i = 1; i < profile.size(); ++i) {
    const auto a = profile[i - 1], b = profile[i];
    if (b.x > a.x && x > a.x && x < b.x)
      return a.z + (b.z - a.z) * (x - a.x) / (b.x - a.x);
  }
  return result;
}
Ray absent() {
  Ray r;
  r.valid = false;
  r.time = std::numeric_limits<double>::quiet_NaN();
  r.wall = 0;
  return r;
}
} // namespace
bool Engine::containsSteel(Point p) const {
  if (general) return containsPolygon(p);
  if (datasetKind == "full")
    return p.z >= 0 && p.z <= material.h;
  return p.z >= surfaceDepth(surfaceProfile, p.x) - tolerance &&
         p.z <= material.h + tolerance;
}
bool Engine::visibleInSteel(Point a, Point b) const {
  if (general) return visibleInPolygon(a,b);
  if (!containsSteel(a) || !containsSteel(b))
    return false;
  if (datasetKind == "full")
    return true;
  const double dx = b.x - a.x;
  if (std::abs(dx) < 1e-15)
    return true;
  // Both the ray and each exposed-surface segment are linear. Endpoints and
  // every surface breakpoint therefore suffice, including both vertical-wall
  // ends.
  for (const auto &knot : surfaceProfile) {
    const double t = (knot.x - a.x) / dx;
    if (t > 1e-12 && t < 1 - 1e-12 &&
        a.z + t * (b.z - a.z) < knot.z - tolerance)
      return false;
  }
  return true;
}
double Engine::sourceDistance(int s, int r) const {
  const double dx = receivers[r].x - sources[s].x;
  return datasetKind == "full" ? std::abs(dx)
                               : std::hypot(dx, receivers[r].z - sources[s].z);
}
Ray Engine::trace(Point a, Point b, const QString &modes, bool outgoing,
                  bool points, const QString &via) const {
  if (general) return polygonRay(a,b,modes,via,points);
  if (!via.isEmpty())
    return grooveRay(a, b, modes, via, outgoing, points);
  if (modes == "SSSS") {
    auto returns = shearReverberations(a, b);
    if (returns.empty())
      return absent();
    auto result = *std::min_element(
        returns.begin(), returns.end(),
        [](const Ray &a, const Ray &b) { return a.time < b.time; });
    if (!points)
      result.points.clear();
    return result;
  }
  if (datasetKind == "full")
    return ray(a, b, modes, material, outgoing, points);
  if (!containsSteel(a) || !containsSteel(b))
    return absent();
  const int legs = modes.size();
  if (legs == 1 && std::abs(a.z - b.z) < 1e-14) {
    if (!visibleInSteel(a, b) || std::abs(a.x - b.x) < 1e-15)
      return absent();
    Ray result;
    result.time =
        std::abs(b.x - a.x) / (modes[0] == 'L' ? material.cl : material.cs);
    result.first = result.last = (b.x > a.x ? 1 : -1) * Pi / 2;
    if (points)
      result.points = {a, b};
    return result;
  }
  Ray best = absent();
  std::vector<double> tops{0};
  if (legs == 3) {
    for (size_t i = 1; i < surfaceProfile.size(); ++i) {
      auto p = surfaceProfile[i - 1], q = surfaceProfile[i];
      if (q.x > p.x && std::abs(p.z - q.z) < 1e-12 &&
          std::none_of(tops.begin(), tops.end(),
                       [&](double z) { return std::abs(z - p.z) < 1e-12; }))
        tops.push_back(p.z);
    }
  }
  for (double top : tops) {
    if (legs == 3 && (outgoing ? a.z : b.z) <= top + 1e-14)
      continue;
    auto candidate = ray(a, b, modes, material, outgoing, true, top);
    bool valid = true;
    if (legs == 3) {
      const auto reflection = candidate.points[outgoing ? 1 : 2];
      valid = std::abs(surfaceDepth(surfaceProfile, reflection.x) - top) <=
              tolerance;
      // A flat-face reflection must lie on a horizontal exposed segment, not
      // on a sloped groove wall or beyond the edge of a finite weld surface.
      bool onFlat = false;
      for (size_t i = 1; i < surfaceProfile.size(); ++i) {
        const auto p = surfaceProfile[i - 1], q = surfaceProfile[i];
        onFlat |= q.x > p.x && std::abs(p.z - top) < 1e-12 &&
                  std::abs(q.z - top) < 1e-12 &&
                  reflection.x >= p.x - tolerance &&
                  reflection.x <= q.x + tolerance;
      }
      valid &= onFlat;
    }
    for (size_t i = 1; valid && i < candidate.points.size(); ++i)
      valid = visibleInSteel(candidate.points[i - 1], candidate.points[i]);
    if (valid && (!best.valid || candidate.time < best.time))
      best = std::move(candidate);
  }
  if (!points)
    best.points.clear();
  return best;
}

std::vector<Ray> Engine::shearReverberations(Point a, Point b) const {
  // This is a specimen-only source -> receiver return, with no image pixel.
  // Retain every exposed flat-top branch: unequal source/receiver depths can
  // admit distinct OD and fill-top returns. Never put a reflector in the air.
  if (!containsSteel(a) || !containsSteel(b))
    return {};
  std::vector<double> tops{0};
  if (datasetKind == "partial")
    for (size_t i = 1; i < surfaceProfile.size(); ++i) {
      auto p = surfaceProfile[i - 1], q = surfaceProfile[i];
      if (q.x > p.x && std::abs(p.z - q.z) < 1e-12 &&
          std::none_of(tops.begin(), tops.end(),
                       [&](double z) { return std::abs(z - p.z) < 1e-12; }))
        tops.push_back(p.z);
    }
  std::vector<Ray> result;
  for (double top : tops) {
    auto candidate = ray(a, b, "SSSS", material, false, true, top);
    if (datasetKind == "partial") {
      const auto bounce = candidate.points[2];
      bool onFlat = false;
      for (size_t i = 1; i < surfaceProfile.size(); ++i) {
        auto p = surfaceProfile[i - 1], q = surfaceProfile[i];
        onFlat |= q.x > p.x && std::abs(p.z - top) < 1e-12 &&
                  std::abs(q.z - top) < 1e-12 && bounce.x >= p.x - tolerance &&
                  bounce.x <= q.x + tolerance;
      }
      if (!onFlat ||
          std::abs(surfaceDepth(surfaceProfile, bounce.x) - top) > tolerance)
        continue;
      bool visible = true;
      for (size_t i = 1; i < candidate.points.size(); ++i)
        visible &= visibleInSteel(candidate.points[i - 1], candidate.points[i]);
      if (!visible)
        continue;
    }
    result.push_back(std::move(candidate));
  }
  return result;
}
} // namespace saft
