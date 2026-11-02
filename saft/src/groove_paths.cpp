#include "saft.hpp"
#include <algorithm>
#include <cmath>
#include <limits>

namespace saft {
namespace {
constexpr double eps = 1e-10; // metres, matches exposed-surface checks
Point sub(Point a, Point b) { return {a.x-b.x, a.z-b.z}; }
Point add(Point a, Point b) { return {a.x+b.x, a.z+b.z}; }
Point mul(Point a, double s) { return {a.x*s, a.z*s}; }
double dot(Point a, Point b) { return a.x*b.x+a.z*b.z; }
double cross(Point a, Point b) { return a.x*b.z-a.z*b.x; }
double length(Point a) { return std::hypot(a.x,a.z); }
Ray missing() {
  Ray r; r.valid=false; r.time=std::numeric_limits<double>::quiet_NaN();
  r.wall=0; return r;
}
struct Face { Point a,b; };
Point mirror(Point p, Face f) {
  const auto d=sub(f.b,f.a);
  const auto foot=add(f.a,mul(d,dot(sub(p,f.a),d)/dot(d,d)));
  return sub(mul(foot,2),p);
}
// Image-source construction for a prescribed ordered set of finite faces.
// Reject endpoint hits: a specular reflection is not a corner diffraction.
bool specular(Point a, Point b, const std::vector<Face> &faces,
              std::vector<Point> &vertices) {
  vertices={a};
  for (size_t i=0;i<faces.size();++i) {
    auto unfolded=b;
    for (size_t j=faces.size();j>i;--j) unfolded=mirror(unfolded,faces[j-1]);
    auto direction=sub(unfolded,vertices.back());
    auto edge=sub(faces[i].b,faces[i].a);
    auto offset=sub(faces[i].a,vertices.back());
    double det=cross(direction,edge);
    if (std::abs(det)<=1e-12*length(direction)*length(edge)) return false;
    double t=cross(offset,edge)/det, u=cross(offset,direction)/det;
    double margin=eps/length(edge);
    if (!(t>0 && t<1 && u>margin && u<1-margin)) return false;
    vertices.push_back(add(faces[i].a,mul(edge,u)));
  }
  vertices.push_back(b);
  return true;
}
struct Feature { QString key,label,family; int index; };
std::vector<Feature> features(const Engine &e) {
  std::vector<Feature> result;
  if (e.datasetKind!="partial") return result;
  for (int i=1;i+1<int(e.surfaceProfile.size());++i) {
    const auto p=e.surfaceProfile[i];
    const QString side=p.x<0?"left":"right";
    const bool floor=std::abs(p.z-std::max_element(e.surfaceProfile.begin(),e.surfaceProfile.end(),
                                [](Point a,Point b){return a.z<b.z;})->z)<eps;
    const QString name=floor?"floor corner":std::abs(p.z)<eps?"lip corner":"bevel knee";
    result.push_back({QString("corner:%1").arg(i),side+" "+name,"corner",i});
  }
  for (int i=0;i+1<int(e.surfaceProfile.size());++i) {
    const auto p=e.surfaceProfile[i],q=e.surfaceProfile[i+1];
    if (q.x-p.x>eps && std::abs(q.z-p.z)>eps)
      result.push_back({QString("wall:%1").arg(i),(p.x<0?QString("left"):QString("right"))+" bevel","wall",i});
  }
  return result;
}
}

Ray Engine::grooveRay(Point a,Point b,const QString &modes,const QString &via,
                     bool outgoing,bool retainPoints) const {
  if (datasetKind!="partial" || !containsSteel(a) || !containsSteel(b)) return missing();
  const auto tokens=via.split(':');
  if (tokens.size()<2 || tokens.size()>3) return missing();
  bool ok=false; int index=tokens[1].toInt(&ok);
  const QString order=tokens.size()==3?tokens[2]:QString();
  if (!ok || (!order.isEmpty() && order!="id-first" && order!="id-last")) return missing();
  if (modes.size()!=(order.isEmpty()?2:3) ||
      (modes.front()!='L' && modes.front()!='S') ||
      std::any_of(modes.begin(),modes.end(),[&](QChar c){return c!=modes.front();})) return missing();
  const char mode=modes.front().toLatin1();
  const double speed=mode=='L'?material.cl:material.cs;
  const Face back{{surfaceProfile.front().x,material.h},
                  {surfaceProfile.back().x,material.h}};
  std::vector<Point> vertices;
  std::vector<std::pair<size_t,Face>> reflecting;
  if (tokens[0]=="wall") {
    if (index<0 || index+1>=int(surfaceProfile.size())) return missing();
    Face face{surfaceProfile[index],surfaceProfile[index+1]};
    if (face.b.x-face.a.x<=eps || std::abs(face.b.z-face.a.z)<=eps) return missing();
    std::vector<Face> faces;
    if (order=="id-first") faces.push_back(back);
    faces.push_back(face);
    if (order=="id-last") faces.push_back(back);
    if (!specular(a,b,faces,vertices)) return missing();
    for (size_t i=0;i<faces.size();++i) reflecting.push_back({i+1,faces[i]});
  } else if (tokens[0]=="corner") {
    if (index<1 || index+1>=int(surfaceProfile.size())) return missing();
    const auto corner=surfaceProfile[index];
    if (order=="id-first") {
      if (!specular(a,corner,{back},vertices)) return missing();
      reflecting.push_back({1,back}); vertices.push_back(b);
    } else if (order=="id-last") {
      std::vector<Point> tail;
      if (!specular(corner,b,{back},tail)) return missing();
      vertices={a};vertices.insert(vertices.end(),tail.begin(),tail.end());
      reflecting.push_back({2,back});
    } else vertices={a,corner,b};
  } else return missing();
  Ray r;
  for (size_t i=1;i<vertices.size();++i) {
    const double leg=length(sub(vertices[i],vertices[i-1]));
    if (leg<=eps || !visibleInSteel(vertices[i-1],vertices[i])) return missing();
    r.time+=leg/speed;
  }
  for (const auto &entry:reflecting) {
    const size_t i=entry.first; const auto face=entry.second;
    auto tangent=mul(sub(face.b,face.a),1/length(sub(face.b,face.a)));
    Point normal{-tangent.z,tangent.x};
    if (std::abs(face.a.z-material.h)<eps && std::abs(face.b.z-material.h)<eps)
      normal=mul(normal,-1);
    auto incident=mul(sub(vertices[i],vertices[i-1]),1/length(sub(vertices[i],vertices[i-1])));
    auto reflected=mul(sub(vertices[i+1],vertices[i]),1/length(sub(vertices[i+1],vertices[i])));
    if (dot(incident,normal)>=-1e-9 || dot(reflected,normal)<=1e-9) return missing();
    const auto expected=sub(incident,mul(normal,2*dot(incident,normal)));
    if (length(sub(expected,reflected))>1e-7) return missing();
    const double angle=std::acos(std::clamp(-dot(incident,normal),0.,1.));
    r.wall*=reflect(angle,mode,material).amplitudes[mode=='L'?0:1];
  }
  auto first=sub(vertices[1],vertices[0]);
  auto last=sub(vertices.back(),vertices[vertices.size()-2]);
  // Generation must leave a horizontal exposed surface into steel;
  // reception must approach it from steel. Grazing endpoints retain 90 deg.
  if ((!outgoing && first.z < -eps) || (outgoing && last.z > eps)) return missing();
  r.first=std::atan2(first.x,std::abs(first.z));
  r.last=std::atan2(last.x,std::abs(last.z));
  // Corner encounters contribute travel time, not a claimed elastic wedge
  // coefficient. Only actual specular events contribute to r.wall.
  if (retainPoints) r.points=std::move(vertices);
  return r;
}

void Engine::addGroovePaths() {
  for (const auto &feature:features(*this))
    for (QChar mode:QString("SL"))
      for (bool outgoing:{false,true})
        for (int extraBackwall:{0,1,2})
          for (int otherLegs:{1,2}) {
            const bool idFirst = extraBackwall == (outgoing ? 2 : 1);
            const QString via=feature.key+(extraBackwall?(idFirst?":id-first":":id-last"):"");
            const QString branch(extraBackwall?3:2,mode), other(otherLegs,mode);
            QString branchName=QString(mode)+"-"+feature.label+"-"+QString(mode);
            if (extraBackwall)
              branchName=idFirst?QString(mode)+"-ID-"+branchName:branchName+"-ID-"+QString(mode);
            const QString id=QString("exp_%1_%2_%3_%4_%5_%6")
                .arg(feature.family).arg(feature.index).arg(mode)
                .arg(outgoing?"rx":"tx").arg(extraBackwall==1?"id":extraBackwall==2?"pixel-id":"direct").arg(otherLegs);
            paths.push_back({id,outgoing?other:branch,outgoing?branch:other,
                "Exp: "+(outgoing?other+" | "+branchName:branchName+" | "+other),
                outgoing?QString():via,outgoing?via:QString(),feature.family});
          }
}
QJsonObject Engine::grooveCatalogue() const {
  QJsonArray fs,ps;
  for (const auto &f:features(*this)) {
    QJsonObject value{{"key",f.key},{"label",f.label},{"family",f.family},{"index",f.index}};
    QJsonArray coords;
    const auto a=surfaceProfile[f.index];coords.append(QJsonArray{a.x,a.z});
    if (f.family=="wall") {const auto b=surfaceProfile[f.index+1];coords.append(QJsonArray{b.x,b.z});}
    value["points_m"]=coords;fs.append(value);
  }
  for (const auto &p:paths) ps.append(QJsonObject{{"id",p.id},{"name",p.name},
      {"before",p.before},{"after",p.after},{"before_via",p.beforeVia},
      {"after_via",p.afterVia},{"family",p.family}});
  return {{"features",fs},{"paths",ps},
      {"scope","Experimental same-mode L/S paths: one exposed bevel reflection or one fixed groove-corner timing encounter; optional backwall before or after it. Corner diffraction amplitude/phase is not modelled."}};
}
} // namespace saft
