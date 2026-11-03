#include "saft.hpp"
#include <algorithm>
#include <cmath>

namespace saft {
namespace {
constexpr double eps=1e-10;
Point sub(Point a,Point b){return {a.x-b.x,a.z-b.z};}
Point add(Point a,Point b){return {a.x+b.x,a.z+b.z};}
Point mul(Point a,double s){return {a.x*s,a.z*s};}
double dot(Point a,Point b){return a.x*b.x+a.z*b.z;}
double cross(Point a,Point b){return a.x*b.z-a.z*b.x;}
double length(Point a){return std::hypot(a.x,a.z);}
bool on(Point p,Point a,Point b) {
  auto v=sub(b,a),w=sub(p,a);
  return std::abs(cross(v,w))<=eps*length(v)&&dot(w,sub(p,b))<=eps*eps;
}
// 0 outside, 1 interior, 2 boundary. Loop orientation is irrelevant here.
int inside(Point p,const std::vector<Point> &poly) {
  bool found=false;
  for(size_t i=0,j=poly.size()-1;i<poly.size();j=i++) {
    auto a=poly[j],b=poly[i];
    if(on(p,a,b)) return 2;
    if((a.z>p.z)!=(b.z>p.z) && p.x<(b.x-a.x)*(p.z-a.z)/(b.z-a.z)+a.x) found=!found;
  }
  return found?1:0;
}
Ray missing() {Ray r;r.valid=false;r.time=std::nan("");r.wall=0;return r;}
}

bool Engine::containsPolygon(Point p) const {
  if(!std::isfinite(p.x)||!std::isfinite(p.z)||polygons.empty()||!inside(p,polygons[0])) return false;
  for(size_t i=1;i<polygons.size();++i) if(inside(p,polygons[i])==1) return false;
  return true;
}
bool Engine::visibleInPolygon(Point a,Point b) const {
  if(!containsPolygon(a)||!containsPolygon(b)) return false;
  Point direction=sub(b,a);double distance=length(direction);
  if(distance<eps) return true;
  Vec cuts{0.,1.};
  for(const auto &edge:boundaries) {
    Point v=sub(edge.b,edge.a),w=sub(edge.a,a);double det=cross(direction,v);
    if(std::abs(det)>1e-20) {
      double t=cross(w,v)/det,u=cross(w,direction)/det;
      if(t>0&&t<1&&u>=-1e-10&&u<=1+1e-10) cuts.push_back(t);
    } else if(std::abs(cross(direction,w))<eps*distance) {
      for(auto p:{edge.a,edge.b}) {
        double t=dot(sub(p,a),direction)/dot(direction,direction);
        if(t>0&&t<1) cuts.push_back(t);
      }
    }
  }
  std::sort(cuts.begin(),cuts.end());
  for(size_t i=1;i<cuts.size();++i)
    if(!containsPolygon(add(a,mul(direction,(cuts[i-1]+cuts[i])*.5)))) return false;
  return true;
}

Ray Engine::polygonRay(Point a,Point b,const QString &modes,const QString &via,bool retain) const {
  if(!containsPolygon(a)||!containsPolygon(b)) return missing();
  QStringList names=via.isEmpty()?QStringList{}:via.split('|');
  if(modes.size()!=names.size()+1||modes.size()>5||!QRegularExpression("^[LS]+$").match(modes).hasMatch())
    throw std::runtime_error("General 2D ray branches require explicit boundaries for every reflection");
  std::vector<const Boundary*> faces;
  for(const auto &name:names) {
    auto it=std::find_if(boundaries.begin(),boundaries.end(),[&](const Boundary &f){return f.id==name;});
    if(it==boundaries.end()) throw std::runtime_error(("Unknown reflecting boundary: "+name).toStdString());
    faces.push_back(&*it);
  }
  Vec speeds;for(auto mode:modes) speeds.push_back(mode=='L'?material.cl:material.cs);
  std::vector<Point> vertices{a};
  for(auto face:faces) vertices.push_back(mul(add(face->a,face->b),.5));
  vertices.push_back(b);
  Vec positions(faces.size(),.5);
  // Fermat travel time is a convex sum of segment lengths over the ordered
  // finite faces. Coordinate minimization enforces tangential slowness (Snell)
  // for both same-mode and converted reflections. Corner optima are rejected:
  // they are diffraction hypotheses, not specular reflections.
  bool converged=faces.empty();
  for(int iteration=0;iteration<256&&!faces.empty();++iteration) {
    double change=0;
    for(size_t i=0;i<faces.size();++i) {
      Point origin=faces[i]->a,edge=sub(faces[i]->b,origin);
      auto derivative=[&](double u) {
        Point p=add(origin,mul(edge,u)),incoming=sub(p,vertices[i]),outgoing=sub(p,vertices[i+2]);
        return dot(edge,incoming)/(std::max(length(incoming),1e-30)*speeds[i])+
               dot(edge,outgoing)/(std::max(length(outgoing),1e-30)*speeds[i+1]);
      };
      double lo=0,hi=1;
      if(derivative(0)>=0) hi=0;
      else if(derivative(1)<=0) lo=1;
      else for(int k=0;k<48;++k) {double mid=(lo+hi)*.5;if(derivative(mid)<0)lo=mid;else hi=mid;}
      double u=(lo+hi)*.5;change=std::max(change,std::abs(u-positions[i]));positions[i]=u;
      vertices[i+1]=add(origin,mul(edge,u));
    }
    if(change<1e-12) {converged=true;break;}
  }
  if(!converged) return missing();
  Ray result;
  for(size_t i=1;i<vertices.size();++i) {
    double distance=length(sub(vertices[i],vertices[i-1]));
    if(distance<=eps||!visibleInPolygon(vertices[i-1],vertices[i])) return missing();
    result.time+=distance/speeds[i-1];
  }
  for(size_t i=0;i<faces.size();++i) {
    if(positions[i]<=1e-9||positions[i]>=1-1e-9) return missing();
    auto incoming=mul(sub(vertices[i+1],vertices[i]),1/length(sub(vertices[i+1],vertices[i])));
    auto outgoing=mul(sub(vertices[i+2],vertices[i+1]),1/length(sub(vertices[i+2],vertices[i+1])));
    Point normal=faces[i]->normal,tangent{normal.z,-normal.x};
    double ci=-dot(incoming,normal),co=dot(outgoing,normal);
    if(ci<=1e-10||co<=1e-10) return missing();
    double pi=dot(incoming,tangent)/speeds[i],po=dot(outgoing,tangent)/speeds[i+1];
    if(std::abs(pi-po)>1e-9/std::max(speeds[i],speeds[i+1])) return missing();
    auto response=reflect(std::atan2(dot(incoming,tangent),ci),modes[int(i)].toLatin1(),material);
    result.wall*=response.amplitudes[modes[int(i+1)]=='L'?0:1]*std::sqrt(speeds[i+1]/speeds[i]*co/ci);
  }
  result.departure=mul(sub(vertices[1],a),1/length(sub(vertices[1],a)));
  result.arrival=mul(sub(b,vertices[vertices.size()-2]),1/length(sub(b,vertices[vertices.size()-2])));
  result.first=std::atan2(result.departure.x,result.departure.z);
  result.last=std::atan2(result.arrival.x,-result.arrival.z);
  if(retain) result.points=vertices;
  return result;
}
} // namespace saft
