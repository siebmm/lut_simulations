#include "saft.hpp"
#include <algorithm>
#include <charconv>
#include <cmath>
#include <cstring>
#include <numeric>

namespace saft {
namespace {
void check(bool ok, const QString &message) {
  if (!ok) throw std::runtime_error(message.toStdString());
}
double number(QJsonObject o, const QString &key) {
  check(o[key].isDouble() && std::isfinite(o[key].toDouble()), "Expected finite number: " + key);
  return o[key].toDouble();
}
int integer(QJsonObject o, const QString &key, int minimum) {
  double value = number(o,key);
  check(value >= minimum && value <= 100000000 && value == std::floor(value), "Invalid integer: " + key);
  return int(value);
}
QString id(QJsonValue v) {
  QString s=v.toString();
  check(QRegularExpression("^[A-Za-z][A-Za-z0-9_-]*$").match(s).hasMatch(), "Invalid identifier: " + s);
  return s;
}
Point point(QJsonValue value) {
  auto a=value.toArray();
  check(a.size()==2 && a[0].isDouble() && a[1].isDouble() && std::isfinite(a[0].toDouble()) && std::isfinite(a[1].toDouble()), "Expected two finite coordinates");
  return {a[0].toDouble(),a[1].toDouble()};
}
Point normal(QJsonObject o, bool required) {
  check(!required || o.contains("inward_normal"), "Laser response requires inward_normal at every source and receiver");
  Point n=o.contains("inward_normal") ? point(o["inward_normal"]) : Point{0,1};
  check(std::abs(std::hypot(n.x,n.z)-1)<1e-8, "inward_normal must be a unit vector");
  return n;
}
QString localFile(const QString &root, QJsonValue value) {
  QString name=value.toString();
  check(!name.isEmpty() && !QDir::isAbsolutePath(name) && !name.contains('\\'), "Use a relative forward-slash input path");
  QString path=QFileInfo(root+"/"+name).canonicalFilePath();
  check(!path.isEmpty() && path.startsWith(QFileInfo(root).canonicalFilePath()+"/"), "Input path is missing or outside the dataset: " + name);
  return QDir::cleanPath(name);
}
}

void Engine::loadGeneral(const QJsonObject &d) {
  check(d["kind"]=="SAFT acquisition dataset", "Not a SAFT acquisition dataset");
  general=bundled=true; datasetKind="general_2d";
  datasetName=d["name"].toString();
  check(!datasetName.isEmpty(), "Dataset name is required");
  manifestHash=shaFile(data+"/dataset.json");
  QString geometryFile=localFile(data,d["geometry"]), acquisitionFile=localFile(data,d["acquisition"]);
  geometry=readJson(data+"/"+geometryFile);
  acquisition=readJson(data+"/"+acquisitionFile);
  check(geometry["schema"]==2 && geometry["coordinate_system"]=="x_z_m", "Geometry requires schema 2 and coordinate_system x_z_m");
  check(acquisition["schema"]==2 && acquisition["velocity_units"]=="m/s", "Acquisition requires schema 2 and velocity_units m/s");
  inputHashes[geometryFile]=shaFile(data+"/"+geometryFile);
  inputHashes[acquisitionFile]=shaFile(data+"/"+acquisitionFile);
  auto m=geometry["material"].toObject();
  double cl=number(m,"cl_m_s"),cs=number(m,"cs_m_s"),cr=number(m,"cr_m_s");
  check(cl>cs && cs>cr && cr>0 && number(m,"rho_kg_m3")>0, "Require cl > cs > cr > 0 and positive density");
  auto loops=geometry["loops"].toArray();
  check(!loops.empty(), "Geometry needs an outer polygon");
  std::set<QString> boundaryIds;
  std::vector<int> edgeLoops, edgeIndices;
  double xmin=INFINITY,xmax=-INFINITY,zmin=INFINITY,zmax=-INFINITY;
  for(int i=0;i<loops.size();++i) {
    auto loop=loops[i].toObject();
    check(loop["role"]==(i==0?"outer":"hole"), "First loop must be outer; remaining loops must be holes");
    auto vertices=loop["vertices_m"].toArray(), names=loop["boundary_ids"].toArray();
    check(vertices.size()>=3 && names.size()==vertices.size(), "Supply one boundary id per polygon edge; omit the repeated closing vertex");
    std::vector<Point> poly;
    for(auto v:vertices) {
      auto p=point(v); poly.push_back(p);
      if(i==0) { xmin=std::min(xmin,p.x);xmax=std::max(xmax,p.x);zmin=std::min(zmin,p.z);zmax=std::max(zmax,p.z); }
    }
    double area=0;
    for(size_t j=0;j<poly.size();++j) {
      auto a=poly[j],b=poly[(j+1)%poly.size()];
      double length=std::hypot(b.x-a.x,b.z-a.z);
      check(length>1e-10, "Zero-length or repeated polygon edge");
      area+=a.x*b.z-a.z*b.x;
      QString name=id(names[int(j)]);
      check(boundaryIds.insert(name).second, "Duplicate boundary id: " + name);
      boundaries.push_back({name,a,b,{-(b.z-a.z)/length,(b.x-a.x)/length}});
      edgeLoops.push_back(i); edgeIndices.push_back(int(j));
    }
    check(i==0?area>1e-20:area<-1e-20, "Outer loop must have positive signed area; holes negative signed area in x,z coordinates");
    polygons.push_back(std::move(poly));
  }
  // Reject crossings, overlapping edges, and touching loops. Adjacent edges may
  // meet only at their shared vertex. The Python validator performs the same checks.
  auto cross=[](Point a,Point b){return a.x*b.z-a.z*b.x;};
  auto sub=[](Point a,Point b){return Point{a.x-b.x,a.z-b.z};};
  for(size_t i=0;i<boundaries.size();++i) for(size_t j=i+1;j<boundaries.size();++j) {
    auto a=boundaries[i],b=boundaries[j];
    Point u=sub(a.b,a.a),v=sub(b.b,b.a),w=sub(b.a,a.a);
    double det=cross(u,v);
    bool shared=std::hypot(a.a.x-b.a.x,a.a.z-b.a.z)<1e-12 || std::hypot(a.a.x-b.b.x,a.a.z-b.b.z)<1e-12 ||
                std::hypot(a.b.x-b.a.x,a.b.z-b.a.z)<1e-12 || std::hypot(a.b.x-b.b.x,a.b.z-b.b.z)<1e-12;
    int apart=std::abs(edgeIndices[i]-edgeIndices[j]);
    bool adjacent=edgeLoops[i]==edgeLoops[j] && (apart==1 || apart==int(polygons[edgeLoops[i]].size())-1);
    check(!shared || adjacent, "Nonadjacent polygon boundaries touch");
    if(std::abs(det)>1e-20) {
      double t=cross(w,v)/det,s=cross(w,u)/det;
      check(!(t>=-1e-10 && t<=1+1e-10 && s>=-1e-10 && s<=1+1e-10 && !shared), "Polygon boundaries cross or touch");
    } else if(std::abs(cross(w,u))<1e-20) {
      double uu=u.x*u.x+u.z*u.z;
      double t0=(w.x*u.x+w.z*u.z)/uu,t1=t0+(v.x*u.x+v.z*u.z)/uu;
      check(std::min(1.,std::max(t0,t1))-std::max(0.,std::min(t0,t1))<1e-10, "Polygon edges overlap");
    }
  }
  for(size_t i=1;i<polygons.size();++i) for(auto p:polygons[i])
    check(containsPolygon(p), "Hole lies outside the specimen or inside another hole");
  auto sampling=acquisition["sampling"].toObject();
  samples=integer(sampling,"samples",3); dt=number(sampling,"dt_s"); startTime=number(sampling,"start_s");
  check(dt>0 && std::isfinite(startTime+dt*(samples-1)), "Invalid time grid");
  filterHalf=integer(acquisition,"filter_half_samples",1);
  check(2*filterHalf+2<=samples && filterHalf<=1000000, "Filter support must fit inside the record");
  auto response=acquisition["response"].toObject();
  QString responseModel=response["model"].toString();
  check(responseModel=="isotropic" || responseModel=="laser_free_surface", "Unknown response model");
  isotropic=responseModel=="isotropic";
  double width=isotropic ? 0 : number(response,"line_half_width_m");
  double pulse=isotropic ? 0 : number(response,"pulse_fwhm_s");
  check(isotropic || (width>0 && pulse>0), "Laser width and pulse duration must be positive");
  auto location=[&](QJsonObject o) {
    auto p=point(o["position_m"]); auto n=normal(o,!isotropic);
    check(containsPolygon(p), "Source or receiver is outside the specimen");
    if(!isotropic) {
      bool on=false;
      for(const auto &edge:boundaries) {
        double along=(p.x-edge.a.x)*(edge.b.x-edge.a.x)+(p.z-edge.a.z)*(edge.b.z-edge.a.z);
        double size=std::pow(edge.b.x-edge.a.x,2)+std::pow(edge.b.z-edge.a.z,2);
        if(std::abs((p.x-edge.a.x)*edge.normal.x+(p.z-edge.a.z)*edge.normal.z)<1e-10 && along>1e-16 && along<size-1e-16 &&
           std::hypot(n.x-edge.normal.x,n.z-edge.normal.z)<1e-8) on=true;
      }
      check(on, "Laser endpoints need a face-interior position and its inward normal");
    }
    return std::make_pair(p,n);
  };
  auto sourceRows=acquisition["sources"].toArray(),receiverRows=acquisition["receivers"].toArray();
  check(!sourceRows.empty() && !receiverRows.empty(), "Supply sources and receivers");
  std::set<QString> sourceIds,receiverLabels;
  for(auto v:sourceRows) {
    auto o=v.toObject(); QString name=id(o["id"]);
    check(sourceIds.insert(name).second,"Duplicate source id");
    auto [p,n]=location(o);
    sources.push_back({name,name,p.x,number(o,"trigger_s"),p.z,n});
  }
  groupNames.clear();
  for(auto v:receiverRows) {
    auto o=v.toObject(); QString name=id(o["id"]),group=id(o["group"]);
    check(receiverLabels.insert(name).second,"Duplicate receiver id");
    if(!groupNames.contains(group)) groupNames.append(group);
    auto [p,n]=location(o);
    receivers.push_back({int(receivers.size()),int(groupNames.indexOf(group)),0,p.x,p.z,n,name});
  }
  signalCases.clear();
  auto cases=acquisition["cases"].toArray(),components=acquisition["components"].toArray();
  QStringList componentNames;
  for(auto v:cases) { QString name=id(v); check(name!="difference"&&!signalCases.contains(name),"Duplicate or reserved case name"); signalCases.append(name); }
  for(auto v:components) { QString name=v.toString(); check((name=="vx"||name=="vz")&&!componentNames.contains(name),"Components must be unique vx and/or vz");componentNames.append(name); }
  check(!signalCases.empty()&&!componentNames.empty(),"Cases and components are required");
  primaryCase=signalCases.front(); referenceCase.clear();
  if(acquisition.contains("difference")) {
    auto diff=acquisition["difference"].toObject(); primaryCase=diff["minuend"].toString();referenceCase=diff["subtrahend"].toString();
    check(primaryCase!=referenceCase && signalCases.contains(primaryCase)&&signalCases.contains(referenceCase),"Difference must name two available cases");
  }
  std::set<QString> keys;
  for(auto v:acquisition["recordings"].toArray()) {
    auto o=v.toObject(); QString source=o["source"].toString(),signal=o["case"].toString(),component=o["component"].toString();
    check(sourceIds.count(source)&&signalCases.contains(signal)&&componentNames.contains(component),"Unknown recording source, case, or component");
    QString key=source+"/"+signal+"/"+component;
    check(keys.insert(key).second,"Duplicate recording entry");
    int si=int(std::find_if(sources.begin(),sources.end(),[&](const Source &s){return s.id==source;})-sources.begin());
    QString file=localFile(data,o["path"]);
    check(!inputHashes.count(file),"Recording paths must be unique");
    inputHashes[file]=shaFile(data+"/"+file);
    recordings.push_back({file,si,signal,component});
  }
  check(recordings.size()==sources.size()*size_t(signalCases.size())*size_t(componentNames.size()),"Missing recording: require every source x case x component");
  auto screens=geometry["screening_boundaries"].toObject();
  backwallBoundary=screens["backwall"].toString();surfaceBoundary=screens["surface"].toString();
  for(auto name:{backwallBoundary,surfaceBoundary}) check(name.isEmpty()||boundaryIds.count(name),"Unknown screening boundary");
  auto branch=[&](QJsonObject o) {
    QString modes=o["modes"].toString();QStringList faces;
    for(auto v:o["boundaries"].toArray()) {QString name=v.toString();check(boundaryIds.count(name),"Unknown path boundary: "+name);faces.append(name);}
    check(QRegularExpression("^[LS]{1,5}$").match(modes).hasMatch()&&modes.size()==faces.size()+1,"A branch needs one L/S mode per leg and one boundary per reflection (up to four)");
    for(int i=1;i<faces.size();++i) check(faces[i]!=faces[i-1],"Adjacent reflections cannot use the same boundary");
    return std::make_pair(modes,faces.join('|'));
  };
  for(auto before:{QString("L"),QString("S")}) for(auto after:{QString("L"),QString("S")})
    paths.push_back({"1-1_"+before+after,before,after,before+" | "+after,{},{},"direct"});
  for(auto v:geometry["paths"].toArray()) {
    auto o=v.toObject();QString name=id(o["id"]);
    check(std::none_of(paths.begin(),paths.end(),[&](const Path &p){return p.id==name;}),"Duplicate path id");
    auto before=branch(o["source"].toObject()),after=branch(o["receiver"].toObject());
    paths.push_back({name,before.first,after.first,name,before.second,after.second,"boundary"});
  }
  QByteArray identity=QJsonDocument(d).toJson(QJsonDocument::Compact);
  for(const auto &[p,h]:inputHashes) identity+=p.toUtf8()+":"+h.toUtf8()+"\n";
  fingerprint=sha(identity);
  audit={{"status","passed"},{"input_fingerprint",fingerprint},{"sources",sourceRows},{"receivers",receiverRows},
         {"source_parameters",QJsonObject{{"line_half_width_m",width},{"pulse_fwhm_s",pulse},{"absorber_width_m",0.}}}};
  m["thickness_m"]=zmax-zmin;
  arrivalAudit={{"baseline_input_fingerprint",fingerprint},{"material",m},{"model_inputs",QJsonArray{}},{"model_evidence",QJsonArray{}}};
  geometry["bounds_m"]=QJsonArray{xmin,xmax,zmin,zmax};
}

Config Engine::config(QJsonObject selection) const {
  if(!general) return Config::fromJson(selection,receiverCount(),sourceCount(),groupNames.size());
  QJsonObject defaults{{"component",acquisition["components"].toArray().first()},
    {"cases",QJsonArray{primaryCase}}, {"x_min_mm",material.xMin*1000},{"x_max_mm",material.xMax*1000},
    {"depth_min_mm",material.zMin*1000},{"depth_max_mm",material.zMax*1000},
    {"spacing_mm",std::max(material.xMax-material.xMin,material.zMax-material.zMin)*1000/100},
    {"low_MHz",std::min(.5,.05e-6/dt)},{"high_MHz",std::min(5.,.4e-6/dt)},
    {"minimum_sources",std::min(2,sourceCount())},{"minimum_receivers",std::min(8,receiverCount())},
    {"minimum_effective_receivers",std::min(8,receiverCount())}, {"min_offset_mm",0.}};
  for(auto it=selection.begin();it!=selection.end();++it) defaults[it.key()]=it.value();
  auto result=Config::fromJson(defaults,receiverCount(),sourceCount(),groupNames.size());
  check(acquisition["components"].toArray().contains(result.component),"Requested velocity component is absent");
  bool back=result.mutes[3]||result.mutes[4]||result.mutes[5]||std::any_of(result.reflectedScreen.begin(),result.reflectedScreen.end(),[](bool b){return b;});
  check(!back||!backwallBoundary.isEmpty(),"Reflected nuisance suppression requires geometry.screening_boundaries.backwall");
  check(!result.reflectedScreen[3]||!surfaceBoundary.isEmpty(),"Repeated-reflection suppression requires geometry.screening_boundaries.surface");
  return result;
}

std::shared_ptr<Vec> Engine::rawGeneral(int si,const QString &signal,const QString &component) {
  check(si>=0&&si<sourceCount(),"Invalid source index");
  QString key=QString::number(si)+"/"+signal+"/"+component;
  if(rawCache.count(key)) return rawCache[key];
  auto result=std::make_shared<Vec>(size_t(receiverCount())*samples);
  if(signal=="difference") {
    check(!referenceCase.isEmpty(),"No paired difference defined in acquisition.json");
    auto a=rawGeneral(si,primaryCase,component),b=rawGeneral(si,referenceCase,component);
    for(size_t i=0;i<result->size();++i) (*result)[i]=(*a)[i]-(*b)[i];
  } else {
    auto entry=std::find_if(recordings.begin(),recordings.end(),[&](const Recording &r){return r.source==si&&r.signal==signal&&r.component==component;});
    check(entry!=recordings.end(),"Requested recording is absent: "+key);
    QFile file(data+"/"+entry->file);check(file.open(QIODevice::ReadOnly),"Cannot read recording: "+entry->file);
    QByteArray bytes=file.readAll();check(sha(bytes)==inputHashes.at(entry->file),"Recording changed since dataset was opened: "+entry->file);
    const char *p=bytes.constData(),*end=p+bytes.size(),*line=std::find(p,end,'\n');
    check(line!=end,"Missing CSV header");QByteArray header(p,line-p);header=header.trimmed();
    if(header.startsWith("\xef\xbb\xbf")) header=header.mid(3);
    QByteArray expected="time_s";for(const auto &r:receivers) expected+=","+r.label.toUtf8();
    check(header==expected,"CSV header must be time_s followed by receiver ids in acquisition order: "+entry->file);
    p=line+1;
    double previousTime=-INFINITY;
    for(int t=0;t<samples;++t) {
      for(int r=-1;r<receiverCount();++r) {
        while(p<end&&(*p==' '||*p=='\t')) ++p;
        double value=0;auto parsed=std::from_chars(p,end,value);
        check(parsed.ec==std::errc()&&std::isfinite(value),"Missing or nonfinite CSV sample: "+entry->file);p=parsed.ptr;
        while(p<end&&(*p==' '||*p=='\t')) ++p;
        if(r<0) {
          check(value>previousTime && std::abs(value-(startTime+t*dt))<=std::max(1e-15,dt*1e-7),"CSV time grid disagrees with acquisition.sampling: "+entry->file);
          previousTime=value;
        }
        else (*result)[size_t(r)*samples+t]=value;
        if(r+1<receiverCount()) {check(p<end&&*p==',',"Incorrect CSV column count");++p;}
      }
      if(p<end&&*p=='\r') ++p;
      if(t+1<samples) {check(p<end&&*p=='\n',"Incorrect CSV row length");++p;}
    }
    while(p<end&&std::isspace(static_cast<unsigned char>(*p))) ++p;
    check(p==end,"Extra CSV samples or columns");++csvParses;++loads[entry->file];
  }
  rawCache[key]=result;return result;
}

double Engine::sourceAngle(const Ray &r,int s) const {
  if(!general) return r.first;
  auto n=sources[s].normal;return std::atan2(r.departure.x*n.z-r.departure.z*n.x,r.departure.x*n.x+r.departure.z*n.z);
}
double Engine::receiverAngle(const Ray &r,int receiver) const {
  if(!general) return r.last;
  auto n=receivers[receiver].normal;return std::atan2(r.arrival.x*n.z-r.arrival.z*n.x,-r.arrival.x*n.x-r.arrival.z*n.z);
}
Complex Engine::receiverResponse(const Ray &r,int receiver,const Config &c,const Curves &curves,char mode) const {
  if(isotropic) return 1.;
  int column=c.component=="vx"?0:1,mi=mode=='S';
  auto surface=reflect(receiverAngle(r,receiver),mode,material).surface;
  if(!general) return (column?-1.:1.)*surface[column]/curves.norm[mi][column];
  auto n=receivers[receiver].normal;
  // Rotate local tangential/normal surface motion into the CSV's global x,z axes.
  double tangent=column?-n.x:n.z,normalComponent=column?n.z:n.x;
  double scale=std::abs(tangent)*curves.norm[mi][0]+std::abs(normalComponent)*curves.norm[mi][1];
  return scale>0?(tangent*surface[0]+normalComponent*surface[1])/scale:Complex(0.);
}
} // namespace saft
