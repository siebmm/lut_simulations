#include "saft.hpp"
#include <iostream>
using namespace saft;
template<class T> void binary(const QString &path,const std::vector<T> &v) {
  QSaveFile f(path);
  if(!f.open(QIODevice::WriteOnly) || f.write(reinterpret_cast<const char*>(v.data()),qint64(v.size()*sizeof(T)))!=qint64(v.size()*sizeof(T)) || !f.commit())
    throw std::runtime_error("Cannot write review statistics");
}
int main(int argc,char **argv) {
 QCoreApplication app(argc,argv);
 try {
  if(argc!=5) throw std::runtime_error("Usage: SAFTReview ROOT DATASET REQUEST OUTPUT");
  Engine engine(argv[1],argv[2]);const auto request=readJson(argv[3]);const QString out=argv[4];QDir().mkpath(out);
  writeJson(out+"/catalogue.json",engine.grooveCatalogue());
  for(auto item:request["requests"].toArray()) {
   const auto job=item.toObject();auto c=Config::fromJson(job["selection"].toObject(),engine.receiverCount());
   auto r=engine.reconstruct(c);const QString folder=out+"/"+job["id"].toString();QDir().mkpath(folder);
   binary(folder+"/points.f64",r->points);binary(folder+"/x.f64",r->x);binary(folder+"/z.f64",r->z);
   for(const auto &[signal,planes]:r->planes) for(int variant=2;variant<=2;++variant) {
    const auto &p=planes[variant];const QString pre=folder+"/"+signal+"_"+QString::number(variant);
    binary(pre+"_amplitude.f64",p.amplitude);binary(pre+"_count.f64",p.count);
    binary(pre+"_means.c128",p.means);binary(pre+"_valid.u8",p.valid);
    Vec stats;CVec sums;
    for(const auto &v:p.stats) {
     for(double x:{v.count,v.eligible,v.weight,v.weight2,v.power,v.candidates,v.screened})stats.push_back(x);
     sums.push_back(v.sum);
    }
    binary(pre+"_stats.f64",stats);binary(pre+"_sums.c128",sums);
   }
   writeJson(folder+"/request.json",job);writeJson(folder+"/timing.json",r->timing());
   std::cout<<job["id"].toString().toStdString()<<" "<<r->seconds<<" s"<<std::endl;
  }
  writeJson(out+"/provenance.json",engine.provenance());
 }catch(const std::exception &e){std::cerr<<e.what()<<std::endl;return 1;}
}
