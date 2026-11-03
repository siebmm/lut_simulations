#include "saft.hpp"
#include <iostream>

using namespace saft;
template<class T> void binary(const QString &name,const std::vector<T> &values) {
  QSaveFile file(name);
  if (!file.open(QIODevice::WriteOnly)) throw std::runtime_error("Cannot write admission evidence");
  if (file.write(reinterpret_cast<const char *>(values.data()),qint64(values.size()*sizeof(T)))!=qint64(values.size()*sizeof(T)))
    throw std::runtime_error("Short admission evidence write");
  if (!file.commit()) throw std::runtime_error("Cannot commit admission evidence");
}
int main(int argc,char **argv) {
  QCoreApplication app(argc,argv);
  try {
    if(argc!=5) throw std::runtime_error("Usage: SAFTAdmissionAudit ROOT DATASET REQUEST OUTPUT");
    Engine engine(argv[1],argv[2]);
    const auto request=readJson(argv[3]); const QString out=argv[4];
    QDir().mkpath(out);
    for(auto value:request["requests"].toArray()) {
      const auto job=value.toObject();
      auto config=engine.config(job["selection"].toObject());
      auto result=engine.reconstruct(config);
      const QString folder=out+"/"+job["id"].toString();QDir().mkpath(folder);
      binary(folder+"/points.f64",result->points);
      binary(folder+"/x.f64",result->x);binary(folder+"/z.f64",result->z);
      for(const auto &[signal,planes]:result->planes) {
        const auto &p=planes[2];const QString prefix=folder+"/"+signal;
        binary(prefix+"_amplitude.f64",p.amplitude);
        binary(prefix+"_source_count.f64",p.count);
        binary(prefix+"_means.c128",p.means);
        binary(prefix+"_valid.u8",p.valid);
        binary(prefix+"_effective.f64",p.effective);
        Vec accepted,eligible,candidates,screened;
        for(const auto &stats:p.stats) {
          accepted.push_back(stats.count);eligible.push_back(stats.eligible);
          candidates.push_back(stats.candidates);screened.push_back(stats.screened);
        }
        binary(prefix+"_accepted.f64",accepted);binary(prefix+"_eligible.f64",eligible);
        binary(prefix+"_complete.f64",candidates);binary(prefix+"_screened.f64",screened);
      }
      writeJson(folder+"/request.json",job);writeJson(folder+"/timing.json",result->timing());
      std::cout<<job["id"].toString().toStdString()<<" "<<result->seconds<<" s"<<std::endl;
    }
    writeJson(out+"/provenance.json",engine.provenance());
  } catch(const std::exception &e) {std::cerr<<e.what()<<std::endl;return 1;}
}
