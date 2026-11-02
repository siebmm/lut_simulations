#include "saft.hpp"
#include <iostream>
#include <cstring>

using namespace saft;
template<class T> void binary(const QString &path, const std::vector<T> &v) {
  QSaveFile file(path);
  if (!file.open(QIODevice::WriteOnly)) throw std::runtime_error("Cannot write array");
  file.write(reinterpret_cast<const char *>(v.data()), qint64(v.size()*sizeof(T)));
  if (!file.commit()) throw std::runtime_error("Cannot commit array");
}
std::shared_ptr<CVec> bank(const QString &folder, int source, QString signal, int count) {
  auto read = [&](QString name) {
    QFile f(folder+"/"+name+"_"+QString::number(source)+".c128");
    if (!f.open(QIODevice::ReadOnly) || f.size()!=qint64(count)*Samples*sizeof(Complex))
      throw std::runtime_error("Missing/wrong-sized processed signal bank");
    auto result=std::make_shared<CVec>(size_t(count)*Samples);
    if (f.read(reinterpret_cast<char *>(result->data()),f.size())!=f.size())
      throw std::runtime_error("Short signal bank read");
    return result;
  };
  if (signal!="difference") return read(signal);
  auto a=read("lof"),b=read("reference");
  for (size_t i=0;i<a->size();++i) (*a)[i]-=(*b)[i];
  return a;
}
int main(int argc,char **argv) {
  QCoreApplication app(argc,argv);
  try {
    if (argc!=5) throw std::runtime_error("Usage: SAFTFilterStudy ROOT DATASET REQUEST OUTPUT");
    Engine engine(argv[1],argv[2]);
    const auto request=readJson(argv[3]); const QString out=argv[4];
    QDir().mkpath(out);
    if (request["dump_raw"].toBool()) {
      writeJson(out+"/input_verification.json",engine.verifyInputs({}));
      for (auto signal:{QString("lof"),QString("reference")})
        for(int s=0;s<5;++s) binary(out+"/"+signal+"_"+QString::number(s)+".f64",*engine.raw(s,signal,"vy"));
      std::cout<<"raw export verified"<<std::endl; return 0;
    }
    QJsonArray timing;
    for (auto value:request["requests"].toArray()) {
      const auto job=value.toObject(); QString id=job["id"].toString();
      auto config=Config::fromJson(job["selection"].toObject(),engine.receiverCount());
      const auto processing=job["signal_processing"].toObject();
      std::map<QString,std::shared_ptr<CVec>> cache;
      if (processing.isEmpty()) engine.setSignalProvider({});
      else {
        auto folder=processing["bank"].toString();
        const auto meta=readJson(folder+"/bank.json");
        if(meta["dataset_fingerprint"].toString()!=engine.fingerprint)
          throw std::runtime_error("Signal bank belongs to a different acquisition");
        engine.setSignalProvider([&,folder](int s,QString signal,const Config &){
          QString key=QString::number(s)+"/"+signal;
          if(!cache.count(key)) cache[key]=bank(folder,s,signal,engine.receiverCount());
          return cache[key];
        },processing["support_half_samples"].toInt(Half),processing["known_initial_rest"].toBool());
      }
      std::vector<Point> points;
      for(auto item:job["points_m"].toArray()) {
        auto p=item.toArray();points.push_back({p[0].toDouble(),p[1].toDouble()});
      }
      auto result=engine.reconstruct(config,{},points.empty()?nullptr:&points);
      QString directory=out+"/"+id;QDir().mkpath(directory);
      binary(directory+"/points.f64",result->points);
      binary(directory+"/x.f64",result->x);binary(directory+"/z.f64",result->z);
      for(const auto &[signal,planes]:result->planes) {
        for(int variant=1;variant<=2;++variant) {
          QString base=directory+"/"+signal+"_"+QString::number(variant);
          const auto &p=planes[variant];
          binary(base+"_amplitude.f64",p.amplitude);
          binary(base+"_count.f64",p.count);binary(base+"_cf.f64",p.cf);
          binary(base+"_retained.f64",p.retained);binary(base+"_screened.f64",p.screened);
          binary(base+"_means.c128",p.means);binary(base+"_valid.u8",p.valid);
          binary(base+"_receiver_usage.u64",p.receiverUsage);
        }
      }
      writeJson(directory+"/request.json",job);
      writeJson(directory+"/timing.json",result->timing());
      if(job["export_windows"].toBool())
        for(const auto &[signal,screen]:result->screening)
          if(screen) writeJson(directory+"/"+signal+"_windows.json",screen->json());
      timing.append(QJsonObject{{"id",id},{"timing",result->timing()}});
      std::cout<<id.toStdString()<<" "<<result->seconds<<" s"<<std::endl;
    }
    writeJson(out+"/benchmark.json",{{"jobs",timing},{"provenance",engine.provenance()}});
  } catch(const std::exception &e) { std::cerr<<e.what()<<std::endl; return 1; }
}
