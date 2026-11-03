// Headless version of the original Workbench --compute exporter.
#include "saft.hpp"
#include <iostream>
int main(int argc, char **argv) {
  QCoreApplication app(argc, argv);
  QCommandLineParser parser;
  parser.addHelpOption();
  parser.addOptions({{"repository", "Repository/cache root", "path"},
                     {"dataset", "Acquisition folder", "path"},
                     {"compute", "Request JSON", "json"},
                     {"output", "Output directory", "path"}});
  parser.process(app);
  try {
    if (!parser.isSet("compute") || !parser.isSet("dataset") || !parser.isSet("output"))
      throw std::runtime_error("Use --dataset DATASET --compute REQUEST --output OUTPUT");
    QString root = parser.isSet("repository") ? QDir(parser.value("repository")).absolutePath() : QDir::currentPath();
      QElapsedTimer init;
      init.start();
      saft::Engine engine(root, parser.value("dataset"));
      double startup = init.nsecsElapsed() / 1e9;
      auto request = saft::readJson(parser.value("compute"));
      QJsonArray geometry;
      for (auto v : request["ray_checks"].toArray()) {
        auto q = v.toObject();
        auto a = q["a_m"].toArray(), b = q["b_m"].toArray();
        const auto r =
            engine.trace({a[0].toDouble(), a[1].toDouble()},
                         {b[0].toDouble(), b[1].toDouble()},
                         q["modes"].toString(), q["outgoing"].toBool(), true, q["via"].toString());
        QJsonArray vertices;
        for (auto p : r.points)
          vertices.append(QJsonArray{p.x, p.z});
        geometry.append(QJsonObject{{"valid", r.valid},
                                    {"time_s", r.time},
                                    {"points_m", vertices},
                                    {"first", r.first},
                                    {"last", r.last},
                                    {"wall_real", r.wall.real()}, {"wall_imag", r.wall.imag()}});
      }
      QJsonArray batch = request["requests"].toArray();
      if (batch.empty())
        batch.append(request);
      QJsonArray reports;
      QString output = parser.value("output");
      if (output.isEmpty())
        output = root + "/studies/laser/build/saft_cpp_computation";
      QDir().mkpath(output);
      saft::writeJson(output + "/path_catalogue.json", engine.grooveCatalogue());
      if (!geometry.empty())
        saft::writeJson(output + "/geometry.json", {{"rays", geometry}});
      if (request["verify_inputs"].toBool())
        saft::writeJson(output + "/input_verification.json",
                        engine.verifyInputs({}));
      for (int i = 0; i < batch.size(); ++i) {
        auto job = batch[i].toObject();
        auto c = saft::Config::fromJson(job["selection"].toObject(),
                                        engine.receiverCount());
        std::vector<saft::Point> points;
        for (auto v : job["points_m"].toArray()) {
          auto a = v.toArray();
          points.push_back({a[0].toDouble(), a[1].toDouble()});
        }
        if (job["clear_memory"].toBool())
          engine.clearMemory();
        if (job["forbid_reference"].toBool()) {
          for (auto it = engine.inputHashes.begin();
               it != engine.inputHashes.end();)
            if (it->first.contains("/reference/"))
              it = engine.inputHashes.erase(it);
            else
              ++it;
        }
        auto result =
            engine.reconstruct(c, {}, points.empty() ? nullptr : &points);
        QString folder = output + "/" + QString::number(i);
        engine.exportResult(folder, *result);
        QJsonObject r = result->timing();
        r["path"] = c.path;
        r["index"] = i;
        r["pixels"] = qint64(result->points.size());
        r["csv_parses"] = engine.csvParses;
        r["cache_reads"] = engine.cacheReads;
        r["filter_builds"] = engine.filterBuilds;
        r["geometry_builds"] = engine.geometryBuilds;
        if (job["inspect"].isObject()) {
          auto q = job["inspect"].toObject();
          auto inspection = engine.inspect(
              c, q["source"].toInt(), q["receiver"].toInt(94),
              q["signal"].toString("lof"), q["competitor"].toString("1-2_SSS"),
              q["has_point"].toBool(),
              {q["x_mm"].toDouble(2.5), q["depth_mm"].toDouble(20)});
          engine.exportInspection(folder, *inspection);
          QJsonArray times;
          for (double t : inspection->times)
            times.append(t);
          r["arrival_us"] = times;
          QJsonArray background;
          for (const auto &family : inspection->arrivals) {
            QJsonArray values;
            for (double t : family)
              values.append(t);
            background.append(values);
          }
          r["background_arrival_us"] = background;
        }
        reports.append(r);
      }
      saft::writeJson(output + "/benchmark.json",
                      {{"status", "passed"},
                       {"engine_initialization_s", startup},
                       {"runs", reports},
                       {"provenance", engine.provenance()}});
      return 0;

  } catch (const std::exception &e) { std::cerr << e.what() << std::endl; return 1; }
}
