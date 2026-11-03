import com.comsol.model.*;
import com.comsol.model.util.*;
import java.util.*;

/** Solve each case and export traces, preserving an intermediate thermal solution. */
public class RunLaserCalibration {
 static void writeRunParameters(Model m, boolean complete) throws Exception {
  java.nio.file.Path path=java.nio.file.Paths.get("results/thermoelastic_source_calibration/exports/calibration_parameters.json");
  java.nio.file.Files.createDirectories(path.getParent());
  String json="{\n  \"complete\": "+complete+",\n  \"absorbed_fluence_J_m2\": "+m.param().evaluate("F_abs")+"\n}\n";
  java.nio.file.Files.write(path,json.getBytes(java.nio.charset.StandardCharsets.UTF_8));
 }
 static void configureExports(Model m) {
  String[] tags={"arcT","arcD","surfT","surfD"};
  String[] names={"calibration_thermal_arc.txt","calibration_dipole_arc.txt","calibration_thermal_surface.txt","calibration_dipole_surface.txt"};
  for(int i=0;i<tags.length;i++)m.result().export(tags[i]).set("filename",ImportLaserIntoWeld.ROOT+"/results/thermoelastic_source_calibration/exports/"+names[i]);
 }
 public static Model run() throws Exception { return run(false); }
 public static Model run(boolean useCuDss) throws Exception {
  Model m=ModelUtil.load("Calibration","models/laser_source_calibration.mph");
  configureExports(m);
  writeRunParameters(m,false);
  for(String sol:new String[]{"sol1","sol2"}){
   m.sol(sol).feature("t1").set("maxstepconstraintgenalpha","expr");
   m.sol(sol).feature("t1").set("maxstepexpressiongenalpha","if(t<60[ns],0.5[ns],max_dt)");
   for(String field:m.sol(sol).feature("v1").feature().tags()){
    System.out.println("DEPENDENT "+sol+" "+field);
    if(field.equals("comp1_u")||field.equals("comp1_udisp")){
     m.sol(sol).feature("v1").feature(field).set("scalemethod","manual");
     m.sol(sol).feature("v1").feature(field).set("scaleval","1e-10");
    }
   }
  }
  if(useCuDss)LaserLinearSolver.configureCuDss(m,"sol1","sol2");
  m.save("models/laser_source_calibration.mph");
  System.out.println("SOLVING_THERMOELASTIC");
  m.sol("sol1").runAll();
  try {m.save("models/laser_source_calibration_solved.mph");}catch(Exception e){throw new RuntimeException(e);}
  m.result().export("arcT").run();
  m.result().export("surfT").run();
  System.out.println("THERMOELASTIC_EXPORTED");
  System.out.println("SOLVING_DIPOLE");
  m.sol("sol2").runAll();
  m.save("models/laser_source_calibration_solved.mph");
  m.result().export("arcD").run();
  m.result().export("surfD").run();
  try {m.save("models/laser_source_calibration_solved.mph");}catch(Exception e){throw new RuntimeException(e);}
  writeRunParameters(m,true);
  System.out.println("CALIBRATION_SOLVE_OK");
  return m;
 }
 public static void main(String[] args) throws Exception {run();}
}
