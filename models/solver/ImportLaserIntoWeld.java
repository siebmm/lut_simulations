import com.comsol.model.*;
import com.comsol.model.util.*;
import java.util.*;

/** Import the calibrated thermoelastic equivalent traction into copies of the weld pair. */
public class ImportLaserIntoWeld {
 static final String ROOT="C:\\Users\\siebm\\Documents\\COMSOL\\NDT\\studies\\laser";
 static void par(Model m,String name,String value,String description){m.param().set(name,value,description);}
 public static void main(String[] args) throws Exception {run();}
 public static Model run() throws Exception {
  Model last=null;
  for(String kind:new String[]{"reference","lof"}){
   Model m=ModelUtil.load("WeldLaser",ROOT+"/models/original_weld_geometry/baseline/weld_"+kind+"_v3.mph");
   m.modelPath(ROOT);
   m.label("weld_"+kind+"_laser.mph");
   m.comments("Equivalent thermoelastic line-source traction imported from the numerical calibration. Absorbed fluence is a normalization, not measured laser energy. Elastic-only, 2D plane strain; quantitative interpretation limited to the validated 0-5 MHz band. See docs/WELD_LASER_IMPORT.md. No defect-detection or experimental SNR claim.");
   par(m,"laser_F_abs","10[J/m^2]","Peak absorbed fluence; normalized small-signal benchmark, NOT incident energy");
   par(m,"laser_Cp","480[J/(kg*K)]","X65 benchmark heat capacity; calibration documentation gives source");
   par(m,"laser_alpha","11e-6[1/K]","X65 benchmark linear expansion coefficient");
   par(m,"laser_wavelength","1064[nm]","Literature metadata; optics are not solved");
   par(m,"laser_fwhm","6.5[ns]","Gaussian power FWHM; literature-informed assumption");
   par(m,"laser_sigma","laser_fwhm/(2*sqrt(2*log(2)))","Gaussian power standard deviation");
   par(m,"laser_t0","30[ns]","Pulse centre / trigger offset; use this timing in reconstruction");
   par(m,"laser_fit","1.0844342935972646","Fitted scalar from results/thermoelastic_source_calibration/calibration_comparison.json, 0-5 MHz");
   par(m,"laser_A","laser_fit*E0*laser_alpha*laser_F_abs/((1-nu0)*rho0*laser_Cp*w0)","Calibrated dipole coefficient, Pa; not the maximum traction");
   par(m,"laser_Eline","laser_F_abs*w0*sqrt(pi/2)","Full line-source absorbed energy per out-of-plane length, J/m");
   par(m,"laser_dt_cap","0.5[ns]","Upper cell-time-scale cap while resolving source rise");
   // Preserve obsolete benchmark parameters for provenance, but label them clearly.
   par(m,"p0",m.param().get("p0"),"UNUSED: original normal-pressure benchmark amplitude");
   for(String p:new String[]{"pulse_t0","pulse_sigma","pulse_duration"})
    par(m,p,m.param().get(p),"UNUSED: original pressure pulse; active source uses laser_* parameters");
   m.component("comp1").variable().create("laserSource");
   m.component("comp1").variable("laserSource").label("Calibrated thermoelastic equivalent line source");
   m.component("comp1").variable("laserSource").set("laser_xi","(x-x_source)/w0","Position relative to laser centre");
   m.component("comp1").variable("laserSource").set("laser_H","0.5*(1+erf((t-laser_t0)/(sqrt(2)*laser_sigma)))","Integrated laser power; retained step over the microsecond record");
   m.component("comp1").variable("laserSource").set("laser_tx","laser_A*4*laser_xi*exp(-2*laser_xi^2)*laser_H","Signed outward tangential dipole traction, Pa");
   m.component("comp1").physics("elte").feature("load1").set("forceType","ForceArea");
   // ForceArea reads forceReferenceArea in COMSOL 6.4; legacy force is inactive.
   m.component("comp1").physics("elte").feature("load1").set("force",new String[]{"0","0","0"});
   m.component("comp1").physics("elte").feature("load1").set("forceReferenceArea_src","userdef");
   m.component("comp1").physics("elte").feature("load1").set("forceReferenceArea",new String[]{"laser_tx","0","0"});
   m.component("comp1").physics("elte").feature("load1").label("Laser: calibrated thermoelastic tangential dipole");
   m.component("comp1").selection("source").label("Flat exposed weld surface: thermoelastic line source");
   // Independently identify the intended flat boundary before accepting the existing selection.
   m.component("comp1").selection().create("laserSurfaceCheck","Box");
   m.component("comp1").selection("laserSurfaceCheck").set("entitydim",1);
   m.component("comp1").selection("laserSurfaceCheck").set("xmin","-groove_floor_half-1e-8[m]");
   m.component("comp1").selection("laserSurfaceCheck").set("xmax","groove_floor_half+1e-8[m]");
   m.component("comp1").selection("laserSurfaceCheck").set("ymin","weld_y-1e-8[m]");
   m.component("comp1").selection("laserSurfaceCheck").set("ymax","weld_y+1e-8[m]");
   m.component("comp1").selection("laserSurfaceCheck").set("condition","inside");
   int[] source=m.component("comp1").selection("source").entities(1);
   int[] check=m.component("comp1").selection("laserSurfaceCheck").entities(1);
   if(source.length!=1||!Arrays.equals(source,check))throw new IllegalStateException("Source boundary mismatch: "+Arrays.toString(source)+" vs "+Arrays.toString(check));
   if(m.component("comp1").physics("elte").feature("frac1").isActive()!=kind.equals("lof"))throw new IllegalStateException("Incorrect fracture state");
   if(Math.abs(m.param().evaluate("x_source"))+4*m.param().evaluate("w0")>=m.param().evaluate("groove_floor_half"))throw new IllegalStateException("Source too close to groove corner for calibrated flat source");
   // Keep the original stability bound, with a tighter bound through the source rise.
   m.sol("sol1").feature("t1").set("exprs",new String[]{"min(root.comp1.elte.wtc,if(t<laser_t0+8*laser_sigma,laser_dt_cap,1[s]))"});
   m.sol("sol1").feature("t1").set("updtstep","manual");
   m.sol("sol1").feature("t1").set("ntimestepsupdate",10);
   m.sol("sol1").feature("t1").set("exprsgporderupdt",new int[]{2});
   // Retain all 113 observation points and the original displacement/velocity channels.
   // Compute Displacement lives on boundaries: domain interpolation cannot evaluate it.
   for(String tag:new String[]{"left","right","weld"}){
    m.result().dataset(tag).set("bndsnap",true);
   }
   for(String tag:m.result().export().tags()){
    if(!m.result().export(tag).getType().equals("Data"))continue;
    m.result().export(tag).set("filename",ROOT+"/results/full_record_receiver_exports/"+kind+"_"+tag+".txt");
    m.result().export(tag).set("header",true);
    m.result().export(tag).set("fullprec",true);
   }
   m.study("std1").label("Elastic propagation: calibrated thermoelastic equivalent source");
   m.save(ROOT+"/models/weld_"+kind+"_laser.mph");
   System.out.println("LASER_IMPORT_OK "+kind+" boundary="+Arrays.toString(source)+" A_Pa="+m.param().evaluate("laser_A")+" fracture="+m.component("comp1").physics("elte").feature("frac1").isActive());
   last=m;
   if(kind.equals("reference"))ModelUtil.remove("WeldLaser");
  }
  return last;
 }
}
