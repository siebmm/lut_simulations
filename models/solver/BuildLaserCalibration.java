import com.comsol.model.*;
import com.comsol.model.util.*;
import java.util.*;

/** Local, symmetric plane-strain source benchmark. All lengths in SI units. */
public class BuildLaserCalibration {
 static Model m;
 static void par(String n,String v,String d){m.param().set(n,v,d);}
 static void box(String tag,int dim,String xmin,String xmax,String ymin,String ymax){
  m.component("comp1").selection().create(tag,"Box");
  m.component("comp1").selection(tag).set("entitydim",dim);
  m.component("comp1").selection(tag).set("xmin",xmin);
  m.component("comp1").selection(tag).set("xmax",xmax);
  m.component("comp1").selection(tag).set("ymin",ymin);
  m.component("comp1").selection(tag).set("ymax",ymax);
  m.component("comp1").selection(tag).set("condition","inside");
 }
 static void rect(String tag,String bottom,String height){
  m.component("comp1").geom("geom1").create(tag,"Rectangle");
  m.component("comp1").geom("geom1").feature(tag).set("size",new String[]{"half_width",height});
  m.component("comp1").geom("geom1").feature(tag).set("pos",new String[]{"0",bottom});
  m.component("comp1").geom("geom1").feature(tag).set("selresult",true);
 }
 static void dataset(String tag,String solution,String[] x,String[] y){
  m.result().dataset().create(tag,"CutPoint2D");
  m.result().dataset(tag).set("data",solution);
  m.result().dataset(tag).set("pointx",x);
  m.result().dataset(tag).set("pointy",y);
 }
 static void export(String tag,String data,String[] expr,String[] units,String file){
  m.result().export().create(tag,"Data");
  m.result().export(tag).set("data",data);
  m.result().export(tag).set("expr",expr);
  m.result().export(tag).set("unit",units);
  m.result().export(tag).set("filename",file);
  m.result().export(tag).set("header",true);
  m.result().export(tag).set("fullprec",true);
 }
 static void study(String tag,boolean thermal){
  m.study().create(tag);
  m.study(tag).label(thermal?"1 - Heat input and thermoelastic wave generation":"2 - Equivalent dipole, elastic only");
  m.study(tag).create("time","Transient");
  m.study(tag).feature("time").set("tlist","range(0,dt_out,t_end)");
  m.study(tag).feature("time").activate("ht",thermal);
  m.study(tag).feature("time").activate("solid",thermal);
  m.study(tag).feature("time").activate("dip",!thermal);
  m.study(tag).createAutoSequences("all");
 }
 public static Model run() {
  m=ModelUtil.create("LaserCalibration");
  m.modelPath("C:\\Users\\siebm\\Documents\\COMSOL\\NDT\\studies\\laser");
  m.label("laser_source_calibration.mph");
  m.comments("Literature-informed, normalized thermoelastic line-source calibration. 6.5 ns duration from Sensors 2021,21,7105 treated as Gaussian FWHM; X65 thermal data from doi:10.1002/srin.202300309 Table 4. Absorbed fluence and line geometry are benchmark assumptions. Symmetric local flat patch; not a reproduction of a pipe experiment. See docs/LASER_CALIBRATION.md.");
  par("E0","212442565550.3707[Pa]","Inherited weld-model elastic modulus");
  par("nu0","0.2889775930054178","Inherited weld-model Poisson ratio");
  par("rho0","7859.97504274123[kg/m^3]","Inherited weld-model density");
  par("k0","52[W/(m*K)]","X65 thermal conductivity; Dolz et al. Table 4");
  par("Cp0","480[J/(kg*K)]","X65 heat capacity; Table 4, converted from mJ/(tonne K)");
  par("alpha0","11e-6[1/K]","X65 linear expansion; Dolz et al. Table 4");
  par("T0","293.15[K]","Ambient and strain reference temperature");
  par("lambda_laser","1064[nm]","Literature laser wavelength; optical field not solved");
  par("pulse_fwhm","6.5[ns]","Reported pulse duration interpreted as Gaussian FWHM");
  par("sigma_t","pulse_fwhm/(2*sqrt(2*log(2)))","Gaussian temporal standard deviation");
  par("t0","30[ns]","Pulse centre; computational trigger offset");
  par("w0","0.5[mm]","1/e^2 line half-width inherited from weld benchmark; not pipe-paper beam diameter");
  par("fluence_scale","10","Absorbed fluence multiplier relative to the solved 10 J/m^2 baseline");
  par("F_abs","fluence_scale*10[J/m^2]","Peak absorbed fluence; default predicts about 78 K rise with constant properties");
  par("Eline_abs","F_abs*w0*sqrt(pi/2)","Full symmetric line: absorbed energy per out-of-plane length");
  par("A_guess","E0*alpha0*F_abs/((1-nu0)*rho0*Cp0*w0)","Thin heated-layer tangential dipole scale; validate against thermal model");
  par("dipole_scale","1","Multiplier to fit to thermoelastic receiver waveforms");
  par("half_width","4[mm]","Half-model x extent; mirror across x=0");
  par("depth","4[mm]","Local calibration depth; not the pipe wall thickness");
  par("h_bulk","0.08[mm]","Quadratic bulk mesh target for analysis up to 5 MHz");
  par("n_strip","80","Mapped surface strip elements across the half-width");
  par("dt_out","2[ns]","Uniform receiver output spacing");
  par("t_end","0.9[us]","Early comparison window, before outer-boundary returns");
  par("max_dt","2[ns]","Internal step cap after pulse generation");
  par("r_obs","2[mm]","Bulk observation arc radius");
  par("cP","sqrt(E0*(1-nu0)/(rho0*(1+nu0)*(1-2*nu0)))","Longitudinal speed");
  par("cS","sqrt(E0/(2*rho0*(1+nu0)))","Shear speed");
  m.component().create("comp1",true);
  m.component("comp1").label("Local laser source; symmetric plane strain");
  m.component("comp1").geom().create("geom1",2);
  m.component("comp1").geom("geom1").lengthUnit("m");
  m.component("comp1").geom("geom1").repairTolType("absolute");
  m.component("comp1").geom("geom1").absRepairTol(1e-12);
  double d=0;
  for(int i=0;i<9;i++){
   double h=5e-8*Math.pow(2,i);
   rect("layer"+i,Double.toString(-d-h),Double.toString(h));
   d+=h;
  }
  par("thermal_depth",Double.toString(d)+"[m]","Graded thermal skin: 50 nm top layer to 12.8 um bottom layer");
  rect("bulk","-depth","depth-thermal_depth");
  m.component("comp1").geom("geom1").feature("fin").set("repairtoltype","absolute");
  m.component("comp1").geom("geom1").feature("fin").set("absrepairtol","1e-12");
  m.component("comp1").geom("geom1").run();
  box("top",1,"-1e-8","half_width+1e-8","-1e-8","1e-8");
  box("sym",1,"-1e-8","1e-8","-depth-1e-8","1e-8");
  box("skin",2,"-1e-8","half_width+1e-8","-thermal_depth-1e-8","1e-8");
  System.out.println("TOP="+Arrays.toString(m.component("comp1").selection("top").entities(1)));
  System.out.println("SYM="+Arrays.toString(m.component("comp1").selection("sym").entities(1)));
  if(m.component("comp1").selection("top").entities(1).length!=1 || m.component("comp1").selection("sym").entities(1).length!=10)
   throw new RuntimeException("Source/symmetry selection failed validation");
  m.component("comp1").variable().create("src");
  m.component("comp1").variable("src").label("Normalized laser and equivalent dipole");
  m.component("comp1").variable("src").set("gtime","exp(-0.5*((t-t0)/sigma_t)^2)/(sqrt(2*pi)*sigma_t)");
  m.component("comp1").variable("src").set("Htime","0.5*(1+erf((t-t0)/(sqrt(2)*sigma_t)))");
  m.component("comp1").variable("src").set("q_laser","F_abs*exp(-2*(x/w0)^2)*gtime");
  m.component("comp1").variable("src").set("tx_dipole","dipole_scale*A_guess*(4*x/w0)*exp(-2*(x/w0)^2)*Htime");
  m.component("comp1").material().create("mat1","Common");
  m.component("comp1").material("mat1").label("X65 benchmark: weld elasticity plus literature thermal data");
  m.component("comp1").material("mat1").propertyGroup("def").set("density","rho0");
  m.component("comp1").material("mat1").propertyGroup("def").set("youngsmodulus","E0");
  m.component("comp1").material("mat1").propertyGroup("def").set("poissonsratio","nu0");
  m.component("comp1").material("mat1").propertyGroup("def").set("thermalconductivity",new String[]{"k0","0","0","0","k0","0","0","0","k0"});
  m.component("comp1").material("mat1").propertyGroup("def").set("heatcapacity","Cp0");
  m.component("comp1").material("mat1").propertyGroup("def").set("thermalexpansioncoefficient",new String[]{"alpha0","0","0","0","alpha0","0","0","0","alpha0"});
  m.component("comp1").physics().create("ht","HeatTransfer","geom1");
  m.component("comp1").physics("ht").feature("init1").set("Tinit","T0");
  m.component("comp1").physics("ht").create("laser","HeatFluxBoundary",1);
  m.component("comp1").physics("ht").feature("laser").selection().named("top");
  m.component("comp1").physics("ht").feature("laser").set("q0_input","q_laser");
  m.component("comp1").physics("ht").feature("laser").label("Absorbed Gaussian laser heat flux");
  for(String tag:new String[]{"solid","dip"}){
   m.component("comp1").physics().create(tag,"SolidMechanics","geom1");
   m.component("comp1").physics(tag).label(tag.equals("solid")?"Thermoelastic reference":"Equivalent dipole comparison");
   if(tag.equals("dip")){
    m.component("comp1").physics(tag).field("displacement").field("udisp");
    m.component("comp1").physics(tag).field("displacement").component(new String[]{"ud","vd","wd"});
   }
   m.component("comp1").physics(tag).create("sym1","SymmetrySolid",1);
   m.component("comp1").physics(tag).feature("sym1").selection().named("sym");
  }
  m.component("comp1").physics("solid").feature("lemm1").create("te1","ThermalExpansion",2);
  m.component("comp1").physics("solid").feature("lemm1").feature("te1").set("minput_temperature_src","userdef");
  m.component("comp1").physics("solid").feature("lemm1").feature("te1").set("minput_temperature","T");
  m.component("comp1").physics("solid").feature("lemm1").feature("te1").set("minput_strainreferencetemperature_src","userdef");
  m.component("comp1").physics("solid").feature("lemm1").feature("te1").set("minput_strainreferencetemperature","T0");
  m.component("comp1").physics("dip").create("load1","BoundaryLoad",1);
  m.component("comp1").physics("dip").feature("load1").selection().named("top");
  m.component("comp1").physics("dip").feature("load1").set("FperArea",new String[]{"tx_dipole","0","0"});
  m.component("comp1").mesh().create("mesh1");
  m.component("comp1").mesh("mesh1").feature("size").set("custom",true);
  m.component("comp1").mesh("mesh1").feature("size").set("hmax","h_bulk");
  m.component("comp1").mesh("mesh1").feature("size").set("hmin","h_bulk/4");
  d=0;
  for(int i=0;i<9;i++){
   double h=5e-8*Math.pow(2,i);
   String map="map"+i;
   box("v"+i,1,"-1e-8","1e-8",Double.toString(-d-h-1e-8),Double.toString(-d+1e-8));
   box("h"+i,1,"-1e-8","half_width+1e-8",Double.toString(-d-1e-8),Double.toString(-d+1e-8));
   if(m.component("comp1").selection("v"+i).entities(1).length!=1 || m.component("comp1").selection("h"+i).entities(1).length!=1)
    throw new RuntimeException("Mapped strip edge selection failed: "+i);
   m.component("comp1").mesh("mesh1").create(map,"Map");
   m.component("comp1").mesh("mesh1").feature(map).selection().named("geom1_layer"+i+"_dom");
   m.component("comp1").mesh("mesh1").feature(map).create("dv","Distribution");
   m.component("comp1").mesh("mesh1").feature(map).feature("dv").selection().named("v"+i);
   m.component("comp1").mesh("mesh1").feature(map).feature("dv").set("numelem",1);
   m.component("comp1").mesh("mesh1").feature(map).create("dh","Distribution");
   m.component("comp1").mesh("mesh1").feature(map).feature("dh").selection().named("h"+i);
   m.component("comp1").mesh("mesh1").feature(map).feature("dh").set("numelem","n_strip");
   d+=h;
  }
  m.component("comp1").mesh("mesh1").create("tri","FreeTri");
  m.component("comp1").mesh("mesh1").feature("tri").selection().named("geom1_bulk_dom");
  m.component("comp1").mesh("mesh1").run();
  study("std1",true);
  study("std2",false);
  for(String sol:new String[]{"sol1","sol2"}){
   m.sol(sol).feature("t1").set("timemethod","genalpha");
   m.sol(sol).feature("t1").set("tstepsgenalpha","strict");
   m.sol(sol).feature("t1").set("initialstepgenalphaactive",true);
   m.sol(sol).feature("t1").set("initialstepgenalpha","0.1[ns]");
   m.sol(sol).feature("t1").set("maxstepconstraintgenalpha","expr");
   m.sol(sol).feature("t1").set("maxstepexpressiongenalpha","if(t<60[ns],0.5[ns],max_dt)");
   m.sol(sol).feature("t1").set("rtol",1e-4);
   for(String field:m.sol(sol).feature("v1").feature().tags()){
    if(field.equals("comp1_u")||field.equals("comp1_udisp")){
     m.sol(sol).feature("v1").feature(field).set("scalemethod","manual");
     m.sol(sol).feature("v1").feature(field).set("scaleval","1e-10");
    }
   }
  }
  String[] ax=new String[7],ay=new String[7];
  for(int i=0;i<7;i++){
   double theta=i*15*Math.PI/180;
   ax[i]=Double.toString(0.002*Math.sin(theta));
   ay[i]=Double.toString(-0.002*Math.cos(theta));
  }
  ay[6]="0";
  dataset("arcThermal","dset1",ax,ay);
  dataset("arcDipole","dset2",ax,ay);
  dataset("surfaceThermal","dset1",new String[]{"0","0.0005","0.001","0.0015","0.002"},new String[]{"0","0","0","0","0"});
  dataset("surfaceDipole","dset2",new String[]{"0","0.0005","0.001","0.0015","0.002"},new String[]{"0","0","0","0","0"});
  export("arcT","arcThermal",new String[]{"u","v","solid.u_tX","solid.u_tY"},new String[]{"m","m","m/s","m/s"},"results/thermoelastic_source_calibration/exports/calibration_thermal_arc.txt");
  export("arcD","arcDipole",new String[]{"ud","vd","solid2.u_tX","solid2.u_tY"},new String[]{"m","m","m/s","m/s"},"results/thermoelastic_source_calibration/exports/calibration_dipole_arc.txt");
  export("surfT","surfaceThermal",new String[]{"u","v","solid.u_tX","solid.u_tY","T"},new String[]{"m","m","m/s","m/s","K"},"results/thermoelastic_source_calibration/exports/calibration_thermal_surface.txt");
  export("surfD","surfaceDipole",new String[]{"ud","vd","solid2.u_tX","solid2.u_tY"},new String[]{"m","m","m/s","m/s"},"results/thermoelastic_source_calibration/exports/calibration_dipole_surface.txt");
  for(String tag:new String[]{"pgT","pgU","pgD"})m.result().create(tag,"PlotGroup2D");
  m.result("pgT").label("Temperature rise near laser (Study 1)");
  m.result("pgT").set("data","dset1");
  m.result("pgT").create("surf1","Surface");
  m.result("pgT").feature("surf1").set("expr","T-T0");
  m.result("pgU").label("Vertical displacement - thermoelastic");
  m.result("pgU").set("data","dset1");
  m.result("pgU").create("surf1","Surface");
  m.result("pgU").feature("surf1").set("expr","v");
  m.result("pgD").label("Vertical displacement - equivalent dipole");
  m.result("pgD").set("data","dset2");
  m.result("pgD").create("surf1","Surface");
  m.result("pgD").feature("surf1").set("expr","vd");
  try {m.save("models/laser_source_calibration.mph");}catch(Exception e){throw new RuntimeException(e);}
  System.out.println("CALIBRATION_BUILD_OK");
  return m;
 }
 public static void main(String[] args){run();}
}
