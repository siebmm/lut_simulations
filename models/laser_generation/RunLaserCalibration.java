import com.comsol.model.*;
import com.comsol.model.util.*;
import java.nio.charset.StandardCharsets;
import java.nio.file.*;
import java.util.*;

/** Solve the original local laser benchmark and export the quantities used in its calibration.
 * Adapted from the NDT study's RunLaserCalibration and CheckLaserCalibration.
 * All writes are relative to a separate run directory, never the original model.
 */
public class RunLaserCalibration {
    static final Path WORK = Paths.get("").toAbsolutePath();
    static final Path EXPORTS = WORK.resolve("exports");

    static void parameters(Model m, boolean complete) throws Exception {
        String[][] names = {
            {"absorbed_fluence_J_m2", "F_abs"}, {"pulse_fwhm_s", "pulse_fwhm"},
            {"pulse_centre_s", "t0"}, {"line_half_width_m", "w0"},
            {"rho_kg_m3", "rho0"}, {"cp_J_kg_K", "Cp0"}, {"k_W_m_K", "k0"},
            {"ambient_K", "T0"}, {"E_Pa", "E0"}, {"nu", "nu0"},
            {"alpha_K_inv", "alpha0"}, {"dipole_scale", "dipole_scale"}
        };
        StringBuilder json = new StringBuilder("{\n  \"complete\": " + complete);
        for (String[] entry : names)
            json.append(",\n  \"").append(entry[0]).append("\": ").append(m.param().evaluate(entry[1]));
        json.append("\n}\n");
        Files.write(EXPORTS.resolve("calibration_parameters.json"), json.toString().getBytes(StandardCharsets.UTF_8));
    }

    static void configureExports(Model m) {
        String[] tags = {"arcT", "arcD", "surfT", "surfD"};
        String[] names = {"thermal_arc", "dipole_arc", "thermal_surface", "dipole_surface"};
        String[][] expressions = {
            {"u", "v", "solid.u_tX", "solid.u_tY"},
            {"ud", "vd", "solid2.u_tX", "solid2.u_tY"},
            {"u", "v", "solid.u_tX", "solid.u_tY", "T"},
            {"ud", "vd", "solid2.u_tX", "solid2.u_tY"}
        };
        for (int i = 0; i < tags.length; i++) {
            m.result().export(tags[i]).set("filename", EXPORTS.resolve("calibration_" + names[i] + ".txt").toString());
            // Older saved models used d(u,t), which cannot be evaluated on cut-point datasets.
            m.result().export(tags[i]).set("expr", expressions[i]);
            m.result().export(tags[i]).set("unit", i == 2 ?
                new String[]{"m", "m", "m/s", "m/s", "K"} : new String[]{"m", "m", "m/s", "m/s"});
            m.result().export(tags[i]).set("header", true);
            m.result().export(tags[i]).set("fullprec", true);
        }
    }

    static void energyCheck(Model m) throws Exception {
        if (Arrays.asList(m.result().numerical().tags()).contains("energyCheck"))
            m.result().numerical().remove("energyCheck");
        m.result().numerical().create("energyCheck", "IntSurface");
        m.result().numerical("energyCheck").set("data", "dset1");
        m.result().numerical("energyCheck").selection().all();
        m.result().numerical("energyCheck").set("expr", new String[]{
            "2*rho0*Cp0*(T-T0)", "2*rho0*Cp0*(T-T0)/Eline_abs"});
        m.result().numerical("energyCheck").set("unit", new String[]{"J/m", "1"});
        if (!Arrays.asList(m.result().table().tags()).contains("tblEnergy"))
            m.result().table().create("tblEnergy", "Table");
        m.result().numerical("energyCheck").set("table", "tblEnergy");
        m.result().numerical("energyCheck").setResult();
        m.result().table("tblEnergy").save(EXPORTS.resolve("calibration_energy_check.txt").toString());
    }

    public static Model run() throws Exception {
        Properties settings = new Properties();
        try (java.io.Reader reader = Files.newBufferedReader(WORK.resolve("run.properties"))) {
            settings.load(reader);
        }
        double fluence = Double.parseDouble(settings.getProperty("fluence_J_m2"));
        if (!Double.isFinite(fluence) || fluence <= 0) throw new IllegalArgumentException("Fluence must be positive");
        Files.createDirectories(EXPORTS);
        Model m = ModelUtil.load("Calibration", WORK.resolve("models/laser_source_calibration.mph").toString());
        m.modelPath(WORK.toString());
        m.param().set("F_abs", Double.toString(fluence) + "[J/m^2]");
        m.param().set("dipole_scale", "1"); // Fit the multiplier after solving both source descriptions.
        configureExports(m);
        for (String sol : new String[]{"sol1", "sol2"}) {
            m.sol(sol).feature("t1").set("maxstepconstraintgenalpha", "expr");
            m.sol(sol).feature("t1").set("maxstepexpressiongenalpha", "if(t<60[ns],0.5[ns],max_dt)");
            for (String field : m.sol(sol).feature("v1").feature().tags()) {
                if (field.equals("comp1_u") || field.equals("comp1_udisp")) {
                    m.sol(sol).feature("v1").feature(field).set("scalemethod", "manual");
                    m.sol(sol).feature("v1").feature(field).set("scaleval", "1e-10");
                }
            }
        }
        parameters(m, false);
        m.save(WORK.resolve("models/laser_source_calibration.mph").toString());
        if (Boolean.parseBoolean(settings.getProperty("prepare_only", "false"))) {
            System.out.println("LASER_MODEL_PREPARED");
            return m;
        }
        System.out.println("Solving heat conduction and thermoelastic wave generation");
        m.sol("sol1").runAll();
        m.save(WORK.resolve("models/laser_source_calibration_solved.mph").toString());
        m.result().export("arcT").run();
        m.result().export("surfT").run();
        System.out.println("Solving the equivalent elastic dipole");
        m.sol("sol2").runAll();
        m.result().export("arcD").run();
        m.result().export("surfD").run();
        energyCheck(m);
        m.save(WORK.resolve("models/laser_source_calibration_solved.mph").toString());
        parameters(m, true);
        System.out.println("LASER_CALIBRATION_COMPLETE");
        return m;
    }

    public static void main(String[] args) throws Exception { run(); }
}
