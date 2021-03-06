package danfulea.phys.resources;

import java.util.ListResourceBundle;

/**
 * Resources for DoseSimCore
 * 
 * @author Dan Fulea, 10 APR. 2005
 */

public class DoseSimCoreResources extends ListResourceBundle {
	public Object[][] getContents() {
		return CONTENTS;
	}

	static final Object[][] CONTENTS = {

			// Coeficienti de atenuare pentru tesuturi in Rx energy range
			// <=0.150MeV!!
			// based on XCOM v.3.1 software developped by NIST (National
			// Institute of Standards
			// and Technology)-->Martin J. Berger,Stephen M. Seltzer and John H.
			// Hubbell
			// last update of software and database:June 1999
			// Partial Interaction Coefficients and Total Attenuation
			// Coefficients
			// energy->MeV and coefficients ->cm2/g. BASED ON MIRD MATH.
			// PHANTOM-TISSUE COMPOSITION
			{
					"commonTissue.diagnostic.attenuationCoef",
					new double[][]
					// energy, scatt_coh, scatt_incoh, PhotoelAbs,
					// Total(+scatter_coh),Total(-scatter_coh)
					{
							{ 1.000E-03, 1.292E+00, 1.390E-02, 3.521E+03,
									3.523E+03, 3.521E+03 },
							{ 2.000E-03, 1.064E+00, 4.342E-02, 5.273E+02,
									5.284E+02, 5.273E+02 },
							{ 3.000E-03, 8.299E-01, 7.259E-02, 1.689E+02,
									1.698E+02, 1.690E+02 },
							{ 4.000E-03, 6.430E-01, 9.593E-02, 7.393E+01,
									7.467E+01, 7.403E+01 },
							{ 5.000E-03, 5.063E-01, 1.135E-01, 3.806E+01,
									3.868E+01, 3.817E+01 },
							{ 6.000E-03, 4.085E-01, 1.265E-01, 2.192E+01,
									2.245E+01, 2.204E+01 },
							{ 8.000E-03, 2.848E-01, 1.438E-01, 9.090E+00,
									9.518E+00, 9.233E+00 },
							{ 1.000E-02, 2.135E-01, 1.545E-01, 4.555E+00,
									4.923E+00, 4.709E+00 },
							{ 1.500E-02, 1.245E-01, 1.694E-01, 1.271E+00,
									1.565E+00, 1.441E+00 },
							{ 2.000E-02, 8.265E-02, 1.770E-01, 5.090E-01,
									7.686E-01, 6.860E-01 },
							{ 3.000E-02, 4.365E-02, 1.823E-01, 1.376E-01,
									3.635E-01, 3.199E-01 },
							{ 4.000E-02, 2.671E-02, 1.820E-01, 5.395E-02,
									2.626E-01, 2.359E-01 },
							{ 5.000E-02, 1.799E-02, 1.794E-01, 2.600E-02,
									2.234E-01, 2.054E-01 },
							{ 6.000E-02, 1.293E-02, 1.761E-01, 1.430E-02,
									2.033E-01, 1.904E-01 },
							{ 8.000E-02, 7.577E-03, 1.687E-01, 5.557E-03,
									1.818E-01, 1.743E-01 },
							{ 1.000E-01, 4.961E-03, 1.616E-01, 2.672E-03,
									1.692E-01, 1.643E-01 },
							{ 1.500E-01, 2.263E-03, 1.464E-01, 7.117E-04,
									1.494E-01, 1.471E-01 }, } },
			// BASED ON MIRD MATH. PHANTOM-TISSUE COMPOSITION
			{ "skeletal_boneWithMarrow.diagnostic.attenuationCoef",
					new double[][]
					// energy, scatt_coh, scatt_incoh, PhotoelAbs,
					// Total(+scatter_coh),Total(-scatter_coh)
					{
							{ 1.000E-03, 1.584E+00, 1.352E-02, 3.475E+03,
									3.477E+03, 3.475E+03 },
							{ 2.000E-03, 1.300E+00, 4.036E-02, 5.303E+02,
									5.316E+02, 5.303E+02 },
							{ 3.000E-03, 1.030E+00, 6.661E-02, 2.220E+02,
									2.231E+02, 2.221E+02 },
							{ 4.000E-03, 8.164E-01, 8.782E-02, 9.909E+01,
									9.999E+01, 9.917E+01 },
							{ 5.000E-03, 6.584E-01, 1.039E-01, 1.066E+02,
									1.074E+02, 1.067E+02 },
							{ 6.000E-03, 5.433E-01, 1.161E-01, 6.437E+01,
									6.503E+01, 6.448E+01 },
							{ 8.000E-03, 3.933E-01, 1.327E-01, 2.869E+01,
									2.922E+01, 2.883E+01 },
							{ 1.000E-02, 3.027E-01, 1.435E-01, 1.513E+01,
									1.557E+01, 1.527E+01 },
							{ 1.500E-02, 1.815E-01, 1.591E-01, 4.601E+00,
									4.942E+00, 4.761E+00 },
							{ 2.000E-02, 1.209E-01, 1.674E-01, 1.946E+00,
									2.234E+00, 2.113E+00 },
							{ 3.000E-02, 6.429E-02, 1.738E-01, 5.645E-01,
									8.026E-01, 7.383E-01 },
							{ 4.000E-02, 3.989E-02, 1.742E-01, 2.315E-01,
									4.456E-01, 4.057E-01 },
							{ 5.000E-02, 2.719E-02, 1.722E-01, 1.152E-01,
									3.147E-01, 2.875E-01 },
							{ 6.000E-02, 1.972E-02, 1.694E-01, 6.492E-02,
									2.540E-01, 2.343E-01 },
							{ 8.000E-02, 1.169E-02, 1.627E-01, 2.613E-02,
									2.005E-01, 1.889E-01 },
							{ 1.000E-01, 7.710E-03, 1.561E-01, 1.292E-02,
									1.767E-01, 1.690E-01 },
							{ 1.500E-01, 3.557E-03, 1.417E-01, 3.580E-03,
									1.488E-01, 1.453E-01 }, } },
			// BASED ON MIRD MATH. PHANTOM-TISSUE COMPOSITION
			{
					"lung.diagnostic.attenuationCoef",
					new double[][]
					// energy, scatt_coh, scatt_incoh, PhotoelAbs,
					// Total(+scatter_coh),Total(-scatter_coh)
					{
							{ 1.000E-03, 1.349E+00, 1.327E-02, 3.825E+03,
									3.826E+03, 3.825E+03 },
							{ 2.000E-03, 1.122E+00, 4.179E-02, 5.775E+02,
									5.786E+02, 5.775E+02 },
							{ 3.000E-03, 8.817E-01, 7.040E-02, 1.865E+02,
									1.875E+02, 1.866E+02 },
							{ 4.000E-03, 6.859E-01, 9.359E-02, 8.156E+01,
									8.234E+01, 8.165E+01 },
							{ 5.000E-03, 5.407E-01, 1.112E-01, 4.195E+01,
									4.261E+01, 4.206E+01 },
							{ 6.000E-03, 4.361E-01, 1.245E-01, 2.418E+01,
									2.474E+01, 2.431E+01 },
							{ 8.000E-03, 3.031E-01, 1.421E-01, 1.013E+01,
									1.057E+01, 1.027E+01 },
							{ 1.000E-02, 2.264E-01, 1.530E-01, 5.081E+00,
									5.460E+00, 5.234E+00 },
							{ 1.500E-02, 1.315E-01, 1.680E-01, 1.423E+00,
									1.722E+00, 1.591E+00 },
							{ 2.000E-02, 8.738E-02, 1.756E-01, 5.700E-01,
									8.330E-01, 7.456E-01 },
							{ 3.000E-02, 4.625E-02, 1.811E-01, 1.544E-01,
									3.818E-01, 3.356E-01 },
							{ 4.000E-02, 2.833E-02, 1.810E-01, 6.064E-02,
									2.700E-01, 2.416E-01 },
							{ 5.000E-02, 1.909E-02, 1.786E-01, 2.926E-02,
									2.269E-01, 2.079E-01 },
							{ 6.000E-02, 1.373E-02, 1.753E-01, 1.610E-02,
									2.052E-01, 1.914E-01 },
							{ 8.000E-02, 8.052E-03, 1.681E-01, 6.266E-03,
									1.824E-01, 1.743E-01 },
							{ 1.000E-01, 5.276E-03, 1.610E-01, 3.015E-03,
									1.693E-01, 1.641E-01 },
							{ 1.500E-01, 2.409E-03, 1.459E-01, 8.040E-04,
									1.491E-01, 1.467E-01 }, } },
			// ---------------MIRD NEWBORN
			// ---------------------------------------
			{
					"commonTissue.newborn.diagnostic.attenuationCoef",
					new double[][]
					// energy, scatt_coh, scatt_incoh, PhotoelAbs,
					// Total(+scatter_coh),Total(-scatter_coh)
					{
							{ 1.000E-03, 1.325E+00, 1.364E-02, 3.708E+03,
									3.709E+03, 3.708E+03 },
							{ 2.000E-03, 1.098E+00, 4.279E-02, 5.569E+02,
									5.580E+02, 5.569E+02 },
							{ 3.000E-03, 8.606E-01, 7.180E-02, 1.786E+02,
									1.795E+02, 1.787E+02 },
							{ 4.000E-03, 6.685E-01, 9.517E-02, 7.892E+01,
									7.969E+01, 7.902E+01 },
							{ 5.000E-03, 5.267E-01, 1.128E-01, 4.057E+01,
									4.121E+01, 4.068E+01 },
							{ 6.000E-03, 4.248E-01, 1.260E-01, 2.338E+01,
									2.393E+01, 2.351E+01 },
							{ 8.000E-03, 2.956E-01, 1.436E-01, 9.708E+00,
									1.015E+01, 9.851E+00 },
							{ 1.000E-02, 2.210E-01, 1.544E-01, 4.863E+00,
									5.238E+00, 5.017E+00 },
							{ 1.500E-02, 1.286E-01, 1.694E-01, 1.359E+00,
									1.657E+00, 1.528E+00 },
							{ 2.000E-02, 8.538E-02, 1.769E-01, 5.433E-01,
									8.055E-01, 7.201E-01 },
							{ 3.000E-02, 4.515E-02, 1.823E-01, 1.469E-01,
									3.743E-01, 3.292E-01 },
							{ 4.000E-02, 2.765E-02, 1.820E-01, 5.760E-02,
									2.673E-01, 2.396E-01 },
							{ 5.000E-02, 1.863E-02, 1.796E-01, 2.776E-02,
									2.259E-01, 2.073E-01 },
							{ 6.000E-02, 1.339E-02, 1.762E-01, 1.526E-02,
									2.049E-01, 1.915E-01 },
							{ 8.000E-02, 7.853E-03, 1.689E-01, 5.931E-03,
									1.827E-01, 1.748E-01 },
							{ 1.000E-01, 5.143E-03, 1.618E-01, 2.851E-03,
									1.698E-01, 1.647E-01 },
							{ 1.500E-01, 2.348E-03, 1.466E-01, 7.591E-04,
									1.497E-01, 1.473E-01 }, } },

			{ "skeletal_boneWithMarrow.newborn.diagnostic.attenuationCoef",
					new double[][]
					// energy, scatt_coh, scatt_incoh, PhotoelAbs,
					// Total(+scatter_coh),Total(-scatter_coh)
					{
							{ 1.000E-03, 1.584E+00, 1.297E-02, 3.852E+03,
									3.853E+03, 3.852E+03 },
							{ 2.000E-03, 1.313E+00, 3.945E-02, 5.910E+02,
									5.924E+02, 5.910E+02 },
							{ 3.000E-03, 1.045E+00, 6.586E-02, 2.287E+02,
									2.298E+02, 2.288E+02 },
							{ 4.000E-03, 8.272E-01, 8.749E-02, 1.014E+02,
									1.023E+02, 1.015E+02 },
							{ 5.000E-03, 6.641E-01, 1.042E-01, 9.580E+01,
									9.656E+01, 9.590E+01 },
							{ 6.000E-03, 5.447E-01, 1.168E-01, 5.750E+01,
									5.816E+01, 5.762E+01 },
							{ 8.000E-03, 3.897E-01, 1.340E-01, 2.542E+01,
									2.594E+01, 2.555E+01 },
							{ 1.000E-02, 2.972E-01, 1.450E-01, 1.331E+01,
									1.376E+01, 1.346E+01 },
							{ 1.500E-02, 1.765E-01, 1.605E-01, 4.012E+00,
									4.349E+00, 4.173E+00 },
							{ 2.000E-02, 1.175E-01, 1.687E-01, 1.685E+00,
									1.971E+00, 1.853E+00 },
							{ 3.000E-02, 6.256E-02, 1.750E-01, 4.852E-01,
									7.228E-01, 6.602E-01 },
							{ 4.000E-02, 3.873E-02, 1.754E-01, 1.981E-01,
									4.123E-01, 3.736E-01 },
							{ 5.000E-02, 2.635E-02, 1.735E-01, 9.829E-02,
									2.981E-01, 2.718E-01 },
							{ 6.000E-02, 1.908E-02, 1.706E-01, 5.526E-02,
									2.449E-01, 2.258E-01 },
							{ 8.000E-02, 1.130E-02, 1.639E-01, 2.217E-02,
									1.973E-01, 1.860E-01 },
							{ 1.000E-01, 7.443E-03, 1.572E-01, 1.090E-02,
									1.755E-01, 1.681E-01 },
							{ 1.500E-01, 3.428E-03, 1.427E-01, 3.007E-03,
									1.491E-01, 1.457E-01 }, } },

			{
					"lung.newborn.diagnostic.attenuationCoef",
					new double[][]
					// energy, scatt_coh, scatt_incoh, PhotoelAbs,
					// Total(+scatter_coh),Total(-scatter_coh)
					{
							{ 1.000E-03, 1.349E+00, 1.327E-02, 3.825E+03,
									3.827E+03, 3.825E+03 },
							{ 2.000E-03, 1.121E+00, 4.179E-02, 5.773E+02,
									5.785E+02, 5.773E+02 },
							{ 3.000E-03, 8.816E-01, 7.040E-02, 1.864E+02,
									1.874E+02, 1.865E+02 },
							{ 4.000E-03, 6.858E-01, 9.360E-02, 8.152E+01,
									8.230E+01, 8.162E+01 },
							{ 5.000E-03, 5.406E-01, 1.112E-01, 4.194E+01,
									4.259E+01, 4.205E+01 },
							{ 6.000E-03, 4.360E-01, 1.245E-01, 2.417E+01,
									2.473E+01, 2.430E+01 },
							{ 8.000E-03, 3.030E-01, 1.421E-01, 1.012E+01,
									1.057E+01, 1.027E+01 },
							{ 1.000E-02, 2.263E-01, 1.530E-01, 5.076E+00,
									5.455E+00, 5.229E+00 },
							{ 1.500E-02, 1.315E-01, 1.680E-01, 1.421E+00,
									1.721E+00, 1.589E+00 },
							{ 2.000E-02, 8.736E-02, 1.756E-01, 5.688E-01,
									8.318E-01, 7.444E-01 },
							{ 3.000E-02, 4.624E-02, 1.812E-01, 1.541E-01,
									3.815E-01, 3.352E-01 },
							{ 4.000E-02, 2.832E-02, 1.810E-01, 6.048E-02,
									2.698E-01, 2.415E-01 },
							{ 5.000E-02, 1.908E-02, 1.786E-01, 2.918E-02,
									2.269E-01, 2.078E-01 },
							{ 6.000E-02, 1.372E-02, 1.753E-01, 1.605E-02,
									2.051E-01, 1.914E-01 },
							{ 8.000E-02, 8.050E-03, 1.681E-01, 6.245E-03,
									1.824E-01, 1.743E-01 },
							{ 1.000E-01, 5.274E-03, 1.610E-01, 3.004E-03,
									1.693E-01, 1.640E-01 },
							{ 1.500E-01, 2.408E-03, 1.459E-01, 8.010E-04,
									1.491E-01, 1.467E-01 }, } },
			// -------------------------------------------------------------------
			// NRPB-1985--composition data
			{
					"breast.diagnostic.attenuationCoef",
					new double[][]
					// energy, scatt_coh, scatt_incoh, PhotoelAbs,
					// Total(+scatter_coh),Total(-scatter_coh)
					{
							{ 1.000E-03, 1.205E+00, 1.496E-02, 3.144E+03,
									3.145E+03, 3.144E+03 },
							{ 2.000E-03, 9.797E-01, 4.629E-02, 4.629E+02,
									4.639E+02, 4.630E+02 },
							{ 3.000E-03, 7.557E-01, 7.662E-02, 1.425E+02,
									1.433E+02, 1.426E+02 },
							{ 4.000E-03, 5.815E-01, 1.004E-01, 6.048E+01,
									6.116E+01, 6.058E+01 },
							{ 5.000E-03, 4.564E-01, 1.180E-01, 3.078E+01,
									3.135E+01, 3.090E+01 },
							{ 6.000E-03, 3.680E-01, 1.308E-01, 1.761E+01,
									1.811E+01, 1.774E+01 },
							{ 8.000E-03, 2.571E-01, 1.478E-01, 7.216E+00,
									7.621E+00, 7.364E+00 },
							{ 1.000E-02, 1.933E-01, 1.583E-01, 3.582E+00,
									3.934E+00, 3.741E+00 },
							{ 1.500E-02, 1.132E-01, 1.731E-01, 9.859E-01,
									1.272E+00, 1.159E+00 },
							{ 2.000E-02, 7.508E-02, 1.804E-01, 3.900E-01,
									6.455E-01, 5.704E-01 },
							{ 3.000E-02, 3.950E-02, 1.854E-01, 1.040E-01,
									3.289E-01, 2.894E-01 },
							{ 4.000E-02, 2.412E-02, 1.847E-01, 4.042E-02,
									2.492E-01, 2.251E-01 },
							{ 5.000E-02, 1.622E-02, 1.819E-01, 1.935E-02,
									2.175E-01, 2.013E-01 },
							{ 6.000E-02, 1.164E-02, 1.784E-01, 1.059E-02,
									2.006E-01, 1.890E-01 },
							{ 8.000E-02, 6.812E-03, 1.708E-01, 4.084E-03,
									1.817E-01, 1.749E-01 },
							{ 1.000E-01, 4.454E-03, 1.635E-01, 1.952E-03,
									1.699E-01, 1.655E-01 },
							{ 1.500E-01, 2.028E-03, 1.480E-01, 5.156E-04,
									1.506E-01, 1.486E-01 }, } },
			// ROSENSTEIN BUT NEVER USED!!!
			{
					"activeBoneMarrow.diagnostic.attenuationCoef",
					new double[][]
					// energy, scatt_coh, scatt_incoh, PhotoelAbs,
					// Total(+scatter_coh),Total(-scatter_coh)
					{
							{ 1.000E-03, 1.174E+00, 1.503E-02, 2.896E+03,
									2.897E+03, 2.896E+03 },
							{ 2.000E-03, 9.448E-01, 4.625E-02, 4.216E+02,
									4.226E+02, 4.217E+02 },
							{ 3.000E-03, 7.227E-01, 7.630E-02, 1.327E+02,
									1.335E+02, 1.328E+02 },
							{ 4.000E-03, 5.539E-01, 9.979E-02, 5.781E+01,
									5.846E+01, 5.791E+01 },
							{ 5.000E-03, 4.347E-01, 1.170E-01, 2.952E+01,
									3.007E+01, 2.964E+01 },
							{ 6.000E-03, 3.511E-01, 1.295E-01, 1.694E+01,
									1.742E+01, 1.707E+01 },
							{ 8.000E-03, 2.469E-01, 1.460E-01, 7.001E+00,
									7.394E+00, 7.147E+00 },
							{ 1.000E-02, 1.867E-01, 1.564E-01, 3.490E+00,
									3.833E+00, 3.647E+00 },
							{ 1.500E-02, 1.100E-01, 1.712E-01, 9.680E-01,
									1.249E+00, 1.139E+00 },
							{ 2.000E-02, 7.294E-02, 1.787E-01, 3.851E-01,
									6.368E-01, 5.638E-01 },
							{ 3.000E-02, 3.829E-02, 1.837E-01, 1.036E-01,
									3.255E-01, 2.872E-01 },
							{ 4.000E-02, 2.338E-02, 1.830E-01, 4.047E-02,
									2.468E-01, 2.234E-01 },
							{ 5.000E-02, 1.572E-02, 1.802E-01, 1.946E-02,
									2.154E-01, 1.997E-01 },
							{ 6.000E-02, 1.128E-02, 1.767E-01, 1.068E-02,
									1.986E-01, 1.874E-01 },
							{ 8.000E-02, 6.598E-03, 1.692E-01, 4.141E-03,
									1.799E-01, 1.733E-01 },
							{ 1.000E-01, 4.313E-03, 1.619E-01, 1.987E-03,
									1.682E-01, 1.639E-01 },
							{ 1.500E-01, 1.963E-03, 1.466E-01, 5.282E-04,
									1.491E-01, 1.471E-01 }, } },

			{ "commonTissue.diagnostic.density", new Double(1.04) },// g/cm3-MIRD
			{ "skeletal_boneWithMarrow.diagnostic.density", new Double(1.4) },// g/cm3-MIRD
			{ "activeBoneMarrow.diagnostic.density", new Double(1.25052223) },// g/cm3-from
																				// composition
			{ "lung.diagnostic.density", new Double(0.296) },// g/cm3-MIRD
			{ "breast.diagnostic.density", new Double(1.006808731) },// g/cm3-from
																		// composition
			{ "commonTissue.newborn.diagnostic.density", new Double(1.04) },// g/cm3-MIRD
			{ "skeletal_boneWithMarrow.newborn.diagnostic.density",
					new Double(1.22) },// g/cm3-MIRD
			{ "lung.newborn.diagnostic.density", new Double(0.296) },// g/cm3-MIRD

			// factori de pondere tisulari doza(mGy)->doza(mSv)
			{ "breast.wt", new Double(0.05) },

			// organs
			{ "mamo.organname", new String[] { "San" } },
			{
					"rad.organname",
					new String[] { "san", "testicule", "schelet",
							"maduva activa", "glanda suprarenala", "creier",
							"vezica biliara", "stomac", "intestin subtire",
							"intestin gros sup.", "intestin gros inf.",
							"inima", "rinichi", "ficat", "plamani", "ovare",
							"pancreas", "splina", "timus", "tiroida",
							"vezica urinara", "uter", "restul" } },
			{
					"wt.organname",
					new double[] { 0.05, 0.20, 0.13, 0.12, 0.05, 0.05, 0.05,
							0.05, 0.12, 0.12, 0.12, 0.05, 0.05, 0.05, 0.12,
							0.20, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05 } },
			// -------------Separated Active Bone Marrow-supposetobe the same
			// for all phantoms
			{
					"abm.sk.fraction",
					new double[][] {
							// MeV, raport(coeff abs!!)
							{ 1.00E-03, 0.82932417 },
							{ 2.00E-03, 0.790885221 },
							{ 3.00E-03, 0.549214227 },
							{ 4.00E-03, 0.534225092 },
							{ 5.00E-03, 0.268965517 },
							{ 6.00E-03, 0.256383298 },
							{ 8.00E-03, 0.240396905 },
							{ 1.00E-02, 0.231849968 },
							{ 1.50E-02, 0.232876712 },
							{ 2.00E-02, 0.260656496 },
							{ 3.00E-02, 0.382270731 },
							{ 4.00E-02, 0.544346979 },
							{ 5.00E-02, 0.690286899 },
							{ 6.00E-02, 0.797446809 },
							{ 8.00E-02, 0.917902542 },
							{ 1.00E-01, 0.971547125 },
							{ 1.50E-01, 1.015883978 }, } },
			// MasaABM/MasaSkelet
			{
					"abm.weight",
					new double[] {
							// newBorn, age1, age 5,age10, age15, adult
							0.133903134, 0.131578947, 0.118081181, 0.13174946,
							0.137254902, 0.112 } },
			// -------------------------------------------
			// for r=0-0.39
			{
					"kn.wielopolski.r1",
					new double[][]
					// B1, B2, B3, B4., B5, B6, B7
					{
							{ -0.01238, -1.5687, -0.42888, -28.63753, 31.96466,
									6.46298, 0 },
							{ 0.02902, 3.8898, 3.73901, 65.38702, -70.51987,
									-8.63537, 0 },
							{ 0.00185, 0.2423, 0.20125, 3.75526, -4.0744,
									0.83815, 0 },
							{ -0.000085, -0.0113, -0.01234, -0.14929, 0.16095,
									0.03104, 0 }, } },

			// for r=0.39-0.7
			{
					"kn.wielopolski.r2",
					new double[][]
					// B1, B2, B3, B4., B5, B6, B7
					{
							{ 5.78087, -22.76518, 52.21054, -29.24388,
									-2.30266, 19.38525, 0.6351 },
							{ -11.74909, 58.75565, -114.7284, 60.70503,
									5.64759, 19.24324, 0.64313 },
							{ -2.28752, 8.55102, -56.57779, 44.52396, 8.55989,
									4.68978, 0.7508 },
							{ 0.03073, -0.16078, 0.47973, -0.34037, -0.03003,
									12.27394, 0.65229 }, } },

			// for r=0.7-1.0
			{
					"kn.wielopolski.r3",
					new double[][]
					// B1, B2, B3, B4., B5, B6, B7
					{
							{ 3.12358, -5.62131, -51.9096, -20.89746,
									191.18823, 10.07783, 0 },
							{ 0.04248, 13.69255, 127.2337, -1.01555, -411.95,
									-23.63232, 0 },
							{ 0.001942, 0.5932, 5.73861, 6.38184, -28.01964,
									-1.71261, 0 },
							{ -0.000045, -0.013426, -0.15548, -0.45873,
									1.31546, 0.04947, 0 }, } },

	// -------------------------------------------
	};
}
