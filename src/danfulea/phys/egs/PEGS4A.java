package danfulea.phys.egs;

import java.io.FileWriter;
import java.util.Calendar;
import java.util.Vector;

/**
 * This class is designed to initialize and perform PEGS specific tasks in order to prepare medium input data for EGS (Electron Gamma Shower). 
 * Based on work done by SLAC (Stanford Linear Accelerator Center) and EGSnrc developers.
 * @author Dan Fulea, 16 AUG. 2005 
 * 
 */

public class PEGS4A {
	private static String medium = "";
	public static boolean mediumB = true;
	public static String mediumFilename = "";
	// private static int imixt=0;
	private static int ne = 0;
	private static double rho = 0.;
	private static String[] asymetry;
	private static double[] rhoze;
	private static double[] pze;
	private static boolean densityB = false;
	private static boolean densityFullB = false;
	private static double gasPressure = 0.0;

	private static String logFile = "pegs4.log";
	private static String intDir = "interactiv";
	private static String extension = ".pegs4dat";
	private static String phpData = "pgs4pepr.dat";// default
	private static String formData = "pgs4form.dat";// default
	private static String densData = "";
	private static boolean densB = false;// true if is compound or mixture,
											// false if is element
	private static boolean cuttB = false;// flag to: true if AE,AP...is set
	private static double g24 = 1.6605402;// g24*10**-24 is u.a.m
	public static String func = "";
	public static String funci = "";
	public static String funceval = "";
	// -----------------------------------------------------
	public static double REM = 0.;
	public static int NJj = 0;
	public static int NKP = 3;
	// -------------------------------------------------------
	public static int $MXEKE = PEGS4.$MXEKE;// 300;//300-WORKS
	public static int NEL = 0;
	public static int NALE = $MXEKE;
	public static double EPE = 0.01;// 1% error
	// public static double EPE=0.05;//5% error=WORKS
	public static double[] ZTHRE = { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
			0., 0., 0., 0., 0., 0., };
	public static double[] ZEPE = { 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,
			0., 0., 0., 0., 0., };
	public static int NIPE = 20;
	public static double AXE = 0.;//
	public static double BXE = 0.;//
	public static double[][] AFE = new double[$MXEKE][8];
	public static double[][] BFE = new double[$MXEKE][8];

	public static int $MXGE = PEGS4.$MXGE;// 300;//NOT WORKS because Photte is
											// with Interp->could lead to
											// multiple singularities
	public static int NGL = 0;
	public static int NALG = $MXGE;
	public static double EPG = 0.01;
	// public static double EPG=0.08;//8% error
	public static double[] ZTHRG = { 0.0, 0.1, 0.0, 0.0 };
	public static double[] ZEPG = { 0.0, 0.01, 0.0, 0.0 };
	public static int NIPG = 20;
	public static double AXG = 0.;//
	public static double BXG = 0.;//
	public static double[][] AFG = new double[$MXGE][4];
	public static double[][] BFG = new double[$MXGE][4];

	public static int $MXRL = PEGS4.$MXRL;// 300;//300-WORKS;
	public static int NGR = 0;
	public static int NALR = $MXRL;
	public static double EPR = 0.01;
	// public static double EPR=0.08;//8% error,WORKS
	public static double[] ZTHRR = { 0.0, 0.0 };
	public static double[] ZEPR = { 0.0, 0.0 };
	public static int NIPR = 20;
	public static double AXR = 0.;//
	public static double BXR = 0.;//
	public static double[][] AFR = new double[$MXRL][1];
	public static double[][] BFR = new double[$MXRL][1];
	// -------------plot---------------------------------
	public static int interv = 200;
	public static double[] energ;
	public static double[] photo;
	public static double[] compt;
	public static double[] pair;
	public static double[] cohe;
	public static double[] totg;
	// --------------ep
	public static int interve = 200;
	public static double[] energe;
	public static double[] brem;
	public static double[] moll;
	public static double[] tote;
	public static double[] bhab;
	public static double[] anni;
	public static double[] totp;

	// set Photoelectric and Pair production data file
	/**
	 * Set Photoelectric and Pair production data file
	 * @param s s
	 */
	public static void setPhpData(String s) {
		phpData = s;
	}

	// set form factor data file
	/**
	 * Set form factor data file
	 * @param s s
	 */
	public static void setFormData(String s) {
		formData = s;
	}

	// set density data file
	/**
	 * Set density data file for mixture (true) or element (false)
	 * @param s s
	 * @param b b
	 */
	public static void setDensData(String s, boolean b) {
		densData = s;
		densB = b;
		densityB = true;
	}

	// set density data file
	/**
	 * Set density data file knowing its full path.
	 * @param s s
	 */
	public static void setDensData(String s) {
		densData = s;
		densityB = true;
		densityFullB = true;
	}

	/**
	 * Set gas pressure
	 * @param d d
	 */
	public static void setGasp(double d) {
		gasPressure = d;
	}

	/**
	 * Reset global variables for re-use
	 */
	public static void reset() {
		densData = "";
		densityB = false;
		densityFullB = false;
		gasPressure = 0.0;
		// PEGS4 global variabiles
		PEGS4.zdepv = new Vector<String>();
		PEGS4.npairdelta = 0;
		PEGS4.jzet = 0;
		PEGS4.jrhozet = 0;
		PEGS4.jenergy = 0;
		PEGS4.jdelta = 0;
		PEGS4.medium = "";
		PEGS4.MTYP = "";
		PEGS4.densPath = "";
		PEGS4.crossPath = "";
		PEGS4.formPath = "";
		PEGS4.GASP = 0.0;
		PEGS4.NE = 0;
		PEGS4.RHO = 0.0;
		PEGS4.WM = 0.0;
		PEGS4.ZC = 0.0;
		PEGS4.TPZ = 0.0;
		PEGS4.EZ = 0.0;
		PEGS4.ZT = 0.0;
		PEGS4.ZA = 0.0;
		PEGS4.ZB = 0.0;
		PEGS4.ZF = 0.0;
		PEGS4.ZG = 0.0;
		PEGS4.ZP = 0.0;
		PEGS4.ZV = 0.0;
		PEGS4.ZU = 0.0;
		PEGS4.ZAB = 0.0;
		PEGS4.ZS = 0.0;
		PEGS4.ZE = 0.0;
		PEGS4.ZX = 0.0;
		PEGS4.TEFF0 = 0.0;
		PEGS4.XCC = 0.0;
		PEGS4.XR0 = 0.0;
		PEGS4.BLCC = 0.0;
		PEGS4.EDEN = 0.0;
		PEGS4.RLC = 0.0;
		PEGS4.IMIXT = 0;
		PEGS4.ZTBL = 0.0;
		PEGS4.IRAYL = 1;
		PEGS4.IUNRST = 0;
		PEGS4.EPSTFL = 1;
		PEGS4.IEPST = 1;
		PEGS4.IAPRIM = 1;
		PEGS4.IAPRFL = 0;
		PEGS4.AFACT = 0.0;
		PEGS4.SK = 0.0;
		PEGS4.X0 = 0.0;
		PEGS4.X1 = 0.0;
		PEGS4.CBAR = 0.0;
		PEGS4.IEV = 0.0;
		PEGS4.ISSB = 0;
		PEGS4.LMED = 1;
		PEGS4.SPC1 = 0.0;
		PEGS4.SPC2 = 0.0;
		PEGS4.NEPST = 0;
		PEGS4.NELEPS = 0;
		PEGS4.AP = 0.;
		PEGS4.AE = 0.;
		PEGS4.UP = 0.;
		PEGS4.UE = 0.;
		PEGS4.AL2 = 0.;
		PEGS4.DELCM = 0.;
		PEGS4.TE = 0.0;
		PEGS4.TET2 = 0.0;
		PEGS4.TEM = 0.0;
		PEGS4.THBREM = 0.0;
		PEGS4.THMOLL = 0.0;
		PEGS4.EBINDA = 0.0;
		PEGS4.NEL = 0;
		PEGS4.AXE = 0.;
		PEGS4.BXE = 0.;
		PEGS4.EPE = 0.01;
		PEGS4.NIPE = 20;
		PEGS4.DELC = 0.;
		PEGS4.CONST = 0.;
		PEGS4.XLNZ = 0.;
		PEGS4.DELTAM = 0.;
		PEGS4.NAPRZ = 0;
		PEGS4.NAPRE = 0;
		PEGS4.E = 0.;
		PEGS4.K = 0.;
		PEGS4.funcname = "";
		// --------------------------
	}

	// set AE,AP,UE,UP
	/**
	 * Set default thresholds, AE,AP,UE and UP (in MeV)
	 */
	public static void setELimit() {
		// PEGS4.defaultAEAPUEUP();
		PEGS4.setAE(0.513);
		PEGS4.setUE(20.);
		PEGS4.setAP(0.010);
		PEGS4.setUP(100.);// MEV!
		cuttB = false;
	}

	/**
	 * Set thresholds in MeV
	 * @param AE AE
	 * @param AP AP
	 * @param UE UE
	 * @param UP UP
	 */
	public static void setELimit(double AE, double AP, double UE, double UP) {
		PEGS4.setAE(AE);
		PEGS4.setUE(UE);
		PEGS4.setAP(AP);
		PEGS4.setUP(UP);// MEV!
		cuttB = true;
	}

	/**
	 * Build data required to plot electron-positron cross sections in range [low energy, high energy].
	 * @param low low energy in MeV
	 * @param high high energy in MeV
	 */
	public static void plotEP(double low, double high) {
		if (low >= high)
			return;
		double w = (high - low) / interve;
		energe = new double[interve];
		brem = new double[interve];
		moll = new double[interve];
		bhab = new double[interve];
		anni = new double[interve];
		tote = new double[interve];
		totp = new double[interve];
		for (int i = 0; i < interve; i++) {
			energe[i] = low + i * w;// do not care if exceed UE
			double[] da = getEPCS(energe[i]);
			brem[i] = da[0];
			moll[i] = da[1];
			bhab[i] = da[2];
			anni[i] = da[3];
			tote[i] = da[4];
			totp[i] = da[5];
		}
	}

	/**
	 * Build data required to plot photon cross sections in range [low energy, high energy].
	 * @param low low energy in MeV
	 * @param high high energy in MeV
	 */
	public static void plotGamma(double low, double high) {
		if (low >= high)
			return;
		double w = (high - low) / interv;
		energ = new double[interv];
		photo = new double[interv];
		compt = new double[interv];
		pair = new double[interv];
		cohe = new double[interv];
		totg = new double[interv];
		for (int i = 0; i < interv; i++) {
			energ[i] = low + i * w;// do not care if exceed UP
			double[] da = getGammaCS(energ[i]);
			photo[i] = da[0];
			compt[i] = da[1];
			pair[i] = da[2];
			cohe[i] = da[3];
			totg[i] = da[4];
		}
	}

	/**
	 * Called by plotGamma. Given energy E, it computes cross sections for specific processes: Photo-electric effect, 
	 * Compton scattering, Pair production, Coherent (Rayleigh) scattering and return the array of cross sections in cm^2/g. 
	 * @param E E in MeV
	 * @return the result
	 */
	public static double[] getGammaCS(double E) {
		double[] result = new double[5];
		// Total Photoelectric cross section
		double d = PEGS4.PHOTTE(E);// in rad length unit
		d = d / PEGS4.RLC;// in cm**-1
		d = d / (PEGS4.RHO);// in cm**2/g
		result[0] = d;
		// Total Compton (incoh) cross section
		d = PEGS4.COMPTM(E);// in rad length unit
		d = d / PEGS4.RLC;// in cm**-1
		d = d / (PEGS4.RHO);// in cm**2/g
		result[1] = d;
		// Total Pair production cross section
		d = PEGS4.PAIRTU(E);// in rad length unit
		d = d / PEGS4.RLC;// in cm**-1
		d = d / (PEGS4.RHO);// in cm**2/g
		result[2] = d;
		// Total Rayleigh (coh) cross section
		d = PEGS4.COHETM(E);// in rad length unit
		d = d / PEGS4.RLC;// in cm**-1
		d = d / (PEGS4.RHO);// in cm**2/g
		result[3] = d;
		// total
		result[4] = result[0] + result[1] + result[2] + result[3];
		// System.out.println("ge: "+E+" ph "+result[0]+" c  "+result[1]+" p "+result[2]+" coh "+result[3]);
		return result;

	}

	/**
	 * Called by plotEP. Given energy E, it computes cross sections for specific processes: Bremsstrahlung production, 
	 * Moller interaction, Bhabha interaction, Annihilation process and return the array of cross sections in cm^2/g. 
	 * @param E E in MeV
	 * @return the result
	 */
	public static double[] getEPCS(double E)// electron positron
	{
		double[] result = new double[6];
		// Total BREMS CS
		double d = PEGS4.BREMTM(E);// in rad length unit
		d = d / PEGS4.RLC;// in cm**-1
		d = d / (PEGS4.RHO);// in cm**2/g
		result[0] = d;
		// Total Moller cross section
		d = PEGS4.AMOLTM(E);// System.out.println("rl: "+d);//in rad length unit
		d = d / PEGS4.RLC;// System.out.println("cm: "+d);//in cm**-1
		d = d / (PEGS4.RHO);// System.out.println("g: "+d);//in cm**2/g
		result[1] = d;
		// Total Bhabha cross section
		d = PEGS4.BHABTM(E);// in rad length unit
		d = d / PEGS4.RLC;// in cm**-1
		d = d / (PEGS4.RHO);// in cm**2/g
		result[2] = d;
		// Total annihilation cross section
		d = PEGS4.ANIHTM(E);// in rad length unit
		d = d / PEGS4.RLC;// in cm**-1
		d = d / (PEGS4.RHO);// in cm**2/g
		result[3] = d;
		// totre
		result[4] = result[0] + result[1];
		// totp
		result[5] = result[0] + result[2] + result[3];
		// System.out.println("en: "+E+" b "+result[0]+" m  "+result[1]+" bh "+result[2]+" a "+result[3]);
		return result;
	}

	// retrieve total cross section for electrons,positrons and electrons at a
	// specific energy
	// it also give the appropriate branching ratios
	/**
	 * Retrieve total cross sections for electrons, positrons and photons at a specific energy. 
	 * It also give the appropriate branching ratios
	 * @param E E
	 * @return the result
	 */
	public static double[] getTotalCrossSections(double E) {
		double[] result = new double[57];
		// electron and positron
		// Total Brems cross section
		double d = PEGS4.BREMTM(E);// in rad length unit
		result[0] = d;
		d = d / PEGS4.RLC;// in cm**-1
		result[1] = d;
		d = d / (PEGS4.RHO);// in cm**2/g
		result[2] = d;
		// 1.6605402E-24 is u.a.m. in grams
		d = d * PEGS4.WM * g24 * PEGS4.NE / (PEGS4.TPZ);// in barns/atoms
		result[3] = d;
		// Total Moller cross section
		d = PEGS4.AMOLTM(E);// in rad length unit
		result[4] = d;
		d = d / PEGS4.RLC;// in cm**-1
		result[5] = d;
		d = d / (PEGS4.RHO);// in cm**2/g
		result[6] = d;
		d = d * PEGS4.WM * g24 * PEGS4.NE / (PEGS4.TPZ);// in barns/atoms
		result[7] = d;
		// "***BRANCHING RATIO BREMS/(BREMS+MOLLER) AND TOTAL ELECTRON CROSS
		// SECTION
		double TEBR = result[0] + result[4];// "TOTAL ELECTRON CROSS-SECTION" in
											// radlength
		result[8] = TEBR;
		d = result[1] + result[5];// in cm**-1
		result[9] = d;
		d = result[2] + result[6];// in cm**2/g
		result[10] = d;
		d = result[3] + result[7];// in barns/atoms
		result[11] = d;
		if (TEBR > 0.0) {
			d = result[0] / TEBR;
		} else {
			// "BELOW THRESHOLD FOR BOTH BREMS AND MOLLER. USE THE BRANCHING"
			// "RATIO THAT EXISTED WHEN CROSS SECTION APPROACHED ZERO"
			if (PEGS4.THBREM <= PEGS4.THMOLL) {
				d = 1.0;
			} else {
				d = 0.0;
			}
		}
		result[12] = d;
		// Total Bhabha cross section
		d = PEGS4.BHABTM(E);// in rad length unit
		result[13] = d;
		d = d / PEGS4.RLC;// in cm**-1
		result[14] = d;
		d = d / (PEGS4.RHO);// in cm**2/g
		result[15] = d;
		d = d * PEGS4.WM * g24 * PEGS4.NE / (PEGS4.TPZ);// in barns/atoms
		result[16] = d;
		// Total annihilation cross section
		d = PEGS4.ANIHTM(E);// in rad length unit
		result[17] = d;
		d = d / PEGS4.RLC;// in cm**-1
		result[18] = d;
		d = d / (PEGS4.RHO);// in cm**2/g
		result[19] = d;
		d = d * PEGS4.WM * g24 * PEGS4.NE / (PEGS4.TPZ);// in barns/atoms
		result[20] = d;
		// "***BRANCHING RATIOS AND TOTAL POSITRON CROSS SECTION
		double PSIG = result[0] + result[13] + result[17];// "TOTAL POSITRON CROSS SECTION"
															// in radlength
		result[21] = PSIG;
		d = result[1] + result[14] + result[18];// //in cm**-1
		result[22] = d;
		d = result[2] + result[15] + result[19];// //in cm**2/g
		result[23] = d;
		d = result[3] + result[16] + result[20];// //in barns/atoms
		result[24] = d;
		// "PBR1=BREM/(BREM+BHABA+ANNIH)"
		result[25] = result[0] / PSIG;
		// "PBR2=(BREM+BHABA)/(PSIG)"
		result[26] = (result[0] + result[13]) / PSIG;
		// //PHOTON
		// Total Photoelectric cross section
		d = PEGS4.PHOTTE(E);// in rad length unit
		result[27] = d;
		d = d / PEGS4.RLC;// in cm**-1
		result[28] = d;
		d = d / (PEGS4.RHO);// in cm**2/g
		result[29] = d;
		d = d * PEGS4.WM * g24 * PEGS4.NE / (PEGS4.TPZ);// in barns/atoms
		result[30] = d;
		// Total Compton (incoh) cross section
		d = PEGS4.COMPTM(E);// in rad length unit
		result[31] = d;
		d = d / PEGS4.RLC;// in cm**-1
		result[32] = d;
		d = d / (PEGS4.RHO);// in cm**2/g
		result[33] = d;
		d = d * PEGS4.WM * g24 * PEGS4.NE / (PEGS4.TPZ);// in barns/atoms
		result[34] = d;
		// Total Pair production cross section
		d = PEGS4.PAIRTU(E);// in rad length unit
		result[35] = d;
		d = d / PEGS4.RLC;// in cm**-1
		result[36] = d;
		d = d / (PEGS4.RHO);// in cm**2/g
		result[37] = d;
		d = d * PEGS4.WM * g24 * PEGS4.NE / (PEGS4.TPZ);// in barns/atoms
		result[38] = d;
		// Total Rayleigh (coh) cross section
		d = PEGS4.COHETM(E);// in rad length unit
		result[39] = d;
		d = d / PEGS4.RLC;// in cm**-1
		result[40] = d;
		d = d / (PEGS4.RHO);// in cm**2/g
		result[41] = d;
		d = d * PEGS4.WM * g24 * PEGS4.NE / (PEGS4.TPZ);// in barns/atoms
		result[42] = d;
		// Total without coh
		double TSANSC = result[27] + result[31] + result[35];
		result[43] = TSANSC;// radlength
		result[44] = result[28] + result[32] + result[36];// in cm**-1
		result[45] = result[29] + result[33] + result[37];// in cm**2/g
		result[46] = result[30] + result[34] + result[38];// in barns/atoms
		// Total + coh
		double TC = result[27] + result[31] + result[35] + result[39];
		result[47] = TC;
		result[48] = result[28] + result[32] + result[36] + result[40];// in
																		// cm**-1
		result[49] = result[29] + result[33] + result[37] + result[41];// in
																		// cm**2/g
		result[50] = result[30] + result[34] + result[38] + result[42];// in
																		// barns/atoms
		// gbr1=PAIR/TSANSC
		result[51] = result[35] / TSANSC;
		result[52] = result[35] / TC;
		// gbr2=PAIR+COMP/TSANSC
		result[53] = (result[35] + result[31]) / TSANSC;
		result[54] = (result[35] + result[31]) / TC;
		// CRATIO=TSANSC/(TSANSC+COHR);
		result[55] = (TSANSC) / (TSANSC + result[39]);
		result[56] = (result[39]) / TC;

		return result;

	}

	// retrieve electron,positron,photon and rayleigh vfuns as are defined in
	// pegs4!!
	/**
	 * Retrieve electron, positron, photon and rayleigh functions vfuns as are defined in PEGS4.
	 * @param E E
	 * @return the result
	 */
	public static double[] getVfuns(double E) {
		double[] result = new double[13];
		double[] epos = new double[8];
		PEGS4.EFUNS(E, epos);
		double[] gamma = new double[4];
		PEGS4.GFUNS(E, gamma);
		double[] ray = new double[1];
		PEGS4.RFUNS(E, ray);
		// depends on PEGS4.IUNRST
		result[0] = epos[0];// esig-total cs without elastic for electrons
		result[1] = epos[1];// psig-total cs without elastic for positrons
		result[2] = epos[2];// // "TOTAL ELECTRON STOPPING POWER"
		result[3] = epos[3];// // "TOTAL POSITRON STOPPING POWER"
		result[4] = epos[4];// BREM/ESIG=brems branching ratio for electrons
		result[5] = epos[5];// BREM/PSIG=brems branching ratio for positrons
		result[6] = epos[6];// (BREM+BHABA)/(PSIG) branching ratio for positrons
		result[7] = epos[7];// MAXIMUM ALLOWED TRANSPORT STEP, FROM MULTIPLE
							// SCATTERING
		result[8] = gamma[0];// GAMMA MEAN FREE PATH=1/(PAIR+COMP+PHOT)
		result[9] = gamma[1];// //PAIR*GMFP branching ratio
		result[10] = gamma[2];// //(PAIR+COMP)*GMFP branching ratio
		result[11] = gamma[3];// //COHERENT RATIO
		result[12] = ray[0];

		return result;

	}

	// public static void preinit(String med,int imix,int nelem,double
	// rhoe,String[] asy,
	/**
	 * Pre-initialization method to be called first for "manual" medium selection.
	 * @param med medium name
	 * @param nelem nelem number of elements in this medium (mixture)
	 * @param rhoe rhoe medium density in g/cm^3
	 * @param asy is the name of each individual element (match the ASSYM table, see PEGS4)
	 * @param rhozz rhozz is the partial density of each individual element (fraction by weight)
	 * @param pzee this is just 0 here, we use rhozz instead. 
	 */
	public static void preinit(String med, int nelem, double rhoe,
			String[] asy, double[] rhozz, double[] pzee) {
		medium = med;
		// imixt=imix;
		ne = nelem;
		rho = rhoe;
		asymetry = new String[asy.length];
		rhoze = new double[rhozz.length];
		pze = new double[pzee.length];
		for (int i = 0; i < ne; i++) {
			asymetry[i] = asy[i];// System.out.println(asymetry[i]);
			rhoze[i] = rhozz[i];
			pze[i] = pzee[i];
		}
	}

	/**
	 * Initialize PEGS4 for material creation.
	 */
	public static void init_pegs4() {
		// "COMPUTE PHYSICAL AND MATHEMATICAL CONSTANTS"
		PEGS4.PMDCON();// constant initialization
		// setting gasp state
		PEGS4.setGASP(gasPressure);
		// "READ IN PHOTOELECTRIC AND PAIR DATA FROM FILE"
		PEGS4.readPhotoelectricPairData(phpData);
		// "READ IN ATOMIC FORM FACTOR DATA FROM A FILE"
		PEGS4.readFormFactorData(formData);
		// NOT ALLWAYS READ IN DENSITY FILE
		String nuls = "";
		if ((densData.compareTo(nuls) == 0) && (!densityB)) {// invalid density
																// file and not
																// file read
			if (ne == 0)// no preinitialization
				return;
			PEGS4.setMedium(medium);// 1=mixture,0=elem or compound
			// PEGS4.setIMIXT(imixt);//1=mixture,0=elem or compound
			PEGS4.setNE(ne);// number of elements
			PEGS4.setRHO(rho);// nai in g/cm3
			PEGS4.setNEVAR2(asymetry, pze, rhoze);
		} else {
			if (!densityFullB)
				PEGS4.readDensityFile(densData, densB);
			else
				// densityFullB=true;
				PEGS4.readDensityFile(densData);
		}
		// set energy cutoffs
		if (!cuttB)
			setELimit();// default
		else {
			if (PEGS4.AE < PEGS4.RM)// could not be lower than electron rest
									// mass MEV!!
				return;
			if (PEGS4.AE >= PEGS4.UE)// could not be higher than upper limit!!
				return;
			if (PEGS4.AP >= PEGS4.UP)// could not be higher than upper limit!!
				return;
		}
		// "NOW WE HAVE DEFINED THE MIXTURE, SO COMPUTE Z-DEPENDENT PARAMETER"
		PEGS4.MIX();
		PEGS4.SPINIT();// default read density correction from file
		// "     COMPUTE DIFFERENTIAL SAMPLING(BREMPR) CONSTANTS.                 "
		PEGS4.DIFFER();
		// SET ENERGY LIMITS FOR ELECTRONS AND PHOTONS"
		PEGS4.setENER();
		// COMPUTE MATERIAL INDEPENDENT MULTIPLE SCATTERING DATA.
		PEGS4.MOLIER();
		// FIRST FIT ELECTRON(AND POSITRON)"
		PEGS4.FITEP();
		func = "ALKE";
		funci = "ALKEI";
		funceval = "EFUNS";
		String nis = "NEL";
		String axs = "AXE";
		String bxs = "BXE";
		String afes = "AFE";
		String bfes = "BFE";
		PWLF1(nis, NALE, PEGS4.AE, PEGS4.UE, PEGS4.THMOLL, EPE, ZTHRE, ZEPE,
				NIPE, func, funci, axs, bxs, $MXEKE, 8, afes, bfes, funceval);
		// " NEXT FIT GAMMA
		func = "ALOG";
		funci = "EXP";
		funceval = "GFUNS";
		nis = "NGL";
		axs = "AXG";
		bxs = "BXG";
		afes = "AFG";
		bfes = "BFG";
		PWLF1(nis, NALG, PEGS4.AP, PEGS4.UP, PEGS4.RMT2, EPG, ZTHRG, ZEPG,
				NIPG, func, funci, axs, bxs, $MXGE, 4, afes, bfes, funceval);
		if (PEGS4.IRAYL == 1) {
			PEGS4.Rayleigh();
			func = "ALIN";
			funci = "ALINI";
			funceval = "RFUNS";
			nis = "NGR";
			axs = "AXR";
			bxs = "BXR";
			afes = "AFR";
			bfes = "BFR";
			PWLF1(nis, NALR, 0.0, 1.0, 0.0, EPR, ZTHRR, ZEPR, NIPR, func,
					funci, axs, bxs, $MXRL, 1, afes, bfes, funceval);
		}
		// HERE FOR OPTION TO PRODUCE A DECK---PASS THE BUCK TO A SUBROUTINE"
		LAY();
		/*
		 * //-------------MIX OUTPUT-------------------
		 * System.out.println("-----NE----------"+PEGS4.NE);
		 * System.out.println("-----RHO----------"+PEGS4.RHO);
		 * System.out.println("-----IEV----------"+PEGS4.IEV); for (int i=0;
		 * i<PEGS4.NE; i++) {
		 * System.out.println("Z: "+PEGS4.Z[i]+"   RHOZ: "+PEGS4
		 * .RHOZ[i]+"  WA:  "+PEGS4.WA[i] +"  PZ:  "+PEGS4.PZ[i]); }
		 * System.out.println("----------------------------------------------");
		 * double wm=PEGS4.WM;System.out.println("WM: "+wm); double
		 * zc=PEGS4.ZC;System.out.println("ZC: "+zc); double
		 * zt=PEGS4.ZT;System.out.println("ZT: "+zt); double
		 * za=PEGS4.ZA;System.out.println("ZA: "+za); double
		 * zb=PEGS4.ZB;System.out.println("ZB: "+zb); double
		 * zab=PEGS4.ZAB;System.out.println("ZAB: "+zab);
		 * System.out.println("----------------------------------------------");
		 * double zf=PEGS4.ZF;System.out.println("ZF: "+zf); double
		 * zg=PEGS4.ZG;System.out.println("ZG: "+zg); double
		 * zp=PEGS4.ZP;System.out.println("ZP: "+zp); double
		 * zv=PEGS4.ZV;System.out.println("ZV: "+zv); double
		 * zu=PEGS4.ZU;System.out.println("ZU: "+zu); double
		 * zs=PEGS4.ZS;System.out.println("ZS: "+zs);
		 * System.out.println("----------------------------------------------");
		 * double zee=PEGS4.ZE;System.out.println("ZE: "+zee); double
		 * zx=PEGS4.ZX;System.out.println("ZX: "+zx); double
		 * rlc=PEGS4.RLC;System.out.println("RLC: "+rlc); double
		 * eden=PEGS4.EDEN;System.out.println("EDEN: "+eden);
		 * System.out.println("----------------------------------------------");
		 * for (int i=0; i<PEGS4.NE; i++) { int j=i+1;
		 * System.out.println("Elem: "
		 * +j+" XSI: "+PEGS4.XSI[i]+" ZZX: "+PEGS4.ZZX[i]+" FZC: "+PEGS4.FZC[i]+
		 * " FCOUL: "+PEGS4.FCOUL[i]+" ZZ: "+PEGS4.ZZ[i]);
		 * System.out.println(" PZ@: "
		 * +PEGS4.PZ[i]+" Z@: "+PEGS4.Z[i]+" WA@: "+PEGS4.WA[i]); }
		 * System.out.println("----------------------------------------------");
		 * double BLCC=PEGS4.BLCC;System.out.println("BLCC: "+BLCC); double
		 * XCC=PEGS4.XCC;System.out.println("XCC: "+XCC); double
		 * TEFF0=PEGS4.TEFF0;System.out.println("TEFF0: "+TEFF0); double
		 * XR0=PEGS4.XR0;System.out.println("XR0: "+XR0); System.out.println(
		 * "-----------------------------------------------------------------");
		 * int ll=PEGS4.DL1.length; for (int i=0;i<ll;i++) {
		 * System.out.println("dl1["+i+"]= "+PEGS4.DL1[i]);
		 * System.out.println("dl2["+i+"]= "+PEGS4.DL2[i]);
		 * System.out.println("dl3["+i+"]= "+PEGS4.DL3[i]);
		 * System.out.println("dl4["+i+"]= "+PEGS4.DL4[i]);
		 * System.out.println("dl5["+i+"]= "+PEGS4.DL5[i]);
		 * System.out.println("dl6["+i+"]= "+PEGS4.DL6[i]); }
		 * System.out.println(
		 * "ALPHI[0]= "+PEGS4.ALPHI[0]+"  ALPHI[1]= "+PEGS4.ALPHI[1]);
		 * System.out
		 * .println("BPAR[0]= "+PEGS4.BPAR[0]+"  BPAR[1]= "+PEGS4.BPAR[1]);
		 * System.out.println("DELCM= "+PEGS4.DELCM);
		 * System.out.println("DELPOS[0]= "
		 * +PEGS4.DELPOS[0]+"  DELPOS[1]= "+PEGS4.DELPOS[1]);
		 * System.out.println(
		 * "-----------------------------------------------------------------");
		 * System
		 * .out.println("AE= "+PEGS4.AE+"  UE= "+PEGS4.UE+"  AP= "+PEGS4.AP
		 * +"  UP= "+PEGS4.UP);
		 * System.out.println("TE= "+PEGS4.TE+"  TET2= "+PEGS4
		 * .TET2+"  TEM= "+PEGS4
		 * .TEM+"  THBREM= "+PEGS4.THBREM+"  THMOLL= "+PEGS4.THMOLL);
		 * System.out.println(
		 * "-----------------------------------------------------------------");
		 * System.out.println(" AVERAGE K-IONIZATION ENERGY(MEV) "+
		 * PEGS4.EBINDA);//WORKS
		 */
	}

	// no lin fit just exact calculations of electrons,positrons and gamma
	// functions
	// is less rapid then egs4 initial evaluations but more exact
	/**
	 * Initialize PEGS4 for computing cross section and plotting charts. It performs exact calculations of 
	 * electrons, positrons and gamma functions.
	 */
	public static void init_pegs4_exact() {
		// "COMPUTE PHYSICAL AND MATHEMATICAL CONSTANTS"
		PEGS4.PMDCON();// constant initialization
		// setting gasp state
		PEGS4.setGASP(gasPressure);
		// "READ IN PHOTOELECTRIC AND PAIR DATA FROM FILE"
		PEGS4.readPhotoelectricPairData(phpData);
		// "READ IN ATOMIC FORM FACTOR DATA FROM A FILE"
		PEGS4.readFormFactorData(formData);
		// NOT ALLWAYS READ IN DENSITY FILE
		String nuls = "";
		if ((densData.compareTo(nuls) == 0) && (!densityB)) {// invalid density
																// file and not
																// file read
			if (ne == 0)// no preinitialization
				return;
			PEGS4.setMedium(medium);// 1=mixture,0=elem or compound
			// PEGS4.setIMIXT(imixt);//1=mixture,0=elem or compound
			PEGS4.setNE(ne);// number of elements
			PEGS4.setRHO(rho);// nai in g/cm3
			PEGS4.setNEVAR2(asymetry, pze, rhoze);
			// System.out.println("enter");
		} else {
			if (!densityFullB) {
				PEGS4.readDensityFile(densData, densB);// System.out.println("enter");
			} else {
				PEGS4.readDensityFile(densData);// System.out.println(densData);
			}
		}
		// set energy cutoffs
		if (!cuttB)
			setELimit();// default
		else {
			if (PEGS4.AE < PEGS4.RM)// could not be lower than electron rest
									// mass MEV!!
				return;
			if (PEGS4.AE >= PEGS4.UE)// could not be higher than upper limit!!
				return;
			if (PEGS4.AP >= PEGS4.UP)// could not be higher than upper limit!!
				return;
		}
		// "NOW WE HAVE DEFINED THE MIXTURE, SO COMPUTE Z-DEPENDENT PARAMETER"
		PEGS4.MIX();
		PEGS4.SPINIT();// default read density correction from file
		// "     COMPUTE DIFFERENTIAL SAMPLING(BREMPR) CONSTANTS.                 "
		PEGS4.DIFFER();
		// SET ENERGY LIMITS FOR ELECTRONS AND PHOTONS"
		PEGS4.setENER();
		// COMPUTE MATERIAL INDEPENDENT MULTIPLE SCATTERING DATA.
		PEGS4.MOLIER();
		// FIRST FIT ELECTRON(AND POSITRON)"
		PEGS4.FITEP();
		// from egs4 just call getVfuns(double E)
		/*
		 * //-------------MIX OUTPUT-------------------
		 * System.out.println("-----NE----------"+PEGS4.NE);
		 * System.out.println("-----RHO----------"+PEGS4.RHO);
		 * System.out.println("-----IEV----------"+PEGS4.IEV); for (int i=0;
		 * i<PEGS4.NE; i++) {
		 * System.out.println("Z: "+PEGS4.Z[i]+"   RHOZ: "+PEGS4
		 * .RHOZ[i]+"  WA:  "+PEGS4.WA[i] +"  PZ:  "+PEGS4.PZ[i]); }
		 * System.out.println("----------------------------------------------");
		 * double wm=PEGS4.WM;System.out.println("WM: "+wm); double
		 * zc=PEGS4.ZC;System.out.println("ZC: "+zc); double
		 * zt=PEGS4.ZT;System.out.println("ZT: "+zt); double
		 * za=PEGS4.ZA;System.out.println("ZA: "+za); double
		 * zb=PEGS4.ZB;System.out.println("ZB: "+zb); double
		 * zab=PEGS4.ZAB;System.out.println("ZAB: "+zab);
		 * System.out.println("----------------------------------------------");
		 * double zf=PEGS4.ZF;System.out.println("ZF: "+zf); double
		 * zg=PEGS4.ZG;System.out.println("ZG: "+zg); double
		 * zp=PEGS4.ZP;System.out.println("ZP: "+zp); double
		 * zv=PEGS4.ZV;System.out.println("ZV: "+zv); double
		 * zu=PEGS4.ZU;System.out.println("ZU: "+zu); double
		 * zs=PEGS4.ZS;System.out.println("ZS: "+zs);
		 * System.out.println("----------------------------------------------");
		 * double zee=PEGS4.ZE;System.out.println("ZE: "+zee); double
		 * zx=PEGS4.ZX;System.out.println("ZX: "+zx); double
		 * rlc=PEGS4.RLC;System.out.println("RLC: "+rlc); double
		 * eden=PEGS4.EDEN;System.out.println("EDEN: "+eden);
		 * System.out.println("----------------------------------------------");
		 * for (int i=0; i<PEGS4.NE; i++) { int j=i+1;
		 * System.out.println("Elem: "
		 * +j+" XSI: "+PEGS4.XSI[i]+" ZZX: "+PEGS4.ZZX[i]+" FZC: "+PEGS4.FZC[i]+
		 * " FCOUL: "+PEGS4.FCOUL[i]+" ZZ: "+PEGS4.ZZ[i]);
		 * System.out.println(" PZ@: "
		 * +PEGS4.PZ[i]+" Z@: "+PEGS4.Z[i]+" WA@: "+PEGS4.WA[i]); }
		 * System.out.println("----------------------------------------------");
		 * double BLCC=PEGS4.BLCC;System.out.println("BLCC: "+BLCC); double
		 * XCC=PEGS4.XCC;System.out.println("XCC: "+XCC); double
		 * TEFF0=PEGS4.TEFF0;System.out.println("TEFF0: "+TEFF0); double
		 * XR0=PEGS4.XR0;System.out.println("XR0: "+XR0); System.out.println(
		 * "-----------------------------------------------------------------");
		 * int ll=PEGS4.DL1.length; for (int i=0;i<ll;i++) {
		 * System.out.println("dl1["+i+"]= "+PEGS4.DL1[i]);
		 * System.out.println("dl2["+i+"]= "+PEGS4.DL2[i]);
		 * System.out.println("dl3["+i+"]= "+PEGS4.DL3[i]);
		 * System.out.println("dl4["+i+"]= "+PEGS4.DL4[i]);
		 * System.out.println("dl5["+i+"]= "+PEGS4.DL5[i]);
		 * System.out.println("dl6["+i+"]= "+PEGS4.DL6[i]); }
		 * System.out.println(
		 * "ALPHI[0]= "+PEGS4.ALPHI[0]+"  ALPHI[1]= "+PEGS4.ALPHI[1]);
		 * System.out
		 * .println("BPAR[0]= "+PEGS4.BPAR[0]+"  BPAR[1]= "+PEGS4.BPAR[1]);
		 * System.out.println("DELCM= "+PEGS4.DELCM);
		 * System.out.println("DELPOS[0]= "
		 * +PEGS4.DELPOS[0]+"  DELPOS[1]= "+PEGS4.DELPOS[1]);
		 * System.out.println(
		 * "-----------------------------------------------------------------");
		 * System
		 * .out.println("AE= "+PEGS4.AE+"  UE= "+PEGS4.UE+"  AP= "+PEGS4.AP
		 * +"  UP= "+PEGS4.UP);
		 * System.out.println("TE= "+PEGS4.TE+"  TET2= "+PEGS4
		 * .TET2+"  TEM= "+PEGS4
		 * .TEM+"  THBREM= "+PEGS4.THBREM+"  THMOLL= "+PEGS4.THMOLL);
		 * System.out.println(
		 * "-----------------------------------------------------------------");
		 * System.out.println(" AVERAGE K-IONIZATION ENERGY(MEV) "+
		 * PEGS4.EBINDA);//WORKS
		 */
	}

	// "PRODUCES DECK OF MATERIAL DEPENDENT DATA.
	/**
	 * Called by init_pegs4. PRODUCES DECK OF MATERIAL DEPENDENT DATA.
	 */
	public static void LAY() {
		int NSGE = 0;
		int NSEKE = 0;
		int NLEKE = 0;
		int NCMFP = 0;
		int NRANGE = 0;
		//String tab = "   ";
		String sp = " ";
		String filename = "";
		if (mediumB)
			filename = PEGS4.datas +PEGS4.file_sep+PEGS4.egsData +PEGS4.file_sep + intDir + PEGS4.file_sep
					+ PEGS4.medium + extension;
		else
			filename = PEGS4.datas +PEGS4.file_sep+PEGS4.egsData+ PEGS4.file_sep + intDir + PEGS4.file_sep
					+ mediumFilename + extension;
		// "   PUT OUT HEADING, AND COMPOSITION CARDS"
		String s = " MEDIUM=" + PEGS4.medium + ",STERNCID=" + PEGS4.medium
				+ "\n";
		if (PEGS4.GASP != 0.0) {// "GASES"
			s = s + " " + PEGS4.MTYP + ",RHO=" + PEGS4.RHO + ",NE=" + PEGS4.NE
					+ ",GASP=" + PEGS4.GASP + ", IUNRST=" + PEGS4.IUNRST
					+ ", EPSTFL=" + PEGS4.EPSTFL + ", IAPRIM=" + PEGS4.IAPRIM
					+ "\n";
		}
		// "   NRC MOD  ========================================="
		else {
			s = s + " " + PEGS4.MTYP + ",RHO=" + PEGS4.RHO + ",NE=" + PEGS4.NE
					+ ", IUNRST=" + PEGS4.IUNRST + ", EPSTFL=" + PEGS4.EPSTFL
					+ ", IAPRIM=" + PEGS4.IAPRIM + "\n";
		}
		for (int IE = 0; IE < PEGS4.NE; IE++) {
			s = s + " ASYM=" + PEGS4.ASYM[IE] + ",Z=" + PEGS4.Z[IE] + ",A="
					+ PEGS4.WA[IE] + ",PZ=" + PEGS4.PZ[IE] + ",RHOZ="
					+ PEGS4.RHOZ[IE] + "\n";
		}
		// "NOW COMES THE DATA PROPER"
		// RLC,AE,AP,UE,UP;
		// s=s+" "+tab+PEGS4.RLC+tab+PEGS4.AE+tab+PEGS4.AP+tab+PEGS4.UE+tab+PEGS4.UP+"\n";
		s = s + sp + PEGS4.RLC + sp + PEGS4.AE + sp + PEGS4.AP + sp + PEGS4.UE
				+ sp + PEGS4.UP + "\n";
		// "FAKE SOME PARAMETERS FOR NOW"
		int NGE = NGL;
		int NEKE = NEL; // "CHANGE NAMES OF SOME VARIABLES"
		s = s + sp + NSGE + sp + NGE + sp + NSEKE + sp + NEKE + sp + NLEKE + sp
				+ NCMFP + sp + NRANGE + sp + PEGS4.IRAYL + sp + PEGS4.IUNRST
				+ "\n";
		// $ECHOWRITE(IP,:FLT:)(DL1(I),DL2(I),DL3(I),DL4(I),DL5(I),DL6(I),I=1,6);
		for (int I = 0; I < 6; I++) {
			s = s + sp + PEGS4.DL1[I] + sp + PEGS4.DL2[I] + sp + PEGS4.DL3[I]
					+ sp + PEGS4.DL4[I] + sp + PEGS4.DL5[I] + sp + PEGS4.DL6[I]
					+ "\n";
		}
		// $ECHOWRITE(IP,:FLT:) DELCM,(ALPHI(I),BPAR(I),DELPOS(I),I=1,2);
		s = s + sp + PEGS4.DELCM + "\n";
		for (int I = 0; I < 2; I++) {
			s = s + sp + PEGS4.ALPHI[I] + sp + PEGS4.BPAR[I] + sp
					+ PEGS4.DELPOS[I] + "\n";
		}
		// $ECHOWRITE(IP,:FLT:) XR0,TEFF0,BLCC,XCC;
		s = s + sp + PEGS4.XR0 + sp + PEGS4.TEFF0 + sp + PEGS4.BLCC + sp
				+ PEGS4.XCC + "\n";
		// $ECHOWRITE(IP,:FLT:) BXE,AXE;
		s = s + sp + BXE + sp + AXE + "\n";
		// $ECHOWRITE(IP,:FLT:) ((BFE(I,IFUN),AFE(I,IFUN),IFUN=1,8),I=1,NEKE);
		for (int i = 0; i < NEKE; i++) {
			for (int IFUN = 0; IFUN < 8; IFUN++) {
				s = s + sp + BFE[i][IFUN] + sp + AFE[i][IFUN] + "\n";
			}
		}
		// $ECHOWRITE(IP,:FLT:) EBINDA,BXG,AXG;
		s = s + sp + PEGS4.EBINDA + sp + BXG + sp + AXG + "\n";
		// $ECHOWRITE(IP,:FLT:) ((BFG(I,IFUN),AFG(I,IFUN),IFUN=1,3),I=1,NGE);
		for (int i = 0; i < NGE; i++) {
			for (int IFUN = 0; IFUN < 3; IFUN++) {
				s = s + sp + BFG[i][IFUN] + sp + AFG[i][IFUN] + "\n";
			}
		}
		if (PEGS4.IRAYL != 0) {
			s = s + sp + NGR + "\n";
			s = s + sp + BXR + sp + AXR + "\n";
			for (int I = 0; I < NGR; I++) {
				// $ECHOWRITE(IP,:FLT:) (BFR(I),AFR(I),I=1,NGR);
				s = s + sp + BFR[I][0] + sp + AFR[I][0] + "\n";
			}
			for (int I = 0; I < NGE; I++) {
				// $ECHOWRITE(IP,:FLT:) (BFG(I,4),AFG(I,4),I=1,NGE);
				s = s + sp + BFG[I][3] + sp + AFG[I][3] + "\n";
			}
		}

		// FileOutputStream sigfos = new FileOutputStream(filename);
		try {
			FileWriter sigfos = new FileWriter(filename);
			sigfos.write(s);
			sigfos.close();
		} catch (Exception ex) {

		}
		// --log file
		filename = PEGS4.datas+PEGS4.file_sep+PEGS4.egsData + PEGS4.file_sep + intDir + PEGS4.file_sep
				+ logFile;
		// --------------------------------------------------------------------
		Calendar cal = Calendar.getInstance();
		String data = " dd: "
				+ PEGS4.intToString(cal.get(Calendar.DAY_OF_MONTH)) + " mm: "
				+ PEGS4.intToString(cal.get(Calendar.MONTH)) + " yyyy: "
				+ PEGS4.intToString(cal.get(Calendar.YEAR)) + " ora: "
				+ PEGS4.intToString(cal.get(Calendar.HOUR_OF_DAY)) + ":"
				+ PEGS4.intToString(cal.get(Calendar.MINUTE)) + ":"
				+ PEGS4.intToString(cal.get(Calendar.SECOND));
		s = "---------------LOG FILE---------------------------------------"
				+ "\n";
		s = s + "medium:  " + PEGS4.medium + "\n";
		s = s + "data generated on:  " + data + "\n";
		s = s + "type=" + PEGS4.MTYP;
		s = s + " rho=" + PEGS4.RHO;
		s = s + " ne=" + PEGS4.NE;
		s = s + " iunrst=" + PEGS4.IUNRST;
		s = s + " epstfl=" + PEGS4.EPSTFL;
		s = s + " iaprim=" + PEGS4.IAPRIM + "\n";
		if (PEGS4.IRAYL != 0) {
			s = s + "Rayleigh data included" + "\n";
		} else {
			s = s + "Rayleigh data not included" + "\n";
		}
		if (PEGS4.EPSTFL == 1) {
			s = s + "density correction file: " + PEGS4.densPath + "\n";
		}
		s = s + "cross section data from: " + PEGS4.crossPath + "\n";
		s = s + "form factor data from: " + PEGS4.formPath + "\n";
		s = s + "ae=" + PEGS4.AE + " ap=" + PEGS4.AP + " ue=" + PEGS4.UE
				+ " up=" + PEGS4.UP + "\n";
		s = s + "data written to: " + PEGS4.userdirs + PEGS4.file_sep
				+ filename;

		try {
			FileWriter sigfos = new FileWriter(filename);
			sigfos.write(s);
			sigfos.close();
		} catch (Exception ex) {

		}
		// "THAT'S ALL FOR NOW"
		return;
	}

	/**
	 * Called by init_pegs4. PWLF1 IS A ROUTINE WHICH WILL FIT UP TO 10(CURRENTLY) FUNCTIONS 
	 * SIMULTANEOUSLY ON AN INTERVAL (XL,XU) OF THE INDEPENDENT VARIABLE X OF THE FUNCTIONS. THE FIT IS A PIECEWISE LINEAR FUNCTION OF 
	 * XFUN(X).  XFI IS THE INVERSE FUNCTION OF XFUN.  THE SUBINTERVALS ARE CHOSEN OF UNIFORM WIDTH IN XFUN(X) AND SUFFICIENT OF THEM 
	 * ARE CHOSEN SO THAT THE FIT GIVES A RELATIVE ERROR[EP FOR ALL THE FUNCTIONS OVER ALL THE SUBINTERVALS. 
	 * @param NI ON RETURN IS NUMBER OF SUBINTERVALS USED FOR THE FIT. 
	 * @param NIMX NIMX
	 * @param XL LOWER LIMIT OF INTERVAL ON WHICH TO FIT THE FUNCTIONS.
	 * @param XU UPPER LIMIT 
	 * @param XR VALUE OF X WHICH WILL BE FORCED TO BE A SUBINTERVAL BOUNDARY, IF ONE WANTS AN EXACT FIT OF THE FUNCTIONS 
	 * AT A PARTICULAR POINT, XR SHOULD BE SET TO THAT VALUE.
	 * @param EP THE MAXIMUM RELATIVE ERROR ALLOWED THE FIT. 
	 * @param ZTHR ZTHR
	 * @param ZEP ZEP
	 * @param NIP THE MINUIMUM NUMBER OF POINTS INTERIOR TO (XL,XU) AT WHICH THE FIT IS TO BE TESTED FOR RELATIVE ERROR VS. EP.
	 * @param XFUN A FUNCTION OF X OVER WHICH IT IS HOPED THE FUNCTIONS TO BE FIT ARE MORE LINEAR THAN OVER X.  XFUN IS EXPECTED TO BE MONOTONICALLY INCREASING IN X.
	 * @param XFI THE INVERSE OF XFUN. THAT IS XFI(XFUN(X))=X. 
	 * @param AX COEFFICIENT USED TO DETERMINE WHICH SUBINTERVAL A VALUE OF X IS IN.
	 * @param BX COEFFICIENT USED TO DETERMINE WHICH SUBINTERVAL A VALUE OF X IS IN.
	 * @param NALM THE MAXIMUM NUMBER OF SUBINTERVALS FOR WHICH ARRAY SPACE HAS BEEN ALLOCATED.
	 * @param NFUN NUMBER OF FUNCTIONS TO BE FITTED(SIMULTANEOUSLY,I.E. ALL FUNCTIONS HAVE THE SAME XFUN AND SUBINTERVALS, AND ALL 
	 * ARE REQUIRED TO BE FIT WITH MAX REL ERR[EP]
	 * @param AF ARRAYS OF COEFFICIENTS USED TO GET VALUES OF THE FUNCS.
	 * @param BF ARRAYS OF COEFFICIENTS USED TO GET VALUES OF THE FUNCS.
	 * @param VFUNS SUBROUTINE TO FILL AN ARRAY WITH THE VALUES OF THE FUNCTIONS TO BE FITTED.
	 */
	public static void PWLF1(String NI, int NIMX, double XL, double XU,
			double XR, double EP, double[] ZTHR, double[] ZEP, int NIP,
			String XFUN, String XFI, String AX, String BX, int NALM, int NFUN,
			String AF, String BF, String VFUNS) {
		// "***PWLF1 IS A ROUTINE WHICH WILL FIT UP TO 10(CURRENTLY) FUNCTIONS  "
		// "   SIMULTANEOUSLY ON AN INTERVAL (XL,XU) OF THE INDEPENDENT VARIABLE"
		// "   X OF THE FUNCTIONS. THE FIT IS A PIECEWISE LINEAR FUNCTION OF    "
		// "   XFUN(X).  XFI IS THE INVERSE FUNCTION OF XFUN.  THE SUBINTERVALS "
		// "   ARE CHOSEN OF UNIFORM WIDTH IN XFUN(X) AND SUFFICIENT OF THEM    "
		// "   ARE CHOSEN SO THAT THE FIT GIVES A RELATIVE ERROR[EP FOR ALL     "
		// "   THE FUNCTIONS OVER ALL THE SUBINTERVALS.                         "
		// "   QFIT IS AN AUXILIARY FUNCTION.       "
		// "   EXPLANATION OF THE ARGUMENTS:                                    "
		// "   NI   ON RETURN IS NUMBER OF SUBINTERVALS USED FOR THE FIT.       "
		// "   XL   LOWER LIMIT OF INTERVAL ON WHICH TO FIT THE FUNCTIONS.      "
		// "   XU   UPPER LIMIT                                                 "
		// "   XR   VALUE OF X WHICH WILL BE FORCED TO BE A SUBINTERVAL BOUNDARY"
		// "        THE SIGNIFICANCE OF THIS IS THAT THE STRAIGHT LINES ON THE  "
		// "        SUBINTERVALS ARE CHOSEN TO FIT EXACTLY AT THE SUBINTERVAL   "
		// "        BOUNDARIES, THUS IF ONE WANTS AN EXACT FIT OF THE FUNCTIONS "
		// "        AT A PARTICULAR POINT, XR SHOULD BE SET TO THAT VALUE.      "
		// "        OTHERWISE XR SHOULD BE SET TO XH.  ANOTHER REQUIREMENT      "
		// "        IS THAT XU SHOULD BE LARGER THAN XL.                        "
		// "   EP   THE MAXIMUM RELATIVE ERROR ALLOWED THE FIT.                 "
		// "   NIP  THE MINUIMUM NUMBER OF POINTS INTERIOR TO (XL,XU) AT WHICH  "
		// "        THE FIT IS TO BE TESTED FOR RELATIVE ERROR VS. EP.          "
		// "   XFUN A FUNCTION OF X OVER WHICH IT IS HOPED THE FUNCTIONS TO BE  "
		// "        FIT ARE MORE LINEAR THAN OVER X.  XFUN IS EXPECTED TO BE    "
		// "        MONOTONICALLY INCREASING IN X.                              "
		// "   XFI  THE INVERSE OF XFUN. THAT IS XFI(XFUN(X))=X.                "
		// "   AX,BX ARE COEFFICIENTS USED AS SHOWN BELOW TO DETERMINE WHICH    "
		// "        SUBINTERVAL A VALUE OF X IS IN.                             "
		// "   AF,BF ARE ARRAYS OF COEFFICIENTS USED TO GET VALUES OF THE FUNCS."
		// "   THE PROCEDURE FOR FINDING THE FIT VALUE OF FUNCTION IFUN IS:     "
		// "   INTERV=AX*XFUN(X)+BX                                             "
		// "   VALUE=AF(INTERV,IFUN)*XFUN(X)+BF(INTERV,IFUN)                    "
		// "   NALM  IS THE MAXIMUM NUMBER OF SUBINTERVALS FOR WHICH ARRAY SPACE"
		// "         HAS BEEN ALLOCATED.                                        "
		// "   NFUN  IS THE NUMBER OF FUNCTIONS TO BE FITTED(SIMULTANEOUSLY,I.E."
		// "         ALL FUNCTIONS HAVE THE SAME XFUN AND SUBINTERVALS, AND ALL "
		// "         ARE REQUIRED TO BE FIT WITH MAX REL ERR[EP)                "
		// "   VFUNS IS A SUBROUTINE TO FILL AN ARRAY WITH THE VALUES OF THE    "
		// "   FUNCTIONS TO BE FITTED.                                          "

		// "   QFIT IS A LOGICAL FUNCTION WHICH IS TRUE IF THE STATED NUMBER    "
		// "   OF INTERVALS GIVES A SUFFICIENTLY CLOSE FIT.                     "

		// "   FIND # OF INTERVALS REQUIRED.                                    "
		int NL = 0;
		int NU = 1;
		//int IPRN = 0;
		int NJ = 0;
		int NK = 0;// double REM=0.;
		// LOOP
		while (true) {
			NJ = Math.min(NU, NIMX);// nimx=150
			if (QFIT(NJ, XL, XU, XR, EP, ZTHR, ZEP, REM, NIP, XFUN, XFI, AX,
					BX, NALM, NFUN, AF, BF, VFUNS, 0))// here is set NJj=NJ!!
				break;// exit the loop
			if (NU >= NIMX) {
				// OUTPUT NIMX,EP;
				// (' NUMBER OF ALLOCATED INTERVALS(=',I5,') WAS INSUFFICIENT'
				// ,/ ,' TO GET MAXIMUM RELATIVE ERROR LESS THAN ',1P,G14.6);
				// NI=NJ;
				System.out
						.println("NUMBER OF ALLOCATED INTERVALS WAS INSUFFICIENT"
								+ " TO GET MAXIMUM RELATIVE ERROR LESS THAN MAXIMUM ALLOWED");
				// PLWF_NI_set(NI,NJ);
				PLWF_NI_set(NI, NJj);
				return;
			}
			NL = NU;
			NU = NU * 2;
		}// System.out.println("enter");
			// "   WE NOW HAVE AN UPPER AND LOWER LIMIT ON NI, REFINE IT.           "
			// NU=NJ ; //"SAVE SUCCESSFUL INDEX"
		NU = NJj; // "SAVE SUCCESSFUL INDEX"
		while (NU > NL + 1) {// "LOOP UNTIL CONVERGENCE"
			NJ = (NL + NU) / 2;
			NK = NJ;
			// NK=NJ;
			// //"THIS IS NECESSARY BECAUSE QFIT MAY LOWER NJ,BUT NEED ORIGINAL"
			// " FOR SETTING NL OR MAY GET INTO INFINITE LOOP."
			if (QFIT(NJ, XL, XU, XR, EP, ZTHR, ZEP, REM, NIP, XFUN, XFI, AX,
					BX, NALM, NFUN, AF, BF, VFUNS, 0)) {
				// NU=NJ;
				NU = NJj;
			} else {
				NL = NK;// original NJ!!
			}
		}
		// "     NU IS NOW THE SMALLEST NI WHICH FITS OK.                       "
		// NI=NU;
		PLWF_NI_set(NI, NU);
		// if (NI==NJ)
		// if (NU==NJ)
		if (NU == NJj)
			return;// "LAST TEST WAS SUCCESS"
		// "     CALL IT ONCE MORE TO GET THE FITS.                             "
		// if (!QFIT(NI,XL,XU,XR,EP,ZTHR,ZEP,REM,NIP,XFUN,XFI,
		if (!QFIT(NU, XL, XU, XR, EP, ZTHR, ZEP, REM, NIP, XFUN, XFI, AX, BX,
				NALM, NFUN, AF, BF, VFUNS, 0)) {
			System.out.println("CATASTROPHE---DOES NOT FIT WHEN IT SHOULD");
			// OUTPUT NI;
			// (' CATASTROPHE---DOES NOT FIT WHEN IT SHOULD,NI=',I5);
		}
		return;
	}

	// el->NJ=150; NALM=150; NFUN=8; NJP=20
	/**
	 * Called by PWLF1. CONSTRUCT THE INTERVAL COEFFICIENTS. XR SHOULD BE IN THE INTERVAL (XL,XH).  IF NOT IT WILL BE SE TO THE 
	 * NEAREST LIMIT.  SUBINTERVALS WILL BE ARRANGED SO THAT XR IS ALWAYS ON A SUBINTERVAL BOUNDARY. THE PURPOSE OF THIS 
	 * FEATURE IS TO MORE EASILY FIT FUNCTIONS WHICH HAVE AN INTERIOR DISCONTINUITY IN SLOPE.  EXAMPLES ARE THE MOLLER AND PAIR 
	 * CROSS SECTIONS WHICH CONTRIBUTE DISCONTINUITIES IN SLOPE TO THE ELECTRON AND PHOTON INTERACTION PROBABILITIES IN 
	 * THE INTERIOR OF THE ENERGY RANGES FOR THESE PARTICLES. <p>
	 * IN ABLE TO GIVE SOME VALUE FOR X'S WHICH MAY LIE SLIGHTLY OUTSIDE THE INTERVAL (XL,XH) AN EXTRA SUBINTERVAL ON EACH SIDE OF (XL,XH) 
	 * IS PROVIDED, WHICH USED THE SAME STRAIGHT LINES AS THE ADJACENT INCLUDED SUBINTERVAL.  NJ IS THE TOTAL NUMBER OF SUBINTERVALS 
	 * AND NI IS DEFINED TO BE THE NUMBER OF INTERNAL SUBINTERVALS=NJ-2. NJP IS THE MINUMUM NUMBER OF POINTS INTERIOR TO THE INTERVAL 
	 * (XL,XH) AT WHICH THE FIT IS TO BE TESTED.  A NUMBER NIP WILL BE CHOSEN AS THE NUMBER OF INTERIOR POINTS WITHIN EACH SUBINTERVAL 
	 * AT WHICH TO TEST SO THAT NIP*NI=NJP.
	 * @param NJ NJ
	 * @param XL XL
	 * @param XH XH
	 * @param XR XR
	 * @param EP EP
	 * @param ZTHR ZTHR
	 * @param ZEP ZEP
	 * @param REM REM
	 * @param NJP NJP
	 * @param XFUN XFUN
	 * @param XFI XFI
	 * @param AX AX
	 * @param BX BX
	 * @param NALM NALM
	 * @param NFUN NFUN
	 * @param AF AF
	 * @param BF BF
	 * @param VFUNS VFUNS
	 * @param IPRN IPRN
	 * @return the result
	 */
	public static boolean QFIT(int NJ, double XL, double XH, double XR,
			double EP, double[] ZTHR, double[] ZEP, double REM, int NJP,
			String XFUN, String XFI, String AX, String BX, int NALM, int NFUN,
			String AF, String BF, String VFUNS, int IPRN) {
		// $REAL AF(NALM,NFUN),BF(NALM,NFUN),ZTHR(NFUN),ZEP(NFUN);
		// "     CONSTRUCT THE INTERVAL COEFFICIENTS.                             "
		// "     XR SHOULD BE IN THE INTERVAL (XL,XH).  IF NOT IT WILL BE SE TO TH"
		// "     NEAREST LIMIT.  SUBINTERVALS WILL BE ARRANGED SO                 "
		// "     THAT XR IS ALWAYS ON A SUBINTERVAL BOUNDARY.   THE PURPOSE OF THI"
		// "     FEATURE IS TO MORE EASILY FIT FUNCTIONS WHICH HAVE AN INTERIOR   "
		// "     DISCONTINUITY IN SLOPE.  EXAMPLES ARE THE MOLLER AND PAIR        "
		// "     CROSS SECTIONS WHICH CONTRIBUTE DISCONTINUITIES IN SLOPE         "
		// "     TO THE ELECTRON AND PHOTON INTERACTION PROBABILITIES IN          "
		// "     THE INTERIOR OF THE ENERGY RANGES FOR THESE PARTICLES.           "
		// "     IN ABLE TO GIVE SOME VALUE FOR X'S WHICH MAY LIE SLIGHTLY OUTSIDE"
		// "     THE INTERVAL (XL,XH) AN EXTRA SUBINTERVAL ON EACH SIDE OF (XL,XH)"
		// "     IS PROVIDED, WHICH USED THE SAME STRAIGHT LINES AS THE ADJACENT  "
		// "     INCLUDED SUBINTERVAL.  NJ IS THE TOTAL NUMBER OF SUBINTERVALS    "
		// "     AND NI IS DEFINED TO BE THE NUMBER OF INTERNAL SUBINTERVALS=NJ-2."
		// "     NJP IS THE MINUMUM NUMBER OF POINTS INTERIOR TO THE INTERVAL     "
		// "     (XL,XH) AT WHICH THE FIT IS TO BE TESTED.  A NUMBER NIP WILL BE  "
		// "     CHOSEN AS THE NUMBER OF INTERIOR POINTS WITHIN EACH SUBINTERVAL  "
		// "     AT WHICH TO TEST SO THAT NIP*NI]=NJP.                        "

		boolean QFIT = false;
		func = XFUN;
		funci = XFI;
		funceval = VFUNS;
		double d = 0.;
		// --------------------------
		double XFL = 0.;
		double XLL = 0.;
		int JSUB = 0;
		double SXFH = 0.;
		double DSXF = 0.;
		double WIP = 0.;
		double SXFIP = 0.;
		double XIP = 0.;
		// ------------------------------
		// int NKP=3;
		if (XH <= XL) {
			// OUTPUT XL,XH;(' QFIT ERROR:XL SHOULD BE < XH. XL,XH=',2G14.6);
			QFIT = false;
			return QFIT;
		}
		double XS = Math.max(XL, Math.min(XH, XR));// allways XR if: XL<=XR<=XH
		// else if XR<XL=>XL or if XR>XH=>XH!!
		// "     GET NUMBER OF INTERNAL SUBINTERVALS ARE ALLOWED AND CHECK        "
		int NI = NJ - 2;// NJ IS THE TOTAL NUMBER OF SUBINTERVALS maxim 150=>148
		// "     AT LEAST 2 SUBINTERVALS ARE NEEDED IF XR(XS) IS NOT AN END POINT."
		if ((((XS == XL) || (XS == XH)) && (NI >= 1)) || (NI >= 2)) {
			XFL = XFUN(XL);
		} else {
			QFIT = false;
			return QFIT;
		}
		double XFH = XFUN(XH);
		double XFS = XFUN(XS);
		// "     SET SUBINTERVAL WIDTH.                                           "
		double XM = Math.max(XFH - XFS, XFS - XFL);
		double DX = XFH - XFL;
		double W = XM / Math.max(1., AINT(NI * XM / DX));// usually is DX/NI but
															// put XS on
															// boundary!!
		// "     RESET NI TO HOW MANY WE'RE ACTUALLY GOING TO USE.                "
		NI = NI - AINT(NI - DX / W);
		// "     COMPUTE HOW MANY INTERIOR POINTS TO SAMPLE IN EACH SUBINTERVAL.  "
		// NJP IS THE MINUMUM NUMBER OF POINTS INTERIOR TO THE INTERVAL "
		// (XL,XH) AT WHICH THE FIT IS TO BE TESTED
		// NIP=NUMBER OF INTERIOR POINTS WITHIN EACH SUBINTERVAL
		// 20+148-1/148=>3
		int NIP = Math.max(NKP, (NJP + NI - 1) / NI);// minim
														// 3@@@@@@@@@@@@@@@@@@@@@@@
		// int NIP=Math.max(NKP,(NJP)/NI);//minim 3@@@@@@@@@@@@@@@@@@@@@@@
		// "     MAKE NIP ODD                                                     "
		NIP = (NIP / 2) * 2 + 1;// =>3
		// "     SET ACTUAL LOWER LIMIT OF INTERVAL.                              "
		if ((XFH - XFS) <= (XFS - XFL)) {
			XLL = XFL;// System.out.println("equal");
		} else {
			XLL = XFH - NI * W;// System.out.println(XLL+" vs XFL "+XFL);
		}
		// "     COEFICIENTS FOR USER TO COMPUTE WHICH SUBINTERVAL TO USE.        "
		// "     ISUBINT=AX*XFUN(X)+BX                                            "
		// AX=1./W;
		double axd = 1. / W;
		PLWF_AX_set(AX, axd);
		// BX=2.-XLL*AX;
		double bxd = 2. - XLL * axd;
		PLWF_BX_set(BX, bxd);
		// for a given value vx let xvx be the XFUN(vx), so ax*xvx+bx is the
		// same with
		// (xvx-xll)Ni/Dx+2, Ni/Dx = 1/w=ax , so when xvx is xll =>interv=2 and
		// in java we will
		// take interv-1 =1.The int value of (ax*xvx+bx) make all xvx>xll and
		// xvx<xll+w
		// to generate a value of 2=>interval=1.when xvx=xfh=> could be
		// Ni+2 so interv=Ni+1!Note Ni=Nj-2-> so we take the last [Nj-1] array 0
		// biased!!!
		// in Fortran for a given I (eg 2)=>
		// "     THE RIGHT BOUNDARY OF SUBINTERVAL I IS XFI(XLL+W*(I-1))          "
		// "     NOW COMPUTE THE FIT COEFFICIENTS FOR THE SUBINTERVALS            "
		// "     AND FIND MAXIMUM RELATIVE ERROR(REM).                            "
		// in Java a interval I = aint(ax*xvx+bx)-1.so the right boundary is
		// XFI(XLL+W*I)!!!
		REM = 0.0;
		QFIT = true;
		// "     LOOP OVER SUBINTERVALS                                           "
		// "     INITIALIZE LOWER BOUNDARY AND VALUE.                             "
		double SXFL = Math.max(XLL, XFL);
		// **************************doing so, if XLL<XFL could missed XS
		// point!!
		// double SXFL=0.;//XLL;
		// if (XLL>=XFL)
		// {
		// SXFL=XLL;
		// }
		// else
		// {
		// while(XLL<XFL)
		// {
		// XLL=XLL+W;//System.out.println("fit "+W);
		// }
		// SXFL=XLL;
		// }
		// double axd=1./W;
		// PLWF_AX_set(AX, axd);
		// double bxd=2.-XLL*axd;//see below ##!!
		// PLWF_BX_set(BX, bxd);
		// ********the above are not necessary!!!

		int ISUB = 0;
		double XSXF = XFI(SXFL);
		double[] FSXL = new double[10];
		double[] FSXH = new double[10];
		double[] FIP = new double[10];
		double[] FFIP = new double[10];
		double[] AFIP = new double[10];
		double[] RE = new double[10];
		double[] AER = new double[10];
		VFUNS(XSXF, FSXL);
		// if (IPRN!=0)
		// WRITE(6,:FMT:) ISUB,SXFL,XSXF,(FSXL(IFUN),IFUN=1,NFUN);
		// :FMT: FORMAT(' QFIT:ISUB,SXF,XSXF,FSX()=',I4,1P,9G11.4/(1X,12G11.4));
		for (ISUB = 1; ISUB <= NI; ISUB++)// ->148
		{
			// "     ALLOW FOR EXTRA SUBINTERVAL OUTSIDE THE MAIN INTERVAL            "
			JSUB = ISUB + 1;// ->2 to 149 for instance####################
			// start from 2->AF[2-1] 0 biased with [0] extrapolation,
			// so JSUB=2 is minimum and XFUN=XLL=>we want interval=AX*XFUN+BX
			// and then
			// value=AF[interval-1][]*XFUN+BF....so 2=ax*xll+bx=>bx=2.-xll*ax
			// so interval=AINT(AX*XFUN+BX);!!!!!!!!!!!!!
			// So for interv=2 we go for AF[1] and so on!!!
			SXFH = Math.min(XLL + W * ISUB, XH);// right boundary <=>
												// w*ISUB=w*(JSUB-1)
			// System.out.println("actual "+SXFH+" sing= " +XFS);//OK
			XSXF = XFI(SXFH);
			// CALL VFUNS(XSXF,FSXH);
			VFUNS(XSXF, FSXH);
			if (IPRN != 0) {
				// WRITE(6,:FMT:)ISUB,SXFH,XSXF,(FSXH(IFUN),IFUN=1,NFUN);
			}
			DSXF = SXFH - SXFL;
			// DO IFUN=1,NFUN
			for (int IFUN = 0; IFUN < NFUN; IFUN++) {
				d = (FSXH[IFUN] - FSXL[IFUN]) / DSXF;
				PLWF_AF_set(AF, d, JSUB - 1, IFUN);// min 1,right boundary <=>
													// w*ISUB=w*(JSUB-1)
				// AF(JSUB,IFUN)=(FSXH(IFUN)-FSXL(IFUN))/DSXF;
				d = (FSXL[IFUN] * SXFH - FSXH[IFUN] * SXFL) / DSXF;
				PLWF_BF_set(BF, d, JSUB - 1, IFUN);// min 1
				// BF(JSUB,IFUN)=(FSXL(IFUN)*SXFH-FSXH(IFUN)*SXFL)/DSXF;
			} // "END OF IFUN"
				// "     LOOP OVER INTERIOR POINTS TO LOOK FOR MAX. REL.ERROR             "
				// "     COMPUTE INTERIOR POINT SPACING.                                  "
			WIP = DSXF / (NIP + 1);// ok
			// DO IP=1,NIP[
			for (int IP = 1; IP <= NIP; IP++) {
				// VALUE OF XFUN AT THE INTERIOR POINT OF THIS SUBINTERVAL. "
				SXFIP = SXFL + IP * WIP;// ok
				XIP = XFI(SXFIP);
				// COMPUTE FUNCTION AT INTERIOR POINT "
				// CALL VFUNS(XIP,FIP);
				VFUNS(XIP, FIP);
				// "COMPUTE FITTED VALUES.                                           "
				// DO IFUN=1,NFUN[
				for (int IFUN = 0; IFUN < NFUN; IFUN++) {
					// FFIP(IFUN)=AF(JSUB,IFUN)*SXFIP+BF(JSUB,IFUN);
					FFIP[IFUN] = PLWF_AF_value(AF, JSUB - 1, IFUN) * SXFIP
							+ PLWF_BF_value(BF, JSUB - 1, IFUN);
					AFIP[IFUN] = Math.abs(FIP[IFUN]);
					AER[IFUN] = Math.abs(FFIP[IFUN] - FIP[IFUN]);
					RE[IFUN] = 0.0;
					if (FIP[IFUN] != 0.0) {
						RE[IFUN] = AER[IFUN] / AFIP[IFUN];
						// if(RE[IFUN]>EP)
						// System.out.println(" isub "+ISUB+" i "+IFUN+"  er "+RE[IFUN]+"  vs "+EP+" ni= "+NI);
					}
					if (AFIP[IFUN] >= ZTHR[IFUN]) {
						REM = Math.max(REM, RE[IFUN]);
					} else if (AER[IFUN] > ZEP[IFUN]) {
						QFIT = false;// System.out.println("enter");
					}
				}// "END OF IFUN"
					// "*****WRITE OUT SO WE CAN SEE HOW WE ARE DOING.                        "
				if (IPRN != 0) {
					// OUTPUT ISUB,IP,SXFIP,XIP,REM,QFIT,(FIP(IFUN),FFIP(IFUN),
					// RE(IFUN),AER(IFUN),IFUN=1,NFUN);
					// (1X,2I4,1P,2G12.5,6P,F12.0,L2,1P,2G11.4,6P,F11.0,1P,G11.4/
					// (1X,3(1P,2G11.4,6P,F11.0,1P,G11.4)));
				}
			} // "END OF IP"
				// "     SAVE RIGHT BOUNDARY AND VALUE FOR NEXT SUBINTERVAL.              "
			SXFL = SXFH;
			// DO IFUN=1,NFUN [FSXL(IFUN)=FSXH(IFUN);]
			for (int IFUN = 0; IFUN < NFUN; IFUN++) {
				FSXL[IFUN] = FSXH[IFUN];
			}
		}// "END OF ISUB"
			// "     SET UP SKIRTING SUBINTERVALS(VIA EXTRAPOLATION)                  "
			// DO IFUN=1,NFUN
		for (int IFUN = 0; IFUN < NFUN; IFUN++) {
			d = PLWF_AF_value(AF, 1, IFUN);
			PLWF_AF_set(AF, d, 0, IFUN);// fill 0
			// AF(1,IFUN)=AF(2,IFUN);
			d = PLWF_BF_value(BF, 1, IFUN);
			PLWF_BF_set(BF, d, 0, IFUN);// fill 0
			// BF(1,IFUN)=BF(2,IFUN);
			d = PLWF_AF_value(AF, NI, IFUN);
			PLWF_AF_set(AF, d, NI + 1, IFUN);// fill last
			// AF(NI+2,IFUN)=AF(NI+1,IFUN);
			d = PLWF_BF_value(BF, NI, IFUN);
			PLWF_BF_set(BF, d, NI + 1, IFUN);// fill last
			// BF(NI+2,IFUN)=BF(NI+1,IFUN);
		} // "END OF IFUN"
			// QFIT=QFIT&(REM<=EP);
		if ((QFIT) && (REM <= EP))
			QFIT = true;
		else
			QFIT = false;

		// NJ=NI+2; //"TELL ACTUAL NO. OF SUBINTERVALS USED."
		NJj = NI + 2;
		return QFIT;
	}

	/**
	 * Internally used by QFIT
	 * @param AFES AFES
	 * @param i i
	 * @param j j
	 * @return the result
	 */
	private static double PLWF_AF_value(String AFES, int i, int j) {
		double r = 0.;
		if (AFES.compareTo("AFE") == 0)
			r = AFE[i][j];
		if (AFES.compareTo("AFG") == 0)
			r = AFG[i][j];
		if (AFES.compareTo("AFR") == 0)
			r = AFR[i][j];

		return r;
	}

	/**
	 * Internally used by QFIT
	 * @param BFES BFES
	 * @param i i
	 * @param j j
	 * @return the result
	 */
	private static double PLWF_BF_value(String BFES, int i, int j) {
		double r = 0.;
		if (BFES.compareTo("BFE") == 0)
			r = BFE[i][j];
		if (BFES.compareTo("BFG") == 0)
			r = BFG[i][j];
		if (BFES.compareTo("BFR") == 0)
			r = BFR[i][j];

		return r;
	}

	/**
	 * Internally used by QFIT
	 * @param AFES AFES
	 * @param value value
	 * @param i i
	 * @param j j 
	 
	 */
	private static void PLWF_AF_set(String AFES, double value, int i, int j) {
		if (AFES.compareTo("AFE") == 0)
			AFE[i][j] = value;
		if (AFES.compareTo("AFG") == 0)
			AFG[i][j] = value;
		if (AFES.compareTo("AFR") == 0)
			AFR[i][j] = value;

	}

	/**
	 * Internally used by QFIT
	 * @param BFES BFES
	 * @param value value
	 * @param i i
	 * @param j j 
	
	 */
	private static void PLWF_BF_set(String BFES, double value, int i, int j) {
		if (BFES.compareTo("BFE") == 0)
			BFE[i][j] = value;
		if (BFES.compareTo("BFG") == 0)
			BFG[i][j] = value;
		if (BFES.compareTo("BFR") == 0)
			BFR[i][j] = value;

	}

	/**
	 * Internally used by PWLF1
	 * @param NIS NIS
	 * @param value value
	 */
	private static void PLWF_NI_set(String NIS, int value) {
		if (NIS.compareTo("NEL") == 0)
			NEL = value;
		if (NIS.compareTo("NGL") == 0)
			NGL = value;
		if (NIS.compareTo("NGR") == 0)
			NGR = value;

	}

	/**
	 * Internally used by QFIT
	 * @param AXS AXS
	 * @param value value
	 */
	private static void PLWF_AX_set(String AXS, double value) {
		if (AXS.compareTo("AXE") == 0)
			AXE = value;
		if (AXS.compareTo("AXG") == 0)
			AXG = value;
		if (AXS.compareTo("AXR") == 0)
			AXR = value;

	}

	/**
	 * Internally used by QFIT
	 * @param BXS BXS
	 * @param value value
	 */
	private static void PLWF_BX_set(String BXS, double value) {
		if (BXS.compareTo("BXE") == 0)
			BXE = value;
		if (BXS.compareTo("BXG") == 0)
			BXG = value;
		if (BXS.compareTo("BXR") == 0)
			BXR = value;

	}

	/**
	 * Convert a double to an int
	 * @param d d
	 * @return the result
	 */
	public static int AINT(double d) {
		Double r = new Double(d);
		return r.intValue();
	}

	/**
	 * Handle XFUN for PWLF1
	 * @param X X
	 * @return the result
	 */
	public static double XFUN(double X) {
		double Y = 0.;
		if (func.compareTo("ALKE") == 0)
			Y = PEGS4.ALKE(X);
		else if (func.compareTo("ALOG") == 0)
			Y = Math.log(X);
		else if (func.compareTo("ALIN") == 0)
			Y = PEGS4.ALIN(X);

		return Y;
	}

	/**
	 * Handle XFI for PWLF1
	 * @param X X
	 * @return the result
	 */
	public static double XFI(double X) {
		double Y = 0.;
		if (funci.compareTo("ALKEI") == 0)
			Y = PEGS4.ALKEI(X);
		else if (funci.compareTo("EXP") == 0)
			Y = Math.exp(X);
		else if (funci.compareTo("ALINI") == 0)
			Y = PEGS4.ALINI(X);

		return Y;
	}

	// public static double[] VFUNS(double X, double [] Z)
	/**
	 * Computes some functions for PWLF1.
	 * @param X X
	 * @param Z Z
	 */
	public static void VFUNS(double X, double[] Z) {
		// double [] Y=new double[10];
		if (funceval.compareTo("EFUNS") == 0) {
			// Y=new double[10];//in fact 8
			// PEGS4.EFUNS(X,Y);
			PEGS4.EFUNS(X, Z);
		}
		if (funceval.compareTo("GFUNS") == 0) {
			// Y=new double[10];//in fact 4
			// PEGS4.GFUNS(X,Y);
			PEGS4.GFUNS(X, Z);
		}
		if (funceval.compareTo("RFUNS") == 0) {
			// Y=new double[10];//in fact 1
			// PEGS4.RFUNS(X,Y);
			PEGS4.RFUNS(X, Z);
		}
		// return Y;
	}
}
