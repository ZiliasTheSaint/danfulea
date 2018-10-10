package danfulea.phys.egs;

import java.io.FileInputStream;
import java.math.BigInteger;
import java.util.Vector;

/**
 * 
 * Prepare medium input data for EGS (Electron Gamma Shower). 
 * Based on EGS4/EGSnrc routines developed by SLAC(Stanford Linear Accelerator Center)/NRC (National Research Council of Canada).
 * DESCRIPTION: PEGS4 ('Preprocessor for EGS4') is a stand-alone program whose purpose is to generate material data for the EGS4 
 * code, and to provide other services for the user who is studying the simulation of electromagnetic interactions.
 * @author Dan Fulea 28 JUL. 2005 
              
 */

public class PEGS4 {
	public static final String datas = "Data";
	public static final String egsData = "egs";
	public static final String compoundss = "compounds";
	public static final String elementss = "elements";

	public static String userdirs = System.getProperty("user.dir");
	public static String file_sep = System.getProperty("file.separator");
	public static Vector<String> zdepv = new Vector<String>();
	public static int npairdelta = 0;
	public static double[] edelta;// energy
	public static double[] delta;// delta correction
	public static int jzet = 0;
	public static int jrhozet = 0;
	public static int jenergy = 0;
	public static int jdelta = 0;
	public static String medium = "";
	public static String MTYP = "";
	public static String densPath = "";
	public static String crossPath = "";
	public static String formPath = "";
	// ------------------------------------------------------------------
	public static double $RERR = 1.E-5; // "RERR-VALUE NEEDED BY DCADRE"
	public static double $AERR = 1.E-16; // "AERR-VALUE NEEDED BY DCADRE"
	// --------------------------READ------------------------------------
	// -----PHYSICAL AND MATHEMATICAL
	// CONSTANTS------------------------------------------------
	public static final double PI = Math.PI;// =CIRCUMFERANCE/DIAMETER.
	public static final double C = 2.997925E+10;// SPEED OF LIGHT(CM/SEC)
	public static final double RME = 9.1091E-28;// ELECTRON REST MASS(GRAMS)
	public static final double HBAR = 1.05450E-27;// PLANCK'S CONSTANT/(2*PI)
													// (ERG SEC)
	public static final double ECGS = 4.80298E-10;// ELECTRON CHARGE (ESU or
													// u.e.s. C.G.S)
	public static final double EMKS = 1.60210E-19;// ELECTRON CHARGE (COULOMB)
	public static final double AN = 6.02252E+23;// AVOGADRO'S NUMBER-mole**-1
	// ----------------------------------------------------------------------------------------
	// --------DERIVED CONSTANTS
	// RESOLVED-----------------------------------------------------
	public static double RADDEG = 0.0;// ONE RADIAN IN DEGREES
	public static double FSC = 0.0;// FINE STRUCTURE CONSTANT ~1/137 =0.00729927
	public static double ERGMEV = 0.0;// ONE MILLION ELECTRON VOLTS EXPRESSED IN
										// ERGS
	public static double R0 = 0.0;// CLASSICAL ELECTRON RADIUS IN CM
	public static double RM = 0.0;// ELECTRON REST ENERGY IN MEV
	public static double RMT2 = 0.0;// 2*RM
	public static double RMSQ = 0.0;// RM*RM
	public static double A22P9 = 0.0;// NAGEL HAD THE ANGLE 22.9 DEGREES IN HIS
										// EXPRESSION FOR THE "
										// CHARACTERISTIC ANGLE,CHIC
	public static double A6680 = 0.0;// IN ADDITION NAGEL USES THE NUMBER 6680
										// IN HIS EXPRESSION FOR "
										// EXP(B) AS OBTAINED BY BETHE(1953),PG
										// 1259,EQN.229
	// ----------------------------------------------------------------------------------------
	// ----------------- DATA FOR ALRAD AND ALRADP FOR ELEMENT Z.LE.4--(Lower or
	// equal)---------
	// ---------- TAKEN FROM TABLE B.2 IN Y.TSAI REV.MOD.PHYS.
	// 46,815(1974)--------------------
	public static double[] ALRAD = { 5.31, 4.79, 4.74, 4.71 };
	public static double[] ALRADP = { 6.144, 5.621, 5.805, 5.924 };
	public static double A1440 = 1194.0;
	public static double A183 = 184.15;
	// -----------------------------------------------------------------------------------------
	// "------------ MIXDAT --line 456 in
	// mortran-----------------------------------------------
	public static double GASP = 0.0;// GASP=GAS PRESSURE=0.0 MEANS 'NON-GAS'
									// STATE"
	// Defines state of mixture: zero (default) for solid or liquid, otherwise
	// value gives gas pressure (atm).
	// -------------MOLECULAR VARIABLES
	// -------------------------------------------------------
	// @@@@@@@@@@@@@@@@@@@@@@@EXTERNAL SET
	public static int NE = 0;// NUMBER OF DIFFERENT TYPES OF ATOMS IN THE
								// MATERIAL
	public static double[] PZ;// PROPORTION OF ELEMENT OF TYPE I. IF A
								// COUMPOUND, "
	// " THEN PZ(I) WILL BE THE NUMBER OF ATOMS OF TYPE I IN THE MOLECULE."
	// " IF A MIXTURE,SUCH AS CONCRETE, PZ(I) COULD BE THE PER CENT OF    "
	// " THE ATOMS WHICH ARE OF TYPE I.
	public static double[] Z;// PERIODIC NUMBER OF ATOMS OF TYPE I
	public static double[] WA;// ATOMIC WEIGHT FOR ATOMS OF TYPE I
	public static double[] RHOZ;// PARTIAL DENSITY DUE TO ATOMS OF TYPE I.
								// (GM/CM**3) "
	// " ELECTRON DENSITY VARIABLE
	public static double RHO = 0.0;// DENSITY OF THE MATERIAL. (IN GRAMS/CM**3)
	// @@@@@@@@@@@@@@@@@@@@END EXTERNAL SET-----------
	public static double WM = 0.0;// SUM(PZ(I)*WA(I)) = MOLECULAR WEIGHT IF A
									// COUMPOUND "
	// " OR A 'MIXTURE WEIGHT' IF A MIXTURE.
	public static double ZC = 0.0;// =SUM(PZ(I)*Z(I)) = NUMBER OF
									// ELECTRONS/MOLECULE "
	// "          BREMSSTRAHLUNG AND PAIR PRODUCTION VARIABLES ARE WEIGHTED"
	// " BY PZ(I)*Z(I)**2 FOR THE NUCLEUS, AND BY PZ(I)*Z(I)*XSI(I) FOR   "
	// " ATOMIC ELECTRONS.
	public static double TPZ = 0.0;// =SUM(PZ(I))
	// @@@ array computed in MIX:
	public static double[] XSI;// =LOG(A1440/Z(I)**(2./3.))/(LOG(A183/Z(I)**(1./3.))-
								// FCOUL(Z(I)) )
	public static double[] ZZX;// =PZ(I)*Z(I)*(Z(I)+XSI(I)) = BREMS AND PAIR
								// WEIGHTS
	public static double[] FZC;//
	public static double[] FCOUL;//
	// @@end array computed in MIX
	public static double EZ = 0.0;// =ZC/TPZ EFFECTIVE Z
	public static double ZT = 0.0;// =SUM(ZZX(I))
	public static double ZA = 0.0;// =LOG(A183)*ZT BUTCHER AND MESSELS L.C.'A'
									// (1960)P.18
	public static double ZB = 0.0;// =SUM(ZZX(I)*LOG(Z(I)**(-1./3.) B&M'S
									// L.C.'B' IBID.
	public static double ZF = 0.0;// =SUM(ZZX(I)*FCOUL(Z(I))),WHERE FCOUL IS THE
									// COULOMB "
									// " CORRECTION FUNCTION.
	// ----END MOLECULAR VARIABLES ---------------------------------------
	// ---------------------RATIOS---------------------------------------------------------
	public static double ZG = 0.0;// = ZB/ZT ,EXP(ZG)=WEIGHTED GEOMETRIC MEAN OF
									// Z**(-1/3)
	public static double ZP = 0.0;// = ZB/ZA , B&M IBID.P18 L.C.'P'
	public static double ZV = 0.0;// = (ZB-ZF)/ZT
	public static double ZU = 0.0;// = (ZB-ZF)/ZA
	// --------------END
	// RATIOS--------------------------------------------------------------
	// -----------------------MULTIPLE SCATTERING
	// VARIABLES--------------------------------
	// " ACCORDING TO MESSEL AND CRAWFORD(1970), MOST OF THE MULTIPLE     "
	// " SCATTERING DUE TO THE FIELD OF THE ATOMIC ELECTRONS IS ALREADY   "
	// " ACCOUNTED FOR BY THE DISCRETE MOLLER SCATTERING. HENCE,THE       "
	// " FOLLOWING VARIABLES ARE ONLY WEIGHTED BY PZ(I)*Z(I)**2           "
	// " HOWEVER I HAVE NOT JUSTIFIED THE ABOVE ASSERTION THEORETICALLY   "
	// " THEORETICALLY, AND COMPARISON OF EXPERIMENTS WITH EGS HAVE       "
	// " SHOWN EGS PHOTON SPECTRUM DEFICIENT IN THE BACKWARD DIRECTION.   "
	// " THEREFORE, I WILL EXPERIMENT WITH RESTORING THE ELECTRON SCATTER-"
	// " ING TERM.  ITS CONTRIBUTION WILL BE REPRESENTED BY THE MACRO VAR-"
	// " IABLE, $FUDGEMS, WHICH WILL HAVE THE VALUE 0.0 FOR NO ELECTRON   "
	// " SCATTERING AND 1.0 FOR FULL ELECTRON SCATTERING.  AND WE NOW HAVE"
	// " WEIGHTING BY PZ(I)*Z(I)*(Z(I)+$FUDGEMS).                         "
	public static double ZAB = 0.0;
	public static double $FUDGEMS = 1.0; // "FULL MS OFF ATOMIC ELECTRONS"
	public static double[] ZZ;// = PZ(I)*Z(I)*(Z(I)+$FUDGEMS) "
	public static double ZS = 0.0;// = SUM(ZZ(I))
	public static double ZE = 0.0;// = SUM(ZZ(I)*LOG(Z(I)**(-2./3.)))
	public static double ZX = 0.0;// = SUM(ZZ(I)*LOG(1.+3.34*(FSC*Z(I))**2))
	public static double TEFF0 = 0.0;
	public static double XCC = 0.0;
	public static double XR0 = 0.0;
	public static double BLCC = 0.0;
	// ----------------------END MULTIPLE SCATTERING
	// VARIABLES--------------------------------
	// ---------------------ELECTON
	// DENSITY(ELECTRONS/CM**3)-------------------------------
	public static double EDEN = 0.0;// =AN*RHO/WM*ZC
	// --------------------END ELECTON
	// DENSITY(ELECTRONS/CM**3)---------------------------
	// RADIATION LENGTH see eq:2.7.22 SLAC 265 -verified error missed *4.0!! "
	// " USEFUL FOR GAUGING THE STEP SIZE, EVEN IF IT IS NOT USED AS THE  "
	// " UNIT OF DISTANCE.                                                "
	// "  1./RLC =(AN*RHO/WM)*4.0*FSC*R0**2*                              "
	// "    SUM( Z(I)*(Z(I)+XSI(I))*(LOG(A183*Z(I)**(-1./3.))-FCOUL(Z(I)))) "
	// =(AN*RHO/WM)*4.*FSC*R0**2*(ZAB-ZF)
	public static double RLC = 0.0;
	public static String[] ASYM;
	// ------------END ALL
	// MIX-------------------------------------------------------------
	public static int IMIXT = 0;// handle if we have elem,compound or mixture!!
	// 0 for elem and compound and 1 for mixture.
	// "*****TABLES OF ATOMIC SYMBOLS, WEIGHTS, DENSITIES, AND"
	// "MEAN EXCITATION ENERGIES FOR Z=1 TO NET(=100)."
	// "DATA FOR COMMON BLOCK ELEMTB"
	public static int NET = 100;
	public static double ZTBL = 0.0;
	public static String[] ASYMT = { "H", "HE", "LI", "BE", "B", "C", "N", "O",
			"F", "NE", "NA", "MG", "AL", "SI", "P", "S", "CL", "AR", "K", "CA",
			"SC", "TI", "V", "CR", "MN", "FE", "CO", "NI", "CU", "ZN", "GA",
			"GE", "AS", "SE", "BR", "KR", "RB", "SR", "Y", "ZR", "NB", "MO",
			"TC", "RU", "RH", "PD", "AG", "CD", "IN", "SN", "SB", "TE", "I",
			"XE", "CS", "BA", "LA", "CE", "PR", "ND", "PM", "SM", "EU", "GD",
			"TB", "DY", "HO", "ER", "TM", "YB", "LU", "HF", "TA", "W", "RE",
			"OS", "IR", "PT", "AU", "HG", "TL", "PB", "BI", "PO", "AT", "RN",
			"FR", "RA", "AC", "TH", "PA", "U", "NP", "PU", "AM", "CM", "BK",
			"CF", "ES", "FM" };
	public static double[] WATBL = { 1.00797, 4.0026, 6.939, 9.0122, 10.811,
			12.01115, 14.0067, 15.9994, 18.9984, 20.183, 22.9898, 24.312,
			26.9815, 28.088, 30.9738, 32.064, 35.453, 39.948, 39.102, 40.08,
			44.956, 47.90, 50.942, 51.998, 54.9380, 55.847, 58.9332, 58.71,
			63.54, 65.37, 69.72, 72.59, 74.9216, 78.96, 79.808, 83.80, 85.47,
			87.62, 88.905, 91.22, 92.906, 95.94, 99.0, 101.07, 102.905, 106.4,
			107.87, 112.4, 114.82, 118.69, 121.75, 127.60, 126.9044, 131.30,
			132.905, 137.34, 138.91, 140.12, 140.907, 144.24, 147., 150.35,
			151.98, 157.25, 158.924, 162.50, 164.930, 167.26, 168.934, 173.04,
			174.97, 178.49, 180.948, 183.85, 186.2, 190.2, 192.2, 195.08,
			196.987, 200.59, 204.37, 207.19, 208.980, 210., 210., 222., 223.,
			226., 227., 232.036, 231., 238.03, 237., 242., 243., 247., 247.,
			248., 254., 253. };
	public static double[] RHOTBL = { 0.0808, 0.19, 0.534, 1.85, 2.5, 2.26,
			1.14, 1.568, 1.5, 1.0, 0.9712, 1.74, 2.702, 2.4, 1.82, 2.07, 2.2,
			1.65, 0.86, 1.55, 3.02, 4.54, 5.87, 7.14, 7.3, 7.86, 8.71, 8.90,
			8.9333, 7.140, 5.91, 5.36, 5.73, 4.80, 4.2, 3.4, 1.53, 2.6, 4.47,
			6.4, 8.57, 9.01, 11.50, 12.20, 12.50, 12., 10.5, 8.65, 7.30, 7.31,
			6.684, 6.24, 4.93, 2.7, 1.873, 3.5, 6.15, 6.90, 6.769, 7.007, 1.,
			7.54, 5.17, 7.87, 8.25, 8.56, 8.80, 9.06, 9.32, 6.96, 9.85, 11.40,
			16.60, 19.30, 20.53, 22.48, 22.42, 21.45, 19.30, 14.19, 11.85,
			11.34, 9.78, 9.30, 1., 4., 1., 5., 1., 11.0, 15.37, 18.90, 20.5,
			19.737, 11.7, 7., 1., 1., 1., 1. };
	public static double[] ITBL = { 19.2, 41.8, 40., 63.7, 76.0, 78.0, 82.0,
			95.0, 115., 137., 149., 156., 166., 173., 173., 180., 174., 188.,
			190., 191., 216., 233., 245., 257., 272., 286., 297., 311., 322.,
			330., 334., 350., 347., 348., 357., 352., 363., 366., 379., 393.,
			417., 424., 428., 441., 449., 470., 470., 469., 488., 488., 487.,
			485., 491., 482., 488., 491., 501., 523., 535., 546., 560., 574.,
			580., 591., 614., 628., 650., 658., 674., 684., 694., 705., 718.,
			727., 736., 746., 757., 790., 790., 800., 810., 823., 823., 830.,
			825., 794., 827., 826., 841., 847., 878., 890., 902., 921., 934.,
			939., 952., 966., 980., 994. };
	public static int[] ISTATB = { 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0,
			0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
			0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
	// ----------------"DATA FOR COMMON BLOCK RAYLEI"-----------------------
	public static int IRAYL = 1;// 0 MEANS NOT INCLUDE RAYLEIGH SCATTERING"
								// 1 MEANS INCLUDE RAYLEIGH SCATTERING"
	// ----------------END DATA FOR COMMON BLOCK RAYLEI"-----------------------
	// ---------------------"DATA FOR COMMON BLOCK THRESH"-----------------------
	public static int IUNRST = 0; // "IUNRST=0 MEANS RESTRICTED STOPPING POWER"
	// "         IUNRST=1 MEANS UNRESTRICTED COLLISION STOPPING POWER"
	// "         IUNRST=2 MEANS CSDA DATA - TOTAL UNRESTRICTED STOPPING POWER"
	// "                  AND NO DISCRETE INTERACTIONS"
	// "         IUNRST=3 MEANS ALLOW BREM EVENTS BUT NO MOLLER INTERACTIONS"
	// "         IUNRST=4 MEANS ALLOW MOLLER EVENTS BUT NO BREM - NOTE THIS"
	// "                  GIVES COMPLETE GARBAGE IN A RUN SINCE ALL BREM IS "
	// "                  DUMPED ON SPOT"
	// "         IUNRST=5 MEANS UNRESTRICTED RADIATIVE STOPPING POWER"
	// ---------------------END DATA FOR COMMON BLOCK
	// THRESH"-----------------------
	// -----------------"DATA FOR COMMON BLOCK MIMSD, MULTIPLE SCATTER"-------------------------
	// PARAMETER $MXEL=20; "MAXIMUM NUMBER OF ELEMENTS IN MATERIAL"
	public static int $MXVRT1 = 1000; // "NUMBER OF REPRESENTATIVE ANGLES IN VERT1"
	public static int $MXJREFF = 200; // "SIZE OF MULTIPLE SCATTERING JREFF MAP"
	public static int $MSSTEPS = 16; // "NUMBER OF MULTIPLE SCATTERING STEP SIZES."
	public static int $MXVRT2 = 200; // "DISTRIBUTIONS OF NONOVERLAPPING PARTS OF VERT"
	// PARAMETER $MXSGE=400; "GAMMA SMALL ENERGY INTERVALS"
	public static int $MXGE = 200; // "GAMMA MAPPED ENERGY INTERVALS"
	// PARAMETER $MXSEKE=300; "ELECTRON SMALL ENERGY INTERVALS"
	public static int $MXEKE = 150; // "ELECTRON MAPPED ENERGY INTERVALS"
	// PARAMETER $MXLEKE=100; "ELECTRON ENERGY INTERVALS BELOW EKELIM"
	// PARAMETER $MXCMFP=100; "CUMULATIVE ELECTRON MEAN FREE PATH"
	// PARAMETER $MXRANGE=100; "ELECTRON RANGE"
	public static int $MXRL = 100; // "RAYLEIGH SCATTERING SAMPLING INTERVAL"
	public static int $MXBLC = 20; // "MOLIERE'S LOWER CASE B"
	// PARAMETER $MXRNTH=20; "RANDOM NUMBER FOR SUBMEDIAN ANGLES"
	// PARAMETER $MXRNTHI=20; "RANDOM NUMBER FOR SUPERMEDIAN ANGLES"
	public static double $F029 = 1.18E-38; // "LAST F0-VALUE"
	public static double BMIN = 4.5;// see MIX routine
	public static int MSTEPS = $MSSTEPS;
	public static int JRMAX = $MXJREFF;
	public static double[] FSTEP = { 1., 2., 3., 4., 6., 8., 10., 15., 20.,
			30., 40., 60., 80., 100., 150., 200. };
	public static int[] MSMAP = new int[$MXJREFF];
	public static double[] FSQR = new double[$MSSTEPS];
	public static int MXV1 = 0;
	public static int MXV2 = 0;
	public static double[] VERT1 = new double[$MXVRT1];
	public static double[][] VERT2 = new double[$MXVRT2][$MSSTEPS];
	public static double[] BLCA = new double[$MSSTEPS];
	public static double[] BA = new double[$MSSTEPS];

	// -----------------END DATA FOR COMMON BLOCK
	// MIMSD"-------------------------------
	// -----------------------"DATA FOR EPSTAR"------------------------------------
	public static int EPSTFL = 1;// EPSTFL default=read density corr from
									// file!!!
	public static int IEPST = 1;
	public static int IAPRIM = 1;// 1=read aprime from file
	public static int IAPRFL = 0;// "IAPRFL is a flag to say if APRIM file has
									// been initialized
	// -----------------------END DATA FOR
	// EPSTAR"------------------------------------
	// -------------------"****STERNHEIMER-SELTZER-BERGER (SSB) DATA                         "
	// "   STDATA IS THE ACTUAL STERNHEIMER DATA                 "
	// " AFACT,SK,X0,X1,IEV,CBAR
	// "     MEDTLB CONTAIN THE IDENTIFIERS FOR THE MEDIUM"
	// "DATA FOR COMMON BLOCK LSPION                                      "
	public static double AFACT = 0.0;
	public static double SK = 0.0;
	public static double X0 = 0.0;
	public static double X1 = 0.0;
	public static double CBAR = 0.0;
	public static double IEV = 0.0;
	public static double ISSB = 0;// no user supply data!!!
	// "DATA FOR COMMON BLOCK SPCOMM "
	public static int $MXSTC = 73; // "NUMBER OF MEDIA WITH STERNHEIMER COEFFICIENTS"
	public static int LMED = 1;// default 1 medium
	public static String[] IDSTRN = new String[LMED];// Sternheimer ID
	public static int NMED = $MXSTC;// number of medium!!
	public static String BLKW = "";// empty string!
	// v. tab 2.13.2 SLAC 265-Default Sternheimer Density Effect Coefficients in
	// PEGS4
	// public static int MREC123=20;
	// public static int MREC4=13;
	public static String[] MEDTB = {// 73 rec
	"H2-GAS", "H2-LIQUID", "HE-GAS", "LI", "BE", "C-2.265 G/CM**3",
			"C-1.70 G/CM**3", "N2-GAS'", "O2-GAS", "NE-GAS", "NA", "MG", "AL",
			"SI", "AR-GAS", "K", "CA", "TI", "V", "MN", "FE", "CO", "NI", "CU",
			"ZN", "GE", "SE", "KR-GAS", "RB", "MO", "AG", "CD", "IN", "SN",
			"XE-GAS", "CS", "GD", "TA", "W", "PT", "AU", "HG", "PB", "RN-GAS",
			"U", "AIR-GAS", "CO2-GAS", "POLYETHYLENE", "POLYPROPYLENE",
			"XYLENE", "TOLUENE", "NYLON", "VINYLTOLUENE", "A150-PLASTIC",
			"STILBENE", "POLYSTYRENE", "ANTHRACENE", "LEXAN", "LUCITE", "H2O",
			"MYLAR", "KAPTON", "LIF", "POLYVINYL-CL", "PYREX-GLASS", "SIO2",
			"CAF2", "PHOTOEMULSION", "AGCL", "NAI", "LII", "AGBR", "CSI" };
	// "   STDATA IS THE ACTUAL STERNHEIMER DATA                 "
	// " AFACT, SK, X0, X1, IEV, CBAR
	public static double[][] STDATA = {// //73 rec and 6 columns
	{ 0.03535, 6.790, 1.864, 3.5, 19.2, 9.584 },
			{ 0.09179, 5.831, 0.476, 2.0, 21.8, 3.263 },
			{ 0.0114, 7.625, 2.202, 4.0, 41.8, 11.139 },
			{ 0.3492, 3.233, 0.0966, 2.0, 40.0, 3.122 },
			{ 0.3518, 3.034, -0.0089, 2.0, 63.7, 2.785 },
			{ 0.5848, 2.360, -0.0089, 2.0, 78.0, 2.868 },
			{ 0.7154, 2.191, -0.0089, 2.0, 78.0, 3.155 },
			{ 0.2120, 3.041, 1.738, 4.0, 82.0, 10.540 },
			{ 0.2666, 2.825, 1.754, 4.0, 95.0, 10.700 },
			{ 0.1202, 3.357, 2.073, 4.5, 137.0, 11.904 },
			{ 0.2204, 3.103, 0.4515, 2.8, 149.0, 5.053 },
			{ 0.1714, 3.223, 0.2386, 2.8, 156.0, 4.530 },
			{ 0.3346, 2.795, 0.0966, 2.5, 166.0, 4.239 },
			{ 0.3755, 2.720, 0.0966, 2.5, 173.0, 4.435 },
			{ 0.1902, 2.982, 1.764, 4.5, 188.0, 11.948 },
			{ 0.3041, 2.674, 0.2386, 3.0, 190.0, 5.642 },
			{ 0.2177, 2.874, 0.1751, 3.0, 191.0, 5.040 },
			{ 0.1782, 2.946, 0.0485, 3.0, 233.0, 4.445 },
			{ 0.1737, 2.935, -0.0089, 3.0, 245.0, 4.266 },
			{ 0.1996, 2.812, -0.0089, 3.0, 272.0, 4.270 },
			{ 0.2101, 2.771, -0.0089, 3.0, 286.0, 4.291 },
			{ 0.2229, 2.713, -0.0089, 3.0, 297.0, 4.260 },
			{ 0.2504, 2.619, -0.0089, 3.0, 311.0, 4.312 },
			{ 0.2557, 2.613, -0.0089, 3.0, 322.0, 4.419 },
			{ 0.3163, 2.468, 0.0485, 3.0, 330.0, 4.691 },
			{ 0.2809, 2.647, 0.2386, 3.0, 350.0, 5.141 },
			{ 0.2979, 2.635, 0.2386, 3.0, 348.0, 5.321 },
			{ 0.1519, 3.030, 1.716, 4.8, 352.0, 12.512 },
			{ 0.1450, 3.078, 0.4515, 3.5, 363.0, 6.478 },
			{ 0.2228, 2.824, 0.1751, 3.0, 424.0, 4.879 },
			{ 0.3091, 2.563, -0.0089, 3.0, 470.0, 5.063 },
			{ 0.1853, 2.819, 0.0485, 3.3, 469.0, 5.273 },
			{ 0.2004, 2.790, 0.1751, 3.3, 487.0, 5.517 },
			{ 0.1898, 2.839, 0.2386, 3.3, 488.0, 5.534 },
			{ 0.1329, 3.020, 1.563, 5.0, 482.0, 12.728 },
			{ 0.2214, 2.784, 0.4515, 3.5, 488.0, 6.914 },
			{ 0.2068, 2.686, 0.0485, 3.5, 591.0, 5.874 },
			{ 0.1663, 2.805, 0.1751, 3.5, 718.0, 5.526 },
			{ 0.1499, 2.870, 0.1751, 3.5, 727.0, 5.406 },
			{ 0.1465, 2.903, 0.0966, 3.5, 790.0, 5.473 },
			{ 0.1533, 2.881, 0.0966, 3.5, 790.0, 5.575 },
			{ 0.1824, 2.798, 0.2386, 3.5, 800.0, 5.961 },
			{ 0.1861, 2.814, 0.2386, 3.5, 823.0, 6.202 },
			{ 0.1130, 3.023, 1.537, 5.3, 794.0, 13.284 },
			{ 0.1362, 3.034, 0.2386, 3.5, 890.0, 5.869 },
			{ 0.2466, 2.879, 1.742, 4.0, 85.7, 10.595 },
			{ 0.1999, 3.022, 1.648, 4.0, 88.7, 10.239 },
			{ 0.4875, 2.544, 0.1379, 2.0, 57.4, 3.002 },
			{ 0.2493, 2.975, 0.1537, 2.3, 59.2, 3.126 },
			{ 0.2755, 2.911, 0.1695, 2.3, 61.8, 3.270 },
			{ 0.2830, 2.890, 0.1722, 2.3, 62.5, 3.303 },
			{ 0.5345, 2.439, 0.1336, 2.0, 63.9, 3.063 },
			{ 0.3495, 2.749, 0.1467, 2.2, 64.7, 3.201 },
			{ 0.5462, 2.435, 0.1329, 2.0, 65.1, 3.110 },
			{ 0.2989, 2.851, 0.1731, 2.3, 67.7, 3.367 },
			{ 0.3670, 2.724, 0.1647, 2.2, 68.7, 3.300 },
			{ 0.5858, 2.364, 0.1146, 2.0, 69.5, 3.151 },
			{ 0.3865, 2.664, 0.1608, 2.2, 73.1, 3.321 },
			{ 0.3996, 2.606, 0.1824, 2.2, 74.0, 3.330 },
			{ 0.2065, 3.007, 0.2400, 2.5, 75.0, 3.502 },
			{ 0.3124, 2.782, 0.1561, 2.3, 78.7, 3.326 },
			{ 0.4061, 2.614, 0.1492, 2.2, 79.3, 3.342 },
			{ 0.1308, 3.476, 0.0171, 2.5, 94.0, 3.167 },
			{ 0.1873, 2.962, 0.1558, 2.8, 108.2, 4.053 },
			{ 0.2988, 2.805, 0.1479, 2.5, 134.0, 3.971 },
			{ 0.1440, 3.220, 0.1385, 2.8, 139.2, 4.003 },
			{ 0.3750, 2.592, 0.0676, 2.5, 166.0, 4.065 },
			{ 0.3416, 2.496, 0.1009, 3.0, 331.0, 5.332 },
			{ 0.1243, 3.002, -0.0138, 3.5, 398.4, 5.344 },
			{ 0.1560, 2.926, 0.1203, 3.5, 452.0, 6.057 },
			{ 0.1785, 2.845, 0.0892, 3.5, 485.1, 6.267 },
			{ 0.1351, 2.976, 0.0358, 3.5, 487.2, 5.616 },
			{ 0.1796, 2.840, 0.0395, 3.5, 553.1, 6.281 } };

	// -------------------END****STERNHEIMER-SELTZER-BERGER (SSB) DATA "
	// --------------------SPINIT---------------------------
	// Initializes stopping power functions for a particular medium.
	public static double TOLN10 = 2.0 * Math.log(10.0);
	public static double SPC1 = 0.0;
	public static double SPC2 = 0.0;
	// "NELEPS IS THE NUMBER OF ELEMENTS IN THE MATERIAL             "
	// "ZEPST(I) (INTEGER), WEPST(I) ARE THE Z VALUE AND FRACTION    "
	// "   BY WEIGHT OF THE I-TH ELEMENT IN THE DATA FILE.           "
	// "IAPRIM is a flag to tell which correction to the brem        "
	// "       cross section to use:      0 =>  old Koch and Motz    "
	// " =1(defaulat) read in new data file; =2, use no corrections. "
	// "IAPRFL is a flag to say if APRIM file has been initialized   "
	// READ(20,:A:)EPSTTL;:A: FORMAT(A);//is title e.g. sodium iodide from
	// .density file
	public static int NEPST = 0;// is npairdelta
	public static int NELEPS = 0;// is NR. of z elements=NE
	public static double[] ZEPST;// is Z[i]
	public static double[] WEPST;// is RHOZ[i]
	public static double[] EPSTEN;// is energy
	public static double[] EPSTD;// is delta
	public static double AP = 0.;// low energy threshold for soft bremsstrahlung
									// production
	public static double AE = 0.;// low energy threshold for knock-on electron
									// production
	public static double UP = 0.;// idem but upper--->restricted by available
									// data !!
	public static double UE = 0.;// idem but upper--->restricted by available
									// data !!
	// -----------------------------------------------------
	// -------------DIFFER brems and pair constants--------------------------
	public static double[] ALPHI = new double[2];
	public static double AL2 = 0.;// ln(2)-------in Differ init!!!
	public static double[] ALFP1 = new double[2];
	public static double[] ALFP2 = new double[2];
	public static double[] BPAR = new double[2];
	public static double DELCM = 0.;
	public static double[] DELPOS = new double[2];
	public static double[] DL1 = new double[6];
	public static double[] DL2 = new double[6];
	public static double[] DL3 = new double[6];
	public static double[] DL4 = new double[6];
	public static double[] DL5 = new double[6];
	public static double[] DL6 = new double[6];
	// -----------------------------------------------------------------------
	// -------------------TRESHOLD-------------------------------------------
	public static double TE = 0.0;
	public static double TET2 = 0.0;
	public static double TEM = 0.0;
	public static double THBREM = 0.0;
	public static double THMOLL = 0.0;
	// -----------------------------------------------------------------------
	// ---------------------------Fit electron, fotoel,pair and form factor
	// data-------------
	public static double EBINDA = 0.0;
	// photelectric and pair data
	public static double[] EKEDGE = new double[100];// each elements!!!"
	// EKEDGE IS THE K IONIZATION ENERGY IN KEV "
	public static int[] NPHE = new int[100];// number of data available for each
											// element
	public static int $mxpaire = 17;
	public static int $mxphote = 61;
	public static double[][] PHE = new double[$mxphote][100];// energies for
																// each element
	public static double[][] PHD = new double[$mxphote][100];
	// "  PHD(IE,Z) IS THE PHOTO CROSS-SECTION AT K=EPP(IE,Z)  "
	public static double[] PRE = new double[$mxpaire];// energies common for
														// each element
	public static double[][] PRD = new double[$mxpaire][100];
	// "  PRD(IE,Z) IS THE PAIR PROC CROSS-SECTION AT K=EPP(IE,Z)     "
	// "  THE CROSS-SECTIONS ARE IN UNITS OF BARNS/ATOM.    AND          "
	// "   NEPP(Z) IS THE # OF ENERGIES AT WHICH DATA IS GIVEN              "
	// "   EPP(IE,Z) IS THE IE'TH ENERGY FOR ELEMENT Z                      "
	public static double[][] COHE = new double[$mxphote][100];
	// ---------------------------------------------------------FormFactor------
	public static double[] XVAL = new double[100];
	public static double[][] AFAC = new double[100][100];
	// ------------------------------------------------------------------
	// --------------------------RESULTS PFWL-----------------------------
	public static int NEL = 0;
	public static double AXE = 0.;
	public static double BXE = 0.;
	public static double[][] AFE = new double[$MXEKE][8];
	public static double[][] BFE = new double[$MXEKE][8];
	public static int NALE = $MXEKE;// 150
	public static double EPE = 0.01;
	public static double[] ZTHRE = new double[8];
	public static double[] ZEPE = new double[8];
	public static int NIPE = 20;
	public static double DELC = 0.;
	public static double CONST = 0.;
	public static double XLNZ = 0.;
	public static double DELTAM = 0.;
	// /aprime.data
	public static int NAPRZ = 0;
	public static int NAPRE = 0;
	public static double[] EPRIM;// (IE),IE=1,NAPRE
	public static double[] ZPRIM;
	public static double[][] APRIMD;
	// ----------------
	public static double E = 0.;
	public static double K = 0.;
	public static String funcname = "";
	// "     DATA FROM BETHE'S TABLE - MOLIERE Functions!! -  "
	public static double[] TH = { .05, .2, .4, .6, .8, 1., 1.2, 1.4, 1.6, 1.8,
			2., 2.2, 2.4, 2.6, 2.8, 3., 3.2, 3.4, 3.6, 3.8, 4.07, 4.5, 5., 5.5,
			6.13, 7., 8., 9., 9.75 };
	public static double[] DTH = { .1, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2,
			0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.35, 0.5,
			0.5, 0.5, 0.75, 1.0, 1.0, 1.0, 0.5 };
	public static double[] F0 = { 2., 1.9216, 1.7214, 1.4094, 1.0546, .7338,
			.4738, .2817, .1546, .0783, .0366, .01581, .0063, .00232, 7.9E-4,
			2.5E-4, 7.3E-5, 1.9E-5, 4.7E-6, 1.1E-6, 2.3E-7, 3.E-9, 2.E-11,
			2.E-13, 5.E-16, 1.E-21, 3.E-28, 1.E-35, $F029 };
	public static double[] F1 = { .8456, .7038, .3437, -0.0777, -0.3981,
			-0.5285, -0.4770, -.3183, -.1396, -.0006, +0.0782, .1054, .1008,
			.08262, .06247, .0455, .03288, .02402, .01791, .01366, .010638,
			.00614, .003831, .002527, .001739, .000908, .0005211, .0003208,
			.0002084 };
	public static double[] F2 = { 2.4929, 2.0694, 1.0488, -.0044, -.6068,
			-.6359, -.3086, .0525, .2423, .2386, .1316, .0196, -.0467, -.0649,
			-.0546, -.03568, -.01923, -.00847, -.00264, 5.E-5, .0010741,
			.0012294, .0008326, .0005368, .0003495, .0001584, 7.83E-5, 4.17E-5,
			2.37E-5 };
	public static int $NAPRE = 115;// " Maximum number of energies ( > 18 )      "
	public static int $NAPRZ = 14;// " Maximum number of elements ( > 5 )       "
	public static double[][] APRIMD1 = new double[$NAPRE][$NAPRZ];
	public static double[] EPRIM1 = new double[$NAPRE];
	public static double[] ZPRIM1 = new double[$NAPRZ];
	public static double[] APRIMZ = new double[$NAPRE];
	public static boolean first_time = true;
	public static double SAFETY = 0.8;
	public static double TABSMX = 10.0;
	// ---------DCADRE--------------------------------------------------
	private static double DIFFq = 0.;
	private static boolean H2CONVq = false;
	private static boolean AITKENq = false;
	private static int LM1q = 0;
	private static int N2q = 0;
	private static double FNq = 0.;
	private static int ISTEPq = 0;
	private static int IIq = 0;
	private static int IIIq = 0;
	private static double HOVNq = 0.;
	private static double FIq = 0.;
	private static int ISTEP2q = 0;
	private static double SUMq = 0.;
	private static double SUMABSq = 0.;
	@SuppressWarnings("unused")
	private static double ABSIq = 0.;
	private static int ITq = 0;
	private static double TABTLMq = 0.;
	private static double ERGLq = 0.;
	private static double ERGOALq = 0.;
	private static double FEXTRPq = 0.;
	private static double ERRERq = 0.;
	private static boolean RIGHTq = false;
	private static double STEPq = 0.;
	private static double ASTEPq = 0.;
	private static double TABSq = 0.;
	private static int Lq = 0;
	private static int Nq = 0;
	private static double[] TSq = new double[2049];
	private static double[][] Tq = new double[10][10];
	private static double ZEROq = 0.0;
	private static double P1q = 0.1;
	private static double HALFq = 0.5;
	private static double ONEq = 1.0;
	private static double TWOq = 2.0;
	private static double FOURq = 4.0;
	private static double FOURP5q = 4.5;
	private static double TENq = 10.0;
	private static double HUNq = 100.0;
	private static double AITLOWq = 1.1;
	private static double H2TOLq = 0.15;
	private static double AITTOLq = 0.1;
	private static double JUMPTLq = 0.01;
	private static int MAXTSq = 2049;
	private static int MAXTBLq = 10;
	private static int MXSTGEq = 30;
	private static double SLOPEq = 0.;
	private static double FBEG2q = 0.;
	private static double[] RNq = new double[4];
	private static boolean[] REGLSVq = new boolean[30];
	private static double[] BEGINq = new double[30];
	private static double[] FINISq = new double[30];
	private static double[] ESTq = new double[30];
	private static int[] IBEGSq = new int[30];
	private static int NNLEFTq = 0;
	private static double ALG4O2q = 0.;
	private static double CADREq = 0.;
	private static double CURESTq = 0.;
	private static double VINTq = 0.;
	private static double DCADREq = 0.;
	private static double LENGTHq = 0.;
	private static double ERRRq = 0.;
	private static double ERRAq = 0.;
	private static double STEPMNq = 0.;
	private static double STEPNMq = 0.;
	private static double STAGEq = 0.;
	private static int ISTAGEq = 0;
	private static double FNSIZEq = 0.;
	private static double PREVERq = 0.;
	private static boolean REGLARq = false;
	private static double BEGq = 0.;
	private static double RVALq = 0.;
	private static double FBEGq = 0.;// @@
	private static int IBEGq = 0;
	private static double ENDq = 0.;
	private static double FENDq = 0.;
	private static int IENDq = 0;
	private static double Aq = 0.;
	private static double Bq = 0.;
	private static double AERRq = 0.;
	private static double RERRq = 0.;
	private static double ERRORq = 0.;
	private static int IERq = 0;
	private static int Iq = 0;
	private static double SINGq = 0.;
	private static double FEXTM1q = 0.;
	private static double[] Rq = new double[10];
	private static double[] AITq = new double[10];
	private static double[] DIFq = new double[10];
	@SuppressWarnings("unused")
	private static double ALPHAq = 0.;
	private static double H2NXTq = 0.;
	private static double SINGNXq = 0.;
	private static double ERRETq = 0.;
	private static double H2TFEXq = 0.;
	// --------------------------------------------
	public static double[] AFAC2 = new double[100];
	public static double[] AFFI = new double[100];
	public static double CMOLL = 0.;
	public static double C1 = 0.;
	public static double C2 = 0.;
	public static double T0 = 0.;
	public static double A = 0.;
	public static double CBHAB = 0.;
	public static double B1 = 0.;
	public static double B2 = 0.;
	public static double B3 = 0.;
	public static double B4 = 0.;
	public static double BETASI = 0.;
	public static int LA = 0;
	public static int LB = 0;
	public static int LD = 0;
	public static double CCOMP = 0.;
	public static double C3 = 0.;
	public static double K0 = 0.;
	public static int LC = 0;
	public static int LE = 0;
	private static boolean exitq = false;

	// ------------------------------------------------------------------
	/**
	 * Set the medium name
	 * @param s s
	 */
	public static void setMedium(String s) {
		medium = s;
	}

	/**
	 * Set IAPRIM flag.<p>
	 * IAPRIM is a flag to tell which correction to the brem cross section to use: <p>
	 * 0 = old Koch and Motz; 1(default) = read in new data file; 2 = use no corrections.
	 * @param n n
	 */
	public static void setIAPRIM(int n) {
		// "   IAPRIM = 0  equivalent to old PEGS4 (default)                      "
		// "            1  reads in values from unit 22                           "
		// " 2 sets APRIM to 1.0
		IAPRIM = n;
	}

	/**
	 * Set IUNRIST flag.<p>
	 * IUNRST=0 (default) MEANS RESTRICTED STOPPING POWER<p>
	 * IUNRST=1 MEANS UNRESTRICTED COLLISION STOPPING POWER.<p>
	 * IUNRST=2 MEANS CSDA DATA - TOTAL UNRESTRICTED STOPPING POWER AND NO DISCRETE INTERACTIONS<p>
	 * IUNRST=3 MEANS ALLOW BREM EVENTS BUT NO MOLLER INTERACTIONS<p>
	 * IUNRST=4 MEANS ALLOW MOLLER EVENTS BUT NO BREM - NOTE THIS GIVES COMPLETE GARBAGE IN A RUN SINCE ALL BREM IS DUMPED ON SPOT<p>
	 * IUNRST=5 MEANS UNRESTRICTED RADIATIVE STOPPING POWER
	 * @param n n
	 */
	public static void setIUNRST(int n)// default=0
	{
		// IUNRST=0 MEANS RESTRICTED STOPPING POWER"
		// "         IUNRST=1 MEANS UNRESTRICTED COLLISION STOPPING POWER"
		// "         IUNRST=2 MEANS CSDA DATA - TOTAL UNRESTRICTED STOPPING POWER"
		// "                  AND NO DISCRETE INTERACTIONS"
		// "         IUNRST=3 MEANS ALLOW BREM EVENTS BUT NO MOLLER INTERACTIONS"
		// "         IUNRST=4 MEANS ALLOW MOLLER EVENTS BUT NO BREM - NOTE THIS"
		// "                  GIVES COMPLETE GARBAGE IN A RUN SINCE ALL BREM IS "
		// "                  DUMPED ON SPOT"
		// "         IUNRST=5 MEANS UNRESTRICTED RADIATIVE STOPPING POWER"
		IUNRST = n;
	}

	/**
	 * Set default values for AE, AP, UE and UP. Must be called before SPINIT and after readDensity...
	 * AP  = low energy threshold for soft bremsstrahlung production; <p>
	 * AE =  low energy threshold for knock-on electron production <p>
	 * UE, UP = same as above but refer to upper limit!
	 */
	public static void defaultAEAPUEUP()// must be called before SPINIT and
										// after readDensity...
	{
		AE = EPSTEN[0] + RM;
		AP = EPSTEN[0];
		UE = EPSTEN[NEPST - 1] + RM;
		UP = EPSTEN[NEPST - 1];
	}

	/**
	 * Set UE, the upper energy threshold for knock-on electron production.
	 * @param n n
	 */
	public static void setUE(double n) {
		UE = n;
	}

	/**
	 * Set UP, the upper energy threshold for soft bremsstrahlung production;.
	 * @param n n
	 */
	public static void setUP(double n) {
		UP = n;
	}

	/**
	 * Set AE, the low energy threshold for knock-on electron production.
	 * @param n n
	 */
	public static void setAE(double n) {
		AE = n;
	}

	/**
	 * Set AP, the low energy threshold for soft bremsstrahlung production;.
	 * @param n n
	 */
	public static void setAP(double n) {
		AP = n;
	}

	// //Defines state of mixture: zero (default) for solid or liquid, otherwise
	// value gives gas pressure (atm).
	/**
	 * Defines state of mixture: zero (default) for solid or liquid, otherwise value gives gas pressure (atm).
	 * @param n n
	 */
	public static void setGASP(double n) {
		GASP = n;
	}

	/**
	 * EPSTFL = density correction flag. If 0, then it is calculated, if 1, then it is read from file. 
	 * Default is 1(read density correction from file).
	 * @param n n
	 */
	public static void setEPSTFL(int n) {
		EPSTFL = n;// default is 1(read density correction from file)
		// only 0 and 1 is available
	}

	/**
	 * Set Sternheimer ID.
	 * @param s s
	 */
	public static void setIDSTRN(String[] s)
	{
		LMED = s.length;
		IDSTRN = new String[LMED];// Sternheimer ID
		for (int I = 0; I < LMED; I++) {
			IDSTRN[I] = s[I];
		}
	}

	/**
	 * IMIXT flag is 0 for element and compound and 1 for mixture.
	 * @param n n
	 */
	public static void setIMIXT(int n)
	{
		IMIXT = n;
	}

	/**
	 * Converts ASCII int value to a String.
	 * 
	 * @param i
	 *            the ASCII integer
	 * @return the string representation
	 */
	public static String asciiToStr(int i) {
		char a[] = new char[1];
		a[0] = (char) i;
		return (new String(a)); // caracterul sub forma de sir caracterizat de
								// codul ASCII
	}

	/**
	 * Return last blank character index from a string s
	 * @param s s
	 * @return the result
	 */
	public static int lsblnk(String s)// last blank character
	{
		int result = 0;
		for (int i = s.length(); i > 0; i--) {
			char chr = s.charAt(i - 1);
			if (Character.isWhitespace(chr)) {
				result = i - 1;
				return result;
			}
		}
		return result;
	}

	/**
	 * Internally used by readDensity routine.
	 * @param s s
	 */
	public static void solveLine2(String s) {
		String result = "";
		Vector<String> v = new Vector<String>();
		String ss = s;// .trim();//remove first and last whitespace characters!
		StringBuffer desc = new StringBuffer();
		for (int i = 0; i < ss.length(); i++) {
			char chr = ss.charAt(i);

			if (!Character.isWhitespace(chr)) {
				desc.append(chr);
				if (i == ss.length() - 1)// last char!!
				{
					result = desc.toString();
					v.addElement(result);
				}
			} else {
				if ((desc.toString()).compareTo("") != 0)// not null
				{
					result = desc.toString();
					v.addElement(result);
				}
				desc = new StringBuffer();
			}
		}
		// for (int i=0;i<10;i++)
		// v.addElement(new Integer(i));//OK!!!!!!
		// ---------------do something with ....
		String npairdeltas = (String) v.elementAt(0);
		// npairdelta=Convertor.stringToInt(npairdeltas);//System.out.println(npairdelta);
		npairdelta = stringToInt(npairdeltas);// System.out.println(npairdelta);
		NEPST = npairdelta;
		edelta = new double[npairdelta];
		delta = new double[npairdelta];
		EPSTEN = new double[NEPST];// is energy
		EPSTD = new double[NEPST];// is delta

		String ievs = (String) v.elementAt(1);
		// double iev=Convertor.stringToDouble(ievs);//System.out.println(iev);
		double iev = stringToDouble(ievs);// System.out.println(iev);
		IEV = iev;// solved from file
		// ---------------------------------
		String rhos = (String) v.elementAt(2);
		// double rhoe=Convertor.stringToDouble(rhos);
		double rhoe = stringToDouble(rhos);
		setRHO(rhoe);
		String nes = (String) v.elementAt(3);
		// int nee=Convertor.stringToInt(nes);
		int nee = stringToInt(nes);
		setNE(nee);
		NELEPS = nee;// solve nee!!
		ZEPST = new double[NELEPS];// is Z[i]
		WEPST = new double[NELEPS];// is RHOZ[i]
		ASYM = new String[nee];
	}

	/**
	 * Internally used by readDensity routine.
	 * @param s s
	 */
	public static void solveZDependencies(String s) {
		String result = "";// System.out.println(zdepv.toString());
		StringBuffer desc = new StringBuffer();
		int ctrl = 0;
		int append = 0;
		int vsize = 0;
		for (int i = 0; i < s.length(); i++) {
			char chr = s.charAt(i);
			if (!Character.isWhitespace(chr)) {
				desc.append(chr);
				if (i == s.length() - 1)// last char!!
				{
					result = desc.toString();
					zdepv.addElement(result);
					append++;
				}
			} else {
				if ((desc.toString()).compareTo("") != 0)// not null
				{
					result = desc.toString();
					zdepv.addElement(result);
					append++;
				}
				desc = new StringBuffer();
			}
			// directly solving the Z params!
			if (append != 0) {// System.out.println(ctrl);
				if (ctrl == 0)// first=Z
				{
					vsize = zdepv.size();
					String ss = (String) zdepv.elementAt(vsize - 1);// System.out.println(ss);
					// Z[jzet]=Convertor.stringToDouble(ss);
					Z[jzet] = stringToDouble(ss);
					ZEPST[jzet] = Z[jzet];// is Z[i]
					Double d = new Double(Z[jzet]);
					int I01 = d.intValue() - 1;// 1,H->index 0 in table!!
					if (WA[jzet] == 0.)
						WA[jzet] = WATBL[I01];// "USE OUR TABLE (DEFAULT)"
					ASYM[jzet] = ASYMT[I01];
					jzet++;
					ctrl++;
				} else// RHOZ or PZ NEXT TO Z!!!
				{
					if (NE == 1) {
						MTYP = "ELEM";
						setIMIXT(0);// elements
						vsize = zdepv.size();
						String ss = (String) zdepv.elementAt(vsize - 1);
						// PZ[jrhozet]=Convertor.stringToDouble(ss);
						PZ[jrhozet] = stringToDouble(ss);
						RHOZ[jrhozet] = PZ[jrhozet] * WA[jrhozet];
						WEPST[jrhozet] = RHOZ[jrhozet];
						jrhozet++;
					} else// compounds or mixture, but from data allways mixture
					{
						MTYP = "MIXT";
						setIMIXT(1);
						vsize = zdepv.size();
						String ss = (String) zdepv.elementAt(vsize - 1);
						// RHOZ[jrhozet]=Convertor.stringToDouble(ss);
						RHOZ[jrhozet] = stringToDouble(ss);
						WEPST[jrhozet] = RHOZ[jrhozet];
						PZ[jrhozet] = RHOZ[jrhozet] / WA[jrhozet];
						jrhozet++;
					}
					ctrl--;
				}
				append--;
			}
		}
		// if Z and RHOZ(and also PZ) are complet->clear vector
		if ((Z[NE - 1] != 0.) && (RHOZ[NE - 1] != 0.))
			zdepv = new Vector<String>();// clear
	}

	/**
	 * Internally used by readDensity routine.
	 * @param s s
	 */
	public static void solveDelta(String s) {
		String result = "";
		char comma = ',';
		StringBuffer desc = new StringBuffer();
		int ctrl = 0;
		int append = 0;
		int vsize = 0;
		for (int i = 0; i < s.length(); i++) {
			char chr = s.charAt(i);
			if ((!Character.isWhitespace(chr)) && (chr != comma)) {
				// if (chr==comma)//final=never comma!!
				desc.append(chr);
				if (i == s.length() - 1)// last char!!
				{
					result = desc.toString();
					zdepv.addElement(result);
					append++;
				}
			} else {
				if ((desc.toString()).compareTo("") != 0)// not null
				{
					result = desc.toString();
					zdepv.addElement(result);
					append++;
				}
				desc = new StringBuffer();
			}
			// directly solving the Z params!
			if (append != 0) {
				if (ctrl == 0)// first=edelta
				{
					vsize = zdepv.size();
					String ss = (String) zdepv.elementAt(vsize - 1);
					// edelta[jenergy]=Convertor.stringToDouble(ss);
					edelta[jenergy] = stringToDouble(ss);
					EPSTEN[jenergy] = edelta[jenergy];// is energy
					jenergy++;
					ctrl++;
				} else// delta
				{
					vsize = zdepv.size();
					String ss = (String) zdepv.elementAt(vsize - 1);
					// delta[jdelta]=Convertor.stringToDouble(ss);
					delta[jdelta] = stringToDouble(ss);
					EPSTD[jdelta] = delta[jdelta];// is delta
					jdelta++;
					ctrl--;
				}
				append--;
			}
		}

	}

	/**
	 * Read density file and initialize global variables. Distinction is made between 
	 * mixture and element.
	 * @param filename filename
	 * @param comp true if mixture, false otherwise
	 */
	public static void readDensityFile(String filename, boolean comp) {
		zdepv = new Vector<String>();// clear
		jzet = 0;
		jrhozet = 0;
		jenergy = 0;
		jdelta = 0;

		medium = stripFileExtension(filename);
		String nume = "";
		if (comp)// mixture
			nume = userdirs + file_sep + datas + file_sep + egsData + file_sep
					+ compoundss + file_sep + filename;
		else
			nume = userdirs + file_sep + datas + file_sep + egsData + file_sep
					+ elementss + file_sep + filename;
		densPath = nume;
		int i = 0;
		int lnr = 0;// line number
		@SuppressWarnings("unused")
		int chri = 0;// line chr number
		StringBuffer desc = new StringBuffer();
		String line = "";
		String temps = "";
		char lineSep = System.getProperty("line.separator").charAt(0);
		try {
			FileInputStream in = new FileInputStream(nume);

			while ((i = in.read()) != -1) {
				if ((char) i != lineSep) {
					desc.append((char) i);
					chri++;
				} else {
					in.read();// step over to not print 2x times EOL (1 from
								// println, 1 from text!)
					lnr++;
					line = desc.toString();
					temps = line;
					if (lnr == 2) {
						solveLine2(temps);
					}
					if (lnr > 2) {
						if (Z[NE - 1] == 0.)// Z dependent data
						{
							solveZDependencies(temps);
						} else// delta data
						{
							solveDelta(temps);
						}
					}
					desc = new StringBuffer();
					// System.out.println(line+" chr= "+chri);//OK
					// System.out.print(line+" chr= "+chri);//YES!
					chri = 0;
				}
			}
			in.close();
		} catch (Exception e) {

		}
		// System.out.println(zdepv.toString());
	}

	/**
	 * Read density file and initialize global variables and we already know if it is 
	 * an element or a mixture and provide its file path. 
	 * @param filenamepath filenamepath
	 */
	public static void readDensityFile(String filenamepath) {
		zdepv = new Vector<String>();// clear
		jzet = 0;
		jrhozet = 0;
		jenergy = 0;
		jdelta = 0;

		medium = stripFileExtension(filenamepath);
		medium = getFileName(medium);// System.out.println(medium);
		String nume = filenamepath;
		densPath = nume;
		int i = 0;
		int lnr = 0;// line number
		@SuppressWarnings("unused")
		int chri = 0;// line chr number
		StringBuffer desc = new StringBuffer();
		String line = "";
		String temps = "";
		char lineSep = System.getProperty("line.separator").charAt(0);
		try {
			FileInputStream in = new FileInputStream(nume);

			while ((i = in.read()) != -1) {
				if ((char) i != lineSep) {
					desc.append((char) i);
					chri++;
				} else {
					in.read();// step over to not print 2x times EOL (1 from
								// println, 1 from text!)
					lnr++;
					line = desc.toString();
					temps = line;
					if (lnr == 2) {
						solveLine2(temps);
					}
					if (lnr > 2) {
						if (Z[NE - 1] == 0.)// Z dependent data
						{
							solveZDependencies(temps);
						} else// delta data
						{
							solveDelta(temps);
						}
					}
					desc = new StringBuffer();
					// System.out.println(line+" chr= "+chri);//OK
					// System.out.print(line+" chr= "+chri);//YES!
					chri = 0;
				}
			}
			in.close();
		} catch (Exception e) {

		}
		// System.out.println(zdepv.toString());
	}

	/**
	 * Read A' factor (correction to the brem cross section) from a file.
	 * @param filename filename
	 */
	public static void readAprimeData(String filename) {
		String nume = "";
		char comma = ',';
		nume = userdirs + file_sep + datas + file_sep + egsData + file_sep
				+ filename;
		boolean haveData = false;
		int i = 0;
		int lnr = 0;// data number
		int ne = 0;// contor for EPRIM
		int nz = 0;// contor for APRIMD[][nz]
		int na = 0;// contor for APRIMD[na][nz]
		int nzz = 0;// contor for ZPRIM

		StringBuffer desc = new StringBuffer();
		try {
			FileInputStream in = new FileInputStream(nume);

			while ((i = in.read()) != -1) {
				if (!Character.isWhitespace((char) i) && ((char) i != comma)) {
					desc.append((char) i);
					haveData = true;
				} else {
					if (haveData)// we have data
					{
						haveData = false;// reset
						lnr++;
						if (lnr == 1)// NAPRZ
						{
							String NAPRZs = desc.toString();
							// NAPRZ=Convertor.stringToInt(NAPRZs);//System.out.println(NAPRZ);
							NAPRZ = stringToInt(NAPRZs);// System.out.println(NAPRZ);
							ZPRIM = new double[NAPRZ];
						}
						if (lnr == 2)// NAPRE
						{
							String NAPREs = desc.toString();
							// NAPRE=Convertor.stringToInt(NAPREs);//System.out.println(NAPRE);
							NAPRE = stringToInt(NAPREs);// System.out.println(NAPRE);
							EPRIM = new double[NAPRE];
							APRIMD = new double[NAPRE][NAPRZ];
						}

						if (lnr > 2)// else
						{
							if (ne < NAPRE)// read EPRIM
							{
								String s = desc.toString();
								// EPRIM[ne]=Convertor.stringToDouble(s);//System.out.println(EPRIM[ne]);
								EPRIM[ne] = stringToDouble(s);// System.out.println("eprim  "
																// +EPRIM[ne]+"   "+ne);
								ne++;
							} else if ((na < NAPRE) && (nz < NAPRZ)) {
								String s = desc.toString();// System.out.println(desc.toString());
								if (nzz == nz)// read for ZPRIM
								{
									// ZPRIM[nzz]=Convertor.stringToDouble(s);//System.out.println(ZPRIM[nzz]);
									ZPRIM[nzz] = stringToDouble(s);// System.out.println("zprim "+
																	// ZPRIM[nzz]+"    "+nzz);
									nzz++;
								} else// read for APRIMD
								{
									// APRIMD[na][nz]=Convertor.stringToDouble(s);//System.out.println(APRIMD[na][nz]);
									APRIMD[na][nz] = stringToDouble(s);// System.out.println(APRIMD[na][nz]+" na "+na+" nz "+nz);
									na++;
									if ((na == NAPRE) && (nz != NAPRZ)) {
										nz++;// System.out.println(nzd);
										na = 0;
									}
								}
							}
						}
					}
					desc = new StringBuffer();
				}
			}
			in.close();
		} catch (Exception e) {

		}

	}

	// if mixture only rhoz matter, else only pz!!!!
	/**
	 * Set NE related variables from external sources.
	 * NE = NUMBER OF DIFFERENT TYPES OF ATOMS IN THE MATERIAL
	 * PZ = PROPORTION OF ELEMENT OF TYPE I. IF A COUMPOUND, 
	 * THEN PZ(I) WILL BE THE NUMBER OF ATOMS OF TYPE I IN THE MOLECULE. 
	 * IF A MIXTURE,SUCH AS CONCRETE, PZ(I) COULD BE THE PER CENT OF THE ATOMS WHICH ARE OF TYPE I. 
	 * Z = PERIODIC NUMBER OF ATOMS OF TYPE I; 
	 * WA = ATOMIC WEIGHT FOR ATOMS OF TYPE I; 
	 * RHOZ = PARTIAL DENSITY DUE TO ATOMS OF TYPE I (GM/CM**3).  
	 * If mixture only rhoz matter, else only pz!
	 * @param ASYMe ASYMe, default it is link to ASYMT table
	 * @param pze pze
	 * @param rhoze rhoze
	 */
	public static void setNEVAR2// (double[] pze,double[] ze, double[] wae,
								// double [] rhoze)
	(String[] ASYMe, double[] pze, double[] rhoze) {
		if (NE == 1)// ELEM->IMIXT=0
		{
			MTYP = "ELEM";
			IMIXT = 0;
			PZ[NE - 1] = 1.;// allways override!! PZ[0] for NE=1 must be 1!!
		} else {
			MTYP = "MIXT";
			IMIXT = 1;
		}
		ASYM = new String[NE];
		setEPSTFL(0);// not read from file
		for (int I = 0; I < NE; I++) {
			if (Z[I] == 0.)
				Z[I] = ZTBL(ASYMe[I]);// default use our table
			Double d = new Double(Z[I]);
			int I01 = d.intValue() - 1;// 1,H->index 0 in table!!
			if (WA[I] == 0.)
				WA[I] = WATBL[I01];// "USE OUR TABLE (DEFAULT)"
			ASYM[I] = ASYMT[I01];
			if (IMIXT != 0) {
				RHOZ[I] = rhoze[I];
				if (NE != 1)// mixture,never!!!
					PZ[I] = RHOZ[I] / WA[I];
			} else {
				if (NE != 1)
					PZ[I] = pze[I];// e from external
				RHOZ[I] = PZ[I] * WA[I];
			}
		}
		setIDSTRN(ASYM);
		if ((NE == 1) && (RHO == 0.)) {
			Double d = new Double(Z[0]);
			int I01 = d.intValue() - 1;
			RHO = RHOTBL[I01] * (1.E-03);
		}

	}

	/**
	 * Set NE and initialize PZ, Z, WA, RHOZ.
	 * @param n n
	 */
	public static void setNE(int n) {
		NE = n;
		PZ = new double[NE];
		Z = new double[NE];
		WA = new double[NE];
		RHOZ = new double[NE];
	}

	/**
	 * Set RHO
	 * @param rhoe rhoe
	 */
	public static void setRHO(double rhoe) {
		RHO = rhoe;
	}

	// length of pze,ze,wae should always be equal with NE!!!!
	// SHOULD NEVER BE USED=OLD!!!!
	/**
	 * Set NE related variable from external sources. setNevar2 is more convenient to be used though.
	 * @param pze pze
	 * @param ze ze
	 * @param wae wae
	 * @param rhoze rhoze
	 */
	public static void setNEVAR(double[] pze, double[] ze, double[] wae,
			double[] rhoze) {
		for (int i = 0; i < NE; i++) {
			PZ[i] = pze[i];// e from external
			Z[i] = ze[i];// e from external
			WA[i] = wae[i];// e from external
			RHOZ[i] = rhoze[i];// e from external
		}
	}

	/**
	 * READ IN PHOTOELECTRIC AND PAIR DATA FROM FILE.
	 * @param filename filename
	 */
	public static void readPhotoelectricPairData(String filename) {
		String nume = "";
		nume = userdirs + file_sep + datas + file_sep + egsData + file_sep
				+ filename;// System.out.println(nume);
		crossPath = nume;
		int i = 0;
		int nphel = NPHE.length;// 100
		int nphec = 0;// contor
		int nphe = 0;// contor=nr of phe read for a specific Z!!
		int nz = 0;// contor=nr of Z read
		int nphd = 0;// contor=nr of phd read for a specific Z!!
		int nzd = 0;// contor=nr of corresponding Z read
		int nekedge = 0;// contor
		int npre = 0;// contor=nr of phe read for all Z!!
		int nprd = 0;
		int nzprd = 0;
		int nphc = 0;// contor=nr of cohe read for a specific Z!!
		int nzc = 0;// contor=nr of Z read
		StringBuffer desc = new StringBuffer();
		String descs = null;
		boolean haveData = false;

		try {
			FileInputStream in = new FileInputStream(nume);

			while ((i = in.read()) != -1) {
				if (!Character.isWhitespace((char) i)) {
					desc.append((char) i);
					haveData = true;
				} else {
					if (haveData)// we have data
					{
						haveData = false;// reset
						if (nphec < nphel)// read nphe
						{
							descs = desc.toString();// System.out.println(desc.toString());
							// NPHE[nphec]=Convertor.stringToInt(descs);//System.out.println(NPHE[nphec]);
							NPHE[nphec] = stringToInt(descs);// System.out.println(NPHE[nphec]);
							nphec++;// increment contor
						} else if ((nphe < $mxphote) && (nz < nphel))
						// PHE data: PHE= new double[$mxphote][nphel];
						{
							descs = desc.toString();// System.out.println(desc.toString());
							// PHE[nphe][nz]=Convertor.stringToDouble(descs);
							PHE[nphe][nz] = stringToDouble(descs);
							nphe++;
							if ((nphe == $mxphote) && (nz != nphel)) {
								nz++;// System.out.println(nz);
								nphe = 0;
							}
						} else if ((nphd < $mxphote) && (nzd < nphel))
						// PHD data PHD= new double[$mxphote][100];
						{
							descs = desc.toString();// System.out.println(desc.toString());
							// PHD[nphd][nzd]=Convertor.stringToDouble(descs);
							PHD[nphd][nzd] = stringToDouble(descs);
							nphd++;
							if ((nphd == $mxphote) && (nzd != nphel)) {
								nzd++;// System.out.println(nzd);
								nphd = 0;
							}
						} else if (nekedge < nphel)
						// EKEDGE=new double[100];nekedge
						{
							descs = desc.toString();// System.out.println(desc.toString());
							// EKEDGE[nekedge]=Convertor.stringToDouble(descs);//System.out.println(EKEDGE[nekedge]);
							EKEDGE[nekedge] = stringToDouble(descs);// System.out.println(EKEDGE[nekedge]);
							nekedge++;// increment contor
						} else if (npre < $mxpaire)
						// PRE= new double[$mxpaire];
						{
							descs = desc.toString();// System.out.println(desc.toString());
							// PRE[npre]=Convertor.stringToDouble(descs);//System.out.println(PRE[npre]);
							PRE[npre] = stringToDouble(descs);// System.out.println(PRE[npre]);
							npre++;// increment contor
						} else if ((nprd < $mxpaire) && (nzprd < nphel))
						// PRD= new double[$mxpaire][100];
						{
							descs = desc.toString();// System.out.println(desc.toString());
							// PRD[nprd][nzprd]=Convertor.stringToDouble(descs);
							PRD[nprd][nzprd] = stringToDouble(descs);
							nprd++;
							if ((nprd == $mxpaire) && (nzprd != nphel)) {
								nzprd++;// System.out.println(nzd);
								nprd = 0;
							}
						} else if ((nphc < $mxphote) && (nzc < nphel))
						// COHE= new double[$mxphote][100];RAYLEIGH data
						{
							descs = desc.toString();// System.out.println(desc.toString());
							// COHE[nphc][nzc]=Convertor.stringToDouble(descs);
							COHE[nphc][nzc] = stringToDouble(descs);
							nphc++;
							if ((nphc == $mxphote) && (nzc != nphel)) {
								nzc++;// System.out.println(nzd);
								nphc = 0;
							}
						} else {
							break;// if it is the case
						}
					}
					desc = new StringBuffer();
				}
			}
			in.close();
		} catch (Exception e) {

		}
	}

	/**
	 * READ IN ATOMIC FORM FACTOR DATA FROM A FILE.
	 * @param filename filename
	 */
	public static void readFormFactorData(String filename) {
		String nume = "";
		nume = userdirs + file_sep + datas + file_sep + egsData + file_sep
				+ filename;// System.out.println(nume);
		formPath = nume;
		int i = 0;
		int nx = 100;
		int ny = 100;
		int nphec = 0;// contor
		int nphe = 0;// contor=nr of XVAL read for a specific Z!!
		int nz = 0;// contor=nr of Z read
		StringBuffer desc = new StringBuffer();
		String descs = null;
		boolean haveData = false;

		try {
			FileInputStream in = new FileInputStream(nume);

			while ((i = in.read()) != -1) {
				if (!Character.isWhitespace((char) i)) {
					desc.append((char) i);
					haveData = true;
				} else {
					if (haveData)// we have data
					{
						haveData = false;// reset
						if (nphec < nx)// read nphe
						{
							descs = desc.toString();// System.out.println(desc.toString());
							// XVAL[nphec]=Convertor.stringToDouble(descs);//System.out.println(NPHE[nphec]);
							XVAL[nphec] = stringToDouble(descs);// System.out.println(NPHE[nphec]);
							nphec++;// increment contor
						} else if ((nphe < nx) && (nz < ny))
						// AFAC data:
						{
							descs = desc.toString();// System.out.println(desc.toString());
							// AFAC[nphe][nz]=Convertor.stringToDouble(descs);
							AFAC[nphe][nz] = stringToDouble(descs);
							nphe++;
							if ((nphe == nx) && (nz != ny)) {
								nz++;// System.out.println(nz);
								nphe = 0;
							}
						} else {
							break;// if it is the case
						}
					}
					desc = new StringBuffer();
				}
			}
			in.close();
		} catch (Exception e) {

		}
	}

	// ------------------------------------------------------------------------------------------
	/**
	 * Retrieve Z from table.
	 * @param IASYM IASYM
	 * @return the result
	 */
	public static double ZTBL(String IASYM) // retrieve Z from table
	{
		ZTBL = 0.0;
		for (int IE = 0; IE < NET; IE++) {
			if (IASYM.compareTo(ASYMT[IE]) == 0) {
				ZTBL = IE + 1.;// 0 is 1 for H!!!
			}
		}
		return ZTBL;
	}

	// ***SUBROUTINE TO ARRIVE AT PHYSICAL, MATHEMATICAL, AND DERIVED
	// CONSTANTS IN A VERY MNEMONIC WAY.
	/**
	 * COMPUTE PHYSICAL AND MATHEMATICAL CONSTANTS.
	 */
	public static void PMDCON() {
		// double FSCI = 0.0;
		RADDEG = 180. / PI;
		FSC = ECGS * ECGS / (HBAR * C);
		// FSCI = 1. / FSC;
		// 1.E+7 IS THE NUMBER OF ERGS PER JOULE
		ERGMEV = (1.E+6) * (EMKS * 1.E+7);
		R0 = (ECGS * ECGS) / (RME * C * C);
		RM = RME * C * C / ERGMEV;
		RMT2 = RM * 2.0;
		RMSQ = RM * RM;
		// CALCULATION OF SOME CONSTANTS
		A22P9 = RADDEG * Math.sqrt(4. * PI * AN) * ECGS * ECGS / ERGMEV;
		A6680 = 4.0 * PI * AN * (HBAR / (RME * C)) * (HBAR / (RME * C))
				* (0.885 * 0.885 / (1.167 * 1.13));
	}

	// for brems and differential pair production cross section, fc(Z) Z) is the
	// Coulomb correction,
	// for pair is for k (MeV) of initial photon>50MeV
	// Z is the atomic number
	// EVALUATES EQUATION 2.7.17 OF SLAC-265"
	/**
	 * For brems and differential pair production cross section, fc(Z) is the 
	 * Coulomb correction. fc(Z) is computed by this routine.
	 * @param Z Z
	 * @return the result
	 */
	public static double FCOULC(double Z) {
		double ASQ, FCOULC = 0.0;
		ASQ = (FSC * Z) * (FSC * Z);
		FCOULC = ASQ
				* (1.0 / (1.0 + ASQ) + 0.20206 + ASQ
						* (-0.0369 + ASQ * (0.0083 + ASQ * (-0.002))));
		return FCOULC;
	}

	// XSIF is a function which is used to take into account bremsstrahlung and
	// pair production in the field of the atomic electrons.
	// EVALUATES EQUATION 2.7.21a OF SLAC-265 on the basis of formulae: 2.7.19 &
	// 2.7.20b of SLAC-265!
	// see also Mustafa proof!!!
	/**
	 * XSIF is a function which is used to take into account bremsstrahlung and 
	 * pair production in the field of the atomic electrons. 
	 * EVALUATES EQUATION 2.7.21a OF SLAC-265 report.
	 * @param Z Z
	 * @return the result
	 */
	public static double XSIF(double Z) {
		double XSIF = 0.0;
		Double d = new Double(Z);// System.out.println(d);
		int IZ = d.intValue() - 1;// 0 biased for array!!!
		if (Z <= 4.0) {
			XSIF = ALRADP[IZ] / (ALRAD[IZ] - FCOULC(Z));
		} else
			XSIF = Math.log(A1440 * Math.pow(Z, -2. / 3.))
					/ (Math.log(A183 * Math.pow(Z, -1. / 3.)) - FCOULC(Z));
		return XSIF;
	}

	/**
	 * COMPUTE Z-DEPENDENT PARAMETERS FOR A MIXTURE.
	 */
	public static void MIX() {
		int IZZ = 0;// System.out.println(NE);
		double AL183 = 0.0;
		// "SCALE RHO FOR GASES BY THE GAS PRESSURE"
		if (GASP != 0.0) {
			RHO = GASP * RHO;
		}
		// " FIND VARIOUS SUMS AND SET OTHER VARIABLES
		AL183 = Math.log(A183);
		// /TPZ,WM,ZC,ZT,ZB,ZF,ZS,ZE,ZX,ZAB/=0.0;
		FZC = new double[NE];
		FCOUL = new double[NE];
		XSI = new double[NE];
		ZZX = new double[NE];
		ZZ = new double[NE];
		for (int I = 0; I < NE; I++) { // System.out.println(Z[I]);
			TPZ = TPZ + PZ[I];
			WM = WM + PZ[I] * WA[I];
			ZC = ZC + PZ[I] * Z[I];
			FZC[I] = (FSC * Z[I]) * (FSC * Z[I]);
			FCOUL[I] = FCOULC(Z[I]);
			XSI[I] = XSIF(Z[I]);
			ZZX[I] = PZ[I] * Z[I] * (Z[I] + XSI[I]);
			if (Z[I] <= 4.0) {
				Double d = new Double(Z[I]);
				IZZ = d.intValue() - 1;// 0 biased for array!!!;
				ZAB = ZAB + ZZX[I] * ALRAD[IZZ];// "IN THE CASE OF Z.LE.4 "
			} else {
				ZAB = ZAB + ZZX[I]
						* (AL183 + Math.log(Math.pow(Z[I], -1. / 3.)));// "Z.GT.4"
			}
			ZT = ZT + ZZX[I];
			ZB = ZB + ZZX[I] * Math.log(Math.pow(Z[I], -1. / 3.));
			ZF = ZF + ZZX[I] * FCOUL[I];
			ZZ[I] = PZ[I] * Z[I] * (Z[I] + $FUDGEMS);
			ZS = ZS + ZZ[I];// //eq 2.14.1 SLAC
			ZE = ZE + ZZ[I] * ((-2. / 3.) * Math.log(Z[I]));
			// same as the:
			// ZE = ZE + ZZ[I]*(Math.log(Math.pow(Z[I],-2./3.)));//eq 2.14.2
			// SLAC
			ZX = ZX + ZZ[I] * Math.log(1. + 3.34 * FZC[I]);// "END OF I LOOP"--eq
															// 2.14.3 SLAC
		}
		EZ = ZC / TPZ;
		ZA = AL183 * ZT;
		ZG = ZB / ZT;
		ZP = ZB / ZA;
		ZV = (ZB - ZF) / ZT;
		ZU = (ZB - ZF) / ZA;
		EDEN = AN * (RHO / WM) * ZC;
		RLC = 1. / ((AN * RHO / WM) * 4.0 * FSC * R0 * R0 * (ZAB - ZF));
		// "NOW COMPUTE MATERIAL-DEPENDENT CONSTANTS FOR MULTIPLE SCATTERING"
		// "LET B BE MOLIERE'S UPPER CASE B, AND LET BLC BE MOLIERE'S LOWER CASE B"
		// "THEN USING SCOTT'S NOTATION, A MEASURE OF THE NUMBER OF SCATTERS IS"
		// "OMEGA0=EXP(BLC).  NOW B AND BLC ARE RELATED BY" eq.2.14.17 SLAC-265
		// "BLC=B-ALOG(B);, OR EXP(BLC)=EXP(B)/B; " eq.2.14.16 SLAC-265
		// "NOW LET T BE THE TRANSPORT DISTANCE IN RL. THEN,"
		// "OMEGA0=EXP(BLC)"eq.2.14.18 - 2.4.19 SLAC-265
		// "      =( A6680*RHO*ZS*EXP(ZE/ZS)*RLC/(WM*EXP(ZX/ZS)) )*T/BETA**2"
		// "WHERE BETA IS V/C FOR THE PARTICLE"
		// "NOW SUPPOSE THAT XC IS CHI-SUB-C, THE CHARACTERISTIC ANGLE OF"
		// "MOLIERE'S THEORY.  THEN IT IS GIVEN BY" eq.2.14.23 SLAC-265
		// "XC= (A22P9/RADDEG)*SQRT( ZS*RHO*RLC/WM ) * SQRT(T) /(E*BETA**2) "
		// " WHERE E IS THE ENERGY OF THE PARTICLE IN MEV."
		// "NOW IN MOLIERE'S THEORY, LET XRM BE XC*SQRT(B).  THEN"
		// " (MOLIERE'S REDUCE ANGLE)=(REAL ANGLE)/XRM."
		// "FOR THE VARIABLE SAMPLING DONE IN EGS, A DIFFERENT REDUCED"
		// "ANGLE IS USED:"
		// " (REDUCED ANGLE FOR VARIABLE SAMPLING)=(REAL ANGLE)/XRV "
		// " WHERE NOW XRV=XC*SQRT(BLC)."
		// "THIS IS DONE BECAUSE BLC IS MORE SIMPLY RELATED TO THE TRANSPORT"
		// "DISTANCE THAN IS B."

		// "WITH THIS BACKGROUND WE NOW PROCEED TO COMPUTE"
		// "XR0, TEFF0, BLCC, AND XCC, WHICH ARE PARAMETERS NEEDED IN EGS"
		// "THESE PARAMETERS ARE USED AS FOLLOWS:(IN OUR NOTATION)"
		// "XRM=XR0*SQRT(T*B/(T0*BMIN))/(E*BETA) "
		// "  WHERE B IS EVALUATED FOR DISTANCE T, AND T0 IS THE"
		// "  DISTANCE THAT WOULD GIVE B=BMIN. T0 IS A FUNCTION OF BETA, AND"
		// "  IS GIVEN BY THE EXPRESSION  "
		// "T0=TEFF0*BETA**2;  THIS IS ALSO THE DEFINING RELATION FOR TEFF0"eq.2.14.52
		// SLAC-265
		// "BMIN IS THE MINIMUM VALUE OF B FOR WHICH MOLIERE CONSIDERED HIS THEORY"
		// "VALID.  WE CURRENTLY GIVE BMIN A VALUE OF 4.5."
		// "OMEGA0=BLCC*T/BETA**2  IS THE USEAGE OF BLCC" eq.2.14.50 SLAC-265
		// "XRV= XCC*SQRT( T*BLC )/(E*BETA**2) IS THE USAGE OF XCC"

		// "WE NOW PROCEED TO SOLVE FOR THE ABOVE QUANTITIES."
		// "COMPARING THE TWO EXPRESSIONS FOR OMEGA0, WE CONCLUDE THAT:"
		// see eq: 2.14.19 SLAC 265
		BLCC = A6680 * RHO * ZS * Math.exp(ZE / ZS) * RLC
				/ (WM * Math.exp(ZX / ZS));
		// "NOW USING THE EXPRESSION FOR OMEGA0 AT B=BMIN,WE HAVE"
		// " OMEGA0=EXP(BLCMIN)"
		// "       =EXP(BMIN)/BMIN"//ok see definition
		// O0=exp(blcmin)=exp(bmin-ln(bmin))
		// "       =BLCC*T0/BETA**2"//eq.2.14.50 SLAC-265
		// "       =BLCC*TEFF0    ,     THEREFORE"

		TEFF0 = (Math.exp(BMIN) / BMIN) / BLCC;// ok

		// "COMPARING THE TWO EXPRESSION FOR XRV WE CONCLUDE THAT"
		// "XC= XCC * SQRT(T) /(E*BETA**2)"
		// "COMPARING THIS WITH OUR PREVIOUS EXPRESSION FOR XC, WE OBTAIN"
		// see eq: 2.14.24 SLAC 265
		XCC = (A22P9 / RADDEG) * Math.sqrt(ZS * RHO * RLC / WM);

		// "NOW COMPARING THE TWO EXPRESSIONS FOR XRM WE OBTAIN"
		// "XRM=XC*SQRT(B)       USING DEFINITION OF XRM"
		// "   =(XCC*SQRT(T)/(E*BETA**2)) *SQRT(B)  AFTER SUBSTITUTING FOR XC"
		// eq.2.14.23 SLAC-265
		// "   =XR0*SQRT((T*B)/(T0*BMIN))/(E*BETA)  USING EXPRESSION DEFINING XR0"
		// "   =XR0*SQRT((T*B)/(TEFF0*BETA**2*BMIN))/(E*BETA)  EXPANDING T0"
		// " THUS AFTER SOME CANCELLATION AND SOLVING FOR XR0 WE OBTAIN"

		XR0 = XCC * Math.sqrt(TEFF0 * BMIN);// ok!

		// "THIS COMPLETE THE MS CALCULATIONS"
	}

	/**
	 * Initializes stopping power functions for a particular medium.
	 */
	public static void SPINIT() {
		boolean b = false;
		int IM = -100;// control var
		int IZ = 0;
		double IMEV = 0.0;
		double VPLASM = 0.0;
		double ALGASP = 0.0;
		double ALIADG = 0.0;
		double EDENL = 0.0;
		// double EPSTRH = 0.;// is RHO
		// double TLRNCE = 0.;
		// int ICHECK = 0;
		if ((EPSTFL < 0) || (EPSTFL > 1)) {
			EPSTFL = 0;// "ERROR ON INPUT, IGNORE", if 1 read density corr from
						// file
		}
		if (EPSTFL == 0) {
			if (ISSB == 0)// no user input data for stern coeff!!
			{
				// "CHECK TO SEE IF MATERIAL IS A 'STERNHEIMER-SELTZER-BERGER' (SSB)"
				// "MATERIAL, FOR WHICH THE DENSITY EFFECT PARAMETERS ARE ALREADY"
				// "SPECIFIED.  IF IT IS NOT, THEN CALCULATE THEM USING THE GENERAL"
				// "FORMULA BY STERNHEIMER-PEIERLS (S-P)."see
				// SLAC-265-eq:2.13.19
				// lookup for materials in default 73 elements table!!
				// "STERNHEIMER-SELTZER-BERGER (SSB) LOOKUP TABLE SECTION:"
				// GO TO :SSB-PARAMETERS-DEFINED:;]
				for (int I = 0; I < NMED; I++) {
					for (int J = 0; J < LMED; J++) {
						if (IDSTRN[J].compareTo(MEDTB[I]) == 0)// found!!
						{
							// "CALCULATION FOLLOWS IF A MATCH IS FOUND"
							AFACT = STDATA[I][0];// System.out.println("enter");
							SK = STDATA[I][1];
							X0 = STDATA[I][2];
							X1 = STDATA[I][3];
							IEV = STDATA[I][4];
							CBAR = STDATA[I][5];
							IMEV = IEV * 1.0E-6; // "EV TO MEV"
							VPLASM = Math.sqrt(EDEN * R0 * C * C / PI);// plasma
																		// frequency
																		// eq.2.13.22
																		// SLAC-265
							// :SSB-PARAMETERS-DEFINED:
							// "GAS PRESSURE CORRECTION COMES NEXT"eq.2.13.22
							// next SLAC-265
							if (GASP != 0.0) {
								ALGASP = Math.log(GASP);
								CBAR = CBAR - ALGASP;
								X0 = X0 - ALGASP / TOLN10;
								X1 = X1 - ALGASP / TOLN10;
							}
							if (IM == 0) {
								AFACT = (CBAR - TOLN10 * X0)
										/ Math.pow((X1 - X0), SK);
							}
							SPC1 = 2. * PI * R0 * R0 * RM * EDEN * RLC;
							SPC2 = Math.log((IMEV / RM) * (IMEV / RM) / 2.0);// eq
																				// 2.13.5
																				// SLAC-265

							return;

						} else {
							b = true;// and remain true!!!->not found at least
										// one component
						}
					}

				}
				if (b)// not found data in default 73 table, so we compute it!!!
				{
					IM = 0;// System.out.println("enter");
					// "DETERMINE THE MEAN EXCITATION ENERGY, IMEV (IN MEV)"
					if (NE == 1)// "ELEMENT"
					{
						Double d = new Double(Z[0]);
						IZ = d.intValue() - 1;// 0 biased Z=1->H->0 index in
												// ITBL
						// "I.E., DIATOMIC MOLECULE"
						// ' ELEMENT (H, N, OR O) CAN ONLY EXIST AS A DIATOMIC
						// MOLECULE.',/,
						// ' REMEDY: USE COMP OPTION FOR H2, N2, OR O2 WITH
						// NE=2,PZ=1,1'/,
						// ' AND, IN THE CASE OF A GAS, DEFINE STERNHEIMER
						// ID',/,
						// ' (I.E., IDSTRN) LIKE H2-GAS');

						IEV = ITBL[IZ]; // "EV"
					} else// "COMPOUND/MIXTURE---USE BRAGG ADDITIVITY RULE"
					{
						ALIADG = 0.0;
						for (int IE = 0; IE < NE; IE++) {
							Double d = new Double(Z[IE]);
							IZ = d.intValue() - 1;// 0 biased Z=1->H->0 index in
													// ITBL
							if (IZ == 0) {
								IEV = 19.2;
							} // "EV" H!!
							else if (IZ == 5) {
								if (GASP == 0.0) {
									IEV = 81.0;
								} else {
									IEV = 70.0;
								}
							} else if (IZ == 6) {
								IEV = 82.0;
							} else if (IZ == 7) {
								if (GASP == 0.0) {
									IEV = 106.0;
								} else {
									IEV = 97.0;
								}
							} else if (IZ == 8) {
								IEV = 112.0;
							} else if (IZ == 16) {
								IEV = 180.0;
							} else {
								IEV = 1.13 * ITBL[IZ];
							}
							// "NRCC comment - above 7 lines reflect table 6 in ref 59 of SLAC-265"
							// "       Berger and Seltzer's fudge to get better agreement with expt"
							ALIADG = ALIADG + PZ[IE] * Z[IE] * Math.log(IEV);
						}// end for
						ALIADG = ALIADG / ZC;
						IEV = Math.exp(ALIADG); // "EV"
					}// end compound
					IMEV = IEV * 1.0E-6;// "EV TO MEV"
					// "COMPUTE VARIOUS STERNHEIMER CONSTANTS"
					if (GASP == 0.0) {
						EDENL = EDEN;
					} else {
						// "VPLASM MUST BE FOR NTP FOR A GAS, AND EDEN HAS"
						// "      BEEN DEFINED IN MIX FOR THE ACTUAL PRESSURE"
						EDENL = EDEN / GASP;
					}
					VPLASM = Math.sqrt(EDENL * R0 * C * C / PI);
					// "ABOVE PATCHED JAN 9,1989 TO REFLECT ERROR POINTED OUT BY"
					// "PROF KAMAE, TOKYO UNIVERSITY, VIA HIDEO HIRAYAMA"
					CBAR = 1. + 2. * Math.log(IMEV
							/ (HBAR * 2 * PI * VPLASM / ERGMEV));
					// IFIX means return
					// Z[I]!!!!!!!!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
					if ((NE == 1) && ((Z[0] == 2.) && (GASP != 0.0))) {
						// "SPECIAL HE(GAS) CASE"
						X0 = 2.191;
						X1 = 3.0;
						SK = 3.297;
					} else if ((NE == 2) && ((Z[0] == 1.) && (Z[1] == 1.))) {
						if (GASP == 0.0) {
							// "SPECIAL H2(LIQUID) CASE"
							X0 = 0.425;
							X1 = 2.0;
							SK = 5.949;
						} else {
							// "SPECIAL H2(GAS) CASE"
							X0 = 1.837;
							X1 = 3.0;
							SK = 4.754;
						}
					} else {
						// "REGULAR CASES"
						SK = 3.0;
						if (GASP == 0.0) {
							// "SOLIDS AND LIQUIDS"
							if (IEV < 100.0) {
								if (CBAR < 3.681) {
									X0 = 0.2;
									X1 = 2.0;
								} else {
									X0 = 0.326 * CBAR - 1.0;
									X1 = 2.0;
								}
							} else {
								// "IEV GE 100.0"
								if (CBAR < 5.215) {
									X0 = 0.2;
									X1 = 3.0;
								} else {
									X0 = 0.326 * CBAR - 1.5;
									X1 = 3.0;
								}
							}
							// allways must be X0<X1
						} else {
							// "GASES---AT NTP AT THIS STAGE"
							if (CBAR < 10.0) {
								X0 = 1.6;
								X1 = 4.0;
							} else if (CBAR < 10.5) {
								X0 = 1.7;
								X1 = 4.0;
							} else if (CBAR < 11.0) {
								X0 = 1.8;
								X1 = 4.0;
							} else if (CBAR < 11.5) {
								X0 = 1.9;
								X1 = 4.0;
							} else if (CBAR < 12.25) {
								X0 = 2.0;
								X1 = 4.0;
							} else if (CBAR < 13.804) {
								X0 = 2.0;
								X1 = 5.0;
							} else {
								X0 = 0.326 * CBAR - 2.5;
								X1 = 5.0;
							}
						}// end GASES---AT NTP AT THIS STAGE"
							// ------------------------------------------
					}// regular case
						// "GAS PRESSURE CORRECTION COMES NEXT"
					if (GASP != 0.0) {
						ALGASP = Math.log(GASP);
						CBAR = CBAR - ALGASP;
						X0 = X0 - ALGASP / TOLN10;
						X1 = X1 - ALGASP / TOLN10;
					}

				}// end b not found
			}// ISSB==0)//no user input data for stern coeff!!
			if (IM == 0) {
				AFACT = (CBAR - TOLN10 * X0) / Math.pow((X1 - X0), SK);
			}
		}// "END OF EPSTFL=0 BLOCK"
		else// 1, read from file EPSTFL=1 BLOCK I.E. READ IN THE INPUT FROM THE
			// DENSITY"
			// "INPUT FILE"
		{
			// EPSTRH = RHO;
			// "CONVERT TO TOTAL ENERGIES"
			for (int I = 0; I < NEPST; I++) {
				EPSTEN[I] = EPSTEN[I] + RM;
			}
			IMEV = IEV * 1.E-06;
			// "CHECK THAT WE HAVE COVERED ENERGY RANGES NEEDED"
			if (AE < EPSTEN[0]) {
				// LOWEST ENERGY INPUT FOR DENSITY EFFECT IS HIGHER THAN THE
				// VALUE OF AE
				EPSTEN[0] = AE;// IT HAS BEEN SET TO AE
			}

			if (UE > EPSTEN[NEPST - 1]) {
				// HIGHEST ENERGY INPUT FOR DENSITY EFFECT IS LOWER THAN THE
				// VALUE OF UE
				EPSTEN[NEPST - 1] = UE;// IT HAS BEEN SET TO UE
			}
			// "DO A CHECK ON THE COMPOSITION AND DENSITY TO INSURE THE RIGHT DELTA"
			// "HAS BEEN PICKED UP. ALLOW A TOLERANCE OF 1 PERCENT ERROR ON THE"
			// "COMPOSITION BY WEIGHT."
			// ICHECK = 0; //
			// "FLAG GETS SET TO UNITY IF THE COMPOSITION DOES NOT MATCH"
			// TLRNCE = 0.01; // "TOLERANCE ALLOWED ON FRACTION BY WEIGHT"
			// "FIRST CHECK THAT THE NUMBER OF ELEMENTS ARE THE SAME"
			if (NELEPS != NE) {
				// ICHECK = 1;// never happen!!!
			}
			// check density but is RHO!!!!!!ETC....we stop now here all is
			// OK!!!
		}// "END OF EPSTFL=1 BLOCK"
		SPC1 = 2. * PI * R0 * R0 * RM * EDEN * RLC;
		SPC2 = Math.log((IMEV / RM) * (IMEV / RM) / 2.0);// eq 2.13.5 SLAC-265
	}

	// COMPUTE DIFFERENTIAL SAMPLING(BREMPR) CONSTANTS.
	/**
	 * COMPUTE DIFFERENTIAL SAMPLING(BREMPR) CONSTANTS.
	 */
	public static void DIFFER() {
		double AL183 = 0.;
		double F10 = 0.;
		double F20 = 0.;
		double A1DEN = 0.;
		double A2DEN = 0.;
		double B1DEN = 0.;
		double B2DEN = 0.;
		double C1DEN = 0.;
		double C2DEN = 0.;
		// "**********************************************************************"
		// "***THE DIFFERENTIAL CROSS-SECTIONS OF BREMSSTRAHLUNG AND PAIR       "
		// "   PRODUCTION ARE Z-DEPENDENT.  BUTCHER AND MESSEL USE A SAMPLING   "
		// "   TECHNIQUE FOR BREMSSTRAHLUNG WHICH ALSO DEPENDS ON THE LARGEST   "
		// "   AND SMALLEST ALLOWED ENERGY VALUES.                              "
		// "   THIS IS NOW TAKEN INTO ACCOUNT BY THE FUNCTION ILOG2(E/AP) WHICH "
		// "   GIVES THE NUMBER OF SUBDISTRIBUTIONS NEED TO PRODUCE PHOTONS DOWN"
		// "   TO THE LIMIT AP.                                                 "
		// "   THE DIFFERENTAL CROSS SECTIONS USED FOR BREMS AND PAIR ARE--     "
		// "   (THEY ARE CROSS SECTION PER RADIATION LENGTH)                    "

		// "***BELOW 50MEV (BETHE-HEITLER)--                                      "
		// "         BREMS--                                                      "
		// "     PHI1=(LOG(2)*(4/3+1/(9*LOG(A183)*(1+ZP)))* (1/LOG(2)*(1-EPS)/EPS) eq 2.7.63-2.7.64 SLAC 265"
		// "    1 * A(DELTAPRIME) + (1/2) * (2*EPS) * B(DELTAPRIME)               "
		// "         PAIR --                                                      "
		// "     PHI2=(2/3 - 1/(36*LOG(A183)*(1+ZP)))*(1) * C(DELTAPRIME)    eq 2.7.98-2.7.102 SLAC 265     "
		// "    1  + (1/12*(4/3+1/(9*LOG(A183)*(1+ZP))))* (12*(EPS-1/2)**2)       "
		// "    2 * A(DELTAPRIME)                                                 "
		// "       WHERE A,B,C ARE SCREENING REJECTION FUNCTIONS GIVEN BY         "
		// "       A(D)=(3*F1(D)-F2(D)+8*ZG)/(3*F1(0)-F2(0)+8*ZG) eq 2.7.60 SLAC 265 "
		// "       B(D)=(F1(D)+4*ZG)/(F1(0)+4*ZG)                 eq 2.7.61 SLAC 265 "
		// "       C(D)=(3*F1(D)+F2(D)+16*ZG)/(3*F1(0)+F2(0)+16*ZG)  eq 2.7.100 SLAC 265 -2.7.111"
		// "       AND WHERE                                                      "
		// "       DELTAPRIME= 136.*EXP(ZG)*RM*DEL= DELCM *DEL eq 2.7.37 SLAC 265 "
		// "       AND                                                            "
		// "       DEL =  EPS/(E*(1-EPS))      BREMSSTRAHLUNG  eq 2.7.36 SLAC 265 "
		// "           =  1/(E*EPS*(1-EPS))    PAIR PRODUCTION eq 2.7.96 SLAC 265 "
		// "       F1(D) AND F2(D) ARE SCREENING FUNCTIONS GIVEN APPROXIMATELY BY "
		// "             IF D.LE.1 ,THEN                                          "
		// "       F1(D) = 20.867 - 3.242*D + 0.625*D**2      eq 2.7.14 SLAC 265  "
		// "       F2(D) = 20.209 - 1.930*D - 0.086*D**2      eq 2.7.15 SLAC 265  "
		// "         BUT IF D.GT.1 ,THEN                                          "
		// "       F1(D)=F2(D)= 21.12 - 4.184*ALOG(D+0.952)                       "
		// "         IN ADDITION WE HAVE THAT                                     "
		// "       F1(0)= 4.*ALOG(A183)             eq 2.7.44 SLAC 265            "
		// "       F2(0)= F1(0) - 2./3.             eq 2.7.45 SLAC 265            "
		// "***ABOVE 50 MEV (COULOMB CORRECTED BETHE HEITLER)                     "
		// "         BREMS--                                                      "
		// "     PHI1=(LOG(2)*(4/3+1/(9*LOG(A183)*(1+ZU)))* (1/LOG(2)*(1-EPS)/EPS) eq 2.7.63-2.7.64 SLAC 265"
		// "    1 * A(DELTAPRIME) + (1/2) * (2*EPS) * B(DELTAPRIME)               "
		// "         PAIR --                                                      "
		// "     PHI2=(2/3-1/(36*LOG(A183)*(1+ZU)))* 1 * C(DELTAPRIME)         eq 2.7.98-2.7.102 SLAC 265   "
		// "    1  + (1/12*(4/3+1/(9*LOG(A183)*(1+ZU))))* (12*(EPS-1/2)**2)       "
		// "    2  * A(DELTAPRIME)                                                "
		// "       WHERE A,B,AND C ARE NOW GIVEN BY                               "
		// "       A(D) = (3*F1(D)-F2(D)+8*ZV)/(3*F1(0)-F2(0)+8*ZV)   eq 2.7.60 SLAC 265"
		// "       B(D) = (F1(D)+4*ZV)/(F1(0)+4*ZV)                   eq 2.7.61 SLAC 265"
		// "       C(D) = (3*F1(D)+F2(D)+16*ZV)/(3*F1(0)+F2(0)+16*ZV) eq 2.7.100 -2.7.111 SLAC 265"
		// "       AND DELTAPRIME,DEL,F1, AND F2 ARE THE SAME AS BELOW 50 MEV.    "
		// "   BREMSSTRAHLUNG ALPHA(I), I=1,...,N.  ALPHA(N+1)=0.5                "
		AL2 = Math.log(2.);
		AL183 = Math.log(A183);
		ALPHI[0] = AL2 * (4. / 3. + 1. / (9. * AL183 * (1. + ZP)));// eq 2.7.64
		ALPHI[1] = AL2 * (4. / 3. + 1. / (9. * AL183 * (1. + ZU)));
		// "*****PAIR PRODUCTION DIFFERENTIAL CROSS-SECTION NUMBERS               "
		ALFP1[0] = 2. / 3. - 1. / (36. * AL183 * (1. + ZP));// eq 2.7.102
		ALFP1[1] = 2. / 3. - 1. / (36. * AL183 * (1. + ZU));
		ALFP2[0] = (1. / 12.) * (4. / 3. + 1. / (9. * AL183 * (1 + ZP)));// eq
																			// 2.7.105
		ALFP2[1] = (1. / 12.) * (4. / 3. + 1. / (9. * AL183 * (1 + ZU)));
		// "*****BRANCHING RATIO FOR PAIR PRODUCTION SAMPLING                     "
		BPAR[0] = ALFP1[0] / (ALFP1[0] + ALFP2[0]);
		BPAR[1] = ALFP1[1] / (ALFP1[1] + ALFP2[1]);
		// "*****WE MUST ALWAYS HAVE A,B,C(DELTA) POSITIVE. SHOWER ASSURES THIS BY"
		// "     REQUIRING DEL.LT.DELPOS(),AND HENCE THAT DELTA.LT.DELCM *DELPOS()"
		// "     THE CROSS-OVER POINT IS IN THE DELTA.GT.1 REGION,AND THE CROSSING"
		// "     CONDITION REDUCES TO --                                          "
		// "       F1(D) + 4*ZG = 0     IF  E.LT.50 MEV                           "
		// "       F1(D) + 4*ZV = 0     IF  E.GT.50 MEV                           "
		// "     OR --                                                            "
		// "       21.12 - 4.184*ALOG(DELCM *DELPOS(1)+0.952) +4.*ZG =0     E.LT.5"
		// "       21.12 - 4.184*ALOG(DELCM *DELPOS(2)+0.952) +4.*ZV =0     E.GT.5"
		// "       THUS DELPOS IN THE TWO ENERGY REGIONS IS GIVEN BY              "
		DELCM = 136.0 * Math.exp(ZG) * RM;
		DELPOS[0] = (Math.exp((21.12 + 4. * ZG) / 4.184) - 0.952) / DELCM;// eq
																			// 2.7.31
																			// SLAC
																			// 265
		DELPOS[1] = (Math.exp((21.12 + 4. * ZV) / 4.184) - 0.952) / DELCM;// ZG<->ln(Z^-1/3)
		// "     CALCULATION OF THE BETHE-HEITLER(B-H) A,B,C(DELTA).AND THE       "
		// "     COULUMB CORRECTED(CC) A,B,C(DELTA).  THE SHOWER PROGRAM COMPUTES "
		// "     THESE SIX SCREENING REJECTION FUNCTION WITH THE EXPRESSIONS--    "
		// "       DL1(LVL)+DELTA*(DL2(LVL)+DELTA*DL3(LVL))   IF DELTA.LT.1 ,AND B"
		// "       DL4(LVL)+DL5(LV)*ALOG(DELTA+DL6(LVL))      IF DELTA.GE.1       "
		// "       WHERE LVL IS USED TO SELECT THE FUNCTION--                     "
		// "       LVL=1   B-H A(DELTA)                                           "
		// "          =2   B-H B(DELTA)                                           "
		// "          =3   B-H C(DELTA)                                           "
		// "          =4   CC  A(DELTA)                                           "
		// "          =5   CC  B(DELTA)                                           "
		// "          =6   CC  C(DELTA)                                           "
		// "     FIRST COMPUTE THE DENOMINATORS                                   "
		F10 = 4. * AL183;
		F20 = F10 - 2. / 3.;
		A1DEN = 3.0 * F10 - F20 + 8.0 * ZG;
		A2DEN = 3.0 * F10 - F20 + 8.0 * ZV;
		B1DEN = F10 + 4.0 * ZG;
		B2DEN = F10 + 4.0 * ZV;
		C1DEN = 3.0 * F10 + F20 + 16.0 * ZG;
		C2DEN = 3.0 * F10 + F20 + 16.0 * ZV;
		// "     LVL=1, B-H A(DELTA)                                              "
		DL1[0] = (3.0 * 20.867 - 20.209 + 8.0 * ZG) / A1DEN;
		DL2[0] = (3.0 * (-3.242) - (-1.930)) / A1DEN;
		// DL3[0]= (3.0*(0.625)-(0.086))/A1DEN;
		DL3[0] = (3.0 * (0.625) - (-0.086)) / A1DEN;// see SLAC!!!-CORRECT
		DL4[0] = (2.0 * 21.12 + 8.0 * ZG) / A1DEN;
		DL5[0] = 2.0 * (-4.184) / A1DEN;
		DL6[0] = 0.952;
		// "     LVL=4, CC  A(DELTA)                                              "
		DL1[3] = (3.0 * 20.867 - 20.209 + 8.0 * ZV) / A2DEN;
		DL2[3] = (3.0 * (-3.242) - (-1.930)) / A2DEN;
		// DL3[3]= (3.0*(0.625)-(0.086))/A2DEN;
		DL3[3] = (3.0 * (0.625) - (-0.086)) / A2DEN;// see SLAC!!!-CORRECT
		DL4[3] = (2.0 * 21.12 + 8.0 * ZV) / A2DEN;
		DL5[3] = 2.0 * (-4.184) / A2DEN;
		DL6[3] = 0.952;
		// "     LVL=2, B-H B(DELTA)   CORRECT                                           "
		DL1[1] = (20.867 + 4.0 * ZG) / B1DEN;
		DL2[1] = -3.242 / B1DEN;
		DL3[1] = 0.625 / B1DEN;
		DL4[1] = (21.12 + 4.0 * ZG) / B1DEN;
		DL5[1] = -4.184 / B1DEN;
		DL6[1] = 0.952;
		// "     LVL=5, CC B(DELTA)    CORRECT                                           "
		DL1[4] = (20.867 + 4.0 * ZV) / B2DEN;
		DL2[4] = -3.242 / B2DEN;
		DL3[4] = 0.625 / B2DEN;
		DL4[4] = (21.12 + 4.0 * ZV) / B2DEN;
		DL5[4] = -4.184 / B2DEN;
		DL6[4] = 0.952;
		// "     LVL=3, B-H C(DELTA)           CORRECT!!                                   "
		DL1[2] = (3.0 * 20.867 + 20.209 + 16.0 * ZG) / C1DEN;
		DL2[2] = (3.0 * (-3.242) + (-1.930)) / C1DEN;
		DL3[2] = (3.0 * 0.625 + (-0.086)) / C1DEN;
		DL4[2] = (4.0 * 21.12 + 16.0 * ZG) / C1DEN;
		DL5[2] = 4.0 * (-4.184) / C1DEN;
		DL6[2] = 0.952;
		// "     LVL=6, CC  C(DELTA)           CORRECT!!                                   "
		DL1[5] = (3.0 * 20.867 + 20.209 + 16.0 * ZV) / C2DEN;
		DL2[5] = (3.0 * (-3.242) + (-1.930)) / C2DEN;
		DL3[5] = (3.0 * 0.625 + (-0.086)) / C2DEN;
		DL4[5] = (4.0 * 21.12 + 16.0 * ZV) / C2DEN;
		DL5[5] = 4.0 * (-4.184) / C2DEN;
		DL6[5] = 0.952;

	}

	// "*****HERE FOR OPTION TO SET ENERGY LIMITS FOR ELECTRONS AND PHOTONS"
	/**
	 * SET ENERGY LIMITS FOR ELECTRONS AND PHOTONS.
	 */
	public static void setENER() {
		if (AE < 0)
			AE = -AE * RM;
		if (UE < 0)
			UE = -UE * RM;
		if (AP < 0)
			AP = -AP * RM;
		if (UP < 0)
			UP = -UP * RM;
		TE = AE - RM;// incoming cuttoff electron energy
		TET2 = TE * 2.0;
		TEM = TE / RM;
		THBREM = RM + AP;// "THRESHOLD FOR DISCRETE BREMSSTRAHLUNG" CORRECT!!
		THMOLL = AE + TE;// see eq 2.10.2 SLAC-265
	}

	// "COMPUTE MATERIAL INDEPENDENT MOLIERE DATA
	/**
	 * COMPUTE MATERIAL INDEPENDENT MOLIERE DATA
	 */
	public static void MOLIER() {
		// "     DATA FROM BETHE'S TABLE - MOLIERE Functions!! -  "
		// double[] TH={.05,.2,.4,.6,.8,1.,1.2,1.4,1.6,1.8,2.,2.2,2.4,2.6,2.8,
		// 3.,3.2,3.4,3.6,3.8,4.07,4.5,5.,5.5,6.13,7.,8.,9.,9.75};
		// double[] DTH={.1,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2,
		// 0.2,0.2,0.2,0.2,0.2,0.2,0.35,0.5,0.5,0.5,0.75,1.0,1.0,1.0,0.5};
		// double[] F0={2.,1.9216,1.7214,1.4094,1.0546,.7338,.4738,.2817,
		// .1546,.0783,.0366,.01581,.0063,.00232,7.9E-4,2.5E-4,7.3E-5,
		// 1.9E-5,4.7E-6,1.1E-6,2.3E-7,3.E-9,2.E-11,2.E-13,5.E-16,1.E-21,
		// 3.E-28,1.E-35,$F029};
		// double[] F1={.8456,.7038,.3437,-0.0777,-0.3981,-0.5285,-0.4770,
		// -.3183,-.1396,-.0006,+0.0782,.1054,.1008,.08262,.06247,.0455,
		// .03288,.02402,.01791,.01366,.010638,.00614,.003831,.002527,
		// .001739,.000908,.0005211,.0003208,.0002084};
		// double[] F2={2.4929,2.0694,1.0488,-.0044,-.6068,-.6359,-.3086,.0525,
		// .2423,.2386,.1316,.0196,-.0467,-.0649,-.0546,-.03568,-.01923,
		// -.00847,-.00264,5.E-5,.0010741,.0012294,.0008326,.0005368,
		// .0003495,.0001584,7.83E-5,4.17E-5,2.37E-5};

		double BLCMIN = 0.;
		double BLC = 0.;
		double B = 0.;
		double BOLD = 0.;
		// double [] BLCA=new double[$MSSTEPS];
		// double [] BA=new double[$MSSTEPS];
		double PTOT = 0.;
		double B1 = 0.;
		double[][] P = new double[29][$MSSTEPS];
		double[][] Q = new double[29][$MSSTEPS];
		// int I=0;
		int L = 0;
		double PPP = 0.;
		double PP = 0.;
		int ITOT = 0;
		int[][] IP1 = new int[29][$MSSTEPS];
		int IDIF = 0;
		boolean b = false;
		int[] IP2 = new int[29];
		int N = 0;
		int I01 = 0;
		int I02 = 0;
		int ISWP = 0;
		int IFLG = 0;
		int IDA = 0;
		int INC = 0;
		int[] IALL = new int[29];
		int II = 0;
		int[][] IXTR = new int[29][$MSSTEPS];

		// System.out.println(F0.length+"    "+F0[28]);
		// "***THE PURPOSE OF THIS ROUTINE IS TO CALCULATE THE DISTRIBUTION OF  "
		// "   THE MOLIERE REDUCED ANGLES FOR USE IN MULTIPLE SCATTERING        "
		// "   THIS IS DONE FOR A SET OF $MSSTEPS STEPSIZES.  EACH STEPSIZE IS A"
		// "   FIXED MULTIPLE OF T0,THE STEP SIZE WHICH GIVES MOLIERE'S UPPER CASE"
		// "   B VARIABLE, B, ITS MINIMUM VALID VALUE BMIN.  BLC IS MOLIERE'S"
		// "   LOWER CASE B, WHICH THEN ALSO HAS A MINIMUM VALUE BLCMIN, WHICH IS"
		// "   RELATED TO BMIN BY THE TRANSCENDENTAL EQUATION LINKING B AND BLC,"
		// "   I.E., BLC=B-ALOG(B);  AND BLCMIN=BMIN-ALOG(BMIN);"
		// "   NOW BLC IS RELATED TO THE TRANSPORT DISTANCE THROUGH THE EQUATION"
		// "   EXP(BLC) = BLCC * T/BETA**2;  (SEE NOTES IN SUBROUTINE MIX)"
		// "   THUS IF WE DEFINE AN EFFECTIVE DISTANCE "
		// "   TEFF=T/BETA**2;  AND ALSO  TEFF0=T0/BETA**2;"
		// "   THEN THE ENERGY DEPENDENCE IF REMOVED FROM THE EQUATION,AND WE HAVE"
		// "   EXP(BLC)=BLCC*TEFF;    FOR THE MINIMAL CONDITION THIS BECOMES"
		// "   EXP(BLCMIN)=BLCC*TEFF0;"
		// "   NOW SUPPOSE WE HAVE A TRANSPORT DISTANCE T WHOSE RATIO TO T0 IS"
		// "   R = T/T0 = TEFF/TEFF0;     WE THEN HAVE"
		// "   EXP(BLC)=BLCC*TEFF=BLCC*TEFF0*(TEFF/TEFF0)"
		// "          =(BLCC*TEFF0)*R"
		// "          = EXP(BLCMIN)*R;     HENCE WE OBTAIN"
		// "   BLC = BLCMIN + ALOG(R); "
		// "   AND THE CORESPONDING VALUE OF B CAN BE DETERMINED BY SOLVING THE"
		// "   TRANSCENDENTAL EQUATION."
		// "   IN WHAT FOLLOWS, FSTEP(IS) CONTAINS THE RATIO OF THE SIZE OF THE"
		// "   IS THE STEPSIZE TO T0.  ITS VALUES ARE ALL INTEGERS.  WHEN USED BY"
		// "   EGS THE RATIO R OF THE DESIRED STEP SIZE TO T0 IS COMPUTED,AND THEN"
		// "   IT IS DESIRED TO FIND THE LARGEST FIXED STEP SIZE THAT IS LESS THAN"
		// "   R*T0.  TO FACILITATE THIS, WE DEFINE AN ARRAY MSMAP, OF SIZE"
		// "   $MXJREFF,SUCH THAT MSMAP(JR) IS THE INDEX OF THE LARGEST FIXED STEP"
		// "   LESS THAN OR EQUAL TO JR.  THUS IN EGS THE SELECTION OF A STEP IS:"
		// "   R=T/T0; JR=R; IS=MSMAP(JR) ; T=FSTEP(IS)*T0; "
		// "   IN ADDITION IT IS REQUIRED TO KNOW THE RATIO OF MOLIERE REDUCING"
		// "   ANGLE XRM=XC*SQRT(B) TO ITS MINIMUM VALUE XRMMIN=XCMIN*SQRT(BMIN);"
		// "   SINCE AT FIXED ENERGY XC IS PROPORTIONAL TO T, THIS RATIO IS"
		// "   FSQR(IS)= R*SQRT(B/BMIN);  THEN"@@@@@@@IN FACT XC~SQRT(T) see
		// MIX!!
		// "   XRM=FSQR*XRMMIN=FSQR*XCMIN*SQRT(BMIN)"
		// "    =FSQR* (XCC*SQRT(T0)/(E*BETA**2)) * SQRT(BMIN) "
		// "    =FSQR* (XCC*SQRT(TEFF0*BETA**2)/(E*BETA**2)) *SQRT(BMIN) "
		// "    = (XCC*SQRT(TEFF0*BMIN)) * FSQR / (E*BETA) "
		// "    = XR0 * FSQR / (E*BETA)"//@@@@@@@@@CORECT!!!
		// "   THIS LAST LINE IS THE FORM USED BY EGS.  XR0 CONTAINS THE MATERIAL"
		// "   DEPENDENCE, FSQR THE STEP DEPENDENCE, AND 1/(E*BETA)"
		// "   THE ENERGY DEPENDENCE."
		// "   IN THE CASE WHEN R.LT.1., THE PROPOSED TRANSPORT DISTANCE IS"
		// "   THAN THE SMALLEST VALID MOLIERE TRANSPORT DISTANCE, AND STRICTLY"
		// "   SPEAKING MOLIERE'S THEORY IS NO LONGER APPLICABLE.  SINCE IN MOST"
		// "   THE SCATTERING OVER THIS SHORT DISTANCE WILL NOT BE TOO SIGIFICANT,"
		// "   WE CURRENTLY USE THE REDUCED ANGLE DISTRIBUTION FOR T0 IN THIS"
		// "   CASE, AND LET FSQR=SQRT(R), I.E. IGNORE CHANGE IN B FROM BMIN."
		// @@@@@@@@@@@@@@XC~SQRT(T) see MIX!!
		// "***NOW FILL UP MSMAP."
		// HERE: MSMAP(JR) IS THE INDEX OF THE LARGEST FIXED STEP"
		// " LESS THAN OR EQUAL TO
		// JR+1@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@!!
		// THUS IN EGS THE SELECTION OF A STEP IS:"
		// "   R=T/T0; JR=R-1@@@@@@@@@@!!!!!!!!; IS=MSMAP(JR) ; T=FSTEP(IS)*T0; "
		for (int IS = 0; IS < MSTEPS - 1; IS++) {
			Double d1 = new Double(FSTEP[IS] - 1);// FSTEP[0]=1->0
			Double d2 = new Double(FSTEP[IS + 1]);
			for (int J = d1.intValue(); J < d2.intValue() - 1; J++) {
				// index is wanted
				MSMAP[J] = IS;// System.out.println(MSMAP[J]);
			}
		}
		// last index!![15]->in FSTEP[15]=200!!that is why is -1 here
		MSMAP[JRMAX - 1] = MSTEPS - 1;// System.out.println(MSMAP[JRMAX-1]);
		BLCMIN = BMIN - Math.log(BMIN);
		// "LOOP OVER STEP SIZES" :STEP:
		for (int IS = 0; IS < MSTEPS; IS++) {
			b = false;
			BLC = BLCMIN + Math.log(FSTEP[IS]);
			// "***NOW SOLVE FOR B"
			B = BLC + Math.log(BLC); // "THIS IS FIRST GUESS",BLC lower case is
										// b, B is B!!
			BOLD = B;// initialize
			B = BOLD - (BOLD - Math.log(BOLD) - BLC) / (1.0 - 1.0 / BOLD);
			while (Math.abs((B - BOLD) / BOLD) > 1.E-5)// UNTIL NEWTON'S METHOD
														// CONVERGES
			{
				BOLD = B;
				B = BOLD - (BOLD - Math.log(BOLD) - BLC) / (1.0 - 1.0 / BOLD);
			}
			// "SAVE BLC AND B FOR PRINTOUT"
			BLCA[IS] = BLC;
			BA[IS] = B;
			// "GET XRM RATIO"
			FSQR[IS] = Math.sqrt(FSTEP[IS] * B / BMIN);// XC~SQRT(T) see MIX!!
			// "***CALCULATE PROBABILITIES OF INPUT REDUCED ANGLES                  "
			B1 = 1.0 / B;// eq 2.14.7 SLAC 265
			PTOT = 0.0;
			for (int I = 0; I < 29; I++) {
				P[I][IS] = TH[I] * DTH[I] * (F0[I] + B1 * (F1[I] + B1 * F2[I]));// //eq
																				// 2.14.7
																				// SLAC
																				// 265
				// "***THESE ANGLES CORRESPOND TO THE PROBABILITY ENTRIES               "
				PTOT = PTOT + P[I][IS];
			}

			for (int I = 0; I < 29; I++) {
				P[I][IS] = P[I][IS] / PTOT;// "NORMALIZE"
			}
			// "***MAKE ALL PROBABILITIES EITHER 0.0 OR GREATER THAN 0.001 IN THE   "
			// "   MANNER NAGEL SUGGESTS.  NOTE THE PROBABILITY IS LUMPED INTO      "
			// "   THE ENTRY FOR THE LARGEST REDUCED ANGLE WHEN INTERVALS ARE MERGED"
			for (int I = 0; I < 29; I++) {
				Q[I][IS] = P[I][IS];// "DON'T CHANGE EXACT PROBABILITIES"
				// System.out.println(Q[I][IS]);
			}

			int I = 29 - 1;// 0 to 28 means 1 to 29 in fortran!!!
			while (I >= 0) {
				L = 1;// same as above
				while ((Q[I][IS] < 0.001) && (I >= L))// I>L is general=> so
														// I-L-1>=0
				{
					Q[I][IS] = Q[I][IS] + Q[I - L][IS];// [0] is also computed
					Q[I - L][IS] = 0.0;// [0] is also computed
					L = L + 1;// increment and asure the salt to previous
								// unresolved line!!
				}
				// I=I-L-1;//same as above
				I = I - L;// same as above
			}

			// for(I=0;I<29;I++)//---------WORKS!!!
			// {
			// System.out.println(" after "+Q[I][IS]);
			// }
			// System.out.println("----------------------------------------------");
			// "***NOW TRY TO MAKE UP 1000 INTEGRAL PROBABILITIES FOR VERT TABLES"
			PPP = 0.5;
			PP = 0.5; // "INITIALIZE ROUNDING PARAMETERS"
			for (int JLR = 1; JLR <= 10; JLR++)// <=10--test cu 4
			{ // "TRY FUDGING ROUNDING DOWN TO 2**-10"
				ITOT = 0;// "INITIALIZE TOTAL"
				for (I = 0; I < 29; I++) {// "ROUND PROBABILITIES FOR THIS ANGLE"
					IP1[I][IS] = new Double(Q[I][IS] * 1000.0 + PP).intValue();
					ITOT = ITOT + IP1[I][IS];
				}
				IDIF = ITOT - 1000;// System.out.println(IDIF);
				if (IDIF == 0) {
					b = true;
					break;
				} else {
					PPP = PPP * 0.5;
					if (IDIF < 0) {
						PP = PP + PPP;
					} else {
						PP = PP - PPP;
					}
				}
			}// end fudging
				// for(I=0;I<29;I++)//---------WORKS!!!
				// {
				// System.out.println(" before "+IP1[I][IS]);
				// }

			if (!b) {
				// System.out.println("enter");la JLR=10 never enter
				// "***GET 1000 ENTRIES WHEN ROUNDING ADJUSTMENT FAILS.  WE ADD OR      "
				// "   SUBTRACT ONE FROM AS MANY OF THE LARGEST PROBABILITIES AS NEEDED "
				// " BUBBLE SORT THE PROBABILITIES TO RANK THEM BY SIZE
				for (I = 0; I < 29; I++) {
					// IP2[I]=1;//error in PEGS4, MOLIER subroutine
					IP2[I] = I;
				}
				N = 29;
				IFLG = 1;// to enter in loop
				while (IFLG != 0)// "LOOP UNTIL PASS HAS NO SWAPS"
				{
					N = N - 1;
					IFLG = 0; // "PREPARE FOR ONE PASS THOUGH BUBBLE SORT"
					for (int J = 0; J < N; J++) {
						I01 = IP2[J];
						I02 = IP2[J + 1];
						if (IP1[I01][IS] < IP1[I02][IS]) {// "MUST SWAP"-look
															// for the LARGEST
															// PROBABILITIES
							ISWP = IP2[J];
							IP2[J] = IP2[J + 1];
							IP2[J + 1] = ISWP;
							IFLG = 1;// "SWAP&SET FLAG"
						}
					}// "END OF PASS"
				}

				if (IDIF < 0) {
					IDA = -IDIF;
					INC = 1;
				} else {
					IDA = IDIF;
					INC = -1;
				}
				for (I = 0; I < IDA; I++) {
					I01 = IP2[I];
					IP1[I01][IS] = IP1[I01][IS] + INC;
				}
			}// b

			// for(I=0;I<29;I++)//---------WORKS!!!
			// {
			// System.out.println(" after "+IP1[I][IS]);
			// }
			// System.out.println("----------------------------------------------");
		}
		// "***NOW FIND MINIMUM FREQUENCY FOR EACH ANGLE OVER ALL STEPS"
		MXV1 = 0;// "MXV1 WILL BE SUM OF MINIMUMS"
		for (int I = 0; I < 29; I++) {
			IALL[I] = IP1[I][0];
			for (int IS = 1; IS < MSTEPS; IS++) {
				IALL[I] = Math.min(IALL[I], IP1[I][IS]);
			}
			MXV1 = MXV1 + IALL[I];// "ADD IN MINIMUM FREQUENCY FOR THIS ANGLE"
		}// "END LOOP OVER ANGLES"
			// System.out.println(MXV1);
		MXV2 = 1000 - MXV1;// "THE VERT2 ARRAY HAS THE LEFT-OVERS"
		// System.out.println(MXV2);
		// "***&&&&&&&&&&LATER ONLY PUT IALL,IXTR TO EGS"
		// "   AND LET EGS FILL IN VERT1 AND VERT2. THIS REQUIRES LESS STORAGE"
		// "***NOW FILL UP THE VERT ARRAYS."
		II = 0;
		for (int I = 0; I < 29; I++) {
			for (int J = 1; J <= IALL[I]; J++) {
				// II=II+1;
				VERT1[II] = TH[I];// System.out.println(VERT1[II]);
				II = II + 1;
			}
		}
		// "***NOW DO VERT2"
		for (int IS = 0; IS < MSTEPS; IS++) {
			II = 0;
			for (int I = 0; I < 29; I++) {
				IXTR[I][IS] = IP1[I][IS] - IALL[I];// "GET NUMBER OF EXTRAS"
				for (int J = 1; J <= IXTR[I][IS]; J++) {
					// II=II+1;
					VERT2[II][IS] = TH[I];// System.out.println("ii "+II+" is "+IS+"   "+VERT2[II][IS]);
					II = II + 1;
				}
			}// "END ANGLE DO LOOP"
		} // "END STEP DO LOOP"

		// "***WE HAVE NOW COMPUTED EVERYTHING THE USER NEEDS, SO NOW LIST IT OUT"
	}

	// "*****FIRST FIT ELECTRON(AND POSITRON)"
	/**
	 * FIT ELECTRON(AND POSITRON)
	 */
	public static void FITEP() {
		EBINDA = EBIND(AP);
	}

	// "***FUNCTION TO GET AN AVERAGE PHOTOELECTRIC BINDING ENERGY.  
	/**
	 * FUNCTION TO GET AN AVERAGE PHOTOELECTRIC BINDING ENERGY.
	 * @param E E
	 * @return the result
	 */
	public static double EBIND(double E) {
		double EBIND = 0.;
		double STOT = 0.;
		int J = 0;
		for (int I = 0; I < NE; I++) {
			Double d = new Double(Z[I]);
			J = d.intValue() - 1;// 0 biased
			EBIND = EBIND + PZ[I] * PHOTTZ(Z[I], E) * EKEDGE[J] * 0.001;// EKEDGE
																		// in
																		// keV->MEV
		}
		STOT = PHOTTE(E);
		if (STOT != 0.0)
			EBIND = EBIND / STOT;

		return EBIND;
	}

	// "***'E' FOR EMPIRICAL OR EXPERIMENTAL PHOTOEFFECT CROSS-SECTION.     "
	// "   THIS FUNCTION GIVE THE PROPER MIX OF PHOTTZ'S.                   "
	/**
	 * THIS FUNCTION GIVE THE PROPER MIX OF PHOTTZ'S. 
	 * 'E' STANDS FOR EMPIRICAL OR EXPERIMENTAL PHOTOEFFECT CROSS-SECTION.
	 * @param K K
	 * @return the result
	 */
	public static double PHOTTE(double K) {
		double PHOTTE = 0.0;
		for (int I = 0; I < NE; I++) {
			PHOTTE = PHOTTE + PZ[I] * PHOTTZ(Z[I], K);
		}
		return PHOTTE;
	}

	// "***FUNCTION TO GIVE INTERPOLATED TOTAL PHOTO EFFECT CROSS-SECTION   "
	// "   FROM TABULATED DATA USING LOG-LOG INTERPOLATION.                 "
	// "   THIS CROSS-SECTION WILL BE THE CONTRIBUTION THIS Z ELEMENT       "
	// "   WOULD MAKE IF THERE WERE ONE/MOLECULE,I.E., TO GET THE ACTUAL    "
	// "   CONTRIBUTION ONE MUST MULTIPLY BY PZ(I).                         "
	// "   COMMON WITH DATA FOR PHOTO AND PAIR.                             "
	/**
	 * FUNCTION TO GIVE INTERPOLATED TOTAL PHOTO EFFECT CROSS-SECTION 
	 * FROM TABULATED DATA USING LOG-LOG INTERPOLATION. THIS CROSS-SECTION WILL BE THE CONTRIBUTION THIS Z ELEMENT 
	 * WOULD MAKE IF THERE WERE ONE/MOLECULE,I.E., TO GET THE ACTUAL CONTRIBUTION ONE MUST MULTIPLY BY PZ(I). 
	 * COMMON WITH DATA FOR PHOTO AND PAIR.
	 * @param Z Z
	 * @param K K
	 * @return the result
	 */
	public static double PHOTTZ(double Z, double K) {
		double PHOTTZ = 0.;
		// " PCON CONVERTS CROSS SECTION IN BARNS TO CROSS SECTION IN RLC**-1 "
		double PCON = 1.E-24 * (AN * RHO / WM) * RLC;
		Double d = new Double(Z);
		int IZ = d.intValue() - 1;// 0 biased
		// construct the doubles!!
		double[] xa = new double[NPHE[IZ]];
		double[] ya = new double[NPHE[IZ]];
		for (int i = 0; i < NPHE[IZ]; i++) {
			xa[i] = PHE[i][IZ];
			ya[i] = PHD[i][IZ];// System.out.println("Z: "+Z+" xa=  "+xa[i]+" ya=  "+ya[i]);
		}
		PHOTTZ = PCON * AINTP(K, xa, NPHE[IZ], ya, 1, true, true);
		return PHOTTZ;
	}

	/**
	 * Linear interpolation based on EGSnrc fortran/mortran code. 
	 * @param X the desired x-value
	 * @param XA the known array of x-values
	 * @param NX its dimension (length of the array)
	 * @param YA the array of y-values.
	 * @param ISK dummy variable having no use!!
	 * @param XLOG true if we want logarithmic values for x-values
	 * @param YLOG true if we want logarithmic values for y-values 
	 * @return the desired y-value
	 */
	public static double AINTP(double X, double[] XA, int NX, double[] YA,
			int ISK, boolean XLOG, boolean YLOG) {
		double AINTP = 0.;
		boolean XLOGL = XLOG;// "SET LOCAL VARIABLE"
		int I = 0;
		int J = 0;
		boolean ok = false;
		double XI = 0.;
		double XJ = 0.;
		double XV = 0.;
		double YI = 0.;
		double YJ = 0.;
		// " FIND INTERVAL FOR X INTERPOLATION
		for (J = 1; J < NX; J++) {
			if (X < XA[J]) {
				I = J - 1;// looking for the nearest I=min, J=max values
				ok = true;
				break;
			}
		}
		if (!ok) {
			J = NX - 1;
			I = J - 1;
		}

		if (XA[I] <= 0.0) {
			XLOGL = false;
		}
		if (!XLOGL) {
			XI = XA[I];
			XJ = XA[J];
			XV = X;
		} else {
			XI = Math.log(XA[I]);
			XJ = Math.log(XA[J]);
			XV = Math.log(X);
		}
		if ((YLOG) && ((YA[I] == 0.0) || (YA[J] == 0.0))) {
			AINTP = 0.0;
		} else {
			if (YLOG) {
				YI = Math.log(YA[I]);
				YJ = Math.log(YA[J]);
				if (XJ == XI) {
					AINTP = YI;
				} else {
					AINTP = (YI * (XJ - XV) + YJ * (XV - XI)) / (XJ - XI);
				}
				AINTP = Math.exp(AINTP);
			} else {
				YI = YA[I];
				YJ = YA[J];
				if (XJ == XI) {
					AINTP = YI;
				} else {
					AINTP = (YI * (XJ - XV) + YJ * (XV - XI)) / (XJ - XI);
				}
			}
		}
		// System.out.println("AINTP  "+AINTP);
		return AINTP;
	}

	// "     FUNCTION TO CALL IF LINEAR DISTRIBUTION FUNCTION IS DESIRED      "
	/**
	 * FUNCTION TO CALL IF LINEAR DISTRIBUTION FUNCTION IS DESIRED. See PEGS4A init_pegs4
	 * @param X X
	 * @return the result
	 */
	public static double ALIN(double X) {
		double ALIN = X;
		return ALIN;
	}

	// "INVERSE OF ALIN"
	/**
	 * INVERSE OF ALIN. See PEGS4A init_pegs4
	 * @param X X
	 * @return the result
	 */
	public static double ALINI(double X) {
		double ALINI = X;
		return ALINI;
	}

	/**
	 * Handle Rayleigh scattering related data.
	 */
	public static void Rayleigh() {
		double AX = 0.;
		double BX = 0.;
		if (IRAYL == 1) {
			// "PWLF FOR DISTRIBUTION FUNCTION FOR RAYLEIGH"
			// "NOTE THAT TOTAL X-SECTION IS FIT WITH REST OF PHOTON CROSS SECTIONS"
			for (int I = 0; I < 100; I++) {
				AFAC2[I] = 0.0; // "ZERO FOR SUMMATION"
				for (int IN = 0; IN < NE; IN++) {
					Double d = new Double(Z[IN]);
					int IZZ = d.intValue() - 1; // "INTEGER Z VALUE"
					AFAC2[I] = AFAC2[I] + PZ[IN] * AFAC[I][IZZ] * AFAC[I][IZZ];
				}
			}// "END OF LOOP ON I"
				// "THE ABOVE CALCULATES THE COHERENT SCATTERING FORM FACTOR FOR THE"
				// "MOLECULE ASSUMING EACH ATOM ACTS INDEPENDENTLY.  THIS ASSUMPTION"
				// "IS ALSO MADE WHEN CALCULATING THE TOTAL COHERENT CROSS SECTION."
				// "THE ASSUMPTION IS KNOWN TO BE VERY WRONG IN SOME SITUATIONS - MOST"
				// "NOTABLY WATER - E.G., SEE L.R. MORIN, J. PHYS. CHEM REF. DATA VOL 11"
				// "(1982) P1091 AND JOHNS AND YAFFE, MEDICAL PHYSICS, VOL. 10 (1983) P40."
			for (int I = 0; I < 100; I++) {
				XVAL[I] = XVAL[I] * XVAL[I];
			}
			AFFI[0] = 0.0;
			for (int I = 1; I < 87; I++) {
				AX = XVAL[I - 1];
				BX = XVAL[I];
				funcname = "AFFACT";
				AFFI[I] = QD(funcname, AX, BX, "AFFACT");// "INTEGRATE AFAC2 FROM AX TO BX"
			}
			for (int I = 1; I < 87; I++) {
				AFFI[I] = AFFI[I] + AFFI[I - 1];
			}
			for (int I = 0; I < 87; I++) {
				AFFI[I] = AFFI[I] / AFFI[86];
			}
			// CALL PWLF1(NGR,NALR,0.0,1.0,0.0,EPR,ZTHRR,ZEPR,NIPR,ALIN,ALINI,
			// AXR,BXR,$MXRL,1,AFR,BFR,RFUNS);
		} // "END OF RAYLEIGH SCATTERING PWLF"
	}

	/**
	 * Internally used.
	 * @param X X
	 * @return the result
	 */
	public static double AFFACT(double X) {
		// construct the doubles!!
		// double [] xa = new double[100];
		// double [] ya = new double[100];
		// for (int i=0; i<100;i++)
		// {
		// xa[i]=XVAL[i];
		// ya[i]=AFAC2[i];
		// }

		double AFFACT = AINTP(X, XVAL, 100, AFAC2, 1, true, true);
		return AFFACT;
	}

	// "*****SUBROUTINE TO COMPUTE RAYLEIGH FUNCTIONS TO BE FIT           "
	// "     IN A WAY THAT AVOIDS REPETITION.                             "
	/**
	 * COMPUTE RAYLEIGH FUNCTIONS TO BE FIT IN A WAY THAT AVOIDS REPETITION. 
	 * @param E E
	 * @param V V
	 */
	public static void RFUNS(double E, double[] V) {
		V[0] = AINTP(E, AFFI, 87, XVAL, 1, true, true);
		return;
	}

	// "*****SUBROUTINE TO COMPUTE PHOTON FUNCTIONS TO BE FIT             "
	// "     IN A WAY THAT AVOIDS REPETITION.                             "
	/**
	 * COMPUTE PHOTON FUNCTIONS TO BE FIT IN A WAY THAT AVOIDS REPETITION.  
	 * @param E E
	 * @param V V
	 */
	public static void GFUNS(double E, double[] V) {
		double PAIR = PAIRTU(E);
		double COMP = COMPTM(E);
		double PHOT = PHOTTE(E);
		double COHR = COHETM(E);
		double TSANSC = PAIR + COMP + PHOT;
		double GMFP = 1.0 / TSANSC;

		// "V(1)=GMFP (GAMMA MEAN FREE PATH)                                  "
		V[0] = GMFP;// TSANSC;//GMFP;
		// "V(2)=GBR1 (GAMMA BRANCHING RATIO NUMBER 1)                        "
		V[1] = PAIR * GMFP;// PAIR/TSANSC;//PAIR*GMFP;
		// "V(3)=GBR2 (GAMMA BRANCHING RATIO NUMBER 2)                        "
		V[2] = (PAIR + COMP) * GMFP;// (PAIR+COMP)/TSANSC;//(PAIR+COMP)*GMFP;
		// "V(4)=CRATIO (COHERENT RATIO)                                      "
		V[3] = TSANSC / (TSANSC + COHR);
		return;
	}

	// Pair production total cross section actually used. Same as PAIRTE for
	// primary energy lese than 50 MeV;otherwise, same as PAIRTM.
	/**
	 * Pair production total cross section actually used. Same as PAIRTE for 
	 * primary energy less than 50 MeV; otherwise, same as PAIRTM.
	 * @param K K
	 * @return the result
	 */
	public static double PAIRTU(double K) {
		// "***THIS IS THE TOTAL PAIR ROUTINE ACTUALLY USED FOR RUNNING.        "
		// "   IT IS THE EMPIRICAL ROUTINE BELOW 50MEV AND THE THEORETICAL ABOVE"
		double PAIRTU = 0.;
		if (K < 50.0) {
			PAIRTU = PAIRTE(K);
		} else {
			PAIRTU = PAIRTM(K);
		}

		return PAIRTU;
	}

	// Empirical total pair production production cross section for a
	// mixture (=sum(Pz(I)*PAIRTZ(z(I))).
	/**
	 * Empirical total pair production production cross section for a mixture.
	 * @param K K
	 * @return the result
	 */
	public static double PAIRTE(double K) {
		// "***'E' FOR EMPIRICAL OR EXPERIMENTAL PAIR PROD CROSS-SECTION.       "
		// "   THIS FUNCTION GIVE THE PROPER MIX OF PAIRTZ'S.                   "
		double PAIRTE = 0.;
		if (K <= 2.0 * RM)
			return PAIRTE;
		// DO I=1,NE [PAIRTE=PAIRTE+PZ(I)*PAIRTZ(Z(I),K);]
		for (int I = 0; I < NE; I++) {
			PAIRTE = PAIRTE + PZ[I] * PAIRTZ(Z[I], K);
		}

		return PAIRTE;
	}

	// Computes contribution to empirical pair production total cross section
	// for an element assuming one atom per molecule. It is obtained by
	// log-linear
	// interpolation of Israel-Storm data.
	/**
	 * Computes contribution to empirical pair production total cross section 
	 * for an element assuming one atom per molecule. It is obtained by log-linear 
	 * interpolation of Israel-Storm data.
	 * @param Z Z
	 * @param K K
	 * @return the result
	 */
	public static double PAIRTZ(double Z, double K) {
		// "***FUNCTION TO GIVE INTERPOLATED TOTAL PAIR PROD CROSS-SECTION      "
		// "   FROM TABULATED DATA USING LOG-LOG INTERPOLATION.                 "
		// "   THIS CROSS-SECTION WILL BE THE CONTRIBUTION THIS Z ELEMENT       "
		// "   WOULD MAKE IF THERE WERE ONE/MOLECULE,I.E., TO GET THE ACTUAL    "
		// "   CONTRIBUTION ONE MUST MULTIPLY BY PZ(I).                         "
		double PAIRTZ = 0.;
		if (K <= RMT2) {
			PAIRTZ = 0.0;
		} else {
			// " PCON CONVERTS CROSS SECTION IN BARNS TO CROSS SECTION IN RLC**-1 "
			double PCON = 1.E-24 * (AN * RHO / WM) * RLC;
			Double d = new Double(Z);
			int IZ = d.intValue() - 1;// 0 biased

			// construct the doubles!!
			double[] xa = new double[17];
			double[] ya = new double[17];// mxpaire
			for (int i = 0; i < 17; i++) {
				xa[i] = PRE[i];
				ya[i] = PRD[i][IZ];
			}
			PAIRTZ = PCON * AINTP(K, PRE, 17, ya, 1, true, false);
		}

		return PAIRTZ;
	}

	// Pair production total cross section for a mixture of elemente, obtained
	// by numerical integration of differential cross section.
	/**
	 * Pair production total cross section for a mixture of elements, obtained 
	 * by numerical integration of differential cross section.
	 * @param K0 K0
	 * @return the result
	 */
	public static double PAIRTM(double K0) {
		double PAIRTM = 0.;
		if (K0 <= 2. * RM) {
			PAIRTM = 0.0;
		} else {
			PAIRTM = PAIRRM(K0, RM, K0 - RM);
		}

		return PAIRTM;
	}

	// Pair production cross section, integrated over some energy range, for a
	// mixture of elements.
	/**
	 * Pair production cross section, integrated over some energy range, for a 
	 * mixture of elements.
	 * @param K K
	 * @param E1 E1
	 * @param E2 E2
	 * @return the result
	 */
	public static double PAIRRM(double K, double E1, double E2) {
		double PAIRRM = 0.;
		// DO I=1,NE [PAIRRM=PAIRRM+PZ(I)*PAIRRZ(Z(I),K,E1,E2);]
		for (int I = 0; I < NE; I++) {
			PAIRRM = PAIRRM + PZ[I] * PAIRRZ(Z[I], K, E1, E2);
		}
		return PAIRRM;
	}

	// Pair production cross section, integrated over some energy range, for an
	// element.
	/**
	 * Pair production cross section, integrated over some energy range, for an element.
	 * @param Z Z
	 * @param K K
	 * @param E1 E1
	 * @param E2 E2
	 * @return the result
	 */
	public static double PAIRRZ(double Z, double K, double E1, double E2) {
		@SuppressWarnings("unused")
		double DUMMY = PAIRDZ(Z, K, E1);
		funcname = "PAIRFZ";
		double PAIRRZ = QD(funcname, E1, E2, "PAIRFZ");
		return PAIRRZ;
	}

	// Pair production differential cross section for an element.
	/**
	 * Pair production differential cross section for an element.
	 * @param Z Z
	 * @param KA KA
	 * @param E E
	 * @return the result
	 */
	public static double PAIRDZ(double Z, double KA, double E) {
		// "***ALL ENTRIES TO THIS FUNCTION GIVE THE CONTRIBUTION THAT ELEMENT Z"
		// "   WOULD HAVE IF THERE WERE ONE PER MOLECULE.                      "
		// "   ENTRIES STARTING WITH D DO THEIR OWN INITIALIZATION.             "
		// "   ENTRIES STARTING WITH F RELY ON PREVIOUS D FOR INITIALIZATION.   "
		// eq.2.7.8 SLAC 265
		K = KA;
		DELC = 136. * Math.pow(Z, -1. / 3.) * RM / K;
		CONST = (AN * RHO / WM) * R0 * R0 * FSC * Z * (Z + XSIF(Z)) * RLC
				/ (K * K * K);
		XLNZ = 4. / 3. * Math.log(Z);
		if (K >= 50.0)
			XLNZ = XLNZ + 4. * FCOULC(Z);
		// ".....DELTAM IS THE DELTA AT WHICH THE SQUARE BRACKETS GO TO ZERO      "
		DELTAM = Math.exp((21.12 - XLNZ) / 4.184) - 0.952;// eq. 2.7.31 SLAC-265
		double PAIRDZ = PAIRFZ(E);
		return PAIRDZ;
	}

	// One argument form of PAIRDZ.
	/**
	 * Called by PAIRDZ.
	 * @param E E
	 * @return the result
	 */
	public static double PAIRFZ(double E) {
		double EPS = E / K;// eq 2.7.95 SLAC 265
		double ONEEPS = 1. - EPS;
		double PAIRFZ = 0.;
		double SB1 = 0.;
		double SB2 = 0.;

		if (ONEEPS == 0.0) {
			ONEEPS = 1.18E-38;
		}
		double DELTA = DELC / (EPS * ONEEPS);// eq 2.7.96 SLAC 265
		if (DELTA >= DELTAM) {
			PAIRFZ = 0.0;
		} else {
			if (DELTA <= 1.) {
				SB1 = 20.867 + DELTA * (-3.242 + DELTA * 0.625) - XLNZ;// eq
																		// 2.7.14
																		// SLAC
																		// 265
				SB2 = 20.209 + DELTA * (-1.930 + DELTA * (-0.086)) - XLNZ;// eq
																			// 2.7.15
																			// SLAC
																			// 265
			} else {
				SB1 = 21.12 - 4.184 * Math.log(DELTA + 0.952) - XLNZ;// eq
																		// 2.7.14
																		// SLAC
																		// 265
				SB2 = SB1;// eq 2.7.15 SLAC 265
			}
			double EPLUS = K - E;
			PAIRFZ = CONST
					* ((E * E + EPLUS * EPLUS) * SB1 + 0.666667 * E * EPLUS
							* SB2);
		}
		return PAIRFZ;
	}

	/**
	 * TOTAL CROSS SECTION FOR COMPTON SCATTERING WITH INCIDENT PHOTON ENERGY K0.
	 * @param K0 K0
	 * @return the result
	 */
	public static double COMPTM(double K0) {
		// "***TOTAL CROSS SECTION FOR COMPTON SCATTERING WITH INCIDENT PHOTON"
		// "   ENERGY K0."
		double COMPTM = 0.;
		double K1 = K0 * RM / (RM + 2. * K0);
		COMPTM = COMPRM(K0, K1, K0); // "EQ. 2.9.7"

		return COMPTM;
	}

	/**
	 * COMPTON CROSS SECTION FOR INCIDENT PHOTON OF ENERGY K0 TO SCATTER 
	 * INTO THE ENERGY RANGE K1 TO K2.
	 * @param K0 K0
	 * @param K1 K1
	 * @param K2 K2
	 * @return the result
	 */
	public static double COMPRM(double K0, double K1, double K2) {
		// "***COMPTON CROSS SECTION FOR INCIDENT PHOTON OF ENERGY K0 TO SCATTER"
		// "   INTO THE ENERGY RANGE K1 TO K2."
		double K0P = K0 / RM;// "FOR VERY LOW ENERGY PHOTONS"
		double CCOMP2 = RLC * EDEN * PI * R0 * R0 / K0P; // "CONSTANT FACTOR IN EQ.2.9.2"
		double C1 = 1. / (K0P * K0P);
		double C2 = 1. - (2. + 2. * K0P) / (K0P * K0P);
		double C3 = (1. + 2. * K0P) / (K0P * K0P);
		double EPS1 = K1 / K0;
		double EPS2 = K2 / K0;
		double COMPRM = CCOMP2
				* (C1 * (1. / EPS1 - 1. / EPS2) + C2 * Math.log(EPS2 / EPS1)
						+ EPS2 * (C3 + 0.5 * EPS2) - EPS1 * (C3 + 0.5 * EPS1)); // "EQ. 2.9.2"
		return COMPRM;
	}

	/**
	 * CROSS SECTION IN 1/RLC FOR COHERENT SCATTERING FOR A MIXTURE.
	 * @param K K
	 * @return the result
	 */
	public static double COHETM(double K) {
		double COHETM = 0.;
		// DO I=1,NE [
		for (int I = 0; I < NE; I++) {
			COHETM = COHETM + PZ[I] * COHETZ(Z[I], K);
		}

		return COHETM;
	}

	/**
	 * CONVERTS STORM AND ISRAEL CROSS SECTION IN BARNS/ATOM 
	 * INTO 1/RLC FOR COHERENT SCATTERING
	 * @param Z Z
	 * @param K K
	 * @return the result
	 */
	public static double COHETZ(double Z, double K) {
		// "PCON CONVERTS STORM AND ISRAEL CROSS SECTION IN BARNS/ATOM"
		// "INTO 1/RLC FOR COHERENT SCATTERING "
		double PCON = 1.E-24 * (AN * RHO / WM) * RLC;
		Double d = new Double(Z);
		int IZ = d.intValue() - 1;// 0 biased
		// construct the doubles!!
		double[] xa = new double[NPHE[IZ]];
		double[] ya = new double[NPHE[IZ]];
		for (int i = 0; i < NPHE[IZ]; i++) {
			xa[i] = PHE[i][IZ];
			ya[i] = COHE[i][IZ];
		}
		double COHETZ = PCON * AINTP(K, xa, NPHE[IZ], ya, 1, true, true);
		return COHETZ;
	}

	// "LOG OF KINETIC ENERGY"
	/**
	 * LOG OF KINETIC ENERGY
	 * @param E E
	 * @return the result
	 */
	public static double ALKE(double E) {
		double ALKE = Math.log(E - RM);
		return ALKE;
	}

	// "INVERSE OF LOG OF KINETIC ENERGY"
	/**
	 * INVERSE OF LOG OF KINETIC ENERGY
	 * @param X X
	 * @return the result
	 */
	public static double ALKEI(double X) {
		double ALKEI = Math.exp(X) + RM;// CORRECT!!ALKEI(ALKE(E))=E
		return ALKEI;
	}

	// $REAL E,V(8);
	// "*****SUBROUTINE TO COMPUTE ELECTRON FUNCTIONS TO BE FIT           "
	// " IN A WAY THAT AVOIDS REPETITION.
	/**
	 * COMPUTE ELECTRON FUNCTIONS TO BE FIT IN A WAY THAT AVOIDS REPETITION.
	 * @param E E
	 * @param V V
	 */
	public static void EFUNS(double E, double[] V) {
		double BREM = 0.;
		double AMOLL = 0.;
		double BHAB = 0.;
		double ANNIH = 0.;
		double ESIG = 0.;
		double PSIG = 0.;
		if (((IUNRST == 0) || (IUNRST == 1)) || (IUNRST == 5)) {
			// "   REGULAR DATA SET OR"
			// "   UNRESTRICTED COLLISIONAL(1) OR RADIATIVE(5) STOPPING POWERS"
			BREM = BREMTM(E);
			AMOLL = AMOLTM(E);
			BHAB = BHABTM(E);
			ANNIH = ANIHTM(E);
			ESIG = BREM + AMOLL; // "TOTAL ELECTRON CROSS-SECTION"
			V[0] = ESIG;// instead of 1...etc
			PSIG = BREM + BHAB + ANNIH;// "TOTAL POSITRON CROSS SECTION"
			V[1] = PSIG;
			V[2] = SPTOTE(E, AE, AP);// "TOTAL ELECTRON STOPPING POWER"
			V[3] = SPTOTP(E, AE, AP);// "TOTAL POSITRON STOPPING POWER"
										// "EBR1=BREM/(BREM+AMOLL)"
			if (ESIG > 0.0) {
				V[4] = BREM / ESIG;
			} else { // "BELOW THRESHOLD FOR BOTH BREMS AND MOLLER. USE THE BRANCHING"
						// "RATIO THAT EXISTED WHEN CROSS SECTION APPROACHED ZERO"
				if (THBREM <= THMOLL) {
					V[4] = 1.0;
				} else {
					V[4] = 0.0;
				}
			}
			V[5] = BREM / PSIG;// "PBR1=BREM/(BREM+BHABA+ANNIH)"
			V[6] = (BREM + BHAB) / PSIG;// "PBR2=(BREM+BHABA)/(PSIG)"
			// "MAXIMUM ALLOWED TRANSPORT STEP, FROM MULTIPLE SCATTERING"
			V[7] = TMXS(E);
		}
		// "  THE FOLLOWING ARE UNDOCUMENTED ADDITIONS"
		else if (IUNRST == 2) {// "FULL CSDA DATA SET WITH NO DISCRETE INTERACTIONS"
								// CSDA=continuous slowing down approximation
			V[0] = 0.0;
			V[1] = 0.0;
			V[4] = 0.0;
			V[5] = 0.0;
			V[6] = 0.0;
			// "ZERO TOTAL CROSS SECTION FOR EL & POS, AND ZERO BRANCHING RATIOS"
			V[2] = SPTOTE(E, E, E);// "  TOTAL UNRESTRICTED STOPPING POWER"
			V[3] = SPTOTP(E, E, E);
			V[7] = TMXS(E);
		} else if (IUNRST == 3) {// "CONSIDER BREM AND ANNIHILATION IN FLIGHT AS"
									// DISCRETE EVENTS BUT TREAT DELTAS IN CSDA"
			BREM = BREMTM(E);
			ANNIH = ANIHTM(E);
			V[0] = BREM;// "TOTAL X-SECTION IS JUST BREM"
			V[1] = BREM + ANNIH;// "POSITRONS ALSO HAVE ANNIHILATION IN FLIGHT"
			V[2] = SPTOTE(E, E, AP);// "UNRESTRICTED COLLISIONAL+RESTRICTED RADIATIVE"
			V[3] = SPTOTP(E, E, AP);// "    ''                             ''    "
			V[4] = 1.0;// "ALL ELECTRON EVENTS ARE BREM EVENTS"
			V[5] = BREM / V[1];// "FRACTION FOR POSITRONS WHICH IS BREM"
			V[6] = V[5];// "FRACTION WHICH IS BREM + COLLISION(=0)"
			V[7] = TMXS(E);
		} else if (IUNRST == 4) {// "CREATE SECONDARIES BUT HAVE NO DISCRETE BREM OR"
									// "ANNIHILATION IN FLIGHT"
			V[0] = AMOLTM(E);// "ONLY MOLLERS FOR ELECTRONS"
			V[1] = BHABTM(E);// "ONLY BHABHA FOR POSITRONS"
			V[2] = SPTOTE(E, AE, E);// "RESTRICTED COLLISIONAL + UNRESTRICTED RADIATIVE"
			V[3] = SPTOTP(E, AE, E);// "         ''                         ''       "
			V[4] = 0.0;// "I.E. NEVER BREMS"
			V[5] = 0.0;// "I.E. NEVER BREMS"
			V[6] = 1.0;// "ALL BHABHA - NO ANNIHILATION"
			V[7] = TMXS(E);
		} else {// "IUNRST=6 OR 7 NOT ALLOWED HERE"

		}
		return;
	}

	// "MAXIMUM STEP SIZE VALID FOR MULTIPLE SCATTERING"
	/**
	 * MAXIMUM STEP SIZE VALID FOR MULTIPLE SCATTERING
	 * @param E E
	 * @return the result
	 */
	public static double TMXS(double E) {
		// double SAFETY=0.8;
		// double TABSMX=10.0;
		double TMXS = 0.;
		TMXS = Math.min(TMXB(E) * SAFETY, TABSMX);
		// "THE FACTORE 'SAFETY' IS TO KEEP SOMEWHAT BELOW BETHE'S LIMIT"
		// "TABSMX IS AN ABSOLUTE LIMIT TO SIZE OF ELECTRON TRANSPORT,"
		// "  INDEPENDENT OF THE MULTIPLE SCATTERING LIMIT"
		return TMXS;
	}

	/**
	 * FINDS THE TRANSPORT DISTANCE WHICH AT THIS ENERGY IS THE LARGEST CONSISTENT WITH BETHE'S CRITERION, NAMELY 
	 * XC**2*B.LE.1;   SINCE XC ANB B ARE INCREASING FUNCTIONS OF THE TRANSPORT DISTANCE, THE CRITERION FOR TMXB IS THEN 
	 * XC**2*B=1;
	 * @param E E
	 * @return the result
	 */
	public static double TMXB(double E) {
		double TMXB = 0.;
		// "THIS FUNCTION FINDS THE TRANSPORT DISTANCE WHICH AT THIS ENERGY"
		// "IS THE LARGEST CONSISTENT WITH BETHE'S CRITERION, NAMELY"
		// "  XC**2*B.LE.1;   SINCE XC ANB B ARE INCREASING FUNCTIONS OF T"eq
		// 2.14.53 SLAC265
		// " THE TRANSPORT DISTANCE, THE CRITERION FOR TMXB IS THEN"
		// "  XC**2*B=1;      OTHER RELATIONS USED IN THE DERIVATION ARE"
		// "EXPLAINED IN SUBROUTINES MIX AND MOLIER.  THEY ARE:"
		// "  XC=XCC*SQRT(T)/(E*BETA**2); "
		// "  EXP(B)/B = BLCC*T/BETA**2;  "
		// "FROM THESE IS DERIVED THE EQUATION THIS FUNCTION IS BASED ON:"
		// "TMXB=(E**2*BETA**2/XCC**2)*BETA**2/ALOG(BLCC*(E**2*BETA**2/XCC**2));"
		double ESQ = E * E;
		double BETA2 = 1.0 - RMSQ / ESQ;
		double PX2 = ESQ * BETA2 / (XCC * XCC);
		TMXB = PX2 * BETA2 / Math.log(BLCC * PX2);

		return TMXB;
	}

	// Calculates the total stopping power(ionization+soft brems) for
	// positrons(+) for specified cutoffs.
	/**
	 * Calculates the total stopping power(ionization+soft brems) for positrons(+) for specified cutoffs.
	 * @param E0 E0
	 * @param EE EE
	 * @param EG EG
	 * @return the result
	 */
	public static double SPTOTP(double E0, double EE, double EG) {
		double SPTOTP = 0.;
		if (IUNRST == 0) {// "RESTRICTED TOTAL STOPPING POWER I.E. NORMAL"
			SPTOTP = SPIONP(E0, EE) + BRMSTM(E0, EG);
		} else if (IUNRST == 1) {// "UNRESTRICTED COLLISION"
			SPTOTP = SPIONP(E0, E0);
		} else if (IUNRST == 2) {// "UNRESTRICTED COLLISION AND RADIATIVE        "
			SPTOTP = SPIONP(E0, E0) + BRMSTM(E0, E0);
		} else if (IUNRST == 3) {// "UNRESTRICTED COLLISION +RESTRICTED RADIATIVE"
			SPTOTP = SPIONP(E0, E0) + BRMSTM(E0, EG);
		} else if (IUNRST == 4) {// "RESTRICTED COLLISION +UNRESTRICTED RADIATIVE"
			SPTOTP = SPIONP(E0, EE) + BRMSTM(E0, E0);
		} else if (IUNRST == 5) {// "UNRESTRICTED RADIATIVE"
			SPTOTP = BRMSTM(E0, E0);
		} else if (IUNRST == 6) {// "RESTRICTED RADIATIVE  "
			SPTOTP = BRMSTM(E0, EG);
		} else if (IUNRST == 7) {// "RESTRICTED COLLISON   "
			SPTOTP = SPIONP(E0, EE);
		}
		return SPTOTP;

	}

	// Calculates the stopping power due to ionizations for positrons (+)
	/**
	 * Calculates the stopping power due to ionizations for positrons (+)
	 * @param E0 E0
	 * @param EE EE
	 * @return the result
	 */
	public static double SPIONP(double E0, double EE) {
		double SPIONP = SPIONB(E0, EE, true);
		return SPIONP;
	}

	// Calculates the total stopping power(ionization+soft brems) for
	// electrons(-) for specified cutoffs.
	/**
	 * Calculates the total stopping power(ionization+soft brems) for electrons(-) for specified cutoffs.
	 * @param E0 E0
	 * @param EE EE
	 * @param EG EG
	 * @return the result
	 */
	public static double SPTOTE(double E0, double EE, double EG) {
		double SPTOTE = 0.;
		if (IUNRST == 0) {// "RESTRICTED TOTAL STOPPING POWER I.E. NORMAL"
			SPTOTE = SPIONE(E0, EE) + BRMSTM(E0, EG);
		} else if (IUNRST == 1) {// "UNRESTRICTED COLLISION"
			SPTOTE = SPIONE(E0, E0);
		} else if (IUNRST == 2) {// "UNRESTRICTED COLLISION AND RADIATIVE        "
			SPTOTE = SPIONE(E0, E0) + BRMSTM(E0, E0);
		} else if (IUNRST == 3) {// "UNRESTRICTED COLLISION +RESTRICTED RADIATIVE"
			SPTOTE = SPIONE(E0, E0) + BRMSTM(E0, EG);
		} else if (IUNRST == 4) {// "RESTRICTED COLLISION +UNRESTRICTED RADIATIVE"
			SPTOTE = SPIONE(E0, EE) + BRMSTM(E0, E0);
		} else if (IUNRST == 5) {// "UNRESTRICTED RADIATIVE"
			SPTOTE = BRMSTM(E0, E0);
		} else if (IUNRST == 6) {// "RESTRICTED RADIATIVE  "
			SPTOTE = BRMSTM(E0, EG);
		} else if (IUNRST == 7) {// "RESTRICTED COLLISON   "
			SPTOTE = SPIONE(E0, EE);
		}

		return SPTOTE;
	}

	// Calculates the stopping power due to ionizations for electrons (-)
	/**
	 * Calculates the stopping power due to ionizations for electrons (-)
	 * @param E0 E0
	 * @param EE EE
	 * @return the result
	 */
	public static double SPIONE(double E0, double EE) {
		double SPIONE = 0.;
		SPIONE = SPIONB(E0, EE, false);
		return SPIONE;

	}

	// Does the work for SPIONE and SPIONP One argument tells whether to compute
	// stopping power for electron or positron.
	/**
	 * Called by SPIONE and SPIONP. One argument tells whether to compute 
	 * stopping power for electron or positron.
	 * @param E0 E0
	 * @param EE EE
	 * @param POSITR POSITR
	 * @return the result
	 */
	public static double SPIONB(double E0, double EE, boolean POSITR) {
		// evaluates eq.2.13.5 SLAC report
		// "***STOPPING POWER FOR AN ELECTRON.  THIS FUNCTION ALSO HAS OTHER    "
		// "   ENTRY POINTS. . .                                                "
		// "   SPIONP(E0,EE) - STOPPING POWER FOR A POSITRON                    "
		// "   SPINIT(MEDIUM) - DOES INITIALIZATION FOR THE OTHER ENTRY POINTS  "
		// "   WHENEVER THE MEDIUM CHANGES.                                     "
		// "   THIS FUNCTION IS FOR STOPPING POWER DUE TO COLLISIONS WITH LESS  "
		// "   THAN EE-RM ENERGY TRANSFER AND DOES NOT INCLUDE SOFT BREMS LOSS. "
		// "   WE USE BERGER AND SELTZER'S FORMULATION.                         "
		// "   STOPPING POWER IS RETURNED IN UNITS OF MEV/R.L.                  "
		double SPIONB = 0.;
		double G = E0 / RM;// gamma
		double EEM = EE / RM - 1.;
		// "     T IS BERGER'S TAU                                                "
		double T = G - 1;// eq 2.13.10 SLAC265
		double ETA2 = T * (G + 1.);
		double BETA2 = ETA2 / (G * G);
		double ALETA2 = Math.log(ETA2);
		double X = 0.21715 * ALETA2;// 2.13.20 SLAC 265
		double D = 0.;
		double FTERM = 0.;
		double TP2 = 0.;
		double D2 = 0.;
		double D3 = 0.;
		double D4 = 0.;
		double DELTA = 0.;
		// "     0.21715=ALOG10(E)/2.   THIS FACTOR IS BECAUSE THE DEFINITION OF  "
		// "     X IS ALOG10(P/(MC)) AND ETA2=ETA**2=(P/MC)**2   OK!!                 "
		if (!POSITR) {// "THIS IS ELECTRON CASE"
						// "     COMPUTE F-TERM FOR ELECTRON.  MAXIMUM TRANSFER IS T/2            "
						// "     D IS BERGER'S CAPITOL DELTA.                                     "
			D = Math.min(EEM, 0.5 * T);// eq 2.13.13 SLAC265+2.13.12
			// "     EEM IS DEFINED AS EE/RM-1 IS ENERGY TRANSFER THRESHOLD FOR       "
			// "     DISCRETE MOLLER AND BHABHA SCATTERING(IN UNITS OF RM.)           "
			FTERM = -1. - BETA2 + Math.log((T - D) * D) + T / (T - D)
					+ (D * D / 2. + (2. * T + 1.) * Math.log(1. - D / T))
					/ (G * G);// eq 2.13.16 SLAC265
			// "     COMPUTE F-TERM FOR POSITRON.  MAXIMUM TRANSFER IS T.             "
		} else {// "THIS IS POSITRON CASE"
			D = Math.min(EEM, T);// eq 2.13.13 SLAC265+2.13.12
			TP2 = T + 2.;// =1/y->//eq 2.13.11 SLAC265
			D2 = D * D;
			D3 = D * D2;
			D4 = D * D3;
			FTERM = Math.log(T * D)
					- (BETA2 / T)
					* (T + 2. * D - (3. * D2 / 2.) / TP2 - (D - D3 / 3.)
							/ (TP2 * TP2) - (D2 / 2. - T * D3 / 3. + D4 / 4.)
							/ (TP2 * TP2 * TP2));// eq 2.13.17 SLAC265
		}

		// "NOW COMPUTE THE DENSITY CORRECTION TERM."
		if (EPSTFL == 0) {// "USE STANDARD PEGS4 METHOD"//eq 2.13.19 SLAC265
			if (X <= X0) {
				DELTA = 0.0;
			} else if (X < X1) {
				DELTA = TOLN10 * X - CBAR + AFACT * Math.pow((X1 - X), SK);
			} else {
				DELTA = TOLN10 * X - CBAR;
			}
		} else {// "USE LINEAR INTERPOLATION OF USER SUPPLIED INPUT TABLE"
				// "IEPST IS A POINTER SUCH THAT EPSTEN(IEPST) <= E0 < EPSTEN(IEPST+1)              "
				// "IEPST IS INITIALIZED IN BLOCK DATA TO 1. WE START FROM    "
				// "THE PREVIOUS VALUE OF THE POINTER SINCE WE ASSUME THAT    "
				// "THE CODE IS WORKING UP OR DOWN A GRID.                    "
				// "THIS CODING IS FAR FROM OPTIMAL                         "
			if (E0 >= EPSTEN[IEPST - 1])// IEPST is 1,2,...NEPST!!!
			{// "AT OR ABOVE PREVIOUS ENTRY"
				if (E0 == EPSTEN[IEPST - 1]) {// "FOUND ENTRY, INCLUDING THE POSSIBILITY"
												// "THAT WE ARE AT THE TOP OF THE TABLE"
												// GO TO :END-SEARCH:
												// "NOW JUST INTERPOLATE LINEARLY IN THE ENERGY"
					if (IEPST < NEPST) {
						DELTA = EPSTD[IEPST - 1] + (E0 - EPSTEN[IEPST - 1])
								/ (EPSTEN[IEPST] - EPSTEN[IEPST - 1])
								* (EPSTD[IEPST] - EPSTD[IEPST - 1]);
					} else {
						DELTA = EPSTD[NEPST - 1];
					}
					SPIONB = (SPC1 / BETA2)
							* (Math.log(T + 2.) - SPC2 + FTERM - DELTA);// ok!!2.13.5
																		// SLAC265
					return SPIONB;
				}
				for (int I = IEPST - 1; I < NEPST - 1; I++) {
					if (E0 < EPSTEN[I + 1]) {// "WE FOUND IT"
						IEPST = I + 1;
						// GO TO :END-SEARCH:
						// "NOW JUST INTERPOLATE LINEARLY IN THE ENERGY"
						if (IEPST < NEPST) {
							DELTA = EPSTD[IEPST - 1] + (E0 - EPSTEN[IEPST - 1])
									/ (EPSTEN[IEPST] - EPSTEN[IEPST - 1])
									* (EPSTD[IEPST] - EPSTD[IEPST - 1]);
						} else {
							DELTA = EPSTD[NEPST - 1];
						}
						SPIONB = (SPC1 / BETA2)
								* (Math.log(T + 2.) - SPC2 + FTERM - DELTA);// ok!!2.13.5
																			// SLAC265
						return SPIONB;
					}
				}
				// "IF WE FALL THRU TO HERE, WE MUST BE AT UPPER ENERGY"
				IEPST = NEPST;
				// GO TO :END-SEARCH:;
				// "NOW JUST INTERPOLATE LINEARLY IN THE ENERGY"
				if (IEPST < NEPST) {
					DELTA = EPSTD[IEPST - 1] + (E0 - EPSTEN[IEPST - 1])
							/ (EPSTEN[IEPST] - EPSTEN[IEPST - 1])
							* (EPSTD[IEPST] - EPSTD[IEPST - 1]);
				} else {
					DELTA = EPSTD[NEPST - 1];
				}
				SPIONB = (SPC1 / BETA2)
						* (Math.log(T + 2.) - SPC2 + FTERM - DELTA);// ok!!2.13.5
																	// SLAC265
				return SPIONB;
			}// "END OF BLOCK E0>EPSTEN(IEPST)"
			else {// "E0<EPSTEN[IEPST-1]"
					// for (int I = IEPST-1;I<2;I--)
				for (int I = IEPST - 1; I >= 2; I--) {
					if (E0 >= EPSTEN[I - 1]) {
						IEPST = I;
						// GO TO :END-SEARCH:;
						// "NOW JUST INTERPOLATE LINEARLY IN THE ENERGY"
						if (IEPST < NEPST) {
							DELTA = EPSTD[IEPST - 1] + (E0 - EPSTEN[IEPST - 1])
									/ (EPSTEN[IEPST] - EPSTEN[IEPST - 1])
									* (EPSTD[IEPST] - EPSTD[IEPST - 1]);
						} else {
							DELTA = EPSTD[NEPST - 1];
						}
						SPIONB = (SPC1 / BETA2)
								* (Math.log(T + 2.) - SPC2 + FTERM - DELTA);// ok!!2.13.5
																			// SLAC265
						return SPIONB;
					}
				}
				// "   IF WE GET HERE WE MUST BE IN THE FIRST REGION"
				IEPST = 1;
				// "NOW JUST INTERPOLATE LINEARLY IN THE ENERGY"
				if (IEPST < NEPST) {
					DELTA = EPSTD[IEPST - 1] + (E0 - EPSTEN[IEPST - 1])
							/ (EPSTEN[IEPST] - EPSTEN[IEPST - 1])
							* (EPSTD[IEPST] - EPSTD[IEPST - 1]);
				} else {
					DELTA = EPSTD[NEPST - 1];
				}
				SPIONB = (SPC1 / BETA2)
						* (Math.log(T + 2.) - SPC2 + FTERM - DELTA);// ok!!2.13.5
																	// SLAC265
				return SPIONB;
			}
		}// /"END OF EPSTFL NON-ZERO BLOCK"
			// "     NOW PUT IT ALL TOGETHER                                          "
		SPIONB = (SPC1 / BETA2) * (Math.log(T + 2.) - SPC2 + FTERM - DELTA);// ok!!2.13.5
																			// SLAC265
		return SPIONB;
	}

	// Soft bremsatrahlung total cross for a mixture of elements.
	/**
	 * Soft bremsstrahlung total cross for a mixture of elements.
	 * @param E0 E0
	 * @param EG EG
	 * @return the result
	 */
	public static double BRMSTM(double E0, double EG) {
		double BRMSTM = 0.;
		double zero = 0.;
		if (E0 <= RM) {
			BRMSTM = 0.;
		} else {
			double AU = Math.min(EG, E0 - RM);
			BRMSTM = BRMSRM(E0, zero, AU);// eq 2.7.112 SLAC!!
		}

		return BRMSTM;
	}

	// Soft bremeetrahlung cross section,integrated over some energy range
	// for a mixture of elements.
	/**
	 * Soft bremsstrahlung cross section,integrated over some energy range 
	 * for a mixture of elements.
	 * @param E E
	 * @param K1 K1
	 * @param K2 K2
	 * @return the result
	 */
	public static double BRMSRM(double E, double K1, double K2) {
		double BRMSRM = 0.;
		// DO I=1,NE[BRMSRM=BRMSRM+PZ(I)*BRMSRZ(Z(I),E,K1,K2);]
		for (int I = 0; I < NE; I++) {
			BRMSRM = BRMSRM + PZ[I] * BRMSRZ(Z[I], E, K1, K2);
		}
		return BRMSRM;
	}

	// Soft bremeetrahlung cross section,integrated over some energy range
	// for a element.
	/**
	 * Soft bremsstrahlung cross section,integrated over some energy range 
	 * for a element.
	 * @param Z Z
	 * @param E E
	 * @param K1 K1
	 * @param K2 K2
	 * @return the result
	 */
	public static double BRMSRZ(double Z, double E, double K1, double K2) {
		double BRMSRZ = 0.;
		@SuppressWarnings("unused")
		double DUMMY = BRMSDZ(Z, E, K1);// initialize
		funcname = "BRMSFZ";
		BRMSRZ = QD(funcname, K1, K2, "BRMSFZ");

		return BRMSRZ;
	}

	// Annihilation total cross section for a mixture of elements.
	/**
	 * Annihilation total cross section for a mixture of elements.
	 * @param E0 E0
	 * @return the result
	 */
	public static double ANIHTM(double E0) {
		// "***TOTAL CROSS SECTION FOR TWO-PHOTON POSITRON-ELECTRON ANNIHILATION"
		// "   WITH INCIDENT POSITRON ENERGY(TOTAL) OF E0."
		double ANIHTM = 0.;
		double GAM = E0 / RM; // "EQ.2.12.3"
		double P0P2 = GAM * GAM - 1.0;
		double P0P = Math.sqrt(P0P2); // "EQ.2.12.6"
		double CANIH = RLC * EDEN * PI * R0 * R0 / (GAM + 1.); // "CONSTANT FACTOR IN EQ.2.12.14"
		ANIHTM = CANIH
				* ((GAM * GAM + 4. * GAM + 1.) / P0P2 * Math.log(GAM + P0P) - (GAM + 3.)
						/ P0P); // "EQ.2.12.14"

		return ANIHTM;
	}

	// Bhabha total cross section for a mixture of elements.
	/**
	 * Bhabha total cross section for a mixture of elements.
	 * @param E0 E0
	 * @return the result
	 */
	public static double BHABTM(double E0) {
		// "***TOTAL CROSS SECTION FOR BHABHA SCATTERING WITH INCIDENT POSITRON"
		// "   ENERGY(TOTAL) OF E0."
		double BHABTM = 0.;
		if (E0 <= AE) {
			BHABTM = 0.;
		} else {
			BHABTM = BHABRM(E0, AE, E0);
		} // "EQ.2.11.4"

		return BHABTM;
	}

	// Bhabha cross section, integrated over some energy range, for a mixture of
	// elements.
	/**
	 * Bhabha cross section, integrated over some energy range, for a mixture of elements.
	 * @param EN0 EN0
	 * @param EN1 EN1
	 * @param EN2 EN2
	 * @return the result
	 */
	public static double BHABRM(double EN0, double EN1, double EN2) {
		// "***BHABHA CROSS SECTION FOR INCIDENT POSITRON OF TOTAL ENERGY EN0 TO"
		// "   PRODUCE SECONDARY ELECTRON IN THE TOTAL ENERGY RANGE EN1 TO EN2."
		double BHABRM = 0.;
		double T0 = EN0 - RM;
		double T1 = EN1 - RM;
		double T2 = EN2 - RM;
		double TM = T0 / RM;
		double EM = TM + 1.;// =gamma
		double Y = 1. / (TM + 2.);
		double BETASI = 1. / (1. - 1. / (EM * EM));// gamma=EM, beta*beta
		double CBHAB2 = RLC * EDEN * 2. * PI * R0 * R0 / TM; // "CONSTANT FACTOR IN EQ.2.11.2"
		double B1 = 2. - Y * Y;
		double B2 = 3. - Y * (6. - Y * (1. - Y * 2.));
		double B3 = 2. - Y * (10. - Y * (16. - Y * 8.));
		double B4 = 1. - Y * (6. - Y * (12. - Y * 8.));
		double EPS1 = T1 / T0;
		double EPS2 = T2 / T0;
		BHABRM = CBHAB2
				* (BETASI * (1. / EPS1 - 1. / EPS2) - B1
						* Math.log(EPS2 / EPS1) + B2 * (EPS2 - EPS1) + EPS2
						* EPS2 * (EPS2 * B4 / 3. - 0.5 * B3) - EPS1 * EPS1
						* (EPS1 * B4 / 3. - 0.5 * B3)); // "EQ.2.11.2"

		return BHABRM;
	}

	// Moller total cross section for a mixture of elements.
	/**
	 * Moller total cross section for a mixture of elements.
	 * @param E0 E0
	 * @return the result
	 */
	public static double AMOLTM(double E0) {
		// "***TOTAL CROSS SECTION FOR MOLLER SCATTERING WITH INCIDENT ELECTRON"
		// "   ENERGY(TOTAL) OF E0."
		double AMOLTM = 0.;
		double T0 = 0.;
		if (E0 <= THMOLL) {
			AMOLTM = 0.;
		} else {
			T0 = E0 - RM;
			AMOLTM = AMOLRM(E0, AE, T0 * 0.5 + RM); // "EQ.2.10.6 SLAC 265
		}
		return AMOLTM;
	}

	// Moller cross section, integrated over some energy range, for a mixture of
	// elements.
	/**
	 * Moller cross section, integrated over some energy range, for a mixture of elements.
	 * @param EN0 EN0
	 * @param EN1 EN1
	 * @param EN2 EN2
	 * @return the result
	 */
	public static double AMOLRM(double EN0, double EN1, double EN2) {
		// "***MOLLER CROSS SECTION FOR INCIDENT ELECTRON OF TOTAL ENERGY EN0 TO"
		// "   PRODUCE SECONDARY ELECTRON IN THE TOTAL ENERGY RANGE EN1 TO EN2."
		double AMOLRM = 0.;
		double T0 = EN0 - RM;
		double T1 = EN1 - RM;
		double T2 = EN2 - RM;
		double TM = T0 / RM;// =gamma-1, gamma=EN0/RM!!!
		double EM = TM + 1.;// gamma
		double C1 = (TM / EM) * (TM / EM);
		double C2 = (2. * TM + 1.) / (EM * EM);
		double BETASQ = 1. - 1. / (EM * EM);
		double CMOLL2 = RLC * EDEN * 2. * PI * R0 * R0 / (BETASQ * TM); // "CONSTANT FACTOR IN EQ.2.10.3"
		double EPS1 = T1 / T0;
		double EPSP1 = 1. - EPS1;
		double EPS2 = T2 / T0;
		double EPSP2 = 1. - EPS2;
		AMOLRM = CMOLL2
				* (C1 * (EPS2 - EPS1) + 1. / EPS1 - 1. / EPS2 + 1. / EPSP2 - 1.
						/ EPSP1 - C2 * Math.log(EPS2 * EPSP1 / (EPS1 * EPSP2))); // "EQ.2.10.3"

		return AMOLRM;
	}

	// Bremeatrahlung total cross section for a mixture of elements.
	// E0 is electorn energy
	/**
	 * Bremsstrahlung total cross section for a mixture of elements.
	 * @param E0 E0
	 * @return the result
	 */
	public static double BREMTM(double E0) {
		double BREMTM = 0.;
		if (E0 < AP + RM) {
			BREMTM = 0.;// System.out.println("enter");
		} else {
			// electron energy, lower photon cuttof, maximum photon energy
			BREMTM = BREMRM(E0, AP, E0 - RM);// System.out.println("enter");
		}
		return BREMTM;
	}

	// BREMRM is the bremsstrahlung cross section, integrated over some energy
	// range,
	// for a mixture of elements.
	/**
	 * Bremsstrahlung cross section, integrated over some energy range for a mixture of elements.
	 * @param E E
	 * @param K1 K1
	 * @param K2 K2
	 * @return the result
	 */
	public static double BREMRM(double E, double K1, double K2) {
		double BREMRM = 0.;
		for (int I = 0; I < NE; I++) {
			// E=electron energy, K1,K2=photon energy range
			BREMRM = BREMRM + PZ[I] * BREMRZ(Z[I], E, K1, K2);
		}
		return BREMRM;
	}

	// BREMRZ is Bremastrahlung cross section, integrated over some energy
	// range, for an
	// element.
	/**
	 * Bremsstrahlung cross section, integrated over some energy range, for an element.
	 * @param Z Z
	 * @param E E
	 * @param K1 K1
	 * @param K2 K2
	 * @return the result
	 */
	public static double BREMRZ(double Z, double E, double K1, double K2) {
		double BREMRZ = 0;// System.out.println("enter");
		@SuppressWarnings("unused")
		double DUMMY = BREMDZ(Z, E, K1);// guess for some further
										// initializations!!!
		funcname = "BREMFZ";// System.out.println("bremdz= "+DUMMY);
		BREMRZ = QD(funcname, K1, K2, "BREMFZ");// System.out.println("bremrz= "+BREMRZ);
		return BREMRZ;
	}

	// BREMDZ is Bremsatrahlung differential cross section for an element
	/**
	 * Bremsstrahlung differential cross section for an element.
	 * @param Z Z
	 * @param E E
	 * @param K K
	 * @return the result
	 */
	public static double BREMDZ(double Z, double E, double K) {
		// "***ALL ENTRIES TO THIS FUNCTION GIVE THE CONTRIBUTION THAT ELEMENT Z"
		// "   WOULD HAVE IF THERE WERE ONE PER MOLECULE.                      "
		// "   ENTRIES STARTING WITH D DO THEIR OWN INITIALIZATION.             "
		// "   ENTRIES STARTING WITH F RELY ON PREVIOUS D FOR INITIALIZATION.   "
		// "   BREMDZ.. D-SIGMA/D-K FOR BREMS IN Z                              "
		// "   BRMSDZ.. K*(D-SIGMA/D-K) FOR SOFT ENERGY LOSS FROM BREMS IN Z    "
		// "EVALUATES EQUATION 2.7.108 IN SLAC-265"NO 2.7.7 in SLAC!!!
		double BREMDZ = BRMSDZ(Z, E, K) / K;// is /K from 2.7.7 SLAC!!!OK!!!THE
											// MORE EXACT
		return BREMDZ;
	}

	// BRMSDZ is Soft bremestrahlung differential cross section for an element.
	/**
	 * Soft bremsstrahlung differential cross section for an element.
	 * @param Z Z
	 * @param EA EA
	 * @param K K
	 * @return the result
	 */
	public static double BRMSDZ(double Z, double EA, double K) {
		double BRMSDZ = 0.;
		E = EA;// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@GLOBAL@@@@@@@@@@@@@@@@@@@@@
		DELC = 136. * Math.pow(Z, -1. / 3.) * RM / E;
		CONST = APRIM(Z, E) * (AN * RHO / WM) * R0 * R0 * FSC * Z
				* (Z + XSIF(Z)) * RLC;// System.out.println("bremsfz= "+APRIM(Z,E));
		XLNZ = (4. / 3.) * Math.log(Z);
		if (E >= 50)
			XLNZ = XLNZ + 4. * FCOULC(Z);
		// ".....DELTAM IS THE DELTA AT WHICH THE SQUARE BRACKETS GO TO ZERO      "
		DELTAM = Math.exp((21.12 - XLNZ) / 4.184) - 0.952;// eq. 2.7.31 SLAC-265
		BRMSDZ = BRMSFZ(K);
		return BRMSDZ;
	}

	// One argument form of BRMSDZ.
	/**
	 * Called by BRMSDZ.
	 * @param K K
	 * @return the result
	 */
	public static double BRMSFZ(double K) {
		double BRMSFZ = 0.;
		// //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@GLOBAL@@@@@@@@@@@@@@@@@@@@@
		double EMKLOC = E - K;// electron energy after brems!!!
		double DELTA = 0.;
		double SB1 = 0.;
		double SB2 = 0.;
		double EE = 0.;
		if (EMKLOC == 0.0) {
			EMKLOC = 1.E-25;
		}
		DELTA = DELC * K / EMKLOC;// eq.2.7.9 SLAC 265
		if (DELTA >= DELTAM) {
			BRMSFZ = 0.0;// System.out.println("bremsfz= "+BRMSFZ);
		} else {
			if (DELTA <= 1.) {
				SB1 = 20.867 + DELTA * (-3.242 + DELTA * 0.625) - XLNZ;// eq
																		// 2.7.14
																		// SLAC
																		// 265
				SB2 = 20.209 + DELTA * (-1.930 + DELTA * (-0.086)) - XLNZ;// eq
																			// 2.7.15
																			// SLAC
																			// 265
			} else {
				SB1 = 21.12 - 4.184 * Math.log(DELTA + 0.952) - XLNZ;// eq
																		// 2.7.14
																		// SLAC
																		// 265
				SB2 = SB1;// eq 2.7.15 SLAC 265
			}
			EE = EMKLOC / E;
			BRMSFZ = CONST * ((1. + EE * EE) * SB1 - 0.666667 * EE * SB2);// System.out.println("bremsfz= "+CONST);
		}// System.out.println("bremsfz= "+BRMSFZ);
		return BRMSFZ;
	}

	/**
	 * This is just BRMSFZ(K) / K.
	 * @param K K
	 * @return the result
	 */
	public static double BREMFZ(double K) {
		double BREMFZ = BRMSFZ(K) / K;// System.out.println("en");
		return BREMFZ;
	}

	// function test for some verification with DCADRE
	/**
	 * Test function for DCADRE.
	 * @param K K
	 * @return the result
	 */
	public static double test(double K) {
		double test = 10.0 * K + 1.0;
		return test;
	}

	// "     EMPIRICAL CORRECTION FACTOR TO BREMS CROSS SECTION               "
	/**
	 * EMPIRICAL CORRECTION FACTOR TO BREMS CROSS SECTION
	 * @param Z Z
	 * @param E E
	 * @return the result
	 */
	public static double APRIM(double Z, double E) {
		// " This version can be switched to use different values:                "
		// "   IAPRIM = 0  equivalent to old PEGS4 (default)                      "
		// "            1  reads in values from unit 22                           "
		// "            2  sets APRIM to 1.0                                      "
		// " Future changes can be accommodated by reading in                     "
		// " different data on unit 22 and if necessary changing the array sizes: "
		double APRIM = 0.;
		// int $NAPRE=115;//" Maximum number of energies ( > 18 )      "
		// int $NAPRZ=14;// " Maximum number of elements ( > 5 )       "
		int $NAPR1 = $NAPRE - 18;
		// int $NAPR2 = $NAPRZ - 5;
		// int $NAPR3 = $NAPRE * $NAPR2;
		// double[][] APRIMD1=new double[$NAPRE][$NAPRZ];
		APRIMD1[0][0] = 1.32;
		APRIMD1[1][0] = 1.26;
		APRIMD1[2][0] = 1.18;
		APRIMD1[3][0] = 1.13;
		APRIMD1[4][0] = 1.09;
		APRIMD1[5][0] = 1.07;
		APRIMD1[6][0] = 1.05;
		APRIMD1[7][0] = 1.04;
		APRIMD1[8][0] = 1.03;
		APRIMD1[9][0] = 1.02;
		for (int i = 1; i <= 8; i++)
			APRIMD1[9 + i][0] = 1.00;
		for (int i = 1; i <= $NAPR1; i++)
			APRIMD1[9 + 8 + i][0] = 0.00;
		APRIMD1[0][1] = 1.34;
		APRIMD1[1][1] = 1.27;
		APRIMD1[2][1] = 1.19;
		APRIMD1[3][1] = 1.13;
		APRIMD1[4][1] = 1.09;
		APRIMD1[5][1] = 1.07;
		APRIMD1[6][1] = 1.05;
		APRIMD1[7][1] = 1.04;
		APRIMD1[8][1] = 1.03;
		APRIMD1[9][1] = 1.02;
		for (int i = 1; i <= 8; i++)
			APRIMD1[9 + i][1] = 1.00;
		for (int i = 1; i <= $NAPR1; i++)
			APRIMD1[9 + 8 + i][1] = 0.00;
		APRIMD1[0][2] = 1.39;
		APRIMD1[1][2] = 1.30;
		APRIMD1[2][2] = 1.21;
		APRIMD1[3][2] = 1.14;
		APRIMD1[4][2] = 1.10;
		APRIMD1[5][2] = 1.07;
		APRIMD1[6][2] = 1.05;
		APRIMD1[7][2] = 1.04;
		APRIMD1[8][2] = 1.03;
		APRIMD1[9][2] = 1.02;
		APRIMD1[10][2] = 0.994;
		APRIMD1[11][2] = 0.991;
		APRIMD1[12][2] = 0.991;
		APRIMD1[13][2] = 0.990;
		APRIMD1[14][2] = 0.989;
		APRIMD1[15][2] = 0.989;
		APRIMD1[16][2] = 0.988;
		APRIMD1[17][2] = 0.988;
		for (int i = 1; i <= $NAPR1; i++)
			APRIMD1[17 + i][2] = 0.00;
		APRIMD1[0][3] = 1.46;
		APRIMD1[1][3] = 1.34;
		APRIMD1[2][3] = 1.23;
		APRIMD1[3][3] = 1.15;
		APRIMD1[4][3] = 1.11;
		APRIMD1[5][3] = 1.08;
		APRIMD1[6][3] = 1.06;
		APRIMD1[7][3] = 1.05;
		APRIMD1[8][3] = 1.03;
		APRIMD1[9][3] = 1.02;
		APRIMD1[10][3] = 0.989;
		APRIMD1[11][3] = 0.973;
		APRIMD1[12][3] = 0.971;
		APRIMD1[13][3] = 0.969;
		APRIMD1[14][3] = 0.967;
		APRIMD1[15][3] = 0.965;
		APRIMD1[16][3] = 0.963;
		APRIMD1[17][3] = 0.963;
		for (int i = 1; i <= $NAPR1; i++)
			APRIMD1[17 + i][3] = 0.00;
		APRIMD1[0][4] = 1.55;
		APRIMD1[1][4] = 1.40;
		APRIMD1[2][4] = 1.26;
		APRIMD1[3][4] = 1.17;
		APRIMD1[4][4] = 1.12;
		APRIMD1[5][4] = 1.09;
		APRIMD1[6][4] = 1.07;
		APRIMD1[7][4] = 1.05;
		APRIMD1[8][4] = 1.03;
		APRIMD1[9][4] = 1.02;
		APRIMD1[10][4] = 0.955;
		APRIMD1[11][4] = 0.935;
		APRIMD1[12][4] = 0.930;
		APRIMD1[13][4] = 0.925;
		APRIMD1[14][4] = 0.920;
		APRIMD1[15][4] = 0.915;
		APRIMD1[16][4] = 0.911;
		APRIMD1[17][4] = 0.911;
		for (int i = 1; i <= $NAPR1; i++)
			APRIMD1[17 + i][4] = 0.00;
		for (int j = 5; j < $NAPRZ; j++) {
			for (int i = 0; i < $NAPRE; i++) {
				APRIMD1[i][j] = 0.00;
			}
		}
		// double[] EPRIM1=new double[$NAPRE];
		EPRIM1[0] = 2.;
		EPRIM1[1] = 3.;
		EPRIM1[2] = 4.;
		EPRIM1[3] = 5.;
		EPRIM1[4] = 6.;
		EPRIM1[5] = 7.;
		EPRIM1[6] = 8.;
		EPRIM1[7] = 9.;
		EPRIM1[8] = 10.;
		EPRIM1[9] = 11.;
		EPRIM1[10] = 21.;
		EPRIM1[11] = 31.;
		EPRIM1[12] = 41.;
		EPRIM1[13] = 51.;
		EPRIM1[14] = 61.;
		EPRIM1[15] = 71.;
		EPRIM1[16] = 81.;
		EPRIM1[17] = 91.;
		for (int i = 1; i <= $NAPR1; i++)
			EPRIM1[17 + i] = 0.00;
		// double[] ZPRIM1=new double[$NAPRZ];
		ZPRIM1[0] = 6.;
		ZPRIM1[1] = 13.;
		ZPRIM1[2] = 29.;
		ZPRIM1[3] = 50.;
		ZPRIM1[4] = 79.;
		ZPRIM1[5] = 0.;
		ZPRIM1[6] = 0.;
		ZPRIM1[7] = 0.;
		ZPRIM1[8] = 0.;
		ZPRIM1[9] = 0.;
		ZPRIM1[10] = 0.;
		ZPRIM1[11] = 0.;
		ZPRIM1[12] = 0.;
		ZPRIM1[13] = 0.;
		// double [] APRIMZ=new double[$NAPRE];
		double EM = 0.;

		if (IAPRIM == 0) {// " PEGS4 default APRIM"
							// if(IAPRFL == 0)
							// {
							// IAPRFL=1;
							// '0IAPRIM=0, i.e. uses KOCH AND MOTZ empirical
							// corrections to',
							// ' brem cross section'/;
			// }
			if (E >= 50) {
				APRIM = 1.;
			} else {
				// " INTERPOLATE APRIM OVER Z "
				EM = E / RM;
				for (int IE = 0; IE < 18; IE++) {
					// double [] xa = new double[NPHE[IZ]];
					double[] ya = new double[$NAPRZ];
					for (int i = 0; i < $NAPRZ; i++) {
						// xa[i]=PHE[i][IZ];
						ya[i] = APRIMD1[IE][i];
					}
					// ------------------------------this($NAPRE)=no use in
					// AINTP
					// first 5 columns and 18 rows in ya is !=0
					APRIMZ[IE] = AINTP(Z, ZPRIM1, 5, ya, $NAPRE, false, false);
				} // " Z INTERPOLATION IS NOW COMPLETE. NOW DO ENERGY "
					// 18=to be consistent with above!!
				APRIM = AINTP(EM, EPRIM1, 18, APRIMZ, 1, false, false);
			}
		} else if (IAPRIM == 1) {
			String filename = "aprime.data";// System.out.println("enter");
			// readAprimeData(filename);
			if (IAPRFL == 0) {// " read in data from APRIME.DATA"
								// uses NRC(based on NIST/ICRU) corrections to
								// brem cross section
				readAprimeData(filename);
				IAPRFL = 1;
				for (int IE = 0; IE < NAPRE; IE++) {
					EPRIM[IE] = 1. + EPRIM[IE] / RM;// System.out.println(EPRIM[IE]);
				}
			}//
			EM = E / RM;
			for (int IE = 0; IE < NAPRE; IE++) {// " INTERPOLATE APRIM OVER
												// LOG(Z)
				double[] ya = new double[$NAPRZ];
				for (int i = 0; i < $NAPRZ; i++) {
					// xa[i]=PHE[i][IZ];
					ya[i] = APRIMD[IE][i];// System.out.println(ya[i]);
				}
				APRIMZ[IE] = AINTP(Z, ZPRIM, NAPRZ, ya, $NAPRE, true, false);
			} // " NOW DO ENERGY INTERPOLATION     "
			APRIM = AINTP(EM, EPRIM, NAPRE, APRIMZ, 1, false, false);
		} else if (IAPRIM == 2)// uses NO corrections to brem
		{
			if (IAPRFL == 0) {
				IAPRFL = 1;
			}
			APRIM = 1.0;
		} else// ILLEGAL VALUE FOR IAPRIM
		{
			// never happen
		}

		return APRIM;
	}

	// DCADRE input data
	/**
	 * Integration using EGSnrc routine DCADRE. Integrates F(x) from A to B using 
	 * cautious adaptive Romberg extrapolation.
	 * @param F the function to be integrated
	 * @param A A
	 * @param B B
	 * @param MSG MSG, not used
	 * @return the result
	 */
	public static double QD(String F, double A, double B, String MSG) {
		double QD = 0.;
		// boolean first_time=true;
		double ADUM = A;
		double BDUM = B;
		double ERRDUM = 0.;
		int IER = 0;
		QD = DCADRE(F, ADUM, BDUM, $AERR, $RERR, ERRDUM, IER);
		IER = IERq;
		ERRDUM = ERRORq;
		if (IER > 66) {
			if (first_time) {
				first_time = false;
				// output_file(lenfn-7:lenfn) = 'pegs4err';
				// open(10,file=output_file,status='unknown');
			}
			// $UOUTPUT(10) IER,MSG,A,B,QD,ERRDUM;
			// (' DCADRE CODE=',I4,' FOR INTEGRAL ',A6,' FROM ',1P,G14.6,' TO
			// ',G14.6,
			// ',QD=',G14.6,'+-',G14.6);
		}
		return QD;
	}

	/**
	 * Used by QD
	 * @param F F
	 * @param A A
	 * @param B B
	 * @param AERR AERR
	 * @param RERR RERR
	 * @param ERROR ERROR
	 * @param IER IER
	 * @return the result
	 */
	public static double DCADRE(String F, double A, double B, double AERR,
			double RERR, double ERROR, int IER) {
		// "------------------------------------------------------------------"
		// "-DCADRE--------D-------LIBRARY 1----------------------------------"
		// "------------------------------------------------------------------"
		// "                                                                  "
		// "FUNCTION:          - INTEGRATE F(X) FROM A TO B, USING CAUTIOUS   "
		// "                     ADAPTIVE ROMBERG EXTRAPOLATION.              "
		// "                                                                  "
		// "USAGE:             - FUNCTION DCADRE(F,A,B,AERR,RERR,ERROR,IER)   "
		// "                                                                  "
		// "PARAMETERS: DCADRE - ESTIMATE OF THE INTEGRAL OF F(X) FROM A TO B."
		// "                                                                  "
		// "            F      - A SINGLE-ARGUMENT REAL FUNCTION SUBPROGRAM   "
		// "                     SUPPLIED BY THE USER.  F MUST BE DECLARED    "
		// "                     EXTERNAL IN THE CALLING PROGRAM.             "
		// "                                                                  "
		// "            A,B    - THE TWO ENDPOINTS OF THE INTERVAL OF         "
		// "                     INTEGRATION (INPUT).                         "
		// "                                                                  "
		// "            AERR   - DESIRED ABSOLUTE ERROR IN THE ANSWER (INPUT)."
		// "                                                                  "
		// "            RERR   - DESIRED RELATIVE ERROR IN THE ANSWER (INPUT)."
		// "                                                                  "
		// "            ERROR  - ESTIMATED BOUND ON THE ABSOLUTE ERROR OF     "
		// "                     THE OUTPUT NUMBER, DCADRE.                   "
		// "                                                                  "
		// "            IER    - ERROR PARAMETER                              "
		// "                                                                  "
		// "                     WARNING ERROR(WITH FIX) = 64 + N             "
		// "                                                                  "
		// "                       N = 1 IMPLIES THAT ONE OR MORE SINGULAR-   "
		// "                             ITIES WERE SUCCESSFULLY HANDLED.     "
		// "                                                                  "
		// "                       N = 2 IMPLIES THAT, IN SOME SUBINTERVAL(S),"
		// "                             THE ESTIMATE OF THE INTEGRAL WAS     "
		// "                             ACCEPTED MERELY BECAUSE THE ESTIMATED"
		// "                             ERROR WAS SMALL, EVEN THOUGH NO REG- "
		// "                             ULAR BEHAVIOR WAS RECOGNIZED.        "
		// "                                                                  "
		// "                     TERMINAL ERROR = 128 + N                     "
		// "                                                                  "
		// "                       N = 3 FAILURE DUE TO INSUFFICIENT INTERNAL "
		// "                             WORKING STORAGE.                     "
		// "                                                                  "
		// "                       N = 4 FAILURE.  THIS MAY BE DUE TO TOO MUCH"
		// "                             NOISE IN THE FUNCTION (RELATIVE TO   "
		// "                             THE GIVEN ERROR REQUIREMENTS) OR DUE "
		// "                             TO AN ILL-BEHAVED INTEGRAND.         "
		// "                                                                  "
		// "                       N = 5 INDICATES THAT RERR IS GREATER THAN  "
		// "                             0.1, OR RERR IS LESS THAN 0.0, OR    "
		// "                             RERR IS TOO SMALL FOR THE PRECISION  "
		// "                             OF THE MACHINE.                      "
		// "                                                                  "
		// "------------------------------------------------------------------"
		// "VERSION DATE:      - 8 OCTOBER 1974                               "
		// "                                                                  "
		// "MORTRAN VERSION    - 4 OCTOBER 1984/1545 (W. R. NELSON)           "
		// "------------------------------------------------------------------"

		Aq = A;
		Bq = B;
		AERRq = AERR;
		RERRq = RERR;
		ERRORq = ERROR;
		IERq = IER;
		// --------------
		resetDCADRE();// reset global values!!
		// test---------------------------------------
		// DCADRE = F(23);//WORKS
		// --------------------------------------------
		if (LENGTHq == ZEROq)// GO TO 215;
		{
			goto215q();
			return DCADREq;// 9005 RETURN;
		}
		if ((RERR > P1q) || (RERR < ZEROq))// GO TO 210;
		{
			IERq = 133;// 210 IER=133;
			goto215q();
			return DCADREq;// 9005 RETURN;
		}
		if ((AERR == ZEROq) && ((RERR + HUNq) <= HUNq)) // GO TO 210;
		{
			IERq = 133;// 210 IER=133;
			goto215q();
			return DCADREq;// 9005 RETURN;
		}

		RIGHTq = false;// 5
						// #########################################################RIGHT=.FALSE.;
		// " INVESTIGATION OF A PARTICULAR       "
		// " SUBINTERVAL BEGINS AT THIS POINT. "
		STEPq = ENDq - BEGq;// 10#################################################10
							// STEP=END - BEG;
		ASTEPq = DABS(STEPq);
		if (ASTEPq < STEPMNq)// GO TO 205;
		{
			IERq = 132;// 205 IER=132;
			goto215q();
			return DCADREq;
		}
		if ((STEPNMq + ASTEPq) == STEPNMq)// GO TO 205;
		{
			IERq = 132;// 205 IER=132;
			goto215q();
			return DCADREq;
		}
		Tq[0][0] = FBEGq + FENDq;// T(1,1)=FBEG + FEND;
		TABSq = DABS(FBEGq) + DABS(FENDq);
		Lq = 1;
		Nq = 1;
		H2CONVq = false;
		AITKENq = false;
		LM1q = Lq;// 15########################################################################
					// LM1=L;
		Lq = Lq + 1;
		// "  CALCULATE THE NEXT TRAPEZOID SUM,   "
		// "  T(L,1), WHICH IS BASED ON *N2* + 1"
		// "  EQUISPACED POINTS. HERE,          "
		// "  N2 = N*2 = 2**(L-1).              "
		N2q = Nq + Nq;
		FNq = N2q;
		ISTEPq = (IENDq - IBEGq) / Nq;
		if (ISTEPq <= 1) // IF(ISTEP>1) GO TO 25;
		{
			notgotto25q();
			if (exitq) {
				IERq = 131;// 200 IER=131;
				return DCADREq;
			}
		}
		ISTEP2q = IBEGq + ISTEPq / 2;// #25 ISTEP2=IBEG + ISTEP/2;
		SUMq = ZEROq;
		SUMABSq = ZEROq;
		for (Iq = ISTEP2q - 1; Iq < IENDq; Iq = Iq + ISTEPq) {
			SUMq = SUMq + TSq[Iq];
			SUMABSq = SUMABSq + DABS(TSq[Iq]);
		}
		Tq[Lq - 1][0] = Tq[Lq - 2][0] * HALFq + SUMq / FNq;// T(L,1)=T(L-1,1)*HALF+SUM/FN;
		TABSq = TABSq * HALFq + SUMABSq / FNq;
		ABSIq = ASTEPq * TABSq;
		Nq = N2q;
		// "     GET PRELIMINARY VALUE FOR *VINT*    "
		// "     FROM LAST TRAPEZOID SUM AND UPDATE"
		// "     THE ERROR REQUIREMENT *ERGOAL*    "
		// "     FOR THIS SUBINTERVAL.             "
		ITq = 1;
		VINTq = STEPq * Tq[Lq - 1][0];// STEP*T(L,1);
		TABTLMq = TABSq * TENq;
		FNSIZEq = DMAX1(FNSIZEq, DABS(Tq[Lq - 1][0]));// DMAX1(FNSIZE,DABS(T(L,1)));
		ERGLq = ASTEPq * FNSIZEq * TENq;
		ERGOALq = STAGEq * DMAX1(ERRAq, ERRRq * DABS(CURESTq + VINTq));
		// "  COMPLETE ROW L AND COLUMN L OF *T*  "
		// "    ARRAY.                            "
		FEXTRPq = ONEq;
		for (Iq = 0; Iq < LM1q; Iq++) {
			FEXTRPq = FEXTRPq * FOURq;
			Tq[Iq][Lq - 1] = Tq[Lq - 1][Iq] - Tq[Lq - 2][Iq];// T(I,L)=T(L,I) -
																// T(L-1,I);
			Tq[Lq - 1][Iq + 1] = Tq[Lq - 1][Iq] + Tq[Iq][Lq - 1]
					/ (FEXTRPq - ONEq);// T(L,I+1)=T(L,I) + T(I,L)/(FEXTRP-ONE);
		}
		ERRERq = ASTEPq * DABS(Tq[0][Lq - 1]);// ERRER=ASTEP*DABS(T(1,L));
		// "  PRELIMINARY DECISION PROCEDURE      "
		// "  IF L = 2 AND T(2,1) = T(1,1),     "
		// "  GO TO 135 TO FOLLOW UP THE        "
		// "  IMPRESSION THAT INTERGRAND IS     "
		// "  STRAIGHT LINE.                    "
		if (Lq > 2)// if (L>2)// GO TO 40;
		{
			goto40q();
			if (exitq) {
				return DCADREq;
			}
		}

		if ((TABSq + P1q * DABS(Tq[0][1])) == TABSq) // IF(TABS+P1*DABS(T(1,2)).EQ.TABS)
														// GO TO 135;
		{
			goto135q();
			if (exitq) {
				return DCADREq;
			}
		}
		// "                              CACULATE NEXT RATIOS FOR            "
		// "                                COLUMNS 1,...,L-2 OF T-TABLE      "
		// "                                RATIO IS SET TO ZERO IF DIFFERENCE"
		// "                                IN LAST TWO ENTRIES OF COLUMN IS  "
		// "                                ABOUT ZERO                        "
		goto15q();
		if (exitq)
			return DCADREq;

		return DCADREq;
	}

	/**
	 * Used by DCADRE
	 */
	private static void goto5q() {
		RIGHTq = false;// 5
						// #########################################################RIGHT=.FALSE.;
		// " INVESTIGATION OF A PARTICULAR       "
		// " SUBINTERVAL BEGINS AT THIS POINT. "
		STEPq = ENDq - BEGq;// 10#################################################10
							// STEP=END - BEG;
		ASTEPq = DABS(STEPq);
		if (ASTEPq < STEPMNq)// GO TO 205;
		{
			IERq = 132;// 205 IER=132;
			goto215q();
			exitq = true;
			return;
		}
		if ((STEPNMq + ASTEPq) == STEPNMq)// GO TO 205;
		{
			IERq = 132;// 205 IER=132;
			goto215q();
			exitq = true;
			return;
		}
		Tq[0][0] = FBEGq + FENDq;// T(1,1)=FBEG + FEND;
		TABSq = DABS(FBEGq) + DABS(FENDq);
		Lq = 1;
		Nq = 1;
		H2CONVq = false;
		AITKENq = false;
		LM1q = Lq;// 15########################################################################
					// LM1=L;
		Lq = Lq + 1;
		// "  CALCULATE THE NEXT TRAPEZOID SUM,   "
		// "  T(L,1), WHICH IS BASED ON *N2* + 1"
		// "  EQUISPACED POINTS. HERE,          "
		// "  N2 = N*2 = 2**(L-1).              "
		N2q = Nq + Nq;
		FNq = N2q;
		ISTEPq = (IENDq - IBEGq) / Nq;
		if (ISTEPq <= 1) // IF(ISTEP>1) GO TO 25;
		{
			notgotto25q();
			if (exitq) {
				IERq = 131;// 200 IER=131;
				return;
			}
		}
		ISTEP2q = IBEGq + ISTEPq / 2;// #25 ISTEP2=IBEG + ISTEP/2;
		SUMq = ZEROq;
		SUMABSq = ZEROq;
		for (Iq = ISTEP2q - 1; Iq < IENDq; Iq = Iq + ISTEPq) {
			SUMq = SUMq + TSq[Iq];
			SUMABSq = SUMABSq + DABS(TSq[Iq]);
		}
		Tq[Lq - 1][0] = Tq[Lq - 2][0] * HALFq + SUMq / FNq;// T(L,1)=T(L-1,1)*HALF+SUM/FN;
		TABSq = TABSq * HALFq + SUMABSq / FNq;
		ABSIq = ASTEPq * TABSq;
		Nq = N2q;
		// "     GET PRELIMINARY VALUE FOR *VINT*    "
		// "     FROM LAST TRAPEZOID SUM AND UPDATE"
		// "     THE ERROR REQUIREMENT *ERGOAL*    "
		// "     FOR THIS SUBINTERVAL.             "
		ITq = 1;
		VINTq = STEPq * Tq[Lq - 1][0];// STEP*T(L,1);
		TABTLMq = TABSq * TENq;
		FNSIZEq = DMAX1(FNSIZEq, DABS(Tq[Lq - 1][0]));// DMAX1(FNSIZE,DABS(T(L,1)));
		ERGLq = ASTEPq * FNSIZEq * TENq;
		ERGOALq = STAGEq * DMAX1(ERRAq, ERRRq * DABS(CURESTq + VINTq));
		// "  COMPLETE ROW L AND COLUMN L OF *T*  "
		// "    ARRAY.                            "
		FEXTRPq = ONEq;
		for (Iq = 0; Iq < LM1q; Iq++) {
			FEXTRPq = FEXTRPq * FOURq;
			Tq[Iq][Lq - 1] = Tq[Lq - 1][Iq] - Tq[Lq - 2][Iq];// T(I,L)=T(L,I) -
																// T(L-1,I);
			Tq[Lq - 1][Iq + 1] = Tq[Lq - 1][Iq] + Tq[Iq][Lq - 1]
					/ (FEXTRPq - ONEq);// T(L,I+1)=T(L,I) + T(I,L)/(FEXTRP-ONE);
		}
		ERRERq = ASTEPq * DABS(Tq[0][Lq - 1]);// ERRER=ASTEP*DABS(T(1,L));
		// "  PRELIMINARY DECISION PROCEDURE      "
		// "  IF L = 2 AND T(2,1) = T(1,1),     "
		// "  GO TO 135 TO FOLLOW UP THE        "
		// "  IMPRESSION THAT INTERGRAND IS     "
		// "  STRAIGHT LINE.                    "
		if (Lq > 2)// if (L>2)// GO TO 40;
		{
			goto40q();
			if (exitq) {
				return;
			}
		}

		if ((TABSq + P1q * DABS(Tq[0][1])) == TABSq) // IF(TABS+P1*DABS(T(1,2)).EQ.TABS)
														// GO TO 135;
		{
			goto135q();
			if (exitq) {
				return;
			}
		}
		goto15q();
		if (exitq)
			return;
	}

	/**
	 * Used by DCADRE
	 */
	private static void goto10q() {
		STEPq = ENDq - BEGq;// 10#################################################10
							// STEP=END - BEG;
		ASTEPq = DABS(STEPq);
		if (ASTEPq < STEPMNq)// GO TO 205;
		{
			IERq = 132;// 205 IER=132;
			goto215q();
			exitq = true;
			return;
		}
		if ((STEPNMq + ASTEPq) == STEPNMq)// GO TO 205;
		{
			IERq = 132;// 205 IER=132;
			goto215q();
			exitq = true;
			return;
		}
		Tq[0][0] = FBEGq + FENDq;// T(1,1)=FBEG + FEND;
		TABSq = DABS(FBEGq) + DABS(FENDq);
		Lq = 1;
		Nq = 1;
		H2CONVq = false;
		AITKENq = false;
		LM1q = Lq;// 15########################################################################
					// LM1=L;
		Lq = Lq + 1;
		// "  CALCULATE THE NEXT TRAPEZOID SUM,   "
		// "  T(L,1), WHICH IS BASED ON *N2* + 1"
		// "  EQUISPACED POINTS. HERE,          "
		// "  N2 = N*2 = 2**(L-1).              "
		N2q = Nq + Nq;
		FNq = N2q;
		ISTEPq = (IENDq - IBEGq) / Nq;
		if (ISTEPq <= 1) // IF(ISTEP>1) GO TO 25;
		{
			notgotto25q();
			if (exitq) {
				IERq = 131;// 200 IER=131;
				return;
			}
		}
		ISTEP2q = IBEGq + ISTEPq / 2;// #25 ISTEP2=IBEG + ISTEP/2;
		SUMq = ZEROq;
		SUMABSq = ZEROq;
		for (Iq = ISTEP2q - 1; Iq < IENDq; Iq = Iq + ISTEPq) {
			SUMq = SUMq + TSq[Iq];
			SUMABSq = SUMABSq + DABS(TSq[Iq]);
		}
		Tq[Lq - 1][0] = Tq[Lq - 2][0] * HALFq + SUMq / FNq;// T(L,1)=T(L-1,1)*HALF+SUM/FN;
		TABSq = TABSq * HALFq + SUMABSq / FNq;
		ABSIq = ASTEPq * TABSq;
		Nq = N2q;
		// "     GET PRELIMINARY VALUE FOR *VINT*    "
		// "     FROM LAST TRAPEZOID SUM AND UPDATE"
		// "     THE ERROR REQUIREMENT *ERGOAL*    "
		// "     FOR THIS SUBINTERVAL.             "
		ITq = 1;
		VINTq = STEPq * Tq[Lq - 1][0];// STEP*T(L,1);
		TABTLMq = TABSq * TENq;
		FNSIZEq = DMAX1(FNSIZEq, DABS(Tq[Lq - 1][0]));// DMAX1(FNSIZE,DABS(T(L,1)));
		ERGLq = ASTEPq * FNSIZEq * TENq;
		ERGOALq = STAGEq * DMAX1(ERRAq, ERRRq * DABS(CURESTq + VINTq));
		// "  COMPLETE ROW L AND COLUMN L OF *T*  "
		// "    ARRAY.                            "
		FEXTRPq = ONEq;
		for (Iq = 0; Iq < LM1q; Iq++) {
			FEXTRPq = FEXTRPq * FOURq;
			Tq[Iq][Lq - 1] = Tq[Lq - 1][Iq] - Tq[Lq - 2][Iq];// T(I,L)=T(L,I) -
																// T(L-1,I);
			Tq[Lq - 1][Iq + 1] = Tq[Lq - 1][Iq] + Tq[Iq][Lq - 1]
					/ (FEXTRPq - ONEq);// T(L,I+1)=T(L,I) + T(I,L)/(FEXTRP-ONE);
		}
		ERRERq = ASTEPq * DABS(Tq[0][Lq - 1]);// ERRER=ASTEP*DABS(T(1,L));
		// "  PRELIMINARY DECISION PROCEDURE      "
		// "  IF L = 2 AND T(2,1) = T(1,1),     "
		// "  GO TO 135 TO FOLLOW UP THE        "
		// "  IMPRESSION THAT INTERGRAND IS     "
		// "  STRAIGHT LINE.                    "
		if (Lq > 2)// if (L>2)// GO TO 40;
		{
			goto40q();
			if (exitq) {
				return;
			}
		}

		if ((TABSq + P1q * DABS(Tq[0][1])) == TABSq) // IF(TABS+P1*DABS(T(1,2)).EQ.TABS)
														// GO TO 135;
		{
			goto135q();
			if (exitq) {
				return;
			}
		}
		goto15q();
		if (exitq)
			return;
	}

	/**
	 * Used by DCADRE
	 */
	private static void goto15q() {
		LM1q = Lq;// 15########################################################################
					// LM1=L;
		Lq = Lq + 1;
		// "  CALCULATE THE NEXT TRAPEZOID SUM,   "
		// "  T(L,1), WHICH IS BASED ON *N2* + 1"
		// "  EQUISPACED POINTS. HERE,          "
		// "  N2 = N*2 = 2**(L-1).              "
		N2q = Nq + Nq;
		FNq = N2q;
		ISTEPq = (IENDq - IBEGq) / Nq;
		if (ISTEPq <= 1) // IF(ISTEP>1) GO TO 25;
		{
			notgotto25q();
			if (exitq) {
				IERq = 131;// 200 IER=131;
				return;
			}
		}
		ISTEP2q = IBEGq + ISTEPq / 2;// #25 ISTEP2=IBEG + ISTEP/2;
		SUMq = ZEROq;
		SUMABSq = ZEROq;
		for (Iq = ISTEP2q - 1; Iq < IENDq; Iq = Iq + ISTEPq) {
			SUMq = SUMq + TSq[Iq];
			SUMABSq = SUMABSq + DABS(TSq[Iq]);
		}
		Tq[Lq - 1][0] = Tq[Lq - 2][0] * HALFq + SUMq / FNq;// T(L,1)=T(L-1,1)*HALF+SUM/FN;
		TABSq = TABSq * HALFq + SUMABSq / FNq;
		ABSIq = ASTEPq * TABSq;
		Nq = N2q;
		// "     GET PRELIMINARY VALUE FOR *VINT*    "
		// "     FROM LAST TRAPEZOID SUM AND UPDATE"
		// "     THE ERROR REQUIREMENT *ERGOAL*    "
		// "     FOR THIS SUBINTERVAL.             "
		ITq = 1;
		VINTq = STEPq * Tq[Lq - 1][0];// STEP*T(L,1);
		TABTLMq = TABSq * TENq;
		FNSIZEq = DMAX1(FNSIZEq, DABS(Tq[Lq - 1][0]));// DMAX1(FNSIZE,DABS(T(L,1)));
		ERGLq = ASTEPq * FNSIZEq * TENq;
		ERGOALq = STAGEq * DMAX1(ERRAq, ERRRq * DABS(CURESTq + VINTq));
		// "  COMPLETE ROW L AND COLUMN L OF *T*  "
		// "    ARRAY.                            "
		FEXTRPq = ONEq;
		for (Iq = 0; Iq < LM1q; Iq++) {
			FEXTRPq = FEXTRPq * FOURq;
			Tq[Iq][Lq - 1] = Tq[Lq - 1][Iq] - Tq[Lq - 2][Iq];// T(I,L)=T(L,I) -
																// T(L-1,I);
			Tq[Lq - 1][Iq + 1] = Tq[Lq - 1][Iq] + Tq[Iq][Lq - 1]
					/ (FEXTRPq - ONEq);// T(L,I+1)=T(L,I) + T(I,L)/(FEXTRP-ONE);
		}
		ERRERq = ASTEPq * DABS(Tq[0][Lq - 1]);// ERRER=ASTEP*DABS(T(1,L));
		// "  PRELIMINARY DECISION PROCEDURE      "
		// "  IF L = 2 AND T(2,1) = T(1,1),     "
		// "  GO TO 135 TO FOLLOW UP THE        "
		// "  IMPRESSION THAT INTERGRAND IS     "
		// "  STRAIGHT LINE.                    "
		if (Lq > 2)// if (L>2)// GO TO 40;
		{
			goto40q();
			if (exitq) {
				return;
			}
		}

		if ((TABSq + P1q * DABS(Tq[0][1])) == TABSq) // IF(TABS+P1*DABS(T(1,2)).EQ.TABS)
														// GO TO 135;
		{
			goto135q();
			if (exitq) {
				return;
			}
		}
		goto15q();
		if (exitq)
			return;
	}

	/**
	 * Used by DCADRE
	 */
	private static void notgotto25q() {
		IIq = IENDq;
		IENDq = IENDq + Nq;
		if (IENDq > MAXTSq)// IF(IEND.GT.MAXTS) GO TO 200;
		{
			goto215q();
			exitq = true;// return DCADREq;
			return;
		}
		HOVNq = STEPq / FNq;
		IIIq = IENDq;
		FIq = ONEq;
		for (Iq = 0; Iq < N2q; Iq = Iq + 2) {
			TSq[IIIq - 1] = TSq[IIq - 1];// TS[III]=TS(II);
			RVALq = ENDq - FIq * HOVNq;
			TSq[IIIq - 2] = F(RVALq);// TS(III-1)=F(RVAL);@@@
			FIq = FIq + TWOq;
			IIIq = IIIq - 2;
			IIq = IIq - 1;
		}
		ISTEPq = 2;
	}

	/**
	 * Used by DCADRE
	 */
	private static void goto40q() {
		for (Iq = 1; Iq < LM1q; Iq++)// 40 DO 45 I=2,LM1;pana la 45
		{
			DIFFq = ZEROq;
			if (TABTLMq + DABS(Tq[Iq - 1][Lq - 1]) != TABTLMq)
				DIFFq = Tq[Iq - 1][LM1q - 1] / Tq[Iq - 1][Lq - 1];
			// IF(TABTLM+DABS(T(I-1,L)).NE.TABTLM) DIFF=T(I-1,LM1)/T(I-1,L);
			Tq[Iq - 1][LM1q - 1] = DIFFq;// T(I-1,LM1)=DIFF;
		}
		// 45 CONTINUE;
		// if (DABS(FOURq-Tq(1,LM1)).LE.H2TOL) GO TO 60;
		if (DABS(FOURq - Tq[0][LM1q - 1]) <= H2TOLq) {
			goto60q();
			if (exitq)
				return;
		}
		// IF(T(1,LM1).EQ.ZERO) GO TO 55;
		if (Tq[0][LM1q - 1] == ZEROq) {
			goto55q();
			if (exitq)
				return;
		}
		// IF(DABS(TWO-DABS(T(1,LM1))).LT.JUMPTL) GO TO 130;
		if (DABS(TWOq - DABS(Tq[0][LM1q - 1])) < JUMPTLq) {
			goto130q();
			if (exitq)
				return;
		}
		// IF(L.EQ.3) GO TO 15;
		if (Lq == 3) {
			goto15q();
			if (exitq)
				return;
		}
		H2CONVq = false;
		// IF(DABS((T(1,LM1)-T(1,L-2))/T(1,LM1)).LE.AITTOL) GO TO 75;
		if (DABS((Tq[0][LM1q - 1] - Tq[0][Lq - 3]) / Tq[0][LM1q - 1]) <= AITTOLq) {
			goto75q();
			if (exitq)
				return;
		}
		// 50 IF(REGLAR) GO TO 55;
		if (REGLARq) {
			goto55q();
			if (exitq)
				return;
		}
		// IF(L.EQ.4) GO TO 15;
		if (Lq == 4) {
			goto15q();
			if (exitq)
				return;
		}
		// 55 IF(ERRER.GT.ERGOAL.AND.(ERGL+ERRER).NE.ERGL) GO TO 175;
		if ((ERRERq > ERGOALq) && (ERGLq + ERRERq) != ERGLq) {
			goto175q();
			if (exitq)
				return;
		}
		goto145q();
		if (exitq)
			return;
	}

	/**
	 * Used by DCADRE
	 */
	private static void goto55q() {
		// 55 IF(ERRER.GT.ERGOAL.AND.(ERGL+ERRER).NE.ERGL) GO TO 175;
		if ((ERRERq > ERGOALq) && ((ERGLq + ERRERq) != ERGLq)) {
			goto175q();
			if (exitq)
				return;
		}
		goto145q();
		if (exitq)
			return;
	}

	/**
	 * Used by DCADRE
	 */
	private static void goto60q() {
		// 60 IF(H2CONV) GO TO 65;
		if (H2CONVq) {
			goto65q();
			if (exitq)
				return;
		}
		AITKENq = false;
		H2CONVq = true;
		// 65 FEXTRP=FOUR;
		FEXTRPq = FOURq;
		// 70 IT=IT + 1;
		ITq = ITq + 1;
		// VINTq=STEPq*Tq(L,IT);
		VINTq = STEPq * Tq[Lq - 1][ITq - 1];
		// ERRERq=DABS(STEPq/(FEXTRPq-ONEq)*T(IT-1,L));
		ERRERq = DABS(STEPq / (FEXTRPq - ONEq) * Tq[ITq - 2][Lq - 1]);
		// IF(ERRER.LE.ERGOAL) GO TO 160;
		if (ERRERq <= ERGOALq) {
			goto160q();
			if (exitq)
				return;
		}
		// IF(ERGL+ERRER.EQ.ERGL) GO TO 160;
		if ((ERGLq + ERRERq) == ERGLq) {
			goto160q();
			if (exitq)
				return;
		}
		// IF(IT.EQ.LM1) GO TO 125;
		if (ITq == LM1q) {
			goto125q();
			if (exitq)
				return;
		}
		// IF(T(IT,LM1).EQ.ZERO) GO TO 70;
		if (Tq[ITq - 1][LM1q - 1] == ZEROq) {
			goto70q();
			if (exitq)
				return;
		}
		// IF(T(IT,LM1).LE.FEXTRP) GO TO 125;
		if (Tq[ITq - 1][LM1q - 1] <= FEXTRPq) {
			goto125q();
			if (exitq)
				return;
		}
		// IF(DABS(T(IT,LM1)/FOUR-FEXTRP)/FEXTRP.LT.AITTOL)
		if (DABS(Tq[ITq - 1][LM1q - 1] / FOURq - FEXTRPq) / FEXTRPq < AITTOLq)
			FEXTRPq = FEXTRPq * FOURq;
		goto70q();
		if (exitq)
			return;
	}

	/**
	 * Used by DCADRE
	 */
	private static void goto65q() {
		// 65 FEXTRP=FOUR;
		FEXTRPq = FOURq;
		// 70 IT=IT + 1;
		ITq = ITq + 1;
		// VINTq=STEPq*Tq(L,IT);
		VINTq = STEPq * Tq[Lq - 1][ITq - 1];
		// ERRERq=DABS(STEPq/(FEXTRPq-ONEq)*T(IT-1,L));
		ERRERq = DABS(STEPq / (FEXTRPq - ONEq) * Tq[ITq - 2][Lq - 1]);
		// IF(ERRER.LE.ERGOAL) GO TO 160;
		if (ERRERq <= ERGOALq) {
			goto160q();
			if (exitq)
				return;
		}
		// IF(ERGL+ERRER.EQ.ERGL) GO TO 160;
		if ((ERGLq + ERRERq) == ERGLq) {
			goto160q();
			if (exitq)
				return;
		}
		// IF(IT.EQ.LM1) GO TO 125;
		if (ITq == LM1q) {
			goto125q();
			if (exitq)
				return;
		}
		// IF(T(IT,LM1).EQ.ZERO) GO TO 70;
		if (Tq[ITq - 1][LM1q - 1] == ZEROq) {
			goto70q();
			if (exitq)
				return;
		}
		// IF(T(IT,LM1).LE.FEXTRP) GO TO 125;
		if (Tq[ITq - 1][LM1q - 1] <= FEXTRPq) {
			goto125q();
			if (exitq)
				return;
		}
		// IF(DABS(T(IT,LM1)/FOUR-FEXTRP)/FEXTRP.LT.AITTOL)
		if (DABS(Tq[ITq - 1][LM1q - 1] / FOURq - FEXTRPq) / FEXTRPq < AITTOLq)
			FEXTRPq = FEXTRPq * FOURq;
		goto70q();
		if (exitq)
			return;
	}

	/**
	 * Used by DCADRE
	 */
	private static void goto70q() {
		// 70 IT=IT + 1;
		ITq = ITq + 1;
		// VINTq=STEPq*Tq(L,IT);
		VINTq = STEPq * Tq[Lq - 1][ITq - 1];
		// ERRERq=DABS(STEPq/(FEXTRPq-ONEq)*T(IT-1,L));
		ERRERq = DABS(STEPq / (FEXTRPq - ONEq) * Tq[ITq - 2][Lq - 1]);
		// IF(ERRER.LE.ERGOAL) GO TO 160;
		if (ERRERq <= ERGOALq) {
			goto160q();
			if (exitq)
				return;
		}
		// IF(ERGL+ERRER.EQ.ERGL) GO TO 160;
		if ((ERGLq + ERRERq) == ERGLq) {
			goto160q();
			if (exitq)
				return;
		}
		// IF(IT.EQ.LM1) GO TO 125;
		if (ITq == LM1q) {
			goto125q();
			if (exitq)
				return;
		}
		// IF(T(IT,LM1).EQ.ZERO) GO TO 70;
		if (Tq[ITq - 1][LM1q - 1] == ZEROq) {
			goto70q();
			if (exitq)
				return;
		}
		// IF(T(IT,LM1).LE.FEXTRP) GO TO 125;
		if (Tq[ITq - 1][LM1q - 1] <= FEXTRPq) {
			goto125q();
			if (exitq)
				return;
		}
		// IF(DABS(T(IT,LM1)/FOUR-FEXTRP)/FEXTRP.LT.AITTOL)
		if (DABS(Tq[ITq - 1][LM1q - 1] / FOURq - FEXTRPq) / FEXTRPq < AITTOLq)
			FEXTRPq = FEXTRPq * FOURq;
		goto70q();
		if (exitq)
			return;
	}

	/**
	 * Used by DCADRE
	 */
	private static void goto75q() {
		// "                              INTEGRAND MAY HAVE X**ALPHA TYPE    "
		// "                                SINGULARITY                       "
		// "                                RESULTING IN A RATIO OF *SING*  = "
		// "                                2**(ALPHA + 1)                    "
		// IF(T(1,LM1).LT.AITLOW) GO TO 175;
		if (Tq[0][LM1q - 1] < AITLOWq) {
			goto175q();
			if (exitq)
				return;
		}
		if (AITKENq) {
			goto80q();
			if (exitq)
				return;
		}
		H2CONVq = false;
		AITKENq = true;
		FEXTRPq = Tq[Lq - 3][LM1q - 1];// 80 FEXTRP=T(L-2,LM1);
		if (FEXTRPq > FOURP5q) {
			goto65q();
			if (exitq)
				return;
		}
		if (FEXTRPq < AITLOWq) {
			goto175q();
			if (exitq)
				return;
		}
		// if(DABS(FEXTRP-T(L-3,LM1))/T(1,LM1).GT.H2TOL) GO TO 175;
		if (DABS(FEXTRPq - Tq[Lq - 4][LM1q - 1]) / Tq[0][LM1q - 1] > H2TOLq) {
			goto175q();
			if (exitq)
				return;
		}
		SINGq = FEXTRPq;
		FEXTM1q = ONEq / (FEXTRPq - ONEq);
		AITq[0] = ZEROq;// AIT(1)=ZERO;
		for (Iq = 1; Iq < Lq; Iq++)// DO 85 I=2,L;
		{
			// AIT(I)=T(I,1) + (T(I,1)-T(I-1,1))*FEXTM1;
			AITq[Iq] = Tq[Iq][0] + (Tq[Iq][0] - Tq[Iq - 1][0]) * FEXTM1q;
			// Rq[Iq]=Tq(1,I-1);
			Rq[Iq] = Tq[0][Iq - 1];
			// DIF(I)=AIT(I) - AIT(I-1);
			DIFq[Iq] = AITq[Iq] - AITq[Iq - 1];
		}
		// 85 CONTINUE;
		ITq = 2;
		VINTq = STEPq * AITq[Lq - 1];// 90 VINT=STEP*AIT(L);
		ERRERq = ERRERq * FEXTM1q;
		// IF(ERRER.GT.ERGOAL.AND.(ERGL+ERRER).NE.ERGL) GO TO 95;
		if ((ERRERq > ERGOALq) && ((ERGLq + ERRERq) != ERGLq)) {
			goto95q();
			if (exitq)
				return;
		}
		ALPHAq = DLOG10(SINGq) / ALG4O2q - ONEq;
		IERq = Math.max(IERq, 65);
		goto160q();
		if (exitq)
			return;
	}

	/**
	 * Used by DCADRE
	 */
	private static void goto80q() {
		// 80 FEXTRP=T(L-2,LM1);
		FEXTRPq = Tq[Lq - 3][LM1q - 1];
		if (FEXTRPq > FOURP5q) {
			goto65q();
			if (exitq)
				return;
		}
		// IF(FEXTRP.LT.AITLOW) GO TO 175;
		if (FEXTRPq < AITLOWq) {
			goto175q();
			if (exitq)
				return;
		}
		// IF(DABS(FEXTRP-T(L-3,LM1))/T(1,LM1).GT.H2TOL) GO TO 175;
		if (DABS(FEXTRPq - Tq[Lq - 4][LM1q - 1]) / Tq[0][LM1q - 1] > H2TOLq) {
			goto175q();
			if (exitq)
				return;
		}
		SINGq = FEXTRPq;
		FEXTM1q = ONEq / (FEXTRPq - ONEq);
		// AIT(1)=ZERO;
		AITq[0] = ZEROq;
		for (Iq = 1; Iq < Lq; Iq++)// DO 85 I=2,L;
		{
			AITq[Iq] = Tq[Iq][0] + (Tq[Iq][0] - Tq[Iq - 1][0]) * FEXTM1q;
			// AIT(I)=T(I,1) + (T(I,1)-T(I-1,1))*FEXTM1;
			Rq[Iq] = Tq[0][Iq - 1];
			// R(I)=T(1,I-1);
			DIFq[Iq] = AITq[Iq] - AITq[Iq - 1];
			// DIF(I)=AIT(I) - AIT(I-1);
		}
		// 85 CONTINUE;
		ITq = 2;
		// 90 VINT=STEP*AIT(L);
		VINTq = STEPq * AITq[Lq - 1];
		ERRERq = ERRERq * FEXTM1q;
		// IF(ERRER.GT.ERGOAL.AND.(ERGL+ERRER).NE.ERGL) GO TO 95;
		if ((ERRERq > ERGOALq) && ((ERGLq + ERRERq) != ERGLq)) {
			goto95q();
			if (exitq)
				return;
		}
		ALPHAq = DLOG10(SINGq) / ALG4O2q - ONEq;
		IERq = Math.max(IERq, 65);
		goto160q();
		if (exitq)
			return;
	}

	/**
	 * Used by DCADRE
	 */
	private static void goto90q() {
		// 90 VINT=STEP*AIT(L);
		VINTq = STEPq * AITq[Lq - 1];
		ERRERq = ERRERq * FEXTM1q;
		// IF(ERRER.GT.ERGOAL.AND.(ERGL+ERRER).NE.ERGL) GO TO 95;
		if ((ERRERq > ERGOALq) && ((ERGLq + ERRERq) != ERGLq)) {
			goto95q();
			if (exitq)
				return;
		}
		ALPHAq = DLOG10(SINGq) / ALG4O2q - ONEq;
		IERq = Math.max(IERq, 65);
		goto160q();
		if (exitq)
			return;
	}

	/**
	 * Used by DCADRE
	 */
	private static void goto95q() {
		ITq = ITq + 1;// 95 IT=IT + 1;
		// IF(IT.EQ.LM1) GO TO 125;
		if (ITq == LM1q) {
			goto125q();
			if (exitq)
				return;
		}
		if (ITq > 3) {
			goto100q();
			if (exitq)
				return;
		}
		H2NXTq = FOURq;
		SINGNXq = SINGq + SINGq;
		// 100 IF(H2NXT.LT.SINGNX) GO TO 105;
		if (H2NXTq < SINGNXq) {
			goto105q();
			if (exitq)
				return;
		}
		FEXTRPq = SINGNXq;
		SINGNXq = SINGNXq + SINGNXq;
		goto110q();
		if (exitq)
			return;
	}

	/**
	 * Used by DCADRE
	 */
	private static void goto100q() {
		// 100 IF(H2NXT.LT.SINGNX) GO TO 105;
		if (H2NXTq < SINGNXq) {
			goto105q();
			if (exitq)
				return;
		}
		FEXTRPq = SINGNXq;
		SINGNXq = SINGNXq + SINGNXq;
		goto110q();
		if (exitq)
			return;
	}

	/**
	 * Used by DCADRE
	 */
	private static void goto105q() {
		FEXTRPq = H2NXTq;
		H2NXTq = FOURq * H2NXTq;
		// 110 DO 115 I=IT,LM1;
		for (Iq = ITq - 1; Iq < LM1q; Iq++) {
			Rq[Iq + 1] = ZEROq;
			// IF(TABTLM+DABS(DIF(I+1)).NE.TABTLM) R(I+1)=DIF(I)/DIF(I+1);
			if (TABTLMq + DABS(DIFq[Iq + 1]) != TABTLMq) {
				Rq[Iq + 1] = DIFq[Iq] / DIFq[Iq + 1];
			}
		}
		// 115 CONTINUE;
		H2TFEXq = -H2TOLq * FEXTRPq;
		// IF(R(L)-FEXTRP.LT.H2TFEX) GO TO 125;
		if ((Rq[Lq - 1] - FEXTRPq) < H2TFEXq) {
			goto125q();
			if (exitq)
				return;
		}
		// IF(R(L-1)-FEXTRP.LT.H2TFEX) GO TO 125;
		if (Rq[Lq - 2] - FEXTRPq < H2TFEXq) {
			goto125q();
			if (exitq)
				return;
		}
		ERRERq = ASTEPq * DABS(DIFq[Lq - 1]);
		FEXTM1q = ONEq / (FEXTRPq - ONEq);
		// DO 120 I=IT,L;
		for (Iq = ITq - 1; Iq < Lq; Iq++) {
			// AIT(I)=AIT(I) + DIF(I)*FEXTM1;
			AITq[Iq] = AITq[Iq] + DIFq[Iq] * FEXTM1q;
			// DIF(I)=AIT(I) - AIT(I-1);
			DIFq[Iq] = AITq[Iq] - AITq[Iq - 1];
		}
		// 120 CONTINUE;
		goto90q();
		if (exitq)
			return;
	}

	/**
	 * Used by DCADRE
	 */
	private static void goto110q() {
		// 110 DO 115 I=IT,LM1;
		for (Iq = ITq - 1; Iq < LM1q; Iq++) {
			Rq[Iq + 1] = ZEROq;
			// IF(TABTLM+DABS(DIF(I+1)).NE.TABTLM) R(I+1)=DIF(I)/DIF(I+1);
			if (TABTLMq + DABS(DIFq[Iq + 1]) != TABTLMq) {
				Rq[Iq + 1] = DIFq[Iq] / DIFq[Iq + 1];
			}
		}
		// 115 CONTINUE;
		H2TFEXq = -H2TOLq * FEXTRPq;
		// IF(R(L)-FEXTRP.LT.H2TFEX) GO TO 125;
		if ((Rq[Lq - 1] - FEXTRPq) < H2TFEXq) {
			goto125q();
			if (exitq)
				return;
		}
		// IF(R(L-1)-FEXTRP.LT.H2TFEX) GO TO 125;
		if (Rq[Lq - 2] - FEXTRPq < H2TFEXq) {
			goto125q();
			if (exitq)
				return;
		}
		ERRERq = ASTEPq * DABS(DIFq[Lq - 1]);
		FEXTM1q = ONEq / (FEXTRPq - ONEq);
		// DO 120 I=IT,L;
		for (Iq = ITq - 1; Iq < Lq; Iq++) {
			// AIT(I)=AIT(I) + DIF(I)*FEXTM1;
			AITq[Iq] = AITq[Iq] + DIFq[Iq] * FEXTM1q;
			// DIF(I)=AIT(I) - AIT(I-1);
			DIFq[Iq] = AITq[Iq] - AITq[Iq - 1];
		}
		// 120 CONTINUE;
		goto90q();
		if (exitq)
			return;
	}

	/**
	 * Used by DCADRE
	 */
	private static void goto125q() {
		// "                              CURRENT TRAPEZOID SUM AND RESULTING "
		// "                                EXTRAPOLATED VALUES DID NOT GIVE  "
		// "                                A SMALL ENOUGH *ERRER*.           "
		// "                                NOTE -- HAVING PREVER .LT. ERRER  "
		// "                                IS AN ALMOST CERTAIN SIGN OF      "
		// "                                BEGINNING TROUBLE WITH IN THE FUNC"
		// "                                TION VALUES. HENCE, A WATCH FOR,  "
		// "                                AND CONTROL OF, NOISE SHOULD      "
		// "                                BEGIN HERE.                       "

		// 125 FEXTRP=DMAX1(PREVER/ERRER,AITLOW);
		FEXTRPq = DMAX1(PREVERq / ERRERq, AITLOWq);
		PREVERq = ERRERq;
		// IF(L.LT.5) GO TO 15;
		if (Lq < 5) {
			goto15q();
			if (exitq)
				return;
		}
		// IF(L-IT.GT.2.AND.ISTAGE.LT.MXSTGE) GO TO 170;
		if (((Lq - ITq) > 2) && (ISTAGEq < MXSTGEq)) {
			goto170q();
			if (exitq)
				return;
		}
		ERRETq = ERRERq / Math.pow(FEXTRPq, (MAXTBLq - Lq));
		// IF(ERRET.GT.ERGOAL.AND.(ERGL+ERRET).NE.ERGL) GO TO 170;
		if ((ERRETq > ERGOALq) && ((ERGLq + ERRETq) != ERGLq)) {
			goto170q();
			if (exitq)
				return;
		}
		goto15q();
		if (exitq)
			return;

	}

	/**
	 * Used by DCADRE
	 */
	private static void goto130q() {
		// 130 IF(ERRER.GT.ERGOAL.AND.(ERGL+ERRER).NE.ERGL) GO TO 170;
		// "   NOTE THAT  2*FN=2**L              "
		if ((ERRERq > ERGOALq) && ((ERGLq + ERRERq) != ERGLq)) {
			goto170q();
			if (exitq)
				return;

		}
		DIFFq = DABS(Tq[0][Lq - 1]) * (FNq + FNq);// DIFFq=DABS(T(1,L))*(FN+FN);
		goto160q();
		if (exitq)
			return;
	}

	/**
	 * Used by DCADRE
	 */
	private static void goto135q() {
		// "                              INTEGRAND IS STRAIGHT LINE          "
		// "                                TEST THIS ASSUMPTION BY COMPARING "
		// "                                THE VALUE OF THE INTEGRAND AT     "
		// "                                FOUR *RANDOMLY CHOSEN* POINTS WITH"
		// "                                THE VALUE OF THE STRAIGHT LINE    "
		// "                                INTERPOLATING THE INTEGRAND AT THE"
		// "                                TWO END POINTS OF THE SUB-INTERVAL"
		// "                                IF TEST IS PASSED, ACCEPT *VINT*  "

		SLOPEq = (FENDq - FBEGq) * TWOq;// 135 SLOPE=(FEND-FBEG)*TWO;
		FBEG2q = FBEGq + FBEGq;
		for (Iq = 0; Iq < 4; Iq++)// DO 140 I=1,4;Fa pana la 140!!
		{
			RVALq = BEGq + RNq[Iq] * STEPq;
			DIFFq = DABS(F(RVALq) - FBEG2q - RNq[Iq] * SLOPEq);
			if ((TABTLMq + DIFFq) != TABTLMq)// GO TO 155;
			{
				goto155q();
				if (exitq)
					return;
			}
		}
		// 140 CONTINUE;
		goto160q();
		if (exitq)
			return;
	}

	/**
	 * Used by DCADRE
	 */
	private static void goto145q() {
		// "   NOISE MAY BE DOMINANT FEATURE       "
		// "   ESTIMATE NOISE LEVEL BY COMPARING "
		// "   THE VALUE OF THE INTEGRAND AT     "
		// "   FOUR *RANDOMLY CHOSEN* POINTS WITH"
		// "   THE VALUE OF THE STRAIGHT LINE    "
		// "   INTERPOLATING THE INTEGRAND AT THE"
		// "   TWO ENDPOINTS. IF SMALL ENOUGH,   "
		// "   ACCEPT *VINT*                     "
		// "   INTERGRATION OVER CURRENT SUB-      "
		// "   INTERVAL SUCCESSFUL               "
		// "   ADD *VINT* TO *DCADRE* AND *ERRER*"
		// "   TO *ERROR*, THEN SET UP NEXT SUB- "
		// "   INTERVAL, IF ANY.                 "

		SLOPEq = (FENDq - FBEGq) * TWOq;// 145 SLOPE=(FEND-FBEG)*TWO;
		FBEG2q = FBEGq + FBEGq;
		Iq = 0;// Iq=1;
		RVALq = BEGq + RNq[Iq] * STEPq;// 150
		DIFFq = DABS(F(RVALq) - FBEG2q - RNq[Iq] * SLOPEq);// 155 next
		ERRERq = DMAX1(ERRERq, ASTEPq * DIFFq);// 155
												// ERRER=DMAX1(ERRER,ASTEP*DIFF);
		if ((ERRERq > ERGOALq) && ((ERGLq + ERRERq) != ERGLq)) {
			goto175q();// GO TO 175;
			if (exitq)
				return;
		}
		Iq = Iq + 1;
		if (Iq <= 3) // IF(I<=4) GO TO 150;
		{
			goto150q();
			if (exitq)
				return;
		}
		IERq = 66;
		// "                              INTERGRATION OVER CURRENT SUB-      "
		// "                                INTERVAL SUCCESSFUL               "
		// "                                ADD *VINT* TO *DCADRE* AND *ERRER*"
		// "                                TO *ERROR*, THEN SET UP NEXT SUB- "
		// "                                INTERVAL, IF ANY.                 "
		CADREq = CADREq + VINTq;// 160 CADRE=CADRE + VINT;
		ERRORq = ERRORq + ERRERq;
		if (RIGHTq)// IF(RIGHT) GO TO 165;
		{
			goto165q();
			if (exitq)
				return;
		}
		ISTAGEq = ISTAGEq - 1;
		if (ISTAGEq == 0) // IF(ISTAGE.EQ.0) GO TO 220;
		{
			IERq = 131; // 200 IER=131;GO TO 215;
			goto215q();
			exitq = true;
			return;
		}
		REGLARq = REGLSVq[ISTAGEq - 1];// REGLAR=REGLSV[ISTAGE];
		BEGq = BEGINq[ISTAGEq - 1];// BEG=BEGIN(ISTAGE);
		ENDq = FINISq[ISTAGEq - 1];// END=FINIS(ISTAGE);
		CURESTq = CURESTq - ESTq[ISTAGEq] + VINTq;// CUREST=CUREST -
													// EST(ISTAGE+1) + VINT;
		IENDq = IBEGq - 1;
		FENDq = TSq[IENDq - 1];// FEND=TS(IEND);
		IBEGq = IBEGSq[ISTAGEq - 1];// IBEG=IBEGS(ISTAGE);
		goto180q();
		if (exitq)
			return;
	}

	/**
	 * Used by DCADRE
	 */
	private static void goto150q() {
		RVALq = BEGq + RNq[Iq] * STEPq;
		DIFFq = DABS(F(RVALq) - FBEG2q - RNq[Iq] * SLOPEq);// 155 next
		ERRERq = DMAX1(ERRERq, ASTEPq * DIFFq);// 155
												// ERRER=DMAX1(ERRER,ASTEP*DIFF);
		if ((ERRERq > ERGOALq) && ((ERGLq + ERRERq) != ERGLq)) {
			goto175q();// GO TO 175;
			if (exitq)
				return;
		}
		Iq = Iq + 1;
		if (Iq <= 3) // IF(I<=4) GO TO 150;
		{
			goto150q();
			if (exitq)
				return;
		}
		IERq = 66;
		// "                              INTERGRATION OVER CURRENT SUB-      "
		// "                                INTERVAL SUCCESSFUL               "
		// "                                ADD *VINT* TO *DCADRE* AND *ERRER*"
		// "                                TO *ERROR*, THEN SET UP NEXT SUB- "
		// "                                INTERVAL, IF ANY.                 "
		CADREq = CADREq + VINTq;// 160 CADRE=CADRE + VINT;
		ERRORq = ERRORq + ERRERq;
		if (RIGHTq)// IF(RIGHT) GO TO 165;
		{
			goto165q();
			if (exitq)
				return;
		}
		ISTAGEq = ISTAGEq - 1;
		if (ISTAGEq == 0) // IF(ISTAGE.EQ.0) GO TO 220;
		{
			IERq = 131; // 200 IER=131;GO TO 215;
			goto215q();
			exitq = true;
			return;
		}
		REGLARq = REGLSVq[ISTAGEq - 1];// REGLAR=REGLSV[ISTAGE];
		BEGq = BEGINq[ISTAGEq - 1];// BEG=BEGIN(ISTAGE);
		ENDq = FINISq[ISTAGEq - 1];// END=FINIS(ISTAGE);
		CURESTq = CURESTq - ESTq[ISTAGEq] + VINTq;// CUREST=CUREST -
													// EST(ISTAGE+1) + VINT;
		IENDq = IBEGq - 1;
		FENDq = TSq[IENDq - 1];// FEND=TS(IEND);
		IBEGq = IBEGSq[ISTAGEq - 1];// IBEG=IBEGS(ISTAGE);
		goto180q();
		if (exitq)
			return;
	}

	/**
	 * Used by DCADRE
	 */
	private static void goto155q() {
		ERRERq = DMAX1(ERRERq, ASTEPq * DIFFq);// 155
												// ERRER=DMAX1(ERRER,ASTEP*DIFF);
		if ((ERRERq > ERGOALq) && ((ERGLq + ERRERq) != ERGLq)) {
			goto175q();// GO TO 175;
			if (exitq)
				return;
		}
		Iq = Iq + 1;
		if (Iq <= 3) // IF(I<=4) GO TO 150;
		{
			goto150q();
			if (exitq)
				return;
		}
		IERq = 66;
		// "                              INTERGRATION OVER CURRENT SUB-      "
		// "                                INTERVAL SUCCESSFUL               "
		// "                                ADD *VINT* TO *DCADRE* AND *ERRER*"
		// "                                TO *ERROR*, THEN SET UP NEXT SUB- "
		// "                                INTERVAL, IF ANY.                 "
		CADREq = CADREq + VINTq;// 160 CADRE=CADRE + VINT;
		ERRORq = ERRORq + ERRERq;
		if (RIGHTq)// IF(RIGHT) GO TO 165;
		{
			goto165q();
			if (exitq)
				return;
		}
		ISTAGEq = ISTAGEq - 1;
		if (ISTAGEq == 0) // IF(ISTAGE.EQ.0) GO TO 220;
		{
			IERq = 131; // 200 IER=131;GO TO 215;
			goto215q();
			exitq = true;
			return;
		}
		REGLARq = REGLSVq[ISTAGEq - 1];// REGLAR=REGLSV[ISTAGE];
		BEGq = BEGINq[ISTAGEq - 1];// BEG=BEGIN(ISTAGE);
		ENDq = FINISq[ISTAGEq - 1];// END=FINIS(ISTAGE);
		CURESTq = CURESTq - ESTq[ISTAGEq] + VINTq;// CUREST=CUREST -
													// EST(ISTAGE+1) + VINT;
		IENDq = IBEGq - 1;
		FENDq = TSq[IENDq - 1];// FEND=TS(IEND);
		IBEGq = IBEGSq[ISTAGEq - 1];// IBEG=IBEGS(ISTAGE);
		goto180q();
		if (exitq)
			return;
	}

	/**
	 * Used by DCADRE
	 */
	private static void goto160q() {
		// "                              INTERGRATION OVER CURRENT SUB-      "
		// "                                INTERVAL SUCCESSFUL               "
		// "                                ADD *VINT* TO *DCADRE* AND *ERRER*"
		// "                                TO *ERROR*, THEN SET UP NEXT SUB- "
		// "                                INTERVAL, IF ANY.                 "
		CADREq = CADREq + VINTq;// 160 CADRE=CADRE + VINT;
		ERRORq = ERRORq + ERRERq;
		if (RIGHTq) {
			goto165q();
			if (exitq)
				return;
		}
		ISTAGEq = ISTAGEq - 1;
		if (ISTAGEq == 0) // IF(ISTAGE.EQ.0) GO TO 220;
		{
			IERq = 131; // 200 IER=131;GO TO 215;
			goto215q();
			exitq = true;
			return;
		}
		REGLARq = REGLSVq[ISTAGEq - 1];// REGLAR=REGLSV[ISTAGE];
		BEGq = BEGINq[ISTAGEq - 1];// BEG=BEGIN(ISTAGE);
		ENDq = FINISq[ISTAGEq - 1];// END=FINIS(ISTAGE);
		CURESTq = CURESTq - ESTq[ISTAGEq] + VINTq;// CUREST=CUREST -
													// EST(ISTAGE+1) + VINT;
		IENDq = IBEGq - 1;
		FENDq = TSq[IENDq - 1];// FEND=TS(IEND);
		IBEGq = IBEGSq[ISTAGEq - 1];// IBEG=IBEGS(ISTAGE);
		goto180q();
		if (exitq)
			return;
	}

	/**
	 * Used by DCADRE
	 */
	private static void goto165q() {
		CURESTq = CURESTq + VINTq;// 165 CUREST=CUREST + VINT;
		STAGEq = STAGEq + STAGEq;
		IENDq = IBEGq;
		IBEGq = IBEGSq[ISTAGEq - 1];// IBEG=IBEGS(ISTAGE);
		ENDq = BEGq;
		BEGq = BEGINq[ISTAGEq - 1];// BEG=BEGIN(ISTAGE);
		FENDq = FBEGq;
		FBEGq = TSq[IBEGq - 1];// FBEG=TS(IBEG);
		goto5q();// @@@@@
		if (exitq)
			return;
	}

	/**
	 * Used by DCADRE
	 */
	private static void goto170q() {
		// "                              INTEGRATION OVER CURRENT SUBINTERVAL"
		// "                                IS UNSUCCESSFUL. MARK SUBINTERVAL "
		// "                                FOR FURTHER SUBDIVISION. SET UP   "
		// "                                NEXT SUBINTERVAL.                 "
		REGLARq = true;
		if (ISTAGEq == MXSTGEq)// 175 IF(ISTAGE.EQ.MXSTGE) GO TO 205;
		{
			IERq = 132;// 205 IER=132;
			goto215q();
			exitq = true;
			return;
		}
		if (RIGHTq) {
			goto185q();
			if (exitq)
				return;
		}
		REGLSVq[ISTAGEq] = REGLARq;// REGLSV(ISTAGE+1)=REGLAR;
		BEGINq[ISTAGEq - 1] = BEGq;// BEGIN(ISTAGE)=BEG;
		IBEGSq[ISTAGEq - 1] = IBEGq;// IBEGS(ISTAGE)=IBEG;
		STAGEq = STAGEq * HALFq;
		RIGHTq = true;// 180 RIGHT=.TRUE.;
		BEGq = (BEGq + ENDq) * HALFq;
		IBEGq = (IBEGq + IENDq) / 2;
		TSq[IBEGq - 1] = TSq[IBEGq - 1] * HALFq;// TS(IBEG)=TS(IBEG)*HALF;
		FBEGq = TSq[IBEGq - 1];// FBEG=TS(IBEG);
		goto10q();// @@;
		if (exitq)
			return;
	}

	/**
	 * Used by DCADRE
	 */
	private static void goto175q() {
		if (ISTAGEq == MXSTGEq)// 175 IF(ISTAGE.EQ.MXSTGE) GO TO 205;
		{
			IERq = 132;// 205 IER=132;
			goto215q();
			exitq = true;
			return;
		}
		if (RIGHTq) {
			goto185q();
			if (exitq)
				return;
		}
		REGLSVq[ISTAGEq] = REGLARq;// REGLSV(ISTAGE+1)=REGLAR;
		BEGINq[ISTAGEq - 1] = BEGq;// BEGIN(ISTAGE)=BEG;
		IBEGSq[ISTAGEq - 1] = IBEGq;// IBEGS(ISTAGE)=IBEG;
		STAGEq = STAGEq * HALFq;
		RIGHTq = true;// 180 RIGHT=.TRUE.;
		BEGq = (BEGq + ENDq) * HALFq;
		IBEGq = (IBEGq + IENDq) / 2;
		TSq[IBEGq - 1] = TSq[IBEGq - 1] * HALFq;// TS(IBEG)=TS(IBEG)*HALF;
		FBEGq = TSq[IBEGq - 1];// FBEG=TS(IBEG);
		goto10q();// @@;
		if (exitq)
			return;
	}

	/**
	 * Used by DCADRE
	 */
	private static void goto180q() {
		RIGHTq = true;// 180 RIGHT=.TRUE.;
		BEGq = (BEGq + ENDq) * HALFq;
		IBEGq = (IBEGq + IENDq) / 2;
		TSq[IBEGq - 1] = TSq[IBEGq - 1] * HALFq;// TS(IBEG)=TS(IBEG)*HALF;
		FBEGq = TSq[IBEGq - 1];// FBEG=TS(IBEG);
		goto10q();// @@@@
		if (exitq)
			return;
	}

	/**
	 * Used by DCADRE
	 */
	private static void goto185q() {
		NNLEFTq = IBEGq - IBEGSq[ISTAGEq - 1];// 185 NNLEFT=IBEG -
												// IBEGS(ISTAGE);
		if ((IENDq + NNLEFTq) >= MAXTSq)// IF(IEND+NNLEFT.GE.MAXTS) GO TO 200;
		{
			IERq = 131;// 200 IER=131;GO TO 215;
			goto215q();
			exitq = true;
			return;// 9005 RETURN;
		}
		IIIq = IBEGSq[ISTAGEq - 1];// III=IBEGS(ISTAGE);
		IIq = IENDq;
		for (int I1 = IIIq - 1; I1 < IBEGq; I1++)// DO 190
													// I=III,IBEG;^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
		{
			IIq = IIq + 1;
			TSq[IIq - 1] = TSq[I1];// TS(II)=TS(I);
		}// 190 CONTINUE;
		for (int I1 = IBEGq - 1; I1 < IIq; I1++)// DO 195
												// I=IBEG,II;II^^^^^^^^^^^^^OK
		{
			TSq[IIIq - 1] = TSq[I1];// TS(III)=TS(I);
			IIIq = IIIq + 1;
		}// 195 CONTINUE;
		IENDq = IENDq + 1;
		IBEGq = IENDq - NNLEFTq;
		FENDq = FBEGq;
		FBEGq = TSq[IBEGq - 1];// FBEG=TS(IBEG);
		FINISq[ISTAGEq - 1] = ENDq;// FINIS(ISTAGE)=END;
		ENDq = BEGq;
		BEGq = BEGINq[ISTAGEq - 1];// BEG=BEGIN(ISTAGE);
		BEGINq[ISTAGEq - 1] = ENDq;// BEGIN(ISTAGE)=END;
		REGLSVq[ISTAGEq - 1] = REGLARq;// REGLSV(ISTAGE)=REGLAR;
		ISTAGEq = ISTAGEq + 1;
		REGLARq = REGLSVq[ISTAGEq - 1];// REGLAR=REGLSV(ISTAGE);
		ESTq[ISTAGEq - 1] = VINTq;// EST(ISTAGE)=VINT;
		CURESTq = CURESTq + ESTq[ISTAGEq - 1];// CUREST=CUREST + EST(ISTAGE);
		goto5q();// @@;
		if (exitq)
			return;
	}

	/**
	 * Used by DCADRE
	 */
	private static void goto215q() {
		CADREq = CURESTq + VINTq;// 215 CADRE=CUREST + VINT;
		DCADREq = CADREq;// //220 DCADRE=CADRE;+//9000 CONTINUE;
	}

	/**
	 * Used by DCADRE
	 */
	private static void resetDCADRE() {
		DIFFq = 0.;
		H2CONVq = false;
		AITKENq = false;
		LM1q = 0;
		N2q = 0;
		FNq = 0.;
		ISTEPq = 0;
		IIq = 0;
		IIIq = 0;
		HOVNq = 0.;
		FIq = 0.;
		ISTEP2q = 0;
		SUMq = 0.;
		SUMABSq = 0.;
		ABSIq = 0.;
		ITq = 0;
		TABTLMq = 0.;
		ERGLq = 0.;
		ERGOALq = 0.;
		FEXTRPq = 0.;
		ERRERq = 0.;
		RIGHTq = false;
		STEPq = 0.;
		ASTEPq = 0.;
		TABSq = 0.;
		Lq = 0;
		Nq = 0;
		TSq = new double[2049];
		Tq = new double[10][10];
		ZEROq = 0.0;
		P1q = 0.1;
		HALFq = 0.5;
		ONEq = 1.0;
		TWOq = 2.0;
		FOURq = 4.0;
		FOURP5q = 4.5;
		TENq = 10.0;
		HUNq = 100.0;
		AITLOWq = 1.1;
		H2TOLq = 0.15;
		AITTOLq = 0.1;
		JUMPTLq = 0.01;
		MAXTSq = 2049;
		MAXTBLq = 10;
		MXSTGEq = 30;
		SLOPEq = 0.;
		FBEG2q = 0.;
		RNq = new double[4];
		RNq[0] = 0.7142005;
		RNq[1] = 0.3466282;
		RNq[2] = 0.843751;
		RNq[3] = 0.1263305;
		REGLSVq = new boolean[30];
		BEGINq = new double[30];
		FINISq = new double[30];
		ESTq = new double[30];
		IBEGSq = new int[30];
		NNLEFTq = 0;
		Iq = 0;
		// ---------------------------------
		DCADREq = 0.;
		ALG4O2q = DLOG10(TWOq);
		CADREq = ZEROq;
		CURESTq = ZEROq;
		VINTq = ZEROq;
		LENGTHq = 0.;
		ERRRq = 0.;
		ERRAq = 0.;
		STEPMNq = 0.;
		STEPNMq = 0.;
		BEGq = 0.;
		RVALq = 0.;
		FBEGq = 0.;
		STAGEq = HALFq;
		ISTAGEq = 1;
		FNSIZEq = ZEROq;
		PREVERq = ZEROq;
		REGLARq = false;
		FBEGq = 0.;
		IBEGq = 1;
		ENDq = 0.;
		FENDq = 0.;
		IENDq = 2;
		ERRORq = ZEROq;
		IERq = 0;
		LENGTHq = DABS(Bq - Aq);
		ERRRq = RERRq;
		ERRAq = DABS(AERRq);
		STEPMNq = (LENGTHq / Math.pow(2.0, MXSTGEq));
		STEPNMq = DMAX1(LENGTHq, DABS(Aq), DABS(Bq)) * TENq;
		// " THE GIVEN INTERVAL OF INTEGRATION   "
		// " IS THE FIRST INTERVAL CONSIDERED. "
		BEGq = Aq;
		RVALq = BEGq;
		FBEGq = F(RVALq) * HALFq;// @@
		TSq[0] = FBEGq;// (1)
		ENDq = Bq;
		RVALq = ENDq;
		FENDq = F(RVALq) * HALFq;// @@
		TSq[1] = FENDq;// (2)
		SINGq = 0.;
		FEXTM1q = 0.;
		Rq = new double[10];
		AITq = new double[10];
		DIFq = new double[10];
		ALPHAq = 0.;
		H2NXTq = 0.;
		SINGNXq = 0.;
		ERRETq = 0.;
		H2TFEXq = 0.;
		// -----------------------
		exitq = false;
	}

	// max from 3 var
	/**
	 * Maximum of 3 numbers
	 * @param a a
	 * @param b b
	 * @param c c
	 * @return the result
	 */
	public static double DMAX1(double a, double b, double c) {
		double r = Math.max(a, b);
		r = Math.max(r, c);
		return r;
	}

	// max from 2 var
	/**
	 * Maximum of 2 numbers
	 * @param a a
	 * @param b b
	 * @return the result
	 */
	public static double DMAX1(double a, double b) {
		double r = Math.max(a, b);
		return r;
	}

	// compute log in 10 base
	// a**(LOGaX)=X due to definition of logaritm.
	// the logarithm of a power equals the exponent multiplied with the
	// logarithm of the base
	// LOGbX=LOGb[a**(LOGaX)]=LOGaX*LOGba!!!
	// so LOG10X=LnX/Ln10
	/**
	 * Logarithm in base 10 of a number x
	 * @param x x
	 * @return the result
	 */
	public static double DLOG10(double x) {
		double result = Math.log(x);
		result = result / Math.log(10);
		return result;
	}

	// absolute value!!
	/**
	 * Absolute value of a number
	 * @param x x
	 * @return the result
	 */
	public static double DABS(double x) {
		double result = Math.abs(x);
		return result;
	}

	// return double -better than float- of int
	/**
	 * Return double representation of an int.
	 * @param x x
	 * @return the result
	 */
	public static double FLOAT(int x) {
		double result = x;
		return result;
	}

	/**
	 * Return the function (of X) to be integrated by QD.
	 * @param X X
	 * @return the result
	 */
	public static double F(double X) {
		double Y = 0.;
		if (funcname.compareTo("test") == 0)
			Y = test(X);
		else if (funcname.compareTo("BREMFZ") == 0)
			Y = BREMFZ(X);
		else if (funcname.compareTo("BRMSFZ") == 0)
			Y = BRMSFZ(X);
		else if (funcname.compareTo("PAIRFZ") == 0)
			Y = PAIRFZ(X);
		else if (funcname.compareTo("AFFACT") == 0)
			Y = AFFACT(X);
		else if (funcname.compareTo("BREMFR") == 0)
			Y = BREMFR(X);
		else if (funcname.compareTo("PAIRFR") == 0)
			Y = PAIRFR(X);

		return Y;
	}

	// END PWLF TO ELECTRON DATA SETS
	// --------------convertors to be independent with other class
	// file--------------------------
	/**
	 * Converts an int to another base
	 * @param value the string representation of an int (base 10 that is)
	 * @param base the desired base for conversion
	 * @return the String representation of the number in desired base
	 * @throws NumberFormatException can throw this exception
	 */
	public static String convertIntToOtherBase(String value, int base)
			throws NumberFormatException// 3a3f15f7=format Hexazecimal ex.
	{
		BigInteger bi = new BigInteger(value);
		return bi.toString(base);
	}

	/**
	 * Converts a number from other base to base 10
	 * @param value the String representation of a number
	 * @param base the base of the number
	 * @return the String representation of number in base 10
	 * @throws NumberFormatException can throw this exception
	 */
	public static String convertOtherBaseToInt(String value, int base)
			throws NumberFormatException {
		BigInteger bi = new BigInteger(value, base);
		return bi.toString(10);
	}

	/**
	 * Converts a string to int value.
	 * 
	 * @param value
	 *            the string
	 * @return the int value of input string
	 * @throws NumberFormatException
	 * can throws this exception
	 */
	public static int stringToInt(String value) throws NumberFormatException {
		value = value.trim();
		if (value.length() == 0) {
			return 0;
		} else {
			return Integer.parseInt(value);
		}
	}

	/**
	 * Converts an int number to string.
	 * 
	 * @param i
	 *            the int value
	 * @return the string representation
	 */
	public static String intToString(int i) {
		return Integer.toString(i);
	}

	/**
	 * Converts a string to float value.
	 * 
	 * @param value
	 *            the string
	 * @return the float value of input string
	 * @throws NumberFormatException
	 * can throws this exception
	 */
	public static float stringToFloat(String value)
			throws NumberFormatException {
		value = value.trim();
		if (value.length() == 0) {
			return 0;
		} else {
			return Float.parseFloat(value);
		}
	}

	/**
	 * Converts a float number to string.
	 * 
	 * @param f
	 *            the float value
	 * @return the string representation
	 */
	public static String floatToString(float f) {
		return Float.toString(f);
	}

	/**
	 * Converts a string to double value.
	 * 
	 * @param value
	 *            the string
	 * @return the double value of input string
	 * @throws NumberFormatException
	 * can throws this exception
	 */
	public static double stringToDouble(String value)
			throws NumberFormatException {
		value = value.trim();
		if (value.length() == 0) {
			return 0;
		} else {
			return Double.parseDouble(value);
		}
	}

	/**
	 * Converts a double number to string.
	 * 
	 * @param d
	 *            the double value
	 * @return the string representation
	 */
	public static String doubleToString(double d) {
		return Double.toString(d);
	}

	// --------------END convertors to be independent with other class
	// file--------------------------
	/**
	 * Return the filename (without extension) 
	 * @param file file
	 * @return the result
	 */
	public static String stripFileExtension(String file) {
		int idx = file.lastIndexOf(".");
		// handles unix hidden files and files without an extension.
		if (idx < 1) {
			return file;
		}
		return file.substring(0, idx);
	}

	/**
	 * Return the filename from file path.
	 * @param file the path
	 * @return the result
	 */
	public static String getFileName(String file) {
		int idx = file.lastIndexOf(file_sep) + 1;
		// handles unix hidden files and files without an extension.
		if (idx < 1) {
			return "";
		}
		return file.substring(idx);
	}

	// ----------------NOT USED YET FOUNCTION------------------
	/**
	 * DERIVATE OF ADFMOL-PROPORTIONAL TO APPROXIMATE P.D.F.
	 * @param X X
	 * @return the result
	 */
	public static double ADDMOL(double X) {
		// "     DERIVATE OF ADFMOL-PROPORTIONAL TO APPROXIMATE P.D.F.          "
		double ADDMOL = 1.0 / ((X - RM) * (X - RM));
		return ADDMOL;
	}

	/**
	 * APPROXIMATE C.D.F. FOR MOLLER AND BHABHA
	 * @param X X
	 * @return the result
	 */
	public static double ADFMOL(double X) {
		// "     APPROXIMATE C.D.F. FOR MOLLER AND BHABHA                       "
		double ADFMOL = -1.0 / (X - RM);
		return ADFMOL;
	}

	/**
	 * INVERSE OF ADFMOL
	 * @param X X
	 * @return the result
	 */
	public static double ADIMOL(double X) {
		// "     INVERSE OF ADFMOL                                              "
		double ADIMOL = -1.0 / X + RM;
		return ADIMOL;
	}

	/**
	 * FUNCTION TO CALL WHEN DERIVATIVE OF LOG IS DESIRED
	 * @param X X
	 * @return the result
	 */
	public static double AREC(double X) {
		// "     FUNCTION TO CALL WHEN DERIVATIVE OF LOG IS DESIRED              "
		double AREC = 1.0 / X;
		return AREC;
	}

	/**
	 * DIFFERENTIAL MOLLER CROSS SECTION FOR INCIDENT ELECTRON OF TOTAL 
	 * ENERGY EN0 TO PRODUCE SCATTERED ELECTRON OF TOTAL ENERGY EN.
	 * @param EN0 EN0
	 * @param EN EN
	 * @return the result
	 */
	public static double AMOLDM(double EN0, double EN) {
		// "***DIFFERENTIAL MOLLER CROSS SECTION FOR INCIDENT ELECTRON OF TOTAL"
		// "   ENERGY EN0 TO PRODUCE SCATTERED ELECTRON OF TOTAL ENERGY EN."
		T0 = EN0 - RM;
		double TM = T0 / RM;
		double EM = TM + 1.;// gamma
		C1 = (TM / EM) * (TM / EM);
		C2 = (2. * TM + 1.) / (EM * EM);
		double BETASQ = 1. - 1. / (EM * EM);
		CMOLL = RLC * EDEN * 2. * PI * R0 * R0 / (BETASQ * T0 * TM);// "CONSTANT FACTOR IN EQ.2.10.1"
		double AMOLDM = AMOLFM(EN);
		return AMOLDM;
	}

	/**
	 * Called by AMOLDM.
	 * @param EN EN
	 * @return the result
	 */
	public static double AMOLFM(double EN) {
		double T = EN - RM;
		double EPS = T / T0;
		double EPSP = 1. - EPS;
		double EPSI = 1. / EPS;
		double EPSPI = 1. / EPSP;
		double AMOLFM = CMOLL
				* (C1 + EPSI * (EPSI - C2) + EPSPI * (EPSPI - C2)); // "EQ.2.10.1"
		return AMOLFM;
	}

	/**
	 * DIFFERENTIAL TWO-PHOTON ANNIHILATION CROSS SECTION FOR POSITRON OF 
	 * INCIDENT TOTAL ENERGY E0 TO PRODUCE SECONDARY PHOTON OF ENERGY K.
	 * @param E0 E0
	 * @param K K
	 * @return the result
	 */
	public static double ANIHDM(double E0, double K) {
		// "***DIFFERENTIAL TWO-PHOTON ANNIHILATION CROSS SECTION FOR POSITRON OF"
		// "   INCIDENT TOTAL ENERGY E0 TO PRODUCE SECONDARY PHOTON OF ENERGY K."
		double GAM = E0 / RM; // "EQ.2.12.3"
		A = GAM + 1.; // "EQ.2.12.4"
		double T0P = GAM - 1.; // "EQ.2.12.5"
		C1 = RLC * EDEN * PI * R0 * R0 / (A * T0P * RM); // "EQ.2.12.9"
		C2 = A + 2.0 * GAM / A; // "2.12.10"
		double ANIHDM = ANIHFM(K);
		return ANIHDM;
	}

	/**
	 * Called by ANIHFM
	 * @param X X
	 * @return the result
	 */
	public static double S1(double X) {
		double S1 = C1 * (-1. + (C2 - 1.0 / X) / X); // "STATEMENT FUNCTION---EQ.2.12.8"
		return S1;
	}

	/**
	 * Called by ANIHRM
	 * @param X X
	 * @return the result
	 */
	public static double S2(double X) {
		double S2 = RM * C1 * (-X + C2 * Math.log(X) + 1.0 / X); // "STATEMENT FUNCTION---EQ.2.12.13"
		return S2;
	}

	/**
	 * Called by ANIHDM
	 * @param K K
	 * @return the result
	 */
	public static double ANIHFM(double K) {
		double KP = K / RM; // "EQ.2.12.7"
		double ANIHFM = S1(KP) + S1(A - KP); // "EQ.2.12.1"
		return ANIHFM;
	}

	/**
	 * TWO-PHOTON ANNIHILATION CROSS SECTION FOR INCIDENT POSITRON OF 
	 * TOTAL ENERGY E0 TO PRODUCE SECONDARY PHOTON IN THE ENERGY RANGE 
	 * K1 TO K2.
	 * @param E0 E0
	 * @param K1 K1
	 * @param K2 K3
	 * @return the result
	 */
	public static double ANIHRM(double E0, double K1, double K2) {
		// "***TWO-PHOTON ANNIHILATION CROSS SECTION FOR INCIDENT POSITRON OF"
		// "   TOTAL ENERGY E0 TO PRODUCE SECONDARY PHOTON IN THE ENERGY RANGE"
		// "   K1 TO K2."
		double GAM = E0 / RM; // "EQ.2.12.3"
		A = GAM + 1.; // "EQ.2.12.4"
		double T0P = GAM - 1.;// "EQ.2.12.5"
		C1 = RLC * EDEN * PI * R0 * R0 / (A * T0P * RM); // "EQ.2.12.9"
		C2 = A + 2. * GAM / A; // "EQ.2.12.10"
		// GAM=E0/RM; //"EQ.2.12.3"
		double KP1 = K1 / RM; // "EQ.2.12.12"
		double KP2 = K2 / RM; // "EQ.2.12.12"
		// A=GAM+1.; // "EQ.2.12.4"
		// double T0P=GAM-1.;// "EQ.2.12.5"
		// C1=RLC*EDEN*PI*R0**2/(A*T0P*RM); //"EQ.2.12.9"
		// C2=A+2.*GAM/A; //"EQ.2.12.10"
		double ANIHRM = S2(KP2) - S2(KP1) + S2(A - KP1) - S2(A - KP2); // "EQ.2.12.11"
		return ANIHRM;
	}

	/**
	 * DIFFERENTIAL BHABHA CROSS SECTION FOR INCIDENT POSITRON OF TOTAL 
	 * ENERGY EN0 TO PRODUCE SCATTERED ELECTRON OF TOTAL ENERGY EN.
	 * @param EN0 EN0
	 * @param EN EN
	 * @return the result
	 */
	public static double BHABDM(double EN0, double EN) {
		// "***DIFFERENTIAL BHABHA CROSS SECTION FOR INCIDENT POSITRON OF TOTAL"
		// "   ENERGY EN0 TO PRODUCE SCATTERED ELECTRON OF TOTAL ENERGY EN."
		double T0 = EN0 - RM;
		double TM = T0 / RM;
		double EM = TM + 1.;
		double Y = 1. / (TM + 2.);
		BETASI = 1. / (1. - 1. / (EM * EM));
		CBHAB = RLC * EDEN * 2. * PI * R0 * R0 / (T0 * TM); // "CONSTANT FACTOR IN EQ.2.11.1"
		B1 = 2. - Y * Y;
		B2 = 3. - Y * (6. - Y * (1. - Y * 2.));
		B3 = 2. - Y * (10. - Y * (16. - Y * 8.));
		B4 = 1. - Y * (6. - Y * (12. - Y * 8.));
		double BHABDM = BHABFM(EN);
		return BHABDM;

	}

	/**
	 * Called by BHABDM
	 * @param EN EN
	 * @return the result
	 */
	public static double BHABFM(double EN) {
		double T = EN - RM;
		double EPS = T / T0;
		double EPSI = 1. / EPS;
		double BHABFM = CBHAB
				* (EPSI * (EPSI * BETASI - B1) + B2 + EPS * (EPS * B4 - B3));// "EQ.2.11.1"
		return BHABFM;
	}

	/**
	 * EVALUATES EQUATION 2.7.108 OF SLAC-265
	 * @param EA EA
	 * @param K K
	 * @return the result
	 */
	public static double BREMDR(double EA, double K)// "EVALUATES EQUATION 2.7.108 OF SLAC-265"
	{
		E = EA;
		int LS = 0;
		if (E >= 50.) {
			LD = 2;
			LS = 3;
		} else {
			LD = 1;
			LS = 0;
		}
		LA = LS + 1;
		LB = LS + 2;
		double BREMDR = BREMFR(K);

		return BREMDR;
	}

	/**
	 * Called by BREMDR.
	 * @param K K
	 * @return the result
	 */
	public static double BREMFR(double K) {
		double EPS = K / E;
		double DEL = EPS / (E * (1 - EPS));
		double BREMFR = 0.;
		double A = 0.;
		double B = 0.;
		if (DEL > DELPOS[LD - 1]) {
			BREMFR = 0.0;
			return BREMFR;
		}
		double DELTA = DELCM * DEL;
		if (DELTA <= 1.) {
			A = DL1[LA - 1] + DELTA * (DL2[LA - 1] + DELTA * DL3[LA - 1]);
			B = DL1[LB - 1] + DELTA * (DL2[LB - 1] + DELTA * DL3[LB - 1]);
		} else {
			A = DL4[LA - 1] + DL5[LA - 1] * Math.log(DELTA + DL6[LA - 1]);
			B = DL4[LB - 1] + DL5[LB - 1] * Math.log(DELTA + DL6[LB - 1]);
		}
		BREMFR = (ALPHI[LD - 1] * ((1. - EPS) / (EPS * AL2)) * A + 0.5
				* (2. * EPS) * B)
				/ E;// divided by E-see eq 2.7.108
		// PEGS may be used to plot these for comparison with the more exact
		// BREMDZ
		return BREMFR;
	}

	/**
	 * Not used.
	 * @param E0 E0
	 * @return the result
	 */
	public static double BREMTR(double E0) {
		double BREMTR = 0.;
		if (E0 <= AP + RM) {
			BREMTR = 0.;
		} else {
			BREMTR = BREMRR(E0, AP, E0 - RM);
		}
		return BREMTR;
	}

	/**
	 * Called by BREMTR
	 * @param E E
	 * @param K1 K1
	 * @param K2 K2
	 * @return the result
	 */
	public static double BREMRR(double E, double K1, double K2) {
		@SuppressWarnings("unused")
		double DUMMY = BREMDR(E, K1);
		funcname = "BREMFR";
		double BREMRR = QD(funcname, K1, K2, "BREMFR");
		return BREMRR;
	}

	// "***DIFFERENTIAL CROSS SECTION FOR INCIDENT PHOTON OF ENERGY K0A TO"
	// "   COMPTON SCATTER TO SECONDARY ENERGY K."
	/**
	 * DIFFERENTIAL CROSS SECTION FOR INCIDENT PHOTON OF ENERGY K0A TO COMPTON SCATTER TO SECONDARY ENERGY K.
	 * @param K0A K0A
	 * @param K K
	 * @return the result
	 */
	public static double COMPDM(double K0A, double K) {
		K0 = K0A;
		double K0P = K0 / RM;
		CCOMP = RLC * EDEN * PI * R0 * R0 / (K0 * K0P); // "CONSTANT FACTOR IN EQ. 2.9.1"
		C1 = 1. / (K0P * K0P);
		C2 = 1. - (2. + 2. * K0P) / (K0P * K0P);
		C3 = (1. + 2. * K0P) / (K0P * K0P);
		double COMPDM = COMPFM(K);
		return COMPDM;
	}

	/**
	 * Called by COMPDM
	 * @param K K
	 * @return the result
	 */
	public static double COMPFM(double K) {
		double EPS = K / K0;
		double EPSI = 1. / EPS;
		double COMPFM = CCOMP * ((C1 * EPSI + C2) * EPSI + C3 + EPS); // "EQ.2.9.1"
		return COMPFM;
	}

	/**
	 * Not used.
	 * @param K0 K0
	 * @return the result
	 */
	public static double PAIRTR(double K0) {
		double PAIRTR = 0.;
		if (K0 <= 2. * RM) {
			PAIRTR = 0.0;
		} else {
			PAIRTR = PAIRRR(K0, RM, K0 - RM);
		}
		return PAIRTR;
	}

	/**
	 * Called by PAIRTR
	 * @param K K
	 * @param E1 E1
	 * @param E2 E2
	 * @return the result
	 */
	public static double PAIRRR(double K, double E1, double E2) {
		@SuppressWarnings("unused")
		double DUMMY = PAIRDR(K, E1);
		funcname = "PAIRFR";
		double PAIRRR = QD(funcname, E1, E2, "PAIRFR");
		return PAIRRR;
	}

	/**
	 * EVALUATES EQUATION 2.7.109 OF SLAC-265
	 * @param KA KA
	 * @param E E
	 * @return the result
	 */
	public static double PAIRDR(double KA, double E)// "EVALUATES EQUATION 2.7.109 OF SLAC-265"
	{
		K = KA;
		int LS = 0;
		if (K < 50.) {
			LE = 1;
			LS = 0;
		} else {
			LE = 2;
			LS = 3;
		}
		LA = LS + 1;
		LC = LS + 3;
		double PAIRDR = PAIRFR(E);
		return PAIRDR;
	}

	/**
	 * Called by PAIRDR
	 * @param E E
	 * @return the result
	 */
	public static double PAIRFR(double E) {
		double PAIRFR = 0.;
		double A = 0.;
		double CC = 0.;
		double DELTA = 0.;
		double EPS = E / K;
		double DEL = 1. / (K * EPS * (1. - EPS));
		if (DEL > DELPOS[LE - 1]) {
			PAIRFR = 0.0;
		} else {
			DELTA = DELCM * DEL;
			if (DELTA <= 1.) {
				A = DL1[LA - 1] + DELTA * (DL2[LA - 1] + DELTA * DL3[LA - 1]);
				CC = DL1[LC - 1] + DELTA * (DL2[LC - 1] + DELTA * DL3[LC - 1]);
			} else {
				A = DL4[LA - 1] + DL5[LA - 1] * Math.log(DELTA + DL6[LA - 1]);
				CC = DL4[LC - 1] + DL5[LC - 1] * Math.log(DELTA + DL6[LC - 1]);
			}
			PAIRFR = (ALFP1[LE - 1] * CC + ALFP2[LE - 1] * 12. * (E / K - 0.5)
					* (E / K - 0.5) * A)
					/ K;
		}
		return PAIRFR;
	}

	/**
	 * Not used.
	 * @param E E
	 * @return the result
	 */
	public static double TMXDE2(double E) {
		double ESQ = E * E;
		double BETASQ = 1.0 - RMSQ / ESQ;
		double TMXDE2 = TMXB(E) / (ESQ * BETASQ * BETASQ);
		return TMXDE2;
	}

	/**
	 * RAYLEIGH CORRECTION FACTOR (I.E., COHERENT RATIO)
	 * @param E E
	 * @return the result
	 */
	public static double CRATIO(double E) {
		// "***RAYLEIGH CORRECTION FACTOR (I.E., COHERENT RATIO)"
		double TOT = PAIRTU(E) + COMPTM(E) + PHOTTE(E);
		double CRATIO = TOT / (TOT + COHETM(E));
		return CRATIO;
	}

	/**
	 * BRANCHING RATIO BREMS/(BREMS+MOLLER)
	 * @param E E
	 * @return the result
	 */
	public static double EBR1(double E) {
		// "***BRANCHING RATIO BREMS/(BREMS+MOLLER)
		double EBR1 = 0.;
		double BREM = BREMTM(E);
		double TEBR = BREM + AMOLTM(E);
		if (TEBR > 0.0) {
			EBR1 = BREM / TEBR;
		} else {
			EBR1 = 0.0;
		}
		return EBR1;
	}

	/**
	 * ELECTRON ENERGY LOST
	 * @param E E
	 * @return the result
	 */
	public static double EDEDX(double E) {
		// "***ELECTRON ENERGY LOST IN EDIST(2)
		double EDEDX = SPTOTE(E, AE, AP);
		return EDEDX;
	}

	/**
	 * ELECTRON(-) CROSS SECTION(1/RL). 1/radiation length = RLC
	 * @param E E
	 * @return the result
	 */
	public static double ESIG(double E) {
		// "ELECTRON(-) CROSS SECTION(1/RL)"----------------1/radiation length
		// RLC
		double ESIG = BREMTM(E) + AMOLTM(E);
		return ESIG;
	}

	/**
	 * BRANCHING RATIO PAIR/(PAIR+COMPTON+PHOTO)
	 * @param E E
	 * @return the result
	 */
	public static double GBR1(double E) {
		// "***BRANCHING RATIO PAIR/(PAIR+COMPTON+PHOTO)
		double PAIR = PAIRTU(E);
		double GBR1 = PAIR / (PAIR + COMPTM(E) + PHOTTE(E));
		return GBR1;
	}

	/**
	 * BRANCHING RATIO (PAIR+COMPTON)/(PAIR+COMPTON+PHOTO)
	 * @param E E
	 * @return the result
	 */
	public static double GBR2(double E) {
		// "***BRANCHING RATIO (PAIR+COMPTON)/(PAIR+COMPTON+PHOTO)              "
		double PRCO = PAIRTU(E) + COMPTM(E);
		double GBR2 = PRCO / (PRCO + PHOTTE(E));
		return GBR2;
	}

	/**
	 * GAMMA MEAN FREE PATH
	 * @param E E
	 * @return the result
	 */
	public static double GMFP(double E) {
		// "***GAMMA MEAN FREE PATH
		double GMFP = 1.0 / (PAIRTU(E) + COMPTM(E) + PHOTTE(E));
		return GMFP;
	}

	/**
	 * BRANCHING RATIO BREMS/(BREMS+BHABHA+ANNIHILATION)
	 * @param E E
	 * @return the result
	 */
	public static double PBR1(double E) {
		// "***BRANCHING RATIO BREMS/(BREMS+BHABHA+ANNIHILATION)
		double BREM = BREMTM(E);
		double PBR1 = BREM / (BREM + BHABTM(E) + ANIHTM(E));
		return PBR1;
	}

	/**
	 * BRANCHING RATIO (BREMS+BHABHA)/(BREMS+BHABHA+ANNIHILATION)
	 * @param E E
	 * @return the result
	 */
	public static double PBR2(double E) {
		// "***BRANCHING RATIO (BREMS+BHABHA)/(BREMS+BHABHA+ANNIHILATION)
		double BRBH = BREMTM(E) + BHABTM(E);
		double PBR2 = BRBH / (BRBH + ANIHTM(E));
		return PBR2;
	}

	/**
	 * POSITRON ENERGY LOSS
	 * @param E E
	 * @return the result
	 */
	public static double PDEDX(double E) {
		// "***POSITRON ENERGY LOSS IN EDIST(2)
		double PDEDX = SPTOTP(E, AE, AP);
		return PDEDX;
	}

	/**
	 * POSITRON CROSS SECTION(1/RL)
	 * @param E E
	 * @return the result
	 */
	public static double PSIG(double E) {
		// "POSITRON CROSS SECTION(1/RL)"
		double PSIG = BREMTM(E) + BHABTM(E) + ANIHTM(E);
		return PSIG;
	}
}
