package danfulea.phys.egsOutput;

import java.io.FileWriter;
import java.util.Calendar;
import java.util.Date;
import java.util.ResourceBundle;

import javax.swing.JTextArea;

import danfulea.phys.egs.EGS4;
import danfulea.phys.egs.EGS4Core;
import danfulea.phys.egs.EGS4Geom;
import danfulea.phys.egs.EGS4Grid;
import danfulea.phys.egs.EGS4Macro;
import danfulea.phys.egs.EGS4SrcEns;
import danfulea.phys.egs.EgsQuestion;

/**
 * Cavity application<br>
 * SIMULATES THE PASSAGE OF AN ELECTRON OR PHOTON BEAM IN A FINITE, RIGHT CYLINDRICAL GEOMETRY. 
 * IT IS INTENDED FOR USE IN CALCULATING QUANTITIES OF INTEREST FOR THICK-WALLED ION CHAMBERS 
 * EXPOSED TO PHOTON BEAMS. IT ALSO MAY BE USED SIMPLY TO SCORE DOSE IN A CYLINDRICAL GEOMETRY.<br>
 * Start as minor modification of GammaDet regarding variables used in handling geometry. 
 * End up all future improvements were made on this class. Bottom line: use this class instead of GammaDet. 
 * @author Dan Fulea, 05 DEC. 2005
 *
 */
public class GammaDetA implements EgsQuestion {
	public static JTextArea jta;
	public static boolean systemOutB = true;

	private static int nearestposition = 0;
	private static final String BASE_RESOURCE_CLASS = "danfulea.phys.egsOutput.GamaDetSimResources";
	private static ResourceBundle resources = ResourceBundle
			.getBundle(BASE_RESOURCE_CLASS);
	public static boolean sourceatt = true;
	public static int sourcecode = 2;// 1-H2O;2-NaCl;3=soil;4=ciment
	public static double[][] sourcecoeftable;
	public static double[] sourcedensity;
	public static double[] energy;
	public static double[] satt_coeff;
	// the source total linear attenuation coefficient
	public static double smiu = 0.0;
	public static double smius = 0.0;
	public static double source_parcurs = 0.0;
	private static boolean indet = true;
	public static double asource = 0.0;
	public static double hsource = 0.0;
	public static double hdettot = 0.0;
	public static double adettot = 0.0;
	public static double ethick = 0.0;// source envelope thickness
	private static boolean cylSupB = false;
	public static double hsourceup = 0.0;
	public static double bsource = 0.0;

	public static int SOURCE = -1;
	public static final int SOURCE_POINT = 0;// means point source on symmetry
												// axis
	public static final int SOURCE_SARPAGAN = 1;// means volumic
												// SARPAGAN(cylinder) source on
												// top of the detecor
	public static final int SOURCE_MARINELLI = 2;// means volumic MARINELLI
													// source surrounding the
													// detecor
	public static final int SOURCE_BEAM = 3;
	public static double dist_source_det = 0.0;
	public static double beam_angle = 0.0;// degree
	public static double beam_radius = 0.0;

	public static int $NSWTCH = 8;// "# OF NRC SWITCHES FOR CONTROLLING SCATTERING "
	public static int $NBATCH = 10;// "OUTPUT BATCHES                             "
	public static int JCASE = 0;; // "no. of histories per batch"
	public static int $NCASEMIN = 100;// "min. no. of histories                        "
	public static int $MXDATA = 1040;// "MAXIMUM DATA POINTS FOR ANALYSIS (i.e.($MXREG-1))"
	public static int $MAXIT = 4;// "MAX # OF PARAMETERS TO BE SCORED             "
	// "                                (1) PRIMARY EDEP TO GAS                       "
	// "                                (2) SECONDARY EDEP TO GAS                     "
	// "                                (3) UNATTENUATED PRIMARY EDEP TO GAS          "
	// "                                     - PRIMARY EDEP TO GAS                    "
	// "                                (4) (R**2-R0**2)(UNATTENUATED PRIMARY         "
	// "                                     EDEP TO GAS)                             "
	// "                                (5) UNSOURCED UNATTENUATED PRIMARY EDEP TO    "
	// "                                     GAS WITH GAS MATERIAL REPLACED BY WALL   "
	// "                                     MATERIAL                                 "
	// "                                (6) UNSOURCED UNATTENUATED PRIMARY EDEP TO    "
	// "                                     WALL MATERIAL WITH GAS MATERIAL REPLACED "
	// "                                     BY WALL MATERIAL                         "
	public static int $MAXCMPTS = $MAXIT;// "FOR THE GRID OUTPUTS"
	public static int $MAX_SC_PLANES = 1;// "required to use phase space macros"

	public static int NCASE = 0;
	public static int NCASEO = 0;
	public static int NCASET = 0;
	public static double[][] AMASS;// ($MAXZREG,$MAXRADII),
	public static double TMCPUO = 0.0;
	public static double TIMMAX = 0.0;
	public static double STATLM = 0.0;// EIN,
	public static int[] MEDSAV;// ($MXREG),
	public static double[] RHOSAV;// ($MXREG),
	public static int IDAT = 0;
	public static int IDOPES = 0;
	public static int IRESTART = 0;
	public static int NNREADO = 0;
	public static int datcount = 0;
	// "AMASS(IZ,IX) MASS OF ZONE WITH COORDINATES (IZ,IX)
	// "TMCPUO CPU TIME USED IN PREVIOUS SESSIONS
	// "TIMMAX MAXIMUM ALLOWED CPU HOURS FOR A GIVEN CALCULATION
	// "STATLM TARGET STATISTICS IN CAVITY USED FOR AN EARLY EXIT
	// "EIN KINETIC ENERGY OF THE EXTERNAL BEAM
	// "ISUMCV(NREG) THE ARRAY OF ZONES COMPRISING THE CAVITY REGION
	// "MEDSAV(NREG) SAVES MEDIUM NUMBERS FOR CORRELATION SCORING
	// "RHOSAV(NREG) SAVES DENSITIES FOR CORRELATION SCORING
	// "IDAT = 0 STORE DATA ARRAYS FOR RE-USE
	// " = 1 DON'T STORE THEM
	// "IDOPES = 1 INCLUDES PHOTOELECTRON ANGLE SELECTION
	// " = 0 DOES NOT INCLUDE PHOTOELECTRON ANGLE SELECTION
	// "NCASE NUMBER OF HISTORIES REMAINING TO BE DONE
	// "NCASEO NUMBER OF HISTORIES DONE IN PREVIOUS SESSIONS
	// "NCASET NUMBER OF HISTORIES ALREADY DONE
	// "NNREADO TOTAL NO. OF PARTICLES READ FROM PHSP SOURCE IN PREVIOUS
	// " RUNS
	// "IRESTART = 0 => INITIAL RUN
	// " = 1 => RESTARTED RUN
	// " = 3 => DATA ANALYSIS ONLY
	// " = 4 => READ STARTING RANDOM NUMBERS FROM A FILE
	// " = 5 => analyse previous parallel runs

	public static double RRZ = 0.0;//
	public static double RRCUT = 0.0;
	public static boolean RUSROU = false;
	// "RRZ COORDINATE OF PLANE AT WHICH RUSSIAN ROULETTE IS PLAYED
	// "RRCUT SURVIVAL PROBABILITY AFTER CROSSING THE PLANE
	// "RUSROU = .FALSE. => RUSSIAN ROULETTE WILL NOT BE PLAYED
	// " = .TRUE. => RUSSIAN ROULETTE WILL BE PLAYED
	public static double cav_dose = 0.0;// => total dose in cavity
	public static double cav_dose0 = 0.0;// => primary dose in cavity
	public static double cav_dose1 = 0.0;// => primary dose corrected for
											// attenuation in cavity
	public static double cav_dose2 = 0.0;// => secondary dose in cavity
	public static double cav2_dose = 0.0;// => total dose squared in cavity
	public static double cav2_dose0 = 0.0;// => primary dose squared in cavity
	public static double cav2_dose1 = 0.0;// => unattenuated primary dose
											// squared in cavity
	public static double cav2_dose2 = 0.0;// => secondary dose squared in cavity
	// " => eventually, they hold the uncertainty in their respective doses
	public static double cav_dosec = 0.0;// => correlation total -
											// primary(unatttenuated) in cavity
	public static double cav_dosec01 = 0.0;// => correlation primary - primary
											// (unattenuated) in cavity
	public static double cav_dosec02 = 0.0;// => correlation secondary - primary
											// in cavity
	public static double SCSTP = 0.0;// => total no. of charged particle steps
	public static double SCSTP2 = 0.0;// => total no. of charged particle steps
										// squared--eventually holds
	// " uncertainty in SCSTP
	public static double SCCSTP = 0.0;// => no. of charged particle steps in
										// cavity
	public static double SCCSTP2 = 0.0;// => no. of charged particle steps in
										// cavity squared--eventually
	// " holds uncertainty in SCCSTP
	public static double[][][] SCDOSE;// (IZ,IX,IT) => dose in cavity voxel
										// IZ,IX:
	// " IT=1 -- total
	// " IT=2 -- primary
	// " IT=3 -- primary corrected for attenuation
	// " IT=4 -- secondary
	public static double[][][] SCDOSE2;// (IZ,IX,IT) => dose in cavity voxel
										// IZ,IX squared. IT same as above.
	// " Eventually holds uncertainties in respective doses.
	public static double[][][] SCDOSE_COV;// (IZ,IX,IT) => correlation in cavity
											// voxel IZ,IX for:
	// " IT=1 -- total and primary (unattenuated)
	// " IT=2 -- primary and primary (unattenuated)
	// " IT=3 -- secondary and primary

	public static double PIISTP = 0.0;// => no. of PRESTA-II steps from previous
										// runs
	// "last_case => last primary history to score dose in cavity
	public static int SCSTP_LAST = 0;// => last primary history to score charged
										// particle step
	public static int SCCSTP_LAST = 0;// => last primary history to score
										// charged particle step
	// " in cavity
	public static int[][] SCDOSE_LAST;// (IZ,IX) => last primary history to
										// score dose in cavity
	// " voxel IZ,IX.
	public static double tmp_dose = 0.0;// => temporay arrays for scoring the
										// different
	public static double tmp_dose0 = 0.0;// dose components in the cavity
	public static double tmp_dose1 = 0.0;//
	public static double tmp_dose2 = 0.0;//
	public static double corr_02 = 0.0;// => correlation primary - secondary
										// cavity doses
	public static double SCSTP_TMP = 0.0;// => temp. variable for scoring total
											// no. of charged
	// " particle steps
	public static double SCCSTP_TMP = 0.0;// => temp. variable for scoring total
											// no. of charged
	// " particle steps in cavity
	public static double[][][] SCDOSE_TMP;// (IZ,IX,IT) => temp. variable for
											// scoring dose in cavity voxel
	// " IZ,IX. IT same as for SCDOSE.
	// "cs_enhance => cross-section enhancement factor

	public static int MXNP = 0;// MAXIMUM LEVEL TO WHICH THE STACK OF DAUGHTER
								// PARTICLES FROM AN
	// " INCIDENT PARTICLE RISES (STACK MAY INCLUDE INCIDENT PARTICLE)
	// "IFULL = 0 JUST CALCULATE TOTAL DOSE AND THAT DUE TO STOPPERS
	// " AND DISCARDS (THE DEFAULT)
	// " = 1 ABOVE PLUS Aatt, Ascat
	public static int ISTORE = 0;// = 0 DO NOT STORE THE INITIAL RANDOM NUMBERS
									// (THE DEFAULT)
	// " = 1 STORE THE INITIAL RANDOM NUMBER FOR THE LAST HISTORY
	// " = 2 STORE THE INITIAL RANDOM NUMBER FOR ALL HISTORIES
	// " THAT DEPOSIT ENERGY IN THE CAVITY
	// " = 3 STORE ALL THE INITIAL RANDOM NUMBERS
	public static int IWATCH = 0;// = 0 FOR NORMAL OUTPUT (THE DEFAULT)
	// " = 1 OUTPUT ON EVERY DISCRETE INTERACTION
	// " = 2 OUTPUT ON EVERY ELECTRON/PHOTON STEP AS WELL
	// " = 3 PRINTS OUT ONLY WHEN ENERGY IS DEPOSITED
	// " = 4 PRINTS OUT FILE FOR GRAPHICS
	public static int IOOPTN = 0;// = 0 SHORT OUTPUT (THE DEFAULT) -JUST CAVITY
									// SUMMARY
	// " AND THE MATERIAL GRID
	// " = 1 ABOVE PLUS OUTPUT GRID
	public static int IOUTSP = 0;// = 0 NO SPECTRUM DATA IN OUTPUT SUMMARY
	// " = 1 INCLUDE SPECTRUM DATA IN OUTPUT SUMMARY
	public static int IFANO = 0;// = 0 NO PHOTON REGENERATION
	// " = 1 PHOTONS REGENERATED AFTER THEY HAVE INTERACTED
	// " = 2 NO PHOTON REGENERATION, ELECTRONS FROM CAVITY WALL ARE ELIMINATED

	public static final String NULLs = "";
	// -------------INPUTS-------------------------------------------
	public static String TITLEs = "";
	// @ IWATCH
	public static final int IWATCH_OFF = 0;
	public static final int IWATCH_INTERACTIONS = 1;
	public static final int IWATCH_STEPS = 2;
	public static final int IWATCH_DEPOSITED = 3;
	public static final int IWATCH_GRAPH = 4;
	// @ STORE INITIAL RANDOM NUMBERS
	public static final int ISTORE_NO = 0;
	public static final int ISTORE_LAST = 1;
	public static final int ISTORE_ALL_DEPOSITED = 2;
	public static final int ISTORE_ALL = 3;
	// @ IRESTART
	public static final int IRESTART_FIRST = 0;
	public static final int IRESTART_RESTART = 1;
	public static final int IRESTART_ANALYZE = 2;
	public static final int IRESTART_START_RNS = 3;
	public static final int IRESTART_PARALLEL = 4;
	// @ OUTPUT OPTIONS
	public static final int IOOPTN_SHORT = 0;
	public static final int IOOPTN_DETAIL = 1;
	// @ STORE DATA ARRAYS
	public static final int IDAT_YES = 0;
	public static final int IDAT_NO = 1;
	// @ NUMBER OF HISTORIES
	public static final int NCASE_MIN = 1;
	public static final int NCASE_MAX = 461168600;// 4.611686e18;
	public static final int NCASE_DEFAULT = 20000;
	// @ MAX CPU HOURS ALLOWED
	public static final double TIMMAX_MIN = 0.0;
	public static final double TIMMAX_MAX = 1000.0;
	public static final double TIMMAX_DEFAULT = 999.0;
	// @ IFULL
	public static final int IFULL_DOSE_AND_STOPPERS = 0;
	public static final int IFULL_AATT_AND_ASCAT = 1;
	// @ STATISTICAL ACCURACY SOUGHT
	public static final double STATLM_MIN = 0.0;
	public static final double STATLM_MAX = 100.0;
	public static final double STATLM_DEFAULT = 0.0;
	// @ PHOTON REGENERATION
	public static final int IFANO_NO = 0;
	public static final int IFANO_YES = 1;
	public static final int IFANO_NO_ELECTRONS_FROM_WALL = 2;
	// @ INITIAL RANDOM NO. SEEDS
	public static final int RANLUX_LEVEL_MIN = 0;
	public static final int RANLUX_LEVEL_MAX = 4;
	public static final int RANLUX_LEVEL_DEFAULT = EGS4.$DEFAULT_LL;// see egs4
	public static final int RANLUX_SEED_MIN = 1;
	public static final int RANLUX_SEED_MAX = 1073741824;// default seed is set
															// in
															// egs4.setdefaults!!
	public static final int RANLUX_SEED_DEFAULT = 999999;
	public static final int RANMAR_SEED_MIN = 1;
	public static final int RANMAR_SEED_MAX = 30081;
	public static final int RANMAR_SEED_DEFAULT = 9373;
	public static int jrng1 = 0;
	public static int jrng2 = 0;
	// @TRANSPORT PARAMS:
	public static double ecut = 0.0;
	public static double pcut = 0.0;
	public static double smax = 0.0;
	public static boolean setEcutRegion = false;
	public static boolean setPcutRegion = false;
	public static boolean setSmaxRegion = false;
	public static final double ECUT_MIN = 0.0;// " ECUT "
	public static final double ECUT_MAX = 1.E15;
	public static final double PCUT_MIN = 0.0;// " PCUT "
	public static final double PCUT_MAX = 1.E15;
	public static final double SMAXIR_MIN = 0.0;// " SMAX "
	public static final double SMAXIR_MAX = 1.E15;
	public static int nEcut = 1;// number of data
	public static double[] Ecut;// =new double[EGS4.$MXREG];
	public static int[] startEcutRegion;
	public static int[] stopEcutRegion;
	public static int nPcut = 1;// number of data
	public static double[] Pcut;// =new double[EGS4.$MXREG];
	public static int[] startPcutRegion;
	public static int[] stopPcutRegion;
	public static int nSmax = 1;// number of data
	public static double[] Smax;// =new double[EGS4.$MXREG];
	public static int[] startSmaxRegion;
	public static int[] stopSmaxRegion;

	public static boolean setIncohRegion = false;
	public static boolean setCohRegion = false;
	public static boolean setRelaxRegion = false;
	public static boolean setPeRegion = false;
	public static int nIncoh = 1;// number of data
	public static int[] Incoh;// =new int[EGS4.$MXREG];
	public static int[] startIncohRegion;
	public static int[] stopIncohRegion;
	public static int nCoh = 1;// number of data
	public static int[] Coh;// =new int[EGS4.$MXREG];
	public static int[] startCohRegion;
	public static int[] stopCohRegion;
	public static int nRelax = 1;// number of data
	public static int[] Relax;// =new int[EGS4.$MXREG];
	public static int[] startRelaxRegion;
	public static int[] stopRelaxRegion;
	public static int nPe = 1;// number of data
	public static int[] Pe;// =new int[EGS4.$MXREG];
	public static int[] startPeRegion;
	public static int[] stopPeRegion;

	// " Incoherent (Compton) scattering "
	public static int incoh = 0;
	public static final int incoh_OFF = 0;
	public static final int incoh_ON = 1;
	public static final int incoh_ON_IN_REGIONS = 2;
	public static final int incoh_OFF_IN_REGIONS = 3;
	// " Radiative corrections for Compton scattering "
	public static final int radc_OFF = 0;
	public static final int radc_ON = 1;
	// " Coherent (Rayleigh) scattering "
	public static int coh = 0;
	public static final int coh_OFF = 0;
	public static final int coh_ON = 1;
	public static final int coh_ON_IN_REGIONS = 2;
	public static final int coh_OFF_IN_REGIONS = 3;
	// " Atomic Relaxations "
	public static int relax = 0;
	public static final int relax_OFF = 0;
	public static final int relax_ON = 1;
	public static final int relax_ON_IN_REGIONS = 2;
	public static final int relax_OFF_IN_REGIONS = 3;
	// " Photoelectron angular sampling "
	public static int pe = 0;
	public static final int pe_ang_OFF = 0;
	public static final int pe_ang_ON = 1;
	public static final int pe_ang_ON_IN_REGIONS = 2;
	public static final int pe_ang_OFF_IN_REGIONS = 3;
	// " Bremsstrahlung angular sampling "
	public static final int brems_ang_SIMPLE = 0;
	public static final int brems_ang_KM = 1;
	// " Bremsstrahlung cross sections "
	public static final int brems_cross_BH = 0;
	public static final int brems_cross_NIST = 1;
	// " Pair angular sampling "
	public static final int pair_ang_OFF = 0;
	public static final int pair_ang_SIMPLE = 1;
	public static final int pair_ang_KM = 2;
	public static final int pair_ang_UNIFORM = 3;
	public static final int pair_ang_BLEND = 4;
	// " Pair cross sections "
	public static final int pair_cross_BH = 0;
	public static final int pair_cross_NRC = 1;
	// " Triplet production "
	public static final int triplet_OFF = 0;
	public static final int triplet_ON = 1;
	// " Spin effects  		"
	public static int ispin = 0;
	public static final int spin_OFF = 0;
	public static final int spin_ON = 1;
	// " Electron impact ionization "
	public static final int eii_OFF = 0;
	public static final int eii_ON = 1;
	public static final int eii_casnati = 2;
	public static final int eii_kolbenstvedt = 3;
	public static final int eii_gryzinski = 4;

	public static String photon_xsection = "";// default

	public static final double ESTEPE_MIN = 1.0E-5;// " ESTEPE "
	public static final double ESTEPE_MAX = 1.0;
	public static final double XIMAX_MIN = 0.0;// " XIMAX "
	public static final double XIMAX_MAX = 1.0;
	// " BCA "
	public static final int BCA_EXACT = 0;
	public static final int BCA_PRESTA_I = 1;

	public static final double Skindepth_MIN = -1.0;// " Skindepth "
	public static final double Skindepth_MAX = 1.0e15;
	// " Electron-step algorithm "
	public static final int estep_alg_PRESTA_I = 1;
	public static final int estep_alg_PRESTA_II = 0;

	// Variance reduction
	public static final int irejct_OFF = 0;
	public static final int irejct_ON = 1;
	public static double ESAVEIN = 0.0;
	public static final double ESAVEIN_MIN = 0.0;
	public static final double ESAVEIN_MAX = 1.0e30;
	public static final double ESAVEIN_DEFAULT = 0.0;
	public static final double cs_enhance_MIN = 0.0;
	public static final double cs_enhance_MAX = 1.0e6;
	public static final double cs_enhance_DEFAULT = 0.5;
	public static final double RRDEPTH_MIN = -1.0e30;
	public static final double RRDEPTH_MAX = 1.0e30;
	public static final double RRDEPTH_DEFAULT = 0.0;
	public static final double RRFRACTION_MIN = -1.0e30;
	public static final double RRFRACTION_MAX = 1.0e30;
	public static final double RRFRACTION_DEFAULT = 0.0;
	public static final double EXPC_MIN = -1.0e30;
	public static final double EXPC_MAX = 1.0e30;
	public static final double EXPC_DEFAULT = 0.0;
	public static int IFARCE = 0;
	public static final int IFARCE_ON = 1;
	public static final int IFARCE_OFF = 0;
	public static final int NFMIN_MIN = 0;
	public static final int NFMIN_DEFAULT = 1;
	public static int NFMIN_MAX = 0;
	public static final int NFMAX_MIN = 0;
	public static final int NFMAX_DEFAULT = 1;
	public static int NFMAX_MAX = 0;
	public static int phsplitt = 1;
	public static final int phsplitt_MIN = 1;
	public static final int phsplitt_DEFAULT = 1;
	public static int phsplitt_MAX = 0;

	// ###########################################
	public static int IHSTRY = 0;// => counter for total no. of histories
	public static double EI = 0.0;
	public static double EKMAX = 0.0;
	public static double DEPTH = 0.0;
	public static double VOLUME = 0.0;
	public static double RLOW2 = 0.0;
	public static int ISUMX = 0;

	public static double FMASSC = 0.0;
	public static double FMASS = 0.0;
	public static double AINFLU_CURRENT = 0.0; // -------
	public static int last_case = 0;
	public static int n_photons = 0;
	public static int n_electrons = 0;
	public static double sumE_photons = 0.0;
	public static double sumE2_photons = 0.0;
	public static double sumE_electrons = 0.0;
	public static double sumE2_electrons = 0.0;
	public static int IBTCH = 0;
	public static double dcav_current = 0.0;
	public static double dcavun = 0.0;
	public static int IDECAV = 0;
	public static int ipass = 0;// not usded
	public static double SCORE_NORM_NUM = 0.0;
	public static double SCORE_TEMP = 0.0;

	public static int IRL = 0;
	public static int MEDNUM = 0;

	// ========================
	public static int IPRINT = 2;// default summary;
	public static boolean createOutputFile = false;
	public static boolean putInFile = false;// internal var defining when and
											// what to print
	private String filename = "";
	FileWriter sigfos;

	public static int ieff = 0;
	public static double edep = 0.0;
	public static double edep2 = 0.0;
	public static double edeptmp = 0.0;
	public static int[] iieff;// =0;
	public static double dedep = 0.0;
	public static double edep_current = 0.0;
	public static double pond = 0.0;
	public static double pondt = 0.0;
	public static double[] pondt_batch;

	public static double ETHRESHOLD = 0.0;
	private boolean ekinokB = true;
	// -------------end INPUTS-------------------------------------------
	// REPLACE {$INVDIM} WITH {1000} "DIMENSION CONTROLS GRID SIZE FOR INVERSE"
	// CDFINV($INVDIM,2)

	public static boolean is_finished = false;

	/**
	 * Constructor
	 */
	public GammaDetA() {
		// createOutputFile=false;
		// createOutputFile=true;
		putInFile = true;
		Calendar cal = Calendar.getInstance();
		String fs = cal.get(Calendar.YEAR) + "_" + cal.get(Calendar.MONTH) + "_"
				+ cal.get(Calendar.DAY_OF_MONTH) + "_" + cal.get(Calendar.HOUR) + "_"
				+ cal.get(Calendar.MINUTE) + "_" + cal.get(Calendar.SECOND) + ".txt";
		filename = fs;// will be 2005_10_25_14_30_56.txt
		// =========================
		sourcecoeftable = (double[][]) resources
				.getObject("source.attenuationCoef");
		sourcedensity = (double[]) resources.getObject("source.densities");// scod-1
																			// will
																			// go
		// the source total linear attenuation coefficient
		smiu = 0.;
		smius = sourcedensity[sourcecode - 1];// *mass total
												// ....@@@@@@@@@@@@@@@@@@@@@
		satt_coeff = new double[sourcecoeftable.length];
		energy = new double[sourcecoeftable.length];
		for (int i = 0; i < satt_coeff.length; i++) {
			energy[i] = sourcecoeftable[i][0];
			satt_coeff[i] = sourcecoeftable[i][sourcecode];
		}
		// =========================
		ekinokB = true;// ETHRESHOLD affect in "hidden mode" only spectrum!!but
						// here no eff calculation is done!!
		// //mono: directly ein <ethresh=>no shower!!
		init();
	}

	/**
	 * Reset global variables for re-use
	 */
	public static void reset() {
		nearestposition = 0;
		sourceatt = true;
		sourcecode = 2;
		smiu = 0.0;
		smius = 0.0;
		source_parcurs = 0.0;
		indet = true;
		asource = 0.0;
		hsource = 0.0;
		cylSupB = false;
		hsourceup = 0.0;
		bsource = 0.0;
		SOURCE = -1;
		$NSWTCH = 8;
		$NBATCH = 10;
		JCASE = 0;
		$NCASEMIN = 100;
		$MXDATA = 1040;
		$MAXIT = 4;
		$MAXCMPTS = $MAXIT;
		$MAX_SC_PLANES = 1;
		NCASE = 0;
		NCASEO = 0;
		NCASET = 0;
		TMCPUO = 0.0;
		TIMMAX = 0.0;
		STATLM = 0.0;
		IDAT = 0;
		IDOPES = 0;
		IRESTART = 0;
		NNREADO = 0;
		datcount = 0;
		RRZ = 0.0;
		RRCUT = 0.0;
		RUSROU = false;
		cav_dose = 0.0;
		cav_dose0 = 0.0;
		cav_dose1 = 0.0;
		cav_dose2 = 0.0;
		cav2_dose = 0.0;
		cav2_dose0 = 0.0;
		cav2_dose1 = 0.0;
		cav2_dose2 = 0.0;
		cav_dosec = 0.0;
		cav_dosec01 = 0.0;
		cav_dosec02 = 0.0;
		SCSTP = 0.0;
		SCSTP2 = 0.0;
		SCCSTP = 0.0;
		SCCSTP2 = 0.0;
		PIISTP = 0.0;
		SCSTP_LAST = 0;
		SCCSTP_LAST = 0;
		tmp_dose = 0.0;
		tmp_dose0 = 0.0;
		tmp_dose1 = 0.0;
		tmp_dose2 = 0.0;
		corr_02 = 0.0;
		SCSTP_TMP = 0.0;
		SCCSTP_TMP = 0.0;
		MXNP = 0;
		ISTORE = 0;
		IWATCH = 0;
		IOOPTN = 0;
		IOUTSP = 0;
		IFANO = 0;
		TITLEs = "";
		jrng1 = 0;
		jrng2 = 0;
		ecut = 0.0;
		pcut = 0.0;
		smax = 0.0;
		setEcutRegion = false;
		setPcutRegion = false;
		setSmaxRegion = false;
		nEcut = 1;
		nPcut = 1;
		nSmax = 1;
		setIncohRegion = false;
		setCohRegion = false;
		setRelaxRegion = false;
		setPeRegion = false;
		nIncoh = 1;
		nCoh = 1;
		nRelax = 1;
		nPe = 1;
		incoh = 0;
		coh = 0;
		relax = 0;
		pe = 0;
		ispin = 0;
		photon_xsection = "";
		ESAVEIN = 0.0;
		IFARCE = 0;
		NFMIN_MAX = 0;
		NFMAX_MAX = 0;
		phsplitt = 1;
		phsplitt_MAX = 0;
		IHSTRY = 0;
		EI = 0.0;
		EKMAX = 0.0;
		DEPTH = 0.0;
		VOLUME = 0.0;
		RLOW2 = 0.0;
		ISUMX = 0;
		FMASSC = 0.0;
		FMASS = 0.0;
		AINFLU_CURRENT = 0.0;
		last_case = 0;
		n_photons = 0;
		n_electrons = 0;
		sumE_photons = 0.0;
		sumE2_photons = 0.0;
		sumE_electrons = 0.0;
		sumE2_electrons = 0.0;
		IBTCH = 0;
		dcav_current = 0.0;
		dcavun = 0.0;
		IDECAV = 0;
		ipass = 0;
		SCORE_NORM_NUM = 0.0;
		SCORE_TEMP = 0.0;
		IRL = 0;
		MEDNUM = 0;
		IPRINT = 2;// default summary;
		ieff = 0;
		edep = 0.0;
		edep2 = 0.0;
		edeptmp = 0.0;
		dedep = 0.0;
		edep_current = 0.0;
		pond = 0.0;
		pondt = 0.0;
		ETHRESHOLD = 0.0;
		is_finished = false;
	}

	/**
	 * Reset global array variables for re-use.
	 */
	private void destroyArrays() {
		AMASS = new double[0][0];// ($MAXZREG,$MAXRADII),
		MEDSAV = new int[0];// ($MXREG),
		RHOSAV = new double[0];// ($MXREG),

		/*
		 * SCPDST= new double[0];SCPDST2= new double[0];SCPCUM= new double[0];
		 * SCPCUM2= new double[0];ECUM= new double[0];ECUM2= new double[0];
		 * SCDFBK=new double[0];SCDFBK2=new double[0]; SCDFEP=new
		 * double[0];SCDFEP2=new double[0]; SCDFDIFF=new double[0];SCDFDIFF2=new
		 * double[0]; DFEN=new double[0][0];//(4,4), BINTOP= new
		 * double[0];//($EBIN), CDFINV=new double[0][0]; REGSVOL=new
		 * int[0];//=0; TOPEBIN=new double[0]; NENHLO=new int[0];NENHHI=new
		 * int[0]; IPRTSP=new int[0]; EBINW=new double[0][0];
		 * 
		 * //===================from init: IPRTSP=new int[0];
		 */
		SCDOSE = new double[0][0][0];
		SCDOSE2 = new double[0][0][0];
		SCDOSE_COV = new double[0][0][0];
		SCDOSE_LAST = new int[0][0];
		SCDOSE_TMP = new double[0][0][0];
		// ($MAXZREG,$MAXRADII,$MAXIT)
		// EGS4Grid.CDSTBL=new String[0];
		// EGS4Grid.CTRTBL=new String[0];
		// EGS4Grid.CABSRB=new String[0];
		EGS4Grid.CAVTRACK = new String[0];

		Ecut = new double[0];
		startEcutRegion = new int[0];
		stopEcutRegion = new int[0];
		Pcut = new double[0];
		startPcutRegion = new int[0];
		stopPcutRegion = new int[0];
		Smax = new double[0];
		startSmaxRegion = new int[0];
		stopSmaxRegion = new int[0];
		Incoh = new int[0];
		startIncohRegion = new int[0];
		stopIncohRegion = new int[0];
		Coh = new int[0];
		startCohRegion = new int[0];
		stopCohRegion = new int[0];
		Relax = new int[0];
		startRelaxRegion = new int[0];
		stopRelaxRegion = new int[0];
		Pe = new int[0];
		startPeRegion = new int[0];
		stopPeRegion = new int[0];
		EGS4Grid.CDSTBL = new String[0];// ////////////here at init
		EGS4Grid.CTRTBL = new String[0];
		EGS4Grid.CABSRB = new String[0];
		// ===============GRID reset==================================
		EGS4Grid.RADIAL_BINS = new double[0];
		EGS4Grid.DEPTH_BINS = new double[0];
		EGS4Grid.RESULTS = new double[0][0][0];
		EGS4Grid.UNCRT = new double[0][0][0];// ($MAXZREG, $MAXRADII,
												// $MAXCMPTS),
		EGS4Grid.LABELS = new String[0];
		EGS4Grid.EXPLANATIONS = new String[0];
		// =============SrcEns reset===================================
		EGS4SrcEns.ENSRCD = new double[0];
		EGS4SrcEns.SRCPDF = new double[0];
		EGS4SrcEns.srcpdf_at = new double[0];
		EGS4SrcEns.srcbin_at = new int[0];
		EGS4SrcEns.RCDFIN = new double[0][0];
		EGS4SrcEns.RDISTF = new double[0];
		EGS4SrcEns.RPDF = new double[0];
		EGS4SrcEns.RCDF = new double[0];
		EGS4SrcEns.source_option = new double[0];
		EGS4SrcEns.OMEGIS = new double[0];
		// =============EGS4Macro reset===================================
		EGS4Macro.iefl = new int[0];
		EGS4Macro.rangerr0 = new double[0];
		EGS4Macro.rangerr1 = new double[0];
		EGS4Macro.NP_INC = new int[0];
		EGS4Macro.GWAITf = new double[0];
		EGS4Macro.GWATE = new double[0];
		// =============EGS4Geom reset===================================
		EGS4Geom.ZPLANE = new double[0];
		EGS4Geom.RCYL = new double[0];
		EGS4Geom.RSPH = new double[0];
		EGS4Geom.CYRAD2 = new double[0];
		EGS4Geom.ntrack = new int[0];
		EGS4Geom.cosalp = new double[0];
		EGS4Geom.TANALP = new double[0];
		EGS4Geom.TANAL2 = new double[0];
		EGS4Geom.SINALP = new double[0];
		EGS4Geom.ALPHA = new double[0];
		EGS4Geom.RSPH2 = new double[0];
		EGS4Geom.NZHI = new int[0];
		EGS4Geom.NZLO = new int[0];
		EGS4Geom.NRHI = new int[0];
		EGS4Geom.NRLO = new int[0];
		EGS4Geom.NSLAB = new int[0];
		EGS4Geom.DELTAZ = new double[0];
		EGS4Geom.MEDNUM = new int[0];
		EGS4Geom.RHOR = new double[0];
		EGS4Geom.NREGLO = new int[0];
		EGS4Geom.NREGHI = new int[0];
		EGS4Geom.NRADIUS = new int[0];
		EGS4Geom.NCONES = new int[0];
		EGS4Geom.ANGLES = new double[0];
		EGS4Geom.RADII = new double[0];
		EGS4Geom.cavreg = new int[0];
		// =============EGS4 reset===================================
		EGS4.ECUT = new double[0];
		EGS4.PCUT = new double[0];
		EGS4.E = new double[0];
		EGS4.X = new double[0];
		EGS4.Y = new double[0];
		EGS4.Z = new double[0];
		EGS4.U = new double[0];
		EGS4.V = new double[0];
		EGS4.W = new double[0];
		EGS4.DNEAR = new double[0];
		EGS4.WT = new double[0];
		EGS4.IQ = new int[0];
		EGS4.IR = new int[0];
		EGS4.LATCH = new int[0];
		EGS4.RHOR = new double[0];
		EGS4.MED = new int[0];
		EGS4.IRAYLR = new int[0];
		EGS4.MEDIA = new String[0];
		EGS4.RLC = new double[0];
		EGS4.RLDU = new double[0];
		EGS4.RHO = new double[0];
		EGS4.MSGE = new int[0];
		EGS4.MGE = new int[0];
		EGS4.MSEKE = new int[0];
		EGS4.MEKE = new int[0];
		EGS4.MLEKE = new int[0];
		EGS4.MCMFP = new int[0];
		EGS4.MRANGE = new int[0];
		EGS4.IRAYLM = new int[0];
		EGS4.iz_array = new int[0];
		EGS4.be_array = new double[0];
		EGS4.Jo_array = new double[0];
		EGS4.erfJo_array = new double[0];
		EGS4.ne_array = new int[0];
		EGS4.shn_array = new int[0];
		EGS4.shell_array = new int[0][0];
		EGS4.eno_array = new double[0][0];
		EGS4.n_shell = new int[0];
		EGS4.ibcmp = new int[0];
		EGS4.binding_energies = new double[0][0];
		EGS4.interaction_prob = new double[0][0];
		EGS4.relaxation_prob = new double[0][0];
		EGS4.edge_energies = new double[0][0];
		EGS4.edge_number = new int[0];
		EGS4.edge_a = new double[0][0];
		EGS4.edge_b = new double[0][0];
		EGS4.edge_c = new double[0][0];
		EGS4.edge_d = new double[0][0];
		EGS4.iedgfl = new int[0];
		EGS4.iphter = new int[0];
		EGS4.eii_xsection_a = new double[0];
		EGS4.eii_xsection_b = new double[0];
		EGS4.eii_z = new int[0];
		EGS4.eii_sh = new int[0];
		EGS4.eii_a = new double[0];
		EGS4.eii_b = new double[0];
		EGS4.eii_nshells = new int[0];
		EGS4.eii_nsh = new int[0];
		EGS4.eii_cons = new double[0];
		EGS4.eii_first = new int[0][0];
		EGS4.eii_no = new int[0][0];
		EGS4.SMAXIR = new double[0];
		EGS4.e_max_rr = new double[0];
		EGS4.i_do_rr = new int[0];
		EGS4.AP = new double[0];
		EGS4.AE = new double[0];
		EGS4.UP = new double[0];
		EGS4.UE = new double[0];
		EGS4.TE = new double[0];
		EGS4.THMOLL = new double[0];
		EGS4.DL1 = new double[0][0];
		EGS4.DL2 = new double[0][0];
		EGS4.DL3 = new double[0][0];
		EGS4.DL4 = new double[0][0];
		EGS4.DL5 = new double[0][0];
		EGS4.DL6 = new double[0][0];
		EGS4.ALPHI = new double[0][0];
		EGS4.BPAR = new double[0][0];
		EGS4.DELPOS = new double[0][0];
		EGS4.WA = new double[0][0];
		EGS4.PZ = new double[0][0];
		EGS4.ZELEM = new double[0][0];
		EGS4.RHOZ = new double[0][0];
		EGS4.PWR2I = new double[0];
		EGS4.DELCM = new double[0];
		EGS4.ZBRANG = new double[0];
		EGS4.LZBRANG = new double[0];
		EGS4.NNE = new int[0];
		EGS4.ASYM = new String[0][0];
		EGS4.nb_fdata = new double[0][0][0];
		EGS4.nb_xdata = new double[0][0][0];
		EGS4.nb_wdata = new double[0][0][0];
		EGS4.nb_idata = new int[0][0][0];
		EGS4.nb_emin = new double[0];
		EGS4.nb_emax = new double[0];
		EGS4.nb_lemin = new double[0];
		EGS4.nb_lemax = new double[0];
		EGS4.nb_dle = new double[0];
		EGS4.nb_dlei = new double[0];
		EGS4.log_ap = new double[0];
		EGS4.iausfl = new int[0];
		EGS4.ums_array = new double[0][0][0];
		EGS4.fms_array = new double[0][0][0];
		EGS4.wms_array = new double[0][0][0];
		EGS4.ims_array = new int[0][0][0];
		EGS4.spin_rej = new double[0][0][0][0][0];
		EGS4.SIN0 = new double[0];
		EGS4.SIN1 = new double[0];
		EGS4.esig_e = new double[0];
		EGS4.psig_e = new double[0];
		EGS4.range_ep = new double[0][0][0];
		EGS4.E_array = new double[0][0];
		EGS4.etae_ms0 = new double[0][0];
		EGS4.etae_ms1 = new double[0][0];
		EGS4.etap_ms0 = new double[0][0];
		EGS4.etap_ms1 = new double[0][0];
		EGS4.q1ce_ms0 = new double[0][0];
		EGS4.q1ce_ms1 = new double[0][0];
		EGS4.q1cp_ms0 = new double[0][0];
		EGS4.q1cp_ms1 = new double[0][0];
		EGS4.q2ce_ms0 = new double[0][0];
		EGS4.q2ce_ms1 = new double[0][0];
		EGS4.q2cp_ms0 = new double[0][0];
		EGS4.q2cp_ms1 = new double[0][0];
		EGS4.blcce0 = new double[0][0];
		EGS4.blcce1 = new double[0][0];
		EGS4.expeke1 = new double[0];
		EGS4.EKE0 = new double[0];
		EGS4.EKE1 = new double[0];
		EGS4.XR0 = new double[0];
		EGS4.TEFF0 = new double[0];
		EGS4.BLCC = new double[0];
		EGS4.XCC = new double[0];
		EGS4.ESIG0 = new double[0][0];
		EGS4.ESIG1 = new double[0][0];
		EGS4.PSIG0 = new double[0][0];
		EGS4.PSIG1 = new double[0][0];
		EGS4.EDEDX0 = new double[0][0];
		EGS4.EDEDX1 = new double[0][0];
		EGS4.PDEDX0 = new double[0][0];
		EGS4.PDEDX1 = new double[0][0];
		EGS4.EBR10 = new double[0][0];
		EGS4.EBR11 = new double[0][0];
		EGS4.PBR10 = new double[0][0];
		EGS4.PBR11 = new double[0][0];
		EGS4.PBR20 = new double[0][0];
		EGS4.PBR21 = new double[0][0];
		EGS4.TMXS0 = new double[0][0];
		EGS4.TMXS1 = new double[0][0];
		EGS4.IUNRST = new int[0];
		EGS4.EPSTFL = new int[0];
		EGS4.IAPRIM = new int[0];
		EGS4.EBINDA = new double[0];
		EGS4.GE0 = new double[0];
		EGS4.GE1 = new double[0];
		EGS4.GMFP0 = new double[0][0];
		EGS4.GMFP1 = new double[0][0];
		EGS4.GBR10 = new double[0][0];
		EGS4.GBR11 = new double[0][0];
		EGS4.GBR20 = new double[0][0];
		EGS4.GBR21 = new double[0][0];
		EGS4.RCO0 = new double[0];
		EGS4.RCO1 = new double[0];
		EGS4.RSCT0 = new double[0][0];
		EGS4.RSCT1 = new double[0][0];
		EGS4.COHE0 = new double[0][0];
		EGS4.COHE1 = new double[0][0];
		EGS4.MPGEM = new int[0][0];
		EGS4.NGR = new int[0];
		EGS4.rng_array = new double[0];
		EGS4.rng_array1 = new double[0];
		EGS4.seeds = new int[0];
		EGS4.state = new int[0];
		EGS4.next = new int[0];
		EGS4.urndm = new int[0];
		EGS4.DATA = new double[0][0];
		EGS4.radc_sigs = new double[0];
		EGS4.radc_sigd = new double[0];
		EGS4.radc_frej = new double[0][0];
		EGS4.radc_x = new double[0];
		EGS4.radc_fdat = new double[0];
		EGS4.radc_Smax = new double[0];
		EGS4.radc_bins = new int[0];
		EGS4.radc_ixmin1 = new int[0];
		EGS4.radc_ixmax1 = new int[0];
		EGS4.radc_ixmin2 = new int[0];
		EGS4.radc_ixmax2 = new int[0];
		EGS4.radc_ixmin3 = new int[0];
		EGS4.radc_ixmax3 = new int[0];
		EGS4.radc_ixmin4 = new int[0];
		EGS4.radc_ixmax4 = new int[0];
		EGS4.radc_startx = new int[0];
		EGS4.radc_startb = new int[0];
		EGS4.a_triplet = new double[0][0];
		EGS4.b_triplet = new double[0][0];
		EGS4.nrcp_idata = new int[0][0][0];
		EGS4.nrcp_xdata = new double[0];
		EGS4.nrcp_fdata = new double[0][0][0];
		EGS4.nrcp_wdata = new double[0][0][0];
		EGS4.sig_ismonotone = new boolean[0][0];

	}

	/**
	 * Perform basic initialization and RUN the Monte Carlo engine.
	 */
	private void init() {
		if (createOutputFile) {
			try {
				sigfos = new FileWriter(filename);
			} catch (Exception ex) {
			}
		}

		EGS4.startSimulationTime = System.currentTimeMillis();

		// EGS4.setMXMED(5);//"MAX # OF MEDIA                               "
		// EGS4.setMXREG(50);//"#REGIONS, $MAXRADII*($MAXZPLANE-1)+1(VAC)    "
		// EGS4.setMXSTACK(4000);//"NEED HUGE STACK FOR CORRELATIONS+splitting   "
		// EGS4.setMXRANGE(500); //"for range arrays used in range_rejection()"
		// --variable init
		EGS4Geom.$MAXZREG = 20;// "MAX # OF DOSE SCORING PLANAR ZONES           "
		EGS4Geom.$MAXRADII = 8;// "MAX # OF DOSE SCORING PLANAR ZONES           "
		EGS4SrcEns.$MXRDIST = 1000;
		EGS4SrcEns.$NENSRC = 300;
		EGS4Geom.$NVALUE = 100;
		AMASS = new double[EGS4Geom.$MAXZREG][EGS4Geom.$MAXRADII];// ($MAXZREG,$MAXRADII),
		MEDSAV = new int[EGS4.$MXREG];
		;// ($MXREG),
		RHOSAV = new double[EGS4.$MXREG];
		;// ($MXREG),
		Ecut = new double[EGS4.$MXREG];
		startEcutRegion = new int[EGS4.$MXREG];
		stopEcutRegion = new int[EGS4.$MXREG];
		Pcut = new double[EGS4.$MXREG];
		startPcutRegion = new int[EGS4.$MXREG];
		stopPcutRegion = new int[EGS4.$MXREG];
		Smax = new double[EGS4.$MXREG];
		startSmaxRegion = new int[EGS4.$MXREG];
		stopSmaxRegion = new int[EGS4.$MXREG];
		Incoh = new int[EGS4.$MXREG];
		startIncohRegion = new int[EGS4.$MXREG];
		stopIncohRegion = new int[EGS4.$MXREG];
		Coh = new int[EGS4.$MXREG];
		startCohRegion = new int[EGS4.$MXREG];
		stopCohRegion = new int[EGS4.$MXREG];
		Relax = new int[EGS4.$MXREG];
		startRelaxRegion = new int[EGS4.$MXREG];
		stopRelaxRegion = new int[EGS4.$MXREG];
		Pe = new int[EGS4.$MXREG];
		startPeRegion = new int[EGS4.$MXREG];
		stopPeRegion = new int[EGS4.$MXREG];
		// ---------LOCAL VAR---------------
		SCDOSE = new double[EGS4Geom.$MAXZREG][EGS4Geom.$MAXRADII][$MAXIT];
		SCDOSE2 = new double[EGS4Geom.$MAXZREG][EGS4Geom.$MAXRADII][$MAXIT];
		SCDOSE_COV = new double[EGS4Geom.$MAXZREG][EGS4Geom.$MAXRADII][3];
		SCDOSE_LAST = new int[EGS4Geom.$MAXZREG][EGS4Geom.$MAXRADII];
		SCDOSE_TMP = new double[EGS4Geom.$MAXZREG][EGS4Geom.$MAXRADII][$MAXIT];
		// ($MAXZREG,$MAXRADII,$MAXIT)
		EGS4Grid.CDSTBL = new String[EGS4.$MXREG];
		EGS4Grid.CTRTBL = new String[EGS4.$MXREG];
		EGS4Grid.CABSRB = new String[EGS4.$MXREG];
		EGS4Grid.CAVTRACK = new String[EGS4.$MXREG];
		// --interface->LINK
		EGS4.eq = this;// pass the printing mode
		EGS4Core.eq = this;// pass the printing mode
		EGS4Macro.eq = this;// pass the printing mode
		EGS4SrcEns.eq = this;// pass the printing mode
		EGS4Geom.eq = this;// pass the printing mode
		EGS4Grid.eq = this;// pass the printing mode
		// --Macro and defaults param:
		// EGS4.egs_set_defaults();//first default then:
		// EGS4.RandomUse=1;//use EGS ranlux or ranmar generators
		// EGS4.ranluxB=true;//NOT use ranmar!!
		EGS4.ispmfp = EGS4.iCavity;// select photon mfp
		EGS4.iurd = EGS4.iCavity;// user range discard
		EGS4.iraycorr = EGS4.iCavity;// rayleigh correction
		EGS4.isemfp = EGS4.iCavity;// $SELECT-ELECTRON-MFP
		EGS4.iGeom = EGS4.iCavity;// 1=RZ GEOM
		EGS4Macro.ismfpfpb = EGS4.iCavity;// select mfp parallel beam
		EGS4Macro.irange_rej = EGS4.iCavity;// range rejection

		EGS4Grid.$MAXCMPTS = $MAXCMPTS;
		// -----------------
		EGS4.USER_CONTROLS_TSTEP_RECURSION = EGS4.iCavity;// no effect here
		EGS4.hatchindex = EGS4.iCavity;// no effect here
		// -----------------
		EGS4.iprint = IPRINT;// summary
		// inputs
		is_finished = true;

		// "******************************************************************************
		// "
		// " *** SECTION 1 ***
		// "
		// "------------------------------------------------------------------------------
		// "
		// "READ INPUTS AND CALCULATE ONE-TIME ONLY CONSTANTS
		// "
		// "------------------------------------------------------------------------------

		Calendar cal = Calendar.getInstance();
		Date d = cal.getTime();
		EGS4.seqStr = "=================================================================================";
		EGS4.seqStr = " ********************CAVITY: APPLICATION*****************************************";
		if (EGS4.iprint > 1)
			printSequence(EGS4.seqStr);
		EGS4.seqStr = " SIMULATES THE PASSAGE OF AN ELECTRON OR PHOTON BEAM IN A FINITE, RIGHT CYLINDRICAL GEOMETRY.";
		if (EGS4.iprint > 1)
			printSequence(EGS4.seqStr);
		EGS4.seqStr = " IT IS INTENDED FOR USE IN CALCULATING QUANTITIES OF INTEREST FOR THICK-WALLED ION CHAMBERS ";
		if (EGS4.iprint > 1)
			printSequence(EGS4.seqStr);
		EGS4.seqStr = " EXPOSED TO PHOTON BEAMS. IT ALSO MAY BE USED SIMPLY TO SCORE DOSE IN A CYLINDRICAL GEOMETRY.";
		if (EGS4.iprint > 1)
			printSequence(EGS4.seqStr);
		EGS4.seqStr = " 					" + d.toString();
		if (EGS4.iprint > 1)
			printSequence(EGS4.seqStr);
		EGS4.seqStr = " ********************************************************************************";
		if (EGS4.iprint > 1)
			printSequence(EGS4.seqStr);
		inputs();
		if (EGS4.STOPPROGRAM) {
			closeFile();
			return;
		}
		if (SOURCE == SOURCE_POINT || SOURCE == SOURCE_SARPAGAN
				|| SOURCE == SOURCE_MARINELLI || SOURCE == SOURCE_BEAM) {
			fixSRCOTO();
		}

		// "OPEN A FILE FOR STORING OR READING RANDOM NUMBERS"
		if (ISTORE > 0)// NOT ALLOWED
		{// "We want to store the rng state in a file"
			// rng_unit = egs_open_file(2,0,1,'.egsrns');
		}
		/*
		 * ELSE IF( irestart = 4 )//NOT ALLOWED [ rng_unit =
		 * egs_open_datfile(2,0,1,'.egsrns'); ]
		 */
		/*
		 * IF(IRESTART.EQ.5)//NOT ALLOWED [
		 * 
		 * call egs_combine_runs(combine_results,'.egsdat');
		 * 
		 * NBATCH=0; "DON'T WANT IT TO RUN ANY HISTORIES" NCASET=NCASEO;
		 * "To prevent a wrong normalization if some of the "
		 * "parallel runs not available, IK, Jan 21 1999" ]
		 * "end of IRESTART = 5, DISTRIBUTED POST-PROCESSING" ELSE [
		 */
		if (NCASE / $NBATCH == 0) {
			NCASE = $NBATCH;
		}
		JCASE = NCASE / $NBATCH;
		NCASE = JCASE * $NBATCH;// "NUMBER OF HISTORIES PER BATCH
		// ]
		MXNP = 0; // "reset the maximum stack indicator"
		IHSTRY = NCASEO; // "reset the number of histories counter"

		// "set up the broad parallel beam defaults"
		if (EGS4SrcEns.ISOURC == 2) {
			EGS4Geom.NR = 1;
			EGS4Geom.RCYL[1] = 1000.;
			EGS4Geom.nreg = EGS4Geom.NZ + 1;
			EGS4Geom.CYRAD2[0] = EGS4Geom.RCYL[1] * EGS4Geom.RCYL[1];
		}
		// "set up ausgab calls"
		for (int J = 1; J <= 5; J++) {
			EGS4.iausfl[J - 1] = 1;// IAUSFL(J)=1;
		}
		for (int J = 6; J <= 25; J++) {
			EGS4.iausfl[J - 1] = 0;
		} // "NORMAL EXECUTION"

		if (EGS4Macro.IFULL == 1) {
			// "these flags are the minimum set needed to identify primary and secondary"
			EGS4.iausfl[7] = 1; // "After BREMSSTRAHLUNG"
			EGS4.iausfl[13] = 1; // "After ANNIHILATION IN FLIGHT"
			EGS4.iausfl[14] = 1; // "After ANNIHILATION AT REST"
			EGS4.iausfl[18] = 1; // "After COMPTON"
			EGS4.iausfl[20] = 1; // "After Photo"
			EGS4.iausfl[24] = 1; // "After Rayleigh"
		}
		if (IFANO == 1) {
			// "AUSGAB will be responsible for making sure that the beam is not"
			// "attenuated and getting rid of the scattered photons."
			EGS4.iausfl[15] = 1; // "Before pair"
			EGS4.iausfl[17] = 1; // "Before Compton"
			EGS4.iausfl[18] = 1; // "After Compton"
			EGS4.iausfl[19] = 1; // "Before photoelectric"
			EGS4.iausfl[20] = 1; // "After photoelectric"
			EGS4.iausfl[23] = 1; // "Before Rayleigh"
			EGS4.iausfl[24] = 1; // "After Rayleigh"

			// "AUSGAB will be responsible for throwing away any photons resulting"
			// "from a primary electron. ie. True equilibtrium requires that all"
			// "energy deposition be local to the primary interaction site."
			EGS4.iausfl[7] = 1; // "After bremsstrahlung"
			EGS4.iausfl[13] = 1; // "A positron has annihilated in-flight"
			EGS4.iausfl[14] = 1; // "A positron has annihilated at rest"
		}

		if (EGS4Macro.n_split > 1) {
			EGS4.iausfl[7] = 1; // "After bremsstrahlung"
			EGS4.iausfl[13] = 1; // "A positron has annihilated in-flight"
			EGS4.iausfl[14] = 1; // "A positron has annihilated at rest"

			// "With n_split > 1, we don't need the following calls even if ifano = 1"
			EGS4.iausfl[15] = 0;
			EGS4.iausfl[17] = 0;
			EGS4.iausfl[19] = 0;
			EGS4.iausfl[23] = 0;
			EGS4.iausfl[24] = 0;
		}

		if (IFANO == 2) {

			// "AUSGAB will be responsible for discarding energy due to electrons set"
			// "in motion in the wall"
			EGS4.iausfl[16] = 1; // "After pair"
			EGS4.iausfl[18] = 1; // "After Compton"
			EGS4.iausfl[20] = 1; // "After photoelectric"
		}

		if (EGS4Macro.use_enhance) {
			EGS4.iausfl[15] = 1; // "Before pair"
			EGS4.iausfl[17] = 1; // "Before Compton"
			EGS4.iausfl[18] = 1; // "After Compton"
			EGS4.iausfl[19] = 1; // "Before photoelectric"
			EGS4.iausfl[20] = 1; // "After photoelectric"
			EGS4.iausfl[23] = 1; // "Before Rayleigh"
			EGS4.iausfl[24] = 1; // "After Rayleigh"
			EGS4.iausfl[7] = 1; // "After bremsstrahlung"
			EGS4.iausfl[13] = 1; // "A positron has annihilated in-flight"
			EGS4.iausfl[14] = 1; // "A positron has annihilated at rest"
		}

		for (int j = 1; j <= 28; j++) {
			if (EGS4.iausfl[j - 1] != 0) {// write(6,'(i3,$)') j;
				int jj = j - 1;
				EGS4.seqStr = " AUSGAB index call: " + jj;
				if (EGS4.iprint > 2)
					printSequence(EGS4.seqStr);
			}
		}

		// "HATCH CALL PREPARATION AND EXECUTION"
		EGS4.DUNIT = 1; // "SET LENGTH UNITS TO CMS"
		EGS4.iprint = 0;
		HATCH();
		if (EGS4.STOPPROGRAM) {
			closeFile();
			return;
		}
		EGS4.iprint = IPRINT;

		if (EGS4Macro.irejct == 1) {
			EGS4Macro.initialize_range_rejection();
		}

		if (EGS4SrcEns.MONOEN == 0 && EGS4SrcEns.ISOURC != 21
				&& EGS4SrcEns.ISOURC != 22) {// "MONOENERGETIC INPUT BEAM"
			if (EGS4SrcEns.iqin == 0) {
				EI = EGS4SrcEns.ein;
			} else {
				EI = EGS4SrcEns.ein + EGS4.RM;
			}
			EKMAX = EGS4SrcEns.ein; // "MAXIMUM KINETIC ENERGY"
		} else if (EGS4SrcEns.MONOEN == 1) {// "ENERGY SPECTRUM"
			EGS4SrcEns.ENSRC1();// "NORMALIZE THE ENERGY DISTRIBUTION"
			EKMAX = EGS4SrcEns.ENSRCD[EGS4SrcEns.NENSRC];// "MAXIMUM KINETIC ENERGY IN THE SPECTRUM"
		} else if (EGS4SrcEns.ISOURC == 21 || EGS4SrcEns.ISOURC == 22) {// "phase space input"
																		// EKMAX=EKSRCM;//NOT
																		// ALLOWED
		} else {
			EKMAX = 0;
		}// " <------------ fixme"

		// "CHECK THAT THE DATA FILE HAD DATA OVER THE ENERGY RANGE REQUIRED"
		for (int I = 1; I <= EGS4.NMED; I++) {
			if ((EKMAX > EGS4.UP[I - 1]) || (EKMAX > EGS4.UE[I - 1] - EGS4.RM)) {
				EGS4.STOPPROGRAM = true;
				EGS4.seqStr = " FOR MEDIUM: " + I + "  INCIDENT ENERGY="
						+ EGS4.format(EKMAX, 10, true) + " MeV";
				printSequence(EGS4.seqStr);
				EGS4.seqStr = " IS GREATER THAN COVERED BY DATA FILE WHERE UP,UE="
						+ EGS4.format(EGS4.UP[I - 1], 10, true)
						+ " "
						+ EGS4.format(EGS4.UE[I - 1], 10, true) + " MeV";
				printSequence(EGS4.seqStr);
				EGS4.seqStr = " EXECUTION WILL BE TERMINATED!";
				printSequence(EGS4.seqStr);

				// return;
			}
		}// "END OF LOOP OVER MEDIA"
		if (EGS4.STOPPROGRAM) {
			closeFile();
			return;
		}

		// "CALCULATE THE MASS OF EACH ZONE (AREAL MASS FOR ISOURC=2 OR 4)
		for (int IZ = 1; IZ <= EGS4Geom.NZ; IZ++) {
			DEPTH = EGS4Geom.ZPLANE[IZ] - EGS4Geom.ZPLANE[IZ - 1];
			for (int IX = 1; IX <= EGS4Geom.NR; IX++) {
				// $GET-IRL(IZ,IX);
				IRL = EGS4Geom.GET_IRL(IZ, IX);
				MEDNUM = EGS4.MED[IRL - 1];
				if (MEDNUM != 0) {
					if ((EGS4SrcEns.ISOURC == 2) || (EGS4SrcEns.ISOURC == 4)) {
						VOLUME = DEPTH;
					} else {
						if (IX == 1) {
							RLOW2 = 0.0;
						} else {
							RLOW2 = EGS4Geom.CYRAD2[IX - 2];
						}
						VOLUME = Math.PI * DEPTH
								* (EGS4Geom.CYRAD2[IX - 1] - RLOW2);
					}
					AMASS[IZ - 1][IX - 1] = EGS4.RHOR[IRL - 1] * VOLUME;
				} else {
					AMASS[IZ - 1][IX - 1] = 0.0;
				}
			}// "END OF IX LOOP"
		}// "END OF IZ LOOP"

		// "CALCULATE ONE-TIME-ONLY CONSTANTS FOR SOURCE"
		if (SOURCE == -1)
			EGS4SrcEns.SRCOTO();// (WEIGHT);
		if ((EGS4SrcEns.IFPB == 0) && (EGS4Macro.IFORCE != 0)
				&& (EGS4SrcEns.iqin == 0) && (EGS4SrcEns.MONOEN == 0)) {
			EGS4Macro.SELECT_MEAN_FREE_PATHS_FOR_FRONTAL_PARALLEL_BEAM();
		}

		// "INITIALIZE DATA ARRAYS FOR FLUORESCENT X-RAYS IF NEEDED"
		ISUMX = 0;
		for (int JJ = 1; JJ <= EGS4Geom.nreg; JJ++) {
			ISUMX = ISUMX + EGS4.iedgfl[JJ - 1]; // "NON-ZERO IF X-RAYS ANYWHERE"
		}
		if (ISUMX != 0) {
			EGS4.EDGSET(2, EGS4Geom.nreg);
		}
		// "NOTE THE ABOVE WILL PRODUCE LOTS OF EXTRA OUTPUT AND SHOULD BE"
		// "CLEANED UP"
		ISUMRY(); // "PRINT THE SUMMARY OF INPUTS"

		// "******************************************************************************
		// "
		// " *** SECTION 2 ***
		// "
		// "------------------------------------------------------------------------------
		// "
		// "LOOP THROUGH THE NUMBER OF HISTORIES. CALCULATE CONSTANTS THAT MAY
		// CHANGE FOR
		// "EACH HISTORY AND DO THE SIMULATION
		// "
		// "------------------------------------------------------------------------------

		// "WRITE THE HEADER"
		// WRITE(IOUT,100) TITLE; call egs_fdate(iout); write(iout,*);
		EGS4.seqStr = "=========================================================================";
		if (EGS4.iprint > 1)
			printSequence(EGS4.seqStr);
		EGS4.seqStr = "                   EXECUTION INFORMATION";
		if (EGS4.iprint > 1)
			printSequence(EGS4.seqStr);
		EGS4.seqStr = "=========================================================================";
		if (EGS4.iprint > 1)
			printSequence(EGS4.seqStr);

		// WRITE(IOUT,200);WRITE(6,200); "PRINT HEADER FOR EXECUTION MESSAGES"

		// "PRINT EXECUTION MODE"
		// if(IRESTART == 0){WRITE(6,201);WRITE(IOUT,201);}
		// ELSEIF(IRESTART.EQ.1)[
		// WRITE(6,202) NCASE,NCASEO;
		// write(6,'(21x,a,$)') 'New RNG state: ';
		// $SHOW-RNG-STATE(6); write(6,*);
		// write(iout,'(21x,a,$)') 'New RNG state: ';
		// $SHOW-RNG-STATE(iout); write(iout,*);
		// ]
		// ELSEIF(IRESTART.EQ.3)[WRITE(6,204);WRITE(IOUT,204);GO TO :END-SIM:;]
		// ELSEIF(IRESTART.EQ.4)[WRITE(6,205);WRITE(IOUT,205);]
		// ELSEIF(IRESTART.EQ.5)[WRITE(6,206);WRITE(IOUT,206);GO TO :END-SIM:;]

		// "Initialize IWATCH routine"
		if (IWATCH != 0)
			EGS4.WATCH(-99, IWATCH);

		// "SET CLOCK AT THE BEGINNING OF SIMULATIONS"
		// $INITIALIZE_ELAPSED_CPU_TIME;
		// $SET_ELAPSED_CPUTIME(CPUT1);
		// $INITIALIZE_ELAPSED_TOTAL_TIME;
		// ETIMETOT=0;
		// TIMEB=0;
		// NETADJ=0;

		// "dcav_old = 0.0;"

		// "Calculate the cavity mass for in-flight dose display"
		FMASSC = 0.0; // "TOTAL CAVITY MASS"
		for (int IX = 1; IX <= EGS4Geom.NR; IX++) {
			for (int IZ = 1; IZ <= EGS4Geom.NZ; IZ++) {
				// $GET-IRL(IZ,IX);
				IRL = EGS4Geom.GET_IRL(IZ, IX);
				if (EGS4Geom.ntrack[IRL - 1] == 1) {
					FMASS = AMASS[IZ - 1][IX - 1];
					FMASSC = FMASSC + FMASS; // "SUM THE CAVITY MASS USED LATER"
				}
			}
		}
		// "end of calculation of cavity mass"

		// "Initialize variables for the cs_enhance scoring "
		tmp_dose = 0;
		tmp_dose1 = 0;
		last_case = 0;
		EGS4SrcEns.NHSTRY = 0;

		// "Open file for data storage, if requested "
		// "The file is opened in the temporary working directory"
		// IF( idat = 0 ) data_unit = egs_open_file(4,0,1,'.egsdat');

		n_photons = 0;
		n_electrons = 0;
		sumE_photons = 0.0;
		sumE2_photons = 0.0;
		sumE_electrons = 0.0;
		sumE2_electrons = 0.0;

		// "Output batches. Statistical analysis is done after each batch.
		// Execution
		// "stops if the desired statistical accuracy is obtained or there is
		// not enough
		// "time to do another batch.
		boolean ENDSIM = false;
		EGS4SrcEns.AINFLU = 0.0;
		// iieff=new int[$NBATCH];pondt_batch=new double[$NBATCH];
		for (int IBATCH = 1; IBATCH <= $NBATCH; IBATCH++) {
			ENDSIM = false;
			long startTime = System.currentTimeMillis();
			IBTCH = IBATCH;
			if (IBATCH == 1) {
				EGS4.seqStr = " BATCH" + EGS4.format("", 8) + "ELAPSED"
						+ EGS4.format("", 3) + "time of day"
						+ EGS4.format("", 2) + "cavity stats(%)"
						+ EGS4.format("", 2) + "cav.dose(Gy.cm^2)";
				if (EGS4.iprint > 0)
					printSequence(EGS4.seqStr);
			} else {// " not first batch"
			}// " end of before batch ne 1 block"

			for (int ICASE = 1; ICASE <= JCASE; ICASE++) {// "now fill each IS bin"

				if (EGS4SrcEns.ISOURC != 23)
					IHSTRY = IHSTRY + 1; // "increment history counter"

				EGS4Macro.NFTIME = 0; // "reset the photon forced interaction counter"

				// CALL SRCHST(XIN,YIN,ZIN,UIN,VIN,WIN,IRIN,WEIGHT,NRCFLG);
				if (SOURCE == -1)
					EGS4SrcEns.SRCHST();
				// "calculate the initial energy if a distribution is to be used"
				if (EGS4SrcEns.MONOEN != 0) {// "if equal to 0, it is monoenergetic"
												// EGS4EnsSrc.ENSRCH(EIN);
												// "returns K.E. from distribution"
					EGS4SrcEns.ein = EGS4SrcEns.ENSRCH();
					if (EGS4SrcEns.iqin == 0) {
						EI = EGS4SrcEns.ein;
					} else {
						EI = EGS4SrcEns.ein + EGS4.RM;
					}// "total energy"
					// " there was a check that the data file had data over the
					// energy
					// "range required, the location of it will eventually be in
					// "ESRCIN.MOR
				}
				// "Set photon weights if gamma interactions are to be forced in
				// the
				// "target in the frontal parallel beam case if monoenergetic
				if ((EGS4SrcEns.MONOEN == 0) && (EGS4SrcEns.iqin == 0)
						&& (EGS4Macro.IFORCE == 1) && (EGS4SrcEns.IFPB == 0)
						&& (EGS4SrcEns.ISOURC != 21)
						&& (EGS4SrcEns.ISOURC != 22)) {
					int IX = (EGS4SrcEns.irin - 2) / EGS4Geom.NZ + 1;
					EGS4Macro.GWAIT = EGS4Macro.GWATE[IX - 1];
					EGS4SrcEns.WEIGHT = EGS4Macro.GWAIT;
				}

				// "FOR AN INPUT ENERGY SPECTRUM, DETAILED FORCING MACRO IS USED"

				EGS4.LATCHI = 0;
				// ****************************FIX WEIGHT
				// *****************************************************
				if (SOURCE == SOURCE_POINT || SOURCE == SOURCE_SARPAGAN
						|| SOURCE == SOURCE_MARINELLI || SOURCE == SOURCE_BEAM) {
					fixEmAll();
				}
				// ********************************************************************************************

				if ((IWATCH != 0) && (IWATCH != 4)) {
					EGS4.seqStr = EGS4.format(1, 5)
							+ EGS4.format(EGS4SrcEns.ein, 9, true)
							+ EGS4.format(EGS4SrcEns.iqin, 4)
							+ EGS4.format(EGS4SrcEns.irin, 4)
							+ EGS4.format(EGS4SrcEns.xin, 8, true)
							+ EGS4.format(EGS4SrcEns.yin, 8, true)
							+ EGS4.format(EGS4SrcEns.zin, 8, true)
							+ EGS4.format(EGS4SrcEns.uin, 8, true)
							+ EGS4.format(EGS4SrcEns.vin, 8, true)
							+ EGS4.format(EGS4SrcEns.win, 8, true)
							+ EGS4.format(EGS4.LATCHI, 10)
							+ EGS4.format(EGS4SrcEns.WEIGHT, 10, false);
					if (EGS4.iprint > 1)
						printSequence(EGS4.seqStr);
					// OUTPUT
					// 1,EIN,IQIN,IRIN,XIN,YIN,ZIN,UIN,VIN,WIN,LATCHI,WEIGHT;
					// (/' INITIAL SHOWER VALUES',T36,':',
					// I5,F9.3,2I4,3F8.3,3F7.3,I10,1PE10.3);
				}
				// "ALL INITIAL SHOWER VARIABLES ARE SET, CALL THE SHOWER ROUTINE"
				if (EGS4SrcEns.iqin == 0) {
					n_photons = n_photons + 1;
					sumE_photons = sumE_photons + EI;
					sumE2_photons = sumE2_photons + EI * EI;
				} else {
					n_electrons = n_electrons + 1;
					sumE_electrons = sumE_electrons + EI - EGS4.RM;
					sumE2_electrons = sumE2_electrons + (EI - EGS4.RM)
							* (EI - EGS4.RM);
				}

				// CALL SHOWER(IQIN,EI,XIN,YIN,ZIN,UIN,VIN,WIN,IRIN,WEIGHT);
				if (EGS4SrcEns.iqin != 0) {
					if (EGS4SrcEns.ein > ETHRESHOLD) {
						ekinokB = true;
					} else {
						ekinokB = false;
					}
				}
				indet = true;// force
				if (indet && ekinokB)
					SHOWER();
				if (EGS4.STOPPROGRAM) {
					closeFile();
					return;
				}

				if (IWATCH > 0)
					EGS4.WATCH(-1, IWATCH);

			}// "END OF THE ICASE LOOP"

			// "******************"
			// "history by history"
			// "    EMH March 2002"
			// "******************"
			if (EGS4SrcEns.ISOURC == 21 || EGS4SrcEns.ISOURC == 22) {
				// AINFLU_CURRENT=
				// EGS4SrcEns.dble(NNREAD+NRCYCL*(NNREAD-IHSTRY))/EGS4SrcEns.dble(NCASE_PHSP)*NINCSRC;
			} else if (EGS4SrcEns.ISOURC == 23) {
				AINFLU_CURRENT = IHSTRY;
			} else {
				AINFLU_CURRENT = EGS4SrcEns.AINFLU * EGS4SrcEns.dble(IHSTRY)
						/ EGS4SrcEns.dble(EGS4SrcEns.NCASET);
			}
			dcav_current = cav_dose * 1.602E-10 / (FMASSC * AINFLU_CURRENT);
			// double mevtojoule=1.60218E-13;//1MeV=1,60218 · 10-13 J
			// and MASS is in g not KG!!!=>E-10!!

			dcavun = (cav2_dose * IHSTRY - cav_dose * cav_dose) / (IHSTRY - 1);
			if (dcavun > 0) {
				dcavun = Math.sqrt(dcavun);
			}
			;
			dcavun = dcavun * 1.602E-10 / (FMASSC * AINFLU_CURRENT);
			dcavun = 100 * dcavun / dcav_current;

			// #####################################################################################
			String timePerBatch = EGS4.timeElapsedShort(startTime);

			Calendar call = Calendar.getInstance();
			String timeday = call.get(Calendar.HOUR) + ":" + call.get(Calendar.MINUTE)
					+ ":" + call.get(Calendar.SECOND);

			EGS4.seqStr = EGS4.format("", 2) + EGS4.format(IBATCH, 3)
					+ EGS4.format("", 2) + EGS4.format(timePerBatch, 14)
					+ EGS4.format("", 3) + EGS4.format(timeday, 9)
					+ EGS4.format("", 6) + EGS4.format(dcavun, 7, true)
					+ EGS4.format("", 10)
					+ EGS4.format(dcav_current, 11, false);
			if (EGS4.iprint > 0)
				printSequence(EGS4.seqStr);

			if (dcavun <= STATLM && STATLM > 0.) {
				// "we have reached the desired statistics, print a message and exit"
				EGS4.seqStr = " DESIRED STATISTICAL ACCURACY OBTAINED.";
				if (EGS4.iprint > 1)
					printSequence(EGS4.seqStr);
				EGS4.seqStr = " STATS IN CAVITY= "
						+ EGS4.format(dcavun, 5, true) + "%" + " AFTER "
						+ EGS4.format(IBTCH, 2) + " BATCHES";
				if (EGS4.iprint > 1)
					printSequence(EGS4.seqStr);

				// WRITE(6,230)dcavun,IBTCH;WRITE(IOUT,230)dcavun,IBTCH;
				// GO TO :END-SIM:;
				ENDSIM = true;
				break;
			}
		}// "END OF SIMULATIONS"

		// "PRINT INSUFFICIENT STATS WARNING"
		if (!ENDSIM) {
			EGS4.seqStr = " DESIRED STATISTICAL ACCURACY OF "
					+ EGS4.format(STATLM, 5, true) + "%" + " NOT REACHED!";
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);
			EGS4.seqStr = " STATS IN CAVITY= " + EGS4.format(dcavun, 5, true)
					+ " %" + " AFTER " + EGS4.format(IBTCH, 2) + " BATCHES";
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);

			// WRITE(IOUT,240) STATLM,DCAVUN,IBTCH;WRITE(6,240)
			// STATLM,DCAVUN,IBTCH;
		}

		// :END-SIM:;

		SCSTP = SCSTP + SCSTP_TMP;
		SCSTP2 = SCSTP2 + SCSTP_TMP * SCSTP_TMP;
		SCCSTP = SCCSTP + SCCSTP_TMP;
		SCCSTP2 = SCCSTP2 + SCCSTP_TMP * SCCSTP_TMP;

		// "******************"
		// "history by history"
		// "    EMH March 2002"
		// "******************"
		// "we have to add the temporary scoring"
		// "variables for the last history      "
		cav_dose = cav_dose + tmp_dose;
		cav2_dose = cav2_dose + tmp_dose * tmp_dose;

		cav_dose0 = cav_dose0 + tmp_dose0;
		cav2_dose0 = cav2_dose0 + tmp_dose0 * tmp_dose0;

		cav_dose1 = cav_dose1 + tmp_dose1;
		cav2_dose1 = cav2_dose1 + tmp_dose1 * tmp_dose1;

		cav_dose2 = cav_dose2 + tmp_dose2;
		cav2_dose2 = cav2_dose2 + tmp_dose2 * tmp_dose2;

		cav_dosec = cav_dosec + tmp_dose * tmp_dose1;
		cav_dosec01 = cav_dosec01 + tmp_dose0 * tmp_dose1;
		cav_dosec02 = cav_dosec02 + tmp_dose0 * tmp_dose2;

		if (EGS4Geom.NSUMCV > 1) {// "store data for individual cavity regions"

			for (int IZ = 1; IZ <= EGS4Geom.NZ; IZ++) {
				for (int IX = 1; IX <= EGS4Geom.NR; IX++) {
					// $GET-IRL(IZ,IX);
					IRL = EGS4Geom.GET_IRL(IZ, IX);

					if (EGS4Geom.ntrack[IRL - 1] == 1) {
						for (int IT = 1; IT <= 4; IT++) {
							SCDOSE[IZ - 1][IX - 1][IT - 1] = SCDOSE[IZ - 1][IX - 1][IT - 1]
									+ SCDOSE_TMP[IZ - 1][IX - 1][IT - 1];
							SCDOSE2[IZ - 1][IX - 1][IT - 1] = SCDOSE2[IZ - 1][IX - 1][IT - 1]
									+ SCDOSE_TMP[IZ - 1][IX - 1][IT - 1]
									* SCDOSE_TMP[IZ - 1][IX - 1][IT - 1];
						}
						SCDOSE_COV[IZ - 1][IX - 1][0] = SCDOSE_COV[IZ - 1][IX - 1][0]
								+ SCDOSE_TMP[IZ - 1][IX - 1][0]
								* SCDOSE_TMP[IZ - 1][IX - 1][2];
						SCDOSE_COV[IZ - 1][IX - 1][1] = SCDOSE_COV[IZ - 1][IX - 1][1]
								+ SCDOSE_TMP[IZ - 1][IX - 1][1]
								* SCDOSE_TMP[IZ - 1][IX - 1][2];
						SCDOSE_COV[IZ - 1][IX - 1][2] = SCDOSE_COV[IZ - 1][IX - 1][2]
								+ SCDOSE_TMP[IZ - 1][IX - 1][1]
								* SCDOSE_TMP[IZ - 1][IX - 1][3];
					}
				}
			}
		}
		if (SOURCE == -1) {
			EGS4.seqStr = " FINAL RANDOM NUMBER STATE:";
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);
			EGS4.SHOW_RNG_STATE();
		}
		// "******************************************************************************
		// "
		// " *** SECTION 3 ***
		// "
		// "------------------------------------------------------------------------------
		//
		// "STATISTICAL AND OTHER DATA HANDLING AND CALL THE OUTPUT SUMMARY ROUTINE"
		//
		// "------------------------------------------------------------------------------

		// :STATS-ANAL:;

		FMASSC = 0.0; // "total cavity mass"
		for (int IX = 1; IX <= EGS4Geom.NR; IX++) {
			for (int IZ = 1; IZ <= EGS4Geom.NZ; IZ++)// DO IZ=1,NZ
			{
				// $GET-IRL(IZ,IX);
				IRL = EGS4Geom.GET_IRL(IZ, IX);
				if (EGS4Geom.ntrack[IRL - 1] == 1) {
					FMASS = AMASS[IZ - 1][IX - 1];
					FMASSC = FMASSC + FMASS;
				}
			}
		}

		EGS4SrcEns.AINFLU = EGS4SrcEns.AINFLU * EGS4SrcEns.dble(IHSTRY)
				/ EGS4SrcEns.dble(EGS4SrcEns.NCASET);
		SCORE_NORM_NUM = EGS4SrcEns.dble(IHSTRY);

		// "******************"
		// "history by history"
		// "    EMH March 2002"
		// "******************"

		// $ANALYZE(SCOMEG, :SCORE_NORM_NUM);
		SCORE_TEMP = EGS4SrcEns.SCOMEG / SCORE_NORM_NUM;
		EGS4SrcEns.SCOMEG2 = EGS4SrcEns.SCOMEG2 / SCORE_NORM_NUM;
		EGS4SrcEns.SCOMEG2 = (EGS4SrcEns.SCOMEG2 - SCORE_TEMP * SCORE_TEMP)
				/ (SCORE_NORM_NUM - 1);
		if (EGS4SrcEns.SCOMEG2 >= 0.)
			EGS4SrcEns.SCOMEG2 = Math.sqrt(EGS4SrcEns.SCOMEG2);
		if (SCORE_TEMP != 0.) {
			EGS4SrcEns.SCOMEG2 = Math.min(EGS4SrcEns.SCOMEG2 / SCORE_TEMP
					* 100., 99.9);
		} else {
			EGS4SrcEns.SCOMEG2 = 99.9;
		}

		EGS4SrcEns.SCOMEG = EGS4SrcEns.SCOMEG / EGS4SrcEns.dble(IHSTRY);// "Corrected, IK May 4 1999"

		// OUTPUT SCOMEG,SCOMEG2;(/' OMEG =',1PE12.3,'(',0PF5.1,'%)'/);
		if (SOURCE == -1) {
			EGS4.seqStr = " OMEG =" + EGS4.format(EGS4SrcEns.SCOMEG, 12, false)
					+ "(" + EGS4.format(EGS4SrcEns.SCOMEG2, 5) + "%)";
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);
		}

		// $ANALYZE(SCSTP, :SCORE_NORM_NUM);
		SCORE_TEMP = SCSTP / SCORE_NORM_NUM;
		SCSTP2 = SCSTP2 / SCORE_NORM_NUM;
		SCSTP2 = (SCSTP2 - SCORE_TEMP * SCORE_TEMP) / (SCORE_NORM_NUM - 1);
		if (SCSTP2 >= 0.)
			SCSTP2 = Math.sqrt(SCSTP2);
		if (SCORE_TEMP != 0.) {
			SCSTP2 = Math.min(SCSTP2 / SCORE_TEMP * 100., 99.9);
		} else {
			SCSTP2 = 99.9;
		}

		// $ANALYZE(SCCSTP, :SCORE_NORM_NUM);
		SCORE_TEMP = SCCSTP / SCORE_NORM_NUM;
		SCCSTP2 = SCCSTP2 / SCORE_NORM_NUM;
		SCCSTP2 = (SCCSTP2 - SCORE_TEMP * SCORE_TEMP) / (SCORE_NORM_NUM - 1);
		if (SCCSTP2 >= 0.)
			SCCSTP2 = Math.sqrt(SCCSTP2);
		if (SCORE_TEMP != 0.) {
			SCCSTP2 = Math.min(SCCSTP2 / SCORE_TEMP * 100., 99.9);
		} else {
			SCCSTP2 = 99.9;
		}

		// "first, estimate uncertainties for entire cavity"
		cav2_dose = (cav2_dose * SCORE_NORM_NUM - cav_dose * cav_dose)
				/ (SCORE_NORM_NUM - 1);
		if (cav2_dose > 0)
			cav2_dose = Math.sqrt(cav2_dose);

		if (EGS4Macro.IFULL == 1) {

			cav2_dose0 = (cav2_dose0 * SCORE_NORM_NUM - cav_dose0 * cav_dose0)
					/ (SCORE_NORM_NUM - 1);
			if (cav2_dose0 > 0)
				cav2_dose0 = Math.sqrt(cav2_dose0);
			cav2_dose1 = (cav2_dose1 * SCORE_NORM_NUM - cav_dose1 * cav_dose1)
					/ (SCORE_NORM_NUM - 1);
			if (cav2_dose1 > 0)
				cav2_dose1 = Math.sqrt(cav2_dose1);
			cav2_dose2 = (cav2_dose2 * SCORE_NORM_NUM - cav_dose2 * cav_dose2)
					/ (SCORE_NORM_NUM - 1);
			if (cav2_dose2 > 0)
				cav2_dose2 = Math.sqrt(cav2_dose2);

			corr_02 = (cav_dosec02 * SCORE_NORM_NUM - cav_dose0 * cav_dose2)
					/ (SCORE_NORM_NUM - 1);
			corr_02 = cav2_dose0 * cav2_dose0 + cav2_dose2 * cav2_dose2 + 2
					* corr_02;
			if (corr_02 > 0)
				corr_02 = Math.sqrt(corr_02);

			cav_dosec = (cav_dosec * SCORE_NORM_NUM - cav_dose * cav_dose1)
					/ (SCORE_NORM_NUM - 1);
			cav_dosec01 = (cav_dosec01 * SCORE_NORM_NUM - cav_dose0 * cav_dose1)
					/ (SCORE_NORM_NUM - 1);
			cav_dosec02 = (cav_dosec02 * SCORE_NORM_NUM - cav_dose0 * cav_dose2)
					/ (SCORE_NORM_NUM - 1);
		}

		cav_dose = cav_dose * 1.602e-10 / (EGS4SrcEns.AINFLU * FMASSC);

		cav_dose0 = cav_dose0 * 1.602e-10 / (EGS4SrcEns.AINFLU * FMASSC);
		cav_dose1 = cav_dose1 * 1.602e-10 / (EGS4SrcEns.AINFLU * FMASSC);
		cav_dose2 = cav_dose2 * 1.602e-10 / (EGS4SrcEns.AINFLU * FMASSC);

		cav2_dose = cav2_dose * 1.602e-10 / (EGS4SrcEns.AINFLU * FMASSC);

		cav2_dose0 = cav2_dose0 * 1.602e-10 / (EGS4SrcEns.AINFLU * FMASSC);
		cav2_dose1 = cav2_dose1 * 1.602e-10 / (EGS4SrcEns.AINFLU * FMASSC);
		cav2_dose2 = cav2_dose2 * 1.602e-10 / (EGS4SrcEns.AINFLU * FMASSC);

		corr_02 = corr_02 * 1.602e-10 / (EGS4SrcEns.AINFLU * FMASSC);

		cav_dosec = cav_dosec * (1.602e-10 / (EGS4SrcEns.AINFLU * FMASSC))
				* (1.602e-10 / (EGS4SrcEns.AINFLU * FMASSC));
		cav_dosec01 = cav_dosec01 * (1.602e-10 / (EGS4SrcEns.AINFLU * FMASSC))
				* (1.602e-10 / (EGS4SrcEns.AINFLU * FMASSC));
		cav_dosec02 = cav_dosec02 * (1.602e-10 / (EGS4SrcEns.AINFLU * FMASSC))
				* (1.602e-10 / (EGS4SrcEns.AINFLU * FMASSC));

		if (EGS4Geom.NSUMCV > 1) {// "multiple cavity regions, analyze quantities in each region"

			// "FOR ISOURC=4 WE NEED THE DATA FOR CIRCLES, NOT RINGS, SO ADD IT UP"
			// "THIS SHOULD ONLY BE USED IF THE CAVITY HAS AN INFINITE DIAMETER"
			if ((EGS4SrcEns.ISOURC == 4) && (EGS4Geom.NR > 1)) {
				for (int IX = 1; IX <= EGS4Geom.NR; IX++)// DO IX=1,NR[
				{
					for (int IZ = 1; IZ <= EGS4Geom.NZ; IZ++)// DO IZ=1,NZ[
					{
						// $GET-IRL(IZ,IX);
						IRL = EGS4Geom.GET_IRL(IZ, IX);
						if (EGS4Geom.ntrack[IRL - 1] == 1) {
							for (int IT = 1; IT <= 4; IT++) {
								// "IK: this must be a bug. for ix=1 ix-1=0 and the"
								// "    array is not defined"
								SCDOSE[IZ - 1][IX - 1][IT - 1] = SCDOSE[IZ - 1][IX - 1][IT - 1]
										+ SCDOSE[IZ - 1][IX - 2][IT - 1];
								SCDOSE2[IZ - 1][IX - 1][IT - 1] = SCDOSE2[IZ - 1][IX - 1][IT - 1]
										+ SCDOSE2[IZ - 1][IX - 2][IT - 1];
								if (IT < 4) {
									SCDOSE_COV[IZ - 1][IX - 1][IT - 1] = SCDOSE_COV[IZ - 1][IX - 1][IT - 1]
											+ SCDOSE_COV[IZ - 1][IX - 2][IT - 1];
								}
							}
						}
					}
				}
			}

			for (int IX = 1; IX <= EGS4Geom.NR; IX++)// DO IX=1,NR[
			{
				for (int IZ = 1; IZ <= EGS4Geom.NZ; IZ++)// DO IZ=1,NZ[
				{
					// $GET-IRL(IZ,IX);
					IRL = EGS4Geom.GET_IRL(IZ, IX);
					if (EGS4Geom.ntrack[IRL - 1] == 1)// NTRACK(IRL).EQ.1)[
					{
						SCDOSE2[IZ - 1][IX - 1][0] = (SCDOSE2[IZ - 1][IX - 1][0]
								* SCORE_NORM_NUM - SCDOSE[IZ - 1][IX - 1][0]
								* SCDOSE[IZ - 1][IX - 1][0])
								/ (SCORE_NORM_NUM - 1);
						if (SCDOSE2[IZ - 1][IX - 1][0] > 0)
							SCDOSE2[IZ - 1][IX - 1][0] = Math
									.sqrt(SCDOSE2[IZ - 1][IX - 1][0]);
						if (EGS4Macro.IFULL == 1) {
							SCDOSE2[IZ - 1][IX - 1][1] = (SCDOSE2[IZ - 1][IX - 1][1]
									* SCORE_NORM_NUM - SCDOSE[IZ - 1][IX - 1][1]
									* SCDOSE[IZ - 1][IX - 1][1])
									/ (SCORE_NORM_NUM - 1);
							if (SCDOSE2[IZ - 1][IX - 1][1] > 0)
								SCDOSE2[IZ - 1][IX - 1][1] = Math
										.sqrt(SCDOSE2[IZ - 1][IX - 1][1]);
							SCDOSE2[IZ - 1][IX - 1][2] = (SCDOSE2[IZ - 1][IX - 1][2]
									* SCORE_NORM_NUM - SCDOSE[IZ - 1][IX - 1][2]
									* SCDOSE[IZ - 1][IX - 1][2])
									/ (SCORE_NORM_NUM - 1);
							if (SCDOSE2[IZ - 1][IX - 1][2] > 0)
								SCDOSE2[IZ - 1][IX - 1][2] = Math
										.sqrt(SCDOSE2[IZ - 1][IX - 1][2]);
							SCDOSE2[IZ - 1][IX - 1][3] = (SCDOSE2[IZ - 1][IX - 1][3]
									* SCORE_NORM_NUM - SCDOSE[IZ - 1][IX - 1][3]
									* SCDOSE[IZ - 1][IX - 1][3])
									/ (SCORE_NORM_NUM - 1);
							if (SCDOSE2[IZ - 1][IX - 1][3] > 0)
								SCDOSE2[IZ - 1][IX - 1][3] = Math
										.sqrt(SCDOSE2[IZ - 1][IX - 1][3]);
							// "now calculate the covariances"

							SCDOSE_COV[IZ - 1][IX - 1][0] = (SCDOSE_COV[IZ - 1][IX - 1][0]
									* SCORE_NORM_NUM - SCDOSE[IZ - 1][IX - 1][0]
									* SCDOSE[IZ - 1][IX - 1][2])
									/ (SCORE_NORM_NUM - 1);

							SCDOSE_COV[IZ - 1][IX - 1][1] = (SCDOSE_COV[IZ - 1][IX - 1][1]
									* SCORE_NORM_NUM - SCDOSE[IZ - 1][IX - 1][1]
									* SCDOSE[IZ - 1][IX - 1][2])
									/ (SCORE_NORM_NUM - 1);

							SCDOSE_COV[IZ - 1][IX - 1][2] = (SCDOSE_COV[IZ - 1][IX - 1][2]
									* SCORE_NORM_NUM - SCDOSE[IZ - 1][IX - 1][1]
									* SCDOSE[IZ - 1][IX - 1][3])
									/ (SCORE_NORM_NUM - 1);
						}

						// "now normalize quantities and convert to dose"
						FMASS = AMASS[IZ - 1][IX - 1];

						SCDOSE[IZ - 1][IX - 1][0] = SCDOSE[IZ - 1][IX - 1][0]
								* 1.602e-10 / (EGS4SrcEns.AINFLU * FMASS);
						SCDOSE[IZ - 1][IX - 1][1] = SCDOSE[IZ - 1][IX - 1][1]
								* 1.602e-10 / (EGS4SrcEns.AINFLU * FMASS);
						SCDOSE[IZ - 1][IX - 1][2] = SCDOSE[IZ - 1][IX - 1][2]
								* 1.602e-10 / (EGS4SrcEns.AINFLU * FMASS);
						SCDOSE[IZ - 1][IX - 1][3] = SCDOSE[IZ - 1][IX - 1][3]
								* 1.602e-10 / (EGS4SrcEns.AINFLU * FMASS);

						SCDOSE2[IZ - 1][IX - 1][0] = SCDOSE2[IZ - 1][IX - 1][0]
								* 1.602e-10 / (EGS4SrcEns.AINFLU * FMASS);
						SCDOSE2[IZ - 1][IX - 1][1] = SCDOSE2[IZ - 1][IX - 1][1]
								* 1.602e-10 / (EGS4SrcEns.AINFLU * FMASS);
						SCDOSE2[IZ - 1][IX - 1][2] = SCDOSE2[IZ - 1][IX - 1][2]
								* 1.602e-10 / (EGS4SrcEns.AINFLU * FMASS);
						SCDOSE2[IZ - 1][IX - 1][3] = SCDOSE2[IZ - 1][IX - 1][3]
								* 1.602e-10 / (EGS4SrcEns.AINFLU * FMASS);

						SCDOSE_COV[IZ - 1][IX - 1][0] = SCDOSE_COV[IZ - 1][IX - 1][0]
								* (1.602e-10 / (EGS4SrcEns.AINFLU * FMASS))
								* (1.602e-10 / (EGS4SrcEns.AINFLU * FMASS));
						SCDOSE_COV[IZ - 1][IX - 1][1] = SCDOSE_COV[IZ - 1][IX - 1][1]
								* (1.602e-10 / (EGS4SrcEns.AINFLU * FMASS))
								* (1.602e-10 / (EGS4SrcEns.AINFLU * FMASS));
						SCDOSE_COV[IZ - 1][IX - 1][2] = SCDOSE_COV[IZ - 1][IX - 1][2]
								* (1.602e-10 / (EGS4SrcEns.AINFLU * FMASS))
								* (1.602e-10 / (EGS4SrcEns.AINFLU * FMASS));
					}
				}
			}
		}

		OSUMRY(); // "PRINT THE OUTPUT SUMRY"

		// :END-OF-RUN:;

		// ;"******************************************************************************
		// "
		// " *** SECTION 4 ***
		// "
		// "------------------------------------------------------------------------------
		// "
		// "THE CONCLUSION"
		// "
		// "------------------------------------------------------------------------------

		// :END:;
		// OUTPUT; (//'END OF RUN',9X,' ',$); call egs_fdate(6);
		// OUTPUT; (//);
		// write(iout,'(/a,$)') 'END OF RUN '; call egs_fdate(iout);
		// write(iout,'(////)');

		// call egs_finish;

		// #ifdef HAVE_C_COMPILER;
		// ;
		// IF( n_parallel > 0 & ~is_finished ) [
		// call egs_pjob_finish(n_job);
		// / IF( n_job = 0 ) [
		// is_finished = .true.;
		// call egs_combine_runs(combine_results,'.egsdat');
		// NCASET=NCASEO; IHSTRY=NCASET;
		// CALL SRCOTO(WEIGHT);
		// goto :STATS-ANAL:;
		// ]
		// ]
		// #endif;

		EGS4SrcEns.SRCEND();// do nothing

		// $CALL_EXIT(0);

		// "FORTRAN FORMAT STATEMENTS. FORMAT STATEMENT N## IS FIRST USED IN SECTION N."
		// %I0
		// 100 FORMAT(80A1//'Calculation using CAVRZnrc(EGSnrc) '$VERSION' ',
		// /' ON '$MACHINE' ',T55,' ',$);
		// 200 FORMAT(//,79('*')/
		// // ,T20,'EXECUTION INFORMATION AND WARNING MESSAGES'/
		// // ,79('*')/
		// //'USING CAVRZnrc(EGSnrc) '$VERSION' ON '$MACHINE' ');
		// 201 FORMAT(/'********* NEW INPUT FILE *********'/);
		// 202 FORMAT(/'********* RESTARTED INPUT FILE ********* '/
		// ' ',10X,I12,' NEW + ',I12,' OLD HISTORIES');
		// 204 FORMAT(/'********* DATA ANALYSIS ONLY *********'/);
		// 205 FORMAT(/'********* RANDOM NUMBERS READ FROM FILE *********'/);
		// 206 FORMAT(/' ********* ANALYZING RESULTS FROM PARALLEL RUNS
		// *******'/);
		// 210 FORMAT(/'********* NOT ENOUGH TIME TO FINISH WITHIN',
		// ' LIMIT OF',F8.2,' HOURS',I5,' BATCHES USED********'/
		// ' ',I12,' HISTORIES RUN, ',I12,' HISTORIES ANALYZED'//);
		// 230 FORMAT(/'DESIRED STATISTICAL ACCURACY OBTAINED.'/
		// ' STATS IN CAVITY= ',F5.2,'%',
		// ' AFTER ',I2,' BATCHES');
		// 240 FORMAT(/'*********DESIRED STATISTICAL ACCURACY OF ',F5.2,'%',
		// ' NOT REACHED*********'/
		// ' STATS IN CAVITY= ',F5.2,' % AFTER ',I2,' BATCHES');
		// 250 FORMAT(/' FOR OLD RUN:'/
		// ' ----------- '/
		// ' Total cputime =',F8.1,'s (=',F5.2,' hr)');
		// 255 FORMAT(/' FOR PARALLEL RUNS:'/
		// ' ----------------- '/
		// ' On ',I5,' machines '/
		// ' Total cputime =',F8.1,'s (=',F8.2,' hr), cputime/machine =',
		// F8.1,'s');
		// 260 FORMAT(/'Finished simulations:'/' time elapsed,cputime',
		// ',ratio= ',2F8.1,'(=',F5.2,'hr)',F8.2);
		// 261 FORMAT(/' Finished: time elapsed this run', F10.1/
		// ' CPUtime total run ', F10.1,'(=',F8.2,'hr)'/
		// ' Ratio ELAPSED/CPU this run:', F8.3);
		// 280 FORMAT(/' CPUtime/history=',F10.5,' sec. Histories/hour=',F12.0);
		//
		// END; "END OF MAIN ROUTINE-CAVRZnrc"

		destroyArrays();
		EGS4.timeElapsed();
		Calendar call = Calendar.getInstance();
		Date da = call.getTime();
		EGS4.seqStr = "End of run:          " + da.toString();
		if (EGS4.iprint > 0)
			printSequence(EGS4.seqStr);

		// ===============================================================
		closeFile();
	}

	/**
	 * Close the file containing the simulation result.
	 */
	public void closeFile() {
		if (createOutputFile) {
			try {
				sigfos.close();
			} catch (Exception ex) {
			}

			putInFile = false;
			EGS4.seqStr = "Check current directory for results, in text file:"
					+ filename;
			if (EGS4.iprint > 0)
				printSequence(EGS4.seqStr);
		}
	}

	/**
	 * Where to print runtime information. Interface method.
	 * @param s the String to be printed
	 */
	public void printSequence(String s) {
		// write file?
		if (createOutputFile && putInFile) {
			try {
				sigfos.write(s + " \n");
			} catch (Exception ex) {
			}
		}

		// output data to console
		if (systemOutB)
			System.out.println(s);
		else {
			jta.append(s + " \n");
			jta.selectAll();
			// int lc=jta.getRows();
			// jta.setSelectionEnd(lc);
		}

	}

	// "******************************************************************************
	// "
	// "
	// " **********
	// " * *
	// " * AUSGAB *
	// " * *
	// " **********
	// "
	// "
	// " An AUSGAB routine to be used with cavrznrc.mortran
	// "
	// " This routine scores the dose and other ionisation cavity parameters
	// " in a finite, azimuthally symmetric cylindrical geometry which the
	// " user defines via plane and radius coordinates. The user must specify
	// " both the target geometry as well as the planes and radii between
	// " which the quantities are to be scored. All the geometrical checks for
	// " crossing 'geometrical' or 'dose' regions are handled by the subroutine
	// " HOWFAR.
	// "
	// " FOR IT = 1 the total primary dose
	// " = 2 the total dose dose less the total primary dose
	// " (i.e. the scatter fraction for calculating Ascat)
	// " = 3 the total primary unattenuated dose less the
	// " total primary dose (for calculating Aatt)
	// " = 4 the total primary unattenuated dose with the source
	// " distribution factored out less the total primary
	// " unattenuated dose (for calculating Apn)
	// " = 5 the total primary unattenuated dose to the cavity
	// " gas with the source distribution factored out with
	// " the electron tranport taking place in the chamber
	// " with the cavity filled with wall material (for
	// " calculating Afl and the stopping power ratio)
	// " = 6 the total primary unattenuated dose to the wall
	// " material with the source distribution factored out with
	// " the electron tranport taking place in the chamber
	// " with the cavity filled with wall material (for
	// " calculating the stopping power ratio)
	// "
	// "
	// ;"******************************************************************************

	/**
	 * In general, AUSGAB is a routine which is called under a series 
	 * of well defined conditions specified by the value of IARG. Interface method.
	 * @param IARG the process code
	 */
	public void AUSGAB(int IARG) {

		// $IMPLICIT-NONE;

		double FTMP = 0.;
		int ip = 0;

		// $INTEGER IRL,IX,IZ,IQL,LATCHL,IARG,IDUMMY;
		// $REAL WTL,FDUMMY,xsi;

		// ;COMIN/
		// ELECIN,EPCONT,GEOM,MEDIA,PHOTIN,RUSROU,SCORE,SOURCE,STACK,USEFUL,USER,
		// RANDOM,BOUNDS/;
		double xsi = 0.0;
		// "STACK OVERFLOW CHECK"
		MXNP = Math.max(MXNP, EGS4.NP);// "keep track of how deep stack is"
		// "MXNP is not output but it should be"

		if (IWATCH > 0) {
			EGS4.WATCH(IARG, IWATCH);
		} // "signal watch routine if active"]

		// "check if particle is leaving the transport geometry"
		int IRL = EGS4.IR[EGS4.NP - 1]; // "local region number"
		if (IRL == 1)
			return; // "outside the chamber"

		// "OBTAIN FREQUENTLY USED LOCAL VARIABLES"
		// $GET-IX-IZ(IRL); //"local plane and radius numbers"
		int IX = EGS4Geom.GET_IX(IRL);
		int IZ = EGS4Geom.GET_IZC(IRL);

		int IQL = EGS4.IQ[EGS4.NP - 1]; // "local charge variable"
		double WTL = EGS4.WT[EGS4.NP - 1]; // "local weight variable"
		int LATCHL = EGS4.LATCH[EGS4.NP - 1]; // "LATCHL=0 for primaries, 1 otherwise"
		// #########################################################
		if (EGS4Macro.use_enhance || EGS4Macro.n_split > 1) {
			if (IARG < 5) {
				if (EGS4.EDEP > 0 && WTL > 0 && EGS4Geom.ntrack[IRL - 1] == 1) {
					// "If we use cross section enhancement or photon splitting,"
					// "all energy scoring is done here and the rest of ausgab is ignored"
					// "We use the technique proposed by the PENELOPE group for scoring "
					// "the energy deposition. This results in a much better estimate   "
					// "of the uncertainty                                              "

					FTMP = WTL * EGS4.EDEP;
					if (EGS4SrcEns.NHSTRY == last_case) {
						// " Still the same history scoring into the cavity => update    "
						// " temporary variables                                         "
						if (LATCHL != 2)
							tmp_dose = tmp_dose + FTMP;
						if (LATCHL != 3)
							tmp_dose1 = tmp_dose1 + FTMP;
					} else {
						// " A new history scoring into the cavity. "
						last_case = EGS4SrcEns.NHSTRY;
						cav_dose = cav_dose + tmp_dose;
						cav2_dose = cav2_dose + tmp_dose * tmp_dose;
						cav_dose1 = cav_dose1 + tmp_dose1;
						cav2_dose1 = cav2_dose1 + tmp_dose1 * tmp_dose1;
						cav_dosec = cav_dosec + tmp_dose * tmp_dose1;
						if (LATCHL != 2) {
							tmp_dose = FTMP;
						} else {
							tmp_dose = 0.;
						}
						if (LATCHL != 3) {
							tmp_dose1 = FTMP;
						} else {
							tmp_dose1 = 0.;
						}
					}
				}
				return;
			}
		}
		// #########################################################
		if (EGS4Macro.use_enhance) {// "If we use cross section enhancement, all scoring "
									// " is done here and the rest of ausgab is ignored  "

			if (IARG == 15 || IARG == 17 || IARG == 19 || IARG == 23) {
				// "A pair/Compton/photoelectric/Rayleigh event is about to take place"
				// "As we have increased the photon cross section by a factor of      "
				// "cs_enhance, we must split the photon into a scattering portion    "
				// "(1/cs_enhance) and a nor-scattering portion (1-1/cs_enhance)      "
				// "Start with placing an identical photon on the stack               "
				EGS4.NP = EGS4.NP + 1;
				if (EGS4.NP + 1 > EGS4.$MXSTACK) {
					EGS4.STOPPROGRAM = true;
					EGS4.seqStr = " ***************************************************"
							+ "  \n"
							+ " Calculation with CS-enhancement: unable to boost stack."
							+ "  \n"
							+ " ***************************************************";
					// if(EGS4.iprint>0)
					printSequence(EGS4.seqStr);
					return;
					// OUTPUT;
					// ( ' Calculation with CS-enhancement: unable to boost
					// stack.'/
					// ' Stopping.'/ 1x,80('*')/);
					// stop;
				}
				// $TRANSFER PROPERTIES TO (np) FROM (np - 1);
				EGS4.X[EGS4.NP - 1] = EGS4.X[EGS4.NP - 2];
				EGS4.Y[EGS4.NP - 1] = EGS4.Y[EGS4.NP - 2];
				EGS4.Z[EGS4.NP - 1] = EGS4.Z[EGS4.NP - 2];
				EGS4.IR[EGS4.NP - 1] = EGS4.IR[EGS4.NP - 2];
				EGS4.WT[EGS4.NP - 1] = EGS4.WT[EGS4.NP - 2];
				EGS4.DNEAR[EGS4.NP - 1] = EGS4.DNEAR[EGS4.NP - 2];
				EGS4.LATCH[EGS4.NP - 1] = EGS4.LATCH[EGS4.NP - 2];

				EGS4.E[EGS4.NP - 1] = EGS4.E[EGS4.NP - 2];
				EGS4.U[EGS4.NP - 1] = EGS4.U[EGS4.NP - 2];
				EGS4.V[EGS4.NP - 1] = EGS4.V[EGS4.NP - 2];
				EGS4.W[EGS4.NP - 1] = EGS4.W[EGS4.NP - 2];
				EGS4.IQ[EGS4.NP - 1] = EGS4.IQ[EGS4.NP - 2];
				if (EGS4.LATCH[EGS4.NP - 2] != 2) {
					// " This is either a primary photon that has not yet been attenuated "
					// " away or a scattered photon. Let's decide what to do with the     "
					// " unscattered fraction of that photon (which is at np-1)           "
					xsi = EGS4.random01();
					if (EGS4Macro.cs_enhance * xsi < 1.) {// " The photon doesn't survive. "
						if (EGS4.LATCH[EGS4.NP - 2] == 3) {// " It is a scattered photon => kill it"
							EGS4.WT[EGS4.NP - 2] = 0.0;
							EGS4.DNEAR[EGS4.NP - 2] = -1.;
						} else {// " This is a primary => mark it as attenuated.         "
								// " From now on, all descendents of this photon will    "
								// " only contribute to the cavity dose with attenuation "
								// " and scatter removed                                 "
							EGS4.LATCH[EGS4.NP - 2] = 2;
						}
					}
				}
				// "Adjust the weight of to be scattered photon"
				EGS4.WT[EGS4.NP - 1] = EGS4.WT[EGS4.NP - 1]
						/ EGS4Macro.cs_enhance;
				return;
			}

			if (IARG == 18 || IARG == 20 || IARG == 24 ||
			// " A Compton/photo-absorption/Rayleigh event just occured"
					IARG == 7 || IARG == 13 || IARG == 14) {
				// " A bremas/annihilation event just occured"
				// " All scattered photons and photons originating in brems/annihilation"
				// " events contribute to the scattered dose. But because all of them   "
				// " have a small weight (initial weight/cs_enhance), we will play      "
				// " Russian Roulette with them, using 1/cs_enhance as a sirvivng       "
				// " probability. If they survive, their weight will become equal to the"
				// " intial photon weight. In addition, we have to set their latch to 3 "
				// " so that they and rheir descendents only contribute to the scattered"
				// " dose.                                                              "
				xsi = EGS4.random01();
				xsi = xsi * EGS4Macro.cs_enhance;
				for (ip = EGS4.NPold; ip <= EGS4.NP; ip++) {
					if (EGS4.IQ[ip - 1] == 0) {
						if (EGS4.LATCH[ip - 1] == 2) {// "that's a descendent of a photon that"
														// "has been attenuated away => kill it"
							EGS4.WT[ip - 1] = 0.0;
							EGS4.DNEAR[ip - 1] = -1.;
						} else {
							if (EGS4.E[ip - 1] >= EGS4.PCUT[IRL - 1]) {
								if (xsi < 1.) {
									EGS4.LATCH[ip - 1] = 3;
									EGS4.WT[ip - 1] = EGS4.WT[ip - 1]
											* EGS4Macro.cs_enhance;
								} else {
									EGS4.WT[ip - 1] = 0.;
									EGS4.DNEAR[ip - 1] = -1.;
								}
							} else {
								EGS4.LATCH[ip - 1] = 3;
							}
							// " i.e. we don't need the Russian Roulette for photons below"
							// " threshold because they will be discarded and their energy"
							// " deposited locally anyway                                 "
						}
					}
				}
				return;
			}
			return;
		}
		// #########################################################

		if (EGS4Macro.n_split > 1) {

			if (IARG == 7 || IARG == 13 || IARG == 14) {
				if (EGS4Macro.iifano == 1) {
					for (ip = EGS4.NPold; ip <= EGS4.NP; ip++) {
						if (EGS4.IQ[ip - 1] == 0) {
							EGS4.WT[ip - 1] = 0.;
							EGS4.E[ip - 1] = 0.;
						}
					}
					return;
				}
				for (ip = EGS4.NPold; ip <= EGS4.NP; ip++)// DO ip=NPold,NP [
				{
					if (EGS4.IQ[ip - 1] == 0) {
						xsi = EGS4.random01();
						if (xsi * EGS4Macro.n_split > 1) {
							EGS4.WT[ip - 1] = 0.;
							EGS4.E[ip - 1] = 0.;
						} else {
							EGS4.LATCH[ip - 1] = 3;
							EGS4.WT[ip - 1] = EGS4.WT[ip - 1]
									* EGS4Macro.n_split;
						}
					}
				}
				return;
			}
		}
		// #########################################################

		// REPLACE {$SCORE(#,#:#)} WITH {;
		// "Scoring macro used in AUSGAB for quantities other than DOSE and KERMA"
		// "{P1}{P2}=scoring array (eg SCSTP)"
		// "{P3}=quantity to be scored (eg 1)"

		// "If the (primary) history number, NHSTRY, is the same as the history"
		// "that last scored in this array, {P1}_LAST{P2}, then {P3} is added"
		// "to a temporary array, {P1}_TMP{P2}.  Otherwise, we add"
		// "{P1}_TMP{P2} to {P1}{P2}, {P1}_TMP{P2}*{P1}_TMP{P2} to {P1}2{P2},"
		// "set {P1}_TMP{P2}={P3}, and set {P1}_LAST{P2}=NHSTRY."
		// "This scoring method allows us to calculate  uncorrelated value"
		// "of {P1}2{P2} which is then used to calculate the uncertainty"
		// "in {P1}{P2}.  In cavrznrc, this macro is only used for counting"
		// "no. of charged particle steps."

		// IF(NHSTRY={P1}_LAST{P2})[
		// {P1}_TMP{P2}={P1}_TMP{P2} + {P3};
		// ]
		// ELSE[
		// {P1}{P2}={P1}{P2}+{P1}_TMP{P2};
		// {P1}2{P2}={P1}2{P2} + {P1}_TMP{P2}*{P1}_TMP{P2};
		// {P1}_TMP{P2}={P3};
		// {P1}_LAST{P2}=NHSTRY;
		// ]
		// ;
		// }

		if (IARG == 0) {// "ABOUT TO TRANSPORT A PARTICLE"
			if (IQL != 0) {
				if (LATCHL == 0) {// "COUNT PRIMARY CHARGED PARTICLES ONLY"
									// $SCORE(SCSTP,
									// :1);"COUNT CHARGED PARTICLE STEPS TAKEN"
					if (EGS4SrcEns.NHSTRY == SCSTP_LAST)// IF(NHSTRY={P1}_LAST{P2})[
					{
						// {P1}_TMP{P2}={P1}_TMP{P2} + {P3};
						SCSTP_TMP = SCSTP_TMP + 1.;
					} else {
						// {P1}{P2}={P1}{P2}+{P1}_TMP{P2};
						SCSTP = SCSTP + SCSTP_TMP;
						// {P1}2{P2}={P1}2{P2} + {P1}_TMP{P2}*{P1}_TMP{P2};
						SCSTP2 = SCSTP2 + SCSTP_TMP * SCSTP_TMP;
						// {P1}_TMP{P2}={P3};
						SCSTP_TMP = 1.;
						// {P1}_LAST{P2}=NHSTRY;
						SCSTP_LAST = EGS4SrcEns.NHSTRY;
					}

					if (EGS4Geom.ntrack[IRL - 1] == 1) {
						// $SCORE(SCCSTP, :1);
						if (EGS4SrcEns.NHSTRY == SCCSTP_LAST)// IF(NHSTRY={P1}_LAST{P2})[
						{
							// {P1}_TMP{P2}={P1}_TMP{P2} + {P3};
							SCCSTP_TMP = SCCSTP_TMP + 1.;
						} else {
							// {P1}{P2}={P1}{P2}+{P1}_TMP{P2};
							SCCSTP = SCCSTP + SCCSTP_TMP;
							// {P1}2{P2}={P1}2{P2} + {P1}_TMP{P2}*{P1}_TMP{P2};
							SCCSTP2 = SCCSTP2 + SCCSTP_TMP * SCCSTP_TMP;
							// {P1}_TMP{P2}={P3};
							SCCSTP_TMP = 1.;
							// {P1}_LAST{P2}=NHSTRY;
							SCCSTP_LAST = EGS4SrcEns.NHSTRY;
						}

						// "WRITE(*,*)' sccstp ',SCCSTP;"
					}
				}// "COUNT STEPS IN CAVITY REGION"
			} else {// "PHOTON STEP - PLAY RUSSIAN ROULETTE?"
				if (RUSROU && (EGS4.W[EGS4.NP - 1] > 0.0)) {// "YES, PLAY IF CROSSES RRZ "
					if ((EGS4.Z[EGS4.NP - 1] <= RRZ)
							&& (EGS4.Z[EGS4.NP - 1] + EGS4.USTEP
									* EGS4.W[EGS4.NP - 1] >= RRZ)) {// "CROSSES"
						xsi = EGS4.random01();
						if (xsi < RRCUT) {// "PARTICLE SURVIVES"
							EGS4.WT[EGS4.NP - 1] = WTL / RRCUT;
						} else {// "DISCARD PARTICLE ON NEXT CALL TO HOWFAR"
							EGS4.WT[EGS4.NP - 1] = 0.0;
						}
					} // "END TEST IF CROSSES RUSSIAN ROULETTE PLANE"
				} // "END TEST FOR PLAYING RUSSIAN ROULETTE"
			}// "END TEST FOR PHOTON STEP"
		}// "END TEST FOR IARG = 0"

		if (IFANO == 1) {
			if (IARG == 15 || IARG == 17 || IARG == 19 || IARG == 23) {
				// "A pair/Compton/photoelectric/Rayleigh event is about to take place"
				EGS4.NP = EGS4.NP + 1; // "Boost the stack"
				if (EGS4.NP + 1 > EGS4.$MXSTACK) {
					EGS4.STOPPROGRAM = true;
					EGS4.seqStr = " ***************************************************"
							+ "  \n"
							+ " Fano calculation unable to boost stack."
							+ "  \n"
							+ " ***************************************************";
					// if(EGS4.iprint>0)
					printSequence(EGS4.seqStr);
					return;
				}
				// "Create an identical photon"
				// $TRANSFER PROPERTIES TO (np) FROM (np - 1);
				EGS4.X[EGS4.NP - 1] = EGS4.X[EGS4.NP - 2];
				EGS4.Y[EGS4.NP - 1] = EGS4.Y[EGS4.NP - 2];
				EGS4.Z[EGS4.NP - 1] = EGS4.Z[EGS4.NP - 2];
				EGS4.IR[EGS4.NP - 1] = EGS4.IR[EGS4.NP - 2];
				EGS4.WT[EGS4.NP - 1] = EGS4.WT[EGS4.NP - 2];
				EGS4.DNEAR[EGS4.NP - 1] = EGS4.DNEAR[EGS4.NP - 2];
				EGS4.LATCH[EGS4.NP - 1] = EGS4.LATCH[EGS4.NP - 2];

				EGS4.E[EGS4.NP - 1] = EGS4.E[EGS4.NP - 2];
				EGS4.U[EGS4.NP - 1] = EGS4.U[EGS4.NP - 2];
				EGS4.V[EGS4.NP - 1] = EGS4.V[EGS4.NP - 2];
				EGS4.W[EGS4.NP - 1] = EGS4.W[EGS4.NP - 2];
				EGS4.IQ[EGS4.NP - 1] = EGS4.IQ[EGS4.NP - 2];

				return;
			}

			// "Throw away any scattered photons from the primary interaction site."
			// " Now there is a stack pointer NPold which points to the particle "
			// " befor the last discrete interaction. This change was necessary "
			// " for the implementation of atomic relaxations. So, check all particles"
			// " between NPold and NP and discard photons "
			// "IF ( (iarg = 18 & NP > NPold)" " Compton has occured"
			if (IARG == 18// " Compton has occured"
					|| IARG == 20 // " After photo-absorption "
					|| IARG == 24 // " After Rayleigh "
					|| IARG == 7 // " After brems "
					|| IARG == 13 // " After annihilation "
					|| IARG == 14)// " After annihilation at rest "
			{
				for (ip = EGS4.NPold; ip <= EGS4.NP; ip++) {
					if (EGS4.IQ[ip - 1] == 0) {
						EGS4.WT[ip - 1] = 0.;
						EGS4.E[ip - 1] = 0.;
					}
				}
			}
		}
		if (IFANO == 2) {
			if (IARG == 16 // " After pair production "
					|| IARG == 18 // " After Compton "
					|| IARG == 20) // " After photo absorption "
			{
				if (EGS4Geom.ntrack[EGS4.IR[EGS4.NP - 1] - 1] == 0) {
					for (ip = EGS4.NPold; ip <= EGS4.NP; ip++)// DO ip=NPold,NP
					{
						if (EGS4.IQ[ip - 1] == 0)// if( iq(ip) ~= 0 )
						{
							EGS4.WT[ip - 1] = 0.;
							EGS4.E[ip - 1] = 0.;
						}// { wt(ip) = 0; e(ip) = 0; }
					}
				}
			}
		}

		// "SCORE THE ENERGY AS REQUIRED FOR THE DIFFERENT MODES"
		// "****************************************************"

		if ((EGS4.EDEP != 0.0) && (WTL > 0.0) && (IARG < 5)
				&& (EGS4Geom.ntrack[IRL - 1] == 1)) {
			// "ENERGY HAS BEEN DEPOSITED IN THE CAVITY REGION"
			// "SCORE PRIMARY AND SECONDARY ENERGY DEPOSITED"
			FTMP = EGS4.WT[EGS4.NP - 1] * EGS4.EDEP;
			// "***************************************************************************"
			// "                                                                           "
			// " Implementation of a history by history scoring scheme for the cavity dose "
			// " taken from Iwan's splitting implementation, EMH, March 2002               "
			// "                                                                           "
			// " For the meaning of the variables see the COMIN/SCORE block                "
			// "                                                                           "
			// "         expmfp = exp(di) - 1 (see macro $SELECT-PHOTON-MFP)               "
			// "                                                                           "
			// "***************************************************************************"
			if (EGS4SrcEns.NHSTRY == last_case) {
				// " Still the same history scoring into the cavity => update    "
				// " temporary variables                                         "
				tmp_dose = tmp_dose + FTMP;
				if (LATCHL == 0) {
					tmp_dose0 = tmp_dose0 + FTMP;
					tmp_dose1 = tmp_dose1 + FTMP * (EGS4Macro.EXPMFP + 1);
				} else {
					tmp_dose2 = tmp_dose2 + FTMP;
				}
			} else {
				// " A new history scoring into the cavity. "
				last_case = EGS4SrcEns.NHSTRY;

				cav_dose = cav_dose + tmp_dose;
				cav2_dose = cav2_dose + tmp_dose * tmp_dose;

				cav_dose0 = cav_dose0 + tmp_dose0;
				cav2_dose0 = cav2_dose0 + tmp_dose0 * tmp_dose0;

				cav_dose1 = cav_dose1 + tmp_dose1;
				cav2_dose1 = cav2_dose1 + tmp_dose1 * tmp_dose1;

				cav_dose2 = cav_dose2 + tmp_dose2;
				cav2_dose2 = cav2_dose2 + tmp_dose2 * tmp_dose2;

				cav_dosec = cav_dosec + tmp_dose * tmp_dose1;
				cav_dosec01 = cav_dosec01 + tmp_dose0 * tmp_dose1;
				cav_dosec02 = cav_dosec02 + tmp_dose0 * tmp_dose2;

				tmp_dose = FTMP;

				if (LATCHL == 0) {
					tmp_dose0 = FTMP;
					tmp_dose1 = FTMP * (EGS4Macro.EXPMFP + 1);
					tmp_dose2 = 0.0;
				} else {
					tmp_dose0 = 0.0;
					tmp_dose1 = 0.0;
					tmp_dose2 = FTMP;
				}
			}

			IDECAV = 1;

			if (EGS4Geom.NSUMCV > 1) {// "calculate quantities in individual cavity regions"
										// "do it the same as for the overall cavity above"

				if (EGS4SrcEns.NHSTRY == SCDOSE_LAST[IZ - 1][IX - 1]) {
					SCDOSE_TMP[IZ - 1][IX - 1][0] = SCDOSE_TMP[IZ - 1][IX - 1][0]
							+ FTMP;
					if (LATCHL == 0) {// "primary dose"
						SCDOSE_TMP[IZ - 1][IX - 1][1] = SCDOSE_TMP[IZ - 1][IX - 1][1]
								+ FTMP;
						SCDOSE_TMP[IZ - 1][IX - 1][2] = SCDOSE_TMP[IZ - 1][IX - 1][2]
								+ FTMP * (1 + EGS4Macro.EXPMFP);
					} else {// "secondary dose"
						SCDOSE_TMP[IZ - 1][IX - 1][3] = SCDOSE_TMP[IZ - 1][IX - 1][3]
								+ FTMP;
					}
				} else {
					SCDOSE_LAST[IZ - 1][IX - 1] = EGS4SrcEns.NHSTRY;

					SCDOSE[IZ - 1][IX - 1][0] = SCDOSE[IZ - 1][IX - 1][0]
							+ SCDOSE_TMP[IZ - 1][IX - 1][0];
					SCDOSE2[IZ - 1][IX - 1][0] = SCDOSE2[IZ - 1][IX - 1][0]
							+ SCDOSE_TMP[IZ - 1][IX - 1][0]
							* SCDOSE_TMP[IZ - 1][IX - 1][0];

					SCDOSE[IZ - 1][IX - 1][1] = SCDOSE[IZ - 1][IX - 1][1]
							+ SCDOSE_TMP[IZ - 1][IX - 1][1];
					SCDOSE2[IZ - 1][IX - 1][1] = SCDOSE2[IZ - 1][IX - 1][1]
							+ SCDOSE_TMP[IZ - 1][IX - 1][1]
							* SCDOSE_TMP[IZ - 1][IX - 1][1];

					SCDOSE[IZ - 1][IX - 1][2] = SCDOSE[IZ - 1][IX - 1][2]
							+ SCDOSE_TMP[IZ - 1][IX - 1][2];
					SCDOSE2[IZ - 1][IX - 1][2] = SCDOSE2[IZ - 1][IX - 1][2]
							+ SCDOSE_TMP[IZ - 1][IX - 1][2]
							* SCDOSE_TMP[IZ - 1][IX - 1][2];

					SCDOSE[IZ - 1][IX - 1][3] = SCDOSE[IZ - 1][IX - 1][3]
							+ SCDOSE_TMP[IZ - 1][IX - 1][3];
					SCDOSE2[IZ - 1][IX - 1][3] = SCDOSE2[IZ - 1][IX - 1][3]
							+ SCDOSE_TMP[IZ - 1][IX - 1][3]
							* SCDOSE_TMP[IZ - 1][IX - 1][3];

					SCDOSE_COV[IZ - 1][IX - 1][0] = SCDOSE_COV[IZ - 1][IX - 1][0]
							+ SCDOSE_TMP[IZ - 1][IX - 1][0]
							* SCDOSE_TMP[IZ - 1][IX - 1][2];
					SCDOSE_COV[IZ - 1][IX - 1][1] = SCDOSE_COV[IZ - 1][IX - 1][1]
							+ SCDOSE_TMP[IZ - 1][IX - 1][1]
							* SCDOSE_TMP[IZ - 1][IX - 1][2];
					SCDOSE_COV[IZ - 1][IX - 1][2] = SCDOSE_COV[IZ - 1][IX - 1][2]
							+ SCDOSE_TMP[IZ - 1][IX - 1][1]
							* SCDOSE_TMP[IZ - 1][IX - 1][3];

					SCDOSE_TMP[IZ - 1][IX - 1][0] = FTMP;

					if (LATCHL == 0) {
						SCDOSE_TMP[IZ - 1][IX - 1][1] = FTMP;
						SCDOSE_TMP[IZ - 1][IX - 1][2] = FTMP
								* (1 + EGS4Macro.EXPMFP);
						SCDOSE_TMP[IZ - 1][IX - 1][3] = 0.0;
					} else {
						SCDOSE_TMP[IZ - 1][IX - 1][1] = 0.0;
						SCDOSE_TMP[IZ - 1][IX - 1][2] = 0.0;
						SCDOSE_TMP[IZ - 1][IX - 1][3] = FTMP;
					}
				}

			}

		} // "END OF ENERGY DEPOSITED IN THE CAVITY"

		// "SET FLAG FOR SECONDARY INTERACTIONS"
		// "***********************************"

		if ((EGS4Macro.IFULL > 0) && (IARG > 5) && (LATCHL == 0)) {
			// "ONLY IF PRIMARY PARTICLES HAVE INTERACTED DISCRETELY                 "
			// "IF A SECONDARY PARTICLE IS CREATED ON THE SECOND PASS, GIVE IT A ZERO"
			// "WEIGHT SO THAT HOWFAR WILL DISCARD IT.                               "

			if (IARG == 7) {// "brem has occured"

				if (IQL == 0) {
					EXCHANGE_STACK(EGS4.NP, EGS4.NP - 1);
				}
				EGS4.LATCH[EGS4.NP - 2] = 1; // " Flag the photon as a secondary"
				if (ipass >= 1) {
					EGS4.WT[EGS4.NP - 2] = 0.;
				}// "To save time in correlation runs"

			} else if (IARG == 18) {// "Compton has occured, with binding effects"
									// "taken into account, 0, 1, or more particles"
									// "may have resulted"

				for (ip = EGS4.NPold; ip <= EGS4.NP; ip++)// DO ip=NPold,NP [
				{
					if (EGS4.IQ[ip - 1] == 0) {
						EGS4.LATCH[ip - 1] = 1;
						if (ipass >= 1) {
							EGS4.WT[ip - 1] = 0.;
						}// "To save time"
					}
				}
				// " NP = NPold means the interactiopn has been rejected and thus "
				// " the emerging (unscattered) photon  is still a primary "
			} else if (IARG == 9) {// "Moller has occured. For now there is only"
									// "one secondary. When impact ionization is implemented"
									// "the following should be changed"
				if (EGS4.E[EGS4.NP - 1] < EGS4.E[EGS4.NP - 2]) {
					EXCHANGE_STACK(EGS4.NP, EGS4.NP - 1);
				}
			} else if (IARG == 13 || IARG == 14) {// "Annihilation, flag the photons"
				EGS4.LATCH[EGS4.NP - 1] = 1;
				EGS4.LATCH[EGS4.NP - 2] = 1;
				if (ipass >= 1) {
					EGS4.WT[EGS4.NP - 1] = 0.;
					EGS4.WT[EGS4.NP - 2] = 0.;
				}
			} else if (IARG == 20) {
				for (ip = EGS4.NPold; ip <= EGS4.NP; ip++)// DO ip=NPold,NP [
				{
					if (EGS4.IQ[ip - 1] == 0) {
						EGS4.LATCH[ip - 1] = 1;
						if (ipass >= 1) {
							EGS4.WT[ip - 1] = 0.;
							EGS4.E[ip - 1] = 0.;
						}
					}
				}
			} else if (IARG == 24) {
				EGS4.LATCH[EGS4.NP - 1] = 1;
				if (ipass >= 1) {
					EGS4.WT[EGS4.NP - 1] = 0.;
				}
			}
		}

		return;
	}// "END OF AUSGAB"

	// "The following is the $CALL-HOWNEAR macro for PRESTA-II

	// REPLACE {$CALL-HOWNEAR(#);} WITH {
	// ;
	// "write(6,'(2i3,4e15.8)') np,ir(np),dnear(np), "
	// "   sqrt(x(np)*x(np)+y(np)*y(np)),z(np),tustep; "
	// IF( dnear(np) < tustep ) [
	// call hownear({P1},x(np),y(np),z(np),ir(np));
	// " write(6,*) ' --> new dnear: ',{P1}; "
	// ]
	// ELSE [ {P1} = dnear(np); ]
	// }
	// "*********************************************************************"
	/**
	 * The following is a general specification of HOWNEAR: 
	 * Given a particle at (x,y,z) in region irl, HOWNEAR answers the 
	 * question, What is the distance tperp to the closest boundary? Interface method.
	 */
	public void HOWNEAR() {
		// "Subroutine arguments
		// $REAL
		// tperp, "nearest distance to any boundary (output)
		// tustep,
		// x, "x-position of the particle (input)
		// y, "y-position of the particle (input)
		// z; "z-position of the particle (input)

		// $INTEGER
		// ir "region number of the particle

		// "Local variables
		double r = 0.0;
		int ix = 0;// "current cylindrical radius number
		int iz = 0;// "current planar slab number

		double z = EGS4.Z[EGS4.NP - 1];
		double y = EGS4.Y[EGS4.NP - 1];
		double x = EGS4.X[EGS4.NP - 1];
		int ir = EGS4.IR[EGS4.NP - 1];
		// FROM MACRO=>ELSE [ {P1} = dnear(np); ]->not necessary, kind of
		// varred(skip calc.)!@@@@@@@@@@@
		if (EGS4.DNEAR[EGS4.NP - 1] >= EGS4.TUSTEP) {
			EGS4.tperp = EGS4.DNEAR[EGS4.NP - 1];
			return;
		}
		// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
		ix = (ir - 2) / EGS4Geom.NZ + 1;// NZ=planar zones=NPLANE-1!!!NR=nr of
										// Radii
		iz = ir - 1 - EGS4Geom.NZ * (ix - 1);
		r = Math.sqrt(x * x + y * y);
		EGS4.tperp = EGS4.min(z - EGS4Geom.ZPLANE[iz - 1], EGS4Geom.ZPLANE[iz]
				- z, EGS4Geom.RCYL[ix] - r);
		if (ix != 1) {
			EGS4.tperp = Math.min(EGS4.tperp, r - EGS4Geom.RCYL[ix - 1]);
		}

		// "IF( tperp < -1e-6 ) [
		// " OUTPUT IHSTRY, tperp; (/' Error in HOWNEAR: IHSTRY=',I12,
		// " ' tperp negative =',1PE10.3);
		// " OUTPUT ir,x,y,z,r; (' ir x y z r=', I5, 4F15.7);
		// " OUTPUT iz,ix;(' depth,radial regions iz,ix=',2I10);
		// " OUTPUT zplane(iz),zplane(iz+1);( ' upper/lower z planes=',2F15.7);
		// " IF( ix > 1) [
		// " OUTPUT rcyl(ix-1),rcyl(ix); (' Radial boundaries=',2F15.7);
		// " ]
		// " ELSE [
		// " OUTPUT rcyl(ix); (' Radial boundary of 1st region =', F15.7);
		// " ]
		// "]

		return;

	}// "end of subroutine HOWNEAR"

	// ;"******************************************************************************
	// "
	// " **********
	// " * *
	// " * HOWFAR *
	// " * *
	// " **********
	// "
	// " A GENERAL PURPOSE CYLINDRICAL GEOMETRY ROUTINE FOR USE WITH THE EGS4
	// " CODE SYSTEM ADAPTED FOR USE WITH CAVRZnrc.
	// "
	// " FOR PARTICLE NP ON THE STACK IN REGION IR(NP), THIS ROUTINE
	// " DETERMINES IF THE PARTICLE CAN GO A DISTANCE USTEP WITHOUT CHANGING
	// " ZONES. IF USTEP CAUSES A ZONE CROSSING, IT IS REDUCED TO PLACE IT ON
	// " THE BOUNDRY AND IRNEW IS SET TO THE ZONE NUMBER ON THE FAR SIDE OF
	// " THE BOUNDARY. IF IR(NP) IS 1 THEN THE PARTICLE HAS ESCAPED THE REGION
	// " OF INTEREST AND THE HISTORY IS TERMINATED.(IDISC IS SET TO 1.)
	// "
	// "
	// "
	// "
	// " SOME VARIABLES
	// " ==============
	// ;"
	// "OUTEND = .TRUE. => PARTICLE MAY TRANSMIT OR BACKSCATTER OUT ENDS
	// " = .FALSE. => PARTICLE STAYS WITHIN THE END BOUNDARIES
	// "OUTSID = .TRUE. => PARTICLE MAY TRANSMIT OUT THE SIDES
	// " = .FALSE. => PARTICLE STAYS WITHIN THE SIDE BOUNDARY
	// "IRL = STARTING REGION NUMBER THE PARTICLE IS IN
	// "IZ = STARTING PLANAR ZONE NUMBER THE PARTICLE IS IN.
	// " THE PARTICLE IS BETWEEN ZPLANE(IZ) AND ZPLANE(IZ+1).
	// "IX = STARTING CYLINDRICAL ZONE NUMBER THE PARTICLE IS IN.
	// " THE PARTICLE IS BETWEEN RCYL(IX-1) AND RCYL(IX).
	// "
	// " COMMON/GEOM/
	// " ZPLANE(IZ) Z VALUES OF PLANES
	// " 1<=IZ<=NZ+1
	// " RCYL(IRR) RADII OF CYLINDERS
	// " 1<=IRR<=NR
	// " CYRAD2(IRR) =RCYL(IRR)**2
	// " NZ # PLANAR GEOMETRICAL ZONES (NPLANE-1)
	// " ZONE(I) IS BETWEEN ZPLANE(I) AND ZPLANE(I+1)
	// " NR # CYLINDRICAL GEOMETRICAL ZONES
	// " ZONE(I) IS BETWEEN RCYL(I-1) AND RCYL(I)
	// " NREG TOTAL # GEOMETRICAL ZONES =NR*NZ +1
	// " +1 FOR VACUUM ENVELOPE
	// %E "cavrznrc.mortran"
	// " DEFINITIONS OF REGION NUMBER, PLANAR ZONE, CYLINDRICAL ZONE
	// " ===========================================================
	// " Z AXIS RUNS ACROSS PAGE SHOWN AS .......
	// "
	// "
	// " 1
	// /"
	// "
	// " --------------------------------------------------------- RCYL(NR)
	// " |(NR-1) |(NR-1) |(NR-1) | . . . . | NR*NZ | NR*NZ | IX=NR
	// " | *NZ+2 | *NZ+3 | *NZ+4 | | | +1 |
	// " --------------------------------------------------------- RCYL(NR-1)
	// " | . | . | . | | . | . |
	// " | . | . | . | | . | . |
	// " | . | . | . | | . | . |
	// " --------------------------------------------------------- RCYL(2)
	// " | NZ+2 | NZ+3 | NZ+4 | . . . . | 2NZ | 2NZ+1 | IX=2
	// " --------------------------------------------------------- RCYL(1)
	// "..1....|...2...|...3...|...4...|...............|...NZ..|..NZ+1.|....IX=1..1..
	// ;" ---------------------------------------------------------
	// " | | | | . . . . | | |
	// " ---------------------------------------------------------
	// " | . | . | . | | . | . |
	// " | . | . | . | | . | . |
	// " | . | . | . | | . | . |
	// " ---------------------------------------------------------
	// " | | | | . . . . | | |
	// " | | | | | | |
	// " ---------------------------------------------------------
	// " IZ=1 IZ=2 IZ=3 IZ=NZ-1 IZ=NZ
	// "
	// " 1
	// "
	// "
	// "
	// " VERSION 1 ADAPTED FROM CAVITY HOWFAR 06/84 ERIC FOX
	// " VERSION 2 THE SUBROUTINE CALLS TO PLANES AND 10/87 AFB
	// " CYLINDER HAVE BEEN REPLACED BY MACROS
	// " TO SPEED THINGS UP
	// "
	// "
	// "******************************************************************************
	/**
	 * The following is a general specification of HOWFAR: 
	 * Given a particle at (X,Y,Z) in region IR and going in direction 
	 * (U,V,W), this routine answers the question, can the particle go 
	 * a distance USTEP without crossing a boundary? If yes, it merely returns; 
	 * If no, it sets USTEP=distance to boundary in the current 
	 * direction and sets IRNEW to the region number on the far side of the boundary. 
	 * Interface method.
	 */
	public void HOWFAR() {
		// $IMPLICIT-NONE;

		// "MACRO USED LOCALLY TO CHANGE REGIONS, ADJUST USTEP, AND EXIT"
		// REPLACE {$SET NEW REGION(#,#);} WITH
		// {
		// IF({P1}.LE.USTEP)[USTEP={P1};IRNEW={P2};]RETURN;
		// }

		// "Debug macro, usually not used, to verify that region is correct on entry"
		// "    It must be used after IX and IZ have been determined"
		// "    REAL Rdebug must be declared"
		// REPLACE {$Check_region;} WITH
		// {
		// Rdebug = SQRT(X(NP)**2 + Y(NP)**2);
		// IF(IX = 1)["we are in inner radial region"
		// IF(Rdebug > (RCYL(IX)+1.E-6) | (Z(NP) < ZPLANE(IZ) & Z(NP) >
		// ZPLANE(IZ+1)))[
		// OUTPUT IRL, IX, IZ, X(NP), Y(NP), Z(NP), Rdebug, RCYL(IX),
		// ZPLANE(IZ),ZPLANE(IZ+1);
		// (/' Error in HOWFAR, particle in wrong region'/
		// ' IRL, IX,IZ =', 3I7/
		// ' X,Y,Z,R=', 4F15.7/
		// ' RCYL(1)=', F15.7/
		// ' ZPLANE(IZ),ZPLANE(IZ+1)=', 2F15.7);
		// ]
		// ]"end of IX=1 block"
		// ELSE [
		// IF((Rdebug > (RCYL(IX)+1E-6) & Rdebug < (RCYL(IX-1) -1.E-6)) |
		// (Z(NP) < ZPLANE(IZ) & Z(NP) > ZPLANE(IZ+1)))[
		// OUTPUT IRL, IX, IZ, X(NP), Y(NP), Z(NP), Rdebug, RCYL(IX-1),
		// RCYL(IX),
		// ZPLANE(IZ),ZPLANE(IZ+1);
		// (/' Error in HOWFAR, particle in wrong region'/
		// ' IRL, IX,IZ =', 3I7/
		// ' X,Y,Z,R=', 4F15.7/
		// ' RCYL(IX-1), RCYL(IX)=', 2F15.7/
		// ' ZPLANE(IZ),ZPLANE(IZ+1)=', 2F15.7);
		// ]
		// ]
		// ;}

		// LOGICAL OUTEND,OUTSID;

		// $INTEGER IRL,IX,IZ,IHITP,IHITC,IZNEW,IXNEW;
		// $REAL WL,TPLANE,U1,V1,A,TCYL,X1,Y1,B,B2,C,COUT,CIN,RAD;
		// "$REAL Rdebug; " "Uncomment if using $Check_region;"

		// "First set idisc and irnew "
		EGS4.IDISC = 0;
		EGS4.IRNEW = EGS4.IR[EGS4.NP - 1];

		// "DISCARD ZERO WEIGHT PARTICLES"
		if (EGS4.WT[EGS4.NP - 1] == 0.0) {
			EGS4.IDISC = 1;
			return;
		}

		// "INITIALLY ASSUME PARTICLE STAYS IN THE TARGET"
		//boolean OUTEND = false;
		//boolean OUTSID = false;

		int IRL = EGS4.IR[EGS4.NP - 1];// "LOCAL REGION NUMBER"

		// "DISCARD IF PARTICLE WANTS TO LEAVE THE GEOMETRY"
		if (IRL == 1) {
			EGS4.IDISC = 1;
			return;
		}

		// $GET-IX-IZ(IRL); //"GET PLANAR AND CYLINDRICAL ZONES NUMBERS"
		int IX = EGS4Geom.GET_IX(IRL);
		int IZ = EGS4Geom.GET_IZC(IRL);

		// "Following commented out usually"
		// "$Check_region;"
		// "write(6,*);"
		// "write(6,*) 'howfar: ',ix,iz,ustep,nz,nr;"

		// $PLANES(IZ,IZ+1,IHITP,TPLANE,ustep);"GET DISTANCE TO PLANE"
		// "IHITP  =  1 => HITS GREATER Z PLANE"
		// "       =  0 => MISSES BOTH PLANES"
		// "       = -1 => HITS LESSER Z PLANE"
		EGS4Geom.ustep = EGS4.USTEP;
		EGS4Geom.PLANES(IZ, IZ + 1);
		EGS4.USTEP = EGS4Geom.ustep;
		// $CYLNDR(IX,IHITC,TCYL,ustep);"GET DISTANCE TO CYLINDER"
		// "       IHITC   =  1 => HITS OUTER CYLINDER"
		// "               =  0 => MISSES BOTH CYLINDERS"
		// "               = -1 => HITS INNER CYLINDER"
		EGS4Geom.ustep = EGS4.USTEP;
		EGS4Geom.CYLNDR(IX);
		EGS4.USTEP = EGS4Geom.ustep;

		return;
	}// "END OF SUBROUTINE HOWFAR"

	/**
	 * Gather all media data required for this simulation.
	 */
	private void HATCH() {
		Calendar cal = Calendar.getInstance();
		Date d = cal.getTime();
		EGS4.seqStr = "Start of run:         " + d.toString();
		if (EGS4.iprint > 0)
			printSequence(EGS4.seqStr);

		EGS4.HATCH();
	}

	/**
	 * Start the shower, i.e. the actual simulation for electron-photon transport
	 */
	private void SHOWER() {
		// //CALL SHOWER(IQIN,EI,XIN,YIN,ZIN,UIN,VIN,WIN,IRIN,WEIGHT);
		EGS4Core.SHOWER(
				EGS4SrcEns.iqin,
				EI,// EGS4SrcEns.ein,
				EGS4SrcEns.xin, EGS4SrcEns.yin, EGS4SrcEns.zin, EGS4SrcEns.uin,
				EGS4SrcEns.vin, EGS4SrcEns.win, EGS4SrcEns.irin,
				EGS4SrcEns.WEIGHT);
	}

	/**
	 * Setup input variables.
	 */
	private void inputs() {
		// @ TITLE
		TITLEs = "cavrznrc_template: 1.25 MeV on graphite pancake chamber";
		// @ IWATCH
		// #off,interactions,steps,deposited,graph;
		// #debug output with increasing detail, graph outputs .gph file for
		// EGS_Windows
		// #if not "off" use very few histories
		// IWATCH=IWATCH_OFF;
		// @ STORE INITIAL RANDOM NUMBERS
		// #no,last,all deposited,all;
		// #last: store initial random numbers for last history in .egsrns
		// #all deposited: store initial random numbers for each history that
		// deposits energy
		// in the cavity
		// #all: store initial random numbers for each history
		// ISTORE=ISTORE_NO;
		// @ IRESTART
		// #first,restart,make,analyze,for graphics,parallel;
		// #first: first run
		// #restart: restart of old run (requires .egsdat file)
		// #make: just create an input file and exit
		// #analyze: read in data from .egsdat file and do statistical analysis
		// and output results
		// # to .egslst
		// #for graphics: read starting random numbers from .egsrns--eg for
		// output to graphics package
		// #parallel: read .egsdat files from parallel jobs (named
		// inputfile_w*), do statistical
		// # analysis and output to .egslst
		// IRESTART=IRESTART_FIRST;
		// @ OUTPUT OPTIONS
		// #short,cavity details;
		// #short: output cavity summary + dose grid
		// #cavity details: above plus details for every zone in cavity
		// IOOPTN=IOOPTN_SHORT;
		// @ STORE DATA ARRAYS
		// #yes,no;
		// #yes: output .egsdat file for restarts, parallel post-processing, etc
		// IDAT=IDAT_YES;
		// @ NUMBER OF HISTORIES
		// #splits into $STAT statistical batches
		// #must be >=$STAT**2 if IWATCH= Off
		// #can have less than this if IWATCH set to another option
		// NCASE=10000;
		// @ MAX CPU HOURS ALLOWED
		// #Will shut down cleanly prior to exceeding this limit, as long as one
		// #batch has completed.
		// TIMMAX=90.000;
		// @ IFULL
		// #dose and stoppers,Aatt and Ascat,Ap,Afl and <s>g/w;
		// #dose and stoppers: output total dose plus that due to stoppers and
		// discards
		// #Aatt and Ascat: above plus Aatt, Ascat
		// #Ap: above plus Ap
		// #Afl and <s>g/w: above plus Afl and stopping power ratio gas/water
		// EGS4Macro.IFULL=0;//IFULL_AATT_AND_ASCAT;
		// @ STATISTICAL ACCURACY SOUGHT
		// #If 0, goes until number of histories or CPU limit exceeded.
		// #If not zero goes until this uncertainty (in %) is achieved in the
		// peak dose region
		// STATLM=0.0;
		// @ PHOTON REGENERATION
		// #no,yes,no electrons from wall;
		// #no: normal calculation
		// #yes: regenerate parent photon after interaction (used for FANO
		// calculations)
		// #no electrons from wall: photons not regenerated,
		// # secondary electrons from cavity wall are eliminated
		// IFANO=IFANO_NO;
		// @ INITIAL RANDOM NO. SEEDS
		// #With ranmar: these must be between 1 and 30081 (default to 9373)
		// #With ranlux: 1st is luxury level 0->4 allowed but should not be 0
		// # 2nd is seed 1 -> 1073741824
		// jrng1=1;jrng2=3;
		// right here, initialize random generator!!!
		// init_random_generator();
		// ###########################################################################################
		// @ METHOD OF INPUT
		// #groups,individual,cavity information:
		// #group: input groups of slabs of equal thickness
		// #individual: input Z of bottom of every slab
		// #cavity information: generate simple geometry from cavity info input
		// in section below.
		// # If you use this, there are no more inputs in this section
		// EGS4Geom.iterseindex=EGS4Geom.iINDIVIDUAL;
		if (EGS4Geom.iterseindex == EGS4Geom.iCAVITY_INFORMATION)// 2
		{
			EGS4Geom.WALLTH = 0.5;// WALL THICKNESS
			EGS4Geom.CAVRAD = 1.0;// CAVITY OUTER RADIUS
			EGS4Geom.CAVLNG = 2.0;// CAVITY LENGTH
			EGS4Geom.ELERAD = 0.01;// ELECTRODE RADIUS
			EGS4Geom.SLENGHT = "H2O_fortran";// WALL MATERIAL
			EGS4Geom.airs = "AIR521ICRU_fortran";// ELERAD=0.0
			EGS4Geom.electrods = "AL521ICRU_fortran";// ELERAD!=0.0;
		} else// //ITERSE->0 sau 1
		{
			// EGS4Geom.Z_OF_FRONT_FACE=0.0;//#Beginning of first slab
			if (EGS4Geom.iterseindex == EGS4Geom.iGROUPS)// 0
			{
				// NSLAB
				EGS4Geom.nNSLAB = 2;// 2 value
				// #Define a group of 10 slabs with thickness 1 cm
				// #followed by 10 slabs with thickness 2 cm
				EGS4Geom.NSLAB[0] = 10;
				EGS4Geom.NSLAB[1] = 10;
				// SLAB THICKNESS
				EGS4Geom.DELTAZ[0] = 1.0;
				EGS4Geom.DELTAZ[1] = 2.0;
			}

			if (EGS4Geom.iterseindex == EGS4Geom.iINDIVIDUAL)// 1
			{
				// NSLAB
				// EGS4Geom.nNSLAB=3;
				// DEPTH BOUNDARIES
				// EGS4Geom.ZPLANE[1]=0.05;//0.3;
				// EGS4Geom.ZPLANE[2]=6.3;//0.5;//+0.95
				// EGS4Geom.ZPLANE[3]=7.25;//0.8;
			}

			// EGS4Geom.nCyl=2;//"number of radial cylinders input"
			// #Radii of cylinders
			// EGS4Geom.RCYL[1]=6.3/2;//1.0;
			// EGS4Geom.RCYL[2]=7.3/2;//1.3;

			// MEDIA=the media in the problem. These must match exactly,
			// including case, one
			// of the media names in the pegs4 data set being used in the
			// problem.
			// #Next we specify which media are in which geometric regions
			// #note that by default all regions contain
			// #medium 1 and which medium to input as 1 should be selected with
			// this in mind.
			// EGS4Geom.nMEDIA=2;
			// EGS4.MEDIA[0]="AL521ICRU_fortran";//"170C521ICRU_fortran";
			// EGS4.MEDIA[1]="NAI_Fortran";//"AIR521ICRU_fortran";

			// DESCRIPTION BY:#planes,regions;
			// #planes: use slab and cylinder no.'s to define what medium goes
			// where
			// #regions: use region numbers to define this (region numbers start
			// at 2 and
			// number from top to bottom of geometry and innermost radius to
			// outermost radius)
			// EGS4Geom.DESCRIBE=EGS4Geom.DESCRIBE_REGIONS;//##INPUT DATA

			// EGS4Geom.nMEDNUM=1;
			// EGS4Geom.MEDNUM[0]=2;//MEDNUM->AIR

			if (EGS4Geom.DESCRIBE == EGS4Geom.DESCRIBE_REGIONS_DENSITY
					|| EGS4Geom.DESCRIBE == EGS4Geom.DESCRIBE_PLANES_DENSITY) {
				EGS4Geom.nRHOR = 1;
				EGS4Geom.RHOR[0] = 0.;
			}
			if (EGS4Geom.DESCRIBE == EGS4Geom.DESCRIBE_REGIONS
					|| EGS4Geom.DESCRIBE == EGS4Geom.DESCRIBE_REGIONS_DENSITY) {
				// EGS4Geom.nNREGLO=1;
				// EGS4Geom.nNREGHI=1;
				// EGS4Geom.NREGLO[0]=3;//START REGION
				// EGS4Geom.NREGHI[0]=3;//STOP REGION
			}
			if (EGS4Geom.DESCRIBE == EGS4Geom.DESCRIBE_PLANES
					|| EGS4Geom.DESCRIBE == EGS4Geom.DESCRIBE_PLANES_DENSITY) {
				EGS4Geom.nNZLO = 2;
				EGS4Geom.nNZHI = 2;
				EGS4Geom.nNRLO = 2;
				EGS4Geom.nNRHI = 2;
				// #This puts MEDIUM 1 everywhere and then #inserts a small
				// column of MEDIUM 2
				// on the central #axis with radius 1cm and going from Z=10cm
				// #to Z=12cm
				EGS4Geom.NZLO[0] = 1;// START ZSLAB
				EGS4Geom.NZLO[1] = 10;// START ZSLAB
				EGS4Geom.NZHI[0] = 20;// STOP ZSLAB
				EGS4Geom.NZHI[1] = 12;// STOP ZSLAB
				EGS4Geom.NRLO[0] = 1;// START RING
				EGS4Geom.NRLO[1] = 1;// START RING
				EGS4Geom.NRHI[0] = 2;// STOP RING
				EGS4Geom.NRHI[1] = 1;// STOP RING
			}
		}// ITERSE->0 sau 1
			// cavity inputs:
			// EGS4Geom.NSUMCV=1;//NUMBER OF CAVITY REGIONS
		// #this defines the small cylinder of
		// #air in the geometry to be the cavity
		// EGS4Geom.ISUMCV[0]=3;//REGION NUMBERS OF THE CAVITY= 3

		// CALL GEOMRZ!!
		EGS4.seqStr = " *** Reading geometrical info ... ***";
		if (EGS4.iprint > 1)
			printSequence(EGS4.seqStr);
		EGS4.iprint = 0;
		EGS4Geom.GEOMRZ();
		if (EGS4.STOPPROGRAM) {
			closeFile();
			return;
		}
		EGS4.iprint = IPRINT;
		if (SOURCE == SOURCE_POINT) {
			EGS4.seqStr = " Point source on symmetry z-axis ";
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);
			EGS4.seqStr = " Distance to entry front face (cm) = "
					+ EGS4.format(EGS4SrcEns.source_option[0], 8, true);
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);
			EGS4.seqStr = " Radius of hitting surface (cm) = "
					+ EGS4.format(EGS4Geom.RCYL[EGS4Geom.NR], 8, true);
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);
		} else if (SOURCE == SOURCE_SARPAGAN) {
			EGS4.seqStr = " Volumic (Sarpagan) source on top of the detector house ";
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);
			EGS4.seqStr = " Source height (cm)= "
					+ EGS4.format(hsource, 8, true);
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);
			EGS4.seqStr = " Source diameter (cm)= "
					+ EGS4.format(2 * asource, 8, true);
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);

		} else if (SOURCE == SOURCE_MARINELLI) {
			EGS4.seqStr = " Volumic (Marinelli) source surrounding the detector house ";
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);
			EGS4.seqStr = " Source full height (cm)= "
					+ EGS4.format(hsource, 8, true);
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);
			EGS4.seqStr = " Source height above the detector top plane(cm)= "
					+ EGS4.format(hsourceup, 8, true);
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);
			EGS4.seqStr = " Source external diameter (cm)= "
					+ EGS4.format(2 * asource, 8, true);
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);
			EGS4.seqStr = " Source internal diameter (cm)= "
					+ EGS4.format(2 * bsource, 8, true);
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);

		} else if (SOURCE == SOURCE_BEAM) {
			EGS4.seqStr = " Beam source at an angle ";
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);
			EGS4.seqStr = " Distance to entry front face (cm) = "
					+ EGS4.format(EGS4SrcEns.source_option[0], 8, true);
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);
			EGS4.seqStr = " Beam radius (cm) = "
					+ EGS4.format(beam_radius, 8, true);
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);
			EGS4.seqStr = " Beam angle (degrees) = "
					+ EGS4.format(beam_angle, 8, true);
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);

		}

		// INCIDENT PARTICLE= photon #electron,photon,positron,all,charged;
		// #all & charged: only for phase space sources
		// ###########################################################################################
		// EGS4SrcEns.ipart=EGS4SrcEns.ipart_photon;
		// SOURCE NUMBER= 1 #0,1,2,3,4,10,11,12,13,14,15,16,20,21,22,23
		// EGS4SrcEns.ISOURC=EGS4SrcEns.point_source_on_axis_incident_from_the_front;//1;
		// SOURCE OPTIONS= 100., 1.3, 0, 0
		// #for source 1: SSD of beam, radius of # beam on front surface
		if (EGS4SrcEns.ISOURC != EGS4SrcEns.parallel_beam_incident_from_the_front_with_radial_distribution)// 20)
		{
			// EGS4SrcEns.source_option[0]=1.5;//100.;
			// EGS4SrcEns.source_option[1]=7.3/2;//1.3;
			// EGS4SrcEns.source_option[2]=0.;
			// EGS4SrcEns.source_option[3]=0.;

			if (EGS4SrcEns.ISOURC == EGS4SrcEns.full_phase_space_from_any_angle)// 22)
			{
				EGS4SrcEns.source_option[4] = 0.;
				EGS4SrcEns.source_option[5] = 0.;
				EGS4SrcEns.source_option[6] = 0.;
				EGS4SrcEns.source_option[7] = 0.;
				EGS4SrcEns.source_option[8] = 0.;
			} else if (EGS4SrcEns.ISOURC == EGS4SrcEns.BEAM_treatment_head_simulation_from_any_angle)// 23
			{
				EGS4SrcEns.source_option[4] = 0.;
			}
		}
		if (EGS4SrcEns.ISOURC == EGS4SrcEns.parallel_beam_incident_from_the_front_with_radial_distribution)// 20)
		{
			// source 20 Ex:
			EGS4SrcEns.MODEIN = EGS4SrcEns.MODEIN_LOCAL;
			EGS4SrcEns.NRDIST = 1;
			EGS4SrcEns.RDISTF[0] = 1.0;
			EGS4SrcEns.RPDF[0] = 1.0;
			EGS4SrcEns.RDIST_IOUTSP = EGS4SrcEns.RDIST_IOUTSP_NONE;
		}
		// CALL SRCRZ->"Get source data"
		if (SOURCE == -1) {
			EGS4.seqStr = " *** Reading radiation source info ... ***";
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);
		}
		if (SOURCE == -1)
			EGS4SrcEns.SRCRZ();
		if (EGS4.STOPPROGRAM) {
			closeFile();
			return;
		}
		// source initialization
		if (SOURCE == -1)
			EGS4SrcEns.SRCINI();
		if (EGS4.STOPPROGRAM) {
			closeFile();
			return;
		}

		if (EGS4SrcEns.ISOURC == 21 || EGS4SrcEns.ISOURC == 22
				|| EGS4SrcEns.ISOURC == 23) {
			EGS4SrcEns.MONOEN = 0;
		}
		// "no need to input monoen for source 21,22"
		else {
			// INCIDENT ENERGY= monoenergetic #monoenergetic, spectrum;
			// EGS4SrcEns.monoindex=EGS4SrcEns.iMONOENERGETIC;
			// INCIDENT KINETIC ENERGY(MEV)= 1.25 #only use for "monoenergetic"
			// EGS4SrcEns.ikemev=0.662;//1.25;

			EGS4SrcEns.ENSRC();
			if (EGS4.STOPPROGRAM) {
				closeFile();
				return;
			}
		}// "Get data re-source energies"
			// ###########################################################################################
		// get_transport_parameter();
		setEcutRegion = false;
		setPcutRegion = false;
		setSmaxRegion = false;
		ecut = 0.521;
		pcut = 0.004;
		smax = 1.e10;

		if (!setEcutRegion)
			EGS4.ECUT[0] = ecut;// #Electron cutoff for transport
		if (!setPcutRegion)
			EGS4.PCUT[0] = pcut;// #Photon cutoff for transport
		if (!setSmaxRegion)
			EGS4.SMAXIR[0] = smax;// #Maximum step size in cm (not needed
									// #unless old PRESTA algorithm used)
		if (setEcutRegion) {
			nEcut = 1;// number of data
			Ecut[0] = 1.;
			startEcutRegion[0] = 1;
			stopEcutRegion[0] = 1;
		}
		if (setPcutRegion) {
			nPcut = 1;// number of data
			Pcut[0] = 1.;
			startPcutRegion[0] = 1;
			stopPcutRegion[0] = 1;
		}
		if (setSmaxRegion) {
			nSmax = 1;// number of data
			Smax[0] = 1.;
			startSmaxRegion[0] = 1;
			stopSmaxRegion[0] = 1;
		}
		// Bound Compton scattering= On or Off->IBCMP
		// #######//On means 1 and Off means 0
		// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		// #Off: Klein-Nishina used for compton
		// # scattering
		// #On: Impulse approximation used for
		// # scattering
		// #It has not been established that the
		// #the scoring routines work with this
		// #option ON.
		// incoh=incoh_ON;//incoh_ON;//incoh_OFF;
		if (incoh == incoh_OFF || incoh == incoh_ON)
			setIncohRegion = false;
		else
			setIncohRegion = true;
		if (setIncohRegion) {
			nIncoh = 1;// number of data
			Incoh[0] = incoh_OFF;// 0 or 1
			startIncohRegion[0] = 1;
			stopIncohRegion[0] = 1;
		}
		// Rayleigh scattering= On->IRAYLR
		// #Off: no coherent scattering
		// #On: simulates coherent scattering
		// coh=coh_ON;//coh_ON;//coh_OFF;
		if (coh == coh_OFF || coh == coh_ON)
			setCohRegion = false;
		else
			setCohRegion = true;
		if (setCohRegion) {
			nCoh = 1;// number of data
			Coh[0] = coh_OFF;// 0 or 1
			startCohRegion[0] = 1;
			stopCohRegion[0] = 1;
		}
		// Atomic relaxations= On->IEDGFL
		// #Note, it has not been verified that
		// #AUSGAB works with this option ON
		// #On: use correct cross section
		// # for p.e. events and shell vacancies
		// # for Compton & p.e. events are relaxed
		// # via emission of fluorescent X-Rays,
		// # Auger and Koster-Cronig electrons
		// relax=relax_ON;//relax_ON;//relax_OFF;
		if (relax == relax_OFF || relax == relax_ON)
			setRelaxRegion = false;
		else
			setRelaxRegion = true;
		if (setRelaxRegion) {
			nRelax = 1;// number of data
			Relax[0] = relax_OFF;// 0 or 1
			startRelaxRegion[0] = 1;
			stopRelaxRegion[0] = 1;
		}
		// Photoelectron angular sampling= On->IPHTER
		// #Off: Photoelectrons get direction of
		// # photon that creates them
		// #On: Sauter's formula is used
		// pe=pe_ang_ON;//pe_ang_ON;//pe_ang_OFF;
		if (pe == pe_ang_OFF || pe == pe_ang_ON)
			setPeRegion = false;
		else
			setPeRegion = true;
		if (setPeRegion) {
			nPe = 1;// number of data
			Pe[0] = pe_ang_OFF;// 0 or 1
			startPeRegion[0] = 1;
			stopPeRegion[0] = 1;
		}
		// Brems angular sampling= On ->IBRDST
		// #Simple: leading term of Koch-Motz
		// # dist'n used to determine angle
		// # of bremsstrahlung photons
		// #KM: Koch-Motz distribution used to
		// # determine angle -> 2BS(modified)
		// EGS4.ibrdst = brems_ang_KM;
		// Pair angular sampling= On->IPRDST
		// #Simple: use leading term of K-M
		// # dist'n
		// #KM: use complete Koch and Motz dist'n
		// #Off: angle of pairs is m/E--like old EGS4
		// EGS4.iprdst = pair_ang_SIMPLE;
		// ESTEPE= 0.25->ESTEPE
		// #Max fractional continuous energy loss
		// #per step. Use 0.25 unless using
		// #PRESTA-I
		EGS4.estepe = 0.25;
		// XIMAX= 0.5->XIMAX
		// #Max first elastic scattering moment
		// #per step. Using default.
		EGS4.ximax = 0.5;
		// Boundary crossing algorithm= exact->bca_algorithm, exact_bca
		// #exact: cross boundaries in single scattering
		// # mode (distance at which to go into
		// # single scattering mode determined by
		// # "Skin depth for BCA"
		// #PRESTA-I: cross boundaries with lateral
		// # correlations off and force multiple
		// # scattering mode
		// EGS4.bca_algorithm = BCA_EXACT;//exact means
		// =EGS4.$BCA_ALGORITHM_DEFAULT=0;
		// Skin depth for BCA= 3->skindepth_for_bca
		// #Distance from a boundary (in elastic
		// #MFP) at which the algorithm will go
		// #into single scattering mode (using
		// #default here)
		EGS4.skindepth_for_bca = 3.0;
		// Electron-step algorithm= default->transport_algorithm
		// PRESTA-II (the default),$PRESTA_II = 0;
		// #Determines the algorithm used to take
		// #into account lateral and longitudinal
		// #correlations in a condensed history
		// #step
		// EGS4.transport_algorithm = estep_alg_PRESTA_II;//0;
		// Spin effects= On->spin_effects
		// #Turns off/on spin effects for electron
		// #elastic scattering. Spin On is
		// #ABSOLUTELY necessary for good
		// #backscattering calculations. Will
		// #make a difference even in `well
		// #conditioned' situations (e.g. depth
		// #dose curves).
		// ispin = spin_ON;
		// Brems cross sections= BH
		// #BH: Bethe-Heitler cross-sections used
		// #NIST: NIST cross-sections used
		// EGS4.ibr_nist=brems_cross_BH;
		// " Pair cross sections "= BH
		// EGS4.pair_nrc = pair_cross_BH;
		// Electron impact ionization
		// EGS4.eii_flag=eii_ON;
		// Triplet
		// EGS4.itriplet=triplet_ON;
		// Radiative compton correction
		// EGS4.radc_flag=radc_ON;
		// /###########################################################################################
		// VARIANCE REDUCTION:
		// ELECTRON RANGE REJECTION
		// EGS4Macro.irejct=irejct_OFF;//irejct_OFF;//irejct_ON;
		// #On: if charged particle energy is below ESAVEIN
		// # and it cannot get out of current region
		// # with energy > ECUT, the particle is
		// # terminated
		// #also terminates all electrons which cannot
		// #reach the cavity under conservative assumptions.
		// ESAVEIN=2.0;//#total energy below which range rejection is considered
		// EGS4Macro.cs_enhance=1.0;// #Photon cross section scaling factors
		// RUSSIAN ROULETTE DEPTH= 0.0000 #play Russian Roulette with photons
		// once they
		// #cross this Z plane
		// RRZ=0.000;
		// RUSSIAN ROULETTE FRACTION= 0.0000 #probability of photon survival--if
		// this
		// #and #RUSSIAN ROULETTE DEPTH both 0, then
		// #photon Russian Roulette is not played
		// RRCUT=0.000;
		// #exponential pathlength biasing can be
		// #used. See Rogers&Bielajew 1990 review for
		// #discussion. C<0 => pathlength shortening
		// # >0 => pathlength stretching
		// # along z axis both cases
		// #CAVRZnrc allows for having the photon cross
		// #section scaled to enhance interactions.
		// # If this input is missing or set to <= 1, it
		// # has no effect on the simulation. But if
		// # the enhancement factor is set to > 1, the
		// # effect is dramatic: all other user input
		// # concerning photon forcing, splitting, exp.
		// # transform, etc., is ignored. In addition,
		// # the calculation result corresponds ALWAYS
		// # to 'Aatt and Ascat', no matter what the
		// # user requested (but only Awall is calculated,
		// # not the individual Ascat and Aatt).
		// # The algorithm employed is implemented via
		// # $RAYLEIGH-CORRECTION and appropriate calls to
		// # AUSGAB. For more detail see the manual and the
		// # header of cavrznrc.mortran
		// EGS4Macro.CEXPTR=0.000;
		// PHOTON FORCING= On #Off (default),On;
		// #On: force photons to interact according to
		// # START FORCING and STOP FORCING AFTER inputs
		// IFARCE=IFARCE_OFF;//IFARCE_ON;
		// EGS4Macro.NFMIN=1;//#Start forcing at this interaction number
		// EGS4Macro.NFMAX=1;//#Number of photon interactions after which
		// #to stop forcing photon interactions
		// PHOTON SPLITTING= 1 #no. of times to split a photon
		// #if < 2-->normal transport
		// #overrides PHOTON FORCING if >= 2
		// #can only be >= 2 if IFULL= dose and stoppers
		// # or if IFULL= Aatt and Ascat
		// phsplitt=1;

		test_inputs();
		if (EGS4.STOPPROGRAM) {
			closeFile();
			return;
		}
		// ##########################################################################################
		// ##########################################################################################
		// ##########################################################################################
		// ##########################################################################################
		// ##########################################################################################
		// " SCORING ARRAY INITIALISATION
		// " ****************************
		NCASEO = 0;
		NCASET = 0;
		TMCPUO = 0;
		NNREADO = 0; // "SET PREVIOUS RUN COUNTERS"
		if (IRESTART == 0 || IRESTART == 5) {// "FRESH START, SET EVERYTHING TO ZERO"
			EGS4SrcEns.NNREAD = 0;
			SCSTP = 0.;
			SCSTP2 = 0.;
			SCSTP_TMP = 0.;
			SCSTP_LAST = 0;
			SCCSTP = 0.;
			SCCSTP2 = 0.;
			SCCSTP_TMP = 0.;
			SCCSTP_LAST = 0;
			cav_dose = 0;
			cav2_dose = 0;
			cav_dose0 = 0;
			cav2_dose0 = 0;
			cav_dose1 = 0;
			cav2_dose1 = 0;
			cav_dose2 = 0;
			cav2_dose2 = 0;

			cav_dosec = 0;
			cav_dosec01 = 0;
			cav_dosec02 = 0;

			tmp_dose = 0.0;
			tmp_dose0 = 0.0;
			tmp_dose1 = 0.0;
			tmp_dose2 = 0.0;

			PIISTP = 0.;
			for (int IZ = 1; IZ <= EGS4Geom.NZ; IZ++) {
				for (int IX = 1; IX <= EGS4Geom.NR; IX++)// DO IX=1,NR[
				{
					for (int IT = 1; IT <= 4; IT++) {
						SCDOSE[IZ - 1][IX - 1][IT - 1] = 0.;
						SCDOSE2[IZ - 1][IX - 1][IT - 1] = 0.;
						SCDOSE_TMP[IZ - 1][IX - 1][IT - 1] = 0.;
						if (IT < 4)
							SCDOSE_COV[IZ - 1][IX - 1][IT - 1] = 0.;
					}
					SCDOSE_LAST[IZ - 1][IX - 1] = 0;
				}
			}
		} else if (IRESTART != 4)// NOT ALLOWED HERE
		{
			// "Restart or stats analysis only, read old data from unit 4"
			// "open unit 4 as an old file"
			// OUTPUT;(' About to read the previous .egsdat file');
		}

		if (IRESTART == 3) {
			NCASE = 0;
		}

		NCASET = NCASE + NCASEO;
		EGS4SrcEns.NCASET = NCASET;

		EGS4.seqStr = " ********* SUCCESSFUL INPUT ACCOMPLISHED *********";
		if (EGS4.iprint > 1)
			printSequence(EGS4.seqStr);

	}// inputs

	/**
	 * Validate more inputs. Called by inputs routine.
	 */
	private void test_inputs() {
		if (IWATCH < 0 || IWATCH > 4) {
			IWATCH = 0;
		}// default
		EGS4Macro.IWATCH = IWATCH;
		EGS4Geom.IWATCH = IWATCH;
		// if(ISTORE<0 || ISTORE>3){ISTORE=0;}//default
		if (ISTORE < 0 || ISTORE > 0) {
			ISTORE = 0;
		}// default
		// if(IRESTART<0 || IRESTART>4){IRESTART=0;}//default
		if (IRESTART < 0 || IRESTART > 0) {
			IRESTART = 0;
		}// default
		if (IOOPTN < 0 || IOOPTN > 1) {
			IOOPTN = 0;
		}// default
		if (IDAT < 0 || IDAT > 1) {
			IDAT = 1;
		}// default
		if (IDAT < 1 || IDAT > 1) {
			IDAT = 1;
		}// default->NO FILE STORAGE
		if (IRESTART == 4)// no parrallelllllllllllll!!
		{
			IDAT = 1; // "do not store output in this case to avoid biasing"
			ISTORE = 0; // "do not store the starting random numbers either"
		}
		if (NCASE < NCASE_MIN || NCASE > NCASE_MAX) {
			NCASE = NCASE_DEFAULT;
		}
		if (TIMMAX < TIMMAX_MIN || TIMMAX > TIMMAX_MAX) {
			TIMMAX = TIMMAX_DEFAULT;
		}
		if (EGS4Macro.IFULL < 0 || EGS4Macro.IFULL > 1) {
			EGS4Macro.IFULL = 0;
		}// default
		if (STATLM < STATLM_MIN || STATLM > STATLM_MAX) {
			STATLM = STATLM_DEFAULT;
		}
		if (IFANO < 0 || IFANO > 2) {
			IFANO = 0;
		}// default
		if (IWATCH == 0 && NCASE < $NCASEMIN) {
			NCASE = $NCASEMIN;
		}
		if (IFANO == 1) {// "With ifano option turned on, it is a waste of time to"
							// "have RAYLEIGH turned on (scattered photon will be killed,"
							// "original photon re-created) => turn Rayleigh off"
			for (int j = 1; j <= EGS4.$MXREG; j++) {
				EGS4.IRAYLR[j - 1] = 0;
			}
			EGS4.seqStr = " ******** ifano set => turning off Rayleigh! **** ";
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);

			// OUTPUT; (//' ******** ifano set => turning off Rayleigh! ****
			// '//);
		}
		if (EGS4Geom.iterseindex < 0 || EGS4Geom.iterseindex > 2) {
			EGS4.STOPPROGRAM = true;
			EGS4.seqStr = " ******** ERROR in geometrical inputs! **** ";
			printSequence(EGS4.seqStr);
			return;
		}
		if (EGS4Geom.iterseindex == EGS4Geom.iCAVITY_INFORMATION)// 2
		{
			if ((EGS4Geom.SLENGHT).compareTo(NULLs) == 0
					|| (EGS4Geom.airs).compareTo(NULLs) == 0) {
				EGS4.STOPPROGRAM = true;
				EGS4.seqStr = " ******** ERROR in geometrical inputs! **** ";
				printSequence(EGS4.seqStr);
				return;
			}
			if (EGS4Geom.ELERAD != 0.0
					&& (EGS4Geom.electrods).compareTo(NULLs) == 0) {
				EGS4.STOPPROGRAM = true;
				EGS4.seqStr = " ******** ERROR in geometrical inputs! **** ";
				printSequence(EGS4.seqStr);
				return;
			}
		}

		if (!setEcutRegion) {
			if (EGS4.ECUT[0] < ECUT_MIN || EGS4.ECUT[0] > ECUT_MAX)
				EGS4.ECUT[0] = EGS4.$GLOBAL_ECUT;
			for (int I = 2; I <= EGS4.$MXREG; I++) {
				EGS4.ECUT[I - 1] = EGS4.ECUT[0];
			}

		}
		if (!setPcutRegion) {
			if (EGS4.PCUT[0] < PCUT_MIN || EGS4.PCUT[0] > PCUT_MAX)
				EGS4.PCUT[0] = EGS4.$GLOBAL_PCUT;
			// "Now set ecut and pcut to the values input by the user"
			for (int I = 2; I <= EGS4.$MXREG; I++) {
				EGS4.PCUT[I - 1] = EGS4.PCUT[0];
			}
		}
		if (!setSmaxRegion) {
			if (EGS4.SMAXIR[0] < SMAXIR_MIN || EGS4.SMAXIR[0] > SMAXIR_MAX)
				EGS4.SMAXIR[0] = EGS4.$MAX_SMAX;
			for (int I = 2; I <= EGS4.$MXREG; I++) {
				EGS4.SMAXIR[I - 1] = EGS4.SMAXIR[0];
			}
		}

		if (setEcutRegion) {
			for (int i = 1; i <= nEcut; i++) {
				if (Ecut[i - 1] < ECUT_MIN || Ecut[i - 1] > ECUT_MAX)
					Ecut[i - 1] = EGS4.$GLOBAL_ECUT;
				if (startEcutRegion[i - 1] < 1
						|| startEcutRegion[i - 1] > EGS4.$MXREG)
					startEcutRegion[i - 1] = 1;
				if (stopEcutRegion[i - 1] < 1
						|| stopEcutRegion[i - 1] > EGS4.$MXREG)
					stopEcutRegion[i - 1] = 1;

				for (int j = startEcutRegion[i - 1]; j <= stopEcutRegion[i - 1]; j++)
					EGS4.ECUT[j - 1] = Ecut[i - 1];
			}
		}
		if (setPcutRegion) {
			for (int i = 1; i <= nPcut; i++) {
				if (Pcut[i - 1] < PCUT_MIN || Pcut[i - 1] > PCUT_MAX)
					Pcut[i - 1] = EGS4.$GLOBAL_PCUT;
				if (startPcutRegion[i - 1] < 1
						|| startPcutRegion[i - 1] > EGS4.$MXREG)
					startPcutRegion[i - 1] = 1;
				if (stopPcutRegion[i - 1] < 1
						|| stopPcutRegion[i - 1] > EGS4.$MXREG)
					stopPcutRegion[i - 1] = 1;

				for (int j = startPcutRegion[i - 1]; j <= stopPcutRegion[i - 1]; j++)
					EGS4.PCUT[j - 1] = Pcut[i - 1];
			}
		}
		if (setSmaxRegion) {
			for (int i = 1; i <= nSmax; i++) {
				if (Smax[i - 1] < SMAXIR_MIN || Smax[i - 1] > SMAXIR_MAX)
					Smax[i - 1] = EGS4.$MAX_SMAX;
				if (startSmaxRegion[i - 1] < 1
						|| startSmaxRegion[i - 1] > EGS4.$MXREG)
					startSmaxRegion[i - 1] = 1;
				if (stopSmaxRegion[i - 1] < 1
						|| stopSmaxRegion[i - 1] > EGS4.$MXREG)
					stopSmaxRegion[i - 1] = 1;

				for (int j = startSmaxRegion[i - 1]; j <= stopSmaxRegion[i - 1]; j++)
					EGS4.SMAXIR[j - 1] = Smax[i - 1];
			}
		}

		if (!setIncohRegion) {
			EGS4.ibcmp[0] = incoh;
			if (EGS4.ibcmp[0] < 0 || EGS4.ibcmp[0] > 3) {
				EGS4.ibcmp[0] = EGS4.$IBCMP_DEFAULT;
				incoh = EGS4.$IBCMP_DEFAULT;
			}
			// "Now set ibcmp for all regions to the value input by the user"
			for (int I = 2; I <= EGS4.$MXREG; I++) {
				EGS4.ibcmp[I - 1] = EGS4.ibcmp[0];
			}
		}
		if (setIncohRegion) {
			for (int i = 1; i <= nIncoh; i++) {
				if (Incoh[i - 1] < 0 || Incoh[i - 1] > 1)// only 2 values
					Incoh[i - 1] = EGS4.$IBCMP_DEFAULT;
				if (startIncohRegion[i - 1] < 1
						|| startIncohRegion[i - 1] > EGS4.$MXREG)
					startIncohRegion[i - 1] = 1;
				if (stopIncohRegion[i - 1] < 1
						|| stopIncohRegion[i - 1] > EGS4.$MXREG)
					stopIncohRegion[i - 1] = 1;

				for (int j = startIncohRegion[i - 1]; j <= stopIncohRegion[i - 1]; j++)
					EGS4.ibcmp[j - 1] = Incoh[i - 1];
			}
		}
		if (!setCohRegion) {
			EGS4.IRAYLR[0] = coh;
			if (EGS4.IRAYLR[0] < 0 || EGS4.IRAYLR[0] > 3) {
				EGS4.IRAYLR[0] = EGS4.$IRAYLR_DEFAULT;
				coh = EGS4.$IRAYLR_DEFAULT;
			}
			// "Now set iraylr for all regions to the value input by the user"
			for (int I = 2; I <= EGS4.$MXREG; I++) {
				EGS4.IRAYLR[I - 1] = EGS4.IRAYLR[0];
			}
		}
		if (setCohRegion) {
			for (int i = 1; i <= nCoh; i++) {
				if (Coh[i - 1] < 0 || Coh[i - 1] > 1)// only 2 values
					Coh[i - 1] = EGS4.$IRAYLR_DEFAULT;
				if (startCohRegion[i - 1] < 1
						|| startCohRegion[i - 1] > EGS4.$MXREG)
					startCohRegion[i - 1] = 1;
				if (stopCohRegion[i - 1] < 1
						|| stopCohRegion[i - 1] > EGS4.$MXREG)
					stopCohRegion[i - 1] = 1;

				for (int j = startCohRegion[i - 1]; j <= stopCohRegion[i - 1]; j++)
					EGS4.IRAYLR[j - 1] = Coh[i - 1];
			}
		}
		if (!setRelaxRegion) {
			EGS4.iedgfl[0] = relax;
			if (EGS4.iedgfl[0] < 0 || EGS4.iedgfl[0] > 3) {
				EGS4.iedgfl[0] = EGS4.$IEDGFL_DEFAULT;
				relax = EGS4.$IEDGFL_DEFAULT;
			}
			// "Now set iedgfl for all regions to the value input by the user"
			for (int I = 2; I <= EGS4.$MXREG; I++) {
				EGS4.iedgfl[I - 1] = EGS4.iedgfl[0];
			}
		}
		if (setRelaxRegion) {
			for (int i = 1; i <= nRelax; i++) {
				if (Relax[i - 1] < 0 || Relax[i - 1] > 1)// only 2 values
					Relax[i - 1] = EGS4.$IEDGFL_DEFAULT;
				if (startRelaxRegion[i - 1] < 1
						|| startRelaxRegion[i - 1] > EGS4.$MXREG)
					startRelaxRegion[i - 1] = 1;
				if (stopRelaxRegion[i - 1] < 1
						|| stopRelaxRegion[i - 1] > EGS4.$MXREG)
					stopRelaxRegion[i - 1] = 1;

				for (int j = startRelaxRegion[i - 1]; j <= stopRelaxRegion[i - 1]; j++)
					EGS4.iedgfl[j - 1] = Relax[i - 1];
			}
		}
		if (!setPeRegion) {
			EGS4.iphter[0] = pe;
			if (EGS4.iphter[0] < 0 || EGS4.iphter[0] > 3) {
				EGS4.iphter[0] = EGS4.$IPHTER_DEFAULT;
				pe = EGS4.$IPHTER_DEFAULT;
			}
			// "Now set iphter for all regions to the value input by the user"
			for (int I = 2; I <= EGS4.$MXREG; I++) {
				EGS4.iphter[I - 1] = EGS4.iphter[0];
			}
		}
		if (setPeRegion) {
			for (int i = 1; i <= nPe; i++) {
				if (Pe[i - 1] < 0 || Pe[i - 1] > 1)// only 2 values
					Pe[i - 1] = EGS4.$IPHTER_DEFAULT;
				if (startPeRegion[i - 1] < 1
						|| startPeRegion[i - 1] > EGS4.$MXREG)
					startPeRegion[i - 1] = 1;
				if (stopPeRegion[i - 1] < 1
						|| stopPeRegion[i - 1] > EGS4.$MXREG)
					stopPeRegion[i - 1] = 1;

				for (int j = startPeRegion[i - 1]; j <= stopPeRegion[i - 1]; j++)
					EGS4.iphter[j - 1] = Pe[i - 1];
			}
		}
		if ((EGS4.ibrdst < 0) || (EGS4.ibrdst > 1)) {
			EGS4.ibrdst = EGS4.$IBRDST_DEFAULT;
		}
		if ((EGS4.iprdst < 0) || (EGS4.iprdst > 4)) {
			EGS4.iprdst = EGS4.$IPRDST_DEFAULT;
		}
		if ((EGS4.estepe < ESTEPE_MIN) || (EGS4.estepe > ESTEPE_MAX)) {
			EGS4.estepe = EGS4.$MAX_ELOSS; // "$MAX-ELOSS is defined in egsnrc.macros at 0.25"
		}
		if ((EGS4.ximax < XIMAX_MIN) || (EGS4.ximax > XIMAX_MAX)) {
			EGS4.ximax = EGS4.$EXACT_BCA_XIMAX; // "$EXACT-BCA-XIMAX set to 0.5 in egsnrc.macros"
		}
		if ((EGS4.bca_algorithm < 0) || (EGS4.bca_algorithm > 1)) {
			EGS4.bca_algorithm = EGS4.$BCA_ALGORITHM_DEFAULT;
		}
		if ((EGS4.skindepth_for_bca < Skindepth_MIN)
				|| (EGS4.skindepth_for_bca > Skindepth_MAX))
			EGS4.skindepth_for_bca = EGS4.$SKIN_DEPTH_FOR_BCA;
		if (EGS4.bca_algorithm == BCA_EXACT) {
			if (EGS4.skindepth_for_bca <= 0.0)
				EGS4.skindepth_for_bca = EGS4.$SKIN_DEPTH_FOR_BCA;
		}
		if ((EGS4.transport_algorithm < 0) || (EGS4.transport_algorithm > 1)) {
			EGS4.transport_algorithm = EGS4.$TRANSPORT_ALGORITHM_DEFAULT;
		}
		if (ispin == spin_ON)
			EGS4.spin_effects = true;// so=>
		else
			EGS4.spin_effects = false;// so=>
		if ((EGS4.ibr_nist < 0) || (EGS4.ibr_nist > 1)) {
			EGS4.ibr_nist = EGS4.$IBR_NIST_DEFAULT;
		}
		if ((EGS4.pair_nrc < 0) || (EGS4.pair_nrc > 1)) {
			EGS4.pair_nrc = EGS4.$PAIR_NRC_DEFAULT;
		}
		if ((EGS4.eii_flag < 0) || (EGS4.eii_flag > 4)) {
			EGS4.eii_flag = eii_OFF;// default
		}

		if (EGS4.eii_flag == eii_ON)
			EGS4.eiifile = "eii_ik";
		else if (EGS4.eii_flag == eii_casnati)
			EGS4.eiifile = "eii_casnati";
		else if (EGS4.eii_flag == eii_kolbenstvedt)
			EGS4.eiifile = "eii_kolbenstvedt";
		else if (EGS4.eii_flag == eii_gryzinski)
			EGS4.eiifile = "eii_gryzinski";

		if (photon_xsection.compareTo("") != 0)
			if (photon_xsection.compareTo(EGS4.photon_xsections_epdl) != 0)
				if (photon_xsection.compareTo(EGS4.photon_xsections_xcom) != 0)
					photon_xsection = "";// reset to default

		EGS4.photon_xsections = photon_xsection;

		if ((EGS4.itriplet < 0) || (EGS4.itriplet > 1)) {
			EGS4.itriplet = EGS4.$TRIPLET_DEFAULT;
		}
		if ((EGS4.radc_flag < 0) || (EGS4.radc_flag > 1)) {
			EGS4.radc_flag = radc_OFF;// default
		}

		if (EGS4Macro.irejct < 0 || EGS4Macro.irejct > 1) {
			EGS4Macro.irejct = irejct_OFF;
		}// default
		if (ESAVEIN < ESAVEIN_MIN || ESAVEIN > ESAVEIN_MAX) {
			ESAVEIN = ESAVEIN_DEFAULT;
		}
		if (EGS4Macro.cs_enhance < cs_enhance_MIN
				|| EGS4Macro.cs_enhance > cs_enhance_MAX) {
			EGS4Macro.cs_enhance = cs_enhance_DEFAULT;
		}
		if (EGS4Macro.cs_enhance > 1.) {
			EGS4Macro.use_enhance = true;
		} else {
			EGS4Macro.use_enhance = false;
		}

		EGS4.seqStr = " *** Reading variance reduction inputs ... ***";
		if (EGS4.iprint > 1)
			printSequence(EGS4.seqStr);

		EGS4.seqStr = " Range rejection is On(1) or Off(0):" + EGS4Macro.irejct;
		if (EGS4.iprint > 1)
			printSequence(EGS4.seqStr);

		if (EGS4Macro.irejct > 0) {
			EGS4.seqStr = " ESAVEIN cutoff value(total) for range rejection:"
					+ EGS4.format(ESAVEIN, 10, true) + " MeV";
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);

			if (ESAVEIN == 0.0) {
				EGS4.seqStr = " WARNING: Have asked for range rejection but left ESAVEIN=0.0!!";
				if (EGS4.iprint > 1)
					printSequence(EGS4.seqStr);
				EGS4.seqStr = " There will be no range rejection!";
				if (EGS4.iprint > 1)
					printSequence(EGS4.seqStr);
			}
			for (int i = 1; i <= EGS4Geom.nreg; i++) {
				EGS4.i_do_rr[i - 1] = 1;
				EGS4.e_max_rr[i - 1] = ESAVEIN;
			}
			// "note  e_max_r is total energy"
			// "above two arrays needed for each region for EGSnrc RANGE-DISCARD macro"
		}
		if (RRZ < RRDEPTH_MIN || RRZ > RRDEPTH_MAX) {
			RRZ = RRDEPTH_DEFAULT;
		}
		if (RRCUT < RRFRACTION_MIN || RRCUT > RRFRACTION_MAX) {
			RRCUT = RRFRACTION_DEFAULT;
		}
		if (EGS4Macro.CEXPTR < EXPC_MIN || EGS4Macro.CEXPTR > EXPC_MAX) {
			EGS4Macro.CEXPTR = EXPC_DEFAULT;
		}
		RUSROU = false;
		if (RRZ + RRCUT != 0.0)
			RUSROU = true;
		if (RUSROU) {
			EGS4.seqStr = " RUSSIAN ROULETTE WILL BE PLAYED";
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);
			EGS4.seqStr = " RUSSIAN ROULETTE PLANE:"
					+ EGS4.format(RRZ, 14, false);
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);
			EGS4.seqStr = " SURVIVAL PROBABILITY:"
					+ EGS4.format(RRCUT, 14, false);
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);
		} else {
			EGS4.seqStr = " RUSSIAN ROULETTE WILL NOT BE PLAYED";
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);
		}
		if (EGS4Macro.CEXPTR == 0.) {
			EGS4.seqStr = " NO PATHLENGTH BIASING TO BE DONE";
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);
		} else {
			EGS4.seqStr = " CEXPTR PARAMATER:"
					+ EGS4.format(EGS4Macro.CEXPTR, 14, false);
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);
		}
		if (IFARCE < 0 || IFARCE > 1) {
			IFARCE = IFARCE_OFF;
		}// default
		NFMIN_MAX = EGS4Geom.nreg;
		NFMAX_MAX = EGS4Geom.nreg + 1;
		if (EGS4Macro.NFMIN < NFMIN_MIN || EGS4Macro.NFMIN > NFMIN_MAX) {
			EGS4Macro.NFMIN = NFMIN_DEFAULT;
		}
		if (EGS4Macro.NFMAX < NFMAX_MIN || EGS4Macro.NFMAX > NFMAX_MAX) {
			EGS4Macro.NFMAX = NFMAX_DEFAULT;
		}

		if (IFARCE == 0) {
			EGS4Macro.IFORCE = 0;
			EGS4Macro.NFMIN = 0;
			EGS4Macro.NFMAX = 0;

			EGS4.seqStr = " NO INTERACTION FORCING IS IN EFFECT";
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);
		} else if (IFARCE == 1) {
			EGS4Macro.IFORCE = 1;
			if (EGS4Macro.NFMAX < EGS4Macro.NFMIN)
				EGS4Macro.NFMAX = EGS4Macro.NFMIN;

			EGS4.seqStr = " FORCED PHOTON INTERACTIONS IN EFFECT FROM "
					+ EGS4Macro.NFMIN + " TO " + EGS4Macro.NFMAX
					+ " INTERACTIONS ";
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);
		}
		EGS4Macro.IQINC = EGS4SrcEns.iqin; // "NEEDED TO TURN OFF FASTSTEP FOR INCIDENT ELECTRONS"
		// "WHEN FORCING INTERACTIONS"
		phsplitt_MAX = EGS4.$MXSTACK - 2;
		if (phsplitt < phsplitt_MIN || phsplitt > phsplitt_MAX) {
			phsplitt = phsplitt_DEFAULT;
		}
		EGS4Macro.n_split = phsplitt;
		if (EGS4Macro.n_split > 1) {
			if (EGS4Macro.IFULL > 1) {
				EGS4.seqStr = " IGNORING INPUT: Photon splitting only for ifull = 0,1! ";
				if (EGS4.iprint > 1)
					printSequence(EGS4.seqStr);

				EGS4Macro.n_split = 1;
			} else {
				EGS4.seqStr = " Calculation with photon splitting, n_split = "
						+ EGS4Macro.n_split;
				if (EGS4.iprint > 1)
					printSequence(EGS4.seqStr);

				EGS4Macro.iifano = IFANO;
			}
		}

		if (EGS4Macro.use_enhance) {
			EGS4.seqStr = "  Calculation with CS enhancement  ";
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);
			EGS4.seqStr = "  photon forcing, exp. transform, etc. input will be ignored, IFULL will be set to 1! (i.e. Ascat and Aatt) !";
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);
			EGS4.seqStr = "  Using cs_enhance = " + EGS4Macro.cs_enhance;
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);

			EGS4Macro.IFULL = 1;
			EGS4Macro.IFORCE = 0;
			EGS4Macro.NFMIN = 0;
			EGS4Macro.NFMAX = 0;
			EGS4Macro.n_split = 1;
		}

		// "if we are calculating Aatt and Ascatt, we must have n_split>1,"
		// "cs_enhance on, or photon forcing on"
		if (EGS4Macro.IFULL == 1 && EGS4Macro.cs_enhance <= 1
				&& EGS4Macro.IFORCE == 0 && EGS4Macro.n_split <= 1) {
			EGS4Macro.IFORCE = 1;
			EGS4Macro.NFMIN = 1;
			EGS4Macro.NFMAX = 1;

			EGS4.seqStr = " If you are calculating Aatt, Ascat (IFULL=1), you must be using";
			if (EGS4.iprint > 0)
				printSequence(EGS4.seqStr);
			EGS4.seqStr = " photon forcing, photon splitting, or cross-section enhancement.";
			if (EGS4.iprint > 0)
				printSequence(EGS4.seqStr);
			EGS4.seqStr = " Currently, none of these are being used.  Will continue run with";
			if (EGS4.iprint > 0)
				printSequence(EGS4.seqStr);
			EGS4.seqStr = " photon forcing on and one interaction forced (NFMIN=NFMAX=1).";
			if (EGS4.iprint > 0)
				printSequence(EGS4.seqStr);
		}

		if (IOOPTN == 1 && (EGS4Macro.n_split > 1 || EGS4Macro.use_enhance)) {
			IOOPTN = 0;
			EGS4.seqStr = " You cannot have cross-section enhancement or photon splitting on ";
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);
			EGS4.seqStr = " and still output detailed results for each cavity zone.  IOOPTN reset to 0.";
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);
		}

	}

	@SuppressWarnings("unused")
	/**
	 * Initialize random number generator.
	 */
	private void init_random_generator() {
		if (EGS4.ranluxB) {
			if (jrng1 < RANLUX_LEVEL_MIN || jrng1 > RANLUX_LEVEL_MAX) {
				jrng1 = RANLUX_LEVEL_DEFAULT;
			}
			if (jrng2 < RANLUX_SEED_MIN || jrng2 > RANLUX_SEED_MAX) {
				jrng2 = RANLUX_SEED_DEFAULT;
			}
		} else {
			if (jrng1 < RANMAR_SEED_MIN || jrng1 > RANMAR_SEED_MAX) {
				jrng1 = RANMAR_SEED_DEFAULT;
			}
			if (jrng2 < RANMAR_SEED_MIN || jrng2 > RANMAR_SEED_MAX) {
				jrng2 = RANMAR_SEED_DEFAULT;
			}
		}
		// IF( i_parallel > 0 ) jrng2 = jrng2 - 1 + i_parallel;//NO PARALLEL JOB
		// USED!!!
		// $INITIALIZE RNG USING jrng1 AND jrng2;
		if (EGS4.ranluxB) {
			EGS4.init_ranlux(jrng1, jrng2);
			EGS4.ranlux(EGS4.rng_array);
			EGS4.rng_seed = 1;
		} else {
			EGS4.ixx = jrng1;
			EGS4.jxx = jrng2;
			EGS4.init_ranmar();
		}
	}

	/**
	 * Print input summary.
	 */
	private void ISUMRY() {
		String s = "";
		String s1 = "";
		int ll = 0;
		int ioff = 0;
		EGS4.seqStr = "************************INPUT SUMMARY:*****************************************";
		if (EGS4.iprint > 1)
			printSequence(EGS4.seqStr);

		EGS4.seqStr = "=========================================================================";
		if (EGS4.iprint > 1)
			printSequence(EGS4.seqStr);
		EGS4.seqStr = "                   Electron/Photon transport parameter";
		if (EGS4.iprint > 1)
			printSequence(EGS4.seqStr);
		EGS4.seqStr = "=========================================================================";
		if (EGS4.iprint > 1)
			printSequence(EGS4.seqStr);
		s = " Photon transport cutoff(MeV)";
		ioff = 32;
		if (setPcutRegion) {
			s = s + EGS4.format("", ioff) + "Set in regions";
		} else {
			if (EGS4.PCUT[0] > 1.e-4) {
				ioff = 30;
				s = s + EGS4.format("", ioff)
						+ EGS4.format(EGS4.PCUT[0], 8, true);
			} else {
				s = s + EGS4.format("", ioff) + "AP(medium)";
			}
		}
		EGS4.seqStr = s;
		if (EGS4.iprint > 1)
			printSequence(EGS4.seqStr);
		s = " Pair angular sampling";
		ioff = 39;
		if (EGS4.iprdst == pair_ang_OFF) {
			s = s + EGS4.format("", ioff) + "OFF";
		} else if (EGS4.iprdst == pair_ang_SIMPLE) {
			s = s + EGS4.format("", ioff) + "SIMPLE";
		} else if (EGS4.iprdst == pair_ang_KM) {
			s = s + EGS4.format("", ioff) + "KM";
		} else if (EGS4.iprdst == pair_ang_UNIFORM) {
			s = s + EGS4.format("", ioff) + "UNIFORM";
		} else if (EGS4.iprdst == pair_ang_BLEND) {
			s = s + EGS4.format("", ioff) + "BLEND";
		}
		EGS4.seqStr = s;
		if (EGS4.iprint > 1)
			printSequence(EGS4.seqStr);
		s = " Pair cross sections";
		ioff = 41;
		if (EGS4.pair_nrc == pair_cross_BH) {
			s = s + EGS4.format("", ioff) + "BH";
		} else if (EGS4.pair_nrc == pair_cross_NRC) {
			s = s + EGS4.format("", ioff) + "NRC";
		}
		EGS4.seqStr = s;
		if (EGS4.iprint > 1)
			printSequence(EGS4.seqStr);
		s = " Triplet production";
		ioff = 42;
		if (EGS4.itriplet == triplet_OFF) {
			s = s + EGS4.format("", ioff) + "OFF";
		} else if (EGS4.itriplet == triplet_ON) {
			s = s + EGS4.format("", ioff) + "ON";
		}
		EGS4.seqStr = s;
		if (EGS4.iprint > 1)
			printSequence(EGS4.seqStr);
		s = " Bound Compton scattering";
		ioff = 36;
		if (incoh == incoh_OFF) {
			s = s + EGS4.format("", ioff) + "OFF";
		} else if (incoh == incoh_ON) {
			s = s + EGS4.format("", ioff) + "ON";
		} else if (incoh == incoh_ON_IN_REGIONS) {
			s = s + EGS4.format("", ioff) + "ON IN REGIONS";
		} else if (incoh == incoh_OFF_IN_REGIONS) {
			s = s + EGS4.format("", ioff) + "OFF IN REGIONS";
		}
		EGS4.seqStr = s;
		if (EGS4.iprint > 1)
			printSequence(EGS4.seqStr);
		s = " Radiative Compton corrections";
		ioff = 31;
		if (EGS4.radc_flag == radc_OFF) {
			s = s + EGS4.format("", ioff) + "OFF";
		} else if (EGS4.radc_flag == radc_ON) {
			s = s + EGS4.format("", ioff) + "ON";
		}
		EGS4.seqStr = s;
		if (EGS4.iprint > 1)
			printSequence(EGS4.seqStr);
		s = " Rayleigh scattering";
		ioff = 41;
		if (coh == coh_OFF) {
			s = s + EGS4.format("", ioff) + "OFF";
		} else if (coh == coh_ON) {
			s = s + EGS4.format("", ioff) + "ON";
		} else if (coh == coh_ON_IN_REGIONS) {
			s = s + EGS4.format("", ioff) + "ON IN REGIONS";
		} else if (coh == coh_OFF_IN_REGIONS) {
			s = s + EGS4.format("", ioff) + "OFF IN REGIONS";
		}
		EGS4.seqStr = s;
		if (EGS4.iprint > 1)
			printSequence(EGS4.seqStr);
		s = " Atomic relaxations";
		ioff = 42;
		if (relax == relax_OFF) {
			s = s + EGS4.format("", ioff) + "OFF";
		} else if (relax == relax_ON) {
			s = s + EGS4.format("", ioff) + "ON";
		} else if (relax == relax_ON_IN_REGIONS) {
			s = s + EGS4.format("", ioff) + "ON IN REGIONS";
		} else if (relax == relax_OFF_IN_REGIONS) {
			s = s + EGS4.format("", ioff) + "OFF IN REGIONS";
		}
		EGS4.seqStr = s;
		if (EGS4.iprint > 1)
			printSequence(EGS4.seqStr);
		s = " Photoelectron angular sampling";
		ioff = 30;
		if (pe == pe_ang_OFF) {
			s = s + EGS4.format("", ioff) + "OFF";
		} else if (pe == pe_ang_ON) {
			s = s + EGS4.format("", ioff) + "ON";
		} else if (pe == pe_ang_ON_IN_REGIONS) {
			s = s + EGS4.format("", ioff) + "ON IN REGIONS";
		} else if (pe == pe_ang_OFF_IN_REGIONS) {
			s = s + EGS4.format("", ioff) + "OFF IN REGIONS";
		}
		EGS4.seqStr = s;
		if (EGS4.iprint > 1)
			printSequence(EGS4.seqStr);
		EGS4.seqStr = "";
		if (EGS4.iprint > 1)
			printSequence(EGS4.seqStr);
		s = " Electron transport cutoff(MeV)";
		ioff = 30;
		if (setEcutRegion) {
			s = s + EGS4.format("", ioff) + "Set in regions";
		} else {
			if (EGS4.ECUT[0] > 1.e-4) {
				ioff = 28;
				s = s + EGS4.format("", ioff)
						+ EGS4.format(EGS4.ECUT[0], 8, true);
			} else {
				s = s + EGS4.format("", ioff) + "AE(medium)";
			}
		}
		EGS4.seqStr = s;
		if (EGS4.iprint > 1)
			printSequence(EGS4.seqStr);
		s = " Bremsstrahlung angular sampling";
		ioff = 29;
		if (EGS4.ibrdst == brems_ang_SIMPLE) {
			s = s + EGS4.format("", ioff) + "SIMPLE";
		} else if (EGS4.ibrdst == brems_ang_KM) {
			s = s + EGS4.format("", ioff) + "KM";
		}
		EGS4.seqStr = s;
		if (EGS4.iprint > 1)
			printSequence(EGS4.seqStr);
		s = " Bremsstrahlung cross sections";
		ioff = 31;
		if (EGS4.ibr_nist == brems_cross_BH) {
			s = s + EGS4.format("", ioff) + "BH";
		} else if (EGS4.ibr_nist == brems_cross_NIST) {
			s = s + EGS4.format("", ioff) + "NIST";
		}
		EGS4.seqStr = s;
		if (EGS4.iprint > 1)
			printSequence(EGS4.seqStr);
		s = " Spin effects";
		ioff = 48;
		if (EGS4.spin_effects) {
			s = s + EGS4.format("", ioff) + "ON";
		} else {
			s = s + EGS4.format("", ioff) + "OFF";
		}
		EGS4.seqStr = s;
		if (EGS4.iprint > 1)
			printSequence(EGS4.seqStr);
		s = " Electron Impact Ionization";
		ioff = 34;
		if (EGS4.eii_flag == eii_OFF) {
			s = s + EGS4.format("", ioff) + "OFF";
		} else if (EGS4.eii_flag == eii_ON) {
			s = s + EGS4.format("", ioff) + "ON";
		}
		EGS4.seqStr = s;
		if (EGS4.iprint > 1)
			printSequence(EGS4.seqStr);
		s = " Maxium electron step in cm (SMAX)";
		ioff = 27;
		if (setSmaxRegion) {
			s = s + EGS4.format("", ioff) + "Set in regions";
		} else {
			if (EGS4.SMAXIR[0] > 1.e-4) {
				ioff = 26;
				s = s + EGS4.format("", ioff)
						+ EGS4.format(EGS4.SMAXIR[0], 6, false);
			} else {
				s = s + EGS4.format("", ioff) + "Restriction is off";
			}
		}
		EGS4.seqStr = s;
		if (EGS4.iprint > 1)
			printSequence(EGS4.seqStr);
		s = " Maximum fractional energy loss/step (ESTEPE)";
		ioff = 16;
		s = s + EGS4.format("", ioff) + EGS4.format(EGS4.estepe, 6, true);
		EGS4.seqStr = s;
		if (EGS4.iprint > 1)
			printSequence(EGS4.seqStr);
		s = " Maximum 1st elastic moment/step (XIMAX)";
		ioff = 21;
		s = s + EGS4.format("", ioff) + EGS4.format(EGS4.ximax, 6, true);
		EGS4.seqStr = s;
		if (EGS4.iprint > 1)
			printSequence(EGS4.seqStr);
		s = " Boundary crossing algorithm";
		ioff = 33;
		if (EGS4.bca_algorithm == BCA_EXACT) {
			s = s + EGS4.format("", ioff) + "EXACT";
		} else if (EGS4.bca_algorithm == BCA_PRESTA_I) {
			s = s + EGS4.format("", ioff) + "PRESTA I";
		}
		EGS4.seqStr = s;
		if (EGS4.iprint > 1)
			printSequence(EGS4.seqStr);
		s = " Skin-depth for boundary crossing (MFP)";
		ioff = 22;
		s = s + EGS4.format("", ioff)
				+ EGS4.format(EGS4.skindepth_for_bca, 6, true);
		EGS4.seqStr = s;
		if (EGS4.iprint > 1)
			printSequence(EGS4.seqStr);
		s = " Electron-step algorithm";
		ioff = 37;
		if (EGS4.transport_algorithm == estep_alg_PRESTA_II) {
			s = s + EGS4.format("", ioff) + "PRESTA II";
		} else if (EGS4.transport_algorithm == estep_alg_PRESTA_I) {
			s = s + EGS4.format("", ioff) + "PRESTA I";
		}
		EGS4.seqStr = s;
		if (EGS4.iprint > 1)
			printSequence(EGS4.seqStr);
		EGS4.seqStr = "=========================================================================";
		if (EGS4.iprint > 1)
			printSequence(EGS4.seqStr);
		EGS4.seqStr = "                   MONTE CARLO, TRANSPORT, AND SCATTER CONTROLS";
		if (EGS4.iprint > 1)
			printSequence(EGS4.seqStr);
		EGS4.seqStr = "=========================================================================";
		if (EGS4.iprint > 1)
			printSequence(EGS4.seqStr);
		s = " Max. # of histories to RUN";
		ll = s.length();
		ll = 54 - ll;
		s = s + EGS4.format("", ll);
		EGS4.seqStr = s + EGS4.format(NCASE, 12);
		if (EGS4.iprint > 1)
			printSequence(EGS4.seqStr);
		s = " Max # of histories to ANALYZE";
		ll = s.length();
		ll = 54 - ll;
		s = s + EGS4.format("", ll);
		EGS4.seqStr = s + EGS4.format(NCASET, 12);
		if (EGS4.iprint > 1)
			printSequence(EGS4.seqStr);
		s = " Incident Charge";// EGS4.format("",20)+"Incident Charge";
		s1 = "";
		ll = s.length();
		ll = 68 - ll;
		if (EGS4SrcEns.iqin == 0)
			s1 = "photons";
		if (EGS4SrcEns.iqin == -1)
			s1 = "electrons";
		if (EGS4SrcEns.iqin == 1)
			s1 = "positrons";
		if (EGS4SrcEns.iqin == 2)
			s1 = "all";
		if (EGS4SrcEns.iqin == 3)
			s1 = "e- & e+";
		s = s + EGS4.format(s1, ll);
		EGS4.seqStr = s;
		if (EGS4.iprint > 1)
			printSequence(EGS4.seqStr);

		if (EGS4SrcEns.MONOEN == 0 && EGS4SrcEns.ISOURC != 21
				&& EGS4SrcEns.ISOURC != 22 && EGS4SrcEns.ISOURC != 23) {
			s = " Incident kinetic energy:";
			ll = s.length();
			ll = 57 - ll;
			s = s + EGS4.format("", ll);
			EGS4.seqStr = s + EGS4.format(EGS4SrcEns.ein, 9, true) + " MeV";
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);

		} else {
			EGS4SrcEns.ENSRCO();
		}
		boolean gotopast = false;
		for (int I = 2; I <= EGS4Geom.nreg; I++) {
			gotopast = false;
			if ((EGS4.ECUT[I - 1] != EGS4.ECUT[1])
					|| (EGS4.PCUT[I - 1] != EGS4.PCUT[1])) {
				// "we failed at least one test, so this means there really are"
				// "varying ECUTs and these will be printed in the grid if we want them"
				// "print the first 12 ECUT & PCUT just to be sure"
				int j = Math.min(12, EGS4Geom.nreg);

				for (int JJ = 2; JJ <= j; JJ++) {
					EGS4.seqStr = "First ECUTs (MeV): "
							+ EGS4.format(EGS4.ECUT[JJ - 1], 12, true) + "  ,"
							+ "First PCUTs (MeV): "
							+ EGS4.format(EGS4.PCUT[JJ - 1], 12, true);
					if (EGS4.iprint > 1)
						printSequence(EGS4.seqStr);
				}
				// GO TO :past:
				gotopast = true;
				break;
			}
		}
		// "if we get here, they were all the same"
		if (!gotopast) {
			s = " GLOBAL ELECTRON TRANSPORT CUT-OFF";
			ll = s.length();
			ll = 57 - ll;
			s = s + EGS4.format("", ll);
			EGS4.seqStr = s + EGS4.format(EGS4.ECUT[1], 9, true) + " MeV";
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);
			s = " GLOBAL PHOTON TRANSPORT CUT-OFF";
			ll = s.length();
			ll = 57 - ll;
			s = s + EGS4.format("", ll);
			EGS4.seqStr = s + EGS4.format(EGS4.PCUT[1], 9, true) + " MeV";
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);

		}
		// :past:
		if (EGS4Macro.IFORCE != 0) {
			s = " Min/max photon step forced";
			ll = s.length();
			ll = 56 - ll;
			s = s + EGS4.format("", ll);
			EGS4.seqStr = s + EGS4.format(EGS4Macro.NFMIN, 6) + "/"
					+ EGS4.format(EGS4Macro.NFMAX, 6);
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);

		} else {
			s = " Photon force interaction switch";
			ll = s.length();
			ll = 56 - ll;
			s = s + EGS4.format("", ll);
			EGS4.seqStr = s + "OFF";
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);
		}
		if (EGS4Macro.irejct > 0) {
			s = " Range rejection on a region by region basis";
			// ll=s.length();ll=56-ll;s=s+EGS4.format("",ll);
			EGS4.seqStr = s;
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);
			s = " Also globally to region between z="
					+ EGS4.format(EGS4Macro.z_cavity_min, 10, true) + " &"
					+ EGS4.format(EGS4Macro.z_cavity_max, 10, true) + " cm";
			// ll=s.length();ll=56-ll;s=s+EGS4.format("",ll);
			EGS4.seqStr = s;
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);
			s = "      and inside radius="
					+ EGS4.format(EGS4Macro.r_cavity_max, 12, true) + " cm";
			// ll=s.length();ll=56-ll;s=s+EGS4.format("",ll);
			EGS4.seqStr = s;
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);
			s = " Range rejection only for electrons < ESAVEIN="
					+ EGS4.format(ESAVEIN, 10, true) + " MeV";
			// ll=s.length();ll=56-ll;s=s+EGS4.format("",ll);
			EGS4.seqStr = s;
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);
		} else {
			s = " RANGE REJECTION SWITCH";
			ll = s.length();
			ll = 61 - ll;
			s = s + EGS4.format("", ll);
			EGS4.seqStr = s + "OFF";
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);
		}
		// WRITE(IOUT,260) TIMMAX,STATLM;
		s = " Maximum cputime allowed";
		ll = s.length();
		ll = 60 - ll;
		s = s + EGS4.format("", ll);
		EGS4.seqStr = s + EGS4.format(TIMMAX, 6, true) + " hrs";
		if (EGS4.iprint > 1)
			printSequence(EGS4.seqStr);
		s = " Stats in cavity objective";
		ll = s.length();
		ll = 61 - ll;
		s = s + EGS4.format("", ll);
		EGS4.seqStr = s + EGS4.format(STATLM, 6, true) + " %";
		if (EGS4.iprint > 1)
			printSequence(EGS4.seqStr);

		if (SOURCE == 1) {
			EGS4.seqStr = " Initial RNG state:";
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);
			EGS4.SHOW_RNG_STATE();
		}
		if (RUSROU) {
			// WRITE(IOUT,265)RRZ,RRCUT;
			s = " RUS ROU FOR PHOTONS CROSSING Z = ";
			ll = s.length();
			ll = 56 - ll;
			s = s + EGS4.format("", ll);
			EGS4.seqStr = s + EGS4.format(RRZ, 10, true) + " cm";
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);
			s = " WITH PROBABILITY OF SURVIVAL:";
			ll = s.length();
			ll = 59 - ll;
			s = s + EGS4.format("", ll);
			EGS4.seqStr = s + EGS4.format(RRCUT, 7, true);
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);
		}

		if (EGS4Macro.CEXPTR != 0) {
			s = " PATHLENGTH EXPONENTIAL TRANSFORMATION";
			EGS4.seqStr = s;
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);
			s = "  VARIABLE FOR FORWARD GOING PHOTNS: ";
			ll = s.length();
			ll = 56 - ll;
			s = s + EGS4.format("", ll);
			EGS4.seqStr = s + EGS4.format(EGS4Macro.CEXPTR, 10, true);
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);
		}

		if (IFANO == 1) {
			s = "  *** REGENERATION REQUESTED (IFANO SET TO 1) ! *** ";
			EGS4.seqStr = s;
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);
		} else if (IFANO == 2) {
			s = "  *** ELECTRONS SET IN MOTION IN WALL WILL BE ELIMINATED (IFANO SET TO 2) ! *** ";
			EGS4.seqStr = s;
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);
		} else {
			s = "  *** NO REGENERATION REQUESTED (IFANO SET TO 0) ! *** ";
			EGS4.seqStr = s;
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);

		}

		// "EK0=EIN;"
		// "$PRESTA-INPUT-SUMMARY; OUTPUT THE PRESTA INPUT VARIABLES"
		// "taken out input-summary at upgrade to PRESTA-II"

		// "MATERIAL INPUT SUMMARY"
		// "====================="
		EGS4.seqStr = "=========================================================================";
		if (EGS4.iprint > 1)
			printSequence(EGS4.seqStr);
		EGS4.seqStr = "                   MATERIAL SUMMARY  " + EGS4.NMED
				+ " MATERIALS USED";
		if (EGS4.iprint > 1)
			printSequence(EGS4.seqStr);
		EGS4.seqStr = "=========================================================================";
		if (EGS4.iprint > 1)
			printSequence(EGS4.seqStr);
		EGS4.seqStr = " # MATERIAL  DENSITY(g/cm**3)" + EGS4.format("", 6)
				+ "AE(MeV)" + EGS4.format("", 4) + "AP(MeV)"
				+ EGS4.format("", 9) + "UE(MeV)" + EGS4.format("", 4)
				+ "UP(MeV)";
		if (EGS4.iprint > 1)
			printSequence(EGS4.seqStr);
		EGS4.seqStr = "  - --------  ----------------" + EGS4.format("", 6)
				+ "-------" + EGS4.format("", 4) + "-------"
				+ EGS4.format("", 9) + "-------" + EGS4.format("", 4)
				+ "-------";
		if (EGS4.iprint > 1)
			printSequence(EGS4.seqStr);
		for (int I = 1; I <= EGS4.NMED; I++) {
			// WRITE(IOUT,310)
			// I,(MEDIA(J,I),J=1,6),RHO(I),AE(I),AP(I),UE(I),UP(I);
			// 310 FORMAT(' ',I1,3X,6A1,4X,1PE10.3,2(7X,0PF9.3,2X,F9.3));
			String meds = "";
			if (EGS4.MEDIA[I - 1].length() > 6) {
				meds = EGS4.MEDIA[I - 1].substring(0, 6);
			} else {
				meds = EGS4.MEDIA[I - 1];
			}
			EGS4.seqStr = "  " + EGS4.format(I,3) + EGS4.format("", 3) + meds
					+ EGS4.format("", 4)
					+ EGS4.format(EGS4.RHO[I - 1], 10, false)
					+ EGS4.format("", 6) + EGS4.format(EGS4.AE[I - 1], 10)
					+ EGS4.format("", 1) + EGS4.format(EGS4.AP[I - 1], 10)
					+ EGS4.format("", 6) + EGS4.format(EGS4.UE[I - 1], 10)
					+ EGS4.format("", 1) + EGS4.format(EGS4.UP[I - 1], 10);
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);
		}
		EGS4.iprint = 0;
		EGS4Geom.GEOMRZ_ISUMRY();
		EGS4.iprint = IPRINT;
		if (SOURCE == -1)
			EGS4SrcEns.SRCOUT();

		// "       PRINT A GRID OF THE ZONE DEPENDENT VARIABLES"
		// "       ============================================"

		// "Set defaults for non-used variables"
		// CDSTBL(1)='0';CTRTBL(1)='0';CABSRB(1)='0';
		EGS4Grid.CDSTBL[0] = " ";
		EGS4Grid.CTRTBL[0] = "0";
		EGS4Grid.CABSRB[0] = "0";
		// "For identifying a cavity region"
		// REPLACE {$IRL} WITH {IZ+1+NZ*(IX-1)}
		for (int IRL = 2; IRL <= EGS4Geom.NR * EGS4Geom.NZ + 1; IRL++) {
			if (EGS4Geom.ntrack[IRL - 1] == 1) {
				// EGS4Grid.CAVTRACK[IRL-1]="C";
				EGS4Grid.CDSTBL[IRL - 1] = "C";
			} else {
				// EGS4Grid.CAVTRACK[IRL-1]=" ";
				EGS4Grid.CDSTBL[IRL - 1] = " ";
			}
		}
		// "Make the material grid"
		EGS4Grid.MATERIALGRID(EGS4Geom.NR, EGS4Geom.NZ, AMASS, 1, EGS4.ECUT,
				EGS4.PCUT, EGS4Geom.RCYL, EGS4Geom.ZPLANE, EGS4.MED, EGS4.MEDIA);
		// ,EGS4Grid.CAVTRACK, EGS4Grid.CTRTBL, EGS4Grid.CABSRB);

	}

	/**
	 * Called by AUSGAB. Exchange particles P1 and P2.
	 * @param P1 P1
	 * @param P2 P2
	 */
	private void EXCHANGE_STACK(int P1, int P2) {
		double FDUMMY = 0.0;
		int IDUMMY = 0;
		FDUMMY = EGS4.U[P2 - 1];
		EGS4.U[P2 - 1] = EGS4.U[P1 - 1];
		EGS4.U[P1 - 1] = FDUMMY;
		FDUMMY = EGS4.V[P2 - 1];
		EGS4.V[P2 - 1] = EGS4.V[P1 - 1];
		EGS4.V[P1 - 1] = FDUMMY;
		FDUMMY = EGS4.W[P2 - 1];
		EGS4.W[P2 - 1] = EGS4.W[P1 - 1];
		EGS4.W[P1 - 1] = FDUMMY;
		FDUMMY = EGS4.E[P2 - 1];
		EGS4.E[P2 - 1] = EGS4.E[P1 - 1];
		EGS4.E[P1 - 1] = FDUMMY;
		FDUMMY = EGS4.WT[P2 - 1];
		EGS4.WT[P2 - 1] = EGS4.WT[P1 - 1];
		EGS4.WT[P1 - 1] = FDUMMY;
		IDUMMY = EGS4.IQ[P2 - 1];
		EGS4.IQ[P2 - 1] = EGS4.IQ[P1 - 1];
		EGS4.IQ[P1 - 1] = IDUMMY;
		// "LATCH IS NOW STANDARD"
		IDUMMY = EGS4.LATCH[P2 - 1];
		EGS4.LATCH[P2 - 1] = EGS4.LATCH[P1 - 1];
		EGS4.LATCH[P1 - 1] = IDUMMY;
	}

	/**
	 * Print output summary
	 */
	private void OSUMRY() {
		String s = "";
		int ll = 0;

		// $IMPLICIT-NONE;

		// COMIN/GEOM,IODAT1,IODAT2,PRINTC,SCORE,SOURCE,USER/;
		// COMIN/CH-Steps/;

		int I = 0;
		int IRL = 0;
		int IX = 0;
		int IZ = 0;
		double ASCT = 0.0;
		double ASCTUN = 0.0;
		double AATT = 0.0;
		double AATTUN = 0.0;
		double AWLL = 0.0;
		double AWLLUN = 0.0;
		double TDAW = 0.0;
		double TDAWUN = 0.0;
		double KSCT = 0.0;
		double KATT = 0.0;
		double KWLL = 0.0;
		double KSCTUN = 0.0;
		double KATTUN = 0.0;
		double KWLLUN = 0.0;

		// "SET UP THE PRINTER"
		// ICHPIN=12; "12 CHARACTERS/INCH"
		// ILPIN=6; "6 LINES/INCH"
		// IPAGE=0; "NO PAGE THROW"
		// "CALL PRNTER(ICHPIN,ILPIN,IOUT,IPAGE);"

		// "WRITE(IOUT,100)TITLE,DATEN,TIMEN; HEADER"

		// IF(ISOURC=21 | ISOURC=22)[
		// WRITE(IOUT,200) SCSTP,SCSTP2,
		// SCSTP/(dble(NNREAD+NRCYCL*(NNREAD-IHSTRY))/dble(NCASE_PHSP)*
		// NINCSRC),SCSTP2,(count_pII_steps+PIISTP)/SCSTP,SCSTP2;
		// ]
		// ELSE[

		s = "                    # primary charged particle steps";
		ll = s.length();
		ll = 58 - ll;
		s = s + EGS4.format("", ll);

		EGS4.seqStr = s + EGS4.format(SCSTP, 10, false) + " +/- "
				+ EGS4.format(SCSTP2, 6) + "%";
		if (EGS4.iprint > 0)
			printSequence(EGS4.seqStr);

		s = "     # primary charged particle steps/initial history";
		ll = s.length();
		ll = 58 - ll;
		s = s + EGS4.format("", ll);

		EGS4.seqStr = s
				+ EGS4.format(SCSTP / EGS4SrcEns.dble(IHSTRY), 10, false)
				+ " +/- " + EGS4.format(SCSTP2, 6) + "%";
		if (EGS4.iprint > 0)
			printSequence(EGS4.seqStr);

		s = "# of presta-II steps/# primary charged particle steps";
		ll = s.length();
		ll = 58 - ll;
		s = s + EGS4.format("", ll);

		EGS4.seqStr = s
				+ EGS4.format((EGS4.count_pII_steps + PIISTP) / SCSTP, 10,
						false) + " +/- " + EGS4.format(SCSTP2, 6) + "%";
		if (EGS4.iprint > 0)
			printSequence(EGS4.seqStr);

		// WRITE(IOUT,200) SCSTP,SCSTP2,SCSTP/dble(IHSTRY),
		// SCSTP2,(count_pII_steps+PIISTP)/SCSTP,SCSTP2;
		// ]
		// 200 FORMAT(//' ' ,' # primary charged particle steps',T58,
		// 1PE10.3,' +/- ',0PF6.3,'%'/
		// ' ',' # primary charged particle steps/initial history',T58,
		// 1PE10.3,' +/- ',0PF6.3,'%'/
		// ' ','# of presta-II steps/# primary charged particle steps',
		// T58,F10.3,' +/- ',0PF6.3,'%');

		// "PRINT # CHARGED PARTICLE STEPS IN cavity REGION, etc"
		// IF(ISOURC=21 | ISOURC=22)[
		// WRITE(IOUT,210) SCCSTP,SCCSTP2,
		// SCCSTP/(dble(NNREAD+NRCYCL*(NNREAD-IHSTRY))/dble(NCASE_PHSP)*
		// NINCSRC),SCCSTP2;
		// ]
		// ELSE[
		s = "   # primary charged particle steps in cavity region";
		ll = s.length();
		ll = 58 - ll;
		s = s + EGS4.format("", ll);

		EGS4.seqStr = s + EGS4.format(SCCSTP, 10, false) + " +/- "
				+ EGS4.format(SCCSTP2, 6) + "%";
		if (EGS4.iprint > 0)
			printSequence(EGS4.seqStr);

		s = "    # primary steps in cavity region/initial history";
		ll = s.length();
		ll = 58 - ll;
		s = s + EGS4.format("", ll);

		EGS4.seqStr = s
				+ EGS4.format(SCCSTP / EGS4SrcEns.dble(IHSTRY), 10, false)
				+ " +/- " + EGS4.format(SCCSTP2, 6) + "%";
		if (EGS4.iprint > 0)
			printSequence(EGS4.seqStr);

		// WRITE(IOUT,210) SCCSTP,SCCSTP2,SCCSTP/dble(IHSTRY),SCCSTP2;
		// ]
		// 210 FORMAT(//' ',' # primary charged particle steps in cavity region'
		// ,T58,1PE10.3,' +/- ',0PF6.3,'%'/
		// ' ',' # primary steps in cavity region/initial history'
		// ,T58,1PE10.3,' +/- ',0PF6.3,'%');

		if (EGS4SrcEns.ISOURC == 15)
			EGS4SrcEns.src15_out();// (iout);

		s = "******************************************************************************************";
		EGS4.seqStr = s;
		if (EGS4.iprint > 0)
			printSequence(EGS4.seqStr);

		s = "  Incident fluence (from source): ";
		ll = s.length();
		ll = 58 - ll;
		s = s + EGS4.format("", ll);

		EGS4.seqStr = s + EGS4.format(EGS4SrcEns.AINFLU, 10, false)
				+ " 1/cm^2 ";
		if (EGS4.iprint > 0)
			printSequence(EGS4.seqStr);

		s = "******************************************************************************************";
		EGS4.seqStr = s;
		if (EGS4.iprint > 0)
			printSequence(EGS4.seqStr);

		// "THE CAVITY SUMMARY"
		// "******************"

		if (EGS4Geom.NSUMCV == 1) {
			EGS4.seqStr = "                   SUM OF RESULTS FOR THE CAVITY: 1 REGION";
			if (EGS4.iprint > 0)
				printSequence(EGS4.seqStr);
			EGS4.seqStr = "                   ***************************************";
			if (EGS4.iprint > 0)
				printSequence(EGS4.seqStr);
			// WRITE(IOUT,300);
		} else {
			EGS4.seqStr = "                   SUM OF RESULTS FOR THE CAVITY: "
					+ EGS4Geom.NSUMCV + " REGION";
			if (EGS4.iprint > 0)
				printSequence(EGS4.seqStr);
			EGS4.seqStr = "                   ***************************************";
			if (EGS4.iprint > 0)
				printSequence(EGS4.seqStr);

			// WRITE(IOUT,301)NSUMCV;
		}

		// 300 FORMAT(// ,T20,'SUM OF RESULTS FOR THE CAVITY: 1 REGION'/
		// ' ',T20,'***************************************');
		// 301 FORMAT(// ,T20,'SUM OF RESULTS FOR THE CAVITY: ',I2,' REGIONS'/
		// ' ',T20,'*****************************************');
		// "******************"
		// "history by history"
		// "    EMH March 2002"
		// "******************"
		cav2_dose = cav2_dose / cav_dose;
		// edep2=edep2/edep;
		// System.out.println("@@ "+EGS4.format(edep,11,false)+" +/- "+EGS4.format(100*edep2,6)+"%");
		if (EGS4Macro.iifano == 1) {
			EGS4.seqStr = " This calculation was performed using regeneration ";
			if (EGS4.iprint > 0)
				printSequence(EGS4.seqStr);
			EGS4.seqStr = " ================================================= ";
			if (EGS4.iprint > 0)
				printSequence(EGS4.seqStr);
			EGS4.seqStr = " D=Total dose; Awall=D/Dose with attenuation and scatter removed ";
			if (EGS4.iprint > 0)
				printSequence(EGS4.seqStr);
			s = " D/Awall (grays/incident fluence): ";
			ll = s.length();
			ll = 50 - ll;
			s = s + EGS4.format("", ll);
			EGS4.seqStr = s + EGS4.format(cav_dose, 11, false) + " +/- "
					+ EGS4.format(100 * cav2_dose, 6) + "%";
			if (EGS4.iprint > 0)
				printSequence(EGS4.seqStr);

			// write(iout,'(a,t50,1PE11.4,a,0PF6.3,a)')
			// 'D/Awall (grays/incident fluence): ',cav_dose,' +/- ',
			// 100*cav2_dose,'%';
			return;
		}

		if (EGS4Macro.use_enhance) {
			EGS4.seqStr = " This calculation was performed using CS enhancement ";
			if (EGS4.iprint > 0)
				printSequence(EGS4.seqStr);
			EGS4.seqStr = "  enhancement factor was "
					+ EGS4.format(EGS4Macro.cs_enhance, 10, true);
			if (EGS4.iprint > 0)
				printSequence(EGS4.seqStr);
			EGS4.seqStr = " ================================================= ";
			if (EGS4.iprint > 0)
				printSequence(EGS4.seqStr);
		} else if (EGS4Macro.n_split > 1) {
			EGS4.seqStr = " This calculation was performed using photon splitting ";
			if (EGS4.iprint > 0)
				printSequence(EGS4.seqStr);
			EGS4.seqStr = "  splitting number was "
					+ EGS4.format(EGS4Macro.n_split, 6);
			if (EGS4.iprint > 0)
				printSequence(EGS4.seqStr);
			EGS4.seqStr = " ================================================= ";
			if (EGS4.iprint > 0)
				printSequence(EGS4.seqStr);
		}
		if (IFANO == 1) {
			EGS4.seqStr = " This calculation was performed using regeneration ";
			if (EGS4.iprint > 0)
				printSequence(EGS4.seqStr);

		}
		if (IFANO == 2) {
			EGS4.seqStr = " This calculation was performed eliminating electrons originating in the cavity wall.";
			if (EGS4.iprint > 0)
				printSequence(EGS4.seqStr);
		}
		// "THE TOTAL DOSE"
		if (EGS4SrcEns.ISOURC == 3 || EGS4SrcEns.ISOURC == 21
				|| EGS4SrcEns.ISOURC == 22) {
			s = " TOTAL DOSE (GRAYS/INCIDENT PARTICLE):";
			ll = s.length();
			ll = 50 - ll;
			s = s + EGS4.format("", ll);

			EGS4.seqStr = s + EGS4.format(cav_dose, 11, false) + " +/- "
					+ EGS4.format(100 * cav2_dose, 5) + "%";
			if (EGS4.iprint > 0)
				printSequence(EGS4.seqStr);

			// WRITE(IOUT,311)cav_dose,100*cav2_dose;
		} else {
			EGS4.seqStr = "Awall=Total dose/Primary dose with attenuation corrected";
			if (EGS4.iprint > 0)
				printSequence(EGS4.seqStr);
			EGS4.seqStr = "Aatt=Primary dose/Primary dose with attenuation corrected";
			if (EGS4.iprint > 0)
				printSequence(EGS4.seqStr);
			EGS4.seqStr = "Ascat=(Primary dose+Secondary dose)/Primary dose";
			if (EGS4.iprint > 0)
				printSequence(EGS4.seqStr);
			EGS4.seqStr = "Kwall=1/Awall";
			if (EGS4.iprint > 0)
				printSequence(EGS4.seqStr);
			EGS4.seqStr = "Katt=1/Aatt";
			if (EGS4.iprint > 0)
				printSequence(EGS4.seqStr);
			EGS4.seqStr = "Kscat=1/Ascat";
			if (EGS4.iprint > 0)
				printSequence(EGS4.seqStr);

			s = " TOTAL DOSE (GRAYS/INCIDENT FLUENCE):";
			ll = s.length();
			ll = 50 - ll;
			s = s + EGS4.format("", ll);

			EGS4.seqStr = s + EGS4.format(cav_dose, 11, false) + " +/- "
					+ EGS4.format(100 * cav2_dose, 5) + "%";
			if (EGS4.iprint > 0)
				printSequence(EGS4.seqStr);

			// WRITE(IOUT,310)cav_dose,100*cav2_dose;
		}
		// 310 FORMAT(//'TOTAL DOSE (GRAYS/INCIDENT FLUENCE):',
		// T50,1PE11.4,' +/- ',0PF5.2,'%');
		// 311 FORMAT(//'TOTAL DOSE (GRAYS/INCIDENT PARTICLE):',
		// T50,1PE11.4,' +/- ',0PF5.2,'%');

		if (EGS4Macro.IFULL == 1) {
			// "CALCULATE Ascat, Aatt, Awall, DOSE/Awall"

			cav2_dose1 = cav2_dose1 / cav_dose1;
			cav_dosec = cav_dosec / cav_dose / cav_dose1;

			// "total dose/Awall"
			TDAW = cav_dose1;
			TDAWUN = cav2_dose1;

			s = " TOTAL DOSE/Awall:";
			ll = s.length();
			ll = 50 - ll;
			s = s + EGS4.format("", ll);
			EGS4.seqStr = s + EGS4.format(TDAW, 11, false) + " +/- "
					+ EGS4.format(100 * TDAWUN, 5) + "%";
			if (EGS4.iprint > 0)
				printSequence(EGS4.seqStr);

			// WRITE(IOUT,350)TDAW,TDAWUN*100.;
			// 350 FORMAT(' TOTAL DOSE/Awall:',
			// T50,1PE11.4,' +/- ',0PF5.2,'%');

			cav_dosec = cav2_dose * cav2_dose + cav2_dose1 * cav2_dose1 - 2
					* cav_dosec;
			if (cav_dosec > 0)
				cav_dosec = Math.sqrt(cav_dosec);
			// "Awall"
			AWLL = cav_dose / cav_dose1;
			AWLLUN = cav_dosec;
			// WRITE(IOUT,340)AWLL,AWLLUN*100.;
			// 340 FORMAT(' Awall:',T50,F8.5,' +/- ',F5.2,'%');
			s = " Awall:";
			ll = s.length();
			ll = 50 - ll;
			s = s + EGS4.format("", ll);
			EGS4.seqStr = s + EGS4.format(AWLL, 8, true) + "  +/- "
					+ EGS4.format(100 * AWLLUN, 5, true) + "%";
			if (EGS4.iprint > 0)
				printSequence(EGS4.seqStr);

			if (!EGS4Macro.use_enhance && EGS4Macro.n_split <= 1) {
				corr_02 = corr_02 / (cav_dose0 + cav_dose2);

				cav2_dose0 = cav2_dose0 / cav_dose0;
				cav2_dose2 = cav2_dose2 / cav_dose2;

				cav_dosec01 = cav_dosec01 / cav_dose0 / cav_dose1;
				cav_dosec02 = cav_dosec02 / cav_dose0 / cav_dose2;

				cav_dosec01 = cav2_dose0 * cav2_dose0 + cav2_dose1 * cav2_dose1
						- 2 * cav_dosec01;
				cav_dosec02 = cav2_dose0 * cav2_dose0 + cav2_dose2 * cav2_dose2
						- 2 * cav_dosec02;

				if (cav_dosec01 > 0)
					cav_dosec01 = Math.sqrt(cav_dosec01);
				if (cav_dosec02 > 0)
					cav_dosec02 = Math.sqrt(cav_dosec02);

				// "... we had rel.error(x), get now rel.error(1+x)"
				cav_dosec02 = cav_dosec02 / (1 + cav_dose0 / cav_dose2);

				// "Ascatt"
				ASCT = 1 + cav_dose2 / cav_dose0;
				ASCTUN = cav_dosec02;
				// WRITE(IOUT,320)ASCT,ASCTUN*100.;
				// / 320 FORMAT(' Ascat:',T50,F8.5,' +/- ',F5.2,'%');
				s = " Ascat:";
				ll = s.length();
				ll = 50 - ll;
				s = s + EGS4.format("", ll);
				EGS4.seqStr = s + EGS4.format(ASCT, 8, true) + "  +/- "
						+ EGS4.format(100 * ASCTUN, 5, true) + "%";
				if (EGS4.iprint > 0)
					printSequence(EGS4.seqStr);

				// "Aatt"
				AATT = cav_dose0 / cav_dose1;
				AATTUN = cav_dosec01;
				// WRITE(IOUT,330)AATT,AATTUN*100.;
				// 330 FORMAT(' Aatt :',T50,F8.5,' +/- ',F5.2,'%');
				s = " Aatt :";
				ll = s.length();
				ll = 50 - ll;
				s = s + EGS4.format("", ll);
				EGS4.seqStr = s + EGS4.format(AATT, 8, true) + "  +/- "
						+ EGS4.format(100 * AATTUN, 5, true) + "%";
				if (EGS4.iprint > 0)
					printSequence(EGS4.seqStr);

				// "Dpr + Dsec"
				// 345 FORMAT(' Dprimary + Dsecondary:',T50,1PE11.4,' +/-
				// ',0PF5.2,'%');
				// WRITE(IOUT,345)cav_dose2+cav_dose0,corr_02*100.;
				s = " Dprimary + Dsecondary:";
				ll = s.length();
				ll = 50 - ll;
				s = s + EGS4.format("", ll);
				EGS4.seqStr = s + EGS4.format(cav_dose2 + cav_dose0, 11, false)
						+ "  +/- " + EGS4.format(corr_02 * 100., 5, true) + "%";
				if (EGS4.iprint > 0)
					printSequence(EGS4.seqStr);

			}
		} // "END OF IFULL = 1"

		// "OUTPUT Kscat,Katt,Kwall,Kpn,Kfl THE INVERSES OF THE Ai's"
		if (EGS4Macro.IFULL > 0) {
			KWLL = 1. / AWLL;
			KWLLUN = AWLLUN / (AWLL * AWLL);
			s = " Kwall:";
			ll = s.length();
			ll = 50 - ll;
			s = s + EGS4.format("", ll);
			EGS4.seqStr = s + EGS4.format(KWLL, 8, true) + "  +/- "
					+ EGS4.format(100 * KWLLUN, 5, true) + "%";
			if (EGS4.iprint > 0)
				printSequence(EGS4.seqStr);

			// WRITE(IOUT,341)KWLL,KWLLUN*100.;
			// 341 FORMAT(' Kwall:',T50,F8.5,' +/- ',F5.2,'%');
			if (!EGS4Macro.use_enhance && EGS4Macro.n_split <= 1) {
				KSCT = 1. / ASCT;
				KSCTUN = ASCTUN / (ASCT * ASCT);
				s = " Kscat:";
				ll = s.length();
				ll = 50 - ll;
				s = s + EGS4.format("", ll);
				EGS4.seqStr = s + EGS4.format(KSCT, 8, true) + "  +/- "
						+ EGS4.format(100 * KSCTUN, 5, true) + "%";
				if (EGS4.iprint > 0)
					printSequence(EGS4.seqStr);

				// WRITE(IOUT,321)KSCT,KSCTUN*100.;
				// / 321 FORMAT(' Kscat:',T50,F8.5,' +/- ',F5.2,'%');
				KATT = 1. / AATT;
				KATTUN = AATTUN / (AATT * AATT);
				// WRITE(IOUT,331)KATT,KATTUN*100.;
				// 331 FORMAT(' Katt :',T50,F8.5,' +/- ',F5.2,'%');
				s = " Katt :";
				ll = s.length();
				ll = 50 - ll;
				s = s + EGS4.format("", ll);
				EGS4.seqStr = s + EGS4.format(KATT, 8, true) + "  +/- "
						+ EGS4.format(100 * KATTUN, 5, true) + "%";
				if (EGS4.iprint > 0)
					printSequence(EGS4.seqStr);

			}
		}

		// "THE DETAILED OUTPUT FOR EACH CAVITY ZONE"
		// "****************************************"
		if (IOOPTN == 1 && EGS4Geom.NSUMCV > 1) {
			// "ONLY IF REQUESTED AND MORE THAN ONE CAVITY ZONE"

			// WRITE(IOUT,400) NSUMCV;
			// 400 FORMAT(// ,T15,'DETAILED RESULTS FOR EACH OF THE ',I4,'
			// CAVITY REGIONS'/
			// /
			// ' ',T15,'****************************************************');
			EGS4.seqStr = "   DETAILED RESULTS FOR EACH OF THE "
					+ EGS4.format(EGS4Geom.NSUMCV, 4) + " CAVITY REGIONS";
			if (EGS4.iprint > 0)
				printSequence(EGS4.seqStr);

			// "PRINT THE TABLE HEADER"
			if (EGS4Macro.IFULL == 0) {
				// WRITE(IOUT,410);
				// 410 FORMAT(//'Z# P# C# Total Dose '/
				// / ' -- -- -- -------------------');
				EGS4.seqStr = " Z# P# C#     Total Dose     ";
				if (EGS4.iprint > 0)
					printSequence(EGS4.seqStr);
				EGS4.seqStr = " -- -- -- -------------------";
				if (EGS4.iprint > 0)
					printSequence(EGS4.seqStr);

			} else if (EGS4Macro.IFULL > 0) {
				// "SET UP THE PRINTER"
				// ICHPIN=16; "16 CHARACTERS/INCH"
				// ILPIN=6; "6 LINES/INCH"
				// IPAGE=0; "NO PAGE THROW"
				// "CALL PRNTER(ICHPIN,ILPIN,IOUT,IPAGE);"

				// WRITE(IOUT,420);
				// / 420 FORMAT(//'Z# P# C#',
				// ' Total Dose ',
				// ' Ascat ',
				// ' Aatt ',
				// ' Awall ',
				// ' Total Dose/Awall'/
				// ' -- -- --',
				// ' ---------- ',
				// ' ----- ',
				// ' ---- ',
				// ' ----- ',
				// ' ----------------');
				EGS4.seqStr = " Z# P# C#" + "     Total Dose   "
						+ "      Ascat      " + "       Aatt      "
						+ "      Awall      " + "  Total Dose/Awall";
				if (EGS4.iprint > 0)
					printSequence(EGS4.seqStr);
				EGS4.seqStr = " -- -- --" + "     ----------   "
						+ "      -----      " + "       ----      "
						+ "      -----      " + "  ----------------";
				if (EGS4.iprint > 0)
					printSequence(EGS4.seqStr);

			}

			// "LOOP OVER THE CAVITY ZONES"
			for (I = 1; I <= EGS4Geom.NSUMCV; I++) {
				IRL = EGS4Geom.ISUMCV[I - 1];
				// $GET-IX-IZ(IRL);
				IX = EGS4Geom.GET_IX(IRL);
				IZ = EGS4Geom.GET_IZC(IRL);

				// "total dose"
				if (SCDOSE[IZ - 1][IX - 1][0] > 0.)
					SCDOSE2[IZ - 1][IX - 1][0] = SCDOSE2[IZ - 1][IX - 1][0]
							/ SCDOSE[IZ - 1][IX - 1][0];
				if (EGS4Macro.IFULL == 1) {
					// "CALCULATE Ascat, Aatt, Awall, DOSE/Awall"

					if (SCDOSE[IZ - 1][IX - 1][2] > 0.) {
						SCDOSE2[IZ - 1][IX - 1][2] = SCDOSE2[IZ - 1][IX - 1][2]
								/ SCDOSE[IZ - 1][IX - 1][2];
						SCDOSE_COV[IZ - 1][IX - 1][0] = SCDOSE_COV[IZ - 1][IX - 1][0]
								/ SCDOSE[IZ - 1][IX - 1][0]
								/ SCDOSE[IZ - 1][IX - 1][2];
					}

					// "TOTAL DOSE/Awall"
					TDAW = SCDOSE[IZ - 1][IX - 1][2];
					TDAWUN = SCDOSE2[IZ - 1][IX - 1][2];

					SCDOSE_COV[IZ - 1][IX - 1][0] = SCDOSE2[IZ - 1][IX - 1][0]
							* SCDOSE2[IZ - 1][IX - 1][0]
							+ SCDOSE2[IZ - 1][IX - 1][2]
							* SCDOSE2[IZ - 1][IX - 1][2] - 2
							* SCDOSE_COV[IZ - 1][IX - 1][0];
					if (SCDOSE_COV[IZ - 1][IX - 1][0] > 0.)
						SCDOSE_COV[IZ - 1][IX - 1][0] = Math
								.sqrt(SCDOSE_COV[IZ - 1][IX - 1][0]);
					// "Awall"
					AWLL = SCDOSE[IZ - 1][IX - 1][0]
							/ SCDOSE[IZ - 1][IX - 1][2];
					AWLLUN = SCDOSE_COV[IZ - 1][IX - 1][0];

					if (!EGS4Macro.use_enhance && EGS4Macro.n_split <= 1) { // "this condition not strictly necessary because"
																			// "IOOPTN set = 0 if either of these options on"

						if (SCDOSE[IZ - 1][IX - 1][1] > 0.) {
							SCDOSE2[IZ - 1][IX - 1][1] = SCDOSE2[IZ - 1][IX - 1][1]
									/ SCDOSE[IZ - 1][IX - 1][1];
							SCDOSE_COV[IZ - 1][IX - 1][1] = SCDOSE_COV[IZ - 1][IX - 1][1]
									/ SCDOSE[IZ - 1][IX - 1][1]
									/ SCDOSE[IZ - 1][IX - 1][2];
							if (SCDOSE[IZ - 1][IX - 1][3] > 0.) {
								SCDOSE_COV[IZ - 1][IX - 1][2] = SCDOSE_COV[IZ - 1][IX - 1][2]
										/ SCDOSE[IZ - 1][IX - 1][1]
										/ SCDOSE[IZ - 1][IX - 1][3];
							}
						}
						if (SCDOSE[IZ - 1][IX - 1][3] > 0.)
							SCDOSE2[IZ - 1][IX - 1][3] = SCDOSE2[IZ - 1][IX - 1][3]
									/ SCDOSE[IZ - 1][IX - 1][3];

						SCDOSE_COV[IZ - 1][IX - 1][1] = SCDOSE2[IZ - 1][IX - 1][1]
								* SCDOSE2[IZ - 1][IX - 1][1]
								+ SCDOSE2[IZ - 1][IX - 1][2]
								* SCDOSE2[IZ - 1][IX - 1][2]
								- 2
								* SCDOSE_COV[IZ - 1][IX - 1][1];
						SCDOSE_COV[IZ - 1][IX - 1][2] = SCDOSE2[IZ - 1][IX - 1][1]
								* SCDOSE2[IZ - 1][IX - 1][1]
								+ SCDOSE2[IZ - 1][IX - 1][3]
								* SCDOSE2[IZ - 1][IX - 1][3]
								- 2
								* SCDOSE_COV[IZ - 1][IX - 1][2];

						if (SCDOSE_COV[IZ - 1][IX - 1][1] > 0)
							SCDOSE_COV[IZ - 1][IX - 1][1] = Math
									.sqrt(SCDOSE_COV[IZ - 1][IX - 1][1]);
						if (SCDOSE_COV[IZ - 1][IX - 1][2] > 0)
							SCDOSE_COV[IZ - 1][IX - 1][2] = Math
									.sqrt(SCDOSE_COV[IZ - 1][IX - 1][2]);

						// "... we had rel.error(x), get now rel.error(1+x)"
						if (SCDOSE[IZ - 1][IX - 1][3] > 0.)
							SCDOSE_COV[IZ - 1][IX - 1][2] = SCDOSE_COV[IZ - 1][IX - 1][2]
									/ (1 + SCDOSE[IZ - 1][IX - 1][1]
											/ SCDOSE[IZ - 1][IX - 1][3]);

						// "Ascatt"
						ASCT = 1 + SCDOSE[IZ - 1][IX - 1][3]
								/ SCDOSE[IZ - 1][IX - 1][1];
						ASCTUN = SCDOSE_COV[IZ - 1][IX - 1][2];

						// "Aatt"
						AATT = SCDOSE[IZ - 1][IX - 1][1]
								/ SCDOSE[IZ - 1][IX - 1][2];
						AATTUN = SCDOSE_COV[IZ - 1][IX - 1][1];

					}
				}
				if (EGS4Macro.IFULL == 0) {
					EGS4.seqStr = " " + EGS4.format(IRL, 2) + " "
							+ EGS4.format(IZ, 2) + " " + EGS4.format(IX, 2)
							+ " "
							+ EGS4.format(SCDOSE[IZ - 1][IX - 1][0], 11, false)
							+ " +/-"
							+ EGS4.format(100 * SCDOSE2[IZ - 1][IX - 1][0], 5)
							+ "%";
					if (EGS4.iprint > 0)
						printSequence(EGS4.seqStr);

					// WRITE(IOUT,430)
					// IRL,IZ,IX, "POSITION"
					// SCDOSE(IZ,IX,1),SCDOSE2(IZ,IX,1)*100.; "TOTAL DOSE"
				}

				// 430 FORMAT(' ',I2,2(1X,I2),1X,1PE11.4,' +/-',0PF5.2,'%');
				// 440 FORMAT(' ',I2,2(1X,I2),
				// 1X,1PE11.4,'(',0PF5.2,'%)',
				// 3(F8.5,'(',F7.5,')'),
				// 1X,1PE11.4,'(',0PF5.2,'%)');

				else if (EGS4Macro.IFULL == 1) {
					EGS4.seqStr = " " + EGS4.format(IRL, 2) + " "
							+ EGS4.format(IZ, 2) + " " + EGS4.format(IX, 2)
							+ " "
							+ EGS4.format(SCDOSE[IZ - 1][IX - 1][0], 11, false)
							+ "("
							+ EGS4.format(100 * SCDOSE2[IZ - 1][IX - 1][0], 5)
							+ "%)" + EGS4.format(ASCT, 8, true) + "("
							+ EGS4.format(ASCTUN, 7, true) + ")"
							+ EGS4.format(AATT, 8, true) + "("
							+ EGS4.format(AATTUN, 7, true) + ")"
							+ EGS4.format(AWLL, 8, true) + "("
							+ EGS4.format(AWLLUN, 7, true) + ")"
							+ EGS4.format(TDAW, 11, false) + "("
							+ EGS4.format(100 * TDAWUN, 5) + "%)";
					if (EGS4.iprint > 0)
						printSequence(EGS4.seqStr);

					// WRITE(IOUT,440)
					// IRL,IZ,IX, "POSITION"
					// SCDOSE(IZ,IX,1),SCDOSE2(IZ,IX,1)*100., "TOTAL DOSE"
					// ASCT,ASCTUN, "Ascat"
					// AATT,AATTUN, "Aatt"
					// AWLL,AWLLUN, "Awall"
					// TDAW,TDAWUN*100; "TOTAL DOSE/Awall"
				}

			}// "END OF LOOP OVER CAVITY ZONES"
		}// "END OF DETAILED SUMMARY"

		// "RESET UP THE PRINTER"
		// ICHPIN=12; "12 CHARACTERS/INCH"
		// ILPIN=6; "6 LINES/INCH"
		// IPAGE=0; "NO PAGE THROW"
		// "CALL PRNTER(ICHPIN,ILPIN,IOUT,IPAGE);"

		// RETURN;
		// %I0
		// "FORMATS"
		// 200 FORMAT(//' ' ,' # primary charged particle steps',T58,
		// 1PE10.3,' +/- ',0PF6.3,'%'/
		// / ' ',' # primary charged particle steps/initial history',T58,
		// 1PE10.3,' +/- ',0PF6.3,'%'/
		// / ' ','# of presta-II steps/# primary charged particle steps',
		// // T58,F10.3,' +/- ',0PF6.3,'%');
		// / 202 FORMAT(//' ' ,' # charged particle steps in run',T58,
		// / 1PE10.3,/
		// / ' ',' # charged particle steps in run/initial history',T58,
		// 1PE10.3/
		// ' ','# of presta-II steps/# primary charged particle steps',
		// / T58,F10.3,' +/- ',0PF6.3,'%');
		// 210 FORMAT(//' ',' # primary charged particle steps in cavity region'
		// ,T58,1PE10.3,' +/- ',0PF6.3,'%'/
		// ' ',' # primary steps in cavity region/initial history'
		// ,T58,1PE10.3,' +/- ',0PF6.3,'%');
		// 220 FORMAT(// ,T8,'STEP COUNTING RESULTS FOR WALL MATERIAL IN THE
		// CAVITY'/
		// / // ,T8,'# primary charged particle steps',T51,
		// / I12,' +/- ',0PF5.2,'%'/
		// ' ',T8,'# OF Times mscat switched off',T51,
		// // I12,' +/- ',0PF5.2,'%'/
		// / ' ',T8,'RATIO',T54,F7.3,' +/- ',0PF5.2,'%');
		// 230 FORMAT(// ,T8,'# Primary charged particle steps in cavity region'
		// / ,T51,I12,' +/- ',0PF5.2,'%'/
		// ' ',T8,'# Times mscat switched off in cavity region.'
		// ,T51,I12,' +/- ',0PF5.2,'%'/
		// ' ',T8,'Ratio',T56,F7.3,' +/- ',0PF5.2,'%');
		// / 300 FORMAT(// ,T20,'SUM OF RESULTS FOR THE CAVITY: 1 REGION'/
		// / ' ',T20,'***************************************');
		// 301 FORMAT(// ,T20,'SUM OF RESULTS FOR THE CAVITY: ',I2,' REGIONS'/
		// ' ',T20,'*****************************************');
		// 310 FORMAT(//'TOTAL DOSE (GRAYS/INCIDENT FLUENCE):',
		// T50,1PE11.4,' +/- ',0PF5.2,'%');
		// 311 FORMAT(//'TOTAL DOSE (GRAYS/INCIDENT PARTICLE):',
		// T50,1PE11.4,' +/- ',0PF5.2,'%');
		// / 320 FORMAT(' Ascat:',T50,F8.5,' +/- ',F5.2,'%');
		// / 321 FORMAT(' Kscat:',T50,F8.5,' +/- ',F5.2,'%');
		// 330 FORMAT(' Aatt :',T50,F8.5,' +/- ',F5.2,'%');
		// 331 FORMAT(' Katt :',T50,F8.5,' +/- ',F5.2,'%');
		// 340 FORMAT(' Awall:',T50,F8.5,' +/- ',F5.2,'%');
		// 341 FORMAT(' Kwall:',T50,F8.5,' +/- ',F5.2,'%');
		// 345 FORMAT(' Dprimary + Dsecondary:',T50,1PE11.4,' +/- ',0PF5.2,'%');
		// 350 FORMAT(' TOTAL DOSE/Awall:',
		// T50,1PE11.4,' +/- ',0PF5.2,'%');
		// / 360 FORMAT(//'Apn :',T50,F8.5,' +/- ',F8.5);
		// 361 FORMAT(' Kpn :',T50,F8.5,' +/- ',F8.5);
		// 370 FORMAT(' TOTAL DOSE/([Apn]*Awall):',
		// T50,1PE11.4,' +/- ',0PF5.2,'%');
		// 380 FORMAT(//'Afl :',T50,F8.5,' +/- ',F8.5);
		// 381 FORMAT(' Kfl :',T50,F8.5,' +/- ',F8.5);
		// 390 FORMAT(' TOTAL DOSE/(Afl*[Apn]*Awall):',
		// T50,1PE11.4,' +/- ',0PF5.2,'%');
		// 395 FORMAT(//'<s>g,w :',T50,F8.5,' +/- ',F8.5);
		// / 396 FORMAT(' TOTAL DOSE/(Afl*[Apn]*Awall*<s>g,w):',
		// / T50,1PE11.4,' +/- ',0PF5.2,'%');
		// 400 FORMAT(// ,T15,'DETAILED RESULTS FOR EACH OF THE ',I4,' CAVITY
		// REGIONS'/
		// / ' ',T15,'****************************************************');
		// 410 FORMAT(//'Z# P# C# Total Dose '/
		// / ' -- -- -- -------------------');
		// / 420 FORMAT(//'Z# P# C#',
		// ' Total Dose ',
		// ' Ascat ',
		// ' Aatt ',
		// ' Awall ',
		// ' Total Dose/Awall'/
		// ' -- -- --',
		// ' ---------- ',
		// ' ----- ',
		// ' ---- ',
		// ' ----- ',
		// ' ----------------');
		// 430 FORMAT(' ',I2,2(1X,I2),1X,1PE11.4,' +/-',0PF5.2,'%');
		// 440 FORMAT(' ',I2,2(1X,I2),
		// 1X,1PE11.4,'(',0PF5.2,'%)',
		// 3(F8.5,'(',F7.5,')'),
		// 1X,1PE11.4,'(',0PF5.2,'%)');
		// 445 FORMAT(' ',I2,2(1X,I2),
		// 1X,1PE11.4,'(',0PF5.2,'%)',
		// 3(F8.5,'(',F7.5,')'),
		// 1X,1PE11.4,'(',0PF5.2,'%)',
		// F8.5,'(',F7.5,')',
		// 1X,1PE11.4,'(',0PF5.2,'%)');
		// / 450 FORMAT(' ',I2,2(1X,I2),
		// 1X,1PE11.4,
		// 3(1X,F8.5,1X),
		// 1X,1PE11.4,
		// 1X,F8.5,1X,
		// 1X,1PE11.4,
		// 1X,F8.5,1X,
		// 1X,1PE11.4,
		// 1X,F8.5,1X,
		// / 1X,1PE11.4/
		// // ' ','ERRORS = ',
		// 4X,'(',0PF5.2,'%)',
		// / 3(1X,'(',F7.5,')'),
		// 4X,'(',0PF5.2,'%)',
		// / 1X,'(',F7.5,')',
		// 4X,'(',0PF5.2,'%)',
		// 1X,'(',F7.5,')',
		// / 4X,'(',0PF5.2,'%)',
		// / 1X,'(',F7.5,')',
		// 4X,'(',0PF5.2,'%)');
		//
		// / END; "END OF SUBROUTINE OSUMRY"
		//

	}

	/**
	 * Fix some input variables in special cases. Called by init.
	 */
	private void fixSRCOTO() {
		// srcrz
		EGS4SrcEns.ENFLAG = 0;
		for (int i = 1; i <= EGS4SrcEns.nsource_option; i++) {
			if (EGS4SrcEns.source_option[i - 1] < EGS4SrcEns.source_option_min
					|| EGS4SrcEns.source_option[i - 1] > EGS4SrcEns.source_option_max)
				EGS4SrcEns.source_option[i - 1] = EGS4SrcEns.source_option_default;
		}

		EGS4SrcEns.TEMP1 = EGS4SrcEns.source_option[0];
		EGS4SrcEns.TEMP2 = EGS4SrcEns.source_option[1];
		EGS4SrcEns.TEMP3 = EGS4SrcEns.source_option[2];
		EGS4SrcEns.TEMP4 = EGS4SrcEns.source_option[3];

		EGS4SrcEns.SVTMP1 = EGS4SrcEns.TEMP1;
		EGS4SrcEns.SVTMP2 = EGS4SrcEns.TEMP2;
		EGS4SrcEns.SVTMP3 = EGS4SrcEns.TEMP3;
		EGS4SrcEns.SVTMP4 = EGS4SrcEns.TEMP4;
		// srcini
		EGS4SrcEns.iqin = EGS4SrcEns.ipart;// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
		EGS4SrcEns.NHSTRY = 0;
		EGS4SrcEns.last_nhstry = 0;

		// srcoto
		if (SOURCE != SOURCE_BEAM)
			EGS4SrcEns.IFPB = 1;
		else
			EGS4SrcEns.IFPB = 0;
		EGS4SrcEns.pi = Math.PI;
	}

	/**
	 * Fix some source and energy related variables. Called by init. Position, direction 
	 * cosines, particle weight and incident fluence are all fixed here.
	 */
	private void fixEmAll() {
		//double D = 0.0;
		double R2 = 0.0;
		int IXIN = 0;
		// NHSTRY
		EGS4SrcEns.NHSTRY = EGS4SrcEns.NHSTRY + 1;// @@@@@@@@@@@@@@@@@@@@@@@@@

		if (SOURCE == SOURCE_POINT) {
			EGS4SrcEns.DISTZ = EGS4SrcEns.TEMP1;// @@@@@@@@@@@@@@@@@@@@@@@@@@@
			EGS4SrcEns.RBEAM = EGS4SrcEns.TEMP2;// @@@@@@@@@@@@@@@@@@@@@@@@@@@
			EGS4SrcEns.RBEAM = EGS4SrcEns.RBEAM * EGS4SrcEns.$ONE_EPS;// @@@@@@@

			EGS4SrcEns.RBEAM2 = EGS4SrcEns.RBEAM * EGS4SrcEns.RBEAM;
			// ########FIX ainflu
			EGS4SrcEns.AINFLU = EGS4SrcEns.dble(EGS4SrcEns.NCASET)
					/ (4 * Math.PI * EGS4SrcEns.DISTZ * EGS4SrcEns.DISTZ);
			/*
			 * while(true) { EGS4SrcEns.xin=EGS4.random01();
			 * EGS4SrcEns.xin=(2.0*EGS4SrcEns.xin-1.0)*EGS4SrcEns.RBEAM;//@@@@
			 * 
			 * EGS4SrcEns.yin=EGS4.random01();
			 * EGS4SrcEns.yin=(2.0*EGS4SrcEns.yin-1.0)*EGS4SrcEns.RBEAM;//@@@@
			 * 
			 * R2=EGS4SrcEns.xin*EGS4SrcEns.xin+EGS4SrcEns.yin*EGS4SrcEns.yin;
			 * 
			 * if(R2<=EGS4SrcEns.RBEAM2) break; }// UNTIL R2.LE.RBEAM2;
			 */
			double diam = EGS4Geom.RCYL[EGS4Geom.NR] * EGS4SrcEns.$ONE_EPS;
			double costmax = EGS4SrcEns.DISTZ
					/ Math.sqrt(EGS4SrcEns.DISTZ * EGS4SrcEns.DISTZ + diam
							* diam);
			// costet evaluation-->polar angle
			double dom = (1 - costmax) / 2;// >0 and <1/2, costmax <1 and >0,
											// thetamax<90!
			double r = EGS4.random01();
			r = r * dom;
			double costet = 1 - 2 * r;// 2*r-1;//<0 always--- negativ z axis!!
			double sintet = Math.sqrt(1 - costet * costet);
			double tgtet = sintet / costet;
			r = EGS4.random01();

			double phi2 = 2 * Math.PI * r;// azimutal angle
			double teta = Math.abs(Math.atan(tgtet));// -pi/2,pi/2
			// initial we have formal directional cosines u,v,w = 0,0,-1 and
			// then we have teta and phi2:
			double u = Math.sin(teta) * Math.cos(phi2);// @@@@@@@@@@@@@@@@@@@@@@@@@@@
			double v = Math.sin(teta) * Math.sin(phi2);// @@@@@@@@@@@@@@@@@@@@@@@@@@@
			double w = costet;// <0//@@@@@@@@@@@@@@@@@@@@@@@@@@@
			source_parcurs = EGS4SrcEns.DISTZ / costet;
			source_parcurs = Math.abs(source_parcurs);// >0
	//		double x0 = 0.0;
	//		double y0 = 0.0;
	//		double x = x0 + u * source_parcurs;// @@@@@@@@@@@@@@@@@@@@@@@@@@@
	//		double y = y0 + v * source_parcurs;// @@@@@@@@@@@@@@@@@@@@@@@@@@@
	//		double z = EGS4Geom.ZPLANE[0];// 0.0;//@@@@@@@@@@@@@@@@@@@@@@@@@@@
			EGS4SrcEns.xin = 0.0;// x;
			EGS4SrcEns.yin = 0.0;// y;

			for (int IX = 1; IX <= EGS4Geom.NR; IX++) {
				IXIN = IX;
				if (R2 <= EGS4Geom.CYRAD2[IX - 1])
					break;
			}
			IXIN = 1;// !!!!!!!!!!
			EGS4SrcEns.irin = 2 + (IXIN - 1) * EGS4Geom.NZ;// @@@@@@@@@@@@@@@@@@@@@@@@
			EGS4SrcEns.zin = EGS4Geom.ZPLANE[0];// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@

			// D=Math.sqrt(R2+EGS4SrcEns.DISTZ*EGS4SrcEns.DISTZ);//point source
			// on symmetry z axis!!

			EGS4SrcEns.uin = u;// EGS4SrcEns.xin/D;//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
			EGS4SrcEns.vin = v;// EGS4SrcEns.yin/D;//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
			EGS4SrcEns.win = w;// EGS4SrcEns.DISTZ/D;//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
			EGS4SrcEns.NRCFLG = 10;
			// ##################fix weight
			double us1 = 2
					* Math.PI
					* (1 - EGS4SrcEns.DISTZ
							/ Math.sqrt(EGS4SrcEns.DISTZ * EGS4SrcEns.DISTZ +
							// EGS4SrcEns.RBEAM*EGS4SrcEns.RBEAM));
									diam * diam));
			EGS4SrcEns.WEIGHT = us1 / (4 * Math.PI);// FIXED
													// WEIGHT!!!@@@@@@@@@@@@@@
		} else if (SOURCE == SOURCE_SARPAGAN) {
			findNearestValue(energy, EGS4SrcEns.ein, true);
			int index_low = nearestposition;
			int index_high = 0;
			if (index_low < energy.length - 1)
				index_high = index_low + 1;
			else
				index_high = index_low;
			double sattc_interp = linInt(energy[index_high],
					satt_coeff[index_high], energy[index_low],
					satt_coeff[index_low], EGS4SrcEns.ein);
			smiu = smius * sattc_interp;// smiu=smiu*sattc_interp;
			// ========================
			// getCylinderRandom2(deltazup,winDetType,winThickness,e_incident);
			getCylinderRandom();
			// indet=atSurfaceDet();
			// ==========================
			// decreasing weight due to attenuation in sorce
			// if (sourceatt)
			// EGS4SrcEns.WEIGHT=EGS4SrcEns.WEIGHT*Math.exp(-smiu*source_parcurs);

		} else if (SOURCE == SOURCE_MARINELLI) {
			findNearestValue(energy, EGS4SrcEns.ein, true);
			int index_low = nearestposition;
			int index_high = 0;
			if (index_low < energy.length - 1)
				index_high = index_low + 1;
			else
				index_high = index_low;
			double sattc_interp = linInt(energy[index_high],
					satt_coeff[index_high], energy[index_low],
					satt_coeff[index_low], EGS4SrcEns.ein);
			smiu = smius * sattc_interp;// smiu=smiu*sattc_interp;
			// =====================================================================
			double r = EGS4.random01();
			if (EGS4.random01() < 0.5)
				r = 1.0 - r;
			double z1 = hsource * r;// -hsourceup;//0=>-hup;1=h-hup
			boolean infcyl = false;// test if inf cylinder exists
			if (hsource > hsourceup + hdettot)// EGS4Geom.ZPLANE[EGS4Geom.NPLANE-1])//hdet)
				infcyl = true;
			// if (z1<0)
			if (z1 < (hsourceup - ethick)) {
				getMCylSup(z1);
			} else if (z1 > hsourceup) {
				if (infcyl) {
					if (z1 > hsourceup + hdettot)// EGS4Geom.ZPLANE[EGS4Geom.NPLANE-1])
					{
						getMCylInf(z1);
					} else// middle region
					{
						getMMidInf(z1);
					}
				} else// middle region
				{
					getMMidInf(z1);
				}

			} else {
				// we are in sup plastic!!! adjust it!!
				while (true) {
					r = EGS4.random01();
					if (EGS4.random01() < 0.5)
						r = 1.0 - r;
					if (r != 1.0)
						break;
				}
				z1 = (hsourceup - ethick) * r;
				getMCylSup(z1);
			}

			// indet=atSurfaceDet();
			// ==========================
			// decreasing weight due to attenuation in sorce
			// if (sourceatt)
			// EGS4SrcEns.WEIGHT=EGS4SrcEns.WEIGHT*Math.exp(-smiu*source_parcurs);
		}
		if (SOURCE == SOURCE_BEAM) {
			// public static double dist_source_det=0.0;
			// public static double beam_angle=0.0;//degree
			// public static double beam_radius=0.0;
			// Sarpagan,Marinelli and source Point=>always u,v, independent of
			// xin,yin due to
			// 4 PI random emision!!
			// Here, Weight=1 and xin,yin are closed related to u and v!!!
			// ########FIX ainflu
			EGS4SrcEns.AINFLU = EGS4SrcEns.dble(EGS4SrcEns.NCASET)
					/ (Math.PI * beam_radius * beam_radius);
			// ##################
			double teta = Math.abs(Math.PI * beam_angle / 180.0);// with z
																	// axis!!
			double w = Math.cos(teta);
			double u = 0.0;
			double v = 0.0;
			R2 = 0.0;// double R2=0.0;
			double rbeam_at_frontFace = beam_radius;// /Math.cos(teta);//>beam_radius!!
			while (true) {// [
				EGS4SrcEns.xin = EGS4.random01();
				EGS4SrcEns.xin = (2.0 * EGS4SrcEns.xin - 1.0)
						* rbeam_at_frontFace;
				// $RANDOMSET YIN;YIN=(2.0*YIN-1.0)*RBEAM;
				EGS4SrcEns.yin = EGS4.random01();
				EGS4SrcEns.yin = (2.0 * EGS4SrcEns.yin - 1.0)
						* rbeam_at_frontFace;
				R2 = EGS4SrcEns.xin * EGS4SrcEns.xin + EGS4SrcEns.yin
						* EGS4SrcEns.yin;

				if (R2 <= rbeam_at_frontFace * rbeam_at_frontFace) {
					double cosphi2 = EGS4SrcEns.xin / R2;// 0.0;
					double sinphi2 = EGS4SrcEns.yin / R2;
					u = Math.sin(teta) * cosphi2;// @@@@@@@@@@@@@@@@@@@@@@@@@@@
					v = Math.sin(teta) * sinphi2;// @@@@@@@@@@@@@@@@@@@@@@@@@@@
					break;
				}
			}

			for (int IX = 1; IX <= EGS4Geom.NR; IX++) {
				IXIN = IX;
				if (R2 <= EGS4Geom.CYRAD2[IX - 1])
					break;
			}
			// System.out.println("ixin "+IXIN);
			EGS4SrcEns.irin = 2 + (IXIN - 1) * EGS4Geom.NZ;// @@@@@@@@@@@@@@@@@@@@@@@@
			EGS4SrcEns.zin = EGS4Geom.ZPLANE[0];// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@

			EGS4SrcEns.uin = u;
			EGS4SrcEns.vin = v;
			EGS4SrcEns.win = w;
			if (w != 1.0)
				EGS4SrcEns.IFPB = 0;// !!!!!!!!!
			EGS4SrcEns.NRCFLG = 10;
			// ##################fix weight
			EGS4SrcEns.WEIGHT = 1.0;
		}

	}

	/**
	 * Called by fixEmAll. It fixes Sarpagan geometry. 
	 */
	private static void getCylinderRandom() {
		cylSupB = false;
		double r = 0.0;
		while (true) {
			r = EGS4.random01();
			if (EGS4.random01() < 0.5)
				r = 1.0 - r;
			if (r != 1.0)
				break;
		}
		double z1 = (hsource - ethick) * r;// in source:0->hs-ethick//double
											// z1=hsource*r;
		EGS4SrcEns.DISTZ = hsource - z1;// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
		// EGS4SrcEns.AINFLU=EGS4SrcEns.dble(EGS4SrcEns.NCASET)/(4*Math.PI*EGS4SrcEns.DISTZ*EGS4SrcEns.DISTZ);//@@@@@@
		EGS4SrcEns.AINFLU = EGS4SrcEns.AINFLU + 1.0
				/ (4.0 * Math.PI * EGS4SrcEns.DISTZ * EGS4SrcEns.DISTZ);
		r = EGS4.random01();
		double ro1 = asource * Math.sqrt(r);// dist to z axis ro1:
		r = EGS4.random01();
		double phi1 = 2 * Math.PI * r;// random angle for x0,y0 eval. from phi1
										// and ro1!!
		// double d=hsource/2+z1+hwin;//distance point 0 to detector----is >0
		// always
		// double dso=hsource/2+z1;//distance point to out of source
		// for MC variance reduction we impose that all photons emerged from
		// point 0 will go into
		// an solid angle us1 having a distance d and a radius=2*detradius or
		// 2*sourceradius in order
		// to make sure that this allways contain the detector!!
		double diam = 0.0;
		if (asource <= EGS4Geom.RCYL[EGS4Geom.NR]) {
			diam = 2 * EGS4Geom.RCYL[EGS4Geom.NR];
		} else {
			diam = 2 * asource;
		}
		double us1 = 2
				* Math.PI
				* (1 - EGS4SrcEns.DISTZ
						/ Math.sqrt(EGS4SrcEns.DISTZ * EGS4SrcEns.DISTZ + diam
								* diam));
		// cosines of theta max.
		double costmax = EGS4SrcEns.DISTZ
				/ Math.sqrt(EGS4SrcEns.DISTZ * EGS4SrcEns.DISTZ + diam * diam);
		// costet evaluation-->polar angle
		double dom = (1 - costmax) / 2;// >0 and <1/2, costmax <1 and >0,
										// thetamax<90!
		r = EGS4.random01();
		r = r * dom;
		double costet = 1 - 2 * r;// 2*r-1;//<0 always--- negativ z axis!!
		double sintet = Math.sqrt(1 - costet * costet);
		double tgtet = sintet / costet;
		r = EGS4.random01();
		double phi2 = 2 * Math.PI * r;// azimutal angle
		EGS4SrcEns.WEIGHT = us1 / (4 * Math.PI);// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
		double teta = Math.abs(Math.atan(tgtet));// -pi/2,pi/2
		// initial we have formal directional cosines u,v,w = 0,0,-1 and then we
		// have teta and phi2:
		double u = Math.sin(teta) * Math.cos(phi2);// @@@@@@@@@@@@@@@@@@@@@@@@@@@
		double v = Math.sin(teta) * Math.sin(phi2);// @@@@@@@@@@@@@@@@@@@@@@@@@@@
		double w = costet;// <0//@@@@@@@@@@@@@@@@@@@@@@@@@@@

		source_parcurs = EGS4SrcEns.DISTZ / costet;
		source_parcurs = Math.abs(source_parcurs);// >0
		double x0 = ro1 * Math.cos(phi1);
		double y0 = ro1 * Math.sin(phi1);
		double x = x0;// +u*source_parcurs;//@@@@@@@@@@@@@@@@@@@@@@@@@@@
		double y = y0;// +v*source_parcurs;//@@@@@@@@@@@@@@@@@@@@@@@@@@@
		double z = z1;// EGS4Geom.ZPLANE[0];//0.0;//@@@@@@@@@@@@@@@@@@@@@@@@@@@
		// -----------------------------------
//		double ux = u;
//		double uy = v;
//		double uz = w;
		// strictly in source volume
	//	double l1 = 0.0;
	//	double zer = (ux * x0 + uy * y0) * (ux * x0 + uy * y0)
	//			+ (ux * ux + uy * uy) * (asource * asource - x0 * x0 - y0 * y0);
		/*
		 * if (zer>0.) { double
		 * s1=(-(ux*x0+uy*y0)-Math.sqrt((ux*x0+uy*y0)*(ux*x0+uy*y0)+
		 * (ux*ux+uy*uy)*(asource*asource-x0*x0-y0*y0)))/(ux*ux+uy*uy); double
		 * s2=(-(ux*x0+uy*y0)+Math.sqrt((ux*x0+uy*y0)*(ux*x0+uy*y0)+
		 * (ux*ux+uy*uy)*(asource*asource-x0*x0-y0*y0)))/(ux*ux+uy*uy); double
		 * s=0.; if ((s1<0.) && (s2>0)){ s=s2;} else
		 * {s=Math.min(Math.abs(s1),Math.abs(s2));}
		 * l1=s;//System.out.println(l1); } else
		 * l1=source_parcurs;//2;//something wrong is happened!! if
		 * (l1<source_parcurs)//2) //it also fly in air (neglected)
		 * source_parcurs=l1;//source_parcurs2=l1;
		 * //ro1>EGS4Geom.RCYL[EGS4Geom.NR] CASE!!!!!!!!!!!!
		 * if(ro1>EGS4Geom.RCYL[EGS4Geom.NR]) { //strictly in detector volume?
		 * double hitd=1.0e30;//infinity double l2=1.0e30;//infinity double
		 * zer2=(ux*x0+uy*y0)*(ux*x0+uy*y0)+
		 * (ux*ux+uy*uy)*(EGS4Geom.RCYL[EGS4Geom
		 * .NR]*EGS4Geom.RCYL[EGS4Geom.NR]-x0*x0-y0*y0); if (zer2>0.)//can hit
		 * the detector external surface { double
		 * s1=(-(ux*x0+uy*y0)-Math.sqrt((ux*x0+uy*y0)*(ux*x0+uy*y0)+
		 * (ux*ux+uy*uy
		 * )*(EGS4Geom.RCYL[EGS4Geom.NR]*EGS4Geom.RCYL[EGS4Geom.NR]-x0
		 * *x0-y0*y0)))/(ux*ux+uy*uy); double
		 * s2=(-(ux*x0+uy*y0)+Math.sqrt((ux*x0+uy*y0)*(ux*x0+uy*y0)+
		 * (ux*ux+uy*uy
		 * )*(EGS4Geom.RCYL[EGS4Geom.NR]*EGS4Geom.RCYL[EGS4Geom.NR]-x0
		 * *x0-y0*y0)))/(ux*ux+uy*uy); double s=0.; if ((s1<0.) && (s2>0)){
		 * s=s2;} else {s=Math.min(Math.abs(s1),Math.abs(s2));}
		 * l2=s;//System.out.println(l1); //test if hit the top plane of
		 * detector placed at DISTZ distance from point source!! double
		 * l22=Math.abs(EGS4SrcEns.DISTZ/uz);//(dist z initial/uz)just in case
		 * //hitd=Math.min(l2,l22); z=-EGS4SrcEns.DISTZ+w*l2;
		 * x=x0+u*l2;x=x*0.9999; y=y0+v*l2;y=y*0.9999;
		 * if(z>EGS4Geom.ZPLANE[EGS4Geom.NPLANE-1])//not hit {
		 * cylSupB=true;return; } if(z<EGS4Geom.ZPLANE[0])// || ) {
		 * z=EGS4Geom.ZPLANE[0];//test if hit superior plane x=x0+u*l22;
		 * y=y0+v*l22;
		 * 
		 * EGS4SrcEns.xin=x;EGS4SrcEns.yin=y;EGS4SrcEns.zin=z;
		 * EGS4SrcEns.uin=u;EGS4SrcEns.vin=v;EGS4SrcEns.win=w;
		 * 
		 * EGS4SrcEns.NRCFLG=10;//ONLY FRONT FACE!! int IXIN=0; double
		 * R2=EGS4SrcEns.xin*EGS4SrcEns.xin+EGS4SrcEns.yin*EGS4SrcEns.yin;
		 * for(int IX=1;IX<=EGS4Geom.NR;IX++) { IXIN=IX;
		 * if(R2<=EGS4Geom.CYRAD2[IX-1])break; }
		 * 
		 * EGS4SrcEns.irin=2+(IXIN-1)*EGS4Geom.NZ;//@@@@@@@@@@@@@@@@@@@@@@@@
		 * 
		 * return; } EGS4SrcEns.xin=x;EGS4SrcEns.yin=y;EGS4SrcEns.zin=z;
		 * EGS4SrcEns.uin=u;EGS4SrcEns.vin=v;EGS4SrcEns.win=w;
		 * 
		 * EGS4SrcEns.NRCFLG=20;//side!! int IZ1=0; for(int
		 * IZ=2;IZ<=EGS4Geom.NPLANE;IZ++) { IZ1=IZ;
		 * if(EGS4SrcEns.zin<=EGS4Geom.ZPLANE[IZ-1]) break; }
		 * EGS4SrcEns.irin=(EGS4Geom.NR-1)*EGS4Geom.NZ+IZ1;
		 * 
		 * 
		 * return;
		 * 
		 * } else { cylSupB=true;return; }
		 * 
		 * }//if(ro1>EGS4Geom.RCYL[EGS4Geom.NR])
		 */
		// END ro1>EGS4Geom.RCYL[EGS4Geom.NR] CASE!!!!!!!!!!!!
		EGS4SrcEns.xin = x;
		EGS4SrcEns.yin = y;
		EGS4SrcEns.zin = z;
		EGS4SrcEns.uin = u;
		EGS4SrcEns.vin = v;
		EGS4SrcEns.win = w;

		EGS4SrcEns.NRCFLG = 10;// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@TEMPORARY
								// ONLY FRONT FACE!!
		int IXIN = 0;
		double R2 = EGS4SrcEns.xin * EGS4SrcEns.xin + EGS4SrcEns.yin
				* EGS4SrcEns.yin;
		for (int IX = 1; IX <= EGS4Geom.NR; IX++) {
			IXIN = IX;
			if (R2 <= EGS4Geom.CYRAD2[IX - 1])
				break;
		}

		EGS4SrcEns.irin = 2 + (IXIN - 1) * EGS4Geom.NZ;// @@@@@@@@@@@@@@@@@@@@@@@@

	}

	@SuppressWarnings("unused")
	/**
	 * Not used here....algorithm changed!
	 * @return the result
	 */
	private static boolean atSurfaceDet() {
		boolean b = true;

		if (cylSupB)
			return false;
		double d2 = EGS4SrcEns.xin * EGS4SrcEns.xin + EGS4SrcEns.yin
				* EGS4SrcEns.yin;
		if (d2 > EGS4Geom.RCYL[EGS4Geom.NR] * EGS4Geom.RCYL[EGS4Geom.NR])
			b = false;

		if (EGS4SrcEns.zin < EGS4Geom.ZPLANE[0]
				|| EGS4SrcEns.zin > EGS4Geom.ZPLANE[EGS4Geom.NPLANE - 1])
			b = false;

		return b;
	}

	// flag=true->valoarea din sir sa fie mai mica--LOWER THAN VALUE
	// flag=false->valoarea din sir sa fie mai mare--NOT LOWER THAN VALUE
	/**
	 * Called by fixEmAll.
	 * Find the nearest value from an array relative to the given value.
	 * If flag is true then the array value should be lower than the given value. 
	 * If flag is set to false then the array value should be greater than the given value. 
	 * The input array is supposed to be already sorted.
	 * @param a the input array
	 * @param value the given value
	 * @param flag the flag
	 * @return the value from array
	 */
	public static double findNearestValue(double[] a, double value, boolean flag) {
		boolean b = true;
		int ip = 0;

		if (a.length > 1) {
			// double[] a1 = newQSort(a);
			while (b) {
				if (flag) {
					if ((a[ip] <= value) && (a[ip + 1] > value)) {
						break;
					}
				} else {
					if (ip > 0)
						if ((a[ip] >= value) && (a[ip - 1] < value)) {
							break;
						}
				}

				ip++;
				if (ip == a.length - 1) {
					b = false;
					break;
				}
			}
			nearestposition = ip;// ----------------
			return a[ip];
		} else {
			nearestposition = 0;// ---------------
			return a[0];
		}
	}

	/**
	 * Linear interpolation
	 * @param x1 first point x-value
	 * @param y1 first point y-value
	 * @param x2 second point x-value
	 * @param y2 second point y-value
	 * @param x desire point x-value
	 * @return desire point y-value
	 */
	public static double linInt(double x1, double y1, double x2, double y2,
			double x) {
		double result = -1.0;
		double[] mn = new double[2];
		// insucces
		mn[0] = -1.0;// m
		mn[1] = -1.0;// n
		double num = x1 - x2;
		if (num != 0.0) {
			mn[0] = (y1 - y2) / num;
			mn[1] = (x1 * y2 - y1 * x2) / num;
			result = mn[0] * x + mn[1];
		}
		return result;
	}

	/**
	 * Called by fixEmAll. It fixes Marinelli geometry. This refers to upper cylinder. 
	 * @param z1 z1
	 */
	private void getMCylSup(double z1) {
		cylSupB = false;
		double r = EGS4.random01();
		if (EGS4.random01() < 0.5)
			r = 1.0 - r;
		EGS4SrcEns.DISTZ = hsourceup - Math.abs(z1);// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
		// EGS4SrcEns.AINFLU=EGS4SrcEns.dble(EGS4SrcEns.NCASET)/(4*Math.PI*EGS4SrcEns.DISTZ*EGS4SrcEns.DISTZ);//@@@@@@
		EGS4SrcEns.AINFLU = EGS4SrcEns.AINFLU + 1.0
				/ (4.0 * Math.PI * EGS4SrcEns.DISTZ * EGS4SrcEns.DISTZ);
		r = EGS4.random01();
		double ro1 = asource * Math.sqrt(r);// dist to z axis ro1:
		r = EGS4.random01();
		double phi1 = 2 * Math.PI * r;// random angle for x0,y0 eval. from phi1
										// and ro1!!
		// double d=hsource/2+z1+hwin;//distance point 0 to detector----is >0
		// always
		// double dso=hsource/2+z1;//distance point to out of source
		// for MC variance reduction we impose that all photons emerged from
		// point 0 will go into
		// an solid angle us1 having a distance d and a radius=2*detradius or
		// 2*sourceradius in order
		// to make sure that this allways contain the detector!!
		double diam = 0.0;
		if (asource <= EGS4Geom.RCYL[EGS4Geom.NR]) {
			diam = 2 * EGS4Geom.RCYL[EGS4Geom.NR];
		} else {
			diam = 2 * asource;
		}
		double us1 = 2
				* Math.PI
				* (1 - EGS4SrcEns.DISTZ
						/ Math.sqrt(EGS4SrcEns.DISTZ * EGS4SrcEns.DISTZ + diam
								* diam));
		// cosines of theta max.
		double costmax = EGS4SrcEns.DISTZ
				/ Math.sqrt(EGS4SrcEns.DISTZ * EGS4SrcEns.DISTZ + diam * diam);
		// costet evaluation-->polar angle
		double dom = (1 - costmax) / 2;// >0 and <1/2, costmax <1 and >0,
										// thetamax<90!
		r = EGS4.random01();
		r = r * dom;
		double costet = 1 - 2 * r;// 2*r-1;//<0 always--- negativ z axis!!
		double sintet = Math.sqrt(1 - costet * costet);
		double tgtet = sintet / costet;
		r = EGS4.random01();
		double phi2 = 2 * Math.PI * r;// azimutal angle
		EGS4SrcEns.WEIGHT = us1 / (4 * Math.PI);// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
		double teta = Math.abs(Math.atan(tgtet));// -pi/2,pi/2
		// initial we have formal directional cosines u,v,w = 0,0,-1 and then we
		// have teta and phi2:
		double u = Math.sin(teta) * Math.cos(phi2);// @@@@@@@@@@@@@@@@@@@@@@@@@@@
		double v = Math.sin(teta) * Math.sin(phi2);// @@@@@@@@@@@@@@@@@@@@@@@@@@@
		double w = costet;// <0//@@@@@@@@@@@@@@@@@@@@@@@@@@@
		source_parcurs = EGS4SrcEns.DISTZ / costet;
		source_parcurs = Math.abs(source_parcurs);// >0
		double x0 = ro1 * Math.cos(phi1);
		double y0 = ro1 * Math.sin(phi1);
		double x = x0;// +u*source_parcurs;//@@@@@@@@@@@@@@@@@@@@@@@@@@@
		double y = y0;// +v*source_parcurs;//@@@@@@@@@@@@@@@@@@@@@@@@@@@
		double z = z1;// EGS4Geom.ZPLANE[0];//0.0;//@@@@@@@@@@@@@@@@@@@@@@@@@@@
		// -----------------------------------
	//	double ux = u;
	//	double uy = v;
	//	double uz = w;
		// strictly in source volume
	//	double l1 = 0.0;
	//	double zer = (ux * x0 + uy * y0) * (ux * x0 + uy * y0)
	//			+ (ux * ux + uy * uy) * (asource * asource - x0 * x0 - y0 * y0);
		/*
		 * if (zer>0.) { double
		 * s1=(-(ux*x0+uy*y0)-Math.sqrt((ux*x0+uy*y0)*(ux*x0+uy*y0)+
		 * (ux*ux+uy*uy)*(asource*asource-x0*x0-y0*y0)))/(ux*ux+uy*uy); double
		 * s2=(-(ux*x0+uy*y0)+Math.sqrt((ux*x0+uy*y0)*(ux*x0+uy*y0)+
		 * (ux*ux+uy*uy)*(asource*asource-x0*x0-y0*y0)))/(ux*ux+uy*uy); double
		 * s=0.; if ((s1<0.) && (s2>0)){ s=s2;} else
		 * {s=Math.min(Math.abs(s1),Math.abs(s2));}
		 * l1=s;//System.out.println(l1); } else
		 * l1=source_parcurs;//2;//something wrong is happened!! if
		 * (l1<source_parcurs)//2) //it also fly in air (neglected)
		 * source_parcurs=l1;//source_parcurs2=l1;
		 * //ro1>EGS4Geom.RCYL[EGS4Geom.NR] CASE!!!!!!!!!!!!
		 * if(ro1>EGS4Geom.RCYL[EGS4Geom.NR]) { //strictly in detector volume?
		 * double hitd=1.0e30;//infinity double l2=1.0e30;//infinity double
		 * zer2=(ux*x0+uy*y0)*(ux*x0+uy*y0)+
		 * (ux*ux+uy*uy)*(EGS4Geom.RCYL[EGS4Geom
		 * .NR]*EGS4Geom.RCYL[EGS4Geom.NR]-x0*x0-y0*y0); if (zer2>0.)//can hit
		 * the detector external surface { double
		 * s1=(-(ux*x0+uy*y0)-Math.sqrt((ux*x0+uy*y0)*(ux*x0+uy*y0)+
		 * (ux*ux+uy*uy
		 * )*(EGS4Geom.RCYL[EGS4Geom.NR]*EGS4Geom.RCYL[EGS4Geom.NR]-x0
		 * *x0-y0*y0)))/(ux*ux+uy*uy); double
		 * s2=(-(ux*x0+uy*y0)+Math.sqrt((ux*x0+uy*y0)*(ux*x0+uy*y0)+
		 * (ux*ux+uy*uy
		 * )*(EGS4Geom.RCYL[EGS4Geom.NR]*EGS4Geom.RCYL[EGS4Geom.NR]-x0
		 * *x0-y0*y0)))/(ux*ux+uy*uy); double s=0.; if ((s1<0.) && (s2>0)){
		 * s=s2;} else {s=Math.min(Math.abs(s1),Math.abs(s2));}
		 * l2=s;//System.out.println(l1); //test if hit the top plane of
		 * detector placed at DISTZ distance from point source!! double
		 * l22=Math.abs(EGS4SrcEns.DISTZ/uz);//(dist z initial/uz)just in case
		 * //hitd=Math.min(l2,l22); z=-EGS4SrcEns.DISTZ+w*l2;
		 * x=x0+u*l2;x=x*0.9999; y=y0+v*l2;y=y*0.9999;
		 * if(z>EGS4Geom.ZPLANE[EGS4Geom.NPLANE-1])//not hit {
		 * cylSupB=true;return; } if(z<EGS4Geom.ZPLANE[0])// || ) {
		 * z=EGS4Geom.ZPLANE[0];//test if hit superior plane x=x0+u*l22;
		 * y=y0+v*l22;
		 * 
		 * EGS4SrcEns.xin=x;EGS4SrcEns.yin=y;EGS4SrcEns.zin=z;
		 * EGS4SrcEns.uin=u;EGS4SrcEns.vin=v;EGS4SrcEns.win=w;
		 * 
		 * EGS4SrcEns.NRCFLG=10;//ONLY FRONT FACE!! int IXIN=0; double
		 * R2=EGS4SrcEns.xin*EGS4SrcEns.xin+EGS4SrcEns.yin*EGS4SrcEns.yin;
		 * for(int IX=1;IX<=EGS4Geom.NR;IX++) { IXIN=IX;
		 * if(R2<=EGS4Geom.CYRAD2[IX-1])break; }
		 * 
		 * EGS4SrcEns.irin=2+(IXIN-1)*EGS4Geom.NZ;//@@@@@@@@@@@@@@@@@@@@@@@@
		 * 
		 * return; } EGS4SrcEns.xin=x;EGS4SrcEns.yin=y;EGS4SrcEns.zin=z;
		 * EGS4SrcEns.uin=u;EGS4SrcEns.vin=v;EGS4SrcEns.win=w;
		 * 
		 * EGS4SrcEns.NRCFLG=20;//side!! int IZ1=0; for(int
		 * IZ=2;IZ<=EGS4Geom.NPLANE;IZ++) { IZ1=IZ;
		 * if(EGS4SrcEns.zin<=EGS4Geom.ZPLANE[IZ-1]) break; }
		 * EGS4SrcEns.irin=(EGS4Geom.NR-1)*EGS4Geom.NZ+IZ1;
		 * 
		 * 
		 * return;
		 * 
		 * } else { cylSupB=true;return; }
		 * 
		 * }//if(ro1>EGS4Geom.RCYL[EGS4Geom.NR])
		 */
		// END ro1>EGS4Geom.RCYL[EGS4Geom.NR] CASE!!!!!!!!!!!!
		EGS4SrcEns.xin = x;
		EGS4SrcEns.yin = y;
		EGS4SrcEns.zin = z;
		EGS4SrcEns.uin = u;
		EGS4SrcEns.vin = v;
		EGS4SrcEns.win = w;

		EGS4SrcEns.NRCFLG = 10;// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@TEMPORARY
								// ONLY FRONT FACE!!
		int IXIN = 0;
		double R2 = EGS4SrcEns.xin * EGS4SrcEns.xin + EGS4SrcEns.yin
				* EGS4SrcEns.yin;
		for (int IX = 1; IX <= EGS4Geom.NR; IX++) {
			IXIN = IX;
			if (R2 <= EGS4Geom.CYRAD2[IX - 1])
				break;
		}

		EGS4SrcEns.irin = 2 + (IXIN - 1) * EGS4Geom.NZ;// @@@@@@@@@@@@@@@@@@@@@@@@

	}

	/**
	 * Called by fixEmAll. It fixes Marinelli geometry. This refers to lower "cylinder", 
	 * i.e. the part of Marinelli baker below the detector (if it is the case). 
	 * @param z1 z1
	 */
	private void getMCylInf(double z1) {
		cylSupB = false;
		double r = EGS4.random01();
		double ro1 = Math.sqrt((asource * asource - (bsource + ethick)
				* (bsource + ethick))
				* r + (bsource + ethick) * (bsource + ethick));// dist to z axis
																// ro1:
		r = EGS4.random01();
		double phi1 = 2 * Math.PI * r;// random angle for x0,y0 eval. from phi1
										// and ro1!!
		// double d=Math.abs(z1+(hsource/2-hsourceup-hdet));//distance point 0
		// to detector>0
		double d = Math.abs(z1 - (hsourceup + hdettot));// -EGS4Geom.ZPLANE[EGS4Geom.NPLANE-1]);
		EGS4SrcEns.DISTZ = d;// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
		// EGS4SrcEns.AINFLU=EGS4SrcEns.dble(EGS4SrcEns.NCASET)/(4*Math.PI*EGS4SrcEns.DISTZ*EGS4SrcEns.DISTZ);//@@@@@@
		EGS4SrcEns.AINFLU = EGS4SrcEns.AINFLU + 1.0
				/ (4.0 * Math.PI * EGS4SrcEns.DISTZ * EGS4SrcEns.DISTZ);
		// for MC variance reduction we impose that all photons emerged from
		// point 0 will go into
		// an solid angle us1 having a distance d and a radius=2*sourceradius in
		// order
		// to make sure that this allways contain the detector!!
		double diam = 2 * asource;
		double us1 = 2 * Math.PI * (1 - d / Math.sqrt(d * d + diam * diam));
		double costmax = d / Math.sqrt(d * d + diam * diam);// cosines of theta
															// max.
		// costet evaluation-->polar angle
		double dom = (1 - costmax) / 2;// >0 and <1/2, costmax <1 and >0,
										// thetamax<90!
		r = EGS4.random01();
		r = r * dom;
		double costet = 2 * r - 1;// modif=><0 always--- neg z
									// axis!!!!!!!!!!!!!!
		double sintet = Math.sqrt(1 - costet * costet);
		double tgtet = sintet / costet;
		r = EGS4.random01();
		double phi2 = 2 * Math.PI * r;// azimutal angle
		EGS4SrcEns.WEIGHT = us1 / (4 * Math.PI);// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
		double teta = Math.abs(Math.atan(tgtet));// -pi/2,pi/2
		// initial we have formal directional cosines u,v,w = 0,0,-1 and then we
		// have teta and phi2:
		double u = Math.sin(teta) * Math.cos(phi2);
		double v = Math.sin(teta) * Math.sin(phi2);
		double w = costet;// <0 here
		source_parcurs = d / costet;
		source_parcurs = Math.abs(source_parcurs);// >0
		double x0 = ro1 * Math.cos(phi1);
		double y0 = ro1 * Math.sin(phi1);
		double x = x0;// +u*source_parcurs;
		double y = y0;// +v*source_parcurs;
		double z = z1;// EGS4Geom.ZPLANE[EGS4Geom.NPLANE-1];
		// -----------------------------------
	//	double ux = u;
	//	double uy = v;
	//	double uz = w;
		// strictly in source volume
	//	double l1 = 0.0;
	////	double zer = (ux * x0 + uy * y0) * (ux * x0 + uy * y0)
	//			+ (ux * ux + uy * uy) * (bsource * bsource - x0 * x0 - y0 * y0);
		/*
		 * if (zer>0.) { double
		 * s1=(-(ux*x0+uy*y0)-Math.sqrt((ux*x0+uy*y0)*(ux*x0+uy*y0)+
		 * (ux*ux+uy*uy)*(bsource*bsource-x0*x0-y0*y0)))/(ux*ux+uy*uy); double
		 * s2=(-(ux*x0+uy*y0)+Math.sqrt((ux*x0+uy*y0)*(ux*x0+uy*y0)+
		 * (ux*ux+uy*uy)*(bsource*bsource-x0*x0-y0*y0)))/(ux*ux+uy*uy); double
		 * s=0.; if ((s1<0.) && (s2>0)){ s=s2;} else
		 * {s=Math.min(Math.abs(s1),Math.abs(s2));}
		 * l1=s;//System.out.println(l1); } else l1=source_parcurs;//something
		 * wrong is happened!! if (l1<source_parcurs) //it also fly in air
		 * (neglected) source_parcurs=l1; //ro1>EGS4Geom.RCYL[EGS4Geom.NR]
		 * CASE!!!!!!!!!!!! if(ro1>EGS4Geom.RCYL[EGS4Geom.NR]) { //strictly in
		 * detector volume? double hitd=1.0e30;//infinity double
		 * l2=1.0e30;//infinity double zer2=(ux*x0+uy*y0)*(ux*x0+uy*y0)+
		 * (ux*ux+uy
		 * *uy)*(EGS4Geom.RCYL[EGS4Geom.NR]*EGS4Geom.RCYL[EGS4Geom.NR]-x0
		 * *x0-y0*y0); if (zer2>0.)//can hit the detector external surface {
		 * double s1=(-(ux*x0+uy*y0)-Math.sqrt((ux*x0+uy*y0)*(ux*x0+uy*y0)+
		 * (ux*ux
		 * +uy*uy)*(EGS4Geom.RCYL[EGS4Geom.NR]*EGS4Geom.RCYL[EGS4Geom.NR]-x0
		 * *x0-y0*y0)))/(ux*ux+uy*uy); double
		 * s2=(-(ux*x0+uy*y0)+Math.sqrt((ux*x0+uy*y0)*(ux*x0+uy*y0)+
		 * (ux*ux+uy*uy
		 * )*(EGS4Geom.RCYL[EGS4Geom.NR]*EGS4Geom.RCYL[EGS4Geom.NR]-x0
		 * *x0-y0*y0)))/(ux*ux+uy*uy); double s=0.; if ((s1<0.) && (s2>0)){
		 * s=s2;} else {s=Math.min(Math.abs(s1),Math.abs(s2));}
		 * l2=s;//System.out.println(l1); //test if hit the top plane of
		 * detector placed at DISTZ distance from point source!! double
		 * l22=Math.abs(EGS4SrcEns.DISTZ/uz);//(dist z initial/uz)just in case
		 * //hitd=Math.min(l2,l22); z=EGS4SrcEns.DISTZ+w*l2;
		 * x=x0+u*l2;x=x*0.9999; y=y0+v*l2;y=y*0.9999;
		 * if(z<EGS4Geom.ZPLANE[0])//>EGS4Geom.ZPLANE[EGS4Geom.NPLANE-1])//not
		 * hit { cylSupB=true;return; }
		 * if(z>EGS4Geom.ZPLANE[EGS4Geom.NPLANE-1])// || ) {
		 * z=EGS4Geom.ZPLANE[EGS4Geom.NPLANE-1]; x=x0+u*l22; y=y0+v*l22;
		 * 
		 * EGS4SrcEns.xin=x;EGS4SrcEns.yin=y;EGS4SrcEns.zin=z;
		 * EGS4SrcEns.uin=u;EGS4SrcEns.vin=v;EGS4SrcEns.win=w;
		 * //irin=1+IXIN*EGS4Geom.NZ;
		 * 
		 * // if(DISTRH==0.0) // { // D=Math.sqrt(R2+DISTB2); // } // else // {
		 * // D=Math.sqrt(R2+DISTRH*(DISTRH-2.0*yin)+DISTB2); // } //
		 * uin=xin/D;vin=(yin-DISTRH)/D;win=-DISTB/D; EGS4SrcEns.NRCFLG=30; int
		 * IXIN=0; double
		 * R2=EGS4SrcEns.xin*EGS4SrcEns.xin+EGS4SrcEns.yin*EGS4SrcEns.yin;
		 * for(int IX=1;IX<=EGS4Geom.NR;IX++) { IXIN=IX;
		 * if(R2<=EGS4Geom.CYRAD2[IX-1])break; }
		 * 
		 * EGS4SrcEns.irin=1+IXIN*EGS4Geom.NZ;//@@@@@@@@@@@@@@@@@@@@@@@@
		 * 
		 * return; } EGS4SrcEns.xin=x;EGS4SrcEns.yin=y;EGS4SrcEns.zin=z;
		 * EGS4SrcEns.uin=u;EGS4SrcEns.vin=v;EGS4SrcEns.win=w;
		 * 
		 * EGS4SrcEns.NRCFLG=20;//side!! int IZ1=0; for(int
		 * IZ=2;IZ<=EGS4Geom.NPLANE;IZ++) { IZ1=IZ;
		 * if(EGS4SrcEns.zin<=EGS4Geom.ZPLANE[IZ-1]) break; }
		 * EGS4SrcEns.irin=(EGS4Geom.NR-1)*EGS4Geom.NZ+IZ1;
		 * 
		 * 
		 * return;
		 * 
		 * } else { cylSupB=true;return; }
		 * 
		 * }//if(ro1>EGS4Geom.RCYL[EGS4Geom.NR])
		 */
		// END ro1>EGS4Geom.RCYL[EGS4Geom.NR] CASE!!!!!!!!!!!!
		EGS4SrcEns.xin = x;
		EGS4SrcEns.yin = y;
		EGS4SrcEns.zin = z;
		EGS4SrcEns.uin = u;
		EGS4SrcEns.vin = v;
		EGS4SrcEns.win = w;
		// irin=1+IXIN*EGS4Geom.NZ;

		// if(DISTRH==0.0)
		// {
		// D=Math.sqrt(R2+DISTB2);
		// }
		// else
		// {
		// D=Math.sqrt(R2+DISTRH*(DISTRH-2.0*yin)+DISTB2);
		// }
		// uin=xin/D;vin=(yin-DISTRH)/D;win=-DISTB/D;
		EGS4SrcEns.NRCFLG = 30;
		int IXIN = 0;
		double R2 = EGS4SrcEns.xin * EGS4SrcEns.xin + EGS4SrcEns.yin
				* EGS4SrcEns.yin;
		for (int IX = 1; IX <= EGS4Geom.NR; IX++) {
			IXIN = IX;
			if (R2 <= EGS4Geom.CYRAD2[IX - 1])
				break;
		}

		EGS4SrcEns.irin = 1 + IXIN * EGS4Geom.NZ;// @@@@@@@@@@@@@@@@@@@@@@@@

	}

	/**
	 * Called by fixEmAll. It fixes Marinelli geometry. This refers to side "cylinder", 
	 * i.e. the part of Marinelli baker surounding the detector at sides. 
	 * @param z1 z1
	 */
	private void getMMidInf(double z1) {
		double r = EGS4.random01();
		double ro1 = Math.sqrt((asource * asource - (bsource + ethick)
				* (bsource + ethick))
				* r + (bsource + ethick) * (bsource + ethick));// dist to z axis
																// ro1:
		r = EGS4.random01();
		double phi1 = 2 * Math.PI * r;// ///////for x0 and y0 eval
										// ???????????????????/
		// double
		// rdet=Math.sqrt(EGS4Geom.RCYL[EGS4Geom.NR]*EGS4Geom.RCYL[EGS4Geom.NR]+
		// EGS4Geom.ZPLANE[EGS4Geom.NPLANE-1]*EGS4Geom.ZPLANE[EGS4Geom.NPLANE-1]/4);//for
		// us eval!!
		double rdet = Math.sqrt(adettot * adettot + hdettot * hdettot / 4.0);// for
																				// us
																				// eval!!
		double d = ro1;// distance for us eval!!
		EGS4SrcEns.DISTZ = d;// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
		// EGS4SrcEns.AINFLU=EGS4SrcEns.dble(EGS4SrcEns.NCASET)/(4*Math.PI*EGS4SrcEns.DISTZ*EGS4SrcEns.DISTZ);//@@@@@@
		EGS4SrcEns.AINFLU = EGS4SrcEns.AINFLU + 1.0
				/ (4.0 * Math.PI * EGS4SrcEns.DISTZ * EGS4SrcEns.DISTZ);
		// to make sure that this allways contain the detector!!
		double diam = 2 * rdet;// for us eval!!
		double us1 = 2 * Math.PI * (1 - d / Math.sqrt(d * d + diam * diam));// OK!
		double costmax = d / Math.sqrt(d * d + diam * diam);// cosines of theta
															// max.
		// costet evaluation-->polar angle
		double dom = (1 - costmax) / 2;// >0 and <1/2, costmax <1 and >0,
										// thetamax<90!
		r = EGS4.random01();
		r = r * dom;
		double costet = -2 * r + 1;// >0 for instance
		double sintet = Math.sqrt(1 - costet * costet);
		double tgtet = sintet / costet;
		r = EGS4.random01();
		// for simplicity and due to simmetry of problem we choose that ***
		double phi2 = 2 * Math.PI * r;// phi1;//2*Math.PI*r;//azimutal angle
		EGS4SrcEns.WEIGHT = us1 / (4 * Math.PI);// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
		double teta = Math.abs(Math.atan(tgtet));// -pi/2,pi/2
		// initial we have formal directional cosines u,v,w and then we have
		// teta and phi2:
		double xsi = Math.PI / 2;// shift with Pi/2 relative to z axis!!!
		double u = Math.sin(xsi) * Math.cos(phi1 + Math.PI);// -changing
															// direction
		double v = Math.sin(xsi) * Math.sin(phi1 + Math.PI);// -changing
															// direction
		double w = Math.cos(xsi);// 0!
		// "new " direction
		double ux = 0.0;
		double uy = 0.0;
		double uz = 0.0;
		if (Math.abs(w) > 0.99999)// miuz->normal incident
		{
			ux = Math.sin(teta) * Math.cos(phi2);// new miux
			uy = Math.sin(teta) * Math.sin(phi2);// new miuy
			if (w < 0)
				uz = -Math.cos(teta);// new miuz--never here---random
			else
				uz = Math.cos(teta);// new miuz
		} else {
			double temp = Math.sqrt(1.0 - w * w);
			ux = Math.sin(teta) * (u * w * Math.cos(phi2) - v * Math.sin(phi2))
					/ temp + u * Math.cos(teta);
			uy = Math.sin(teta) * (v * w * Math.cos(phi2) + u * Math.sin(phi2))
					/ temp + v * Math.cos(teta);
			uz = -Math.sin(teta) * Math.cos(phi2) * temp + w * Math.cos(teta);
		}
		// --------------------------------------
		double x0 = ro1 * Math.cos(phi1);
		double y0 = ro1 * Math.sin(phi1);
		// -----------------------------------
	//	double zer = (ux * x0 + uy * y0)
	//			* (ux * x0 + uy * y0)
	//			+ (ux * ux + uy * uy)
	//			* (EGS4Geom.RCYL[EGS4Geom.NR] * EGS4Geom.RCYL[EGS4Geom.NR] - x0
	//					* x0 - y0 * y0);
		/*
		 * if (zer>0.) { //source_parcurs=-(ux*x0+uy*y0);//System.out.println(
		 * source_parcurs);>0 double
		 * s1=(-(ux*x0+uy*y0)-Math.sqrt((ux*x0+uy*y0)*(ux*x0+uy*y0)+
		 * (ux*ux+uy*uy
		 * )*(EGS4Geom.RCYL[EGS4Geom.NR]*EGS4Geom.RCYL[EGS4Geom.NR]-x0
		 * *x0-y0*y0)))/(ux*ux+uy*uy); double
		 * s2=(-(ux*x0+uy*y0)+Math.sqrt((ux*x0+uy*y0)*(ux*x0+uy*y0)+
		 * (ux*ux+uy*uy
		 * )*(EGS4Geom.RCYL[EGS4Geom.NR]*EGS4Geom.RCYL[EGS4Geom.NR]-x0
		 * *x0-y0*y0)))/(ux*ux+uy*uy);double s=0.; if ((s1<0.) && (s2>0)){
		 * s=s2;} else {s=Math.min(Math.abs(s1),Math.abs(s2));}
		 * source_parcurs=s+0.1;//System.out.println(source_parcurs); } else
		 * source_parcurs=0.;//not hit
		 */
		// ----------------------------------
		double x = x0;// +ux*source_parcurs;
		double y = y0;// +uy*source_parcurs;
		double z = z1;// +uz*source_parcurs;
		// attenuation is always in source:
		EGS4SrcEns.xin = x;
		EGS4SrcEns.yin = y;
		EGS4SrcEns.zin = z;
		EGS4SrcEns.uin = ux;
		EGS4SrcEns.vin = uy;
		EGS4SrcEns.win = uz;
		EGS4SrcEns.NRCFLG = 20;// side!!
		int IZ1 = 0;
		for (int IZ = 2; IZ <= EGS4Geom.NPLANE; IZ++) {
			IZ1 = IZ;
			if (EGS4SrcEns.zin <= EGS4Geom.ZPLANE[IZ - 1])
				break;
		}
		EGS4SrcEns.irin = (EGS4Geom.NR - 1) * EGS4Geom.NZ + IZ1;

	}
}
