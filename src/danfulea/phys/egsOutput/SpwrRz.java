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
 * SPWR APPLICATION <br>
 * SIMULATES THE PASSAGE OF AN ELECTRON OR PHOTON BEAM IN A FINITE, RIGHT CYLINDRICAL GEOMETRY. 
 * STOPPING POWER RATIO are computed in the volume of interest. SPW=(de/dx)/RHO. 
 * Based on EGS4/EGSnrc routines developed by SLAC(Stanford Linear Accelerator Center)/NRC (National Research Council of Canada).
 * @author Dan Fulea, 13 APR. 2007
 */

public class SpwrRz implements EgsQuestion {
	// public static double NORM_XRAY=1.0E12;
	public static int imoi = 5;// in main!~!!!!!!
	public static String mois = "";// NaI=>in main!~!!!!!!
	public static int ndowntotal = 0;
	public static int ndownscatt = 0;
	public static boolean scattCaseB = false;
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
	// ====================================================
	// REPLACE {$GOODSTATS} WITH {20}
	public static int $GOODSTATS = 20;
	// REPLACE {$MAXDOS} WITH {3} "NUMBER OF DOSE RELATED COMPONENTS "
	// REPLACE {$NSPR} WITH {2} "STOPPING POWER RATIO INDEX       "
	// "                                (1) STANDARD SPR --EDPN IT=1/2 "
	// "                                (2) SPR SANS STOPPERS  EDPN IT=3/4 "
	public static int $MAXDOS = 3;
	public static int $NSPR = 2;
	public static int ISPRREG = 0;
	public static int[] ISPR;// =new int[EGS4.$MXREG];
	/*
	 * "SCEDPN(IZ,IX,IT) = scores energy deposited in voxel IZ,IX for: " IT=1 --
	 * local material excluding bremsstrahlung created " below PCUT, atomic
	 * relaxations, e- created below " ECUT but including e- slowing across ECUT
	 * " IT=2 -- second material, same details as IT=1 " IT=3 -- local material
	 * with no stoppers " IT=4 -- second material with no stoppers " After
	 * analysis, SCEDPN(IZ,IX,1) is used to store the " stopping power ratio in
	 * voxel IZ,IX without stoppers and " SCEDPN(IZ,IX,3) stores the stopping
	 * power ratio in voxel " IZ,IX including stoppers
	 */
	public static double[][][] SCEDPN;// ($MAXZPLANE,$MAXRADII,$MAXIT),
	public static double[][][] SCEDPN2;// ($MAXZPLANE,$MAXRADII,$MAXIT),
	// SCDOSE($MAXZPLANE,$MAXRADII,$MAXDOS),
	// SCDOSE2($MAXZPLANE,$MAXRADII,$MAXDOS),
	public static double[][][] SCDOSE_COV;// ($MAXZPLANE,$MAXRADII,$MAXDOS-1),
	public static double[][][] SCSPR_COV;// ($MAXZPLANE,$MAXRADII,$NSPR),
	public static double DELTA = 0.0;
	public static double DELTAT = 0.0;
	public static double DELTLG = 0.0;
	public static int[][] SCDOSE_LAST;// ($MAXZPLANE,$MAXRADII),
	public static int[][] SCNUM;// ($MAXZPLANE,$MAXRADII);
	public static double[][][] SCEDPN_TMP;// ($MAXZPLANE,$MAXRADII,$MAXIT),
	public static double[][] SCDOSE_TMP3;// ($MAXZPLANE,$MAXRADII);
	public static int ICROSS;
	public static int ITRACKE;
	public static int[][] BSTAT;// ($MAXZPLANE,$NSPR),
	public static int[] BDSTAT;// ($MAXZPLANE);
	public static int NUM_STARTR = 0;// LIST FLUENCE START REGION
	public static int NUM_STOPR = 0;// LIST FLUENCE STOP REGION

	public static double ECUTMX = 0.0;
	public static double PCUTMX = 0.0;
	public static int IFLDIF = 0;
	public static int MAXIB = 0;
	public static int INEXT = 0;
	public static int ILOW = 0;
	public static int IHIGH = 0;

	/*
	 * public static double ECUTMX=0.0; public static double PCUTMX = 0.0;
	 * public static int IFLDIF = 0; public static int $MAXRZ=0;// public static
	 * int MAXIB=0; public static int IPRNTP=0;//PRINT FLUENCE SPECTRA public
	 * static final int IPRNTP_ALL=0; public static final int
	 * IPRNTP_CONDENSED=1; public static final int IPRNTP_SPECIFIED=2; public
	 * static final int IPRNTP_NONE=3; public static final int
	 * IPRIMARY_TOTAL_FLUENCE=0;//score total fluence spectra only public static
	 * final int IPRIMARY_ELECTRON_PRIMARIES=1;//score electron primaries
	 * //separately as well include brem generated photons as primaries public
	 * static final int IPRIMARY_INCLUDE_BREM_SECONDARIES=2;//score electron
	 * primaries but include //those generated by brem as secondaries public
	 * static final int IPRIMARY_PHOTON_PRIMARIES=3;//primaries are all first
	 * generation //photons,including brem.secondaries include all scattered
	 * photons and annihilations photons. public static final int
	 * IPRIMARY_ELECTRON_SECONDARIES=4;//score electron secondaries. public
	 * static int IPRIMARY=-1;//not fluence public static int
	 * NUM_STARTR=0;//LIST FLUENCE START REGION public static int
	 * NUM_STOPR=0;//LIST FLUENCE STOP REGION public static int INEXT=0; public
	 * static int[] IPRTSP;//IPRTSP($MXREG) public static int ILOW=0; public
	 * static int IHIGH=0; //stores fluence for particles of type IQL in energy
	 * bin IB and voxel IZ,IX. IT=1 for primaries. IT=2 for secondaries. public
	 * static double[][][][][] SCFLU;//(IB,IZ,IX,IQL,IT) //stores fluence
	 * squared for particles of type IQL in bin IB and voxel IZ,IX. Eventually
	 * holds //" uncertainty on SCFLU(IB,IZ,IX,IQL,IT). public static
	 * double[][][][][] SCFLU2; public static int[][][][][] SCFLU_LAST; //stores
	 * energy of particles IQL crossing voxel IZ,IX. IT same as SCFLU. public
	 * static double[][][][] SCEAVE;//SCEAVE(IZ,IX,IQL,IT) //stores energy
	 * squared of particles IQL crossing voxel IZ,IX. Eventually holds
	 * uncertainty in //" SCEAVE(IZ,IX,IQL,IT) public static double[][][][]
	 * SCEAVE2; public static double[][][][] SCEAVE_COV; //stores total fluence
	 * of particle type IQL in voxel IZ,IX. IT same as above. public static
	 * double[][][][] SCTFLU; public static double[][][][] SCTFLU2; public
	 * static int[][][][] SCTFLU_LAST; //stores fluence of charged particles in
	 * energy bin IB in voxel IZ,IX. IT same as above. public static
	 * double[][][][] SCFLEP;//SCFLEP(IB,IZ,IX,IT)=0.; public static
	 * double[][][][] SCFLEP2; public static int[][][][] SCFLEP_LAST; //stores
	 * energy of charged particles crossing voxel IZ,IX. public static
	 * double[][][] SCELEP;//SCELEP(IZ,IX,IT)=0.; public static double[][][]
	 * SCELEP2; public static double[][][] SCELEP_COV; //stores total fluence of
	 * charged particles crossing voxel IZ,IX. IT same as above. public static
	 * double[][][] SCTLEP; public static double[][][] SCTLEP2; public static
	 * int[][][] SCTLEP_LAST; //"SCFLU_TMP(IB,IZ,IX,IQL,IT) = temporary variable
	 * for SCFLU(IB,IZ,IX,IQL,IT) public static double[][][][][] SCFLU_TMP;
	 * //"SCEAVE_TMP(IZ,IX,IQL,IT) = temporary variable for SCEAVE(IZ,IX,IQL,IT)
	 * public static double[][][][] SCEAVE_TMP; //"SCSTP_TMP
	 * //"SCELEP_TMP(IZ,IX,IT) public static double[][][] SCELEP_TMP;
	 * //"SCFLEP_TMP(IB,IZ,IX,IT) more temporary variables public static
	 * double[][][][] SCFLEP_TMP; //"SCTFLU_TMP(IZ,IX,IQL,IT) public static
	 * double[][][][] SCTFLU_TMP; //"SCTLEP_TMP(IZ,IX,IT) public static
	 * double[][][] SCTLEP_TMP;
	 */
	// ===============================================

	public static int $EBIN = 100;// 500;//"MAX NUMBER OF BINS IN PULSE HEIGHT DISTN     "
									// Note if you use large values of EBIN so
									// that "
									// 2*EBIN > MXDATA (below), then use 2*EBIN
									// "

	// public static double[][] EBINW=new double[$EBIN][3];
	public static int $NSWTCH = 8;// "# OF NRC SWITCHES FOR CONTROLLING SCATTERING "

	public static int $NBATCH = 10;// "OUTPUT BATCHES                             "
	public static int JCASE = 0;; // "no. of histories per batch"
	public static int $NCASEMIN = 100;// "min. no. of histories                        "
	public static int $MXDATA;// =1040;//
								// "MAXIMUM DATA POINTS FOR ANALYSIS (i.e.($MXREG-1))"
	public static int $MAXIT = 4;// 7;//
									// "MAX # OF PARAMETERS TO BE SCORED             "
	// "                                (1) TOTAL DOSE                                "
	// "                                (2) STOPPERS AND DISCARDS DOSE                "
	// "                                (3) TOTAL DOSE FROM FRONT WALL                "
	// "                                (4) TOTAL DOSE FROM SIDE WALL                 "
	// "                                (5) TOTAL DOSE FROM BACK WALL                 "
	// "                                (6) TOTAL DOSE FROM INSIDE WALL               "
	// "                                (7) TOTAL DOSE SCORED IN REGION NSRCRG        "
	// "                                    DUE TO PARTICLES CREATED IN THE           "
	// "                                    REGION (ONLY WHEN ISOURC = 3)             "
	// "MAX # OF FLUENCE COMPONENTS TO BE SCORED"(1) TOTAL
	// FLUENCE"(2) PRIMARY FLUENCE"
	public static int $MAXCMPTS = 12;// $MAXIT;// "FOR THE GRID OUTPUTS"
	public static int $MAX_SC_PLANES = 1;// "required to use phase space macros"
	public static int $MAXBRSPLIT = 200;// "MAX BREM SPLITTING NUMBER"

	public static int NCASE = 0;
	public static int NCASEO = 0;
	public static int NCASET = 0;
	public static double[][] AMASS;// ($MAXZREG,$MAXRADII),
	public static double TMCPUO = 0.0;
	public static double TIMMAX = 0.0;
	public static double STATLM = 0.0;// EIN,
	public static int IDAT = 0;
	public static int IRESTART = 0;
	public static int DATCOUNT = 0;
	// "AMASS(IZ,IX) MASS OF ZONE WITH COORDINATES (IZ,IX)
	// "TMCPUO CPU TIME USED IN PREVIOUS SESSIONS
	// "TIMMAX MAXIMUM ALLOWED CPU HOURS FOR A GIVEN CALCULATION
	// "STATLM TARGET STATISTICS IN CAVITY USED FOR AN EARLY EXIT
	// "EIN KINETIC ENERGY OF THE EXTERNAL BEAM
	// "ISUMCV(NREG) THE ARRAY OF ZONES COMPRISING THE CAVITY REGION
	// "IDAT = 0 STORE DATA ARRAYS FOR RE-USE
	// " = 1 DON'T STORE THEM
	// " = 0 DOES NOT INCLUDE PHOTOELECTRON ANGLE SELECTION
	// "NCASE NUMBER OF HISTORIES REMAINING TO BE DONE
	// "NCASEO NUMBER OF HISTORIES DONE IN PREVIOUS SESSIONS
	// "NCASET NUMBER OF HISTORIES ALREADY DONE
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

	/*
	 * "SCDOSE(IZ,IX,IT) = scores energy deposited (for output as dose) in voxel
	 * " IZ,IX for: " IT=1 -- fraction of total dose due to stoppers " that
	 * exclude bremsstrahlung created " below PCUT, atomic relaxations and e- "
	 * created below ECUT but include e- slowing " across ECUT " IT=2 --
	 * fraction of total dose due to " bremsstrahlung created below PCUT, "
	 * atomic relaxations and e- created below " ECUT " IT=3 -- total dose
	 */
	public static double[][][] SCDOSE;
	public static double[][][] SCDOSE2;
	// ACCUMULATES ENERGY OF ELECTRONS CREATED (ENERGY OF ELECTRONS
	// " CREATED^2), EVENTUALLY HOLDS KERMA (UNCERTAINTY IN KERMA)
	public static double[][][] SCKERMA;
	public static double[][][] SCKERMA2;
	// ACCUMULATES (ENERGY DEPOSITED)*(ENERGY OF ELECTRONS CREATED)
	// " SUMMED OVER ALL PRIMARY HISTORIES, EVENTUALLY HOLDS UNCERTAINTY
	// " IN DOSE/KERMA
	public static double[][][] SCDOSEtoKERMA2;
	// ACCUMULATES COUNTS (COUNTS^2) IN EACH ENERGY BIN OF PULSE
	// " HEIGHT SENSITIVE REGION, EVENTUALLY HOLDS PULSE HEIGHT DISTN
	// " (UNCERTAINTY IN PULSE HEIGHT DISTN)
	public static double[] SCPDST = new double[$EBIN];
	public static double[] SCPDST2 = new double[$EBIN];
	// ACCUMULATES SUM OF COUNTS (COUNTS^2) IN ALL BINS <= CURRENT BIN
	// " IN PULSE HEIGHT REGION, EVENTUALLY HOLDS CUMULATIVE PULSE
	// " HEIGHT DISTN (UNCERTAINTY IN CUMULATIVE PULSE HEIGHT DISTN)
	public static double[] SCPCUM = new double[$EBIN];
	public static double[] SCPCUM2 = new double[$EBIN];
	public static double[] ECUM = new double[$EBIN];
	public static double[] ECUM2 = new double[$EBIN];
	public static double ETHRESHOLD = 0.0;
	public static double eefficiency = 0.0;
	public static double eefficiency_error = 0.0;
	// ACCUMULATES ALL COUNTS (COUNTS^2) IN ALL BINS IN PULSE HEIGHT
	// " REGION, EVENTUALLY HOLDS TOTAL OF ALL COUNTS (UNCERTAINTY ON
	// " TOTAL) AND IS USED TO NORMALIZE SCPDST, SCPCUM AND SCDFEP
	public static double SCPTOT = 0.0;
	public static double SCPTOT2 = 0.0;
	// ACCUMULATES TOTAL ENERGY (ENERGY^2) DEPOSITED IN PULSE HEIGHT
	// " REGION, EVENTUALL HOLDS ENERGY DEPOSITED IN PULSE HEIGHT REGION
	// " (UNCERTAINTY ON ENERGY IN PULSE HEIGHT REGION) PER INCIDENT
	// " PARTICLE
	public static double SCPHEN = 0.0;
	public static double SCPHEN2 = 0.0;
	// ACCUMULATES BACKGROUND COUNTS (COUNTS^2) FOR EACH PEAK IN PULSE
	// " HEIGHT DISTRIBUTION, WILL BE SUBTRACTED FROM SCDFEP BELOW
	public static double[] SCDFBK = new double[4];
	public static double[] SCDFBK2 = new double[4];
	// ACCUMULATES COUNTS (COUNTS^2) IN EACH ENERGY PEAK IN PULSE
	// " HEIGHT DISTRIBUTION
	public static double[] SCDFEP = new double[4];
	public static double[] SCDFEP2 = new double[4];
	// ACCUMULATES DIFFERENCE (DIFFERENCE^2) BETWEEN PEAK AND
	// " BACKGROUND COUNTS FOR EACH PEAK IN PULSE HEIGHT DISTRIBUTION,
	// " GETS USED TO DETERMINE COVARIANCE BETWEEN SCPTOT AND
	// " (SCDFEP-SCDFBK) SO WE CAN CALCULATE THE UNCERTAINTY ON THE
	// " FINAL VALUE OF SCDFEP
	public static double[] SCDFDIFF = new double[4];
	public static double[] SCDFDIFF2 = new double[4];
	// ACCUMULATES TOTAL PARTICLE STEPS (PARTICLE STEPS^2) IN PHANTOM,
	// " EVENTUALLY, SCSTP2 HOLDS UNCERTAINTY ON THIS NUMBER
	public static double SCSTP = 0.0;
	public static double SCSTP2 = 0.0;
	// ACCUMULATES TOTAL PARTICLE STEPS (PARTICLE STEPS^2) IN DOSE
	// " REGION OF PHANTOM, EVENTUALLY, SCDSTP2 HOLDS UNCERTAINTY ON
	// " THIS NUMBER
	public static double SCDSTP = 0.0;
	public static double SCDSTP2 = 0.0;
	public static double PIISTP = 0.0;// HOLDS NO. OF PRESTA-II STEPS FROM
										// PREVIOUS RUNS
	public static int SCSTP_LAST = 0;// LAST PRIMARY HISTORY TO SCORE A PARTICLE
										// STEP IN PHANTOM
	// LAST PRIMARY HISTORY TO SCORE A STEP IN THE DOSE REGION OF THE
	// " PHANTOM
	public static int SCDSTP_LAST = 0;
	// LAST PRIMARY HISTORY TO DEPOSIT ENERGY IN EACH DOSE REGION OF
	// " THE PHANTOM
	// public static int [][][] SCDOSE_LAST;//($MAXZREG,$MAXRADII,$MAXIT),
	// LAST PRIMARY HISTORY TO DEPOSIT KERMA IN EACH REGION OF
	// " PHANTOM
	public static int[][][] SCKERMA_LAST;// ($MAXZREG,$MAXRADII,$MAXIT),
	// SECOND LAST PRIMARY HISTORY TO DEPOSIT KERMA IN EACH REGION
	// " OF PHANTOM. USED TO CALCULATE SCDOSE*SCKERMA FOR COVARIANCE
	// " IN UNCERTAINTY ON DOSE/KERMA
	public static int[][][] SCKERMA_LASTOLD;// ($MAXZREG,$MAXRADII,$MAXIT),
	// LAST PRIMARY HISTORY TO SCORE ENERGY IN PULSE HEIGHT REGION
	public static int SCPDST_LAST = 0;
	// COUNTER FOR TOTAL NUMBER OF HISTORIES SUCCESSFULLY SIMULATED
	public static int IHSTRY = 0;
	// ENERGY LIMITS FOR FOUR PEAK AREA AND BACKGROUNDS
	public static double[][] DFEN = new double[4][4];// (4,4),
	// ENERGY DEPOSITED IN SENSITIVE VOLUME FOR CURRENT HISTORY ONLY
	public static double PHENER = 0.0;
	// stores value of WT(1) from last history for scoring phd
	public static double WT1OLD = 0.0;
	// TOPS OF ENERGY BINS FOR PULSE HEIGHT DISTRIBUTION
	public static double[] BINTOP = new double[$EBIN];// ($EBIN),
	// >0, WIDTH OF PULSE HEIGHT DISTRIBUTION ENERGY BINS
	// " ELSE FLAG TO USE BINTOP
	public static double SLOTE = 0.0;
	// WIDTH OF ENERGY BINS USED FOR PEAK EFFICIENCIES
	public static double DELTAE = 0.0;
	// ACCUMULATES ENERGY DEPOSITED IN EACH REGION BY CURRENT PRIMARY
	// " HISTORY
	public static double[][][] SCDOSE_TMP;// ($MAXZREG,$MAXRADII,$MAXIT),
	// ACCUMULATES ENERGY OF ELECTRONS CREATED IN EACH REGION BY
	// " CURRENT PRIMARY HISTORY
	public static double[][][] SCKERMA_TMP;// ($MAXZREG,$MAXRADII,$MAXIT),
	// ENERGY OF ELECTRONS CREATED IN EACH REGION BY LAST PRIMARY
	// " HISTORY. USED TO CALCULATE SCDOSE*SCKERMA FOR COVARIANCE
	// " IN UNCERTAINTY ON DOSE/KERMA
	public static double[][][] SCKERMA_TMPOLD;// ($MAXZREG,$MAXRADII,$MAXIT),
	// ACCUMULATES CHARGED PARTICLE STEPS IN CURRENT PRIMARY HISTORY
	public static double SCSTP_TMP = 0.0;
	// ACCUMULATES CHARGED PARTICLE STEPS IN DOSE REGION IN CURRENT
	// " PRIMARY HISTORY
	public static double SCDSTP_TMP = 0.0;
	// MAXIMUM LEVEL TO WHICH THE STACK OF DAUGHTER PARTICLES FROM AN
	// " INCIDENT PARTICLE RISES (STACK MAY INCLUDE INCIDENT PARTICLE)
	public static int MXNP = 0;
	// "IFULL = 0 JUST CALCULATE TOTAL DOSE AND THAT DUE TO STOPPERS
	// " AND DISCARDS (THE DEFAULT)
	// " = 1 ABOVE ANALYSE WHERE THE DOSE IS COMING FROM
	// " = 2 IFULL = 0 SCORING PLUS PULSE HEIGHT DISTRIBUTIONS
	// " = 3 SCORE THE SCATTER FRACTION INSTEAD OF STOPPERS

	// "ISTORE = 0 DO NOT STORE THE INITIAL RANDOM NUMBERS (THE DEFAULT)
	// " = 1 STORE THE INITIAL RANDOM NUMBER FOR LAST HISTORY
	// " = 2 STORE INITIAL RANDOM NUMBERS FOR ALL HISTORIES
	public static int ISTORE = 0;
	// "IKERMA = 0 DO NOT SCORE KERMA
	// " = 1 SCORE KERMA
	public static int IKERMA = 0;
	// "IWATCH = 0 FOR NORMAL OUTPUT (THE DEFAULT)
	// " = 1 OUTPUT ON EVERY DISCRETE INTERACTION
	// " = 2 OUTPUT ON EVERY ELECTRON/PHOTON STEP AS WELL
	// " = 3 PRINTS OUT ONLY WHEN ENERGY IS DEPOSITED
	// " = 4 PRINTS OUT FILE FOR GRAPHICS
	public static int IWATCH = 0;
	// "IOOPTN = 0 SHORT OUTPUT (THE DEFAULT) -JUST DOSE GRID(DG)
	// " = 1 OUTPUT DOSE SUMMARY ONLY (DS)
	// " = 2 OUTPUT MATERIAL SUMMARY GRID(MG) + DG
	// " = 3 OUTPUT MG + DS
	// " = 4 OUTPUT MG + DS + DG
	public static int IOOPTN = 0;
	// "IOUTSP = 0 NO ENERGY INPUT SPECTRUM DATA IN OUTPUT SUMMARY
	// " = 1 INCLUDE ENERGY INOUT SPECTRUM DATA IN OUTPUT SUMMARY
	public static int IOUTSP = 0;
	// FLAG ARRAY FOR EACH REGION, NON-ZERO ONLY IF PULSE HEIGHT
	// " DISTRIBUTION WANTED IN THIS GEOMETRIC REGION
	public static int[] IPHR;// ($MXREG),
	// NUMBER OF ENERGY BINS FOR PULSE HEIGHT DISTRIBUTION
	public static int MAXBIN = 0;
	public static int NCOMPT = 0;
	public static int IFANO = 0;// = 0 NO PHOTON REGENERATION
	// " = 1 PHOTONS REGENERATED AFTER THEY HAVE INTERACTED
	// " = 2 NO PHOTON REGENERATION, ELECTRONS FROM CAVITY WALL ARE ELIMINATED

	// "CDFINV INVERSE OF THE CUMULATIVE ENERGY PROBABILITY DISTRIBUTION
	// FUNCTION
	public static int $INVDIM = 1000;
	public static double[][] CDFINV = new double[$INVDIM][2];
	public static final String NULLs = "";
	public static final String BLANK = " ";
	public static final String DCHAR = "D";
	public static final String TCHAR = "T";
	public static final String ACHAR = "A";
	public static final String ASTER = "*";
	public static final String CCHAR = "C";
	public static final String BCHAR = "B";
	public static final String ICHAR = "I";
	public static final String MCHAR = "M";
	public static final String RCHAR = "R";
	public static final String SCHAR = "S";
	public static final String UCHAR = "U";
	public static final String VCHAR = "V";
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
	// public static final int ISTORE_ALL_DEPOSITED=2;
	public static final int ISTORE_ALL = 2;
	// @ IRESTART
	public static final int IRESTART_FIRST = 0;
	public static final int IRESTART_RESTART = 1;
	public static final int IRESTART_ANALYZE = 2;
	public static final int IRESTART_START_RNS = 3;
	public static final int IRESTART_PARALLEL = 4;
	// @ OUTPUT OPTIONS
	public static final int IOOPTN_SHORT = 0;
	public static final int IOOPTN_DOSE_SUMMARY = 1;
	public static final int IOOPTN_MATERIAL_SUMMARY = 2;
	public static final int IOOPTN_MATERIAL_AND_DOSE_SUMMARY = 3;
	public static final int IOOPTN_LONG = 4;
	// @ ELECTRON TRANSPORT
	public static final int ETRANS_NORMAL = 0;
	public static final int ETRANS_NO_INTERACTION = 1;
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
	// public static final int IFULL_AATT_AND_ASCAT=1;
	public static final int IFULL_ENTRANCE_REGIONS = 1;
	public static final int IFULL_PULSE_HEIGHT_DISTRIBUTION = 2;
	public static final int IFULL_SCATTER_FRACTION = 3;
	public static final int IFULL_OFMET_Fricke = 4;
	// @ STATISTICAL ACCURACY SOUGHT
	public static final double STATLM_MIN = 0.0;
	public static final double STATLM_MAX = 100.0;
	public static final double STATLM_DEFAULT = 0.0;
	// @ PHOTON REGENERATION
	public static final int IFANO_NO = 0;
	public static final int IFANO_YES = 1;
	public static final int IFANO_NO_ELECTRONS_FROM_WALL = 2;
	// @ SCORE KERMA
	public static final int KERMA_YES = 1;
	public static final int KERMA_NO = 0;
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
	public static double EI = 0.0;
	public static double EKMAX = 0.0;
	public static double DEPTH = 0.0;
	// public static double[][] VOLUME;//=0.0;
	public static double VOLUME;// =0.0;
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
	// ================================dose
	public static int ITMAX = 0;
	// "DECISION = 1 if ustep = 0 and region number changes--this ensures that
	// " LATCH bit gets set to reflect NEWNRC before dose is deposited
	public static int DECISION = 0;
	public static int NDOSE = 0;
	public static int NZDOSE = 0;
	public static int NRDOSE = 0;
	public static int MINZ = 0;
	public static int MAXZ = 0;
	public static int MINR = 0;
	public static int MAXR = 0;
	public static int IZD = 0;
	public static int IXD = 0;
	public static int NZDMIN = 0;
	public static int NZDMAX = 0;
	public static int NRDMIN = 0;
	public static int NRDMAX = 0;
	// "IDSTBL(IRL,1)=DOSE PLANAR SLAB NUMBER IF IRL IS IN DOSE REGION
	// "IDSTBL(IRL,2)=DOSE CYLINDRICAL REGION NUMBER IF IRL IS IN DOSE REGION
	public static int[][] IDSTBL;// ($MXREG,2);
	public static double EK0 = 0.0;
	public static double TDSMAX = 0.0;
	public static int IDSMAX = 0;
	public static double TDOS = 0.0;
	public static double TDOS2 = 0.0;
	public static int IB = 0;
	public static double SCORE_TEMP2 = 0.0;
	public static int IGBUG1 = 0;
	public static int IGBUG2 = 0;

	public static int nREGSOLV = 0;
	public static int[] REGSVOL = new int[1000];// =0;
	public static int nTOPEBIN = 0;
	public static double[] TOPEBIN = new double[1000];
	public static final double SLOTE_DEFAULT = 1.25;
	public static final double SLOTE_MIN = -10000.0;
	public static final double SLOTE_MAX = 10000.0;
	public static final double DELTAE_DEFAULT = 0.005;
	public static final double DELTAE_MIN = 1.e-20;
	public static final double DELTAE_MAX = 100000.0;
	public static final double TOPEBIN_DEFAULT = 1.25;
	public static final double TOPEBIN_MIN = 1.e-20;
	public static final double TOPEBIN_MAX = 100000.0;
	public static int nENHREG = 0;
	public static int[] NENHLO = new int[1000];
	public static int[] NENHHI = new int[1000];
	public static int ics_enhance = 0;
	public static double efficiency = 0.0;
	public static double efficiency_error = 0.0;
	public static double photontotal_efficiency = 0.0;
	public static double photontotal_efficiency_error = 0.0;

	// ========================
	public static int IPRINT = 2;
	public static boolean createOutputFile = false;
	public static boolean putInFile = false;// internal var defining when and
											// what to print
	private String filename = "";
	FileWriter sigfos;

	private boolean ekinokB = true;
	// -------------end INPUTS-------------------------------------------
	// REPLACE {$INVDIM} WITH {1000} "DIMENSION CONTROLS GRID SIZE FOR INVERSE"
	// CDFINV($INVDIM,2)

	public static boolean is_finished = false;

	/**
	 * Constructor
	 */
	public SpwrRz() {
		// createOutputFile=false;
		// createOutputFile=true;
		ndowntotal = 0;
		ndownscatt = 0;
		// FluenceRz.NCASE=100;
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
		// ETHRESHOLD=ETHRESHOLD+EGS4.RM;//kinetic
		ekinokB = true;// ETHRESHOLD affect in "hidden mode" only spectrum!!but
						// here no eff calculation is done!!
		// mono: directly ein <ethresh=>no shower!!
		// ==========================
		init();
	}

	/**
	 * Reset global variables for re-use
	 */
	public static void reset() {
		scattCaseB = false;
		nearestposition = 0;
		sourceatt = true;
		sourcecode = 2;// 1-H2O;2-NaCl;3=soil;4=ciment
		source_parcurs = 0.0;
		indet = true;
		asource = 0.0;
		hsource = 0.0;
		cylSupB = false;
		hsourceup = 0.0;
		bsource = 0.0;
		SOURCE = -1;// $EBIN=500;$NSWTCH=8;
		// $NBATCH=10;
		JCASE = 0;
		// $NCASEMIN=100;$MXDATA=1040;
		// $MAXIT=7;
		$MAXCMPTS = $MAXIT;// this case
		// $MAX_SC_PLANES=1;$MAXBRSPLIT=200;
		NCASE = 0;
		NCASEO = 0;
		NCASET = 0;
		TMCPUO = 0.0;
		TIMMAX = 0.0;
		STATLM = 0.0;
		IDAT = 0;
		IRESTART = 0;
		DATCOUNT = 0;
		RRZ = 0.0;
		RRCUT = 0.0;
		RUSROU = false;

		ETHRESHOLD = 0.0;
		eefficiency = 0.0;
		eefficiency_error = 0.0;
		SCPTOT = 0.0;
		SCPTOT2 = 0.0;
		SCPHEN = 0.0;
		SCPHEN2 = 0.0;

		SCSTP = 0.0;
		SCSTP2 = 0.0;
		SCDSTP = 0.0;
		SCDSTP2 = 0.0;
		PIISTP = 0.0;
		SCSTP_LAST = 0;
		SCDSTP_LAST = 0;
		SCPDST_LAST = 0;
		IHSTRY = 0;

		PHENER = 0.0;
		WT1OLD = 0.0;

		SLOTE = 0.0;
		DELTAE = 0.0;
		SCSTP_TMP = 0.0;
		SCDSTP_TMP = 0.0;
		MXNP = 0;
		ISTORE = 0;
		IKERMA = 0;
		IWATCH = 0;
		IOOPTN = 0;
		IOUTSP = 0;
		MAXBIN = 0;
		NCOMPT = 0;// $INVDIM=1000;

		SCPDST = new double[$EBIN];
		SCPDST2 = new double[$EBIN];
		SCPCUM = new double[$EBIN];
		SCPCUM2 = new double[$EBIN];
		ECUM = new double[$EBIN];
		ECUM2 = new double[$EBIN];
		SCDFBK = new double[4];
		SCDFBK2 = new double[4];
		SCDFEP = new double[4];
		SCDFEP2 = new double[4];
		SCDFDIFF = new double[4];
		SCDFDIFF2 = new double[4];
		DFEN = new double[4][4];// (4,4),
		BINTOP = new double[$EBIN];// ($EBIN),
		CDFINV = new double[$INVDIM][2];
		REGSVOL = new int[1000];// =0;
		TOPEBIN = new double[1000];
		NENHLO = new int[1000];
		NENHHI = new int[1000];

		// IPRTSP=new int[EGS4.$MXREG];
		// EBINW=new double[$EBIN][3];

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
		ispin = 0;// photon_xsection="";
		ESAVEIN = 0.0;
		IFARCE = 0;
		phsplitt = 1;
		phsplitt_MAX = 0;
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
		ITMAX = 0;
		DECISION = 0;
		NDOSE = 0;
		NZDOSE = 0;
		NRDOSE = 0;
		MINZ = 0;
		MAXZ = 0;
		MINR = 0;
		MAXR = 0;
		IZD = 0;
		IXD = 0;
		NZDMIN = 0;
		NZDMAX = 0;
		NRDMIN = 0;
		NRDMAX = 0;
		EK0 = 0.0;
		TDSMAX = 0.0;
		IDSMAX = 0;
		TDOS = 0.0;
		TDOS2 = 0.0;
		IB = 0;
		SCORE_TEMP2 = 0.0;
		IGBUG1 = 0;
		IGBUG2 = 0;
		nREGSOLV = 0;
		nTOPEBIN = 0;
		nENHREG = 0;
		ics_enhance = 0;
		efficiency = 0.0;
		efficiency_error = 0.0;
		photontotal_efficiency = 0.0;
		photontotal_efficiency_error = 0.0;
		IPRINT = 2;
		is_finished = false;

		ECUTMX = 0.0;
		PCUTMX = 0.0;
		IFLDIF = 0;// $MAXRZ=0;//
		MAXIB = 0;// IPRNTP=0;//PRINT FLUENCE SPECTRA ==DEFAULT
		// IPRIMARY=-1;//not fluence
		NUM_STARTR = 0;// LIST FLUENCE START REGION
		NUM_STOPR = 0;// LIST FLUENCE STOP REGION

		INEXT = 0;
		ILOW = 0;
		IHIGH = 0;

		ICROSS = 0;
		ITRACKE = 0;
		ISPRREG = 0;
		DELTA = 0.0;
		DELTAT = 0.0;
		DELTLG = 0.0;

	}

	/**
	 * Reset global array variables for re-use
	 */
	public static void destroyArrays() {

		ISPR = new int[0];// =new int[EGS4.$MXREG];
		SCEDPN = new double[0][0][0];// ($MAXZPLANE,$MAXRADII,$MAXIT),
		SCEDPN2 = new double[0][0][0];// ($MAXZPLANE,$MAXRADII,$MAXIT),
		SCDOSE = new double[0][0][0];// ($MAXZPLANE,$MAXRADII,$MAXDOS),
		SCDOSE2 = new double[0][0][0];// ($MAXZPLANE,$MAXRADII,$MAXDOS),
		SCDOSE_COV = new double[0][0][0];// ($MAXZPLANE,$MAXRADII,$MAXDOS-1),
		SCSPR_COV = new double[0][0][0];// ($MAXZPLANE,$MAXRADII,$NSPR),
		SCDOSE_LAST = new int[0][0];// ($MAXZPLANE,$MAXRADII),
		SCNUM = new int[0][0];// //($MAXZPLANE,$MAXRADII);
		SCEDPN_TMP = new double[0][0][0];// ($MAXZPLANE,$MAXRADII,$MAXIT),
		SCDOSE_TMP3 = new double[0][0];// ($MAXZPLANE,$MAXRADII);
		BSTAT = new int[0][0];// ($MAXZPLANE,$NSPR),
		BDSTAT = new int[0];// ($MAXZPLANE);

		/*
		 * SCFLU=new double[0][0][0][0][0];//DESTROY SCFLU2=new
		 * double[0][0][0][0][0];//DESTROY SCFLU_LAST=new
		 * int[0][0][0][0][0];//DESTROY SCFLU_TMP=new
		 * double[0][0][0][0][0];//DESTROY SCEAVE=new
		 * double[0][0][0][0];//DESTROY SCEAVE2=new double[0][0][0][0];//DESTROY
		 * SCEAVE_COV=new double[0][0][0][0];//DESTROY SCEAVE_TMP=new
		 * double[0][0][0][0];//DESTROY SCTFLU=new double[0][0][0][0];//DESTROY
		 * SCTFLU2=new double[0][0][0][0];//DESTROY SCTFLU_LAST=new
		 * int[0][0][0][0];//DESTROY SCTFLU_TMP=new double[0][0][0][0];//DESTROY
		 * SCFLEP=new double[0][0][0][0];//DESTROY SCFLEP2=new
		 * double[0][0][0][0];//DESTROY SCFLEP_LAST=new
		 * int[0][0][0][0];//DESTROY SCFLEP_TMP=new double[0][0][0][0];//DESTROY
		 * SCELEP=new double[0][0][0];//DESTROY SCELEP2=new
		 * double[0][0][0];//DESTROY SCELEP_COV=new double[0][0][0];//DESTROY
		 * SCELEP_TMP=new double[0][0][0];//DESTROY SCTLEP=new
		 * double[0][0][0];//DESTROY SCTLEP2=new double[0][0][0];//DESTROY
		 * SCTLEP_LAST=new int[0][0][0];//DESTROY SCTLEP_TMP=new
		 * double[0][0][0];//DESTROY
		 */

		// VOLUME=new double[0][0];//at init
		AMASS = new double[0][0];// at init
		SCPDST = new double[0];
		SCPDST2 = new double[0];
		SCPCUM = new double[0];
		SCPCUM2 = new double[0];
		ECUM = new double[0];
		ECUM2 = new double[0];
		SCDFBK = new double[0];
		SCDFBK2 = new double[0];
		SCDFEP = new double[0];
		SCDFEP2 = new double[0];
		SCDFDIFF = new double[0];
		SCDFDIFF2 = new double[0];
		DFEN = new double[0][0];// (4,4),
		BINTOP = new double[0];// ($EBIN),
		CDFINV = new double[0][0];
		REGSVOL = new int[0];// =0;
		TOPEBIN = new double[0];
		NENHLO = new int[0];
		NENHHI = new int[0];
		// IPRTSP=new int[0];
		// EBINW=new double[0][0];

		// ===================from init:
		// IPRTSP=new int[0];
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
		// System.out.println("init");

		EGS4.startSimulationTime = System.currentTimeMillis();
		// --variable init
		// ======================================================================================
		// IPRTSP=new int[EGS4.$MXREG];
		// ======================================================================================
		$MXDATA = $MAXIT * $EBIN * EGS4Geom.$MAXRADII * EGS4Geom.$MAXZREG;// COMPUTE
																			// $MXREG-1
		// REPLACE {$MXDATA} WITH {{COMPUTE
		// $MAXIT*$EBIN*$MAXRADII*($MAXZPLANE-1) }}
		// REPLACE {$MAXZPLANE} WITH {{COMPUTE $MAXZREG+1}}=at geom.reset(

		// EGS4SrcEns.$MXRDIST=1000;//$MXRDIST->4
		// EGS4SrcEns.$NENSRC=500;//300;
		EGS4Geom.$NVALUE = 100;
		// ======================
		nREGSOLV = GammaDetEffA.nREGSOLV;
		for (int i = 1; i <= nREGSOLV; i++)// at reset GammaDetEffA, then
											// construct in main
		{
			REGSVOL[i - 1] = GammaDetEffA.REGSVOL[i - 1];
		}
		// =======================

		AMASS = new double[EGS4Geom.$MAXZREG][EGS4Geom.$MAXRADII];// ($MAXZREG,$MAXRADII),

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
		// ====================================================================================
		ISPR = new int[EGS4.$MXREG];// =new int[EGS4.$MXREG];
		SCEDPN = new double[EGS4Geom.$MAXZPLANE][EGS4Geom.$MAXRADII][$MAXIT];// ($MAXZPLANE,$MAXRADII,$MAXIT),
		SCEDPN2 = new double[EGS4Geom.$MAXZPLANE][EGS4Geom.$MAXRADII][$MAXIT];// ($MAXZPLANE,$MAXRADII,$MAXIT),
		SCDOSE = new double[EGS4Geom.$MAXZPLANE][EGS4Geom.$MAXRADII][$MAXDOS];// ($MAXZPLANE,$MAXRADII,$MAXDOS),
		SCDOSE2 = new double[EGS4Geom.$MAXZPLANE][EGS4Geom.$MAXRADII][$MAXDOS];// ($MAXZPLANE,$MAXRADII,$MAXDOS),
		SCDOSE_COV = new double[EGS4Geom.$MAXZPLANE][EGS4Geom.$MAXRADII][$MAXDOS - 1];// ($MAXZPLANE,$MAXRADII,$MAXDOS-1),
		SCSPR_COV = new double[EGS4Geom.$MAXZPLANE][EGS4Geom.$MAXRADII][$NSPR];// ($MAXZPLANE,$MAXRADII,$NSPR),
		SCDOSE_LAST = new int[EGS4Geom.$MAXZPLANE][EGS4Geom.$MAXRADII];// ($MAXZPLANE,$MAXRADII),
		SCNUM = new int[EGS4Geom.$MAXZPLANE][EGS4Geom.$MAXRADII];// //($MAXZPLANE,$MAXRADII);
		SCEDPN_TMP = new double[EGS4Geom.$MAXZPLANE][EGS4Geom.$MAXRADII][$MAXIT];// ($MAXZPLANE,$MAXRADII,$MAXIT),
		SCDOSE_TMP3 = new double[EGS4Geom.$MAXZPLANE][EGS4Geom.$MAXRADII];// ($MAXZPLANE,$MAXRADII);
		BSTAT = new int[EGS4Geom.$MAXZPLANE][$NSPR];// ($MAXZPLANE,$NSPR),
		BDSTAT = new int[EGS4Geom.$MAXZPLANE];// ($MAXZPLANE);

		/*
		 * SCFLU=new
		 * double[$EBIN][EGS4Geom.$MAXZPLANE][EGS4Geom.$MAXRADII][3][$MAXIT
		 * ];//($EBIN,$MAXZPLANE,$MAXRADII,3,$MAXIT) SCFLU2=new
		 * double[$EBIN][EGS4Geom.$MAXZPLANE][EGS4Geom.$MAXRADII][3][$MAXIT];
		 * SCFLU_LAST=new
		 * int[$EBIN][EGS4Geom.$MAXZPLANE][EGS4Geom.$MAXRADII][3][$MAXIT];
		 * SCFLU_TMP=new
		 * double[$EBIN][EGS4Geom.$MAXZPLANE][EGS4Geom.$MAXRADII][3][$MAXIT];
		 * SCEAVE=new
		 * double[EGS4Geom.$MAXZPLANE][EGS4Geom.$MAXRADII][3][$MAXIT];
		 * //SCEAVE($MAXZPLANE,$MAXRADII,3,$MAXIT), SCEAVE2=new
		 * double[EGS4Geom.$MAXZPLANE][EGS4Geom.$MAXRADII][3][$MAXIT];
		 * SCEAVE_COV=new
		 * double[EGS4Geom.$MAXZPLANE][EGS4Geom.$MAXRADII][3][$MAXIT];
		 * SCEAVE_TMP=new
		 * double[EGS4Geom.$MAXZPLANE][EGS4Geom.$MAXRADII][3][$MAXIT];
		 * SCTFLU=new
		 * double[EGS4Geom.$MAXZPLANE][EGS4Geom.$MAXRADII][3][$MAXIT];
		 * //SCTFLU($MAXZPLANE,$MAXRADII,3,$MAXIT), SCTFLU2=new
		 * double[EGS4Geom.$MAXZPLANE][EGS4Geom.$MAXRADII][3][$MAXIT];
		 * SCTFLU_LAST=new
		 * int[EGS4Geom.$MAXZPLANE][EGS4Geom.$MAXRADII][3][$MAXIT];
		 * SCTFLU_TMP=new
		 * double[EGS4Geom.$MAXZPLANE][EGS4Geom.$MAXRADII][3][$MAXIT];
		 * SCFLEP=new
		 * double[$EBIN][EGS4Geom.$MAXZPLANE][EGS4Geom.$MAXRADII][$MAXIT
		 * ];//SCFLEP($EBIN,$MAXZPLANE,$MAXRADII,$MAXIT), SCFLEP2=new
		 * double[$EBIN][EGS4Geom.$MAXZPLANE][EGS4Geom.$MAXRADII][$MAXIT];
		 * SCFLEP_LAST=new
		 * int[$EBIN][EGS4Geom.$MAXZPLANE][EGS4Geom.$MAXRADII][$MAXIT];
		 * SCFLEP_TMP=new
		 * double[$EBIN][EGS4Geom.$MAXZPLANE][EGS4Geom.$MAXRADII][$MAXIT];
		 * SCELEP=new
		 * double[EGS4Geom.$MAXZPLANE][EGS4Geom.$MAXRADII][$MAXIT];//SCELEP
		 * ($MAXZPLANE,$MAXRADII,$MAXIT), SCELEP2=new
		 * double[EGS4Geom.$MAXZPLANE][EGS4Geom.$MAXRADII][$MAXIT];
		 * SCELEP_COV=new
		 * double[EGS4Geom.$MAXZPLANE][EGS4Geom.$MAXRADII][$MAXIT];
		 * SCELEP_TMP=new
		 * double[EGS4Geom.$MAXZPLANE][EGS4Geom.$MAXRADII][$MAXIT]; SCTLEP=new
		 * double
		 * [EGS4Geom.$MAXZPLANE][EGS4Geom.$MAXRADII][$MAXIT];//SCTLEP($MAXZPLANE
		 * ,$MAXRADII,$MAXIT), SCTLEP2=new
		 * double[EGS4Geom.$MAXZPLANE][EGS4Geom.$MAXRADII][$MAXIT];
		 * SCTLEP_LAST=new int[EGS4Geom.$MAXZPLANE][EGS4Geom.$MAXRADII][$MAXIT];
		 * SCTLEP_TMP=new
		 * double[EGS4Geom.$MAXZPLANE][EGS4Geom.$MAXRADII][$MAXIT];
		 */
		// VOLUME=new double[EGS4Geom.$MAXZREG][EGS4Geom.$MAXRADII];
		AMASS = new double[EGS4Geom.$MAXZREG][EGS4Geom.$MAXRADII];
		// ===================================================================================
		/*
		 * SCDOSE=new double[EGS4Geom.$MAXZREG][EGS4Geom.$MAXRADII][$MAXIT];
		 * SCDOSE2=new double[EGS4Geom.$MAXZREG][EGS4Geom.$MAXRADII][$MAXIT];
		 * SCKERMA=new double[EGS4Geom.$MAXZREG][EGS4Geom.$MAXRADII][$MAXIT];
		 * SCKERMA2=new double[EGS4Geom.$MAXZREG][EGS4Geom.$MAXRADII][$MAXIT];
		 * SCDOSEtoKERMA2=new
		 * double[EGS4Geom.$MAXZREG][EGS4Geom.$MAXRADII][$MAXIT];
		 * SCDOSE_LAST=new
		 * int[EGS4Geom.$MAXZREG][EGS4Geom.$MAXRADII][$MAXIT];//(
		 * $MAXZREG,$MAXRADII,$MAXIT), SCKERMA_LAST=new
		 * int[EGS4Geom.$MAXZREG][EGS4Geom.$MAXRADII][$MAXIT];
		 * SCKERMA_LASTOLD=new
		 * int[EGS4Geom.$MAXZREG][EGS4Geom.$MAXRADII][$MAXIT]; SCDOSE_TMP=new
		 * double[EGS4Geom.$MAXZREG][EGS4Geom.$MAXRADII][$MAXIT];
		 * SCKERMA_TMP=new
		 * double[EGS4Geom.$MAXZREG][EGS4Geom.$MAXRADII][$MAXIT];
		 * SCKERMA_TMPOLD=new
		 * double[EGS4Geom.$MAXZREG][EGS4Geom.$MAXRADII][$MAXIT];
		 * //($MAXZREG,$MAXRADII,$MAXIT) IPHR=new int[EGS4.$MXREG]; IDSTBL=new
		 * int[EGS4.$MXREG][2];//($MXREG,2);
		 */// @@@@@@@@@@@@@@@@@@@@
		EGS4Grid.CDSTBL = new String[EGS4.$MXREG];
		EGS4Grid.CTRTBL = new String[EGS4.$MXREG];
		EGS4Grid.CABSRB = new String[EGS4.$MXREG];
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
		// EGS4.ranluxB=false;//use ranmar!!
		EGS4.ispmfp = EGS4.iSpr;// EGS4.ispmfp=EGS4.iFlur;//EGS4.iDose;//select
								// photon mfp
		EGS4.iurd = 0;// no user range discard
		// EGS4.iraycorr=EGS4.iDose;//rayleigh correction NO USED
		EGS4.isemfp = EGS4.iDose;// $SELECT-ELECTRON-MFPCOINCIDE@@@@@@@@@@@@
		EGS4.iGeom = 1;// EGS4.iCavity;//1=RZ GEOM@@@@@@@@@@@@@@YES
		EGS4Macro.ismfpfpb = EGS4.iCavity;// select mfp parallel beam
											// COINCIDE@@@@@@@@@@@@
		// EGS4Macro.irange_rej=EGS4.iCavity;//not used range rejection

		EGS4Grid.$MAXCMPTS = $MAXCMPTS;
		// -----------------
		// EGS4.USER_CONTROLS_TSTEP_RECURSION=EGS4.iDose;//no effect here
		EGS4.USER_CONTROLS_TSTEP_RECURSION = EGS4.iSpr;
		EGS4.hatchindex = EGS4.iDose;// no effect here
		// -----------------
		EGS4.iprint = IPRINT;// EGS4.iprint=2;//summary
		// inputs
		// is_finished = true;
		// "******************************************************************************
		// "
		// " *** SECTION 1 ***
		// "
		// "------------------------------------------------------------------------------
		// "
		// "READ INPUTS AND CALCULATE ONE-TIME ONLY CONSTANTS
		// "
		// "------------------------------------------------------------------------------

		//Calendar cal = Calendar.getInstance();
		//Date d = cal.getTime();
		EGS4.seqStr = "=================================================================================";
		EGS4.seqStr = " **********************SPWR: APPLICATION****************************************";
		if (EGS4.iprint > 1)
			printSequence(EGS4.seqStr);
		EGS4.seqStr = " SIMULATES THE PASSAGE OF AN ELECTRON OR PHOTON BEAM IN A FINITE, RIGHT CYLINDRICAL GEOMETRY.";
		if (EGS4.iprint > 1)
			printSequence(EGS4.seqStr);
		EGS4.seqStr = "STOPPING POWER RATIO are computed in the volume of interest. SPW=(de/dx)/RHO";
		if (EGS4.iprint > 1)
			printSequence(EGS4.seqStr);
		EGS4.seqStr = " ********************************************************************************";
		if (EGS4.iprint > 1)
			printSequence(EGS4.seqStr);
		inputs();
		if (EGS4.STOPPROGRAM) {
			closeFile();
			destroyArrays();
			return;
		}

		if (SOURCE == SOURCE_POINT || SOURCE == SOURCE_SARPAGAN
				|| SOURCE == SOURCE_MARINELLI || SOURCE == SOURCE_BEAM) {
			fixSRCOTO();
		}
		// " This is copied directly from below.  Don't know why it needs "
		// " repeating, but it does.  -- JT "
		// if((EGS4Macro.IFULL==0)||(EGS4Macro.IFULL==2)||(EGS4Macro.IFULL==3)){ITMAX=2;}
		// else{ITMAX=$MAXIT;}
		// if(IPRIMARY == 0){ITMAX=1;}else{ITMAX=$MAXIT;}
		// //"# OF FLUENCE COMPONENTS"
		if (NCASE / $NBATCH == 0) {
			NCASE = $NBATCH;
		}
		JCASE = NCASE / $NBATCH;
		NCASE = JCASE * $NBATCH;// "NUMBER OF HISTORIES PER BATCH

		// DECISION=0; //"this is for howfar and ausgab when IFULL=1"

		MXNP = 0; // "reset the maximum stack indicator"
		IHSTRY = NCASEO; // "reset the number of histories counter"
		// =================================
		// EGS4Geom.NPLANE=EGS4Geom.NZ+1; //"NUMBER OF PLANAR BOUNDARIES"=>GEOM
		// if(EGS4SrcEns.iqin ==
		// 0){EI=EGS4SrcEns.ein;}else{EI=EGS4SrcEns.ein+EGS4.RM;}//
		// "TOTAL ENERGY"
		// LATER!!!!!!!!!!!!!!!!!!!!!
		// ===============================
		EGS4SrcEns.NHSTRY = 0; // "start the no. of primary histories counter at zero"

		// "set up ausgab calls"
		// DO J=1,5[IAUSFL(J)=1;]
		for (int J = 1; J <= 5; J++) {
			EGS4.iausfl[J - 1] = 1;
		}
		// DO J=6,25[IAUSFL(J)=0;]
		for (int J = 6; J <= 25; J++) {
			EGS4.iausfl[J - 1] = 0;
		}
		// "now set up to call AUSGAB under various circumstances so the flags"
		// "can be set                                                         "
		// " IAUSFL(8)=1; " "after brem created"
		// " IAUSFL(10)=1; " "after a Moller interaction occurs"
		// " IAUSFL(12)=1; " "after a Bhabha"
		// " IAUSFL(14)=1; " "after ANNIH"
		// " IAUSFL(15)=1; " "after annihilation at rest"
		EGS4.iausfl[16] = 1;// (17)=1; "after a pair event"
		EGS4.iausfl[18] = 1;// (19)=1; "after a Compton scatter"
		EGS4.iausfl[20] = 1;// (21)=1; "after a p.e. event"

		if (IFANO == 1) {
			// "AUSGAB will be responsible for making sure that the beam is not"
			// "attenuated and getting rid of the scattered photons."
			EGS4.iausfl[15] = 1;// (16) = 1; "Before pair"
			EGS4.iausfl[16] = 1;// (17) = 1; "After pair"
			EGS4.iausfl[17] = 1;// (18) = 1; "Before Compton"
			EGS4.iausfl[18] = 1;// (19) = 1; "After Compton"
			EGS4.iausfl[19] = 1;// (20) = 1; "Before photoelectric"
			EGS4.iausfl[20] = 1;// (21) = 1; "After photoelectric"
			EGS4.iausfl[23] = 1;// (24) = 1; "Before Rayleigh"
			EGS4.iausfl[24] = 1;// (25) = 1; "After Rayleigh"

			// "AUSGAB will be responsible for throwing away any photons resulting"
			// "from a primary electron. ie. True equilibtrium requires that all"
			// "energy deposition be local to the primary interaction site."
			EGS4.iausfl[7] = 1;// (8) = 1; "After bremsstrahlung"
			EGS4.iausfl[13] = 1;// (14) = 1;
								// "A positron has annihilated in-flight"
			EGS4.iausfl[14] = 1;// (15) = 1;
								// "A positron has annihilated at rest"
		}

		// "HATCH CALL PREPARATION AND EXECUTION"
		// System.out.println("before h");
		EGS4.DUNIT = 1; // "SET LENGTH UNITS TO CMS"
		EGS4.iprint = 0;
		HATCH();
		if (EGS4.STOPPROGRAM) {
			closeFile();
			destroyArrays();
			return;
		}
		// System.out.println("after h");
		EGS4.iprint = IPRINT;
		if (EGS4SrcEns.MONOEN == 0 && EGS4SrcEns.ISOURC != 21
				&& EGS4SrcEns.ISOURC != 22 && EGS4SrcEns.ISOURC != 23) {// "MONOENERGETIC INPUT BEAM"
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
		}
		// else { EKMAX=0; }// " <------------ fixme"

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
			destroyArrays();
			return;
		}

		DELTA = EGS4.AE[0] - EGS4.PRM + 1.E-7;
		;// AE(1) - PRM+1.D-7;
		DELTAT = EGS4.AE[0] + 1.E-7;// AE(1)+1.D-7;
		DELTLG = Math.log(DELTA);// "all passed in SCORE"
		// "The 1E-7 is needed to avoid roundoff problems since AE is single prec"
		// "CALCULATE THE MASS OF EACH ZONE (AREAL MASS FOR ISOURC=2 OR 4)
		// ==============================================
		for (int IZ = 1; IZ <= EGS4Geom.NZ; IZ++)// DO IZ=1,NZ
		{
			DEPTH = EGS4Geom.ZPLANE[IZ] - EGS4Geom.ZPLANE[IZ - 1];// DEPTH=ZPLANE(IZ+1)-ZPLANE(IZ);
			for (int IX = 1; IX <= EGS4Geom.NR; IX++)// DO IX=1,NR [
			{
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

			}
		}
		// ==================================================================================
		// "CALCULATE ONE-TIME-ONLY CONSTANTS FOR SOURCE"
		if (SOURCE == -1)
			EGS4SrcEns.SRCOTO();// (WEIGHT);

		if ((EGS4SrcEns.IFPB == 0) && (EGS4Macro.IFORCE != 0)
				&& (EGS4SrcEns.iqin == 0) && (EGS4SrcEns.MONOEN == 0)) {
			EGS4Macro.SELECT_MEAN_FREE_PATHS_FOR_FRONTAL_PARALLEL_BEAM();
			// -------------
			EGS4Macro.do_fast_step = true;
		}
		// -------------
		else {
			EGS4Macro.do_fast_step = false;
		}
		// =============================
		// ECUTMN = 1.E10;
		// DO I=2,NREG [
		// IF(MED(I)~= 0)[
		// ECUTMN=MIN(ECUTMN,ECUT(I));
		// ]
		// } //"NEED MINIMUM ECUT FOR PRESTA"
		// EK0=EIN; //"NOTE EIN IS SET AS MAXIMUM ENERGY FOR A SPECTRUM TOO"
		// $PRESTA-INPUTS;=>REPLACE {$PRESTA-INPUTS;} WITH {;}

		// =============================
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

		EK0 = EGS4SrcEns.ein;

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

		// "Initialize IWATCH routine"
		if (IWATCH != 0)
			EGS4.WATCH(-99, IWATCH);

		// "Output batches. Statistical analysis is done after each batch.
		// Execution
		// "stops if the desired statistical accuracy is obtained or there is
		// not enough
		// "time to do another batch.
		//boolean ENDSIM = false;
		EGS4SrcEns.AINFLU = 0.0;
		for (int IBATCH = 1; IBATCH <= $NBATCH; IBATCH++) {
			//ENDSIM = false;
			long startTime = System.currentTimeMillis();
			IBTCH = IBATCH;
			if (IBATCH == 1) {
				EGS4.seqStr = " BATCH" + EGS4.format("", 8) + "ELAPSED"
						+ EGS4.format("", 3) + "time of day";// +EGS4.format("",2)+"peak region stats(%)";
				if (EGS4.iprint > 1)
					printSequence(EGS4.seqStr);
			} else {// " not first batch"
			}// " end of before batch ne 1 block"

			for (int ICASE = 1; ICASE <= JCASE; ICASE++) {// "now fill each IS bin"

				if (EGS4SrcEns.ISOURC != 23)
					IHSTRY = IHSTRY + 1; // "increment history counter"

				EGS4Macro.NFTIME = 0; // "reset the photon forced interaction counter"
				if (SOURCE == -1)
					EGS4SrcEns.SRCHST();
				// "calculate the initial energy if a distribution is to be used"
				if (EGS4SrcEns.MONOEN != 0 && EGS4SrcEns.ISOURC != 21
						&& EGS4SrcEns.ISOURC != 22 && EGS4SrcEns.ISOURC != 23) {// "if equal to 0, it is monoenergetic"
					EGS4SrcEns.ein = EGS4SrcEns.ENSRCH();
					if (EGS4SrcEns.iqin == 0) {
						EI = EGS4SrcEns.ein;
					} else {
						EI = EGS4SrcEns.ein + EGS4.RM;
					}// "total energy"
				}

				// IFPB always 1!!!
				if (EGS4Macro.IFORCE == 1 && EGS4SrcEns.IFPB == 0
						&& EGS4SrcEns.MONOEN == 0 && EGS4SrcEns.iqin == 0) {
					// EGS4Macro.MAFORC=1;//flag for fast calculations in photon
					// forcing CMa"
					// GWAIT(1)=GWATE((IRIN-2)/NZ+1);WEIGHT=GWAIT(1);
					int IX = (EGS4SrcEns.irin - 2) / EGS4Geom.NZ + 1;
					EGS4Macro.GWAIT = EGS4Macro.GWATE[IX - 1];
					EGS4SrcEns.WEIGHT = EGS4Macro.GWAIT;
				}
				// else{EGS4Macro.MAFORC=0;}//"set the flag to 0 so that we won't do it for e-,e+ beams"]

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
				}
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
					destroyArrays();
					return;
				}

				if (IWATCH > 0)
					EGS4.WATCH(-1, IWATCH);

			}// "END OF THE ICASE LOOP"
			/*
			 * ============NEVER IF(IDAT = 0)[
			 * 
			 * ]"END OF CONDITIONAL DATA STORAGE"
			 * 
			 * $SET_ELAPSED_CPUTIME(CPUT2);
			 * TIMCPU=$CONVERSION_TO_SECONDS*(CPUT2-CPUT1)+TMCPUO; IF(IDAT = 0)[
			 * $PUT RNG STATE ON UNIT data_unit; WRITE(data_unit,*)
			 * IHSTRY,TIMCPU,NNREAD,PIISTP+count_pII_steps; ]
			 */
			String timePerBatch = EGS4.timeElapsedShort(startTime);
			Calendar call = Calendar.getInstance();
			String timeday = call.get(Calendar.HOUR) + ":" + call.get(Calendar.MINUTE)
					+ ":" + call.get(Calendar.SECOND);

			EGS4.seqStr = EGS4.format("", 2) + EGS4.format(IBATCH, 3)
					+ EGS4.format("", 2) + EGS4.format(timePerBatch, 14)
					+ EGS4.format("", 3) + EGS4.format(timeday, 9);// +
			// EGS4.format("",6)+EGS4.format(TDOS2,7,true);
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);

		}// "END OF SIMULATIONS"

		// "add unscored portions of _TMP arrays"

		SCSTP = SCSTP + SCSTP_TMP;
		SCSTP2 = SCSTP2 + SCSTP_TMP * SCSTP_TMP;
		for (int IX = 1; IX <= EGS4Geom.NR; IX++) {
			for (int IZ = 1; IZ <= EGS4Geom.NZ; IZ++) {
				for (int IT = 1; IT <= $MAXIT; IT++) {
					SCEDPN[IZ - 1][IX - 1][IT - 1] = SCEDPN[IZ - 1][IX - 1][IT - 1]
							+ SCEDPN_TMP[IZ - 1][IX - 1][IT - 1];
					SCEDPN2[IZ - 1][IX - 1][IT - 1] = SCEDPN2[IZ - 1][IX - 1][IT - 1]
							+ SCEDPN_TMP[IZ - 1][IX - 1][IT - 1]
							* SCEDPN_TMP[IZ - 1][IX - 1][IT - 1];
				}
				SCDOSE[IZ - 1][IX - 1][0] = SCDOSE[IZ - 1][IX - 1][0]
						+ SCEDPN_TMP[IZ - 1][IX - 1][0]
						- SCEDPN_TMP[IZ - 1][IX - 1][2];
				SCDOSE2[IZ - 1][IX - 1][0] = SCDOSE2[IZ - 1][IX - 1][0]
						+ (SCEDPN_TMP[IZ - 1][IX - 1][0] - SCEDPN_TMP[IZ - 1][IX - 1][2])
						* (SCEDPN_TMP[IZ - 1][IX - 1][0] - SCEDPN_TMP[IZ - 1][IX - 1][2]);
				SCDOSE_COV[IZ - 1][IX - 1][0] = SCDOSE_COV[IZ - 1][IX - 1][0]
						+ (SCEDPN_TMP[IZ - 1][IX - 1][0] - SCEDPN_TMP[IZ - 1][IX - 1][2])
						* SCDOSE_TMP3[IZ - 1][IX - 1];

				SCDOSE[IZ - 1][IX - 1][1] = SCDOSE[IZ - 1][IX - 1][1]
						+ SCDOSE_TMP3[IZ - 1][IX - 1]
						- SCEDPN_TMP[IZ - 1][IX - 1][0];
				SCDOSE2[IZ - 1][IX - 1][1] = SCDOSE2[IZ - 1][IX - 1][1]
						+ (SCDOSE_TMP3[IZ - 1][IX - 1] - SCEDPN_TMP[IZ - 1][IX - 1][0])
						* (SCDOSE_TMP3[IZ - 1][IX - 1] - SCEDPN_TMP[IZ - 1][IX - 1][0]);
				SCDOSE_COV[IZ - 1][IX - 1][1] = SCDOSE_COV[IZ - 1][IX - 1][1]
						+ (SCDOSE_TMP3[IZ - 1][IX - 1] - SCEDPN_TMP[IZ - 1][IX - 1][0])
						* SCDOSE_TMP3[IZ - 1][IX - 1];

				SCDOSE[IZ - 1][IX - 1][2] = SCDOSE[IZ - 1][IX - 1][2]
						+ SCDOSE_TMP3[IZ - 1][IX - 1];
				SCDOSE2[IZ - 1][IX - 1][2] = SCDOSE2[IZ - 1][IX - 1][2]
						+ SCDOSE_TMP3[IZ - 1][IX - 1]
						* SCDOSE_TMP3[IZ - 1][IX - 1];
				SCSPR_COV[IZ - 1][IX - 1][0] = SCSPR_COV[IZ - 1][IX - 1][0]
						+ SCEDPN_TMP[IZ - 1][IX - 1][0]
						* SCEDPN_TMP[IZ - 1][IX - 1][1];
				SCSPR_COV[IZ - 1][IX - 1][1] = SCSPR_COV[IZ - 1][IX - 1][1]
						+ SCEDPN_TMP[IZ - 1][IX - 1][2]
						* SCEDPN_TMP[IZ - 1][IX - 1][3];
				if (SCDOSE_TMP3[IZ - 1][IX - 1] > 0.)
					SCNUM[IZ - 1][IX - 1] = SCNUM[IZ - 1][IX - 1] + 1;

			}
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
		if (EGS4SrcEns.ISOURC == 21 || EGS4SrcEns.ISOURC == 22) {
			// SCORE_NORM_NUM=EGS4SrcEns.dble(NNREAD+NRCYCL*(NNREAD-IHSTRY))/
			// EGS4SrcEns.dble(NCASE_PHSP)*NINCSRC;
		} else {
			SCORE_NORM_NUM = EGS4SrcEns.dble(IHSTRY);
		}

		// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
		// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

		// "STATISTICAL ANALYSES ON THE RAW DATA and do normalizations"
		/*
		 * REPLACE {$ANALYZE(#,#:#)} WITH {;
		 * 
		 * "Macro to analyze uncertainty:"
		 * "{P1}{P2}=scoring array (eg SCEDPN(IZ,IX,IT))"
		 * "{P3}=quantity to normalize by (eg incident no. of particles)"
		 * 
		 * "Calculates the uncertainty on {P1}{P2}/{P3}.  The "
		 * "uncertainty is stored in {P1}2{P2} and is expressed as a percentage of"
		 * "{P1}{P2}/{P3} (max 99.9%).  Note that you must define the REAL variable"
		 * "SCORE_TEMP in any subroutine where this macro is used.  This macro"
		 * "is only used in the analysis of no. of steps."
		 * 
		 * SCORE_TEMP={P1}{P2}/{P3}; {P1}2{P2}={P1}2{P2}/{P3};
		 * {P1}2{P2}=({P1}2{P2}-SCORE_TEMP*SCORE_TEMP)/({P3}-1);
		 * IF({P1}2{P2}>=0.) {P1}2{P2}= SQRT({P1}2{P2}); IF(SCORE_TEMP~=0.)[
		 * {P1}2{P2}={P1}2{P2}/SCORE_TEMP*100; IF( {P1}2{P2} > 99.9 ) [
		 * {P1}2{P2} = 99.9; ] ] ELSE[ {P1}2{P2}=99.9; ] ; }
		 */

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

		imoi = REGSVOL[0];// where to score IRL of interest!!
		// EGS4Core:=>EGS4.MEDIUM=EGS4.MED[irl-1];//eg. med[1]=1 mediul 1
		// imoi==regoin in wich we have medium of interest
		imoi = EGS4.MED[imoi - 1];// medium of interest!!!
		// System.out.println("imoi= "+imoi+" media "+EGS4.MEDIA[imoi-1]);
		for (int IX = 1; IX <= EGS4Geom.NR; IX++) {
			for (int IZ = 1; IZ <= EGS4Geom.NZ; IZ++) {
				// $GET-IRL(IZ,IX); MEDNUM = MED(IRL);
				IRL = EGS4Geom.GET_IRL(IZ, IX);
				MEDNUM = EGS4.MED[IRL - 1];
				if (SCNUM[IZ - 1][IX - 1] < $GOODSTATS) {
					// "mark this as bad statistics"
					BSTAT[IZ - 1][0] = -1;
					BSTAT[IZ - 1][1] = -1;
				}// ]

				if (SCEDPN[IZ - 1][IX - 1][1] == 0.) {
					SCEDPN[IZ - 1][IX - 1][0] = 0.;
					SCEDPN2[IZ - 1][IX - 1][0] = 99.9;
					BSTAT[IZ - 1][0] = -1; // "FLAG FOR BAD STATISTICS"
				} else {// " CALCULATE STOPPING POWER RATIO FOR REGULAR EDPN"
						// $ANALYZE(SCEDPN,(IZ,IX,1):SCORE_NORM_NUM);
					SCORE_TEMP = SCEDPN[IZ - 1][IX - 1][0] / SCORE_NORM_NUM;
					SCEDPN2[IZ - 1][IX - 1][0] = SCEDPN2[IZ - 1][IX - 1][0]
							/ SCORE_NORM_NUM;
					SCEDPN2[IZ - 1][IX - 1][0] = (SCEDPN2[IZ - 1][IX - 1][0] - SCORE_TEMP
							* SCORE_TEMP)
							/ (SCORE_NORM_NUM - 1);
					if (SCEDPN2[IZ - 1][IX - 1][0] >= 0.)
						SCEDPN2[IZ - 1][IX - 1][0] = Math
								.sqrt(SCEDPN2[IZ - 1][IX - 1][0]);
					if (SCORE_TEMP != 0.) {
						SCEDPN2[IZ - 1][IX - 1][0] = Math.min(
								SCEDPN2[IZ - 1][IX - 1][0] / SCORE_TEMP * 100.,
								99.9);
					} else {
						SCEDPN2[IZ - 1][IX - 1][0] = 99.9;
					}

					// $ANALYZE(SCEDPN,(IZ,IX,2):SCORE_NORM_NUM);
					SCORE_TEMP = SCEDPN[IZ - 1][IX - 1][1] / SCORE_NORM_NUM;
					SCEDPN2[IZ - 1][IX - 1][1] = SCEDPN2[IZ - 1][IX - 1][1]
							/ SCORE_NORM_NUM;
					SCEDPN2[IZ - 1][IX - 1][1] = (SCEDPN2[IZ - 1][IX - 1][1] - SCORE_TEMP
							* SCORE_TEMP)
							/ (SCORE_NORM_NUM - 1);
					if (SCEDPN2[IZ - 1][IX - 1][1] >= 0.)
						SCEDPN2[IZ - 1][IX - 1][1] = Math
								.sqrt(SCEDPN2[IZ - 1][IX - 1][1]);
					if (SCORE_TEMP != 0.) {
						SCEDPN2[IZ - 1][IX - 1][1] = Math.min(
								SCEDPN2[IZ - 1][IX - 1][1] / SCORE_TEMP * 100.,
								99.9);
					} else {
						SCEDPN2[IZ - 1][IX - 1][1] = 99.9;
					}

					if (SCEDPN2[IZ - 1][IX - 1][0] == 99.9
							|| SCEDPN2[IZ - 1][IX - 1][1] == 99.9) {
						BSTAT[IZ - 1][0] = -1;
						SCSPR_COV[IZ - 1][IX - 1][0] = 0.;
					}// ]
					else {
						SCSPR_COV[IZ - 1][IX - 1][0] = SCSPR_COV[IZ - 1][IX - 1][0]
								/ SCORE_NORM_NUM
								- SCEDPN[IZ - 1][IX - 1][0]
								* SCEDPN[IZ - 1][IX - 1][1]
								/ (SCORE_NORM_NUM * SCORE_NORM_NUM);

						SCSPR_COV[IZ - 1][IX - 1][0] = SCSPR_COV[IZ - 1][IX - 1][0]
								/ (SCEDPN[IZ - 1][IX - 1][0]
										* SCEDPN[IZ - 1][IX - 1][1]
										/ (SCORE_NORM_NUM * SCORE_NORM_NUM) * (SCORE_NORM_NUM - 1));
					}
					SCEDPN2[IZ - 1][IX - 1][0] = (SCEDPN2[IZ - 1][IX - 1][0] / 100.)
							* (SCEDPN2[IZ - 1][IX - 1][0] / 100.)
							+ (SCEDPN2[IZ - 1][IX - 1][1] / 100.)
							* (SCEDPN2[IZ - 1][IX - 1][1] / 100.)
							- 2.0
							* SCSPR_COV[IZ - 1][IX - 1][0];
					if (SCEDPN2[IZ - 1][IX - 1][0] > 0.) {
						SCEDPN2[IZ - 1][IX - 1][0] = Math
								.sqrt(SCEDPN2[IZ - 1][IX - 1][0]) * 100;
					}
					if (SCEDPN2[IZ - 1][IX - 1][0] > 99.9)
						SCEDPN2[IZ - 1][IX - 1][0] = 99.9;
					SCEDPN[IZ - 1][IX - 1][0] = (SCEDPN[IZ - 1][IX - 1][0] / EGS4.RHOR[IRL - 1])
							/ (SCEDPN[IZ - 1][IX - 1][1] / EGS4.RHO[imoi - 1]);// /RHO(2));
					// public static double[] RHOR=new double[$MXREG];
					// public static double[] RHO=new double[$MXMED];!!!!!!!!
				}
				if (SCEDPN[IZ - 1][IX - 1][3] == 0.) {// [
					SCEDPN[IZ - 1][IX - 1][2] = 0.;
					SCEDPN2[IZ - 1][IX - 1][2] = 99.9;
					BSTAT[IZ - 1][1] = -1; // "FLAG FOR BAD STATISTICS"
				}// ]
				else {// [" CALCULATE STOPPING POWER RATIO FOR REGULAR EDPN"
						// $ANALYZE(SCEDPN,(IZ,IX,3):SCORE_NORM_NUM);
					SCORE_TEMP = SCEDPN[IZ - 1][IX - 1][2] / SCORE_NORM_NUM;
					SCEDPN2[IZ - 1][IX - 1][2] = SCEDPN2[IZ - 1][IX - 1][2]
							/ SCORE_NORM_NUM;
					SCEDPN2[IZ - 1][IX - 1][2] = (SCEDPN2[IZ - 1][IX - 1][2] - SCORE_TEMP
							* SCORE_TEMP)
							/ (SCORE_NORM_NUM - 1);
					if (SCEDPN2[IZ - 1][IX - 1][2] >= 0.)
						SCEDPN2[IZ - 1][IX - 1][2] = Math
								.sqrt(SCEDPN2[IZ - 1][IX - 1][2]);
					if (SCORE_TEMP != 0.) {
						SCEDPN2[IZ - 1][IX - 1][2] = Math.min(
								SCEDPN2[IZ - 1][IX - 1][2] / SCORE_TEMP * 100.,
								99.9);
					} else {
						SCEDPN2[IZ - 1][IX - 1][2] = 99.9;
					}

					// $ANALYZE(SCEDPN,(IZ,IX,4):SCORE_NORM_NUM);
					SCORE_TEMP = SCEDPN[IZ - 1][IX - 1][3] / SCORE_NORM_NUM;
					SCEDPN2[IZ - 1][IX - 1][3] = SCEDPN2[IZ - 1][IX - 1][3]
							/ SCORE_NORM_NUM;
					SCEDPN2[IZ - 1][IX - 1][3] = (SCEDPN2[IZ - 1][IX - 1][3] - SCORE_TEMP
							* SCORE_TEMP)
							/ (SCORE_NORM_NUM - 1);
					if (SCEDPN2[IZ - 1][IX - 1][3] >= 0.)
						SCEDPN2[IZ - 1][IX - 1][3] = Math
								.sqrt(SCEDPN2[IZ - 1][IX - 1][3]);
					if (SCORE_TEMP != 0.) {
						SCEDPN2[IZ - 1][IX - 1][3] = Math.min(
								SCEDPN2[IZ - 1][IX - 1][3] / SCORE_TEMP * 100.,
								99.9);
					} else {
						SCEDPN2[IZ - 1][IX - 1][3] = 99.9;
					}

					if (SCEDPN2[IZ - 1][IX - 1][2] == 99.9
							|| SCEDPN2[IZ - 1][IX - 1][3] == 99.9) {
						BSTAT[IZ - 1][1] = -1;
						SCSPR_COV[IZ - 1][IX - 1][1] = 0.;
					}// ]
					else {// [
						SCSPR_COV[IZ - 1][IX - 1][1] = SCSPR_COV[IZ - 1][IX - 1][1]
								/ SCORE_NORM_NUM
								- SCEDPN[IZ - 1][IX - 1][2]
								* SCEDPN[IZ - 1][IX - 1][3]
								/ (SCORE_NORM_NUM * SCORE_NORM_NUM);
						SCSPR_COV[IZ - 1][IX - 1][1] = SCSPR_COV[IZ - 1][IX - 1][1]
								/ (SCEDPN[IZ - 1][IX - 1][2]
										* SCEDPN[IZ - 1][IX - 1][3]
										/ (SCORE_NORM_NUM * SCORE_NORM_NUM) * (SCORE_NORM_NUM - 1));
					}// ]
					SCEDPN2[IZ - 1][IX - 1][2] = (SCEDPN2[IZ - 1][IX - 1][2] / 100.)
							* (SCEDPN2[IZ - 1][IX - 1][2] / 100.)
							+ (SCEDPN2[IZ - 1][IX - 1][3] / 100.)
							* (SCEDPN2[IZ - 1][IX - 1][3] / 100.)
							- 2.0
							* SCSPR_COV[IZ - 1][IX - 1][1];
					if (SCEDPN2[IZ - 1][IX - 1][2] > 0.) {
						SCEDPN2[IZ - 1][IX - 1][2] = Math
								.sqrt(SCEDPN2[IZ - 1][IX - 1][2]) * 100;
					}// ]
					if (SCEDPN2[IZ - 1][IX - 1][2] > 99.9)
						SCEDPN2[IZ - 1][IX - 1][2] = 99.9;
					SCEDPN[IZ - 1][IX - 1][2] = (SCEDPN[IZ - 1][IX - 1][2] / EGS4.RHOR[IRL - 1])
							/ (SCEDPN[IZ - 1][IX - 1][3] / EGS4.RHO[imoi - 1]);// RHO(2));
				}// ]
					// "analyze uncertainty on doses and normalize doses"
					// $ANALYZE(SCDOSE,(IZ,IX,1):SCORE_NORM_NUM);
				SCORE_TEMP = SCDOSE[IZ - 1][IX - 1][0] / SCORE_NORM_NUM;
				SCDOSE2[IZ - 1][IX - 1][0] = SCDOSE2[IZ - 1][IX - 1][0]
						/ SCORE_NORM_NUM;
				SCDOSE2[IZ - 1][IX - 1][0] = (SCDOSE2[IZ - 1][IX - 1][0] - SCORE_TEMP
						* SCORE_TEMP)
						/ (SCORE_NORM_NUM - 1);
				if (SCDOSE2[IZ - 1][IX - 1][0] >= 0.)
					SCDOSE2[IZ - 1][IX - 1][0] = Math
							.sqrt(SCDOSE2[IZ - 1][IX - 1][0]);
				if (SCORE_TEMP != 0.) {
					SCDOSE2[IZ - 1][IX - 1][0] = Math.min(
							SCDOSE2[IZ - 1][IX - 1][0] / SCORE_TEMP * 100.,
							99.9);
				} else {
					SCDOSE2[IZ - 1][IX - 1][0] = 99.9;
				}

				// $ANALYZE(SCDOSE,(IZ,IX,2):SCORE_NORM_NUM);
				SCORE_TEMP = SCDOSE[IZ - 1][IX - 1][1] / SCORE_NORM_NUM;
				SCDOSE2[IZ - 1][IX - 1][1] = SCDOSE2[IZ - 1][IX - 1][1]
						/ SCORE_NORM_NUM;
				SCDOSE2[IZ - 1][IX - 1][1] = (SCDOSE2[IZ - 1][IX - 1][1] - SCORE_TEMP
						* SCORE_TEMP)
						/ (SCORE_NORM_NUM - 1);
				if (SCDOSE2[IZ - 1][IX - 1][1] >= 0.)
					SCDOSE2[IZ - 1][IX - 1][1] = Math
							.sqrt(SCDOSE2[IZ - 1][IX - 1][1]);
				if (SCORE_TEMP != 0.) {
					SCDOSE2[IZ - 1][IX - 1][1] = Math.min(
							SCDOSE2[IZ - 1][IX - 1][1] / SCORE_TEMP * 100.,
							99.9);
				} else {
					SCDOSE2[IZ - 1][IX - 1][1] = 99.9;
				}

				// $ANALYZE(SCDOSE,(IZ,IX,3):SCORE_NORM_NUM);
				SCORE_TEMP = SCDOSE[IZ - 1][IX - 1][2] / SCORE_NORM_NUM;
				SCDOSE2[IZ - 1][IX - 1][2] = SCDOSE2[IZ - 1][IX - 1][2]
						/ SCORE_NORM_NUM;
				SCDOSE2[IZ - 1][IX - 1][2] = (SCDOSE2[IZ - 1][IX - 1][2] - SCORE_TEMP
						* SCORE_TEMP)
						/ (SCORE_NORM_NUM - 1);
				if (SCDOSE2[IZ - 1][IX - 1][2] >= 0.)
					SCDOSE2[IZ - 1][IX - 1][2] = Math
							.sqrt(SCDOSE2[IZ - 1][IX - 1][2]);
				if (SCORE_TEMP != 0.) {
					SCDOSE2[IZ - 1][IX - 1][2] = Math.min(
							SCDOSE2[IZ - 1][IX - 1][2] / SCORE_TEMP * 100.,
							99.9);
				} else {
					SCDOSE2[IZ - 1][IX - 1][2] = 99.9;
				}

				if (SCDOSE2[IZ - 1][IX - 1][0] == 99.9
						|| SCDOSE2[IZ - 1][IX - 1][2] == 99.9) {
					SCDOSE_COV[IZ - 1][IX - 1][0] = 0.0;
					// if(BSTAT[0][IZ-1]!=-1)BDSTAT[IZ-1]=-1;//
					// "only set if BSTAT(1,IZ)"
					// "is not already set"
					if (BSTAT[IZ - 1][0] != -1)
						BDSTAT[IZ - 1] = -1;
				}// ]
				else {
					SCDOSE_COV[IZ - 1][IX - 1][0] = SCDOSE_COV[IZ - 1][IX - 1][0]
							/ SCORE_NORM_NUM
							- SCDOSE[IZ - 1][IX - 1][0]
							* SCDOSE[IZ - 1][IX - 1][2]
							/ (SCORE_NORM_NUM * SCORE_NORM_NUM);
					SCDOSE_COV[IZ - 1][IX - 1][0] = SCDOSE_COV[IZ - 1][IX - 1][0]
							/ (SCDOSE[IZ - 1][IX - 1][0]
									* SCDOSE[IZ - 1][IX - 1][2]
									/ (SCORE_NORM_NUM * SCORE_NORM_NUM) * (SCORE_NORM_NUM - 1));
				}
				SCDOSE2[IZ - 1][IX - 1][0] = (SCDOSE2[IZ - 1][IX - 1][0] / 100.)
						* (SCDOSE2[IZ - 1][IX - 1][0] / 100.)
						+ (SCDOSE2[IZ - 1][IX - 1][2] / 100.)
						* (SCDOSE2[IZ - 1][IX - 1][2] / 100.)
						- 2.0
						* SCDOSE_COV[IZ - 1][IX - 1][0];

				if (SCDOSE2[IZ - 1][IX - 1][0] > 0.) {
					SCDOSE2[IZ - 1][IX - 1][0] = Math
							.sqrt(SCDOSE2[IZ - 1][IX - 1][0]) * 100;
				}
				if (SCDOSE2[IZ - 1][IX - 1][0] > 99.9)
					SCDOSE2[IZ - 1][IX - 1][0] = 99.9;

				if (SCDOSE[IZ - 1][IX - 1][2] != 0.) {
					SCDOSE[IZ - 1][IX - 1][0] = SCDOSE[IZ - 1][IX - 1][0]
							/ SCDOSE[IZ - 1][IX - 1][2];
				}

				if (SCDOSE2[IZ - 1][IX - 1][1] == 99.9
						|| SCDOSE2[IZ - 1][IX - 1][2] == 99.9) {
					SCDOSE_COV[IZ - 1][IX - 1][1] = 0.0;
					// if(BSTAT[0][IZ-1]!=-1)BDSTAT[IZ-1]=-1;
					if (BSTAT[IZ - 1][0] != -1)
						BDSTAT[IZ - 1] = -1;
				} else {
					SCDOSE_COV[IZ - 1][IX - 1][1] = SCDOSE_COV[IZ - 1][IX - 1][1]
							/ SCORE_NORM_NUM
							- SCDOSE[IZ - 1][IX - 1][1]
							* SCDOSE[IZ - 1][IX - 1][2]
							/ (SCORE_NORM_NUM * SCORE_NORM_NUM);
					SCDOSE_COV[IZ - 1][IX - 1][1] = SCDOSE_COV[IZ - 1][IX - 1][1]
							/ (SCDOSE[IZ - 1][IX - 1][1]
									* SCDOSE[IZ - 1][IX - 1][2]
									/ (SCORE_NORM_NUM * SCORE_NORM_NUM) * (SCORE_NORM_NUM - 1));
				}
				SCDOSE2[IZ - 1][IX - 1][1] = (SCDOSE2[IZ - 1][IX - 1][1] / 100.)
						* (SCDOSE2[IZ - 1][IX - 1][1] / 100.)
						+ (SCDOSE2[IZ - 1][IX - 1][2] / 100.)
						* (SCDOSE2[IZ - 1][IX - 1][2] / 100.)
						- 2.0
						* SCDOSE_COV[IZ - 1][IX - 1][1];
				if (SCDOSE2[IZ - 1][IX - 1][1] > 0.) {
					SCDOSE2[IZ - 1][IX - 1][1] = Math
							.sqrt(SCDOSE2[IZ - 1][IX - 1][1]) * 100;
				}
				if (SCDOSE2[IZ - 1][IX - 1][1] > 99.9)
					SCDOSE2[IZ - 1][IX - 1][1] = 99.9;
				if (SCDOSE[IZ - 1][IX - 1][2] != 0.) {
					SCDOSE[IZ - 1][IX - 1][1] = SCDOSE[IZ - 1][IX - 1][1]
							/ SCDOSE[IZ - 1][IX - 1][2];
				}

				// "now, normalize total dose"
				if (SCDOSE[IZ - 1][IX - 1][2] != 0.) {
					FMASS = AMASS[IZ - 1][IX - 1];
					if (FMASS == 0.0)
						FMASS = 1.0; // "AVOIDS /0 FOR VACUUM"
					SCDOSE[IZ - 1][IX - 1][2] = SCDOSE[IZ - 1][IX - 1][2]
							* 1.602e-10 / (FMASS * EGS4SrcEns.AINFLU);
				}// ]"end non-zero block"

			}
		}

		OSUMRY();// CALL OSUMRY; "PRINT THE OUTPUT SUMRY"
		// finalizations
		destroyArrays();
		// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
		// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
		/*
		 * IF(ITRACKE ~= 0) [WRITE(IOUT,:TEMP:) dble(ICROSS)/dble(ITRACKE);]
		 * :TEMP: FORMAT(///' Fraction of tracks treated as trackends which
		 * were', /T20,' really steps crossing DELTA =', F10.3/' Usually close
		 * to 0.5 because', ' each crosser turns into a stopper'/);
		 * :END-OF-RUN:;
		 */
		// :END-OF-RUN:;

		// "******************************************************************************
		// "
		// " *** SECTION 4 ***
		// "
		// "------------------------------------------------------------------------------
		// "
		// "THE CONCLUSION"
		// "
		// "------------------------------------------------------------------------------

		// :END:;
		// " ******* IK: using wrapper subroutines for date and time.
		// " Was:
		// "CALL DATE(DATEN);CALL TIME(TIMEN);
		// "OUTPUT DATEN,TIMEN;(/' END OF FILE
		// READ',10X,$DATEN_FORMAT,1X,$TIMEN_FORMAT//);
		// "OUTPUT DATEN,TIMEN;(/' END OF
		// RUN',10X,$DATEN_FORMAT,1X,$TIMEN_FORMAT//);
		// OUTPUT; (/' END OF RUN',9X,' ',$); call egs_fdate(6);
		// write(6,'(//)');

		// IF (ERROR_FLAG=1) [
		// OUTPUT ;(//' ********** ERRORS AND WARNING MESSAGES ********** '//);
		// ]

		// REWIND (15);
		// 9998 CONTINUE;
		// READ (15, END=9999, FMT='(A80)') TEXT;
		// WRITE (*,'(A79)') TEXT;
		// GOTO 9998;
		// 9999 CONTINUE;
		// " ******* IK: using wrapper subroutines for date and time.
		// " Was:
		// " WRITE(15 , 400) DATEN,TIMEN;
		// " WRITE(IOUT,400) DATEN,TIMEN;
		// write(15,400); call egs_fdate(15); write(15,*);
		// write(iout,400); call egs_fdate(iout); write(iout,*);

		// call egs_finish;

		// #ifdef HAVE_C_COMPILER;
		// ;
		// IF( n_parallel > 0 & ~is_finished ) [
		// call egs_pjob_finish(n_job);
		// IF( n_job = 0 ) [
		// is_finished = .true.;
		// call egs_combine_runs(combine_results,'.egsdat');
		// NCASET=NCASEO; IHSTRY=NCASET;
		// CALL SRCOTO(WEIGHT);
		// goto :STATS-ANAL:;
		// ]
		// ]
		// #endif;

		// call SRCEND;

		// $CALL_EXIT(0);
		/*
		 * 200 FORMAT(' ',79('*')// T20,'EXECUTION INFORMATION AND WARNING
		 * MESSAGES'// 1x,79('*')//' Using FLURZnrc(EGSnrc) '$VERSION' ') ;
		 * 
		 * " ******* IK: using wrapper subroutines for date and time. " Was: 400
		 * FORMAT(/' End of run',10X,$DATEN_FORMAT,1X,$TIMEN_FORMAT/); 400
		 * FORMAT(/' End of run',10X,$);
		 * 
		 * "FORTRAN format statements. format statement N## is first used in section N."
		 * %F %I0 201 FORMAT(/' ********* NEW INPUT FILE *********'/) 202
		 * FORMAT(/' ********* RESTARTED INPUT FILE ********* '/ 1' ',10X,I12,'
		 * NEW + ',I12,' OLD HISTORIES'/ 2' ',20X,'NEW RN SEEDS=',$) 204
		 * FORMAT(/' ********* DATA ANALYSIS ONLY *********'/) 205 FORMAT(/'
		 * ********* POST-PROCESSING PARALLEL RUNS **********'/) 206
		 * FORMAT(/'********* RANDOM NUMBERS READ FROM UNIT 2 *********'/) 210
		 * FORMAT(/' ********* NOT ENOUGH TIME to finish within', 1' limit
		 * of',F8.2,' hours',I5,' batches used********',/10X, 2' ',I12,'
		 * Histories run, ',I12,' Histories analyzed'//) 220 FORMAT(I12,' ',I12)
		 * 250 FORMAT(/' ********* FINAL RANDOM NUMBER SEED: ',2I12,'
		 * *********') 255 FORMAT(/' FOR OLD RUN:'/ 1' ----------- '/ 2' TOTAL
		 * CPUTIME =',F8.1,'s (=',F5.2,' hr)') 256 FORMAT(/' FOR PARALLEL
		 * RUNS:'/ 1' ----------------- '/ 2' ON ',I5,' MACHINES '/ 3' TOTAL
		 * CPUTIME =',F8.1,'s (=',F8.2,' hr), CPUTIME/MACHINE =', 4F8.1,'s') 260
		 * FORMAT(/'FINISHED: TIME ELSAPSED, CPUTIME', 1',RATIO=
		 * ',2F8.1,'(=',F5.2,'hr)',F8.2) 270 FORMAT(/' FINISHED SIMULATIONS:
		 * CPUTIME=',F8.1,'(=',F5.2,' hr)') 280 FORMAT(/' CPUTIME / history
		 * =',F10.5,' sec. # histories', 1' / hour =',F12.0/) 500 FORMAT(//'
		 * NORMALIZATION PROBLEMS IN FIRST BIN(/MeV) SINCE', X' DIFFERENT ECUTS
		 * ARE USED IN DIFFERENT REGIONS '/ X' MAXIMUM ECUT=',F10.3,' USED
		 * EVERYWHERE***'//) %M %I4
		 * 
		 * END; "end of main routine-FLURZnrc"
		 */

		EGS4SrcEns.SRCEND();// do nothing

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
	// "
	// "
	// "******************************************************************************
	/*
	 * " An AUSGAB routine to be used with sprrznrc.mortran " " Energy
	 * deposition is scored in SCEDPN(IZ,IX,IT) and doses are " scored in
	 * SCDOSE(IZ,IX,IT). " " For SCEDPN: " " IT=1 EDPN to local material
	 * excluding brem created below PCUT, atomic " relaxation and electrons
	 * created below ECUT (compton, pe or pair) " but includes electrons slowing
	 * across ECUT " IT=2 EDPN to second material with transport in first
	 * material " IT=3 EDPN to local material without any stoppers " IT=4 like
	 * IT=3 to 2nd material with transport in first material " IT=5 EDPN for
	 * everything in local material " " For SCDOSE: " " IT=1 Fraction of total
	 * dose due to stoppers including " electrons/positrons slowing down across
	 * delta. " IT=2 Fraction of total dose due to e- created below delta and
	 * photons " created below AP " IT=3 Total dose " " SPR IS FOR current
	 * MEDIUM TO MEDIUM 2. " " SPR(1) = EDPN(IT=1)/RHOR(IRL) / EDPN(IT=2)/RHO(2)
	 * " SPR(2) = EDPN(IT=3)/RHOR(IRL) / EDPN(IT=4)/RHO(2) i.e. no stoppers " "
	 * Note: SPR(1) is actually placed in SCEDPN(IZ,IX,1) and SPR(2) is " placed
	 * in SCEDPN(IZ,IX,2) at the end of the run. " " RHOs are needed since we
	 * use linear not mass stopping powers. " " VERSION 01 adapted from
	 * CAVITY(EGS3) AUSGAB 10/87 AFB " VERSION 02 major revision Jan 1990, DWOR
	 * " removed unused code March 1990 DR " VERSION 03 another revision with
	 * DELTA etc introduced " Mar 30, 1990 DR " " Variables: " DELTA set in
	 * main, ke of AE, comin SCORE " DELTLG log of DELTA " DELTAT AE = DELTA +
	 * RM " ERATIO(I) ratio of medium 2 to medium I stopping power at DELTA "
	 * set once for use with track ends - must save for F77 " ERATIO(2) = 1 once
	 * initialized, and used to test if set " PRATIO(I) same for positrons
	 * 
	 * OBS=>REPLACE {$MXMED} WITH {3}; //1 vacuum; 2 cyl; 3 exit
	 */
	/**
	 * In general, AUSGAB is a routine which is called under a series 
	 * of well defined conditions specified by the value of IARG. Interface method.
	 * @param IARG the process code
	 */
	public void AUSGAB(int IARG) {

		// $IMPLICIT-NONE;
		int IRL = 0;
		int IQL = 0;
		//int IBTOP = 0;
		int IX = 0;
		int IZ = 0;// int IB=0;int NBIN=0;
		int I = 0;
		//int LATCHL = 0;
		int MEDL = 0;
		int LELKE = 0;
		int ip = 0;
		int IFLAG = 0;
		//int IT = 0;

		double WTL = 0.0;
		double DEDXME = 0.0;
		double DEDXMP = 0.0;
		double DEDXMI = 0.0;
		double DEDXMJ = 0.0;
		double XSI = 0.0;
		double FTMP = 0.0;
		double EKEtmp = 0.0;
		double ELKEtmp = 0.0;
		double DEDXM1 = 0.0;
		double DEDXM2 = 0.0;
		double SRATIO = 0.0;
		double[] ERATIO = new double[EGS4.$MXMED];// ($MXMED),
		double[] PRATIO = new double[EGS4.$MXMED];// ($MXMED);

		// "STACK OVERFLOW CHECK"
		MXNP = Math.max(MXNP, EGS4.NP);// "keep track of how deep stack is"
		// "MXNP is not output but it should be"

		// "check if particle is leaving the transport geometry"
		IRL = EGS4.IR[EGS4.NP - 1]; // "local region number"
		if (IRL == 1) {// "outside the chamber IDISC will have been set in HOWFAR"
			if (IWATCH > 0)
				EGS4.WATCH(IARG, IWATCH);// CALL WATCH(IARG,IWATCH);
			return;
		}// ]
			// $GET-IX-IZ(IRL); "local plane and radius numbers"
		IX = EGS4Geom.GET_IX(IRL);
		IZ = EGS4Geom.GET_IZC(IRL);

		IQL = EGS4.IQ[EGS4.NP - 1]; // "local charge variable"
		WTL = EGS4.WT[EGS4.NP - 1]; // "local weight variable"
		//LATCHL = EGS4.LATCH[EGS4.NP - 1]; // "LATCHL=0 except = 99 for e created below delta"
		MEDL = EGS4.MED[IRL - 1];// "Current medium for this region - comin MISC"
		// ++++++++++++++++++++++++++++++++++++++++++++++
		int MEDIUM = 0;// @@@@@@@@@@
		imoi = REGSVOL[0];// where to score IRL of interest!!
		// EGS4Core:=>EGS4.MEDIUM=EGS4.MED[irl-1];//eg. med[1]=1 mediul 1
		imoi = EGS4.MED[imoi - 1];// medium of interest!!!
		// +++++++++++++++++++++++++++++++++++++++++++++++
		// System.out.println(ERATIO[imoi-1]+ "  ");
		if (ERATIO[imoi - 1] == 0.0)// (2) =
									// 0.0)["ratio of stopping powers at E=DELTA are "
		// "                   not initialized yet - do it now"
		{

			EKEtmp = EGS4.EKE;
			ELKEtmp = EGS4.ELKE;
			MEDIUM = imoi;// 2; //"set to second medium" OR MEDIUM OF INTEREST
							// RISKY TO USE EGS4.MEDIUM!!
			EGS4.EKE = DELTA;
			EGS4.ELKE = DELTLG; // "DELTLG is calculated in inputs and passed in SCORE"
			// $SET INTERVAL ELKE, EKE; "prepare to interpolate in  data tables"
			Double dbl = new Double(EGS4.EKE1[MEDIUM - 1] * EGS4.ELKE
					+ EGS4.EKE0[MEDIUM - 1]);
			LELKE = dbl.intValue();

			// $EVALUATE DEDXME USING EDEDX(ELKE); "sp for electrons in med2"
			// $EVALUATE dedx0 USING ededx(elke);
			DEDXME = EGS4.EDEDX1[LELKE - 1][MEDIUM - 1] * EGS4.ELKE
					+ EGS4.EDEDX0[LELKE - 1][MEDIUM - 1];
			// $EVALUATE DEDXMP USING PDEDX(ELKE); "sp for positrons in med2"
			DEDXMP = EGS4.PDEDX1[LELKE - 1][MEDIUM - 1] * EGS4.ELKE
					+ EGS4.PDEDX0[LELKE - 1][MEDIUM - 1];
			// write(6,*) ' Initializing track end terms ERATIO and PRATIO: ';
			for (I = 1; I <= EGS4.NMED; I++) {
				if (I == imoi)// (I=2)
				{
					// ERATIO(2),PRATIO(2)/ = 1.0; "from definition"
					ERATIO[imoi - 1] = 1.0;
					PRATIO[imoi - 1] = 1.0;
				} else {
					MEDIUM = I; // "set to current medium"
					EGS4.EKE = DELTA;
					EGS4.ELKE = DELTLG;
					// $SET INTERVAL ELKE, EKE; "prepare to interpolate tables"
					dbl = // Double dbl=
					new Double(EGS4.EKE1[MEDIUM - 1] * EGS4.ELKE
							+ EGS4.EKE0[MEDIUM - 1]);
					LELKE = dbl.intValue();

					// $EVALUATE DEDXMI USING EDEDX(ELKE);
					// $EVALUATE DEDXMJ USING PDEDX(ELKE);
					DEDXMI = EGS4.EDEDX1[LELKE - 1][MEDIUM - 1] * EGS4.ELKE
							+ EGS4.EDEDX0[LELKE - 1][MEDIUM - 1];
					DEDXMJ = EGS4.PDEDX1[LELKE - 1][MEDIUM - 1] * EGS4.ELKE
							+ EGS4.PDEDX0[LELKE - 1][MEDIUM - 1];
					ERATIO[I - 1] = DEDXME / DEDXMI;// ERATIO(I)=DEDXME/DEDXMI;
					PRATIO[I - 1] = DEDXMP / DEDXMJ;// PRATIO(I)=DEDXMP/DEDXMJ;
					// write(6,*) ' Medium ',i,' e-: ',ERATIO(I),' e+:
					// ',PRATIO(I);
				}// ]"end of else"
			}// ]"end loop on media"
			EGS4.EKE = EKEtmp;
			EGS4.ELKE = ELKEtmp;// reseeeeeeeeeeeeeeeeeeeeeeeeeeeeet
			// System.out.println(ERATIO[imoi-1]+ "  ");
		}// ]"end of initialization block"
			// System.out.println(ERATIO[imoi-1]+ "  ");
		if (IARG == 0) {// ["about to transport a particle"
			if (IQL != 0) {// ["count charged particle steps taken"
				if (EGS4SrcEns.NHSTRY == SCSTP_LAST) {// [
					SCSTP_TMP = SCSTP_TMP + 1;
				}// ]
				else {// [
					SCSTP = SCSTP + SCSTP_TMP;
					SCSTP2 = SCSTP2 + SCSTP_TMP * SCSTP_TMP;
					SCSTP_TMP = 1;
					SCSTP_LAST = EGS4SrcEns.NHSTRY;
				}// ]
			}// ]
			else {// ["photon step - play russian roulette?"
				if (RUSROU && (EGS4.W[EGS4.NP - 1] > 0.0)) {// ["yes, play if crosses rrz "
					if ((EGS4.Z[EGS4.NP - 1] <= RRZ)
							&& (EGS4.Z[EGS4.NP - 1] + EGS4.USTEP
									* EGS4.W[EGS4.NP - 1] >= RRZ)) {// "crosses"
																	// $RANDOMSET
																	// XSI;
						XSI = EGS4.random01();
						if (XSI < RRCUT) {// [ "particle survives"
							EGS4.WT[EGS4.NP - 1] = WTL / RRCUT;
						} else {
							// "discard particle on next call to HOWFAR"
							EGS4.WT[EGS4.NP - 1] = 0.0;
						}
					}// ] "end test if crosses russian roulette plane"
				}// ] "end test for playing russian roulette"
			}// ] "end test for photon step"
		}// ] "end test for iarg = 0"

		// "    Set flags for electrons created below DELTA"

		if (IFANO == 1) {// [ "option to regenerate photons after interaction"
							// "and throw away scattered particles"
			if (IARG == 15 || IARG == 17 || IARG == 19 || IARG == 23) {
				// "A pair/Compton/photoelectric/pair event is about to take place"
				EGS4.NP = EGS4.NP + 1;// "Boost the stack"
				if (EGS4.NP + 1 > EGS4.$MXSTACK) {// "check there is space on STACK"
													// OUTPUT; ( ' Fano
													// calculation unable to
													// boost stack.'/
													// ' Stopping.'/
													// 1x,80('*')/);
													// stop;
					EGS4.STOPPROGRAM = true;
					EGS4.seqStr = " ***************************************************"
							+ "  \n"
							+ " Fano calculation unable to boost stack.Stopping."
							+ "  \n"
							+ " ***************************************************";
					// if(EGS4.iprint>0)
					printSequence(EGS4.seqStr);
					return;
				}// ]
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

				// e(np) = e(np - 1);
				// u(np) = u(np - 1);
				// v(np) = v(np - 1);
				// w(np) = w(np - 1);
				// iq(np) = iq(np - 1);
				return;
			}// ]

			// "Throw away any scattered photons from the primary interaction site."
			// " NPold points to the particle before the last discrete interaction."
			// " This change was necessary  for the implementation of atomic"
			// " relaxations. So, check all particles"
			// " between NPold and NP and discard photons "
			if (IARG == 18// " Compton has occured"
					|| IARG == 20 // " After photo-absorption "
					|| IARG == 24 // " After Rayleigh "
					|| IARG == 7 // " After brems "
					|| IARG == 13 // " After annihilation "
					|| IARG == 14)// " After annihilation at rest "
			{
				for (ip = EGS4.NPold; ip <= EGS4.NP; ip++)// DO ip=NPold,NP [
				{
					if (EGS4.IQ[ip - 1] == 0)// if( iq(ip) = 0 )
					{
						// wt(ip) = 0; e(ip) = 0;
						EGS4.WT[ip - 1] = 0.;
						EGS4.E[ip - 1] = 0.;
					}
				}// ]
					// "Note in the case of Bound Compton event not happening, there"
					// "are two identical photons on the stack because of the boosting"
					// "done just above. In this case NP=NPOLD but we want to delete the"
					// "photon at NP, even though it is not scattered!"
			}// ] "end of block for when scattered photon created"
		}// ] " end of ifano = 1 block"

		// "in the following we make use of the knowledge that internal Russian"
		// "Roulette is NOT being used"
		if (IARG == 18 || IARG == 20) {// "after compton or pe event"
										// "flag e- created below DELTA"
			for (ip = EGS4.NPold; ip <= EGS4.NP; ip++)// DO ip=NPold, NP [
			{
				if (EGS4.IQ[ip - 1] == -1 && EGS4.E[EGS4.NP - 1] < DELTAT) {
					EGS4.LATCH[ip - 1] = 99;
				}
			}// ]
		}// ]"end of iarg = 18 or 20 block"
		else if (IARG == 16) {// " after a pair event, lower energy particle is NP"
			if (EGS4.E[EGS4.NP - 1] < DELTAT) {
				EGS4.LATCH[EGS4.NP - 1] = 99;
			}
			// "at most one particle can be below DELTAT, except for DELTA > 500 keV"
		}// ]
			// "end of flag setting"

		if (IWATCH > 0)
			EGS4.WATCH(IARG, IWATCH);// CALL WATCH(IARG,IWATCH);
										// "signal watch routine if active"

		// "Score the energy as required for the different modes"
		// "****************************************************"

		IFLAG = 0; // "left zero unless discards a brem below AP, fluorescent discard"
		// "         or compton, PE or pair  created e- below ECUT -LATCH=99"

		if ((EGS4.EDEP != 0.0) && (WTL > 0.0) && (IARG < 5)) {// "energy deposited"

			// "score energy deposited"
			FTMP = WTL * EGS4.EDEP;

			// "Note, EKE is not set after a compton event"
			// "if the scattered electron is about to be thrown out"
			// "so do not use internal EKE at this point, calculate it here"
			// "We can also get to this point when an electron has fallen below"
			// "AE, and in fact, due to numerical round-off it might get"
			// "below PRM. We don't need EKE for these track ends but except for"
			// "a test - really we just don't need the log"
			if (IQL != 0) {// [
				EGS4.EKE = EGS4.E[EGS4.NP - 1] - EGS4.PRM;
				if (EGS4.EKE <= 0.0) {
					EGS4.EKE = 1.E-5;
					EGS4.ELKE = 1E30;
				}// "big so we notice if use it!" ]
				else {
					EGS4.ELKE = Math.log(EGS4.EKE);
				}
			}// ]

			// "store k.e. for start of step FOR UNIFORM USE BELOW"
			EKEtmp = EGS4.EKE;
			ELKEtmp = EGS4.ELKE;

			if (IARG == 0) {// "step to occur so this is not a track end"
							// "and is not a photon since edep ~=0      "
							// "so adjust k.e. to mid-step              "
				EGS4.EKE = EKEtmp - 0.5 * EGS4.EDEP;
				EGS4.ELKE = Math.log(EGS4.EKE);
			}// ]
			else if (IQL != 0 && EGS4.i_do_rr[IRL - 1] == 0)// EGS4.i_do_rr(IRL)==0)
			{
				// "this must be a charged particle track end because IARG ~=0"
				// "it can also be  an IARG=4 deposit of energy below binding"
				// "if relaxation is on - and in this case, the energy is that of"
				// "the electron we are following - withou relaxation it must"
				// "have always been the photon on the top of stack?"
				if (EGS4.EKE >= DELTA && (IARG != 4)) {// [
														// "not  possible - this proves it"
														// "note that DELTA is a kinetic energy "
														// "IARG=4 is discard of energy below binding"
														// OUTPUT EKE;(' AUSGAB:
														// EKE=',F15.9,' is
														// impossible***'/);
														// "recall ECUT = AE = DELTAT, and we have not allowed"
														// "range rejection. The only other discard is out of"
														// "geometry and that has no edep so doesn't get to here"
														// OUTPUT IARG, IQL,
														// DELTA,E(NP),AE(1),PRM,RM,EDEP;
														// (' IARG, IQL, DELTA,
														// E(NP), AE(1), PRM,
														// RM, EDEP:'/
														// 2I5, 6F11.7);
														// STOP;
					EGS4.STOPPROGRAM = true;
					EGS4.seqStr = " ***************************************************"
							+ "  \n"
							+ " Impossible kinetic energy EKE in charged Roussian Roulette range rejection.Stopping."
							+ "  \n"
							+ " ***************************************************";
					// if(EGS4.iprint>0)
					printSequence(EGS4.seqStr);
					return;

				}// ]
			}// ]"end of charged particle track end block";
				// "Why do we have the above anymore? We know we never get these messages"
				// "to catch logic problems as with the relaxation!"

			if (IQL != 0 && IARG != 4) {// "e- or e+, either taking step or being discarded"
										// "but we must explicitly exclude IARG=4 discards"
										// "which might happen with an electron on the top of stack"
				if (EGS4.LATCH[EGS4.NP - 1] == 99) {// "this must be a discarded electron which "
													// "was created below delta"
					IFLAG = 1;// "we do not include these in spr"

					// IF(IWATCH > 0 & IWATCH ~=4 )[OUTPUT IARG,IRL,EKEtmp;
					// (' AUSGAB: discard e- or e+ created below ECUT:
					// IARG,IRL,'
					// ,'EKEtmp:',2I5,F15.8);
					// ]
					EGS4.LATCH[EGS4.NP - 1] = 0; // "just in case"
				}// ]"end of block for discarded e- or e+ created below delta"

				else if (EKEtmp - EGS4.EDEP > DELTA) {// "a straight transport step "
														// "the energy is always > DELTA"
					MEDIUM = imoi;// 2; //"set to second medium to get "
					// $SET INTERVAL ELKE,EKE;
					// "prepare to interpolate in  data tables"
					Double dbl = new Double(EGS4.EKE1[MEDIUM - 1] * EGS4.ELKE
							+ EGS4.EKE0[MEDIUM - 1]);
					LELKE = dbl.intValue();

					if (IQL == -1) {
						// $EVALUATE DEDXM2 USING EDEDX(ELKE);
						DEDXM2 = EGS4.EDEDX1[LELKE - 1][MEDIUM - 1] * EGS4.ELKE
								+ EGS4.EDEDX0[LELKE - 1][MEDIUM - 1];
					} else {
						// $EVALUATE DEDXM2 USING PDEDX(ELKE);
						DEDXM2 = EGS4.PDEDX1[LELKE - 1][MEDIUM - 1] * EGS4.ELKE
								+ EGS4.PDEDX0[LELKE - 1][MEDIUM - 1];
					}
					MEDIUM = MEDL;// "set back to current medium"
					// $SET INTERVAL ELKE,EKE;
					// "prepare to interpolate in  data tables"
					dbl = new Double(EGS4.EKE1[MEDIUM - 1] * EGS4.ELKE
							+ EGS4.EKE0[MEDIUM - 1]);
					LELKE = dbl.intValue();

					if (IQL == -1) {
						// $EVALUATE DEDXM1 USING EDEDX(ELKE);
						DEDXM1 = EGS4.EDEDX1[LELKE - 1][MEDIUM - 1] * EGS4.ELKE
								+ EGS4.EDEDX0[LELKE - 1][MEDIUM - 1];
					} else {
						// $EVALUATE DEDXM1 USING PDEDX(ELKE);
						DEDXM1 = EGS4.PDEDX1[LELKE - 1][MEDIUM - 1] * EGS4.ELKE
								+ EGS4.PDEDX0[LELKE - 1][MEDIUM - 1];
					}
					SRATIO = DEDXM2 / DEDXM1; // "ratio of energy depositions"
				}// ]"end of straight transport step"
				else {// "These are not electrons created below DELTA"
						// "and they are not electrons with energy always above DELTA"
						// "this is mostly for track ends, i.e. particles with "
						// "IARG = 1,2,3,4"
						// "but also has e- or e+ crossing AE in their IARG=0 step"

					// "this is an approximation but can we do better? Yes, we can"
					// "i.e.  we could be on a straight transport which crosses AE"
					// "and has EKE < DELTA.  We treat as if at AE"

					if (IARG != 0) {// "stoppers  IARG = 1,2,3,4"
						if (IQL == -1) {
							SRATIO = ERATIO[MEDL - 1];
						} else {
							SRATIO = PRATIO[MEDL - 1];
						}// PRATIO(MEDL);]
						// "SRATIO is just ratio at DELTA"
					}// ]
					else {// "This is a step crossing delta, we set therefore"
							// "sratio to ( (e-delta)*Sratio((e+delta)/2)+     "
							// "            (edep-e+delta)*Sratio(delta))/edep "
							// "i.e. we take energy above delta and use ratio at"
							// "mid-point of that part of step and energy below"
							// "DELTA as being the ration at DELTA"
						EGS4.EKE = 0.5 * (EKEtmp + DELTA);
						EGS4.ELKE = Math.log(EGS4.EKE); // "mid-point above delta"
						MEDIUM = imoi;// 2; "set to second medium"
						// $SET INTERVAL ELKE,EKE;
						// "prepare to interpolate in  data tables"
						Double dbl = new Double(EGS4.EKE1[MEDIUM - 1]
								* EGS4.ELKE + EGS4.EKE0[MEDIUM - 1]);
						LELKE = dbl.intValue();

						if (IQL == -1) {
							// $EVALUATE DEDXM2 USING EDEDX(ELKE);
							DEDXM2 = EGS4.EDEDX1[LELKE - 1][MEDIUM - 1]
									* EGS4.ELKE
									+ EGS4.EDEDX0[LELKE - 1][MEDIUM - 1];
						} else {
							// $EVALUATE DEDXM2 USING PDEDX(ELKE);
							DEDXM2 = EGS4.PDEDX1[LELKE - 1][MEDIUM - 1]
									* EGS4.ELKE
									+ EGS4.PDEDX0[LELKE - 1][MEDIUM - 1];
						}
						MEDIUM = MEDL; // "set to current medium"
						// $SET INTERVAL ELKE,EKE;
						// "prepare to interpolate in  data tables"
						dbl = new Double(EGS4.EKE1[MEDIUM - 1] * EGS4.ELKE
								+ EGS4.EKE0[MEDIUM - 1]);
						LELKE = dbl.intValue();

						if (IQL == -1) {
							// $EVALUATE DEDXM1 USING EDEDX(ELKE);
							DEDXM1 = EGS4.EDEDX1[LELKE - 1][MEDIUM - 1]
									* EGS4.ELKE
									+ EGS4.EDEDX0[LELKE - 1][MEDIUM - 1];
							SRATIO = DEDXM2 / DEDXM1 * (EKEtmp - DELTA)
									+ ERATIO[MEDL - 1]
									* (EGS4.EDEP - EKEtmp + DELTA);
						}// ]
						else {// [
								// $EVALUATE DEDXM1 USING PDEDX(ELKE);
							DEDXM1 = EGS4.PDEDX1[LELKE - 1][MEDIUM - 1]
									* EGS4.ELKE
									+ EGS4.PDEDX0[LELKE - 1][MEDIUM - 1];
							SRATIO = DEDXM2 / DEDXM1 * (EKEtmp - DELTA)
									+ PRATIO[MEDL - 1]
									* (EGS4.EDEP - EKEtmp + DELTA);
						}// ]
						SRATIO = SRATIO / EGS4.EDEP;
						ICROSS = ICROSS + 1; // "count the crossers"
					}// ] "end of crossers block"

					// "check how often things are crossing rather than track ends"
					// "results passed back in SCORE"
					ITRACKE = ITRACKE + 1;
				}// ]"end of track-end block -excludes if created with energy below delta"
			}// ]"end of e- and e+ block"

			else if (IQL == 0 || IARG == 4) {// "this is a photon created below PCUT(BREM)"
												// " or is a fluorescent discard"
				IFLAG = 1;// "we do not include this at present in the spr"
			}// ]

			if (IARG == 0) {
				// "recover original energy of particle"
				EGS4.EKE = EKEtmp;
				EGS4.ELKE = ELKEtmp;
			}// ]

			if (EGS4SrcEns.NHSTRY == SCDOSE_LAST[IZ - 1][IX - 1])// (IZ,IX))
			{// [
				SCDOSE_TMP3[IZ - 1][IX - 1] = SCDOSE_TMP3[IZ - 1][IX - 1]
						+ FTMP;
				if (IFLAG == 0) {// "do not include bremsstrahlung photons created below PCUT"
									// "or discarded fluorescent energy for IT=1 and 2"
									// SCEDPN_TMP(IZ,IX,1)=SCEDPN_TMP(IZ,IX,1)+FTMP;
					SCEDPN_TMP[IZ - 1][IX - 1][0] = SCEDPN_TMP[IZ - 1][IX - 1][0]
							+ FTMP;
					SCEDPN_TMP[IZ - 1][IX - 1][1] = SCEDPN_TMP[IZ - 1][IX - 1][1]
							+ FTMP * SRATIO;
				}// ]
				if (IARG == 0) {// "SCORE SANS ANY TRACK ENDS -recall we are in the EDEP block"
								// SCEDPN_TMP(IZ,IX,3)=SCEDPN_TMP(IZ,IX,3)+FTMP;
					SCEDPN_TMP[IZ - 1][IX - 1][2] = SCEDPN_TMP[IZ - 1][IX - 1][2]
							+ FTMP;
					SCEDPN_TMP[IZ - 1][IX - 1][3] = SCEDPN_TMP[IZ - 1][IX - 1][3]
							+ FTMP * SRATIO;
				}// ]
			}// ]
			else {
				// "add to IT=1,2,3,4 regardless of whether there is a contribution or"
				// "not"
				SCEDPN[IZ - 1][IX - 1][0] = SCEDPN[IZ - 1][IX - 1][0]
						+ SCEDPN_TMP[IZ - 1][IX - 1][0];
				SCEDPN2[IZ - 1][IX - 1][0] = SCEDPN2[IZ - 1][IX - 1][0]
						+ SCEDPN_TMP[IZ - 1][IX - 1][0]
						* SCEDPN_TMP[IZ - 1][IX - 1][0];

				SCEDPN[IZ - 1][IX - 1][1] = SCEDPN[IZ - 1][IX - 1][1]
						+ SCEDPN_TMP[IZ - 1][IX - 1][1];
				SCEDPN2[IZ - 1][IX - 1][1] = SCEDPN2[IZ - 1][IX - 1][1]
						+ SCEDPN_TMP[IZ - 1][IX - 1][1]
						* SCEDPN_TMP[IZ - 1][IX - 1][1];

				SCSPR_COV[IZ - 1][IX - 1][0] = SCSPR_COV[IZ - 1][IX - 1][0]
						+ SCEDPN_TMP[IZ - 1][IX - 1][1]
						* SCEDPN_TMP[IZ - 1][IX - 1][0];

				SCEDPN[IZ - 1][IX - 1][2] = SCEDPN[IZ - 1][IX - 1][2]
						+ SCEDPN_TMP[IZ - 1][IX - 1][2];
				SCEDPN2[IZ - 1][IX - 1][2] = SCEDPN2[IZ - 1][IX - 1][2]
						+ SCEDPN_TMP[IZ - 1][IX - 1][2]
						* SCEDPN_TMP[IZ - 1][IX - 1][2];

				SCEDPN[IZ - 1][IX - 1][3] = SCEDPN[IZ - 1][IX - 1][3]
						+ SCEDPN_TMP[IZ - 1][IX - 1][3];
				SCEDPN2[IZ - 1][IX - 1][3] = SCEDPN2[IZ - 1][IX - 1][3]
						+ SCEDPN_TMP[IZ - 1][IX - 1][3]
						* SCEDPN_TMP[IZ - 1][IX - 1][3];

				SCSPR_COV[IZ - 1][IX - 1][1] = SCSPR_COV[IZ - 1][IX - 1][1]
						+ SCEDPN_TMP[IZ - 1][IX - 1][3]
						* SCEDPN_TMP[IZ - 1][IX - 1][2];

				// "now add contributions to doses"
				SCDOSE[IZ - 1][IX - 1][0] = SCDOSE[IZ - 1][IX - 1][0]
						+ SCEDPN_TMP[IZ - 1][IX - 1][0]
						- SCEDPN_TMP[IZ - 1][IX - 1][2];
				SCDOSE2[IZ - 1][IX - 1][0] = SCDOSE2[IZ - 1][IX - 1][0]
						+ (SCEDPN_TMP[IZ - 1][IX - 1][0] - SCEDPN_TMP[IZ - 1][IX - 1][2])
						* (SCEDPN_TMP[IZ - 1][IX - 1][0] - SCEDPN_TMP[IZ - 1][IX - 1][2]);

				SCDOSE_COV[IZ - 1][IX - 1][0] = SCDOSE_COV[IZ - 1][IX - 1][0]
						+ SCDOSE_TMP3[IZ - 1][IX - 1]
						* (SCEDPN_TMP[IZ - 1][IX - 1][0] - SCEDPN_TMP[IZ - 1][IX - 1][2]);

				SCDOSE[IZ - 1][IX - 1][1] = SCDOSE[IZ - 1][IX - 1][1]
						+ SCDOSE_TMP3[IZ - 1][IX - 1]
						- SCEDPN_TMP[IZ - 1][IX - 1][0];
				SCDOSE2[IZ - 1][IX - 1][1] = SCDOSE2[IZ - 1][IX - 1][1]
						+ (SCDOSE_TMP3[IZ - 1][IX - 1] - SCEDPN_TMP[IZ - 1][IX - 1][0])
						* (SCDOSE_TMP3[IZ - 1][IX - 1] - SCEDPN_TMP[IZ - 1][IX - 1][0]);

				SCDOSE_COV[IZ - 1][IX - 1][1] = SCDOSE_COV[IZ - 1][IX - 1][1]
						+ SCDOSE_TMP3[IZ - 1][IX - 1]
						* (SCDOSE_TMP3[IZ - 1][IX - 1] - SCEDPN_TMP[IZ - 1][IX - 1][0]);

				SCDOSE[IZ - 1][IX - 1][2] = SCDOSE[IZ - 1][IX - 1][2]
						+ SCDOSE_TMP3[IZ - 1][IX - 1];
				SCDOSE2[IZ - 1][IX - 1][2] = SCDOSE2[IZ - 1][IX - 1][2]
						+ SCDOSE_TMP3[IZ - 1][IX - 1]
						* SCDOSE_TMP3[IZ - 1][IX - 1];

				if (IFLAG == 0) {
					SCEDPN_TMP[IZ - 1][IX - 1][0] = FTMP;
					SCEDPN_TMP[IZ - 1][IX - 1][1] = FTMP * SRATIO;
				}// ]
				else {// [
					SCEDPN_TMP[IZ - 1][IX - 1][0] = 0.;
					SCEDPN_TMP[IZ - 1][IX - 1][1] = 0.;
				}// ]
				if (IARG == 0) {
					SCEDPN_TMP[IZ - 1][IX - 1][2] = FTMP;
					SCEDPN_TMP[IZ - 1][IX - 1][3] = FTMP * SRATIO;
				}// ]
				else {// [
					SCEDPN_TMP[IZ - 1][IX - 1][2] = 0.;
					SCEDPN_TMP[IZ - 1][IX - 1][3] = 0.;
				}// ]
				SCDOSE_TMP3[IZ - 1][IX - 1] = FTMP;
				SCDOSE_LAST[IZ - 1][IX - 1] = EGS4SrcEns.NHSTRY;
				SCNUM[IZ - 1][IX - 1] = SCNUM[IZ - 1][IX - 1] + 1;
			}// ]

			// IF((IWATCH > 1) & (IWATCH ~= 4))[
			// OUTPUT FTMP/RHOR(IRL),IRL,IARG;
			// (5x,' ***EDEP scored = ',E14.7,' MeV. IRL= ',I3,' IARG= ',I3);
			// OUTPUT FTMP*SRATIO/RHO(2),IRL,IARG;
			// (5x,' ***EDEP*SRATIO/RHO(2) scored = ',
			// E14.7,' MeV. IRL= ',I3,' IARG= ',I3);
			// ]
		}// ] "end of energy deposited block"

		// "score the total curved pathlength in the volume"

		return;
	}// "END OF AUSGAB"

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
		ix = (ir - 2) / EGS4Geom.NZ + 1;// NZ=planar zones=NPLANE-1!!!NR=nr of
										// Radii
		iz = ir - 1 - EGS4Geom.NZ * (ix - 1);
		r = Math.sqrt(x * x + y * y);
		EGS4.tperp = EGS4.min(z - EGS4Geom.ZPLANE[iz - 1], EGS4Geom.ZPLANE[iz]
				- z, EGS4Geom.RCYL[ix] - r);
		if (ix != 1) {
			EGS4.tperp = Math.min(EGS4.tperp, r - EGS4Geom.RCYL[ix - 1]);
		}

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

		// "DISCARD ZERO WEIGHT PARTICLES"
		if (EGS4.WT[EGS4.NP - 1] == 0.0) {
			EGS4.IDISC = 1;
			return;
		}

		// "INITIALLY ASSUME PARTICLE STAYS IN THE TARGET"
		boolean OUTEND = false;
		boolean OUTSID = false;
		//int IQL = EGS4.IQ[EGS4.NP - 1];
		int IRL = EGS4.IR[EGS4.NP - 1];// "LOCAL REGION NUMBER"

		// "DISCARD IF PARTICLE WANTS TO LEAVE THE GEOMETRY OR OF THE REGION IS TOTALLY"
		// "ABSORBING"
		if ((IRL == 1))// ||(EGS4Grid.CABSRB[IRL-1].compareTo(ACHAR)==0))
		{
			EGS4.IDISC = 1;
			return;
		}

		// $GET-IX-IZ(IRL); //"GET PLANAR AND CYLINDRICAL ZONES NUMBERS"
		int IX = EGS4Geom.GET_IX(IRL);
		int IZ = EGS4Geom.GET_IZC(IRL);

		EGS4Geom.PLANES2(IZ, IZ + 1);
		EGS4Geom.CYLNDR2(IX);

		/*
		 * REPLACE {$SET NEW REGION(#,#);} WITH
		 * {;IF(#1.LE.USTEP)[USTEP=#1;IRNEW=#2;RETURN;]}
		 */
		// "CALCULATE DNEAR"
		double ZL = 0.0;
		double RL = 0.0;
		double ZLEFT = 0.0;
		double ZRIGHT = 0.0;
		double ROUT = 0.0;
		double DNTMP = 0.0;
		double RIN = 0.0;
		if ((EGS4.IQ[EGS4.NP - 1] != 0)) {
			ZL = EGS4.Z[EGS4.NP - 1];
			RL = Math.sqrt(EGS4.X[EGS4.NP - 1] * EGS4.X[EGS4.NP - 1]
					+ EGS4.Y[EGS4.NP - 1] * EGS4.Y[EGS4.NP - 1]);
			ZLEFT = ZL - EGS4Geom.ZPLANE[IZ - 1];
			ZRIGHT = EGS4Geom.ZPLANE[IZ] - ZL;// +1)-ZL;
			ROUT = EGS4Geom.RCYL[IX] - RL;
			DNTMP = EGS4.min(ZLEFT, ZRIGHT, ROUT);
			if (IX != 1) {
				RIN = RL - EGS4Geom.RCYL[IX - 1];
				DNTMP = Math.min(DNTMP, RIN);
			}
			EGS4.DNEAR[EGS4.NP - 1] = DNTMP;
		}

		int IZNEW = IZ + EGS4Geom.ihitp;// IHITP; "GET NEW PLANAR REGION"
		if ((IZNEW < 1) || (IZNEW > EGS4Geom.NZ))
			OUTEND = true; // "FLAG IF LEAVES BY THE ENDS"

		int IXNEW = IX + EGS4Geom.ihitc;// IHITC; "GET NEW CYLINDRICAL REGION"
		if (IXNEW > EGS4Geom.NR)
			OUTSID = true; // "FLAG IF LEAVES BY THE SIDES"

		// int NWNRCL=0;int IMSOFF=0;
		// "DO MOST PROBABLE CASE FIRST WHERE A PLANE AND A CYLINDER CAN BE HIT"
		if ((EGS4Geom.ihitp != 0) && (EGS4Geom.ihitc != 0)) {
			if (EGS4Geom.tplane < EGS4Geom.tcyl) {// "HITS PLANE FIRST"
				if (OUTEND) {
					// NWNRCL=0;
					// $SET NEW REGION(TPLANE,1);
					// if({P1}.LE.USTEP)[
					if (EGS4Geom.tplane <= EGS4.USTEP) {
						EGS4.USTEP = EGS4Geom.tplane;// {P1};
						EGS4.IRNEW = 1;// {P2};
						return;// IMSOFF=0;
						// EGS4Macro.NEWNRC=NWNRCL;
						// if(EGS4.USTEP==0 && EGS4.IRNEW!=EGS4.IROLD)
						// DECISION=1;
					}
					// return;
				} else {
					// NWNRCL=20-10*EGS4Geom.ihitp;//IHITP;
					// $SET NEW REGION(TPLANE,IRL+IHITP);
					// if({P1}.LE.USTEP)[
					if (EGS4Geom.tplane <= EGS4.USTEP) {
						EGS4.USTEP = EGS4Geom.tplane;// {P1};
						EGS4.IRNEW = IRL + EGS4Geom.ihitp;// {P2};
						return;// IMSOFF=0;
						// EGS4Macro.NEWNRC=NWNRCL;
						// if(EGS4.USTEP==0 && EGS4.IRNEW!=EGS4.IROLD)
						// DECISION=1;
					}
					// return;
				}
			} else if (EGS4Geom.tcyl < EGS4Geom.tplane) {// "HITS CYLINDER FIRST"
				if (OUTSID) {
					// NWNRCL=0;
					// $SET NEW REGION(TCYL,1);
					// if({P1}.LE.USTEP)[
					if (EGS4Geom.tcyl <= EGS4.USTEP) {
						EGS4.USTEP = EGS4Geom.tcyl;// {P1};
						EGS4.IRNEW = 1;// {P2};
						return;// IMSOFF=0;
						// EGS4Macro.NEWNRC=NWNRCL;
						// if(EGS4.USTEP==0 && EGS4.IRNEW!=EGS4.IROLD)
						// DECISION=1;
					}
					// return;
				} else {
					// NWNRCL=30+10*EGS4Geom.ihitc;
					// $SET NEW REGION(TCYL,IRL+NZ*IHITC);
					// if({P1}.LE.USTEP)[
					if (EGS4Geom.tcyl <= EGS4.USTEP) {
						EGS4.USTEP = EGS4Geom.tcyl;// {P1};
						EGS4.IRNEW = IRL + EGS4Geom.NZ * EGS4Geom.ihitc;// {P2};
						return;// IMSOFF=0;
						// EGS4Macro.NEWNRC=NWNRCL;
						// if(EGS4.USTEP==0 && EGS4.IRNEW!=EGS4.IROLD)
						// DECISION=1;
					}
					// return;
				}
			} else {// "ODD CASE TCYL=TPLANE:HITS PLANE AND CYLINDER TOGETHER"
				if (OUTEND || OUTSID) {
					// NWNRCL=0;
					// $SET NEW REGION(TCYL,1);
					// if({P1}.LE.USTEP)[
					if (EGS4Geom.tcyl <= EGS4.USTEP) {
						EGS4.USTEP = EGS4Geom.tcyl;// {P1};
						EGS4.IRNEW = 1;// {P2};
						return;// IMSOFF=0;
						// EGS4Macro.NEWNRC=NWNRCL;
						// if(EGS4.USTEP==0 && EGS4.IRNEW!=EGS4.IROLD)
						// DECISION=1;
					}
					// return;
				} else {
					// NWNRCL=20-10*EGS4Geom.ihitp;//IHITP;
					// $SET NEW REGION(TCYL,IRL+IHITP+NZ*IHITC);
					// if({P1}.LE.USTEP)[
					if (EGS4Geom.tcyl <= EGS4.USTEP) {
						EGS4.USTEP = EGS4Geom.tcyl;// {P1};
						EGS4.IRNEW = IRL + EGS4Geom.ihitp + EGS4Geom.NZ
								* EGS4Geom.ihitc;// {P2};
						return;// IMSOFF=0;
						// EGS4Macro.NEWNRC=NWNRCL;
						// if(EGS4.USTEP==0 && EGS4.IRNEW!=EGS4.IROLD)
						// DECISION=1;
					}
					// return;
				}
			}
		}

		// "DO ODD CASE-PARTICLE CAN HIT PLANE BUT NOT CYLINDER"
		else if (EGS4Geom.ihitp != 0) {
			if (OUTEND) {
				// NWNRCL=0;
				// $SET NEW REGION(TPLANE,1);
				// if({P1}.LE.USTEP)[
				if (EGS4Geom.tplane <= EGS4.USTEP) {
					EGS4.USTEP = EGS4Geom.tplane;// {P1};
					EGS4.IRNEW = 1;// {P2};
					return;// IMSOFF=0;
					// EGS4Macro.NEWNRC=NWNRCL;
					// if(EGS4.USTEP==0 && EGS4.IRNEW!=EGS4.IROLD) DECISION=1;
				}
				// return;
			} else {
				// NWNRCL=20-10*EGS4Geom.ihitp;
				// $SET NEW REGION(TPLANE,IRL+IHITP);
				// if({P1}.LE.USTEP)[
				if (EGS4Geom.tplane <= EGS4.USTEP) {
					EGS4.USTEP = EGS4Geom.tplane;// {P1};
					EGS4.IRNEW = IRL + EGS4Geom.ihitp;// {P2};
					return;// IMSOFF=0;
					// EGS4Macro.NEWNRC=NWNRCL;
					// if(EGS4.USTEP==0 && EGS4.IRNEW!=EGS4.IROLD) DECISION=1;
				}
				// return;
			}
		}

		// "DO ODD CASE-PARTICLE CAN HIT CYLINDER BUT NOT PLANE"
		else {
			if (OUTSID) {
				// NWNRCL=0;
				// $SET NEW REGION(TCYL,1);
				// if({P1}.LE.USTEP)[
				if (EGS4Geom.tcyl <= EGS4.USTEP) {
					EGS4.USTEP = EGS4Geom.tcyl;// {P1};
					EGS4.IRNEW = 1;// {P2};
					return;// IMSOFF=0;
					// EGS4Macro.NEWNRC=NWNRCL;
					// if(EGS4.USTEP==0 && EGS4.IRNEW!=EGS4.IROLD) DECISION=1;
				}
				// return;

			} else {
				// NWNRCL=30+10*EGS4Geom.ihitc;//IHITC;
				// $SET NEW REGION(TCYL,IRL+NZ*IHITC);
				// if({P1}.LE.USTEP)[
				if (EGS4Geom.tcyl <= EGS4.USTEP) {
					EGS4.USTEP = EGS4Geom.tcyl;// {P1};
					EGS4.IRNEW = IRL + EGS4Geom.NZ * EGS4Geom.ihitc;// {P2};
					return;// IMSOFF=0;
					// EGS4Macro.NEWNRC=NWNRCL;
					// if(EGS4.USTEP==0 && EGS4.IRNEW!=EGS4.IROLD) DECISION=1;
				}
				// return;
			}
		}

		// "AT THIS STAGE ALL GEOMETRICAL POSSIBILITIES HAVE BEEN CHECKED AND CONTROL"
		// "HAS ALREADY BEEN TRANSFERRED TO EGS"

		// return;
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
		TITLEs = "dosrznrc_template--depth dose in H2O due to Cobalt beam";
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
		// IOOPTN=IOOPTN_MATERIAL_SUMMARY;//IOOPTN_SHORT;
		// @ ELECTRON TRANSPORT
		// " = normal (0) normal electron transport (discrete interactions)
		// " = no interactions (1) no discrete interactions (used for CDSA
		// " calculations but note that special data
		// " sets are also needed to do a full CSDA
		// " calculation. All turning off
		// " interactions does is just that. See use
		// " of IUNRST=2,3,4 PEGS4 data sets for real CSDA)
		// " [ ICSDA]
		// EGS4Macro.ICSDA=ETRANS_NORMAL;
		// ==============================
		// " DOSE ZBOUND MIN (I) Minimum plane # defining dose region
		// (default=1)
		// " [NZDMIN]
		// " DOSE ZBOUND MAX (I) Maximum plane # defining dose region
		// " [NZDMAX]
		// " DOSE RBOUND MIN (I) Minimum cylinder # defining dose region
		// (default=0)
		// " [NRDMIN]
		// " DOSE RBOUND MAX (I) Maximum cylinder # defining dose region
		// " [NRDMAX]
		// ============================================
		// IPRNTP=IPRNTP_ALL;
		// IPRIMARY=IPRIMARY_TOTAL_FLUENCE
		// NUM_STARTR=2;
		// NUM_STOPR=EGS4.$MXREG;
		if (IDAT == 1) {
			INEXT = 0;
		} else {
			INEXT = 1;
		}
		if (IRESTART == 4) {
			IDAT = 1; // "do not store output in this case to avoid biasing"
			ISTORE = 0; // "do not store the starting random numbers either"
		}
		// "Region numbers for printing in listing files is done below GEOMRZ."
		ISPRREG = 0;// always by region
		NUM_STARTR = 2;
		NUM_STOPR = EGS4.$MXREG;
		// =================================
		// NZDMIN=1;
		// NZDMAX=61;
		// NRDMIN=0;
		// NRDMAX=60;
		// ===========================
		// @ STORE DATA ARRAYS
		// #yes,no;
		// #yes: output .egsdat file for restarts, parallel post-processing, etc
		// IDAT=IDAT_NO;
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
		// EGS4Macro.IFULL=IFULL_PULSE_HEIGHT_DISTRIBUTION;//IFULL_AATT_AND_ASCAT;
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
		// @SCORE KERMA
		// IKERMA=KERMA_YES;
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
				// EGS4Geom.ZPLANE[1]=0.05;
				// EGS4Geom.ZPLANE[2]=6.35;
				// EGS4Geom.ZPLANE[3]=7.30;
			}

			// EGS4Geom.nCyl=2;//"number of radial cylinders input"
			// #Radii of cylinders
			// EGS4Geom.RCYL[1]=3.15;//6.3/2
			// EGS4Geom.RCYL[2]=3.65;//7.3/2

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

			// EGS4Geom.nMEDNUM=2;
			// EGS4Geom.MEDNUM[0]=1;//Al
			// EGS4Geom.MEDNUM[1]=2;//NAI

			if (EGS4Geom.DESCRIBE == EGS4Geom.DESCRIBE_REGIONS_DENSITY
					|| EGS4Geom.DESCRIBE == EGS4Geom.DESCRIBE_PLANES_DENSITY) {
				EGS4Geom.nRHOR = 1;
				EGS4Geom.RHOR[0] = 0.;
			}
			if (EGS4Geom.DESCRIBE == EGS4Geom.DESCRIBE_REGIONS
					|| EGS4Geom.DESCRIBE == EGS4Geom.DESCRIBE_REGIONS_DENSITY) {
				// EGS4Geom.nNREGLO=2;
				// EGS4Geom.nNREGHI=2;
				// EGS4Geom.NREGLO[0]=2;//START REGION
				// EGS4Geom.NREGHI[0]=7;//STOP REGION
				// EGS4Geom.NREGLO[1]=3;//START REGION
				// EGS4Geom.NREGHI[1]=3;//STOP REGION
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
			destroyArrays();
			return;
		}
		// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
		// "This is done here because NREG is needed from geomrz..."
		// if(IPRNTP==2)
		// {
		// OUTPUT;(/'***** SPECTRUM OUTPUT IN LISTING FILE *****');
		// for(int IRL=2;IRL<=EGS4Geom.nreg;IRL++) {IPRTSP[IRL-1]=1;}
		// //"DEFAULT PRINT NOTHING"//EGS4Geom.nreg
		for (int J = 1; J <= EGS4Geom.nreg; J++) {
			ISPR[J - 1] = 0;
		}
		// for(int J=1;J<= NUM_STARTR;J++)
		for (int J = 1; J <= 1; J++)// only one region allowed here!!!
		{
			ILOW = NUM_STARTR;// VALUE(NUM_STARTR,J);
			IHIGH = NUM_STOPR;// VALUE(NUM_STOPR,J);
			if (ILOW > EGS4Geom.nreg) {
				ILOW = EGS4Geom.nreg;
			}
			if (IHIGH > EGS4Geom.nreg) {
				IHIGH = EGS4Geom.nreg;
			}
			if (ILOW > IHIGH) {
				IHIGH = ILOW;
			}
			if (ILOW == IHIGH) {
				// EGS4.seqStr="  Print spectrum in region "+ILOW;
				// if(EGS4.iprint>1)
				// printSequence(EGS4.seqStr);

				// OUTPUT ILOW;(' Print spectrum in region',I3);
				ISPR[ILOW - 1] = 1;// 0;
			} else {
				for (int IRL = ILOW; IRL <= IHIGH; IRL++) {
					ISPR[IRL - 1] = 1;
				}
				// OUTPUT ILOW,IHIGH;(' Print spectra in regions',I4,' TO',I4);
				// EGS4.seqStr=" Print spectra in regions "+ILOW+" TO "+IHIGH;
				// if(EGS4.iprint>1)
				// printSequence(EGS4.seqStr);

			}
		}
		// }

		// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
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
			if (scattCaseB) {
				EGS4.seqStr = " Radius of hitting surface-colimatted beam- (cm) = "
						+ EGS4.format(EGS4Geom.RCYL[1], 8, true);
			}
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
		// #######################################PULSE HEIHHT
		// @REGION OF SENSITIVE VOLUME
		// nREGSOLV=1;
		// REGSVOL[0]=3;
		// SLOTE=0.01;
		// DELTAE=0.005;

		if (SLOTE < 0) {
			nTOPEBIN = 2;
			TOPEBIN[0] = 0.12;
			TOPEBIN[1] = 0.02;
		}
		// ###########################################################################################
		// EGS4SrcEns.ipart=EGS4SrcEns.ipart_photon;
		// SOURCE NUMBER= 1 #0,1,2,3,4,10,11,12,13,14,15,16,20,21,22,23
		// EGS4SrcEns.ISOURC=EGS4SrcEns.point_source_on_axis_incident_from_the_front;//1;
		// SOURCE OPTIONS= 100., 1.3, 0, 0
		// #for source 1: SSD of beam, radius of # beam on front surface
		if (EGS4SrcEns.ISOURC != EGS4SrcEns.parallel_beam_incident_from_the_front_with_radial_distribution)// 20)
		{
			// EGS4SrcEns.source_option[0]=3.5;
			// EGS4SrcEns.source_option[1]=7.3;
			// EGS4SrcEns.source_option[2]=0.;
			// EGS4SrcEns.source_option[3]=1.;

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
			destroyArrays();
			return;
		}
		// source initialization
		if (SOURCE == -1)
			EGS4SrcEns.SRCINI();
		if (EGS4.STOPPROGRAM) {
			closeFile();
			destroyArrays();
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
			// EGS4SrcEns.ikemev=0.662;

			EGS4SrcEns.ENSRC();
			if (EGS4.STOPPROGRAM) {
				closeFile();
				destroyArrays();
				return;
			}
		}// "Get data re-source energies"
			// ###########################################################################################
		// get_transport_parameter();
		setEcutRegion = false;
		setPcutRegion = false;
		setSmaxRegion = false;
		ecut = 0.521;
		pcut = 0.001;
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
		// incoh=incoh_ON;
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
		// coh=coh_OFF;
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
		// relax=relax_ON;
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
		// pe=pe_ang_ON;
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
		// EGS4.eii_flag=eii_OFF;
		// Triplet
		// EGS4.itriplet=triplet_OFF;
		// Radiative compton correction
		// EGS4.radc_flag=radc_OFF;
		// ###########################################################################################
		// VARIANCE REDUCTION:
		// @BREM SPLITTING
		// EGS4.nbr_split=1;//$MAXBRSPLIT
		// EGS4.i_play_RR = 0;
		// ELECTRON RANGE REJECTION
		// EGS4Macro.irejct=irejct_OFF;
		// #On: if charged particle energy is below ESAVEIN
		// # and it cannot get out of current region
		// # with energy > ECUT, the particle is
		// # terminated
		// #also terminates all electrons which cannot
		// #reach the cavity under conservative assumptions.
		// ESAVEIN=2.0;//#total energy below which range rejection is considered

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
		// IFARCE=IFARCE_OFF;
		// EGS4Macro.NFMIN=1;//#Start forcing at this interaction number
		// EGS4Macro.NFMAX=1;//#Number of photon interactions after which
		// #to stop forcing photon interactions
		// PHOTON SPLITTING= 1 #no. of times to split a photon
		// #if < 2-->normal transport
		// #overrides PHOTON FORCING if >= 2
		// #can only be >= 2 if IFULL= dose and stoppers
		// # or if IFULL= Aatt and Ascat
		// phsplitt=1;

		// EGS4Macro.cs_enhance=1.0;// #Photon cross section scaling factors
		// nENHREG=1;
		// NENHLO[0]=1;
		// NENHHI[0]=EGS4.$MXREG;//all geometry
		// NENHLO[1]=1;
		// NENHHI[1]=1;
		// CS ENHANCEMENT FACTOR= 1
		// CS ENHANCEMENT START REGION= 1, 1
		// CS ENHANCEMENT STOP REGION= 1, 1

		// PLOTSN();

		test_inputs();
		if (EGS4.STOPPROGRAM) {
			closeFile();
			destroyArrays();
			return;
		}
		// ##########################################################################################
		// ##########################################################################################

		// "CALCULATE THE NUMBER OF DOSE COMPONENTS"
		// if((EGS4Macro.IFULL==0)||(EGS4Macro.IFULL==2)||(EGS4Macro.IFULL==3)){ITMAX=2;}
		// else{ITMAX=$MAXIT;}
		// if(IPRIMARY == 0){ITMAX=1;}else{ITMAX=$MAXIT;}
		ICROSS = 0;
		ITRACKE = 0;
		NCASEO = 0;
		NCASET = 0;
		TMCPUO = 0;
		SCSTP_TMP = 0;
		if (IRESTART == 0 || IRESTART == 5)// ALWAYS 0!!!!
		{// "0-FRESH START, SET EVERYTHING TO"
			// "ZERO OR 5- INITIALIZE ALL ARRAYS TO 0 FOR PARALLEL POST-PROCESSING"

			PIISTP = 0;// SCSTP=0.;SCSTP2=0.;SCSTP_LAST=0;
			SCSTP = 0;
			SCSTP2 = 0;// SCSTP_TMP=0;
			SCSTP_LAST = 0;
			// SCDSTP=0; SCDSTP2=0; SCDSTP_TMP=0; SCDSTP_LAST=0;

			for (int IX = 1; IX <= EGS4Geom.NR; IX++)// ( IX=1,NR[
			{
				for (int IZ = 1; IZ <= EGS4Geom.NZ; IZ++)// DO IZ=1,NZ[
				{
					for (int IT = 1; IT <= $MAXIT; IT++)// DO IT=1,ITMAX[
					{
						SCEDPN[IZ - 1][IX - 1][IT - 1] = 0.;
						SCEDPN2[IZ - 1][IX - 1][IT - 1] = 0.;
						SCEDPN_TMP[IZ - 1][IX - 1][IT - 1] = 0.;
					}
					SCDOSE_TMP3[IZ - 1][IX - 1] = 0.;
					for (int IT = 1; IT <= $MAXDOS; IT++) {
						SCDOSE[IZ - 1][IX - 1][IT - 1] = 0.;
						SCDOSE2[IZ - 1][IX - 1][IT - 1] = 0.;
						if (IT < $MAXDOS)
							SCDOSE_COV[IZ - 1][IX - 1][IT - 1] = 0.;
					}
					for (int IT = 1; IT <= $NSPR; IT++) {
						SCSPR_COV[IZ - 1][IX - 1][IT - 1] = 0.;
					}
					SCDOSE_LAST[IZ - 1][IX - 1] = 0;
					SCNUM[IZ - 1][IX - 1] = 0;

				}
			}

			/*
			 * for(int IT=1;IT<=ITMAX;IT++) { for(int IX=1;IX<=EGS4Geom.NR;IX++)
			 * { for(int IZ=1;IZ<=EGS4Geom.NZ;IZ++) {
			 * SCDOSE[IZ-1][IX-1][IT-1]=0.0; SCDOSE2[IZ-1][IX-1][IT-1]=0.0;
			 * SCDOSE_TMP[IZ-1][IX-1][IT-1]=0.0;
			 * SCDOSE_LAST[IZ-1][IX-1][IT-1]=0; if(IKERMA==1) {
			 * SCKERMA[IZ-1][IX-1][IT-1]=0.0; SCKERMA2[IZ-1][IX-1][IT-1]=0.0;
			 * SCKERMA_TMP[IZ-1][IX-1][IT-1]=0.0;
			 * SCKERMA_LAST[IZ-1][IX-1][IT-1]=0;
			 * SCDOSEtoKERMA2[IZ-1][IX-1][IT-1]=0.0; } } } }
			 * if(EGS4Macro.IFULL==2) { for(int IB=1;IB<=MAXBIN;IB++) {
			 * SCPDST[IB-1]=0.0; SCPDST2[IB-1]=0.0; SCPCUM[IB-1]=0.0;
			 * SCPCUM2[IB-1]=0.0; //********************** ECUM[IB-1]=0.0;
			 * ECUM2[IB-1]=0.0; //************************** } for(int
			 * IPK=1;IPK<=4;IPK++) { SCDFEP[IPK-1]=0.0; SCDFEP2[IPK-1]=0.0;
			 * SCDFBK[IPK-1]=0.0; SCDFBK2[IPK-1]=0.0; SCDFDIFF[IPK-1]=0.0;
			 * SCDFDIFF2[IPK-1]=0.0; } SCPTOT=0.0; SCPTOT2=0.0; SCPDST_LAST=0; }
			 */
		}// "END OF IRESTART =0 OR 5"
		// NCASET=NCASE+NCASEO;
		NCASET = NCASE + NCASEO;
		EGS4SrcEns.NCASET = NCASET;

		// :FINISHED: CONTINUE;
		// "************************"
		// "* Check for any errors *"
		// "************************"
		EGS4.seqStr = " ********* SUCCESSFUL INPUT ACCOMPLISHED *********";
		if (EGS4.iprint > 1)
			printSequence(EGS4.seqStr);

	}// inputs

	/**
	 * Validate more inputs. Called by inputs routine.
	 */
	private void test_inputs() {
		// =================================
		// ==================================
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
		if (IOOPTN < 0 || IOOPTN > 4) {
			IOOPTN = 0;
		}// default
		if (EGS4Macro.ICSDA < 0 || EGS4Macro.ICSDA > 1) {
			EGS4Macro.ICSDA = 0;
		}// default
		// ================
		// if(NZDMIN<1 || NZDMIN>EGS4Geom.$MAXZREG){NZDMIN=1;}
		// if(NZDMAX<2 ||
		// NZDMAX>EGS4Geom.$MAXZREG+1){NZDMAX=EGS4Geom.$MAXZREG+1;}
		// if(NRDMIN<0 || NRDMIN>EGS4Geom.$MAXRADII-1){NRDMIN=0;}
		// if(NRDMAX<1 || NRDMAX>EGS4Geom.$MAXRADII){NRDMAX=EGS4Geom.$MAXRADII;}
		// NZDOSE=NZDMAX-NZDMIN;
		// NRDOSE=NRDMAX-NRDMIN;
		// ==============
		if (IDAT < 0 || IDAT > 1) {
			IDAT = 1;
		}// default
		if (IDAT < 1 || IDAT > 1) {
			IDAT = 1;
		}// default->NO FILE STORAGE
		// =============================================================
		// IF(IDAT.EQ.1)[INEXT=0;]ELSE[INEXT=1;]
		if (IFANO == 1) {
			// OUTPUT; (/' ********* Calculation WITH photon regeneration
			// ******** '/);
			EGS4.seqStr = " ********* Calculation WITH photon regeneration ******** ";
			if (EGS4.iprint > 2)
				printSequence(EGS4.seqStr);
		} else {
			// OUTPUT; (/' ********* Calculation WITHOUT photon regeneration
			// ******* '/);
			EGS4.seqStr = " ********* Calculation WITHOUT photon regeneration ******* ";
			if (EGS4.iprint > 2)
				printSequence(EGS4.seqStr);
		}
		// for(int
		// J=1;J<=EGS4Geom.nreg;J++){ISPR[J-1]=0;}//]"default to not output/plot stopping-power ratio "
		// =============================================================
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
		if (EGS4Macro.IFULL < 0 || EGS4Macro.IFULL > 4) {
			EGS4Macro.IFULL = 0;
		}// default
		if (STATLM < STATLM_MIN || STATLM > STATLM_MAX) {
			STATLM = STATLM_DEFAULT;
		}
		if (IFANO < 0 || IFANO > 2) {
			IFANO = 0;
		}// default
		EGS4Macro.iifano = IFANO;
		if (IWATCH == 0 && NCASE < $NCASEMIN) {
			NCASE = $NCASEMIN;
		}
		// if(IKERMA<0 || IKERMA>1){IKERMA=0;}//default

		if (EGS4Geom.iterseindex < 0 || EGS4Geom.iterseindex > 2) {
			EGS4.STOPPROGRAM = true;
			EGS4.seqStr = " ******** ERROR in geometrical inputs! **** ";
			printSequence(EGS4.seqStr);
			destroyArrays();
			return;
		}
		if (EGS4Geom.iterseindex == EGS4Geom.iCAVITY_INFORMATION)// 2
		{
			if ((EGS4Geom.SLENGHT).compareTo(NULLs) == 0
					|| (EGS4Geom.airs).compareTo(NULLs) == 0) {
				EGS4.STOPPROGRAM = true;
				EGS4.seqStr = " ******** ERROR in geometrical inputs! **** ";
				printSequence(EGS4.seqStr);
				destroyArrays();
				return;
			}
			if (EGS4Geom.ELERAD != 0.0
					&& (EGS4Geom.electrods).compareTo(NULLs) == 0) {
				EGS4.STOPPROGRAM = true;
				EGS4.seqStr = " ******** ERROR in geometrical inputs! **** ";
				printSequence(EGS4.seqStr);
				destroyArrays();
				return;
			}
		}
		/*
		 * if ( NZDMAX > EGS4Geom.NZ+1 ) { NZDMAX=EGS4Geom.NZ+1; //OUTPUT
		 * NZDMAX;(/'===> MAX. SCORING PLANE # RESET TO: ',I6);
		 * EGS4.seqStr=" ===> MAX. SCORING PLANE # RESET TO: "
		 * +EGS4.format(NZDMAX,6); if (EGS4.iprint>2)
		 * printSequence(EGS4.seqStr);
		 * 
		 * } if ( NZDMIN >= NZDMAX ) { NZDMIN = NZDMAX -1; //OUTPUT
		 * NZDMIN;(/'===> MIN. SCORING PLANE # RESET TO: ',I6);
		 * EGS4.seqStr=" ===> MIN. SCORING PLANE # RESET TO: "
		 * +EGS4.format(NZDMIN,6); if (EGS4.iprint>2)
		 * printSequence(EGS4.seqStr);
		 * 
		 * }
		 * 
		 * if ( NRDMAX > EGS4Geom.NR ) { NRDMAX=EGS4Geom.NR; // OUTPUT
		 * NRDMAX;(/'===> MAX. SCORING CYLINDER # RESET TO: ',I6);
		 * EGS4.seqStr=" ===> MAX. SCORING CYLINDER # RESET TO: "
		 * +EGS4.format(NRDMAX,6); if (EGS4.iprint>2)
		 * printSequence(EGS4.seqStr);
		 * 
		 * } if ( NRDMIN >= NRDMAX ) { NRDMIN = NRDMAX -1; // OUTPUT
		 * NRDMIN;(/'===> MIN. SCORING CYLINDER # RESET TO: ',I6);
		 * EGS4.seqStr=" ===> MIN. SCORING CYLINDER # RESET TO: "
		 * +EGS4.format(NRDMIN,6); if (EGS4.iprint>2)
		 * printSequence(EGS4.seqStr);
		 * 
		 * } NZDOSE=NZDMAX-NZDMIN; NRDOSE=NRDMAX-NRDMIN;
		 * //"-------------------------------------------------------------"
		 * //PULSE HEIGHT=============
		 * if(EGS4Macro.IFULL==IFULL_PULSE_HEIGHT_DISTRIBUTION) { //OUTPUT;(/'
		 * INPUT FOR PULSE HEIGHT DISTRIBUTION'/);
		 * //"INITIALIZE FLAGS TO NO PULSE HEIGHT DISTRIBUTION IN EACH REGION"
		 * for(int J=1;J<=EGS4Geom.nreg;J++) {IPHR[J-1]=0;} }
		 * if(EGS4Macro.IFULL==IFULL_PULSE_HEIGHT_DISTRIBUTION) { for(int
		 * I=1;I<=nREGSOLV;I++) if(REGSVOL[I-1]<0 ||
		 * REGSVOL[I-1]>EGS4Geom.nreg){REGSVOL[I-1]=EGS4Geom.nreg;}//default
		 * if(SLOTE<SLOTE_MIN || SLOTE>SLOTE_MAX){SLOTE=SLOTE_DEFAULT;}//default
		 * if(DELTAE<DELTAE_MIN ||
		 * DELTAE>DELTAE_MAX){DELTAE=DELTAE_DEFAULT;}//default
		 * 
		 * for(int J=1;J<= nREGSOLV;J++) { //REGNUM=VALUE(NUM_REGSVOL,J);
		 * IPHR[REGSVOL[J-1]-1]=1; //OUTPUT REGNUM,MED(REGNUM); (/T10,'
		 * REGION',I4,' HAS MEDIUM',I3);
		 * EGS4.seqStr=" 	 REGION"+EGS4.format(REGSVOL[J-1],4)+"  HAS MEDIUM"+
		 * EGS4.format(EGS4.MED[REGSVOL[J-1]-1],3); if (EGS4.iprint>2)
		 * printSequence(EGS4.seqStr);
		 * 
		 * }
		 */
		if (SLOTE < SLOTE_MIN || SLOTE > SLOTE_MAX) {
			SLOTE = SLOTE_DEFAULT;
		}// default
		if (DELTAE < DELTAE_MIN || DELTAE > DELTAE_MAX) {
			DELTAE = DELTAE_DEFAULT;
		}// default

		if (SLOTE > 0.0) {
			// OUTPUT SLOTE,SLOTE*dble($EBIN);
			// (/' EQUAL BINS OF',F10.4,' MeV WILL COVER UP TO',F10.3,' MeV');
			EGS4.seqStr = " EQUAL BINS OF" + EGS4.format(SLOTE, 10, true)
					+ " MeV WILL COVER UP TO"
					+ EGS4.format(SLOTE * EGS4SrcEns.dble($EBIN), 10, true)
					+ " MeV";
			if (EGS4.iprint > 2)
				printSequence(EGS4.seqStr);

			// "NOTE THAT LATER, WHEN WE KNOW THE MAXIMUM ENERGY IN THE"
			// "INPUT SPECTRUM, WE WILL INCREASE SLOTE TO MAKE SURE IT WORKS"
		} else {// "SLOTE <= 0.0 means we input energy bins"

			for (int I = 1; I <= nTOPEBIN; I++)
				if (TOPEBIN[I - 1] < TOPEBIN_MIN
						|| TOPEBIN[I - 1] > TOPEBIN_MAX) {
					TOPEBIN[I - 1] = TOPEBIN_DEFAULT;
				}// default

			// $GET_INPUT(NUM_TOPEBIN);
			// DO J=1, NVALUE(NUM_TOPEBIN) [BINTOP(J)=VALUE(NUM_TOPEBIN,J);]
			for (int J = 1; J <= nTOPEBIN; J++) {
				BINTOP[J - 1] = TOPEBIN[J - 1];
			}
			MAXBIN = nTOPEBIN;// NVALUE(NUM_TOPEBIN);

			if (nTOPEBIN > $EBIN + 1) // (NVALUE(NUM_TOPEBIN) > $EBIN+1)
			{
				EGS4.STOPPROGRAM = true;
				EGS4.seqStr = " Tried to use more than max number of energy bins";
				printSequence(EGS4.seqStr);
				EGS4.seqStr = " Either increase $EBIN or reduce number of bins";
				printSequence(EGS4.seqStr);
				destroyArrays();

				return;

				// OUTPUT $EBIN;
				// (/' ****Tried to use more than max number of energy bins
				// ',I5/
				// ' Either increase $EBIN or reduce number of bins');
				// STOP;
			}
			// OUTPUT NVALUE(NUM_TOPEBIN),(BINTOP(J),J=1,NVALUE(NUM_TOPEBIN));
			// (/' READ A TOTAL OF',I4,' ENERGY BINS'/ (5F15.4) );
		}

		// }
		/*
		 * if(EGS4Macro.IFULL == 3 && EGS4SrcEns.iqin!= 0) {//
		 * "only allow option 3 for pure photon beam" //
		 * "                           could modify for full phase space if needed"
		 * //OUTPUT; //(// 1x,70('@')/' Changed IFULL to 0 from 3 since photon
		 * beam not input'); EGS4Macro.IFULL = 0;
		 * EGS4.seqStr="  Changed IFULL to 0 from 3 since photon beam not input"
		 * ; if (EGS4.iprint>1) printSequence(EGS4.seqStr);
		 * 
		 * }
		 */
		double EPHTOP = 0.0;
		if (EGS4Macro.IFULL == IFULL_PULSE_HEIGHT_DISTRIBUTION) {
			// "DO SOME CHECKS ON PULSE HEIGHT DISTRIBUTION BINS"
			if (EGS4SrcEns.iqin == 1) {
				EPHTOP = EGS4SrcEns.ein + 1.022;// "INCLUDE ANNIHILATION FOR POSITRONS IN"
			} else {
				EPHTOP = EGS4SrcEns.ein;
			}
			if (SLOTE > 0.0) {
				while (true) {
					if (SLOTE * EGS4SrcEns.dble($EBIN) > 1.05 * EPHTOP)
						break;
					SLOTE = SLOTE * 2.;
					EGS4.seqStr = "  HAVE DOUBLED SLOTE TO"
							+ EGS4.format(SLOTE, 12, true)
							+ " MeV TO REACH MAXIMUM INPUT ENERGY";
					if (EGS4.iprint > 1)
						printSequence(EGS4.seqStr);
				}
				Double dbb = new Double(1.05 * EPHTOP / SLOTE + 1.0);
				MAXBIN = dbb.intValue();// IFIX(1.05*EPHTOP/SLOTE+0.999);
				EGS4.seqStr = " MAXBIN SET TO " + EGS4.format(MAXBIN, 3)
						+ " TO COVER SPECTRUM";
				if (EGS4.iprint > 1)
					printSequence(EGS4.seqStr);

				// OUTPUT MAXBIN;(/' MAXBIN SET TO',I3,' TO COVER SPECTRUM'/);
				// ================
				// "set up BINTOP array needed within AUSGAB to allocate across bins"
				// "recall MAXIB was set so just greater than max input covered."
				// DO IB=1,MAXIB [BINTOP(IB) =SLOTE*dble(IB);]
				for (int I = 1; I <= MAXBIN; I++) {
					BINTOP[I - 1] = SLOTE * I;// dble(IB);
				}
				// ================
			}// "END OF SLOTE>0 BLOCK"
			else {
				if (BINTOP[MAXBIN - 1] <= EPHTOP) {
					EGS4.seqStr = " CHANGED BINTOP(" + EGS4.format(MAXBIN, 3)
							+ "-1)" + " TO REACH"
							+ EGS4.format(EPHTOP, 10, true) + " MeV";
					if (EGS4.iprint > 1)
						printSequence(EGS4.seqStr);

					// OUTPUT MAXBIN,EPHTOP;
					// (/' ***CHANGED BINTOP(',I3,') TO REACH', F10.3,'
					// MeV***');
					BINTOP[MAXBIN - 1] = EPHTOP;
				}
			}// "END SLOTE<=0 BLOCK"
				// =========================
			MAXIB = MAXBIN;
			// ============================
			/*
			 * DFEN[0][1]=EPHTOP-DELTAE;DFEN[1][1]=DFEN[0][1]-0.511;
			 * DFEN[2][1]=DFEN[0][1]-1.022;DFEN[3][1]=0.511-DELTAE;
			 * //"I.E. WE SET LOWER ENERGIES FOR FULL ENERGY, SINGLE ESCAPE, DOUBLE"
			 * //"ESCAPE AND THE 511 KEV LINE" for(int IPK = 1;IPK<=4;IPK++) {
			 * DFEN[IPK-1][0] = DFEN[IPK-1][1] - DELTAE; DFEN[IPK-1][2] =
			 * DFEN[IPK-1][1] +2*DELTAE; DFEN[IPK-1][3] = DFEN[IPK-1][2]
			 * +DELTAE; }
			 */

		}

		// =========================

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
		// =================================
		if ((EGS4.nbr_split > $MAXBRSPLIT))
			EGS4.nbr_split = $MAXBRSPLIT;
		// CHARGED PARTICLE RUSSIAN ROULETTE
		if (EGS4.i_play_RR == 1) {
			EGS4.prob_RR = 1. / EGS4SrcEns.dble(EGS4.nbr_split);
		} else {
			EGS4.prob_RR = 1.;
		}
		/*
		 * if(EGS4Macro.IFULL==2 && EGS4.nbr_split > 1) {// "cannot have this"
		 * EGS4.seqStr=" You cannot calculate a pulse height distribution with";
		 * if (EGS4.iprint>1) printSequence(EGS4.seqStr);
		 * EGS4.seqStr=" bremsstrahlung splitting on.  Will run simulation with"
		 * ; if (EGS4.iprint>1) printSequence(EGS4.seqStr);
		 * EGS4.seqStr=" IFULL= dose and stoppers"; if (EGS4.iprint>1)
		 * printSequence(EGS4.seqStr);
		 * 
		 * // WRITE(1,'(//'' ****WARNING****''/ // ''You cannot calculate a
		 * pulse height distribution with''/ // '' bremsstrahlung splitting on.
		 * Will run simulation with''/ // '' IFULL= dose and stoppers''//)');
		 * EGS4Macro.IFULL=0; }
		 */

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
		// if( EGS4Macro.cs_enhance > 1. ) { EGS4Macro.use_enhance = true; }
		// else { EGS4Macro.use_enhance = false; }

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
		/*
		 * if(EGS4Macro.IFULL==2 && RUSROU) {
		 * EGS4.seqStr=" You cannot calculate a pulse height distribution with";
		 * if (EGS4.iprint>1) printSequence(EGS4.seqStr);
		 * EGS4.seqStr=" Russian Roulette on.  Will run simulation with"; if
		 * (EGS4.iprint>1) printSequence(EGS4.seqStr);
		 * EGS4.seqStr=" IFULL= dose and stoppers"; if (EGS4.iprint>1)
		 * printSequence(EGS4.seqStr);
		 * 
		 * EGS4Macro.IFULL=0; } if(EGS4Macro.IFULL==2 && EGS4Macro.CEXPTR!=0) {
		 * EGS4.seqStr=" You cannot calculate a pulse height distribution with";
		 * if (EGS4.iprint>1) printSequence(EGS4.seqStr);
		 * EGS4.seqStr=" pathlength biasing.  Will run simulation with"; if
		 * (EGS4.iprint>1) printSequence(EGS4.seqStr);
		 * EGS4.seqStr=" IFULL= dose and stoppers"; if (EGS4.iprint>1)
		 * printSequence(EGS4.seqStr);
		 * 
		 * EGS4Macro.IFULL=0; }
		 */

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

		/*
		 * if(EGS4Macro.IFULL==2 && (EGS4Macro.NFMAX>EGS4Macro.NFMIN ||
		 * EGS4Macro.NFMIN>1)) {
		 * EGS4.seqStr=" You cannot calculate a pulse height distribution with";
		 * if (EGS4.iprint>1) printSequence(EGS4.seqStr); EGS4.seqStr=
		 * " more than 1 interaction forced or if the forced interaction"; if
		 * (EGS4.iprint>1) printSequence(EGS4.seqStr);
		 * EGS4.seqStr=" is > the first interaction.  Will run simulation with";
		 * if (EGS4.iprint>1) printSequence(EGS4.seqStr);
		 * EGS4.seqStr=" IFULL= dose and stoppers"; if (EGS4.iprint>1)
		 * printSequence(EGS4.seqStr);
		 * 
		 * EGS4Macro.IFULL=0; }
		 */
		ics_enhance = 0;
		for (int jj = 1; jj <= EGS4Geom.nreg; jj++) {
			EGS4Macro.iefl[jj - 1] = 0;
		}
		for (int ii = 1; ii <= nENHREG; ii++) {
			int ics_start = NENHLO[ii - 1];
			int ics_stop = NENHHI[ii - 1];
			for (int jj = ics_start; jj <= ics_stop; jj++) {
				EGS4Macro.iefl[jj - 1] = 1;
			}
		}

		int COUNT = 0;
		for (int jj = 2; jj <= EGS4Geom.nreg; jj++) {
			COUNT = COUNT + EGS4Macro.iefl[jj - 1];
		}
		// "We don't care about region 1 since outside geometry"
		if (COUNT > 0 && (EGS4Macro.cs_enhance > 1.0001)) {// "there is enhancement somewhere"
			EGS4.seqStr = " Cross section enhancement in regions ";
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);
			for (int jj = 1; jj <= EGS4Geom.nreg; jj++) {
				EGS4.seqStr = " 		" + EGS4.format(EGS4Macro.iefl[jj - 1], 4);
				if (EGS4.iprint > 1)
					printSequence(EGS4.seqStr);
			}
			EGS4.seqStr = " Cross section enhancement factor: "
					+ EGS4.format(EGS4Macro.cs_enhance, 6, true);
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);

			// OUTPUT; ('Cross section enhancement in regions ');
			// OUTPUT (iefl(jj),jj=1,NREG); (20 I4);
			// OUTPUT cs_enhance;
			// ( ' Cross section enhancement factor: ',T60,F6.1/);
			ics_enhance = 1;
		} else {
			EGS4.seqStr = " No cross section enhancement";
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);

			ics_enhance = 0;
		}

		/*
		 * if(EGS4Macro.IFULL==2 && ics_enhance==1) {
		 * EGS4.seqStr=" You cannot calculate a pulse height distribution with";
		 * if (EGS4.iprint>1) printSequence(EGS4.seqStr);
		 * EGS4.seqStr=" cross section enhancement.  Will run simulation with";
		 * if (EGS4.iprint>1) printSequence(EGS4.seqStr);
		 * EGS4.seqStr=" IFULL= dose and stoppers"; if (EGS4.iprint>1)
		 * printSequence(EGS4.seqStr); EGS4Macro.IFULL=0; }
		 */
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
			ll = 61 - ll;
			s = s + EGS4.format("", ll);
			EGS4.seqStr = s + "OFF";
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);
		}
		if (ics_enhance == 1) {
			// WRITE(IOUT,235) cs_enhance,(j,iefl(j),j=2,NREG);
			// 235 FORMAT(T20,'Cross section enhancement factor of',T56,F8.1/
			// T20,'In regions with a 1:'/
			// (T10, 10('(',I3,',',I1,')')));
			s = " Cross section enhancement factor of";
			ll = s.length();
			ll = 57 - ll;
			s = s + EGS4.format("", ll);
			EGS4.seqStr = s + EGS4.format(EGS4Macro.cs_enhance, 8, true);
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);

		} else {
			s = " No cross section enhancement used";
			// ll=s.length();ll=56-ll;s=s+EGS4.format("",ll);
			EGS4.seqStr = s;// +"OFF";
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);

			// WRITE(IOUT,236);236 FORMAT(T20,'No cross section enhancement
			// used');
		}

		if (EGS4Macro.irejct > 0) {
			s = " RANGE REJECTION SWITCH";
			ll = s.length();
			ll = 61 - ll;
			s = s + EGS4.format("", ll);
			EGS4.seqStr = s + "ON";
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);
			s = "  Range rejection for energy  <";
			ll = s.length();
			ll = 57 - ll;
			s = s + EGS4.format("", ll);
			EGS4.seqStr = s + EGS4.format(ESAVEIN, 9, true) + " MeV";
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);
			s = "  Ranges determined internally from stopping powers";
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
		for (int I = 1; I <= EGS4Geom.nreg; I++) {
			if (EGS4.IRAYLR[I - 1] == 1) {
				// WRITE(IOUT,244);
				s = " RAYLEIGH SCATTERING INCLUDED";
				// ll=s.length();ll=61-ll;s=s+EGS4.format("",ll);
				EGS4.seqStr = s;// +"ON";
				if (EGS4.iprint > 1)
					printSequence(EGS4.seqStr);

				break;// EXIT;
			}
		}

		// WRITE(IOUT,260) TIMMAX,STATLM;
		s = " Maximum cputime allowed";
		ll = s.length();
		ll = 60 - ll;
		s = s + EGS4.format("", ll);
		EGS4.seqStr = s + EGS4.format(TIMMAX, 6, true) + " hrs";
		if (EGS4.iprint > 1)
			printSequence(EGS4.seqStr);
		s = " Stats in PEAK REGION objective";
		ll = s.length();
		ll = 61 - ll;
		s = s + EGS4.format("", ll);
		EGS4.seqStr = s + EGS4.format(STATLM, 6, true) + " %";
		if (EGS4.iprint > 1)
			printSequence(EGS4.seqStr);

		if (SOURCE == -1) {
			EGS4.seqStr = " Initial RNG state:";
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);
			EGS4.SHOW_RNG_STATE();
		}
		if (EGS4.nbr_split == 1) {
			s = " Bremsstrahlung splitting";
			ll = s.length();
			ll = 61 - ll;
			s = s + EGS4.format("", ll);
			EGS4.seqStr = s + "OFF";
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);

			// WRITE(IOUT,312);
		} else {
			// WRITE(IOUT,313) nbr_split;
			s = " Bremsstrahlung splitting";
			ll = s.length();
			ll = 61 - ll;
			s = s + EGS4.format("", ll);
			EGS4.seqStr = s + "ON";
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);
			s = " Initially, each bremsstrahlung photon split into ";
			EGS4.seqStr = s + EGS4.format(EGS4.nbr_split, 3) + " photons";
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);

		}
		if (EGS4.i_play_RR == 0) {
			// WRITE(IOUT,314);
			s = " Charged particle Russian Roulette";
			ll = s.length();
			ll = 61 - ll;
			s = s + EGS4.format("", ll);
			EGS4.seqStr = s + "OFF";
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);

		} else {
			// WRITE(IOUT,315) EGS4.prob_RR;
			s = " Charged particle Russian Roulette";
			ll = s.length();
			ll = 61 - ll;
			s = s + EGS4.format("", ll);
			EGS4.seqStr = s + "ON";
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);
			s = " With probability of survival = ";
			ll = s.length();
			ll = 61 - ll;
			s = s + EGS4.format("", ll);
			EGS4.seqStr = s + EGS4.format(EGS4.prob_RR, 9);
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);
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
		if (EGS4Macro.ICSDA == 0) {
			// WRITE(IOUT,289);
		} else {
			// WRITE(IOUT,290);
		}
		EK0 = EGS4SrcEns.ein;
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
		if (SOURCE == -1)
			EGS4Geom.GEOMRZ_ISUMRY();

		if (SOURCE == -1)
			EGS4SrcEns.SRCOUT();
		// "       PRINT A GRID OF THE ZONE DEPENDENT VARIABLES"
		// "       ============================================"
		// CDSTBL(1)='0';CTRTBL(1)='0';CABSRB(1)='0';
		// CDSTBL(1)='0';CTRTBL(1)='0';CABSRB(1)='0';
		EGS4Grid.CABSRB[0] = "0";
		EGS4Grid.CDSTBL[0] = "0";
		EGS4Grid.CTRTBL[0] = "0";
		// "Make the material grid"
		EGS4Grid.MATERIALGRID(EGS4Geom.NR, EGS4Geom.NZ, AMASS, 1, EGS4.ECUT,
				EGS4.PCUT, EGS4Geom.RCYL, EGS4Geom.ZPLANE, EGS4.MED, EGS4.MEDIA);

		// EGS4Grid.MATERIALGRID(EGS4Geom.NR, EGS4Geom.NZ, VOLUME, 2, EGS4.ECUT,
		// EGS4.PCUT, EGS4Geom.RCYL, EGS4Geom.ZPLANE, EGS4.MED, EGS4.MEDIA);

		/*
		 * if(EGS4Macro.IFULL==2) {
		 * s=" PULSE HEIGHT DISTRIBUTION WILL BE SCORED IN DETECTOR VOLUME!";//
		 * IN THOSE REGIONS DENOTED WITH A 1"; EGS4.seqStr=s; if(EGS4.iprint>1)
		 * printSequence(EGS4.seqStr);
		 * 
		 * }
		 */// "END IFULL=2 BLOCK"

	}

	@SuppressWarnings("unused")
	/**
	 * Not used here.
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

		//int I = 0;
		int IRL = 0;// int IX=0;int IZ=0;
		//double ASCT = 0.0;
		//double ASCTUN = 0.0;
		//double AATT = 0.0;
		//double AATTUN = 0.0;
		//double AWLL = 0.0;
	//double AWLLUN = 0.0;
	//	double TDAW = 0.0;
	//	double TDAWUN = 0.0;
	//	double KSCT = 0.0;
	//	double KATT = 0.0;
	//	double KWLL = 0.0;
	//	double KSCTUN = 0.0;
	//	double KATTUN = 0.0;
	//	double KWLLUN = 0.0;

	//	int IPLTUNX = 0;
	//	int SPRRZ_COUNT = 0;
		// int BIN_CUT,
	//	int ICOUNT = 0;
	//	IPLTUNX = 23; // "UNIT number used for optput file *.plotdat"
	//	SPRRZ_COUNT = 0;// "keeps count of curvenumber in stopping power calculations"

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
		// "SCALE DOSE FRACTIONS"
	//	double TEMP = 0.0;
		// FRACS($MAXZREG,$MAXRADII,3:6);
	//	double[][][] FRACS = new double[EGS4Geom.$MAXZREG][EGS4Geom.$MAXRADII][6];
		// "PRINT A SUMMARY OF THE DOSE REGION RESULTS"
		IRL = 0;
		// =================================================================
		EGS4Grid.RESULTS = new double[EGS4Geom.$MAXZREG][EGS4Geom.$MAXRADII][$MAXCMPTS];
		EGS4Grid.UNCRT = new double[EGS4Geom.$MAXZREG][EGS4Geom.$MAXRADII][$MAXCMPTS];// ($MAXZREG,
																						// $MAXRADII,
																						// $MAXCMPTS),
		EGS4Grid.RADIAL_BINS = new double[EGS4Geom.$MAXRADII];
		EGS4Grid.DEPTH_BINS = new double[EGS4Geom.$MAXZPLANE];

		// if(EGS4Geom.nreg > 2)//if((IOOPTN==0)||(IOOPTN==2)||(IOOPTN==4))
		// {
		// " OUTPUT GRID ONLY IF REQUESTED"
		// "Grid routine by Aaron Merovitz, 1998"

		// "1) Set up the arrays"
		for (int IX = 1; IX <= EGS4Geom.NR; IX++) {
			for (int IZ = 1; IZ <= EGS4Geom.NZ; IZ++) {
				// $GET-IRL(IZD,IXD);
				IRL = EGS4Geom.GET_IRL(IZ, IX);
				if (EGS4.i_do_rr[IRL - 1] == 1) {
					SCEDPN[IZ - 1][IX - 1][0] = 0.;
					SCEDPN2[IZ - 1][IX - 1][0] = 99.9;
					SCEDPN[IZ - 1][IX - 1][2] = 0.;
					SCEDPN2[IZ - 1][IX - 1][2] = 99.9;
				}
				EGS4Grid.RESULTS[IZ - 1][IX - 1][0] = SCEDPN[IZ - 1][IX - 1][0];
				EGS4Grid.UNCRT[IZ - 1][IX - 1][0] = SCEDPN2[IZ - 1][IX - 1][0];
				EGS4Grid.RESULTS[IZ - 1][IX - 1][1] = SCEDPN[IZ - 1][IX - 1][2];
				EGS4Grid.UNCRT[IZ - 1][IX - 1][1] = SCEDPN2[IZ - 1][IX - 1][2];

			}
		}
		// "2) Determine the number of components"
		// IF(IPRIM ~= 0) [NCOMP = 12;] ELSE [NCOMP = 6;]
		// if (IPRIMARY!=0){EGS4Grid.NCOMP=12;}else{EGS4Grid.NCOMP=6;}
		EGS4Grid.NCOMP = 2;
		// "3) Set up the bin indicators"
		// DO IX=1,NR+1 [RADIAL_BINS(IX)=RCYL(IX-1);]
		// DO IZ=1,NZ+1 [DEPTH_BINS(IZ)=ZPLANE(IZ);]
		for (int IX = 1; IX <= EGS4Geom.NR + 1; IX++) {
			EGS4Grid.RADIAL_BINS[IX - 1] = EGS4Geom.RCYL[IX - 1];
		}
		for (int IZ = 1; IZ <= EGS4Geom.NZ + 1; IZ++) {
			EGS4Grid.DEPTH_BINS[IZ - 1] = EGS4Geom.ZPLANE[IZ - 1];
		}

		// "4) Set up the labels and the explanations"
		// IF(ISOURC=3|ISOURC=21|ISOURC=22|ISOURC=23)[
		// EXPLANATIONS(1)='Total dose (Gray/incident particle)';
		// ]
		// ELSE[
		EGS4Grid.EXPLANATIONS[0] = "spr local medium/medium (DETECTOR) (with Track End)";
		// ]
		EGS4Grid.LABELS[0] = "spr:";

		EGS4Grid.EXPLANATIONS[1] = "spr sans TRACK ENDS";
		EGS4Grid.LABELS[1] = "ste:";

		// "5) Make the grid"
		// EGS4Grid.ZONEGRID(NRDOSE, NZDOSE, NRDMIN, NZDMIN, EGS4Geom.NZ,
		// RESULTS,
		// UNCRTY, NCOMP, RADIAL_BINS, DEPTH_BINS);//, LABELS, EXPLANATIONS);
		// System.out.println(NCOMP);
		if (EGS4Macro.irejct == 1) {
			// WRITE(IOUT,'(//'' ***Range rejection is ON:''/
			// '' stopping power ratios in those regions where the user''/
			// '' did not request output are zeroed below''/)');
			s = "Range rejection is ON. stopping power ratios in those regions are zeroed below";
			EGS4.seqStr = s;
			if (EGS4.iprint > 0)
				printSequence(EGS4.seqStr);

		}// ]

		EGS4Grid.ZONEGRID(EGS4Geom.NR, EGS4Geom.NZ, 0, 1, EGS4Geom.NZ);// ,
																		// EGS4Grid.NCOMP);

		// ==========================================================
		// write(iout,100) title; call egs_fdate(iout); write(iout,*);
		// " **** IK: using wrapper routines for date and time functions."
		// " Was: WRITE(IOUT,100) TITLE,DATEN,TIMEN;
		// WRITE(IOUT,320);
		s = "***** spr with track ends *****";
		EGS4.seqStr = s;
		if (EGS4.iprint > 0)
			printSequence(EGS4.seqStr);

		// WRITE(IOUT,305);
		/*
		 * 305 FORMAT(/' Region',T19,'depth (cm)',T35, ' spr (medium /
		 * med2)',T58,'dose/inc flu. Gy.cm**2'/);
		 */
		s = " Region";
		ll = s.length();
		ll = 19 - ll;
		s = s + EGS4.format("", 19);
		s = s + "depth (cm)";
		// ll=s.length();ll=35-ll;s=s+EGS4.format("",35);
		s = s + "	spr (medium / medDet)";
		// ll=s.length();ll=60-ll;s=s+EGS4.format("",60);
		s = s + "	dose/inc flu. Gy.cm**2";
		EGS4.seqStr = s;
		if (EGS4.iprint > 0)
			printSequence(EGS4.seqStr);

		for (int IX = 1; IX <= EGS4Geom.NR; IX++)// DO IX=1,NR[
		{
			for (int IZ = 1; IZ <= EGS4Geom.NZ; IZ++)// DO IZ=1,NZ[
			{
				// WRITE(IOUT,310) IZ,IX,ZPLANE(IZ),ZPLANE(IZ+1),
				// SCEDPN(IZ,IX,1),(SCEDPN(IZ,IX,1)*SCEDPN2(IZ,IX,1)/100.),
				// SCDOSE(IZ,IX,3), SCDOSE2(IZ,IX,3) ;
				/*
				 * 310 FORMAT(' NZ,NR=',I3,',',I3,F8.3,' -',F8.3,T33, F10.5,'
				 * (+/-',F8.5,')', T60,1PE10.3,'(+/-',0PF4.1,'%)');
				 */
				s = " NZ,NR= "
						+ EGS4.format(IZ, 3)
						+ ","
						+ EGS4.format(IX, 3)
						+ "  "
						+ EGS4.format(EGS4Geom.ZPLANE[IZ - 1], 8, true)
						+ " - "
						+ EGS4.format(EGS4Geom.ZPLANE[IZ], 8, true)
						+ "  "
						+ EGS4.format(SCEDPN[IZ - 1][IX - 1][0], 10, true)
						+ " (+/-"
						+ EGS4.format(SCEDPN[IZ - 1][IX - 1][0]
								* SCEDPN2[IZ - 1][IX - 1][0] / 100., 8, true)
						+ ")  "
						+ EGS4.format(SCDOSE[IZ - 1][IX - 1][2], 10, false)
						+ " (+/- "
						+ EGS4.format(SCDOSE2[IZ - 1][IX - 1][2], 6, true)
						+ "%)";
				EGS4.seqStr = s;
				if (EGS4.iprint > 0)
					printSequence(EGS4.seqStr);

				if (BSTAT[IZ - 1][0] == -1)// (BSTAT(IZ,1) = -1)
				{
					// WRITE(IOUT,315); **** WARNING **** BAD STATISTICS!
					s = "**** WARNING **** BAD STATISTICS!";
					EGS4.seqStr = s;
					if (EGS4.iprint > 0)
						printSequence(EGS4.seqStr);

				}
			}// ]
		}// ]
			// IF( (NR*NZ) > 23)["CALL PRNTER(ICHPIN,ILPIN,IOUT,IPAGE);"
			// "WRITE THE HEADER"
			// write(iout,100) title; call egs_fdate(iout); write(iout,*);
			// " **** IK: using wrapper routines for date and time functions."
			// " Was: WRITE(IOUT,100) TITLE,DATEN,TIMEN;
			// ]
			// ELSE [WRITE(IOUT,'(/)');]
			// WRITE(IOUT,325);***** spr without track ends *****
		s = "***** spr without track ends *****";
		EGS4.seqStr = s;
		if (EGS4.iprint > 0)
			printSequence(EGS4.seqStr);

		// WRITE(IOUT,305);
		s = " Region";
		ll = s.length();
		ll = 19 - ll;
		s = s + EGS4.format("", 19);
		s = s + "depth (cm)";
		// ll=s.length();ll=35-ll;s=s+EGS4.format("",35);
		s = s + "	 spr (medium / medDet)";
		// ll=s.length();ll=60-ll;s=s+EGS4.format("",60);
		s = s + "	dose/inc flu. Gy.cm**2";
		EGS4.seqStr = s;
		if (EGS4.iprint > 0)
			printSequence(EGS4.seqStr);

		// DO IX=1,NR[
		// DO IZ=1,NZ[
		// WRITE(IOUT,310) IZ,IX,ZPLANE(IZ),ZPLANE(IZ+1),
		// SCEDPN(IZ,IX,3),(SCEDPN(IZ,IX,3)*SCEDPN2(IZ,IX,3)/100.);
		// IF (BSTAT(IZ,2) = -1) [WRITE(IOUT,315);]
		// ]
		// ]
		for (int IX = 1; IX <= EGS4Geom.NR; IX++)// DO IX=1,NR[
		{
			for (int IZ = 1; IZ <= EGS4Geom.NZ; IZ++)// DO IZ=1,NZ[
			{
				// WRITE(IOUT,310) IZ,IX,ZPLANE(IZ),ZPLANE(IZ+1),
				// SCEDPN(IZ,IX,1),(SCEDPN(IZ,IX,1)*SCEDPN2(IZ,IX,1)/100.),
				// SCDOSE(IZ,IX,3), SCDOSE2(IZ,IX,3) ;
				/*
				 * 310 FORMAT(' NZ,NR=',I3,',',I3,F8.3,' -',F8.3,T33, F10.5,'
				 * (+/-',F8.5,')', T60,1PE10.3,'(+/-',0PF4.1,'%)');
				 */
				s = " NZ,NR= "
						+ EGS4.format(IZ, 3)
						+ ","
						+ EGS4.format(IX, 3)
						+ "  "
						+ EGS4.format(EGS4Geom.ZPLANE[IZ - 1], 8, true)
						+ " - "
						+ EGS4.format(EGS4Geom.ZPLANE[IZ], 8, true)
						+ "  "
						+ EGS4.format(SCEDPN[IZ - 1][IX - 1][2], 10, true)
						+ " (+/-"
						+ EGS4.format(SCEDPN[IZ - 1][IX - 1][2]
								* SCEDPN2[IZ - 1][IX - 1][2] / 100., 8, true)
						+ ")  "
				// EGS4.format(SCDOSE[IZ-1][IX-1][2],10,false)+" (+/- "+
				// EGS4.format(SCDOSE2[IZ-1][IX-1][2],6,true)+"%)"
				;
				EGS4.seqStr = s;
				if (EGS4.iprint > 0)
					printSequence(EGS4.seqStr);

				if (BSTAT[IZ - 1][1] == -1)// (BSTAT(IZ,1) = -1)
				{
					// WRITE(IOUT,315); **** WARNING **** BAD STATISTICS!
					s = "**** WARNING **** BAD STATISTICS!";
					EGS4.seqStr = s;
					if (EGS4.iprint > 0)
						printSequence(EGS4.seqStr);

				}
			}// ]
		}// ]

		// "CALL PRNTER(ICHPIN,ILPIN,IOUT,IPAGE);"
		// "WRITE THE HEADER"
		// write(iout,100) title; call egs_fdate(iout); write(iout,*);
		// " **** IK: using wrapper routines for date and time functions."
		// " Was: WRITE(IOUT,100) TITLE,DATEN,TIMEN;
		// WRITE(IOUT,350);
		s = "**** percentage of dose in local medium from stoppers ****";
		EGS4.seqStr = s;
		if (EGS4.iprint > 0)
			printSequence(EGS4.seqStr);
		s = " Included fraction of stoppers includes those stoppers that come from";
		EGS4.seqStr = s;
		if (EGS4.iprint > 0)
			printSequence(EGS4.seqStr);
		s = "                   electrons/positrons slowing down across delta.";
		EGS4.seqStr = s;
		if (EGS4.iprint > 0)
			printSequence(EGS4.seqStr);
		s = " Excluded fraction of stoppers are those electrons created below delta";
		EGS4.seqStr = s;
		if (EGS4.iprint > 0)
			printSequence(EGS4.seqStr);
		s = "                   and photons created below AP";
		EGS4.seqStr = s;
		if (EGS4.iprint > 0)
			printSequence(EGS4.seqStr);
		s = " Region";
		ll = s.length();
		ll = 19 - ll;
		s = s + EGS4.format("", 19);
		s = s + "depth (cm)";
		// ll=s.length();ll=20-ll;s=s+EGS4.format("",20);
		s = s + "		 included 	        excluded";
		EGS4.seqStr = s;
		if (EGS4.iprint > 0)
			printSequence(EGS4.seqStr);

		// DO IX=1,NR[
		// DO IZ=1,NZ[
		// WRITE(IOUT,360) IZ,IX,ZPLANE(IZ),ZPLANE(IZ+1),
		// (100.* SCDOSE(IZ,IX,I),SCDOSE2(IZ,IX,I)*SCDOSE(IZ,IX,I),I=1,2);
		// "100 TO GET PERCENT, AND ERROR WAS PERCENT SO NEEDED TO"
		// "DIVIDE BY 100 AND HENCE CANCELED"
		// IF (BSTAT(IZ,1) = -1 | BDSTAT(IZ)= -1) [WRITE(IOUT,315);]
		// ]
		// ]
		for (int IX = 1; IX <= EGS4Geom.NR; IX++)// DO IX=1,NR[
		{
			for (int IZ = 1; IZ <= EGS4Geom.NZ; IZ++)// DO IZ=1,NZ[
			{
				// WRITE(IOUT,310) IZ,IX,ZPLANE(IZ),ZPLANE(IZ+1),
				// SCEDPN(IZ,IX,1),(SCEDPN(IZ,IX,1)*SCEDPN2(IZ,IX,1)/100.),
				// SCDOSE(IZ,IX,3), SCDOSE2(IZ,IX,3) ;
				/*
				 * 360 FORMAT(' NZ,NR=',I3,',',I3,F8.2,' -',F8.2,2X,2(F7.3,'
				 * (+/- ',F6.3,' )')); %I4
				 */
				s = " NZ,NR= "
						+ EGS4.format(IZ, 3)
						+ ","
						+ EGS4.format(IX, 3)
						+ "  "
						+ EGS4.format(EGS4Geom.ZPLANE[IZ - 1], 8, true)
						+ " - "
						+ EGS4.format(EGS4Geom.ZPLANE[IZ], 8, true)
						+ "  "
						+
						// EGS4.format(SCEDPN[IZ-1][IX-1][2],10,true)+" (+/-"+
						// EGS4.format(SCEDPN[IZ-1][IX-1][2]*SCEDPN2[IZ-1][IX-1][2]/100.,8,true)+")  "+
						EGS4.format(100. * SCDOSE[IZ - 1][IX - 1][0], 8, true)
						+ " (+/- "
						+ EGS4.format(SCDOSE2[IZ - 1][IX - 1][0]
								* SCDOSE[IZ - 1][IX - 1][0], 7, true)
						+ "  "
						+ EGS4.format(100. * SCDOSE[IZ - 1][IX - 1][1], 8, true)
						+ " (+/- "
						+ EGS4.format(SCDOSE2[IZ - 1][IX - 1][1]
								* SCDOSE[IZ - 1][IX - 1][1], 7, true);
				EGS4.seqStr = s;
				if (EGS4.iprint > 0)
					printSequence(EGS4.seqStr);

				if (BSTAT[IZ - 1][0] == -1 || BDSTAT[IZ - 1] == -1)// (BSTAT(IZ,1)
																	// = -1)
				{
					// WRITE(IOUT,315); **** WARNING **** BAD STATISTICS!
					s = "**** WARNING **** BAD STATISTICS!";
					EGS4.seqStr = s;
					if (EGS4.iprint > 0)
						printSequence(EGS4.seqStr);

				}
			}// ]
		}// ]

		return;

	}

	/**
	 * Return base 10 logarithm.
	 * @param x x
	 * @return the result
	 */
	public static double DLOG10(double x) {
		double result = Math.log(x);
		result = result / Math.log(10);
		return result;
	}

	// CALL PLOTPH(TITLE,SCPDST,SCPDST2,SCPCUM,SCPCUM2,
	// SCPTOT,SCPTOT2,SCDFEP,SCDFEP2,MAXBIN,SLOTE,BINTOP,
	// IHSTRY,SCOMEG,SCOMEG2,SCPHEN,SCPHEN2);

	// SUBROUTINE PLOTPH(TITLE,PDST,PDSTUN,PCUM,PCUMUN,PTOT,PTOTUN,DFEP,DFEPUN,
	// MAXBIN,SLOTE,BINTOP,IHSTRY,OMEG,OMEGUN,PHEN,PHENUN);
	/**
	 * Plot some useful information (pulse height distribution).
	 */
	public void PLOTPH() {
		String s = "";
		int ll = 0;
		// "GET PEAK IN SPECTRUM"
		double sfac = 0.0;
		double EB = 0.0;
		String SPACE = " ";// /,BAR/'|'/,SYMBOL/$S'*+$-#@'/;
		String BAR = "|";
		String[] SYMBOL = { "*", "+", "$", "-", "#", "@" };
		String[] LINE = new String[61];
		int NHIST = 61;
		double HIST = 61.;
	//	int IOUT = 1;
		int ILEV = 0;

		for (int IB = 1; IB <= MAXBIN; IB++) {
			if (SCPDST[IB - 1] > sfac) {
				sfac = SCPDST[IB - 1];
			}
		}
		// WRITE(IOUT,4) IHSTRY;
		s = " " + EGS4.format(IHSTRY, 12) + " HISTORIES ANALYSED";
		EGS4.seqStr = s;// \n";
		if (EGS4.iprint > 0)
			printSequence(EGS4.seqStr);

		// WRITE(IOUT,5); "HEADER FOR PLOT"
		s = EGS4.format("", 65) + "EBIN";
		EGS4.seqStr = s;// \n";
		if (EGS4.iprint > 0)
			printSequence(EGS4.seqStr);
		s = " 3" + EGS4.format("", 19) + "2" + EGS4.format("", 19) + "1"
				+ EGS4.format("", 7) + "-LOG" + EGS4.format("", 8) + "0    TOP"
				+ EGS4.format("", 9) + "PDST              CUMULATIVE";
		EGS4.seqStr = s;// \n";
		if (EGS4.iprint > 0)
			printSequence(EGS4.seqStr);

		for (int IB = 1; IB <= MAXBIN; IB++) {
			if (SLOTE > 0.0) {
				EB = IB * SLOTE;
			} else {
				EB = BINTOP[IB - 1];
			}
			for (int J = 1; J <= NHIST; J++) {
				LINE[J - 1] = SPACE;
				// "BLANK ENTIRE LINE EACH TIME"
			}
			ILEV = 0;
			if (SCPDST[IB - 1] > 0.0) {
				Double dbl = new Double(HIST + DLOG10(SCPDST[IB - 1] / sfac)
						* 20. + 0.5);
				ILEV = dbl.intValue();
				// "RUNS FROM 1 TO 61 FOR 0.0 TO 1.0"
			}
			if ((ILEV > 0) && (ILEV <= NHIST)) {
				LINE[ILEV - 1] = SYMBOL[0];
			}
			// "ADD REFERENCE BARS"
			for (int J = 1; J <= 4; J++) {
				ILEV = 20 * J - 19;
				if (LINE[ILEV - 1].compareTo(SPACE) == 0)
					LINE[ILEV - 1] = BAR;
			}
			// WRITE(IOUT,10) LINE,EB,PDST(IB),PDSTUN(IB),PCUM(IB),PCUMUN(IB);
			// 10 FORMAT(1X,61A1,F8.4,1X,2(0PF11.4,' (',F6.3,'%)'));
			s = " ";// +EGS4.format(IHSTRY,12)+" HISTORIES ANALYSED";
			for (int i = 1; i <= 61; i++) {
				s = s + LINE[i - 1];
			}
			s = s + EGS4.format(EB, 8, true) + " "
					+ EGS4.format(SCPDST[IB - 1], 11, true) + " ("
					+ EGS4.format(SCPDST2[IB - 1], 6, true) + "%)"
					+ EGS4.format(SCPCUM[IB - 1], 11, true) + " ("
					+ EGS4.format(SCPCUM2[IB - 1], 6, true) + "%)";
			EGS4.seqStr = s;// \n";
			if (EGS4.iprint > 0)
				printSequence(EGS4.seqStr);

		}// "END OF IB LOOP"

		// 4 FORMAT(/I12,' HISTORIES ANALYSED'/);
		// 5 FORMAT(/65X,'EBIN'/' 3',19X,'2',19X,'1',7X,
		// '-LOG',8X,'0 TOP',9X,'PDST CUMULATIVE');
		// WRITE(IOUT,30);
		// 30 FORMAT(/' PLOT NORMALIZED TO PEAK OF ONE, ',
		// 'PDST IS NORMALIZED TO UNIT AREA'/);
		s = " PLOT NORMALIZED TO PEAK OF ONE, "
				+ "PDST IS NORMALIZED TO UNIT AREA";
		EGS4.seqStr = s;// \n";
		if (EGS4.iprint > 0)
			printSequence(EGS4.seqStr);

		double SUM = 0.0;
		if (EGS4SrcEns.iqin == 0
				&& EGS4SrcEns.monoindex == EGS4SrcEns.iMONOENERGETIC) {
			for (int IPK = 1; IPK <= 4; IPK++) {
				SUM = SUM + SCDFEP[IPK - 1];
			}
			if (SUM >= 0.00005) {
				// "THERE ARE SOME COUNTS IN THE PEAKS"
				// WRITE(IOUT,40) (DFEP(IPK),DFEPUN(IPK),IPK=1,4);
				// 40 FORMAT(' PEAK EFFICIENCIES PER COUNT IN SPECTRUM'//
				// ' FULL ENERGY PEAK',T30,F10.4,'(+-',F6.3,'%)'/
				// ' SINGLE ESCAPE PK',T30,F10.4,'(+-',F6.3,'%)'/
				// ' DOUBLE ESCAPE PK',T30,F10.4,'(+-',F6.3,'%)'/
				// ' 511 KEV PK' ,T30,F10.4,'(+-',F6.3,'%)'/);
				s = " PEAK EFFICIENCIES PER COUNT IN SPECTRUM";
				EGS4.seqStr = s;// \n";
				if (EGS4.iprint > 0)
					printSequence(EGS4.seqStr);
				s = " FULL ENERGY PEAK";
				ll = s.length();
				ll = 30 - ll;
				s = s + EGS4.format("", ll);
				EGS4.seqStr = s + EGS4.format(SCDFEP[0], 10, true) + "(+-"
						+ EGS4.format(SCDFEP2[0], 6, true) + "%)";
				if (EGS4.iprint > 0)
					printSequence(EGS4.seqStr);
				s = " SINGLE ESCAPE PK";
				ll = s.length();
				ll = 30 - ll;
				s = s + EGS4.format("", ll);
				EGS4.seqStr = s + EGS4.format(SCDFEP[1], 10, true) + "(+-"
						+ EGS4.format(SCDFEP2[1], 6, true) + "%)";
				if (EGS4.iprint > 0)
					printSequence(EGS4.seqStr);
				s = " DOUBLE ESCAPE PK";
				ll = s.length();
				ll = 30 - ll;
				s = s + EGS4.format("", ll);
				EGS4.seqStr = s + EGS4.format(SCDFEP[2], 10, true) + "(+-"
						+ EGS4.format(SCDFEP2[2], 6, true) + "%)";
				if (EGS4.iprint > 0)
					printSequence(EGS4.seqStr);
				s = " 511 KEV PK";
				ll = s.length();
				ll = 30 - ll;
				s = s + EGS4.format("", ll);
				EGS4.seqStr = s + EGS4.format(SCDFEP[3], 10, true) + "(+-"
						+ EGS4.format(SCDFEP2[3], 6, true) + "%)";
				if (EGS4.iprint > 0)
					printSequence(EGS4.seqStr);

			}
		} else {
			if (EGS4SrcEns.monoindex == EGS4SrcEns.iMONOENERGETIC) {
				s = " ELECTRON KINETIC THRESHOLD ENERGY="
						+ EGS4.format(ETHRESHOLD, 10, true) + " MeV";
				EGS4.seqStr = s;// \n";
				if (EGS4.iprint > 0)
					printSequence(EGS4.seqStr);
			}
		}
		double UNCERT = 0.0;
		if (EGS4SrcEns.SCOMEG == 0.0) {
			// "THIS IS A PARALLEL EXTERNAL BEAM OR INTERNAL SOURCE CASE"
			// WRITE(IOUT,50) PTOT,PTOTUN;
			// 50 FORMAT(/' FRACTION OF INITIAL PARTICLES (HITTING DETECTOR
			// HOUSING ',
			// '(if parallel beam) WHICH CAUSE PULSE = ',0PF14.4,'(',F6.3,'%)');
			// MODIF weight=realsolidangle /4 pi, so
			if (SOURCE == -1) {
				s = "  FRACTION OF INITIAL PARTICLES (HITTING DETECTOR HOUSING "
						+ "(if parallel beam) WHICH CAUSE PULSE = "
						+ EGS4.format(SCPTOT * 4 * Math.PI, 14, true)
						+ "("
						+ EGS4.format(SCPTOT2, 6, true) + "%)";
				// ll=s.length();ll=30-ll;s=s+EGS4.format("",ll);
				EGS4.seqStr = s;
				if (EGS4.iprint > 0)
					printSequence(EGS4.seqStr);
			}
		} else {
			// "EXTERNAL POINT SOURCE GEOMETRY"
			UNCERT = Math.sqrt(SCPTOT2 * SCPTOT2 + EGS4SrcEns.SCOMEG2
					* EGS4SrcEns.SCOMEG2); // "UNCERTAINTY ON RATIO"
			// WRITE(IOUT,60)
			// OMEG,OMEGUN,PTOT/(4.*3.141593),PTOTUN,PTOT/OMEG,UNCERT;
			// 60 FORMAT(/' SOLID ANGLE SUBTENDED BY DETECTOR HOUSING = ',
			// 1PE13.3,'(',0PF5.1,'%)'/
			// /' FRACTION OF PARTICLES INTO 4-PI WHICH CAUSE A PULSE = ',
			// 1PE13.3,'(',0PF5.1,'%)'/
			// /' FRACTION OF PARTICLES INCIDENT ON HOUSING WHICH CAUSE A PULSE
			// =',
			// 1PE13.3,'(',0PF5.1,'%)');
			if (SOURCE == -1) {
				s = " SOLID ANGLE SUBTENDED BY DETECTOR HOUSING =~                   "
						+ EGS4.format(EGS4SrcEns.SCOMEG, 13, true)
						+ "("
						+ EGS4.format(EGS4SrcEns.SCOMEG2, 5, true) + "%)";
				// ll=s.length();ll=30-ll;s=s+EGS4.format("",ll);
				EGS4.seqStr = s;
				if (EGS4.iprint > 0)
					printSequence(EGS4.seqStr);
			}
			s = " FRACTION OF PARTICLES INTO 4-PI WHICH CAUSE A PULSE =          "
					+ EGS4.format(SCPTOT, 13, true)
					+ "("
					+ EGS4.format(SCPTOT2, 5, true) + "%)";
			// ll=s.length();ll=30-ll;s=s+EGS4.format("",ll);
			EGS4.seqStr = s;
			if (EGS4.iprint > 0)
				printSequence(EGS4.seqStr);
			if (SOURCE == -1) {
				s = " FRACTION OF PARTICLES INCIDENT ON HOUSING WHICH CAUSE A PULSE ="
						+ EGS4.format(SCPTOT * 4 * Math.PI / EGS4SrcEns.SCOMEG,
								13, true)
						+ "("
						+ EGS4.format(UNCERT, 5, true)
						+ "%)";
				// ll=s.length();ll=30-ll;s=s+EGS4.format("",ll);
				EGS4.seqStr = s;
				if (EGS4.iprint > 0)
					printSequence(EGS4.seqStr);
			}
		}
		// WRITE(IOUT,70) PHEN,PHENUN;//=>SCPHEN,SCPHEN2
		// 70 FORMAT(/' ENERGY DEPOSITED IN DETECTOR PER INITIAL PARTICLE:',
		// 0PF14.5,' MeV (+-',F6.3,'%)');
		s = "  ENERGY DEPOSITED IN DETECTOR PER INITIAL PARTICLE:"
				+ EGS4.format(SCPHEN, 14, true) + "(+-"
				+ EGS4.format(SCPHEN2, 6, true) + "%)";
		// ll=s.length();ll=30-ll;s=s+EGS4.format("",ll);
		EGS4.seqStr = s;
		if (EGS4.iprint > 0)
			printSequence(EGS4.seqStr);

		s = "*******************************************************************************************";
		String PHOTONs = "PHOTON";
		if (EGS4SrcEns.monoindex == EGS4SrcEns.iMONOENERGETIC) {
			if (EGS4SrcEns.iqin == -1) {
				PHOTONs = "ELECTRON";
			} else if (EGS4SrcEns.iqin == 1) {
				PHOTONs = "POSITRON";
			}
			EGS4.seqStr = s;
			if (EGS4.iprint > 0)
				printSequence(EGS4.seqStr);
			if (EGS4SrcEns.iqin == 0) {
				s = " REAL DETECTOR EFFICIENCY (PEAK EFF.) FOR INCIDENT "
						+ PHOTONs + " ENERGY (%)= "
						+ EGS4.format(efficiency, 14, true) + "(+-"
						+ EGS4.format(efficiency_error, 6, true) + "%)";
				// ll=s.length();ll=30-ll;s=s+EGS4.format("",ll);
				EGS4.seqStr = s;
				if (EGS4.iprint > 0)
					printSequence(EGS4.seqStr);

				s = " DETECTOR TOTAL EFFICIENCY (WHOLE SPECTRUM) FOR INCIDENT "
						+ PHOTONs + " ENERGY (%)= "
						+ EGS4.format(photontotal_efficiency, 14, true) + "(+-"
						+ EGS4.format(photontotal_efficiency_error, 6, true)
						+ "%)";
				// ll=s.length();ll=30-ll;s=s+EGS4.format("",ll);
				EGS4.seqStr = s;
				if (EGS4.iprint > 0)
					printSequence(EGS4.seqStr);

			}
			// *************************
			else {
				s = " REAL DETECTOR EFFICIENCY FOR INCIDENT " + PHOTONs
						+ " ENERGY (%)= " + EGS4.format(eefficiency, 14, true)
						+ "(+-" + EGS4.format(eefficiency_error, 6, true)
						+ "%)";
				// ll=s.length();ll=30-ll;s=s+EGS4.format("",ll);
				EGS4.seqStr = s;
				if (EGS4.iprint > 0)
					printSequence(EGS4.seqStr);
			}
			// *************************
			s = "******************************************************************************************";
			EGS4.seqStr = s;
			if (EGS4.iprint > 0)
				printSequence(EGS4.seqStr);
		}
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
		// EGS4SrcEns.IFPB=1;
	//	double D = 0.0;
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

			if (scattCaseB)
				diam = EGS4Geom.RCYL[1] * EGS4SrcEns.$ONE_EPS;// inner!!

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
			// ---------------------------------
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
					break;// gresit R2 always 0=>ixin in old case=1!!
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
			double z1 = hsource * r;// -hsourceup;//0=>-hup;1=h-hup!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			boolean infcyl = false;// test if inf cylinder exists
			if (hsource > hsourceup + hdettot)// +EGS4Geom.ZPLANE[EGS4Geom.NPLANE-1])//hdet)
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
			// x and y does not change with angle.Anyway, in angle case, circle
			// projection is not a circle!!
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
		double z1 = (hsource - ethick) * r;// in source:0->hs-ethick
		EGS4SrcEns.DISTZ = hsource - z1;// +plasticul@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
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

	@SuppressWarnings("unused")
	/**
	 * Not used here.
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
	///	double uz = w;
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
		double d = Math.abs(z1 - (hsourceup + hdettot));// EGS4Geom.ZPLANE[EGS4Geom.NPLANE-1]);//yes
														// ca-i ABS!!!
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

	public static boolean NEGVAL = false;
	public static int CHECK = 0;

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

	// REPLACE {$CONVERT_INT(#)_TO_CHAR(#);} WITH {
	// ;"This MACRO changes integers of three digits"
	// "or less into a character string."
	// "In each case, 48 must be added to the digit we are"
	// "converting into ASCII since '0' is ASCII # 48."
	// "Aaron Merovitz, Jan 1998"
	// "{P1} is the integer value < 1000"
	// "{P2} is the character string returned"
	/**
	 * Changes integers of three digits or less into a character string.
	 * @param P1 P1
	 * @return the result
	 */
	public static String CONVERT_INT(int P1) {
		String P2 = "";
		NEGVAL = false;
		if (P1 < 0) {
			NEGVAL = true;
			P1 = Math.abs(P1);
		}
		if (P1 == P1 % 1000)// if ((P1 = MOD({P1},1000)))
		{
			CHECK = 1;
			if (P1 == P1 % 10)// if({P1} = MOD({P1},10))
			{
				if (NEGVAL) {
					// P2="-"//CHAR({P1}+48);// // = Concatenation
					P2 = "-" + asciiToStr(P1 + 48);//
				} else {
					P2 = asciiToStr(P1 + 48);// {P2}=CHAR({P1}+48);
				}
				CHECK = 0;
			}
			if ((CHECK != 0) && (P1 == P1 % 100))// ({P1} = MOD({P1},100)))
			{
				String a1 = asciiToStr(P1 / 10 + 48);// CHAR(({P1}/10)+48);
				int integ = (P1 % 10) + 48;// (MOD({P1},10)+48);
				String a2 = asciiToStr(integ);// a(2)=CHAR(int);
				if (NEGVAL) {
					P2 = "-" + a1 + a2;// {P2}='-'//a(1)//a(2);
				} else {
					P2 = a1 + a2;// {P2}=a(1)//a(2);
				}
				CHECK = 0;
			}
			if (CHECK != 0) {
				String a1 = asciiToStr(P1 / 100 + 48);
				;// a(1)=CHAR(({P1}/100)+48);
				int integ = (P1 % 100) / 10 + 48;// (MOD({P1},100)/10+48);
				String a2 = asciiToStr(integ);// a(2)=CHAR(int);
				integ = (P1 % 10) + 48;// MOD({P1},10)+48;
				String a3 = asciiToStr(integ);// CHAR(int);
				if (NEGVAL) {
					P2 = "-" + a1 + a2 + a3;// {P2}='-'//a(1)//a(2)//a(3);
				} else {
					// {P2}=a(1)//a(2)//a(3);
					P2 = a1 + a2 + a3;
				}
			}
		} else {
			P2 = "3max";// {P2}='3max';
		}

		return P2;
	}

	/*
	 * ;"THE FOLLOWING MACROS MAKE THE STATISTICAL ANALYSIS EASY. THEY ARE
	 * "CODED IN GENERAL FASHION AND MAY BE USEFUL FOR OTHER USER CODES.
	 * 
	 * REPLACE {$SCORE(#,#:#)} WITH {;
	 * 
	 * "Scoring macro used in AUSGAB for SCSTP,SCFLU,SCFLEP"
	 * "{P1}{P2}=scoring array (eg SCFLU(IB,IZ,IX,IQL,IT))"
	 * "{P3}=quantity to be scored (eg BINFR(IB)*WTTVST)"
	 * 
	 * "If the (primary) history number, NHSTRY, is the same as the history"
	 * "that last scored the quantity of interest in this voxel, {P1}_LAST{P2},"
	 * "then {P3} is added to a temporary array, {P1}_TMP{P2}.  Otherwise, we add"
	 * "{P1}_TMP{P2} to {P1}{P2}, {P1}_TMP{P2}*{P1}_TMP{P2} to {P1}2{P2},"
	 * "set {P1}_TMP{P2}={P3}, and set {P1}_LAST{P2}=NHSTRY."
	 * 
	 * IF(NHSTRY={P1}_LAST{P2})[ {P1}_TMP{P2}={P1}_TMP{P2} + {P3}; ] ELSE[
	 * {P1}{P2}={P1}{P2}+{P1}_TMP{P2}; {P1}2{P2}={P1}2{P2} +
	 * {P1}_TMP{P2}*{P1}_TMP{P2}; {P1}_TMP{P2}={P3}; {P1}_LAST{P2}=NHSTRY; ] ; }
	 * //$SCORE(SCSTP, :1); "counting number of charged particle steps"
	 * if(EGS4SrcEns.NHSTRY==SCSTP_LAST) { SCSTP_TMP=SCSTP_TMP + 1; } else {
	 * SCSTP=SCSTP+SCSTP_TMP; SCSTP2=SCSTP2 + SCSTP_TMP*SCSTP_TMP; SCSTP_TMP=1;
	 * SCSTP_LAST=EGS4SrcEns.NHSTRY; }
	 * 
	 * REPLACE {$SCORET(#,#,#:#:#)} WITH {;
	 * 
	 * "Scoring macro used in AUSGAB for SCTFLU,SCEAVE,SCTLEP and SCELEP"
	 * "{P1}{P3} and {P2}{P3}= scoring arrays of interest (eg SCTFLU(IZ,IX,IQL,IT)"
	 * "                       and SCEAVE(IZ,IX,IQL,IT)"
	 * "{P4}=quantity to be scored in {P1}{P3} (ie WTTVST)"
	 * "{P5}=quantity to be scored in {P2}{P3} (ie WTEKIN)"
	 * 
	 * "If the (primary) history number, NHSTRY, is the same as the history"
	 * "that last scored the quantities of interest in this voxel, {P1}_LAST{P3},"
	 * "then {P4} is added to {P1}_TMP{P3} and {P5} is added to {P2}_TMP{P3}."
	 * "Otherwise, we add {P1}_TMP{P3} to {P1}{P3}, {P1}_TMP{P3}^2 to"
	 * "{P1}2{P3}, add {P2}_TMP{P3} to {P2}{P3}, {P2}_TMP{P3}^2"
	 * "to {P2}2{P3}, and {P1}_TMP{P3}*{P2}_TMP{P3} to"
	 * "{P2}_COV{P3}.  Then, we set {P1}_TMP{P3}={P4}, {P2}_TMP{P3}={P5}"
	 * "and {P1}_LAST{P3}=NHSTRY."
	 * 
	 * IF(NHSTRY={P1}_LAST{P3})[ {P1}_TMP{P3}={P1}_TMP{P3} + {P4};
	 * {P2}_TMP{P3}={P2}_TMP{P3} + {P5}; ] ELSE[ {P1}{P3}={P1}{P3}+{P1}_TMP{P3};
	 * {P1}2{P3}={P1}2{P3} + {P1}_TMP{P3}*{P1}_TMP{P3};
	 * {P2}{P3}={P2}{P3}+{P2}_TMP{P3}; {P2}2{P3}={P2}2{P3} +
	 * {P2}_TMP{P3}*{P2}_TMP{P3};
	 * {P2}_COV{P3}={P2}_COV{P3}+{P1}_TMP{P3}*{P2}_TMP{P3}; {P1}_TMP{P3}={P4};
	 * {P2}_TMP{P3}={P5}; {P1}_LAST{P3}=NHSTRY; ] ; }
	 * 
	 * REPLACE {$ANALYZE(#,#:#)} WITH {;
	 * 
	 * "Macro to analyze uncertainty:"
	 * "{P1}{P2}=scoring array (eg SCFLU(IB,IZ,IX,IQL,IT))"
	 * "{P3}=quantity to normalize by (eg incident no. of particles)"
	 * 
	 * "Calculates the uncertainty on {P1}{P2}/{P3}.  The "
	 * "uncertainty is stored in {P1}2{P2} and is expressed as a percentage of"
	 * "{P1}{P2}/{P3} (max 99.9%).  Note that you must define the REAL*8 variable"
	 * "SCORE_TEMP in any subroutine where this macro is used.  This macro"
	 * "is only used in the analysis of no. of steps."
	 * 
	 * SCORE_TEMP={P1}{P2}/{P3}; {P1}2{P2}={P1}2{P2}/{P3};
	 * {P1}2{P2}=({P1}2{P2}-SCORE_TEMP*SCORE_TEMP)/({P3}-1); IF({P1}2{P2}>=0.)
	 * {P1}2{P2}= SQRT({P1}2{P2}); IF(SCORE_TEMP~=0.)[
	 * {P1}2{P2}=100*{P1}2{P2}/SCORE_TEMP; IF( {P1}2{P2} > 99.9 ) [ {P1}2{P2} =
	 * 99.9; ] ] ELSE[ {P1}2{P2}=99.9; ] ; }
	 */
}
