package danfulea.phys.egsOutput;

import danfulea.phys.Phantom;
import danfulea.phys.XRaySpectrum;
import danfulea.phys.egs.EgsQuestion;
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

/**
 * Class for computing dose involved in mammography where the breast can be 
 * considered a cylinder due to the compression device. Thus,RZ geometry is suitable for 
 * Monte-Carlo simulation using EGSnrc toolkit. <p>
 * For general human body, we must know all 
 * TPERP to all boundaries and exact boundary crossing distance USTEP which is annoying 
 * and very time-consuming. Therefore for general body (i.e. radiography, CT etc) one must 
 * consider different MC engines. DoseSimCore from phys package works well but runs in KERMA 
 * approximation which is not good enough for very accurate results especially at organs interfaces. 
 * Geant4 simulation toolkit is a very good alternative and maybe there are other. For me, GEANT4 is 
 * the answer for all MC based simulation. Here are some reasons: It is reliable (CERN); 
 * it is constantly updated when new "physics" is born; It integrates geometry-core engine in a very 
 * simple and convenient way so one can simulate anything, human body, space ships, with relatively ease.    
 * 
 * @author Dan Fulea, 05 DEC. 2005
 */

public class DoseMamo implements EgsQuestion {
	private static final String BASE_RESOURCE_CLASS = "danfulea.phys.resources.DoseSimCoreResources";
	private static ResourceBundle resources = ResourceBundle
			.getBundle(BASE_RESOURCE_CLASS);
	public static Phantom phantom;// =new Phantom(Phantom.MAMO_INDEX);
	public static XRaySpectrum xrs;
	public static double kvp = 30.;
	public static double uanod = 17.;
	public static double filtration = 0.5;
	public static double airExposure = 7.52;// mGy
	public static double phantomThickness = 5.0;
	public static double xSkinSize = 7.0;

	private static double[] ENSRCD;
	private static double[] SRCPDF;
	public static int $NENSRC = 300;
	public static double[] srcpdf_at = new double[$NENSRC];// ($NENSRC),
	public static int[] srcbin_at = new int[$NENSRC];// ($NENSRC)

	public static int imin = 0;
	public static double breastThickness = 1.0;// 4.2;
	public static double breastRadius = 7.0;
	public static double fsd = 55.0;
	public static double effdose = 0.0;
	public static double effdose2 = 0.0;
	public static double[] dosem = new double[1];// a single organ!!
	public static double[] dosem2 = new double[1];// a single organ!!
	public static double[] density = new double[1];// a single organ!!
	public static double[] volume = new double[1];// a single organ!!
	public static double[] mass = new double[1];// a single organ!!
	public static double[] wt = new double[1];
	public static double[] dose = new double[1];// a single organ!!
	public static double[] dose_tmp = new double[1];// a single organ!!
	public static double[] dose2 = new double[1];// a single organ!!
	public static double[] kerma = new double[1];// a single organ!!
	public static double[] kerma_tmp = new double[1];// a single organ!!
	public static double[] kerma2 = new double[1];// a single organ!!
	public static int[] dose_last = new int[1];
	public static boolean indet = false;
	public static int IHSTRY = 0;
	public static int NHSTRY = 0;
	public static int nreg = 2;// number of regions
	public static double ein = 0.0;
	public static double E_inc = 10.0;// total energy in MeV
	public static int iqin = 0;
	public static double xin = 0.0;
	public static double yin = 0.0;
	public static double zin = 0.0;
	public static double uin = 0.0;
	public static double vin = 0.0;
	public static double win = 0.0;
	public static int irin = 2;
	public static double WEIGHT = 1.0;
	public static double AINFLU = 0.0;
	public static int MONOEN = 0;
	// LAST PRIMARY HISTORY TO DEPOSIT KERMA
	public static int[] kerma_last = new int[1];
	public static int[] kerma_lastold = new int[1];
	public static double[] kerma_tmpold = new double[1];// a single organ!!
	public static double[] dosetokerma2 = new double[1];

	public static JTextArea jta;
	public static boolean systemOutB = true;

	public static int $NBATCH = 10;// "OUTPUT BATCHES                             "
	public static int JCASE = 0;; // "no. of histories per batch"
	public static int $NCASEMIN = 100;// "min. no. of histories                        "

	public static int $MAXBRSPLIT = 200;// "MAX BREM SPLITTING NUMBER"

	public static int NCASE = 0;
	public static int NCASEO = 0;
	public static int NCASET = 0;

	public static double STATLM = 0.0;// EIN,

	// -------------INPUTS-------------------------------------------
	// @ ELECTRON TRANSPORT
	public static final int ETRANS_NORMAL = 0;
	public static final int ETRANS_NO_INTERACTION = 1;
	// @ NUMBER OF HISTORIES
	public static final int NCASE_MIN = 1;
	public static final int NCASE_MAX = 461168600;// 4.611686e18;
	public static final int NCASE_DEFAULT = 20000;
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

	public static String photon_xsection = "";// default

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

	public static int IBTCH = 0;
	public static double SCORE_NORM_NUM = 0.0;
	public static double SCORE_TEMP = 0.0;

	public static int IRL = 0;
	public static int MEDNUM = 0;
	// ================================dose
	public static double TDSMAX = 0.0;
	public static int IDSMAX = 0;
	public static double TDOS = 0.0;
	public static double TDOS2 = 0.0;

	public static int nENHREG = 0;
	public static int[] NENHLO = new int[1000];
	public static int[] NENHHI = new int[1000];
	public static int ics_enhance = 0;
	public static double efficiency = 0.0;
	public static double efficiency_error = 0.0;
	// ========================
	public static int IPRINT = 2;
	public static boolean createOutputFile = false;
	public static boolean putInFile = false;// internal var defining when and
											// what to print
	private String filename = "";
	FileWriter sigfos;

	// -------------end INPUTS-------------------------------------------
	/**
	 * Constructor.
	 */
	public DoseMamo() {
		createOutputFile = false;// allways

		putInFile = true;
		Calendar cal = Calendar.getInstance();
		String fs = cal.get(Calendar.YEAR) + "_" + cal.get(Calendar.MONTH) + "_"
				+ cal.get(Calendar.DAY_OF_MONTH) + "_" + cal.get(Calendar.HOUR) + "_"
				+ cal.get(Calendar.MINUTE) + "_" + cal.get(Calendar.SECOND) + ".txt";
		filename = fs;// will be 2005_10_25_14_30_56.txt
		// ====================================================================================

		// phantomThickness=5.0;
		// xSkinSize=7.0;
		// kvp=30.;
		// uanod=17.;
		// filtration=0.5;
		// airExposure=7.52;//mGy
		// fsd=55.0;

		phantom = new Phantom(Phantom.MAMO_INDEX);

		phantom.setOrganThickness(phantomThickness);// in AEIA 4-6cm
		phantom.setMaximumOrganDimension(xSkinSize);

		// phantom.setBounds();

		breastThickness = phantom.getOrganThickness();
		breastRadius = phantom.getMaximumOrganDimension();

		density[0] = ((Double) resources.getObject("breast.diagnostic.density"))
				.doubleValue();
		volume[0] = Math.PI * breastRadius * breastRadius * breastThickness;// cm3
		mass[0] = density[0] * volume[0];// grams
		wt[0] = ((Double) resources.getObject("breast.wt")).doubleValue();// 0.05;
		// ====================================================================================
		init();
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

		EGS4.setMXMED(10);// allways//"MAX # OF MEDIA 		"
		EGS4.setMXREG(200);// allways
		EGS4.setMXSTACK(4000);// allways//"NEED HUGE STACK FOR CORRELATIONS+splitting   "
		EGS4.setMXRANGE(500);// allways
								// //"for range arrays used in range_rejection()"
		EGS4.$MXMDSH = 100;// allways
		// flush
		EGS4.reset();// allways
		EGS4.startSimulationTime = System.currentTimeMillis();

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
		// --interface->LINK
		EGS4.eq = this;// pass the printing mode
		EGS4Core.eq = this;// pass the printing mode
		EGS4Macro.eq = this;// pass the printing mode
		EGS4SrcEns.eq = this;// pass the printing mode
		EGS4Geom.eq = this;// pass the printing mode
		EGS4Grid.eq = this;// pass the printing mode
		// --Macro and defaults param:

		EGS4.egs_set_defaults();// first default then://allways
		EGS4.RandomUse = 1;// use EGS ranlux or ranmar generators//allways
		EGS4.ranluxB = true;// use ranmar!!//allways

		EGS4.ispmfp = EGS4.iDose;// select photon mfp
		EGS4.iurd = 0;// no user range discard
		EGS4.iraycorr = EGS4.iDose;// rayleigh correction
		EGS4.isemfp = EGS4.iDose;// $SELECT-ELECTRON-MFP
		// EGS4.iGeom=1;//EGS4.iCavity;//1=RZ GEOM
		EGS4Macro.ismfpfpb = EGS4.iDose;// select mfp parallel beam=NOT USED see
										// gammadetEff
		// EGS4Macro.irange_rej=EGS4.iCavity;//not used range rejection
		EGS4.USER_CONTROLS_TSTEP_RECURSION = 0;// no effect here
		EGS4.hatchindex = 0;// no effect here
		// -----------------
		// iqin=0;//incident particle allways photon
		// MONOEN=1;//allways
		// E_inc=10.0;//total in MeV!!!!!!!
		IPRINT = 1;// allways
		EGS4.iprint = IPRINT;// EGS4.iprint=2;//summary
		// inputs
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
		EGS4.seqStr = " **DOSE: MAMMOGRAPHY**";
		if (EGS4.iprint > 0)
			printSequence(EGS4.seqStr);

		inputs();
		if (EGS4.STOPPROGRAM) {
			closeFile();
			return;
		}

		ENSRC();

		NHSTRY = 0;

		if (NCASE / $NBATCH == 0) {
			NCASE = $NBATCH;
		}
		JCASE = NCASE / $NBATCH;
		NCASE = JCASE * $NBATCH;// "NUMBER OF HISTORIES PER BATCH

		IHSTRY = NCASEO; // "reset the number of histories counter"

		NHSTRY = 0;

		// "set up ausgab calls"
		for (int J = 1; J <= 5; J++) {
			EGS4.iausfl[J - 1] = 1;// IAUSFL(J)=1;
		}
		for (int J = 6; J <= 25; J++) {
			EGS4.iausfl[J - 1] = 0;
		} // "NORMAL EXECUTION"

		EGS4.iausfl[5] = 1; // "AFTER TRANSPORT"
		EGS4.iausfl[16] = 1; // "after pair
		EGS4.iausfl[18] = 1; // "After COMPTON"
		EGS4.iausfl[20] = 1; // "After Photo"

		if (EGS4Macro.cs_enhance > 1.0001) {
			// write(6,*) 'flagged all photon intereaction types';
			EGS4.iausfl[15] = 1; // "Before pair"
			EGS4.iausfl[17] = 1; // "Before Compton"
			EGS4.iausfl[18] = 1; // "After Compton"
			EGS4.iausfl[19] = 1; // "Before photoelectric"
			EGS4.iausfl[20] = 1; // "After photoelectric"
			EGS4.iausfl[23] = 1; // "Before Rayleigh"
			EGS4.iausfl[24] = 1; // "After Rayleigh"
		} else {

			for (int j = 1; j <= EGS4.$MXREG; j++) {
				EGS4Macro.iefl[j - 1] = 0;
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

		if (MONOEN == 0) {// "MONOENERGETIC INPUT BEAM"
			if (iqin == 0) {
				EI = ein;
			} else {
				EI = ein + EGS4.RM;
			}
			EKMAX = ein; // "MAXIMUM KINETIC ENERGY"
		} else if (MONOEN == 1) {// "ENERGY SPECTRUM"
			ENSRC1();// "NORMALIZE THE ENERGY DISTRIBUTION"
			EKMAX = ENSRCD[ENSRCD.length - 1] / 1000.;// "MAXIMUM KINETIC ENERGY IN THE SPECTRUM"
		}

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

		// "******************************************************************************
		// "
		// " *** SECTION 2 ***
		// "
		// "------------------------------------------------------------------------------
		// "Output batches. Statistical analysis is done after each batch.
		// Execution
		// "stops if the desired statistical accuracy is obtained or there is
		// not enough
		// "time to do another batch.
		boolean ENDSIM = false;
		for (int IBATCH = 1; IBATCH <= $NBATCH; IBATCH++) {
			ENDSIM = false;
			long startTime = System.currentTimeMillis();
			IBTCH = IBATCH;
			if (IBATCH == 1) {
				EGS4.seqStr = " BATCH" + EGS4.format("", 8) + "ELAPSED"
						+ EGS4.format("", 3) + "time of day"
						+ EGS4.format("", 2) + " dose stats(%)";
				if (EGS4.iprint > 0)
					printSequence(EGS4.seqStr);
			} else {// " not first batch"
			}// " end of before batch ne 1 block"

			for (int ICASE = 1; ICASE <= JCASE; ICASE++) {// "now fill each IS bin"

				IHSTRY = IHSTRY + 1;

				EGS4Macro.NFTIME = 0;

				// "calculate the initial energy if a distribution is to be used"
				if (MONOEN != 0)// allways photons with kev
				{// "if equal to 0, it is monoenergetic"
					ein = ENSRCH() / 1000.;
					if (iqin == 0) {
						EI = ein;
					} else {
						EI = ein + EGS4.RM;
					}// "total energy"
				}

				// "FOR AN INPUT ENERGY SPECTRUM, DETAILED FORCING MACRO IS USED"

				EGS4.LATCHI = 0;
				// "SET INITIAL DOSE COMPONENTS"
				// ****************************FIX WEIGHT
				// *****************************************************
				fixEmAll();
				// ********************************************************************************************
				if (indet)
					SHOWER();
				if (EGS4.STOPPROGRAM) {
					closeFile();
					return;
				}

			}// "END OF THE ICASE LOOP"

			// "DO STATISTICAL ANALYSIS ON THE PEAK DOSE REGION TO SEE IF WE QUIT EARLY"
			TDSMAX = 0.0;
			for (IRL = 2; IRL <= nreg; IRL++) {
				if ((dose[IRL - 2] + dose_tmp[IRL - 2]) / mass[IRL - 2] > TDSMAX) {
					TDSMAX = (dose[IRL - 2] + dose_tmp[IRL - 2])
							/ mass[IRL - 2];
					IDSMAX = IRL;
				}
			}

			// "NOW DO STATS ON THE PEAK REGION"
			if (TDSMAX > 0.0) {
				TDOS = dose[IDSMAX - 2] + dose_tmp[IDSMAX - 2];
				TDOS2 = dose2[IDSMAX - 2] + dose_tmp[IDSMAX - 2]
						* dose_tmp[IDSMAX - 2];
				SCORE_NORM_NUM = dble(IHSTRY);

				if (SCORE_NORM_NUM > 1.) {
					TDOS = TDOS / SCORE_NORM_NUM;
					TDOS2 = TDOS2 / SCORE_NORM_NUM;
					TDOS2 = (TDOS2 - TDOS * TDOS) / (SCORE_NORM_NUM - 1);
					if (TDOS2 >= 0.)
						TDOS2 = Math.sqrt(TDOS2);
					TDOS2 = Math.min(TDOS2 / TDOS * 100., 99.9);
					if ((TDOS2 <= STATLM) && (STATLM != 0.0)) {
						// "REACHED OBJECTIVE - PRINT MESSAGE AND JUMP OUT OF SIMULATION LOOP"
						EGS4.seqStr = " DESIRED STATISTICAL ACCURACY OBTAINED.";
						if (EGS4.iprint > 1)
							printSequence(EGS4.seqStr);
						EGS4.seqStr = " STATS IN PEAK DOSE REGION (REGION "
								+ EGS4.format(IDSMAX, 3)
								+ EGS4.format(TDOS2, 6, true) + "%" + " AFTER "
								+ EGS4.format(IBTCH, 2) + " BATCHES";
						if (EGS4.iprint > 1)
							printSequence(EGS4.seqStr);

						AINFLU = AINFLU * dble(IHSTRY) / dble(NCASET);
						// GO TO :END-SIM:;
						ENDSIM = true;
						break;
					}
				}
			}

			String timePerBatch = EGS4.timeElapsedShort(startTime);
			Calendar call = Calendar.getInstance();
			String timeday = call.get(Calendar.HOUR) + ":" + call.get(Calendar.MINUTE)
					+ ":" + call.get(Calendar.SECOND);

			EGS4.seqStr = EGS4.format("", 2) + EGS4.format(IBATCH, 3)
					+ EGS4.format("", 2) + EGS4.format(timePerBatch, 14)
					+ EGS4.format("", 3) + EGS4.format(timeday, 9)
					+ EGS4.format("", 6) + EGS4.format(TDOS2, 7, true);
			if (EGS4.iprint > 0)
				printSequence(EGS4.seqStr);

		}// "END OF SIMULATIONS"

		if (!ENDSIM) {
			EGS4.seqStr = " DESIRED STATISTICAL ACCURACY OF "
					+ EGS4.format(STATLM, 6, true) + "%" + " NOT REACHED!";
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);
			EGS4.seqStr = " STATS IN PEAK DOSE REGION= "
					+ EGS4.format(TDOS2, 6, true) + " %" + " AFTER "
					+ EGS4.format(IBTCH, 2) + " BATCHES";
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);

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
		SCORE_NORM_NUM = EGS4SrcEns.dble(IHSTRY);

		for (IRL = 2; IRL <= nreg; IRL++) {
			dose[IRL - 2] = dose[IRL - 2] + dose_tmp[IRL - 2];
			dose2[IRL - 2] = dose2[IRL - 2] + dose_tmp[IRL - 2]
					* dose_tmp[IRL - 2];

			kerma[IRL - 2] = kerma[IRL - 2] + kerma_tmp[IRL - 2];
			kerma2[IRL - 2] = kerma2[IRL - 2] + kerma_tmp[IRL - 2]
					* kerma_tmp[IRL - 2];

			if (dose_last[IRL - 2] == kerma_last[IRL - 2]) {
				dosetokerma2[IRL - 2] = dosetokerma2[IRL - 2]
						+ dose_tmp[IRL - 2] * kerma_tmp[IRL - 2];
			}
		}

		for (IRL = 2; IRL <= nreg; IRL++) {
			SCORE_TEMP = dose[IRL - 2] / SCORE_NORM_NUM;
			dose2[IRL - 2] = dose2[IRL - 2] / SCORE_NORM_NUM;
			dose2[IRL - 2] = (dose2[IRL - 2] - SCORE_TEMP * SCORE_TEMP)
					/ (SCORE_NORM_NUM - 1);
			if (dose2[IRL - 2] >= 0.)
				dose2[IRL - 2] = Math.sqrt(dose2[IRL - 2]);
			if (SCORE_TEMP != 0.) {
				dose2[IRL - 2] = Math.min(dose2[IRL - 2] / SCORE_TEMP * 100.,
						99.9);
			} else {
				dose2[IRL - 2] = 99.9;
			}

			SCORE_TEMP = kerma[IRL - 2] / SCORE_NORM_NUM;
			kerma2[IRL - 2] = kerma2[IRL - 2] / SCORE_NORM_NUM;
			kerma2[IRL - 2] = (kerma2[IRL - 2] - SCORE_TEMP * SCORE_TEMP)
					/ (SCORE_NORM_NUM - 1);
			if (kerma2[IRL - 2] >= 0.)
				kerma2[IRL - 2] = Math.sqrt(kerma2[IRL - 2]);
			if (SCORE_TEMP != 0.) {
				kerma2[IRL - 2] = Math.min(kerma2[IRL - 2] / SCORE_TEMP * 100.,
						99.9);
			} else {
				kerma2[IRL - 2] = 99.9;
			}

			// "now analyze the uncertainty on the dose/kerma ratio"
			dosetokerma2[IRL - 2] = dosetokerma2[IRL - 2] / SCORE_NORM_NUM
					- dose[IRL - 2] * kerma[IRL - 2]
					/ (SCORE_NORM_NUM * SCORE_NORM_NUM);
			dosetokerma2[IRL - 2] = dosetokerma2[IRL - 2]
					/ (dose[IRL - 2] * kerma[IRL - 2] / (SCORE_NORM_NUM * SCORE_NORM_NUM));
			dosetokerma2[IRL - 2] = dosetokerma2[IRL - 2]
					/ (SCORE_NORM_NUM - 1);
			// "now estimate the uncertainty on dose/fluence"
			dosetokerma2[IRL - 2] = (dose2[IRL - 2] / 100.)
					* (dose2[IRL - 2] / 100.) + (kerma2[IRL - 2] / 100.)
					* (kerma2[IRL - 2] / 100.) - 2 * dosetokerma2[IRL - 2];
			if (dosetokerma2[IRL - 2] > 0.) {
				dosetokerma2[IRL - 2] = 100 * Math.sqrt(dosetokerma2[IRL - 2]);
			}
			if (dosetokerma2[IRL - 2] > 99.9) {
				dosetokerma2[IRL - 2] = 99.9;
			}
		}

		effdose = 0.0;
		for (IRL = 2; IRL <= nreg; IRL++) {
			dose[IRL - 2] = dose[IRL - 2] * 1.602E-10
					/ (mass[IRL - 2] * AINFLU);
			// if(kerma[IRL-2] != 0)
			// {
			kerma[IRL - 2] = kerma[IRL - 2] * 1.602E-10
					/ (mass[IRL - 2] * AINFLU);
			// }
			if (MONOEN != 0) {
				double exposure = airExposure * 1000;// mGy=>microGy
				double factor = 100 * exposure * xrs.getPhotonFlux()
						/ xrs.getAirKerma();
				// 100=1cm2=100mm2
				// flux=photons/(mm2*mAs);airKerma=microGy/mAs=>factor is
				// photons/cm2 so

				dosem[IRL - 2] = 1000 * dose[IRL - 2] * factor;// mGy
				dosem2[IRL - 2] = dose2[IRL - 2];

				effdose = effdose + wt[IRL - 2] * dosem[IRL - 2];
				effdose2 = effdose2
						+ Math.pow(wt[IRL - 2] * dosem[IRL - 2]
								* dosem2[IRL - 2] / 100., 2);
			}
		}

		if (MONOEN != 0) {
			if (effdose != 0) {
				effdose2 = (100 * Math.sqrt(effdose2)) / effdose;
			} else {
				effdose2 = 99.9;
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
		}

	}

	/**
	 * In general, AUSGAB is a routine which is called under a series 
	 * of well defined conditions specified by the value of IARG. Interface method.
	 * @param IARG the process code
	 */
	public void AUSGAB(int IARG) {

		// $IMPLICIT-NONE;

		double FTMP = 0.;
		int ip = 0;

		//double xsi = 0.0;
		double R1 = 0.0;
		double aux1 = 0.0;
		int IRL = 0;
		// "STACK OVERFLOW CHECK"
		if (EGS4Macro.ienhance == 1) {// " Option to enhance photon cross section in some region"
										// "write(6,*) 'in ausgab to recreate photon';"
			if (IARG == 15 || IARG == 17 || IARG == 19 || IARG == 23) {
				// "A pair/Compton/photoelectric/pair event is about to take place"
				EGS4.NP = EGS4.NP + 1; // "Boost the stack"
				if (EGS4.NP > EGS4.$MXSTACK) {
					EGS4.STOPPROGRAM = true;
					EGS4.seqStr = " ***************************************************"
							+ "  \n"
							+ " Calculation with CS-enhancement: unable to boost stack."
							+ "  \n"
							+ " ***************************************************";
					if (EGS4.iprint > 0)
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

				R1 = EGS4.random01();
				aux1 = 1. - 1. / EGS4Macro.cs_enhance_current;
				if (R1 > aux1) {
					EGS4.WT[EGS4.NP - 2] = 0.;
					EGS4.E[EGS4.NP - 2] = 0.;
					EGS4.DNEAR[EGS4.NP - 2] = -1.;
				}
				EGS4.WT[EGS4.NP - 1] = EGS4.WT[EGS4.NP - 1] * 1.0
						/ EGS4Macro.cs_enhance_current;
				// "write(6,*) ' Creating an unscattered photon!
				// ',ir(np),wt(np),wt(np-1);
				return;
			}// "end of block after photon events"

			// " Play Russian Roulette with scattered photons to avoid transport
			// " of many low weight particles.
			// " Now there is a stack pointer NPold which points to the particle "
			// " before the last discrete interaction. This change was necessary "
			// " for the implementation of atomic relaxations. So, check all particles"
			// " between NPold and NP for RR"
			if ((IARG == 18 // " Compton has occured"
					|| IARG == 20 // " After photo-absorption "
			|| IARG == 24) // " After Rayleigh "
					&& EGS4Macro.ienhance == 1) {
				EGS4Macro.ienhance = 0;
				// "write(6,*) ' iarg = ',iarg,' NP NPold: ',NP,NPold;

				for (ip = EGS4.NPold; ip <= EGS4.NP; ip++) {
					if (EGS4.IQ[ip - 1] == 0) {
						R1 = EGS4.random01();
						if (R1 < 1.0 / EGS4Macro.cs_enhance_current) {
							// "PARTICLE SURVIVES"
							EGS4.WT[ip - 1] = EGS4.WT[ip - 1]
									* EGS4Macro.cs_enhance_current;
							// "write(6,*) ' iarg = ',iarg,' particle has
							// survived ',
							// " wt(ip),ip;
						} else {
							EGS4.WT[ip - 1] = 0.;
							EGS4.E[ip - 1] = 0.;
							EGS4.DNEAR[ip - 1] = -1.;
							// "write(6,*) ' Have killed particle ',ip,
							// " ' after iarg = ',iarg;
						}
					}
				}
				return;
			}

		}

		// "Check if particle is leaving the transport geometry"
		IRL = EGS4.IR[EGS4.NP - 1]; // "local region number"
		if (IRL == 1) {
			return; // "outside the chamber, howfar will discard"
		}

		// "obtain frequently used local variables"
		IRL = EGS4.IR[EGS4.NP - 1];// IR(NP);
		if (IRL == 1)
			return; // "outside the chamber"
		if (IARG == 4 && !EGS4.BTEST_LATCH(EGS4.NP, 8)) {// "local energy deposition"
															// "include deposited energy as kerma"
			if (NHSTRY == kerma_last[IRL - 2]) {
				kerma_tmp[IRL - 2] = kerma_tmp[IRL - 2] + EGS4.WT[EGS4.NP - 1]
						* EGS4.EDEP;
			} else {
				kerma[IRL - 2] = kerma[IRL - 2] + kerma_tmp[IRL - 2];
				kerma2[IRL - 2] = kerma2[IRL - 2] + kerma_tmp[IRL - 2]
						* kerma_tmp[IRL - 2];

				kerma_tmpold[IRL - 2] = kerma_tmp[IRL - 2];
				kerma_lastold[IRL - 2] = kerma_last[IRL - 2];

				kerma_tmp[IRL - 2] = EGS4.WT[EGS4.NP - 1] * EGS4.EDEP;
				kerma_last[IRL - 2] = NHSTRY;
			}
		}
		if (IARG == 16) {// "pair event just occured"
			if (EGS4.NP > EGS4.NPold || EGS4.i_survived_RR > 0) {
				for (ip = EGS4.NPold; ip <= EGS4.NP; ip++) {
					if (EGS4.IQ[ip - 1] != 0 && !EGS4.BTEST_LATCH(ip, 8)) {
						if (NHSTRY == kerma_last[IRL - 2]) {
							kerma_tmp[IRL - 2] = kerma_tmp[IRL - 2]
									+ EGS4.WT[ip - 1]
									* (EGS4.E[ip - 1] - EGS4.PRM);
						} else {
							kerma[IRL - 2] = kerma[IRL - 2]
									+ kerma_tmp[IRL - 2];
							kerma2[IRL - 2] = kerma2[IRL - 2]
									+ kerma_tmp[IRL - 2] * kerma_tmp[IRL - 2];

							kerma_tmpold[IRL - 2] = kerma_tmp[IRL - 2];
							kerma_lastold[IRL - 2] = kerma_last[IRL - 2];

							kerma_tmp[IRL - 2] = EGS4.WT[ip - 1]
									* (EGS4.E[ip - 1] - EGS4.PRM);
							kerma_last[IRL - 2] = NHSTRY;
						}

						EGS4.LATCH[ip - 1] = EGS4.IBSET_LATCH(ip, 8);
					}
				}
			}
		}// "end of pair case"
		if (IARG == 18) {// "compton event just occured"
							// "must score kerma for all resultant electrons"
			if (EGS4.NP > EGS4.NPold) {// "compton occurred and we have not cleared the stack with"
										// "russian roulette"
				for (ip = EGS4.NPold; ip <= EGS4.NP; ip++) {
					if (EGS4.IQ[ip - 1] != 0 && !EGS4.BTEST_LATCH(ip, 8)) {
						if (NHSTRY == kerma_last[IRL - 2]) {
							kerma_tmp[IRL - 2] = kerma_tmp[IRL - 2]
									+ EGS4.WT[ip - 1]
									* (EGS4.E[ip - 1] - EGS4.PRM);
						} else {
							kerma[IRL - 2] = kerma[IRL - 2]
									+ kerma_tmp[IRL - 2];
							kerma2[IRL - 2] = kerma2[IRL - 2]
									+ kerma_tmp[IRL - 2] * kerma_tmp[IRL - 2];

							kerma_tmpold[IRL - 2] = kerma_tmp[IRL - 2];
							kerma_lastold[IRL - 2] = kerma_last[IRL - 2];

							kerma_tmp[IRL - 2] = EGS4.WT[ip - 1]
									* (EGS4.E[ip - 1] - EGS4.PRM);
							kerma_last[IRL - 2] = NHSTRY;
						}
						EGS4.LATCH[ip - 1] = EGS4.IBSET_LATCH(ip, 8);
					}
				}
			}
		}// "end of compton case"
		if (IARG == 20) {// "photoelectric event just occured"
			for (ip = EGS4.NPold; ip <= EGS4.NP; ip++) {
				if (EGS4.IQ[ip - 1] != 0 && !EGS4.BTEST_LATCH(ip, 8)) {
					if (NHSTRY == kerma_last[IRL - 2]) {
						kerma_tmp[IRL - 2] = kerma_tmp[IRL - 2]
								+ EGS4.WT[ip - 1] * (EGS4.E[ip - 1] - EGS4.PRM);
					} else {
						kerma[IRL - 2] = kerma[IRL - 2] + kerma_tmp[IRL - 2];
						kerma2[IRL - 2] = kerma2[IRL - 2] + kerma_tmp[IRL - 2]
								* kerma_tmp[IRL - 2];

						kerma_tmpold[IRL - 2] = kerma_tmp[IRL - 2];
						kerma_lastold[IRL - 2] = kerma_last[IRL - 2];

						kerma_tmp[IRL - 2] = EGS4.WT[ip - 1]
								* (EGS4.E[ip - 1] - EGS4.PRM);
						kerma_last[IRL - 2] = NHSTRY;
					}
					EGS4.LATCH[ip - 1] = EGS4.IBSET_LATCH(ip, 8);
				}
			}
		}// "end of photoelectric case"

		// "do some basic checks to see if scoring is needed"
		if (IARG >= 5 || EGS4.EDEP == 0)
			return;

		// "score total energy deposited"
		// "============================="

		FTMP = EGS4.WT[EGS4.NP - 1] * EGS4.EDEP;// System.out.println(FTMP);
		// score dose every where
		if (NHSTRY == dose_last[IRL - 2]) {
			dose_tmp[IRL - 2] = dose_tmp[IRL - 2] + FTMP;
		} else {
			dose[IRL - 2] = dose[IRL - 2] + dose_tmp[IRL - 2];
			dose2[IRL - 2] = dose2[IRL - 2] + dose_tmp[IRL - 2]
					* dose_tmp[IRL - 2];

			if (dose_last[IRL - 2] == kerma_lastold[IRL - 2]) {
				dosetokerma2[IRL - 2] = dosetokerma2[IRL - 2]
						+ dose_tmp[IRL - 2] * kerma_tmpold[IRL - 2];
			}
			dose_tmp[IRL - 2] = FTMP;
			dose_last[IRL - 2] = NHSTRY;
		}

		return;
	}// "END OF AUSGAB"
	// "*********************************************************************"
	/**
	 * The following is a general specification of HOWNEAR: 
	 * Given a particle at (x,y,z) in region irl, HOWNEAR answers the 
	 * question, What is the distance tperp to the closest boundary? Interface method.
	 */
	public void HOWNEAR() {
		// "Local variables
		// double r=0.0;
		double z = EGS4.Z[EGS4.NP - 1];
		double y = EGS4.Y[EGS4.NP - 1];
		double x = EGS4.X[EGS4.NP - 1];
		int IRL = EGS4.IR[EGS4.NP - 1];

		// r = Math.sqrt(x*x + y*y);
		// EGS4.tperp = EGS4.min(z-0.0,breastThickness-z,breastRadius-r);//or

		// new box aprox method:NO USED.IMPROVE IN interaction sampling-EGS-
		// lack in
		// geomtry approximation=>must be exact otherwise old sim is used!!!
		EGS4.tperp = phantom.getTperp(x, y, z, IRL);

		return;

	}// "end of subroutine HOWNEAR"

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
		double z = EGS4.Z[EGS4.NP - 1];
		double y = EGS4.Y[EGS4.NP - 1];
		double x = EGS4.X[EGS4.NP - 1];
		double u = EGS4.U[EGS4.NP - 1];
		double v = EGS4.V[EGS4.NP - 1];
		double w = EGS4.W[EGS4.NP - 1];

		// "DISCARD ZERO WEIGHT PARTICLES"
		if (EGS4.WT[EGS4.NP - 1] == 0.0) {
			EGS4.IDISC = 1;
			return;
		}

		int IRL = EGS4.IR[EGS4.NP - 1];// "LOCAL REGION NUMBER"
		// "DISCARD IF PARTICLE WANTS TO LEAVE THE GEOMETRY OR OF THE REGION IS TOTALLY"
		// "ABSORBING"
		if ((IRL == 1)) {
			EGS4.IDISC = 1;
			return;
		}

		double tplane = phantom.getMamoTz(w, z);
		double tcyl = phantom.getMamoTxy(u, v, x, y);
		// "DO MOST PROBABLE CASE FIRST WHERE A PLANE AND A CYLINDER CAN BE HIT"
		if ((Phantom.ihitp != 0) && (Phantom.ihitc != 0)) {
			if (tplane < tcyl) {// "HITS PLANE FIRST"
				if (tplane <= EGS4.USTEP) {
					EGS4.USTEP = tplane;
					EGS4.IRNEW = 1;
				}
				return;
			} else if (tcyl < tplane) {// "HITS CYLINDER FIRST"
				if (tcyl <= EGS4.USTEP) {
					EGS4.USTEP = tcyl;
					EGS4.IRNEW = 1;
				}
				return;
			} else {// "ODD CASE TCYL=TPLANE:HITS PLANE AND CYLINDER TOGETHER"
				if (tcyl <= EGS4.USTEP) {
					EGS4.USTEP = tcyl;
					EGS4.IRNEW = 1;
				}
				return;
			}
		}
		// "DO ODD CASE-PARTICLE CAN HIT PLANE BUT NOT CYLINDER"
		else if (Phantom.ihitp != 0) {
			if (tplane <= EGS4.USTEP) {
				EGS4.USTEP = tplane;
				EGS4.IRNEW = 1;
			}
			return;
		}
		// "DO ODD CASE-PARTICLE CAN HIT CYLINDER BUT NOT PLANE"
		else {
			if (tcyl <= EGS4.USTEP) {
				EGS4.USTEP = tcyl;
				EGS4.IRNEW = 1;
			}
			return;
		}
		// "AT THIS STAGE ALL GEOMETRICAL POSSIBILITIES HAVE BEEN CHECKED AND CONTROL"
		// "HAS ALREADY BEEN TRANSFERRED TO EGS"
	}// "END OF SUBROUTINE HOWFAR"

	/**
	 * Minimum of 3 doubles.
	 * @param a a
	 * @param b b
	 * @param c c
	 * @return the result
	 */
	public static double min(double a, double b, double c) {
		double r = Math.min(a, b);
		if (a < b)
			imin = 0;// first
		else
			imin = 1;// second
		if (r > c)
			imin = 2;// third
		r = Math.min(r, c);
		return r;
	}

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
		EGS4Core.SHOWER(iqin, ein, xin, yin, zin, uin, vin, win, irin, WEIGHT);
	}

	/**
	 * Setup input variables.
	 */
	private void inputs() {
		EGS4Macro.ICSDA = ETRANS_NORMAL;

		// NCASE=10000;

		EGS4Macro.IFULL = 0;
		STATLM = 0.0;

		// ######################################################
		nreg = 2;// 2 regions//allways
		EGS4.MEDIA[0] = "tissue_soft_icrp";// allways
		int[] MEDNUM = new int[EGS4.$MXMED];// allways
		MEDNUM[0] = 1;// first medium which is tissue//allways

		EGS4.MED[0] = 0;// allways
		for (int I = 2; I <= nreg; I++) {
			EGS4.MED[I - 1] = 1;
		} // "defaults"//allways

		EGS4.MED[1] = MEDNUM[0];// med is related to regions//allways

		EGS4.seqStr = " X-ray source on symmetry z-axis ";
		if (EGS4.iprint > 1)
			printSequence(EGS4.seqStr);
		EGS4.seqStr = " Distance to entry front face (cm) = "
				+ EGS4.format(fsd, 8, true);
		if (EGS4.iprint > 1)
			printSequence(EGS4.seqStr);
		EGS4.seqStr = " Radius of hitting surface (cm) = "
				+ EGS4.format(breastRadius, 8, true);
		if (EGS4.iprint > 1)
			printSequence(EGS4.seqStr);

		// get_transport_parameter();
		setEcutRegion = false;
		setPcutRegion = false;
		setSmaxRegion = false;
		ecut = 0.521;// allways
		pcut = 0.001;// allways
		smax = 1.e10;// allways

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

		incoh = incoh_ON;// allways
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
		coh = coh_ON;// allways
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
		relax = relax_ON;// allways
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

		pe = pe_ang_ON;// allways
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
		EGS4.ibrdst = brems_ang_KM;// allways
		EGS4.iprdst = pair_ang_SIMPLE;// allways

		EGS4.estepe = 0.25;// allways
		EGS4.ximax = 0.5;// allways

		EGS4.bca_algorithm = BCA_EXACT;// exact means
										// =EGS4.$BCA_ALGORITHM_DEFAULT=0;//allways
		EGS4.skindepth_for_bca = 3.0;
		EGS4.transport_algorithm = estep_alg_PRESTA_II;// 0;//allways
		ispin = spin_ON;// allways
		EGS4.ibr_nist = brems_cross_BH;// allways
		EGS4.pair_nrc = pair_cross_BH;// allways
		EGS4.eii_flag = eii_ON;// allways so=>
		EGS4.eiifile = "eii_ik";
		EGS4.photon_xsections = photon_xsection;// "" default no EPDL or XCOM
		// ==================================
		EGS4.itriplet = triplet_ON;// allways
		EGS4.radc_flag = radc_ON;// allways

		// VARIANCE REDUCTION:
		// @BREM SPLITTING
		EGS4.nbr_split = 1;// $MAXBRSPLIT//allways
		EGS4.i_play_RR = 0;// allways
		// ELECTRON RANGE REJECTION
		EGS4Macro.irejct = irejct_OFF;// allways
		ESAVEIN = 2.0;// #total energy below which range rejection is
						// considered//allways
		EGS4Macro.CEXPTR = 0.000;// allways
		IFARCE = IFARCE_OFF;// allways
		EGS4Macro.NFMIN = 1;// #Start forcing at this interaction
							// number//allways
		EGS4Macro.NFMAX = 1;// #Number of photon interactions after
							// which//allways
		phsplitt = 1;// allways

		EGS4Macro.cs_enhance = 1.0;// #Photon cross section scaling
									// factors//allways
		nENHREG = 1;// allways
		NENHLO[0] = 1;// allways
		NENHHI[0] = EGS4.$MXREG;// all geometry//allways

		test_inputs();
		if (EGS4.STOPPROGRAM) {
			closeFile();
			return;
		}

		NCASEO = 0;
		NCASET = 0;

		for (int IRL = 2; IRL <= nreg; IRL++) {
			dose[IRL - 2] = 0.0;
			dose2[IRL - 2] = 0.0;
			dose_tmp[IRL - 2] = 0.0;
			dose_last[IRL - 2] = 0;
			kerma[IRL - 2] = 0.0;
			kerma2[IRL - 2] = 0.0;
			kerma_tmp[IRL - 2] = 0.0;
			kerma_last[IRL - 2] = 0;
			dosetokerma2[IRL - 2] = 0.0;
		}

		NCASET = NCASE + NCASEO;

		EGS4.seqStr = " ********* SUCCESSFUL INPUT ACCOMPLISHED *********";
		if (EGS4.iprint > 1)
			printSequence(EGS4.seqStr);

	}// inputs

	/**
	 * Validate more inputs. Called by inputs routine.
	 */
	private void test_inputs() {
		if (EGS4Macro.ICSDA < 0 || EGS4Macro.ICSDA > 1) {
			EGS4Macro.ICSDA = 0;
		}// default

		if (NCASE < NCASE_MIN || NCASE > NCASE_MAX) {
			NCASE = NCASE_DEFAULT;
		}

		if (EGS4Macro.IFULL < 0 || EGS4Macro.IFULL > 0) {
			EGS4Macro.IFULL = 0;
		}// default
		if (STATLM < STATLM_MIN || STATLM > STATLM_MAX) {
			STATLM = STATLM_DEFAULT;
		}

		if (NCASE < $NCASEMIN) {
			NCASE = $NCASEMIN;
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
		// =================================
		if ((EGS4.nbr_split > $MAXBRSPLIT))
			EGS4.nbr_split = $MAXBRSPLIT;
		// CHARGED PARTICLE RUSSIAN ROULETTE
		if (EGS4.i_play_RR == 1) {
			EGS4.prob_RR = 1. / dble(EGS4.nbr_split);
		} else {
			EGS4.prob_RR = 1.;
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
			for (int i = 1; i <= nreg; i++) {
				EGS4.i_do_rr[i - 1] = 1;
				EGS4.e_max_rr[i - 1] = ESAVEIN;
			}
			// "note  e_max_r is total energy"
			// "above two arrays needed for each region for EGSnrc RANGE-DISCARD macro"
		}

		if (EGS4Macro.CEXPTR < EXPC_MIN || EGS4Macro.CEXPTR > EXPC_MAX) {
			EGS4Macro.CEXPTR = EXPC_DEFAULT;
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
		NFMIN_MAX = nreg;
		NFMAX_MAX = nreg + 1;
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

		ics_enhance = 0;
		for (int jj = 1; jj <= nreg; jj++) {
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
		for (int jj = 2; jj <= nreg; jj++) {
			COUNT = COUNT + EGS4Macro.iefl[jj - 1];
		}
		// "We don't care about region 1 since outside geometry"
		if (COUNT > 0 && (EGS4Macro.cs_enhance > 1.0001)) {// "there is enhancement somewhere"
			EGS4.seqStr = " Cross section enhancement in regions ";
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);
			for (int jj = 1; jj <= nreg; jj++) {
				EGS4.seqStr = " 		" + EGS4.format(EGS4Macro.iefl[jj - 1], 4);
				if (EGS4.iprint > 1)
					printSequence(EGS4.seqStr);
			}
			EGS4.seqStr = " Cross section enhancement factor: "
					+ EGS4.format(EGS4Macro.cs_enhance, 6, true);
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);
			ics_enhance = 1;
		} else {
			EGS4.seqStr = " No cross section enhancement";
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);

			ics_enhance = 0;
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
	 * Print output summary
	 */
	private void OSUMRY() {
		String s = "";
		int ll = 0;

		int IRL = 0;

		IRL = 2;

		s = " Total dose (Gray/(incident fluence)):";
		ll = s.length();
		ll = 50 - ll;
		s = s + EGS4.format("", ll);
		EGS4.seqStr = s + EGS4.format(dose[IRL - 2], 10, false) + " +/- "
				+ EGS4.format(dose2[IRL - 2], 6, true) + "%" + " \n";
		if (EGS4.iprint > 0)
			printSequence(EGS4.seqStr);
		if (MONOEN != 0) {
			s = " Total dose includes stoppers and discards so it can be greater than kerma.";
			EGS4.seqStr = s;
			if (EGS4.iprint > 0)
				printSequence(EGS4.seqStr);
			s = " This is not important for our purpose.";
			EGS4.seqStr = s;
			if (EGS4.iprint > 0)
				printSequence(EGS4.seqStr);
		}

		s = " Total kerma (Gray/(incident fluence)):";
		ll = s.length();
		ll = 50 - ll;
		s = s + EGS4.format("", ll);
		EGS4.seqStr = s + EGS4.format(kerma[IRL - 2], 10, false) + " +/- "
				+ EGS4.format(kerma2[IRL - 2], 6, true) + "%" + " \n";
		if (EGS4.iprint > 0)
			printSequence(EGS4.seqStr);

		s = " Dose/kerma: ";
		ll = s.length();
		ll = 50 - ll;
		s = s + EGS4.format("", ll);
		double rap = 0.0;
		if (kerma[IRL - 2] != 0.0) {
			rap = dose[IRL - 2] / kerma[IRL - 2];
		}
		EGS4.seqStr = s + EGS4.format(rap, 10, false) + " +/- "
				+ EGS4.format(dosetokerma2[IRL - 2], 6, true) + "%" + " \n";
		if (EGS4.iprint > 0)
			printSequence(EGS4.seqStr);

		// 100=1cm2=100mm2
		// flux=photons/(mm2*mAs);airKerma=microGy/mAs=>factor is photons/cm2 so
		if (MONOEN != 0) {
			s = " Total absorbed dose (mGy):";
			ll = s.length();
			ll = 50 - ll;
			s = s + EGS4.format("", ll);
			EGS4.seqStr = s + EGS4.format(dosem[IRL - 2], 10, false) + " +/- "
					+ EGS4.format(dosem2[IRL - 2], 6, true) + "%" + " \n";
			if (EGS4.iprint > 0)
				printSequence(EGS4.seqStr);

			s = " Total effective dose (mSv):";
			ll = s.length();
			ll = 50 - ll;
			s = s + EGS4.format("", ll);
			EGS4.seqStr = s + EGS4.format(effdose, 10, false) + " +/- "
					+ EGS4.format(effdose2, 6, true) + "%" + " \n";
			if (EGS4.iprint > 0)
				printSequence(EGS4.seqStr);
		}
		return;

	}

	/**
	 * Fix some source and energy related variables. Called by init. Position, direction 
	 * cosines, particle weight and incident fluence are all fixed here.
	 */
	private void fixEmAll() {
	//	double D = 0.0;
	//	double R2 = 0.0;
	//	int IXIN = 0;

		NHSTRY = NHSTRY + 1;// @@@@@@@@@@@@@@@@@@@@@@@@@

		// AINFLU=dble(NCASET)/(4*Math.PI*fsd*fsd);
		AINFLU = dble(NCASET) / (Math.PI * breastRadius * breastRadius);// real
																		// fluence
																		// NOT
																		// isotropic
																		// point
																		// source

		double diam = breastRadius * 0.999;
		double costmax = fsd / Math.sqrt(fsd * fsd + diam * diam);
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
		// initial we have formal directional cosines u,v,w = 0,0,-1 and then we
		// have teta and phi2:
		double u = Math.sin(teta) * Math.cos(phi2);// @@@@@@@@@@@@@@@@@@@@@@@@@@@
		double v = Math.sin(teta) * Math.sin(phi2);// @@@@@@@@@@@@@@@@@@@@@@@@@@@
		double w = costet;// <0//@@@@@@@@@@@@@@@@@@@@@@@@@@@
		double parcurs = fsd / costet;
		parcurs = Math.abs(parcurs);// >0
		double x0 = 0.0;
		double y0 = 0.0;
		double x = x0 + u * parcurs;// @@@@@@@@@@@@@@@@@@@@@@@@@@@
		double y = y0 + v * parcurs;// @@@@@@@@@@@@@@@@@@@@@@@@@@@
		double z = 0.0;// @@@@@@@@@@@@@@@@@@@@@@@@@@@
		xin = x;
		yin = y;

		irin = 2;// allways
		zin = z;

		uin = u;// EGS4SrcEns.xin/D;//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
		vin = v;// EGS4SrcEns.yin/D;//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
		win = w;// EGS4SrcEns.DISTZ/D;//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

		// ##################fix weight
	//	double us1 = 2 * Math.PI
	//			* (1 - fsd / Math.sqrt(fsd * fsd + diam * diam));

		// WEIGHT=us1/(4*Math.PI);//FIXED WEIGHT!!!@@@@@@@@@@@@@@
		WEIGHT = 1.0;// flux collimated beam!!
		// hit the detector allways
		indet = true;
	}

	/**
	 * Convert an int to a double to avoid integer division.
	 * @param n n
	 * @return the result
	 */
	public static double dble(int n) {
		Integer in = new Integer(n);
		return in.doubleValue();
	}

	/**
	 * Called by init. Initialize variables related to XRay spectrum.
	 */
	private void ENSRC() {
		// ##########################################################################################
		if (MONOEN == 0) {
			ein = E_inc;// 0.03;//30kev
			return;
		}
		// ###########################################################################################

		xrs = new XRaySpectrum(kvp, filtration, uanod);
		// ======================================
		ENSRCD = xrs.getXRayEnergies();
		SRCPDF = xrs.getXRayIntensities();
		double norm = xrs.getNormalizedValue();

		double sume1 = 0.0;
		double sume = 0.0;
		int NENSRC = ENSRCD.length - 1;

		for (int i = 0; i < ENSRCD.length; i++) {
			SRCPDF[i] = SRCPDF[i] / norm;
			// System.out.println("E "+ENSRCD[i]+"   "+"w "+SRCPDF[i]);
		}
		// ======================================
		for (int ib = 1; ib <= NENSRC; ib++) {
			sume = sume + SRCPDF[ib - 1];
			sume1 = sume1 + SRCPDF[ib - 1] * (ENSRCD[ib] + ENSRCD[ib - 1])
					/ 2.0;
		}
		// System.out.println("  @@@@@@@  "+sume);
		double rap = sume1 / sume;
		EGS4.seqStr = " Average spectrum energy is "
				+ EGS4.format(rap / 1000., 10, true) + " MeV";
		if (EGS4.iprint > 1)
			printSequence(EGS4.seqStr);

		ein = ENSRCD[ENSRCD.length - 1] / 1000.0;// "SET TO MAX ENERGY FOR SOME CHECKS"
	}

	/**
	 * Called by init. Prepare alias table for sampling incident energy.
	 */
	public static void ENSRC1() {
		if (MONOEN == 0)
			return;

		// prepare_alias_sampling(nensrc,srcpdf,srcpdf_at,srcbin_at);
		// 300 (srcpdf_at) is sufficient to cover spectrum energy in
		// radiodiagnostic
		// EGS4.prepare_alias_sampling(ENSRCD.length-1,SRCPDF,srcpdf_at,srcbin_at);

		EGS4.prepare_alias_table(ENSRCD.length - 1, ENSRCD, SRCPDF, srcpdf_at,
				srcbin_at);
	}

	/**
	 * Called by init in SHOWER loop. Incident energy is sampled here.
	 * @return the result
	 */
	public static double ENSRCH() {
		double result = 0.0;

		if (MONOEN == 0) {
			result = ein;
			return result;
		}

		// result =
		// EGS4.alias_sample(ENSRCD.length-1,ENSRCD,srcpdf_at,srcbin_at);

		result = EGS4.alias_sample1(ENSRCD.length - 1, ENSRCD, SRCPDF,
				srcpdf_at, srcbin_at);
		// System.out.println(" E: "+result);
		return result;
	}
}
