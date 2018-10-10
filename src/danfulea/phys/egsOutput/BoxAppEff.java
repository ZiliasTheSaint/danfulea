package danfulea.phys.egsOutput;

import java.io.FileWriter;
import java.util.Calendar;
import java.util.Date;

import javax.swing.JTextArea;

import danfulea.phys.egs.EGS4;
import danfulea.phys.egs.EGS4Core;
import danfulea.phys.egs.EGS4Geom;
import danfulea.phys.egs.EGS4Grid;
import danfulea.phys.egs.EGS4Macro;
import danfulea.phys.egs.EGS4SrcEns;
import danfulea.phys.egs.EgsQuestion;

/**
 * Class for quick evaluation of an ideal surface detector efficiency (XDET,YDET). 
 * This detection area is on a X,Y BOX surface. The box is uniformly filled with air,  
 * and the beta nuclides are uniformly distributed in X,Y,Z Box volume. Optional, the 
 * source of electrons can also be a point source placed inside the box.
 * @author Dan Fulea, 03 AUG. 2006
 */
// a quick evaluation of an ideal surface detector efficiency (XDET,YDET)
// this detection area is on a X,Y BOX surface. The box is uniformely filled
// with AIR,
// and the beta nuclides are uniformely distributed in X,Y,Z Box volume.
public class BoxAppEff implements EgsQuestion {
	public static JTextArea jta;
	public static boolean systemOutB = true;

	public static int SOURCE = -1;
	public static final int SOURCE_POINT = 0;// means point source on symmetry
												// axis
	public static final int SOURCE_SARPAGAN = 1;// means volumic
												// SARPAGAN(cylinder) source on
												// top of the detecor

	public static int $NBATCH = 10;// "OUTPUT BATCHES                             "
	public static int JCASE = 0;; // "no. of histories per batch"

	public static int $MAXBRSPLIT = 200;// "MAX BREM SPLITTING NUMBER"

	public static int NCASE = 0;
	public static int NCASEO = 0;
	public static int NCASET = 0;

	public static double ETHRESHOLD = 0.0;
	public static double eefficiency = 0.0;
	public static double eefficiency_error = 0.0;
	// COUNTER FOR TOTAL NUMBER OF HISTORIES SUCCESSFULLY SIMULATED
	public static int IHSTRY = 0;
	public static double WT1OLD = 0.0;
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
	public static final int IFULL_ENTRANCE_REGIONS = 1;
	public static final int IFULL_PULSE_HEIGHT_DISTRIBUTION = 2;
	public static final int IFULL_SCATTER_FRACTION = 3;
	public static final int IFULL_OFMET_Fricke = 4;
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

	public static int ISUMX = 0;

	public static double SCORE_NORM_NUM = 0.0;
	public static double SCORE_TEMP = 0.0;
	public static double SCORE_TEMP2 = 0.0;

	public static int IRL = 0;

	public static int ics_enhance = 0;
	// ========================
	public static int IPRINT = 2;
	public static boolean createOutputFile = false;
	public static boolean putInFile = false;// internal var defining when and
											// what to print
	private String filename = "";
	FileWriter sigfos;

	private boolean ekinokB = true;
	// -------------end INPUTS-------------------------------------------

	// =================================================================
	public static double ZBOUND = 0.0;// ZBOX
	public static double YBOUND = 0.0;// YBOX
	public static double XBOUND = 0.0;// XBOX
	public static double XDET = 0.0;// XDET<XBOX
	public static double YDET = 0.0;// YDET<YBOX
	public static double YDETCENTER = 0.0;// YDET<YBOX
	public static double XDETCENTER = 0.0;
	public static double ps_distance = 0.0;// point source distance on z axis to
											// the CENTERDET
	private static double pondt = 0.0;
	private static double pondt2 = 0.0;
	public static String airmaterial = "";

	/**
	 * Constructor.
	 */
	public BoxAppEff() {
		// ===DE MUTAT IN EXTERN===================================
		/*
		 * BoxAppEff.reset(); BoxAppEff.systemOutB=true;
		 * BoxAppEff.ZBOUND=12.0;//cm BoxAppEff.XBOUND=12.0;//cm
		 * BoxAppEff.YBOUND=12.0;//cm //ZCENTERDET=6.0;//cm,at 1/2
		 * //DET.CENTER=AT SIDE CENTER (XBOUND/2, YBOUND/2)
		 * BoxAppEff.XDET=6.0;//cm BoxAppEff.YDET=6.0;//cm
		 * 
		 * BoxAppEff.XDETCENTER=0;//deviation from center XBOUND/2.;
		 * BoxAppEff.YDETCENTER=0;//deviation from center YBOUND/2.;
		 * 
		 * BoxAppEff.ps_distance=6.2;//cm
		 * 
		 * BoxAppEff.NCASE=20000; BoxAppEff.airmaterial="air_dry_nearsealevel";
		 * BoxAppEff.ETHRESHOLD=0.0;
		 * 
		 * //BoxAppEff.SOURCE=BoxAppEff.SOURCE_SARPAGAN;
		 * BoxAppEff.SOURCE=BoxAppEff.SOURCE_POINT;
		 * 
		 * EGS4SrcEns.ipart=EGS4SrcEns.ipart_electron;//or
		 * //EGS4SrcEns.ipart=EGS4SrcEns.ipart_positron;
		 * 
		 * EGS4.setMXMED(1);//"MAX # OF MEDIA 		" EGS4.setMXREG(3);
		 * EGS4.setMXSTACK
		 * (4000);//"NEED HUGE STACK FOR CORRELATIONS+splitting   "
		 * EGS4.setMXRANGE(500); //"for range arrays used in range_rejection()"
		 * EGS4.$MXMDSH=100;//allways
		 * 
		 * EGS4.reset(); EGS4Core.reset(); EGS4Geom.reset(); EGS4Grid.reset();
		 * EGS4Macro.reset(); EGS4SrcEns.reset();
		 * EGS4.egs_set_defaults();//first default then:
		 * 
		 * //EGS4SrcEns.monoindex=EGS4SrcEns.iMONOENERGETIC;
		 * //EGS4SrcEns.ikemev=2.662; EGS4SrcEns.monoindex=EGS4SrcEns.iSPECTRUM;
		 * EGS4SrcEns.enerFilename="sr90y90";//spectrumTf.getText();
		 */
		// ===============================
		EGS4.MEDIA[0] = airmaterial;
		EGS4.MED[0] = 0;
		EGS4.MED[2] = 0;
		EGS4.MED[1] = 1;

		BoxAppEff.IPRINT = 2;
		BoxAppEff.createOutputFile = false;
		EGS4.RandomUse = 1;
		EGS4.ranluxB = true;

		EGS4Macro.ICSDA = BoxAppEff.ETRANS_NORMAL;

		EGS4Macro.IFULL = BoxAppEff.IFULL_PULSE_HEIGHT_DISTRIBUTION;

		if (EGS4SrcEns.ipart == 4) {
			EGS4SrcEns.ipart = -1;
		}

		EGS4SrcEns.ISOURC = 0;

		BoxAppEff.incoh = BoxAppEff.incoh_ON;
		BoxAppEff.coh = BoxAppEff.coh_ON;
		BoxAppEff.relax = BoxAppEff.relax_ON;
		BoxAppEff.pe = BoxAppEff.pe_ang_ON;
		EGS4.ibrdst = BoxAppEff.brems_ang_KM;
		EGS4.iprdst = BoxAppEff.pair_ang_SIMPLE;
		EGS4.ibr_nist = BoxAppEff.brems_cross_BH;
		EGS4.pair_nrc = BoxAppEff.pair_cross_BH;
		BoxAppEff.ispin = BoxAppEff.spin_ON;
		EGS4.eii_flag = BoxAppEff.eii_ON;
		BoxAppEff.photon_xsection = "";
		EGS4.itriplet = BoxAppEff.triplet_ON;
		EGS4.radc_flag = BoxAppEff.radc_ON;
		EGS4.bca_algorithm = BoxAppEff.BCA_EXACT;
		EGS4.transport_algorithm = BoxAppEff.estep_alg_PRESTA_II;
		EGS4.nbr_split = 1;
		EGS4.i_play_RR = 0;
		EGS4Macro.irejct = BoxAppEff.irejct_OFF;
		BoxAppEff.ESAVEIN = 2.0;
		EGS4Macro.cs_enhance = 1.0;
		EGS4Macro.CEXPTR = 0.000;
		BoxAppEff.IFARCE = BoxAppEff.IFARCE_OFF;
		// ==============================
		// =====================================THIS MUST BE:============
		pondt = 0.0;

		putInFile = true;
		Calendar cal = Calendar.getInstance();
		String fs = cal.get(Calendar.YEAR) + "_" + cal.get(Calendar.MONTH) + "_"
				+ cal.get(Calendar.DAY_OF_MONTH) + "_" + cal.get(Calendar.HOUR) + "_"
				+ cal.get(Calendar.MINUTE) + "_" + cal.get(Calendar.SECOND) + ".txt";
		filename = fs;// will be 2005_10_25_14_30_56.txt

		ekinokB = true;
		// ==========================
		init();
	}

	/**
	 * Reset global variables for re-use.
	 */
	public static void reset() {
		SOURCE = -1;
		$NBATCH = 10;
		JCASE = 0;

		$MAXBRSPLIT = 200;
		NCASE = 0;
		NCASEO = 0;
		NCASET = 0;

		ETHRESHOLD = 0.0;
		eefficiency = 0.0;
		eefficiency_error = 0.0;
		IHSTRY = 0;
		WT1OLD = 0.0;

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
		phsplitt = 1;
		phsplitt_MAX = 0;
		EI = 0.0;
		EKMAX = 0.0;
		ISUMX = 0;

		SCORE_NORM_NUM = 0.0;
		SCORE_TEMP = 0.0;
		IRL = 0;

		SCORE_TEMP2 = 0.0;

		ics_enhance = 0;
		IPRINT = 2;

		pondt = 0.0;
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
		// --variable init
		EGS4SrcEns.$MXRDIST = 1000;
		EGS4SrcEns.$NENSRC = 300;
		EGS4Geom.$NVALUE = 100;

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

		// @@@@@@@@@@@@@@@@@@@@
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
		EGS4.ispmfp = EGS4.iDose;// select photon mfp
		EGS4.iurd = 0;// no user range discard
		EGS4.iraycorr = EGS4.iDose;// rayleigh correction
		EGS4.isemfp = EGS4.iDose;// $SELECT-ELECTRON-MFP
		EGS4.iGeom = 1;// EGS4.iCavity;//1=RZ GEOM
		EGS4Macro.ismfpfpb = EGS4.iDose;// select mfp parallel beam

		// -----------------
		EGS4.USER_CONTROLS_TSTEP_RECURSION = EGS4.iDose;// no effect here
		EGS4.hatchindex = EGS4.iDose;// no effect here
		// -----------------
		EGS4.iprint = IPRINT;// EGS4.iprint=2;//summary
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
		EGS4.seqStr = " **********************BOXGeometry: APPLICATION****************************************";
		if (EGS4.iprint > 1)
			printSequence(EGS4.seqStr);
		EGS4.seqStr = " COMPUTES THE IDEAL SURFACE DETECTOR EFFICIENCY FOR AN ELECTRON SOURCE HAVING A FINITE BOX GEOMETRY OR BEING A POINT SOURCE.";
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

		if (SOURCE == SOURCE_POINT || SOURCE == SOURCE_SARPAGAN) {
			fixSRCOTO();
		}

		if (NCASE / $NBATCH == 0) {
			NCASE = $NBATCH;
		}
		JCASE = NCASE / $NBATCH;
		NCASE = JCASE * $NBATCH;// "NUMBER OF HISTORIES PER BATCH

		IHSTRY = NCASEO; // "reset the number of histories counter"

		EGS4SrcEns.NHSTRY = 0; // "start the no. of primary histories counter at zero"

		// "set up ausgab calls"
		for (int J = 1; J <= 5; J++) {
			EGS4.iausfl[J - 1] = 1;// IAUSFL(J)=1;
		}
		for (int J = 6; J <= 25; J++) {
			EGS4.iausfl[J - 1] = 0;
		} // "NORMAL EXECUTION"
		EGS4.iausfl[5] = 1; // "AFTER TRANSPORT"

		for (int j = 1; j <= EGS4.$MXREG; j++) {
			EGS4Macro.iefl[j - 1] = 0;
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

		if (EGS4SrcEns.MONOEN == 0) {// "MONOENERGETIC INPUT BEAM"
			if (EGS4SrcEns.iqin == 0) {
				EI = EGS4SrcEns.ein;
			} else {
				EI = EGS4SrcEns.ein + EGS4.RM;
			}
			EKMAX = EGS4SrcEns.ein; // "MAXIMUM KINETIC ENERGY"
		} else if (EGS4SrcEns.MONOEN == 1) {// "ENERGY SPECTRUM"
			EGS4SrcEns.ENSRC1();// "NORMALIZE THE ENERGY DISTRIBUTION"
			EKMAX = EGS4SrcEns.ENSRCD[EGS4SrcEns.NENSRC];// "MAXIMUM KINETIC ENERGY IN THE SPECTRUM"
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

		EGS4Macro.do_fast_step = false;

		// "INITIALIZE DATA ARRAYS FOR FLUORESCENT X-RAYS IF NEEDED"
		ISUMX = 0;
		for (int JJ = 1; JJ <= EGS4Geom.nreg; JJ++) {
			ISUMX = ISUMX + EGS4.iedgfl[JJ - 1]; // "NON-ZERO IF X-RAYS ANYWHERE"
		}
		if (ISUMX != 0) {
			EGS4.EDGSET(2, EGS4Geom.nreg);
		}

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

			if (IBATCH == 1) {
				EGS4.seqStr = " BATCH" + EGS4.format("", 8) + "ELAPSED"
						+ EGS4.format("", 3) + "time of day";
				if (EGS4.iprint > 1)
					printSequence(EGS4.seqStr);
			}

			for (int ICASE = 1; ICASE <= JCASE; ICASE++) {// "now fill each IS bin"

				if (EGS4SrcEns.ISOURC != 23)
					IHSTRY = IHSTRY + 1; // "increment history counter"

				EGS4Macro.NFTIME = 0; // "reset the photon forced interaction counter"
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

				// "FOR AN INPUT ENERGY SPECTRUM, DETAILED FORCING MACRO IS USED"
				EGS4.LATCHI = 0;
				// ****************************FIX WEIGHT
				// *****************************************************
				if (SOURCE == SOURCE_POINT || SOURCE == SOURCE_SARPAGAN) {
					fixEmAll();
				}
				// ********************************************************************************************
				if (EGS4SrcEns.iqin != 0) {
					if (EGS4SrcEns.ein > ETHRESHOLD) {
						ekinokB = true;
					} else {
						ekinokB = false;
					}
				}

				if (ekinokB)
					SHOWER();
				if (EGS4.STOPPROGRAM) {
					closeFile();
					return;
				}

			}// "END OF THE ICASE LOOP"

			// "NOW DO STATS ON THE PEAK REGION"
			SCORE_NORM_NUM = EGS4SrcEns.dble(IHSTRY);

			String timePerBatch = EGS4.timeElapsedShort(startTime);
			Calendar call = Calendar.getInstance();
			String timeday = call.get(Calendar.HOUR) + ":" + call.get(Calendar.MINUTE)
					+ ":" + call.get(Calendar.SECOND);

			EGS4.seqStr = EGS4.format("", 2) + EGS4.format(IBATCH, 3)
					+ EGS4.format("", 2) + EGS4.format(timePerBatch, 14)
					+ EGS4.format("", 3) + EGS4.format(timeday, 9)
					+ EGS4.format("", 6);
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);

		}// "END OF SIMULATIONS"

		// "******************************************************************************
		// "
		// " *** SECTION 3 ***
		// "
		// "------------------------------------------------------------------------------
		//
		// "STATISTICAL AND OTHER DATA HANDLING AND CALL THE OUTPUT SUMMARY ROUTINE"
		//
		// "------------------------------------------------------------------------------
		SCORE_NORM_NUM = EGS4SrcEns.dble(IHSTRY);

		SCORE_TEMP2 = pondt2;
		SCORE_TEMP = pondt / SCORE_NORM_NUM;
		pondt2 = pondt2 / SCORE_NORM_NUM;
		pondt2 = (pondt2 - SCORE_TEMP * SCORE_TEMP) / (SCORE_NORM_NUM - 1);
		if (pondt2 >= 0.)
			pondt2 = Math.sqrt(pondt2);
		if (SCORE_TEMP != 0.) {
			pondt2 = Math.min(pondt2 / SCORE_TEMP * 100., 99.9);
		} else {
			pondt2 = 99.9;
		}
		eefficiency = pondt * 100 / SCORE_NORM_NUM;// %
		eefficiency_error = pondt2;
		// =================================================================================
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
	// "******************************************************************************

	/**
	 * In general, AUSGAB is a routine which is called under a series 
	 * of well defined conditions specified by the value of IARG. Interface method.
	 * @param IARG the process code
	 */
	public void AUSGAB(int IARG) {
		IRL = EGS4.IR[EGS4.NP - 1]; // "local region number"
		if (IRL == 1) {
			return; // "outside the chamber, howfar will discard"
		}

		if (IRL == 3 && IARG == 5
				&& (EGS4.IQ[EGS4.NP - 1] == -1 || EGS4.IQ[EGS4.NP - 1] == 1))// after
																				// transport
																				// electrons
																				// exit
																				// on
																				// the
																				// detector
																				// side
		{
			double y = EGS4.Y[EGS4.NP - 1];
			double x = EGS4.X[EGS4.NP - 1];
			WT1OLD = EGS4.WT[0];// primary weight
			if (x > XDETCENTER - XDET / 2. && x < XDETCENTER + XDET / 2
					&& y > YDETCENTER - YDET / 2. && x < YDETCENTER + YDET / 2) {
				// hit detector surface=>score it don't care about energy and if
				// the same history generates more electrons!
				// just all electrons hitting the ideal detector are scored!!
				pondt = pondt + WT1OLD;
				pondt2 = pondt2 + WT1OLD * WT1OLD;
			}

		}
	}// "END OF AUSGAB"

	// "*********************************************************************"
	// "
	// "
	// " ***********
	// " * *
	// " * HOWNEAR *
	// " * *
	// " ***********
	// "
	// "
	// "
	// "******************************************************************************

	/**
	 * The following is a general specification of HOWNEAR: 
	 * Given a particle at (x,y,z) in region irl, HOWNEAR answers the 
	 * question, What is the distance tperp to the closest boundary? Interface method.
	 */
	public void HOWNEAR() {
		// =====================================================================
		double z = EGS4.Z[EGS4.NP - 1];
		double y = EGS4.Y[EGS4.NP - 1];
		double x = EGS4.X[EGS4.NP - 1];
		int IRL = EGS4.IR[EGS4.NP - 1];// "LOCAL REGION NUMBER"

		if (SOURCE == SOURCE_POINT) {
			// "   REGION 1        |   REGION 2    |       REGION 3                  "
			// "                   |               |                                 "
			// "   e- =========>   |               | e- or photon ====>              "
			// "                   |               |                                 "
			// "   vacuum          |     AIR       |       vacuum                    "
			// "                                                                     "
			if (IRL == 3) {
				return;
			} else if (IRL == 2) {// "We are in the Ta plate - check the geometry"
				EGS4.tperp = Math.min(z, (ps_distance - z));
			} else if (IRL == 1) {
				return;
			}
		} else if (SOURCE == SOURCE_SARPAGAN) {
			if (IRL == 2)// inside the box
			{
				double t1 = x + 0.5 * XBOUND;
				double t2 = 0.5 * XBOUND - x;
				if (t1 < t2)
					EGS4.tperp = t1;
				else
					EGS4.tperp = t2;

				t1 = y + 0.5 * YBOUND;
				if (t1 < EGS4.tperp)
					EGS4.tperp = t1;
				t1 = 0.5 * YBOUND - y;
				if (t1 < EGS4.tperp)
					EGS4.tperp = t1;
				t1 = z + 0.5 * ZBOUND;
				if (t1 < EGS4.tperp)
					EGS4.tperp = t1;
				t1 = 0.5 * ZBOUND - z;
				if (t1 < EGS4.tperp)
					EGS4.tperp = t1;
				return;
			}
			double s1 = 0;
			double s2 = 0;
			int nout = 0;

			if (2. * x + XBOUND < 0) {
				double t = -0.5 * XBOUND - x;
				s1 += t;
				s2 += t * t;
				nout++;
			} else if (2. * x - XBOUND > 0) {
				double t = x - 0.5 * XBOUND;
				s1 += t;
				s2 += t * t;
				nout++;
			}
			if (2. * y + YBOUND < 0) {
				double t = -0.5 * YBOUND - y;
				s1 += t;
				s2 += t * t;
				nout++;
			} else if (2. * y - YBOUND > 0) {
				double t = y - 0.5 * YBOUND;
				s1 += t;
				s2 += t * t;
				nout++;
			}
			if (2. * z + ZBOUND < 0) {
				double t = -0.5 * ZBOUND - z;
				s1 += t;
				s2 += t * t;
				nout++;
			} else if (2. * z - ZBOUND > 0) {
				double t = z - 0.5 * ZBOUND;
				s1 += t;
				s2 += t * t;
				nout++;
			}
			if (nout < 2) {
				EGS4.tperp = s1;
				return;
			}
			EGS4.tperp = Math.sqrt(s2);
		}

	}// "end of subroutine HOWNEAR"

	// ;"******************************************************************************
	// "
	// " **********
	// " * *
	// " * HOWFAR *
	// " * *
	// " **********
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
		double z = EGS4.Z[EGS4.NP - 1];
		double y = EGS4.Y[EGS4.NP - 1];
		double x = EGS4.X[EGS4.NP - 1];
		double w = EGS4.W[EGS4.NP - 1];
		double v = EGS4.V[EGS4.NP - 1];
		double u = EGS4.U[EGS4.NP - 1];
		//int IQL = EGS4.IQ[EGS4.NP - 1];
		int IRL = EGS4.IR[EGS4.NP - 1];// "LOCAL REGION NUMBER"
		// "DISCARD IF PARTICLE WANTS TO LEAVE THE GEOMETRY OR OF THE REGION IS TOTALLY"
		// "ABSORBING"

		if (SOURCE == SOURCE_POINT) {
			double TVAL = 0.0;

			if (EGS4.IR[EGS4.NP - 1] == 3) {
				EGS4.IDISC = 1;
				return;// "terminate this history: it is past the plate"
			} else if (EGS4.IR[EGS4.NP - 1] == 2) {// "We are in the AIR plate - check the geometry"

				if (EGS4.W[EGS4.NP - 1] > 0.0) {
					// "going forward - consider first since  most frequent"
					TVAL = (ps_distance - EGS4.Z[EGS4.NP - 1])
							/ EGS4.W[EGS4.NP - 1];// "TVAL is dist to boundary"
					// "                          in this direction"
					if (TVAL > EGS4.USTEP) {
						return;// "can take currently requested step"
					} else {
						EGS4.USTEP = TVAL;
						EGS4.IRNEW = 3;
						return;
					}
				}// "END OF W(NP)>0 CASE"
				else if (EGS4.W[EGS4.NP - 1] < 0.0) {// "going back towards origin"
					TVAL = -EGS4.Z[EGS4.NP - 1] / EGS4.W[EGS4.NP - 1]; // "distance to plane at origin"
					if (TVAL > EGS4.USTEP) {
						return;// "can take currently requested step"
					} else {
						EGS4.USTEP = TVAL;
						EGS4.IRNEW = 1;
						return;
					}
				}// "END W(NP)<0 CASE"
				else if (EGS4.W[EGS4.NP - 1] == 0.0) {
					// "cannot hit boundary"
					return;
				}
			}// "end of region 2 case"
			else if (EGS4.IR[EGS4.NP - 1] == 1) {// "in region with source"
				if (EGS4.W[EGS4.NP - 1] > 0.0) {// "this must be a source particle on z=0 boundary"
					EGS4.USTEP = 0.0;
					EGS4.IRNEW = 2;
					return;
				} else {// "it must be a reflected particle-discard it"
					EGS4.IDISC = 1;
					return;
				}
			}// "end region 1 case"
		} else if (SOURCE == SOURCE_SARPAGAN) {
			if (IRL == 1) {
				EGS4.IDISC = 1;
				return;
			}

			if (IRL == 3) {
				EGS4.IDISC = 1;
				return;// "terminate this history: it is past the plate"
			}

			if (IRL == 2)// inside
			{
				double t1 = 1.e30;
				//int inew = 0;
				if (u > 0) {
					t1 = (XBOUND - 2. * x) / (2. * u);
				} else if (u < 0) {
					t1 = -(XBOUND + 2. * x) / (2. * u);
				}
				if (t1 < EGS4.USTEP) {
					EGS4.USTEP = t1;
					EGS4.IRNEW = 1;
					// return;
				}
				t1 = 1.e30;

				if (v > 0)
					t1 = (YBOUND - 2. * y) / (2. * v);
				else if (v < 0)
					t1 = -(YBOUND + 2. * y) / (2. * v);
				if (t1 < EGS4.USTEP) {
					EGS4.USTEP = t1;
					EGS4.IRNEW = 1;
				}

				t1 = 1.e30;
				if (w > 0) {
					t1 = (ZBOUND - 2. * z) / (2. * w);
					if (t1 < EGS4.USTEP) {
						EGS4.USTEP = t1;
						EGS4.IRNEW = 3;// detector side
					}
				} else if (w < 0) {
					t1 = -(ZBOUND + 2. * z) / (2. * w);
					if (t1 < EGS4.USTEP) {
						EGS4.USTEP = t1;
						EGS4.IRNEW = 1;
					}
				}

				return;
			}

			double t1 = 1.e30;
			if (2. * x + XBOUND < 0 && u > 0)
				t1 = -(2. * x + XBOUND) / (2. * u);
			else if (2. * x - XBOUND > 0 && u < 0)
				t1 = -(2 * x - XBOUND) / (2. * u);
			if (t1 < EGS4.USTEP) {
				double y1 = y + v * t1;
				double z1 = z + w * t1;
				if (2. * y1 + YBOUND >= 0 && 2. * y1 - YBOUND <= 0
						&& 2. * z1 + ZBOUND >= 0 && 2. * z1 - ZBOUND <= 0) {
					EGS4.USTEP = t1;
					EGS4.IRNEW = 2;// go inside the XYZBOUND box
					return;
				}
			}

			t1 = 1.e30;
			if (2. * y + YBOUND < 0 && v > 0)
				t1 = -(2. * y + YBOUND) / (2. * v);
			else if (2. * y - YBOUND > 0 && v < 0)
				t1 = -(2. * y - YBOUND) / (2. * v);
			if (t1 < EGS4.USTEP) {
				double x1 = x + u * t1;
				double z1 = z + w * t1;
				if (2. * x1 + XBOUND >= 0 && 2. * x1 - XBOUND <= 0
						&& 2. * z1 + ZBOUND >= 0 && 2. * z1 - ZBOUND <= 0) {
					EGS4.USTEP = t1;
					EGS4.IRNEW = 2;// go inside the XYZBOUND box
					return;
				}
			}

			t1 = 1.e30;
			if (2. * z + ZBOUND < 0 && w > 0)
				t1 = -(2. * z + ZBOUND) / (2. * w);
			else if (2. * z - ZBOUND > 0 && w < 0)
				t1 = -(2. * z - ZBOUND) / (2. * w);
			if (t1 < EGS4.USTEP) {
				double x1 = x + u * t1;
				double y1 = y + v * t1;
				if (2. * x1 + XBOUND >= 0 && 2. * x1 - XBOUND <= 0
						&& 2. * y1 + YBOUND >= 0 && 2. * y1 - YBOUND <= 0) {
					EGS4.USTEP = t1;
					EGS4.IRNEW = 2;// go inside the XYZBOUND box
					return;
				}
			}

			EGS4.IRNEW = 1;
			return;
		}// sarpagan
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
		EGS4Core.SHOWER(EGS4SrcEns.iqin, EI, EGS4SrcEns.xin, EGS4SrcEns.yin,
				EGS4SrcEns.zin, EGS4SrcEns.uin, EGS4SrcEns.vin, EGS4SrcEns.win,
				EGS4SrcEns.irin, EGS4SrcEns.WEIGHT);
	}

	/**
	 * Setup input variables.
	 */
	private void inputs() {

		if (SOURCE == SOURCE_POINT) {
			EGS4.seqStr = " Point source on symmetry z-axis ";
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);
			EGS4.seqStr = " Distance to entry front face (cm) = "
					+ EGS4.format(ps_distance, 8, true);
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);
			EGS4.seqStr = " Ideal surface Detector -X dimension (cm) = "
					+ EGS4.format(XDET, 8, true);
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);
			EGS4.seqStr = " Ideal surface Detector -Y dimension (cm) = "
					+ EGS4.format(YDET, 8, true);
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);

		} else if (SOURCE == SOURCE_SARPAGAN) {
			EGS4.seqStr = " Volumic (BOX) source on top of the detector house ";
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);
			EGS4.seqStr = " BOX source -X dimension (cm) = "
					+ EGS4.format(XBOUND, 8, true);
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);
			EGS4.seqStr = " BOX source -Y dimension (cm) = "
					+ EGS4.format(YBOUND, 8, true);
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);
			EGS4.seqStr = " BOX source -Z dimension (cm) = "
					+ EGS4.format(ZBOUND, 8, true);
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);
			EGS4.seqStr = " Ideal surface Detector -X dimension (cm) = "
					+ EGS4.format(XDET, 8, true);
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);
			EGS4.seqStr = " Ideal surface Detector -Y dimension (cm) = "
					+ EGS4.format(YDET, 8, true);
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);

		}

		EGS4SrcEns.ENSRC();
		if (EGS4.STOPPROGRAM) {
			closeFile();
			return;
		}
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
		EGS4.estepe = 0.25;
		EGS4.ximax = 0.5;
		EGS4.skindepth_for_bca = 3.0;

		test_inputs();
		if (EGS4.STOPPROGRAM) {
			closeFile();
			return;
		}

		NCASEO = 0;
		NCASET = 0;

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
		if (EGS4Macro.ICSDA < 0 || EGS4Macro.ICSDA > 1) {
			EGS4Macro.ICSDA = 0;
		}// default

		if (NCASE < NCASE_MIN || NCASE > NCASE_MAX) {
			NCASE = NCASE_DEFAULT;
		}

		if (EGS4Macro.IFULL < 0 || EGS4Macro.IFULL > 4) {
			EGS4Macro.IFULL = 0;
		}// default

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

		if (EGS4Macro.IFULL == 2 && EGS4.nbr_split > 1) {// "cannot have this"
			EGS4.seqStr = " You cannot calculate a pulse height distribution with";
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);
			EGS4.seqStr = " bremsstrahlung splitting on.  Will run simulation with";
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);
			EGS4.seqStr = " IFULL= dose and stoppers";
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);

			EGS4Macro.IFULL = 0;
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

		if (EGS4Macro.IFULL == 2 && EGS4Macro.CEXPTR != 0) {
			EGS4.seqStr = " You cannot calculate a pulse height distribution with";
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);
			EGS4.seqStr = " pathlength biasing.  Will run simulation with";
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);
			EGS4.seqStr = " IFULL= dose and stoppers";
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);

			EGS4Macro.IFULL = 0;
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

		if (EGS4Macro.IFULL == 2
				&& (EGS4Macro.NFMAX > EGS4Macro.NFMIN || EGS4Macro.NFMIN > 1)) {
			EGS4.seqStr = " You cannot calculate a pulse height distribution with";
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);
			EGS4.seqStr = " more than 1 interaction forced or if the forced interaction";
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);
			EGS4.seqStr = " is > the first interaction.  Will run simulation with";
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);
			EGS4.seqStr = " IFULL= dose and stoppers";
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);

			EGS4Macro.IFULL = 0;
		}

		ics_enhance = 0;
		for (int jj = 1; jj <= EGS4Geom.nreg; jj++) {
			EGS4Macro.iefl[jj - 1] = 0;
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
			ics_enhance = 1;
		} else {
			EGS4.seqStr = " No cross section enhancement";
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);

			ics_enhance = 0;
		}

		if (EGS4Macro.IFULL == 2 && ics_enhance == 1) {
			EGS4.seqStr = " You cannot calculate a pulse height distribution with";
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);
			EGS4.seqStr = " cross section enhancement.  Will run simulation with";
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);
			EGS4.seqStr = " IFULL= dose and stoppers";
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);
			EGS4Macro.IFULL = 0;
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
				s = " RAYLEIGH SCATTERING INCLUDED";
				EGS4.seqStr = s;// +"ON";
				if (EGS4.iprint > 1)
					printSequence(EGS4.seqStr);

				break;// EXIT;
			}
		}

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
		} else {
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
			s = " Charged particle Russian Roulette";
			ll = s.length();
			ll = 61 - ll;
			s = s + EGS4.format("", ll);
			EGS4.seqStr = s + "OFF";
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);
		} else {
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

	}

	/**
	 * Print output summary
	 */
	private void OSUMRY() {
		PLOTPH();
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

	/**
	 * Called by OSUMRY
	 */
	public void PLOTPH() {
		String s = "";
		//int ll = 0;

		s = " " + EGS4.format(IHSTRY, 12) + " HISTORIES ANALYSED";
		EGS4.seqStr = s;// \n";
		if (EGS4.iprint > 0)
			printSequence(EGS4.seqStr);

		s = "*******************************************************************************************";
		String PHOTONs = "PHOTON";
		// if(EGS4SrcEns.monoindex==EGS4SrcEns.iMONOENERGETIC)
		// {
		if (EGS4SrcEns.iqin == -1) {
			PHOTONs = "ELECTRON";
		} else if (EGS4SrcEns.iqin == 1) {
			PHOTONs = "POSITRON";
		}
		EGS4.seqStr = s;
		if (EGS4.iprint > 0)
			printSequence(EGS4.seqStr);
		if (EGS4SrcEns.iqin == 0) {

		} else {
			s = " DETECTOR EFFICIENCY FOR INCIDENT " + PHOTONs + " (%)= "
					+ EGS4.format(eefficiency, 14, true) + "(+-"
					+ EGS4.format(eefficiency_error, 6, true) + "%)";

			EGS4.seqStr = s;
			if (EGS4.iprint > 0)
				printSequence(EGS4.seqStr);
		}

		s = "******************************************************************************************";
		EGS4.seqStr = s;
		if (EGS4.iprint > 0)
			printSequence(EGS4.seqStr);
		// }
	}

	/**
	 * Fix some input variables in special cases. Called by init.
	 */
	private void fixSRCOTO() {
		EGS4SrcEns.iqin = EGS4SrcEns.ipart;// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
		EGS4SrcEns.NHSTRY = 0;
		EGS4SrcEns.last_nhstry = 0;

		EGS4SrcEns.pi = Math.PI;
	}

	/**
	 * Fix some source and energy related variables. Called by init. Position, direction 
	 * cosines, particle weight and incident fluence are all fixed here.
	 */
	private void fixEmAll() {
		//double D = 0.0;
		//double R2 = 0.0;
		//int IXIN = 0;
		// NHSTRY
		EGS4SrcEns.NHSTRY = EGS4SrcEns.NHSTRY + 1;// @@@@@@@@@@@@@@@@@@@@@@@@@

		if (SOURCE == SOURCE_POINT) {
			EGS4SrcEns.DISTZ = ps_distance;
			// ########FIX ainflu
			EGS4SrcEns.AINFLU = EGS4SrcEns.dble(EGS4SrcEns.NCASET)
					/ (4 * Math.PI * EGS4SrcEns.DISTZ * EGS4SrcEns.DISTZ);

			double diam = Math.sqrt(XDET * XDET + YDET * YDET) / 2.0;// diag==max

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

			double x0 = 0.0;
			double y0 = 0.0;
			EGS4SrcEns.xin = x0;
			EGS4SrcEns.yin = y0;
			EGS4SrcEns.irin = 2;
			EGS4SrcEns.zin = 0.0;

			EGS4SrcEns.uin = u;
			EGS4SrcEns.vin = v;
			EGS4SrcEns.win = w;
			// ##################fix weight
			double us1 = 2
					* Math.PI
					* (1 - EGS4SrcEns.DISTZ
							/ Math.sqrt(EGS4SrcEns.DISTZ * EGS4SrcEns.DISTZ
									+ diam * diam));
			EGS4SrcEns.WEIGHT = us1 / (4 * Math.PI);// FIXED
													// WEIGHT!!!@@@@@@@@@@@@@@
		} else if (SOURCE == SOURCE_SARPAGAN) {
			getCylinderRandom();
		}
	}

	/**
	 * Improper name of method. It handles BOX geometry and it is called by fixEmAll.
	 */
	private static void getCylinderRandom() {
		// BOX with XBOUND,YBOUND,ZBOUND and SC in the MASS CENTER!!
		double r = EGS4.random01();
		if (EGS4.random01() < 0.5)
			r = 1.0 - r;
		double z1 = (2. * r - 1.) * ZBOUND / 2.;
		r = EGS4.random01();
		if (EGS4.random01() < 0.5)
			r = 1.0 - r;
		double y1 = (2. * r - 1.) * YBOUND / 2.;
		r = EGS4.random01();
		if (EGS4.random01() < 0.5)
			r = 1.0 - r;
		double x1 = (2. * r - 1.) * XBOUND / 2.;

		EGS4SrcEns.DISTZ = z1;// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
		// EGS4SrcEns.AINFLU=EGS4SrcEns.dble(EGS4SrcEns.NCASET)/(4*Math.PI*EGS4SrcEns.DISTZ*EGS4SrcEns.DISTZ);//@@@@@@
		EGS4SrcEns.AINFLU = EGS4SrcEns.AINFLU + 1.0
				/ (4.0 * Math.PI * EGS4SrcEns.DISTZ * EGS4SrcEns.DISTZ);
		// for MC variance reduction we impose that all particles emerged from
		// point 0 will go into
		// an solid angle us1 having a distance d and a radius=2*detradius or
		// 2*sourceradius in order
		// to make sure that this allways contain the detector!!
		// here,detradius or sourceradius=sqrt(x*x+y*y)diag
		double diam = Math.sqrt(XBOUND * XBOUND + YBOUND * YBOUND);// new
																	// radius!!

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
		// -----------------------------------
		//double ux = u;
		//double uy = v;
		//double uz = w;

		EGS4SrcEns.xin = x1;
		EGS4SrcEns.yin = y1;
		EGS4SrcEns.zin = z1;
		EGS4SrcEns.uin = u;
		EGS4SrcEns.vin = v;
		EGS4SrcEns.win = w;

		EGS4SrcEns.irin = 2;// inside!!
	}

}
