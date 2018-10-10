package danfulea.phys.egsOutput;

import java.io.FileInputStream;
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
import danfulea.phys.egs.Graph;

/**
 * Class for XRAy spectra simulations. <p>
 * It calculates the emerging photon fluence distribution from material of interest. 
 * It is a rough approximation of real X-ray spectra due to limitation of RZ geometry. 
 * For simplicity, it is assumed a 4PI photon angle distribution when they are ejected from material. 
 * 
 * @author Dan Fulea 25 MAR. 2007
 */
public class BremsApp2 implements EgsQuestion {
	public static JTextArea jta;
	public static boolean systemOutB = true;

	public static int nREGSOLV = 0;
	public static int[] REGSVOL = new int[1000];// =0;
	public static double beam_angle = 0.0;// degree
	public static double beam_radius = 0.0;
	public static double[][] VOLUME;
	// @TRANSPORT
	// PARAMS:=============================================================
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
	public static int $MAXBRSPLIT = 200;// "MAX BREM SPLITTING NUMBER"
	// =============================================================================
	public static int $MXENE = 301;// REPLACE {$MXENE} WITH {200};
	public static int ncase = 0;
	public static int calc_type = 0;
	public static double e_tot = 0.0;
	public static double e_tot2 = 0.0;
	public static double e_brem = 0.0;
	public static double e_brem2 = 0.0;
	public static double e_bremc = 0.0;
	public static double e_rad = 0.0;
	public static double e_rad2 = 0.0;
	public static double e_radc = 0.0;
	public static double etot_tmp = 0.0;
	public static double ebrem_tmp = 0.0;
	public static double erad_tmp = 0.0;
	public static double my_gle = 0.0;
	public static double accu = 0.0;

	public static int neis = 0;
	public static double[] eis;// =new double[$MXENE];
	// ==========================================================
	public static int nbatch = 0;
	public static int nperbatch = 0;
	public static int ibatch = 0;
	public static int itimes = 0;
	public static int ntimes = 0;
	public static double aux = 0.0;
	public static double aux2 = 0.0;
	public static double total = 0.0;
	public static double anorm = 0.0;
	public static double sum = 0.0;
	public static double sum2 = 0.0;
	public static double xtest = 0.0;
	public static double xtest2 = 0.0;
	public static int ncasei = 0;
	public static int npgi = 0;
	public static int npei = 0;
	public static int nspliti = 0;
	public static int nbini = 0;
	public static double de_pulsei = 0.0;
	public static double Eave = 0.0;
	public static double err_Eave = 0.0;
	public static double factor = 0.0;

	public static double gmfp = 0.0;
	public static int lmy_gle = 0;
	public static double e_mutr = 0.0;
	public static double e_mutr2 = 0.0;
	public static double e_muen = 0.0;
	public static double e_muen2 = 0.0;
	public static double mutr = 0.0;
	public static double mutr2 = 0.0;
	public static double muen = 0.0;
	public static double muen2 = 0.0;
	public static double ert = 0.0;
	public static double ert2 = 0.0;
	public static double ett = 0.0;
	public static double ett2 = 0.0;

	// "For interacting with the source routine and shower"
	public static int iqin = 0;
	public static int irin = 0;
	public static int icase = 0;
	public static int ip = 0;
	public static double ein = 0.0;
	public static double uin = 0.0;
	public static double vin = 0.0;
	public static double win = 0.0;
	public static double wtin = 0.0;
	public static double xin = 0.0;
	public static double yin = 0.0;
	public static double zin = 0.0;
	public static double ecut_ask = 0.0;
	public static double pcut_ask = 0.0;
	public static double gbr1 = 0.0;
	public static double gbr2 = 0.0;
	public static double rnno = 0.0;
	public static int nmutr = 0;
	public static int nmuen = 0;
	// ============================================
	public static int source_type = 0;
	public static double angle = 0.0;
	public static double rbeam = 0.0;
	public static double distance = 0.0;
	public static double emax = 0.0;
	public static int mono = 0;

	public static String enerFilename = "";
	private static final String datas = "Data";
	private static final String egsDatas = "egs";
	private static final String dataspec = "spectra";
	private static String file_sep = System.getProperty("file.separator");
	private static String defaultext = ".spectrum";
	public static int $NENSRC = 300;// "MAX # OF POINTS IN ENERGY DISTRIBUTION  "
	public static int NENSRC = 0;
	public static double[] ENSRCD = new double[$NENSRC + 1];// (0:$NENSRC)
	public static double[] SRCPDF = new double[$NENSRC];// ($NENSRC)
	public static double[] srcpdf_at = new double[$NENSRC];// ($NENSRC),
	public static int[] srcbin_at = new int[$NENSRC];// ($NENSRC)
	public static double enmin = 0.0;
	public static double enmax = 0.0;
	public static int mode = 0;
	public static double eik = 0.0;
	public static double esum = 0.0;
	public static double esum2 = 0.0;
	public static double ecount = 0.0;
	public static double es = 0.0;
	public static double des = 0.0;

	public static int $EBIN = 500;
	public static double[] SCPDST = new double[$EBIN];
	public static double[] SCPDST2 = new double[$EBIN];
	public static double[] SCPCUM = new double[$EBIN];
	public static double[] SCPCUM2 = new double[$EBIN];
	public static double[] ECUM = new double[$EBIN];
	public static double[] ECUM2 = new double[$EBIN];
	public static double SCPTOT = 0.0;
	public static double SCPTOT2 = 0.0;
	public static double PHENER = 0.0;// brems photo energy
	public static double SLOTE = 0.0005;// delta E=>0.5keV
	public static int MAXBIN = 0;
	public static double EPHTOP = 0.100;// EGS4SrcEns.ein; MeV=====100kev
	public static int IB = 0;
	public static double[] BINTOP = new double[$EBIN];

	/**
	 * Constructor.
	 */
	public BremsApp2() {
		init();

	}

	/**
	 * Reset global variables for re-use.
	 */
	public static void reset() {
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

		ncase = 0;
		calc_type = 0;
		e_tot = 0.0;
		e_tot2 = 0.0;
		e_brem = 0.0;
		e_brem2 = 0.0;
		e_bremc = 0.0;
		e_rad = 0.0;
		e_rad2 = 0.0;
		e_radc = 0.0;
		etot_tmp = 0.0;
		ebrem_tmp = 0.0;
		erad_tmp = 0.0;
		my_gle = 0.0;
		accu = 0.0;
		neis = 0;

		SCPDST = new double[$EBIN];
		SCPDST2 = new double[$EBIN];
		SCPCUM = new double[$EBIN];
		SCPCUM2 = new double[$EBIN];
		ECUM = new double[$EBIN];
		ECUM2 = new double[$EBIN];
		REGSVOL = new int[1000];// =0;
		eis = new double[$MXENE];// mutrs=new double[$MXENE];
		// muens=new double[$MXENE];gs=new double[$MXENE];

		nbatch = 0;
		nperbatch = 0;
		ibatch = 0;
		itimes = 0;
		ntimes = 0;
		aux = 0.0;
		aux2 = 0.0;
		total = 0.0;
		anorm = 0.0;
		sum = 0.0;
		sum2 = 0.0;
		xtest = 0.0;
		xtest2 = 0.0;
		ncasei = 0;
		npgi = 0;
		npei = 0;
		nspliti = 0;
		nbini = 0;
		de_pulsei = 0.0;
		Eave = 0.0;
		err_Eave = 0.0;
		factor = 0.0;

		gmfp = 0.0;
		lmy_gle = 0;
		e_mutr = 0.0;
		e_mutr2 = 0.0;
		e_muen = 0.0;
		e_muen2 = 0.0;
		mutr = 0.0;
		mutr2 = 0.0;
		muen = 0.0;
		muen2 = 0.0;
		ert = 0.0;
		ert2 = 0.0;
		ett = 0.0;
		ett2 = 0.0;
		iqin = 0;
		irin = 0;
		icase = 0;
		ip = 0;
		ein = 0.0;
		uin = 0.0;
		vin = 0.0;
		win = 0.0;
		wtin = 0.0;
		xin = 0.0;
		yin = 0.0;
		zin = 0.0;
		ecut_ask = 0.0;
		pcut_ask = 0.0;
		gbr1 = 0.0;
		gbr2 = 0.0;
		rnno = 0.0;
		nmutr = 0;
		nmuen = 0;

		source_type = 0;
		angle = 0.0;
		rbeam = 0.0;
		distance = 0.0;
		emax = 0.0;
		mono = 0;
		enerFilename = "";
		NENSRC = 0;
		ENSRCD = new double[$NENSRC + 1];
		SRCPDF = new double[$NENSRC];
		srcpdf_at = new double[$NENSRC];
		srcbin_at = new int[$NENSRC];
		enmin = 0.0;
		enmax = 0.0;
		mode = 0;
		eik = 0.0;
		esum = 0.0;
		esum2 = 0.0;
		ecount = 0.0;
		es = 0.0;
		des = 0.0;
	}

	/**
	 * Perform basic initialization and RUN the Monte Carlo engine.
	 */
	private void init() {
		EGS4.startSimulationTime = System.currentTimeMillis();

		nREGSOLV = GammaDetEffA.nREGSOLV;
		for (int i = 1; i <= nREGSOLV; i++)// at reset GammaDetEffA, then
											// construct in main
		{
			REGSVOL[i - 1] = GammaDetEffA.REGSVOL[i - 1];
		}
		VOLUME = new double[EGS4Geom.$MAXZREG][EGS4Geom.$MAXRADII];
		EGS4Grid.CDSTBL = new String[EGS4.$MXREG];
		EGS4Grid.CTRTBL = new String[EGS4.$MXREG];
		EGS4Grid.CABSRB = new String[EGS4.$MXREG];
		// --variable init
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
		// ----------------------------------------
		// --interface->LINK
		EGS4.eq = this;// pass the printing mode
		EGS4Core.eq = this;// pass the printing mode
		EGS4Macro.eq = this;// pass the printing mode
		EGS4SrcEns.eq = this;// pass the printing mode
		EGS4Geom.eq = this;// pass the printing mode
		EGS4Grid.eq = this;// pass the printing mode
		// --Macro and defaults param:
		// EGS4.ispmfp=EGS4.iFlur;//EGS4.iGe;//select photon mfp
		EGS4.iGeom = 1;// 1=RZ GEOM@@@@@@@@@@@@@@YES
		EGS4.iprint = 2;// summary
		// =========================================inputs
		// Calendar cal = Calendar.getInstance();
		// Date d = cal.getTime();
		EGS4.seqStr = "=================================================================================";
		EGS4.seqStr = " ********************brems App: ************************************";
		if (EGS4.iprint > 1)
			printSequence(EGS4.seqStr);
		EGS4.seqStr = " It calculates the emerging photon fluence distribution from material of interest!";
		if (EGS4.iprint > 1)
			printSequence(EGS4.seqStr);
		EGS4.seqStr = " It is a rough approximation of real X-ray spectra due to limitation of RZ geometry!";
		if (EGS4.iprint > 1)
			printSequence(EGS4.seqStr);
		EGS4.seqStr = " For simplicity, it is assumed a 4PI photon angle distribution when they are ejected from material!";
		if (EGS4.iprint > 1)
			printSequence(EGS4.seqStr);

		EGS4.seqStr = " ********************************************************************************";
		if (EGS4.iprint > 1)
			printSequence(EGS4.seqStr);

		Double dbb = new Double(1.05 * EPHTOP / SLOTE + 1.0);
		MAXBIN = dbb.intValue();
		System.out.println("" + MAXBIN);

		inputs();
		if (EGS4.STOPPROGRAM) {
			return;
		}
		esum = 0.0;
		esum2 = 0.0;
		ecount = 0.0;

		fixSRCOTO();

		ecut_ask = EGS4.ECUT[0];
		pcut_ask = EGS4.PCUT[0];

		EGS4.seqStr = "  HATCHING media to get cross-section data...";
		if (EGS4.iprint > 1)
			printSequence(EGS4.seqStr);

		EGS4.iprint = 0;
		HATCH();
		if (EGS4.STOPPROGRAM) {
			closeFile();
			return;
		}
		EGS4.iprint = 2;

		if (ecut_ask < EGS4.AE[0] || pcut_ask < EGS4.AP[0]) {
			EGS4.seqStr = "  There is a mismatch between asked for and available cut-offs  ";
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);
			EGS4.seqStr = "  Asked for ECUT of " + EGS4.format(ecut_ask, 10)
					+ " MeV and have AE of " + EGS4.format(EGS4.AE[0], 10)
					+ " MeV";
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);
			EGS4.seqStr = "  Asked for PCUT of " + EGS4.format(pcut_ask, 10)
					+ " MeV and have AP of " + EGS4.format(EGS4.AP[0], 10)
					+ " MeV";
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);
			EGS4.seqStr = "  Exiting until there is a match.";
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);

			return;// STOP;
		}

		nbatch = 10;
		nperbatch = ncase / nbatch;
		if (nperbatch == 0)
			nperbatch = 1;
		// source_sumry();
		factor = 160.2176462; // "e_mutr and e_muen scored as MeV cm^2/g"
								// "*1.602176462E-13 J/MeV *1000 g/kg =160.2176E-12 Gy cm^2"
								// "Note, this is used elsewhere too"
		// boolean theEnd = false;
		if (neis > 1) {
			ntimes = neis;
			// write(6,'(//A)') '=====================================';
			EGS4.seqStr = "=====================================";
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);
			// write(6,'(A,I4,A)') 'Will loop through ',ntimes,' energies';
			EGS4.seqStr = "Will loop through " + ntimes + " energies";
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);
			// write(6,'(A//)') '=====================================';
			EGS4.seqStr = "=====================================";
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);
		} else {
			ntimes = 1;
		}
		// System.out.println("NZ "+EGS4Geom.NZ+" NR "+EGS4Geom.NR);
		for (int IZ = 1; IZ <= EGS4Geom.NZ; IZ++)// DO IZ=1,NZ
		{
			double DEPTH = EGS4Geom.ZPLANE[IZ] - EGS4Geom.ZPLANE[IZ - 1];// DEPTH=ZPLANE(IZ+1)-ZPLANE(IZ);
			for (int IX = 1; IX <= EGS4Geom.NR; IX++)// DO IX=1,NR [
			{
				double RLOW2 = 0.0;
				if (IX == 1) {
					RLOW2 = 0.0;
				} else {
					RLOW2 = EGS4Geom.CYRAD2[IX - 2];// CYRAD2(IX-1);
				}
				// VOLUME(IZ,IX)=Math.PI*DEPTH*(CYRAD2(IX)-RLOW2);
				VOLUME[IZ - 1][IX - 1] = Math.PI * DEPTH
						* (EGS4Geom.CYRAD2[IX - 1] - RLOW2);
				// System.out.println("depth "+DEPTH+" RLOW "+RLOW2+"  cyrad "+EGS4Geom.CYRAD2[IX-1]+" vol "+VOLUME[IZ-1][IX-1]);
			}
		}

		source_sumry();

		for (itimes = 1; itimes <= ntimes; itimes++) {
			e_mutr = 0.0;
			e_mutr2 = 0.0;
			e_muen = 0.0;
			e_muen2 = 0.0;
			mutr = 0.0;
			mutr2 = 0.0;
			muen = 0.0;
			muen2 = 0.0;
			nmutr = 0;
			nmuen = 0;

			source_switch_energy(itimes);

			for (icase = 1; icase < ncase; icase++) {
				source_sample();
				fixEmAll();
				// "============================================================================"
				SHOWER();
				// "============================================================================"
				ibatch = icase / nperbatch;
				if (ibatch * nperbatch == icase) {
					// write(6,'(a,i2,a,i2,a,f9.2,a)')
					// '+ Finished part ',ibatch,' out of ',nbatch,
					// ', cpu time = ',etime(time_array)-cpu,' sec.';
					// String timePerBatch=EGS4.timeElapsedShort(startTime);
					Calendar call = Calendar.getInstance();
					String timeday = call.get(Calendar.HOUR) + ":"
							+ call.get(Calendar.MINUTE) + ":"
							+ call.get(Calendar.SECOND);

					EGS4.seqStr = " Finished part " + ibatch + "  out of  "
							+ nbatch + " cpu time: " + EGS4.format(timeday, 9);// +//+EGS4.format(100*e_rad2*(e_tot/e_rad-1),12)+" %";
					if (EGS4.iprint > 1)
						printSequence(EGS4.seqStr);

				}

			}// for
				// "-----------------------------------------------------------------"
				// "STEP 8   OUTPUT-OF-RESULTS                                       "
				// "-----------------------------------------------------------------"

			/*
			 * if(PHENER>0.) { //"EQUAL ENERGY BINS CASE" Double dbll=new
			 * Double(PHENER/SLOTE+1.0);//NOT ZEROES
			 * //IB=Math.min(IFIX(PHENER/SLOTE+0.999),$EBIN);
			 * IB=Math.min(dbll.intValue(),$EBIN);
			 * 
			 * //"ACCUMULATE THE PULSE HEIGHT DISTRIBUTION"
			 * SCPDST[IB-1]=SCPDST[IB-1]+1.0;//WT1OLD;
			 * SCPDST2[IB-1]=SCPDST2[IB-1]+1.0;//WT1OLD*WT1OLD;
			 * //"also add to cumulative distn" for(int
			 * ICUM=IB;ICUM<=MAXBIN;ICUM++) {
			 * SCPCUM[ICUM-1]=SCPCUM[ICUM-1]+1.0;//WT1OLD;
			 * SCPCUM2[ICUM-1]=SCPCUM2[ICUM-1]+1.0;//WT1OLD*WT1OLD;
			 * 
			 * }
			 * 
			 * PHENER=0.0;//reset }
			 */

			// EGS4.seqStr=" Average spectrum energy: "+EGS4.format(Eave,12);//+
			// " cm2/g";
			// if(EGS4.iprint>1)
			// printSequence(EGS4.seqStr);

		}// ]"end of itimes loop"

		int SCORE_NORM_NUM = ncase;
		double ech = 1.60217646E-19;// C
		double totfluence = 0.0;
		for (IB = 1; IB <= MAXBIN; IB++) {
			SCPTOT = SCPTOT + SCPDST[IB - 1];
			SCPTOT2 = SCPTOT2 + SCPDST2[IB - 1];

			// double SCORE_TEMP2 = SCPDST2[IB - 1];
			// $ANALYZE(SCPDST,(IB):SCORE_NORM_NUM);
			double SCORE_TEMP = SCPDST[IB - 1] / SCORE_NORM_NUM;
			SCPDST2[IB - 1] = SCPDST2[IB - 1] / SCORE_NORM_NUM;
			SCPDST2[IB - 1] = (SCPDST2[IB - 1] - SCORE_TEMP * SCORE_TEMP)
					/ (SCORE_NORM_NUM - 1);
			if (SCPDST2[IB - 1] >= 0.)
				SCPDST2[IB - 1] = Math.sqrt(SCPDST2[IB - 1]);
			if (SCORE_TEMP != 0.) {
				SCPDST2[IB - 1] = Math.min(SCPDST2[IB - 1] / SCORE_TEMP * 100.,
						99.9);
			} else {
				SCPDST2[IB - 1] = 99.9;
			}
			// ======================
			SCPDST[IB - 1] = SCPDST[IB - 1] * 1.0
					/ (Math.PI * 4.0 * 750.0 * 750.0);// fluence [ph/mm2] at 750
														// mm
			SCPDST[IB - 1] = (SCPDST[IB - 1] / ncase) * 1.0E-03 / ech;
			// is x/n * N at 1mAs, N1mAs=10-3As/echarge
			// SCPDST=>photons/(mm2 x mAs) at 750 mm
			// =========================
			totfluence = totfluence + SCPDST[IB - 1];
		}

		double SCORE_TEMP = SCPTOT / SCORE_NORM_NUM;
		SCPTOT2 = SCPTOT2 / SCORE_NORM_NUM;
		SCPTOT2 = (SCPTOT2 - SCORE_TEMP * SCORE_TEMP) / (SCORE_NORM_NUM - 1);
		if (SCPTOT2 >= 0.)
			SCPTOT2 = Math.sqrt(SCPTOT2);
		if (SCORE_TEMP != 0.) {
			SCPTOT2 = Math.min(SCPTOT2 / SCORE_TEMP * 100., 99.9);
		} else {
			SCPTOT2 = 99.9;
		}

		PLOTPH();

		EGS4.seqStr = "total emerging photon fluence [photons/(mm2*mAs) at 750mm]: "
				+ totfluence + "  +/- " + SCPTOT2 + "  %";
		if (EGS4.iprint > 0)
			printSequence(EGS4.seqStr);

		int ndat = MAXBIN;
		Graph.yS = "photon fluence/mAs at 750mm";
		Graph.titleS = "Photons spectrum (Brems application)";

		Graph.interv2 = ndat;
		Graph.xarray2 = new double[ndat];
		Graph.yarray2 = new double[ndat];
		for (int i = 1; i <= ndat; i++) {
			Graph.xarray2[i - 1] = 1000.0 * i * SLOTE;// i;//KeV
			Graph.yarray2[i - 1] = SCPDST[i - 1];
		}
		new Graph();

		// =======================================================================
		EGS4.timeElapsed();
		Calendar call = Calendar.getInstance();
		Date da = call.getTime();
		EGS4.seqStr = "End of run:          " + da.toString();
		if (EGS4.iprint > 0)
			printSequence(EGS4.seqStr);

		// ===============================================================
	}

	/**
	 * Display input data summary.
	 */
	private void source_sumry() {
		String s1 = "";
		String s = "";
		int ll = 0;

		EGS4.seqStr = "=========================================================================";
		if (EGS4.iprint > 1)
			printSequence(EGS4.seqStr);
		EGS4.seqStr = "                   MONTE CARLO, ENERGY AND SOURCE";
		if (EGS4.iprint > 1)
			printSequence(EGS4.seqStr);
		EGS4.seqStr = "=========================================================================";
		if (EGS4.iprint > 1)
			printSequence(EGS4.seqStr);
		s = " Max. # of histories to RUN";
		ll = s.length();
		ll = 54 - ll;
		s = s + EGS4.format("", ll);
		EGS4.seqStr = s + EGS4.format(ncase, 12);
		if (EGS4.iprint > 1)
			printSequence(EGS4.seqStr);

		s = " Incident Charge";// EGS4.format("",20)+"Incident Charge";
		s1 = "";
		ll = s.length();
		ll = 68 - ll;
		if (iqin == 0)
			s1 = "photons";
		if (iqin == -1)
			s1 = "electrons";
		if (iqin == 1)
			s1 = "positrons";
		s = s + EGS4.format(s1, ll);
		EGS4.seqStr = s;
		if (EGS4.iprint > 1)
			printSequence(EGS4.seqStr);

		if (mono == 0 || mono == 2) {
			EGS4.seqStr = " Incident spectrum:"
					+ EGS4.format("MONO-ENERGY", ll);
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);

			for (int i = 1; i <= neis; i++) {
				s = " Incident kinetic energy:";
				ll = s.length();
				ll = 57 - ll;
				s = s + EGS4.format("", ll);

				EGS4.seqStr = s
						+ EGS4.format(eis[i - 1] - EGS4.RM * Math.abs(iqin), 9,
								true) + " MeV";
				if (EGS4.iprint > 1)
					printSequence(EGS4.seqStr);
			}

		} else {
			EGS4.seqStr = " Incident spectrum:" + EGS4.format("SPECTRUM", ll);
			;// "SPECTRUM";
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);

			s = " Average spectrum energy       : ";
			ll = s.length();
			ll = 57 - ll;
			s = s + EGS4.format("", ll);
			EGS4.seqStr = s + EGS4.format(Eave, 9, true) + " MeV";// ensrcd(nensrc)
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);
			s = " Maximum spectrum energy       : ";
			ll = s.length();
			ll = 57 - ll;
			s = s + EGS4.format("", ll);
			EGS4.seqStr = s + EGS4.format(enmax, 9, true) + " MeV";// ensrcd(nensrc)
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);

		}

		/*
		 * s=" Source type                  : ";
		 * ll=s.length();ll=57-ll;s=s+EGS4.format("",ll);
		 * EGS4.seqStr=s+EGS4.format(source_type,9);//ensrcd(nensrc)
		 * if(EGS4.iprint>1) printSequence(EGS4.seqStr);
		 * 
		 * if (source_type==0) { s=" Incident angle               : ";
		 * ll=s.length();ll=57-ll;s=s+EGS4.format("",ll);
		 * EGS4.seqStr=s+EGS4.format(angle,9,true)+" degrees";//ensrcd(nensrc)
		 * if(EGS4.iprint>1) printSequence(EGS4.seqStr);
		 * 
		 * } else if (source_type==1) { s=" Beam size on front face      : ";
		 * ll=s.length();ll=57-ll;s=s+EGS4.format("",ll);
		 * EGS4.seqStr=s+EGS4.format(rbeam,9,true)+" cm";//ensrcd(nensrc)
		 * if(EGS4.iprint>1) printSequence(EGS4.seqStr);
		 * s=" Source-face distance         : ";
		 * ll=s.length();ll=57-ll;s=s+EGS4.format("",ll);
		 * EGS4.seqStr=s+EGS4.format(distance,9,true)+" cm";//ensrcd(nensrc)
		 * if(EGS4.iprint>1) printSequence(EGS4.seqStr);
		 * 
		 * }
		 */
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
		EGS4.seqStr = "=========================================================================";
		if (EGS4.iprint > 1)
			printSequence(EGS4.seqStr);
		EGS4.seqStr = "                   Variance reduction and electron thresholds have no effect here!!";
		if (EGS4.iprint > 1)
			printSequence(EGS4.seqStr);
		EGS4.seqStr = "                   Valid top energy is 150keV!";
		if (EGS4.iprint > 1)
			printSequence(EGS4.seqStr);

		EGS4.seqStr = "=========================================================================";
		if (EGS4.iprint > 1)
			printSequence(EGS4.seqStr);

		EGS4Grid.CABSRB[0] = "0";
		EGS4Grid.CDSTBL[0] = "0";
		EGS4Grid.CTRTBL[0] = "0";
		EGS4Grid.MATERIALGRID(EGS4Geom.NR, EGS4Geom.NZ, VOLUME, 2, EGS4.ECUT,
				EGS4.PCUT, EGS4Geom.RCYL, EGS4Geom.ZPLANE, EGS4.MED, EGS4.MEDIA);

	}

	/**
	 * Setup input variables.
	 */
	public void inputs() {
		EGS4.seqStr = " *** Reading geometrical info ... ***";
		if (EGS4.iprint > 1)
			printSequence(EGS4.seqStr);

		EGS4.iprint = 0;
		EGS4Geom.GEOMRZ();
		if (EGS4.STOPPROGRAM) {
			closeFile();
			return;
		}
		EGS4.iprint = 2;
		// call get_transport_parameter(6);
		// "read transport parameter settings and"
		setEcutRegion = false;
		setPcutRegion = false;
		setSmaxRegion = false;
		ecut = 0.521;
		// pcut=0.001;
		pcut = 0.01;
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
		// Pair angular sampling= On->IPRDST
		EGS4.estepe = 0.25;
		// XIMAX= 0.5->XIMAX
		EGS4.ximax = 0.5;
		// Boundary crossing algorithm= exact->bca_algorithm, exact_bca
		// Skin depth for BCA= 3->skindepth_for_bca
		EGS4.skindepth_for_bca = 3.0;
		// Electron-step algorithm= default->transport_algorithm
		// Spin effects= On->spin_effects
		// Brems cross sections= BH
		// " Pair cross sections "= BH
		// Electron impact ionization
		// Triplet
		// Radiative compton correction

		// " Media input "
		// ====OUTSIDE==================
		// ==============================
		EGS4.DUNIT = 1; // "SET LENGTH UNITS TO CMS"

		if (neis > $MXENE) {
			// too many energies requested, increase MXENE
			EGS4.STOPPROGRAM = true;
			EGS4.seqStr = "too many energies requested, max. allowed= "
					+ $MXENE;
			printSequence(EGS4.seqStr);
		}

		energy_inputs();
		source_init();

		for (int j = 1; j <= 25; j++) {
			EGS4.iausfl[j - 1] = 0;
		}

		// EGS4.iausfl[7] = 1; //"after brems"
		// EGS4.iausfl[25] = 1;// "after Xcharacteritic"

		test_inputs();
		if (EGS4.STOPPROGRAM) {
			closeFile();
			return;
		}
	}

	/**
	 * Validate more inputs. Called by inputs routine.
	 */
	private void test_inputs() {
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
			EGS4.prob_RR = 1. / EGS4.nbr_split;
		} else {
			EGS4.prob_RR = 1.;
		}

	}

	/**
	 * Dummy routine
	 */
	public void closeFile() {
	}

	/**
	 * Where to print runtime information. Interface method.
	 * @param s the String to be printed
	 */
	public void printSequence(String s) {
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
		// DO NOT NEED ANY CALL OF AUSGAB!!!anything is scored in HOWFAR in
		// checkBrems
		/*
		 * if( IARG == 7 ) {//
		 * "after a call to brems - add up the energy of the photon" //
		 * "and then discard the photon" if( EGS4.IQ[EGS4.NP-1] == 0 )//if(
		 * iq(np) = 0 ) {//"photon is top of stack" //ebrem_tmp = ebrem_tmp +
		 * EGS4.E[EGS4.NP-1];//e(np); //erad_tmp = erad_tmp +
		 * EGS4.E[EGS4.NP-1];//e(np); //IF( $DEBUGIT ) [ write(18,*) ' brems:
		 * ',e(np),wt(np); ] //wt(np) = 0; e(np) = 0 //EGS4.WT[EGS4.NP-1] = 0.0;
		 * //EGS4.E[EGS4.NP-1] = 0.0; PHENER=EGS4.E[EGS4.NP-1];//score photon!!
		 * if (PHENER>0.0) { Double dbll=new Double(PHENER/SLOTE+1.0);//NOT
		 * ZEROES IB=Math.min(dbll.intValue(),$EBIN);
		 * SCPDST[IB-1]=SCPDST[IB-1]+1.0;//WT1OLD;
		 * SCPDST2[IB-1]=SCPDST2[IB-1]+1.0;//WT1OLD*WT1OLD;
		 * //"also add this to the cumulative pulse height distn" for(int
		 * ICUM=IB;ICUM<=MAXBIN;ICUM++) {
		 * SCPCUM[ICUM-1]=SCPCUM[ICUM-1]+1.0;//WT1OLD;
		 * SCPCUM2[ICUM-1]=SCPCUM2[ICUM-1]+1.0;//WT1OLD*WT1OLD; } }
		 * EGS4.WT[EGS4.NP-1] = 0.0; } else {//
		 * "electron was top of stack=> photon is np-1" //ebrem_tmp = ebrem_tmp
		 * + EGS4.E[EGS4.NP-2];//e(np-1); //erad_tmp = erad_tmp +
		 * EGS4.E[EGS4.NP-2];//e(np-1); //IF( $DEBUGIT ) [ write(18,*) ' brems:
		 * ',e(np-1),wt(np-1); ] //wt(np-1) = 0; e(np-1) = 0 EGS4.WT[EGS4.NP-2]
		 * = 0.0; EGS4.E[EGS4.NP-2] = 0.0;
		 * 
		 * PHENER=EGS4.E[EGS4.NP-2];//score photon!! if (PHENER>0.0) { Double
		 * dbll=new Double(PHENER/SLOTE+1.0);//NOT ZEROES
		 * IB=Math.min(dbll.intValue(),$EBIN);
		 * SCPDST[IB-1]=SCPDST[IB-1]+1.0;//WT1OLD;
		 * SCPDST2[IB-1]=SCPDST2[IB-1]+1.0;//WT1OLD*WT1OLD;
		 * //"also add this to the cumulative pulse height distn" for(int
		 * ICUM=IB;ICUM<=MAXBIN;ICUM++) {
		 * SCPCUM[ICUM-1]=SCPCUM[ICUM-1]+1.0;//WT1OLD;
		 * SCPCUM2[ICUM-1]=SCPCUM2[ICUM-1]+1.0;//WT1OLD*WT1OLD; } }
		 * EGS4.WT[EGS4.NP-1] = 0.0;///terminate history, no more transport }
		 * return; } if( IARG == 25 )//fluoro xray {
		 * PHENER=EGS4.E[EGS4.NP-1];//score photon!! if (PHENER>0.0) { Double
		 * dbll=new Double(PHENER/SLOTE+1.0);//NOT ZEROES
		 * IB=Math.min(dbll.intValue(),$EBIN);
		 * SCPDST[IB-1]=SCPDST[IB-1]+1.0;//WT1OLD;
		 * SCPDST2[IB-1]=SCPDST2[IB-1]+1.0;//WT1OLD*WT1OLD;
		 * //"also add this to the cumulative pulse height distn" for(int
		 * ICUM=IB;ICUM<=MAXBIN;ICUM++) {
		 * SCPCUM[ICUM-1]=SCPCUM[ICUM-1]+1.0;//WT1OLD;
		 * SCPCUM2[ICUM-1]=SCPCUM2[ICUM-1]+1.0;//WT1OLD*WT1OLD; } }
		 * EGS4.WT[EGS4.NP-1] = 0.0; } return;
		 */
	}// "END OF AUSGAB"

	/**
	 * The following is a general specification of HOWNEAR: 
	 * Given a particle at (x,y,z) in region irl, HOWNEAR answers the 
	 * question, What is the distance tperp to the closest boundary? Interface method.
	 */
	public void HOWNEAR() {
		/* EGS4.tperp=1.0e10; */
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

	/**
	 * Replaces AUSGAB. Everything is scored here. Called by HOWFAR!
	 */
	public void checkBrems() {
		// System.out.println("enter");
		int IQL = EGS4.IQ[EGS4.NP - 1];// curent
		int IRL = EGS4.IR[EGS4.NP - 1];// "LOCAL REGION NUMBER"
		// EGS4.IRNEW=>changed by HOWFAR!!
		int imoi = REGSVOL[0];// where to score IRL of interest!!
		if (imoi != IRL || IQL != 0 || EGS4.IRNEW == IRL)
			return;
		// study only for photons emerging from imoi at any angle=>4PI
		// distribution!!!!!
		// ANY PHOTONS HERE!!!
		double PHENER = EGS4.E[EGS4.NP - 1];
		int IB = 0;
		int MAXIB = MAXBIN;// !!!!!!!!!!
		if (PHENER > 0.0)// SLOTE always>0!!
		{
			Double dbll = new Double(PHENER / SLOTE + 1.0);// NOT ZEROES
			IB = Math.min(dbll.intValue(), MAXIB);// MAXIB=>IB=Math.min(dbll.intValue(),MAXIB);
			SCPDST[IB - 1] = SCPDST[IB - 1] + 1.0;// WT1OLD~for paralel beam
													// always 1!!;
			SCPDST2[IB - 1] = SCPDST2[IB - 1] + 1.0;// WT1OLD*WT1OLD;
			// "also add this to the cumulative pulse height distn"
			for (int ICUM = IB; ICUM <= MAXIB; ICUM++) {
				SCPCUM[ICUM - 1] = SCPCUM[ICUM - 1] + 1.0;// WT1OLD;
				SCPCUM2[ICUM - 1] = SCPCUM2[ICUM - 1] + 1.0;// WT1OLD*WT1OLD;
			}
		}

	}

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
		// System.out.println("enter");
		if (EGS4.WT[EGS4.NP - 1] == 0.0) {
			EGS4.IDISC = 1;
			return;
		}

		// "INITIALLY ASSUME PARTICLE STAYS IN THE TARGET"
		boolean OUTEND = false;
		boolean OUTSID = false;
		// int IQL = EGS4.IQ[EGS4.NP - 1];
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

		int IZNEW = IZ + EGS4Geom.ihitp;// IHITP; "GET NEW PLANAR REGION"
		if ((IZNEW < 1) || (IZNEW > EGS4Geom.NZ))
			OUTEND = true; // "FLAG IF LEAVES BY THE ENDS"

		int IXNEW = IX + EGS4Geom.ihitc;// IHITC; "GET NEW CYLINDRICAL REGION"
		if (IXNEW > EGS4Geom.NR)
			OUTSID = true; // "FLAG IF LEAVES BY THE SIDES"

		if ((EGS4Geom.ihitp != 0) && (EGS4Geom.ihitc != 0)) {
			if (EGS4Geom.tplane < EGS4Geom.tcyl) {// "HITS PLANE FIRST"
				if (OUTEND) {
					// NWNRCL=0;
					// $SET NEW REGION(TPLANE,1);
					// if({P1}.LE.USTEP)[
					if (EGS4Geom.tplane <= EGS4.USTEP) {
						EGS4.USTEP = EGS4Geom.tplane;// {P1};
						EGS4.IRNEW = 1;// {P2};
						checkBrems();
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
						checkBrems();
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
						checkBrems();
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
						checkBrems();
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
						checkBrems();
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
						checkBrems();
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
					checkBrems();
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
					checkBrems();
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
					checkBrems();
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
					checkBrems();
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
		// System.out.println("enter "+irin+"  ein "+ein+ " wtin "+wtin );
		// call shower(iqin,ein,xin,yin,zin,uin,vin,win,irin,wtin);
		EGS4Core.SHOWER(iqin, ein, xin, yin, zin, uin, vin, win, irin, wtin);
	}

	/**
	 * Setup energy inputs. Called by inputs routine.
	 */
	private void energy_inputs() {
		if (mono == 1) {
			readSpectrum(enerFilename);
			EGS4.seqStr = " Spectrum file:  " + enerFilename;
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);

			EGS4.seqStr = "    HAVE READ" + EGS4.format(NENSRC, 5)
					+ " INPUT ENERGY BINS FROM FILE";
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);

			if (mode == 0) {
				EGS4.seqStr = "      Counts/bin assumed";
				if (EGS4.iprint > 1)
					printSequence(EGS4.seqStr);
			} else if (mode == 1) {
				EGS4.seqStr = "      Counts/MeV assumed";
				if (EGS4.iprint > 1)
					printSequence(EGS4.seqStr);

				for (int IB = 1; IB <= NENSRC; IB++) {
					SRCPDF[IB - 1] = SRCPDF[IB - 1]
							* (ENSRCD[IB] - ENSRCD[IB - 1]);
				}
			}// "end mode = 1 block"
			else {
				EGS4.seqStr = "*****MODE not 0 or 1 in spectrum file? **";
				if (EGS4.iprint > 1)
					printSequence(EGS4.seqStr);
			}

			ein = ENSRCD[NENSRC];// "SET TO MAX ENERGY FOR SOME CHECKS"

			EGS4.seqStr = "    ENERGY RANGES FROM"
					+ EGS4.format(enmin, 10, true) + " MeV TO"
					+ EGS4.format(ein, 12, true) + " MeV";
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);

			enmax = ein;
			Eave = 0;
			sum = 0;
			for (int i = 1; i <= NENSRC; i++) {
				sum = sum + SRCPDF[i - 1];
				Eave = Eave + 0.5 * SRCPDF[i - 1] * (ENSRCD[i] + ENSRCD[i - 1]);
			}
			Eave = Eave / sum;
			// call prepare_alias_sampling(nensrc,srcpdf,srcpdf_at,srcbin_at);
			EGS4.prepare_alias_sampling(NENSRC, SRCPDF, srcpdf_at, srcbin_at);
		}
	}

	/**
	 * Plot some useful information.
	 */
	public void PLOTPH() {
		String s = "";
		// int ll = 0;
		// "GET PEAK IN SPECTRUM"
		double sfac = 0.0;
		double EB = 0.0;
		String SPACE = " ";// /,BAR/'|'/,SYMBOL/$S'*+$-#@'/;
		String BAR = "|";
		String[] SYMBOL = { "*", "+", "$", "-", "#", "@" };
		String[] LINE = new String[61];
		int NHIST = 61;
		double HIST = 61.;
		// int IOUT = 1;
		int ILEV = 0;

		for (int IB = 1; IB <= MAXBIN; IB++) {
			if (SCPDST[IB - 1] > sfac) {
				sfac = SCPDST[IB - 1];
			}
		}

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
				EB = 1000.0 * IB * SLOTE;
			} else {
				EB = BINTOP[IB - 1];
			}// kev
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
		s = " PLOT NORMALIZED TO PEAK OF ONE, "
				+ "PDST IS NORMALIZED TO UNIT AREA";
		EGS4.seqStr = s;// \n";
		if (EGS4.iprint > 0)
			printSequence(EGS4.seqStr);
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
	 * Dummy method.
	 */
	private void source_init() {
		// if (source_type==0)
		// {
		// win = angle/180*Math.PI; win = Math.cos(win);
		// uin = Math.sqrt((1-win)*(1+win)); vin = 0;
		// }
		// else if (source_type==1)
		// {
		// read rbeam, distance=>outside
		// }
	}

	/**
	 * Read spectrum from a file.
	 * @param enerfile the spectrum file
	 */
	private void readSpectrum(String enerfile) {
		String filename = datas + file_sep + egsDatas + file_sep + dataspec
				+ file_sep + enerfile + defaultext;
		int iread = 0;
		int lnr = 0;// data number
		int lnrr = 0;// line number
		int indx = 1;
		StringBuffer desc = new StringBuffer();
		boolean haveData = false;
		// String equals = "=";
		char comma = ',';
		char lineSep = '\n';// System.getProperty("line.separator").charAt(0);

		boolean enB = false;
		boolean srB = false;

		try {
			FileInputStream in = new FileInputStream(filename);

			while ((iread = in.read()) != -1) {
				if (!Character.isWhitespace((char) iread)
						&& ((char) iread != comma)) {
					desc.append((char) iread);
					haveData = true;
				} else {
					if (haveData)// we have data
					{
						// System.out.println("ok");
						haveData = false;// reset
						if (lnrr != 0)// for skip first line
						{
							lnr++;
							// READ(9,*) nensrc,ensrcd(0),mode;
							if (lnr == 1) {
								String s = desc.toString();
								NENSRC = EGS4.stringToInt(s);// System.out.println(NENSRC);
							}
							if (lnr == 2) {
								String s = desc.toString();
								ENSRCD[0] = EGS4.stringToDouble(s);// System.out.println(ENSRCD[0]);
							}
							if (lnr == 3) {
								String s = desc.toString();
								mode = EGS4.stringToInt(s);// System.out.println(mode);
								enmin = ENSRCD[0];
								if (NENSRC > $NENSRC) {
									// OUTPUT NENSRC,$NENSRC;
									// (//' ********** Asked for too many energy
									// bins******'/
									// ' NENSRC =',I4, ' reduced to max allowed
									// =',I4/1x,30('*')//);
									EGS4.seqStr = " ********** Asked for too many energy bins******";
									if (EGS4.iprint > 1)
										printSequence(EGS4.seqStr);
									EGS4.seqStr = " NENSRC =" + NENSRC
											+ " reduced to max allowed ="
											+ $NENSRC;
									// EGS4.STOPPROGRAM=true;
									if (EGS4.iprint > 1)
										printSequence(EGS4.seqStr);

									NENSRC = $NENSRC;
								}

							} else if (lnr > 3) {
								// READ(9,*)(ENSRCD(IB),SRCPDF(IB),IB=1,NENSRC);
								String s = desc.toString();
								if (!enB) {
									ENSRCD[indx] = EGS4.stringToDouble(s);
									enB = true;
								} else if (!srB) {
									SRCPDF[indx - 1] = EGS4.stringToDouble(s);
									// System.out.println(ENSRCD[indx]+"    "+SRCPDF[indx-1]);

									indx++;
									if (indx == NENSRC + 1)
										break;
									srB = true;
									enB = false;
									srB = false;
								}
							}

						}

					}// have data
					if ((char) iread == lineSep)
						lnrr++;
					desc = new StringBuffer();
				}// else
			}// main while
			in.close();
		}// try
		catch (Exception exc) {

		}

	}

	// private void source_sample(int iqin,int irin,double ein,double xin,double
	// yin,
	// double zin,double uin,double vin,double win,double wtin);
	/**
	 * Internally used. Called by init.
	 * 
	 */
	private void source_sample() {
		// "=============================================================="
		// iqin = iqi;
		// irin = 1;
		if (mono == 0 || mono == 2)// | mono = 3 )
		{
			// ein = ei;
			wtin = 1;
		} else {
			// ein = alias_sample(nensrc,ensrcd,srcpdf_at,srcbin_at);
			ein = EGS4.alias_sample(NENSRC, ENSRCD, srcpdf_at, srcbin_at);
			wtin = 1;
			eik = ein;
			if (iqin != 0) {
				ein = ein + EGS4.RM;
			}// total energy
		}
		ecount = ecount + 1;
		esum = esum + eik;
		esum2 = esum2 + eik * eik;

		// if( source_type == 0 )
		// {
		// xin = 0; uin = 0; zin = 0;
		// uin = ui; vin = vi; win = wi;
		// }
		// else
		// {
		// double r=EGS4.random01();//$RANDOMSET r;
		// r = rbeam*Math.sqrt(r);
		// double phi=EGS4.random01();//$RANDOMSET phi;
		// phi = 2*phi*Math.PI;
		// xin = r*Math.cos(phi); yin = r*Math.sin(phi);
		// aux = 1/Math.sqrt(xin*xin + yin*yin + distance*distance);
		// uin = xin*aux; vin = yin*aux; win = distance*aux;
		// zin = 0;
		// }

		return;
	}

	/**
	 * Internally used. Called by init.
	 * @param itimes itimes
	 */
	private void source_switch_energy(int itimes) {
		// System.out.println("enter "+itimes);
		// "========================="
		if (mono == 1) {
			return;
		}
		ein = eis[itimes - 1];
		eik = ein - EGS4.RM * Math.abs(iqin);// ei=eis(itimes); eik = ei;
		ecount = 0;
		esum = 0;
		esum2 = 0;
		if (neis > 1) {
			// write(6,'(a,$)')
			// '===> Incident kinetic energy (MeV): ';
			// write(6,'(f15.8/)') ei-abs(iqi)*rm;
			EGS4.seqStr = " ===> Incident kinetic energy (MeV): "
					+ EGS4.format(ein - Math.abs(iqin) * EGS4.RM, 12, true);
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);

		}
		// System.out.println("ein "+ein);
		return;
		// end; "of subroutine source"
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
		// if(SOURCE!=SOURCE_BEAM)
		EGS4SrcEns.IFPB = 1;// does not matter here!!no var red will be used!!!
		// else
		// EGS4SrcEns.IFPB=0;
		EGS4SrcEns.pi = Math.PI;
	}

	/**
	 * Fix some source and energy related variables. Called by init. Position, direction 
	 * cosines, particle weight and incident fluence are all fixed here.
	 */
	private void fixEmAll() {
		// double D = 0.0;
		double R2 = 0.0;
		int IXIN = 0;
		// NHSTRY
		// EGS4SrcEns.NHSTRY=EGS4SrcEns.NHSTRY+1;//@@@@@@@@@@@@@@@@@@@@@@@@@

		// if(SOURCE==SOURCE_BEAM)
		// {
		// public static double dist_source_det=0.0;
		// public static double beam_angle=0.0;//degree
		// public static double beam_radius=0.0;
		// Sarpagan,Marinelli and source Point=>always u,v, independent of
		// xin,yin due to
		// 4 PI random emision!!
		// Here, Weight=1 and xin,yin are closed related to u and v!!!
		// ########FIX ainflu
		// EGS4SrcEns.AINFLU=EGS4SrcEns.dble(EGS4SrcEns.NCASET)/(Math.PI*beam_radius*beam_radius);
		EGS4SrcEns.AINFLU = EGS4SrcEns.dble(ncase)
				/ (Math.PI * beam_radius * beam_radius);
		// ##################
		double teta = Math.abs(Math.PI * beam_angle / 180.0);// with z axis!!
		double w = Math.cos(teta);
		double u = 0.0;
		double v = 0.0;
		R2 = 0.0;// double R2=0.0;
		double rbeam_at_frontFace = beam_radius;// /Math.cos(teta);//>beam_radius!!
		// deplasarea fata de centrul cercului=identic cu deplasarea pe 'cercul'
		// de pe
		// marginea exterioara superioara a geometriei!! fascicul cu rbeam ff
		// mic!!
		// si angle de maxim 20!!=>rezultate corecte pt fascicule inguste
		// relativ la raza det!!!
		while (true) {// [
			xin = EGS4.random01();
			xin = (2.0 * xin - 1.0) * rbeam_at_frontFace;
			// $RANDOMSET YIN;YIN=(2.0*YIN-1.0)*RBEAM;
			yin = EGS4.random01();
			yin = (2.0 * yin - 1.0) * rbeam_at_frontFace;
			R2 = xin * xin + yin * yin;

			if (R2 <= rbeam_at_frontFace * rbeam_at_frontFace) {
				double cosphi2 = xin / R2;// 0.0;
				double sinphi2 = yin / R2;
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
		irin = 2 + (IXIN - 1) * EGS4Geom.NZ;// @@@@@@@@@@@@@@@@@@@@@@@@
		zin = EGS4Geom.ZPLANE[0];// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@

		uin = u;
		vin = v;
		win = w;
		// if(w!=1.0)EGS4SrcEns.IFPB=0;//!!!!!!!!!
		// EGS4SrcEns.NRCFLG=10;
		// ##################fix weight
		// EGS4SrcEns.WEIGHT=1.0;
		wtin = 1.0;
		// }
		// other sample not good here fascicul e- paralel!!!!!!:
		// double r=EGS4.random01();//$RANDOMSET r;
		// eventual daca distance=
		// r = rbeam*Math.sqrt(r);
		// double phi=EGS4.random01();//$RANDOMSET phi;
		// phi = 2*phi*Math.PI;
		// xin = r*Math.cos(phi); yin = r*Math.sin(phi);
		// aux = 1/Math.sqrt(xin*xin + yin*yin + distance*distance);
		// uin = xin*aux; vin = yin*aux; win = distance*aux;

	}
}
