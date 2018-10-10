package danfulea.phys.egs;

import java.io.FileInputStream;
import java.util.Calendar;
import java.util.Date;

import javax.swing.JTextArea;

/**
 * Class for calculation of the energy fraction lost to radiation when electrons slow down, if the incident beam is photons, 
 * or the radiative yield, if the incident beam is electrons. For photons, the code also calculates mu_tr (proportional to KERMA) and mu_en (proportional to absorbed dose). 
 * The energy fraction lost to radiation, the g-factor, is used to compute the absorbed dose (D) from kerma (K): D = K x (1-g). 
 * Based on EGS4/EGSnrc routines developed by SLAC(Stanford Linear Accelerator Center)/NRC (National Research Council of Canada).
 * 
 * @author Dan Fulea, 25 MAR. 2007
 */
public class G implements EgsQuestion {
	
	// "Variables definitions
	// "
	// "g e_brem/e_tot and e_rad/e_tot
	// " from NIST web page:the average fraction of the kinetic energy
	// " of secondary charged particles (produced in all the types of
	// " interactions) that is subsequently lost in radiative (photon-
	// " emitting) energy-loss processes as the particles slow to rest
	// " in the medium.
	// " http://physics.nist.gov/PhysRefData/XrayMassCoef/chap3.html
	// "e_brem average energy lost to bremsstrahlung
	// "ebrem_tmp
	// "e_rad average energy lost to bremsstrahlung & kinetic energy
	// " transferred in annihilation radiation (ie only annihilation
	// " in flight. Note this used to include fluorescent photons, but
	// " they should be excluded from e_tot as well.
	// "erad_tmp
	// "etot_tmp includes subthreshold energy (iarg=4)
	// "
	// "e_tot average energy released per particle
	// "e_radc = sum (etot_tmp*erad_tmp)
	// "e_bremc= sum (etot_tmp*ebrem_tmp)
	// "anorm
	// "npgi
	// "npei
	// "E_ave
	// "factor 1.602176462E-13 J/MeV *1000 g/kg =160.2176E-12 Gy cm^2
	// " to convert e_mutr and e_muen scored as MeV cm^2/g
	// "de_pulsei
	// "e_mutr
	// "e_muen
	// "mutr
	// "muen
	// "Eave average energy of incident spectrum: either actual or sampled

	public static JTextArea jta;
	public static boolean systemOutB = true;

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
	// public static double[] mutrs;
	// public static double[] muens;
	// public static double[] gs;

	/*
	 * REPLACE {$ANALYZE#(#);} WITH {; {P1} = {P1}/{P2}; {P1}2 = {P1}2/{P2};
	 * {P1}2 = {P1}2 - {P1}*{P1}; IF( {P1}2 > 0 ) [ {P1}2 =
	 * sqrt({P1}2/({P2}-1)); ] };
	 */
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

	/**
	 * Constructor. This application calculates the energy fraction lost to radiation when 
	 * electrons slow down if the incident beam is made of photons, or the radiative yield 
	 * if the incident beam is made of electrons. For photons, the code also calculates the 
	 * mass-energy transfer and absorbtion coefficients (mu_tr/rho,mu_en/rho).
	 */
	public G() {
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
		// --Macro and defaults param:
		EGS4.ispmfp = EGS4.iGe;// select photon mfp
		EGS4.iprint = 2;// summary
		// =========================================inputs
		Calendar cal = Calendar.getInstance();
		Date d = cal.getTime();
		EGS4.seqStr = "=================================================================================";
		EGS4.seqStr = " ********************g Kerma, mu_en && mu_tr: ************************************";
		if (EGS4.iprint > 1)
			printSequence(EGS4.seqStr);
		EGS4.seqStr = " It calculates the energy fraction lost to radiation when electrons slow down, ";
		if (EGS4.iprint > 1)
			printSequence(EGS4.seqStr);
		EGS4.seqStr = " if the incident beam is photons, or the radiative yield, if the incident beam is electrons.";
		if (EGS4.iprint > 1)
			printSequence(EGS4.seqStr);
		EGS4.seqStr = " For photons, the code also calculates the mass-energy transfer and absortion coefficients (mu_tr/rho,mu_en/rho)";
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
			return;
		}
		esum = 0.0;
		esum2 = 0.0;
		ecount = 0.0;

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

		EGS4.seqStr = "  Start g values calculations...";
		if (EGS4.iprint > 1)
			printSequence(EGS4.seqStr);

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

		// "call test_brems;
		// "call test_compton;

		nbatch = 10;
		nperbatch = ncase / nbatch;
		if (nperbatch == 0)
			nperbatch = 1;
		// "Note that nbatch and nperbatch are not used for statistical analysis"
		// "They are merely for printing information about the progress of the"
		// "simulation"

		source_sumry();

		factor = 160.2176462; // "e_mutr and e_muen scored as MeV cm^2/g"
								// "*1.602176462E-13 J/MeV *1000 g/kg =160.2176E-12 Gy cm^2"
								// "Note, this is used elsewhere too"

		boolean theEnd = false;
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

			// call source_switch_energy(itimes);
			source_switch_energy(itimes);

			if (calc_type == 1) {
				// write(6,*); write(6,*) 'First calculating mutr only';
				// write(6,*);
				// x(1),y(1),z(1)/=0; wt(1)=1; ir(1)=1; medium=1;
				EGS4.X[0] = 0.0;
				EGS4.Y[0] = 0.0;
				EGS4.Z[0] = 0.0;
				EGS4.WT[0] = 1.0;
				EGS4.IR[0] = 1;
				EGS4.MEDIUM = 1;
				for (icase = 1; icase <= ncase; icase++) {
					// call
					// source_sample(iqin,irin,ein,xin,yin,zin,uin,vin,win,wtin);
					source_sample();
					if (iqin != 0) {
						// write(6,*) 'type=1 calculation only works for
						// photons!';
						// $CALL_EXIT(1);
						EGS4.STOPPROGRAM = true;
						EGS4.seqStr = "type=1 calculation only works for photons!";
						printSequence(EGS4.seqStr);
						return;

					}// ]
					my_gle = Math.log(ein);
					EGS4.GLE = my_gle;
					// $SET INTERVAL my_gle,ge;
					Double dbl = new Double(EGS4.GE1[EGS4.MEDIUM - 1] * my_gle
							+ EGS4.GE0[EGS4.MEDIUM - 1]);
					lmy_gle = dbl.intValue();
					// $EVALUATE gmfp USING gmfp(my_gle);
					gmfp = EGS4.GMFP1[lmy_gle - 1][EGS4.MEDIUM - 1] * my_gle
							+ EGS4.GMFP0[lmy_gle - 1][EGS4.MEDIUM - 1];
					// $EVALUATE gbr1 USING gbr1(my_gle);
					gbr1 = EGS4.GBR11[lmy_gle - 1][EGS4.MEDIUM - 1] * my_gle
							+ EGS4.GBR10[lmy_gle - 1][EGS4.MEDIUM - 1];
					// $EVALUATE gbr2 USING gbr2(my_gle);
					gbr2 = EGS4.GBR21[lmy_gle - 1][EGS4.MEDIUM - 1] * my_gle
							+ EGS4.GBR20[lmy_gle - 1][EGS4.MEDIUM - 1];

					// np = 1; e(np) = ein; u(np)=0; v(np)=0; w(np)=1; iq(np) =
					// 0;
					EGS4.NP = 1;
					EGS4.E[EGS4.NP - 1] = ein;
					EGS4.U[EGS4.NP - 1] = 0.0;
					EGS4.V[EGS4.NP - 1] = 0.0;
					EGS4.W[EGS4.NP - 1] = 1.0;
					EGS4.IQ[EGS4.NP - 1] = 0;
					// $RANDOMSET rnno; edep = 0;
					rnno = EGS4.random01();
					EGS4.EDEP = 0;
					if (rnno < gbr1 && ein > EGS4.RMT2) {
						EGS4Core.PAIR();
					}// { call pair; }
					else if (rnno < gbr2) {
						EGS4Core.COMPT();
					}// [ call compt; ]
					else {
						EGS4Core.PHOTO();
					}// [ call photo; ]
					for (ip = 1; ip <= EGS4.NP; ip++)// np [
					{
						// if( iq(ip) ~= 0 ) edep = edep + e(ip)-prm;
						if (EGS4.IQ[ip - 1] != 0)
							EGS4.EDEP = EGS4.EDEP + EGS4.E[ip - 1] - EGS4.PRM;
					}
					nmutr = nmutr + 1;
					aux = EGS4.EDEP / gmfp / EGS4.RHO[0];// rho(1);
					e_mutr = e_mutr + aux;
					e_mutr2 = e_mutr2 + aux * aux;

					// if( mod(icase,10000) = 0 & accu > 0 )
					if (icase % 10000 == 0)// & accu > 0 )
					{
						xtest = e_mutr;
						xtest2 = e_mutr2;
						// $ANALYZE xtest(nmutr);
						// REPLACE {$ANALYZE#(#);} WITH {;
						// {P1} = {P1}/{P2}; {P1}2 = {P1}2/{P2};
						// {P1}2 = {P1}2 - {P1}*{P1};
						// IF( {P1}2 > 0 ) [ {P1}2 = sqrt({P1}2/({P2}-1)); ]
						// };
						xtest = xtest / nmutr;
						xtest2 = xtest2 / nmutr;
						xtest2 = xtest2 - xtest * xtest;
						if (xtest2 > 0.0) {
							xtest2 = Math.sqrt(xtest2 / (nmutr - 1));
						}
						// "write(6,*) icase,xtest2/xtest;"
						// IF( xtest2/xtest < accu ) EXIT;
					}
				}// for
				xtest = e_mutr;
				xtest2 = e_mutr2; // $ANALYZE xtest(nmutr);
				xtest = xtest / nmutr;
				xtest2 = xtest2 / nmutr;
				xtest2 = xtest2 - xtest * xtest;
				if (xtest2 > 0.0) {
					xtest2 = Math.sqrt(xtest2 / (nmutr - 1));
				}

				// IF( accu > 0 & xtest2/xtest < accu ) [
				// write(6,*) 'Converged after ',nmutr,' histories'; write(6,*);
				// ]
				// write(6,*) ' <E*mu_tr/rho> = ',factor*xtest,' 10^-12 Gy cm^2
				// +/- ',
				// 100*xtest2/xtest,'%';
				// write(6,*);
				EGS4.seqStr = "***   " + "<E*mu_tr/rho>     = "
						+ EGS4.format(factor * xtest, 12)
						+ " 10^-12 Gy cm^2  +/- "
						+ EGS4.format(100 * xtest2 / xtest, 12) + " %";
				// if(EGS4.iprint>1)
				printSequence(EGS4.seqStr);
				EGS4.seqStr = "***   " + "EDEP     = "
						+ EGS4.format(EGS4.EDEP, 12) + " MEV";
				// if(EGS4.iprint>1)
				printSequence(EGS4.seqStr);
				EGS4.seqStr = "***   " + "ein      = " + EGS4.format(ein, 12)
						+ " MEV";
				// if(EGS4.iprint>1)
				printSequence(EGS4.seqStr);
				// ein=eis[itimes-1];
				// =================
				double einc = ein;// 1.602176462E-13 J/MeV
				einc = ein * 1.602176462E-13;// J
				double miutr = factor * xtest * 1.0E-12;// Gycm2=Jcm2/kg
				miutr = miutr / einc / 1000.0;// //1kg=1000g;cm2/g
				EGS4.seqStr = "***   " + "miu_tr1  = " + EGS4.format(miutr, 12)
						+ " cm2/g";
				// if(EGS4.iprint>1)
				printSequence(EGS4.seqStr);
				// OK VERIF.............ASTA E miu_tr=> KERMA?? NU IN
				// NIST=miu_en!!!!

				/*
				 * einc=EGS4.EDEP; einc=EGS4.EDEP*1.602176462E-13;
				 * miutr=factor*xtest*1.0E-12;
				 * miutr=miutr/einc/1000.0;////1kg=1000g;cm2/g
				 * EGS4.seqStr="***   "+"miu_tr2  = "+EGS4.format(miutr,12)+
				 * " cm2/g"; //if(EGS4.iprint>1) printSequence(EGS4.seqStr);
				 */

				// ================

				for (icase = 1; icase <= ncase; icase++) {
					// call
					// source_sample(iqin,irin,ein,xin,yin,zin,uin,vin,win,wtin);
					source_sample();
					ebrem_tmp = 0.0;
					erad_tmp = 0.0;
					etot_tmp = 0.0;

					// call shower(iqin,ein,xin,yin,zin,uin,vin,win,irin,wtin);
					SHOWER();

					// $SET INTERVAL my_gle,ge;
					Double dbl = new Double(EGS4.GE1[EGS4.MEDIUM - 1] * my_gle
							+ EGS4.GE0[EGS4.MEDIUM - 1]);
					lmy_gle = dbl.intValue();

					// $EVALUATE gmfp USING gmfp(my_gle);
					gmfp = EGS4.GMFP1[lmy_gle - 1][EGS4.MEDIUM - 1] * my_gle
							+ EGS4.GMFP0[lmy_gle - 1][EGS4.MEDIUM - 1];

					aux = etot_tmp / gmfp / EGS4.RHO[0];// rho(1);
					e_mutr = e_mutr + aux;
					e_mutr2 = e_mutr2 + aux * aux;
					nmutr = nmutr + 1;
					nmuen = nmuen + 1;
					e_tot = e_tot + etot_tmp;
					e_tot2 = e_tot2 + etot_tmp * etot_tmp;
					e_brem = e_brem + ebrem_tmp;
					e_brem2 = e_brem2 + ebrem_tmp * ebrem_tmp;
					e_rad = e_rad + erad_tmp;
					e_rad2 = e_rad2 + erad_tmp * erad_tmp;
					e_radc = e_radc + etot_tmp * erad_tmp;

					// if( mod(icase,10000) = 0 & accu > 0 )
					if (icase % 10000 == 0) {
						if (e_rad > 0) {
							ert = e_rad;
							ert2 = e_rad2;
							// $ANALYZE ert(nmuen);
							ert = ert / nmuen;
							ert2 = ert2 / nmuen;
							ert2 = ert2 - ert * ert;
							if (ert2 > 0.0) {
								ert2 = Math.sqrt(ert2 / (nmuen - 1));
							}

							ett = e_tot;
							ett2 = e_tot2;
							// $ANALYZE ett(nmuen);
							ett = ett / nmuen;
							ett2 = ett2 / nmuen;
							ett2 = ett2 - ett * ett;
							if (ett2 > 0.0) {
								ett2 = Math.sqrt(ett2 / (nmuen - 1));
							}

							ert2 = ert
									* Math.sqrt((ert2 / ert) * (ert2 / ert)
											+ (ett2 / ett) * (ett2 / ett))
									/ (ett - ert);
							// IF( ert2 < accu/2 ) EXIT;
						} else {
							if (icase > 100000) {
								ert2 = 0; // EXIT;
								break;
							}
						}
					}
				}
				// IF( accu > 0 & ert2 < accu ) [
				// write(6,*) 'Converged after ',nmuen,' histories'; write(6,*);
				// ]

				// $ANALYZE e_mutr(nmutr);
				e_mutr = e_mutr / nmutr;
				e_mutr2 = e_mutr2 / nmutr;
				e_mutr2 = e_mutr2 - e_mutr * e_mutr;
				if (e_mutr2 > 0.0) {
					e_mutr2 = Math.sqrt(e_mutr2 / (nmutr - 1));
				}
				// $ANALYZE e_tot(nmuen);
				e_tot = e_tot / nmuen;
				e_tot2 = e_tot2 / nmuen;
				e_tot2 = e_tot2 - e_tot * e_tot;
				if (e_tot2 > 0.0) {
					e_tot2 = Math.sqrt(e_tot2 / (nmuen - 1));
				}
				// $ANALYZE e_rad(nmuen);
				e_rad = e_rad / nmuen;
				e_rad2 = e_rad2 / nmuen;
				e_rad2 = e_rad2 - e_rad * e_rad;
				if (e_rad2 > 0.0) {
					e_rad2 = Math.sqrt(e_rad2 / (nmuen - 1));
				}

				if (e_rad < 1.0e-15 * e_tot)
					e_rad = 1.0e-15 * e_tot;

				// write(6,*) 'Final:';
				EGS4.seqStr = "==================Final:===========================";
				// if(EGS4.iprint>1)
				printSequence(EGS4.seqStr);

				// write(6,*) ' <E*mu_tr/rho> = ',factor*e_mutr,' 10^-12 Gy cm^2
				// +/- ',
				// 100*e_mutr2/e_mutr,'%';
				EGS4.seqStr = "***   " + "<E*mu_tr/rho>     = "
						+ EGS4.format(factor * e_mutr, 12)
						+ " 10^-12 Gy cm^2  +/- "
						+ EGS4.format(100 * e_mutr2 / e_mutr, 12) + " %";
				// if(EGS4.iprint>1)
				printSequence(EGS4.seqStr);

				einc = ein * 1.602176462E-13;// J
				miutr = factor * e_mutr * 1.0E-12;// Gycm2=Jcm2/kg
				miutr = miutr / einc / 1000.0;// //1kg=1000g;cm2/g
				EGS4.seqStr = "***   " + "miu_tr1  = " + EGS4.format(miutr, 12)
						+ " cm2/g";
				// if(EGS4.iprint>1)
				printSequence(EGS4.seqStr);

				e_rad2 = e_rad
						* Math.sqrt((e_rad2 / e_rad) * (e_rad2 / e_rad)
								+ (e_tot2 / e_tot) * (e_tot2 / e_tot))
						/ (e_tot - e_rad);

				// write(6,*) ' <E*mu_en/rho> = ',factor*e_mutr*(1-e_rad/e_tot),
				// ' 10^-12 Gy cm^2 +/-
				// ',100*sqrt((e_mutr2/e_mutr)**2+e_rad2**2),'%';
				EGS4.seqStr = "***   "
						+ "<E*mu_en/rho>     = "
						+ EGS4.format(factor * e_mutr * (1 - e_rad / e_tot), 12)
						+ " 10^-12 Gy cm^2  +/- "
						+ EGS4.format(
								100 * Math.sqrt((e_mutr2 / e_mutr)
										* (e_mutr2 / e_mutr) + e_rad2 * e_rad2),
								12) + " %";
				// if(EGS4.iprint>1)
				printSequence(EGS4.seqStr);
				einc = ein * 1.602176462E-13;// J
				double miuen = factor * e_mutr * (1 - e_rad / e_tot) * 1.0E-12;// Gycm2=Jcm2/kg
				miuen = miuen / einc / 1000.0;// //1kg=1000g;cm2/g
				EGS4.seqStr = "***   " + "miu_en  = " + EGS4.format(miuen, 12)
						+ " cm2/g";
				// if(EGS4.iprint>1)
				printSequence(EGS4.seqStr);

				// write(6,*) ' g = ',e_rad/e_tot,' +/- ',
				// 100*e_rad2*(e_tot/e_rad-1),'%';
				// write(6,*);
				EGS4.seqStr = "***   " + "g         = "
						+ EGS4.format(e_rad / e_tot, 12) + " +/- "
						+ EGS4.format(100 * e_rad2 * (e_tot / e_rad - 1), 12)
						+ " %";
				// if(EGS4.iprint>1)
				printSequence(EGS4.seqStr);

				// //ein=eis[itimes-1];

				theEnd = true;
				// goto :the_end:;
			}// calc=1

			if (!theEnd) {
				for (icase = 1; icase < ncase; icase++) {
					// call
					// source_sample(iqin,irin,ein,xin,yin,zin,uin,vin,win,wtin);
					source_sample();
					// IF( $DEBUGIT ) [
					// write(18,*) ' ******* new shower, e = ',ein,' iq =
					// ',iqin;
					// ]
					ebrem_tmp = 0.0; // "keeps track of ke of charged particeles lost"
					// "to brem and annihilation events"
					erad_tmp = 0.0; // "loses just to annihilation"
					// "initialize etot_tmp differently for photons and electrons"
					if (iqin == -1) {
						etot_tmp = (ein - EGS4.RM);
					} else {
						etot_tmp = 0.0;
					}
					// "============================================================================"
					// call shower(iqin,ein,xin,yin,zin,uin,vin,win,irin,wtin);
					SHOWER();
					// "============================================================================"
					// IF( $DEBUGIT ) [
					// write(18,*) ' energy released: ',etot_tmp,' lost to
					// brems: ',ebrem_tmp,
					// ' erad: ',erad_tmp;
					// ]
					if (iqin == 0) {
						my_gle = Math.log(ein);
						// $SET INTERVAL my_gle,ge;
						Double dbl = new Double(EGS4.GE1[EGS4.MEDIUM - 1]
								* my_gle + EGS4.GE0[EGS4.MEDIUM - 1]);
						lmy_gle = dbl.intValue();
						// $EVALUATE gmfp USING gmfp(my_gle);
						gmfp = EGS4.GMFP1[lmy_gle - 1][EGS4.MEDIUM - 1]
								* my_gle
								+ EGS4.GMFP0[lmy_gle - 1][EGS4.MEDIUM - 1];
						// "when doing bound compton, the above is not the correct cross-section"
						// "but we also throw out the uncollided photons, so it works out"
						// aux = etot_tmp/gmfp/rho(1);
						aux = etot_tmp / gmfp / EGS4.RHO[0];
						e_mutr = e_mutr + aux;
						e_mutr2 = e_mutr2 + aux * aux;
						mutr = mutr + aux / ein;
						mutr2 = mutr2 + aux * aux / (ein * ein);
						// aux = (etot_tmp - erad_tmp)/gmfp/rho(1);
						aux = (etot_tmp - erad_tmp) / gmfp / EGS4.RHO[0];
						e_muen = e_muen + aux;
						e_muen2 = e_muen2 + aux * aux;
						muen = muen + aux / ein;
						muen2 = muen2 + aux * aux / (ein * ein);
					}
					e_tot = e_tot + etot_tmp;
					e_tot2 = e_tot2 + etot_tmp * etot_tmp;
					e_brem = e_brem + ebrem_tmp;
					e_brem2 = e_brem2 + ebrem_tmp * ebrem_tmp;
					e_rad = e_rad + erad_tmp;
					e_rad2 = e_rad2 + erad_tmp * erad_tmp;
					e_radc = e_radc + etot_tmp * erad_tmp;
					e_bremc = e_bremc + etot_tmp * ebrem_tmp;
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
								+ nbatch + " cpu time: "
								+ EGS4.format(timeday, 9);// +//+EGS4.format(100*e_rad2*(e_tot/e_rad-1),12)+" %";
						if (EGS4.iprint > 1)
							printSequence(EGS4.seqStr);

					}

				}// for
					// "-----------------------------------------------------------------"
					// "STEP 8   OUTPUT-OF-RESULTS                                       "
					// "-----------------------------------------------------------------"

				// call source_get_eave(Eave);
				Eave = source_get_eave();
				// write(6,*);
				// write(6,*) ' Average spectrum energy: ',Eave;
				EGS4.seqStr = " Average spectrum energy: "
						+ EGS4.format(Eave, 12);// + " cm2/g";
				// if(EGS4.iprint>1)
				printSequence(EGS4.seqStr);

				// call source_get_samplede(Eave,err_Eave);
				source_get_samplede();
				Eave = es;
				err_Eave = des;
				// write(6,*) ' Average sampled energy : ',Eave,' +/-
				// ',err_Eave;
				EGS4.seqStr = " Average sampled energy : "
						+ EGS4.format(Eave, 12) + " +/- "
						+ EGS4.format(err_Eave, 12);// + " cm2/g";
				// if(EGS4.iprint>1)
				printSequence(EGS4.seqStr);

				// "The correct definition of mu/rho for a spectrum is that given by"
				// "Attix = int[psi mu/rho dE]/int[psi dE]  This is equivalent to the"
				// "ICRU60 defn as long as one sums the total energy transferred from"
				// "the spectrum and then divides by the average energy incident"
				// "For a while we scored <mu/rho> directly, but this can differ"
				// "substantially from the correct results and amounts to averaging"
				// "over the fluence, not the energy fluence"
				// "These incorrect values are now commented out"
				// $ANALYZE e_tot(ncase);
				e_tot = e_tot / ncase;
				e_tot2 = e_tot2 / ncase;
				e_tot2 = e_tot2 - e_tot * e_tot;
				if (e_tot2 > 0.0) {
					e_tot2 = Math.sqrt(e_tot2 / (ncase - 1));
				}
				// $ANALYZE e_brem(ncase);
				e_brem = e_brem / ncase;
				e_brem2 = e_brem2 / ncase;
				e_brem2 = e_brem2 - e_brem * e_brem;
				if (e_brem2 > 0.0) {
					e_brem2 = Math.sqrt(e_brem2 / (ncase - 1));
				}
				// $ANALYZE e_rad(ncase);
				e_rad = e_rad / ncase;
				e_rad2 = e_rad2 / ncase;
				e_rad2 = e_rad2 - e_rad * e_rad;
				if (e_rad2 > 0.0) {
					e_rad2 = Math.sqrt(e_rad2 / (ncase - 1));
				}
				factor = 160.2176462; // "e_mutr and e_muen scored as MeV cm^2/g"
				// "*1.602176462E-13 J/MeV *1000 g/kg =160.2176E-12 Gy cm^2"
				// "Note, this is used elsewhere too"
				if (e_mutr > 0) {
					// write(6,'(/)');
					// $ANALYZE e_mutr(ncase);
					e_mutr = e_mutr / ncase;
					e_mutr2 = e_mutr2 / ncase;
					e_mutr2 = e_mutr2 - e_mutr * e_mutr;
					if (e_mutr2 > 0.0) {
						e_mutr2 = Math.sqrt(e_mutr2 / (ncase - 1));
					}

					// write(6,*) ' <E*mu_tr/rho> = ',factor*e_mutr,' 10^-12 Gy
					// cm^2 +/- ',
					// 100*e_mutr2/e_mutr,'%';
					EGS4.seqStr = "		   " + "<E*mu_tr/rho>     = "
							+ EGS4.format(factor * e_mutr, 12, true)
							+ " 10^-12 Gy cm^2  +/- "
							+ EGS4.format(100 * e_mutr2 / e_mutr, 12) + " %";
					// if(EGS4.iprint>1)
					printSequence(EGS4.seqStr);

					// "$ANALYZE mutr(ncase);
					mutr = mutr / ncase;
					mutr2 = mutr2 / ncase;
					mutr2 = mutr2 - mutr * mutr;
					if (mutr2 > 0.0) {
						mutr2 = Math.sqrt(mutr2 / (ncase - 1));
					}

					// "write(6,*) ' <mu_tr/rho> = ',mutr,' cm^2/g +/- ',
					// " 100*mutr2/mutr,'%';
					// String pattern="0.###E0";
					// DecimalFormatSymbols dfs=new
					// DecimalFormatSymbols(Locale.US);
					// DecimalFormat nff = new DecimalFormat(pattern,dfs);

					// EGS4.seqStr="***   "+"<mu_tr/rho>     = "+EGS4.format(mutr,12,false)+
					// " 	cm^2/g	+/- "+EGS4.format(100*mutr2/mutr,12)+" %";
					// if(EGS4.iprint>1)
					// printSequence(EGS4.seqStr);

					// write(6,*) '<E*mu_tr/rho>/Eave = ',e_mutr/Eave,' cm^2/g';
					EGS4.seqStr = "		   "
							+ "mu_tr/rho = <E*mu_tr/rho>/Eave     = "
							+ EGS4.format(e_mutr / Eave, 12, false)
							+ " 	cm^2/g " + "	+/- "
							+ EGS4.format(100 * e_mutr2 / e_mutr, 12) + " %";
					// if(EGS4.iprint>1)
					printSequence(EGS4.seqStr);

					// write(6,*) 'The above is the spectrum averaged
					// coefficient';
					EGS4.seqStr = "		   The above is the spectrum averaged coefficient";
					// if(EGS4.iprint>1)
					printSequence(EGS4.seqStr);

					// "write(6,*) ' difference = ',(e_mutr/Eave -
					// mutr)/mutr*100.,'%';
					// EGS4.seqStr=" difference        = "+EGS4.format((e_mutr/Eave
					// - mutr)/mutr*100.,12)+
					// " 	% ";
					// if(EGS4.iprint>1)
					// printSequence(EGS4.seqStr);

				}// if( e_mutr > 0 )
				if (e_muen > 0) {
					// $ANALYZE e_muen(ncase);
					e_muen = e_muen / ncase;
					e_muen2 = e_muen2 / ncase;
					e_muen2 = e_muen2 - e_muen * e_muen;
					if (e_muen2 > 0.0) {
						e_muen2 = Math.sqrt(e_muen2 / (ncase - 1));
					}

					// write(6,*);
					// write(6,*) ' <E*mu_en/rho> = ',factor*e_muen,' 10^-12 Gy
					// cm^2 +/- ',
					// 100*e_muen2/e_muen,'%';
					EGS4.seqStr = "		   " + "<E*mu_en/rho>     = "
							+ EGS4.format(factor * e_muen, 12, true)
							+ " 10^-12 Gy cm^2  +/- "
							+ EGS4.format(100 * e_muen2 / e_muen, 12) + " %";
					// if(EGS4.iprint>1)
					printSequence(EGS4.seqStr);

					// "$ANALYZE muen(ncase);
					muen = muen / ncase;
					muen2 = muen2 / ncase;
					muen2 = muen2 - muen * muen;
					if (muen2 > 0.0) {
						muen2 = Math.sqrt(muen2 / (ncase - 1));
					}

					// "write(6,*) ' <mu_en/rho> = ',muen,' cm^2/g +/- ',
					// " 100*muen2/muen,'%';
					// EGS4.seqStr="***   "+"<mu_en/rho>     = "+EGS4.format(muen,12,false)+
					// " 	cm^2/g	+/- "+EGS4.format(100*muen2/muen,12)+" %";
					// if(EGS4.iprint>1)
					// printSequence(EGS4.seqStr);

					// write(6,*) '<E*mu_en/rho>/Eave = ',e_muen/Eave,' cm^2/g';
					// write(6,*) 'The above is the spectrum averaged
					// coefficient';
					EGS4.seqStr = "		   "
							+ "mu_en/rho = <E*mu_en/rho>/Eave     = "
							+ EGS4.format(e_muen / Eave, 12, false)
							+ " 	cm^2/g " + "	+/- "
							+ EGS4.format(100 * e_muen2 / e_muen, 12) + " %";
					// if(EGS4.iprint>1)
					printSequence(EGS4.seqStr);

					// write(6,*) 'The above is the spectrum averaged
					// coefficient';
					EGS4.seqStr = "		   The above is the spectrum averaged coefficient";
					// if(EGS4.iprint>1)
					printSequence(EGS4.seqStr);

					// "write(6,*) ' difference = ',(e_muen/Eave -
					// muen)/muen*100.,'%';
					// write(6,'(/)');
					// EGS4.seqStr=" difference        = "+EGS4.format((e_muen/Eave
					// - muen)/muen*100.,12)+
					// " 	% ";
					// if(EGS4.iprint>1)
					// printSequence(EGS4.seqStr);

				}// if( e_muen > 0 )

				// write(6,*) ' Average energy released per particle: ',e_tot, '
				// +/- ',e_tot2;
				EGS4.seqStr = " Average energy released per particle:  "
						+ EGS4.format(e_tot, 12, false) + " 	+/-  "
						+ EGS4.format(e_tot2, 12, false);
				// if(EGS4.iprint>1)
				printSequence(EGS4.seqStr);

				// write(6,*) ' Average energy lost to bremsstrahlung:
				// ',e_brem,' +/- ',e_brem2;
				EGS4.seqStr = " Average energy lost to bremsstrahlung: "
						+ EGS4.format(e_brem, 12, false) + " 	+/-  "
						+ EGS4.format(e_brem2, 12, false);
				// if(EGS4.iprint>1)
				printSequence(EGS4.seqStr);

				// write(6,*) ' Average energy lost to all radiation: ',e_rad, '
				// +/- ',e_rad2;
				EGS4.seqStr = " Average energy lost to all radiation:  "
						+ EGS4.format(e_rad, 12, false) + " 	+/-  "
						+ EGS4.format(e_rad2, 12, false);
				// if(EGS4.iprint>1)
				printSequence(EGS4.seqStr);

				// write(6,*) '\n fractions g(brem) = ',e_brem/e_tot,' +/- ',
				// 100*sqrt((e_brem2/e_brem)**2+(e_tot2/e_tot)**2),'%';
				if (e_brem > 0 && e_tot > 0) {
					EGS4.seqStr = "		   fractions   g(brem)  = "
							+ EGS4.format(e_brem / e_tot, 10, false)
							+ " 	+/-  "
							+ EGS4.format(
									100 * Math.sqrt((e_brem2 / e_brem)
											* (e_brem2 / e_brem)
											+ (e_tot2 / e_tot)
											* (e_tot2 / e_tot)), 12) + " %";
					// if(EGS4.iprint>1)
					printSequence(EGS4.seqStr);
				}

				// write(6,*) ' g(all rad) = ',e_rad/e_tot,' +/- ',
				// 100*sqrt((e_rad2/e_rad)**2+(e_tot2/e_tot)**2),'%';
				if (e_rad > 0 && e_tot > 0) {
					EGS4.seqStr = "		             g(all rad) = "
							+ EGS4.format(e_rad / e_tot, 10, false)
							+ " 	+/-  "
							+ EGS4.format(
									100 * Math.sqrt((e_rad2 / e_rad)
											* (e_rad2 / e_rad)
											+ (e_tot2 / e_tot)
											* (e_tot2 / e_tot)), 12) + " %";
					// if(EGS4.iprint>1)
					printSequence(EGS4.seqStr);

					// write(6,*);
					// write(6,*)
					EGS4.seqStr = "		    the above fraction error estimates are made ignoring correlations ";
					// if(EGS4.iprint>1)
					printSequence(EGS4.seqStr);
					EGS4.seqStr = "		    between energy released and energy lost to radiation ";
					// if(EGS4.iprint>1)
					printSequence(EGS4.seqStr);
				}
				// ' the above fraction error estimates are made ignoring
				// correlations ';
				// write(6,*)
				// ' between energy released and energy lost to radiation ';

				if (e_rad > 0 && e_tot > 0 && e_brem > 0) {
					e_radc = e_radc / ncase;
					e_bremc = e_bremc / ncase;
					if (e_rad > 0) {
						e_radc = (e_radc / e_tot / e_rad - 1) / (ncase - 1);
					}
					if (e_brem > 0) {
						e_bremc = (e_bremc / e_tot / e_brem - 1) / (ncase - 1);
					}
					e_radc = (e_rad2 / e_rad) * (e_rad2 / e_rad)
							+ (e_tot2 / e_tot) * (e_tot2 / e_tot) - 2.0
							* e_radc;
					if (e_radc > 0) {
						e_radc = Math.sqrt(e_radc);
					}
					e_bremc = (e_brem2 / e_brem) * (e_brem2 / e_brem)
							+ (e_tot2 / e_tot) * (e_tot2 / e_tot) - 2.0
							* e_bremc;
					if (e_bremc > 0) {
						e_bremc = Math.sqrt(e_bremc);
					}
					// write(6,*) ' fraction errors with correlations are: ';
					EGS4.seqStr = "		    fraction errors with correlations are: ";
					// if(EGS4.iprint>1)
					printSequence(EGS4.seqStr);

					// write(6,*) ' brems: ',100*e_bremc,'%';
					// write(6,*) ' all radiative: ',100*e_radc,'%';
					// write(6,'(//)');
					EGS4.seqStr = "		      brems:         "
							+ EGS4.format(100 * e_bremc, 12) + " %";
					// if(EGS4.iprint>1)
					printSequence(EGS4.seqStr);
					EGS4.seqStr = "		      all radiative: "
							+ EGS4.format(100 * e_bremc, 12) + " %";
					// if(EGS4.iprint>1)
					printSequence(EGS4.seqStr);
				}
				// $FLUSH_UNIT(6);

				// :the_end:;
			}// if(!theEnd)
		}// ]"end of itimes loop"

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

		s = " Source type                  : ";
		ll = s.length();
		ll = 57 - ll;
		s = s + EGS4.format("", ll);
		EGS4.seqStr = s + EGS4.format(source_type, 9);// ensrcd(nensrc)
		if (EGS4.iprint > 1)
			printSequence(EGS4.seqStr);

		if (source_type == 0) {
			s = " Incident angle               : ";
			ll = s.length();
			ll = 57 - ll;
			s = s + EGS4.format("", ll);
			EGS4.seqStr = s + EGS4.format(angle, 9, true) + " degrees";// ensrcd(nensrc)
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);

		} else if (source_type == 1) {
			s = " Beam size on front face      : ";
			ll = s.length();
			ll = 57 - ll;
			s = s + EGS4.format("", ll);
			EGS4.seqStr = s + EGS4.format(rbeam, 9, true) + " cm";// ensrcd(nensrc)
			if (EGS4.iprint > 1)
				printSequence(EGS4.seqStr);
			s = " Source-face distance         : ";
			ll = s.length();
			ll = 57 - ll;
			s = s + EGS4.format("", ll);
			EGS4.seqStr = s + EGS4.format(distance, 9, true) + " cm";// ensrcd(nensrc)
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
			EGS4.seqStr = "  " + I + EGS4.format("", 3) + meds
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
	 * Setup input variables.
	 */
	public void inputs() {
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
		EGS4.iausfl[4] = 1; // "in order to score sub-threshold energy after relaxations"
		EGS4.iausfl[7] = 1; // "after brems"
		EGS4.iausfl[13] = 1;// "after annih"
		EGS4.iausfl[14] = 1;// "after annih at rest"
		EGS4.iausfl[16] = 1;// "after pair"
		EGS4.iausfl[18] = 1;// "after compt"
		EGS4.iausfl[20] = 1;// "after photo"
		EGS4.iausfl[9] = 1; // "after Moller (to count radiative losses due to EII)"

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
	 * Dummy method.
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

		// "Note that once impact ionization is included in EGSnrc, we need
		// "to change this since, g (e_rad) needs to include fluorescence
		// emission
		// "as a result of electron- and positron-impact ionization,
		// "************************************************************************""
		// $IMPLICIT-NONE;
		// $INTEGER iarg,irl,iql,jp;
		//int irl = 0;
		//int iql = 0;
		//int jp = 0;
		//double aux = 0.0;
		//double usave = 0.0;
		//double vsave = 0.0;
		//double wsave = 0.0;
		//double wtsave = 0.0;
		//double eesave = 0.0;
		//double rnno = 0.0;
		//double ee = 0.0;
		//int isplit = 0;
		//int ip = 0;

		// ;COMIN/SCORE, "to get the scoring arrays and iscore"
		// EPCONT, "to get EDEP"
		// STACK, "to get the region number"
		// USEFUL
		// /;

		if (IARG == 4) {// "depositing sub-threshold energy"
						// IF( $DEBUGIT ) [ write(19,*) ' Depositing iarg 4
						// ',edep; ]
			etot_tmp = etot_tmp + EGS4.EDEP;// edep; return;
			return;
		}
		if (IARG == 16 || IARG == 18 || IARG == 20) {// "after a photon interaction"
														// IF( $DEBUGIT ) [
														// write(18,*) ' iarg =
														// ',iarg,' np = ',np; ]

			for (ip = EGS4.NPold; ip <= EGS4.NP; ip++) {// "loop over all particles on stack after interaction"
														// IF( $DEBUGIT ) [
														// write(18,*) '
														// ',ip,iq(ip),wt(ip),e(ip);
														// ]
				if (EGS4.IQ[ip - 1] != 0)// iq(ip) ~= 0 )
				{
					// etot_tmp = etot_tmp + e(ip) - prm;
					etot_tmp = etot_tmp + EGS4.E[ip - 1] - EGS4.PRM;
				} else {// "photon"
						// wt(ip) = 0; e(ip) = 0;
					EGS4.WT[ip - 1] = 0.0;
					EGS4.E[ip - 1] = 0.0;
					// "terminate the history-no more transport"
					// "note, even a rejected bound Compton event is terminated"
				}
			}// ]"end block over particles on stack after interaction"
			return;
		}// "end photon interaction block"

		if (IARG == 7) {// "after a call to brems - add up the energy of the photon"
						// "and then discard the photon"
			if (EGS4.IQ[EGS4.NP - 1] == 0)// if( iq(np) = 0 )
			{// "photon is top of stack"
				ebrem_tmp = ebrem_tmp + EGS4.E[EGS4.NP - 1];// e(np);
				erad_tmp = erad_tmp + EGS4.E[EGS4.NP - 1];// e(np);
				// IF( $DEBUGIT ) [ write(18,*) ' brems: ',e(np),wt(np); ]
				// wt(np) = 0; e(np) = 0
				EGS4.WT[EGS4.NP - 1] = 0.0;
				EGS4.E[EGS4.NP - 1] = 0.0;
			} else {// "electron was top of stack=> photon is np-1"
				ebrem_tmp = ebrem_tmp + EGS4.E[EGS4.NP - 2];// e(np-1);
				erad_tmp = erad_tmp + EGS4.E[EGS4.NP - 2];// e(np-1);
				// IF( $DEBUGIT ) [ write(18,*) ' brems: ',e(np-1),wt(np-1); ]
				// wt(np-1) = 0; e(np-1) = 0
				EGS4.WT[EGS4.NP - 2] = 0.0;
				EGS4.E[EGS4.NP - 2] = 0.0;
			}
			return;
		}

		if (IARG == 13 || IARG == 14) {// "after annihilation at rest or in flight"
										// IF( $DEBUGIT ) [ write(18,*) '
										// annihilation '; ]
			if (IARG == 13) {// "in flight - we only want to score any k.e."
								// "transferred to photons"
								// erad_tmp = erad_tmp + e(np) + e(np-1) -
								// 2*prm;
				erad_tmp = erad_tmp + EGS4.E[EGS4.NP - 1] + EGS4.E[EGS4.NP - 2]
						- 2.0 * EGS4.PRM;
			}
			// e(np) = 0; wt(np) = 0;
			EGS4.WT[EGS4.NP - 1] = 0.0;
			EGS4.E[EGS4.NP - 1] = 0.0;
			// e(np-1) = 0; wt(np-1) = 0;
			EGS4.WT[EGS4.NP - 2] = 0.0;
			EGS4.E[EGS4.NP - 2] = 0.0;
			return;
		}

		if (IARG == 9) {// "after moller"
						// IF( $DEBUGIT ) [ write(18,*) ' Moller '; ]
			for (ip = 1; ip <= EGS4.NP; ip++)// np [
			{
				// if( iq(ip) = 0 )
				if (EGS4.IQ[ip - 1] == 0) {
					erad_tmp = erad_tmp + EGS4.E[ip - 1];// e(ip);
					// wt(ip) = 0; e(ip) = 0;
					EGS4.WT[ip - 1] = 0.0;
					EGS4.E[ip - 1] = 0.0;
				}
			}
			return;
		}

		// OUTPUT;('We should not get here!!!!');
		return;
	}// "END OF AUSGAB"

	// REPLACE {$CALL-HOWNEAR(#);} WITH {;
	// {P1} = 1e10;
	// }=>$CALL-HOWNEAR(tperp);in ELECTR egsnrc.mortran
	/**
	 * The following is a general specification of HOWNEAR: 
	 * Given a particle at (x,y,z) in region irl, HOWNEAR answers the 
	 * question, What is the distance tperp to the closest boundary? Interface method.
	 */
	public void HOWNEAR() {
		EGS4.tperp = 1.0e10;
	}// "end of subroutine HOWNEAR"

	// "Don't need any calls to howfar"
	// REPLACE {$CALL-HOWFAR-IN-ELECTR;} WITH {;}
	// REPLACE {$CALL-HOWFAR-IN-PHOTON;} WITH {;
	// IF( wt(np) <= 0 ) idisc = 1;
	// };
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
		if (EGS4.IQ[EGS4.NP - 1] == 0)// photon
		{
			if (EGS4.WT[EGS4.NP - 1] <= 0)
				EGS4.IDISC = 1;
		}
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
	 * Initialize directional cosines. Called by inputs routine.
	 */
	private void source_init() {
		if (source_type == 0) {
			win = angle / 180 * Math.PI;
			win = Math.cos(win);
			uin = Math.sqrt((1 - win) * (1 + win));
			vin = 0;
		} else if (source_type == 1) {
			// read rbeam, distance=>outside
		}
	}

	/**
	 * Read spectrum from a file.
	 * @param enerfile the spectrum file
	 */
	private void readSpectrum(String enerfile) {
		String filename = datas+file_sep+egsDatas + file_sep + dataspec + file_sep + enerfile
				+ defaultext;
		int iread = 0;
		int lnr = 0;// data number
		int lnrr = 0;// line number
		int indx = 1;
		StringBuffer desc = new StringBuffer();
		boolean haveData = false;
		//String equals = "=";
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

	@SuppressWarnings("unused")
	/**
	 * Internally used.
	 * @return the result
	 */
	private double source_get_emax()// emax);
	// "=========================="
	{
		double emax = 0.0;
		if (mono == 0 || mono == 2)// | mono = 3 ) [
		{
			// if( iqi = 0 ) [ emax = ei; ]
			if (iqin == 0) {
				emax = ein;
			} else {
				emax = ein - EGS4.RM;
			}
		} else {
			emax = ENSRCD[NENSRC];// ensrcd(nensrc);
		}
		return emax;
	}

	/**
	 * Internally used.
	 * @return the result
	 */
	private double source_get_eave()// (e);
	{
		// "========================"
		double e = 0.0;
		if (mono == 0 || mono == 2)// | mono = 3 ) [
		{
			// IF( iqi = 0 ) [ e = ei; ]
			if (iqin == 0) {
				e = ein;
			} else {
				e = ein - EGS4.RM;
			}
		} else {
			e = Eave;
		}
		return e;
	}

	/**
	 * Internally used.
	 * 
	 */
	private void source_get_samplede()// es,des);
	// "================================="
	{
		if (ecount <= 2)
			return;
		aux = esum / ecount;
		aux2 = esum2 / ecount;
		aux2 = (aux2 - aux * aux) / (ecount - 1);
		if (aux2 > 0) {
			aux2 = Math.sqrt(aux2);
		}
		es = aux;
		des = aux2;
		return;
	}

	// entry source_set_energy(ein);
	// ei = ein;
	// return;

	// private void source_sample(int iqin,int irin,double ein,double xin,double
	// yin,
	// double zin,double uin,double vin,double win,double wtin);
	/**
	 * Internally used.
	 * 
	 */
	private void source_sample() {
		// "=============================================================="
		// iqin = iqi;
		irin = 1;
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

		if (source_type == 0) {
			xin = 0;
			uin = 0;
			zin = 0;
			// uin = ui; vin = vi; win = wi;
		} else {
			double r = EGS4.random01();// $RANDOMSET r;
			r = rbeam * Math.sqrt(r);
			double phi = EGS4.random01();// $RANDOMSET phi;
			phi = 2 * phi * Math.PI;
			xin = r * Math.cos(phi);
			yin = r * Math.sin(phi);
			aux = 1 / Math.sqrt(xin * xin + yin * yin + distance * distance);
			uin = xin * aux;
			vin = yin * aux;
			win = distance * aux;
			zin = 0;
		}

		return;
	}

	/**
	 * Internally used.
	 * @param itimes itimes
	 */
	private void source_switch_energy(int itimes) {
		// "========================="
		if (mono == 1) {
			return;
		}
		ein = eis[itimes - 1];
		eik = ein;// ei=eis(itimes); eik = ei;
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
		return;
		// end; "of subroutine source"
	}

	/*
	 * "************************************************************************""
	 * subroutine test_brems;
	 * "************************************************************************"
	 * " implicit none; ;COMIN/ELECIN,EPCONT,SCORE,STACK,USEFUL/;
	 * 
	 * $INTEGER icase,iqin,irin,ip,lelke; $REAL ein, xin,yin,zin,
	 * uin,vin,win,wtin; $REAL dedx,ebr1,sig,sig_brem,eb; real*8 sum,sum2;
	 * 
	 * call source_sample(iqin,irin,ein,xin,yin,zin,uin,vin,win,wtin); np=1;
	 * e(np)=ein; iq(np)=-1; ir(np)=1; wt(np)=1; /x(np),y(np),z(np)/=0; u(np)=0;
	 * v(np)=0; w(np)=1; elke = log(ein-prm); medium=1; $SET INTERVAL elke,eke;
	 * $EVALUATE dedx USING ededx(elke); $EVALUATE ebr1 USING ebr1(elke);
	 * $EVALUATE sig USING esig(elke); sig_brem = sig*ebr1; write(6,*) '
	 * Incident energy: ',ein-prm,elke; write(6,*) ' dedx: ',dedx; write(6,*) '
	 * sig: ',sig; write(6,*) ' ebr1: ',ebr1; write(6,*) ' Brems cross section:
	 * ',sig_brem;
	 * 
	 * DO icase=1,ncase [
	 * 
	 * np=1; e(np)=ein; call brems; eb=0; DO ip=1,np [ IF(iq(ip) = 0) [ eb = eb
	 * + e(ip); ] ] sum = sum + eb; sum2 = sum2 + eb*eb; ]
	 * 
	 * sum = sum/ncase; sum2 = sum2/ncase; sum2 = sum2 - sum*sum; IF( sum2 > 0 )
	 * [ sum2 = sqrt(sum2/(ncase-1)); ] write(6,*) ' Average energy per brem:
	 * ',sum,' +/- ',sum2; write(6,*) ' Rad. stopping power: ',sum*sig_brem,'
	 * +/- ',sum2*sig_brem; stop; "return; end;
	 * 
	 * "************************************************************************""
	 * subroutine test_compton;
	 * "************************************************************************"
	 * "
	 * 
	 * REPLACE {$ANALYZE1#;} WITH {; sum{P1} = sum{P1}/ncase; sum{P1}2 =
	 * sum{P1}2/ncase; sum{P1}2 = sum{P1}2 - sum{P1}*sum{P1}; IF( sum{P1}2 > 0 )
	 * [ sum{P1}2 = sqrt(sum{P1}2/(ncase-1)); ] };
	 * 
	 * implicit none; ;COMIN/PHOTIN,EPCONT,MEDIA,SCORE,STACK,USEFUL/;
	 * 
	 * $INTEGER icase,iqin,irin,ip,lgle; $REAL ein, xin,yin,zin,
	 * uin,vin,win,wtin; $REAL gmfp,gbr1,gbr2; real*8
	 * sumr,sumr2,sume,sume2,sumtr,sumtr2; $REAL es,des,ee,factor;
	 * 
	 * medium = 1; DO icase=1,ncase [ call
	 * source_sample(iqin,irin,ein,xin,yin,zin,uin,vin,win,wtin); IF( iqin ~= 0
	 * ) [ write(6,*) ' test_compton: only works for photons!'; stop; ] gle =
	 * log(ein); $SET INTERVAL gle,ge; $EVALUATE gmfp USING gmfp(gle); $EVALUATE
	 * gbr1 USING GBR1(GLE); $EVALUATE GBR2 USING GBR2(GLE); np=1; e(np)=ein;
	 * iq(np)=0; ir(np)=1; wt(np)=1; /x(np),y(np),z(np)/=0; u(np)=0; v(np)=0;
	 * w(np)=1; call compt; IF( np > 1 ) [ IF( iq(1) = 0 ) [ ee = ein - e(1); ]
	 * ELSE [ ee = ein - e(2); ] sume = sume + ee; sume2 = sume2 + ee*ee; ee =
	 * ee*(gbr2-gbr1) + (1-gbr2)*ein; IF( ein > 2*prm & gbr1 > 0 ) [ ee = ee +
	 * gbr1*(ein-2*prm); ] ee = ee/gmfp/rho(medium); sumtr = sumtr + ee; sumtr2
	 * = sumtr2 + ee*ee; ] ELSE [ sumr = sumr + 1; sumr2 = sumr2 + 1; ] ]
	 * 
	 * $ANALYZE1 r; $ANALYZE1 e; $ANALYZE1 tr; call source_get_samplede(es,des);
	 * write(6,*) ' Average sampled energy: ',es,' +/- ',des; write(6,*) '
	 * Average number of Compton rejections: ',sumr,' +/- ',sumr2; write(6,*) '
	 * Average energy released in Compton events: ',sume,' +/- ',sume2; factor =
	 * 160.2176462; "e_mutr and e_muen scored as MeV cm^2/g"
	 * "*1.602176462E-13 J/MeV *1000 g/kg =160.2176462E-12 Gy cm^2"
	 * "Note, this is used elsewhere too" write(6,*) ' mutr: ',factor*sumtr,'
	 * +/- ',factor*sumtr2; stop; return; end;
	 */
}
