package danfulea.phys.egs;

/**
 * Core class for Monte Carlo simulation of photon-electorn transport. All physics behind particle interaction with matter are sampled here. 
 * Based on EGS4/EGSnrc routines developed by SLAC(Stanford Linear Accelerator Center)/NRC (National Research Council of Canada).
 * 
 * @author Dan Fulea, 02 OCT. 2005
 */
public class EGS4Core {
	public static int IRCODE = 0;
	// REPLACE {$EPSGMFP} WITH {1.E-5} "SMALLEST GAMMA MFP VALUE"
	public static double $EPSGMFP = 1.E-5;
	// "THRESHOLD BELOW WHICH ONLY LOWEST ORDER ANGULAR DISTRIBUTION OF THE  "
	// "PAIR ANGLE IS EMPLOYED. SCALE IS ENERGY (MeV).                       "
	// "USERS MAY OVERRIDE THIS WITH A HIGHER VALUE BUT A LOWER VALUE WILL   "
	// "CAUSE NON-PHYSICAL SAMPLING                                          "
	// "                                                                     "
	// REPLACE {$BHPAIR} WITH {4.14}
	public static double $BHPAIR = 4.14;
	// REPLACE {$RELAX-CUTOFF} WITH {0.001}
	public static double $RELAX_CUTOFF = 0.001;
	// REPLACE {$MXVAC} WITH {50}
	public static int $MXVAC = 50;
	// REPLACE {$EPSEMFP} WITH {1.E-5} "SMALLEST ELECTRON MFP VALUE"
	public static double $EPSEMFP = 1.E-5;
	// REPLACE {$RANDOMIZE-TUSTEP} WITH {.false.}
	public static boolean $RANDOMIZE_TUSTEP = false;
	// *********EMF*************************
	// REPLACE {$EMULMT} WITH {0.02}
	public static double $EMULMT = 0.02;
	// REPLACE {$EMELMT} WITH {0.02}
	public static double $EMELMT = 0.02;
	// REPLACE {$EMFLMT} WITH {0.02}
	public static double $EMFLMT = 0.02;
	// REPLACE {$EMMLMT} WITH {0.20}
	public static double $EMMLMT = 0.02;
	// ***************************************

	private static double uscat = 0.0;// "x-axis direction cosine for scattering"
	private static double vscat = 0.0;// "y-axis direction cosine for scattering"
	private static double wscat = 0.0;// "z-axis direction cosine for scattering"
	private static double xtrans = 0.0;// "final x-axis position after transport"
	private static double ytrans = 0.0;// "final y-axis position after transport"
	private static double ztrans = 0.0;// "final z-axis position after transport"
	private static boolean spin_index = false;
	private static boolean find_index = false;

	private static int sscatcallindex = 0;// for passig parameters in sscat.
											// Note in Java the args are passing
	// only by value and does not affect external variables!!!!
	private static int mscatcallindex = 0;
	private static double w1 = 0.0;
	private static double sint1 = 0.0;// 1
	private static double w2 = 0.0;
	private static double sint2 = 0.0;// 2
	private static double ws = 0.0;
	private static double sint = 0.0;// 3
	private static int msdist1call = 0;
	private static int msdist2call = 0;

	private static int ierust = 0;// "error counter for negative ustep errors"

	private static double CTHET = 0.0;// "5/2*pi-THETA, used to evaluate cos(THETA) using the sine table"
	//private static double PHI = 0.0;// "azimuthal scattering angle"
	//private static double CPHI = 0.0;// "5/2*pi-PHI"
	private static double A = 0.0;
	private static double B = 0.0;
	private static double C = 0.0;// "direction cosines before rotation"
	private static double SINPS2 = 0.0;// "SINPS2=A*A+B*B"
	private static double SINPSI = 0.0;// "Sqrt(SINPS2)"
	private static double US = 0.0;
	private static double VS = 0.0;// "x- and y- component of scattering vector"
	private static double SINDEL = 0.0;
	private static double COSDEL = 0.0;
	private static int i = 0;
	private static int j = 0;
	private static double omega2 = 0.0;
	private static int isrj = 0;
	private static int jsrj = 0;
	public static int[] first_transition = { 1, 20, 27, 33, 38 };// ($MXINTER),//$MXINTER=5
	public static int[] last_transition = { 19, 26, 32, 37, 39 };// ($MXINTER);
	// "first and last transition for a given "
	// "shell in the list of all possible     "
	// "transitions                           "
	public static int[] final_state = { 4, 3, 5, 6, // " K-shell fluorescence    "
			202, 302, 402, 404, 403, 303, // " K-shell Auger           "
			502, 503, 504, 602, 603, 604, // " K-shell Auger           "
			505, 605, 606, // " K-shell Auger           "
			13, 14, // " L1 Coster-Kronig        "
			5, 6, // " L1 fluorescence         "
			505, 605, 606, // " L1 Auger                "
			14, // " L2 Coster-Kronig        "
			5, 6, // " L2 fluorescence         "
			505, 605, 606, // " L2 Auger                "
			5, 6, // " L3 fluorescence         "
			505, 605, 606, // " L3 Auger                "
			6, // " M  fluorescence         "
			606 };// ($MXTRANS);//=39
					// " final_state(i) is the final atomic state                "
					// " after transition i coded as follows:                    "
					// "   * fluorescence - final_state is the shell number      "
					// "                    of the new vacancy                   "
					// "   * Coster-Kronig - final_state is the shell number     "
					// "                     of the new vacancy + 10             "
					// "   * Auger - final_state is n1 + 100*n2 where n1 and n2  "
					// "             are the shell numbers of the 2 new vacancies"
	public static int n_warning = 0;// "a warning counter"
	// triplket
	private static boolean is_initialized = false;
	private static double[] fmax_array = new double[EGS4.$MAX_TRIPLET];
	private static double[] eta_p_array = new double[EGS4.$MAX_TRIPLET];// ($MAX_TRIPLET),
	private static double[] eta_Ep_array = new double[EGS4.$MAX_TRIPLET];// ($MAX_TRIPLET),
	private static double[] eta_costp_array = new double[EGS4.$MAX_TRIPLET];// ($MAX_TRIPLET),
	private static double[] eta_costm_array = new double[EGS4.$MAX_TRIPLET];// ($MAX_TRIPLET),
	private static double[] ebin_array = new double[EGS4.$MAX_TRIPLET];// ($MAX_TRIPLET),
	private static double[] wp_array = new double[EGS4.$MAX_TRIPLET];// ($MAX_TRIPLET),
	private static double[] qmin_array = new double[EGS4.$MAX_TRIPLET];// ($MAX_TRIPLET);
	private static double kmin = 0.0;
	private static double kmax = 0.0;
	private static double dlogki = 0.0;
	private static double alogkm = 0.0;
	private static double prmi = 0.0;
	private static double tiny_eta = 0.0;

	// private int IQI=0;
	// private double EI=0.0;
	// private double XI=0.0;
	// private double YI=0.0;
	// private double ZI=0.0;
	// private double UI=0.0;
	// private double VI=0.0;
	// private double WI=0.0;
	// private int IRI=0;
	// private double WTI=0.0;
	public static EgsQuestion eq;

	// SEE EGS4
	// public EGS4Core(int IQI,double EI,double XI,double YI,double ZI,double
	// UI,double VI,double WI,
	// int IRI,double WTI, EgsQuestion eq)
	// {
	// this.IQI=IQI;
	// this.EI=EI;
	// this.XI=XI;
	// this.YI=YI;
	// this.ZI=ZI;
	// this.UI=UI;
	// this.VI=VI;
	// this.WI=WI;
	// this.IRI=IRI;
	// this.WTI=WTI;
	// this.eq=eq;
	// SHOWER(IQI,EI,XI,YI,ZI,UI,VI,WI,IRI,WTI);

	// }

	/**
	 * Starts the "shower", the coupled photon-electron simulation.
	 * @param IQI the incident particle charge
	 * @param EI the incident particle energy (total, i.e. kinetic + restMass if electrons)
	 * @param XI the incident particle X coordinate
	 * @param YI the incident particle Y coordinate
	 * @param ZI the incident particle Z coordinate
	 * @param UI the incident particle U directional cosine
	 * @param VI the incident particle V directional cosine
	 * @param WI the incident particle W directional cosine
	 * @param IRI the incident particle region index
	 * @param WTI the incident particle weight
	 */
	public static void SHOWER(int IQI, double EI, double XI, double YI,
			double ZI, double UI, double VI, double WI, int IRI, double WTI) {
		// ; Copyright NRC;

		// $COMIN-SHOWER; "DEFAULT REPLACEMENT PRODUCES THE FOLLOWING:
		// "COMIN/DEBUG,STACK,UPHIOT,RANDOM/;

		// "Input variables"
		// $REAL EI, "initial shower energy"
		// XI,YI,ZI,"initial co-ordinates"
		// UI,VI,WI,"initial direction cosines"
		// WTI; "initial weight"

		// $INTEGER
		// IQI, "initial particle charge"
		// IRI; "initial region number"

		// "Local variables"
		double DEG = 0.0;// "energy for pi-zero option"
		double DPGL = 0.0;// "angle factor for pi-zero option"
		double DEI = 0.0;// "incident energy for pi-zero option"
		double DPI = 0.0;// "intermediate factor for pi-zero option"
		double DCSTH = 0.0;// "random number for pi-zero option"
		double DCOSTH = 0.0;// "cos(theta) for pi-zero option"
		double PI0MSQ = 0.0;// "pi-zero mass squared (in MeV**2)"

		double DNEARI = 0.0; // "initial distance to closest boundary"
		double CSTH = 0.0; // "random number for pi-zero option"

		// int IRCODE=0; //"status returned by ELECTR or PHOTON"
		EGS4.STOPPROGRAM = false;// make it false if here
		IRCODE = 0;

		PI0MSQ = 1.8215416E4; // "PI-ZERO MASS (MEV) SQUARED"

		EGS4.NP = 1;
		EGS4.NPold = EGS4.NP; // "Set the old stack counter"
		DNEARI = 0.0;
		// EGS4.IQ[1]=IQI; EGS4.E[1]=EI; EGS4.U[1]=UI; EGS4.V[1]=VI;
		// EGS4.W[1]=WI;
		EGS4.IQ[0] = IQI;
		EGS4.E[0] = EI;
		EGS4.U[0] = UI;
		EGS4.V[0] = VI;
		EGS4.W[0] = WI;
		// $TRANSFER PROPERTIES TO (1) FROM I;
		// "MACRO TO ALLOW USER TO ADD TO PROPERTIES THAT ARE"
		// "PASSED TO NEW PARTICLES"
		// REPLACE {$TRANSFERPROPERTIESTO#FROM#;} WITH {
		// X{P1}=X{P2};Y{P1}=Y{P2};Z{P1}=Z{P2};IR{P1}=IR{P2};
		// WT{P1}=WT{P2};DNEAR{P1}=DNEAR{P2};LATCH{P1}=LATCH{P2}
		// X(1)=XI;Y(1)=YI;Z(1)=ZI;IR(1)=IRI;
		// WT(1)=WTI;DNEAR(1)=DNEARI;LATCH(1)=LATCHI
		EGS4.X[0] = XI;
		EGS4.Y[0] = YI;
		EGS4.Z[0] = ZI;
		EGS4.IR[0] = IRI;
		EGS4.WT[0] = WTI;
		EGS4.DNEAR[0] = DNEARI;
		EGS4.LATCH[0] = EGS4.LATCHI;

		// ;}
		// "IN THE USAGE OF THE ABOVE MACRO, '(IP)' WILL REFER TO THE"
		// "PARTICLE WHOSE STACK INDEX IS 'IP'. 'NP' IS THE TOP OF THE STACK."
		// "JUST PLAIN 'I' WILL REFER TO THE INITIAL VALUES SUPPLIED"
		// "AS ARGUMENTS TO SHOWER, OR SUPPLIED IN COMMON BLOCKS OR"
		// "DATA STATEMENTS IN SHOWER."

		if (IQI == 2) {// "PI-ZERO OPTION"
						// "IF(EI <= PI0MSQ) [OUTPUT EI;    corrected Oct 24 1995 e-mail Hideo H "
						// " noted by Dr. Muroyama at Nagoya University
			if (EI * EI <= PI0MSQ) {
				// OUTPUT EI;
				EGS4.STOPPROGRAM = true;
				EGS4.seqStr = " Stopped in subroutine SHOWER---PI-ZERO option invoked"
						+ " but the total energy was too small (EI="
						+ EI
						+ " MeV)";
				// if(EGS4.iprint>2)
				eq.printSequence(EGS4.seqStr);
				return;
			}
			// $RANDOMSET CSTH;
			CSTH = EGS4.random01();
			DCSTH = CSTH;
			DEI = EI;
			DPI = Math.sqrt(DEI * DEI - PI0MSQ);
			DEG = DEI + DPI * DCSTH;
			DPGL = DPI + DEI * DCSTH;
			DCOSTH = DPGL / DEG;
			EGS4.COSTHE = DCOSTH;
			EGS4.SINTHE = Math.sqrt(1.0 - DCOSTH * DCOSTH);
			// IQ(1)=0; E(1)=DEG/2.;
			EGS4.IQ[0] = 0;
			EGS4.E[0] = DEG / 2.;
			UPHI(2, 1);// $SELECT-AZIMUTHAL-ANGLE and OLD-PARTICLE:
			EGS4.NP = 2;
			DEG = DEI - DPI * DCSTH;
			DPGL = DPI - DEI * DCSTH;
			DCOSTH = DPGL / DEG;
			EGS4.COSTHE = DCOSTH;
			EGS4.SINTHE = -Math.sqrt(1.0 - DCOSTH * DCOSTH);
			// IQ(2)=0; E(2)=DEG/2.;
			EGS4.IQ[1] = 0;
			EGS4.E[1] = DEG / 2.;
			UPHI(3, 2);// NEW-PARTICLE
		}// "end of pi-zero option"

		// "The following convoluted logic is difficult to follow"
		// "when one modifies the outcome of certain interactions"
		// "using nbr_split, Russian Roulette, or one of the     "
		// "particle selection macros. I'm simplifying it        "
		// "so that ircode becomes irrelevant. IK, August 2002   "
		// ":TOPSTACK:"
		// "LOOP["
		// "   $KERMA-INSERT;"
		// "   IF(IQ(NP) = 0) GO TO :PHOTON:;"
		// "   LOOP ["
		// "      :ELECTRON:"
		// "      CALL ELECTR(IRCODE);"
		// "      IF(IRCODE.EQ.2) EXIT; "
		// "      :PHOTON:"
		// "      CALL PHOTON(IRCODE);"
		// "      IF(IRCODE.EQ.2) EXIT;"
		// "   ]REPEAT"
		// "   IF(NP <= 0) EXIT;"
		// "]REPEAT "

		// LOOP [
		while (true) {
			// ;
			if (EGS4.NP <= 0)
				break;
			// $KERMA-INSERT; //" DEFAULT FOR $KERMA-INSERT; IS ; (NULL)"
			// IF( iq(np) = 0 ) [ call photon(ircode); ]
			// ELSE [ call electr(ircode); ]
			if (EGS4.IQ[EGS4.NP - 1] == 0) {
				PHOTON();// IRCODE);//ircode);
			} else {
				ELECTR();// IRCODE);//(ircode);
			}

			if (EGS4.STOPPROGRAM) {
				return;
			}
		}

		return;
		// "end of subroutine shower" END;

	}

	// "******************************************************************"
	// "                               NATIONAL RESEARCH COUNCIL OF CANADA"
	/**
	 * Called by SHOWER. Handles electron interaction with matter.
	 */
	protected static void ELECTR()// (int ircode)
	{
		// "******************************************************************"
		// "   This subroutine has been almost completely recoded to include  "
		// "   the EGSnrc enhancements.                                       "
		// "                                                                  "
		// "   Version 1.0   Iwan Kawrakow       Complete recoding            "
		// "   Version 1.1   Iwan Kawrakow       Corrected implementation of  "
		// "                                     fictitious method (important "
		// "                                     for low energy transport     "
		// "******************************************************************"

		// ; Copyright NRC;

		// $COMIN-ELECTR; "default replacement produces the following:
		// "COMIN/DEBUG,BOUNDS,EGS-VARIANCE-REDUCTION, ELECIN,EPCONT,"
		// "ET-Control,MEDIA,MISC,STACK,THRESH,UPHIIN,"
		// "UPHIOT,USEFUL,USER,RANDOM/;"
		// ;COMIN/EII-DATA/;

		// $DEFINE-LOCAL-VARIABLES-ELECTR;
		// REPLACE {$DEFINE-LOCAL-VARIABLES-ELECTR;} WITH
		// {;
		// " Local ELECTR variables"
		// $ENERGY PRECISION "($ENERGY PRECISION means double precision)"
		double demfp = 0.0;// "differential electron mean free path"
		double peie = 0.0;// "precise energy of incident electron"
		double total_tstep = 0.0;// "total path-length to next discrete interaction"
		double total_de = 0.0;// "total energy loss to next discrete interaction"

		double ekems = 0.0;// "kinetic energy used to sample MS angle (normally midpoint)"
		double elkems = 0.0;// "Log(ekems)"
		double chia2 = 0.0;// "Multiple scattering screening angle"
		double etap = 0.0;// "correction to Moliere screening angle from PWA cross sections"
		double lambda = 0.0;// "number of mean free paths (elastic scattering cross section)"
		double blccl = 0.0;// "blcc(medium)*rhof"
		double xccl = 0.0;// "xcc(medium)*rhof"
		double xi = 0.0;// "used for PLC calculations (first GS moment times path-length)"
		double xi_corr = 0.0;// "correction to xi due to spin effects"
		double ms_corr = 0.0;
		double p2 = 0.0;// "electron momentum times c, squared"
		double beta2 = 0.0;// "electron speed in units of c, squared"
		double de = 0.0;// "energy loss to dedx"
		double save_de = 0.0;// "de saved before $DE-FLUCTUATION"
		double dedx = 0.0;// "stopping power after density scaling"
		double dedx0 = 0.0;// "stopping power before density scaling"
		double dedxmid = 0.0;// "stopping power at mid-step before density scaling"
		double ekei = 0.0;// "used in $CALCULATE-TSTEP-FROM-DEMFP;"
		double elkei = 0.0;// "Log(ekei), used in $CALCULATE-TSTEP-FROM-DEMFP;"
		double aux = 0.0;// "aux. variable"
		double ebr1 = 0.0;// "e- branching ratio into brem"
		double eie = 0.0;// "energy of incident electron"
		double ekef = 0.0;// "kinetic energy after a step"
		double elkef = 0.0;// "Log(ekef)"
		//double ekeold = 0.0;// "kinetic energy before a step"
		double eketmp = 0.0;// "used to evaluate average kinetic energy of a step"
		double elktmp = 0.0;// "log(eketmp)"
		double fedep = 0.0;// "fractional energy loss used in stopping power calculation"
		double tuss = 0.0;// "sampled path-length to a single scattering event"
		double pbr1 = 0.0;// "e+ branching ratio into brem"
		double pbr2 = 0.0;// "e+ branching ratio into brem or Bhabha"
		double range = 0.0;// "electron range"
		double rfict = 0.0;// "rejection function for fictitious cross section"
		double rnne1 = 0.0;// "random number"
		double rnno24 = 0.0;// "random number"
		double rnno25 = 0.0;// "random number"
		double rnnotu = 0.0;// "random number"
		double rnnoss = 0.0;// "random number"
		double sig = 0.0;// "cross section after density scaling but before a step"
		double sig0 = 0.0;// "cross section before density scaling but before a step"
		double sigf = 0.0;// "cross section before density scaling but after a step"
		double skindepth = 0.0;// "skin depth employed for PRESTA-II boundary crossing"
		double ssmfp = 0.0;// "distance of one single elastic scattering mean free path"
		double tmxs = 0.0;// "electron step-size restriction"
		double tperp = 0.0;// "perpendicular distance to the closest boundary"
		double ustep0 = 0.0;// "temporary storage for ustep"
		// double uscat=0.0;// "x-axis direction cosine for scattering"
		// double vscat=0.0;// "y-axis direction cosine for scattering"
		// double wscat=0.0;// "z-axis direction cosine for scattering"
		// double xtrans=0.0;// "final x-axis position after transport"
		// double ytrans=0.0;// "final y-axis position after transport"
		// double ztrans=0.0;// "final z-axis position after transport"
		//double cphi = 0.0;
		//double sphi = 0.0;// "for azimuthal angle selection for annih at rest"

		// $DEFINE-VARIABLES-FOR-SELECT-AZIMUTHAL-ANGLE;
		// REPLACE {$DEFINE-VARIABLES-FOR-SELECT-AZIMUTHAL-ANGLE;} WITH {;
		//double xphi = 0.0;
		//double xphi2 = 0.0;
		//double yphi = 0.0;
		//double yphi2 = 0.0;
		//double rhophi2 = 0.0;
		// };
		// THESE MACROS ARE USED BY ELECTR FOR CALCULATIONS WITH AN EXTERNAL em
		// FIELD. "
		double x0 = 0.0;
		double y0 = 0.0;
		double z0 = 0.0;
		//double tustp0 = 0.0;
		//double fbtemp = 0.0;
		double Ex0 = 0.0;
		double Ey0 = 0.0;
		double Ez0 = 0.0;
		double Bx0 = 0.0;
		double By0 = 0.0;
		double Bz0 = 0.0;
		double ekin0 = 0.0;
		double beta20 = 0.0;
		double beta0 = 0.0;
		double gamma0 = 0.0;
		double fnorm = 0.0;
		double fnormb = 0.0;
		double u0 = 0.0;
		double v0 = 0.0;
		double w0 = 0.0;
		double edotu = 0.0;
		double dperpx = 0.0;
		double dperpy = 0.0;
		double dperpz = 0.0;
		double dperp2 = 0.0;
		double bperpx = 0.0;
		double bperpy = 0.0;
		double bperpz = 0.0;
		double bperp2 = 0.0;
		double eperpx = 0.0;
		double eperpy = 0.0;
		double eperpz = 0.0;
		double eperp2 = 0.0;
		double eperp = 0.0;
		double xf = 0.0;
		double yf = 0.0;
		double zf = 0.0;
		double test = 0.0;
		double x1 = 0.0;
		double y1 = 0.0;
		double z1 = 0.0;
		double efx0 = 0.0;
		double efy0 = 0.0;
		double efz0 = 0.0;
		//double efxf = 0.0;
		//double efyf = 0.0;
		//double efzf = 0.0;
		//double bfx0 = 0.0;
		//double bfy0 = 0.0;
		//double bfz0 = 0.0;
		double pot1 = 0.0;
		double pot2 = 0.0;
		double potdif = 0.0;
		double ufx = 0.0;
		double ufy = 0.0;
		double ufz = 0.0;
		double bdotub = 0.0;
		double bsqrd = 0.0;
		double ufxpar = 0.0;
		double ufypar = 0.0;
		double ufzpar = 0.0;
		double ufxprp = 0.0;
		double ufyprp = 0.0;
		double ufzprp = 0.0;
		boolean EMFb = EGS4.inEMF;
		int iemf = EGS4.i_EMF;

		int iarg = 0;// "calling code for eq.AUSGAB"
		int idr = 0;// "calling code for eq.AUSGAB"
		// int ierust=0;// "error counter for negative ustep errors"
		int irl = 0;// "region number"
		int lelec = 0;// "charge of electron"
		int qel = 0;// " = 0 for electrons, = 1 for positrons "
		int lelke = 0;// "index into the energy grid of tabulated functions"
		int lelkems = 0;// "index into the energy grid of tabulated functions"
		int lelkef = 0;// "index into the energy grid of tabulated functions"
		int lelktmp = 0;// "index into the energy grid of tabulated functions"
		//int ibr = 0;// "a loop variable"

		// "BCA = boundary crossing algorithm"
		boolean callHOWFAR = false;// "= .true.  => BCA requires a call to eq.HOWFAR"
		// "= .false. => BCA does not require a call to eq.HOWFAR"
		boolean domultiple = false;// "= .true.  => inexact BCA requires multiple scattering"
		boolean dosingle = false;// "= .true.  => exact BCA requires single scattering"
		// "= .false. => exact BCA requires no single scattering"
		boolean callmsdist = false;// "= .true.  => normal condensed-history transport"
		// "= .false. => one of the BCA's will be invoked"
		//boolean findindex = false;// "used for mscat"
		// boolean spin_index=false;// "used for mscat with spin effects"
		boolean compute_tstep = false;
		// ;
		// }

		double lambda_max = 0.0;
		double sigratio = 0.0;
		boolean random_tustep = false;
		boolean NEWELECTRONb = false;
		boolean EBREMSb = false;
		boolean ECUTb = false;
		boolean USERb = false;
		boolean NEXTb = false;

		// ierust=0; // "To count negative ustep's"
		// save ierust;

		// $CALL-USER-ELECTRON;REPLACE {$CALL-USER-ELECTRON} WITH {;}

		IRCODE = 1; // "Set up normal return-which means there is a photon
					// "with less available energy than the lowest energy
					// electron,
					// "so return to shower so it can call photon to follow it.
					// "(For efficiency's sake, we like to stay in this routine
					// " as long as there are electrons to process. That's why
					// this
					// " apparently convoluted scheme of STACK contro is
					// effected.)

		EGS4.IROLD = EGS4.IR[EGS4.NP - 1]; // "Initialize previous region
		// "(ir() is an integer that is attached to the particle's
		// " phase space. It contains the region
		// " number that the current particle is in.
		// " Np is the stack pointer, it points to where on the
		// " stack the current particle is.)
		irl = EGS4.IROLD; // "region number in local variable

		// $start_new_particle;
		// REPLACE {$start_new_particle;} WITH { medium = med(irl); };
		EGS4.MEDIUM = EGS4.MED[irl - 1];// eg. med[1]=1 mediul 1

		// " Default replacement for the above is medium = med(irl); "
		// " This is made a macro so that it can be replaced with a call to a "
		// " user provided function start_new_particle(); for the C/C++ interface "
		BREMSLOOP: while (true)// go to :NEWELECTRON loop for outside main loop
		{
			NEWELECTRON: while (true)// main loop
			{

				NEWELECTRONb = false;
				EBREMSb = false;
				USERb = false;

				// "Go once through this loop for each 'new' electron whose
				// charge and
				// "energy has not been checked

				lelec = EGS4.IQ[EGS4.NP - 1]; // "Save charge in local variable
				// "(iq = -1 for electrons, 0 for photons and 1 for positrons)
				qel = (1 + lelec) / 2; // " = 0 for electrons, = 1 for positrons "
				peie = EGS4.E[EGS4.NP - 1]; // "precise energy of incident
											// electron (double precision)
				eie = peie; // "energy incident electron (conversion to single)

				if (eie <= EGS4.ECUT[irl - 1])// ecut(irl))
				{
					// go to :ECUT-DISCARD:;
					ECUTb = true;
					break BREMSLOOP;// NEWELECTRON;
				}
				// "(Ecut is the lower transport threshold.)

				// "medium = med(irl);" "(This renders the above assignment
				// redundant!)
				// "The above assignment is unnecessary, IK, June 2003"

				if (EGS4.WT[EGS4.NP - 1] == 0.0) {
					// go to :USER-ELECTRON-DISCARD:;
					USERb = true;
					break BREMSLOOP;// NEWELECTRON;
				} // "added May 01"

				TSTEP: while (true)// LOOP
				{
					NEXTb = false;// reset
					// "Go through this loop each time we recompute distance to
					// an interaction
					compute_tstep = true; // "MFP resampled => calculate
											// distance to the
											// "interaction in the USTEP loop
					EGS4.EKE = eie - EGS4.RM; // "moved here so that kinetic energy will be known"
					// "to user even for a vacuum step, IK January 2000"
					if (EGS4.MEDIUM != 0) {
						// "Not vacuum. Must sample to see how far to next
						// interaction.
						if (EGS4.isemfp == 0) {
							// $SELECT-ELECTRON-MFP;
							// " Default FOR $SELECT-ELECTRON-MFP; is:
							// $RANDOMSET rnne1;
							// " demfp = -log(rnne1);
							// "($RANDOMSET is a macro'ed random number
							// generator)
							// "(demfp = differential electron mean free path)
							// REPLACE {$SELECT-ELECTRON-MFP;} WITH {
							rnne1 = EGS4.random01();
							if (rnne1 == 0.0) {
								rnne1 = 1.E-30;
							}
							demfp = Math.max(-Math.log(rnne1), $EPSEMFP);
							// }
						} else {
							// pass var==============
							EGS4Macro.demfp = demfp;
							// end pass var===========
							EGS4Macro.SELECT_ELECTRON_MFP();
							// get var==================
							demfp = EGS4Macro.demfp;
							// end get var========================
						}

						EGS4.ELKE = Math.log(EGS4.EKE);
						// "(eke = kinetic energy, rm = rest mass, all in units
						// of MeV)
						// $SET INTERVAL elke,eke; //"Prepare to approximate
						// cross section
						Double dbl = new Double(EGS4.EKE1[EGS4.MEDIUM - 1]
								* EGS4.ELKE + EGS4.EKE0[EGS4.MEDIUM - 1]);
						lelke = dbl.intValue();

						// $EVALUATE-SIG0;
						// REPLACE {$EVALUATE-SIG0;} WITH
						// "        ==============="
						// {;
						// if( sig_ismonotone(qel,medium) ) [
						if (EGS4.sig_ismonotone[qel][EGS4.MEDIUM - 1]) {
							// $EVALUATE-SIGF;
							// REPLACE {$EVALUATE-SIGF;} WITH
							// "        ==============="
							// {;
							if (lelec < 0) {
								// $EVALUATE sigf USING esig(elke);
								sigf = EGS4.ESIG1[lelke - 1][EGS4.MEDIUM - 1]
										* EGS4.ELKE
										+ EGS4.ESIG0[lelke - 1][EGS4.MEDIUM - 1];
								// $EVALUATE dedx0 USING ededx(elke);
								dedx0 = EGS4.EDEDX1[lelke - 1][EGS4.MEDIUM - 1]
										* EGS4.ELKE
										+ EGS4.EDEDX0[lelke - 1][EGS4.MEDIUM - 1];
								sigf = sigf / dedx0;
							} else {
								// $EVALUATE sigf USING psig(elke);
								sigf = EGS4.PSIG1[lelke - 1][EGS4.MEDIUM - 1]
										* EGS4.ELKE
										+ EGS4.PSIG0[lelke - 1][EGS4.MEDIUM - 1];
								// $EVALUATE dedx0 USING pdedx(elke);
								dedx0 = EGS4.PDEDX1[lelke - 1][EGS4.MEDIUM - 1]
										* EGS4.ELKE
										+ EGS4.PDEDX0[lelke - 1][EGS4.MEDIUM - 1];
								sigf = sigf / dedx0;
							}
							// }

							sig0 = sigf;
						} else {
							if (lelec < 0) {
								sig0 = EGS4.esig_e[EGS4.MEDIUM - 1];
							} else {
								sig0 = EGS4.psig_e[EGS4.MEDIUM - 1];
							}
						}
						// }

						// "The fix up of the fictitious method uses cross section per"
						// "energy loss. Therefore, demfp/sig is sub-threshold energy loss"
						// "until the next discrete interaction occures (see below)"
						// "As this quantity is a single constant for a material,"
						// "$SET INTERVAL is not necessary at this point. However, to not"
						// "completely alter the logic of the TSTEP and USTEP loops,"
						// "this is left for now"

					}// "end non-vacuum test

					USTEP: while (true)// LOOP
					{
						// "Here for each check with user geometry.
						// "Compute size of maximum acceptable step, which is
						// limited
						// "by multiple scattering or other approximations.
						if (EGS4.MEDIUM == 0) {
							// "vacuum
							// tstep = vacdst; ustep = tstep; tustep = ustep;
							EGS4.TSTEP = EGS4.VACDST;
							EGS4.USTEP = EGS4.TSTEP;
							EGS4.TUSTEP = EGS4.USTEP;
							callHOWFAR = true; // "Always call eq.HOWFAR for vacuum steps!"

							// "(Important definitions:
							// " tstep = total pathlength to the next discrete
							// interaction
							// " vacdst = infinity (actually 10^8)
							// " tustep = total pathlength of the electron step
							// " ustep = projected transport distance in the
							// " direction of motion at the start of the step
							// " Note that tustep and ustep are modified below.
							// " The above provide defaults.)
						} else {
							// "non-vacuum
							// $SET-RHOF; //"density ratio scaling template
							// "EGS allows the density to vary
							// "continuously (user option)
							// REPLACE {$SET-RHOF;} WITH
							// {RHOF=RHOR(IRL)/RHO(MEDIUM);} //"DEFAULT"
							EGS4.RHOF = EGS4.RHOR[irl - 1]
									/ EGS4.RHO[EGS4.MEDIUM - 1];

							// $SCALE-SIG0;
							// " Because the cross section is interactions per energy loss, no "
							// " rhof-scaling is required "
							// REPLACE {$SCALE-SIG0;} WITH
							// "        ============"
							// {
							sig = sig0;
							// }
							;
							if (sig <= 0.0) {
								// "This can happen if the threshold for brems,
								// "(ap + rm), is greater than ae. Moller
								// threshold is
								// "2*ae - rm. If sig is zero, we are below the
								// "thresholds for both bremsstrahlung and
								// Moller.
								// "In this case we will just lose energy by
								// "ionization loss until we go below cut-off.
								// Do not
								// "assume range is available, so just ask for
								// step
								// "same as vacuum. Electron transport will
								// reduce
								// "into little steps.
								// "(Note: ae is the lower threshold for
								// creation of a
								// " secondary Moller electron, ap is the lower
								// " threshold for creation of a brem.)
								EGS4.TSTEP = EGS4.VACDST;// tstep = vacdst;
								sig0 = 1.E-15;
							} else {
								// $CALCULATE-TSTEP-FROM-DEMFP;
								// " Once the sub-threshold processes energy loss to the next discrete "
								// " interaction is determined, the corresponding path-length has to be"
								// " calculated. This is done by the macro below. This macro           "
								// " assumes the energy at the begining to be eke, the logarithm of it "
								// " elke, lelke - the corresponding interpolation index and makes     "
								// " use of $COMPUTE-DRANGE(#,#,#,#,#,#)                               "

								// REPLACE {$CALCULATE-TSTEP-FROM-DEMFP;} WITH
								// "        ============================"
								// {;
								if (compute_tstep) {
									total_de = demfp / sig;
									fedep = total_de;
									ekef = EGS4.EKE - fedep;
									// "write(6,*) ' CALCULATE-TSTEP-FROM-DEMFP: eke ekef = ',eke,ekef;"
									if (ekef <= EGS4.E_array[0][EGS4.MEDIUM - 1])// (1,medium)
																					// )
									{
										EGS4.TSTEP = EGS4.VACDST;// tstep =
																	// vacdst;
									} else {
										elkef = Math.log(ekef);
										// $SET INTERVAL elkef,eke;
										Double dbll = new Double(
												EGS4.EKE1[EGS4.MEDIUM - 1]
														* elkef
														+ EGS4.EKE0[EGS4.MEDIUM - 1]);
										lelkef = dbll.intValue();

										if (lelkef == lelke) {
											// " initial and final energy are in the same interpolation bin "
											// "write(6,*) ' CALCULATE-TSTEP-FROM-DEMFP: same bin ',lelkef,lelke;"
											// $COMPUTE-DRANGE(eke,ekef,lelke,elke,elkef,tstep);
											// ###REPLACE
											// {$COMPUTE-DRANGE(#,#,#,#,#,#);}
											// WITH
											// {
											// "write(6,*) ' COMPUTE-DRANGE: ',{P1},{P2},{P3},{P4},{P5};"
											fedep = 1.0 - ekef / EGS4.EKE;// 1 -
																			// {P2}/{P1};
											// elktmp =
											// 0.5*({P4}+{P5}+0.25*fedep*fedep*(1+fedep*(1+0.875*fedep)));
											elktmp = 0.5 * (EGS4.ELKE + elkef + 0.25
													* fedep
													* fedep
													* (1.0 + fedep
															* (1.0 + 0.875 * fedep)));
											// " the above evaluates the logarithm of the midpoint energy"
											// "write(6,*) ' COMPUTE-DRANGE: fedep elktmp = ',fedep,elktmp;"
											lelktmp = lelke;// lelktmp = {P3};
											if (lelec < 0) {
												// $EVALUATE dedxmid USING
												// ededx(elktmp);
												dedxmid = EGS4.EDEDX1[lelktmp - 1][EGS4.MEDIUM - 1]
														* elktmp
														+ EGS4.EDEDX0[lelktmp - 1][EGS4.MEDIUM - 1];

												aux = EGS4.EDEDX1[lelktmp - 1][EGS4.MEDIUM - 1]
														/ dedxmid;
												// (lelktmp,medium)/dedxmid;
											} else {
												// $EVALUATE dedxmid USING
												// pdedx(elktmp);
												dedxmid = EGS4.PDEDX1[lelktmp - 1][EGS4.MEDIUM - 1]
														* elktmp
														+ EGS4.PDEDX0[lelktmp - 1][EGS4.MEDIUM - 1];

												aux = EGS4.PDEDX1[lelktmp - 1][EGS4.MEDIUM - 1]
														/ dedxmid;
												// (lelktmp,medium)/dedxmid;
											}
											aux = aux * (1.0 + 2.0 * aux)
													* (fedep / (2.0 - fedep))
													* (fedep / (2.0 - fedep))
													/ 6.0;
											// {P6} =
											// fedep*{P1}/dedxmid*(1+aux);
											EGS4.TSTEP = fedep * EGS4.EKE
													/ dedxmid * (1.0 + aux);
											// "write(6,*) ' COMPUTE-DRANGE: aux {P6} = ',aux,{P6};"
											// }####
										} else { // " initial and final energy are in different interpolation bins, "
													// " calc range from ekef to E(lelkef+1) and from E(lelke) to eke  "
													// " and add the pre-calculated range from E(lelkef+1) to E(lelke) "
													// "write(6,*) ' CALCULATE-TSTEP-FROM-DEMFP: diff bins: ',lelkef,lelke;"
											ekei = EGS4.E_array[lelke - 1][EGS4.MEDIUM - 1];// (lelke,medium);
											// "write(6,*) ' CALCULATE-TSTEP-FROM-DEMFP dr1: ',eke,ekei;"
											elkei = (lelke - EGS4.EKE0[EGS4.MEDIUM - 1])
													/ EGS4.EKE1[EGS4.MEDIUM - 1];
											// $COMPUTE-DRANGE(eke,ekei,lelke,elke,elkei,tuss);

											fedep = 1.0 - ekei / EGS4.EKE;// 1 -
																			// {P2}/{P1};
											// elktmp =
											// 0.5*({P4}+{P5}+0.25*fedep*fedep*(1+fedep*(1+0.875*fedep)));
											elktmp = 0.5 * (EGS4.ELKE + elkei + 0.25
													* fedep
													* fedep
													* (1.0 + fedep
															* (1.0 + 0.875 * fedep)));
											// " the above evaluates the logarithm of the midpoint energy"
											// "write(6,*) ' COMPUTE-DRANGE: fedep elktmp = ',fedep,elktmp;"
											lelktmp = lelke;// lelktmp = {P3};
											if (lelec < 0) {
												// $EVALUATE dedxmid USING
												// ededx(elktmp);
												dedxmid = EGS4.EDEDX1[lelktmp - 1][EGS4.MEDIUM - 1]
														* elktmp
														+ EGS4.EDEDX0[lelktmp - 1][EGS4.MEDIUM - 1];

												aux = EGS4.EDEDX1[lelktmp - 1][EGS4.MEDIUM - 1]
														/ dedxmid;
												// (lelktmp,medium)/dedxmid;
											} else {
												// $EVALUATE dedxmid USING
												// pdedx(elktmp);
												dedxmid = EGS4.PDEDX1[lelktmp - 1][EGS4.MEDIUM - 1]
														* elktmp
														+ EGS4.PDEDX0[lelktmp - 1][EGS4.MEDIUM - 1];

												aux = EGS4.PDEDX1[lelktmp - 1][EGS4.MEDIUM - 1]
														/ dedxmid;
												// (lelktmp,medium)/dedxmid;
											}
											aux = aux * (1.0 + 2.0 * aux)
													* (fedep / (2.0 - fedep))
													* (fedep / (2.0 - fedep))
													/ 6.0;
											// {P6} =
											// fedep*{P1}/dedxmid*(1+aux);
											tuss = fedep * EGS4.EKE / dedxmid
													* (1.0 + aux);
											// "write(6,*) ' COMPUTE-DRANGE: aux {P6} = ',aux,{P6};"

											ekei = EGS4.E_array[lelkef][EGS4.MEDIUM - 1];// (lelkef+1,medium);
											// "write(6,*) ' CALCULATE-TSTEP-FROM-DEMFP dr2: ',ekei,ekef;"
											elkei = (lelkef + 1.0 - EGS4.EKE0[EGS4.MEDIUM - 1])
													/ EGS4.EKE1[EGS4.MEDIUM - 1];
											// $COMPUTE-DRANGE(ekei,ekef,lelkef,elkei,elkef,tstep);

											fedep = 1.0 - ekef / ekei;// 1 -
																		// {P2}/{P1};
											// elktmp =
											// 0.5*({P4}+{P5}+0.25*fedep*fedep*(1+fedep*(1+0.875*fedep)));
											elktmp = 0.5 * (elkei + elkef + 0.25
													* fedep
													* fedep
													* (1.0 + fedep
															* (1.0 + 0.875 * fedep)));
											// " the above evaluates the logarithm of the midpoint energy"
											// "write(6,*) ' COMPUTE-DRANGE: fedep elktmp = ',fedep,elktmp;"
											lelktmp = lelkef;// lelktmp = {P3};
											if (lelec < 0) {
												// $EVALUATE dedxmid USING
												// ededx(elktmp);
												dedxmid = EGS4.EDEDX1[lelktmp - 1][EGS4.MEDIUM - 1]
														* elktmp
														+ EGS4.EDEDX0[lelktmp - 1][EGS4.MEDIUM - 1];

												aux = EGS4.EDEDX1[lelktmp - 1][EGS4.MEDIUM - 1]
														/ dedxmid;
												// (lelktmp,medium)/dedxmid;
											} else {
												// $EVALUATE dedxmid USING
												// pdedx(elktmp);
												dedxmid = EGS4.PDEDX1[lelktmp - 1][EGS4.MEDIUM - 1]
														* elktmp
														+ EGS4.PDEDX0[lelktmp - 1][EGS4.MEDIUM - 1];

												aux = EGS4.PDEDX1[lelktmp - 1][EGS4.MEDIUM - 1]
														/ dedxmid;
												// (lelktmp,medium)/dedxmid;
											}
											aux = aux * (1.0 + 2.0 * aux)
													* (fedep / (2.0 - fedep))
													* (fedep / (2.0 - fedep))
													/ 6.0;
											// {P6} =
											// fedep*{P1}/dedxmid*(1+aux);
											EGS4.TSTEP = fedep * ekei / dedxmid
													* (1.0 + aux);
											// "write(6,*) ' COMPUTE-DRANGE: aux {P6} = ',aux,{P6};"

											EGS4.TSTEP = EGS4.TSTEP
													+ tuss
													+ // tstep=tstep+tuss+
													EGS4.range_ep[qel][lelke - 1][EGS4.MEDIUM - 1]
													- // (qel,lelke,medium)-
													EGS4.range_ep[qel][lelkef][EGS4.MEDIUM - 1];// (qel,lelkef+1,medium);
										}
										// "write(6,*) ' CALCULATE-TSTEP-FROM-DEMFP: tstep = ',tstep;"
									}
									total_tstep = EGS4.TSTEP;// total_tstep =
																// tstep;
									compute_tstep = false;
								}
								EGS4.TSTEP = total_tstep / EGS4.RHOF; // " non-default density scaling "
								// tstep = total_tstep/rhof;
								// }

							}// "end sig if-else

							// "calculate stopping power"
							if (lelec < 0) {
								// $EVALUATE dedx0 USING ededx(elke);
								dedx0 = EGS4.EDEDX1[lelke - 1][EGS4.MEDIUM - 1]
										* EGS4.ELKE
										+ EGS4.EDEDX0[lelke - 1][EGS4.MEDIUM - 1];
							} // "e-"
							else {
								// $EVALUATE dedx0 USING pdedx(elke);
								dedx0 = EGS4.PDEDX1[lelke - 1][EGS4.MEDIUM - 1]
										* EGS4.ELKE
										+ EGS4.PDEDX0[lelke - 1][EGS4.MEDIUM - 1];
							}// "e+"
							dedx = EGS4.RHOF * dedx0;

							// "Determine maximum step-size (Formerly
							// $SET-TUSTEP)
							// $EVALUATE tmxs USING tmxs(elke);
							tmxs = EGS4.TMXS1[lelke - 1][EGS4.MEDIUM - 1]
									* EGS4.ELKE
									+ EGS4.TMXS0[lelke - 1][EGS4.MEDIUM - 1];

							tmxs = tmxs / EGS4.RHOF;

							// "Compute the range to E_min(medium) (e_min is the
							// first
							// "energy in the table). Do not go more than range.
							// "Don't replace this macro and don't override
							// range, because
							// "the energy loss evaluation below relies on the
							// accurate
							// "(and self-consistent) evaluation of range!
							// $COMPUTE-RANGE;
							// REPLACE {$COMPUTE-RANGE;} WITH
							// {
							ekei = EGS4.E_array[lelke - 1][EGS4.MEDIUM - 1];// (lelke,medium);
							elkei = (lelke - EGS4.EKE0[EGS4.MEDIUM - 1])
									/ EGS4.EKE1[EGS4.MEDIUM - 1];
							// $COMPUTE-DRANGE(eke,ekei,lelke,elke,elkei,range);

							fedep = 1.0 - ekei / EGS4.EKE;// 1 - {P2}/{P1};
							// elktmp =
							// 0.5*({P4}+{P5}+0.25*fedep*fedep*(1+fedep*(1+0.875*fedep)));
							elktmp = 0.5 * (EGS4.ELKE + elkei + 0.25 * fedep
									* fedep
									* (1.0 + fedep * (1.0 + 0.875 * fedep)));
							// " the above evaluates the logarithm of the midpoint energy"
							// "write(6,*) ' COMPUTE-DRANGE: fedep elktmp = ',fedep,elktmp;"
							lelktmp = lelke;// lelktmp = {P3};
							if (lelec < 0) {
								// $EVALUATE dedxmid USING ededx(elktmp);
								dedxmid = EGS4.EDEDX1[lelktmp - 1][EGS4.MEDIUM - 1]
										* elktmp
										+ EGS4.EDEDX0[lelktmp - 1][EGS4.MEDIUM - 1];

								aux = EGS4.EDEDX1[lelktmp - 1][EGS4.MEDIUM - 1]
										/ dedxmid;
								// (lelktmp,medium)/dedxmid;
							} else {
								// $EVALUATE dedxmid USING pdedx(elktmp);
								dedxmid = EGS4.PDEDX1[lelktmp - 1][EGS4.MEDIUM - 1]
										* elktmp
										+ EGS4.PDEDX0[lelktmp - 1][EGS4.MEDIUM - 1];

								aux = EGS4.PDEDX1[lelktmp - 1][EGS4.MEDIUM - 1]
										/ dedxmid;
								// (lelktmp,medium)/dedxmid;
							}
							aux = aux * (1.0 + 2.0 * aux)
									* (fedep / (2.0 - fedep))
									* (fedep / (2.0 - fedep)) / 6.0;
							// {P6} = fedep*{P1}/dedxmid*(1+aux);
							range = fedep * EGS4.EKE / dedxmid * (1.0 + aux);
							// "write(6,*) ' COMPUTE-DRANGE: aux {P6} = ',aux,{P6};"

							range = (range + EGS4.range_ep[qel][lelke - 1][EGS4.MEDIUM - 1])
									/ EGS4.RHOF;// (qel,lelke,medium))/rhof;
							// }

							// "The RANDOMIZE-TUSTEP option as coded by AFB
							// forced the
							// "electrons to approach discrete events
							// (Moller,brems etc.)
							// "only in a single scattering mode => waste of CPU
							// time.
							// "Moved here and changed by IK Oct 22 1997
							random_tustep = $RANDOMIZE_TUSTEP;
							if (random_tustep) {
								rnnotu = EGS4.random01();
								tmxs = rnnotu
										* Math.min(tmxs, EGS4.SMAXIR[irl - 1]);
							} else {
								tmxs = Math.min(tmxs, EGS4.SMAXIR[irl - 1]);
							}
							EGS4.TUSTEP = EGS4.min(EGS4.TSTEP, tmxs, range);// tustep
																			// =
																			// min(tstep,tmxs,range);
							// ###############################
							// $SET-TUSTEP-EM-FIELD; //"optional tustep
							// restriction in EM field
							// REPLACE {$SET-TUSTEP-EM-FIELD;} WITH {;}-->more
							// ABSTRACT METHODS IN TIME
							if (EMFb) {
								// "THIS MACRO SETS LIMITS ON THE ELECTRON STEP-SIZE FOR EM-FIELD ELECTRON       "
								// "TRANSPORT. THE VARIOUS CONSTRAINTS ARE PUT IN MACRO FORM SO THAT THEY MAY    "
								// "BE EXCLUDED BY THE USER. THEIR EVEALUATION MAY BE QUITE TIME CONSUMING.      "
								// "IF THE FIELDS ARE WEAK ENOUGH THEY SHOULDN'T AFFECT THE TRANSPORT.           "
								// REPLACE {$SET-TUSTEP-EM-FIELD;} WITH {;
								//tustp0 = EGS4.TUSTEP;
								x0 = EGS4.X[EGS4.NP - 1];
								y0 = EGS4.Y[EGS4.NP - 1];
								z0 = EGS4.Z[EGS4.NP - 1];
								EGS4emf.GET_EM_FIELD("E", "B", x0, y0, z0, iemf);
								// -----------------------------------
								Ex0 = EGS4emf.Ex0;
								Ey0 = EGS4emf.Ey0;
								Ez0 = EGS4emf.Ez0;
								Bx0 = EGS4emf.Bx0;
								By0 = EGS4emf.By0;
								Bz0 = EGS4emf.Bz0;
								// -----------------------------------
								ekin0 = eie - EGS4.RM;
								beta20 = Math.max(1.E-8, ekin0
										* (eie + EGS4.RM) / (eie * eie));
								beta0 = Math.sqrt(beta20);
								gamma0 = 1. / Math.sqrt((1. + beta0)
										* (1. - beta0));
								fnorm = lelec / (beta20 * gamma0);
								Ex0 = fnorm * Ex0;
								Ey0 = fnorm * Ey0;
								Ez0 = fnorm * Ez0;
								fnormb = beta0 * fnorm;
								Bx0 = fnormb * Bx0;
								By0 = fnormb * By0;
								Bz0 = fnormb * Bz0;
								u0 = EGS4.U[EGS4.NP - 1];
								v0 = EGS4.V[EGS4.NP - 1];
								w0 = EGS4.W[EGS4.NP - 1];
								edotu = Ex0 * u0 + Ey0 * v0 + Ez0 * w0;
								dperpx = Ex0 - u0 * edotu;
								dperpy = Ey0 - v0 * edotu;
								dperpz = Ez0 - w0 * edotu;
								dperp2 = dperpx * dperpx + dperpy * dperpy
										+ dperpz * dperpz;
								bperpx = v0 * Bz0 - w0 * By0;
								bperpy = w0 * Bx0 - u0 * Bz0;
								bperpz = u0 * By0 - v0 * Bx0;
								bperp2 = bperpx * bperpx + bperpy * bperpy
										+ bperpz * bperpz;
								eperpx = dperpx + bperpx;
								eperpy = dperpy + bperpy;
								eperpz = dperpz + bperpz;
								eperp2 = eperpx * eperpx + eperpy * eperpy
										+ eperpz * eperpz;
								// $SET-TUSTEP-DIRECTION-VECTOR-CHANGE-EM-FIELD;
								// "THESE MACROS PERFORM THE ACTUAL STEP-SIZE SHORTENING TO SATISFY THE          "
								// "VALIDITY CRITERIA OF TRANSPORT IN EXTERNAL EM FIELDS                         "
								// "                                                                             "
								// REPLACE
								// {$SET-TUSTEP-DIRECTION-VECTOR-CHANGE-EM-FIELD;}
								// WITH {;
								if (eperp2 != 0.0) {
									eperp = Math.sqrt(eperp2);
									EGS4.TUSTEP = Math.min(EGS4.TUSTEP, $EMULMT
											/ eperp);
								}
								// }
								// $SET-TUSTEP-ENERGY-CHANGE-EM-FIELD;
								// REPLACE {$SET-TUSTEP-ENERGY-CHANGE-EM-FIELD;}
								// WITH {;
								if (edotu != 0.0)
									EGS4.TUSTEP = Math.min(
											EGS4.TUSTEP,
											Math.abs($EMELMT
													* (1. + 1. / gamma0)
													/ edotu));
								// }
								// $SET-TUSTEP-CHANGE-OF-EM-FIELD;
								// REPLACE {$SET-TUSTEP-CHANGE-OF-EM-FIELD;}
								// WITH {;
								if (EGS4.TUSTEP != 0.0) {
									xf = x0 + u0 * EGS4.TUSTEP;
									yf = y0 + v0 * EGS4.TUSTEP;
									zf = z0 + w0 * EGS4.TUSTEP;
									EGS4emf.GET_EM_FIELD("ef", "bf", xf, yf,
											zf, iemf);
									// -----------------------------------
									efx0 = EGS4emf.efx0;
									efy0 = EGS4emf.efy0;
									efz0 = EGS4emf.efz0;
									//bfx0 = EGS4emf.bfx0;
									//bfy0 = EGS4emf.bfy0;
									//bfz0 = EGS4emf.bfz0;
									// -----------------------------------
									test = Ex0 * Ex0 + Ey0 * Ey0 + Ez0 * Ez0;
									if (test != 0.0) {
										// test=((EFXF-EX0)**2+(EFYF-EY0)**2+(EFZF-EZ0)**2)/TEST;
										test = ((efx0 - Ex0) * (efx0 - Ex0)
												+ (efy0 - Ey0) * (efy0 - Ey0) + (efz0 - Ez0)
												* (efz0 - Ez0))
												/ test;
										if ((test != 0.0) && (test < 1.0))
											EGS4.TUSTEP = Math.min(EGS4.TUSTEP,
													$EMFLMT * EGS4.TUSTEP
															/ Math.sqrt(test));
									}
								}
								// }
								// $SET-TUSTEP-DIRECTION-VECTOR-CHANGE-BY-MULTIPLE-SCATTERING;
								// Replace with null for now, IK, July 2004.
								// / *
								// Replace with null for now, IK, July 2004.
								// REPLACE
								// {$SET-TUSTEP-DIRECTION-VECTOR-CHANGE-BY-MULTIPLE-SCATTERING;}
								// WITH {;
								// IF(MEDIUM.NE.0)[
								// AMSPLC=BLCC(MEDIUM)/BETA20;
								// OMEGA0=AMSPLC*TUSTEP*RHOF;
								// IF(OMEGA0.LE.2.718282)[B=1;]
								// ELSE[
								// BLC=ALOG(OMEGA0);
								// IF(BLC.LT.1.306853)[B=-10.27666+BLC*(17.82596-6.468813*BLC);]
								// ELSE[IB=B0BGB+BLC*B1BGB;
								// IF(IB.GT.NBGB)[OUTPUT IB;(' NBGB[IB=',I5);]
								// B=BGB0(IB)+BLC*(BGB1(IB)+BLC*BGB2(IB));
								// ]
								// ]
								// GMSPLC=RHOF*(XCC(MEDIUM)/(EIE*BETA20))**2;
								// XCC2BT=GMSPLC*B;
								// IF(XCC2BT.NE.0.0)
								// TUSTEP=AMIN1(TUSTEP,$EMMLMT/XCC2BT);
								// ]
								// }
								// * /
								// REPLACE
								// {$SET-TUSTEP-DIRECTION-VECTOR-CHANGE-BY-MULTIPLE-SCATTERING;}
								// WITH {;}
								// }
							}// if(EMFb)
								// ##########################################
								// $CALL-eq.HOWNEAR(tperp);
								// ---------------------------------
							EGS4.irl = irl;// pass the irl
							eq.HOWNEAR();
							tperp = EGS4.tperp;
							// ---------------------------------
							EGS4.DNEAR[EGS4.NP - 1] = tperp;

							// $RANGE-DISCARD;
							// //"optional regional range rejection for"
							// "particles below e_max_rr if i_do_rr set"
							// REPLACE {$RANGE-DISCARD;} WITH {
							if ((EGS4.i_do_rr[irl - 1] == 1)
									&& (EGS4.E[EGS4.NP - 1] < EGS4.e_max_rr[irl - 1])) {
								if (tperp >= range) {// "particle cannot escape local region"
									EGS4.IDISC = 50 + 49 * EGS4.IQ[EGS4.NP - 1];
									// "1 for electrons, 99 for positrons"
									// go to :USER-ELECTRON-DISCARD: ;
									USERb = true;
									break BREMSLOOP;// NEWELECTRON;
								}
							}
							// };
							// $USER-RANGE-DISCARD;
							// //"default is ;, but user may implement"
							// pass the var========================
							EGS4Macro.range = range;
							EGS4Macro.irl = irl;
							// end pass var========================
							if (EGS4Macro.USER_RANGE_DISCARD())// the abstract
																// function can
																// return true
																// or false;
							{
								USERb = true;
								break BREMSLOOP;// NEWELECTRON;
							}
							// getvar=======================

							// end getvar=====================

							// $SET-SKINDEPTH(eke,elke);
							// "This macro sets the minimum step size for a condensed"
							// "history (CH) step. When the exact BCA is used, the minimum"
							// "CH step is determined by efficiency considerations only"
							// "At about 3 elastic MFP's single scattering becomes more"
							// "efficient than CH and so the algorithm switches off CH"
							// "If one of the various inexact BCA's is invoked, this macro"
							// "provides a simple way to include more sophisticated"
							// "decisions about the maximum acceptable approximated CH step"
							// REPLACE {$SET-SKINDEPTH(#,#);} WITH
							// {
							// $CALCULATE-ELASTIC-SCATTERING-MFP(ssmfp,{P1},{P2});//eke,elke
							// "This macro calculates the elastic scattering MFP"
							// "If spin_effects is .false., the screened Rutherford cross section"
							// "is used, else the the elastic MFP is based on PWA cross sections"
							// REPLACE
							// {$CALCULATE-ELASTIC-SCATTERING-MFP(#,#,#);} WITH
							// {
							blccl = EGS4.RHOF * EGS4.BLCC[EGS4.MEDIUM - 1];
							xccl = EGS4.RHOF * EGS4.XCC[EGS4.MEDIUM - 1];
							// p2 = {P2}*({P2}+rmt2);
							p2 = EGS4.EKE * (EGS4.EKE + EGS4.RMT2);
							beta2 = p2 / (p2 + EGS4.RMSQ);
							if (EGS4.spin_effects) {
								if (lelec < 0) {
									// $EVALUATE etap USING etae_ms({P3});
									etap = EGS4.etae_ms1[lelke - 1][EGS4.MEDIUM - 1]
											* EGS4.ELKE
											+ EGS4.etae_ms0[lelke - 1][EGS4.MEDIUM - 1];
								} else {
									// $EVALUATE etap USING etap_ms({P3});
									etap = EGS4.etap_ms1[lelke - 1][EGS4.MEDIUM - 1]
											* EGS4.ELKE
											+ EGS4.etap_ms0[lelke - 1][EGS4.MEDIUM - 1];
								}
								// $EVALUATE ms_corr USING blcce({P3});//elke
								// here
								ms_corr = EGS4.blcce1[lelke - 1][EGS4.MEDIUM - 1]
										* EGS4.ELKE
										+ EGS4.blcce0[lelke - 1][EGS4.MEDIUM - 1];

								blccl = blccl
										/ etap
										/ (1.0 + 0.25 * etap * xccl / blccl
												/ p2) * ms_corr;
							}
							ssmfp = beta2 / blccl;// {P1}=beta2/blccl;
							// }

							skindepth = EGS4.skindepth_for_bca * ssmfp;
							// }

							EGS4.TUSTEP = Math.min(EGS4.TUSTEP,
									Math.max(tperp, skindepth));
							// tustep = min(tustep,max(tperp,skindepth));
							// "The transport logic below is determined by the
							// logical
							// "variables callheq.HOWFAR, domultiple and
							// dosingle
							// "
							// "There are the following possibilities:
							// "
							// " calleq.HOWFAR = .false. This indicates that the
							// " ==================== intended step is shorter
							// than tperp
							// " independent of BCA used
							// " - domultiple = .false. dosingle = .false. and
							// " callmsdist = .true.
							// " ==> everything has been done in msdist
							// " - domultiple = .true. and dosingle = .false.
							// " ==> should happen only if exact_bca = .false.
							// " indicates that MS remains to be done
							// " - domultiple = .false. and dosingle = .true.
							// " ==> should happen only if exact_bca = .true.
							// " sampled distance to a single scattering event
							// is
							// / " shorter than tperp ==> do single scattering
							// at the
							// " end of the step
							// " - domultiple = .true. and dosingle = .true.
							// " ==> error condition, something with the logic
							// is wrong!
							// "
							// " calleq.HOWFAR = .true. This indicates that the
							// intended step
							// " =================== is longer than tperp and
							// forces a
							// " call to hawfar which returns the
							// " straight line distance to the boundary
							// " in the initial direction of motion
							// " (via a modification of ustep)
							// " - domultiple = .false. and dosingle = .false.
							// " ==> should happen only of exact_bca=.true.
							// " simply put the particle on the boundary
							// " - domultiple = .false. and dosingle = .true.
							// " ==> should happen only of exact_bca=.true.
							// " single elastic scattering has to be done
							// " - domultiple = .true. and dosingle = .false.
							// " ==> should happen only of exact_bca=.false.
							// " indicates that MS remains to be done
							// " - domultiple = .true. and dosingle = .true.
							// " ==> error condition, something with the logic
							// is wrong!

							// "IF(tustep <= tperp & tustep > skindepth)"
							// "This statement changed to be consistent with PRESTA-I"
							EGS4.count_all_steps = EGS4.count_all_steps + 1;
							EGS4.is_ch_step = false;
							if ((EGS4.TUSTEP <= tperp)
									&& ((!EGS4.exact_bca) || (EGS4.TUSTEP > skindepth))) {
								// "We are further way from a boundary than a
								// skindepth, so
								// "perform a normal condensed-history step
								callHOWFAR = false; // "Do not call HAWFAR
								domultiple = false; // "Multiple scattering done
													// here
								dosingle = false; // "MS => no single scattering
								callmsdist = true; // "Remember that msdist has
													// been called

								// "Fourth order technique for de
								// $COMPUTE-ELOSS-G(tustep,eke,elke,lelke,de);
								// REPLACE {$COMPUTE-ELOSS-G(#,#,#,#,#);} WITH
								// {
								tuss = range
										- EGS4.range_ep[qel][lelke - 1][EGS4.MEDIUM - 1]
										/ EGS4.RHOF;
								// (qel,{P4},medium)/rhof;
								// " here tuss is the range between the initial energy and the next lower "
								// " energy on the interpolation grid "
								// if( tuss >= {P1} )
								if (tuss >= EGS4.TUSTEP) {// " Final energy is in the same interpolation bin "
															// $COMPUTE-ELOSS({P1},{P2},{P3},{P4},{P5});
															// HERE(tustep,eke,elke,lelke,de)the
															// same;
															// REPLACE
															// {$COMPUTE-ELOSS(#,#,#,#,#);}
															// WITH
															// {;
									if (lelec < 0) {
										// $EVALUATE dedxmid USING ededx({P3});
										dedxmid = EGS4.EDEDX1[lelke - 1][EGS4.MEDIUM - 1]
												* EGS4.ELKE
												+ EGS4.EDEDX0[lelke - 1][EGS4.MEDIUM - 1];

										// aux = ededx1({P4},medium)/dedxmid;
										aux = EGS4.EDEDX1[lelke - 1][EGS4.MEDIUM - 1]
												/ dedxmid;
									} else {
										// $EVALUATE dedxmid USING pdedx({P3});
										dedxmid = EGS4.PDEDX1[lelke - 1][EGS4.MEDIUM - 1]
												* EGS4.ELKE
												+ EGS4.PDEDX0[lelke - 1][EGS4.MEDIUM - 1];

										// aux = pdedx1({P4},medium)/dedxmid;
										aux = EGS4.PDEDX1[lelke - 1][EGS4.MEDIUM - 1]
												/ dedxmid;
									}
									// {P5} = dedxmid*{P1};
									// //" Energy loss using stopping power at the beginning "
									// de = dedxmid*EGS4.TUSTEP;
									de = dedxmid * EGS4.TUSTEP * EGS4.RHOF;
									/*
									 * {P5} = dedxmid*{P1}*rhof;
									 * "IK: rhof scaling bug, June 9 2006"
									 * "rhof scaling must be done here and NOT in "
									 * "$COMPUTE-ELOSS-G below!"
									 */
									fedep = de / EGS4.EKE;// fedep = {P5}/{P2};
									// {P5} =
									// {P5}*(1-0.5*fedep*aux*(1-0.333333*fedep*(aux-1-
									de = de
											* (1.0 - 0.5
													* fedep
													* aux
													* (1.0 - 0.333333
															* fedep
															* (aux - 1.0 - 0.25
																	* fedep
																	* (2.0 - aux
																			* (4.0 - aux)))));
									// }

									// {P5} = {P5}*rhof; //"IK, rhof bug"
									// de =
									// de*EGS4.RHOF;//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
									// "IK: rhof scaling bug, June 9 2006. rhof scaling is done in "
								} else {// " Must find first the table index where the step ends using "
										// " pre-calculated ranges                                     "
									lelktmp = lelke;// lelktmp = {P4};
									// tuss = (range - {P1})*rhof;
									tuss = (range - EGS4.TUSTEP) * EGS4.RHOF;
									// " now tuss is the range of the final energy electron "
									// " scaled to the default mass density from PEGS4      "
									if (tuss <= 0.0) {
										// {P5} = {P2};- TE(medium)*0.99;
										de = EGS4.EKE
												- EGS4.TE[EGS4.MEDIUM - 1]
												* 0.99;
									}
									// " i.e., if the step we intend to take is longer than the particle "
									// //" range, the particle loses its entire energy (which is in {P2})  "
									// " range, the particle energy goes down to the threshold "
									// "({P2} is the initial particle energy)  "
									// "originally the entire energy was lost, but msdist_xxx is not prepared"
									// "to deal with such large eloss fractions => changed July 2005."
									else {
										while (tuss < EGS4.range_ep[qel][lelktmp - 1][EGS4.MEDIUM - 1])// (qel,lelktmp,medium)
																										// )
										{
											lelktmp = lelktmp - 1;
										}
										elktmp = (lelktmp + 1 - EGS4.EKE0[EGS4.MEDIUM - 1])
												/ EGS4.EKE1[EGS4.MEDIUM - 1];
										eketmp = EGS4.E_array[lelktmp][EGS4.MEDIUM - 1];// (lelktmp+1,medium);
										// tuss =
										// EGS4.range_ep[qel][lelktmp][EGS4.MEDIUM-1]
										// - tuss;//(qel,lelktmp+1,medium) -
										// tuss;
										/*
										 * "IK: rhof scaling bug, June 9 2006: because of the change in "
										 * "    $COMPUTE-ELOSS above, we must scale tuss by rhof        "
										 * tuss =
										 * (range_ep(qel,lelktmp+1,medium) -
										 * tuss)/rhof;
										 */
										tuss = (EGS4.range_ep[qel][lelktmp][EGS4.MEDIUM - 1] - tuss)
												/ EGS4.RHOF;
										// $COMPUTE-ELOSS(tuss,eketmp,elktmp,lelktmp,{P5});

										// REPLACE {$COMPUTE-ELOSS(#,#,#,#,#);}
										// WITH
										// {;
										if (lelec < 0) {
											// $EVALUATE dedxmid USING
											// ededx({P3});
											dedxmid = EGS4.EDEDX1[lelktmp - 1][EGS4.MEDIUM - 1]
													* elktmp
													+ EGS4.EDEDX0[lelktmp - 1][EGS4.MEDIUM - 1];

											// aux =
											// ededx1({P4},medium)/dedxmid;
											aux = EGS4.EDEDX1[lelktmp - 1][EGS4.MEDIUM - 1]
													/ dedxmid;
										} else {
											// $EVALUATE dedxmid USING
											// pdedx({P3});
											dedxmid = EGS4.PDEDX1[lelktmp - 1][EGS4.MEDIUM - 1]
													* elktmp
													+ EGS4.PDEDX0[lelktmp - 1][EGS4.MEDIUM - 1];

											// aux =
											// pdedx1({P4},medium)/dedxmid;
											aux = EGS4.PDEDX1[lelktmp - 1][EGS4.MEDIUM - 1]
													/ dedxmid;
										}
										// {P5} = dedxmid*{P1};
										// " Energy loss using stopping power at the beginning "
										// de = dedxmid*tuss;
										de = dedxmid * tuss * EGS4.RHOF;
										fedep = de / eketmp;// fedep =
															// {P5}/{P2};
										de = de
												* (1.0 - 0.5
														* fedep
														* aux
														* (1.0 - 0.333333
																* fedep
																* (aux - 1.0 - 0.25
																		* fedep
																		* (2.0 - aux
																				* (4.0 - aux)))));
										// }

										// {P5} = {P5} + {P2} - eketmp;
										de = de + EGS4.EKE - eketmp;
									}
								}
								// }

								EGS4.TVSTEP = EGS4.TUSTEP;// tvstep = tustep;
								EGS4.is_ch_step = true;
								if (EGS4.transport_algorithm == EGS4.$PRESTA_II) {
									msdist2call = 0;// added
									msdist_pII(
											// "Inputs
											EGS4.EKE, de, EGS4.TUSTEP,
											EGS4.RHOF, EGS4.MEDIUM, qel,
											EGS4.spin_effects,
											EGS4.U[EGS4.NP - 1],
											EGS4.V[EGS4.NP - 1],
											EGS4.W[EGS4.NP - 1],
											EGS4.X[EGS4.NP - 1],
											EGS4.Y[EGS4.NP - 1],
											EGS4.Z[EGS4.NP - 1]
									// "Outputs
									// ,uscat,vscat,wscat,xtrans,ytrans,ztrans,ustep
									);
								} else {
									msdist1call = 0;// added
									msdist_pI(
											// "Inputs
											EGS4.EKE, de, EGS4.TUSTEP,
											EGS4.RHOF, EGS4.MEDIUM, qel,
											EGS4.spin_effects,
											EGS4.U[EGS4.NP - 1],
											EGS4.V[EGS4.NP - 1],
											EGS4.W[EGS4.NP - 1],
											EGS4.X[EGS4.NP - 1],
											EGS4.Y[EGS4.NP - 1],
											EGS4.Z[EGS4.NP - 1]
									// "Outputs
									// ,uscat,vscat,wscat,xtrans,ytrans,ztrans,ustep
									);
								}
							} else {
								// "We are within a skindepth from a boundary,
								// invoke
								// "one of the various boundary-crossing
								// algorithms
								callmsdist = false;
								// "Remember that msdist has not been called
								if (EGS4.exact_bca) {
									// "Cross the boundary in a single
									// scattering mode
									domultiple = false; // "Do not do multiple
														// scattering
									// "Sample the distance to a single
									// scattering event
									rnnoss = EGS4.random01();
									lambda = -Math.log(1.0 - rnnoss);
									lambda_max = 0.5 * blccl * EGS4.RM / dedx
											* (EGS4.EKE / EGS4.RM + 1.0)
											* (EGS4.EKE / EGS4.RM + 1.0)
											* (EGS4.EKE / EGS4.RM + 1.0);
									if ((lambda >= 0.0) && (lambda_max > 0.0)) {
										if (lambda < lambda_max) {
											tuss = lambda
													* ssmfp
													* (1.0 - 0.5 * lambda
															/ lambda_max);
										} else {
											tuss = 0.5 * lambda * ssmfp;
										}

										if (tuss < EGS4.TUSTEP) {
											EGS4.TUSTEP = tuss;
											dosingle = true;
										} else {
											dosingle = false;
										}
									} else {
										// write(6,*) ' lambda > lambda_max:
										// ',lambda,lambda_max;
										// write(6,*) ' eke dedx: ',eke,dedx;
										// write(6,*) ' medium blcc:
										// ',medium,blcc(medium);
										dosingle = false;
										EGS4.NP = EGS4.NP - 1;
										return;
									}
									EGS4.USTEP = EGS4.TUSTEP;// ustep = tustep;
								} else {
									// "Boundary crossing a la EGS4/PRESTA-I but
									// using
									// "exact PLC
									dosingle = false;
									domultiple = true;
									// $SET-USTEP;
									// REPLACE {$SET-USTEP;} WITH
									// {
									ekems = EGS4.EKE - 0.5 * EGS4.TUSTEP * dedx; // "Use mid-point energy to calculate"
																					// "energy dependent quantities"
									// $CALCULATE-XI(tustep);
									// REPLACE {$CALCULATE-XI(#);} WITH
									// {
									p2 = ekems * (ekems + EGS4.RMT2);
									beta2 = p2 / (p2 + EGS4.RMSQ);
									chia2 = xccl / (4.0 * blccl * p2);
									// "Note that our chia2 is Moliere chia2/4"
									// "Note also that xcc is now old egs xcc**2"
									// xi = 0.5*xccl/p2/beta2*{P1};
									xi = 0.5 * xccl / p2 / beta2 * EGS4.TUSTEP;
									if (EGS4.spin_effects) {
										elkems = Math.log(ekems);
										// $SET INTERVAL elkems,eke;
										Double dbl = new Double(
												EGS4.EKE1[EGS4.MEDIUM - 1]
														* elkems
														+ EGS4.EKE0[EGS4.MEDIUM - 1]);
										lelkems = dbl.intValue();

										if (lelec < 0) {
											// $EVALUATE etap USING
											// etae_ms(elkems);
											etap = EGS4.etae_ms1[lelkems - 1][EGS4.MEDIUM - 1]
													* elkems
													+ EGS4.etae_ms0[lelkems - 1][EGS4.MEDIUM - 1];

											// $EVALUATE xi_corr USING
											// q1ce_ms(elkems);
											xi_corr = EGS4.q1ce_ms1[lelkems - 1][EGS4.MEDIUM - 1]
													* elkems
													+ EGS4.q1ce_ms0[lelkems - 1][EGS4.MEDIUM - 1];

										} else {
											// $EVALUATE etap USING
											// etap_ms(elkems);
											etap = EGS4.etap_ms1[lelkems - 1][EGS4.MEDIUM - 1]
													* elkems
													+ EGS4.etap_ms0[lelkems - 1][EGS4.MEDIUM - 1];

											// $EVALUATE xi_corr USING
											// q1cp_ms(elkems);
											xi_corr = EGS4.q1cp_ms1[lelkems - 1][EGS4.MEDIUM - 1]
													* elkems
													+ EGS4.q1cp_ms0[lelkems - 1][EGS4.MEDIUM - 1];

										}
										chia2 = chia2 * etap;
										xi = xi * xi_corr;
										// $EVALUATE ms_corr USING
										// blcce(elkems);
										ms_corr = EGS4.blcce1[lelkems - 1][EGS4.MEDIUM - 1]
												* elkems
												+ EGS4.blcce0[lelkems - 1][EGS4.MEDIUM - 1];

										blccl = blccl * ms_corr;
									} else {
										xi_corr = 1.0;
										etap = 1.0;
									}
									xi = xi
											* (Math.log(1.0 + 1. / chia2) - 1.0 / (1.0 + chia2));
									// }

									if (xi < 0.1) {
										EGS4.USTEP = EGS4.TUSTEP
												* (1.0 - xi
														* (0.5 - xi * 0.166667));
									} else {
										EGS4.USTEP = EGS4.TUSTEP
												* (1.0 - Math.exp(-xi)) / xi;
									}
									// }

								}
								if (EGS4.USTEP < tperp)// ustep < tperp
								{
									callHOWFAR = false;
								} else {
									callHOWFAR = true;
								}
							}

						} // "^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^end
							// non-vacuum test

						// $SET-USTEP-EM-FIELD; //"additional ustep restriction
						// in em field
						// "default for $SET-USTEP-EM-FIELD; is ;(null)
						if (EMFb) {
							// "*****************************************************************************"
							// "                                                                             "
							// "THIS MACRO RESETS THE STRAIGHT LINE DISTANCE TO BE TRANSPORTED BEFORE        "
							// "THE CALL TO eq.HOWFAR AND IT ACCOUNTS FOR THE STRAIGHT LINE LENGTHENING OR      "
							// "TRUNCATION THAT MAY OCCUR IN AN EM FIELD. FOR NOW IT DOES NOT DO ANYTHING    "
							// "ADDITIONAL BECAUSE IT IS A SECOND ORDER EFFECT WHICH HAS NOT YET BEEN        "
							// "DISCOVERED.                                                                  "
							// REPLACE {$SET-USTEP-EM-FIELD;} WITH {;}
							// "IF(MEDIUM.EQ.0)[USTEP=TUSTEP;]"
							// "ELSEIF(TUSTP0.NE.TUSTEP) [$SET-USTEP;]"
							// "}"

						}

						// irnew = ir(np); "default new region is old region
						EGS4.IRNEW = EGS4.IR[EGS4.NP - 1];
						EGS4.IDISC = 0; // "default is no discard (this flag is
										// initialized here)
						ustep0 = EGS4.USTEP;// ustep; "Save the intended ustep."

						// "IF(calleq.HOWFAR) [ call eq.HOWFAR; ]"
						// $CALL-eq.HOWFAR-IN-ELECTR;
						// //"The above is the default replacement"
						// REPLACE {$CALL-eq.HOWFAR-IN-ELECTR;} WITH {;
						if ((callHOWFAR) || (EGS4.WT[EGS4.NP - 1] <= 0.0)) {
							eq.HOWFAR();
						}
						// };

						// "Now see if user requested discard
						if (EGS4.IDISC > 0) // "(idisc is returned by eq.HOWFAR)
						{
							// "User requested immediate discard
							// go to :USER-ELECTRON-DISCARD:;
							USERb = true;
							break BREMSLOOP;// NEWELECTRON;

						}
						// 22222222222222222222222222222222222222222222222222
						// $CHECK-NEGATIVE-USTEP
						// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
						if (EGS4.USTEP <= 0.0)// (ustep <= 0)
						{
							// "Negative ustep---probable truncation problem at
							// a
							// "boundary, which means we are not in the region
							// we think
							// "we are in. The default macro assumes that user
							// has set
							// "irnew to the region we are really most likely to
							// be
							// "in. A message is written out whenever ustep is
							// less than
							// "-1.e-4
							if (EGS4.USTEP < -1.e-4) {
								ierust = ierust + 1;
								// OUTPUT
								// ierust,ustep,dedx,e(np)-prm,
								// ir(np),irnew,irold,x(np),y(np),z(np)
								// ;
								// (
								// i4,' Negative ustep = ',e12.5,' dedx=',F8.4,'
								// ke=',F8.4,
								// ' ir,irnew,irold =',3i4,' x,y,z =',4e10.3
								// );
								EGS4.seqStr = " " + ierust
										+ " Negative ustep = " + EGS4.USTEP
										+ " dedx= " + dedx + "  totalen="
										+ EGS4.E[EGS4.NP - 1] + "  ir="
										+ EGS4.IR[EGS4.NP - 1] + "  irnew="
										+ EGS4.IRNEW + "  irold=" + EGS4.IROLD
										+ "  x=" + EGS4.X[EGS4.NP - 1] + "  y="
										+ EGS4.Y[EGS4.NP - 1] + " z="
										+ EGS4.Z[EGS4.NP - 1];// +" \n";
								if (EGS4.iprint > 2)
									eq.printSequence(EGS4.seqStr);

								if (ierust > 1000) {
									// OUTPUT;(////' Called exit---too many
									// ustep errors'///);
									EGS4.STOPPROGRAM = true;
									EGS4.seqStr = "Called exit---too many ustep errors";
									// if(EGS4.iprint>2)
									eq.printSequence(EGS4.seqStr);

									return;
									// stop;
								}
							}
							EGS4.USTEP = 0.0;
						}

						if ((EGS4.USTEP == 0.0) || (EGS4.MEDIUM == 0)) {
							// "Do fast step
							if (EGS4.USTEP != 0.0) {
								// "Step in vacuum
								EGS4.VSTEP = EGS4.USTEP;// vstep = ustep;
								EGS4.TVSTEP = EGS4.VSTEP;// tvstep = vstep;
								// "( vstep is ustep truncated (possibly) by
								// eq.HOWFAR
								// " tvstep is the total curved path associated
								// with vstep)
								EGS4.EDEP = EGS4.PZERO;// edep = pzero; "no
														// energy loss in vacuum
								// $VACUUM-TRANSPORT-EM-FIELD;
								// "additional vacuum transport in em field
								if (EMFb) {
									// "THIS MACRO IS USED TO SET THE ANGLES AND PATHLENGTH CORRECTION IN VACUUM     "
									// REPLACE {$VACUUM-TRANSPORT-EM-FIELD;}
									// WITH {;
									if (EGS4.MEDIUM == 0) {
										// $SET-TVSTEP-EM-FIELD;
										// "THIS MACRO ADJUSTS THE PATHLENGTH FOR BENDING IN THE EM FIELD.               "
										// "THIS CORRECTION HAS NOT YET BEEN DISCOVERED. IT IS HIGHER ORDER              "
										// "IN THE PERTURBATION AND MIXES WITH THE MULTIPLE SCATTERING.                  "
										// REPLACE {$SET-TVSTEP-EM-FIELD;} WITH
										// {;}
										de = 0.0;
										// $ADD-WORK-EM-FIELD;
										// "THIS MACRO ADJUSTS DE TO INCORPORATE THE POTENTIAL CHANGE IN THE E FIELD     "
										// REPLACE {$ADD-WORK-EM-FIELD;} WITH {;
										EGS4emf.GET_POTENTIAL("pot1", x0, y0,
												z0, iemf);
										// -----------------
										pot1 = EGS4emf.pot1;
										// -----------------
										x1 = x0 + EGS4.USTEP
												* EGS4.U[EGS4.NP - 1];
										y1 = y0 + EGS4.USTEP
												* EGS4.V[EGS4.NP - 1];
										z1 = z0 + EGS4.USTEP
												* EGS4.W[EGS4.NP - 1];
										EGS4emf.GET_POTENTIAL("pot2", x1, y1,
												z1, iemf);
										// -----------------
										pot2 = EGS4emf.pot2;
										// -----------------
										potdif = pot2 - pot1;
										de = de + potdif;
										// }

										peie = peie - de;
										eie = peie;
										EGS4.E[EGS4.NP - 1] = peie;
									}

								}
								EGS4.E_RANGE = EGS4.VACDST;// e_range = vacdst;
								// $AUSCALL($TRANAUSB);
								iarg = EGS4.$TRANAUSB;
								if (EGS4.iausfl[iarg] != 0) {
									eq.AUSGAB(iarg);
								}

								// "Transport the particle
								EGS4.X[EGS4.NP - 1] = EGS4.X[EGS4.NP - 1]
										+ EGS4.U[EGS4.NP - 1] * EGS4.VSTEP;// vstep;
								EGS4.Y[EGS4.NP - 1] = EGS4.Y[EGS4.NP - 1]
										+ EGS4.V[EGS4.NP - 1] * EGS4.VSTEP;
								EGS4.Z[EGS4.NP - 1] = EGS4.Z[EGS4.NP - 1]
										+ EGS4.W[EGS4.NP - 1] * EGS4.VSTEP;
								EGS4.DNEAR[EGS4.NP - 1] = EGS4.DNEAR[EGS4.NP - 1]
										- EGS4.VSTEP;
								// "(dnear is distance to the nearest boundary
								// " that goes along with particle stack and
								// " which the user's eq.HOWFAR can supply
								// (option)
								EGS4.IROLD = EGS4.IR[EGS4.NP - 1]; // "save
																	// previous
																	// region
								// $SET-ANGLES-EM-FIELD;
								if (EMFb) {
									// "THIS MACRO SETS THE NEW DIRECTION COSINES IN THE PRESENCE OF THE EM FIELD"
									// REPLACE {$SET-ANGLES-EM-FIELD;} WITH {;
									if (EGS4.TVSTEP != 0.0) {
										if (bperp2 != 0.0) {
											ufx = u0 + EGS4.TVSTEP * bperpx;
											ufy = v0 + EGS4.TVSTEP * bperpy;
											ufz = w0 + EGS4.TVSTEP * bperpz;
											bdotub = Bx0 * ufx + By0 * ufy
													+ Bz0 * ufz;
											bsqrd = Bx0 * Bx0 + By0 * By0 + Bz0
													* Bz0;
											bdotub = bdotub / bsqrd;
											ufxpar = Bx0 * bdotub;
											ufypar = By0 * bdotub;
											ufzpar = Bz0 * bdotub;
											ufxprp = ufx - ufxpar;
											ufyprp = ufy - ufypar;
											ufzprp = ufz - ufzpar;
											fnorm = ufxprp * ufxprp + ufyprp
													* ufyprp + ufzprp * ufzprp;
											fnorm = fnorm
													/ (1.0 - (Bx0 * u0 + By0
															* v0 + Bz0 * w0)
															* (Bx0 * u0 + By0
																	* v0 + Bz0
																	* w0)
															/ bsqrd);
											fnorm = Math.sqrt(fnorm);
											EGS4.U[EGS4.NP - 1] = ufxpar
													+ ufxprp / fnorm;
											EGS4.V[EGS4.NP - 1] = ufypar
													+ ufyprp / fnorm;
											EGS4.W[EGS4.NP - 1] = ufzpar
													+ ufzprp / fnorm;
										}
										if (dperp2 != 0.0) {
											test = 0.5 * dperp2 * EGS4.TVSTEP
													* EGS4.TVSTEP;
											fnorm = 1. + test
													* (1. - 0.5 * test);
											EGS4.U[EGS4.NP - 1] = (EGS4.U[EGS4.NP - 1] + EGS4.TVSTEP
													* dperpx)
													/ fnorm;
											EGS4.V[EGS4.NP - 1] = (EGS4.V[EGS4.NP - 1] + EGS4.TVSTEP
													* dperpy)
													/ fnorm;
											EGS4.W[EGS4.NP - 1] = (EGS4.W[EGS4.NP - 1] + EGS4.TVSTEP
													* dperpz)
													/ fnorm;
										}
									}
									// }

								}
								// "default for $SET-ANGLES-EM-FIELD; is ;
								// (null)
								// "(allows for EM field deflection
							} // "end of vacuum step

							// $electron_region_change;=>NOTHING DEFAULT
							// REPLACE {$electron_region_change;} WITH {
							// ir(np) = irnew; irl = irnew; medium = med(irl);
							EGS4.IR[EGS4.NP - 1] = EGS4.IRNEW;
							irl = EGS4.IRNEW;
							EGS4.MEDIUM = EGS4.MED[irl - 1];

							if (EGS4.USTEP != 0.0) {
								// $AUSCALL($TRANAUSA);
								iarg = EGS4.$TRANAUSA;
								if (EGS4.iausfl[iarg] != 0) {
									eq.AUSGAB(iarg);
								}

							}
							if (eie <= EGS4.ECUT[irl - 1])// ecut(irl))
							{
								// go to :ECUT-DISCARD:;
								ECUTb = true;
								break BREMSLOOP;// NEWELECTRON;
							}
							if ((EGS4.USTEP != 0.0) && (EGS4.IDISC < 0)) {
								// go to :USER-ELECTRON-DISCARD:;
								USERb = true;
								break BREMSLOOP;// NEWELECTRON;
							}
							// NEXT :TSTEP: ; //"(Start again at :TSTEP:)
							NEXTb = true;// particle transported see above!!!
							break USTEP;

						} // "Go try another big step in (possibly) new medium
						// @@@@@@@@@@@@@@@@@@@@@@@@THE KEY@@@@
						EGS4.VSTEP = EGS4.USTEP;// vstep = ustep;//vstep =
												// ustep;
						// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
						// ENERGY LOSS:
						if (callHOWFAR) {
							if (EGS4.exact_bca) {
								// "If calleq.HOWFAR=.true. and exact_bca=.true. we are"
								// "in a single scattering mode"
								EGS4.TVSTEP = EGS4.VSTEP;// tvstep = vstep;
								if (EGS4.TVSTEP != EGS4.TUSTEP)// (tvstep ~=
																// tustep)
								{
									// "Boundary was crossed. Shut off single scattering"
									dosingle = false;
								}
							} else {
								// "calleq.HOWFAR=.true. and exact_bca=.false."
								// "=>we are doing an approximate CH step"
								// "calculate the average curved path-length corresponding"
								// "to vstep"
								// $SET-TVSTEP;
								// REPLACE {$SET-TVSTEP;} WITH
								// {
								if (EGS4.VSTEP < ustep0) // ( vstep < ustep0 )
								{
									ekems = EGS4.EKE - 0.5 * EGS4.TUSTEP
											* EGS4.VSTEP / ustep0 * dedx;
									// "This estimates the energy loss to the boundary."
									// "tustep was the intended curved path-length,"
									// "ustep0 is the average transport distance in the initial direction"
									// "       resulting from tustep"
									// "vstep = ustep is the reduced average transport distance in the "
									// "              initial direction due to boundary crossing"
									// $CALCULATE-XI(vstep);
									// REPLACE {$CALCULATE-XI(#);} WITH
									// {
									p2 = ekems * (ekems + EGS4.RMT2);
									beta2 = p2 / (p2 + EGS4.RMSQ);
									chia2 = xccl / (4.0 * blccl * p2);
									// "Note that our chia2 is Moliere chia2/4"
									// "Note also that xcc is now old egs xcc**2"
									// xi = 0.5*xccl/p2/beta2*{P1};
									xi = 0.5 * xccl / p2 / beta2 * EGS4.VSTEP;
									if (EGS4.spin_effects) {
										elkems = Math.log(ekems);
										// $SET INTERVAL elkems,eke;
										Double dbl = new Double(
												EGS4.EKE1[EGS4.MEDIUM - 1]
														* elkems
														+ EGS4.EKE0[EGS4.MEDIUM - 1]);
										lelkems = dbl.intValue();

										if (lelec < 0) {
											// $EVALUATE etap USING
											// etae_ms(elkems);
											etap = EGS4.etae_ms1[lelkems - 1][EGS4.MEDIUM - 1]
													* elkems
													+ EGS4.etae_ms0[lelkems - 1][EGS4.MEDIUM - 1];

											// $EVALUATE xi_corr USING
											// q1ce_ms(elkems);
											xi_corr = EGS4.q1ce_ms1[lelkems - 1][EGS4.MEDIUM - 1]
													* elkems
													+ EGS4.q1ce_ms0[lelkems - 1][EGS4.MEDIUM - 1];

										} else {
											// $EVALUATE etap USING
											// etap_ms(elkems);
											etap = EGS4.etap_ms1[lelkems - 1][EGS4.MEDIUM - 1]
													* elkems
													+ EGS4.etap_ms0[lelkems - 1][EGS4.MEDIUM - 1];

											// $EVALUATE xi_corr USING
											// q1cp_ms(elkems);
											xi_corr = EGS4.q1cp_ms1[lelkems - 1][EGS4.MEDIUM - 1]
													* elkems
													+ EGS4.q1cp_ms0[lelkems - 1][EGS4.MEDIUM - 1];

										}
										chia2 = chia2 * etap;
										xi = xi * xi_corr;
										// $EVALUATE ms_corr USING
										// blcce(elkems);
										ms_corr = EGS4.blcce1[lelkems - 1][EGS4.MEDIUM - 1]
												* elkems
												+ EGS4.blcce0[lelkems - 1][EGS4.MEDIUM - 1];

										blccl = blccl * ms_corr;
									} else {
										xi_corr = 1.0;
										etap = 1.0;
									}
									xi = xi
											* (Math.log(1.0 + 1. / chia2) - 1.0 / (1.0 + chia2));
									// }

									if (xi < 0.1) {
										EGS4.TVSTEP = EGS4.VSTEP
												* (1.0 + xi
														* (0.5 + xi * 0.333333));
									} else {

										if (xi < 0.999999) {
											EGS4.TVSTEP = -EGS4.VSTEP
													* Math.log(1.0 - xi) / xi;
										} else {
											// "This is an error condition because the average transition "
											// "in the initial direction of motion is always smaller than 1/Q1"
											EGS4.STOPPROGRAM = true;
											EGS4.seqStr = " Stoped in SET-TVSTEP because xi > 1!  "
													+ "  \n"
													+ " Medium: "
													+ EGS4.MEDIUM
													+ "  \n"
													+ " Initial energy: "
													+ EGS4.EKE
													+ "  \n"
													+ " Average step energy:  "
													+ ekems
													+ "  \n"
													+ " tustep: "
													+ EGS4.TUSTEP
													+ "  \n"
													+ " ustep0  "
													+ ustep0
													+ "  \n"
													+ " vstep:  "
													+ EGS4.VSTEP
													+ "  \n"
													+ "  ==> xi = "
													+ xi;// +"  \n";
											// if(EGS4.iprint>2)
											eq.printSequence(EGS4.seqStr);
											return;// stop;
										}
									}
								} else {
									EGS4.TVSTEP = EGS4.TUSTEP;// tvstep =
																// tustep;
								}
								// }

							}
							// "Fourth order technique for dedx
							// "Must be done for an approx. CH step or a
							// "single scattering step.
							// $COMPUTE-ELOSS-G(tvstep,eke,elke,lelke,de);
							// REPLACE {$COMPUTE-ELOSS-G(#,#,#,#,#);} WITH
							// {
							tuss = range
									- EGS4.range_ep[qel][lelke - 1][EGS4.MEDIUM - 1]
									/ EGS4.RHOF;
							// (qel,{P4},medium)/rhof;
							// " here tuss is the range between the initial energy and the next lower "
							// " energy on the interpolation grid "
							// if( tuss >= {P1} )
							if (tuss >= EGS4.TVSTEP) {// " Final energy is in the same interpolation bin "
														// $COMPUTE-ELOSS({P1},{P2},{P3},{P4},{P5});
														// HERE(tustep,eke,elke,lelke,de)the
														// same;
														// REPLACE
														// {$COMPUTE-ELOSS(#,#,#,#,#);}
														// WITH
														// {;
								if (lelec < 0) {
									// $EVALUATE dedxmid USING ededx({P3});
									dedxmid = EGS4.EDEDX1[lelke - 1][EGS4.MEDIUM - 1]
											* EGS4.ELKE
											+ EGS4.EDEDX0[lelke - 1][EGS4.MEDIUM - 1];

									// aux = ededx1({P4},medium)/dedxmid;
									aux = EGS4.EDEDX1[lelke - 1][EGS4.MEDIUM - 1]
											/ dedxmid;
								} else {
									// $EVALUATE dedxmid USING pdedx({P3});
									dedxmid = EGS4.PDEDX1[lelke - 1][EGS4.MEDIUM - 1]
											* EGS4.ELKE
											+ EGS4.PDEDX0[lelke - 1][EGS4.MEDIUM - 1];

									// aux = pdedx1({P4},medium)/dedxmid;
									aux = EGS4.PDEDX1[lelke - 1][EGS4.MEDIUM - 1]
											/ dedxmid;
								}
								// {P5} = dedxmid*{P1};
								// //" Energy loss using stopping power at the beginning "
								// de = dedxmid*EGS4.TVSTEP;
								de = dedxmid * EGS4.TVSTEP * EGS4.RHOF;
								fedep = de / EGS4.EKE;// fedep = {P5}/{P2};
								// {P5} =
								// {P5}*(1-0.5*fedep*aux*(1-0.333333*fedep*(aux-1-
								de = de
										* (1.0 - 0.5
												* fedep
												* aux
												* (1.0 - 0.333333
														* fedep
														* (aux - 1.0 - 0.25
																* fedep
																* (2.0 - aux
																		* (4.0 - aux)))));
								// }

								// {P5} = {P5}*rhof; //"IK, rhof bug"
								// de = de*EGS4.RHOF;
							} else {// " Must find first the table index where the step ends using "
									// " pre-calculated ranges                                     "
								lelktmp = lelke;// lelktmp = {P4};
								// tuss = (range - {P1})*rhof;
								tuss = (range - EGS4.TVSTEP) * EGS4.RHOF;
								// " now tuss is the range of the final energy electron "
								// " scaled to the default mass density from PEGS4      "
								if (tuss <= 0.0) {
									// {P5} = {P2}; - TE(medium)*0.99;
									de = EGS4.EKE - EGS4.TE[EGS4.MEDIUM - 1]
											* 0.99;
								}
								// " i.e., if the step we intend to take is longer than the particle "
								// " range, the particle loses its entire energy (which is in {P2})  "
								else {
									while (tuss < EGS4.range_ep[qel][lelktmp - 1][EGS4.MEDIUM - 1])// (qel,lelktmp,medium)
																									// )
									{
										lelktmp = lelktmp - 1;
									}
									elktmp = (lelktmp + 1.0 - EGS4.EKE0[EGS4.MEDIUM - 1])
											/ EGS4.EKE1[EGS4.MEDIUM - 1];
									eketmp = EGS4.E_array[lelktmp][EGS4.MEDIUM - 1];// (lelktmp+1,medium);
									// tuss =
									// EGS4.range_ep[qel][lelktmp][EGS4.MEDIUM-1]
									// - tuss;//(qel,lelktmp+1,medium) - tuss;
									tuss = (EGS4.range_ep[qel][lelktmp][EGS4.MEDIUM - 1] - tuss)
											/ EGS4.RHOF;
									// $COMPUTE-ELOSS(tuss,eketmp,elktmp,lelktmp,{P5});

									// REPLACE {$COMPUTE-ELOSS(#,#,#,#,#);} WITH
									// {;
									if (lelec < 0) {
										// $EVALUATE dedxmid USING ededx({P3});
										dedxmid = EGS4.EDEDX1[lelktmp - 1][EGS4.MEDIUM - 1]
												* elktmp
												+ EGS4.EDEDX0[lelktmp - 1][EGS4.MEDIUM - 1];

										// aux = ededx1({P4},medium)/dedxmid;
										aux = EGS4.EDEDX1[lelktmp - 1][EGS4.MEDIUM - 1]
												/ dedxmid;
									} else {
										// $EVALUATE dedxmid USING pdedx({P3});
										dedxmid = EGS4.PDEDX1[lelktmp - 1][EGS4.MEDIUM - 1]
												* elktmp
												+ EGS4.PDEDX0[lelktmp - 1][EGS4.MEDIUM - 1];

										// aux = pdedx1({P4},medium)/dedxmid;
										aux = EGS4.PDEDX1[lelktmp - 1][EGS4.MEDIUM - 1]
												/ dedxmid;
									}
									// {P5} = dedxmid*{P1};
									// " Energy loss using stopping power at the beginning "
									// de = dedxmid*tuss;
									de = dedxmid * tuss * EGS4.RHOF;
									fedep = de / eketmp;// fedep = {P5}/{P2};
									de = de
											* (1.0 - 0.5
													* fedep
													* aux
													* (1.0 - 0.333333
															* fedep
															* (aux - 1.0 - 0.25
																	* fedep
																	* (2.0 - aux
																			* (4.0 - aux)))));
									// }

									// {P5} = {P5} + {P2} - eketmp;
									de = de + EGS4.EKE - eketmp;
								}
							}
							// }

						} else {
							// "calleq.HOWFAR=.false. => step has not been
							// reduced due to
							// " boundaries
							EGS4.TVSTEP = EGS4.TUSTEP;// tvstep = tustep;
							if (!callmsdist)// ( ~callmsdist )
							{
								// "Second order technique for dedx
								// "Already done in a normal CH step with call
								// to msdist
								// $COMPUTE-ELOSS-G(tvstep,eke,elke,lelke,de);
								// REPLACE {$COMPUTE-ELOSS-G(#,#,#,#,#);} WITH
								// {
								tuss = range
										- EGS4.range_ep[qel][lelke - 1][EGS4.MEDIUM - 1]
										/ EGS4.RHOF;
								// (qel,{P4},medium)/rhof;
								// " here tuss is the range between the initial energy and the next lower "
								// " energy on the interpolation grid "
								// if( tuss >= {P1} )
								if (tuss >= EGS4.TVSTEP) {// " Final energy is in the same interpolation bin "
															// $COMPUTE-ELOSS({P1},{P2},{P3},{P4},{P5});
															// HERE(tustep,eke,elke,lelke,de)the
															// same;
															// REPLACE
															// {$COMPUTE-ELOSS(#,#,#,#,#);}
															// WITH
															// {;
									if (lelec < 0) {
										// $EVALUATE dedxmid USING ededx({P3});
										dedxmid = EGS4.EDEDX1[lelke - 1][EGS4.MEDIUM - 1]
												* EGS4.ELKE
												+ EGS4.EDEDX0[lelke - 1][EGS4.MEDIUM - 1];

										// aux = ededx1({P4},medium)/dedxmid;
										aux = EGS4.EDEDX1[lelke - 1][EGS4.MEDIUM - 1]
												/ dedxmid;
									} else {
										// $EVALUATE dedxmid USING pdedx({P3});
										dedxmid = EGS4.PDEDX1[lelke - 1][EGS4.MEDIUM - 1]
												* EGS4.ELKE
												+ EGS4.PDEDX0[lelke - 1][EGS4.MEDIUM - 1];

										// aux = pdedx1({P4},medium)/dedxmid;
										aux = EGS4.PDEDX1[lelke - 1][EGS4.MEDIUM - 1]
												/ dedxmid;
									}
									// {P5} = dedxmid*{P1};
									// //" Energy loss using stopping power at the beginning "
									// de = dedxmid*EGS4.TVSTEP;
									de = dedxmid * EGS4.TVSTEP * EGS4.RHOF;
									fedep = de / EGS4.EKE;// fedep = {P5}/{P2};
									// {P5} =
									// {P5}*(1-0.5*fedep*aux*(1-0.333333*fedep*(aux-1-
									de = de
											* (1.0 - 0.5
													* fedep
													* aux
													* (1.0 - 0.333333
															* fedep
															* (aux - 1.0 - 0.25
																	* fedep
																	* (2.0 - aux
																			* (4.0 - aux)))));
									// }

									// {P5} = {P5}*rhof; //"IK, rhof bug"
									// de = de*EGS4.RHOF;
								} else {// " Must find first the table index where the step ends using "
										// " pre-calculated ranges                                     "
									lelktmp = lelke;// lelktmp = {P4};
									// tuss = (range - {P1})*rhof;
									tuss = (range - EGS4.TVSTEP) * EGS4.RHOF;
									// " now tuss is the range of the final energy electron "
									// " scaled to the default mass density from PEGS4      "
									if (tuss <= 0.0) {
										// {P5} = {P2};-
										// EGS4.TE[EGS4.MEDIUM-1]*0.99;
										de = EGS4.EKE
												- EGS4.TE[EGS4.MEDIUM - 1]
												* 0.99;
									}
									// " i.e., if the step we intend to take is longer than the particle "
									// " range, the particle loses its entire energy (which is in {P2})  "
									else {
										while (tuss < EGS4.range_ep[qel][lelktmp - 1][EGS4.MEDIUM - 1])// (qel,lelktmp,medium)
																										// )
										{
											lelktmp = lelktmp - 1;
										}
										elktmp = (lelktmp + 1 - EGS4.EKE0[EGS4.MEDIUM - 1])
												/ EGS4.EKE1[EGS4.MEDIUM - 1];
										eketmp = EGS4.E_array[lelktmp][EGS4.MEDIUM - 1];// (lelktmp+1,medium);
										// tuss =
										// EGS4.range_ep[qel][lelktmp][EGS4.MEDIUM-1]
										// - tuss;//(qel,lelktmp+1,medium) -
										// tuss;
										tuss = (EGS4.range_ep[qel][lelktmp][EGS4.MEDIUM - 1] - tuss)
												/ EGS4.RHOF;
										// $COMPUTE-ELOSS(tuss,eketmp,elktmp,lelktmp,{P5});

										// REPLACE {$COMPUTE-ELOSS(#,#,#,#,#);}
										// WITH
										// {;
										if (lelec < 0) {
											// $EVALUATE dedxmid USING
											// ededx({P3});
											dedxmid = EGS4.EDEDX1[lelktmp - 1][EGS4.MEDIUM - 1]
													* elktmp
													+ EGS4.EDEDX0[lelktmp - 1][EGS4.MEDIUM - 1];

											// aux =
											// ededx1({P4},medium)/dedxmid;
											aux = EGS4.EDEDX1[lelktmp - 1][EGS4.MEDIUM - 1]
													/ dedxmid;
										} else {
											// $EVALUATE dedxmid USING
											// pdedx({P3});
											dedxmid = EGS4.PDEDX1[lelktmp - 1][EGS4.MEDIUM - 1]
													* elktmp
													+ EGS4.PDEDX0[lelktmp - 1][EGS4.MEDIUM - 1];

											// aux =
											// pdedx1({P4},medium)/dedxmid;
											aux = EGS4.PDEDX1[lelktmp - 1][EGS4.MEDIUM - 1]
													/ dedxmid;
										}
										// {P5} = dedxmid*{P1};
										// " Energy loss using stopping power at the beginning "
										// de = dedxmid*tuss;
										de = dedxmid * tuss * EGS4.RHOF;
										fedep = de / eketmp;// fedep =
															// {P5}/{P2};
										de = de
												* (1.0 - 0.5
														* fedep
														* aux
														* (1.0 - 0.333333
																* fedep
																* (aux - 1.0 - 0.25
																		* fedep
																		* (2.0 - aux
																				* (4.0 - aux)))));
										// }

										// {P5} = {P5} + {P2} - eketmp;
										de = de + EGS4.EKE - eketmp;
									}
								}
								// }

							}
						}

						// $SET-TVSTEP-EM-FIELD; "additional path length
						// correction in em field
						// "( Calculates tvstep given vstep
						// " default for $SET-TVSTEP-EM-FIELD; is ; (null)
						if (EMFb) {
							// "THIS MACRO ADJUSTS THE PATHLENGTH FOR BENDING IN THE EM FIELD.               "
							// "THIS CORRECTION HAS NOT YET BEEN DISCOVERED. IT IS HIGHER ORDER              "
							// "IN THE PERTURBATION AND MIXES WITH THE MULTIPLE SCATTERING.                  "
							// REPLACE {$SET-TVSTEP-EM-FIELD;} WITH {;}

						}
						// END ENERGY LOSS
						save_de = de; // "the energy loss is used to calculate
										// the number
										// "of MFP gone up to now. If energy
										// loss
										// "fluctuations are implemented, de
										// will be
										// "changed in $DE-FLUCTUATION; => save

						// "The following macro template allows the user to
						// change the
						// "ionization loss.
						// "(Provides a user hook for Landau/Vavilov processes)
						// $DE-FLUCTUATION;
						// REPLACE {$DE-FLUCTUATION;} WITH {;}
						// "default for $DE-FLUCTUATION; is ; (null)
						EGS4.EDEP = de; // "energy deposition variable for user
						// $ADD-WORK-EM-FIELD; "e-loss or gain in em field
						// "Default for $ADD-WORK-EM-FIELD; is ; (null)
						if (EMFb) {
							// "THIS MACRO ADJUSTS DE TO INCORPORATE THE POTENTIAL CHANGE IN THE E FIELD     "
							// REPLACE {$ADD-WORK-EM-FIELD;} WITH {;
							EGS4emf.GET_POTENTIAL("pot1", x0, y0, z0, iemf);
							// -----------------
							pot1 = EGS4emf.pot1;
							// -----------------
							x1 = x0 + EGS4.USTEP * EGS4.U[EGS4.NP - 1];
							y1 = y0 + EGS4.USTEP * EGS4.V[EGS4.NP - 1];
							z1 = z0 + EGS4.USTEP * EGS4.W[EGS4.NP - 1];
							EGS4emf.GET_POTENTIAL("pot2", x1, y1, z1, iemf);
							// -----------------
							pot2 = EGS4emf.pot2;
							// -----------------
							potdif = pot2 - pot1;
							de = de + potdif;
							// }

						}
						ekef = EGS4.EKE - de; // "(final kinetic energy)
						EGS4.EOLD = eie; // "save old value
						EGS4.ENEW = EGS4.EOLD - de; // "energy at end of
													// transport

						// "Now do multiple scattering->COSTHE and SINTHE
						if (!callmsdist) // "everything done if callmsdist =
											// .true.
						{
							if (domultiple) {
								// "Approximated CH step => do multiple
								// scattering
								// "
								// "ekems, elkems, beta2 have been set in either
								// $SET-TUSTEP
								// "or $SET-TVSTEP if spin_effects is .true.,
								// they are
								// "not needed if spin_effects is .false.
								// "
								// "chia2,etap,xi,xi_corr are also set in the
								// above macros
								// "
								// "qel (0 for e-, 1 for e+) and medium are now
								// also required
								// "(for the spin rejection loop)
								// "
								lambda = blccl * EGS4.TVSTEP / beta2 / etap
										/ (1.0 + chia2);
								xi = xi / xi_corr;
								//findindex = true;
								spin_index = true;
								mscatcallindex = 0;// added
								mscat(lambda, chia2, xi, elkems, beta2, qel,
										EGS4.MEDIUM, EGS4.spin_effects);// ,findindex);
								// ,spin_index,EGS4.COSTHE,EGS4.SINTHE);
							} else {
								if (dosingle) {
									// "Single scattering

									ekems = Math.max(ekef, EGS4.ECUT[irl - 1]
											- EGS4.RM);
									p2 = ekems * (ekems + EGS4.RMT2);
									beta2 = p2 / (p2 + EGS4.RMSQ);
									chia2 = EGS4.XCC[EGS4.MEDIUM - 1]
											/ (4.0 * EGS4.BLCC[EGS4.MEDIUM - 1] * p2);
									if (EGS4.spin_effects) {
										elkems = Math.log(ekems);
										// $SET INTERVAL elkems,eke;
										Double dbl = new Double(
												EGS4.EKE1[EGS4.MEDIUM - 1]
														* elkems
														+ EGS4.EKE0[EGS4.MEDIUM - 1]);
										lelkems = dbl.intValue();

										if (lelec < 0) {
											// $EVALUATE etap USING
											// etae_ms(elkems);
											etap = EGS4.etae_ms1[lelkems - 1][EGS4.MEDIUM - 1]
													* elkems
													+ EGS4.etae_ms0[lelkems - 1][EGS4.MEDIUM - 1];
										} else {
											// $EVALUATE etap USING
											// etap_ms(elkems);
											etap = EGS4.etap_ms1[lelkems - 1][EGS4.MEDIUM - 1]
													* elkems
													+ EGS4.etap_ms0[lelkems - 1][EGS4.MEDIUM - 1];
										}
										chia2 = chia2 * etap;
									}
									sscatcallindex = 0;// added
									sscat(chia2, elkems, beta2, qel,
											EGS4.MEDIUM, EGS4.spin_effects);
									// ,EGS4.COSTHE,EGS4.SINTHE);
								} else {
									EGS4.THETA = 0.0; // "No deflection in
														// single scattering
														// model
									EGS4.SINTHE = 0.0;
									EGS4.COSTHE = 1.0;
								}
							}
						}

						// "We now know distance and amount of energy loss for
						// this step,
						// "and the angle by which the electron will be
						// scattered. Hence,
						// "it is time to call the user and inform him of this
						// transport,
						// "after which we will do it.

						// "Now transport, deduct energy loss, and do multiple
						// scatter.
						// e_range = range;
						EGS4.E_RANGE = range;
						// $AUSCALL($TRANAUSB);
						iarg = EGS4.$TRANAUSB;
						if (EGS4.iausfl[iarg] != 0) {
							eq.AUSGAB(iarg);
						}

						// "Transport the particle
						if (callmsdist) {
							// "Deflection and scattering have been
							// calculated/sampled in msdist
							EGS4.U[EGS4.NP - 1] = uscat;
							EGS4.V[EGS4.NP - 1] = vscat;
							EGS4.W[EGS4.NP - 1] = wscat;
							EGS4.X[EGS4.NP - 1] = xtrans;
							EGS4.Y[EGS4.NP - 1] = ytrans;
							EGS4.Z[EGS4.NP - 1] = ztrans;
						} else {
							EGS4.X[EGS4.NP - 1] = EGS4.X[EGS4.NP - 1]
									+ EGS4.U[EGS4.NP - 1] * EGS4.VSTEP;
							EGS4.Y[EGS4.NP - 1] = EGS4.Y[EGS4.NP - 1]
									+ EGS4.V[EGS4.NP - 1] * EGS4.VSTEP;
							EGS4.Z[EGS4.NP - 1] = EGS4.Z[EGS4.NP - 1]
									+ EGS4.W[EGS4.NP - 1] * EGS4.VSTEP;
							if (domultiple || dosingle) {
								// //$SELECT-AZIMUTHAL-ANGLE and
								// OLD-PARTICLE:->COSTHE and SINTHE
								UPHI(2, 1); // "Apply the deflection, save call
											// to uphi if
											// "no deflection in a single
											// scattering mode
							}
						}

						EGS4.DNEAR[EGS4.NP - 1] = EGS4.DNEAR[EGS4.NP - 1]
								- EGS4.VSTEP;
						EGS4.IROLD = EGS4.IR[EGS4.NP - 1]; // "save previous
															// region
						// $SET-ANGLES-EM-FIELD;
						// "Default for $SET-ANGLES-EM-FIELD; is ; (null)
						if (EMFb) {
							// "THIS MACRO SETS THE NEW DIRECTION COSINES IN THE PRESENCE OF THE EM FIELD"
							// REPLACE {$SET-ANGLES-EM-FIELD;} WITH {;
							if (EGS4.TVSTEP != 0.0) {
								if (bperp2 != 0.0) {
									ufx = u0 + EGS4.TVSTEP * bperpx;
									ufy = v0 + EGS4.TVSTEP * bperpy;
									ufz = w0 + EGS4.TVSTEP * bperpz;
									bdotub = Bx0 * ufx + By0 * ufy + Bz0 * ufz;
									bsqrd = Bx0 * Bx0 + By0 * By0 + Bz0 * Bz0;
									bdotub = bdotub / bsqrd;
									ufxpar = Bx0 * bdotub;
									ufypar = By0 * bdotub;
									ufzpar = Bz0 * bdotub;
									ufxprp = ufx - ufxpar;
									ufyprp = ufy - ufypar;
									ufzprp = ufz - ufzpar;
									fnorm = ufxprp * ufxprp + ufyprp * ufyprp
											+ ufzprp * ufzprp;
									fnorm = fnorm
											/ (1.0 - (Bx0 * u0 + By0 * v0 + Bz0
													* w0)
													* (Bx0 * u0 + By0 * v0 + Bz0
															* w0) / bsqrd);
									fnorm = Math.sqrt(fnorm);
									EGS4.U[EGS4.NP - 1] = ufxpar + ufxprp
											/ fnorm;
									EGS4.V[EGS4.NP - 1] = ufypar + ufyprp
											/ fnorm;
									EGS4.W[EGS4.NP - 1] = ufzpar + ufzprp
											/ fnorm;
								}
								if (dperp2 != 0.0) {
									test = 0.5 * dperp2 * EGS4.TVSTEP
											* EGS4.TVSTEP;
									fnorm = 1. + test * (1. - 0.5 * test);
									EGS4.U[EGS4.NP - 1] = (EGS4.U[EGS4.NP - 1] + EGS4.TVSTEP
											* dperpx)
											/ fnorm;
									EGS4.V[EGS4.NP - 1] = (EGS4.V[EGS4.NP - 1] + EGS4.TVSTEP
											* dperpy)
											/ fnorm;
									EGS4.W[EGS4.NP - 1] = (EGS4.W[EGS4.NP - 1] + EGS4.TVSTEP
											* dperpz)
											/ fnorm;
								}
							}
							// }
						}

						// "Now done with multiple scattering,
						// "update energy and see if below cut
						peie = peie - EGS4.EDEP;
						eie = peie;
						EGS4.E[EGS4.NP - 1] = peie;

						// "IF( irnew ~= irl & eie <= ecut(irl)) [
						// "IK: the above is clearly a bug. If the particle energy falls "
						// "    below ecut, but the particle is actually entering a new "
						// "    region, the discard will happen in the current region "
						// "    instead the next. If the particle is a positron, all "
						// "    resulting annihilation photons will have the new position "
						// "    but the old region => confusion in the geometry routine "
						// "    is very likely.      Jan 27 2004 "
						if ((EGS4.IRNEW == irl) && (eie <= EGS4.ECUT[irl - 1])) {
							// go to :ECUT-DISCARD:;
							ECUTb = true;
							break BREMSLOOP;// NEWELECTRON;

						}

						EGS4.MEDOLD = EGS4.MEDIUM;// medold = medium;
						if (EGS4.MEDIUM != 0) {
							//ekeold = EGS4.EKE;
							EGS4.EKE = eie - EGS4.RM; // "update kinetic energy
							EGS4.ELKE = Math.log(EGS4.EKE);
							// $SET INTERVAL elke,eke; //"Get updated interval
							Double dbl = new Double(EGS4.EKE1[EGS4.MEDIUM - 1]
									* EGS4.ELKE + EGS4.EKE0[EGS4.MEDIUM - 1]);
							lelke = dbl.intValue();

						}

						if (EGS4.IRNEW != EGS4.IROLD) {
							// $electron_region_change;
							// REPLACE {$electron_region_change;} WITH {
							// ir(np) = irnew; irl = irnew; medium = med(irl);
							EGS4.IR[EGS4.NP - 1] = EGS4.IRNEW;
							irl = EGS4.IRNEW;
							EGS4.MEDIUM = EGS4.MED[irl - 1];

						}

						// "After transport call to user scoring routine
						// $AUSCALL($TRANAUSA);
						iarg = EGS4.$TRANAUSA;
						if (EGS4.iausfl[iarg] != 0) {
							eq.AUSGAB(iarg);
						}

						if (eie <= EGS4.ECUT[irl - 1]) {
							// go to :ECUT-DISCARD:;
							ECUTb = true;
							break BREMSLOOP;// NEWELECTRON;

						}

						// "Now check for deferred discard request. May have
						// been set
						// "by either eq.HOWFAR, or one of the transport
						// eq.AUSGAB calls
						if (EGS4.IDISC < 0) {
							// go to :USER-ELECTRON-DISCARD:;
							USERb = true;
							break BREMSLOOP;// NEWELECTRON;
						}

						if (EGS4.MEDIUM != EGS4.MEDOLD)// medold)
						{
							// NEXT :TSTEP:;
							NEXTb = true;// ok, particle was transported!!
							break USTEP;

						}

						// $USER_CONTROLS_TSTEP_RECURSION;
						// REPLACE {$USER_CONTROLS_TSTEP_RECURSION;} WITH {;}
						// "NRCC update 87/12/08--default is null
						EGS4Macro.USER_CONTROLS_TSTEP_RECURSION();
						if (EGS4.nextB) {
							NEXTb = true;// ok, particle was transported!!
							break USTEP;
						}

						// $UPDATE-DEMFP;
						// " The following macro updates demfp. As energy loss is used as the  "
						// " 'path-length' variable (see above), it just substracts the energy "
						// " loss for the step.                                                "
						// REPLACE {$UPDATE-DEMFP;} WITH
						// {
						demfp = demfp - save_de * sig;
						total_de = total_de - save_de;
						total_tstep = total_tstep - EGS4.TVSTEP * EGS4.RHOF;// tvstep*rhof;
						if (total_tstep < 1.e-9) {
							demfp = 0.0;
						}
						// }

						// ---------------------
						if (demfp < $EPSEMFP) {
							break USTEP;
						}
						// --------------------

					} // UNTIL(demfp < $EPSEMFP); //"end ustep loop

					// "Compute final sigma to see if resample is needed.
					// "this will take the energy variation of the sigma into
					// "account using the fictitious sigma method.
					if (!NEXTb) {
						// $EVALUATE-SIGF;
						// REPLACE {$EVALUATE-SIGF;} WITH
						// "        ==============="
						// {;
						if (lelec < 0) {
							// $EVALUATE sigf USING esig(elke);
							sigf = EGS4.ESIG1[lelke - 1][EGS4.MEDIUM - 1]
									* EGS4.ELKE
									+ EGS4.ESIG0[lelke - 1][EGS4.MEDIUM - 1];
							// $EVALUATE dedx0 USING ededx(elke);
							dedx0 = EGS4.EDEDX1[lelke - 1][EGS4.MEDIUM - 1]
									* EGS4.ELKE
									+ EGS4.EDEDX0[lelke - 1][EGS4.MEDIUM - 1];
							sigf = sigf / dedx0;
						} else {
							// $EVALUATE sigf USING psig(elke);
							sigf = EGS4.PSIG1[lelke - 1][EGS4.MEDIUM - 1]
									* EGS4.ELKE
									+ EGS4.PSIG0[lelke - 1][EGS4.MEDIUM - 1];
							// $EVALUATE dedx0 USING pdedx(elke);
							dedx0 = EGS4.PDEDX1[lelke - 1][EGS4.MEDIUM - 1]
									* EGS4.ELKE
									+ EGS4.PDEDX0[lelke - 1][EGS4.MEDIUM - 1];
							sigf = sigf / dedx0;
						}
						// }

						sigratio = sigf / sig0;

						rfict = EGS4.random01();
						// ------------------------------
						if (rfict <= sigratio)
							break TSTEP;
						// -----------------------------------
					}// if (!NEXTb)
				}// UNTIL (rfict <= sigratio) ; "end tstep loop

				// " Now sample electron interaction

				if (lelec < 0) {
					// "e-,check branching ratio
					// $EVALUATE ebr1 USING ebr1(elke); "ebr1 = brems/total
					// {P1}={P2}1(L{P3},MEDIUM)*{P3}+{P2}0(L{P3},MEDIUM)
					ebr1 = EGS4.EBR11[lelke - 1][EGS4.MEDIUM - 1] * EGS4.ELKE
							+ EGS4.EBR10[lelke - 1][EGS4.MEDIUM - 1];

					rnno24 = EGS4.random01();
					if (rnno24 <= ebr1) {
						// "It was bremsstrahlung
						// go to :EBREMS:;
						EBREMSb = true;
						break NEWELECTRON;
					} else {
						// "It was Moller, but first check the kinematics.
						// "However, if EII is on, we should still permit an
						// interaction
						// "even if E<moller threshold as EII interactions go
						// down to
						// "the ionization threshold which may be less than
						// thmoll.
						if ((EGS4.E[EGS4.NP - 1] <= EGS4.THMOLL[EGS4.MEDIUM - 1])
								&& (EGS4.eii_flag == 0))
						// "(thmoll = lower Moller threshold)
						{
							// "Not enough energy for Moller, so
							// "force it to be a bremsstrahlung---provided ok
							// kinematically.
							if (ebr1 <= 0.0) {
								break NEWELECTRON;// it also repeat if
													// newelectron is set to
													// false
								// go to :NEWELECTRON:;
							}
							// "Brems not allowed either.
							// go to :EBREMS:;
							EBREMSb = true;
							break NEWELECTRON;

						}
						// $AUSCALL($MOLLAUSB);
						iarg = EGS4.$MOLLAUSB;
						if (EGS4.iausfl[iarg] != 0) {
							eq.AUSGAB(iarg);
						}

						MOLLER();
						// "The following macro template allows the user to
						// change the
						// "particle selection scheme (e.g., adding importance
						// sampling
						// "such as splitting, leading particle selection,
						// etc.).
						// "(Default macro is template
						// '$PARTICLE-SELECTION-ELECTR'
						// "which in turn has the 'null' replacement ';')
						// $PARTICLE-SELECTION-MOLLER;
						// $AUSCALL($MOLLAUSA);
						iarg = EGS4.$MOLLAUSA;
						if (EGS4.iausfl[iarg] != 0) {
							eq.AUSGAB(iarg);
						}
						if (EGS4.IQ[EGS4.NP - 1] == 0)
							return;
					}

					// go to :NEWELECTRON:; //"Electron is lowest energy-follow
					// it
					break NEWELECTRON;
				}

				// "e+ interaction. pbr1 = brems/(brems + bhabha + annih)
				// $EVALUATE pbr1 USING pbr1(elke);
				pbr1 = EGS4.PBR11[lelke - 1][EGS4.MEDIUM - 1] * EGS4.ELKE
						+ EGS4.PBR10[lelke - 1][EGS4.MEDIUM - 1];

				rnno25 = EGS4.random01();
				if (rnno25 < pbr1) {
					// go to :EBREMS:;
					EBREMSb = true;
					break NEWELECTRON;

				}// "It was bremsstrahlung
					// "Decide between bhabha and annihilation
					// "pbr2 is (brems + bhabha)/(brems + bhabha + annih)
					// $EVALUATE pbr2 USING pbr2(elke);
				pbr2 = EGS4.PBR21[lelke - 1][EGS4.MEDIUM - 1] * EGS4.ELKE
						+ EGS4.PBR20[lelke - 1][EGS4.MEDIUM - 1];

				if (rnno25 < pbr2) {
					// "It is bhabha
					// $AUSCALL($BHABAUSB);
					iarg = EGS4.$BHABAUSB;
					if (EGS4.iausfl[iarg] != 0) {
						eq.AUSGAB(iarg);
					}

					BHABHA();
					// "The following macro template allows the user to change
					// the
					// "particle selection scheme (e.g., adding importance
					// sampling
					// "such as splitting, leading particle selection, etc.).
					// (default
					// "macro is template '$PARTICLE-SELECTION-ELECTR' which in
					// turn
					// "has the 'null' replacement ';')
					// $PARTICLE-SELECTION-BHABHA;
					// $AUSCALL($BHABAUSA);
					iarg = EGS4.$BHABAUSA;
					if (EGS4.iausfl[iarg] != 0) {
						eq.AUSGAB(iarg);
					}
					if (EGS4.IQ[EGS4.NP - 1] == 0)
						return;
				} else {
					// "It is in-flight annihilation
					// $AUSCALL($ANNIHFAUSB);
					iarg = EGS4.$ANNIHFAUSB;
					if (EGS4.iausfl[iarg] != 0) {
						eq.AUSGAB(iarg);
					}

					ANNIH();
					// "The following macro template allows the user to change
					// the
					// "particle selection scheme (e.g., adding importance
					// sampling
					// "such as splitting, leading particle selection, etc.).
					// (default
					// "macro is template '$PARTICLE-SELECTION-ELECTR' which in
					// turn
					// "has the 'null' replacement ';')
					// $PARTICLE-SELECTION-ANNIH;
					// $AUSCALL($ANNIHFAUSA);
					iarg = EGS4.$ANNIHFAUSA;
					if (EGS4.iausfl[iarg] != 0) {
						eq.AUSGAB(iarg);
					}

					// EXIT :NEWELECTRON:; "i.e., in order to return to shower
					NEWELECTRONb = true;
					break NEWELECTRON;
					// "After annihilation the gammas are bound to be the lowest
					// energy
					// "particles, so return and follow them.
				} // "end pbr2 else

			} // REPEAT "NEWELECTRON=======>//main loop

			// return; //"i.e., return to shower

			// "---------------------------------------------
			// "Bremsstrahlung-call section
			// "---------------------------------------------
			// :EBREMS:
			// $AUSCALL($BREMAUSB); //EBREMSb=true;
			if (EBREMSb) {
				iarg = EGS4.$BREMAUSB;
				if (EGS4.iausfl[iarg] != 0) {
					eq.AUSGAB(iarg);
				}

				BREMS();
				// "The following macro template allows the user to change the
				// particle
				// "selection scheme (e.g., adding importance sampling such as
				// splitting,
				// "leading particle selection, etc.). (default macro is
				// template
				// "'$PARTICLE-SELECTION-ELECTR' which in turn has the 'null'
				// replacement ';')
				// $PARTICLE-SELECTION-BREMS;
				// REPLACE {$PARTICLE-SELECTION-ELECTR;} WITH {;}
				// REPLACE {$PARTICLE-SELECTION-ANNIH;} WITH {
				// $PARTICLE-SELECTION-ELECTR;}
				// REPLACE {$PARTICLE-SELECTION-ANNIHREST;} WITH {
				// $PARTICLE-SELECTION-ELECTR;}
				// REPLACE {$PARTICLE-SELECTION-BHABHA;} WITH {
				// $PARTICLE-SELECTION-ELECTR;}
				// REPLACE {$PARTICLE-SELECTION-BREMS;} WITH {
				// $PARTICLE-SELECTION-ELECTR;}
				// REPLACE {$PARTICLE-SELECTION-MOLLER;}
				// WITH {$PARTICLE-SELECTION-ELECTR;}

				// $AUSCALL($BREMAUSA);
				iarg = EGS4.$BREMAUSA;
				if (EGS4.iausfl[iarg] != 0) {
					eq.AUSGAB(iarg);
				}

				if (EGS4.IQ[EGS4.NP - 1] == 0) {
					// "Photon was selected.
					return;
					// "i.e., return to shower
				}// //otherwise loop again
					// else
					// {
					// "Electron was selected
					// go to :NEWELECTRON:;
				// }
			}// if(EBREMSb)
			else {
				if (NEWELECTRONb)// otherwise loop again
					return;
			}

		}// go to NEWELECTRON loop for outside main
			// loop;#############################

		// "---------------------------------------------
		// "Electron cutoff energy discard section
		// "---------------------------------------------
		// :ECUT-DISCARD:
		if (ECUTb) {
			if (EGS4.MEDIUM > 0) {
				if (eie > EGS4.AE[EGS4.MEDIUM - 1]) {
					idr = EGS4.$EGSCUTAUS;
					if (lelec < 0) {
						EGS4.EDEP = EGS4.E[EGS4.NP - 1] - EGS4.PRM;
					} else {
						// $POSITRON-ECUT-DISCARD;
						// REPLACE {$POSITRON-ECUT-DISCARD;} WITH
						// {EDEP=PEIE-PRM;}
						EGS4.EDEP = peie - EGS4.PRM;
					}
				} else {
					idr = EGS4.$PEGSCUTAUS;
					EGS4.EDEP = EGS4.E[EGS4.NP - 1] - EGS4.PRM;
				}
			} else {
				idr = EGS4.$EGSCUTAUS;
				EGS4.EDEP = EGS4.E[EGS4.NP - 1] - EGS4.PRM;// edep = e(np) -
															// prm;
			}

			// $ELECTRON-TRACK-END;
			// "The default replacement for this macros is "
			// "          $AUSCALL(idr);                   "
			// "Use this macro if you wish to modify the   "
			// "treatment of track ends                    "
			// REPLACE {$ELECTRON-TRACK-END;} WITH {; $AUSCALL(idr); }
			iarg = idr;
			if (EGS4.iausfl[iarg] != 0) {
				eq.AUSGAB(iarg);
			}

			// :POSITRON-ANNIHILATION:; "NRCC extension 86/9/12
			if (lelec > 0) {
				// "It's a positron. Produce annihilation gammas if edep < peie
				if (EGS4.EDEP < peie) {
					// $AUSCALL($ANNIHRAUSB);
					iarg = EGS4.$ANNIHRAUSB;
					if (EGS4.iausfl[iarg] != 0) {
						eq.AUSGAB(iarg);
					}

					ANNIH_AT_REST();
					// "The following macro template allows the user to change
					// the
					// "particle selection scheme (e.g., adding importance
					// sampling
					// "such as splitting, leading particle selection, etc.).
					// (Default
					// "macro is template '$PARTICLE-SELECTION-ELECTR' which in
					// turn
					// "has the 'null' replacement ';')

					// $PARTICLE-SELECTION-ANNIHREST;see above
					// $AUSCALL($ANNIHRAUSA);
					iarg = EGS4.$ANNIHRAUSA;
					if (EGS4.iausfl[iarg] != 0) {
						eq.AUSGAB(iarg);
					}

					// "Now discard the positron and take normal return to
					// follow
					// "the annihilation gammas.
					return; // "i.e., return to shower
				}
			} // "end of positron block

			EGS4.NP = EGS4.NP - 1;
			IRCODE = 2;// ircode = 2; "tell shower an e- or un-annihilated
			// "e+ has been discarded

			return; // "i.e., return to shower"
		}// if(ECUTb)
			// "---------------------------------------------
			// "User requested electron discard section
			// "---------------------------------------------
			// :USER-ELECTRON-DISCARD:
		if (USERb) {
			EGS4.IDISC = Math.abs(EGS4.IDISC);

			if ((lelec < 0) || (EGS4.IDISC == 99)) {
				EGS4.EDEP = EGS4.E[EGS4.NP - 1] - EGS4.PRM;
			} else {
				EGS4.EDEP = EGS4.E[EGS4.NP - 1] + EGS4.PRM;
			}

			// $AUSCALL($USERDAUS);
			iarg = EGS4.$USERDAUS;
			if (EGS4.iausfl[iarg] != 0) {
				eq.AUSGAB(iarg);
			}

			if (EGS4.IDISC == 99) {
				// goto :POSITRON-ANNIHILATION:;
				// :POSITRON-ANNIHILATION:; "NRCC extension 86/9/12

				if (lelec > 0) {
					// "It's a positron. Produce annihilation gammas if edep <
					// peie
					if (EGS4.EDEP < peie) {
						// $AUSCALL($ANNIHRAUSB);
						iarg = EGS4.$ANNIHRAUSB;
						if (EGS4.iausfl[iarg] != 0) {
							eq.AUSGAB(iarg);
						}

						ANNIH_AT_REST();
						// "The following macro template allows the user to
						// change the
						// "particle selection scheme (e.g., adding importance
						// sampling
						// "such as splitting, leading particle selection,
						// etc.). (Default
						// "macro is template '$PARTICLE-SELECTION-ELECTR' which
						// in turn
						// "has the 'null' replacement ';')

						// $PARTICLE-SELECTION-ANNIHREST;see above
						// $AUSCALL($ANNIHRAUSA);
						iarg = EGS4.$ANNIHRAUSA;
						if (EGS4.iausfl[iarg] != 0) {
							eq.AUSGAB(iarg);
						}

						// "Now discard the positron and take normal return to
						// follow
						// "the annihilation gammas.
						return; // "i.e., return to shower
					}
				} // "end of positron block

			}

			// np = np - 1; ircode = 2;
			EGS4.NP = EGS4.NP - 1;
			IRCODE = 2;

			return; // "i.e., return to shower
		}// if(USERb)
			// end; "End of subroutine electr
			// "*******************************************************************************

	}

	// "******************************************************************"
	// "                               National Research Council of Canada"
	/**
	 * Called by ELECTR. Samples bremsstrahlung energy using Coulomb corrected Bethe-Heitler above 50 MeV and 
	 * Bethe-Heitler below 50 MeV.<br>
	 * if ibr_nist = 0, the NIST bremsstrahlung cross section data base (prepared in a form of an alias table for rapid sampling) is used. <br>
	 * To use bremsstrahlung splitting, set nbr_split to the desired number greater than 1 (1 is the default). Be aware that event-by-event energy conservation 
	 * is NOT guaranteed, so don't use for calculations where this is important (e.g. calculation of detector response functions). The result will be nbr_split photons, 
	 * all with the weight wt(npold)/nbr_split, and an electron with the original weight and energy given by the incident energy - energy of last photon.
	 * 
	 */
	protected static void BREMS() {
		// "                                                                  "
		// "******************************************************************"
		// "   Samples bremsstrahlung energy using                            "
		// "    - Coulomb corrected Bethe-Heitler above 50 MeV                "
		// "    - Bethe-Heitler below 50 MeV                                  "
		// "   if ibr_nist = 0, or                                            "
		// "    - the NIST bremsstrahlung cross section data base             "
		// "      (prepared in a form of an alias table for rapid sampling)   "
		// "   if ibr_nist = 1                                                "
		// "   and direction using                                            "
		// "    - formula 2BS from from Koch and Motz if IBRDST=1             "
		// "    - leading term of the brems angular dsstr. if IBRDST=0        "
		// "    - photon direction = electron direction if IBRDST<0           "
		// "                                                                  "
		// "   This version replaces the original EGS4 implementation         "
		// "   because of a bug discovered in the EGS4 brems routine          "
		// "   In order to work properly, the parameter DL1,..,DL6            "
		// "   are re-calculated in subroutine fix_brems which is called      "
		// "   from HATCH                                                     "
		// "   In addition, this version has the internal capability of       "
		// "   bremsstrahlung splitting.                                      "
		// "   To use bremsstrahlung splitting, set nbr_split (COMON/BREMPR/) "
		// "   to the desired number > 1 (1 is the default)                   "
		// "   Be aware that event-by-event energy conservation is NOT        "
		// "   guaranteed, so don't use for calculations where this is        "
		// "   important (e.g. calculation of detector response functions)    "
		// "   The result will be nbr_split photons, all with the weight      "
		// "   wt(npold)/nbr_split, and an electron with the original weight  "
		// "   and energy given by the incident energy - energy of last photon"
		// "                                                                  "
		// " I. Kawrakow, January 2000                                        "
		// "                                                                  "
		// "******************************************************************"

		// ; Copyright NRC;

		// $COMIN-BREMS; "DEFAULT REPLACEMENT PRODUCES THE FOLLOWING:
		// "COMIN/DEBUG,BREMPR,EGS-VARIANCE-REDUCTION,   "
		// "STACK,THRESH,UPHIOT,USEFUL,RANDOM/;"

		// $DEFINE-LOCAL-VARIABLES-BREMS;
		// REPLACE {$DEFINE-LOCAL-VARIABLES-BREMS;} WITH
		// {;
		// "Local variables in order of appearance"
		double PEIE = 0.0;// "precise incident electron energy"
		double PESG = 0.0;// "presice energy of emitted photon"
		double PESE = 0.0;// "precise total energy of scattered electron"
		double EIE = 0.0;// "total incident electron energy"
		double EKIN = 0.0;// "kinetic incident energy"
		double brmin = 0.0;// " ap(medium)/ekin"
		double waux = 0.0;// "for faster sampling of 1/br"
		double aux = 0.0;// "ese/eie"
		double r1 = 0.0;// "a random number"
		double ajj = 0.0;// "for energy bin determination if alias sampling is employed"
		//double alias_sample1 = 0.0;//
		double RNNO06 = 0.0;// "random number"
		double RNNO07 = 0.0;// "random number"
		double BR = 0.0;// "energy fraction of secondary photon"
		double ESG = 0.0;// "energy of secondary photon"
		double ESE = 0.0;// "total energy of secondary electron"
		double DELTA = 0.0;// "scaled momentum transfer"
		double phi1 = 0.0;// "screening function"
		double phi2 = 0.0;// "screening function"
		double REJF = 0.0;// "screening rejection function"

		// "Brems angle selection variables"
		double a = 0.0;//
		double b = 0.0;//
		double c = 0.0;// "direction cosines of incident `electron'"
		double sinpsi = 0.0;//
		double sindel = 0.0;//
		double cosdel = 0.0;//
		double us = 0.0;//
		double vs = 0.0;//
		// "all used for rotations"
		double ztarg = 0.0;// "(Zeff**1/3/111)**2, used for 2BS angle sampling"
		double tteie = 0.0;// "total energy in units of rest energy"
		double beta = 0.0;// "electron velocity in units of speed of light"
		double y2max = 0.0;// "maximum possible scaled angle"
		double y2maxi = 0.0;// "inverse of the above"
		double ttese = 0.0;// "new electron energy in units of rm"
		double rjarg1 = 0.0;//
		double rjarg2 = 0.0;//
		double rjarg3 = 0.0;//
		//double rejmin = 0.0;//
		//double rejmid = 0.0;//
		double rejmax = 0.0;//
		//double rejtop = 0.0;//
		double rejtst = 0.0;//
		// "all of them used for angle rejection function calcs"
		double esedei = 0.0;// "new total energy over old total energy"
		double y2tst = 0.0;// "scaled angle, costhe = 1 - 2*y2tst/y2max"
		double y2tst1 = 0.0;//
		double rtest = 0.0;// "random number for rejection"
		double xphi = 0.0;//
		double yphi = 0.0;//
		double xphi2 = 0.0;//
		double yphi2 = 0.0;//
		double rhophi2 = 0.0;//
		double cphi = 0.0;//
		double sphi = 0.0;//
		// "all of the above is for azimuthal angle sampling"

		int L = 0;
		int L1 = 0;
		int ibr = 0;
		int jj = 0;
		int j = 0;

		// }

		double z2max = 0.0;
		double z2maxi = 0.0;
		double aux1 = 0.0;
		double aux3 = 0.0;
		double aux4 = 0.0;
		double aux5 = 0.0;
		double aux2 = 0.0;

		if (EGS4.nbr_split < 1)
			return; // "i.e. the user can turn off brems production"
		// "by setting nbr_split to zero!"

		EGS4.NPold = EGS4.NP;// "Set the old stack counter"
		PEIE = EGS4.E[EGS4.NP - 1]; // "PRECISE ENERGY OF INCIDENT 'ELECTRON'"
		EIE = PEIE; // "ENERGY OF INCIDENT 'ELECTRON'"
		// @=======================================
		// weight = wt(np)/nbr_split;
		double weight = EGS4.WT[EGS4.NP - 1] / EGS4.nbr_split;
		// ============================
		// "   DECIDE WHICH DISTRIBUTION TO USE (B-H COULOMB CORRECTED IS     "
		// "   USED FROM 50 TO 20000 MEV, B-H IS USED 1.5 TO 50 MEV)          "
		if (EIE < 50.0) {
			L = 1;
		} else {
			L = 3;
		}
		L1 = L + 1;

		EKIN = PEIE - EGS4.PRM;
		brmin = EGS4.AP[EGS4.MEDIUM - 1] / EKIN;
		// "waux = -log(brmin);"
		waux = EGS4.ELKE - EGS4.log_ap[EGS4.MEDIUM - 1]; // "this saves the time consuming log evaluation"
		// "log_ap = log(ap(medium)) is calculated in   "
		// "fix_brems for each medium, elke is needed   "
		// "in electr to calculate the branching ratios "
		// "and therefore it must be known at this point"

		if (EGS4.ibrdst >= 0) {// "inrdst >=0 means we will sample the photon emmision"
								// "angle from KM-2BS (ibrdst=1) or from the leading"
								// "term (ibrdst=0). If nbr_split > 1, we can re-use"
								// "the following quantities several time"

			a = EGS4.U[EGS4.NP - 1];
			b = EGS4.V[EGS4.NP - 1];
			c = EGS4.W[EGS4.NP - 1];
			sinpsi = a * a + b * b;
			if (sinpsi > 1.e-20) {
				sinpsi = Math.sqrt(sinpsi);
				sindel = b / sinpsi;
				cosdel = a / sinpsi;
			}

			ztarg = EGS4.ZBRANG[EGS4.MEDIUM - 1];// (medium);
			tteie = EIE / EGS4.RM;
			beta = Math.sqrt((tteie - 1.0) * (tteie + 1.0)) / tteie;
			y2max = 2.0 * beta * (1.0 + beta) * tteie * tteie;
			y2maxi = 1.0 / y2max;
			if (EGS4.ibrdst == 1) {
				z2max = y2max + 1.0;
				z2maxi = Math.sqrt(z2max);
			}

		}

		if (EGS4.ibr_nist == 1) {
			ajj = 1.0
					+ (waux + EGS4.log_ap[EGS4.MEDIUM - 1] - EGS4.nb_lemin[EGS4.MEDIUM - 1])
					* EGS4.nb_dlei[EGS4.MEDIUM - 1];
			Double dbl = new Double(ajj);
			jj = dbl.intValue();
			ajj = ajj - jj;
			if (jj > EGS4.$MXBRES) {
				jj = EGS4.$MXBRES;
				ajj = -1.0;
			}
		}

		// DO ibr = 1,nbr_split [
		for (ibr = 1; ibr <= EGS4.nbr_split; ibr++) {
			if (EGS4.ibr_nist == 1) {// "use the NIST bremsstrahlung cross section"
										// "data base"
				if (EKIN > EGS4.nb_emin[EGS4.MEDIUM - 1]) {
					r1 = EGS4.random01();
					if (r1 < ajj) {
						j = jj + 1;
					} else {
						j = jj;
					}
					// ----------------------
					// nb_xdata=new double[$MXBRXS+1][$MXBRES][$MXMED];
					double[] patx = new double[EGS4.$MXBRXS + 1];
					double[] patf = new double[EGS4.$MXBRXS + 1];
					double[] patw = new double[EGS4.$MXBRXS];
					int[] pati = new int[EGS4.$MXBRXS];
					for (int ll = 0; ll <= EGS4.$MXBRXS; ll++) {
						patx[ll] = EGS4.nb_xdata[ll][j - 1][EGS4.MEDIUM - 1];
						patf[ll] = EGS4.nb_fdata[ll][j - 1][EGS4.MEDIUM - 1];
					}
					for (int ll = 0; ll <= EGS4.$MXBRXS - 1; ll++) {
						patw[ll] = EGS4.nb_wdata[ll][j - 1][EGS4.MEDIUM - 1];
						pati[ll] = EGS4.nb_idata[ll][j - 1][EGS4.MEDIUM - 1];
					}

					// ------------------
					// br =
					// EGS4.alias_sample1(EGS4.$MXBRXS,nb_xdata(0,j,medium),
					// nb_fdata(0,j,medium),
					// nb_wdata(1,j,medium),nb_idata(1,j,medium));
					BR = EGS4.alias_sample1(EGS4.$MXBRXS, patx, patf, patw,
							pati);
					// System.out.println(BR);
				} else {
					BR = EGS4.random01();
				}
				ESG = EGS4.AP[EGS4.MEDIUM - 1] * Math.exp(BR * waux);
				PESG = ESG;
				PESE = PEIE - PESG;
				ESE = PESE;
			} else {

				// LOOP
				while (true) {// "User wants to use Bethe-Heitler"

					RNNO06 = EGS4.random01();
					RNNO07 = EGS4.random01();
					BR = brmin * Math.exp(RNNO06 * waux);
					ESG = EKIN * BR;
					PESG = ESG;
					PESE = PEIE - PESG;
					ESE = PESE;// pese = peie - pesg; ese = pese;
					DELTA = ESG / EIE / ESE * EGS4.DELCM[EGS4.MEDIUM - 1];
					aux = ESE / EIE;
					if (DELTA < 1.0) {
						phi1 = EGS4.DL1[L - 1][EGS4.MEDIUM - 1]
								+ DELTA
								* (EGS4.DL2[L - 1][EGS4.MEDIUM - 1] + DELTA
										* EGS4.DL3[L - 1][EGS4.MEDIUM - 1]);
						phi2 = EGS4.DL1[L1 - 1][EGS4.MEDIUM - 1]
								+ DELTA
								* (EGS4.DL2[L1 - 1][EGS4.MEDIUM - 1] + DELTA
										* EGS4.DL3[L1 - 1][EGS4.MEDIUM - 1]);
					} else {
						phi1 = EGS4.DL4[L - 1][EGS4.MEDIUM - 1]
								+ EGS4.DL5[L - 1][EGS4.MEDIUM - 1]
								* Math.log(DELTA
										+ EGS4.DL6[L - 1][EGS4.MEDIUM - 1]);
						phi2 = phi1;
					}
					REJF = (1.0 + aux * aux) * phi1 - 2.0 * aux * phi2 / 3.0;

					if (RNNO07 < REJF)
						break;
				}// UNTIL (rnno07 < rejf);
			}

			// "   SET UP THE NEW PHOTON                                          "
			EGS4.NP = EGS4.NP + 1;
			if (EGS4.NP > EGS4.$MXSTACK) {
				// OUTPUT np; (//' Stack overflow in BREMS! np = ',i6,
				// ' Increase $MXSTACK and try again'//);
				// stop;
				EGS4.STOPPROGRAM = true;
				EGS4.seqStr = " ***************************************************"
						+ "  \n"
						+ " Stack overflow in BREMS! np = "
						+ EGS4.NP
						+ "  \n" + " Increase $MXSTACK and try again ";// +"  \n";
				// if(EGS4.iprint>2)
				eq.printSequence(EGS4.seqStr);

				return;// stop;
			}
			EGS4.E[EGS4.NP - 1] = PESG;// //pesg;
			EGS4.IQ[EGS4.NP - 1] = 0;
			// $TRANSFER PROPERTIES TO (np) FROM (NPold);
			// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
			// $TRANSFER PROPERTIES TO (np) FROM (np-1);
			// ==================================================
			EGS4.X[EGS4.NP - 1] = EGS4.X[EGS4.NP - 2];
			EGS4.Y[EGS4.NP - 1] = EGS4.Y[EGS4.NP - 2];
			EGS4.Z[EGS4.NP - 1] = EGS4.Z[EGS4.NP - 2];
			EGS4.IR[EGS4.NP - 1] = EGS4.IR[EGS4.NP - 2];
			EGS4.WT[EGS4.NP - 1] = EGS4.WT[EGS4.NP - 2];
			EGS4.DNEAR[EGS4.NP - 1] = EGS4.DNEAR[EGS4.NP - 2];
			EGS4.LATCH[EGS4.NP - 1] = EGS4.LATCH[EGS4.NP - 2];

			// EGS4.WT[EGS4.NP-1] = EGS4.WT[EGS4.NP-1]/EGS4.nbr_split;
			EGS4.WT[EGS4.NP - 1] = weight;// =EGS4.WT[EGS4.NP-1]/EGS4.nbr_split

			if (EGS4.ibrdst < 0) {// "The photon will inherit the direction from "
									// "the electron. This option is given so that "
									// "the user can implement their own brems angle "
									// "schemes via a call to eq.AUSGAB"
				EGS4.U[EGS4.NP - 1] = EGS4.U[EGS4.NPold - 1];
				EGS4.V[EGS4.NP - 1] = EGS4.V[EGS4.NPold - 1];
				EGS4.W[EGS4.NP - 1] = EGS4.W[EGS4.NPold - 1];
			} else {
				if (EGS4.ibrdst == 1) {
					// "/ *"
					// This is the original implementation
					// suggested by Alex Bielajew. Commented out as
					// the implementation below is way more efficient.
					// IK, Sep. 2004.
					// ttese = ese/rm; esedei = ttese/tteie;
					// rjarg1 = 1+esedei*esedei;
					// rjarg2 = 3*rjarg1 - 2*esedei;
					// rjarg3 = ((1-esedei)/(2*tteie*esedei))**2;
					// $SET-BREM-REJECTION-FUNCTION(REJMIN,0.0);
					// $SET-BREM-REJECTION-FUNCTION(REJMID,1.0);
					// $SET-BREM-REJECTION-FUNCTION(REJMAX,y2max);
					// rejtop = max(rejmin,rejmid,rejmax);
					// LOOP [
					// $RANDOMSET y2tst; y2tst = y2tst/(1-y2tst+y2maxi);
					// $SET-BREM-REJECTION-FUNCTION(REJTST,Y2TST);
					// $RANDOMSET rtest;
					// ] UNTIL (rtest*rejtop <= REJTST);
					// " * //
					ttese = ESE / EGS4.RM;
					esedei = ttese / tteie;
					rjarg1 = 1.0 + esedei * esedei;
					rjarg2 = rjarg1 + 2 * esedei;
					aux = 2.0 * ESE * tteie / ESG;
					aux = aux * aux;
					aux1 = aux * ztarg;
					if (aux1 > 10.0) {
						rjarg3 = EGS4.LZBRANG[EGS4.MEDIUM - 1] + (1.0 - aux1)
								/ (aux1 * aux1);
					} else {
						rjarg3 = Math.log(aux / (1.0 + aux1));
					}
					rejmax = rjarg1 * rjarg3 - rjarg2;
					// LOOP [
					while (true) {
						y2tst = EGS4.random01();
						rtest = EGS4.random01();
						aux3 = z2maxi / (y2tst + (1.0 - y2tst) * z2maxi);
						rtest = rtest * aux3 * rejmax;
						y2tst = aux3 * aux3 - 1.0;
						y2tst1 = esedei * y2tst / Math.pow(aux3, 4.0);
						aux4 = 16.0 * y2tst1 - rjarg2;
						aux5 = rjarg1 - 4.0 * y2tst1;
						if (rtest < aux4 + aux5 * rjarg3)
							break;
						aux2 = Math.log(aux
								/ (1.0 + aux1 / Math.pow(aux3, 4.0)));
						rejtst = aux4 + aux5 * aux2;

						if (rtest < rejtst)
							break;
					}// UNTIL (rtest < rejtst );

				} else {
					y2tst = EGS4.random01();
					y2tst = y2tst / (1.0 - y2tst + y2maxi);
				}
				EGS4.COSTHE = 1.0 - 2.0 * y2tst * y2maxi;
				EGS4.SINTHE = Math.sqrt(Math.max((1.0 - EGS4.COSTHE)
						* (1.0 + EGS4.COSTHE), 0.0));
				// $SELECT-AZIMUTHAL-ANGLE(cphi,sphi);
				do {
					// $RANDOMSET xphi;
					xphi = EGS4.random01();
					xphi = 2.0 * xphi - 1.0;
					xphi2 = xphi * xphi;
					// $RANDOMSET yphi;
					yphi = EGS4.random01();
					yphi2 = yphi * yphi;
					rhophi2 = xphi2 + yphi2;
				} while (rhophi2 > 1.0);
				rhophi2 = 1.0 / rhophi2;
				cphi = (xphi2 - yphi2) * rhophi2;
				sphi = 2.0 * xphi * yphi * rhophi2;

				if (sinpsi >= 1.e-10) {
					us = EGS4.SINTHE * cphi;
					vs = EGS4.SINTHE * sphi;
					EGS4.U[EGS4.NP - 1] = c * cosdel * us - sindel * vs + a
							* EGS4.COSTHE;
					EGS4.V[EGS4.NP - 1] = c * sindel * us + cosdel * vs + b
							* EGS4.COSTHE;
					EGS4.W[EGS4.NP - 1] = c * EGS4.COSTHE - sinpsi * us;
				} else {
					EGS4.U[EGS4.NP - 1] = EGS4.SINTHE * cphi;
					EGS4.V[EGS4.NP - 1] = EGS4.SINTHE * sphi;
					EGS4.W[EGS4.NP - 1] = EGS4.COSTHE;
				}
			}
		}

		EGS4.E[EGS4.NPold - 1] = PESE;// e(npold) = pese;

		return;
		// "END OF SUBROUTINE BREMS" END;

	}

	// "******************************************************************"
	// "                               National Research Council of Canada"
	/**
	 * Called by ELECTR. DISCRETE MOLLER SCATTERING HAS BEEN ARBITRARILY DEFINED AND CALCULATED 
	 * TO MEAN MOLLER SCATTERINGS WHICH IMPART TO THE SECONDARY ELECTRON SUFFICIENT ENERGY THAT 
	 * IT BE TRANSPORTED DISCRETELY.  THE THRESHOLD TO TRANSPORT AN ELECTRON DISCRETELY IS A TOTAL ENERGY OF AE OR A KINETIC ENERGY 
	 * OF TE=AE-RM.  SINCE THE KINETIC ENERGY TRANSFER IS ALWAYS, BY DEFINITION, LESS THAN HALF OF THE INCIDENT KINETIC ENERGY, THIS 
	 * IMPLIES THAT THE INCIDENT ENERGY, EIE, MUST BE LARGER THAN THMOLL=TE*2+RM.  THE REST OF THE COLLISION CONTRIBUTION IS 
	 * SUBTRACTED CONTINUOUSLY FROM THE ELECTRON AS IONIZATION LOSS DURING TRANSPORT.
	 * 
	*/
	protected static void MOLLER() {
		// "                                                                  "
		// "******************************************************************"
		// "   DISCRETE MOLLER SCATTERING (A CALL TO THIS ROUTINE) HAS BEEN   "
		// "   ARBITRARILY DEFINED AND CALCULATED TO MEAN MOLLER SCATTERINGS  "
		// "   WHICH IMPART TO THE SECONDARY ELECTRON SUFFICIENT ENERGY THAT  "
		// "   IT BE TRANSPORTED DISCRETELY.  THE THRESHOLD TO TRANSPORT AN   "
		// "   ELECTRON DISCRETELY IS A TOTAL ENERGY OF AE OR A KINETIC ENERGY"
		// "   OF TE=AE-RM.  SINCE THE KINETIC ENERGY TRANSFER IS ALWAYS, BY  "
		// "   DEFINITION, LESS THAN HALF OF THE INCIDENT KINETIC ENERGY, THIS"
		// "   IMPLIES THAT THE INCIDENT ENERGY, EIE, MUST BE LARGER THAN     "
		// "   THMOLL=TE*2+RM.  THE REST OF THE COLLISION CONTRIBUTION IS     "
		// "   SUBTRACTED CONTINUOUSLY FROM THE ELECTRON AS IONIZATION        "
		// "   LOSS DURING TRANSPORT.                                         "
		// "******************************************************************"

		// ; Copyright NRC;

		// $COMIN-MOLLER; "DEFAULT REPLACEMENT PRODUCES THE FOLLOWING:
		// "  COMIN/EGS-VARIANCE-REDUCTION, DEBUG,STACK,THRESH,"
		// "UPHIOT,USEFUL,RANDOM/;"
		// ;COMIN/EII-DATA,ELECIN,EPCONT,EDGE,BREMPR/;

		// $DEFINE-LOCAL-VARIABLES-MOLLER;
		// REPLACE {$DEFINE-LOCAL-VARIABLES-MOLLER;} WITH
		// {;
		// "Local MOLLER variables in order of their appearance"

		double PEIE = 0.0;// "precise total energy of incident electron"
		double PEKSE2 = 0.0;// "precise kinetic energy of 2nd secondary electron"
		double PESE1 = 0.0;// "precise total energy of 1st secondary electron"
		double PESE2 = 0.0;// "precise total energy of 2nd secondary electron"
		double PEKIN = 0.0;// "precise kinetic energy of incident electron"
		double H1 = 0.0;// "used for polar scattering angle calculation"
		double DCOSTH = 0.0;// "polar scattering angle squared"
		double EIE = 0.0;// "total energy of incident electron"
		double EKIN = 0.0;// "kinetic energy of incident electron"
		double T0 = 0.0;// "kinetic energy of incident electron in units of RM"
		double E0 = 0.0;// "total energy of incident electron in units of RM"
		double EXTRAE = 0.0;// "energy above the Moller threshold"
		double E02 = 0.0;// "E0**2"
		//double EP0 = 0.0;// "minimum alowed kinetic energy fraction"
		double G2 = 0.0;//
		double G3 = 0.0;// "used for rejection function calculation"
		double GMAX = 0.0;// "maximum value of the rejection function"
		double BR = 0.0;// "kinetic energy fraction to lowew energy electron"
		double R = 0.0;// "(1-BR)/BR"
		double REJF4 = 0.0;// "rejection function"
		double RNNO27 = 0.0;// "random number for BR sampling"
		double RNNO28 = 0.0;// "random number for rejection"
		//double ESE1 = 0.0;// "energy of 1st secondary electron"
		//double ESE2 = 0.0;// "energy of 2nd secondary electron"
		// }

		double sigm = 0.0;
		double pbrem = 0.0;
		double rsh = 0.0;
		double Uj = 0.0;
		double sig_j = 0.0;//
		int lelke = 0;
		int iele = 0;
		int ish = 0;
		int nsh = 0;
		int ifirst = 0;
		int i = 0;
		int jj = 0;
		int iZ = 0;

		// "IRCODE=1;  appears to be unused, IK Oct 97"
		EGS4.NPold = EGS4.NP;// "Set the old stack counter"
		PEIE = EGS4.E[EGS4.NP - 1];// "PRECISE ENERGY OF INCIDENT ELECTRON"
		EIE = PEIE; // "ENERGY OF INCIDENT ELECTRON"
		PEKIN = PEIE - EGS4.PRM; // "PRECISE K.E. OF INCIDENT ELECTRON"
		EKIN = PEKIN;

		if ((EGS4.eii_flag > 0) && (EGS4.eii_nsh[EGS4.MEDIUM - 1] > 0)) {
			// "The EII flag is set and this medium has shells for which we want to"
			// "simulate EII => sample if the interaction is with a EII shell"
			// $SET INTERVAL elke,eke;
			Double dbl = new Double(EGS4.EKE1[EGS4.MEDIUM - 1] * EGS4.ELKE
					+ EGS4.EKE0[EGS4.MEDIUM - 1]);
			lelke = dbl.intValue();

			// $EVALUATE sigm USING esig(elke);
			sigm = EGS4.ESIG1[lelke - 1][EGS4.MEDIUM - 1] * EGS4.ELKE
					+ EGS4.ESIG0[lelke - 1][EGS4.MEDIUM - 1];
			// $EVALUATE pbrem USING ebr1(elke);
			pbrem = EGS4.EBR11[lelke - 1][EGS4.MEDIUM - 1] * EGS4.ELKE
					+ EGS4.EBR10[lelke - 1][EGS4.MEDIUM - 1];
			sigm = sigm * (1.0 - pbrem);
			rsh = EGS4.random01();
			rsh = sigm * rsh;
			for (iele = 1; iele <= EGS4.NNE[EGS4.MEDIUM - 1]; iele++) {
				Double dbl1 = new Double(
						EGS4.ZELEM[EGS4.MEDIUM - 1][iele - 1] + 0.5);
				iZ = dbl1.intValue();// int(zelem(medium,iele)+0.5);
				nsh = EGS4.eii_no[EGS4.MEDIUM - 1][iele - 1];// (medium,iele);
				if (nsh > 0) {
					ifirst = EGS4.eii_first[EGS4.MEDIUM - 1][iele - 1];// (medium,iele);
					for (ish = 1; ish <= nsh; ish++) {
						Uj = EGS4.binding_energies[ish - 1][iZ - 1];// (ish,iZ);
						if ((EKIN > Uj)
								&& ((Uj > EGS4.TE[EGS4.MEDIUM - 1]) || (Uj > EGS4.AP[EGS4.MEDIUM - 1]))) {
							jj = ifirst + ish - 1;
							Double dbll = new Double(EGS4.eii_a[jj - 1]
									* EGS4.ELKE + EGS4.eii_b[jj - 1] + (jj - 1)
									* EGS4.$N_EII_BINS);
							i = dbll.intValue();
							sig_j = EGS4.eii_xsection_a[i - 1] * EGS4.ELKE
									+ EGS4.eii_xsection_b[i - 1];
							sig_j = sig_j * EGS4.PZ[EGS4.MEDIUM - 1][iele - 1]
									* EGS4.eii_cons[EGS4.MEDIUM - 1];
							rsh = rsh - sig_j;
							if (rsh < 0) {
								eii_sample(ish, iZ, Uj);
								return;
							}
						}
					}
				}
			}
		}

		if (EKIN <= 2.0 * EGS4.TE[EGS4.MEDIUM - 1])
			return;
		T0 = EKIN / EGS4.RM;
		E0 = T0 + 1.0;
		EXTRAE = EIE - EGS4.THMOLL[EGS4.MEDIUM - 1];
		E02 = E0 * E0;
		// "BETAI2=E02/(E02-1.0); "
		// "BLIF 96/2/1 -- not needed for Moller fix-up"
		//EP0 = EGS4.TE[EGS4.MEDIUM - 1] / EKIN;
		// "G1=(1.-2.*EP0)*BETAI2;"
		// "BLIF 96/2/1 -- not needed for Moller fix-up"
		G2 = T0 * T0 / E02;
		G3 = (2. * T0 + 1.) / E02;
		// "   H.H.NAGEL HAS CONSTRUCTED A FACTORIZATION OF THE FREQUENCY"
		// "   DISTRIBUTION FUNCTION FOR THE MOLLER DIFFERENTIAL CROSS"
		// "   SECTION USED AS SUGGESTED BY BUTCHER AND MESSEL."
		// "   (H.H.NAGEL, OP.CIT., P. 53-55)                                 "
		// "   HOWEVER, A MUCH SIMPLER SAMPLING METHOD WHICH DOES NOT BECOME  "
		// "   VERY INEFFICIENT NEAR THMOLL IS THE FOLLOWING. . .             "
		// "   LET BR=EKS/EKIN,  WHERE EKS IS KINETIC ENERGY TRANSFERED TO THE"
		// "   SECONDARY ELECTRON AND EKIN IS THE INCIDENT KINETIC ENERGY.    "

		// "   MODIFIED (7 FEB 1974) TO USE THE TRUE MOLLER CROSS SECTION."
		// "   THAT IS, INSTEAD OF THE E+ E- AVERAGE GIVEN IN THE ROSSI"
		// "   FORMULA USED BY NAGEL.  THE SAMPLING SCHEME IS THAT"
		// "   USED BY MESSEL AND CRAWFORD (EPSDF 1970 P.13)"
		// "   FIRST SAMPLE (1/BR**2) OVER (TE/EKIN,1/2) . . .            "

		GMAX = (1. + 1.25 * G2); // "BLIF 96/2/1 -- Moller fix-up"
		while (true)// LOOP[" TO RETRY IF REJECTED"
		{
			RNNO27 = EGS4.random01();
			BR = EGS4.TE[EGS4.MEDIUM - 1] / (EKIN - EXTRAE * RNNO27);

			// "   USE MESSEL AND CRAWFORDS REJECTION FUNCTION."
			R = BR / (1. - BR);
			RNNO28 = EGS4.random01();
			REJF4 = (1. + G2 * BR * BR + R * (R - G3));// "G1*"(1.+G2*BR*BR+R*(R-G3));
														// //"BLIF 96/2/1 -- Moller fix-up"
			RNNO28 = GMAX * RNNO28; // "BLIF 96/2/1 -- Moller fix-up"

			if (RNNO28 <= REJF4)
				break;
		}// UNTIL RNNO28.LE.REJF4; "TRY UNTIL ACCEPTED. END REJECTION LOOP"

		PEKSE2 = BR * EKIN; // "PRECISE KINETIC ENERGY OF SECONDARY ELECTRON #2"
		PESE1 = PEIE - PEKSE2; // "PRECISE ENERGY OF SECONDARY ELECTRON #1"
		PESE2 = PEKSE2 + EGS4.PRM; // "PRECISE ENERGY OF SECONDARY ELECTRON #2"
		//ESE1 = PESE1; // "ENERGY OF SECONDARY ELECTRON 1"
		//ESE2 = PESE2; // "ENERGY OF SECONDARY ELECTRON 2"
		EGS4.E[EGS4.NP - 1] = PESE1;
		// $CHECK-STACK(np+1,'MOLLER');
		if (EGS4.NP + 1 > EGS4.$MXSTACK) {
			EGS4.STOPPROGRAM = true;
			EGS4.seqStr = " ***************************************************"
					+ "  \n"
					+ " In subroutine "
					+ "MOLLER"
					+ " stack size exceeded!"
					+ "  \n"
					+ " $MXSTACK = "
					+ EGS4.$MXSTACK
					+ " np = "
					+ EGS4.NP
					+ "  \n"
					+ " Increase $MXSTACK and try again "
					+ "  \n"
					+ " Terminating execution "
					+ "  \n"
					+ " ***************************************************";// +"  \n";
																				// ;
			// if(EGS4.iprint>2)
			eq.printSequence(EGS4.seqStr);

			return;// stop;
		}

		EGS4.E[EGS4.NP] = PESE2;// E(NP+1)=PESE2;
		// "   SINCE BR.LE.0.5, E(NP+1) MUST BE .LE. E(NP). "
		// "   MOLLER ANGLES ARE UNIQUELY DETERMINED BY KINEMATICS            "

		// " One possible way of dealing with double counting of angular      "
		// " deflections in inelastic scattering would be to                  "
		// " not deflect the 'old' electron as these deflections are          "
		// " already taken into account in the multiple elastic scattering    "
		// " This approach has the disadvantage of loosing correlations       "
		// " between big energy losses and strong angular deflections         "
		// " The advantage of such an approach is its simplicity.             "
		// " If spin effects for multiple elastic scattering are turned on,   "
		// " the double counting is taken into account by the appropriate     "
		// " modification of the scattering power (which depends on AE)       "
		// "                                                                  "
		// "                                                                  "
		// " IK, June 1999                                                    "

		H1 = (PEIE + EGS4.PRM) / PEKIN;
		// "   DIRECTION COSINE CHANGE FOR 'OLD' ELECTRON                     "
		DCOSTH = H1 * (PESE1 - EGS4.PRM) / (PESE1 + EGS4.PRM);
		EGS4.SINTHE = Math.sqrt(1.0 - DCOSTH);
		EGS4.COSTHE = Math.sqrt(DCOSTH);

		// "sinthe = 0; costhe = 1; <- this will turn off the Moller ang. deflections"

		UPHI(2, 1);

		// "   RELATED CHANGE AND (X,Y,Z) SETUP FOR 'NEW' ELECTRON            "
		EGS4.NP = EGS4.NP + 1;
		EGS4.IQ[EGS4.NP - 1] = -1;
		DCOSTH = H1 * (PESE2 - EGS4.PRM) / (PESE2 + EGS4.PRM);
		EGS4.SINTHE = -Math.sqrt(1.0 - DCOSTH);
		EGS4.COSTHE = Math.sqrt(DCOSTH);
		UPHI(3, 2);

		return;
	}// "END OF SUBROUTINE MOLLER" END;

	/**
	 * Called by MOLLER. Sample EII (Electron Impact Ionization).
	 * @param ish ish
	 * @param iZ iZ
	 * @param Uj Uj
	 */
	protected static void eii_sample(int ish, int iZ, double Uj) {
		// $INTEGER ish,iZ;
		// $REAL Uj;

		// $COMIN-EII-SAMPLE;

		double T = 0.0;
		double tau = 0.0;
		double tau1 = 0.0;
		double tau12 = 0.0;
		double tau2 = 0.0;
		double p2 = 0.0;
		double beta2 = 0.0;
		double c1 = 0.0;
		double c2 = 0.0;
		double Wmax = 0.0;
		double xmax = 0.0;
		double fm_s = 0.0;
		double fm_h = 0.0;
		double prob_s = 0.0;
		double prob = 0.0;
		double r1 = 0.0;
		double r2 = 0.0;
		double r3 = 0.0;
		double wx = 0.0;
		double wxx = 0.0;
		double aux = 0.0;
		double frej = 0.0;
		double peie = 0.0;
		double pese1 = 0.0;
		double pese2 = 0.0;
		double dcosth = 0.0;
		double h1 = 0.0;
		int iarg = 0;
		//double eta = 0.0;
		//double cphi = 0.0;
		//double sphi = 0.0;
		//int np_save = 0;
		//int ip = 0;
		//int j = 0;
		// $DEFINE-VARIABLES-FOR-SELECT-AZIMUTHAL-ANGLE;
		// REPLACE {$DEFINE-VARIABLES-FOR-SELECT-AZIMUTHAL-ANGLE;} WITH {;
		//double xphi = 0.0;
		//double xphi2 = 0.0;
		//double yphi = 0.0;
		//double yphi2 = 0.0;
		//double rhophi2 = 0.0;

		// / *
		// real*8 sum_es, sum_es2, sum_eh, sum_eh2;
		// $REAL ss_0, sh_0, ss_1, sh_1, uwm;
		// $INTEGER count_s,iloop;
		// $LOGICAL is_soft;
		// * /

		// calculate some useful constants //
		peie = EGS4.E[EGS4.NP - 1];
		T = peie - EGS4.RM;
		tau = T / EGS4.RM;
		tau1 = tau + 1.0;
		tau12 = tau1 * tau1;
		tau2 = tau * tau;
		p2 = tau2 + 2.0 * tau;
		beta2 = p2 / tau12;
		// "c1 = tau2/tau12; "
		Wmax = 0.5 * (T + Uj);
		xmax = Uj / Wmax;
		c1 = (Wmax / peie) * (Wmax / peie);
		c2 = (2.0 * tau + 1.0) / tau12;
		fm_s = Math.log(EGS4.RMT2 * p2 / Uj) - beta2 - 0.5;
		prob_s = 0.66666667 * fm_s * (1 + xmax + xmax * xmax);
		// "fm_h = 1 + c1 - c2;"
		fm_h = 2.0 + c1 - c2;
		if (fm_h < 1.0)
			fm_h = 1.0;
		prob = fm_h + prob_s;

		// / *
		// /sum_es,sum_es2,sum_eh,sum_eh2,count_s/=0;
		// * /

		// LOOP [
		while (true) {
			r1 = EGS4.random01();
			r2 = EGS4.random01();
			r3 = EGS4.random01();
			if (r1 * prob < fm_h) {// "Use the hard collision cross section "
				wx = 1.0 / (r2 * xmax + 1.0 - r2);
				wxx = wx * xmax;
				aux = wxx / (2.0 - wxx);
				frej = (1.0 + aux * (aux - c2) + c1 * wxx * wxx) / fm_h;
			} else {// "Use the soft collision cross section "
				wx = 1.0 / Math.pow(r2 * xmax * xmax * xmax + 1.0 - r2,
						0.333333333);
				frej = 1.0 - Math.log(wx) / fm_s;
			}
			if (r3 < frej)
				break;
		} // UNTIL ( r3 < frej );

		wx = wx * Uj;
		// / *
		// IF( is_soft ) [
		// sum_es = sum_es + wx; sum_es2 = sum_es2 + wx*wx; count_s = count_s +
		// 1;
		// ]
		// ELSE [
		// sum_eh = sum_eh + wx; sum_eh2 = sum_eh2 + wx*wx;
		// ]
		// * /

		// / *
		// iloop = 10000000;
		// sum_es = sum_es/count_s; sum_es2 = sum_es2/count_s;
		// sum_es2 = sum_es2 - sum_es*sum_es;
		// IF( sum_es2 > 0 ) sum_es2 = sqrt(sum_es2/(count_s-1));
		// sum_eh = sum_eh/(iloop-count_s); sum_eh2 = sum_eh2/(iloop-count_s);
		// sum_eh2 = sum_eh2 - sum_eh*sum_eh;
		// IF( sum_eh2 > 0 ) sum_eh2 = sqrt(sum_eh2/(iloop-count_s-1));
		// write(6,*) 'soft: ',count_s,dble(count_s)/dble(iloop);
		// write(6,*) '<Es>: ',sum_es,' +/- ',sum_es2;
		// write(6,*) '<Eh>: ',sum_eh,' +/- ',sum_eh2;
		// p2 = 2*rm*tau*(tau+2); uwm = Uj/Wmax;
		// ss_0 = 2*(log(p2/Uj)-uwm**3*log(p2/Wmax)-
		// (beta2+0.833333)*(1-uwm**3))/3/Uj;
		// sh_0 = ((1-uwm)*(1+uwm/(2-uwm))+Uj*(Wmax-Uj)/e(np)**2
		// - (2*tau+1)/(tau+1)**2*uwm/2*log((2-uwm)/uwm))/Uj;
		// ss_1 = log(p2/Uj)-uwm**2*log(p2/Wmax)-
		// (beta2+1)*(1-uwm**2);
		// sh_1 = log(Wmax/Uj/(2-uwm))+2*(Wmax-Uj)/(2*Wmax-Uj)
		// +(Wmax**2-Uj**2)/e(np)**2/2
		// -(2*tau+1)/(tau+1)**2*log((2*Wmax-Uj)/Wmax);
		// write(6,*) 'expected: ',ss_0/(ss_0 + sh_0);
		// write(6,*) '<Es>: ',ss_1/ss_0;
		// write(6,*) '<Eh>: ',sh_1/sh_0;

		// stop;
		// *//

		// set-up new particles /
		h1 = (peie + EGS4.PRM) / T;
		pese1 = peie - wx;
		EGS4.E[EGS4.NP - 1] = pese1;
		dcosth = h1 * (pese1 - EGS4.PRM) / (pese1 + EGS4.PRM);
		EGS4.SINTHE = Math.sqrt(1.0 - dcosth);
		EGS4.COSTHE = Math.sqrt(dcosth);
		UPHI(2, 1);

		pese2 = wx - Uj + EGS4.PRM;
		if (pese2 > EGS4.AE[EGS4.MEDIUM - 1]) {
			// $CHECK-STACK(np+1,'eii_sample');
			if (EGS4.NP + 1 > EGS4.$MXSTACK) {
				EGS4.STOPPROGRAM = true;
				EGS4.seqStr = " ***************************************************"
						+ "  \n"
						+ " In subroutine "
						+ "eii_sample"
						+ " stack size exceeded!"
						+ "  \n"
						+ " $MXSTACK = "
						+ EGS4.$MXSTACK
						+ " np = "
						+ EGS4.NP
						+ "  \n"
						+ " Increase $MXSTACK and try again "
						+ "  \n"
						+ " Terminating execution "
						+ "  \n"
						+ " ***************************************************";// +"  \n";;
				// if(EGS4.iprint>2)
				eq.printSequence(EGS4.seqStr);

				return;// stop;
			}

			EGS4.NP = EGS4.NP + 1;
			EGS4.E[EGS4.NP - 1] = pese2;
			dcosth = h1 * (pese2 - EGS4.PRM) / (pese2 + EGS4.PRM);
			EGS4.SINTHE = -Math.sqrt(1.0 - dcosth);
			EGS4.COSTHE = Math.sqrt(dcosth);
			EGS4.IQ[EGS4.NP - 1] = -1;
			UPHI(3, 2);
			EGS4.EDEP = 0.0;
		} else {
			EGS4.EDEP = wx - Uj;
		}

		RELAX(Uj, ish, iZ);

		if (EGS4.EDEP > 0.0) {
			// $AUSCALL($PHOTXAUS);
			iarg = EGS4.$PHOTXAUS;
			if (EGS4.iausfl[iarg] != 0) {
				eq.AUSGAB(iarg);
			}

		}
		/*
		 * if( EGS4.nbr_split > 1 ) { np_save = EGS4.NP;
		 * for(ip=EGS4.NPold+1;ip<=np_save;ip++) { if( EGS4.IQ[ip-1] == 0 ) {
		 * EGS4.WT[ip-1] = EGS4.WT[ip-1]/EGS4.nbr_split; for(j =
		 * 1;j<=EGS4.nbr_split-1;j++) { EGS4.NP = EGS4.NP + 1;
		 * //$CHECK-STACK(np,'eii_sample'); if( EGS4.NP > EGS4.$MXSTACK ) {
		 * EGS4.STOPPROGRAM=true;
		 * EGS4.STOPPROGRAMs=" ***************************************************"
		 * +"  \n"; EGS4.STOPPROGRAMs=EGS4.STOPPROGRAMs+
		 * " In subroutine "+"eii_sample"+" stack size exceeded!"+"  \n";
		 * EGS4.STOPPROGRAMs=EGS4.STOPPROGRAMs+
		 * " $MXSTACK = "+EGS4.$MXSTACK+" np = "+EGS4.NP+"  \n";
		 * EGS4.STOPPROGRAMs=EGS4.STOPPROGRAMs+
		 * " Increase $MXSTACK and try again "+"  \n";
		 * EGS4.STOPPROGRAMs=EGS4.STOPPROGRAMs+
		 * " Terminating execution "+"  \n";
		 * EGS4.STOPPROGRAMs=EGS4.STOPPROGRAMs+
		 * " ***************************************************"+"  \n";
		 * //System.out.println(EGS4.STOPPROGRAMs); return;//stop; }
		 * 
		 * EGS4.IQ[EGS4.NP-1] = 0; EGS4.E[EGS4.NP-1] = EGS4.E[ip-1]; //$TRANSFER
		 * PROPERTIES TO (np) FROM (ip); EGS4.X[EGS4.NP-1]=EGS4.X[ip-1];
		 * EGS4.Y[EGS4.NP-1]=EGS4.Y[ip-1]; EGS4.Z[EGS4.NP-1]=EGS4.Z[ip-1];
		 * EGS4.IR[EGS4.NP-1]=EGS4.IR[ip-1]; EGS4.WT[EGS4.NP-1]=EGS4.WT[ip-1];
		 * EGS4.DNEAR[EGS4.NP-1]=EGS4.DNEAR[ip-1];
		 * EGS4.LATCH[EGS4.NP-1]=EGS4.LATCH[ip-1];
		 * 
		 * eta=EGS4.random01(); eta = 2.0*eta - 1.0; EGS4.W[EGS4.NP-1] = eta;
		 * eta = (1.0-eta)*(1.0+eta); if( eta > 1.e-20 ) { eta = Math.sqrt(eta);
		 * //$SELECT-AZIMUTHAL-ANGLE(cphi,sphi); do { //$RANDOMSET xphi;
		 * xphi=EGS4.random01(); xphi = 2.0*xphi - 1.0; xphi2 = xphi*xphi;
		 * //$RANDOMSET yphi; yphi=EGS4.random01(); yphi2 = yphi*yphi; rhophi2 =
		 * xphi2 + yphi2; } while(rhophi2 > 1.0); rhophi2 = 1.0/rhophi2; cphi =
		 * (xphi2 - yphi2)*rhophi2; sphi = 2.0*xphi*yphi*rhophi2;
		 * 
		 * EGS4.U[EGS4.NP-1] = eta*cphi; EGS4.V[EGS4.NP-1] = eta*sphi; } else {
		 * EGS4.U[EGS4.NP-1] = 0.0; EGS4.V[EGS4.NP-1] = 0.0; EGS4.W[EGS4.NP-1] =
		 * 1.0; } } } } }
		 */
		return;// end;
	}

	// "******************************************************************"
	// "                               National Research Council of Canada"
	/**
	 * Called by ELECTR. DISCRETE BHABHA SCATTERING HAS BEEN ARBITRARILY DEFINED AND CALCULATED TO MEAN BHABHA SCATTERINGS 
	 * WHICH IMPART TO THE SECONDARY ELECTRON SUFFICIENT ENERGY THAT IT BE TRANSPORTED DISCRETELY, I.E. E=AE OR T=TE. 
	 * IT IS NOT GUARANTEED THAT THE FINAL POSITRON WILL HAVE THIS MUCH ENERGY HOWEVER.  THE EXACT BHABHA DIFFERENTIAL CROSS SECTION IS USED.
	 */
	protected static void BHABHA() {
		// "                                                                  "
		// "******************************************************************"
		// "   DISCRETE BHABHA SCATTERING (A CALL TO THIS ROUTINE) HAS BEEN   "
		// "   ARBITRARILY DEFINED AND CALCULATED TO MEAN BHABHA SCATTERINGS  "
		// "   WHICH IMPART TO THE SECONDARY ELECTRON SUFFICIENT ENERGY THAT  "
		// "   IT BE TRANSPORTED DISCRETELY, I.E. E=AE OR T=TE.  IT IS NOT    "
		// "   GUARANTEED THAT THE FINAL POSITRON WILL HAVE THIS MUCH ENERGY  "
		// "   HOWEVER.  THE EXACT BHABHA DIFFERENTIAL CROSS SECTION IS USED. "
		// "******************************************************************"

		// ; Copyright NRC;

		// $COMIN-BHABHA; "DEFAULT REPLACEMENT PRODUCES THE FOLLOWING:      "
		// "COMIN/DEBUG,EGS-VARIANCE-REDUCTION,STACK,"
		// "THRESH,UPHIOT,USEFUL,RANDOM/;"

		// $DEFINE-LOCAL-VARIABLES-BHABHA;
		// REPLACE {$DEFINE-LOCAL-VARIABLES-BHABHA;} WITH
		// {;
		// "Local variables in order of appearance"
		double PEIP = 0.0;// "precise total energy of incident positron"
		double PEKIN = 0.0;// "precise kinetic energy of incident positron"
		double PEKSE2 = 0.0;// "precise kinetic energy of second 'electron'"
		double PESE1 = 0.0;// "precise total energy of first 'electron'"
		double PESE2 = 0.0;// "precise total energy of second 'electron'"
		double H1 = 0.0;// "used in direction cosine calculations"
		double DCOSTH = 0.0;// "polar scattering angle for more energetic 'electron'"
		//double EIP = 0.0;// "total energy of incident positron"
		double EKIN = 0.0;// "kinetic energy of incident positron"
		double T0 = 0.0;// "kinetic energy of incident positron in units of RM"
		double E0 = 0.0;// "total energy of incident positron in units of RM"
		double E02 = 0.0;// "E0**2"
		double YY = 0.0;// "1/(T0+2)"
		double Y2 = 0.0;//
		double YP = 0.0;//
		double YP2 = 0.0;// "various functions of YY"
		double BETA2 = 0.0;// "incident positron velocity in units of c"
		double EP0 = 0.0;// "minimum fractional energy of a secondary 'electron'"
		double EP0C = 0.0;// "1-EP0"
		double B1 = 0.0;//
		double B2 = 0.0;//
		double B3 = 0.0;//
		double B4 = 0.0;// "used in rejection function calculation"
		double RNNO03 = 0.0;//
		double RNNO04 = 0.0;// "random numbers"
		double BR = 0.0;// "kinetic energy fraction of the 2nd 'electron'"
		double REJF2 = 0.0;// "rejection function"
		//double ESE1 = 0.0;// "total energy of 1st 'electron'"
		//double ESE2 = 0.0;// "total energy of 2nd 'electron'"
		// }

		EGS4.NPold = EGS4.NP; // "Set the old stack counter"
		PEIP = EGS4.E[EGS4.NP - 1]; // "PRECISE ENERGY OF INCIDENT POSITRON"
		//EIP = PEIP; // "ENERGY OF INCIDENT POSITRON"
		PEKIN = PEIP - EGS4.PRM;// "PRECISE K.E. OF INCIDENT POSITRON"
		EKIN = PEKIN;
		T0 = EKIN / EGS4.RM;
		E0 = T0 + 1.;
		YY = 1. / (T0 + 2.);
		E02 = E0 * E0;
		// "BETAI2=E02/(E02-1.);" "BLIF 96/2/1 -- not needed for Bhabha fix-up"
		BETA2 = (E02 - 1.) / E02; // "BLIF 96/2/1 -- needed for Bhabha fix-up"
		EP0 = EGS4.TE[EGS4.MEDIUM - 1] / EKIN;
		EP0C = 1. - EP0;
		Y2 = YY * YY;
		YP = 1. - 2. * YY;
		YP2 = YP * YP;
		B4 = YP2 * YP;
		B3 = B4 + YP2;
		B2 = YP * (3. + Y2);
		B1 = 2. - Y2;
		// "   SAMPLE BR FROM MINIMUM(EP0) TO 1."
		// LOOP[
		while (true) {
			RNNO03 = EGS4.random01();
			BR = EP0 / (1. - EP0C * RNNO03);
			// "   APPLY REJECTION FUNCTION"
			RNNO04 = EGS4.random01();
			// "REJF2=EP0C*(BETAI2-BR*(B1-BR*(B2-BR*(B3-BR*B4))));BLIF 96/2/1 -- Bhabha fix-up"
			REJF2 = (1.0 - BETA2 * BR * (B1 - BR * (B2 - BR * (B3 - BR * B4)))); // "BLIF 96/2/1 -- Bhabha fix-up"
			if (RNNO04 <= REJF2)
				break;
		}// UNTIL RNNO04.LE.REJF2 ;
			// "   IF E- GOT MORE THAN E+, MOVE THE E+ POINTER AND REFLECT B"
			// $CHECK-STACK(np+1,'BHABHA');
		if (EGS4.NP + 1 > EGS4.$MXSTACK) {
			EGS4.STOPPROGRAM = true;

			EGS4.seqStr = " ***************************************************"
					+ "  \n"
					+ " In subroutine "
					+ "BHABHA"
					+ " stack size exceeded!"
					+ "  \n"
					+ " $MXSTACK = "
					+ EGS4.$MXSTACK
					+ " np = "
					+ EGS4.NP
					+ "  \n"
					+ " Increase $MXSTACK and try again "
					+ "  \n"
					+ " Terminating execution "
					+ "  \n"
					+ " ***************************************************";// +"  \n";
			// if(EGS4.iprint>2)
			eq.printSequence(EGS4.seqStr);

			return;// stop;
		}

		if (BR < 0.5) {
			EGS4.IQ[EGS4.NP] = -1;// IQ(NP+1)=-1;
		} else {
			EGS4.IQ[EGS4.NP - 1] = -1;// IQ(NP)=-1;
			EGS4.IQ[EGS4.NP] = 1;// IQ(NP+1)=1;
			BR = 1. - BR;
		}
		// "THE ABOVE PUTS E+ ON TOP OF STACK IF IT HAS LESS ENERGY"
		// "   DIVIDE UP THE ENERGY"
		BR = Math.max(BR, 0.0); // "AVOIDS POSSIBLE NEGATIVE NUMBER DUE TO ROUND-OFF"
		PEKSE2 = BR * EKIN; // "PRECISE KINETIC ENERGY OF SECONDARY 'ELECTRON' 2"
		PESE1 = PEIP - PEKSE2;// "PRECISE ENERGY OF SECONDARY 'ELECTRON' 1"
		PESE2 = PEKSE2 + EGS4.PRM; // "PRECISE ENERGY OF SECONDARY 'ELECTRON' 2"
		//ESE1 = PESE1;
		//ESE2 = PESE2;
		EGS4.E[EGS4.NP - 1] = PESE1;
		EGS4.E[EGS4.NP] = PESE2;
		// "   BHABHA ANGLES ARE UNIQUELY DETERMINED BY KINEMATICS"
		H1 = (PEIP + EGS4.PRM) / PEKIN;
		// "   DIRECTION COSINE CHANGE FOR 'OLD' ELECTRON"

		// "AFB modified the following statement 92/10/28 to avoid"
		// "numerical difficulties"
		// "DCOSTH=H1*(PESE1-PRM)/(PESE1+PRM);"
		DCOSTH = Math.min(1.0, H1 * (PESE1 - EGS4.PRM) / (PESE1 + EGS4.PRM));

		EGS4.SINTHE = Math.sqrt(1.0 - DCOSTH);
		EGS4.COSTHE = Math.sqrt(DCOSTH);
		UPHI(2, 1);
		EGS4.NP = EGS4.NP + 1;
		DCOSTH = H1 * (PESE2 - EGS4.PRM) / (PESE2 + EGS4.PRM);
		EGS4.SINTHE = -Math.sqrt(1.0 - DCOSTH);
		EGS4.COSTHE = Math.sqrt(DCOSTH);
		UPHI(3, 2);

		return;
	} // "END OF SUBROUTINE BHABHA" END;

	// "******************************************************************"
	// "                               National Research Council of Canada"
	/**
	 * Called by ELECTR. Handle two gamma from electron-positron anihhilation. 
	 * If the user requests radiative splitting (via nbr_split greater than 1), this routine produces 2*nbr_split annihilation photons at once, 
	 * each carying the fraction 1/nbr_split of the weight of the incident positron.
	 */
	protected static void ANNIH() {
		// "                                                                  "
		// "******************************************************************"
		// "   GAMMA SPECTRUM FOR TWO GAMMA IN-FLIGHT POSITRON ANNIHILATION.  "
		// "   USING SCHEME BASED ON HEITLER'S P269-27O FORMULAE.             "
		// "                                                                  "
		// "   If the user requests radiative splitting (via nbr_split > 1),  "
		// "   this routine produces 2*nbr_split annihilation photons at once,"
		// "   each carying the fraction 1/nbr_split of the weight of the     "
		// "   incident positron.                                             "
		// "                                                                  "
		// "   Except for taking out the calculation of                       "
		// "   LOG((1.0-EP0)/EP0) out of the sampling loop and using a        "
		// "   rejection function normalized to its maximum, the sampling     "
		// "   technique is the same as the original EGS4 implementation.     "
		// "                                                                  "
		// "   I. Kawrakow, January 2000                                      "
		// "                                                                  "
		// "******************************************************************"
		// ; Copyright NRC;

		// $COMIN-ANNIH; "DEFAULT REPLACEMENT PRODUCES THE FOLLOWING:      "
		// "COMIN/DEBUG,STACK,UPHIOT,USEFUL,RANDOM,          "
		// "EGS-VARIANCE-REDUCTION/;                         "

		// $DEFINE-LOCAL-VARIABLES-ANNIH;
		// REPLACE {$DEFINE-LOCAL-VARIABLES-ANNIH;} WITH
		// {;
		// "Local variables in order of appearance"
		double PAVIP = 0.0;// "precise total energy in the laboratory frame"
		double PESG1 = 0.0;// "precise energy of 1st annihilation photon"
		double PESG2 = 0.0;// "precise energy of 2nd annihilation photon"
		double AVIP = 0.0;// "total energy in the laboratory frame"
		double A = 0.0;// "total energy in units of the electron's rest energy"
		double G = 0.0;
		double T = 0.0;
		double P = 0.0;// "energy, kinetic energy and momentum in units of RM"
		double POT = 0.0;// "P/T"
		double EP0 = 0.0;// "minimum fractional energy"
		double WSAMP = 0.0;// "to avoid un-necessary calc. of Log((1-ep0)/ep0)"
		double RNNO01 = 0.0;// "random numbers"
		double RNNO02 = 0.0;//
		double EP = 0.0;// "fractional energy of the more energetic photon"
		double REJF = 0.0;// "rejection function"
		double ESG1 = 0.0;// "energy of the more energetic photon"
		double ESG2 = 0.0;// "energy of the less energetic photon"
		double aa = 0.0;
		double bb = 0.0;
		double cc = 0.0;
		double sinpsi = 0.0;
		double sindel = 0.0;
		double cosdel = 0.0;
		double us = 0.0;
		double vs = 0.0;
		double cphi = 0.0;
		double sphi = 0.0;
		// "for inline rotations"
		int ibr = 0;
		// $DEFINE-VARIABLES-FOR-SELECT-AZIMUTHAL-ANGLE;
		// $REAL xphi,xphi2,yphi,yphi2,rhophi2;
		double xphi = 0.0;
		double xphi2 = 0.0;
		double yphi = 0.0;
		double yphi2 = 0.0;
		double rhophi2 = 0.0;
		// }

		EGS4.NPold = EGS4.NP;// "Set the old stack counter"
		if (EGS4.nbr_split <= 0) {
			return;
		}
		PAVIP = EGS4.E[EGS4.NP - 1] + EGS4.PRM; // "PRECISE AVAILABLE ENERGY OF
												// INCIDENT POSITRON,
		// "i.e. electron assumed to be at rest
		AVIP = PAVIP; // "AVAILABLE ENERGY OF INCIDENT POSITRON"
		A = AVIP / EGS4.RM;
		// "AI=1.0/A;  AI not necessary, IK Oct 97"
		G = A - 1.0;
		T = G - 1.0;
		P = Math.sqrt(A * T);
		POT = P / T;
		EP0 = 1.0 / (A + P);
		// "   SAMPLE 1/EP FROM EP=EP0 TO 1.0-EP0"
		// "Take the calculation of the logarithm out of the loop, IK Oct 97"
		WSAMP = Math.log((1.0 - EP0) / EP0);

		aa = EGS4.U[EGS4.NP - 1];
		bb = EGS4.V[EGS4.NP - 1];
		cc = EGS4.W[EGS4.NP - 1];
		sinpsi = aa * aa + bb * bb;
		if (sinpsi > 1.e-20) {
			sinpsi = Math.sqrt(sinpsi);
			sindel = bb / sinpsi;
			cosdel = aa / sinpsi;
		}

		if (EGS4.nbr_split > 1) {
			EGS4.WT[EGS4.NP - 1] = EGS4.WT[EGS4.NP - 1] / EGS4.nbr_split;
		}

		// DO ibr = 1,nbr_split
		for (ibr = 1; ibr <= EGS4.nbr_split; ibr++) {// "nbr_split > 1 means we want splitting for any"
														// "radiative event                              "

			if (EGS4.NP + 1 > EGS4.$MXSTACK) {
				// OUTPUT np+1; (//' Stack overflow in ANNIH! np = ',i6,
				// ' Increase $MXSTACK and try again'//);
				// stop;
				EGS4.STOPPROGRAM = true;
				EGS4.seqStr = " Stack overflow in ANNIH at rest! np = "
						+ EGS4.NP + "  \n" + " Increase $MXSTACK and try again";// +"  \n";;
				// if(EGS4.iprint>2)
				eq.printSequence(EGS4.seqStr);

				return;

			}

			while (true)// LOOP[
			{
				RNNO01 = EGS4.random01();
				EP = EP0 * Math.exp(RNNO01 * WSAMP);
				// "   NOW DECIDE WHETHER TO ACCEPT"
				RNNO02 = EGS4.random01();
				// "REJF=1.0-EP+AI*AI*(2.0*G-1.0/EP);"
				// "The above rejection function has a maximum = 1 - 2/A**2"
				// "For efficiency, it is better to divide by the maximum value, IK Oct 97"
				REJF = 1.0 - (EP * A - 1.0) * (EP * A - 1.0)
						/ (EP * (A * A - 2.0));

				if (RNNO02 <= REJF)
					break;
			}// UNTIL (RNNO02 <= REJF);

			// "   SET UP ENERGIES"
			ESG1 = AVIP * EP; // "ENERGY OF SECONDARY GAMMA 1"
			PESG1 = ESG1; // "PRECISE ENERGY OF SECONDARY GAMMA 1"
			EGS4.E[EGS4.NP - 1] = PESG1;
			EGS4.IQ[EGS4.NP - 1] = 0;
			// $TRANSFER PROPERTIES TO (np) FROM (NPold);
			// @=======================================================
			int ip = 0;
			// IF( ibr = 1 ) [ ip = npold; ] ELSE [ ip = np-1; ]
			if (ibr == 1) {
				ip = EGS4.NPold;
			} else {
				ip = EGS4.NP - 1;
			}
			// $TRANSFER PROPERTIES TO (np) FROM (ip);
			// ======================================================
			EGS4.X[EGS4.NP - 1] = EGS4.X[ip - 1];
			EGS4.Y[EGS4.NP - 1] = EGS4.Y[ip - 1];
			EGS4.Z[EGS4.NP - 1] = EGS4.Z[ip - 1];
			EGS4.IR[EGS4.NP - 1] = EGS4.IR[ip - 1];
			EGS4.WT[EGS4.NP - 1] = EGS4.WT[ip - 1];
			EGS4.DNEAR[EGS4.NP - 1] = EGS4.DNEAR[ip - 1];
			EGS4.LATCH[EGS4.NP - 1] = EGS4.LATCH[ip - 1];

			EGS4.COSTHE = Math.min(1.0, (ESG1 - EGS4.RM) * POT / ESG1);
			EGS4.SINTHE = Math.sqrt(1.0 - EGS4.COSTHE * EGS4.COSTHE);
			// $SELECT-AZIMUTHAL-ANGLE(cphi,sphi);
			do {
				// $RANDOMSET xphi;
				xphi = EGS4.random01();
				xphi = 2.0 * xphi - 1.0;
				xphi2 = xphi * xphi;
				// $RANDOMSET yphi;
				yphi = EGS4.random01();
				yphi2 = yphi * yphi;
				rhophi2 = xphi2 + yphi2;
			} while (rhophi2 > 1.0);
			rhophi2 = 1.0 / rhophi2;
			cphi = (xphi2 - yphi2) * rhophi2;
			sphi = 2.0 * xphi * yphi * rhophi2;

			if (sinpsi >= 1.e-10) {
				us = EGS4.SINTHE * cphi;
				vs = EGS4.SINTHE * sphi;
				EGS4.U[EGS4.NP - 1] = cc * cosdel * us - sindel * vs + aa
						* EGS4.COSTHE;
				EGS4.V[EGS4.NP - 1] = cc * sindel * us + cosdel * vs + bb
						* EGS4.COSTHE;
				EGS4.W[EGS4.NP - 1] = cc * EGS4.COSTHE - sinpsi * us;
			} else {
				EGS4.U[EGS4.NP - 1] = EGS4.SINTHE * cphi;
				EGS4.V[EGS4.NP - 1] = EGS4.SINTHE * sphi;
				EGS4.W[EGS4.NP - 1] = EGS4.COSTHE;
			}
			EGS4.NP = EGS4.NP + 1;
			PESG2 = PAVIP - PESG1;
			ESG2 = PESG2;// esg2 = pesg2;
			EGS4.E[EGS4.NP - 1] = PESG2;// pesg2;
			EGS4.IQ[EGS4.NP - 1] = 0;
			// $TRANSFER PROPERTIES TO (np) FROM (NPold);
			// @====================================================
			// $TRANSFER PROPERTIES TO (np) FROM (np-1);
			// ======================================================
			EGS4.X[EGS4.NP - 1] = EGS4.X[EGS4.NP - 2];
			EGS4.Y[EGS4.NP - 1] = EGS4.Y[EGS4.NP - 2];
			EGS4.Z[EGS4.NP - 1] = EGS4.Z[EGS4.NP - 2];
			EGS4.IR[EGS4.NP - 1] = EGS4.IR[EGS4.NP - 2];
			EGS4.WT[EGS4.NP - 1] = EGS4.WT[EGS4.NP - 2];
			EGS4.DNEAR[EGS4.NP - 1] = EGS4.DNEAR[EGS4.NP - 2];
			EGS4.LATCH[EGS4.NP - 1] = EGS4.LATCH[EGS4.NP - 2];

			EGS4.COSTHE = Math.min(1.0, (ESG2 - EGS4.RM) * POT / ESG2);
			EGS4.SINTHE = -Math.sqrt(1.0 - EGS4.COSTHE * EGS4.COSTHE);
			if (sinpsi >= 1.e-10) {
				us = EGS4.SINTHE * cphi;
				vs = EGS4.SINTHE * sphi;
				EGS4.U[EGS4.NP - 1] = cc * cosdel * us - sindel * vs + aa
						* EGS4.COSTHE;
				EGS4.V[EGS4.NP - 1] = cc * sindel * us + cosdel * vs + bb
						* EGS4.COSTHE;
				EGS4.W[EGS4.NP - 1] = cc * EGS4.COSTHE - sinpsi * us;
			} else {
				EGS4.U[EGS4.NP - 1] = EGS4.SINTHE * cphi;
				EGS4.V[EGS4.NP - 1] = EGS4.SINTHE * sphi;
				EGS4.W[EGS4.NP - 1] = EGS4.COSTHE;
			}
			EGS4.NP = EGS4.NP + 1;

		}
		EGS4.NP = EGS4.NP - 1;

		return;
		// "END OF SUBROUTINE ANNIH" END;
	}

	// "======================================================================"
	// "                 subroutine msdist_pII                                "
	// "                 =====================                                "
	// "                                                                      "
	// "  This subroutine models multiple elastic scattering and spatial      "
	// "  deflections for a given path-length tustep.                         "
	// "  For description of input and output variables see below             "
	// "                                                                      "
	// "  September 1996      Iwan Kawrakow        Initial coding (in fortran)"
	// "  March 1997          Alex Bielajew        Adaption for EGS4          "
	// "  April/Mai 1997      Iwan Kawrakow        Debuging of the EGS4       "
	// "                                           mortran version by Bielajew"
	// "  June 1997           Iwan Kawrakow        Improved energy loss       "
	// "                                           corrections                "
	// "  June 1999           Iwan Kawrakow        spin effects, removed      "
	// "                                           $SUBSTEP-ELOSS-EVALUATION  "
	// "                                                                      "
	// "======================================================================"
	// "                                                                      "
	/**
	 * Called by ELECTR. This subroutine models multiple elastic scattering and spatial deflections for a given path-length tustep. 
	 * @param e0 e0, electron kinetic energy at the beginning of step
	 * @param eloss eloss, energy loss for this step
	 * @param tustep tustep, total pathlength of the step
	 * @param rhof rhof, density scaling template (as in EGS4)
	 * @param medium medium, medium number
	 * @param qel qel, qel =0 for e-, =1 for e+, needed for spin effects
	 * @param spin_effects spin_effects
	 * @param u0 u0, x-direction cosine before scattering
	 * @param v0 v0, y-direction cosine before scattering
	 * @param w0 w0, z-direction cosine before scattering
	 * @param x0 x0, initial x-position
	 * @param y0 y0, initial y-position
	 * @param z0 z0, initial z-position
	 */
	protected static void msdist_pII(double e0, double eloss, double tustep,
			double rhof, int medium, int qel, boolean spin_effects, double u0,
			double v0, double w0, double x0, double y0, double z0 // "Inputs
	// ,us,vs,ws,xf,yf,zf,ustep //"Outputs
	) {
		// ; Copyright NRC;

		// " Input variables
		// " ===============
		// $REAL
		// e0, "electron kinetic energy at the beginning of step
		// eloss, "energy loss for this step
		// rhof, "density scaling template (as in EGS)
		// tustep, "total pathlength of the step,
		// u0, "x-direction cosine before scattering
		// v0, "y-direction cosine before scattering
		// w0, "z-direction cosine before scattering
		// x0, "initial x-position
		// y0, "initial y-position
		// z0 "initial z-position
		// $INTEGER
		// medium,"medium number
		// qel "=0 for e-, =1 for e+, needed for spin effects

		// $LOGICAL spin_effects;

		// " Output variables
		// " ================
		// $REAL
		double us = 0.0;// "x-direction cosine after scattering
		double vs = 0.0;// "y-direction cosine after scattering
		// ws, "z-direction cosine after scattering
		double xf = 0.0;// "final x-position after transport
		double yf = 0.0;// "final y-position after transport
		double zf = 0.0;// "final z-position after transport
		double ustep = 0.0;// "straight line distance between the initial and
							// final position

		// " Local variables
		// " ===============
		double b = 0.0;// "substep transport distance,
		double blccc = 0.0;// "multiple scattering parameter
		double xcccc = 0.0;// "multiple scattering parameter
		double c = 0.0;// "substep transport distance,
		double eta = 0.0;//
		double eta1 = 0.0;// "randomization of the substep transport distances
		double chia2 = 0.0;// "screening angle, note: our chia2 is Moliere's
							// chia2/4
		double chilog = 0.0;// "log(1+1/chia2)
		double cphi0 = 0.0;// "cosine of the azimuthal angle of the initial
							// particle relative
		// "to its coordinates
		double cphi1 = 0.0;// "cosine of the first azimuthal angle
		double cphi2 = 0.0;// "cosine of the second azimuthal angle
		// /////////////////w1, "cosine of the first substep polar scattering
		// angle
		// ////////////////w2, "cosine of the second substep polar scattering
		// angle
		double w1v2 = 0.0;// "w1*v2;
		double delta = 0.0;// "transport parameter (see paper)
		double e = 0.0;// "average kinetic energy over the step
		double elke = 0.0;// "Log(e)"
		double beta2 = 0.0;// "speed at e in units of c, squared"
		double etap = 0.0;// "correction to the screening parameter derived from
							// PWA
		double xi_corr = 0.0;// "correction to the first MS moments due to spin
		double ms_corr = 0.0;//
		double tau = 0.0;// "average kinetic energy over the step divided by
							// electron mass
		double tau2 = 0.0;// "tau squared
		double epsilon = 0.0;// "fractional energy loss
		double epsilonp = 0.0;// "fractional energy loss
		double temp = 0.0;//
		double temp1 = 0.0;// "auxilarity variables for energy loss corrections
		double temp2 = 0.0;// "
		//double factor = 0.0;// "intermediate factor employed in the energy-loss
							// calculations
		double gamma = 0.0;// "q2/q1
		double lambda = 0.0;// "distance in number of elastic scattering mean
							// free paths
		// "for each sample of the multiple scattering angle
		double p2 = 0.0;// "average momentum over the step
		//double p2i = 0.0;// "inverse of ap2
		double q1 = 0.0;// "first moment of the single scattering cross section
		double rhophi2 = 0.0;// "xphi**2 + yphi**2 or its inverse
		double sint0 = 0.0;// "sine of the initial particle relative to its
							// coordinates
		double sint02 = 0.0;// "sint0**2
		double sint0i = 0.0;// "1/sint0
		// /////////////sint1, "sine of the first substep polar scattering angle
		// /////////////sint2, "sine of the second substep polar scattering
		// angle
		double sphi0 = 0.0;// "sine of the azimuthal angle of the initial
							// particle relative
		// "to its coordinates
		double sphi1 = 0.0;// "sine of the first azimuthal angle
		double sphi2 = 0.0;// "sine of the second azimuthal angle
		double u2p = 0.0;// "intermediate scatter or transport direction cosine
		double u2 = 0.0;// "sint2*cphi2;
		double v2 = 0.0;// "sint2*sphi2;
		double ut = 0.0;// "x-direction cosine for transport
		double vt = 0.0;// "y-direction cosine for transport
		double wt = 0.0;// "z-direction cosine for transport
		double xi = 0.0;// "first GS - moment
		double xphi = 0.0;// "x - used to calculated azimuthal angles
		double xphi2 = 0.0;// "xphi**2
		double yphi = 0.0;// "y - used to calculated azimuthal angles
		double yphi2 = 0.0;// "yphi**2

		// $LOGICAL
		// //////////find_index,
		// "needed to save locating the q2 index in the 2. call to mscat"
		// //spin_index
		// "saves locating the spin rejection index in 2. call to mscat"

		int lelke = 0;

		// $declare_max_medium;//REPLACE {$declare_max_medium;} WITH {;};
		// ;COMIN/ELECIN,THRESH,UPHIOT,RANDOM,CH-Steps/;

		EGS4.count_pII_steps = EGS4.count_pII_steps + 1;

		blccc = EGS4.BLCC[medium - 1];
		xcccc = EGS4.XCC[medium - 1];

		// "Commonly used factors
		e = e0 - 0.5 * eloss;
		tau = e / 0.5110034;
		tau2 = tau * tau;
		epsilon = eloss / e0;
		epsilonp = eloss / e;
		// "e = e * (1 -
		// epsilonp*epsilonp*((6+tau*(10+5*tau))/(tau+1)/(tau+2))/24);
		e = e
				* (1.0 - epsilonp * epsilonp * (6.0 + 10.0 * tau + 5.0 * tau2)
						/ (24.0 * tau2 + 48.0 * tau + 72.0));
		p2 = e * (e + EGS4.RMT2);
		// "p2i = 1/p2;
		beta2 = p2 / (p2 + EGS4.RMSQ);
		// "chia2 = xcccc*p2i/(4*blccc);
		chia2 = xcccc / (4.0 * p2 * blccc);
		lambda = 0.5 * tustep * rhof * blccc / beta2; // "The 0.5 implies a
														// half-step

		temp2 = 0.166666
				* (4.0 + tau * (6.0 + tau * (7.0 + tau * (4.0 + tau))))
				* (epsilonp / (tau + 1.0) / (tau + 2.0))
				* (epsilonp / (tau + 1.0) / (tau + 2.0));
		lambda = lambda * (1.0 - temp2);

		if (spin_effects) {
			elke = Math.log(e);
			// $SET INTERVAL elke,eke;
			Double dbl = new Double(EGS4.EKE1[medium - 1] * elke
					+ EGS4.EKE0[medium - 1]);
			lelke = dbl.intValue();

			if (lelke < 1) {// "This should normally not happen"
				lelke = 1;
				elke = (1.0 - EGS4.EKE0[medium - 1]) / EGS4.EKE1[medium - 1];
			}
			// "write(6,*) ' spin_effects: ',e,elke,lelke;
			if (qel == 0) {
				// $EVALUATE etap USING etae_ms(elke);
				etap = EGS4.etae_ms1[lelke - 1][medium - 1] * elke
						+ EGS4.etae_ms0[lelke - 1][medium - 1];

				// $EVALUATE xi_corr USING q1ce_ms(elke);
				xi_corr = EGS4.q1ce_ms1[lelke - 1][medium - 1] * elke
						+ EGS4.q1ce_ms0[lelke - 1][medium - 1];

				// $EVALUATE gamma USING q2ce_ms(elke);
				gamma = EGS4.q2ce_ms1[lelke - 1][medium - 1] * elke
						+ EGS4.q2ce_ms0[lelke - 1][medium - 1];

			} else {
				// $EVALUATE etap USING etap_ms(elke);
				etap = EGS4.etap_ms1[lelke - 1][medium - 1] * elke
						+ EGS4.etap_ms0[lelke - 1][medium - 1];

				// $EVALUATE xi_corr USING q1cp_ms(elke);
				xi_corr = EGS4.q1cp_ms1[lelke - 1][medium - 1] * elke
						+ EGS4.q1cp_ms0[lelke - 1][medium - 1];

				// $EVALUATE gamma USING q2cp_ms(elke);
				gamma = EGS4.q2cp_ms1[lelke - 1][medium - 1] * elke
						+ EGS4.q2cp_ms0[lelke - 1][medium - 1];

			}
			// $EVALUATE ms_corr USING blcce(elke);
			ms_corr = EGS4.blcce1[lelke - 1][medium - 1] * elke
					+ EGS4.blcce0[lelke - 1][medium - 1];

		} else {
			etap = 1.0;
			xi_corr = 1.0;
			gamma = 1.0;
			ms_corr = 1.0;
		}

		chia2 = chia2 * etap;
		lambda = lambda / etap / (1.0 + chia2) * ms_corr;
		chilog = Math.log(1.0 + 1.0 / chia2);
		q1 = 2.0 * chia2 * (chilog * (1.0 + chia2) - 1.0);
		gamma = 6.0 * chia2 * (1.0 + chia2)
				* (chilog * (1.0 + 2.0 * chia2) - 2.0) / q1 * gamma;
		xi = q1 * lambda;

		// "Sample first substep scattering angle
		find_index = true;
		spin_index = true;

		mscatcallindex = 1;// /added
		mscat(lambda, chia2, xi, elke, beta2, qel, medium, spin_effects);// ,find_index,spin_index,
		// w1,sint1);

		// $SELECT-AZIMUTHAL-ANGLE(cphi1,sphi1);
		do {
			// $RANDOMSET xphi;
			xphi = EGS4.random01();
			xphi = 2.0 * xphi - 1.0;
			xphi2 = xphi * xphi;
			// $RANDOMSET yphi;
			yphi = EGS4.random01();
			yphi2 = yphi * yphi;
			rhophi2 = xphi2 + yphi2;
		} while (rhophi2 > 1.0);
		rhophi2 = 1.0 / rhophi2;
		cphi1 = (xphi2 - yphi2) * rhophi2;
		sphi1 = 2.0 * xphi * yphi * rhophi2;

		// "Sample second substep scattering angle

		mscatcallindex = 2;
		mscat(lambda, chia2, xi, elke, beta2, qel, medium, spin_effects);// ,find_index,spin_index,
		// w2,sint2);

		// $SELECT-AZIMUTHAL-ANGLE(cphi2,sphi2);
		do {
			// $RANDOMSET xphi;
			xphi = EGS4.random01();
			xphi = 2.0 * xphi - 1.0;
			xphi2 = xphi * xphi;
			// $RANDOMSET yphi;
			yphi = EGS4.random01();
			yphi2 = yphi * yphi;
			rhophi2 = xphi2 + yphi2;
		} while (rhophi2 > 1.0);
		rhophi2 = 1.0 / rhophi2;
		cphi2 = (xphi2 - yphi2) * rhophi2;
		sphi2 = 2.0 * xphi * yphi * rhophi2;

		// "Final direction of motion, relative to z-axis motion
		u2 = sint2 * cphi2;
		v2 = sint2 * sphi2;
		u2p = w1 * u2 + sint1 * w2;
		us = u2p * cphi1 - v2 * sphi1;
		vs = u2p * sphi1 + v2 * cphi1;
		ws = w1 * w2 - sint1 * u2;

		// "Calculate delta, b, c

		xi = 2.0 * xi * xi_corr; // "xi was for half step, xi_corr corrects for
									// spin effects

		eta = EGS4.random01();
		eta = Math.sqrt(eta);
		eta1 = 0.5 * (1.0 - eta);
		delta = 0.9082483 - (0.1020621 - 0.0263747 * gamma) * xi;

		// "Correct the coefficients for energy loss
		temp1 = 2.0 + tau;
		temp = (2.0 + tau * temp1) / (tau + 1.0) / temp1;
		// "Take logarithmic dependence into account as well
		temp = temp - (tau + 1.0) / (tau + 2.0)
				/ (chilog * (1.0 + chia2) - 1.0);
		temp = temp * epsilonp;
		temp1 = 1.0 - temp;
		delta = delta
				+ 0.40824829
				* (epsilon * (tau + 1.0) / (tau + 2.0)
						/ (chilog * (1.0 + chia2) - 1.0)
						/ (chilog * (1.0 + 2.0 * chia2) - 2.0) - 0.25 * temp
						* temp);
		// "0.40824829 is 1/Sqrt(6)"
		b = eta * delta;
		c = eta * (1.0 - delta);

		// "Calculate transport direction cosines
		w1v2 = w1 * v2;
		ut = b * sint1 * cphi1 + c * (cphi1 * u2 - sphi1 * w1v2) + eta1 * us
				* temp1;
		vt = b * sint1 * sphi1 + c * (sphi1 * u2 + cphi1 * w1v2) + eta1 * vs
				* temp1;
		wt = eta1 * (1.0 + temp) + b * w1 + c * w2 + eta1 * ws * temp1;

		// "Calculate transport distance
		ustep = tustep * Math.sqrt(ut * ut + vt * vt + wt * wt);

		// "Rotate into the final direction of motion and transport
		// "relative to original direction of motion
		sint02 = u0 * u0 + v0 * v0;
		if (sint02 > 1.e-20) {
			sint0 = Math.sqrt(sint02);
			sint0i = 1.0 / sint0;
			cphi0 = sint0i * u0;
			sphi0 = sint0i * v0;

			// "Scattering angles
			u2p = w0 * us + sint0 * ws;
			ws = w0 * ws - sint0 * us;
			us = u2p * cphi0 - vs * sphi0;
			vs = u2p * sphi0 + vs * cphi0;

			// "Transport angles
			u2p = w0 * ut + sint0 * wt;
			wt = w0 * wt - sint0 * ut;
			ut = u2p * cphi0 - vt * sphi0;
			vt = u2p * sphi0 + vt * cphi0;
		}
		// ELSE [ wt = w0*wt; ws = w0*ws; ]
		else {
			wt = w0 * wt;
			ws = w0 * ws;
		}

		// "Transport-----------------------if msdist2call==0
		xf = x0 + tustep * ut;
		yf = y0 + tustep * vt;
		zf = z0 + tustep * wt;

		if (msdist2call == 0) {
			uscat = us;// "x-axis direction cosine for scattering"
			vscat = vs;// "y-axis direction cosine for scattering"
			wscat = ws;// "z-axis direction cosine for scattering"
			xtrans = xf;// "final x-axis position after transport"
			ytrans = yf;// "final y-axis position after transport"
			ztrans = zf;// "final z-axis position after transport"
			EGS4.USTEP = ustep;
		}

		return;
	}// end;

	// "                 subroutine msdist_pI                                 "
	// "                 ====================                                 "
	// "                                                                      "
	// "  This subroutine models multiple elastic scattering and spatial      "
	// "  deflections for a given path-length tustep                          "
	// "  resampling PRESTA-I behaviour.                                      "
	// "                                                                      "
	// "  October 1997        Iwan Kawrakow        Initial coding             "
	// "  June    1999        Iwan Kawrakow        spin effects               "
	// "                                                                      "
	// "======================================================================"
	// "                                                                      "
	/**
	 * Called by ELECTR. This subroutine models multiple elastic scattering and spatial deflections for a given path-length tustep 
	 * resampling (old) PRESTA-I behaviour.
	 * @param e0 e0, electron kinetic energy at the beginning of step
	 * @param eloss eloss, energy loss for this step
	 * @param tustep tustep, total pathlength of the step
	 * @param rhof rhof, density scaling template (as in EGS4)
	 * @param medium medium, medium number
	 * @param qel qel, qel =0 for e-, =1 for e+, needed for spin effects
	 * @param spin_effects spin_effects
	 * @param u0 u0, x-direction cosine before scattering
	 * @param v0 v0, y-direction cosine before scattering
	 * @param w0 w0, z-direction cosine before scattering
	 * @param x0 x0, initial x-position
	 * @param y0 y0, initial y-position
	 * @param z0 z0, initial z-position	
	 */
	protected static void msdist_pI(double e0, double eloss, double tustep,
			double rhof, int medium, int qel, boolean spin_effects, double u0,
			double v0, double w0, double x0, double y0, double z0 // "Inputs
	// ,us,vs,ws,xf,yf,zf,ustep // "Outputs
	) {

		// ; Copyright NRC;

		// " Input variables
		// " ===============
		// double
		// e0, "electron kinetic energy at the beginning of step
		// eloss, "energy loss for this step
		// rhof, "density scaling template (as in EGS)
		// tustep, "total pathlength of the step,
		// u0, "x-direction cosine before scattering
		// v0, "y-direction cosine before scattering
		// w0, "z-direction cosine before scattering
		// x0, "initial x-position
		// y0, "initial y-position
		// z0 "initial z-position
		// $INTEGER
		// medium,"medium number
		// qel "=0 for e-, =1 for e+, needed for spin effects
		// $LOGICAL spin_effects

		// " Output variables
		// " ================
		double us = 0.0;// "x-direction cosine after scattering
		double vs = 0.0;// "y-direction cosine after scattering
		// //////////////ws, //"z-direction cosine after scattering
		double xf = 0.0;// , "final x-position after transport
		double yf = 0.0; // "final y-position after transport
		double zf = 0.0; // "final z-position after transport
		double ustep = 0.0; // "straight line distance between the initial and
							// final position

		// " Local variables
		// " ===============
		// $REAL
		double blccc = 0.0; // "multiple scattering parameter
		double xcccc = 0.0;// "multiple scattering parameter
		double z = 0.0;
		double r = 0.0;
		double z2 = 0.0;
		double r2 = 0.0;
		// "used to calculate PLC and lateral deflection a la PRESTA-I
		double r2max = 0.0;
		double chia2 = 0.0;// "screening angle, note: our chia2 is Moliere's
							// chia2/4
		double chilog = 0.0;// "log(1+1/chia2)
		double cphi0 = 0.0;// "cosine of the azimuthal angle of the initial
							// particle relative
		// "to its coordinates
		double cphi = 0.0;// "cosine of the azimuthal scattering angle
		double sphi = 0.0;// "sine of the azimuthal scattering angle
		double e = 0.0;// "average kinetic energy over the step
		double elke = 0.0;// "Log(e)
		double beta2 = 0.0;// "speed at e in units of c, squared
		double etap = 0.0;// "correction to the screening angle derived from PWA
		double xi_corr = 0.0;// "correction to the first MS moment due to spin
		double ms_corr = 0.0;//
		double epsilon = 0.0;// "fractional energy loss
		double temp = 0.0;// "auxilarity variable for energy loss corrections
		double factor = 0.0;// "intermediate factor employed in the energy-loss
							// calculations
		double lambda = 0.0;// "distance in number of elastic scattering mean
							// free paths
		double p2 = 0.0;// "average momentum over the step
		double p2i = 0.0;// "inverse of p2
		double q1 = 0.0;// "first moment of the single scattering cross section
		double rhophi2 = 0.0;// "xphi**2 + yphi**2 or its inverse
		// ///////////////sint, "sine of the MS angle
		double sint0 = 0.0;// "sine of the initial particle relative to its
							// coordinates
		double sint02 = 0.0;// "sint0**2
		double sint0i = 0.0;// "1/sint0
		double sphi0 = 0.0;// "sine of the azimuthal angle of the initial
							// particle relative
		// "to its coordinates
		double u2p = 0.0;// "intermediate scatter or transport direction cosine
		double ut = 0.0;// "x-direction cosine for transport
		double vt = 0.0;// "y-direction cosine for transport
		double wt = 0.0;// "z-direction cosine for transport
		double xi = 0.0;// "first GS - moment
		double xphi = 0.0;// "x - used to calculated azimuthal angles
		double xphi2 = 0.0;// "xphi**2
		double yphi = 0.0;// "y - used to calculated azimuthal angles
		double yphi2 = 0.0;// "yphi**2

		// $LOGICAL
		// //find_index,
		// "needed to save locating the q2 index in the 2. call to mscat"
		// //spin_index

		int lelke = 0;

		// $declare_max_medium;//REPLACE {$declare_max_medium;} WITH {;};
		// ;COMIN/ELECIN,THRESH,UPHIOT,RANDOM/;

		blccc = EGS4.BLCC[medium - 1];
		xcccc = EGS4.XCC[medium - 1];

		e = e0 - 0.5 * eloss;
		p2 = e * (e + EGS4.RMT2);
		p2i = 1.0 / p2;
		chia2 = xcccc * p2i / (4.0 * blccc);
		beta2 = p2 / (p2 + EGS4.RMSQ);
		lambda = tustep * rhof * blccc / beta2;

		// "Account for energy loss in the MS distribution
		factor = 1.0 / (1.0 + 0.9784671 * e); // "0.9784671 = 1/(2*rm)
		epsilon = eloss / e0;
		epsilon = epsilon / (1.0 - 0.5 * epsilon);
		temp = 0.25 * (1.0 - factor * (1.0 - 0.333333 * factor)) * epsilon
				* epsilon;
		lambda = lambda * (1.0 + temp);

		if (spin_effects) {
			elke = Math.log(e);
			// $SET INTERVAL elke,eke;
			Double dbl = new Double(EGS4.EKE1[medium - 1] * elke
					+ EGS4.EKE0[medium - 1]);
			lelke = dbl.intValue();

			if (lelke < 1) {// "This should normally not happen"
				lelke = 1;
				elke = (1.0 - EGS4.EKE0[medium - 1]) / EGS4.EKE1[medium - 1];
			}
			if (qel == 0) {
				// $EVALUATE etap USING etae_ms(elke);
				etap = EGS4.etae_ms1[lelke - 1][medium - 1] * elke
						+ EGS4.etae_ms0[lelke - 1][medium - 1];

				// $EVALUATE xi_corr USING q1ce_ms(elke);
				xi_corr = EGS4.q1ce_ms1[lelke - 1][medium - 1] * elke
						+ EGS4.q1ce_ms0[lelke - 1][medium - 1];

			} else {
				// $EVALUATE etap USING etap_ms(elke);
				etap = EGS4.etap_ms1[lelke - 1][medium - 1] * elke
						+ EGS4.etap_ms0[lelke - 1][medium - 1];

				// $EVALUATE xi_corr USING q1cp_ms(elke);
				xi_corr = EGS4.q1cp_ms1[lelke - 1][medium - 1] * elke
						+ EGS4.q1cp_ms0[lelke - 1][medium - 1];

			}
			// $EVALUATE ms_corr USING blcce(elke);
			ms_corr = EGS4.blcce1[lelke - 1][medium - 1] * elke
					+ EGS4.blcce0[lelke - 1][medium - 1];

		} else {
			etap = 1.0;
			xi_corr = 1.0;
			ms_corr = 1.0;
		}

		chia2 = xcccc * p2i / (4.0 * blccc) * etap;
		lambda = lambda / etap / (1.0 + chia2) * ms_corr;
		chilog = Math.log(1.0 + 1.0 / chia2);
		q1 = 2.0 * chia2 * (chilog * (1.0 + chia2) - 1.0);
		xi = q1 * lambda;

		// "Sample multiple scattering angle
		find_index = true;
		spin_index = true;

		mscatcallindex = 3;
		mscat(lambda, chia2, xi, elke, beta2, qel, medium, spin_effects);// ,find_index,spin_index,
		// ws,sint);

		// $SELECT-AZIMUTHAL-ANGLE(cphi,sphi);
		do {
			// $RANDOMSET xphi;
			xphi = EGS4.random01();
			xphi = 2.0 * xphi - 1.0;
			xphi2 = xphi * xphi;
			// $RANDOMSET yphi;
			yphi = EGS4.random01();
			yphi2 = yphi * yphi;
			rhophi2 = xphi2 + yphi2;
		} while (rhophi2 > 1.0);
		rhophi2 = 1.0 / rhophi2;
		cphi = (xphi2 - yphi2) * rhophi2;
		sphi = 2.0 * xphi * yphi * rhophi2;

		us = sint * cphi;
		vs = sint * sphi;

		// "Correct xi used for the PLC calc. for spin effects
		xi = xi * xi_corr;

		// "Calculate PLC and lateral transport a la PRESTA-I
		// "Note that we use here the exact expression for <z>
		// "because it is much simpler and faster than the original PRESTA-I
		// formulas
		// "(which are also second order approximations)
		if (xi < 0.1) {
			z = 1.0 - xi * (0.5 - xi * (0.166666667 - 0.041666667 * xi));
		} else {
			z = (1.0 - Math.exp(-xi)) / xi;
		}
		r = 0.5 * sint;
		r2 = r * r;
		z2 = z * z;
		r2max = 1.0 - z2;
		if (r2max < r2) {
			r2 = r2max;
			r = Math.sqrt(r2);
		}

		// "Calculate final position vector
		ut = r * cphi;
		vt = r * sphi;
		wt = z;

		// "Calculate transport distance
		ustep = Math.sqrt(z2 + r2) * tustep;// ustep = Sqrt(z2 + r2)*tustep;

		// "Rotate into the final direction of motion and transport
		// "relative to original direction of motion
		sint02 = u0 * u0 + v0 * v0;
		if (sint02 > 1.e-20) {
			sint0 = Math.sqrt(sint02);
			sint0i = 1.0 / sint0;
			cphi0 = sint0i * u0;
			sphi0 = sint0i * v0;

			// "Scattering angles
			u2p = w0 * us + sint0 * ws;
			ws = w0 * ws - sint0 * us;
			us = u2p * cphi0 - vs * sphi0;
			vs = u2p * sphi0 + vs * cphi0;

			// "Transport angles
			u2p = w0 * ut + sint0 * wt;
			wt = w0 * wt - sint0 * ut;
			ut = u2p * cphi0 - vt * sphi0;
			vt = u2p * sphi0 + vt * cphi0;
		}
		// ELSE [ wt = w0*wt; ws = w0*ws; ]
		else {
			wt = w0 * wt;
			ws = w0 * ws;
		}

		// "Transport-----------------------if msdist1call==0
		xf = x0 + tustep * ut;
		yf = y0 + tustep * vt;
		zf = z0 + tustep * wt;

		if (msdist1call == 0) {
			uscat = us;// "x-axis direction cosine for scattering"
			vscat = vs;// "y-axis direction cosine for scattering"
			wscat = ws;// "z-axis direction cosine for scattering"
			xtrans = xf;// "final x-axis position after transport"
			ytrans = yf;// "final y-axis position after transport"
			ztrans = zf;// "final z-axis position after transport"
			EGS4.USTEP = ustep;
		}

		return;
	}// end;

	/**
	 * Internally used by ELECTR, msdist_pII and msdist_pI. Subroutine to sample multiple electron scattering angles from the exact 
	 * distribution resulting from elastic scattering described by the screened Rutherford cross section (spin_effects=false) or by 
	 * the screened Rutherford cross times Mott correction (spin_effects=true). 
	 * @param lambda lambda
	 * @param chia2 chia2
	 * @param q1 q1
	 * @param elke elke
	 * @param beta2 beta2
	 * @param qel qe1
	 * @param medium medium
	 * @param spin_effects spin_effects
	 */
	protected static void mscat(double lambda, double chia2, double q1,
			double elke, double beta2, int qel, int medium, boolean spin_effects)// ,
																					// boolean
																					// find_index)
	// ,boolean spin_index,double cost,double sint)
	{
		// "========================================================================="
		// "                                                                         "
		// " Subroutine to sample multiple electron scattering angles from the exact "
		// " distribution resulting from elastic scattering described by the screened"
		// " Rutherford cross section (spin_effects=.false.) or by the screened      "
		// " Rutherford cross times Mott correction (spin_effects=.true.)            "
		// "                                                                         "
		// " I.Kawrakow, NRC                                                         "
		// "========================================================================="

		// ; Copyright NRC;

		// $REAL lambda, chia2,q1,elke,beta2,cost,sint;
		// $INTEGER qel,medium;
		// $LOGICAL spin_effects,find_index,spin_index;

		// COMIN/MS-Data,RANDOM/;

		double sprob = 0.0;
		double explambda = 0.0;
		double wsum = 0.0;
		double wprob = 0.0;
		double xi = 0.0;
		double rejf = 0.0;
		//double spin_rejection = 0.0;
		double cosz = 0.0;
		double sinz = 0.0;
		double phi = 0.0;
		// double omega2=0.0;
		double llmbda = 0.0;
		double ai = 0.0;
		double aj = 0.0;
		double ak = 0.0;
		double a = 0.0;
		double u = 0.0;
		double du = 0.0;
		double x1 = 0.0;
		double rnno = 0.0;

		int icount = 0;
		// int i=0;
		// int j=0;double omega2=0.0;
		int k = 0;
		double costt = 0.0;
		double sintt = 0.0;
		// save i,j,omega2;

		if (lambda <= 13.8) {
			// "Test only for lambda = 13.8 implies a 1e-6 error, ie
			// large-lambda cases
			// "that contribute to the forward no-scattering amplitude.
			sprob = EGS4.random01();
			explambda = Math.exp(-lambda);
			if (sprob < explambda) {
				// "It was a no scattering event
				costt = 1.0;
				sintt = 0.0;
				if (mscatcallindex == 0) {
					EGS4.COSTHE = costt;
					EGS4.SINTHE = sintt;
				} else if (mscatcallindex == 1) {
					w1 = costt;
					sint1 = sintt;// 1
				} else if (mscatcallindex == 2) {
					w2 = costt;
					sint2 = sintt;// 2
				}
				if (mscatcallindex == 3) {
					ws = costt;
					sint = sintt;// 3
				}
				return;
			}
			wsum = (1.0 + lambda) * explambda;
			if (sprob < wsum) {
				// :RETRY_1:;
				while (true)// ^^^^^^^^^^^^^^^^
				{// ^^^^^^^^^^^^^^^^^^^
					xi = EGS4.random01();
					xi = 2.0 * chia2 * xi / (1.0 - xi + chia2);
					costt = 1.0 - xi;
					if (spin_effects) {
						rejf = spin_rejection(qel, medium, elke, beta2, q1,
								costt, false);
						// spin_index,false);
						rnno = EGS4.random01();
						if (rnno <= rejf)// if( rnno > rejf )
						{
							break;// GOTO :RETRY_1:;
						}
					} else
						// no spin effect
						break;
				}// ^^^^^^^^^^^
				sintt = Math.sqrt(xi * (2.0 - xi));
				if (mscatcallindex == 0) {
					EGS4.COSTHE = costt;
					EGS4.SINTHE = sintt;
				} else if (mscatcallindex == 1) {
					w1 = costt;
					sint1 = sintt;// 1
				} else if (mscatcallindex == 2) {
					w2 = costt;
					sint2 = sintt;// 2
				}
				if (mscatcallindex == 3) {
					ws = costt;
					sint = sintt;// 3
				}

				return;
			}
			if (lambda <= 1.0) // "IK introduced this portion because with
								// "alternative BCAs mscat can be called with
								// " lambda < 1 where there are no
								// pre-calculated
								// "data
			{
				wprob = explambda;
				wsum = explambda;
				costt = 1.0;
				sintt = 0.0;
				icount = 0;
				while (true)// LOOP
				{
					icount = icount + 1;
					if (icount > 20)
						break;// EXIT; //"To avoid underflow if sprob very close
								// to 1
					wprob = wprob * lambda / icount;
					wsum = wsum + wprob;
					// :RETRY_2:;
					while (true)// ^^^^^^^^^
					{// ^^^^^^^^^
						xi = EGS4.random01();
						xi = 2.0 * chia2 * xi / (1.0 - xi + chia2);
						cosz = 1.0 - xi;
						if (spin_effects) {
							rejf = spin_rejection(qel, medium, elke, beta2, q1,
									cosz, false);
							// spin_index,.false.);
							rnno = EGS4.random01();
							if (rnno <= rejf)// ( rnno > rejf )
							{
								break;
								// GOTO :RETRY_2:;
							}
						} else
							break;
					}// //^^^^^^^^^
					sinz = xi * (2.0 - xi);
					if (sinz > 1.e-20) {
						sinz = Math.sqrt(sinz);
						xi = EGS4.random01();
						phi = xi * 6.2831853;
						costt = costt * cosz - sintt * sinz * Math.cos(phi);
						sintt = Math.sqrt(Math.max(0.0, (1.0 - costt)
								* (1.0 + costt)));
					}
					// ----------------------------
					if (wsum > sprob)
						break;
				} // UNTIL ( wsum > sprob);
				if (mscatcallindex == 0) {
					EGS4.COSTHE = costt;
					EGS4.SINTHE = sintt;
				} else if (mscatcallindex == 1) {
					w1 = costt;
					sint1 = sintt;// 1
				} else if (mscatcallindex == 2) {
					w2 = costt;
					sint2 = sintt;// 2
				}
				if (mscatcallindex == 3) {
					ws = costt;
					sint = sintt;// 3
				}

				return;
			}
		}

		// "It was a multiple scattering event
		// "Sample the angle from the q^(2+) surface

		if (lambda <= EGS4.$LAMBMAX_MS) {

			if (find_index) {
				llmbda = Math.log(lambda);

				// " First fix lambda bin
				ai = llmbda * EGS4.dllambi;
				Double dbl = new Double(ai);
				i = dbl.intValue();
				ai = ai - i;
				xi = EGS4.random01();
				if (xi < ai)
					i = i + 1;

				// " fix now q1 bin
				if (q1 < EGS4.$QMIN_MS) {
					j = 0;
				} else if (q1 < EGS4.$QMAX_MS) {
					aj = q1 * EGS4.dqmsi;
					Double dbll = new Double(aj);
					j = dbll.intValue();
					aj = aj - j;
					xi = EGS4.random01();
					if (xi < aj)
						j = j + 1;
				} else {
					j = EGS4.$MAXQ_MS;
				}

				// " Calculate omega2 "
				if (llmbda < 2.2299) {
					omega2 = chia2
							* (lambda + 4.0)
							* (1.347006 + llmbda
									* (0.209364 - llmbda
											* (0.45525 - llmbda
													* (0.50142 - 0.081234 * llmbda))));
				} else {
					omega2 = chia2
							* (lambda + 4.0)
							* (-2.77164 + llmbda
									* (2.94874 - llmbda
											* (0.1535754 - llmbda * 0.00552888)));
				}

				find_index = false;
			}
			// "If this is a re-iteration with the same lambda, then omega2, i,
			// and k
			// "should have been defined in the previous iteration

			// RETRY_3:;
			while (true)// ^^^^^^^^^
			{// ^^^^^^^^^
				xi = EGS4.random01();
				ak = xi * EGS4.$MAXU_MS;
				Double dblll = new Double(ak);
				k = dblll.intValue();
				ak = ak - k;
				if (ak > EGS4.wms_array[i][j][k])// see egs is 0 biased!!!
					k = EGS4.ims_array[i][j][k];
				a = EGS4.fms_array[i][j][k];
				u = EGS4.ums_array[i][j][k];
				du = EGS4.ums_array[i][j][k + 1] - u;
				xi = EGS4.random01();
				if (Math.abs(a) < 0.2) {
					x1 = 0.5 * (1.0 - xi) * a;
					u = u + xi * du * (1.0 + x1 * (1.0 - xi * a));
				} else {
					u = u - du / a
							* (1.0 - Math.sqrt(1.0 + xi * a * (2.0 + a)));
				}

				xi = omega2 * u / (1.0 + 0.5 * omega2 - u);
				if (xi > 1.99999) {
					xi = 1.99999;
				}
				// "some machines have trouble when xi is very close to 2 in subsequent"
				// "calculations. IK, April 25 2002"
				costt = 1.0 - xi;
				if (spin_effects) {
					rejf = spin_rejection(qel, medium, elke, beta2, q1, costt,
							false);//
					// spin_index,.false.);
					rnno = EGS4.random01();
					if (rnno <= rejf)// ( rnno > rejf )
					{
						break;// GOTO :RETRY_3:;
					}
				} else
					break;
			}
			sintt = Math.sqrt(xi * (2.0 - xi));
			if (mscatcallindex == 0) {
				EGS4.COSTHE = costt;
				EGS4.SINTHE = sintt;
			} else if (mscatcallindex == 1) {
				w1 = costt;
				sint1 = sintt;// 1
			} else if (mscatcallindex == 2) {
				w2 = costt;
				sint2 = sintt;// 2
			}
			if (mscatcallindex == 3) {
				ws = costt;
				sint = sintt;// 3
			}

			return;
		}

		// "This is an error condition
		EGS4.STOPPROGRAM = true;

		EGS4.seqStr = "*************************************"
				+ "Maximum step size in mscat exceeded!"
				+ "Maximum step size initialized: 100000" + "Present lambda: "
				+ lambda + "chia2: " + chia2 + "q1, elke,  beta2: " + q1
				+ " , " + elke + " , " + beta2 + "medium: " + medium
				+ "Stopping execution "
				+ "*************************************";// +" \n";
		// if(EGS4.iprint>2)
		eq.printSequence(EGS4.seqStr);

		return;

		// end;
	}

	/**
	 * Called by ELECTR. Sample the single elastic scattering. 
	 * @param chia2 chia2
	 * @param elke elke
	 * @param beta2 beta2
	 * @param qel qe1
	 * @param medium medium
	 * @param spin_effects spin_effects
	 */
	protected static void sscat(double chia2, double elke, double beta2,
			int qel, int medium, boolean spin_effects)
	// ,double cost,double sint)
	{
		// "============================================================================"
		// "                                                                            "
		// " single elastic scattering                                                  "
		// "                                                                            "
		// " I.Kawrakow, NRC                                                            "
		// "============================================================================"
		// ; Copyright NRC;

		// $REAL chia2,elke,beta2,cost,sint;
		// $INTEGER qel,medium;
		// $LOGICAL spin_effects;

		// COMIN/RANDOM/;

		double xi = 0.0;
		double rnno = 0.0;
		double rejf = 0.0;
		//double spin_rejection = 0.0;
		double qzero = 0.0;
		// $LOGICAL spin_index;
		boolean retry_spin = false;
		spin_index = true;
		// :RETRY-SPIN:;
		while (true) {
			retry_spin = false;
			xi = EGS4.random01();
			xi = 2.0 * chia2 * xi / (1.0 - xi + chia2);
			double cost = 1.0 - xi;
			// EGS4.COSTHE=cost;
			if (spin_effects) {
				qzero = 0.0;
				// rejf =
				// spin_rejection(qel,medium,elke,beta2,qzero,cost,spin_index,.true.);
				rejf = spin_rejection(qel, medium, elke, beta2, qzero, cost,
						true);
				rnno = EGS4.random01();
				if (rnno > rejf)
					retry_spin = true;// goto :RETRY-SPIN:;
			}
			if (!retry_spin) {
				double sint = Math.sqrt(xi * (2.0 - xi));// ------>added
				if (sscatcallindex == 0) {
					EGS4.SINTHE = sint;
					EGS4.COSTHE = cost;
				}
				return;
			}
		}// while

	}

	/**
	 * Used internally in mscat and sscat. Determines the rejection function due to spin effects.
	 * @param qel charge qel (=0 for e-, =1 for e+)  
	 * @param medium medium
	 * @param elke log(energy)
	 * @param beta2 beta2
	 * @param q1 first MS moment
	 * @param cost cos(theta)
	 * @param is_single is_single
	 * @return the result
	 */ 
	protected static double spin_rejection(int qel, int medium, double elke,
			double beta2, double q1, double cost, boolean is_single)
	// boolean spin_index,boolean is_single)
	// "============================================================================="
	// "                                                                             "
	// " Determines the rejection function due to spin effects for                   "
	// "   charge        qel (=0 for e-, =1 for e+)                                  "
	// "   log(energy)   elke                                                        "
	// "   speed         beta2                                                       "
	// "   1. MS moment  q1                                                          "
	// "   cos(theta)    cost                                                        "
	// "                                                                             "
	// " I.Kawrakow, NRC                                                             "
	// "============================================================================="
	// ; Copyright NRC;
	{
		double spin_rejection = 0.0;

		// $REAL elke,beta2,q1,cost;
		// $INTEGER qel,medium;
		// $LOGICAL spin_index,is_single;

		// $declare_max_medium;//REPLACE {$declare_max_medium;} WITH {;};
		// COMIN/Spin-Data,RANDOM/;

		double rnno = 0.0;
		double ai = 0.0;
		double qq1 = 0.0;
		double aj = 0.0;
		double xi = 0.0;
		double ak = 0.0;
		// int i=0;int j=0;
		int k = 0;

		// save i,j;-->srj local saved var!!!

		if (spin_index) {
			// "Determine the energy and q1 index
			spin_index = false;
			if (beta2 >= EGS4.b2spin_min) {
				ai = (beta2 - EGS4.b2spin_min) * EGS4.dbeta2i;
				Double dbl = new Double(ai);
				isrj = dbl.intValue();
				ai = ai - isrj;
				isrj = isrj + EGS4.$MAXE_SPIN + 1;
			} else if (elke > EGS4.espml) {
				ai = (elke - EGS4.espml) * EGS4.dleneri;
				Double dbl = new Double(ai);
				isrj = dbl.intValue();
				ai = ai - isrj;
			} else {
				isrj = 0;
				ai = -1.0;
			}
			rnno = EGS4.random01();
			if (rnno < ai)
				isrj = isrj + 1;
			if (is_single) {
				jsrj = 0;
			} else {
				qq1 = 2.0 * q1;
				qq1 = qq1 / (1.0 + qq1);
				aj = qq1 * EGS4.dqq1i;
				Double dbl = new Double(aj);
				jsrj = dbl.intValue();
				if (jsrj >= EGS4.$MAXQ_SPIN) {
					jsrj = EGS4.$MAXQ_SPIN;
				} else {
					aj = aj - jsrj;
					rnno = EGS4.random01();
					if (rnno < aj)
						jsrj = jsrj + 1;
				}
			}
		}
		xi = Math.sqrt(0.5 * (1.0 - cost));
		ak = xi * EGS4.$MAXU_SPIN;
		Double dbl = new Double(ak);
		k = dbl.intValue();
		ak = ak - k;
		// spin_rej($MXMED,0:1,$0-MAXE_SPI1,$0-MAXQ_SPIN,$0-MAXU_SPIN)
		spin_rejection = (1.0 - ak)
				* EGS4.spin_rej[medium - 1][qel][isrj][jsrj][k] + // (medium,qel,i,j,k)
																	// +
				ak * EGS4.spin_rej[medium - 1][qel][isrj][jsrj][k + 1];// (medium,qel,i,j,k+1);

		return spin_rejection;
	}

	// "******************************************************************"
	// "                               National Research Council of Canada"
	/**
	 * Called by SHOWER. Handles photon interaction with matter.
	 */
	protected static void PHOTON()// (int ircode)
	{
		// "                                                                  "
		// "******************************************************************"
		// ; Copyright NRC;
		// $INTEGER IRCODE; "1 => normal return"
		// $COMIN-PHOTON; "default replacement produces the following:
		// "COMIN/DEBUG,BOUNDS,MEDIA,MISC,EPCONT,PHOTIN,STACK,THRESH,"
		// "  UPHIOT,USEFUL,USER,RANDOM,EGS-VARIANCE-REDUCTION/;"

		// $DEFINE-LOCAL-VARIABLES-PHOTON;
		double PEIG = 0.0;// "precise photon energy"
		double EIG = 0.0;// "photon energy"
		double RNNO35 = 0.0;// "random number for default MFP selection"
		double GMFPR0 = 0.0;// photon MFP before density scaling and coherent
							// correction
		double GMFP = 0.0;// "photon MFP after density scaling"@mean free path
		double COHFAC = 0.0;// "Rayleigh scattering correction"
		double RNNO37 = 0.0;// "random number for Rayleigh scattering selection"
		double XXX = 0.0; // "random number for momentum transfer sampling in Rayleigh"
		double X2 = 0.0; // "scaled momentum transfer in Rayleigh scattering event"
		double Q2 = 0.0; // "momentum transfer squared in Rayleigh scattering event"
		double CSQTHE = 0.0;// "COSTHE**2"
		double REJF = 0.0; // "Rayleigh scattering rejection function"
		double RNNORJ = 0.0;// "random number for rejection in Rayleigh scattering"
		double RNNO36 = 0.0;// "random number for interaction branching"
		double GBR1 = 0.0; // "probability for pair production"
		double GBR2 = 0.0; // "probability for pair + compton"
		//double T = 0.0; // "used for particle exchange on the stack"
		int IARG = 0;// "parameter for eq.AUSGAB"
		int IDR = 0;// "parameter for eq.AUSGAB"
		int IRL = 0;// "region number"
		int LGLE = 0;// "index for GMFP interpolation"
		int LXXX = 0;// "index for Rayleigh scattering cummulative distribution int."
		// -------------------
		boolean PCUT_DISCARD = false;
		boolean USER_PHOTON_DISCARD = false;
		boolean PNEWENERGYb = false;
		boolean SAMPLING_LOOPb = false;
		boolean PAIR_ELECTRONS_KILLED = false;
		boolean electronB = false;
		// -------------------

		IRCODE = 1;// "set up normal return"
		PEIG = EGS4.E[EGS4.NP - 1];// System.out.println(PEIG);
		EIG = PEIG; // "energy of incident gamma"
		IRL = EGS4.IR[EGS4.NP - 1];// eg: region 2 correspond to med[1]

		// $start_new_particle;
		// REPLACE {$start_new_particle;} WITH { medium = med(irl); };
		EGS4.MEDIUM = EGS4.MED[IRL - 1];// eg. med[1]=1 mediul 1

		if (EIG <= EGS4.PCUT[IRL - 1]) {
			// GO TO :PCUT-DISCARD:
			PCUT_DISCARD = true;
		} else {
			PCUT_DISCARD = false;
		}

		if (!PCUT_DISCARD) {
			PNEWENERGY: while (true) {// "enter this loop for each photon with new energy"
										// exit is by command: break PNEWENERGY;
				PAIR_ELECTRONS_KILLED = false;// @@@@@@@@@@@@@@@@@@@
				PNEWENERGYb = false;// @@@@@@@@@@@@@@@@@22
				electronB = false;
				if (EGS4.WT[EGS4.NP - 1] == 0.0) {
					// go to :USER-PHOTON-DISCARD:;
					USER_PHOTON_DISCARD = true;
					electronB = true;
					break PNEWENERGY;
				}// "added May 01"
					// @@@@@@@test
					// USER_PHOTON_DISCARD=true;
					// break PNEWENERGY;
					// @@@@@test
				EGS4.GLE = Math.log(EIG);// "GLE IS GAMMA LOG ENERGY"
				// "   here to sample no. mfp to transport before interacting"
				if (EGS4.ispmfp == 0) {
					// $SELECT-PHOTON-MFP;
					// "MACRO FOR SELECTION OF THE PHOTON MEAN-FREE-PATH"
					// REPLACE {$SELECT-PHOTON-MFP;} WITH {
					RNNO35 = EGS4.random01();
					if (RNNO35 == 0.0) {
						RNNO35 = 1.E-30;
					}
					EGS4.DPMFP = -Math.log(RNNO35);
					// }
				} else {
					// pass var===========================
					// EGS4Macro.eq =eq;
					EGS4Macro.peig = PEIG;
					EGS4Macro.eig = EIG;
					EGS4Macro.irl = IRL;
					EGS4Macro.LGLE = LGLE;
					EGS4Macro.GMFPR0 = GMFPR0;
					EGS4Macro.GMFP = GMFP;
					EGS4Macro.COHFAC = COHFAC;
					EGS4Macro.GBR1 = GBR1;
					EGS4Macro.GBR2 = GBR2;
					// end //pass var===========================
					EGS4Macro.SELECT_PHOTON_MFP();
					// output var==================================
					PEIG = EGS4Macro.peig;
					EIG = EGS4Macro.eig;
					IRL = EGS4Macro.irl;
					LGLE = EGS4Macro.LGLE;
					GMFPR0 = EGS4Macro.GMFPR0;
					GMFP = EGS4Macro.GMFP;
					COHFAC = EGS4Macro.COHFAC;
					GBR1 = EGS4Macro.GBR1;
					GBR2 = EGS4Macro.GBR2;

					if (EGS4Macro.returnB) {
						electronB = false;
						break PNEWENERGY;
					}
					if (EGS4Macro.PCUT_DISCARD) {
						PCUT_DISCARD = true;
						electronB = true;
						break PNEWENERGY;
					}

					// end output var==================================
				}
				// " DEFAULT FOR $SELECT-PHOTON-MFP; IS:  $RANDOMSET RNNO35;"
				// "                                      DPMFP=-LOG(RNNO35);"
				// "NOTE:  THIS TEMPLATE CAN ALSO BE OVER-RIDDEN BY OTHER SCHEMES,"
				// "       SUCH AS THE 'EXPONENTIAL TRANSFORM' TECHNIQUE."

				EGS4.IROLD = EGS4.IR[EGS4.NP - 1];// "INITIALIZE PREVIOUS REGION"

				PNEWMEDIUM: while (true) {// "HERE EACH TIME WE CHANGE MEDIUM DURING PHOTON TRANSPORT"
					if (EGS4.MEDIUM != 0) {
						// $SET INTERVAL GLE,GE;//"SET PWLF INTERVAL"
						// else->L{P1}={P2}1(MEDIUM)*{P1}+{P2}0(MEDIUM)
						Double dbl = new Double(EGS4.GE1[EGS4.MEDIUM - 1]
								* EGS4.GLE + EGS4.GE0[EGS4.MEDIUM - 1]);
						LGLE = dbl.intValue();
						// $EVALUATE GMFPR0 USING GMFP(GLE);
						// {P1}={P2}1(L{P3},MEDIUM)*{P3}+{P2}0(L{P3},MEDIUM)
						GMFPR0 = EGS4.GMFP1[LGLE - 1][EGS4.MEDIUM - 1]
								* EGS4.GLE
								+ EGS4.GMFP0[LGLE - 1][EGS4.MEDIUM - 1];
					}

					PTRANS: while (true) {// "PHOTON TRANSPORT LOOP"
						if (EGS4.MEDIUM == 0) {
							EGS4.TSTEP = EGS4.VACDST;
						} else {
							// $SET-RHOF; //"DENSITY RATIO SCALING TEMPLATE"
							// REPLACE {$SET-RHOF;} WITH
							// {RHOF=RHOR(IRL)/RHO(MEDIUM);} //"DEFAULT"
							EGS4.RHOF = EGS4.RHOR[IRL - 1]
									/ EGS4.RHO[EGS4.MEDIUM - 1];
							GMFP = GMFPR0 / EGS4.RHOF;
							if (EGS4.iraycorr == 0) {
								// $RAYLEIGH-CORRECTION;//
								// "A RAYLEIGH SCATTERING TEMPLATE"
								// REPLACE {$RAYLEIGH-CORRECTION;} WITH {
								if (EGS4.IRAYLR[IRL - 1] == 1) {
									// $EVALUATE COHFAC USING COHE(GLE);
									COHFAC = EGS4.COHE1[LGLE - 1][EGS4.MEDIUM - 1]
											* EGS4.GLE
											+ EGS4.COHE0[LGLE - 1][EGS4.MEDIUM - 1];
									GMFP = GMFP * COHFAC;
								}
								// }
							} else {
								// pass var====================
								EGS4Macro.LGLE = LGLE;
								EGS4Macro.COHFAC = COHFAC;
								EGS4Macro.GMFP = GMFP;
								EGS4Macro.irl = IRL;
								// =============================
								EGS4Macro.RAYLEIGH_CORRECTION();
								// get var====================
								COHFAC = EGS4Macro.COHFAC;
								GMFP = EGS4Macro.GMFP;
								// =============================

							}

							EGS4.TSTEP = GMFP * EGS4.DPMFP;
						}
						// "   SET DEFAULT VALUES FOR FLAGS SENT BACK FROM USER"
						EGS4.IRNEW = EGS4.IR[EGS4.NP - 1];// "SET DEFAULT NEW REGION NUMBER"
						EGS4.IDISC = 0;// "ASSUME PHOTON NOT DISCARDED"
						EGS4.USTEP = EGS4.TSTEP;// "TRANSFER TRANSPORT DISTANCE TO USER VARIABLE"
						EGS4.TUSTEP = EGS4.USTEP;
						// "IF (USTEP.GT.DNEAR(NP)) [;CALL eq.HOWFAR;]"
						// $CALL-eq.HOWFAR-IN-PHOTON;
						// //"The above is the default replacement"
						// REPLACE {$CALL-eq.HOWFAR-IN-PHOTON;} WITH {;
						if ((EGS4.USTEP > EGS4.DNEAR[EGS4.NP - 1])
								|| (EGS4.WT[EGS4.NP - 1] <= 0)) {
							eq.HOWFAR();
						}
						// };
						// "   NOW CHECK FOR USER DISCARD REQUEST"
						if (EGS4.IDISC > 0) {// "USER REQUESTED IMMEDIATE DISCARD"
												// GO TO :USER-PHOTON-DISCARD:;
							USER_PHOTON_DISCARD = true;
							electronB = true;
							break PNEWENERGY;
						}

						EGS4.VSTEP = EGS4.USTEP; // "SET VARIABLE FOR OUTPUT CODE"
						EGS4.TVSTEP = EGS4.VSTEP;
						EGS4.EDEP = EGS4.PZERO; // "NO ENERGY DEPOSITION ON PHOTON TRANSPORT"

						// $AUSCALL($TRANAUSB);
						// REPLACE {$AUSCALL(#);} WITH
						// {IARG={P1} ; IF (IAUSFL(IARG+1).NE.0) [CALL
						// eq.AUSGAB(IARG);]} ;
						IARG = EGS4.$TRANAUSB;
						if (EGS4.iausfl[IARG] != 0) {
							eq.AUSGAB(IARG);
						}

						// "   TRANSPORT THE PHOTON"
						EGS4.X[EGS4.NP - 1] = EGS4.X[EGS4.NP - 1]
								+ EGS4.U[EGS4.NP - 1] * EGS4.USTEP;
						EGS4.Y[EGS4.NP - 1] = EGS4.Y[EGS4.NP - 1]
								+ EGS4.V[EGS4.NP - 1] * EGS4.USTEP;
						EGS4.Z[EGS4.NP - 1] = EGS4.Z[EGS4.NP - 1]
								+ EGS4.W[EGS4.NP - 1] * EGS4.USTEP;
						EGS4.DNEAR[EGS4.NP - 1] = EGS4.DNEAR[EGS4.NP - 1]
								- EGS4.USTEP;// "DEDUCT FROM DISTANCE TO NEAREST BOUNDARY"
						if (EGS4.MEDIUM != 0) {
							EGS4.DPMFP = Math.max(0., EGS4.DPMFP - EGS4.USTEP
									/ GMFP);
						}// "DEDUCT MFP'S"
						EGS4.IROLD = EGS4.IR[EGS4.NP - 1];// "SAVE PREVIOUS REGION"

						EGS4.MEDOLD = EGS4.MEDIUM;
						if (EGS4.IRNEW != EGS4.IROLD) {// "REGION CHANGE"
														// $photon_region_change;
														// REPLACE
														// {$photon_region_change;}
														// WITH {
														// $electron_region_change;
														// }
														// REPLACE
														// {$electron_region_change;}
														// WITH {
														// ir(np) = irnew; irl =
														// irnew; medium =
														// med(irl);
							EGS4.IR[EGS4.NP - 1] = EGS4.IRNEW;
							IRL = EGS4.IRNEW;
							EGS4.MEDIUM = EGS4.MED[IRL - 1];
							// };

						}

						// "   AFTER TRANSPORT CALL TO USER"
						// $AUSCALL($TRANAUSA);
						IARG = EGS4.$TRANAUSA;
						if (EGS4.iausfl[IARG] != 0) {
							eq.AUSGAB(IARG);
						}
						// "oct 31 bug found by C Ma. PCUT discard now after eq.AUSGAB call"
						if (EIG <= EGS4.PCUT[IRL - 1]) {
							// GO TO :PCUT-DISCARD:
							PCUT_DISCARD = true;
							electronB = true;
							break PNEWENERGY;
						}

						// "   NOW CHECK FOR DEFERRED DISCARD REQUEST.  MAY HAVE BEEN SET"
						// "   BY EITHER eq.HOWFAR, OR ONE OF THE TRANSPORT eq.AUSGAB CALLS"
						if (EGS4.IDISC < 0) {
							// GO TO :USER-PHOTON-DISCARD:;
							USER_PHOTON_DISCARD = true;
							electronB = true;
							break PNEWENERGY;
						}

						if (EGS4.MEDIUM != EGS4.MEDOLD) {
							// EXIT :PTRANS:;
							break PTRANS;
						}

						if ((EGS4.MEDIUM != 0) && (EGS4.DPMFP <= $EPSGMFP)) {// "TIME FOR AN INTERACTION"
																				// EXIT
																				// :PNEWMEDIUM:;
							break PNEWMEDIUM;
						}

					}// PTRANS:-->]REPEAT ":PTRANS: LOOP"

				}// PNEWMEDIUM:-->]REPEAT ":PNEWMEDIUM: LOOP"

				// "   IT IS FINALLY TIME TO INTERACT."
				// System.out.println("Interaction: ");
				// "   THE FOLLOWING MACRO ALLOWS ONE TO INTRODUCE RAYLEIGH SCATTERING"
				// $RAYLEIGH-SCATTERING;
				// REPLACE {$RAYLEIGH-SCATTERING;} WITH {
				if (EGS4.IRAYLR[IRL - 1] == 1) {
					RNNO37 = EGS4.random01();
					if (RNNO37 <= (1.0 - COHFAC)) {
						// $AUSCALL($RAYLAUSB);
						// System.out.println("Ray");
						IARG = EGS4.$RAYLAUSB;
						if (EGS4.iausfl[IARG] != 0) {
							eq.AUSGAB(IARG);
						}
						EGS4.NPold = EGS4.NP;

						SAMPLING_LOOP: while (true)// RNNORJ > REJF)->PIRS eq
													// 2.4.1 - 2.4.5
						{
							SAMPLING_LOOPb = false;
							XXX = EGS4.random01();
							// $SET INTERVAL XXX,RCO;
							// else->L{P1}={P2}1(MEDIUM)*{P1}+{P2}0(MEDIUM)
							Double dbl1 = new Double(EGS4.RCO1[EGS4.MEDIUM - 1]
									* XXX + EGS4.RCO0[EGS4.MEDIUM - 1]);
							LXXX = dbl1.intValue();
							// $EVALUATE X2 USING RSCT(XXX);
							// //{P1}={P2}1(L{P3},MEDIUM)*{P3}+{P2}0(L{P3},MEDIUM)
							X2 = EGS4.RSCT1[LXXX - 1][EGS4.MEDIUM - 1] * XXX
									+ EGS4.RSCT0[LXXX - 1][EGS4.MEDIUM - 1];
							Q2 = X2 * EGS4.RMSQ / (20.60744 * 20.60744);
							EGS4.COSTHE = 1.
									- Q2
									/ (2. * EGS4.E[EGS4.NP - 1] * EGS4.E[EGS4.NP - 1]);// OK
							// System.out.println("cos  vechi"+EGS4.COSTHE);
							// EGS4.COSTHE=1.-2.*Q2/(EGS4.E[EGS4.NP-1]*EGS4.E[EGS4.NP-1]);//OK
							// System.out.println("cos  nou"+EGS4.COSTHE);error
							// in pirs correct here
							// q2=p2+p2-2*p*p*costhe=>q2=2p2(1-costhe)!!!;p<->k
							if (Math.abs(EGS4.COSTHE) > 1.0) {
								// GO TO :SAMPLING-LOOP:;
								SAMPLING_LOOPb = true;
							}

							if (!SAMPLING_LOOPb) {
								CSQTHE = EGS4.COSTHE * EGS4.COSTHE;
								REJF = (1.0 + CSQTHE) / 2.0;
								RNNORJ = EGS4.random01();
								// ------------------------------------
								if (RNNORJ <= REJF) {
									break SAMPLING_LOOP;
								}
							}
						}// UNTIL (RNNORJ <= REJF);
						EGS4.SINTHE = Math.sqrt(1.0 - CSQTHE);
						UPHI(2, 1);// //$SELECT-AZIMUTHAL-ANGLE and
									// OLD-PARTICLE:
						// $AUSCALL($RAYLAUSA);
						IARG = EGS4.$RAYLAUSA;
						if (EGS4.iausfl[IARG] != 0) {
							eq.AUSGAB(IARG);
						}

						// GOTO :PNEWENERGY:;
						PNEWENERGYb = true;// repeat photon transport
					}
				}// end //REPLACE {$RAYLEIGH-SCATTERING;} WITH {
					// }//REPLACE
				if (!PNEWENERGYb) {
					RNNO36 = EGS4.random01(); // "THIS RANDOM NUMBER DETERMINES WHICH INTERACTION"
					// "   GBR1=PAIR/(PAIR+COMPTON+PHOTO)=PAIR/GTOTAL"
					// $EVALUATE GBR1 USING GBR1(GLE);
					// {P1}={P2}1(L{P3},MEDIUM)*{P3}+{P2}0(L{P3},MEDIUM)
					GBR1 = EGS4.GBR11[LGLE - 1][EGS4.MEDIUM - 1] * EGS4.GLE
							+ EGS4.GBR10[LGLE - 1][EGS4.MEDIUM - 1];

					if ((RNNO36 <= GBR1) && (EGS4.E[EGS4.NP - 1] > EGS4.RMT2)) {// "IT WAS A PAIR PRODUCTION"
																				// System.out.println("Pair");
						// $AUSCALL($PAIRAUSB);
						IARG = EGS4.$PAIRAUSB;
						if (EGS4.iausfl[IARG] != 0) {
							eq.AUSGAB(IARG);
						}

						PAIR();
						// "   THE FOLLOWING MACRO ALLOWS THE USER TO CHANGE THE PARTICLE"
						// "   SELECTION SCHEME (E.G., ADDING IMPORTANCE SAMPLING (SPLITTING, "
						// "   LEADING PARTICLE SELECTION, ETC.))."
						// "   (DEFAULT MACRO IS TEMPLATE '$PARTICLE-SELECTION-PHOTON' "
						// "   WHICH IN TURN HAS THE 'NULL' REPLACEMENT ';') "
						// $PARTICLE-SELECTION-PAIR;
						// REPLACE {$PARTICLE-SELECTION-PHOTON;} WITH {;}
						// REPLACE {$PARTICLE-SELECTION-COMPT;} WITH {
						// $PARTICLE-SELECTION-PHOTON;}
						// REPLACE {$PARTICLE-SELECTION-PAIR;} WITH {
						// $PARTICLE-SELECTION-PHOTON;}
						// REPLACE {$PARTICLE-SELECTION-PHOTO;} WITH {
						// $PARTICLE-SELECTION-PHOTON;}
						// $AUSCALL($PAIRAUSA);
						IARG = EGS4.$PAIRAUSA;
						if (EGS4.iausfl[IARG] != 0) {
							eq.AUSGAB(IARG);
						}
						if (EGS4.IQ[EGS4.NP - 1] != 0) {
							// EXIT :PNEWENERGY:;
							break PNEWENERGY;// is electron so electronB is set
												// to false!!
							// terminate PHOTON
						} else {// "this may happen if pair electrons killed via Russian Roulette"
								// goto :PAIR-ELECTRONS-KILLED:;
							PAIR_ELECTRONS_KILLED = true;
						}
					}// //"IT WAS A PAIR PRODUCTION"

					if (!PAIR_ELECTRONS_KILLED)// the rest
					{
						// "GBR2=(PAIR+COMPTON)/GTOTAL"
						// $EVALUATE GBR2 USING GBR2(GLE);
						// {P1}={P2}1(L{P3},MEDIUM)*{P3}+{P2}0(L{P3},MEDIUM)
						GBR2 = EGS4.GBR21[LGLE - 1][EGS4.MEDIUM - 1] * EGS4.GLE
								+ EGS4.GBR20[LGLE - 1][EGS4.MEDIUM - 1];
						if (RNNO36 < GBR2) {// "IT WAS A COMPTON"
											// System.out.println("Compton");
							// $AUSCALL($COMPAUSB);
							IARG = EGS4.$COMPAUSB;
							if (EGS4.iausfl[IARG] != 0) {
								eq.AUSGAB(IARG);
							}

							COMPT();
							// "   THE FOLLOWING MACRO ALLOWS THE USER TO CHANGE THE PARTICLE"
							// "   SELECTION SCHEME (E.G., ADDING IMPORTANCE SAMPLING (SPLITTING, "
							// "   LEADING PARTICLE SELECTION, ETC.))."
							// "   (DEFAULT MACRO IS TEMPLATE '$PARTICLE-SELECTION-PHOTON' "
							// "   WHICH IN TURN HAS THE 'NULL' REPLACEMENT ';') "
							// $PARTICLE-SELECTION-COMPT;
							// $AUSCALL($COMPAUSA);
							IARG = EGS4.$COMPAUSA;
							if (EGS4.iausfl[IARG] != 0) {
								eq.AUSGAB(IARG);
							}

							if (EGS4.IQ[EGS4.NP - 1] != 0)// "NOT PHOTON"
							{
								// EXIT:PNEWENERGY:;
								break PNEWENERGY;// is electron so electronB is
													// set to false!!
								// terminate PHOTON
							}
						} else {// "IT WAS PHOTOELECTRIC EFFECT"
								// System.out.println("Photo");
							// $AUSCALL($PHOTOAUSB);
							IARG = EGS4.$PHOTOAUSB;
							if (EGS4.iausfl[IARG] != 0) {
								eq.AUSGAB(IARG);
							}

							PHOTO();
							// "   THE FOLLOWING MACRO ALLOWS THE USER TO CHANGE THE PARTICLE"
							// "   SELECTION SCHEME (E.G., ADDING IMPORTANCE SAMPLING (SPLITTING, "
							// "   LEADING PARTICLE SELECTION, ETC.))."
							// "   (DEFAULT MACRO IS TEMPLATE '$PARTICLE-SELECTION-PHOTON' "
							// "   WHICH IN TURN HAS THE 'NULL' REPLACEMENT ';') "
							// $PARTICLE-SELECTION-PHOTO;

							// if (EGS4.NP==0)
							// {
							// IRCODE=2;
							// return;
							// }// "FOR SPECIAL PHOTO SUBPROGRAM"
							/*
							 * IF (NP = 0 | NP < NPOLD ) [RETURN;]
							 * "The above may happen if Russian Roulette is on"
							 * "NP<NPOLD means that only electrons were created in the interaction"
							 * "and that all of them were killed. Hence, the top particle on the "
							 * "stack is from a previous interaction and may be in another region"
							 * "To avoid problems with the :PNEWENERGY: loop logic, we simply force"
							 * "a return to shower so that ELECTR or PHOTON are properly re-entered."
							 * "Changed by IK Dec. 21 2006 after D. Rogers and R. Taylor found a"
							 * "wrong dose with brems splitting and Russian Roulette on in a low "
							 * "energy calculation."
							 */
							if (EGS4.NP == 0 || EGS4.NP < EGS4.NPold) {
								return;
							}
							// "    WHERE STACK BECOMES EMPTY (I.E., "
							// "    FOR FOLLOWING FLUORESCENT PHOTONS)"
							// $AUSCALL($PHOTOAUSA);
							IARG = EGS4.$PHOTOAUSA;
							if (EGS4.iausfl[IARG] != 0) {
								eq.AUSGAB(IARG);
							}

							if (EGS4.IQ[EGS4.NP - 1] != 0) {
								// EXIT :PNEWENERGY:;
								break PNEWENERGY;// is electron so electronB is
													// set to false!!
								// //terminate PHOTON
							}
						}// "END OF PHOTO ELECTRIC BLOCK"
					}
					// :PAIR-ELECTRONS-KILLED:
					// "   IF HERE, THEN GAMMA IS LOWEST ENERGY PARTICLE."
					PEIG = EGS4.E[EGS4.NP - 1];
					EIG = PEIG;
					if (EIG < EGS4.PCUT[IRL - 1]) {
						// GO TO :PCUT-DISCARD:;
						PCUT_DISCARD = true;
						electronB = true;
						break PNEWENERGY;
					}// otherwise->repeat photon transport
				}// if(!PNEWENERGYb)---from RAyleigh sampling

			}// main while=PNEWENERGY->REPEAT ":PNEWENERGY: LOOP"

			// "   IF HERE, MEANS ELECTRON TO BE TRANSPORTED NEXT"
			if (!electronB)
				return;

		}// if(!PCUT_DISCARD)

		// "---------------------------------------------"
		// "PHOTON CUTOFF ENERGY DISCARD SECTION         "
		// "---------------------------------------------"
		// :PCUT-DISCARD:
		if (PCUT_DISCARD) {
			if (EGS4.MEDIUM > 0) {
				if (EIG > EGS4.AP[EGS4.MEDIUM - 1]) {
					IDR = EGS4.$EGSCUTAUS;
				} else {
					IDR = EGS4.$PEGSCUTAUS;
				}
			} else {
				IDR = EGS4.$EGSCUTAUS;
			}
			EGS4.EDEP = PEIG;// "GET ENERGY DEPOSITION FOR USER"
			// $PHOTON-TRACK-END;
			// REPLACE {$PHOTON-TRACK-END;} WITH {; $AUSCALL(IDR); }
			IARG = IDR;
			if (EGS4.iausfl[IARG] != 0) {
				eq.AUSGAB(IARG);
			}
			IRCODE = 2;
			EGS4.NP = EGS4.NP - 1;// dec NP!!
			return;
		}

		// "---------------------------------------------"
		// "User requested photon discard section        "
		// "---------------------------------------------"
		// :USER-PHOTON-DISCARD:
		if (USER_PHOTON_DISCARD) {
			EGS4.EDEP = PEIG;
			// $AUSCALL($USERDAUS);
			IARG = EGS4.$USERDAUS;
			if (EGS4.iausfl[IARG] != 0) {
				eq.AUSGAB(IARG);
			}
			IRCODE = 2;
			EGS4.NP = EGS4.NP - 1;// dec NP!!
			// System.out.println("exit from user ph discard");
			return;
		}

	}// "END OF SUBROUTINE PHOTON" END;

	// "******************************************************************"
	// "                               National Research Council of Canada"
	/**
	 * Called by PHOTON. For a photon energy below 2.1 MeV, the energies of the pair particles are uniformly distributed in the allowed range 
	 * via the default replacement for SELECT-LOW-ENERGY-PAIR-PRODICTION; For a photon energy between 2.1 and 50 MeV the Bethe-Heitler 
	 * cross section is employed, above 50 MeV the Coulomb-corrected Bethe-Heitler is used.
	 */
	protected static void PAIR() {
		// "                                                                  "
		// "******************************************************************"
		// "   For a photon energy below 2.1 MeV, the energies of the pair    "
		// "   particles are uniformly distributed in the allowed range via   "
		// "   the default replacement for $SELECT-LOW-ENERGY-PAIR-PRODICTION;"
		// "   If the user has a better approach, modify this macro.          "
		// "   For a photon energy between 2.1 and 50 MeV the Bethe-Heitler   "
		// "   cross section is employed, above 50 MeV the Coulomb-corrected  "
		// "   Bethe-Heitler is used.                                         "
		// "   Modified from its original version to make compatible with the "
		// "   changes made in BREMS.                                         "
		// "                                                                  "
		// "   I. Kawrakow                                                    "
		// "******************************************************************"
		// ; Copyright NRC;

		// $COMIN-PAIR; "DEFAULT REPLACEMENT PRODUCES THE FOLLOWING:
		// "COMIN/DEBUG,BREMPR,EGS-VARIANCE-REDUCTION,STACK,"
		// "THRESH,UPHIOT,USEFUL,RANDOM/;"

		// $DEFINE-LOCAL-VARIABLES-PAIR;
		// REPLACE {$DEFINE-LOCAL-VARIABLES-PAIR;} WITH
		// {;
		// "Local PAIR variables in order of their appearance"
		double PEIG = 0.0;// "precise energy of incident photon"
		double PESE1 = 0.0;// "precise energy of 1st 'electron'"
		double PESE2 = 0.0;// "precise energy of 2nd 'electron'"

		double EIG = 0.0;// "energy of incident photon"
		double ESE2 = 0.0;// "total energy of lower energy 'electron'"
		double rnno30 = 0.0;
		double rnno31 = 0.0;
		double rnno32 = 0.0;
		double rnno33 = 0.0;
		double rnno34 = 0.0;// "random numbers"

		double delta = 0.0;// "scaled momentum transfer"
		double rejf = 0.0; // "screening rejection function"
		double rejmax = 0.0; // "the maximum of rejf"
		double aux1 = 0.0;
		double aux2 = 0.0; // "auxilary variables"
		double Amax = 0.0; // "Maximum of the screening function used with (br-1/2)**2"
		double Bmax = 0.0; // "Maximum of the screening function used with the uniform part"
		double del0 = 0.0; // "delcm*eig"
		double br = 0.0; // "fraction of the available energy (eig-rmt2) going to the"
		// "lower energy `electron'"
		double Eminus = 0.0;
		double Eplus = 0.0;
		double Eavail = 0.0;
		double rnno_RR = 0.0;

		int L = 0;
		int L1 = 0; // "flags for high/low energy distributions"
		// }

		// $DEFINE-VARIABLES-FOR-SET-PAIR-ANGLE;
		// REPLACE {$DEFINE-VARIABLES-FOR-SET-PAIR-ANGLE;} WITH
		// {;
		double ESE = 0.0;// "total energy of one of the 'electrons'"
		double PSE = 0.0;// "momentum corresponding to ESE"
		double ZTARG = 0.0;// "( (1/111)*Zeff**(1/3) )**2"
		double TTEIG = 0.0;// "incident photon energy in units of RM"
		double TTESE = 0.0;// "energy of one of the 'electrons' in units of RM"
		//double TTPSE = 0.0;// "momentum of one of the 'electrons' in units of RM"
		double ESEDEI = 0.0;// "TTESE/(TTEIG-TTESE) = ratio of secondary electron energies"
		double ESEDER = 0.0;// "1/ESEDEI"
		double XIMIN = 0.0;// "1st argument where rejection function might have a maximum"
		double XIMID = 0.0;// "2nd argument where rejection function might have a maximum"
		double REJMIN = 0.0;// "rejection function at XIMIN"
		double REJMID = 0.0;// "rejection function at XIMID"
		double REJTOP = 0.0;// "max(REJMIN,REJMID)"
		double YA = 0.0;
		double XITRY = 0.0;
		double GALPHA = 0.0;
		double GBETA = 0.0;
		// "aux. variables for XIMID calculation"
		double XITST = 0.0;// "random number for pair angle sampling"
		double REJTST_on_REJTOP = 0.0;// "ratio for rejection test"
		double REJTST = 0.0;// "rejection function at XITST"
		double RTEST = 0.0;// "random number for rejection"

		int ICHRG = 0;// "loop variable"
		// }
		// 33333333333333333333333
		double k = 0.0;
		double xx = 0.0;
		double abin = 0.0;
		double rbin = 0.0;
		//double alias_sample1 = 0.0;
		int ibin = 0;
		int iq1 = 0;
		int iq2 = 0;
		int iprdst_use = 0;
		boolean do_nrc_pair = false;
		int itrip = 0;
		double ftrip = 0.0;
		// 333333333333333333333333

		EGS4.NPold = EGS4.NP; // "Set the old stack counter"

		if (EGS4.i_play_RR == 1) {// " The user wants to play Russian Roulette. For pair "
									// " it is much more efficient to do it BEFORE the "
									// " actual sampling "
			EGS4.i_survived_RR = 0; // "flag they all survive inititally"
			if (EGS4.prob_RR <= 0.0) {
				if (EGS4.n_RR_warning < EGS4.$MAX_RR_WARNING) {
					EGS4.n_RR_warning = EGS4.n_RR_warning + 1;

					EGS4.seqStr = "**** Warning, attempt to play Russian Roulette with prob_RR<0! "
							+ EGS4.prob_RR;// +" \n";
					if (EGS4.iprint > 2)
						eq.printSequence(EGS4.seqStr);

				}
			} else {
				rnno_RR = EGS4.random01();
				if (rnno_RR > EGS4.prob_RR) { // "The pair was killed "
					EGS4.i_survived_RR = 2; // "flag both particles eliminated"
					if (EGS4.NP > 1) {
						EGS4.NP = EGS4.NP - 1;// dec np photon dispatched and
												// also the pair
					} else { // " We have just one photon left on the stack. In order to  "
								// " get a proper exit from PHOTO, we have to leave at least "
								// " one particle on the stack                               "
						EGS4.WT[EGS4.NP - 1] = 0.0;
						EGS4.E[EGS4.NP - 1] = 0.0;
					}
					return;// TERMINATE PAIR!!
				} else {// "The pair survived, increase the weight"
					EGS4.WT[EGS4.NP - 1] = EGS4.WT[EGS4.NP - 1] / EGS4.prob_RR;
				}
			}
		}
		// $CHECK-STACK(np+1,'PAIR');
		// REPLACE {$CHECK-STACK(#,#);} WITH {;
		if (EGS4.NP + 1 > EGS4.$MXSTACK) {
			EGS4.STOPPROGRAM = true;

			EGS4.seqStr = " ***************************************************"
					+ "  \n"
					+ " In subroutine "
					+ "PAIR"
					+ " stack size exceeded!"
					+ "  \n"
					+ " $MXSTACK = "
					+ EGS4.$MXSTACK
					+ " np = "
					+ EGS4.NP
					+ "  \n"
					+ " Increase $MXSTACK and try again "
					+ "  \n"
					+ " Terminating execution "
					+ "  \n"
					+ " ***************************************************";// +"  \n";
			// if(EGS4.iprint>2)
			eq.printSequence(EGS4.seqStr);

			return;// stop;
		}
		// };

		PEIG = EGS4.E[EGS4.NP - 1];// "PRECISE ENERGY OF INCIDENT GAMMA"
		EIG = PEIG; // "ENERGY OF INCIDENT GAMMA"
		// ##################
		do_nrc_pair = false;

		if ((EGS4.itriplet > 0) && (EIG > 4.0 * EGS4.RM)) {
			Double tripd = new Double(EGS4.dli_triplet * EGS4.GLE
					+ EGS4.bli_triplet);
			itrip = tripd.intValue();
			ftrip = EGS4.a_triplet[itrip - 1][EGS4.MEDIUM - 1] * EGS4.GLE
					+ EGS4.b_triplet[itrip - 1][EGS4.MEDIUM - 1];// (itrip,medium);
			rnno34 = EGS4.random01();
			if (rnno34 < ftrip) { // " Triplet production "
				sample_triplet();
				return;
			}
		}

		if (EGS4.pair_nrc == 1) {// "Sample from the NRC pair cross section data base"
									// "(privided the energy is within the available range)"
			k = EIG / EGS4.RM;
			if (k < EGS4.nrcp_emax) {
				do_nrc_pair = true;
				if (k <= EGS4.nrcp_emin) {
					ibin = 1;
				} else {
					abin = 1.0 + Math.log((k - 2.0) / (EGS4.nrcp_emin - 2.0))
							* EGS4.nrcp_dlei;
					Double abind = new Double(abin);
					ibin = abind.intValue();
					abin = abin - ibin;
					rbin = EGS4.random01();
					if (rbin < abin)
						ibin = ibin + 1;
				}
				// ----------------------
				// nrcp_xdata=new double[$NRC_PAIR_NXX];65,-1=>64
				double[] patx = new double[EGS4.$NRC_PAIR_NXX];
				double[] patf = new double[EGS4.$NRC_PAIR_NXX];
				double[] patw = new double[EGS4.$NRC_PAIR_NXX];
				int[] pati = new int[EGS4.$NRC_PAIR_NXX];
				for (int ll = 0; ll < EGS4.$NRC_PAIR_NXX; ll++) {
					patx[ll] = EGS4.nrcp_xdata[ll];
					patf[ll] = EGS4.nrcp_fdata[ll][ibin - 1][EGS4.MEDIUM - 1];
				}
				for (int ll = 0; ll < EGS4.$NRC_PAIR_NXX; ll++) {
					patw[ll] = EGS4.nrcp_wdata[ll][ibin - 1][EGS4.MEDIUM - 1];
					pati[ll] = EGS4.nrcp_idata[ll][ibin - 1][EGS4.MEDIUM - 1];
				}
				// ------------------
				// xx = EGS4.alias_sample1(EGS4.$NRC_PAIR_NX_1,nrcp_xdata,
				// nrcp_fdata(1,ibin,medium),nrcp_wdata(1,ibin,medium),
				// nrcp_idata(1,ibin,medium));
				xx = EGS4.alias_sample1(EGS4.$NRC_PAIR_NX_1, patx, patf, patw,
						pati);

				// " The above returns the energy fraction of the positron "
				if (xx > 0.5) {
					PESE1 = EGS4.PRM * (1.0 + xx * (k - 2.0));
					iq1 = 1;
					PESE2 = PEIG - PESE1;
					iq2 = -1;
				} else {
					PESE2 = EGS4.PRM * (1.0 + xx * (k - 2.0));
					iq2 = 1;
					PESE1 = PEIG - PESE2;
					iq1 = -1;
				}
			}
		}
		// ##################
		if (!do_nrc_pair) {
			if (EIG <= 2.1) {
				// "   BELOW 2.1,USE APPROXIMATION"
				// $SELECT-LOW-ENERGY-PAIR-PRODICTION;
				// REPLACE {$SELECT-LOW-ENERGY-PAIR-PRODICTION;} WITH
				// {
				rnno30 = EGS4.random01();
				rnno34 = EGS4.random01();
				PESE2 = EGS4.PRM + 0.5 * rnno30 * (PEIG - 2.0 * EGS4.PRM);
				PESE1 = PEIG - PESE2;
				if (rnno34 < 0.5) {
					iq1 = -1;
					iq2 = 1;
				} else {
					iq1 = 1;
					iq2 = -1;
				}
				// see PIRS701 after eq 2.1.10
				// }
				// " IK introduced this macro because uniform energy distribution"
				// " is probably a better approximation than a zero energy 'electron'"
				// " for low energy pair production"
			} else { // "ABOVE 2.1, MUST SAMPLE"

				// " DECIDE WHETHER TO USE BETHE-HEITLER or BH
				// " COULOMB CORRECTED

				if (EIG < 50.) {// "Use BH without Coulomb correction"

					L = 5;
					L1 = L + 1;

					// "Find the actual rejection maximum for this photon energy"
					delta = 4.0 * EGS4.DELCM[EGS4.MEDIUM - 1] / EIG;// 4*delcm(medium)/eig;see
																	// PEgs and
																	// pirs eq
																	// 2.1.14
					if (delta < 1.0) {
						Amax = EGS4.DL1[L - 1][EGS4.MEDIUM - 1]
								+ delta
								* (EGS4.DL2[L - 1][EGS4.MEDIUM - 1] + delta
										* EGS4.DL3[L - 1][EGS4.MEDIUM - 1]);
						Bmax = EGS4.DL1[L1 - 1][EGS4.MEDIUM - 1]
								+ delta
								* (EGS4.DL2[L1 - 1][EGS4.MEDIUM - 1] + delta
										* EGS4.DL3[L1 - 1][EGS4.MEDIUM - 1]);
					} else {
						aux2 = Math.log(delta
								+ EGS4.DL6[L - 1][EGS4.MEDIUM - 1]);
						Amax = EGS4.DL4[L - 1][EGS4.MEDIUM - 1]
								+ EGS4.DL5[L - 1][EGS4.MEDIUM - 1] * aux2;
						Bmax = EGS4.DL4[L1 - 1][EGS4.MEDIUM - 1]
								+ EGS4.DL5[L1 - 1][EGS4.MEDIUM - 1] * aux2;
					}
					// "and then calculate the probability for sampling from (br-1/2)**2"
					aux1 = 1.0 - EGS4.RMT2 / EIG;
					aux1 = aux1 * aux1;
					aux1 = aux1 * Amax / 3.0;
					aux1 = aux1 / (Bmax + aux1);// pirs 701 eq 2.1.15
				} else {
					// "Use BH Coulomb-corrected"
					L = 7;
					// "The absolute maxima are close to the actual maxima at high energies"
					// "=>use the absolute maxima to save time"
					Amax = EGS4.DL1[L - 1][EGS4.MEDIUM - 1];
					Bmax = EGS4.DL1[L][EGS4.MEDIUM - 1];// dl1(l+1,medium);
					aux1 = EGS4.BPAR[1][EGS4.MEDIUM - 1]
							* (1.0 - EGS4.BPAR[0][EGS4.MEDIUM - 1] * EGS4.RM
									/ EIG);
				}

				del0 = EIG * EGS4.DELCM[EGS4.MEDIUM - 1];// delcm(medium);
				Eavail = EIG - EGS4.RMT2;

				// LOOP
				while (true) {

					rnno30 = EGS4.random01();
					rnno31 = EGS4.random01();
					rnno34 = EGS4.random01();
					if (rnno30 > aux1)// aux1=1-alfa (pirs701) so probab is the
										// same!!
					{ // "use the uniform part"
						br = 0.5 * rnno31;
						rejmax = Bmax;
						L1 = L + 1;// 2.1.17 pirs
					} else { // "use the (br-1/2)**2 part of the distribution"
						rnno32 = EGS4.random01();
						rnno33 = EGS4.random01();
						br = 0.5 * (1.0 - EGS4.max(rnno31, rnno32, rnno33));// //2.1.16
																			// pirs
						rejmax = Amax;
						L1 = L;
					}
					Eminus = br * Eavail + EGS4.RM;
					Eplus = EIG - Eminus;
					delta = del0 / (Eminus * Eplus);
					if (delta < 1.0) {
						rejf = EGS4.DL1[L1 - 1][EGS4.MEDIUM - 1]
								+ delta
								* (EGS4.DL2[L1 - 1][EGS4.MEDIUM - 1] + delta
										* EGS4.DL3[L1 - 1][EGS4.MEDIUM - 1]);
					} else {
						rejf = EGS4.DL4[L1 - 1][EGS4.MEDIUM - 1]
								+ EGS4.DL5[L1 - 1][EGS4.MEDIUM - 1]
								* Math.log(delta
										+ EGS4.DL6[L1 - 1][EGS4.MEDIUM - 1]);
					}

					if (rnno34 * rejmax <= rejf)
						break;
				} // UNTIL ( rnno34*rejmax <= rejf );

				PESE2 = Eminus;
				PESE1 = PEIG - PESE2;
				rnno34 = EGS4.random01();
				if (rnno34 < 0.5) {
					iq1 = -1;
					iq2 = 1;
				} else {
					iq1 = 1;
					iq2 = -1;
				}

			}
		}
		// "   ENERGY GOING TO LOWER SECONDARY HAS NOW BEEN DETERMINED"
		ESE2 = PESE2;
		EGS4.E[EGS4.NP - 1] = PESE1;
		EGS4.E[EGS4.NP] = PESE2;
		// "   THIS AVERAGE ANGLE OF EMISSION FOR BOTH PAIR PRODUCTION AND"
		// "   BREMSSTRAHLUNG IS MUCH SMALLER THAN THE AVERAGE ANGLE OF"
		// "   MULTIPLE SCATTERING FOR DELTA T TRANSPORT=0.01 R.L."
		// "   THE INITIAL AND FINAL MOMENTA ARE COPLANAR "
		// "   SET UP A NEW 'ELECTRON'  "
		// $SET-PAIR-ANGLE;===================================
		// "THE FOLLOWING REPLACES THE EGS4 DEFAULT $SET-PAIR-ANGLE MACRO    "
		// "IT'S USE REQUIRES AN ASSOCIATE MACRO $SET-PAIR-REJECTION-FUNCTION"
		// "DEFINED BELOW                                                    "
		// "                                                                 "
		// "USAGE: IPRDST=0 => EGS4 DEFAULT ANGLE SELECTION                  "
		// "       IPRDST=1 => LOWEST ORDER ANGULAR DISTRIBUTION             "
		// "                                                                 "
		// "              d(Probability)            sin(theta)               "
		// "              -------------- = -------------------------------   "
		// "                 d(theta)      2*P*[E_total - P*cos(theta)]**2   "
		// "                                                                 "
		// "       IPRDST=2 => MOTZ, OLSEN AND KOCH (1969) EQ. 3D-2003       "
		// "                   IF IPRDST IS NON-ZERO AND E_PHOTON < $BHPAIR  "
		// "                   THE IPRDST=1 DISTRIBUTION IS USED             "
		// "                                                                 "
		// REPLACE {$SET-PAIR-ANGLE;} WITH {
		// if((EGS4.iprdst==1)||((EGS4.iprdst==2)&&(EIG<$BHPAIR)))
		if (EGS4.iprdst > 0) {
			if (EGS4.iprdst == 4) {
				RTEST = EGS4.random01();
				// "gbeta = (1-rmt2/eig)**8;"
				GBETA = PESE1 / (PESE1 + 10);
				if (RTEST < GBETA) {
					iprdst_use = 1;
				} else {
					iprdst_use = 4;
				}
			} else if ((EGS4.iprdst == 2) && (EIG < $BHPAIR)) {
				iprdst_use = 1;
			} else {
				iprdst_use = EGS4.iprdst;
			}
			// DO ICHRG=1,2[
			for (ICHRG = 1; ICHRG <= 2; ICHRG++) {
				if (ICHRG == 1) {
					ESE = PESE1;
				} else {
					ESE = ESE2;
					if (EGS4.iprdst == 4) {
						GBETA = ESE / (ESE + 10.0);
						RTEST = EGS4.random01();
						if (RTEST < GBETA) {
							iprdst_use = 1;
						} else {
							iprdst_use = 4;
						}
					}
				}
				if (iprdst_use == 1) {
					PSE = Math.sqrt(Math.max(0.0, (ESE - EGS4.RM)
							* (ESE + EGS4.RM)));
					EGS4.COSTHE = EGS4.random01();
					EGS4.COSTHE = 1.0 - 2.0 * EGS4.COSTHE;
					EGS4.SINTHE = EGS4.RM
							* Math.sqrt((1.0 - EGS4.COSTHE)
									* (1.0 + EGS4.COSTHE))
							/ (PSE * EGS4.COSTHE + ESE);
					EGS4.COSTHE = (ESE * EGS4.COSTHE + PSE)
							/ (PSE * EGS4.COSTHE + ESE);
				} else if (iprdst_use == 2) {
					// "ZBRANG=( (1/111)*Zeff**(1/3) )**2"
					ZTARG = EGS4.ZBRANG[EGS4.MEDIUM - 1];
					// "TTEIG = TOTAL INITIAL PHOTON ENERGY IN ELECTRON REST MASS UNITS"
					TTEIG = EIG / EGS4.RM;
					// "TTESE = TOTAL FINAL ELECTRON ENERGY IN ELECTRON REST MASS UNITS"
					TTESE = ESE / EGS4.RM;
					// "TTPSE = TOTAL FINAL ELECTRON MOMENTUM IN ELECTRON REST MASS UNITS"
					//TTPSE = Math.sqrt((TTESE - 1.0) * (TTESE + 1.0));
					// "THIS IS THE RATIO (r IN PIRS0287)"
					ESEDEI = TTESE / (TTEIG - TTESE);
					ESEDER = 1.0 / ESEDEI;
					// "DETERMINE THE NORMALIZATION "
					XIMIN = 1.0 / (1.0 + (3.141593 * TTESE)
							* (3.141593 * TTESE));
					// $SET-PAIR-REJECTION-FUNCTION(REJMIN,XIMIN);
					REJMIN = 2.0
							+ 3.0
							* (ESEDEI + ESEDER)
							- 4.00
							* (ESEDEI + ESEDER + 1.0 - 4.0 * (XIMIN - 0.5)
									* (XIMIN - 0.5))
							* (1.0 + 0.25 * Math
									.log(((1.0 + ESEDER) * (1.0 + ESEDEI) / (2. * TTEIG))
											* ((1.0 + ESEDER) * (1.0 + ESEDEI) / (2. * TTEIG))
											+ ZTARG * XIMIN * XIMIN));
					YA = (2.0 / TTEIG) * (2.0 / TTEIG);
					XITRY = Math.max(
							0.01,
							Math.max(XIMIN,
									Math.min(0.5, Math.sqrt(YA / ZTARG))));
					GALPHA = 1.0 + 0.25 * Math.log(YA + ZTARG * XITRY * XITRY);
					GBETA = 0.5 * ZTARG * XITRY / (YA + ZTARG * XITRY * XITRY);
					GALPHA = GALPHA - GBETA * (XITRY - 0.5);
					XIMID = GALPHA / (3.0 * GBETA);
					if (GALPHA >= 0.0) {
						XIMID = 0.5 - XIMID + Math.sqrt(XIMID * XIMID + 0.25);
					} else {
						XIMID = 0.5 - XIMID - Math.sqrt(XIMID * XIMID + 0.25);
					}
					XIMID = Math.max(0.01,
							Math.max(XIMIN, Math.min(0.5, XIMID)));
					// $SET-PAIR-REJECTION-FUNCTION(REJMID,XIMID);
					REJMID = 2.0
							+ 3.0
							* (ESEDEI + ESEDER)
							- 4.00
							* (ESEDEI + ESEDER + 1.0 - 4.0 * (XIMID - 0.5)
									* (XIMID - 0.5))
							* (1.0 + 0.25 * Math
									.log(((1.0 + ESEDER) * (1.0 + ESEDEI) / (2. * TTEIG))
											* ((1.0 + ESEDER) * (1.0 + ESEDEI) / (2. * TTEIG))
											+ ZTARG * XIMID * XIMID));
					// "ESTIMATE MAXIMUM OF THE REJECTION FUNCTION"
					// "FOR LATER USE BY THE REJECTION TECHNIQUE  "
					REJTOP = 1.02 * Math.max(REJMIN, REJMID);
					while (true) {
						XITST = EGS4.random01();
						// $SET-PAIR-REJECTION-FUNCTION(REJTST,XITST);
						REJTST = 2.0
								+ 3.0
								* (ESEDEI + ESEDER)
								- 4.00
								* (ESEDEI + ESEDER + 1.0 - 4.0 * (XITST - 0.5)
										* (XITST - 0.5))
								* (1.0 + 0.25 * Math
										.log(((1.0 + ESEDER) * (1.0 + ESEDEI) / (2. * TTEIG))
												* ((1.0 + ESEDER)
														* (1.0 + ESEDEI) / (2. * TTEIG))
												+ ZTARG * XITST * XITST));

						RTEST = EGS4.random01();
						// "CONVERT THE SUCCESSFUL CANDIDATE XITST TO AN ANGLE"
						EGS4.THETA = Math.sqrt(1.0 / XITST - 1.0) / TTESE;
						// "LOOP UNTIL REJECTION TECHNIQUE ACCEPTS XITST"
						REJTST_on_REJTOP = REJTST / REJTOP;

						if ((RTEST <= REJTST_on_REJTOP)
								&& (EGS4.THETA < EGS4.PI))
							break;
					}// UNTIL((RTEST <= REJTST_on_REJTOP) & (THETA < PI) );
					EGS4.SINTHE = Math.sin(EGS4.THETA);
					EGS4.COSTHE = Math.cos(EGS4.THETA);
				} else if (iprdst_use == 3) {
					EGS4.COSTHE = EGS4.random01();
					EGS4.COSTHE = 1.0 - 2.0 * EGS4.COSTHE;
					EGS4.SINTHE = (1.0 - EGS4.COSTHE) * (1.0 + EGS4.COSTHE);
					if (EGS4.SINTHE > 0) {
						EGS4.SINTHE = Math.sqrt(EGS4.SINTHE);
					} else {
						EGS4.SINTHE = 0.0;
					}
				} else {
					// "PSE=SQRT(MAX(1e-10,(ESE-RM)*(ESE+RM)));"
					// "$RANDOMSET costhe;"
					// "costhe=(ese-(ese+pse)*exp(-2*costhe*log((ese+pse)/rm)))/pse;"
					EGS4.COSTHE = EGS4.random01();
					EGS4.COSTHE = 1.0 - 2.0 * Math.sqrt(EGS4.COSTHE);
					EGS4.SINTHE = (1 - EGS4.COSTHE) * (1 + EGS4.COSTHE);
					if (EGS4.SINTHE > 0) {
						EGS4.SINTHE = Math.sqrt(EGS4.SINTHE);
					} else {
						EGS4.SINTHE = 0.0;
					}
				}

				if (ICHRG == 1) {
					UPHI(2, 1);// $SELECT-AZIMUTHAL-ANGLE and OLD-PARTICLE:
				} else {
					EGS4.NP = EGS4.NP + 1;
					EGS4.SINTHE = -EGS4.SINTHE;
					UPHI(3, 2);// //NEW-PARTICLE
				}
			}// for
			EGS4.IQ[EGS4.NP - 1] = iq2;
			EGS4.IQ[EGS4.NP - 2] = iq1;
			return;// exit PAIR!!!!!!!!!!!!!!!!
		}// if( EGS4.iprdst > 0 )
		else {
			EGS4.THETA = 0.0;// EGS4.RM/EIG;
		}
		// }
		// ==================================================
		// " DEFAULT FOR $SET-PAIR-ANGLE; is to select the angle from the leading term"
		// " of the angular distribution "
		UPHI(1, 1);// theta + azimutal + :OLD-PARTICLE
		// "   SET UP A NEW 'ELECTRON' "
		EGS4.NP = EGS4.NP + 1;
		EGS4.SINTHE = -EGS4.SINTHE;
		UPHI(3, 2);// //NEW-PARTICLE

		EGS4.IQ[EGS4.NP - 1] = iq2;
		EGS4.IQ[EGS4.NP - 2] = iq1;
		return;
		// "END OF SUBROUTINE PAIR" END;
	}

	// "******************************************************************"
	// "                                                          NRCC    "
	/**
	 * Called by PHOTON. Subroutine for sampling incoherent (Compton) scattering. If the flag ibcmp(region) is zero, Klein-Nishina is used. 
	 * Otherwise scattering is modelled in the impuls approximation (see R.Ribberfors and K.F.Berggren, Phys.Rev.A26 (1982) 3325). 
	 * As the total cross section from PEGS4 is not modified (and thus calculated using Klein-Nishina), all rejections leed to an 
	 * unscattered photon and a zero energy electron. If a K,L1,L2,L3,M or N vacancy is created, the subsequent 
	 * atomic relaxation is treated in RELAX. This has as a consequence that more than one particle can be created as a result of an incoherent scattering. 
	 */
	protected static void COMPT() {
		// "                                VERSION 1.00  --  12 JAN 1999     "
		// "******************************************************************"
		// /"                                                                  "
		// "   Subroutine for sampling incoherent (Compton) scattering        "
		// "   If the flag ibcmp(region) is zero, Klein-Nishina is used.      "
		// "   Otherwise scattering is modelled in the impuls approximation   "
		// "    (see R.Ribberfors and K.F.Berggren, Phys.Rev.A26 (1982) 3325) "
		// "   As the total cross section from PEGS4 is not modified (and thus"
		// "   calculated using Klein-Nishina), all rejections leed to an     "
		// "   unscattered photon and a zero energy electron.                 "
		// "   If a K,L1,L2,L3,M or N vacancy is created, the subsequent      "
		// "   atomic relaxation is treated in RELAX. This has as a           "
		// "   consequence that more than one particle can be created as a    "
		// "   result of an incoherent scattering. The user should therefore  "
		// "   check their user codes for possible inconsistencies.           "
		// "                                                                  "
		// "   I.Kawrakow, January 1999                                       "
		// "******************************************************************"
		// ; Copyright NRC;

		// $COMIN-COMPT; "DEFAULT REPLACEMENT PRODUCES THE FOLLOWING:
		// "COMIN/COMPTON-DATA,DEBUG,STACK,THRESH,UPHIOT,USEFUL,RANDOM/;

		// $DEFINE-LOCAL-VARIABLES-COMPT;
		// append to $DEFINE-LOCAL-VARIABLES-COMPT;
		double acheck = 0.0;
		double acheck1 = 0.0;
		double sig_sc = 0.0;
		double sig_dc = 0.0;
		double beta_rad = 0.0;
		double alpha_rad = 0.0;
		double ux = 0.0;
		double au = 0.0;
		int icheck = 0;
		int iu = 0;

		// REPLACE {$DEFINE-LOCAL-VARIABLES-COMPT;} WITH
		// {;
		// "Local COMPT variables in order of appearance"
		double PEIG = 0.0;// "precise energy of incident photon"
		double PESG = 0.0;// "precise energy of scattered photon"
		double PESE = 0.0;// "precise total energy of compton electron"
		double ko = 0.0;// "energy of incident photon in units of RM"
		double broi = 0.0;// "1+2*ko"
		double broi2 = 0.0;// "broi*broi"
		double bro = 0.0;// "1/broi"
		double bro1 = 0.0;// "1-bro"
		double alph1 = 0.0;// "probability for the 1/BR part"
		double alph2 = 0.0;// "probability for the BR part"
		double alpha = 0.0;// "alpha1/(alph1+alph2)"
		double rnno15 = 0.0;
		double rnno16 = 0.0;
		double rnno17 = 0.0;
		double rnno18 = 0.0;
		double rnno19 = 0.0; // "random numbers"
		double br = 0.0;// "scattered photon energy fraction"
		double temp = 0.0;// "aux. variable for polar angle calculation"
		double rejf3 = 0.0;// "rejection function"
		double rejmax = 0.0;// "max. of rejf3 in thge case of uniform sampling"
		double Uj = 0.0;// "binding energy of the selected shell"
		double Jo = 0.0;// "the Compton profile parameter"
		double br2 = 0.0;// "br*br"
		double fpz = 0.0;
		double fpz1 = 0.0;// "used for limited pz-range rejection"
		double qc = 0.0;// "momentum transfer corresponding to the Compton line energy"
		double qc2 = 0.0;// "qc squared"
		double af = 0.0;// "for calculating F"
		double Fmax = 0.0;// "maximum of F"
		double frej = 0.0;// "used for F-rejection"
		double eta_incoh = 0.0;
		double eta = 0.0;// "random numbers"
		double aux = 0.0;
		double aux1 = 0.0;
		double aux2 = 0.0;
		double aux3 = 0.0;
		double aux4 = 0.0;// "aux. variables"
		double pzmax = 0.0;// "max. possible z-component of the initial electron momentum"
		double pz = 0.0;// "initial electron momentum projection"
		double pz2 = 0.0;// "pz*pz"
		double rnno_RR = 0.0;// "for playing Russian Roulette"

		int irl = 0;// "local region number"
		int i = 0;// "loop variable for shell sampling (and then shell sampled)"
		int j = 0;// "pointer to the shell in the shell data list"
		int iarg = 0;// "argument for eq.AUSGAB call"
		int ip = 0;// "a loop variable"
		// }

		boolean finishedComptonSampling = false;
		boolean RETRY_PZ = false;
		boolean RESAMPLEb = false;

		EGS4.NPold = EGS4.NP;// "Set the old stack counter"
		PEIG = EGS4.E[EGS4.NP - 1];// "PRECISE ENERGY OF INCIDENT GAMMA"
		ko = PEIG / EGS4.RM;// "Gamma energy in units of electron rest energy"
		broi = 1.0 + 2.0 * ko; // "Needed for scattering angle sampling"
		// ####################################
		// $RADC_CHECK;
		// *****************************************************************************
		// Macro to decide if single or double Compton event occures.
		// ****************************************************************************/
		// REPLACE {$RADC_CHECK;} WITH {;
		if ((EGS4.radc_flag == 1)
				&& ((ko > EGS4.radc_emin) && (ko < EGS4.radc_emax))) {
			acheck = EGS4.radc_dlei * EGS4.GLE + EGS4.radc_le1;
			Double ddh = new Double(acheck);
			icheck = ddh.intValue();
			acheck = acheck - icheck;
			acheck1 = 1.0 - acheck;
			// //radc_sigs=new double[$RADC_NE+1];//($0-RADC_NE)
			sig_sc = EGS4.radc_sigs[icheck] * acheck1
					+ EGS4.radc_sigs[icheck + 1] * acheck;
			// radc_sigd=new double[$RADC_NE+1];//($0-RADC_NE),
			sig_dc = EGS4.radc_sigd[icheck] * acheck1
					+ EGS4.radc_sigd[icheck + 1] * acheck;
			rnno15 = EGS4.random01();
			// / *
			// IF( rnno15 > sig_sc + sig_dc ) return;
			// "i.e. reject the interaction to take into account the cross "
			// "section reduction compared to Klein-Nishina"
			// * /
			rnno16 = EGS4.random01();
			if (rnno16 < acheck)
				icheck = icheck + 1;
			if (rnno15 * (sig_sc + sig_dc) < sig_dc) {
				sample_double_compton(ko, icheck);
				return;
			}
			// "If here, we are doing single Compton with radiative corrections"
			// "Calculate a few useful constants."
			beta_rad = ko * (1.0 + ko);
			alpha_rad = 1.0 / Math.log(1.0 + 2.0 * beta_rad);
		}
		// };
		// ####################################
		irl = EGS4.IR[EGS4.NP - 1];
		if (EGS4.ibcmp[irl - 1] == 1) {// "User wants to take into account binding effects"
										// "=>first sample the shell and see whether an    "
										// "  interaction is possible                      "
			rnno17 = EGS4.random01();
			// DO i=1,n_shell(medium) [
			for (i = 1; i <= EGS4.n_shell[EGS4.MEDIUM - 1]; i++) {
				rnno17 = rnno17 - EGS4.eno_array[i - 1][EGS4.MEDIUM - 1];// (i,medium);
				if (rnno17 <= 0)
					break;// EXIT;
			}
			// "j is the shell number in the data list"
			j = EGS4.shell_array[i - 1][EGS4.MEDIUM - 1];// (i,medium);
			// "Uj is the binding energy in units of rm"
			Uj = EGS4.be_array[j - 1];// (j);

			// " Binding energy rejection "
			if (ko <= Uj) {
				// goto :INTERACTION-REJECTED:;
				// :INTERACTION-REJECTED:
				// " Create here a zero energy electron if required (check user codes) "
				return;
			}
			// "Jo is the Compton profile parameter"
			Jo = EGS4.Jo_array[j - 1];// (j);

		}

		// " We always sample the scattering angle from Klein-Nishina"
		while (true)// added
		{// added
			RESAMPLEb = false;
			// :RESAMPLE:
			if (ko > 2.0) {// "At high energies the original EGS4 method is most efficient"
				broi2 = broi * broi;
				alph1 = Math.log(broi);
				alph2 = ko * (broi + 1.0) / broi2;
				alpha = alph1 / (alph1 + alph2);
				// LOOP [
				while (true) {
					rnno15 = EGS4.random01();
					rnno16 = EGS4.random01();
					if (rnno15 < alpha) {// "Use 1/br part"
						br = Math.exp(alph1 * rnno16) / broi;
					} else {// "Use the br part."
						br = Math.sqrt(rnno16 + (1.0 - rnno16) / broi2);
					}
					temp = (1.0 - br) / ko / br;
					EGS4.SINTHE = Math.max(0., temp * (2.0 - temp));
					rejf3 = 1.0 - br * EGS4.SINTHE / (1.0 + br * br);
					rnno19 = EGS4.random01();
					// ------------------------------------------
					if (rnno19 <= rejf3) {
						break;
					}
					// ---------------------------------------------
				}// UNTIL rnno19.le.rejf3;
			} else {// "At low energies it is faster to sample br uniformely"
				bro = 1. / broi;
				bro1 = 1.0 - bro;
				rejmax = broi + bro;
				while (true) {
					rnno15 = EGS4.random01();
					rnno16 = EGS4.random01();
					br = bro + bro1 * rnno15;
					temp = (1.0 - br) / ko / br;
					EGS4.SINTHE = Math.max(0., temp * (2.0 - temp));
					rejf3 = (br + 1. / br - EGS4.SINTHE) / rejmax;
					// ------------------------------------------
					if (rnno16 <= rejf3) {
						break;
					}
					// --------------------------------------
				}// UNTIL rnno16.le.rejf3;
			}

			if ((br < 1. / broi) || (br > 1.0))
			// if((br >= 1./broi) && (br <= 1.0))
			{
				// if( (br < 0.99999/broi) || (br > 1.00001 ))
				// {
				// write(6,*) ' sampled br outside of allowed range!
				// ',ko,1./broi,br;
				// }
				// goto :RESAMPLE: ;
				RESAMPLEb = true;
				// break;
			}
			if (!RESAMPLEb) {
				// $RADC_REJECTION
				// *****************************************************************************
				// Rejection loop macro when sampling single Compton events
				// ****************************************************************************/
				// REPLACE {$RADC_REJECTION;} WITH {;
				if ((EGS4.radc_flag == 1)
						&& ((ko > EGS4.radc_emin) && (ko < EGS4.radc_emax))) {
					ux = Math.log(1.0 + beta_rad * temp) * alpha_rad;
					au = ux * EGS4.$RADC_NU;
					Double aud = new Double(au);
					iu = aud.intValue();// au;
					au = au - iu;// radc_frej=new
									// double[$RADC_NE+1][$RADC_NU+1];
					rejf3 = EGS4.radc_frej[icheck][iu] * (1.0 - au)
							+ EGS4.radc_frej[icheck][iu + 1] * au;
					rnno15 = EGS4.random01();
					if (rnno15 > rejf3) {
						// goto :RESAMPLE:;
						RESAMPLEb = true;
					}
				}
				// };

			}

			if (!RESAMPLEb)
				break;

		}// main loop added
		EGS4.COSTHE = 1.0 - temp;// costhe = 1 - temp;

		if (EGS4.ibcmp[irl - 1] == 0) {// "User wants to use Klein-Nishina, so we are done"
			Uj = 0.0;
			// goto :FINISHED-COMPTON-SAMPLING:;
			finishedComptonSampling = true;
		}

		if (!finishedComptonSampling) {
			// " Check for rejection due to the limited range of pzmax "
			br2 = br * br;
			aux = ko * (ko - Uj) * temp;
			pzmax = (aux - Uj) / Math.sqrt(2.0 * aux + Uj * Uj);// pirs 2.2.17
			if (pzmax <= -1.0) {
				// goto :INTERACTION-REJECTED:;
				// :INTERACTION-REJECTED:
				// " Create here a zero energy electron if required (check user codes) "
				return;
			}
			qc2 = 1.0 + br * br - 2.0 * br * EGS4.COSTHE;
			qc = Math.sqrt(qc2);// ~pirs 2.2.4

			if (pzmax > 1.0) {
				pzmax = 1.0;
				af = 0.0;
				Fmax = 1.0;
				fpz = 1.0;
				// goto :RETRY-PZ:;
				RETRY_PZ = true;
			}
			if (!RETRY_PZ) {
				aux3 = 1.0 + 2.0 * Jo * Math.abs(pzmax);
				aux4 = 0.5 * (1 - aux3 * aux3);// =-b
				fpz = 0.5 * Math.exp(aux4);// si
				af = qc * (1.0 + br * (br - EGS4.COSTHE) / qc2);// 2.2.21;22,23

				if (af < 0.0) {
					if (pzmax > 0.0)
						fpz = 1.0 - fpz;
					eta_incoh = EGS4.random01();
					if (eta_incoh > fpz) {
						// goto :INTERACTION-REJECTED:;
						return;
					}
					af = 0.0;
					Fmax = 1.0;
					// goto :RETRY-PZ:;
					RETRY_PZ = true;
				}
			}// if(!RETRY_PZ)
			if (!RETRY_PZ) {
				if (pzmax < -0.15) {
					Fmax = 1.0 - af * 0.15;// pirs 2.2.20
					fpz1 = fpz * Fmax;
				} else if (pzmax < 0.15) {
					Fmax = 1.0 + af * pzmax;
					aux3 = 1.0 / (1 + 0.33267252734 * aux3);
					// "0.33267252734 is p/sqrt(2), p is the parameter from Eq. 7.1.25"
					// "of Abramowitz and Stegun, needed for approximating Erf        "
					aux4 = fpz
							* aux3
							* (0.3480242 + aux3
									* (-0.0958798 + aux3 * 0.7478556))
							+ EGS4.erfJo_array[j - 1];// (j);
					if (pzmax > 0.0) {
						// "fpz1 = 1 - Fmax*fpz - 0.31332853433*af/Jo_array(j)*aux4;"
						// "missing factor 1/2 in the above found by Cerneliu Costescu"
						// "0.62665706866 is sqrt(Pi/8)"
						fpz1 = 1.0 - Fmax * fpz - 0.62665706866 * af
								/ EGS4.Jo_array[j - 1] * aux4;
						fpz = 1.0 - fpz;
					} else {
						// "fpz1 = Fmax*fpz - 0.31332853433*af/Jo_array(j)*aux4;"
						// "missing factor 1/2 in the above found by Cerneliu Costescu"
						// "0.62665706866 is sqrt(Pi/8)"
						fpz1 = Fmax * fpz - 0.62665706866 * af
								/ EGS4.Jo_array[j - 1] * aux4;
					}
				} else {
					Fmax = 1.0 + af * 0.15;
					fpz1 = 1.0 - Fmax * fpz;
					fpz = 1.0 - fpz;
				}
				eta_incoh = EGS4.random01();
				if (eta_incoh > fpz1)
					return;// goto :INTERACTION-REJECTED:;
			}// if(!RETRY_PZ)
				// "At this point, all rejections are handled, now we need to sample pz "
				// "between -1 and pzmax using the Compton profile of the selected shell"
				// "and F(pz,cos(theta)) as a rejection function                        "

			// :RETRY-PZ:;
			// RETRY_PZ=false;
			while (true) {
				RETRY_PZ = false;
				rnno18 = EGS4.random01();
				rnno18 = rnno18 * fpz;
				if (rnno18 < 0.5) {
					rnno18 = Math.max(1.e-30, 2.0 * rnno18);
					pz = 0.5 * (1.0 - Math.sqrt(1.0 - 2.0 * Math.log(rnno18)))
							/ Jo;
				} else {
					rnno18 = 2.0 * (1.0 - rnno18);
					pz = 0.5 * (Math.sqrt(1.0 - 2.0 * Math.log(rnno18)) - 1.0)
							/ Jo;
				}
				if (Math.abs(pz) > 1.0) {
					// goto :RETRY-PZ:;
					RETRY_PZ = true;
				}
				// "Due to the non-relativistic approximation"
				// "for pz, it has to be between -1 and 1    "
				if (!RETRY_PZ) {
					if (pz < 0.15) {
						if (pz < -0.15) {
							frej = (1.0 - af * 0.15) / Fmax;
						} else {
							frej = (1.0 + af * pz) / Fmax;
						}
						eta = EGS4.random01();
						if (eta > frej) {
							// goto :RETRY-PZ:;
							RETRY_PZ = true;
						}
					}
				}// if(!RETRY_PZ)
				if (!RETRY_PZ)
					break;
			}// while true
				// "If pz > 0.15, F is always 1 => no need for rejection"
				// " Calculate energy of scattered photon "
			pz2 = pz * pz;
			if (Math.abs(pz) < 0.01) {
				br = br * (1.0 + pz * (qc + (br2 - EGS4.COSTHE) * pz));
			} else {
				aux = 1.0 - pz2 * br * EGS4.COSTHE;
				aux1 = 1.0 - pz2 * br2;
				aux2 = qc2 - br2 * pz2 * EGS4.SINTHE;
				if (aux2 > 1.e-10) {
					br = br / aux1 * (aux + pz * Math.sqrt(aux2));
				}
			}
			Uj = Uj * EGS4.PRM;
		}// :FINISHED-COMPTON-SAMPLING:
			// :FINISHED-COMPTON-SAMPLING:
		PESG = br * PEIG;
		PESE = PEIG - PESG - Uj + EGS4.PRM;
		EGS4.SINTHE = Math.sqrt(EGS4.SINTHE);
		UPHI(2, 1);// //$SELECT-AZIMUTHAL-ANGLE and OLD-PARTICLE:
		EGS4.E[EGS4.NP - 1] = PESG;
		aux = 1.0 + br * br - 2.0 * br * EGS4.COSTHE;
		if (aux > 1.e-8) {
			EGS4.COSTHE = (1.0 - br * EGS4.COSTHE) / Math.sqrt(aux);
			EGS4.SINTHE = (1.0 - EGS4.COSTHE) * (1.0 + EGS4.COSTHE);
			if (EGS4.SINTHE > 0.0) {
				EGS4.SINTHE = -Math.sqrt(EGS4.SINTHE);
			} else {
				EGS4.SINTHE = 0.0;
			}
		} else {
			EGS4.COSTHE = 0.0;
			EGS4.SINTHE = -1.0;
		}
		EGS4.NP = EGS4.NP + 1;// inc NP
		// $CHECK-STACK(np,'COMPT');
		// REPLACE {$CHECK-STACK(#,#);} WITH {;
		if (EGS4.NP > EGS4.$MXSTACK) {
			EGS4.STOPPROGRAM = true;

			EGS4.seqStr = " ***************************************************"
					+ "  \n"
					+ " In subroutine "
					+ "COMPT"
					+ " stack size exceeded!"
					+ "  \n"
					+ " $MXSTACK = "
					+ EGS4.$MXSTACK
					+ " np = "
					+ EGS4.NP
					+ "  \n"
					+ " Increase $MXSTACK and try again "
					+ "  \n"
					+ " Terminating execution "
					+ "  \n"
					+ " ***************************************************";// +"  \n";
			// if(EGS4.iprint>2)
			eq.printSequence(EGS4.seqStr);

			return;// stop;
		}
		// };

		UPHI(3, 2);// NEW-PARTICLE
		EGS4.E[EGS4.NP - 1] = PESE;
		EGS4.IQ[EGS4.NP - 1] = -1;

		if (EGS4.ibcmp[irl - 1] == 1) {
			// " Shell vacancy "
			if (Uj > 1.e-3) {
				EGS4.EDEP = 0.0;
				// relax(Uj,shn_array(j),iz_array(j));
				RELAX(Uj, EGS4.shn_array[j - 1], EGS4.iz_array[j - 1]);
				// "relax will put all particles with energies above ecut,pcut on the "
				// "stack, the remaining energy will be scored in edep and deposited  "
				// "localy (via the call to eq.AUSGAB below)                             "
			} else {
				EGS4.EDEP = Uj;
			}
			// $AUSCALL($PHOTXAUS); "generates IARG = 4 call"
			if (EGS4.EDEP > 0.0) {
				iarg = EGS4.$PHOTXAUS;
				if (EGS4.iausfl[iarg] != 0) {
					eq.AUSGAB(iarg);
				}
			}

		}

		// " Now play Russian Roulette with resulting electrons if the user asked for it"
		// $PLAY RUSSIAN ROULETTE WITH ELECTRONS FROM NPold+1; //"TO NP"
		// "This macro implements Russian Roulette (most useful  with brems splitting)"
		// "It is more efficient than having the user do it via eq.AUSGAB since it avoids"
		// "considerable handling of the particles by ELECTR"
		// "The user must set i_play_RR (defaults to 0) and prob_RR"
		// "Both are in COMIN EGS-VARIANCE-REDUCTION"
		// ""
		// "Note that this macro is called as $PLAY RUSSIAN ROULETTE WITH ELECTRONS..."
		// "Note also that subroutine pair has its own, internal version"

		// REPLACE {$PLAYRUSSIANROULETTEWITHELECTRONSFROM#;} WITH {;

		EGS4.i_survived_RR = 0; // "flag all survive"
		if (EGS4.i_play_RR == 1) {
			if (EGS4.prob_RR <= 0.0) {
				if (EGS4.n_RR_warning < EGS4.$MAX_RR_WARNING) {
					EGS4.n_RR_warning = EGS4.n_RR_warning + 1;
					// OUTPUT prob_RR;
					EGS4.seqStr = "**** Warning, attempt to play Russian Roulette with prob_RR<0! "
							+ EGS4.prob_RR;// +" \n";
					if (EGS4.iprint > 2)
						eq.printSequence(EGS4.seqStr);

				}
			} else {
				ip = EGS4.NPold + 1;
				while (true) {// "handle all particles from p1 to np"
					if (EGS4.IQ[ip - 1] != 0) {// "i.e. charged particles"
						rnno_RR = EGS4.random01();
						if (rnno_RR < EGS4.prob_RR) {// "particle survives"
							EGS4.WT[ip - 1] = EGS4.WT[ip - 1] / EGS4.prob_RR;
							ip = ip + 1; // "increase local pointer"
						} else {// "particle killed"
							EGS4.i_survived_RR = EGS4.i_survived_RR + 1;
							if (ip < EGS4.NP) {// "=>replace it with last particle"
												// "on stack"
								EGS4.E[ip - 1] = EGS4.E[EGS4.NP - 1];
								EGS4.IQ[ip - 1] = EGS4.IQ[EGS4.NP - 1];
								EGS4.WT[ip - 1] = EGS4.WT[EGS4.NP - 1];
								EGS4.U[ip - 1] = EGS4.U[EGS4.NP - 1];
								EGS4.V[ip - 1] = EGS4.V[EGS4.NP - 1];
								EGS4.W[ip - 1] = EGS4.W[EGS4.NP - 1];
							}
							EGS4.NP = EGS4.NP - 1; // "reduce stack by one=> particle gone"
						}// "end of kill particle block"
					} else {// "this is a photon, leave it. Change pointer"
						ip = ip + 1;
					}

					if (ip > EGS4.NP)
						break;
				} // UNTIL (ip > np);
					// "loops until either np is decreased to ip, or ip increased to np"
				if (EGS4.NP == 0) {// " we need at least one particle on the stack "
									// " so that the transport routines can exit properly"
					EGS4.NP = 1;
					EGS4.E[EGS4.NP - 1] = 0.0;
					EGS4.IQ[EGS4.NP - 1] = 0;
					EGS4.WT[EGS4.NP - 1] = 0.0;
				}
			} // "end of russian roulette block"
		} // "end of flag set block"
		// };

		return;

		// :INTERACTION-REJECTED:
		// " Create here a zero energy electron if required (check user codes) "
		// return;
	}

	// "******************************************************************"
	/**
	 * Called by PHOTON. Sample the photoelectric effect. The emission of Fluorescent X-rays, Auger, 
	 * Coster-Kronig electrons are treated in RELAX.
	 */
	protected static void PHOTO() {
		// "******************************************************************"
		// " Programmers:  I. Kawrakow, complete recoding,                    "
		// "                            Fluorescent X-rays, Auger,            "
		// "                            Coster-Kronig treated in RELAX        "
		// "               A.F. Bielajew (NRC) photoelectric angular distn    "
		// "******************************************************************"

		// ; Copyright NRC;

		// $COMIN-PHOTO; "default replacement is:
		// "COMIN/BOUNDS,DEBUG,EDGE,EGS-VARIANCE-REDUCTION,EPCONT,"
		// "MEDIA,PHOTIN,RANDOM,STACK,UPHIOT,USEFUL/"

		// $DEFINE-VARIABLES-FOR-SELECT-PHOTOELECTRON-DIRECTION;
		// REPLACE {$DEFINE-VARIABLES-FOR-SELECT-PHOTOELECTRON-DIRECTION;} WITH
		// {;
		// "Photo-electron angle selection variables"
		double EELEC = 0.0;// "total energy of photo-electron"
		double BETA = 0.0;// "velocity of electron in units of c"
		double GAMMA = 0.0;// "total energy of photo-electron in units of RM"
		double ALPHA = 0.0;// "kinematic factor"
		double RATIO = 0.0;// "=BETA/ALPHA"
		double RNPHT = 0.0;// "random number"
		double FKAPPA = 0.0;// "aux. variable for COSTHE calculation"
		double XI = 0.0;// "used in rejection function calculation"
		double SINTH2 = 0.0;// "SINTHE**2"
		double RNPHT2 = 0.0;// "random number for rejection"
		// }

		// $DEFINE-LOCAL-VARIABLES-PHOTO;
		// REPLACE {$DEFINE-LOCAL-VARIABLES-PHOTO;} WITH
		// {;
		// "Local PHOTO variables in order of their appearance"

		double PEIG = 0.0;// "precise energy of incident photon"
		double BR = 0.0;// "random number"
		double sigma = 0.0;// "elemental cross section"
		double aux = 0.0;
		double aux1 = 0.0;// "aux. variables"
		double[] probs = new double[EGS4.$MXEL];// "probability for an interaction with a given element"
		double sigtot = 0.0;// "total cross section"
		double e_vac = 0.0;// "shell binding energy"
		double rnno_RR = 0.0;// "for playing Russian Roulette"
		int IARG = 0;// "eq.AUSGAB calling switch"
		int iZ = 0;// "Atomic number of the element the photon is "
					// "interactiong with"
		int irl = 0;// "local region number"
		int[] ints = new int[EGS4.$MXEL];// "energy interval number for a given element"
		int j = 0;
		int ip = 0;// "loop variables"
		int k = 0;// "shell number"

		boolean do_relax = false;
		// save n_warning;
		// }

		// int n_warning=0;//"a warning counter"

		EGS4.NPold = EGS4.NP; // "Set the old stack counter"
		PEIG = EGS4.E[EGS4.NP - 1];
		irl = EGS4.IR[EGS4.NP - 1];
		if (PEIG < EGS4.edge_energies[1][0]) {// edge_energies(2,1) )[

			if (n_warning < 100) {
				n_warning = n_warning + 1;
				EGS4.seqStr = " Subroutine PHOTO called with E = " + PEIG
						+ " which is below the current min. energy of 1 keV! "
						+ " \n"
						+ " Converting now this photon to an electron, "
						+ " but you should check your code! ";// +" \n";
				if (EGS4.iprint > 2)
					eq.printSequence(EGS4.seqStr);

			}
			EGS4.IQ[EGS4.NP - 1] = -1;
			EGS4.E[EGS4.NP - 1] = PEIG + EGS4.PRM;
			return;
		}

		iZ = EGS4.iedgfl[irl - 1];
		do_relax = false;
		EGS4.EDEP = EGS4.PZERO;// pzero;
		if (EGS4.iedgfl[irl - 1] != 0) {// " User requested atomic relaxations "
										// " first sample the element "
			if (EGS4.NNE[EGS4.MEDIUM - 1] == 1) {
				Double dbl = new Double(EGS4.ZELEM[EGS4.MEDIUM - 1][0] + 0.5);
				iZ = dbl.intValue();// int( zelem(medium,1) + 0.5 );
				for (j = 1; j <= EGS4.edge_number[iZ - 1]; j++) {
					if (PEIG >= EGS4.edge_energies[j - 1][iZ - 1])
						break;// EXIT;
				}

			} else {
				aux = PEIG * PEIG;
				aux1 = aux * PEIG;
				aux = aux * Math.sqrt(PEIG);
				sigtot = 0.0;

				// "write(1,*) nne(medium),' elements ';
				for (k = 1; k <= EGS4.NNE[EGS4.MEDIUM - 1]; k++) {
					Double dbl = new Double(
							EGS4.ZELEM[EGS4.MEDIUM - 1][k - 1] + 0.5);
					iZ = dbl.intValue();// int( zelem(medium,k) + 0.5 );
					if ((iZ < 1) || (iZ > EGS4.$MXELEMENT)) {
						EGS4.STOPPROGRAM = true;

						EGS4.seqStr = " Error in PHOTO: " + " \n"
								+ "   Atomic number of element " + k
								+ " in medium " + EGS4.MEDIUM
								+ " is not between 1 and " + EGS4.$MXELEMENT;// +" \n";
						// if(EGS4.iprint>2)
						eq.printSequence(EGS4.seqStr);

						// write(6,*) ' Error in PHOTO: ';
						// write(6,*) ' Atomic number of element ',k,
						// ' in medium ',medium,' is not between 1 and ',
						// $MXELEMENT;
						return;
						// stop;
					}
					if (PEIG > EGS4.edge_energies[0][iZ - 1])// (1,iZ) )
					{
						j = 1;
						sigma = (EGS4.edge_a[0][iZ - 1]
								+ EGS4.edge_b[0][iZ - 1] / PEIG
								+ EGS4.edge_c[0][iZ - 1] / aux + EGS4.edge_d[0][iZ - 1]
								/ aux1)
								/ PEIG;
					} else {
						for (j = 2; j <= EGS4.edge_number[iZ - 1]; j++) {
							if (PEIG >= EGS4.edge_energies[j - 1][iZ - 1])// (j,iZ)
																			// )
								break;
						}
						sigma = EGS4.edge_a[j - 1][iZ - 1]
								+ EGS4.GLE
								* (EGS4.edge_b[j - 1][iZ - 1] + EGS4.GLE
										* (EGS4.edge_c[j - 1][iZ - 1] + EGS4.GLE
												* EGS4.edge_d[j - 1][iZ - 1]));
						sigma = Math.exp(sigma);
					}// see 2.3.2 PIRS

					sigma = sigma * EGS4.PZ[EGS4.MEDIUM - 1][k - 1];// pz(medium,k);
					sigtot = sigtot + sigma;
					probs[k - 1] = sigma;
					ints[k - 1] = j;
				}
				BR = EGS4.random01();

				for (k = 1; k <= EGS4.NNE[EGS4.MEDIUM - 1]; k++) {
					sigma = probs[k - 1] / sigtot;// see 2.3.1 PIRS
					BR = BR - sigma;
					if (BR <= 0.0)
						break;
				}
				Double dbb = new Double(
						EGS4.ZELEM[EGS4.MEDIUM - 1][k - 1] + 0.5);
				iZ = dbb.intValue();// int( zelem(medium,k) + 0.5 );
				j = ints[k - 1];

			}
			// " Now we know the atomic number (iZ) and the energy interval the "
			// " photon energy is in (j). It is time to sample the shell the photon "
			// " is interacting with. "
			// " left for now as before, to be changed!!! "
			if (PEIG <= EGS4.binding_energies[EGS4.$MXSHELL - 1][iZ - 1])// ($MXSHELL,iZ)
																			// )
			{ // "Below  N-shell -> local energy deposition "
				EGS4.EDEP = PEIG;
				EGS4.E[EGS4.NP - 1] = EGS4.PZERO;
				EGS4.WT[EGS4.NP - 1] = 0.0;
			} else {// "Above  N-shell -> sample the shell the photon is interacting with"
				BR = EGS4.random01();
				for (k = 1; k <= EGS4.$MXINTER; k++) {
					if (PEIG > EGS4.binding_energies[k - 1][iZ - 1])// (k,iZ) )
					{
						if (BR < EGS4.interaction_prob[k - 1][iZ - 1])// (k,iZ)
																		// )
							break;
						BR = (BR - EGS4.interaction_prob[k - 1][iZ - 1])
								/ (1.0 - EGS4.interaction_prob[k - 1][iZ - 1]);
					}
				}
				e_vac = EGS4.binding_energies[k - 1][iZ - 1];
				EGS4.E[EGS4.NP - 1] = PEIG - e_vac + EGS4.PRM;
				do_relax = true;
				// "IF( e_vac > 1e-3 ) [ e(np) = peig - e_vac + prm; do_relax =
				// .true.; ]
				// "ELSE [ e(np) = peig + prm; e_vac = 0; ]
				EGS4.IQ[EGS4.NP - 1] = -1;
			}
		} else {
			EGS4.E[EGS4.NP - 1] = PEIG + EGS4.PRM;
			EGS4.IQ[EGS4.NP - 1] = -1;// electron, e-
		}

		if (EGS4.IQ[EGS4.NP - 1] == -1) {
			// $SELECT-PHOTOELECTRON-DIRECTION;
			// //"Samples photo-electron direction"
			// "--------------------------------------------------------------"
			// "                                                              "
			// "           PHOTOELECTRON ANGLE SELECTION                      "
			// "           =============================                      "
			// "                                                              "
			// "--------------------------------------------------------------"
			// "This macro can be used to select the photoelectron direction  "

			// REPLACE {$SELECT-PHOTOELECTRON-DIRECTION;} WITH {
			// "        ================================"
			if (EGS4.iphter[EGS4.IR[EGS4.NP - 1] - 1] == 1) {
				EELEC = EGS4.E[EGS4.NP - 1];
				if (EELEC > EGS4.ECUT[EGS4.IR[EGS4.NP - 1] - 1]) {
					BETA = Math.sqrt((EELEC - EGS4.RM) * (EELEC + EGS4.RM))
							/ EELEC;
					GAMMA = EELEC / EGS4.RM;
					ALPHA = 0.5 * GAMMA - 0.5 + 1. / GAMMA;
					RATIO = BETA / ALPHA;
					while (true) {
						RNPHT = EGS4.random01();
						RNPHT = 2. * RNPHT - 1.;
						if (RATIO <= 0.2) {
							FKAPPA = RNPHT + 0.5 * RATIO * (1. - RNPHT)
									* (1. + RNPHT);
							if (GAMMA < 100.0) {
								EGS4.COSTHE = (BETA + FKAPPA)
										/ (1. + BETA * FKAPPA);
							} else {
								if (FKAPPA > 0.0) {
									EGS4.COSTHE = 1.0
											- (1.0 - FKAPPA)
											* (GAMMA - 3.0)
											/ (2.0 * (1.0 + FKAPPA)
													* (GAMMA - 1.0)
													* (GAMMA - 1.0) * (GAMMA - 1.0));
								} else {
									EGS4.COSTHE = (BETA + FKAPPA)
											/ (1. + BETA * FKAPPA);
								}
							}
							// "XI=1./(1.-BETA*COSTHE); <-- this numerically problematic "
							// "                            at high energies, IK"
							XI = (1.0 + BETA * FKAPPA) * GAMMA * GAMMA;
						} else {
							XI = GAMMA
									* GAMMA
									* (1. + ALPHA
											* (Math.sqrt(1. + RATIO
													* (2. * RNPHT + RATIO)) - 1.));
							EGS4.COSTHE = (1. - 1. / XI) / BETA;
						}
						SINTH2 = Math.max(0., (1. - EGS4.COSTHE)
								* (1. + EGS4.COSTHE));
						RNPHT2 = EGS4.random01();

						if (RNPHT2 <= 0.5 * (1. + GAMMA) * SINTH2 * XI / GAMMA)
							break;
					}// WHILE(RNPHT2.GT.0.5*(1.+GAMMA)*SINTH2*XI/GAMMA);
					EGS4.SINTHE = Math.sqrt(SINTH2);
					UPHI(2, 1);
				}// if(EELEC>EGS4.ECUT[EGS4.IR[EGS4.NP-1]-1])
			}// if(EGS4.IPHTER[EGS4.IR[EGS4.NP-1]-1]==1)
				// }

		}

		if (do_relax) {
			// "write(1,*) ' Calling relax, e_vac = ',e_vac;
			RELAX(e_vac, k, iZ);
		}

		if (EGS4.EDEP > 0.0) {
			// $AUSCALL($PHOTXAUS);
			IARG = EGS4.$PHOTXAUS;
			if (EGS4.iausfl[IARG] != 0) {
				eq.AUSGAB(IARG);
			}

		}// "generates IARG = 4 call"

		// $PLAY RUSSIAN ROULETTE WITH ELECTRONS FROM NPold; " TO NP;"
		EGS4.i_survived_RR = 0; // "flag all survive"
		if (EGS4.i_play_RR == 1) {
			if (EGS4.prob_RR <= 0.0) {
				if (EGS4.n_RR_warning < EGS4.$MAX_RR_WARNING) {
					EGS4.n_RR_warning = EGS4.n_RR_warning + 1;
					// OUTPUT prob_RR;
					EGS4.seqStr = "**** Warning, attempt to play Russian Roulette with prob_RR<0! "
							+ EGS4.prob_RR;// +" \n";
					if (EGS4.iprint > 2)
						eq.printSequence(EGS4.seqStr);

				}
			} else {
				ip = EGS4.NPold;
				while (true) {// "handle all particles from p1 to np"
					if (EGS4.IQ[ip - 1] != 0) {// "i.e. charged particles"
						rnno_RR = EGS4.random01();
						if (rnno_RR < EGS4.prob_RR) {// "particle survives"
							EGS4.WT[ip - 1] = EGS4.WT[ip - 1] / EGS4.prob_RR;
							ip = ip + 1; // "increase local pointer"
						} else {// "particle killed"
							EGS4.i_survived_RR = EGS4.i_survived_RR + 1;
							if (ip < EGS4.NP) {// "=>replace it with last particle"
												// "on stack"
								EGS4.E[ip - 1] = EGS4.E[EGS4.NP - 1];
								EGS4.IQ[ip - 1] = EGS4.IQ[EGS4.NP - 1];
								EGS4.WT[ip - 1] = EGS4.WT[EGS4.NP - 1];
								EGS4.U[ip - 1] = EGS4.U[EGS4.NP - 1];
								EGS4.V[ip - 1] = EGS4.V[EGS4.NP - 1];
								EGS4.W[ip - 1] = EGS4.W[EGS4.NP - 1];
							}
							EGS4.NP = EGS4.NP - 1; // "reduce stack by one=> particle gone"
						}// "end of kill particle block"
					} else {// "this is a photon, leave it. Change pointer"
						ip = ip + 1;
					}

					if (ip > EGS4.NP)
						break;
				} // UNTIL (ip > np);
					// "loops until either np is decreased to ip, or ip increased to np"
				if (EGS4.NP == 0) {// " we need at least one particle on the stack "
									// " so that the transport routines can exit properly"
					EGS4.NP = 1;
					EGS4.E[EGS4.NP - 1] = 0.0;
					EGS4.IQ[EGS4.NP - 1] = 0;
					EGS4.WT[EGS4.NP - 1] = 0.0;
				}
			} // "end of russian roulette block"
		} // "end of flag set block"
		// };

		return;

	}

	// "******************************************************************"
	/**
	 * Called by eii_sample, COMPT, PHOTO. Sample atomic relaxation. Subroutine to fill a vacancy in shell n, element iZ 
	 * by emitting fluorescent X-rays, Auger and Coster-Kronig electrons. Transitions between K,L1,L2,L3,average M,average N are taken into 
	 * account. Particles with energies above the transport cut-offs (ECUT and PCUT) are placed on the stack, energy of sub-threshold 
	 * particles is stored in EDEP. In this version a global cut-off of 1 keV applies, i.e. if ECUT-RM or PCUT is below 1 keV, binding energies below 
	 * 1 keV will still be absorbed locally (due to lack of data).
	 * @param energy energy
	 * @param n n
	 * @param iZ iZ 
	 */
	protected static void RELAX(double energy, int n, int iZ) {
		// "******************************************************************"
		// " Subroutine to fill a vacancy in shell n, element iZ              "
		// " by emitting fluorescent X-rays, Auger and Coster-Kronig electrons"
		// " Transitions between K,L1,L2,L3,average M,average N are taken into"
		// " account. Particles with energies above the transport cut-offs    "
		// " (ECUT and PCUT) are placed on the stack, energy of sub-threshold "
		// " particles is stored in EDEP.                                     "
		// " In this version a global cut-off of 1 keV applies                "
		// "  i.e. if ECUT-RM or PCUT is below 1 keV, binding energies below  "
		// "  1 keV will still be absorbed locally (due to lack of data)      "
		// "                                                                  "
		// " Version 1:  I. Kawrakow, December 1998                           "
		// "******************************************************************"

		// REPLACE {$RELAX-CUTOFF} WITH {0.001}
		$RELAX_CUTOFF = 0.001;
		// ; Copyright NRC;

		// " Input variables "
		// "================="
		// $INTEGER n,iZ;
		// $REAL energy; "

		// " Local variables "
		// "================="
		// REPLACE {$MXVAC} WITH {50}
		// "Maximum number of vacancies during the"
		// "relaxation cascade                    "
		$MXVAC = 50;

		int[] vac_array = new int[$MXVAC];// "array with shell vacancies            "
		int n_vac = 0;// "current number of vacancies           "
		int shell = 0;// "current shell                         "
		int finalL = 0;
		int finala = 0;// "code of final state                   "
		int final1 = 0;
		int final2 = 0;// "two vacancies in the final state in   "
						// "the case of Auger transitions         "
		int iql = 0;// "particle charge                       "
		int irl = 0;// "present region"
		// int[] first_transition={1,20,27,33,38};//($MXINTER),//$MXINTER=5
		// int[] last_transition={19,26,32,37,39};//($MXINTER);
		// "first and last transition for a given "
		// "shell in the list of all possible     "
		// "transitions                           "
		// int[] final_state={
		// 4,3,5,6, //" K-shell fluorescence    "
		// 202,302,402,404,403,303, //" K-shell Auger           "
		// 502,503,504,602,603,604, //" K-shell Auger           "
		// 505,605,606, //" K-shell Auger           "
		// 13,14, //" L1 Coster-Kronig        "
		// 5,6, //" L1 fluorescence         "
		// 505,605,606, //" L1 Auger                "
		// 14, //" L2 Coster-Kronig        "
		// 5,6, //" L2 fluorescence         "
		// 505,605,606, //" L2 Auger                "
		// 5,6, //" L3 fluorescence         "
		// 505,605,606, //" L3 Auger                "
		// 6, //" M  fluorescence         "
		// 606};//($MXTRANS);//=39
		// " final_state(i) is the final atomic state                "
		// " after transition i coded as follows:                    "
		// "   * fluorescence - final_state is the shell number      "
		// "                    of the new vacancy                   "
		// "   * Coster-Kronig - final_state is the shell number     "
		// "                     of the new vacancy + 10             "
		// "   * Auger - final_state is n1 + 100*n2 where n1 and n2  "
		// "             are the shell numbers of the 2 new vacancies"

		int k = 0;
		//int np_old = 0;
		//int ip = 0;
		int iarg = 0;
		double[] e_array = new double[$MXVAC];// "array with vacancy energies           "
		double Ei = 0.0;
		double Ef = 0.0;// "initial,final binding energies        "
		double Ex = 0.0;// "kinetic energy of emitted particle    "
		double eta = 0.0;// "a random number                       "
		//double e_check = 0.0;// "energy conservation check             "
		double e_cut = 0.0;
		double ekcut = 0.0;
		double pkcut = 0.0;
		double elcut = 0.0; // "cut-off energies                "

		double xphi = 0.0;
		double yphi = 0.0;
		double xphi2 = 0.0;
		double yphi2 = 0.0;
		double rhophi2 = 0.0;
		double cphi = 0.0;
		double sphi = 0.0; // "for azimuthal angle selection"

		boolean repeatLOOP = false;
		// " Global EGS4 variables "
		// "======================="
		// $COMIN-RELAX;
		// ;COMIN/RELAX-USER/;

		if ((n < 1) || (n > EGS4.$MXSHELL)) {
			return;
		}// "unknown vacancy"

		EGS4.iZ_relax = iZ;
		irl = EGS4.IR[EGS4.NP - 1];
		ekcut = EGS4.ECUT[irl - 1] - EGS4.RM;
		pkcut = EGS4.PCUT[irl - 1];
		e_cut = Math.min(ekcut, pkcut);
		e_cut = Math.max($RELAX_CUTOFF, e_cut);

		if (energy <= e_cut) {
			EGS4.EDEP = EGS4.EDEP + energy; // "We assume that edep is zeroed "
			// "(or set to the appropriate value in the routine "
			// "calling RELAX "
			return;
		}

		// " Set-up the array of vacancies for the relaxation cascade "
		n_vac = 1;
		vac_array[n_vac - 1] = n;
		//np_old = EGS4.NP;
		//e_check = 0.0;
		e_array[n_vac - 1] = energy;

		// :START:
		// LOOP[ "Until no  >N-shell vacancies"
		while (true) {// "Until no  >N-shell vacancies"
			repeatLOOP = false;// reset
			shell = vac_array[n_vac - 1];
			Ei = e_array[n_vac - 1];
			n_vac = n_vac - 1;
			if (Ei <= e_cut) {// " Below cut-off -> local absorption "
				EGS4.EDEP = EGS4.EDEP + Ei;
				if (n_vac > 0)
					repeatLOOP = true;// goto :START: ;
				else
					// added
					break;
			}
			// "Set the relax_user common block variables, IK March 22 2004"
			if (!repeatLOOP) {
				EGS4.ish_relax = shell;
				EGS4.u_relax = Ei;
				if (shell == EGS4.$MXSHELL) {// "This is N-shell vacancy -> just produce Auger"
					if (Ei > e_cut) {
						EGS4.NP = EGS4.NP + 1;
						// $CHECK-STACK(np,'RELAX');
						if (EGS4.NP > EGS4.$MXSTACK) {
							EGS4.STOPPROGRAM = true;

							EGS4.seqStr = " ***************************************************"
									+ "  \n"
									+ " In subroutine "
									+ "RELAX"
									+ " stack size exceeded!"
									+ "  \n"
									+ " $MXSTACK = "
									+ EGS4.$MXSTACK
									+ " np = "
									+ EGS4.NP
									+ "  \n"
									+ " Increase $MXSTACK and try again "
									+ "  \n"
									+ " Terminating execution "
									+ "  \n"
									+ " ***************************************************";// +"  \n";
							// if(EGS4.iprint>2)
							eq.printSequence(EGS4.seqStr);

							return;// stop;
						}

						EGS4.E[EGS4.NP - 1] = Ei + EGS4.PRM;
						EGS4.IQ[EGS4.NP - 1] = -1;
						// $TRANSFER PROPERTIES TO (np) FROM (np_old);
						// $TRANSFER PROPERTIES TO (np) FROM (np-1);
						EGS4.X[EGS4.NP - 1] = EGS4.X[EGS4.NP - 2];
						EGS4.Y[EGS4.NP - 1] = EGS4.Y[EGS4.NP - 2];
						EGS4.Z[EGS4.NP - 1] = EGS4.Z[EGS4.NP - 2];
						EGS4.IR[EGS4.NP - 1] = EGS4.IR[EGS4.NP - 2];
						EGS4.WT[EGS4.NP - 1] = EGS4.WT[EGS4.NP - 2];
						EGS4.DNEAR[EGS4.NP - 1] = EGS4.DNEAR[EGS4.NP - 2];
						EGS4.LATCH[EGS4.NP - 1] = EGS4.LATCH[EGS4.NP - 2];

						eta = EGS4.random01();
						eta = 2.0 * eta - 1.0;
						EGS4.W[EGS4.NP - 1] = eta;
						eta = (1.0 - eta) * (1.0 + eta);
						if (eta > 1.e-20) {
							eta = Math.sqrt(eta);
							// $SELECT-AZIMUTHAL-ANGLE(cphi,sphi);
							// while(rhophi2 <= 1)//LOOP
							do {
								// $RANDOMSET xphi;
								xphi = EGS4.random01();
								xphi = 2.0 * xphi - 1.0;
								xphi2 = xphi * xphi;
								// $RANDOMSET yphi;
								yphi = EGS4.random01();
								yphi2 = yphi * yphi;
								rhophi2 = xphi2 + yphi2;
							} while (rhophi2 > 1.0);
							rhophi2 = 1.0 / rhophi2;
							cphi = (xphi2 - yphi2) * rhophi2;
							sphi = 2.0 * xphi * yphi * rhophi2;

							EGS4.U[EGS4.NP - 1] = eta * cphi;
							EGS4.V[EGS4.NP - 1] = eta * sphi;
						} else {
							EGS4.U[EGS4.NP - 1] = 0.0;
							EGS4.V[EGS4.NP - 1] = 0.0;
							EGS4.W[EGS4.NP - 1] = 1.0;
						}
						// $AUSCALL($AUGERTRA);
						iarg = EGS4.$AUGERTRA;
						if (EGS4.iausfl[iarg] != 0) {
							eq.AUSGAB(iarg);
						}
					} else {
						EGS4.EDEP = EGS4.EDEP + Ei;
					}
					if (n_vac > 0)
						repeatLOOP = true;// goto :START: ;repeatLOOP
					else
						// added
						break;// EXIT;
				}
			}// if(!repeatLOOP)
				// " Sample transition number for this vacancy "
			if (!repeatLOOP) {
				eta = EGS4.random01();
				// DO k=first_transition(shell),last_transition(shell)-1 [
				for (k = first_transition[shell - 1]; k <= last_transition[shell - 1] - 1; k++) {
					eta = eta - EGS4.relaxation_prob[k - 1][iZ - 1];// (k,iZ);
					if (eta <= 0.0)
						break;
				}
				finalL = final_state[k - 1];// (k);
				finala = finalL;
				if (finalL < 100) {
					if (finalL < 10) {// "fuorescence"
						iql = 0;
						elcut = pkcut;
					} else {// "Coster-Kronig"
						finalL = finalL - 10;
						iql = -1;
						elcut = ekcut;
					}
					Ef = EGS4.binding_energies[finalL - 1][iZ - 1];// (finalL,iZ);
					Ex = Ei - Ef;
					n_vac = n_vac + 1;
					vac_array[n_vac - 1] = finalL;
					e_array[n_vac - 1] = Ef;
				} else {// "Auger"
					final1 = finalL / 100;
					final2 = finalL - final1 * 100;
					n_vac = n_vac + 1;
					vac_array[n_vac - 1] = final1;
					e_array[n_vac - 1] = EGS4.binding_energies[final1 - 1][iZ - 1];// (final1,iZ);
					n_vac = n_vac + 1;
					vac_array[n_vac - 1] = final2;
					e_array[n_vac - 1] = EGS4.binding_energies[final2 - 1][iZ - 1];// (final2,iZ);
					iql = -1;
					Ex = Ei - e_array[n_vac - 1] - e_array[n_vac - 2];
					elcut = ekcut;
				}
				if (Ex <= elcut) {// "Below cut-off"
					EGS4.EDEP = EGS4.EDEP + Ex;
				} else {
					EGS4.NP = EGS4.NP + 1;
					// $CHECK-STACK(np,'RELAX');
					if (EGS4.NP > EGS4.$MXSTACK) {
						EGS4.STOPPROGRAM = true;

						EGS4.seqStr = " ***************************************************"
								+ "  \n"
								+ " In subroutine "
								+ "RELAX"
								+ " stack size exceeded!"
								+ "  \n"
								+ " $MXSTACK = "
								+ EGS4.$MXSTACK
								+ " np = "
								+ EGS4.NP
								+ "  \n"
								+ " Increase $MXSTACK and try again "
								+ "  \n"
								+ " Terminating execution "
								+ "  \n"
								+ " ***************************************************";// +"  \n";
						// if(EGS4.iprint>2)
						eq.printSequence(EGS4.seqStr);

						return;// stop;
					}

					EGS4.IQ[EGS4.NP - 1] = iql;
					if (iql == 0) {
						EGS4.E[EGS4.NP - 1] = Ex;
					} else {
						EGS4.E[EGS4.NP - 1] = Ex + EGS4.RM;
					}
					// $TRANSFER PROPERTIES TO (np) FROM (np_old);
					// $TRANSFER PROPERTIES TO (np) FROM (np-1);
					EGS4.X[EGS4.NP - 1] = EGS4.X[EGS4.NP - 2];
					EGS4.Y[EGS4.NP - 1] = EGS4.Y[EGS4.NP - 2];
					EGS4.Z[EGS4.NP - 1] = EGS4.Z[EGS4.NP - 2];
					EGS4.IR[EGS4.NP - 1] = EGS4.IR[EGS4.NP - 2];
					EGS4.WT[EGS4.NP - 1] = EGS4.WT[EGS4.NP - 2];
					EGS4.DNEAR[EGS4.NP - 1] = EGS4.DNEAR[EGS4.NP - 2];
					EGS4.LATCH[EGS4.NP - 1] = EGS4.LATCH[EGS4.NP - 2];

					eta = EGS4.random01();
					eta = 2.0 * eta - 1.0;
					EGS4.W[EGS4.NP - 1] = eta;
					eta = (1.0 - eta) * (1.0 + eta);
					if (eta > 1.e-20) {
						eta = Math.sqrt(eta);
						// $SELECT-AZIMUTHAL-ANGLE(cphi,sphi);
						do {
							// $RANDOMSET xphi;
							xphi = EGS4.random01();
							xphi = 2.0 * xphi - 1.0;
							xphi2 = xphi * xphi;
							// $RANDOMSET yphi;
							yphi = EGS4.random01();
							yphi2 = yphi * yphi;
							rhophi2 = xphi2 + yphi2;
						} while (rhophi2 > 1.0);
						rhophi2 = 1.0 / rhophi2;
						cphi = (xphi2 - yphi2) * rhophi2;
						sphi = 2.0 * xphi * yphi * rhophi2;

						EGS4.U[EGS4.NP - 1] = eta * cphi;
						EGS4.V[EGS4.NP - 1] = eta * sphi;
					} else {
						EGS4.U[EGS4.NP - 1] = 0.0;
						EGS4.V[EGS4.NP - 1] = 0.0;
						EGS4.W[EGS4.NP - 1] = 1.0;
					}
					if (finala < 10) {
						// $AUSCALL($FLUORTRA);
						iarg = EGS4.$FLUORTRA;
						if (EGS4.iausfl[iarg] != 0) {
							eq.AUSGAB(iarg);
						}

					} else if (finala < 100) {
						// $AUSCALL($COSKROTRA);
						iarg = EGS4.$COSKROTRA;
						if (EGS4.iausfl[iarg] != 0) {
							eq.AUSGAB(iarg);
						}

					} else {
						// $AUSCALL($AUGERTRA);
						iarg = EGS4.$AUGERTRA;
						if (EGS4.iausfl[iarg] != 0) {
							eq.AUSGAB(iarg);
						}

					}
				}
			}// if(!repeatLOOP)
		}// loop

		return;
	}

	// "******************************************************************"
	// "                               National Research Council of Canada"
	/**
	 * Internally used. UPHI STANDS FOR 'UNIFORM PHI DISTRIBUTION'. SET COORDINATES FOR NEW PARTICLE OR RESET DIRECTION COSINES OF 
	 * OLD ONE.  GENERATE RANDOM AZIMUTH SELECTION AND REPLACE THE DIRECTION COSINES WITH THEIR NEW VALUES.
	 * @param IENTRY IENTRY
	 * @param LVL LVL
	 */
	protected static void UPHI(int IENTRY, int LVL) {
		// "                                                                  "
		// "******************************************************************"
		// "   UPHI STANDS FOR 'UNIFORM PHI DISTRIBUTION'.                    "
		// "   SET COORDINATES FOR NEW PARTICLE OR RESET DIRECTION COSINES OF "
		// "   OLD ONE.  GENERATE RANDOM AZIMUTH SELECTION AND REPLACE THE    "
		// "   DIRECTION COSINES WITH THEIR NEW VALUES.                       "
		// "******************************************************************"
		// ; Copyright NRC;

		// $COMIN-UPHI; "DEFAULT REPLACEMENT PRODUCES THE FOLLOWING:
		// "COMIN/DEBUG,EPCONT,STACK,UPHIIN,UPHIOT,RANDOM/;

		// "Input variables"
		// integer IENTRY,LVL; "entry switches"

		// "Local variables"
		//double RNNO38 = 0.0;// "random number for azimuthal angle selection"
		// double CTHET=0.0;//
		// "5/2*pi-THETA, used to evaluate cos(THETA) using the sine table"
		// double PHI=0.0;// "azimuthal scattering angle"
		// double CPHI=0.0;// "5/2*pi-PHI"
		// double A=0.0;
		// double B=0.0;
		// double C=0.0;// "direction cosines before rotation"
		// double SINPS2=0.0;// "SINPS2=A*A+B*B"
		// double SINPSI=0.0;// "Sqrt(SINPS2)"
		// double US=0.0;double VS=0.0;//
		// "x- and y- component of scattering vector"
		// double SINDEL=0.0;double COSDEL=0.0;
		// "aux. variables for the rotation algorithm"

		// $INTEGER
		int IARG = 0;// "index for eq.AUSGAB"
		//int LPHI = 0;
		//int LTHETA = 0;
		//int LCTHET = 0;
		//int LCPHI = 0;
		// "indeces for sine table"
		boolean ie3b = false;
		// $DEFINE-VARIABLES-FOR-SELECT-AZIMUTHAL-ANGLE;
		// REPLACE {$DEFINE-VARIABLES-FOR-SELECT-AZIMUTHAL-ANGLE;} WITH {;
		// $REAL xphi,xphi2,yphi,yphi2,rhophi2;
		double xphi = 0.0;
		double xphi2 = 0.0;
		double yphi = 0.0;
		double yphi2 = 0.0;
		double rhophi2 = 0.0;
		// };

		// save CTHET,PHI,CPHI,A,B,C,SINPS2,SINPSI,US,VS,SINDEL,COSDEL;

		// $AUSCALL($UPHIAUSB);
		// REPLACE {$AUSCALL(#);} WITH{
		// IARG={P1} ; IF (IAUSFL(IARG+1).NE.0) [CALL eq.AUSGAB(IARG);]
		IARG = EGS4.$UPHIAUSB;
		if (EGS4.iausfl[IARG] != 0) {
			eq.AUSGAB(IARG);
		}
		// } ;
		// GO TO (:UPHI:,:UPHI2:,:NRK:),IENTRY;
		if (IENTRY == 1) {// UPHI
							// "NOTE: AFB 88/12/12 ADDED SEMI-COLON, ELSE BUG WHEN OVERRIDING SIN"
							// "TABLE LOOK-UP"
							// $SET INTERVAL THETA,SINC;
							// NEW=> REPLACE {$EVALUATE#USING SIN(#);} WITH
							// {{P1}=sin({P2});}
			// NEW=> REPLACE {$SET INTERVAL#,SINC;} WITH {;}

			// Double dth=new
			// Double(EGS4.SINC1*EGS4.THETA+EGS4.SINC0);//System.out.println(dth);
			// LTHETA=dth.intValue();

			// $EVALUATE SINTHE USING SIN(THETA);
			// EGS4.SINTHE=EGS4.SIN1[LTHETA-1]*EGS4.THETA+EGS4.SIN0[LTHETA-1];
			EGS4.SINTHE = Math.sin(EGS4.THETA);

			CTHET = EGS4.PI5D2 - EGS4.THETA;
			// $SET INTERVAL CTHET,SINC;
			// Double dcth=new Double(EGS4.SINC1*CTHET+EGS4.SINC0);
			// LCTHET=dcth.intValue();
			// $EVALUATE COSTHE USING SIN(CTHET);
			// EGS4.COSTHE=EGS4.SIN1[LCTHET-1]*CTHET+EGS4.SIN0[LCTHET-1];
			EGS4.COSTHE = Math.sin(CTHET);

			// "   USE THE FOLLOWING ENTRY IF SINTHE AND COSTHE ARE ALREADY KNOWN."
			// "   SELECT PHI UNIFORMLY OVER THE INTERVAL (0,TWO PI). THEN USE    "
			// "   PWLF OF SIN FUNCTION TO GET SIN(PHI) AND COS(PHI).  THE COSINE "
			// "   IS GOTTEN BY COS(PHI)=SIN(9*PI/4 - PHI).                       "
			// UPHI2 is next
			do {
				// $RANDOMSET xphi;
				xphi = EGS4.random01();
				xphi = 2.0 * xphi - 1.0;
				xphi2 = xphi * xphi;
				// $RANDOMSET yphi;
				yphi = EGS4.random01();
				yphi2 = yphi * yphi;
				rhophi2 = xphi2 + yphi2;
			} while (rhophi2 > 1.0);
			rhophi2 = 1.0 / rhophi2;
			EGS4.COSPHI = (xphi2 - yphi2) * rhophi2;
			EGS4.SINPHI = 2.0 * xphi * yphi * rhophi2;

			ie3b = true;
		} else if (IENTRY == 2) {// UPHI2
									// " It is much faster to use the box method for azimuthal angle selection"
									// " than the following                                                   "
									// " $RANDOMSET RNNO38;
									// " PHI=RNNO38*TWOPI;$SET INTERVAL
									// PHI,SINC;
									// " $EVALUATE SINPHI USING SIN(PHI);
									// " CPHI=PI5D2-PHI;$SET INTERVAL CPHI,SINC;
									// " $EVALUATE COSPHI USING SIN(CPHI);

			// $SELECT-AZIMUTHAL-ANGLE(cosphi,sinphi);
			// "Macro for azimuthal angle selection
			// "using a sampling within a box method
			// "Choose a point randomly within a box such that
			// "-1 <= x <= 1 and 0 <= y < = 1
			// "Reject the set if it lies without the inscribed unit semicircle
			// centered
			// "at (x,y) = (0,0)
			// "once out of the loop, use the trigonimetric relations (TeX
			// notation)
			// "\cos 2\phi = (x^2 - y^2)/(x^2 + y^2)
			// "\sin 2\phi = 2xy/(x^2 + y^2)
			// REPLACE {$SELECT-AZIMUTHAL-ANGLE(#,#);} WITH
			// {
			// ;
			// while(rhophi2 <= 1)//LOOP
			do {
				// $RANDOMSET xphi;
				xphi = EGS4.random01();
				xphi = 2.0 * xphi - 1.0;
				xphi2 = xphi * xphi;
				// $RANDOMSET yphi;
				yphi = EGS4.random01();
				yphi2 = yphi * yphi;
				rhophi2 = xphi2 + yphi2;
			} while (rhophi2 > 1.0);
			rhophi2 = 1.0 / rhophi2;
			EGS4.COSPHI = (xphi2 - yphi2) * rhophi2;
			EGS4.SINPHI = 2.0 * xphi * yphi * rhophi2;
			// }
			// ;
			// "   USE THE FOLLOWING ENTRY FOR THE SECOND OF TWO PARTICLES WHEN WE"
			// "   KNOW TWO PARTICLES HAVE A RELATIONSHIP IN THEIR CORRECTIONS.   "
			// "   NOTE: SINTHE AND COSTHE CAN BE CHANGED OUTSIDE THROUGH COMMON. "
			// "   LVL IS A PARAMETER TELLING WHICH PARTICLES TO WORK WITH.       "
			// "   THETA (SINTHE AND COSTHE) ARE ALWAYS RELATIVE TO THE DIRECTION "
			// "   OF THE INCIDENT PARTICLE BEFORE ITS DIRECTION WAS ADJUSTED.    "
			// "   THUS WHEN TWO PARTICLES NEED TO HAVE THEIR DIRECTIONS COMPUTED,"
			// "   THE ORIGINAL INCIDENT DIRECTION IS SAVED IN THE VARIABLE A,B,C "
			// "   SO THAT IT CAN BE USED ON BOTH CALLS."

			// "   LVL=1 -- OLD PARTICLE, SAVE ITS DIRECTION AND ADJUST IT"
			// "   LVL=2 -- NEW PARTICLE. ADJUST DIRECTION USING SAVED A,B,C"
			// "   LVL=3 -- BREMSSTRAHLUNG GAMMA.  SAVE ELECTRON DIRECTION (NEXT  "
			// "   TO TOP OF STACK), AND THEN ADJUST GAMMA DIRECTION."

			ie3b = true;
		} else if (IENTRY == 3) {// NRK
									// GO TO
									// (:OLD-PARTICLE:,:NEW-PARTICLE:,:BREMS-GAMMA:),LVL;
									// "   LVL OUT-OF-BOUNDS IF HERE" GO TO
									// :ERROR:;
			ie3b = true;
		} else {// ERROR
				// "REACH THIS POINT IF EITHER IENTRY OR LVL NE 1,2, OR 3"

			EGS4.STOPPROGRAM = true;

			EGS4.seqStr = " STOPPED IN UPHI WITH IENTRY,LVL=" + IENTRY + "  "
					+ LVL;
			// if(EGS4.iprint>2)
			eq.printSequence(EGS4.seqStr);

			return;// stop
		}

		if (ie3b) {
			if (LVL == 1) {// :OLD-PARTICLE:
							// A=U(NP);B=V(NP);C=W(NP);
				A = EGS4.U[EGS4.NP - 1];
				B = EGS4.V[EGS4.NP - 1];
				C = EGS4.W[EGS4.NP - 1];
				// GO TO :ADJUST:;

			} else if (LVL == 2) {// :NEW-PARTICLE
									// "   SEE H.H. NAGEL DISSERTATION FOR COORDINATE SYSTEM DESCRIPTION. "
									// "   A ROTATION IS PERFORMED TO TRANSFORM DIRECTION COSINES OF THE  "
									// "   PARTICLE BACK TO THE PHYSICAL FRAME (FROM THE TRANSPORT FRAME) "
									// $TRANSFER PROPERTIES TO (NP) FROM (NP-1);
				// $TRANSFER PROPERTIES TO (1) FROM I;
				// X(1)=XI;Y(1)=YI;Z(1)=ZI;IR(1)=IRI;
				// WT(1)=WTI;DNEAR(1)=DNEARI;LATCH(1)=LATCHI
				EGS4.X[EGS4.NP - 1] = EGS4.X[EGS4.NP - 2];
				EGS4.Y[EGS4.NP - 1] = EGS4.Y[EGS4.NP - 2];
				EGS4.Z[EGS4.NP - 1] = EGS4.Z[EGS4.NP - 2];
				EGS4.IR[EGS4.NP - 1] = EGS4.IR[EGS4.NP - 2];
				EGS4.WT[EGS4.NP - 1] = EGS4.WT[EGS4.NP - 2];
				EGS4.DNEAR[EGS4.NP - 1] = EGS4.DNEAR[EGS4.NP - 2];
				EGS4.LATCH[EGS4.NP - 1] = EGS4.LATCH[EGS4.NP - 2];

			} else if (LVL == 3) {// :BREMS-GAMMA:
									// A=U(NP-1);B=V(NP-1);C=W(NP-1);
				A = EGS4.U[EGS4.NP - 2];
				B = EGS4.V[EGS4.NP - 2];
				C = EGS4.W[EGS4.NP - 2];
				// $TRANSFER PROPERTIES TO (NP) FROM (NP-1);
				EGS4.X[EGS4.NP - 1] = EGS4.X[EGS4.NP - 2];
				EGS4.Y[EGS4.NP - 1] = EGS4.Y[EGS4.NP - 2];
				EGS4.Z[EGS4.NP - 1] = EGS4.Z[EGS4.NP - 2];
				EGS4.IR[EGS4.NP - 1] = EGS4.IR[EGS4.NP - 2];
				EGS4.WT[EGS4.NP - 1] = EGS4.WT[EGS4.NP - 2];
				EGS4.DNEAR[EGS4.NP - 1] = EGS4.DNEAR[EGS4.NP - 2];
				EGS4.LATCH[EGS4.NP - 1] = EGS4.LATCH[EGS4.NP - 2];
			} else {
				// ERROR
				// "REACH THIS POINT IF EITHER IENTRY OR LVL NE 1,2, OR 3"
				EGS4.STOPPROGRAM = true;

				EGS4.seqStr = " STOPPED IN UPHI WITH IENTRY,LVL=" + IENTRY
						+ "  " + LVL;
				// if(EGS4.iprint>2)
				eq.printSequence(EGS4.seqStr);
				return;// stop
			}

			// adjust
			SINPS2 = A * A + B * B;
			// "   If SINPS2 is small, no rotation is needed    "
			if (SINPS2 < 1.0E-20) {// "small polar angle case"
				EGS4.U[EGS4.NP - 1] = EGS4.SINTHE * EGS4.COSPHI;
				EGS4.V[EGS4.NP - 1] = EGS4.SINTHE * EGS4.SINPHI;
				EGS4.W[EGS4.NP - 1] = C * EGS4.COSTHE; // "fixed March 2001 from =COSTHE"
			} // "end small polar angle case"
			else {// "large polar angle case"
				SINPSI = Math.sqrt(SINPS2);
				US = EGS4.SINTHE * EGS4.COSPHI;
				VS = EGS4.SINTHE * EGS4.SINPHI;
				SINDEL = B / SINPSI;
				COSDEL = A / SINPSI;
				EGS4.U[EGS4.NP - 1] = C * COSDEL * US - SINDEL * VS + A
						* EGS4.COSTHE;
				EGS4.V[EGS4.NP - 1] = C * SINDEL * US + COSDEL * VS + B
						* EGS4.COSTHE;
				EGS4.W[EGS4.NP - 1] = -SINPSI * US + C * EGS4.COSTHE;
			}// "end large polar angle case"

			// $AUSCALL($UPHIAUSA);
			IARG = EGS4.$UPHIAUSA;
			if (EGS4.iausfl[IARG] != 0) {
				eq.AUSGAB(IARG);
			}

		}
	}// "END OF SUBROUTINE UPHI" END;

	// "******************************************************************"
	// "                               National Research Council of Canada"
	/**
	 * Called by ELECTR. Handle annihilation at rest.
	 */
	protected static void ANNIH_AT_REST() {
		// "                                                                  "
		// " It is handy to be able to initiate annihilation at rest from     "
		// " places other than the electron discard section (e.g. eq.AUSGAB)     "
		// " Annihilation at rest takes a sufficent amount of time to not     "
		// " have any real benefit from this code being inline in the         "
		// " ELECTR subroutine.                                               "
		// " I. Kawrakow, June 2005.                                          "
		// "******************************************************************"
		// ; Copyright NRC;
		// $IMPLICIT-NONE;
		// $COMIN-ANNIH-ATREST;
		// $REAL costhe,sinthe,cphi,sphi;
		double cphi = 0.0;
		double sphi = 0.0;
		// $INTEGER ibr;
		int ibr = 0;
		// $DEFINE-VARIABLES-FOR-SELECT-AZIMUTHAL-ANGLE;
		double xphi = 0.0;
		double xphi2 = 0.0;
		double yphi = 0.0;
		double yphi2 = 0.0;
		double rhophi2 = 0.0;

		EGS4.NPold = EGS4.NP;
		// $CHECK-STACK(np+2*nbr_split-1,'ANNIH_AT_REST');
		// $CHECK-STACK(np+1,'MOLLER');
		if (EGS4.NP + 2 * EGS4.nbr_split - 1 > EGS4.$MXSTACK) {
			EGS4.STOPPROGRAM = true;

			EGS4.seqStr = " ***************************************************"
					+ "  \n"
					+ " In subroutine "
					+ "ANNIH_AT_REST"
					+ " stack size exceeded!"
					+ "  \n"
					+ " $MXSTACK = "
					+ EGS4.$MXSTACK
					+ " np = "
					+ EGS4.NP
					+ "  \n"
					+ " Increase $MXSTACK and try again "
					+ "  \n"
					+ " Terminating execution "
					+ "  \n"
					+ " ***************************************************";// +"  \n";
			// if(EGS4.iprint>2)
			eq.printSequence(EGS4.seqStr);

			return;// stop;
		}

		if (EGS4.nbr_split > 1) {
			EGS4.WT[EGS4.NP - 1] = EGS4.WT[EGS4.NP - 1] / EGS4.nbr_split;
		}
		// " nbr_split > 1 means user wants to use radiative "
		// " splitting => produce 2*nbr_split annihilation   "
		// " photons at once                                 "
		for (ibr = 1; ibr <= EGS4.nbr_split; ibr++) {
			// "Pick random direction for first gamma
			// $RANDOMSET costhe;
			EGS4.COSTHE = EGS4.random01();
			EGS4.COSTHE = 2.0 * EGS4.COSTHE - 1.0;
			EGS4.SINTHE = Math.sqrt(Math.max(0.0, (1.0 - EGS4.COSTHE)
					* (1.0 + EGS4.COSTHE)));
			// $SELECT-AZIMUTHAL-ANGLE(cphi,sphi);
			do {
				xphi = EGS4.random01();
				xphi = 2.0 * xphi - 1.0;
				xphi2 = xphi * xphi;
				yphi = EGS4.random01();
				yphi2 = yphi * yphi;
				rhophi2 = xphi2 + yphi2;
			} while (rhophi2 > 1.0);
			rhophi2 = 1.0 / rhophi2;
			cphi = (xphi2 - yphi2) * rhophi2;
			sphi = 2.0 * xphi * yphi * rhophi2;

			EGS4.E[EGS4.NP - 1] = EGS4.PRM;
			EGS4.IQ[EGS4.NP - 1] = 0;
			// @=========================================================
			int ip = 0;
			// IF( ibr = 1 ) [ ip = npold; ] ELSE [ ip = np-1; ]
			if (ibr == 1) {
				ip = EGS4.NPold;
			} else {
				ip = EGS4.NP - 1;
			}
			// ===========================================================
			// $TRANSFER PROPERTIES TO (np) FROM (npold);
			// $TRANSFER PROPERTIES TO (np) FROM (ip);
			EGS4.X[EGS4.NP - 1] = EGS4.X[ip - 1];
			EGS4.Y[EGS4.NP - 1] = EGS4.Y[ip - 1];
			EGS4.Z[EGS4.NP - 1] = EGS4.Z[ip - 1];
			EGS4.IR[EGS4.NP - 1] = EGS4.IR[ip - 1];
			EGS4.WT[EGS4.NP - 1] = EGS4.WT[ip - 1];
			EGS4.DNEAR[EGS4.NP - 1] = EGS4.DNEAR[ip - 1];
			EGS4.LATCH[EGS4.NP - 1] = EGS4.LATCH[ip - 1];

			EGS4.U[EGS4.NP - 1] = EGS4.SINTHE * cphi;
			EGS4.V[EGS4.NP - 1] = EGS4.SINTHE * sphi;
			EGS4.W[EGS4.NP - 1] = EGS4.COSTHE;
			EGS4.NP = EGS4.NP + 1;
			// e(np) = prm; iq(np) = 0;
			EGS4.E[EGS4.NP - 1] = EGS4.PRM;
			EGS4.IQ[EGS4.NP - 1] = 0;
			// $TRANSFER PROPERTIES TO (np) FROM (npold);
			// @========================================================
			// $TRANSFER PROPERTIES TO (np) FROM (np-1);
			// @========================================================
			EGS4.X[EGS4.NP - 1] = EGS4.X[EGS4.NP - 2];
			EGS4.Y[EGS4.NP - 1] = EGS4.Y[EGS4.NP - 2];
			EGS4.Z[EGS4.NP - 1] = EGS4.Z[EGS4.NP - 2];
			EGS4.IR[EGS4.NP - 1] = EGS4.IR[EGS4.NP - 2];
			EGS4.WT[EGS4.NP - 1] = EGS4.WT[EGS4.NP - 2];
			EGS4.DNEAR[EGS4.NP - 1] = EGS4.DNEAR[EGS4.NP - 2];
			EGS4.LATCH[EGS4.NP - 1] = EGS4.LATCH[EGS4.NP - 2];

			EGS4.U[EGS4.NP - 1] = -EGS4.U[EGS4.NP - 2];
			EGS4.V[EGS4.NP - 1] = -EGS4.V[EGS4.NP - 2];// (np-1);
			EGS4.W[EGS4.NP - 1] = -EGS4.W[EGS4.NP - 2];
			EGS4.NP = EGS4.NP + 1;
		}
		EGS4.NP = EGS4.NP - 1;
	}// return; end;

	// "****************************************************************************"
	/**
	 * Called by COMPT. Samples a double Compton scattering event and puts the resulting particles (2 photons and 1 electron) on the stack. 
	 * wo is the incident photon energy in units of the electron rest energy; ie is the energy index that was needed and determined in COMPT in 
	 * order to decide if single or double Compton scattering is to be simulated.
	 * @param wo wo
	 * @param ie ie
	 */
	protected static void sample_double_compton(double wo, int ie) {
		// "****************************************************************************"
		// "
		// " Samples a double Compton scattering event and puts the resulting
		// particles
		// " (2 photons and 1 electron) on the stack.
		// " wo is the incident photon energy in units of the electron rest
		// energy
		// " ie is the energy index that was needed and determined in COMPT in
		// order
		// " to decide if single or double Compton scattering is to be
		// simulated.
		// "
		// " I Kawrakow, September 2005.
		// "============================================================================"
		// implicit none;
		// $REAL wo;
		// $INTEGER ie;
		// $declare_max_medium;
		// REPLACE {$declare_max_medium;} WITH {;};

		double rnno = 0.0;
		double y1 = 0.0;
		double y2 = 0.0;
		double yw = 0.0;
		double s_max = 0.0;
		double asamp = 0.0;
		double alpha1 = 0.0;
		double yb1 = 0.0;
		double yb2 = 0.0;
		double yb3 = 0.0;
		double yb4 = 0.0;
		double dy1 = 0.0;
		double dy2 = 0.0;
		double dy3 = 0.0;
		double dy4 = 0.0;
		double wo_save = 0.0;
		double Vol = 0.0;
		double cost1 = 0.0;
		double cost2 = 0.0;
		double w1 = 0.0;
		double w2 = 0.0;
		double cost12 = 0.0;
		double cphi = 0.0;
		double acphi = 0.0;
		double a1 = 0.0;
		double b1 = 0.0;
		double z1 = 0.0;
		double z2 = 0.0;
		double facct = 0.0;
		double aux = 0.0;
		double ax = 0.0;
		double bx = 0.0;
		double w1_max = 0.0;
		double w1t = 0.0;
		double ww1tot = 0.0;
		double zz = 0.0;
		double facw1 = 0.0;
		double k1 = 0.0;
		double k2 = 0.0;
		double k3 = 0.0;
		double k1i = 0.0;
		double k2i = 0.0;
		double k3i = 0.0;
		double k1p = 0.0;
		double k2p = 0.0;
		double k3p = 0.0;
		double k1pi = 0.0;
		double k2pi = 0.0;
		double k3pi = 0.0;
		double a = 0.0;
		double b = 0.0;
		double c = 0.0;
		double Ac = 0.0;
		double Bc = 0.0;
		double rho = 0.0;
		double s = 0.0;
		double xx = 0.0;
		double Xc = 0.0;
		double px1 = 0.0;
		double py1 = 0.0;
		double pz1 = 0.0;
		double px2 = 0.0;
		double py2 = 0.0;
		double pz2 = 0.0;
		double pxe = 0.0;
		double pye = 0.0;
		double pze = 0.0;
		double pp = 0.0;
		double Ep = 0.0;
		double sindel = 0.0;
		double cosdel = 0.0;
		double sinpsi = 0.0;
		double sphi = 0.0;
		double sint1 = 0.0;
		double sint2 = 0.0;
		// $DEFINE-VARIABLES-FOR-SELECT-AZIMUTHAL-ANGLE;
		// REPLACE {$DEFINE-VARIABLES-FOR-SELECT-AZIMUTHAL-ANGLE;} WITH {;
		double xphi = 0.0;
		double xphi2 = 0.0;
		double yphi = 0.0;
		double yphi2 = 0.0;
		double rhophi2 = 0.0;
		int ibin = 0;
		int nbox = 0;
		int ix = 0;
		double woo = wo;
		boolean dcloop = false;
		boolean dcbegin = false;
		// "****************************************************************************"
		// $CHECK-STACK(np+2,'sample_double_compton');
		if (EGS4.NP + 2 > EGS4.$MXSTACK) {
			EGS4.STOPPROGRAM = true;

			EGS4.seqStr = " ***************************************************"
					+ "  \n"
					+ " In subroutine "
					+ "sample_double_compton"
					+ " stack size exceeded!"
					+ "  \n"
					+ " $MXSTACK = "
					+ EGS4.$MXSTACK
					+ " np = "
					+ EGS4.NP
					+ "  \n"
					+ " Increase $MXSTACK and try again "
					+ "  \n"
					+ " Terminating execution "
					+ "  \n"
					+ " ***************************************************";// +"  \n";
			// if(EGS4.iprint>2)
			eq.printSequence(EGS4.seqStr);

			return;// stop;
		}
		// *******************
		// Pick the cell (box) within the final particles phase-spece from which
		// to
		// sample energies and directions.
		// radc_startb(ie) is the first cell belonging to this energy.
		// Number of cells for this energy is nbox =
		// radc_startb(ie+1)-radc_startb(ie)
		// The data is in form of an alias table, i.e. we pick a cell randomly
		// and then based on another random number stay in the cell or go to the
		// cell that has filled up the probability for the alias table.
		// *******************/
		rnno = EGS4.random01();
		// radc_startb=new int[$RADC_NE+1];//$0-RADC_NE);
		nbox = EGS4.radc_startb[ie + 1] - EGS4.radc_startb[ie];
		Double ddd = new Double(EGS4.radc_startb[ie] + rnno * nbox);
		ibin = ddd.intValue();
		rnno = EGS4.random01();
		// if( rnno > radc_fdat(ibin) ) ibin = radc_startb(ie) +
		// radc_bins(ibin);
		if (rnno > EGS4.radc_fdat[ibin - 1])
			ibin = EGS4.radc_startb[ie] + EGS4.radc_bins[ibin - 1];
		// ******************
		// Set the cell boundaries and the maximum value of the cross section
		// within this cell
		// ******************/
		s_max = EGS4.radc_Smax[ibin - 1];// radc_Smax(ibin);
		ix = EGS4.radc_startx[ie];// radc_startx=new
									// int[$RADC_NE+1];//$0-RADC_NE);
		// radc_x=new double[$RADC_NX];radc_ixmin1=new int[$RADC_NBOX];
		yb1 = EGS4.radc_x[ix + EGS4.radc_ixmin1[ibin - 1] - 1];
		dy1 = EGS4.radc_x[ix + EGS4.radc_ixmax1[ibin - 1] - 1] - yb1;
		yb2 = EGS4.radc_x[ix + EGS4.radc_ixmin2[ibin - 1] - 1];
		dy2 = EGS4.radc_x[ix + EGS4.radc_ixmax2[ibin - 1] - 1] - yb2;
		yb3 = EGS4.radc_x[ix + EGS4.radc_ixmin3[ibin - 1] - 1];
		dy3 = EGS4.radc_x[ix + EGS4.radc_ixmax3[ibin - 1] - 1] - yb3;
		yb4 = EGS4.radc_x[ix + EGS4.radc_ixmin4[ibin - 1] - 1];
		dy4 = EGS4.radc_x[ix + EGS4.radc_ixmax4[ibin - 1] - 1] - yb4;
		Vol = dy1 * dy2 * dy3 * dy4;

		// :DC-BEGIN:;
		while (true) {
			dcbegin = false;
			// *******************
			// We use the energy from the energy grid on which the data was
			// created
			// to sample the 4 random numbers needed to obtain the final
			// energies and
			// directions. Once we have these, we will use the actual energy wo
			// to
			// determine energies and directions.
			// *******************/
			wo_save = wo;
			woo = wo;
			woo = EGS4.radc_emin * Math.exp(EGS4.radc_dle * ie);
			// wo = EGS4.radc_emin*Math.exp(EGS4.radc_dle*ie);
			asamp = 0.25 * woo * (1.0 + woo);
			// asamp = 0.25*wo*(1.0+wo);
			alpha1 = Math.log(0.5 * (1.0 + Math.sqrt(1.0 + 4.0 * asamp)));

			// :DC-LOOP:LOOP
			while (true) {
				dcloop = false;
				// *************************
				// y1 selects the cosine of the polar angle of the lower energy
				// photon
				// y2 selects the cosine of the polar angle of the higher energy
				// photon
				// yw selects the energy of the lower energy photon
				// yphi selects the azimuthal angle between the 2 photon
				// directions
				// all other kinematic variables follow from energy-momentum
				// conservation
				// We pick y1,y2,yw and yphi within the selected cell.
				// *************************/
				rnno = EGS4.random01();
				y1 = yb1 + dy1 * rnno;
				facct = 2.0 * y1;
				y1 = y1 * y1;
				rnno = EGS4.random01();
				y2 = yb2 + dy2 * rnno;
				rnno = EGS4.random01();
				yw = yb3 + dy3 * rnno;
				rnno = EGS4.random01();
				yphi = yb4 + dy4 * rnno;

				// **********************
				// Set cost1, cost2 and the corresponding transformation factor
				// facct
				// **********************/
				if (y2 > 0.0) {
					z2 = 2.0 * Math.exp(y2 * alpha1) - 1.0;
					cost2 = 1.0 - (z2 * z2 - 1.0) / (2.0 * asamp);
					z1 = z2 / (1.0 + y1 * (z2 - 1.0));
					cost1 = 1.0 - 2.0 * (z1 * z1 - 1.0) / (z2 * z2 - 1.0);
				} else {
					z2 = 1.0;
					z1 = 1.0;
					cost2 = 1.0;
					cost1 = 2.0 * y1 - 1.0;
				}
				facct = facct * alpha1 / asamp * z1 * z1 * z1;

				// *********************
				// Set the azimuthal angle and the corresponding transformation
				// factor
				// and calculate the cosine of the angle between the 2 photons
				// cost12
				// **********************/
				acphi = Math
						.sqrt((1.0 - cost1 * cost1) * (1.0 - cost2 * cost2));
				a1 = 1.0 + woo * (1.0 + woo) * (1.0 - cost1) * (1.0 - cost2);
				// a1 = 1.0 + wo*(1.0+wo)*(1.0-cost1)*(1.0-cost2);
				b1 = woo * acphi;
				// b1 = wo*acphi;
				if (Math.abs(yphi - 0.5) > 1.e-4) {
					aux = Math.tan(3.14159265358979323846 * yphi);
					aux = aux * aux;
					cphi = (a1 - b1 - (a1 + b1) * aux)
							/ (a1 - b1 + (a1 + b1) * aux);
				} else {
					cphi = -1.0;
				}
				facct = facct * (a1 + b1 * cphi) / Math.sqrt(a1 * a1 - b1 * b1);
				cost12 = cost1 * cost2 + acphi * cphi;

				// ******************
				// The maximum possible energy w1m for the lower energy photon
				// *******************/
				ax = 2.0 + woo * (1.0 - cost1) + woo * (1.0 - cost2);
				bx = woo * (1.0 - cost12) / (ax * ax);
				// ax = 2.0 + wo*(1.0-cost1) + wo*(1.0-cost2);
				// bx = wo*(1.0-cost12)/(ax*ax);

				if (bx > 1.e-3) {
					w1_max = woo * (1.0 - Math.sqrt(1.0 - 4.0 * bx))
							/ (2.0 * bx * ax);
				} else {
					w1_max = woo / ax
							* (1.0 + bx * (1.0 + bx * (2.0 + 5.0 * bx)));
				}
				w1t = woo / (1.0 + woo * (1.0 - cost1));
				// if( bx > 1.e-3 ) { w1_max =
				// wo*(1.0-Math.sqrt(1.0-4.0*bx))/(2.0*bx*ax); }
				// else { w1_max = wo/ax*(1.0 + bx*(1.0 + bx*(2.0 + 5.0*bx))); }
				// w1t = wo/(1.0 + wo*(1.0-cost1));

				if (w1_max <= EGS4.radc_dw) {// goto :DC-LOOP:;
					dcloop = true;
				}
				if (!dcloop) {
					// ******************
					// Set the energy w1 of the lower energy photon and the
					// corresponding transformation factor
					// *******************/
					ww1tot = Math.log(w1_max * (w1t - EGS4.radc_dw)
							/ (EGS4.radc_dw * (w1t - w1_max)));
					zz = Math.exp(yw * ww1tot);
					w1 = zz * w1t * EGS4.radc_dw
							/ (w1t + (zz - 1.0) * EGS4.radc_dw);
					facw1 = w1 * (w1t - w1) * ww1tot / w1t;

					// ******************
					// Calculate the energy of the second photon
					// /******************
					w2 = (woo - w1 * (1.0 + woo * (1.0 - cost1)))
							/ (1.0 + woo * (1.0 - cost2) - w1 * (1.0 - cost12));
					// w2 = (wo -
					// w1*(1.0+wo*(1.0-cost1)))/(1.0+wo*(1.0-cost2)-w1*(1.0-cost12));
					if (w1 > w2) // goto :DC-LOOP:;
					{
						dcloop = true;
					}
					if (!dcloop) {
						// ******************
						// Now calculate the cross section
						// *****************/
						k1 = w1;
						k2 = w2;
						// k3 = -wo;
						k3 = -woo;
						k1p = -w1
								* (1.0 + woo * (1.0 - cost1) - w2
										* (1.0 - cost12));
						k2p = -w2
								* (1.0 + woo * (1.0 - cost2) - w1
										* (1.0 - cost12));
						k3p = woo
								* (1.0 - w1 * (1.0 - cost1) - w2
										* (1.0 - cost2));
						// k1p = -w1*(1.0 + wo*(1.0-cost1) - w2*(1.0-cost12));
						// k2p = -w2*(1.0 + wo*(1.0-cost2) - w1*(1.0-cost12));
						// k3p = wo*(1.0 - w1*(1.0-cost1) - w2*(1.0-cost2));

						k1i = 1.0 / k1;
						k2i = 1.0 / k2;
						k3i = 1.0 / k3;
						k1pi = 1.0 / k1p;
						k2pi = 1.0 / k2p;
						k3pi = 1.0 / k3p;
						a = k1i + k2i + k3i;
						b = k1pi + k2pi + k3pi;
						c = k1i * k1pi + k2i * k2pi + k3i * k3pi;
						xx = k1 + k2 + k3;
						zz = k1 * k1p + k2 * k2p + k3 * k3p;
						Ac = k1 * k2 * k3;
						Bc = k1p * k2p * k3p;
						rho = k1 * k1pi + k1p * k1i + k2 * k2pi + k2p * k2i
								+ k3 * k3pi + k3p * k3i;
						Xc = 2.0
								* (a * b - c)
								* ((a + b) * (xx + 2.0) - (a * b - c) - 8.0)
								- 2.0
								* xx
								* (a * a + b * b)
								- 8.0
								* c
								+ 4.0
								* xx
								/ (Ac * Bc)
								* ((Ac + Bc) * (xx + 1.0) - (a * Ac + b * Bc)
										* (2.0 + zz * (1.0 - xx) / xx) + xx
										* xx * (1.0 - zz) + 2.0 * zz) - 2.0
								* rho * (a * b + c * (1.0 - xx));
						s = Xc
								* Vol
								* w1
								* w2
								* facct
								* facw1
								/ (woo * (1.0 + woo * (1.0 - cost2) - w1
										* (1.0 - cost12)));
						// s =
						// Xc*Vol*w1*w2*facct*facw1/(wo*(1.0+wo*(1.0-cost2)-w1*(1.0-cost12)));
						rnno = EGS4.random01();
						if (rnno * s_max <= s)
							break;
					}// if(!dcloop)
				}// if(!dcloop)
			}// WHILE ( rnno*s_max > s );

			// *****************
			// OK, now the y1, y2, yw and yphi have been accepted.
			// Calculate the energies and directions of the final particles
			// from the actual incident photon energy wo.
			// ****************/
			woo = wo_save;// and woo=wo so the following are correct
			// wo = wo_save;
			asamp = 0.25 * woo * (1.0 + woo);
			// asamp = 0.25*wo*(1.0+wo);
			alpha1 = Math.log(0.5 * (1.0 + Math.sqrt(1.0 + 4.0 * asamp)));

			z2 = 2.0 * Math.exp(y2 * alpha1) - 1.0;
			cost2 = 1.0 - (z2 * z2 - 1.0) / (2.0 * asamp);

			z1 = z2 / (1.0 + y1 * (z2 - 1.0));
			cost1 = 1.0 - 2.0 * (z1 * z1 - 1.0) / (z2 * z2 - 1.0);

			acphi = Math.sqrt((1.0 - cost1 * cost1) * (1.0 - cost2 * cost2));
			a1 = 1 + woo * (1.0 + woo) * (1.0 - cost1) * (1.0 - cost2);
			b1 = woo * acphi;
			// a1 = 1 + wo*(1.0+wo)*(1.0-cost1)*(1.0-cost2);
			// b1 = wo*acphi;

			if (Math.abs(yphi - 0.5) > 1.e-4) {
				aux = Math.tan(3.14159265358979323846 * yphi);
				aux = aux * aux;
				cphi = (a1 - b1 - (a1 + b1) * aux)
						/ (a1 - b1 + (a1 + b1) * aux);
			} else {
				cphi = -1.0;
			}
			cost12 = cost1 * cost2 + acphi * cphi;

			ax = 2.0 + woo * (1.0 - cost1) + woo * (1.0 - cost2);
			bx = woo * (1.0 - cost12) / (ax * ax);
			if (bx > 1.e-3) {
				w1_max = woo * (1.0 - Math.sqrt(1.0 - 4.0 * bx))
						/ (2.0 * bx * ax);
			} else {
				w1_max = woo / ax * (1.0 + bx * (1.0 + bx * (2.0 + 5.0 * bx)));
			}
			w1t = woo / (1 + woo * (1.0 - cost1));
			// ax = 2.0 + wo*(1.0-cost1) + wo*(1.0-cost2);
			// bx = wo*(1.0-cost12)/(ax*ax);
			// if( bx > 1.e-3 ) { w1_max =
			// wo*(1.0-Math.sqrt(1.0-4.0*bx))/(2.0*bx*ax); }
			// else { w1_max = wo/ax*(1.0 + bx*(1.0 + bx*(2.0 + 5.0*bx))); }
			// w1t = wo/(1 + wo*(1.0-cost1));

			if (w1_max > EGS4.radc_dw) {
				ww1tot = Math.log(w1_max * (w1t - EGS4.radc_dw)
						/ (EGS4.radc_dw * (w1t - w1_max)));
				zz = Math.exp(yw * ww1tot);
				w1 = zz * w1t * EGS4.radc_dw
						/ (w1t + (zz - 1.0) * EGS4.radc_dw);
				w2 = (woo - w1 * (1.0 + woo * (1.0 - cost1)))
						/ (1.0 + woo * (1.0 - cost2) - w1 * (1.0 - cost12));
				// w2 = (wo -
				// w1*(1.0+wo*(1.0-cost1)))/(1.0+wo*(1.0-cost2)-w1*(1.0-cost12));
				if (w1 < w2) {
					sphi = Math.sqrt(1.0 - cphi * cphi);
					sint1 = Math.sqrt(1.0 - cost1 * cost1);
					sint2 = Math.sqrt(1.0 - cost2 * cost2);
				} else {
					// goto :DC-BEGIN:;
					dcbegin = true;
				}
			} else {
				// goto :DC-BEGIN:;
				dcbegin = true;
			}

			// ################
			if (!dcbegin)
				break;
			// #################

		}// dcbegin
			// ************************************
			// Now set the final momenta in the frame where the
			// incident photon moves along the positive z-axis.
			// ************************************/
		px1 = w1 * sint1 * cphi;
		py1 = w1 * sint1 * sphi;
		pz1 = w1 * cost1;
		px2 = w2 * sint2;
		py2 = 0.0;
		pz2 = w2 * cost2;
		// "Pick now another azimuthal angle and rotate the above"
		// $SELECT-AZIMUTHAL-ANGLE(cphi,sphi);
		do {
			// $RANDOMSET xphi;
			xphi = EGS4.random01();
			xphi = 2.0 * xphi - 1.0;
			xphi2 = xphi * xphi;
			// $RANDOMSET yphi;
			yphi = EGS4.random01();
			yphi2 = yphi * yphi;
			rhophi2 = xphi2 + yphi2;
		} while (rhophi2 > 1.0);
		rhophi2 = 1.0 / rhophi2;
		cphi = (xphi2 - yphi2) * rhophi2;
		sphi = 2.0 * xphi * yphi * rhophi2;

		aux = px1 * sphi;
		px1 = px1 * cphi - py1 * sphi;
		py1 = aux + py1 * cphi;
		py2 = sphi * px2;
		px2 = cphi * px2;
		// "Get the electron momentum from momentum conservation"
		pxe = -px1 - px2;
		pye = -py1 - py2;
		pze = woo - pz1 - pz2;
		Ep = woo - w1 - w2 + 1.0;
		pp = 1.0 / Math.sqrt(pxe * pxe + pye * pye + pze * pze);
		// pze = wo - pz1 - pz2;
		// Ep = wo - w1 - w2 + 1.0; pp = 1.0/Math.sqrt(pxe*pxe + pye*pye +
		// pze*pze);

		// ******************************
		// Set up particles on the stack and rotate back to the lab frame
		// ******************************/
		EGS4.NPold = EGS4.NP;// np;
		// $TRANSFER PROPERTIES TO (np+1) FROM (np);
		EGS4.X[EGS4.NP] = EGS4.X[EGS4.NP - 1];
		EGS4.Y[EGS4.NP] = EGS4.Y[EGS4.NP - 1];
		EGS4.Z[EGS4.NP] = EGS4.Z[EGS4.NP - 1];
		EGS4.IR[EGS4.NP] = EGS4.IR[EGS4.NP - 1];
		EGS4.WT[EGS4.NP] = EGS4.WT[EGS4.NP - 1];
		EGS4.DNEAR[EGS4.NP] = EGS4.DNEAR[EGS4.NP - 1];
		EGS4.LATCH[EGS4.NP] = EGS4.LATCH[EGS4.NP - 1];

		// $TRANSFER PROPERTIES TO (np+2) FROM (np);
		EGS4.X[EGS4.NP + 1] = EGS4.X[EGS4.NP - 1];
		EGS4.Y[EGS4.NP + 1] = EGS4.Y[EGS4.NP - 1];
		EGS4.Z[EGS4.NP + 1] = EGS4.Z[EGS4.NP - 1];
		EGS4.IR[EGS4.NP + 1] = EGS4.IR[EGS4.NP - 1];
		EGS4.WT[EGS4.NP + 1] = EGS4.WT[EGS4.NP - 1];
		EGS4.DNEAR[EGS4.NP + 1] = EGS4.DNEAR[EGS4.NP - 1];
		EGS4.LATCH[EGS4.NP + 1] = EGS4.LATCH[EGS4.NP - 1];

		a = EGS4.U[EGS4.NP - 1];
		b = EGS4.V[EGS4.NP - 1];
		c = EGS4.W[EGS4.NP - 1];
		sinpsi = a * a + b * b;
		if (sinpsi > 1.e-20) {
			sinpsi = Math.sqrt(sinpsi);
			sindel = b / sinpsi;
			cosdel = a / sinpsi;
			EGS4.U[EGS4.NP - 1] = (c * cosdel * px2 - sindel * py2 + a * pz2)
					/ w2;
			EGS4.V[EGS4.NP - 1] = (c * sindel * px2 + cosdel * py2 + b * pz2)
					/ w2;
			EGS4.W[EGS4.NP - 1] = (c * pz2 - sinpsi * px2) / w2;
			EGS4.IQ[EGS4.NP - 1] = 0;
			EGS4.E[EGS4.NP - 1] = w2 * EGS4.PRM;
			EGS4.NP = EGS4.NP + 1;
			EGS4.U[EGS4.NP - 1] = (c * cosdel * px1 - sindel * py1 + a * pz1)
					/ w1;
			EGS4.V[EGS4.NP - 1] = (c * sindel * px1 + cosdel * py1 + b * pz1)
					/ w1;
			EGS4.W[EGS4.NP - 1] = (c * pz1 - sinpsi * px1) / w1;
			EGS4.IQ[EGS4.NP - 1] = 0;
			EGS4.E[EGS4.NP - 1] = w1 * EGS4.PRM;
			EGS4.NP = EGS4.NP + 1;
			EGS4.U[EGS4.NP - 1] = (c * cosdel * pxe - sindel * pye + a * pze)
					* pp;
			EGS4.V[EGS4.NP - 1] = (c * sindel * pxe + cosdel * pye + b * pze)
					* pp;
			EGS4.W[EGS4.NP - 1] = (c * pze - sinpsi * pxe) * pp;
			EGS4.IQ[EGS4.NP - 1] = -1;
			EGS4.E[EGS4.NP - 1] = Ep * EGS4.PRM;
		} else {
			EGS4.U[EGS4.NP - 1] = px2 / w2;
			EGS4.V[EGS4.NP - 1] = py2 / w2;
			EGS4.W[EGS4.NP - 1] = c * pz2 / w2;
			EGS4.IQ[EGS4.NP - 1] = 0;
			EGS4.E[EGS4.NP - 1] = w2 * EGS4.PRM;
			EGS4.NP = EGS4.NP + 1;
			EGS4.U[EGS4.NP - 1] = px1 / w1;
			EGS4.V[EGS4.NP - 1] = py1 / w1;
			EGS4.W[EGS4.NP - 1] = c * pz1 / w1;
			EGS4.IQ[EGS4.NP - 1] = 0;
			EGS4.E[EGS4.NP - 1] = w1 * EGS4.PRM;
			EGS4.NP = EGS4.NP + 1;
			EGS4.U[EGS4.NP - 1] = pxe * pp;
			EGS4.V[EGS4.NP - 1] = pye * pp;
			EGS4.W[EGS4.NP - 1] = c * pze * pp;
			EGS4.IQ[EGS4.NP - 1] = -1;
			EGS4.E[EGS4.NP - 1] = Ep * EGS4.PRM;
		}

	}

	// #############################################################################################
	// "***************************************************************************"
	// "                                                                           "
	// " Sampling of triplet production events.                                    "
	// "                                                                           "
	// " The treatment is based on Borsellino's first Born approximation           "
	// " result (see Eq. 4B-3002 in the pair article of Motz, Olsen & Koch)        "
	// " As the kinematic of the process is already complicated enough and the     "
	// " cross section itself is not simple either, a Markov-chain method is used  "
	// " to sample triplet events from the Borsellino equation without any         "
	// " additional approximations (other then the use of the first Born           "
	// " approximation and the assumption of free electrons implied by             "
	// " Borsellino's derivation)                                                  "
	// "                                                                           "
	// " Iwan Kawrakow, April 2005.                                                "
	// "***************************************************************************"
	/**
	 * Called by PAIR. Handle triple production. The treatment is based on Borsellino's first Born approximation 
	 * result (see Eq. 4B-3002 in the pair article of Motz, Olsen and Koch). As the kinematic of the process is already complicated enough and the 
	 * cross section itself is not simple either, a Markov-chain method is used to sample triplet events from the Borsellino equation without any 
	 * additional approximations (other then the use of the first Born approximation and the assumption of free electrons implied by Borsellino's derivation).
	 */
	protected static void sample_triplet() {

		// implicit none;
		// $declare_max_medium;
		// ;COMIN/EPCONT,STACK,MEDIA,THRESH,USEFUL,RANDOM,USER/;
		// $declare_write_buffer;

		// " We use double precision throughout as in many cases the kinematically "
		// " permitted angular interval is too small to be resolved accurately enough "
		// " in single precision "

		// real*8 fmax_array($MAX_TRIPLET), eta_p_array($MAX_TRIPLET),
		// eta_Ep_array($MAX_TRIPLET), eta_costp_array($MAX_TRIPLET),
		// eta_costm_array($MAX_TRIPLET), ebin_array($MAX_TRIPLET),
		// wp_array($MAX_TRIPLET), qmin_array($MAX_TRIPLET);
		// real*8 kmin, kmax, dlogki, alogkm, prmi, tiny_eta;

		double ai = 0.0;
		double rnno = 0.0;
		double k = 0.0;
		double qmin = 0.0;
		double qmax = 0.0;
		double aux = 0.0;
		double a1 = 0.0;
		double a2 = 0.0;
		double a3 = 0.0;
		double D = 0.0;
		double px1 = 0.0;
		double px2 = 0.0;
		double pp_min = 0.0;
		double pp_max = 0.0;
		double Ep_min = 0.0;
		double Ep_max = 0.0;
		double k2p2 = 0.0;
		double k2p2x = 0.0;
		double peig = 0.0;
		double b = 0.0;
		double aux1 = 0.0;
		double aux12 = 0.0;
		double D1 = 0.0;
		double aux3 = 0.0;
		double xmin = 0.0;
		double xmax = 0.0;
		double aux6 = 0.0;
		double aux7 = 0.0;
		double uu = 0.0;
		double cphi = 0.0;
		double sphi = 0.0;
		double cphi_factor = 0.0;
		double aux5 = 0.0;
		double phi = 0.0;
		double tmp = 0.0;
		double Er = 0.0;
		double pr = 0.0;
		double pr2 = 0.0;
		double eta_pr = 0.0;
		double Ep = 0.0;
		double pp = 0.0;
		double pp2 = 0.0;
		double wEp = 0.0;
		double cost_p = 0.0;
		double sint_p = 0.0;
		double eta_Ep = 0.0;
		double mup_min = 0.0;
		double wmup = 0.0;
		double eta_costp = 0.0;
		double Epp = 0.0;
		double pp_sintp = 0.0;
		double pp_sntp2 = 0.0;
		double Em = 0.0;
		double pm = 0.0;
		double pm2 = 0.0;
		double cost_m = 0.0;
		double sint_m = 0.0;
		//double Emm = 0.0;
		double wmum = 0.0;
		double pm_sintm = 0.0;
		double eta_costm = 0.0;
		double k2 = 0.0;
		double k3 = 0.0;
		double s2 = 0.0;
		double s3 = 0.0;
		double k2k3i = 0.0;
		double k22 = 0.0;
		double k32 = 0.0;
		double q2 = 0.0;
		double aux4 = 0.0;
		double S_1 = 0.0;
		double S_2 = 0.0;
		double sigma = 0.0;
		double ppx = 0.0;
		double ppy = 0.0;
		double ppz = 0.0;
		double pmx = 0.0;
		double pmy = 0.0;
		double pmz = 0.0;
		double prx = 0.0;
		double pry = 0.0;
		double prz = 0.0;
		double a = 0.0;
		double c = 0.0;
		double sindel = 0.0;
		double cosdel = 0.0;
		double sinpsi = 0.0;

		int i = 0;
		boolean use_it = false;
		//int iscore = 0; // " needed for BEAM "

		// $LOGICAL is_initialized;
		// data is_initialized/.false./;
		// save
		// is_initialized,fmax_array,eta_p_array,eta_Ep_array,eta_costp_array,
		// eta_costm_array,ebin_array,wp_array,qmin_array,
		// kmin,kmax,dlogki,alogkm,prmi,tiny_eta;
		boolean retry_triplet = false;

		if (!is_initialized) {
			is_initialized = true;
			tiny_eta = 1.e-6;
			// " Set current cross section value to -1 in each energy bin "
			// DO i=1,$MAX_TRIPLET [ fmax_array(i) = -1; ]
			for (i = 1; i <= EGS4.$MAX_TRIPLET; i++) {
				fmax_array[i - 1] = -1;
			}
			// " Find the maximum energy of the cross section data "
			kmax = 0;
			kmin = 4.1 * EGS4.PRM;// prm;
			// DO i=1,nmed [ IF( up(i) > kmax ) kmax = UP(i); ]
			for (i = 1; i <= EGS4.NMED; i++) {
				if (EGS4.UP[i - 1] > kmax)
					kmax = EGS4.UP[i - 1];
			}
			if (kmax <= kmin)
				return;
			dlogki = EGS4.$MAX_TRIPLET - 1.0;
			dlogki = dlogki / Math.log(kmax / kmin);
			alogkm = 1.0 - dlogki * Math.log(kmin);
			prmi = 1.0 / EGS4.PRM;
			for (i = 1; i <= EGS4.$MAX_TRIPLET; i++)// DO i=1,$MAX_TRIPLET [
			{
				k = 4.1 * Math.exp((i - 1.) / dlogki);
				ebin_array[i - 1] = k;
				qmin = 4.0
						* k
						/ (k * (k - 1.0) + (k + 1.0) * Math.sqrt(k * (k - 4.0)));
				qmax = (k * (k - 1.0) + (k + 1.0) * Math.sqrt(k * (k - 4.0)))
						/ (2.0 * k + 1.0);
				qmin_array[i - 1] = qmin;
				wp_array[i - 1] = Math.log(qmax / qmin);
			}
		}

		peig = EGS4.E[EGS4.NP - 1];
		if (peig <= 4.0 * EGS4.PRM)
			return;
		// $CHECK-STACK(np+2,'sample_triplet');
		if (EGS4.NP + 2 > EGS4.$MXSTACK) {
			EGS4.STOPPROGRAM = true;

			EGS4.seqStr = " ***************************************************"
					+ "  \n"
					+ " In subroutine "
					+ "sample_triplet"
					+ " stack size exceeded!"
					+ "  \n"
					+ " $MXSTACK = "
					+ EGS4.$MXSTACK
					+ " np = "
					+ EGS4.NP
					+ "  \n"
					+ " Increase $MXSTACK and try again "
					+ "  \n"
					+ " Terminating execution "
					+ "  \n"
					+ " ***************************************************";// +"  \n";
			// if(EGS4.iprint>2)
			eq.printSequence(EGS4.seqStr);

			return;// stop;
		}

		// " Determine energy bin "
		if (peig <= kmin) {
			i = 1;
		} else if (peig >= kmax) {
			i = EGS4.$MAX_TRIPLET;
		} else {
			ai = alogkm + dlogki * EGS4.GLE;
			Double aid = new Double(ai);
			i = aid.intValue();
			ai = ai - i;
			rnno = EGS4.random01();
			if (rnno < ai) {
				i = i + 1;
			}
		}

		// " First use the bin energy to sample the random numbers "
		// " that determine recoil momentum and electron/postron angles "
		k = ebin_array[i - 1];// (i);

		// In the following: k is incident photon energy in units of m*c^2
		// (all energies are in units of m*c^2, momenta in
		// units of m*c)
		// Er,pr is energy, momentum of the recoil electron
		// Ep,pp is energy, momentum of the pair positron
		// Em,pm is energy, momentum of the pair electron
		// cost_p, sint_p is cos, sin of the positron angle
		// with respect to k
		// cost_m, sint_m same but for the electron
		// cphi is cos of azimuthal angle between positron
		// and pair electron directions.

		// :retry_triplet:;
		while (true) {
			retry_triplet = false;
			// " Pick the recoil electron momentum from 1/p.
			eta_pr = EGS4.random01();
			if (eta_pr < tiny_eta)
				eta_pr = tiny_eta;
			pr = qmin_array[i - 1] * Math.exp(eta_pr * wp_array[i - 1]);
			pr2 = pr * pr;
			Er = Math.sqrt(1.0 + pr2);

			// " Determine min./max. kinematically permitted postron energy for "
			// " this k and p "
			aux = Er - pr - 1.0;
			a1 = (k - pr) * (1.0 - Er - k * aux);
			a2 = 1.0 + k - Er;
			a3 = 1.0 / (aux * (pr + Er - 2.0 * k - 1.0));
			D = a2
					* Math.sqrt(aux
							* (2.0 * k * Er + k * k * aux - pr
									* (Er + pr + 1.0) / 2.0));
			px1 = (a1 + D) * a3;
			px2 = (a1 - D) * a3;
			if (px1 < px2) {
				pp_min = px1;
				pp_max = px2;
			} else {
				pp_min = px2;
				pp_max = px1;
			}
			Ep_min = Math.sqrt(1.0 + pp_min * pp_min);
			Ep_max = Math.sqrt(1.0 + pp_max * pp_max);

			// " Pick the positron energy "
			eta_Ep = EGS4.random01();
			if (eta_Ep < tiny_eta)
				eta_Ep = tiny_eta;
			wEp = Ep_max - Ep_min;
			Ep = Ep_min + eta_Ep * wEp;
			pp2 = Ep * Ep - 1.0;
			pp = Math.sqrt(pp2);
			k2p2 = k * k + pp2;

			// " Now we can determine the pair electron energy from energy conservation "
			Em = k + 1.0 - Er - Ep;
			pm2 = Em * Em - 1.0;
			pm = Math.sqrt(pm2);

			// " The minimum cosine of the positron angle follows from the kinematics. "
			mup_min = (k2p2 - (pr + pm) * (pr + pm)) / (2.0 * k * pp);

			// " Now pick the positron direction from 1/(Ep-pp*cost_p) "
			eta_costp = EGS4.random01();
			if (eta_costp < tiny_eta)
				eta_costp = tiny_eta;
			Epp = Ep / pp;
			wmup = Math.log((Epp - 1.0) / (Epp - mup_min));
			cost_p = Epp - (Epp - mup_min) * Math.exp(wmup * eta_costp);
			wmup = wmup * (cost_p - Epp);
			sint_p = 1.0 - cost_p * cost_p;
			if (sint_p > 1.e-20) {
				sint_p = Math.sqrt(sint_p);
			} else {
				sint_p = 1.e-10;
			}
			k2p2x = k2p2 - 2.0 * k * pp * cost_p;

			// " The minimum amd maximum cosine of the pair electron angle follows from "
			// " the kinematics "
			b = pr2 - k2p2x - pm2;
			aux1 = k - pp * cost_p;
			aux12 = aux1 * aux1;
			pp_sintp = pp * sint_p;
			pp_sntp2 = pp_sintp * pp_sintp;
			D1 = pm2 * (aux12 + pp_sntp2) - b * b / 4.0;
			if (D1 <= 0) {
				// goto :retry_triplet:;
				retry_triplet = true;
			}
			if (!retry_triplet)// 1
			{
				D = 2.0 * pp_sintp * Math.sqrt(D1);
				aux3 = 0.5 / (aux12 + pp_sntp2);
				xmin = (-b * aux1 - D) * aux3;
				xmax = (-b * aux1 + D) * aux3;

				// " Now pick the electron direction from "
				// "  1/(Em-pm*cost_m)/sqrt((cost_m_max-cost_m)*(cost_m-cost_m_min)) "
				// " We have to take into account the "
				// " 1/sqrt((cost_m_max-cost_m)*(cost_m-cost_m_min)) factor in the sampling "
				// " otherwise we end up with 1/sqrt() singularities near the ends of the "
				// " allowed cost_m range                                                 "
				eta_costm = EGS4.random01();
				if (eta_costm < tiny_eta)
					eta_costm = tiny_eta;
				aux6 = Math.sqrt((Em - xmin) / (Em - xmax));
				aux7 = aux6 * Math.tan(1.570796326794897 * eta_costm);
				uu = (aux7 - 1.0) / (aux7 + 1.0);
				cost_m = 0.5 * (xmax + xmin + 2.0 * uu * (xmax - xmin)
						/ (1.0 + uu * uu));
				wmum = Math.sqrt((xmax - cost_m) * (cost_m - xmin));
				wmum = wmum * aux6 * (Em - cost_m) / (Em - xmin);
				cost_m = cost_m / pm;
				sint_m = Math.sqrt(1.0 - cost_m * cost_m);
				pm_sintm = pm * sint_m;

				// " Now we have selected all independent kinematic variables. "
				// " Determine the azimuthal angle between the pair electrons "
				cphi = (b + 2 * pm * cost_m * aux1) / (2 * pp_sintp * pm_sintm);
				if (Math.abs(cphi) >= 1) {
					// goto :retry_triplet:;
					retry_triplet = true;
				}
				if (!retry_triplet)// 2
				{

					sphi = Math.sqrt(1.0 - cphi * cphi);

					// " And now evaluate the Borsellino cross section "
					k3 = k * (pp * cost_p - Ep);
					k2 = k * (pm * cost_m - Em);
					k22 = k2 * k2;
					k32 = k3 * k3;
					k2k3i = 1.0 / (k2 * k3);
					s2 = pp * pm * (cost_p * cost_m + sint_p * sint_m * cphi)
							- Ep * Em;
					s3 = k2 - Em + 1.0 - s2;
					q2 = 2.0 * (Er - 1.0);
					S_1 = k32 + k22 + (q2 - 2.0) * s2 - (1.0 - q2 / 2)
							* (k32 + k22) * k2k3i;
					aux4 = k3 * Ep - k2 * Em;
					S_2 = -q2 * (Ep * Ep + Em * Em) + 2.0 * s2
							- (2.0 * aux4 * aux4 - k22 - k32) * k2k3i;
					sigma = Math
							.abs(pp
									* pm2
									* pm
									* k2k3i
									/ (q2 * q2 * (Em * s3 + Er))
									* (S_1 * (1.0 - q2 / 4.0) + S_2
											* (1.0 + q2 / 4.0)));

					// " We get the following factor due to the transformation from phi to "
					// " the recoil momentum pr "
					cphi_factor = Math.abs(2.0 * Er * pm2 - Em
							* (k2p2x - pr2 - pm2))
							/ (2.0 * pp_sintp * pm_sintm * pm2 * sphi);

					// " We have to also multiply by the various factors from the sampling of "
					// " pr, Ep, cost_p and cost_m "
					sigma = sigma * cphi_factor * wEp * wmup * wmum
							* wp_array[i - 1] * pr2 / Er;
					if (sigma < 0) {
						// $egs_warning(*,'In triplet sigma < 0 ? ',sigma);
						EGS4.seqStr = " In triplet sigma < 0 ?  " + sigma;// +" \n";
						if (EGS4.iprint > 2)
							eq.printSequence(EGS4.seqStr);

					}

					// " Now determine if we accept this new event "
					use_it = true;
					if (sigma < fmax_array[i - 1]) {
						rnno = EGS4.random01();
						if (sigma < fmax_array[i - 1] * rnno) {
							use_it = false;
						}
					}
					if (use_it) { // " Yes, event accepted "
						fmax_array[i - 1] = sigma;
						eta_p_array[i - 1] = eta_pr;
						eta_Ep_array[i - 1] = eta_Ep;
						eta_costp_array[i - 1] = eta_costp;
						eta_costm_array[i - 1] = eta_costm;
					} else { // " Nop, event rejected => use last accepted "
						eta_pr = eta_p_array[i - 1];
						eta_Ep = eta_Ep_array[i - 1];
						eta_costp = eta_costp_array[i - 1];
						eta_costm = eta_costm_array[i - 1];
					}

					// " We now have a set of random number accepted for sampling around "
					// " the i'th bin energy. We need to recalculate all variables using "
					// " the actual photon energy "

					k = peig * prmi;
					aux5 = k * (k - 1.0) + (k + 1.0) * Math.sqrt(k * (k - 4.0));
					qmin = 4.0 * k / aux5;
					qmax = aux5 / (2.0 * k + 1.0);
					pr = qmin * Math.exp(eta_pr * Math.log(qmax / qmin));
					pr2 = pr * pr;
					Er = Math.sqrt(1.0 + pr2);

					aux = Er - pr - 1.0;
					a1 = (k - pr) * (1.0 - Er - k * aux);
					a2 = 1.0 + k - Er;
					a3 = 1.0 / (aux * (pr + Er - 2.0 * k - 1.0));
					D = a2
							* Math.sqrt(aux
									* (2.0 * k * Er + k * k * aux - pr
											* (Er + pr + 1.0) / 2.0));
					px1 = (a1 + D) * a3;
					px2 = (a1 - D) * a3;
					if (px1 < px2) {
						pp_min = px1;
						pp_max = px2;
					} else {
						pp_min = px2;
						pp_max = px1;
					}
					Ep_min = Math.sqrt(1.0 + pp_min * pp_min);
					Ep_max = Math.sqrt(1.0 + pp_max * pp_max);

					wEp = Ep_max - Ep_min;
					Ep = Ep_min + eta_Ep * wEp;
					pp2 = Ep * Ep - 1.0;
					pp = Math.sqrt(pp2);
					k2p2 = k * k + pp2;
					Em = k + 1.0 - Er - Ep;
					pm2 = Em * Em - 1;
					pm = Math.sqrt(pm2);

					mup_min = (k2p2 - (pr + pm) * (pr + pm)) / (2.0 * k * pp);
					Epp = Ep / pp;
					wmup = Math.log((Epp - 1.0) / (Epp - mup_min));
					cost_p = Epp - (Epp - mup_min) * Math.exp(wmup * eta_costp);
					sint_p = Math.sqrt(1.0 - cost_p * cost_p);
					k2p2x = k2p2 - 2.0 * k * pp * cost_p;

					b = pr2 - k2p2x - pm2;
					aux1 = k - pp * cost_p;
					aux12 = aux1 * aux1;
					pp_sintp = pp * sint_p;
					pp_sntp2 = pp_sintp * pp_sintp;
					D1 = pm2 * (aux12 + pp_sntp2) - b * b / 4.0;
					if (D1 <= 0) {
						// goto :retry_triplet:;
						retry_triplet = true;
					}
					if (!retry_triplet)// 3
					{
						D = 2.0 * pp_sintp * Math.sqrt(D1);
						aux3 = 0.5 / (aux12 + pp_sntp2);
						xmin = (-b * aux1 - D) * aux3;
						xmax = (-b * aux1 + D) * aux3;
						aux6 = Math.sqrt((Em - xmin) / (Em - xmax));
						aux7 = aux6 * Math.tan(1.570796326794897 * eta_costm);
						uu = (aux7 - 1.0) / (aux7 + 1.0);
						cost_m = 0.5
								* (xmax + xmin + 2.0 * uu * (xmax - xmin)
										/ (1.0 + uu * uu)) / pm;
						sint_m = Math.sqrt(1.0 - cost_m * cost_m);
						pm_sintm = pm * sint_m;

						cphi = (b + 2.0 * pm * cost_m * aux1)
								/ (2.0 * pp_sintp * pm_sintm);
						if (Math.abs(cphi) >= 1) {
							// goto :retry_triplet:;
							retry_triplet = true;
						}
						if (!retry_triplet)// 4
						{
							sphi = Math.sqrt(1.0 - cphi * cphi);
							// @@@@@@@@@@@@@@@@@@@@@@
							// if here is ok,so:
							break;
							// @@@@@@@@@@@@@@@@@@@@@@
						}// if(!retry_triplet)//4
					}// if(!retry_triplet)//3
				}// if(!retry_triplet)//2
			}// if(!retry_triplet)//1
		}// while

		// OK, now the final momenta are
		// Positron: (pp*sint_p, 0, pp*cost_p)
		// Electron: (pm*sint_m*cphi,pm*sint_m*sphi, pm*cost_m)
		// Recoil electron: k - pp - pm
		// This is in a frame where the photon is moving along the z axis.
		// We have to pick another azimuthal angle randomly, rotate the
		// x- and y- components of pp and pm by that, determine the recoil
		// momentum from momentum conservation and then rotate all three
		// momenta back into the lab frame.

		phi = EGS4.random01();
		phi = phi * 6.283185307179586;
		ppx = pp * sint_p;
		ppy = 0.0;
		pmx = pm * sint_m * cphi;
		pmy = pm * sint_m * sphi;
		cphi = Math.cos(phi);
		sphi = Math.sin(phi);
		tmp = ppx * sphi;
		ppx = ppx * cphi - ppy * sphi;
		ppy = tmp + ppy * cphi;
		tmp = pmx * sphi;
		pmx = pmx * cphi - pmy * sphi;
		pmy = tmp + pmy * cphi;
		ppz = pp * cost_p;
		pmz = pm * cost_m;
		prx = -ppx - pmx;
		pry = -ppy - pmy;
		prz = k - ppz - pmz;
		// " Set up particles on the stack ";
		// " We always put the recoil electron on top (even if its energy is higher "
		// " then the energies of the pair particles) because                       "
		// "   - that way, we know which particle is the recoil  electron in case   "
		// "     we want to score some quantity related to it                       "
		// "   - its energy is, on average, lower than the pair particle energies   "

		EGS4.NPold = EGS4.NP;// np;
		// $TRANSFER PROPERTIES TO (np) FROM (np);?????????????????
		// $TRANSFER PROPERTIES TO (np+1) FROM (np);
		EGS4.X[EGS4.NP] = EGS4.X[EGS4.NP - 1];
		EGS4.Y[EGS4.NP] = EGS4.Y[EGS4.NP - 1];
		EGS4.Z[EGS4.NP] = EGS4.Z[EGS4.NP - 1];
		EGS4.IR[EGS4.NP] = EGS4.IR[EGS4.NP - 1];
		EGS4.WT[EGS4.NP] = EGS4.WT[EGS4.NP - 1];
		EGS4.DNEAR[EGS4.NP] = EGS4.DNEAR[EGS4.NP - 1];
		EGS4.LATCH[EGS4.NP] = EGS4.LATCH[EGS4.NP - 1];

		// $TRANSFER PROPERTIES TO (np+2) FROM (np); replace with
		// $TRANSFER PROPERTIES TO (np+2) FROM (np+1);
		EGS4.X[EGS4.NP + 1] = EGS4.X[EGS4.NP];
		EGS4.Y[EGS4.NP + 1] = EGS4.Y[EGS4.NP];
		EGS4.Z[EGS4.NP + 1] = EGS4.Z[EGS4.NP];
		EGS4.IR[EGS4.NP + 1] = EGS4.IR[EGS4.NP];
		EGS4.WT[EGS4.NP + 1] = EGS4.WT[EGS4.NP];
		EGS4.DNEAR[EGS4.NP + 1] = EGS4.DNEAR[EGS4.NP];
		EGS4.LATCH[EGS4.NP + 1] = EGS4.LATCH[EGS4.NP];

		pp = 1.0 / pp;
		pm = 1.0 / pm;
		pr = 1.0 / pr;
		a = EGS4.U[EGS4.NP - 1];
		b = EGS4.V[EGS4.NP - 1];
		c = EGS4.W[EGS4.NP - 1];
		sinpsi = a * a + b * b;
		if (sinpsi > 1.e-20) {
			sinpsi = Math.sqrt(sinpsi);
			sindel = b / sinpsi;
			cosdel = a / sinpsi;
			if (Ep > Em) {
				EGS4.U[EGS4.NP - 1] = pp
						* (c * cosdel * ppx - sindel * ppy + a * ppz);
				EGS4.V[EGS4.NP - 1] = pp
						* (c * sindel * ppx + cosdel * ppy + b * ppz);
				EGS4.W[EGS4.NP - 1] = pp * (c * ppz - sinpsi * ppx);
				EGS4.IQ[EGS4.NP - 1] = 1;
				EGS4.E[EGS4.NP - 1] = Ep * EGS4.PRM;// prm;
				EGS4.U[EGS4.NP] = pm
						* (c * cosdel * pmx - sindel * pmy + a * pmz);
				EGS4.V[EGS4.NP] = pm
						* (c * sindel * pmx + cosdel * pmy + b * pmz);
				EGS4.W[EGS4.NP] = pm * (c * pmz - sinpsi * pmx);
				EGS4.IQ[EGS4.NP] = -1;
				EGS4.E[EGS4.NP] = Em * EGS4.PRM;
			} else {
				EGS4.U[EGS4.NP] = pp
						* (c * cosdel * ppx - sindel * ppy + a * ppz);
				EGS4.V[EGS4.NP] = pp
						* (c * sindel * ppx + cosdel * ppy + b * ppz);
				EGS4.W[EGS4.NP] = pp * (c * ppz - sinpsi * ppx);
				EGS4.IQ[EGS4.NP] = 1;
				EGS4.E[EGS4.NP] = Ep * EGS4.PRM;
				EGS4.U[EGS4.NP - 1] = pm
						* (c * cosdel * pmx - sindel * pmy + a * pmz);
				EGS4.V[EGS4.NP - 1] = pm
						* (c * sindel * pmx + cosdel * pmy + b * pmz);
				EGS4.W[EGS4.NP - 1] = pm * (c * pmz - sinpsi * pmx);
				EGS4.IQ[EGS4.NP - 1] = -1;
				EGS4.E[EGS4.NP - 1] = Em * EGS4.PRM;
			}
			EGS4.NP = EGS4.NP + 2;
			EGS4.U[EGS4.NP - 1] = pr
					* (c * cosdel * prx - sindel * pry + a * prz);
			EGS4.V[EGS4.NP - 1] = pr
					* (c * sindel * prx + cosdel * pry + b * prz);
			EGS4.W[EGS4.NP - 1] = pr * (c * prz - sinpsi * prx);
			EGS4.IQ[EGS4.NP - 1] = -1;
			EGS4.E[EGS4.NP - 1] = Er * EGS4.PRM;
		} else {
			if (Ep > Em) {
				EGS4.U[EGS4.NP - 1] = pp * ppx;
				EGS4.V[EGS4.NP - 1] = pp * ppy;
				EGS4.W[EGS4.NP - 1] = c * pp * ppz;
				EGS4.IQ[EGS4.NP - 1] = 1;
				EGS4.E[EGS4.NP - 1] = Ep * EGS4.PRM;
				EGS4.U[EGS4.NP] = pm * pmx;
				EGS4.V[EGS4.NP] = pm * pmy;
				EGS4.W[EGS4.NP] = c * pm * pmz;
				EGS4.IQ[EGS4.NP] = -1;
				EGS4.E[EGS4.NP] = Em * EGS4.PRM;
			} else {
				EGS4.U[EGS4.NP] = pp * ppx;
				EGS4.V[EGS4.NP] = pp * ppy;
				EGS4.W[EGS4.NP] = c * pp * ppz;
				EGS4.IQ[EGS4.NP] = 1;
				EGS4.E[EGS4.NP] = Ep * EGS4.PRM;
				EGS4.U[EGS4.NP - 1] = pm * pmx;
				EGS4.V[EGS4.NP - 1] = pm * pmy;
				EGS4.W[EGS4.NP - 1] = c * pm * pmz;
				EGS4.IQ[EGS4.NP - 1] = -1;
				EGS4.E[EGS4.NP - 1] = Em * EGS4.PRM;
			}
			EGS4.NP = EGS4.NP + 2;
			EGS4.U[EGS4.NP - 1] = pr * prx;
			EGS4.V[EGS4.NP - 1] = pr * pry;
			EGS4.W[EGS4.NP - 1] = c * pr * prz;
			EGS4.IQ[EGS4.NP - 1] = -1;
			EGS4.E[EGS4.NP - 1] = Er * EGS4.PRM;
		}

		return;
	}

	// see EGS4 reset!!!
	/**
	 * Reset global variables for re-use.
	 */
	public static void reset() {
		IRCODE = 0;
		uscat = 0.0;
		vscat = 0.0;
		wscat = 0.0;
		xtrans = 0.0;
		ytrans = 0.0;
		ztrans = 0.0;
		spin_index = false;
		find_index = false;
		sscatcallindex = 0;
		mscatcallindex = 0;
		w1 = 0.0;
		sint1 = 0.0;
		w2 = 0.0;
		sint2 = 0.0;
		ws = 0.0;
		sint = 0.0;
		msdist1call = 0;
		msdist2call = 0;
		ierust = 0;
		CTHET = 0.0;
		//PHI = 0.0;
		//CPHI = 0.0;
		A = 0.0;
		B = 0.0;
		C = 0.0;
		SINPS2 = 0.0;
		SINPSI = 0.0;
		US = 0.0;
		VS = 0.0;
		SINDEL = 0.0;
		COSDEL = 0.0;
		i = 0;
		j = 0;
		omega2 = 0.0;
		isrj = 0;
		jsrj = 0;
		n_warning = 0;
		is_initialized = false;
		fmax_array = new double[EGS4.$MAX_TRIPLET];
		eta_p_array = new double[EGS4.$MAX_TRIPLET];
		eta_Ep_array = new double[EGS4.$MAX_TRIPLET];
		eta_costp_array = new double[EGS4.$MAX_TRIPLET];
		eta_costm_array = new double[EGS4.$MAX_TRIPLET];
		ebin_array = new double[EGS4.$MAX_TRIPLET];
		wp_array = new double[EGS4.$MAX_TRIPLET];
		qmin_array = new double[EGS4.$MAX_TRIPLET];
		kmin = 0.0;
		kmax = 0.0;
		dlogki = 0.0;
		alogkm = 0.0;
		prmi = 0.0;
		tiny_eta = 0.0;

		// IQI=0;EI=0.0;XI=0.0;YI=0.0;ZI=0.0;UI=0.0;VI=0.0;WI=0.0;IRI=0;WTI=0.0;

	}

	// "============================================================================"
	/**
	 * Curently not used.
	 * @param imed imed
	 * @param fac fac
	 * @param which which
	 */
	public static void egs_scale_photon_xsection(int imed, double fac, int which) {
		// "============================================================================"
		// "
		// " Scale the photon cross section 'which' by factor fac for medium
		// imed.
		// " which = 0 for all cross sections
		// " = 1 for Rayleigh scattering
		// " = 2 for Compton scattering
		// " = 3 for Pair production
		// " = 4 for photo-absorption
		// " If imed = 0, scaling is done for all media.
		// "============================================================================"
		// implicit none;
		// $INTEGER imed,which;
		// $REAL fac;

		// $declare_max_medium;//REPLACE {$declare_max_medium;} WITH {;};
		// ;COMIN/MEDIA,PHOTIN/;

		int ifirst = 0;
		int ilast = 0;
		int medium = 0;
		int j = 0;
		boolean has_r = false;
		double gle = 0.0;
		double gmfp = 0.0;
		double gbr1 = 0.0;
		double gbr2 = 0.0;
		double cohfac = 0.0;
		double aux = 0.0;
		double gmfp_old = 0.0;
		double gbr1_old = 0.0;
		double gbr2_old = 0.0;
		double cohfac_old = 0.0;
		// character*8 strings(5);
		String[] strings = { "photon", "Rayleigh", "Compton", "pair", "photo" };

		if (which < 0 || which > 4) {
			return;
		}
		if (imed > 0 && imed <= EGS4.NMED) {
			ifirst = imed;
			ilast = imed;
		} else {
			ifirst = 1;
			ilast = EGS4.NMED;
		}
		if (which == 1) {
			has_r = false;
			for (medium = ifirst; medium <= ilast; medium++) {
				if (EGS4.IRAYLM[medium - 1] == 1) {
					has_r = true;
				}
			}
			if (!has_r)
				return;
		}
		// $egs_info(*,' ');
		for (medium = ifirst; medium <= ilast; medium++) {

			// $egs_info('(a,a,a,i3,a,f9.5)',
			// 'Scaling ',strings(which+1),' x-section data for medium',
			// medium,' with ',fac);
			EGS4.seqStr = "Scaling " + strings[which]
					+ " x-section data for medium " + medium + " with "
					+ EGS4.format(fac, 9);
			if (EGS4.iprint > 1)
				eq.printSequence(EGS4.seqStr);

			for (j = 1; j <= EGS4.MGE[medium - 1]; j++) {
				gle = (j - EGS4.GE0[medium - 1]) / EGS4.GE1[medium - 1];
				gmfp = EGS4.GMFP0[j - 1][medium - 1]
						+ EGS4.GMFP1[j - 1][medium - 1] * gle;
				gbr1 = EGS4.GBR10[j - 1][medium - 1]
						+ EGS4.GBR11[j - 1][medium - 1] * gle;
				gbr2 = EGS4.GBR20[j - 1][medium - 1]
						+ EGS4.GBR21[j - 1][medium - 1] * gle;
				if (EGS4.IRAYLM[medium - 1] == 1) {
					cohfac = EGS4.COHE0[j - 1][medium - 1]
							+ EGS4.COHE1[j - 1][medium - 1] * gle;
				} else {
					cohfac = 1.;
				}
				if (which == 0) {
					gmfp = gmfp / fac;
				} else if (which == 1) {
					cohfac = cohfac / (fac * (1.0 - cohfac) + cohfac);
				} else {
					if (which == 2) {
						aux = fac * (gbr2 - gbr1) + gbr1 + 1.0 - gbr2;
						gbr2 = (gbr1 + fac * (gbr2 - gbr1)) / aux;
						gbr1 = gbr1 / aux;
					} else if (which == 3) {
						aux = fac * gbr1 + 1.0 - gbr1;
						gbr2 = (fac * gbr1 + gbr2 - gbr1) / aux;
						gbr1 = fac * gbr1 / aux;
					} else {
						aux = gbr2 + fac * (1.0 - gbr2);
						gbr1 = gbr1 / aux;
						gbr2 = gbr2 / aux;
					}
					gmfp = gmfp / aux;
					cohfac = cohfac * aux / (aux * cohfac + 1.0 - cohfac);
				}
				if (j > 1) {
					EGS4.GMFP1[j - 2][medium - 1] = (gmfp - gmfp_old)
							* EGS4.GE1[medium - 1];
					EGS4.GMFP0[j - 2][medium - 1] = gmfp
							- EGS4.GMFP1[j - 2][medium - 1] * gle;
					EGS4.GBR11[j - 2][medium - 1] = (gbr1 - gbr1_old)
							* EGS4.GE1[medium - 1];
					EGS4.GBR10[j - 2][medium - 1] = gbr1
							- EGS4.GBR11[j - 2][medium - 1] * gle;
					EGS4.GBR21[j - 2][medium - 1] = (gbr2 - gbr2_old)
							* EGS4.GE1[medium - 1];
					EGS4.GBR20[j - 2][medium - 1] = gbr2
							- EGS4.GBR21[j - 2][medium - 1] * gle;
					EGS4.COHE1[j - 2][medium - 1] = (cohfac - cohfac_old)
							* EGS4.GE1[medium - 1];
					EGS4.COHE0[j - 2][medium - 1] = cohfac
							- EGS4.COHE1[j - 2][medium - 1] * gle;
				}
				gmfp_old = gmfp;
				gbr1_old = gbr1;
				gbr2_old = gbr2;
				cohfac_old = cohfac;

			}
			EGS4.GMFP1[EGS4.MGE[medium - 1] - 1][medium - 1] = EGS4.GMFP1[EGS4.MGE[medium - 1] - 2][medium - 1];
			EGS4.GMFP0[EGS4.MGE[medium - 1] - 1][medium - 1] = EGS4.GMFP0[EGS4.MGE[medium - 1] - 2][medium - 1];
			EGS4.GBR11[EGS4.MGE[medium - 1] - 1][medium - 1] = EGS4.GBR11[EGS4.MGE[medium - 1] - 2][medium - 1];
			EGS4.GBR10[EGS4.MGE[medium - 1] - 1][medium - 1] = EGS4.GBR10[EGS4.MGE[medium - 1] - 2][medium - 1];
			EGS4.GBR21[EGS4.MGE[medium - 1] - 1][medium - 1] = EGS4.GBR21[EGS4.MGE[medium - 1] - 2][medium - 1];
			EGS4.GBR20[EGS4.MGE[medium - 1] - 1][medium - 1] = EGS4.GBR20[EGS4.MGE[medium - 1] - 2][medium - 1];
			EGS4.COHE1[EGS4.MGE[medium - 1] - 1][medium - 1] = EGS4.COHE1[EGS4.MGE[medium - 1] - 2][medium - 1];
			EGS4.COHE0[EGS4.MGE[medium - 1] - 1][medium - 1] = EGS4.COHE0[EGS4.MGE[medium - 1] - 2][medium - 1];
		}

		return; // end;
	}
}
