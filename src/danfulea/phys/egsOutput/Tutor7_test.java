package danfulea.phys.egsOutput;

import java.util.Calendar;
import java.util.Date;

import danfulea.phys.egs.EGS4;
import danfulea.phys.egs.EGS4Core;
import danfulea.phys.egs.EGS4Macro;
import danfulea.phys.egs.EgsQuestion;

/**
 * Demo class
 * Based on EGS4/EGSnrc routines developed by SLAC(Stanford Linear Accelerator Center)/NRC (National Research Council of Canada).
 * @author Dan Fulea, 17 OCT. 2005
 */
public class Tutor7_test implements EgsQuestion {
	// I DO NOT WANT TO READ INPUTS FROM FILE BUT DIRECTLY FROM
	// GUI!!!!!!!!!!!!!!!!!!!!!!!!

	// " An EGSnrc user code. It calculates reflected, deposited and         "
	// " transmitted energy for electron and photon beams incident on        "
	// " a slab geometry.                                                    "
	// "*******************************************************************************
	// " MC TRANSPORT PARAMETER
	// " **********************
	// "
	// " Global ECUT= Global (in all regions) electron transport cut
	// " off energy (in MeV). If this imput is missing,
	// " AE(medium) will be used.
	// " [ ECUT ]
	// " Global PCUT= Global (in all regions) photon transport cut
	// " off energy (in MeV). If this imput is missing,
	// " AP(medium) will be used.
	// " [ PCUT ]
	// " Global SMAX= Global (in all regions) maximum step-size
	// " restriction for electron transport (in cm).
	// " If missing, no geometrical step-size restrictions
	// " will be employed. Note that if you use the default
	// " EGSnrc electron-step algorithm, no SMAX-restriction
	// " is necessary. Option is useful for transport in low
	// " density materials (air) when PRESTA behaviour is
	// " turned on (see below)
	// " [ SMAXIR ]
	// " ESTEPE= Maximum fractional energy loss per step.
	// " Note that this is a global option only, no
	// " region-by-region setting is possible. If missing,
	// " the defualt is 0.25 (25%).
	// " [ ESTEPE ]
	// " XImax= Maximum first elastic scattering moment per step.
	// " Default is 0.5, NEVER use value greater than 1 as
	// " this is beyond the range of MS data available.
	// " [ XIMAX ]
	// " Boundary crossing algorithm=
	// " There are two selections possible: EXACT, means
	// " the algorithm will cross boundaries in a single
	// " scattering (SS) mode, the distance from a boundary
	// " at which the transition to SS mode is made is
	// " determined by 'Skin depth for BCA' (see below).
	// " The second option is PRESTA-I, if selected boundaries
	// " will be crossed a la PRESTA, i.e. with lateral
	// " correlations turned off and MS forced at boundaries.
	// " Default is EXACT.
	// " [ bca_algorithm, exact_bca ]
	// " Skin depth for BCA=
	// " Determines the distance from a boundary (in elastic
	// " MFP) at which the algorithm will go into single
	// " scattering mode (if EXACT boundary crossing) or
	// " swith off lateral correlations (if PRESTA-I boundary
	// " crossing). Default value is 3 for EXACT or
	// " exp(BLCMIN)/BLCMIN for PRESTA-I (see the PRESTA paper
	// " for a definition of BLCMIN). Note that if you choose
	// " EXACT boundary crossing and set Skin depth for BCA
	// " to a very large number (e.g. 1e10), the entire
	// " calculation will be in SS mode. If you choose
	// " PRESTA-I boundary crossing and make Skin depth for BCA
	// " large, you will get default EGS4 behavious (no PRESTA)
	// " [ skindepth_for_bca ]
	// " Electron-step algorithm=
	// " PRESTA-II (the default), the name is
	// " used for historical reasons
	// " or PRESTA-I
	// " Determines the algorithm used to take into account
	// " lateral and longitudinal correlations in a
	// " condensed history step.
	// " [ transport_algorithm ]
	// " Spin effects= Off, On, default is On
	// " Turns off/on spin effects for electron elastic
	// " scattering. Spin On is ABSOLUTELY necessary for
	// " good backscattering calculations. Will make a
	// " difference even in `well conditioned' situations
	// " (e.g. depth dose curves for RTP energy range
	// " electrons).
	// " [ spin_effects ]
	// " Brems angular sampling= Simple, KM, default is KM
	// " If Simple, use only the leading term of the Koch-Motz
	// " distribution to determine the emission angle of
	// " bremsstrahlung photons. If On, complete
	// " modified Koch-Motz 2BS is used (modifications
	// " concern proper handling of kinematics at low energies,
	// " makes 2BS almost the same as 2BN at low energies).
	// " [ IBRDST ]
	// " Brems cross sections= BH, NIST, default is BH
	// " If BH is selected, the Bethe-Heitler bremsstrahlung
	// " cross sections (Coulomb corrected above 50 MeV)
	// " will be used. If NIST is selected, the NIST brems
	// " cross section data base (which is the basis for
	// " the ICRU radiative stopping powers) will be employed.
	// " Differences are negligible for E > ,say, 10 MeV,
	// " but signifficant in the keV energy range.
	// " Bound Compton scattering= On or Off
	// " If Off, Compton scattering will be treated with
	// " Klein-Nishina, with On Compton scattering is
	// " treated in the Impulse approximation. Default is On
	// " Make sure to turn on for low energy applications,
	// " not necessary above, say, 1 MeV.
	// " [ IBCMP ]
	// " Pair angular sampling= Off, Simple or KM
	// " If off, pairs are set in motion at an angle m/E
	// " relative to the photon direction (m is electron rest
	// " energy, E the photon energy). Simple turns on
	// " the leading term of the angular distribution
	// " (this is sufficient for most applications),
	// " KM (comes from Koch and Motz) turns on using 2BS
	// " from the article by Koch and Motz.
	// " Default is Simple, make sure you always use
	// " Simple or KM
	// " [ IPRDST ]
	// " Photoelectron angular sampling= Off or On
	// " If Off, photo-electrons get the direction of the
	// " `mother' photon, with On, Sauter's furmula is
	// " used (which is, striktly speaking, valid only for
	// " K-shell photo-absorption).
	// " If the user has a better approach, replace the macro
	// " $SELECT-PHOTOELECTRON-DIRECTION;
	// " The only application that
	// " I encountered until now where this option made a
	// " small difference was a big ion chamber (cavity size
	// " comparable with electron range) with high-Z walls
	// " in a low energy photon beam.
	// " Default is On
	// " [ IPHTER ]
	// " Rayleigh scattering= Off, On
	// " If On, turnd on coherent (Rayleigh) scattering.
	// " Default is Off. Should be turned on for low energy
	// " applications. Not set to On by default because
	// " On requires a sperial PEGS4 data set
	// " [ IRAYLR ]
	// " Atomic relaxations= Off, On
	// " Default is On. The effect of using On is twofold:
	// " - In photo-electric absorption events, the element
	// " (if material is mixture) and the shell the photon
	// " is interacting with are sampled from the appropriate
	// " cross seections
	// " - Shell vacancies created in photo-absorption events
	// " are relaxed via emission of fluorescent X-Rays,
	// " Auger and Koster-Cronig electrons.
	// " Make sure to turn this option on for low energy
	// " applications.
	// " [ IEDGFL ]
	// "
	// " Atomic relaxations, Rayleigh scattering,
	// " Photoelectron angular sampling and Bound Compton scattering
	// " can also be turned On/Off on a region-by-region
	// " basis. To do so, put e.g.
	// "
	// " Atomic relaxations= On in Regions or
	// " Atomic relaxations= Off in regions
	// "
	// " in your input file. Then use
	// "
	// " Bound Compton start region=
	// " Bound Compton stop region=
	// " or
	// " Rayleigh start region=
	// " Rayleigh stop region=
	// " or
	// " Relaxations start region=
	// " Relaxations stop region=
	// " or
	// " PE sampling start region=
	// " PE sampling stop region=
	// "
	// " each followed by a lost of of one or more
	// " start and stop regions separated by commas.
	// " Example:
	// " Atomic relaxations= On in Regions
	// " Relaxations start region= 1, 40
	// " Relaxations stop region= 10, 99
	// " will first turn off relaxations everywhere and
	// " then turn on in regions 1-10 and 40-99.
	// " Note that input is checked against min. and max.
	// " region number and ignored if
	// " start region < 1 or stop_region > $MXREG or
	// " start region > stop region.
	// "
	// " ECUT, PCUT and SMAX can also be set on a
	// " region-by-region basis. To do so, iclude
	// " in your input file
	// "
	// " Set XXXX= f_value1, f_value2, ...
	// " Set XXXX start region= i_value1, i_value2, ...
	// " Set XXXX stop region= j_value1, j_value2, ...
	// "
	// " where XXXX is ECUT, PCUT or SMAX ,
	// " f_value1, f_value2,... are the desired values for XXXX
	// " and i_value_i and j_value_i are the start and
	// " stop regions.
	// "
	// "*******************************************************************************
	private int I = 0;
	private int J = 0;
	private int IQIN = 0;
	private int IRIN = 0;
	private int NCASE = 0;
	private double XIN = 0.;
	private double YIN = 0.;
	private double ZIN = 0.;
	private double EIN = 0.;
	private double WTIN = 0.;
	private double UIN = 0.;
	private double VIN = 0.;
	private double WIN = 0.;

	private double[] sc_array; // "for scoring energy deposited in all regions"
	private double[] sc_array2;// "for scoring energy squared on a history-by-history basis"
	private double[] sc_tmp;
	private int[] sc_last;
	private double[][] sc_pulse_height;// ($MXEBIN,$MXREG),
										// "for pulse height dstn"
	private double de_pulse = 0.0;
	// 0=Nai,1=HPGe,etc
	public static int detType = 0;
	private int $MXEBIN = 0;// 400;//200;
	private int icase = 0;
	private int ieff = 0;
	private double[] edep;
	private double[] edep2;
	private double[] edeptmp;
	private int[] iieff;// =0;
	private int irejct = 0;// "turn range rejection on (1) or off (0)"
	private double esave = 0.;// Don't use range rejection for electrons with
								// E>esave"

	private double[] zbound;// "array with z-plane co-ordinates"
	int nzb = 0;

	int nbatch = 0;
	int nperbatch = 0;
	int ibatch = 0;
	double aux = 0.0;
	double aux2 = 0.0;
	double total = 0.0;
	double anorm = 0.0;
	double sum = 0.0;
	double sum2 = 0.0;

	public int iterse = 0;// groups of planes or individual

	public int n_groups = 0;// number of groups of slabs
	public int[] num_slab;// number of slabs per group
	public double[] slabthickness;
	public double[] zfront;
	public double[] deltaz;
	public int ndepthboundaries = 0;
	public double[] planesz;
	private boolean de_set = false;

	public Tutor7_test()// auto init
	{
		de_set = false;
		init();
	}

	public Tutor7_test(double deltae) {
		de_pulse = deltae;
		de_set = true;
		init();
	}

	// default set de_pulse as fwhm, taking as reference the QUANTUM ASSAYER
	// user manual
	private void setBinWidth() {
		if (detType == 0) {
			double einc = EIN * 1000.0;// kev
			de_pulse = 53.0 * Math.sqrt(einc / 661.66);
			de_pulse = de_pulse / 1000.0;// MeV
		} else if (detType == 1) {
			double einc = EIN * 1000.0;// kev
			de_pulse = 0.6097 + 0.0011 * einc;
			de_pulse = de_pulse / 1000.0;// MeV
		}
	}

	private void init() {
		EGS4.startSimulationTime = System.currentTimeMillis();// grab start time
		// "---------------------------------------------------------------------"
		// "@@@@@STEP 1: USER-OVERRIDE-OF-EGSnrc-MACROS
		// "---------------------------------------------------------------------"
		EGS4.setMXMED(5);// "up to 5 media in the problem(default 10)
		EGS4.setMXREG(102);// "up to 102 geometric regions (default 2000)
		EGS4.setMXSTACK(50);// "less than 50 particles on stack at once"

		zbound = new double[EGS4.$MXREG];
		sc_array = new double[EGS4.$MXREG];
		sc_array2 = new double[EGS4.$MXREG];
		sc_tmp = new double[EGS4.$MXREG];
		sc_last = new int[EGS4.$MXREG];

		edep = new double[EGS4.$MXREG];
		edep2 = new double[EGS4.$MXREG];
		edeptmp = new double[EGS4.$MXREG];

		EGS4.egs_set_defaults();// EGS4.RandomUse=2;
		EGS4.eq = this;
		EGS4Core.eq = this;
		// " "Read" the input file "
		inputs();

		if (EGS4.STOPPROGRAM) {
			return;
		}
		// $RNG-INITIALIZATION;
		// REPLACE {$RNG-INITIALIZATION;} WITH {;
		// call init_ranlux($DEFAULT-LL,0);
		// call ranlux(rng_array); rng_seed = 1;
		// }
		EGS4.init_ranlux(EGS4.$DEFAULT_LL, 0);
		EGS4.ranlux(EGS4.rng_array);
		EGS4.rng_seed = 1;

		HATCH();
		if (EGS4.STOPPROGRAM) {
			return;
		}

		double xae = EGS4.AE[0] - 0.511;
		double xae2 = EGS4.AP[0];
		// ;OUTPUT AE(1)-0.511, AP(1);
		EGS4.seqStr = " knock-on electrons can be created and any electron followed down to "
				+ EGS4.format(xae, 8, true) + " MeV kinetic energy";// +" \n";
		if (EGS4.iprint > 1)
			printSequence(EGS4.seqStr);
		EGS4.seqStr = "brem photons can be created and any photon followed down to      "
				+ EGS4.format(xae2, 8, true) + " MeV";// +" \n";
		if (EGS4.iprint > 1)
			printSequence(EGS4.seqStr);
		// "Compton events can create electrons and photons below these cutoffs"
		// "---------------------------------------------------------------------"
		// "@@@@@STEP 7   SHOWER-CALL                                                 "
		// "---------------------------------------------------------------------"
		// "initiate the shower ncase times"
		// ;OUTPUT;(' Starting shower simulation ...');
		EGS4.seqStr = " Starting shower simulation ...";
		if (EGS4.iprint > 0)
			printSequence(EGS4.seqStr);

		nbatch = 10;
		nperbatch = NCASE / nbatch;
		if (nperbatch == 0)
			nperbatch = 1;
		// "Note that nbatch and nperbatch are not used for statistical analysis"
		// "They are merely for printing information about the progress of the"
		// "simulation"
		ieff = 0;
		iieff = new int[nbatch];
		long startTime = System.currentTimeMillis();
		for (icase = 1; icase <= NCASE; icase++)// "Note the use of icase as the shower counter"
		// "icase is in common/score/ and is used for"
		// "shower-by-shower scoring"
		{
			SHOWER();
			// if insucces:
			if (EGS4.STOPPROGRAM) {
				return;
			}

			ibatch = icase / nperbatch;
			if (ibatch * nperbatch == icase) {// "print every batch end"
				String timePerBatch = EGS4.timeElapsed(startTime);
				EGS4.seqStr = " Finished batch" + EGS4.format(ibatch, 7)
						+ " out of " + EGS4.format(nbatch, 7) + " time: "
						+ timePerBatch;
				startTime = System.currentTimeMillis();
				if (EGS4.iprint > 0)
					printSequence(EGS4.seqStr);
			}

			for (I = 1; I <= nzb + 1; I++) {
				edep[I - 1] = edep[I - 1] + edeptmp[I - 1];
				edep2[I - 1] = edep2[I - 1] + edeptmp[I - 1] * edeptmp[I - 1];

				if (I > 1 && I <= nzb) {
					if (edeptmp[I - 1] == EIN)// loss all energy
					{
						ieff++;
						if (ibatch == nbatch)// transfer to nbatch-1
							iieff[ibatch - 1]++;
						else
							iieff[ibatch]++;
					}
				}

				edeptmp[I - 1] = 0.0;// reset
			}

		}
		// ;OUTPUT ncase;(/' Finished shower simulation with', I10,' cases. '/);
		EGS4.seqStr = " Finished shower simulation with"
				+ EGS4.format(NCASE, 10) + " cases";
		if (EGS4.iprint > 0)
			printSequence(EGS4.seqStr);
		// "-----------------------------------------------------------------"
		// "@@@@STEP 8   OUTPUT-OF-RESULTS                                       "
		// "-----------------------------------------------------------------"
		sum = 0.0;
		sum2 = 0.0;
		double sume = 0.0;
		double sume2 = 0.0;
		double totale = 0.0;
		total = 0.0;
		Integer i3 = new Integer(IQIN);
		anorm = 1.0 / (EIN + i3.doubleValue() * EGS4.RM); // "for e+ add 2*rm to k.e."OK!!!
		// IF( iqin = 0 ) [ anorm = 1/ein; ] ELSE [ anorm = 1/(ein-rm); ]
		double auxe = 0.0;
		double auxe2 = 0.0;
		for (I = 1; I <= nzb + 1; I++) {
			// "first put non-scored energy portions into sc_array and sc_array2"
			aux = sc_tmp[I - 1];
			aux2 = aux * aux;
			sc_array[I - 1] = sc_array[I - 1] + aux;
			sc_array2[I - 1] = sc_array2[I - 1] + aux2;

			// if(I>1 && I<=nzb)
			// {
			// if(aux==EIN)//loss all energy
			// {
			// ieff++;
			// }
			// }

			aux = sc_array[I - 1] / NCASE;
			aux2 = sc_array2[I - 1] / NCASE;
			aux2 = (aux2 - aux * aux) / (NCASE - 1);
			if (aux2 > 0)
				aux2 = Math.sqrt(aux2);
			// errors->(x2med-xmed2)/(n-1)!!OK
			aux = aux * anorm;
			aux2 = aux2 * anorm;
			sc_array[I - 1] = aux;
			sc_array2[I - 1] = aux2;

			auxe = edep[I - 1] / NCASE;
			auxe2 = edep2[I - 1] / NCASE;
			auxe2 = (auxe2 - auxe * auxe) / (NCASE - 1);
			if (auxe2 > 0)
				auxe2 = Math.sqrt(auxe2);
			// errors->(x2med-xmed2)/(n-1)!!OK
			auxe = auxe * anorm;
			auxe2 = auxe2 * anorm;
			edep[I - 1] = auxe;
			edep2[I - 1] = auxe2;

			if ((I > 1) && (I <= nzb)) {
				sume = sume + auxe;
				sume2 = sume2 + auxe2 * auxe2;
			}

			totale = totale + auxe;

			if ((I > 1) && (I <= nzb)) {
				sum = sum + aux;
				sum2 = sum2 + aux2 * aux2;
			}

			total = total + aux;

		}
		sum2 = Math.sqrt(sum2);
		sume2 = Math.sqrt(sume2);
		EGS4.seqStr = " ********************************************** ";
		if (EGS4.iprint > 0)
			printSequence(EGS4.seqStr);

		EGS4.seqStr = "   Reflected energy fraction: "
				+ EGS4.format(sc_array[0], 10, true) + " +/- "
				+ EGS4.format(sc_array2[0], 10, true);
		if (EGS4.iprint > 0)
			printSequence(EGS4.seqStr);
		EGS4.seqStr = "   @Reflected energy fraction: "
				+ EGS4.format(edep[0], 10, true) + " +/- "
				+ EGS4.format(edep2[0], 10, true);
		if (EGS4.iprint > 0)
			printSequence(EGS4.seqStr);

		EGS4.seqStr = "   Deposited energy fraction: "
				+ EGS4.format(sum, 10, true) + " +/- "
				+ EGS4.format(sum2, 10, true);
		if (EGS4.iprint > 0)
			printSequence(EGS4.seqStr);
		EGS4.seqStr = "   @Deposited energy fraction: "
				+ EGS4.format(sume, 10, true) + " +/- "
				+ EGS4.format(sume2, 10, true);
		if (EGS4.iprint > 0)
			printSequence(EGS4.seqStr);

		EGS4.seqStr = " Transmitted energy fraction: "
				+ EGS4.format(sc_array[nzb], 10, true) + " +/- "
				+ EGS4.format(sc_array2[nzb], 10, true);
		if (EGS4.iprint > 0)
			printSequence(EGS4.seqStr);
		EGS4.seqStr = " @Transmitted energy fraction: "
				+ EGS4.format(edep[nzb], 10, true) + " +/- "
				+ EGS4.format(edep2[nzb], 10, true);
		if (EGS4.iprint > 0)
			printSequence(EGS4.seqStr);

		EGS4.seqStr = " -------------------------------------------------------------";
		if (EGS4.iprint > 0)
			printSequence(EGS4.seqStr);

		EGS4.seqStr = "                       total: "
				+ EGS4.format(total, 10, true);
		if (EGS4.iprint > 0)
			printSequence(EGS4.seqStr);
		EGS4.seqStr = "                       @total: "
				+ EGS4.format(totale, 10, true);
		if (EGS4.iprint > 0)
			printSequence(EGS4.seqStr);

		// "print detailed result in the list file"
		EGS4.seqStr = " ------------------------DETAILED-------------------------------------";
		if (EGS4.iprint > 0)
			printSequence(EGS4.seqStr);

		EGS4.seqStr = "   depth " + "    " + " deposited energy fraction "
				+ "    " + " error ";
		if (EGS4.iprint > 0)
			printSequence(EGS4.seqStr);

		EGS4.seqStr = " -------------------------------------------------------------";
		if (EGS4.iprint > 0)
			printSequence(EGS4.seqStr);

		EGS4.setDigits(6);
		for (I = 2; I <= nzb; I++) {
			double dbl = 0.5 * (zbound[I - 2] + zbound[I - 1]);
			EGS4.seqStr = EGS4.format(dbl, 10, true) + "        "
					+ EGS4.format(sc_array[I - 1], 11, true) + "            "
					+ EGS4.format(sc_array2[I - 1], 11, true);
			if (EGS4.iprint > 0)
				printSequence(EGS4.seqStr);

		}

		EGS4.seqStr = " -------------------------------------------------------------";
		if (EGS4.iprint > 0)
			printSequence(EGS4.seqStr);

		EGS4.seqStr = " ---------------- Response function--------------";
		if (EGS4.iprint > 0)
			printSequence(EGS4.seqStr);

		EGS4.seqStr = " =============================================================";
		if (EGS4.iprint > 0)
			printSequence(EGS4.seqStr);

		EGS4.seqStr = " Energy	Counts/incident photon     error";
		if (EGS4.iprint > 0)
			printSequence(EGS4.seqStr);

		EGS4.seqStr = " -------------------------------------------------------------";
		if (EGS4.iprint > 0)
			printSequence(EGS4.seqStr);

		for (I = 1; I <= nzb + 1; I++) {
			EGS4.seqStr = " ********* Region " + I;
			if (EGS4.iprint > 0)
				printSequence(EGS4.seqStr);

			// "Score also the remaining pulses"
			aux = sc_tmp[I - 1];// (i);
			Double db = new Double(1.0 + aux / de_pulse);
			J = db.intValue();
			if ((J > 0) && (J <= $MXEBIN)) {
				sc_pulse_height[J - 1][I - 1] = sc_pulse_height[J - 1][I - 1] + 1.0;
			}

			for (J = 1; J <= $MXEBIN; J++) {
				aux = sc_pulse_height[J - 1][I - 1] / NCASE;
				aux2 = Math.sqrt(aux * (1 - aux) / (NCASE - 1));
				// aux = aux/de_pulse; aux2 = aux2/de_pulse;

				double dd = (J) * de_pulse;
				EGS4.seqStr = EGS4.format(dd, 9, true) + EGS4.format(aux, 15)
						+ EGS4.format(aux2, 15);
				if (EGS4.iprint > 0)
					printSequence(EGS4.seqStr);

			}

		}

		EGS4.seqStr = " ********************************************** ";
		if (EGS4.iprint > 0)
			printSequence(EGS4.seqStr);
		Integer intg = new Integer(NCASE);
		double dbl2 = intg.doubleValue();
		double dbl = 100 * ieff / dbl2;
		EGS4.seqStr = "Eff (%)=" + EGS4.format(dbl, 9, true);// WORKS
		if (EGS4.iprint > 0)
			printSequence(EGS4.seqStr);
		double[] effi = new double[nbatch];
		double[] effi2 = new double[nbatch];
		Integer itg = new Integer(nperbatch);
		double nd = itg.doubleValue();
		sume = 0.0;
		sume2 = 0.0;
		for (int i = 1; i <= nbatch; i++) {
			effi[i - 1] = iieff[i - 1] * 100 / nd;
			effi2[i - 1] = effi[i - 1] * effi[i - 1];
			sume = sume + effi[i - 1];
			sume2 = sume2 + effi2[i - 1];
		}
		double effm = sume / nbatch;
		double seffm = sume2 / nbatch;
		seffm = (seffm - effm * effm) / (nbatch - 1);
		if (seffm > 0)
			seffm = Math.sqrt(seffm);
		// errors->(x2med-xmed2)/(n-1)!!OK
		EGS4.seqStr = "@Eff (%)=" + EGS4.format(effm, 9, true) + "  +/-  "
				+ EGS4.format(seffm, 9, true);// WORKS
		if (EGS4.iprint > 0)
			printSequence(EGS4.seqStr);

		EGS4.timeElapsed();
		Calendar cal = Calendar.getInstance();
		Date d = cal.getTime();
		EGS4.seqStr = "End of run:          " + d.toString();
		if (EGS4.iprint > 0)
			printSequence(EGS4.seqStr);

		// "-----------------------------------------------------------------"
		// "@@@@@STEP 9   finish run                                              "
		// Java can handle finishing by itself.It auto-runs the Garbage
		// collector!!!
		// "-----------------------------------------------------------------"

	}

	private void HATCH() {
		Calendar cal = Calendar.getInstance();
		Date d = cal.getTime();
		EGS4.seqStr = "Start of run:         " + d.toString();
		if (EGS4.iprint > 0)
			printSequence(EGS4.seqStr);

		EGS4.HATCH();
	}

	private void SHOWER() {
		EGS4Core.SHOWER(IQIN, EIN, XIN, YIN, ZIN, UIN, VIN, WIN, IRIN, WTIN);
	}

	public void printSequence(String s) {

		System.out.println(s);
	}

	public void AUSGAB(int iarg) {
		// ;Copyright NRC;
		int irl = 0;
		int jp = 0;

		if (EGS4.NP >= EGS4.$MXSTACK) {// "STACK is as deep as allowed"
										// OUTPUT NP,$MXSTACK;(//' In AUSGAB,
										// NP=',I5,' >= Maximum value allowed=',
										// I3/' Adjust $MXSTACK =',I5,',
										// recompile and try again
										// '/1X,80('*')/);
										// STOP;
			EGS4.STOPPROGRAM = true;
			EGS4.seqStr = " ***************************************************"
					+ "  \n"
					+ " In AUSGAB, NP="
					+ EGS4.NP
					+ " >= Maximum value allowed="
					+ EGS4.$MXSTACK
					+ "  \n"
					+ "  Adjust $MXSTACK, recompile and try again "
					+ "  \n"
					+ " ***************************************************";
			// if(EGS4.iprint>0)
			printSequence(EGS4.seqStr);

			return;// stop;
		}
		// "Note the above is not foolproof"
		if (iarg < 5) {// "energy is being deposited"
			irl = EGS4.IR[EGS4.NP - 1];// ir(np);
			edeptmp[irl - 1] = edeptmp[irl - 1] + EGS4.EDEP;
			if (icase == sc_last[irl - 1]) {// "still the same shower that deposited energy"
											// "last time in this region"
				sc_tmp[irl - 1] = sc_tmp[irl - 1] + // edep*wt(np);
						EGS4.EDEP;// *EGS4.WT[EGS4.NP-1];
				// "OUTPUT iarg,edep,e(np),iq(np); "
				// "(' scoring (1) ',f10.3,f10.3,f10.3,f10.3); "
			} else {// "we have the next shower depositing energy into region irl"
					// " => put sc_tmp into  the scoring arrays and set sc_last"
				aux = sc_tmp[irl - 1];
				sc_array[irl - 1] = sc_array[irl - 1] + aux;
				sc_array2[irl - 1] = sc_array2[irl - 1] + aux * aux;
				sc_tmp[irl - 1] = EGS4.EDEP;// *EGS4.WT[EGS4.NP-1];//edep*wt(np);
				sc_last[irl - 1] = icase;

				Double dbl = new Double(1.0 + aux / de_pulse);
				jp = dbl.intValue();
				if ((jp > 0) && (jp <= $MXEBIN)) {
					sc_pulse_height[jp - 1][irl - 1] = sc_pulse_height[jp - 1][irl - 1] + 1.0;
				}

				// *****************************
				// if ((irl > 1)&&(irl <= nzb))
				// {
				// if(aux==EIN)//loss all energy
				// ieff++;
				// }
				// *****************************

			}
		}

	}// "END OF AUSGAB"

	// /"*********************************************************************"
	public void HOWFAR() {
		// ;Copyright NRC;

		double tval = 0.0;
		int irl = EGS4.IR[EGS4.NP - 1];// ir(np);//"region number and direction into local variables"
		double wl = EGS4.W[EGS4.NP - 1];// w(np);

		if (irl > nzb)// "past the geometry ?"
		{
			if (irl > nzb + 1) {// "This should not happen";
								// OUTPUT irl; (' irl > nzb+1 !? ',f10.3);
								// EGS4.OUTPUTs=EGS4.OUTPUTs+" irl > nzb+1 !? "+irl;
			}
			if (wl > 0.0) {// " yes, terminate the history"
				EGS4.IDISC = 1;
			} else {
				// "No. This should not happen for exact boundary crossing but "
				// "possible to happen with boundary crossing a la PRESTA"
				// "(particle reflected at the boundary)"
				EGS4.USTEP = 0.0;
				EGS4.IRNEW = nzb;// ustep = 0; irnew = nzb;
			}

			return;
		} else if (irl > 1) {// "in the geometry, do transport checks"
			if (wl > 0.0) {
				// "going forward
				tval = (zbound[irl - 1] - EGS4.Z[EGS4.NP - 1]) / wl;
				// "TVAL is dist to boundary in this direction"
				if (tval > EGS4.USTEP) {// "can take currently requested step"
					EGS4.IRNEW = irl;
				} else {// "requested step longer than distance to boundary => adjust"
					EGS4.USTEP = tval;
					EGS4.IRNEW = irl + 1;
				}
			} else if (wl < 0.0) {// "going back towards origin"
									// tval = (zbound(irl-1)-z(np))/wl;
				tval = (zbound[irl - 2] - EGS4.Z[EGS4.NP - 1]) / wl;
				if (tval > EGS4.USTEP) {// "can take currently requested step"
					EGS4.IRNEW = irl;
				} else {// "requested step longer than distance to boundary => adjust"
					EGS4.USTEP = tval;
					EGS4.IRNEW = irl - 1;
				}
			} else {// "going parallel to the boundary"
				EGS4.IRNEW = irl;
			}
			return;
		}
		// "at this point it is clear that the particle is in front of the geometry"
		if (EGS4.W[EGS4.NP - 1] < 0)// w(np) < 0)
		{// "This is a backscattered particle, discard it"
			EGS4.IDISC = 1;
		} else {// "this is either a particle reflected on the boundary (possible "
				// "for PRESTA) or a particle with an incorrectly initialized entry region"
			EGS4.USTEP = 0.0;// ustep = 0;
			EGS4.IRNEW = 2;// irnew = 2;
		}

	}// "END OF SUBROUTINE HOWFAR"

	// "*********************************************************************"
	// CALL HOWNEAR({P1},X(NP),Y(NP),Z(NP),IRL);
	public void HOWNEAR()// double tperp, double x,double y,double z,int irl)
	{
		// " The following is a general specification of HOWNEAR                 "
		// "   Given a particle at (x,y,z) in region irl, HOWNEAR answers the    "
		// "   question, What is the distance tperp to the closest boundary?     "
		// "  In general this can be a complex subroutine.                       "
		// "*********************************************************************"
		// ;Copyright NRC;
		// ;COMIN/GEOM/; // "       COMMON GEOM contains ZBOUND"
		// ####################---------------------------------
		double z = EGS4.Z[EGS4.NP - 1];
	//	double y = EGS4.Y[EGS4.NP - 1];
	//	double x = EGS4.X[EGS4.NP - 1];
		// ####################---------------------------------
		int irl = EGS4.irl;

		if ((irl > 1) || (irl <= nzb)) {// "particle is in the geometry"
										// tperp = min(z - zbound(irl-1),
										// zbound(irl) - z);
			EGS4.tperp = Math.min(z - zbound[irl - 2], zbound[irl - 1] - z);
		} else {
			EGS4.tperp = 0.0;
		}

	}// "end of subroutine HOWNEAR"

	private void inputs() {
		// ###########################################call
		// get_transport_parameter(6);#######
		// Global ECUT ->ECUT
		// Global PCUT ->PCUT
		EGS4.ECUT[0] = 0.711;
		EGS4.PCUT[0] = 0.010;
		// "Now set ecut and pcut to the values input by the user"
		for (I = 2; I <= EGS4.$MXREG; I++) {
			EGS4.ECUT[I - 1] = EGS4.ECUT[0];
			EGS4.PCUT[I - 1] = EGS4.PCUT[0];
		}
		// Global SMAX ->SMAXIR
		EGS4.SMAXIR[0] = 1.e10;
		for (I = 2; I <= EGS4.$MXREG; I++) {
			EGS4.SMAXIR[I - 1] = EGS4.SMAXIR[0];
		}
		// Bound Compton scattering= On or Off->IBCMP
		// #######//On means 1 and Off means 0
		// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		EGS4.ibcmp[0] = 1;
		// "Now set ibcmp for all regions to the value input by the user"
		for (I = 2; I <= EGS4.$MXREG; I++) {
			EGS4.ibcmp[I - 1] = EGS4.ibcmp[0];
		}
		// Rayleigh scattering= On->IRAYLR
		EGS4.IRAYLR[0] = 1;
		// "Now set iraylr for all regions to the value input by the user"
		for (I = 2; I <= EGS4.$MXREG; I++) {
			EGS4.IRAYLR[I - 1] = EGS4.IRAYLR[0];
		}
		// Atomic relaxations= On->IEDGFL
		EGS4.iedgfl[0] = 1;
		// "Now set iedgfl for all regions to the value input by the user"
		for (I = 2; I <= EGS4.$MXREG; I++) {
			EGS4.iedgfl[I - 1] = EGS4.iedgfl[0];
		}
		// Photoelectron angular sampling= On->IPHTER
		EGS4.iphter[0] = 1;
		// "Now set iphter for all regions to the value input by the user"
		for (I = 2; I <= EGS4.$MXREG; I++) {
			EGS4.iphter[I - 1] = EGS4.iphter[0];
		}
		// Brems angular sampling= On ->IBRDST
		// ' 0: leading term of Koch and Motz distn'/
		// ' 1: Koch and Motz 2BS(modified): '/
		EGS4.ibrdst = 1;
		if ((EGS4.ibrdst < 0) || (EGS4.ibrdst > 1)) {
			EGS4.ibrdst = 1;
		}
		// Pair angular sampling= On->IPRDST
		// ' 0: fixed pair angle (EGS4 default) '/
		// ' 1: leading term of the distribution '/
		// ' 2: Koch and Motz ');
		EGS4.iprdst = 1;
		if ((EGS4.iprdst < 0) || (EGS4.iprdst > 2)) {
			EGS4.iprdst = 1;
		}
		// ESTEPE= 0.25->ESTEPE
		EGS4.estepe = 0.25;
		if ((EGS4.estepe <= 0) || (EGS4.estepe >= 1)) {
			EGS4.estepe = EGS4.$MAX_ELOSS; // "$MAX-ELOSS is defined in egsnrc.macros at 0.25"
		}
		// XIMAX= 0.5->XIMAX
		EGS4.ximax = 0.5;
		if ((EGS4.ximax <= 0) || (EGS4.ximax >= 1)) {
			EGS4.ximax = EGS4.$EXACT_BCA_XIMAX; // "$EXACT-BCA-XIMAX set to 0.5 in egsnrc.macros"
		}
		// Boundary crossing algorithm= exact->bca_algorithm, exact_bca
		EGS4.bca_algorithm = 0;// exact means =EGS4.$BCA_ALGORITHM_DEFAULT=0;
		if ((EGS4.bca_algorithm < 0) || (EGS4.bca_algorithm > 1)) {
			EGS4.bca_algorithm = 0;
		}
		// Skin depth for BCA= 3->skindepth_for_bca
		EGS4.skindepth_for_bca = 3;
		if (EGS4.bca_algorithm == 0) {
			// INPUT skindepth_for_bca; (F10.0);
			if (EGS4.skindepth_for_bca <= 0)
				EGS4.skindepth_for_bca = 3;
		}
		// Electron-step algorithm= default->transport_algorithm
		// PRESTA-II (the default),$PRESTA_II = 0;
		EGS4.transport_algorithm = 0;
		if ((EGS4.transport_algorithm < 0) || (EGS4.transport_algorithm > 1)) {
			EGS4.transport_algorithm = 0;
		}
		// Spin effects= On->spin_effects
		EGS4.spin_effects = true;
		// ##########################END TRANSPORT
		// PARAM#####################################
		// ######################GEOMETRY###################################################
		iterse = 0;
		// "User wants to input groups of planes "
		if (iterse == 0)// //"User wants to input groups of planes "
		{
			// %%%%%%%%%%%%%INPUT%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			n_groups = 1;
			// %%%%%%%%%%%%%END INPUT%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			num_slab = new int[n_groups];
			slabthickness = new double[n_groups];
			zfront = new double[n_groups];
			deltaz = new double[n_groups];
			// %%%%%%%%%%%%%INPUT%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			num_slab[0] = 1;// 'NUMBER OF SLABS'; "The code word"
			zfront[0] = 0.0;// value(num_zfront,1);==='Z OF FRONT FACE'
							// "The code word"
			deltaz[0] = 6.3;// 1.0;////slab thickness= 0.01;value(num_deltaz,i)
			// %%%%%%%%%%%%%END INPUT%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

			zbound[0] = zfront[0];
			nzb = 1;
			for (I = 1; I <= n_groups; I++) {
				for (J = 1; J <= num_slab[I - 1]; J++) {
					nzb = nzb + 1;
					if (nzb > EGS4.$MXREG) {
						EGS4.STOPPROGRAM = true;
						EGS4.seqStr = " ***************************************************"
								+ "  \n"
								+ " You are asking for too many planes! "
								+ " maximum allowed is "
								+ EGS4.$MXREG
								+ "  \n"
								+ " ***************************************************";
						// if(EGS4.iprint>0)
						printSequence(EGS4.seqStr);
						return;
					}
					zbound[nzb - 1] = zbound[nzb - 2] + deltaz[I - 1];
				}
			}
		} else {
			// %%%%%%%%%%%%%INPUT%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			ndepthboundaries = 2;// for instance
			// %%%%%%%%%%%%%END INPUT%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			planesz = new double[ndepthboundaries];
			// %%%%%%%%%%%%%INPUT%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			planesz[0] = 0.0;
			planesz[1] = 0.5;
			// %%%%%%%%%%%%%END INPUT%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			nzb = ndepthboundaries;
			for (I = 1; I <= nzb; I++) {
				zbound[I - 1] = planesz[I - 1];
			}

		}
		// ######################END
		// GEOMETRY##############################################
		// #####################" Media input "##########################################
		EGS4.MEDIA[0] = "NAI_fortran";// "GE_fortran";//"TA_fortran";
		// "Set medium 1 everywhere"
		for (I = 1; I <= nzb; I++) {
			EGS4.MED[I - 1] = 1;
		}
		// " and the vacuum in region 1 and nzb+1 "
		EGS4.MED[0] = 0;
		EGS4.MED[nzb] = 0;// med(1),med(nzb+1)/ = 0;
		// #####################" END Media input "##########################################
		// override for user media-no case here
		// med(j) = value(num_set,i);
		// ----------------------------

		EGS4.DUNIT = 1.0; // "i.e. we work in cm"
		// "Define range rejection parameter. Although not directly related "
		// "to ausgab, range rejection is an `user' variance reduction technique "
		// "and so, this is the most appropriate place to initialize it"

		// Range Rejection= On
		EGS4.iurd = EGS4.iTut7;
		irejct = 0;// WE MAKE IT OFF
		EGS4Macro.irejct = irejct;
		// Esave= 1
		if (irejct == 1) {
			esave = 1.0;
			EGS4Macro.esave = esave;
		}
		// " Source inputs "
		IQIN = 0;// 'INCIDENT CHARGE'
		EIN = 0.662;// 'INCIDENT KINETIC ENERGY'
		if (IQIN != 0)
			EIN = EIN + EGS4.RM;
		NCASE = 40000;// 5000000;//'NUMBER OF HISTORIES'
		WIN = 0.0;// 'INCIDENT ANGLE'
		WIN = WIN / 180 * Math.PI;
		WIN = Math.cos(WIN);
		// choosing lets say VIN=0, then uin is equal with sin!!
		UIN = Math.sqrt(Math.max(0.0, (1 - WIN) * (1 + WIN)));
		VIN = 0.0;
		IRIN = 2; // "starts in region 2, could be 1"
		WTIN = 1; // "statistical weight is 1"
		XIN = 0.0;
		YIN = 0.0;
		ZIN = zbound[0];

		if (!de_set)
			setBinWidth();

		Double ddd = new Double(1.0 + EIN / de_pulse);
		$MXEBIN = ddd.intValue();

		sc_pulse_height = new double[$MXEBIN][EGS4.$MXREG];
		// "Set all scoring arrays to zero. This could be avoided if"
		// "the compiler being used has a `initialize to zero' option"
		// "It is a good coding habit to not rely on variables being"
		// "automatically zeroed"
		for (I = 1; I <= EGS4.$MXREG; I++) {
			sc_array[I - 1] = 0.0;
			sc_array2[I - 1] = 0.0;
			sc_tmp[I - 1] = 0.0;
			sc_last[I - 1] = 0;
			for (J = 1; J <= $MXEBIN; J++) {
				sc_pulse_height[J - 1][I - 1] = 0.0;
			}
		}
		return;
	}
}
