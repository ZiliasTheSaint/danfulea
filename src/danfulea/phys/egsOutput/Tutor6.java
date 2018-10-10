package danfulea.phys.egsOutput;

import java.util.Calendar;
import java.util.Date;

import danfulea.phys.egs.EGS4;
import danfulea.phys.egs.EGS4Core;
import danfulea.phys.egs.EgsQuestion;

/**
 * Demo class
 * Based on EGS4/EGSnrc routines developed by SLAC(Stanford Linear Accelerator Center)/NRC (National Research Council of Canada). 
 * @author Dan Fulea, 17 OCT. 2005
 */

// " An EGSnrc user code. It calculates reflected, deposited and         "
// " transmitted energy for electron and photon beams incident on        "
// " a slab geometry.                                                    "

public class Tutor6 implements EgsQuestion {
	private int I = 0;
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
	private double[] sc_array2;// "for scoring energy squared on a history-by-"
								// "history basis"
	private double[] sc_tmp;
	private int[] sc_last;
	private int ispin = 0;
	private int icase = 0;
	private int irejct = 0;
	private double esave = 0.;
	private double[] zbound;
	int nzb = 0;

	int nbatch = 0;
	int nperbatch = 0;
	int ibatch = 0;
	double aux = 0.0;
	double aux2 = 0.0;
	double total = 0.0;
	double anorm = 0.0;

	public Tutor6() {
		init();
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
		// REPLACE {$CALL-HOWNEAR(#);} WITH {
		// ;call hownear({P1},z(np),ir(np));
		// "---------------------------------------------------------------------"
		// "@@@@@STEP 2 PRE-HATCH-CALL-INITIALIZATION
		// "---------------------------------------------------------------------"
		EGS4.egs_set_defaults();
		EGS4.eq = this;
		EGS4Core.eq = this;
		// " Read the input file "
		inputs();
		// "---------------------------------------------------------------------"
		// "@@@@@STEP 3   HATCH-CALL                                                  "
		// "---------------------------------------------------------------------"
		// ' CALL HATCH to get cross-section data'/)
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
			if (icase % nperbatch == 0) {// "print every batch end"
											// OUTPUT
											// float(100*icase)/float(ncase);
											// (' Finished ',F7.1,'% of cases');
				Integer i1 = new Integer(icase);
				Integer i2 = new Integer(NCASE);
				double dbl = 100. * i1.doubleValue() / i2.doubleValue();
				EGS4.seqStr = " Finished " + EGS4.format(dbl, 7, true)
						+ " % of cases";
				if (EGS4.iprint > 0)
					printSequence(EGS4.seqStr);
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
		total = 0.0;
		Integer i3 = new Integer(IQIN);
		anorm = 1.0 / (EIN + i3.doubleValue() * EGS4.RM); // "for e+ add 2*rm to k.e."OK!!!

		for (I = 1; I <= nzb + 1; I++) {
			// "first put non-scored energy portions into sc_array and sc_array2"
			aux = sc_tmp[I - 1];
			aux2 = aux * aux;
			sc_array[I - 1] = sc_array[I - 1] + aux;
			sc_array2[I - 1] = sc_array2[I - 1] + aux2;
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
			// "OUTPUT i,aux,aux2;"
			// "  (' region ',i3,' deposited fraction ',f10.6,' +/- ',f10.6);"
			total = total + aux;
		}
		EGS4.seqStr = " ********************************************** ";
		if (EGS4.iprint > 0)
			printSequence(EGS4.seqStr);

		EGS4.seqStr = "   Reflected energy fraction: "
				+ EGS4.format(sc_array[0], 10, true) + " +/- "
				+ EGS4.format(sc_array2[0], 10, true);
		if (EGS4.iprint > 0)
			printSequence(EGS4.seqStr);

		EGS4.seqStr = "   Deposited energy fraction: "
				+ EGS4.format(sc_array[1], 10, true) + " +/- "
				+ EGS4.format(sc_array2[1], 10, true);
		if (EGS4.iprint > 0)
			printSequence(EGS4.seqStr);

		EGS4.seqStr = " Transmitted energy fraction: "
				+ EGS4.format(sc_array[2], 10, true) + " +/- "
				+ EGS4.format(sc_array2[2], 10, true);
		if (EGS4.iprint > 0)
			printSequence(EGS4.seqStr);

		EGS4.seqStr = " -------------------------------------------------------------";
		if (EGS4.iprint > 0)
			printSequence(EGS4.seqStr);

		EGS4.seqStr = "                       total: "
				+ EGS4.format(total, 10, true);
		if (EGS4.iprint > 0)
			printSequence(EGS4.seqStr);

		// OUTPUT sc_array(1),sc_array2(1);
		// (' Reflected energy fraction: ',f10.6,' +/- ',f10.6);
		// OUTPUT sc_array(2),sc_array2(2);
		// (' Deposited energy fraction: ',f10.6,' +/- ',f10.6);
		// OUTPUT sc_array(3),sc_array2(3);
		// (' Transmitted energy fraction: ',f10.6,' +/- ',f10.6);
		// OUTPUT;
		// ('-------------------------------------------------------------');
		// OUTPUT total; (' total: ',f10.6///);

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

		if (EGS4.NP >= EGS4.$MXSTACK) {// "STACK is as deep as allowed"
										// OUTPUT NP,$MXSTACK;(//' In AUSGAB,
										// NP=',I5,' >= Maximum value allowed=',
										// I3/' Adjust $MXSTACK =',I5,',
										// recompile and try again
										// '/1X,80('*')/);
										// STOP;
			EGS4.STOPPROGRAM = true;
			EGS4.seqStr = " ********************************************** ";
			if (EGS4.iprint > 0)
				printSequence(EGS4.seqStr);
			EGS4.seqStr = " In AUSGAB, NP=" + EGS4.NP
					+ " >= Maximum value allowed=" + EGS4.$MXSTACK;
			if (EGS4.iprint > 0)
				printSequence(EGS4.seqStr);
			EGS4.seqStr = "  Adjust $MXSTACK, recompile and try again ";
			if (EGS4.iprint > 0)
				printSequence(EGS4.seqStr);
			EGS4.seqStr = " ********************************************** ";
			if (EGS4.iprint > 0)
				printSequence(EGS4.seqStr);

			return;// stop;
		}
		// "Note the above is not foolproof"
		if (iarg < 5) {// "energy is being deposited"
			irl = EGS4.IR[EGS4.NP - 1];// ir(np);
			if (icase == sc_last[irl - 1]) {// "still the same shower that deposited energy"
											// "last time in this region"
				sc_tmp[irl - 1] = sc_tmp[irl - 1] + // edep*wt(np);
						EGS4.EDEP * EGS4.WT[EGS4.NP - 1];
				// "OUTPUT iarg,edep,e(np),iq(np); "
				// "(' scoring (1) ',f10.3,f10.3,f10.3,f10.3); "
			} else {// "we have the next shower depositing energy into region irl"
					// " => put sc_tmp into  the scoring arrays and set sc_last"
				aux = sc_tmp[irl - 1];
				sc_array[irl - 1] = sc_array[irl - 1] + aux;
				sc_array2[irl - 1] = sc_array2[irl - 1] + aux * aux;
				sc_tmp[irl - 1] = EGS4.EDEP * EGS4.WT[EGS4.NP - 1];// edep*wt(np);
				sc_last[irl - 1] = icase;
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
		//double y = EGS4.Y[EGS4.NP - 1];
		//double x = EGS4.X[EGS4.NP - 1];
		// ####################---------------------------------
		int irl = EGS4.IR[EGS4.NP - 1];// EGS4.irl;

		if ((irl > 1) || (irl <= nzb)) {// "particle is in the geometry"
										// tperp = min(z - zbound(irl-1),
										// zbound(irl) - z);
			EGS4.tperp = Math.min(z - zbound[irl - 2], zbound[irl - 1] - z);
		} else {
			EGS4.tperp = 0.0;
		}

	}// "end of subroutine HOWNEAR"

	private void inputs() {
		// " This routine handles the input of all necessary parameter and passes"
		// " control to the EGSnrc system by setting various flags               "
		// "                                                                     "
		// " It handles the following standard ``steps'' in an EGSnrc user code  "
		// "                                                                     "
		// " Steps 2 (pre-hatch initializations),                                "
		// "       4 (initializations for howfar and hownear)                    "
		// "       5 (initializations for ausgab)                                "
		// "       6 (determination of incident particle parameters)             "
		// " are done is subroutine inputs                                       "
		// "                                                                     "
		// " Copyright National Research Council of Canada 2000                  "
		// "*********************************************************************"
		// ;Copyright NRC;

		// "The following are EGSnrc internal (private) common blocks"
		// "They are included in order to get access to various switches, "
		// "material array, cut-off energies, etc."

		// ;COMIN/BOUNDS, "to get access to ecut and pcut"
		// BREMPR, "to get access to ibrdst, iprdst and ibr_nist"
		// COMPTON-DATA, "to get access to the array ibcmp"
		// EDGE, "to get access to the arrays iedgfl and iphter"
		// EGS-VARIANCE-REDUCTION, "to get access to i_do"
		// ET-Control, "to get access to estepe,ximax,"
		// "                 skindepth_for_bca,transport_algorithm,"
		// "                 bca_algorithm,exact_bca,spin_effects"
		// MEDIA, "to get access to nmed and media"
		// MISC, "to get access to the arrays med and iraylr"
		// RANDOM, "to give access to random number seeds"
		// UPHIOT, "to get PI"
		// USEFUL, "to get electron rest energy RM"
		// "                                             "
		// " The following are user-defined common blocks"
		// "                                             "
		// GEOM,
		// SCORE,
		// SOURCE,
		// USER
		// /;
		// character*60 medium_name;
		// $INTEGER i,j,ispin,irejct,luxury_level,iseed;
		// $REAL esave;
		// "---------------------------------------------------------------------"
		// "STEP 2 PRE-HATCH-CALL-INITIALIZATION                                 "
		// "---------------------------------------------------------------------"

		EGS4.NMED = 1;// nmed = 1;
						// //"in this version of the tutor6 code, use just one medium"
		EGS4.DUNIT = 1;// dunit = 1; //"i.e. we work in cm"

		// DO i=1,nmed [
		// OUTPUT i; (/' Input name of medium ',I3,':',$);
		// INPUT medium_name; (A);
		// DO j=1,24 [ media(j,1) = medium_name(j:j); ]
		// ]
		// try
		// {
		// FileInputStream in = new FileInputStream(nume);

		// EGS4.OUTPUTs=EGS4.OUTPUTs+" Input name of medium "+1+"\n";
		// OUTPUTs=OUTPUTs+" Input name of medium "+1+"\n";
		EGS4.MEDIA[0] = "TA_fortran";// "tantalum";//"tantalum";//"TA_fortran";
		// "Set medium 1 everywhere"
		for (I = 1; I <= EGS4.$MXREG; I++) {
			EGS4.MED[I - 1] = 1;
		}

		// OUTPUT; (/' Input minimum electron transport energy (total, MeV):
		// ',$);
		// EGS4.OUTPUTs=EGS4.OUTPUTs+" Input minimum electron transport energy (total, MeV): "+"\n";
		// OUTPUTs=OUTPUTs+" Input minimum electron transport energy (total, MeV): "+"\n";
		// INPUT ecut(1); (F10.0);
		// OUTPUT; (/' Input minimum photon transport energy (MeV): ',$);
		// INPUT pcut(1); (F10.0);
		EGS4.ECUT[0] = 1.0;
		EGS4.PCUT[0] = 0.01;
		// "Now set ecut and pcut to the values input by the user"
		for (I = 2; I <= EGS4.$MXREG; I++) {
			EGS4.ECUT[I - 1] = EGS4.ECUT[0];
			EGS4.PCUT[I - 1] = EGS4.PCUT[0];
		}

		// "Rayleigh switch, must be input in step 2 so that HATCH can check whether"
		// "Rayleight data is available"
		// while(true)
		// {
		// OUTPUT; (/' Rayleigh scattering on (1) or off (0)? ',$);
		// INPUT iraylr(1); (I5);
		// }// UNTIL ( iraylr(1) = 1 | iraylr(1) = 0 );
		EGS4.IRAYLR[0] = 0;
		// "Now set iraylr for all regions to the value input by the user"
		for (I = 2; I <= EGS4.$MXREG; I++) {
			EGS4.IRAYLR[I - 1] = EGS4.IRAYLR[0];
		}

		// "Relaxations switch, must be done before HATCH so that the necessary"
		// "additional data can be read in in HATCH if the user requsted relaxations"
		// LOOP [
		// OUTPUT; (/' Atomic relaxations on (1) or off (0)? ',$);
		// INPUT iedgfl(1); (I5);
		// ] UNTIL ( iedgfl(1) = 1 | iedgfl(1) = 0 );
		EGS4.iedgfl[0] = 1;
		// "Now set iedgfl for all regions to the value input by the user"
		for (I = 2; I <= EGS4.$MXREG; I++) {
			EGS4.iedgfl[I - 1] = EGS4.iedgfl[0];
		}

		// "Photo-electron angular distribution switch. It does not need to be "
		// "before HATCH, we do it here because this is the most logical place"
		// LOOP [
		// OUTPUT; (/' Photo-electron angular distribution on (1) or off (0)?
		// ',$);
		// INPUT iphter(1); (I5);
		// ] UNTIL ( iphter(1) = 1 | iphter(1) = 0 );
		EGS4.iphter[0] = 0;
		// "Now set iphter for all regions to the value input by the user"
		for (I = 2; I <= EGS4.$MXREG; I++) {
			EGS4.iphter[I - 1] = EGS4.iphter[0];
		}

		// "Bound Compton scattering switch"
		// "Must be done before HATCH in order to get the data for bound Compton"
		// LOOP [
		// OUTPUT; (/' Binding effects for Compton scattering on (1) or off (0)?
		// ',$);
		// INPUT ibcmp(1); (I5);
		// ] UNTIL ( ibcmp(1) = 1 | ibcmp(1) = 0 );
		EGS4.ibcmp[0] = 1;
		// "Now set ibcmp for all regions to the value input by the user"
		for (I = 2; I <= EGS4.$MXREG; I++) {
			EGS4.ibcmp[I - 1] = EGS4.ibcmp[0];
		}

		// "Pair production angular distribution switch"
		// "Must be done before HATCH in order to get the material composition needed"
		// "for some initializations"
		// OUTPUT; (/' Pair angular distribution: the following choices are
		// available: '/
		// ' 0: fixed pair angle (EGS4 default) '/
		// ' 1: leading term of the distribution '/
		// ' 2: Koch and Motz ');
		// OUTPUT; (' your choice: ',$);
		// INPUT iprdst; (I5);
		// IF( iprdst < 0 | iprdst > 2 ) [
		// OUTPUT; (' wrong input, iprdst set to 1');
		// iprdst = 1;
		// ]
		EGS4.iprdst = 1;
		if ((EGS4.iprdst < 0) || (EGS4.iprdst > 2)) {
			// OUTPUT; (' wrong input, iprdst set to 1');
			EGS4.iprdst = 1;
		}

		// "Bremsstrahlung angular distribution switch"
		// "Must be done before HATCH in order to get the material composition needed"
		// "for some initializations"
		// OUTPUT; (/' Bremsstrahlung angular distribution,'/
		// ' 0: leading term of Koch and Motz distn'/
		// ' 1: Koch and Motz 2BS(modified): '/
		// ' your choice: ',$);
		// INPUT ibrdst; (I5);
		// IF( ibrdst < 0 | ibrdst > 1 ) [
		// OUTPUT; (' wrong input, ibrdst set to 1');
		// ibrdst = 1;
		// ]
		EGS4.ibrdst = 1;
		if ((EGS4.ibrdst < 0) || (EGS4.ibrdst > 1)) {
			// OUTPUT; (' wrong input, ibrdst set to 1');
			EGS4.ibrdst = 1;
		}

		// "Bremsstrahlung photon differential cross secton switch"
		// "Must be selected before call to HATCH if non-zero value being used"
		// OUTPUT;(/' Bremsstrahlung differential photon cross section to
		// sample,'/
		// ' 0: use Bethe-Heitler distribution as in EGS4'/
		// ' 1: use NIST/ICRU 37 distributions'/
		// ' your choice: ',$)
		// INPUT ibr_nist;(I5);
		// IF( ibr_nist <0 | ibr_nist > 1) [
		// OUTPUT; (' wrong input, ibr_nist set to 0');
		// ibr_nist = 1;
		// ]
		EGS4.ibr_nist = 1;// 0;//1;
		if ((EGS4.ibr_nist < 0) || (EGS4.ibr_nist > 1)) {
			// OUTPUT; (' wrong input, ibr_nist set to 0');
			EGS4.ibr_nist = 0;
		}

		// "Spin effects for electron/positron elastic scattering"
		// "Must be done before hatch in order to get the data for the spin "
		// "rejection loop"
		// LOOP [
		// OUTPUT; (/' Spin effects on (1) or off (0)? ',$);
		// INPUT ispin; (I5);
		// ] UNTIL ( ispin = 1 | ispin = 0 );
		// IF( ispin = 0 ) [ spin_effects = .false.; ]
		// ELSE [ spin_effects = .true.; ]
		// 1 ,ISPIN SO:
		ispin = 1;
		if (ispin == 0) {
			EGS4.spin_effects = false;
		} else {
			EGS4.spin_effects = true;
		}

		// "estepe: maximum fractional energy loss per step"
		// "Must be done before HATCH so that step-lengths are properly initialized"
		// OUTPUT; (/' Input maximum fractional energy loss per step (estepe):
		// ',$);
		// INPUT estepe; (F10.0);
		// IF( estepe <= 0 | estepe >= 1 ) [
		// estepe = $MAX-ELOSS; "$MAX-ELOSS is defined in egsnrc.macros at 0.25"
		// OUTPUT estepe; (' using default value: ',f10.3);
		// ]
		// 0 ,ESTEPE so
		EGS4.estepe = 0;
		if ((EGS4.estepe <= 0) || (EGS4.estepe >= 1)) {
			EGS4.estepe = EGS4.$MAX_ELOSS; // "$MAX-ELOSS is defined in egsnrc.macros at 0.25"
			// OUTPUT estepe; (' using default value: ',f10.3);
		}

		// "ximax: maximum first elastic scattering moment per step"
		// "used to determine step-length together with estepe"
		// "maximum possible value restricted by the maximum value for which MS data"
		// "is available"
		// "Must be done before HATCH so that step-lengths are properly initialized"
		// OUTPUT; (/' Input maximum 1st elastic scattering moment per step:
		// ',$);
		// INPUT ximax; (F10.0);
		// IF( ximax <= 0 | ximax >= 1 ) [
		// ximax = $EXACT-BCA-XIMAX;
		// "$EXACT-BCA-XIMAX set to 0.5 in egsnrc.macros"
		// OUTPUT ximax; (' using default, ximax = ',f10.3);
		// ]
		// 0 ,XIMAX so:
		EGS4.ximax = 0;
		if ((EGS4.ximax <= 0) || (EGS4.ximax >= 1)) {
			EGS4.ximax = EGS4.$EXACT_BCA_XIMAX; // "$EXACT-BCA-XIMAX set to 0.5 in egsnrc.macros"
			// OUTPUT ximax; (' using default, ximax = ',f10.3);
		}

		// "Transport algorithm"
		// OUTPUT; (/' Electron-step algorithm: EGSnrc default (0) or PRESTA
		// (1)? ',$);
		// INPUT transport_algorithm; (I5);
		// IF( transport_algorithm < 0 | transport_algorithm > 1 ) [
		// transport_algorithm = 0;
		// OUTPUT; (' using EGSnrc default');
		// ]
		// 0 ,transport_algorithm so
		EGS4.transport_algorithm = 0;
		if ((EGS4.transport_algorithm < 0) || (EGS4.transport_algorithm > 1)) {
			EGS4.transport_algorithm = 0;
			// OUTPUT; (' using EGSnrc default');
		}

		// "Boundary crossing algorithm"
		// OUTPUT; (/' Boundary crossing algorithm: exact (0) or PRESTA (1)?
		// ',$);
		// INPUT bca_algorithm; (I5);
		// IF( bca_algorithm < 0 | bca_algorithm > 1 ) [
		// bca_algorithm = 0;
		// OUTPUT; (' using exact boundary crossing ');
		// ]
		// 0 ,bca_algorithm so
		EGS4.bca_algorithm = 0;
		if ((EGS4.bca_algorithm < 0) || (EGS4.bca_algorithm > 1)) {
			EGS4.bca_algorithm = 0;
			// OUTPUT; (' using exact boundary crossing ');
		}

		// "skin-depth for BCA"
		// OUTPUT; (/' Skin-depth for BCA: ');
		// IF( bca_algorithm = 0 ) [
		// OUTPUT; (' this is the distance from a boundary ' /
		// ' (measured in elastic mean-free-paths) at which the '/
		// ' simulation switches to single scattering mode '/
		// ' Best choice for efficiency is 3. If you set this '/
		// ' parameter to a very large number (e.g. 1e10), you '/
		// ' can force single scattering simulation in the entire '/
		// ' geometry (this is very slow) '/
		// ' your choice: ',$);
		// INPUT skindepth_for_bca; (F10.0);
		// IF( skindepth_for_bca <= 0 ) skindepth_for_bca = 3;
		// ]
		// ELSE [
		// OUTPUT; (' this is the distance from a boundary '/
		// ' (measured in elastic mean-free-paths) at which '/
		// ' lateral deflections will be turned off. If you select '/
		// ' a very large number (e.g. 1e10), standard EGS4 '/
		// ' behaviour (no PRESTA) will result. If you input '/
		// ' a number < 1, this parameter will be determined in the'/
		// ' way it was with PRESTA (i.e. depending on ECUT)'/
		// ' your choice: ',$);
		// INPUT skindepth_for_bca; (F10.0);
		// ]
		// our case:
		EGS4.skindepth_for_bca = 3;
		if (EGS4.bca_algorithm == 0) {
			// INPUT skindepth_for_bca; (F10.0);
			if (EGS4.skindepth_for_bca <= 0)
				EGS4.skindepth_for_bca = 3;
		}

		// LOOP [
		// OUTPUT;
		// (/' Input random number luxury level (0-4) & seed (>0) (0 defaults
		// OK): ',$);
		// INPUT luxury_level,iseed; (2I5);
		// ] UNTIL ((luxury_level <= 4 & luxury_level >= 0) & (iseed >= 0));
		// $INITIALIZE RNG USING luxury_level AND iseed;
		// REPLACE {$INITIALIZE RNG USING # AND #; } WITH
		// {
		// call init_ranlux({P1},{P2});
		// call ranlux(rng_array); rng_seed = 1;
		// }
		int llev = 1;
		int isd = 1;
		EGS4.init_ranlux(llev, isd);
		EGS4.ranlux(EGS4.rng_array);
		EGS4.rng_seed = 1;

		// " At this point we could perform step 3 (call to HATCH), "
		// " we will defer the HATCH call to be performed in the main routine"
		// "---------------------------------------------------------------------"
		// "STEP 4  INITIALIZATION-FOR-HOWFAR and HOWNEAR                        "
		// "---------------------------------------------------------------------"

		// "In this simplified version of the code we will allow for a single    "
		// "slab of material only. (we need 2 planes to define 1 region)         "
		nzb = 2;

		// "and set region 1 and 3 to vacuum"
		EGS4.MED[0] = 0;
		EGS4.MED[2] = 0;// med(1) = 0; med(3) = 0;

		zbound[0] = 0.0;// zbound(1) = 0;
		// OUTPUT; (/' Input slab thickness: ',$);
		// INPUT zbound(2); (F10.0);
		zbound[1] = 0.1;

		// "---------------------------------------------------------------------"
		// "STEP 5  INITIALIZATION-FOR-AUSGAB                                    "
		// "---------------------------------------------------------------------"

		// "Set all scoring arrays to zero. This could be avoided if"
		// "the compiler being used has a `initialize to zero' option"
		// "It is a good coding habit to not rely on variables being"
		// "automatically zeroed"

		for (I = 1; I <= EGS4.$MXREG; I++) {
			sc_array[I - 1] = 0.0;
			sc_array2[I - 1] = 0.0;
			sc_tmp[I - 1] = 0.0;
			sc_last[I - 1] = 0;
		}

		// "Define range rejection parameter. Although not directly related "
		// "to ausgab, range rejection is an `user' variance reduction technique "
		// "and so, this is the most appropriate place to initialize it"

		// "---------------------------------------------------------------------"
		// "STEP 5b  INITIALIZATION-FOR-Variance-Reduction                       "
		// "---------------------------------------------------------------------"

		// OUTPUT; (/' Use (1) or do not use (0) electron range rejection? ',$);
		// INPUT irejct; (I5);
		irejct = 1;
		if (irejct == 1) {
			for (I = 1; I <= EGS4.$MXREG; I++) {
				EGS4.i_do_rr[I - 1] = 1;
			} // "initialize for all regions"
			// OUTPUT; (/' Input the maximum energy to apply range rejection:
			// ',$);
			// INPUT esave; (F10.0);
			esave = 5.0;
			for (I = 1; I <= EGS4.$MXREG; I++) {
				EGS4.e_max_rr[I - 1] = esave;
			} // "initialize for all regions"
		}

		// OUTPUT; (/' How many brem photons to create per event (0=>just 1):
		// ',$);
		// INPUT nbr_split;(I10);
		EGS4.nbr_split = 1;
		if (EGS4.nbr_split <= 0) {
			// OUTPUT; (' Negative or zero value of nbr_split made 1=> no
			// splitting');
			EGS4.nbr_split = 1;
		}

		// OUTPUT;
		// (/' Russian Roulette all secondary charged particles(yes=1,no=0)?
		// ',$);
		// INPUT i_play_RR;(I3);
		EGS4.i_play_RR = 0;
		if (EGS4.i_play_RR != 0) {
			EGS4.i_play_RR = 1;
			Integer itg = new Integer(EGS4.nbr_split);
			EGS4.prob_RR = 1. / itg.doubleValue();// "We are assuming the purpose is to"
			// "keep a natural number of charged"
			// "particles, even when using splitting"
		}
		// "---------------------------------------------------------------------"
		// "STEP 6   DETERMINATION-OF-INICIDENT-PARTICLE-PARAMETERS              "
		// "---------------------------------------------------------------------"
		// LOOP [
		// OUTPUT; (/' Input incident charge(-1, 0, +1): ',$);
		// INPUT iqin; (I5);
		// ] UNTIL ( iqin >= -1 & iqin <= 1 );
		IQIN = -1;

		// OUTPUT; (/' Input incident kinetic energy (MeV): ',$);
		// INPUT ein; (F10.0);
		EIN = 5.0;
		if (IQIN != 0)
			EIN = EIN + EGS4.RM; // "add rest energy for electrons and positrons"

		// LOOP [
		// OUTPUT; (/' Input incident angle (degrees): ',$);
		// INPUT win; (F10.0);
		// ] UNTIL ( win >= 0 & win <= 90 );
		WIN = 0.0;
		// "Convert to direction cosine and calculate other 2 direction cosines"
		WIN = WIN / 180 * Math.PI;
		WIN = Math.cos(WIN);
		// choosing lets say VIN=0, then uin is equal with sin!!
		UIN = Math.sqrt(Math.max(0.0, (1 - WIN) * (1 + WIN)));
		VIN = 0.0;
		IRIN = 2; // "starts in region 2, could be 1"
		WTIN = 1; // "statistical weight is 1"

		// OUTPUT; (/' Input number of showers to be simulated: ',$);
		// INPUT ncase; (I12);
		NCASE = 1000;

		return;
		// end; "end of subroutine inputs for tutor6.mortran"

	}
}
