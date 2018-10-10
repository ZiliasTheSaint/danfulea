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

// " An EGSnrc user code which scores the number and average energy of   "
// " primary, Rayleigh scattered and Compton scattered photons passing   "
// " through a  5 mm thick slab of water when a 50 keV pencil beam of    "
// " photons is incident normally                                        "

public class Tutor5 implements EgsQuestion {
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
	private double ANORM = 0.;
	private double[] COUNT = new double[3];
	private double[] ENTOT = new double[3];
	// "  in COUNT(1),(2),(3) AUSGAB will count the number of transmitted  "
	// "  primaries, rayleigh scattered or only compton scattered photons  "
	// "  ENTOT adds up  the total energy of each of these components      "

	private double ZBOUND = 0.0;// thickness of medium: e.g 0.1 cm of Tl plate

	public Tutor5() {
		init();
	}

	private void init() {
		EGS4.startSimulationTime = System.currentTimeMillis();// grab start time
		// "---------------------------------------------------------------------"
		// "@@@@@STEP 1: USER-OVERRIDE-OF-EGSnrc-MACROS
		// "---------------------------------------------------------------------"
		EGS4.setMXMED(1);// "only 1 medium in the problem(default 10)"-->2 for 2
							// media
		EGS4.setMXREG(3);// "only 3 geometric regions (default 2000)"-->4 for 4
							// regions
		EGS4.setMXSTACK(15);// "less than 15 particles on stack at once"
		// REPLACE {$EBIN} WITH {25} "user parameter -# bins in energy spectrum"
		// REPLACE {$CALL-HOWNEAR(#);} WITH {
		// ;CALL HOWNEAR({P1},X(NP),Y(NP),Z(NP),IRL);
		// "---------------------------------------------------------------------"
		// "@@@@@STEP 2 PRE-HATCH-CALL-INITIALIZATION
		// "---------------------------------------------------------------------"
		EGS4.egs_set_defaults();// EGS4.RandomUse=2;
		EGS4.eq = this;
		EGS4Core.eq = this;
		// EGS4.MEDIA[0]="water_liquid";//"TA_fortran";//only 1 medium and
		// pegs4file is TA_fortran.pegs4dat!!!!!!!!!!!!!!!!
		EGS4.MEDIA[0] = "H2O_fortran";// "sodiumiodide";//"NAI_fortran";
		// /MED(1),MED(3)/=0;MED(2)=1;"vacuum in regions 1 and 3, H2O in region 2"
		EGS4.MED[0] = 0;
		EGS4.MED[2] = 0;
		EGS4.MED[1] = 1;
		EGS4.ECUT[1] = 1.5;// "   terminate electron histories at 1.5 MeV in the plate"
		EGS4.PCUT[1] = 0.010;// "   terminate   photon histories at 0.01 MeV in the slab"
		// "               only needed for region 2 since no transport elsewhere"
		// "               ECUT is total energy = 0.989   MeV kinetic energy"

		EGS4.IRAYLR[1] = 1;// "     turn on Rayleigh scattering in the slab "
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
		// "@@@@@STEP 4  INITIALIZATION-FOR-HOWFAR and HOWNEAR                        "
		// "---------------------------------------------------------------------"
		ZBOUND = 0.5;// "     plate is 5 mm thick"
		// "---------------------------------------------------------------------"
		// "@@@@@STEP 5  INITIALIZATION-FOR-AUSGAB                                    "
		// "---------------------------------------------------------------------"
		for (I = 1; I <= 3; I++) {
			COUNT[I - 1] = 0.0;
			ENTOT[I - 1] = 0.0;
		}// "zero scoring array before starting"

		// "  We want to set flags in AUSGAB every time a Rayleigh scattering   "
		// "  or compton scattering occurs. Set the flags in IAUSFL(COMIN      "
		// "  EPCONT) to signal the  EGS system to make the appropriate calls  "
		EGS4.iausfl[17] = 1;// [18]=1;
		EGS4.iausfl[23] = 1;// [24]=1;
		// "---------------------------------------------------------------------"
		// "@@@@@STEP 6   DETERMINATION-OF-INICIDENT-PARTICLE-PARAMETERS              "
		// "---------------------------------------------------------------------"
		// "Define initial variables for 5 MeV beam of photons incident"
		// "perpendicular to the slab"
		IQIN = 0;// "               incident charge - electrons"
		EIN = 0.050;// "           MeV kinetic energy"
		XIN = 0.0;
		YIN = 0.0;
		ZIN = 0.0;// "     incident at origin"
		UIN = 0.0;
		VIN = 0.0;
		WIN = 1.0;// " moving along Z axis"
		IRIN = 2;// "                starts in region 2, could be 1"
		WTIN = 1.0;// "              weight = 1 since no variance reduction used"
		EGS4.LATCHI = 0;// "              LATCH set to zero at start of each history"
		// "---------------------------------------------------------------------"
		// "@@@@@STEP 7   SHOWER-CALL                                                 "
		// "---------------------------------------------------------------------"
		// "initiate the shower 10 times"
		NCASE = 10000; // "INITIATE THE SHOWER NCASE TIMES"
		for (I = 1; I <= NCASE; I++) {
			SHOWER();
			// if insucces:
			if (EGS4.STOPPROGRAM) {
				return;
			}

		}
		// "-----------------------------------------------------------------"
		// "@@@@STEP 8   OUTPUT-OF-RESULTS                                       "
		// "-----------------------------------------------------------------"
		// "note output is at the end of each history in subroutine ausgab"
		ANORM = 100. / NCASE;// (new Integer(NCASE)).doubleValue();
		// "normalize to % of total input energy"
		for (I = 1; I <= 3; I++) {
			if (COUNT[I - 1] != 0)
				ENTOT[I - 1] = ENTOT[I - 1] / COUNT[I - 1];// "get average energies"
		}
		// OUTPUT EIN*1000.,ZBOUND, PCUT(2), (ANORM*COUNT(I),ENTOT(I),I=1,3);
		// (//' For',F6.1,' keV photons incident on',F4.1,' cm OF H2O'/
		// ' with PCUT=',F5.3,' MeV'
		// //' Transmitted primaries=',T40,F8.2,'% ave energy=',F10.3,' MeV'//
		// ' Fraction rayleigh scattering=',T40,F8.2,'% ave energy=',F10.3,'
		// MeV'
		// //' Fraction compton scattering only=',T40,F8.2,'% ave
		// energy=',F10.3,
		// ' MeV'//)
		double einc = EIN * 1000.;

		EGS4.seqStr = " ********************************************** ";
		if (EGS4.iprint > 0)
			printSequence(EGS4.seqStr);

		EGS4.seqStr = " For " + einc + " keV photons incident on"
				+ EGS4.format(ZBOUND, 4, true) + " cm OF H2O" + " \n"
				+ " with PCUT=" + EGS4.format(EGS4.PCUT[1], 5, true);
		if (EGS4.iprint > 0)
			printSequence(EGS4.seqStr);

		// for(I=1;I<=3;I++)
		// {
		double f = ANORM * COUNT[0];
		EGS4.seqStr = " Transmitted primaries=" + EGS4.format(f, 8, true)
				+ "%  ave energy=" + EGS4.format(ENTOT[0], 8, true) + " MeV";
		if (EGS4.iprint > 0)
			printSequence(EGS4.seqStr);

		f = ANORM * COUNT[1];
		EGS4.seqStr = " Fraction rayleigh scattering="
				+ EGS4.format(f, 8, true) + "%  ave energy="
				+ EGS4.format(ENTOT[1], 8, true) + " MeV";
		if (EGS4.iprint > 0)
			printSequence(EGS4.seqStr);

		f = ANORM * COUNT[2];
		EGS4.seqStr = " Fraction compton scattering only="
				+ EGS4.format(f, 8, true) + "%  ave energy="
				+ EGS4.format(ENTOT[2], 8, true) + " MeV";
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

	public void AUSGAB(int IARG) {
		// " In this AUSGAB routine for tutor5 we both set flags whenever there is"
		// " a scattering event and then count histories when they have come      "
		// " through the slab , according to what kind of scattering they have    "
		// " undergone.                                                           "
		// " The logic is as follows                                              "
		// " set FLAG1 if a compton event occurs                                  "
		// " set FLAG2 if a Rayleigh event occurs                                 "
		// " The FLAGS are the units and thousands digits in the parameter LATCH  "
		// "                                                                      "
		// " When a history is terminated, increment various counters according   "
		// " to whether no flags are set - i.e. its a primary, FLAG2 is set,      "
		// " i.e. it has Rayleigh scattered or FLAG1 is set and FLAG2 is not set  "
		// /" i.e. only  compton scattering has occured.                           "
		// "                                                                      "
		// " First a few macros are defined to make the logic simpler to read and "
		// " therefore verify                                                     "
		// "                                                                      "
		// "********************************************************************"
		// ;Copyright NRC;
		int JJ = 0;
		// REPLACE {$SET-FLAG1;} WITH {LATCH(NP)=LATCH(NP)+1;}
		// REPLACE {$SET-FLAG2;} WITH {LATCH(NP)=LATCH(NP)+1000;}
		// REPLACE {$FLAG1} WITH {MOD(LATCH(NP),100)}"i.e. units digit of LATCH"
		// REPLACE {$FLAG2} WITH {MOD(LATCH(NP),10000)-$FLAG1} "thousands digit"
		int $FLAG1 = EGS4.LATCH[EGS4.NP - 1] % 100;
		int $FLAG2 = EGS4.LATCH[EGS4.NP - 1] % 10000 - $FLAG1;
		// ;COMIN/SCORE,STACK/; "we use IR(NP) from STACK"

		// "   first set flags when scattering events occur - IAUSFL was set "
		// "   in step 5 of main to ensure AUSGAB was called at these points "

		if (IARG == 17) { // "a compton scatter is about to occur"
							// $SET-FLAG1;
			EGS4.LATCH[EGS4.NP - 1] = EGS4.LATCH[EGS4.NP - 1] + 1;
		} else if (IARG == 23) {// "a rayleigh scatter is about to occur"
								// $SET-FLAG2;
			EGS4.LATCH[EGS4.NP - 1] = EGS4.LATCH[EGS4.NP - 1] + 1000;
		}
		// " if a history has terminated because leaving the slab, score it"
		else if (IARG == 3) {// "particle has left slab"
			if ((EGS4.IR[EGS4.NP - 1] == 3) || (EGS4.IR[EGS4.NP - 1] == 1)) {
				// "it is transmitted or reflected"
				JJ = 0;
				if (EGS4.LATCH[EGS4.NP - 1] == 0) {// "no scattering - a primary"
					JJ = 1;
				} else if ($FLAG2 != 0) {// "at least one rayleigh scatter"
					JJ = 2;
				} else if ($FLAG1 != 0) {// ">=1 compton scatter, no Rayleigh"
					JJ = 3;
				} else {// "debug";
				}

				if (JJ != 0) {
					COUNT[JJ - 1] = COUNT[JJ - 1] + 1.;
					ENTOT[JJ - 1] = ENTOT[JJ - 1] + EGS4.E[EGS4.NP - 1];
				}
			}// "end region 3 or 1 block"
		}// "end IARG 3 block"

	}// "END OF AUSGAB"

	// /"*********************************************************************"
	public void HOWFAR() {
		// " The following is a general specification of HOWFAR                  "
		// "   Given a particle at (X,Y,Z) in region IR and going in direction   "
		// "   (U,V,W), this routine answers the question, can the particle go   "
		// "   a distance USTEP without crossing a boundary?                     "
		// "           If yes, it merely returns                                 "
		// "           If no, it sets USTEP=distance to boundary in the current  "
		// "           direction and sets IRNEW to the region number   on the    "
		// "           far side of the boundary (this can be messy in general!)  "
		// "                                                                     "
		// "   The user can terminate a history by setting IDISC>0. Here we      "
		// "   terminate all histories which enter region 3 or are going         "
		// "   backwards in region 1                                             "
		// "                                                                     "
		// "                   |               |                                 "
		// "   REGION 1        |   REGION 2    |       REGION 3                  "
		// "                   |               |                                 "
		// "   e- =========>   |               | e- or photon ====>              "
		// "                   |               |                                 "
		// "   vacuum          |     Ta        |       vacuum                    "
		// "                                                                     "
		// "*********************************************************************"
		// ;Copyright NRC;

		double TVAL = 0.0;

		if (EGS4.IR[EGS4.NP - 1] == 3) {
			EGS4.IDISC = 1;
			return;// "terminate this history: it is past the plate"
		} else if (EGS4.IR[EGS4.NP - 1] == 2) {// "We are in the Ta plate - check the geometry"

			if (EGS4.W[EGS4.NP - 1] > 0.0) {
				// "going forward - consider first since  most frequent"
				TVAL = (ZBOUND - EGS4.Z[EGS4.NP - 1]) / EGS4.W[EGS4.NP - 1];// "TVAL is dist to boundary"
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
		if (EGS4.irl == 3) {
			return;
		} else if (EGS4.irl == 2) {// "We are in the Ta plate - check the geometry"
			EGS4.tperp = Math.min(z, (ZBOUND - z));
		} else if (EGS4.irl == 1) {
			return;
		}

	}// "end of subroutine HOWNEAR"
}
