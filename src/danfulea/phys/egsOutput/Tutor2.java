package danfulea.phys.egsOutput;

import java.util.Calendar;
import java.util.Date;

import danfulea.phys.egs.EGS4;
import danfulea.phys.egs.EGS4Core;
import danfulea.phys.egs.EgsQuestion;

/**
 * Demo class
 * Based on EGS4/EGSnrc routines developed by SLAC(Stanford Linear Accelerator Center)/NRC (National Research Council of Canada).
 * @author Dan Fulea, 02 OCT. 2005
 */

// "  An EGSnrc user code. It scores total amount of energy reflected,  "
// " deposited and transmitted when a 20 MeV beam of electrons is        "
// " incident on a 1mm slab of Ta.

public class Tutor2 implements EgsQuestion {
	private int I = 0;
	private int IQIN = 0;
	private int IRIN = 0;
	private double XIN = 0.;
	private double YIN = 0.;
	private double ZIN = 0.;
	private double EIN = 0.;
	private double WTIN = 0.;
	private double UIN = 0.;
	private double VIN = 0.;
	private double WIN = 0.;
	private double[] ESCORE = new double[3];
	private int NCASE = 0;
	private double ANORM = 0.0;
	private double TOTAL = 0.0;

	private double ZBOUND = 0.1;// thickness of medium: e.g 0.1 cm of Tl plate

	public Tutor2() {
		init();
	}

	private void init() {
		EGS4.startSimulationTime = System.currentTimeMillis();
		// "---------------------------------------------------------------------"
		// "@@@@@STEP 1: USER-OVERRIDE-OF-EGSnrc-MACROS
		// "---------------------------------------------------------------------"
		EGS4.setMXMED(1);// "only 1 medium in the problem(default 10)"-->2 for 2
							// media
		EGS4.setMXREG(3);// "only 3 geometric regions (default 2000)"-->4 for 4
							// regions
		EGS4.setMXSTACK(15);// "less than 15 particles on stack at once"
		// "---------------------------------------------------------------------"
		// "@@@@@STEP 2 PRE-HATCH-CALL-INITIALIZATION
		// "---------------------------------------------------------------------"
		EGS4.eq = this;
		EGS4Core.eq = this;
		EGS4.egs_set_defaults();
		EGS4.MEDIA[0] = "TA_fortran";// only 1 medium and pegs4file is
										// TA_fortran.pegs4dat!!!!!!!!!!!!!!!!
		// EGS4.MEDIA[0]="tantalum";
		EGS4.MED[0] = 0;
		EGS4.MED[2] = 0;
		EGS4.MED[1] = 1;
		EGS4.ECUT[1] = 1.5;
		EGS4.PCUT[1] = 0.1;
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
		ZBOUND = 0.1;// "     plate is 1 mm thick"
		// "---------------------------------------------------------------------"
		// "@@@@@STEP 5  INITIALIZATION-FOR-AUSGAB                                    "
		// "---------------------------------------------------------------------"
		for (I = 1; I <= 3; I++) {
			ESCORE[I - 1] = 0.0;
		}// "zero scoring array before starting"
		// "---------------------------------------------------------------------"
		// "@@@@@STEP 6   DETERMINATION-OF-INICIDENT-PARTICLE-PARAMETERS              "
		// "---------------------------------------------------------------------"
		// "Define initial variables for 20 MeV beam of electrons incident"
		// "perpendicular to the slab"
		IQIN = -1;// "               incident charge - electrons"
		EIN = 20.511;// "            20 MeV kinetic energy"
		XIN = 0.0;
		YIN = 0.0;
		ZIN = 0.0;// "     incident at origin"
		UIN = 0.0;
		VIN = 0.0;
		WIN = 1.0;// " moving along Z axis"
		IRIN = 2;// "                starts in region 2, could be 1"
		WTIN = 1.0;// "              weight = 1 since no variance reduction used"
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
		ANORM = 100. / ((EIN + IQIN * 0.511) * NCASE);
		// "normalize to % of total input energy"
		TOTAL = 0.0;
		for (I = 1; I <= 3; I++) {
			TOTAL = TOTAL + ESCORE[I - 1];
		}

		EGS4.seqStr = " ********************************************** "
				+ " \n";
		if (EGS4.iprint > 0)
			printSequence(EGS4.seqStr);
		double f = ANORM * ESCORE[0];
		EGS4.seqStr = " Fraction of energy reflected from plate= "
				+ EGS4.format(f, 10, true) + " %" + " \n";
		if (EGS4.iprint > 0)
			printSequence(EGS4.seqStr);
		f = ANORM * ESCORE[1];
		EGS4.seqStr = " Fraction of energy deposited in plate= "
				+ EGS4.format(f, 10, true) + " %" + " \n";
		if (EGS4.iprint > 0)
			printSequence(EGS4.seqStr);
		f = ANORM * ESCORE[2];
		EGS4.seqStr = " Fraction of energy transmitted through plate= "
				+ EGS4.format(f, 10, true) + " %" + " \n";
		if (EGS4.iprint > 0)
			printSequence(EGS4.seqStr);
		f = TOTAL * ANORM;
		EGS4.seqStr = " Total fraction of energy accounted for= "
				+ EGS4.format(f, 10, true) + " \n";
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
		// output data to console
		System.out.println(s);
	}

	public void AUSGAB(int IARG) {
		// "  In general, AUSGAB is a routine which is called under a series    "
		// "  of well defined conditions specified by the value of IARG (see the"
		// "  EGSnrc manual for the list). This is a particularily simple AUSGAB"
		// "  Whenever this routine is called with IARG=3 , a particle has      "
		// "  been discarded by the user in HOWFAR                              "
		// "  We get AUSGAB to print the required information at that point     "
		// "********************************************************************"
		// ;Copyright NRC;
		// "we use EDEP from EPCONT,IR(NP) from STACK"
		// "ESCORE is passed in user defined COMIN SCORE"
		if (IARG <= 4) {
			int IRL = EGS4.IR[EGS4.NP - 1];// " pick up current region number"
			ESCORE[IRL - 1] = ESCORE[IRL - 1] + EGS4.EDEP;
		}
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
	//	double y = EGS4.Y[EGS4.NP - 1];
	//	double x = EGS4.X[EGS4.NP - 1];
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
