package danfulea.phys.egs;

/**
 * 
 * Utility class for handling the effect of external electromagnetic field (EMF) upon electron interaction with matter. 
 * Based on EGS4/EGSnrc routines developed by SLAC(Stanford Linear Accelerator Center)/NRC (National Research Council of Canada).
 * @author Dan Fulea, 08 OCT. 2005
 */

public class EGS4emf {
	public static double Ex0 = 0.0;
	public static double Ey0 = 0.0;
	public static double Ez0 = 0.0;
	public static double Bx0 = 0.0;
	public static double By0 = 0.0;
	public static double Bz0 = 0.0;
	public static double efx0 = 0.0;
	public static double efy0 = 0.0;
	public static double efz0 = 0.0;
	public static double bfx0 = 0.0;
	public static double bfy0 = 0.0;
	public static double bfz0 = 0.0;
	public static double fbtemp = 0.0;
	public static double pot1 = 0.0;
	public static double pot2 = 0.0;

	/**
	 * Configure EMF. See default EMF.
	 * @param e e
	 * @param b b
	 * @param x0 x0
	 * @param y0 y0
	 * @param z0 z0
	 * @param iemf iemf
	 */
	public static void GET_EM_FIELD(String e, String b, double x0, double y0,
			double z0, int iemf) {
		if (iemf == 0)// default example
		{
			defaultEMF(e, b, x0, y0, z0);
		} else if (iemf == 1) {
			// further implementation
		}
	}

	/**
	 * Get the potential for electric field. See defaultPot.
	 * @param P P
	 * @param X0 X0
	 * @param Y0 Y0
	 * @param Z0 Z)
	 * @param iemf iemf
	 */
	public static void GET_POTENTIAL(String P, double X0, double Y0, double Z0,
			int iemf) {
		if (iemf == 0)// default example
		{
			defaultPot(P, X0, Y0, Z0);
		} else if (iemf == 1) {
			// further implementation
		}
	}

	// "*****************************************************************************"
	// "                                                                             "
	// "THIS IS THE DEFAULT FIELD CONFIGURATION. THE ELECTRIC FIELD IS A CONSTANT IN "
	// "THE Y-DIRECTION AND THE MAGNETIC FIELD IS ZERO.                              "
	// "THE USER MUST SPECIFY THE ELECTRIC AND MAGNETIC FIELD AT P3=X,P4=Y,P5=Z      "
	// "IN THE UNITS 1/(cm).                                                         "
	// "THIS IS DONE BY MULTIPLYING THE ELECTRIC FIELD IN volts/cm BY THE UNIT CHARGE"
	// "e(positive) AND DIVIDING BY THE REST MASS OF THE ELECTRON IN electron-volts, "
	// "AND BY MULTIPLYING THE MAGNETIC FIELD IN volts-sec BY THE UNIT CHARGE        "
	// "e(positive) AND THE SPEED OF LIGHT IN cm/sec AND DIVIDING BY THE REST MASS   "
	// "OF THE ELECTRON IN electron-volts.                                           "
	// "FOR EXAMPLE, FOR A CONSTANT E-FIELD, P1X0=0.0;P1Y0=1.0;P1Z0=0.0;             "
	// "                                     P2X0=0.0;P2Y0=0.0;P2Z0=0.0;             "
	// "             FOR A CONSTANT B-FIELD, P1X0=0.0;P1Y0=0.0;P1Z0=0.0;             "
	// "                                     P2X0=0.0;P2Y0=0.0;P2Z0=1.0;             "
	// "             FOR A POINT SOURCE      P1X0=P3/(SQRT(P3**2+P4**2+P5**2))**3;   "
	// "                                     P1Y0=P4/(SQRT(P3**2+P4**2+P5**2))**3;   "
	// "                                     P1Z0=P5/(SQRT(P3**2+P4**2+P5**2))**3;   "
	// "                                     P2X0=0.0;P2Y0=0.0;P2Z0=1.0;             "
	// "THE FOLLOWING EXAMPLE IS FOR A CONDUCTING CYLINDER WITH RADIUS 0.5 CM        "
	// //(E,B,X0,Y0,Z0) or (EF,BF,XF,YF,ZF)
	/**
	 * THIS IS THE DEFAULT FIELD CONFIGURATION. THE ELECTRIC FIELD IS A CONSTANT IN THE Y-DIRECTION AND THE MAGNETIC FIELD IS ZERO. 
	 * THE USER MUST SPECIFY THE ELECTRIC AND MAGNETIC FIELD AT P3=X,P4=Y,P5=Z IN THE UNITS 1/(cm). 
	 * THIS IS DONE BY MULTIPLYING THE ELECTRIC FIELD IN volts/cm BY THE UNIT CHARGE e(positive) AND DIVIDING BY THE REST MASS OF 
	 * THE ELECTRON IN electron-volts, AND BY MULTIPLYING THE MAGNETIC FIELD IN volts-sec BY THE UNIT CHARGE 
	 * e(positive) AND THE SPEED OF LIGHT IN cm/sec AND DIVIDING BY THE REST MASS OF THE ELECTRON IN electron-volts. <p>
	 * FOR EXAMPLE, FOR A CONSTANT E-FIELD, P1X0=0.0;P1Y0=1.0;P1Z0=0.0; P2X0=0.0;P2Y0=0.0;P2Z0=0.0; <p>
	 * FOR A CONSTANT B-FIELD, P1X0=0.0;P1Y0=0.0;P1Z0=0.0; P2X0=0.0;P2Y0=0.0;P2Z0=1.0; <p>
	 * FOR A POINT SOURCE      P1X0=P3/(SQRT(P3**2+P4**2+P5**2))**3; P1Y0=P4/(SQRT(P3**2+P4**2+P5**2))**3; P1Z0=P5/(SQRT(P3**2+P4**2+P5**2))**3; 
	 * P2X0=0.0;P2Y0=0.0;P2Z0=1.0;  <p>
	 * THE FOLLOWING EXAMPLE IS FOR A CONDUCTING CYLINDER WITH RADIUS 0.5 CM (E,B,X0,Y0,Z0) or (EF,BF,XF,YF,ZF).
	 * @param e E or ef
	 * @param b B or bf
	 * @param x0 x0
	 * @param y0 y0
	 * @param z0 z0
	 */
	private static void defaultEMF(String e, String b, double x0, double y0,
			double z0) {
		// REPLACE {$GET-EM-FIELD(#,#,#,#,#);} WITH {;
		fbtemp = y0 * y0 + z0 * z0;// fbtemp={P4}**2+{P5}**2;
		if (fbtemp <= 0.25) {
			if (e.compareTo("E") == 0) {
				Ex0 = 0.0;
				Ey0 = 0.0;
				Ez0 = 0.0;
			} else if (e.compareTo("ef") == 0) {
				efx0 = 0.0;
				efy0 = 0.0;
				efz0 = 0.0;
			}
		} else {
			fbtemp = 12.43 / fbtemp;
			if (e.compareTo("E") == 0) {
				Ex0 = 0.0;
				Ey0 = y0 * fbtemp;// {P4}*fbtemp;
				Ez0 = z0 * fbtemp;// {P5}*fbtemp;
			} else if (e.compareTo("ef") == 0) {
				efx0 = 0.0;
				efy0 = y0 * fbtemp;// {P4}*fbtemp;
				efz0 = z0 * fbtemp;// {P5}*fbtemp;
			}
		}
		if (b.compareTo("B") == 0) {
			Bx0 = 0.0;
			By0 = 0.0;
			Bz0 = 0.0;
		} else if (b.compareTo("bf") == 0) {
			bfx0 = 0.0;
			bfy0 = 0.0;
			bfz0 = 0.0;
		}
		// }
	}

	// "*****************************************************************************"
	// "                                                                             "
	// "THIS MACRO RETURNS THE POTENTIAL FOR THE ELECTRIC FIELD                      "
	// "THE USER MUST RETURN THE POTENTIAL P1 IN MeV AT X=P2,Y=P3,Z=P4.              "
	// "FOR EXAMPLE, FOR A CONSTANT E-FIELD, P1=-P3*IQ(NP)*RM;                       "
	// "             FOR A CONSTANT B-FIELD, P1=0.0;                                 "
	// "             FOR A POINT SOURCE P1=-IQ(NP)*RM/SQRT(P2**2+P3**2+P4**2);       "
	// "THE FOLLOWING EXAMPLE IS FOR A CONDUCTING CYLINDER WITH RADIUS 0.5 CM        "
	/**
	 * RETURNS THE POTENTIAL FOR THE ELECTRIC FIELD. THE USER MUST RETURN THE POTENTIAL P1 IN MeV AT X=P2,Y=P3,Z=P4. <p>
	 * FOR EXAMPLE, FOR A CONSTANT E-FIELD, P1=-P3*IQ(NP)*RM; <p>
	 * FOR A CONSTANT B-FIELD, P1=0.0; <p>
	 * FOR A POINT SOURCE P1=-IQ(NP)*RM/SQRT(P2**2+P3**2+P4**2); <p>
	 * THE FOLLOWING EXAMPLE IS FOR A CONDUCTING CYLINDER WITH RADIUS 0.5 CM 
	 * @param P pot1 or pot2
	 * @param X0 X0
	 * @param Y0 Y0
	 * @param Z0 Z0
	 */
	private static void defaultPot(String P, double X0, double Y0, double Z0) {
		// REPLACE {$GET-POTENTIAL(#,#,#,#);} WITH {;
		fbtemp = Y0 * Y0 + Z0 * Z0;
		if (fbtemp <= 0.25) {
			if (P.compareTo("pot1") == 0)
				pot1 = 0.0;
			else if (P.compareTo("pot2") == 0)
				pot2 = 0.0;
		} else {
			if (P.compareTo("pot1") == 0)
				pot1 = -EGS4.IQ[EGS4.NP - 1] * EGS4.RM * 6.215
						* Math.log(4. * (Y0 * Y0 + Z0 * Z0));
			else if (P.compareTo("pot2") == 0)
				pot2 = -EGS4.IQ[EGS4.NP - 1] * EGS4.RM * 6.215
						* Math.log(4. * (Y0 * Y0 + Z0 * Z0));
		}
		// }

	}
}
