package danfulea.phys.egs;

/**
 * Utility class for output results. 
 * Based on EGS4/EGSnrc routines developed by SLAC(Stanford Linear Accelerator Center)/NRC (National Research Council of Canada).
 * 
 * @author Dan Fulea, 22 NOV. 2005
 */
public class EGS4Grid {
	public static EgsQuestion eq;
	public static String[] CDSTBL;// =new String[EGS4.$MXREG];
	public static String[] CAVTRACK;// =new String[EGS4.$MXREG];
	public static String[] CTRTBL;// =new String[EGS4.$MXREG];
	public static String[] CABSRB;// =new String[EGS4.$MXREG];

	private static int RADNUM = 0;
	private static int HSET = 0;
	private static int FMT = 0;
	private static int $MAXRZ = EGS4Geom.$MAXZREG;// =$MAXZREG}
													// "MAX(MAXRADII,MAXZREG)"
	public static double[] RADIAL_BINS = new double[$MAXRZ];
	public static double[] DEPTH_BINS = new double[$MAXRZ];

	private static int PGTHROW = 0;
	private static int DEEPNUM = 0;
	private static boolean ROT = false;
	private static int DLYPT = 0;
	private static int REGNUM = 0;
	private static int MNUM1 = 0;
	private static int MNUM2 = 0;
	private static int MNUM3 = 0;
	private static String[] MED_NAME1 = new String[11];
	private static String[] MED_NAME2 = new String[11];
	private static String[] MED_NAME3 = new String[11];
	public static int NCOMP = 0;
	public static int $MAXCMPTS = 14;
	public static double[][][] RESULTS = new double[EGS4Geom.$MAXZREG][EGS4Geom.$MAXRADII][$MAXCMPTS];
	public static double[][][] UNCRT = new double[EGS4Geom.$MAXZREG][EGS4Geom.$MAXRADII][$MAXCMPTS];// ($MAXZREG,
																									// $MAXRADII,
																									// $MAXCMPTS),
	public static String[] LABELS = new String[$MAXCMPTS];
	private static double[][][] VALUES = new double[EGS4Geom.$MAXZREG][EGS4Geom.$MAXRADII][$MAXCMPTS];
	// next could be locals!!
	//private static boolean ESTEPSON = false;
	private static boolean ECUTON = false;
	private static boolean PCUTON = false;
	private static double[][][] TMP2 = new double[EGS4Geom.$MAXZREG][EGS4Geom.$MAXRADII][$MAXCMPTS];
	private static double[] TMP3 = new double[$MAXRZ];
	private static double[] CHINDEX = new double[$MAXRZ];
	private static int TMP1 = 0;
	public static String[] EXPLANATIONS = new String[$MAXCMPTS];
	private static int IRL = 0;

	private static boolean SFIG = false;
	private static int COUNT = 0;

	// for different formats:
	//private static boolean zoneB = false;

	// "***************************************************************************
	// "
	// " VERSION 1 MATERIALGRID.MORTRAN VERSION 1
	// " \ ================== /
	// " \ /
	// " By Aaron Merovitz and D.W.O. Rogers, Feb 1998
	// "
	/**
	 * This routine provides a material grid output.
	 * @param NRADIAL NRADIAL : number of radial results to display
	 * @param NDEPTH NDEPTH : number of depth results to display
	 * @param MASSVOL MASSVOL : array of the masses, or volumes for each region
	 * @param MORV MORV: 1 for mass, 2 for volume.
	 * @param ECUTS the 'ECUT' array
	 * @param PCUTS the 'PCUT' arrays
	 * @param RCYL RCYL : array of outer boundaries of each radial bin (1st radius=0)
	 * @param ZPLANE ZPLANE : array of depth boundaries (1st plane not always=0)
	 * @param MED MED : array of the medium numbers for each IRL.
	 * @param MEDIA MEDIA : array of medium names.
	 */
	public static void MATERIALGRID(int NRADIAL, int NDEPTH,
			double[][] MASSVOL, int MORV, double[] ECUTS, double[] PCUTS,
			double[] RCYL, double[] ZPLANE, int[] MED, String[] MEDIA)
	// ,String[] CDSTBL,String[] CTRTBL,String[] CABSRB)
	{
		// "
		// " DESCRIPTION : routine provides a material grid output.
		// " OF : the grid will be rotated if nradial>3
		// " THE : and ndepth<=3. the grid contains the IRL, the type
		// " SUBROUTINE : of the region and material aswell as it's mass.
		// "
		// " \ one must add the macro '$maxcmpts'(usually=$maxit)
		// " FOR INSTALLING \ to the program in which this subroutine is called.
		// " THIS > be careful that all arrays passed to the subroutine
		// " SUBROUTINE / (especially 'MASSVOL()' and 'MEDIA()') have
		// " / exactly the same dimensions.

		// " DESCRIPTION OF THE ARGUMENTS:
		// " NRADIAL : number of radial results to display
		// " NDEPTH : number of depth results to display
		// " MASSVOL : array of the masses, or volumes for each region
		// " MORV : 1 for mass, 2 for volume.
		// " ECUTS, PCUTS: the 'ECUT' and 'PCUT' arrays.
		// " RCYL : array of outer boundaries of each radial bin(1st
		// " radius=0)
		// " ZPLANE : array of depth boundaries (1st plane not always=0)
		// " MED : array of the medium numbers for each IRL.
		// " MEDIA : 2-d array of the masse for each medium number.
		// " CDSTBL : array that contains a 'D' for each dose scoring
		// " region(for each IRL). default is 0.
		// " CTRTBL : contains a 'T' for each tracking region. default 0
		// " CABSRB : contains an 'A' for each tatally absorbing region.
		// " 0 default means these values are not calculated.
		// "***************************************************************************
		// ;IMPLICIT NONE;

		String s = "";
		//int ll = 0;
		//zoneB = false;

		ROT = false; // "For a rotated grid (.FALSE. = off)"
		//ESTEPSON = false;
		ECUTON = false;
		PCUTON = false; // "default these options off"

		NCOMP = 1; // "start the grid with one component(minimum)"
		PGTHROW = 14 + NCOMP; // "For printing: throw page at 66 lines(15 1st for title)"

		// IOUT=1; //"Designates output to file fort.1, the *.egs4lst file"
		// "Set up the bin number indicators"
		DEEPNUM = NDEPTH;
		RADNUM = NRADIAL;
		// "Set up the bin indicators"
		if (RCYL[1] == 0.)// 0 biased
		{
			for (int IX = 1; IX <= RADNUM + 1; IX++) {
				RADIAL_BINS[IX - 1] = RCYL[IX];
			}
		} else {
			for (int IX = 1; IX <= RADNUM + 1; IX++) {
				RADIAL_BINS[IX] = RCYL[IX];
			}
		}
		for (int IZ = 1; IZ <= DEEPNUM + 1; IZ++) {
			DEPTH_BINS[IZ - 1] = ZPLANE[IZ - 1];
		}
		// "Load the array of VALUES to be output"
		if (MORV == 1) {
			LABELS[NCOMP - 1] = "MASS";
			EXPLANATIONS[NCOMP - 1] = "MASS OF EACH REGION IN GRAMS";
		}
		if (MORV == 2) {
			LABELS[NCOMP - 1] = "VOL ";
			EXPLANATIONS[NCOMP - 1] = "VOLUME OF EACH REGION IN cm^3";
		}
		for (int IX = 1; IX <= RADNUM; IX++) {
			for (int IZ = 1; IZ <= DEEPNUM; IZ++) {
				VALUES[IZ - 1][IX - 1][NCOMP - 1] = MASSVOL[IZ - 1][IX - 1];
				IRL = IZ + DEEPNUM * (IX - 1) + 1;// int
													// lll=CDSTBL.length;System.out.println(lll+"  "+IRL);
				if (CDSTBL[0].compareTo("0") == 0) {
					CDSTBL[IRL - 1] = " ";
				}// could be pass by param (array!!)
				if (CTRTBL[0].compareTo("0") == 0) {
					CTRTBL[IRL - 1] = " ";
				}
				if (CABSRB[0].compareTo("0") == 0) {
					CABSRB[IRL - 1] = " ";
				}
				if (ECUTS[IRL - 1] != ECUTS[1]) {
					ECUTON = true;
				}
				if (PCUTS[IRL - 1] != PCUTS[1]) {
					PCUTON = true;
				}
			}
		}
		if (ECUTON) {
			NCOMP = NCOMP + 1;
			LABELS[NCOMP - 1] = "ECUT";
			EXPLANATIONS[NCOMP - 1] = "ECUT (PRINTED BECAUSE DIFFERENT FROM GLOBAL)";
			for (int IX = 1; IX <= RADNUM; IX++)// DO IX=1, RADNUM [
			{
				for (int IZ = 1; IZ <= DEEPNUM; IZ++)// DO IZ=1, DEEPNUM [
				{
					IRL = IZ + DEEPNUM * (IX - 1) + 1;
					VALUES[IZ - 1][IX - 1][NCOMP - 1] = ECUTS[IRL - 1];
					// VALUES(IZ, IX, NCOMP)=ECUTS(IRL);
				}
			}
		}
		if (PCUTON) {
			NCOMP = NCOMP + 1;
			LABELS[NCOMP - 1] = "PCUT";
			EXPLANATIONS[NCOMP - 1] = "PCUT (PRINTED BECAUSE DIFFERENT FROM GLOBAL)";
			for (int IX = 1; IX <= RADNUM; IX++)// DO IX=1, RADNUM [
			{
				for (int IZ = 1; IZ <= DEEPNUM; IZ++)// DO IZ=1, DEEPNUM [
				{
					IRL = IZ + DEEPNUM * (IX - 1) + 1;
					VALUES[IZ - 1][IX - 1][NCOMP - 1] = PCUTS[IRL - 1];
					// VALUES(IZ, IX, NCOMP)=PCUTS(IRL);
				}
			}
		}
		PGTHROW = 14 + NCOMP; // "For printing: throw page at 66 lines(for title)"
		// "Check for rotation of the grid"
		if ((DEEPNUM <= 3) && (RADNUM > 3)) {
			ROT = true;
			CHINDEX[RADNUM] = RADIAL_BINS[RADNUM];// CHINDEX(RADNUM+1)=RADIAL_BINS(RADNUM+1);
			for (int IX = 1; IX <= RADNUM; IX++)// DO IX=1, RADNUM [
			{
				TMP3[IX - 1] = DEPTH_BINS[IX - 1];
				CHINDEX[IX - 1] = RADIAL_BINS[IX - 1];
				DEPTH_BINS[IX - 1] = CHINDEX[IX - 1];
				RADIAL_BINS[IX - 1] = TMP3[IX - 1];

				for (int IZ = 1; IZ <= DEEPNUM; IZ++)// DO IZ=1, DEEPNUM [
				{
					for (int ICOMP = 1; ICOMP <= NCOMP; ICOMP++) {
						TMP2[IZ - 1][IX - 1][ICOMP - 1] = VALUES[IZ - 1][IX - 1][ICOMP - 1];
						// TMP2(IZ, IX, ICOMP)=VALUES(IZ, IX, ICOMP);
					} // "end do ICOMP"
				} // " end DO IZ"
			} // "end DO IX"

			for (int IX = 1; IX <= RADNUM; IX++)// DO IX=1, RADNUM [
			{
				for (int IZ = 1; IZ <= DEEPNUM; IZ++)// DO IZ=1, DEEPNUM [
				{
					for (int ICOMP = 1; ICOMP <= NCOMP; ICOMP++)// DO ICOMP=1,
																// NCOMP [
					{
						VALUES[IX - 1][IZ - 1][ICOMP - 1] = TMP2[IZ - 1][IX - 1][ICOMP - 1];
						// VALUES(IX, IZ, ICOMP)=TMP2(IZ, IX, ICOMP);
					}
				}
			}

			DEPTH_BINS[RADNUM] = CHINDEX[RADNUM];// DEPTH_BINS(RADNUM+1)=CHINDEX(RADNUM+1);
			TMP1 = RADNUM;
			RADNUM = DEEPNUM;
			DEEPNUM = TMP1;
		} // "end check for rotation"
		// WRITE (IOUT, *) '\f'; "Page break"
		// WRITE(IOUT, 400) ' '; call egs_fdate(iout); write(iout,'(//)');
		if (ROT) {
			s = EGS4.format("", 19);
			EGS4.seqStr = s + "ZONAL MATERIAL GRID: ROTATED";
			if (EGS4.iprint > 1)
				eq.printSequence(EGS4.seqStr);
			s = EGS4.format("", 21);
			EGS4.seqStr = s + "**************************";
			if (EGS4.iprint > 1)
				eq.printSequence(EGS4.seqStr);
			EGS4.seqStr = "\n";
			if (EGS4.iprint > 1)
				eq.printSequence(EGS4.seqStr);
			s = EGS4.format("", 4);
			EGS4.seqStr = s
					+ "/X/Y/Z/MED :  X = ' ' IS DEFAULT: OPTION NOT USED";
			if (EGS4.iprint > 1)
				eq.printSequence(EGS4.seqStr);
			s = EGS4.format("", 18);
			EGS4.seqStr = s + "X = 'D' IF DOSE SCORING REGION";
			if (EGS4.iprint > 1)
				eq.printSequence(EGS4.seqStr);
			s = EGS4.format("", 18);
			EGS4.seqStr = s + "X = 'C' IF CAVITY REGION";
			if (EGS4.iprint > 1)
				eq.printSequence(EGS4.seqStr);
			s = EGS4.format("", 18);
			EGS4.seqStr = s + "X = 'S' IF SPR SCORING REGION";
			if (EGS4.iprint > 1)
				eq.printSequence(EGS4.seqStr);
			s = EGS4.format("", 18);
			EGS4.seqStr = s + "Y = 'T' IF TRACKING REGION";
			if (EGS4.iprint > 1)
				eq.printSequence(EGS4.seqStr);
			s = EGS4.format("", 18);
			EGS4.seqStr = s + "Z = 'A' IF TOTALLY ABSORBING REGION";
			if (EGS4.iprint > 1)
				eq.printSequence(EGS4.seqStr);
			s = EGS4.format("", 16);
			EGS4.seqStr = s + "MED = MEDIUM NAME, 11 CHARACTER ABREVIATION";
			if (EGS4.iprint > 1)
				eq.printSequence(EGS4.seqStr);

			// WRITE(IOUT, 93) TITLE;
		} else {
			s = EGS4.format("", 19);
			EGS4.seqStr = s + "ZONAL MATERIAL GRID: NON-ROTATED";
			if (EGS4.iprint > 1)
				eq.printSequence(EGS4.seqStr);
			s = EGS4.format("", 21);
			EGS4.seqStr = s + "******************************";
			if (EGS4.iprint > 1)
				eq.printSequence(EGS4.seqStr);
			EGS4.seqStr = "\n";
			if (EGS4.iprint > 1)
				eq.printSequence(EGS4.seqStr);
			s = EGS4.format("", 4);
			EGS4.seqStr = s
					+ "/X/Y/Z/MED :  X = ' ' IS DEFAULT: OPTION NOT USED";
			if (EGS4.iprint > 1)
				eq.printSequence(EGS4.seqStr);
			s = EGS4.format("", 18);
			EGS4.seqStr = s + "  = 'D' IF DOSE SCORING REGION";
			if (EGS4.iprint > 1)
				eq.printSequence(EGS4.seqStr);
			s = EGS4.format("", 18);
			EGS4.seqStr = s + "  = 'C' IF CAVITY REGION";
			if (EGS4.iprint > 1)
				eq.printSequence(EGS4.seqStr);
			s = EGS4.format("", 18);
			EGS4.seqStr = s + "  = 'S' IF SPR SCORING REGION";
			if (EGS4.iprint > 1)
				eq.printSequence(EGS4.seqStr);
			s = EGS4.format("", 18);
			EGS4.seqStr = s + "Y = 'T' IF TRACKING REGION";
			if (EGS4.iprint > 1)
				eq.printSequence(EGS4.seqStr);
			s = EGS4.format("", 18);
			EGS4.seqStr = s + "Z = 'A' IF TOTALLY ABSORBING REGION";
			if (EGS4.iprint > 1)
				eq.printSequence(EGS4.seqStr);
			s = EGS4.format("", 16);
			EGS4.seqStr = s + "MED = MEDIUM NAME, 11 CHARACTER ABREVIATION";
			if (EGS4.iprint > 1)
				eq.printSequence(EGS4.seqStr);

			// WRITE(IOUT, 94) TITLE;
		}

		for (int ICOMP = 1; ICOMP <= NCOMP; ICOMP++) {
			s = EGS4.format("", 10) + EGS4.format(LABELS[ICOMP - 1], 4) + " = "
					+ EGS4.format(EXPLANATIONS[ICOMP - 1], 60);
			EGS4.seqStr = s;
			if (EGS4.iprint > 1)
				eq.printSequence(EGS4.seqStr);
		}
		// WRITE(IOUT, 95) (LABELS(ICOMP), EXPLANATIONS(ICOMP), ICOMP=1, NCOMP);
		// //95 FORMAT (T10, A4, ' = ', A60);
		// "Make the grids"
		// DO HSET=1,RADNUM,3 [$MKGRID;] "Feed in horiz. sets of three"
		for (HSET = 1; HSET <= RADNUM; HSET = HSET + 3) {
			MKGRID(MED, MEDIA, CDSTBL, CTRTBL, CABSRB);
			// $MKGRID;
		} // "Feed in horiz. sets of three"
		// "FORMATS"
		// 1 FORMAT (T11, '|', 2X, A4, 2X, 1PE10.3, 2X, ' |');
		// 2 FORMAT (T11, '|', 2 (2X, A4, 2X, 1PE10.3, 2X, ' |'));
		// 3 FORMAT (T11, '|', 3 (2X, A4, 2X, 1PE10.3, 2X, ' |'));
		// 4 FORMAT (T11, '|', 1X, '/', 3 (A1, '/'), 11A1, ' |');
		// 5 FORMAT (T11, '|', 2 (1X, '/', 3 (A1, '/'), 11A1, ' |'));
		// 6 FORMAT (T11, '|', 3 (1X, '/', 3 (A1, '/'), 11A1, ' |'));
		// 10 FORMAT (1X, F9.4, T11, 23 ('-'));
		// 11 FORMAT (1X, F9.4, T11, 45 ('-'));
		// 12 FORMAT (1X, F9.4, T11, 67 ('-'));
		// 13 FORMAT (T11, '|',
		// 3 ('IRL', I3, 1X, 'IZ ', I3, 1X, 'IX ', I3, 1X, '|'));
		// 14 FORMAT (T11, '|',
		// 2 ('IRL', I3, 1X, 'IZ ', I3, 1X, 'IX ', I3, 1X, '|'));
		// 15 FORMAT (T11, '|',
		// 'IRL', I3, 1X, 'IZ ', I3, 1X, 'IX ', I3, 1X, '|');
		// 91 FORMAT (/ T7, F9.4, T30, F9.4, T53, F9.4, T70, F9.4);
		// 93 FORMAT (' ',79A1 //
		// T19, 29H ZONAL MATERIAL GRID: ROTATED /
		// T21, '**************************'/,
		// /T4 , '/X/Y/Z/MED : X = " " IS DEFAULT: OPTION NOT USED'/,
		// T18, ' = "D" IF DOSE SCORING REGION'/,
		// T18, 'X = "C" IF CAVITY REGION'/,
		// T18, 'X = "S" IF SPR SCORING REGION'/,
		// T18, 'Y = "T" IF TRACKING REGION'/,
		// T18, 'Z = "A" IF TOTALLY ABSORBING REGION'/,
		// T16, 'MED = MEDIUM NAME, 11 CHARACTER ABREVIATION');
		// 94 FORMAT (' ',79A1 //
		// T19, 33H ZONAL MATERIAL GRID: NON-ROTATED /
		// T21, '******************************'/
		// /T4 , '/X/Y/Z/MED : X = " " IS DEFAULT: OPTION NOT USED'/,
		// T18, ' = "D" IF DOSE SCORING REGION'/,
		// T18, ' = "C" IF CAVITY REGION'/,
		// T18, ' = "S" IF SPR SCORING REGION'/,
		// T18, 'Y = "T" IF TRACKING REGION'/,
		// T18, 'Z = "A" IF TOTALLY ABSORBING REGION'/,
		// T16, 'MED = MEDIUM NAME, 11 CHARACTER ABREVIATION');
		// 95 FORMAT (T10, A4, ' = ', A60);
		// 400 FORMAT (/T54,a1,$);

		return;// RETURN;
	}// END;
		// "    end of grids.mortran (Rev 1.2 last edited  2002-05-06 10:09:54-04)"

	// "This MACRO outputs a grid of values."
	// "There is room for three sets of values(blocks or boxes) in"
	// "80 columns(standard window size), so the grids each"
	// "contain sets of three values, e.g. radii(if non-rotated)."
	/**
	 * Internally used. It outputs a grid of values. Called by MATERIALGRID.
	 * @param MED MED
	 * @param MEDIA MEDIA
	 * @param CDSTBL CDSTBL
	 * @param CTRTBL CTRTBL
	 * @param CABSRB CABSRB
	 */
	private static void MKGRID(int[] MED, String[] MEDIA, String[] CDSTBL,
			String[] CTRTBL, String[] CABSRB) {
		String s = "";
		int ll = 0;
		String hs = "";

		// REPLACE {$MKGRID;} WITH {
		// ;//------------->l=jxx % 169;//l = mod(jxx, 169) ;
		// "Set up the formatting indicator FMT"
		if (RADNUM - HSET > 1) {
			FMT = 3;
		} else {
			// if (MOD(RADNUM,3) = 1) [FMT=1;]
			// if (MOD(RADNUM,3) = 2) [FMT=2;]

			if (RADNUM % 3 == 1) {
				FMT = 1;
			}
			if (RADNUM % 3 == 2) {
				FMT = 2;
			}

		}
		// "Horizontal axis indicators"
		// WRITE (IOUT, 91) (RADIAL_BINS(IX), IX=HSET,HSET+FMT);
		// 91 FORMAT (/ T7, F9.4, T30, F9.4, T53, F9.4, T70, F9.4);
		for (int IX = HSET; IX <= HSET + FMT; IX++) {
			if (IX == HSET) {
				s = EGS4.format("", 7);
			} else if (IX == HSET + 1) {
				ll = 30 - ll;
				if (ll > 0) {
					s = EGS4.format("", ll);
				} else {
					s = EGS4.format("", 30);
				}
			} else if (IX == HSET + 2) {
				ll = 53 - ll;
				if (ll > 0) {
					s = EGS4.format("", ll);
				} else {
					s = EGS4.format("", 53);
				}
			} else if (IX == HSET + 3) {
				ll = 70 - ll;
				if (ll > 0) {
					s = EGS4.format("", ll);
				} else {
					s = EGS4.format("", 70);
				}
			}

			// ll=s.length();ll=60-ll;s=s+EGS4.format("",ll);
			// EGS4.seqStr=s+EGS4.format(RADIAL_BINS[IX-1],9,true);
			// if(EGS4.iprint>1)
			// eq.printSequence(EGS4.seqStr);
			if (IX == HSET) {
				EGS4.seqStr = s + EGS4.format(RADIAL_BINS[IX - 1], 9, true);
				ll = EGS4.seqStr.length();
			} else {
				EGS4.seqStr = EGS4.seqStr + s
						+ EGS4.format(RADIAL_BINS[IX - 1], 9, true);
				ll = EGS4.seqStr.length();
			}
		}
		// EGS4.seqStr=s+EGS4.format(RADIAL_BINS[IX-1],9,true);
		if (EGS4.iprint > 1)
			eq.printSequence(EGS4.seqStr);// auto salt at new line

		PGTHROW = PGTHROW + 1;
		for (int IZ = 1; IZ <= DEEPNUM; IZ++) {
			// "This is for IRL, IZ and IX"
			if (ROT) {
				REGNUM = 2 + (IZ - 1) * FMT + (HSET - 1) * RADNUM;
			} else {
				REGNUM = (IZ + 1 + ((HSET - 1) * DEEPNUM));
			}
			// "Setting up medium name in each region"
			if (ROT) {
				MNUM1 = MED[REGNUM - 1];
				MNUM2 = MED[REGNUM];
				MNUM3 = MED[REGNUM + 1];
			} else {
				MNUM1 = MED[REGNUM - 1];
				MNUM2 = MED[REGNUM + DEEPNUM - 1];
				MNUM3 = MED[REGNUM + DEEPNUM * 2 - 1];
			}
			if (MNUM1 == 0) {
				MED_NAME1[0] = "V";
				MED_NAME1[1] = "A";
				MED_NAME1[2] = "C";
				MED_NAME1[3] = "U";
				MED_NAME1[4] = "U";
				MED_NAME1[5] = "M";
				for (int J = 7; J <= 11; J++) {
					MED_NAME1[J - 1] = " ";
				}
			} else {
				for (int COUNT = 1; COUNT <= 11; COUNT++) {
					String meds = "";
					if (MEDIA[MNUM1 - 1].length() >= COUNT) {
						meds = MEDIA[MNUM1 - 1].substring(COUNT - 1, COUNT);
					} else {
						meds = " ";
					}

					MED_NAME1[COUNT - 1] = meds;// MEDIA(COUNT, MNUM1);
				}
			}
			if (MNUM2 == 0) {
				MED_NAME2[0] = "V";
				MED_NAME2[1] = "A";
				MED_NAME2[2] = "C";
				MED_NAME2[3] = "U";
				MED_NAME2[4] = "U";
				MED_NAME2[5] = "M";
				for (int J = 7; J <= 11; J++) {
					MED_NAME2[J - 1] = " ";
				}
			} else {
				// DO COUNT=1, 11 [MED_NAME2(COUNT) = MEDIA(COUNT,MNUM2);]
				for (int COUNT = 1; COUNT <= 11; COUNT++) {
					String meds = "";
					if (MEDIA[MNUM2 - 1].length() >= COUNT) {
						meds = MEDIA[MNUM2 - 1].substring(COUNT - 1, COUNT);
					} else {
						meds = " ";
					}

					MED_NAME2[COUNT - 1] = meds;
				}
			}
			if (MNUM3 == 0) {
				MED_NAME3[0] = "V";
				MED_NAME3[1] = "A";
				MED_NAME3[2] = "C";
				MED_NAME3[3] = "U";
				MED_NAME3[4] = "U";
				MED_NAME3[5] = "M";
				for (int J = 7; J <= 11; J++) {
					MED_NAME3[J - 1] = " ";
				}
			} else {
				// DO COUNT=1, 11 [MED_NAME3(COUNT) = MEDIA(COUNT,MNUM3);]
				for (int COUNT = 1; COUNT <= 11; COUNT++) {
					String meds = "";
					if (MEDIA[MNUM3 - 1].length() >= COUNT) {
						meds = MEDIA[MNUM3 - 1].substring(COUNT - 1, COUNT);
					} else {
						meds = " ";
					}

					MED_NAME3[COUNT - 1] = meds;
				}
			}
			if (FMT == 1) {
				// WRITE (IOUT, 10) DEPTH_BINS(IZ);
				hs = "";
				for (int i = 1; i <= 23; i++) {
					hs = hs + "-";
				}
				s = EGS4.format("", 1)
						+ EGS4.format(DEPTH_BINS[IZ - 1], 9, true);
				ll = s.length();
				ll = 11 - ll;
				if (ll > 0) {
					s = s + EGS4.format("", ll);
				}
				s = s + hs;// EGS4.format("-",23);
				EGS4.seqStr = s;
				if (EGS4.iprint > 1)
					eq.printSequence(EGS4.seqStr);

				if (ROT) {
					s = EGS4.format("", 11) + "|" + "IRL"
							+ EGS4.format(REGNUM, 3) + EGS4.format("", 1)
							+ "IZ " + EGS4.format(HSET, 3) + EGS4.format("", 1)
							+ "IX " + EGS4.format(IZ, 3) + EGS4.format("", 1)
							+ "|";
					EGS4.seqStr = s;
					if (EGS4.iprint > 1)
						eq.printSequence(EGS4.seqStr);

					// WRITE (IOUT, 15) REGNUM,HSET,IZ;
					// 15 FORMAT (T11, '|',
					// 'IRL', I3, 1X, 'IZ ', I3, 1X, 'IX ', I3, 1X, '|');
				} else {
					s = EGS4.format("", 11) + "|" + "IRL"
							+ EGS4.format(REGNUM, 3) + EGS4.format("", 1)
							+ "IZ " + EGS4.format(IZ, 3) + EGS4.format("", 1)
							+ "IX " + EGS4.format(HSET, 3) + EGS4.format("", 1)
							+ "|";
					EGS4.seqStr = s;
					if (EGS4.iprint > 1)
						eq.printSequence(EGS4.seqStr);

					// WRITE (IOUT, 15) REGNUM,IZ,HSET;
				}
				String mds = "";
				for (int J = 1; J <= 11; J++) {
					mds = mds + MED_NAME1[J - 1];
				}
				s = EGS4.format("", 11) + "|" + EGS4.format("", 1) + "/"
						+ EGS4.format(CDSTBL[REGNUM - 1], 1) + "/"
						+ EGS4.format(CTRTBL[REGNUM - 1], 1) + "/"
						+ EGS4.format(CABSRB[REGNUM - 1], 1) + "/"
						+ EGS4.format(mds, 11) + "  |";
				EGS4.seqStr = s;
				if (EGS4.iprint > 1)
					eq.printSequence(EGS4.seqStr);

				// WRITE (IOUT, 4) CDSTBL(REGNUM),CTRTBL(REGNUM),
				// CABSRB(REGNUM), (MED_NAME1(J),J=1,11);
				// 4 FORMAT (T11, '|', 1X, '/', 3 (A1, '/'), 11A1, ' |');
			}
			if (FMT == 2) {
				// WRITE (IOUT, 11) DEPTH_BINS(IZ);
				hs = "";
				for (int i = 1; i <= 45; i++) {
					hs = hs + "-";
				}
				s = EGS4.format("", 1)
						+ EGS4.format(DEPTH_BINS[IZ - 1], 9, true);
				ll = s.length();
				ll = 11 - ll;
				if (ll > 0) {
					s = s + EGS4.format("", ll);
				}
				s = s + hs;// EGS4.format("-",45);
				EGS4.seqStr = s;
				if (EGS4.iprint > 1)
					eq.printSequence(EGS4.seqStr);

				if (ROT) {
					s = EGS4.format("", 11) + "|" + "IRL"
							+ EGS4.format(REGNUM, 3) + EGS4.format("", 1)
							+ "IZ " + EGS4.format(HSET, 3) + EGS4.format("", 1)
							+ "IX " + EGS4.format(IZ, 3) + EGS4.format("", 1)
							+ "|" + "IRL" + EGS4.format(REGNUM + 1, 3)
							+ EGS4.format("", 1) + "IZ "
							+ EGS4.format(HSET + 1, 3) + EGS4.format("", 1)
							+ "IX " + EGS4.format(IZ, 3) + EGS4.format("", 1)
							+ "|";

					EGS4.seqStr = s;
					if (EGS4.iprint > 1)
						eq.printSequence(EGS4.seqStr);

					// WRITE(IOUT,14) REGNUM,HSET,IZ,
					// REGNUM+1,HSET+1,IZ;
					// 14 FORMAT (T11, '|',
					// 2 ('IRL', I3, 1X, 'IZ ', I3, 1X, 'IX ', I3, 1X, '|'));

					String mds2 = "";
					String mds1 = "";
					for (int J = 1; J <= 11; J++) {
						mds1 = mds1 + MED_NAME1[J - 1];
						mds2 = mds2 + MED_NAME2[J - 1];
					}
					s = EGS4.format("", 11) + "|" + EGS4.format("", 1) + "/"
							+ EGS4.format(CDSTBL[REGNUM - 1], 1) + "/"
							+ EGS4.format(CTRTBL[REGNUM - 1], 1) + "/"
							+ EGS4.format(CABSRB[REGNUM - 1], 1) + "/"
							+ EGS4.format(mds1, 11) + "  |"
							+ EGS4.format("", 1) + "/"
							+ EGS4.format(CDSTBL[REGNUM], 1) + "/"
							+ EGS4.format(CTRTBL[REGNUM], 1) + "/"
							+ EGS4.format(CABSRB[REGNUM], 1) + "/"
							+ EGS4.format(mds2, 11) + "  |";
					EGS4.seqStr = s;
					if (EGS4.iprint > 1)
						eq.printSequence(EGS4.seqStr);

					// WRITE(IOUT, 5) CDSTBL(REGNUM),CTRTBL(REGNUM),
					// CABSRB(REGNUM), (MED_NAME1(J),J=1,11),
					// CDSTBL(REGNUM+1),CTRTBL(REGNUM+1),
					// CABSRB(REGNUM+1), (MED_NAME2(J),J=1,11);
					// 5 FORMAT (T11, '|', 2 (1X, '/', 3 (A1, '/'), 11A1, '
					// |'));
				} else {
					s = EGS4.format("", 11) + "|" + "IRL"
							+ EGS4.format(REGNUM, 3) + EGS4.format("", 1)
							+ "IZ " + EGS4.format(IZ, 3) + EGS4.format("", 1)
							+ "IX " + EGS4.format(HSET, 3) + EGS4.format("", 1)
							+ "|" + "IRL" + EGS4.format(REGNUM + DEEPNUM, 3)
							+ EGS4.format("", 1) + "IZ " + EGS4.format(IZ, 3)
							+ EGS4.format("", 1) + "IX "
							+ EGS4.format(HSET + 1, 3) + EGS4.format("", 1)
							+ "|";

					EGS4.seqStr = s;
					if (EGS4.iprint > 1)
						eq.printSequence(EGS4.seqStr);

					// WRITE(IOUT,14) REGNUM,IZ,HSET,
					// REGNUM+DEEPNUM,IZ,HSET+1;

					String mds2 = "";
					String mds1 = "";
					for (int J = 1; J <= 11; J++) {
						mds1 = mds1 + MED_NAME1[J - 1];
						mds2 = mds2 + MED_NAME2[J - 1];
					}
					s = EGS4.format("", 11) + "|" + EGS4.format("", 1) + "/"
							+ EGS4.format(CDSTBL[REGNUM - 1], 1) + "/"
							+ EGS4.format(CTRTBL[REGNUM - 1], 1) + "/"
							+ EGS4.format(CABSRB[REGNUM - 1], 1) + "/"
							+ EGS4.format(mds1, 11) + "  |"
							+ EGS4.format("", 1) + "/"
							+ EGS4.format(CDSTBL[REGNUM + DEEPNUM - 1], 1)
							+ "/"
							+ EGS4.format(CTRTBL[REGNUM + DEEPNUM - 1], 1)
							+ "/"
							+ EGS4.format(CABSRB[REGNUM + DEEPNUM - 1], 1)
							+ "/" + EGS4.format(mds2, 11) + "  |";
					EGS4.seqStr = s;
					if (EGS4.iprint > 1)
						eq.printSequence(EGS4.seqStr);

					// WRITE(IOUT, 5) CDSTBL(REGNUM),CTRTBL(REGNUM),
					// CABSRB(REGNUM), (MED_NAME1(J),J=1,11),
					// CDSTBL(REGNUM+DEEPNUM),CTRTBL(REGNUM+DEEPNUM),
					// CABSRB(REGNUM+DEEPNUM),(MED_NAME2(J),J=1,11);
				}
			}
			if (FMT == 3) {
				// WRITE (IOUT, 12) DEPTH_BINS(IZ);
				hs = "";
				for (int i = 1; i <= 67; i++) {
					hs = hs + "-";
				}
				s = EGS4.format("", 1)
						+ EGS4.format(DEPTH_BINS[IZ - 1], 9, true);
				ll = s.length();
				ll = 11 - ll;
				if (ll > 0) {
					s = s + EGS4.format("", ll);
				}
				s = s + hs;// EGS4.format("-",67);
				EGS4.seqStr = s;
				if (EGS4.iprint > 1)
					eq.printSequence(EGS4.seqStr);

				if (ROT) {
					s = EGS4.format("", 11) + "|" + "IRL"
							+ EGS4.format(REGNUM, 3) + EGS4.format("", 1)
							+ "IZ " + EGS4.format(HSET, 3) + EGS4.format("", 1)
							+ "IX " + EGS4.format(IZ, 3) + EGS4.format("", 1)
							+ "|" + "IRL" + EGS4.format(REGNUM + 1, 3)
							+ EGS4.format("", 1) + "IZ "
							+ EGS4.format(HSET + 1, 3) + EGS4.format("", 1)
							+ "IX " + EGS4.format(IZ, 3) + EGS4.format("", 1)
							+ "|" + "IRL" + EGS4.format(REGNUM + 2, 3)
							+ EGS4.format("", 1) + "IZ "
							+ EGS4.format(HSET + 2, 3) + EGS4.format("", 1)
							+ "IX " + EGS4.format(IZ, 3) + EGS4.format("", 1)
							+ "|";

					EGS4.seqStr = s;
					if (EGS4.iprint > 1)
						eq.printSequence(EGS4.seqStr);

					// WRITE (IOUT, 13)
					// REGNUM,HSET,IZ,REGNUM+1,
					// HSET+1,IZ,REGNUM+2,HSET+2,IZ;
					// 13 FORMAT (T11, '|',
					// 3 ('IRL', I3, 1X, 'IZ ', I3, 1X, 'IX ', I3, 1X, '|'));

					String mds2 = "";
					String mds1 = "";
					String mds3 = "";
					for (int J = 1; J <= 11; J++) {
						mds1 = mds1 + MED_NAME1[J - 1];
						mds2 = mds2 + MED_NAME2[J - 1];
						mds3 = mds3 + MED_NAME3[J - 1];
					}
					s = EGS4.format("", 11) + "|" + EGS4.format("", 1) + "/"
							+ EGS4.format(CDSTBL[REGNUM - 1], 1) + "/"
							+ EGS4.format(CTRTBL[REGNUM - 1], 1) + "/"
							+ EGS4.format(CABSRB[REGNUM - 1], 1) + "/"
							+ EGS4.format(mds1, 11) + "  |"
							+ EGS4.format("", 1) + "/"
							+ EGS4.format(CDSTBL[REGNUM], 1) + "/"
							+ EGS4.format(CTRTBL[REGNUM], 1) + "/"
							+ EGS4.format(CABSRB[REGNUM], 1) + "/"
							+ EGS4.format(mds2, 11) + "  |"
							+ EGS4.format("", 1) + "/"
							+ EGS4.format(CDSTBL[REGNUM + 1], 1) + "/"
							+ EGS4.format(CTRTBL[REGNUM + 1], 1) + "/"
							+ EGS4.format(CABSRB[REGNUM + 1], 1) + "/"
							+ EGS4.format(mds3, 11) + "  |";
					EGS4.seqStr = s;
					if (EGS4.iprint > 1)
						eq.printSequence(EGS4.seqStr);

					// WRITE (IOUT, 6) CDSTBL(REGNUM),CTRTBL(REGNUM),
					// CABSRB(REGNUM), (MED_NAME1(J),J=1,11),
					// CDSTBL(REGNUM+1),CTRTBL(REGNUM+1),
					// CABSRB(REGNUM+1), (MED_NAME2(J),J=1,11),
					// CDSTBL(REGNUM+2),CTRTBL(REGNUM+2),
					// CABSRB(REGNUM+2), (MED_NAME3(J),J=1,11);
					// 6 FORMAT (T11, '|', 3 (1X, '/', 3 (A1, '/'), 11A1, '
					// |'));
				} else {
					s = EGS4.format("", 11) + "|" + "IRL"
							+ EGS4.format(REGNUM, 3) + EGS4.format("", 1)
							+ "IZ " + EGS4.format(IZ, 3) + EGS4.format("", 1)
							+ "IX " + EGS4.format(HSET, 3) + EGS4.format("", 1)
							+ "|" + "IRL" + EGS4.format(REGNUM + DEEPNUM, 3)
							+ EGS4.format("", 1) + "IZ " + EGS4.format(IZ, 3)
							+ EGS4.format("", 1) + "IX "
							+ EGS4.format(HSET + 1, 3) + EGS4.format("", 1)
							+ "|" + "IRL"
							+ EGS4.format(REGNUM + DEEPNUM * 2, 3)
							+ EGS4.format("", 1) + "IZ " + EGS4.format(IZ, 3)
							+ EGS4.format("", 1) + "IX "
							+ EGS4.format(HSET + 2, 3) + EGS4.format("", 1)
							+ "|";

					EGS4.seqStr = s;
					if (EGS4.iprint > 1)
						eq.printSequence(EGS4.seqStr);

					// WRITE (IOUT, 13)
					// REGNUM,IZ,HSET,REGNUM+DEEPNUM,
					// IZ,HSET+1,REGNUM+DEEPNUM*2,IZ,HSET+2;

					String mds2 = "";
					String mds1 = "";
					String mds3 = "";
					for (int J = 1; J <= 11; J++) {
						mds1 = mds1 + MED_NAME1[J - 1];
						mds2 = mds2 + MED_NAME2[J - 1];
						mds3 = mds3 + MED_NAME3[J - 1];
					}
					s = EGS4.format("", 11) + "|" + EGS4.format("", 1) + "/"
							+ EGS4.format(CDSTBL[REGNUM - 1], 1) + "/"
							+ EGS4.format(CTRTBL[REGNUM - 1], 1) + "/"
							+ EGS4.format(CABSRB[REGNUM - 1], 1) + "/"
							+ EGS4.format(mds1, 11) + "  |"
							+ EGS4.format("", 1) + "/"
							+ EGS4.format(CDSTBL[REGNUM + DEEPNUM - 1], 1)
							+ "/"
							+ EGS4.format(CTRTBL[REGNUM + DEEPNUM - 1], 1)
							+ "/"
							+ EGS4.format(CABSRB[REGNUM + DEEPNUM - 1], 1)
							+ "/" + EGS4.format(mds2, 11) + "  |"
							+ EGS4.format("", 1) + "/"
							+ EGS4.format(CDSTBL[REGNUM + DEEPNUM * 2 - 1], 1)
							+ "/"
							+ EGS4.format(CTRTBL[REGNUM + DEEPNUM * 2 - 1], 1)
							+ "/"
							+ EGS4.format(CABSRB[REGNUM + DEEPNUM * 2 - 1], 1)
							+ "/" + EGS4.format(mds3, 11) + "  |";
					EGS4.seqStr = s;
					if (EGS4.iprint > 1)
						eq.printSequence(EGS4.seqStr);

					// WRITE (IOUT, 6) CDSTBL(REGNUM),CTRTBL(REGNUM),
					// CABSRB(REGNUM), (MED_NAME1(J),J=1,11),
					// CDSTBL(REGNUM+DEEPNUM),CTRTBL(REGNUM+DEEPNUM),
					// CABSRB(REGNUM+DEEPNUM), (MED_NAME2(J),J=1,11),
					// CDSTBL(REGNUM+DEEPNUM*2),CTRTBL(REGNUM+DEEPNUM*2),
					// CABSRB(REGNUM+DEEPNUM*2), (MED_NAME3(J),J=1,11);
				}
			}

			PGTHROW = PGTHROW + 3;
			// "Main part of the grid"
			for (int ICOMP = 1; ICOMP <= NCOMP; ICOMP++) {
				if (FMT == 3) {
					s = EGS4.format("", 11)
							+ "|"
							+ EGS4.format("", 2)
							+ EGS4.format(LABELS[ICOMP - 1], 4)
							+ EGS4.format("", 2)
							+ EGS4.format(VALUES[IZ - 1][HSET - 1][ICOMP - 1],
									10, false)
							+ EGS4.format("", 2)
							+ " |"
							+ EGS4.format("", 2)
							+ EGS4.format(LABELS[ICOMP - 1], 4)
							+ EGS4.format("", 2)
							+ EGS4.format(VALUES[IZ - 1][HSET][ICOMP - 1], 10,
									false)
							+ EGS4.format("", 2)
							+ " |"
							+ EGS4.format("", 2)
							+ EGS4.format(LABELS[ICOMP - 1], 4)
							+ EGS4.format("", 2)
							+ EGS4.format(VALUES[IZ - 1][HSET + 1][ICOMP - 1],
									10, false) + EGS4.format("", 2) + " |";
					EGS4.seqStr = s;
					if (EGS4.iprint > 1)
						eq.printSequence(EGS4.seqStr);

					// WRITE(IOUT, 3) LABELS(ICOMP),VALUES(IZ,HSET,ICOMP),
					// LABELS(ICOMP),VALUES(IZ,HSET+1,ICOMP),
					// LABELS(ICOMP),VALUES(IZ,HSET+2,ICOMP);
					// 3 FORMAT (T11, '|', 3 (2X, A4, 2X, 1PE10.3, 2X, ' |'));
				}
				if (FMT == 2) {
					s = EGS4.format("", 11)
							+ "|"
							+ EGS4.format("", 2)
							+ EGS4.format(LABELS[ICOMP - 1], 4)
							+ EGS4.format("", 2)
							+ EGS4.format(VALUES[IZ - 1][HSET - 1][ICOMP - 1],
									10, false)
							+ EGS4.format("", 2)
							+ " |"
							+ EGS4.format("", 2)
							+ EGS4.format(LABELS[ICOMP - 1], 4)
							+ EGS4.format("", 2)
							+ EGS4.format(VALUES[IZ - 1][HSET][ICOMP - 1], 10,
									false) + EGS4.format("", 2) + " |";
					EGS4.seqStr = s;
					if (EGS4.iprint > 1)
						eq.printSequence(EGS4.seqStr);

					// WRITE(IOUT, 2) LABELS(ICOMP),VALUES(IZ,HSET,ICOMP),
					// LABELS(ICOMP),VALUES(IZ,HSET+1,ICOMP);
					// 2 FORMAT (T11, '|', 2 (2X, A4, 2X, 1PE10.3, 2X, ' |'));
				}
				if (FMT == 1) {
					s = EGS4.format("", 11)
							+ "|"
							+ EGS4.format("", 2)
							+ EGS4.format(LABELS[ICOMP - 1], 4)
							+ EGS4.format("", 2)
							+ EGS4.format(VALUES[IZ - 1][HSET - 1][ICOMP - 1],
									10, false) + EGS4.format("", 2) + " |";
					EGS4.seqStr = s;
					if (EGS4.iprint > 1)
						eq.printSequence(EGS4.seqStr);

					// WRITE(IOUT, 1) LABELS(ICOMP),VALUES(IZ,HSET,ICOMP);
					// 1 FORMAT (T11, '|', 2X, A4, 2X, 1PE10.3, 2X, ' |');
				}
				PGTHROW = PGTHROW + 1;
			} // "end DO ICOMP"
			// "This piece of code surveys PGTHROW, and throws the page when"
			// "PGTHROW is near 65 lines."
			if ((PGTHROW % 65) > (60 - NCOMP)) {
				// "This is for the last vertical bin"
				if (IZ != DEEPNUM) {// "only for grid cut by PGTHROW"
					if (FMT == 1) {
						// WRITE (IOUT, 10) DEPTH_BINS(IZ+1);
						hs = "";
						for (int i = 1; i <= 23; i++) {
							hs = hs + "-";
						}
						s = EGS4.format("", 1)
								+ EGS4.format(DEPTH_BINS[IZ], 9, true);
						ll = s.length();
						ll = 11 - ll;
						if (ll > 0) {
							s = s + EGS4.format("", ll);
						}
						s = s + hs;// EGS4.format("-",23);
						EGS4.seqStr = s;
						if (EGS4.iprint > 1)
							eq.printSequence(EGS4.seqStr);

					}
					if (FMT == 2) {
						// WRITE (IOUT, 11) DEPTH_BINS(IZ+1);
						hs = "";
						for (int i = 1; i <= 45; i++) {
							hs = hs + "-";
						}
						s = EGS4.format("", 1)
								+ EGS4.format(DEPTH_BINS[IZ], 9, true);
						ll = s.length();
						ll = 11 - ll;
						if (ll > 0) {
							s = s + EGS4.format("", ll);
						}
						s = s + hs;// EGS4.format("-",45);
						EGS4.seqStr = s;
						if (EGS4.iprint > 1)
							eq.printSequence(EGS4.seqStr);

					}
					if (FMT == 3) {
						// WRITE (IOUT, 12) DEPTH_BINS(IZ+1);
						hs = "";
						for (int i = 1; i <= 67; i++) {
							hs = hs + "-";
						}
						s = EGS4.format("", 1)
								+ EGS4.format(DEPTH_BINS[IZ], 9, true);
						ll = s.length();
						ll = 11 - ll;
						if (ll > 0) {
							s = s + EGS4.format("", ll);
						}
						s = s + hs;// EGS4.format("-",67);
						EGS4.seqStr = s;
						if (EGS4.iprint > 1)
							eq.printSequence(EGS4.seqStr);

					}
					// "Start a new page"
					// WRITE (IOUT, *) '\f'; "page break"
					PGTHROW = 15;
					// WRITE(IOUT, 400) ' '; call egs_fdate(iout);
					// write(iout,'(//)');
					if (ROT) {
						// WRITE(IOUT, 93) TITLE;
						s = EGS4.format("", 19);
						EGS4.seqStr = s + "ZONAL MATERIAL GRID: ROTATED";
						if (EGS4.iprint > 1)
							eq.printSequence(EGS4.seqStr);
						s = EGS4.format("", 21);
						EGS4.seqStr = s + "**************************";
						if (EGS4.iprint > 1)
							eq.printSequence(EGS4.seqStr);
						EGS4.seqStr = "\n";
						if (EGS4.iprint > 1)
							eq.printSequence(EGS4.seqStr);
						s = EGS4.format("", 4);
						EGS4.seqStr = s
								+ "/X/Y/Z/MED :  X = ' ' IS DEFAULT: OPTION NOT USED";
						if (EGS4.iprint > 1)
							eq.printSequence(EGS4.seqStr);
						s = EGS4.format("", 18);
						EGS4.seqStr = s + "X = 'D' IF DOSE SCORING REGION";
						if (EGS4.iprint > 1)
							eq.printSequence(EGS4.seqStr);
						s = EGS4.format("", 18);
						EGS4.seqStr = s + "X = 'C' IF CAVITY REGION";
						if (EGS4.iprint > 1)
							eq.printSequence(EGS4.seqStr);
						s = EGS4.format("", 18);
						EGS4.seqStr = s + "X = 'S' IF SPR SCORING REGION";
						if (EGS4.iprint > 1)
							eq.printSequence(EGS4.seqStr);
						s = EGS4.format("", 18);
						EGS4.seqStr = s + "Y = 'T' IF TRACKING REGION";
						if (EGS4.iprint > 1)
							eq.printSequence(EGS4.seqStr);
						s = EGS4.format("", 18);
						EGS4.seqStr = s + "Z = 'A' IF TOTALLY ABSORBING REGION";
						if (EGS4.iprint > 1)
							eq.printSequence(EGS4.seqStr);
						s = EGS4.format("", 16);
						EGS4.seqStr = s
								+ "MED = MEDIUM NAME, 11 CHARACTER ABREVIATION";
						if (EGS4.iprint > 1)
							eq.printSequence(EGS4.seqStr);

					} else {
						// WRITE(IOUT, 94) TITLE;
						s = EGS4.format("", 19);
						EGS4.seqStr = s + "ZONAL MATERIAL GRID: NON-ROTATED";
						if (EGS4.iprint > 1)
							eq.printSequence(EGS4.seqStr);
						s = EGS4.format("", 21);
						EGS4.seqStr = s + "******************************";
						if (EGS4.iprint > 1)
							eq.printSequence(EGS4.seqStr);
						EGS4.seqStr = "\n";
						if (EGS4.iprint > 1)
							eq.printSequence(EGS4.seqStr);
						s = EGS4.format("", 4);
						EGS4.seqStr = s
								+ "/X/Y/Z/MED :  X = ' ' IS DEFAULT: OPTION NOT USED";
						if (EGS4.iprint > 1)
							eq.printSequence(EGS4.seqStr);
						s = EGS4.format("", 18);
						EGS4.seqStr = s + "  = 'D' IF DOSE SCORING REGION";
						if (EGS4.iprint > 1)
							eq.printSequence(EGS4.seqStr);
						s = EGS4.format("", 18);
						EGS4.seqStr = s + "  = 'C' IF CAVITY REGION";
						if (EGS4.iprint > 1)
							eq.printSequence(EGS4.seqStr);
						s = EGS4.format("", 18);
						EGS4.seqStr = s + "  = 'S' IF SPR SCORING REGION";
						if (EGS4.iprint > 1)
							eq.printSequence(EGS4.seqStr);
						s = EGS4.format("", 18);
						EGS4.seqStr = s + "Y = 'T' IF TRACKING REGION";
						if (EGS4.iprint > 1)
							eq.printSequence(EGS4.seqStr);
						s = EGS4.format("", 18);
						EGS4.seqStr = s + "Z = 'A' IF TOTALLY ABSORBING REGION";
						if (EGS4.iprint > 1)
							eq.printSequence(EGS4.seqStr);
						s = EGS4.format("", 16);
						EGS4.seqStr = s
								+ "MED = MEDIUM NAME, 11 CHARACTER ABREVIATION";
						if (EGS4.iprint > 1)
							eq.printSequence(EGS4.seqStr);

					}

					// WRITE (IOUT, 91) (RADIAL_BINS(IX), IX=HSET,HSET+FMT);
					for (int IX = HSET; IX <= HSET + FMT; IX++) {
						if (IX == HSET) {
							s = EGS4.format("", 7);
						} else if (IX == HSET + 1) {
							ll = 30 - ll;
							if (ll > 0) {
								s = EGS4.format("", ll);
							} else {
								s = EGS4.format("", 30);
							}
						} else if (IX == HSET + 2) {
							ll = 53 - ll;
							if (ll > 0) {
								s = EGS4.format("", ll);
							} else {
								s = EGS4.format("", 53);
							}
						} else if (IX == HSET + 3) {
							ll = 70 - ll;
							if (ll > 0) {
								s = EGS4.format("", ll);
							} else {
								s = EGS4.format("", 70);
							}
						}

						// ll=s.length();ll=60-ll;s=s+EGS4.format("",ll);
						// EGS4.seqStr=s+EGS4.format(RADIAL_BINS[IX-1],9,true);
						// if(EGS4.iprint>1)
						// eq.printSequence(EGS4.seqStr);
						if (IX == HSET) {
							EGS4.seqStr = s
									+ EGS4.format(RADIAL_BINS[IX - 1], 9, true);
							ll = EGS4.seqStr.length();
						} else {
							EGS4.seqStr = EGS4.seqStr + s
									+ EGS4.format(RADIAL_BINS[IX - 1], 9, true);
							ll = EGS4.seqStr.length();
						}
					}
					// EGS4.seqStr=s+EGS4.format(RADIAL_BINS[IX-1],9,true);
					if (EGS4.iprint > 1)
						eq.printSequence(EGS4.seqStr);// auto salt at new line

				}
				// "The following statement delays the page throw until the "
				// "last line of the grid can be printed"
				else {
					DLYPT = HSET;
				}
			}

		} // "end DO IZ"
			// "This is for the last vertical bin"
		if (FMT == 1) {
			hs = "";
			for (int i = 1; i <= 23; i++) {
				hs = hs + "-";
			}
			s = EGS4.format("", 1) + EGS4.format(DEPTH_BINS[DEEPNUM], 9, true);
			ll = s.length();
			ll = 11 - ll;
			if (ll > 0) {
				s = s + EGS4.format("", ll);
			}
			s = s + hs;// EGS4.format("-",23);
			EGS4.seqStr = s;
			if (EGS4.iprint > 1)
				eq.printSequence(EGS4.seqStr);

			// WRITE (IOUT, 10) DEPTH_BINS(DEEPNUM+1);
			// 10 FORMAT (1X, F9.4, T11, 23 ('-'));
		}
		if (FMT == 2) {
			hs = "";
			for (int i = 1; i <= 45; i++) {
				hs = hs + "-";
			}
			s = EGS4.format("", 1) + EGS4.format(DEPTH_BINS[DEEPNUM], 9, true);
			ll = s.length();
			ll = 11 - ll;
			if (ll > 0) {
				s = s + EGS4.format("", ll);
			}
			s = s + hs;// EGS4.format("-",45);
			EGS4.seqStr = s;
			if (EGS4.iprint > 1)
				eq.printSequence(EGS4.seqStr);

			// WRITE (IOUT, 11) DEPTH_BINS(DEEPNUM+1);
			// 11 FORMAT (1X, F9.4, T11, 45 ('-'));

		}
		if (FMT == 3) {
			hs = "";
			for (int i = 1; i <= 67; i++) {
				hs = hs + "-";
			}
			s = EGS4.format("", 1) + EGS4.format(DEPTH_BINS[DEEPNUM], 9, true);
			ll = s.length();
			ll = 11 - ll;
			if (ll > 0) {
				s = s + EGS4.format("", ll);
			}
			s = s + hs;// EGS4.format("-",67);
			EGS4.seqStr = s;
			if (EGS4.iprint > 1)
				eq.printSequence(EGS4.seqStr);

			// WRITE (IOUT, 12) DEPTH_BINS(DEEPNUM+1);
			// 12 FORMAT (1X, F9.4, T11, 67 ('-'));
		}
		// "For the delayed page throw"
		if ((DLYPT == HSET) && (RADNUM - HSET > 2)) {
			// "Start a new page"->I DO NOT CARE!
			// WRITE (IOUT, *) '\f'; "page break"
			PGTHROW = 17;
			// WRITE(IOUT, 400) ' '; call egs_fdate(iout); write(iout,'(//)');
			if (ROT) {
				s = EGS4.format("", 19);
				EGS4.seqStr = s + "ZONAL MATERIAL GRID: ROTATED";
				if (EGS4.iprint > 1)
					eq.printSequence(EGS4.seqStr);
				s = EGS4.format("", 21);
				EGS4.seqStr = s + "**************************";
				if (EGS4.iprint > 1)
					eq.printSequence(EGS4.seqStr);
				EGS4.seqStr = "\n";
				if (EGS4.iprint > 1)
					eq.printSequence(EGS4.seqStr);
				s = EGS4.format("", 4);
				EGS4.seqStr = s
						+ "/X/Y/Z/MED :  X = ' ' IS DEFAULT: OPTION NOT USED";
				if (EGS4.iprint > 1)
					eq.printSequence(EGS4.seqStr);
				s = EGS4.format("", 18);
				EGS4.seqStr = s + "X = 'D' IF DOSE SCORING REGION";
				if (EGS4.iprint > 1)
					eq.printSequence(EGS4.seqStr);
				s = EGS4.format("", 18);
				EGS4.seqStr = s + "X = 'C' IF CAVITY REGION";
				if (EGS4.iprint > 1)
					eq.printSequence(EGS4.seqStr);
				s = EGS4.format("", 18);
				EGS4.seqStr = s + "X = 'S' IF SPR SCORING REGION";
				if (EGS4.iprint > 1)
					eq.printSequence(EGS4.seqStr);
				s = EGS4.format("", 18);
				EGS4.seqStr = s + "Y = 'T' IF TRACKING REGION";
				if (EGS4.iprint > 1)
					eq.printSequence(EGS4.seqStr);
				s = EGS4.format("", 18);
				EGS4.seqStr = s + "Z = 'A' IF TOTALLY ABSORBING REGION";
				if (EGS4.iprint > 1)
					eq.printSequence(EGS4.seqStr);
				s = EGS4.format("", 16);
				EGS4.seqStr = s + "MED = MEDIUM NAME, 11 CHARACTER ABREVIATION";
				if (EGS4.iprint > 1)
					eq.printSequence(EGS4.seqStr);

				// WRITE(IOUT, 93) TITLE;
			} else {
				s = EGS4.format("", 19);
				EGS4.seqStr = s + "ZONAL MATERIAL GRID: NON-ROTATED";
				if (EGS4.iprint > 1)
					eq.printSequence(EGS4.seqStr);
				s = EGS4.format("", 21);
				EGS4.seqStr = s + "******************************";
				if (EGS4.iprint > 1)
					eq.printSequence(EGS4.seqStr);
				EGS4.seqStr = "\n";
				if (EGS4.iprint > 1)
					eq.printSequence(EGS4.seqStr);
				s = EGS4.format("", 4);
				EGS4.seqStr = s
						+ "/X/Y/Z/MED :  X = ' ' IS DEFAULT: OPTION NOT USED";
				if (EGS4.iprint > 1)
					eq.printSequence(EGS4.seqStr);
				s = EGS4.format("", 18);
				EGS4.seqStr = s + "  = 'D' IF DOSE SCORING REGION";
				if (EGS4.iprint > 1)
					eq.printSequence(EGS4.seqStr);
				s = EGS4.format("", 18);
				EGS4.seqStr = s + "  = 'C' IF CAVITY REGION";
				if (EGS4.iprint > 1)
					eq.printSequence(EGS4.seqStr);
				s = EGS4.format("", 18);
				EGS4.seqStr = s + "  = 'S' IF SPR SCORING REGION";
				if (EGS4.iprint > 1)
					eq.printSequence(EGS4.seqStr);
				s = EGS4.format("", 18);
				EGS4.seqStr = s + "Y = 'T' IF TRACKING REGION";
				if (EGS4.iprint > 1)
					eq.printSequence(EGS4.seqStr);
				s = EGS4.format("", 18);
				EGS4.seqStr = s + "Z = 'A' IF TOTALLY ABSORBING REGION";
				if (EGS4.iprint > 1)
					eq.printSequence(EGS4.seqStr);
				s = EGS4.format("", 16);
				EGS4.seqStr = s + "MED = MEDIUM NAME, 11 CHARACTER ABREVIATION";
				if (EGS4.iprint > 1)
					eq.printSequence(EGS4.seqStr);

				// WRITE(IOUT, 94) TITLE;
			}
		}
		PGTHROW = PGTHROW + 1;
		// }

	}

	// public static void ZONEGRID(int NRADIAL,int NDEPTH,int NRMIN,int
	// NZMIN,int NZ,
	// double[][][] RESULTS,double[][][]UNCRT,int NCOMP,double[] RADIAL_BINS,
	// double[] DEPTH_BINS),String[] LABELS,String[] EXPLANATIONS)
	/**
	 * This routine provides a generalized grid output. 
	 * @param NRADIAL NRADIAL : number of radial results to display
	 * @param NDEPTH NDEPTH : number of depth results to display
	 * @param NRMIN NRMIN : index of min radius of results zone
	 * @param NZMIN NZMIN : index of min Z plane of results zone
	 * @param NZ NZ : total number of slabs in the geometry
	 */
	public static void ZONEGRID(int NRADIAL, int NDEPTH, int NRMIN, int NZMIN,
			int NZ) {
		// "
		// " DESCRIPTION : Routine provides a generalized grid output of
		// " OF : values. The grid will be rotated if radnum>3
		// " THE : and deepnum<=3. Also, if 10% of all the uncertainties
		// " SUBROUTINE : are less than 0.3% (this is the value of the
		// " parameter 'CUTOFF') an extra significant figure will
		// " appear in the grid output.
		// "
		// " Ensure that the macro $MAXRZ has been declared.
		// " \ One must add the macro '$MAXCMPTS'(usually=$MAXIT)
		// /" FOR INSTALLING \ to the program in which this subroutine is
		// called.
		// " THIS > Be careful that all arrays passed to the subroutine
		// /" SUBROUTINE / (especially 'RESULTS()' and 'UNCRT()') have exactly
		// " / the same dimensions as used here.
		// "
		// " Description of the arguments:
		// " NRADIAL : number of radial results to display
		// " NDEPTH : number of depth results to display
		// " NRMIN : index of min radius of results zone
		// " NZMIN : index of min Z plane of results zone
		// " NZ : total number of slabs in the geometry
		// " RESULTS : array of results
		// " UNCRT : corresponding array of uncertainties(max at 99.9%)
		// " NCOMP : number of components of the results
		// " RADIAL_BINS : array of outer boundaries of each radial bin(1st
		// " radius=0)
		// " DEPTH_BINS : array of depth boundaries (1st plane not always 0)
		// " LABELS : array of length ncomp, the label for each component
		// " EXPLANATIONS: explanations for the labels
		// "
		// "***************************************************************************
		// ;IMPLICIT NONE;
		//zoneB = true;
		String s = "";

		// $LOGICAL SFIG; "controls significant figures on UNCRT()"
		// "true => use 1 significant figure"
		// "false => use 2 significant figures"
		double CUTOFF = 0.3; // "changes SFIG if 10% of uncertainties are less"
		// // "than 0.2%"
		// $REAL RESULTS($MAXZREG,$MAXRADII,$MAXCMPTS),
		// UNCRT($MAXZREG,$MAXRADII,$MAXCMPTS),
		// RADIAL_BINS($MAXRZ),DEPTH_BINS($MAXRZ),
		double TMP2 = 0.0;
		double TMP3 = 0.0;
		;
		// CHARACTER*60 EXPLANATIONS($MAXCMPTS);
		// CHARACTER*4 LABELS($MAXCMPTS);
		// $INTEGER IOUT, ICOMP, IX, IZ, HSET, PGTHROW, DLYPT,
		// COUNT,NZMIN,NRMIN,NZ;

		ROT = false; // "For a rotated grid (.FALSE. = off)"
		SFIG = true; // "Default sig. fig. option off"
		if (NCOMP > 4) {
			PGTHROW = 20;
		} else {
			PGTHROW = 17;
		}
		// "For printing: throw page at 66 lines(16 1st for title)"
		// IOUT=1; "Designates output to file fort.1, the *.egs4lst file"
		// "Change variable names so any changes are not passed back to MAIN prog"
		RADNUM = NRADIAL;
		DEEPNUM = NDEPTH;
		COUNT = 0;
		for (int IX = 1; IX <= RADNUM; IX++) {// "Check for needed precision of uncertainties"
			for (int IZ = 1; IZ <= DEEPNUM; IZ++) {
				for (int ICOMP = 1; ICOMP <= NCOMP; ICOMP++) {
					if (UNCRT[IZ - 1][IX - 1][ICOMP - 1] < CUTOFF) {
						COUNT = COUNT + 1;
					}
				}
			}
		}
		// "IF better precision in UNCRT() is needed, change the grid format."
		if (COUNT / (RADNUM * DEEPNUM * NCOMP) > 0.1) {
			SFIG = false;
		}
		// "Check for rotation of the grid"
		if ((DEEPNUM <= 3) && (RADNUM > 3)) {
			ROT = true;
			for (int IX = 1; IX <= RADNUM; IX++)// DO IX=1, RADNUM [
			{
				if (IX <= DEEPNUM + 1) {
					TMP3 = RADIAL_BINS[IX - 1];
					RADIAL_BINS[IX - 1] = DEPTH_BINS[IX - 1];
					DEPTH_BINS[IX - 1] = TMP3;
				} else {
					DEPTH_BINS[IX - 1] = RADIAL_BINS[IX - 1];
				}
				if (IX <= DEEPNUM) {
					for (int IZ = IX + 1; IZ <= DEEPNUM; IZ++)// DO
																// IZ=IX+1,DEEPNUM[
					{
						for (int ICOMP = 1; ICOMP <= NCOMP; ICOMP++)// DO
																	// ICOMP=1,NCOMP[
						{
							TMP2 = RESULTS[IZ - 1][IX - 1][ICOMP - 1];
							RESULTS[IZ - 1][IX - 1][ICOMP - 1] = RESULTS[IX - 1][IZ - 1][ICOMP - 1];
							RESULTS[IX - 1][IZ - 1][ICOMP - 1] = TMP2;
							TMP2 = UNCRT[IZ - 1][IX - 1][ICOMP - 1];
							UNCRT[IZ - 1][IX - 1][ICOMP - 1] = UNCRT[IX - 1][IZ - 1][ICOMP - 1];
							UNCRT[IX - 1][IZ - 1][ICOMP - 1] = TMP2;
						}
					}
				} else {
					for (int IZ = 1; IZ <= DEEPNUM; IZ++)// DO IZ=1,DEEPNUM[
					{
						for (int ICOMP = 1; ICOMP <= NCOMP; ICOMP++)// DO
																	// ICOMP=1,NCOMP[
						{
							RESULTS[IX - 1][IZ - 1][ICOMP - 1] = RESULTS[IZ - 1][IX - 1][ICOMP - 1];
							UNCRT[IX - 1][IZ - 1][ICOMP - 1] = UNCRT[IZ - 1][IX - 1][ICOMP - 1];
						}
					}
				}
			}// "end DO IX"
			DEPTH_BINS[RADNUM] = RADIAL_BINS[RADNUM];
			TMP1 = RADNUM;
			RADNUM = DEEPNUM;
			DEEPNUM = TMP1;
		}// "end check for rotation"
			// WRITE (IOUT, *) '\f'; "Page break"
			// WRITE(IOUT, 400) ' '; call egs_fdate(iout); write(iout,'(//)');
		if (ROT) {
			// WRITE(IOUT, 93) TITLE;
			// 93 FORMAT (' ',79A1 //
			// T19, 27H ZONAL OUTPUT GRID: ROTATED /
			// T20, '**************************');

			s = EGS4.format("", 19);
			EGS4.seqStr = s + "ZONAL OUTPUT GRID: ROTATED";
			if (EGS4.iprint > 0)
				eq.printSequence(EGS4.seqStr);
			s = EGS4.format("", 20);
			EGS4.seqStr = s + "**************************";
			if (EGS4.iprint > 0)
				eq.printSequence(EGS4.seqStr);
		} else {
			s = EGS4.format("", 19);
			EGS4.seqStr = s + "ZONAL OUTPUT GRID: NON-ROTATED";
			if (EGS4.iprint > 0)
				eq.printSequence(EGS4.seqStr);
			s = EGS4.format("", 20);
			EGS4.seqStr = s + "******************************";
			if (EGS4.iprint > 0)
				eq.printSequence(EGS4.seqStr);

			// WRITE(IOUT, 94) TITLE;
		}
		for (int ICOMP = 1; ICOMP <= NCOMP; ICOMP++) {
			s = EGS4.format("", 10) + EGS4.format(LABELS[ICOMP - 1], 4) + "  "
					+ EGS4.format(EXPLANATIONS[ICOMP - 1], 60);
			EGS4.seqStr = s;
			if (EGS4.iprint > 0)
				eq.printSequence(EGS4.seqStr);
			// WRITE(IOUT, 95) LABELS(ICOMP), EXPLANATIONS(ICOMP);
		}
		// "Make the grids"
		for (HSET = 1; HSET <= RADNUM; HSET = HSET + 3) {
			// System.out.println("@@@ "+NCOMP);
			// $MKGRID;
			MKGRIDZ(NRMIN, NZMIN, NZ, RESULTS, UNCRT);
		}// "Feed in horiz. sets of three"
			// "FORMATS"
			// 1 FORMAT (T11, '|',
			// A4, 1PE10.3, '+-', 0PF4.1, '%', '|');
			// 2 FORMAT (T11, '|',
			// 1 (A4, 1PE10.3, '+-', 0PF4.1, '%', ' |'),
			// A4, 1PE10.3, '+-', 0PF4.1, '%', '|');
			// 3 FORMAT (T11, '|',
			// 2 (A4, 1PE10.3, '+-', 0PF4.1, '%', ' |'),
			// A4, 1PE10.3, '+-', 0PF4.1, '%', '|');
			// 4 FORMAT (T11, '|',
			// A4, 1PE10.3, '+-', 0PF5.2, '%', '|');
			// 5 FORMAT (T11, '|',
			// 1 (A4, 1PE10.3, '+-', 0PF5.2, '%', ' |'),
			// A4, 1PE10.3, '+-', 0PF5.2, '%', '|');
			// 6 FORMAT (T10, '|',
			// 2 (A4, 1PE10.3, '+-', 0PF5.2, '%', '|'),
			// A4, 1PE10.3, '+-', 0PF5.2, '%', '|');
			// 10 FORMAT (1X, F9.4, T11, 1 (23 ('-')));
			// 11 FORMAT (1X, F9.4, T11, 2 (23 ('-')));
			// 12 FORMAT (1X, F9.4, T11, 3 (23 ('-')));
			// 13 FORMAT (T11, '|',
			// 2 ('IRL',I3,1X,'IZ ',I3,1X,'IX ',I3,1X,' |'),
			// 'IRL',I3,1X,'IZ ',I3,1X,'IX ',I3,1X,'|');
			// 14 FORMAT (T11, '|',
			// 'IRL',I3,1X,'IZ ',I3,1X,'IX ',I3,1X,' |',
			// 'IRL',I3,1X,'IZ ',I3,1X,'IX ',I3,1X,'|');
			// 15 FORMAT (T11, '|',
			// 'IRL',I3,1X,'IZ ',I3,1X,'IX ',I3,1X,'|');
			// 16 FORMAT (T10, '|',
			// 2 ('IRL',I3,1X,'IZ ',I3,1X,'IX ',I3,1X,' |'),
			// 'IRL',I3,1X,'IZ ',I3,1X,'IX ',I3,1X,' |');
			// 17 FORMAT (T11, '|',
			// 1 ('IRL',I3,1X,'IZ ',I3,1X,'IX ',I3,1X,' |'),
			// 'IRL',I3,1X,'IZ ',I3,1X,'IX ',I3,1X,' |');
			// 18 FORMAT (T11, '|',
			// 'IRL',I3,1X,'IZ ',I3,1X,'IX ',I3,1X,' |');
			// 91 FORMAT (/ T7,F9.4,T30,F9.4,T53,F9.4,T70,F9.4);
			// 93 FORMAT (' ',79A1 //
			// T19, 27H ZONAL OUTPUT GRID: ROTATED /
			// T20, '**************************');
			// 94 FORMAT (' ',79A1 //
			// T19, 31H ZONAL OUTPUT GRID: NON-ROTATED /
			// T20, '******************************');
			// 95 FORMAT (T14, A4, T19, A60);
			// 400 FORMAT (T54,a1,$);

		return;
	}

	/**
	 * Internally used. Called by ZONEGRID.
	 * @param NRMIN NRMIN
	 * @param NZMIN NZMIN
	 * @param NZ NZ
	 * @param RESULTS RESULTS
	 * @param UNCRT UNCRT
	 */
	private static void MKGRIDZ(int NRMIN, int NZMIN, int NZ,
			double[][][] RESULTS, double[][][] UNCRT) {
		String s = "";
		String hs = "";
		int ll = 0;
		// "Set up the formatting indicator FMT"
		if (RADNUM - HSET > 1) {
			FMT = 3;
		} else {
			if (RADNUM % 3 == 1) {
				FMT = 1;
			}
			if (RADNUM % 3 == 2) {
				FMT = 2;
			}

			// IF (MOD(RADNUM,3) = 1) {FMT=1;}
			// IF (MOD(RADNUM,3) = 2) [FMT=2;}
		}
		// "Horizontal axis indicators"
		// WRITE (IOUT, 91) (RADIAL_BINS(IX), IX=HSET,HSET+FMT);
		// 91 FORMAT (/ T7,F9.4,T30,F9.4,T53,F9.4,T70,F9.4);
		for (int IX = HSET; IX <= HSET + FMT; IX++) {
			if (IX == HSET) {
				s = EGS4.format("", 7);
			} else if (IX == HSET + 1) {
				ll = 30 - ll;
				if (ll > 0) {
					s = EGS4.format("", ll);
				} else {
					s = EGS4.format("", 30);
				}
			} else if (IX == HSET + 2) {
				ll = 53 - ll;
				if (ll > 0) {
					s = EGS4.format("", ll);
				} else {
					s = EGS4.format("", 53);
				}
			} else if (IX == HSET + 3) {
				ll = 70 - ll;
				if (ll > 0) {
					s = EGS4.format("", ll);
				} else {
					s = EGS4.format("", 70);
				}
			}

			// ll=s.length();ll=60-ll;s=s+EGS4.format("",ll);
			// EGS4.seqStr=s+EGS4.format(RADIAL_BINS[IX-1],9,true);
			// if(EGS4.iprint>1)
			// eq.printSequence(EGS4.seqStr);
			if (IX == HSET) {
				EGS4.seqStr = s + EGS4.format(RADIAL_BINS[IX - 1], 9, true);
				ll = EGS4.seqStr.length();
			} else {
				EGS4.seqStr = EGS4.seqStr + s
						+ EGS4.format(RADIAL_BINS[IX - 1], 9, true);
				ll = EGS4.seqStr.length();
			}
		}
		// EGS4.seqStr=s+EGS4.format(RADIAL_BINS[IX-1],9,true);
		if (EGS4.iprint > 0)
			eq.printSequence(EGS4.seqStr);// auto salt at new line

		PGTHROW = PGTHROW + 1;
		for (int IZ = 1; IZ <= DEEPNUM; IZ++) {
			// "This is for IRL, IZ and IX"
			if (ROT) {
				REGNUM = (IZ + NRMIN - 1) * NZ + HSET + NZMIN;
			} else {
				REGNUM = (HSET + NRMIN - 1) * NZ + IZ + NZMIN;
			}
			if (FMT == 1) {
				// WRITE (IOUT, 10) DEPTH_BINS(IZ)
				// 10 FORMAT (1X, F9.4, T11, 1 (23 ('-')));
				hs = "";
				for (int i = 1; i <= 23; i++) {
					hs = hs + "-";
				}
				s = EGS4.format("", 1)
						+ EGS4.format(DEPTH_BINS[IZ - 1], 9, true);
				ll = s.length();
				ll = 11 - ll;
				if (ll > 0) {
					s = s + EGS4.format("", ll);
				}
				s = s + hs;// EGS4.format("-",23);
				EGS4.seqStr = s;
				if (EGS4.iprint > 0)
					eq.printSequence(EGS4.seqStr);

				if (SFIG) {
					if (ROT) {
						s = EGS4.format("", 11) + "|" + "IRL"
								+ EGS4.format(REGNUM, 3) + EGS4.format("", 1)
								+ "IZ " + EGS4.format(HSET + NZMIN - 1, 3)
								+ EGS4.format("", 1) + "IX "
								+ EGS4.format(IZ + NRMIN, 3)
								+ EGS4.format("", 1) + "|";
						EGS4.seqStr = s;
						if (EGS4.iprint > 0)
							eq.printSequence(EGS4.seqStr);

						// WRITE (IOUT, 15) REGNUM,HSET+NZMIN-1,IZ+NRMIN;
						// 15 FORMAT (T11, '|',
						// 'IRL',I3,1X,'IZ ',I3,1X,'IX ',I3,1X,'|');
					} else {
						s = EGS4.format("", 11) + "|" + "IRL"
								+ EGS4.format(REGNUM, 3) + EGS4.format("", 1)
								+ "IZ " + EGS4.format(IZ + NZMIN - 1, 3)
								+ EGS4.format("", 1) + "IX "
								+ EGS4.format(HSET + NRMIN, 3)
								+ EGS4.format("", 1) + "|";
						EGS4.seqStr = s;
						if (EGS4.iprint > 0)
							eq.printSequence(EGS4.seqStr);

						// WRITE (IOUT, 15) REGNUM,IZ+NZMIN-1,HSET+NRMIN;
					}
				} else {
					if (ROT) {
						s = EGS4.format("", 11) + "|" + "IRL"
								+ EGS4.format(REGNUM, 3) + EGS4.format("", 1)
								+ "IZ " + EGS4.format(HSET + NZMIN - 1, 3)
								+ EGS4.format("", 1) + "IX "
								+ EGS4.format(IZ + NRMIN, 3)
								+ EGS4.format("", 1) + " |";
						EGS4.seqStr = s;
						if (EGS4.iprint > 0)
							eq.printSequence(EGS4.seqStr);

						// WRITE (IOUT, 18) REGNUM,HSET+NZMIN-1,IZ+NRMIN;
						// 18 FORMAT (T11, '|',
						// 'IRL',I3,1X,'IZ ',I3,1X,'IX ',I3,1X,' |');

					} else {
						s = EGS4.format("", 11) + "|" + "IRL"
								+ EGS4.format(REGNUM, 3) + EGS4.format("", 1)
								+ "IZ " + EGS4.format(HSET + NZMIN - 1, 3)
								+ EGS4.format("", 1) + "IX "
								+ EGS4.format(IZ + NRMIN, 3)
								+ EGS4.format("", 1) + " |";
						EGS4.seqStr = s;
						if (EGS4.iprint > 0)
							eq.printSequence(EGS4.seqStr);

						// WRITE (IOUT, 18) REGNUM,IZ+NZMIN-1,HSET+NRMIN;
					}
				}
			}
			if (FMT == 2) {
				// WRITE (IOUT, 11) DEPTH_BINS(IZ);
				// 11 FORMAT (1X, F9.4, T11, 2 (23 ('-')));
				hs = "";
				for (int i = 1; i <= 46; i++) {
					hs = hs + "-";
				}
				s = EGS4.format("", 1)
						+ EGS4.format(DEPTH_BINS[IZ - 1], 9, true);
				ll = s.length();
				ll = 11 - ll;
				if (ll > 0) {
					s = s + EGS4.format("", ll);
				}
				s = s + hs;// EGS4.format("-",45);
				EGS4.seqStr = s;
				if (EGS4.iprint > 0)
					eq.printSequence(EGS4.seqStr);

				if (SFIG) {
					if (ROT) {
						s = EGS4.format("", 11) + "|" + "IRL"
								+ EGS4.format(REGNUM, 3) + EGS4.format("", 1)
								+ "IZ " + EGS4.format(HSET + NZMIN - 1, 3)
								+ EGS4.format("", 1) + "IX "
								+ EGS4.format(IZ + NRMIN, 3)
								+ EGS4.format("", 1) + " |" + "IRL"
								+ EGS4.format(REGNUM + 1, 3)
								+ EGS4.format("", 1) + "IZ "
								+ EGS4.format(HSET + NZMIN + 1, 3)
								+ EGS4.format("", 1) + "IX "
								+ EGS4.format(IZ + NRMIN, 3)
								+ EGS4.format("", 1) + "|";

						EGS4.seqStr = s;
						if (EGS4.iprint > 0)
							eq.printSequence(EGS4.seqStr);

						// WRITE(IOUT,14) REGNUM,HSET+NZMIN-1,IZ+NRMIN,
						// REGNUM+1,HSET+NZMIN,IZ+NRMIN;
						// 14 FORMAT (T11, '|',
						// 'IRL',I3,1X,'IZ ',I3,1X,'IX ',I3,1X,' |',
						// 'IRL',I3,1X,'IZ ',I3,1X,'IX ',I3,1X,'|');

					} else {
						s = EGS4.format("", 11) + "|" + "IRL"
								+ EGS4.format(REGNUM, 3) + EGS4.format("", 1)
								+ "IZ " + EGS4.format(IZ + NZMIN - 1, 3)
								+ EGS4.format("", 1) + "IX "
								+ EGS4.format(HSET + NRMIN, 3)
								+ EGS4.format("", 1) + " |" + "IRL"
								+ EGS4.format(REGNUM + NZ, 3)
								+ EGS4.format("", 1) + "IZ "
								+ EGS4.format(IZ + NZMIN - 1, 3)
								+ EGS4.format("", 1) + "IX "
								+ EGS4.format(HSET + NRMIN + 1, 3)
								+ EGS4.format("", 1) + "|";

						EGS4.seqStr = s;
						if (EGS4.iprint > 0)
							eq.printSequence(EGS4.seqStr);

						// WRITE(IOUT,14) REGNUM,IZ+NZMIN-1,HSET+NRMIN,
						// REGNUM+NZ,IZ+NZMIN-1,HSET+NRMIN+1;
					}
				} else {
					if (ROT) {
						s = EGS4.format("", 11) + "|" + "IRL"
								+ EGS4.format(REGNUM, 3) + EGS4.format("", 1)
								+ "IZ " + EGS4.format(HSET + NZMIN - 1, 3)
								+ EGS4.format("", 1) + "IX "
								+ EGS4.format(IZ + NRMIN, 3)
								+ EGS4.format("", 1) + "  |" + "IRL"
								+ EGS4.format(REGNUM + 1, 3)
								+ EGS4.format("", 1) + "IZ "
								+ EGS4.format(HSET + NZMIN, 3)
								+ EGS4.format("", 1) + "IX "
								+ EGS4.format(IZ + NRMIN, 3)
								+ EGS4.format("", 1) + " |";

						EGS4.seqStr = s;
						if (EGS4.iprint > 0)
							eq.printSequence(EGS4.seqStr);

						// WRITE(IOUT,17) REGNUM,HSET+NZMIN-1,IZ+NRMIN,
						// REGNUM+1,HSET+NZMIN, IZ+NRMIN;
						// 17 FORMAT (T11, '|',
						// 1 ('IRL',I3,1X,'IZ ',I3,1X,'IX ',I3,1X,' |'),
						// 'IRL',I3,1X,'IZ ',I3,1X,'IX ',I3,1X,' |');

					} else {
						s = EGS4.format("", 11) + "|" + "IRL"
								+ EGS4.format(REGNUM, 3) + EGS4.format("", 1)
								+ "IZ " + EGS4.format(IZ + NZMIN - 1, 3)
								+ EGS4.format("", 1) + "IX "
								+ EGS4.format(HSET + NRMIN, 3)
								+ EGS4.format("", 1) + "  |" + "IRL"
								+ EGS4.format(REGNUM + NZ, 3)
								+ EGS4.format("", 1) + "IZ "
								+ EGS4.format(IZ + NZMIN - 1, 3)
								+ EGS4.format("", 1) + "IX "
								+ EGS4.format(HSET + NRMIN + 1, 3)
								+ EGS4.format("", 1) + " |";

						EGS4.seqStr = s;
						if (EGS4.iprint > 0)
							eq.printSequence(EGS4.seqStr);

						// WRITE(IOUT,17) REGNUM,IZ+NZMIN-1,HSET+NRMIN,
						// REGNUM+NZ,IZ+NZMIN-1,HSET+NRMIN+1;
					}
				}
			}
			if (FMT == 3) {
				hs = "";
				for (int i = 1; i <= 69; i++) {
					hs = hs + "-";
				}
				s = EGS4.format("", 1)
						+ EGS4.format(DEPTH_BINS[IZ - 1], 9, true);
				ll = s.length();
				ll = 11 - ll;
				if (ll > 0) {
					s = s + EGS4.format("", ll);
				}
				s = s + hs;// EGS4.format("-",67);
				EGS4.seqStr = s;
				if (EGS4.iprint > 0)
					eq.printSequence(EGS4.seqStr);

				// WRITE (IOUT, 12) DEPTH_BINS(IZ);
				// 12 FORMAT (1X, F9.4, T11, 3 (23 ('-')));
				if (SFIG) {
					if (ROT) {
						s = EGS4.format("", 11) + "|" + "IRL"
								+ EGS4.format(REGNUM, 3) + EGS4.format("", 1)
								+ "IZ " + EGS4.format(HSET + NZMIN - 1, 3)
								+ EGS4.format("", 1) + "IX "
								+ EGS4.format(IZ + NRMIN, 3)
								+ EGS4.format("", 1) + " |" + "IRL"
								+ EGS4.format(REGNUM + 1, 3)
								+ EGS4.format("", 1) + "IZ "
								+ EGS4.format(HSET + NZMIN, 3)
								+ EGS4.format("", 1) + "IX "
								+ EGS4.format(IZ + NRMIN, 3)
								+ EGS4.format("", 1) + " |" + "IRL"
								+ EGS4.format(REGNUM + 2, 3)
								+ EGS4.format("", 1) + "IZ "
								+ EGS4.format(HSET + NZMIN + 1, 3)
								+ EGS4.format("", 1) + "IX "
								+ EGS4.format(IZ + NRMIN, 3)
								+ EGS4.format("", 1) + "|";

						EGS4.seqStr = s;
						if (EGS4.iprint > 0)
							eq.printSequence(EGS4.seqStr);

						// WRITE (IOUT,13)REGNUM,HSET+NZMIN-1,IZ+NRMIN,
						// REGNUM+1,HSET+NZMIN,IZ+NRMIN,
						// REGNUM+2,HSET+NZMIN+1,IZ+NRMIN;
						// 13 FORMAT (T11, '|',
						// 2 ('IRL',I3,1X,'IZ ',I3,1X,'IX ',I3,1X,' |'),
						// 'IRL',I3,1X,'IZ ',I3,1X,'IX ',I3,1X,'|');

					} else {
						s = EGS4.format("", 11) + "|" + "IRL"
								+ EGS4.format(REGNUM, 3) + EGS4.format("", 1)
								+ "IZ " + EGS4.format(IZ + NZMIN - 1, 3)
								+ EGS4.format("", 1) + "IX "
								+ EGS4.format(HSET + NRMIN, 3)
								+ EGS4.format("", 1) + " |" + "IRL"
								+ EGS4.format(REGNUM + NZ, 3)
								+ EGS4.format("", 1) + "IZ "
								+ EGS4.format(IZ + NZMIN - 1, 3)
								+ EGS4.format("", 1) + "IX "
								+ EGS4.format(HSET + NRMIN + 1, 3)
								+ EGS4.format("", 1) + " |" + "IRL"
								+ EGS4.format(REGNUM + NZ * 2, 3)
								+ EGS4.format("", 1) + "IZ "
								+ EGS4.format(IZ + NZMIN - 1, 3)
								+ EGS4.format("", 1) + "IX "
								+ EGS4.format(HSET + NRMIN + 2, 3)
								+ EGS4.format("", 1) + "|";

						EGS4.seqStr = s;
						if (EGS4.iprint > 0)
							eq.printSequence(EGS4.seqStr);

						// WRITE (IOUT, 13) REGNUM,IZ+NZMIN-1,HSET+NRMIN,
						// REGNUM+NZ,IZ+NZMIN-1,HSET+NRMIN+1,
						// REGNUM+NZ*2,IZ+NZMIN-1,HSET+NRMIN+2;
					}
				}// "end SFIG true block"
				else {// "not SFIG"
					if (ROT) {
						s = EGS4.format("", 10) + "|" + "IRL"
								+ EGS4.format(REGNUM, 3) + EGS4.format("", 1)
								+ "IZ " + EGS4.format(HSET + NZMIN - 1, 3)
								+ EGS4.format("", 1) + "IX "
								+ EGS4.format(IZ + NRMIN, 3)
								+ EGS4.format("", 1) + " |" + "IRL"
								+ EGS4.format(REGNUM + 1, 3)
								+ EGS4.format("", 1) + "IZ "
								+ EGS4.format(HSET + NZMIN, 3)
								+ EGS4.format("", 1) + "IX "
								+ EGS4.format(IZ + NRMIN, 3)
								+ EGS4.format("", 1) + " |" + "IRL"
								+ EGS4.format(REGNUM + 2, 3)
								+ EGS4.format("", 1) + "IZ "
								+ EGS4.format(HSET + NZMIN + 1, 3)
								+ EGS4.format("", 1) + "IX "
								+ EGS4.format(IZ + NRMIN, 3)
								+ EGS4.format("", 1) + " |";

						EGS4.seqStr = s;
						if (EGS4.iprint > 0)
							eq.printSequence(EGS4.seqStr);

						// WRITE (IOUT, 16)REGNUM,HSET+NZMIN-1,IZ+NRMIN,
						// REGNUM+1,HSET+NZMIN,IZ+NRMIN,
						// REGNUM+2,HSET+NZMIN+1,IZ+NRMIN;
						// 16 FORMAT (T10, '|',
						// 2 ('IRL',I3,1X,'IZ ',I3,1X,'IX ',I3,1X,' |'),
						// 'IRL',I3,1X,'IZ ',I3,1X,'IX ',I3,1X,' |');

					} else {
						s = EGS4.format("", 10) + "|" + "IRL"
								+ EGS4.format(REGNUM, 3) + EGS4.format("", 1)
								+ "IZ " + EGS4.format(IZ + NZMIN - 1, 3)
								+ EGS4.format("", 1) + "IX "
								+ EGS4.format(HSET + NRMIN, 3)
								+ EGS4.format("", 1) + " |" + "IRL"
								+ EGS4.format(REGNUM + NZ, 3)
								+ EGS4.format("", 1) + "IZ "
								+ EGS4.format(IZ + NZMIN - 1, 3)
								+ EGS4.format("", 1) + "IX "
								+ EGS4.format(HSET + NRMIN + 1, 3)
								+ EGS4.format("", 1) + " |" + "IRL"
								+ EGS4.format(REGNUM + NZ * 2, 3)
								+ EGS4.format("", 1) + "IZ "
								+ EGS4.format(IZ + NZMIN - 1, 3)
								+ EGS4.format("", 1) + "IX "
								+ EGS4.format(HSET + NRMIN + 2, 3)
								+ EGS4.format("", 1) + " |";

						EGS4.seqStr = s;
						if (EGS4.iprint > 0)
							eq.printSequence(EGS4.seqStr);

						// WRITE (IOUT, 16)REGNUM,IZ+NZMIN-1,HSET+NRMIN,
						// REGNUM+NZ,IZ+NZMIN-1,HSET+NRMIN+1,
						// REGNUM+NZ*2,IZ+NZMIN-1,HSET+NRMIN+2;

					}
				}
			}// "end FMT=3 block"
			PGTHROW = PGTHROW + 2;
			// "Main part of the grid"
			for (int ICOMP = 1; ICOMP <= NCOMP; ICOMP++) {
				if (SFIG) {
					if (FMT == 3) {
						s = EGS4.format("", 11)
								+ "|"
								+ EGS4.format(LABELS[ICOMP - 1], 4)
								+ EGS4.format(
										RESULTS[IZ - 1][HSET - 1][ICOMP - 1],
										10)
								+ "+-"
								+ EGS4.format1(
										UNCRT[IZ - 1][HSET - 1][ICOMP - 1], 4,
										true)
								+ "%"
								+ " |"
								+ EGS4.format(LABELS[ICOMP - 1], 4)
								+ EGS4.format(RESULTS[IZ - 1][HSET][ICOMP - 1],
										10)
								+ "+-"
								+ EGS4.format1(UNCRT[IZ - 1][HSET][ICOMP - 1],
										4, true)
								+ "%"
								+ " |"
								+ EGS4.format(LABELS[ICOMP - 1], 4)
								+ EGS4.format(
										RESULTS[IZ - 1][HSET + 1][ICOMP - 1],
										10)
								+ "+-"
								+ EGS4.format1(
										UNCRT[IZ - 1][HSET + 1][ICOMP - 1], 4,
										true) + "%" + "|";

						EGS4.seqStr = s;
						if (EGS4.iprint > 0)
							eq.printSequence(EGS4.seqStr);

						// WRITE(IOUT, 3) LABELS(ICOMP),RESULTS(IZ,HSET,ICOMP),
						// UNCRT(IZ,HSET,ICOMP),LABELS(ICOMP),
						// RESULTS(IZ,HSET+1,ICOMP),UNCRT(IZ,HSET+1,ICOMP),
						// LABELS(ICOMP),RESULTS(IZ,HSET+2,ICOMP),
						// UNCRT(IZ,HSET+2,ICOMP);
						// 3 FORMAT (T11, '|',
						// 2 (A4, 1PE10.3, '+-', 0PF4.1, '%', ' |'),
						// A4, 1PE10.3, '+-', 0PF4.1, '%', '|');

					}
					if (FMT == 2) {
						s = EGS4.format("", 11)
								+ "|"
								+ EGS4.format(LABELS[ICOMP - 1], 4)
								+ EGS4.format(
										RESULTS[IZ - 1][HSET - 1][ICOMP - 1],
										10)
								+ "+-"
								+ EGS4.format1(
										UNCRT[IZ - 1][HSET - 1][ICOMP - 1], 4,
										true)
								+ "%"
								+ " |"
								+ EGS4.format(LABELS[ICOMP - 1], 4)
								+ EGS4.format(RESULTS[IZ - 1][HSET][ICOMP - 1],
										10)
								+ "+-"
								+ EGS4.format1(UNCRT[IZ - 1][HSET][ICOMP - 1],
										4, true) + "%" + "|";

						EGS4.seqStr = s;
						if (EGS4.iprint > 0)
							eq.printSequence(EGS4.seqStr);

						// WRITE(IOUT, 2) LABELS(ICOMP),RESULTS(IZ,HSET,ICOMP),
						// UNCRT(IZ,HSET,ICOMP),LABELS(ICOMP),
						// RESULTS(IZ,HSET+1,ICOMP),UNCRT(IZ,HSET+1,ICOMP);
						// 2 FORMAT (T11, '|',
						// 1 (A4, 1PE10.3, '+-', 0PF4.1, '%', ' |'),
						// A4, 1PE10.3, '+-', 0PF4.1, '%', '|');

					}
					if (FMT == 1) {
						s = EGS4.format("", 11)
								+ "|"
								+ EGS4.format(LABELS[ICOMP - 1], 4)
								+ EGS4.format(
										RESULTS[IZ - 1][HSET - 1][ICOMP - 1],
										10)
								+ "+-"
								+ EGS4.format1(
										UNCRT[IZ - 1][HSET - 1][ICOMP - 1], 4,
										true) + "%" + "|";

						EGS4.seqStr = s;
						if (EGS4.iprint > 0)
							eq.printSequence(EGS4.seqStr);

						// WRITE(IOUT, 1) LABELS(ICOMP),RESULTS(IZ,HSET,ICOMP),
						// UNCRT(IZ,HSET,ICOMP);
						// 1 FORMAT (T11, '|',
						// A4, 1PE10.3, '+-', 0PF4.1, '%', '|');

					}
				}// ] "end IF (SFIG)"
				else {
					if (FMT == 3) {
						s = EGS4.format("", 10)
								+ "|"
								+ EGS4.format(LABELS[ICOMP - 1], 4)
								+ EGS4.format(
										RESULTS[IZ - 1][HSET - 1][ICOMP - 1],
										10)
								+ "+-"
								+ EGS4.format1(
										UNCRT[IZ - 1][HSET - 1][ICOMP - 1], 5,
										true)
								+ "%"
								+ "|"
								+ EGS4.format(LABELS[ICOMP - 1], 4)
								+ EGS4.format(RESULTS[IZ - 1][HSET][ICOMP - 1],
										10)
								+ "+-"
								+ EGS4.format1(UNCRT[IZ - 1][HSET][ICOMP - 1],
										5, true)
								+ "%"
								+ "|"
								+ EGS4.format(LABELS[ICOMP - 1], 4)
								+ EGS4.format(
										RESULTS[IZ - 1][HSET + 1][ICOMP - 1],
										10)
								+ "+-"
								+ EGS4.format1(
										UNCRT[IZ - 1][HSET + 1][ICOMP - 1], 5,
										true) + "%" + "|";

						EGS4.seqStr = s;
						if (EGS4.iprint > 0)
							eq.printSequence(EGS4.seqStr);

						// WRITE(IOUT, 6) LABELS(ICOMP),RESULTS(IZ,HSET,ICOMP),
						// UNCRT(IZ,HSET,ICOMP),LABELS(ICOMP),
						// RESULTS(IZ,HSET+1,ICOMP),UNCRT(IZ,HSET+1,ICOMP),
						// LABELS(ICOMP),RESULTS(IZ,HSET+2,ICOMP),
						// UNCRT(IZ,HSET+2,ICOMP);
						// 6 FORMAT (T10, '|',
						// 2 (A4, 1PE10.3, '+-', 0PF5.2, '%', '|'),
						// A4, 1PE10.3, '+-', 0PF5.2, '%', '|');

					}// ]
					if (FMT == 2) {
						s = EGS4.format("", 11)
								+ "|"
								+ EGS4.format(LABELS[ICOMP - 1], 4)
								+ EGS4.format(
										RESULTS[IZ - 1][HSET - 1][ICOMP - 1],
										10)
								+ "+-"
								+ EGS4.format1(
										UNCRT[IZ - 1][HSET - 1][ICOMP - 1], 5,
										true)
								+ "%"
								+ " |"
								+ EGS4.format(LABELS[ICOMP - 1], 4)
								+ EGS4.format(RESULTS[IZ - 1][HSET][ICOMP - 1],
										10)
								+ "+-"
								+ EGS4.format1(UNCRT[IZ - 1][HSET][ICOMP - 1],
										5, true) + "%" + "|";

						EGS4.seqStr = s;
						if (EGS4.iprint > 0)
							eq.printSequence(EGS4.seqStr);

						// WRITE(IOUT, 5) LABELS(ICOMP),RESULTS(IZ,HSET,ICOMP),
						// UNCRT(IZ,HSET,ICOMP),LABELS(ICOMP),
						// RESULTS(IZ,HSET+1,ICOMP),UNCRT(IZ,HSET+1,ICOMP);
						// 5 FORMAT (T11, '|',
						// 1 (A4, 1PE10.3, '+-', 0PF5.2, '%', ' |'),
						// A4, 1PE10.3, '+-', 0PF5.2, '%', '|');

					}
					if (FMT == 1) {
						s = EGS4.format("", 11)
								+ "|"
								+ EGS4.format(LABELS[ICOMP - 1], 4)
								+ EGS4.format(
										RESULTS[IZ - 1][HSET - 1][ICOMP - 1],
										10)
								+ "+-"
								+ EGS4.format1(
										UNCRT[IZ - 1][HSET - 1][ICOMP - 1], 5,
										true) + "%" + "|";

						EGS4.seqStr = s;
						if (EGS4.iprint > 0)
							eq.printSequence(EGS4.seqStr);

						// WRITE(IOUT, 4) LABELS(ICOMP),RESULTS(IZ,HSET,ICOMP),
						// UNCRT(IZ,HSET,ICOMP);
						// 4 FORMAT (T11, '|',
						// A4, 1PE10.3, '+-', 0PF5.2, '%', '|');

					}
				}// ] "end ELSE (SFIG)"
				PGTHROW = PGTHROW + 1;
			}// ] "end DO ICOMP"
				// "This piece of code surveys PGTHROW, and throws the page when"
				// "PGTHROW is near 65 lines."
				// if (MOD(PGTHROW,65)>(61-NCOMP)) [if ((PGTHROW %
				// 65)>(60-NCOMP))
			if ((PGTHROW % 65) > (61 - NCOMP)) {
				// "This is for the last vertical bin"
				if (IZ != DEEPNUM) {// "only for grid cut by PGTHROW"
					if (FMT == 1) {
						// WRITE (IOUT, 10) DEPTH_BINS(IZ+1);
						hs = "";
						for (int i = 1; i <= 23; i++) {
							hs = hs + "-";
						}
						s = EGS4.format("", 1)
								+ EGS4.format(DEPTH_BINS[IZ], 9, true);
						ll = s.length();
						ll = 11 - ll;
						if (ll > 0) {
							s = s + EGS4.format("", ll);
						}
						s = s + hs;// EGS4.format("-",23);
						EGS4.seqStr = s;
						if (EGS4.iprint > 0)
							eq.printSequence(EGS4.seqStr);

					}
					if (FMT == 2) {
						hs = "";
						for (int i = 1; i <= 46; i++) {
							hs = hs + "-";
						}
						s = EGS4.format("", 1)
								+ EGS4.format(DEPTH_BINS[IZ], 9, true);
						ll = s.length();
						ll = 11 - ll;
						if (ll > 0) {
							s = s + EGS4.format("", ll);
						}
						s = s + hs;// EGS4.format("-",45);
						EGS4.seqStr = s;
						if (EGS4.iprint > 0)
							eq.printSequence(EGS4.seqStr);

						// WRITE (IOUT, 11) DEPTH_BINS(IZ+1);
					}
					if (FMT == 3) {
						hs = "";
						for (int i = 1; i <= 69; i++) {
							hs = hs + "-";
						}
						s = EGS4.format("", 1)
								+ EGS4.format(DEPTH_BINS[IZ], 9, true);
						ll = s.length();
						ll = 11 - ll;
						if (ll > 0) {
							s = s + EGS4.format("", ll);
						}
						s = s + hs;// EGS4.format("-",67);
						EGS4.seqStr = s;
						if (EGS4.iprint > 0)
							eq.printSequence(EGS4.seqStr);

						// WRITE (IOUT, 12) DEPTH_BINS(IZ+1);
					}
					// "Start a new page"
					// WRITE (IOUT, *) '\f'; "page break"
					PGTHROW = 10;
					// WRITE(IOUT, 400) ' '; call egs_fdate(iout);
					// write(iout,'(//)');
					if (ROT) {
						s = EGS4.format("", 19);
						EGS4.seqStr = s + "ZONAL OUTPUT GRID: ROTATED";
						if (EGS4.iprint > 0)
							eq.printSequence(EGS4.seqStr);
						s = EGS4.format("", 20);
						EGS4.seqStr = s + "**************************";
						if (EGS4.iprint > 0)
							eq.printSequence(EGS4.seqStr);

						// WRITE(IOUT, 93) TITLE;
					} else {
						s = EGS4.format("", 19);
						EGS4.seqStr = s + "ZONAL OUTPUT GRID: NON-ROTATED";
						if (EGS4.iprint > 0)
							eq.printSequence(EGS4.seqStr);
						s = EGS4.format("", 20);
						EGS4.seqStr = s + "******************************";
						if (EGS4.iprint > 0)
							eq.printSequence(EGS4.seqStr);

						// WRITE(IOUT, 94) TITLE;
					}

					// WRITE (IOUT, 91) (RADIAL_BINS(IX), IX=HSET,HSET+FMT);]
					for (int IX = HSET; IX <= HSET + FMT; IX++) {
						if (IX == HSET) {
							s = EGS4.format("", 7);
						} else if (IX == HSET + 1) {
							ll = 30 - ll;
							if (ll > 0) {
								s = EGS4.format("", ll);
							} else {
								s = EGS4.format("", 30);
							}
						} else if (IX == HSET + 2) {
							ll = 53 - ll;
							if (ll > 0) {
								s = EGS4.format("", ll);
							} else {
								s = EGS4.format("", 53);
							}
						} else if (IX == HSET + 3) {
							ll = 70 - ll;
							if (ll > 0) {
								s = EGS4.format("", ll);
							} else {
								s = EGS4.format("", 70);
							}
						}

						// ll=s.length();ll=60-ll;s=s+EGS4.format("",ll);
						// EGS4.seqStr=s+EGS4.format(RADIAL_BINS[IX-1],9,true);
						// if(EGS4.iprint>1)
						// eq.printSequence(EGS4.seqStr);
						if (IX == HSET) {
							EGS4.seqStr = s
									+ EGS4.format(RADIAL_BINS[IX - 1], 9, true);
							ll = EGS4.seqStr.length();
						} else {
							EGS4.seqStr = EGS4.seqStr + s
									+ EGS4.format(RADIAL_BINS[IX - 1], 9, true);
							ll = EGS4.seqStr.length();
						}
					}
					// EGS4.seqStr=s+EGS4.format(RADIAL_BINS[IX-1],9,true);
					if (EGS4.iprint > 0)
						eq.printSequence(EGS4.seqStr);// auto salt at new line
				}// see wrtite
					// "The following statement delays the page throw until the "
					// "last line of the grid can be printed"
				else {
					DLYPT = HSET;
				}
			}// "end IF mod(pgthrow)",
		}// ] "end DO IZ"
			// "This is for the last vertical bin"
		if (FMT == 1) {
			hs = "";
			for (int i = 1; i <= 23; i++) {
				hs = hs + "-";
			}
			s = EGS4.format("", 1) + EGS4.format(DEPTH_BINS[DEEPNUM], 9, true);
			ll = s.length();
			ll = 11 - ll;
			if (ll > 0) {
				s = s + EGS4.format("", ll);
			}
			s = s + hs;// EGS4.format("-",23);
			EGS4.seqStr = s;
			if (EGS4.iprint > 0)
				eq.printSequence(EGS4.seqStr);

			// WRITE (IOUT, 10) DEPTH_BINS(DEEPNUM+1);
		}
		if (FMT == 2) {
			hs = "";
			for (int i = 1; i <= 46; i++) {
				hs = hs + "-";
			}
			s = EGS4.format("", 1) + EGS4.format(DEPTH_BINS[DEEPNUM], 9, true);
			ll = s.length();
			ll = 11 - ll;
			if (ll > 0) {
				s = s + EGS4.format("", ll);
			}
			s = s + hs;// EGS4.format("-",45);
			EGS4.seqStr = s;
			if (EGS4.iprint > 0)
				eq.printSequence(EGS4.seqStr);

			// WRITE (IOUT, 11) DEPTH_BINS(DEEPNUM+1);
			// 11 FORMAT (1X, F9.4, T11, 2 (23 ('-')));
		}
		if (FMT == 3) {
			hs = "";
			for (int i = 1; i <= 69; i++) {
				hs = hs + "-";
			}
			s = EGS4.format("", 1) + EGS4.format(DEPTH_BINS[DEEPNUM], 9, true);
			ll = s.length();
			ll = 11 - ll;
			if (ll > 0) {
				s = s + EGS4.format("", ll);
			}
			s = s + hs;// EGS4.format("-",67);
			EGS4.seqStr = s;
			if (EGS4.iprint > 0)
				eq.printSequence(EGS4.seqStr);

			// WRITE (IOUT, 12) DEPTH_BINS(DEEPNUM+1);
			// 12 FORMAT (1X, F9.4, T11, 3 (23 ('-')));
		}

		// "For the delayed page throw"
		if ((DLYPT == HSET) && (RADNUM - HSET > 2)) {
			// "Start a new page"
			// WRITE (IOUT, *) '\f'; "page break"
			PGTHROW = 10;
			// WRITE(IOUT, 400) ' '; call egs_fdate(iout); write(iout,'(//)');
			if (ROT) {
				s = EGS4.format("", 19);
				EGS4.seqStr = s + "ZONAL OUTPUT GRID: ROTATED";
				if (EGS4.iprint > 0)
					eq.printSequence(EGS4.seqStr);
				s = EGS4.format("", 20);
				EGS4.seqStr = s + "**************************";
				if (EGS4.iprint > 0)
					eq.printSequence(EGS4.seqStr);

				// WRITE(IOUT, 93) TITLE;
			} else {
				s = EGS4.format("", 19);
				EGS4.seqStr = s + "ZONAL OUTPUT GRID: NON-ROTATED";
				if (EGS4.iprint > 0)
					eq.printSequence(EGS4.seqStr);
				s = EGS4.format("", 20);
				EGS4.seqStr = s + "******************************";
				if (EGS4.iprint > 0)
					eq.printSequence(EGS4.seqStr);

				// WRITE(IOUT, 94) TITLE;
			}
		}
		PGTHROW = PGTHROW + 1;

	}

	/**
	 * Reset global variables for re-use.
	 */
	public static void reset() {
		RADNUM = 0;
		HSET = 0;
		FMT = 0;
		$MAXRZ = EGS4Geom.$MAXZREG;// =$MAXZREG} "MAX(MAXRADII,MAXZREG)"
		RADIAL_BINS = new double[$MAXRZ];
		DEPTH_BINS = new double[$MAXRZ];

		PGTHROW = 0;
		DEEPNUM = 0;
		ROT = false;
		DLYPT = 0;
		REGNUM = 0;
		MNUM1 = 0;
		MNUM2 = 0;
		MNUM3 = 0;
		MED_NAME1 = new String[11];
		MED_NAME2 = new String[11];
		MED_NAME3 = new String[11];
		NCOMP = 0;
		$MAXCMPTS = 14;// MUSAI AICI E MAXIMUM!!!FLUR are max de 12
		RESULTS = new double[EGS4Geom.$MAXZREG][EGS4Geom.$MAXRADII][$MAXCMPTS];
		UNCRT = new double[EGS4Geom.$MAXZREG][EGS4Geom.$MAXRADII][$MAXCMPTS];// ($MAXZREG,
																				// $MAXRADII,
																				// $MAXCMPTS),
		LABELS = new String[$MAXCMPTS];
		VALUES = new double[EGS4Geom.$MAXZREG][EGS4Geom.$MAXRADII][$MAXCMPTS];
		// next could be locals!!
		//ESTEPSON = false;
		ECUTON = false;
		PCUTON = false;
		TMP2 = new double[EGS4Geom.$MAXZREG][EGS4Geom.$MAXRADII][$MAXCMPTS];
		TMP3 = new double[$MAXRZ];
		CHINDEX = new double[$MAXRZ];
		TMP1 = 0;
		EXPLANATIONS = new String[$MAXCMPTS];
		IRL = 0;
		SFIG = false;
		COUNT = 0;
		//zoneB = false;
	}
}
