package danfulea.phys.egs;

/**
 * 
 * Utility class for handling basic geometry (cylindrical RZ or spherical). 
 * Based on EGS4/EGSnrc routines developed by SLAC(Stanford Linear Accelerator Center)/NRC (National Research Council of Canada).
 * @author Dan Fulea, 07 NOV. 2005
 */
public class EGS4Geom {
	// ;Copyright NRC;
	// ;"******************************************************************************
	// "
	// " ******************
	// " * *
	// " * geomrz.mortran *
	// " * *
	// " ******************
	// "
	// "
	// " This subroutine is used to do all the I/O associated with having
	// " a cylindrical geometry input.
	// "
	// " VERSION 1 A. Merovitz 03/98 for DWOR

	// "*******************************************************************************
	// " CYLINDRICAL GEOMETRY INPUT
	// " **************************
	// "*******************************************************************************
	// "
	// " GEOMRZ DELIMETERS: :start geometrical inputs:
	// " :stop geometrical inputs:
	// "
	// " METHOD OF INPUT
	// " = Groups (0) input groups of slabs of equal thickness
	// " = Individual (1) detailed input of the geometry and media.
	// " [ITERSE]
	// "
	// "-------------------------------------------------------------------------------
	// "
	// " Information defining depth boundaries along z axis (all dimensions cm)
	// "
	// " Only if METHOD OF INPUT= Groups
	// "
	// " Z OF FRONT FACE (R) start of first slab (real)
	// " NSLAB (M) # planar slabs in a group (integers)
	// " SLAB THICKNESS (M) thickness of each slab in the group (reals)
	// "
	// "
	// " Only if METHOD OF INPUT= Individual
	// "
	// " Z OF FRONT FACE (R) start of first slab (real)
	// " DEPTH BOUNDARIES (M) geometrical z-plane coordinates (reals)
	// "
	// "-------------------------------------------------------------------------------
	// "
	// " Information defining radial boundaries
	// "
	// " RADII (M) radii of cylinders defining the geometry (reals)
	// "
	// "-------------------------------------------------------------------------------
	// " MATERIAL INPUT
	// " **************
	// "
	// " MEDIA (M) material name which must match that in the
	// " pegs4 data set EXACTLY, including case.
	// " 24 characters max per medium, ended by , or ;
	// "
	// " Define which media in which regions, numbering in order given above.
	// "
	// " DESCRIPTION BY= Regions use the individual geometric region
	// " numbers
	// " = Planes use the IX, IZ values
	// " = Regions + Density same as Regions but specify medium
	// " density as well
	// " = Planes + Density same as Planes but specify medium
	// " density as well
	// " [DESCRIBE]
	// "
	// " If DESCRIPTION BY= Regions
	// "
	// " MEDNUM (M) the material number (integers)
	// " (MEDNUM=0 => vacuum)
	// " RHOR (only if DESCRIPTION= Regions + Density)
	// " (M) material density if different from the default
	// " (real) (if 0 then assumed to be default)
	// " START REGION (M) initial geometrical zone(irl) (integers) for
	// " this medium [NREGLO]
	// " STOP REGION (M) final geometrical zone(irl) (integers) for
	// " this medium.[NREGHI]
	// " ( >NREGLO to input more than one zone)
	// " DEFAULTS: MEDNUM=0 FOR REGION=1 (i.e. VACUUM)
	// " MEDNUM=1 FOR REGION=2,NREG
	// "
	// " These inputs should be thought of as triplets
	// " (quadruplets if RHOR is also specified) of
	// " MEDNUM, (RHOR,) START and STOP REGIONs which are used
	// " to specify the medium numbers for all regions where
	// " the medium is not the default (medium 1).
	// "
	// " If DESCRIPTION BY= Planes
	// "
	// " MEDNUM (M) the material number (integers)
	// " (MEDNUM=0 => vacuum)
	// " RHOR (only if DESCRIPTION= Planes + Density)
	// " (M) material density if different from the default
	// " (real) (if 0 then assumed to be default)
	// " START ZSLAB (M) initial zslab(iz) (integers)
	// " STOP ZSLAB (M) final zslab(iz) (integers)
	// " START RING (M) initial radial ring (ix) (integers)
	// " STOP RING (M) final radial ring (ix) (integers)
	// " DEFAULTS: MEDNUM=0 FOR REGION=1 (i.e. VACUUM)
	// " MEDNUM=1 FOR REGION=2,NREG
	// " These inputs shpuld be thought of as quintuples
	// " (sextuples if RHOR is specified) of numbers
	// " which specify the medium numbers and density by
	// " planar - radial regions
	// "
	// " One must use one type of input or the other so you must decide
	// " which is more convenient for any given case.
	// "
	// "***************************************************************************"
	public static EgsQuestion eq;

	public static int NR = 0;
	public static int NZ = 0;
	public static int NC = 0;
	public static int $MAXRADII = 8;// "MAX # OF DOSE SCORING RADIAL ZONES           "
	public static int $MAXZREG = 20;// "MAX # OF DOSE SCORING PLANAR ZONES           "
	public static int $MAXZPLANE = $MAXZREG + 1;
	public static double[] ZPLANE = new double[$MAXZPLANE];
	public static double[] RCYL = new double[$MAXRADII + 1];// (0:$MAXRADII)
	public static double[] RSPH = new double[$MAXRADII + 1];// (0:$MAXRADII),
	public static double[] CYRAD2 = new double[$MAXRADII + 1];// 1 biased
																// ($MAXRADII+1)
	public static int[] ntrack = new int[EGS4.$MXREG];// cav
	public static int nreg = 0;
	public static int NPLANE = 0;
	public static int IDNEAR = 0;

	public static int $MAXCANGLE = 48;// default==="MAX # OF CONICAL BOUNDARIES (includes 180o )"

	public static int[] NRADIUS = new int[$MAXRADII + 1];// (0:$MAXRADII)
	public static int[] NCONES = new int[$MAXCANGLE + 1];// (0:$MAXCANGLE)
	public static double[] ANGLES = new double[$MAXCANGLE + 1];// (0:$MAXCANGLE)
	public static int nCONES = 0;
	public static int nANGLES = 0;
	public static int nRADIUS = 0;
	public static int nRADII = 0;
	public static double[] RADII = new double[$MAXRADII + 1];// (0:$MAXRADII)
	public static final double RADII_MIN = 0.0000001;
	public static final double RADII_MAX = 999999.9;
	public static final double RADII_DEFAULT = 1.0;
	public static final int nCONES_MAX = 1000000;
	public static final int nCONES_MIN = 0;
	public static final int nCONES_DEFAULT = 1;
	public static final double nANGLES_MAX = 180;
	public static final double nANGLES_MIN = 0.0;
	public static final double nANGLES_DEFAULT = 0.0;
	public static int nRADIUS_MAX = $MAXRADII;
	public static final int nRADIUS_MIN = 0;
	public static int nRADIUS_DEFAULT = $MAXRADII;
	public static int numcavreg = 0;
	public static int[] cavreg = new int[EGS4.$MXREG];
	public static final int cavreg_MIN = 0;
	public static final int cavreg_MAX = 999999;
	public static final int cavreg_DEFAULT = 1;

	public static double[] cosalp = new double[$MAXCANGLE + 1];// (0:$MAXCANGLE)
	public static double[] TANALP = new double[$MAXCANGLE + 1];// (0:$MAXCANGLE),
	public static double[] TANAL2 = new double[$MAXCANGLE + 1];// (0:$MAXCANGLE),
	public static double[] SINALP = new double[$MAXCANGLE + 1];// (0:$MAXCANGLE),
	public static double[] ALPHA = new double[$MAXCANGLE + 1];// (0:$MAXCANGLE),
	public static double[] RSPH2 = new double[$MAXRADII + 1];// 1
																// biased->($MAXRADII+1),
	public static int NPLAN1 = 0;
	public static int NPLAN2 = 0;
	// "ZPLANE(IZ) CONTAINS THE REAL VALUED COORDINATE OF THE IZ'TH PLANE
	// "RCYL(IX) CONTAINS THE REAL VALUED COORDINATE OF THE IX'TH CYL.
	// "CYRAD2(IX) =RCYL(IX)**2
	// "NTRACK(IRL) =1 IF CAVITY REGION ELSE =0
	// "NZ/NR NUMBER OF PLANAR SLABS/CYLINDRICAL REGIONS DEFINING THE TARGET
	// "NREG =NZ*NR+1 (+1 FOR THE SURROUNDING VACUUM), NPLANE=NZ+1
	// "NPLANE NUMBER OF PLANES DEFINING THE GEOMETRY
	public static int NSUMCV = 0;// = NO. OF REGIONS MAKING UP THE CAVITY
	// THE ARRAY OF ZONES COMPRISING THE CAVITY REGION:
	public static int[] ISUMCV = new int[EGS4.$MXREG];
	public static int[] NREGLO = new int[EGS4.$MXMED];// =0;
	public static int[] NREGHI = new int[EGS4.$MXMED];// =0;// "material input"
	public static int[] NZHI = new int[EGS4.$MXMED];// =0;=0;
	public static int[] NZLO = new int[EGS4.$MXMED];// =0;=0;
	public static int[] NRHI = new int[EGS4.$MXMED];// =0;=0;
	public static int[] NRLO = new int[EGS4.$MXMED];// =0;=0; //
													// "material input"
	public static int ITERSE = 0;
	public static int $NVALUE = 100;// "max number of values per input"
	public static int[] NSLAB = new int[$NVALUE];// ($NVALUE),
	public static int nNSLAB = 1;
	public static int nCyl = 1;
	public static int DESCRIBE = 0;// "cylindrical geom. input"
	public static double[] DELTAZ = new double[$NVALUE];// ($NVALUE),
	public static double RHORI = 0.0;// "cylindrical geom. input"

	public static int iterseindex = 0;
	public static final int iGROUPS = 0;
	public static final int iINDIVIDUAL = 1;
	public static final int iCAVITY_INFORMATION = 2;

	public static double WALLTH = 0.0;
	public static double CAVRAD = 0.0;
	public static double CAVLNG = 0.0;
	public static double ELERAD = 0.0;
	public static String SLENGHT = "";
	public static int medset = 0;
	public static String airs = "";
	public static String electrods = "";
	public static int nMEDIA = 0;
	public static final double WALLTH_MIN = 1.e-10;
	public static final double WALLTH_MAX = 999999;
	public static final double WALLTH_DEFAULT = 0.273;
	public static final double CAVRAD_MIN = 1.e-10;
	public static final double CAVRAD_MAX = 999999;
	public static final double CAVRAD_DEFAULT = 1.0;
	public static final double CAVLNG_MIN = 1.e-10;
	public static final double CAVLNG_MAX = 999999;
	public static final double CAVLNG_DEFAULT = 0.2;
	public static final double ELERAD_MIN = 0.0;
	public static final double ELERAD_MAX = 999999;
	public static final double ELERAD_DEFAULT = 0.0;

	public static double Z_OF_FRONT_FACE = 0.0;
	public static final double Z_OF_FRONT_FACE_default = 0.0;
	public static final double Z_OF_FRONT_FACE_min = -999999.0;
	public static final double Z_OF_FRONT_FACE_max = 999999.0;
	public static final double ZPLANE_MIN = -999999.0;
	public static final double ZPLANE_MAX = 999999.0;
	public static final double ZPLANE_DEFAULT = 1.0;
	public static final int NSLAB_MIN = 0;
	public static final int NSLAB_MAX = 999999;
	public static final int NSLAB_DEFAULT = 1;
	public static final double DELTAZ_MIN = 0.0;
	public static final double DELTAZ_MAX = 999999.0;
	public static final double DELTAZ_DEFAULT = 1.0;
	public static final double RCYL_MIN = 0.0000001;
	public static final double RCYL_MAX = 999999.9;
	public static final double RCYL_DEFAULT = 1000.0;
	public static final int DESCRIBE_REGIONS = 0;
	public static final int DESCRIBE_PLANES = 1;
	public static final int DESCRIBE_REGIONS_DENSITY = 2;
	public static final int DESCRIBE_PLANES_DENSITY = 3;
	public static int nMEDNUM = 1;
	public static final int MEDNUM_DEFAULT = 1;
	public static final int MEDNUM_MIN = 0;
	public static int MEDNUM_MAX = EGS4.NMED;
	public static final double RHOR_MIN = 0.0;
	public static final double RHOR_MAX = 999999.0;
	public static final double RHOR_DEFAULT = 0.0;
	public static int nNREGLO = 0;
	public static final int NREGLO_MIN = 0;
	public static final int NREGLO_DEFAULT = 1;
	public static int NREGLO_MAX = 0;
	public static int nNREGHI = 0;
	public static final int NREGHI_MIN = 0;
	public static final int NREGHI_DEFAULT = 1;
	public static int NREGHI_MAX = 0;
	public static int nRHOR = 0;
	public static double[] RHOR = new double[EGS4.$MXMED];// =0;
	public static int[] MEDNUM = new int[EGS4.$MXMED];
	public static int nNZLO = 0;
	public static int nNZHI = 0;
	public static int nNRLO = 0;
	public static int nNRHI = 0;
	public static final int NZLO_MIN = 1;
	public static final int NZLO_DEFAULT = 1;
	public static int NZLO_MAX = 0;
	public static final int NZHI_MIN = 1;
	public static final int NZHI_DEFAULT = 1;
	public static int NZHI_MAX = 0;
	public static final int NRLO_MIN = 1;
	public static final int NRLO_DEFAULT = 1;
	public static int NRLO_MAX = 0;
	public static final int NRHI_MIN = 1;
	public static final int NRHI_DEFAULT = 1;
	public static int NRHI_MAX = 0;
	public static final int NSUMCV_MIN = 0;
	public static final int NSUMCV_DEFAULT = 1;
	public static int NSUMCV_MAX = 0;
	public static final int ISUMCV_MIN = 2;
	public static int ISUMCV_DEFAULT = 0;
	public static int ISUMCV_MAX = 0;

	/**
	 * Reset all global variables for re-use.
	 */
	public static void reset() {
		NR = 0;
		NZ = 0;
		NC = 0;// $MAXRADII=8;$MAXZREG=20;
		$MAXZPLANE = $MAXZREG + 1;
		ZPLANE = new double[$MAXZPLANE];
		RCYL = new double[$MAXRADII + 1];
		RSPH = new double[$MAXRADII + 1];
		CYRAD2 = new double[$MAXRADII + 1];
		ntrack = new int[EGS4.$MXREG];
		nreg = 0;
		NPLANE = 0;
		IDNEAR = 0;
		// $MAXCANGLE=48;
		cosalp = new double[$MAXCANGLE + 1];
		TANALP = new double[$MAXCANGLE + 1];
		TANAL2 = new double[$MAXCANGLE + 1];
		SINALP = new double[$MAXCANGLE + 1];
		ALPHA = new double[$MAXCANGLE + 1];
		RSPH2 = new double[$MAXRADII + 1];
		NPLAN1 = 0;
		NPLAN2 = 0;
		NZHI = new int[EGS4.$MXMED];
		NZLO = new int[EGS4.$MXMED];
		NRHI = new int[EGS4.$MXMED];
		NRLO = new int[EGS4.$MXMED];
		ITERSE = 0;// $NVALUE=100;
		NSLAB = new int[$NVALUE];
		DESCRIBE = 0;
		DELTAZ = new double[$NVALUE];
		RHORI = 0.0;
		iterseindex = 0;

		Z_OF_FRONT_FACE = 0.0;
		WALLTH = 0.0;
		CAVRAD = 0.0;
		CAVLNG = 0.0;
		ELERAD = 0.0;
		SLENGHT = "";
		medset = 0;
		airs = "";
		electrods = "";
		nMEDIA = 1;
		nNSLAB = 1;
		nCyl = 1;
		nMEDNUM = 1;
		MEDNUM_MAX = EGS4.NMED;
		nRHOR = 0;
		nNREGLO = 0;
		nNREGHI = 0;
		MEDNUM = new int[EGS4.$MXMED];
		RHOR = new double[EGS4.$MXMED];
		NREGLO = new int[EGS4.$MXMED];
		NREGHI = new int[EGS4.$MXMED];
		NREGLO_MAX = 0;
		NREGHI_MAX = 0;
		nNZLO = 0;
		nNZHI = 0;
		nNRLO = 0;
		nNRHI = 0;
		NZLO_MAX = 0;
		NZHI_MAX = 0;
		NRLO_MAX = 0;
		NRHI_MAX = 0;

		NSUMCV_MAX = 0;
		ISUMCV_DEFAULT = 0;
		ISUMCV_MAX = 0;
		NRADIUS = new int[$MAXRADII + 1];
		NCONES = new int[$MAXCANGLE + 1];
		nCONES = 0;
		nANGLES = 0;
		ANGLES = new double[$MAXCANGLE + 1];
		nRADIUS_MAX = $MAXRADII;
		nRADIUS_DEFAULT = $MAXRADII;
		nRADIUS = 0;
		nRADII = 0;
		RADII = new double[$MAXRADII + 1];
		cavreg = new int[EGS4.$MXREG];
		numcavreg = 0;
	}

	/**
	 * Build cylindrical RZ geometry.
	 */
	public static void GEOMRZ() {
		// iterse must be known.
		ITERSE = iterseindex;
		if (ITERSE == iCAVITY_INFORMATION) {
			// WALL THICKNESS
			if (WALLTH < WALLTH_MIN || WALLTH > WALLTH_MAX) {
				WALLTH = WALLTH_DEFAULT;
			}
			// CAVITY RADIUS
			if (CAVRAD < CAVRAD_MIN || CAVRAD > CAVRAD_MAX) {
				CAVRAD = CAVRAD_DEFAULT;
			}
			// CAVITY LENGTH
			if (CAVLNG < CAVLNG_MIN || CAVLNG > CAVLNG_MAX) {
				CAVLNG = CAVLNG_DEFAULT;
			}
			// ELECTRODE RADIUS
			if (ELERAD < ELERAD_MIN || ELERAD > ELERAD_MAX) {
				ELERAD = ELERAD_DEFAULT;
			}
			// WALL MATERIAL-->.pegs4dat!!!!!!!!!!!
			EGS4.MEDIA[0] = SLENGHT;
			medset++;// score the actual set!!
			if (ELERAD == 0.0) {// NO ELECTRODE"
				EGS4.seqStr = " WALL THICKNESS "
						+ EGS4.format(WALLTH, 10, true) + " cms "
						+ " WALL MATERIAL: " + SLENGHT + " GAS CAVITY RADIUS "
						+ EGS4.format(CAVRAD, 10, true) + " cms "
						+ " GAS CAVITY LENGTH " + EGS4.format(CAVLNG, 10, true)
						+ " cms ";
				if (EGS4.iprint > 0)// WATCH is next->allways show
					eq.printSequence(EGS4.seqStr);

				// "SECOND MATERIAL IS GAS (AIR)"
				// MEDIA(1,2)='A';MEDIA(2,2)='I';MEDIA(3,2)='R';
				// DO J=4,24 [MEDIA(J,2)=' ';]
				EGS4.MEDIA[medset] = airs;
				medset++;// 1

				// "FILL THE GEOMETRY WITH THE MATERIALS"
				EGS4.NMED = 2; // "TWO MEDIA"
				EGS4.MED[0] = 0; // "VACUUM ENVELOPE"
				EGS4.MED[2] = 2; // "THE GAS CAVITY"
				EGS4.MED[1] = 1;
				for (int J = 4; J <= 6; J++) {
					EGS4.MED[J - 1] = 1;
				} // "THE WALLS"

				// "DEFINE THE RADII OF THE GEOMETRY"
				NR = 2; // "TWO RADII"
				RCYL[1] = CAVRAD; // "FIRST IS THE CAVITY RADIUS" 0 BIASED
				RCYL[2] = CAVRAD + WALLTH; // "SECOND IS THE CHAMBER RADIUS"
				for (int J = 1; J <= 2; J++) {
					CYRAD2[J - 1] = RCYL[J] * RCYL[J];
				}
				;

				// "DEFINE THE PLANES OF THE GEOMETRY"
				NZ = 3;
				NPLANE = 4; // "THREE PLANAR ZONES, FOUR DEFINING PLANES"
				ZPLANE[0] = 0.; // "FIRST PLANE AT 0"
				ZPLANE[1] = WALLTH;
				ZPLANE[2] = WALLTH + CAVLNG;
				ZPLANE[3] = 2. * WALLTH + CAVLNG;
				nreg = NR * NZ + 1; // "TOTAL NUMBER OF GEOMETRICAL REGIONS"

				// "DEFINE THE CAVITY REGION"
				NSUMCV = 1;
				ISUMCV[0] = 3;
				for (int I = 2; I <= nreg; I++) {
					ISUMCV[I - 1] = 0;
				}
			} else {// WITH AN ELECTRODE"
					// "SECOND MATERIAL IS GAS"
				EGS4.MEDIA[medset] = airs;
				medset++;// 1
				// "READ IN THE ELECTRODE MATERIAL"
				EGS4.MEDIA[medset] = electrods;
				medset++;// 2
				// DO J=SLENGHT+1, 24 [MEDIA(J,3)=' ';]

				EGS4.seqStr = " WALL THICKNESS "
						+ EGS4.format(WALLTH, 10, true) + " cms "
						+ " WALL MATERIAL: " + SLENGHT + " GAS CAVITY RADIUS "
						+ EGS4.format(CAVRAD, 10, true) + " cms "
						+ " GAS CAVITY LENGTH " + EGS4.format(CAVLNG, 10, true)
						+ " cms " + " ELECTRODE RADIUS"
						+ EGS4.format(ELERAD, 10, true) + " cms "
						+ " ELECTRODE MATERIAL: " + electrods;
				if (EGS4.iprint > 0)// WATCH is next->allways show
					eq.printSequence(EGS4.seqStr);

				// "FILL THE GEOMETRY WITH THE MATERIALS"
				EGS4.MED[0] = 0; // "VACUUM ENVELOPE"

				// "CHECK ALL THE CHARACTERS TO SEE IF THE WALL MATERIAL"
				// "AND THE ELECTRODE MATERIAL ARE THE SAME"->ISSAME = 0=>same;
				// 1 otherwise

				if (SLENGHT.compareTo(electrods) == 0) {// "WALL AND ELECTRODE ARE THE SAME MATERIAL"
					EGS4.NMED = 2;
					EGS4.MED[2] = 1;// MED(3)=1; //"THE ELECTRODE"
					EGS4.MED[5] = 2;// MED(6)=2; //"THE GAS CAVITY"
					// MED(5)=1; DO J=7,10 [MED(J)=1;] //"THE WALLS"
					EGS4.MED[4] = 1;
					for (int J = 7; J <= 10; J++) {
						EGS4.MED[J - 1] = 1;
					} // "THE WALLS"
				} else {// "WALL AND ELECTRODE HAVE DIFFERENT MATERIALS"
					EGS4.NMED = 3;
					EGS4.MED[2] = 3;// MED(3)=3; //"THE ELECTRODE"
					EGS4.MED[5] = 2;// MED(6)=2; //"THE GAS CAVITY"
					// MED(5)=1;DO J=7,10[MED(J)=1]; //"THE WALLS"
					EGS4.MED[4] = 1;
					for (int J = 7; J <= 10; J++) {
						EGS4.MED[J - 1] = 1;
					} // "THE WALLS"
				}

				// "DEFINE THE RADII OF THE GEOMETRY"
				NR = 3; // "THREE RADII"
				RCYL[1] = ELERAD; // "FIRST IS THE ELECTRODE RADIUS"
				RCYL[2] = CAVRAD; // "SECOND IS THE OUTER CAVITY RADIUS"
				RCYL[3] = CAVRAD + WALLTH; // "LAST IS THE CHAMBER RADIUS"
				for (int J = 1; J <= 3; J++) {
					CYRAD2[J - 1] = RCYL[J] * RCYL[J];
				}

				// "DEFINE THE PLANES OF THE GEOMETRY"
				NZ = 3;
				NPLANE = 4; // "THREE PLANAR ZONES, FOUR DEFINING PLANES"
				ZPLANE[0] = 0.; // "FIRST PLANE AT 0"
				ZPLANE[1] = WALLTH;
				ZPLANE[2] = WALLTH + CAVLNG;
				ZPLANE[3] = 2. * WALLTH + CAVLNG;
				nreg = NR * NZ + 1; // "TOTAL NUMBER OF GEOMETRICAL REGIONS"

				// "DEFINE THE CAVITY REGION"
				NSUMCV = 1;
				ISUMCV[0] = 6;
				for (int J = 2; J <= nreg; J++) {
					ISUMCV[J - 1] = 0;
				}
			}

			if (EGS4.NMED > EGS4.$MXMED) {
				EGS4.STOPPROGRAM = true;
				EGS4.seqStr = " ERROR: No. of media > max. no. of media allowed";
				eq.printSequence(EGS4.seqStr);
				EGS4.seqStr = " Increase $MXMED!!";
				eq.printSequence(EGS4.seqStr);
				return;
			}

		} else// ITERSE->0 sau 1
		{
			if (Z_OF_FRONT_FACE < Z_OF_FRONT_FACE_min
					|| Z_OF_FRONT_FACE > Z_OF_FRONT_FACE_max)
				Z_OF_FRONT_FACE = Z_OF_FRONT_FACE_default;
			ZPLANE[0] = Z_OF_FRONT_FACE;

			if (ITERSE == 0) {
				if (nNSLAB < 1) {
					EGS4.STOPPROGRAM = true;
					EGS4.seqStr = " ******** ERROR in geometrical inputs! **** ";
					eq.printSequence(EGS4.seqStr);
					return;
				}
				// NSLAB
				for (int i = 1; i <= nNSLAB; i++) {
					if (NSLAB[i - 1] < NSLAB_MIN || NSLAB[i - 1] > NSLAB_MAX)
						NSLAB[i - 1] = NSLAB_DEFAULT;
					// SLAB THICKNESS
					if (DELTAZ[i - 1] < DELTAZ_MIN
							|| DELTAZ[i - 1] > DELTAZ_MAX)
						DELTAZ[i - 1] = DELTAZ_DEFAULT;
				}

				int COUNT = 1;
				double ADDING = 0.0;
				for (int PLN = 1; PLN <= nNSLAB; PLN++) {
					if (PLN == 1) {
						ADDING = 0.0;
					}
					// else {ADDING=ADDING+DELTAZ[PLN-1]*NSLAB[PLN-1];}
					else {
						ADDING = ADDING + DELTAZ[PLN - 2] * NSLAB[PLN - 2];
					}

					for (int k = 1; k <= NSLAB[PLN - 1]; k++) {
						COUNT = COUNT + 1;
						ZPLANE[COUNT - 1] = ZPLANE[0] + ADDING
								+ DELTAZ[PLN - 1] * k;
					}

				}
				NZ = COUNT - 1;
				// ==================
				EGS4.seqStr = " NUMBER OF PLANAR SLABS: " + NZ;
				if (EGS4.iprint > 1)
					eq.printSequence(EGS4.seqStr);

				NPLANE = NZ + 1;

				EGS4.seqStr = " GEOMETRICAL Z-PLANE COORDINATES:";
				if (EGS4.iprint > 1)
					eq.printSequence(EGS4.seqStr);
				for (int I = 1; I <= NZ + 1; I++) {
					EGS4.seqStr = EGS4.format(ZPLANE[I - 1], 10, true);
					if (EGS4.iprint > 1)
						eq.printSequence(EGS4.seqStr);
				}

			}

			if (ITERSE == 1) {
				NZ = nNSLAB;
				if (NZ - 1 > $MAXZREG) {
					EGS4.STOPPROGRAM = true;
					EGS4.seqStr = " Number of depth regions is greater than $MAXZREG!!";
					eq.printSequence(EGS4.seqStr);
					return;
				}
				// DO I=1, NZ [ZPLANE(I+1)=VALUE(NUM_ZPLANES,I);]
				for (int i = 1; i <= NZ; i++) {
					// DEPTH BOUNDARIES
					if (ZPLANE[i] < ZPLANE_MIN || ZPLANE[i] > ZPLANE_MAX)
						ZPLANE[i] = ZPLANE_DEFAULT;
				}

				EGS4.seqStr = " NUMBER OF PLANAR SLABS: " + NZ;
				if (EGS4.iprint > 1)
					eq.printSequence(EGS4.seqStr);

				NPLANE = NZ + 1;

				EGS4.seqStr = " GEOMETRICAL Z-PLANE COORDINATES:";
				if (EGS4.iprint > 1)
					eq.printSequence(EGS4.seqStr);
				for (int I = 1; I <= NPLANE; I++) {
					EGS4.seqStr = EGS4.format(ZPLANE[I - 1], 10, true);
					if (EGS4.iprint > 1)
						eq.printSequence(EGS4.seqStr);
				}

			}

			NR = nCyl;
			if (NR > $MAXRADII) {
				EGS4.STOPPROGRAM = true;
				EGS4.seqStr = " Number of radial regions is greater than $MAXRADII!!";
				eq.printSequence(EGS4.seqStr);
				return;
			}
			nreg = NR * NZ + 1; // "define the number of regions"
			if (nreg > EGS4.$MXREG) {
				EGS4.STOPPROGRAM = true;
				EGS4.seqStr = " Total number of regions is greater than $MXREG!!";
				eq.printSequence(EGS4.seqStr);
				return;
			}

			EGS4.seqStr = " Number of concentric rings: " + NR;
			if (EGS4.iprint > 1)
				eq.printSequence(EGS4.seqStr);

			RCYL[0] = 0.0;
			for (int i = 1; i <= NR; i++) {
				if (RCYL[i] < RCYL_MIN || RCYL[i] > RCYL_MAX)
					RCYL[i] = RCYL_DEFAULT;

				CYRAD2[i - 1] = RCYL[i] * RCYL[i];

				EGS4.seqStr = " Ring radius: " + i + " "
						+ EGS4.format(RCYL[i], 10, true);
				if (EGS4.iprint > 1)
					eq.printSequence(EGS4.seqStr);
			}
			// MEDIA
			EGS4.NMED = nMEDIA;
			if (EGS4.NMED > EGS4.$MXMED) {
				EGS4.STOPPROGRAM = true;
				EGS4.seqStr = " ERROR: No. of media > max. no. of media allowed";
				eq.printSequence(EGS4.seqStr);
				EGS4.seqStr = " Increase $MXMED!!";
				eq.printSequence(EGS4.seqStr);
				return;
			}

			EGS4.seqStr = " Number of media: " + EGS4.NMED;
			if (EGS4.iprint > 1)
				eq.printSequence(EGS4.seqStr);

			for (int i = 1; i <= EGS4.NMED; i++) {
				EGS4.seqStr = " MEDIUM: " + i + "  " + EGS4.MEDIA[i - 1];
				if (EGS4.iprint > 1)
					eq.printSequence(EGS4.seqStr);
			}

			EGS4.MED[0] = 0;
			for (int I = 2; I <= nreg; I++) {
				EGS4.MED[I - 1] = 1;
			} // "defaults"
			EGS4.seqStr = " Number of geometrical zones = " + nreg
					+ " ,vacuum in first region";
			if (EGS4.iprint > 1)
				eq.printSequence(EGS4.seqStr);

			// DESCRIPTION BY
			if (DESCRIBE < DESCRIBE_REGIONS
					|| DESCRIBE > DESCRIBE_PLANES_DENSITY) {
				EGS4.STOPPROGRAM = true;
				EGS4.seqStr = " ERROR in describing regions!!!";
				eq.printSequence(EGS4.seqStr);
				return;
			}

			MEDNUM_MAX = EGS4.NMED;
			for (int i = 1; i <= nMEDNUM; i++) {
				if (MEDNUM[i - 1] < MEDNUM_MIN || MEDNUM[i - 1] > MEDNUM_MAX)
					MEDNUM[i - 1] = MEDNUM_DEFAULT;
			}

			if (DESCRIBE == DESCRIBE_REGIONS_DENSITY
					|| DESCRIBE == DESCRIBE_PLANES_DENSITY) {
				if (nRHOR != nMEDNUM) {
					EGS4.STOPPROGRAM = true;
					EGS4.seqStr = " MEDNUM AND RHOR MUST BOTH HAVE THE SAME NUMBER OF VALUES!";
					eq.printSequence(EGS4.seqStr);
					return;
				}

				for (int i = 1; i <= nMEDNUM; i++) {
					if (RHOR[i - 1] < RHOR_MIN || RHOR[i - 1] > RHOR_MAX)
						RHOR[i - 1] = RHOR_DEFAULT;
				}
			}

			NREGLO_MAX = nreg;
			NREGHI_MAX = nreg;
			if (DESCRIBE == DESCRIBE_REGIONS
					|| DESCRIBE == DESCRIBE_REGIONS_DENSITY) {
				if (nNREGLO != nMEDNUM || nNREGHI != nMEDNUM) {
					EGS4.STOPPROGRAM = true;
					EGS4.seqStr = " MEDNUM AND START AND STOP REGIONS MUST HAVE THE SAME NUMBER OF VALUES!";
					eq.printSequence(EGS4.seqStr);
					return;
				}
				for (int i = 1; i <= nMEDNUM; i++) {
					if (NREGLO[i - 1] < NREGLO_MIN
							|| NREGLO[i - 1] > NREGLO_MAX)
						NREGLO[i - 1] = NREGLO_DEFAULT;
					if (NREGHI[i - 1] < NREGHI_MIN
							|| NREGHI[i - 1] > NREGHI_MAX)
						NREGHI[i - 1] = NREGHI_DEFAULT;
				}
				// ---------------------------------------
				// "In the following, we allow for vacuum input"
				if (MEDNUM[0] >= 0)// (VALUE(NUM_MEDNUM,1) >= 0)
				{
					for (int i = 1; i <= nMEDNUM; i++) {
						RHORI = 0.0;
						if (DESCRIBE == DESCRIBE_REGIONS_DENSITY)
							RHORI = RHOR[i - 1];// VALUE(NUM_RHOR,I);
						// ->MEDNUM[i-1]
						if (NREGHI[i - 1] <= NREGLO[i - 1]) {
							if (NREGLO[i - 1] > 0
									&& NREGLO[i - 1] <= EGS4.$MXREG) {
								// MED(NREGLO)=MEDNUM;
								EGS4.MED[NREGLO[i - 1] - 1] = MEDNUM[i - 1];
								if (RHORI > 0.0)
									EGS4.RHOR[NREGLO[i - 1] - 1] = RHORI;
								if (DESCRIBE == DESCRIBE_REGIONS_DENSITY) {
									if (RHORI > 0.0) {
										EGS4.seqStr = " REGION: "
												+ NREGLO[i - 1]
												+ " = MATERIAL: "
												+ MEDNUM[i - 1] + " RHOR ="
												+ EGS4.format(RHORI, 8, true);
										if (EGS4.iprint > 1)
											eq.printSequence(EGS4.seqStr);
									} else {
										EGS4.seqStr = " REGION: "
												+ NREGLO[i - 1]
												+ " = MATERIAL: "
												+ MEDNUM[i - 1]
												+ " RHOR = DEFAULT";
										if (EGS4.iprint > 1)
											eq.printSequence(EGS4.seqStr);
									}
								} else {
									EGS4.seqStr = " REGION: " + NREGLO[i - 1]
											+ " = MATERIAL: " + MEDNUM[i - 1];
									if (EGS4.iprint > 1)
										eq.printSequence(EGS4.seqStr);
								}
							}

						} else {
							for (int K = NREGLO[i - 1]; K <= NREGHI[i - 1]; K++) {
								if (K > 0 && K <= EGS4.$MXREG) {
									EGS4.MED[K - 1] = MEDNUM[i - 1];
									if (RHORI > 0.0)
										EGS4.RHOR[K - 1] = RHORI;
								}
							}
							if (DESCRIBE == DESCRIBE_REGIONS_DENSITY) {
								if (RHORI > 0.0) {
									EGS4.seqStr = " REGION: " + NREGLO[i - 1]
											+ " TO REGION " + NREGHI[i - 1]
											+ " = MATERIAL: " + MEDNUM[i - 1]
											+ " RHOR ="
											+ EGS4.format(RHORI, 8, true);
									if (EGS4.iprint > 1)
										eq.printSequence(EGS4.seqStr);
								} else {
									EGS4.seqStr = " REGION: " + NREGLO[i - 1]
											+ " TO REGION " + NREGHI[i - 1]
											+ " = MATERIAL: " + MEDNUM[i - 1]
											+ " RHOR = DEFAULT";
									if (EGS4.iprint > 1)
										eq.printSequence(EGS4.seqStr);
								}
							} else {
								EGS4.seqStr = " REGION: " + NREGLO[i - 1]
										+ " TO REGION " + NREGHI[i - 1]
										+ " = MATERIAL: " + MEDNUM[i - 1];
								if (EGS4.iprint > 1)
									eq.printSequence(EGS4.seqStr);
							}

						}

					}// for
				}
			}

			NZLO_MAX = NZ;
			NZHI_MAX = NZ;
			NRLO_MAX = NR;
			NRHI_MAX = NR;
			if (DESCRIBE == DESCRIBE_PLANES
					|| DESCRIBE == DESCRIBE_PLANES_DENSITY) {
				if ((nNZLO != nMEDNUM || nNZHI != nMEDNUM)
						|| (nNRLO != nMEDNUM || nNRHI != nMEDNUM)) {
					EGS4.STOPPROGRAM = true;
					EGS4.seqStr = " MEDNUM AND START AND STOP REGIONS MUST HAVE THE SAME NUMBER OF VALUES!";
					eq.printSequence(EGS4.seqStr);
					return;
				}

				for (int i = 1; i <= nMEDNUM; i++) {
					if (NZLO[i - 1] < NZLO_MIN || NZLO[i - 1] > NZLO_MAX)
						NZLO[i - 1] = NZLO_DEFAULT;
					if (NZHI[i - 1] < NZHI_MIN || NZHI[i - 1] > NZHI_MAX)
						NZHI[i - 1] = NZHI_DEFAULT;
					if (NRLO[i - 1] < NRLO_MIN || NRLO[i - 1] > NRLO_MAX)
						NRLO[i - 1] = NRLO_DEFAULT;
					if (NRHI[i - 1] < NRHI_MIN || NRHI[i - 1] > NRHI_MAX)
						NRHI[i - 1] = NRHI_DEFAULT;
				}

				// ---------------------------------------
				// "In the following, we allow for vacuum input"
				if (MEDNUM[0] >= 0)// (VALUE(NUM_MEDNUM,1) >= 0)
				{
					for (int i = 1; i <= nMEDNUM; i++) {
						RHORI = 0.0;
						if (DESCRIBE == DESCRIBE_PLANES_DENSITY)
							RHORI = RHOR[i - 1];
						if (NZLO[i - 1] > NZHI[i - 1]) {
							NZHI[i - 1] = NZLO[i - 1];
						}
						if (NRLO[i - 1] > NRHI[i - 1]) {
							NRHI[i - 1] = NRLO[i - 1];
						}

						for (int IZ = NZLO[i - 1]; IZ <= NZHI[i - 1]; IZ++) {
							for (int IX = NRLO[i - 1]; IX <= NRHI[i - 1]; IX++) {
								int REGNUM = IZ + NZ * (IX - 1) + 1;
								if (REGNUM > 0 && REGNUM <= EGS4.$MXREG) {
									EGS4.MED[REGNUM - 1] = MEDNUM[i - 1];
									if (RHORI > 0.0)
										EGS4.RHOR[REGNUM - 1] = RHORI;
									if (DESCRIBE == DESCRIBE_PLANES_DENSITY) {
										if (RHORI > 0.0) {
											EGS4.seqStr = " ZPLANE: "
													+ IZ
													+ " ,RADIUS: "
													+ IX
													+ " = MATERIAL: "
													+ MEDNUM[i - 1]
													+ " RHOR="
													+ EGS4.format(RHORI, 8,
															true);
											if (EGS4.iprint > 1)
												eq.printSequence(EGS4.seqStr);
										} else {
											EGS4.seqStr = " ZPLANE: " + IZ
													+ " ,RADIUS: " + IX
													+ " = MATERIAL: "
													+ MEDNUM[i - 1]
													+ " RHOR= DEFAULT";
											if (EGS4.iprint > 1)
												eq.printSequence(EGS4.seqStr);
										}
									} else {
										EGS4.seqStr = " ZPLANE: " + IZ
												+ " ,RADIUS: " + IX
												+ " = MATERIAL: "
												+ MEDNUM[i - 1];
										if (EGS4.iprint > 1)
											eq.printSequence(EGS4.seqStr);
									}
								}
							}
						}

					}
				}// if (MEDNUM[0]>=0)

			}

			NSUMCV_MAX = nreg - 1;
			ISUMCV_MAX = nreg;
			ISUMCV_DEFAULT = nreg;

			if (NSUMCV < NSUMCV_MIN || NSUMCV > NSUMCV_MAX)
				NSUMCV = NSUMCV_DEFAULT;// 'NUMBER OF CAVITY REGIONS';
			for (int i = 1; i <= NSUMCV; i++) {
				if (ISUMCV[i - 1] < ISUMCV_MIN || ISUMCV[i - 1] > ISUMCV_MAX)
					ISUMCV[i - 1] = ISUMCV_DEFAULT;// 'REGION NUMBERS OF THE
													// CAVITY'
			}

			if (NSUMCV == nreg - 1) {// "WHOLE GEOMETRY IS A CAVITY ZONE, SCORE EVERYWHERE"
				for (int J = 2; J <= nreg; J++) {
					ISUMCV[J - 2] = J;
				}
			}
			if (NSUMCV == 1) {
				EGS4.seqStr = " THE CAVITY ZONE REGION NUMBER IS: " + ISUMCV[0];
				if (EGS4.iprint > 1)
					eq.printSequence(EGS4.seqStr);
			} else {
				EGS4.seqStr = " THE " + NSUMCV + " CAVITY ZONES ARE: ";
				if (EGS4.iprint > 1)
					eq.printSequence(EGS4.seqStr);

				for (int I = 1; I <= NSUMCV; I++)
					EGS4.seqStr = "   " + ISUMCV[I - 1];
				if (EGS4.iprint > 1)
					eq.printSequence(EGS4.seqStr);
			}

		}// ITERSE->0 sau 1

		// "SET THE GEOMETRY FLAG FOR THE CAVITY"
		for (int J = 1; J <= nreg; J++) {
			ntrack[J - 1] = 0;
		}// "ASSUME IT IS NON-CAVITY EVERYWHERE"

		// "SET THE CAVITY FLAG FOR THE CAVITY REGIONS"
		for (int J = 1; J <= NSUMCV; J++) {
			ntrack[ISUMCV[J - 1] - 1] = 1;
		}

	}// GEOMRZ()

	// ;"******************************************************************************
	// "
	// " *******************
	// " * *
	// " * geomsph.mortran *
	// " * *
	// " *******************
	// "
	// "
	// " THIS SUBROUTINE IS USED TO DO ALL THE WORK ASSOCIATED WITH HAVING
	// " A SPHERICAL GEOMETRY INPUT.
	// "
	// " VERSION 1.1 J. Treurniet 05/99
	// "
	// "
	// "===============================================================================
	// "
	// " CONES NOW SUPPORTED. SINCE 90 DEGREES CONE IS NEEDED BY THE CONE
	// GEOMETRY
	// " CHECKING MACRO, IF THE USER DOES NOT INCLUDE THE 90 DEGREES CONE, THE
	// NEXT
	// " ANGLE GREATER THAN 90 IS AUTOMATICALLY FORCED TO BE 90.
	// "
	// " ANGLES AND RADII CAN BE ENTERED INDIVIDUALLY OR IN GROUPS.
	// "
	// "
	// " VERSION 2.0 E. Mainegra-Hing 10/10/2003
	// "
	// "*******************************************************************************
	// "*******************************************************************************
	// "* (1) * Integers between parentheses show the value of the internal *
	// "* * variable corresponding to this input. These are for reference only.
	// *
	// "* * *
	// "* (M) * The 'M' indicates that the variable at hand has multiple *
	// "* * input capability. One may assign an arbitrary number of *
	// "* * values to that input. *
	// "* * E.g.: NSLAB= 2, 10, 2, 16... *
	// "* * *
	// "* (M2) * The 'M' with an integer beside it means that the variable *
	// "* * has that number of inputs. *
	// "* * E.g.: RANDOM NUMBER SEEDS= 97, 33 *
	// "* * *
	// "* (I) * Regular (one number) integer input value. *
	// "* * E.g.: SOURCE NUMBER= 0 *
	// "* * *
	// "* (R) * Regular (one number) real input value. *
	// "* * E.g.: Z OF FRONT FACE= -1000.0 *
	// "* * *
	// "* (C) * Regular (one string) character input value(no ; or #). *
	// "* * E.g.: TITLE= NRCC EGS4 simulation *
	// "* * *
	// "*******************************************************************************
	// ;
	// "*******************************************************************************
	// " SPHERICAL GEOMETRY INPUT
	// " **************************
	// "*******************************************************************************
	// "
	// " GEOMSPH DELIMETERS: :start geometrical inputs:
	// " :stop geometrical inputs:
	// "
	// "
	// " NUMBER OF CONES (M) number of cones (individual or by group)
	// " If omitted or ZERO, pure spherical geometry
	// " assumed.
	// "
	// " ANGLES (M) ANGLES defining the geometry (reals)
	// " No needed in pure spherical geometries.
	// "
	// " For group input there must be as many entries
	// " as for the NUMBER OF CONES, i.e. :
	// " NCON1,NCON2,...,NCONn
	// " DANG1,DANG2,...,DANGn
	// "
	// " For individual input, ncones must be equal
	// " to the number of entries, i.e.:
	// " ncones
	// " DANG1, DANG2,...,DANGncones
	// "
	// " NUMBER OF SPHERES (M) number of spheres (individual or by group)
	// "
	// " For individual inputs, number of spheres
	// " can be omitted
	// "
	// "
	// " RADII (M) radii of spheres defining the geometry (reals)
	// "
	// " For group input there must be as many entries
	// " as for the NUMBER OF SPHERES, i.e. :
	// " NSPH1,NSPH2,...,NSPHn
	// " DRAD1,DRAD2,...,DRADn
	// "
	// " CAVITY ZONES (M) geometrical zone numbers in the cavity (reals)
	// "
	// "-------------------------------------------------------------------------------
	// " MATERIAL INPUT
	// " **************
	// "
	// " MEDIA (M) material name which must match that in the
	// " pegs4 data set EXACTLY, including case.
	// " 24 characters max per medium, ended by , or ;
	// "
	// " MEDNUM (M) the material number (integers)
	// " (MEDNUM=0 => vacuum)
	// " START REGION (M) initial geometrical zone(irl) (integers) for
	// " this medium [NREGLO]
	// " STOP REGION (M) final geometrical zone(irl) (integers) for
	// " this medium.[NREGHI]
	// " ( >NREGLO to input more than one zone)
	// " DEFAULTS: MEDNUM=0 FOR REGION=1 (i.e. VACUUM)
	// " MEDNUM=1 FOR REGION=2,NREG
	// "
	// " These inputs should be thought of as triplets of
	// " MEDNUM,START and STOP REGIONs which are used
	// " to specify the medium numbers for all regions where
	// " the medium is not the default (medium 1).
	// "
	// "***************************************************************************"

	/**
	 * Build spherical geometry.
	 */
	public static void GEOMSPH() {
		// 'NUMBER OF CONES'
		// OBS must supply NCONES[1],..,NCONES[nCONES]!!!!!!!!!!!!!!!!!!!!!
		// idem ANGLES[1].and NRADIUS[1]..RADII[1]
		for (int i = 1; i <= nCONES; i++) {
			if (NCONES[i] < nCONES_MIN || NCONES[i] > nCONES_MAX)
				NCONES[i] = nCONES_DEFAULT;
		}

		ALPHA[0] = 0.0;
		int NCONE = 0;// "get total number of cones"
		for (int i = 1; i <= nCONES; i++)// DO I = 1, NVALUE(NUM_CONES)[
		{
			NCONE = NCONE + NCONES[i];
			if (NCONE > $MAXCANGLE) {
				EGS4.STOPPROGRAM = true;
				EGS4.seqStr = "ERROR: Number of cones " + NCONE
						+ " is greater than $MAXCANGLE=" + $MAXCANGLE;
				eq.printSequence(EGS4.seqStr);
				return;
			}
		}

		EGS4.seqStr = "NUMBER OF ANGLES TO BE INPUT NCONE = " + NCONE;
		if (EGS4.iprint > 1)
			eq.printSequence(EGS4.seqStr);

		NC = NCONE + 1;// "number of conical sections = total number of cones + 1"
		if (NC == 1) {// "PURE SPHERICAL GEOMETRY"
			ALPHA[1] = 180.0; // "NONE CONES"
			NPLAN1 = 1;
			NPLAN2 = 2;
		} else if (NC == 2) {// "HEMISPHERICAL GEOMETRY"
			ALPHA[1] = 90.0; // "ONE CONE ONLY AT 90 DEGREES"
			NPLAN1 = 1;
			NPLAN2 = 2;
		} else if (NC > 2) {// "USER DEFINES NC-2 CONES"
							// 'ANGLES'
			for (int i = 1; i <= nANGLES; i++) {
				if (ANGLES[i] < nANGLES_MIN || ANGLES[i] > nANGLES_MAX)
					ANGLES[i] = nANGLES_DEFAULT;
			}

			if (nANGLES == nCONES) { // "group input: NCON1,NCON2,...,NCONn"
										// "DANG1,DANG2,...,DANGn"
				NCONES[0] = 0;
				int K = 0;
				for (int i = 1; i <= nCONES; i++)// DO I = 1, NVALUE(NUM_CONES)[
				{
					K = K + NCONES[i - 1];
					for (int J = 1; J <= NCONES[i]; i++)// DO J = 1,NCONES(I)[
					{
						ALPHA[K + J] = ALPHA[K + J - 1] + ANGLES[i];
					}
				}
			} else if (nANGLES == NCONE) {
				// "individual input: NCONE"
				// "                  ANG1,ANG2,...,ANGncone"
				for (int I = 1; I <= NCONE; I++) {
					ALPHA[I] = ANGLES[I];// VALUE(NUM_ANGLE,I);
				}
			} else {// "input error for cones"
				EGS4.STOPPROGRAM = true;
				EGS4.seqStr = "ERROR IN GEOMETRICAL INPUT FOR CONES";
				eq.printSequence(EGS4.seqStr);
				EGS4.seqStr = "MISMATCH BETWEEN NUMBER OF CONES AND NUMBER OF ANGULAR ENTRIES";
				eq.printSequence(EGS4.seqStr);

				return;
			}

			boolean NoNinety = true;// "assume initially 90 degree not included"
			int IC = 0; // "LOOP INDEX"
			while (true) {// "CHECK IF USER PUT IN 90 DEGREE CONE"
				IC = IC + 1;
				if (ALPHA[IC] == 90.0) {// "user included 90 degree, great!!!"
					NPLAN1 = IC;
					NPLAN2 = NPLAN1 + 1;
					NoNinety = false;

					EGS4.seqStr = "WARNING: 90 DEGREE CONE INPUT";
					if (EGS4.iprint > 1)
						eq.printSequence(EGS4.seqStr);

					break;// EXIT;
				}

				if (IC == NCONE)
					break;
			}// UNTIL (IC = NCONE);

			if (NoNinety) {
				if (NC > $MAXCANGLE) {// "remeber NC = NCONE + 1"
					EGS4.STOPPROGRAM = true;
					EGS4.seqStr = "ERROR IN GEOMETRICAL INPUT FOR CONES";
					eq.printSequence(EGS4.seqStr);
					EGS4.seqStr = "ABOUT TO ADD 90 DEGREES, BUT DIMENSIONS EXCEEDED MAXIMUM NUMBER OF ANGLES.";
					eq.printSequence(EGS4.seqStr);

					return;
				}
				IC = NCONE; // "LOOP INDEX"
				while (true) {// "PUT IN THE 90 DEGREE CONE"
					if (ALPHA[IC] > 90.0) {// "SHIFT UPWARDS
						ALPHA[IC + 1] = ALPHA[IC];
					} else {// "PUT IN THE 90 DEGREE CONE RIGHT HERE"
						ALPHA[IC + 1] = 90.0;
						NPLAN1 = IC + 1;
						NPLAN2 = NPLAN1 + 1;
						NCONE = NCONE + 1;
						NC = NCONE + 1;

						EGS4.seqStr = " FORCING ANGLE # " + NPLAN1
								+ " TO BE 90 DEGREES";
						if (EGS4.iprint > 1)
							eq.printSequence(EGS4.seqStr);
						break;
					}
					IC = IC - 1;

					if (IC == 0)
						break;
				}// UNTIL (IC = 0);

				if (IC == 0) {// "FELL THROUGH THE LOOP => ALL INPUT CONES > 90"
					ALPHA[1] = 90.0;
					NPLAN1 = 1;
					NPLAN2 = 2;

					EGS4.seqStr = " FORCING ANGLE # " + NPLAN1
							+ " TO BE 90 DEGREES";
					if (EGS4.iprint > 1)
						eq.printSequence(EGS4.seqStr);
				}
			}
		} // "END OF ANGLE INPUT FOR NC>2"

		for (int IC = 1; IC <= NCONE; IC++) {
			if (ALPHA[IC] > 180.0) {

				EGS4.STOPPROGRAM = true;
				EGS4.seqStr = "ERROR IN GEOMETRICAL INPUT FOR CONES";
				eq.printSequence(EGS4.seqStr);
				EGS4.seqStr = " POLAR ANGLES CAN NOT BE GREATER THAN 180 DEGREES";
				eq.printSequence(EGS4.seqStr);

				return;
			}
		}

		// "0 and 180 degrees are the boundaries(i.e. 0th and NCth), "
		// "if those values are entered by input let's use them."
		// "ALPHA(0) =   0.0; set above"

		if (ALPHA[1] == 0.0) {// "shift down angle array since 0 degree is 0th element"
			for (int IC = 1; IC <= NCONE - 1; IC++)// DO IC=1,NCONE-1
			{
				ALPHA[IC] = ALPHA[IC + 1];
			}
			NC = NCONE;
			NCONE = NCONE - 1;
		}
		if (ALPHA[NCONE] == 180.0) {// "user did the work, let's use it !!!"
			NC = NCONE;
			NCONE = NCONE - 1;
		}
		// "should put here a check whether NC > $MAXCANGLE!!!!!!!!!!!"
		ALPHA[NC] = 180.0;// "just in case it didn't fall in the IF"

		double ANGRAD = 0.0;
		if (NC == 1) {// "PURE SPHERICAL GEOMETRY"
			EGS4.seqStr = " ==> This is a pure spherical geometry !";
			if (EGS4.iprint > 1)
				eq.printSequence(EGS4.seqStr);

			for (int IC = 0; IC <= NC; IC++)// DO IC=0,NC
			{// "O and 180 degree included"
				ANGRAD = (Math.PI / 180.) * ALPHA[IC];
				cosalp[IC] = Math.cos(ANGRAD);// 0 biased
				SINALP[IC] = Math.sin(ANGRAD);
				TANALP[IC] = Math.tan(ANGRAD);
				TANAL2[IC] = TANALP[IC] * TANALP[IC];
			}
		} else if (NC > 1) {
			EGS4.seqStr = " CONE OPENING ANGLES:";
			if (EGS4.iprint > 1)
				eq.printSequence(EGS4.seqStr);
			for (int IC = 1; IC <= NC; IC++) {
				EGS4.seqStr = EGS4.format(ALPHA[IC], 10, true);
				if (EGS4.iprint > 1)
					eq.printSequence(EGS4.seqStr);
			}

			NPLAN2 = NPLAN1 + 1;
			for (int IC = 0; IC <= NC; IC++)// DO IC=0,NC
			{// "O and 180 degree included"
				ANGRAD = (Math.PI / 180.) * ALPHA[IC];
				cosalp[IC] = Math.cos(ANGRAD);
				SINALP[IC] = Math.sin(ANGRAD);
				TANALP[IC] = Math.tan(ANGRAD);
				TANAL2[IC] = TANALP[IC] * TANALP[IC];
			}
		}

		EGS4.seqStr = "NUMBER OF CONICAL REGIONS NC = " + NC;
		if (EGS4.iprint > 1)
			eq.printSequence(EGS4.seqStr);

		EGS4.seqStr = "REGION WHERE 90 deg IS UPPER CONE :" + NPLAN1;
		if (EGS4.iprint > 1)
			eq.printSequence(EGS4.seqStr);

		// NUMBER OF SPHERES
		for (int i = 1; i <= nRADIUS; i++) {
			if (NRADIUS[i] < nRADIUS_MIN || NRADIUS[i] > nRADIUS_MAX)
				NRADIUS[i] = nRADIUS_DEFAULT;
		}

		NR = 0;// "get total number of spheres"
		for (int I = 1; I <= nRADIUS; I++)// DO I = 1, NVALUE(NUM_RADII)
		{
			// NRADIUS[I] = VALUE(NUM_RADII,I);
			NR = NR + NRADIUS[I];
			if (NR > $MAXRADII) {
				EGS4.STOPPROGRAM = true;
				EGS4.seqStr = "ERROR:Number of spheres is greater than $MAXRADII!";
				eq.printSequence(EGS4.seqStr);
				return;
			}
		}

		EGS4.seqStr = NR + " spheres in the problem ...";
		if (EGS4.iprint > 1)
			eq.printSequence(EGS4.seqStr);
		// RADII
		for (int i = 1; i <= nRADII; i++) {
			if (RADII[i] < RADII_MIN || RADII[i] > RADII_MAX)
				RADII[i] = RADII_DEFAULT;
		}

		if (NR == 0) {
			NR = nRADII;// NVALUE(NUM_RSPH);

			EGS4.seqStr = NR + " spheres in the problem ...";
			if (EGS4.iprint > 1)
				eq.printSequence(EGS4.seqStr);
		}

		RSPH[0] = 0.0;

		if (nRADIUS == nRADII) { // "group input: NRAD1,NCON2,...,NCONn"
									// "DRAD1,DRAD2,...,DRADn"
			NRADIUS[0] = 0;
			int K = 0;
			for (int I = 1; I <= nRADIUS; I++)// NVALUE(NUM_RADII)[
			{
				K = K + NRADIUS[I - 1];
				for (int J = 1; J <= NRADIUS[I]; J++) {
					RSPH[K + J] = RSPH[K + J - 1] + RADII[I];// VALUE(NUM_RSPH,I);
					RSPH2[K + J - 1] = RSPH[K + J] * RSPH[K + J];
				}
			}
		} else if (nRADII == NR) { // "individual input: NR"
									// "              RAD1,RAD2,...,RADnr"
			for (int I = 1; I <= NR; I++) {
				RSPH[I] = RADII[I];// VALUE(NUM_RSPH,I);
				RSPH2[I - 1] = RSPH[I] * RSPH[I];
			}
		} else {// "input error for spheres"
			EGS4.STOPPROGRAM = true;
			EGS4.seqStr = "ERROR IN GEOMETRICAL INPUT FOR SPHERES !!";
			eq.printSequence(EGS4.seqStr);
			EGS4.seqStr = " MISMATCH BETWEEN NUMBER OF SPHERES AND NUMBER OF RADIAL ENTRIES";
			eq.printSequence(EGS4.seqStr);

			return;
		}

		int I = Math.min(100, NR);
		for (int IX = 1; IX <= I; IX++) {
			if (RSPH[IX] == 0.0) {
				EGS4.seqStr = " IMPROPER INPUT. RADIUS OF 0.0 NOT ALLOWED";
				if (EGS4.iprint > 1)
					eq.printSequence(EGS4.seqStr);
			}

			EGS4.seqStr = " RING RADIUS nr." + IX + " :"
					+ EGS4.format(RSPH[IX], 12, true) + " cms";
			if (EGS4.iprint > 1)
				eq.printSequence(EGS4.seqStr);
		}

		nreg = NR * NC + 1;

		// "define cavity zones"
		// 'CAVITY ZONES'
		for (int i = 1; i <= numcavreg; i++) {
			if (cavreg[i - 1] < cavreg_MIN || cavreg[i - 1] > cavreg_MAX)
				cavreg[i - 1] = cavreg_DEFAULT;
		}

		EGS4.seqStr = " number of cavity regions: " + numcavreg;
		if (EGS4.iprint > 1)
			eq.printSequence(EGS4.seqStr);

		// "CHECK THAT THE NUMBER OF CAVITY ZONES DOES NOT EXCEED ITS BOUNDS"
		if (numcavreg <= 0) {// "no cavity defined, just dose in all the regions desired"
			for (I = 1; I <= nreg; I++) {
				cavreg[I - 1] = 0;
			}

			EGS4.seqStr = " no cavity regions defined !!";
			if (EGS4.iprint > 1)
				eq.printSequence(EGS4.seqStr);
		} else if (numcavreg >= nreg) {
			// "IF IT DOES, REVERT TO THE STANDARD CHAMBER CONFIGURATION"
			numcavreg = 1;
			cavreg[0] = 3;
			for (I = 2; I <= nreg; I++) {
				cavreg[I - 1] = 0;
			}
			// DO I=2,NREG[cavreg(I)=0;]
			EGS4.seqStr = " TOO MANY CAVITY ZONES, REVERTING TO STANDARD CHAMBER";
			if (EGS4.iprint > 1)
				eq.printSequence(EGS4.seqStr);
		} else if (numcavreg == (nreg - 1)) {
			// "WHOLE GEOMETRY IS A CAVITY ZONE, SCORE EVERYWHERE"
			EGS4.seqStr = " WHOLE GEOMETRY IS A CAVITY ZONE, SCORING EVERYWHERE";
			if (EGS4.iprint > 1)
				eq.printSequence(EGS4.seqStr);
			numcavreg = nreg - 1;// DO I=2,NREG[cavreg(I-1)=I;]
			for (I = 2; I <= nreg; I++) {
				cavreg[I - 2] = I;
			}
		} else {// "check cavity regions are within the right limits"
			I = 0;
			while (true) {
				I = I + 1;
				if ((cavreg[I - 1] <= 1) || (cavreg[I - 1] > nreg)) {// "wrong cavity region number"
					EGS4.seqStr = " wrong cavity region number: "
							+ cavreg[I - 1]
							+ " INAPPROPRIATE CAVITY ZONES, REVERTING TO STANDARD CHAMBER";
					if (EGS4.iprint > 1)
						eq.printSequence(EGS4.seqStr);

					numcavreg = 1;
					cavreg[0] = 3;
					for (I = 2; I <= nreg; I++) {
						cavreg[I - 1] = 0;
					}
					// DO I=2,NREG[cavreg(I)=0;]
					break;// "THEN EXIT THE LOOP"
				} else {
					EGS4.seqStr = "  ===> cavity region: "
							+ I
							+ " is region : "
							+ cavreg[I - 1]
							+ " INAPPROPRIATE CAVITY ZONES, REVERTING TO STANDARD CHAMBER";
					if (EGS4.iprint > 1)
						eq.printSequence(EGS4.seqStr);
				}

				if (I >= numcavreg)
					break;

			}// WHILE(I.LT.numcavreg);
		}// "CAVITY ZONES ARE DEFINED"

		// " MATERIAL INPUT
		// " **************
		EGS4.NMED = nMEDIA;
		for (int i = 1; i <= EGS4.NMED; i++) {
			EGS4.seqStr = " MEDIUM: " + i + "  " + EGS4.MEDIA[i - 1];
			if (EGS4.iprint > 1)
				eq.printSequence(EGS4.seqStr);
		}

		if (EGS4.NMED < 1 || EGS4.NMED > EGS4.$MXMED) {
			EGS4.NMED = 1;

			EGS4.seqStr = "NO MEDIUM OR TOO MANY MEDIA. RESET TO ONE MEDIUM INPUT!";
			if (EGS4.iprint > 1)
				eq.printSequence(EGS4.seqStr);

		}

		MEDNUM_MAX = EGS4.NMED;
		for (int i = 1; i <= nMEDNUM; i++) {
			if (MEDNUM[i - 1] < MEDNUM_MIN || MEDNUM[i - 1] > MEDNUM_MAX)
				MEDNUM[i - 1] = MEDNUM_DEFAULT;
		}

		EGS4.MED[0] = 0;
		for (I = 2; I <= nreg; I++) {
			EGS4.MED[I - 1] = 1;
		} // "defaults"
		EGS4.seqStr = " Number of geometrical zones = " + nreg
				+ " ,vacuum in first region";
		if (EGS4.iprint > 1)
			eq.printSequence(EGS4.seqStr);

		NREGLO_MAX = nreg;
		NREGHI_MAX = nreg;

		if (nNREGLO != nMEDNUM || nNREGHI != nMEDNUM) {
			EGS4.STOPPROGRAM = true;
			EGS4.seqStr = " MEDNUM AND START AND STOP REGIONS MUST HAVE THE SAME NUMBER OF VALUES!";
			eq.printSequence(EGS4.seqStr);
			return;
		}
		for (int i = 1; i <= nMEDNUM; i++) {
			if (NREGLO[i - 1] < NREGLO_MIN || NREGLO[i - 1] > NREGLO_MAX)
				NREGLO[i - 1] = NREGLO_DEFAULT;
			if (NREGHI[i - 1] < NREGHI_MIN || NREGHI[i - 1] > NREGHI_MAX)
				NREGHI[i - 1] = NREGHI_DEFAULT;
		}
		// ---------------------------------------
		// "In the following, we allow for vacuum input"
		if (MEDNUM[0] >= 0)// (VALUE(NUM_MEDNUM,1) >= 0)
		{
			for (int i = 1; i <= nMEDNUM; i++) {
				RHORI = 0.0;
				if (NREGHI[i - 1] <= NREGLO[i - 1]) {
					if (NREGLO[i - 1] > 0 && NREGLO[i - 1] <= EGS4.$MXREG) {
						// MED(NREGLO)=MEDNUM;
						EGS4.MED[NREGLO[i - 1] - 1] = MEDNUM[i - 1];
						EGS4.seqStr = " REGION: " + NREGLO[i - 1]
								+ " = MATERIAL: " + MEDNUM[i - 1];
						if (EGS4.iprint > 1)
							eq.printSequence(EGS4.seqStr);
					}

				} else {
					for (int K = NREGLO[i - 1]; K <= NREGHI[i - 1]; K++) {
						if (K > 0 && K <= EGS4.$MXREG) {
							EGS4.MED[K - 1] = MEDNUM[i - 1];
						}
					}

					EGS4.seqStr = " REGION: " + NREGLO[i - 1] + " TO REGION "
							+ NREGHI[i - 1] + " = MATERIAL: " + MEDNUM[i - 1];
					if (EGS4.iprint > 1)
						eq.printSequence(EGS4.seqStr);
				}

			}// for
		}
		// ================
		// "SET THE GEOMETRY FLAG FOR THE CAVITY"
		for (int J = 1; J <= nreg; J++) {
			ntrack[J - 1] = 0;
		} // "ASSUME IT IS NON-CAVITY EVERYWHERE"

		// "SET THE CAVITY FLAG FOR THE CAVITY REGIONS"
		// DO J=1, numcavreg [NTRACK(cavreg(J))=1;]
		for (int J = 1; J <= numcavreg; J++) {
			ntrack[cavreg[J - 1] - 1] = 1;
		}

	}// END;

	// "MACRO THAT GETS CONE AND RADIUS NUMBERS FROM THE REGION NUMBER
	// REPLACE {$GET-IX-IC(#);} WITH {;IX=({P1}-2)/NC+1; IC={P1}-1-NC*(IX-1);}
	// "MACRO THAT GETS PLANE AND RADIUS NUMBERS FROM THE REGION NUMBER
	// REPLACE {$GET-IX-IZ(#);} WITH {;IX=({P1}-2)/NZ+1; IZ={P1}-1-NZ*(IX-1);}
	/**
	 * GETS CONE AND RADIUS NUMBERS FROM THE REGION NUMBER
	 * @param in in
	 * @return the result
	 */
	public static int GET_IX(int in) {
		int result = 0;
		if (EGS4.iGeom == EGS4.iCavity)// 1=RZ geom
			result = (in - 2) / NZ + 1;
		else if (EGS4.iGeom == EGS4.iCavitySPH)// 2=SPH Geom
			result = (in - 2) / NC + 1;
		return result;

	}

	/**
	 * GETS PLANE AND RADIUS NUMBERS FROM THE REGION NUMBER
	 * @param in in
	 * @return the result
	 */
	public static int GET_IZC(int in) {
		int ix = GET_IX(in);
		int result = 0;
		if (EGS4.iGeom == EGS4.iCavity)
			result = in - 1 - NZ * (ix - 1);
		else if (EGS4.iGeom == EGS4.iCavitySPH)
			result = in - 1 - NC * (ix - 1);
		return result;

	}

	// "MACRO THAT GETS THE GEOMETRY NUMBER FROM THE PLANAR AND RADIAL ZONES
	// REPLACE {$GET-IRL(#,#);} WITH {;IRL={P1}+NZ*({P2}-1)+1;}
	// "MACRO THAT GETS THE GEOMETRY NUMBER FROM THE CONICAL AND RADIAL ZONES
	// REPLACE {$GET-IRL(#,#);} WITH {;IRL={P1}+NC*({P2}-1)+1;}
	/**
	 * GETS THE GEOMETRY REGION NUMBER FROM THE PLANAR/CONICAL AND RADIAL ZONES
	 * @param in1 in1
	 * @param in2 in2
	 * @return the result
	 */
	public static int GET_IRL(int in1, int in2) {
		int result = 0;
		if (EGS4.iGeom == EGS4.iCavity)
			result = in1 + NZ * (in2 - 1) + 1;
		else if (EGS4.iGeom == EGS4.iCavitySPH)
			result = in1 + NC * (in2 - 1) + 1;
		return result;

	}

	/**
	 * Print the geometry summary.
	 */
	public static void GEOMRZ_ISUMRY() {
		if (DESCRIBE == 2 || DESCRIBE == 3) {
			EGS4.seqStr = "=========================================================================";
			if (EGS4.iprint > 1)
				eq.printSequence(EGS4.seqStr);
			EGS4.seqStr = "                   NON-DEFAULT DENSITIES";
			if (EGS4.iprint > 1)
				eq.printSequence(EGS4.seqStr);
			EGS4.seqStr = "=========================================================================";
			if (EGS4.iprint > 1)
				eq.printSequence(EGS4.seqStr);

			EGS4.seqStr = "    Region #         Medium             Rhor  ";
			if (EGS4.iprint > 1)
				eq.printSequence(EGS4.seqStr);
			EGS4.seqStr = "from        to          #             (g/cm**3)";
			if (EGS4.iprint > 1)
				eq.printSequence(EGS4.seqStr);

			RHORI = 0.0;
			int NREGHI = 0;
			int NREGLO = 1;
			for (int I = 2; I <= EGS4Geom.nreg; I++) {// "loop through the geometrical regions"
				if (RHORI > 0.0
						&& (EGS4.RHOR[I - 1] != RHORI || EGS4.MED[I - 1] != EGS4.MED[I - 2])) {
					NREGHI = I - 1;

					EGS4.seqStr = " " + EGS4.format(NREGLO, 3)
							+ EGS4.format(NREGHI, 11)
							+ EGS4.format(EGS4.MED[NREGLO - 1], 11)
							+ EGS4.format(RHORI, 20, true);
					if (EGS4.iprint > 1)
						eq.printSequence(EGS4.seqStr);

					RHORI = 0.0;
				}
				if (EGS4.RHOR[I - 1] > 0.0
						&& EGS4.RHOR[I - 1] != EGS4.RHO[EGS4.MED[I - 1] - 1]
						&& RHORI == 0.0) {
					NREGLO = I;
					RHORI = EGS4.RHOR[I - 1];
				}
				if (I == EGS4Geom.nreg && RHORI > 0.0) {// "needed if ir=NREG has non-default density"
					NREGHI = I;
					EGS4.seqStr = " " + EGS4.format(NREGLO, 3)
							+ EGS4.format(NREGHI, 11)
							+ EGS4.format(EGS4.MED[NREGLO - 1], 11)
							+ EGS4.format(RHORI, 20, true);
					if (EGS4.iprint > 1)
						eq.printSequence(EGS4.seqStr);
				}
			}
		}
	}

	// ###############UTILS for
	// HOWFAR######################################################
	// "THE MACRO REPLACING THE CALL TO CYLNDR
	// "*****************************************************************************
	// "
	// " ***********
	// " * *
	// " * $CYLNDR *
	// " * *
	// " ***********
	// "
	// "MACRO TO BE CALLED BY SUBROUTINE HOWFAR IN THE EGS CODE SYSTEM
	// "A FULLY DOCUMENTED SUBROUTINE VERSION IS CONTAINED IN NRCCAUX.MOR, PART
	// "OF THE STANDARD NRCC DISTRIBUTION
	// "
	// "FOR A PARTICLE TRAVELLING INSIDE TWO CONCENTRIC, INFINITE, RIGHT
	// CYLINDERS,
	// "THIS SUBROUTINE DETERMINES THE MINIMUM DISTANCE IT MUST GO TO HIT A
	// CYLINDER.
	// "THE CYLINDERS ARE ASSUMED TO BE ALIGNED AND CENTERED ALONG THE Z-AXIS.
	// "
	// " SOME VARIABLES
	// " ==============
	// "
	// "{P1} = ICYL = THE NUMBER OF THE OUTER CYLINDER
	// "{P2} = IHITC = 1 => PARTICLE HITS THE OUTER SURFACE
	// " = 0 => PARTICLE MISSES THE SURFACES
	// " =-1 => PARTICLE HITS THE INNER SURFACE
	// "{P3} = TCYL = DISTANCE TO SURFACE IF IT HITS
	// "CYRAD2(ICYL) = RADIUS**2 OF THE OUTER CYLINDER
	// "
	// "THIS CODE IS OPTIMIZED FOR SPEED, NOT SIZE.
	// "
	// "FOR PARTICLES NEAR THE SURFACE, A FIRST ORDER APPROXIMATION IS MADE
	// "FOR EXPRESSIONS LIKE X-SQRT(X**2+EPSILON), WHICH SOMETIMES SETS TCYL=0.
	// "
	// "MACRO VERSION 1 A.F.BIELAJEW NRCC 87/10
	// "
	// "
	// "******************************************************************************
	public static int ihitc = 0;// P2
	public static double tcyl = 0.0;// P3
	public static double ustep = 0.0;// P4
	public static int IWATCH = 0;

	public static int ihitp = 0;
	public static double tplane = 0.0;

	public static int geobug = 0;
	public static int IRL = 0;
	public static int IXX = 0;
	public static int ihits = 0;
	public static double tsph = 0.0;
	public static int ICC = 0;
	public static int ihitco = 0;
	public static double tcone = 0.0;

	/**
	 * TO BE CALLED BY SUBROUTINE HOWFAR. FOR A PARTICLE TRAVELLING INSIDE TWO CONCENTRIC, INFINITE, RIGHT CYLINDERS, 
	 * THIS SUBROUTINE DETERMINES THE MINIMUM DISTANCE IT MUST GO TO HIT A CYLINDER. 
	 * THE CYLINDERS ARE ASSUMED TO BE ALIGNED AND CENTERED ALONG THE Z-AXIS. <p>
	 * IF IHITC = 1 - PARTICLE HITS THE OUTER SURFACE; 0 - PARTICLE MISSES THE SURFACES; -1 - PARTICLE HITS THE INNER SURFACE;<p>
	 * TCYL = DISTANCE TO SURFACE IF IT HITS; CYRAD2(ICYL) = RADIUS**2 OF THE OUTER CYLINDER.
	 * CYLINDR IS BASED ON GIVEN ustep!
	 * @param icyl THE NUMBER OF THE OUTER CYLINDER
	 */
	public static void CYLNDR(int icyl)// $CYLNDR(IX,IHITC,TCYL,ustep)
	{
		int irl = EGS4.IR[EGS4.NP - 1]; // "LOCAL REGION NUMBER"
		// REPLACE {$CYLNDR(#,#,#,#);} WITH
		// {
		double U1 = EGS4.U[EGS4.NP - 1];
		double V1 = EGS4.V[EGS4.NP - 1];
		double A = U1 * U1 + V1 * V1;

		double X1 = 0.0;
		double Y1 = 0.0;
		double B = 0.0;
		double B2 = 0.0;
		double C = 0.0;
		double COUT = 0.0;
		double CIN = 0.0;
		double RAD = 0.0;

		if (A == 0.0) {
			// {P2}=0;{P3}=1.0E30;
			ihitc = 0;
			tcyl = 1.0E30;
		} else {
			X1 = EGS4.X[EGS4.NP - 1];
			Y1 = EGS4.Y[EGS4.NP - 1];
			B = X1 * U1 + Y1 * V1;
			B2 = B * B;
			C = X1 * X1 + Y1 * Y1;
			// COUT=C-CYRAD2({P1});
			COUT = C - CYRAD2[icyl - 1];// 1 biased
			if (COUT > 0.0) {
				// if(IWATCH>0){OUTPUT COUT;(' COUT=',E11.3);}
				COUT = 0.0;
			}
			if (B > 0.0) {
				// {P2}=1;
				ihitc = 1;
				if (COUT / B2 > -1.0E-5) {
					// {P3}=-0.5*COUT/B;
					tcyl = -0.5 * COUT / B;
				} else {
					// {P3}=-COUT/(SQRT(B2-A*COUT)+B);
					tcyl = -COUT / (Math.sqrt(B2 - A * COUT) + B);
				}
				// if( {P3} <= {P4} )
				if (tcyl <= ustep) {
					// {P4} = {P3};
					ustep = tcyl;
					// IF( {P1}+1 <= nr )
					if (icyl + 1 <= NR) {
						EGS4.IRNEW = irl + NZ;
					} else {
						EGS4.IRNEW = 1;
					}
				}
			} else if (B < 0.0) {
				// {P2}=1;
				ihitc = 1;
				if (COUT / B2 > -1.0E-5) {
					// {P3}=-2.*B/A*(1.-0.25*A*COUT/B2);
					tcyl = -2. * B / A * (1. - 0.25 * A * COUT / B2);
				} else {
					// {P3}=(SQRT(B2-A*COUT)-B)/A;
					tcyl = (Math.sqrt(B2 - A * COUT) - B) / A;
				}
				// if( {P3} <= {P4} )
				if (tcyl <= ustep) {
					// {P4} = {P3};
					ustep = tcyl;
					// if( {P1}+1 <= nr )
					if (icyl + 1 <= NR) {
						EGS4.IRNEW = irl + NZ;
					} else {
						EGS4.IRNEW = 1;
					}
				}
				// if({P1}.NE.1)
				if (icyl != 1) {
					// CIN=C-CYRAD2({P1}-1);
					CIN = C - CYRAD2[icyl - 2];
					if (CIN < 0.0) {
						// IF(IWATCH.GT.0){OUTPUT CIN;(' CIN=',E11.3);}
						CIN = 0.0;
					}
					RAD = B2 - A * CIN;
					if (RAD >= 0.0) {
						if (CIN / B2 < 1.0E-5) {
							// {P3}=-0.5*CIN/B;
							tcyl = -0.5 * CIN / B;
						} else {
							tcyl = CIN / (Math.sqrt(RAD) - B);
						}
						// if( {P3} <= {P4} )
						if (tcyl <= ustep) {
							// {P4} = {P3};
							ustep = tcyl;
							EGS4.IRNEW = irl - NZ;
						}
					}
				}
			}
		}
		// ;
		// }->REPLACE
		IRL = irl;// just in case!!!

	}

	// $CYLNDR(IX,IHITC,TCYL);"GET DISTANCE TO CYLINDER"
	// "       IHITC   =  1 => HITS OUTER CYLINDER"
	// "               =  0 => MISSES BOTH CYLINDERS"
	// "               = -1 => HITS INNER CYLINDER"
	// ICYL,IHITC,TCYL
	/**
	 * TO BE CALLED BY SUBROUTINE HOWFAR. FOR A PARTICLE TRAVELLING INSIDE TWO CONCENTRIC, INFINITE, RIGHT CYLINDERS, 
	 * THIS SUBROUTINE DETERMINES THE MINIMUM DISTANCE IT MUST GO TO HIT A CYLINDER. 
	 * THE CYLINDERS ARE ASSUMED TO BE ALIGNED AND CENTERED ALONG THE Z-AXIS. <p>
	 * IF IHITC = 1 - PARTICLE HITS THE OUTER SURFACE; 0 - PARTICLE MISSES THE SURFACES; -1 - PARTICLE HITS THE INNER SURFACE;<p>
	 * TCYL = DISTANCE TO SURFACE IF IT HITS; CYRAD2(ICYL) = RADIUS**2 OF THE OUTER CYLINDER.
	 * 
	 * @param icyl THE NUMBER OF THE OUTER CYLINDER
	 */
	public static void CYLNDR2(int icyl) {
		// REPLACE {$CYLNDR(#,#,#,#);} WITH
		// {
		double U1 = EGS4.U[EGS4.NP - 1];
		double V1 = EGS4.V[EGS4.NP - 1];
		double A = U1 * U1 + V1 * V1;

		double X1 = 0.0;
		double Y1 = 0.0;
		double B = 0.0;
		double B2 = 0.0;
		double C = 0.0;
		double COUT = 0.0;
		double CIN = 0.0;
		double RAD = 0.0;

		if (A == 0.0) {
			// {P2}=0;{P3}=1.0E30;
			ihitc = 0;
			tcyl = 1.0E30;
		} else {
			X1 = EGS4.X[EGS4.NP - 1];
			Y1 = EGS4.Y[EGS4.NP - 1];
			B = X1 * U1 + Y1 * V1;
			B2 = B * B;
			C = X1 * X1 + Y1 * Y1;
			// COUT=C-CYRAD2({P1});
			COUT = C - CYRAD2[icyl - 1];// 1 biased
			if (COUT > 0.0) {
				// if(IWATCH>0){OUTPUT COUT;(' COUT=',E11.3);}
				COUT = 0.0;
			}
			if (B > 0.0) {
				// {P2}=1;
				ihitc = 1;
				if (COUT / B2 > -1.0E-3) {
					// {P3}=-0.5*COUT/B;
					tcyl = -0.5 * COUT / B;
				} else {
					// {P3}=-COUT/(SQRT(B2-A*COUT)+B);
					tcyl = -COUT / (Math.sqrt(B2 - A * COUT) + B);
				}
			} else if (B < 0.0) {
				// {P2}=1;
				ihitc = 1;
				if (COUT / B2 > -1.0E-3) {
					// {P3}=-2.*B/A*(1.-0.25*A*COUT/B2);
					tcyl = -2. * B / A * (1. - 0.25 * A * COUT / B2);
				} else {
					// {P3}=(SQRT(B2-A*COUT)-B)/A;
					tcyl = (Math.sqrt(B2 - A * COUT) - B) / A;
				}
				// if({P1}.NE.1)
				if (icyl != 1) {
					// CIN=C-CYRAD2({P1}-1);
					CIN = C - CYRAD2[icyl - 2];
					if (CIN < 0.0) {
						// IF(IWATCH.GT.0){OUTPUT CIN;(' CIN=',E11.3);}
						CIN = 0.0;
					}
					RAD = B2 - A * CIN;
					if (RAD >= 0.0) {
						ihitc = -1;// {P2}=-1;
						if (CIN / B2 < 1.0E-3) {
							// {P3}=-0.5*CIN/B;
							tcyl = -0.5 * CIN / B;
						} else {
							tcyl = CIN / (Math.sqrt(RAD) - B);
						}
					}
				}
			} else {
				ihitc = 1;// {P2}=1;
				tcyl = Math.sqrt(-COUT / A);// {P3}=SQRT(-COUT/A);
			}

		}
		// ;
		// }->REPLACE

	}

	// "THE MACRO REPLACING THE CALL TO PLANES
	// "******************************************************************************
	// "
	// " ***********
	// " * *
	// " * $PLANES *
	// " * *
	// " ***********
	// "
	// "PROGRAM TO BE CALLED BY HOWFAR IN THE EGS CODE SYSTEM
	// "A FULLY DOCUMENTED SUBROUTINE VERSION IS CONTAINED IN NRCCAUX.MOR, PART
	// "OF THE STANDARD NRCC DISTRIBUTION
	// "
	// "FOR A PARTICLE TRAVELLING INSIDE TWO PARALLEL, INFINITE PLANES WITH
	// NORMALS
	// "ALIGNED ALONG THE Z-AXIS, THIS ROUTINE CALCULATES THE STRAIGHT LINES
	// DISTANCE
	// "IT MUST GO TO HIT ONE OF THE PLANES
	// "
	// " SOME VARIABLES
	// " ==============
	// "
	// "{P1} = THE NUMBER OF THE INNER PLANE (LESSER Z-COORDINATE)
	// " IT MUST BE .GE. 1
	// "{P2} = THE NUMBER OF THE OUTER PLANE (GREATER Z-COORDINATE)
	// " IT MUST BE .GE. 2
	// "{P3} = IHITP = 1 => PARTICLE HITS THE OUTER PLANE
	// " = 0 => PARTICLE MISSES BOTH PLANES
	// " = -1 => PARTICLE HITS THE INNER PLANE
	// "{P4} = TPLANE = DISTANCE TO THE PLANE THAT IT HITS
	// "ZPLANE(IPLANE) = Z-COORDINATE OF THE OUTER PLANE
	// "
	// "MACRO VERSION 1 A.F.BIELAJEW NRCC 87/10
	// "
	// "
	// "******************************************************************************
	// CRYSTAL CLEAR!!->I am worry about handling geometry for simulation of
	// human body,
	// where different organs are ellipsoids cutted by planes etc!!!
	/**
	 * TO BE CALLED BY HOWFAR. FOR A PARTICLE TRAVELLING INSIDE TWO PARALLEL, INFINITE PLANES WITH 
	 * NORMALS ALIGNED ALONG THE Z-AXIS, THIS ROUTINE CALCULATES THE STRAIGHT LINES DISTANCE IT MUST GO TO HIT ONE OF THE PLANES. <p>
	 * IF IHITP = 1 - PARTICLE HITS THE OUTER PLANE; 0 - PARTICLE MISSES BOTH PLANES;  -1 - PARTICLE HITS THE INNER PLANE; <p>
	 * TPLANE = DISTANCE TO THE PLANE THAT IT HITS; ZPLANE(IPLANE) = Z-COORDINATE OF THE OUTER PLANE. 
	 * PLANES IS BASED ON USTEP.
	 * @param iz1 iz1, THE NUMBER OF THE INNER PLANE (LESSER Z-COORDINATE), IT MUST BE .GE. 1
	 * @param iz2 iz2, THE NUMBER OF THE OUTER PLANE (GREATER Z-COORDINATE), IT MUST BE .GE. 2
	 */
	public static void PLANES(int iz1, int iz2)// (IZ,IZ+1,IHITP,TPLANE,ustep)
	{
		int irl = EGS4.IR[EGS4.NP - 1]; // "LOCAL REGION NUMBER"

		// REPLACE {$PLANES(#,#,#,#,#);} WITH {
		double WL = EGS4.W[EGS4.NP - 1];
		if (WL > 0.0) {
			// {P4}=(ZPLANE({P2})-Z(NP))/WL;
			tplane = (ZPLANE[iz2 - 1] - EGS4.Z[EGS4.NP - 1]) / WL;
			// IF( {P4} <= {P5} )
			if (tplane <= ustep) {
				// {P5} = {P4};
				ustep = tplane;
				// if( {P1}+1 <= nz )
				if (iz1 + 1 <= NZ) {
					EGS4.IRNEW = irl + 1;
				} else {
					EGS4.IRNEW = 1;
				}
			}
		} else if (WL < 0.0) {
			// {P4}=(ZPLANE({P1})-Z(NP))/WL;
			tplane = (ZPLANE[iz1 - 1] - EGS4.Z[EGS4.NP - 1]) / WL;
			// IF( {P4} <= {P5} ) [
			if (tplane <= ustep) {
				// {P5} = {P4};
				ustep = tplane;
				// if( {P1}-1 >= 1 )
				if (iz1 - 1 >= 1) {
					EGS4.IRNEW = irl - 1;
				} else {
					EGS4.IRNEW = 1;
				}
			}
		}
		// }
		IRL = irl;// just in case!!!
	}

	/**
	 * TO BE CALLED BY HOWFAR. FOR A PARTICLE TRAVELLING INSIDE TWO PARALLEL, INFINITE PLANES WITH 
	 * NORMALS ALIGNED ALONG THE Z-AXIS, THIS ROUTINE CALCULATES THE STRAIGHT LINES DISTANCE IT MUST GO TO HIT ONE OF THE PLANES. <p>
	 * IF IHITP = 1 - PARTICLE HITS THE OUTER PLANE; 0 - PARTICLE MISSES BOTH PLANES;  -1 - PARTICLE HITS THE INNER PLANE; <p>
	 * TPLANE = DISTANCE TO THE PLANE THAT IT HITS; ZPLANE(IPLANE) = Z-COORDINATE OF THE OUTER PLANE. <p>
	 * Comments: <p>
	 * To do these routines compatible with human body organs where we have complex geometries (ellipsoids, planes cut etc.and not simple RZ geometry) is a real pain. Therefore, 
	 * EGSGeom is not well suited for these tasks, so we must use different toolkit such as GEANT4 where geometries and geometry-MC engine integrations are handled in a much more convenient (and simple) way.    
	 * One major drawback for EGS is its cumbersome design regarding geometry integration (EGS4Geom-EGS4Core by implementing HOWFAR and HOWNEAR). HOWFAR and HOWNEAR are ok for 
	 * simple geometries such as RZ, spherical or pure planar but for complex geometries is quite hard to implement them. For instance write HOWFAR and HOWNEAR for a half-elipsoid 
	 * cut by a plane at given intersection points...good luck with that...it is possible but it is hard, annoying and error prone. 
	 * @param iz1 iz1, THE NUMBER OF THE INNER PLANE (LESSER Z-COORDINATE), IT MUST BE .GE. 1
	 * @param iz2 iz2, THE NUMBER OF THE OUTER PLANE (GREATER Z-COORDINATE), IT MUST BE .GE. 2
	 */
	public static void PLANES2(int iz1, int iz2)// IZ,IZ+1,IHITP,TPLANE
	{
		// REPLACE {$PLANES(#,#,#,#,#);} WITH {
		double WL = EGS4.W[EGS4.NP - 1];
		if (WL > 0.0) {
			ihitp = 1;// {P3}=1;
			// {P4}=(ZPLANE({P2})-Z(NP))/WL;
			tplane = (ZPLANE[iz2 - 1] - EGS4.Z[EGS4.NP - 1]) / WL;
		} else if (WL < 0.0) {
			ihitp = -1;// {P3}=-1
			// {P4}=(ZPLANE({P1})-Z(NP))/WL;
			tplane = (ZPLANE[iz1 - 1] - EGS4.Z[EGS4.NP - 1]) / WL;
		} else {
			ihitp = 0;// {P3}=0;
			tplane = 1.0E30;// {P4}=1.0E30;
		}
		// }
	}

	// "*****************************************************************************
	// "
	// " ***********
	// " * *
	// " * $SPHERE *
	// " * *
	// " ***********
	// "
	// "MACRO TO BE CALLED BY SUBROUTINE HOWFAR IN THE EGS CODE SYSTEM
	// "
	// "FOR A PARTICLE TRAVELLING INSIDE TWO CONCENTRIC, INFINITE, RIGHT
	// SPHERES,
	// "THIS SUBROUTINE DETERMINES THE MINIMUM DISTANCE IT MUST GO TO HIT A
	// SPHERE.
	// "THE SPHERES ARE ASSUMED TO BE CENTERED AT THE ORIGIN.
	// "
	// " SOME VARIABLES
	// " ==============
	// "
	// "{P1} = IX = THE NUMBER OF THE OUTER SPHERE
	// "{P2} = IHITS = 1 => PARTICLE HITS THE OUTER SURFACE
	// " = 0 => PARTICLE MISSES THE SURFACES
	// " =-1 => PARTICLE HITS THE INNER SURFACE
	// "{P3} = TSPH = DISTANCE TO SURFACE IF IT HITS
	// "RSPH2(IX) = RADIUS**2 OF SPHERE
	// "
	// "FOR PARTICLES NEAR THE SURFACE, A FIRST ORDER APPROXIMATION IS MADE
	// "FOR EXPRESSIONS LIKE X-SQRT(X**2+EPSILON), WHICH SOMETIMES SETS TSPH=0.
	// "
	// "MACRO VERSION 1 A.F.BIELAJEW NRCC 88/03
	// " ADAPTED FROM SCASPH (EGS3)
	// "
	// "******************************************************************************
	// "
	/**
	 * TO BE CALLED BY SUBROUTINE HOWFAR. FOR A PARTICLE TRAVELLING INSIDE TWO CONCENTRIC, INFINITE, RIGHT SPHERES, 
	 * THIS SUBROUTINE DETERMINES THE MINIMUM DISTANCE IT MUST GO TO HIT A SPHERE. THE SPHERES ARE ASSUMED TO BE CENTERED AT THE ORIGIN. <p>
	 * IF IHITS = 1 - PARTICLE HITS THE OUTER SURFACE; 0 - PARTICLE MISSES THE SURFACES; -1 - PARTICLE HITS THE INNER SURFACE; TSPH - DISTANCE TO SURFACE IF IT HITS
	 * @param ix ix, THE NUMBER OF THE OUTER SPHERE
	 */
	public static void SPHERE(int ix) {
		IXX = ix;
		boolean startS = false;
		//boolean endS = false;
		IRL = EGS4.IR[EGS4.NP - 1];// "LOCAL REGION NUMBER"
		// "LOCAL VARIABLES"
		double UL = EGS4.U[EGS4.NP - 1];
		double VL = EGS4.V[EGS4.NP - 1];
		double WL = EGS4.W[EGS4.NP - 1];
		double XL = EGS4.X[EGS4.NP - 1];
		double YL = EGS4.Y[EGS4.NP - 1];
		double ZL = EGS4.Z[EGS4.NP - 1];

		double CIN = 0.0;
		double TEST = 0.0;
		double BS = 0.0;
		double COUT = 0.0;
		// REPLACE {$SPHERE(#,#,#);} WITH {
		double B = UL * XL + VL * YL + WL * ZL;// "projection of particle radio vector"
		double B2 = B * B; // "on particle direction"
		double C = XL * XL + YL * YL + ZL * ZL;// "distance of particle to origin"

		// :START-SPHERE:;
		while (true) {
			startS = false;
			//endS = false;
			if (IXX != 1) {
				// "not inner most sphere"
				// CIN=C-RSPH2({P1}-1);
				CIN = C - RSPH2[IXX - 2];// 1 biased
				if (CIN < 0.0) {// "distance to particle smaller than inner radius"
					IRL = IRL - NC;// "reset everything ..."
					EGS4.IR[EGS4.NP - 1] = IRL;// "reset to the corresponding zone"
					// {P1}={P1}-1;//"reset to correspponding spheres"
					IXX = IXX - 1;// "reset to correspponding spheres"
					geobug = 1;// GEOBUG=1; //"signal a bug"
					// GOTO :START-SPHERE:;//"try again"
					startS = true;
				}
				if (!startS) {
					if (B2 < CIN) {
						ihits = 0;
					} // {{P2}=0;}
					else if (B >= 0.) {
						ihits = 0;
					}// {{P2}=0;}
					else {// "B is negative"
						ihits = -1;// {P2}=-1;
						TEST = CIN / B2;
						if (Math.abs(TEST) < 0.001) {
							BS = 2. * B * (1. - 0.25 * TEST);
						} else if (TEST >= 1.) {
							BS = B;
						} else {
							BS = B * (1. + Math.sqrt(1. - TEST));
						}
						// "{P3}=-CIN/BS; since B is < 0.0"
						if (BS != 0.0) {
							tsph = -CIN / BS;// {P3}=-CIN/BS;
						} else if (BS == 0.0) {
							EGS4.STOPPROGRAM = true;
							EGS4.seqStr = " ERROR: caught in sphere macro. BS= "
									+ BS;
							eq.printSequence(EGS4.seqStr);
							return;
						} else {
							EGS4.STOPPROGRAM = true;
							EGS4.seqStr = " ERROR: BS is undefined!! BS= " + BS;
							eq.printSequence(EGS4.seqStr);
							return;
						}
						// GOTO :END-SPHERE:;
						return;// endS=true;
					}
				}// if(!startS)
			}// if(IXX!=1)
			if (!startS) {
				// COUT=C-RSPH2({P1});
				COUT = C - RSPH2[IXX - 1];
				if (COUT > 0.0) {// "distance to particle larger than outer radius"
					IRL = IRL + NC;// "reset everything ..."
					EGS4.IR[EGS4.NP - 1] = IRL;// "update zone"
					// {P1}={P1}+1;//""
					IXX = IXX + 1;// ""
					geobug = 1;// GEOBUG=1;
					// GOTO :START-SPHERE:;
					startS = true;
				}
			}// if(!startS)

			if (!startS)
				break;
		}// while
			// {P2}=1;
		ihits = 1;
		if (B == 0.) {
			tsph = Math.sqrt(-COUT);
		}// [{P3}=SQRT(-COUT);]
		else {
			TEST = COUT / B2;
			if (Math.abs(TEST) < 0.001) {
				BS = 2. * B * (1. - 0.25 * TEST);
			} else if (TEST >= 1.) {
				BS = B;
			} else {
				BS = B * (1. + Math.sqrt(1. - TEST));
			}
			if (B <= 0.) {
				tsph = -BS;
			}// [{P3}=-BS;]
			else {
				if ((BS != 0.0)) {
					tsph = -COUT / BS;// {P3}=-COUT/BS;
				} else if (BS == 0.0) {
					EGS4.STOPPROGRAM = true;
					EGS4.seqStr = " ERROR: caught in sphere macro. BS= " + BS;
					eq.printSequence(EGS4.seqStr);
					return;

					// write(6,*) ' exception caught in sphere macro';
					// write(6,*) ' BS= ', BS;
				} else {
					EGS4.STOPPROGRAM = true;
					EGS4.seqStr = " ERROR: BS is undefined!! BS= " + BS;
					eq.printSequence(EGS4.seqStr);
					return;

					// write(6,*) 'BS is undefined!!! BS= ', BS;
					// write(6,*) 'TEST=', TEST, ' B2=', B2, ' COUT=',COUT;
					// stop;
				}
			}
		}
		// :END-SPHERE:;
		// }
	}

	// "THE MACRO REPLACING THE CALL TO CONES
	// "******************************************************************************
	// "
	// " **********
	// " * *
	// " * $CONES *
	// " * *
	// " **********
	// "
	// "PROGRAM TO BE CALLED BY HOWFAR IN THE EGS CODE SYSTEM
	// "
	// "FOR A PARTICLE TRAVELLING INSIDE TWO CONES WITH AXES
	// "ALIGNED ALONG THE Z-AXIS, THIS ROUTINE CALCULATES THE STRAIGHT LINES
	// DISTANCE
	// "IT MUST GO TO HIT ONE OF THE CONES
	// "
	// " SOME VARIABLES
	// " ==============
	// "
	// "{P1} = IC = THE CONE NUMBER OF THE REGION THE PARTICLE IS IN
	// "{P2} = IHITC = 1 => PARTICLE HITS THE OUTER CONE
	// " = 0 => PARTICLE MISSES BOTH CONES
	// " = -1 => PARTICLE HITS THE INNER CONE
	// "{P3} = TCONE = DISTANCE TO THE CONE THAT IT HITS
	// "TANAL2(IX) = TANGENT**2 OF CONE OPENING ANGLE
	// "
	// "MACRO VERSION 1 A.F.BIELAJEW NRCC 88/03
	// " ADAPTED FROM SCASPH (EGS3)
	// "
	// "******************************************************************************
	// "
	// "------------------------------------------------------------------------------"
	// "                                                                              "
	// "   macro below replaces GOTO statements between IF-BLOCKS avoiding compiler   "
	// "   warnings.                                                                  "
	// "                                                                              "
	// " {P1} => as above, cone number of the region the particle is in               "
	// " {P2} => +/- depending whether region number is to be increased or decreased  "
	// "                                                          EMH, June 5, 2002   "
	// "------------------------------------------------------------------------------"
	// ;
	// REPLACE {$RESET REGION # #;} WITH {;
	// IRL = IRL {P2} 1;
	// IR(NP) = IRL;
	// {P1} = {P1} {P2} 1;
	// GEOBUG = 1;
	// GOTO :START-CONE:;
	// }

	/**
	 * TO BE CALLED BY HOWFAR IN THE EGS. FOR A PARTICLE TRAVELLING INSIDE TWO CONES WITH AXES ALIGNED ALONG THE Z-AXIS, THIS ROUTINE CALCULATES THE STRAIGHT LINES 
	 * DISTANCE IT MUST GO TO HIT ONE OF THE CONES.<p>
	 * IF IHITC = 1 - PARTICLE HITS THE OUTER CONE; 0 - PARTICLE MISSES BOTH CONES; -1 - PARTICLE HITS THE INNER CONE; TCONE - DISTANCE TO THE CONE THAT IT HITS.
	 * @param ic ic, THE CONE NUMBER OF THE REGION THE PARTICLE IS IN
	 */
	public static void CONES(int ic) {
		ICC = ic;
		IRL = EGS4.IR[EGS4.NP - 1];
		int IHITCO = 0;
		int IHITCI = 0;
		double TCONEI = 0.0;
		double TCONEO = 0.0;
		double TANA2 = 0.0;
		double A = 0.0;
		double B = 0.0;
		double C = 0.0;
		double AC = 0.0;
		double B2 = 0.0;
		double TEST = 0.0;
		boolean startC = false;
		double BS = 0.0;
		// "LOCAL VARIABLES"
		double UL = EGS4.U[EGS4.NP - 1];
		double VL = EGS4.V[EGS4.NP - 1];
		double WL = EGS4.W[EGS4.NP - 1];
		double XL = EGS4.X[EGS4.NP - 1];
		double YL = EGS4.Y[EGS4.NP - 1];
		double ZL = EGS4.Z[EGS4.NP - 1];
		// REPLACE {$CONES(#,#,#);} WITH {;

		// :START-CONE:;
		while (true) {
			startC = false;
			IHITCI = 0;
			IHITCO = 0;
			// if({P1}.EQ.NPLAN2)[
			if (ICC == NPLAN2) {
				if (ZL > 0.) {
					// $RESET REGION {P1} -;
					IRL = IRL - 1;// IRL = IRL {P2} 1;
					EGS4.IR[EGS4.NP - 1] = IRL;
					ICC = ICC - 1;// {P1} = {P1} {P2} 1;
					geobug = 1;
					// GOTO :START-CONE:;
					startC = true;
				} // "GOTO :BUG1:;"
				else if (WL <= 0.) {
				} else {
					IHITCI = -1;
					TCONEI = -ZL / WL;
				}
			} else if (ICC != 1)// ({P1}.NE.1)
			{
				TANA2 = TANAL2[ICC - 1];// TANAL2({P1}-1);//0-biased
				A = UL * UL + VL * VL - TANA2 * WL * WL;
				B = UL * XL + VL * YL - TANA2 * WL * ZL;
				C = XL * XL + YL * YL - TANA2 * ZL * ZL;
				AC = A * C;
				B2 = B * B;
				if (C >= 0) {
					if (ZL < 0.0) {
						if (ICC < NPLAN2) {
							// $RESET REGION {P1} +;
							IRL = IRL + 1;// IRL = IRL {P2} 1;
							EGS4.IR[EGS4.NP - 1] = IRL;
							ICC = ICC + 1;// {P1} = {P1} {P2} 1;
							geobug = 1;
							// GOTO :START-CONE:;
							startC = true;
						}
						if (!startC) {
							// $RESET REGION {P1} -;
							IRL = IRL - 1;// IRL = IRL {P2} 1;
							EGS4.IR[EGS4.NP - 1] = IRL;
							ICC = ICC - 1;// {P1} = {P1} {P2} 1;
							geobug = 1;
							// GOTO :START-CONE:;
							startC = true;
						}
					}
					if (!startC) {
						if (B2 < AC) {
						} else if ((A >= 0.) && (B >= 0.)) {
						} else {
							if (B == 0.) {
								TCONEI = Math.sqrt(-C / A);
							} else {
								TEST = AC / B2;
								if (Math.abs(TEST) < 0.001) {
									BS = 2. * B * (1. - 0.25 * TEST);
								} else if (TEST >= 1.) {
									BS = B;
								} else {
									BS = B * (1. + Math.sqrt(1. - TEST));
								}
								if (B > 0.) {
									TCONEI = -BS / A;
								} else {
									TCONEI = -C / BS;
								}
							}
							// TEST=TANALP({P1}-1)*(ZL+WL*TCONEI);//0 biased
							TEST = TANALP[ICC - 1] * (ZL + WL * TCONEI);
							if (TEST < 0.) {
							} else {
								IHITCI = -1;
							}
						}
					}// if(!startC)
				}
				// else if(TANALP({P1}-1)*ZL.GE.0.)[
				else if (TANALP[ICC - 1] * ZL >= 0.) {
					if (ZL > 0.0) {
						// $RESET REGION {P1} -;
						IRL = IRL - 1;// IRL = IRL {P2} 1;
						EGS4.IR[EGS4.NP - 1] = IRL;
						ICC = ICC - 1;// {P1} = {P1} {P2} 1;
						geobug = 1;
						// GOTO :START-CONE:;
						startC = true;
					}
					if (!startC) {
						if ((A <= 0.) && (B <= 0.)) {
						} else {
							if (B == 0.) {
								TCONEI = Math.sqrt(-C / A);
							} else {
								TEST = AC / B2;
								if (Math.abs(TEST) < 0.001) {
									BS = 2. * B * (1. - 0.25 * TEST);
								} else if (TEST >= 1.) {
									BS = B;
								} else {
									BS = B * (1. + Math.sqrt(1. - TEST));
								}
								if (B <= 0.) {
									TCONEI = -BS / A;
								} else {
									TCONEI = -C / BS;
								}
							}
							IHITCI = -1;
						}
					}// if(!startC)
				} else {
					// GOTO :SKIP-THIS:;
					// if((A.LT.0.).AND.(B.GT.0.))[
					// TEST=AC/B2;
					// IF(ABS(TEST).LT.0.001)[BS=2.*B*(1.-0.25*TEST);]
					// ELSEIF(TEST.GE.1.)[BS=B;]
					// ELSE[BS=B*(1.+SQRT(1.-TEST));]
					// TCONEI=-BS/A;
					// IHITCI=-1;
					// ]
					// ELSE[;]
					// :SKIP-THIS:;
				}
			}
			// IF({P1}.EQ.NPLAN1)[
			if (!startC)// 1
			{
				if (ICC == NPLAN1) {
					if (ZL < 0.) {
						// $RESET REGION {P1} +;//here auto loop
						IRL = IRL + 1;// IRL = IRL {P2} 1;
						EGS4.IR[EGS4.NP - 1] = IRL;
						ICC = ICC + 1;// {P1} = {P1} {P2} 1;
						geobug = 1;
						// GOTO :START-CONE:;
						startC = true;

					} else if (WL >= 0.) {
					} else {
						IHITCO = 1;
						TCONEO = -ZL / WL;
					}
				}
				// else if({P1}.NE.NC)[
				else if (ICC != NC) {
					TANA2 = TANAL2[ICC];// TANAL2({P1});
					A = UL * UL + VL * VL - TANA2 * WL * WL;
					B = UL * XL + VL * YL - TANA2 * WL * ZL;
					C = XL * XL + YL * YL - TANA2 * ZL * ZL;
					AC = A * C;
					B2 = B * B;
					if (C >= 0) {
						if (ZL > 0.0) {
							// if({P1}.GT.NPLAN1)
							if (ICC > NPLAN1) {
								// $RESET REGION {P1} -;
								IRL = IRL - 1;// IRL = IRL {P2} 1;
								EGS4.IR[EGS4.NP - 1] = IRL;
								ICC = ICC - 1;// {P1} = {P1} {P2} 1;
								geobug = 1;
								// GOTO :START-CONE:;
								startC = true;

							}
							if (!startC) {
								// $RESET REGION {P1} +;
								IRL = IRL + 1;// IRL = IRL {P2} 1;
								EGS4.IR[EGS4.NP - 1] = IRL;
								ICC = ICC + 1;// {P1} = {P1} {P2} 1;
								geobug = 1;
								// GOTO :START-CONE:;
								startC = true;
							}
						}
						if (!startC) {
							if (B2 < AC) {
							} else if ((A >= 0.) && (B >= 0.)) {
							} else {
								if (B == 0.) {
									TCONEO = Math.sqrt(-C / A);
								} else {
									TEST = AC / B2;
									if (Math.abs(TEST) < 0.001) {
										BS = 2. * B * (1. - 0.25 * TEST);
									} else if (TEST >= 1.) {
										BS = B;
									} else {
										BS = B * (1. + Math.sqrt(1. - TEST));
									}
									if (B > 0.) {
										TCONEO = -BS / A;
									} else {
										TCONEO = -C / BS;
									}
								}
								// TEST=TANALP({P1})*(ZL+WL*TCONEO);
								TEST = TANALP[ICC] * (ZL + WL * TCONEO);
								if (TEST < 0.) {
								} else {
									IHITCO = 1;
								}
							}
						}// if(!startC)
					}
					// ELSEIF(TANALP({P1})*ZL.GE.0.)
					else if (TANALP[ICC] * ZL >= 0.) {
						if (ZL < 0.0) {
							// $RESET REGION {P1} +;
							IRL = IRL + 1;// IRL = IRL {P2} 1;
							EGS4.IR[EGS4.NP - 1] = IRL;
							ICC = ICC + 1;// {P1} = {P1} {P2} 1;
							geobug = 1;
							// GOTO :START-CONE:;
							startC = true;
						}
						if (!startC) {
							if ((A <= 0.) && (B <= 0.)) {
							} else {
								if (B == 0.) {
									TCONEO = Math.sqrt(-C / A);
								} else {
									TEST = AC / B2;
									if (Math.abs(TEST) < 0.001) {
										BS = 2. * B * (1. - 0.25 * TEST);
									} else if (TEST >= 1.) {
										BS = B;
									} else {
										BS = B * (1. + Math.sqrt(1. - TEST));
									}
									if (B <= 0.) {
										TCONEO = -BS / A;
									} else {
										TCONEO = -C / BS;
									}
								}
								IHITCO = 1;
							}
						}// if(!startC)
					} else {
						// GOTO :END-CONE:;
						// IF((A.LT.0.).AND.(B.GT.0.))
						// [
						// TEST=AC/B2;
						// IF(ABS(TEST).LT.0.001)[BS=2.*B*(1.-0.25*TEST);]
						// ELSEIF(TEST.GE.1.)[BS=B;]
						// ELSE[BS=B*(1.+SQRT(1.-TEST));]
						// TCONEO=-BS/A;
						// IHITCO=1;
						// ]
						// ELSE[;]
					}
				}
			}// if(!startC)//1
			if (!startC)
				break;
		}// while
			// :END-CONE:;

		if ((IHITCI == -1) && (IHITCO == 1)) {
			if (TCONEI <= TCONEO)
			// {{P2}=IHITCI;{P3}=TCONEI;}
			{
				ihitco = IHITCI;
				tcone = TCONEI;
			} else
			// {{P2}=IHITCO;{P3}=TCONEO;}
			{
				ihitco = IHITCO;
				tcone = TCONEO;
			}
		} else if (IHITCI == -1)
		// [{P2}=IHITCI;{P3}=TCONEI;]
		{
			ihitco = IHITCI;
			tcone = TCONEI;
		} else if (IHITCO == 1)
		// [{P2}=IHITCO;{P3}=TCONEO;]
		{
			ihitco = IHITCO;
			tcone = TCONEO;
		} else {
			ihitco = 0;
		}// [{P2}=0;]
		// }

	}

	// ###############UTILS for
	// HOWFAR######################################################
}
