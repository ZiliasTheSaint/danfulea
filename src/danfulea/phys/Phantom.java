package danfulea.phys;

import java.util.ResourceBundle;

/**
 * This class simulates the human body using MIRD 5 mathematical phantom. Provides 
 * geometry to Monte Carlo engine (e.g. DoseSimCore).
 *
 * @author Dan Fulea, 10 APR. 2005
 */
public class Phantom {
	private static final String BASE_RESOURCE_CLASS = "danfulea.phys.resources.PhantomResources";
	private static ResourceBundle resources = ResourceBundle
			.getBundle(BASE_RESOURCE_CLASS);

	// --standard dimension of compressed breast->compression device
	private static double BREAST_RADIUS = 7.0; // cm
	private static double BREAST_THICKNESS = 4.2; // cm
	private static double BREAST_MASS = 300.0; // g??-------NEVER USED
	// ---------------------
	public static final int MAMO_INDEX = 0;// breast phantom
	public static final int NEWBORN_INDEX = 1;// newborn human phantom
	public static final int AGE1_INDEX = 2;// age 1 human phantom
	public static final int AGE5_INDEX = 3;// age 5 human phantom
	public static final int AGE10_INDEX = 4;// age 10 human phantom
	public static final int AGE15_INDEX = 5;// age 15 human phantom
	public static final int ADULT_INDEX = 6;// adult human phantom
	private int mode = 0;
	// -----------------------------------
	public static final int NO_ORGAN = 0;
	public static final int BREASTS = 1;
	public static final int TESTES = 2;
	public static final int SKELETON = 3;
	public static final int ACTVEBONEMARROW = 4;
	public static final int ADRENALS = 5;
	public static final int BRAIN = 6;
	public static final int GALLBLADDER = 7;
	public static final int STOMACH = 8;
	public static final int SMALLINTESTINE = 9;
	public static final int UPPERLARGEINTESTINE = 10;
	public static final int LOWERLARGEINTESTINE = 11;
	public static final int HEART = 12;
	public static final int KIDNEYS = 13;
	public static final int LIVER = 14;
	public static final int LUNGS = 15;
	public static final int OVARIES = 16;
	public static final int PANCREAS = 17;
	public static final int SPLEEN = 18;
	public static final int THYMUS = 19;
	public static final int THYROID = 20;
	public static final int URINARYBLADDER = 21;
	public static final int UTERUS = 22;
	public static final int REMAINDER = 23;

	// if flag is true->it is used for manual setting several features
	/**
	 * Constructor. Modes available are: MAMO_INDEX (for simulating the breast in mammography), NEWBORN_INDEX, AGE1_INDEX, 
	 * AGE5_INDEX, AGE10_INDEX, AGE15_INDEX and ADULT_INDEX.
	 * @param mode the mode index
	 */
	public Phantom(int mode) {
		this.mode = mode;
	}

	/**
	 * 
	 * @return the mode index
	 */
	public int getIndex() {
		return this.mode;
	}

	// ==**********************************************************************************************
	/**
	 * Find the minimum of 3 doubles.
	 * @param a a
	 * @param b b
	 * @param c c
	 * @return the result
	 */
	public static double min(double a, double b, double c) {
		double r = Math.min(a, b);
		r = Math.min(r, c);
		return r;
	}

	/**
	 * The minimum distance from a point inside the cylinder to the cylinder surface.
	 * @param x x
	 * @param y y
	 * @param z z
	 * @param irl irl not used though
	 * @return the result
	 */
	public double getTperp(double x, double y, double z, int irl) {
		double result = 0.0;
		if (mode == MAMO_INDEX)// only one organ
		{
			double r = Math.sqrt(x * x + y * y);
			result = min(z - 0.0, BREAST_THICKNESS - z, BREAST_RADIUS - r);
		}

		return result;
	}

	public static int ihitp = 0;

	/**
	 * Using directional cosine w and z coordinate, computes the distance from the point 
	 * to the base surface (top or bottom) exit (distance traveled along w direction).
	 * @param w w
	 * @param z z
	 * @return the result
	 */
	public double getMamoTz(double w, double z)// w is the directional cosine
												// relative of z axis
	{
		// w is EGS4.W[EGS4.NP-1];
		double tplane = 1.0E30;
		ihitp = 0;
		if (w > 0.0) {
			ihitp = 1;
			tplane = (BREAST_THICKNESS - z) / w;
		} else if (w < 0.0) {
			ihitp = -1;
			tplane = (0.0 - z) / w;
		}

		return tplane;
	}

	/*
	 * //With emotions,( :)) ), found the following correct for our case!!
	 * double l1=0.0; double
	 * zer=(ux*x0+uy*y0)*(ux*x0+uy*y0)+(ux*ux+uy*uy)*(asource
	 * *asource-x0*x0-y0*y0); if (zer>0.) { double
	 * s1=(-(ux*x0+uy*y0)-Math.sqrt((ux*x0+uy*y0)*(ux*x0+uy*y0)+
	 * (ux*ux+uy*uy)*(asource*asource-x0*x0-y0*y0)))/(ux*ux+uy*uy); double
	 * s2=(-(ux*x0+uy*y0)+Math.sqrt((ux*x0+uy*y0)*(ux*x0+uy*y0)+
	 * (ux*ux+uy*uy)*(asource*asource-x0*x0-y0*y0)))/(ux*ux+uy*uy); double s=0.;
	 * if ((s1<0.) && (s2>0)){ s=s2;} else
	 * {s=Math.min(Math.abs(s1),Math.abs(s2));} l1=s;//System.out.println(l1); }
	 * else l1=source_parcurs;//2;//something wrong is happened!! if
	 * (l1<source_parcurs)//2) //it also fly in air (neglected)
	 * source_parcurs=l1;//source_parcurs2=l1;
	 */
	public static int ihitc = 0;

	/**
	 * Same as getMamoTz but here the distance traveled is until it exists the cylinder 
	 * lateral surface. It uses directional cosines u, v and its associated coordinates x, y. 
	 * @param u u
	 * @param v v
	 * @param x x
	 * @param y y
	 * @return the result
	 */
	public double getMamoTxy(double u, double v, double x, double y)// hitting
																	// cylindersurface?
	{
		double txy = 1.0E30;
		ihitc = 0;

		double A = u * u + v * v;// see numitor in gammadeteff, getCylRandom,
									// zer
		if (A != 0) {
			double B = u * x + v * y;
			double B2 = B * B;
			double C = x * x + y * y;
			double COUT = C - BREAST_RADIUS * BREAST_RADIUS;

			double zer = B2 - A * COUT;
			if (zer > 0.)// this case, COUT<0 (inside cyl)=>sqrt(zer)>abs(B).
			{
				ihitc = 1;
				double s = 0.0;
				double s1 = (-B - Math.sqrt(zer)) / A;
				double s2 = (-B + Math.sqrt(zer)) / A;
				if ((s1 <= 0.0) && (s2 >= 0.0)) {
					s = s2;
				}// allways happen!! because:
				// B>0=>s2>0,s1<0
				// B<0=>s1<0 (sqrt(zer)>abs(B)),s2>0
				else {
					s = Math.min(Math.abs(s1), Math.abs(s2));
				}// just in case
				txy = s;
			} else {
				// something weird for point inside cylinder case=> not hit!!
			}

		}

		return txy;
	}

	// ==**********************************************************************************************

	// here, Organ means organism!!!@@@@@@@@@@@@@@@@@@@@@@@ONLY In MAMO USED
	/**
	 * Set breast radius in cm. Used for mammography simulation.
	 * @param organDim organDim
	 */
	public void setMaximumOrganDimension(double organDim) {
		if (mode == MAMO_INDEX)
			BREAST_RADIUS = organDim;
	}

	// here, Organ means organism!!!@@@@@@@@@@@@@@@@@@@@@@@ONLY In MAMO USED
	/**
	 * Set breast thickness in cm. Used for mammography simulation.
	 * @param organTh organTh
	 */
	public void setOrganThickness(double organTh) {
		if (mode == MAMO_INDEX)
			BREAST_THICKNESS = organTh;
	}

	// here, Organ means organism!!!@@@@@@@@@@@@@@@@@@@@@@@ONLY In MAMO USED
	/**
	 * Used for mammography simulation.
	 * @return breast thickness in cm.
	 */
	public double getOrganThickness() {
		double result = 0.0;
		if (mode == MAMO_INDEX)
			result = BREAST_THICKNESS;
		return result;
	}

	// used in computation of theta max for opening tube
	// here, Organ means organism!!!
	/**
	 * Used for mammography simulation.
	 * @return breast radius in cm. 
	 */
	public double getMaximumOrganDimension()// @@@@@@@@@@@@@@@@@@@@@@ONLY In
											// MAMO USED
	{
		double result = 0.0;
		if (mode == MAMO_INDEX)
			result = BREAST_RADIUS;
		return result;
	}

	/**
	 * Used for mammography simulation.
	 * @return breast mass in g. 
	 */
	@Deprecated
	public double getOrganMass()// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ONLY In MAMO
								// USED
	{
		double result = 0.0;
		if (mode == MAMO_INDEX)
			result = BREAST_MASS;
		return result;
	}

	// --The entrance area of XRay Field->general x*y--@@@@@@@@@@@@@@ONLY In
	// MAMO USED
	/**
	 * Return the entrance traversal area, i.e. the breast traversal area in cm^2. 
	 * Used for mammography simulation.
	 * @return the result
	 */
	public double getEntranceTransversalArea() {
		double result = 0.0;
		if (mode == MAMO_INDEX)
			result = Math.PI * BREAST_RADIUS * BREAST_RADIUS;
		return result;
	}

	// -------------------------------------------------------------------------------------
	/**
	 * Return the MIRD5 phantom vertex z-coordinate in cm (top of the head).
	 * @return the result
	 */
	public double getVertexCoordonate() {
		double result = 0.0;
		if (mode == MAMO_INDEX)
			result = 0.0;// BREAST_THICKNESS;
		else // ADULT_INDEX=6 and resources =5
		{
			double[][] headthicknesstab = (double[][]) resources
					.getObject("thickness.head");
			result = headthicknesstab[mode - 1][2];
		}
		return result;
	}

	/**
	 * Given the projection index and z-coordinate, return the thickness of body-part in cm. 
	 * The projection code is: none for mammography (use any int here); 0 for AP, 
	 * 2 for PA, 1 for LLAT and 3 for RLAT.
	 * @param proj proj
	 * @param z z
	 * @return the result.
	 */
	public double getOrganProjThickness(int proj, double z) {
		double result = 0.0;
		if (mode == MAMO_INDEX)
			result = BREAST_THICKNESS;
		else // ADULT_INDEX=6 and resources =5
		{
			double[][] headthicknesstab = (double[][]) resources
					.getObject("thickness.head");
			double[][] trunkthicknesstab = (double[][]) resources
					.getObject("thickness.trunk");
			if (z >= trunkthicknesstab[mode - 1][2]
					&& z <= headthicknesstab[mode - 1][2])// head
			{
				if (proj == 0 || proj == 2)// AP,PA
				{
					result = headthicknesstab[mode - 1][0];
				} else// LL,RL-->1,3
				{
					result = headthicknesstab[mode - 1][1];
				}
			} else// trunk or legs
			{
				if (proj == 0 || proj == 2)// AP,PA
				{
					result = trunkthicknesstab[mode - 1][0];
				} else// LL,RL-->1,3
				{
					result = trunkthicknesstab[mode - 1][1];
				}
			}
			// legs average~trunk
		}

		return result;
	}

	// ---------------------------------------------------------------------------------------------
	/**
	 * Given the projection index and z-coordinate, return the width of body-part in cm. 
	 * The projection code is: none for mammography (use any int here); 0 for AP, 
	 * 2 for PA, 1 for LLAT and 3 for RLAT.
	 * @param proj proj
	 * @param z z
	 * @return the result.
	 */
	public double getOrganProjWidth(int proj, double z) {
		double result = 0.0;
		if (mode == MAMO_INDEX)
			result = 2 * BREAST_RADIUS;
		else // ADULT_INDEX=6 and resources =5
		{
			double[][] headthicknesstab = (double[][]) resources
					.getObject("thickness.head");
			double[][] trunkthicknesstab = (double[][]) resources
					.getObject("thickness.trunk");
			if (z >= trunkthicknesstab[mode - 1][2]
					&& z <= headthicknesstab[mode - 1][2])// head
			{
				if (proj == 0 || proj == 2)// AP,PA
				{
					result = headthicknesstab[mode - 1][1];
				} else// LL,RL-->1,3
				{
					result = headthicknesstab[mode - 1][0];
				}
			} else// trunk or legs
			{
				if (proj == 0 || proj == 2)// AP,PA
				{
					result = trunkthicknesstab[mode - 1][1];
				} else// LL,RL-->1,3
				{
					result = trunkthicknesstab[mode - 1][0];
				}
			}
			// legs average~trunk
		}

		return result;
	}

	// ---------------------------------------------------------------------------------------------
	// proj 0-360 degres!!xy represents sqrt(x*x+y*y);
	/**
	 * Return the thickness of bodypart in cm based on projection angle (0-360) and the 
	 * z coordinates in cm.
	 * @param proj proj in degrees (0-360)
	 * @param z z
	 * @return the result
	 */
	public double getGeneralOrganProjThickness(double proj, double z) {
		double result = 0.0;
		// -----------------
		// lambda=proj in radians=>AP=PI/2;LL=0.0;PA=3*PI/2;RL=PI
		// -------------------
		double lambda = proj * Math.PI / 180;
		double t = 0.0;// angle from parametric eq of elipsis
		double[][] tab = (double[][]) resources.getObject("trunk.constant");
		double[] abc = tab[mode - 1];
		double at = abc[0];
		double bt = abc[1];
		//double ct = abc[2];
		double[][] tabh = (double[][]) resources.getObject("head.constant");
		double[] abch = tabh[mode - 1];
		//double rh = abch[0];
		double ah = abch[1];
		double bh = abch[2];
		//double ch0 = abch[3];
		//double ch1 = abch[4];
		//double ch2 = abch[5];
		// we are going for maximum in order to further fitting!!!!
		double xv = 0.0;// *Math.sin(lambda);//virtual
		double yv = 0.0;// *Math.cos(lambda);//virtual

		if (mode == MAMO_INDEX)
			result = BREAST_THICKNESS;
		else // ADULT_INDEX=6 and resources =5
		{
			double[][] headthicknesstab = (double[][]) resources
					.getObject("thickness.head");
			double[][] trunkthicknesstab = (double[][]) resources
					.getObject("thickness.trunk");
			if (z >= trunkthicknesstab[mode - 1][2]
					&& z <= headthicknesstab[mode - 1][2])// head
			{
				if (lambda != Math.PI / 2 && lambda != 3 * Math.PI / 2)
					t = Math.atan(ah * Math.tan(lambda) / bh);
				else
					t = Math.PI / 2;
				xv = ah * Math.cos(t);
				yv = bh * Math.sin(t);

				result = 2 * Math.sqrt(xv * xv + yv * yv);
			} else// trunk or legs
			{

				if (lambda != Math.PI / 2 && lambda != 3 * Math.PI / 2)
					t = Math.atan(at * Math.tan(lambda) / bt);
				else
					t = Math.PI / 2;
				// System.out.println(t);
				xv = at * Math.cos(t);
				yv = bt * Math.sin(t);

				result = 2 * Math.sqrt(xv * xv + yv * yv);
			}
			// legs average~trunk
		}

		return result;
	}

	// proj 0-360 degres!!xy represents sqrt(x*x+y*y);
	/**
	 * Return the width of bodypart in cm based on projection angle (0-360) and the 
	 * z coordinates in cm.
	 * @param proj proj in degrees (0-360)
	 * @param z z
	 * @return the result
	 */
	public double getGeneralOrganProjWidth(double proj, double z) {
		double result = 0.0;
		// -----------------
		// lambda=proj in radians=>AP=PI/2;LL=0.0;PA=3*PI/2;RL=PI
		// -------------------
		double lambda = proj * Math.PI / 180;
		double t = 0.0;// angle from parametric eq of elipsis
		double[][] tab = (double[][]) resources.getObject("trunk.constant");
		double[] abc = tab[mode - 1];
		double at = abc[0];
		double bt = abc[1];
		//double ct = abc[2];
		double[][] tabh = (double[][]) resources.getObject("head.constant");
		double[] abch = tabh[mode - 1];
		//double rh = abch[0];
		double ah = abch[1];
		double bh = abch[2];
		//double ch0 = abch[3];
		//double ch1 = abch[4];
		//double ch2 = abch[5];
		// we are going for maximum in order to further fitting!!!!
		double xv = 0.0;// *Math.sin(lambda);//virtual
		double yv = 0.0;// *Math.cos(lambda);//virtual

		if (mode == MAMO_INDEX)
			result = 2 * BREAST_RADIUS;
		else // ADULT_INDEX=6 and resources =5
		{
			double[][] headthicknesstab = (double[][]) resources
					.getObject("thickness.head");
			double[][] trunkthicknesstab = (double[][]) resources
					.getObject("thickness.trunk");
			if (z >= trunkthicknesstab[mode - 1][2]
					&& z <= headthicknesstab[mode - 1][2])// head
			{
				if (lambda != Math.PI / 2 && lambda != 3 * Math.PI / 2)
					t = Math.PI / 2 + Math.atan(ah * Math.tan(lambda) / bh);
				else
					t = 0.0;
				xv = ah * Math.cos(t);
				yv = bh * Math.sin(t);

				result = 2 * Math.sqrt(xv * xv + yv * yv);
			} else// trunk or legs
			{

				if (lambda != Math.PI / 2 && lambda != 3 * Math.PI / 2)
					t = Math.PI / 2 + Math.atan(at * Math.tan(lambda) / bt);
				else
					t = 0.0;
				// System.out.println(t);
				xv = at * Math.cos(t);
				yv = bt * Math.sin(t);

				result = 2 * Math.sqrt(xv * xv + yv * yv);
			}
			// legs average~trunk
		}

		return result;
	}

	/**
	 * Get the MIRD5 organ mass. For mammography the result is a dummy value. 
	 * @param indx indx of the organ
	 * @return the result
	 */
	public double getOrganMass(int indx)// ONLY In MAMO USED
	{
		double result = 0.0;
		double[][] om = (double[][]) resources.getObject("organs.mass");
		if (mode == MAMO_INDEX)
			result = BREAST_MASS;
		else
			result = om[mode - 1][indx];
		return result;
	}

	// function to test if the foton reach an organ or not
	// for mammography the main coordonate system is at base of cylinder, in
	// center, with
	// z axis upward (heading to focus), phantom "looks" in positive yaxis and
	// the "left hand"
	// is in positive x axis
	// double x, double y, double z-->0,1,2
	/**
	 * Based on input coordinates, this routine find in what organ the particle is.
	 * @param coord coord
	 * @return the organ index.
	 */
	public int inWhatOrgan(double[] coord) {

		if (mode == MAMO_INDEX) {
			/*
			 * if(coord[0]>=-BREAST_RADIUS && coord[0]<=BREAST_RADIUS) {
			 * if(coord[1]>=-BREAST_RADIUS && coord[1]<=BREAST_RADIUS) {
			 * if(coord[2]>=0 && coord[2]<=BREAST_THICKNESS) { return REMAINDER;
			 * } } }
			 */
			if (coord[0] * coord[0] + coord[1] * coord[1] <= BREAST_RADIUS
					* BREAST_RADIUS) {
				if (coord[2] >= 0 && coord[2] <= BREAST_THICKNESS) {
					return REMAINDER;
				}
			}

		}
		// -----------------------------------
		else {
			if (!inBody(coord)) {
				return NO_ORGAN;
			} else {
				if (inBreasts(coord))
					return BREASTS;

				if (inTestes(coord))
					return TESTES;

				if (inThymus(coord))
					return THYMUS;

				if (inThyroid(coord))
					return THYROID;

				if (inSkeleton(coord))
					return SKELETON;

				if (inAdrenals(coord))
					return ADRENALS;

				if (inBrain(coord))
					return BRAIN;

				if (inGallBladder(coord))
					return GALLBLADDER;

				if (inStomach(coord))
					return STOMACH;

				if (inSmallIntestine(coord))
					return SMALLINTESTINE;

				if (inUpperLargeIntestine(coord))
					return UPPERLARGEINTESTINE;

				if (inLowerLargeIntestine(coord))
					return LOWERLARGEINTESTINE;

				if (inHeart(coord))
					return HEART;

				if (inKidneys(coord))
					return KIDNEYS;

				if (inLiver(coord))
					return LIVER;

				if (inLungs(coord))
					return LUNGS;

				if (inOvaries(coord))
					return OVARIES;

				if (inPancreas(coord))
					return PANCREAS;

				if (inSpleen(coord))
					return SPLEEN;

				if (inUrinaryBladder(coord))
					return URINARYBLADDER;

				if (inUterus(coord))
					return UTERUS;

				return REMAINDER;
			}
		}

		return NO_ORGAN;
	}

	// ------------------------------------
	// Fitting x or y to hit the body->INITIAL APROXIMATION
	// ONLY FOR RAD SHOULD BE CALLED!!!-->MULAJ PE FANTOM TO HIT!!!!!
	/**
	 * Return the required initial coordinate making sure the MIRD5 phantom is hit.
	 * @param proj the projection index
	 * @param zinit initial z coordinate
	 * @param xy the initial x or y
	 * @return the required y or x
	 */
	public double getFittedCoord(int proj, double zinit, double xy) {
		double result = 100000.0;// being sure that with this value we will miss
									// the body!!!!!!

		double[][] tab = (double[][]) resources.getObject("trunk.constant");
		double[] abc = tab[mode - 1];
		double at = abc[0];
		double bt = abc[1];
		double ct = abc[2];
		// ///////////////////////////////////////
		// BREASTS---REGION---AP,LL,RL ONLY
		double[][] tabb = (double[][]) resources.getObject("breasts.constant");
		double[] cnst = tabb[mode - 1];
		double a = cnst[0];
		double b = cnst[1];
		double cbrst = cnst[2];
		double x0 = cnst[3];
		double z0 = cnst[4];
		double y0 = -bt * Math.sqrt(1 - Math.pow(x0 / at, 2));
		if (proj == 0)// AP-->here xy is x!!!
		{
			double test = 1 - Math.pow((xy - x0) / a, 2)
					- Math.pow((zinit - z0) / cbrst, 2);// left breast
			if (test >= 0) {
				double igrec = -(y0 + b * Math.sqrt(test)) + 0.001;
				double tst = Math.pow(xy / at, 2) + Math.pow(igrec / bt, 2);
				if (tst > 1)
					return igrec;
			}
			double test1 = 1 - Math.pow((xy + x0) / a, 2)
					- Math.pow((zinit - z0) / cbrst, 2);// right brst
			if (test1 >= 0) {
				double igrec = -(y0 + b * Math.sqrt(test1)) + 0.001;
				double tst = Math.pow(xy / at, 2) + Math.pow(igrec / bt, 2);
				if (tst > 1)
					return igrec;
			}

		}
		if (proj == 1)// LL-->here xy is y!!!->xnegativ.ONLY RIGHT
		{
			double test = 1 - Math.pow((xy - y0) / b, 2)
					- Math.pow((zinit - z0) / cbrst, 2);// left breast
			if (test >= 0) {
				double igrec = -(-x0 + a * Math.sqrt(test)) + 0.001;
				double tst = Math.pow(igrec / at, 2) + Math.pow(xy / bt, 2);
				if (tst > 1)
					return igrec;
			}

		}
		if (proj == 3)// RL-->here xy is y!!!->xnegativ.ONLY LEFT
		{
			double test = 1 - Math.pow((xy - y0) / b, 2)
					- Math.pow((zinit - z0) / cbrst, 2);// left breast
			if (test >= 0) {
				double igrec = (x0 + a * Math.sqrt(test)) - 0.001;
				double tst = Math.pow(igrec / at, 2) + Math.pow(xy / bt, 2);
				if (tst > 1)
					return igrec;
			}

		}
		// System.out.println("tapa");
		// //////////////////////////////////////
		// test TRUNK REGION (INITIAL GO FOR BREASTS!!)
		if (0 <= zinit && zinit <= ct) {
			if (proj == 0)// AP-->here xy is x!!!
			{
				double test = 1 - Math.pow(xy / at, 2);
				if (test >= 0) {
					result = -bt * Math.sqrt(test) + 0.001;// for be sure to hit
					return result;
				}
			}

			if (proj == 2)// PA-->here xy is x!!!
			{
				double test = 1 - Math.pow(xy / at, 2);
				if (test >= 0) {
					result = bt * Math.sqrt(test) - 0.001;// for be sure to hit
					return result;
				}
			}

			if (proj == 1)// LL-->here xy is y!!!
			{
				double test = 1 - Math.pow(xy / bt, 2);
				if (test >= 0) {
					result = -at * Math.sqrt(test) + 0.001;// for be sure to hit
					return result;
				}
			}

			if (proj == 3)// RL-->here xy is y!!!
			{
				double test = 1 - Math.pow(xy / bt, 2);
				if (test >= 0) {
					result = at * Math.sqrt(test) - 0.001;// for be sure to hit
					return result;
				}
			}
		}
		// end trunk
		// test HEAD+NECK REGION
		double[][] tabh = (double[][]) resources.getObject("head.constant");
		double[] abch = tabh[mode - 1];
		double rh = abch[0];
		double ah = abch[1];
		double bh = abch[2];
		double ch0 = abch[3];
		double ch1 = abch[4];
		double ch2 = abch[5];
		// neck region----------------------------------------
		if (ct <= zinit && zinit <= ct + ch0) {
			if (proj == 0 || proj == 1)// AP-->here xy is x!!!
			{
				double test = rh * rh - Math.pow(xy, 2);// System.out.println("neck  "+test);
				if (test >= 0) {
					result = -Math.sqrt(test) + 0.001;// for be sure to hit
					return result;
				}
			}

			if (proj == 2 || proj == 3)// PA-->here xy is x!!!
			{
				double test = rh * rh - Math.pow(xy, 2);
				if (test >= 0) {
					result = Math.sqrt(test) - 0.001;// for be sure to hit
					return result;
				}
			}
		}
		// ---------------------end neck
		// head main part
		if (ct + ch0 <= zinit && zinit <= ct + ch0 + ch1) {
			if (proj == 0)// AP-->here xy is x!!!
			{
				double test = 1 - Math.pow(xy / ah, 2);
				if (test >= 0) {
					result = -bh * Math.sqrt(test) + 0.001;// for be sure to hit
					return result;
				}
			}

			if (proj == 2)// PA-->here xy is x!!!
			{
				double test = 1 - Math.pow(xy / ah, 2);
				if (test >= 0) {
					result = bh * Math.sqrt(test) - 0.001;// for be sure to hit
					return result;
				}
			}

			if (proj == 1)// LL-->here xy is y!!!
			{
				double test = 1 - Math.pow(xy / bh, 2);
				if (test >= 0) {
					result = -ah * Math.sqrt(test) + 0.001;// for be sure to hit
					return result;
				}
			}

			if (proj == 3)// RL-->here xy is y!!!
			{
				double test = 1 - Math.pow(xy / bh, 2);
				if (test >= 0) {
					result = ah * Math.sqrt(test) - 0.001;// for be sure to hit
					return result;
				}
			}
		}
		// end head main part
		// top head
		if (zinit > ct + ch0 + ch1) {
			if (proj == 0)// AP-->here xy is x!!!
			{
				double test = 1 - Math.pow(xy / ah, 2)
						- Math.pow((zinit - ct - ch0 - ch1) / (ch2), 2);
				if (test >= 0) {
					result = -bh * Math.sqrt(test) + 0.001;// for be sure to hit
					return result;
				}
			}

			if (proj == 2)// PA-->here xy is x!!!
			{
				double test = 1 - Math.pow(xy / ah, 2)
						- Math.pow((zinit - ct - ch0 - ch1) / (ch2), 2);
				if (test >= 0) {
					result = bh * Math.sqrt(test) - 0.001;// for be sure to hit
					return result;
				}
			}

			if (proj == 1)// LL-->here xy is y!!!
			{
				double test = 1 - Math.pow(xy / bh, 2)
						- Math.pow((zinit - ct - ch0 - ch1) / (ch2), 2);
				if (test >= 0) {
					result = -ah * Math.sqrt(test) + 0.001;// for be sure to hit
					return result;
				}
			}

			if (proj == 3)// RL-->here xy is y!!!
			{
				double test = 1 - Math.pow(xy / bh, 2)
						- Math.pow((zinit - ct - ch0 - ch1) / (ch2), 2);
				if (test >= 0) {
					result = ah * Math.sqrt(test) - 0.001;// for be sure to hit
					return result;
				}
			}
		}

		// end top head
		// end HEAD+NECK
		// MALEGENITALA---AP ONLY
		double[][] tabl = (double[][]) resources.getObject("legs.constant");
		double[] abcl = tabl[mode - 1];
		double cl = abcl[0];
		double clprime = abcl[1];
		// ===============================================================
		double[][] tabs = (double[][]) resources.getObject("skin.constant");
		double[] abcs = tabs[mode - 1];
		double s = abcs[0];
		double[][] tabtst = (double[][]) resources.getObject("testes.constant");
		double[] abctst = tabtst[mode - 1];
		double c = abctst[2];
		double r = 0.5 * at * (1 + zinit / clprime);
		double z1 = -(2 * c + s);
		if (z1 <= zinit && zinit <= 0) {
			if (proj == 0)// AP-->here xy is x!!!
			{
				if (xy <= r && xy >= -r) {
					result = -r + 0.001;
					return result;
				}
			}
		}
		// end MALEGENITALA
		// LEGS REGION
		double legsconst = at + at * zinit / cl;
		if (-cl <= zinit && zinit <= 0) {
			if (proj == 0)// AP-->here xy is x!!!
			{
				double test = xy * legsconst - Math.pow(xy, 2);// left leg
				if (test >= 0) {
					result = -Math.sqrt(test) + 0.001;// for be sure to hit
					return result;
				}
				double test1 = -xy * legsconst - Math.pow(xy, 2);// right leg
				if (test1 >= 0) {
					result = -Math.sqrt(test1) + 0.001;// for be sure to hit
					return result;
				}

			}

			if (proj == 2)// PA-->here xy is x!!!
			{
				double test = xy * legsconst - Math.pow(xy, 2);// left leg
				if (test >= 0) {
					result = Math.sqrt(test) - 0.001;// for be sure to hit
					return result;
				}
				double test1 = -xy * legsconst - Math.pow(xy, 2);// right leg
				if (test1 >= 0) {
					result = Math.sqrt(test1) + 0.001;// for be sure to hit
					return result;
				}

			}
			// ------------------------
			if (proj == 1)// LL-->here xy is y!!!ONLY RIGHT
			{
				double test = legsconst * legsconst - 4 * Math.pow(xy, 2);
				if (test >= 0) {
					result = 0.001 + (-legsconst - Math.sqrt(test)) / 2;// for
																		// be
																		// sure
																		// to
																		// hit
					return result;
				}
			}

			if (proj == 3)// RL-->here xy is y!!!ONLY LEFT
			{
				double test = legsconst * legsconst - 4 * Math.pow(xy, 2);
				if (test >= 0) {
					result = -0.001 + (legsconst + Math.sqrt(test)) / 2;// for
																		// be
																		// sure
																		// to
																		// hit
					return result;
				}
			}
		}

		// END LEGS REGION
		// -----------------------END----------------------------------------------------------
		result = 100000.0;// being sure that with this value we will miss the
							// body!!!!!!
		return result;
	}

	// ------------------------------------
	// Fitting x or y to hit the body->INITIAL APROXIMATION
	// ONLY FOR RAD SHOULD BE CALLED!!!-->MULAJ PE FANTOM TO HIT!!!!!
	// proj 0-360 degres!!xy represents sqrt(x*x+y*y);
	/** 
	 * Return the required initial coordinate making sure the MIRD5 phantom is hit.
	 * @param proj the projection angle in degrees (0-360)
	 * @param zinit initial z coordinate
	 * @param xy the initial x or y
	 * @return the required y or x
	 */
	public double getGeneralFittedCoord(double proj, double zinit, double xy) {
		double result = 100000.0;// being sure that with this value we will miss
									// the body!!!!!!

		// -----------------
		// lambda=proj in radians=>AP=PI/2;LL=0.0;PA=3*PI/2;RL=PI
		// -------------------
		double lambda = proj * Math.PI / 180;
		// double xv=xy*Math.sin(lambda);//virtual
		// double yv=xy*Math.cos(lambda);//virtual
		// -------------------------------------------------------
		double xv = xy;
		double yv = xy;
		// -------------------------------------------------------
		double[][] tab = (double[][]) resources.getObject("trunk.constant");
		double[] abc = tab[mode - 1];
		double at = abc[0];
		double bt = abc[1];
		double ct = abc[2];
		// ///////////////////////////////////////
		// BREASTS---REGION---AP,LL,RL ONLY
		double[][] tabb = (double[][]) resources.getObject("breasts.constant");
		double[] cnst = tabb[mode - 1];
		double a = cnst[0];
		double b = cnst[1];
		double cbrst = cnst[2];
		double x0 = cnst[3];
		double z0 = cnst[4];
		double y0 = -bt * Math.sqrt(1 - Math.pow(x0 / at, 2));
		// if ((lambda>=Math.PI/2 && lambda<Math.PI))// || (lambda>=3*Math.PI/2
		// && lambda<2*Math.PI))
		if ((lambda >= Math.PI / 4 && lambda <= 3 * Math.PI / 4))// only AP
																	// HERE!!!!
		// AP or PA in extrem-->here xy is xv!!!
		// only AP HERE!!!!
		{
			double test = 1 - Math.pow((xv - x0) / a, 2)
					- Math.pow((zinit - z0) / cbrst, 2);// left breast
			if (test >= 0) {
				double igrec = 0.0;
				if (lambda >= Math.PI / 2 && lambda < Math.PI)
					igrec = -((y0 + b * Math.sqrt(test)) - 0.001);
				else
					igrec = ((y0 + b * Math.sqrt(test)) - 0.001);
				double tst = Math.pow(xv / at, 2) + Math.pow(igrec / bt, 2);
				if (tst > 1)
					return igrec;
			}
			double test1 = 1 - Math.pow((xv + x0) / a, 2)
					- Math.pow((zinit - z0) / cbrst, 2);// right brst
			if (test1 >= 0) {
				double igrec = 0.0;
				if (lambda >= Math.PI / 2 && lambda < Math.PI)
					igrec = -((y0 + b * Math.sqrt(test1)) - 0.001);
				else
					igrec = ((y0 + b * Math.sqrt(test1)) - 0.001);
				double tst = Math.pow(xv / at, 2) + Math.pow(igrec / bt, 2);
				if (tst > 1)
					return igrec;
			}

		}
		// if ((lambda>=0.0 && lambda<Math.PI/2) || (lambda>=Math.PI &&
		// lambda<3*Math.PI/2))
		if (((lambda >= 0.0 && lambda < Math.PI / 4) || (lambda > 7 * Math.PI / 4 && lambda <= 2 * Math.PI))// LL
				|| (lambda > 3 * Math.PI / 4 && lambda < 5 * Math.PI / 4))// RL
		// LL OR RL AT EXTREEM-->xy is yv
		{
			double test = 1 - Math.pow((yv - y0) / b, 2)
					- Math.pow((zinit - z0) / cbrst, 2);// left breast
			if (test >= 0) {
				double igrec = 0.0;
				if (((lambda >= 0.0 && lambda < Math.PI / 4) || (lambda > 7 * Math.PI / 4 && lambda <= 2 * Math.PI)))
					igrec = -(-x0 + a * Math.sqrt(test) - 0.001);
				else
					igrec = (x0 + a * Math.sqrt(test)) - 0.001;
				double tst = Math.pow(igrec / at, 2) + Math.pow(yv / bt, 2);
				if (tst > 1)
					return igrec;
			}

		}
		// //////////////////////////////////////
		// test TRUNK REGION (INITIAL GO FOR BREASTS!!)
		if (0 <= zinit && zinit <= ct) {
			// if ((lambda>=Math.PI/2 && lambda<Math.PI) || (lambda>=3*Math.PI/2
			// && lambda<2*Math.PI))
			if ((lambda >= Math.PI / 4 && lambda <= 3 * Math.PI / 4)
					|| (lambda >= 5 * Math.PI / 4 && lambda <= 7 * Math.PI / 4))
			// here xy is xv!!!
			{
				double test = 1 - Math.pow(xv / at, 2);
				if (test >= 0) {
					if (lambda >= Math.PI / 4 && lambda <= 3 * Math.PI / 4)// @@@@@@@@@@@@@@@@@
						result = -(bt * Math.sqrt(test) - 0.001);// for be sure
																	// to hit
					else
						result = bt * Math.sqrt(test) - 0.001;// for be sure to
																// hit
					return result;
				}
			}

			// if ((lambda>=0.0 && lambda<Math.PI/2) || (lambda>=Math.PI &&
			// lambda<3*Math.PI/2))
			if (((lambda >= 0.0 && lambda < Math.PI / 4) || (lambda > 7 * Math.PI / 4 && lambda <= 2 * Math.PI))// LL
					|| (lambda > 3 * Math.PI / 4 && lambda < 5 * Math.PI / 4))// RL
			// here xy is yv!!!
			{
				double test = 1 - Math.pow(yv / bt, 2);
				if (test >= 0) {
					if (((lambda >= 0.0 && lambda < Math.PI / 4) || (lambda > 7 * Math.PI / 4 && lambda <= 2 * Math.PI)))
						result = -at * Math.sqrt(test) + 0.001;// for be sure to
																// hit
					else
						result = at * Math.sqrt(test) - 0.001;// for be sure to
																// hit
					return result;
				}
			}
		}
		// end trunk
		// test HEAD+NECK REGION
		double[][] tabh = (double[][]) resources.getObject("head.constant");
		double[] abch = tabh[mode - 1];
		double rh = abch[0];
		double ah = abch[1];
		double bh = abch[2];
		double ch0 = abch[3];
		double ch1 = abch[4];
		double ch2 = abch[5];
		// neck region----------------------------------------
		if (ct <= zinit && zinit <= ct + ch0) {
			// if ((lambda>=Math.PI/2 && lambda<Math.PI) || (lambda>=3*Math.PI/2
			// && lambda<2*Math.PI))
			if ((lambda >= Math.PI / 4 && lambda <= 3 * Math.PI / 4)
					|| (lambda >= 5 * Math.PI / 4 && lambda <= 7 * Math.PI / 4))
			// AP-->here xy is x!!!
			{
				double test = rh * rh - Math.pow(xv, 2);// System.out.println("neck  "+test);
				if (test >= 0) {
					if (lambda >= Math.PI / 4 && lambda <= 3 * Math.PI / 4)
						result = -Math.sqrt(test) + 0.001;// for be sure to hit
					else
						result = Math.sqrt(test) - 0.001;// for be sure to hit
					return result;
				}
			}

			// if ((lambda>=0.0 && lambda<Math.PI/2) || (lambda>=Math.PI &&
			// lambda<3*Math.PI/2))
			if (((lambda >= 0.0 && lambda < Math.PI / 4) || (lambda > 7 * Math.PI / 4 && lambda <= 2 * Math.PI))// LL
					|| (lambda > 3 * Math.PI / 4 && lambda < 5 * Math.PI / 4))// RL
			{
				double test = rh * rh - Math.pow(yv, 2);
				if (test >= 0) {
					if (((lambda >= 0.0 && lambda < Math.PI / 4) || (lambda > 7 * Math.PI / 4 && lambda <= 2 * Math.PI)))
						result = -Math.sqrt(test) + 0.001;// for be sure to hit
					else
						result = Math.sqrt(test) - 0.001;// for be sure to hit
					return result;
				}
			}
		}
		// ---------------------end neck
		// head main part
		if (ct + ch0 <= zinit && zinit <= ct + ch0 + ch1) {
			// if ((lambda>=Math.PI/2 && lambda<Math.PI) || (lambda>=3*Math.PI/2
			// && lambda<2*Math.PI))
			if ((lambda >= Math.PI / 4 && lambda <= 3 * Math.PI / 4)
					|| (lambda >= 5 * Math.PI / 4 && lambda <= 7 * Math.PI / 4))
			// AP-->here xy is x!!!
			{
				double test = 1 - Math.pow(xv / ah, 2);
				if (test >= 0) {
					if (lambda >= Math.PI / 4 && lambda <= 3 * Math.PI / 4)
						result = -bh * Math.sqrt(test) + 0.001;// for be sure to
																// hit
					else
						result = bh * Math.sqrt(test) - 0.001;// for be sure to
																// hit
					return result;
				}
			}

			// if ((lambda>=0.0 && lambda<Math.PI/2) || (lambda>=Math.PI &&
			// lambda<3*Math.PI/2))
			if (((lambda >= 0.0 && lambda < Math.PI / 4) || (lambda > 7 * Math.PI / 4 && lambda <= 2 * Math.PI))// LL
					|| (lambda > 3 * Math.PI / 4 && lambda < 5 * Math.PI / 4))// RL
			// LL-->here xy is y!!!
			{
				double test = 1 - Math.pow(yv / bh, 2);
				if (test >= 0) {
					if (((lambda >= 0.0 && lambda < Math.PI / 4) || (lambda > 7 * Math.PI / 4 && lambda <= 2 * Math.PI)))
						result = -ah * Math.sqrt(test) + 0.001;// for be sure to
																// hit
					else
						result = ah * Math.sqrt(test) - 0.001;// for be sure to
																// hit
					return result;
				}
			}

		}
		// end head main part
		// top head
		if (zinit > ct + ch0 + ch1) {
			// if ((lambda>=Math.PI/2 && lambda<Math.PI) || (lambda>=3*Math.PI/2
			// && lambda<2*Math.PI))
			if ((lambda >= Math.PI / 4 && lambda <= 3 * Math.PI / 4)
					|| (lambda >= 5 * Math.PI / 4 && lambda <= 7 * Math.PI / 4))
			// AP-->here xy is x!!!
			{
				double test = 1 - Math.pow(xv / ah, 2)
						- Math.pow((zinit - ct - ch0 - ch1) / (ch2), 2);
				if (test >= 0) {
					if (lambda >= Math.PI / 4 && lambda <= 3 * Math.PI / 4)
						result = -bh * Math.sqrt(test) + 0.001;// for be sure to
																// hit
					else
						result = bh * Math.sqrt(test) - 0.001;// for be sure to
																// hit
					return result;
				}
			}

			if (((lambda >= 0.0 && lambda < Math.PI / 4) || (lambda > 7 * Math.PI / 4 && lambda <= 2 * Math.PI))// LL
					|| (lambda > 3 * Math.PI / 4 && lambda < 5 * Math.PI / 4))// RL
			// LL-->here xy is y!!!
			{
				double test = 1 - Math.pow(yv / bh, 2)
						- Math.pow((zinit - ct - ch0 - ch1) / (ch2), 2);
				if (test >= 0) {
					if (((lambda >= 0.0 && lambda < Math.PI / 4) || (lambda > 7 * Math.PI / 4 && lambda <= 2 * Math.PI)))
						result = -ah * Math.sqrt(test) + 0.001;// for be sure to
																// hit
					else
						result = ah * Math.sqrt(test) - 0.001;// for be sure to
																// hit
					return result;
				}
			}

		}

		// end top head
		// end HEAD+NECK
		// MALEGENITALA---AP ONLY
		double[][] tabl = (double[][]) resources.getObject("legs.constant");
		double[] abcl = tabl[mode - 1];
		double cl = abcl[0];
		double clprime = abcl[1];
		// ===============================================================
		double[][] tabs = (double[][]) resources.getObject("skin.constant");
		double[] abcs = tabs[mode - 1];
		double s = abcs[0];
		double[][] tabtst = (double[][]) resources.getObject("testes.constant");
		double[] abctst = tabtst[mode - 1];
		double c = abctst[2];
		double r = 0.5 * at * (1 + zinit / clprime);
		double z1 = -(2 * c + s);
		if (z1 <= zinit && zinit <= 0) {
			// if ((lambda>=Math.PI/2 && lambda<Math.PI))// ||
			// (lambda>=3*Math.PI/2 && lambda<2*Math.PI))
			if ((lambda >= Math.PI / 4 && lambda <= 3 * Math.PI / 4))// only AP
																		// HERE!!!!
			// AP-->here xy is x!!!
			{
				if (xv <= r && xv >= -r) {
					result = -r + 0.001;
					return result;
				}
			}
		}
		// end MALEGENITALA
		// LEGS REGION
		double legsconst = at + at * zinit / cl;
		if (-cl <= zinit && zinit <= 0) {
			// if ((lambda>=Math.PI/2 && lambda<Math.PI) || (lambda>=3*Math.PI/2
			// && lambda<2*Math.PI))
			if ((lambda >= Math.PI / 4 && lambda <= 3 * Math.PI / 4)
					|| (lambda >= 5 * Math.PI / 4 && lambda <= 7 * Math.PI / 4))
			// AP-->here xy is x!!!
			{
				double test = xv * legsconst - Math.pow(xv, 2);// left leg
				if (test >= 0) {
					if (lambda >= Math.PI / 4 && lambda <= 3 * Math.PI / 4)
						result = -Math.sqrt(test) + 0.001;// for be sure to hit
					else
						result = Math.sqrt(test) - 0.001;// for be sure to hit
					return result;
				}
				double test1 = -xv * legsconst - Math.pow(xv, 2);// right leg
				if (test1 >= 0) {
					if (lambda >= Math.PI / 2 && lambda < Math.PI)
						result = -Math.sqrt(test1) + 0.001;// for be sure to hit
					else
						result = Math.sqrt(test1) + 0.001;// for be sure to hit
					return result;
				}

			}
			// ------------------------
			// if ((lambda>=0.0 && lambda<Math.PI/2) || (lambda>=Math.PI &&
			// lambda<3*Math.PI/2))
			if (((lambda >= 0.0 && lambda < Math.PI / 4) || (lambda > 7 * Math.PI / 4 && lambda <= 2 * Math.PI))// LL
					|| (lambda > 3 * Math.PI / 4 && lambda < 5 * Math.PI / 4))// RL
			// LL-->here xy is y!!!ONLY RIGHT
			{
				double test = legsconst * legsconst - 4 * Math.pow(yv, 2);
				if (test >= 0) {
					if (((lambda >= 0.0 && lambda < Math.PI / 4) || (lambda > 7 * Math.PI / 4 && lambda <= 2 * Math.PI)))
						result = 0.001 + (-legsconst - Math.sqrt(test)) / 2;// for
																			// be
																			// sure
																			// to
																			// hit
					else
						result = -0.001 + (legsconst + Math.sqrt(test)) / 2;// for
																			// be
																			// sure
																			// to
																			// hit
					return result;
				}
			}

		}

		// END LEGS REGION
		// -----------------------END----------------------------------------------------------
		result = 100000.0;// being sure that with this value we will miss the
							// body!!!!!!
		return result;
	}

	// return true only if is in organism
	/**
	 * Based on input coordinates, find if particle is in this bodypart
	 * @param coord coord
	 * @return true or false
	 */
	public boolean inBody(double[] coord) {
		if (testBreastsRegion(coord)) // useful in LL,RL proj
			return true;

		if (testTrunk(coord))
			return true;

		if (testHead(coord))
			return true;

		if (testMaleGenitala(coord))
			return true;

		if (testLegs(coord))
			return true;

		return false;
	}

	// -TEST TRUNK--REGION
	/**
	 * Based on input coordinates, find if particle is in this bodypart
	 * @param coord coord
	 * @return true or false
	 */
	public boolean testTrunk(double[] coord) {
		boolean b = false;
		double[][] tab = (double[][]) resources.getObject("trunk.constant");
		double[] abc = tab[mode - 1];
		if (0 <= coord[2] && coord[2] <= abc[2])
			b = true;
		if (b) {
			double d = Math.pow(coord[0] / abc[0], 2)
					+ Math.pow(coord[1] / abc[1], 2);
			if (d <= 1)
				return true;
		}
		// else is false
		return false;
	}

	// -TEST HEAD--REGION
	/**
	 * Based on input coordinates, find if particle is in this bodypart
	 * @param coord coord
	 * @return true or false
	 */
	public boolean testHead(double[] coord) {
		boolean bneck = false;
		boolean bhead = false;
		double[][] tab = (double[][]) resources.getObject("trunk.constant");
		double[] abc = tab[mode - 1];
		double ct = abc[2];
		double[][] tabh = (double[][]) resources.getObject("head.constant");
		double[] cnst = tabh[mode - 1];
		double ch0 = cnst[3];
		double rh = cnst[0];
		double ah = cnst[1];
		double bh = cnst[2];
		double ch1 = cnst[4];
		double ch2 = cnst[5];
		// neck region
		if (ct <= coord[2] && coord[2] <= ct + ch0)
			bneck = true;
		if (bneck) {
			double d = Math.pow(coord[0], 2) + Math.pow(coord[1], 2);// System.out.println("neck  "+d);
			if (d <= rh * rh)
				bneck = true;
			else
				bneck = false;
		}
		if (bneck)
			return true;
		// end neck
		// head region
		if (ct + ch0 <= coord[2] && coord[2] <= ct + ch0 + ch1)
			bhead = true;
		if (bhead) {
			double d = Math.pow(coord[0] / ah, 2) + Math.pow(coord[1] / bh, 2);
			if (d <= 1)
				bhead = true;
			else
				bhead = false;
		}
		if (bhead)
			return true;
		// end part1
		if (coord[2] > ct + ch0 + ch1)
			bhead = true;
		if (bhead) {
			double d = Math.pow(coord[0] / ah, 2) + Math.pow(coord[1] / bh, 2)
					+ Math.pow((coord[2] - ct - ch0 - ch1) / (ch2), 2);
			if (d <= 1)
				bhead = true;
			else
				bhead = false;
		}
		if (bhead)
			return true;

		return false;
	}

	// -TEST LEGS--REGION
	/**
	 * Based on input coordinates, find if particle is in this bodypart
	 * @param coord coord
	 * @return true or false
	 */
	public boolean testLegs(double[] coord) {
		boolean b = false;
		double[][] tab = (double[][]) resources.getObject("trunk.constant");
		double[] abc = tab[mode - 1];
		double at = abc[0];
		double[][] tabl = (double[][]) resources.getObject("legs.constant");
		double[] cnst = tabl[mode - 1];
		double cl = cnst[0];
		if (-cl <= coord[2] && coord[2] <= 0)
			b = true;
		if (b) {
			double d1 = Math.pow(coord[0], 2) + Math.pow(coord[1], 2);
			double d2 = coord[0] * (at + at * coord[2] / cl);
			// left leg
			if (d1 <= d2)
				return true;
			// right leg
			if (d1 <= -d2)
				return true;
		}

		return false;
	}

	// -TEST MALEGENITAKA--REGION
	/**
	 * Based on input coordinates, find if particle is in this bodypart
	 * @param coord coord
	 * @return true or false
	 */
	public boolean testMaleGenitala(double[] coord) {
		boolean b = false;
		double[][] tab = (double[][]) resources.getObject("trunk.constant");
		double[] abc = tab[mode - 1];
		double at = abc[0];
		double[][] tabl = (double[][]) resources.getObject("legs.constant");
		double[] cnst = tabl[mode - 1];
		double clprime = cnst[1];
		double[][] tabs = (double[][]) resources.getObject("skin.constant");
		double[] cnsts = tabs[mode - 1];
		double s = cnsts[0];
		double[][] tabt = (double[][]) resources.getObject("testes.constant");
		double[] cnstt = tabt[mode - 1];
		double c = cnstt[2];
		// ------------------
		double r = 0.5 * at * (1 + coord[2] / clprime);
		double z1 = -(2 * c + s);
		// ------------------------
		if (z1 <= coord[2] && coord[2] <= 0)
			b = true;
		if (b) {
			if (-r <= coord[0] && coord[0] <= r)
				b = true;
			else
				b = false;
			if (b) {
				if (-r <= coord[1] && coord[1] <= 0)
					b = true;
				else
					b = false;
				if (b) {
					double d = Math.pow(coord[0] + r, 2)
							+ Math.pow(coord[1], 2);
					if (d >= r * r)
						b = true;
					else
						b = false;
					if (b) {
						double d1 = Math.pow(coord[0] - r, 2)
								+ Math.pow(coord[1], 2);
						if (d1 >= r * r)
							return true;
					}
				}

			}
		}

		return false;
	}

	// -TEST BREASTS--REGION--WITH SKIN (WITHOUT -S at numitor!!!!)
	/**
	 * Based on input coordinates, find if particle is in this bodypart
	 * @param coord coord
	 * @return true or false
	 */
	public boolean testBreastsRegion(double[] coord) {
		boolean bool = false;
		double[][] tab = (double[][]) resources.getObject("trunk.constant");
		double[] abc = tab[mode - 1];
		double at = abc[0];
		double bt = abc[1];
		double[][] tabb = (double[][]) resources.getObject("breasts.constant");
		double[] cnst = tabb[mode - 1];
		double a = cnst[0];
		double b = cnst[1];
		double c = cnst[2];
		double x0 = cnst[3];
		double z0 = cnst[4];
		double y0 = -bt * Math.sqrt(1 - Math.pow(x0 / at, 2));
		// -------------
		double d = Math.pow(coord[0] / at, 2) + Math.pow(coord[1] / bt, 2);
		if (d > 1)
			bool = true;
		if (bool) {
			double d1 = Math.pow((coord[0] - x0) / a, 2)
					+ Math.pow((coord[1] - y0) / b, 2)
					+ Math.pow((coord[2] - z0) / c, 2);
			double d2 = Math.pow((coord[0] + x0) / a, 2)
					+ Math.pow((coord[1] - y0) / b, 2)
					+ Math.pow((coord[2] - z0) / c, 2);
			// left breast
			if (d1 <= 1)
				return true;
			// right breast
			if (d2 <= 1)
				return true;
		}

		return false;
	}

	// ------------------END REGIONS
	// ---------------ORGANS------------------------------------------------------------------------------
	// ----------SKELETON PARTIAL
	/**
	 * Based on input coordinates, find if particle is in this bodypart
	 * @param coord coord
	 * @return true or false
	 */
	public boolean inLegBones(double[] coord) {
		double[][] tab = (double[][]) resources.getObject("trunk.constant");
		double[] abc = tab[mode - 1];
		double at = abc[0];
		double[][] tabl = (double[][]) resources.getObject("legs.constant");
		double[] cnst = tabl[mode - 1];
		double cl = cnst[0];
		double clprime = cnst[1];
		double[][] tabs = (double[][]) resources.getObject("skin.constant");
		double[] cnss = tabs[mode - 1];
		double s = cnss[0];
		// ------------------------
		double k = at * (1 - (clprime - cl) / clprime) / 2;
		double r1 = 0.175 * at;
		double r2 = at * ((clprime - cl) / clprime) / 4;
		if (coord[2] >= -(cl - s) && coord[2] <= 0) {
			// left
			double dd = Math.pow(coord[0]
					- ((at / 2) + k * coord[2] / (cl - s)), 2)
					+ Math.pow(coord[1], 2);
			double ddd = Math.pow(r1 + coord[2] * ((r1 - r2) / (cl - s)), 2);
			if (dd <= ddd)
				return true;
			// right
			dd = Math.pow(coord[0] + ((at / 2) + k * coord[2] / (cl - s)), 2)
					+ Math.pow(coord[1], 2);
			ddd = Math.pow(r1 + coord[2] * ((r1 - r2) / (cl - s)), 2);
			if (dd <= ddd)
				return true;
		}

		return false;
	}

	/**
	 * Based on input coordinates, find if particle is in this bodypart
	 * @param coord coord
	 * @return true or false
	 */
	public boolean inArmBones(double[] coord) {
		double[][] tabl = (double[][]) resources.getObject("armBones.constant");
		double[] cnst = tabl[mode - 1];
		double a = cnst[0];
		double b = cnst[1];
		double x0 = cnst[2];
		double z2 = cnst[3];

		if (coord[2] >= 0 && coord[2] <= z2) {
			// /left
			double dd = Math.pow((coord[0] - x0 + (coord[2] - z2) * a
					/ (2 * z2))
					/ a, 2)
					+ Math.pow(coord[1] / b, 2);
			double ddd = Math.pow((2 * z2 + coord[2] - z2) / (2 * z2), 2);
			if (dd <= ddd)
				return true;
			// right
			dd = Math.pow((coord[0] + x0 + (coord[2] - z2) * a / (2 * z2)) / a,
					2) + Math.pow(coord[1] / b, 2);
			ddd = Math.pow((2 * z2 + coord[2] - z2) / (2 * z2), 2);
			if (dd <= ddd)
				return true;

		}

		return false;
	}

	/**
	 * Based on input coordinates, find if particle is in this bodypart
	 * @param coord coord
	 * @return true or false
	 */
	public boolean inPelvis(double[] coord) {
		double[][] tabl = (double[][]) resources.getObject("pelvis.constant");
		double[] cnst = tabl[mode - 1];
		double a1 = cnst[0];
		double b1 = cnst[1];
		double a2 = cnst[2];
		double b2 = cnst[3];
		double y01 = cnst[4];
		double y02 = cnst[5];
		double y1 = cnst[6];
		double z1 = cnst[7];
		double z2 = cnst[8];

		if ((coord[2] >= 0 && coord[2] <= z2) && coord[1] >= y02) {
			if (coord[2] <= z1) {
				if (coord[1] <= y1) {
					double dd = Math.pow((coord[0]) / a2, 2)
							+ Math.pow((coord[1] - y02) / b2, 2);
					if (dd <= 1) {
						double ddd = Math.pow((coord[0]) / a1, 2)
								+ Math.pow((coord[1] - y01) / b1, 2);
						if (ddd >= 1)
							return true;
					}
				}
			} else {
				double dd = Math.pow((coord[0]) / a2, 2)
						+ Math.pow((coord[1] - y02) / b2, 2);
				if (dd <= 1) {
					double ddd = Math.pow((coord[0]) / a1, 2)
							+ Math.pow((coord[1] - y01) / b1, 2);
					if (ddd >= 1)
						return true;
				}
			}
		}
		return false;
	}

	/**
	 * Based on input coordinates, find if particle is in this bodypart
	 * @param coord coord
	 * @return true or false
	 */
	public boolean inSpine(double[] coord) {
		double[][] tabl = (double[][]) resources.getObject("spine.constant");
		double[] cnst = tabl[mode - 1];
		double a = cnst[0];
		double b = cnst[1];
		double y0 = cnst[2];
		double z1 = cnst[3];
		//double z2 = cnst[4];
		//double z3 = cnst[5];
		double z4 = cnst[6];

		if (coord[2] >= z1 && coord[2] <= z4) {
			double dd = Math.pow((coord[0]) / a, 2)
					+ Math.pow((coord[1] - y0) / b, 2);
			if (dd <= 1)
				return true;
		}
		return false;
	}

	/**
	 * Based on input coordinates, find if particle is in this bodypart
	 * @param coord coord
	 * @return true or false
	 */
	public boolean inSkull(double[] coord) {
		double[][] tab = (double[][]) resources.getObject("brain.constant");
		double[] cnst = tab[mode - 1];
		double[][] tabc = (double[][]) resources.getObject("cranium.constant");
		double[] cnstc = tabc[mode - 1];
		double[][] tabf = (double[][]) resources.getObject("facial.constant");
		double[] cnstf = tabf[mode - 1];
		double[][] tabt = (double[][]) resources.getObject("trunk.constant");
		double[] cnstt = tabt[mode - 1];
		double[][] tabh = (double[][]) resources.getObject("head.constant");
		double[] cnsth = tabh[mode - 1];
		double ct = cnstt[2];
		double ch1 = cnsth[4];
		// ------------------------
		double d = cnstc[0];
		double a = cnst[0];
		double b = cnst[1];
		double c = cnst[2];
		// cranium
		double dd = Math.pow((coord[0]) / a, 2) + Math.pow((coord[1]) / b, 2)
				+ Math.pow((coord[2] - ct - ch1) / c, 2);
		double ddd = Math.pow((coord[0]) / (a + d), 2)
				+ Math.pow((coord[1]) / (b + d), 2)
				+ Math.pow((coord[2] - ct - ch1) / (c + d), 2);
		if (dd >= 1 && ddd <= 1)
			return true;
		// -----------
		double a2 = a + d;
		double b2 = b + d;
		double c2 = c + d;
		double a1 = cnstf[0];
		double b1 = cnstf[1];
		d = cnstf[2];
		double z1 = cnstf[3];
		double z5 = cnstf[4];
		if (coord[1] <= 0) {
			if (coord[2] >= z1 + ct && coord[2] <= z5 + ct) {
				double dd0 = Math.pow((coord[0]) / a1, 2)
						+ Math.pow((coord[1]) / b1, 2);
				if (dd0 <= 1) {
					double dd1 = Math.pow((coord[0]) / (a1 - d), 2)
							+ Math.pow((coord[1]) / (b1 - d), 2);
					if (dd1 >= 1) {
						double dd2 = Math.pow((coord[0]) / (a2), 2)
								+ Math.pow((coord[1]) / (b2), 2)
								+ Math.pow((coord[2] - ct - ch1) / (c2), 2);
						if (dd2 > 1)
							return true;
					}
				}

			}
		}

		return false;
	}

	/**
	 * Based on input coordinates, find if particle is in this bodypart
	 * @param coord coord
	 * @return true or false
	 */
	public boolean inRibCage(double[] coord) {
		double[][] tab = (double[][]) resources.getObject("ribcage.constant");
		double[] cnst = tab[mode - 1];
		double a = cnst[0];
		double b = cnst[1];
		double d = cnst[2];
		double z1 = cnst[3];
		double z2 = cnst[4];
		double c = cnst[5];
		// -------------------
		double testint = (coord[2] - z1) / c;
		Double testd = new Double(testint);
		int testintmodif = testd.intValue();// integer part of double
		// /test if is even
		int intst = testintmodif % 2;
		if (intst == 0)// is even
		{
			if (coord[2] >= z1 && coord[2] <= z2) {
				double dd1 = Math.pow((coord[0]) / (a), 2)
						+ Math.pow((coord[1]) / (b), 2);
				if (dd1 <= 1) {
					double dd2 = Math.pow((coord[0]) / (a - d), 2)
							+ Math.pow((coord[1]) / (b - d), 2);
					if (dd2 >= 1)
						return true;
				}
			}
		}

		return false;
	}

	/**
	 * Based on input coordinates, find if particle is in this bodypart
	 * @param coord coord
	 * @return true or false
	 */
	public boolean inClavicles(double[] coord) {
		double[][] tab = (double[][]) resources.getObject("clavicles.constant");
		double[] cnst = tab[mode - 1];
		double y0 = cnst[0];
		double z1 = cnst[1];
		double rr = cnst[2];
		double r = cnst[3];
		double cot1 = cnst[4];
		double cot2 = cnst[5];
		// -------------------
		if (coord[1] < 0) {
			double cot = (y0 - coord[1]) / Math.abs(coord[0]);
			if (cot >= cot2 && cot <= cot1) {
				double dd = Math.pow((coord[2] - z1), 2)
						+ Math.pow(
								(rr - Math.sqrt(coord[0] * coord[0]
										+ (coord[1] - y0) * (coord[1] - y0))),
								2);
				if (dd <= r * r)
					return true;
			}
		}

		return false;
	}

	/**
	 * Based on input coordinates, find if particle is in this bodypart
	 * @param coord coord
	 * @return true or false
	 */
	public boolean inScapulae(double[] coord) {
		double[][] tab = (double[][]) resources.getObject("scapulae.constant");
		double[] cnst = tab[mode - 1];
		double a1 = cnst[0];
		double a2 = cnst[1];
		double b = cnst[2];
		double z1 = cnst[3];
		double z2 = cnst[4];
		double m1 = cnst[5];
		double m2 = cnst[6];
		// ----------------------------------
		if (coord[2] >= z1 && coord[2] <= z2) {
			if (coord[1] > 0) {
				double m = coord[1] / Math.abs(coord[0]);
				if (m > m1 && m < m2) {
					double dd = Math.pow((coord[0]) / (a2), 2)
							+ Math.pow((coord[1]) / (b), 2);
					double dd1 = Math.pow((coord[0]) / (a1), 2)
							+ Math.pow((coord[1]) / (b), 2);
					if (dd <= 1 && dd1 > 1) {
						return true;
					}
				}
			}
		}
		return false;
	}

	// -------------------------------------------------------------------------------------------------
	/**
	 * Based on input coordinates, find if particle is in this bodypart
	 * @param coord coord
	 * @return true or false
	 */
	public boolean inBreasts(double[] coord) {
		boolean bool = false;
		double[][] tabs = (double[][]) resources.getObject("skin.constant");
		double[] cnsts = tabs[mode - 1];
		double s = cnsts[0];
		double[][] tab = (double[][]) resources.getObject("trunk.constant");
		double[] abc = tab[mode - 1];
		double at = abc[0];
		double bt = abc[1];
		double[][] tabb = (double[][]) resources.getObject("breasts.constant");
		double[] cnst = tabb[mode - 1];
		double a = cnst[0];
		double b = cnst[1];
		double c = cnst[2];
		double x0 = cnst[3];
		double z0 = cnst[4];
		double y0 = -bt * Math.sqrt(1 - Math.pow(x0 / at, 2));
		// -------------
		double d = Math.pow(coord[0] / at, 2) + Math.pow(coord[1] / bt, 2);
		if (d > 1)
			bool = true;
		if (bool) {
			double d1 = Math.pow((coord[0] - x0) / (a - s), 2)
					+ Math.pow((coord[1] - y0) / (b - s), 2)
					+ Math.pow((coord[2] - z0) / (c - s), 2);
			double d2 = Math.pow((coord[0] + x0) / (a - s), 2)
					+ Math.pow((coord[1] - y0) / (b - s), 2)
					+ Math.pow((coord[2] - z0) / (c - s), 2);
			// left breast
			if (d1 <= 1)
				return true;
			// right breast
			if (d2 <= 1)
				return true;
		}

		return false;
	}

	/**
	 * Based on input coordinates, find if particle is in this bodypart
	 * @param coord coord
	 * @return true or false
	 */
	public boolean inTestes(double[] coord) {
		double[][] tabt = (double[][]) resources.getObject("testes.constant");
		double[] cnstt = tabt[mode - 1];
		double a = cnstt[0];
		double b = cnstt[1];
		double c = cnstt[2];
		double y0 = cnstt[3];

		double d = Math.pow((coord[0] - a) / a, 2)
				+ Math.pow((coord[1] - y0) / b, 2)
				+ Math.pow((coord[2] + c) / c, 2);// left testes
		double d1 = Math.pow((coord[0] + a) / a, 2)
				+ Math.pow((coord[1] - y0) / b, 2)
				+ Math.pow((coord[2] + c) / c, 2);// right testes
		if (d <= 1)
			return true;
		if (d1 <= 1)
			return true;

		return false;
	}

	/**
	 * Based on input coordinates, find if particle is in this bodypart
	 * @param coord coord
	 * @return true or false
	 */
	public boolean inSkeleton(double[] coord) {
		if (inLegBones(coord))
			return true;
		if (inArmBones(coord))
			return true;
		if (inPelvis(coord))
			return true;
		if (inSpine(coord))
			return true;
		if (inSkull(coord))
			return true;
		if (inRibCage(coord))
			return true;
		if (inClavicles(coord))
			return true;
		if (inScapulae(coord))
			return true;

		return false;
	}

	/**
	 * Based on input coordinates, find if particle is in this bodypart
	 * @param coord coord
	 * @return true or false
	 */
	public boolean inAdrenals(double[] coord) {
		boolean bool = false;
		double[][] tab = (double[][]) resources.getObject("adrenals.constant");
		double[] abc = tab[mode - 1];
		double a = abc[0];
		double b = abc[1];
		double c = abc[2];
		double x0 = abc[3];// +/- for both adrenals
		double x00 = -abc[3];// +/- for both adrenals
		double y0 = abc[4];
		double z0 = abc[5];
		double theta = abc[6];// +/- for both adrenals
		double theta0 = -abc[6];// +/- for both adrenals
		double x1 = (coord[0] - x0) * Math.cos(theta) + (coord[1] - y0)
				* Math.sin(theta);
		double y1 = -(coord[0] - x0) * Math.sin(theta) + (coord[1] - y0)
				* Math.cos(theta);
		double x10 = (coord[0] - x00) * Math.cos(theta0) + (coord[1] - y0)
				* Math.sin(theta0);
		double y10 = -(coord[0] - x00) * Math.sin(theta0) + (coord[1] - y0)
				* Math.cos(theta0);
		double z1 = coord[2] - z0;

		if (z1 >= 0)
			bool = true;
		if (bool) {
			double d1 = Math.pow(x1 / a, 2) + Math.pow(y1 / b, 2)
					+ Math.pow(z1 / c, 2);
			if (d1 <= 1)
				return true;
			// --else we go here
			double d10 = Math.pow(x10 / a, 2) + Math.pow(y10 / b, 2)
					+ Math.pow(z1 / c, 2);
			if (d10 <= 1)
				return true;
		}

		return false;
	}

	/**
	 * Based on input coordinates, find if particle is in this bodypart
	 * @param coord coord
	 * @return true or false
	 */
	public boolean inBrain(double[] coord) {
		double[][] tabt = (double[][]) resources.getObject("brain.constant");
		double[] cnstt = tabt[mode - 1];
		double a = cnstt[0];
		double b = cnstt[1];
		double c = cnstt[2];
		double[][] tab = (double[][]) resources.getObject("trunk.constant");
		double[] cnst = tab[mode - 1];
		double ct = cnst[2];
		double[][] tabh = (double[][]) resources.getObject("head.constant");
		double[] cnsth = tabh[mode - 1];
		double ch1 = cnsth[4];

		double dd = Math.pow((coord[0]) / a, 2) + Math.pow((coord[1]) / b, 2)
				+ Math.pow((coord[2] - ct - ch1) / c, 2);
		if (dd <= 1)
			return true;
		return false;
	}

	/**
	 * Based on input coordinates, find if particle is in this bodypart
	 * @param coord coord
	 * @return true or false
	 */
	public boolean inGallBladder(double[] coord) {
		double[][] tabt = (double[][]) resources
				.getObject("gallbladder.constant");
		double[] cnstt = tabt[mode - 1];
		double alfa1 = cnstt[0];
		double beta1 = cnstt[1];
		double gama1 = cnstt[2];
		double alfa2 = cnstt[3];
		double beta2 = cnstt[4];
		double gama2 = cnstt[5];
		double alfa3 = cnstt[6];
		double beta3 = cnstt[7];
		double gama3 = cnstt[8];
		//double r1 = cnstt[9];
		double r2 = cnstt[10];
		double s = cnstt[11];
		double h = cnstt[12];
		double x0 = cnstt[13];
		double y0 = cnstt[14];
		double z0 = cnstt[15];
		// ------------------------
		double x1 = alfa1 * (coord[0] - x0) + beta1 * (coord[1] - y0) + gama1
				* (coord[2] - z0);
		double y1 = alfa2 * (coord[0] - x0) + beta2 * (coord[1] - y0) + gama2
				* (coord[2] - z0);
		double z1 = alfa3 * (coord[0] - x0) + beta3 * (coord[1] - y0) + gama3
				* (coord[2] - z0);
		// hemispherical part
		if (z1 < 0) {
			double dd = Math.pow((x1), 2) + Math.pow((y1), 2)
					+ Math.pow((z1), 2);
			if (dd <= r2 * r2)
				return true;
		}
		// conical part
		if (z1 >= 0 && z1 <= h) {
			double dd = Math.pow((x1), 2) + Math.pow((y1), 2);
			double ddd = Math.pow((r2 - s * z1), 2);
			if (dd <= ddd)
				return true;
		}

		return false;
	}

	/**
	 * Based on input coordinates, find if particle is in this bodypart
	 * @param coord coord
	 * @return true or false
	 */
	public boolean inStomach(double[] coord) {
		double[][] tabt = (double[][]) resources.getObject("stomach.constant");
		double[] cnstt = tabt[mode - 1];
		double a = cnstt[0];
		double b = cnstt[1];
		double c = cnstt[2];
		//double d = cnstt[3];
		double x0 = cnstt[4];
		double y0 = cnstt[5];
		double z0 = cnstt[6];

		double dd = Math.pow((coord[0] - x0) / a, 2)
				+ Math.pow((coord[1] - y0) / b, 2)
				+ Math.pow((coord[2] - z0) / c, 2);
		if (dd <= 1)
			return true;

		return false;
	}

	/**
	 * Based on input coordinates, find if particle is in this bodypart
	 * @param coord coord
	 * @return true or false
	 */
	public boolean inSmallIntestine(double[] coord) {
		double[][] tabt = (double[][]) resources
				.getObject("smallintestine.constant");
		double[] cnstt = tabt[mode - 1];
		double a = cnstt[0];
		double b = cnstt[1];
		double y0 = cnstt[2];
		double y1 = cnstt[3];
		double y2 = cnstt[4];
		double z1 = cnstt[5];
		double z2 = cnstt[6];
		if (coord[1] >= y1 && coord[1] <= y2) {
			if (coord[2] >= z1 && coord[2] <= z2) {
				double dd = Math.pow((coord[0]) / a, 2)
						+ Math.pow((coord[1] - y0) / b, 2);
				if (dd <= 1)
					return true;
			}
		}

		return false;
	}

	/**
	 * Based on input coordinates, find if particle is in this bodypart
	 * @param coord coord
	 * @return true or false
	 */
	public boolean inUpperLargeIntestine(double[] coord) {
		double[][] tabt = (double[][]) resources
				.getObject("upperlargeintestine.asccolon.constant");
		double[] cnstt = tabt[mode - 1];
		double[][] tab = (double[][]) resources
				.getObject("upperlargeintestine.transcolon.constant");
		double[] cnst = tab[mode - 1];
		double a = cnstt[0];
		double b = cnstt[1];
		//double d = cnstt[2];
		double x0 = cnstt[3];
		double y0 = cnstt[4];
		double z1 = cnstt[5];
		double z2 = cnstt[6];
		// ascending colon TOTAL---------------------
		if (coord[2] >= z1 && coord[2] <= z2) {
			double dd = Math.pow((coord[0] - x0) / a, 2)
					+ Math.pow((coord[1] - y0) / b, 2);
			if (dd <= 1)
				return true;

		}
		// ---------------transverse colon
		b = cnst[0];
		double c = cnst[1];
		//d = cnst[2];
		y0 = cnst[3];
		double z0 = cnst[4];
		double x1 = cnst[5];
		if (coord[0] >= -x1 && coord[0] <= x1) {
			double dd = Math.pow((coord[1] - y0) / b, 2)
					+ Math.pow((coord[2] - z0) / c, 2);
			if (dd <= 1)
				return true;

		}

		return false;
	}

	/**
	 * Based on input coordinates, find if particle is in this bodypart
	 * @param coord coord
	 * @return true or false
	 */
	public boolean inLowerLargeIntestine(double[] coord) {
		double[][] tabt = (double[][]) resources
				.getObject("lowerlargeintestine.desccolon.constant");
		double[] cnstt = tabt[mode - 1];
		double[][] tab = (double[][]) resources
				.getObject("lowerlargeintestine.sigcolon.constant");
		double[] cnst = tab[mode - 1];

		double a = cnstt[0];
		double b = cnstt[1];
		//double d = cnstt[2];
		double x1 = cnstt[3];
		double mx = cnstt[4];
		double my = cnstt[5];
		double z1 = cnstt[6];
		double z2 = cnstt[7];
		// descending colon TOTAL---------------------
		double x0 = x1 + mx * (coord[2] - z2) / (z2 - z1);
		double y0 = my * (z1 - coord[2]) / (z2 - z1);
		if (coord[2] >= z1 && coord[2] <= z2) {
			double dd = Math.pow((coord[0] - x0) / a, 2)
					+ Math.pow((coord[1] - y0) / b, 2);
			if (dd <= 1)
				return true;

		}
		// ---------------sigmoid colon
		a = cnst[0];
		b = cnst[1];
		//d = cnst[2];
		x0 = cnst[3];
		double z0 = cnst[4];
		double r1 = cnst[5];
		double r2 = cnst[6];
		// upper part
		if (coord[0] >= x0 && coord[2] <= z0) {
			double dd = Math.pow(
					(Math.sqrt((coord[0] - x0) * (coord[0] - x0)
							+ (coord[2] - z0) * (coord[2] - z0)) - r1)
							/ a, 2)
					+ Math.pow((coord[1]) / b, 2);
			if (dd <= 1)
				return true;
		}
		// lower part
		if (coord[0] <= x0 && coord[2] >= 0) {
			double dd = Math.pow(
					(Math.sqrt((coord[0] - x0) * (coord[0] - x0) + (coord[2])
							* (coord[2])) - r2)
							/ a, 2)
					+ Math.pow((coord[1]) / b, 2);
			if (dd <= 1)
				return true;
		}

		return false;
	}

	/**
	 * Based on input coordinates, find if particle is in this bodypart
	 * @param coord coord
	 * @return true or false
	 */
	public boolean inHeart(double[] coord) {
		double[][] tabt = (double[][]) resources.getObject("heart.constant");
		double[] cnstt = tabt[mode - 1];
		double alfa1 = cnstt[0];
		double beta1 = cnstt[1];
		double gama1 = cnstt[2];
		double alfa2 = cnstt[3];
		double beta2 = cnstt[4];
		double gama2 = cnstt[5];
		double alfa3 = cnstt[6];
		double beta3 = cnstt[7];
		double gama3 = cnstt[8];
		double x0 = cnstt[17];
		double y0 = cnstt[18];
		double z0 = cnstt[19];
		double vx = cnstt[9];
		double avy = cnstt[10];
		double lavz = cnstt[11];
		double ravz = cnstt[12];
		double ax = cnstt[13];
		double tlvw = cnstt[14];
		//double trvw = cnstt[15];
		double taw = cnstt[16];
		// ------------------------
		double x1 = alfa1 * (coord[0] - x0) + beta1 * (coord[1] - y0) + gama1
				* (coord[2] - z0);
		double y1 = alfa2 * (coord[0] - x0) + beta2 * (coord[1] - y0) + gama2
				* (coord[2] - z0);
		double z1 = alfa3 * (coord[0] - x0) + beta3 * (coord[1] - y0) + gama3
				* (coord[2] - z0);
		// System.out.println("  "+x1+"  "+y1+"   "+z1);
		// left ventricle-WALL+CONTENTS
		if (x1 >= 0) {
			double dd = Math.pow((x1) / vx, 2) + Math.pow((y1) / avy, 2)
					+ Math.pow((z1) / lavz, 2);
			if (dd <= 1)
				return true;
		}
		// right ventricle -WALL+CONTENTS and dd previous>1---is taken into
		// account
		if (x1 >= 0 && z1 < 0) {
			double dd = Math.pow((x1) / vx, 2) + Math.pow((y1) / avy, 2)
					+ Math.pow((z1) / ravz, 2);
			if (dd <= 1)
				return true;
		}
		// left atrium-WALL+CONTENTS
		// part 1
		if (x1 < 0 && z1 >= 0) {
			double dd = Math.pow((x1) / ax, 2) + Math.pow((y1) / avy, 2)
					+ Math.pow((z1) / lavz, 2);
			if (dd <= 1)
				return true;
		}
		// part2
		if (x1 < 0 && z1 < 0) {
			double dd = Math.pow((x1) / ax, 2) + Math.pow((y1) / avy, 2)
					+ Math.pow((z1) / (lavz - tlvw + taw), 2);
			if (dd <= 1)
				return true;
		}
		// right atrium-WALL+CONTENTS dd previous>1
		if (x1 < 0 && z1 < 0) {
			double dd = Math.pow((x1) / ax, 2) + Math.pow((y1) / avy, 2)
					+ Math.pow((z1) / ravz, 2);
			if (dd <= 1)
				return true;
		}

		return false;
	}

	/**
	 * Based on input coordinates, find if particle is in this bodypart
	 * @param coord coord
	 * @return true or false
	 */
	public boolean inKidneys(double[] coord) {

		double[][] tabt = (double[][]) resources.getObject("kidneys.constant");
		double[] cnstt = tabt[mode - 1];
		double a = cnstt[0];
		double b = cnstt[1];
		double c = cnstt[2];
		double x0 = cnstt[3];
		double y0 = cnstt[4];
		double z0 = cnstt[5];
		double x1 = cnstt[6];

		if (Math.abs(coord[0]) >= x1) {
			// left
			double dd = Math.pow((coord[0] - x0) / a, 2)
					+ Math.pow((coord[1] - y0) / b, 2)
					+ Math.pow((coord[2] - z0) / c, 2);
			if (dd <= 1)
				return true;
			// right
			double dd1 = Math.pow((coord[0] + x0) / a, 2)
					+ Math.pow((coord[1] - y0) / b, 2)
					+ Math.pow((coord[2] - z0) / c, 2);
			if (dd1 <= 1)
				return true;

		}
		return false;
	}

	/**
	 * Based on input coordinates, find if particle is in this bodypart
	 * @param coord coord
	 * @return true or false
	 */
	public boolean inLiver(double[] coord) {
		double[][] tabt = (double[][]) resources.getObject("liver.constant");
		double[] cnstt = tabt[mode - 1];
		double a = cnstt[0];
		double b = cnstt[1];
		double xm = cnstt[2];
		double ym = cnstt[3];
		double zm = cnstt[4];
		double z1 = cnstt[5];
		double z2 = cnstt[6];

		if (coord[2] >= z1 && coord[2] <= z2) {
			double dd = Math.pow((coord[0]) / a, 2)
					+ Math.pow((coord[1]) / b, 2);
			if (dd <= 1) {
				double ddd = (coord[0] / xm) + (coord[1] / ym)
						- (coord[2] / zm);
				if (ddd <= -1)
					return true;
			}
		}

		return false;
	}

	/**
	 * Based on input coordinates, find if particle is in this bodypart
	 * @param coord coord
	 * @return true or false
	 */
	public boolean inLungs(double[] coord) {
		double[][] tabt = (double[][]) resources.getObject("lungs.constant");
		double[] cnstt = tabt[mode - 1];
		double a = cnstt[0];
		double b = cnstt[1];
		double c = cnstt[2];
		double x0 = cnstt[3];
		double z0 = cnstt[4];
		double x1r = cnstt[5];
		double y1r = cnstt[6];
		double z1r = cnstt[7];
		double z2r = cnstt[8];
		double x1l = cnstt[9];
		double y1l = cnstt[10];
		double z2l = cnstt[11];
		// right lung
		if (coord[2] >= z0) {
			if ((coord[2] >= z1r && coord[2] <= z2r) && (coord[1] < y1r)) {
				double dd = Math.pow((coord[0] + x0) / a, 2)
						+ Math.pow((coord[1]) / b, 2)
						+ Math.pow((coord[2] - z0) / c, 2);
				if ((dd <= 1) && (coord[0] <= x1r))
					return true;
			} else {
				double dd = Math.pow((coord[0] + x0) / a, 2)
						+ Math.pow((coord[1]) / b, 2)
						+ Math.pow((coord[2] - z0) / c, 2);
				if (dd <= 1)
					return true;
			}
		}
		// left
		if (coord[2] >= z0) {
			if ((coord[2] >= z0 && coord[2] <= z2l) && (coord[1] < y1l)) {
				double dd = Math.pow((coord[0] - x0) / a, 2)
						+ Math.pow((coord[1]) / b, 2)
						+ Math.pow((coord[2] - z0) / c, 2);
				if ((dd <= 1) && (coord[0] >= x1l))
					return true;
			} else {
				double dd = Math.pow((coord[0] - x0) / a, 2)
						+ Math.pow((coord[1]) / b, 2)
						+ Math.pow((coord[2] - z0) / c, 2);
				if (dd <= 1)
					return true;
			}
		}

		return false;
	}

	/**
	 * Based on input coordinates, find if particle is in this bodypart
	 * @param coord coord
	 * @return true or false
	 */
	public boolean inOvaries(double[] coord) {
		double[][] tabt = (double[][]) resources.getObject("ovaries.constant");
		double[] cnstt = tabt[mode - 1];
		double a = cnstt[0];
		double b = cnstt[1];
		double c = cnstt[2];
		double x0 = cnstt[3];
		double z0 = cnstt[4];
		// left ovaries
		double dd = Math.pow((coord[0] - x0) / a, 2)
				+ Math.pow((coord[1]) / b, 2)
				+ Math.pow((coord[2] - z0) / c, 2);
		if (dd <= 1)
			return true;
		// right
		double dd1 = Math.pow((coord[0] + x0) / a, 2)
				+ Math.pow((coord[1]) / b, 2)
				+ Math.pow((coord[2] - z0) / c, 2);
		if (dd1 <= 1)
			return true;

		return false;
	}

	/**
	 * Based on input coordinates, find if particle is in this bodypart
	 * @param coord coord
	 * @return true or false
	 */
	public boolean inPancreas(double[] coord) {
		double[][] tabt = (double[][]) resources.getObject("pancreas.constant");
		double[] cnstt = tabt[mode - 1];
		double a = cnstt[0];
		double b = cnstt[1];
		double c = cnstt[2];
		double x0 = cnstt[3];
		double z0 = cnstt[4];
		double x1 = cnstt[5];

		if (coord[0] >= x0) {
			if (coord[0] >= x1) {
				if (coord[2] >= z0) {
					double dd = Math.pow((coord[0] - x0) / a, 2)
							+ Math.pow((coord[1]) / b, 2)
							+ Math.pow((coord[2] - z0) / c, 2);
					if (dd <= 1)
						return true;
				}
			} else {
				double dd = Math.pow((coord[0] - x0) / a, 2)
						+ Math.pow((coord[1]) / b, 2)
						+ Math.pow((coord[2] - z0) / c, 2);
				if (dd <= 1)
					return true;

			}
		}
		return false;
	}

	/**
	 * Based on input coordinates, find if particle is in this bodypart
	 * @param coord coord
	 * @return true or false
	 */
	public boolean inSpleen(double[] coord) {
		double[][] tabt = (double[][]) resources.getObject("spleen.constant");
		double[] cnstt = tabt[mode - 1];
		double a = cnstt[0];
		double b = cnstt[1];
		double c = cnstt[2];
		double x0 = cnstt[3];
		double y0 = cnstt[4];
		double z0 = cnstt[5];

		double dd = Math.pow((coord[0] - x0) / a, 2)
				+ Math.pow((coord[1] - y0) / b, 2)
				+ Math.pow((coord[2] - z0) / c, 2);
		if (dd <= 1)
			return true;

		return false;
	}

	/**
	 * Based on input coordinates, find if particle is in this bodypart
	 * @param coord coord
	 * @return true or false
	 */
	public boolean inThymus(double[] coord) {
		double[][] tabt = (double[][]) resources.getObject("thymus.constant");
		double[] cnstt = tabt[mode - 1];
		double a = cnstt[0];
		double b = cnstt[1];
		double c = cnstt[2];
		double y0 = cnstt[3];
		double z0 = cnstt[4];

		double dd = Math.pow((coord[0]) / a, 2)
				+ Math.pow((coord[1] - y0) / b, 2)
				+ Math.pow((coord[2] - z0) / c, 2);
		if (dd <= 1)
			return true;

		return false;
	}

	/**
	 * Based on input coordinates, find if particle is in this bodypart
	 * @param coord coord
	 * @return true or false
	 */
	public boolean inThyroid(double[] coord) {
		double[][] tabt = (double[][]) resources.getObject("thyroid.constant");
		double[] cnstt = tabt[mode - 1];
		double rr = cnstt[0];// R
		double r = cnstt[1];// r
		double c = cnstt[2];
		double y0 = cnstt[3];
		double[][] tab = (double[][]) resources.getObject("trunk.constant");
		double[] abc = tab[mode - 1];
		double ct = abc[2];
		double tau = 0.0;

		if (coord[2] - ct >= 0 && coord[2] - ct <= c) {
			if (coord[2] - ct >= 0 && coord[2] - ct <= 0.25 * c) {
				tau = 1 + ((Math.sqrt(2) - 2) / 2)
						* ((coord[2] - ct) / (0.25 * c));// System.out.println("tau "+tau);
			} else {
				if (coord[2] - ct >= 0.25 * c && coord[2] - ct <= c) {
					tau = ((-1 + 2 * Math.sqrt(2)) / 3)
							+ ((2 - Math.sqrt(2)) / 2)
							* ((coord[2] - ct) / (0.75 * c));
				}
			}
			// -----------test------------------------------
			double d = Math.pow((coord[0]), 2) + Math.pow((coord[1] - y0), 2);
			if (d <= rr * rr) {
				if (d >= r * r) {
					if (coord[1] <= y0) {
						double d1 = Math.pow(
								(coord[1] - y0 - Math.abs(coord[0])), 2);
						double d2 = 2
								* tau
								* tau
								* (Math.pow((coord[0]), 2) + Math.pow(
										(coord[1] - y0), 2));
						if (d1 >= d2)
							return true;
					}
				}
			}
		}
		return false;
	}

	/**
	 * Based on input coordinates, find if particle is in this bodypart
	 * @param coord coord
	 * @return true or false
	 */
	public boolean inUrinaryBladder(double[] coord) {
		double[][] tabt = (double[][]) resources
				.getObject("urinarybladder.constant");
		double[] cnstt = tabt[mode - 1];
		double a = cnstt[0];
		double b = cnstt[1];
		double c = cnstt[2];
		//double d = cnstt[3];
		double y0 = cnstt[4];
		double z0 = cnstt[5];
		// consider wall+contents together
		double dd = Math.pow((coord[0]) / a, 2)
				+ Math.pow((coord[1] - y0) / b, 2)
				+ Math.pow((coord[2] - z0) / c, 2);// left testes
		if (dd <= 1)
			return true;

		return false;
	}

	/**
	 * Based on input coordinates, find if particle is in this bodypart
	 * @param coord coord
	 * @return true or false
	 */
	public boolean inUterus(double[] coord) {
		double[][] tabt = (double[][]) resources.getObject("uterus.constant");
		double[] cnstt = tabt[mode - 1];
		double a = cnstt[0];
		double b = cnstt[1];
		double c = cnstt[2];
		double y0 = cnstt[3];
		double z0 = cnstt[4];
		double y1 = cnstt[5];

		if (coord[1] >= y1) {
			double d = Math.pow((coord[0]) / a, 2)
					+ Math.pow((coord[1] - y0) / b, 2)
					+ Math.pow((coord[2] - z0) / c, 2);// left testes
			if (d <= 1)
				return true;
		}

		return false;
	}

	// ----------------------------------------------------------------------------------------------------
}
