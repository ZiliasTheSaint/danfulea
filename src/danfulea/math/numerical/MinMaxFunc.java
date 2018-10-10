package danfulea.math.numerical;

/**
 * Min. and Max. of functions
 * Based on literature, including Numerical Recipes in C (Cambridge Univ.).
 * @author Dan Fulea, 11 OCT. 2006
 */
public class MinMaxFunc {
	private Function func;
	public static boolean failB = false;
	public static String failS = "";

	public static double GOLD = 1.618034;
	public static double GLIMIT = 100.0;
	public static double TINY = 1.0e-20;
	// Here GOLD is the default ratio by which successive intervals are
	// magnified; GLIMIT is the
	// maximum magnification allowed for a parabolic-fit step.
	// #define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);
	// public static double a_SHFT=0.0;
	// public static double b_SHFT=0.0;
	// public static double c_SHFT=0.0;
	// public static double d_SHFT=0.0;

	public static double ax_mnbrak = 0.0;
	public static double bx_mnbrak = 0.0;
	public static double cx_mnbrak = 0.0;
	public static double fa_mnbrak = 0.0;
	public static double fb_mnbrak = 0.0;
	public static double fc_mnbrak = 0.0;

	public static double R = 0.61803399; // The golden ratios.
	public static double C = (1.0 - R);
	public static double xmin_golden = 0.0;

	public static int ITMAX = 100;
	public static double CGOLD = 0.3819660;
	public static double ZEPS = 1.0e-10;
	public static double xmin_brent = 0.0;

	public static double xmin_dbrent = 0.0;
	// ============================================
	public static double TINY1 = 1.0e-10;// A small number.
	public static int NMAX = 5000;
	public static int nfunk_amoeba = 0;

	public static double TINY2 = 1.0e-25; // A small number.
	public static int ITMAX2 = 200;
	public static double TOL = 2.0e-4;// Tolerance passed to brent.
	public static int ncom = 0; // Global variables communicate with f1dim.
	// float *pcom,*xicom,(*nrfunc)(float []);
	public static double[] pcom;
	public static double[] xicom;
	public static double fret = 0.0;
	public static int iter_powel = 0;

	// #define ITMAX 200=>2
	public static double EPS = 1.0e-10;
	public static int iter_fr = 0;
	public static double fret_fr = 0.0;
	public static double fret_dl = 0.0;

	// #define ITMAX 200 Maximum allowed number of iterations.
	public static double EPS2 = 3.0e-8;// Machine precision.
	public static double TOLX = 4.0 * EPS2;// ) Convergence criterion on x
											// values.
	public static double STPMX = 100.0;
	public static double ALF = 1.0e-4;// Ensures sufficient decrease in function
										// value.
	// public static double TOLX =1.0e-7;// Convergence criterion on ?x.
	public static int check_ln = 0;
	public static double f_ln = 0.0;

	public static int iter_df = 0;
	public static double fret_df = 0.0;

	/**
	 * Constuctor. A user supplied function is an object whose class implements the Function interface
	 * @param func the function
	 */
	public MinMaxFunc(Function func) {
		this.func = func;
	}

	/**
	 * Change the sign of a if b is negative!
	 * @param a a
	 * @param b b
	 * @return a or -a
	 */
	public static double SIGN(double a, double b) {
		return b >= 0.0 ? a : -a;
	}

	/*
	 * public static void SWAP(double a, double b) { a_swap=a;b_swap=b; double
	 * temp=a_swap; a_swap=b_swap; b_swap=temp; }
	 * 
	 * public static void SHFT(double a, double b, double c, double d) {
	 * a_SHFT=a;b_SHFT=b;c_SHFT=c;d_SHFT=d;
	 * 
	 * a_SHFT=b_SHFT; b_SHFT=c_SHFT; c_SHFT=d_SHFT; }
	 */
	/*
	 * Routine for Initially Bracketing a Minimum
	 * 
	 * The preceding discussion has assumed that you are able to bracket the
	 * minimum in the first place. We consider this initial bracketing to be an
	 * essential part of any one-dimensional minimization. There are some
	 * one-dimensional algorithms that do not require a rigorous initial
	 * bracketing. However, we would never trade the secure feeling of knowing
	 * that a minimum is “in there somewhere” for the dubious reduction of
	 * function evaluations that these nonbracketing routines may promise.
	 * Please bracket your minima (or, for that matter, your zeros) before
	 * isolating them! There is not much theory as to how to do this bracketing.
	 * Obviously you want to step downhill. But how far? We like to take larger
	 * and larger steps, starting with some (wild?) initial guess and then
	 * increasing the stepsize at each step either by a constant factor, or else
	 * by the result of a parabolic extrapolation of the preceding points that
	 * is designed to take us to the extrapolated turning point. It doesn’t much
	 * matter if the steps get big. After all, we are stepping downhill, so we
	 * already have the left and middle points of the bracketing triplet. We
	 * just need to take a big enough step to stop the downhill trend and get a
	 * high third point.
	 */
	// void mnbrak(float *ax, float *bx, float *cx, float *fa, float *fb, float
	// *fc,
	// float (*func)(float))
	/**
	 * Given a function func, and given distinct initial points ax and bx, this routine searches in 
	 * the downhill direction (defined by the function as evaluated at the initial points) and returns 
	 * new points ax_mnbrak, bx_mnbrak, cx_mnbrak that bracket a minimum of the function. Also 
	 * returned are the function values at the three points, fa_mnbrak, fb_mnbrak, and fc_mnbrak.
	 * @param ax ax
	 * @param bx bx
	 */
	public void mnbrak(double ax, double bx)// , float *cx, float *fa, float
											// *fb, float *fc,
	// Given a function func, and given distinct initial points ax and bx, this
	// routine searches in
	// the downhill direction (defined by the function as evaluated at the
	// initial points) and returns
	// new points ax, bx, cx that bracket a minimum of the function. Also
	// returned are the function
	// values at the three points, fa, fb, and fc.
	{
		ax_mnbrak = ax;
		bx_mnbrak = bx;

		double ulim = 0.0;
		double u = 0.0;
		double r = 0.0;
		double q = 0.0;
		double fu = 0.0;
		double dum = 0.0;
		fa_mnbrak = func.F(ax_mnbrak);// *fa=(*func)(*ax);
		fb_mnbrak = func.F(bx_mnbrak);// *fb=(*func)(*bx);
		if (fb_mnbrak > fa_mnbrak) {// Switch roles of a and b so that we can go
									// downhill in the direction from a to b.
									// SHFT(dum,*ax,*bx,dum)
			dum = ax_mnbrak;
			ax_mnbrak = bx_mnbrak;
			bx_mnbrak = dum;
			// SHFT(dum,*fb,*fa,dum)
			dum = fb_mnbrak;
			fb_mnbrak = fa_mnbrak;
			fa_mnbrak = dum;
		}
		cx_mnbrak = bx_mnbrak + GOLD * (bx_mnbrak - ax_mnbrak);// First guess
																// for c.
		fc_mnbrak = func.F(cx_mnbrak);
		while (fb_mnbrak > fc_mnbrak) {// Keep returning here until we bracket.
			r = (bx_mnbrak - ax_mnbrak) * (fb_mnbrak - fc_mnbrak);// Compute u
																	// by
																	// parabolic
																	// extrapolation
																	// from
			// a, b, c. TINY is used to prevent any possible division by zero.
			q = (bx_mnbrak - cx_mnbrak) * (fb_mnbrak - fa_mnbrak);
			u = (bx_mnbrak)
					- ((bx_mnbrak - cx_mnbrak) * q - (bx_mnbrak - ax_mnbrak)
							* r)
					/ (2.0 * SIGN(Math.max(Math.abs(q - r), TINY), q - r));
			ulim = (bx_mnbrak) + GLIMIT * (cx_mnbrak - bx_mnbrak);
			// We won’t go farther than this. Test various possibilities:
			if ((bx_mnbrak - u) * (u - cx_mnbrak) > 0.0) {// Parabolic u is
															// between b and c:
															// try it.
				fu = func.F(u);
				if (fu < fc_mnbrak) {// Got a minimum between b and c.
					ax_mnbrak = bx_mnbrak;
					bx_mnbrak = u;
					fa_mnbrak = fb_mnbrak;
					fb_mnbrak = fu;
					return;
				} else if (fu > fb_mnbrak) {// Got a minimum between between a
											// and u.
					cx_mnbrak = u;
					fc_mnbrak = fu;
					return;
				}
				u = (cx_mnbrak) + GOLD * (cx_mnbrak - bx_mnbrak);// Parabolic
																	// fit was
																	// no use.
																	// Use
																	// default
																	// magnification.
				fu = func.F(u);
			} else if ((cx_mnbrak - u) * (u - ulim) > 0.0) {// Parabolic fit is
															// between c and its
															// allowed limit.
				fu = func.F(u);
				if (fu < fc_mnbrak) {
					// SHFT(*bx,*cx,u,*cx+GOLD*(*cx-*bx))
					bx_mnbrak = cx_mnbrak;
					cx_mnbrak = u;
					u = cx_mnbrak + GOLD * (cx_mnbrak - bx_mnbrak);
					// SHFT(*fb,*fc,fu,(*func)(u))
					fb_mnbrak = fc_mnbrak;
					fc_mnbrak = fu;
					fu = func.F(u);
				}
			} else if ((u - ulim) * (ulim - cx_mnbrak) >= 0.0) {// Limit
																// parabolic u
																// to maximum
																// allowed
																// value.
				u = ulim;
				fu = func.F(u);
			} else {// Reject parabolic u, use default magnification.
				u = cx_mnbrak + GOLD * (cx_mnbrak - bx_mnbrak);
				fu = func.F(u);
			}
			// SHFT(*ax,*bx,*cx,u) Eliminate oldest point and continue.
			ax_mnbrak = bx_mnbrak;
			bx_mnbrak = cx_mnbrak;
			cx_mnbrak = u;
			// SHFT(*fa,*fb,*fc,fu)
			fa_mnbrak = fb_mnbrak;
			fb_mnbrak = fc_mnbrak;
			fc_mnbrak = fu;
		}
	}

	/**
	 * Same as above but this is related to a local function and it is used internally.
	 * @param ax ax
	 * @param bx bx
	 * @param namS the name of local function to be used.
	 */
	public void mnbrak1(double ax, double bx, String namS)// , float *cx, float
															// *fa, float *fb,
															// float *fc,
	// Given a function func, and given distinct initial points ax and bx, this
	// routine searches in
	// the downhill direction (defined by the function as evaluated at the
	// initial points) and returns
	// new points ax, bx, cx that bracket a minimum of the function. Also
	// returned are the function
	// values at the three points, fa, fb, and fc.
	{
		ax_mnbrak = ax;
		bx_mnbrak = bx;

		double ulim = 0.0;
		double u = 0.0;
		double r = 0.0;
		double q = 0.0;
		double fu = 0.0;
		double dum = 0.0;
		fa_mnbrak = localFunc(ax_mnbrak, namS);// *fa=(*func)(*ax);
		fb_mnbrak = localFunc(bx_mnbrak, namS);// *fb=(*func)(*bx);
		if (fb_mnbrak > fa_mnbrak) {// Switch roles of a and b so that we can go
									// downhill in the direction from a to b.
									// SHFT(dum,*ax,*bx,dum)
			dum = ax_mnbrak;
			ax_mnbrak = bx_mnbrak;
			bx_mnbrak = dum;
			// SHFT(dum,*fb,*fa,dum)
			dum = fb_mnbrak;
			fb_mnbrak = fa_mnbrak;
			fa_mnbrak = dum;
		}
		cx_mnbrak = bx_mnbrak + GOLD * (bx_mnbrak - ax_mnbrak);// First guess
																// for c.
		fc_mnbrak = localFunc(cx_mnbrak, namS);
		while (fb_mnbrak > fc_mnbrak) {// Keep returning here until we bracket.
			r = (bx_mnbrak - ax_mnbrak) * (fb_mnbrak - fc_mnbrak);// Compute u
																	// by
																	// parabolic
																	// extrapolation
																	// from
			// a, b, c. TINY is used to prevent any possible division by zero.
			q = (bx_mnbrak - cx_mnbrak) * (fb_mnbrak - fa_mnbrak);
			u = (bx_mnbrak)
					- ((bx_mnbrak - cx_mnbrak) * q - (bx_mnbrak - ax_mnbrak)
							* r)
					/ (2.0 * SIGN(Math.max(Math.abs(q - r), TINY), q - r));
			ulim = (bx_mnbrak) + GLIMIT * (cx_mnbrak - bx_mnbrak);
			// We won’t go farther than this. Test various possibilities:
			if ((bx_mnbrak - u) * (u - cx_mnbrak) > 0.0) {// Parabolic u is
															// between b and c:
															// try it.
				fu = localFunc(u, namS);
				if (fu < fc_mnbrak) {// Got a minimum between b and c.
					ax_mnbrak = bx_mnbrak;
					bx_mnbrak = u;
					fa_mnbrak = fb_mnbrak;
					fb_mnbrak = fu;
					return;
				} else if (fu > fb_mnbrak) {// Got a minimum between between a
											// and u.
					cx_mnbrak = u;
					fc_mnbrak = fu;
					return;
				}
				u = (cx_mnbrak) + GOLD * (cx_mnbrak - bx_mnbrak);// Parabolic
																	// fit was
																	// no use.
																	// Use
																	// default
																	// magnification.
				fu = localFunc(u, namS);
			} else if ((cx_mnbrak - u) * (u - ulim) > 0.0) {// Parabolic fit is
															// between c and its
															// allowed limit.
				fu = localFunc(u, namS);
				if (fu < fc_mnbrak) {
					// SHFT(*bx,*cx,u,*cx+GOLD*(*cx-*bx))
					bx_mnbrak = cx_mnbrak;
					cx_mnbrak = u;
					u = cx_mnbrak + GOLD * (cx_mnbrak - bx_mnbrak);
					// SHFT(*fb,*fc,fu,(*func)(u))
					fb_mnbrak = fc_mnbrak;
					fc_mnbrak = fu;
					fu = localFunc(u, namS);
				}
			} else if ((u - ulim) * (ulim - cx_mnbrak) >= 0.0) {// Limit
																// parabolic u
																// to maximum
																// allowed
																// value.
				u = ulim;
				fu = localFunc(u, namS);
			} else {// Reject parabolic u, use default magnification.
				u = cx_mnbrak + GOLD * (cx_mnbrak - bx_mnbrak);
				fu = localFunc(u, namS);
			}
			// SHFT(*ax,*bx,*cx,u) Eliminate oldest point and continue.
			ax_mnbrak = bx_mnbrak;
			bx_mnbrak = cx_mnbrak;
			cx_mnbrak = u;
			// SHFT(*fa,*fb,*fc,fu)
			fa_mnbrak = fb_mnbrak;
			fb_mnbrak = fc_mnbrak;
			fc_mnbrak = fu;
		}
	}

	/*
	 * (Because of the housekeeping involved in moving around three or four
	 * points and their function values, the above program ends up looking
	 * deceptively formidable. That is true of several other programs in this
	 * chapter as well. The underlying ideas, however, are quite simple.)
	 * Routine for Golden Section Search
	 */
	// float golden(float ax, float bx, float cx, float (*f)(float), float tol,
	// float *xmin)
	/**
	 * Given a function f, and given a bracketing triplet of abscissas ax, bx, cx (such that bx is 
	 * between ax and cx, and f(bx) is less than both f(ax) and f(cx)), this routine performs a 
	 * golden section search for the minimum, isolating it to a fractional precision of about tol. The 
	 * abscissa of the minimum is returned as xmin_golden, and the minimum function value is returned as 
	 * golden, the returned function value.
	 * @param ax ax
	 * @param bx bx
	 * @param cx cx
	 * @param tol tol
	 * @return the function value
	 */
	public double golden(double ax, double bx, double cx, double tol)
	// Given a function f, and given a bracketing triplet of abscissas ax, bx,
	// cx (such that bx is
	// between ax and cx, and f(bx) is less than both f(ax) and f(cx)), this
	// routine performs a
	// golden section search for the minimum, isolating it to a fractional
	// precision of about tol. The
	// abscissa of the minimum is returned as xmin, and the minimum function
	// value is returned as
	// golden, the returned function value.
	{
		double f1 = 0.0;
		double f2 = 0.0;
		double x0 = 0.0;
		double x1 = 0.0;
		double x2 = 0.0;
		double x3 = 0.0;

		x0 = ax; // At any given time we will keep track of four
					// points, x0,x1,x2,x3.
		x3 = cx;
		if (Math.abs(cx - bx) > Math.abs(bx - ax)) {// Make x0 to x1 the smaller
													// segment,
			x1 = bx;
			x2 = bx + C * (cx - bx);// and fill in the new point to be tried.
		} else {
			x2 = bx;
			x1 = bx - C * (bx - ax);
		}
		f1 = func.F(x1); // The initial function evaluations. Note that
		// we never need to evaluate the function at the original endpoints.
		f2 = func.F(x2);
		while (Math.abs(x3 - x0) > tol * (Math.abs(x1) + Math.abs(x2))) {
			if (f2 < f1) {// One possible outcome,
							// SHFT3(x0,x1,x2,R*x1+C*x3) its housekeeping,
				x0 = x1;
				x1 = x2;
				x2 = R * x1 + C * x3;
				// SHFT2(f1,f2,(*f)(x2)) and a new function evaluation.
				f1 = f2;
				f2 = func.F(x2);
			} else {// The other outcome,
					// SHFT3(x3,x2,x1,R*x2+C*x0)
				x3 = x2;
				x2 = x1;
				x1 = R * x2 + C * x0;
				// SHFT2(f2,f1,(*f)(x1)) and its new function evaluation.
				f2 = f1;
				f1 = func.F(x1);
			}
		}// Back to see if we are done.
		if (f1 < f2) {// We are done. Output the best of the two current values.
			xmin_golden = x1;
			return f1;
		} else {
			xmin_golden = x2;
			return f2;
		}
	}

	/*
	 * Parabolic Interpolation and Brent’s Method We already tipped our hand
	 * about the desirability of parabolic interpolation in the previous
	 * section’s mnbrak routine, but it is now time to be more explicit. A
	 * golden section search is designed to handle, in effect, the worst
	 * possible case of function minimization, with the uncooperative minimum
	 * hunted down and cornered like a scared rabbit. But why assume the worst?
	 * If the function is nicely parabolic near to the minimum— surely the
	 * generic case for sufficiently smooth functions— then the parabola fitted
	 * through any three points ought to take us in a single leap to the
	 * minimum, or at least very near to it The exacting task is to invent a
	 * scheme that relies on a sure-but-slow technique, like golden section
	 * search, when the function is not cooperative, but that switches over to
	 * (10.2.1) when the function allows.
	 */
	// public double brent(float ax, float bx, float cx, float (*f)(float),
	// float tol,
	// float *xmin)
	/**
	 * Given a function f, and given a bracketing triplet of abscissas ax, bx, cx (such that bx is 
	 * between ax and cx, and f(bx) is less than both f(ax) and f(cx)), this routine isolates 
	 * the minimum to a fractional precision of about tol using Brent’s method. The abscissa of 
	 * the minimum is returned as xmin_brent, and the minimum function value is returned as brent, the 
	 * returned function value.
	 * @param ax ax
	 * @param bx bx
	 * @param cx cx
	 * @param tol tol
	 * @return the function value
	 */
	public double brent(double ax, double bx, double cx, double tol)
	// Given a function f, and given a bracketing triplet of abscissas ax, bx,
	// cx (such that bx is
	// between ax and cx, and f(bx) is less than both f(ax) and f(cx)), this
	// routine isolates
	// the minimum to a fractional precision of about tol using Brent’s method.
	// The abscissa of
	// the minimum is returned as xmin, and the minimum function value is
	// returned as brent, the
	// returned function value.
	{
		failB = false;
		int iter = 0;
		double a = 0.0;
		double b = 0.0;
		double d = 0.0;
		double etemp = 0.0;
		double fu = 0.0;
		double fv = 0.0;
		double fw = 0.0;
		double fx = 0.0;
		double p = 0.0;
		double q = 0.0;
		double r = 0.0;
		double tol1 = 0.0;
		double tol2 = 0.0;
		double u = 0.0;
		double v = 0.0;
		double w = 0.0;
		double x = 0.0;
		double xm = 0.0;
		double e = 0.0; // This will be the distance moved on the step before
						// last.

		a = (ax < cx ? ax : cx); // a and b must be in ascending order,
		// but input abscissas need not be.
		b = (ax > cx ? ax : cx);
		x = w = v = bx; // Initializations...
		fw = fv = fx = func.F(x);
		for (iter = 1; iter <= ITMAX; iter++) {// Main program loop.
			xm = 0.5 * (a + b);
			tol2 = 2.0 * (tol1 = tol * Math.abs(x) + ZEPS);
			if (Math.abs(x - xm) <= (tol2 - 0.5 * (b - a))) {// Test for done
																// here.
				xmin_brent = x;
				return fx;
			}
			if (Math.abs(e) > tol1) {// Construct a trial parabolic fit.
				r = (x - w) * (fx - fv);
				q = (x - v) * (fx - fw);
				p = (x - v) * q - (x - w) * r;
				q = 2.0 * (q - r);
				if (q > 0.0)
					p = -p;
				q = Math.abs(q);
				etemp = e;
				e = d;
				if (Math.abs(p) >= Math.abs(0.5 * q * etemp)
						|| p <= q * (a - x) || p >= q * (b - x))
					d = CGOLD * (e = (x >= xm ? a - x : b - x));
				// The above conditions determine the acceptability of the
				// parabolic fit. Here we
				// take the golden section step into the larger of the two
				// segments.
				else {
					d = p / q; // Take the parabolic step.
					u = x + d;
					if (u - a < tol2 || b - u < tol2)
						d = SIGN(tol1, xm - x);
				}
			} else {
				d = CGOLD * (e = (x >= xm ? a - x : b - x));
			}
			u = (Math.abs(d) >= tol1 ? x + d : x + SIGN(tol1, d));
			fu = func.F(u);
			// This is the one function evaluation per iteration.
			if (fu <= fx) {// Now decide what to do with our function
							// evaluation.
				if (u >= x)
					a = x;
				else
					b = x;
				// SHFT(v,w,x,u) Housekeeping follows:
				v = w;
				w = x;
				x = u;
				// SHFT(fv,fw,fx,fu)
				fv = fw;
				fw = fx;
				fx = fu;
			} else {
				if (u < x)
					a = u;
				else
					b = u;
				if (fu <= fw || w == x) {
					v = w;
					w = u;
					fv = fw;
					fw = fu;
				} else if (fu <= fv || v == x || v == w) {
					v = u;
					fv = fu;
				}
			}// Done with housekeeping. Back for another iteration.
		}
		// nrerror("Too many iterations in brent");
		failB = true;
		failS = "Too many iterations in brent";
		xmin_brent = x;// Never get here.
		return fx;
	}

	/**
	 * Same as above but this is related to a local function and it is used internally.
	 * @param ax ax
	 * @param bx bx
	 * @param cx cx
	 * @param tol tol
	 * @param namS the name of local function to be used.
	 * @return the result
	 */
	public double brent1(double ax, double bx, double cx, double tol,
			String namS)
	// Given a function f, and given a bracketing triplet of abscissas ax, bx,
	// cx (such that bx is
	// between ax and cx, and f(bx) is less than both f(ax) and f(cx)), this
	// routine isolates
	// the minimum to a fractional precision of about tol using Brent’s method.
	// The abscissa of
	// the minimum is returned as xmin, and the minimum function value is
	// returned as brent, the
	// returned function value.
	{
		failB = false;
		int iter = 0;
		double a = 0.0;
		double b = 0.0;
		double d = 0.0;
		double etemp = 0.0;
		double fu = 0.0;
		double fv = 0.0;
		double fw = 0.0;
		double fx = 0.0;
		double p = 0.0;
		double q = 0.0;
		double r = 0.0;
		double tol1 = 0.0;
		double tol2 = 0.0;
		double u = 0.0;
		double v = 0.0;
		double w = 0.0;
		double x = 0.0;
		double xm = 0.0;
		double e = 0.0; // This will be the distance moved on the step before
						// last.

		a = (ax < cx ? ax : cx); // a and b must be in ascending order,
		// but input abscissas need not be.
		b = (ax > cx ? ax : cx);
		x = w = v = bx; // Initializations...
		fw = fv = fx = localFunc(x, namS);
		for (iter = 1; iter <= ITMAX; iter++) {// Main program loop.
			xm = 0.5 * (a + b);
			tol2 = 2.0 * (tol1 = tol * Math.abs(x) + ZEPS);
			if (Math.abs(x - xm) <= (tol2 - 0.5 * (b - a))) {// Test for done
																// here.
				xmin_brent = x;
				return fx;
			}
			if (Math.abs(e) > tol1) {// Construct a trial parabolic fit.
				r = (x - w) * (fx - fv);
				q = (x - v) * (fx - fw);
				p = (x - v) * q - (x - w) * r;
				q = 2.0 * (q - r);
				if (q > 0.0)
					p = -p;
				q = Math.abs(q);
				etemp = e;
				e = d;
				if (Math.abs(p) >= Math.abs(0.5 * q * etemp)
						|| p <= q * (a - x) || p >= q * (b - x))
					d = CGOLD * (e = (x >= xm ? a - x : b - x));
				// The above conditions determine the acceptability of the
				// parabolic fit. Here we
				// take the golden section step into the larger of the two
				// segments.
				else {
					d = p / q; // Take the parabolic step.
					u = x + d;
					if (u - a < tol2 || b - u < tol2)
						d = SIGN(tol1, xm - x);
				}
			} else {
				d = CGOLD * (e = (x >= xm ? a - x : b - x));
			}
			u = (Math.abs(d) >= tol1 ? x + d : x + SIGN(tol1, d));
			fu = localFunc(u, namS);
			// This is the one function evaluation per iteration.
			if (fu <= fx) {// Now decide what to do with our function
							// evaluation.
				if (u >= x)
					a = x;
				else
					b = x;
				// SHFT(v,w,x,u) Housekeeping follows:
				v = w;
				w = x;
				x = u;
				// SHFT(fv,fw,fx,fu)
				fv = fw;
				fw = fx;
				fx = fu;
			} else {
				if (u < x)
					a = u;
				else
					b = u;
				if (fu <= fw || w == x) {
					v = w;
					w = u;
					fv = fw;
					fw = fu;
				} else if (fu <= fv || v == x || v == w) {
					v = u;
					fv = fu;
				}
			}// Done with housekeeping. Back for another iteration.
		}
		// nrerror("Too many iterations in brent");
		failB = true;
		failS = "Too many iterations in brent";
		xmin_brent = x;// Never get here.
		return fx;
	}

	/*
	 * Here we want to accomplish precisely the same goal as in the previous
	 * section, namely to isolate a functional minimum that is bracketed by the
	 * triplet of abscissas (a, b, c), but utilizing an additional capability to
	 * compute the function’s first derivative as well as its value. In
	 * principle, we might simply search for a zero of the derivative, ignoring
	 * the function value information, using a root finder like rtflsp or zbrent
	 * (§§9.2–9.3). It doesn’t take long to reject that idea: How do we
	 * distinguish maxima from minima? Where do we go from initial conditions
	 * where the derivatives on one or both of the outer bracketing points
	 * indicate that “downhill” is in the direction out of the bracketed
	 * interval?
	 */
	// float dbrent(float ax, float bx, float cx, float (*f)(float),
	// float (*df)(float), float tol, float *xmin)
	/**
	 * Given a function f and its derivative function df, and given a bracketing triplet of abscissas ax, 
	 * bx, cx [such that bx is between ax and cx, and f(bx) is less than both f(ax) and f(cx)], 
	 * this routine isolates the minimum to a fractional precision of about tol using a modification of 
	 * Brent’s method that uses derivatives. The abscissa of the minimum is returned as xmin_dbrent, and 
	 * the minimum function value is returned as dbrent, the returned function value.
	 * @param ax ax
	 * @param bx bx
	 * @param cx cx
	 * @param tol tol
	 * @return the function value
	 */
	public double dbrent(double ax, double bx, double cx, double tol)// , float
																		// *xmin)
	// Given a function f and its derivative function df, and given a bracketing
	// triplet of abscissas ax,
	// bx, cx [such that bx is between ax and cx, and f(bx) is less than both
	// f(ax) and f(cx)],
	// this routine isolates the minimum to a fractional precision of about tol
	// using a modification of
	// Brent’s method that uses derivatives. The abscissa of the minimum is
	// returned as xmin, and
	// the minimum function value is returned as dbrent, the returned function
	// value.
	{
		failB = false;
		int iter = 0;
		int ok1 = 0;
		int ok2 = 0; // Will be used as flags for whether proposed
		// steps are acceptable or not.
		double a = 0.0;
		double b = 0.0;
		double d = 0.0;
		double d1 = 0.0;
		double d2 = 0.0;
		double du = 0.0;
		double dv = 0.0;
		double dw = 0.0;
		double dx = 0.0;
		double e = 0.0;
		double fu = 0.0;
		double fv = 0.0;
		double fw = 0.0;
		double fx = 0.0;
		double olde = 0.0;
		double tol1 = 0.0;
		double tol2 = 0.0;
		double u = 0.0;
		double u1 = 0.0;
		double u2 = 0.0;
		double v = 0.0;
		double w = 0.0;
		double x = 0.0;
		double xm = 0.0;
		// Comments following will point out only differences from the routine
		// brent. Read that
		// routine first.
		a = (ax < cx ? ax : cx);
		b = (ax > cx ? ax : cx);
		x = w = v = bx;

		// fw=fv=fx=(*f)(x);
		double[] ffd = func.FD(x);
		fw = fv = fx = ffd[0];

		// dw=dv=dx=(*df)(x); All our housekeeping chores are doubled by the
		// necessity of moving
		// derivative values around as well as function values.
		dw = dv = dx = ffd[1];

		for (iter = 1; iter <= ITMAX; iter++) {
			xm = 0.5 * (a + b);
			tol1 = tol * Math.abs(x) + ZEPS;
			tol2 = 2.0 * tol1;
			if (Math.abs(x - xm) <= (tol2 - 0.5 * (b - a))) {
				xmin_dbrent = x;
				return fx;
			}
			if (Math.abs(e) > tol1) {
				d1 = 2.0 * (b - a); // Initialize these d’s to an out-of-bracket
									// value.
				d2 = d1;
				if (dw != dx)
					d1 = (w - x) * dx / (dx - dw); // Secant method with one
													// point.
				if (dv != dx)
					d2 = (v - x) * dx / (dx - dv);// And the other.
				// Which of these two estimates of d shall we take? We will
				// insist that they be within
				// the bracket, and on the side pointed to by the derivative at
				// x:
				u1 = x + d1;
				u2 = x + d2;// if (sum!=0.0) ii=i;//if (sum) ii=i;
				if ((a - u1) * (u1 - b) > 0.0 && dx * d1 <= 0.0)
					ok1 = 1;// true
				else
					ok1 = 0;// false
				// ok1 = (a-u1)*(u1-b) > 0.0 && dx*d1 <= 0.0;
				// ok2 = (a-u2)*(u2-b) > 0.0 && dx*d2 <= 0.0;
				if ((a - u2) * (u2 - b) > 0.0 && dx * d2 <= 0.0)
					ok2 = 1;// true
				else
					ok2 = 0;// false
				olde = e; // Movement on the step before last.
				e = d;
				if (ok1 == 1 || ok2 == 1)// if (ok1 || ok2)
				{// Take only an acceptable d, and if both are acceptable, then
					// take
					// the smallest one.
					if (ok1 == 1 && ok2 == 1)// if (ok1 && ok2)
						d = (Math.abs(d1) < Math.abs(d2) ? d1 : d2);
					else if (ok1 == 1)// if (ok1)
						d = d1;
					else
						d = d2;
					if (Math.abs(d) <= Math.abs(0.5 * olde)) {
						u = x + d;
						if (u - a < tol2 || b - u < tol2)
							d = SIGN(tol1, xm - x);
					} else {// Bisect, not golden section.
						d = 0.5 * (e = (dx >= 0.0 ? a - x : b - x));
						// Decide which segment by the sign of the derivative.
					}
				} else {
					d = 0.5 * (e = (dx >= 0.0 ? a - x : b - x));
				}
			} else {
				d = 0.5 * (e = (dx >= 0.0 ? a - x : b - x));
			}
			if (Math.abs(d) >= tol1) {
				u = x + d;
				ffd = func.FD(u);
				fu = ffd[0];// fu=(*f)(u);//@@@@@@@@@
			} else {
				u = x + SIGN(tol1, d);
				ffd = func.FD(u);
				fu = ffd[0];// fu=(*f)(u);//@@@@@@@@@
				if (fu > fx) {// If the minimum step in the downhill direction
								// takes us uphill, then
								// we are done.
					xmin_dbrent = x;
					return fx;
				}
			}
			du = ffd[1];// du=(*df)(u); //Now all the housekeeping,
						// sigh.//@@@@@@@@@
			if (fu <= fx) {
				if (u >= x)
					a = x;
				else
					b = x;
				// MOV3(v,fv,dv, w,fw,dw)
				v = w;
				fv = fw;
				dv = dw;
				// MOV3(w,fw,dw, x,fx,dx)
				w = x;
				fw = fx;
				dw = dx;
				// MOV3(x,fx,dx, u,fu,du)
				x = u;
				fx = fu;
				dx = du;
			} else {
				if (u < x)
					a = u;
				else
					b = u;
				if (fu <= fw || w == x) {
					// MOV3(v,fv,dv, w,fw,dw)
					v = w;
					fv = fw;
					dv = dw;
					// MOV3(w,fw,dw, u,fu,du)
					w = u;
					fw = fu;
					dw = du;
				} else if (fu < fv || v == x || v == w) {
					// MOV3(v,fv,dv, u,fu,du)
					v = u;
					fv = fu;
					dv = du;
				}
			}
		}
		// nrerror("Too many iterations in routine dbrent");
		failB = true;
		failS = "Too many iterations in dbrent";
		return 0.0;// Never get here.
	}

	/**
	 * Same as above but with local functions used internally.
	 * @param ax ax
	 * @param bx bx
	 * @param cx cx
	 * @param tol tol
	 * @param namS the name of local function to be used.
	 * @param namS2 the name of local function to be used.
	 * @return the result
	 */
	public double dbrent1(double ax, double bx, double cx, double tol,
			String namS, String namS2)// , float *xmin)
	// Given a function f and its derivative function df, and given a bracketing
	// triplet of abscissas ax,
	// bx, cx [such that bx is between ax and cx, and f(bx) is less than both
	// f(ax) and f(cx)],
	// this routine isolates the minimum to a fractional precision of about tol
	// using a modification of
	// Brent’s method that uses derivatives. The abscissa of the minimum is
	// returned as xmin, and
	// the minimum function value is returned as dbrent, the returned function
	// value.
	{
		failB = false;
		int iter = 0;
		int ok1 = 0;
		int ok2 = 0; // Will be used as flags for whether proposed
		// steps are acceptable or not.
		double a = 0.0;
		double b = 0.0;
		double d = 0.0;
		double d1 = 0.0;
		double d2 = 0.0;
		double du = 0.0;
		double dv = 0.0;
		double dw = 0.0;
		double dx = 0.0;
		double e = 0.0;
		double fu = 0.0;
		double fv = 0.0;
		double fw = 0.0;
		double fx = 0.0;
		double olde = 0.0;
		double tol1 = 0.0;
		double tol2 = 0.0;
		double u = 0.0;
		double u1 = 0.0;
		double u2 = 0.0;
		double v = 0.0;
		double w = 0.0;
		double x = 0.0;
		double xm = 0.0;
		// Comments following will point out only differences from the routine
		// brent. Read that
		// routine first.
		a = (ax < cx ? ax : cx);
		b = (ax > cx ? ax : cx);
		x = w = v = bx;

		// fw=fv=fx=(*f)(x);
		// double[] ffd=func.FD(x);fw=fv=fx=ffd[0];
		fw = fv = fx = localFunc(x, namS);

		// dw=dv=dx=(*df)(x); All our housekeeping chores are doubled by the
		// necessity of moving
		// derivative values around as well as function values.
		dw = dv = dx = dlocalFunc(x, namS2);// ffd[1];

		for (iter = 1; iter <= ITMAX; iter++) {
			xm = 0.5 * (a + b);
			tol1 = tol * Math.abs(x) + ZEPS;
			tol2 = 2.0 * tol1;
			if (Math.abs(x - xm) <= (tol2 - 0.5 * (b - a))) {
				xmin_dbrent = x;
				return fx;
			}
			if (Math.abs(e) > tol1) {
				d1 = 2.0 * (b - a); // Initialize these d’s to an out-of-bracket
									// value.
				d2 = d1;
				if (dw != dx)
					d1 = (w - x) * dx / (dx - dw); // Secant method with one
													// point.
				if (dv != dx)
					d2 = (v - x) * dx / (dx - dv);// And the other.
				// Which of these two estimates of d shall we take? We will
				// insist that they be within
				// the bracket, and on the side pointed to by the derivative at
				// x:
				u1 = x + d1;
				u2 = x + d2;// if (sum!=0.0) ii=i;//if (sum) ii=i;
				if ((a - u1) * (u1 - b) > 0.0 && dx * d1 <= 0.0)
					ok1 = 1;// true
				else
					ok1 = 0;// false
				// ok1 = (a-u1)*(u1-b) > 0.0 && dx*d1 <= 0.0;
				// ok2 = (a-u2)*(u2-b) > 0.0 && dx*d2 <= 0.0;
				if ((a - u2) * (u2 - b) > 0.0 && dx * d2 <= 0.0)
					ok2 = 1;// true
				else
					ok2 = 0;// false
				olde = e; // Movement on the step before last.
				e = d;
				if (ok1 == 1 || ok2 == 1)// if (ok1 || ok2)
				{// Take only an acceptable d, and if both are acceptable, then
					// take
					// the smallest one.
					if (ok1 == 1 && ok2 == 1)// if (ok1 && ok2)
						d = (Math.abs(d1) < Math.abs(d2) ? d1 : d2);
					else if (ok1 == 1)// if (ok1)
						d = d1;
					else
						d = d2;
					if (Math.abs(d) <= Math.abs(0.5 * olde)) {
						u = x + d;
						if (u - a < tol2 || b - u < tol2)
							d = SIGN(tol1, xm - x);
					} else {// Bisect, not golden section.
						d = 0.5 * (e = (dx >= 0.0 ? a - x : b - x));
						// Decide which segment by the sign of the derivative.
					}
				} else {
					d = 0.5 * (e = (dx >= 0.0 ? a - x : b - x));
				}
			} else {
				d = 0.5 * (e = (dx >= 0.0 ? a - x : b - x));
			}
			if (Math.abs(d) >= tol1) {
				u = x + d;
				// ffd=func.FD(u);fu=ffd[0];//fu=(*f)(u);//@@@@@@@@@
				fu = localFunc(u, namS);
			} else {
				u = x + SIGN(tol1, d);
				// ffd=func.FD(u);fu=ffd[0];//fu=(*f)(u);//@@@@@@@@@
				fu = localFunc(u, namS);
				if (fu > fx) {// If the minimum step in the downhill direction
								// takes us uphill, then
								// we are done.
					xmin_dbrent = x;
					return fx;
				}
			}
			// du=ffd[1];//du=(*df)(u); //Now all the housekeeping,
			// sigh.//@@@@@@@@@
			du = dlocalFunc(u, namS2);
			if (fu <= fx) {
				if (u >= x)
					a = x;
				else
					b = x;
				// MOV3(v,fv,dv, w,fw,dw)
				v = w;
				fv = fw;
				dv = dw;
				// MOV3(w,fw,dw, x,fx,dx)
				w = x;
				fw = fx;
				dw = dx;
				// MOV3(x,fx,dx, u,fu,du)
				x = u;
				fx = fu;
				dx = du;
			} else {
				if (u < x)
					a = u;
				else
					b = u;
				if (fu <= fw || w == x) {
					// MOV3(v,fv,dv, w,fw,dw)
					v = w;
					fv = fw;
					dv = dw;
					// MOV3(w,fw,dw, u,fu,du)
					w = u;
					fw = fu;
					dw = du;
				} else if (fu < fv || v == x || v == w) {
					// MOV3(v,fv,dv, u,fu,du)
					v = u;
					fv = fu;
					dv = du;
				}
			}
		}
		// nrerror("Too many iterations in routine dbrent");
		failB = true;
		failS = "Too many iterations in dbrent";
		return 0.0;// Never get here.
	}

	/*
	 * With this section we begin consideration of multidimensional
	 * minimization, that is, finding the minimum of a function of more than one
	 * independent variable. This section stands apart from those which follow,
	 * however: All of the algorithms after this section will make explicit use
	 * of a one-dimensional minimization algorithm as a part of their
	 * computational strategy. This section implements an entirely
	 * self-contained strategy, in which one-dimensional minimization does not
	 * figure. The downhill simplex method is due to Nelder and Mead [1]. The
	 * method requires only function evaluations, not derivatives. It is not
	 * very efficient in terms of the number of function evaluations that it
	 * requires. Powell’s method (§10.5) is almost surely faster in all likely
	 * applications. However, the downhill simplex method may frequently be the
	 * best method to use if the figure of merit is “get something working
	 * quickly” for a problem whose computational burden is small.
	 */
	/**
	 * Multidimensional minimization of the function funk(x) where x[1..ndim] is a vector in ndim 
	 * dimensions, by the downhill simplex method of Nelder and Mead. The matrix p[1..ndim+1][1..ndim] is input. 
	 * Its ndim+1 rows are ndim-dimensional vectors which are the vertices of the starting simplex. Also input is the vector y[1..ndim+1], whose 
	 * components must be preinitialized to the values of funk evaluated at the ndim+1 vertices (rows) of p; and ftol the 
	 * fractional convergence tolerance to be achieved in the function value. On output, p and y will have been reset to ndim+1 new points all within ftol of a minimum 
	 * function value, and nfunk gives the number of function evaluations taken.
	 * @param p p
	 * @param y y
	 * @param ndim ndim
	 * @param ftol ftol
	 */
	public void amoeba(double[][] p, double[] y, int ndim, double ftol)// ,
	// float (*funk)(float []), int *nfunk)
	// Multidimensional minimization of the function funk(x) where x[1..ndim] is
	// a vector in ndim
	// dimensions, by the downhill simplex method of Nelder and Mead. The matrix
	// p[1..ndim+1]
	// [1..ndim] is input. Its ndim+1 rows are ndim-dimensional vectors which
	// are the vertices of
	// the starting simplex. Also input is the vector y[1..ndim+1], whose
	// components must be preinitialized
	// to the values of funk evaluated at the ndim+1 vertices (rows) of p; and
	// ftol the
	// fractional convergence tolerance to be achieved in the function value
	// (n.b.!). On output, p and
	// y will have been reset to ndim+1 new points all within ftol of a minimum
	// function value, and
	// nfunk gives the number of function evaluations taken.
	{
		failB = false;
		// float amotry(float **p, float y[], float psum[], int ndim,
		// float (*funk)(float []), int ihi, float fac);
		int i = 0;
		int ihi = 0;
		int ilo = 0;
		int inhi = 0;
		int j = 0;
		int mpts = ndim + 1;
		double rtol = 0.0;
		double sum = 0.0;
		double swap = 0.0;
		double ysave = 0.0;
		double ytry = 0.0;
		double[] psum = new double[ndim];
		// psum=vector(1,ndim);
		nfunk_amoeba = 0;// *nfunk=0;
		// GET_PSUM
		for (j = 1; j <= ndim; j++) {// \
			for (sum = 0.0, i = 1; i <= mpts; i++)
				sum += p[i - 1][j - 1];// sum += p[i][j];//\
			psum[j - 1] = sum;// psum[j]=sum;
		}

		for (;;) {
			ilo = 1;
			// First we must determine which point is the highest (worst),
			// next-highest, and lowest
			// (best), by looping over the points in the simplex.
			// ihi = y[1]>y[2] ? (inhi=2,1) : (inhi=1,2);
			// ihi = y[0]>y[1] ? (inhi=2,1) : (inhi=1,2);
			if (y[0] > y[1]) {
				ihi = 1;
				inhi = 2;
			} else {
				ihi = 2;
				inhi = 1;
			}

			for (i = 1; i <= mpts; i++) {
				if (y[i - 1] <= y[ilo - 1])
					ilo = i;// if (y[i] <= y[ilo]) ilo=i;
				if (y[i - 1] > y[ihi - 1])// if (y[i] > y[ihi])
				{
					inhi = ihi;
					ihi = i;
				} else if (y[i - 1] > y[inhi - 1] && i != ihi)
					inhi = i;// (y[i] > y[inhi] && i != ihi) inhi=i;
			}
			// rtol=2.0*Math.abs(y[ihi]-y[ilo])/(fabs(y[ihi])+fabs(y[ilo])+TINY1);
			rtol = 2.0 * Math.abs(y[ihi - 1] - y[ilo - 1])
					/ (Math.abs(y[ihi - 1]) + Math.abs(y[ilo - 1]) + TINY1);
			// Compute the fractional range from highest to lowest and return if
			// satisfactory.
			if (rtol < ftol) {// If returning, put best point and value in slot
								// 1.
								// SWAP(y[1],y[ilo])
				swap = y[0];
				y[0] = y[ilo - 1];
				y[ilo - 1] = swap;
				for (i = 1; i <= ndim; i++) // SWAP(p[1][i],p[ilo][i])
				{
					swap = p[0][i - 1];
					p[0][i - 1] = p[ilo - 1][i - 1];
					p[ilo - 1][i - 1] = swap;
				}
				break;
			}
			if (nfunk_amoeba >= NMAX) {
				// nrerror("NMAX exceeded");
				failB = true;
				failS = "NMAX exceeded";
				return;
			}
			nfunk_amoeba += 2;
			// Begin a new iteration. First extrapolate by a factor -1 through
			// the face of the simplex
			// across from the high point, i.e., reflect the simplex from the
			// high point.
			ytry = amotry(p, y, psum, ndim, ihi, -1.0);// ytry=amotry(p,y,psum,ndim,funk,ihi,-1.0);
			if (ytry <= y[ilo - 1])// if (ytry <= y[ilo])
				// Gives a result better than the best point, so try an
				// additional extrapolation by a factor 2.
				ytry = amotry(p, y, psum, ndim, ihi, 2.0);// ytry=amotry(p,y,psum,ndim,funk,ihi,2.0);
			else if (ytry >= y[inhi - 1])// else if (ytry >= y[inhi])
			{
				// The reflected point is worse than the second-highest, so look
				// for an intermediate
				// lower point, i.e., do a one-dimensional contraction.
				ysave = y[ihi - 1];// ysave=y[ihi];
				ytry = amotry(p, y, psum, ndim, ihi, 0.5);// ytry=amotry(p,y,psum,ndim,funk,ihi,0.5);
				if (ytry >= ysave) {// Can’t seem to get rid of that high point.
									// Better
									// contract around the lowest (best) point.
					for (i = 1; i <= mpts; i++) {
						if (i != ilo) {
							for (j = 1; j <= ndim; j++)
								// p[i][j]=psum[j]=0.5*(p[i][j]+p[ilo][j]);
								p[i - 1][j - 1] = psum[j - 1] = 0.5 * (p[i - 1][j - 1] + p[ilo - 1][j - 1]);
							y[i - 1] = func.MF(psum);// y[i]=(*funk)(psum);
						}
					}
					nfunk_amoeba += ndim;// *nfunk += ndim; Keep track of
											// function evaluations.
					// GET_PSUM Recompute psum.
					for (j = 1; j <= ndim; j++) {// \
						for (sum = 0.0, i = 1; i <= mpts; i++)
							sum += p[i - 1][j - 1];// sum += p[i][j];//\
						psum[j - 1] = sum;// psum[j]=sum;
					}
				}
			} else
				--nfunk_amoeba;// (*nfunk); Correct the evaluation count.
		}// Go back for the test of doneness and the next
			// iteration. free_vector(psum,1,ndim);
	}

	// float amotry(float **p, float y[], float psum[], int ndim,
	// float (*funk)(float []), int ihi, float fac)
	/**
	 * Extrapolates by a factor fac through the face of the simplex across from the high point, tries it, and replaces the high point if the new point is better.
	 * Used internally by amoeba. 
	 * @param p p
	 * @param y y
	 * @param psum psum
	 * @param ndim ndim
	 * @param ihi ihi
	 * @param fac fac
	 * @return the result
	 */
	public double amotry(double[][] p, double[] y, double[] psum, int ndim,
			int ihi, double fac)
	// Extrapolates by a factor fac through the face of the simplex across from
	// the high point, tries
	// it, and replaces the high point if the new point is better.
	{
		int j = 0;
		double fac1 = 0.0;
		double fac2 = 0.0;
		double ytry = 0.0;
		double[] ptry = new double[ndim];
		// ptry=vector(1,ndim);
		fac1 = (1.0 - fac) / ndim;
		fac2 = fac1 - fac;
		for (j = 1; j <= ndim; j++)
			ptry[j - 1] = psum[j - 1] * fac1 - p[ihi - 1][j - 1] * fac2;// ptry[j]=psum[j]*fac1-p[ihi][j]*fac2;
		ytry = func.MF(ptry);// Evaluate the function at the trial point.
		if (ytry < y[ihi - 1])// if (ytry < y[ihi])
		{// If it’s better than the highest, then replace the highest.
			y[ihi - 1] = ytry;// y[ihi]=ytry;
			for (j = 1; j <= ndim; j++) {
				psum[j - 1] += ptry[j - 1] - p[ihi - 1][j - 1];// psum[j] +=
																// ptry[j]-p[ihi][j];
				p[ihi - 1][j - 1] = ptry[j - 1];// p[ihi][j]=ptry[j];
			}
		}
		// free_vector(ptry,1,ndim);
		return ytry;
	}

	/**
	 * Local function.
	 * @param x x
	 * @param name name
	 * @return the result
	 */
	private double localFunc(double x, String name) {
		double res = 0.0;
		if (name.compareTo("f1dim") == 0) {
			res = f1dim(x);
		}
		return res;
	}

	/**
	 * Local function
	 * @param x x
	 * @param name name
	 * @return the result
	 */
	private double dlocalFunc(double x, String name) {
		double res = 0.0;
		if (name.compareTo("df1dim") == 0) {
			res = df1dim(x);
		}
		return res;
	}

	// Direction Set (Powell’s) Methods
	/**
	 * Used internally.
	 * @param x x
	 * @return the result
	 */
	public double f1dim(double x)
	// Must accompany linmin.
	{
		int j = 0;
		double f = 0.0;
		double[] xt = new double[ncom];
		// xt=vector(1,ncom);
		for (j = 1; j <= ncom; j++)
			xt[j - 1] = pcom[j - 1] + x * xicom[j - 1];// xt[j]=pcom[j]+x*xicom[j];
		f = func.MF(xt);// f=(*nrfunc)(xt);
		// free_vector(xt,1,ncom);
		return f;
	}

	/**
	 * Given an n-dimensional point p[1..n] and an n-dimensional direction xi[1..n], moves and 
	 * resets p to where the function func(p) takes on a minimum along the direction xi from p, 
	 * and replaces xi by the actual vector displacement that p was moved. Also returns as fret 
	 * the value of func at the returned location p. This is actually all accomplished by calling the 
	 * routines mnbrak and brent.
	 * @param p p
	 * @param xi xi
	 * @param n n
	 */
	public void linmin(double[] p, double[] xi, int n)// , float *fret, float
														// (*func)(float []))
	// Given an n-dimensional point p[1..n] and an n-dimensional direction
	// xi[1..n], moves and
	// resets p to where the function func(p) takes on a minimum along the
	// direction xi from p,
	// and replaces xi by the actual vector displacement that p was moved. Also
	// returns as fret
	// the value of func at the returned location p. This is actually all
	// accomplished by calling the
	// routines mnbrak and brent.
	{
		// float brent(float ax, float bx, float cx,float (*f)(float),float tol,
		// float *xmin);
		// float f1dim(float x);
		// void mnbrak(float *ax, float *bx, float *cx, float *fa, float
		// *fb,float *fc, float (*func)(float));
		int j = 0;
		double xx = 0.0;
		double xmin = 0.0;
		// double fx=0.0;double fb=0.0;
		// double fa=0.0;
		double bx = 0.0;
		double ax = 0.0;
		ncom = n; // Define the global variables.
		pcom = new double[n];// vector(1,n);
		xicom = new double[n];// vector(1,n);
		// nrfunc=func;
		for (j = 1; j <= n; j++) {
			pcom[j - 1] = p[j - 1];// pcom[j]=p[j];
			xicom[j - 1] = xi[j - 1];// xicom[j]=xi[j];
		}
		ax = 0.0;// Initial guess for brackets.
		xx = 1.0;
		// mnbrak(&ax,&xx,&bx,&fa,&fx,&fb,f1dim);//!!!!!!!!!
		mnbrak1(ax, xx, "f1dim");
		ax = ax_mnbrak;
		xx = bx_mnbrak;
		bx = cx_mnbrak;// fa=fa_mnbrak;fx=fb_mnbrak;fb=fc_mnbrak;
		// *fret=brent(ax,xx,bx,f1dim,TOL,&xmin);//!!!!!!!!!!!!
		fret = brent1(ax, xx, bx, TOL, "f1dim");// ,&xmin);
		xmin = xmin_brent;
		for (j = 1; j <= n; j++) {// Construct the vector results to return.
			xi[j - 1] *= xmin;// xi[j] *= xmin;
			p[j - 1] += xi[j - 1];// p[j] += xi[j];
		}
		// free_vector(xicom,1,n);
		// free_vector(pcom,1,n);
	}

	/**
	 * Minimization of a function func of n variables. Input consists of an initial starting point 
	 * p[1..n]; an initial matrix xi[1..n][1..n], whose columns contain the initial set of directions 
	 * (usually the n unit vectors); and ftol, the fractional tolerance in the function value 
	 * such that failure to decrease by more than this amount on one iteration signals doneness. On 
	 * output, p is set to the best point found, xi is the then-current direction set, fret is the returned 
	 * function value at p, and iter is the number of iterations taken. The routine linmin is used.
	 * @param p p
	 * @param xi xi
	 * @param n n
	 * @param ftol ftol
	 */
	public void powell(double[] p, double[][] xi, int n, double ftol)// , int
																		// *iter,
																		// float
																		// *fret,
	// float (*func)(float []))
	// Minimization of a function func of n variables. Input consists of an
	// initial starting point
	// p[1..n]; an initial matrix xi[1..n][1..n], whose columns contain the
	// initial set of directions
	// (usually the n unit vectors); and ftol, the fractional tolerance in the
	// function value
	// such that failure to decrease by more than this amount on one iteration
	// signals doneness. On
	// output, p is set to the best point found, xi is the then-current
	// direction set, fret is the returned
	// function value at p, and iter is the number of iterations taken. The
	// routine linmin is used.
	{
		// void linmin(float p[], float xi[], int n, float *fret,
		// float (*func)(float []));
		failB = false;

		int i = 0;
		int ibig = 0;
		int j = 0;
		;
		double del = 0.0;
		double fp = 0.0;
		double fptt = 0.0;
		double t = 0.0;
		double[] pt = new double[n];
		double[] ptt = new double[n];
		double[] xit = new double[n];
		// pt=vector(1,n);
		// ptt=vector(1,n);
		// xit=vector(1,n);
		fret = func.MF(p);// *fret=(*func)(p)
		for (j = 1; j <= n; j++)
			pt[j - 1] = p[j - 1];// pt[j]=p[j]; Save the initial point.
		for (iter_powel = 1;; ++(iter_powel)) {
			fp = fret;
			ibig = 0;
			del = 0.0;// Will be the biggest function decrease.
			for (i = 1; i <= n; i++) {// In each iteration, loop over all
										// directions in the set.
				for (j = 1; j <= n; j++)
					xit[j - 1] = xi[j - 1][i - 1];// xit[j]=xi[j][i]; Copy the
													// direction,
				fptt = fret;
				linmin(p, xit, n);// linmin(p,xit,n,fret,func);// minimize along
									// it,
				if (fptt - fret > del) {// and record it if it is the largest
										// decrease so far.
					del = fptt - fret;
					ibig = i;
				}
			}
			if (2.0 * (fp - fret) <= ftol * (Math.abs(fp) + Math.abs(fret))
					+ TINY2) {
				// free_vector(xit,1,n); Termination criterion.
				// free_vector(ptt,1,n);
				// free_vector(pt,1,n);
				return;
			}
			if (iter_powel == ITMAX2) {
				// nrerror("powell exceeding maximum iterations.");
				failB = true;
				failS = "powell exceeding maximum iterations.";
				return;
			}
			for (j = 1; j <= n; j++) {// Construct the extrapolated point and
										// the average direction moved. Save the
										// old starting point.
				ptt[j - 1] = 2.0 * p[j - 1] - pt[j - 1];// ptt[j]=2.0*p[j]-pt[j];
				xit[j - 1] = p[j - 1] - pt[j - 1];// xit[j]=p[j]-pt[j];
				pt[j - 1] = p[j - 1];// pt[j]=p[j];
			}
			fptt = func.MF(ptt); // Function value at extrapolated point.
			if (fptt < fp) {
				// t=2.0*(fp-2.0*(fret)+fptt)*SQR(fp-(*fret)-del)-del*SQR(fp-fptt);
				t = 2.0 * (fp - 2.0 * (fret) + fptt) * (fp - (fret) - del)
						* (fp - (fret) - del) - del * (fp - fptt) * (fp - fptt);
				if (t < 0.0) {
					linmin(p, xit, n);
					// linmin(p,xit,n,fret,func); Move to the minimum of the new
					// direction,
					// and save the new direction.
					for (j = 1; j <= n; j++) {
						xi[j - 1][ibig - 1] = xi[j - 1][n - 1];// xi[j][ibig]=xi[j][n];
						xi[j - 1][n - 1] = xit[j - 1];// xi[j][n]=xit[j];
					}
				}
			}
		}// Back for another iteration.
	}

	/*
	 * Conjugate Gradient Methods in Multidimensions We consider now the case
	 * where you are able to calculate, at a given Ndimensional point P, not
	 * just the value of a function f(P) but also the gradient (vector of first
	 * partial derivatives) ?f(P).
	 * 
	 * The following routine implements the Polak-Ribiere variant, which we
	 * recommend; but changing one program line, as shown, will give you
	 * Fletcher-Reeves. The routine presumes the existence of a function
	 * func(p), where p[1..n] is a vector of length n, and also presumes the
	 * existence of a function dfunc(p,df) that sets the vector gradient
	 * df[1..n] evaluated at the input point p.
	 */
	// public void frprmn(float p[], int n, float ftol, int *iter, float *fret,
	// float (*func)(float []), void (*dfunc)(float [], float []))

	/**
	 * Given a starting point p[1..n], Fletcher-Reeves-Polak-Ribiere minimization is performed on a 
	 * function func, using its gradient as calculated by a routine dfunc. The convergence tolerance 
	 * on the function value is input as ftol. Returned quantities are p (the location of the minimum), 
	 * iter (the number of iterations that were performed), and fret (the minimum value of the 
	 * function). The routine linmin is called to perform line minimizations.
	 * @param p p
	 * @param n n
	 * @param ftol ftol
	 */
	public void frprmn(double[] p, int n, double ftol)// , int *iter, float
														// *fret,
	// float (*func)(float []), void (*dfunc)(float [], float []))

	// Given a starting point p[1..n], Fletcher-Reeves-Polak-Ribiere
	// minimization is performed on a
	// /function func, using its gradient as calculated by a routine dfunc. The
	// convergence tolerance
	// on the function value is input as ftol. Returned quantities are p (the
	// location of the minimum),
	// iter (the number of iterations that were performed), and fret (the
	// minimum value of the
	// function). The routine linmin is called to perform line minimizations.
	{
		failB = false;
		// void linmin(float p[], float xi[], int n, float *fret, float
		// (*func)(float []));
		int j = 0;
		int its = 0;
		double gg = 0.0;
		double gam = 0.0;
		double fp = 0.0;
		double dgg = 0.0;
		double[] g = new double[n];
		double[] h = new double[n];
		double[] xi = new double[n];
		// g=vector(1,n);
		// h=vector(1,n);
		// xi=vector(1,n);
		fp = func.MF(p);// *fret=(*func)(p)//fp=(*func)(p); Initializations.
		xi = func.DMF(p);// (*dfunc)(p,xi);
		for (j = 1; j <= n; j++) {
			g[j - 1] = -xi[j - 1];// g[j] = -xi[j];
			xi[j - 1] = h[j - 1] = g[j - 1];// xi[j]=h[j]=g[j];
		}
		for (its = 1; its <= ITMAX2; its++) {// Loop over iterations.
			iter_fr = its;
			// linmin(p,xi,n,fret,func); Next statement is the normal return:
			linmin(p, xi, n);
			fret_fr = fret;// @@@@@@@@@@@@@@@@@@@@@@@@@@@@
			if (2.0 * Math.abs(fret_fr - fp) <= ftol
					* (Math.abs(fret_fr) + Math.abs(fp) + EPS)) {
				// FREEALL
				return;
			}
			fp = fret_fr;
			xi = func.DMF(p);// (*dfunc)(p,xi);
			dgg = gg = 0.0;
			for (j = 1; j <= n; j++) {
				gg += g[j - 1] * g[j - 1];// gg += g[j]*g[j];
				/* //dgg += xi[j]*xi[j]; */// This statement for
											// Fletcher-Reeves.
				/* dgg += xi[j-1]*xi[j-1]; */// This statement for
												// Fletcher-Reeves.
				// dgg += (xi[j]+g[j])*xi[j]; This statement for Polak-Ribiere.
				dgg += (xi[j - 1] + g[j - 1]) * xi[j - 1]; // This statement for
															// Polak-Ribiere.
			}
			if (gg == 0.0) {// Unlikely. If gradient is exactly zero then we are
							// already done. FREEALL
				return;
			}
			gam = dgg / gg;
			for (j = 1; j <= n; j++) {
				g[j - 1] = -xi[j - 1];// g[j] = -xi[j];
				xi[j - 1] = h[j - 1] = g[j - 1] + gam * h[j - 1];// xi[j]=h[j]=g[j]+gam*h[j];
			}
		}
		// nrerror("Too many iterations in frprmn");
		failB = true;
		failS = "Too many iterations in frprmn";
		return;
	}

	/*
	 * Kindly reread the last part of §10.5. We here want to do the same thing,
	 * but using derivative information in performing the line minimization. The
	 * modified version of linmin, called dlinmin, and its required companion
	 * routine df1dim follow:
	 */
	/**
	 * Given an n-dimensional point p[1..n] and an n-dimensional direction xi[1..n], moves and 
	 * resets p to where the function func(p) takes on a minimum along the direction xi from p, 
	 * and replaces xi by the actual vector displacement that p was moved. Also returns as fret 
	 * the value of func at the returned location p. This is actually all accomplished by calling the 
	 * routines mnbrak and dbrent.
	 * @param p p
	 * @param xi xi
	 * @param n n
	 */
	public void dlinmin(double p[], double xi[], int n)// , float *fret, float
														// (*func)(float []),
	// void (*dfunc)(float [], float []))
	// Given an n-dimensional point p[1..n] and an n-dimensional direction
	// xi[1..n], moves and
	// resets p to where the function func(p) takes on a minimum along the
	// direction xi from p,
	// and replaces xi by the actual vector displacement that p was moved. Also
	// returns as fret
	// the value of func at the returned location p. This is actually all
	// accomplished by calling the
	// routines mnbrak and dbrent.
	{
		// float dbrent(float ax, float bx, float cx, float (*f)(float), float
		// (*df)(float), float tol, float *xmin);
		// float f1dim(float x);
		// float df1dim(float x);
		// void mnbrak(float *ax, float *bx, float *cx, float *fa, float
		// *fb,float *fc, float (*func)(float));

		int j = 0;
		double xx = 0.0;
		double xmin = 0.0;
		// double fx=0.0;
		// double fb=0.0;double fa=0.0;
		double bx = 0.0;
		double ax = 0.0;

		ncom = n;// Define the global variables.
		pcom = new double[n];// vector(1,n);
		xicom = new double[n];// =vector(1,n);
		// nrfunc=func;
		// nrdfun=dfunc;
		for (j = 1; j <= n; j++) {
			pcom[j - 1] = p[j - 1];// pcom[j]=p[j];
			xicom[j - 1] = xi[j - 1];// xicom[j]=xi[j];
		}
		ax = 0.0;// Initial guess for brackets.
		xx = 1.0;

		// mnbrak(&ax,&xx,&bx,&fa,&fx,&fb,f1dim);
		mnbrak1(ax, xx, "f1dim");
		ax = ax_mnbrak;
		xx = bx_mnbrak;
		bx = cx_mnbrak;// fa=fa_mnbrak;fx=fb_mnbrak;fb=fc_mnbrak;

		// *fret=dbrent(ax,xx,bx,f1dim,df1dim,TOL,&xmin);
		fret_dl = dbrent1(ax, xx, bx, TOL, "f1dim", "df1dim");
		xmin = xmin_brent;

		for (j = 1; j <= n; j++) {// Construct the vector results to return.
			xi[j - 1] *= xmin;// xi[j] *= xmin;
			p[j - 1] += xi[j - 1];// p[j] += xi[j];
		}
		// free_vector(xicom,1,n);
		// free_vector(pcom,1,n);
	}

	/**
	 * Used internally
	 * @param x x
	 * @return the result
	 */
	public double df1dim(double x) {
		int j = 0;
		double df1 = 0.0;
		double[] xt = new double[ncom];
		double[] df = new double[ncom];
		// xt=vector(1,ncom);
		// df=vector(1,ncom);
		for (j = 1; j <= ncom; j++)
			xt[j - 1] = pcom[j - 1] + x * xicom[j - 1];// xt[j]=pcom[j]+x*xicom[j];
		// (*nrdfun)(xt,df);
		df = func.DMF(xt);// gradient??
		for (j = 1; j <= ncom; j++)
			df1 += df[j - 1] * xicom[j - 1];// df1 += df[j]*xicom[j];
		// free_vector(df,1,ncom);
		// free_vector(xt,1,ncom);
		return df1;
	}

	/*
	 * Variable Metric Methods in Multidimensions The goal of variable metric
	 * methods, which are sometimes called quasi-Newton methods, is not
	 * different from the goal of conjugate gradient methods: to accumulate
	 * information from successive line minimizations so that N such line
	 * minimizations lead to the exact minimum of a quadratic form in N
	 * dimensions. In that case, the method will also be quadratically
	 * convergent for more general smooth functions. Both variable metric and
	 * conjugate gradient methods require that you are able to compute your
	 * function’s gradient, or first partial derivatives, at arbitrary points. ×
	 * N. Generally, for any moderate N, this is an entirely trivial
	 * disadvantage.
	 * 
	 * On the other hand, there is not, as far as we know, any overwhelming
	 * advantage that the variable metric methods hold over the conjugate
	 * gradient techniques, except perhaps a historical one. Developed somewhat
	 * earlier, and more widely propagated, the variable metric methods have by
	 * now developed a wider constituency of satisfied users. Likewise, some
	 * fancier implementations of variable metric methods (going beyond the
	 * scope of this book, see below) have been developed to a greater level of
	 * sophistication on issues like the minimization of roundoff error,
	 * handling of special conditions, and so on. We tend to use variable metric
	 * rather than conjugate gradient, but we have no reason to urge this habit
	 * on you.
	 */
	/**
	 * Given an n-dimensional point xold[1..n], the value of the function and gradient there, fold 
	 * and g[1..n], and a direction p[1..n], finds a new point x[1..n] along the direction p from 
	 * xold where the function func has decreased “sufficiently.” The new function value is returned 
	 * in f. stpmax is an input quantity that limits the length of the steps so that you do not try to 
	 * evaluate the function in regions where it is undefined or subject to overflow. p is usually the 
	 * Newton direction. The output quantity check is false (0) on a normal exit. It is true (1) when 
	 * x is too close to xold. In a minimization algorithm, this usually signals convergence and can 
	 * be ignored. However, in a zero-finding algorithm the calling program should check whether the 
	 * convergence is spurious. 
	 * @param n n
	 * @param xold xold
	 * @param fold fold
	 * @param g g
	 * @param p p
	 * @param x x
	 * @param stpmax stpmax
	 */
	public void lnsrch(int n, double[] xold, double fold, double[] g,
			double[] p, double[] x, double stpmax)// , String namS)//, int
													// *check, float
													// (*func)(float []))
	// Given an n-dimensional point xold[1..n], the value of the function and
	// gradient there, fold
	// and g[1..n], and a direction p[1..n], finds a new point x[1..n] along the
	// direction p from
	// xold where the function func has decreased “sufficiently.” The new
	// function value is returned
	// in f. stpmax is an input quantity that limits the length of the steps so
	// that you do not try to
	// evaluate the function in regions where it is undefined or subject to
	// overflow. p is usually the
	// Newton direction. The output quantity check is false (0) on a normal
	// exit. It is true (1) when
	// x is too close to xold. In a minimization algorithm, this usually signals
	// convergence and can
	// be ignored. However, in a zero-finding algorithm the calling program
	// should check whether the
	// convergence is spurious. Some “difficult” problems may require double
	// precision in this routine.
	{
		int i = 0;
		double a = 0.0;
		double alam = 0.0;
		double alam2 = 0.0;
		double alamin = 0.0;
		double b = 0.0;
		double disc = 0.0;
		double f2 = 0.0;
		double rhs1 = 0.0;
		double rhs2 = 0.0;
		double slope = 0.0;
		double sum = 0.0;
		double temp = 0.0;
		double test = 0.0;
		double tmplam = 0.0;

		failB = false;

		check_ln = 0;
		for (sum = 0.0, i = 1; i <= n; i++)
			sum += p[i - 1] * p[i - 1];// sum += p[i]*p[i];
		sum = Math.sqrt(sum);
		if (sum > stpmax)
			for (i = 1; i <= n; i++)
				p[i - 1] *= stpmax / sum;// p[i] *= stpmax/sum; Scale if
											// attempted step is too big.
		for (slope = 0.0, i = 1; i <= n; i++)
			slope += g[i - 1] * p[i - 1];// slope += g[i]*p[i];
		if (slope >= 0.0) {
			// nrerror("Roundoff problem in lnsrch.");
			failB = true;
			failS = "Roundoff problem in lnsrch.";
			return;
		}
		test = 0.0; // Compute ?min.
		for (i = 1; i <= n; i++) {
			temp = Math.abs(p[i - 1]) / Math.max(Math.abs(xold[i - 1]), 1.0);// fabs(p[i])/FMAX(fabs(xold[i]),1.0);
			if (temp > test)
				test = temp;
		}
		alamin = TOLX / test;
		alam = 1.0;// Always try full Newton step first.
		for (;;) {// Start of iteration loop.
			for (i = 1; i <= n; i++)
				x[i - 1] = xold[i - 1] + alam * p[i - 1];// x[i]=xold[i]+alam*p[i];
			f_ln = func.MF(x);// localFunc(x,namS);//func.vecfunc(n,x);//*f=(*func)(x);
			if (alam < alamin) {// Convergence on ?x. For zero finding,the
								// calling program should verify the
								// convergence.
				for (i = 1; i <= n; i++)
					x[i - 1] = xold[i - 1];// x[i]=xold[i];
				check_ln = 1;
				return;
			} else if (f_ln <= fold + ALF * alam * slope)
				return; // Sufficient function decrease.
			else {// Backtrack.
				if (alam == 1.0)
					tmplam = -slope / (2.0 * (f_ln - fold - slope)); // First
																		// time.
				else {// Subsequent backtracks.
					rhs1 = f_ln - fold - alam * slope;
					rhs2 = f2 - fold - alam2 * slope;
					a = (rhs1 / (alam * alam) - rhs2 / (alam2 * alam2))
							/ (alam - alam2);
					b = (-alam2 * rhs1 / (alam * alam) + alam * rhs2
							/ (alam2 * alam2))
							/ (alam - alam2);
					if (a == 0.0)
						tmplam = -slope / (2.0 * b);
					else {
						disc = b * b - 3.0 * a * slope;
						if (disc < 0.0)
							tmplam = 0.5 * alam;
						else if (b <= 0.0)
							tmplam = (-b + Math.sqrt(disc)) / (3.0 * a);
						else
							tmplam = -slope / (b + Math.sqrt(disc));
					}
					if (tmplam > 0.5 * alam)
						tmplam = 0.5 * alam; // ? = 0.5?1.
				}
			}
			alam2 = alam;
			f2 = f_ln;
			alam = Math.max(tmplam, 0.1 * alam); // ? = 0.1?1.
		}// Try again.
	}

	/**
	 * Given a starting point p[1..n] that is a vector of length n, the Broyden-Fletcher-Goldfarb-
	 * Shanno variant of Davidon-Fletcher-Powell minimization is performed on a function func, using 
	 * its gradient as calculated by a routine dfunc. The convergence requirement on zeroing the 
	 * gradient is input as gtol. Returned quantities are p[1..n] (the location of the minimum), 
	 * iter (the number of iterations that were performed), and fret (the minimum value of the 
	 * function). The routine lnsrch is called to perform approximate line minimizations.
	 * @param p p
	 * @param n n
	 * @param gtol gtol
	 */
	public void dfpmin(double[] p, int n, double gtol)// , int *iter, float
														// *fret,
	// float(*func)(float []), void (*dfunc)(float [], float []))
	// Given a starting point p[1..n] that is a vector of length n, the
	// Broyden-Fletcher-Goldfarb-
	// Shanno variant of Davidon-Fletcher-Powell minimization is performed on a
	// function func, using
	// its gradient as calculated by a routine dfunc. The convergence
	// requirement on zeroing the
	// gradient is input as gtol. Returned quantities are p[1..n] (the location
	// of the minimum),
	// iter (the number of iterations that were performed), and fret (the
	// minimum value of the
	// function). The routine lnsrch is called to perform approximate line
	// minimizations.
	{
		failB = false;
		// void lnsrch(int n, float xold[], float fold, float g[], float p[],
		// float x[], float *f, float stpmax, int *check, float (*func)(float
		// []));
		// int check=0;
		int i = 0;
		int its = 0;
		int j = 0;
		double den = 0.0;
		double fac = 0.0;
		double fad = 0.0;
		double fae = 0.0;
		double fp = 0.0;
		double stpmax = 0.0;
		double sum = 0.0;
		double sumdg = 0.0;
		double sumxi = 0.0;
		double temp = 0.0;
		double test = 0.0;

		double[] dg = new double[n];
		double[] g = new double[n];
		double[] hdg = new double[n];
		double[][] hessin = new double[n][n];
		double[] pnew = new double[n];
		double[] xi = new double[n];
		// dg=vector(1,n);
		// g=vector(1,n);
		// hdg=vector(1,n);
		// hessin=matrix(1,n,1,n);
		// pnew=vector(1,n);
		// xi=vector(1,n);
		fp = func.MF(p);
		// fp=(*func)(p); Calculate starting function value and gradient,
		g = func.DMF(p);// (*dfunc)(p,g);

		for (i = 1; i <= n; i++) { // and initialize the inverse Hessian to the
									// unit matrix.
			for (j = 1; j <= n; j++)
				hessin[i - 1][j - 1] = 0.0;// hessin[i][j]=0.0;
			hessin[i - 1][i - 1] = 1.0;// hessin[i][i]=1.0;
			xi[i - 1] = -g[i - 1];// xi[i] = -g[i]; Initial line direction.
			sum += p[i - 1] * p[i - 1];// sum += p[i]*p[i];
		}
		stpmax = STPMX * Math.max(Math.sqrt(sum), (double) n);
		for (its = 1; its <= ITMAX2; its++) {// Main loop over the iterations.
			iter_df = its;
			// lnsrch(n,p,fp,g,xi,pnew,fret,stpmax,&check,func);
			lnsrch(n, p, fp, g, xi, pnew, stpmax);
			fret_df = f_ln;
			// The new function evaluation occurs in lnsrch; save the function
			// value in fp for the
			// next line search. It is usually safe to ignore the value of
			// check.
			fp = fret_df;
			for (i = 1; i <= n; i++) {
				xi[i - 1] = pnew[i - 1] - p[i - 1];// xi[i]=pnew[i]-p[i]; Update
													// the line direction,
				p[i - 1] = pnew[i - 1];// p[i]=pnew[i]; and the current point.
			}
			test = 0.0; // Test for convergence on ?x.
			for (i = 1; i <= n; i++) {
				// temp=fabs(xi[i])/FMAX(fabs(p[i]),1.0);
				temp = Math.abs(xi[i - 1]) / Math.max(Math.abs(p[i - 1]), 1.0);
				if (temp > test)
					test = temp;
			}
			if (test < TOLX) {
				// FREEALL
				return;
			}
			for (i = 1; i <= n; i++)
				dg[i - 1] = g[i - 1];// dg[i]=g[i]; Save the old gradient,
			// (*dfunc)(p,g); and get the new gradient.
			g = func.DMF(p);
			test = 0.0; // Test for convergence on zero gradient.
			den = Math.max(fret_df, 1.0);
			for (i = 1; i <= n; i++) {
				// temp=fabs(g[i])*FMAX(fabs(p[i]),1.0)/den;
				temp = Math.abs(g[i - 1]) * Math.max(Math.abs(p[i - 1]), 1.0)
						/ den;
				if (temp > test)
					test = temp;
			}
			if (test < gtol) {
				// FREEALL
				return;
			}
			for (i = 1; i <= n; i++)
				dg[i - 1] = g[i - 1] - dg[i - 1];// dg[i]=g[i]-dg[i]; Compute
													// difference of gradients,
			for (i = 1; i <= n; i++) {// and difference times current matrix.
				hdg[i - 1] = 0.0;// hdg[i]=0.0;
				for (j = 1; j <= n; j++)
					hdg[i - 1] += hessin[i - 1][j - 1] * dg[j - 1];// hdg[i] +=
																	// hessin[i][j]*dg[j];
			}
			fac = fae = sumdg = sumxi = 0.0;// Calculate dot products for the
											// denominators.
			for (i = 1; i <= n; i++) {
				fac += dg[i - 1] * xi[i - 1];// fac += dg[i]*xi[i];
				fae += dg[i - 1] * hdg[i - 1];// fae += dg[i]*hdg[i];
				sumdg += dg[i - 1] * dg[i - 1];// sumdg += SQR(dg[i]);
				sumxi += xi[i - 1] * xi[i - 1];// sumxi += SQR(xi[i]);
			}
			if (fac > Math.sqrt(EPS2 * sumdg * sumxi)) {// Skip update if fac
														// not sufficiently
														// positive.
				fac = 1.0 / fac;
				fad = 1.0 / fae;
				// The vector that makes BFGS different from DFP:
				for (i = 1; i <= n; i++)
					dg[i - 1] = fac * xi[i - 1] - fad * hdg[i - 1];// dg[i]=fac*xi[i]-fad*hdg[i];
				for (i = 1; i <= n; i++) {// The BFGS updating formula:
					for (j = i; j <= n; j++) {
						// hessin[i][j] += fac*xi[i]*xi[j]
						hessin[i - 1][j - 1] += fac * xi[i - 1] * xi[j - 1]
								// -fad*hdg[i]*hdg[j]+fae*dg[i]*dg[j];
								- fad * hdg[i - 1] * hdg[j - 1] + fae
								* dg[i - 1] * dg[j - 1];
						hessin[j - 1][i - 1] = hessin[i - 1][j - 1];// hessin[j][i]=hessin[i][j];
					}
				}
			}
			for (i = 1; i <= n; i++) {// Now calculate the next direction to go,
				xi[i - 1] = 0.0;// xi[i]=0.0;
				for (j = 1; j <= n; j++)
					xi[i - 1] -= hessin[i - 1][j - 1] * g[j - 1];// xi[i] -=
																	// hessin[i][j]*g[j];
			}
		}// and go back for another iteration.
		failB = true;
		failS = "too many iterations in dfpmin";
		return;
		// nrerror("too many iterations in dfpmin");
		// FREEALL
	}

	/*
	 * Linear Programming and the Simplex
	 * 
	 * The subject of linear programming, sometimes called linear optimization,
	 * concerns itself with the following problem: ForN independent variables
	 * x1, . . . , xN, maximize the function z = a01x1 + a02x2 + · · · + a0NxN
	 * (10.8.1) subject to the primary constraints x1 >= 0, x2 >= 0, . . . xN
	 * >=0
	 * 
	 * and simultaneously subject to M = m1 + m2 + m3 additional constraints, m1
	 * of them of the form ai1x1 + ai2x2 + · · · + aiNxN = bi (bi = 0) i = 1, .
	 * . .,m1 (10.8.3) m2 of them of the form aj1x1 + aj2x2 + · · · + ajNxN = bj
	 * = 0 j = m1 +1, . . .,m1 +m2 (10.8.4) and m3 of them of the form ak1x1 +
	 * ak2x2 + · · · + akNxN = bk = 0 k = m1 + m2 + 1, . . .,m1 + m2 + m3
	 * 
	 * Here is a specific example of a problem in linear programming, which has
	 * N = 4, m1 = 2, m2 = m3 = 1, hence M = 4: Maximize z = x1 + x2 + 3x3 -
	 * 1/2*x4
	 * 
	 * The following routine is based algorithmically on the implementation
	 * ofKuenzi, Tzschach, and Zehnder [4]. Aside from input values of M, N, m1,
	 * m2, m3, the principal input to the routine is a two-dimensional array a
	 * containing the portion of the tableau (10.8.18) that is contained between
	 * the double lines. This input occupies the M + 1 rows and N + 1 columns of
	 * a[1..m+1][1..n+1]. Note, however, that reference is made internally to
	 * row M + 2 of a (used for the auxiliary objective function, just as in
	 * 10.8.18). Therefore the variable declared as float **a, must point to
	 * allocated memory allowing references in the subrange a[i][k], i= 1. . .
	 * m+2, k = 1. . . n+1
	 * 
	 * #include "nrutil.h" #define EPS 1.0e-6 Here EPS is the absolute
	 * precision, which should be adjusted to the scale of your variables.
	 * #define FREEALL free_ivector(l3,1,m);free_ivector(l1,1,n+1); void
	 * simplx(float **a, int m, int n, int m1, int m2, int m3, int *icase, int
	 * izrov[], int iposv[]) Simplex method for linear programming. Input
	 * parameters a, m, n, mp, np, m1, m2, and m3, and output parameters a,
	 * icase, izrov, and iposv are described above. { void simp1(float **a, int
	 * mm, int ll[], int nll, int iabf, int *kp, float *bmax); void simp2(float
	 * **a, int m, int n, int *ip, int kp); void simp3(float **a, int i1, int
	 * k1, int ip, int kp); int i,ip,is,k,kh,kp,nl1; int *l1,*l3; float q1,bmax;
	 * if (m != (m1+m2+m3)) nrerror("Bad input constraint counts in simplx");
	 * l1=ivector(1,n+1); l3=ivector(1,m); nl1=n; for (k=1;k<=n;k++)
	 * l1[k]=izrov[k]=k; Initialize index list of columns admissible for
	 * exchange, and make all variables initially right-hand. for (i=1;i<=m;i++)
	 * { if (a[i+1][1] < 0.0) nrerror("Bad input tableau in simplx"); Constants
	 * bi must be nonnegative. iposv[i]=n+i; Initial left-hand variables. m1
	 * type constraints are represented by having their slack variable initially
	 * left-hand, with no artificial variable. m2 type constraints have their
	 * slack variable initially left-hand, with a minus sign, and their
	 * artificial variable handled implicitly during their first exchange. m3
	 * type constraints have their artificial variable initially left-hand. } if
	 * (m2+m3) { Origin is not a feasible starting solution: we must do phase
	 * one. for (i=1;i<=m2;i++) l3[i]=1; Initialize list of m2 constraints whose
	 * slack variables have never been exchanged out of the initial basis. for
	 * (k=1;k<=(n+1);k++) { Compute the auxiliary objective function. q1=0.0;
	 * for (i=m1+1;i<=m;i++) q1 += a[i+1][k]; a[m+2][k] = -q1; } for (;;) {
	 * simp1(a,m+1,l1,nl1,0,&kp,&bmax); Find max. coeff. of auxiliary objective
	 * fn. if (bmax <= EPS && a[m+2][1] < -EPS) {icase = -1; Auxiliary objective
	 * function is still negative and can’t be improved, hence no feasible
	 * solution exists. FREEALL return; } else if (bmax <= EPS && a[m+2][1] <=
	 * EPS) { Auxiliary objective function is zero and can’t be improved; we
	 * have a feasible starting vector. Clean out the artificial variables
	 * corresponding to any remaining equality constraints by goto one and then
	 * move on to phase two. for (ip=m1+m2+1;ip<=m;ip++) { if (iposv[ip] ==
	 * (ip+n)) { Found an artificial variable for an equality constraint.
	 * simp1(a,ip,l1,nl1,1,&kp,&bmax); if (bmax > EPS) Exchange with column
	 * corresponding to maximum pivot element in row. goto one; } } for
	 * (i=m1+1;i<=m1+m2;i++) Change sign of row for any m2 constraints still
	 * present from the initial basis. if (l3[i-m1] == 1) for (k=1;k<=n+1;k++)
	 * a[i+1][k] = -a[i+1][k]; break; Go to phase two. } simp2(a,m,n,&ip,kp);
	 * Locate a pivot element (phase one). if (ip == 0) { Maximum of auxiliary
	 * objective function is unbounded, so no feasible solution exists.icase =
	 * -1; FREEALL return; } one: simp3(a,m+1,n,ip,kp); Exchange a left- and a
	 * right-hand variable (phase one), then update lists. if (iposv[ip] >=
	 * (n+m1+m2+1)) { Exchanged out an artificial variable for an equality
	 * constraint. Make sure it stays out by removing it from the l1 list. for
	 * (k=1;k<=nl1;k++) if (l1[k] == kp) break; --nl1; for (is=k;is<=nl1;is++)
	 * l1[is]=l1[is+1]; } else { kh=iposv[ip]-m1-n; if (kh >= 1 && l3[kh]) {
	 * Exchanged out an m2 type constraint for the first time. Correct the pivot
	 * column for the minus sign and the implicit artificial variable. l3[kh]=0;
	 * ++a[m+2][kp+1]; for (i=1;i<=m+2;i++) a[i][kp+1] = -a[i][kp+1]; } }
	 * is=izrov[kp]; Update lists of left- and right-hand variables.
	 * izrov[kp]=iposv[ip]; iposv[ip]=is; } Still in phase one, go back to the
	 * for(;;). } End of phase one code for finding an initial feasible
	 * solution. Now, in phase two, optimize it. for (;;) {
	 * simp1(a,0,l1,nl1,0,&kp,&bmax); Test the z-row for doneness. if (bmax <=
	 * EPS) { Done. Solution found. Return with the good news. *icase=0; FREEALL
	 * return; } simp2(a,m,n,&ip,kp); Locate a pivot element (phase two). if (ip
	 * == 0) { Objective function is unbounded. Report and return. *icase=1;
	 * FREEALL return; } simp3(a,m,n,ip,kp); Exchange a left- and a right-hand
	 * variable (phase two), is=izrov[kp]; izrov[kp]=iposv[ip]; iposv[ip]=is; }
	 * and return for another iteration. }
	 * 
	 * 
	 * 
	 * The preceding routine makes use of the following utility functions. void
	 * simp1(float **a, int mm, int ll[], int nll, int iabf, int *kp, float
	 * *bmax) Determines the maximum of those elements whose index is contained
	 * in the supplied list ll, either with or without taking the absolute
	 * value, as flagged by iabf. { int k; float test; if (nll <= 0) No eligible
	 * columns.bmax=0.0; else {kp=ll[1];bmax=a[mm+1][*kp+1]; for
	 * (k=2;k<=nll;k++) { if (iabf == 0) test=a[mm+1][ll[k]+1]-(*bmax); else
	 * test=fabs(a[mm+1][ll[k]+1])-fabs(*bmax); if (test > 0.0) {
	 * bmax=a[mm+1][ll[k]+1];kp=ll[k]; } } } }
	 * 
	 * 
	 * #define EPS 1.0e-6 void simp2(float **a, int m, int n, int *ip, int kp)
	 * Locate a pivot element, taking degeneracy into account. { int k,i; float
	 * qp,q0,q,q1;ip=0; for (i=1;i<=m;i++) if (a[i+1][kp+1] < -EPS) break; Any
	 * possible pivots? if (i>m) return; q1 = -a[i+1][1]/a[i+1][kp+1];ip=i; for
	 * (i=*ip+1;i<=m;i++) { if (a[i+1][kp+1] < -EPS) { q =
	 * -a[i+1][1]/a[i+1][kp+1]; if (q < q1) {ip=i; q1=q; } else if (q == q1) {
	 * We have a degeneracy. for (k=1;k<=n;k++) { qp =
	 * -a[*ip+1][k+1]/a[*ip+1][kp+1]; q0 = -a[i+1][k+1]/a[i+1][kp+1]; if (q0 !=
	 * qp) break; } if (q0 < qp) *ip=i; } } } }
	 * 
	 * 
	 * void simp3(float **a, int i1, int k1, int ip, int kp) Matrix operations
	 * to exchange a left-hand and right-hand variable (see text). { int kk,ii;
	 * float piv; piv=1.0/a[ip+1][kp+1]; for (ii=1;ii<=i1+1;ii++) if (ii-1 !=
	 * ip) { a[ii][kp+1] *= piv; for (kk=1;kk<=k1+1;kk++) if (kk-1 != kp)
	 * a[ii][kk] -= a[ip+1][kk]*a[ii][kp+1]; } for (kk=1;kk<=k1+1;kk++) if (kk-1
	 * != kp) a[ip+1][kk] *= -piv; a[ip+1][kp+1]=piv; }
	 * 
	 * 
	 * 
	 * A concrete illustration is provided by the traveling salesman problem.
	 * The proverbial seller visits N cities with given positions (xi, yi),
	 * returning finally to his or her city of origin. Each city is to be
	 * visited only once, and the route is to be made as short as possible. This
	 * problem belongs to a class known as NP-complete problems, whose
	 * computation time for an exact solution increases with N as exp(const.×
	 * N), becoming rapidly prohibitive in cost as N increases. The traveling
	 * salesman problem also belongs to a class of minimization problems for
	 * which the objective function E has many local minima. In practical cases,
	 * it is often enough to be able to choose from these a minimum which, even
	 * if not absolute, cannot be significantly improved upon. The annealing
	 * method manages to achieve this, while limiting its calculations to scale
	 * as a small power of N.
	 * 
	 * ANNEALING=>NOT IMPLEMENTED!!
	 */
}
