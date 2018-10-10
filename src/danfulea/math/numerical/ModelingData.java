package danfulea.math.numerical;

/**
 * Modeling of data<br>
 * Based on literature, including Numerical Recipes in C (Cambridge Univ.).
 * @author Dan Fulea, 18 OCT. 2006
 */
public class ModelingData {
	public static Function func;
	public static boolean failB = false;
	public static String failS = "";

	public static double a_fit = 0.0;
	public static double b_fit = 0.0;
	public static double siga_fit = 0.0;
	public static double sigb_fit = 0.0;
	public static double chi2_fit = 0.0;
	public static double q_fit = 0.0;

	public static double POTN = 1.571000;
	public static double BIG = 1.0e30;
	public static double PI = 3.14159265;
	public static double ACC = 1.0e-3;
	public static int nn_fitexy = 0;
	public static double[] xx_fitexy;
	public static double[] yy_fitexy;
	public static double[] sx_fitexy;
	public static double[] sy_fitexy;
	public static double[] ww_fitexy;
	public static double aa_fitexy = 0.0;
	public static double offs_fitexy = 0.0;
	public static double a_fitexy = 0.0;
	public static double b_fitexy = 0.0;
	public static double siga_fitexy = 0.0;
	public static double sigb_fitexy = 0.0;
	public static double chi2_fitexy = 0.0;
	public static double q_fitexy = 0.0;
	public static double GOLD = 1.618034;
	public static double GLIMIT = 100.0;
	public static double TINY = 1.0e-20;
	public static double ax_mnbrak = 0.0;
	public static double bx_mnbrak = 0.0;
	public static double cx_mnbrak = 0.0;
	public static double fa_mnbrak = 0.0;
	public static double fb_mnbrak = 0.0;
	public static double fc_mnbrak = 0.0;
	public static int ITMAX = 100;
	public static double CGOLD = 0.3819660;
	public static double ZEPS = 1.0e-10;
	public static double xmin_brent = 0.0;
	public static double EPS = 3.0e-8;

	public static double chisq_lfit = 0.0;

	public static double TOL = 1.0e-5;
	public static double chisq_svdfit = 0.0;
	public static double chisq_mrqmin = 0.0;
	public static int mfit_mrqmin = 0;
	public static double alamda_mrqmin = 0;
	public static double ochisq_mrqmin = 0.0;
	public static double[] atry_mrqmin;
	public static double[] beta_mrqmin;
	public static double[] da_mrqmin;
	public static double[][] oneda_mrqmin;
	public static boolean convB = false;

	public static int ndatat_medfit = 0;
	public static double[] xt_medfit;
	public static double[] yt_medfit;
	public static double aa_medfit = 0.0;
	public static double abdevt_medfit = 0.0;
	public static double EPS2 = 1.0e-7;
	public static double a_medfit = 0.0;
	public static double b_medfit = 0.0;
	public static double abdev_medfit = 0.0;

	/*
	 * Given a set of observations, one often wants to condense and summarize
	 * the data by fitting it to a “model” that depends on adjustable
	 * parameters. Sometimes the model is simply a convenient class of
	 * functions, such as polynomials or Gaussians, and the fit supplies the
	 * appropriate coefficients. Other times, the model’s parameters come from
	 * some underlying theory that the data are supposed to satisfy; examples
	 * are coefficients of rate equations in a complex network of chemical
	 * reactions, or orbital elements of a binary star. Modeling can also be
	 * used as a kind of constrained interpolation, where you want to extend a
	 * few data points into a continuous function, but with some underlying idea
	 * of what that function should look like. The basic approach in all cases
	 * is usually the same: You choose or design a figure-of-merit function
	 * (“merit function,” for short) that measures the agreement between the
	 * data and the model with a particular choice of parameters. The merit
	 * function is conventionally arranged so that small values represent close
	 * agreement. The parameters of the model are then adjusted to achieve a
	 * minimum in the merit function, yielding best-fit parameters. The
	 * adjustment process is thus a problem in minimization in many dimensions.
	 * This optimization was the subject of Chapter 10; however, there exist
	 * special, more efficient, methods that are specific to modeling, and we
	 * will discuss these in this chapter. There are important issues that go
	 * beyond the mere finding of best-fit parameters. Data are generally not
	 * exact. They are subject to measurement errors (called noise in the
	 * context of signal-processing). Thus, typical data never exactly fit the
	 * model that is being used, even when that model is correct. We need the
	 * means to assess whether or not the model is appropriate, that is, we need
	 * to test the goodness-of-fit against some useful statistical standard. We
	 * usually also need to knowthe accuracy with which parameters are
	 * determined by the data set. In other words, we need to know the likely
	 * errors of the best-fit parameters. Finally, it is not uncommon in fitting
	 * data to discover that the merit function is not unimodal, with a single
	 * minimum. In some cases, we may be interested in global rather than local
	 * questions. Not, “how good is this fit?” but rather, “how sure am I that
	 * there is not a very much better fit in some corner of parameter space?”
	 * As we have seen in Chapter 10, especially §10.9, this kind of problem is
	 * generally quite difficult to solve.
	 * 
	 * Least Squares as a Maximum Likelihood Estimator y(x) = y(x; a1 . . .
	 * aM)=> minimize over a1 . . . aM :SUM i=1,N of [yi - y(xi; a1 . . . aM)]^2
	 * 
	 * Chi-Square Fitting chi2 =SUM i=1,N of {[yi - y(xi; a1 . . . aM)]/sigma i
	 * }^2
	 * 
	 * Fitting Data to a Straight Line
	 */
	
	/**
	 * Change the sign of a if b is negative!
	 * @param a a
	 * @param b b
	 * @return a or -a
	 */
	public static double SIGN(double a, double b) {
		return b >= 0.0 ? a : -a;
	}

	/**
	 * Given a set of data points x[1..ndata],y[1..ndata] with individual standard deviations sig[1..ndata], fit them to a straight 
	 * line y = a + bx by minimizing chi2. Returned are a_fit,b_fit and their respective probable uncertainties siga_fit and sigb_fit, the 
	 * chi-square chi2_fit, and the goodness-of-fit probability q_fit (that the fit would have chi2 this large or larger). If mwt=0 on 
	 * input, then the standard deviations are assumed to be unavailable: q_fit is returned as 1.0 and the normalization of chi2_fit is 
	 * to unit standard deviation on all points.
	 * @param x x
	 * @param y y
	 * @param ndata ndata
	 * @param sig sig
	 * @param mwt mwt
	 */
	public static void fit(double[] x, double[] y, int ndata, double[] sig,
			int mwt)
	// , float *a, float *b, float *siga, float *sigb, float *chi2, float *q)
	// Given a set of data points x[1..ndata],y[1..ndata] with individual
	// standard deviations
	// sig[1..ndata], fit them to a straight line y = a + bx by minimizing chi2.
	// Returned are
	// a,b and their respective probable uncertainties siga and sigb, the
	// chi-square chi2, and the
	// goodness-of-fit probability q (that the fit would have chi2 this large or
	// larger). If mwt=0 on
	// input, then the standard deviations are assumed to be unavailable: q is
	// returned as 1.0 and
	// the normalization of chi2 is to unit standard deviation on all points.
	{
		// float gammq(float a, float x);
		int i = 0;
		double wt = 0.0;
		double t = 0.0;
		double sxoss = 0.0;
		double sx = 0.0;
		double sy = 0.0;
		double st2 = 0.0;
		double ss = 0.0;
		double sigdat = 0.0;

		b_fit = 0.0;
		if (mwt != 0)// if (mwt)
		{// Accumulate sums ...
			ss = 0.0;
			for (i = 1; i <= ndata; i++) {// ...with weights
				wt = 1.0 / ((sig[i - 1]) * (sig[i - 1]));// wt=1.0/SQR(sig[i]);
				ss += wt;
				sx += x[i - 1] * wt;// sx += x[i]*wt;
				sy += y[i - 1] * wt;// sy += y[i]*wt;
			}
		} else {
			for (i = 1; i <= ndata; i++) {// ...or without weights.
				sx += x[i - 1];// sx += x[i];
				sy += y[i - 1];// sy += y[i];
			}
			ss = ndata;
		}
		sxoss = sx / ss;
		if (mwt != 0)// if (mwt)
		{
			for (i = 1; i <= ndata; i++) {
				t = (x[i - 1] - sxoss) / sig[i - 1];// t=(x[i]-sxoss)/sig[i];
				st2 += t * t;
				b_fit += t * y[i - 1] / sig[i - 1];// b_fit += t*y[i]/sig[i];
			}
		} else {
			for (i = 1; i <= ndata; i++) {
				t = x[i - 1] - sxoss;// t=x[i]-sxoss;
				st2 += t * t;
				b_fit += t * y[i - 1];// b_fit += t*y[i];
			}
		}
		b_fit /= st2;// Solve for a, b, ?a, and ?b.
		a_fit = (sy - sx * (b_fit)) / ss;
		siga_fit = Math.sqrt((1.0 + sx * sx / (ss * st2)) / ss);
		sigb_fit = Math.sqrt(1.0 / st2);
		chi2_fit = 0.0;// Calculate ?2.
		q_fit = 1.0;
		if (mwt == 0) {
			for (i = 1; i <= ndata; i++)
				// chi2_fit += SQR(y[i]-(a_fit)-(b_fit)*x[i]);
				chi2_fit += (y[i - 1] - (a_fit) - (b_fit) * x[i - 1])
						* (y[i - 1] - (a_fit) - (b_fit) * x[i - 1]);
			sigdat = Math.sqrt((chi2_fit) / (ndata - 2)); // For unweighted data
															// evaluate typical
			// sig using chi2, and adjust the standard deviations.
			siga_fit *= sigdat;
			sigb_fit *= sigdat;
		} else {
			for (i = 1; i <= ndata; i++)
				// chi2_fit += SQR((y[i]-(a_fit)-(b_fit)*x[i])/sig[i]);
				chi2_fit += ((y[i - 1] - (a_fit) - (b_fit) * x[i - 1]) / sig[i - 1])
						* ((y[i - 1] - (a_fit) - (b_fit) * x[i - 1]) / sig[i - 1]);
			if (ndata > 2)
				q_fit = SpecFunc.gammq(0.5 * (ndata - 2), 0.5 * (chi2_fit));// Equation
																			// (15.2.12).
		}
	}

	/*
	 * Straight-Line Data with Errors in Both Coordinates
	 */
	/**
	 * Using Brent’s method, find the root of a function func known to lie between x1 and x2. The 
	 * root, returned as zbrent, will be refined until its accuracy is tol.
	 * @param x1 x1
	 * @param x2 x2
	 * @param tol tol
	 * @param namS the name of local function to be used.
	 * @return the root
	 */
	public static double zbrent1(double x1, double x2, double tol, String namS)
	// Using Brent’s method, find the root of a function func known to lie
	// between x1 and x2. The
	// root, returned as zbrent, will be refined until its accuracy is tol.
	{
		failB = false;
		int iter = 0;
		double a = x1;
		double b = x2;
		double c = x2;
		double d = 0.0;
		double e = 0.0;
		double min1 = 0.0;
		double min2 = 0.0;

		double fa = localFunc(a, namS);
		double fb = localFunc(b, namS);
		double fc = 0.0;
		double p = 0.0;
		double q = 0.0;
		double r = 0.0;
		double s = 0.0;
		double tol1 = 0.0;
		double xm = 0.0;

		if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0)) {
			// nrerror("Root must be bracketed in zbrent");
			failB = true;
			failS = "Root must be bracketed in zbrent";
			return 0.0;
		}
		fc = fb;
		for (iter = 1; iter <= ITMAX; iter++)// ->100 ok
		{
			if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) {
				c = a;// Rename a, b, c and adjust bounding interval d.
				fc = fa;
				e = d = b - a;
			}
			if (Math.abs(fc) < Math.abs(fb)) {
				a = b;
				b = c;
				c = a;
				fa = fb;
				fb = fc;
				fc = fa;
			}
			tol1 = 2.0 * EPS * Math.abs(b) + 0.5 * tol;// Convergence check.
			xm = 0.5 * (c - b);
			if (Math.abs(xm) <= tol1 || fb == 0.0)
				return b;
			if (Math.abs(e) >= tol1 && Math.abs(fa) > Math.abs(fb)) {
				s = fb / fa;// Attempt inverse quadratic interpolation.
				if (a == c) {
					p = 2.0 * xm * s;
					q = 1.0 - s;
				} else {
					q = fa / fc;
					r = fb / fc;
					p = s * (2.0 * xm * q * (q - r) - (b - a) * (r - 1.0));
					q = (q - 1.0) * (r - 1.0) * (s - 1.0);
				}
				if (p > 0.0)
					q = -q; // Check whether in bounds.
				p = Math.abs(p);
				min1 = 3.0 * xm * q - Math.abs(tol1 * q);
				min2 = Math.abs(e * q);
				if (2.0 * p < (min1 < min2 ? min1 : min2)) {
					e = d;// Accept interpolation.
					d = p / q;
				} else {
					d = xm;// Interpolation failed, use bisection.
					e = d;
				}
			} else {// Bounds decreasing too slowly, use bisection.
				d = xm;
				e = d;
			}
			a = b;// Move last best guess to a.
			fa = fb;
			if (Math.abs(d) > tol1) // Evaluate new trial root.
				b += d;
			else
				b += SIGN(tol1, xm);
			fb = localFunc(b, namS);
		}
		// nrerror("Maximum number of iterations exceeded in zbrent");
		failB = true;
		failS = "Maximum number of iterations exceeded in zbrent";

		return 0.0; // Never get here.
	}

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
	 * @param namS the name of local function to be used.
	 * @return the function value
	 */
	public static double brent1(double ax, double bx, double cx, double tol,
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

	/**
	 * Given a function func, and given distinct initial points ax and bx, this routine searches in 
	 * the downhill direction (defined by the function as evaluated at the initial points) and returns 
	 * new points ax_mnbrak, bx_mnbrak, cx_mnbrak that bracket a minimum of the function. Also 
	 * returned are the function values at the three points, fa_mnbrak, fb_mnbrak, and fc_mnbrak.
	 * @param ax ax
	 * @param bx bx
	 * @param namS the name of local function to be used.
	 */
	public static void mnbrak1(double ax, double bx, String namS)// , float *cx,
																	// float
																	// *fa,
																	// float
																	// *fb,
																	// float
																	// *fc,
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

	/**
	 * Straight-line fit to input data x[1..ndat] and y[1..ndat] with errors in both x and y, the respective 
	 * standard deviations being the input quantities sigx[1..ndat] and sigy[1..ndat]. 
	 * Output quantities are a and b such that y = a + bx minimizes chi2, whose value is returned as chi2_fitexy. 
	 * The chi2 probability is returned as q_fitexy, a small value indicating a poor fit (sometimes indicating underestimated errors). 
	 * Standard errors on a and b are returned as siga_fitexy and sigb_fitexy. These are not meaningful if either (i) the fit is poor, 
	 * or (ii) b is so large that the data are consistent with a vertical (infinite b) line. If siga and sigb are returned as BIG, then the data 
	 * are consistent with all values of b.
	 * @param x x
	 * @param y y
	 * @param ndat ndat
	 * @param sigx sigx
	 * @param sigy sigy
	 */
	public static void fitexy(double[] x, double[] y, int ndat, double[] sigx,
			double[] sigy)// ,
	// float *a, float *b, float *siga, float *sigb, float *chi2, float *q)
	// Straight-line fit to input data x[1..ndat] and y[1..ndat] with errors in
	// both x and y, the respective
	// standard deviations being the input quantities sigx[1..ndat] and
	// sigy[1..ndat].
	// Output quantities are a and b such that y = a + bx minimizes ?2, whose
	// value is returned
	// as chi2. The ?2 probability is returned as q, a small value indicating a
	// poor fit (sometimes
	// indicating underestimated errors). Standard errors on a and b are
	// returned as siga and sigb.
	// These are not meaningful if either (i) the fit is poor, or (ii) b is so
	// large that the data are
	// consistent with a vertical (infinite b) line. If siga and sigb are
	// returned as BIG, then the data
	// are consistent with all values of b.
	{
		// void avevar(float data[], unsigned long n, float *ave, float *var);
		// float brent(float ax, float bx, float cx,
		// float (*f)(float), float tol, float *xmin);
		// float chixy(float bang);
		// void fit(float x[], float y[], int ndata, float sig[], int mwt,
		// float *a, float *b, float *siga, float *sigb, float *chi2, float *q);
		// float gammq(float a, float x);
		// void mnbrak(float *ax, float *bx, float *cx, float *fa, float *fb,
		// float *fc, float (*func)(float));
		// float zbrent(float (*func)(float), float x1, float x2, float tol);

		int j = 0;
		double swap = 0.0;
		double amx = 0.0;
		double amn = 0.0;
		double varx = 0.0;
		double vary = 0.0;
		double[] ang = new double[7];
		double[] ch = new double[7];
		double scale = 0.0;
		double bmn = 0.0;
		double bmx = 0.0;
		double d1 = 0.0;
		double d2 = 0.0;
		double r2 = 0.0;
		//double dum1 = 0.0;
		//double dum2 = 0.0;
		//double dum3 = 0.0;
		//double dum4 = 0.0;
		//double dum5 = 0.0;

		// xx=vector(1,ndat);
		// yy=vector(1,ndat);
		// sx=vector(1,ndat);
		// sy=vector(1,ndat);
		// ww=vector(1,ndat);
		xx_fitexy = new double[ndat];
		yy_fitexy = new double[ndat];
		sx_fitexy = new double[ndat];
		sy_fitexy = new double[ndat];
		ww_fitexy = new double[ndat];

		Stats.avevar(x, ndat);// ,&dum1,&varx); //Find the x and y variances,
								// and scale
		// the data into the global variables for communication with the
		// function chixy.
		//dum1 = Stats.ave_avevar;
		varx = Stats.var_avevar;
		Stats.avevar(y, ndat);// ,&dum1,&vary);
		//dum1 = Stats.ave_avevar;
		vary = Stats.var_avevar;
		scale = Math.sqrt(varx / vary);
		nn_fitexy = ndat;
		for (j = 1; j <= ndat; j++) {
			xx_fitexy[j - 1] = x[j - 1];
			yy_fitexy[j - 1] = y[j - 1] * scale;
			sx_fitexy[j - 1] = sigx[j - 1];
			sy_fitexy[j - 1] = sigy[j - 1] * scale;
			ww_fitexy[j - 1] = Math.sqrt((sx_fitexy[j - 1])
					* (sx_fitexy[j - 1]) + (sy_fitexy[j - 1])
					* (sy_fitexy[j - 1]));
			// Use both x and y weights in first trial fit.
		}
		fit(xx_fitexy, yy_fitexy, nn_fitexy, ww_fitexy, 1);// ,&dum1,b,&dum2,&dum3,&dum4,&dum5);
															// Trial fit for b.
		// float *a, float *b, float *siga, float *sigb, float *chi2, float *q)
		//dum1 = a_fit;
		b_fitexy = b_fit;
		//dum2 = siga_fit;
		//dum3 = sigb_fit;
		//dum4 = chi2_fit;
		//dum5 = q_fit;
		offs_fitexy = ang[0] = 0.0; // ang[1]=0.0;
		// Construct several angles for reference points, and make b an angle.
		ang[1] = Math.atan(b_fitexy);// ang[2]=atan(*b);
		ang[3] = 0.0;// ang[4]=0.0;
		ang[4] = ang[1];// ang[5]=ang[2];
		ang[5] = POTN;// ang[6]=POTN;
		for (j = 4; j <= 6; j++)
			// ch[j]=chixy(ang[j]);
			ch[j - 1] = chixy(ang[j - 1]);
		// MinMaxFunc.mnbrak(&ang[1],&ang[2],&ang[3],&ch[1],&ch[2],&ch[3],chixy);
		mnbrak1(ang[0], ang[1], "chixy");// ,&ang[3],&ch[1],&ch[2],&ch[3],chixy);
		ang[2] = cx_mnbrak;
		ch[0] = fa_mnbrak;
		ch[1] = fb_mnbrak;
		ch[2] = fc_mnbrak;
		// Bracket the ?2 minimum and then locate it with brent.
		// chi2_fitexy=brent1(ang[1],ang[2],ang[3],"chixy",ACC,b);
		chi2_fitexy = brent1(ang[0], ang[1], ang[2], ACC, "chixy");
		b_fitexy = xmin_brent;
		chi2_fitexy = chixy(b_fitexy);
		a_fitexy = aa_fitexy;
		q_fitexy = SpecFunc.gammq(0.5 * (nn_fitexy - 2), chi2_fitexy * 0.5); // Compute
																				// ?2
																				// probability.
		for (r2 = 0.0, j = 1; j <= nn_fitexy; j++)
			r2 += ww_fitexy[j - 1];// r2 += ww[j]; Save the inverse sum of
									// weights at the minimum.
		r2 = 1.0 / r2;
		bmx = BIG; // Now, find standard errors for b as points where ??2 = 1.
		bmn = BIG;
		offs_fitexy = (chi2_fitexy) + 1.0;
		for (j = 1; j <= 6; j++) {// Go through saved values to bracketthe
									// desired roots. Note periodicity
									// in slope angles.
			if (ch[j - 1] > offs_fitexy) // if (ch[j] > offs)
			{
				d1 = Math.abs(ang[j - 1] - (b_fitexy));// d1=Math.abs(ang[j]-(*b));
				while (d1 >= PI)
					d1 -= PI;
				d2 = PI - d1;
				if (ang[j - 1] < b_fitexy)// if (ang[j] < *b)
				{
					swap = d1;
					d1 = d2;
					d2 = swap;
				}
				if (d1 < bmx)
					bmx = d1;
				if (d2 < bmn)
					bmn = d2;
			}
		}
		if (bmx < BIG) {// Call zbrent to find the roots.
						// bmx=zbrent(chixy,*b,*b+bmx,ACC)-(*b);
			bmx = zbrent1(b_fitexy, b_fitexy + bmx, ACC, "chixy") - (b_fitexy);
			amx = aa_fitexy - (a_fitexy);
			// bmn=zbrent(chixy,*b,*b-bmn,ACC)-(*b);
			bmn = zbrent1(b_fitexy, b_fitexy - bmn, ACC, "chixy") - (b_fitexy);
			amn = aa_fitexy - (a_fitexy);
			sigb_fitexy = Math.sqrt(0.5 * (bmx * bmx + bmn * bmn))
					/ (scale * (Math.cos(b_fitexy)) * (Math.cos(b_fitexy)));
			siga_fitexy = Math.sqrt(0.5 * (amx * amx + amn * amn) + r2) / scale; // Error
																					// in
																					// a
																					// has
																					// additional
																					// piece
																					// r2.
		} else
			(sigb_fitexy) = (siga_fitexy) = BIG;
		a_fitexy /= scale;// Unscale the answers.
		b_fitexy = Math.tan(b_fitexy) / scale;
		// free_vector(ww,1,ndat);
		// free_vector(sy,1,ndat);
		// free_vector(sx,1,ndat);
		// free_vector(yy,1,ndat);
		// free_vector(xx,1,ndat);
	}

	/**
	 * Local function.
	 * @param x x
	 * @param name name
	 * @return the result
	 */
	private static double localFunc(double x, String name) {
		double res = 0.0;
		if (name.compareTo("chixy") == 0) {
			res = chixy(x);
		}
		return res;
	}

	/**
	 * Captive function of fitexy, returns the value of (chi2 - offs) for the 
	 * slope b=tan(bang). Scaled data and offs are communicated via the global variables.
	 * @param bang bang
	 * @return the result
	 */
	public static double chixy(double bang)
	// Captive function of fitexy, returns the value of (chi2 - offs) for the
	// slope b=tan(bang).
	// Scaled data and offs are communicated via the global variables.
	{
		int j = 0;
		double ans;
		double avex = 0.0;
		double avey = 0.0;
		double sumw = 0.0;
		double b = 0.0;
		b = Math.tan(bang);
		for (j = 1; j <= nn_fitexy; j++) {
			// ww[j] = SQR(b*sx[j])+SQR(sy[j]);
			ww_fitexy[j - 1] = (b * sx_fitexy[j - 1]) * (b * sx_fitexy[j - 1])
					+ (sy_fitexy[j - 1]) * (sy_fitexy[j - 1]);
			// sumw += (ww[j] = (ww[j] < 1.0/BIG ? BIG : 1.0/ww[j]));
			sumw += (ww_fitexy[j - 1] = (ww_fitexy[j - 1] < 1.0 / BIG ? BIG
					: 1.0 / ww_fitexy[j - 1]));
			avex += ww_fitexy[j - 1] * xx_fitexy[j - 1];// avex += ww[j]*xx[j];
			avey += ww_fitexy[j - 1] * yy_fitexy[j - 1];// avey += ww[j]*yy[j];
		}
		avex /= sumw;
		avey /= sumw;
		aa_fitexy = avey - b * avex;
		for (ans = -offs_fitexy, j = 1; j <= nn_fitexy; j++)
			// ans += ww[j]*SQR(yy[j]-aa-b*xx[j]);
			ans += ww_fitexy[j - 1]
					* (yy_fitexy[j - 1] - aa_fitexy - b * xx_fitexy[j - 1])
					* (yy_fitexy[j - 1] - aa_fitexy - b * xx_fitexy[j - 1]);
		return ans;
	}

	/*
	 * General Linear Least Squares An immediate generalization of §15.2 is to
	 * fit a set of data points (xi, yi) to a model that is not just a linear
	 * combination of 1 and x (namely a + bx), but rather a linear combination
	 * of any M specified functions of x. For example, the functions could be 1,
	 * x, x2, . . . , xM-1, in which case their general linear combination, y(x)
	 * = a1 + a2x + a3x^2 + · · · + aMx^M-1
	 * 
	 * Or, the functions could be sines and cosines, in which case their general
	 * linear combination is a harmonic series. The general form of this kind of
	 * model is y(x) =sum k=1,M of akXk(x) (15.4.2) where X1(x), . . .,XM(x) are
	 * arbitrary fixed functions of x, called the basis functions.
	 * 
	 * We will now give a routine that implements the above formulas for the
	 * general linear least-squares problem, by the method of normal equations.
	 * Since we wish to compute not only the solution vector a but also the
	 * covariance matrix [C], it is most convenient to use Gauss-Jordan
	 * elimination (routine gaussj of §2.1) to perform the linear algebra. The
	 * operation count, in this application, is no larger than that for LU
	 * decomposition. orm the equations, and Gauss-Jordan is quite adequate. We
	 * need to warn you that the solution of a least-squares problem directly
	 * from the normal equations is rather susceptible to roundoff error. An
	 * alternative, and preferred, technique involves QR decomposition (§2.10,
	 * §11.3, and §11.6) of the design matrix A. This is essentially what we did
	 * at the end of §15.2 for fitting data to a straight line, but without
	 * invoking all the machinery of QR to derive the necessary formulas. Later
	 * in this section, we will discuss other difficulties in the least-squares
	 * problem, for which the cure is singular value decomposition (SVD), of
	 * which we give an implementation. It turns out that SVD also fixes the
	 * roundoff problem, so it is our recommended technique for all but “easy”
	 * least-squares problems. It is for these easy problems that the following
	 * routine, which solves the normal equations, is intended.
	 */

	/**
	 * Given a set of data points x[1..ndat], y[1..ndat] with individual standard deviations 
	 * sig[1..ndat], use chi2 minimization to fit for some or all of the coefficients a[1..ma] of 
	 * a function that depends linearly on a, y = sum i a_i × afunc_i(x). The input array ia[1..ma] 
	 * indicates by nonzero entries those components of a that should be fitted for, and by zero entries 
	 * those components that should be held fixed at their input values. The program returns values 
	 * for a[1..ma], chi2 = chisq_lfit, and the covariance matrix covar[1..ma][1..ma]. (Parameters 
	 * held fixed will return zero covariances). The user supplies a routine funcs(x,afunc,ma) that 
	 * returns the ma basis functions evaluated at x = x in the array afunc[1..ma]. Mainly used for 
	 * polynomial fit of data.
	 * @param x x
	 * @param y y
	 * @param sig sig
	 * @param ndat ndat
	 * @param a a
	 * @param ia ia
	 * @param ma ma
	 * @param covar covar
	 */
	public static void lfit(double[] x, double[] y, double[] sig, int ndat,
			double[] a, int[] ia, int ma, double[][] covar)// , float *chisq,
															// void
															// (*funcs)(float,
															// float [], int))
	// Given a set of data points x[1..ndat], y[1..ndat] with individual
	// standard deviations
	// sig[1..ndat], use chi2 minimization to fit for some or all of the
	// coefficients a[1..ma] of
	// a function that depends linearly on a, y = sum i a_i × afunc_i(x). The
	// input array ia[1..ma]
	// indicates by nonzero entries those components of a that should be fitted
	// for, and by zero entries
	// those components that should be held fixed at their input values. The
	// program returns values
	// for a[1..ma], chi2 = chisq, and the covariance matrix
	// covar[1..ma][1..ma]. (Parameters
	// held fixed will return zero covariances.)Th e user supplies a routine
	// funcs(x,afunc,ma) that
	// returns the ma basis functions evaluated at x = x in the array
	// afunc[1..ma].
	{
		// void covsrt(float **covar, int ma, int ia[], int mfit);
		// void gaussj(float **a, int n, float **b, int m);
		int i = 0;
		int j = 0;
		int k = 0;
		int l = 0;
		int m = 0;
		int mfit = 0;
		double ym = 0.0;
		double wt = 0.0;
		double sum = 0.0;
		double sig2i = 0.0;
		double[][] beta = new double[ma][1];
		double[] afunc = new double[ma];
		failB = false;
		// beta=matrix(1,ma,1,1);
		// afunc=vector(1,ma);
		for (j = 1; j <= ma; j++)
			if (ia[j - 1] != 0)
				mfit++;// if (ia[j]) mfit++;
		if (mfit == 0) {
			failB = true;
			failS = "lfit: no parameters to be fitted";
			return;
			// nrerror("lfit: no parameters to be fitted");
		}
		for (j = 1; j <= mfit; j++) {// Initialize the (symmetric)mat rix.
			for (k = 1; k <= mfit; k++)
				covar[j - 1][k - 1] = 0.0;// covar[j][k]=0.0;
			beta[j - 1][0] = 0.0;// beta[j][1]=0.0;
		}
		for (i = 1; i <= ndat; i++) {// Loop over data to accumulate
										// coefficients of the normal equations.
			afunc = func.aF(x[i - 1], ma);// (*funcs)(x[i],afunc,ma);
			ym = y[i - 1];// ym=y[i];
			if (mfit < ma) {// Subtract off dependences on known pieces of the
							// fitting function.
				for (j = 1; j <= ma; j++)
					// if (!ia[j]) ym -= a[j]*afunc[j];
					if (ia[j - 1] == 0)
						ym -= a[j - 1] * afunc[j - 1];
			}
			sig2i = 1.0 / ((sig[i - 1]) * (sig[i - 1]));// sig2i=1.0/SQR(sig[i]);
			for (j = 0, l = 1; l <= ma; l++) {
				if (ia[l - 1] != 0) // if (ia[l])
				{
					wt = afunc[l - 1] * sig2i;// wt=afunc[l]*sig2i;
					for (j++, k = 0, m = 1; m <= l; m++)
						// if (ia[m]) covar[j][++k] += wt*afunc[m];
						if (ia[m - 1] != 0)
							covar[j - 1][++k - 1] += wt * afunc[m - 1];
					beta[j - 1][0] += ym * wt;// beta[j][1] += ym*wt;
				}
			}
		}
		for (j = 2; j <= mfit; j++)
			// Fill in above the diagonal from symmetry.
			for (k = 1; k < j; k++)
				covar[k - 1][j - 1] = covar[j - 1][k - 1];// covar[k][j]=covar[j][k];
		LinAEq.gaussj(covar, mfit, beta, 1); // Matrix solution.
		for (j = 0, l = 1; l <= ma; l++)
			// if (ia[l]) a[l]=beta[++j][1]; Partition solution to appropriate
			// coefficients a.
			if (ia[l - 1] != 0)
				a[l - 1] = beta[++j - 1][0];
		chisq_lfit = 0.0;
		for (i = 1; i <= ndat; i++) {// Evaluate ?2 of the fit.
			afunc = func.aF(x[i - 1], ma);// (*funcs)(x[i],afunc,ma);
			for (sum = 0.0, j = 1; j <= ma; j++)
				sum += a[j - 1] * afunc[j - 1];// sum += a[j]*afunc[j];
			chisq_lfit += ((y[i - 1] - sum) / sig[i - 1])
					* ((y[i - 1] - sum) / sig[i - 1]);// SQR((y[i]-sum)/sig[i]);
		}
		covsrt(covar, ma, ia, mfit);// Sort covariance matrix to true order of
									// fitting coefficients.
		// free_vector(afunc,1,ma);
		// free_matrix(beta,1,ma,1,1);
	}

	// #define SWAP(a,b) {swap=(a);(a)=(b);(b)=swap;}
	/**
	 * Used internally. Expand in storage the covariance matrix covar, so as to take into account 
	 * parameters that are being held fixed. (For the latter, return zero covariances).
	 * @param covar covar
	 * @param ma ma
	 * @param ia ia
	 * @param mfit mfit
	 */
	public static void covsrt(double[][] covar, int ma, int ia[], int mfit)
	// Expand in storage the covariance matrix covar, so as to take into account
	// parameters that are
	// being held fixed. (For the latter, return zero covariances.)
	{
		int i = 0;
		int j = 0;
		int k = 0;
		double swap = 0.0;
		for (i = mfit + 1; i <= ma; i++)
			for (j = 1; j <= i; j++)
				covar[i - 1][j - 1] = covar[j - 1][i - 1] = 0.0;// covar[i][j]=covar[j][i]=0.0;
		k = mfit;
		for (j = ma; j >= 1; j--) {
			if (ia[j - 1] != 0) // if (ia[j])
			{
				for (i = 1; i <= ma; i++) // SWAP(covar[i][k],covar[i][j])
				{
					swap = covar[i - 1][k - 1];
					covar[i - 1][k - 1] = covar[i - 1][j - 1];
					covar[i - 1][j - 1] = swap;
				}
				for (i = 1; i <= ma; i++) // SWAP(covar[k][i],covar[j][i])
				{
					swap = covar[k - 1][i - 1];
					covar[k - 1][i - 1] = covar[j - 1][i - 1];
					covar[j - 1][i - 1] = swap;
				}
				k--;
			}
		}
	}

	/*
	 * In some applications, the normal equations are perfectly adequate for
	 * linear least-squares problems. However, in many cases the normal
	 * equations are very close to singular. A zero pivot element may be
	 * encountered during the solution of the linear equations (e.g., in
	 * gaussj), in which case you get no solution at all. Or a very small pivot
	 * may occur, in which case you typically get fitted parameters a k with
	 * very large magnitudes that are delicately (and unstably) balanced to
	 * cancel out almost precisely when the fitted function is evaluated. Why
	 * does this commonly occur? The reason is that, more often than
	 * experimenters would like to admit, data do not clearly distinguish
	 * between two or more of the basis functions provided. If two such
	 * functions, or two different combinations of functions, happen to fit the
	 * data about equally well — or equally badly — then the matrix [?], unable
	 * to distinguish between them, neatly folds up its tent and becomes
	 * singular. There is a certain mathematical irony in the fact that
	 * least-squares problems are both overdetermined (number of data points
	 * greater than number of parameters) and underdetermined (ambiguous
	 * combinations of parameters exist); but that is how it frequently is. The
	 * ambiguities can be extremely hard to notice a priori in complicated
	 * problems
	 * 
	 * Generally speaking, we recommend that you always use SVD techniques
	 * instead of using the normal equations. SVD’s only significant
	 * disadvantage is that it requires an extra array of size N × M to store
	 * the whole design matrix. This storage is overwritten by the matrix U.
	 * Storage is also required for the M × M matrix V, but this is instead of
	 * the same-sized matrix for the coefficients of the normal equations. SVD
	 * can be significantly slower than solving the normal equations; however,
	 * its great advantage, that it (theoretically) cannot fail, more than makes
	 * up for the speed disadvantage. In the routine that follows, the matrices
	 * u,v and the vector w are input as working space. The logical dimensions
	 * of the problem are ndata data points by ma basis functions (and fitted
	 * parameters). If you care only about the values a of the fitted
	 * parameters, then u,v,w contain no useful information on output. If you
	 * want probable errors for the fitted parameters, read on.
	 */
	/**
	 * Given a set of data points x[1..ndata],y[1..ndata] with individual standard deviations sig[1..ndata], use chi2 minimization to determine the coefficients 
	 * a[1..ma] of the fitting function y = sum i a_i × afunc_i(x). Here we solve the fitting equations using singular 
	 * value decomposition of the ndata by ma matrix. Arrays u[1..ndata][1..ma], v[1..ma][1..ma], and w[1..ma] provide workspace on input; on output they 
	 * define the singular value decomposition, and can be used to obtain the covariance matrix. The program 
	 * returns values for the ma fit parameters a, and chi2, chisq_svdfit. The user supplies a routine 
	 * funcs(x,afunc,ma) that returns the ma basis functions evaluated at x = x in the array afunc[1..ma].
	 * @param x x
	 * @param y y
	 * @param sig sig
	 * @param ndata ndata
	 * @param a a
	 * @param ma ma
	 * @param u u
	 * @param v v
	 * @param w w
	 */
	public static void svdfit(double[] x, double[] y, double[] sig, int ndata,
			double[] a, int ma, double[][] u, double[][] v, double[] w)// ,
																		// float
																		// *chisq,void
																		// (*funcs)(float,
																		// float
																		// [],
																		// int))
	// Given a set of data points x[1..ndata],y[1..ndata] with individual
	// standard deviations
	// sig[1..ndata], use chi2 minimization to determine the coefficients
	// a[1..ma] of the fitting
	// function y = sum i a_i × afunc_i(x). Here we solve the fitting equations
	// using singular
	// value decomposition of the ndata by ma matrix, as in §2.6. Arrays
	// u[1..ndata][1..ma],
	// v[1..ma][1..ma], and w[1..ma] provide workspace on input; on output they
	// define the
	// singular value decomposition, and can be used to obtain the covariance
	// matrix. The program
	// returns values for the ma fit parameters a, and chi2, chisq. The user
	// supplies a routine
	// funcs(x,afunc,ma) that returns the ma basis functions evaluated at x = x
	// in the array
	// afunc[1..ma].
	{
		// void svbksb(float **u, float w[], float **v, int m, int n, float
		// b[],float x[]);
		// void svdcmp(float **a, int m, int n, float w[], float **v);
		int j = 0;
		int i = 0;
		double wmax = 0.0;
		double tmp = 0.0;
		double thresh = 0.0;
		double sum = 0.0;
		double[] b = new double[ndata];
		double[] afunc = new double[ma];
		// b=vector(1,ndata);
		// afunc=vector(1,ma);
		for (i = 1; i <= ndata; i++) {// Accumulate coefficients of the fitting
										// matrix.
			afunc = func.aF(x[i - 1], ma);// (*funcs)(x[i],afunc,ma);
			tmp = 1.0 / sig[i - 1];// tmp=1.0/sig[i];
			for (j = 1; j <= ma; j++)
				u[i - 1][j - 1] = afunc[j - 1] * tmp;// u[i][j]=afunc[j]*tmp;
			b[i - 1] = y[i - 1] * tmp;// b[i]=y[i]*tmp;
		}
		LinAEq.svdcmp(u, ndata, ma, w, v);// Singular value decomposition.
		wmax = 0.0;// Edit the singular values, given TOL from the
		// #define statement, between here ...
		for (j = 1; j <= ma; j++)
			if (w[j - 1] > wmax)
				wmax = w[j - 1];// if (w[j] > wmax) wmax=w[j];
		thresh = TOL * wmax;
		for (j = 1; j <= ma; j++)
			if (w[j - 1] < thresh)
				w[j - 1] = 0.0;// if (w[j] < thresh) w[j]=0.0; ...and here.
		LinAEq.svbksb(u, w, v, ndata, ma, b, a);
		chisq_svdfit = 0.0; // Evaluate chi-square.
		for (i = 1; i <= ndata; i++) {
			afunc = func.aF(x[i - 1], ma);// (*funcs)(x[i],afunc,ma);
			for (sum = 0.0, j = 1; j <= ma; j++)
				sum += a[j - 1] * afunc[j - 1];// sum += a[j]*afunc[j];
			// chisq_svdfit += (tmp=(y[i]-sum)/sig[i],tmp*tmp);
			chisq_svdfit += ((y[i - 1] - sum) / sig[i - 1])
					* ((y[i - 1] - sum) / sig[i - 1]);
		}
		// free_vector(afunc,1,ma);
		// free_vector(b,1,ndata);
	}

	/*
	 * Feeding the matrix v and vector w output by the above program into the
	 * following short routine, you easily obtain variances and covariances of
	 * the fitted parameters a. The square roots of the variances are the
	 * standard deviations of the fitted parameters. The routine
	 * straightforwardly implements equation (15.4.20) above, with the
	 * convention that singular values equal to zero are recognized as having
	 * been edited out of the fit.
	 */

	/**
	 * To evaluate the covariance matrix cvm[1..ma][1..ma] of the fit for ma parameters obtained 
	 * by svdfit, call this routine with matrices v[1..ma][1..ma], w[1..ma] as returned from svdfit.
	 * @param v v
	 * @param ma ma
	 * @param w w
	 * @param cvm cvm
	 */
	public static void svdvar(double[][] v, int ma, double[] w, double[][] cvm)
	// To evaluate the covariance matrix cvm[1..ma][1..ma] of the fit for ma
	// parameters obtained
	// by svdfit, call this routine with matrices v[1..ma][1..ma], w[1..ma] as
	// returned from
	// svdfit.
	{
		int k = 0;
		int j = 0;
		int i = 0;
		double sum = 0.0;
		double[] wti = new double[ma];
		// wti=vector(1,ma);
		for (i = 1; i <= ma; i++) {
			wti[i - 1] = 0.0;// wti[i]=0.0;
			// if (w[i]) wti[i]=1.0/(w[i]*w[i]);
			if (w[i - 1] != 0)
				wti[i - 1] = 1.0 / (w[i - 1] * w[i - 1]);
		}
		for (i = 1; i <= ma; i++) {// Sum contributions to covariance matrix
									// (15.4.20).
			for (j = 1; j <= i; j++) {
				for (sum = 0.0, k = 1; k <= ma; k++)
					// sum += v[i][k]*v[j][k]*wti[k];
					sum += v[i - 1][k - 1] * v[j - 1][k - 1] * wti[k - 1];
				// cvm[j][i]=cvm[i][j]=sum;
				cvm[j - 1][i - 1] = cvm[i - 1][j - 1] = sum;
			}
		}
		// free_vector(wti,1,ma);
	}

	/*
	 * Examples Be aware that some apparently nonlinear problems can be
	 * expressed so that they are linear. For example, an exponential model with
	 * two parameters a and b, y(x) = a exp(-bx) (15.4.21) can be rewritten as
	 * log[y(x)] = c - bx (15.4.22) which is linear in its parameters c and b.
	 * (Of course you must be aware that such transformations do not exactly
	 * take Gaussian errors into Gaussian errors.) Also watch out for
	 * “non-parameters,” as in y(x) = a exp(-bx + d) (15.4.23) Here the
	 * parameters a and d are, in fact, indistinguishable. This is a good
	 * example of where the normal equations will be exactly singular, and where
	 * SVD will find a zero singular value. SVD will then make a “least-squares”
	 * choice for setting a balance between a and d (or, rather, their
	 * equivalents in the linear model derived by taking the logarithms).
	 * However — and this is true whenever SVD gives back a zero singular value
	 * — you are better advised to figure out analytically where the degeneracy
	 * is among your basis functions, and then make appropriate deletions in the
	 * basis set. Here are two examples for user-supplied routines funcs. The
	 * first one is trivial and fits a general polynomial to a set of data:
	 * 
	 * void fpoly(float x, float p[], int np) Fitting routine for a polynomial
	 * of degree np-1, with coefficients in the array p[1..np]. { int j;
	 * p[1]=1.0; for (j=2;j<=np;j++) p[j]=p[j-1]*x; }
	 * 
	 * The second example is slightly less trivial. It is used to fit Legendre
	 * polynomials up to some order nl-1 through a data set. void fleg(float x,
	 * float pl[], int nl) Fitting routine for an expansion with nl Legendre
	 * polynomials pl, evaluated using the recurrence relation as in §5.5. { int
	 * j; float twox,f2,f1,d; pl[1]=1.0; pl[2]=x; if (nl > 2) { twox=2.0*x;
	 * f2=x; d=1.0; for (j=3;j<=nl;j++) { f1=d++; f2 += twox;
	 * pl[j]=(f2*pl[j-1]-f1*pl[j-2])/d; } } }
	 */

	/*
	 * Nonlinear Models We now consider fitting when the model depends
	 * nonlinearly on the set of M unknown parameters ak, k = 1, 2, . . .,M. We
	 * use the same approach as in previous sections, namely to define a ?2
	 * merit function and determine best-fit parameters by its minimization.
	 * With nonlinear dependences, however, the minimization must proceed
	 * iteratively. Given trial values for the parameters, we develop a
	 * procedure that improves the trial solution. The procedure is then
	 * repeated until ?2 stops (or effectively stops) decreasing.
	 * 
	 * The following pair of functions encodes Marquardt’s method for nonlinear
	 * parameter estimation. Much of the organization matches that used in lfit
	 * of §15.4. In particular the array ia[1..ma] must be input with components
	 * one or zero corresponding to whether the respective parameter values
	 * a[1..ma] are to be fitted for or held fixed at their input values,
	 * respectively. The routine mrqmin performs one iteration of Marquardt’s
	 * method. It is first called (once) with alamda < 0, which signals the
	 * routine to initialize. alamda is set on the first and all subsequent
	 * calls to the suggested value of ? for the next iteration; a and chisq are
	 * always given back as the best parameters found so far and their ?2. When
	 * convergence is deemed satisfactory, set alamda to zero before a final
	 * call. The matrices alpha and covar (which were used as workspace in all
	 * previous calls) will then be set to the curvature and covariance matrices
	 * for the converged parameter values. The arguments alpha, a, and chisq
	 * must not be modified between calls, nor should alamda be, except to set
	 * it to zero for the final call. When an uphill step is taken, chisq and a
	 * are given back with their input (best) values, but alamda is set to an
	 * increased value. The routine mrqmin calls the routine mrqcof for the
	 * computation of the matrix [?] (equation 15.5.11) and vector ß (equations
	 * 15.5.6 and 15.5.8). In turn mrqcof calls the user-supplied routine
	 * funcs(x,a,y,dyda), which for input values x ? xi and a ? a calculates the
	 * model function y ? y(xi; a) and the vector of derivatives dyda ? ?y/?ak.
	 */

	/**
	 * Levenberg-Marquardt method, attempting to reduce the value chi2 of a fit between a set of data points x[1..ndata], y[1..ndata] with individual standard deviations 
	 * sig[1..ndata], and a nonlinear function dependent on ma coefficients a[1..ma]. The input array ia[1..ma] indicates by nonzero entries those components of a that should be fitted 
	 * for, and by zero entries those components that should be held fixed at their input values. The program returns 
	 * current best-fit values for the parameters a[1..ma], and chi2 = chisq_mrqmin. The arrays covar[1..ma][1..ma], alpha[1..ma][1..ma] are used as working space during  
	 * most iterations. Supply a routine funcs(x,a,yfit,dyda,ma) that evaluates the fitting function yfit, and its derivatives dyda[1..ma] with respect to the fitting 
	 * parameters a at x. On the first call provide an initial guess for the parameters a, and set alamda less than 0 for initialization 
	 * (which then sets alamda=.001). If a step succeeds chisq becomes smaller and alamda decreases by a factor of 10. If a step fails alamda grows by a factor of 10. You 
	 * must call this routine repeatedly until convergence is achieved. Then, make one final call with alamda=0, so that covar[1..ma][1..ma] returns the covariance matrix, and alpha the 
	 * curvature matrix (parameters held fixed will return zero covariances).
	 * @param x x
	 * @param y y
	 * @param sig sig
	 * @param ndata ndata
	 * @param a a
	 * @param ia ia
	 * @param ma ma
	 * @param covar covar
	 * @param alpha alpha
	 * @param alamda alambda
	 */
	public static void mrqmin(double[] x, double[] y, double[] sig, int ndata,
			double[] a, int[] ia, int ma, double[][] covar, double[][] alpha,
			double alamda)// , float *chisq,
	// void (*funcs)(float, float [], float *, float [], int), float *alamda)
	// Levenberg-Marquardt method, attempting to reduce the value ?2 of a fit
	// between a set of data
	// points x[1..ndata], y[1..ndata] with individual standard deviations
	// sig[1..ndata],
	// and a nonlinear function dependent on ma coefficients a[1..ma]. The input
	// array ia[1..ma]
	// indicates by nonzero entries those components of a that should be fitted
	// for, and by zero
	// entries those components that should be held fixed at their input values.
	// The program returns
	// current best-fit values for the parameters a[1..ma], and ?2 = chisq. The
	// arrays
	// covar[1..ma][1..ma], alpha[1..ma][1..ma] are used as working space during
	// most
	// iterations. Supply a routine funcs(x,a,yfit,dyda,ma) that evaluates the
	// fitting function
	// yfit, and its derivatives dyda[1..ma] with respect to the fitting
	// parameters a at x. On
	// the first call provide an initial guess for the parameters a, and set
	// alamda<0 for initialization
	// (which then sets alamda=.001). If a step succeeds chisq becomes smaller
	// and alamda decreases
	// by a factor of 10. If a step fails alamda grows by a factor of 10. You
	// must call this
	// routine repeatedly until convergence is achieved. Then, make one final
	// call with alamda=0, so
	// that covar[1..ma][1..ma] returns the covariance matrix, and alpha the
	// curvature matrix.
	// (Parameters held fixed will return zero covariances.)
	{
		// void covsrt(float **covar, int ma, int ia[], int mfit);
		// void gaussj(float **a, int n, float **b, int m);
		// void mrqcof(float x[], float y[], float sig[], int ndata, float a[],
		// int ia[], int ma, float **alpha, float beta[], float *chisq,
		// void (*funcs)(float, float [], float *, float [], int));
		int j = 0;
		int k = 0;
		int l = 0;
		alamda_mrqmin = alamda;
		// static int mfit;
		// static float ochisq,*atry,*beta,*da,**oneda;
		if (alamda_mrqmin < 0.0) {// Initialization.
			atry_mrqmin = new double[ma];// vector(1,ma);
			beta_mrqmin = new double[ma];// vector(1,ma);
			da_mrqmin = new double[ma];// vector(1,ma);
			for (mfit_mrqmin = 0, j = 1; j <= ma; j++)
				if (ia[j - 1] != 0)
					mfit_mrqmin++;// if (ia[j]) mfit++;
			oneda_mrqmin = new double[mfit_mrqmin][1];// matrix(1,mfit,1,1);
			alamda_mrqmin = 0.001;
			mrqcof(x, y, sig, ndata, a, ia, ma, alpha, beta_mrqmin);// ,chisq_mrqmin,funcs);
			ochisq_mrqmin = (chisq_mrqmin);
			for (j = 1; j <= ma; j++)
				atry_mrqmin[j - 1] = a[j - 1];// atry[j]=a[j];
		}
		for (j = 1; j <= mfit_mrqmin; j++) {// Alter linearized fitting matrix,
											// by augmenting diagonal elements.
			for (k = 1; k <= mfit_mrqmin; k++)
				covar[j - 1][k - 1] = alpha[j - 1][k - 1];// covar[j][k]=alpha[j][k];
			covar[j - 1][j - 1] = alpha[j - 1][j - 1] * (1.0 + (alamda_mrqmin));// covar[j][j]=alpha[j][j]*(1.0+(*alamda));
			oneda_mrqmin[j - 1][0] = beta_mrqmin[j - 1];// oneda[j][1]=beta[j];
		}
		LinAEq.gaussj(covar, mfit_mrqmin, oneda_mrqmin, 1);// Matrix solution.
		for (j = 1; j <= mfit_mrqmin; j++)
			da_mrqmin[j - 1] = oneda_mrqmin[j - 1][0];// da[j]=oneda[j][1];
		if (alamda_mrqmin == 0.0) {// Once converged, evaluate covariance
									// matrix.
			covsrt(covar, ma, ia, mfit_mrqmin);
			covsrt(alpha, ma, ia, mfit_mrqmin); // Spread out alpha to its full
												// size too.
			// free_matrix(oneda,1,mfit,1,1);
			// free_vector(da,1,ma);
			// free_vector(beta,1,ma);
			// free_vector(atry,1,ma);
			return;
		}
		for (j = 0, l = 1; l <= ma; l++)
			// Did the trial succeed?
			// if (ia[l]) atry[l]=a[l]+da[++j];
			if (ia[l - 1] != 0)
				atry_mrqmin[l - 1] = a[l - 1] + da_mrqmin[++j - 1];
		mrqcof(x, y, sig, ndata, atry_mrqmin, ia, ma, covar, da_mrqmin);// ,chisq,funcs);
		if (chisq_mrqmin < ochisq_mrqmin) {// Success, accept the new solution.
			convB = true;
			alamda_mrqmin *= 0.1;
			ochisq_mrqmin = (chisq_mrqmin);
			for (j = 1; j <= mfit_mrqmin; j++) {
				for (k = 1; k <= mfit_mrqmin; k++)
					alpha[j - 1][k - 1] = covar[j - 1][k - 1];// alpha[j][k]=covar[j][k];
				beta_mrqmin[j - 1] = da_mrqmin[j - 1];// beta[j]=da[j];
			}
			for (l = 1; l <= ma; l++)
				a[l - 1] = atry_mrqmin[l - 1];// a[l]=atry[l];
		} else {// Failure, increase alamda and return.
			alamda_mrqmin *= 10.0;
			chisq_mrqmin = ochisq_mrqmin;
			convB = false;
		}
	}

	/**
	 * Used internally by mrqmin.
	 * @param x x
	 * @param y y
	 * @param sig sig
	 * @param ndata ndata
	 * @param a a
	 * @param ia ia
	 * @param ma ma
	 * @param alpha alpha
	 * @param beta beta
	 */
	public static void mrqcof(double[] x, double[] y, double[] sig, int ndata,
			double[] a, int[] ia, int ma, double[][] alpha, double[] beta)// ,
																			// float
																			// *chisq,
	// void (*funcs)(float, float [], float *, float [], int))
	// Used by mrqmin to evaluate the linearized fitting matrix alpha, and
	// vector beta as in (15.5.8),
	// and calculate ?2.
	{
		int i = 0;
		int j = 0;
		int k = 0;
		int l = 0;
		int m = 0;
		int mfit = 0;
		double ymod = 0.0;
		double wt = 0.0;
		double sig2i = 0.0;
		double dy = 0.0;
		double[] dyda = new double[ma];
		// dyda=vector(1,ma);
		for (j = 1; j <= ma; j++)
			if (ia[j - 1] != 0)
				mfit++;// if (ia[j]) mfit++;
		for (j = 1; j <= mfit; j++) {// Initialize (symmetric) alpha, beta.
			for (k = 1; k <= j; k++)
				alpha[j - 1][k - 1] = 0.0;// alpha[j][k]=0.0;
			beta[j - 1] = 0.0;// beta[j]=0.0;
		}
		chisq_mrqmin = 0.0;
		for (i = 1; i <= ndata; i++) {// Summation loop over all data.
										// (*funcs)(x[i],a,&ymod,dyda,ma);
										// System.out.println("before"
										// +dyda[2]);
			ymod = func.fdf(x[i - 1], a, dyda, ma);
			// System.out.println("after" +dyda[2]);
			sig2i = 1.0 / (sig[i - 1] * sig[i - 1]);// sig2i=1.0/(sig[i]*sig[i]);
			dy = y[i - 1] - ymod;// dy=y[i]-ymod;
			for (j = 0, l = 1; l <= ma; l++) {
				if (ia[l - 1] != 0) // if (ia[l])
				{
					wt = dyda[l - 1] * sig2i;// wt=dyda[l]*sig2i;
					for (j++, k = 0, m = 1; m <= l; m++)
						// if (ia[m]) alpha[j][++k] += wt*dyda[m];
						if (ia[m - 1] != 0)
							alpha[j - 1][++k - 1] += wt * dyda[m - 1];
					beta[j - 1] += dy * wt;// beta[j] += dy*wt;
				}
			}
			chisq_mrqmin += dy * dy * sig2i;// And find ?2.
		}
		for (j = 2; j <= mfit; j++)
			// Fill in the symmetric side.
			for (k = 1; k < j; k++)
				// alpha[k][j]=alpha[j][k];
				alpha[k - 1][j - 1] = alpha[j - 1][k - 1];
		// free_vector(dyda,1,ma);
	}

	/*
	 * Example The following function fgauss is an example of a user-supplied
	 * function funcs. Used with the above routine mrqmin (in turn using mrqcof,
	 * covsrt, and gaussj), it fits for the model y(x) =sum k=1 to K Bk exp [-(x
	 * - Ek)/Gk ]^2 (15.5.16) which is a sum of K Gaussians, each having a
	 * variable position, amplitude, and width. We store the parameters in the
	 * order B1,E1,G1,B2,E2,G2, . . . , BK, EK,GK.
	 * 
	 * void fgauss(float x, float a[], float *y, float dyda[], int na) y(x; a)
	 * is the sum of na/3 Gaussians (15.5.16). The amplitude, center, and width
	 * of the Gaussians are stored in consecutive locations of a: a[i] = Bk,
	 * a[i+1] = Ek, a[i+2] = Gk, k = 1, ..., na/3. The dimensions of the arrays
	 * are a[1..na], dyda[1..na]. { int i; float fac,ex,arg;y=0.0; for
	 * (i=1;i<=na-1;i+=3) { arg=(x-a[i+1])/a[i+2]; ex=exp(-arg*arg);
	 * fac=a[i]*ex*2.0*arg;y += a[i]*ex; dyda[i]=ex; dyda[i+1]=fac/a[i+2];
	 * dyda[i+2]=fac*arg/a[i+2]; } }
	 */

	/*
	 * Robust Estimation
	 * 
	 * The concept of robustness has been mentioned in passing several times
	 * already. In §14.1 we noted that the median was a more robust estimator of
	 * central value than the mean; in §14.6 it was mentioned that rank
	 * correlation is more robust than linear correlation. The concept of
	 * outlier points as exceptions to a Gaussian model for experimental error
	 * was discussed in §15.1. The term “robust” was coined in statistics by
	 * G.E.P. Box in 1953. Various definitions of greater or lesser mathematical
	 * rigor are possible for the term, but in general, referring to a
	 * statistical estimator, it means “insensitive to small departures from the
	 * idealized assumptions for which the estimator is optimized.” [1,2] The
	 * word “small” can have two different interpretations, both important:
	 * either fractionally small departures for all data points, or else
	 * fractionally large departures for a small number of data points. It is
	 * the latter interpretation, leading to the notion of outlier points, that
	 * is generally the most stressful for statistical procedures. Statisticians
	 * have developed various sorts of robust statistical estimators. Many, if
	 * not most, can be grouped in one of three categories. M-estimates follow
	 * from maximum-likelihood arguments very much as equations (15.1.5) and
	 * (15.1.7) followed from equation (15.1.3). M-estimates are usually the
	 * most relevant class for model-fitting, that is, estimation of parameters.
	 * We therefore consider these estimates in some detail below. L-estimates
	 * are “linear combinations of order statistics.” These are most applicable
	 * to estimations of central value and central tendency, though they can
	 * occasionally be applied to some problems in estimation of parameters. Two
	 * “typical” L-estimates will give you the general idea. They are (i) the
	 * median, and (ii) Tukey’s trimean, defined as the weighted average of the
	 * first, second, and third quartile points in a distribution, with weights
	 * 1/4, 1/2, and 1/4, respectively. R-estimates are estimates based on rank
	 * tests. For example, the equality or inequality of two distributions can
	 * be estimated by the Wilcoxon test of computing the mean rank of one
	 * distribution in a combined sample of both distributions.
	 * ............................ Fitting a Line by Minimizing Absolute
	 * Deviation
	 */

	/**
	 * Fits y = a + bx by the criterion of least absolute deviations. The arrays x[1..ndata] and 
	 * y[1..ndata] are the input experimental points. The fitted parameters a_medfit and b_medfit are output, 
	 * along with abdev_medfit, which is the mean absolute deviation (in y) of the experimental points from 
	 * the fitted line. This routine uses the routine rofunc, with communication via global variables.
	 * @param x x
	 * @param y y
	 * @param ndata ndata
	 */
	public static void medfit(double[] x, double[] y, int ndata)// , float *a,
																// float *b,
																// float *abdev)
	// Fits y = a + bx by the criterion of least absolute deviations. The arrays
	// x[1..ndata] and
	// y[1..ndata] are the input experimental points. The fitted parameters a
	// and b are output,
	// along with abdev, which is the mean absolute deviation (in y) of the
	// experimental points from
	// the fitted line. This routine uses the routine rofunc, with communication
	// via global variables.
	{
		// float rofunc(float b);
		int j = 0;
		double bb = 0.0;
		double b1 = 0.0;
		double b2 = 0.0;
		double del = 0.0;
		double f = 0.0;
		double f1 = 0.0;
		double f2 = 0.0;
		double sigb = 0.0;
		double temp = 0.0;
		double sx = 0.0;
		double sy = 0.0;
		double sxy = 0.0;
		double sxx = 0.0;
		double chisq = 0.0;
		ndatat_medfit = ndata;
		xt_medfit = new double[ndata];
		yt_medfit = new double[ndata];
		for (j = 1; j <= ndata; j++) {
			xt_medfit[j - 1] = x[j - 1];// System.out.println("--- "+xt_medfit[j-1]);
			yt_medfit[j - 1] = y[j - 1];
		}
		// xt_medfit=x;
		// yt_medfit=y;
		for (j = 1; j <= ndata; j++) {// As a first guess for a and b, we will
										// find the least-squares fitting line.
			sx += x[j - 1];// x[j];
			sy += y[j - 1];// y[j];
			sxy += x[j - 1] * y[j - 1];// x[j]*y[j];
			sxx += x[j - 1] * x[j - 1];// x[j]*x[j];
		}
		del = ndata * sxx - sx * sx;
		aa_medfit = (sxx * sy - sx * sxy) / del; // Least-squares solutions.
		bb = (ndata * sxy - sx * sy) / del;
		for (j = 1; j <= ndata; j++) {
			// chisq += (temp=y[j]-(aa+bb*x[j]),temp*temp);
			temp = y[j - 1] - (aa_medfit + bb * x[j - 1]);
			chisq += temp * temp;// System.out.println("--- "+chisq);
		}
		sigb = Math.sqrt(chisq / del);// The standard deviation will give some
										// idea of
		// how big an iteration step to take.
		b1 = bb;
		f1 = rofunc(b1);
		if (sigb > 0.0) {
			b2 = bb + SIGN(3.0 * sigb, f1);
			// Guess bracket as 3-sigma away, in the downhill direction known
			// from f1.
			f2 = rofunc(b2);
			if (b2 == b1) {
				a_medfit = aa_medfit;
				b_medfit = bb;
				abdev_medfit = abdevt_medfit / ndata;
				return;
			}
			while (f1 * f2 > 0.0) {// Bracketing.
				bb = b2 + 1.6 * (b2 - b1);
				b1 = b2;
				f1 = f2;
				b2 = bb;
				f2 = rofunc(b2);
			}
			sigb = 0.01 * sigb;// Refine until error a negligible number of
								// standard deviations.
			while (Math.abs(b2 - b1) > sigb) {
				bb = b1 + 0.5 * (b2 - b1); // Bisection.
				if (bb == b1 || bb == b2)
					break;
				f = rofunc(bb);
				if (f * f1 >= 0.0) {
					f1 = f;
					b1 = bb;
				} else {
					f2 = f;
					b2 = bb;
				}
			}
		}
		a_medfit = aa_medfit;
		b_medfit = bb;
		abdev_medfit = abdevt_medfit / ndata;
	}

	// extern int ndatat; Defined in medfit.
	// extern float *xt,*yt,aa,abdevt;

	/**
	 * Used internally by medfit.
	 * @param b b
	 * @return the result
	 */
	public static double rofunc(double b)
	// Evaluates the right-hand side of equation (15.7.16) for a given value of
	// b. Communication with
	// the routine medfit is through global variables.
	{
		// float select(unsigned long k, unsigned long n, float arr[]);Sorting.
		int j = 0;
		double[] arr = new double[ndatat_medfit];
		double d = 0.0;
		double sum = 0.0;
		// arr=vector(1,ndatat);
		for (j = 1; j <= ndatat_medfit; j++)
			arr[j - 1] = yt_medfit[j - 1] - b * xt_medfit[j - 1];// arr[j]=yt[j]-b*xt[j];
		// if (ndatat_medfit & 1)
		if ((ndatat_medfit & 1) != 0) {
			// aa_medfit=select((ndatat+1)>>1,ndatat,arr);
			aa_medfit = Sorting.select((ndatat_medfit + 1) >> 1, ndatat_medfit,
					arr);
		} else {
			j = ndatat_medfit >> 1;
			aa_medfit = 0.5 * (Sorting.select(j, ndatat_medfit, arr) + Sorting
					.select(j + 1, ndatat_medfit, arr));
		}
		abdevt_medfit = 0.0;
		for (j = 1; j <= ndatat_medfit; j++) {
			// d=yt_medfit[j]-(b*xt_medfit[j]+aa_medfit);
			d = yt_medfit[j - 1] - (b * xt_medfit[j - 1] + aa_medfit);
			abdevt_medfit += Math.abs(d);
			// if (yt_medfit[j] != 0.0) d /= Math.abs(yt_medfit[j]);
			if (yt_medfit[j - 1] != 0.0)
				d /= Math.abs(yt_medfit[j - 1]);
			// if (Math.abs(d) > EPS) sum += (d >= 0.0 ? xt_medfit[j] :
			// -xt_medfit[j]);
			if (Math.abs(d) > EPS2)
				sum += (d >= 0.0 ? xt_medfit[j - 1] : -xt_medfit[j - 1]);
		}
		// free_vector(arr,1,ndatat);
		return sum;
	}
}
