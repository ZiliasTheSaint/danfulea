package danfulea.math.numerical;

import java.text.NumberFormat;
import java.util.Locale;

/**
 * Root finding class
 * Based on literature, including Numerical Recipes in C (Cambridge Univ.).
 * @author Dan Fulea, 10 OCT. 2006
 */
public class RootFind {
	private Function func;
	public static boolean failB = false;
	public static String failS = "";

	public static double FACTOR = 1.6;
	public static int NTRY = 50;

	public static double x1_zbrac = 0.0;
	public static double x2_zbrac = 0.0;

	public static int nb_zbrak = 0;

	public static int JMAX = 40;
	public static int MAXIT = 30;
	public static int MAXIT1 = 60;
	public static double UNUSED = -1.11e30;
	public static int ITMAX = 100;// Maximum allowed number of iterations.
	public static double EPS = 3.0e-8;// Machine floating-point precision.
	public static int JMAX1 = 20;
	public static int MAXIT2 = 100;

	public static double EPSS = 1.0e-7;
	public static int MR = 8;
	public static int MT = 10;
	public static int MAXIT3 = MT * MR;
	public static Complex x_laguer;
	public static int its_laguer;

	public static double EPS1 = 2.0e-6;
	public static int MAXM = 100;
	public static double RADIX = 2.0;
	public static int MAXM1 = 50;
	public static int ITMAX1 = 20; // At most ITMAX iterations.
	public static double TINY = 1.0e-6;
	public static double b_qroot = 0.0;
	public static double c_qroot = 0.0;
	private static final double eps = 1e-10;
	public static double EPS2 = 1.0e-4;
	public static int nn_newt = 0;
	public static double[] fvec_newt;// =1.0e-4;

	public static double ALF = 1.0e-4;// Ensures sufficient decrease in function
										// value.
	public static double TOLX = 1.0e-7;// Convergence criterion on ?x.
	public static int check_ln = 0;
	public static double f_ln = 0.0;

	public static int MAXITS = 200;
	public static double TOLF = 1.0e-4;
	public static double TOLMIN = 1.0e-6;
	// #define TOLX 1.0e-7
	public static double STPMX = 100.0;

	// #define MAXITS 200
	// #define EPS 1.0e-7
	// #define TOLF 1.0e-4
	// #define TOLX EPS
	// #define STPMX 100.0
	// #define TOLMIN 1.0e-6

	/**
	 * Constructor. A user supplied function (an object whose class implements the Function interface) is passed here.
	 * @param func the function
	 */
	public RootFind(Function func) {
		this.func = func;
	}

	/**
	 * Dummy constructor.
	 */
	public RootFind()// Function func)
	{
		// for other purpose!
		// this.func=func;
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

	// @#@Bracketing and Bisection
	/*
	 * We will say that a root is bracketed in the interval (a, b) if f(a) and
	 * f(b) have opposite signs. If the function is continuous, then at least
	 * one root must lie in that interval (the intermediate value theorem). If
	 * the function is discontinuous, but bounded, then instead of a root there
	 * might be a step discontinuity which crosses zero (see Figure 9.1.1). For
	 * numerical purposes, that might as well be a root, since the behavior is
	 * indistinguishable from the case of a continuous function whose zero
	 * crossing occurs in between two “adjacent” floating-point numbers in a
	 * machine’s finite-precision representation. Only for functions with
	 * singularities is there the possibility that a bracketed root is not
	 * really there,
	 */
	// public int zbrac(float (*func)(float), float *x1, float *x2)
	/**
	 * Given a function func and an initial guessed range x1 to x2, the routine expands the range 
	 * geometrically until a root is bracketed by the returned values x1_zbrac and x2_zbrac (in which case zbrac 
	 * returns 1) or until the range becomes unacceptably large (in which case zbrac returns 0).
	 * @param x1 x1
	 * @param x2 x2
	 * @return the result
	 */
	public int zbrac(double x1, double x2)
	// Given a function func and an initial guessed range x1 to x2, the routine
	// expands the range
	// geometrically until a root is bracketed by the returned values x1 and x2
	// (in which case zbrac
	// returns 1) or until the range becomes unacceptably large (in which case
	// zbrac returns 0).
	{
		// void nrerror(char error_text[]);
		failB = false;
		int j = 0;
		double f1 = 0.0;
		double f2 = 0.0;

		x1_zbrac = x1;
		x2_zbrac = x2;

		if (x1_zbrac == x2_zbrac) {
			// nrerror("Bad initial range in zbrac");
			failB = true;
			failS = "Bad initial range in zbrac";
			return 0;
		}
		f1 = func.F(x1_zbrac);
		f2 = func.F(x2_zbrac);
		for (j = 1; j <= NTRY; j++) {
			if (f1 * f2 < 0.0)
				return 1;
			if (Math.abs(f1) < Math.abs(f2))
				f1 = func.F(x1_zbrac += FACTOR * (x1_zbrac - x2_zbrac));
			else
				f2 = func.F(x2_zbrac += FACTOR * (x2_zbrac - x1_zbrac));
		}
		return 0;
	}

	/*
	 * Alternatively, you might want to “look inward” on an initial interval,
	 * rather than “look outward” from it, asking if there are any roots of the
	 * function f(x) in the interval from x1 to x2 when a search is carried out
	 * by subdivision into n equal intervals. The following function calculates
	 * brackets for up to nb distinct intervals which each contain one or more
	 * roots.
	 */
	// public void zbrak(float (*fx)(float), float x1, float x2, int n, float
	// xb1[],
	// float xb2[], int *nb)
	/**
	 * Given a function fx defined on the interval from x1-x2 subdivide the interval into n equally 
	 * spaced segments, and search for zero crossings of the function. nb is input as the maximum number 
	 * of roots sought, and is reset to the number of bracketing pairs xb1[1..nb], xb2[1..nb] that are found.
	 * @param x1 x1
	 * @param x2 x2
	 * @param n n
	 * @param xb1 xb1
	 * @param xb2 xb2
	 */
	public void zbrak(double x1, double x2, int n, double[] xb1, double[] xb2)// ,
																				// int
																				// *nb)

	// Given a function fx defined on the interval from x1-x2 subdivide the
	// interval into n equally
	// spaced segments, and search for zero crossings of the function. nb is
	// input as the maximum number
	// of roots sought, and is reset to the number of bracketing pairs
	// xb1[1..nb], xb2[1..nb]
	// that are found.
	{
		int nbb = 0;
		int i = 0;
		double x = 0.0;
		double fp = 0.0;
		double fc = 0.0;
		double dx = 0.0;
		nbb = 0;
		dx = (x2 - x1) / n; // Determine the spacing appropriate to the mesh.
		fp = func.F(x = x1);
		for (i = 1; i <= n; i++) {// Loop over all intervals
			fc = func.F(x += dx);
			if (fc * fp <= 0.0) {// If a sign change occurs then record values
									// for the
									// bounds.
				xb1[++nbb - 1] = x - dx;// xb1[++nbb]=x-dx;
				xb2[nbb - 1] = x;// xb2[nbb]=x;
				if (nb_zbrak == nbb)
					return;
			}
			fp = fc;
		}
		nb_zbrak = nbb;
	}

	// @@@@@@Bisection Method
	/*
	 * Once we know that an interval contains a root, several classical
	 * procedures are available to refine it. The bisection method is one that
	 * cannot fail. It is thus not to be sneered at as a method for otherwise
	 * badly behaved problems. The idea is simple. Over some interval the
	 * function is known to pass through zero because it changes sign. Evaluate
	 * the function at the interval’s midpoint and examine its sign. Use the
	 * midpoint to replace whichever limit has the same sign. After each
	 * iteration the bounds containing the root decrease by a factor of two. If
	 * after n iterations the root is known to be within an interval of size
	 */

	// public double rtbis(float (*func)(float), float x1, float x2, float xacc)
	/**
	 * Using bisection, find the root of a function func known to lie between x1 and x2. The root, 
	 * returned as rtbis, will be refined until its accuracy is ±xacc.
	 * @param x1 x1
	 * @param x2 x2
	 * @param xacc xacc
	 * @return the result
	 */
	public double rtbis(double x1, double x2, double xacc)
	// Using bisection, find the root of a function func known to lie between x1
	// and x2. The root,
	// returned as rtbis, will be refined until its accuracy is ±xacc.
	{
		// void nrerror(char error_text[]);
		failB = false;
		int j = 0;
		double dx = 0.0;
		double f = 0.0;
		double fmid = 0.0;
		double xmid = 0.0;
		double rtb = 0.0;
		f = func.F(x1);
		fmid = func.F(x2);
		if (f * fmid >= 0.0) {
			// nrerror("Root must be bracketed for bisection in rtbis");
			failB = true;
			failS = "Root must be bracketed for bisection in rtbis";
			return 0.0;
		}
		// rtb = f < 0.0 ? (dx=x2-x1,x1) : (dx=x1-x2,x2); //Orient the search so
		// that f>0 lies at x+dx.
		rtb = f < 0.0 ? x1 : x2;
		dx = f < 0.0 ? x2 - x1 : x1 - x2;
		for (j = 1; j <= JMAX; j++) {
			fmid = func.F(xmid = rtb + (dx *= 0.5)); // Bisection loop.
			if (fmid <= 0.0)
				rtb = xmid;
			if (Math.abs(dx) < xacc || fmid == 0.0)
				return rtb;
		}
		// nrerror("Too many bisections in rtbis");
		failB = true;
		failS = "Too many bisections in rtbis";

		return 0.0; // Never get here.
	}

	/*
	 * For functions that are smooth near a root, the methods known respectively
	 * as false position (or regula falsi) and secant method generally converge
	 * faster than bisection. In both of these methods the function is assumed
	 * to be approximately linear in the local region of interest, and the next
	 * improvement in the root is taken as the point where the approximating
	 * line crosses the axis. After each iteration one of the previous boundary
	 * points is discarded in favor of the latest estimate of the root.
	 * Mathematically,the secant method converges more rapidly near a root of a
	 * sufficiently continuous function.
	 * 
	 * The secant method has, however, the disadvantage that the root does not
	 * necessarily remain bracketed. For functions that are not sufficiently
	 * continuous, the algorithm can therefore not be guaranteed to converge:
	 * Local behavior might send it off towards infinity.
	 */
	// public double rtflsp(float (*func)(float), float x1, float x2, float
	// xacc)
	
	/**
	 * Using the false position method, find the root of a function func known to lie between x1 and 
	 * x2. The root, returned as rtflsp, is refined until its accuracy is ±xacc.
	 * @param x1 x1
	 * @param x2 x2
	 * @param xacc xacc
	 * @return the result
	 */
	public double rtflsp(double x1, double x2, double xacc)
	// Using the false position method, find the root of a function func known
	// to lie between x1 and
	// x2. The root, returned as rtflsp, is refined until its accuracy is ±xacc.
	{
		// void nrerror(char error_text[]);
		failB = false;
		int j = 0;
		double fl = 0.0;
		double fh = 0.0;
		double xl = 0.0;
		double xh = 0.0;
		double swap = 0.0;
		double dx = 0.0;
		double del = 0.0;
		double f = 0.0;
		double rtf = 0.0;
		fl = func.F(x1);// fl=(*func)(x1);
		fh = func.F(x2);// Be sure the interval brackets a root.
		if (fl * fh > 0.0) {
			// nrerror("Root must be bracketed in rtflsp");
			failB = true;
			failS = "Root must be bracketed in rtflsp";
			return 0.0;
		}
		if (fl < 0.0) {// Identify the limits so that xl corresponds to the low
						// side.
			xl = x1;
			xh = x2;
		} else {
			xl = x2;
			xh = x1;
			swap = fl;
			fl = fh;
			fh = swap;
		}
		dx = xh - xl;
		for (j = 1; j <= MAXIT; j++) {// False position loop.
			rtf = xl + dx * fl / (fl - fh);// Increment with respect to latest
											// value.
			f = func.F(rtf);
			if (f < 0.0) {// Replace appropriate limit.
				del = xl - rtf;
				xl = rtf;
				fl = f;
			} else {
				del = xh - rtf;
				xh = rtf;
				fh = f;
			}
			dx = xh - xl;
			if (Math.abs(del) < xacc || f == 0.0)
				return rtf; // Convergence.
		}
		failB = true;
		failS = "Maximum number of iterations exceeded in rtflsp";

		// nrerror("Maximum number of iterations exceeded in rtflsp");
		return 0.0;// Never get here.
	}

	// public double rtsec(float (*func)(float), float x1, float x2, float xacc)
	/**
	 * Using the secant method, find the root of a function func thought to lie between x1 and x2. 
	 * The root, returned as rtsec, is refined until its accuracy is ±xacc.
	 * @param x1 x1
	 * @param x2 x2
	 * @param xacc xacc
	 * @return the result
	 */
	public double rtsec(double x1, double x2, double xacc)
	// /Using the secant method, find the root of a function func thought to lie
	// between x1 and x2.
	// The root, returned as rtsec, is refined until its accuracy is ±xacc.
	{
		// void nrerror(char error_text[]);
		failB = false;
		int j = 0;
		double fl = 0.0;
		double f = 0.0;
		double dx = 0.0;
		double swap = 0.0;
		double xl = 0.0;
		double rts = 0.0;
		fl = func.F(x1);
		f = func.F(x2);
		if (Math.abs(fl) < Math.abs(f)) {// Pick the bound with the smaller
											// function value as the most recent
											// guess.
			rts = x1;
			xl = x2;
			swap = fl;
			fl = f;
			f = swap;
		} else {
			xl = x1;
			rts = x2;
		}
		for (j = 1; j <= MAXIT; j++) {// Secant loop.
			dx = (xl - rts) * f / (f - fl); // Increment with respect to latest
											// value.
			xl = rts;
			fl = f;
			rts += dx;
			f = func.F(rts);
			if (Math.abs(dx) < xacc || f == 0.0)
				return rts;// Convergence.
		}
		// nrerror("Maximum number of iterations exceeded in rtsec");
		failB = true;
		failS = "Maximum number of iterations exceeded in rtsec";

		return 0.0; // Never get here.
	}

	/*
	 * A powerful variant on false position is due to Ridders [1]. In both
	 * reliability and speed, Ridders’ method is generally competitive with the
	 * more highly developed and better established (but more complicated)
	 * method ofVanWijngaarden, Dekker, and Brent, which we next discuss.
	 */
	// public double zriddr(float (*func)(float), float x1, float x2, float
	// xacc)
	/**
	 * Using Ridders’ method, return the root of a function func known to lie between x1 and x2. 
	 * The root, returned as zriddr, will be refined to an approximate accuracy xacc.
	 * @param x1 x1
	 * @param x2 x2
	 * @param xacc xacc
	 * @return the result
	 */
	public double zriddr(double x1, double x2, double xacc)
	// Using Ridders’ method, return the root of a function func known to lie
	// between x1 and x2.
	// The root, returned as zriddr, will be refined to an approximate accuracy
	// xacc.
	{
		int j = 0;
		double ans = 0.0;
		double fh = 0.0;
		double fl = 0.0;
		double fm = 0.0;
		double fnew = 0.0;
		double s = 0.0;
		double xh = 0.0;
		double xl = 0.0;
		double xm = 0.0;
		double xnew = 0.0;

		failB = false;

		fl = func.F(x1);
		fh = func.F(x2);
		if ((fl > 0.0 && fh < 0.0) || (fl < 0.0 && fh > 0.0)) {
			xl = x1;
			xh = x2;
			ans = UNUSED;// / Any highly unlikely value, to simplify logic
							// below.
			for (j = 1; j <= MAXIT1; j++) {
				xm = 0.5 * (xl + xh);
				fm = func.F(xm);// First of two function evaluations per
								// iteration.
				s = Math.sqrt(fm * fm - fl * fh);
				if (s == 0.0)
					return ans;
				xnew = xm + (xm - xl) * ((fl >= fh ? 1.0 : -1.0) * fm / s); // Updating
																			// formula.
				if (Math.abs(xnew - ans) <= xacc)
					return ans;
				ans = xnew;
				fnew = func.F(ans); // Second of two function evaluations per
									// iteration.
				if (fnew == 0.0)
					return ans;
				if (SIGN(fm, fnew) != fm) {// Bookkeeping to keep the root
											// bracketed on next iteration.
					xl = xm;
					fl = fm;
					xh = ans;
					fh = fnew;
				} else if (SIGN(fl, fnew) != fl) {
					xh = ans;
					fh = fnew;
				} else if (SIGN(fh, fnew) != fh) {
					xl = ans;
					fl = fnew;
				} else {
					// nrerror("never get here.");
					failB = true;
					failS = "never get here.";
					return 0.0;
				}
				if (Math.abs(xh - xl) <= xacc)
					return ans;
			}
			// nrerror("zriddr exceed maximum iterations");
			failB = true;
			failS = "zriddr exceed maximum iterations";
			return 0.0;
		} else {
			if (fl == 0.0)
				return x1;
			if (fh == 0.0)
				return x2;
			// nrerror("root must be bracketed in zriddr.");
			failB = true;
			failS = "root must be bracketed in zriddr.";
			return 0.0;
		}
		// return 0.0; //Never get here.
	}

	/*
	 * While secant and false position formally converge faster than bisection,
	 * one finds in practice pathological functions for which bisection
	 * converges more rapidly. These can be choppy, discontinuous functions, or
	 * even smooth functions if the second derivative changes sharply near the
	 * root. Bisection always halves the interval, while secant and false
	 * position can sometimes spend many cycles slowly pulling distant bounds
	 * closer to a root. Ridders’ method does a much better job, but it too can
	 * sometimes be fooled. Is there a way to combine superlinear convergence
	 * with the sureness of bisection? Yes. We can keep track of whether a
	 * supposedly superlinear method is actually converging the way it is
	 * supposed to, and, if it is not, we can intersperse bisection steps so as
	 * to guarantee at least linear convergence. This kind of super-strategy
	 * requires attention to bookkeeping detail, and also careful consideration
	 * of how roundoff errors can affect the guiding strategy. Also, we must be
	 * able to determine reliably when convergence has been achieved. An
	 * excellent algorithm that pays close attention to these matters was
	 * developed in the 1960s by van Wijngaarden, Dekker, and others at the
	 * Mathematical Center in Amsterdam, and later improved by Brent [1]. For
	 * brevity, we refer to the final form of the algorithm as Brent’s method.
	 * The method @@@@@@@@@@@@ is guaranteed (by Brent)@@@@@@ to converge, so
	 * long as the function can be evaluated within the initial interval known
	 * to contain a root.
	 */
	// float zbrent(float (*func)(float), float x1, float x2, float tol)
	
	/**
	 * Using Brent’s method, find the root of a function func known to lie between x1 and x2. The 
	 * root, returned as zbrent, will be refined until its accuracy is tol.
	 * @param x1 x1
	 * @param x2 x2
	 * @param tol tol
	 * @return the result
	 */
	public double zbrent(double x1, double x2, double tol)
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

		double fa = func.F(a);
		double fb = func.F(b);
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
		for (iter = 1; iter <= ITMAX; iter++) {
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
			fb = func.F(b);
		}
		// nrerror("Maximum number of iterations exceeded in zbrent");
		failB = true;
		failS = "Maximum number of iterations exceeded in zbrent";

		return 0.0; // Never get here.
	}

	/*
	 * Perhaps the most celebrated of all one-dimensional root-finding routines
	 * is Newton’s method, also called the Newton-Raphson method. This method is
	 * distinguished from the methods of previous sections by the fact that it
	 * requires the evaluation of both the function f(x), and the derivative f
	 * '(x), at arbitrary points x. Far from a root, where the higher-order
	 * terms in the series are important, the Newton-Raphson formula can give
	 * grossly inaccurate, meaningless corrections. For instance, the initial
	 * guess for the root might be so far from the true root as to let the
	 * search interval include a local maximum or minimum of the function. This
	 * can be death to the method Like most powerful tools, Newton-Raphson can
	 * be destructive used in inappropriate circumstances. Why do we call
	 * Newton-Raphson powerful? The answer lies in its rate of convergence: This
	 * very strong convergence property makes Newton-Raphson the method of
	 * choice for any function whose derivative can be evaluated efficiently,
	 * and whose derivative is continuous and nonzero in the neighborhood of a
	 * root. For an efficient realization of Newton-Raphson the user provides a
	 * routine that evaluates both f(x) and its first derivative f (x) at the
	 * point x. The Newton-Raphson formula can also be applied using a numerical
	 * difference to approximate the true local derivative, This is not,
	 * however, a recommended procedure for the following reasons: (i) You are
	 * doing two function evaluations per step, so at best the superlinear order
	 * of convergence will be only ?2. (ii) If you take dx too small you will be
	 * wiped out by roundoff, while if you take it too large your order of
	 * convergence will be only linear, no better than using the initial
	 * evaluation f (x0) for all subsequent steps. Therefore, Newton-Raphson
	 * with numerical derivatives is (in one dimension) always dominated by the
	 * secant method of §9.2. (In multidimensions, where there is a paucity of
	 * available methods, Newton-Raphson with numerical derivatives must be
	 * taken more seriously. The following function calls a user supplied
	 * function funcd(x,fn,df) which supplies the function value as fn and the
	 * derivative as df. We have included input bounds on the root simply to be
	 * consistent with previous root-finding routines: Newton does not adjust
	 * bounds, and works only on local information at the point x. The bounds
	 * are used only to pick the midpoint as the first guess, and to reject the
	 * solution if it wanders outside of the bounds.
	 */
	// float rtnewt(void (*funcd)(float, float *, float *), float x1, float
	// x2,float xacc)
	
	/**
	 * Using the Newton-Raphson method, find the root of a function known to lie in the interval 
	 * [x1, x2]. The root rtnewt will be refined until its accuracy is known within ±xacc. funcd 
	 * is a user-supplied routine that returns both the function value and the first derivative of the 
	 * function at the point x.
	 * @param x1 x1
	 * @param x2 x2
	 * @param xacc xacc
	 * @return the result
	 */
	public double rtnewt(double x1, double x2, double xacc)
	// Using the Newton-Raphson method, find the root of a function known to lie
	// in the interval
	// [x1, x2]. The root rtnewt will be refined until its accuracy is known
	// within ±xacc. funcd
	// is a user-supplied routine that returns both the function value and the
	// first derivative of the
	// function at the point x.
	{
		// void nrerror(char error_text[]);
		failB = false;
		int j = 0;
		double df = 0.0;
		double dx = 0.0;
		double f = 0.0;
		double rtn = 0.0;

		rtn = 0.5 * (x1 + x2); // Initial guess.
		for (j = 1; j <= JMAX1; j++) {
			// (*funcd)(rtn,&f,&df);
			double[] ffd = func.FD(rtn);
			f = ffd[0];
			df = ffd[1];
			dx = f / df;
			rtn -= dx;
			if ((x1 - rtn) * (rtn - x2) < 0.0) {
				// nrerror("Jumped out of brackets in rtnewt");
				failB = true;
				failS = "Jumped out of brackets in rtnewt";

				return 0.0;
			}
			if (Math.abs(dx) < xacc)
				return rtn;// Convergence.
		}
		// nrerror("Maximum number of iterations exceeded in rtnewt");
		failB = true;
		failS = "Maximum number of iterations exceeded in rtnewt";

		return 0.0; // Never get here.
	}

	/*
	 * While Newton-Raphson’s global convergence properties are poor, it is
	 * fairly easy to design a @@@@@fail-safe@@@ routine that utilizes a
	 * combination of bisection and Newton- Raphson. The hybrid algorithm takes
	 * a bisection step whenever Newton-Raphson would take the solution out of
	 * bounds, or whenever Newton-Raphson is not reducing the size of the
	 * brackets rapidly enough.
	 */
	// float rtsafe(void (*funcd)(float, float *, float *), float x1, float
	// x2,float xacc)
	
	/**
	 * Using a combination of Newton-Raphson and bisection, find the root of a function bracketed 
	 * between x1 and x2. The root, returned as the function value rtsafe, will be refined until 
	 * its accuracy is known within ±xacc. funcd is a user-supplied routine that returns both the 
	 * function value and the first derivative of the function.
	 * @param x1 x1
	 * @param x2 x2
	 * @param xacc xacc
	 * @return the result
	 */
	public double rtsafe(double x1, double x2, double xacc)
	// Using a combination of Newton-Raphson and bisection, find the root of a
	// function bracketed
	// between x1 and x2. The root, returned as the function value rtsafe, will
	// be refined until
	// its accuracy is known within ±xacc. funcd is a user-supplied routine that
	// returns both the
	// function value and the first derivative of the function.
	{
		// void nrerror(char error_text[]);
		failB = false;
		int j = 0;
		double df = 0.0;
		double dx = 0.0;
		double dxold = 0.0;
		double f = 0.0;
		double fh = 0.0;
		double fl = 0.0;
		double temp = 0.0;
		double xh = 0.0;
		double xl = 0.0;
		double rts = 0.0;
		// (*funcd)(x1,&fl,&df);
		double[] ffd = func.FD(x1);
		fl = ffd[0];
		df = ffd[1];
		// (*funcd)(x2,&fh,&df);
		ffd = func.FD(x2);
		fh = ffd[0];
		df = ffd[1];
		if ((fl > 0.0 && fh > 0.0) || (fl < 0.0 && fh < 0.0)) {
			// nrerror("Root must be bracketed in rtsafe");
			failB = true;
			failS = "Root must be bracketed in rtsafe";

			return 0.0;
		}
		if (fl == 0.0)
			return x1;
		if (fh == 0.0)
			return x2;
		if (fl < 0.0) {// Orient the search so that f(xl) < 0.
			xl = x1;
			xh = x2;
		} else {
			xh = x1;
			xl = x2;
		}
		rts = 0.5 * (x1 + x2);// Initialize the guess for root,
		dxold = Math.abs(x2 - x1);// the “stepsize before last,”
		dx = dxold; // and the last step.
		// (*funcd)(rts,&f,&df);
		ffd = func.FD(rts);
		f = ffd[0];
		df = ffd[1];
		for (j = 1; j <= MAXIT2; j++) {// Loop over allowed iterations.
			if ((((rts - xh) * df - f) * ((rts - xl) * df - f) > 0.0) // Bisect
																		// if
																		// Newton
																		// out
																		// of
																		// range,
					|| (Math.abs(2.0 * f) > Math.abs(dxold * df))) { // or not
																		// decreasing
																		// fast
																		// enough.
				dxold = dx;
				dx = 0.5 * (xh - xl);
				rts = xl + dx;
				if (xl == rts)
					return rts;// Change in root is negligible.
			} else {// Newton step acceptable. Take it.
				dxold = dx;
				dx = f / df;
				temp = rts;
				rts -= dx;
				if (temp == rts)
					return rts;
			}
			if (Math.abs(dx) < xacc)
				return rts;// Convergence criterion.
			// (*funcd)(rts,&f,&df);
			ffd = func.FD(rts);
			f = ffd[0];
			df = ffd[1];
			// The one new function evaluation per iteration.
			if (f < 0.0) // Maintain the bracket on the root.
				xl = rts;
			else
				xh = rts;
		}
		// nrerror("Maximum number of iterations exceeded in rtsafe");
		failB = true;
		failS = "Maximum number of iterations exceeded in rtsafe";

		return 0.0;// Never get here.
	}

	/*
	 * Roots of Polynomials Newton-Raphson may work, but slowly, since large
	 * roundoff errors can occur. When a root is known in advance to be
	 * multiple, then special methods of attack are readily devised. Problems
	 * arise when (as is generally the case) we do not know in advance what
	 * pathology a root will display. When seeking several or all roots of a
	 * polynomial, the total effort can be significantly reduced by the use of
	 * deflation. There are two schools of thought about how to proceed when
	 * faced with a polynomial of real coefficients. One school says to go after
	 * the easiest quarry, the real, distinct roots, by the same kinds of
	 * methods that we have discussed in previous sections for general
	 * functions, i.e., trial-and-error bracketing followed by a safe
	 * Newton-Raphson as in rtsafe. Sometimes you are only interested in real
	 * roots, in which case the strategy is complete. Otherwise, you then go
	 * after quadratic factors of the form (9.5.1) by any of a variety of
	 * methods. One such is Bairstow’s method, which we will discuss below in
	 * the context of root polishing. Another is Muller’s method, which we here
	 * briefly discuss.
	 * 
	 * Laguerre’s method is by far the most straightforward of these general,
	 * complex methods. It does require complex arithmetic, even while
	 * converging to real roots; however, for polynomials with all real roots,
	 * it is guaranteed to converge to a root from any starting point. For
	 * polynomials with some complex roots, little is theoretically proved about
	 * the method’s convergence.
	 */

	// public void laguer(fcomplex a[], int m, fcomplex *x, int *its)
	/**
	 * Given the degree m and the m+1 complex coefficients a[0..m] of the polynomial sum from 
	 * i=0 to m of a[i]x^i, and given a complex value x, this routine improves x by Laguerre’s method 
	 * until it converges, within the achievable roundoff limit, to a root of the given polynomial. 
	 * The number of iterations taken is returned as its.
	 * @param a a
	 * @param m m
	 * @param x x
	 */
	public void laguer(Complex[] a, int m, Complex x)// , int *its)
	// Given the degree m and the m+1 complex coefficients a[0..m] of the
	// polynomial sum from
	// i=0 to m of a[i]xi,
	// and given a complex value x, this routine improves x by Laguerre’s method
	// until it converges,
	// within the achievable roundoff limit, to a root of the given polynomial.
	// The number of iterations taken is returned as its.
	{
		int iter = 0;
		int j = 0;
		double abx = 0.0;
		double abp = 0.0;
		double abm = 0.0;
		double err = 0.0;
		Complex dx, x1, b, d, f, g, h, sq, gp, gm, g2;
		// static float frac[MR+1] =
		// {0.0,0.5,0.25,0.75,0.13,0.38,0.62,0.88,1.0};
		double[] frac = { 0.0, 0.5, 0.25, 0.75, 0.13, 0.38, 0.62, 0.88, 1.0 };
		// Fractions used to break a limit cycle.
		failB = false;
		x_laguer = x;

		for (iter = 1; iter <= MAXIT3; iter++) {// Loop over iterations up to
												// allowed maximum.
			its_laguer = iter;
			b = a[m];
			err = Complex.Cabs(b);
			d = f = new Complex(0.0, 0.0);
			abx = Complex.Cabs(x_laguer);
			for (j = m - 1; j >= 0; j--) {// Efficient computation of the
											// polynomial and
											// its first two derivatives. f
											// stores P''/2.
				f = Complex.Cadd(Complex.Cmul(x_laguer, f), d);
				d = Complex.Cadd(Complex.Cmul(x_laguer, d), b);
				b = Complex.Cadd(Complex.Cmul(x_laguer, b), a[j]);
				err = Complex.Cabs(b) + abx * err;
			}
			err *= EPSS;
			// Estimate of roundoff error in evaluating polynomial.
			if (Complex.Cabs(b) <= err)
				return;// We are on the root.
			g = Complex.Cdiv(d, b);// The generic case: use Laguerre’s formula.
			g2 = Complex.Cmul(g, g);
			h = Complex.Csub(g2, Complex.RCmul(2.0, Complex.Cdiv(f, b)));
			sq = Complex.Csqrt(Complex.RCmul((double) (m - 1),
					Complex.Csub(Complex.RCmul((double) m, h), g2)));
			gp = Complex.Cadd(g, sq);
			gm = Complex.Csub(g, sq);
			abp = Complex.Cabs(gp);
			abm = Complex.Cabs(gm);
			if (abp < abm)
				gp = gm;
			dx = ((Math.max(abp, abm) > 0.0 ? Complex.Cdiv(new Complex(
					(double) m, 0.0), gp) : Complex.RCmul(1 + abx, new Complex(
					Math.cos((double) iter), Math.sin((double) iter)))));
			x1 = Complex.Csub(x_laguer, dx);
			// if (x->r == x1.r && x->i == x1.i) return; //Converged.
			if (x_laguer.r == x1.r && x_laguer.i == x1.i)
				return;
			// if (iter % MT) x_laguer=x1;
			if (iter % MT != 0)
				x_laguer = x1;
			else
				x_laguer = Complex.Csub(x_laguer,
						Complex.RCmul(frac[iter / MT], dx));
			// Every so often we take a fractional step, to break any limit
			// cycle (itself a rare occurrence).
		}
		// nrerror("too many iterations in laguer");
		// Very unusual — can occur only for complex roots. Try a different
		// starting guess for the
		// root.
		failB = true;
		failS = "too many iterations in laguer";

		return;
	}

	/*
	 * Here is a driver routine that calls laguer in succession for each root,
	 * performs the deflation, optionally polishes the roots by the same
	 * Laguerre method — if you are not going to polish in some other way — and
	 * finally sorts the roots by their real parts. (We will use this routine in
	 * Chapter 13.)
	 */
	/**
	 * Given the degree m and the m+1 complex coefficients a[0..m] of the polynomial 
	 * sum from i=0 to m of a(i)x^i, this routine successively calls laguer and finds all m complex roots in 
	 * roots[1..m]. The boolean variable polish should be input as true (1) if polishing (also by 
	 * Laguerre’s method) is desired, false (0) if the roots will be subsequently polished by other 
	 * means.
	 * @param a a
	 * @param m m
	 * @param roots roots
	 * @param polish polish
	 */
	public void zroots(Complex[] a, int m, Complex[] roots, int polish)
	// Given the degree m and the m+1 complex coefficients a[0..m] of the
	// polynomial
	// sum from i=0 to m of a(i)xi,
	// this routine successively calls laguer and finds all m complex roots in
	// roots[1..m]. The
	// boolean variable polish should be input as true (1) if polishing (also by
	// Laguerre’s method)
	// is desired, false (0) if the roots will be subsequently polished by other
	// means.
	{
		// void laguer(fcomplex a[], int m, fcomplex *x, int *its);
		int i = 0;
		// int its = 0;
		int j = 0;
		int jj = 0;
		Complex x, b, c;
		Complex[] ad = new Complex[MAXM];

		for (j = 0; j <= m; j++)
			ad[j] = a[j];// Copy of coefficients for successive deflation.
		for (j = m; j >= 1; j--) {// Loop over each root to be found.
			x = new Complex(0.0, 0.0);// Start at zero to favor convergence to
										// smallest
										// remaining root, and find the root.
			laguer(ad, j, x);// ,&x,&its);
			x = x_laguer;
			// its = its_laguer;
			if (Math.abs(x.i) <= 2.0 * EPS1 * Math.abs(x.r))
				x.i = 0.0;
			roots[j - 1] = x;// roots[j]=x;
			b = ad[j];// Forward deflation.
			for (jj = j - 1; jj >= 0; jj--) {
				c = ad[jj];
				ad[jj] = b;
				b = Complex.Cadd(Complex.Cmul(x, b), c);
			}
		}
		if (polish == 1)// (polish)
			for (j = 1; j <= m; j++)// Polish the roots using the undeflated
									// coefficients.
			{
				laguer(a, m, roots[j - 1]);// ,&roots[j],&its);
				roots[j - 1] = x_laguer;
				// its = its_laguer;
			}

		for (j = 2; j <= m; j++) {// Sort roots by their real parts by straight
									// insertion.
			x = roots[j - 1];// x=roots[j];
			for (i = j - 1; i >= 1; i--) {
				if (roots[i - 1].r <= x.r)
					break;// if (roots[i].r <= x.r) break;
				roots[i] = roots[i - 1];// roots[i+1]=roots[i];
			}
			roots[i] = x;// roots[i+1]=x;
		}
	}

	// ====================================
	/**
	 * Given a matrix a[1..n][1..n], this routine replaces it by a balanced matrix with identical 
	 * eigenvalues. A symmetric matrix is already balanced and is unaffected by this procedure. The 
	 * parameter RADIX should be the machine’s floating-point radix.
	 * @param a a
	 * @param n n
	 */
	public void balanc(double[][] a, int n)
	// Given a matrix a[1..n][1..n], this routine replaces it by a balanced
	// matrix with identical
	// eigenvalues. A symmetric matrix is already balanced and is unaffected by
	// this procedure. The
	// parameter RADIX should be the machine’s floating-point radix.
	{
		int last = 0;
		int j = 0;
		int i = 0;
		double s = 0.0;
		double r = 0.0;
		double g = 0.0;
		double f = 0.0;
		double c = 0.0;
		double sqrdx = 0.0;

		sqrdx = RADIX * RADIX;
		last = 0;
		while (last == 0) {
			last = 1;
			for (i = 1; i <= n; i++) {// Calculate row and column norms.
				r = c = 0.0;
				for (j = 1; j <= n; j++)
					if (j != i) {
						c += Math.abs(a[j - 1][i - 1]);// c += fabs(a[j][i]);
						r += Math.abs(a[i - 1][j - 1]);// r += fabs(a[i][j]);
					}
				if (c != 0.0 && r != 0.0)// if (c && r)
				{// If both are nonzero,
					g = r / RADIX;
					f = 1.0;
					s = c + r;
					while (c < g) {// find the integer power of the machine
									// radix that
									// comes closest to balancing the matrix.
						f *= RADIX;
						c *= sqrdx;
					}
					g = r * RADIX;
					while (c > g) {
						f /= RADIX;
						c /= sqrdx;
					}
					if ((c + r) / f < 0.95 * s) {
						last = 0;
						g = 1.0 / f;
						for (j = 1; j <= n; j++)
							a[i - 1][j - 1] *= g;// a[i][j] *= g; //Apply
													// similarity
													// transformation.
						for (j = 1; j <= n; j++)
							a[j - 1][i - 1] *= f;// a[j][i] *= f;
					}
				}
			}
		}
	}

	/**
	 * Finds all eigenvalues of an upper Hessenberg matrix a[1..n][1..n]. On input a can be 
	 * exactly as output from elmhes (see Eigensystem), on output it is destroyed. The real and imaginary parts 
	 * of the eigenvalues are returned in wr[1..n] and wi[1..n], respectively.
	 * @param a a
	 * @param n n
	 * @param wr wr
	 * @param wi wi
	 */
	public void hqr(double[][] a, int n, double wr[], double wi[])
	// Finds all eigenvalues of an upper Hessenberg matrix a[1..n][1..n]. On
	// input a can be
	// exactly as output from elmhes §11.5; on output it is destroyed. The real
	// and imaginary parts
	// of the eigenvalues are returned in wr[1..n] and wi[1..n], respectively.
	{
		int nn = 0;
		int m = 0;
		int l = 0;
		int k = 0;
		int j = 0;
		int its = 0;
		int i = 0;
		int mmin = 0;
		double z = 0.0;
		double y = 0.0;
		double x = 0.0;
		double w = 0.0;
		double v = 0.0;
		double u = 0.0;
		double t = 0.0;
		double s = 0.0;
		double r = 0.0;
		double q = 0.0;
		double p = 0.0;
		double anorm = 0.0;

		failB = false;

		anorm = 0.0;// Compute matrix norm for possible use in locating
					// single small subdiagonal element.
		for (i = 1; i <= n; i++)
			for (j = Math.max(i - 1, 1); j <= n; j++)
				anorm += Math.abs(a[i - 1][j - 1]);// anorm += fabs(a[i][j]);
		nn = n;
		t = 0.0; // Gets changed only by an exceptional shift.
		while (nn >= 1) {// Begin search for next eigenvalue.
			its = 0;
			do {
				for (l = nn; l >= 2; l--) {// Begin iteration: look for single
											// small subdiagonal element.
					s = Math.abs(a[l - 2][l - 2]) + Math.abs(a[l - 1][l - 1]);// s=fabs(a[l-1][l-1])+fabs(a[l][l]);
					if (s == 0.0)
						s = anorm;
					if ((double) (Math.abs(a[l - 1][l - 2]) + s) == s)// ((double)(Math.abs(a[l][l-1])
																		// + s)
																		// == s)
					{
						a[l - 1][l - 2] = 0.0;// a[l][l-1]=0.0;
						break;
					}
				}
				x = a[nn - 1][nn - 1];// x=a[nn][nn];
				if (l == nn) {// One root found.
					wr[nn - 1] = x + t;// wr[nn]=x+t;
					wi[nn-- - 1] = 0.0;// wi[nn--]=0.0;
				} else {
					y = a[nn - 2][nn - 2];// y=a[nn-1][nn-1];
					w = a[nn - 1][nn - 2] * a[nn - 2][nn - 1];// w=a[nn][nn-1]*a[nn-1][nn];
					if (l == (nn - 1)) {// Two roots found...
						p = 0.5 * (y - x);
						q = p * p + w;
						z = Math.sqrt(Math.abs(q));
						x += t;
						if (q >= 0.0) {// ...a real pair.
							z = p + SIGN(z, p);
							wr[nn - 2] = wr[nn - 1] = x + z;// wr[nn-1]=wr[nn]=x+z;
							if (z != 0.0)
								wr[nn - 1] = x - w / z;// wr[nn]=x-w/z;
							wi[nn - 2] = wi[nn - 1] = 0.0;// wi[nn-1]=wi[nn]=0.0;
						} else {// ...a complex pair.
							wr[nn - 2] = wr[nn - 1] = x + p;// wr[nn-1]=wr[nn]=x+p;
							wi[nn - 2] = -(wi[nn - 1] = z);// wi[nn-1]=
															// -(wi[nn]=z);
						}
						nn -= 2;
					} else {// No roots found. Continue iteration.
						if (its == 30) {
							// nrerror("Too many iterations in hqr");
							failB = true;
							failS = "Too many iterations in hqr";

							return;
						}
						if (its == 10 || its == 20) {// Form exceptional shift.
							t += x;
							for (i = 1; i <= nn; i++)
								a[i - 1][i - 1] -= x;// a[i][i] -= x;
							// s=fabs(a[nn][nn-1])+fabs(a[nn-1][nn-2]);
							s = Math.abs(a[nn - 1][nn - 2])
									+ Math.abs(a[nn - 2][nn - 3]);
							y = x = 0.75 * s;
							w = -0.4375 * s * s;
						}
						++its;
						for (m = (nn - 2); m >= l; m--) {// Form shift and then
															// look for 2
															// consecutive small
															// subdiagonal
															// elements.
							z = a[m - 1][m - 1];// z=a[m][m];
							r = x - z;
							s = y - z;
							// p=(r*s-w)/a[m+1][m]+a[m][m+1]; Equation
							// (11.6.23).
							p = (r * s - w) / a[m][m - 1] + a[m - 1][m];
							q = a[m][m] - z - r - s;// q=a[m+1][m+1]-z-r-s;
							r = a[m + 1][m];// r=a[m+2][m+1];
							s = Math.abs(p) + Math.abs(q) + Math.abs(r); // Scale
																			// to
																			// prevent
																			// overflow
																			// or
																			// underflow.
							p /= s;
							q /= s;
							r /= s;
							if (m == l)
								break;
							// u=fabs(a[m][m-1])*(fabs(q)+fabs(r));
							u = Math.abs(a[m - 1][m - 2])
									* (Math.abs(q) + Math.abs(r));
							// v=fabs(p)*(fabs(a[m-1][m-1])+fabs(z)+fabs(a[m+1][m+1]));
							v = Math.abs(p)
									* (Math.abs(a[m - 2][m - 2]) + Math.abs(z) + Math
											.abs(a[m][m]));
							if ((double) (u + v) == v)
								break; // Equation (11.6.26).
						}
						for (i = m + 2; i <= nn; i++) {
							a[i - 1][i - 3] = 0.0;// a[i][i-2]=0.0;
							if (i != (m + 2))
								a[i - 1][i - 4] = 0.0;// a[i][i-3]=0.0;
						}
						for (k = m; k <= nn - 1; k++) {
							// Double QR step on rows l to nn and columns m to
							// nn.
							if (k != m) {
								p = a[k - 1][k - 2];// p=a[k][k-1]; Begin setup
													// of Householder vector.
								q = a[k][k - 2];// q=a[k+1][k-1];
								r = 0.0;
								if (k != (nn - 1))
									r = a[k + 1][k - 2];// r=a[k+2][k-1];
								if ((x = Math.abs(p) + Math.abs(q)
										+ Math.abs(r)) != 0.0) {
									p /= x; // Scale to prevent overflow or
											// underflow.
									q /= x;
									r /= x;
								}
							}
							if ((s = SIGN(Math.sqrt(p * p + q * q + r * r), p)) != 0.0) {
								if (k == m) {
									if (l != m)
										a[k - 1][k - 2] = -a[k - 1][k - 2];// a[k][k-1]
																			// =
																			// -a[k][k-1];
								} else
									a[k - 1][k - 2] = -s * x;// a[k][k-1] =
																// -s*x;
								p += s;// Equations (11.6.24).
								x = p / s;
								y = q / s;
								z = r / s;
								q /= p;
								r /= p;
								for (j = k; j <= nn; j++) {// Row modification.
									p = a[k - 1][j - 1] + q * a[k][j - 1];// p=a[k][j]+q*a[k+1][j];
									if (k != (nn - 1)) {
										p += r * a[k + 1][j - 1];// p +=
																	// r*a[k+2][j];
										a[k + 1][j] -= p * z;// a[k+2][j] -=
																// p*z;
									}
									a[k][j - 1] -= p * y;// a[k+1][j] -= p*y;
									a[k - 1][j - 1] -= p * x;// a[k][j] -= p*x;
								}
								mmin = nn < k + 3 ? nn : k + 3;
								for (i = l; i <= mmin; i++) {// Column
																// modification.
									p = x * a[i - 1][k - 1] + y * a[i - 1][k];// p=x*a[i][k]+y*a[i][k+1];
									if (k != (nn - 1)) {
										p += z * a[i - 1][k + 1];// p +=
																	// z*a[i][k+2];
										a[i - 1][k + 1] -= p * r;// a[i][k+2] -=
																	// p*r;
									}
									a[i - 1][k] -= p * q;// a[i][k+1] -= p*q;
									a[i - 1][k - 1] -= p;// a[i][k] -= p;
								}
							}
						}
					}
				}
			} while (l < nn - 1);
		}
	}

	/*
	 * This method, implemented in the routine zrhqr following, is typically
	 * about a factor 2 slower than zroots (above). However, for some classes of
	 * polynomials, it is a more robust technique, largely because of the fairly
	 * sophisticated convergence methods embodied in hqr. If your polynomial has
	 * real coefficients, and you are having trouble with zroots, then zrhqr is
	 * a recommended alternative.
	 */
	/**
	 * Find all the roots of a polynomial with real coefficients, sum from i=0 to m of a(i)x^i, given the degree m 
	 * and the coefficients a[0..m]. The method is to construct an upper Hessenberg matrix whose 
	 * eigenvalues are the desired roots, and then use the routines balanc and hqr. The real and 
	 * imaginary parts of the roots are returned in rtr[1..m] and rti[1..m], respectively.
	 * @param a a
	 * @param m m
	 * @param rtr rtr
	 * @param rti rti
	 */
	public void zrhqr(double[] a, int m, double[] rtr, double[] rti)
	// Find all the roots of a polynomial with real coefficients,
	// sum from i=0 to m of a(i)x^i, given the degree m
	// and the coefficients a[0..m]. The method is to construct an upper
	// Hessenberg matrix whose
	// eigenvalues are the desired roots, and then use the routines balanc and
	// hqr. The real and
	// imaginary parts of the roots are returned in rtr[1..m] and rti[1..m],
	// respectively.
	{
		// void balanc(float **a, int n);
		// void hqr(float **a, int n, float wr[], float wi[]);
		int j = 0;
		int k = 0;
		double[][] hess = new double[MAXM1][MAXM1];
		double xr = 0.0;
		double xi = 0.0;

		failB = false;
		// hess=matrix(1,MAXM,1,MAXM);
		if (m > MAXM1 || a[m] == 0.0) {
			// nrerror("bad args in zrhqr");
			failB = true;
			failS = "bad args in zrhqr";
			return;
		}
		for (k = 1; k <= m; k++) {// Construct the matrix.
			hess[0][k - 1] = -a[m - k] / a[m];// hess[1][k] = -a[m-k]/a[m];
			for (j = 2; j <= m; j++)
				hess[j - 1][k - 1] = 0.0;// hess[j][k]=0.0;
			if (k != m)
				hess[k][k - 1] = 1.0;// hess[k+1][k]=1.0;
		}
		balanc(hess, m); // Find its eigenvalues.
		hqr(hess, m, rtr, rti);
		if (failB)
			return;// @@@@@@@@@@@@@@@@@@@@@@@@
		for (j = 2; j <= m; j++) {// Sort roots by their real parts by straight
									// insertion.
			xr = rtr[j - 1];// rtr[j];
			xi = rti[j - 1];// rti[j];
			for (k = j - 1; k >= 1; k--) {
				if (rtr[k - 1] <= xr)
					break;// if (rtr[k] <= xr) break;
				rtr[k] = rtr[k - 1];// rtr[k+1]=rtr[k];
				rti[k] = rti[k - 1];// rti[k+1]=rti[k];
			}
			rtr[k] = xr;// rtr[k+1]=xr;
			rti[k] = xi;// rti[k+1]=xi;
		}
		// free_matrix(hess,1,MAXM,1,MAXM);
	}

	// pe baza coeficientilor functiei de evaluat (polinom de interpolare) (a)
	// se returneza solutia in tabloul double[].
	// depinzand de coeficienti,si de ordinul polinomului de interpolare se
	// construieste xx!!
	/**
	 * Given the coefficients of the polynomial (usually used in interpolation), 
	 * this routine calculates the roots using Birge-Viette method.
	 * @param a a
	 * @return the result
	 */
	public static double[] birge(double[] a) {
		double x = 0.0;
		double d = 0.0;
		double p = 0.0;
		failB = false;// true;
		int nx = 0;
		int i = 0;
		int n = a.length - 1;// grad polinom
		// daca am 3 coeficienti<-> polinomul de interpolare e de grad 2 si am
		// 3-1 solutii!!!
		double[] xx = new double[n];
		int imax = 100;
		if (n > 1) {
			for (int m = n; m >= 2; --m) {
				i = 0;
				while (!failB)// (zero)
				{
					i++;
					d = a[0];
					p = a[0] * x + a[1];
					for (int j = 2; j <= m; j++) {
						d = d * x + p;
						p = p * x + a[j];
					}
					if (d != 0)
						d = p / d;
					else
						d = p;
					x = x - d;// noua aproximare a radacinii
					if (x != 0)
						d = d / x;
					if (Math.abs(d) <= eps || i >= imax)
						break;
				}
				if (i >= imax) {
					xx[nx] = Double.NaN;
					return trunc(xx);// mereu sunt gasite solutiile reale!!!
				}
				nx++;
				xx[nx - 1] = x;
				for (int j = 1; j <= m - 1; j++)
					a[j] = a[j - 1] * x + a[j];
			}
		}
		nx = n;
		xx[n - 1] = -a[1] / a[0];

		return xx;
	}

	/**
	 * Used internally.
	 * @param xx xx
	 * @return the result
	 */
	private static double[] trunc(double[] xx) {
		int j = -1;
		for (int i = 0; i < xx.length; i++) {
			j++;
			if (Double.toString(xx[i]).compareTo("NaN") == 0)
				break;
		}
		double[] x = new double[j];
		for (int i = 0; i < j; i++)
			x[i] = xx[i];

		return x;
	}

	/*
	 * Techniques for Root-Polishing
	 * 
	 * Newton-Raphson works very well for real roots once the neighborhood of a
	 * root has been identified. The polynomial and its derivative can be
	 * efficiently simultaneously evaluated as in §5.3. For a polynomial of
	 * degree n with coefficients c[0]...c[n], the following segment of code
	 * embodies one cycle of Newton- Raphson: p=c[n]*x+c[n-1]; p1=c[n];
	 * for(i=n-2;i>=0;i--) { p1=p+p1*x; p=c[i]+p*x; } if (p1 == 0.0)
	 * nrerror("derivative should not vanish"); x -= p/p1
	 * 
	 * Bairstow’s method now consists of using Newton-Raphson in two dimensions
	 * (which is actually the subject of the next section) to find a
	 * simultaneous zero of R and S. Synthetic division is used twice per cycle
	 * to evaluate R, S and their partial derivatives with respect to B,C. Like
	 * one-dimensional Newton-Raphson, the method works well in the vicinity of
	 * a root pair (real or complex), but it can fail miserably when started at
	 * a random point. We therefore recommend it only in the context of
	 * polishing tentative complex roots.
	 */
	/**
	 * Given n+1 coefficients p[0..n] of a polynomial of degree n, and trial values for the coefficients 
	 * of a quadratic factor x*x+b*x+c, improve the solution until the coefficients b,c change by less 
	 * than eps. The routine poldiv (see EvalFunc class) is used.
	 * @param p p
	 * @param n n
	 * @param b b
	 * @param c c
	 * @param eps eps
	 */
	public void qroot(double[] p, int n, double b, double c, double eps)
	// Given n+1 coefficients p[0..n] of a polynomial of degree n, and trial
	// values for the coefficients
	// of a quadratic factor x*x+b*x+c, improve the solution until the
	// coefficients b,c change by less
	// than eps. The routine poldiv §5.3 is used.
	{
		// void poldiv(float u[], int n, float v[], int nv, float q[], float
		// r[]);EvalFunc
		int iter = 0;
		double sc = 0.0;
		double sb = 0.0;
		double s = 0.0;
		double rc = 0.0;
		double rb = 0.0;
		double r = 0.0;
		double dv = 0.0;
		double delc = 0.0;
		double delb = 0.0;
		double[] q = new double[n + 1];
		double[] qq = new double[n + 1];
		double[] rem = new double[n + 1];
		double[] d = new double[3];// float d[3];
		// q=vector(0,n);
		// qq=vector(0,n);
		// rem=vector(0,n);
		b_qroot = b;
		c_qroot = c;
		failB = false;

		d[2] = 1.0;
		for (iter = 1; iter <= ITMAX1; iter++) {
			d[1] = b_qroot;
			d[0] = c_qroot;
			EvalFunc.poldiv(p, n, d, 2, q, rem);
			s = rem[0];// First division r,s.
			r = rem[1];
			EvalFunc.poldiv(q, (n - 1), d, 2, qq, rem);
			sb = -c_qroot * (rc = -rem[1]);// Second division partial r,s with
											// respect to c.
			rb = -b_qroot * rc + (sc = -rem[0]);
			dv = 1.0 / (sb * rc - sc * rb);// Solve 2x2 equation.
			delb = (r * sc - s * rc) * dv;
			delc = (-r * sb + s * rb) * dv;
			b_qroot += (delb = (r * sc - s * rc) * dv);
			c_qroot += (delc = (-r * sb + s * rb) * dv);
			if ((Math.abs(delb) <= eps * Math.abs(b_qroot) || Math.abs(b_qroot) < TINY)
					&& (Math.abs(delc) <= eps * Math.abs(c_qroot) || Math
							.abs(c_qroot) < TINY)) {
				// free_vector(rem,0,n);// Coefficients converged.
				// free_vector(qq,0,n);
				// free_vector(q,0,n);
				return;
			}
		}
		// nrerror("Too many iterations in routine qroot");
		failB = true;
		failS = "Too many iterations in routine qroot";
		return;
	}

	// Newton-Raphson Method for Nonlinear Systems of Equations
	/*
	 * We make an extreme, but wholly defensible, statement: There are no good,
	 * general methods for solving systems of more than one nonlinear equation.
	 * Furthermore, it is not hard to see why (very likely) there never will be
	 * any good, general methods: Consider the case of two dimensions, where we
	 * want to solve simultaneously f(x, y) = 0 g(x, y) = 0 The functions f and
	 * g are two arbitrary functions, each of which has zero contour lines that
	 * divide the (x, y) plane into regions where their respective function is
	 * positive or negative. These zero contour boundaries are of interest to
	 * us. The solutions that we seek are those points (if any) that are common
	 * to the zero contours of f and g (see Figure 9.6.1). Unfortunately, the
	 * functions f and g have, in general, no relation to each other at all!
	 * There is nothing special about a common point from either f’s point of
	 * view, or from g’s. In this section we will discuss the simplest
	 * multidimensional root finding method, Newton-Raphson. This method gives
	 * you a very efficient means of converging to a root, if you have a
	 * sufficiently good initial guess. It can also spectacularly fail to
	 * converge, indicating (though not proving) that your putative root does
	 * not exist nearby. In §9.7 we discuss more sophisticated implementations
	 * of the Newton-Raphson method, which try to improve on Newton-Raphson’s
	 * poor global convergence. A multidimensional generalization of the secant
	 * method, called Broyden’s method, is also discussed in §9.7.
	 * 
	 * A typical problem givesN functional relations to be zeroed, involving
	 * variables xi, i = 1, 2, . . .,N: Fi(x1, x2, . . . , xN) = 0 i = 1, 2, . .
	 * .,N. (9.6.2) We let x denote the entire vector of values xi and F denote
	 * the entire vector of functions Fi.
	 * 
	 * The following routine mnewt performs ntrial iterations starting from an
	 * initial guess at the solution vector x[1..n]. Iteration stops if either
	 * the sum of the magnitudes of the functions Fi is less than some tolerance
	 * tolf, or the sum of the absolute values of the corrections to dxi is less
	 * than some tolerance tolx. mnewt calls a user supplied function usrfun
	 * which must provide the function values F and the Jacobian matrix J. If J
	 * is difficult to compute analytically, you can try having usrfun call the
	 * routine fdjac of §9.7 to compute the partial derivatives by finite
	 * differences. You should not make ntrial too big; rather inspect to see
	 * what is happening before continuing for some further iterations.
	 * 
	 * In this section we will discuss the simplest multidimensional root
	 * finding method, Newton-Raphson. This method gives you a very efficient
	 * means of converging to a root, if you have a sufficiently good initial
	 * guess. It can also spectacularly fail to converge, indicating (though not
	 * proving) that your putative root does not exist nearby. In §9.7 we
	 * discuss more sophisticated implementations of the Newton-Raphson method,
	 * which try to improve on Newton-Raphson’s poor global convergence. A
	 * multidimensional generalization of the secant method, called Broyden’s
	 * method, is also discussed in §9.7.
	 */
	// public void fdjac(int n, float x[], float fvec[], float **df,
	// void (*vecfunc)(int, float [], float []))
	/**
	 * Computes forward-difference approximation to Jacobian. On input, x[1..n] is the point at 
	 * which the Jacobian is to be evaluated, fvec[1..n] is the vector of function values at the 
	 * point, and vecfunc(n,x,f) is a user-supplied routine that returns the vector of functions at 
	 * x. On output, df[1..n][1..n] is the Jacobian array.
	 * @param n n
	 * @param x x
	 * @param fvec fvec
	 * @return the result
	 */
	public double[][] fdjac(int n, double[] x, double[] fvec)// , double[][] df)
	// Computes forward-difference approximation to Jacobian. On input, x[1..n]
	// is the point at
	// which the Jacobian is to be evaluated, fvec[1..n] is the vector of
	// function values at the
	// point, and vecfunc(n,x,f) is a user-supplied routine that returns the
	// vector of functions at
	// x. On output, df[1..n][1..n] is the Jacobian array.
	{
		double[][] df = new double[n][n];

		int i = 0;
		int j = 0;
		double h = 0.0;
		double temp = 0.0;
		double[] f = new double[n];
		// f=vector(1,n);
		for (j = 1; j <= n; j++) {
			temp = x[j - 1];// temp=x[j];
			h = EPS2 * Math.abs(temp);
			if (h == 0.0)
				h = EPS2;
			x[j - 1] = temp + h;// x[j]=temp+h; Trick to reduce finite precision
								// error.
			h = x[j - 1] - temp;// h=x[j]-temp;
			f = func.vecfunc(n, x);// (*vecfunc)(n,x,f);
			x[j - 1] = temp;// x[j]=temp;
			for (i = 1; i <= n; i++)
				// df[i][j]=(f[i]-fvec[i])/h; Forward difference formula.
				df[i - 1][j - 1] = (f[i - 1] - fvec[i - 1]) / h;
		}
		// free_vector(f,1,n);
		return df;
	}

	/**
	 * Returns f = 1/2 F · F at x. The vecfunc is a routine that returns the 
	 * vector of functions at x. It is set to point to a user-supplied routine in the calling program. 
	 * Global variables also communicate the function values back to the calling program.
	 * @param x x
	 * @return the result
	 */
	public double fmin(double x[])
	// Returns f = 1/2 F · F at x. The global pointer *nrfuncv points to a
	// routine that returns the
	// vector of functions at x. It is set to point to a user-supplied routine
	// in the calling program.
	// Global variables also communicate the function values back to the calling
	// program.
	{
		int i = 0;
		double sum = 0.0;
		// (*nrfuncv)(nn,x,fvec);
		fvec_newt = func.vecfunc(nn_newt, x);
		for (sum = 0.0, i = 1; i <= nn_newt; i++)
			// sum += (fvec_newt[i])*(fvec_newt[i]);
			sum += (fvec_newt[i - 1]) * (fvec_newt[i - 1]);
		return 0.5 * sum;
	}

	// void usrfun(float *x,int n,float *fvec,float **fjac);

	/**
	 * Given an initial guess x[1..n] for a root in n dimensions, take ntrial Newton-Raphson steps 
	 * to improve the root. Stop if the root converges in either summed absolute variable increments 
	 * tolx or summed absolute function values tolf.
	 * @param ntrial ntrial
	 * @param x x
	 * @param n n
	 * @param tolx tolx
	 * @param tolf tolf
	 */
	public void mnewt(int ntrial, double[] x, int n, double tolx, double tolf)
	// Given an initial guess x[1..n] for a root in n dimensions, take ntrial
	// Newton-Raphson steps
	// to improve the root. Stop if the root converges in either summed absolute
	// variable increments
	// tolx or summed absolute function values tolf.
	{
		// void lubksb(float **a, int n, int *indx, float b[]);
		// void ludcmp(float **a, int n, int *indx, float *d);
		int k = 0;
		int i = 0;
		int[] indx = new int[n];
		double errx = 0.0;
		double errf = 0.0;
		// double d = 0.0;
		double[] fvec = new double[n];
		double[][] fjac = new double[n][n];
		double[] p = new double[n];
		// indx=ivector(1,n);
		// p=vector(1,n);
		// fvec=vector(1,n);
		// fjac=matrix(1,n,1,n);
		for (k = 1; k <= ntrial; k++) {
			// /usrfun(x,n,fvec,fjac); User function supplies function values at
			// x in
			// fvec and Jacobian matrix in fjac.
			fvec = func.vecfunc(n, x);
			fjac = fdjac(n, x, fvec);

			errf = 0.0;
			for (i = 1; i <= n; i++)
				errf += Math.abs(fvec[i - 1]);// errf += fabs(fvec[i]); Check
												// function convergence.
			if (errf <= tolf)
				return;// FREERETURN
			for (i = 1; i <= n; i++)
				p[i - 1] = -fvec[i - 1];// p[i] = -fvec[i]; Right-hand side of
										// linear equations.
			// ludcmp(fjac,n,indx,&d); //Solve linear equations using LU
			// decomposition.
			LinAEq.ludcmp(fjac, n, indx);
			LinAEq.lubksb(fjac, n, indx, p);
			errx = 0.0;// Check root convergence.
			for (i = 1; i <= n; i++) {// Update solution.
				errx += Math.abs(p[i - 1]);// (p[i]);
				x[i - 1] += p[i - 1];// x[i] += p[i];
			}
			if (errx <= tolx)
				return;// FREERETURN
		}
		// FREERETURN
	}

	/**
	 * Used internally.
	 * @param x x
	 * @param s s
	 * @return the result
	 */
	public double localFunc(double[] x, String s) {
		double res = 0.0;
		if (s.compareTo("fmin") == 0)
			res = fmin(x);
		return res;
	}

	// Globally Convergent Methods for Nonlinear Systems of Equations

	// public void lnsrch(int n, float xold[], float fold, float g[], float p[],
	// float x[],
	// float *f, float stpmax, int *check, float (*func)(float []))

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
	 * @param namS namS pointing to fmin via localFunc
	 */
	public void lnsrch(int n, double[] xold, double fold, double[] g,
			double[] p, double[] x, double stpmax, String namS)// , int *check,
																// float
																// (*func)(float
																// []))
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
			f_ln = localFunc(x, namS);// func.vecfunc(n,x);//*f=(*func)(x);
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

	// Here MAXITS is the maximum number of iterations; TOLF sets the
	// convergence criterion on
	// function values; TOLMIN sets the criterion for deciding whether spurious
	// convergence to a
	// minimum of fmin has occurred; TOLX is the convergence criterion on dx;
	// STPMX is the scaled
	// maximum step length allowed in line searches.
	// int nn; Global variables to communicate with fmin.
	// float *fvec;
	// void (*nrfuncv)(int n, float v[], float f[]);
	// #define FREERETURN {free_vector(fvec,1,n);free_vector(xold,1,n);\
	// free_vector(p,1,n);free_vector(g,1,n);free_matrix(fjac,1,n,1,n);\
	// free_ivector(indx,1,n);return;}

	// public void newt(float x[], int n, int *check,void (*vecfunc)(int, float
	// [], float []))
	/**
	 * Given an initial guess x[1..n] for a root in n dimensions, find the root 
	 * by a globally convergent Newton’s method. The vector of functions to be zeroed, called fvec_newt
	 * in the routine below, is returned by the user-supplied routine. The 
	 * output quantity check is false (0) on a normal return and true (1) if the routine has 
	 * converged to a local minimum of the function fmin defined below. In this case try restarting 
	 * from a different initial guess.
	 * @param x x
	 * @param n n
	 */
	public void newt(double[] x, int n)
	
	// Given an initial guess x[1..n] for a root in n dimensions, find the root
	// by a globally convergent
	// Newton’s method. The vector of functions to be zeroed, called fvec[1..n]
	// in the routine
	// below, is returned by the user-supplied routine vecfunc(n,x,fvec). The
	// output quantity
	// check is false (0) on a normal return and true (1) if the routine has
	// converged to a local
	// minimum of the function fmin defined below. In this case try restarting
	// from a different initial
	// guess.
	{
		// void fdjac(int n, float x[], float fvec[], float **df,void
		// (*vecfunc)(int, float [], float []));
		// float fmin(float x[]);
		// void lnsrch(int n, float xold[], float fold, float g[], float p[],
		// float x[],float *f, float stpmax, int *check, float (*func)(float
		// []));
		// void lubksb(float **a, int n, int *indx, float b[]);
		// void ludcmp(float **a, int n, int *indx, float *d);

		int i = 0;
		int its = 0;
		int j = 0;
		int[] indx = new int[n];
		// double d = 0.0;
		double den = 0.0;
		double f = 0.0;
		double fold = 0.0;
		double stpmax = 0.0;
		double sum = 0.0;
		double temp = 0.0;
		double test = 0.0;
		double[][] fjac = new double[n][n];
		double[] g = new double[n];
		double[] p = new double[n];
		double[] xold = new double[n];

		failB = false;
		// indx=ivector(1,n);
		// fjac=matrix(1,n,1,n);
		// g=vector(1,n);
		// p=vector(1,n);
		// xold=vector(1,n);
		// fvec=vector(1,n); Define global variables.

		nn_newt = n;
		// nrfuncv=vecfunc;
		f = fmin(x); // fvec is also computed by this call.
		test = 0.0; // Test for initial guess being a root. Use more stringent
					// test than simply TOLF.
		for (i = 1; i <= n; i++)
			// if (Math.abs(fvec[i]) > test) test=fabs(fvec[i]);
			if (Math.abs(fvec_newt[i - 1]) > test)
				test = Math.abs(fvec_newt[i - 1]);
		if (test < 0.01 * TOLF) {
			check_ln = 0;
			return;// FREERETURN
		}
		for (sum = 0.0, i = 1; i <= n; i++)
			// sum += SQR(x[i]); Calculate stpmax for line searches.
			sum += x[i - 1] * x[i - 1];
		stpmax = STPMX * Math.max(Math.sqrt(sum), (double) n);
		for (its = 1; its <= MAXITS; its++) {// Start of iteration loop.
			fjac = fdjac(n, x, fvec_newt);// fdjac(n,x,fvec,fjac,vecfunc);
			// If analytic Jacobian is available, you can replace the routine
			// fdjac below with your
			// own routine.
			for (i = 1; i <= n; i++) {// Compute gradientf for the line search.
				for (sum = 0.0, j = 1; j <= n; j++)
					sum += fjac[j - 1][i - 1] * fvec_newt[j - 1];// sum +=
																	// fjac[j][i]*fvec[j];
				g[i - 1] = sum;// g[i]=sum;
			}
			for (i = 1; i <= n; i++)
				xold[i - 1] = x[i - 1];// xold[i]=x[i]; Store x,
			fold = f;// and f.
			for (i = 1; i <= n; i++)
				p[i - 1] = -fvec_newt[i - 1];// p[i] = -fvec[i]; Right-hand side
												// for linear equations.
			LinAEq.ludcmp(fjac, n, indx);// ,&d);// Solve linear equations by LU
											// decomposition.
			LinAEq.lubksb(fjac, n, indx, p);
			// lnsrch(n,xold,fold,g,p,x,&f,stpmax,check,fmin);
			lnsrch(n, xold, fold, g, p, x, stpmax, "fmin");
			f = f_ln;
			// lnsrch returns new x and f. It also calculates fvec at the new x
			// when it calls fmin.
			test = 0.0; // Test for convergence on function values.
			for (i = 1; i <= n; i++)
				// if (fabs(fvec[i]) > test) test=fabs(fvec[i]);
				if (Math.abs(fvec_newt[i - 1]) > test)
					test = Math.abs(fvec_newt[i - 1]);
			if (test < TOLF) {
				check_ln = 0;
				return;// FREERETURN
			}
			if (check_ln != 0)// (*check)
			{// Check for gradient of f zero, i.e., spurious convergence.
				test = 0.0;
				den = Math.max(f, 0.5 * n);
				for (i = 1; i <= n; i++) {
					// temp=fabs(g[i])*FMAX(fabs(x[i]),1.0)/den;
					temp = Math.abs(g[i - 1])
							* Math.max(Math.abs(x[i - 1]), 1.0) / den;
					if (temp > test)
						test = temp;
				}
				check_ln = (test < TOLMIN ? 1 : 0);
				return;// FREERETURN
			}
			test = 0.0;// Test for convergence on dx.
			for (i = 1; i <= n; i++) {
				// temp=(fabs(x[i]-xold[i]))/FMAX(fabs(x[i]),1.0);
				temp = (Math.abs(x[i - 1] - xold[i - 1]))
						/ Math.max(Math.abs(x[i - 1]), 1.0);
				if (temp > test)
					test = temp;
			}
			if (test < TOLX)
				return;// FREERETURN
		}
		// nrerror("MAXITS exceeded in newt");
		failB = true;
		failS = "MAXITS exceeded in newt";
	}

	// Newton’s method as implemented above is quite powerful, but it still has
	// several
	// disadvantages. One drawback is that the Jacobian matrix is needed. In
	// many problems
	// analytic derivatives are unavailable. If function evaluation is
	// expensive, then the cost of
	// finite-difference determination of the Jacobian can be prohibitive.
	// Just as the quasi-Newton methods to be discussed in §10.7 provide cheap
	// approximations
	// for the Hessian matrix in minimization algorithms, there are quasi-Newton
	// methods that
	// provide cheap approximations to the Jacobian for zero finding. These
	// methods are often called
	// secant methods, since they reduce to the secant method (§9.2) in one
	// dimension (see, e.g., [1]).
	// The best of these methods still seems to be the first one introduced,
	// Broyden’s method [2].

	// Here MAXITS is the maximum number of iterations; EPS is a number close to
	// the machine
	// precision; TOLF is the convergence criterion on function values; TOLX is
	// the convergence criterion
	// on dx; STPMX is the scaled maximum step length allowed in line searches;
	// TOLMIN is used to
	// decide whether spurious convergence to a minimum of fmin has occurred.
	// #define FREERETURN {free_vector(fvec,1,n);free_vector(xold,1,n);\
	// free_vector(w,1,n);free_vector(t,1,n);free_vector(s,1,n);\
	// free_matrix(r,1,n,1,n);free_matrix(qt,1,n,1,n);free_vector(p,1,n);\
	// free_vector(g,1,n);free_vector(fvcold,1,n);free_vector(d,1,n);\
	// free_vector(c,1,n);return;}
	// int nn; Global variables to communicate with fmin.
	// float *fvec;
	// void (*nrfuncv)(int n, float v[], float f[]);

	// void broydn(float x[], int n, int *check,
	// void (*vecfunc)(int, float [], float []))
	/**
	 * Given an initial guess x[1..n] for a root in n dimensions, find the root by Broyden’s method 
	 * embedded in a globally convergent strategy. The vector of functions to be zeroed, called 
	 * fvec_newt in the routine below, is returned by the user-supplied routine. 
	 * The routine fdjac and the function fmin from newt are used. The output quantity check 
	 * is false (0) on a normal return and true (1) if the routine has converged to a local minimum 
	 * of the function fmin or if Broyden’s method can make no further progress. In this case try 
	 * restarting from a different initial guess.
	 * @param x x
	 * @param n n
	 */
	public void broydn(double[] x, int n)
	// /Given an initial guess x[1..n] for a root in n dimensions, find the root
	// by Broyden’s method
	// embedded in a globally convergent strategy. The vector of functions to be
	// zeroed, called
	// fvec[1..n] in the routine below, is returned by the user-supplied routine
	// vecfunc(n,x,fvec).
	// The routine fdjac and the function fmin from newt are used. The output
	// quantity check
	// is false (0) on a normal return and true (1) if the routine has converged
	// to a local minimum
	// of the function fmin or if Broyden’s method can make no further progress.
	// In this case try
	// restarting from a different initial guess.
	{
		// void fdjac(int n, float x[], float fvec[], float **df,void
		// (*vecfunc)(int, float [], float []));
		// float fmin(float x[]);
		// void lnsrch(int n, float xold[], float fold, float g[], float p[],
		// float x[],float *f, float stpmax, int *check, float (*func)(float
		// []));
		// void qrdcmp(float **a, int n, float *c, float *d, int *sing);
		// void qrupdt(float **r, float **qt, int n, float u[], float v[]);
		// void rsolv(float **a, int n, float d[], float b[]);

		int i = 0;
		int its = 0;
		int j = 0;
		int k = 0;
		int restrt = 0;
		// int sing = 0;
		int skip = 0;
		double den = 0.0;
		double f = 0.0;
		double fold = 0.0;
		double stpmax = 0.0;
		double sum = 0.0;
		double temp = 0.0;
		double test = 0.0;
		double[] c = new double[n];
		double[] d = new double[n];
		double[] fvcold = new double[n];
		double[] g = new double[n];
		double[] p = new double[n];
		double[][] qt = new double[n][n];
		double[][] r = new double[n][n];
		double[] s = new double[n];
		double[] t = new double[n];
		double[] w = new double[n];
		double[] xold = new double[n];
		// c=vector(1,n);
		// d=vector(1,n);
		// fvcold=vector(1,n);
		// g=vector(1,n);
		// p=vector(1,n);
		// qt=matrix(1,n,1,n);
		// /r=matrix(1,n,1,n);
		// s=vector(1,n);
		// t=vector(1,n);
		// w=vector(1,n);
		// xold=vector(1,n);
		// fvec=vector(1,n); Define global variables.
		nn_newt = n;
		failB = false;
		// nrfuncv=vecfunc;
		f = fmin(x); // The vector fvec is also computed by this call.
		test = 0.0;
		for (i = 1; i <= n; i++)
			// Test for initial guess being a root. Use more
			// stringent test than simply TOLF.
			// if (fabs(fvec[i]) > test)test=fabs(fvec[i]);
			if (Math.abs(fvec_newt[i - 1]) > test)
				test = Math.abs(fvec_newt[i - 1]);
		if (test < 0.01 * TOLF) {
			check_ln = 0;
			return;// FREERETURN
		}
		for (sum = 0.0, i = 1; i <= n; i++)
			// sum += SQR(x[i]); //Calculate stpmax for line searches.
			sum += x[i - 1] * x[i - 1];
		stpmax = STPMX * Math.max(Math.sqrt(sum), (double) n);
		restrt = 1;// Ensure initial Jacobian gets computed.
		for (its = 1; its <= MAXITS; its++) {// Start of iteration loop.
			if (restrt != 0) // if (restrt)
			{
				// fdjac(n,x,fvec,r,vecfunc); //Initialize or reinitialize
				// Jacobian in r.
				r = fdjac(n, x, fvec_newt);
				LinAEq.qrdcmp(r, n, c, d);// qrdcmp(r,n,c,d,&sing); //QR
											// decomposition of Jacobian.
				if (LinAEq.sing != 0) // if (sing)
				{
					failB = true;
					failS = "singular Jacobian in broydn";
					return;
					// nrerror("singular Jacobian in broydn");
				}
				for (i = 1; i <= n; i++) { // Form QT explicitly.
					for (j = 1; j <= n; j++)
						qt[i - 1][j - 1] = 0.0;// qt[i][j]=0.0;
					qt[i - 1][i - 1] = 1.0;// qt[i][i]=1.0;
				}
				for (k = 1; k < n; k++) {
					if (c[k - 1] != 0.0) // if (c[k])
					{
						for (j = 1; j <= n; j++) {
							sum = 0.0;
							for (i = k; i <= n; i++)
								// sum += r[i][k]*qt[i][j];
								sum += r[i - 1][k - 1] * qt[i - 1][j - 1];
							sum /= c[k - 1];// sum /= c[k];
							for (i = k; i <= n; i++)
								qt[i - 1][j - 1] -= sum * r[i - 1][k - 1];// qt[i][j]
																			// -=
																			// sum*r[i][k];
						}
					}
				}
				for (i = 1; i <= n; i++) {// Form R explicitly.
					r[i - 1][i - 1] = d[i - 1];// r[i][i]=d[i];
					for (j = 1; j < i; j++)
						r[i - 1][j - 1] = 0.0;// r[i][j]=0.0;
				}
			} else {// Carry out Broyden update.
				for (i = 1; i <= n; i++)
					s[i - 1] = x[i - 1] - xold[i - 1]; // s =
														// dx.//s[i]=x[i]-xold[i];
														// s = dx.
				for (i = 1; i <= n; i++) {// t = R · s.
					for (sum = 0.0, j = i; j <= n; j++)
						sum += r[i - 1][j - 1] * s[j - 1];// sum +=
															// r[i][j]*s[j];
					t[i - 1] = sum;// t[i]=sum;
				}
				skip = 1;
				for (i = 1; i <= n; i++) {// w = dF - B · s.
					for (sum = 0.0, j = 1; j <= n; j++)
						sum += qt[j - 1][i - 1] * t[j - 1];// sum +=
															// qt[j][i]*t[j];
					w[i - 1] = fvec_newt[i - 1] - fvcold[i - 1] - sum;// w[i]=fvec[i]-fvcold[i]-sum;
					// if (fabs(w[i]) >= EPS*(fabs(fvec[i])+fabs(fvcold[i])))
					// skip=0;
					if (Math.abs(w[i - 1]) >= EPSS
							* (Math.abs(fvec_newt[i - 1]) + Math
									.abs(fvcold[i - 1])))
						skip = 0;
					// Don’t update with noisy components of w.
					else
						w[i - 1] = 0.0;// w[i]=0.0;
				}
				if (skip == 0) // if (!skip)
				{
					for (i = 1; i <= n; i++) {// t = QT · w.
						for (sum = 0.0, j = 1; j <= n; j++)
							sum += qt[i - 1][j - 1] * w[j - 1];// sum +=
																// qt[i][j]*w[j];
						t[i - 1] = sum; // t[i]=sum;
					}
					for (den = 0.0, i = 1; i <= n; i++)
						den += s[i - 1] * s[i - 1];// den += SQR(s[i]);
					for (i = 1; i <= n; i++)
						s[i - 1] /= den;// s[i] /= den; Store s/(s · s) in s.
					LinAEq.qrupdt(r, qt, n, t, s); // Update R and QT .
					for (i = 1; i <= n; i++) {
						// if (r[i][i] == 0.0) nrerror("r singular in broydn");
						if (r[i - 1][i - 1] == 0.0) {
							// nrerror("r singular in broydn");
							failB = true;
							failS = "r singular in broydn";
							return;
						}
						d[i - 1] = r[i - 1][i - 1];// d[i]=r[i][i]; Diagonal of
													// R stored in d.
					}
				}
			}
			for (i = 1; i <= n; i++) {// Right-hand side for linear equations is
										// -QT · F.
				for (sum = 0.0, j = 1; j <= n; j++)
					sum += qt[i - 1][j - 1] * fvec_newt[j - 1];// sum +=
																// qt[i][j]*fvec[j];
				p[i - 1] = -sum;// p[i] = -sum;
			}
			for (i = n; i >= 1; i--) {// Compute ?f ? (Q·R)T · F for the line
										// search.
				for (sum = 0.0, j = 1; j <= i; j++)
					sum -= r[j - 1][i - 1] * p[j - 1];// sum -= r[j][i]*p[j];
				g[i - 1] = sum;// g[i]=sum;
			}
			for (i = 1; i <= n; i++) {// Store x and F.
				xold[i - 1] = x[i - 1];// xold[i]=x[i];
				fvcold[i - 1] = fvec_newt[i - 1];// fvcold[i]=fvec[i];
			}
			fold = f;// Store f.
			LinAEq.rsolv(r, n, d, p); // Solve linear equations.
			// lnsrch(n,xold,fold,g,p,x,&f,stpmax,check,fmin);
			lnsrch(n, xold, fold, g, p, x, stpmax, "fmin");
			f = f_ln;
			// lnsrch returns new x and f. It also calculates fvec at the new x
			// when it calls fmin.
			test = 0.0; // Test for convergence on function values.
			for (i = 1; i <= n; i++)
				// if (Math.abs(fvec[i]) > test) test=fabs(fvec[i]);
				if (Math.abs(fvec_newt[i - 1]) > test)
					test = Math.abs(fvec_newt[i - 1]);
			if (test < TOLF) {
				check_ln = 0;
				return;// FREERETURN
			}
			if (check_ln != 0)// if (*check)
			{// True if line search failed to find a new x.
				if (restrt != 0)
					return;// if (restrt) //FREERETURN Failure; already tried
							// reinitializing the Jacobian.
				else {
					test = 0.0;// Check for gradient of f zero, i.e., spurious
								// convergence.
					den = Math.max(f, 0.5 * n);
					for (i = 1; i <= n; i++) {
						// temp=fabs(g[i])*FMAX(fabs(x[i]),1.0)/den;
						temp = Math.abs(g[i - 1])
								* Math.max(Math.abs(x[i - 1]), 1.0) / den;
						if (temp > test)
							test = temp;
					}
					if (test < TOLMIN)
						return;// FREERETURN
					else
						restrt = 1; // Try reinitializing the Jacobian.
				}
			} else {// Successful step; will use Broyden update for next step.
				restrt = 0;
				test = 0.0; // Test for convergence on dx.
				for (i = 1; i <= n; i++) {
					// temp=(fabs(x[i]-xold[i]))/FMAX(fabs(x[i]),1.0);
					temp = (Math.abs(x[i - 1] - xold[i - 1]))
							/ Math.max(Math.abs(x[i - 1]), 1.0);
					if (temp > test)
						test = temp;
				}
				if (test < TOLX)
					return;// FREERETURN
			}
		}
		failB = true;
		failS = "MAXITS exceeded in broydn";
		return;

		// nrerror("MAXITS exceeded in broydn");
		// FREERETURN
	}
//=======================
	private static boolean zero = false;
	//private static final double eps = 1e-10;

	private NumberFormat nf = NumberFormat.getInstance(Locale.US);
    //private String pattern="0.###E0";
    //private DecimalFormatSymbols dfs=new DecimalFormatSymbols(Locale.US);
	//private DecimalFormat nff =  new DecimalFormat(pattern,dfs);
	private int idigits=3;
	
	/**
	 * 3 digits number format.
	 */
	public void set_defaults()
	{
	    idigits=3;
		nf.setMinimumFractionDigits(idigits);//3);//default is 2!!
	    nf.setMaximumFractionDigits(idigits);//(3);//default is 2!!
	    nf.setGroupingUsed(false);//no 4,568.02 but 4568.02
	}

	/**
	 * Check if the root finding method succeeds.
	 * @return the result
	 */
	public boolean validate()
	{
		return zero;
	}

	/**
	 * Another implementation for secant method. Find the root of user supplied function in 
	 * [a,b] interval.
	 * @param a a
	 * @param b b
	 * @return the result
	 */
    public double Secant(double a, double b)
	{
        zero=true;
        double x=a;
        double dx=1;
        double fx =1;
        //---particularizare
        double fa=func.F(a);//functia de prelucrat
        //------------------
        if (Math.abs(fa)<=eps)
           return x;
        x=b;
        double fb=func.F(b);
        if (Math.abs(fb)<=eps)
           return x;

		if (fa*fb<0)
        {
		   while(zero)
		   {

              x=(a*fb-b*fa)/(fb-fa);
	          fx=func.F(x);
	          if (fa*fx>0)
	          {
	            a=x;
	            fa=fx;
		      }
	          else
	          {
	            b=x;
	            fb=fx;
		      }
	          dx=b-a;
	          if (x!=0)
                dx=dx/x;

              if(Math.abs(dx)<=eps || Math.abs(fx)<=eps)
                 break;
		   }
	    }
		else
	    {
	 	  zero=false;
	    }

		if(zero)
		{
			func.printSequence("==========advanced Secant evaluation=================");
			func.printSequence("Results= "+nf.format(x));
			func.printSequence("=====================================================");

	 	    return x;
	    }
	    else
	      return 0.0;

	}

    /**
     * Successive iterations for solution refinement.
     * @param x x, the gross solution or initial guess.
     * @return the result
     */
    public double Iter(double x)
	{
         int imax = 100;
         zero =true;
         int i=0;
         double test,dx = 1.0;

         while(zero)
         {
			 i++;
			 dx=func.F(x);//functia de prelucrat
			 test=1/dx;//de control al neconvergentei----+/-infinit
			 if (test==0.0)
		     {
				 //cand sirul e neconvergent
				 zero=false;
				 break;
			 }

			 x=x-dx;
			 if(x!=0)
			   dx=dx/x;

			 if(Math.abs(dx)<=eps || i>=imax)
                 break;
		 }

		 if (i>imax)
		    zero=false;


		  if (zero)
		  {

			func.printSequence("==========advanced Iter evaluation=================");
			func.printSequence("Results= "+nf.format(x));
			func.printSequence("===================================================");

		     return x;

 		  }
		  else
		     return 0.0;
	}

	/*//pe baza coeficientilor functiei de evaluat (polinom de interpolare) (a)
    //se returneza solutia in tabloul double[].
    //depinzand de coeficienti,si de ordinul polinomului de interpolare se construieste xx!!
	public static double[] birge(double[] a)
	{
		double x=0.0; double d=0.0; double p=0.0;
        zero=true;
		int nx=0; int i=0;
		int n = a.length-1;//grad polinom
		//daca am 3 coeficienti<-> polinomul de interpolare e de grad 2 si am 3-1 solutii!!!
        double [] xx = new double[n];
		int imax=100;
        if(n>1)
        {
           for (int m=n; m>=2; --m)
           {
			  i=0;
			  while(zero)
			  {
			      i++;
			      d=a[0];
			      p=a[0]*x+a[1];
			      for (int j=2; j<=m; j++)
			      {
				     d=d*x+p;
                     p=p*x+a[j];
				  }
				  if(d!=0)
				     d=p/d;
				  else
				     d=p;
				  x=x-d;//noua aproximare a radacinii
				  if(x!=0)
				     d=d/x;
				  if(Math.abs(d)<=eps || i>=imax)
                      break;
		      }
		      if (i>=imax)
		      {
				 xx[nx]=Double.NaN;
			     return trunc(xx);//mereu sunt gasite solutiile reale!!!
			  }
			  nx++;
			  xx[nx-1]=x;
			  for (int j=1; j<= m-1; j++)
                  a[j]=a[j-1]*x+a[j];
			}
		}
		nx =n;
		xx[n-1]=-a[1]/a[0];

		return xx;
	}

	private static double[] trunc(double[] xx)
	{
		 int j=-1;
	     for(int i=0 ; i<xx.length; i++)
	     {
            j++;
            //if (Convertor.doubleToString(xx[i]).compareTo("NaN")==0)
            if(Double.toString(xx[i]).compareTo("NaN")==0)
              break;
         }
         double [] x = new double[j];
         for (int i=0; i<j; i++)
            x[i]=xx[i];

         return x;
	}
	//pentru birge tre inversati coeficientii de la polinomul dat de interpolare
	//care este a0 + a1x + a2x^2 + a3x^3 + ... + anx^n cu tabloul a0,,...,an
	//la polinomul an,...a0!!!!*/
}
