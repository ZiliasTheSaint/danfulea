package danfulea.math.numerical;

import java.text.NumberFormat;
import java.util.Locale;

/**
 * Ordinary differential equations.
 * Based on literature, including Numerical Recipes in C (Cambridge Univ.).
 * @author Dan Fulea, 23 OCT. 2006
 */

public class OrdinaryDiffEq {
	public static Function func;
	public static boolean failB = false;
	public static String failS = "";

	public static double[][] y_rkdumb;
	public static double[] xx_rkdumb;

	public static double SAFETY = 0.9;
	public static double PGROW = -0.2;
	public static double PSHRNK = -0.25;
	public static double ERRCON = 1.89e-4;

	public static double MAXSTP = 10000;
	public static double TINY = 1.0e-30;
	public static int kmax_odeint = 0;
	public static int kount_odeint = 0;
	public static double[] xp_odeint;
	public static double[][] yp_odeint;
	public static double dxsav_odeint = 0.0;
	public static double x_rkqs = 0.0;
	public static double hdid_rkqs = 0.0;
	public static double hnext_rkqs = 0.0;

	public static int nok_odeint = 0;
	public static int nbad_odeint = 0;

	public static int KMAXX = 8;// Maximum row number used in the extrapolation.
	public static int IMAXX = (KMAXX + 1);
	public static double SAFE1 = 0.25;// Safety factors.
	public static double SAFE2 = 0.7;
	public static double REDMAX = 1.0e-5;// Maximum factor for stepsize
											// reduction.
	public static double REDMIN = 0.7;// Minimum factor for stepsize reduction.
	// public static double TINY =1.0e-30;// Prevents division by zero.
	public static double SCALMX = 0.1;// 1/SCALMX is the maximum factor by which
										// a stepsize can be increased.
	public static double[][] d_bsstep;// d=matrix(1,nv,1,KMAXX);
	public static double[] x_bsstep;// x=vector(1,KMAXX);
	// Pointers to matrix and vector used by pzextr or rzextr.
	public static double hdid_bsstep = 0.0;
	public static double hnext_bsstep = 0.0;
	public static double xx_bsstep = 0.0;

	public static int KMAXX1 = 7;// Maximum row number used in the
									// extrapolation.
	public static int IMAXX1 = (KMAXX1 + 1);
	/*
	 * #define KMAXX 7 #define IMAXX (KMAXX+1) #define SAFE1 0.25 #define SAFE2
	 * 0.7 #define REDMAX 1.0e-5 #define REDMIN 0.7 #define TINY 1.0e-30 #define
	 * SCALMX 0.1
	 */
	// float **d,*x;
	public static double[][] d_stifbs;// d=matrix(1,nv,1,KMAXX);
	public static double[] x_stifbs;// x=vector(1,KMAXX);
	public static double hdid_stifbs = 0.0;
	public static double hnext_stifbs = 0.0;
	public static double xx_stifbs = 0.0;
	// ------------------
	public static int first_stifbs = 1;
	public static int kmax_stifbs = 0;
	public static int kopt_stifbs = 0;
	public static int nvold_stifbs = -1;
	public static double epsold_stifbs = -1.0;
	public static double xnew_stifbs = 0.0;
	public static double[] a_stifbs = new double[IMAXX1 + 1];
	public static double[][] alf_stifbs = new double[KMAXX1 + 1][KMAXX1 + 1];
	public static int[] nseq_stifbs = { 0, 2, 6, 10, 14, 22, 34, 50, 70 }; // [IMAXX+1]=>1=>[9]
	// -----------
	public static int first_bsstep = 1;
	public static int kmax_bsstep = 0;
	public static int kopt_bsstep = 0;
	public static double epsold_bsstep = -1.0;
	public static double xnew_bsstep = 0.0;
	public static double[] a_bsstep = new double[IMAXX + 1];
	public static double[][] alf_bsstep = new double[KMAXX + 1][KMAXX + 1];
	// public static int nseq_bsstep[IMAXX+1]={0,2,4,6,8,10,12,14,16,18};
	public static int[] nseq_bsstep = { 0, 2, 4, 6, 8, 10, 12, 14, 16, 18 };

	public static double SAFETY_stiff = 0.9;
	public static double GROW_stiff = 1.5;
	public static double PGROW_stiff = -0.25;
	public static double SHRNK_stiff = 0.5;
	public static double PSHRNK_stiff = (-1.0 / 3.0);
	public static double ERRCON_stiff = 0.1296;
	public static int MAXTRY_stiff = 40;
	// Here NMAX is the maximum value of n; GROW and SHRNK are the largest and
	// smallest factors
	// by which stepsize can change in one step; ERRCON equals (GROW/SAFETY)
	// raised to the power
	// (1/PGROW) and handles the case when errmax  0.
	public static double GAM_stiff = (1.0 / 2.0);
	public static double A21_stiff = 2.0;
	public static double A31_stiff = (48.0 / 25.0);
	public static double A32_stiff = (6.0 / 25.0);
	public static double C21_stiff = -8.0;
	public static double C31_stiff = (372.0 / 25.0);
	public static double C32_stiff = (12.0 / 5.0);
	public static double C41_stiff = (-112.0 / 125.0);
	public static double C42_stiff = (-54.0 / 125.0);
	public static double C43_stiff = (-2.0 / 5.0);
	public static double B1_stiff = (19.0 / 9.0);
	public static double B2_stiff = (1.0 / 2.0);
	public static double B3_stiff = (25.0 / 108.0);
	public static double B4_stiff = (125.0 / 108.0);
	public static double E1_stiff = (17.0 / 54.0);
	public static double E2_stiff = (7.0 / 36.0);
	public static double E3_stiff = 0.0;
	public static double E4_stiff = (125.0 / 108.0);
	public static double C1X_stiff = (1.0 / 2.0);
	public static double C2X_stiff = (-3.0 / 2.0);
	public static double C3X_stiff = (121.0 / 50.0);
	public static double C4X_stiff = (29.0 / 250.0);
	public static double A2X_stiff = 1.0;
	public static double A3X_stiff = (3.0 / 5.0);
	// Here are the Kaps-Rentrop parameters, which can be substituted for those
	// of Shampine
	// simply by replacing the #define statements:

	public static double GAM_stiff1 = 0.231;
	public static double A21_stiff1 = 2.0;
	public static double A31_stiff1 = 4.52470820736;
	public static double A32_stiff1 = 4.16352878860;
	public static double C21_stiff1 = -5.07167533877;
	public static double C31_stiff1 = 6.02015272865;
	public static double C32_stiff1 = 0.159750684673;
	public static double C41_stiff1 = -1.856343618677;
	public static double C42_stiff1 = -8.50538085819;
	public static double C43_stiff1 = -2.08407513602;
	public static double B1_stiff1 = 3.95750374663;
	public static double B2_stiff1 = 4.62489238836;
	public static double B3_stiff1 = 0.617477263873;
	public static double B4_stiff1 = 1.282612945268;
	public static double E1_stiff1 = -2.30215540292;
	public static double E2_stiff1 = -3.07363448539;
	public static double E3_stiff1 = 0.873280801802;
	public static double E4_stiff1 = 1.282612945268;
	public static double C1X_stiff1 = GAM_stiff1;
	public static double C2X_stiff1 = -0.396296677520e-01;
	public static double C3X_stiff1 = 0.550778939579;
	public static double C4X_stiff1 = -0.553509845700e-01;
	public static double A2X_stiff1 = 0.462;
	public static double A3X_stiff1 = 0.880208333333;

	public static double hdid_stiff = 0.0;
	public static double hnext_stiff = 0.0;
	public static double x_stiff = 0.0;

	public static double SIGN(double a, double b) {
		return b >= 0.0 ? a : -a;
	}

	/*
	 * User storage for intermediate results. Preset kmax and dxsav in the
	 * calling program. If kmax != 0 results are stored at approximate intervals
	 * dxsav in the arrays xp[1..kount], yp[1..nvar] [1..kount], where kount is
	 * output by odeint. Defining declarations for these variables, with memory
	 * allocations xp[1..kmax] and yp[1..nvar][1..kmax] for the arrays, should
	 * be in the calling program.
	 */

	/*
	 * Problems involving ordinary differential equations (ODEs) can always be
	 * reduced to the study of sets of first-order differential equations. For
	 * example the second-order equation d2y/dx2 + q(x) dy/dx = r(x) (16.0.1)
	 * can be rewritten as two first-order equations dy/dx = z(x) dz/dx = r(x) -
	 * q(x)z(x) (16.0.2) where z is a new variable. This exemplifies the
	 * procedure for an arbitrary ODE. The usual choice for the new variables is
	 * to let them be just derivatives of each other (and of the original
	 * variable). Occasionally, it is useful to incorporate into their
	 * definition some other factors in the equation, or some powers of the
	 * independent variable, for the purpose of mitigating singular behavior
	 * that could result in overflows or increased roundoff error. Let common
	 * sense be your guide: If you find that the original variables are smooth
	 * in a solution, while your auxiliary variables are doing crazy things,
	 * then figure out why and choose different auxiliary variables.
	 * 
	 * The generic problem in ordinary differential equations is thus reduced to
	 * the study of a set of N coupled first-order differential equations for
	 * the functions yi, i = 1, 2, . . .,N, having the general form dyi(x)/dx =
	 * fi(x, y1, . . . , yN), i= 1, . . .,N (16.0.3) where the functions fi on
	 * the right-hand side are known. A problem involving ODEs is not completely
	 * specified by its equations. Even more crucial in determining how to
	 * attack the problem numerically is the nature of the problem’s boundary
	 * conditions. Boundary conditions are algebraic conditions on the values of
	 * the functions yi in (16.0.3). In general they can be satisfied at
	 * discrete specified points, but do not hold between those points, i.e.,
	 * are not preserved automatically by the differential equations. Boundary
	 * conditions can be as simple as requiring that certain variables have
	 * certain numerical values, or as complicated as a set of nonlinear
	 * algebraic equations among the variables. Usually, it is the nature of the
	 * boundary conditions that determines which numerical methods will be
	 * feasible. Boundary conditions divide into two broad categories. • In
	 * initial value problems all the yi are given at some starting value xs,
	 * and it is desired to find the yi’s at some final point xf , or at some
	 * discrete list of points (for example, at tabulated intervals). • In
	 * two-point boundary value problems, on the other hand, boundary conditions
	 * are specified at more than one x. Typically, some of the conditions will
	 * be specified at xs and the remainder at xf .
	 * 
	 * This chapter will consider exclusively the initial value problem,
	 * deferring twopoint boundary value problems, which are generally more
	 * difficult, to Chapter 17. The underlying idea of any routine for solving
	 * the initial value problem is always this: Rewrite the dy’s and dx’s in
	 * (16.0.3) as finite steps?y and?x, and multiply the equations by?x. This
	 * gives algebraic formulas for the change in the functions when the
	 * independent variable x is “stepped” by one “stepsize” ?x. In the limit of
	 * making the stepsize very small, a good approximation to the underlying
	 * differential equation is achieved. Literal implementation of this
	 * procedure results in Euler’s method (16.1.1, below), which is, however,
	 * not recommended for any practical use. Euler’s method is conceptually
	 * important, however; oneway or another, practical methods all come down to
	 * this same idea: Add small increments to your functions corresponding to
	 * derivatives (right-hand sides of the equations) multiplied by stepsizes.
	 * In this chapter we consider three major types of practical numerical
	 * methods for solving initial value problems for ODEs: • Runge-Kutta
	 * methods • Richardson extrapolation and its particular implementation as
	 * the Bulirsch- Stoer method • predictor-corrector methods. A brief
	 * description of each of these types follows. 1. Runge-Kutta methods
	 * propagate a solution over an interval by combining the information from
	 * several Euler-style steps (each involving one evaluation of the
	 * right-hand f’s), and then using the information obtained to match a
	 * Taylor series expansion up to some higher order. 2. Richardson
	 * extrapolation uses the powerful idea of extrapolating a computed result
	 * to the value that would have been obtained if the stepsize had been very
	 * much smaller than it actually was. In particular, extrapolation to zero
	 * stepsize is the desired goal. The first practical ODE integrator that
	 * implemented this idea was developed by Bulirsch and Stoer, and so
	 * extrapolation methods are often called Bulirsch-Stoer methods. 3.
	 * Predictor-corrector methods store the solution along the way, and use
	 * those results to extrapolate the solution one step advanced; they then
	 * correct the extrapolation using derivative information at the new point.
	 * These are best for very smooth functions. Runge-Kutta is what you use
	 * when (i) you don’t know any better, or (ii) you have an intransigent
	 * problem where Bulirsch-Stoer is failing, or (iii) you have a trivial
	 * problem where computational efficiency is of no concern. Runge-Kutta
	 * succeeds virtually always; but it is not usually fastest, except when
	 * evaluating fi is cheap and moderate accuracy (<~ 10-5) is required.
	 * Predictor-corrector methods, since they use past information, are
	 * somewhat more difficult to start up, but, for many smooth problems, they
	 * are computationally more efficient than Runge-Kutta. In recent years
	 * Bulirsch-Stoer has been replacing predictor-corrector in many
	 * applications, but it is too soon to say that predictor-corrector is
	 * dominated in all cases. However, it appears that only rather
	 * sophisticated predictor-corrector routines are competitive. Accordingly,
	 * we have chosen not to give an implementation of predictor-corrector in
	 * this book. We discuss predictor-corrector further in §16.7, so that you
	 * can use a canned routine should you encounter a suitable problem. In our
	 * experience, the relatively simple Runge-Kutta and Bulirsch-Stoer routines
	 * we give are adequate for most problems. Each of the three types of
	 * methods can be organized to monitor internal consistency. This allows
	 * numerical errors which are inevitably introduced into the solution to be
	 * controlled by automatic, (adaptive) changing of the fundamental stepsize.
	 * We always recommend that adaptive stepsize control be implemented, and we
	 * will do so below. In general, all three types of methods can be applied
	 * to any initial value problem. Each comes with its own set of debits and
	 * credits that must be understood before it is used.
	 * 
	 * Runge-Kutta Method Fourth-order Runge-Kutta method. You input the values
	 * of the independent variables, and you get out new values which are
	 * stepped by a stepsize h (which can be positive or negative). You will
	 * notice that the routine requires you to supply not only function derivs
	 * for calculating the right-hand side, but also values of the derivatives
	 * at the starting point. Why not let the routine call derivs for this first
	 * value? The answer will become clear only in the next section, but in
	 * brief is this: This call may not be your only one with these starting
	 * conditions. You may have taken a previous step with too large a stepsize,
	 * and this is your replacement. In that case, you do not want to call
	 * derivs unnecessarily at the start. Note that the routine that follows
	 * has, therefore, only three calls to derivs.
	 */

	/**
	 * Given values for the variables y[1..n] and their derivatives dydx[1..n] known at x, use the 
	 * fourth-order Runge-Kutta method to advance the solution over an interval h and return the 
	 * incremented variables as yout[1..n], which need not be a distinct array from y. The user 
	 * supplies the routine derivs(x,y,dydx), which returns derivatives dydx at x.
	 * @param y y
	 * @param dydx dydx
	 * @param n n
	 * @param x x
	 * @param h h
	 * @param yout yout
	 */
	public static void rk4(double[] y, double[] dydx, int n, double x,
			double h, double[] yout)// ,
	// void (*derivs)(float, float [], float []))
	// Given values for the variables y[1..n] and their derivatives dydx[1..n]
	// known at x, use the
	// fourth-order Runge-Kutta method to advance the solution over an interval
	// h and return the
	// incremented variables as yout[1..n], which need not be a distinct array
	// from y. The user
	// supplies the routine derivs(x,y,dydx), which returns derivatives dydx at
	// x.
	{
		int i = 0;
		double xh = 0.0;
		double hh = 0.0;
		double h6 = 0.0;
		double[] dym = new double[n];
		double[] dyt = new double[n];
		double[] yt = new double[n];
		// dym=vector(1,n);
		// dyt=vector(1,n);
		// yt=vector(1,n);
		hh = h * 0.5;
		h6 = h / 6.0;
		xh = x + hh;
		for (i = 1; i <= n; i++)
			yt[i - 1] = y[i - 1] + hh * dydx[i - 1];// yt[i]=y[i]+hh*dydx[i];
													// First step.
		// (*derivs)(xh,yt,dyt); Second step.
		dyt = func.derivF(xh, yt);
		for (i = 1; i <= n; i++)
			yt[i - 1] = y[i - 1] + hh * dyt[i - 1];// yt[i]=y[i]+hh*dyt[i];
		// (*derivs)(xh,yt,dym); Third step.
		dym = func.derivF(xh, yt);
		for (i = 1; i <= n; i++) {
			yt[i - 1] = y[i - 1] + h * dym[i - 1];// yt[i]=y[i]+h*dym[i];
			dym[i - 1] += dyt[i - 1];// dym[i] += dyt[i];
		}
		// (*derivs)(x+h,yt,dyt); Fourth step.
		dyt = func.derivF(x + h, yt);
		for (i = 1; i <= n; i++)
			// Accumulate increments with proper weights.
			// yout[i]=y[i]+h6*(dydx[i]+dyt[i]+2.0*dym[i]);
			yout[i - 1] = y[i - 1] + h6
					* (dydx[i - 1] + dyt[i - 1] + 2.0 * dym[i - 1]);
		// free_vector(yt,1,n);
		// free_vector(dyt,1,n);
		// free_vector(dym,1,n);
	}

	/**
	 * Initialize rkdumb.
	 * @param nvar nvar
	 * @param nstep nstep
	 */
	public static void rkdumb_init(int nvar, int nstep) {
		y_rkdumb = new double[nvar][nstep + 1];
		xx_rkdumb = new double[nstep + 1];
	}

	// float **y,*xx; For communication back to main.
	/**
	 * Starting from initial values vstart[1..nvar] known at x1 use fourth-order Runge-Kutta 
	 * to advance nstep equal increments to x2. The user-supplied routine derivs(x,v,dvdx) 
	 * evaluates derivatives. Results are stored in the global variables y_rkdumb[1..nvar][1..nstep+1] 
	 * and xx_rkdumb[1..nstep+1].
	 * @param vstart vstar
	 * @param nvar nvar
	 * @param x1 x1
	 * @param x2 x2
	 * @param nstep nstep
	 */
	public static void rkdumb(double[] vstart, int nvar, double x1, double x2,
			int nstep)// ,
	// void (*derivs)(float, float [], float []))
	// Starting from initial values vstart[1..nvar] known at x1 use fourth-order
	// Runge-Kutta
	// to advance nstep equal increments to x2. The user-supplied routine
	// derivs(x,v,dvdx)
	// evaluates derivatives. Results are stored in the global variables
	// y[1..nvar][1..nstep+1]
	// and xx[1..nstep+1].
	{
		failB = false;
		// void rk4(float y[], float dydx[], int n, float x, float h, float
		// yout[],
		// void (*derivs)(float, float [], float []));
		int i = 0;
		int k = 0;
		double x = 0.0;
		double h = 0.0;

		double[] v = new double[nvar];
		double[] vout = new double[nvar];
		double[] dv = new double[nvar];

		// v=vector(1,nvar);
		// vout=vector(1,nvar);
		// dv=vector(1,nvar);
		for (i = 1; i <= nvar; i++) {// Load starting values.
			v[i - 1] = vstart[i - 1];// v[i]=vstart[i];
			y_rkdumb[i - 1][0] = v[i - 1];// y[i][1]=v[i];
		}
		xx_rkdumb[0] = x1;// xx[1]=x1;
		x = x1;
		h = (x2 - x1) / nstep;
		for (k = 1; k <= nstep; k++) {// Take nstep steps.
										// (*derivs)(x,v,dv);
			dv = func.derivF(x, v);
			rk4(v, dv, nvar, x, h, vout);// ,derivs);
			if ((double) (x + h) == x) {
				// nrerror("Step size too small in routine rkdumb");
				failB = true;
				failS = "Step size too small in routine rkdumb";
				return;
			}
			x += h;
			xx_rkdumb[k] = x;// xx[k+1]=x; Store intermediate steps.
			for (i = 1; i <= nvar; i++) {
				v[i - 1] = vout[i - 1];// v[i]=vout[i];
				y_rkdumb[i - 1][k] = v[i - 1];// y_rkdumb[i][k+1]=v[i];
			}
		}
		// free_vector(dv,1,nvar);
		// free_vector(vout,1,nvar);
		// free_vector(v,1,nvar);
	}

	/*
	 * Adaptive Stepsize Control for Runge-Kutta
	 * 
	 * AgoodODEintegrator should exert some adaptive control over its own
	 * progress, making frequent changes in its stepsize. Usually the purpose of
	 * this adaptive stepsize control is to achieve some predetermined accuracy
	 * in the solution with minimum computational effort. Many small steps
	 * should tiptoe through treacherous terrain, while a few great strides
	 * should speed through smooth uninteresting countryside. The resulting
	 * gains in efficiency are not mere tens of percents or factors of two; they
	 * can sometimes be factors of ten, a hundred, or more. Sometimes accuracy
	 * may be demanded not directly in the solution itself, but in some related
	 * conserved quantity that can be monitored. Implementation of adaptive
	 * stepsize control requires that the stepping algorithm signal information
	 * about its performance,most important, an estimate of its truncation
	 * error. In this section we will learn how such information can be
	 * obtained.
	 */

	// public static void rkqs(double[] y, double[] dydx, int n, float *x, float
	// htry, float eps,
	// float yscal[], float *hdid, float *hnext,
	// void (*derivs)(float, float [], float []))
	/**
	 * Fifth-order Runge-Kutta step with monitoring of local truncation error to ensure accuracy and 
	 * adjust stepsize. Input are the dependent variable vector y[1..n] and its derivative dydx[1..n] 
	 * at the starting value of the independent variable x. Also input are the stepsize to be attempted 
	 * htry, the required accuracy eps, and the vector yscal[1..n] against which the error is 
	 * scaled. On output, y and x are replaced by their new values, hdid is the stepsize that was 
	 * actually accomplished, and hnext is the estimated next stepsize. derivs is the user-supplied 
	 * routine that computes the right-hand side derivatives.
	 * @param y y
	 * @param dydx dydx
	 * @param n n
	 * @param x x
	 * @param htry htry
	 * @param eps eps
	 * @param yscal yscal
	 */
	public static void rkqs(double[] y, double[] dydx, int n, double x,
			double htry, double eps, double[] yscal)// , float *hdid, float
													// *hnext,
	// void (*derivs)(float, float [], float []))
	// Fifth-order Runge-Kutta step with monitoring of local truncation error to
	// ensure accuracya nd
	// adjust stepsize. Input are the dependent variable vector y[1..n] and its
	// derivative dydx[1..n]
	// at the starting value of the independent variable x. Also input are the
	// stepsize to be attempted
	// htry, the required accuracy eps, and the vector yscal[1..n] against which
	// the error is
	// scaled. On output, y and x are replaced bythei r new values, hdid is the
	// stepsize that was
	// actuallyac complished, and hnext is the estimated next stepsize. derivs
	// is the user-supplied
	// routine that computes the right-hand side derivatives.
	{
		// void rkck(float y[], float dydx[], int n, float x, float h,
		// float yout[], float yerr[], void (*derivs)(float, float [], float
		// []));
		x_rkqs = x;
		failB = false;

		int i = 0;
		double errmax = 0.0;
		double h = 0.0;
		double htemp = 0.0;
		double xnew = 0.0;
		double[] yerr = new double[n];
		double[] ytemp = new double[n];
		// yerr=vector(1,n);
		// ytemp=vector(1,n);
		h = htry; // Set stepsize to the initial trial value.
		for (;;) {
			rkck(y, dydx, n, x_rkqs, h, ytemp, yerr);// ,derivs); Take a step.
			errmax = 0.0; // Evaluate accuracy.
			for (i = 1; i <= n; i++)
				// errmax=Math.max(errmax,Math.abs(yerr[i]/yscal[i]));
				errmax = Math.max(errmax, Math.abs(yerr[i - 1] / yscal[i - 1]));
			errmax /= eps;// Scale relative to required tolerance.
			if (errmax <= 1.0)
				break; // Step succeeded. Compute size of next step.
			htemp = SAFETY * h * Math.pow(errmax, PSHRNK);
			// Truncation error too large, reduce stepsize.
			h = (h >= 0.0 ? Math.max(htemp, 0.1 * h) : Math.min(htemp, 0.1 * h));
			// No more than a factor of 10.
			xnew = (x_rkqs) + h;
			if (xnew == x_rkqs) {
				failB = true;
				failS = "stepsize underflow in rkqs";
				return;
				// nrerror("stepsize underflow in rkqs");
			}
		}
		if (errmax > ERRCON)
			hnext_rkqs = SAFETY * h * Math.pow(errmax, PGROW);
		else
			hnext_rkqs = 5.0 * h; // No more than a factor of 5 increase.
		x_rkqs += (hdid_rkqs = h);
		for (i = 1; i <= n; i++)
			y[i - 1] = ytemp[i - 1];// y[i]=ytemp[i];
		// free_vector(ytemp,1,n);
		// free_vector(yerr,1,n);
	}

	/**
	 * Given values for n variables y[1..n] and their derivatives dydx[1..n] known at x, use 
	 * the fifth-order Cash-Karp Runge-Kutta method to advance the solution over an interval h 
	 * and return the incremented variables as yout[1..n]. Also return an estimate of the local 
	 * truncation error in yout using the embedded fourth-order method. The user supplies the routine 
	 * derivs(x,y,dydx), which returns derivatives dydx at x.
	 * @param y y
	 * @param dydx dydx
	 * @param n n
	 * @param x x
	 * @param h h
	 * @param yout yout
	 * @param yerr yerr
	 */
	public static void rkck(double[] y, double[] dydx, int n, double x,
			double h, double[] yout, double[] yerr)// , void (*derivs)(float,
													// float [], float []))
	// Given values for n variables y[1..n] and their derivatives dydx[1..n]
	// known at x, use
	// the fifth-order Cash-Karp Runge-Kutta method to advance the solution over
	// an interval h
	// and return the incremented variables as yout[1..n]. Also return an
	// estimate of the local
	// truncation error in yout using the embedded fourth-order method. The user
	// supplies the routine
	// derivs(x,y,dydx), which returns derivatives dydx at x.
	{
		int i = 0;
		/*
		 * static float a2=0.2,a3=0.3,a4=0.6,a5=1.0,a6=0.875,b21=0.2,
		 * b31=3.0/40.0,b32=9.0/40.0,b41=0.3,b42 = -0.9,b43=1.2, b51 =
		 * -11.0/54.0, b52=2.5,b53 = -70.0/27.0,b54=35.0/27.0,
		 * b61=1631.0/55296.0,b62=175.0/512.0,b63=575.0/13824.0,
		 * b64=44275.0/110592.0,b65=253.0/4096.0,c1=37.0/378.0,
		 * c3=250.0/621.0,c4=125.0/594.0,c6=512.0/1771.0, dc5 = -277.00/14336.0;
		 */
		double a2 = 0.2;
		double a3 = 0.3;
		double a4 = 0.6;
		double a5 = 1.0;
		double a6 = 0.875;
		double b21 = 0.2;
		double b31 = 3.0 / 40.0;
		double b32 = 9.0 / 40.0;
		double b41 = 0.3;
		double b42 = -0.9;
		double b43 = 1.2;
		double b51 = -11.0 / 54.0;
		double b52 = 2.5;
		double b53 = -70.0 / 27.0;
		double b54 = 35.0 / 27.0;
		double b61 = 1631.0 / 55296.0;
		double b62 = 175.0 / 512.0;
		double b63 = 575.0 / 13824.0;
		double b64 = 44275.0 / 110592.0;
		double b65 = 253.0 / 4096.0;
		double c1 = 37.0 / 378.0;
		double c3 = 250.0 / 621.0;
		double c4 = 125.0 / 594.0;
		double c6 = 512.0 / 1771.0;
		double dc5 = -277.00 / 14336.0;

		double dc1 = c1 - 2825.0 / 27648.0;
		double dc3 = c3 - 18575.0 / 48384.0;
		double dc4 = c4 - 13525.0 / 55296.0;
		double dc6 = c6 - 0.25;

		double[] ak2 = new double[n];
		double[] ak3 = new double[n];
		double[] ak4 = new double[n];
		double[] ak5 = new double[n];
		double[] ak6 = new double[n];
		double[] ytemp = new double[n];
		// ak2=vector(1,n);
		// ak3=vector(1,n);
		// ak4=vector(1,n);
		// ak5=vector(1,n);
		// ak6=vector(1,n);
		// ytemp=vector(1,n);
		for (i = 1; i <= n; i++)
			// First step.
			// ytemp[i]=y[i]+b21*h*dydx[i];
			ytemp[i - 1] = y[i - 1] + b21 * h * dydx[i - 1];
		// (*derivs)(x+a2*h,ytemp,ak2); Second step.
		ak2 = func.derivF(x + a2 * h, ytemp);
		for (i = 1; i <= n; i++)
			// ytemp[i]=y[i]+h*(b31*dydx[i]+b32*ak2[i]);
			ytemp[i - 1] = y[i - 1] + h
					* (b31 * dydx[i - 1] + b32 * ak2[i - 1]);
		// (*derivs)(x+a3*h,ytemp,ak3); Third step.
		ak3 = func.derivF(x + a3 * h, ytemp);
		for (i = 1; i <= n; i++)
			// ytemp[i]=y[i]+h*(b41*dydx[i]+b42*ak2[i]+b43*ak3[i]);
			ytemp[i - 1] = y[i - 1] + h
					* (b41 * dydx[i - 1] + b42 * ak2[i - 1] + b43 * ak3[i - 1]);
		// (*derivs)(x+a4*h,ytemp,ak4); Fourth step.
		ak4 = func.derivF(x + a4 * h, ytemp);
		for (i = 1; i <= n; i++)
			// ytemp[i]=y[i]+h*(b51*dydx[i]+b52*ak2[i]+b53*ak3[i]+b54*ak4[i]);
			ytemp[i - 1] = y[i - 1]
					+ h
					* (b51 * dydx[i - 1] + b52 * ak2[i - 1] + b53 * ak3[i - 1] + b54
							* ak4[i - 1]);
		// (*derivs)(x+a5*h,ytemp,ak5); Fifth step.
		ak5 = func.derivF(x + a5 * h, ytemp);
		for (i = 1; i <= n; i++)
			// ytemp[i]=y[i]+h*(b61*dydx[i]+b62*ak2[i]+b63*ak3[i]+b64*ak4[i]+b65*ak5[i]);
			ytemp[i - 1] = y[i - 1]
					+ h
					* (b61 * dydx[i - 1] + b62 * ak2[i - 1] + b63 * ak3[i - 1]
							+ b64 * ak4[i - 1] + b65 * ak5[i - 1]);
		// (*derivs)(x+a6*h,ytemp,ak6); Sixth step.
		ak6 = func.derivF(x + a6 * h, ytemp);
		for (i = 1; i <= n; i++)
			// Accumulate increments with proper weights.
			// yout[i]=y[i]+h*(c1*dydx[i]+c3*ak3[i]+c4*ak4[i]+c6*ak6[i]);
			yout[i - 1] = y[i - 1]
					+ h
					* (c1 * dydx[i - 1] + c3 * ak3[i - 1] + c4 * ak4[i - 1] + c6
							* ak6[i - 1]);
		for (i = 1; i <= n; i++)
			// yerr[i]=h*(dc1*dydx[i]+dc3*ak3[i]+dc4*ak4[i]+dc5*ak5[i]+dc6*ak6[i]);
			yerr[i - 1] = h
					* (dc1 * dydx[i - 1] + dc3 * ak3[i - 1] + dc4 * ak4[i - 1]
							+ dc5 * ak5[i - 1] + dc6 * ak6[i - 1]);
		// Estimate error as difference between fourth and fifth order methods.
		// free_vector(ytemp,1,n);
		// free_vector(ak6,1,n);
		// free_vector(ak5,1,n);
		// free_vector(ak4,1,n);
		// free_vector(ak3,1,n);
		// free_vector(ak2,1,n);
	}

	/**
	 * Runge-Kutta driver with adaptive stepsize control. Integrate starting values ystart[1..nvar] 
	 * from x1 to x2 with accuracy eps, storing intermediate results in global variables. h1 should 
	 * be set as a guessed first stepsize, hmin as the minimum allowed stepsize (can be zero). On 
	 * output nok and nbad are the number of good and bad (but retried and fixed) steps taken, 
	 * and ystart is replaced by values at the end of the integration interval. 
	 * derivs is the user-supplied routine for calculating the right-hand side derivative, while rkqs is the 
	 * name of the stepper routine to be used.
	 * @param ystart ystart
	 * @param nvar nvar
	 * @param x1 x1
	 * @param x2 x2
	 * @param eps eps
	 * @param h1 h1
	 * @param hmin hmin
	 */
	public static void odeint(double[] ystart, int nvar, double x1, double x2,
			double eps, double h1, double hmin)// , int *nok, int *nbad,
	// void (*derivs)(float, float [], float []),
	// void (*rkqs)(float [], float [], int, float *, float, float, float [],
	// float *, float *, void (*)(float, float [], float [])))
	// Runge-Kutta driver with adaptive stepsize control. Integrate starting
	// values ystart[1..nvar]
	// from x1 to x2 with accuracy eps, storing intermediate results in global
	// variables. h1 should
	// be set as a guessed first stepsize, hmin as the minimum allowed stepsize
	// (can be zero). On
	// output nok and nbad are the number of good and bad (but retried and
	// fixed) steps taken, and
	// ystart is replaced byv alues at the end of the integration interval.
	// derivs is the user-supplied
	// routine for calculating the right-hand side derivative, while rkqs is the
	// name of the stepper
	// routine to be used.
	{
		failB = false;

		int nstp = 0;
		int i = 0;
		double xsav = 0.0;
		double x = 0.0;
		double hnext = 0.0;// double hdid=0.0;
		double h = 0.0;
		double[] yscal = new double[nvar];
		double[] y = new double[nvar];
		double[] dydx = new double[nvar];
		// yscal=vector(1,nvar);
		// y=vector(1,nvar);
		// dydx=vector(1,nvar);
		x = x1;
		h = SIGN(h1, x2 - x1);
		nok_odeint = (nbad_odeint) = kount_odeint = 0;
		for (i = 1; i <= nvar; i++)
			y[i - 1] = ystart[i - 1];// y[i]=ystart[i];
		if (kmax_odeint > 0)
			xsav = x - dxsav_odeint * 2.0; // Assures storage of first step.
		for (nstp = 1; nstp <= MAXSTP; nstp++) {// Take at most MAXSTP steps.
												// (*derivs)(x,y,dydx);
			dydx = func.derivF(x, y);
			for (i = 1; i <= nvar; i++)
				// Scaling used to monitor accuracy. This general-purpose choice
				// can be modified if need be.
				// yscal[i]=fabs(y[i])+fabs(dydx[i]*h)+TINY;
				yscal[i - 1] = Math.abs(y[i - 1]) + Math.abs(dydx[i - 1] * h)
						+ TINY;
			if (kmax_odeint > 0 && kount_odeint < kmax_odeint - 1
					&& Math.abs(x - xsav) > Math.abs(dxsav_odeint)) {
				xp_odeint[++kount_odeint - 1] = x;// xp_odeint[++kount_odeint]=x;
													// //Store intermediate
													// results.
				for (i = 1; i <= nvar; i++)
					yp_odeint[i - 1][kount_odeint - 1] = y[i - 1];// yp_odeint[i][kount]=y[i];
				xsav = x;
			}
			if ((x + h - x2) * (x + h - x1) > 0.0)
				h = x2 - x; // If stepsize can overshoot, decrease.
			// (*rkqs)(y,dydx,nvar,&x,h,eps,yscal,&hdid,&hnext,derivs);
			rkqs(y, dydx, nvar, x, h, eps, yscal);// ,&hdid,&hnext,derivs);
			x = x_rkqs;// hdid=hdid_rkqs;
			hnext = hnext_rkqs;
			if (failB)
				return;// @@@@@@@@@@@@@@@@@
			if (hdid_rkqs == h)
				++(nok_odeint);
			else
				++(nbad_odeint);
			if ((x - x2) * (x2 - x1) >= 0.0) {// Are we done?
				for (i = 1; i <= nvar; i++)
					ystart[i - 1] = y[i - 1];// ystart[i]=y[i];
				if (kmax_odeint != 0) // if (kmax)
				{
					xp_odeint[++kount_odeint - 1] = x;// xp_odeint[++kount_odeint]=x;
														// //Save final step.
					for (i = 1; i <= nvar; i++)
						yp_odeint[i - 1][kount_odeint - 1] = y[i];// yp[i][kount]=y[i];
				}
				// free_vector(dydx,1,nvar);
				// free_vector(y,1,nvar);
				// free_vector(yscal,1,nvar);
				return;// Normal exit.
			}
			if (Math.abs(hnext) <= hmin) {
				// nrerror("Step size too small in odeint");
				failB = true;
				failS = "Step size too small in odeint";
				return;
			}
			h = hnext;
		}
		// nrerror("Too many steps in routine odeint");
		failB = true;
		failS = "Too many steps in routine odeint";
		return;
	}

	/*
	 * Modified Midpoint Method
	 * 
	 * The modified midpoint method is a second-order method, like (16.1.2), but
	 * with the advantage of requiring (asymptotically for large n) only one
	 * derivative evaluation per step h instead of the two required by
	 * second-order Runge-Kutta. Perhaps there are applications where the
	 * simplicity of (16.3.2), easily coded in-line in some other program,
	 * recommends it. In general, however, use of the modified midpoint method
	 * by itself will be dominated by the embedded Runge-Kutta method with
	 * adaptive stepsize control, as implemented in the preceding section.
	 */
	/**
	 * Modified midpoint step. At xs, input the dependent variable vector y[1..nvar] and its derivative 
	 * vector dydx[1..nvar]. Also input is htot, the total step to be made, and nstep, the 
	 * number of substeps to be used. The output is returned as yout[1..nvar], which need not 
	 * be a distinct array from y; if it is distinct, however, then y and dydx are returned undamaged.
	 * @param y y
	 * @param dydx dydx
	 * @param nvar nvar
	 * @param xs xs
	 * @param htot htot
	 * @param nstep nstep
	 * @param yout yout
	 */
	public static void mmid(double[] y, double[] dydx, int nvar, double xs,
			double htot, int nstep, double[] yout)// , void (*derivs)(float,
													// float[], float[]))
	// Modified midpoint step. At xs, input the dependent variable vector
	// y[1..nvar] and its derivative
	// vector dydx[1..nvar]. Also input is htot, the total step to be made, and
	// nstep, the
	// number of substeps to be used. The output is returned as yout[1..nvar],
	// which need not
	// be a distinct array from y; if it is distinct, however, then y and dydx
	// are returned undamaged.
	{
		int n = 0;
		int i = 0;
		double x = 0.0;
		double swap = 0.0;
		double h2 = 0.0;
		double h = 0.0;
		double[] ym = new double[nvar];
		double[] yn = new double[nvar];
		// ==============================================
		double[] yout1 = new double[nvar];
		// for (i=1;i<=nvar;i++) yout1[i-1]=yout[i-1];//0
		// ============================
		// ym=vector(1,nvar);
		// yn=vector(1,nvar);
		h = htot / nstep; // Stepsize this trip.
		for (i = 1; i <= nvar; i++) {
			ym[i - 1] = y[i - 1];// ym[i]=y[i];
			yn[i - 1] = y[i - 1] + h * dydx[i - 1];// yn[i]=y[i]+h*dydx[i];
													// First step.
		}
		x = xs + h;
		// (*derivs)(x,yn,yout); Will use yout for temporary storage of
		// derivatives.
		yout1 = func.derivF(x, yn);// System.out.println("yout= "+yout[0]);@@@@@@@@@@@@@@@@@@@@@@
		// yout=func.derivF(x,yn);
		h2 = 2.0 * h;
		for (n = 2; n <= nstep; n++) {// General step.
			for (i = 1; i <= nvar; i++) {
				swap = ym[i - 1] + h2 * yout1[i - 1];// swap=ym[i]+h2*yout[i];@@@@@@@@@@@@@@@@@@@@@@@@@@@@
				// swap=ym[i-1]+h2*yout[i-1];
				ym[i - 1] = yn[i - 1];// ym[i]=yn[i];
				yn[i - 1] = swap;// yn[i]=swap;
			}
			x += h;
			// (*derivs)(x,yn,yout);
			yout1 = func.derivF(x, yn);// System.out.println("yout= "+yout[0]);@@@@@@@@@@@@@@@@
			// yout=func.derivF(x,yn);
		}// System.out.println("yout= "+yout[0]);
		for (i = 1; i <= nvar; i++)
			// Last step.
			// yout[i]=0.5*(ym[i]+yn[i]+h*yout[i]);@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
			yout1[i - 1] = 0.5 * (ym[i - 1] + yn[i - 1] + h * yout1[i - 1]);
		// yout[i-1]=0.5*(ym[i-1]+yn[i-1]+h*yout[i-1]);

		for (i = 1; i <= nvar; i++)
			yout[i - 1] = yout1[i - 1];// 0
		// System.out.println("yout= "+yout1[0]);
		// System.out.println("yout= "+yout[0]);
		// free_vector(yn,1,nvar);
		// free_vector(ym,1,nvar);
	}

	/*
	 * Richardson Extrapolation and the Bulirsch-Stoer Method
	 * 
	 * The techniques described in this section are not for differential
	 * equations containing nonsmooth functions. For example, you might have a
	 * differential equation whose right-hand side involves a function that is
	 * evaluated by table look-up and interpolation. If so, go back to
	 * Runge-Kutta with adaptive stepsize choice: That method does an excellent
	 * job of feeling its way through rocky or discontinuous terrain. It is also
	 * an excellent choice for quick-and-dirty, low-accuracy solution of a set
	 * of equations. A second warning is that the techniques in this section are
	 * not particularly good for differential equations that have singular
	 * points inside the interval of integration. A regular solution must tiptoe
	 * very carefully across such points. Runge-Kutta with adaptive stepsize can
	 * sometimes effect this; more generally, there are special techniques
	 * available for such problems, beyond our scope here. Apart from those two
	 * caveats, we believe that the Bulirsch-Stoer method, discussed in this
	 * section, is the best known way to obtain high-accuracy solutions to
	 * ordinary differential equations with minimal computational effort
	 */
	// //float *xx, float htry, float eps,
	/**
	 * Bulirsch-Stoer step with monitoring of local truncation error to ensure accuracy and adjust 
	 * stepsize. Input are the dependent variable vector y[1..nv] and its derivative dydx[1..nv] 
	 * at the starting value of the independent variable x. Also input are the stepsize to be attempted 
	 * htry, the required accuracy eps, and the vector yscal[1..nv] against which the error is 
	 * scaled. On output, y and x are replaced by their new values, hdid is the stepsize that was 
	 * actually accomplished, and hnext is the estimated next stepsize. derivs is the user-supplied 
	 * routine that computes the right-hand side derivatives. Be sure to set htry on successive steps 
	 * to the value of hnext returned from the previous step, as is the case if the routine is called by odeint.
	 * @param y y
	 * @param dydx dydx
	 * @param nv nv
	 * @param xx xx
	 * @param htry htry
	 * @param eps eps
	 * @param yscal yscal
	 */
	public static void bsstep(double[] y, double[] dydx, int nv, double xx,
			double htry, double eps, double[] yscal)// , float *hdid, float
													// *hnext,
	// void (*derivs)(float, float [], float []))
	// Bulirsch-Stoer step with monitoring of local truncation error to ensure
	// accuracy and adjust
	// stepsize. Input are the dependent variable vector y[1..nv] and its
	// derivative dydx[1..nv]
	// at the starting value of the independent variable x. Also input are the
	// stepsize to be attempted
	// htry, the required accuracy eps, and the vector yscal[1..nv] against
	// which the error is
	// scaled. On output, y and x are replaced by their new values, hdid is the
	// stepsize that was
	// actually accomplished, and hnext is the estimated next stepsize. derivs
	// is the user-supplied
	// routine that computes the right-hand side derivatives. Be sure to set
	// htry on successive steps
	// to the value of hnext returned from the previous step, as is the case if
	// the routine is called
	// by odeint.
	{
		// void mmid(float y[], float dydx[], int nvar, float xs, float htot,
		// int nstep, float yout[], void (*derivs)(float, float[], float[]));
		// void pzextr(int iest, float xest, float yest[], float yz[], float
		// dy[],
		// int nv);
		int i = 0;
		int iq = 0;
		int k = 0;
		int kk = 0;
		int km = 0;
		// static int first=1,kmax,kopt;
		// static float epsold = -1.0,xnew;
		// static float a[IMAXX+1];
		// static float alf[KMAXX+1][KMAXX+1];
		// static int nseq[IMAXX+1]={0,2,4,6,8,10,12,14,16,18};

		double eps1 = 0.0;
		double errmax = 0.0;
		double fact = 0.0;
		double h = 0.0;
		double red = 0.0;
		double scale = 0.0;
		double work = 0.0;
		double wrkmin = 0.0;
		double xest = 0.0;

		double[] err = new double[KMAXX];
		double[] yerr = new double[nv];
		double[] ysav = new double[nv];
		double[] yseq = new double[nv];

		xx_bsstep = xx;
		failB = false;

		int reduct = 0;
		int exitflag = 0;

		d_bsstep = new double[nv][KMAXX];// matrix(1,nv,1,KMAXX);
		// err=vector(1,KMAXX);
		x_bsstep = new double[KMAXX];// x=vector(1,KMAXX);
		// yerr=vector(1,nv);
		// ysav=vector(1,nv);
		// yseq=vector(1,nv);
		if (eps != epsold_bsstep) {// A new tolerance, so reinitialize.
			hnext_bsstep = xnew_bsstep = -1.0e29; // “Impossible” values.
			eps1 = SAFE1 * eps;
			a_bsstep[0] = nseq_bsstep[1] + 1;// a[1]=nseq[1]+1; Compute work
												// coefficients Ak.
			for (k = 1; k <= KMAXX; k++)
				a_bsstep[k] = a_bsstep[k - 1] + nseq_bsstep[k + 1];// a[k+1]=a[k]+nseq[k+1];
			for (iq = 2; iq <= KMAXX; iq++) {// Compute ?(k, q).
				for (k = 1; k < iq; k++)
					// alf[k][iq]=pow(eps1,(a[k+1]-a[iq+1])/((a[iq+1]-a[1]+1.0)*(2*k+1)));
					alf_bsstep[k - 1][iq - 1] = Math
							.pow(eps1,
									(a_bsstep[k] - a_bsstep[iq])
											/ ((a_bsstep[iq] - a_bsstep[0] + 1.0) * (2 * k + 1)));
			}
			epsold_bsstep = eps;
			for (kopt_bsstep = 2; kopt_bsstep < KMAXX; kopt_bsstep++)
				// Determine optimal row number for convergence.
				// if (a[kopt+1] > a[kopt]*alf[kopt-1][kopt]) break;
				if (a_bsstep[kopt_bsstep] > a_bsstep[kopt_bsstep - 1]
						* alf_bsstep[kopt_bsstep - 2][kopt_bsstep - 1])
					break;
			kmax_bsstep = kopt_bsstep;
		}
		h = htry;
		for (i = 1; i <= nv; i++)
			ysav[i - 1] = y[i - 1];// ysav[i]=y[i]; Save the starting values.
		if (xx_bsstep != xnew_bsstep || h != (hnext_bsstep)) {// A new stepsize
																// or a new
																// integration:re-establish
																// the order
																// window.
			first_bsstep = 1;
			kopt_bsstep = kmax_bsstep;
		}
		reduct = 0;
		for (;;) {
			for (k = 1; k <= kmax_bsstep; k++) {// Evaluate the sequence of
												// modified midpoint
												// integrations.
				xnew_bsstep = (xx_bsstep) + h;
				if (xnew_bsstep == (xx_bsstep)) {
					failB = true;
					failS = "step size underflow in bsstep";
					return;
					// nrerror("step size underflow in bsstep");
				}
				// mmid(ysav,dydx,nv,*xx,h,nseq[k],yseq)//,derivs);
				mmid(ysav, dydx, nv, xx_bsstep, h, nseq_bsstep[k], yseq);
				// mmid(ysav,dydx,nv,xnew_bsstep,h,nseq_bsstep[k],yseq);
				// xest=SQR(h/nseq[k]); Squared, since error series is even.
				xest = (h / nseq_bsstep[k]) * (h / nseq_bsstep[k]);
				pzextr(k, xest, yseq, y, yerr, nv); // Perform extrapolation.
				if (k != 1) {// Compute normalized error estimate (k).
					errmax = TINY;
					for (i = 1; i <= nv; i++)
						// errmax=FMAX(errmax,fabs(yerr[i]/yscal[i]));
						errmax = Math.max(errmax,
								Math.abs(yerr[i - 1] / yscal[i - 1]));
					errmax /= eps; // Scale error relative to tolerance.
					km = k - 1;
					// err[km]=pow(errmax/SAFE1,1.0/(2*km+1));
					err[km - 1] = Math.pow(errmax / SAFE1, 1.0 / (2 * km + 1));
				}
				if (k != 1 && (k >= kopt_bsstep - 1 || first_bsstep != 0))// first_bsstep))
				{// In order window.
					if (errmax < 1.0) {// Converged.
						exitflag = 1;
						break;
					}
					if (k == kmax_bsstep || k == kopt_bsstep + 1) {// Check for
																	// possible
																	// stepsize
																	// reduction.
						red = SAFE2 / err[km - 1];// red=SAFE2/err[km];
						break;
					}
					// else if (k == kopt && alf[kopt-1][kopt] < err[km])
					else if (k == kopt_bsstep
							&& alf_bsstep[kopt_bsstep - 2][kopt_bsstep - 1] < err[km - 1]) {
						red = 1.0 / err[km - 1];// red=1.0/err[km];
						break;
					}
					// else if (kopt == kmax && alf[km][kmax-1] < err[km])
					else if (kopt_bsstep == kmax_bsstep
							&& alf_bsstep[km - 1][kmax_bsstep - 2] < err[km - 1]) {
						red = alf_bsstep[km - 1][kmax_bsstep - 2] * SAFE2
								/ err[km - 1];// red=alf[km][kmax-1]*SAFE2/err[km];
						break;
					}
					// else if (alf[km][kopt] < err[km])
					else if (alf_bsstep[km - 1][kopt_bsstep - 1] < err[km - 1]) {
						red = alf_bsstep[km - 1][kopt_bsstep - 2] / err[km - 1];// red=alf[km][kopt-1]/err[km];
						break;
					}
				}
			}
			if (exitflag != 0)
				break;// if (exitflag) break;
			red = Math.min(red, REDMIN);// Reduce stepsize by at least REDMIN
										// and at most REDMAX.
			red = Math.max(red, REDMAX);
			h *= red;
			reduct = 1;
		}// Try again.
		xx_bsstep = xnew_bsstep; // Successful step taken.
		hdid_bsstep = h;
		first_bsstep = 0;
		wrkmin = 1.0e35; // Compute optimal row for convergence and
							// corresponding stepsize.
		for (kk = 1; kk <= km; kk++) {
			fact = Math.max(err[kk - 1], SCALMX);// (err[kk],SCALMX);
			work = fact * a_bsstep[kk];// a_bsstep[kk+1];
			if (work < wrkmin) {
				scale = fact;
				wrkmin = work;
				kopt_bsstep = kk + 1;
			}
		}
		hnext_bsstep = h / scale;
		if (kopt_bsstep >= k && kopt_bsstep != kmax_bsstep && reduct == 0)// !reduct)
		{
			// Check for possible order increase, but not if stepsize was just
			// reduced.
			// fact=FMAX(scale/alf[kopt-1][kopt],SCALMX);
			fact = Math.max(scale
					/ alf_bsstep[kopt_bsstep - 2][kopt_bsstep - 1], SCALMX);
			// if (a[kopt+1]*fact <= wrkmin)
			if (a_bsstep[kopt_bsstep] * fact <= wrkmin) {
				hnext_bsstep = h / fact;
				kopt_bsstep++;
			}
		}

		// System.out.println(" inside "+y[0]);
		// free_vector(yseq,1,nv);
		// free_vector(ysav,1,nv);
		// free_vector(yerr,1,nv);
		// free_vector(x,1,KMAXX);
		// free_vector(err,1,KMAXX);
		// free_matrix(d,1,nv,1,KMAXX);
	}

	/**
	 * Use polynomial extrapolation to evaluate nv functions at x = 0 by fitting a polynomial to a sequence of estimates with progressively smaller values x = xest, and 
	 * corresponding function vectors yest[1..nv]. This call is number iest in the sequence of calls. Extrapolated function 
	 * values are output as yz[1..nv], and their estimated error is output as dy[1..nv].
	 * @param iest iest
	 * @param xest xest
	 * @param yest yest
	 * @param yz yz
	 * @param dy dy
	 * @param nv nv
	 */
	public static void pzextr(int iest, double xest, double[] yest,
			double[] yz, double[] dy, int nv)
	// Use polynomial extrapolation to evaluate nv functions at x = 0 by fitting
	// a polynomial to a
	// sequence of estimates with progressively smaller values x = xest, and
	// corresponding function
	// vectors yest[1..nv]. This call is number iest in the sequence of calls.
	// Extrapolated function
	// values are output as yz[1..nv], and their estimated error is output as
	// dy[1..nv].
	{
		int k1 = 0;
		int j = 0;
		double q = 0.0;
		double f2 = 0.0;
		double f1 = 0.0;
		double delta = 0.0;
		double[] c = new double[nv];
		// c=vector(1,nv);
		x_bsstep[iest - 1] = xest;// x[iest]=xest; Save current independent
									// variable.
		for (j = 1; j <= nv; j++)
			dy[j - 1] = yz[j - 1] = yest[j - 1];// dy[j]=yz[j]=yest[j];
		if (iest == 1) {// Store first estimate in first column.
			for (j = 1; j <= nv; j++)
				d_bsstep[j - 1][0] = yest[j - 1];// d[j][1]=yest[j];
		} else {
			for (j = 1; j <= nv; j++)
				c[j - 1] = yest[j - 1];// c[j]=yest[j];
			for (k1 = 1; k1 < iest; k1++) {
				delta = 1.0 / (x_bsstep[iest - k1 - 1] - xest);// delta=1.0/(x[iest-k1]-xest);
				f1 = xest * delta;
				f2 = x_bsstep[iest - k1 - 1] * delta;// f2=x[iest-k1]*delta;
				for (j = 1; j <= nv; j++) {// Propagate tableau 1 diagonal more.
					q = d_bsstep[j - 1][k1 - 1];// q=d[j][k1];
					d_bsstep[j - 1][k1 - 1] = dy[j - 1];// d[j][k1]=dy[j];
					delta = c[j - 1] - q;// delta=c[j]-q;
					dy[j - 1] = f1 * delta;// dy[j]=f1*delta;
					c[j - 1] = f2 * delta;// c[j]=f2*delta;
					yz[j - 1] += dy[j - 1];// yz[j] += dy[j];
				}
			}
			for (j = 1; j <= nv; j++)
				d_bsstep[j - 1][iest - 1] = dy[j - 1];// d[j][iest]=dy[j];
		}
		// free_vector(c,1,nv);
	}

	/*
	 * Current wisdom favors polynomial extrapolation over rational function
	 * extrapolation in the Bulirsch-Stoer method. However, our feeling is that
	 * this view is guided more by the kinds of problems used for tests than by
	 * one method being actually “better.” Accordingly, we provide the optional
	 * routine rzextr for rational function extrapolation, an exact substitution
	 * for pzextr above.
	 */
	/**
	 * Exact substitute for pzextr, but uses diagonal rational function extrapolation instead of polynomial extrapolation. 
	 * Input data are the same as in pzextr.
	 * @param iest iest
	 * @param xest xest
	 * @param yest yest
	 * @param yz yz
	 * @param dy dy
	 * @param nv nv
	 */
	public static void rzextr(int iest, double xest, double[] yest,
			double[] yz, double[] dy, int nv)
	// Exact substitute for pzextr, but uses diagonal rational function
	// extrapolation instead of polynomial
	// extrapolation.
	{
		int k = 0;
		int j = 0;
		double yy = 0.0;
		double v = 0.0;
		double ddy = 0.0;
		double c = 0.0;
		double b1 = 0.0;
		double b = 0.0;
		double[] fx = new double[iest];
		// fx=vector(1,iest);
		x_bsstep[iest - 1] = xest;// x[iest]=xest; Save current independent
									// variable.
		if (iest == 1)
			for (j = 1; j <= nv; j++) {
				yz[j - 1] = yest[j - 1];// yz[j]=yest[j];
				d_bsstep[j - 1][0] = yest[j - 1];// d[j][1]=yest[j];
				dy[j - 1] = yest[j - 1];// dy[j]=yest[j];
			}
		else {
			for (k = 1; k < iest; k++)
				// fx[k+1]=x_bsstep[iest-k]/xest;
				fx[k] = x_bsstep[iest - k - 1] / xest;
			for (j = 1; j <= nv; j++) {// Evaluate next diagonal in tableau.
				v = d_bsstep[j - 1][0];// v=d_bsstep[j][1];
				d_bsstep[j - 1][0] = yy = c = yest[j - 1];// d_bsstep[j][1]=yy=c=yest[j];
				for (k = 2; k <= iest; k++) {
					b1 = fx[k - 1] * v;// b1=fx[k]*v;
					b = b1 - c;
					if (b != 0.0)// if (b)
					{
						b = (c - v) / b;
						ddy = c * b;
						c = b1 * b;
					} else
						// Care needed to avoid division by 0.
						ddy = v;
					if (k != iest)
						v = d_bsstep[j - 1][k - 1];// v=d_bsstep[j][k];
					d_bsstep[j - 1][k - 1] = ddy;// d_bsstep[j][k]=ddy;
					yy += ddy;
				}
				dy[j - 1] = ddy;// dy[j]=ddy;
				yz[j - 1] = yy;// yz[j]=yy;
			}
		}
		// free_vector(fx,1,iest);
	}

	/*
	 * Second-Order Conservative Equations
	 * 
	 * Usually when you have a system of high-order differential equations to
	 * solve it is best to reformulate them as a system of first-order
	 * equations, as discussed in §16.0. There is a particular class of
	 * equations that occurs quite frequently in practice where you can gain
	 * about a factor of two in efficiency by differencing the equations
	 * directly. The equations are second-order systems where the derivative
	 * does not appear on the right-hand side: y'' = f(x, y), y(x0) = y0, y'(x0)
	 * = z0 (16.5.1) As usual, y can denote a vector of values.
	 */

	/**
	 * Stoermer’s rule for integrating y'' = f(x, y) for a system of n = nv/2 equations. On input 
	 * y[1..nv] contains y in its first n elements and y' in its second n elements, all evaluated at 
	 * xs. d2y[1..nv] contains the right-hand side function f (also evaluated at xs) in its first n 
	 * elements. Its second n elements are not referenced. Also input is htot, the total step to be 
	 * taken, and nstep, the number of substeps to be used. The output is returned as yout[1..nv], 
	 * with the same storage arrangement as y. derivs is the user-supplied routine that calculates f.
	 * @param y y
	 * @param d2y d2y
	 * @param nv nv
	 * @param xs xs
	 * @param htot htot
	 * @param nstep nstep
	 * @param yout yout
	 */
	public static void stoerm(double[] y, double[] d2y, int nv, double xs,
			double htot, int nstep, double[] yout)// , void (*derivs)(float,
													// float [], float []))
	// Stoermer’s rule for integrating y'' = f(x, y) for a system of n = nv/2
	// equations. On input
	// y[1..nv] contains y in its first n elements and y' in its second n
	// elements, all evaluated at
	// xs. d2y[1..nv] contains the right-hand side function f (also evaluated at
	// xs) in its first n
	// elements. Its second n elements are not referenced. Also input is htot,
	// the total step to be
	// taken, and nstep, the number of substeps to be used. The output is
	// returned as yout[1..nv],
	// with the same storage arrangement as y. derivs is the user-supplied
	// routine that calculates f.
	{
		int i = 0;
		int n = 0;
		int neqns = 0;
		int nn = 0;
		double h = 0.0;
		double h2 = 0.0;
		double halfh = 0.0;
		double x = 0.0;
		double[] ytemp = new double[nv];
		double[] yout1 = new double[nv];// @@@@@@@@@@@@@@@!!!!!!!!!!!!!!!!!!!
		// ytemp=vector(1,nv);
		h = htot / nstep; // Stepsize this trip.
		halfh = 0.5 * h;
		neqns = nv / 2; // Number of equations.
		for (i = 1; i <= neqns; i++) {// First step.
			n = neqns + i;
			// ytemp[i]=y[i]+(ytemp[n]=h*(y[n]+halfh*d2y[i]));
			ytemp[i - 1] = y[i - 1]
					+ (ytemp[n - 1] = h * (y[n - 1] + halfh * d2y[i - 1]));
		}
		x = xs + h;
		// (*derivs)(x,ytemp,yout); Use yout for temporary storage of
		// derivatives.
		yout1 = func.derivF(x, ytemp);
		h2 = h * h;
		for (nn = 2; nn <= nstep; nn++) {// General step.
			for (i = 1; i <= neqns; i++) {
				// ytemp[i] += (ytemp[(n=neqns+i)] += h2*yout[i]);
				n = neqns + i;
				ytemp[n - 1] += h2 * yout1[i - 1];
				ytemp[i - 1] += ytemp[n - 1];
			}
			x += h;
			// (*derivs)(x,ytemp,yout);
			yout1 = func.derivF(x, ytemp);
		}
		for (i = 1; i <= neqns; i++) {// Last step.
			n = neqns + i;
			// yout[n]=ytemp[n]/h+halfh*yout[i];
			yout1[n - 1] = ytemp[n - 1] / h + halfh * yout1[i - 1];
			yout1[i - 1] = ytemp[i - 1];
		}
		// free_vector(ytemp,1,nv);
		for (i = 1; i <= nv; i++)
			yout[i - 1] = yout1[i - 1];// 0@@@@@@@@@@@@@@@@@@@@@@@@@
	}

	/*
	 * Stiff Sets of Equations
	 * 
	 * As soon as one deals with more than one first-order differential
	 * equation, the possibility of a stiff set of equations arises. Stiffness
	 * occurs in a problem where there are two or more very different scales of
	 * the independent variable on which the dependent variables are changing.
	 * For example, consider the following set of equations [1]: u' = 998u +
	 * 1998v v' = -999u - 1999v (16.6.1) with boundary conditions u(0) = 1 v(0)
	 * = 0 (16.6.2) By means of the transformation u = 2y -z v= -y + z (16.6.3)
	 * we find the solution u = 2e-x - e-1000x v = -e-x + e-1000x (16.6.4) If we
	 * integrated the system (16.6.1) with any of the methods given so far in
	 * this chapter, the presence of the e-1000x term would require a stepsize h
	 * < 1/1000 for the method to be stable (the reason for this is explained
	 * below).
	 * 
	 * So far we have dealt only with implicit methods that are first-order
	 * accurate. While these are very robust, most problems will benefit from
	 * higher-order methods. There are three important classes of higher-order
	 * methods for stiff systems: • Generalizations of the Runge-Kutta method,
	 * of which the most useful are the Rosenbrock methods. The first practical
	 * implementation of these ideas was by Kaps and Rentrop, and so these
	 * methods are also called Kaps-Rentrop methods. • Generalizations of the
	 * Bulirsch-Stoer method, in particular a semi-implicit extrapolation method
	 * due to Bader and Deuflhard. • Predictor-corrector methods, most of which
	 * are descendants of Gear’s backward differentiation method. We shall give
	 * implementations of the first two methods.
	 * 
	 * Rosenbrock Methods
	 */
	/**
	 * Fourth-order Rosenbrock step for integrating stiffo. d.e.’s, with monitoring of local truncation 
	 * error to adjust stepsize. Input are the dependent variable vector y[1..n] and its derivative 
	 * dydx[1..n] at the starting value of the independent variable x. Also input are the stepsize to 
	 * be attempted htry, the required accuracy eps, and the vector yscal[1..n] against which 
	 * the error is scaled. On output, y and x are replaced by their new values, hdid is the stepsize 
	 * that was actually accomplished, and hnext is the estimated next stepsize. derivs is a usersupplied 
	 * routine that computes the derivatives of the right-hand side with respect to x, while 
	 * jacobn (a fixed name) is a user-supplied routine that computes the Jacobi matrix of derivatives 
	 * of the right-hand side with respect to the components of y.
	 * @param y y
	 * @param dydx dydx
	 * @param n n
	 * @param x x
	 * @param htry htry
	 * @param eps eps
	 * @param yscal yscal
	 * @param indxjac indxjac
	 */
	public static void stiff(double[] y, double[] dydx, int n, double x,
			double htry, double eps, double[] yscal, int indxjac)// [], float
																	// *hdid,
																	// float
																	// *hnext,
	// void (*derivs)(float, float [], float []))
	// Fourth-order Rosenbrock step for integrating stiffo. d.e.’s, with
	// monitoring of local truncation
	// error to adjust stepsize. Input are the dependent variable vector y[1..n]
	// and its derivative
	// dydx[1..n] at the starting value of the independent variable x. Also
	// input are the stepsize to
	// be attempted htry, the required accuracy eps, and the vector yscal[1..n]
	// against which
	// the error is scaled. On output, y and x are replaced by their new values,
	// hdid is the stepsize
	// that was actually accomplished, and hnext is the estimated next stepsize.
	// derivs is a usersupplied
	// routine that computes the derivatives of the right-hand side with respect
	// to x, while
	// jacobn (a fixed name) is a user-supplied routine that computes the Jacobi
	// matrix of derivatives
	// of the right-hand side with respect to the components of y.
	{
		x_stiff = x;
		failB = false;
		// void jacobn(float x, float y[], float dfdx[], float **dfdy, int n);
		// void lubksb(float **a, int n, int *indx, float b[]);
		// void ludcmp(float **a, int n, int *indx, float *d);LinAEq.//d
		int i = 0;
		int j = 0;
		int jtry = 0;
		int[] indx = new int[n];
		// double d=0.0;
		double errmax = 0.0;
		double h = 0.0;
		double xsav = 0.0;
		double[][] a = new double[n][n];
		double[] dfdx = new double[n];
		double[][] dfdy = new double[n][n];
		double[] dysav = new double[n];
		double[] err = new double[n];
		double[] g1 = new double[n];
		double[] g2 = new double[n];
		double[] g3 = new double[n];
		double[] g4 = new double[n];
		double[] ysav = new double[n];

		double[] dydx1 = new double[n];// @@@@@@@@@
		// indx=ivector(1,n);
		// a=matrix(1,n,1,n);
		// dfdx=vector(1,n);
		// dfdy=matrix(1,n,1,n);
		// dysav=vector(1,n);
		// err=vector(1,n);
		// g1=vector(1,n);
		// g2=vector(1,n);
		// g3=vector(1,n);
		// g4=vector(1,n);
		// ysav=vector(1,n);

		xsav = (x_stiff); // Save initial values.
		for (i = 1; i <= n; i++) {
			ysav[i - 1] = y[i - 1];// ysav[i]=y[i];
			dysav[i - 1] = dydx[i - 1];// dysav[i]=dydx[i];
		}
		// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
		jacobn(xsav, ysav, dfdx, dfdy, n, indxjac);// jacobn(xsav,ysav,dfdx,dfdy,n);
		// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
		// The user must supply this routine to return the n-by-n matrix dfdy
		// and the vector dfdx.
		h = htry; // Set stepsize to the initial trial value.
		for (jtry = 1; jtry <= MAXTRY_stiff; jtry++) {
			for (i = 1; i <= n; i++) {// Set up the matrix 1 - ?hf.
				for (j = 1; j <= n; j++)
					a[i - 1][j - 1] = -dfdy[i - 1][j - 1];// a[i][j] =
															// -dfdy[i][j];
				a[i - 1][i - 1] += 1.0 / (GAM_stiff * h);// a[i][i] +=
															// 1.0/(GAM_stiff*h);
			}
			LinAEq.ludcmp(a, n, indx);// ,&d); LU decomposition of the matrix.
			for (i = 1; i <= n; i++)
				// Set up right-hand side for g1.
				// g1[i]=dysav[i]+h*C1X_stiff*dfdx[i];
				g1[i - 1] = dysav[i - 1] + h * C1X_stiff * dfdx[i - 1];
			LinAEq.lubksb(a, n, indx, g1); // Solve for g1.
			for (i = 1; i <= n; i++)
				// Compute intermediate values of y and x.
				// y[i]=ysav[i]+A21_stiff*g1[i];
				y[i - 1] = ysav[i - 1] + A21_stiff * g1[i - 1];
			x_stiff = xsav + A2X_stiff * h;
			// (*derivs)(*x,y,dydx); Compute dydx at the intermediate values.
			// dydx=func.derivF(x_stiff,y);
			dydx1 = func.derivF(x_stiff, y);
			for (i = 1; i <= n; i++)
				// Set up right-hand side for g2.
				// g2[i]=dydx[i]+h*C2X_stiff*dfdx[i]+C21_stiff*g1[i]/h;
				g2[i - 1] = dydx1[i - 1] + h * C2X_stiff * dfdx[i - 1]
						+ C21_stiff * g1[i - 1] / h;
			LinAEq.lubksb(a, n, indx, g2);// Solve for g2.
			for (i = 1; i <= n; i++)
				// Compute intermediate values of y and x.
				// y[i]=ysav[i]+A31_stiff*g1[i]+A32_stiff*g2[i];
				y[i - 1] = ysav[i - 1] + A31_stiff * g1[i - 1] + A32_stiff
						* g2[i - 1];
			x_stiff = xsav + A3X_stiff * h;
			// (*derivs)(*x,y,dydx); Compute dydx at the intermediate values.
			// dydx=func.derivF(x_stiff,y);
			dydx1 = func.derivF(x_stiff, y);
			for (i = 1; i <= n; i++)
				// Set up right-hand side for g3.
				// g3[i]=dydx[i]+h*C3X_stiff*dfdx[i]+(C31_stiff*g1[i]+C32_stiff*g2[i])/h;
				g3[i - 1] = dydx1[i - 1] + h * C3X_stiff * dfdx[i - 1]
						+ (C31_stiff * g1[i - 1] + C32_stiff * g2[i - 1]) / h;
			LinAEq.lubksb(a, n, indx, g3);// Solve for g3.
			for (i = 1; i <= n; i++)
				// Set up right-hand side for g4.
				// g4[i]=dydx[i]+h*C4X_stiff*dfdx[i]+(C41_stiff*g1[i]+C42_stiff*g2[i]+C43_stiff*g3[i])/h;
				g4[i - 1] = dydx1[i - 1]
						+ h
						* C4X_stiff
						* dfdx[i - 1]
						+ (C41_stiff * g1[i - 1] + C42_stiff * g2[i - 1] + C43_stiff
								* g3[i - 1]) / h;
			LinAEq.lubksb(a, n, indx, g4); // Solve for g4.
			for (i = 1; i <= n; i++) {// Get fourth-order estimate of y and
										// error estimate.
										// y[i]=ysav[i]+B1_stiff*g1[i]+B2_stiff*g2[i]+B3_stiff*g3[i]+B4_stiff*g4[i];
				y[i - 1] = ysav[i - 1] + B1_stiff * g1[i - 1] + B2_stiff
						* g2[i - 1] + B3_stiff * g3[i - 1] + B4_stiff
						* g4[i - 1];
				// err[i]=E1_stiff*g1[i]+E2_stiff*g2[i]+E3_stiff*g3[i]+E4_stiff*g4[i];
				err[i - 1] = E1_stiff * g1[i - 1] + E2_stiff * g2[i - 1]
						+ E3_stiff * g3[i - 1] + E4_stiff * g4[i - 1];
			}
			x_stiff = xsav + h;
			if (x_stiff == xsav) {
				failB = true;
				failS = "stepsize not significant in stiff";
				return;
				// nrerror("stepsize not significant in stiff");
			}
			errmax = 0.0; // Evaluate accuracy.
			for (i = 1; i <= n; i++)
				// errmax=Math.max(errmax,Math.abs(err[i]/yscal[i]));
				errmax = Math.max(errmax, Math.abs(err[i - 1] / yscal[i - 1]));
			errmax /= eps; // Scale relative to required tolerance.
			if (errmax <= 1.0) {// Step succeeded. Compute size of next step and
								// return.
				hdid_stiff = h;
				hnext_stiff = (errmax > ERRCON_stiff ? SAFETY_stiff * h
						* Math.pow(errmax, PGROW_stiff) : GROW_stiff * h);
				// free_vector(ysav,1,n);
				// free_vector(g4,1,n);
				// free_vector(g3,1,n);
				// free_vector(g2,1,n);
				// free_vector(g1,1,n);
				// free_vector(err,1,n);
				// free_vector(dysav,1,n);
				// free_matrix(dfdy,1,n,1,n);
				// free_vector(dfdx,1,n);
				// free_matrix(a,1,n,1,n);
				// free_ivector(indx,1,n);
				for (i = 1; i <= n; i++)
					dydx[i - 1] = dydx1[i - 1];// 0@@@@@@@@@@@@@@@@@@@@@@@@@
				return;
			} else {// Truncation error too large, reduce stepsize.
				hnext_stiff = SAFETY_stiff * h * Math.pow(errmax, PSHRNK_stiff);
				h = (h >= 0.0 ? Math.max(hnext_stiff, SHRNK_stiff * h) : Math
						.min(hnext_stiff, SHRNK_stiff * h));
			}
		}// Go back and re-try step.
			// nrerror("exceeded MAXTRY in stiff");
		failB = true;
		failS = "exceeded MAXTRY in stiff";
		return;
	}

	/*
	 * As an example of how stiff is used, one can solve the system y1 =
	 * -.013y1 - 1000y1y3 y2 = -2500y2y3 y3 = -.013y1 - 1000y1y3 - 2500y2y3
	 * (16.6.27) with initial conditions y1(0) = 1, y2(0) = 1, y3(0) = 0
	 * (16.6.28) (This is test problem D4 in [4].) We integrate the system up to
	 * x = 50 with an initial stepsize of h = 2.9 × 10-4 using odeint. The
	 * components of C in (16.6.20) are all set to unity. The routines derivs
	 * and jacobn for this problem are given below. Even though the ratio of
	 * largest to smallest decay constants for this problem is around 106, stiff
	 * succeeds in integrating this set in only 29 steps with  = 10-4. By
	 * contrast, the Runge-Kutta routine rkqs requires 51,012 steps!
	 */

	// public static void jacobn(double x, double[] y, double[] dfdx, float
	// **dfdy, int n)
	/**
	 * Test/Example of jacobian.
	 * @param x x
	 * @param y y
	 * @param dfdx dfdx
	 * @param dfdy dfdy
	 * @param n n
	 * @param index index
	 */
	public static void jacobn(double x, double[] y, double[] dfdx,
			double[][] dfdy, int n, int index) {
		int i = 0;
		if (index == 0) {
			for (i = 1; i <= n; i++)
				dfdx[i - 1] = 0.0;// dfdx[i]=0.0;
			dfdy[0][0] = -0.013 - 1000.0 * y[2];// dfdy[1][1] =
												// -0.013-1000.0*y[3];
			dfdy[0][1] = 0.0;// dfdy[1][2]=0.0;
			dfdy[0][2] = -1000.0 * y[0];// dfdy[1][3] = -1000.0*y[1];
			dfdy[1][0] = 0.0;// dfdy[2][1]=0.0;
			dfdy[1][1] = -2500.0 * y[2];// dfdy[2][2] = -2500.0*y[3];
			dfdy[1][2] = -2500.0 * y[1];// dfdy[2][3] = -2500.0*y[2];
			dfdy[2][0] = -0.013 - 1000.0 * y[2];// dfdy[3][1] =
												// -0.013-1000.0*y[3];
			dfdy[2][1] = -2500.0 * y[2];// dfdy[3][2] = -2500.0*y[3];
			dfdy[2][2] = -1000.0 * y[0] - 2500.0 * y[1];// dfdy[3][3] =
														// -1000.0*y[1]-2500.0*y[2];
		}
		if (index == 1) {
			for (i = 1; i <= n; i++)
				dfdx[i - 1] = 0.0;
			// f[0]=y[1];
			// f[1]=-y[0];
			dfdy[0][0] = 0.0;
			dfdy[0][1] = 1.0;
			dfdy[1][0] = -1.0;
			dfdy[1][1] = 0.0;
		}

	}

	/*
	 * Semi-implicit Extrapolation Method
	 * 
	 * The Bulirsch-Stoer method, which discretizes the differential equation
	 * using the modified midpoint rule, does not work for stiff problems.Bader
	 * and Deuflhard [5] discovered a semiimplicit discretization that works
	 * very well and that lends itself to extrapolation exactly as in the
	 * original Bulirsch-Stoer method.
	 */

	// public static void simpr(double[] y, double[] dydx, double[] dfdx, float
	// **dfdy, int n,
	/**
	 * Performs one step of semi-implicit midpoint rule. Input are the dependent variable y[1..n], its 
	 * derivative dydx[1..n], the derivative of the right-hand side with respect to x, dfdx[1..n], 
	 * and the Jacobian dfdy[1..n][1..n] at xs. Also input are htot, the total step to be taken, and nstep, the number of substeps to be used. The output is returned as 
	 * yout[1..n]. derivs is the user-supplied routine that calculates dydx.
	 * @param y y
	 * @param dydx dydx
	 * @param dfdx dfdx
	 * @param dfdy dfdy
	 * @param n n
	 * @param xs xs
	 * @param htot htot
	 * @param nstep nstep
	 * @param yout yout
	 */
	public static void simpr(double[] y, double[] dydx, double[] dfdx,
			double[][] dfdy, int n, double xs, double htot, int nstep,
			double[] yout)// ,
	// void (*derivs)(float, float [], float []))
	// Performs one step of semi-implicit midpoint rule. Input are the dependent
	// variable y[1..n], its
	// derivative dydx[1..n], the derivative of the right-hand side with respect
	// to x, dfdx[1..n],
	// and the Jacobian dfdy[1..n][1..n] at xs. Also input are htot, the total
	// step to be taken,
	// and nstep, the number of substeps to be used. The output is returned as
	// yout[1..n].
	// derivs is the user-supplied routine that calculates dydx.
	{
		// void lubksb(float **a, int n, int *indx, float b[]);
		// void ludcmp(float **a, int n, int *indx, float *d);
		int i = 0;
		int j = 0;
		int nn = 0;
		int[] indx = new int[n];
		// double d=0.0;
		double h = 0.0;
		double x = 0.0;
		double[][] a = new double[n][n];
		double[] del = new double[n];
		double[] ytemp = new double[n];

		double[] yout1 = new double[n];
		// indx=ivector(1,n);
		// a=matrix(1,n,1,n);
		// del=vector(1,n);
		// ytemp=vector(1,n);
		h = htot / nstep; // Stepsize this trip.
		for (i = 1; i <= n; i++) {// Set up the matrix 1 - hf.
			for (j = 1; j <= n; j++)
				a[i - 1][j - 1] = -h * dfdy[i - 1][j - 1];// a[i][j] =
															// -h*dfdy[i][j];
			++a[i - 1][i - 1];// ++a[i][i];
		}
		// ludcmp(a,n,indx,&d); LU decomposition of the matrix.
		LinAEq.ludcmp(a, n, indx);
		for (i = 1; i <= n; i++)
			// Set up right-hand side for first step. Use yout for temporary
			// storage.
			// yout[i]=h*(dydx[i]+h*dfdx[i]);
			yout[i - 1] = h * (dydx[i - 1] + h * dfdx[i - 1]);
		LinAEq.lubksb(a, n, indx, yout);
		for (i = 1; i <= n; i++)
			// First step.
			ytemp[i - 1] = y[i - 1] + (del[i - 1] = yout[i - 1]);// ytemp[i]=y[i]+(del[i]=yout[i]);
		x = xs + h;
		// (*derivs)(x,ytemp,yout); Use yout for temporary storage of
		// derivatives.
		yout1 = func.derivF(x, ytemp);
		for (nn = 2; nn <= nstep; nn++) {// General step.
			for (i = 1; i <= n; i++)
				// Set up right-hand side for general step.
				// yout[i]=h*yout[i]-del[i];
				yout1[i - 1] = h * yout1[i - 1] - del[i - 1];
			LinAEq.lubksb(a, n, indx, yout1);// ,yout);
			for (i = 1; i <= n; i++)
				// ytemp[i] += (del[i] += 2.0*yout[i]);
				ytemp[i - 1] += (del[i - 1] += 2.0 * yout1[i - 1]);
			x += h;
			// (*derivs)(x,ytemp,yout);
			yout1 = func.derivF(x, ytemp);
		}
		for (i = 1; i <= n; i++)
			// Set up right-hand side for last step.
			// yout[i]=h*yout[i]-del[i];
			yout1[i - 1] = h * yout1[i - 1] - del[i - 1];
		LinAEq.lubksb(a, n, indx, yout1);// ,yout);
		for (i = 1; i <= n; i++)
			// Take last step.
			// yout[i] += ytemp[i];
			yout1[i - 1] += ytemp[i - 1];

		for (i = 1; i <= n; i++)
			yout[i - 1] = yout1[i - 1];// @@@@@@@@@@@@@@@@@
		// free_vector(ytemp,1,n);
		// free_vector(del,1,n);
		// free_matrix(a,1,n,1,n);
		// free_ivector(indx,1,n);
	}

	/**
	 * Use polynomial extrapolation to evaluate nv functions at x = 0 by fitting a polynomial to a 
	 * sequence of estimates with progressively smaller values x = xest, and corresponding function 
	 * vectors yest[1..nv]. This call is number iest in the sequence of calls. Extrapolated function 
	 * values are output as yz[1..nv], and their estimated error is output as dy[1..nv].
	 * @param iest iest
	 * @param xest xest
	 * @param yest yest
	 * @param yz yz
	 * @param dy dy
	 * @param nv nv
	 */
	public static void pzextr1(int iest, double xest, double[] yest,
			double[] yz, double[] dy, int nv)
	// Use polynomial extrapolation to evaluate nv functions at x = 0 by fitting
	// a polynomial to a
	// sequence of estimates with progressively smaller values x = xest, and
	// corresponding function
	// vectors yest[1..nv]. This call is number iest in the sequence of calls.
	// Extrapolated function
	// values are output as yz[1..nv], and their estimated error is output as
	// dy[1..nv].
	{
		int k1 = 0;
		int j = 0;
		double q = 0.0;
		double f2 = 0.0;
		double f1 = 0.0;
		double delta = 0.0;
		double[] c = new double[nv];
		// c=vector(1,nv);
		x_stifbs[iest - 1] = xest;// x[iest]=xest; Save current independent
									// variable.
		for (j = 1; j <= nv; j++)
			dy[j - 1] = yz[j - 1] = yest[j - 1];// dy[j]=yz[j]=yest[j];
		if (iest == 1) {// Store first estimate in first column.
			for (j = 1; j <= nv; j++)
				d_stifbs[j - 1][0] = yest[j - 1];// d[j][1]=yest[j];
		} else {
			for (j = 1; j <= nv; j++)
				c[j - 1] = yest[j - 1];// c[j]=yest[j];
			for (k1 = 1; k1 < iest; k1++) {
				delta = 1.0 / (x_stifbs[iest - k1 - 1] - xest);// delta=1.0/(x[iest-k1]-xest);
				f1 = xest * delta;
				f2 = x_stifbs[iest - k1 - 1] * delta;// f2=x[iest-k1]*delta;
				for (j = 1; j <= nv; j++) {// Propagate tableau 1 diagonal more.
					q = d_stifbs[j - 1][k1 - 1];// q=d[j][k1];
					d_stifbs[j - 1][k1 - 1] = dy[j - 1];// d[j][k1]=dy[j];
					delta = c[j - 1] - q;// delta=c[j]-q;
					dy[j - 1] = f1 * delta;// dy[j]=f1*delta;
					c[j - 1] = f2 * delta;// c[j]=f2*delta;
					yz[j - 1] += dy[j - 1];// yz[j] += dy[j];
				}
			}
			for (j = 1; j <= nv; j++)
				d_stifbs[j - 1][iest - 1] = dy[j - 1];// d[j][iest]=dy[j];
		}
		// free_vector(c,1,nv);
	}

	/*
	 * The routine simpr is intended to be used in a routine stifbs that is
	 * almost exactly the same as bsstep. The only differences are: • The
	 * stepsize sequence is n = 2, 6, 10, 14, 22, 34, 50, . . . , (16.6.35)
	 * where each member differs from its predecessor by the smallest multiple
	 * of 4 that makes the ratio of successive terms be ? 5 7 . The parameter
	 * KMAXX is taken to be 7. • The work per unit step now includes the cost of
	 * Jacobian evaluations as well as function evaluations. We count one
	 * Jacobian evaluation as equivalent to N function evaluations, where N is
	 * the number of equations. • Once again the user-supplied routine derivs is
	 * a dummy argument and so can have any name. However, to maintain
	 * “plug-compatibility” with rkqs, bsstep and stiff, the routine jacobn is
	 * not an argument and must have exactly this name. It is called once per
	 * step to return f (dfdy) and ?f/?x (dfdx) as functions of x and y. Here
	 * is the routine, with comments pointing out only the differences from
	 * bsstep:
	 */
	/**
	 * Semi-implicit extrapolation step for integrating stiffo .d.e.’s, with monitoring of local truncation 
	 * error to adjust stepsize. Input are the dependent variable vector y[1..nv] and its derivative 
	 * dydx[1..nv] at the starting value of the independent variable x. Also input are the stepsize 
	 * to be attempted htry, the required accuracy eps, and the vector yscal[1..nv] against which the error is scaled. On output, y and x are replaced by their new 
	 * values, hdid is the stepsize that was actually accomplished, and hnext is the estimated next stepsize. derivs 
	 * is a user-supplied routine that computes the derivatives of the right-hand side with respect to x, while jacobn (a fixed name) is a user-supplied routine that computes 
	 * the Jacobi matrix of derivatives of the right-hand side with respect to the components of y. Be sure to set htry 
	 * on successive steps to the value of hnext returned from the previous step, as is the case if the routine is called by odeint.
	 * @param y y
	 * @param dydx dydx
	 * @param nv nv
	 * @param xx xx
	 * @param htry htry
	 * @param eps eps
	 * @param yscal yscal
	 * @param indxjac indxjac
	 */
	public static void stifbs(double[] y, double[] dydx, int nv, double xx,
			double htry, double eps, double[] yscal, int indxjac)// , float
																	// *hdid,
																	// float
																	// *hnext,
	// void (*derivs)(float, float [], float []))
	// Semi-implicit extrapolation step for integrating stiffo .d.e.’s, with
	// monitoring of local truncation
	// error to adjust stepsize. Input are the dependent variable vector
	// y[1..nv] and its derivative
	// dydx[1..nv] at the starting value of the independent variable x. Also
	// input are the stepsize
	// to be attempted htry, the required accuracy eps, and the vector
	// yscal[1..nv] against
	// which the error is scaled. On output, y and x are replaced by their new
	// values, hdid is the
	// stepsize that was actually accomplished, and hnext is the estimated next
	// stepsize. derivs
	// is a user-supplied routine that computes the derivatives of the
	// right-hand side with respect to
	// x, while jacobn (a fixed name) is a user-supplied routine that computes
	// the Jacobi matrix of
	// derivatives of the right-hand side with respect to the components of y.
	// Be sure to set htry
	// on successive steps to the value of hnext returned from the previous
	// step, as is the case if the
	// routine is called by odeint.
	{
		// void jacobn(float x, float y[], float dfdx[], float **dfdy, int n);
		// void simpr(float y[], float dydx[], float dfdx[], float **dfdy,
		// int n, float xs, float htot, int nstep, float yout[],
		// void (*derivs)(float, float [], float []));
		// void pzextr(int iest, float xest, float yest[], float yz[], float
		// dy[],
		// int nv);
		xx_stifbs = xx;
		failB = false;
		int i = 0;
		int iq = 0;
		int k = 0;
		int kk = 0;
		int km = 0;
		// static int first=1,kmax,kopt,nvold = -1;
		// static float epsold = -1.0,xnew;
		double eps1 = 0.0;
		double errmax = 0.0;
		double fact = 0.0;
		double h = 0.0;
		double red = 0.0;
		double scale = 0.0;
		double work = 0.0;
		double wrkmin = 0.0;
		double xest = 0.0;
		// double[] dfdx,
		// double[][] dfdy,
		// double[] err,
		// double[] yerr,
		// double[] ysav,
		// double[] yseq;
		// static float a[IMAXX+1];
		// static float alf[KMAXX+1][KMAXX+1];
		// static int nseq[IMAXX+1]={0,2,6,10,14,22,34,50,70};
		// Sequence is different from bsstep.
		int reduct = 0;
		int exitflag = 0;
		d_stifbs = new double[nv][KMAXX1];
		x_stifbs = new double[KMAXX1];
		double[] dfdx = new double[nv];
		double[][] dfdy = new double[nv][nv];
		double[] err = new double[KMAXX1];
		double[] yerr = new double[nv];
		double[] ysav = new double[nv];
		double[] yseq = new double[nv];

		// d=matrix(1,nv,1,KMAXX);
		// dfdx=vector(1,nv);
		// dfdy=matrix(1,nv,1,nv);
		// err=vector(1,KMAXX);
		// x=vector(1,KMAXX);
		// yerr=vector(1,nv);
		// ysav=vector(1,nv);
		// yseq=vector(1,nv);
		if (eps != epsold_stifbs || nv != nvold_stifbs) {// Reinitialize also if
															// nv has changed.
			hnext_stifbs = xnew_stifbs = -1.0e29;
			eps1 = SAFE1 * eps;
			a_stifbs[0] = nseq_stifbs[1] + 1;// a[1]=nseq[1]+1;
			for (k = 1; k <= KMAXX1; k++)
				// a[k+1]=a[k]+nseq[k+1];
				a_stifbs[k] = a_stifbs[k - 1] + nseq_stifbs[k + 1];
			for (iq = 2; iq <= KMAXX1; iq++) {
				for (k = 1; k < iq; k++)
					// alf[k][iq]=pow(eps1,((a[k+1]-a[iq+1])/((a[iq+1]-a[1]+1.0)*(2*k+1))));
					alf_stifbs[k - 1][iq - 1] = Math.pow(eps1,
							((a_stifbs[k] - a_stifbs[iq]) / ((a_stifbs[iq]
									- a_stifbs[0] + 1.0) * (2 * k + 1))));
			}
			epsold_stifbs = eps;
			nvold_stifbs = nv;// Save nv.
			a_stifbs[0] += nv;// a[1] += nv; Add cost of Jacobian evaluations to
								// work coefficients.
			for (k = 1; k <= KMAXX1; k++)
				a_stifbs[k] = a_stifbs[k - 1] + nseq_stifbs[k + 1];// a[k+1]=a[k]+nseq[k+1];
			for (kopt_stifbs = 2; kopt_stifbs < KMAXX1; kopt_stifbs++)
				// if (a[kopt+1] > a[kopt]*alf[kopt-1][kopt]) break;
				if (a_stifbs[kopt_stifbs] > a_stifbs[kopt_stifbs - 1]
						* alf_stifbs[kopt_stifbs - 2][kopt_stifbs - 1])
					break;
			kmax_stifbs = kopt_stifbs;
		}
		h = htry;
		for (i = 1; i <= nv; i++)
			ysav[i - 1] = y[i - 1];// ysav[i]=y[i];
		// jacobn(*xx,y,dfdx,dfdy,nv); Evaluate Jacobian.
		jacobn(xx_stifbs, y, dfdx, dfdy, nv, indxjac);
		if (xx_stifbs != xnew_stifbs || h != (hnext_stifbs)) {
			first_stifbs = 1;
			kopt_stifbs = kmax_stifbs;
		}
		reduct = 0;
		for (;;) {
			for (k = 1; k <= kmax_stifbs; k++) {
				xnew_stifbs = (xx_stifbs) + h;
				if (xnew_stifbs == (xx_stifbs)) {
					failB = true;
					failS = "step size underflow in stifbs";
					return;
					// nrerror("step size underflow in stifbs");
				}
				// simpr(ysav,dydx,dfdx,dfdy,nv,*xx,h,nseq[k],yseq,derivs);
				simpr(ysav, dydx, dfdx, dfdy, nv, xx_stifbs, h, nseq_stifbs[k],
						yseq);// ,derivs);
				// Semi-implicit midpoint rule.
				// xest=SQR(h/nseq[k]); The rest of the routine is identical to
				// bsstep.
				xest = (h / nseq_stifbs[k]) * (h / nseq_stifbs[k]);
				pzextr1(k, xest, yseq, y, yerr, nv);
				if (k != 1) {
					errmax = TINY;
					for (i = 1; i <= nv; i++)
						// errmax=FMAX(errmax,fabs(yerr[i]/yscal[i]));
						errmax = Math.max(errmax,
								Math.abs(yerr[i - 1] / yscal[i - 1]));
					errmax /= eps;
					km = k - 1;
					// err[km]=pow(errmax/SAFE1,1.0/(2*km+1));
					err[km - 1] = Math.pow(errmax / SAFE1, 1.0 / (2 * km + 1));
				}
				// if (k != 1 && (k >= kopt_stifbs-1 || first_stifbs))
				if (k != 1 && (k >= kopt_stifbs - 1 || first_stifbs != 0)) {
					if (errmax < 1.0) {
						exitflag = 1;
						break;
					}
					if (k == kmax_stifbs || k == kopt_stifbs + 1) {
						red = SAFE2 / err[km - 1];// red=SAFE2/err[km];
						break;
					}
					// else if (k == kopt && alf[kopt-1][kopt] < err[km])
					else if (k == kopt_stifbs
							&& alf_stifbs[kopt_stifbs - 2][kopt_stifbs - 1] < err[km - 1]) {
						red = 1.0 / err[km - 1];// red=1.0/err[km];
						break;
					}
					// else if (kopt == kmax && alf[km][kmax-1] < err[km])
					else if (kopt_stifbs == kmax_stifbs
							&& alf_stifbs[km - 1][kmax_stifbs - 2] < err[km - 1]) {
						// red=alf[km][kmax-1]*SAFE2/err[km];
						red = alf_stifbs[km - 1][kmax_stifbs - 2] * SAFE2
								/ err[km - 1];
						break;
					}
					// else if (alf[km][kopt] < err[km])
					else if (alf_stifbs[km - 1][kopt_stifbs - 1] < err[km - 1]) {
						// red=alf[km][kopt-1]/err[km];
						red = alf_stifbs[km - 1][kopt_stifbs - 2] / err[km - 1];
						break;
					}
				}
			}
			// if (exitflag) break;
			if (exitflag != 0)
				break;
			red = Math.min(red, REDMIN);
			red = Math.max(red, REDMAX);
			h *= red;
			reduct = 1;
		}
		xx_stifbs = xnew_stifbs;
		hdid_stifbs = h;
		first_stifbs = 0;
		wrkmin = 1.0e35;
		for (kk = 1; kk <= km; kk++) {
			fact = Math.max(err[kk - 1], SCALMX);// fact=Math.max(err[kk],SCALMX);
			work = fact * a_stifbs[kk];// work=fact*a[kk+1];
			if (work < wrkmin) {
				scale = fact;
				wrkmin = work;
				kopt_stifbs = kk + 1;
			}
		}
		hnext_stifbs = h / scale;
		// if (kopt_stifbs >= k && kopt_stifbs != kmax_stifbs && !reduct)
		if (kopt_stifbs >= k && kopt_stifbs != kmax_stifbs && reduct == 0) {
			// fact=Math.max(scale/alf[kopt-1][kopt],SCALMX);
			fact = Math.max(scale
					/ alf_stifbs[kopt_stifbs - 2][kopt_stifbs - 1], SCALMX);
			// if (a[kopt+1]*fact <= wrkmin)
			if (a_stifbs[kopt_stifbs] * fact <= wrkmin) {
				hnext_stifbs = h / fact;
				kopt_stifbs++;
			}
		}
		// free_vector(yseq,1,nv);
		// free_vector(ysav,1,nv);
		// free_vector(yerr,1,nv);
		// free_vector(x,1,KMAXX);
		// free_vector(err,1,KMAXX);
		// free_matrix(dfdy,1,nv,1,nv);
		// free_vector(dfdx,1,nv);
		// free_matrix(d,1,nv,1,KMAXX);
	}
	/*
	 * The routine stifbs is an excellent routine for all stiff problems,
	 * competitive with the best Gear-type routines.
	 */
	//==============
	public static boolean printB=true;
	private NumberFormat nf = NumberFormat.getInstance(Locale.US);
	//metoda Runge-Kutta de ordin 4
	//intra valoarea x pentru evaluarea lui y si a derivatelor sale
	//pasul de propagare a solutiei
	//tabloul conditiilor initiale(y0,y0',....
	//ordinul ecuatiei diferentiale (1,2,....
	/**
	 * Other implementation of Runke kutta (order 4) method. x is used for evaluate y and its derivatives, h is the step to be taken, 
	 * y0 is the array of initial condition (y0, y0',...y0(n)). n is the order of differential equation, the array size.
	 * @param x x
	 * @param h h
	 * @param y0 y0
	 * @param n n
	 * @return the result, i.e. the value y and all derivatives y',....
	 */
	    public double[] RK4(double x,double h, double[] y0, int n)
	    {
	        double[] y=new double[n];
	        double[] yy=new double[n];
			double hh=0.5*h;
			double[] f1=func.derivF(x,y0);
			for (int i=1; i<=n; i++)
				yy[i-1]=y0[i-1]+hh*f1[i-1];
	        double[] f2=func.derivF(x+hh,yy);
			for (int i=1; i<=n; i++)
				yy[i-1]=y0[i-1]+hh*f2[i-1];
	        double[] f3=func.derivF(x+hh,yy);//df.func(x+hh,yy);
			for (int i=1; i<=n; i++)
			{
				yy[i-1]=y0[i-1]+h*f3[i-1];
				f2[i-1]=f2[i-1]+f3[i-1];
			}
			f3=func.derivF(x+h,yy);
			hh=h/6.0;
			for (int i=1; i<=n; i++)
			{
				y[i-1]=y0[i-1]+hh*(f1[i-1]+2*f2[i-1]+f3[i-1]);

				if (printB)
				{
					func.printSequence("============RK4 evaluation===================");
					func.printSequence("Results for '"+new Integer(i-1)+"' derivate= "+nf.format(y[i-1])+
					                 " at x= x0("+nf.format(x)+")+step("+nf.format(h)+")= "+nf.format(x+h));
					func.printSequence("=============================================");
				}

			}
			return y;
		}

	//solutia pe un sir de p puncte echidistante delimitate de intervalul a<x<b
	    /**
	     * Runge-Kutta method of order 4 for an array of p equaly spaced points in the interval [a,b].
	     * @param a a
	     * @param b b
	     * @param p p
	     * @param y0 y0
	     * @param n n
	     * @return the solution, i.e. the y [][0] and all its derivatives [][1]....and for all p points.
	     */
	    public double[][] allRK4(double a, double b, int p, double[] y0, int n)
	    {
			printB=false;
			double h=(b-a)/(p-1);//evaluare pas
			double[] x=new double[p];
			for (int i=0; i<p; i++)
				x[i]=a+h*i;
			//gata cu constructia retelei de puncte
			double[][] y=new double[p][n];
			double[] ytmp=new double[n];
			//[][0]-solutia, [][1]-prima derivata,....[n-1] derivata.
			for (int i=0; i<n; i++)
			{
				y[0][i]=y0[i];//initial conditions at x0=a!!!
				ytmp[i]=y0[i];//to not alter the initial array, y0!!
			}
			for (int j=1; j<p; j++)
			{
				double[] ys=RK4(x[j-1],h,ytmp,n);
				for (int i=0; i<n; i++)
				{
					y[j][i]=ys[i];
					ytmp[i]=ys[i];
				}
			}

			printB=true;
			for (int i=1; i<=n; i++)
			{
				if (printB)
				{
					func.printSequence("============RK4 evaluation===================");
					func.printSequence("Results for '"+new Integer(i-1)+"' derivate= "+nf.format(y[p-1][i-1])+
					                 " at x= "+nf.format(b));
					func.printSequence("=============================================");
				}

			}

			return y;
		}
}
