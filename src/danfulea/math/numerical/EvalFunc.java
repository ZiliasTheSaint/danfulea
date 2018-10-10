package danfulea.math.numerical;

import java.text.NumberFormat;
import java.util.Locale;

/**
 * Evaluation of functions class. 
 * Based on literature, including Numerical Recipes in C (Cambridge Univ.).
 * @author Dan Fulea, 05 OCT. 2006.
 */
public class EvalFunc {
	private NumberFormat nf = NumberFormat.getInstance(Locale.US);
	//private String pattern = "0.###E0";
	//private DecimalFormatSymbols dfs = new DecimalFormatSymbols(Locale.US);
	//private DecimalFormat nff = nff = new DecimalFormat(pattern, dfs);
	private int idigits = 3;

	public static boolean failB = false;
	public static String failS = "";
	public static double sumeul = 0.0;
	public static int nterm = 0;

	private Function func;
	public static double CON = 1.4;// Stepsize is decreased by CON at each
									// iteration.
	public static double CON2 = CON * CON;
	public static double BIG = 1.0e30;
	public static int NTAB = 10;// Sets maximum size of tableau.
	public static double SAFE = 2.0;// Return when error is SAFE worse than the
									// best so
	// far.
	public static double PI = 3.141592653589793;
	public static int NPFAC = 8;
	public static int MAXIT = 5;
	public static double PIO2 = 3.141592653589793 / 2.0;

	/**
	 * Constructor
	 * @param func a class implementing Function interface.
	 */
	public EvalFunc(Function func) {
		this.func = func;
		set_defaults();
	}

	public EvalFunc() {
		// null constructor for predef functions
	}

	/**
	 * Choose the number of decimal digits for number format. E.g. if i=2 then a number 4.21456 would be displayed as 4.21.
	 * @param i the number of digits.
	 */
	public void set_defaults(int i) {
		idigits = i;
		nf.setMinimumFractionDigits(idigits);// 3);//default is 2!!
		nf.setMaximumFractionDigits(idigits);// (3);//default is 2!!
		nf.setGroupingUsed(false);// no 4,568.02 but 4568.02
	}

	/**
	 * Default number of digits, which is 3.
	 */
	public void set_defaults() {
		idigits = 3;
		nf.setMinimumFractionDigits(idigits);// 3);//default is 2!!
		nf.setMaximumFractionDigits(idigits);// (3);//default is 2!!
		nf.setGroupingUsed(false);// no 4,568.02 but 4568.02
	}

	/*
	 * There is an elegant and subtle implementation of Euler’s transformation
	 * due to van Wijngaarden [1]: It incorporates the terms of the original
	 * alternating series one at a time, in order. For each incorporation it
	 * either increases p by 1, equivalent to computing one further difference
	 * (5.1.6), or else retroactively increases n by 1, without having to redo
	 * all the difference calculations based on the old n value! The decision as
	 * to which to increase, n or p, is taken in such a way as to make the
	 * convergence most rapid. Van Wijngaarden’s technique requires only one
	 * vector of saved partial differences.
	 */
	// public static void eulsum(float *sum, float term, int jterm, float
	// wksp[])
	/**
	 * Convergence of a series. Euler transformation. Incorporates into sum the jterm’th term, with value term, of an 
	 * alternating series. sum is input as the previous partial sum, and is output as the new partial sum. The first call to this 
	 * routine, with the first term in the series, should be with jterm=1. On the second call, term should be set to the second term 
	 * of the series, with sign opposite to that of the first call, and jterm should be 2. And so on. 
	 * wksp is a workspace array provided by the calling program, dimensioned at least as large as the maximum number of terms to be 
	 * incorporated.
	 * @param sum sum
	 * @param term term
	 * @param jterm jterm
	 * @param wksp wksp
	 */
	public static void eulsum(double sum, double term, int jterm, double[] wksp)
	// Incorporates into sum the jterm’th term, with value term, of an
	// alternating series. sum is
	// input as the previous partial sum, and is output as the new partial sum.
	// The first call to this
	// routine, with the first term in the series, should be with jterm=1. On
	// the second call, term
	// should be set to the second term of the series, with sign opposite to
	// that of the first call, and
	// jterm should be 2. And so on. wksp is a workspace array provided by the
	// calling program,
	// dimensioned at least as large as the maximum number of terms to be
	// incorporated.
	{
		int j = 0;
		// static int nterm;
		double tmp = 0.0;
		double dum = 0.0;

		sumeul = sum;
		if (jterm == 1) {// Initialize:
			nterm = 1;// Number of saved differences in wksp.
			sumeul = 0.5 * (wksp[0] = term);// *sum=0.5*(wksp[1]=term); //Return
											// first estimate.
		} else {
			tmp = wksp[0];// tmp=wksp[1];
			wksp[0] = term;// wksp[1]=term;
			for (j = 1; j <= nterm - 1; j++) {
				// Update saved quantities by van Wijngaarden’s algorithm.
				dum = wksp[j];// dum=wksp[j+1];
				wksp[j] = 0.5 * (wksp[j - 1] + tmp);// wksp[j+1]=0.5*(wksp[j]+tmp);
				tmp = dum;
			}
			wksp[nterm] = 0.5 * (wksp[nterm - 1] + tmp);// wksp[nterm+1]=0.5*(wksp[nterm]+tmp);
			if (Math.abs(wksp[nterm]) <= Math.abs(wksp[nterm - 1]))// if
																	// (fabs(wksp[nterm+1])
																	// <=
																	// fabs(wksp[nterm]))
																	// //Favorable
																	// to
																	// increase
																	// p,
				sumeul += (0.5 * wksp[++nterm - 1]);// *sum +=
													// (0.5*wksp[++nterm]);
													// //and the table becomes
													// longer.
			else
				// Favorable to increase n,
				sumeul += wksp[nterm];// *sum += wksp[nterm+1]; the table
										// doesn’t become longer.
		}
	}

	// The above is a lit a bit weird!! where nterm is increased?=>see line 117!!. It is incresed in line 138 NEVERMIND!
	/*
	 * A polynomial of degree N is represented numerically as a stored array of
	 * coefficients, c[j] with j= 0, . . .,N. We will always take c[0] to be the
	 * constant term in the polynomial, c[N] the coefficient of xN; but of
	 * course other conventions are possible. There are two kinds of
	 * manipulations that you can do with a polynomial: numerical manipulations
	 * (such as evaluation), where you are given the numerical value of its
	 * argument, or algebraic manipulations, where you want to transform the
	 * coefficient array in some way without choosing any particular argument.
	 * Let’s start with the numerical. We assume that you know enough never to
	 * evaluate a polynomial this way:
	 * p=c[0]+c[1]*x+c[2]*x*x+c[3]*x*x*x+c[4]*x*x*x*x; or (even worse!),
	 * p=c[0]+c[1]*x+c[2]*pow(x,2.0)+c[3]*pow(x,3.0)+c[4]*pow(x,4.0); Come the
	 * (computer) revolution, all persons found guilty of such criminal behavior
	 * will be summarily executed, and their programs won’t be! It is a matter
	 * of taste, however, whether to write
	 * p=c[0]+x*(c[1]+x*(c[2]+x*(c[3]+x*c[4]))); or
	 * p=(((c[4]*x+c[3])*x+c[2])*x+c[1])*x+c[0]; If the number of coefficients
	 * c[0..n] is large, one writes p=c[n]; for(j=n-1;j>=0;j--) p=p*x+c[j]; or
	 * p=c[j=n]; while (j>0) p=p*x+c[--j]; Another useful trick is for
	 * evaluating a polynomial P(x) and its derivative dP(x)/dx simultaneously:
	 * p=c[n]; dp=0.0; for(j=n-1;j>=0;j--) {dp=dp*x+p; p=p*x+c[j];} or p=c[j=n];
	 * dp=0.0; while (j>0) {dp=dp*x+p; p=p*x+c[--j];} which yields the
	 * polynomial as p and its derivative as dp.
	 */

	/**
	 * A polynomial of degree N is represented numerically as a stored array of coefficients, c[j] with j= 0, . . .,N.
	 * We will always take c[0] to be the constant term in the polynomial, c[N] the coefficient of xN; 
	 * Given the nc+1 coefficients of a polynomial of degree nc as an array c[0..nc] with c[0] 
	 * being the constant term, and given a value x, and given a value nd greater than 1, this routine returns the 
	 * polynomial evaluated at x as pd[0] and nd derivatives as pd[1..nd].
	 * @param c c
	 * @param nc nc
	 * @param x x
	 * @param pd pd
	 * @param nd nd
	 */
	public static void ddpoly(double[] c, int nc, double x, double[] pd, int nd)
	// Given the nc+1 coefficients of a polynomial of degree nc as an array
	// c[0..nc] with c[0]
	// being the constant term, and given a value x, and given a value nd>1,
	// this routine returns the
	// polynomial evaluated at x as pd[0] and nd derivatives as pd[1..nd].
	{
		int nnd = 0;
		int j = 0;
		int i = 0;
		double cnst = 1.0;

		pd[0] = c[nc];
		for (j = 1; j <= nd; j++)
			pd[j] = 0.0;
		for (i = nc - 1; i >= 0; i--) {
			nnd = (nd < (nc - i) ? nd : nc - i);
			for (j = nnd; j >= 1; j--)
				pd[j] = pd[j] * x + pd[j - 1];
			pd[0] = pd[0] * x + c[i];
		}
		for (i = 2; i <= nd; i++) {// After the first derivative, factorial
									// constants come in.
			cnst *= i;
			pd[i] *= cnst;
		}
	}

	/**
	 * Given the n+1 coefficients of a polynomial of degree n in u[0..n], and the nv+1 coefficients 
	 * of another polynomial of degree nv in v[0..nv], divide the polynomial u by the polynomial v 
	 * (“u”/“v”)giving a quotient polynomial whose coefficients are returned in q[0..n], and a 
	 * remainder polynomial whose coefficients are returned in r[0..n]. The elements r[nv..n] 
	 * and q[n-nv+1..n] are returned as zero.
	 * @param u u
	 * @param n n
	 * @param v v
	 * @param nv nv
	 * @param q q
	 * @param r r
	 */
	public static void poldiv(double[] u, int n, double[] v, int nv,
			double[] q, double[] r)
	// Given the n+1 coefficients of a polynomial of degree n in u[0..n], and
	// the nv+1 coefficients
	// of another polynomial of degree nv in v[0..nv], divide the polynomial u
	// by the polynomial
	// v (“u”/“v”)giving a quotient polynomial whose coefficients are returned
	// in q[0..n], and a
	// remainder polynomial whose coefficients are returned in r[0..n]. The
	// elements r[nv..n]
	// and q[n-nv+1..n] are returned as zero.
	{
		int k = 0;
		int j = 0;
		for (j = 0; j <= n; j++) {
			r[j] = u[j];
			q[j] = 0.0;
		}
		for (k = n - nv; k >= 0; k--) {
			q[k] = r[nv + k] / v[nv];
			for (j = nv + k - 1; j >= k; j--)
				r[j] -= q[k] * v[j - k];
		}
		for (j = nv; j <= n; j++)
			r[j] = 0.0;
	}

	/**
	 * Given mm, kk, and cof[0..mm+kk], evaluate and return the rational function (cof[0] + 
	 * cof[1]x + · · · + cof[mm]x^mm)/(1 + cof[mm+1]x + · · · + cof[mm+kk]x^kk).
	 * @param x x
	 * @param cof cof
	 * @param mm mm
	 * @param kk kk
	 * @return the result
	 */
	public static double ratval(double x, double[] cof, int mm, int kk)
	// Given mm, kk, and cof[0..mm+kk], evaluate and return the rational
	// function (cof[0] +
	// cof[1]x + · · · + cof[mm]x^mm)/(1 + cof[mm+1]x + · · · + cof[mm+kk]x^kk).
	{
		int j = 0;
		double sumd = 0.0;
		double sumn = 0.0; // Note precision! Change to float if desired.
		for (sumn = cof[mm], j = mm - 1; j >= 0; j--)
			sumn = sumn * x + cof[j];
		for (sumd = 0.0, j = mm + kk; j >= mm + 1; j--)
			sumd = (sumd + cof[j]) * x;
		return sumn / (1.0 + sumd);
	}

	/*
	 * ax2 + bx + c = 0 x = -b ± sqrt(b2 - 4ac)/2a (5.6.2) AND viete=> x =
	 * 2c/(-b ± sqrt(b2 - 4ac))
	 * 
	 * If you use either (5.6.2) or (5.6.3) to get the two roots, you are asking
	 * for trouble: If either a or c (or both) are small, then one of the roots
	 * will involve the subtraction of b from a very nearly equal quantity (the
	 * discriminant); you will get that root very inaccurately. The correct way
	 * to compute the roots is q = -1/2(b + sgn(b)sqrt(b2 - 4ac)) (5.6.4) Then
	 * the two roots are x1 = q/ a and x2 = c/q
	 */

	/**
	 * Solve the quadratic equation ax2 + bx + c = 0.
	 * @param a a
	 * @param b b
	 * @param c c
	 * @return the array containing the solutions.
	 */
	public static double[] quadratic(double a, double b, double c) {
		double[] d = new double[2];
		double delta = Math.pow(b, 2) - 4.0 * a * c;
		// double dbl=(g >= 0.0 ? r : -r);==SIGN(r,g)
		double dbl = (b >= 0.0 ? 1.0 : -1.0);
		double q = 0.0;
		if (a != 0) {
			if (delta >= 0) {
				// d[0] = (-b + Math.pow(delta,0.5))/(2*a);
				// d[1] = (-b - Math.pow(delta,0.5))/(2*a);
				q = -0.5 * (b + dbl * Math.sqrt(delta));
				d[0] = q / a;
				d[1] = c / q;
				failB = false;
			} else
				failB = true;
		} else // bx+c=0
		if (b != 0) {
			d[0] = -c / b;
			d[1] = d[0];
			failB = false;
		} else
			failB = true;

		return d;
	}

	// ax3 +bx2 +cx+d =0-------------->x1,x2,x3 in cel mai bun caz
	// formulele Cardano:
	/*
	 * x3 + ax2 + bx + c = 0 Q = (a^2 - 3b)/ 9 and R = (2a^3 - 9ab + 27c)/54
	 * 
	 * If Q and R are real (always true when a, b, c are real) and R^2 < Q^3,
	 * then the cubic equation has three real roots. Find them by computing
	 * theta = arccos(R/sqrt(Q^3))=>
	 * 
	 * x1 = -2sqrt(Q)cos(theta/3) -a/3 x2 = -2sqrt(Q)cos((theta+2PI)/3) -a/3 x3
	 * = -2sqrt(Q)cos((theta-2PI)/3) -a/3
	 * 
	 * OTHERWISE: If Q and R are both real, equations (5.6.13)–(5.6.14) are
	 * equivalent to A = -sgn(R) [|R| +sqrt(R^2 - Q^3]^1/3 Next compute B = Q/A,
	 * (A != 0) 0 (A = 0)
	 * 
	 * x1 = (A + B) -a/3; the other are complex!!!
	 */
	/**
	 * Solves equation of 3rd order (ax3 +bx2 +cx+d =0) using Cardano formula.
	 * @param a a
	 * @param b b
	 * @param c c
	 * @param d d
	 * @return the array containing the solutions.
	 */
	public static double[] cubic(double a, double b, double c, double d) {
		double[] dd = new double[3];
		/*
		 * double radical3 =1.0/3.0; //prima reducere a ecuatiei<->x3+rx2+sx+t=0
		 * double r=b/a; double s=c/a; double t=d/a; //a doua reducere a
		 * ecuatiei<->x3+px+q=0 double p=s-Math.pow(r,2)/3; double
		 * q=2*Math.pow(r,3)/27-r*s/3+t;
		 * 
		 * double delta = Math.pow(q/2,2)+Math.pow(p/3,3); if (delta>=0)//1
		 * solutie reala { double u1 =
		 * Math.pow(-q/2+Math.sqrt(Math.pow(q/2,2)+Math.pow(p/3,3)),radical3);
		 * double v1 =
		 * Math.pow(-q/2-Math.sqrt(Math.pow(q/2,2)+Math.pow(p/3,3)),radical3);
		 * dd[0]=u1+v1-r/3;//+-r/3 dd[1]=dd[0]; dd[2]=dd[0]; } else//casus
		 * ireductibilus { r=Math.sqrt(-Math.pow(p,3)/27); //unghiul in
		 * intervalul 0-pi; double fi = Math.acos((-q/2)/r);
		 * dd[0]=2*Math.pow(r,radical3)*Math.cos(fi/3)-b/(3*a);
		 * dd[1]=2*Math.pow(r,radical3)*Math.cos(fi/3+2*Math.PI/3)-b/(3*a);
		 * dd[2]=2*Math.pow(r,radical3)*Math.cos(fi/3+4*Math.PI/3)-b/(3*a); }
		 */// ax^3+bx^2+cx+d=0<=>x^3+ax^2+bx+c=0.0;

		double aa = b / a;
		double bb = c / a;
		double cc = d / a;

		double Q = (aa * aa - 3.0 * bb) / 9.0;
		double R = (2.0 * aa * aa * aa - 9.0 * aa * bb + 27.0 * cc) / 54.0;

		if (R * R < Q * Q * Q) {
			double theta = Math.acos(R / Math.sqrt(Q * Q * Q));// 0-PI
			dd[0] = -2.0 * Math.sqrt(Q) * Math.cos(theta / 3.0) - aa / 3.0;
			dd[1] = -2.0 * Math.sqrt(Q)
					* Math.cos((theta + 2.0 * Math.PI) / 3.0) - aa / 3.0;
			dd[2] = -2.0 * Math.sqrt(Q)
					* Math.cos((theta - 2.0 * Math.PI) / 3.0) - aa / 3.0;
		} else {
			double dbl = (R >= 0.0 ? 1.0 : -1.0);
			double A = -dbl
					* Math.pow((Math.abs(R) + Math.sqrt(R * R - Q * Q * Q)),
							1.0 / 3.0);
			double B = 0.0;
			if (A != 0)
				B = Q / A;

			dd[0] = (A + B) - aa / 3.0;
			dd[1] = dd[0];// conjugate imagfinare
			dd[2] = dd[0];// conjugate imagfinare

		}

		return dd;
	}

	// float dfridr(float (*func)(float), float x, float h, float *err)
	public static double errdfridr = 0.0;

	/**
	 * Returns the derivative of a function func at a point x by Ridders’ method of polynomial extrapolation. 
	 * The value h is input as an estimated initial stepsize; it need not be small, but rather should be an increment in x over
	 * which func changes substantially.
	 * @param x s
	 * @param h h
	 * @return the result
	 */
	public double dfridr(double x, double h)// , double *err)
	// Returns the derivative of a function func at a point x by Ridders’ method
	// of polynomial
	// extrapolation. The value h is input as an estimated initial stepsize; it
	// need not be small, but
	// rather should be an increment in x over which func changes substantially.
	// An estimate of the
	// error in the derivative is returned as err.
	{
		int i = 0;
		int j = 0;
		double errt = 0.0;
		double fac = 0.0;
		double hh = 0.0;// ,**a,
		double ans = 0.0;
		double[][] a = new double[NTAB][NTAB];

		failB = false;
		if (h == 0.0) {
			// nrerror("h must be nonzero in dfridr.");
			failS = "h must be nonzero in dfridr.";
			failB = true;

			return 0.0;// Never get here.
		}
		// a=matrix(1,NTAB,1,NTAB);
		hh = h;
		// a[1][1]=((*func)(x+hh)-(*func)(x-hh))/(2.0*hh);
		a[0][0] = (func.F(x + hh) - func.F(x - hh)) / (2.0 * hh);
		errdfridr = BIG;// *err=BIG;
		for (i = 2; i <= NTAB; i++) {
			// Successive columns in the Neville tableau will go to smaller
			// stepsizes and
			// higher orders of extrapolation.
			hh /= CON;
			// a[1][i]=((*func)(x+hh)-(*func)(x-hh))/(2.0*hh); Try new, smaller
			// stepsize.
			a[0][i - 1] = (func.F(x + hh) - func.F(x - hh)) / (2.0 * hh);
			fac = CON2;
			for (j = 2; j <= i; j++) {// Compute extrapolations of various
										// orders, requiring
										// no new function evaluations.
										// a[j][i]=(a[j-1][i]*fac-a[j-1][i-1])/(fac-1.0);
				a[j - 1][i - 1] = (a[j - 2][i - 1] * fac - a[j - 2][i - 2])
						/ (fac - 1.0);
				fac = CON2 * fac;
				// errt=FMAX(fabs(a[j][i]-a[j-1][i]),fabs(a[j][i]-a[j-1][i-1]));
				errt = Math.max(Math.abs(a[j - 1][i - 1] - a[j - 2][i - 1]),
						Math.abs(a[j - 1][i - 1] - a[j - 2][i - 2]));
				// The error strategy is to compare each new extrapolation to
				// one order lower, both
				// at the present stepsize and the previous one.
				if (errt <= errdfridr)// (errt <= *err)
				{ // If error is decreased, save the improved answer.
					errdfridr = errt;// *err=errt;
					ans = a[j - 1][i - 1];// ans=a[j][i];
				}
			}
			// if (fabs(a[i][i]-a[i-1][i-1]) >= SAFE*(*err)) break;
			if (Math.abs(a[i - 1][i - 1] - a[i - 2][i - 2]) >= SAFE
					* (errdfridr))
				break;
			// If higher order is worse by a significant factor SAFE, then quit
			// early.
		}
		// free_matrix(a,1,NTAB,1,NTAB);
		return ans;
	}

	/*
	 * The first of these tasks is straightforward. A generalization of equation
	 * (5.8.7) that is here implemented is to allow the range of approximation
	 * to be between two arbitrary limits a and b, instead of just-1 to 1. This
	 * is effected by a change of variable y =[x - 1/2*(b + a)]/[1/2*(b-a)]
	 * (5.8.10) and by the approximation of f(x) by a Chebyshev polynomial in y.
	 * Tn(x) = cos(n arccos x).
	 */
	/**
	 * Chebyshev fit: Given a function func, lower and upper limits of the interval [a,b], and a 
	 * maximum degree n, this routine computes the n coefficients c[0..n-1] such that func(x) ~ 
	 * [sum from k=0 to n-1 of ckTk(y)]-c0/2. Here, y = (x-0.5(b+a))/(0.5(b-a)) and tk(y)=cos(k arccos(y))
	 * This routine is to be used with moderately large n (e.g., 30 or 50), the array of c’s subsequently to be 
	 * truncated at the smaller value m such that cm and subsequent elements are negligible.
	 * @param a a
	 * @param b b
	 * @param c c
	 * @param n n
	 */
	public void chebft(double a, double b, double[] c, int n)// , float
																// (*func)(float))
	// Chebyshev fit: Given a function func, lower and upper limits of the
	// interval [a,b], and a
	// maximum degree n, this routine computes the n coefficients c[0..n-1] such
	// that func(x) ~
	// [sum from k=0 to n-1 of ckTk(y)]-c0/2, where y and x are related by
	// (5.8.10). This routine is to be used with
	// moderately large n (e.g., 30 or 50), the array of c’s subsequently to be
	// truncated at the smaller
	// value m such that cm and subsequent elements are negligible.
	{
		int k = 0;
		int j = 0;
		double fac = 0.0;
		double bpa = 0.0;
		double bma = 0.0;
		double[] f = new double[n];
		// f=vector(0,n-1);
		bma = 0.5 * (b - a);
		bpa = 0.5 * (b + a);
		for (k = 0; k < n; k++) {// We evaluate the function at the n points
									// required by (5.8.7).
			double y = Math.cos(PI * (k + 0.5) / n);
			f[k] = func.F(y * bma + bpa);// f[k]=(*func)(y*bma+bpa);
		}
		fac = 2.0 / n;
		for (j = 0; j < n; j++) {
			double sum = 0.0; // We will accumulate the sum in double precision,
								// a nicety that you can ignore.
			for (k = 0; k < n; k++)
				sum += f[k] * Math.cos(PI * j * (k + 0.5) / n);
			c[j] = fac * sum;
		}
		// free_vector(f,0,n-1);
	}

	/*
	 * Now that we have the Chebyshev coefficients, how do we evaluate the
	 * approximation? One could use the recurrence relation of equation (5.8.2)
	 * to generate values for Tk(x) from T0 = 1, T1 = x, while also accumulating
	 * the sum of (5.8.9). It is better to use Clenshaw’s recurrence formula
	 * (§5.5), effecting the two processes simultaneously.
	 */
	/**
	 * Chebyshev evaluation: c[0..m-1] is an array of Chebyshev coefficients, the first m elements of c output from chebft 
	 * (which must have been called with the same a and b). The Chebyshev polynomial sum from k=0 to m-1 of [ckTk(y)] - c0/2 
	 * is evaluated at a point y = [x - (b + a)/2]/[(b - a)/2], and the result is returned as the function value.
	 * @param a a
	 * @param b b
	 * @param c c
	 * @param m m
	 * @param x x
	 * @return the result
	 */
	public double chebev(double a, double b, double[] c, int m, double x)
	// Chebyshev evaluation: All arguments are input. c[0..m-1] is an array of
	// Chebyshev coeffi-
	// cients, the first m elements of c output from chebft (which must have
	// been called with the
	// same a and b). The Chebyshev polynomial
	// sum from k=0 to m-1 of [ckTk(y)] - c0/2 is evaluated at a point
	// y = [x - (b + a)/2]/[(b - a)/2], and the result is returned as the
	// function value.
	{
		// void nrerror(char error_text[]);
		double d = 0.0;
		double dd = 0.0;
		double sv = 0.0;
		double y = 0.0;
		double y2 = 0.0;
		int j = 0;

		failB = false;
		if ((x - a) * (x - b) > 0.0) {
			// nrerror("x not in range in routine chebev");
			failS = "x not in range in routine chebev";
			failB = true;

			return 0.0;// Never get here.
		}
		y2 = 2.0 * (y = (2.0 * x - a - b) / (b - a)); // Change of variable.
		for (j = m - 1; j >= 1; j--) {// Clenshaw’s recurrence.
			sv = d;
			d = y2 * d - dd + c[j];
			dd = sv;
		}
		return y * d - dd + 0.5 * c[0]; // Last step is different.
	}

	/*
	 * If you have obtained the Chebyshev coefficients that approximate a
	 * function in a certain range (e.g., from chebft in §5.8), then it is a
	 * simple matter to transform them to Chebyshev coefficients corresponding
	 * to the derivative or integral of the function. Having done this, you can
	 * evaluate the derivative or integral just as if it were a function that
	 * you had Chebyshev-fitted ab initio.
	 */
	/**
	 * Given a,b,c[0..n-1], as output from routine chebft and given n, the desired degree 
	 * of approximation (length of c to be used), this routine returns the array cder[0..n-1], the
	 * Chebyshev coefficients of the derivative of the function whose coefficients are c.
	 * @param a a
	 * @param b b
	 * @param c c
	 * @param cder cder
	 * @param n n
	 */
	public void chder(double a, double b, double[] c, double[] cder, int n)
	// Given a,b,c[0..n-1], as output from routine chebft §5.8, and given n, the
	// desired degree
	// of approximation (length of c to be used), this routine returns the array
	// cder[0..n-1], the
	// Chebyshev coefficients of the derivative of the function whose
	// coefficients are c.
	{
		int j = 0;
		double con = 0.0;
		cder[n - 1] = 0.0; // n-1 and n-2 are special cases.
		cder[n - 2] = 2 * (n - 1) * c[n - 1];
		for (j = n - 3; j >= 0; j--)
			cder[j] = cder[j + 2] + 2 * (j + 1) * c[j + 1];// Equation (5.9.2).
		con = 2.0 / (b - a);
		for (j = 0; j < n; j++)
			// Normalize to the interval b-a.
			cder[j] *= con;
	}

	/*
	 * integral f(x) = Const+g(x) g- is a function of x. Const+g(a)=0.0;=>
	 * Const. Then, Intgral=Const+g(xinput);!!!!
	 */
	/**
	 * Given a,b,c[0..n-1], as output from routine chebft and given n, the desired degree 
	 * of approximation (length of c to be used), this routine returns the array cint[0..n-1], the 
	 * Chebyshev coefficients of the integral of the function whose coefficients are c. The constant of 
	 * integration is set so that the integral vanishes at a.
	 * @param a a
	 * @param b b
	 * @param c c
	 * @param cint cint
	 * @param n n
	 */
	public void chint(double a, double b, double[] c, double[] cint, int n)
	// Given a,b,c[0..n-1], as output from routine chebft §5.8, and given n, the
	// desired degree
	// of approximation (length of c to be used), this routine returns the array
	// cint[0..n-1], the
	// Chebyshev coefficients of the integral of the function whose coefficients
	// are c. The constant of
	// integration is set so that the integral vanishes at a.
	{
		int j = 0;
		double sum = 0.0;
		double fac = 1.0;
		double con = 0.0;
		con = 0.25 * (b - a);// Factor that normalizes to the interval b-a.
		for (j = 1; j <= n - 2; j++) {
			cint[j] = con * (c[j - 1] - c[j + 1]) / j; // Equation (5.9.1).
			sum += fac * cint[j]; // Accumulates the constant of integration.
			fac = -fac; // Will equal ±1.
		}
		cint[n - 1] = con * c[n - 2] / (n - 1);// Special case of (5.9.1) for
												// n-1.
		sum += fac * cint[n - 1];
		cint[0] = 2.0 * sum; // Set the constant of integration.
	}

	/*
	 * You may well ask after reading the preceding two sections, “Must I store
	 * and evaluate my Chebyshev approximation as an array of Chebyshev
	 * coefficients for a transformed variable y? Can’t I convert the ck’s into
	 * actual polynomial coefficients in the original variable x and have an
	 * approximation of the following form?” f(x) = sum from k=0 to m-1 of gkxk
	 */
	/**
	 * Chebyshev polynomial coefficients. Given a coefficient array c[0..n-1], this routine generates 
	 * a coefficient array d[0..n-1] such that sum from k=0 to n-1 of dky^k = sum from k=0 to n-1 of ckTk(y) - c0/2. The method 
	 * is Clenshaw’s recurrence.
	 * @param c c
	 * @param d d
	 * @param n n
	 */
	public void chebpc(double[] c, double[] d, int n)
	// Chebyshev polynomial coefficients. Given a coefficient array c[0..n-1],
	// this routine generates
	// a coefficient array d[0..n-1] such that sum from k=0 to n-1 of dky^k =
	// sum from k=0 to n-1
	// of ckTk(y) - c0/2. The method
	// is Clenshaw’s recurrence (5.8.11), but now applied algebraically rather
	// than arithmetically.
	{
		int k = 0;
		int j = 0;
		double sv = 0.0;
		double[] dd = new double[n];
		// dd=vector(0,n-1);
		for (j = 0; j < n; j++)
			d[j] = dd[j] = 0.0;
		d[0] = c[n - 1];
		for (j = n - 2; j >= 1; j--) {
			for (k = n - j; k >= 1; k--) {
				sv = d[k];
				d[k] = 2.0 * d[k - 1] - dd[k];
				dd[k] = sv;
			}
			sv = d[0];
			d[0] = -dd[0] + c[j];
			dd[0] = sv;
		}
		for (j = n - 1; j >= 1; j--)
			d[j] = d[j - 1] - dd[j];
		d[0] = -dd[0] + 0.5 * c[0];
		// free_vector(dd,0,n-1);
	}

	/**
	 * Polynomial coefficient shift. Given a coefficient array d[0..n-1], this routine generates a 
	 * coefficient array g[0..n-1] such that sum from k=0 to n-1 of dky^k = sum from k=0 to n-1 gkx^k, where x and y are related 
	 * by y = [x - (b + a)/2]/[(b - a)/2], i.e., the interval -1 less y less 1 is mapped to the interval a less x less b. The array 
	 * g is returned in d.
	 * @param a a
	 * @param b b
	 * @param d d
	 * @param n n
	 */
	public void pcshft(double a, double b, double[] d, int n)
	// Polynomial coefficient shift. Given a coefficient array d[0..n-1], this
	// routine generates a
	// coefficient array g[0..n-1] such that sum from k=0 to n-1 of dky^k = sum
	// from k=0 to n-1
	// gkx^k, where x and y are related
	// by (5.8.10), i.e., the interval -1 < y < 1 is mapped to the interval a <
	// x < b. The array
	// g is returned in d.
	{
		int k = 0;
		int j = 0;
		double fac = 0.0;
		double cnst = 0.0;
		cnst = 2.0 / (b - a);
		fac = cnst;
		for (j = 1; j < n; j++) {// First we rescale by the factor const...
			d[j] *= fac;
			fac *= cnst;
		}
		cnst = 0.5 * (a + b); // ...which is then redefined as the desired
								// shift.
		for (j = 0; j <= n - 2; j++)
			// We accomplish the shift by synthetic division. Synthetic
			// division is a miracle of high-school algebra. If you
			// never learned it, go do so. You won’t be sorry.
			for (k = n - 2; k >= j; k--)
				d[k] -= cnst * d[k + 1];
	}

	/**
	 * Inverse of routine chebpc: given an array of polynomial coefficients d[0..n-1], returns an 
	 * equivalent array of Chebyshev coefficients c[0..n-1].
	 * @param d d
	 * @param c c
	 * @param n n
	 */
	public void pccheb(double[] d, double[] c, int n)
	// Inverse of routine chebpc: given an array of polynomial coefficients
	// d[0..n-1], returns an
	// equivalent array of Chebyshev coefficients c[0..n-1].
	{
		int j = 0;
		int jm = 0;
		int jp = 0;
		int k = 0;
		double fac = 0.0;
		double pow = 0.0;
		pow = 1.0; // Will be powers of 2.
		c[0] = 2.0 * d[0];
		for (k = 1; k < n; k++) {// Loopov er orders of x in the polynomial.
			c[k] = 0.0; // Zero corresponding order of Chebyshev.
			fac = d[k] / pow;
			jm = k;
			jp = 1;
			for (j = k; j >= 0; j -= 2, jm--, jp++) {
				// Increment this and lower orders of Chebyshev with the
				// combinatorial coefficent times
				// d[k]; see text for formula.
				c[j] += fac;
				fac *= ((float) jm) / ((float) jp);
			}
			pow += pow;
		}
	}

	/*
	 * The fourth and fifth steps are accomplished by the routines chebpc and
	 * pcshft, respectively. Here is how the procedure looks all together:
	 * #define NFEW .. #define NMANY .. float *c,*d,*e,a,b; Economize NMANY
	 * power series coefficients e[0..NMANY-1] in the range (a, b) into NFEW
	 * coefficients d[0..NFEW-1]. c=vector(0,NMANY-1); d=vector(0,NFEW-1);
	 * e=vector(0,NMANY-1); pcshft((-2.0-b-a)/(b-a),(2.0-b-a)/(b-a),e,NMANY);
	 * pccheb(e,c,NMANY); ... Here one would normally examine the Chebyshev
	 * coefficients c[0..NMANY-1] to decide how small NFEW can be.
	 * chebpc(c,d,NFEW); pcshft(a,b,d,NFEW);
	 */

	// @@@@@@Pad´e Approximants of a infinite series f(x)=sum ckx^k
	public static double resid = 0.0;

	/**
	 * Given cof[0..2*n], the leading terms in the power series expansion of a function, solve the 
	 * linear Pad´e equations to return the coefficients of a diagonal rational function approximation to 
	 * the same function, namely (cof[0] + cof[1]x+· · ·+ cof[n]x^N)/(1 + cof[n+1]x +· · ·+cof[2*n]x^N). 
	 * The value resid is the norm of the residual vector; a small value indicates a well-converged solution.
	 * @param cof cof
	 * @param n n
	 */
	public static void pade(double[] cof, int n)// double[] *resid)
	// Given cof[0..2*n], the leading terms in the power series expansion of a
	// function, solve the
	// linear Pad´e equations to return the coefficients of a diagonal rational
	// function approximation to
	// the same function, namely (cof[0] + cof[1]x+· · ·+ cof[n]x^N)/(1 +
	// cof[n+1]x +· · ·+
	// cof[2*n]x^N). The value resid is the norm of the residual vector; a small
	// value indicates a
	// well-converged solution. Note that cof is double precision for
	// consistency with ratval.
	{
		// void lubksb(float **a, int n, int *indx, float b[]);
		// void ludcmp(float **a, int n, int *indx, float *d);
		// void mprove(float **a, float **alud, int n, int indx[], float b[],
		// float x[]);
		int j = 0;
		int k = 0;
		int[] indx = new int[n];
		//double d = 0.0;
		double rr = 0.0;
		double rrold = 0.0;
		double sum = 0.0;
		double[][] q = new double[n][n];
		double[][] qlu = new double[n][n];
		double[] x = new double[n];
		double[] y = new double[n];
		double[] z = new double[n];
		// indx=ivector(1,n);
		// q=matrix(1,n,1,n);
		// qlu=matrix(1,n,1,n);
		// x=vector(1,n);
		// y=vector(1,n);
		// z=vector(1,n);
		for (j = 1; j <= n; j++) {// Set up matrix for solving.Given cof[0..2*n]
			y[j - 1] = x[j - 1] = cof[n + j];// y[j]=x[j]=cof[n+j];
			for (k = 1; k <= n; k++) {
				q[j - 1][k - 1] = cof[j - k + n];// q[j][k]=cof[j-k+n];
				qlu[j - 1][k - 1] = q[j - 1][k - 1];// qlu[j][k]=q[j][k];
			}
		}
		// =>dlu
		LinAEq.ludcmp(qlu, n, indx);// ,&d); //Solve by LU decomposition and
									// backsubstitution.
		LinAEq.lubksb(qlu, n, indx, x);
		rr = BIG;
		do {// Important to use iterative improvement, since
			// the Pad´e equations tend to be ill-conditioned.
			rrold = rr;
			for (j = 1; j <= n; j++)
				z[j - 1] = x[j - 1];// z[j]=x[j];
			LinAEq.mprove(q, qlu, n, indx, y, x);
			for (rr = 0.0, j = 1; j <= n; j++)
				// Calculate residual.
				rr += (z[j - 1] - x[j - 1]) * (z[j - 1] - x[j - 1]);// rr +=
																	// SQR(z[j]-x[j]);
		} while (rr < rrold); // If it is no longer improving, call it quits.
		resid = Math.sqrt(rrold);
		for (k = 1; k <= n; k++) {// Calculate the remaining coefficients.
			for (sum = cof[k], j = 1; j <= k; j++)
				sum -= z[j - 1] * cof[k - j];// sum -= z[j]*cof[k-j];
			y[k - 1] = sum;// y[k]=sum;
		}// Copy answers to output.
		for (j = 1; j <= n; j++) {
			cof[j] = y[j - 1];// cof[j]=y[j];
			cof[j + n] = -z[j - 1];// cof[j+n] = -z[j];
		}
		// free_vector(z,1,n);
		// free_vector(y,1,n);
		// free_vector(x,1,n);
		// free_matrix(qlu,1,n,1,n);
		// free_matrix(q,1,n,1,n);
		// free_ivector(indx,1,n);
	}

	// @@@Rational Chebyshev Approximation
	public static double dev = 0.0;

	// public static void ratlsq(double (*fn)(double), double a, double b, int
	// mm, int kk,
	/**
	 * Returns in cof[0..mm+kk] the coefficients of a rational function approximation to the function 
	 * fn in the interval (a, b). Input quantities mm and kk specify the order of the numerator and 
	 * denominator, respectively. The maximum absolute deviation of the approximation (insofar as is known) is returned as dev.
	 * @param a a
	 * @param b b
	 * @param mm mm
	 * @param kk kk
	 * @param cof cof
	 */
	public void ratlsq(double a, double b, int mm, int kk, double[] cof)// ,
																		// double
																		// *dev)
	// Returns in cof[0..mm+kk] the coefficients of a rational function
	// approximation to the function
	// fn in the interval (a, b). Input quantities mm and kk specify the order
	// of the numerator and
	// denominator, respectively. The maximum absolute deviation of the
	// approximation (insofar as
	// is known) is returned as dev.
	{
		// double ratval(double x, double cof[], int mm, int kk);
		// void dsvbksb(double **u, double w[], double **v, int m, int n, double
		// b[],
		// double x[]);
		// void dsvdcmp(double **a, int m, int n, double w[], double **v);
		// These are double versions of svdcmp, svbksb.
		int i = 0;
		int it = 0;
		int j = 0;
		int ncof = 0;
		int npt = 0;
		double devmax = 0.0;
		double e = 0.0;
		double hth = 0.0;
		double power = 0.0;
		double sum = 0.0;
		ncof = mm + kk + 1;
		npt = NPFAC * ncof;// Number of points where function is evaluated,
		// i.e., fineness of the mesh.
		double[] bb = new double[npt];
		double[] coff = new double[ncof];// 0-ncof-1!!!!
		double[] ee = new double[npt];
		double[] fs = new double[npt];
		double[][] u = new double[npt][ncof];
		double[][] v = new double[ncof][ncof];
		double[] w = new double[ncof];
		double[] wt = new double[npt];
		double[] xs = new double[npt];

		// bb=dvector(1,npt);
		// coff=dvector(0,ncof-1);
		// ee=dvector(1,npt);
		// fs=dvector(1,npt);
		// u=dmatrix(1,npt,1,ncof);
		// v=dmatrix(1,ncof,1,ncof);
		// w=dvector(1,ncof);
		// wt=dvector(1,npt);
		// xs=dvector(1,npt);

		dev = BIG;
		for (i = 1; i <= npt; i++) {// Fill arrays with mesh abscissas and
									// function values.
			if (i < npt / 2) {
				hth = PIO2 * (i - 1) / (npt - 1.0);// At each end, use formula
													// that minimizes roundoff
													// sensitivity.
				xs[i - 1] = a + (b - a) * Math.sin(hth) * Math.sin(hth);// xs[i]=a+(b-a)*DSQR(sin(hth));
			} else {
				hth = PIO2 * (npt - i) / (npt - 1.0);
				xs[i - 1] = b - (b - a) * Math.sin(hth) * Math.sin(hth);// xs[i]=b-(b-a)*DSQR(sin(hth));
			}
			fs[i - 1] = func.F(xs[i - 1]);// fs[i]=(*fn)(xs[i]);
			wt[i - 1] = 1.0; // In later iterations we will adjust these weights
								// to
								// combat the largest deviations.
			ee[i - 1] = 1.0;
		}
		e = 0.0;
		for (it = 1; it <= MAXIT; it++) {// Loop over iterations.
			for (i = 1; i <= npt; i++) {// Set up the “design matrix” for the
										// least-squares fit.
				power = wt[i - 1];// power=wt[i];
				// double dbl=(g >= 0.0 ? r : -r);==SIGN(r,g)
				double dbl = (ee[i - 1] >= 0.0 ? e : -e);
				// bb[i]=power*(fs[i]+SIGN(e,ee[i]));
				bb[i - 1] = power * (fs[i - 1] + dbl);
				// Key idea here: Fit to fn(x)+e where the deviation is
				// positive, to fn(x)-e where
				// it is negative. Then e is supposed to become an approximation
				// to the equal-ripple
				// deviation.
				for (j = 1; j <= mm + 1; j++) {
					u[i - 1][j - 1] = power;// u[i][j]=power;
					power *= xs[i - 1];// power *= xs[i];
				}
				power = -bb[i - 1];// power = -bb[i];
				for (j = mm + 2; j <= ncof; j++) {
					power *= xs[i - 1];// power *= xs[i];
					u[i - 1][j - 1] = power;// u[i][j]=power;
				}
			}
			// dsvdcmp(u,npt,ncof,w,v); Singular Value Decomposition.
			LinAEq.svdcmp(u, npt, ncof, w, v);
			// In especially singular or difficult cases, one might here edit
			// the singular values w[1..ncof],
			// replacing small values by zero. Note that dsvbksb works with
			// one-based arrays, so we
			// must subtract 1 when we pass it the zero-based array coff.
			// dsvbksb(u,w,v,npt,ncof,bb,coff-1);
			LinAEq.svbksb(u, w, v, npt, ncof, bb, coff);// -1);
			devmax = sum = 0.0;
			for (j = 1; j <= npt; j++) {// Tabulate the deviations and revise
										// the weights.
										// ee[j]=ratval(xs[j],coff,mm,kk)-fs[j];
				ee[j - 1] = ratval(xs[j - 1], coff, mm, kk) - fs[j - 1];
				wt[j - 1] = Math.abs(ee[j - 1]);// wt[j]=fabs(ee[j]);// Use
												// weighting to emphasize most
												// deviant points.
				sum += wt[j - 1];// sum += wt[j];
				// if (wt[j] > devmax) devmax=wt[j];
				if (wt[j - 1] > devmax)
					devmax = wt[j - 1];
			}
			e = sum / npt;// Update e to be the mean absolute deviation.
			if (devmax <= dev) {// Save only the best coefficient set found.
				for (j = 0; j < ncof; j++)
					cof[j] = coff[j];
				dev = devmax;
			}
			// printf(" ratlsq iteration= %2d max error= %10.3e\n",it,devmax);
		}
		// free_dvector(xs,1,npt);
		// free_dvector(wt,1,npt);
		// free_dvector(w,1,ncof);
		// free_dmatrix(v,1,ncof,1,ncof);
		// free_dmatrix(u,1,npt,1,ncof);
		// free_dvector(fs,1,npt);
		// free_dvector(ee,1,npt);
		// free_dvector(coff,0,ncof-1);
		// free_dvector(bb,1,npt);
	}
}
