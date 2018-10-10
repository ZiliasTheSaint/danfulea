package danfulea.math.numerical;

//import java.text.DecimalFormat;
//import java.text.DecimalFormatSymbols;
import java.text.NumberFormat;
import java.util.Locale;



/**
 * Class for solving integrals of functions.
 * Based on literature, including Numerical Recipes in C (Cambridge Univ.).
 * @author Dan Fulea, 02 OCT. 2006
 */
public class Integrator {
	public static boolean failB = false;
	public static String failS = "";

	private Function func;
	private NumberFormat nf = NumberFormat.getInstance(Locale.US);
	// private String pattern="0.###E0";
	// private DecimalFormatSymbols dfs=new DecimalFormatSymbols(Locale.US);
	// private DecimalFormat nff = new DecimalFormat(pattern,dfs);
	private int idigits = 3;

	public static double strap = 0.0;
	public static double smidpnt = 0.0;
	public static double smidinf = 0.0;
	public static double smidsql = 0.0;
	public static double smidsqu = 0.0;
	public static double smidexp = 0.0;
	public static double EPS = 1.0e-5;
	public static int JMAX = 100;
	public static int JMAXP = JMAX + 1;
	public static int K = 5;
	public static double EPS1 = 1.0e-6;
	public static int JMAX1 = 14;
	public static int JMAXP1 = JMAX1 + 1;
	public static double EPS2 = 3.0e-11;// EPS is the relative precision.
	public static double EPS3 = 3.0e-14;// Increase EPS if you don’t have this
										// precision.
	public static int MAXIT = 10;
	public static double PIM4 = 0.7511255444649425;// 1/PI^1/4.

	public static double xsav = 0.0;
	public static double ysav = 0.0;

	// ---------DCADRE--------------------------------------------------
	private double DIFFq = 0.;
	private boolean H2CONVq = false;
	private boolean AITKENq = false;
	private int LM1q = 0;
	private int N2q = 0;
	private double FNq = 0.;
	private int ISTEPq = 0;
	private int IIq = 0;
	private int IIIq = 0;
	private double HOVNq = 0.;
	private double FIq = 0.;
	private int ISTEP2q = 0;
	private double SUMq = 0.;
	private double SUMABSq = 0.;
	@SuppressWarnings("unused")
	private double ABSIq = 0.;
	private int ITq = 0;
	private double TABTLMq = 0.;
	private double ERGLq = 0.;
	private double ERGOALq = 0.;
	private double FEXTRPq = 0.;
	private double ERRERq = 0.;
	private boolean RIGHTq = false;
	private double STEPq = 0.;
	private double ASTEPq = 0.;
	private double TABSq = 0.;
	private int Lq = 0;
	private int Nq = 0;
	private double[] TSq = new double[2049];
	private double[][] Tq = new double[10][10];
	private double ZEROq = 0.0;
	private double P1q = 0.1;
	private double HALFq = 0.5;
	private double ONEq = 1.0;
	private double TWOq = 2.0;
	private double FOURq = 4.0;
	private double FOURP5q = 4.5;
	private double TENq = 10.0;
	private double HUNq = 100.0;
	private double AITLOWq = 1.1;
	private double H2TOLq = 0.15;
	private double AITTOLq = 0.1;
	private double JUMPTLq = 0.01;
	private int MAXTSq = 2049;
	private int MAXTBLq = 10;
	private int MXSTGEq = 30;
	private double SLOPEq = 0.;
	private double FBEG2q = 0.;
	private double[] RNq = new double[4];
	private boolean[] REGLSVq = new boolean[30];
	private double[] BEGINq = new double[30];
	private double[] FINISq = new double[30];
	private double[] ESTq = new double[30];
	private int[] IBEGSq = new int[30];
	private int NNLEFTq = 0;
	private double ALG4O2q = 0.;
	private double CADREq = 0.;
	private double CURESTq = 0.;
	private double VINTq = 0.;
	private double DCADREq = 0.;
	private double LENGTHq = 0.;
	private double ERRRq = 0.;
	private double ERRAq = 0.;
	private double STEPMNq = 0.;
	private double STEPNMq = 0.;
	private double STAGEq = 0.;
	private int ISTAGEq = 0;
	private double FNSIZEq = 0.;
	private double PREVERq = 0.;
	private boolean REGLARq = false;
	private double BEGq = 0.;
	private double RVALq = 0.;
	private double FBEGq = 0.;
	private int IBEGq = 0;
	private double ENDq = 0.;
	private double FENDq = 0.;
	private int IENDq = 0;
	private double Aq = 0.;
	private double Bq = 0.;
	private double AERRq = 0.;
	private double RERRq = 0.;
	private double ERRORq = 0.;
	private int IERq = 0;
	private int Iq = 0;
	private double SINGq = 0.;
	private double FEXTM1q = 0.;
	private double[] Rq = new double[10];
	private double[] AITq = new double[10];
	private double[] DIFq = new double[10];
	@SuppressWarnings("unused")
	private double ALPHAq = 0.;
	private double H2NXTq = 0.;
	private double SINGNXq = 0.;
	private double ERRETq = 0.;
	private double H2TFEXq = 0.;
	// ----------------------------------------------
	private boolean exitq = false;
	public boolean first_time = true;
	public double $RERR = 1.E-5; // "RERR-VALUE NEEDED BY DCADRE"
	public double $AERR = 1.E-16; // "AERR-VALUE NEEDED BY DCADRE"
	public String funcname = "";

	/**
	 * Constructor. Setup function and some defaults.
	 * @param func the class implementing Function interface.
	 */
	public Integrator(Function func) {
		this.func = func;
		set_defaults();
	}

	/**
	 * Set the numbers to be displayed with i digits.
	 * @param i the number of digits
	 */
	public void set_defaults(int i) {
		idigits = i;
		nf.setMinimumFractionDigits(idigits);// 3);//default is 2!!
		nf.setMaximumFractionDigits(idigits);// (3);//default is 2!!
		nf.setGroupingUsed(false);// no 4,568.02 but 4568.02
	}

	/**
	 * Set the numbers to be displayed with 3 digits.
	 */
	public void set_defaults() {
		idigits = 3;
		nf.setMinimumFractionDigits(idigits);// 3);//default is 2!!
		nf.setMaximumFractionDigits(idigits);// (3);//default is 2!!
		nf.setGroupingUsed(false);// no 4,568.02 but 4568.02
	}

	/**
	 * This routine computes the nth stage of refinement of an extended trapezoidal rule. func is 
	 * to be integrated between limits a and b. When called with n=1, the routine returns the crudest estimate of integral from a to b 
	 * f(x)dx. Subsequent calls with n=2,3,..(in that sequential order) will improve the accuracy by adding 2^n-2 
	 * additional interior points.
	 * @param a a
	 * @param b b
	 * @param n n
	 * @return the result
	 */
	public double trapzd(double a, double b, int n)// float (*func)(float),
													// float a, float b, int n)
	// This routine computes the nth stage of refinement of an extended
	// trapezoidal rule. func is input
	// as a pointer to the function to be integrated between limits a and b,
	// also input. When called with
	// n=1, the routine returns the crudest estimate of integral from a to b
	// f(x)dx. Subsequent calls with n=2,3,...
	// (in that sequential order) will improve the accuracy by adding 2^n-2
	// additional interior points.
	{
		double x = 0.0;
		double tnm = 0.0;
		double sum = 0.0;
		double del = 0.0;
		// double s=0.0;
		int it = 0;
		int j = 0;

		if (n == 1) {
			// return
			// (s=0.5*(b-a)*(func.F(a)+func.F(b)));//(s=0.5*(b-a)*(FUNC(a)+FUNC(b)));
			return (strap = 0.5 * (b - a) * (func.F(a) + func.F(b)));
		} else {
			for (it = 1, j = 1; j < n - 1; j++) {
				it <<= 1;// System.out.println("it= "+it);=>2,4,8,16,32....
			}
			tnm = it;
			del = (b - a) / tnm; // This is the spacing of the points to be
									// added.
			x = a + 0.5 * del;
			for (sum = 0.0, j = 1; j <= it; j++, x += del)
				sum += func.F(x);// sum += FUNC(x);
			// s=0.5*(s+(b-a)*sum/tnm); //This replaces s by its refined value.
			strap = 0.5 * (strap + (b - a) * sum / tnm);// s=(s+(b-a)*sum/tnm);
			return strap;
		}
		/*
		 * double eps=1e-8; int in=1; //evaluare optimizata a pasului de
		 * integrare double h=b-a; double s=0.5*h*(func.F(a)+func.F(b)); double
		 * del=0.0; boolean zero=true;//control while(zero) { del=0.0; for (int
		 * i=1; i<=in; i++) del=del+func.F(a+(i-0.5)*h); del=0.5*(h*del-s);
		 * s=s+del; h=0.5*h; in=2*in;//pana la 2^30!! maxvalue=2^31-1!!!
		 * if(in>n) zero=false; } return s;
		 */
	}

	/*
	 * Much better, of course, is to refine the trapezoidal rule until some
	 * specified degree of accuracy has been achieved:
	 */
	// float qtrap(float (*func)(float), float a, float b)
	/**
	 * Returns the integral of the function func from a to b. The parameters EPS can be set to the 
	 * desired fractional accuracy and JMAX so that 2 to the power JMAX-1 is the maximum allowed 
	 * number of steps. Integration is performed by the trapezoidal rule.
	 * @param a a
	 * @param b b
	 * @return the result
	 */
	public double qtrap(double a, double b)
	// Returns the integral of the function func from a to b. The parameters EPS
	// can be set to the
	// desired fractional accuracy and JMAX so that 2 to the power JMAX-1 is the
	// maximum allowed
	// number of steps. Integration is performed by the trapezoidal rule.
	{
		// float trapzd(float (*func)(float), float a, float b, int n);
		// void nrerror(char error_text[]);
		int j = 0;
		double s = 0.0;
		double olds = 0.0; // Initial value of olds is arbitrary.

		failB = false;
		for (j = 1; j <= JMAX; j++) {
			// s=trapzd(func,a,b,j);
			s = trapzd(a, b, j);
			if (j > 5) // Avoid spurious early convergence.
				// if (fabs(s-olds) < EPS*fabs(olds) ||
				if (Math.abs(s - olds) < EPS * Math.abs(olds)
						|| (s == 0.0 && olds == 0.0))
					return s;
			olds = s;
		}
		// nrerror("Too many steps in routine qtrap");
		failS = "Too many steps in routine qtrap";
		failB = true;

		return 0.0;// Never get here.
	}

	// float qsimp(float (*func)(float), float a, float b)
	/**
	 * Returns the integral of the function func from a to b. The parameters EPS can be set to the 
	 * desired fractional accuracy and JMAX so that 2 to the power JMAX-1 is the maximum allowed 
	 * number of steps. Integration is performed by Simpson’s rule.
	 * @param a a
	 * @param b b
	 * @return the result
	 */
	public double qsimp(double a, double b)
	// Returns the integral of the function func from a to b. The parameters EPS
	// can be set to the
	// desired fractional accuracy and JMAX so that 2 to the power JMAX-1 is the
	// maximum allowed
	// number of steps. Integration is performed by Simpson’s rule.
	{
		// float trapzd(float (*func)(float), float a, float b, int n);
		// void nrerror(char error_text[]);
		int j = 0;
		double s = 0.0;
		double st = 0.0;
		double ost = 0.0;
		double os = 0.0;

		failB = false;
		for (j = 1; j <= JMAX; j++) {
			st = trapzd(a, b, j);// trapzd(func,a,b,j);
			s = (4.0 * st - ost) / 3.0; // Compare equation (4.2.4), above.
			if (j > 5) // Avoid spurious early convergence.
				// if (fabs(s-os) < EPS*fabs(os) ||
				if (Math.abs(s - os) < EPS * Math.abs(os)
						|| (s == 0.0 && os == 0.0))
					return s;
			os = s;
			ost = st;
		}

		failS = "Too many steps in routine qsimp";
		failB = true;

		// nrerror("Too many steps in routine qsimp");
		return 0.0; // Never get here.
	}

	// public double qromb(float (*func)(float), float a, float b)
	/*
	 * Quite often you will want to call polint with the dummy arguments xa and
	 * ya replaced by actual arrays with offsets. For example, the construction
	 * polint(&xx[14],&yy[14],4,x,y,dy) performs 4-point interpolation on the
	 * tabulated values xx[15..18], yy[15..18]. For more on this, see the end of
	 * §3.4.
	 */
	/**
	 * Returns the integral of the function func from a to b. Integration is performed by Romberg’s 
	 * method of order 2K, where, e.g., K=2 is Simpson’s rule.
	 * @param a a
	 * @param b b
	 * @return the result
	 */
	public double qromb(double a, double b)
	// Returns the integral of the function func from a to b. Integration is
	// performed by Romberg’s
	// method of order 2K, where, e.g., K=2 is Simpson’s rule.
	{
		// void polint(float xa[], float ya[], int n, float x, float *y, float
		// *dy);
		// float trapzd(float (*func)(float), float a, float b, int n);
		// void nrerror(char error_text[]);
		double ss = 0.0;
		double dss = 0.0;
		double[] s = new double[JMAXP];
		double[] h = new double[JMAXP + 1];
		// These store the successive trapezoidal approximations
		// and their relative stepsizes.
		failB = false;
		int j = 0;
		h[0] = 1.0;// h[1]=1.0;
		for (j = 1; j <= JMAX; j++) {
			// s[j]=trapzd(func,a,b,j);
			s[j - 1] = trapzd(a, b, j);
			// if (j >= K)
			if (j >= K) {
				// polint(&h[j-K],&s[j-K],K,0.0,&ss,&dss);
				// &h[j-K]= pointer to the j-K element of vector=> it was taken
				// the elements from j-K to the end!!!
				double[] hh = new double[K];// [j-K];
				double[] sss = new double[K];// [j-K];
				for (int k = 1; k <= K; k++)// for (int k=1;k<=j-K;k++)
				{
					hh[k - 1] = h[j - K + k - 1];// k-1];
					sss[k - 1] = s[j - K + k - 1];// [k-1];
				}
				// Interpolation.polint(hh,sss,j-K,0.0);//,&ss,&dss);
				Interpolator.polint(hh, sss, K, 0.0);// ,&ss,&dss);
				ss = Interpolator.ypoli;
				dss = Interpolator.dypoli;
				if (Math.abs(dss) <= EPS * Math.abs(ss))
					return ss;
			}
			h[j] = 0.25 * h[j - 1];// h[j+1]=0.25*h[j];
			// This is a key step: The factor is 0.25 even though the stepsize
			// is decreased by only
			// 0.5. This makes the extrapolation a polynomial in h2 as allowed
			// by equation (4.2.1),
			// not just a polynomial in h.
		}

		failS = "Too many steps in routine qromb";
		failB = true;

		// nrerror("Too many steps in routine qromb");
		return 0.0;// Never get here.
	}

	// @@Improper Integrals
	// Second Euler-Maclaurin summation formula,
	// float midpnt(float (*func)(float), float a, float b, int n)
	/**
	 * This routine computes the nth stage of refinement of an extended midpoint rule. func is 
	 * integrated between limits a and b. When called with n=1, the routine returns the crudest estimate of integral from a to b 
	 * f(x)dx. Subsequent calls with n=2,3,...(in that sequential order) will improve the accuracy of s by adding (2/3) 
	 * × 3n-1 additional interior points. s should not be modified between sequential calls.
	 * @param a a
	 * @param b b
	 * @param n n
	 * @return the result
	 */
	public double midpnt(double a, double b, int n)
	// This routine computes the nth stage of refinement of an extended midpoint
	// rule. func is input
	// as a pointer to the function to be integrated between limits a and b,
	// also input. When called with
	// n=1, the routine returns the crudest estimate of integral from a to b
	// f(x)dx.
	// Subsequent calls with n=2,3,...
	// (in that sequential order) will improve the accuracy of s by adding (2/3)
	// × 3n-1 additional
	// interior points. s should not be modified between sequential calls.
	{
		double x = 0.0;
		double tnm = 0.0;
		double sum = 0.0;
		double del = 0.0;
		double ddel = 0.0;
		// double s=0.0;
		int it = 0;
		int j = 0;

		if (n == 1) {
			// return (s=(b-a)*func.F(0.5*(a+b)));
			return (smidpnt = (b - a) * func.F(0.5 * (a + b)));
		} else {
			for (it = 1, j = 1; j < n - 1; j++)
				it *= 3;
			tnm = it;
			del = (b - a) / (3.0 * tnm);
			ddel = del + del; // The added points alternate in spacing between
								// del and ddel.
			x = a + 0.5 * del;
			sum = 0.0;
			for (j = 1; j <= it; j++) {
				sum += func.F(x);
				x += ddel;
				sum += func.F(x);
				x += del;
			}
			// s=(s+(b-a)*sum/tnm)/2.0;///3.0; //The new sum is combined with
			// the old integral
			// to give a refined integral.
			smidpnt = (smidpnt + (b - a) * sum / tnm) / 3.0;
			return smidpnt;
		}
	}

	/*
	 * The routine midpnt can exactly replace trapzd in a driver routine like
	 * qtrap (§4.2); one simply changes trapzd(func,a,b,j) to midpnt(func,a,b,
	 * j), and perhaps also decreases the parameter JMAX since 3JMAX-1 (from
	 * step tripling) is a much larger number than 2JMAX-1 (step doubling). The
	 * open formula implementation analogous to Simpson’s rule (qsimp in §4.2)
	 * substitutes midpnt for trapzd and decreases JMAX as above, but now also
	 * changes the extrapolation step to be s=(9.0*st-ost)/8.0; since, when the
	 * number of steps is tripled, the error decreases to 1/9th its size, not
	 * 1/4th as with step doubling. Either the modified qtrap or the modified
	 * qsimp will fix the first problem on the list at the beginning of this
	 * section. Yet more sophisticated is to generalize Romberg integration in
	 * like manner
	 * 
	 * A typical invocation (integrating the Bessel function Y0(x) from 0 to 2)
	 * is simply #include "nr.h" float answer; ...
	 * answer=qromo(bessy0,0.0,2.0,midpnt);
	 */
	// float qromo(float (*func)(float), float a, float b,
	// float (*choose)(float(*)(float), float, float, int))
	/**
	 * Romberg integration on an open interval. Returns the integral of the function func from a to b, 
	 * using any specified integrating function choose and Romberg’s method. Normally choose will 
	 * be an open formula, not evaluating the function at the endpoints. It is assumed that choose 
	 * triples the number of steps on each call, and that its error series contains only even powers of 
	 * the number of steps. The routines midpnt, midinf, midsql, midsqu, midexp, are possible 
	 * choices for choose. The parameters have the same meaning as in qromb.
	 * @param a a
	 * @param b b
	 * @param funcs funcs of choice 
	 * @return the result
	 */
	public double qromo(double a, double b, String funcs)
	// Romberg integration on an open interval. Returns the integral of the
	// function func from a to b,
	// using any specified integrating function choose and Romberg’s method.
	// Normally choose will
	// be an open formula, not evaluating the function at the endpoints. It is
	// assumed that choose
	// triples the number of steps on each call, and that its error series
	// contains only even powers of
	// the number of steps. The routines midpnt, midinf, midsql, midsqu, midexp,
	// are possible
	// choices for choose. The parameters have the same meaning as in qromb.
	{
		// void polint(float xa[], float ya[], int n, float x, float *y, float
		// *dy);
		// void nrerror(char error_text[]);
		int j = 0;
		double ss = 0.0;
		double dss = 0.0;
		double[] h = new double[JMAXP1 + 1];
		double[] s = new double[JMAXP];

		failB = false;
		h[0] = 1.0;// h[1]=1.0;
		for (j = 1; j <= JMAX1; j++) {
			s[j - 1] = choose(funcs, a, b, j);// s[j]=choose(funcs,a,b,j);//(*choose)(func,a,b,j);
			if (j >= K) {
				// polint(&h[j-K],&s[j-K],K,0.0,&ss,&dss);
				// &h[j-K]= pointer to the j-K element of vector=> it was taken
				// the elements from j-K to the end!!!
				double[] hh = new double[K];
				double[] sss = new double[K];
				for (int k = 1; k <= K; k++) {
					hh[k - 1] = h[j - K + k - 1];
					sss[k - 1] = s[j - K + k - 1];
				}
				Interpolator.polint(hh, sss, K, 0.0);
				ss = Interpolator.ypoli;
				dss = Interpolator.dypoli;

				if (Math.abs(dss) <= EPS1 * Math.abs(ss))
					return ss;
			}
			// h[j+1]=h[j]/9.0; //This is where the assumption of step tripling
			// and an even
			// error series is used.
			h[j] = h[j - 1] / 9.0;
		}

		failS = "Too many steps in routing qromo";
		failB = true;

		// nrerror("Too many steps in routing qromo");
		return 0.0; // Never get here.
	}

	/*
	 * If you need to integrate from a negative lower limit to positive
	 * infinity, you do this by breaking the integral into two pieces at some
	 * positive value, for example,
	 * answer=qromo(funk,-5.0,2.0,midpnt)+qromo(funk,2.0,1.0e30,midinf);
	 * 
	 * #define FUNC(x) ((*funk)(1.0/(x))/((x)*(x)))
	 */
	// float midinf(float (*funk)(float), float aa, float bb, int n)
	/**
	 * Used by qromo.
	 * This routine is an exact replacement for midpnt, i.e., returns the nth stage of refinement of 
	 * the integral of funk from aa to bb, except that the function is evaluated at evenly spaced 
	 * points in 1/x rather than in x. This allows the upper limit bb to be as large and positive as 
	 * the computer allows, or the lower limit aa to be as large and negative, but not both. aa and 
	 * bb must have the same sign.
	 * @param aa aa
	 * @param bb bb
	 * @param n n
	 * @return the result
	 */
	double midinf(double aa, double bb, int n)
	// This routine is an exact replacement for midpnt, i.e., returns the nth
	// stage of refinement of
	// the integral of funk from aa to bb, except that the function is evaluated
	// at evenly spaced
	// points in 1/x rather than in x. This allows the upper limit bb to be as
	// large and positive as
	// the computer allows, or the lower limit aa to be as large and negative,
	// but not both. aa and
	// bb must have the same sign.
	{
		double x = 0.0;
		double tnm = 0.0;
		double sum = 0.0;
		double del = 0.0;
		double ddel = 0.0;
		double b = 0.0;
		double a = 0.0;
		// double s=0.0;
		int it = 0;
		int j = 0;
		b = 1.0 / aa; // These two statements change the limits of integration.
		a = 1.0 / bb;
		if (n == 1) {// From this point on, the routine is identical to midpnt.
						// return (s=(b-a)*func.F(0.5*(a+b)));
						// return
						// (s=(b-a)*func.F(1.0/(0.5*(a+b)))/((0.5*(a+b))*(0.5*(a+b))));
			return (smidinf = (b - a) * func.F(1.0 / (0.5 * (a + b)))
					/ ((0.5 * (a + b)) * (0.5 * (a + b))));
		} else {
			for (it = 1, j = 1; j < n - 1; j++)
				it *= 3;
			tnm = it;
			del = (b - a) / (3.0 * tnm);
			ddel = del + del;
			x = a + 0.5 * del;
			sum = 0.0;
			for (j = 1; j <= it; j++) {
				// sum += func.F(x);
				sum += func.F(1.0 / x) / (x * x);
				x += ddel;
				// sum += func.F(x);
				sum += func.F(1.0 / x) / (x * x);
				x += del;
			}
			// return (s=(s+(b-a)*sum/tnm)/2.0);//3.0);
			return (smidinf = (smidinf + (b - a) * sum / tnm) / 3.0);
		}
	}

	// #define FUNC(x) (2.0*(x)*(*funk)(aa+(x)*(x)))
	/**
	 * Used by qromo
	 * @param aa aa
	 * @param bb bb
	 * @param n n
	 * @return the result
	 */
	double midsql(double aa, double bb, int n) {
		double x = 0.0;
		double tnm = 0.0;
		double sum = 0.0;
		double del = 0.0;
		double ddel = 0.0;
		double b = 0.0;
		double a = 0.0;
		// double s=0.0;
		int it = 0;
		int j = 0;
		b = Math.sqrt(bb - aa);
		a = 0.0;
		if (n == 1) {// From this point on, the routine is identical to midpnt.
						// return (s=(b-a)*func.F(0.5*(a+b)));
						// return
						// (s=(b-a)*2.0*(0.5*(a+b))*func.F(aa+(0.5*(a+b))*(0.5*(a+b))));
			return (smidsql = (b - a) * 2.0 * (0.5 * (a + b))
					* func.F(aa + (0.5 * (a + b)) * (0.5 * (a + b))));
		} else {
			for (it = 1, j = 1; j < n - 1; j++)
				it *= 3;
			tnm = it;
			del = (b - a) / (3.0 * tnm);
			ddel = del + del;
			x = a + 0.5 * del;
			sum = 0.0;
			for (j = 1; j <= it; j++) {
				// sum += func.F(x);
				sum += 2.0 * x * func.F(aa + (x * x));
				x += ddel;
				// sum += func.F(x);
				sum += 2.0 * x * func.F(aa + (x * x));
				x += del;
			}
			// return (s=(s+(b-a)*sum/tnm)/2.0);//3.0);
			return (smidsql = (smidsql + (b - a) * sum / tnm) / 3.0);
		}
	}

	// /#define FUNC(x) (2.0*(x)*(*funk)(bb-(x)*(x)))
	/**
	 * Used by qromo
	 * @param aa aa
	 * @param bb bb
	 * @param n n
	 * @return the result
	 */
	double midsqu(double aa, double bb, int n) {
		double x = 0.0;
		double tnm = 0.0;
		double sum = 0.0;
		double del = 0.0;
		double ddel = 0.0;
		double b = 0.0;
		double a = 0.0;
		// double s=0.0;
		int it = 0;
		int j = 0;
		b = Math.sqrt(bb - aa);
		a = 0.0;
		if (n == 1) {// From this point on, the routine is identical to midpnt.
						// return (s=(b-a)*func.F(0.5*(a+b)));
						// return
						// (s=(b-a)*2.0*(0.5*(a+b))*func.F(bb-(0.5*(a+b))*(0.5*(a+b))));
			return (smidsqu = (b - a) * 2.0 * (0.5 * (a + b))
					* func.F(bb - (0.5 * (a + b)) * (0.5 * (a + b))));
		} else {
			for (it = 1, j = 1; j < n - 1; j++)
				it *= 3;
			tnm = it;
			del = (b - a) / (3.0 * tnm);
			ddel = del + del;
			x = a + 0.5 * del;
			sum = 0.0;
			for (j = 1; j <= it; j++) {
				// sum += func.F(x);
				sum += 2.0 * x * func.F(bb - (x * x));
				x += ddel;
				// sum += func.F(x);
				sum += 2.0 * x * func.F(bb - (x * x));
				x += del;
			}
			// return (s=(s+(b-a)*sum/tnm)/2.0);//3.0);
			return (smidsqu = (smidsqu + (b - a) * sum / tnm) / 3.0);
		}
	}

	// /#define FUNC(x) ((*funk)(-log(x))/(x))
	/**
	 * Used by qromo
	 * @param aa aa
	 * @param bb bb
	 * @param n n
	 * @return the result
	 */
	double midexp(double aa, double bb, int n) {
		double x = 0.0;
		double tnm = 0.0;
		double sum = 0.0;
		double del = 0.0;
		double ddel = 0.0;
		double b = 0.0;
		double a = 0.0;
		// double s=0.0;
		int it = 0;
		int j = 0;
		b = Math.exp(-aa);
		a = 0.0;
		if (n == 1) {// From this point on, the routine is identical to midpnt.
						// return (s=(b-a)*func.F(0.5*(a+b)));
						// return
						// (s=(b-a)*func.F(-Math.log(0.5*(a+b)))/(0.5*(a+b)));
			return (smidexp = (b - a) * func.F(-Math.log(0.5 * (a + b)))
					/ (0.5 * (a + b)));
		} else {
			for (it = 1, j = 1; j < n - 1; j++)
				it *= 3;
			tnm = it;
			del = (b - a) / (3.0 * tnm);
			ddel = del + del;
			x = a + 0.5 * del;
			sum = 0.0;
			for (j = 1; j <= it; j++) {
				// sum += func.F(x);
				sum += func.F(-Math.log(x)) / x;
				x += ddel;
				// sum += func.F(x);
				sum += func.F(-Math.log(x)) / x;
				x += del;
			}
			// return (s=(s+(b-a)*sum/tnm)/2.0);//3.0);
			return (smidexp = (smidexp + (b - a) * sum / tnm) / 3.0);
		}
	}

	/**
	 * Used by qromo
	 * @param s s
	 * @param a a
	 * @param b b
	 * @param j j
	 * @return the result
	 */
	 
	private double choose(String s, double a, double b, int j) {
		double res = 0.0;
		if (s.compareTo("midpnt") == 0)
			res = midpnt(a, b, j);
		else if (s.compareTo("trapzd") == 0)
			res = trapzd(a, b, j);
		else if (s.compareTo("midinf") == 0)
			res = midinf(a, b, j);
		else if (s.compareTo("midsql") == 0)
			res = midsql(a, b, j);
		else if (s.compareTo("midsqu") == 0)
			res = midsqu(a, b, j);
		else if (s.compareTo("midexp") == 0)
			res = midexp(a, b, j);

		return res;
	}

	// @@@@@@Gaussian Quadratures and Orthogonal Polynomials
	// float qgaus(float (*func)(float), float a, float b)
	/**
	 * Returns the integral of the function func between a and b, by ten-point Gauss-Legendre integration: 
	 * the function is evaluated exactly ten times at interior points in the range of integration.
	 * @param a a
	 * @param b b
	 * @return the result
	 */
	public double qgaus(double a, double b)
	// Returns the integral of the function func between a and b, by ten-point
	// Gauss-Legendre integration:
	// the function is evaluated exactly ten times at interior points in the
	// range of integration.
	{
		int j = 0;
		double xr = 0.0;
		double xm = 0.0;
		double dx = 0.0;
		double s = 0.0;
		// The abscissas and weights. First value of each array not used.
		double[] x = { 0.0, 0.1488743389, 0.4333953941, 0.6794095682,
				0.8650633666, 0.9739065285 };
		double[] w = { 0.0, 0.2955242247, 0.2692667193, 0.2190863625,
				0.1494513491, 0.0666713443 };
		xm = 0.5 * (b + a);
		xr = 0.5 * (b - a);
		s = 0;
		// Will be twice the average value of the function, since the
		// ten weights (five numbers above each used twice)sum to 2.
		for (j = 1; j <= 5; j++) {
			dx = xr * x[j];// dx=xr*x[j];
			// s += w[j]*((*func)(xm+dx)+(*func)(xm-dx));
			s += w[j] * (func.F(xm + dx) + func.F(xm - dx));
		}
		return s *= xr; // Scale the answer to the range of integration.

		/*
		 * int n=5; double[]
		 * w={0.23692688,0.47862868,0.56888889,0.47862868,0.23692688}; double[]
		 * x={-0.90617985,-0.53846931,0.0,0.53846931,0.90617985}; b=0.5*(b+a);
		 * a=b-a; double s=0.0; for (int i=1; i<=n; i++)
		 * s=s+w[i-1]*func.F(a*x[i-1]+b); s=a*s; return s;
		 */
	}

	/**
	 * Given the lower and upper limits of integration x1 and x2, and given n, this routine returns 
	 * arrays x[1..n] and w[1..n] of length n, containing the abscissas and weights of the Gauss-
	 * Legendre n-point quadrature formula.
	 * @param x1 x1
	 * @param x2 x2
	 * @param x x
	 * @param w w
	 * @param n n
	 */
	public static void gauleg(double x1, double x2, double[] x, double[] w,
			int n)
	// Given the lower and upper limits of integration x1 and x2, and given n,
	// this routine returns
	// arrays x[1..n] and w[1..n] of length n, containing the abscissas and
	// weights of the Gauss-
	// Legendre n-point quadrature formula.
	{
		int m = 0;
		int j = 0;
		int i = 0;
		double z1 = 0.0;
		double z = 0.0;
		double xm = 0.0;
		double xl = 0.0;
		double pp = 0.0;
		double p3 = 0.0;
		double p2 = 0.0;
		double p1 = 0.0;
		// High precision is a good idea for this routine.
		m = (n + 1) / 2;
		// The roots are symmetric in the interval, so we only have to find half
		// of them.
		xm = 0.5 * (x2 + x1);
		xl = 0.5 * (x2 - x1);
		for (i = 1; i <= m; i++) {// Loop over the desired roots.
			z = Math.cos(3.141592654 * (i - 0.25) / (n + 0.5));
			// Starting with the above approximation to the ith root, we enter
			// the main loop of
			// refinement by Newton’s method.
			do {
				p1 = 1.0;
				p2 = 0.0;
				for (j = 1; j <= n; j++) {
					// Loop up the recurrence relation to get the
					// Legendre polynomial evaluated at z.
					p3 = p2;
					p2 = p1;
					p1 = ((2.0 * j - 1.0) * z * p2 - (j - 1.0) * p3) / j;
				}
				// p1 is now the desired Legendre polynomial. We next compute
				// pp, its derivative,
				// by a standard relation involving also p2, the polynomial of
				// one lower order.
				pp = n * (z * p1 - p2) / (z * z - 1.0);
				z1 = z;
				z = z1 - p1 / pp;// Newton’s method.
			} while (Math.abs(z - z1) > EPS2);
			x[i - 1] = xm - xl * z;// x[i]=xm-xl*z;// Scale the root to the
									// desired interval,
			x[n - i] = xm + xl * z;// x[n+1-i]=xm+xl*z; //and put in its
									// symmetric counterpart.
			w[i - 1] = 2.0 * xl / ((1.0 - z * z) * pp * pp);// w[i]=2.0*xl/((1.0-z*z)*pp*pp);//
															// Compute the
															// weight
			w[n - i] = w[i - 1];// w[n+1-i]=w[i]; and its symmetric counterpart.
		}
	}

	/*
	 * Next we give three routines that use initial approximations for the roots
	 * given by Stroud and Secrest [2]. The first is for Gauss-Laguerre
	 * abscissas and weights, to be used with the integration formula
	 */
	/**
	 * Returns the value ln[GAMMA(xx)] for xx greater than 0.
	 * @param xx xx
	 * @return the result
	 */
	public static double gammln(double xx)
	// Returns the value ln[?(xx)] for xx > 0.
	{
		// Internal arithmetic will be done in double precision, a nicety that
		// you can omit if five-figure
		// accuracy is good enough.
		double x = 0.0;
		double y = 0.0;
		double tmp = 0.0;
		double ser = 0.0;
		// cof->[6]
		double[] cof = { 76.18009172947146, -86.50532032941677,
				24.01409824083091, -1.231739572450155, 0.1208650973866179e-2,
				-0.5395239384953e-5 };
		int j = 0;
		y = x = xx;
		tmp = x + 5.5;
		tmp -= (x + 0.5) * Math.log(tmp);
		ser = 1.000000000190015;
		for (j = 0; j <= 5; j++)
			ser += cof[j] / ++y;
		return -tmp + Math.log(2.5066282746310005 * ser / x);// =>tested OK!!

		/*
		 * //mine old int j; double stp = 2.506628274650; double cof[] = new
		 * double[6]; cof[0]=76.18009173; cof[1]=-86.50532033;
		 * cof[2]=24.01409822; cof[3]=-1.231739516; cof[4]=0.120858003E-02;
		 * cof[5]=-0.536382E-05;
		 * 
		 * double x = xx-1; double tmp = x + 5.5; tmp = (x + 0.5)*Math.log(tmp)
		 * - tmp; double ser = 1.0; for(j=0;j<6;j++) { x++; ser = ser +
		 * cof[j]/x; } double retVal = tmp + Math.log(stp*ser); return retVal;
		 */
	}

	/**
	 * Given alf, the parameter of the Laguerre polynomials, this routine returns arrays x[1..n] 
	 * and w[1..n] containing the abscissas and weights of the n-point Gauss-Laguerre quadrature 
	 * formula. The smallest abscissa is returned in x[1], the largest in x[n].
	 * @param x x
	 * @param w w
	 * @param n n
	 * @param alf alf
	 */
	public static void gaulag(double[] x, double[] w, int n, double alf)
	// Given alf, the parameter of the Laguerre polynomials, this routine
	// returns arrays x[1..n]
	// and w[1..n] containing the abscissas and weights of the n-point
	// Gauss-Laguerre quadrature
	// formula. The smallest abscissa is returned in x[1], the largest in x[n].
	{
		// float gammln(float xx);=>@@@@@@@@@@@@@@NEXT CHAPTERS!!!
		// void nrerror(char error_text[]);
		int i = 0;
		int its = 0;
		int j = 0;
		double ai = 0.0;
		double p1 = 0.0;
		double p2 = 0.0;
		double p3 = 0.0;
		double pp = 0.0;
		double z = 0.0;
		double z1 = 0.0;// High precision is a good idea for this routine.

		failB = false;

		for (i = 1; i <= n; i++) {
			// Loop over the desired roots.
			if (i == 1) {// Initial guess for the smallest root.
				z = (1.0 + alf) * (3.0 + 0.92 * alf)
						/ (1.0 + 2.4 * n + 1.8 * alf);
			} else if (i == 2) {// Initial guess for the second root.
				z += (15.0 + 6.25 * alf) / (1.0 + 0.9 * alf + 2.5 * n);
			} else {// Initial guess for the other roots.
				ai = i - 2.0;
				// z +=
				// ((1.0+2.55*ai)/(1.9*ai)+1.26*ai*alf/(1.0+3.5*ai))*(z-x[i-2])/(1.0+0.3*alf);
				z += ((1.0 + 2.55 * ai) / (1.9 * ai) + 1.26 * ai * alf
						/ (1.0 + 3.5 * ai))
						* (z - x[i - 3]) / (1.0 + 0.3 * alf);
			}

			for (its = 1; its <= MAXIT; its++) {// Refinement by Newton’s
												// method.
				p1 = 1.0;
				p2 = 0.0;
				for (j = 1; j <= n; j++) { // Loop up the recurrence relation to
											// get the Laguerre polynomial
											// evaluated at z.
					p3 = p2;
					p2 = p1;
					p1 = ((2 * j - 1 + alf - z) * p2 - (j - 1 + alf) * p3) / j;
				}
				// p1 is now the desired Laguerre polynomial. We next compute
				// pp, its derivative,
				// by a standard relation involving also p2, the polynomial of
				// one lower order.
				pp = (n * p1 - (n + alf) * p2) / z;
				z1 = z;
				z = z1 - p1 / pp;// Newton’s formula.
				if (Math.abs(z - z1) <= EPS3)
					break;
			}
			if (its > MAXIT) {
				// nrerror("too many iterations in gaulag");
				failS = "too many iterations in gaulag";
				failB = true;
				return;
			}
			x[i - 1] = z;// x[i]=z; Store the root and the weight.
			// w[i] = -exp(gammln(alf+n)-gammln((float)n))/(pp*n*p2);
			w[i - 1] = -Math.exp(gammln(alf + n) - gammln((float) n))
					/ (pp * n * p2);
		}
	}

	/**
	 * Given n, this routine returns arrays x[1..n] and w[1..n] containing the abscissas and weights 
	 * of the n-point Gauss-Hermite quadrature formula. The largest abscissa is returned in x[1], the most negative in x[n].
	 * @param x x
	 * @param w w
	 * @param n n
	 */
	public static void gauher(double[] x, double[] w, int n)
	// Given n, this routine returns arrays x[1..n] and w[1..n] containing the
	// abscissas and weights
	// of the n-point Gauss-Hermite quadrature formula. The largest abscissa is
	// returned in x[1], the
	// most negative in x[n].
	{
		// void nrerror(char error_text[]);
		int i = 0;
		int its = 0;
		int j = 0;
		int m = 0;
		double p1 = 0.0;
		double p2 = 0.0;
		double p3 = 0.0;
		double pp = 0.0;
		double z = 0.0;
		double z1 = 0.0;
		// High precision is a good idea for this routine.

		failB = false;

		m = (n + 1) / 2;
		double nd = n;
		// The roots are symmetric about the origin, so we have to find only
		// half of them.
		for (i = 1; i <= m; i++) {// Loop over the desired roots.
			if (i == 1) {// Initial guess for the largest root.
				z = Math.sqrt((2.0 * nd + 1.0)) - 1.85575
						* Math.pow((2.0 * nd + 1.0), -0.16667);
			} else if (i == 2) {// Initial guess for the second largest root.
				z -= 1.14 * Math.pow(nd, 0.426) / z;
			} else if (i == 3) {// Initial guess for the third largest root.
				z = 1.86 * z - 0.86 * x[0];// z=1.86*z-0.86*x[1];
			} else if (i == 4) {// Initial guess for the fourth largest root.
				z = 1.91 * z - 0.91 * x[1];// z=1.91*z-0.91*x[2];
			} else {// Initial guess for the other roots.
				z = 2.0 * z - x[i - 3];// z=2.0*z-x[i-2];
			}
			for (its = 1; its <= MAXIT; its++) {// Refinement by Newton’s
												// method.
				p1 = PIM4;
				p2 = 0.0;
				for (j = 1; j <= n; j++) {
					// Loop up the recurrence relation to get
					// the Hermite polynomial evaluated at z.
					p3 = p2;
					p2 = p1;
					p1 = z * Math.sqrt(2.0 / j) * p2
							- Math.sqrt(((j - 1.0)) / j) * p3;
				}
				// p1 is now the desired Hermite polynomial. We next compute pp,
				// its derivative, by
				// the relation (4.5.21) using p2, the polynomial of one lower
				// order.
				pp = Math.sqrt(2.0 * nd) * p2;
				z1 = z;
				z = z1 - p1 / pp;// Newton’s formula.
				if (Math.abs(z - z1) <= EPS3)
					break;
			}
			if (its > MAXIT) {// nrerror("too many iterations in gauher");
				failS = "too many iterations in gauher";
				failB = true;
				return;
			}
			x[i - 1] = z;// x[i]=z; //Store the root
			x[n - i] = -z;// x[n+1-i] = -z; and its symmetric counterpart.
			w[i - 1] = 2.0 / (pp * pp);// w[i]=2.0/(pp*pp); Compute the weight
			w[n - i] = w[i - 1];// w[n+1-i]=w[i]; and its symmetric counterpart.
		}
	}

	/**
	 * Given alf and bet, the parameters ALPHA and BETA of the Jacobi polynomials, this routine returns 
	 * arrays x[1..n] and w[1..n] containing the abscissas and weights of the n-point Gauss-Jacobi 
	 * quadrature formula. The largest abscissa is returned in x[1], the smallest in x[n].
	 * @param x x
	 * @param w w
	 * @param n n
	 * @param alf alf
	 * @param bet bet
	 */
	public static void gaujac(double[] x, double[] w, int n, double alf,
			double bet)
	// Given alf and bet, the parameters ALPHA and BETA of the Jacobi polynomials, this
	// routine returns
	// arrays x[1..n] and w[1..n] containing the abscissas and weights of the
	// n-point Gauss-Jacobi
	// quadrature formula. The largest abscissa is returned in x[1], the
	// smallest in x[n].
	{
		// float gammln(float xx);
		// void nrerror(char error_text[]);
		int i = 0;
		int its = 0;
		int j = 0;
		double alfbet = 0.0;
		double an = 0.0;
		double bn = 0.0;
		double r1 = 0.0;
		double r2 = 0.0;
		double r3 = 0.0;
		double a = 0.0;
		double b = 0.0;
		double c = 0.0;
		double p1 = 0.0;
		double p2 = 0.0;
		double p3 = 0.0;
		double pp = 0.0;
		double temp = 0.0;
		double z = 0.0;
		double z1 = 0.0;
		// High precision is a good idea for this routine.

		failB = false;

		for (i = 1; i <= n; i++) {// Loop over the desired roots.
			if (i == 1) {// Initial guess for the largest root.
				an = alf / n;
				bn = bet / n;
				r1 = (1.0 + alf) * (2.78 / (4.0 + n * n) + 0.768 * an / n);
				r2 = 1.0 + 1.48 * an + 0.96 * bn + 0.452 * an * an + 0.83 * an
						* bn;
				z = 1.0 - r1 / r2;
			} else if (i == 2) {// Initial guess for the second largest root.
				r1 = (4.1 + alf) / ((1.0 + alf) * (1.0 + 0.156 * alf));
				r2 = 1.0 + 0.06 * (n - 8.0) * (1.0 + 0.12 * alf) / n;
				r3 = 1.0 + 0.012 * bet * (1.0 + 0.25 * Math.abs(alf)) / n;
				z -= (1.0 - z) * r1 * r2 * r3;
			} else if (i == 3) {// Initial guess for the third largest root.
				r1 = (1.67 + 0.28 * alf) / (1.0 + 0.37 * alf);
				r2 = 1.0 + 0.22 * (n - 8.0) / n;
				r3 = 1.0 + 8.0 * bet / ((6.28 + bet) * n * n);
				z -= (x[0] - z) * r1 * r2 * r3;// z -= (x[1]-z)*r1*r2*r3;
			} else if (i == n - 1) {// Initial guess for the second smallest
									// root.
				r1 = (1.0 + 0.235 * bet) / (0.766 + 0.119 * bet);
				r2 = 1.0 / (1.0 + 0.639 * (n - 4.0) / (1.0 + 0.71 * (n - 4.0)));
				r3 = 1.0 / (1.0 + 20.0 * alf / ((7.5 + alf) * n * n));
				// z += (z-x[n-3])*r1*r2*r3;
				z += (z - x[n - 4]) * r1 * r2 * r3;
			} else if (i == n) {// Initial guess for the smallest root.
				r1 = (1.0 + 0.37 * bet) / (1.67 + 0.28 * bet);
				r2 = 1.0 / (1.0 + 0.22 * (n - 8.0) / n);
				r3 = 1.0 / (1.0 + 8.0 * alf / ((6.28 + alf) * n * n));
				// z += (z-x[n-2])*r1*r2*r3;
				z += (z - x[n - 3]) * r1 * r2 * r3;
			} else {// Initial guess for the other roots.
				z = 3.0 * x[i - 2] - 3.0 * x[i - 3] + x[i - 4];// z=3.0*x[i-1]-3.0*x[i-2]+x[i-3];
			}
			alfbet = alf + bet;
			for (its = 1; its <= MAXIT; its++) {// Refinement by Newton’s
												// method.
				temp = 2.0 + alfbet;
				// Start the recurrence with P0 and P1 to avoid
				// a division by zero when ? + ß = 0 or-1.
				p1 = (alf - bet + temp * z) / 2.0;
				p2 = 1.0;
				for (j = 2; j <= n; j++) {// Loop up the recurrence relation to
											// get the
											// Jacobi polynomial evaluated at z.
					p3 = p2;
					p2 = p1;
					temp = 2.0 * j + alfbet;
					a = 2.0 * j * (j + alfbet) * (temp - 2.0);
					b = (temp - 1.0)
							* (alf * alf - bet * bet + temp * (temp - 2.0) * z);
					c = 2.0 * (j - 1.0 + alf) * (j - 1.0 + bet) * temp;
					p1 = (b * p2 - c * p3) / a;
				}
				pp = (n * (alf - bet - temp * z) * p1 + 2.0 * (n + alf)
						* (n + bet) * p2)
						/ (temp * (1.0 - z * z));
				// p1 is now the desired Jacobi polynomial. We next compute pp,
				// its derivative, by
				// a standard relation involving also p2, the polynomial of one
				// lower order.
				z1 = z;
				z = z1 - p1 / pp;// Newton’s formula.
				if (Math.abs(z - z1) <= EPS3)
					break;
			}
			if (its > MAXIT) {
				// nrerror("too many iterations in gaujac");
				failS = "too many iterations in gaujac";
				failB = true;
				return;
			}
			x[i - 1] = z;// x[i]=z; Store the root and the weight.
			// w[i]=exp(gammln(alf+n)+gammln(bet+n)-gammln(n+1.0)-
			// gammln(n+alfbet+1.0))*temp*pow(2.0,alfbet)/(pp*p2);
			w[i - 1] = Math.exp(gammln(alf + n) + gammln(bet + n)
					- gammln(n + 1.0) - gammln(n + alfbet + 1.0))
					* temp * Math.pow(2.0, alfbet) / (pp * p2);
		}
	}

	/*
	 * ] Case of Known Recurrences Turn now to the case where you do not know
	 * good initial guesses for the zeros of your orthogonal polynomials, but
	 * you do have available the coefficients aj and bj that generate them. As
	 * we have seen, the zeros of pN(x) are the abscissas for the N-point
	 * Gaussian quadrature formula. The most useful computational formula for
	 * the weights is equation (4.5.9) above, since the derivative pN can be
	 * efficiently computed by the derivative of (4.5.6) in the general case, or
	 * by special relations for the classical polynomials. Note that (4.5.9) is
	 * valid as written only for monic polynomials; for other normalizations,
	 * there is an extra factor of ?N/?N-1, where ?N is the coefficient of xN in
	 * pN.
	 */

	/**
	 * Given the eigenvalues d[1..n] and eigenvectors v[1..n][1..n] as output from jacobi 
	 * this routine sorts the eigenvalues into escending order, and rearranges 
	 * the columns of v correspondingly. The method is straight insertion.
	 * @param d d
	 * @param v v
	 * @param n n
	 */
	public static void eigsrt(double[] d, double[][] v, int n)
	// Given the eigenvalues d[1..n] and eigenvectors v[1..n][1..n] as output
	// from jacobi
	// (§11.1) or tqli (§11.3), this routine sorts the eigenvalues into
	// descending order, and rearranges
	// the columns of v correspondingly. The method is straight insertion.
	{
		int k = 0;
		int j = 0;
		int i = 0;
		double p = 0.0;
		for (i = 1; i < n; i++) {
			k = i;
			p = d[i - 1];
			// p=d[k=i];
			for (j = i + 1; j <= n; j++)
				// if (d[j] >= p) p=d[k=j];
				if (d[j - 1] >= p) {
					k = j;
					p = d[j - 1];
				}
			if (k != i) {
				d[k - 1] = d[i - 1];// d[k]=d[i];
				d[i - 1] = p;// d[i]=p;
				for (j = 1; j <= n; j++) {
					p = v[j - 1][i - 1];// p=v[j][i];
					v[j - 1][i - 1] = v[j - 1][k - 1];// v[j][i]=v[j][k];
					v[j - 1][k - 1] = p;// v[j][k]=p;
				}
			}
		}
	}

	/**
	 * QL algorithm with implicit shifts, to determine the eigenvalues and eigenvectors of a real, symmetric, 
	 * tridiagonal matrix, or of a real, symmetric matrix previously reduced by tred2. On 
	 * input, d[1..n] contains the diagonal elements of the tridiagonal matrix. On output, it returns 
	 * the eigenvalues. The vector e[1..n] inputs the subdiagonal elements of the tridiagonal matrix, 
	 * with e[1] arbitrary. On output e is destroyed. When finding only the eigenvalues, several lines 
	 * may be omitted, as noted in the comments. If the eigenvectors of a tridiagonal matrix are desired, 
	 * the matrix z[1..n][1..n] is input as the identity matrix. If the eigenvectors of a matrix 
	 * that has been reduced by tred2 are required, then z is input as the matrix output by tred2. 
	 * In either case, the kth column of z returns the normalized eigenvector corresponding to d[k].
	 * @param d d
	 * @param e e
	 * @param n n
	 * @param z z
	 */
	public static void tqli(double[] d, double[] e, int n, double[][] z)
	// QL algorithm with implicit shifts, to determine the eigenvalues and
	// eigenvectors of a real, symmetric,
	// tridiagonal matrix, or of a real, symmetric matrix previously reduced by
	// tred2 §11.2. On
	// input, d[1..n] contains the diagonal elements of the tridiagonal matrix.
	// On output, it returns
	// the eigenvalues. The vector e[1..n] inputs the subdiagonal elements of
	// the tridiagonal matrix,
	// with e[1] arbitrary. On output e is destroyed. When finding only the
	// eigenvalues, several lines
	// may be omitted, as noted in the comments. If the eigenvectors of a
	// tridiagonal matrix are desired,
	// the matrix z[1..n][1..n] is input as the identity matrix. If the
	// eigenvectors of a matrix
	// that has been reduced by tred2 are required, then z is input as the
	// matrix output by tred2.
	// In either case, the kth column of z returns the normalized eigenvector
	// corresponding to d[k].
	{
		// float pythag(float a, float b);=>LinAEq.pythag(double a, double b)
		int m = 0;
		int l = 0;
		int iter = 0;
		int i = 0;
		int k = 0;
		double s = 0.0;
		double r = 0.0;
		double p = 0.0;
		double g = 0.0;
		double f = 0.0;
		double dd = 0.0;
		double c = 0.0;
		double b = 0.0;

		failB = false;

		for (i = 2; i <= n; i++)
			e[i - 2] = e[i - 1];// e[i-1]=e[i]; //Convenient to renumber the
								// elements of e.
		e[n - 1] = 0.0;// e[n]=0.0;
		for (l = 1; l <= n; l++) {
			iter = 0;
			do {
				for (m = l; m <= n - 1; m++) {// Look for a single small
												// subdiagonal
												// element to split the matrix.
					dd = Math.abs(d[m - 1]) + Math.abs(d[m]);// dd=fabs(d[m])+fabs(d[m+1]);
					if (Math.abs(e[m - 1]) + dd == dd)
						break;// if ((float)(fabs(e[m])+dd) == dd) break;
				}
				if (m != l) {
					if (iter++ == 30) {
						// nrerror("Too many iterations in tqli");
						failS = "too many iterations in gaujac";
						failB = true;
						return;
					}
					g = (d[l] - d[l - 1]) / (2.0 * e[l - 1]);// g=(d[l+1]-d[l])/(2.0*e[l]);
																// Form shift.
					r = LinAEq.pythag(g, 1.0);
					// g=d[m]-d[l]+e[l]/(g+SIGN(r,g)); //This is dm - ks.

					double dbl = (g >= 0.0 ? r : -r);
					g = d[m - 1] - d[l - 1] + e[l - 1] / (g + dbl);
					s = c = 1.0;
					p = 0.0;
					for (i = m - 1; i >= l; i--) {// A plane rotation as in the
													// original
													// QL, followed by Givens
													// rotations to restore
													// tridiagonal form.
						f = s * e[i - 1];// f=s*e[i];
						b = c * e[i - 1];// b=c*e[i];
						e[i] = (r = LinAEq.pythag(f, g));// e[i+1]=(r=pythag(f,g));
						if (r == 0.0) {// Recover from underflow.
							d[i] -= p;// d[i+1] -= p;
							e[m - 1] = 0.0;// e[m]=0.0;
							break;
						}
						s = f / r;
						c = g / r;
						g = d[i] - p;// g=d[i+1]-p;
						r = (d[i - 1] - g) * s + 2.0 * c * b;// r=(d[i]-g)*s+2.0*c*b;
						d[i] = g + (p = s * r);// d[i+1]=g+(p=s*r);
						g = c * r - b;
						// Next loop can be omitted if eigenvectors not wanted
						for (k = 1; k <= n; k++) {// Form eigenvectors.
							f = z[k - 1][i];// f=z[k][i+1];
							z[k - 1][i] = s * z[k - 1][i - 1] + c * f;// z[k][i+1]=s*z[k][i]+c*f;
							z[k - 1][i - 1] = c * z[k - 1][i - 1] - s * f;// z[k][i]=c*z[k][i]-s*f;
						}
					}
					if (r == 0.0 && i >= l)
						continue;
					d[l - 1] -= p;// d[l] -= p;
					e[l - 1] = g;// e[l]=g;
					e[m - 1] = 0.0;// e[m]=0.0;
				}
			} while (m != l);
		}
	}

	/**
	 * Computes the abscissas and weights for a Gaussian quadrature formula from the Jacobi matrix. 
	 * On input, a[1..n] and b[1..n] are the coefficients of the recurrence relation for the set of 
	 * monic orthogonal polynomials. The quantity µ0 =integeral from a to b of W(x) dx is input as amu0. The abscissas 
	 * x[1..n] are returned in descending order, with the corresponding weights in w[1..n]. The 
	 * arrays a and b are modified. Execution can be speeded up by modifying tqli and eigsrt to 
	 * compute only the first component of each eigenvector.
	 * @param n n
	 * @param a a
	 * @param b b
	 * @param amu0 amu0
	 * @param x x
	 * @param w w
	 */
	public static void gaucof(int n, double[] a, double[] b, double amu0,
			double[] x, double[] w)
	// Computes the abscissas and weights for a Gaussian quadrature formula from
	// the Jacobi matrix.
	// On input, a[1..n] and b[1..n] are the coefficients of the recurrence
	// relation for the set of
	// monic orthogonal polynomials. The quantity µ0 =integeral from a to b of
	// W(x) dx is input as amu0. The abscissas
	// x[1..n] are returned in descending order, with the corresponding weights
	// in w[1..n]. The
	// arrays a and b are modified. Execution can be speeded up by modifying
	// tqli and eigsrt to
	// compute only the first component of each eigenvector.
	{
		// void eigsrt(float d[], float **v, int n);
		// void tqli(float d[], float e[], int n, float **z);
		int i = 0;
		int j = 0;
		double[][] z = new double[n][n];
		// z=matrix(1,n,1,n);
		for (i = 1; i <= n; i++) {
			if (i != 1)
				b[i - 1] = Math.sqrt(b[i - 1]);// b[i]=sqrt(b[i]); //Set up
												// superdiagonal of Jacobi
												// matrix.
			for (j = 1; j <= n; j++) // z[i-1][j-1]=(double)(i ==
										// j);//z[i][j]=(float)(i == j);
			{
				if (i == j)
					z[i - 1][j - 1] = 1.0;
				else
					z[i - 1][j - 1] = 0.0;
			}
			// Set up identity matrix for tqli to compute eigenvectors.
		}
		tqli(a, b, n, z);
		eigsrt(a, z, n); // Sort eigenvalues into descending order.
		for (i = 1; i <= n; i++) {
			x[i - 1] = a[i - 1];// x[i]=a[i];
			w[i - 1] = amu0 * z[0][i - 1] * z[0][i - 1];// w[i]=amu0*z[1][i]*z[1][i];
														// //Equation (4.5.27).
		}
		// free_matrix(z,1,n,1,n);
	}

	// @@@@Orthogonal Polynomials with Nonclassical Weights
	/**
	 * Computes the coefficients aj and bj, j = 0, . . .N - 1, of the recurrence 
	 * relation for monic orthogonal polynomials with weight function W(x) by Wheeler’s algorithm. 
	 * On input, the arrays alpha[1..2*n-1] and beta[1..2*n-1] are the coefficients 
	 * ALPHAj and BETAj, j = 0, . . . 2N-2, of the recurrence relation for the chosen basis of orthogonal polynomials. 
	 * The modified moments NIUj are input in anu[1..2*n]. The first n coefficients are returned in 
	 * a[1..n] and b[1..n].
	 * @param n n
	 * @param anu anu
	 * @param alpha alpha
	 * @param beta beta
	 * @param a a
	 * @param b b
	 */
	public static void orthog(int n, double[] anu, double[] alpha,
			double[] beta, double[] a, double[] b)
	// Computes the coefficients aj and bj, j = 0, . . .N - 1, of the recurrence
	// relation for monic
	// orthogonal polynomials with weight function W(x) by Wheeler’s algorithm.
	// On input, the arrays
	// alpha[1..2*n-1] and beta[1..2*n-1] are the coefficients ?j and ßj, j = 0,
	// . . . 2N-2, of
	// the recurrence relation for the chosen basis of orthogonal polynomials.
	// The modified moments
	// ?j are input in anu[1..2*n]. The first n coefficients are returned in
	// a[1..n] and b[1..n].
	{
		int k = 0;
		int l = 0;
		double[][] sig = new double[2 * n + 1][2 * n + 1];
		int looptmp = 0;
		// sig=matrix(1,2*n+1,1,2*n+1);
		looptmp = 2 * n;
		for (l = 3; l <= looptmp; l++)
			sig[0][l - 1] = 0.0;// sig[1][l]=0.0; Initialization, Equation
								// (4.5.33).
		looptmp++;
		for (l = 2; l <= looptmp; l++)
			sig[1][l - 1] = anu[l - 2];// sig[2][l]=anu[l-1];
		a[0] = alpha[0] + anu[1] / anu[0];// a[1]=alpha[1]+anu[2]/anu[1];
		b[0] = 0.0;// b[1]=0.0;
		for (k = 3; k <= n + 1; k++) {// Equation (4.5.34).
			looptmp = 2 * n - k + 3;
			for (l = k; l <= looptmp; l++) {
				// sig[k][l]=sig[k-1][l+1]+(alpha[l-1]-a[k-2])*sig[k-1][l]-
				// b[k-2]*sig[k-2][l]+beta[l-1]*sig[k-1][l-1];
				sig[k - 1][l - 1] = sig[k - 2][l] + (alpha[l - 2] - a[k - 3])
						* sig[k - 2][l - 1] - b[k - 3] * sig[k - 3][l - 1]
						+ beta[l - 2] * sig[k - 2][l - 2];

			}
			// a[k-1]=alpha[k-1]+sig[k][k+1]/sig[k][k]-sig[k-1][k]/sig[k-1][k-1];
			a[k - 2] = alpha[k - 2] + sig[k - 1][k] / sig[k - 1][k - 1]
					- sig[k - 2][k - 1] / sig[k - 2][k - 2];
			// b[k-1]=sig[k][k]/sig[k-1][k-1];
			b[k - 2] = sig[k - 1][k - 1] / sig[k - 2][k - 2];
		}
		// free_matrix(sig,1,2*n+1,1,2*n+1);
	}

	// @@@@@@@Multidimensional Integrals
	/*
	 * I ?   dx dy dzf(x, y, z) =  x2 x1 dx  y2(x) y1(x) dy  z2(x,y)
	 * z1(x,y) dz f(x, y, z)
	 */

	// static float (*nrfunc)(float,float,float);
	// float quad3d(float (*func)(float, float, float), float x1, float x2)
	/**
	 * Returns the integral of a user-supplied function func over a three-dimensional region specified 
	 * by the limits x1, x2, and by the user-supplied functions yy1, yy2, z1, and z2. 
	 * Integration is performed by calling qgaus recursively.
	 * @param x1 x1
	 * @param x2 x2
	 * @return the result
	 */
	public double quad3d(double x1, double x2)
	// Returns the integral of a user-supplied function func over a
	// three-dimensional region specified
	// by the limits x1, x2, and by the user-supplied functions yy1, yy2, z1,
	// and z2, as defined in
	// (4.6.2). (The functions y1 and y2 are here called yy1 and yy2 to avoid
	// conflict with the names
	// of Bessel functions in some C libraries). Integration is performed by
	// calling qgaus recursively.
	{
		// float qgaus(float (*func)(float), float a, float b);
		// float f1(float x);
		// nrfunc=func;
		return qgaus2("f1", x1, x2);
	}

	/**
	 * used in qgaus2
	 * @param x x
	 * @return the result
	 */
	public double f1(double x) // This is H of eq. (4.6.5).
	{
		// float qgaus(float (*func)(float), float a, float b);
		// float f2(float y);
		// float yy1(float),yy2(float);
		xsav = x;
		// return qgaus(f2,yy1(x),yy2(x));
		return qgaus2("f2", func.yy1(x), func.yy2(x));
	}

	/**
	 * used in qgaus2
	 * @param y y
	 * @return the result
	 */
	public double f2(double y) // This is G of eq. (4.6.4).
	{
		// float qgaus(float (*func)(float), float a, float b);
		// float f3(float z);
		// float z1(float,float),z2(float,float);
		ysav = y;
		return qgaus2("f3", func.z1(xsav, y), func.z2(xsav, y));
	}

	/**
	 * used in qgaus2
	 * @param z z
	 * @return the result
	 */
	public double f3(double z) // The integrand f(x, y, z) evaluated at fixed x
								// and y.
	{
		// return (*nrfunc)(xsav,ysav,z);
		return func.F3D(xsav, ysav, z);
	}

	// The necessary user-supplied functions have the following prototypes:
	// float func(float x,float y,float z); The 3-dimensional function to be
	// integrated.
	// float yy1(float x);
	// float yy2(float x);
	// float z1(float x,float y);
	// float z2(float x,float y);
	/**
	 * Returns the integral of the function func between a and b, by ten-point 
	 * Gauss-Legendre integration: the function is evaluated exactly ten times at interior points in the 
	 * range of integration.
	 * @param fs fs
	 * @param a a
	 * @param b b
	 * @return the result
	 */
	public double qgaus2(String fs, double a, double b)
	// Returns the integral of the function func between a and b, by ten-point
	// Gauss-Legendre integration:
	// the function is evaluated exactly ten times at interior points in the
	// range of integration.
	{
		int j = 0;
		double xr = 0.0;
		double xm = 0.0;
		double dx = 0.0;
		double s = 0.0;
		// The abscissas and weights. First value of each array not used.
		double[] x = { 0.0, 0.1488743389, 0.4333953941, 0.6794095682,
				0.8650633666, 0.9739065285 };
		double[] w = { 0.0, 0.2955242247, 0.2692667193, 0.2190863625,
				0.1494513491, 0.0666713443 };
		xm = 0.5 * (b + a);
		xr = 0.5 * (b - a);
		s = 0;
		// Will be twice the average value of the function, since the
		// ten weights (five numbers above each used twice)sum to 2.
		for (j = 1; j <= 5; j++) {
			dx = xr * x[j];// dx=xr*x[j];
			// s += w[j]*((*func)(xm+dx)+(*func)(xm-dx));
			double ff1 = 0.0;
			double ff2 = 0.0;
			if (fs.compareTo("f1") == 0) {
				ff1 = f1(xm + dx);
				ff2 = f1(xm - dx);
			} else if (fs.compareTo("f2") == 0) {
				ff1 = f2(xm + dx);
				ff2 = f2(xm - dx);
			} else if (fs.compareTo("f3") == 0) {
				ff1 = f3(xm + dx);
				ff2 = f3(xm - dx);
			}

			// s += w[j]*(func.F(xm+dx)+func.F(xm-dx));
			s += w[j] * (ff1 + ff2);
		}
		return s *= xr; // Scale the answer to the range of integration.

		/*
		 * int n=5; double[]
		 * w={0.23692688,0.47862868,0.56888889,0.47862868,0.23692688}; double[]
		 * x={-0.90617985,-0.53846931,0.0,0.53846931,0.90617985}; b=0.5*(b+a);
		 * a=b-a; double s=0.0; for (int i=1; i<=n; i++) {
		 * if(fs.compareTo("f1")==0) s=s+w[i-1]*f1(a*x[i-1]+b); else
		 * if(fs.compareTo("f2")==0) s=s+w[i-1]*f2(a*x[i-1]+b); else
		 * if(fs.compareTo("f3")==0) s=s+w[i-1]*f3(a*x[i-1]+b); } s=a*s; return
		 * s;
		 */
	}

	// =======================
	/**
	 * Another integration by QTrapez rule (between a and b).
	 * @param a a
	 * @param b b
	 * @return the result
	 */
	public double qTrapez(double a, double b) {
		double eps = 1e-8;
		int n = 1; // evaluare optimizata a pasului de integrare
		double h = b - a;
		double s = 0.5 * h * (func.F(a) + func.F(b));
		double del = 0.0;
		boolean zero = true;// control
		while (zero) {
			del = 0.0;
			for (int i = 1; i <= n; i++)
				del = del + func.F(a + (i - 0.5) * h);
			del = 0.5 * (h * del - s);
			s = s + del;
			h = 0.5 * h;
			n = 2 * n;// pana la 2^30!! maxvalue=2^31-1!!!
			if (Math.abs(del) <= eps * Math.abs(s)
					|| n >= (1. + Integer.MAX_VALUE) * 0.5)
				zero = false;
		}

		func.printSequence("============QTrapez evaluation===================");
		func.printSequence("Results= " + nf.format(s));
		func.printSequence("=================================================");

		return s;
	}

	/**
	 * Another integration by QSimpson rule (between a and b).
	 * @param a a
	 * @param b b
	 * @return the result
	 */
	public double qSimpson(double a, double b) {
		double eps = 1e-8;
		int n = 1; // evaluare optimizata a pasului de integrare
		double h = b - a;
		double t = 0.5 * h * (func.F(a) + func.F(b));
		double snew = 0.0;
		double s = t;
		double del = 0.0;
		boolean zero = true;// control
		while (zero) {
			del = 0.0;
			for (int i = 1; i <= n; i++)
				del = del + func.F(a + (i - 0.5) * h);
			del = 0.5 * (h * del - t);
			t = t + del;
			snew = t + del / 3;
			del = snew - s;
			s = snew;
			h = 0.5 * h;
			n = 2 * n;// pana la 2^30!! maxvalue=2^31-1!!!
			if (Math.abs(del) <= eps * Math.abs(s)
					|| n >= (1. + Integer.MAX_VALUE) * 0.5)
				zero = false;
		}

		func.printSequence("============QSimpson evaluation==================");
		func.printSequence("Results= " + nf.format(s));
		func.printSequence("=================================================");

		return s;
	}

	/**
	 * Another integration routine using QGauss rule (between a and b).
	 * @param a a
	 * @param b b
	 * @return the result
	 */
	public double qGauss5(double a, double b) {
		int n = 5;
		double[] w = { 0.23692688, 0.47862868, 0.56888889, 0.47862868,
				0.23692688 };
		double[] x = { -0.90617985, -0.53846931, 0.0, 0.53846931, 0.90617985 };
		b = 0.5 * (b + a);
		a = b - a;
		double s = 0.0;
		for (int i = 1; i <= n; i++)
			s = s + w[i - 1] * func.F(a * x[i - 1] + b);
		s = a * s;

		func.printSequence("============QGauss5 evaluation===================");
		func.printSequence("Results= " + nf.format(s));
		func.printSequence("=================================================");

		return s;
	}

	// DCADRE input data
	/**
	 * Integration using EGSnrc routine DCADRE. Integrates F(x) from A to B using 
	 * cautious adaptive Romberg extrapolation.
	 * @param A A
	 * @param B B
	 * @return the result
	 */
	public double QD(double A, double B) {
		double QD = 0.;
		first_time = true;
		double ADUM = A;
		double BDUM = B;
		double ERRDUM = 0.;
		int IER = 0;
		QD = DCADRE(ADUM, BDUM, $AERR, $RERR, ERRDUM, IER);
		IER = IERq;
		ERRDUM = ERRORq;
		if (IER > 66) {
			if (first_time) {
				first_time = false;
			}
		}

		// func.printSequence("==========advanced QD evaluation=================");
		// func.printSequence("Results= "+nf.format(QD));
		// func.printSequence("=================================================");

		return QD;
	}

	/**
	 * Used by QD
	 * @param A A
	 * @param B B
	 * @param AERR AERR
	 * @param RERR RERR
	 * @param ERROR ERROR
	 * @param IER IER
	 * @return the result
	 */
	public double DCADRE(double A, double B, double AERR, double RERR,
			double ERROR, int IER) {
		// "------------------------------------------------------------------"
		// "-DCADRE--------D-------LIBRARY 1----------------------------------"
		// "------------------------------------------------------------------"
		// "                                                                  "
		// "FUNCTION:          - INTEGRATE F(X) FROM A TO B, USING CAUTIOUS   "
		// "                     ADAPTIVE ROMBERG EXTRAPOLATION.              "
		// "                                                                  "
		// "USAGE:             - FUNCTION DCADRE(F,A,B,AERR,RERR,ERROR,IER)   "
		// "                                                                  "
		// "PARAMETERS: DCADRE - ESTIMATE OF THE INTEGRAL OF F(X) FROM A TO B."
		// "                                                                  "
		// "            F      - A SINGLE-ARGUMENT REAL FUNCTION SUBPROGRAM   "
		// "                     SUPPLIED BY THE USER.  F MUST BE DECLARED    "
		// "                     EXTERNAL IN THE CALLING PROGRAM.             "
		// "                                                                  "
		// "            A,B    - THE TWO ENDPOINTS OF THE INTERVAL OF         "
		// "                     INTEGRATION (INPUT).                         "
		// "                                                                  "
		// "            AERR   - DESIRED ABSOLUTE ERROR IN THE ANSWER (INPUT)."
		// "                                                                  "
		// "            RERR   - DESIRED RELATIVE ERROR IN THE ANSWER (INPUT)."
		// "                                                                  "
		// "            ERROR  - ESTIMATED BOUND ON THE ABSOLUTE ERROR OF     "
		// "                     THE OUTPUT NUMBER, DCADRE.                   "
		// "                                                                  "
		// "            IER    - ERROR PARAMETER                              "
		// "                                                                  "
		// "                     WARNING ERROR(WITH FIX) = 64 + N             "
		// "                                                                  "
		// "                       N = 1 IMPLIES THAT ONE OR MORE SINGULAR-   "
		// "                             ITIES WERE SUCCESSFULLY HANDLED.     "
		// "                                                                  "
		// "                       N = 2 IMPLIES THAT, IN SOME SUBINTERVAL(S),"
		// "                             THE ESTIMATE OF THE INTEGRAL WAS     "
		// "                             ACCEPTED MERELY BECAUSE THE ESTIMATED"
		// "                             ERROR WAS SMALL, EVEN THOUGH NO REG- "
		// "                             ULAR BEHAVIOR WAS RECOGNIZED.        "
		// "                                                                  "
		// "                     TERMINAL ERROR = 128 + N                     "
		// "                                                                  "
		// "                       N = 3 FAILURE DUE TO INSUFFICIENT INTERNAL "
		// "                             WORKING STORAGE.                     "
		// "                                                                  "
		// "                       N = 4 FAILURE.  THIS MAY BE DUE TO TOO MUCH"
		// "                             NOISE IN THE FUNCTION (RELATIVE TO   "
		// "                             THE GIVEN ERROR REQUIREMENTS) OR DUE "
		// "                             TO AN ILL-BEHAVED INTEGRAND.         "
		// "                                                                  "
		// "                       N = 5 INDICATES THAT RERR IS GREATER THAN  "
		// "                             0.1, OR RERR IS LESS THAN 0.0, OR    "
		// "                             RERR IS TOO SMALL FOR THE PRECISION  "
		// "                             OF THE MACHINE.                      "
		// "                                                                  "
		// "------------------------------------------------------------------"
		// "VERSION DATE:      - 8 OCTOBER 1974                               "
		// "                                                                  "
		// "MORTRAN VERSION    - 4 OCTOBER 1984/1545 (W. R. NELSON)           "
		// "------------------------------------------------------------------"

		Aq = A;
		Bq = B;
		AERRq = AERR;
		RERRq = RERR;
		ERRORq = ERROR;
		IERq = IER;
		// --------------
		resetDCADRE();// reset global values!!
		// test---------------------------------------
		// DCADRE = F(23);//WORKS
		// --------------------------------------------
		if (LENGTHq == ZEROq)// GO TO 215;
		{
			goto215q();
			return DCADREq;// 9005 RETURN;
		}
		if ((RERR > P1q) || (RERR < ZEROq))// GO TO 210;
		{
			IERq = 133;// 210 IER=133;
			goto215q();
			return DCADREq;// 9005 RETURN;
		}
		if ((AERR == ZEROq) && ((RERR + HUNq) <= HUNq)) // GO TO 210;
		{
			IERq = 133;// 210 IER=133;
			goto215q();
			return DCADREq;// 9005 RETURN;
		}

		RIGHTq = false;// 5
						// #########################################################RIGHT=.FALSE.;
		// " INVESTIGATION OF A PARTICULAR       "
		// " SUBINTERVAL BEGINS AT THIS POINT. "
		STEPq = ENDq - BEGq;// 10#################################################10
							// STEP=END - BEG;
		ASTEPq = DABS(STEPq);
		if (ASTEPq < STEPMNq)// GO TO 205;
		{
			IERq = 132;// 205 IER=132;
			goto215q();
			return DCADREq;
		}
		if ((STEPNMq + ASTEPq) == STEPNMq)// GO TO 205;
		{
			IERq = 132;// 205 IER=132;
			goto215q();
			return DCADREq;
		}
		Tq[0][0] = FBEGq + FENDq;// T(1,1)=FBEG + FEND;
		TABSq = DABS(FBEGq) + DABS(FENDq);
		Lq = 1;
		Nq = 1;
		H2CONVq = false;
		AITKENq = false;
		LM1q = Lq;// 15########################################################################
					// LM1=L;
		Lq = Lq + 1;
		// "  CALCULATE THE NEXT TRAPEZOID SUM,   "
		// "  T(L,1), WHICH IS BASED ON *N2* + 1"
		// "  EQUISPACED POINTS. HERE,          "
		// "  N2 = N*2 = 2**(L-1).              "
		N2q = Nq + Nq;
		FNq = N2q;
		ISTEPq = (IENDq - IBEGq) / Nq;
		if (ISTEPq <= 1) // IF(ISTEP>1) GO TO 25;
		{
			notgotto25q();
			if (exitq) {
				IERq = 131;// 200 IER=131;
				return DCADREq;
			}
		}
		ISTEP2q = IBEGq + ISTEPq / 2;// #25 ISTEP2=IBEG + ISTEP/2;
		SUMq = ZEROq;
		SUMABSq = ZEROq;
		for (Iq = ISTEP2q - 1; Iq < IENDq; Iq = Iq + ISTEPq) {
			SUMq = SUMq + TSq[Iq];
			SUMABSq = SUMABSq + DABS(TSq[Iq]);
		}
		Tq[Lq - 1][0] = Tq[Lq - 2][0] * HALFq + SUMq / FNq;// T(L,1)=T(L-1,1)*HALF+SUM/FN;
		TABSq = TABSq * HALFq + SUMABSq / FNq;
		ABSIq = ASTEPq * TABSq;
		Nq = N2q;
		// "     GET PRELIMINARY VALUE FOR *VINT*    "
		// "     FROM LAST TRAPEZOID SUM AND UPDATE"
		// "     THE ERROR REQUIREMENT *ERGOAL*    "
		// "     FOR THIS SUBINTERVAL.             "
		ITq = 1;
		VINTq = STEPq * Tq[Lq - 1][0];// STEP*T(L,1);
		TABTLMq = TABSq * TENq;
		FNSIZEq = DMAX1(FNSIZEq, DABS(Tq[Lq - 1][0]));// DMAX1(FNSIZE,DABS(T(L,1)));
		ERGLq = ASTEPq * FNSIZEq * TENq;
		ERGOALq = STAGEq * DMAX1(ERRAq, ERRRq * DABS(CURESTq + VINTq));
		// "  COMPLETE ROW L AND COLUMN L OF *T*  "
		// "    ARRAY.                            "
		FEXTRPq = ONEq;
		for (Iq = 0; Iq < LM1q; Iq++) {
			FEXTRPq = FEXTRPq * FOURq;
			Tq[Iq][Lq - 1] = Tq[Lq - 1][Iq] - Tq[Lq - 2][Iq];// T(I,L)=T(L,I) -
																// T(L-1,I);
			Tq[Lq - 1][Iq + 1] = Tq[Lq - 1][Iq] + Tq[Iq][Lq - 1]
					/ (FEXTRPq - ONEq);// T(L,I+1)=T(L,I) + T(I,L)/(FEXTRP-ONE);
		}
		ERRERq = ASTEPq * DABS(Tq[0][Lq - 1]);// ERRER=ASTEP*DABS(T(1,L));
		// "  PRELIMINARY DECISION PROCEDURE      "
		// "  IF L = 2 AND T(2,1) = T(1,1),     "
		// "  GO TO 135 TO FOLLOW UP THE        "
		// "  IMPRESSION THAT INTERGRAND IS     "
		// "  STRAIGHT LINE.                    "
		if (Lq > 2)// if (L>2)// GO TO 40;
		{
			goto40q();
			if (exitq) {
				return DCADREq;
			}
		}

		if ((TABSq + P1q * DABS(Tq[0][1])) == TABSq) // IF(TABS+P1*DABS(T(1,2)).EQ.TABS)
														// GO TO 135;
		{
			goto135q();
			if (exitq) {
				return DCADREq;
			}
		}
		// "                              CACULATE NEXT RATIOS FOR            "
		// "                                COLUMNS 1,...,L-2 OF T-TABLE      "
		// "                                RATIO IS SET TO ZERO IF DIFFERENCE"
		// "                                IN LAST TWO ENTRIES OF COLUMN IS  "
		// "                                ABOUT ZERO                        "
		goto15q();
		if (exitq)
			return DCADREq;

		return DCADREq;
	}

	/**
	 * Used by DCADRE
	 */
	private void goto5q() {
		RIGHTq = false;// 5
						// #########################################################RIGHT=.FALSE.;
		// " INVESTIGATION OF A PARTICULAR       "
		// " SUBINTERVAL BEGINS AT THIS POINT. "
		STEPq = ENDq - BEGq;// 10#################################################10
							// STEP=END - BEG;
		ASTEPq = DABS(STEPq);
		if (ASTEPq < STEPMNq)// GO TO 205;
		{
			IERq = 132;// 205 IER=132;
			goto215q();
			exitq = true;
			return;
		}
		if ((STEPNMq + ASTEPq) == STEPNMq)// GO TO 205;
		{
			IERq = 132;// 205 IER=132;
			goto215q();
			exitq = true;
			return;
		}
		Tq[0][0] = FBEGq + FENDq;// T(1,1)=FBEG + FEND;
		TABSq = DABS(FBEGq) + DABS(FENDq);
		Lq = 1;
		Nq = 1;
		H2CONVq = false;
		AITKENq = false;
		LM1q = Lq;// 15########################################################################
					// LM1=L;
		Lq = Lq + 1;
		// "  CALCULATE THE NEXT TRAPEZOID SUM,   "
		// "  T(L,1), WHICH IS BASED ON *N2* + 1"
		// "  EQUISPACED POINTS. HERE,          "
		// "  N2 = N*2 = 2**(L-1).              "
		N2q = Nq + Nq;
		FNq = N2q;
		ISTEPq = (IENDq - IBEGq) / Nq;
		if (ISTEPq <= 1) // IF(ISTEP>1) GO TO 25;
		{
			notgotto25q();
			if (exitq) {
				IERq = 131;// 200 IER=131;
				return;
			}
		}
		ISTEP2q = IBEGq + ISTEPq / 2;// #25 ISTEP2=IBEG + ISTEP/2;
		SUMq = ZEROq;
		SUMABSq = ZEROq;
		for (Iq = ISTEP2q - 1; Iq < IENDq; Iq = Iq + ISTEPq) {
			SUMq = SUMq + TSq[Iq];
			SUMABSq = SUMABSq + DABS(TSq[Iq]);
		}
		Tq[Lq - 1][0] = Tq[Lq - 2][0] * HALFq + SUMq / FNq;// T(L,1)=T(L-1,1)*HALF+SUM/FN;
		TABSq = TABSq * HALFq + SUMABSq / FNq;
		ABSIq = ASTEPq * TABSq;
		Nq = N2q;
		// "     GET PRELIMINARY VALUE FOR *VINT*    "
		// "     FROM LAST TRAPEZOID SUM AND UPDATE"
		// "     THE ERROR REQUIREMENT *ERGOAL*    "
		// "     FOR THIS SUBINTERVAL.             "
		ITq = 1;
		VINTq = STEPq * Tq[Lq - 1][0];// STEP*T(L,1);
		TABTLMq = TABSq * TENq;
		FNSIZEq = DMAX1(FNSIZEq, DABS(Tq[Lq - 1][0]));// DMAX1(FNSIZE,DABS(T(L,1)));
		ERGLq = ASTEPq * FNSIZEq * TENq;
		ERGOALq = STAGEq * DMAX1(ERRAq, ERRRq * DABS(CURESTq + VINTq));
		// "  COMPLETE ROW L AND COLUMN L OF *T*  "
		// "    ARRAY.                            "
		FEXTRPq = ONEq;
		for (Iq = 0; Iq < LM1q; Iq++) {
			FEXTRPq = FEXTRPq * FOURq;
			Tq[Iq][Lq - 1] = Tq[Lq - 1][Iq] - Tq[Lq - 2][Iq];// T(I,L)=T(L,I) -
																// T(L-1,I);
			Tq[Lq - 1][Iq + 1] = Tq[Lq - 1][Iq] + Tq[Iq][Lq - 1]
					/ (FEXTRPq - ONEq);// T(L,I+1)=T(L,I) + T(I,L)/(FEXTRP-ONE);
		}
		ERRERq = ASTEPq * DABS(Tq[0][Lq - 1]);// ERRER=ASTEP*DABS(T(1,L));
		// "  PRELIMINARY DECISION PROCEDURE      "
		// "  IF L = 2 AND T(2,1) = T(1,1),     "
		// "  GO TO 135 TO FOLLOW UP THE        "
		// "  IMPRESSION THAT INTERGRAND IS     "
		// "  STRAIGHT LINE.                    "
		if (Lq > 2)// if (L>2)// GO TO 40;
		{
			goto40q();
			if (exitq) {
				return;
			}
		}

		if ((TABSq + P1q * DABS(Tq[0][1])) == TABSq) // IF(TABS+P1*DABS(T(1,2)).EQ.TABS)
														// GO TO 135;
		{
			goto135q();
			if (exitq) {
				return;
			}
		}
		goto15q();
		if (exitq)
			return;
	}

	/**
	 * Used by DCADRE
	 */
	private void goto10q() {
		STEPq = ENDq - BEGq;// 10#################################################10
							// STEP=END - BEG;
		ASTEPq = DABS(STEPq);
		if (ASTEPq < STEPMNq)// GO TO 205;
		{
			IERq = 132;// 205 IER=132;
			goto215q();
			exitq = true;
			return;
		}
		if ((STEPNMq + ASTEPq) == STEPNMq)// GO TO 205;
		{
			IERq = 132;// 205 IER=132;
			goto215q();
			exitq = true;
			return;
		}
		Tq[0][0] = FBEGq + FENDq;// T(1,1)=FBEG + FEND;
		TABSq = DABS(FBEGq) + DABS(FENDq);
		Lq = 1;
		Nq = 1;
		H2CONVq = false;
		AITKENq = false;
		LM1q = Lq;// 15########################################################################
					// LM1=L;
		Lq = Lq + 1;
		// "  CALCULATE THE NEXT TRAPEZOID SUM,   "
		// "  T(L,1), WHICH IS BASED ON *N2* + 1"
		// "  EQUISPACED POINTS. HERE,          "
		// "  N2 = N*2 = 2**(L-1).              "
		N2q = Nq + Nq;
		FNq = N2q;
		ISTEPq = (IENDq - IBEGq) / Nq;
		if (ISTEPq <= 1) // IF(ISTEP>1) GO TO 25;
		{
			notgotto25q();
			if (exitq) {
				IERq = 131;// 200 IER=131;
				return;
			}
		}
		ISTEP2q = IBEGq + ISTEPq / 2;// #25 ISTEP2=IBEG + ISTEP/2;
		SUMq = ZEROq;
		SUMABSq = ZEROq;
		for (Iq = ISTEP2q - 1; Iq < IENDq; Iq = Iq + ISTEPq) {
			SUMq = SUMq + TSq[Iq];
			SUMABSq = SUMABSq + DABS(TSq[Iq]);
		}
		Tq[Lq - 1][0] = Tq[Lq - 2][0] * HALFq + SUMq / FNq;// T(L,1)=T(L-1,1)*HALF+SUM/FN;
		TABSq = TABSq * HALFq + SUMABSq / FNq;
		ABSIq = ASTEPq * TABSq;
		Nq = N2q;
		// "     GET PRELIMINARY VALUE FOR *VINT*    "
		// "     FROM LAST TRAPEZOID SUM AND UPDATE"
		// "     THE ERROR REQUIREMENT *ERGOAL*    "
		// "     FOR THIS SUBINTERVAL.             "
		ITq = 1;
		VINTq = STEPq * Tq[Lq - 1][0];// STEP*T(L,1);
		TABTLMq = TABSq * TENq;
		FNSIZEq = DMAX1(FNSIZEq, DABS(Tq[Lq - 1][0]));// DMAX1(FNSIZE,DABS(T(L,1)));
		ERGLq = ASTEPq * FNSIZEq * TENq;
		ERGOALq = STAGEq * DMAX1(ERRAq, ERRRq * DABS(CURESTq + VINTq));
		// "  COMPLETE ROW L AND COLUMN L OF *T*  "
		// "    ARRAY.                            "
		FEXTRPq = ONEq;
		for (Iq = 0; Iq < LM1q; Iq++) {
			FEXTRPq = FEXTRPq * FOURq;
			Tq[Iq][Lq - 1] = Tq[Lq - 1][Iq] - Tq[Lq - 2][Iq];// T(I,L)=T(L,I) -
																// T(L-1,I);
			Tq[Lq - 1][Iq + 1] = Tq[Lq - 1][Iq] + Tq[Iq][Lq - 1]
					/ (FEXTRPq - ONEq);// T(L,I+1)=T(L,I) + T(I,L)/(FEXTRP-ONE);
		}
		ERRERq = ASTEPq * DABS(Tq[0][Lq - 1]);// ERRER=ASTEP*DABS(T(1,L));
		// "  PRELIMINARY DECISION PROCEDURE      "
		// "  IF L = 2 AND T(2,1) = T(1,1),     "
		// "  GO TO 135 TO FOLLOW UP THE        "
		// "  IMPRESSION THAT INTERGRAND IS     "
		// "  STRAIGHT LINE.                    "
		if (Lq > 2)// if (L>2)// GO TO 40;
		{
			goto40q();
			if (exitq) {
				return;
			}
		}

		if ((TABSq + P1q * DABS(Tq[0][1])) == TABSq) // IF(TABS+P1*DABS(T(1,2)).EQ.TABS)
														// GO TO 135;
		{
			goto135q();
			if (exitq) {
				return;
			}
		}
		goto15q();
		if (exitq)
			return;
	}

	/**
	 * Used by DCADRE
	 */
	private void goto15q() {
		LM1q = Lq;// 15########################################################################
					// LM1=L;
		Lq = Lq + 1;
		// "  CALCULATE THE NEXT TRAPEZOID SUM,   "
		// "  T(L,1), WHICH IS BASED ON *N2* + 1"
		// "  EQUISPACED POINTS. HERE,          "
		// "  N2 = N*2 = 2**(L-1).              "
		N2q = Nq + Nq;
		FNq = N2q;
		ISTEPq = (IENDq - IBEGq) / Nq;
		if (ISTEPq <= 1) // IF(ISTEP>1) GO TO 25;
		{
			notgotto25q();
			if (exitq) {
				IERq = 131;// 200 IER=131;
				return;
			}
		}
		ISTEP2q = IBEGq + ISTEPq / 2;// #25 ISTEP2=IBEG + ISTEP/2;
		SUMq = ZEROq;
		SUMABSq = ZEROq;
		for (Iq = ISTEP2q - 1; Iq < IENDq; Iq = Iq + ISTEPq) {
			SUMq = SUMq + TSq[Iq];
			SUMABSq = SUMABSq + DABS(TSq[Iq]);
		}
		Tq[Lq - 1][0] = Tq[Lq - 2][0] * HALFq + SUMq / FNq;// T(L,1)=T(L-1,1)*HALF+SUM/FN;
		TABSq = TABSq * HALFq + SUMABSq / FNq;
		ABSIq = ASTEPq * TABSq;
		Nq = N2q;
		// "     GET PRELIMINARY VALUE FOR *VINT*    "
		// "     FROM LAST TRAPEZOID SUM AND UPDATE"
		// "     THE ERROR REQUIREMENT *ERGOAL*    "
		// "     FOR THIS SUBINTERVAL.             "
		ITq = 1;
		VINTq = STEPq * Tq[Lq - 1][0];// STEP*T(L,1);
		TABTLMq = TABSq * TENq;
		FNSIZEq = DMAX1(FNSIZEq, DABS(Tq[Lq - 1][0]));// DMAX1(FNSIZE,DABS(T(L,1)));
		ERGLq = ASTEPq * FNSIZEq * TENq;
		ERGOALq = STAGEq * DMAX1(ERRAq, ERRRq * DABS(CURESTq + VINTq));
		// "  COMPLETE ROW L AND COLUMN L OF *T*  "
		// "    ARRAY.                            "
		FEXTRPq = ONEq;
		for (Iq = 0; Iq < LM1q; Iq++) {
			FEXTRPq = FEXTRPq * FOURq;
			Tq[Iq][Lq - 1] = Tq[Lq - 1][Iq] - Tq[Lq - 2][Iq];// T(I,L)=T(L,I) -
																// T(L-1,I);
			Tq[Lq - 1][Iq + 1] = Tq[Lq - 1][Iq] + Tq[Iq][Lq - 1]
					/ (FEXTRPq - ONEq);// T(L,I+1)=T(L,I) + T(I,L)/(FEXTRP-ONE);
		}
		ERRERq = ASTEPq * DABS(Tq[0][Lq - 1]);// ERRER=ASTEP*DABS(T(1,L));
		// "  PRELIMINARY DECISION PROCEDURE      "
		// "  IF L = 2 AND T(2,1) = T(1,1),     "
		// "  GO TO 135 TO FOLLOW UP THE        "
		// "  IMPRESSION THAT INTERGRAND IS     "
		// "  STRAIGHT LINE.                    "
		if (Lq > 2)// if (L>2)// GO TO 40;
		{
			goto40q();
			if (exitq) {
				return;
			}
		}

		if ((TABSq + P1q * DABS(Tq[0][1])) == TABSq) // IF(TABS+P1*DABS(T(1,2)).EQ.TABS)
														// GO TO 135;
		{
			goto135q();
			if (exitq) {
				return;
			}
		}
		goto15q();
		if (exitq)
			return;
	}

	/**
	 * Used by DCADRE
	 */
	private void notgotto25q() {
		IIq = IENDq;
		IENDq = IENDq + Nq;
		if (IENDq > MAXTSq)// IF(IEND.GT.MAXTS) GO TO 200;
		{
			goto215q();
			exitq = true;// return DCADREq;
			return;
		}
		HOVNq = STEPq / FNq;
		IIIq = IENDq;
		FIq = ONEq;
		for (Iq = 0; Iq < N2q; Iq = Iq + 2) {
			TSq[IIIq - 1] = TSq[IIq - 1];// TS[III]=TS(II);
			RVALq = ENDq - FIq * HOVNq;
			TSq[IIIq - 2] = func.F(RVALq);// TS(III-1)=F(RVAL);@@@
			FIq = FIq + TWOq;
			IIIq = IIIq - 2;
			IIq = IIq - 1;
		}
		ISTEPq = 2;
	}

	/**
	 * Used by DCADRE
	 */
	private void goto40q() {
		for (Iq = 1; Iq < LM1q; Iq++)// 40 DO 45 I=2,LM1;pana la 45
		{
			DIFFq = ZEROq;
			if (TABTLMq + DABS(Tq[Iq - 1][Lq - 1]) != TABTLMq)
				DIFFq = Tq[Iq - 1][LM1q - 1] / Tq[Iq - 1][Lq - 1];
			// IF(TABTLM+DABS(T(I-1,L)).NE.TABTLM) DIFF=T(I-1,LM1)/T(I-1,L);
			Tq[Iq - 1][LM1q - 1] = DIFFq;// T(I-1,LM1)=DIFF;
		}
		// 45 CONTINUE;
		// if (DABS(FOURq-Tq(1,LM1)).LE.H2TOL) GO TO 60;
		if (DABS(FOURq - Tq[0][LM1q - 1]) <= H2TOLq) {
			goto60q();
			if (exitq)
				return;
		}
		// IF(T(1,LM1).EQ.ZERO) GO TO 55;
		if (Tq[0][LM1q - 1] == ZEROq) {
			goto55q();
			if (exitq)
				return;
		}
		// IF(DABS(TWO-DABS(T(1,LM1))).LT.JUMPTL) GO TO 130;
		if (DABS(TWOq - DABS(Tq[0][LM1q - 1])) < JUMPTLq) {
			goto130q();
			if (exitq)
				return;
		}
		// IF(L.EQ.3) GO TO 15;
		if (Lq == 3) {
			goto15q();
			if (exitq)
				return;
		}
		H2CONVq = false;
		// IF(DABS((T(1,LM1)-T(1,L-2))/T(1,LM1)).LE.AITTOL) GO TO 75;
		if (DABS((Tq[0][LM1q - 1] - Tq[0][Lq - 3]) / Tq[0][LM1q - 1]) <= AITTOLq) {
			goto75q();
			if (exitq)
				return;
		}
		// 50 IF(REGLAR) GO TO 55;
		if (REGLARq) {
			goto55q();
			if (exitq)
				return;
		}
		// IF(L.EQ.4) GO TO 15;
		if (Lq == 4) {
			goto15q();
			if (exitq)
				return;
		}
		// 55 IF(ERRER.GT.ERGOAL.AND.(ERGL+ERRER).NE.ERGL) GO TO 175;
		if ((ERRERq > ERGOALq) && (ERGLq + ERRERq) != ERGLq) {
			goto175q();
			if (exitq)
				return;
		}
		goto145q();
		if (exitq)
			return;
	}

	/**
	 * Used by DCADRE
	 */
	private void goto55q() {
		// 55 IF(ERRER.GT.ERGOAL.AND.(ERGL+ERRER).NE.ERGL) GO TO 175;
		if ((ERRERq > ERGOALq) && ((ERGLq + ERRERq) != ERGLq)) {
			goto175q();
			if (exitq)
				return;
		}
		goto145q();
		if (exitq)
			return;
	}

	/**
	 * Used by DCADRE
	 */
	private void goto60q() {
		// 60 IF(H2CONV) GO TO 65;
		if (H2CONVq) {
			goto65q();
			if (exitq)
				return;
		}
		AITKENq = false;
		H2CONVq = true;
		// 65 FEXTRP=FOUR;
		FEXTRPq = FOURq;
		// 70 IT=IT + 1;
		ITq = ITq + 1;
		// VINTq=STEPq*Tq(L,IT);
		VINTq = STEPq * Tq[Lq - 1][ITq - 1];
		// ERRERq=DABS(STEPq/(FEXTRPq-ONEq)*T(IT-1,L));
		ERRERq = DABS(STEPq / (FEXTRPq - ONEq) * Tq[ITq - 2][Lq - 1]);
		// IF(ERRER.LE.ERGOAL) GO TO 160;
		if (ERRERq <= ERGOALq) {
			goto160q();
			if (exitq)
				return;
		}
		// IF(ERGL+ERRER.EQ.ERGL) GO TO 160;
		if ((ERGLq + ERRERq) == ERGLq) {
			goto160q();
			if (exitq)
				return;
		}
		// IF(IT.EQ.LM1) GO TO 125;
		if (ITq == LM1q) {
			goto125q();
			if (exitq)
				return;
		}
		// IF(T(IT,LM1).EQ.ZERO) GO TO 70;
		if (Tq[ITq - 1][LM1q - 1] == ZEROq) {
			goto70q();
			if (exitq)
				return;
		}
		// IF(T(IT,LM1).LE.FEXTRP) GO TO 125;
		if (Tq[ITq - 1][LM1q - 1] <= FEXTRPq) {
			goto125q();
			if (exitq)
				return;
		}
		// IF(DABS(T(IT,LM1)/FOUR-FEXTRP)/FEXTRP.LT.AITTOL)
		if (DABS(Tq[ITq - 1][LM1q - 1] / FOURq - FEXTRPq) / FEXTRPq < AITTOLq)
			FEXTRPq = FEXTRPq * FOURq;
		goto70q();
		if (exitq)
			return;
	}

	/**
	 * Used by DCADRE
	 */
	private void goto65q() {
		// 65 FEXTRP=FOUR;
		FEXTRPq = FOURq;
		// 70 IT=IT + 1;
		ITq = ITq + 1;
		// VINTq=STEPq*Tq(L,IT);
		VINTq = STEPq * Tq[Lq - 1][ITq - 1];
		// ERRERq=DABS(STEPq/(FEXTRPq-ONEq)*T(IT-1,L));
		ERRERq = DABS(STEPq / (FEXTRPq - ONEq) * Tq[ITq - 2][Lq - 1]);
		// IF(ERRER.LE.ERGOAL) GO TO 160;
		if (ERRERq <= ERGOALq) {
			goto160q();
			if (exitq)
				return;
		}
		// IF(ERGL+ERRER.EQ.ERGL) GO TO 160;
		if ((ERGLq + ERRERq) == ERGLq) {
			goto160q();
			if (exitq)
				return;
		}
		// IF(IT.EQ.LM1) GO TO 125;
		if (ITq == LM1q) {
			goto125q();
			if (exitq)
				return;
		}
		// IF(T(IT,LM1).EQ.ZERO) GO TO 70;
		if (Tq[ITq - 1][LM1q - 1] == ZEROq) {
			goto70q();
			if (exitq)
				return;
		}
		// IF(T(IT,LM1).LE.FEXTRP) GO TO 125;
		if (Tq[ITq - 1][LM1q - 1] <= FEXTRPq) {
			goto125q();
			if (exitq)
				return;
		}
		// IF(DABS(T(IT,LM1)/FOUR-FEXTRP)/FEXTRP.LT.AITTOL)
		if (DABS(Tq[ITq - 1][LM1q - 1] / FOURq - FEXTRPq) / FEXTRPq < AITTOLq)
			FEXTRPq = FEXTRPq * FOURq;
		goto70q();
		if (exitq)
			return;
	}

	/**
	 * Used by DCADRE
	 */
	private void goto70q() {
		// 70 IT=IT + 1;
		ITq = ITq + 1;
		// VINTq=STEPq*Tq(L,IT);
		VINTq = STEPq * Tq[Lq - 1][ITq - 1];
		// ERRERq=DABS(STEPq/(FEXTRPq-ONEq)*T(IT-1,L));
		ERRERq = DABS(STEPq / (FEXTRPq - ONEq) * Tq[ITq - 2][Lq - 1]);
		// IF(ERRER.LE.ERGOAL) GO TO 160;
		if (ERRERq <= ERGOALq) {
			goto160q();
			if (exitq)
				return;
		}
		// IF(ERGL+ERRER.EQ.ERGL) GO TO 160;
		if ((ERGLq + ERRERq) == ERGLq) {
			goto160q();
			if (exitq)
				return;
		}
		// IF(IT.EQ.LM1) GO TO 125;
		if (ITq == LM1q) {
			goto125q();
			if (exitq)
				return;
		}
		// IF(T(IT,LM1).EQ.ZERO) GO TO 70;
		if (Tq[ITq - 1][LM1q - 1] == ZEROq) {
			goto70q();
			if (exitq)
				return;
		}
		// IF(T(IT,LM1).LE.FEXTRP) GO TO 125;
		if (Tq[ITq - 1][LM1q - 1] <= FEXTRPq) {
			goto125q();
			if (exitq)
				return;
		}
		// IF(DABS(T(IT,LM1)/FOUR-FEXTRP)/FEXTRP.LT.AITTOL)
		if (DABS(Tq[ITq - 1][LM1q - 1] / FOURq - FEXTRPq) / FEXTRPq < AITTOLq)
			FEXTRPq = FEXTRPq * FOURq;
		goto70q();
		if (exitq)
			return;
	}

	/**
	 * Used by DCADRE
	 */
	private void goto75q() {
		// "                              INTEGRAND MAY HAVE X**ALPHA TYPE    "
		// "                                SINGULARITY                       "
		// "                                RESULTING IN A RATIO OF *SING*  = "
		// "                                2**(ALPHA + 1)                    "
		// IF(T(1,LM1).LT.AITLOW) GO TO 175;
		if (Tq[0][LM1q - 1] < AITLOWq) {
			goto175q();
			if (exitq)
				return;
		}
		if (AITKENq) {
			goto80q();
			if (exitq)
				return;
		}
		H2CONVq = false;
		AITKENq = true;
		FEXTRPq = Tq[Lq - 3][LM1q - 1];// 80 FEXTRP=T(L-2,LM1);
		if (FEXTRPq > FOURP5q) {
			goto65q();
			if (exitq)
				return;
		}
		if (FEXTRPq < AITLOWq) {
			goto175q();
			if (exitq)
				return;
		}
		// if(DABS(FEXTRP-T(L-3,LM1))/T(1,LM1).GT.H2TOL) GO TO 175;
		if (DABS(FEXTRPq - Tq[Lq - 4][LM1q - 1]) / Tq[0][LM1q - 1] > H2TOLq) {
			goto175q();
			if (exitq)
				return;
		}
		SINGq = FEXTRPq;
		FEXTM1q = ONEq / (FEXTRPq - ONEq);
		AITq[0] = ZEROq;// AIT(1)=ZERO;
		for (Iq = 1; Iq < Lq; Iq++)// DO 85 I=2,L;
		{
			// AIT(I)=T(I,1) + (T(I,1)-T(I-1,1))*FEXTM1;
			AITq[Iq] = Tq[Iq][0] + (Tq[Iq][0] - Tq[Iq - 1][0]) * FEXTM1q;
			// Rq[Iq]=Tq(1,I-1);
			Rq[Iq] = Tq[0][Iq - 1];
			// DIF(I)=AIT(I) - AIT(I-1);
			DIFq[Iq] = AITq[Iq] - AITq[Iq - 1];
		}
		// 85 CONTINUE;
		ITq = 2;
		VINTq = STEPq * AITq[Lq - 1];// 90 VINT=STEP*AIT(L);
		ERRERq = ERRERq * FEXTM1q;
		// IF(ERRER.GT.ERGOAL.AND.(ERGL+ERRER).NE.ERGL) GO TO 95;
		if ((ERRERq > ERGOALq) && ((ERGLq + ERRERq) != ERGLq)) {
			goto95q();
			if (exitq)
				return;
		}
		ALPHAq = DLOG10(SINGq) / ALG4O2q - ONEq;
		IERq = Math.max(IERq, 65);
		goto160q();
		if (exitq)
			return;
	}

	/**
	 * Used by DCADRE
	 */
	private void goto80q() {
		// 80 FEXTRP=T(L-2,LM1);
		FEXTRPq = Tq[Lq - 3][LM1q - 1];
		if (FEXTRPq > FOURP5q) {
			goto65q();
			if (exitq)
				return;
		}
		// IF(FEXTRP.LT.AITLOW) GO TO 175;
		if (FEXTRPq < AITLOWq) {
			goto175q();
			if (exitq)
				return;
		}
		// IF(DABS(FEXTRP-T(L-3,LM1))/T(1,LM1).GT.H2TOL) GO TO 175;
		if (DABS(FEXTRPq - Tq[Lq - 4][LM1q - 1]) / Tq[0][LM1q - 1] > H2TOLq) {
			goto175q();
			if (exitq)
				return;
		}
		SINGq = FEXTRPq;
		FEXTM1q = ONEq / (FEXTRPq - ONEq);
		// AIT(1)=ZERO;
		AITq[0] = ZEROq;
		for (Iq = 1; Iq < Lq; Iq++)// DO 85 I=2,L;
		{
			AITq[Iq] = Tq[Iq][0] + (Tq[Iq][0] - Tq[Iq - 1][0]) * FEXTM1q;
			// AIT(I)=T(I,1) + (T(I,1)-T(I-1,1))*FEXTM1;
			Rq[Iq] = Tq[0][Iq - 1];
			// R(I)=T(1,I-1);
			DIFq[Iq] = AITq[Iq] - AITq[Iq - 1];
			// DIF(I)=AIT(I) - AIT(I-1);
		}
		// 85 CONTINUE;
		ITq = 2;
		// 90 VINT=STEP*AIT(L);
		VINTq = STEPq * AITq[Lq - 1];
		ERRERq = ERRERq * FEXTM1q;
		// IF(ERRER.GT.ERGOAL.AND.(ERGL+ERRER).NE.ERGL) GO TO 95;
		if ((ERRERq > ERGOALq) && ((ERGLq + ERRERq) != ERGLq)) {
			goto95q();
			if (exitq)
				return;
		}
		ALPHAq = DLOG10(SINGq) / ALG4O2q - ONEq;
		IERq = Math.max(IERq, 65);
		goto160q();
		if (exitq)
			return;
	}

	/**
	 * Used by DCADRE
	 */
	private void goto90q() {
		// 90 VINT=STEP*AIT(L);
		VINTq = STEPq * AITq[Lq - 1];
		ERRERq = ERRERq * FEXTM1q;
		// IF(ERRER.GT.ERGOAL.AND.(ERGL+ERRER).NE.ERGL) GO TO 95;
		if ((ERRERq > ERGOALq) && ((ERGLq + ERRERq) != ERGLq)) {
			goto95q();
			if (exitq)
				return;
		}
		ALPHAq = DLOG10(SINGq) / ALG4O2q - ONEq;
		IERq = Math.max(IERq, 65);
		goto160q();
		if (exitq)
			return;
	}

	/**
	 * Used by DCADRE
	 */
	private void goto95q() {
		ITq = ITq + 1;// 95 IT=IT + 1;
		// IF(IT.EQ.LM1) GO TO 125;
		if (ITq == LM1q) {
			goto125q();
			if (exitq)
				return;
		}
		if (ITq > 3) {
			goto100q();
			if (exitq)
				return;
		}
		H2NXTq = FOURq;
		SINGNXq = SINGq + SINGq;
		// 100 IF(H2NXT.LT.SINGNX) GO TO 105;
		if (H2NXTq < SINGNXq) {
			goto105q();
			if (exitq)
				return;
		}
		FEXTRPq = SINGNXq;
		SINGNXq = SINGNXq + SINGNXq;
		goto110q();
		if (exitq)
			return;
	}

	/**
	 * Used by DCADRE
	 */
	private void goto100q() {
		// 100 IF(H2NXT.LT.SINGNX) GO TO 105;
		if (H2NXTq < SINGNXq) {
			goto105q();
			if (exitq)
				return;
		}
		FEXTRPq = SINGNXq;
		SINGNXq = SINGNXq + SINGNXq;
		goto110q();
		if (exitq)
			return;
	}

	/**
	 * Used by DCADRE
	 */
	private void goto105q() {
		FEXTRPq = H2NXTq;
		H2NXTq = FOURq * H2NXTq;
		// 110 DO 115 I=IT,LM1;
		for (Iq = ITq - 1; Iq < LM1q; Iq++) {
			Rq[Iq + 1] = ZEROq;
			// IF(TABTLM+DABS(DIF(I+1)).NE.TABTLM) R(I+1)=DIF(I)/DIF(I+1);
			if (TABTLMq + DABS(DIFq[Iq + 1]) != TABTLMq) {
				Rq[Iq + 1] = DIFq[Iq] / DIFq[Iq + 1];
			}
		}
		// 115 CONTINUE;
		H2TFEXq = -H2TOLq * FEXTRPq;
		// IF(R(L)-FEXTRP.LT.H2TFEX) GO TO 125;
		if ((Rq[Lq - 1] - FEXTRPq) < H2TFEXq) {
			goto125q();
			if (exitq)
				return;
		}
		// IF(R(L-1)-FEXTRP.LT.H2TFEX) GO TO 125;
		if (Rq[Lq - 2] - FEXTRPq < H2TFEXq) {
			goto125q();
			if (exitq)
				return;
		}
		ERRERq = ASTEPq * DABS(DIFq[Lq - 1]);
		FEXTM1q = ONEq / (FEXTRPq - ONEq);
		// DO 120 I=IT,L;
		for (Iq = ITq - 1; Iq < Lq; Iq++) {
			// AIT(I)=AIT(I) + DIF(I)*FEXTM1;
			AITq[Iq] = AITq[Iq] + DIFq[Iq] * FEXTM1q;
			// DIF(I)=AIT(I) - AIT(I-1);
			DIFq[Iq] = AITq[Iq] - AITq[Iq - 1];
		}
		// 120 CONTINUE;
		goto90q();
		if (exitq)
			return;
	}

	/**
	 * Used by DCADRE
	 */
	private void goto110q() {
		// 110 DO 115 I=IT,LM1;
		for (Iq = ITq - 1; Iq < LM1q; Iq++) {
			Rq[Iq + 1] = ZEROq;
			// IF(TABTLM+DABS(DIF(I+1)).NE.TABTLM) R(I+1)=DIF(I)/DIF(I+1);
			if (TABTLMq + DABS(DIFq[Iq + 1]) != TABTLMq) {
				Rq[Iq + 1] = DIFq[Iq] / DIFq[Iq + 1];
			}
		}
		// 115 CONTINUE;
		H2TFEXq = -H2TOLq * FEXTRPq;
		// IF(R(L)-FEXTRP.LT.H2TFEX) GO TO 125;
		if ((Rq[Lq - 1] - FEXTRPq) < H2TFEXq) {
			goto125q();
			if (exitq)
				return;
		}
		// IF(R(L-1)-FEXTRP.LT.H2TFEX) GO TO 125;
		if (Rq[Lq - 2] - FEXTRPq < H2TFEXq) {
			goto125q();
			if (exitq)
				return;
		}
		ERRERq = ASTEPq * DABS(DIFq[Lq - 1]);
		FEXTM1q = ONEq / (FEXTRPq - ONEq);
		// DO 120 I=IT,L;
		for (Iq = ITq - 1; Iq < Lq; Iq++) {
			// AIT(I)=AIT(I) + DIF(I)*FEXTM1;
			AITq[Iq] = AITq[Iq] + DIFq[Iq] * FEXTM1q;
			// DIF(I)=AIT(I) - AIT(I-1);
			DIFq[Iq] = AITq[Iq] - AITq[Iq - 1];
		}
		// 120 CONTINUE;
		goto90q();
		if (exitq)
			return;
	}

	/**
	 * Used by DCADRE
	 */
	private void goto125q() {
		// "                              CURRENT TRAPEZOID SUM AND RESULTING "
		// "                                EXTRAPOLATED VALUES DID NOT GIVE  "
		// "                                A SMALL ENOUGH *ERRER*.           "
		// "                                NOTE -- HAVING PREVER .LT. ERRER  "
		// "                                IS AN ALMOST CERTAIN SIGN OF      "
		// "                                BEGINNING TROUBLE WITH IN THE FUNC"
		// "                                TION VALUES. HENCE, A WATCH FOR,  "
		// "                                AND CONTROL OF, NOISE SHOULD      "
		// "                                BEGIN HERE.                       "

		// 125 FEXTRP=DMAX1(PREVER/ERRER,AITLOW);
		FEXTRPq = DMAX1(PREVERq / ERRERq, AITLOWq);
		PREVERq = ERRERq;
		// IF(L.LT.5) GO TO 15;
		if (Lq < 5) {
			goto15q();
			if (exitq)
				return;
		}
		// IF(L-IT.GT.2.AND.ISTAGE.LT.MXSTGE) GO TO 170;
		if (((Lq - ITq) > 2) && (ISTAGEq < MXSTGEq)) {
			goto170q();
			if (exitq)
				return;
		}
		ERRETq = ERRERq / Math.pow(FEXTRPq, (MAXTBLq - Lq));
		// IF(ERRET.GT.ERGOAL.AND.(ERGL+ERRET).NE.ERGL) GO TO 170;
		if ((ERRETq > ERGOALq) && ((ERGLq + ERRETq) != ERGLq)) {
			goto170q();
			if (exitq)
				return;
		}
		goto15q();
		if (exitq)
			return;

	}

	/**
	 * Used by DCADRE
	 */
	private void goto130q() {
		// 130 IF(ERRER.GT.ERGOAL.AND.(ERGL+ERRER).NE.ERGL) GO TO 170;
		// "   NOTE THAT  2*FN=2**L              "
		if ((ERRERq > ERGOALq) && ((ERGLq + ERRERq) != ERGLq)) {
			goto170q();
			if (exitq)
				return;

		}
		DIFFq = DABS(Tq[0][Lq - 1]) * (FNq + FNq);// DIFFq=DABS(T(1,L))*(FN+FN);
		goto160q();
		if (exitq)
			return;
	}

	/**
	 * Used by DCADRE
	 */
	private void goto135q() {
		// "                              INTEGRAND IS STRAIGHT LINE          "
		// "                                TEST THIS ASSUMPTION BY COMPARING "
		// "                                THE VALUE OF THE INTEGRAND AT     "
		// "                                FOUR *RANDOMLY CHOSEN* POINTS WITH"
		// "                                THE VALUE OF THE STRAIGHT LINE    "
		// "                                INTERPOLATING THE INTEGRAND AT THE"
		// "                                TWO END POINTS OF THE SUB-INTERVAL"
		// "                                IF TEST IS PASSED, ACCEPT *VINT*  "

		SLOPEq = (FENDq - FBEGq) * TWOq;// 135 SLOPE=(FEND-FBEG)*TWO;
		FBEG2q = FBEGq + FBEGq;
		for (Iq = 0; Iq < 4; Iq++)// DO 140 I=1,4;Fa pana la 140!!
		{
			RVALq = BEGq + RNq[Iq] * STEPq;
			DIFFq = DABS(func.F(RVALq) - FBEG2q - RNq[Iq] * SLOPEq);
			if ((TABTLMq + DIFFq) != TABTLMq)// GO TO 155;
			{
				goto155q();
				if (exitq)
					return;
			}
		}
		// 140 CONTINUE;
		goto160q();
		if (exitq)
			return;
	}

	/**
	 * Used by DCADRE
	 */
	private void goto145q() {
		// "   NOISE MAY BE DOMINANT FEATURE       "
		// "   ESTIMATE NOISE LEVEL BY COMPARING "
		// "   THE VALUE OF THE INTEGRAND AT     "
		// "   FOUR *RANDOMLY CHOSEN* POINTS WITH"
		// "   THE VALUE OF THE STRAIGHT LINE    "
		// "   INTERPOLATING THE INTEGRAND AT THE"
		// "   TWO ENDPOINTS. IF SMALL ENOUGH,   "
		// "   ACCEPT *VINT*                     "
		// "   INTERGRATION OVER CURRENT SUB-      "
		// "   INTERVAL SUCCESSFUL               "
		// "   ADD *VINT* TO *DCADRE* AND *ERRER*"
		// "   TO *ERROR*, THEN SET UP NEXT SUB- "
		// "   INTERVAL, IF ANY.                 "

		SLOPEq = (FENDq - FBEGq) * TWOq;// 145 SLOPE=(FEND-FBEG)*TWO;
		FBEG2q = FBEGq + FBEGq;
		Iq = 0;// Iq=1;
		RVALq = BEGq + RNq[Iq] * STEPq;// 150
		DIFFq = DABS(func.F(RVALq) - FBEG2q - RNq[Iq] * SLOPEq);// 155 next
		ERRERq = DMAX1(ERRERq, ASTEPq * DIFFq);// 155
												// ERRER=DMAX1(ERRER,ASTEP*DIFF);
		if ((ERRERq > ERGOALq) && ((ERGLq + ERRERq) != ERGLq)) {
			goto175q();// GO TO 175;
			if (exitq)
				return;
		}
		Iq = Iq + 1;
		if (Iq <= 3) // IF(I<=4) GO TO 150;
		{
			goto150q();
			if (exitq)
				return;
		}
		IERq = 66;
		// "                              INTERGRATION OVER CURRENT SUB-      "
		// "                                INTERVAL SUCCESSFUL               "
		// "                                ADD *VINT* TO *DCADRE* AND *ERRER*"
		// "                                TO *ERROR*, THEN SET UP NEXT SUB- "
		// "                                INTERVAL, IF ANY.                 "
		CADREq = CADREq + VINTq;// 160 CADRE=CADRE + VINT;
		ERRORq = ERRORq + ERRERq;
		if (RIGHTq)// IF(RIGHT) GO TO 165;
		{
			goto165q();
			if (exitq)
				return;
		}
		ISTAGEq = ISTAGEq - 1;
		if (ISTAGEq == 0) // IF(ISTAGE.EQ.0) GO TO 220;
		{
			IERq = 131; // 200 IER=131;GO TO 215;
			goto215q();
			exitq = true;
			return;
		}
		REGLARq = REGLSVq[ISTAGEq - 1];// REGLAR=REGLSV[ISTAGE];
		BEGq = BEGINq[ISTAGEq - 1];// BEG=BEGIN(ISTAGE);
		ENDq = FINISq[ISTAGEq - 1];// END=FINIS(ISTAGE);
		CURESTq = CURESTq - ESTq[ISTAGEq] + VINTq;// CUREST=CUREST -
													// EST(ISTAGE+1) + VINT;
		IENDq = IBEGq - 1;
		FENDq = TSq[IENDq - 1];// FEND=TS(IEND);
		IBEGq = IBEGSq[ISTAGEq - 1];// IBEG=IBEGS(ISTAGE);
		goto180q();
		if (exitq)
			return;
	}

	/**
	 * Used by DCADRE
	 */
	private void goto150q() {
		RVALq = BEGq + RNq[Iq] * STEPq;
		DIFFq = DABS(func.F(RVALq) - FBEG2q - RNq[Iq] * SLOPEq);// 155 next
		ERRERq = DMAX1(ERRERq, ASTEPq * DIFFq);// 155
												// ERRER=DMAX1(ERRER,ASTEP*DIFF);
		if ((ERRERq > ERGOALq) && ((ERGLq + ERRERq) != ERGLq)) {
			goto175q();// GO TO 175;
			if (exitq)
				return;
		}
		Iq = Iq + 1;
		if (Iq <= 3) // IF(I<=4) GO TO 150;
		{
			goto150q();
			if (exitq)
				return;
		}
		IERq = 66;
		// "                              INTERGRATION OVER CURRENT SUB-      "
		// "                                INTERVAL SUCCESSFUL               "
		// "                                ADD *VINT* TO *DCADRE* AND *ERRER*"
		// "                                TO *ERROR*, THEN SET UP NEXT SUB- "
		// "                                INTERVAL, IF ANY.                 "
		CADREq = CADREq + VINTq;// 160 CADRE=CADRE + VINT;
		ERRORq = ERRORq + ERRERq;
		if (RIGHTq)// IF(RIGHT) GO TO 165;
		{
			goto165q();
			if (exitq)
				return;
		}
		ISTAGEq = ISTAGEq - 1;
		if (ISTAGEq == 0) // IF(ISTAGE.EQ.0) GO TO 220;
		{
			IERq = 131; // 200 IER=131;GO TO 215;
			goto215q();
			exitq = true;
			return;
		}
		REGLARq = REGLSVq[ISTAGEq - 1];// REGLAR=REGLSV[ISTAGE];
		BEGq = BEGINq[ISTAGEq - 1];// BEG=BEGIN(ISTAGE);
		ENDq = FINISq[ISTAGEq - 1];// END=FINIS(ISTAGE);
		CURESTq = CURESTq - ESTq[ISTAGEq] + VINTq;// CUREST=CUREST -
													// EST(ISTAGE+1) + VINT;
		IENDq = IBEGq - 1;
		FENDq = TSq[IENDq - 1];// FEND=TS(IEND);
		IBEGq = IBEGSq[ISTAGEq - 1];// IBEG=IBEGS(ISTAGE);
		goto180q();
		if (exitq)
			return;
	}

	/**
	 * Used by DCADRE
	 */
	private void goto155q() {
		ERRERq = DMAX1(ERRERq, ASTEPq * DIFFq);// 155
												// ERRER=DMAX1(ERRER,ASTEP*DIFF);
		if ((ERRERq > ERGOALq) && ((ERGLq + ERRERq) != ERGLq)) {
			goto175q();// GO TO 175;
			if (exitq)
				return;
		}
		Iq = Iq + 1;
		if (Iq <= 3) // IF(I<=4) GO TO 150;
		{
			goto150q();
			if (exitq)
				return;
		}
		IERq = 66;
		// "                              INTERGRATION OVER CURRENT SUB-      "
		// "                                INTERVAL SUCCESSFUL               "
		// "                                ADD *VINT* TO *DCADRE* AND *ERRER*"
		// "                                TO *ERROR*, THEN SET UP NEXT SUB- "
		// "                                INTERVAL, IF ANY.                 "
		CADREq = CADREq + VINTq;// 160 CADRE=CADRE + VINT;
		ERRORq = ERRORq + ERRERq;
		if (RIGHTq)// IF(RIGHT) GO TO 165;
		{
			goto165q();
			if (exitq)
				return;
		}
		ISTAGEq = ISTAGEq - 1;
		if (ISTAGEq == 0) // IF(ISTAGE.EQ.0) GO TO 220;
		{
			IERq = 131; // 200 IER=131;GO TO 215;
			goto215q();
			exitq = true;
			return;
		}
		REGLARq = REGLSVq[ISTAGEq - 1];// REGLAR=REGLSV[ISTAGE];
		BEGq = BEGINq[ISTAGEq - 1];// BEG=BEGIN(ISTAGE);
		ENDq = FINISq[ISTAGEq - 1];// END=FINIS(ISTAGE);
		CURESTq = CURESTq - ESTq[ISTAGEq] + VINTq;// CUREST=CUREST -
													// EST(ISTAGE+1) + VINT;
		IENDq = IBEGq - 1;
		FENDq = TSq[IENDq - 1];// FEND=TS(IEND);
		IBEGq = IBEGSq[ISTAGEq - 1];// IBEG=IBEGS(ISTAGE);
		goto180q();
		if (exitq)
			return;
	}

	/**
	 * Used by DCADRE
	 */
	private void goto160q() {
		// "                              INTERGRATION OVER CURRENT SUB-      "
		// "                                INTERVAL SUCCESSFUL               "
		// "                                ADD *VINT* TO *DCADRE* AND *ERRER*"
		// "                                TO *ERROR*, THEN SET UP NEXT SUB- "
		// "                                INTERVAL, IF ANY.                 "
		CADREq = CADREq + VINTq;// 160 CADRE=CADRE + VINT;
		ERRORq = ERRORq + ERRERq;
		if (RIGHTq) {
			goto165q();
			if (exitq)
				return;
		}
		ISTAGEq = ISTAGEq - 1;
		if (ISTAGEq == 0) // IF(ISTAGE.EQ.0) GO TO 220;
		{
			IERq = 131; // 200 IER=131;GO TO 215;
			goto215q();
			exitq = true;
			return;
		}
		REGLARq = REGLSVq[ISTAGEq - 1];// REGLAR=REGLSV[ISTAGE];
		BEGq = BEGINq[ISTAGEq - 1];// BEG=BEGIN(ISTAGE);
		ENDq = FINISq[ISTAGEq - 1];// END=FINIS(ISTAGE);
		CURESTq = CURESTq - ESTq[ISTAGEq] + VINTq;// CUREST=CUREST -
													// EST(ISTAGE+1) + VINT;
		IENDq = IBEGq - 1;
		FENDq = TSq[IENDq - 1];// FEND=TS(IEND);
		IBEGq = IBEGSq[ISTAGEq - 1];// IBEG=IBEGS(ISTAGE);
		goto180q();
		if (exitq)
			return;
	}

	/**
	 * Used by DCADRE
	 */
	private void goto165q() {
		CURESTq = CURESTq + VINTq;// 165 CUREST=CUREST + VINT;
		STAGEq = STAGEq + STAGEq;
		IENDq = IBEGq;
		IBEGq = IBEGSq[ISTAGEq - 1];// IBEG=IBEGS(ISTAGE);
		ENDq = BEGq;
		BEGq = BEGINq[ISTAGEq - 1];// BEG=BEGIN(ISTAGE);
		FENDq = FBEGq;
		FBEGq = TSq[IBEGq - 1];// FBEG=TS(IBEG);
		goto5q();// @@@@@
		if (exitq)
			return;
	}

	/**
	 * Used by DCADRE
	 */
	private void goto170q() {
		// "                              INTEGRATION OVER CURRENT SUBINTERVAL"
		// "                                IS UNSUCCESSFUL. MARK SUBINTERVAL "
		// "                                FOR FURTHER SUBDIVISION. SET UP   "
		// "                                NEXT SUBINTERVAL.                 "
		REGLARq = true;
		if (ISTAGEq == MXSTGEq)// 175 IF(ISTAGE.EQ.MXSTGE) GO TO 205;
		{
			IERq = 132;// 205 IER=132;
			goto215q();
			exitq = true;
			return;
		}
		if (RIGHTq) {
			goto185q();
			if (exitq)
				return;
		}
		REGLSVq[ISTAGEq] = REGLARq;// REGLSV(ISTAGE+1)=REGLAR;
		BEGINq[ISTAGEq - 1] = BEGq;// BEGIN(ISTAGE)=BEG;
		IBEGSq[ISTAGEq - 1] = IBEGq;// IBEGS(ISTAGE)=IBEG;
		STAGEq = STAGEq * HALFq;
		RIGHTq = true;// 180 RIGHT=.TRUE.;
		BEGq = (BEGq + ENDq) * HALFq;
		IBEGq = (IBEGq + IENDq) / 2;
		TSq[IBEGq - 1] = TSq[IBEGq - 1] * HALFq;// TS(IBEG)=TS(IBEG)*HALF;
		FBEGq = TSq[IBEGq - 1];// FBEG=TS(IBEG);
		goto10q();// @@;
		if (exitq)
			return;
	}

	/**
	 * Used by DCADRE
	 */
	private void goto175q() {
		if (ISTAGEq == MXSTGEq)// 175 IF(ISTAGE.EQ.MXSTGE) GO TO 205;
		{
			IERq = 132;// 205 IER=132;
			goto215q();
			exitq = true;
			return;
		}
		if (RIGHTq) {
			goto185q();
			if (exitq)
				return;
		}
		REGLSVq[ISTAGEq] = REGLARq;// REGLSV(ISTAGE+1)=REGLAR;
		BEGINq[ISTAGEq - 1] = BEGq;// BEGIN(ISTAGE)=BEG;
		IBEGSq[ISTAGEq - 1] = IBEGq;// IBEGS(ISTAGE)=IBEG;
		STAGEq = STAGEq * HALFq;
		RIGHTq = true;// 180 RIGHT=.TRUE.;
		BEGq = (BEGq + ENDq) * HALFq;
		IBEGq = (IBEGq + IENDq) / 2;
		TSq[IBEGq - 1] = TSq[IBEGq - 1] * HALFq;// TS(IBEG)=TS(IBEG)*HALF;
		FBEGq = TSq[IBEGq - 1];// FBEG=TS(IBEG);
		goto10q();// @@;
		if (exitq)
			return;
	}

	/**
	 * Used by DCADRE
	 */
	private void goto180q() {
		RIGHTq = true;// 180 RIGHT=.TRUE.;
		BEGq = (BEGq + ENDq) * HALFq;
		IBEGq = (IBEGq + IENDq) / 2;
		TSq[IBEGq - 1] = TSq[IBEGq - 1] * HALFq;// TS(IBEG)=TS(IBEG)*HALF;
		FBEGq = TSq[IBEGq - 1];// FBEG=TS(IBEG);
		goto10q();// @@@@
		if (exitq)
			return;
	}

	/**
	 * Used by DCADRE
	 */
	private void goto185q() {
		NNLEFTq = IBEGq - IBEGSq[ISTAGEq - 1];// 185 NNLEFT=IBEG -
												// IBEGS(ISTAGE);
		if ((IENDq + NNLEFTq) >= MAXTSq)// IF(IEND+NNLEFT.GE.MAXTS) GO TO 200;
		{
			IERq = 131;// 200 IER=131;GO TO 215;
			goto215q();
			exitq = true;
			return;// 9005 RETURN;
		}
		IIIq = IBEGSq[ISTAGEq - 1];// III=IBEGS(ISTAGE);
		IIq = IENDq;
		for (int I1 = IIIq - 1; I1 < IBEGq; I1++)// DO 190
													// I=III,IBEG;^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
		{
			IIq = IIq + 1;
			TSq[IIq - 1] = TSq[I1];// TS(II)=TS(I);
		}// 190 CONTINUE;
		for (int I1 = IBEGq - 1; I1 < IIq; I1++)// DO 195
												// I=IBEG,II;II^^^^^^^^^^^^^OK
		{
			TSq[IIIq - 1] = TSq[I1];// TS(III)=TS(I);
			IIIq = IIIq + 1;
		}// 195 CONTINUE;
		IENDq = IENDq + 1;
		IBEGq = IENDq - NNLEFTq;
		FENDq = FBEGq;
		FBEGq = TSq[IBEGq - 1];// FBEG=TS(IBEG);
		FINISq[ISTAGEq - 1] = ENDq;// FINIS(ISTAGE)=END;
		ENDq = BEGq;
		BEGq = BEGINq[ISTAGEq - 1];// BEG=BEGIN(ISTAGE);
		BEGINq[ISTAGEq - 1] = ENDq;// BEGIN(ISTAGE)=END;
		REGLSVq[ISTAGEq - 1] = REGLARq;// REGLSV(ISTAGE)=REGLAR;
		ISTAGEq = ISTAGEq + 1;
		REGLARq = REGLSVq[ISTAGEq - 1];// REGLAR=REGLSV(ISTAGE);
		ESTq[ISTAGEq - 1] = VINTq;// EST(ISTAGE)=VINT;
		CURESTq = CURESTq + ESTq[ISTAGEq - 1];// CUREST=CUREST + EST(ISTAGE);
		goto5q();// @@;
		if (exitq)
			return;
	}

	/**
	 * Used by DCADRE
	 */
	private void goto215q() {
		CADREq = CURESTq + VINTq;// 215 CADRE=CUREST + VINT;
		DCADREq = CADREq;// //220 DCADRE=CADRE;+//9000 CONTINUE;
	}

	/**
	 * Used by DCADRE
	 */
	private void resetDCADRE() {
		DIFFq = 0.;
		H2CONVq = false;
		AITKENq = false;
		LM1q = 0;
		N2q = 0;
		FNq = 0.;
		ISTEPq = 0;
		IIq = 0;
		IIIq = 0;
		HOVNq = 0.;
		FIq = 0.;
		ISTEP2q = 0;
		SUMq = 0.;
		SUMABSq = 0.;
		ABSIq = 0.;
		ITq = 0;
		TABTLMq = 0.;
		ERGLq = 0.;
		ERGOALq = 0.;
		FEXTRPq = 0.;
		ERRERq = 0.;
		RIGHTq = false;
		STEPq = 0.;
		ASTEPq = 0.;
		TABSq = 0.;
		Lq = 0;
		Nq = 0;
		TSq = new double[2049];
		Tq = new double[10][10];
		ZEROq = 0.0;
		P1q = 0.1;
		HALFq = 0.5;
		ONEq = 1.0;
		TWOq = 2.0;
		FOURq = 4.0;
		FOURP5q = 4.5;
		TENq = 10.0;
		HUNq = 100.0;
		AITLOWq = 1.1;
		H2TOLq = 0.15;
		AITTOLq = 0.1;
		JUMPTLq = 0.01;
		MAXTSq = 2049;
		MAXTBLq = 10;
		MXSTGEq = 30;
		SLOPEq = 0.;
		FBEG2q = 0.;
		RNq = new double[4];
		RNq[0] = 0.7142005;
		RNq[1] = 0.3466282;
		RNq[2] = 0.843751;
		RNq[3] = 0.1263305;
		REGLSVq = new boolean[30];
		BEGINq = new double[30];
		FINISq = new double[30];
		ESTq = new double[30];
		IBEGSq = new int[30];
		NNLEFTq = 0;
		Iq = 0;
		// ---------------------------------
		DCADREq = 0.;
		ALG4O2q = DLOG10(TWOq);
		CADREq = ZEROq;
		CURESTq = ZEROq;
		VINTq = ZEROq;
		LENGTHq = 0.;
		ERRRq = 0.;
		ERRAq = 0.;
		STEPMNq = 0.;
		STEPNMq = 0.;
		BEGq = 0.;
		RVALq = 0.;
		FBEGq = 0.;
		STAGEq = HALFq;
		ISTAGEq = 1;
		FNSIZEq = ZEROq;
		PREVERq = ZEROq;
		REGLARq = false;
		FBEGq = 0.;
		IBEGq = 1;
		ENDq = 0.;
		FENDq = 0.;
		IENDq = 2;
		ERRORq = ZEROq;
		IERq = 0;
		LENGTHq = DABS(Bq - Aq);
		ERRRq = RERRq;
		ERRAq = DABS(AERRq);
		STEPMNq = (LENGTHq / Math.pow(2.0, MXSTGEq));
		STEPNMq = DMAX1(LENGTHq, DABS(Aq), DABS(Bq)) * TENq;
		// " THE GIVEN INTERVAL OF INTEGRATION   "
		// " IS THE FIRST INTERVAL CONSIDERED. "
		BEGq = Aq;
		RVALq = BEGq;
		FBEGq = func.F(RVALq) * HALFq;// @@
		TSq[0] = FBEGq;// (1)
		ENDq = Bq;
		RVALq = ENDq;
		FENDq = func.F(RVALq) * HALFq;// @@
		TSq[1] = FENDq;// (2)
		SINGq = 0.;
		FEXTM1q = 0.;
		Rq = new double[10];
		AITq = new double[10];
		DIFq = new double[10];
		ALPHAq = 0.;
		H2NXTq = 0.;
		SINGNXq = 0.;
		ERRETq = 0.;
		H2TFEXq = 0.;
		// -----------------------
		exitq = false;
	}

	// max from 3 var

	// compute log in 10 base
	// a**(LOGaX)=X due to definition of logaritm.
	// the logarithm of a power equals the exponent multiplied with the
	// logarithm of the base
	// LOGbX=LOGb[a**(LOGaX)]=LOGaX*LOGba!!!
	// so LOG10X=LnX/Ln10
	/**
	 * Logarithm in base 10 of a number x
	 * @param x x
	 * @return the result
	 */
	public static double DLOG10(double x) {
		double result = Math.log(x);
		result = result / Math.log(10);
		return result;
	}

	/**
	 * Maximum of 3 numbers
	 * @param a a
	 * @param b b
	 * @param c c
	 * @return the result
	 */
	public static double DMAX1(double a, double b, double c) {
		double r = Math.max(a, b);
		r = Math.max(r, c);
		return r;
	}

	/**
	 * Maximum of 2 numbers
	 * @param a a
	 * @param b b
	 * @return the result
	 */
	public static double DMAX1(double a, double b) {
		double r = Math.max(a, b);
		return r;
	}

	/**
	 * Absolute value of a number
	 * @param x x
	 * @return the result
	 */
	public static double DABS(double x) {
		double result = Math.abs(x);
		return result;
	}
	// ==================================================================================
}
