package danfulea.math.numerical;

/**
 * Evaluation of special functions class.
 * Based on literature, including Numerical Recipes in C (Cambridge Univ.).
 * @author Dan Fulea, 06 OCT. 2006
 */
public class SpecFunc {
	public static boolean failB = false;
	public static String failS = "";

	public static int ntop = 4;
	public static double[] afctln = new double[101];
	public static int ITMAX = 100;
	public static double EPS = 3.0e-7;
	public static double FPMIN = 1.0e-30;
	public static double gamser_out = 0.0;
	public static double gln_out_gser = 0.0;
	public static double gammcf_out = 0.0;
	public static double gln_out_gcf = 0.0;

	public static int MAXIT = 100;// Maximum allowed number of iterations.
	public static double EULER = 0.5772156649;// Euler’s constant gamma lower
												// case.
	// #define FPMIN 1.0e-30 Close to smallest representable floating-point
	// number.
	public static double EPS1 = 1.0e-7;// Desired relative error, not smaller
										// than the machine precision.
	public static double EPS2 = 6.0e-8;

	public static double ACC = 40.0;// Make larger to increase accuracy.
	public static double BIGNO = 1.0e10;
	public static double BIGNI = 1.0e-10;

	public static double EPSb = 1.0e-10;
	public static double FPMINb = 1.0e-30;
	public static int MAXITb = 10000;
	public static double XMIN = 2.0;
	public static double PI = 3.141592653589793;
	public static double rjb = 0.0;
	public static double ryb = 0.0;
	public static double rjpb = 0.0;
	public static double rypb = 0.0;
	public static int NUSE1 = 5;
	public static int NUSE2 = 5;
	public static double gam1b = 0.0;
	public static double gam2b = 0.0;
	public static double gamplb = 0.0;
	public static double gammib = 0.0;
	public static double rib = 0.0;
	public static double rkb = 0.0;
	public static double ripb = 0.0;
	public static double rkpb = 0.0;

	public static double THIRD = 1.0 / 3.0;
	public static double TWOTHR = 2.0 * THIRD;
	public static double ONOVRT = 0.57735027;
	public static double ai_airy = 0.0;
	public static double bi_airy = 0.0;
	public static double aip_airy = 0.0;
	public static double bip_airy = 0.0;

	public static double RTPIO2 = 1.2533141;
	public static double sjb = 0.0;
	public static double syb = 0.0;
	public static double sjpb = 0.0;
	public static double sypb = 0.0;

	public static double XMIN2 = 1.5;
	public static double PIBY2 = (PI / 2.0);
	public static int TRUE = 1;
	public static Complex ONE = new Complex(1.0, 0.0);
	public static double s_fresnel = 0.0;
	public static double c_fresnel = 0.0;

	public static double TMIN = 2.0;
	public static double ci_integral = 0.0;
	public static double si_integral = 0.0;

	public static int NMAX = 6;
	public static double H = 0.4;
	public static double A1 = (2.0 / 3.0);
	public static double A2 = 0.4;
	public static double A3 = (2.0 / 7.0);
	public static double[] c = new double[NMAX + 1];
	public static int init = 0;

	public static double ERRTOL = 0.08;
	public static double TINY = 1.5e-38;
	public static double BIG = 3.0e37;
	public static double C1 = (1.0 / 24.0);
	public static double C2 = 0.1;
	public static double C3 = (3.0 / 44.0);
	public static double C4 = (1.0 / 14.0);

	public static double ERRTOL1 = 0.05;
	public static double TINY1 = 1.0e-25;
	public static double BIG1 = 4.5e21;
	public static double C11 = (3.0 / 14.0);
	public static double C21 = (1.0 / 6.0);
	public static double C31 = (9.0 / 22.0);
	public static double C41 = (3.0 / 26.0);
	public static double C51 = (0.25 * C31);
	public static double C61 = (1.5 * C41);

	public static double ERRTOL2 = 0.05;
	public static double TINY2 = 2.5e-13;
	public static double BIG2 = 9.0e11;
	public static double C12 = (3.0 / 14.0);
	public static double C22 = (1.0 / 3.0);
	public static double C32 = (3.0 / 22.0);
	public static double C42 = (3.0 / 26.0);
	public static double C52 = (0.75 * C32);
	public static double C62 = (1.5 * C42);
	public static double C72 = (0.5 * C22);
	public static double C82 = (C32 + C32);

	public static double ERRTOL3 = 0.04;
	public static double TINY3 = 1.69e-38;
	public static double SQRTNY = 1.3e-19;
	public static double BIG3 = 3.e37;
	public static double TNBG = (TINY3 * BIG3);
	public static double COMP1 = (2.236 / SQRTNY);
	public static double COMP2 = (TNBG * TNBG / 25.0);
	public static double C13 = 0.3;
	public static double C23 = (1.0 / 7.0);
	public static double C33 = 0.375;
	public static double C43 = (9.0 / 22.0);

	public static double CA = 0.0003;
	public static double sn_jacob = 0.0;
	public static double cn_jacob = 0.0;
	public static double dn_jacob = 0.0;

	public static double EPS3 = 1.0e-6;
	public static Complex aa, bb, cc, z0, dz;
	public static int kmax, kount;

	// public static double *xp,**yp,dxsav;

	/**
	 * Returns the value ln[GAMMA(xx)] for xx greater than 0. 
	 * @param xx xx
	 * @return the result
	 */
	public static double gammln(double xx)
	// Returns the value ln[GAMMA(xx)] for xx > 0.
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

	/*
	 * How shall we write a routine for the factorial function n!? Generally the
	 * factorial function will be called for small integer values (for large
	 * values it will overflow anyway!), and in most applications the same
	 * integer value will be called for many times. It is a profligate waste of
	 * computer time to call exp(gammln(n+1.0)) for each required factorial.
	 * Better to go back to basics, holding gammln in reserve for unlikely
	 * calls:
	 */
	/**
	 * Returns the value n! as a double number.
	 * @param n n
	 * @return the result
	 */
	public static double factrl(int n)
	// Returns the value n! as a floating-point number.
	{
		// float gammln(float xx);
		// void nrerror(char error_text[]);
		// static int ntop=4;

		failB = false;
		// static float a[33]={1.0,1.0,2.0,6.0,24.0}; //Fill in table only as
		// required.
		double[] a = new double[33];// {1.0,1.0,2.0,6.0,24.0};
		a[0] = 1.0;
		a[1] = 1.0;
		a[2] = 2.0;
		a[3] = 6.0;
		a[4] = 24.0;
		int j = 0;
		if (n < 0) {
			// nrerror("Negative factorial in routine factrl");
			failB = true;
			failS = "Negative factorial in routine factrl";
			return 1.0;
		}
		if (n > 32)
			return Math.exp(gammln(n + 1.0));
		// Larger value than size of table is required. Actually, this big a
		// value is going to overflow
		// on many computers, but no harm in trying.
		while (ntop < n) {// Fill in table up to desired value.
			j = ntop++;// j is still ntop; but ntop becoms ntop+1
			a[ntop] = a[j] * ntop;
		}
		return a[n];
	}

	/**
	 * Returns the binomial coefficient (n,k)=Cn^k as a double number.
	 * @param n n
	 * @param k k
	 * @return the result
	 */
	public static double bico(int n, int k)
	// Returns the binomial coefficient (n,k)=Cnk
	// as a floating-point number.
	{
		// float factln(int n);
		return Math
				.floor(0.5 + Math.exp(factln(n) - factln(k) - factln(n - k)));
		// The floor function cleans up roundoff error for smaller values of n
		// and k.
	}

	/**
	 * Returns ln(n!).
	 * @param n n
	 * @return the result
	 */
	public static double factln(int n)
	// Returns ln(n!).
	{
		// float gammln(float xx);
		// void nrerror(char error_text[]);
		// static float a[101]; A static array is automatically initialized to
		// zero.
		failB = false;
		// double[] a = new double[101];
		if (n < 0) {
			// nrerror("Negative factorial in routine factln");
			failB = true;
			failS = "Negative factorial in routine factln";
			return 1.0;
		}
		if (n <= 1)
			return 0.0;
		if (n <= 100)// return a[n] ? a[n] : (a[n]=gammln(n+1.0)); //In range of
						// table.
			return afctln[n] != 0.0 ? afctln[n] : (afctln[n] = gammln(n + 1.0));
		else
			return gammln(n + 1.0); // Out of range of table.
	}

	/*
	 * B(z,w) =gamma(z)gamma(w)/gamma(z + w)
	 */
	/**
	 * Returns the value of the beta function B(z,w).
	 * @param z z
	 * @param w w
	 * @return the result
	 */
	public static double beta(double z, double w)
	// Returns the value of the beta function B(z,w).
	{
		// float gammln(float xx);
		return Math.exp(gammln(z) + gammln(w) - gammln(z + w));
	}

	// Incomplete Gamma Function,Error Function, Chi-Square Probability
	// Function,
	// Cumulative Poisson Function
	/*
	 * P(a, x) =gammalowercase(a, x)/gamma(a) =[1/gamma(a)] integral 0-x of
	 * exp(-t)t^(a-1)dt;a>0 Q(a, x) =1-P(a,x)=gamma(a, x)/gamma(a) =[1/gamma(a)]
	 * integral x-+inf of exp(-t)t^(a-1)dt;a>0
	 */
	/**
	 * Returns the incomplete gamma function P(a, x) = [1/gamma(a)] Integral (from 0 to x) of exp(-t)t^(a-1)dt; a greater than 0.
	 * @param a a
	 * @param x x
	 * @return the result
	 */
	public static double gammp(double a, double x)
	// Returns the incomplete gamma function P(a, x).
	{
		// void gcf(float *gammcf, float a, float x, float *gln);
		// void gser(float *gamser, float a, float x, float *gln);
		// void nrerror(char error_text[]);
		double gamser = 0.0;
		// double gammcf = 0.0;
		// double gln = 0.0;
		failB = false;

		if (x < 0.0 || a <= 0.0) {
			// nrerror("Invalid arguments in routine gammp");
			failB = true;
			failS = "Invalid arguments in routine gammp";
			return 1.0;
		}
		if (x < (a + 1.0)) {// Use the series representation.
							// gser(&gamser,a,x,&gln);
			gser(a, x);
			gamser = gamser_out;
			return gamser;
		} else {// Use the continued fraction representation
				// gcf(&gammcf,a,x,&gln);
			gcf(a, x);
			// return 1.0-gammcf; //and take its complement.
			return 1.0 - gammcf_out;
		}
	}

	/**
	 * Returns the incomplete gamma function Q(a, x) = 1 - P(a, x).
	 * @param a a
	 * @param x x
	 * @return the result
	 */
	public static double gammq(double a, double x)
	// Returns the incomplete gamma function Q(a, x) = 1 - P(a, x).
	{
		// void gcf(float *gammcf, float a, float x, float *gln);
		// void gser(float *gamser, float a, float x, float *gln);
		// void nrerror(char error_text[]);
		// double gamser = 0.0;
		double gammcf = 0.0;
		// double gln = 0.0;

		failB = false;
		if (x < 0.0 || a <= 0.0) {
			// nrerror("Invalid arguments in routine gammq");
			failB = true;
			failS = "Invalid arguments in routine gammq";
			return 1.0;
		}
		if (x < (a + 1.0)) {// Use the series representation
			gser(a, x);// gser(&gamser,a,x,&gln);
			// return 1.0-gamser; and take its complement.
			return 1.0 - gamser_out;
		} else {// Use the continued fraction representation.
			gcf(a, x);// gcf(&gammcf,a,x,&gln);
			gammcf = gammcf_out;
			return gammcf;
		}
	}

	/**
	 * Returns the incomplete gamma function P(a, x) evaluated by its series representation as gamser_out. 
	 * Also returns ln GAMMA(a) as gln_out_gser.
	 * @param a a
	 * @param x x
	 */
	public static void gser(double a, double x)// (float *gamser, float a, float
												// x, float *gln)
	// Returns the incomplete gamma function P(a, x) evaluated by its series
	// representation as gamser.
	// Also returns ln GAMMA(a) as gln.
	{
		// float gammln(float xx);
		// void nrerror(char error_text[]);
		int n = 0;
		double sum = 0.0;
		double del = 0.0;
		double ap = 0.0;
		failB = false;

		gln_out_gser = gammln(a);// *gln=gammln(a);
		if (x <= 0.0) {
			if (x < 0.0) {
				// nrerror("x less than 0 in routine gser");
				failB = true;
				failS = "x less than 0 in routine gser";
				return;
			}
			gamser_out = 0.0;// *gamser=0.0;
			return;
		} else {
			ap = a;
			del = sum = 1.0 / a;
			for (n = 1; n <= ITMAX; n++) {
				++ap;
				del *= x / ap;
				sum += del;
				if (Math.abs(del) < Math.abs(sum) * EPS)// (fabs(del) <
														// fabs(sum)*EPS)
				{
					// *gamser=sum*exp(-x+a*log(x)-(*gln));
					gamser_out = sum
							* Math.exp(-x + a * Math.log(x) - gln_out_gser);
					return;
				}
			}
			// nrerror("a too large, ITMAX too small in routine gser");
			failB = true;
			failS = "a too large, ITMAX too small in routine gser";

			return;
		}
	}

	/**
	 * Returns the incomplete gamma function Q(a, x) evaluated by its continued fraction representation 
	 * as gammcf_out. Also returns lnGAMMA(a) as gln_out_gcf.
	 * @param a a
	 * @param x x
	 */
	public static void gcf(double a, double x)// (float *gammcf, float a, float
												// x, float *gln)
	// Returns the incomplete gamma function Q(a, x) evaluated by its continued
	// fraction representation
	// as gammcf. Also returns ln?(a) as gln.
	{
		// float gammln(float xx);
		// void nrerror(char error_text[]);
		int i = 0;
		double an = 0.0;
		double b = 0.0;
		double c = 0.0;
		double d = 0.0;
		double del = 0.0;
		double h = 0.0;

		failB = false;

		gln_out_gcf = gammln(a);// *gln=gammln(a);
		b = x + 1.0 - a; // Set up for evaluating continued fraction
							// by modified Lentz’s method (§5.2) with b0 = 0.
		c = 1.0 / FPMIN;
		d = 1.0 / b;
		h = d;
		for (i = 1; i <= ITMAX; i++) {// Iterate to convergence.
			an = -i * (i - a);
			b += 2.0;
			d = an * d + b;
			if (Math.abs(d) < FPMIN)
				d = FPMIN;// if (fabs(d) < FPMIN) d=FPMIN;
			c = b + an / c;
			if (Math.abs(c) < FPMIN)
				c = FPMIN;
			d = 1.0 / d;
			del = d * c;
			h *= del;
			if (Math.abs(del - 1.0) < EPS)
				break;
		}
		if (i > ITMAX) {
			// nrerror("a too large, ITMAX too small in gcf");
			failB = true;
			failS = "a too large, ITMAX too small in routine gcf";

			return;
		}
		// *gammcf=exp(-x+a*log(x)-(*gln))*h; Put factors in front.
		gammcf_out = Math.exp(-x + a * Math.log(x) - gln_out_gcf) * h;
	}

	/*
	 * erf(x) =[2/sqrt(PI)]*Integral 0-x of exp(-t^2dt erfc(x)=1-erf(x)=
	 * [2/sqrt(PI)]*Integral x-+inf of exp(-t^2dt
	 */
	/**
	 * Returns the error function erf(x).
	 * @param x x
	 * @return the result
	 */
	public static double erff(double x)
	// Returns the error function erf(x).
	{
		// float gammp(float a, float x);
		return x < 0.0 ? -gammp(0.5, x * x) : gammp(0.5, x * x);
	}

	/**
	 * Returns the complementary error function erfc(x).
	 * @param x x
	 * @return the result
	 */
	public static double erffc(double x)
	// Returns the complementary error function erfc(x).
	{
		// float gammp(float a, float x);
		// float gammq(float a, float x);
		return x < 0.0 ? 1.0 + gammp(0.5, x * x) : gammq(0.5, x * x);
	}

	/**
	 * Returns the complementary error function erfc(x) with fractional error everywhere less than 1.2 × 10-7.
	 * @param x x
	 * @return the result
	 */
	public static double erfcc(double x)
	// Returns the complementary error function erfc(x) with fractional error
	// everywhere less than
	// 1.2 × 10-7.
	{
		double t = 0.0;
		double z = 0.0;
		double ans = 0.0;
		z = Math.abs(x);
		t = 1.0 / (1.0 + 0.5 * z);
		ans = t
				* Math.exp(-z
						* z
						- 1.26551223
						+ t
						* (1.00002368 + t
								* (0.37409196 + t
										* (0.09678418 + t
												* (-0.18628806 + t
														* (0.27886807 + t
																* (-1.13520398 + t
																		* (1.48851587 + t
																				* (-0.82215223 + t * 0.17087277)))))))));
		return x >= 0.0 ? ans : 2.0 - ans;
	}

	/**
	 * Another algorithm for error function based on EGSnrc
	 * @param X X
	 * @return the result
	 */
	public static double ERF1(double X) {
		double erf1 = 0.0;
		double x = X;

		// double precision A(0:22,2); //" Coefficients in expansion for erf(x)
		// if x<3
		// " (K=1) and for erfc(x) x>3 (K=2)
		// double precision
		// CONST, //" 2/sqrt(pi)
		// BN,BN1,BN2, //" Recursion coefficients B(n),B(n+1),B(n+2)
		// Y,FAC; //" y=x/3 or 3/x and FAC = 2(2y**2-1)
		// $INTEGER N, //" recursion index n
		// K, //" K=1,2 for x <= 3 or x > 3
		// NLIM(2); //" Maximum value of n in sum for K=1,2
		double y = 0.0;
		double FAC = 0.0;
		double BN = 0.0;
		double BN1 = 0.0;
		double BN2 = 0.0;
		int k = 0;
		int n = 0;
		double[][] A = { { 1.0954712997776232, 0.9750834237085559 },
				{ -0.2891754011269890, -0.0240493938504146 },
				{ 0.1104563986337951, 0.0008204522408804 },
				{ -0.0412531882278565, -0.0000434293081303 },
				{ 0.0140828380706516, 0.0000030184470340 },
				{ -0.0043292954474314, -0.0000002544733193 },
				{ 0.0011982719015923, 0.0000000248583530 },
				{ -0.0002999729623532, -0.0000000027317201 },
				{ 0.0000683258603789, 0.0000000003308472 },
				{ -0.0000142469884549, 0.0000000000001464 },
				{ 0.0000027354087728, -0.0000000000000244 },
				{ -0.0000004861912872, 0.0000000000000042 },
				{ 0.0000000803872762, -0.0000000000000008 },
				{ -0.0000000124184183, 0.0000000000000001 },
				{ 0.0000000017995326, 0.0 }, { -0.0000000002454795, 0.0 },
				{ 0.0000000000316251, 0.0 }, { -0.0000000000038590, 0.0 },
				{ 0.0000000000004472, 0.0 }, { -0.0000000000000493, 0.0 },
				{ 0.0000000000000052, 0.0 }, { -0.0000000000000005, 0.0 },
				{ 0.0000000000000001, 0.0 } };
		int[] NLIM = { 22, 16 };
		double CONST = 2 / Math.sqrt(Math.PI);// 1.128379167095513 /;

		if (x > 3) {
			y = 3.0 / x;
			k = 2;
		} else {
			y = x / 3.0;
			k = 1;
		}
		double Y = y;
		int K = k;
		// int N = n;
		// " Calculate sum of Chebyshev polynomials by backwards recursion
		// "
		// " sum { A(n)*T(2n+1;y) : n=0,N } = y * ( B(0) - B(1) )
		// " sum { A(n)*T(2n;y) : n=0,N } = ( B(0) - (2*y**2-1) * B(1) ) / 2
		// " = ( B(0) - B(2) + A(0) ) / 2
		// "
		// " where B(N+2) = B(N+1) = 0
		// " and B(n) = 2*(2*y**2-1)*B(n+1) - B(n+2) + A(n) for
		// n=N,(N-1),...,1,0
		// "
		FAC = 2.0 * (2.0 * Y * Y - 1.0);
		BN1 = 0.0; // " Initialise B(N+2) = 0
		BN = 0.0; // " Initialise B(N+1) = 0

		// for(n = NLIM(K),0,-1 [
		for (n = NLIM[K - 1]; n >= 0; n--) {
			BN2 = BN1;
			BN1 = BN;
			// BN = FAC * BN1 - BN2 + A(N,K)
			BN = FAC * BN1 - BN2 + A[n][K - 1];
		}

		if (k == 1) {
			erf1 = CONST * Y * (BN - BN1);
		} else {
			// erf1 = 1 - CONST * Math.exp(-X*X) * ( BN - BN2 + A(0,K) )/(4.0 *
			// X);
			erf1 = 1.0 - CONST * Math.exp(-X * X) * (BN - BN2 + A[0][K - 1])
					/ (4.0 * X);
		}

		return erf1;
	}

	// ****************************************************************************
	// erf from c
	/**
	 * Another algorithm for complementary error function based on EGSnrc
	 * @param x x
	 * @return th result
	 */
	public static double erfc1(double x)
	// ****************************************************************************
	//
	// Purpose:
	//
	// ERF1 evaluates the error function.
	//
	// Parameters:
	//
	// Input, double *X, the argument.
	//
	// Output, double ERF1, the value of the error function at X.
	//
	{
		double c = .564189583547756e0;
		double[] a = { .771058495001320e-04, -.133733772997339e-02,
				.323076579225834e-01, .479137145607681e-01,
				.128379167095513e+00 };
		double[] b = { .301048631703895e-02, .538971687740286e-01,
				.375795757275549e+00 };
		double[] p = { -1.36864857382717e-07, 5.64195517478974e-01,
				7.21175825088309e+00, 4.31622272220567e+01,
				1.52989285046940e+02, 3.39320816734344e+02,
				4.51918953711873e+02, 3.00459261020162e+02 };
		double[] q = { 1.00000000000000e+00, 1.27827273196294e+01,
				7.70001529352295e+01, 2.77585444743988e+02,
				6.38980264465631e+02, 9.31354094850610e+02,
				7.90950925327898e+02, 3.00459260956983e+02 };
		double[] r = { 2.10144126479064e+00, 2.62370141675169e+01,
				2.13688200555087e+01, 4.65807828718470e+00,
				2.82094791773523e-01 };
		double[] s = { 9.41537750555460e+01, 1.87114811799590e+02,
				9.90191814623914e+01, 1.80124575948747e+01 };
		double erf1, ax, bot, t, top, x2;

		ax = Math.abs(x);
		// if(ax > 0.5e0)
		// {
		// goto S10;
		// }
		if (ax <= 0.5e0) {
			t = x * x;
			top = (((a[0] * t + a[1]) * t + a[2]) * t + a[3]) * t + a[4]
					+ 1.0e0;
			bot = ((b[0] * t + b[1]) * t + b[2]) * t + 1.0e0;
			erf1 = x * (top / bot);
			return erf1;
		}
		// S10:
		// if(ax > 4.0e0) goto S20;
		else if (ax <= 4.0e0) {
			top = ((((((p[0] * ax + p[1]) * ax + p[2]) * ax + p[3]) * ax + p[4])
					* ax + p[5])
					* ax + p[6])
					* ax + p[7];
			bot = ((((((q[0] * ax + q[1]) * ax + q[2]) * ax + q[3]) * ax + q[4])
					* ax + q[5])
					* ax + q[6])
					* ax + q[7];
			erf1 = 0.5e0 + (0.5e0 - Math.exp(-(x * x)) * top / bot);
			if (x < 0.0e0)
				erf1 = -erf1;
			return erf1;
		}
		// S20:
		// if(ax >= 5.8e0) goto S30;
		else if (ax < 5.8e0) {
			x2 = x * x;
			t = 1.0e0 / x2;
			top = (((r[0] * t + r[1]) * t + r[2]) * t + r[3]) * t + r[4];
			bot = (((s[0] * t + s[1]) * t + s[2]) * t + s[3]) * t + 1.0e0;
			erf1 = (c - top / (x2 * bot)) / ax;
			erf1 = 0.5e0 + (0.5e0 - Math.exp(-x2) * erf1);
			if (x < 0.0e0)
				erf1 = -erf1;
			return erf1;
		}
		// S30:
		// erf1 = fifdsign(1.0e0,*x);
		else {
			erf1 = 1.0e0;
			if (x < 0.0e0)
				erf1 = -erf1;
		}
		return erf1;
	}

	/*
	 * Exponential integral En(x) = integral from 1 to +inf of exp(-xt)/t^n dt,
	 * x > 0, n= 0, 1, . . .
	 */
	/**
	 * Evaluates the exponential integral En(x)= integral from 1 to +inf of exp(-xt)/t^n dt.
	 * @param n n
	 * @param x x
	 * @return the result
	 */
	public static double expint(int n, double x)
	// Evaluates the exponential integral En(x).
	{
		// void nrerror(char error_text[]);
		int i = 0;
		int ii = 0;
		int nm1 = 0;
		double a = 0.0;
		double b = 0.0;
		double c = 0.0;
		double d = 0.0;
		double del = 0.0;
		double fact = 0.0;
		double h = 0.0;
		double psi = 0.0;
		double ans = 0.0;

		failB = false;

		nm1 = n - 1;
		if (n < 0 || x < 0.0 || (x == 0.0 && (n == 0 || n == 1))) {
			// nrerror("bad arguments in expint");
			failB = true;
			failS = "bad arguments in expint";

			return 0.0;

		} else// e1
		{
			if (n == 0)
				ans = Math.exp(-x) / x; // Special case.
			else// e2
			{
				if (x == 0.0)
					ans = 1.0 / nm1; // Another special case.
				else// e3
				{
					if (x > 1.0) {// Lentz’s algorithm (§5.2).
						b = x + n;
						c = 1.0 / FPMIN;
						d = 1.0 / b;
						h = d;
						for (i = 1; i <= MAXIT; i++) {
							a = -i * (nm1 + i);
							b += 2.0;
							d = 1.0 / (a * d + b);// Denominators cannot be
													// zero.
							c = b + a / c;
							del = c * d;
							h *= del;
							if (Math.abs(del - 1.0) < EPS1) {
								ans = h * Math.exp(-x);
								return ans;
							}
						}
						// nrerror("continued fraction failed in expint");
						failB = true;
						failS = "continued fraction failed in expint";

						return 0.0;
					} else {// Evaluate series.
						ans = (nm1 != 0 ? 1.0 / nm1 : -Math.log(x) - EULER); // Set
																				// first
																				// term.
						fact = 1.0;
						for (i = 1; i <= MAXIT; i++) {
							fact *= -x / i;
							if (i != nm1)
								del = -fact / (i - nm1);
							else {
								psi = -EULER; // Compute ?(n).
								for (ii = 1; ii <= nm1; ii++)
									psi += 1.0 / ii;
								del = fact * (-Math.log(x) + psi);
							}
							ans += del;
							if (Math.abs(del) < Math.abs(ans) * EPS1)
								return ans;
						}
						// nrerror("series failed in expint");
						failB = true;
						failS = "series failed in expint";

						return 0.0;

					}
				}// //e3
			}// e2
		}// e1
		return ans;
	}

	/*
	 * Ei(x) = -integral from -x to +inf of exp(-t)/t dt = integral form -inf to
	 * x of exp(t)/tdt, x >0
	 */
	/**
	 * Computes the exponential integral Ei(x) = integral form -inf to x of exp(t)/tdt, for x greater than 0.
	 * @param x x
	 * @return the result
	 */
	public static double ei(double x)
	// Computes the exponential integral Ei(x) for x > 0.
	{
		// void nrerror(char error_text[]);
		failB = false;

		int k = 0;
		double fact = 0.0;
		double prev = 0.0;
		double sum = 0.0;
		double term = 0.0;
		if (x <= 0.0) {
			// nrerror("Bad argument in ei");
			failB = true;
			failS = "Bad argument in ei";

			return 0.0;
		}
		if (x < FPMIN)
			return Math.log(x) + EULER;// Special case: avoid failure of
										// convergence
										// test because of underflow.
		if (x <= -Math.log(EPS2)) {
			sum = 0.0; // Use power series.
			fact = 1.0;
			for (k = 1; k <= MAXIT; k++) {
				fact *= x / k;
				term = fact / k;
				sum += term;
				if (term < EPS2 * sum)
					break;
			}
			if (k > MAXIT) {
				// nrerror("Series failed in ei");
				failB = true;
				failS = "Series failed in ei";

				return 0.0;
			}
			return sum + Math.log(x) + EULER;
		} else { // Use asymptotic series.
			sum = 0.0;// Start with second term.
			term = 1.0;
			for (k = 1; k <= MAXIT; k++) {
				prev = term;
				term *= k / x;
				if (term < EPS2)
					break;
				// Since final sum is greater than one, term itself approximates
				// the relative error.
				if (term < prev)
					sum += term;// Still converging: add new term.
				else {
					sum -= prev; // Diverging: subtract previous term and
									// exit.
					break;
				}
			}
			return Math.exp(x) * (1.0 + sum) / x;
		}
	}

	/*
	 * INCOMPLETE BETA FUNCTION Ix(a, b) =Bx(a, b)/B(a, b) =[1/B(a, b)]*
	 * integral from 0 to x of t^(a-1)(1 - t)^(b-1)dt (a, b > 0)
	 */
	/**
	 * Returns the incomplete beta function Ix(a, b) = integral from 0 to x of t^(a-1)(1 - t)^(b-1)dt (a, b positives).
	 * @param a a
	 * @param b b
	 * @param x x
	 * @return the result
	 */
	public static double betai(double a, double b, double x)
	// /Returns the incomplete beta function Ix(a, b).
	{
		// float betacf(float a, float b, float x);
		// float gammln(float xx);
		// void nrerror(char error_text[]);
		failB = false;
		double bt = 0.0;
		if (x < 0.0 || x > 1.0) {
			// nrerror("Bad x in routine betai");
			failB = true;
			failS = "Bad x in routine betai";

			return 0.0;
		}
		if (x == 0.0 || x == 1.0)
			bt = 0.0;
		else
			// Factors in front of the continued fraction.
			bt = Math.exp(gammln(a + b) - gammln(a) - gammln(b) + a
					* Math.log(x) + b * Math.log(1.0 - x));
		if (x < (a + 1.0) / (a + b + 2.0))// Use continued fraction directly.
			return bt * betacf(a, b, x) / a;
		else
			// Use continued fraction after making the symmetry transformation.
			return 1.0 - bt * betacf(b, a, 1.0 - x) / b;
	}

	/**
	 * Used by betai: Evaluates continued fraction for incomplete beta function by modified Lentz’s method.
	 * @param a a
	 * @param b b
	 * @param x x
	 * @return the result
	 */
	public static double betacf(double a, double b, double x)
	// Used by betai: Evaluates continued fraction for incomplete beta function
	// by modified Lentz’s
	// method (§5.2).
	{
		// void nrerror(char error_text[]);
		failB = false;
		int m = 0;
		int m2 = 0;
		double aa = 0.0;
		double c = 0.0;
		double d = 0.0;
		double del = 0.0;
		double h = 0.0;
		double qab = 0.0;
		double qam = 0.0;
		double qap = 0.0;

		qab = a + b; // These q’s will be used in factors that occur in the
						// coefficients (6.4.6).
		qap = a + 1.0;
		qam = a - 1.0;
		c = 1.0;// First step of Lentz’s method.
		d = 1.0 - qab * x / qap;
		if (Math.abs(d) < FPMIN)
			d = FPMIN;
		d = 1.0 / d;
		h = d;
		for (m = 1; m <= MAXIT; m++) {
			m2 = 2 * m;
			aa = m * (b - m) * x / ((qam + m2) * (a + m2));
			d = 1.0 + aa * d; // One step (the even one) of the recurrence.
			if (Math.abs(d) < FPMIN)
				d = FPMIN;
			c = 1.0 + aa / c;
			if (Math.abs(c) < FPMIN)
				c = FPMIN;
			d = 1.0 / d;
			h *= d * c;
			aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2));
			d = 1.0 + aa * d;// Next step of the recurrence (the odd one).
			if (Math.abs(d) < FPMIN)
				d = FPMIN;
			c = 1.0 + aa / c;
			if (Math.abs(c) < FPMIN)
				c = FPMIN;
			d = 1.0 / d;
			del = d * c;
			h *= del;
			if (Math.abs(del - 1.0) < EPS)
				break; // Are we done?
		}
		if (m > MAXIT) {
			// nrerror("a or b too big, or MAXIT too small in betacf");
			failB = true;
			failS = "a or b too big, or MAXIT too small in betacf";

			return 0.0;
		}
		return h;
	}

	// @@@@@@@Student’s Distribution Probability Function
	/*
	 * Student’s distribution, denoted A(t|?), is useful in several statistical
	 * contexts, notably in the test of whether two observed distributions have
	 * the same mean. A(t|?) is the probability, for ? degrees of freedom, that
	 * a certain statistic t (measuring the observed difference of means) would
	 * be smaller than the observed value if the means were in fact the same.
	 * (See Chapter 14 for further details.) Two means are significantly
	 * different if, e.g., A(t|?) > 0.99. In other words, 1 - A(t|?) is the
	 * significance level at which the hypothesis that the means are equal is
	 * disproved. A(t|?) = 1 - I v/v+t^2 (v/2,1/2)=> So, you can use (6.4.9) and
	 * the above routine betai to evaluate the function.
	 */

	/*
	 * F-Distribution Probability Function This function occurs in the
	 * statistical test of whether two observed samples have the same variance.
	 * A certain statistic F, essentially the ratio of the observed dispersion
	 * of the first sample to that of the second one, is calculated. (For
	 * further details, see Chapter 14.) The probability that F would be as
	 * large as it is if the first sample’s underlying distribution actually has
	 * smaller variance than the second’s is denoted Q(F|?1, ?2), where ?1 and
	 * ?2 are the number of degrees of freedom in the first and second samples,
	 * respectively. In otherwords,Q(F|?1, ?2) is the significance level at
	 * which the hypothesis “1 has smaller variance than 2” can be rejected. A
	 * small numerical value implies a very significant rejection, in turn
	 * implying high confidence in the hypothesis “1 has variance greater or
	 * equal to 2.”
	 * 
	 * Q(F|v1, v2) = I v2/v2+v1F (v2/2,1/2)
	 */

	// @@@@@@@@@@Bessel Functions of Integer Order
	/*
	 * This section and the next one present practical algorithms for computing
	 * various kinds of Bessel functions of integer order. In §6.7 we deal with
	 * fractional order. In fact, the more complicated routines for fractional
	 * order work fine for integer order too. For integer order, however, the
	 * routines in this section (and §6.6) are simpler and faster. Their only
	 * drawback is that they are limited by the precision of the underlying
	 * rational approximations. For full double precision, it is best to work
	 * with the routines for fractional order in §6.7.
	 * 
	 * For any real v, the Bessel function Jv(x) can be defined by the series
	 * representation Jv (x) = (1/2 x)^v sum from k=0 to +inf of
	 * (-x^2/4)^k/[k!gama(v + k + 1)] Yv(x) = [Jv(x) cos(vPI) - J-v(x)]/sin(vPI)
	 */
	/**
	 * Returns the Bessel function J0(x) for any real x.
	 * @param x x
	 * @return the result
	 */
	public static double bessj0(double x)
	// Returns the Bessel function J0(x) for any real x.
	{
		double ax = 0.0;
		double z = 0.0;
		double xx = 0.0;
		double y = 0.0;
		double ans = 0.0;
		double ans1 = 0.0;
		double ans2 = 0.0; // Accumulate polynomials in double precision.
		if ((ax = Math.abs(x)) < 8.0) {// Direct rational function fit.
			y = x * x;
			ans1 = 57568490574.0
					+ y
					* (-13362590354.0 + y
							* (651619640.7 + y
									* (-11214424.18 + y
											* (77392.33017 + y * (-184.9052456)))));
			ans2 = 57568490411.0
					+ y
					* (1029532985.0 + y
							* (9494680.718 + y
									* (59272.64853 + y
											* (267.8532712 + y * 1.0))));
			ans = ans1 / ans2;
		} else {// Fitting function (6.5.9).
			z = 8.0 / ax;
			y = z * z;
			xx = ax - 0.785398164;
			ans1 = 1.0
					+ y
					* (-0.1098628627e-2 + y
							* (0.2734510407e-4 + y
									* (-0.2073370639e-5 + y * 0.2093887211e-6)));
			ans2 = -0.1562499995e-1
					+ y
					* (0.1430488765e-3 + y
							* (-0.6911147651e-5 + y
									* (0.7621095161e-6 - y * 0.934945152e-7)));
			ans = Math.sqrt(0.636619772 / ax)
					* (Math.cos(xx) * ans1 - z * Math.sin(xx) * ans2);
		}
		return ans;
	}

	/**
	 * Returns the Bessel function Y0(x) for positive x.
	 * @param x x
	 * @return the result
	 */
	public static double bessy0(double x)
	// Returns the Bessel function Y0(x) for positive x.
	{
		// float bessj0(float x);
		double z = 0.0;
		double xx = 0.0;
		double y = 0.0;
		double ans = 0.0;
		double ans1 = 0.0;
		double ans2 = 0.0; // Accumulate polynomials in double precision.
		if (x < 8.0) {// Rational function approximation of (6.5.8).
			y = x * x;
			ans1 = -2957821389.0
					+ y
					* (7062834065.0 + y
							* (-512359803.6 + y
									* (10879881.29 + y
											* (-86327.92757 + y * 228.4622733))));
			ans2 = 40076544269.0
					+ y
					* (745249964.8 + y
							* (7189466.438 + y
									* (47447.26470 + y
											* (226.1030244 + y * 1.0))));
			ans = (ans1 / ans2) + 0.636619772 * bessj0(x) * Math.log(x);
		} else {// Fitting function (6.5.10).
			z = 8.0 / x;
			y = z * z;
			xx = x - 0.785398164;
			ans1 = 1.0
					+ y
					* (-0.1098628627e-2 + y
							* (0.2734510407e-4 + y
									* (-0.2073370639e-5 + y * 0.2093887211e-6)));
			ans2 = -0.1562499995e-1
					+ y
					* (0.1430488765e-3 + y
							* (-0.6911147651e-5 + y
									* (0.7621095161e-6 + y * (-0.934945152e-7))));
			ans = Math.sqrt(0.636619772 / x)
					* (Math.sin(xx) * ans1 + z * Math.cos(xx) * ans2);
		}
		return ans;
	}

	/**
	 * Returns the Bessel function J1(x) for any real x.
	 * @param x x
	 * @return the result
	 */
	public static double bessj1(double x)
	// Returns the Bessel function J1(x) for any real x.
	{
		double ax = 0.0;
		double z = 0.0;
		double xx = 0.0;
		double y = 0.0;
		double ans = 0.0;
		double ans1 = 0.0;
		double ans2 = 0.0; // Accumulate polynomials in double precision.
		if ((ax = Math.abs(x)) < 8.0) {// Direct rational approximation.
			y = x * x;
			ans1 = x
					* (72362614232.0 + y
							* (-7895059235.0 + y
									* (242396853.1 + y
											* (-2972611.439 + y
													* (15704.48260 + y
															* (-30.16036606))))));
			ans2 = 144725228442.0
					+ y
					* (2300535178.0 + y
							* (18583304.74 + y
									* (99447.43394 + y
											* (376.9991397 + y * 1.0))));
			ans = ans1 / ans2;
		} else {// Fitting function (6.5.9).
			z = 8.0 / ax;
			y = z * z;
			xx = ax - 2.356194491;
			ans1 = 1.0
					+ y
					* (0.183105e-2 + y
							* (-0.3516396496e-4 + y
									* (0.2457520174e-5 + y * (-0.240337019e-6))));
			ans2 = 0.04687499995
					+ y
					* (-0.2002690873e-3 + y
							* (0.8449199096e-5 + y
									* (-0.88228987e-6 + y * 0.105787412e-6)));
			ans = Math.sqrt(0.636619772 / ax)
					* (Math.cos(xx) * ans1 - z * Math.sin(xx) * ans2);
			if (x < 0.0)
				ans = -ans;
		}
		return ans;
	}

	/**
	 * Returns the Bessel function Y1(x) for positive x.
	 * @param x x
	 * @return the result
	 */
	public static double bessy1(double x)
	// Returns the Bessel function Y1(x) for positive x.
	{
		// float bessj1(float x);
		double z = 0.0;
		double xx = 0.0;
		double y = 0.0;
		double ans = 0.0;
		double ans1 = 0.0;
		double ans2 = 0.0; // Accumulate polynomials in double precision.
		if (x < 8.0) {// Rational function approximation of (6.5.8).
			y = x * x;
			ans1 = x
					* (-0.4900604943e13 + y
							* (0.1275274390e13 + y
									* (-0.5153438139e11 + y
											* (0.7349264551e9 + y
													* (-0.4237922726e7 + y * 0.8511937935e4)))));
			ans2 = 0.2499580570e14
					+ y
					* (0.4244419664e12 + y
							* (0.3733650367e10 + y
									* (0.2245904002e8 + y
											* (0.1020426050e6 + y
													* (0.3549632885e3 + y)))));
			ans = (ans1 / ans2) + 0.636619772
					* (bessj1(x) * Math.log(x) - 1.0 / x);
		} else {// Fitting function (6.5.10).
			z = 8.0 / x;
			y = z * z;
			xx = x - 2.356194491;
			ans1 = 1.0
					+ y
					* (0.183105e-2 + y
							* (-0.3516396496e-4 + y
									* (0.2457520174e-5 + y * (-0.240337019e-6))));
			ans2 = 0.04687499995
					+ y
					* (-0.2002690873e-3 + y
							* (0.8449199096e-5 + y
									* (-0.88228987e-6 + y * 0.105787412e-6)));
			ans = Math.sqrt(0.636619772 / x)
					* (Math.sin(xx) * ans1 + z * Math.cos(xx) * ans2);
		}
		return ans;
	}

	/*
	 * We now turn to the second task, namely how to use the recurrence formulas
	 * (6.5.6) and (6.5.7) to get the Bessel functions Jn(x) and Yn(x) for n ?
	 * 2. The latter of these is straightforward, since its upward recurrence is
	 * always stable:
	 */
	/**
	 * Returns the Bessel function Yn(x) for positive x and n greater (or equal) than 2.
	 * @param n n
	 * @param x x
	 * @return the result
	 */
	public static double bessy(int n, double x)
	// Returns the Bessel function Yn(x) for positive x and n ? 2.
	{
		// float bessy0(float x);
		// float bessy1(float x);
		// void nrerror(char error_text[]);
		failB = false;

		int j = 0;
		double by = 0.0;
		double bym = 0.0;
		double byp = 0.0;
		double tox = 0.0;
		if (n < 2) {
			// nrerror("Index n less than 2 in bessy");
			failB = true;
			failS = "Index n less than 2 in bessy";

			return 0.0;
		}
		tox = 2.0 / x;
		by = bessy1(x);// Starting values for the recurrence.
		bym = bessy0(x);
		for (j = 1; j < n; j++) {// Recurrence (6.5.7).
			byp = j * tox * by - bym;
			bym = by;
			by = byp;
		}
		return by;
	}

	/**
	 * Returns the Bessel function Jn(x) for any real x and n greater (or equal) than 2.
	 * @param n n
	 * @param x x
	 * @return the result
	 */
	public static double bessj(int n, double x)
	// Returns the Bessel function Jn(x) for any real x and n ? 2.
	{
		// float bessj0(float x);
		// float bessj1(float x);
		// void nrerror(char error_text[]);
		failB = false;

		int j = 0;
		int jsum = 0;
		int m = 0;
		double ax = 0.0;
		double bj = 0.0;
		double bjm = 0.0;
		double bjp = 0.0;
		double sum = 0.0;
		double tox = 0.0;
		double ans = 0.0;
		if (n < 2) {
			// nrerror("Index n less than 2 in bessj");
			failB = true;
			failS = "Index n less than 2 in bessj";

			return 0.0;
		}
		ax = Math.abs(x);
		if (ax == 0.0)
			return 0.0;
		else if (ax > (double) n) {// if (ax > (float) n) { Upwards recurrence
									// from J0 and J1.
			tox = 2.0 / ax;
			bjm = bessj0(ax);
			bj = bessj1(ax);
			for (j = 1; j < n; j++) {
				bjp = j * tox * bj - bjm;
				bjm = bj;
				bj = bjp;
			}
			ans = bj;
		} else {// Downwards recurrence from an even m here computed.
			tox = 2.0 / ax;
			m = 2 * ((n + (int) Math.sqrt(ACC * n)) / 2);
			jsum = 0;// jsum will alternate between 0 and 1; when it is
			// 1, we accumulate in sum the even terms in (5.5.16).
			bjp = ans = sum = 0.0;
			bj = 1.0;
			for (j = m; j > 0; j--) {// The downward recurrence.
				bjm = j * tox * bj - bjp;
				bjp = bj;
				bj = bjm;
				if (Math.abs(bj) > BIGNO) {// Renormalize to prevent overflows.
					bj *= BIGNI;
					bjp *= BIGNI;
					ans *= BIGNI;
					sum *= BIGNI;
				}
				if (jsum != 0)
					sum += bj;// if (jsum) sum += bj;// Accumulate the sum.
				if (jsum == 0)
					jsum = 1;// jsum=!jsum;// Change 0 to 1 or vice versa.
				else
					jsum = 0;
				if (j == n)
					ans = bjp;// Save the unnormalized answer.
			}
			sum = 2.0 * sum - bj;// Compute (5.5.16)
			ans /= sum; // and use it to normalize the answer.
		}
		// return x < 0.0 && (n & 1) ? -ans : ans;
		return (x < 0.0) && ((n & 1) != 0) ? -ans : ans;
	}

	/*
	 * The modified Bessel functions In(x) and Kn(x) are equivalent to the usual
	 * Bessel functions Jn and Yn evaluated for purely imaginary arguments. In
	 * detail, the relationship is
	 */
	/**
	 * Returns the modified Bessel function I0(x) for any real x.
	 * @param x x
	 * @return the result
	 */
	public static double bessi0(double x)
	// Returns the modified Bessel function I0(x) for any real x.
	{
		double ax = 0.0;
		double ans = 0.0;
		double y = 0.0;// Accumulate polynomials in double precision.
		if ((ax = Math.abs(x)) < 3.75) {// Polynomial fit.
			y = x / 3.75;
			y *= y;
			ans = 1.0
					+ y
					* (3.5156229 + y
							* (3.0899424 + y
									* (1.2067492 + y
											* (0.2659732 + y
													* (0.360768e-1 + y * 0.45813e-2)))));
		} else {
			y = 3.75 / ax;
			ans = (Math.exp(ax) / Math.sqrt(ax))
					* (0.39894228 + y
							* (0.1328592e-1 + y
									* (0.225319e-2 + y
											* (-0.157565e-2 + y
													* (0.916281e-2 + y
															* (-0.2057706e-1 + y
																	* (0.2635537e-1 + y
																			* (-0.1647633e-1 + y * 0.392377e-2))))))));
		}
		return ans;
	}

	/**
	 * Returns the modified Bessel function K0(x) for positive real x.
	 * @param x x
	 * @return the result
	 */
	public static double bessk0(double x)
	// Returns the modified Bessel function K0(x) for positive real x.
	{
		// float bessi0(float x);
		double y = 0.0;
		double ans = 0.0; // Accumulate polynomials in double precision.
		if (x <= 2.0) {// Polynomial fit.
			y = x * x / 4.0;
			ans = (-Math.log(x / 2.0) * bessi0(x))
					+ (-0.57721566 + y
							* (0.42278420 + y
									* (0.23069756 + y
											* (0.3488590e-1 + y
													* (0.262698e-2 + y
															* (0.10750e-3 + y * 0.74e-5))))));
		} else {
			y = 2.0 / x;
			ans = (Math.exp(-x) / Math.sqrt(x))
					* (1.25331414 + y
							* (-0.7832358e-1 + y
									* (0.2189568e-1 + y
											* (-0.1062446e-1 + y
													* (0.587872e-2 + y
															* (-0.251540e-2 + y * 0.53208e-3))))));
		}
		return ans;
	}

	/**
	 * Returns the modified Bessel function I1(x) for any real x.
	 * @param x x
	 * @return the result
	 */
	public static double bessi1(double x)
	// Returns the modified Bessel function I1(x) for any real x.
	{
		double ax = 0.0;
		double ans = 0.0;
		double y = 0.0;// Accumulate polynomials in double precision.
		if ((ax = Math.abs(x)) < 3.75) {// Polynomial fit.
			y = x / 3.75;
			y *= y;
			ans = ax
					* (0.5 + y
							* (0.87890594 + y
									* (0.51498869 + y
											* (0.15084934 + y
													* (0.2658733e-1 + y
															* (0.301532e-2 + y * 0.32411e-3))))));
		} else {
			y = 3.75 / ax;
			ans = 0.2282967e-1 + y
					* (-0.2895312e-1 + y * (0.1787654e-1 - y * 0.420059e-2));
			ans = 0.39894228
					+ y
					* (-0.3988024e-1 + y
							* (-0.362018e-2 + y
									* (0.163801e-2 + y
											* (-0.1031555e-1 + y * ans))));
			ans *= (Math.exp(ax) / Math.sqrt(ax));
		}
		return x < 0.0 ? -ans : ans;
	}

	/**
	 * Returns the modified Bessel function K1(x) for positive real x.
	 * @param x x
	 * @return the result
	 */
	public static double bessk1(double x)
	// Returns the modified Bessel function K1(x) for positive real x.
	{
		// float bessi1(float x);
		double y = 0.0;
		double ans = 0.0;// Accumulate polynomials in double precision.
		if (x <= 2.0) { // Polynomial fit.
			y = x * x / 4.0;
			ans = (Math.log(x / 2.0) * bessi1(x))
					+ (1.0 / x)
					* (1.0 + y
							* (0.15443144 + y
									* (-0.67278579 + y
											* (-0.18156897 + y
													* (-0.1919402e-1 + y
															* (-0.110404e-2 + y
																	* (-0.4686e-4)))))));
		} else {
			y = 2.0 / x;
			ans = (Math.exp(-x) / Math.sqrt(x))
					* (1.25331414 + y
							* (0.23498619 + y
									* (-0.3655620e-1 + y
											* (0.1504268e-1 + y
													* (-0.780353e-2 + y
															* (0.325614e-2 + y
																	* (-0.68245e-3)))))));
		}
		return ans;
	}

	/**
	 * Returns the modified Bessel function Kn(x) for positive x and n greater (or equal) than 2.
	 * @param n n
	 * @param x x
	 * @return the result
	 */
	public static double bessk(int n, double x)
	// Returns the modified Bessel function Kn(x) for positive x and n ? 2.
	{
		// float bessk0(float x);
		// float bessk1(float x);
		// void nrerror(char error_text[]);
		failB = false;
		int j = 0;
		double bk = 0.0;
		double bkm = 0.0;
		double bkp = 0.0;
		double tox = 0.0;
		if (n < 2) {
			// nrerror("Index n less than 2 in bessk");
			failB = true;
			failS = "Index n less than 2 in bessk";

			return 0.0;
		}
		tox = 2.0 / x;
		bkm = bessk0(x);// Upward recurrence for all x...
		bk = bessk1(x);
		for (j = 1; j < n; j++) {// ...and here it is.
			bkp = bkm + j * tox * bk;
			bkm = bk;
			bk = bkp;
		}
		return bk;
	}

	/**
	 * Returns the modified Bessel function In(x) for any real x and n greater (or equal) than 2.
	 * @param n n
	 * @param x x
	 * @return the result
	 */
	public static double bessi(int n, double x)
	// Returns the modified Bessel function In(x) for any real x and n ? 2.
	{
		// float bessi0(float x);
		// void nrerror(char error_text[]);
		failB = false;
		int j = 0;
		double bi = 0.0;
		double bim = 0.0;
		double bip = 0.0;
		double tox = 0.0;
		double ans = 0.0;
		if (n < 2) {
			// nrerror("Index n less than 2 in bessi");
			failB = true;
			failS = "Index n less than 2 in bessi";

			return 0.0;
		}
		if (x == 0.0)
			return 0.0;
		else {
			tox = 2.0 / Math.abs(x);
			bip = ans = 0.0;
			bi = 1.0;
			for (j = 2 * (n + (int) Math.sqrt(ACC * n)); j > 0; j--) {// Downward
																		// recurrence
																		// from
																		// even
																		// m.
				bim = bip + j * tox * bi;
				bip = bi;
				bi = bim;
				if (Math.abs(bi) > BIGNO) {// Renormalize to prevent overflows.
					ans *= BIGNI;
					bi *= BIGNI;
					bip *= BIGNI;
				}
				if (j == n)
					ans = bip;
			}
			ans *= bessi0(x) / bi;// Normalize with bessi0.
			// return x < 0.0 && (n & 1) ? -ans : ans;
			return (x < 0.0) && ((n & 1) != 0) ? -ans : ans;
		}
	}

	// @@@Bessel Functions of Fractional Order, Airy Functions, Spherical Bessel
	// Functions
	// public static void bessjy(float x, float xnu, float *rj, float *ry, float
	// *rjp, float *ryp)
	/**
	 * Returns the Bessel functions rj = Jv, ry = Yv and their derivatives rjp = J'v , ryp = Y'v, for 
	 * positive x and for xnu = v greater than 0. The relative accuracy is within one or two significant digits 
	 * of EPS, except near a zero of one of the functions, where EPS controls its absolute accuracy. 
	 * FPMIN is a number close to the machine’s smallest floating-point number. All internal arithmetic is in double precision.
	 * @param x x
	 * @param xnu xnu
	 */
	public static void bessjy(double x, double xnu)// , float *rj, float *ry,
													// float *rjp, float *ryp)
	// Returns the Bessel functions rj = Jv, ry = Yv and their derivatives rjp =
	// J'v , ryp = Y'v, for
	// positive x and for xnu = v > 0. The relative accuracy is within one or
	// two significant digits
	// of EPS, except near a zero of one of the functions, where EPS controls
	// its absolute accuracy.
	// FPMIN is a number close to the machine’s smallest floating-point number.
	// All internal arithmetic
	// is in double precision. To convert the entire routine to double
	// precision, change the float
	// declarations above to double and decrease EPS to 10-16. Also convert the
	// function beschb.
	{
		// void beschb(double x, double *gam1, double *gam2, double
		// *gampl,double *gammi);
		int i = 0;
		int isign = 0;
		int l = 0;
		int nl = 0;
		double a = 0.0;
		double b = 0.0;
		double br = 0.0;
		double bi = 0.0;
		double c = 0.0;
		double cr = 0.0;
		double ci = 0.0;
		double d = 0.0;
		double del = 0.0;
		double del1 = 0.0;
		double den = 0.0;
		double di = 0.0;
		double dlr = 0.0;
		double dli = 0.0;
		double dr = 0.0;
		double e = 0.0;
		double f = 0.0;
		double fact = 0.0;
		double fact2 = 0.0;
		double fact3 = 0.0;
		double ff = 0.0;
		double gam = 0.0;
		double gam1 = 0.0;
		double gam2 = 0.0;
		double gammi = 0.0;
		double gampl = 0.0;
		double h = 0.0;
		double p = 0.0;
		double pimu = 0.0;
		double pimu2 = 0.0;
		double q = 0.0;
		double r = 0.0;
		double rjl = 0.0;
		double rjl1 = 0.0;
		double rjmu = 0.0;
		double rjp1 = 0.0;
		double rjpl = 0.0;
		double rjtemp = 0.0;
		double ry1 = 0.0;
		double rymu = 0.0;
		double rymup = 0.0;
		double rytemp = 0.0;
		double sum = 0.0;
		double sum1 = 0.0;
		double temp = 0.0;
		double w = 0.0;
		double x2 = 0.0;
		double xi = 0.0;
		double xi2 = 0.0;
		double xmu = 0.0;
		double xmu2 = 0.0;

		failB = false;

		if (x <= 0.0 || xnu < 0.0) {
			// nrerror("bad arguments in bessjy");
			failB = true;
			failS = "bad arguments in bessjy";

			return;
		}
		// nl=(x < XMIN ? (int)(xnu+0.5) : IMAX(0,(int)(xnu-x+1.5)));
		nl = (x < XMIN ? (int) (xnu + 0.5) : Math.max(0, (int) (xnu - x + 1.5)));
		// nl is the number of downward recurrences of the J’s and upward
		// recurrences of Y ’s. xmu
		// lies between -1/2 and 1/2 for x < XMIN, while it is chosen so that x
		// is greater than the
		// turning point for x ? XMIN.
		xmu = xnu - nl;
		xmu2 = xmu * xmu;
		xi = 1.0 / x;
		xi2 = 2.0 * xi;
		w = xi2 / PI; // The Wronskian.
		isign = 1; // Evaluate CF1 by modified Lentz’s method (§5.2).
					// isign keeps track of sign changes in the denominator.
		h = xnu * xi;
		if (h < FPMINb)
			h = FPMINb;
		b = xi2 * xnu;
		d = 0.0;
		c = h;
		for (i = 1; i <= MAXITb; i++) {
			b += xi2;
			d = b - d;
			if (Math.abs(d) < FPMINb)
				d = FPMINb;
			c = b - 1.0 / c;
			if (Math.abs(c) < FPMINb)
				c = FPMINb;
			d = 1.0 / d;
			del = c * d;
			h = del * h;
			if (d < 0.0)
				isign = -isign;
			if (Math.abs(del - 1.0) < EPSb)
				break;
		}
		if (i > MAXITb) {
			// nrerror("x too large in bessjy; try asymptotic expansion");
			failB = true;
			failS = "x too large in bessjy; try asymptotic expansion";

			return;
		}
		rjl = isign * FPMINb;// Initialize Jv and J'v for downward recurrence.
		rjpl = h * rjl;
		rjl1 = rjl; // Store values for later rescaling.
		rjp1 = rjpl;
		fact = xnu * xi;
		for (l = nl; l >= 1; l--) {
			rjtemp = fact * rjl + rjpl;
			fact -= xi;
			rjpl = fact * rjtemp - rjl;
			rjl = rjtemp;
		}
		if (rjl == 0.0)
			rjl = EPSb;
		f = rjpl / rjl; // Now have unnormalized Jµ and J'µ .
		if (x < XMIN) {// Use series.
			x2 = 0.5 * x;
			pimu = PI * xmu;
			fact = (Math.abs(pimu) < EPSb ? 1.0 : pimu / Math.sin(pimu));
			d = -Math.log(x2);
			e = xmu * d;
			fact2 = (Math.abs(e) < EPSb ? 1.0 : sinh(e) / e);
			// beschb(xmu,&gam1,&gam2,&gampl,&gammi); Chebyshev evaluation of ?1
			// and ?2.
			beschb(xmu);
			gam1 = gam1b;
			gam2 = gam2b;
			gampl = gamplb;
			gammi = gammib;
			ff = 2.0 / PI * fact * (gam1 * cosh(e) + gam2 * fact2 * d); // f0.
			e = Math.exp(e);
			p = e / (gampl * PI);// p0.
			q = 1.0 / (e * PI * gammi);// q0.
			pimu2 = 0.5 * pimu;
			fact3 = (Math.abs(pimu2) < EPSb ? 1.0 : Math.sin(pimu2) / pimu2);
			r = PI * pimu2 * fact3 * fact3;
			c = 1.0;
			d = -x2 * x2;
			sum = ff + r * q;
			sum1 = p;
			for (i = 1; i <= MAXITb; i++) {
				ff = (i * ff + p + q) / (i * i - xmu2);
				c *= (d / i);
				p /= (i - xmu);
				q /= (i + xmu);
				del = c * (ff + r * q);
				sum += del;
				del1 = c * p - i * del;
				sum1 += del1;
				if (Math.abs(del) < (1.0 + Math.abs(sum)) * EPSb)
					break;
			}
			if (i > MAXITb) {
				// nrerror("bessy series failed to converge");
				failB = true;
				failS = "bessy series failed to converge";

				return;
			}
			rymu = -sum;
			ry1 = -sum1 * xi2;
			rymup = xmu * xi * rymu - ry1;
			rjmu = w / (rymup - f * rymu); // Equation (6.7.13).
		} else {// Evaluate CF2 by modified Lentz’s method (§5.2).
			a = 0.25 - xmu2;
			p = -0.5 * xi;
			q = 1.0;
			br = 2.0 * x;
			bi = 2.0;
			fact = a * xi / (p * p + q * q);
			cr = br + q * fact;
			ci = bi + p * fact;
			den = br * br + bi * bi;
			dr = br / den;
			di = -bi / den;
			dlr = cr * dr - ci * di;
			dli = cr * di + ci * dr;
			temp = p * dlr - q * dli;
			q = p * dli + q * dlr;
			p = temp;
			for (i = 2; i <= MAXITb; i++) {
				a += 2 * (i - 1);
				bi += 2.0;
				dr = a * dr + br;
				di = a * di + bi;
				if (Math.abs(dr) + Math.abs(di) < FPMINb)
					dr = FPMINb;
				fact = a / (cr * cr + ci * ci);
				cr = br + cr * fact;
				ci = bi - ci * fact;
				if (Math.abs(cr) + Math.abs(ci) < FPMINb)
					cr = FPMINb;
				den = dr * dr + di * di;
				dr /= den;
				di /= -den;
				dlr = cr * dr - ci * di;
				dli = cr * di + ci * dr;
				temp = p * dlr - q * dli;
				q = p * dli + q * dlr;
				p = temp;
				if (Math.abs(dlr - 1.0) + Math.abs(dli) < EPSb)
					break;
			}
			if (i > MAXITb) {
				// nrerror("cf2 failed in bessjy");
				failB = true;
				failS = "cf2 failed in bessjy";

				return;
			}
			gam = (p - f) / q; // Equations (6.7.6) – (6.7.10).
			rjmu = Math.sqrt(w / ((p - f) * gam + q));
			// rjmu=SIGN(rjmu,rjl);//@@@double dbl=(g >= 0.0 ? r :
			// -r);==SIGN(r,g)
			rjmu = (rjl >= 0.0 ? rjmu : -rjmu);
			rymu = rjmu * gam;
			rymup = rymu * (p + q / gam);
			ry1 = xmu * xi * rymu - rymup;
		}
		fact = rjmu / rjl;
		rjb = rjl1 * fact; // Scale original Jv and J'v .
		rjpb = rjp1 * fact;
		for (i = 1; i <= nl; i++) {// Upward recurrence of Y?.
			rytemp = (xmu + i) * xi2 * ry1 - rymu;
			rymu = ry1;
			ry1 = rytemp;
		}
		ryb = rymu;
		rypb = xnu * xi * rymu - ry1;
	}

	// public static void beschb(double x, double *gam1, double *gam2, double
	// *gampl, double *gammi)
	/**
	 * Evaluates gama1 and gama2 by Chebyshev expansion for |x| less or equal 1/2. Also returns 1/gama(1 + x) and 1/gama(1 - x).
	 * @param x x
	 */
	public static void beschb(double x)// , double *gam1, double *gam2, double
										// *gampl, double *gammi)
	// Evaluates gama1 and gama2 by Chebyshev expansion for |x| <= 1/2. Also
	// returns 1/gama(1 + x) and
	// 1/gama(1 - x). If converting to double precision, set NUSE1 = 7, NUSE2 =
	// 8.
	{
		// float chebev(float a, float b, float c[], int m, float x);
		double xx = 0.0;
		double[] c1 = { -1.142022680371168e0, 6.5165112670737e-3,
				3.087090173086e-4, -3.4706269649e-6, 6.9437664e-9, 3.67795e-11,
				-1.356e-13 };
		double[] c2 = { 1.843740587300905e0, -7.68528408447867e-2,
				1.2719271366546e-3, -4.9717367042e-6, -3.31261198e-8,
				2.423096e-10, -1.702e-13, -1.49e-15 };
		xx = 8.0 * x * x - 1.0; // Multiply x by 2 to make range be -1 to 1,
		// and then apply transformation for evaluating even Chebyshev series.
		EvalFunc ef = new EvalFunc();
		gam1b = ef.chebev(-1.0, 1.0, c1, NUSE1, xx);
		gam2b = ef.chebev(-1.0, 1.0, c2, NUSE2, xx);
		gamplb = gam2b - x * (gam1b);
		gammib = gam2b + x * (gam1b);
	}

	/*
	 * // Obtain angle in degrees from user degs = 20d; // Convert degrees to
	 * radian rads = Math.toRadians(degs);
	 * 
	 * // Calculate hyperbolic sine sinHA = (Math.exp(rads) -
	 * Math.exp(-rads))/2; System.out.println("Hyperbolic sine = " + sinHA);
	 * 
	 * // Calculate Hyperbolic cosine cosHA = (Math.exp(rads) +
	 * Math.exp(-rads))/2; System.out.println("Hyperbolic cosine = " + cosHA);
	 * 
	 * // Calculate hyperbolic tangent tanHA = sinHA/ cosHA;
	 * System.out.println("Hyperbolic tangent = " + tanHA);
	 * 
	 * // Calculate hyperbolic arc-sine asinHA = Math.log(sinHA +
	 * Math.sqrt((sinHA * sinHA)+ 1.0)); degs = Math.toDegrees(asinHA);
	 * System.out.println("Arc hyperbolic sine = " + degs);
	 */
	/**
	 * 
	 * @param degs degs
	 * @return the radians
	 */
	public static double degreeToRadian(double degs) {
		return Math.toRadians(degs);
	}

	// "Hyperbolic sine = "
	/**
	 * 
	 * @param rads input double
	 * @return the hyperbolic sine
	 */
	public static double sinh(double rads) {
		return (Math.exp(rads) - Math.exp(-rads)) / 2.0;
	}

	// "Hyperbolic cosine = "
	/**
	 * 
	 * @param rads input double
	 * @return the hyperbolic cosine
	 */
	public static double cosh(double rads) {
		return (Math.exp(rads) + Math.exp(-rads)) / 2.0;
	}

	// "Hyperbolic tangent = "
	/**
	 * 
	 * @param rads input double
	 * @return the hyperbolic tangent
	 */
	public static double tanh(double rads) {
		return sinh(rads) / cosh(rads);
	}

	// "Arc hyperbolic sine = = "
	/**
	 * This is just for testing! The output should be = input.
	 * @param rads input double
	 * @return the arc hyperbolic sine
	 */
	public static double asinh(double rads) {
		return Math
				.log(sinh(rads) + Math.sqrt((sinh(rads) * sinh(rads)) + 1.0));// rads
	}
	
	/**
	 * The real asinh function.
	 * @param rads input double which is a sinh
	 * @return its inverse
	 */
	public static double asinh_ok(double rads) {
		return Math
				.log(rads + Math.sqrt((rads * rads) + 1.0));
	}

	// Modified Bessel Functions
	// public static void bessik(float x, float xnu, float *ri, float *rk, float
	// *rip, float *rkp)
	/**
	 * Returns the modified Bessel functions ri = Iv, rk = Kv and their derivatives rip = I'v , rkp = K'v , for positive x and 
	 * for xnu = v greater (or equal) than 0. The relative accuracy is 
	 * within one or two significant digits of EPS. FPMIN is a number close to the machine’s smallest floating-point number. All internal arithmetic is in double precision. 
	 * 
	 * @param x x
	 * @param xnu xnu
	 */
	public static void bessik(double x, double xnu)// , float *ri, float *rk,
													// float *rip, float *rkp)
	// Returns the modified Bessel functions ri = Iv, rk = Kv and their
	// derivatives rip = I'v ,
	// rkp = K'v , for positive x and for xnu = v >= 0. The relative accuracy is
	// within one or two
	// significant digits of EPS. FPMIN is a number close to the machine’s
	// smallest floating-point
	// number. All internal arithmetic is in double precision. To convert the
	// entire routine to double
	// precision, change the float declarations above to double and decrease EPS
	// to 10-16. Also
	// convert the function beschb.
	{
		// void beschb(double x, double *gam1, double *gam2, double *gampl,
		// double *gammi);
		// void nrerror(char error_text[]);

		failB = false;

		int i = 0;
		int l = 0;
		int nl = 0;
		double a = 0.0;
		double a1 = 0.0;
		double b = 0.0;
		double c = 0.0;
		double d = 0.0;
		double del = 0.0;
		double del1 = 0.0;
		double delh = 0.0;
		double dels = 0.0;
		double e = 0.0;
		double f = 0.0;
		double fact = 0.0;
		double fact2 = 0.0;
		double ff = 0.0;
		double gam1 = 0.0;
		double gam2 = 0.0;
		double gammi = 0.0;
		double gampl = 0.0;
		double h = 0.0;
		double p = 0.0;
		double pimu = 0.0;
		double q = 0.0;
		double q1 = 0.0;
		double q2 = 0.0;
		double qnew = 0.0;
		double ril = 0.0;
		double ril1 = 0.0;
		double rimu = 0.0;
		double rip1 = 0.0;
		double ripl = 0.0;
		double ritemp = 0.0;
		double rk1 = 0.0;
		double rkmu = 0.0;
		double rkmup = 0.0;
		double rktemp = 0.0;
		double s = 0.0;
		double sum = 0.0;
		double sum1 = 0.0;
		double x2 = 0.0;
		double xi = 0.0;
		double xi2 = 0.0;
		double xmu = 0.0;
		double xmu2 = 0.0;

		if (x <= 0.0 || xnu < 0.0) {
			// nrerror("bad arguments in bessik");
			failB = true;
			failS = "bad arguments in bessik";

			return;
		}
		nl = (int) (xnu + 0.5);// nl is the number of downward recurrences of
								// the I’s and upward
								// recurrences of K’s. xmu lies between -1/2 and
								// 1/2.
		xmu = xnu - nl;
		xmu2 = xmu * xmu;
		xi = 1.0 / x;
		xi2 = 2.0 * xi;
		h = xnu * xi;// Evaluate CF1 by modified Lentz’s method (§5.2).
		if (h < FPMINb)
			h = FPMINb;
		b = xi2 * xnu;
		d = 0.0;
		c = h;
		for (i = 1; i <= MAXITb; i++) {
			b += xi2;
			d = 1.0 / (b + d); // Denominators cannot be zero here,so no need
								// for special precautions.
			c = b + 1.0 / c;
			del = c * d;
			h = del * h;
			if (Math.abs(del - 1.0) < EPSb)
				break;
		}
		if (i > MAXITb) {
			// nrerror("x too large in bessik; try asymptotic expansion");
			failB = true;
			failS = "x too large in bessik; try asymptotic expansion";

			return;

		}
		ril = FPMINb;// Initialize I? and I? for downward recurrence.
		ripl = h * ril;
		ril1 = ril;// Store values for later rescaling.
		rip1 = ripl;
		fact = xnu * xi;
		for (l = nl; l >= 1; l--) {
			ritemp = fact * ril + ripl;
			fact -= xi;
			ripl = fact * ritemp + ril;
			ril = ritemp;
		}
		f = ripl / ril;// Now have unnormalized Iµ and I'µ .
		if (x < XMIN) {// Use series.
			x2 = 0.5 * x;
			pimu = PI * xmu;
			fact = (Math.abs(pimu) < EPSb ? 1.0 : pimu / Math.sin(pimu));
			d = -Math.log(x2);
			e = xmu * d;
			fact2 = (Math.abs(e) < EPSb ? 1.0 : sinh(e) / e);
			// beschb(xmu,&gam1,&gam2,&gampl,&gammi); Chebyshev evaluation of ?1
			// and ?2.
			beschb(xmu);
			gam1 = gam1b;
			gam2 = gam2b;
			gampl = gamplb;
			gammi = gammib;
			ff = fact * (gam1 * cosh(e) + gam2 * fact2 * d); // f0.
			sum = ff;
			e = Math.exp(e);
			p = 0.5 * e / gampl;// p0.
			q = 0.5 / (e * gammi); // q0.
			c = 1.0;
			d = x2 * x2;
			sum1 = p;
			for (i = 1; i <= MAXITb; i++) {
				ff = (i * ff + p + q) / (i * i - xmu2);
				c *= (d / i);
				p /= (i - xmu);
				q /= (i + xmu);
				del = c * ff;
				sum += del;
				del1 = c * (p - i * ff);
				sum1 += del1;
				if (Math.abs(del) < Math.abs(sum) * EPSb)
					break;
			}
			if (i > MAXITb) {
				// nrerror("bessk series failed to converge");
				failB = true;
				failS = "bessk series failed to converge";

				return;
			}
			rkmu = sum;
			rk1 = sum1 * xi2;
		} else {// Evaluate CF2 by Steed’s algorithm (§5.2), which is OK because
				// there
				// can be no zero denominators.
			b = 2.0 * (1.0 + x);
			d = 1.0 / b;
			h = delh = d;
			q1 = 0.0; // Initializations for recurrence (6.7.35).
			q2 = 1.0;
			a1 = 0.25 - xmu2;
			q = c = a1;// First term in equation (6.7.34).
			a = -a1;
			s = 1.0 + q * delh;
			for (i = 2; i <= MAXITb; i++) {
				a -= 2 * (i - 1);
				c = -a * c / i;
				qnew = (q1 - b * q2) / a;
				q1 = q2;
				q2 = qnew;
				q += c * qnew;
				b += 2.0;
				d = 1.0 / (b + a * d);
				delh = (b * d - 1.0) * delh;
				h += delh;
				dels = q * delh;
				s += dels;
				if (Math.abs(dels / s) < EPSb)
					break;
				// Need only test convergence of sum since CF2 itself converges
				// more quickly.
			}
			if (i > MAXITb) {
				// nrerror("bessik: failure to converge in cf2");
				failB = true;
				failS = "bessik: failure to converge in cf2";

				return;
			}
			h = a1 * h;
			rkmu = Math.sqrt(PI / (2.0 * x)) * Math.exp(-x) / s; // Omit the
																	// factor
																	// exp(-x)
																	// to scale
			// all the returned functions by exp(x) for x >= XMIN.
			rk1 = rkmu * (xmu + x + 0.5 - h) * xi;
		}
		rkmup = xmu * xi * rkmu - rk1;
		rimu = xi / (f * rkmu - rkmup); // Get Iµ from Wronskian.
		rib = (rimu * ril1) / ril;// Scale original Iv and I'v .
		ripb = (rimu * rip1) / ril;
		for (i = 1; i <= nl; i++) {// Upward recurrence of Kv.
			rktemp = (xmu + i) * xi2 * rk1 + rkmu;
			rkmu = rk1;
			rk1 = rktemp;
		}
		rkb = rkmu;
		rkpb = xnu * xi * rkmu - rk1;
	}

	// @@@@@@@Airy Functions
	// void airy(float x, float *ai, float *bi, float *aip, float *bip)
	// Returns Airy functions Ai(x), Bi(x), and their derivatives Ai'(x),
	// Bi'(x).
	/**
	 * Returns Airy functions Ai(x), Bi(x), and their derivatives Ai'(x), Bi'(x).
	 * @param x x
	 */
	public static void airy(double x)// , float *ai, float *bi, float *aip,
										// float *bip)
	{
		// void bessik(float x, float xnu, float *ri, float *rk, float *rip,
		// float *rkp);
		// void bessjy(float x, float xnu, float *rj, float *ry, float *rjp,
		// float *ryp);
		double absx = 0.0;
		double ri = 0.0;
		// double rip = 0.0;
		double rj = 0.0;
		// double rjp = 0.0;
		double rk = 0.0;
		// double rkp = 0.0;
		double rootx = 0.0;
		double ry = 0.0;
		// double ryp = 0.0;
		double z = 0.0;

		absx = Math.abs(x);
		rootx = Math.sqrt(absx);
		z = TWOTHR * absx * rootx;
		if (x > 0.0) {
			// bessik(z,THIRD,&ri,&rk,&rip,&rkp);
			bessik(z, THIRD);
			ri = rib;
			rk = rkb;
			// rip = ripb;
			// rkp = rkpb;

			ai_airy = rootx * ONOVRT * rk / PI;
			bi_airy = rootx * (rk / PI + 2.0 * ONOVRT * ri);

			// bessik(z,TWOTHR,&ri,&rk,&rip,&rkp);
			bessik(z, TWOTHR);
			ri = rib;
			rk = rkb;
			// rip = ripb;
			// rkp = rkpb;

			aip_airy = -x * ONOVRT * rk / PI;
			bip_airy = x * (rk / PI + 2.0 * ONOVRT * ri);
		} else if (x < 0.0) {
			// bessjy(z,THIRD,&rj,&ry,&rjp,&ryp);
			bessjy(z, THIRD);
			rj = rjb;
			ry = ryb;
			// rjp = rjpb;
			// ryp = rypb;

			ai_airy = 0.5 * rootx * (rj - ONOVRT * ry);
			bi_airy = -0.5 * rootx * (ry + ONOVRT * rj);

			// bessjy(z,TWOTHR,&rj,&ry,&rjp,&ryp);
			bessjy(z, TWOTHR);
			rj = rjb;
			ry = ryb;
			// rjp = rjpb;
			// ryp = rypb;

			aip_airy = 0.5 * absx * (ONOVRT * ry + rj);
			bip_airy = 0.5 * absx * (ONOVRT * rj - ry);
		} else {// Case x = 0.
			ai_airy = 0.35502805;
			bi_airy = (ai_airy) / ONOVRT;
			aip_airy = -0.25881940;
			bip_airy = -(aip_airy) / ONOVRT;
		}
	}

	// Spherical Bessel Functions
	/**
	 * Returns spherical Bessel functions jn(x), yn(x), and their derivatives j'n (x), y'n (x) for integer n.
	 * @param n n
	 * @param x x
	 */
	public static void sphbes(int n, double x)// , float *sj, float *sy, float
												// *sjp, float *syp)
	// Returns spherical Bessel functions jn(x), yn(x), and their derivatives
	// j'n (x), y'n (x) for integer n.
	{
		// void bessjy(float x, float xnu, float *rj, float *ry, float
		// *rjp,float *ryp);
		// void nrerror(char error_text[]);
		failB = false;
		double factor = 0.0;
		double order = 0.0;
		double rj = 0.0;
		double rjp = 0.0;
		double ry = 0.0;
		double ryp = 0.0;
		if (n < 0 || x <= 0.0) {
			// nrerror("bad arguments in sphbes");
			failB = true;
			failS = "bad arguments in sphbes";

			return;
		}
		order = n + 0.5;
		// bessjy(x,order,&rj,&ry,&rjp,&ryp);
		bessjy(x, order);
		rj = rjb;
		ry = ryb;
		rjp = rjpb;
		ryp = rypb;
		factor = RTPIO2 / Math.sqrt(x);
		sjb = factor * rj;
		syb = factor * ry;
		sjpb = factor * rjp - (sjb) / (2.0 * x);
		sypb = factor * ryp - (syb) / (2.0 * x);
	}

	// @@@@@@@@Spherical Harmonics
	/**
	 * Computes the associated Legendre polynomial Plm (x). Here m and l are integers satisfying 
	 * m is in [0,l] while x is in the range [-1,1]. Useful in quantum mechanics when describing atoms (wave equation is solved by 
	 * separation of variables in spherical coordinates).
	 * @param l l
	 * @param m m
	 * @param x x
	 * @return the result
	 */
	public static double plgndr(int l, int m, double x)
	// Computes the associated Legendre polynomial Plm (x). Here m and l are
	// integers satisfying
	// 0 <= m <= l, while x lies in the range -1 <= x <= 1.
	{
		// void nrerror(char error_text[]);
		failB = false;
		double fact = 0.0;
		double pll = 0.0;
		double pmm = 0.0;
		double pmmp1 = 0.0;
		double somx2 = 0.0;
		int i = 0;
		int ll = 0;
		if (m < 0 || m > l || Math.abs(x) > 1.0) {
			// nrerror("Bad arguments in routine plgndr");
			failB = true;
			failS = "Bad arguments in routine plgndr";

			return 0.0;
		}
		pmm = 1.0; // Compute Pmm.
		if (m > 0) {
			somx2 = Math.sqrt((1.0 - x) * (1.0 + x));
			fact = 1.0;
			for (i = 1; i <= m; i++) {
				pmm *= -fact * somx2;
				fact += 2.0;
			}
		}
		if (l == m)
			return pmm;
		else {// Compute Pm m+1.
			pmmp1 = x * (2 * m + 1) * pmm;
			if (l == (m + 1))
				return pmmp1;
			else {// Compute Pm l , l > m+ 1.
				for (ll = m + 2; ll <= l; ll++) {
					pll = (x * (2 * ll - 1) * pmmp1 - (ll + m - 1) * pmm)
							/ (ll - m);
					pmm = pmmp1;
					pmmp1 = pll;
				}
				return pll;
			}
		}
	}

	// @@Fresnel Integrals, Cosine and Sine Integrals
	/*
	 * C(x) = integral from 0 to x of cos(PIt^2/2)dt, S(x) = integral from 0 to
	 * x of sin(PIt^2/2)dt,
	 */
	/**
	 * Computes the Fresnel integrals S(x) and C(x) for all real x. 
	 * @param x x
	 */
	public static void frenel(double x)// , float *s, float *c)
	// Computes the Fresnel integrals S(x) and C(x) for all real x.
	{
		// void nrerror(char error_text[]);
		failB = false;
		int k = 0;
		int n = 0;
		int odd = 0;
		double a = 0.0;
		double ax = 0.0;
		double fact = 0.0;
		double pix2 = 0.0;
		double sign = 0.0;
		double sum = 0.0;
		double sumc = 0.0;
		double sums = 0.0;
		double term = 0.0;
		double test = 0.0;
		Complex b, cc, d, h, del, cs;

		ax = Math.abs(x);
		if (ax < Math.sqrt(FPMIN)) {// Special case: avoid failure of
									// convergence test because of underflow.
			s_fresnel = 0.0;
			c_fresnel = ax;
		} else if (ax <= XMIN2) {// Evaluate both series simultaneously.
			sum = sums = 0.0;
			sumc = ax;
			sign = 1.0;
			fact = PIBY2 * ax * ax;
			odd = TRUE;
			term = ax;
			n = 3;
			for (k = 1; k <= MAXIT; k++) {
				term *= fact / k;
				sum += sign * term / n;
				test = Math.abs(sum) * EPS2;
				if (odd != 0) // if (odd)
				{
					sign = -sign;
					sums = sum;
					sum = sumc;
				} else {
					sumc = sum;
					sum = sums;
				}
				if (term < test)
					break;
				// odd=!odd;
				if (odd == 0)
					odd = 1;
				else
					odd = 0;
				n += 2;
			}
			if (k > MAXIT) {
				// nrerror("series failed in frenel");
				failB = true;
				failS = "series failed in frenel";

				return;
			}
			s_fresnel = sums;
			c_fresnel = sumc;
		} else {// Evaluate continued fraction by modified Lentz’s method
				// (§5.2).
			pix2 = PI * ax * ax;
			b = new Complex(1.0, -pix2);
			cc = new Complex(1.0 / FPMIN, 0.0);
			d = h = Complex.Cdiv(ONE, b);
			n = -1;
			for (k = 2; k <= MAXIT; k++) {
				n += 2;
				a = -n * (n + 1);
				b = Complex.Cadd(b, new Complex(4.0, 0.0));
				d = Complex.Cdiv(ONE, Complex.Cadd(Complex.RCmul(a, d), b));// Denominators
																			// cannot
																			// be
																			// zero.
				cc = Complex.Cadd(b, Complex.Cdiv(new Complex(a, 0.0), cc));
				del = Complex.Cmul(cc, d);
				h = Complex.Cmul(h, del);
				if (Math.abs(del.r - 1.0) + Math.abs(del.i) < EPS2)
					break;
			}
			if (k > MAXIT) {
				// nrerror("cf failed in frenel");
				failB = true;
				failS = "cf failed in frenel";

				return;
			}
			h = Complex.Cmul(new Complex(ax, -ax), h);
			cs = Complex.Cmul(new Complex(0.5, 0.5), Complex.Csub(ONE, Complex
					.Cmul(new Complex(Math.cos(0.5 * pix2), Math
							.sin(0.5 * pix2)), h)));
			c_fresnel = cs.r;
			s_fresnel = cs.i;
		}
		if (x < 0.0) {// Use antisymmetry.
			c_fresnel = -(c_fresnel);
			s_fresnel = -(s_fresnel);
		}
	}

	// @@@Cosine and Sine Integrals
	/*
	 * Ci(x) = gamalowercase + lnx + integral 0-x of [(cos t - 1)/t]dt Si(x) =
	 * integral 0-x of [sint)/t]dt
	 */
	/**
	 * Computes the cosine and sine integrals Ci(x) and Si(x). Ci(0) is returned as a large negative 
	 * number and no error message is generated. For x less than 0 the routine returns Ci(-x) and you must supply the -iPI yourself.
	 * @param x x
	 */
	public static void cisi(double x)// , float *ci, float *si)
	// Computes the cosine and sine integrals Ci(x) and Si(x). Ci(0) is returned
	// as a large negative
	// number and no error message is generated. For x < 0 the routine returns
	// Ci(-x) and you must
	// supply the -ip yourself.
	{
		// void nrerror(char error_text[]);
		failB = false;
		int i = 0;
		int k = 0;
		int odd = 0;
		double a = 0.0;
		double err = 0.0;
		double fact = 0.0;
		double sign = 0.0;
		double sum = 0.0;
		double sumc = 0.0;
		double sums = 0.0;
		double t = 0.0;
		double term = 0.0;
		Complex h, b, c, d, del;

		t = Math.abs(x);
		if (t == 0.0) {// Special case.
			si_integral = 0.0;
			ci_integral = -1.0 / FPMIN;
			return;
		}
		if (t > TMIN) {// Evaluate continued fraction by modified Lentz’s method
						// (§5.2).
			b = new Complex(1.0, t);
			c = new Complex(1.0 / FPMIN, 0.0);
			d = h = Complex.Cdiv(ONE, b);
			for (i = 2; i <= MAXIT; i++) {
				a = -(i - 1) * (i - 1);
				b = Complex.Cadd(b, new Complex(2.0, 0.0));
				d = Complex.Cdiv(ONE, Complex.Cadd(Complex.RCmul(a, d), b)); // Denominators
																				// cannot
																				// be
																				// zero.
				c = Complex.Cadd(b, Complex.Cdiv(new Complex(a, 0.0), c));
				del = Complex.Cmul(c, d);
				h = Complex.Cmul(h, del);
				if (Math.abs(del.r - 1.0) + Math.abs(del.i) < EPS2)
					break;
			}
			if (i > MAXIT) {
				// nrerror("cf failed in cisi");
				failB = true;
				failS = "cf failed in cisi";

				return;
			}
			h = Complex.Cmul(new Complex(Math.cos(t), -Math.sin(t)), h);
			ci_integral = -h.r;
			si_integral = PIBY2 + h.i;
		} else {// Evaluate both series simultaneously.
			if (t < Math.sqrt(FPMIN)) {// Special case: avoid failure of
										// convergence test because of
										// underflow.
				sumc = 0.0;
				sums = t;
			} else {
				sum = sums = sumc = 0.0;
				sign = fact = 1.0;
				odd = TRUE;
				for (k = 1; k <= MAXIT; k++) {
					fact *= t / k;
					term = fact / k;
					sum += sign * term;
					err = term / Math.abs(sum);
					if (odd != 0)// if (odd)<=>if (sum!=0.0) ii=i;//if (sum)
									// ii=i;
					{
						sign = -sign;
						sums = sum;
						sum = sumc;
					} else {
						sumc = sum;
						sum = sums;
					}
					if (err < EPS2)
						break;
					// odd=!odd;
					if (odd == 0)
						odd = 1;
					else
						odd = 0;
				}
				if (k > MAXIT) {
					// nrerror("maxits exceeded in cisi");
					failB = true;
					failS = "maxits exceeded in cisi";

					return;
				}
			}
			si_integral = sums;
			ci_integral = sumc + Math.log(t) + EULER;
		}
		if (x < 0.0)
			si_integral = -(si_integral);
	}

	// @@Dawson’s Integral
	/*
	 * F(x) = exp(-x2) Integral0-x of exp(t^2)dt
	 */
	/**
	 * Returns Dawson’s integral for any real x.
	 * @param x x
	 * @return the result
	 */
	public static double dawson(double x)
	// Returns Dawson’s integral for any real x.
	{
		int i = 0;
		int n0 = 0;
		double d1 = 0.0;
		double d2 = 0.0;
		double e1 = 0.0;
		double e2 = 0.0;
		double sum = 0.0;
		double x2 = 0.0;
		double xp = 0.0;
		double xx = 0.0;
		double ans = 0.0;
		// static float c[NMAX+1];
		// static int init = 0; Flag is 0 if we need to initialize, else 1.
		if (init == 0) {
			init = 1;
			for (i = 1; i <= NMAX; i++)
				c[i] = Math.exp(-((2.0 * i - 1.0) * H) * ((2.0 * i - 1.0) * H));// Math.exp(-SQR((2.0*i-1.0)*H));
		}
		if (Math.abs(x) < 0.2) {// Use series expansion.
			x2 = x * x;
			ans = x * (1.0 - A1 * x2 * (1.0 - A2 * x2 * (1.0 - A3 * x2)));
		} else {// Use sampling theorem representation.
			xx = Math.abs(x);
			n0 = 2 * (int) (0.5 * xx / H + 0.5);
			xp = xx - n0 * H;
			e1 = Math.exp(2.0 * xp * H);
			e2 = e1 * e1;
			d1 = n0 + 1;
			d2 = d1 - 2.0;
			sum = 0.0;
			for (i = 1; i <= NMAX; i++, d1 += 2.0, d2 -= 2.0, e1 *= e2)
				sum += c[i] * (e1 / d1 + 1.0 / (d2 * e1));
			// ans=0.5641895835*SIGN(exp(-xp*xp),x)*sum; Constant is 1/vp.
			// //@@@double dbl=(g >= 0.0 ? r : -r);==SIGN(r,g)
			double dbl = (x >= 0.0 ? Math.exp(-xp * xp) : -Math.exp(-xp * xp));
			ans = 0.5641895835 * dbl * sum;
		}
		return ans;
	}

	// @@Elliptic Integrals and Jacobian Elliptic Functions
	/**
	 * Computes Carlson’s elliptic integral of the first kind, RF (x, y, z). x, y, and z must be nonnegative, 
	 * and at most one can be zero. TINY must be at least 5 times the machine underflow limit, BIG at most one fifth the machine overflow limit. 
	 * @param x x
	 * @param y y
	 * @param z z
	 * @return the result
	 */	
	public static double rf(double x, double y, double z)
	// Computes Carlson’s elliptic integral of the first kind, RF (x, y, z). x,
	// y, and z must be nonnegative,
	// and at most one can be zero. TINY must be at least 5 times the machine
	// underflow limit,
	// BIG at most one fifth the machine overflow limit.
	{
		double alamb = 0.0;
		double ave = 0.0;
		double delx = 0.0;
		double dely = 0.0;
		double delz = 0.0;
		double e2 = 0.0;
		double e3 = 0.0;
		double sqrtx = 0.0;
		double sqrty = 0.0;
		double sqrtz = 0.0;
		double xt = 0.0;
		double yt = 0.0;
		double zt = 0.0;

		failB = false;

		if (Math.min(Math.min(x, y), z) < 0.0
				|| Math.min(Math.min(x + y, x + z), y + z) < TINY
				|| Math.max(Math.max(x, y), z) > BIG) {
			// nrerror("invalid arguments in rf");
			failB = true;
			failS = "invalid arguments in rf";

			return 0.0;
		}
		xt = x;
		yt = y;
		zt = z;
		do {
			sqrtx = Math.sqrt(xt);
			sqrty = Math.sqrt(yt);
			sqrtz = Math.sqrt(zt);
			alamb = sqrtx * (sqrty + sqrtz) + sqrty * sqrtz;
			xt = 0.25 * (xt + alamb);
			yt = 0.25 * (yt + alamb);
			zt = 0.25 * (zt + alamb);
			ave = THIRD * (xt + yt + zt);
			delx = (ave - xt) / ave;
			dely = (ave - yt) / ave;
			delz = (ave - zt) / ave;
		} while (Math.max(Math.max(Math.abs(delx), Math.abs(dely)),
				Math.abs(delz)) > ERRTOL);
		e2 = delx * dely - delz * delz;
		e3 = delx * dely * delz;
		return (1.0 + (C1 * e2 - C2 - C3 * e3) * e2 + C4 * e3) / Math.sqrt(ave);
	}

	/**
	 * Computes Carlson’s elliptic integral of the second kind, RD(x, y, z). x and y must be nonnegative, 
	 * and at most one can be zero. z must be positive. TINY must be at least twice the 
	 * negative 2/3 power of the machine overflow limit. BIG must be at most 0.1 × ERRTOL times the negative 2/3 power of the machine underflow limit.
	 * @param x x
	 * @param y y
	 * @param z z
	 * @return the result
	 */
	public static double rd(double x, double y, double z)
	// Computes Carlson’s elliptic integral of the second kind, RD(x, y, z). x
	// and y must be nonnegative,
	// and at most one can be zero. z must be positive. TINY must be at least
	// twice the
	// negative 2/3 power of the machine overflow limit. BIG must be at most 0.1
	// × ERRTOL times
	// the negative 2/3 power of the machine underflow limit.
	{
		double alamb = 0.0;
		double ave = 0.0;
		double delx = 0.0;
		double dely = 0.0;
		double delz = 0.0;
		double ea = 0.0;
		double eb = 0.0;
		double ec = 0.0;
		double ed = 0.0;
		double ee = 0.0;
		double fac = 0.0;
		double sqrtx = 0.0;
		double sqrty = 0.0;
		double sqrtz = 0.0;
		double sum = 0.0;
		double xt = 0.0;
		double yt = 0.0;
		double zt = 0.0;

		failB = false;

		if (Math.min(x, y) < 0.0 || Math.min(x + y, z) < TINY1
				|| Math.max(Math.max(x, y), z) > BIG1) {
			// nrerror("invalid arguments in rd");
			failB = true;
			failS = "invalid arguments in rd";

			return 0.0;
		}
		xt = x;
		yt = y;
		zt = z;
		sum = 0.0;
		fac = 1.0;
		do {
			sqrtx = Math.sqrt(xt);
			sqrty = Math.sqrt(yt);
			sqrtz = Math.sqrt(zt);
			alamb = sqrtx * (sqrty + sqrtz) + sqrty * sqrtz;
			sum += fac / (sqrtz * (zt + alamb));
			fac = 0.25 * fac;
			xt = 0.25 * (xt + alamb);
			yt = 0.25 * (yt + alamb);
			zt = 0.25 * (zt + alamb);
			ave = 0.2 * (xt + yt + 3.0 * zt);
			delx = (ave - xt) / ave;
			dely = (ave - yt) / ave;
			delz = (ave - zt) / ave;
		} while (Math.max(Math.max(Math.abs(delx), Math.abs(dely)),
				Math.abs(delz)) > ERRTOL1);
		ea = delx * dely;
		eb = delz * delz;
		ec = ea - eb;
		ed = ea - 6.0 * eb;
		ee = ed + ec + ec;
		return 3.0
				* sum
				+ fac
				* (1.0 + ed * (-C11 + C51 * ed - C61 * delz * ee) + delz
						* (C21 * ee + delz * (-C31 * ec + delz * C41 * ea)))
				/ (ave * Math.sqrt(ave));
	}

	/**
	 * Computes Carlson’s elliptic integral of the third kind, RJ (x, y, z, p). x, y, and z must be 
	 * nonnegative, and at most one can be zero. p must be nonzero. If p less than 0, the Cauchy principal value is returned. TINY must be at least twice the cube root of the 
	 * machine underflow limit, BIG at most one fifth the cube root of the machine overflow limit.
	 * @param x x
	 * @param y y
	 * @param z z
	 * @param p p
	 * @return the result
	 */
	public static double rj(double x, double y, double z, double p)
	// Computes Carlson’s elliptic integral of the third kind, RJ (x, y, z, p).
	// x, y, and z must be
	// nonnegative, and at most one can be zero. p must be nonzero. If p < 0,
	// the Cauchy principal
	// value is returned. TINY must be at least twice the cube root of the
	// machine underflow limit,
	// BIG at most one fifth the cube root of the machine overflow limit.
	{
		// float rc(float x, float y);
		// float rf(float x, float y, float z);
		failB = false;
		double a = 0.0;
		double alamb = 0.0;
		double alpha = 0.0;
		double ans = 0.0;
		double ave = 0.0;
		double b = 0.0;
		double beta = 0.0;
		double delp = 0.0;
		double delx = 0.0;
		double dely = 0.0;
		double delz = 0.0;
		double ea = 0.0;
		double eb = 0.0;
		double ec = 0.0;
		double ed = 0.0;
		double ee = 0.0;
		double fac = 0.0;
		double pt = 0.0;
		double rcx = 0.0;
		double rho = 0.0;
		double sqrtx = 0.0;
		double sqrty = 0.0;
		double sqrtz = 0.0;
		double sum = 0.0;
		double tau = 0.0;
		double xt = 0.0;
		double yt = 0.0;
		double zt = 0.0;

		if (Math.min(Math.min(x, y), z) < 0.0
				|| Math.min(Math.min(Math.min(x + y, x + z), y + z),
						Math.abs(p)) < TINY2
				|| Math.max(Math.max(Math.max(x, y), z), Math.abs(p)) > BIG2) {
			// nrerror("invalid arguments in rj");
			failB = true;
			failS = "invalid arguments in rj";

			return 0.0;
		}
		sum = 0.0;
		fac = 1.0;
		if (p > 0.0) {
			xt = x;
			yt = y;
			zt = z;
			pt = p;
		} else {
			xt = Math.min(Math.min(x, y), z);
			zt = Math.max(Math.max(x, y), z);
			yt = x + y + z - xt - zt;
			a = 1.0 / (yt - p);
			b = a * (zt - yt) * (yt - xt);
			pt = yt + b;
			rho = xt * zt / yt;
			tau = p * pt / yt;
			rcx = rc(rho, tau);
		}
		do {
			sqrtx = Math.sqrt(xt);
			sqrty = Math.sqrt(yt);
			sqrtz = Math.sqrt(zt);
			alamb = sqrtx * (sqrty + sqrtz) + sqrty * sqrtz;
			// alpha=SQR(pt*(sqrtx+sqrty+sqrtz)+sqrtx*sqrty*sqrtz);
			alpha = (pt * (sqrtx + sqrty + sqrtz) + sqrtx * sqrty * sqrtz)
					* (pt * (sqrtx + sqrty + sqrtz) + sqrtx * sqrty * sqrtz);
			beta = pt * (pt + alamb) * (pt + alamb);// beta=pt*SQR(pt+alamb);
			sum += fac * rc(alpha, beta);
			fac = 0.25 * fac;
			xt = 0.25 * (xt + alamb);
			yt = 0.25 * (yt + alamb);
			zt = 0.25 * (zt + alamb);
			pt = 0.25 * (pt + alamb);
			ave = 0.2 * (xt + yt + zt + pt + pt);
			delx = (ave - xt) / ave;
			dely = (ave - yt) / ave;
			delz = (ave - zt) / ave;
			delp = (ave - pt) / ave;
		} while (Math.max(
				Math.max(Math.max(Math.abs(delx), Math.abs(dely)),
						Math.abs(delz)), Math.abs(delp)) > ERRTOL2);
		ea = delx * (dely + delz) + dely * delz;
		eb = delx * dely * delz;
		ec = delp * delp;
		ed = ea - 3.0 * ec;
		ee = eb + 2.0 * delp * (ea - ec);
		ans = 3.0
				* sum
				+ fac
				* (1.0 + ed * (-C12 + C52 * ed - C62 * ee) + eb
						* (C72 + delp * (-C82 + delp * C42)) + delp * ea
						* (C22 - delp * C32) - C22 * delp * ec)
				/ (ave * Math.sqrt(ave));
		if (p <= 0.0)
			ans = a * (b * ans + 3.0 * (rcx - rf(xt, yt, zt)));
		return ans;
	}

	/**
	 * Computes Carlson’s degenerate elliptic integral, RC(x, y). x must be nonnegative and y must 
	 * be nonzero. If y less than 0, the Cauchy principal value is returned. TINY must be at least 5 times 
	 * the machine underflow limit, BIG at most one fifth the machine maximum overflow limit.
	 * @param x x
	 * @param y y
	 * @return the result
	 */
	public static double rc(double x, double y)
	// Computes Carlson’s degenerate elliptic integral, RC(x, y). x must be
	// nonnegative and y must
	// be nonzero. If y < 0, the Cauchy principal value is returned. TINY must
	// be at least 5 times
	// the machine underflow limit, BIG at most one fifth the machine maximum
	// overflow limit.
	{
		double alamb = 0.0;
		double ave = 0.0;
		double s = 0.0;
		double w = 0.0;
		double xt = 0.0;
		double yt = 0.0;

		failB = false;

		if (x < 0.0 || y == 0.0 || (x + Math.abs(y)) < TINY3
				|| (x + Math.abs(y)) > BIG3
				|| (y < -COMP1 && x > 0.0 && x < COMP2)) {
			// nrerror("invalid arguments in rc");
			failB = true;
			failS = "invalid arguments in rc";

			return 0.0;
		}
		if (y > 0.0) {
			xt = x;
			yt = y;
			w = 1.0;
		} else {
			xt = x - y;
			yt = -y;
			w = Math.sqrt(x) / Math.sqrt(xt);
		}
		do {
			alamb = 2.0 * Math.sqrt(xt) * Math.sqrt(yt) + yt;
			xt = 0.25 * (xt + alamb);
			yt = 0.25 * (yt + alamb);
			ave = THIRD * (xt + yt + yt);
			s = (yt - ave) / ave;
		} while (Math.abs(s) > ERRTOL3);
		return w * (1.0 + s * s * (C13 + s * (C23 + s * (C33 + s * C43))))
				/ Math.sqrt(ave);
	}

	/**
	 * Legendre elliptic integral of the 1st kind F(phi, k), evaluated using Carlson’s function RF. The 
	 * argument ranges are: PHI in [0, PI/2], k x sin(PHI) in [0,1].
	 * @param phi phi
	 * @param ak ak
	 * @return the result
	 */
	public static double ellf(double phi, double ak)
	// Legendre elliptic integral of the 1st kind F(f, k), evaluated using
	// Carlson’s function RF. The
	// argument ranges are 0 = f = p/2, 0 = k sin f = 1.
	{
		// float rf(float x, float y, float z);
		double s = 0.0;
		s = Math.sin(phi);
		return s
				* rf((Math.cos(phi)) * (Math.cos(phi)), (1.0 - s * ak)
						* (1.0 + s * ak), 1.0);
	}

	/**
	 * Legendre elliptic integral of the 2nd kind E(PHI, k), evaluated using Carlson’s functions RD and 
	 * RF . The argument ranges are PHI in [0, PI/2], k x sin(PHI) in [0,1].
	 * @param phi phi
	 * @param ak ak
	 * @return the result
	 */
	public static double elle(double phi, double ak)
	// Legendre elliptic integral of the 2nd kind E(f, k), evaluated using
	// Carlson’s functions RD and
	// RF . The argument ranges are 0 = f = p/2, 0 = k sin f = 1.
	{
		// float rd(float x, float y, float z);
		// float rf(float x, float y, float z);
		double cc = 0.0;
		double q = 0.0;
		double s = 0.0;
		s = Math.sin(phi);
		cc = (Math.cos(phi)) * (Math.cos(phi));// SQR(cos(phi));
		q = (1.0 - s * ak) * (1.0 + s * ak);
		// return s*(rf(cc,q,1.0)-(SQR(s*ak))*rd(cc,q,1.0)/3.0);
		return s
				* (rf(cc, q, 1.0) - ((s * ak) * (s * ak)) * rd(cc, q, 1.0)
						/ 3.0);
	}

	/**
	 * Legendre elliptic integral of the 3rd kind PI(PHI, n, k), evaluated using Carlson’s functions RJ and 
	 * RF . (Note that the sign convention on n is opposite that of Abramowitz and Stegun.) The 
	 * argument ranges are PHI in [0, PI/2], k x sin(PHI) in [0,1]. 
	 * @param phi phi
	 * @param en en
	 * @param ak ak
	 * @return the result
	 */
	public static double ellpi(double phi, double en, double ak)
	// Legendre elliptic integral of the 3rd kind PI(PHI, n, k), evaluated using
	// Carlson’s functions RJ and
	// RF . (Note that the sign convention on n is opposite that of Abramowitz
	// and Stegun.) The
	// ranges of f and k are 0 = f = p/2, 0 = k sin f = 1.
	{
		// float rf(float x, float y, float z);
		// float rj(float x, float y, float z, float p);
		double cc = 0.0;
		double enss = 0.0;
		double q = 0.0;
		double s = 0.0;
		s = Math.sin(phi);
		enss = en * s * s;
		cc = (Math.cos(phi)) * (Math.cos(phi));// SQR(cos(phi));
		q = (1.0 - s * ak) * (1.0 + s * ak);
		// return s*(rf(cc,q,1.0)-enss*rj(cc,q,1.0,1.0+enss)/3.0);
		return s * (rf(cc, q, 1.0) - enss * rj(cc, q, 1.0, 1.0 + enss) / 3.0);
	}

	/**
	 * Returns the Jacobian elliptic functions sn(u, kc), cn(u, kc), and dn(u, kc). Here uu = u, while emmc = kc^2.
	 * @param uu uu
	 * @param emmc emmc
	 */
	public static void sncndn(double uu, double emmc)// , float *sn, float *cn,
														// float *dn)
	// Returns the Jacobian elliptic functions sn(u, kc), cn(u, kc), and dn(u,
	// kc). Here uu = u, while
	// emmc = kc2
	{
		double a = 0.0;
		double b = 0.0;
		double c = 0.0;
		double d = 0.0;
		double emc = 0.0;
		double u = 0.0;
		double[] em = new double[14];
		double[] en = new double[14];
		int i = 0;
		int ii = 0;
		int l = 0;
		int bo = 0;

		emc = emmc;
		u = uu;
		if (emc != 0.0)// if (emc)//if (odd!=0) //if (odd)
		{
			// bo=(emc < 0.0);
			if (emc < 0.0)
				bo = 1;
			else
				bo = 0;
			if (bo != 0) // if (bo)
			{
				d = 1.0 - emc;
				emc /= -1.0 / d;
				u *= (d = Math.sqrt(d));
			}
			a = 1.0;
			dn_jacob = 1.0;
			for (i = 1; i <= 13; i++) {
				l = i;
				em[i] = a;
				en[i] = (emc = Math.sqrt(emc));
				c = 0.5 * (a + emc);
				if (Math.abs(a - emc) <= CA * a)
					break;
				emc *= a;
				a = c;
			}
			u *= c;
			sn_jacob = Math.sin(u);
			cn_jacob = Math.cos(u);
			if (sn_jacob != 0.0) // if (*sn)
			{
				a = (cn_jacob) / (sn_jacob);
				c *= a;
				for (ii = l; ii >= 1; ii--) {
					b = em[ii];
					a *= c;
					c *= (dn_jacob);
					dn_jacob = (en[ii] + a) / (b + a);
					a = c / b;
				}
				a = 1.0 / Math.sqrt(c * c + 1.0);
				sn_jacob = (sn_jacob >= 0.0 ? a : -a);
				cn_jacob = c * (sn_jacob);
			}
			if (bo != 0)// if (bo)
			{
				a = (dn_jacob);
				dn_jacob = (cn_jacob);
				cn_jacob = a;
				sn_jacob /= d;
			}
		} else {
			cn_jacob = 1.0 / cosh(u);
			dn_jacob = (cn_jacob);
			sn_jacob = tanh(u);
		}
	}

	/*
	 * //NOT IMPLEMENTED DUE TO The further diff eq implementation
	 * 
	 * #define EPS 1.0e-6 Accuracy parameter. fcomplex aa,bb,cc,z0,dz;
	 * Communicates with hypdrv. int kmax,kount; Used by odeint. float
	 * *xp,**yp,dxsav; fcomplex hypgeo(fcomplex a, fcomplex b, fcomplex c,
	 * fcomplex z) Complex hypergeometric function 2F1 for complex a, b, c, and
	 * z, by direct integration of the hypergeometric equation in the complex
	 * plane. The branch cut is taken to lie along the real axis, Re z > 1. {
	 * void bsstep(float y[], float dydx[], int nv, float *xx, float htry, float
	 * eps, float yscal[], float *hdid, float *hnext, void (*derivs)(float,
	 * float [], float [])); void hypdrv(float s, float yy[], float dyyds[]);
	 * void hypser(fcomplex a, fcomplex b, fcomplex c, fcomplex z, fcomplex
	 * *series, fcomplex *deriv); void odeint(float ystart[], int nvar, float
	 * x1, float x2, float eps, float h1, float hmin, int *nok, int *nbad, void
	 * (*derivs)(float, float [], float []), void (*rkqs)(float [], float [],
	 * int, float *, float, float, float [], float *, float *, void (*)(float,
	 * float [], float []))); int nbad,nok; fcomplex ans,y[3]; float *yy;
	 * kmax=0; if (z.r*z.r+z.i*z.i <= 0.25) { Use series...
	 * hypser(a,b,c,z,&ans,&y[2]); return ans; } else if (z.r < 0.0)
	 * z0=Complex(-0.5,0.0); ...or pick a starting point for the path
	 * integration. else if (z.r <= 1.0) z0=Complex(0.5,0.0); else
	 * z0=Complex(0.0,z.i >= 0.0 ? 0.5 : -0.5); aa=a; Load the global variables
	 * to pass parameters “over the head” of odeint to hypdrv. bb=b; cc=c;
	 * dz=Csub(z,z0); hypser(aa,bb,cc,z0,&y[1],&y[2]); Get starting function and
	 * derivative. yy=vector(1,4); yy[1]=y[1].r; yy[2]=y[1].i; yy[3]=y[2].r;
	 * yy[4]=y[2].i;
	 * odeint(yy,4,0.0,1.0,EPS,0.1,0.0001,&nok,&nbad,hypdrv,bsstep); The
	 * arguments to odeint are the vector of independent variables, its length,
	 * the starting and ending values of the dependent variable, the accuracy
	 * parameter, an initial guess for stepsize, a minimum stepsize, the
	 * (returned) number of good and bad steps taken, and the names of the
	 * derivative routine and the (here Bulirsch-Stoer) stepping routine.
	 * y[1]=Complex(yy[1],yy[2]); free_vector(yy,1,4); return y[1]; }
	 * 
	 * void hypser(fcomplex a, fcomplex b, fcomplex c, fcomplex z, fcomplex
	 * *series, fcomplex *deriv) Returns the hypergeometric series 2F1 and its
	 * derivative, iterating to machine accuracy. For |z| = 1/2 convergence is
	 * quite rapid. { void nrerror(char error_text[]); int n; fcomplex
	 * aa,bb,cc,fac,temp; deriv->r=0.0; deriv->i=0.0; fac=Complex(1.0,0.0);
	 * temp=fac; aa=a; bb=b; cc=c; for (n=1;n<=1000;n++) {
	 * fac=Cmul(fac,Cdiv(Cmul(aa,bb),cc)); deriv->r+=fac.r; deriv->i+=fac.i;
	 * fac=Cmul(fac,RCmul(1.0/n,z));series=Cadd(temp,fac); if (series->r ==
	 * temp.r && series->i == temp.i) return; temp= *series; aa=Cadd(aa,ONE);
	 * bb=Cadd(bb,ONE); cc=Cadd(cc,ONE); }
	 * nrerror("convergence failure in hypser"); }
	 * 
	 * extern fcomplex aa,bb,cc,z0,dz; Defined in hypgeo. void hypdrv(float s,
	 * float yy[], float dyyds[]) Computes derivatives for the hypergeometric
	 * equation, see text equation (5.14.4). { fcomplex z,y[3],dyds[3];
	 * y[1]=Complex(yy[1],yy[2]); y[2]=Complex(yy[3],yy[4]);
	 * z=Cadd(z0,RCmul(s,dz)); dyds[1]=Cmul(y[2],dz);
	 * dyds[2]=Cmul(Csub(Cmul(Cmul(aa,bb),y[1]),Cmul(Csub(cc,
	 * Cmul(Cadd(Cadd(aa,bb),ONE),z)),y[2])), Cdiv(dz,Cmul(z,Csub(ONE,z))));
	 * dyyds[1]=dyds[1].r; dyyds[2]=dyds[1].i; dyyds[3]=dyds[2].r;
	 * dyyds[4]=dyds[2].i; }
	 */
}
