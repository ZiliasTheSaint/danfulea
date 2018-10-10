package danfulea.math;

import java.math.BigInteger;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.text.NumberFormat;
import java.util.Locale;

/**
 * Some math utilities such as base conversions, spline interpolations and number formatting.
 * Some routines are based on EGSnrc MORTRAN code (NRC Canada).
 * @author Dan Fulea, 15 APR. 2011
 *
 */
public class MathUtil {
	public static NumberFormat nf = NumberFormat.getInstance(Locale.US);
	public static String pattern = "0.###E0";
	public static DecimalFormatSymbols dfs = new DecimalFormatSymbols(Locale.US);
	public static DecimalFormat nff = new DecimalFormat(pattern, dfs);
	public static int idigits = 3;

	// ==========Sigma
	public static int IERRsgm = 0;
	public static int $MXDATA = 10;// 10 arrays of date to handle
	public static int $STAT = 200;// maximum data available per array
	// ex Data[n][i] means n siruri de cate i date maxim fiecare
	public static double[][] DATA = new double[$MXDATA][$STAT];

	/**
	 * Default number format with 3 digits.
	 */
	public static void set_defaults() {
		idigits = 3;
		nf.setMinimumFractionDigits(idigits);// 3);//default is 2!!
		nf.setMaximumFractionDigits(idigits);// (3);//default is 2!!
		nf.setGroupingUsed(false);// no 4,568.02 but 4568.02
	}

	/**
	 * Linear interpolation based on EGSnrc fortran/mortran code. 
	 * @param X the desired x-value
	 * @param XA the known array of x-values
	 * @param NX its dimension (length of the array)
	 * @param YA the array of y-values.
	 * @param ISK dummy variable having no use!!
	 * @param XLOG true if we want logarithmic values for x-values
	 * @param YLOG true if we want logarithmic values for y-values 
	 * @return the desired y-value
	 */
	public static double AINTP(double X, double[] XA, int NX, double[] YA,
			int ISK, boolean XLOG, boolean YLOG) {
		double AINTP = 0.;
		boolean XLOGL = XLOG;// "SET LOCAL VARIABLE"
		int I = 0;
		int J = 0;
		boolean ok = false;
		double XI = 0.;
		double XJ = 0.;
		double XV = 0.;
		double YI = 0.;
		double YJ = 0.;
		// " FIND INTERVAL FOR X INTERPOLATION
		for (J = 1; J < NX; J++) {
			if (X < XA[J]) {
				I = J - 1;// looking for the nearest I=min, J=max values
				ok = true;
				break;
			}
		}
		if (!ok) {
			J = NX - 1;
			I = J - 1;
		}

		if (XA[I] <= 0.0) {
			XLOGL = false;
		}
		if (!XLOGL) {
			XI = XA[I];
			XJ = XA[J];
			XV = X;
		} else {
			XI = Math.log(XA[I]);
			XJ = Math.log(XA[J]);
			XV = Math.log(X);
		}
		if ((YLOG) && ((YA[I] == 0.0) || (YA[J] == 0.0))) {
			AINTP = 0.0;
		} else {
			if (YLOG) {
				YI = Math.log(YA[I]);
				YJ = Math.log(YA[J]);
				if (XJ == XI) {
					AINTP = YI;
				} else {
					AINTP = (YI * (XJ - XV) + YJ * (XV - XI)) / (XJ - XI);
				}
				AINTP = Math.exp(AINTP);
			} else {
				YI = YA[I];
				YJ = YA[J];
				if (XJ == XI) {
					AINTP = YI;
				} else {
					AINTP = (YI * (XJ - XV) + YJ * (XV - XI)) / (XJ - XI);
				}
			}
		}
		// System.out.println("AINTP  "+AINTP);
		return AINTP;
	}

	/**
	 * Sets cubic spline interpolation coefficients for the data contained in the array f(n) at the abscissas x(n).
	 * Based on EGSnrc fortran/mortran code. 
	 * @param x the x-array data
	 * @param f the y-array data
	 * @param a the a coefficient array of n-1 dimension to be filled by this routine
	 * @param b the b coefficient array of n-1 dimension to be filled by this routine
	 * @param c the c coefficient array of n-1 dimension to be filled by this routine
	 * @param d the d coefficient array of n-1 dimension to be filled by this routine
	 * @param n the dimension of input arrays (x or y)
	 */
	public static void set_spline(double[] x, double[] f, double[] a,
			double[] b, double[] c, double[] d, int n) {
		// "======================================================================"
		// "
		// " Sets cubic spline interpolation coefficients for the data contained  "
		// " in the array f(n) at the abscissas x(n)                              "
		// " original fortran by I.Kawrakow, NRC,                                 "
		// "======================================================================"
		// ; Copyright NRC;

		int m1 = 0;
		int m2 = 0;
		int m = 0;
		int mr = 0;
		double s = 0.0;
		double r = 0.0;
		m1 = 2;
		m2 = n - 1;
		s = 0.0;
		// DO m=1,m2 [
		for (m = 0; m < m2; m++)// 0 - n-2
		{
			d[m] = x[m + 1] - x[m];
			r = (f[m + 1] - f[m]) / d[m];
			c[m] = r - s;
			s = r;
		}
		// /s,r,c(1),c(n)/=0;
		s = 0.0;
		r = 0.0;
		c[0] = 0.0;// reset to null
		c[n - 1] = 0.0;// reset to null
		// DO m=m1,m2 [
		for (m = m1 - 1; m < m2; m++)// 1 - n-2
		{
			c[m] = c[m] + r * c[m - 1];
			b[m] = 2.0 * (x[m - 1] - x[m + 1]) - r * s;
			s = d[m];
			r = s / b[m];
		}
		// mr = m2;
		mr = m2 - 1;// n-2
		// DO m=m1,m2 [
		for (m = m1 - 1; m < m2; m++)// 1 - n-2
		{
			c[mr] = (d[mr] * c[mr + 1] - c[mr]) / b[mr];
			mr = mr - 1;
		}
		// DO m=1,m2 [
		for (m = 0; m < m2; m++)// 0 - n-2
		{
			s = d[m];
			r = c[m + 1] - c[m];
			d[m] = r / s;
			c[m] = 3.0 * c[m];
			b[m] = (f[m + 1] - f[m]) / s - (c[m] + r) * s;
			a[m] = f[m];
		}
	}

	/**
	 * Returns the value of the function at s using the spline coefficients a,b,c,d, which must have been set using set_spline 
	 * @param s the value at which function is evaluated 
	 * @param x the x-array of dimension n
	 * @param a the a coefficient array of n-1 dimension
	 * @param b the b coefficient array of n-1 dimension
	 * @param c the c coefficient array of n-1 dimension
	 * @param d the d coefficient array of n-1 dimension
	 * @param n the dimension of x-array
	 * @return the value of the function at s
	 */
	public static double spline(double s, double[] x, double[] a, double[] b,
			double[] c, double[] d, int n) {
		// "======================================================================"
		// "                                                                      "
		// " Returns the value of the function at s using the spline coefficients "
		// " a,b,c,d, which must have been set using set_spline                   "
		// "                                                                      "
		// " original fortran I.Kawrakow, NRC                                     "
		// "======================================================================"

		// ; Copyright NRC;

		// $INTEGER n;
		// $REAL s,x(n),a(n),b(n),c(n),d(n);

		int m_lower = 0;
		int m_upper = 0;
		int direction = 0;
		int m = 0;
		int ml = 0;
		int mu = 0;
		int mav = 0;
		double q = 0.0;

		// if( x(1) > x(n) ) [ direction = 1; m_lower = n; m_upper = 0; ]
		if (x[0] > x[n - 1]) {
			direction = 1;
			m_lower = n;
			m_upper = 0;
		} else {
			direction = 0;
			m_lower = 0;
			m_upper = n;
		}
		// if ( s >= x(m_upper + direction) ) [
		if (s >= x[m_upper + direction - 1]) {
			m = m_upper + 2 * direction - 1;
		}
		// else if ( s <= x(m_lower+1-direction) ) [
		else if (s <= x[m_lower - direction]) {
			m = m_lower - 2 * direction + 1;
		} else {// " Perform a binary search to find the interval s is in "
			ml = m_lower;
			mu = m_upper;
			// while ( iabs(mu-ml) > 1 ) [
			while (Math.abs(mu - ml) > 1) {
				mav = (ml + mu) / 2;
				// if( s < x(mav) )
				if (s < x[mav - 1]) {
					mu = mav;
				} else {
					ml = mav;
				}
			}
			m = mu + direction - 1;
		}

		q = s - x[m - 1];
		// double spline = a(m) + q*(b(m) + q*(c(m) + q*d(m)));
		double spline = a[m - 1] + q
				* (b[m - 1] + q * (c[m - 1] + q * d[m - 1]));
		return spline;
	}

	/**
	 * Given the lower and upper limit of integration, x1 and x2, and given n, this routine returns arrays x and w,
	 * containing the abscissas and weights of the Gauss-Legendre n - point quadrature formula
	 * @param x1 lower limit of integration
	 * @param x2 uper limit of integration
	 * @param x abscissas to be filled by this routine
	 * @param w weights to be filled by this routine
	 * @param n array dimension
	 */
	public static void gauss_legendre(double x1, double x2, double[] x,
			double[] w, int n) {

		// " Given the lower and upper limit of integration, x1 and x2,
		// " and given n, this routine returns arrays x and w,
		// " containing the abscissas and weights of the Gauss-Legendre
		// " n - point quadrature formula
		// "******************************************************************************
		// ; Copyright NRC;

		double eps = 3.E-14;
		double Pi = Math.PI;
		// parameter (eps = 3.D-14, Pi = 3.141592654D0);

		int i = 0;
		int m = 0;
		int j = 0;
		double xm = 0.0;
		double xl = 0.0;
		double z = 0.0;
		double z1 = 0.0;
		double p1 = 0.0;
		double p2 = 0.0;
		double p3 = 0.0;
		double pp = 0.0;

		m = (n + 1) / 2;
		// xm=0.5d0*(x2+x1); xl=0.5d0*(x2-x1);
		xm = 0.5 * (x2 + x1);
		xl = 0.5 * (x2 - x1);
		for (i = 1; i <= m; i++) {
			// z=cos(Pi*(i-.25d0)/(n+.5d0));
			z = Math.cos(Pi * (i - 0.25) / (n + 0.5));
			while (Math.abs(z - z1) >= eps)
			// LOOP
			{
				// p1=1.d0;
				// p2=0.d0;
				p1 = 1.0;
				p2 = 0.0;
				for (j = 1; j <= n; j++) {
					p3 = p2;
					p2 = p1;
					// p1=((2.d0*j-1.d0)*z*p2-(j-1.d0)*p3)/j;
					p1 = ((2.0 * j - 1.0) * z * p2 - (j - 1.0) * p3) / j;
				}
				// pp=n*(z*p1-p2)/(z*z-1.d0);
				pp = n * (z * p1 - p2) / (z * z - 1.0);
				z1 = z;
				z = z1 - p1 / pp;
			}// UNTIL (abs(z-z1) < eps);
				// x(i)=xm-xl*z; x(n+1-i)=xm+xl*z;
				// w(i)=2.d0*xl/((1.d0-z*z)*pp*pp); w(n+1-i)=w(i);
			x[i - 1] = xm - xl * z;
			x[n - i] = xm + xl * z;
			// w(i-1)=2.d0*xl/((1.d0-z*z)*pp*pp);
			w[i - 1] = 2.0 * xl / ((1.0 - z * z) * pp * pp);
			w[n - i] = w[i];
		}
	}

	// "************************************************************************
	// " an error function routine which is needed since some of
	// " the compiler don't have it as an intrinsic
	// " Originally came from some library somewhere (Harwell I think)
	// " recoded in mortran
	// "************************************************************************
	/**
	 * Error function routine
	 * @param X value where ERF is evaluated
	 * @return the error function value
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
		//int N = n;
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
	 * Complementary error function evaluated at x
	 * @param x x
	 * @return the erfc at x
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

	// --------------convertors to be independent with other class
	// file--------------------------
	/**
	 * Converts an int to another base
	 * @param value the string representation of an int (base 10 that is)
	 * @param base the desired base for conversion
	 * @return the String representation of the number in desired base
	 * @throws NumberFormatException can throw this exception
	 */
	public static String convertIntToOtherBase(String value, int base)
			throws NumberFormatException// 3a3f15f7=format Hexazecimal ex.
	{
		BigInteger bi = new BigInteger(value);
		return bi.toString(base);
	}

	/**
	 * Converts a number from other base to base 10
	 * @param value the String representation of a number
	 * @param base the base of the number
	 * @return the String representation of number in base 10
	 * @throws NumberFormatException can throw this exception
	 */
	public static String convertOtherBaseToInt(String value, int base)
			throws NumberFormatException {
		BigInteger bi = new BigInteger(value, base);
		return bi.toString(10);
	}

	/**
	 * Converts a string to int value.
	 * 
	 * @param value
	 *            the string
	 * @return the int value of input string
	 * @throws NumberFormatException
	 * can throws this exception
	 */
	public static int stringToInt(String value) throws NumberFormatException {
		value = value.trim();
		if (value.length() == 0) {
			return 0;
		} else {
			return Integer.parseInt(value);
		}
	}

	/**
	 * Converts an int number to string.
	 * 
	 * @param i
	 *            the int value
	 * @return the string representation
	 */
	public static String intToString(int i) {
		return Integer.toString(i);
	}

	/**
	 * Converts a string to float value.
	 * 
	 * @param value
	 *            the string
	 * @return the float value of input string
	 * @throws NumberFormatException
	 * can throws this exception
	 */
	public static float stringToFloat(String value)
			throws NumberFormatException {
		value = value.trim();
		if (value.length() == 0) {
			return 0;
		} else {
			return Float.parseFloat(value);
		}
	}

	/**
	 * Converts a float number to string.
	 * 
	 * @param f
	 *            the float value
	 * @return the string representation
	 */
	public static String floatToString(float f) {
		return Float.toString(f);
	}

	/**
	 * Converts a string to double value.
	 * 
	 * @param value
	 *            the string
	 * @return the double value of input string
	 * @throws NumberFormatException
	 * can throws this exception
	 */
	public static double stringToDouble(String value)
			throws NumberFormatException {
		value = value.trim();
		if (value.length() == 0) {
			return 0;
		} else {
			return Double.parseDouble(value);
		}
	}

	/**
	 * Converts a double number to string.
	 * 
	 * @param d
	 *            the double value
	 * @return the string representation
	 */
	public static String doubleToString(double d) {
		return Double.toString(d);
	}

	/**
	 * Format a number d using either number format (true) or scientific format(false) and if its length is less then an offset
	 * value then remaining characters are filled with blank. 
	 * @param d the number as a double
	 * @param offset the offset
	 * @param nrformat true if number format false if scientific
	 * @return the String representation of the formatted number
	 */
	public static String format(double d, int offset, boolean nrformat) {
		String result = "";
		String f = "";
		if (nrformat) {
			f = nf.format(d);
		} else
			f = nff.format(d);
		int k = f.length();
		if (k <= offset) {
			int i = offset - k;
			char[] chs = new char[i];
			for (int j = 0; j < i; j++)
				chs[j] = ' ';// blank
			// f=f+new String(chs);//after
			f = new String(chs) + f;// before
		} else {// 1 character blank by default
				// f=f+" ";//after
			f = " " + f;// before
		}

		result = f;
		return result;
	}

	/**
	 * Format a number d and if its length is less then an offset
	 * value then remaining characters are filled with blank. 
	 * @param d the number as an int
	 * @param offset the offset
	 * @return the String representation of the formatted number
	 */
	public static String format(int d, int offset) {
		String result = "";
		String f = "" + d;
		int k = f.length();
		if (k <= offset) {
			int i = offset - k;
			char[] chs = new char[i];
			for (int j = 0; j < i; j++)
				chs[j] = ' ';// blank
			// f=f+new String(chs);//after
			f = new String(chs) + f;// before
		} else {// 1 character blank by default
				// f=f+" ";//after
			f = " " + f;// before
		}

		result = f;
		return result;
	}

	/**
	 * Try to auto-format the number using either number format or scientific format based on offset.
	 * @param d the number as a double
	 * @param offset the offset
	 * @return the String representation of the formatted number
	 */
	public static String format(double d, int offset) {
		int k = 0;
		String result = "";
		String f = "";
		f = nf.format(d);

		double dbl = stringToDouble(f);
		if (dbl == 0.0) {
			f = nff.format(d);
			k = f.length();
			if (k <= offset) {
				int i = offset - k;
				char[] chs = new char[i];
				for (int j = 0; j < i; j++)
					chs[j] = ' ';// blank
				// f=f+new String(chs);//after
				f = new String(chs) + f;// before
				return f;
			} else {
				// f=f+" ";//after
				f = " " + f;// before
				return f;
			}
		}

		k = f.length();
		if (k <= offset) {
			int i = offset - k;
			char[] chs = new char[i];
			for (int j = 0; j < i; j++)
				chs[j] = ' ';// blank
			// f=f+new String(chs);//after
			f = new String(chs) + f;// before
			return f;
		} else {
			f = nff.format(d);
			k = f.length();
			if (k <= offset) {
				int i = offset - k;
				char[] chs = new char[i];
				for (int j = 0; j < i; j++)
					chs[j] = ' ';// blank
				// f=f+new String(chs);//after
				f = new String(chs) + f;// before
				return f;
			} else {
				// f=f+" ";//after
				f = " " + f;// before
			}
		}

		result = f;
		return result;
	}

	/**
	 * Special case for number format option which is 1 fraction digit.
	 * @param d the number as a double
	 * @param offset the offset
	 * @param nrformat the format, true if number format false if scientific
	 * @return the String representation of formatted number
	 */
	public static String format1(double d, int offset, boolean nrformat) {
		String result = "";
		String f = "";
		NumberFormat nf1 = NumberFormat.getInstance(Locale.US);
		nf1.setMinimumFractionDigits(1);// 3);//default is 2!!
		nf1.setMaximumFractionDigits(1);// (3);//default is 2!!
		nf1.setGroupingUsed(false);// no 4,568.02 but 4568.02

		if (nrformat) {
			f = nf1.format(d);
		} else
			f = nff.format(d);
		int k = f.length();
		if (k <= offset) {
			int i = offset - k;
			char[] chs = new char[i];
			for (int j = 0; j < i; j++)
				chs[j] = ' ';// blank
			// f=f+new String(chs);//after
			f = new String(chs) + f;// before
		}
		// else
		// {//1 character blank by default
		// f=f+" ";//after
		// f=" "+f;//before
		// }

		result = f;
		return result;
	}

	/**
	 * Format a number d and if its length is less then an offset
	 * value then remaining characters are filled with blank. 
	 * @param d the number as a String
	 * @param offset the offset
	 * @return the String representation of the formatted number
	 */
	public static String format(String d, int offset) {
		String result = "";
		String f = "" + d;
		int k = f.length();
		if (k <= offset) {
			int i = offset - k;
			char[] chs = new char[i];
			for (int j = 0; j < i; j++)
				chs[j] = ' ';// blank
			// f=f+new String(chs);//after
			f = new String(chs) + f;// before

		} else {
			// f=f+" ";//after
			f = " " + f;// before

		}

		result = f;
		return result;
	}

	// "*******************************************************************************
	// "
	// "
	// " *****************
	// " * *
	// " * SIGMA.MORTRAN *
	// " * *
	// " *****************
	// "
	// "
	// " SIGMA IS A STATISTICAL ANALYSIS ROUTINE DESIGNED TO BE USED BY EGS
	// " USER PROGRAMS TO GIVE THE TOTALS OR AVERAGES AND THEIR UNCERTAINTIES
	// " OF THE DATA CALCULATED BY THE MONTE CARLO CODE.
	// " THE UNCERTAINTIES ARE RETURNED AS PERCENTS.
	// "
	// " VARIABLES
	// " =========
	// "
	// " DATA(NDATA,ISTAT) THE TWO DIMENSIONAL ARRAY OF DATA TO BE
	// " ANALYZED. ISTAT IS THE NUMBER OF STATISTICAL
	// " BATCHES AND NDATA IS THE NUMBER OF ERRORS TO
	// " BE CALCULATED. AFTER THE END OF THE CALCULATION,
	// " DATA(N,1) CONTAINS THE TOTAL OR AVERAGE AND
	// " DATA(N,2) CONTAINS THE ERROR. NDATA SHOULD
	// " BE < OR = $MAXDATA AND ISTAT SHOULD BE < OR =
	// " $STAT WHCH MUST BE DEFINED IN THE MAIN ROUTINE.
	// " Note $STAT must be 2 or greater, even if istat=1
	// "
	// " MODE = 0 => ANALYSIS ON MEAN VALUES WHERE ZERO DATA IS
	// " IGNORED. (eg. STOPPING POWER RATIO)
	// " = 1 => ANALYSIS ON MEAN VALUES WHERE ZERO DATA IS NOT
	// " IGNORED. (e.g. DOSE)
	// " = 2 => ANALYSIS ON TOTAL VALUES (eg. TOTAL EDEP)
	// "
	// " IERR = 0 => NORMAL COMPLETION.
	// " = 1 => WARNING: MODE OUT OF RANGE, DEFAULTED TO 0
	// " = 10 => ERROR: ONLY ONE BATCH INPUT, QUICK CALCULATION
	// " DONE. ERROR=99.9%
	// " = 11 => ERROR: NO NON-ZERO DATA FOUND IN A GIVEN SET,
	// " ERROR=99.9%
	// " = -1 => FATAL ERROR: NDATA OR ISTAT OUT OF RANGE, NO
	// " CALCULATION DONE.
	// "
	// "
	// " VERSION 1 A.F.B. 83/7/22
	// " Version 2 IK Jan 6 6000 implemented implicit none
	// "
	// "*******************************************************************************

	/**
	 * Statistical routine based on EGSnrc MORTRAN code. Just for testing.
	 * @param NDATA NDATA
	 * @param ISTAT ISTAT
	 * @param MODE MODE
	 */
	public static void SIGMA(int NDATA, int ISTAT, int MODE)// ,int IERR)
	{

		// ;Copyright NRC;

		// $INTEGER NDATA,ISTAT,MODE,IERR;

		// REPLACE {;COMIN/ERROR/;} WITH {
		// ;COMMON/ERROR/DATA($MXDATA,$STAT);
		// $REAL data;
		// }
		// double[][] DATA=new DATA[$MXDATA][$STAT];
		// ;COMIN/ERROR/;

		int N = 0;
		int NON0 = 0;
		int I = 0;
		double STAT = 0.0;
		double SDENOM = 0.0;
		//double emax = 0.0;
		double AVG = 0.0;
		double ERROR = 0.0;
		double DATUM = 0.0;
		double ARGMNT = 0.0;
		// "It is a good idea to use double precision"
		// "in cases with very low stat. uncertainties"

		double EMAX = 99.9;

		IERRsgm = 0; // "ASSUME NORMAL COMPLETION"
		boolean transferb = false;

		// "TEST INPUTS AND SET ERROR CODES AND RETURN IF NEEDED."

		if ((MODE < 0) || (MODE > 2)) {
			MODE = 2;
			IERRsgm = 1;
		}

		if (((NDATA <= 0) || (NDATA > $MXDATA))
				|| ((ISTAT <= 0) || (ISTAT > $STAT))) {
			IERRsgm = -1;
			return;// "FATAL INPUT ERROR, RETURN IMMEDIATELY"
		}
		if (ISTAT == 1) {
			IERRsgm = 10;// "ONLY ONE STATISTICAL BATCH, QUICK CALCULATION"
			// DO N=1,NDATA[DATA(N,2)=EMAX;]
			for (N = 1; N <= NDATA; N++) {
				// DATA(N,2)=EMAX;
				DATA[N - 1][1] = EMAX;
			}
			return;
		}

		// "MOST ANOMALIES HAVE BEEN HANDLED. NOW DO THE ANALYSIS"

		if (MODE != 0) {
			STAT = ISTAT;
			SDENOM = STAT * (STAT - 1.);
		}
		for (N = 1; N <= NDATA; N++) {
			transferb = false;
			NON0 = 0;// "NON-ZERO COUNTER"
			AVG = 0.0;
			ERROR = 0.0;
			for (I = 1; I <= ISTAT; I++) {
				DATUM = DATA[N - 1][I - 1];// DATA(N,I);
				if (DATUM != 0.0) {
					NON0 = NON0 + 1;
					AVG = AVG + DATUM;
					ERROR = ERROR + DATUM * DATUM;
				}
			}
			if (NON0 == 0) {
				IERRsgm = 11;
				ERROR = EMAX;
				// GOTO :TRANSFER:;//"NO NON-ZERO DATA "
				transferb = true;

			} else if ((NON0 == 1) && (MODE == 0)) {
				ERROR = EMAX;
				// GOTO:TRANSFER:;//"ONLY ONE DATUM"
				transferb = true;
			} else {
				if (MODE == 0) {
					STAT = NON0;
					SDENOM = STAT * (STAT - 1.);
				}
			}

			if (!transferb) {
				AVG = AVG / STAT;
				ARGMNT = ERROR - STAT * AVG * AVG;
				// "FLAG -VE SQUARE ROOTS THAT CAN ONLY OCCUR DUE TO ROUND-OFF ERRORS"
				if (ARGMNT < 0.0) {
					// OUTPUT ARGMNT,ERROR,STAT,AVG,SDENOM;
					// (' ***** - SQ RT IN SIGMA.
					// ARGMNT,ERROR,STAT,AVG,SDENOM='/' ',5E12.4);
					ARGMNT = 0.0;
				}
				ERROR = Math.sqrt(ARGMNT / SDENOM);

				if (AVG == 0.) {
					ERROR = EMAX;
				} else {
					ERROR = 100. * ERROR / Math.abs(AVG);
				}

				if (MODE == 2)
					AVG = AVG * STAT;
			}

			// :TRANSFER:;
			// DATA(N,1)=AVG;DATA(N,2)=MIN(EMAX,ERROR);
			DATA[N - 1][0] = AVG;
			DATA[N - 1][1] = Math.min(EMAX, ERROR);
		}// "END OF NDATA LOOP"
		return;
	}// "END OF SIGMA"

	// "============================================================================="
	/**
	 * Heap sort will sort the real array rarray of dimension n in ascending order and at the same time put into the integer array
	 * jarray the original position of the elements, e.g. if rarray was on input (5,14,8,2), it will be after completion
	 * of heap_sort (2,5,8,14) and jarray will be (4,1,3,2). heap_sort uses the heap sort algorithm, the implementation is
	 * based on hpsort from Numerical Recipies with a couple of modifications.
	 * @param n dimension of the array
	 * @param rarray input array which will be sorted
	 * @param jarray original indexes of input array
	 */
	public static void heap_sort(int n, double[] rarray, int[] jarray) {
		// "************************************************************************
		// " egs_heap_sort will sort the real array rarray of dimension n in
		// " ascending order and at the same time put into the integer array
		// " jarray the original position of the elements, e.g.
		// " if rarray was on input (5,14,8,2), it will be after completion
		// " of heap_sort (2,5,8,14) and jarray will be (4,1,3,2).
		// " heap_sort uses the heap sort algorithm, the implementation is
		// " based on hpsort from Numerical Recipies with a couple of
		// " modifications.
		// "
		// " Iwan Kawrakow, NRC, July 2001
		// "*************************************************************************

		// implicit none;

		// $INTEGER n,jarray(*);
		// $REAL rarray(*);

		int i = 0;
		int ir = 0;
		int j = 0;
		int l = 0;
		int ira = 0;
		double rra = 0.0;

		for (i = 1; i <= n; i++) {
			jarray[i - 1] = i;
		}// { jarray(i)=i; }
		if (n < 2)
			return;
		l = n / 2 + 1;
		ir = n;

		while (true) {
			if (l > 1) {
				l = l - 1;
				rra = rarray[l - 1];// rarray(l);
				ira = l;
			} else {
				rra = rarray[ir - 1];// rarray(ir);
				ira = jarray[ir - 1];// jarray(ir);
				rarray[ir - 1] = rarray[0];// rarray(ir)=rarray(1);
				jarray[ir - 1] = jarray[0];// jarray(ir)=jarray(1);
				ir = ir - 1;
				if (ir == 1) {
					rarray[0] = rra;// rarray(1)=rra;
					jarray[0] = ira; // jarray(1)=ira;
					return;
				}
			}
			i = l;
			j = l + l;
			while (true) {
				if (j > ir)
					break;
				if (j < ir) {
					// IF (rarray(j) < rarray(j+1) ) j=j+1;
					if (rarray[j - 1] < rarray[j])
						j = j + 1;
				}
				if (rra < rarray[j - 1])// if (rra < rarray(j))
				{
					// rarray(i)=rarray(j); jarray(i)=jarray(j);
					rarray[i - 1] = rarray[j - 1];
					jarray[i - 1] = jarray[j - 1];
					i = j;
					j = j + j;
				} else {
					j = ir + 1;
				}
			}
			// rarray(i)=rra; jarray(i)=ira;
			rarray[i - 1] = rra;
			jarray[i - 1] = ira;
		}
		// return; //end;
	}

}
