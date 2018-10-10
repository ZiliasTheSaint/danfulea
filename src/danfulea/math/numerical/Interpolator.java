package danfulea.math.numerical;

/**
 * Interpolator for function evaluation. 
 * Based on Numerical Recipes in C (Cambridge Univ.).
 * 
 * @author Dan Fulea, 29 SEP. 2006
 */
public class Interpolator {
	public static double ypoli = 0.0;
	public static double dypoli = 0.0;
	public static boolean failB = false;
	public static String failS = "";

	public static double TINY = 1.0e-25;

	// @@@@@@Neville algorithm much better than the old Lagrange 
	/**
	 * Given arrays xa[1..n] and ya[1..n], and given a value x, this routine returns a value ypoli, and 
	 * an error estimate dypoli. If P(x) is the polynomial of degree N - 1 such that P(xai) = yai, i = 
	 * 1, . . . , n, then the returned value y = P(x). The interpolating polynomial of degree N - 1 through 
	 * the N points y1 = f(x1), y2 = f(x2), . . . , yN = f(xN) is given explicitly by 
	 * Lagrange’s classical formula. 
	 * @param xa xa
	 * @param ya ya
	 * @param n n
	 * @param x x
	 */
	public static void polint(double[] xa, double[] ya, int n, double x)// ,
																		// double
																		// ypoli,
																		// double
																		// dypoli)
	// Given arrays xa[1..n] and ya[1..n], and given a value x, this routine
	// returns a value y, and
	// an error estimate dy. If P(x) is the polynomial of degree N - 1 such that
	// P(xai) = yai, i =
	// 1, . . . , n, then the returned value y = P(x).
	// The interpolating polynomial of degree N - 1 through
	// the N points y1 = f(x1), y2 = f(x2), . . . , yN = f(xN) is given
	// explicitly by
	// Lagrange’s classical formula,=>etc.
	// It is not terribly wrong to implement the Lagrange formula
	// straightforwardly,
	// but it is not terribly right either. The resulting algorithm gives no
	// error estimate, and
	// it is also somewhat awkward to program. A much better algorithm (for
	// constructing
	// the same, unique, interpolating polynomial) is Neville’s algorithm,
	{
		int i = 0;
		int m = 0;
		int ns = 1;
		double den = 0.0;
		double dif = 0.0;
		double dift = 0.0;
		double ho = 0.0;
		double hp = 0.0;
		double w = 0.0;
		double[] c = new double[n];
		double[] d = new double[n];// float *c,*d;

		failB = false;
		dif = Math.abs(x - xa[0]);// fabs(x-xa[1]);
		// c=vector(1,n);
		// d=vector(1,n);
		for (i = 1; i <= n; i++) { // Here we find the index ns of the closest
									// table entry,
			if ((dift = Math.abs(x - xa[i - 1])) < dif) // if (
														// (dift=fabs(x-xa[i]))
														// < dif)
			{
				ns = i;
				dif = dift;
			}
			c[i - 1] = ya[i - 1];// c[i]=ya[i]; and initialize the tableau of
									// c’s and d’s.
			d[i - 1] = ya[i - 1];// d[i]=ya[i];
		}
		ns = ns - 1;// @@@@@@@@@@@@@@for array matching!!!
		ypoli = ya[ns--];// *y=ya[ns--]; This is the initial approximation to y.
		for (m = 1; m < n; m++) { // For each column of the tableau,we loop over
									// the current c’s and d’s and update them.
			for (i = 1; i <= n - m; i++) {
				ho = xa[i - 1] - x;// ho=xa[i]-x;
				hp = xa[i + m - 1] - x;// hp=xa[i+m]-x;
				w = c[i] - d[i - 1];// w=c[i+1]-d[i];
				if ((den = ho - hp) == 0.0) {
					// nrerror("Error in routine polint");
					failS = "Error in routine polint";// System.out.println("Error in routine polint");
					failB = true;
					return;
				}
				// This error can occur only if two input xa’s are (to within
				// roundoff) identical.
				den = w / den;
				d[i - 1] = hp * den;// d[i]=hp*den; Here the c’s and d’s are
									// updated.
				c[i - 1] = ho * den;// c[i]=ho*den;
			}
			// *y += (*dy=(2*ns < (n-m) ? c[ns+1] : d[ns--]));
			ypoli += (dypoli = (2 * (ns + 1) < (n - m) ? c[ns + 1] : d[ns--]));
			// After each column in the tableau is completed, we decide which
			// correction, c or d,
			// we want to add to our accumulating value of y, i.e., which path
			// to take through the
			// tableau—forking up or down. We do this in such a way as to take
			// the most “straight
			// line” route through the tableau to its apex, updating ns
			// accordingly to keep track of
			// where we are. This route keeps the partial approximations
			// centered (insofar as possible)
			// on the target x. The last dy added is thus the error indication.
		}
		// free_vector(d,1,n);
		// free_vector(c,1,n);

	}

	// @@@@@@Rational Function Interpolation
	/*
	 * Some functions are not well approximated by polynomials, but are well
	 * approximated by rational functions, that is quotients of polynomials.
	 */
	public static double yrat = 0.0;
	public static double dyrat = 0.0;

	/**
	 * Given arrays xa[1..n] and ya[1..n], and given a value of x, this routine returns a value of 
	 * yrat and an accuracy estimate dyrat. The value returned is that of the diagonal rational function, 
	 * evaluated at x, which passes through the n points (xai, yai), i = 1...n.
	 * @param xa xa
	 * @param ya ya
	 * @param n n
	 * @param x x
	 */
	public static void ratint(double[] xa, double[] ya, int n, double x)// ,
																		// float
																		// *y,
																		// float
																		// *dy)
	// Given arrays xa[1..n] and ya[1..n], and given a value of x, this routine
	// returns a value of
	// y and an accuracy estimate dy. The value returned is that of the diagonal
	// rational function,
	// evaluated at x, which passes through the n points (xai, yai), i = 1...n.
	{
		int m = 0;
		int i = 0;
		int ns = 1;
		double w = 0.0;
		double t = 0.0;
		double hh = 0.0;
		double h = 0.0;
		double dd = 0.0;
		double[] c = new double[n];
		double[] d = new double[n];

		// c=vector(1,n);
		// d=vector(1,n);
		failB = false;
		hh = Math.abs(x - xa[0]);// fabs(x-xa[1]);
		for (i = 1; i <= n; i++) {
			h = Math.abs(x - xa[i - 1]);// h=fabs(x-xa[i]);
			if (h == 0.0) {
				yrat = ya[i - 1];// ya[i];
				dyrat = 0.0;
				return;// FREERETURN
			} else if (h < hh) {
				ns = i;
				hh = h;
			}
			c[i - 1] = ya[i - 1];// c[i]=ya[i];
			d[i - 1] = ya[i - 1] + TINY;// d[i]=ya[i]+TINY; //The TINY part is
										// needed to prevent a rare
										// zero-over-zero
										// condition.
		}
		ns = ns - 1;// @@@@@@@@@@@@@@for array matching!!!
		yrat = ya[ns--];
		for (m = 1; m < n; m++) {
			for (i = 1; i <= n - m; i++) {
				w = c[i] - d[i - 1];// w=c[i+1]-d[i];
				h = xa[i + m - 1] - x;// h=xa[i+m]-x;// h will never be zero,
										// since this was tested in the
										// initializing
										// loop.
				t = (xa[i - 1] - x) * d[i - 1] / h;// t=(xa[i]-x)*d[i]/h;
				dd = t - c[i];// dd=t-c[i+1];
				if (dd == 0.0) {
					// nrerror("Error in routine ratint");
					failS = "Error in routine ratint";
					failB = true;
					return;
				}
				// This error condition indicates that the interpolating
				// function has a pole at the
				// requested value of x.
				dd = w / dd;
				d[i - 1] = c[i] * dd;// d[i]=c[i+1]*dd;
				c[i - 1] = t * dd;// c[i]=t*dd;
			}
			// *y += (*dy=(2*ns < (n-m) ? c[ns+1] : d[ns--]));
			yrat += (dyrat = (2 * (ns + 1) < (n - m) ? c[ns + 1] : d[ns--]));
		}
		return;// FREERETURN
	}

	// @@@@@@@@@SPLINE
	/**
	 * Given arrays x[1..n] and y[1..n] containing a tabulated function, i.e., yi = f(xi), with 
	 * x1 ... xN in ascending order, and given values yp1 and ypn for the first derivative of the interpolating 
	 * function at points 1 and n, respectively, this routine returns an array y2[1..n] that contains 
	 * the second derivatives of the interpolating function at the tabulated points xi. If yp1 and/or 
	 * ypn are equal to 1 × 1030 or larger, the routine is signaled to set the corresponding boundary 
	 * condition for a natural spline, with zero second derivative on that boundary.
	 
	 * @param x x
	 * @param y y
	 * @param n n
	 * @param yp1 yp1
	 * @param ypn ypn
	 * @param y2 y2
	 */
	public static void spline(double[] x, double[] y, int n, double yp1,
			double ypn, double[] y2)
	// Given arrays x[1..n] and y[1..n] containing a tabulated function, i.e.,
	// yi = f(xi), with
	// x1 < x2 < .. . < xN, and given values yp1 and ypn for the first
	// derivative of the interpolating
	// function at points 1 and n, respectively, this routine returns an array
	// y2[1..n] that contains
	// the second derivatives of the interpolating function at the tabulated
	// points xi. If yp1 and/or
	// ypn are equal to 1 × 1030 or larger, the routine is signaled to set the
	// corresponding boundary
	// condition for a natural spline, with zero second derivative on that
	// boundary.
	// double
	// yp1=(xi[1]-xi[0])/(yi[1]-yi[0])-((xi[2]-xi[1])/(yi[2]-yi[1]))+(xi[2]-xi[0])/(yi[2]-yi[0]);
	// double
	// ypn=-((xi[ni-2]-xi[ni-3])/(yi[ni-2]-yi[ni-3]))+(xi[ni-1]-xi[ni-2])/(yi[ni-1]-yi[ni-2])+(xi[ni-1]-xi[ni-3])/(yi[ni-1]-yi[ni-3]);

	{
		int i = 0;
		int k = 0;
		double p = 0.0;
		double qn = 0.0;
		double sig = 0.0;
		double un = 0.0;
		double[] u = new double[n - 1];
		// u=vector(1,n-1);
		if (yp1 > 0.99e30) // The lower boundary condition is set either to be
							// “natural” or else to have a specified first
							// derivative.
			y2[0] = u[0] = 0.0;// y2[1]=u[1]=0.0;
		else {
			y2[0] = -0.5;// y2[1] = -0.5;
			u[0] = (3.0 / (x[1] - x[0]))
					* ((y[1] - y[0]) / (x[1] - x[0]) - yp1);// u[1]=(3.0/(x[2]-x[1]))*((y[2]-y[1])/(x[2]-x[1])-yp1);
		}
		for (i = 2; i <= n - 1; i++) {
			// This is the decomposition loop of the tridiagonal algorithm.
			// y2 and u are used for temporary
			// storage of the decomposed factors.
			sig = (x[i - 1] - x[i - 2]) / (x[i] - x[i - 2]);// sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);
			p = sig * y2[i - 2] + 2.0;// p=sig*y2[i-1]+2.0;
			y2[i - 1] = (sig - 1.0) / p;// y2[i]=(sig-1.0)/p;
			u[i - 1] = (y[i] - y[i - 1]) / (x[i] - x[i - 1])
					- (y[i - 1] - y[i - 2]) / (x[i - 1] - x[i - 2]);// u[i]=(y[i+1]-y[i])/(x[i+1]-x[i])
																	// -
																	// (y[i]-y[i-1])/(x[i]-x[i-1]);
			u[i - 1] = (6.0 * u[i - 1] / (x[i] - x[i - 2]) - sig * u[i - 2])
					/ p;// u[i]=(6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
		}
		if (ypn > 0.99e30) // the upper boundary condition is set either to be
							// “natural”
			qn = un = 0.0;
		else { // or else to have a specified first derivative.
			qn = 0.5;
			// un=(3.0/(x[n]-x[n-1]))*(ypn-(y[n]-y[n-1])/(x[n]-x[n-1]));
			un = (3.0 / (x[n - 1] - x[n - 2]))
					* (ypn - (y[n - 1] - y[n - 2]) / (x[n - 1] - x[n - 2]));
		}
		y2[n - 1] = (un - qn * u[n - 2]) / (qn * y2[n - 2] + 1.0);// y2[n]=(un-qn*u[n-1])/(qn*y2[n-1]+1.0);
		for (k = n - 1; k >= 1; k--)
			// This is the backsubstitution loop of the tridiagonal algorithm.
			y2[k - 1] = y2[k - 1] * y2[k] + u[k - 1];// y2[k]=y2[k]*y2[k+1]+u[k];
		// free_vector(u,1,n-1);
	}

	/**
	 * Call the spline function using pre-computed yp1 and ypn.
	 * * <p>
	 * yp1=(xi[1]-xi[0])/(yi[1]-yi[0])-((xi[2]-xi[1])/(yi[2]-yi[1]))+(xi[2]-xi[0])/(yi[2]-yi[0]);<p>
	 * ypn=-((xi[ni-2]-xi[ni-3])/(yi[ni-2]-yi[ni-3]))+(xi[ni-1]-xi[ni-2])/(yi[ni-1]-yi[ni-2])+(xi[ni-1]-xi[ni-3])/(yi[ni-1]-yi[ni-3]);
	 * @param xi xi
	 * @param yi yi
	 * @param ni ni
	 * @param y2 y2
	 */
	public static void spline_default(double[] xi, double[] yi, int ni,
			double[] y2) {
		double yyp1 = (xi[1] - xi[0]) / (yi[1] - yi[0])
				- ((xi[2] - xi[1]) / (yi[2] - yi[1])) + (xi[2] - xi[0])
				/ (yi[2] - yi[0]);
		double yypn = -((xi[ni - 2] - xi[ni - 3]) / (yi[ni - 2] - yi[ni - 3]))
				+ (xi[ni - 1] - xi[ni - 2]) / (yi[ni - 1] - yi[ni - 2])
				+ (xi[ni - 1] - xi[ni - 3]) / (yi[ni - 1] - yi[ni - 3]);

		spline(xi, yi, ni, yyp1, yypn, y2);
	}

	/*
	 * It is important to understand that the program spline is called only once
	 * to process an entire tabulated function in arrays xi and yi.
	 */
	public static double yspline = 0.0;

	/**
	 * Given the arrays xa[1..n] and ya[1..n], which tabulate a function (with the xai’s in order), and given the array y2a[1..n], 
	 * which is the output from spline above, and given a value of x, this routine returns a cubic-spline interpolated value yspline.
	 * @param xa xa
	 * @param ya ya
	 * @param y2a y2a
	 * @param n n
	 * @param x x
	 */
	public static void splint(double[] xa, double[] ya, double[] y2a, int n,
			double x)// , float *y)
	// Given the arrays xa[1..n] and ya[1..n], which tabulate a function (with
	// the xai’s in order),
	// and given the array y2a[1..n], which is the output from spline above, and
	// given a value of
	// x, this routine returns a cubic-spline interpolated value y.
	{
		// void nrerror(char error_text[]);
		int klo = 0;
		int khi = 0;
		int k = 0;
		double h = 0.0;
		double b = 0.0;
		double a = 0.0;

		failB = false;
		klo = 1;
		// We will find the right place in the table by means of
		// bisection. This is optimal if sequential calls to this
		// routine are at random values of x. If sequential calls
		// are in order, and closely spaced, one would do better
		// to store previous values of klo and khi and test if
		// they remain appropriate on the next call.
		khi = n;
		while (khi - klo > 1) {
			k = (khi + klo) >> 1;
			if (xa[k - 1] > x)
				khi = k;// if (xa[k] > x) khi=k;
			else
				klo = k;
		} // klo and khi now bracket the input value of x.
		h = xa[khi - 1] - xa[klo - 1];// h=xa[khi]-xa[klo];
		if (h == 0.0) {
			// nrerror("Bad xa input to routine splint"); The xa’s must be
			// distinct.
			failS = "Bad xa input to routine splint";
			failB = true;
		}
		a = (xa[khi - 1] - x) / h;// a=(xa[khi]-x)/h;
		b = (x - xa[klo - 1]) / h;// b=(x-xa[klo])/h; Cubic spline polynomial is
									// now evaluated.
		// *y=a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0;
		yspline = a
				* ya[klo - 1]
				+ b
				* ya[khi - 1]
				+ ((a * a * a - a) * y2a[klo - 1] + (b * b * b - b)
						* y2a[khi - 1]) * (h * h) / 6.0;
	}

	// @@@@@LOCATOR
	/**
	 * Given an array xx[1..n], and given a value x, returns a value j such that x is between xx[j] 
	 * and xx[j+1]. xx must be monotonic, either increasing or decreasing. j=0 or j=n is returned 
	 * to indicate that x is out of range.
	 * @param xx xx
	 * @param n n
	 * @param x x
	 * @return the result
	 */
	public static int locate(double[] xx, int n, double x)// , int *j)
	// Given an array xx[1..n], and given a value x, returns a value j such that
	// x is between xx[j]
	// and xx[j+1]. xx must be monotonic, either increasing or decreasing. j=0
	// or j=n is returned
	// to indicate that x is out of range.
	{
		int j = 0;// result
		int ju = 0;
		int jm = 0;
		int jl = 0;
		boolean ascnd = false;// int ascnd=0;
		jl = 0; // Initialize lower
		ju = n + 1; // and upper limits.
		ascnd = (xx[n - 1] >= xx[0]);// ascnd=(xx[n] >= xx[1]);
		while (ju - jl > 1) { // If we are not yet done,
			jm = (ju + jl) >> 1;// compute a midpoint,
			if (x >= xx[jm - 1] == ascnd)// if (x >= xx[jm] == ascnd)
				jl = jm; // and replace either the lower limit
			else
				ju = jm; // or the upper limit, as appropriate.
		}// Repeat until the test condition is satisfied.
		if (x == xx[0])
			j = 1;// if (x == xx[1]) *j=1; //Then set the output
		else if (x == xx[n - 1])
			j = n - 1;// if(x == xx[n]) *j=n-1;
		else
			j = jl;

		return j - 1;// array 0 start!!
	}

	/*
	 * Sometimes you will be in the situation of searching a large table many
	 * times, and with nearly identical abscissas on consecutive searches. For
	 * example, you may be generating a function that is used on the right-hand
	 * side of a differential equation: Most differential-equation integrators,
	 * as we shall see in Chapter 16, call for right-hand side evaluations at
	 * points that hop back and forth a bit, but whose trend moves slowly in the
	 * direction of the integration.
	 */
	// @@@hunting locators!!!
	/**
	 * Given an array xx[1..n], and given a value x, returns a value j such that x is between 
	 * xx[jlo] and xx[jlo+1]. xx[1..n] must be monotonic, either increasing or decreasing. 
	 * j=0 or j=n is returned to indicate that x is out of range. jlo on input is taken as the 
	 * initial guess for j on output.
	 * @param xx x
	 * @param n n
	 * @param x x
	 * @param jlo jlo
	 * @return the result
	 */
	public static int hunt(double xx[], int n, double x, int jlo)// , unsigned
																	// long
																	// *jlo)
	// Given an array xx[1..n], and given a value x, returns a value jlo such
	// that x is between
	// xx[jlo] and xx[jlo+1]. xx[1..n] must be monotonic, either increasing or
	// decreasing.
	// jlo=0 or jlo=n is returned to indicate that x is out of range. jlo on
	// input is taken as the
	// initial guess for jlo on output.
	{
		int j = jlo;// initial guess in [1,n] interval!!!!!
		int jm = 0;
		int jhi = 0;
		int inc = 0;
		boolean ascnd = false;// int ascnd=0;
		ascnd = (xx[n - 1] >= xx[0]);// ascnd=(xx[n] >= xx[1]); //True if
										// ascending order of table, false
										// otherwise.
		if (j <= 0 || j > n)// if (*jlo <= 0 || *jlo > n) {
		{ // Input guess not useful. Go immediately to bisection.
			j = 0;// *jlo=0;
			jhi = n + 1;
		} else {
			inc = 1; // Set the hunting increment.
			if (x >= xx[j - 1] == ascnd)// if (x >= xx[*jlo] == ascnd)
			{// Hunt up:
				if (j == n)
					return j - 1;// if (*jlo == n) return;
				jhi = (j) + 1;// jhi=(*jlo)+1;
				while (x >= xx[jhi - 1] == ascnd)// while (x >= xx[jhi] ==
													// ascnd)
				{// Not done hunting,
					j = jhi;// *jlo=jhi;
					inc += inc; // so double the increment
					jhi = (j) + inc;// jhi=(*jlo)+inc;
					if (jhi > n) {// Done hunting, since off end of table.
						jhi = n + 1;
						break;
					}// Try again.
				}// Done hunting, value bracketed.
			} else {// Hunt down:
				if (j == 1)// if (*jlo == 1)
				{
					j = 0;// *jlo=0;
					return j - 1;
				}
				jhi = (j)--;// jhi=(*jlo)--;
				while (x < xx[j - 1] == ascnd)// while (x < xx[*jlo] == ascnd)
				{ // Not done hunting,
					jhi = (j);// jhi=(*jlo);
					inc <<= 1; // so double the increment
					if (inc >= jhi) {// Done hunting, since off end of table.
						j = 0;// *jlo=0;
						break;
					} else
						j = jhi - inc;// /*jlo=jhi-inc;
				}// and try again.
			}// Done hunting, value bracketed.
		}// Hunt is done, so begin the final bisection phase:
		while (jhi - (j) != 1)// while (jhi-(*jlo) != 1)
		{
			jm = (jhi + (j)) >> 1;// jm=(jhi+(*jlo)) >> 1;
			if (x >= xx[jm - 1] == ascnd)// if (x >= xx[jm] == ascnd)
				j = jm;// *jlo=jm;
			else
				jhi = jm;
		}
		if (x == xx[n - 1])
			j = n - 1;// if (x == xx[n]) *jlo=n-1;
		if (x == xx[0])
			j = 1;// if (x == xx[1]) *jlo=1;
		return j - 1;
	}

	/*
	 * Occasionally you may wish to knownot the value of the interpolating
	 * polynomial that passes through a (small!) number of points, but the
	 * coefficients of that polynomial. A valid use of the coefficients might
	 * be, for example, to compute simultaneous interpolated values of the
	 * function and of several of its derivatives (see §5.3), or to convolve a
	 * segment of the tabulated function with some other function, where the
	 * moments of that other function (i.e., its convolution with powers of x)
	 * are known analytically. Also, you should not mistake the interpolating
	 * polynomial (and its coefficients) for its cousin, the best fit polynomial
	 * through a data set. Fitting is a smoothing process, since the number of
	 * fitted coefficients is typically much less than the number of data
	 * points. Therefore, fitted coefficients can be accurately and stably
	 * determined even in the presence of statistical errors in the tabulated
	 * values. (See §14.8.) Interpolation, where the number of coefficients and
	 * number of tabulated points are equal, takes the tabulated values as
	 * perfect. If they in fact contain statistical errors, these can be
	 * magnified into oscillations of the interpolating polynomial in between
	 * the tabulated points. As before, we take the tabulated points to be yi ?
	 * y(xi). If the interpolating polynomial is written as y = c0 + c1x + c2x^2
	 * + · · · + cNx^N (3.5.1) then the ci’s are required to satisfy the linear
	 * equation =>This is a Vandermonde matrix, as described in
	 */
	/**
	 * Given arrays x[0..n] and y[0..n] containing a tabulated function yi =f(xi), this routine 
	 * returns an array of coefficients cof[0..n], such that yi = sumj cofjxi^j.
	 * @param x x
	 * @param y y
	 * @param n n
	 * @param cof cof
	 */
	public static void polcoe(double[] x, double[] y, int n, double[] cof)
	// Given arrays x[0..n] and y[0..n] containing a tabulated function yi =
	// f(xi), this routine
	// returns an array of coefficients cof[0..n], such that yi = sumj cofjxi^j.
	{
		int k = 0;
		int j = 0;
		int i = 0;
		double phi = 0.0;
		double ff = 0.0;
		double b = 0.0;// *s;
		double[] s = new double[n + 1];// s=vector(0,n);
		for (i = 0; i <= n; i++)
			s[i] = cof[i] = 0.0;
		s[n] = -x[0];
		for (i = 1; i <= n; i++) {// Coefficients si of the master polynomial
									// P(x) are
									// found by recurrence.
			for (j = n - i; j <= n - 1; j++)
				s[j] -= x[i] * s[j + 1];
			s[n] -= x[i];
		}
		for (j = 0; j <= n; j++) {
			phi = n + 1;
			for (k = n; k >= 1; k--)
				// The quantity phi = prodj!=k (xj - xk) is found as a
				// derivative of P(xj).
				phi = k * s[k] + x[j] * phi;
			ff = y[j] / phi;
			b = 1.0; // Coefficients of polynomials in each term of the Lagrange
						// formula are found by synthetic division of
						// P(x) by (x - xj ). The solution ck is accumulated.
			for (k = n; k >= 0; k--) {
				cof[k] += b * ff;
				b = s[k] + x[j] * b;
			}
		}
		// free_vector(s,0,n);
	}

	/*
	 * Another technique is to make use of the function value interpolation
	 * routine already given (polint). If we interpolate (or extrapolate)
	 * to find the value of the interpolating polynomial at x = 0, then this
	 * value will evidently be c 0. Now we can subtract c0 from the yi’s and
	 * divide each by its corresponding xi. Throwing out one point (the one with
	 * smallest xi is a good candidate), we can repeat the procedure to find c1,
	 * and so on. It is not instantly obvious that this procedure is stable, but
	 * we have generally found it to be somewhat more stable than the routine
	 * immediately preceding.
	 */
	/**
	 * Given arrays xa[0..n] and ya[0..n] containing a tabulated function yai =f(xai), this 
	 * routine returns an array of coefficients cof[0..n] such that yai = sumjcofjxai^j .
	 * @param xa xa
	 * @param ya ya
	 * @param n n
	 * @param cof cof
	 */
	public static void polcof(double[] xa, double[] ya, int n, double[] cof)
	// Given arrays xa[0..n] and ya[0..n] containing a tabulated function yai =
	// f(xai), this
	// routine returns an array of coefficients cof[0..n] such that yai = sumj
	// cofjxai^j .
	{
		// void polint(float xa[], float ya[], int n, float x, float *y, float
		// *dy);
		int k = 0;
		int j = 0;
		int i = 0;
		double xmin = 0.0;
		// double dy = 0.0;// double *x,double *y;
		double[] x = new double[n + 1];// x=vector(0,n);
		double[] y = new double[n + 1];// y=vector(0,n);
		for (j = 0; j <= n; j++) {
			x[j] = xa[j];
			y[j] = ya[j];
		}
		for (j = 0; j <= n; j++) {
			// polint(x-1,y-1,n+1-j,0.0,&cof[j],&dy);
			polint(x, y, n + 1 - j, 0.0);// ,&cof[j],&dy);
			cof[j] = ypoli;
			// dy = dypoli;
			// Subtract 1 from the pointers to x and y because polint uses
			// dimensions [1..n]. We
			// extrapolate to x = 0.
			xmin = 1.0e38;
			k = -1;
			for (i = 0; i <= n - j; i++) { // Find the remaining xi of smallest
											// absolute value,
				if (Math.abs(x[i]) < xmin)// if (fabs(x[i]) < xmin)
				{
					xmin = Math.abs(x[i]);// xmin=fabs(x[i]);
					k = i;
				}
				// if (x[i]) y[i]=(y[i]-cof[j])/x[i]; //(meanwhile reducing all
				// the terms)
				if (x[i] != 0.0)
					y[i] = (y[i] - cof[j]) / x[i];
			}
			for (i = k + 1; i <= n - j; i++) {// and eliminate it.
				y[i - 1] = y[i];
				x[i - 1] = x[i];
			}
		}
		// free_vector(y,0,n);
		// free_vector(x,0,n);
	}

	// @@@@@@@@Interpolation in Two or More Dimensions
	/*
	 * In two dimensions, we imagine that we are given a matrix of functional
	 * values ya[1..m][1..n].We are also given an array x1a[1..m], and an array
	 * x2a[1..n]. The relation of these input quantities to an underlying
	 * function y(x 1, x2) is ya[j][k] = y(x1a[j], x2a[k]) Bilinear
	 * interpolation is frequently “close enough for government work.”
	 */
	public static double ypoli2 = 0.0;
	public static double dypoli2 = 0.0;

	/**
	 * Given arrays x1a[1..m] and x2a[1..n] of independent variables, and a submatrix of function 
	 * values ya[1..m][1..n], tabulated at the grid points defined by x1a and x2a; and given values 
	 * x1 and x2 of the independent variables; this routine returns an interpolated function value ypoli2, 
	 * and an accuracy indication dypoli2 (based only on the interpolation in the x1 direction, however).
	 * @param x1a x1a
	 * @param x2a x2a
	 * @param ya ya
	 * @param m m
	 * @param n n
	 * @param x1 x1
	 * @param x2 x2
	 */
	public static void polin2(double[] x1a, double[] x2a, double[][] ya, int m,
			int n, double x1, double x2)// , float *y, float *dy)
	// Given arrays x1a[1..m] and x2a[1..n] of independent variables, and a
	// submatrix of function
	// values ya[1..m][1..n], tabulated at the grid points defined by x1a and
	// x2a; and given values
	// x1 and x2 of the independent variables; this routine returns an
	// interpolated function value y,
	// and an accuracy indication dy (based only on the interpolation in the x1
	// direction, however).
	{
		// void polint(float xa[], float ya[], int n, float x, float *y, float
		// *dy);
		int j = 0;
		double[] ymtmp = new double[m];
		// ymtmp=vector(1,m);
		for (j = 1; j <= m; j++) { // Loop over rows.Interpolate answer into
									// temporary storage.
			double[] yy = new double[n];
			for (int k = 1; k <= n; k++)
				yy[k - 1] = ya[j - 1][k - 1];
			polint(x2a, yy, n, x2);// polint(x2a,ya[j],n,x2);//,&ymtmp[j],dy);
			ymtmp[j - 1] = ypoli;
			dypoli2 = dypoli;
		}
		polint(x1a, ymtmp, m, x1);// ,y,dy); //Do the final interpolation.
		ypoli2 = ypoli;
		dypoli2 = dypoli;
		// free_vector(ymtmp,1,m);
	}

	/*
	 * We will give two methods that are in common use, and which are themselves
	 * not unrelated. The first is usually called bicubic interpolation. Best of
	 * all is to know the derivatives analytically, or to be able to compute
	 * them accurately by numerical means, at the grid points. Next best is to
	 * determine them by numerical differencing from the functional values
	 * already tabulated on the grid. The relevant code would be something like
	 * this (using centered differencing):
	 * y1a[j][k]=(ya[j+1][k]-ya[j-1][k])/(x1a[j+1]-x1a[j-1]);
	 * y2a[j][k]=(ya[j][k+1]-ya[j][k-1])/(x2a[k+1]-x2a[k-1]);
	 * y12a[j][k]=(ya[j+1][k+1]-ya[j+1][k-1]-ya[j-1][k+1]+ya[j-1][k-1])
	 * /((x1a[j+1]-x1a[j-1])*(x2a[k+1]-x2a[k-1]));
	 */
	/**
	 * Given arrays y[1..4], y1[1..4], y2[1..4], and y12[1..4], containing the function, gradients, 
	 * and cross derivative at the four grid points of a rectangular grid cell (numbered counterclockwise 
	 * from the lower left), and given d1 and d2, the length of the grid cell in the 1- and 2-directions, 
	 * this routine returns the table c[1..4][1..4] that is used by routine bcuint for bicubic interpolation.
	 * @param y y
	 * @param y1 y1
	 * @param y2 y2
	 * @param y12 y12
	 * @param d1 d1
	 * @param d2 d2
	 * @param c c
	 */
	public static void bcucof(double[] y, double[] y1, double[] y2,
			double[] y12, double d1, double d2, double[][] c)
	// Given arrays y[1..4], y1[1..4], y2[1..4], and y12[1..4], containing the
	// function, gradients,
	// and cross derivative at the four grid points of a rectangular grid cell
	// (numbered counterclockwise
	// from the lower left), and given d1 and d2, the length of the grid cell in
	// the 1- and
	// 2-directions, this routine returns the table c[1..4][1..4] that is used
	// by routine bcuint
	// for bicubic interpolation.
	{
		// static int wt[16][16]=
		int[][] wt = // [16][16]=
		{ { 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0 },
				{ -3, 0, 0, 3, 0, 0, 0, 0, -2, 0, 0, -1, 0, 0, 0, 0 },
				{ 2, 0, 0, -2, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 },
				{ 0, 0, 0, 0, -3, 0, 0, 3, 0, 0, 0, 0, -2, 0, 0, -1 },
				{ 0, 0, 0, 0, 2, 0, 0, -2, 0, 0, 0, 0, 1, 0, 0, 1 },
				{ -3, 3, 0, 0, -2, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 0, 0, 0, 0, -3, 3, 0, 0, -2, -1, 0, 0 },
				{ 9, -9, 9, -9, 6, 3, -3, -6, 6, -6, -3, 3, 4, 2, 1, 2 },
				{ -6, 6, -6, 6, -4, -2, 2, 4, -3, 3, 3, -3, -2, -1, -1, -2 },
				{ 2, -2, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 },
				{ 0, 0, 0, 0, 0, 0, 0, 0, 2, -2, 0, 0, 1, 1, 0, 0 },
				{ -6, 6, -6, 6, -3, -3, 3, 3, -4, 4, 2, -2, -2, -2, -1, -1 },
				{ 4, -4, 4, -4, 2, 2, -2, -2, 2, -2, -2, 2, 1, 1, 1, 1 } };
		int l = 0;
		int k = 0;
		int j = 0;
		int i = 0;
		double xx = 0.0;
		double d1d2 = 0.0;
		double[] cl = new double[16];
		double[] x = new double[16];
		d1d2 = d1 * d2;
		for (i = 1; i <= 4; i++) {// Pack a temporary vector x.
			x[i - 1] = y[i - 1];// y[i];
			x[i + 3] = y1[i - 1] * d1;// y1[i]*d1;
			x[i + 7] = y2[i - 1] * d2;// y2[i]*d2;
			x[i + 11] = y12[i - 1] * d1d2;// y12[i]*d1d2;
		}
		for (i = 0; i <= 15; i++) {// Matrix multiply by the stored table.
			xx = 0.0;
			for (k = 0; k <= 15; k++)
				xx += wt[i][k] * x[k];
			cl[i] = xx;
		}
		l = 0;
		for (i = 1; i <= 4; i++)
			// Unpack the result into the output table.
			for (j = 1; j <= 4; j++)
				c[i - 1][j - 1] = cl[l++];// c[i][j]=cl[l++];
	}

	public static double ansy = 0.0;
	public static double ansy1 = 0.0;
	public static double ansy2 = 0.0;

	/**
	 * Bicubic interpolation within a grid square. Input quantities are y,y1,y2,y12 (as described in bcucof); 
	 * x1l and x1u, the lower and upper coordinates of the grid square in the 1-direction; 
	 * x2l and x2u likewise for the 2-direction; and x1,x2, the coordinates of the desired point for 
	 * the interpolation. The interpolated function value is returned as ansy, and the interpolated 
	 * gradient values as ansy1 and ansy2. This routine calls bcucof.
	 * @param y y
	 * @param y1 y1
	 * @param y2 y2
	 * @param y12 y12
	 * @param x1l x1l
	 * @param x1u x1u
	 * @param x2l x2l
	 * @param x2u x2u
	 * @param x1 x1
	 * @param x2 x2
	 */
	public static void bcuint(double[] y, double[] y1, double[] y2,
			double[] y12, double x1l, double x1u, double x2l, double x2u,
			double x1, double x2)// , float *ansy,float *ansy1, float *ansy2)
	// Bicubic interpolation within a grid square. Input quantities are
	// y,y1,y2,y12 (as described in
	// bcucof); x1l and x1u, the lower and upper coordinates of the grid square
	// in the 1-direction;
	// x2l and x2u likewise for the 2-direction; and x1,x2, the coordinates of
	// the desired point for
	// the interpolation. The interpolated function value is returned as ansy,
	// and the interpolated
	// gradient values as ansy1 and ansy2. This routine calls bcucof.
	{
		// void bcucof(float y[], float y1[], float y2[], float y12[], float d1,
		// float d2, float **c);
		int i = 0;
		double t = 0.0;
		double u = 0.0;
		double d1 = 0.0;
		double d2 = 0.0;
		double[][] c = new double[4][4];
		// c=matrix(1,4,1,4);
		failB = false;
		d1 = x1u - x1l;
		d2 = x2u - x2l;
		bcucof(y, y1, y2, y12, d1, d2, c); // Get the c’s.
		if (x1u == x1l || x2u == x2l) {
			// nrerror("Bad input in routine bcuint");
			failS = "Bad input in routine bcuint";
			failB = true;
			return;
		}
		t = (x1 - x1l) / d1; // Equation (3.6.4).
		u = (x2 - x2l) / d2;
		// *ansy=(*ansy2)=(*ansy1)=0.0;
		ansy = (ansy2) = (ansy1) = 0.0;
		for (i = 4; i >= 1; i--) {// Equation (3.6.6).
									// *ansy=t*(*ansy)+((c[i][4]*u+c[i][3])*u+c[i][2])*u+c[i][1];
									// *ansy2=t*(*ansy2)+(3.0*c[i][4]*u+2.0*c[i][3])*u+c[i][2];
									// *ansy1=u*(*ansy1)+(3.0*c[4][i]*t+2.0*c[3][i])*t+c[2][i];

			ansy = t * (ansy)
					+ ((c[i - 1][3] * u + c[i - 1][2]) * u + c[i - 1][1]) * u
					+ c[i - 1][0];
			ansy2 = t * (ansy2) + (3.0 * c[i - 1][3] * u + 2.0 * c[i - 1][2])
					* u + c[i - 1][1];
			ansy1 = u * (ansy1) + (3.0 * c[3][i - 1] * t + 2.0 * c[2][i - 1])
					* t + c[1][i - 1];
		}
		ansy1 /= d1;
		ansy2 /= d2;
		// free_matrix(c,1,4,1,4);
	}

	/**
	 * Given an m by n tabulated function ya[1..m][1..n], and tabulated independent variables x1a[1..n] and
	 * x2a[1..n], this routine constructs one-dimensional natural cubic splines of the rows of ya 
	 * and returns the second-derivatives in the array y2a[1..m][1..n]. (The array x1a[1..m] is 
	 * included in the argument list merely for consistency with routine splin2.)
	 * @param x1a x1a
	 * @param x2a x2a
	 * @param ya ya
	 * @param m m
	 * @param n n
	 * @param y2a y2a
	 */
	public static void splie2(double[] x1a, double[] x2a, double[][] ya, int m,
			int n, double[][] y2a)
	// Given an m by n tabulated function ya[1..m][1..n], and tabulated
	// independent variables
	// x2a[1..n], this routine constructs one-dimensional natural cubic splines
	// of the rows of ya
	// and returns the second-derivatives in the array y2a[1..m][1..n]. (The
	// array x1a[1..m] is
	// included in the argument list merely for consistency with routine
	// splin2.)
	// double[] x, double[] y, int n, double yp1, double ypn, double[] y2)
	// Given arrays x[1..n] and y[1..n] containing a tabulated function, i.e.,
	// yi = f(xi), with
	// x1 < x2 < .. . < xN, and given values yp1 and ypn for the first
	// derivative of the interpolating
	// function at points 1 and n, respectively, this routine returns an array
	// y2[1..n] that contains
	// the second derivatives of the interpolating function at the tabulated
	// points xi.
	{
		// void spline(float x[], float y[], int n, float yp1, float ypn, float
		// y2[]);
		int j = 0;
		for (j = 1; j <= m; j++) {
			double[] yy = new double[n];
			double[] yy2 = new double[n];
			for (int k = 1; k <= n; k++) {
				yy[k - 1] = ya[j - 1][k - 1];
				yy2[k - 1] = 0.0;
			}
			// spline(x2a,ya[j],n,1.0e30,1.0e30,y2a[j]); //Values 1×10^30 signal
			// a natural spline.
			spline(x2a, yy, n, 1.0e30, 1.0e30, yy2);

			for (int k = 1; k <= n; k++) {
				y2a[j - 1][k - 1] = yy2[k - 1];
			}
		}
	}

	// After the above routine has been executed once, any number of bicubic
	// spline
	// interpolations can be performed by successive calls of the following
	// routine:
	/*
	 * public static void splint(double[] xa, double[] ya, double[] y2a, int n,
	 * double x)//, float *y) //Given the arrays xa[1..n] and ya[1..n], which
	 * tabulate a function (with the xai’s in order), //and given the array
	 * y2a[1..n], which is the output from spline above, and given a value of
	 * //x, this routine returns a cubic-spline interpolated value y.
	 */
	public static double yspline2 = 0.0;

	/**
	 * Given x1a, x2a, ya, m, n as described in splie2 and y2a as produced by that routine; and 
	 * given a desired interpolating point x1,x2; this routine returns an interpolated function value yspline2 
	 * by bicubic spline interpolation.
	 * @param x1a x1a
	 * @param x2a x2a
	 * @param ya ya
	 * @param y2a y2a
	 * @param m m
	 * @param n n
	 * @param x1 x1
	 * @param x2 x2
	 */
	public static void splin2(double[] x1a, double[] x2a, double[][] ya,
			double[][] y2a, int m, int n, double x1, double x2)// , double[] y)
	// Given x1a, x2a, ya, m, n as described in splie2 and y2a as produced by
	// that routine; and
	// given a desired interpolating point x1,x2; this routine returns an
	// interpolated function value y
	// by bicubic spline interpolation.
	{
		// void spline(float x[], float y[], int n, float yp1, float ypn, float
		// y2[]);
		// void splint(float xa[], float ya[], float y2a[], int n, float x,
		// float *y);
		int j = 0;
		// float *ytmp,*yytmp;
		double[] ytmp = new double[m];// vector(1,m);
		double[] yytmp = new double[m];// vector(1,m);
		// Perform m evaluations of the row splines constructed by
		// splie2, using the one-dimensional spline evaluator splint.
		for (j = 1; j <= m; j++) {
			double[] yy = new double[n];
			double[] yy2 = new double[n];
			for (int k = 1; k <= n; k++) {
				yy[k - 1] = ya[j - 1][k - 1];
				yy2[k - 1] = y2a[j - 1][k - 1];
			}

			// splint(x2a,ya[j],y2a[j],n,x2,&yytmp[j]);
			splint(x2a, yy, yy2, n, x2);// ,&yytmp[j]);
			yytmp[j - 1] = yspline;
		}

		spline(x1a, yytmp, m, 1.0e30, 1.0e30, ytmp);
		// Construct the one-dimensional column spline and evaluate it.
		splint(x1a, yytmp, ytmp, m, x1);// ,y);
		yspline2 = yspline;
		// free_vector(yytmp,1,m);
		// free_vector(ytmp,1,m);
	}

	/**
	 * Sets cubic spline interpolation coefficients for the data contained in the array f(n) at the abscissas x(n). 
	 * This is EGSnrc implementation of spline method.   
	 * @param x x
	 * @param f f
	 * @param a a
	 * @param b b
	 * @param c c
	 * @param d d
	 * @param n n
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

		// aa=new double[n];bb=new double[n];cc=new double[n];dd=new
		// double[n];xx=new double[n];

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
	 * Returns the value of the function at s using the spline coefficients ,b,c,d, which must have been set using set_spline. 
	 * 
	 * @param s s
	 * @param x x
	 * @param a a
	 * @param b b
	 * @param c c
	 * @param d d
	 * @param n n
	 * @return the result
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
	//====================
	//intra sirurile coordonatelor x si y, ordinul polinomului de interpolare dorit n
	//rezultatul se stocheaza in matricea (de fapt vectorul) tl!!
	/**
	 * Another polynomial method based on least square method. x, y are tabulated arrays y=f(x), n is the polynomial order and tl is a 1 dimensional 
	 * matrix (basically a vector) storing the polynomial coefficients ai:
	 * P(x) = a0 + a1x + a2x^2 + a3x^3 + ... + anx^n.
	 * Least square means minimizing the chi2, that is the quantity S=SUM[y-P(x)]^2. dS/dai=0 yields a system of equations!
	 * @param x x
	 * @param y y
	 * @param n n
	 * @param tl tl
	 */
	public static void polynomial(double[] x, double[] y, int n, double[][] tl)
	{
	//construiesc sirul termenilor liberi ai sistemului de ecuatii obtinut prin metoda celor mai mici patrate!
	    double tls,cf = 0.0;
        for (int j=0; j<=n; j++)
        {
            tls=0.0;
            for (int i =0; i<x.length; i++)
               tls=tls+y[i]*Math.pow(x[i],j);
            tl[j][0]=tls;
		}
	//construiesc sirul elementelor matricii sistemului de ecuatii
	    double[][] coef = new double[n+1][n+1];
	    for (int i=0; i<=n; i++)
	       for (int j=0; j<=n; j++)
	       {
			   cf=0.0;
			   for (int k=0; k<x.length; k++)
			      cf=cf+Math.pow(x[k],i+j);
			   coef[i][j]=cf;
		   }
	//apelez rezolvarea sistemului de ecuatii
	    LinAEq.sysEqGauss(coef,tl,n+1,1);
	    //returneaza pentru ecuatia de interpolare in aceasta ordine:
	    //a0 + a1x + a2x^2 + a3x^3 + ... + anx^n
	//returneaza solutia in tabloul tl!!
	}

    //intra coordonatele xi,yi ale celor ni puncte de retea
    //tabloul x pentru care se doreste evaluarea
    //tabloul y al evaluarii
    //a,b,c,d--tabloul restricttlor pe intervake ca functii:ax^3+bx^2+cx+d
    //dm,dp --valorile derivatelor la capete SAU daca nu se cunosc
    //pot fi evaluate pe baza interpolarii Lagrange
	/**
	 * Another implementation of spline method. xi,yi are the known point coordinates, ni is the number of points. a,b,c,d are 
	 * the spline coefficients which are also evaluated here (the restrictions ax^3+bx^2+cx+d between each successive points), dm,dp are the derivatives 
	 * at both ends which, if are unknown, can be evaluated using Langrange method. x and y are desired points for evaluation.
	 * @param xi xi
	 * @param yi yi
	 * @param ni ni
	 * @param x x
	 * @param y y
	 * @param n n
	 * @param a a
	 * @param b b
	 * @param c c
	 * @param d d
	 * @param dm dm
	 * @param dp dp
	 */
	public static void spline(double[] xi, double[] yi, int ni, double[] x, double[] y, int n,
	double[] a, double[] b, double[] c, double[] d, double dm, double dp)
	{
         double dd=0.0;
         double hh=0.0;
         double hm=0.0;
         double xx=0.0;
         double xp=0.0;
         int ip=0;

         //hh=0;
		 dd=dm;
		 for(int i=0; i<=ni-2; i++)
		 {
		      hm=hh;
		      hh=xi[i+1]-xi[i];
		      dm=dd;
		      dd=(yi[i+1]-yi[i])/hh;
		      a[i]=hm;
		      b[i]=2*(hm+hh);
		      c[i]=hh;
		      d[i]=6*(dd-dm);
		 }
		 a[ni-1]=hh;
		 b[ni-1]=2*hh;
		 c[ni-1]=0.0;
		 d[ni-1]=6*(dp-dd);
		 //rezolv sistemul
		 LinAEq.sysEqTriDiag(a,b,c,d,ni);
		 ///
		 //calculez coeficientii functiilor spline cubice pe fiecare interval din reteaua de puncte data
		 for(int i=0; i<=ni-2; i++)
		 {
		      ip=i+1;
		      xx=xi[i];
		      xp=xi[ip];
		      hh=xp-xx; //pasul
		      dd=d[i];  //aci sunt derivatele calculate cu tridiag
		      dp=d[ip];
		      a[i]=(dp-dd)/(6*hh);
		      b[i]=(dd*xp-dp*xx)/(2*hh);
		      c[i]=(dp*xx*xx-dd*xp*xp)/(2*hh)+(yi[ip]-yi[i])/hh-a[i]*hh*hh;
		      d[i]=(dd*xp*xp*xp-dp*xx*xx*xx)/(6*hh)+(yi[i]*xp-yi[ip]*xx)/hh-b[i]*hh*hh/3;
		}
		///pentru x dat ca vector evaluez spline coresp. si dau y=solutia
		ip=0;
		for(int i=0; i<=n-1; i++)
		{
		      xx=x[i];
		      while (xx>xi[ip+1])  //calculez intii restrictia
		         ip++;
		      y[i]=((a[ip]*xx+b[ip])*xx+c[ip])*xx+d[ip]; //pe baza restrictiei gasite
        }
	}
}
