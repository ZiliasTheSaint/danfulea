package danfulea.math.numerical;

/**
 * Fast Fourier analysis.
 * Based on literature, including Numerical Recipes in C (Cambridge Univ.).
 * @author Dan Fulea, 15 OCT. 2006
 */
public class FFTAnalysis {
	// public static Function func;
	public static boolean failB = false;
	public static String failS = "";
	public static double xms_memcof = 0.0;

	public static int NMAX = 100;// Largest expected value of m.
	public static Complex ZERO = new Complex(0.0, 0.0);
	public static Complex ONE = new Complex(1.0, 0.0);
	public static double TWOPID = 6.2831853071795865;

	public static int nout_period = 0;
	public static int jmax_period = 0;
	public static double prob_period = 0;

	public static double ave_avevar = 0.0;
	public static double var_avevar = 0.0;

	public static int MACC = 4;
	public static int nout_fasper = 0;
	public static int jmax_fasper = 0;
	public static double prob_fasper = 0.0;

	public static double corre_dftcor = 0.0;
	public static double corim_dftcor = 0.0;
	public static double corfac_dftcor = 0.0;

	public static int M = 64;
	public static int NDFT = 1024;
	public static int MPOL = 6;
	public static double TWOPI = (2.0 * 3.14159265);
	private Function func;
	public static double sinint_dftint = 0.0;
	public static double cosint_dftint = 0.0;

	public static int init = 0;
	public static double aold = -1.e30;
	public static double bold = -1.e30;
	public static double delta = 0.0;
	// public static double (*funcold)(float);
	public static double[] data = new double[NDFT + 1];
	public static double[] endpts = new double[9];

	public static double C0 = 0.4829629131445341;
	public static double C1 = 0.8365163037378079;
	public static double C2 = 0.2241438680420134;
	public static double C3 = -0.1294095225512604;

	public static double[] c4r = new double[5];
	public static double[] c12r = new double[13];
	public static double[] c20r = new double[21];

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
	 * Constructor. Here the function is passed by a class implementing Function interface.
	 * @param func the function
	 */
	public FFTAnalysis(Function func) {
		this.func = func;
	}

	/*
	 * Fourier methods have revolutionized fields of science and engineering,
	 * from radio astronomy to medical imaging, from seismology to spectroscopy.
	 * In this chapter, we present some of the basic applications of Fourier and
	 * spectral methods that have made these revolutions possible. Say the word
	 * “Fourier” to a numericist, and the response, as if by Pavlovian
	 * conditioning, will likely be “FFT.” Indeed, the wide application of
	 * Fourier methods must be credited principally to the existence of the fast
	 * Fourier transform. Better mousetraps stand aside: If you speed up any
	 * nontrivial algorithm by a factor of a million or so, the world will beat
	 * a path towards finding useful applications for it. The most direct
	 * applications of the FFT are to the convolution or deconvolution of data
	 * (§13.1), correlation and autocorrelation (§13.2), optimal filtering
	 * (§13.3), power spectrum estimation (§13.4), and the computation of
	 * Fourier integrals (§13.9). As important as they are, however, FFT methods
	 * are not the be-all and end-all of spectral analysis. Section 13.5 is a
	 * brief introduction to the field of time-domain digital filters. In the
	 * spectral domain, one limitation of the FFT is that it always represents a
	 * function’s Fourier transform as a polynomial in z = exp(2PI*i*f*delta)
	 * (cf. equation 12.1.7). Sometimes, processes have spectra whose shapes are
	 * not well represented by this form. An alternative form, which allows the
	 * spectrum to have poles in z, is used in the techniques of linear
	 * prediction (§13.6) and maximum entropy spectral estimation (§13.7).
	 * Another significant limitation of all FFT methods is that they require
	 * the input data to be sampled at evenly spaced intervals. For irregularly
	 * or incompletely sampled data, other (albeit slower) methods are
	 * available, as discussed in §13.8. So-called wavelet methods inhabit a
	 * representation of function space that is neither in the temporal, nor in
	 * the spectral, domain, but rather something in-between. Section 13.10 is
	 * an introduction to this subject. Finally §13.11 is an excursion into
	 * numerical use of the Fourier sampling theorem.
	 * 
	 * 
	 * Convolution and Deconvolution Using FFT
	 * 
	 * Here is our routine for convolution and deconvolution, using the FFT as
	 * implemented in four1 of §12.2. Since the data and response functions are
	 * real, not complex, both of their transforms can be taken simultaneously
	 * using twofft. Note, however, that two calls to realft should be
	 * substituted if data and respns have very different magnitudes, to
	 * minimize roundoff. The data are assumed to be stored in a float array
	 * data[1..n], with n an integer power of two. The response function is
	 * assumed to be stored in wrap-around order in a sub-array respns[1..m] of
	 * the array respns[1..n]. The value of m can be any odd integer less than
	 * or equal to n, since the first thing the program does is to recopy the
	 * response function into the appropriate wrap-around order in respns[1..n].
	 * The answer is provided in ans.
	 */

	// public static void convlv(float data[], unsigned long n, float respns[],
	// unsigned long m,
	// int isign, float ans[])
	/**
	 * Convolves or deconvolves a real data set data[1..n] (including any user-supplied zero padding) 
	 * with a response function respns[1..n]. The response function must be stored in wrap-around 
	 * order in the first m elements of respns, where m is an odd integer = n. Wrap-around order 
	 * means that the first half of the array respns contains the impulse response function at positive 
	 * times, while the second half of the array contains the impulse response function at negative times, 
	 * counting down from the highest element respns[m]. On input isign is +1 for convolution, 
	 * -1 for deconvolution. The answer is returned in the first n components of ans. However, 
	 * ans must be supplied in the calling program with dimensions [1..2*n], for consistency with 
	 * twofft. n MUST be an integer power of two.
	 * @param data data
	 * @param n n
	 * @param respns respns
	 * @param m m
	 * @param isign isign
	 * @param ans ans
	 */
	public static void convlv(double[] data, int n, double[] respns, int m,
			int isign, double[] ans)
	// Convolves or deconvolves a real data set data[1..n] (including any
	// user-supplied zero padding)
	// with a response function respns[1..n]. The response function must be
	// stored in wrap-around
	// order in the first m elements of respns, where m is an odd integer = n.
	// Wrap-around order
	// means that the first half of the array respns contains the impulse
	// response function at positive
	// times, while the second half of the array contains the impulse response
	// function at negative times,
	// counting down from the highest element respns[m]. On input isign is +1
	// for convolution,
	// -1 for deconvolution. The answer is returned in the first n components of
	// ans. However,
	// ans must be supplied in the calling program with dimensions [1..2*n], for
	// consistency with
	// twofft. n MUST be an integer power of two.
	{
		// void realft(float data[], unsigned long n, int
		// isign);FastFourierTransform
		// void twofft(float data1[], float data2[], float fft1[], float fft2[],
		// unsigned long n);
		// unsigned long i,no2;
		failB = false;
		int i = 0;
		int no2 = 0;
		double dum = 0.0;
		double mag2 = 0.0;
		double[] fft = new double[n << 1];
		// fft=vector(1,n<<1);
		for (i = 1; i <= (m - 1) / 2; i++)
			// Put respns in array of length n.
			respns[n - i] = respns[m - i];// respns[n+1-i]=respns[m+1-i];
		for (i = (m + 3) / 2; i <= n - (m - 1) / 2; i++)
			// Pad with zeros.
			respns[i - 1] = 0.0;// respns[i]=0.0;
		FastFourierTransform.twofft(data, respns, fft, ans, n); // FFT both at
																// once.
		no2 = n >> 1;
		for (i = 2; i <= n + 2; i += 2) {
			if (isign == 1) {
				// ans[i-1]=(fft[i-1]*(dum=ans[i-1])-fft[i]*ans[i])/no2;
				// Multiply FFTs to convolve.
				ans[i - 2] = (fft[i - 2] * (dum = ans[i - 2]) - fft[i - 1]
						* ans[i - 1])
						/ no2;
				// ans[i]=(fft[i]*dum+fft[i-1]*ans[i])/no2;
				ans[i - 1] = (fft[i - 1] * dum + fft[i - 2] * ans[i - 1]) / no2;
			} else if (isign == -1) {
				// if ((mag2=SQR(ans[i-1])+SQR(ans[i])) == 0.0)
				if ((mag2 = (ans[i - 2]) * (ans[i - 2]) + (ans[i - 1])
						* (ans[i - 1])) == 0.0) {
					failB = true;
					failS = "Deconvolving at response zero in convlv";
					return;
					// nrerror("Deconvolving at response zero in convlv");
				}
				// ans[i-1]=(fft[i-1]*(dum=ans[i-1])+fft[i]*ans[i])/mag2/no2;Divide
				// FFTs to deconvolve.
				ans[i - 2] = (fft[i - 2] * (dum = ans[i - 2]) + fft[i - 1]
						* ans[i - 1])
						/ mag2 / no2;
				// ans[i]=(fft[i]*dum-fft[i-1]*ans[i])/mag2/no2;
				ans[i - 1] = (fft[i - 1] * dum - fft[i - 2] * ans[i - 1])
						/ mag2 / no2;
			} else {
				// nrerror("No meaning for isign in convlv");
				failB = true;
				failS = "No meaning for isign in convlv";
				return;
			}
		}
		// ans[2]=ans[n+1]; Pack last element with first for realft.
		ans[1] = ans[n];
		FastFourierTransform.realft(ans, n, -1); // Inverse transform back to
													// time domain.
		// free_vector(fft,1,n<<1);
	}

	/*
	 * Correlation and Autocorrelation Using FFT Correlation is the close
	 * mathematical cousin of convolution. It is in some ways simpler, however,
	 * because the two functions that go into a correlation are not as
	 * conceptually distinct as were the data and response functions that
	 * entered into convolution. Rather, in correlation, the functions are
	 * represented by different, but generally similar, data sets. We
	 * investigate their “correlation,” by comparing them both directly
	 * superposed, and with one of them shifted left or right.
	 */
	// public static void correl(double[] data1, double[] data2, unsigned long
	// n, double[] ans)
	/**
	 * Computes the correlation of two real data sets data1[1..n] and data2[1..n] (including any 
	 * user-supplied zero padding). n MUST be an integer power of two. The answer is returned as 
	 * the first n points in ans[1..2*n] stored in wrap-around order, i.e., correlations at increasingly 
	 * negative lags are in ans[n] on down to ans[n/2+1], while correlations at increasingly positive 
	 * lags are in ans[1] (zero lag) on up to ans[n/2]. Note that ans must be supplied in the calling 
	 * program with length at least 2*n, since it is also used as working space. Sign convention of 
	 * this routine: if data1 lags data2, i.e., is shifted to the right of it, then ans will show a peak 
	 * at positive lags.
	 * @param data1 data1
	 * @param data2 data2
	 * @param n n
	 * @param ans ans
	 */
	public static void correl(double[] data1, double[] data2, int n,
			double[] ans)
	// Computes the correlation of two real data sets data1[1..n] and
	// data2[1..n] (including any
	// user-supplied zero padding). n MUST be an integer power of two. The
	// answer is returned as
	// the first n points in ans[1..2*n] stored in wrap-around order, i.e.,
	// correlations at increasingly
	// negative lags are in ans[n] on down to ans[n/2+1], while correlations at
	// increasingly positive
	// lags are in ans[1] (zero lag) on up to ans[n/2]. Note that ans must be
	// supplied in the calling
	// program with length at least 2*n, since it is also used as working space.
	// Sign convention of
	// this routine: if data1 lags data2, i.e., is shifted to the right of it,
	// then ans will show a peak
	// at positive lags.
	{
		// void realft(float data[], unsigned long n, int isign);
		// void twofft(float data1[], float data2[], float fft1[], float fft2[],
		// unsigned long n);
		// unsigned long no2,i;
		int no2 = 0;
		int i = 0;
		double dum = 0.0;
		double[] fft = new double[n << 1];
		// fft=vector(1,n<<1);
		FastFourierTransform.twofft(data1, data2, fft, ans, n);// Transform both
																// data vectors
																// at once.
		no2 = n >> 1; // Normalization for inverse FFT.
		for (i = 2; i <= n + 2; i += 2) {
			// ans[i-1]=(fft[i-1]*(dum=ans[i-1])+fft[i]*ans[i])/no2; //Multiply
			// to find FFT of their correlation.
			ans[i - 2] = (fft[i - 2] * (dum = ans[i - 2]) + fft[i - 1]
					* ans[i - 1])
					/ no2;
			// ans[i]=(fft[i]*dum-fft[i-1]*ans[i])/no2;
			ans[i - 1] = (fft[i - 1] * dum - fft[i - 2] * ans[i - 1]) / no2;
		}
		ans[1] = ans[n];// ans[2]=ans[n+1]; Pack first and last into one
						// element.
		FastFourierTransform.realft(ans, n, -1);// Inverse transform gives
												// correlation.
		// free_vector(fft,1,n<<1);
	}

	/*
	 * Power Spectrum Estimation Using the FFT In the previous section we
	 * “informally” estimated the power spectral density of a function c(t) by
	 * taking the modulus-squared of the discrete Fourier transform of some
	 * finite, sampled stretch of it. In this section we’ll do roughly the same
	 * thing, but with considerably greater attention to details. Our attention
	 * will uncover some surprises. The first detail is power spectrum (also
	 * called a power spectral density or PSD) normalization. In general there
	 * is some relation of proportionality between a measure of the squared
	 * amplitude of the function and a measure of the amplitude of the PSD.
	 * 
	 * #define WINDOW(j,a,b) (1.0-fabs((((j)-1)-(a))*(b))) // Bartlett // //
	 * #define WINDOW(j,a,b) 1.0 // // Square // // #define WINDOW(j,a,b)
	 * (1.0-SQR((((j)-1)-(a))*(b))) // // Welch // void spctrm(FILE *fp, float
	 * p[], int m, int k, int ovrlap) Reads data from input stream specified by
	 * file pointer fp and returns as p[j] the data’s power (mean square
	 * amplitude) at frequency (j-1)/(2*m) cycles per gridpoint, for
	 * j=1,2,...,m, based on (2*k+1)*m data points (if ovrlap is set true (1))
	 * or 4*k*m data points (if ovrlap is set false (0)). The number of segments
	 * of the data is 2*k in both cases: The routine calls four1 k times, each
	 * call with 2 partitions each of 2*m real data points. { void four1(float
	 * data[], unsigned long nn, int isign); int
	 * mm,m44,m43,m4,kk,joffn,joff,j2,j; float
	 * w,facp,facm,*w1,*w2,sumw=0.0,den=0.0; mm=m+m; Useful factors.
	 * m43=(m4=mm+mm)+3; m44=m43+1; w1=vector(1,m4); w2=vector(1,m); facm=m;
	 * facp=1.0/m; for (j=1;j<=mm;j++) sumw += SQR(WINDOW(j,facm,facp));
	 * Accumulate the squared sum of the weights. for (j=1;j<=m;j++) p[j]=0.0;
	 * Initialize the spectrum to zero. if (ovrlap) Initialize the “save”
	 * half-buffer. for (j=1;j<=m;j++) fscanf(fp,"%f",&w2[j]); for
	 * (kk=1;kk<=k;kk++) { Loop over data set segments in groups of two. for
	 * (joff = -1;joff<=0;joff++) { Get two complete segments into workspace. if
	 * (ovrlap) { for (j=1;j<=m;j++) w1[joff+j+j]=w2[j]; for (j=1;j<=m;j++)
	 * fscanf(fp,"%f",&w2[j]); joffn=joff+mm; for (j=1;j<=m;j++)
	 * w1[joffn+j+j]=w2[j]; } else { for (j=joff+2;j<=m4;j+=2)
	 * fscanf(fp,"%f",&w1[j]); } } for (j=1;j<=mm;j++) { Apply the window to the
	 * data. j2=j+j; w=WINDOW(j,facm,facp); w1[j2] *= w; w1[j2-1] *= w; }
	 * four1(w1,mm,1); Fourier transform the windowed data. p[1] +=
	 * (SQR(w1[1])+SQR(w1[2])); Sum results into previous segments. for
	 * (j=2;j<=m;j++) { j2=j+j; p[j] += (SQR(w1[j2])+SQR(w1[j2-1])
	 * +SQR(w1[m44-j2])+SQR(w1[m43-j2])); } den += sumw; } den *= m4; Correct
	 * normalization. for (j=1;j<=m;j++) p[j] /= den; Normalize the output.
	 * free_vector(w2,1,m); free_vector(w1,1,m4); }
	 * 
	 * 
	 * 
	 * 
	 * 
	 * 
	 * 
	 * 
	 * Linear Prediction and Linear Predictive Coding
	 */
	// public static void memcof(float data[], int n, int m, float *xms, float
	// d[])
	/**
	 * Given a real vector of data[1..n], and given m, this routine returns m linear prediction (LP)
	 * coefficients as d[1..m], and returns the mean square discrepancy as xms. 
	 * @param data data
	 * @param n n
	 * @param m m
	 * @param d d
	 */
	public static void memcof(double[] data, int n, int m, double[] d)
	// Given a real vector of data[1..n], and given m, this routine returns m
	// linear prediction coef-
	// ficients as d[1..m], and returns the mean square discrepancy as xms.
	{
		int k = 0;
		int j = 0;
		int i = 0;
		double p = 0.0;
		double[] wk1 = new double[n];
		double[] wk2 = new double[n];
		double[] wkm = new double[m];
		failB = false;
		// wk1=vector(1,n);
		// wk2=vector(1,n);
		// wkm=vector(1,m);
		// for (j=1;j<=n;j++) p += SQR(data[j]);
		for (j = 1; j <= n; j++)
			p += (data[j - 1]) * (data[j - 1]);
		xms_memcof = p / n;
		wk1[0] = data[0];// wk1[1]=data[1];
		wk2[n - 2] = data[n - 1];// wk2[n-1]=data[n];
		for (j = 2; j <= n - 1; j++) {
			wk1[j - 1] = data[j - 1];// wk1[j]=data[j];
			wk2[j - 2] = data[j - 1];// wk2[j-1]=data[j];
		}
		for (k = 1; k <= m; k++) {
			double num = 0.0;
			double denom = 0.0;
			for (j = 1; j <= (n - k); j++) {
				num += wk1[j - 1] * wk2[j - 1];// num += wk1[j]*wk2[j];
				denom += (wk1[j - 1]) * (wk1[j - 1]) + (wk2[j - 1])
						* (wk2[j - 1]);// denom += SQR(wk1[j])+SQR(wk2[j]);
			}
			d[k - 1] = 2.0 * num / denom;// d[k]=2.0*num/denom;
			xms_memcof *= (1.0 - (d[k - 1]) * (d[k - 1]));// *xms *=
															// (1.0-SQR(d[k]));
			for (i = 1; i <= (k - 1); i++)
				d[i - 1] = wkm[i - 1] - d[k - 1] * wkm[k - i - 1];// d[i]=wkm[i]-d[k]*wkm[k-i];
			// The algorithm is recursive, building up the answer for larger and
			// larger values of m
			// until the desired value is reached. At this point in the
			// algorithm, one could return
			// the vector d and scalar xms for a set of LP coefficients with k
			// (rather than m)
			// terms.
			if (k == m) {
				// free_vector(wkm,1,m);
				// free_vector(wk2,1,n);
				// free_vector(wk1,1,n);
				return;
			}
			for (i = 1; i <= k; i++)
				wkm[i - 1] = d[i - 1];// wkm[i]=d[i];
			for (j = 1; j <= (n - k - 1); j++) {
				wk1[j - 1] -= wkm[k - 1] * wk2[j - 1];// wk1[j] -=
														// wkm[k]*wk2[j];
				wk2[j - 1] = wk2[j] - wkm[k - 1] * wk1[j];// wk2[j]=wk2[j+1]-wkm[k]*wk1[j+1];
			}
		}
		// nrerror("never get here in memcof.");
		failB = true;
		failS = "never get here in memcof.";
		return;
	}

	/*
	 * Here are procedures for rendering the LP coefficients stable (if you
	 * choose to do so), and for extrapolating a data set by linear prediction,
	 * using the original or massaged LP coefficients. The routine zroots (§9.5)
	 * is used to find all complex roots of a polynomial.
	 */

	// #define NMAX 100 Largest expected value of m.
	// #define ZERO Complex(0.0,0.0)
	// #define ONE Complex(1.0,0.0)
	/**
	 * Given the LP coefficients d[1..m], this routine finds all roots of the characteristic polynomial, 
	 * reflects any roots that are outside the unit circle back inside, and then returns a modified set of coefficients d[1..m].
	 * @param d d
	 * @param m m
	 */
	public static void fixrts(double[] d, int m)
	// Given the LP coefficients d[1..m], this routine finds all roots of the
	// characteristic polynomial
	// (13.6.14), reflects any roots that are outside the unit circle back
	// inside, and then returns a
	// modified set of coefficients d[1..m].
	{
		// void zroots(fcomplex a[], int m, fcomplex roots[], int
		// polish);RootFind
		int i = 0;
		int j = 0;
		int polish = 0;
		Complex[] a = new Complex[NMAX];
		Complex[] roots = new Complex[NMAX];
		a[m] = ONE;// a[m]=ONE;
		for (j = m - 1; j >= 0; j--)
			// Set up complex coefficients for polynomial root finder.
			a[j] = new Complex(-d[m - j - 1], 0.0);// a[j]=Complex(-d[m-j],0.0);
		polish = 1;
		RootFind rf = new RootFind();
		rf.zroots(a, m, roots, polish);// Find all the roots.
		for (j = 1; j <= m; j++)
			// Look for a...
			// if (Complex.Cabs(roots[j]) > 1.0) root outside the unit circle,
			if (Complex.Cabs(roots[j - 1]) > 1.0)// roots[1...m] in zroots!!!
				// roots[j]=Cdiv(ONE,Conjg(roots[j])); and reflect it back
				// inside.
				roots[j - 1] = Complex.Cdiv(ONE, Complex.Conjg(roots[j - 1]));
		// a[0]=Csub(ZERO,roots[1]); Now reconstruct the polynomial
		// coefficients,
		a[0] = Complex.Csub(ZERO, roots[0]);
		a[1] = ONE;
		for (j = 2; j <= m; j++) {// by looping over the roots
			a[j] = ONE;
			for (i = j; i >= 2; i--)
				// and synthetically multiplying.
				// a[i-1]=Csub(a[i-2],Cmul(roots[j],a[i-1]));
				a[i - 1] = Complex.Csub(a[i - 2],
						Complex.Cmul(roots[j - 1], a[i - 1]));
			// a[0]=Csub(ZERO,Cmul(roots[j],a[0]));
			a[0] = Complex.Csub(ZERO, Complex.Cmul(roots[j - 1], a[0]));
		}
		for (j = 0; j <= m - 1; j++)
			// The polynomial coefficients are guaranteed to be
			// real, so we need only return the real part as
			// new LP coefficients.
			// d[m-j] = -a[j].r;
			d[m - j - 1] = -a[j].r;
	}

	/**
	 * Given data[1..ndata], and given the data’s LP coefficients d[1..m], this routine 
	 * predict the next nfut data points, which it returns in the array future[1..nfut]. 
	 * Note that the routine references only the last m values of data, as initial values for the prediction.
	 * @param data data
	 * @param ndata ndata
	 * @param d d
	 * @param m m
	 * @param future future
	 * @param nfut nfut
	 */
	public static void predic(double[] data, int ndata, double[] d, int m,
			double[] future, int nfut)
	// Given data[1..ndata], and given the data’s LP coefficients d[1..m], this
	// routine applies
	// equation (13.6.11) to predict the next nfut data points, which it returns
	// in the array
	// future[1..nfut]. Note that the routine references only the last m values
	// of data, as initial
	// values for the prediction.
	{
		int k = 0;
		int j = 0;
		double sum = 0.0;
		double discrp = 0.0;
		double[] reg = new double[m];
		// reg=vector(1,m);
		for (j = 1; j <= m; j++)
			reg[j - 1] = data[ndata - j];// reg[j]=data[ndata+1-j];
		for (j = 1; j <= nfut; j++) {
			discrp = 0.0;
			// This is where you would put in a known discrepancy if you were
			// reconstructing a
			// function by linear predictive coding rather than extrapolating a
			// function by linear prediction.
			// See text.
			sum = discrp;
			for (k = 1; k <= m; k++)
				sum += d[k - 1] * reg[k - 1];// sum += d[k]*reg[k];
			for (k = m; k >= 2; k--)
				reg[k - 1] = reg[k - 2];// reg[k]=reg[k-1]; [If you want to
										// implement circular
			// arrays, you can avoid this shifting of coefficients.]
			future[j - 1] = reg[0] = sum;// future[j]=reg[1]=sum;
		}
		// free_vector(reg,1,m);
	}

	/*
	 * Power Spectrum Estimation by the Maximum Entropy (All Poles) Method
	 */
	/**
	 * Given d[1..m], m, xms as returned by memcof, this function returns the power spectrum 
	 * estimate P(f) as a function of fdt = fDELTA.
	 * @param fdt fdt
	 * @param d d
	 * @param m m
	 * @param xms xms
	 * @return the result
	 */
	public static double evlmem(double fdt, double[] d, int m, float xms)
	// Given d[1..m], m, xms as returned by memcof, this function returns the
	// power spectrum
	// estimate P(f) as a function of fdt = f?.
	{
		int i = 0;
		double sumr = 1.0;
		double sumi = 0.0;
		double wr = 1.0;
		double wi = 0.0;
		double wpr = 0.0;
		double wpi = 0.0;
		double wtemp = 0.0;
		double theta = 0.0;// Trig. recurrences in double precision.
		theta = 6.28318530717959 * fdt;
		wpr = Math.cos(theta); // Set up for recurrence relations.
		wpi = Math.sin(theta);
		for (i = 1; i <= m; i++) {// Loop over the terms in the sum.
			wr = (wtemp = wr) * wpr - wi * wpi;
			wi = wi * wpr + wtemp * wpi;
			sumr -= d[i - 1] * wr;// sumr -= d[i]*wr; These accumulate the
									// denominator of (13.7.4).
			sumi -= d[i - 1] * wi;// sumi -= d[i]*wi;
		}
		return xms / (sumr * sumr + sumi * sumi);// Equation (13.7.4).
	}

	/*
	 * Spectral Analysis of Unevenly Sampled Data
	 */
	// public static void avevar(double data[], unsigned long n, float *ave,
	// float *var)
	/**
	 * Given array data[1..n], returns its mean as ave_avevar and its variance as var_avevar.
	 * @param data data
	 * @param n n
	 */
	public static void avevar(double data[], int n)// , float *ave, float *var)
	// Given array data[1..n], returns its mean as ave and its variance as var.
	{
		int j = 0;// unsigned long j;
		double s = 0.0;
		double ep = 0.0;
		for (ave_avevar = 0.0, j = 1; j <= n; j++)
			ave_avevar += data[j - 1];// *ave += data[j];
		ave_avevar /= n;// *ave /= n;
		var_avevar = ep = 0.0;
		for (j = 1; j <= n; j++) {
			s = data[j - 1] - (ave_avevar);// s=data[j]-(*ave);
			ep += s;
			var_avevar += s * s;
		}
		// *var=(*var-ep*ep/n)/(n-1); Corrected two-pass formula (14.1.8).
		var_avevar = (var_avevar - ep * ep / n) / (n - 1);
	}

	/**
	 * Given n data points with abscissas x[1..n] (which need not be equally spaced) and ordinates 
	 * y[1..n], and given a desired oversampling factor ofac (a typical value being 4 or larger), 
	 * this routine fills array px[1..np] with an increasing sequence of frequencies (not angular 
	 * frequencies) up to hifac times the “average” Nyquist frequency, and fills array py[1..np] 
	 * with the values of the Lomb normalized periodogram at those frequencies. 
	 * The arrays x and y are not altered. np, the dimension of px and py, must be large enough to 
	 * contain the output, or an error results. The routine also returns jmax such that py[jmax] is 
	 * the maximum element in py, and prob, an estimate of the significance of that maximum against 
	 * the hypothesis of random noise. A small value of prob indicates that a significant periodic 
	 * signal is present.
	 * @param x x
	 * @param y y
	 * @param n n
	 * @param ofac ofac
	 * @param hifac hifac
	 * @param px px
	 * @param py py
	 * @param np np
	 */
	public static void period(double[] x, double[] y, int n, double ofac,
			double hifac, double[] px, double[] py, int np)// , int *nout, int
															// *jmax, float
															// *prob)
	// Given n data points with abscissas x[1..n] (which need not be equally
	// spaced) and ordinates
	// y[1..n], and given a desired oversampling factor ofac (a typical value
	// being 4 or larger),
	// this routine fills array px[1..np] with an increasing sequence of
	// frequencies (not angular
	// frequencies) up to hifac times the “average” Nyquist frequency, and fills
	// array py[1..np]
	// with the values of the Lomb normalized periodogram at those frequencies.
	// The arrays x and y
	// are not altered. np, the dimension of px and py, must be large enough to
	// contain the output,
	// or an error results. The routine also returns jmax such that py[jmax] is
	// the maximum element
	// in py, and prob, an estimate of the significance of that maximum against
	// the hypothesis of
	// random noise. A small value of prob indicates that a significant periodic
	// signal is present.
	{
		// void avevar(float data[], unsigned long n, float *ave, float *var);
		int i = 0;
		int j = 0;
		double ave = 0.0;
		double c = 0.0;
		double cc = 0.0;
		double cwtau = 0.0;
		double effm = 0.0;
		double expy = 0.0;
		double pnow = 0.0;
		double pymax = 0.0;
		double s = 0.0;
		double ss = 0.0;
		double sumc = 0.0;
		double sumcy = 0.0;
		double sums = 0.0;
		double sumsh = 0.0;
		double sumsy = 0.0;
		double swtau = 0.0;
		double var = 0.0;
		double wtau = 0.0;
		double xave = 0.0;
		double xdif = 0.0;
		double xmax = 0.0;
		double xmin = 0.0;
		double yy = 0.0;
		double arg = 0.0;
		double wtemp = 0.0;

		double[] wi = new double[n];
		double[] wpi = new double[n];
		double[] wpr = new double[n];
		double[] wr = new double[n];
		// wi=dvector(1,n);
		// wpi=dvector(1,n);
		// wpr=dvector(1,n);
		// wr=dvector(1,n);
		failB = false;
		// nout_period=0.5*ofac*hifac*n;
		nout_period = (new Double(0.5 * ofac * hifac * n)).intValue();
		if (nout_period > np) {
			failB = true;
			failS = "output arrays too short in period";
			return;
			// nrerror("output arrays too short in period");
		}
		avevar(y, n);// ,&ave,&var); Get mean and variance of the input data.
		ave = ave_avevar;
		var = var_avevar;
		if (var == 0.0) {
			failB = true;
			failS = "zero variance in period";
			return;
			// nrerror("zero variance in period");
		}
		xmax = xmin = x[0];// xmax=xmin=x[1]; Go through data to get the range
							// of abscissas.
		for (j = 1; j <= n; j++) {
			if (x[j - 1] > xmax)
				xmax = x[j - 1];// if (x[j] > xmax) xmax=x[j];
			if (x[j - 1] < xmin)
				xmin = x[j - 1];// if (x[j] < xmin) xmin=x[j];
		}
		xdif = xmax - xmin;
		xave = 0.5 * (xmax + xmin);
		pymax = 0.0;
		pnow = 1.0 / (xdif * ofac);// Starting frequency.
		for (j = 1; j <= n; j++) {// Initialize values for the trigonometric
									// recurrences
									// at each data point. The recurrences are
									// done in double precision.
			arg = TWOPID * ((x[j - 1] - xave) * pnow);// arg=TWOPID*((x[j]-xave)*pnow);
			wpr[j - 1] = -2.0 * (Math.sin(0.5 * arg)) * (Math.sin(0.5 * arg));// wpr[j]
																				// =
																				// -2.0*SQR(sin(0.5*arg));
			wpi[j - 1] = Math.sin(arg);// wpi[j]=sin(arg);
			wr[j - 1] = Math.cos(arg);// wr[j]=cos(arg);
			wi[j - 1] = wpi[j - 1];// wi[j]=wpi[j];
		}
		for (i = 1; i <= (nout_period); i++) {// Main loop over the frequencies
												// to be evaluated.
			px[i - 1] = pnow;// px[i]=pnow;
			sumsh = sumc = 0.0;// First, loop over the data to get tau and
								// related quantities.
			for (j = 1; j <= n; j++) {
				c = wr[j - 1];// c=wr[j];
				s = wi[j - 1];// s=wi[j];
				sumsh += s * c;
				sumc += (c - s) * (c + s);
			}
			wtau = 0.5 * Math.atan2(2.0 * sumsh, sumc);
			/*
			 * double atan2 ( double y, double x ); Calculate arctangent, 2
			 * parameters. Performs the trigonometric arctangent operation on
			 * y/x and returns an angle in the range from -PI to PI expressed in
			 * radians, using the signs of the parameters to determine the
			 * quadrant. The result is valid even if x is 0 (angle is PI/2 or
			 * -PI/2). In fact this function returns the angle of bidimensional
			 * vector (x,y). returns:Arctangent of y/x.
			 * 
			 * Java: public static double atan2(double y, double x) Converts
			 * rectangular coordinates (x, y) to polar (r, theta). This method
			 * computes the phase theta by computing an arc tangent of y/x in
			 * the range of -pi to pi.
			 */
			swtau = Math.sin(wtau);
			cwtau = Math.cos(wtau);
			sums = sumc = sumsy = sumcy = 0.0;// Then, loop over the data again
												// to get the periodogram value.
			for (j = 1; j <= n; j++) {
				s = wi[j - 1];// s=wi[j];
				c = wr[j - 1];// c=wr[j];
				ss = s * cwtau - c * swtau;
				cc = c * cwtau + s * swtau;
				sums += ss * ss;
				sumc += cc * cc;
				yy = y[j - 1] - ave;// yy=y[j]-ave;
				sumsy += yy * ss;
				sumcy += yy * cc;
				// wr[j]=((wtemp=wr[j])*wpr[j]-wi[j]*wpi[j])+wr[j]; Update the
				// trigonometric recurrences.
				wr[j - 1] = ((wtemp = wr[j - 1]) * wpr[j - 1] - wi[j - 1]
						* wpi[j - 1])
						+ wr[j - 1];
				// wi[j]=(wi[j]*wpr[j]+wtemp*wpi[j])+wi[j];
				wi[j - 1] = (wi[j - 1] * wpr[j - 1] + wtemp * wpi[j - 1])
						+ wi[j - 1];
			}
			// py[i]=0.5*(sumcy*sumcy/sumc+sumsy*sumsy/sums)/var;
			py[i - 1] = 0.5 * (sumcy * sumcy / sumc + sumsy * sumsy / sums)
					/ var;
			// if (py[i] >= pymax) pymax=py[(*jmax=i)];
			if (py[i - 1] >= pymax)
				pymax = py[(jmax_period = i) - 1];
			pnow += 1.0 / (ofac * xdif);// The next frequency.
		}
		expy = Math.exp(-pymax);// Evaluate statistical significance of the
								// maximum.
		effm = 2.0 * (nout_period) / ofac;
		prob_period = effm * expy;
		if (prob_period > 0.01)
			prob_period = 1.0 - Math.pow(1.0 - expy, effm);
		// free_dvector(wr,1,n);
		// free_dvector(wpr,1,n);
		// free_dvector(wpi,1,n);
		// free_dvector(wi,1,n);
	}

	// public static void fasper(double[] x, double[] y, unsigned long n, float
	// ofac, float hifac,
	// float wk1[], float wk2[], unsigned long nwk, unsigned long *nout,
	// unsigned long *jmax, float *prob)
	/**
	 * Given n data points with abscissas x[1..n] (which need not be equally spaced) and ordinates 
	 * y[1..n], and given a desired oversampling factor ofac (a typical value being 4 or larger), this 
	 * routine fills array wk1[1..nwk] with a sequence of nout_fasper increasing frequencies (not angular 
	 * frequencies) up to hifac times the “average” Nyquist frequency, and fills array wk2[1..nwk] 
	 * with the values of the Lomb normalized periodogram at those frequencies. 
	 * The arrays x and y are not altered. nwk, the dimension of wk1 and wk2, must be large enough for intermediate 
	 * work space, or an error results. The routine also returns jmax such that wk2[jmax] is the 
	 * maximum element in wk2, and prob_fasper, an estimate of the significance of that maximum against 
	 * the hypothesis of random noise. A small value of prob_fasper indicates that a significant periodic signal is present.
	 * @param x x
	 * @param y y
	 * @param n n
	 * @param ofac ofac
	 * @param hifac hifac
	 * @param wk1 wk1
	 * @param wk2 wk2
	 * @param nwk nwk
	 */
	public static void fasper(double[] x, double[] y, int n, double ofac,
			double hifac, double[] wk1, double[] wk2, int nwk)// , unsigned long
																// *nout,
	// unsigned long *jmax, float *prob)
	// Given n data points with abscissas x[1..n] (which need not be equally
	// spaced) and ordinates
	// y[1..n], and given a desired oversampling factor ofac (a typical value
	// being 4o r larger), this
	// routine fills array wk1[1..nwk] with a sequence of nout increasing
	// frequencies (not angular
	// frequencies) up to hifac times the “average” Nyquist frequency, and fills
	// array wk2[1..nwk]
	// with the values of the Lomb normalized periodogram at those frequencies.
	// The arrays x and
	// y are not altered. nwk, the dimension of wk1 and wk2, must be large
	// enough for intermediate
	// work space, or an error results. The routine also returns jmax such that
	// wk2[jmax] is the
	// maximum element in wk2, and prob, an estimate of the significance of that
	// maximum against
	// the hypothesis of random noise. A small value of prob indicates that a
	// significant periodic
	// signal is present.
	{
		// void avevar(float data[], unsigned long n, float *ave, float *var);
		// void realft(float data[], unsigned long n, int isign);
		// void spread(float y, float yy[], unsigned long n, float x, int m);
		// unsigned long j,k,ndim,nfreq,nfreqt;
		int j = 0;
		int k = 0;
		int ndim = 0;
		int nfreq = 0;
		int nfreqt = 0;
		double ave = 0.0;
		double ck = 0.0;
		double ckk = 0.0;
		double cterm = 0.0;
		double cwt = 0.0;
		double den = 0.0;
		double df = 0.0;
		double effm = 0.0;
		double expy = 0.0;
		double fac = 0.0;
		double fndim = 0.0;
		double hc2wt = 0.0;
		double hs2wt = 0.0;
		double hypo = 0.0;
		double pmax = 0.0;
		double sterm = 0.0;
		double swt = 0.0;
		double var = 0.0;
		double xdif = 0.0;
		double xmax = 0.0;
		double xmin = 0.0;

		failB = false;

		// nout_fasper=0.5*ofac*hifac*n;
		nout_fasper = (new Double(0.5 * ofac * hifac * n)).intValue();
		// nfreqt=ofac*hifac*n*MACC;// Size the FFT as next power of 2 above
		// nfreqt.
		nfreqt = (new Double(ofac * hifac * n * MACC)).intValue();
		nfreq = 64;
		while (nfreq < nfreqt)
			nfreq <<= 1;
		ndim = nfreq << 1;
		if (ndim > nwk) {
			failB = true;
			failS = "workspaces too small in fasper";
			return;
			// nrerror("workspaces too small in fasper");
		}
		// avevar(y,n,&ave,&var); Compute the mean, variance, and range of the
		// data.
		avevar(y, n);
		ave = ave_avevar;
		var = var_avevar;
		if (var == 0.0) {
			// nrerror("zero variance in fasper");
			failB = true;
			failS = "zero variance in fasper";
			return;
		}
		xmin = x[0];// xmin=x[1];
		xmax = xmin;
		for (j = 2; j <= n; j++) {
			if (x[j - 1] < xmin)
				xmin = x[j - 1];// if (x[j] < xmin) xmin=x[j];
			if (x[j - 1] > xmax)
				xmax = x[j - 1];// if (x[j] > xmax) xmax=x[j];
		}
		xdif = xmax - xmin;
		for (j = 1; j <= ndim; j++)
			wk1[j - 1] = wk2[j - 1] = 0.0;// wk1[j]=wk2[j]=0.0; Zero the
											// workspaces.
		fac = ndim / (xdif * ofac);
		fndim = ndim;
		for (j = 1; j <= n; j++) {// Extirpolate the data into the workspaces.
			ck = (x[j - 1] - xmin) * fac;// ck=(x[j]-xmin)*fac;
			// #define MOD(a,b) while(a >= b) a -= b;
			// MOD(ck,fndim)
			while (ck >= fndim)
				ck -= fndim;
			ckk = 2.0 * (ck++);
			// MOD(ckk,fndim)
			while (ckk >= fndim)
				ckk -= fndim;
			++ckk;
			// spread(y[j]-ave,wk1,ndim,ck,MACC);
			spread(y[j - 1] - ave, wk1, ndim, ck, MACC);
			spread(1.0, wk2, ndim, ckk, MACC);
		}
		FastFourierTransform.realft(wk1, ndim, 1); // Take the Fast Fourier
													// Transforms.
		FastFourierTransform.realft(wk2, ndim, 1);
		df = 1.0 / (xdif * ofac);
		pmax = -1.0;
		for (k = 3, j = 1; j <= (nout_fasper); j++, k += 2) {// Compute the Lomb
																// value for
																// each
																// frequency.
																// hypo=sqrt(wk2[k]*wk2[k]+wk2[k+1]*wk2[k+1]);
			hypo = Math.sqrt(wk2[k - 1] * wk2[k - 1] + wk2[k] * wk2[k]);
			hc2wt = 0.5 * wk2[k - 1] / hypo;// hc2wt=0.5*wk2[k]/hypo;
			hs2wt = 0.5 * wk2[k] / hypo;// hs2wt=0.5*wk2[k+1]/hypo;
			cwt = Math.sqrt(0.5 + hc2wt);
			swt = SIGN(Math.sqrt(0.5 - hc2wt), hs2wt);
			// den=0.5*n+hc2wt*wk2[k]+hs2wt*wk2[k+1];
			den = 0.5 * n + hc2wt * wk2[k - 1] + hs2wt * wk2[k];
			// cterm=SQR(cwt*wk1[k]+swt*wk1[k+1])/den;
			cterm = (cwt * wk1[k - 1] + swt * wk1[k])
					* (cwt * wk1[k - 1] + swt * wk1[k]) / den;
			// sterm=SQR(cwt*wk1[k+1]-swt*wk1[k])/(n-den);
			sterm = (cwt * wk1[k] - swt * wk1[k - 1])
					* (cwt * wk1[k] - swt * wk1[k - 1]) / (n - den);
			wk1[j - 1] = j * df;// wk1[j]=j*df;
			wk2[j - 1] = (cterm + sterm) / (2.0 * var);// wk2[j]=(cterm+sterm)/(2.0*var);
			if (wk2[j] > pmax)
				pmax = wk2[(jmax_fasper = j) - 1];// pmax=wk2[(*jmax=j)];
		}
		expy = Math.exp(-pmax);// Estimate significance of largest peak value.
		effm = 2.0 * (nout_fasper) / ofac;
		prob_fasper = effm * expy;
		if (prob_fasper > 0.01)
			prob_fasper = 1.0 - Math.pow(1.0 - expy, effm);
	}

	// public static void spread(double y, double[] yy, unsigned long n, float
	// x, int m)
	/**
	 * Given an array yy[1..n], extirpolate (spread) a value y into m actual array elements that best 
	 * approximate the “fictional” (i.e., possibly noninteger) array element number x. The weights 
	 * used are coefficients of the Lagrange interpolating polynomial.
	 * @param y y
	 * @param yy yy
	 * @param n n
	 * @param x x
	 * @param m m
	 */
	public static void spread(double y, double[] yy, int n, double x, int m)
	// Given an array yy[1..n], extirpolate (spread) a value y into m actual
	// array elements that best
	// approximate the “fictional” (i.e., possibly noninteger) array element
	// number x. The weights
	// used are coefficients of the Lagrange interpolating polynomial.
	{
		int ihi = 0;
		int ilo = 0;
		int ix = 0;
		int j = 0;
		int nden = 0;
		// static long nfac[11]={0,1,1,2,6,24,120,720,5040,40320,362880};
		int[] nfac = { 0, 1, 1, 2, 6, 24, 120, 720, 5040, 40320, 362880 };
		double fac = 0.0;

		failB = false;
		if (m > 10) {
			failB = true;
			failS = "factorial table too small in spread";
			return;
			// nrerror("factorial table too small in spread");
		}
		ix = (new Double(x)).intValue();// (int)x;
		// if (x == (float)ix) yy[ix] += y;
		if (x == (new Integer(ix)).doubleValue())
			yy[ix - 1] += y;// yy[ix] += y;
		else {
			// ilo=LMIN(LMAX((long)(x-0.5*m+1.0),1),n-m+1);
			ilo = Math.min(
					Math.max((new Double(x - 0.5 * m + 1.0)).intValue(), 1), n
							- m + 1);
			ihi = ilo + m - 1;
			nden = nfac[m];// nden=nfac[m];
			fac = x - ilo;
			for (j = ilo + 1; j <= ihi; j++)
				fac *= (x - j);
			yy[ihi - 1] += y * fac / (nden * (x - ihi));// yy[ihi] +=
														// y*fac/(nden*(x-ihi));
			for (j = ihi - 1; j >= ilo; j--) {
				nden = (nden / (j + 1 - ilo)) * (j - ihi);
				// yy[j] += y*fac/(nden*(x-j));
				yy[j - 1] += y * fac / (nden * (x - j));
			}
		}
	}

	/*
	 * Computing Fourier Integrals Using the FFT
	 */
	/**
	 * For an integral approximated by a discrete Fourier transform, this routine computes the correction 
	 * factor that multiplies the DFT and the endpoint correction to be added. Input is the 
	 * angular frequency w, stepsize delta, lower and upper limits of the integral a and b, while the 
	 * array endpts contains the first 4 and last 4 function values. The correction factor W(THETA) is 
	 * returned as corfac_dftcor, while the real and imaginary parts of the endpoint correction are returned 
	 * as corre_dftcor and corim_dftcor.
	 * @param w w
	 * @param delta delta
	 * @param a a
	 * @param b b
	 * @param endpts endpts
	 */
	public static void dftcor(double w, double delta, double a, double b,
			double endpts[])// ,
	// float *corre, float *corim, float *corfac)
	// For an integral approximated by a discrete Fourier transform, this
	// routine computes the correction
	// factor that multiplies the DFT and the endpoint correction to be added.
	// Input is the
	// angular frequency w, stepsize delta, lower and upper limits of the
	// integral a and b, while the
	// array endpts contains the first 4 and last 4 function values. The
	// correction factor W(?) is
	// returned as corfac, while the real and imaginary parts of the endpoint
	// correction are returned
	// as corre and corim.
	{
		// void nrerror(char error_text[]);
		double a0i = 0.0;
		double a0r = 0.0;
		double a1i = 0.0;
		double a1r = 0.0;
		double a2i = 0.0;
		double a2r = 0.0;
		double a3i = 0.0;
		double a3r = 0.0;
		double arg = 0.0;
		double c = 0.0;
		double cl = 0.0;
		double cr = 0.0;
		double s = 0.0;
		double sl = 0.0;
		double sr = 0.0;
		double t = 0.0;
		double t2 = 0.0;
		double t4 = 0.0;
		double t6 = 0.0;
		double cth = 0.0;
		double ctth = 0.0;
		double spth2 = 0.0;
		double sth = 0.0;
		double sth4i = 0.0;
		double stth = 0.0;
		double th = 0.0;
		double th2 = 0.0;
		double th4 = 0.0;
		double tmth2 = 0.0;
		double tth4i = 0.0;
		failB = false;

		th = w * delta;
		if (a >= b || th < 0.0e0 || th > 3.1416e0) {
			failB = true;
			failS = "bad arguments to dftcor";
			return;
			// nrerror("bad arguments to dftcor");
		}
		if (Math.abs(th) < 5.0e-2) {// Use series.
			t = th;
			t2 = t * t;
			t4 = t2 * t2;
			t6 = t4 * t2;
			corfac_dftcor = 1.0 - (11.0 / 720.0) * t4 + (23.0 / 15120.0) * t6;
			a0r = (-2.0 / 3.0) + t2 / 45.0 + (103.0 / 15120.0) * t4
					- (169.0 / 226800.0) * t6;
			a1r = (7.0 / 24.0) - (7.0 / 180.0) * t2 + (5.0 / 3456.0) * t4
					- (7.0 / 259200.0) * t6;
			a2r = (-1.0 / 6.0) + t2 / 45.0 - (5.0 / 6048.0) * t4 + t6 / 64800.0;
			a3r = (1.0 / 24.0) - t2 / 180.0 + (5.0 / 24192.0) * t4 - t6
					/ 259200.0;
			a0i = t
					* (2.0 / 45.0 + (2.0 / 105.0) * t2 - (8.0 / 2835.0) * t4 + (86.0 / 467775.0)
							* t6);
			a1i = t
					* (7.0 / 72.0 - t2 / 168.0 + (11.0 / 72576.0) * t4 - (13.0 / 5987520.0)
							* t6);
			a2i = t
					* (-7.0 / 90.0 + t2 / 210.0 - (11.0 / 90720.0) * t4 + (13.0 / 7484400.0)
							* t6);
			a3i = t
					* (7.0 / 360.0 - t2 / 840.0 + (11.0 / 362880.0) * t4 - (13.0 / 29937600.0)
							* t6);
		} else {// Use trigonometric formulas in double precision.
			cth = Math.cos(th);
			sth = Math.sin(th);
			ctth = cth * cth - sth * sth;
			stth = 2.0e0 * sth * cth;
			th2 = th * th;
			th4 = th2 * th2;
			tmth2 = 3.0e0 - th2;
			spth2 = 6.0e0 + th2;
			sth4i = 1.0 / (6.0e0 * th4);
			tth4i = 2.0e0 * sth4i;
			corfac_dftcor = tth4i * spth2 * (3.0e0 - 4.0e0 * cth + ctth);
			a0r = sth4i
					* (-42.0e0 + 5.0e0 * th2 + spth2 * (8.0e0 * cth - ctth));
			a0i = sth4i * (th * (-12.0e0 + 6.0e0 * th2) + spth2 * stth);
			a1r = sth4i * (14.0e0 * tmth2 - 7.0e0 * spth2 * cth);
			a1i = sth4i * (30.0e0 * th - 5.0e0 * spth2 * sth);
			a2r = tth4i * (-4.0e0 * tmth2 + 2.0e0 * spth2 * cth);
			a2i = tth4i * (-12.0e0 * th + 2.0e0 * spth2 * sth);
			a3r = sth4i * (2.0e0 * tmth2 - spth2 * cth);
			a3i = sth4i * (6.0e0 * th - spth2 * sth);
		}
		// cl=a0r*endpts[1]+a1r*endpts[2]+a2r*endpts[3]+a3r*endpts[4];
		cl = a0r * endpts[0] + a1r * endpts[1] + a2r * endpts[2] + a3r
				* endpts[3];
		// sl=a0i*endpts[1]+a1i*endpts[2]+a2i*endpts[3]+a3i*endpts[4];
		sl = a0i * endpts[0] + a1i * endpts[1] + a2i * endpts[2] + a3i
				* endpts[3];
		// cr=a0r*endpts[8]+a1r*endpts[7]+a2r*endpts[6]+a3r*endpts[5];
		cr = a0r * endpts[7] + a1r * endpts[6] + a2r * endpts[5] + a3r
				* endpts[4];
		// sr = -a0i*endpts[8]-a1i*endpts[7]-a2i*endpts[6]-a3i*endpts[5];
		sr = -a0i * endpts[7] - a1i * endpts[6] - a2i * endpts[5] - a3i
				* endpts[4];
		arg = w * (b - a);
		c = Math.cos(arg);
		s = Math.sin(arg);
		corre_dftcor = cl + c * cr - s * sr;
		corim_dftcor = sl + s * cr + c * sr;
	}

	// public void dftint(float (*func)(float), float a, float b, float w, float
	// *cosint,
	// float *sinint)float *cosint,
	/**
	 * Example program illustrating how to use the routine dftcor. The user supplies an external 
	 * function func that returns the quantity h(t). The routine then returns integral a-b of cos(omega*t)h(t) dt as 
	 * cosint and integral a-b of sin(omegat)h(t) dt as sinint.
	 * @param a a
	 * @param b b
	 * @param w w
	 */
	public void dftint(double a, double b, double w)// , float *cosint,
	// Example program illustrating how to use the routine dftcor. The user
	// supplies an external
	// function func that returns the quantity h(t). The routine then returns
	// integral a-b
	// of cos(omega*t)h(t) dt as
	// cosint and integral a-b of sin(?t)h(t) dt as sinint.
	{
		// void dftcor(float w, float delta, float a, float b, float endpts[],
		// float *corre, float *corim, float *corfac);
		// void polint(float xa[], float ya[], int n, float x, float *y, float
		// *dy);
		// void realft(float data[], unsigned long n, int isign);
		// static int init=0;
		int j = 0;
		int nn = 0;
		// static float aold = -1.e30,bold = -1.e30,delta,(*funcold)(float);
		// static float data[NDFT+1],endpts[9];
		double c = 0.0;
		double cdft = 0.0;
		// double cerr = 0.0;
		double corfac = 0.0;
		double corim = 0.0;
		double corre = 0.0;
		double en = 0.0;
		double s = 0.0;
		double sdft = 0.0;
		// double serr = 0.0;
		double[] cpol = new double[MPOL];
		double[] spol = new double[MPOL];
		double[] xpol = new double[MPOL];
		// cpol=vector(1,MPOL);
		// spol=vector(1,MPOL);
		// xpol=vector(1,MPOL);
		if (init != 1 || a != aold || b != bold)// || func != funcold)
		{
			// Do we need to initialize, or is only ? changed?
			init = 1;
			aold = a;
			bold = b;
			// funcold=func;
			delta = (b - a) / M;
			// Load the function values into the data array.
			for (j = 1; j <= M + 1; j++)
				data[j - 1] = func.F(a + (j - 1) * delta);// data[j]=see
															// def!!!BUT..for
															// consistency with
															// realft
			for (j = M + 2; j <= NDFT; j++)
				// Zero pad the rest of the data array.
				data[j - 1] = 0.0;// [j]
			for (j = 1; j <= 4; j++) {// Load the endpoints.
				endpts[j - 1] = data[j - 1];// endpts[j]=data[j];
				endpts[j + 4 - 1] = data[M - 3 + j - 1];// endpts[j+4]=data[M-3+j];
			}
			FastFourierTransform.realft(data, NDFT, 1);
			// realft returns the unused value corresponding to ?N/2 in data[2].
			// We actually want
			// this element to contain the imaginary part corresponding to ?0,
			// which is zero.
			data[1] = 0.0;// data[2]=0.0;
		}
		// Now interpolate on the DFT result for the desired frequency. If the
		// frequency is an ?n,
		// i.e., the quantity en is an integer, then cdft=data[2*en-1],
		// sdft=data[2*en], and you
		// could omit the interpolation.
		en = w * delta * NDFT / TWOPI + 1.0;
		nn = Math.min(Math.max((int) (en - 0.5 * MPOL + 1.0), 1), NDFT / 2
				- MPOL + 1);// Leftmost point for the
		// interpolation.
		for (j = 1; j <= MPOL; j++, nn++) {
			cpol[j - 1] = data[2 * nn - 2];// cpol[j]=data[2*nn-1];
			spol[j - 1] = data[2 * nn - 1];// spol[j]=data[2*nn];
			xpol[j - 1] = nn;// xpol[j]=nn;
		}
		Interpolator.polint(xpol, cpol, MPOL, en);// ,&cdft,&cerr);
		cdft = Interpolator.ypoli;
		// cerr = Interpolator.dypoli;
		Interpolator.polint(xpol, spol, MPOL, en);// ,&sdft,&serr);
		sdft = Interpolator.ypoli;
		// serr = Interpolator.dypoli;
		// dftcor(w,delta,a,b,endpts,&corre,&corim,&corfac); Now get the
		// endpoint correction
		// and the multiplicative factor W(theta).
		dftcor(w, delta, a, b, endpts);//
		corre = corre_dftcor;
		corim = corim_dftcor;
		corfac = corfac_dftcor;
		cdft *= corfac;
		sdft *= corfac;
		cdft += corre;
		sdft += corim;
		c = delta * Math.cos(w * a);// Finally multiply by ? and exp(i?a).
		s = delta * Math.sin(w * a);
		cosint_dftint = c * cdft - s * sdft;
		sinint_dftint = s * cdft + c * sdft;
		// free_vector(cpol,1,MPOL);
		// free_vector(spol,1,MPOL);
		// free_vector(xpol,1,MPOL);
	}

	/*
	 * Wavelet Transforms Like the fast Fourier transform (FFT), the discrete
	 * wavelet transform (DWT) is a fast, linear operation that operates on a
	 * data vector whose length is an integer power of two, transforming it into
	 * a numerically different vector of the same length. Also like the FFT, the
	 * wavelet transform is invertible and in fact orthogonal—the inverse
	 * transform, when viewed as a big matrix, is simply the transpose of the
	 * transform. Both FFT and DWT, therefore, can be viewed as a rotation in
	 * function space, from the input space (or time) domain, where the basis
	 * functions are the unit vectors e i, or Dirac delta functions in the
	 * continuum limit, to a different domain. For the FFT, this new domain has
	 * basis functions that are the familiar sines and cosines. In the wavelet
	 * domain, the basis functions are somewhat more complicated and have the
	 * fanciful names “mother functions” and “wavelets.” Of course there are an
	 * infinity of possible bases for function space, almost all of them
	 * uninteresting! What makes the wavelet basis interesting is that, unlike
	 * sines and cosines, individual wavelet functions are quite localized in
	 * space; simultaneously, like sines and cosines, individual wavelet
	 * functions are quite localized in frequency or (more precisely)
	 * characteristic scale. As we will see below, the particular kind of dual
	 * localization achieved by wavelets renders large classes of functions and
	 * operators sparse, or sparse to some high accuracy, when transformed into
	 * the wavelet domain. Analogously with the Fourier domain, where a class of
	 * computations, like convolutions, become computationally fast, there is a
	 * large class of computations
	 */
	/**
	 * Calling daub4 if namS is set to "daub4";
	 * @param a a
	 * @param n n
	 * @param isign isign
	 * @param namS namS
	 */
	public static void wtstep(double a[], int n, int isign, String namS) {
		if (namS.compareTo("daub4") == 0) {
			daub4(a, n, isign);
		}
		// else if ()

	}

	/**
	 * One-dimensional discrete wavelet transform. This routine implements the pyramid algorithm, 
	 * replacing a[1..n] by its wavelet transform (for isign=1), or performing the inverse operation 
	 * (for isign=-1). Note that n MUST be an integer power of 2. The routine wtstep, whose 
	 * actual name must be supplied in calling this routine, is the underlying wavelet filter. Examples 
	 * of wtstep are daub4 and (preceded by pwtset) pwt.
	 * @param a a
	 * @param n n
	 * @param isign isign
	 * @param namS namS
	 */
	public static void wt1(double[] a, int n, int isign, String namS)
	// void (*wtstep)(float [], unsigned long, int))
	// One-dimensional discrete wavelet transform. This routine implements the
	// pyramid algorithm,
	// replacing a[1..n] by its wavelet transform (for isign=1), or performing
	// the inverse operation
	// (for isign=-1). Note that n MUST be an integer power of 2. The routine
	// wtstep, whose
	// actual name must be supplied in calling this routine, is the underlying
	// wavelet filter. Examples
	// of wtstep are daub4 and (preceded by pwtset) pwt.
	{
		int nn = 0;// unsigned long nn;
		if (n < 4)
			return;
		if (isign >= 0) {// Wavelet transform.
			for (nn = n; nn >= 4; nn >>= 1)
				wtstep(a, nn, isign, "daub4");// (*wtstep)(a,nn,isign);
			// Start at largest hierarchy, and work towards smallest.
		} else {// Inverse wavelet transform.
			for (nn = 4; nn <= n; nn <<= 1)
				wtstep(a, nn, isign, "daub4");
			// Start at smallest hierarchy, and work towards largest.
		}
	}

	/**
	 * Applies the Daubechies 4-coefficient wavelet filter to data vector a[1..n] (for isign=1) or 
	 * applies its transpose (for isign=-1). Used hierarchically by routines wt1 and wtn.
	 * @param a a
	 * @param n n
	 * @param isign isign
	 */
	public static void daub4(double a[], int n, int isign)
	// Applies the Daubechies 4-coefficient wavelet filter to data vector
	// a[1..n] (for isign=1) or
	// applies its transpose (for isign=-1). Used hierarchically by routines wt1
	// and wtn.
	{
		double[] wksp = new double[n];
		int nh = 0;
		int nh1 = 0;
		int i = 0;
		int j = 0;
		if (n < 4)
			return;
		// wksp=vector(1,n);
		nh1 = (nh = n >> 1) + 1;
		if (isign >= 0) {// Apply filter.
			for (i = 1, j = 1; j <= n - 3; j += 2, i++) {
				// wksp[i]=C0*a[j]+C1*a[j+1]+C2*a[j+2]+C3*a[j+3];
				wksp[i - 1] = C0 * a[j - 1] + C1 * a[j] + C2 * a[j + 1] + C3
						* a[j + 2];
				// wksp[i+nh] = C3*a[j]-C2*a[j+1]+C1*a[j+2]-C0*a[j+3];
				wksp[i + nh - 1] = C3 * a[j - 1] - C2 * a[j] + C1 * a[j + 1]
						- C0 * a[j + 2];
			}
			// wksp[i]=C0*a[n-1]+C1*a[n]+C2*a[1]+C3*a[2];
			wksp[i - 1] = C0 * a[n - 2] + C1 * a[n - 1] + C2 * a[0] + C3 * a[1];
			// wksp[i+nh] = C3*a[n-1]-C2*a[n]+C1*a[1]-C0*a[2];
			wksp[i + nh - 1] = C3 * a[n - 2] - C2 * a[n - 1] + C1 * a[0] - C0
					* a[1];
		} else {// Apply transpose filter.
				// wksp[1]=C2*a[nh]+C1*a[n]+C0*a[1]+C3*a[nh1];
			wksp[0] = C2 * a[nh - 1] + C1 * a[n - 1] + C0 * a[0] + C3
					* a[nh1 - 1];
			// wksp[2] = C3*a[nh]-C0*a[n]+C1*a[1]-C2*a[nh1];
			wksp[1] = C3 * a[nh - 1] - C0 * a[n - 1] + C1 * a[0] - C2
					* a[nh1 - 1];
			for (i = 1, j = 3; i < nh; i++) {
				// wksp[j++]=C2*a[i]+C1*a[i+nh]+C0*a[i+1]+C3*a[i+nh1];
				wksp[j++ - 1] = C2 * a[i - 1] + C1 * a[i + nh - 1] + C0 * a[i]
						+ C3 * a[i + nh1 - 1];
				// wksp[j++] = C3*a[i]-C0*a[i+nh]+C1*a[i+1]-C2*a[i+nh1];
				wksp[j++ - 1] = C3 * a[i - 1] - C0 * a[i + nh - 1] + C1 * a[i]
						- C2 * a[i + nh1 - 1];
			}
		}
		for (i = 1; i <= n; i++)
			a[i - 1] = wksp[i - 1];// a[i]=wksp[i];
		// free_vector(wksp,1,n);
	}

	// typedef struct {
	// int ncof,ioff,joff;
	// float *cc,*cr;
	// } wavefilt;
	// wavefilt wfilt;

	/**
	 * Initializing routine for pwt, here implementing the Daubechies wavelet filters with 4, 12, and 
	 * 20 coefficients, as selected by the input value n. Further wavelet filters can be included in the 
	 * obvious manner. This routine must be called (once) before the first use of pwt. (For the case 
	 * n=4, the specific routine daub4 is considerably faster than pwt.)
	 * @param n n
	 */
	public static void pwtset(int n)
	// Initializing routine for pwt, here implementing the Daubechies wavelet
	// filters with 4, 12, and
	// 20 coefficients, as selected by the input value n. Further wavelet
	// filters can be included in the
	// obvious manner. This routine must be called (once) before the first use
	// of pwt. (For the case
	// n=4, the specific routine daub4 is considerably faster than pwt.)
	{
		// void nrerror(char error_text[]);
		failB = false;
		int k = 0;
		double sig = -1.0;
		double[] c4 = { 0.0, 0.4829629131445341, 0.8365163037378079,
				0.2241438680420134, -0.1294095225512604 };// [5]
		double[] c12 = { 0.0, 0.111540743350, 0.494623890398, 0.751133908021,
				0.315250351709, -0.226264693965, -0.129766867567,
				0.097501605587, 0.027522865530, -0.031582039318,
				0.000553842201, 0.004777257511, -0.001077301085 };// [13]
		double[] c20 = { 0.0, 0.026670057901, 0.188176800078, 0.527201188932,
				0.688459039454, 0.281172343661, -0.249846424327,
				-0.195946274377, 0.127369340336, 0.093057364604,
				-0.071394147166, -0.029457536822, 0.033212674059,
				0.003606553567, -0.010733175483, 0.001395351747,
				0.001992405295, -0.000685856695, -0.000116466855,
				0.000093588670, -0.000013264203 };// [21]

		// Wavefilt wfilt = new Wavefilt();

		Wavefilt.ncof = n;
		if (n == 4) {
			Wavefilt.cc = c4;
			Wavefilt.cr = c4r;
		} else if (n == 12) {
			Wavefilt.cc = c12;
			Wavefilt.cr = c12r;
		} else if (n == 20) {
			Wavefilt.cc = c20;
			Wavefilt.cr = c20r;
		} else {
			failB = true;
			failS = "unimplemented value n in pwtset";
			return;
			// nrerror("unimplemented value n in pwtset");
		}
		for (k = 1; k <= n; k++) {
			Wavefilt.cr[Wavefilt.ncof + 1 - k] = sig * Wavefilt.cc[k];// ok!!!
			sig = -sig;
		}
		Wavefilt.ioff = Wavefilt.joff = -(n >> 1);
		// These values center the “support” of the wavelets at each level.
		// Alternatively, the “peaks”
		// of the wavelets can be approximately centered by the choices ioff=-2
		// and joff=-n+2.
		// Note that daub4 and pwtset with n=4 use different default centerings.
	}

	// typedef struct {
	// int ncof,ioff,joff;
	// float *cc,*cr;
	// } wavefilt;
	// extern wavefilt wfilt; Defined in pwtset.
	/**
	 * Partial wavelet transform: applies an arbitrary wavelet filter to data vector a[1..n] (for isign = 
	 * 1) or applies its transpose (for isign = -1). Used hierarchically by routines wt1 and wtn. 
	 * The actual filter is determined by a preceding (and required) call to pwtset, which initializes the structure wfilt.
	 * @param a a
	 * @param n n
	 * @param isign isign
	 */
	public static void pwt(double a[], int n, int isign)
	// Partial wavelet transform: applies an arbitrary wavelet filter to data
	// vector a[1..n] (for isign =
	// 1) or applies its transpose (for isign = -1). Used hierarchically by
	// routines wt1 and wtn.
	// The actual filter is determined by a preceding (and required) call to
	// pwtset, which initializes
	// the structure wfilt.
	{
		double ai = 0.0;
		double ai1 = 0.0;
		double[] wksp = new double[n];
		int i = 0;
		int ii = 0;
		int j = 0;
		int jf = 0;
		int jr = 0;
		int k = 0;
		int n1 = 0;
		int ni = 0;
		int nj = 0;
		int nh = 0;
		int nmod = 0;
		if (n < 4)
			return;
		// wksp=vector(1,n);
		// Wavefilt wfilt = new Wavefilt();

		nmod = Wavefilt.ncof * n;// A positive constant equal to zero mod n.
		n1 = n - 1;// Mask of all bits, since n a power of 2.
		nh = n >> 1;
		for (j = 1; j <= n; j++)
			wksp[j - 1] = 0.0;// wksp[j]=0.0;
		if (isign >= 0) {// Apply filter.
			for (ii = 1, i = 1; i <= n; i += 2, ii++) {
				ni = i + nmod + Wavefilt.ioff;// Pointer to be incremented and
												// wrapped-around.
				nj = i + nmod + Wavefilt.joff;
				for (k = 1; k <= Wavefilt.ncof; k++) {
					jf = n1 & (ni + k);// We use bitwise and to wrap-around the
										// pointers.
					jr = n1 & (nj + k);
					// wksp[ii] += wfilt.cc[k]*a[jf+1];
					wksp[ii - 1] += Wavefilt.cc[k] * a[jf];
					// wksp[ii+nh] += wfilt.cr[k]*a[jr+1];
					wksp[ii + nh - 1] += Wavefilt.cr[k] * a[jr];
				}
			}
		} else {// Apply transpose filter.
			for (ii = 1, i = 1; i <= n; i += 2, ii++) {
				ai = a[ii - 1];// ai=a[ii];
				ai1 = a[ii + nh - 1];// ai1=a[ii+nh];
				ni = i + nmod + Wavefilt.ioff; // See comments above.
				nj = i + nmod + Wavefilt.joff;
				for (k = 1; k <= Wavefilt.ncof; k++) {
					jf = (n1 & (ni + k)) + 1;
					jr = (n1 & (nj + k)) + 1;
					wksp[jf - 1] += Wavefilt.cc[k] * ai;// wksp[jf] +=
														// wfilt.cc[k]*ai;
					wksp[jr - 1] += Wavefilt.cr[k] * ai1;// wksp[jr] +=
															// wfilt.cr[k]*ai1;
				}
			}
		}
		for (j = 1; j <= n; j++)
			a[j - 1] = wksp[j - 1];// a[j]=wksp[j]; Copy the results back from
									// workspace.
		// free_vector(wksp,1,n);
	}

	/*
	 * Wavelet Transform in Multidimensions
	 */
	/**
	 * Replaces a by its ndim-dimensional discrete wavelet transform, if isign is input as 1. Here 
	 * nn[1..ndim] is an integer array containing the lengths of each dimension (number of real 
	 * values), which MUST all be powers of 2. a is a real array of length equal to the product of 
	 * these lengths, in which the data are stored as in a multidimensional real array. If isign is input 
	 * as -1, a is replaced by its inverse wavelet transform. The routine wtstep, whose actual name 
	 * must be supplied in calling this routine, is the underlying wavelet filter. Examples of wtstep 
	 * are daub4 and (preceded by pwtset) pwt.
	 * @param a a
	 * @param nn nn
	 * @param ndim ndim
	 * @param isign isign
	 */
	public static void wtn(double a[], int nn[], int ndim, int isign)// ,
	// void (*wtstep)(float [], unsigned long, int))
	// Replaces a by its ndim-dimensional discrete wavelet transform, if isign
	// is input as 1. Here
	// nn[1..ndim] is an integer array containing the lengths of each dimension
	// (number of real
	// values), which MUST all be powers of 2. a is a real array of length equal
	// to the product of
	// these lengths, in which the data are stored as in a multidimensional real
	// array. If isign is input
	// as -1, a is replaced by its inverse wavelet transform. The routine
	// wtstep, whose actual name
	// must be supplied in calling this routine, is the underlying wavelet
	// filter. Examples of wtstep
	// are daub4 and (preceded by pwtset) pwt.
	{
		// unsigned long i1,i2,i3,k,n,nnew,nprev=1,nt,ntot=1;
		int i1 = 0;
		int i2 = 0;
		int i3 = 0;
		int k = 0;
		int n = 0;
		int nnew = 0;
		int nprev = 1;
		int nt = 0;
		int ntot = 1;
		int idim = 0;

		for (idim = 1; idim <= ndim; idim++)
			ntot *= nn[idim - 1];// ntot *= nn[idim];
		double[] wksp = new double[ntot];// wksp=vector(1,ntot);
		for (idim = 1; idim <= ndim; idim++) {// Main loop over the dimensions.
			n = nn[idim - 1];// n=nn[idim];
			nnew = n * nprev;
			if (n > 4) {
				for (i2 = 0; i2 < ntot; i2 += nnew) {
					for (i1 = 1; i1 <= nprev; i1++) {
						for (i3 = i1 + i2, k = 1; k <= n; k++, i3 += nprev)
							wksp[k - 1] = a[i3 - 1];// wksp[k]=a[i3];
						// Copy the relevant row or column or etc. into
						// workspace.
						if (isign >= 0) {// Do one-dimensional wavelet
											// transform.
							for (nt = n; nt >= 4; nt >>= 1)
								wtstep(wksp, nt, isign, "daub4");// (*wtstep)(wksp,nt,isign);
						} else {// Or inverse transform.
							for (nt = 4; nt <= n; nt <<= 1)
								wtstep(wksp, nt, isign, "daub4");// (*wtstep)(wksp,nt,isign);
						}
						for (i3 = i1 + i2, k = 1; k <= n; k++, i3 += nprev)
							a[i3 - 1] = wksp[k - 1];// a[i3]=wksp[k];
						// Copy back from workspace.
					}
				}
			}
			nprev = nnew;
		}
		// free_vector(wksp,1,ntot);
	}
}
