package danfulea.math.numerical;

/**
 * FastFourierTransform class
 * Based on literature, including Numerical Recipes in C (Cambridge Univ.).
 * @author Dan Fulea, 13 OCT. 2006
 */
public class FastFourierTransform {
	// public static Function func;
	public static boolean failB = false;
	public static String failS = "";
	public static double PI = 3.141592653589793;

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
	 * A very large class of important computational problems falls under the
	 * general rubric of “Fourier transform methods” or “spectral methods.” For
	 * some of these problems, the Fourier transform is simply an efficient
	 * computational tool for accomplishing certain common manipulations of
	 * data. In other cases, we have problems for which the Fourier transform
	 * (or the related “power spectrum”) is itself of intrinsic interest. These
	 * two kinds of problems share a common methodology. Largely for historical
	 * reasons the literature on Fourier and spectral methods has been disjoint
	 * from the literature on “classical” numerical analysis. Nowadays there is
	 * no justification for such a split. Fourier methods are commonplace in
	 * research and we shall not treat them as specialized or arcane. At the
	 * same time, we realize that many computer users have had relatively less
	 * experience with this field than with, say, differential equations or
	 * numerical integration. Therefore our summary of analytical results will
	 * be more complete. Numerical algorithms, per se, begin in §12.2. Various
	 * applications of Fourier transform methods are discussed in Chapter 13.
	 * 
	 *  A
	 * physical process can be described either in the time domain, by the
	 * values of some quantity h as a function of time t, e.g., h(t), or else in
	 * the frequency domain, where the process is specified by giving its
	 * amplitude H (generally a complex number indicating phase also) as a
	 * function of frequency f, that is H(f), with -inf < f < inf. For many
	 * purposes it is useful to think of h(t) and H(f) as being two different
	 * representations of the same function. One goes back and forth between
	 * these two representations by means of the Fourier transform equations,
	 * H(f) = integral from -inf to +inf of h(t)exp(2pift)dt h(t) = integral
	 * from -inf to +inf of H(f)exp(-2pift)df If you are trained as a physicist
	 * or mathematician, you are probably more used to using angular frequency
	 * omega, which is given in radians per sec. The relation between omega and
	 * f, H(omega) and H(f) is => H(omega) = integral from -inf to +inf of
	 * h(t)exp(i*omega*t)dt h(t) == integral from -inf to +inf of
	 * H(omega)exp(-i*omega*t)domega
	 * 
	 * From equation (12.0.1) it is evident at once that Fourier transformation
	 * is a linear operation. The transform of the sum of two functions is equal
	 * to the sum of the transforms. The transform of a constant times a
	 * function is that same constant times the transform of the function. In
	 * the time domain, function h(t) may happen to have one or more special
	 * symmetries It might be purely real or purely imaginary or it might be
	 * even, h(t) = h(-t), or odd, h(t) = -h(-t). In the frequency domain, these
	 * symmetries lead to relationships between H(f) and H(-f). The following
	 * table gives the correspondence between symmetries in the two domains: If
	 * . . . then . . . h(t) is real H(-f) = [H(f)]* h(t) is imaginary H(-f) =
	 * -[H(f)]* h(t) is even H(-f) = H(f) [i.e., H(f) is even] h(t) is odd H(-f)
	 * = -H(f) [i.e., H(f) is odd] h(t) is real and even H(f) is real and even
	 * h(t) is real and odd H(f) is imaginary and odd h(t) is imaginary and even
	 * H(f) is imaginary and even h(t) is imaginary and odd H(f) is real and odd
	 * 
	 * 
	 * Fourier Transform of Discretely Sampled Data
	 * 
	 * In the most common situations, function h(t) is sampled (i.e., its value
	 * is recorded) at evenly spaced intervals in time. Let delta denote the
	 * time interval between consecutive samples, so that the sequence of
	 * sampled values is hn = h(ndelta) n = . . . ,-3,-2,-1, 0, 1, 2, 3, . . .
	 * 
	 * We now estimate the Fourier transform of a function from a finite number
	 * of its sampled points. Suppose that we have N consecutive sampled values
	 * The remaining step is to approximate the integral in (12.0.1) by a
	 * discrete sum: H(fn) = integral from -inf to +inf of h(t)exp(2pifnt)dt =
	 * sum from k=0 to N-1 of hk exp(2pifntk)delta = delta * sum from k=0 to N-1
	 * of hk exp(2pikn/N) (12.1.6) Here equations (12.1.4) and (12.1.5) have
	 * been used in the final equality. The final summation in equation (12.1.6)
	 * is called the discrete Fourier transform of the N points hk. Let us
	 * denote it by Hn, Hn =sum from k=0 to N-1 of hk exp(2pikn/N) (12.1.7) The
	 * discrete Fourier transform maps N complex numbers (the hk’s) into N
	 * complex numbers (the Hn’s). It does not depend on any dimensional
	 * parameter, such as the time scale ?. The relation (12.1.6) between the
	 * discrete Fourier transform of a set of numbers and their continuous
	 * Fourier transform when they are viewed as samples of a continuous
	 * function sampled at an interval ? can be rewritten as H(fn) ? ?Hn
	 * (12.1.8) where fn is given by (12.1.5). Up to now we have taken the view
	 * that the index n in
	 * 
	 * 
	 * Fast Fourier Transform (FFT) How much computation is involved in
	 * computing the discrete Fourier transform (12.1.7) of N points? The FFT
	 * routine given below is based on one originally written by N. M. Brenner.
	 * The input quantities are the number of complex data points (nn), the data
	 * array (data[1..2*nn]), and isign, which should be set to either ±1 and is
	 * the sign of i in the exponential of equation (12.1.7). When isign is set
	 * to -1, the routine thus calculates the inverse transform (12.1.9) —
	 * except that it does not multiply by the normalizing factor 1/N that
	 * appears in that equation. You can do that yourself. Notice that the
	 * argument nn is the number of complex data points. The actual length of
	 * the real array (data[1..2*nn]) is 2 times nn, with each complex value
	 * occupying two consecutive locations. In other words, data[1] is the real
	 * part of f0, data[2] is the imaginary part of f0, and so on up to
	 * data[2*nn-1], which is the real part of fN-1, and data[2*nn], which is
	 * the imaginary part of fN-1. The FFT routine gives back the Fn’s packed in
	 * exactly the same fashion, as nn complex numbers. The real and imaginary
	 * parts of the zero frequency component F 0 are in data[1] and data[2]; the
	 * smallest nonzero positive frequency has real and imaginary parts in
	 * data[3] and data[4]; the smallest (in magnitude) nonzero negative
	 * frequency has real and imaginary parts in data[2*nn-1] and data[2*nn].
	 * Positive frequencies increasing in magnitude are stored in the
	 * real-imaginary pairs data[5], data[6] up to data[nn-1], data[nn].
	 * Negative frequencies of increasing magnitude are stored in data[2*nn-3],
	 * data[2*nn-2] down to data[nn+3], data[nn+4]. Finally, the pair
	 * data[nn+1], data[nn+2] contain the real and imaginary parts of the one
	 * aliased point that contains the most positive and the most negative
	 * frequency. You should try to develop a familiarity with this storage
	 * arrangement of complex spectra, also shown in Figure 12.2.2, since it is
	 * the practical standard. This is a good place to remind you that you can
	 * also use a routine like four1 without modification even if your input
	 * data array is zero-offset, that is has the range data[0..2*nn-1]. In this
	 * case, simply decrement the pointer to data by one when four1 is invoked,
	 * e.g., four1(data-1,1024,1);. The real part of f0 will now be returned in
	 * data[0], the imaginary part in data[1], and so on. See §1.2.
	 */

	// #define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr
	// public static void four1(double[] data, unsigned long nn, int isign)
	
	/**
	 * Replaces data[1..2*nn] by its discrete Fourier transform, if isign is input as 1; or replaces
	 * data[1..2*nn] by nn times its inverse discrete Fourier transform, if isign is input as -1.
	 * data is a complex array of length nn or, equivalently, a real array of length 2*nn. nn MUST
	 * be an integer power of 2 (this is not checked for!).
	 * @param data data
	 * @param nn nn
	 * @param isign isign
	 */
	public static void four1(double[] data, int nn, int isign)
	// Replaces data[1..2*nn] by its discrete Fourier transform, if isign is
	// input as 1; or replaces
	// data[1..2*nn] by nn times its inverse discrete Fourier transform, if
	// isign is input as -1.
	// data is a complex array of length nn or, equivalently, a real array of
	// length 2*nn. nn MUST
	// be an integer power of 2 (this is not checked for!).
	{
		// unsigned long n,mmax,m,j,istep,i;
		int n = 0;
		int mmax = 0;
		int m = 0;
		int j = 0;
		int istep = 0;
		int i = 0;
		double wtemp = 0.0;
		double wr = 0.0;
		double wpr = 0.0;
		double wpi = 0.0;
		double wi = 0.0;
		double theta = 0.0;// Double precision for the trigonometric
							// recurrences.
		double tempr = 0.0;
		double tempi = 0.0;

		n = nn << 1;
		j = 1;
		for (i = 1; i < n; i += 2) {// This is the bit-reversal section of the
									// routine.
			if (j > i) {
				// SWAP(data[j],data[i]); Exchange the two complex numbers.
				tempr = data[j - 1];
				data[j - 1] = data[i - 1];
				data[i - 1] = tempr;
				// SWAP(data[j+1],data[i+1]);
				tempr = data[j];
				data[j] = data[i];
				data[i] = tempr;
			}
			m = nn;
			while (m >= 2 && j > m) {
				j -= m;
				m >>= 1;
			}
			j += m;
		}
		// Here begins the Danielson-Lanczos section of the routine.
		mmax = 2;
		while (n > mmax) {// Outer loop executed log2 nn times.
			istep = mmax << 1;
			theta = isign * (6.28318530717959 / mmax);// Initialize the
														// trigonometric
														// recurrence.
			wtemp = Math.sin(0.5 * theta);
			wpr = -2.0 * wtemp * wtemp;
			wpi = Math.sin(theta);
			wr = 1.0;
			wi = 0.0;
			for (m = 1; m < mmax; m += 2) {// Here are the two nested inner
											// loops.
				for (i = m; i <= n; i += istep) {
					j = i + mmax;// This is the Danielson-Lanczos formula:
					// tempr=wr*data[j]-wi*data[j+1];
					tempr = wr * data[j - 1] - wi * data[j];
					// tempi=wr*data[j+1]+wi*data[j];
					tempi = wr * data[j] + wi * data[j - 1];
					// data[j]=data[i]-tempr;
					data[j - 1] = data[i - 1] - tempr;
					// data[j+1]=data[i+1]-tempi;
					data[j] = data[i] - tempi;
					data[i - 1] += tempr;// data[i] += tempr;
					data[i] += tempi;// data[i+1] += tempi;
				}
				wr = (wtemp = wr) * wpr - wi * wpi + wr; // Trigonometric
															// recurrence.
				wi = wi * wpr + wtemp * wpi + wi;
			}
			mmax = istep;
		}
	}

	/*
	 * 12.3 FFT of Real Functions, Sine and Cosine Transforms
	 * 
	 * Transform of Two Real Functions Simultaneously
	 */

	/**
	 * Given two real input arrays data1[1..n] and data2[1..n], this routine calls four1 and 
	 * returns two complex output arrays, fft1[1..2n] and fft2[1..2n], each of complex length 
	 * n (i.e., real length 2*n), which contain the discrete Fourier transforms 
	 * of the respective data arrays. n MUST be an integer power of 2.
	 * @param data1 data1
	 * @param data2 data2
	 * @param fft1 fft1
	 * @param fft2 fft2
	 * @param n n
	 */
	public static void twofft(double[] data1, double[] data2, double[] fft1,
			double[] fft2, int n)// unsigned long n)
	// Given two real input arrays data1[1..n] and data2[1..n], this routine
	// calls four1 and
	// returns two complex output arrays, fft1[1..2n] and fft2[1..2n], each of
	// complex length
	// n (i.e., real length 2*n), which contain the discrete Fourier transforms
	// of the respective data
	// arrays. n MUST be an integer power of 2.
	{
		// void four1(float data[], unsigned long nn, int isign);
		// unsigned long nn3,nn2,jj,j;
		int nn3 = 0;
		int nn2 = 0;
		int jj = 0;
		int j = 0;
		double rep = 0.0;
		double rem = 0.0;
		double aip = 0.0;
		double aim = 0.0;

		nn3 = 1 + (nn2 = 2 + n + n);
		for (j = 1, jj = 2; j <= n; j++, jj += 2) { // Pack the two real arrays
													// into one complex array.
			fft1[jj - 2] = data1[j - 1];// fft1[jj-1]=data1[j];
			fft1[jj - 1] = data2[j - 1];// fft1[jj]=data2[j];
		}
		four1(fft1, n, 1); // Transform the complex array.
		fft2[0] = fft1[1];// fft2[1]=fft1[2];
		fft1[1] = fft2[1] = 0.0;// fft1[2]=fft2[2]=0.0;
		for (j = 3; j <= n + 1; j += 2) {
			// rep=0.5*(fft1[j]+fft1[nn2-j]); Use symmetries to separate the two
			// transforms.
			rep = 0.5 * (fft1[j - 1] + fft1[nn2 - j - 1]);
			rem = 0.5 * (fft1[j - 1] - fft1[nn2 - j - 1]);// rem=0.5*(fft1[j]-fft1[nn2-j]);
			aip = 0.5 * (fft1[j] + fft1[nn3 - j - 1]);// aip=0.5*(fft1[j+1]+fft1[nn3-j]);
			aim = 0.5 * (fft1[j] - fft1[nn3 - j - 1]);// aim=0.5*(fft1[j+1]-fft1[nn3-j]);
			fft1[j - 1] = rep;// fft1[j]=rep; Ship them out in two complex
								// arrays.
			fft1[j] = aim;// fft1[j+1]=aim;
			fft1[nn2 - j - 1] = rep;// fft1[nn2-j]=rep;
			fft1[nn3 - j - 1] = -aim;// fft1[nn3-j] = -aim;
			fft2[j - 1] = aip;// fft2[j]=aip;
			fft2[j] = -rem;// fft2[j+1] = -rem;
			fft2[nn2 - j - 1] = aip;// fft2[nn2-j]=aip;
			fft2[nn3 - j - 1] = rem;// fft2[nn3-j]=rem;
		}
	}

	/*
	 * To implement the second method, which allows us to perform the FFT of a
	 * single real function without redundancy, we split the data set in half,
	 * thereby forming two real arrays of half the size. We can apply the
	 * program above to these two, but of course the result will not be the
	 * transform of the original data. It will be a schizophrenic combination of
	 * two transforms, each of which has half of the information we need.
	 * Fortunately, this schizophrenia is treatable. It works like this: The
	 * right way to split the original data is to take the even-numbered f j as
	 * one data set, and the odd-numbered fj as the other. The beauty of this is
	 * that we can take the original real array and treat it as a complex array
	 * hj of half the length. The first data set is the real part of this array,
	 * and the second is the imaginary part, as prescribed for twofft. No
	 * repacking is required. In other words hj = f2j +if2j+1, j= 0, . .
	 * .,N/2-1. We submit this to four1, and it will give back a complex array
	 * Hn = Fe n + iF o n, n= 0, . . .,N/2 - 1 with................
	 */
	/**
	 * Calculates the Fourier transform of a set of n real-valued data points. 
	 * Replaces this data (which is stored in array data[1..n]) by the positive frequency half of its 
	 * complex Fourier transform. The real-valued first and last components of the complex transform are 
	 * returned as elements data[1] and data[2], respectively. n must be a power of 2. This routine 
	 * also calculates the inverse transform of a complex data array if it is the transform of real 
	 * data. (Result in this case must be multiplied by 2/n.)
	 * @param data data
	 * @param n n
	 * @param isign isign of 1 means the forward transform!
	 */
	public static void realft(double[] data, int n, int isign)// (float data[],
																// unsigned long
																// n, int isign)
	// Calculates the Fourier transform of a set of n real-valued data points.
	// Replaces this data (which
	// is stored in array data[1..n]) by the positive frequency half of its
	// complex Fourier transform.
	// The real-valued first and last components of the complex transform are
	// returned as elements
	// data[1] and data[2], respectively. n must be a power of 2. This routine
	// also calculates the
	// inverse transform of a complex data array if it is the transform of real
	// data. (Result in this case
	// must be multiplied by 2/n.)
	{
		// void four1(float data[], unsigned long nn, int isign);
		// unsigned long i,i1,i2,i3,i4,np3;
		int i = 0;
		int i1 = 0;
		int i2 = 0;
		int i3 = 0;
		int i4 = 0;
		int np3 = 0;
		double c1 = 0.5;
		double c2 = 0.0;
		double h1r = 0.0;
		double h1i = 0.0;
		double h2r = 0.0;
		double h2i = 0.0;
		double wr = 0.0;
		double wi = 0.0;
		double wpr = 0.0;
		double wpi = 0.0;
		double wtemp = 0.0;
		double theta = 0.0;
		// Double precision for the trigonometric recurrences.
		theta = 3.141592653589793 / (double) (n >> 1); // Initialize the
														// recurrence.
		if (isign == 1) {
			c2 = -0.5;
			four1(data, n >> 1, 1); // The forward transform is here.
		} else {
			c2 = 0.5; // Otherwise set up for an inverse transform.
			theta = -theta;
		}
		wtemp = Math.sin(0.5 * theta);
		wpr = -2.0 * wtemp * wtemp;
		wpi = Math.sin(theta);
		wr = 1.0 + wpr;
		wi = wpi;
		np3 = n + 3;
		for (i = 2; i <= (n >> 2); i++) {// Case i=1 done separately below.
			i4 = 1 + (i3 = np3 - (i2 = 1 + (i1 = i + i - 1)));
			// h1r=c1*(data[i1]+data[i3]); //The two separate transforms are
			// separated out of data.
			h1r = c1 * (data[i1 - 1] + data[i3 - 1]);
			h1i = c1 * (data[i2 - 1] - data[i4 - 1]);// h1i=c1*(data[i2]-data[i4]);
			h2r = -c2 * (data[i2 - 1] + data[i4 - 1]);// h2r =
														// -c2*(data[i2]+data[i4]);
			h2i = c2 * (data[i1 - 1] - data[i3 - 1]);// h2i=c2*(data[i1]-data[i3]);
			data[i1 - 1] = h1r + wr * h2r - wi * h2i;// data[i1]=h1r+wr*h2r-wi*h2i;
			// Here they are recombined to form the true transform of the
			// original real data.
			data[i2 - 1] = h1i + wr * h2i + wi * h2r;// data[i2]=h1i+wr*h2i+wi*h2r;
			data[i3 - 1] = h1r - wr * h2r + wi * h2i;// data[i3]=h1r-wr*h2r+wi*h2i;
			data[i4 - 1] = -h1i + wr * h2i + wi * h2r;// data[i4] =
														// -h1i+wr*h2i+wi*h2r;
			wr = (wtemp = wr) * wpr - wi * wpi + wr; // The recurrence.
			wi = wi * wpr + wtemp * wpi + wi;
		}
		if (isign == 1) {
			data[0] = (h1r = data[0]) + data[1];// data[1] =
												// (h1r=data[1])+data[2];
			// Squeeze the first and last data together to get them all within
			// the original array.
			data[1] = h1r - data[1];// data[2] = h1r-data[2];
		} else {
			data[0] = c1 * ((h1r = data[0]) + data[1]);// data[1]=c1*((h1r=data[1])+data[2]);
			data[1] = c1 * (h1r - data[1]);// data[2]=c1*(h1r-data[2]);
			four1(data, n >> 1, -1); // This is the inverse transform for the
										// case isign=-1.
		}
	}

	/*
	 * Fast Sine and Cosine Transforms Among their other uses, the Fourier
	 * transforms of functions can be used to solve differential equations (see
	 * §19.4). The most common boundary conditions for the solutions are 1) they
	 * have the value zero at the boundaries, or 2) their derivatives are zero
	 * at the boundaries. In the first instance, the natural transform to use is
	 * the sine transform, given by....
	 */
	/**
	 * Calculates the sine transform of a set of n real-valued data points stored in array y[1..n]. 
	 * The number n must be a power of 2. On exit y is replaced by its transform. This program, 
	 * without changes, also calculates the inverse sine transform, but in this case the output array 
	 * should be multiplied by 2/n.
	 * @param y y
	 * @param n n
	 */
	public static void sinft(double[] y, int n)
	// Calculates the sine transform of a set of n real-valued data points
	// stored in array y[1..n].
	// The number n must be a power of 2. On exit y is replaced by its
	// transform. This program,
	// without changes, also calculates the inverse sine transform, but in this
	// case the output array
	// should be multiplied by 2/n.
	{
		// void realft(float data[], unsigned long n, int isign);
		int j = 0;
		int n2 = n + 2;
		double sum = 0.0;
		double y1 = 0.0;
		double y2 = 0.0;
		double theta = 0.0;
		double wi = 0.0;
		double wr = 1.0;
		double wpi = 0.0;
		double wpr = 0.0;
		double wtemp = 0.0;
		// Double precision in the trigonometric recurrences.
		theta = 3.14159265358979 / (double) n; // Initialize the recurrence.
		wtemp = Math.sin(0.5 * theta);
		wpr = -2.0 * wtemp * wtemp;
		wpi = Math.sin(theta);
		y[0] = 0.0;// y[1]=0.0;
		for (j = 2; j <= (n >> 1) + 1; j++) {
			wr = (wtemp = wr) * wpr - wi * wpi + wr;// Calculate the sine for
													// the auxiliary array.
			wi = wi * wpr + wtemp * wpi + wi; // The cosine is needed to
												// continue the recurrence.
			y1 = wi * (y[j - 1] + y[n2 - j - 1]);// y1=wi*(y[j]+y[n2-j]);
													// Construct the auxiliary
													// array.
			y2 = 0.5 * (y[j - 1] - y[n2 - j - 1]);// y2=0.5*(y[j]-y[n2-j]);
			y[j - 1] = y1 + y2;// y[j]=y1+y2; Terms j and N - j are related.
			y[n2 - j - 1] = y1 - y2;// y[n2-j]=y1-y2;
		}
		realft(y, n, 1); // Transform the auxiliary array.
		y[0] *= 0.5;// y[1]*=0.5; Initialize the sum used for odd terms below.
		sum = y[1] = 0.0;// sum=y[2]=0.0;
		for (j = 1; j <= n - 1; j += 2) {
			sum += y[j - 1];// sum += y[j];
			y[j - 1] = y[j];// y[j]=y[j+1]; Even terms determined directly.
			y[j] = sum;// y[j+1]=sum; Odd terms determined by this running sum.
		}
	}

	/*
	 * The other common boundary condition for differential equations is that
	 * the derivative of the function is zero at the boundary. In this case the
	 * natural transform is the cosine transform. There are several possible
	 * ways of defining the transform. Each can be thought of as resulting from
	 * a different way of extending a given array to create an even array of
	 * double the length, and/or from whether the extended array contains 2N -
	 * 1, 2N, or some other number of points. In practice, only two of the
	 * numerous possibilities are useful so we will restrict ourselves to just
	 * these two.
	 */

	/**
	 * Calculates the cosine transform of a set y[1..n+1] of real-valued data points. The transformed 
	 * data replace the original data in array y. n must be a power of 2. This program, without 
	 * changes, also calculates the inverse cosine transform, but in this case the output array should 
	 * be multiplied by 2/n.
	 * @param y y
	 * @param n n
	 */
	public static void cosft1(double[] y, int n)
	// Calculates the cosine transform of a set y[1..n+1] of real-valued data
	// points. The transformed
	// data replace the original data in array y. n must be a power of 2. This
	// program, without
	// changes, also calculates the inverse cosine transform, but in this case
	// the output array should
	// be multiplied by 2/n.
	{
		// void realft(float data[], unsigned long n, int isign);
		int j = 0;
		int n2 = 0;
		double sum = 0.0;
		double y1 = 0.0;
		double y2 = 0.0;
		double theta = 0.0;
		double wi = 0.0;
		double wpi = 0.0;
		double wpr = 0.0;
		double wr = 1.0;
		double wtemp = 0.0;
		// Double precision for the trigonometric recurrences.
		theta = PI / n; // Initialize the recurrence.
		wtemp = Math.sin(0.5 * theta);
		wpr = -2.0 * wtemp * wtemp;
		wpi = Math.sin(theta);
		sum = 0.5 * (y[0] - y[n]);// sum=0.5*(y[1]-y[n+1]);
		y[0] = 0.5 * (y[0] + y[n]);// y[1]=0.5*(y[1]+y[n+1]);
		n2 = n + 2;
		for (j = 2; j <= (n >> 1); j++) { // j=n/2+1 unnecessary since y[n/2+1]
											// unchanged.
			wr = (wtemp = wr) * wpr - wi * wpi + wr; // Carry out the
														// recurrence.
			wi = wi * wpr + wtemp * wpi + wi;
			y1 = 0.5 * (y[j - 1] + y[n2 - j - 1]);// y1=0.5*(y[j]+y[n2-j]);
													// Calculate the auxiliary
													// function.
			y2 = (y[j - 1] - y[n2 - j - 1]);// y2=(y[j]-y[n2-j]);
			y[j - 1] = y1 - wi * y2;// y[j]=y1-wi*y2; The values for j and N - j
									// are related.
			y[n2 - j - 1] = y1 + wi * y2;// y[n2-j]=y1+wi*y2;
			sum += wr * y2;// Carry along this sum for later use in unfolding
							// the transform.
		}
		realft(y, n, 1); // Calculate the transform of the auxiliary function.
		y[n] = y[1];// y[n+1]=y[2];
		y[1] = sum;// y[2]=sum; sum is the value of F1 in equation (12.3.21).
		for (j = 4; j <= n; j += 2) {
			sum += y[j - 1];// sum += y[j]; Equation (12.3.20).
			y[j - 1] = sum;// y[j]=sum;
		}
	}

	/*
	 * The second important form of the cosine transform is defined by Fk = N-1
	 * j=0 fj cos pk(j + 1 2 ) N (12.3.22) with inverse fj = 2 N N-1  k=0 Fk
	 * cos pk(j + 1 2 ) N
	 * 
	 * This form of the cosine transform is useful when solving differential
	 * equations on “staggered” grids, where the variables are centered midway
	 * between mesh points. It is also the standard form in the field of data
	 * compression and image processing.
	 */
	/**
	 * Calculates the “staggered” cosine transform of a set y[1..n] of real-valued data points. The 
	 * transformed data replace the original data in array y. n must be a power of 2. Set isign to 
	 * +1 for a transform, and to -1 for an inverse transform. For an inverse transform, the output 
	 * array should be multiplied by 2/n.
	 * @param y y
	 * @param n n
	 * @param isign isign
	 */
	public static void cosft2(double[] y, int n, int isign)
	// Calculates the “staggered” cosine transform of a set y[1..n] of
	// real-valued data points. The
	// transformed data replace the original data in array y. n must be a power
	// of 2. Set isign to
	// +1 for a transform, and to -1 for an inverse transform. For an inverse
	// transform, the output
	// array should be multiplied by 2/n.
	{
		// void realft(float data[], unsigned long n, int isign);
		int i = 0;
		double sum = 0.0;
		double sum1 = 0.0;
		double y1 = 0.0;
		double y2 = 0.0;
		double ytemp = 0.0;
		double theta = 0.0;
		double wi = 0.0;
		double wi1 = 0.0;
		double wpi = 0.0;
		double wpr = 0.0;
		double wr = 1.0;
		double wr1 = 0.0;
		double wtemp = 0.0;
		// Double precision for the trigonometric recurrences.
		theta = 0.5 * PI / n; // Initialize the recurrences.
		wr1 = Math.cos(theta);
		wi1 = Math.sin(theta);
		wpr = -2.0 * wi1 * wi1;
		wpi = Math.sin(2.0 * theta);
		if (isign == 1) {// Forward transform.
			for (i = 1; i <= n / 2; i++) {
				y1 = 0.5 * (y[i - 1] + y[n - i]);// y1=0.5*(y[i]+y[n-i+1]);
													// Calculate the auxiliary
													// function.
				y2 = wi1 * (y[i - 1] - y[n - i]);// y2=wi1*(y[i]-y[n-i+1]);
				y[i - 1] = y1 + y2;// y[i]=y1+y2;
				y[n - i] = y1 - y2;// y[n-i+1]=y1-y2;
				wr1 = (wtemp = wr1) * wpr - wi1 * wpi + wr1; // Carry out the
																// recurrence.
				wi1 = wi1 * wpr + wtemp * wpi + wi1;
			}
			realft(y, n, 1);// Transform the auxiliary function.
			for (i = 3; i <= n; i += 2) {// Even terms.
				wr = (wtemp = wr) * wpr - wi * wpi + wr;
				wi = wi * wpr + wtemp * wpi + wi;
				y1 = y[i - 1] * wr - y[i] * wi;// y1=y[i]*wr-y[i+1]*wi;
				y2 = y[i] * wr + y[i - 1] * wi;// y2=y[i+1]*wr+y[i]*wi;
				y[i - 1] = y1;// y[i]=y1;
				y[i] = y2;// y[i+1]=y2;
			}
			sum = 0.5 * y[1];// sum=0.5*y[2]; Initialize recurrence for odd
								// terms with 1 2RN/2.
			for (i = n; i >= 2; i -= 2) {
				sum1 = sum;// Carry out recurrence for odd terms.
				sum += y[i - 1];// sum += y[i];
				y[i - 1] = sum1;// y[i]=sum1;
			}
		} else if (isign == -1) {// Inverse transform.
			ytemp = y[n - 1];// ytemp=y[n];
			for (i = n; i >= 4; i -= 2)
				// y[i]=y[i-2]-y[i]; Form difference of odd terms.
				y[i - 1] = y[i - 3] - y[i - 1];
			y[1] = 2.0 * ytemp;// y[2]=2.0*ytemp;
			for (i = 3; i <= n; i += 2) {// Calculate Rk and Ik.
				wr = (wtemp = wr) * wpr - wi * wpi + wr;
				wi = wi * wpr + wtemp * wpi + wi;
				y1 = y[i - 1] * wr + y[i] * wi;// y1=y[i]*wr+y[i+1]*wi;
				y2 = y[i] * wr - y[i - 1] * wi;// y2=y[i+1]*wr-y[i]*wi;
				y[i - 1] = y1;// y[i]=y1;
				y[i] = y2;// y[i+1]=y2;
			}
			realft(y, n, -1);
			for (i = 1; i <= n / 2; i++) {// Invert auxiliary array.
				y1 = y[i - 1] + y[n - i];// y1=y[i]+y[n-i+1];
				y2 = (0.5 / wi1) * (y[i - 1] - y[n - i]);// y2=(0.5/wi1)*(y[i]-y[n-i+1]);
				y[i - 1] = 0.5 * (y1 + y2);// y[i]=0.5*(y1+y2);
				y[n - i] = 0.5 * (y1 - y2);// y[n-i+1]=0.5*(y1-y2);
				wr1 = (wtemp = wr1) * wpr - wi1 * wpi + wr1;
				wi1 = wi1 * wpr + wtemp * wpi + wi1;
			}
		}
	}

	/*
	 * FFT in Two or More Dimensions Given a complex function h(k1, k2) defined
	 * over the two-dimensional grid 0 = k1 = N1 - 1, 0 = k2 = N2 - 1, we can
	 * define its two-dimensional discrete Fourier transform as a complex
	 * function H(n1, n2), defined over the same grid, H(n1, n2) = N2-1  k2=0
	 * N1-1  k1=0 exp(2pik2n2/N2) exp(2pik1n1/N1) h(k1, k2)
	 * 
	 * Arrangement of a two-dimensional FFT: data[1] Re data[1] Im,.....pana la
	 * 2N2 float numbers=ROW1
	 * ....................................................=ROW2
	 * ......................................................
	 * ......................................................ROW N1/2
	 * ......................................................ROW N1/2+1
	 * ......................................................ROW N1/2+2
	 * ......................................................
	 * ................................................data[2N1N2]=ROW N1
	 * 
	 * The total number of (real)array elements is 2N1N2. A few words about the
	 * data array: fourn accesses it as a one-dimensional array of real numbers,
	 * that is, data[1..(2N1N2 . . . NL)], of length equal to twice the product
	 * of the lengths of the L dimensions. It assumes that the array represents
	 * an L-dimensional complex array, with individual components ordered as
	 * follows: (i) each complex value occupies two sequential locations, real
	 * part followed by imaginary; (ii) the first subscript changes least
	 * rapidly as one goes through the array; the last subscript changes most
	 * rapidly (that is, “store by rows,” the C norm); (iii) subscripts range
	 * from 1 to their maximum values (N1,N2, . . . , NL, respectively), rather
	 * than from 0 to N1 -1, N2 -1, . . . , NL -1. Almost all failures to get
	 * fourn to work result from improper understanding of the above ordering of
	 * the data array, so take care! (Figure 12.4.1 illustrates the format of
	 * the output array.)
	 */
	// public static void fourn(double[] data, unsigned long nn[], int ndim, int
	// isign)
	/**
	 * Replaces data by its ndim-dimensional discrete Fourier transform, if isign is input as 1. 
	 * nn[1..ndim] is an integer array containing the lengths of each dimension (number of complex values), 
	 * which MUST all be powers of 2. data is a real array of length twice the product of 
	 * these lengths, in which the data are stored as in a multidimensional complex array: real and 
	 * imaginary parts of each element are in consecutive locations, and the rightmost index of the 
	 * array increases most rapidly as one proceeds along data. For a two-dimensional array, this is 
	 * equivalent to storing the array by rows. If isign is input as -1, data is replaced by its inverse 
	 * transform times the product of the lengths of all dimensions.
	 * @param data data
	 * @param nn nn
	 * @param ndim ndim
	 * @param isign isign
	 */
	public static void fourn(double[] data, int[] nn, int ndim, int isign)
	// Replaces data by its ndim-dimensional discrete Fourier transform, if
	// isign is input as 1.
	// nn[1..ndim] is an integer array containing the lengths of each dimension
	// (number of complex
	// values), which MUST all be powers of 2. data is a real array of length
	// twice the product of
	// these lengths, in which the data are stored as in a multidimensional
	// complex array: real and
	// /imaginary parts of each element are in consecutive locations, and the
	// rightmost index of the
	// array increases most rapidly as one proceeds along data. For a
	// two-dimensional array, this is
	// equivalent to storing the array by rows. If isign is input as -1, data is
	// replaced by its inverse
	// transform times the product of the lengths of all dimensions.
	{
		int idim = 0;
		// unsigned long i1,i2,i3,i2rev,i3rev,ip1,ip2,ip3,ifp1,ifp2;
		int i1 = 0;
		int i2 = 0;
		int i3 = 0;
		int i2rev = 0;
		int i3rev = 0;
		int ip1 = 0;
		int ip2 = 0;
		int ip3 = 0;
		int ifp1 = 0;
		int ifp2 = 0;
		// unsigned long ibit,k1,k2,n,nprev,nrem,ntot;
		int ibit = 0;
		int k1 = 0;
		int k2 = 0;
		int n = 0;
		int nprev = 0;
		int nrem = 0;
		int ntot = 0;
		double tempi = 0.0;
		double tempr = 0.0;
		double theta = 0.0;
		double wi = 0.0;
		double wpi = 0.0;
		double wpr = 0.0;
		double wr = 0.0;
		double wtemp = 0.0; // Double precision for trigonometric recurrences.

		for (ntot = 1, idim = 1; idim <= ndim; idim++)
			// Compute total number of complex values.
			ntot *= nn[idim - 1];// ntot *= nn[idim];
		nprev = 1;
		for (idim = ndim; idim >= 1; idim--) {// Main loop over the dimensions.
			n = nn[idim - 1];// n=nn[idim];
			nrem = ntot / (n * nprev);
			/*
			 * Siftarea cu semn la stânga reprezinta o operatie identica cu
			 * înmultirea cu 2 de n ori, unde n este al doilea operand. Siftarea
			 * cu semn la dreapta reprezinta împartirea întreaga
			 */
			ip1 = nprev << 1;// =>init=1
			ip2 = ip1 * n;
			ip3 = ip2 * nrem;
			i2rev = 1;
			for (i2 = 1; i2 <= ip2; i2 += ip1) {// This is the bit-reversal
												// section of the routine.
				if (i2 < i2rev) {
					for (i1 = i2; i1 <= i2 + ip1 - 2; i1 += 2) {
						for (i3 = i1; i3 <= ip3; i3 += ip2) {
							i3rev = i2rev + i3 - i2;
							// SWAP(data[i3],data[i3rev]);tempr=(a);(a)=(b);(b)=tempr
							tempr = data[i3 - 1];
							data[i3 - 1] = data[i3rev - 1];
							data[i3rev - 1] = tempr;
							// SWAP(data[i3+1],data[i3rev+1]);
							tempr = data[i3];
							data[i3] = data[i3rev];
							data[i3rev] = tempr;
						}
					}
				}
				ibit = ip2 >> 1;
				while (ibit >= ip1 && i2rev > ibit) {
					i2rev -= ibit;
					ibit >>= 1;
				}
				i2rev += ibit;
			}
			ifp1 = ip1;// Here begins the Danielson-Lanczos section of the
						// routine.
			while (ifp1 < ip2) {
				ifp2 = ifp1 << 1;
				theta = isign * 6.28318530717959 / (ifp2 / ip1); // Initialize
																	// for the
																	// trig.
																	// recurrence.
				wtemp = Math.sin(0.5 * theta);
				wpr = -2.0 * wtemp * wtemp;
				wpi = Math.sin(theta);
				wr = 1.0;
				wi = 0.0;
				for (i3 = 1; i3 <= ifp1; i3 += ip1) {
					for (i1 = i3; i1 <= i3 + ip1 - 2; i1 += 2) {
						for (i2 = i1; i2 <= ip3; i2 += ifp2) {
							k1 = i2;// Danielson-Lanczos formula:
							k2 = k1 + ifp1;
							// tempr=(float)wr*data[k2]-(float)wi*data[k2+1];
							tempr = (double) wr * data[k2 - 1] - (double) wi
									* data[k2];
							// tempi=(double)wr*data[k2+1]+(double)wi*data[k2];
							tempi = (double) wr * data[k2] + (double) wi
									* data[k2 - 1];
							// data[k2]=data[k1]-tempr;
							data[k2 - 1] = data[k1 - 1] - tempr;
							// data[k2+1]=data[k1+1]-tempi;
							data[k2] = data[k1] - tempi;
							data[k1 - 1] += tempr;// data[k1] += tempr;
							data[k1] += tempi;// data[k1+1] += tempi;
						}
					}
					wr = (wtemp = wr) * wpr - wi * wpi + wr; // Trigonometric
																// recurrence.
					wi = wi * wpr + wtemp * wpi + wi;
				}
				ifp1 = ifp2;
			}
			nprev *= n;
		}
	}

	/*
	 * Fourier Transforms of Real Data in Two and Three Dimensions
	 * 
	 * Two-dimensional FFTs are particularly important in the field of image
	 * processing.
	 * 
	 * @@@@@@@An image is usually represented as a two-dimensional array of
	 * pixel intensities, real (and usually positive) numbers. One commonly
	 * desires to filter high, or low, frequency spatial components from an
	 * image; or to convolve or deconvolve the image with some instrumental
	 * point spread function. Use of the FFT is by far the most efficient
	 * technique,
	 * 
	 * In three dimensions, a common use of the FFT is to solve Poisson’s
	 * equation for a potential (e.g., electromagnetic or gravitational) on a
	 * three-dimensional lattice that represents the discretization of
	 * three-dimensional space. Here the source terms (mass or charge
	 * distribution) and the desired potentials are also real. In two and three
	 * dimensions, with large arrays, memory is often at a premium. It is
	 * therefore important to perform the FFTs, insofar as possible, on the data
	 * “in place.” We want a routine with functionality similar to the
	 * multidimensional FFT routine fourn (§12.4), but which operates on real,
	 * not complex, input data. We give such a routine in this section. The
	 * development is analogous to that of §12.3 leading to the one-dimensional
	 * routine realft. (You might wish to review that material at this point,
	 * particularly equation 12.3.5.)
	 */

	// public static void rlft3(double[][][] data, double[][] speq, unsigned
	// long nn1, unsigned long nn2,
	// unsigned long nn3, int isign)
	/**
	 * Given a three-dimensional real array data[1..nn1][1..nn2][1..nn3] (where nn1 = 1 for 
	 * the case of a logically two-dimensional array), this routine returns (for isign=1) the complex 
	 * fast Fourier transform as two complex arrays: On output, data contains the zero and positive 
	 * frequency values of the third frequency component, while speq[1..nn1][1..2*nn2] contains 
	 * the Nyquist critical frequency values of the third frequency component. 
	 * First (and second) frequency components are stored for zero, positive, and negative frequencies, in standard wraparound 
	 * order. For isign=-1, the inverse transform (times nn1*nn2*nn3/2 as a constant multiplicative 
	 * factor) is performed, with output data (viewed as a real array) deriving from input data 
	 * (viewed as complex) and speq. The dimensions nn1, nn2, nn3 must always be integer powers of 2.
	 
	 * @param data data
	 * @param speq speq
	 * @param nn1 nn1
	 * @param nn2 nn2
	 * @param nn3 nn3
	 * @param isign isign
	 */
	public static void rlft3(double[][][] data, double[][] speq, int nn1,
			int nn2, int nn3, int isign)
	// Given a three-dimensional real array data[1..nn1][1..nn2][1..nn3] (where
	// nn1 = 1 for
	// the case of a logically two-dimensional array), this routine returns (for
	// isign=1) the complex
	// fast Fourier transform as two complex arrays: On output, data contains
	// the zero and positive
	// frequency values of the third frequency component, while
	// speq[1..nn1][1..2*nn2] contains
	// the Nyquist critical frequency values of the third frequency component.
	// First (and second)
	// frequency components are stored for zero, positive, and negative
	// frequencies, in standard wraparound
	// order. See text for description of how complex values are arranged. For
	// isign=-1, the
	// inverse transform (times nn1*nn2*nn3/2 as a constant multiplicative
	// factor) is performed,
	// with output data (viewed as a real array) deriving from input data
	// (viewed as complex) and
	// speq. For inverse transforms on data not generated first by a forward
	// transform, make sure
	// the complex input data array satisfies property (12.5.2). The dimensions
	// nn1, nn2, nn3 must
	// always be integer powers of 2.
	{
		// void fourn(float data[], unsigned long nn[], int ndim, int isign);
		// void nrerror(char error_text[]);
		failB = false;
		// unsigned long i1,i2,i3,j1,j2,j3,nn[4],ii3;
		int i1 = 0;
		int i2 = 0;
		int i3 = 0;
		int j1 = 0;
		int j2 = 0;
		int j3 = 0;
		int[] nn = new int[4];
		int ii3 = 0;
		;
		double theta = 0.0;
		double wi = 0.0;
		double wpi = 0.0;
		double wpr = 0.0;
		double wr = 0.0;
		double wtemp = 0.0;
		double c1 = 0.0;
		double c2 = 0.0;
		double h1r = 0.0;
		double h1i = 0.0;
		double h2r = 0.0;
		double h2i = 0.0;

		// if (1+&data[nn1][nn2][nn3]-&data[1][1][1] != nn1*nn2*nn3)
		// if (1+data[nn1-1][nn2-1][nn3-1]-data[0][0][0] != nn1*nn2*nn3)
		// {
		// nrerror("rlft3: problem with dimensions or contiguity of data array\n");
		// failB=true;
		// failS="rlft3: problem with dimensions or contiguity of data array";
		// return;
		// }
		// data is [n1][n2][n3] or dataN=[n1n2n3]

		double[] dataN = new double[nn1 * nn2 * nn3];// @@
		for (int i = 1; i <= nn1; i++) {
			for (int j = 1; j <= nn2; j++) {
				for (int k = 1; k <= nn3; k++) {
					dataN[(i - 1) * nn2 * nn3 + (j - 1) * nn3 + k - 1] = data[i - 1][j - 1][k - 1];
				}
			}
		}
		// For a two-dimensional array, this is equivalent to storing the array
		// by rows.!!!
		// /===========================================
		c1 = 0.5;
		c2 = -0.5 * isign;
		theta = isign * (6.28318530717959 / nn3);
		wtemp = Math.sin(0.5 * theta);
		wpr = -2.0 * wtemp * wtemp;
		wpi = Math.sin(theta);
		nn[0] = nn1;// nn[1]=nn1;
		nn[1] = nn2;// nn[2]=nn2;
		nn[2] = nn3 >> 1;// nn[3]=nn3 >> 1;IMPARTIT LA 2 e ok!! data is corect
							// definit!!data is [n1][n2][n3]
		if (isign == 1) {// Case of forward transform.
							// ========================
			fourn(dataN, nn, 3, isign);
			for (int i = 1; i <= nn1; i++) {
				for (int j = 1; j <= nn2; j++) {
					for (int k = 1; k <= nn3; k++) {
						data[i - 1][j - 1][k - 1] = dataN[(i - 1) * nn2 * nn3
								+ (j - 1) * nn3 + k - 1];
					}
				}
			}
			// ========================
			// fourn(&data[1][1][1]-1,nn,3,isign);
			// //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			// &h[j-K]= pointer to the j-K element of vector=> it was taken the
			// elements from j-K to the end!!!
			// Here is where most all of the compute time is spent.
			for (i1 = 1; i1 <= nn1; i1++)
				for (i2 = 1, j2 = 0; i2 <= nn2; i2++) {// Extend data
														// periodically into
														// speq.
					speq[i1 - 1][++j2 - 1] = data[i1 - 1][i2 - 1][0];// speq[i1][++j2]=data[i1][i2][1];
					speq[i1 - 1][++j2 - 1] = data[i1 - 1][i2 - 1][1];// speq[i1][++j2]=data[i1][i2][2];
				}
		}
		for (i1 = 1; i1 <= nn1; i1++) {
			j1 = (i1 != 1 ? nn1 - i1 + 2 : 1);
			// Zero frequency is its own reflection, otherwise locate
			// corresponding negative frequency
			// in wrap-around order.
			wr = 1.0; // Initialize trigonometric recurrence.
			wi = 0.0;
			for (ii3 = 1, i3 = 1; i3 <= (nn3 >> 2) + 1; i3++, ii3 += 2) {
				for (i2 = 1; i2 <= nn2; i2++) {
					if (i3 == 1) {// Equation (12.3.5).
						j2 = (i2 != 1 ? ((nn2 - i2) << 1) + 3 : 1);
						// h1r=c1*(data[i1][i2][1]+speq[j1][j2]);
						h1r = c1
								* (data[i1 - 1][i2 - 1][0] + speq[j1 - 1][j2 - 1]);
						// h1i=c1*(data[i1][i2][2]-speq[j1][j2+1]);
						h1i = c1 * (data[i1 - 1][i2 - 1][1] - speq[j1 - 1][j2]);
						// h2i=c2*(data[i1][i2][1]-speq[j1][j2]);
						h2i = c2
								* (data[i1 - 1][i2 - 1][0] - speq[j1 - 1][j2 - 1]);
						// h2r= -c2*(data[i1][i2][2]+speq[j1][j2+1]);
						h2r = -c2
								* (data[i1 - 1][i2 - 1][1] + speq[j1 - 1][j2]);
						data[i1 - 1][i2 - 1][0] = h1r + h2r;// data[i1][i2][1]=h1r+h2r;
						data[i1 - 1][i2 - 1][1] = h1i + h2i;// data[i1][i2][2]=h1i+h2i;
						speq[j1 - 1][j2 - 1] = h1r - h2r;// speq[j1][j2]=h1r-h2r;
						speq[j1 - 1][j2] = h2i - h1i;// speq[j1][j2+1]=h2i-h1i;
					} else {
						j2 = (i2 != 1 ? nn2 - i2 + 2 : 1);
						j3 = nn3 + 3 - (i3 << 1);
						// h1r=c1*(data[i1][i2][ii3]+data[j1][j2][j3]);
						h1r = c1
								* (data[i1 - 1][i2 - 1][ii3 - 1] + data[j1 - 1][j2 - 1][j3 - 1]);
						// h1i=c1*(data[i1][i2][ii3+1]-data[j1][j2][j3+1]);
						h1i = c1
								* (data[i1 - 1][i2 - 1][ii3] - data[j1 - 1][j2 - 1][j3]);
						// h2i=c2*(data[i1][i2][ii3]-data[j1][j2][j3]);
						h2i = c2
								* (data[i1 - 1][i2 - 1][ii3 - 1] - data[j1 - 1][j2 - 1][j3 - 1]);
						// h2r= -c2*(data[i1][i2][ii3+1]+data[j1][j2][j3+1]);
						h2r = -c2
								* (data[i1 - 1][i2 - 1][ii3] + data[j1 - 1][j2 - 1][j3]);
						// data[i1][i2][ii3]=h1r+wr*h2r-wi*h2i;
						data[i1 - 1][i2 - 1][ii3 - 1] = h1r + wr * h2r - wi
								* h2i;
						// data[i1][i2][ii3+1]=h1i+wr*h2i+wi*h2r;
						data[i1 - 1][i2 - 1][ii3] = h1i + wr * h2i + wi * h2r;
						// data[j1][j2][j3]=h1r-wr*h2r+wi*h2i;
						data[j1 - 1][j2 - 1][j3 - 1] = h1r - wr * h2r + wi
								* h2i;
						// data[j1][j2][j3+1]= -h1i+wr*h2i+wi*h2r;
						data[j1 - 1][j2 - 1][j3] = -h1i + wr * h2i + wi * h2r;
					}
				}
				wr = (wtemp = wr) * wpr - wi * wpi + wr; // Do the recurrence.
				wi = wi * wpr + wtemp * wpi + wi;
			}
		}
		if (isign == -1)// Case of reverse transform.
		{
			fourn(dataN, nn, 3, isign);
			for (int i = 1; i <= nn1; i++) {
				for (int j = 1; j <= nn2; j++) {
					for (int k = 1; k <= nn3; k++) {
						data[i - 1][j - 1][k - 1] = dataN[(i - 1) * nn2 * nn3
								+ (j - 1) * nn3 + k - 1];
					}
				}
			}
			// fourn(&data[1][1][1]-1,nn,3,isign);//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		}
	}

	/*
	 * The first program fragment FFTs a two-dimensional data array, allows for
	 * some processing on it, e.g., filtering, and then takes the inverse
	 * transform. Figure 12.5.2 shows an example of the use of this kind of
	 * code: A sharp image becomes blurry when its high-frequency spatial
	 * components are suppressed by the factor (here) max (1 - 6f^2/fc^2, 0).
	 * The second program example illustrates a three-dimensional transform,
	 * where the three dimensions have different lengths. The third program
	 * example is an example of convolution, as it might occur in a program to
	 * compute the potential generated by a three-dimensional distribution of
	 * sources. NOTA:
	 * 
	 * float ***f3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long
	 * ndh) // allocate a float 3tensor with range
	 * t[nrl..nrh][ncl..nch][ndl..ndh]
	 * 
	 * #define N2 256 #define N3 256 Note that the first component must be set
	 * to 1. int main(void) // example1 // This fragment shows how one might
	 * filter a 256 by 256 digital image. { void rlft3(float ***data, float
	 * **speq, unsigned long nn1, unsigned long nn2, unsigned long nn3, int
	 * isign); float ***data, **speq; data=f3tensor(1,1,1,N2,1,N3);
	 * speq=matrix(1,1,1,2*N2); // ...// Here the image would be loaded into
	 * data. rlft3(data,speq,1,N2,N3,1); // ...// Here the arrays data and speq
	 * would be multiplied by a suitable filter function (of frequency).
	 * rlft3(data,speq,1,N2,N3,-1); // ...// Here the filtered image would be
	 * unloaded from data. free_matrix(speq,1,1,1,2*N2);
	 * free_f3tensor(data,1,1,1,N2,1,N3); return 0; }
	 * 
	 * 
	 * #define N1 32 #define N2 64 #define N3 16 int main(void) // example2 //
	 * This fragment shows how one might FFT a real three-dimensional array of
	 * size 32 by 64 by 16. { void rlft3(float ***data, float **speq, unsigned
	 * long nn1, unsigned long nn2, unsigned long nn3, int isign); int j; float
	 * ***data,**speq; data=f3tensor(1,N1,1,N2,1,N3); speq=matrix(1,N1,1,2*N2);
	 * // ...// Here load data. rlft3(data,speq,N1,N2,N3,1); // ...// Here
	 * unload data and speq. free_matrix(speq,1,N1,1,2*N2);
	 * free_f3tensor(data,1,N1,1,N2,1,N3); return 0; }
	 * 
	 * 
	 * #define N 32 int main(void) // example3 // This fragment shows how one
	 * might convolve two real, three-dimensional arrays of size 32 by 32 by 32,
	 * replacing the first array by the result. { void rlft3(float ***data,
	 * float **speq, unsigned long nn1, unsigned long nn2, unsigned long nn3,
	 * int isign); int j; float
	 * fac,r,i,***data1,***data2,**speq1,**speq2,*sp1,*sp2;
	 * data1=f3tensor(1,N,1,N,1,N); data2=f3tensor(1,N,1,N,1,N);
	 * speq1=matrix(1,N,1,2*N); speq2=matrix(1,N,1,2*N); // ...//
	 * rlft3(data1,speq1,N,N,N,1); FFT both input arrays.
	 * rlft3(data2,speq2,N,N,N,1); fac=2.0/(N*N*N); Factor needed to get
	 * normalized inverse. sp1 = &data1[1][1][1]; sp2 = &data2[1][1][1]; for
	 * (j=1;j<=N*N*N/2;j++) { Note how this can be made a single for-loop
	 * instead of three nested ones by using the pointers sp1 and sp2. r =
	 * sp1[0]*sp2[0] - sp1[1]*sp2[1]; i = sp1[0]*sp2[1] + sp1[1]*sp2[0]; sp1[0]
	 * = fac*r; sp1[1] = fac*i; sp1 += 2; sp2 += 2; } sp1 = &speq1[1][1]; sp2 =
	 * &speq2[1][1]; for (j=1;j<=N*N;j++) { r = sp1[0]*sp2[0] - sp1[1]*sp2[1]; i
	 * = sp1[0]*sp2[1] + sp1[1]*sp2[0]; sp1[0] = fac*r; sp1[1] = fac*i; sp1 +=
	 * 2; sp2 += 2; } rlft3(data1,speq1,N,N,N,-1); Inverse FFT the product of
	 * the two FFTs. // ...// free_matrix(speq2,1,N,1,2*N);
	 * free_matrix(speq1,1,N,1,2*N); free_f3tensor(data2,1,N,1,N,1,N);
	 * free_f3tensor(data1,1,N,1,N,1,N); return 0; }
	 * 
	 * 
	 * //======================================== External Storage or
	 * Memory-Local FFTs Sometime in your life, you might have to compute the
	 * Fourier transform of a really large data set, larger than the size of
	 * your computer’s physical memory. In such a case, the data will be stored
	 * on some external medium, such as magnetic or optical tape or disk. Needed
	 * is an algorithm that makes some manageable number of sequential passes
	 * through the external data, processing it on the fly and outputting
	 * intermediate results to other external media, which can be read on
	 * subsequent passes. In fact, an algorithm of just this description was
	 * developed by Singleton [1] very soon after the discovery of the FFT. The
	 * algorithm requires four sequential storage devices, each capable of
	 * holding half of the input data. The first half of the input data is
	 * initially on one device, the second half on another. Singleton’s
	 * algorithm is based on the observation that it is possible to bit-reverse
	 * 2M values by the following sequence of operations: On the first pass,
	 * values are read alternately from the two input devices, and written to a
	 * single output device (until it holds half the data), and then to the
	 * other output device. On the second pass, the output devices become input
	 * devices, and vice versa. Now, we copy two values from the first device,
	 * then two values from the second, writing them (as before) first to fill
	 * one output device, then to fill a second. Subsequent passes read 4, 8,
	 * etc., input values at a time. After completion of pass M - 1, the data
	 * are in bit-reverse order. Singleton’s next observation is that it is
	 * possible to alternate the passes of essentially this bit-reversal
	 * technique with passes that implement one stage of the Danielson-Lanczos
	 * combination formula (12.2.3). The scheme, roughly, is this: One starts as
	 * before with half the input data on one device, half on another. In the
	 * first pass, one complex value is read from each input device. Two
	 * combinations are formed, and one is written to each of two output
	 * devices. After this “computing” pass, the devices are rewound, and a
	 * “permutation” pass is performed, where groups of values are read from the
	 * first input device and alternately written to the first and second output
	 * devices; when the first input device is exhausted, the second is
	 * similarly processed. This sequence of computing and permutation passes is
	 * repeated M - K - 1 times, where 2K is the size of internal buffer
	 * available to the program. The second phase of the computation consists of
	 * a finalK computation passes. What distinguishes the second phase from the
	 * first is that, now, the permutations are local enough to do in place
	 * during the computation. There are thus no separate permutation passes in
	 * the second phase. In all, there are 2M - K - 2 passes through the data.
	 * Here is an implementation of Singleton’s algorithm, based on [1]:
	 * 
	 * 
	 * #define KBF 128 void fourfs(FILE *file[5], unsigned long nn[], int ndim,
	 * int isign) One- or multi-dimensional Fourier transform of a large data
	 * set stored on external media. On input, ndim is the number of dimensions,
	 * and nn[1..ndim] contains the lengths of each dimension (number of real
	 * and imaginary value pairs), which must be powers of two. file[1..4]
	 * contains the stream pointers to 4 temporary files, each large enough to
	 * hold half of the data. The four streams must be opened in the system’s
	 * “binary” (as opposed to “text”) mode. The input data must be in C normal
	 * order, with its first half stored in file file[1], its second half in
	 * file[2], in native floating point form. KBF real numbers are processed
	 * per buffered read or write. isign should be set to 1 for the Fourier
	 * transform, to -1 for its inverse. On output, values in the array file may
	 * have been permuted; the first half of the result is stored in file[3],
	 * the second half in file[4]. N.B.: For ndim > 1, the output is stored by
	 * columns, i.e., not in C normal order; in other words, the output is the
	 * transpose of that which would have been produced by routine fourn. { void
	 * fourew(FILE *file[5], int *na, int *nb, int *nc, int *nd); unsigned long
	 * j,j12,jk,k,kk,n=1,mm,kc=0,kd,ks,kr,nr,ns,nv; int cc,na,nb,nc,nd; float
	 * tempr,tempi,*afa,*afb,*afc; double wr,wi,wpr,wpi,wtemp,theta; static int
	 * mate[5] = {0,2,1,4,3}; afa=vector(1,KBF); afb=vector(1,KBF);
	 * afc=vector(1,KBF); for (j=1;j<=ndim;j++) { n *= nn[j]; if (nn[j] <= 1)
	 * nrerror("invalid float or wrong ndim in fourfs"); } nv=1; jk=nn[nv];
	 * mm=n; ns=n/KBF; nr=ns >> 1; kd=KBF >> 1; ks=n;
	 * fourew(file,&na,&nb,&nc,&nd); The first phase of the transform starts
	 * here. for (;;) { Start of the computing pass.
	 * theta=isign*3.141592653589793/(n/mm); wtemp=sin(0.5*theta); wpr =
	 * -2.0*wtemp*wtemp; wpi=sin(theta); wr=1.0; wi=0.0; mm >>= 1; for
	 * (j12=1;j12<=2;j12++) { kr=0; do {
	 * cc=fread(&afa[1],sizeof(float),KBF,file[na]); if (cc != KBF)
	 * nrerror("read error in fourfs");
	 * cc=fread(&afb[1],sizeof(float),KBF,file[nb]); if (cc != KBF)
	 * nrerror("read error in fourfs"); for (j=1;j<=KBF;j+=2) {
	 * tempr=((float)wr)*afb[j]-((float)wi)*afb[j+1];
	 * tempi=((float)wi)*afb[j]+((float)wr)*afb[j+1]; afb[j]=afa[j]-tempr;
	 * afa[j] += tempr; afb[j+1]=afa[j+1]-tempi; afa[j+1] += tempi; } kc += kd;
	 * if (kc == mm) { kc=0; wr=(wtemp=wr)*wpr-wi*wpi+wr;
	 * wi=wi*wpr+wtemp*wpi+wi; } cc=fwrite(&afa[1],sizeof(float),KBF,file[nc]);
	 * if (cc != KBF) nrerror("write error in fourfs");
	 * cc=fwrite(&afb[1],sizeof(float),KBF,file[nd]); if (cc != KBF)
	 * nrerror("write error in fourfs"); } while (++kr < nr); if (j12 == 1 && ks
	 * != n && ks == KBF) { na=mate[na]; nb=na; } if (nr == 0) break; }
	 * fourew(file,&na,&nb,&nc,&nd); Start of the permutation pass. jk >>= 1;
	 * while (jk == 1) { mm=n; jk=nn[++nv]; } ks >>= 1; if (ks > KBF) { for
	 * (j12=1;j12<=2;j12++) { for (kr=1;kr<=ns;kr+=ks/KBF) { for
	 * (k=1;k<=ks;k+=KBF) { cc=fread(&afa[1],sizeof(float),KBF,file[na]); if (cc
	 * != KBF) nrerror("read error in fourfs");
	 * cc=fwrite(&afa[1],sizeof(float),KBF,file[nc]); if (cc != KBF)
	 * nrerror("write error in fourfs"); } nc=mate[nc]; } na=mate[na]; }
	 * fourew(file,&na,&nb,&nc,&nd); } else if (ks == KBF) nb=na; else break; }
	 * j=1; The second phase of the transform starts here. Now, the remaining
	 * permutations are suf- ficiently local to be done in place. for (;;) {
	 * theta=isign*3.141592653589793/(n/mm); wtemp=sin(0.5*theta); wpr =
	 * -2.0*wtemp*wtemp; wpi=sin(theta); wr=1.0; wi=0.0; mm >>= 1; ks=kd; kd >>=
	 * 1; for (j12=1;j12<=2;j12++) { for (kr=1;kr<=ns;kr++) {
	 * cc=fread(&afc[1],sizeof(float),KBF,file[na]); if (cc != KBF)
	 * nrerror("read error in fourfs"); kk=1; k=ks+1; for (;;) {
	 * tempr=((float)wr)*afc[kk+ks]-((float)wi)*afc[kk+ks+1];
	 * tempi=((float)wi)*afc[kk+ks]+((float)wr)*afc[kk+ks+1];
	 * afa[j]=afc[kk]+tempr; afb[j]=afc[kk]-tempr; afa[++j]=afc[++kk]+tempi;
	 * afb[j++]=afc[kk++]-tempi; if (kk < k) continue; kc += kd; if (kc == mm) {
	 * kc=0; wr=(wtemp=wr)*wpr-wi*wpi+wr; wi=wi*wpr+wtemp*wpi+wi; } kk += ks; if
	 * (kk > KBF) break; else k=kk+ks; } if (j > KBF) {
	 * cc=fwrite(&afa[1],sizeof(float),KBF,file[nc]); if (cc != KBF)
	 * nrerror("write error in fourfs");
	 * cc=fwrite(&afb[1],sizeof(float),KBF,file[nd]); if (cc != KBF)
	 * nrerror("write error in fourfs"); j=1; } } na=mate[na]; }
	 * fourew(file,&na,&nb,&nc,&nd); jk >>= 1; if (jk > 1) continue; mm=n; do {
	 * if (nv < ndim) jk=nn[++nv]; else { free_vector(afc,1,KBF);
	 * free_vector(afb,1,KBF); free_vector(afa,1,KBF); return; } } while (jk ==
	 * 1); } }
	 * 
	 * 
	 * 
	 * 
	 * 
	 * 
	 * 
	 * 
	 * #define SWAP(a,b) ftemp=(a);(a)=(b);(b)=ftemp void fourew(FILE *file[5],
	 * int *na, int *nb, int *nc, int *nd) Utility used by fourfs. Rewinds and
	 * renumbers the four files. { int i; FILE *ftemp; for (i=1;i<=4;i++)
	 * rewind(file[i]); SWAP(file[2],file[4]); SWAP(file[1],file[3]);na=3;nb=4;
	 * nc=1;nd=2; }
	 */
}
