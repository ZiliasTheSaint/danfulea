package danfulea.math.numerical;

/**
 * Random number generators.
 * Based on literature, including Numerical Recipes in C (Cambridge Univ.) and EGSnrc.
 * @author Dan Fulea, 09 OCT. 2006
 */
public class RandomNrs {
	public static Function func;
	public static boolean failB = false;
	public static String failS = "";

	public static int IA = 16807;
	public static int IM = 2147483647;
	public static double AM = (1.0 / IM);
	public static int IQ = 127773;
	public static int IR = 2836;
	public static int MASK = 123459876;

	public static int NTAB = 32;
	public static int NDIV = (1 + (IM - 1) / NTAB);
	public static double EPS = 1.2e-7;
	public static double RNMX = (1.0 - EPS);
	public static long iy = 0;// int iy=0;
	public static long[] iv = new long[NTAB];// int[] iv=new int[NTAB];

	public static int IM1 = 2147483563;
	public static int IM2 = 2147483399;
	public static double AM1 = (1.0 / IM1);// AM (1.0/IM1)
	public static int IMM1 = (IM1 - 1);
	public static int IA1 = 40014;
	public static int IA2 = 40692;
	public static int IQ1 = 53668;
	public static int IQ2 = 52774;
	public static int IR1 = 12211;
	public static int IR2 = 3791;
	// #define NTAB 32
	public static int NDIV1 = (1 + IMM1 / NTAB);
	// public static double EPS =1.2e-7;
	// public static double RNMX =(1.0-EPS);

	public static long idum2 = 123456789;// int idum2=123456789;
	public static long iy1 = 0;// int iy1=0;
	public static long[] iv1 = new long[NTAB];// int[] iv1=new int[NTAB];

	public static int MBIG = 1000000000;
	public static int MSEED = 161803398;
	public static int MZ = 0;
	public static double FAC = (1.0 / MBIG);

	public static int inext = 0;
	public static int inextp = 0;
	public static long[] ma = new long[56]; // The value 56 (range ma[1..55]) is
											// special and
	// should not be modified; see Knuth.
	public static int iff = 0;

	public static long idum_ran0 = 0;
	public static long idum_ran1 = 0;
	public static long idum_ran2 = 0;
	public static long idum_ran3 = 0;

	public static boolean ranluxB = true;// use ranlux if true, ranmar otherwise
	public static int $NRANMAR = 128;
	public static double[] rng_array1 = new double[$NRANMAR];

	public static int rng_seed = 999999;// "current pointer for rng_array            "
	public static double[] rng_array = new double[24];// "containes 24 random numbers              "
	public static int[] seeds = new int[24];// "for storing the rng state                "
	public static int seedin = 0;
	public static int luxury_level = 0;
	public static int[] state = new int[25];
	public static int carry = 0;// " The state of the generator "
	public static int i24 = 0;
	public static int j24 = 0;// " The rng seeds "
	public static int[] next = new int[24];// " for convinience "
	public static int jseed_dflt = 314159265;// 0;
	public static int nskipRnd = 0;
	public static int icon = 2147483563;// 0;
	public static int status = 0;
	public static int jseed = 0;
	public static int[] nskipll = { 0, 24, 73, 199, 365 };// new int[5];//(0:4),
	public static int icarry = 0;
	public static int kRnd = 0;
	public static int jRnd = 0;
	public static boolean not_initialized = true;
	public static int uni = 0;
	public static double twom24 = 0.0;
	public static double twop24 = 0.0;
	public static int $DEFAULT_LL = 1;// REPLACE {$DEFAULT-LL} WITH {1}
	public static int $MXRNGDIM = 5;// REPLACE {$MXRNGDIM} WITH {5}
	// ranmar
	public static int[] urndm = new int[97];
	public static int crndm = 0;
	public static int cdrndm = 0;
	public static int cmrndm = 0;
	public static int i4opt = 0;
	public static int ixx = 0;
	public static int jxx = 0;
	public static int fool_optimizer = 0;

	public static int iset = 0;
	public static double gset = 0.0;

	public static double PI = 3.141592654;
	public static double sq_poidev = 0.0;
	public static double alxm_poidev = 0.0;
	public static double g_poidev = 0.0;
	public static double oldm_poidev = (-1.0); // oldm is a flag for whether xm
												// has changed
	// since last call.

	public static int nold_bnldev = (-1);
	public static double pold_bnldev = (-1.0);
	public static double pc_bnldev = 0.0;
	public static double plog_bnldev = 0.0;
	public static double pclog_bnldev = 0.0;
	public static double en_bnldev = 0.0;
	public static double oldg_bnldev = 0.0;

	public static double ALPH = 1.5;
	public static int NDMX = 50;
	public static int MXDIM = 10;
	public static double TINY = 1.0e-30;
	public static double tgral_vegas = 0.0;
	public static double sd_vegas = 0.0;
	public static double chi2a_vegas = 0.0;
	public static int i_vegas = 0;
	public static int it_vegas = 0;
	public static int j_vegas = 0;
	public static int k_vegas = 0;
	public static int mds_vegas = 0;
	public static int nd_vegas = 0;
	public static int ndo_vegas = 0;
	public static int ng_vegas = 0;
	public static int npg_vegas = 0;
	public static int[] ia_vegas = new int[MXDIM + 1];
	public static int[] kg_vegas = new int[MXDIM + 1];
	public static double calls_vegas = 0.0;
	public static double dv2g_vegas = 0.0;
	public static double dxg_vegas = 0.0;
	public static double f_vegas = 0.0;
	public static double f2_vegas = 0.0;
	public static double f2b_vegas = 0.0;
	public static double fb_vegas = 0.0;
	public static double rc_vegas = 0.0;
	public static double ti_vegas = 0.0;
	public static double tsi_vegas = 0.0;
	public static double wgt_vegas = 0.0;
	public static double xjac_vegas = 0.0;
	public static double xn_vegas = 0.0;
	public static double xnd_vegas = 0.0;
	public static double xo_vegas = 0.0;
	public static double[][] d_vegas = new double[NDMX + 1][MXDIM + 1];
	public static double[][] di_vegas = new double[NDMX + 1][MXDIM + 1];
	public static double[] dt_vegas = new double[MXDIM + 1];
	public static double[] dx_vegas = new double[MXDIM + 1];
	public static double[] r_vegas = new double[NDMX + 1];
	public static double[] x_vegas = new double[MXDIM + 1];
	public static double[][] xi_vegas = new double[MXDIM + 1][NDMX + 1];
	public static double[] xin_vegas = new double[NDMX + 1];
	public static double schi_vegas = 0.0;
	public static double si_vegas = 0.0;
	public static double swgt_vegas = 0.0;

	public static double PFAC = 0.1;
	public static int MNPT = 15;
	public static int MNBS = 60;
	// public static double TINY =1.0e-30;
	public static double BIG = 1.0e30;
	public static int iran = 0;
	public static double ave_miser = 0.0;
	public static double var_miser = 0.0;

	/*
	 * First proposed by Lewis, Goodman, and Miller in 1969, this generator has
	 * in subsequent years passed all new theoretical tests, and (perhaps more
	 * importantly) has accumulated a large amount of successful use. Park and
	 * Miller do not claim that the generator is “perfect” (we will see below
	 * that it is not), but only that it is a good minimal standard against
	 * which other generators should be judged.
	 * 
	 * The period of ran0 is 231 - 2 ? 2.1 × 109. A peculiarity of generators of
	 * the form (7.1.2) is that the value 0 must never be allowed as the initial
	 * seed — it perpetuates itself — and it never occurs for any nonzero
	 * initial seed. Experience has shown that users always manage to call
	 * random number generators with the seed idum=0. That is why ran0 performs
	 * its exclusive-or with an arbitrary constant both on entry and exit. If
	 * you are the first user in history to be proof against human error, you
	 * can remove the two lines with the ^ operation.
	 * 
	 * The routine ran0 is a Minimal Standard, satisfactory for the majority of
	 * applications, but we do not recommend it as the final word on random
	 * number generators. Our reason is precisely the simplicity of the Minimal
	 * Standard. It is not hard to think of situations where successive random
	 * numbers might be used in a way that accidentally conflicts with the
	 * generation algorithm. For example, since successive numbers differ by a
	 * multiple of only 1.6 × 104 out of a modulus of more than 2 × 109, very
	 * small random numbers will tend to be followed by smaller than average
	 * values. One time in 106, for example, there will be a value < 10-6
	 * returned (as there should be), but this will always be followed by a
	 * value less than about 0.0168. One can easily think of applications
	 * involving rare events where this property would lead to wrong results.
	 */
	/**
	 * “Minimal” random number generator of Park and Miller. Returns a uniform random deviate between 0.0 and 1.0.
	 * Set or reset idum to any integer value (except the unlikely value MASK) to initialize the sequence; idum must not be altered between calls for 
	 * successive deviates in a sequence.
	 * @return random number
	 */
	public static double ran0()// int idum)//(long *idum)
	// “Minimal” random number generator of Park and Miller. Returns a uniform
	// random deviate
	// between 0.0 and 1.0. Set or reset idum to any integer value (except the
	// unlikely value MASK)
	// to initialize the sequence; idum must not be altered between calls for
	// successive deviates in
	// a sequence.
	{
		// idum_ran0=idum;
		// int k=0;
		long k = 0;
		double ans = 0.0;
		// *idum ^= MASK; XORing with MASK allows use of zero and other
		// simple bit patterns for idum.
		idum_ran0 ^= MASK;
		k = (idum_ran0) / IQ;
		idum_ran0 = IA * (idum_ran0 - k * IQ) - IR * k; // Compute
														// idum=(IA*idum) % IM
														// without over-
		// flows by Schrage’s method.
		if (idum_ran0 < 0)
			idum_ran0 += IM;
		ans = AM * (idum_ran0);// Convert idum to a floating result.
		idum_ran0 ^= MASK; // Unmask before return.
		return ans;
	}

	/*
	 * The following routine, ran1, uses the Minimal Standard for its random
	 * value, but it shuffles the output to remove low-order serial
	 * correlations. The routine ran1 passes those statistical tests that ran0
	 * is known to fail. In fact, we do not know of any statistical test that
	 * ran1 fails to pass, except when the number of calls starts to become on
	 * the order of the period m, say > 10 ^8 ~ m/20.
	 */
	/**
	 * “Minimal” random number generator of Park and Miller with Bays-Durham shuffle and added safeguards. Returns a uniform random deviate between 0.0 and 1.0 
	 * (exclusive of the endpoint values). Call with idum a negative integer to initialize; thereafter, do not alter idum between 
	 * successive deviates in a sequence. RNMX should approximate the largest floating value that is less than 1.
	 * @return the random number
	 */
	public static double ran1()// long *idum)
	// “Minimal” random number generator of Park and Miller with Bays-Durham
	// shuffle and added
	// safeguards. Returns a uniform random deviate between 0.0 and 1.0
	// (exclusive of the endpoint
	// values). Call with idum a negative integer to initialize; thereafter, do
	// not alter idum between
	// successive deviates in a sequence. RNMX should approximate the largest
	// floating value that is
	// less than 1.
	{
		int j = 0;
		// int k=0;
		long k = 0;
		// static long iy=0;
		// static long iv[NTAB];
		double temp = 0.0;
		if (idum_ran1 <= 0 || iy == 0)// !iy)
		{// Initialize.
			if (-(idum_ran1) < 1)
				idum_ran1 = 1; // Be sure to prevent idum = 0.
			else
				idum_ran1 = -(idum_ran1);
			for (j = NTAB + 7; j >= 0; j--) {// Load the shuffle table (after 8
												// warm-ups).
				k = (idum_ran1) / IQ;
				idum_ran1 = IA * (idum_ran1 - k * IQ) - IR * k;
				if (idum_ran1 < 0)
					idum_ran1 += IM;
				if (j < NTAB)
					iv[j] = idum_ran1;
			}
			iy = iv[0];
		}
		k = (idum_ran1) / IQ;// Start here when not initializing.
		idum_ran1 = IA * (idum_ran1 - k * IQ) - IR * k;// Compute idum=(IA*idum)
														// % IM without over-
		// flows by Schrage’s method.
		if (idum_ran1 < 0)
			idum_ran1 += IM;
		j = (int) (iy / NDIV);// Will be in the range 0..NTAB-1.
		iy = iv[j];// Output previously stored value and refill the
					// shuffle table.
		iv[j] = idum_ran1;
		if ((temp = AM * iy) > RNMX)
			return RNMX; // Because users don’t expect endpoint values.
		else
			return temp;
	}

	/*
	 * For situations when even longer random sequences are needed, L’Ecuyer [6]
	 * has given a good way of combining two different sequences with different
	 * periods so as to obtain a new sequence whose period is the least common
	 * multiple of the two periods. Combining the two generators breaks up
	 * serial correlations to a considerable extent. We nevertheless recommend
	 * the additional shuffle that is implemented in the following routine,
	 * ran2. We think that, within the limits of its floating-point precision,
	 * ran2 provides perfect random numbers; a practical definition of “perfect”
	 * is that we will pay $1000 to the first reader who convinces us otherwise
	 * (by finding a statistical test that ran2 fails in a nontrivial way,
	 * excluding the ordinary limitations of a machine’s floating-point
	 * representation).
	 * 
	 * On balance, we recommend ran1 for general use. It is portable, based on
	 * Park and Miller’s Minimal Standard generator with an additional shuffle,
	 * and has no known (to us) flaws other than period exhaustion.=> ==> If you
	 * are generating more than 100,000,000 random numbers in a single
	 * calculation (that is, more than about 5% of ran1’s period), we recommend
	 * the use of ran2, with its much longer period.
	 */
	// @@@@VERY GOOD they say!!!
	/**
	 * Long period (greater than 2 × 1018) random number generator of L’Ecuyer with Bays-Durham shuffle and added safeguards. 
	 *  Returns a uniform random deviate between 0.0 and 1.0 (exclusive of the endpoint values). Call with idum a negative integer to initialize; 
	 *  thereafter, do not alter idum between successive deviates in a sequence. RNMX should approximate the largest floating value that is less than 1.
	 * @return the random number
	 */
	public static double ran2()// long *idum)
	// Long period (> 2 × 1018) random number generator of L’Ecuyer with
	// Bays-Durham shuffle
	// and added safeguards. Returns a uniform random deviate between 0.0 and
	// 1.0 (exclusive of
	// the endpoint values). Call with idum a negative integer to initialize;
	// thereafter, do not alter
	// idum between successive deviates in a sequence. RNMX should approximate
	// the largest floating
	// value that is less than 1.
	{
		int j = 0;
		// int k=0;
		long k = 0;
		double temp = 0.0;
		if (idum_ran2 <= 0) {// Initialize.
			if (-(idum_ran2) < 1)
				idum_ran2 = 1;// Be sure to prevent idum = 0.
			else
				idum_ran2 = -(idum_ran2);
			idum2 = (idum_ran2);
			for (j = NTAB + 7; j >= 0; j--) {// Load the shuffle table (after 8
												// warm-ups).
				k = (idum_ran2) / IQ1;
				idum_ran2 = IA1 * (idum_ran2 - k * IQ1) - k * IR1;
				if (idum_ran2 < 0)
					idum_ran2 += IM1;
				if (j < NTAB)
					iv1[j] = idum_ran2;
			}
			iy1 = iv1[0];
		}
		k = (idum_ran2) / IQ1; // Start here when not initializing.
		idum_ran2 = IA1 * (idum_ran2 - k * IQ1) - k * IR1;// Compute
															// idum=(IA1*idum) %
															// IM1 without
		// overflows by Schrage’s method.
		if (idum_ran2 < 0)
			idum_ran2 += IM1;
		k = idum2 / IQ2;
		idum2 = IA2 * (idum2 - k * IQ2) - k * IR2; // Compute idum2=(IA2*idum) %
													// IM2 likewise.
		if (idum2 < 0)
			idum2 += IM2;
		j = (int) (iy1 / NDIV1);// Will be in the range 0..NTAB-1.
		iy1 = iv1[j] - idum2; // Here idum is shuffled, idum and idum2 are
		// combined to generate output.
		iv1[j] = idum_ran2;
		if (iy1 < 1)
			iy1 += IMM1;
		if ((temp = AM1 * iy1) > RNMX)
			return RNMX;// Because users don’t expect endpoint values.
		else
			return temp;
	}

	/*
	 * Finally, we give you Knuth’s suggestion [4] for a portable routine, which
	 * we have translated to the present conventions as ran3. This is not based
	 * on the linear congruential method at all, but rather on a subtractive
	 * method (see also [5]). One might hope that its weaknesses, if any, are
	 * therefore of a highly different character from the weaknesses, if any, of
	 * ran1 above. If you ever suspect trouble with one routine, it is a good
	 * idea to try the other in the same application.
	 * 
	 * Knuth’s subtractive routine ran3 seems to be the timing winner among
	 * portable routines. Unfortunately the subtractive method is not so well
	 * studied, and not a standard. We like to keep ran3 in reserve for a
	 * “second opinion,” substituting it when we suspect another generator of
	 * introducing unwanted correlations into a calculation.
	 */
	/**
	 * Returns a uniform random deviate between 0.0 and 1.0. Set idum to any static long ma[56]; The value 56 (range ma[1..55]) is special and 
	 * should not be modified;
	 * @return the random number
	 */
	public static double ran3()// long *idum)
	// Returns a uniform random deviate between 0.0 and 1.0. Set idum to any
	// negative value to
	// initialize or reinitialize the sequence.
	{
		// static int inext,inextp;
		// static long ma[56]; The value 56 (range ma[1..55]) is special and
		// should not be modified; see Knuth. static int iff=0;
		long mj = 0;
		long mk = 0;
		int i = 0;
		int ii = 0;
		int k = 0;
		if (idum_ran3 < 0 || iff == 0) {// Initialization.
			iff = 1;
			mj = Math.abs(MSEED - Math.abs(idum_ran3)); // Initialize ma[55]
														// using the seed idum
														// and the
			// large number MSEED. mj %= MBIG;
			ma[55] = mj;
			mk = 1;
			for (i = 1; i <= 54; i++) {// Now initialize the rest of the table,
				ii = (21 * i) % 55; // in a slightly random order,
				ma[ii] = mk; // with numbers that are not especially random.
				mk = mj - mk;
				if (mk < MZ)
					mk += MBIG;
				mj = ma[ii];
			}
			for (k = 1; k <= 4; k++)
				// We randomize them by “warming upthe generator.”
				for (i = 1; i <= 55; i++) {
					ma[i] -= ma[1 + (i + 30) % 55];
					if (ma[i] < MZ)
						ma[i] += MBIG;
				}
			inext = 0; // Prepare indices for our first generated number.
			inextp = 31; // The constant 31 is special; see Knuth.
			idum_ran3 = 1;
		}
		// Here is where we start, except on initialization.
		if (++inext == 56)
			inext = 1; // Increment inext and inextp, wrapping around
		// 56 to 1.
		if (++inextp == 56)
			inextp = 1;
		mj = ma[inext] - ma[inextp]; // Generate a new random number
										// subtractively.
		if (mj < MZ)
			mj += MBIG; // Be sure that it is in range.
		ma[inext] = mj;// Store it,
		return mj * FAC; // and output the derived uniform deviate.
	}

	/**
	 * The EGSnrc implementation for RNG (RandomNumberGenerator).
	 * @return the random number
	 */
	public static double RANDOMSET() {
		double result = 0.0;
		if (ranluxB) {
			if (rng_seed > 24) {
				ranlux(rng_array);
				rng_seed = 1;
			}

			result = rng_array[rng_seed - 1];
			rng_seed = rng_seed + 1;
		} else// ranmar
		{
			if (rng_seed > $NRANMAR)
				ranmar_get();
			result = rng_array1[rng_seed - 1] * twom24;
			rng_seed = rng_seed + 1;
		}

		return result;
		// " i.e. take the rng_seed'th random number from rng_array, "
		// " if all 24 numbers used, generate a new set of 24        "
	}

	/**
	 * Internal use by RANDOMSET.
	 */
	public static void ranmar_get() {
		int i = 0;
		int iopt = 0;
		if (rng_seed == 999999)
			init_ranmar();
		for (i = 1; i <= $NRANMAR; i++) {
			iopt = urndm[ixx - 1] - urndm[jxx - 1];
			if (iopt < 0)
				iopt = iopt + 16777216;
			urndm[ixx - 1] = iopt;
			ixx = ixx - 1;
			jxx = jxx - 1;
			if (ixx == 0) {
				ixx = 97;
			} else if (jxx == 0) {
				jxx = 97;
			}
			crndm = crndm - cdrndm;
			if (crndm < 0)
				crndm = crndm + cmrndm;
			iopt = iopt - crndm;
			if (iopt < 0)
				iopt = iopt + 16777216;
			rng_array1[i - 1] = iopt;
		}
		rng_seed = 1;
		return;
	}

	/**
	 * Internal use by ranmar_get.
	 */
	public static void init_ranmar() {
		int s = 0;
		int t = 0;
		int i = 0;
		int j = 0;
		int k = 0;
		int l = 0;
		int m = 0;
		int ii = 0;
		int jj = 0;

		if (ixx <= 0 || ixx > 31328)
			ixx = 1802; // "Sets Marsaglia default"
		if (jxx <= 0 || jxx > 30081)
			jxx = 9373; // "sets Marsaglia default"

		i = ixx / 177;
		i = i % 177;
		i = i + 2;// i = mod(ixx/177,177) + 2;
		j = ixx % 177;
		j = j + 2;// j = mod(ixx, 177) + 2;
		k = jxx / 169;
		k = k % 178;
		k = k + 1;// k = mod(jxx/169,178) + 1;
		l = jxx % 169;// l = mod(jxx, 169) ;

		for (ii = 1; ii <= 97; ii++) {

			s = 0;
			t = 8388608;// "t is 2**23 i.e. half of the maximum allowed"
						// "(note that only 24 bits are used)          "

			for (jj = 1; jj <= 24; jj++) {

				// "The if( fool_optimizer ...) statements below are"
				// "to prevent re-arangement of statements for high "
				// "level optimizations and thus different sequences"
				// "on different architectures                      "
				m = i * j;
				m = m % 179;
				m = m * k;
				m = m % 179;
				// m = mod(mod(i*j,179)*k,179);
				// IF( fool_optimizer = 999 ) [ write(6,*) i,j,k,m,s,t; ]
				i = j;
				// IF( fool_optimizer = 999 ) [ write(6,*) i,j,k,m,s,t; ]
				j = k;
				// IF( fool_optimizer = 999 ) [ write(6,*) i,j,k,m,s,t; ]
				k = m;
				// IF( fool_optimizer = 999 ) [ write(6,*) i,j,k,m,s,t; ]
				l = 53 * l + 1;
				l = l % 169;// l = mod(53*l+1,169);
				// IF( fool_optimizer = 999 ) [ write(6,*) i,j,k,m,s,t; ]
				int lm = l * m;
				lm = lm % 64;
				if (lm >= 32)
					s = s + t;// IF(mod(l*m,64) >= 32) s = s + t;
				// IF( fool_optimizer = 999 ) [ write(6,*) i,j,k,m,s,t; ]
				t = t / 2;
				// IF( fool_optimizer = 999 ) [ write(6,*) i,j,k,m,s,t; ]
			}
			urndm[ii - 1] = s;
		}

		crndm = 362436;
		cdrndm = 7654321;
		cmrndm = 16777213;

		twom24 = 1. / 16777216.;

		ixx = 97;
		jxx = 33;

		rng_seed = $NRANMAR + 1;

		return;
	}

	/**
	 * Internal use by RANDOMSET
	 * @param rng_array rng_array
	 */
	public static void ranlux(double[] rng_array) {

		if (not_initialized) {
			not_initialized = false;
			nskipRnd = nskipll[$DEFAULT_LL];
			twom24 = 1.0;
			twop24 = 1.0;
			jseed = jseed_dflt;
			for (jRnd = 1; jRnd <= 24; jRnd++) {
				twom24 = twom24 * 0.5;
				twop24 = twop24 * 2.0;
				kRnd = jseed / 53668;
				jseed = 40014 * (jseed - kRnd * 53668) - kRnd * 12211;
				if (jseed < 0) {
					jseed = jseed + icon;
				}
				seeds[jRnd - 1] = jseed % 16777216;// mod(jseed,16777216);
				next[jRnd - 1] = jRnd - 1;
			}
			next[0] = 24;// next(1) = 24;
			i24 = 24;
			j24 = 10;
			carry = 0;
			// if( seeds(24) = 0 )
			if (seeds[23] == 0) {
				carry = 1;
			}
		}

		for (jRnd = 1; jRnd <= 24; jRnd++) // j=1,24
		{
			uni = seeds[j24 - 1] - seeds[i24 - 1] - carry;// seeds(j24) -
															// seeds(i24) -
															// carry;
			if (uni < 0) {
				uni = uni + 16777216;
				carry = 1;
			} else {
				carry = 0;
			}
			seeds[i24 - 1] = uni;// seeds(i24) = uni;
			// "IF( uni = 0 ) [ uni = twom24*twom24; ]"
			i24 = next[i24 - 1];// i24 = next(i24);
			j24 = next[j24 - 1];// j24 = next(j24);
			if (uni >= 4096) {
				rng_array[jRnd - 1] = uni * twom24;// rng_array(j) = uni*twom24;
			} else {
				// rng_array(j) = uni*twom24 + seeds(j24)*twom24*twom24;
				rng_array[jRnd - 1] = uni * twom24 + seeds[j24 - 1] * twom24
						* twom24;
			}
		}

		if (nskipRnd > 0) {
			for (jRnd = 1; jRnd <= nskipRnd; jRnd++)// DO jRnd=1,nskipRnd
			{
				// uni = seeds(j24) - seeds(i24) - carry;
				uni = seeds[j24 - 1] - seeds[i24 - 1] - carry;
				if (uni < 0) {
					uni = uni + 16777216;
					carry = 1;
				} else {
					carry = 0;
				}
				seeds[i24 - 1] = uni;// seeds(i24) = uni;
				i24 = next[i24 - 1];// i24 = next[i24];
				j24 = next[j24 - 1];// j24 = next[j24];
			}
		}
		return;
	}

	/**
	 * Initialize ranlux.
	 * @param luxury_level1 luxury_level1
	 * @param seedin1 seedin1
	 */
	public static void init_ranlux(int luxury_level1, int seedin1) {
		seedin = seedin1;// $$$$$$$$$$$
		luxury_level = luxury_level1;// $$$$$$$$$$$

		jseed = seedin;
		if (jseed <= 0)
			jseed = jseed_dflt;
		if ((luxury_level < 0) || (luxury_level > 4)) {
			luxury_level = $DEFAULT_LL;
		}

		nskipRnd = nskipll[luxury_level];

		// OUTPUT luxury_level,jseed;
		// System.out.println(" ** RANLUX initialization **"+
		// " luxury level: "+luxury_level+
		// " initial seed: "+jseed+
		// "*******");

		not_initialized = false;
		twom24 = 1;
		twop24 = 1;
		for (jRnd = 1; jRnd <= 24; jRnd++) {
			twom24 = twom24 * 0.5;
			twop24 = twop24 * 2.0;
			kRnd = jseed / 53668;
			jseed = 40014 * (jseed - kRnd * 53668) - kRnd * 12211;
			if (jseed < 0) {
				jseed = jseed + icon;
			}
			seeds[jRnd - 1] = jseed % 16777216;// mod(jseed,16777216);
			next[jRnd - 1] = jRnd - 1;
		}
		next[0] = 24;
		i24 = 24;
		j24 = 10;
		carry = 0;
		if (seeds[23] == 0) {
			carry = 1;
		}

		return;
	}

	/**
	 * Get ranlux state based on seeds. It is stored in input state variable.
	 * @param state state
	 */
	public static void get_ranlux_state(int[] state) {

		for (jRnd = 1; jRnd <= 24; jRnd++) {
			state[jRnd - 1] = seeds[jRnd - 1];
		}
		state[24] = i24 + 100 * (j24 + 100 * nskipRnd);
		if (carry > 0)
			state[24] = -state[24];
		return;
	}

	/**
	 * Set ranlux status.
	 * @param state state
	 */
	public void set_ranlux_state(int[] state)// public static void
												// set_ranlux_state(int[] state)
	{
		twom24 = 1;
		twop24 = 1;
		for (jRnd = 1; jRnd <= 24; jRnd++) {
			twom24 = twom24 * 0.5;
			twop24 = twop24 * 2;
			next[jRnd - 1] = jRnd - 1;
		}
		next[0] = 24;
		for (jRnd = 1; jRnd <= 24; jRnd++) {
			seeds[jRnd - 1] = state[jRnd - 1];
		}
		if (state[24] <= 0) {
			status = -state[24];
			carry = 1;
		} else {
			status = state[24];
			carry = 0;
		}
		nskipRnd = status / 10000;
		status = status - nskipRnd * 10000;
		j24 = status / 100;
		i24 = status - 100 * j24;
		if (((j24 < 1) || (j24 > 24)) || ((i24 < 1) || (i24 > 24))) {
			// stop;
			failB = true;

			failS = "*** Error in set_ranlux_state: seeds outside of allowed range!"
					+ " \n"
					+ "   status = "
					+ state[24]
					+ " \n"
					+ "   nskip = "
					+ nskipRnd
					+ " \n"
					+ "   i24 = "
					+ i24
					+ " \n" + "   j24 = " + j24;// +" \n";

			return;
		}
		not_initialized = false;
		return;
	}

	/**
	 * Show summary of Ranlux seeds
	 */
	public static void show_ranlux_seeds() {
		if (carry > 0) {
			icarry = 1;
		} else {
			icarry = 0;
		}
	}

	/**
	 * Show summary of RNG state.
	 */
	public static void SHOW_RNG_STATE() {
		if (ranluxB) {
			show_ranlux_seeds();
		} else// ranmar
		{
		}
	}

	/**
	 * Returns an exponentially distributed, positive, random deviate of unit mean, using 
	 * ran1(idum) as the source of uniform deviates.
	 * @return the result
	 */
	public static double expdev()// long *idum)
	// Returns an exponentially distributed, positive, random deviate of unit
	// mean, using
	// ran1(idum) as the source of uniform deviates.
	{
		// float ran1(long *idum);
		double dum = 0.0;
		do
			dum = ran1();// idum);
		while (dum == 0.0);
		return -Math.log(dum);
	}

	/**
	 * Returns a normally distributed deviate with zero mean and unit variance, using ran1(idum) 
	 * as the source of uniform deviates.
	 * @return the result
	 */
	public static double gasdev()// long *idum)
	// Returns a normally distributed deviate with zero mean and unit variance,
	// using ran1(idum)
	// as the source of uniform deviates.
	{
		// float ran1(long *idum);
		double fac = 0.0;
		double rsq = 0.0;
		double v1 = 0.0;
		double v2 = 0.0;
		if (idum_ran1 < 0)
			iset = 0; // if (*idum < 0) iset=0; Reinitialize.
		if (iset == 0) {// We don’t have an extra deviate handy, so
			do {
				v1 = 2.0 * ran1() - 1.0;// ran1(idum)-1.0; //pick two uniform
										// numbers in the square extending
				// from -1 to +1 in each direction,
				v2 = 2.0 * ran1() - 1.0;// ran1(idum)-1.0;
				rsq = v1 * v1 + v2 * v2;// see if they are in the unit circle,
			} while (rsq >= 1.0 || rsq == 0.0); // and if they are not, try
												// again.
			fac = Math.sqrt(-2.0 * Math.log(rsq) / rsq);
			// Now make the Box-Muller transformation to get two normal
			// deviates. Return one and
			// save the other for next time.
			gset = v1 * fac;
			iset = 1; // Set flag.
			return v2 * fac;
		} else {// We have an extra deviate handy,
			iset = 0; // so unset the flag,
			return gset;// and return it.
		}
	}

	/*
	 * The gamma distribution of integer order a > 0 is the waiting time to the
	 * ath event in a Poisson random process of unit mean. For example, when a =
	 * 1, it is just the exponential distribution of §7.2, the waiting time to
	 * the first event.
	 * 
	 * A gamma deviate has probability pa(x)dx of occurring with a value between
	 * x and x + dx, where pa(x)dx = x^(a-1)exp(-x)/gama(a) dx x >0 (7.3.1) To
	 * generate deviates of (7.3.1) for small values of a, it is best to add up
	 * a exponentially distributed waiting times, i.e., logarithms of uniform
	 * deviates. Since the sum of logarithms is the logarithm of the product,
	 * one really has only to generate the product of a uniform deviates, then
	 * take the log. For larger values of a, the distribution (7.3.1) has a
	 * typically “bell-shaped” form, with a peak at x = a and a half-width of
	 * about ?a.
	 */
	/**
	 * Returns a deviate distributed as a gamma distribution of integer order ia, i.e., a waiting time 
	 * to the iath event in a Poisson process of unit mean, using ran1(idum) as the source of uniform deviates.
	 * @param ia ia
	 * @return the result
	 */
	public static double gamdev(int ia)// , long *idum)
	// Returns a deviate distributed as a gamma distribution of integer order
	// ia, i.e., a waiting time
	// to the iath event in a Poisson process of unit mean, using ran1(idum) as
	// the source of
	// uniform deviates.
	{
		// float ran1(long *idum);
		// void nrerror(char error_text[]);
		failB = false;

		int j = 0;
		double am = 0.0;
		double e = 0.0;
		double s = 0.0;
		double v1 = 0.0;
		double v2 = 0.0;
		double x = 0.0;
		double y = 0.0;
		if (ia < 1) {
			// nrerror("Error in routine gamdev");
			failB = true;
			failS = "Error in routine gamdev";
			return 0.0;
		}
		if (ia < 6) {// Use direct method, adding waiting times.
			x = 1.0;
			for (j = 1; j <= ia; j++)
				x *= ran1();// idum);
			x = -Math.log(x);
		} else {// Use rejection method.
			do {
				do {
					do {
						// These four lines generate the tangent of a random
						// angle, i.e., they
						// are equivalent to y = tan(PI * ran1(idum)).
						v1 = ran1();// idum);
						v2 = 2.0 * ran1() - 1.0;// idum)-1.0;
					} while (v1 * v1 + v2 * v2 > 1.0);
					y = v2 / v1;
					am = ia - 1;
					s = Math.sqrt(2.0 * am + 1.0);
					x = s * y + am; // We decide whether to reject x:
				} while (x <= 0.0);// Reject in region of zero probability.
				e = (1.0 + y * y) * Math.exp(am * Math.log(x / am) - s * y);// Ratio
																			// of
																			// prob.
																			// fn.
																			// to
																			// comparison
																			// fn.
			} while (ran1() > e);// idum) > e); Reject on basis of a second
									// uniform deviate.
		}
		return x;
	}

	/*
	 * The Poisson distribution is conceptually related to the gamma
	 * distribution. It gives the probability of a certain integer number m of
	 * unit rate Poisson random events occurring in a given interval of time x,
	 * while the gamma distribution was the probability of waiting time between
	 * x and x+dx to themth event. Note thatmtakes on only integer values ? 0,
	 * so that the Poisson distribution, viewed as a continuous distribution
	 * function px(m)dm, is zero everywhere except where m is an integer ? 0. At
	 * such places, it is infinite, such that the integrated probability over a
	 * region containing the integer is some finite number. The total
	 * probability at an integer j is Prob(j) = integral j-eps to j+eps of
	 * px(m)dm = x^jexp(-x)/j!
	 */
	/**
	 * Returns as a floating-point number an integer value that is a random deviate drawn from a 
	 * Poisson distribution of mean xm, using ran1(idum) as a source of uniform random deviates.
	 * @param xm xm
	 * @return the result
	 */
	public static double poidev(double xm)// , long *idum)
	// Returns as a floating-point number an integer value that is a random
	// deviate drawn from a
	// Poisson distribution of mean xm, using ran1(idum) as a source of uniform
	// random deviates.
	{
		// float gammln(float xx);
		// float ran1(long *idum);
		// static float sq,alxm,g,oldm=(-1.0);
		double em = 0.0;
		double t = 0.0;
		double y = 0.0;
		if (xm < 12.0) {// Use direct method.
			if (xm != oldm_poidev) {
				oldm_poidev = xm;
				g_poidev = Math.exp(-xm);// If xm is new, compute the
											// exponential.
			}
			em = -1;
			t = 1.0;
			do {// Instead of adding exponential deviates it is equivalent
				// to multiply uniform deviates. We never actually have to take
				// the log,
				// merely compare to the pre-computed exponential.
				++em;
				t *= ran1();// idum);
			} while (t > g_poidev);
		} else {// Use rejection method.
			if (xm != oldm_poidev) {
				// If xm has changed since the last call, then precompute
				// some functions that occur below.
				oldm_poidev = xm;
				sq_poidev = Math.sqrt(2.0 * xm);
				alxm_poidev = Math.log(xm);
				g_poidev = xm * alxm_poidev - SpecFunc.gammln(xm + 1.0);
				// The function gammln is the natural log of the gamma function,
				// as given in §6.1.
			}
			do {
				do {// y is a deviate from a Lorentzian comparison function.
					y = Math.tan(PI * ran1());// idum));
					em = sq_poidev * y + xm; // em is y, shifted and scaled.
				} while (em < 0.0); // Reject if in regime of zero probability.
				em = Math.floor(em);// The trick for integer-valued
									// distributions.
				t = 0.9
						* (1.0 + y * y)
						* Math.exp(em * alxm_poidev - SpecFunc.gammln(em + 1.0)
								- g_poidev);
				// The ratio of the desired distribution to the comparison
				// function; we accept or
				// reject by comparing it to another uniform deviate. The factor
				// 0.9 is chosen so
				// that t never exceeds 1.
			} while (ran1() > t);// idum) > t);
		}
		return em;
	}

	/*
	 * If an event occurs with probability q, and we make n trials, then the
	 * number of times m that it occurs has the binomial distribution, integral
	 * from j-eps to j+eps of pn,q(m)dm = Comb(n,j)q^j(1 - q)^(n-j) (7.3.7) The
	 * binomial distribution is integer valued, with m taking on possible values
	 * from 0 to n. It depends on two parameters, n and q, so is correspondingly
	 * a bit harder to implement than our previous examples.
	 */
	/**
	 * Returns as a floating-point number an integer value that is a random deviate drawn from 
	 * a binomial distribution of n trials each of probability pp, using ran1(idum) as a source of 
	 * uniform random deviates.
	 * @param pp pp
	 * @param n n
	 * @return the result
	 */
	public static double bnldev(double pp, int n)// , long *idum)
	// Returns as a floating-point number an integer value that is a random
	// deviate drawn from
	// a binomial distribution of n trials each of probability pp, using
	// ran1(idum) as a source of
	// uniform random deviates.
	{
		// float gammln(float xx);
		// float ran1(long *idum);
		int j = 0;
		// static int nold=(-1);
		double am = 0.0;
		double em = 0.0;
		double g = 0.0;
		double angle = 0.0;
		double p = 0.0;
		double bnl = 0.0;
		double sq = 0.0;
		double t = 0.0;
		double y = 0.0;
		// static float pold=(-1.0),pc,plog,pclog,en,oldg;
		p = (pp <= 0.5 ? pp : 1.0 - pp);
		// The binomial distribution is invariant under changing pp to 1-pp, if
		// we also change the
		// answer to n minus itself; we’ll remember to do this below.
		am = n * p; // This is the mean of the deviate to be produced.
		if (n < 25) {// Use the direct method while n is not too large.
						// This can require up to 25 calls to ran1.
			bnl = 0.0;
			for (j = 1; j <= n; j++)
				if (ran1() < p)
					++bnl;
		} else if (am < 1.0) {// If fewer than one event is expected out of 25
								// or more trials, then the distribution is
								// quite
								// accurately Poisson. Use direct Poisson
								// method.
			g = Math.exp(-am);
			t = 1.0;
			for (j = 0; j <= n; j++) {
				t *= ran1();
				if (t < g)
					break;
			}
			bnl = (j <= n ? j : n);
		} else {// Use the rejection method.
			if (n != nold_bnldev) {// If n has changed, then compute useful
									// quantities.
				en_bnldev = n;
				oldg_bnldev = SpecFunc.gammln(en_bnldev + 1.0);
				nold_bnldev = n;
			}
			if (p != pold_bnldev) {// If p has changed, then compute useful
									// quantities.
				pc_bnldev = 1.0 - p;
				plog_bnldev = Math.log(p);
				pclog_bnldev = Math.log(pc_bnldev);
				pold_bnldev = p;
			}
			sq = Math.sqrt(2.0 * am * pc_bnldev);// The following code should by
													// now seem familiar:
			// rejection method with a Lorentzian comparison function.
			do {
				do {
					angle = PI * ran1();
					y = Math.tan(angle);
					em = sq * y + am;
				} while (em < 0.0 || em >= (en_bnldev + 1.0));// Reject.
				em = Math.floor(em); // Trick for integer-valued distribution.
				t = 1.2
						* sq
						* (1.0 + y * y)
						* Math.exp(oldg_bnldev - SpecFunc.gammln(em + 1.0)
								- SpecFunc.gammln(en_bnldev - em + 1.0) + em
								* plog_bnldev + (en_bnldev - em) * pclog_bnldev);
			} while (ran1() > t);// Reject. This happens about 1.5 times per
									// deviate,on average.
			bnl = em;
		}
		if (p != pp)
			bnl = n - bnl;// Remember to undo the symmetry transformation.
		return bnl;
	}

	/*
	 * Generation of Random Bits=NOT IMPLEMENTED YET!! Random Sequences Based on
	 * Data Encryption=NOT IMPLEMENTED YET!! MC&Quasi- (that is, Sub-) Random
	 * Sequences=NOT IMPLEMENTED YET!!
	 */

	/*
	 * The VEGAS algorithm, invented by Peter Lepage [1,2], is widely used for
	 * multidimensional integrals that occur in elementary particle physics.
	 * VEGAS is primarily based on importance sampling, but it also does some
	 * stratified sampling if the dimension d is small enough to avoid Kd
	 * explosion (specifically, if (K/2)d < N/2, with N the number of sample
	 * points). The basic technique for importance sampling in VEGAS is to
	 * construct, adaptively, a multidimensional weight function g that is
	 * separable, The weakness of VEGAS is the obvious one: To the extent that
	 * the projection of the function f onto individual coordinate directions is
	 * uniform, VEGAS gives no concentration of sample points in those
	 * dimensions. The worst case for VEGAS, e.g., is an integrand that is
	 * concentrated close to a body diagonal line, e.g., one from (0, 0, 0, . .
	 * .) to (1, 1, 1, . . .). Since this geometry is completely nonseparable,
	 * VEGAS can give no advantage at all. More generally, VEGAS may not do well
	 * when the integrand is concentrated in one-dimensional (or higher) curved
	 * trajectories (or hypersurfaces), unless these happen to be oriented close
	 * to the coordinate directions Note that the user-supplied integrand
	 * function, fxn, has an argument wgt in addition to the expected evaluation
	 * point x. In most applications @@@@you ignore wgt inside the function@@@@.
	 * Occasionally, however, you may want to integrate some additional function
	 * or functions along with the principal function f. The integral of any
	 * such function g can be estimated by.... where the wi’s and x’s are the
	 * arguments wgt and x, respectively. It is straightforward to accumulate
	 * this sum inside your function fxn, and to pass the answer back to your
	 * main program via global variables. Of course, g(x) had better resemble
	 * the principal function f to some degree, since the sampling will be
	 * optimized for f.
	 * 
	 * The input flag init can be used to advantage. One might have a call with
	 * init=0, ncall=1000, itmx=5 immediately followed by a call with init=1,
	 * ncall=100000, itmx=1. The effect would be to develop a sampling grid over
	 * 5 iterations of a small number of samples, then to do a single high
	 * accuracy integration on the optimized grid.
	 */
	// extern long idum; For random number initialization in main.
	// public void vegas(double[] regn, int ndim, float (*fxn)(float [], float),
	// int init,
	// unsigned long ncall, int itmx, int nprn, float *tgral, float *sd,float
	// *chi2a)
	/**
	 * Performs Monte Carlo integration of a user-supplied ndim-dimensional function fxn over a 
	 * rectangular volume specified by regn[1..2*ndim], a vector consisting of ndim “lower left” 
	 * coordinates of the region followed by ndim “upper right” coordinates. The integration consists 
	 * of itmx iterations, each with approximately ncall calls to the function. After each iteration 
	 * the grid is refined; more than 5 or 10 iterations are rarely useful. The input flag init signals 
	 * whether this call is a new start, or a subsequent call for additional iterations. 
	 * The input flag nprn (normally 0) controls the amount of diagnostic output. Returned 
	 * answers are tgral (the best estimate of the integral), sd (its standard deviation), and chi2a ( 
	 * chi2 per degree of freedom, an indicator of whether consistent results are being obtained).
	 * @param regn regn
	 * @param ndim ndim
	 * @param init init
	 * @param ncall ncall
	 * @param itmx itmx
	 * @param nprn nprn
	 */
	public static void vegas(double[] regn, int ndim, int init, int ncall,
			int itmx, int nprn)// , float *tgral, float *sd,float *chi2a)

	// Performs Monte Carlo integration of a user-supplied ndim-dimensional
	// function fxn over a
	// rectangular volume specified by regn[1..2*ndim], a vector consisting of
	// ndim “lower left”
	// coordinates of the region followed by ndim “upper right” coordinates. The
	// integration consists
	// of itmx iterations, each with approximately ncall calls to the function.
	// After each iteration
	// the grid is refined; more than 5 or 10 iterations are rarely useful. The
	// input flag init signals
	// whether this call is a new start, or a subsequent call for additional
	// iterations (see comments
	// below). The input flag nprn (normally 0) controls the amount of
	// diagnostic output. Returned
	// answers are tgral (the best estimate of the integral), sd (its standard
	// deviation), and chi2a
	// (?2 per degree of freedom, an indicator of whether consistent results are
	// being obtained). See
	// text for further details.
	{
		// float ran2(long *idum);
		// void rebin(float rc, int nd, float r[], float xin[], float xi[]);
		// Best make everything static, allowing restarts.
		if (init <= 0) {// Normal entry. Enter here on a cold start.
			mds_vegas = ndo_vegas = 1; // Change to mds=0 to disable stratified
										// sampling,
			// i.e., use importance sampling only.
			for (j_vegas = 1; j_vegas <= ndim; j_vegas++)
				xi_vegas[j_vegas - 1][0] = 1.0;// xi_vegas[j_vegas][1]=1.0;
		}
		if (init <= 1)
			si_vegas = swgt_vegas = schi_vegas = 0.0;
		// Enter here to inherit the grid from a previous call, but not its
		// answers.
		if (init <= 2) {// Enter here to inherit the previous grid and its
						// answers.
			nd_vegas = NDMX;
			ng_vegas = 1;
			if (mds_vegas != 0) // if (mds)
			{// Set up for stratification.
				ng_vegas = (int) Math.pow(ncall / 2.0 + 0.25, 1.0 / ndim);
				mds_vegas = 1;
				if ((2 * ng_vegas - NDMX) >= 0) {
					mds_vegas = -1;
					npg_vegas = ng_vegas / NDMX + 1;
					nd_vegas = ng_vegas / npg_vegas;
					ng_vegas = npg_vegas * nd_vegas;
				}
			}
			for (k_vegas = 1, i_vegas = 1; i_vegas <= ndim; i_vegas++)
				k_vegas *= ng_vegas;
			npg_vegas = Math.max(ncall / k_vegas, 2);
			calls_vegas = (double) npg_vegas * (double) k_vegas;
			dxg_vegas = 1.0 / ng_vegas;
			for (dv2g_vegas = 1, i_vegas = 1; i_vegas <= ndim; i_vegas++)
				dv2g_vegas *= dxg_vegas;
			dv2g_vegas = (calls_vegas * dv2g_vegas)
					* (calls_vegas * dv2g_vegas) / npg_vegas / npg_vegas
					/ (npg_vegas - 1.0);// SQR(calls*dv2g)/npg/npg/(npg-1.0);
			xnd_vegas = nd_vegas;
			dxg_vegas *= xnd_vegas;
			xjac_vegas = 1.0 / calls_vegas;
			for (j_vegas = 1; j_vegas <= ndim; j_vegas++) {
				// dx_vegas[j_vegas]=regn[j_vegas+ndim]-regn[j_vegas];
				dx_vegas[j_vegas - 1] = regn[j_vegas + ndim - 1]
						- regn[j_vegas - 1];
				xjac_vegas *= dx_vegas[j_vegas - 1];// xjac *= dx[j];
			}
			if (nd_vegas != ndo_vegas) {// Do binning if necessary.
				for (i_vegas = 1; i_vegas <= Math.max(nd_vegas, ndo_vegas); i_vegas++)
					r_vegas[i_vegas - 1] = 1.0;// r[i]=1.0;
				for (j_vegas = 1; j_vegas <= ndim; j_vegas++) {
					// rebin(ndo/xnd,nd,r,xin,xi[j]);
					rebin2(ndo_vegas / xnd_vegas, nd_vegas, r_vegas, xin_vegas,
							xi_vegas, j_vegas);
				}
				ndo_vegas = nd_vegas;
			}
			if (nprn >= 0) {/*
							 * printf("%s: ndim= %3d ncall= %8.0f\n",
							 * " Input parameters for vegas",ndim,calls);
							 * printf("%28s it=%5d itmx=%5d\n"," ",it,itmx);
							 * printf
							 * ("%28s nprn=%3d ALPH=%5.2f\n"," ",nprn,ALPH);
							 * printf("%28s mds=%3d nd=%4d\n"," ",mds,nd); for
							 * (j=1;j<=ndim;j++) {
							 * printf("%30s xl[%2d]= %11.4g xu[%2d]= %11.4g\n",
							 * " ",j,regn[j],j,regn[j+ndim]); }
							 */
			}
		}
		for (it_vegas = 1; it_vegas <= itmx; it_vegas++) {
			// Main iteration loop. Can enter here (init ? 3) to do an
			// additional itmx iterations with
			// all other parameters unchanged.
			ti_vegas = tsi_vegas = 0.0;
			for (j_vegas = 1; j_vegas <= ndim; j_vegas++) {
				kg_vegas[j_vegas - 1] = 1;// kg[j]=1;
				for (i_vegas = 1; i_vegas <= nd_vegas; i_vegas++)
					d_vegas[i_vegas - 1][j_vegas - 1] = di_vegas[i_vegas - 1][j_vegas - 1] = 0.0;// d[i][j]=di[i][j]=0.0;
			}
			for (;;) {
				fb_vegas = f2b_vegas = 0.0;
				for (k_vegas = 1; k_vegas <= npg_vegas; k_vegas++) {
					wgt_vegas = xjac_vegas;
					for (j_vegas = 1; j_vegas <= ndim; j_vegas++) {
						// xn_vegas=(kg_vegas[j_vegas]-ran2(&idum))*dxg_vegas+1.0;
						xn_vegas = (kg_vegas[j_vegas - 1] - ran2()) * dxg_vegas
								+ 1.0;
						// ia_vegas[j_vegas]=Math.max(Math.min((int)(xn_vegas),NDMX),1);
						ia_vegas[j_vegas - 1] = Math.max(
								Math.min((int) (xn_vegas), NDMX), 1);
						if (ia_vegas[j_vegas - 1] > 1) // if (ia[j] > 1)
						{
							// xo_vegas=xi_vegas[j_vegas][ia_vegas[j_vegas]]-
							// xi_vegas[j_vegas][ia_vegas[j_vegas]-1];
							xo_vegas = xi_vegas[j_vegas - 1][ia_vegas[j_vegas - 1] - 1]
									- xi_vegas[j_vegas - 1][ia_vegas[j_vegas - 1] - 2];
							// rc=xi[j][ia[j]-1]+(xn-ia[j])*xo;
							rc_vegas = xi_vegas[j_vegas - 1][ia_vegas[j_vegas - 1] - 2]
									+ (xn_vegas - ia_vegas[j_vegas - 1])
									* xo_vegas;
						} else {
							// xo_vegas=xi_vegas[j_vegas][ia_vegas[j_vegas]];
							xo_vegas = xi_vegas[j_vegas - 1][ia_vegas[j_vegas - 1] - 1];
							// rc_vegas=(xn_vegas-ia_vegas[j_vegas])*xo_vegas;
							rc_vegas = (xn_vegas - ia_vegas[j_vegas - 1])
									* xo_vegas;
						}
						// x_vegas[j_vegas]=regn_vegas[j_vegas]+rc_vegas*dx_vegas[j_vegas];
						x_vegas[j_vegas - 1] = regn[j_vegas - 1] + rc_vegas
								* dx_vegas[j_vegas - 1];
						wgt_vegas *= xo_vegas * xnd_vegas;
					}
					f_vegas = wgt_vegas * func.MF(x_vegas);// f=wgt*(*fxn)(x,wgt);
					f2_vegas = f_vegas * f_vegas;
					fb_vegas += f_vegas;
					f2b_vegas += f2_vegas;
					for (j_vegas = 1; j_vegas <= ndim; j_vegas++) {
						// di_vegas[ia_vegas[j_vegas]][j_vegas] += f_vegas;
						di_vegas[ia_vegas[j_vegas - 1] - 1][j_vegas - 1] += f_vegas;
						// if (mds_vegas >= 0)
						// d_vegas[ia_vegas[j_vegas]][j_vegas] += f2_vegas;
						if (mds_vegas >= 0)
							d_vegas[ia_vegas[j_vegas - 1] - 1][j_vegas - 1] += f2_vegas;
					}
				}
				f2b_vegas = Math.sqrt(f2b_vegas * npg_vegas);
				f2b_vegas = (f2b_vegas - fb_vegas) * (f2b_vegas + fb_vegas);
				if (f2b_vegas <= 0.0)
					f2b_vegas = TINY;
				ti_vegas += fb_vegas;
				tsi_vegas += f2b_vegas;
				if (mds_vegas < 0) {// Use stratified sampling.
					for (j_vegas = 1; j_vegas <= ndim; j_vegas++)
						// d[ia[j]][j] += f2b;
						d_vegas[ia_vegas[j_vegas - 1] - 1][j_vegas - 1] += f2b_vegas;
				}
				for (k_vegas = ndim; k_vegas >= 1; k_vegas--) {
					kg_vegas[k_vegas - 1] %= ng_vegas;// kg[k] %= ng;
					// if (++kg[k] != 1) break;
					if (++kg_vegas[k_vegas - 1] != 1)
						break;
				}
				if (k_vegas < 1)
					break;
			}
			tsi_vegas *= dv2g_vegas; // Compute final results for this
										// iteration.
			wgt_vegas = 1.0 / tsi_vegas;
			si_vegas += wgt_vegas * ti_vegas;
			schi_vegas += wgt_vegas * ti_vegas * ti_vegas;
			swgt_vegas += wgt_vegas;
			tgral_vegas = si_vegas / swgt_vegas;
			chi2a_vegas = (schi_vegas - si_vegas * (tgral_vegas))
					/ (it_vegas - 0.9999);
			if (chi2a_vegas < 0.0)
				chi2a_vegas = 0.0;
			sd_vegas = Math.sqrt(1.0 / swgt_vegas);
			tsi_vegas = Math.sqrt(tsi_vegas);
			if (nprn >= 0) {
				// printf("%s %3d : integral = %14.7g +/- %9.2g\n",
				// " iteration no.",it,ti,tsi);
				// printf("%s integral =%14.7g+/-%9.2g chi**2/IT n = %9.2g\n",
				// " all iterations: ",*tgral,*sd,*chi2a);
				// if (nprn!=0) //if (nprn)
				// {
				// for (j=1;j<=ndim;j++)
				// {
				// printf(" DATA FOR axis %2d\n",j);
				// printf("%6s%13s%11s%13s%11s%13s\n",
				// "X","delta i","X","delta i","X","delta i");
				// for (i=1+nprn/2;i<=nd;i += nprn+2)
				// {
				// printf("%8.5f%12.4g%12.5f%12.4g%12.5f%12.4g\n",
				// xi[j][i],di[i][j],xi[j][i+1],
				// di[i+1][j],xi[j][i+2],di[i+2][j]);
				// }
				// }
				// }
			}
			for (j_vegas = 1; j_vegas <= ndim; j_vegas++) {// Refine the grid.
															// Consult
															// references to
															// understand
															// the subtlety of
															// this procedure.
															// The refinement
															// is damped, to
															// avoid rapid,
															// destabilizing
															// changes, and also
															// compressed in
															// range
															// by the exponent
															// ALPH.
				xo_vegas = d_vegas[0][j_vegas - 1];// xo_vegas=d_vegas[1][j_vegas];
				xn_vegas = d_vegas[1][j_vegas - 1];// xn_vegas=d_vegas[2][j_vegas];
				// d_vegas[1][j_vegas]=(xo_vegas+xn_vegas)/2.0;
				d_vegas[0][j_vegas - 1] = (xo_vegas + xn_vegas) / 2.0;
				dt_vegas[j_vegas - 1] = d_vegas[0][j_vegas - 1];// dt_vegas[j_vegas]=d_vegas[1][j_vegas];
				for (i_vegas = 2; i_vegas < nd_vegas; i_vegas++) {
					rc_vegas = xo_vegas + xn_vegas;
					xo_vegas = xn_vegas;
					// xn_vegas=d_vegas[i_vegas+1][j_vegas];
					xn_vegas = d_vegas[i_vegas][j_vegas - 1];
					// d_vegas[i_vegas][j_vegas] = (rc_vegas+xn_vegas)/3.0;
					d_vegas[i_vegas - 1][j_vegas - 1] = (rc_vegas + xn_vegas) / 3.0;
					// dt_vegas[j_vegas] += d_vegas[i_vegas][j_vegas];
					dt_vegas[j_vegas - 1] += d_vegas[i_vegas - 1][j_vegas - 1];
				}
				// d_vegas[nd_vegas][j_vegas]=(xo_vegas+xn_vegas)/2.0;
				d_vegas[nd_vegas - 1][j_vegas - 1] = (xo_vegas + xn_vegas) / 2.0;
				// dt_vegas[j_vegas] += d_vegas[nd_vegas][j_vegas];
				dt_vegas[j_vegas - 1] += d_vegas[nd_vegas - 1][j_vegas - 1];
			}
			for (j_vegas = 1; j_vegas <= ndim; j_vegas++) {
				rc_vegas = 0.0;
				for (i_vegas = 1; i_vegas <= nd_vegas; i_vegas++) {
					// if (d_vegas[i_vegas][j_vegas] < TINY)
					// d_vegas[i_vegas][j_vegas]=TINY;
					if (d_vegas[i_vegas - 1][j_vegas - 1] < TINY)
						d_vegas[i_vegas - 1][j_vegas - 1] = TINY;
					// r_vegas[i_vegas]=pow((1.0-d_vegas[i_vegas][j_vegas]/dt_vegas[j_vegas])/
					// (log(dt_vegas[j_vegas])-log(d_vegas[i_vegas][j_vegas])),ALPH);
					r_vegas[i_vegas - 1] = Math
							.pow((1.0 - d_vegas[i_vegas - 1][j_vegas - 1]
									/ dt_vegas[j_vegas - 1])
									/ (Math.log(dt_vegas[j_vegas - 1]) - Math
											.log(d_vegas[i_vegas - 1][j_vegas - 1])),
									ALPH);

					rc_vegas += r_vegas[i_vegas - 1];// rc_vegas +=
														// r_vegas[i_vegas];
				}
				// rebin(rc/xnd,nd,r,xin,xi[j]);
				rebin2(rc_vegas / xnd_vegas, nd_vegas, r_vegas, xin_vegas,
						xi_vegas, j_vegas);
			}
		}
	}

	/**
	 * Utility routine used by vegas, to rebin a vector of densities xi into new bins defined by a vector r.
	 * @param rc rc
	 * @param nd nd
	 * @param r r
	 * @param xin xin
	 * @param xi xi
	 */
	public static void rebin(double rc, int nd, double[] r, double[] xin,
			double[] xi)
	// Utility routine used by vegas, to rebin a vector of densities xi into new
	// bins defined by a
	// vector r.
	{
		int i = 0;
		int k = 0;
		double dr = 0.0;
		double xn = 0.0;
		double xo = 0.0;
		for (i = 1; i < nd; i++) {
			while (rc > dr)
				dr += r[++k - 1];// r[++k];
			if (k > 1)
				xo = xi[k - 2];// xi[k-1];
			xn = xi[k - 1];// xi[k];
			dr -= rc;
			// xin[i]=xn-(xn-xo)*dr/r[k];
			xin[i - 1] = xn - (xn - xo) * dr / r[k - 1];
		}
		for (i = 1; i < nd; i++)
			xi[i - 1] = xin[i - 1];// xi[i]=xin[i];
		xi[nd - 1] = 1.0;// xi[nd]=1.0;
	}

	/**
	 * Utility routine used by vegas, to rebin a vector of densities xi into new bins defined by a vector r.
	 * @param rc rc
	 * @param nd nd
	 * @param r r
	 * @param xin xin
	 * @param xi xi
	 * @param jj jj
	 */
	public static void rebin2(double rc, int nd, double[] r, double[] xin,
			double[][] xi, int jj)
	// Utility routine used by vegas, to rebin a vector of densities xi into new
	// bins defined by a
	// vector r.
	{
		int i = 0;
		int k = 0;
		double dr = 0.0;
		double xn = 0.0;
		double xo = 0.0;
		for (i = 1; i < nd; i++) {
			while (rc > dr)
				dr += r[++k - 1];// r[++k];
			if (k > 1)
				xo = xi[jj - 1][k - 2];// xo=xi[k-2][jj-1];//xi[k-1];
			xn = xi[jj - 1][k - 1];// xn=xi[k-1][jj-1];//xi[k];
			dr -= rc;
			// xin[i]=xn-(xn-xo)*dr/r[k];
			xin[i - 1] = xn - (xn - xo) * dr / r[k - 1];
		}
		for (i = 1; i < nd; i++)
			xi[jj - 1][i - 1] = xin[i - 1];// xi[i-1][jj-1]=xin[i-1];//xi[i]=xin[i];
		xi[jj - 1][nd - 1] = 1.0;// xi[nd-1][jj-1]=1.0;//xi[nd]=1.0;
	}

	// and anaother
	// Here PFAC is the fraction of remaining function evaluations used at each
	// stage to explore the
	// variance of func. At least MNPT function evaluations are performed in any
	// terminal subregion;
	// a subregion is further bisected only if at least MNBS function
	// evaluations are available. We take
	// MNBS = 4*MNPT.

	// public void miser(float (*func)(float []), float regn[], int ndim,
	// unsigned long npts,
	// float dith, float *ave, float *var)
	/**
	 * Monte Carlo samples a user-supplied ndim-dimensional function func in a rectangular volume 
	 * specified by regn[1..2*ndim], a vector consisting of ndim “lower-left” coordinates of the 
	 * region followed by ndim “upper-right” coordinates. The function is sampled a total of npts 
	 * times, at locations determined by the method of recursive stratified sampling. The mean value 
	 * of the function in the region is returned as ave; an estimate of the statistical uncertainty of ave 
	 * (square of standard deviation) is returned as var. The input parameter dith should normally 
	 * be set to zero, but can be set to (e.g.) 0.1 if func’s active region falls on the boundary of a 
	 * power-of-two subdivision of region.
	 * @param regn regn
	 * @param ndim ndim
	 * @param npts npts
	 * @param dith dith
	 */
	public static void miser(double[] regn, int ndim, int npts, double dith)// ,
																			// float
																			// *ave,
																			// float
																			// *var)
	// Monte Carlo samples a user-supplied ndim-dimensional function func in a
	// rectangular volume
	// specified by regn[1..2*ndim], a vector consisting of ndim “lower-left”
	// coordinates of the
	// region followed by ndim “upper-right” coordinates. The function is
	// sampled a total of npts
	// times, at locations determined by the method of recursive stratified
	// sampling. The mean value
	// of the function in the region is returned as ave; an estimate of the
	// statistical uncertainty of ave
	// (square of standard deviation) is returned as var. The input parameter
	// dith should normally
	// be set to zero, but can be set to (e.g.) 0.1 if func’s active region
	// falls on the boundary of a
	// power-of-two subdivision of region.
	{
		// void ranpt(float pt[], float regn[], int n);
		// float *regn_temp;
		double[] regn_temp;
		int n = 0;
		int npre = 0;
		int nptl = 0;
		int nptr = 0;// long
		int j = 0;
		int jb = 0;
		double avel = 0.0;
		double varl = 0.0;
		double fracl = 0.0;
		double fval = 0.0;
		double rgl = 0.0;
		double rgm = 0.0;
		double rgr = 0.0;
		double s = 0.0;
		double sigl = 0.0;
		double siglb = 0.0;
		double sigr = 0.0;
		double sigrb = 0.0;
		double sum = 0.0;
		double sumb = 0.0;
		double summ = 0.0;
		double summ2 = 0.0;

		double[] fmaxl;
		double[] fmaxr;
		double[] fminl;
		double[] fminr;
		double[] pt;
		double[] rmid;
		pt = new double[ndim];// pt=vector(1,ndim);

		if (npts < MNBS) {// Too few points to bisect; do straight Monte Carlo.
			summ = summ2 = 0.0;
			for (n = 1; n <= npts; n++) {
				ranpt(pt, regn, ndim);
				fval = func.MF(pt);
				summ += fval;
				summ2 += fval * fval;
			}
			ave_miser = summ / npts;
			var_miser = Math.max(TINY, (summ2 - summ * summ / npts)
					/ (npts * npts));
		} else {// Do the preliminary (uniform) sampling.
			rmid = new double[ndim];// vector(1,ndim);
			npre = Math.max((int) (npts * PFAC), MNPT);
			fmaxl = new double[ndim];// vector(1,ndim);
			fmaxr = new double[ndim];// vector(1,ndim);
			fminl = new double[ndim];// vector(1,ndim);
			fminr = new double[ndim];// vector(1,ndim);
			for (j = 1; j <= ndim; j++) { // Initialize the left and right
											// bounds for each dimension.
				iran = (iran * 2661 + 36979) % 175000;
				s = SIGN(dith, (double) (iran - 87500));
				// rmid[j]=(0.5+s)*regn[j]+(0.5-s)*regn[ndim+j];
				rmid[j - 1] = (0.5 + s) * regn[j - 1] + (0.5 - s)
						* regn[ndim + j - 1];
				fminl[j - 1] = fminr[j - 1] = BIG;// fminl[j]=fminr[j]=BIG;
				fmaxl[j - 1] = fmaxr[j - 1] = -BIG;// fmaxl[j]=fmaxr[j] = -BIG;
			}
			for (n = 1; n <= npre; n++) {// Loop over the points in the sample.
				ranpt(pt, regn, ndim);
				fval = func.MF(pt);
				for (j = 1; j <= ndim; j++) {// Find the left and right bounds
												// for each dimension.
					if (pt[j - 1] <= rmid[j - 1]) // if (pt[j]<=rmid[j])
					{
						fminl[j - 1] = Math.min(fminl[j - 1], fval);// fminl[j]=FMIN(fminl[j],fval);
						fmaxl[j - 1] = Math.max(fmaxl[j - 1], fval);// fmaxl[j]=FMAX(fmaxl[j],fval);
					} else {
						fminr[j - 1] = Math.min(fminr[j - 1], fval);// fminr[j]=FMIN(fminr[j],fval);
						fmaxr[j - 1] = Math.max(fmaxr[j - 1], fval);// fmaxr[j]=FMAX(fmaxr[j],fval);
					}
				}
			}
			sumb = BIG; // Choose which dimension jb to bisect.
			jb = 0;
			siglb = sigrb = 1.0;
			for (j = 1; j <= ndim; j++) {
				// if (fmaxl[j] > fminl[j] && fmaxr[j] > fminr[j])
				if (fmaxl[j - 1] > fminl[j - 1] && fmaxr[j - 1] > fminr[j - 1]) {
					// sigl=FMAX(TINY,pow(fmaxl[j]-fminl[j],2.0/3.0));
					sigl = Math.max(TINY,
							Math.pow(fmaxl[j - 1] - fminl[j - 1], 2.0 / 3.0));
					// sigr=FMAX(TINY,pow(fmaxr[j]-fminr[j],2.0/3.0));
					sigr = Math.max(TINY,
							Math.pow(fmaxr[j - 1] - fminr[j - 1], 2.0 / 3.0));
					sum = sigl + sigr;// Equation (7.8.24), see text.
					if (sum <= sumb) {
						sumb = sum;
						jb = j;
						siglb = sigl;
						sigrb = sigr;
					}
				}
			}
			// free_vector(fminr,1,ndim);
			// free_vector(fminl,1,ndim);
			// free_vector(fmaxr,1,ndim);
			// free_vector(fmaxl,1,ndim);
			// if (!jb) jb=1+(ndim*iran)/175000; //MNPT may be too small.
			if (jb == 0)
				jb = 1 + (ndim * iran) / 175000;
			// rgl=regn[jb]; //Apportion the remaining points between left and
			// right.
			rgl = regn[jb - 1];
			rgm = rmid[jb - 1];// rgm=rmid[jb];
			rgr = regn[ndim + jb - 1];// rgr=regn[ndim+jb];
			fracl = Math.abs((rgm - rgl) / (rgr - rgl));
			nptl = (int) (MNPT + (npts - npre - 2 * MNPT) * fracl * siglb
					/ (fracl * siglb + (1.0 - fracl) * sigrb));// Equation
																// (7.8.23).
			nptr = npts - npre - nptl;
			regn_temp = new double[2 * ndim];// vector(1,2*ndim);
			// Now allocate and integrate the two subregions.
			for (j = 1; j <= ndim; j++) {
				regn_temp[j - 1] = regn[j - 1];// regn_temp[j]=regn[j];
				regn_temp[ndim + j - 1] = regn[ndim + j - 1];// regn_temp[ndim+j]=regn[ndim+j];
			}
			regn_temp[ndim + jb - 1] = rmid[jb - 1];// regn_temp[ndim+jb]=rmid[jb];
			// miser(func,regn_temp,ndim,nptl,dith,&avel,&varl);
			miser(regn_temp, ndim, nptl, dith);
			// regn_temp[jb]=rmid[jb]; Dispatch recursive call; will return back
			// here eventually.
			regn_temp[jb - 1] = rmid[jb - 1];
			// regn_temp[ndim+jb]=regn[ndim+jb];
			regn_temp[ndim + jb - 1] = regn[ndim + jb - 1];
			// miser(func,regn_temp,ndim,nptr,dith,ave,var);
			miser(regn_temp, ndim, nptr, dith);
			// free_vector(regn_temp,1,2*ndim);
			ave_miser = fracl * avel + (1 - fracl) * (ave_miser);
			var_miser = fracl * fracl * varl + (1 - fracl) * (1 - fracl)
					* (var_miser);
			// Combine left and right regions by equation (7.8.11) (1st line).
			// free_vector(rmid,1,ndim);
		}
		// free_vector(pt,1,ndim);
	}

	// The miser routine calls a short function ranpt to get a random point
	// within a specified
	// d-dimensional region. The following version of ranpt makes consecutive
	// calls to a uniform
	// random number generator and does the obvious scaling. One can easily
	// modify ranpt to
	// generate its points via the quasi-random routine sobseq (§7.7). We find
	// that miser with
	// sobseq can be considerably more accurate than miser with uniform random
	// deviates. Since
	// the use of RSS and the use of quasi-random numbers are completely
	// separable, however, we
	// have not made the code given here dependent on sobseq. A similar remark
	// might be made
	// regarding importance sampling, which could in principle be combined with
	// RSS. (One could
	// in principle combine vegas and miser, although the programming would be
	// intricate.)

	// extern long idum;
	/**
	 * Returns a uniformly random point pt in an n-dimensional rectangular region. Used by miser; 
	 * calls ran1 for uniform deviates. Your main program should initialize the global variable idum 
	 * to a negative seed integer.
	 * @param pt pt
	 * @param regn regn
	 * @param n n
	 */
	public static void ranpt(double pt[], double regn[], int n)
	// Returns a uniformly random point pt in an n-dimensional rectangular
	// region. Used by miser;
	// calls ran1 for uniform deviates. Your main program should initialize the
	// global variable idum
	// to a negative seed integer.
	{
		// float ran1(long *idum);
		int j = 0;
		for (j = 1; j <= n; j++)
			// pt[j]=regn[j]+(regn[n+j]-regn[j])*ran1(&idum);
			pt[j - 1] = regn[j - 1] + (regn[n + j - 1] - regn[j - 1]) * ran1();
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
}
