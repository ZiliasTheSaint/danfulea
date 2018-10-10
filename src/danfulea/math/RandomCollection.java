package danfulea.math;

/**
 * Class used for handling several random number methods.
 * 
 * @author Dan Fulea, 14 APR. 2011
 * 
 */
public class RandomCollection {

	public static int RandomUse = 1;// if 0- old Java random
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

	/**
	 * Reseting variables for random number generator.
	 */
	public static void reset() {
		rng_array1 = new double[$NRANMAR];
		rng_seed = 999999;
		jseed_dflt = 314159265;
		icon = 2147483563;
		rng_array = new double[24];
		seeds = new int[24];
		seedin = 0;
		luxury_level = 0;
		state = new int[25];
		carry = 0;
		i24 = 0;
		j24 = 0;
		next = new int[24];
		nskipRnd = 0;
		status = 0;
		jseed = 0;
		icarry = 0;
		kRnd = 0;
		jRnd = 0;
		not_initialized = true;
		uni = 0;
		twom24 = 0.0;
		twop24 = 0.0;
		urndm = new int[97];
		crndm = 0;
		cdrndm = 0;
		cmrndm = 0;
		i4opt = 0;
		ixx = 0;
		jxx = 0;
		fool_optimizer = 0;
	}

	/**
	 * Several methods for random number generator. Default = the EGS based
	 * method
	 * 
	 * @return a random number
	 */
	public static double random01() {
		// To get <0,1> instead of <0,1), we substract the random number from 1
		// in half cases
		double result = 0.0;
		if (RandomUse == 0) {
			double randomNumber = Math.random();
			if (Math.random() < 0.5)
				randomNumber = 1.0 - randomNumber;
			return randomNumber;
		} else if (RandomUse == 1) {
			result = RANDOMSET();
		} else if (RandomUse == 2) {
			result = Math.random();
		}

		return result;
	}

	/**
	 * EGS-based random number generator.
	 * 
	 * @return random number
	 */
	private static double RANDOMSET() {
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
	 * method used by the EGS-based random number generator.
	 * 
	 */
	private static void ranmar_get() {
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
	 * method used by the EGS-based random number generator.
	 * 
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
	 * Method used by the EGS-based random number generator.
	 * @param rng_array the rng array
	 */
	private static void ranlux(double[] rng_array) {

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
	 * Method used by the EGS-based random number generator.
	 * @param luxury_level1 luxury_level
	 * @param seedin1 the seed
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
	 * Method used by the EGS-based random number generator.
	 * @param state the state array
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
	 * Method used by the EGS-based random number generator.
	 * @param state the state array
	 */
	public void set_ranlux_state(int[] state) {
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
			/*
			 * STOPPROGRAM=true;
			 * 
			 * seqStr=
			 * "*** Error in set_ranlux_state: seeds outside of allowed range!"
			 * +" \n"+ "   status = "+state[24]+" \n"+
			 * "   nskip = "+nskipRnd+" \n"+ "   i24 = "+i24+" \n"+
			 * "   j24 = "+j24;//+" \n"; //if(iprint>1)
			 * eq.printSequence(seqStr);
			 */

			return;
		}
		not_initialized = false;
		return;
	}

	/**
	 * method used by the EGS-based random number generator.
	 * 
	 */
	public static void show_ranlux_seeds() {
		if (carry > 0) {
			icarry = 1;
		} else {
			icarry = 0;
		}
		/*
		 * seqStr="' skip = "+format(nskipRnd,4)+",  ix jx = "+format(i24,3)+" ,"
		 * +format(j24,3)+ " carry = "+format(icarry,2); if(iprint>1)
		 * eq.printSequence(seqStr);
		 */

	}

	/**
	 * method used by the EGS-based random number generator.
	 * 
	 */
	public static void SHOW_RNG_STATE() {
		if (ranluxB) {
			show_ranlux_seeds();
		} else// ranmar
		{
			/*
			 * seqStr="  ixx jxx = "+format(ixx,4)+" ,"+format(jxx,4);
			 * if(iprint>1) eq.printSequence(seqStr);
			 */
		}
	}

	// =============
	public static double negExp(double lambda) {
		double y, x;

		lambda = Math.abs(lambda);
		if (lambda != 0.0) {
			do {
				y = Math.random();
			} while (y >= 1.0); /*
								 * Can't accept 1.0 because of ln(1-y) --> would
								 * go to -inf
								 */

			x = (-1.0 * Math.log(1.0 - y)) / lambda;
			return x;
		} else
			return (0.0);
	} // negExp

	/**
	 * Returns a random number having uniform distribution with parameters a and
	 * b.
	 * 
	 * @param a
	 *            the lower limit of the interval &lt;a,b&gt;
	 * @param b
	 *            the upper limit of the interval &lt;a,b&gt;
	 * 
	 * @return a random number having uniform distribution with parameters a and
	 *         b
	 */
	public static double uniform(double a, double b) {
		double xchg;
		double randomNumber;

		if (b < a) {
			xchg = b;
			b = a;
			a = xchg;
		}

		/*
		 * To get <a,b> instead of <a,b), we substract the random number from 1
		 * in half cases
		 */
		randomNumber = Math.random();
		if (Math.random() < 0.5)
			randomNumber = 1.0 - randomNumber;

		return ((randomNumber * (b - a)) + a);
	} // uniform

	/**
	 * Returns true with given probability.
	 * 
	 * @param probability
	 *            probability that true will be returned
	 * 
	 * @return true if a randomly generated number from &lt;0,1) is less than or
	 *         equal to the parameter
	 */
	public static boolean draw(double probability) {
		if (Math.random() > probability)
			return false;
		else
			return true;
	} // draw

	/**
	 * Returns a random number having Gaussian (normal) distribution with
	 * parameters mi and sigma
	 * 
	 * @param mu
	 *            parameter mu (mean value) of the normal-distribution generator
	 * @param sigma
	 *            parameter sigma (variance) of the normal-distribution
	 *            generator
	 * 
	 * @return a random number having normal distribution with parameters mu and
	 *         sigma
	 */
	public static double gauss(double mu, double sigma) {
		int i;
		double sum;

		sum = 0.0;
		for (i = 0; i < 12; i++)
			sum += uniform(0.0, 1.0);

		return (sigma * (sum - 6.0) + mu);
	} // gauss

}
