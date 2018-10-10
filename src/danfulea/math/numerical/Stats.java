package danfulea.math.numerical;

/**
 * Statistical description of data<br>
 * Based on literature, including Numerical Recipes in C (Cambridge Univ.).
 * @author Dan Fulea, 16 OCT. 2006
 */
public class Stats {
	public static Function func;
	public static boolean failB = false;
	public static String failS = "";

	public static double ave_moment = 0.0;
	public static double adev_moment = 0.0;
	public static double sdev_moment = 0.0;
	public static double var_moment = 0.0;
	public static double skew_moment = 0.0;
	public static double curt_moment = 0.0;

	public static double t_ttest = 0.0;
	public static double prob_ttest = 0.0;

	public static double ave_avevar = 0.0;
	public static double var_avevar = 0.0;

	public static double t_tutest = 0.0;
	public static double prob_tutest = 0.0;

	public static double t_tptest = 0.0;
	public static double prob_tptest = 0.0;

	public static double f_ftest = 0.0;
	public static double prob_ftest = 0.0;
	// if probab>0.95 (95%) then the means are different
	public static double prob_reff = 0.95;
	public static double prob_2tail_fordiff_ttest = 0.0;
	public static double prob_1tail_fordiff_ttest = 0.0;
	public static double prob_2tail_fordiff_tutest = 0.0;
	public static double prob_1tail_fordiff_tutest = 0.0;
	public static double prob_2tail_fordiff_tptest = 0.0;
	public static double prob_1tail_fordiff_tptest = 0.0;
	public static String significance_ttest = "";
	public static String significance_tutest = "";
	public static String significance_tptest = "";
	public static double prob_2tail_fordiff_ftest = 0.0;// "2 tail"
	public static double prob_1tail_fordiff_ftest = 0.0;// "1 tail"
	public static String significance_ftest = "";
	public static boolean differentB = false;// true if different statistically,
												// false otherwise!

	public static double df_chsone = 0.0;
	public static double chsq_chsone = 0.0;
	public static double prob_chsone = 0.0;

	public static double df_chstwo = 0.0;
	public static double chsq_chstwo = 0.0;
	public static double prob_chstwo = 0.0;

	public static double d_ksone = 0.0;
	public static double prob_ksone = 0.0;

	public static double d_kstwo = 0.0;
	public static double prob_kstwo = 0.0;

	public static double EPS1 = 0.001;
	public static double EPS2 = 1.0e-8;
	public static double TINY = 1.0e-30;

	public static double chisq_cntab1 = 0.0;
	public static double df_cntab1 = 0.0;
	public static double prob_cntab1 = 0.0;
	public static double cramrv_cntab1 = 0.0;
	public static double ccc_cntab1 = 0.0;

	public static double h_cntab2 = 0.0;
	public static double hx_cntab2 = 0.0;
	public static double hy_cntab2 = 0.0;
	public static double hygx_cntab2 = 0.0;
	public static double hxgy_cntab2 = 0.0;
	public static double uygx_cntab2 = 0.0;
	public static double uxgy_cntab2 = 0.0;
	public static double uxy_cntab2 = 0.0;

	public static double TINY2 = 1.0e-20;
	public static double r_pearsn = 0.0;
	public static double prob_pearsn = 0.0;
	public static double z_pearsn = 0.0;

	public static double d_spear = 0.0;
	public static double zd_spear = 0.0;
	public static double probd_spear = 0.0;
	public static double rs_spear = 0.0;
	public static double probrs_spear = 0.0;
	public static double s_crank = 0.0;

	public static double tau_kendl1 = 0.0;
	public static double z_kendl1 = 0.0;
	public static double prob_kendl1 = 0.0;

	public static double tau_kendl2 = 0.0;
	public static double z_kendl2 = 0.0;
	public static double prob_kendl2 = 0.0;

	public static double fa_quadvl = 0.0;
	public static double fb_quadvl = 0.0;
	public static double fc_quadvl = 0.0;
	public static double fd_quadvl = 0.0;
	// float *fa, float *fb, float *fc, float *fd)
	public static double fa_quadct = 0.0;
	public static double fb_quadct = 0.0;
	public static double fc_quadct = 0.0;
	public static double fd_quadct = 0.0;

	public static double d1_ks2d1s = 0.0;
	public static double prob_ks2d1s = 0.0;
	// float *d, float *prob
	public static double d_ks2d2s = 0.0;
	public static double prob_ks2d2s = 0.0;

	/*
	 * Moments of a Distribution: Mean, Variance, Skewness, and So Forth
	 * 
	 * Here we will be content to note that the N - 1 should be changed to N if
	 * you are ever in the situation of measuring the variance of a distribution
	 * whose mean x is known a priori rather than being estimated from the data.
	 * (We might also comment that if the difference between N and N - 1 ever
	 * matters to you, then you are probably up to no good anyway — e.g., trying
	 * to substantiate a questionable hypothesis with marginal data.)
	 * 
	 * NOTE:var definition corrected: The second sum would be zero if x(mean)
	 * were exact, but otherwise it does a good job of correcting the roundoff
	 * error in the first term.
	 */
	/**
	 * Given an array of data[1..n], this routine returns its mean ave, average deviation adev, 
	 * standard deviation sdev, variance var, skewness skew, and kurtosis curt.
	 * @param data data
	 * @param n n
	 */
	public static void moment(double[] data, int n)// , float *ave, float *adev,
													// float *sdev,
	// float *var, float *skew, float *curt)
	// Given an array of data[1..n], this routine returns its mean ave, average
	// deviation adev,
	// standard deviation sdev, variance var, skewness skew, and kurtosis curt.
	{
		// void nrerror(char error_text[]);
		failB = false;
		int j = 0;
		double ep = 0.0;
		double s = 0.0;
		double p = 0.0;
		if (n <= 1) {
			failB = true;
			failS = "n must be at least 2 in moment";
			return;
			// nrerror("n must be at least 2 in moment");
		}
		s = 0.0;// First pass to get the mean.
		for (j = 1; j <= n; j++)
			s += data[j - 1];// s += data[j];
		ave_moment = s / n;
		adev_moment = (var_moment) = (skew_moment) = (curt_moment) = 0.0;
		// Second pass to get the first (absolute), second, third, and fourth
		// moments of the
		// deviation from the mean.
		for (j = 1; j <= n; j++) {
			adev_moment += Math.abs(s = data[j - 1] - (ave_moment));// Math.abs(s=data[j]-(*ave));
			ep += s;
			var_moment += (p = s * s);
			skew_moment += (p *= s);
			curt_moment += (p *= s);
		}
		adev_moment /= n;
		var_moment = (var_moment - ep * ep / n) / (n - 1); // Corrected two-pass
															// formula.
		sdev_moment = Math.sqrt(var_moment);// Put the pieces together according
											// to the conventional
		// definitions.
		if (var_moment != 0.0)// if (var_moment)
		{
			skew_moment /= (n * (var_moment) * (sdev_moment));
			curt_moment = (curt_moment) / (n * (var_moment) * (var_moment))
					- 3.0;
		} else {
			failB = true;
			failS = "No skew/kurtosis when variance = 0 (in moment)";
			return;
			// nrerror("No skew/kurtosis when variance = 0 (in moment)");
		}
	}

	/*
	 * Do Two Distributions Have the Same Means or Variances?
	 * 
	 * Applying the concept of standard error, the conventional statistic for
	 * measuring the significance of a difference of means is termed Student’s
	 * t. When the two distributions are thought to have the same variance, but
	 * possibly different means, then Student’s t is computed as follows: First,
	 * estimate the standard error of the difference of the means, sD, from the
	 * “pooled variance” by the formula sD = SQRT[[(sum i of (xi-xAmean)^2+sum i
	 * of (xi-xBmean)^2)]*(1/NA+1/NB)/(NA+NB+2)] (14.2.1) where each sum is over
	 * the points in one sample, the first or second, each mean likewise refers
	 * to one sample or the other, and NA and NB are the numbers of points in
	 * the first and second samples, respectively. Second, compute t by t =
	 * [xA(mean) - xB(mean)]/sD (14.2.2) Third, evaluate the significance of
	 * this value of t for Student’s distribution with NA + NB - 2 degrees of
	 * freedom, by equations (6.4.7) and (6.4.9), and by the routine betai
	 * (incomplete beta function) of §6.4. The significance is a number between
	 * zero and one, and is the probability that |t| could be this large or
	 * larger just by chance, for distributions with equal means. Therefore, a
	 * small numerical value of the significance (0.05 or 0.01) means that the
	 * observed difference is “very significant.” The function A(t|?) in
	 * equation (6.4.7) is one minus the significance.
	 */
	/**
	 * Given the arrays data1[1..n1] and data2[1..n2], this routine returns Student’s t as t, 
	 * and its significance as prob, small values of prob indicating that the arrays have significantly 
	 * different means. The data arrays are assumed to be drawn from populations with the same true variance (PROBAB=NUMBER NOT %).
	 * @param data1 data1
	 * @param n1 n1
	 * @param data2 data2
	 * @param n2 n2
	 */
	public static void ttest(double[] data1, int n1, double[] data2, int n2)// ,
	// float *t, float *prob)
	// Given the arrays data1[1..n1] and data2[1..n2], this routine returns
	// Student’s t as t,
	// and its significance as prob, small values of prob indicating that the
	// arrays have significantly
	// different means. The data arrays are assumed to be drawn from populations
	// with the same
	// true variance.
	// PROBAB=NUMBER NOT % !!!
	{
		// void avevar(float data[], unsigned long n, float *ave, float *var);
		// float betai(float a, float b, float x);
		double var1 = 0.0;
		double var2 = 0.0;
		double svar = 0.0;
		double df = 0.0;
		double ave1 = 0.0;
		double ave2 = 0.0;

		avevar(data1, n1);// ,&ave1,&var1);
		ave1 = ave_avevar;
		var1 = var_avevar;
		avevar(data2, n2);// ,&ave2,&var2);
		ave2 = ave_avevar;
		var2 = var_avevar;

		/*
		 * MAB: ttestmediirate function...THE SAME!!! t:=abs(r1-r2);
		 * t:=t*power((f1+1)*(f2+1)/(f1+f2+2),0.5);
		 * t:=t/(power((f1*s1*s1+f2*s2*s2)/(f1+f2),0.5));
		 */
		df = n1 + n2 - 2; // Degrees of freedom.
		svar = ((n1 - 1) * var1 + (n2 - 1) * var2) / df; // Pooled variance.

		t_ttest = (ave1 - ave2) / Math.sqrt(svar * (1.0 / n1 + 1.0 / n2));
		prob_ttest = SpecFunc.betai(0.5 * df, 0.5, df
				/ (df + (t_ttest) * (t_ttest))); // See equation (6.4.9).

		prob_2tail_fordiff_ttest = 1.0 - prob_ttest;// 1-alpha
		prob_1tail_fordiff_ttest = 1.0 - prob_ttest / 2.0;
		if (prob_1tail_fordiff_ttest > prob_reff) {
			differentB = true;
			significance_ttest = "The means are different!";
		} else {
			differentB = false;
			significance_ttest = "The means are NOT different!";
		}
	}

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

	/*
	 * The next case to consider is where the two distributions have
	 * significantly different variances, but we nevertheless want to know if
	 * their means are the same or different. (A treatment for baldness has
	 * caused some patients to lose all their hair and turned others into
	 * werewolves, but we want to know if it helps cure baldness on the
	 * average!) Be suspicious of the unequal-variance t-test: If two
	 * distributions have very different variances, then they may also be
	 * substantially different in shape; in that case, the difference of the
	 * means may not be a particularly useful thing to know. To find out whether
	 * the two data sets have variances that are significantly different, you
	 * use the F-test, described later on in this section. The relevant
	 * statistic for the unequal variance t-test is t = [xA - xB]/[Var(xA)/NA +
	 * Var(xB)/NB]^1/2 (14.2.3) This statistic is distributed approximately as
	 * Student’s t with a number of degrees of freedom equal to [Var(xA)/NA +
	 * Var(xB)/NB ]^2/{[Var(xA)/NA]2/(NA - 1) + [Var(xB)/NB]2/(NB - 1)}
	 */

	/**
	 * Given the arrays data1[1..n1] and data2[1..n2], this routine returns Student’s t as t, and 
	 * its significance as prob, small values of prob indicating that the arrays have significantly different 
	 * means. The data arrays are allowed to be drawn from populations with unequal variances.
	 * @param data1 data1
	 * @param n1 n1
	 * @param data2 data2
	 * @param n2 n2
	 */
	public static void tutest(double[] data1, int n1, double[] data2, int n2)// ,
	// float *t, float *prob)
	// Given the arrays data1[1..n1] and data2[1..n2], this routine returns
	// Student’s t as t, and
	// its significance as prob, small values of prob indicating that the arrays
	// have significantly different
	// means. The data arrays are allowed to be drawn from populations with
	// unequal variances.
	{
		// void avevar(float data[], unsigned long n, float *ave, float *var);
		// float betai(float a, float b, float x);
		double var1 = 0.0;
		double var2 = 0.0;
		double df = 0.0;
		double ave1 = 0.0;
		double ave2 = 0.0;

		avevar(data1, n1);// ,&ave1,&var1);
		ave1 = ave_avevar;
		var1 = var_avevar;
		avevar(data2, n2);// ,&ave2,&var2);
		ave2 = ave_avevar;
		var2 = var_avevar;

		t_tutest = (ave1 - ave2) / Math.sqrt(var1 / n1 + var2 / n2);
		df = (var1 / n1 + var2 / n2)
				* (var1 / n1 + var2 / n2)
				/ ((var1 / n1) * (var1 / n1) / (n1 - 1) + (var2 / n2)
						* (var2 / n2) / (n2 - 1));
		prob_tutest = SpecFunc.betai(0.5 * df, 0.5, df
				/ (df + (t_tutest) * (t_tutest)));

		prob_2tail_fordiff_tutest = 1.0 - prob_tutest;// 1-alpha
		prob_1tail_fordiff_tutest = 1.0 - prob_tutest / 2.0;
		if (prob_1tail_fordiff_tutest > prob_reff) {
			differentB = true;
			significance_tutest = "The means are different!";
		} else {
			differentB = false;
			significance_tutest = "The means are NOT different!";
		}
	}

	/*
	 * Our final example of a Student’s t test is the case of paired samples.
	 * Here we imagine that much of the variance in both samples is due to
	 * effects that are point-by-point identical in the two samples. For
	 * example, we might have two job candidates who have each been rated by the
	 * same ten members of a hiring committee. We want to know if the means of
	 * the ten scores differ significantly. We first try ttest above, and obtain
	 * a value of prob that is not especially significant (e.g., > 0.05). But
	 * perhaps the significance is being washed out by the tendency of some
	 * committee members always to give high scores, others always to give low
	 * scores, which increases the apparent variance and thus decreases the
	 * significance of any difference in the means. We thus try the
	 * paired-sample formulas,
	 * 
	 * Cov(xA, xB) =1/(N - 1) sum i=1 to N of (xAi - xAMed)(xBi - xBMed)
	 * (14.2.5) sD = {[Var(xA) + Var(xB) - 2Cov(xA, xB)]/N }^1/2 (14.2.6) t =
	 * (xA - xB)/sD
	 */

	/**
	 * Given the paired arrays data1[1..n] and data2[1..n], this routine returns Student’s t for 
	 * paired data as t, and its significance as prob, small values of prob indicating a significant difference of means.
	 * @param data1 data1
	 * @param data2 data2
	 * @param n n
	 */
	public static void tptest(double[] data1, double[] data2, int n)// , float
																	// *t,
	// float *prob)
	// Given the paired arrays data1[1..n] and data2[1..n], this routine returns
	// Student’s t for
	// paired data as t, and its significance as prob, small values of prob
	// indicating a significant
	// difference of means.
	{
		// void avevar(float data[], unsigned long n, float *ave, float *var);
		// float betai(float a, float b, float x);
		int j = 0;// unsigned long j;
		double var1 = 0.0;
		double var2 = 0.0;
		double ave1 = 0.0;
		double ave2 = 0.0;
		double sd = 0.0;
		double df = 0.0;
		double cov = 0.0;
		avevar(data1, n);// ,&ave1,&var1);
		ave1 = ave_avevar;
		var1 = var_avevar;
		avevar(data2, n);// ,&ave2,&var2);
		ave2 = ave_avevar;
		var2 = var_avevar;

		for (j = 1; j <= n; j++)
			cov += (data1[j - 1] - ave1) * (data2[j - 1] - ave2);// (data1[j]-ave1)*(data2[j]-ave2);
		cov /= df = n - 1;
		sd = Math.sqrt((var1 + var2 - 2.0 * cov) / n);
		t_tptest = (ave1 - ave2) / sd;
		prob_tptest = SpecFunc.betai(0.5 * df, 0.5, df
				/ (df + (t_tptest) * (t_tptest)));

		prob_2tail_fordiff_tptest = 1.0 - prob_tptest;// 1-alpha
		prob_1tail_fordiff_tptest = 1.0 - prob_tptest / 2.0;
		if (prob_1tail_fordiff_tptest > prob_reff) {
			differentB = true;
			significance_tptest = "The means are different!";
		} else {
			differentB = false;
			significance_tptest = "The means are NOT different!";
		}

	}

	/*
	 * F-Test for Significantly Different Variances
	 * 
	 * The F-test tests the hypothesis that two samples have different variances
	 * by trying to reject the null hypothesis that their variances are actually
	 * consistent. The statistic F is the ratio of one variance to the other, so
	 * values either  1 or  1 will indicate very significant differences. The
	 * distribution of F in the null case is given in equation (6.4.11), which
	 * is evaluated using the routine betai. In the most common case, we are
	 * willing to disprove the null hypothesis (of equal variances) by either
	 * very large or very small values of F, so the correct significance is
	 * two-tailed, the sum of two incomplete beta functions. It turns out, by
	 * equation (6.4.3), that the two tails are always equal; we need compute
	 * only one, and double it. Occasionally, when the null hypothesis is
	 * strongly viable, the identity of the two tails can become confused,
	 * giving an indicated probability greater than one. Changing the
	 * probability to two minus itself correctly exchanges the tails.
	 */
	/**
	 * Given the arrays data1[1..n1] and data2[1..n2], this routine returns the value of f, and 
	 * its significance as prob. Small values of prob indicate that the two arrays have significantly 
	 * different variances.
	 * @param data1 data1
	 * @param n1 n1
	 * @param data2 data2
	 * @param n2 n2
	 */
	public static void ftest(double[] data1, int n1, double[] data2, int n2)// ,
	// float *f, float *prob)
	// Given the arrays data1[1..n1] and data2[1..n2], this routine returns the
	// value of f, and
	// its significance as prob. Small values of prob indicate that the two
	// arrays have significantly
	// different variances.
	{
		// void avevar(float data[], unsigned long n, float *ave, float *var);
		// float betai(float a, float b, float x);
		double var1 = 0.0;
		double var2 = 0.0;
		// double ave1 = 0.0;
		// double ave2 = 0.0;
		double df1 = 0.0;
		double df2 = 0.0;

		avevar(data1, n1);// ,&ave1,&var1);
		// ave1 = ave_avevar;
		var1 = var_avevar;
		avevar(data2, n2);// ,&ave2,&var2);
		// ave2 = ave_avevar;
		var2 = var_avevar;

		if (var1 > var2) {// Make F the ratio of the larger variance to the
							// smaller one.
			f_ftest = var1 / var2;
			df1 = n1 - 1;
			df2 = n2 - 1;
		} else {
			f_ftest = var2 / var1;
			df1 = n2 - 1;
			df2 = n1 - 1;
		}
		// =======test========df1=df2=12;f_ftest=2.69;
		prob_ftest = 2.0 * SpecFunc.betai(0.5 * df2, 0.5 * df1, df2
				/ (df2 + df1 * (f_ftest)));
		if (prob_ftest > 1.0)
			prob_ftest = 2.0 - prob_ftest;

		prob_2tail_fordiff_ftest = 1.0 - prob_ftest;
		prob_1tail_fordiff_ftest = 1.0 - prob_ftest / 2.0;
		if (prob_1tail_fordiff_ftest > prob_reff) {
			differentB = true;
			significance_ftest = "The variances are different!";
		} else {
			differentB = false;
			significance_ftest = "The variances are NOT different!";
		}

	}

	/*
	 * Are Two Distributions Different?
	 * 
	 * Given two sets of data, we can generalize the questions asked in the
	 * previous section and ask the single question: Are the two sets drawn from
	 * the same distribution function, or from different distribution functions?
	 * Equivalently, in proper statistical language, “Can we disprove, to a
	 * certain required level of significance, the null hypothesis that two data
	 * sets are drawn from the same population distribution function?”
	 * Disproving the null hypothesis in effect proves that the data sets are
	 * from different distributions. Failing to disprove the null hypothesis, on
	 * the other hand, only shows that the data sets can be consistent with a
	 * single distribution function. One can never prove that two data sets come
	 * from a single distribution, since (e.g.) no practical amount of data can
	 * distinguish between two distributions which differ only by one part in
	 * 1010. Proving that two distributions are different, or showing that they
	 * are consistent, is a task that comes up all the time in many areas of
	 * research: Are the visible stars distributed uniformly in the sky? (That
	 * is, is the distribution of stars as a function of declination — position
	 * in the sky — the same as the distribution of sky area as a function of
	 * declination?) Are educational patterns the same in Brooklyn as in the
	 * Bronx? (That is, are the distributions of people as a function of
	 * last-grade-attended the same?) Do two brands of fluorescent lights have
	 * the same distribution of burn-out times? Is the incidence of chicken pox
	 * the same for first-born, second-born, third-born children, etc.?
	 * 
	 * Chi-Square Test Suppose that Ni is the number of events observed in the
	 * ith bin, and that ni is the number expected according to some known
	 * distribution. Note that the Ni’s are integers, while the ni’s may not be.
	 * .... where the sum is over all bins. A large value of chi2 indicates that
	 * the null hypothesis (that theNi’s are drawnfrom the population
	 * represented by the ni’s) is rather unlikely. The chi-square probability
	 * function Q(?2|?) is an incomplete gamma function, and was already
	 * discussed in §6.2 (see equation 6.2.18). Strictly speaking Q(?2|?) is the
	 * probability that the sum of the squares of ? random normal variables of
	 * unit variance (and zero mean) will be greater than ?2. The terms in the
	 * sum (14.3.1) are not individually normal. However, if either the number
	 * of bins is large ( 1), or the number of events in each bin is large (
	 * 1), then the chi-square probability function is a good approximation to
	 * the distribution of (14.3.1) in the case of the null hypothesis. Its use
	 * to estimate the significance of the chi-square test is standard. The
	 * appropriate value of ?, the number of degrees of freedom, bears some
	 * additional discussion. If the data are collected with the model ni’s
	 * fixed — that is, not later renormalized to fit the total observed number
	 * of events SNi — then ? equals the number of bins NB. (Note that this is
	 * not the total number of events!) Much more commonly, the ni’s are
	 * normalized after the fact so that their sum equals the sum of the Ni’s.
	 * In this case the correct value for ? is NB - 1, and the model is said to
	 * have one constraint (knstrn=1 in the program below). If the model that
	 * gives the ni’s has additional free parameters that were adjusted after
	 * the fact to agree with the data, then each of these additional “fitted”
	 * parameters decreases ? (and increases knstrn) by one additional unit.
	 */
	/**
	 * Given the array bins[1..nbins] containing the observed numbers of events, and an array 
	 * ebins[1..nbins] containing the expected numbers of events, and given the number of constraints 
	 * knstrn (normally one), this routine returns (trivially) the number of degrees of freedom 
	 * df, and (nontrivially) the chi-square chsq and the significance prob. A small value of prob 
	 * indicates a significant difference between the distributions bins and ebins. Note that bins 
	 * and ebins are both double arrays, although bins will normally contain integer values.
	 * @param bins bins
	 * @param ebins ebins
	 * @param nbins nbins
	 * @param knstrn knstrn
	 */
	public static void chsone(double[] bins, double[] ebins, int nbins,
			int knstrn)//
	// , float *df,float *chsq, float *prob)
	// Given the array bins[1..nbins] containing the observed numbers of events,
	// and an array
	// ebins[1..nbins] containing the expected numbers of events, and given the
	// number of constraints
	// knstrn (normally one), this routine returns (trivially) the number of
	// degrees of freedom
	// df, and (nontrivially) the chi-square chsq and the significance prob. A
	// small value of prob
	// indicates a significant difference between the distributions bins and
	// ebins. Note that bins
	// and ebins are both float arrays, although bins will normally contain
	// integer values.
	{
		// float gammq(float a, float x);
		// void nrerror(char error_text[]);
		failB = false;
		int j = 0;
		double temp = 0.0;
		df_chsone = nbins - knstrn;
		chsq_chsone = 0.0;
		for (j = 1; j <= nbins; j++) {
			if (ebins[j - 1] <= 0.0)// if (ebins[j] <= 0.0)
			{
				failB = true;
				failS = "Bad expected number in chsone";
				return;
				// nrerror("Bad expected number in chsone");
			}
			temp = bins[j - 1] - ebins[j - 1];// temp=bins[j]-ebins[j];
			chsq_chsone += temp * temp / ebins[j - 1];
		}
		prob_chsone = SpecFunc.gammq(0.5 * (df_chsone), 0.5 * (chsq_chsone)); // Chi-square
																				// probability
																				// function.
																				// See
																				// §6.2.
	}

	/**
	 * Given the arrays bins1[1..nbins] and bins2[1..nbins], containing two sets of binned 
	 * data, and given the number of constraints knstrn (normally 1 or 0), this routine returns the 
	 * number of degrees of freedom df, the chi-square chsq, and the significance prob. A small value 
	 * of prob indicates a significant difference between the distributions bins1 and bins2.  Note that 
	 * bins1 and bins2 are both double arrays, although they will normally contain integer values.
	 * @param bins1 bins1
	 * @param bins2 bins2
	 * @param nbins nbins
	 * @param knstrn knstrn
	 */
	public static void chstwo(double[] bins1, double[] bins2, int nbins,
			int knstrn)
	// , float *df,float *chsq, float *prob)
	// Given the arrays bins1[1..nbins] and bins2[1..nbins], containing two sets
	// of binned
	// data, and given the number of constraints knstrn (normally 1 or 0), this
	// routine returns the
	// number of degrees of freedom df, the chi-square chsq, and the
	// significance prob. A small value
	// of prob indicates a significant difference between the distributions
	// bins1 and bins2. Notethat
	// bins1 and bins2 are both float arrays, although they will normally
	// contain integer values.
	{
		// float gammq(float a, float x);
		int j = 0;
		double temp = 0.0;
		df_chstwo = nbins - knstrn;
		chsq_chstwo = 0.0;
		for (j = 1; j <= nbins; j++)
			// if (bins1[j] == 0.0 && bins2[j] == 0.0)
			if (bins1[j - 1] == 0.0 && bins2[j - 1] == 0.0)
				--(df_chstwo); // No data means one less degree of freedom.
			else {
				temp = bins1[j - 1] - bins2[j - 1];// temp=bins1[j]-bins2[j];
				// *chsq += temp*temp/(bins1[j]+bins2[j]);
				chsq_chstwo += temp * temp / (bins1[j - 1] + bins2[j - 1]);
			}
		prob_chstwo = SpecFunc.gammq(0.5 * (df_chstwo), 0.5 * (chsq_chstwo)); // Chi-square
																				// probability
																				// function.
																				// See
																				// §6.2.
	}

	/*
	 * Kolmogorov-Smirnov Test
	 * 
	 * The Kolmogorov-Smirnov (or K–S) test is applicable to unbinned
	 * distributions that are functions of a single independent variable, that
	 * is, to data sets where each data point can be associated with a single
	 * number (lifetime of each lightbulb when it burns out, or declination of
	 * each star). In such cases, the list of data points can be easily
	 * converted to an unbiased estimator SN(x) of the cumulative distribution
	 * function of the probability distribution from which it was drawn: If the
	 * N events are located at values xi, i = 1, . . .,N, then SN(x) is the
	 * function giving the fraction of data points to the left of a given value
	 * x. This function is obviously constant between consecutive (i.e., sorted
	 * into ascending order) xi’s, and jumps by the same constant 1/N at each
	 * xi. (..............
	 */
	/**
	 * Kolmogorov-Smirnov probability function. Used in ksone.
	 * @param alam alam
	 * @return the result
	 */
	public static double probks(double alam)
	// Kolmogorov-Smirnov probability function.
	{
		int j = 0;
		double a2 = 0.0;
		double fac = 2.0;
		double sum = 0.0;
		double term = 0.0;
		double termbf = 0.0;

		a2 = -2.0 * alam * alam;
		for (j = 1; j <= 100; j++) {
			term = fac * Math.exp(a2 * j * j);
			sum += term;
			if (Math.abs(term) <= EPS1 * termbf || Math.abs(term) <= EPS2 * sum)
				return sum;
			fac = -fac; // Alternating signs in sum.
			termbf = Math.abs(term);
		}
		return 1.0; // Get here only by failing to converge.
	}

	/**
	 * Given an array data[1..n], and given a user-supplied function of a single variable func which 
	 * is a cumulative distribution function ranging from 0 (for smallest values of its argument) to 1 
	 * (for largest values of its argument), this routine returns the K–S statistic d, and the significance 
	 * level prob. Small values of prob show that the cumulative distribution function of data is significantly different from func. 
	 * The array data is modified by being sorted into ascending order.
	 * @param data data
	 * @param n n
	 */
	public static void ksone(double[] data, int n)// ,
	// float (*func)(float), float *d, float *prob)
	// Given an array data[1..n], and given a user-supplied function of a single
	// variable func which
	// is a cumulative distribution function ranging from 0 (for smallest values
	// of its argument) to 1
	// (for largest values of its argument), this routine returns the K–S
	// statistic d, and the significance
	// level prob. Small values of prob showtha t the cumulative distribution
	// function of data is
	// significantly different from func. The array data is modified by being
	// sorted into ascending
	// order.
	{
		// float probks(float alam);
		// void sort(unsigned long n, float arr[]);
		int j = 0;// unsigned long j;
		double dt = 0.0;
		double en = 0.0;
		double ff = 0.0;
		double fn = 0.0;
		double fo = 0.0;
		Sorting.sort(n, data);// If the data are already sorted into ascending
								// order, then this call can be omitted.
		en = n;
		d_ksone = 0.0;
		for (j = 1; j <= n; j++) {// Loop over the sorted data points.
			fn = j / en; // Data’s c.d.f. after this step.
			ff = func.F(data[j - 1]);// ff=(*func)(data[j]); Compare to the
										// user-supplied function.
			dt = Math.max(Math.abs(fo - ff), Math.abs(fn - ff)); // Maximum
																	// distance.
			if (dt > d_ksone)
				d_ksone = dt;
			fo = fn;
		}
		en = Math.sqrt(en);
		prob_ksone = probks((en + 0.12 + 0.11 / en) * (d_ksone));// Compute
																	// significance.
	}

	/**
	 * Given an array data1[1..n1], and an array data2[1..n2], this routine returns the K–S statistic d, and the significance 
	 * level prob for the null hypothesis that the data sets are drawn from the same distribution. Small values of prob show that the 
	 * cumulative distribution function of data1 is significantly different from that of data2. The arrays data1 and data2 
	 * are modified by being sorted into ascending order.
	 * @param data1 data1
	 * @param n1 n1
	 * @param data2 data2
	 * @param n2 n2
	 */
	public static void kstwo(double[] data1, int n1, double[] data2, int n2)// ,
	// float *d, float *prob)
	// Given an array data1[1..n1], and an array data2[1..n2], this routine
	// returns the K–
	// S statistic d, and the significance level prob for the null hypothesis
	// that the data sets are
	// drawn from the same distribution. Small values of prob showtha t the
	// cumulative distribution
	// function of data1 is significantly different from that of data2. The
	// arrays data1 and data2
	// are modified by being sorted into ascending order.
	{
		// float probks(float alam);
		// void sort(unsigned long n, float arr[]);
		// unsigned long j1=1,j2=1;
		int j1 = 1;
		int j2 = 1;
		double d1 = 0.0;
		double d2 = 0.0;
		double dt = 0.0;
		double en1 = 0.0;
		double en2 = 0.0;
		double en = 0.0;
		double fn1 = 0.0;
		double fn2 = 0.0;
		Sorting.sort(n1, data1);
		Sorting.sort(n2, data2);
		en1 = n1;
		en2 = n2;
		d_kstwo = 0.0;
		while (j1 <= n1 && j2 <= n2) {// If we are not done...
										// if ((d1=data1[j1]) <= (d2=data2[j2]))
										// fn1=j1++/en1; Next step is in data1.
			if ((d1 = data1[j1 - 1]) <= (d2 = data2[j2 - 1]))
				fn1 = j1++ / en1;
			if (d2 <= d1)
				fn2 = j2++ / en2; // Next step is in data2.
			if ((dt = Math.abs(fn2 - fn1)) > d_kstwo)
				d_kstwo = dt;
		}
		en = Math.sqrt(en1 * en2 / (en1 + en2));
		prob_kstwo = probks((en + 0.12 + 0.11 / en) * (d_kstwo)); // Compute
																	// significance.
	}

	/*
	 * In this section, and the next two sections, we deal with measures of
	 * association for two distributions. The situation is this: Each data point
	 * has two or more different quantities associated with it, and wewant to
	 * knowwhether knowledge of one quantity gives us any demonstrable advantage
	 * in predicting the value of another quantity. In many cases, one variable
	 * will be an “independent” or “control” variable, and another will be a
	 * “dependent” or “measured” variable. Then, we want to know if the latter
	 * variable is in fact dependent on or associated with the former variable.
	 * If it is, we want to have some quantitative measure of the strength of
	 * the association. One often hears this loosely stated as the question of
	 * whether two variables are correlated or uncorrelated, but we will reserve
	 * those terms for a particular kind of association (linear, or at least
	 * monotonic), as discussed in §14.5 and §14.6.
	 * 
	 * Measures of Association Based on Chi-Square
	 */

	/**
	 * Given a two-dimensional contingency table in the form of an integer array nn[1..ni][1..nj], 
	 * this routine returns the chi-square chisq, the number of degrees of freedom df, the significance 
	 * level prob (small values indicating a significant association), and two measures of association, 
	 * Cramer’s V (cramrv) and the contingency coefficient C (ccc).
	 * @param nn nn
	 * @param ni ni
	 * @param nj nj
	 */
	public static void cntab1(int[][] nn, int ni, int nj)
	// , float *chisq, float *df, float *prob, float *cramrv, float *ccc)
	// Given a two-dimensional contingency table in the form of an integer array
	// nn[1..ni][1..nj],
	// this routine returns the chi-square chisq, the number of degrees of
	// freedom df, the significance
	// level prob (small values indicating a significant association), and two
	// measures of association,
	// Cramer’s V (cramrv) and the contingency coefficient C (ccc).
	{
		// float gammq(float a, float x);
		int nnj = 0;
		int nni = 0;
		int j = 0;
		int i = 0;
		int minij = 0;
		double sum = 0.0;
		double expctd = 0.0;
		double[] sumi = new double[ni];
		double[] sumj = new double[nj];
		double temp = 0.0;
		// sumi=vector(1,ni);
		// sumj=vector(1,nj);
		nni = ni;// Number of rows
		nnj = nj;// and columns.
		for (i = 1; i <= ni; i++) {// Get the row totals.
			sumi[i - 1] = 0.0;// sumi[i]=0.0;
			for (j = 1; j <= nj; j++) {
				sumi[i - 1] += nn[i - 1][j - 1];// sumi[i] += nn[i][j];
				sum += nn[i - 1][j - 1];// sum += nn[i][j];
			}
			// if (sumi[i] == 0.0) --nni; Eliminate any zero rows by reducing
			// the number.
			if (sumi[i - 1] == 0.0)
				--nni;
		}
		for (j = 1; j <= nj; j++) {// Get the column totals.
			sumj[j - 1] = 0.0;// sumj[j]=0.0;
			for (i = 1; i <= ni; i++)
				sumj[j - 1] += nn[i - 1][j - 1];// sumj[j] += nn[i][j];
			// if (sumj[j] == 0.0) --nnj; Eliminate any zero columns.
			if (sumj[j - 1] == 0.0)
				--nnj;
		}
		df_cntab1 = nni * nnj - nni - nnj + 1;// Corrected number of degrees of
												// freedom.
		chisq_cntab1 = 0.0;
		for (i = 1; i <= ni; i++) {// Do the chi-square sum.
			for (j = 1; j <= nj; j++) {
				expctd = sumj[j - 1] * sumi[i - 1] / sum;// expctd=sumj[j]*sumi[i]/sum;
				temp = nn[i - 1][j - 1] - expctd;// temp=nn[i][j]-expctd;
				chisq_cntab1 += temp * temp / (expctd + TINY);
				// Here TINY guarantees that any eliminated row or column will
				// not contribute to the sum.
			}
		}
		prob_cntab1 = SpecFunc.gammq(0.5 * (df_cntab1), 0.5 * (chisq_cntab1)); // Chi-square
																				// probability
																				// function.
		minij = nni < nnj ? nni - 1 : nnj - 1;
		cramrv_cntab1 = Math.sqrt(chisq_cntab1 / (sum * minij));
		ccc_cntab1 = Math.sqrt(chisq_cntab1 / (chisq_cntab1 + sum));
		// free_vector(sumj,1,nj);
		// free_vector(sumi,1,ni);
	}

	/**
	 * Given a two-dimensional contingency table in the form of an integer array nn[i][j], where i 
	 * labels the x variable and ranges from 1 to ni, j labels the y variable and ranges from 1 to nj, 
	 * this routine returns the entropy h of the whole table, the entropy hx of the x distribution, the 
	 * entropy hy of the y distribution, the entropy hygx of y given x, the entropy hxgy of x given y, 
	 * the dependency uygx of y on x, the dependency uxgy of x on y, and the symmetrical dependency uxy.
	 * @param nn nn
	 * @param ni ni
	 * @param nj nj
	 */
	public static void cntab2(int[][] nn, int ni, int nj)
	// , float *h, float *hx, float *hy, float *hygx, float *hxgy, float *uygx,
	// float *uxgy, float *uxy)
	// Given a two-dimensional contingency table in the form of an integer array
	// nn[i][j], where i
	// labels the x variable and ranges from 1 to ni, j labels the y variable
	// and ranges from 1 to nj,
	// this routine returns the entropy h of the whole table, the entropy hx of
	// the x distribution, the
	// entropy hy of the y distribution, the entropy hygx of y given x, the
	// entropy hxgy of x given y,
	// the dependency uygx of y on x (eq. 14.4.15), the dependency uxgy of x on
	// y (eq. 14.4.16),
	// and the symmetrical dependency uxy (eq. 14.4.17).
	{
		int i = 0;
		int j = 0;
		double sum = 0.0;
		double p = 0.0;
		double[] sumi = new double[ni];
		double[] sumj = new double[nj];
		// sumi=vector(1,ni);
		// sumj=vector(1,nj);
		for (i = 1; i <= ni; i++) {// Get the row totals.
			sumi[i - 1] = 0.0;// sumi[i]=0.0;
			for (j = 1; j <= nj; j++) {
				sumi[i - 1] += nn[i - 1][j - 1];// sumi[i] += nn[i][j];
				sum += nn[i - 1][j - 1];// sum += nn[i][j];
			}
		}
		for (j = 1; j <= nj; j++) {// Get the column totals.
			sumj[j - 1] = 0.0;// sumj[j]=0.0;
			for (i = 1; i <= ni; i++)
				sumj[j - 1] += nn[i - 1][j - 1];// sumj[j] += nn[i][j];
		}
		hx_cntab2 = 0.0;// Entropy of the x distribution,
		for (i = 1; i <= ni; i++)
			if (sumi[i - 1] != 0.0) // if (sumi[i])
			{
				p = sumi[i - 1] / sum;// p=sumi[i]/sum;
				hx_cntab2 -= p * Math.log(p);
			}
		hy_cntab2 = 0.0;// and of the y distribution.
		for (j = 1; j <= nj; j++)
			if (sumj[j - 1] != 0.0) // if (sumj[j])
			{
				p = sumj[j - 1] / sum;// p=sumj[j]/sum;
				hy_cntab2 -= p * Math.log(p);
			}
		h_cntab2 = 0.0;
		for (i = 1; i <= ni; i++)
			// Total entropy: loop over both x
			for (j = 1; j <= nj; j++)
				// and y.
				if (nn[i - 1][j - 1] != 0.0) // if (nn[i][j])
				{
					p = nn[i - 1][j - 1] / sum;// p=nn[i][j]/sum;
					h_cntab2 -= p * Math.log(p);
				}
		hygx_cntab2 = (h_cntab2) - (hx_cntab2); // Uses equation (14.4.18),
		hxgy_cntab2 = (h_cntab2) - (hy_cntab2); // as does this.
		uygx_cntab2 = (hy_cntab2 - hygx_cntab2) / (hy_cntab2 + TINY);// Equation
																		// (14.4.15).
		uxgy_cntab2 = (hx_cntab2 - hxgy_cntab2) / (hx_cntab2 + TINY);// Equation
																		// (14.4.16).
		uxy_cntab2 = 2.0 * (hx_cntab2 + hy_cntab2 - h_cntab2)
				/ (hx_cntab2 + hy_cntab2 + TINY);// Equation (14.4.17).
		// free_vector(sumj,1,nj);
		// free_vector(sumi,1,ni);
	}

	/**
	 * Given two arrays x[1..n] and y[1..n], this routine computes their correlation coefficient 
	 * r (returned as r_pearsn), the significance level at which the null hypothesis of zero correlation is 
	 * disproved (prob whose small value indicates a significant correlation), and Fisher’s z (returned 
	 * as z_pearsn), whose value can be used in further statistical tests as described above.
	 * @param x x
	 * @param y y
	 * @param n n
	 */
	public static void pearsn(double[] x, double[] y, int n)
	// float *r, float *prob, float *z)
	// Given two arrays x[1..n] and y[1..n], this routine computes their
	// correlation coefficient
	// r (returned as r), the significance level at which the null hypothesis of
	// zero correlation is
	// disproved (prob whose small value indicates a significant correlation),
	// and Fisher’s z (returned
	// as z), whose value can be used in further statistical tests as described
	// above.
	{
		// float betai(float a, float b, float x);
		// float erfcc(float x);
		int j = 0;// unsigned long j;
		double yt = 0.0;
		double xt = 0.0;
		double t = 0.0;
		double df = 0.0;
		double syy = 0.0;
		double sxy = 0.0;
		double sxx = 0.0;
		double ay = 0.0;
		double ax = 0.0;
		for (j = 1; j <= n; j++) {// Find the means.
			ax += x[j - 1];// ax += x[j];
			ay += y[j - 1];// ay += y[j];
		}
		ax /= n;
		ay /= n;
		for (j = 1; j <= n; j++) {// Compute the correlation coefficient.
			xt = x[j - 1] - ax;// xt=x[j]-ax;
			yt = y[j - 1] - ay;// yt=y[j]-ay;
			sxx += xt * xt;
			syy += yt * yt;
			sxy += xt * yt;
		}
		r_pearsn = sxy / (Math.sqrt(sxx * syy) + TINY2);
		z_pearsn = 0.5 * Math.log((1.0 + (r_pearsn) + TINY2)
				/ (1.0 - (r_pearsn) + TINY2)); // Fisher’s z transformation.
		df = n - 2;
		t = (r_pearsn)
				* Math.sqrt(df
						/ ((1.0 - (r_pearsn) + TINY2) * (1.0 + (r_pearsn) + TINY2)));// Equation
																						// (14.5.5).
		prob_pearsn = SpecFunc.betai(0.5 * df, 0.5, df / (df + t * t)); // Student’s
																		// t
																		// probability.
		/* *prob=erfcc(fabs((*z)*sqrt(n-1.0))/1.4142136) */
		// For large n, this easier computation of prob, using the short routine
		// erfcc, would give approximately
		// the same value.
	}

	/**
	 * Given a sorted array w[1..n], replaces the elements by their rank, including midranking of ties, 
	 * and returns as s_crank the sum of f3 - f, w here f is the number of elements in each tie.
	 * @param n n
	 * @param w w
	 */
	public static void crank(int n, double[] w)// , float *s)
	// Given a sorted array w[1..n], replaces the elements by their rank,
	// including midranking of ties,
	// and returns as s the sum of f3 - f, w here f is the number of elements in
	// each tie.
	{
		int j = 1;
		int ji = 0;
		int jt = 0;
		double t = 0.0;
		double rank = 0.0;
		s_crank = 0.0;
		while (j < n) {
			if (w[j] != w[j - 1])// if (w[j+1] != w[j])
			{// Not a tie.
				w[j - 1] = j;// w[j]=j;
				++j;
			} else {// A tie:
					// for (jt=j+1;jt<=n && w[jt]==w[j];jt++); Howf ar does it
					// go?
				for (jt = j + 1; jt <= n && w[jt - 1] == w[j - 1]; jt++)
					;
				rank = 0.5 * (j + jt - 1);// This is the mean rank of the tie,
				for (ji = j; ji <= (jt - 1); ji++)
					w[ji - 1] = rank;// w[ji]=rank; so enter it into all the
										// tied entries,
				t = jt - j;
				s_crank += t * t * t - t; // and update s.
				j = jt;
			}
		}
		if (j == n)
			w[n - 1] = n;// w[n]=n; //If the last element was not tied, this is
							// its rank.
	}

	/**
	 * Given two data arrays, data1[1..n] and data2[1..n], this routine returns their sum-squared 
	 * difference of ranks as d_spear, the number of standard deviations by which d_spear deviates from its null hypothesis 
	 * expected value as zd_spear, the two-sided significance level of this deviation as probd_spear, Spearman’s rank correlation 
	 * rs as rs_spear, and the two-sided significance level of its deviation from zero as probrs_spear. The external routines crank 
	 * and sort2 (see Sorting) are used. A small value of either probd_spear or probrs_spear indicates a significant correlation (rs 
	 * positive) or anticorrelation (rs negative).
	 * @param data1 data1
	 * @param data2 data2
	 * @param n n
	 */
	public static void spear(double[] data1, double[] data2, int n)
	// , float *d, float *zd, float *probd, float *rs, float *probrs)
	// Given two data arrays, data1[1..n] and data2[1..n], this routine returns
	// their sum-squared
	// difference of ranks as D, the number of standard deviations by which D
	// deviates from its nullhypothesis
	// expected value as zd, the two-sided significance level of this deviation
	// as probd,
	// Spearman’s rank correlation rs as rs, and the two-sided significance
	// level of its deviation from
	// zero as probrs. The external routines crank (below) and sort2 (§8.2) are
	// used. A small value
	// of either probd or probrs indicates a significant correlation (rs
	// positive) or anticorrelation
	// (rs negative).
	{
		// float betai(float a, float b, float x);
		// void crank(unsigned long n, float w[], float *s);
		// float erfcc(float x);
		// void sort2(unsigned long n, float arr[], float brr[]);
		int j = 0;// unsigned long j;
		double vard = 0.0;
		double t = 0.0;
		double sg = 0.0;
		double sf = 0.0;
		double fac = 0.0;
		double en3n = 0.0;
		double en = 0.0;
		double df = 0.0;
		double aved = 0.0;
		double[] wksp1 = new double[n];
		double[] wksp2 = new double[n];
		// wksp1=vector(1,n);
		// wksp2=vector(1,n);
		for (j = 1; j <= n; j++) {
			wksp1[j - 1] = data1[j - 1];// wksp1[j]=data1[j];
			wksp2[j - 1] = data2[j - 1];// wksp2[j]=data2[j];
		}
		Sorting.sort2(n, wksp1, wksp2); // Sort each of the data arrays, and
										// convert the entries to
		// ranks. The values sf and sg return the sums (f3k-fk) and (g3m -
		// gm), respectively.
		crank(n, wksp1);// ,&sf);
		sf = s_crank;
		Sorting.sort2(n, wksp2, wksp1);
		crank(n, wksp2);// ,&sg);
		sg = s_crank;

		d_spear = 0.0;
		for (j = 1; j <= n; j++)
			// Sum the squared difference of ranks.
			d_spear += (wksp1[j - 1] - wksp2[j - 1])
					* (wksp1[j - 1] - wksp2[j - 1]);// SQR(wksp1[j]-wksp2[j]);
		en = n;
		en3n = en * en * en - en;
		aved = en3n / 6.0 - (sf + sg) / 12.0;// Expectation value of D,
		fac = (1.0 - sf / en3n) * (1.0 - sg / en3n);
		vard = ((en - 1.0) * en * en * (en + 1.0) * (en + 1.0) / 36.0) * fac; // and
																				// variance
																				// of
																				// D
																				// give
		zd_spear = (d_spear - aved) / Math.sqrt(vard); // number of standard
														// deviations and
														// significance.
		probd_spear = SpecFunc.erfcc(Math.abs(zd_spear) / 1.4142136);
		rs_spear = (1.0 - (6.0 / en3n) * (d_spear + (sf + sg) / 12.0))
				/ Math.sqrt(fac); // Rank correlation coefficient,
		fac = (rs_spear + 1.0) * (1.0 - (rs_spear));
		if (fac > 0.0) {
			t = (rs_spear) * Math.sqrt((en - 2.0) / fac); // and its t value,
			df = en - 2.0;
			probrs_spear = SpecFunc.betai(0.5 * df, 0.5, df / (df + t * t)); // sgive
																				// its
																				// significance.
		} else
			probrs_spear = 0.0;
		// free_vector(wksp2,1,n);
		// free_vector(wksp1,1,n);
	}

	/*
	 * On balance, we prefer r s as being the more straightforward nonparametric
	 * test, but both statistics are in general use. In fact, t and rs are very
	 * strongly correlated and, in most applications, are effectively the same
	 * test.
	 */
	/**
	 * Given data arrays data1[1..n] and data2[1..n], this program returns Kendall’s t as tau_kendl1, 
	 * its number of standard deviations from zero as z_kendl1, and its two-sided significance level as prob_kendl1. 
	 * Small values of prob_kendl1 indicate a significant correlation (tau positive) or anticorrelation (tau negative).
	 * @param data1 data1
	 * @param data2 data2
	 * @param n n
	 */
	public static void kendl1(double[] data1, double[] data2, int n)
	// , float *tau, float *z, float *prob)
	// ..Given data arrays data1[1..n] and data2[1..n], this program returns
	// Kendall’s t as tau,
	// its number of standard deviations from zero as z, and its two-sided
	// significance level as prob.
	// Small values of prob indicate a significant correlation (tau positive) or
	// anticorrelation (tau
	// negative).
	{
		// float erfcc(float x);
		int n2 = 0;
		int n1 = 0;
		int k = 0;
		int j = 0;
		int is = 0;
		double svar = 0.0;
		double aa = 0.0;
		double a2 = 0.0;
		double a1 = 0.0;
		for (j = 1; j < n; j++) {// Loop over first member of pair,
			for (k = (j + 1); k <= n; k++) {// and second member.
				a1 = data1[j - 1] - data1[k - 1];// a1=data1[j]-data1[k];
				a2 = data2[j - 1] - data2[k - 1];// a2=data2[j]-data2[k];
				aa = a1 * a2;
				if (aa != 0.0)// (aa)
				{// Neither array has a tie.
					++n1;
					++n2;
					// aa > 0.0 ? ++is : --is;
					if (aa > 0.0)
						++is;
					else
						--is;
				} else {// One or both arrays have ties.
					if (a1 != 0.0)
						++n1; // An “extra x” event.
					if (a2 != 0.0)
						++n2; // An “extra y” event.
				}
			}
		}
		tau_kendl1 = is / (Math.sqrt((double) n1) * Math.sqrt((double) n2)); // Equation
																				// (14.6.8).
		svar = (4.0 * n + 10.0) / (9.0 * n * (n - 1.0)); // Equation (14.6.9).
		z_kendl1 = (tau_kendl1) / Math.sqrt(svar);
		prob_kendl1 = SpecFunc.erfcc(Math.abs(z_kendl1) / 1.4142136); // Significance.
	}

	/**
	 * Given a two-dimensional table tab[1..i][1..j], such that tab[k][l] contains the number 
	 * of events falling in bin k of one variable and bin l of another, this program returns Kendall’s t 
	 * as tau_kendl2, its number of standard deviations from zero as z_kendl2, and its two-sided significance level as 
	 * prob_kendl2. Small values of prob_kendl2 indicate a significant correlation (tau positive) or anticorrelation 
	 * (tau negative) between the two variables. Although tab is a double array, it will normally contain integral values.
	 * @param tab tab
	 * @param i i
	 * @param j j
	 */
	public static void kendl2(double[][] tab, int i, int j)// , float *tau,
															// float *z, float
															// *prob)
	// Given a two-dimensional table tab[1..i][1..j], such that tab[k][l]
	// contains the number
	// of events falling in bin k of one variable and bin l of another, this
	// program returns Kendall’s t
	// as tau_kendl2, its number of standard deviations from zero as z_kendl2, and its
	// two-sided significance level as
	// prob. Small values of prob indicate a significant correlation (tau
	// positive) or anticorrelation
	// (tau negative) between the two variables. Although tab is a float array,
	// it will normally
	// contain integral values.
	{
		// float erfcc(float x);
		int nn = 0;
		int mm = 0;
		int m2 = 0;
		int m1 = 0;
		int lj = 0;
		int li = 0;
		int l = 0;
		int kj = 0;
		int ki = 0;
		int k = 0;

		double svar = 0.0;
		double s = 0.0;
		double points = 0.0;
		double pairs = 0.0;
		double en2 = 0.0;
		double en1 = 0.0;

		nn = i * j; // Total number of entries in contingency table.
		points = tab[i - 1][j - 1];// points=tab[i][j];
		for (k = 0; k <= nn - 2; k++) {// Loop over entries in table,
			ki = (k / j);// decoding a row,
			kj = k - j * ki;// and a column.
			points += tab[ki][kj];// tab[ki+1][kj+1]; Increment the total count
									// of events.
			for (l = k + 1; l <= nn - 1; l++) {// Loop over other member of the
												// pair,
				li = l / j;// decoding its row
				lj = l - j * li;// and column.
				mm = (m1 = li - ki) * (m2 = lj - kj);
				// pairs=tab[ki+1][kj+1]*tab[li+1][lj+1];
				pairs = tab[ki][kj] * tab[li][lj];
				if (mm != 0)// if (mm)
				{// Not a tie.
					en1 += pairs;
					en2 += pairs;
					s += (mm > 0 ? pairs : -pairs);// Concordant, or discordant.
				} else {
					if (m1 != 0)
						en1 += pairs;
					if (m2 != 0)
						en2 += pairs;
				}
			}
		}
		tau_kendl2 = s / Math.sqrt(en1 * en2);
		svar = (4.0 * points + 10.0) / (9.0 * points * (points - 1.0));
		z_kendl2 = (tau_kendl2) / Math.sqrt(svar);
		prob_kendl2 = SpecFunc.erfcc(Math.abs(z_kendl2) / 1.4142136);
	}

	/*
	 * Do Two-Dimensional Distributions Differ?
	 */
	/**
	 * Two-dimensional Kolmogorov-Smirnov test of one sample against a model. Given the x and y 
	 * coordinates of n1 data points in arrays x1[1..n1] and y1[1..n1], and given a user-supplied 
	 * function quadvl that exemplifies the model, this routine returns the two-dimensional K-S 
	 * statistic as d1_ks2d1s, and its significance level as prob_ks2d1s. Small values of prob_ks2d1s show that the sample 
	 * is significantly different from the model. Note that the test is slightly distribution-dependent, 
	 * so prob_ks2d1s is only an estimate.
	 * @param x1 x1
	 * @param y1 y1
	 * @param n1 n
	 */
	public static void ks2d1s(double[] x1, double[] y1, int n1)// ,
	// void (*quadvl)(float, float, float *, float *, float *, float *),float
	// *d1, float *prob)
	// Two-dimensional Kolmogorov-Smirnov test of one sample against a model.
	// Given the x and y
	// coordinates of n1 data points in arrays x1[1..n1] and y1[1..n1], and
	// given a user-supplied
	// function quadvl that exemplifies the model, this routine returns the
	// two-dimensional K-S
	// statistic as d1, and its significance level as prob. Small values of prob
	// show that the sample
	// is significantly different from the model. Note that the test is slightly
	// distribution-dependent,
	// so prob is only an estimate.
	{
		// void pearsn(float x[], float y[], unsigned long n, float *r, float
		// *prob,float *z);
		// float probks(float alam);
		// void quadct(float x, float y, float xx[], float yy[], unsigned long
		// nn,
		// float *fa, float *fb, float *fc, float *fd);

		int j = 0;
		// double dum = 0.0;
		// double dumm = 0.0;
		double fa = 0.0;
		double fb = 0.0;
		double fc = 0.0;
		double fd = 0.0;
		double ga = 0.0;
		double gb = 0.0;
		double gc = 0.0;
		double gd = 0.0;
		double r1 = 0.0;
		double rr = 0.0;
		double sqen = 0.0;

		d1_ks2d1s = 0.0;
		for (j = 1; j <= n1; j++) {// Loop over the data points.
									// quadct(x1[j],y1[j],x1,y1,n1,&fa,&fb,&fc,&fd);
			quadct(x1[j], y1[j], x1, y1, n1);
			fa = fa_quadct;
			fb = fb_quadct;
			fc = fc_quadct;
			fd = fd_quadct;
			// (*quadvl)(x1[j],y1[j],&ga,&gb,&gc,&gd);
			QDVL(x1[j], y1[j], "quadvl");
			ga = fa_quadvl;
			gb = fb_quadvl;
			gc = fc_quadvl;
			gd = fd_quadvl;

			d1_ks2d1s = Math.max(d1_ks2d1s, Math.abs(fa - ga));
			d1_ks2d1s = Math.max(d1_ks2d1s, Math.abs(fb - gb));
			d1_ks2d1s = Math.max(d1_ks2d1s, Math.abs(fc - gc));
			d1_ks2d1s = Math.max(d1_ks2d1s, Math.abs(fd - gd));
			// For both the sample and the model, the distribution is integrated
			// in each of four
			// quadrants, and the maximum difference is saved.
		}
		// pearsn(x1,y1,n1,&r1,&dum,&dumm); Get the linear correlation
		// coefficient r1.
		pearsn(x1, y1, n1);// ,&r1,&dum,&dumm);
		r1 = r_pearsn;
		// dum = prob_pearsn;
		// dumm = z_pearsn;

		sqen = Math.sqrt((double) n1);
		rr = Math.sqrt(1.0 - r1 * r1);
		// Estimate the probability using the K-S probability function probks.
		prob_ks2d1s = probks(d1_ks2d1s * sqen
				/ (1.0 + rr * (0.25 - 0.75 / sqen)));
	}

	/**
	 * Given an origin (x, y), and an array of nn points with coordinates xx[1..nn] and yy[1..nn], 
	 * count how many of them are in each quadrant around the origin, and return the normalized 
	 * fractions. Quadrants are labeled alphabetically, counterclockwise from the upper right. Used 
	 * by ks2d1s and ks2d2s.
	 * @param x x
	 * @param y y
	 * @param xx xx
	 * @param yy yy
	 * @param nn nn
	 */
	public static void quadct(double x, double y, double[] xx, double[] yy,
			int nn)// ,
	// float *fa, float *fb, float *fc, float *fd)
	// Given an origin (x, y), and an array of nn points with coordinates
	// xx[1..nn] and yy[1..nn],
	// count how many of them are in each quadrant around the origin, and return
	// the normalized
	// fractions. Quadrants are labeled alphabetically, counterclockwise from
	// the upper right. Used
	// by ks2d1s and ks2d2s.
	{
		int k = 0;
		int na = 0;
		int nb = 0;
		int nc = 0;
		int nd = 0;
		double ff = 0.0;

		na = nb = nc = nd = 0;
		for (k = 1; k <= nn; k++) {
			if (yy[k - 1] > y) // if (yy[k] > y)
			{
				// xx[k] > x ? ++na : ++nb;
				if (xx[k - 1] > x)
					++na;
				else
					++nb;
			} else {
				// xx[k] > x ? ++nd : ++nc;
				if (xx[k - 1] > x)
					++nd;
				else
					++nc;
			}
		}
		ff = 1.0 / nn;
		fa_quadct = ff * na;
		fb_quadct = ff * nb;
		fc_quadct = ff * nc;
		fd_quadct = ff * nd;
	}

	/**
	 * Internal use.
	 * @param x x
	 * @param y y
	 * @param namS namS
	 */
	public static void QDVL(double x, double y, String namS) {
		if (namS.compareTo("quadvl") == 0)
			quadvl(x, y);
	}

	/**
	 * User-supplied function for ks2
	 * @param x x
	 * @param y y
	 */
	public static void quadvl(double x, double y)//
	// , float *fa, float *fb, float *fc, float *fd)
	// This is a sample of a user-supplied routine to be used with ks2d1s. In
	// this case, the model
	// distribution is uniform inside the square -1 < x < 1, -1 < y < 1. In
	// general this routine
	// should return, for any point (x, y), the fraction of the total
	// distribution in each of the four
	// quadrants around that point. The fractions, fa, fb, fc, and fd, must add
	// up to 1. Quadrants
	// are alphabetical, counterclockwise from the upper right.
	{
		double qa = 0.0;
		double qb = 0.0;
		double qc = 0.0;
		double qd = 0.0;
		qa = Math.min(2.0, Math.max(0.0, 1.0 - x));
		qb = Math.min(2.0, Math.max(0.0, 1.0 - y));
		qc = Math.min(2.0, Math.max(0.0, x + 1.0));
		qd = Math.min(2.0, Math.max(0.0, y + 1.0));
		fa_quadvl = 0.25 * qa * qb;
		fb_quadvl = 0.25 * qb * qc;
		fc_quadvl = 0.25 * qc * qd;
		fd_quadvl = 0.25 * qd * qa;
	}

	/**
	 * Two-dimensional Kolmogorov-Smirnov test on two samples. Given the x and y coordinates of 
	 * the first sample as n1 values in arrays x1[1..n1] and y1[1..n1], and likewise for the second 
	 * sample, n2 values in arrays x2 and y2, this routine returns the two-dimensional, two-sample 
	 * K-S statistic as d_ks2d2s, and its significance level as prob_ks2d2s. Small values of prob_ks2d2s show that the 
	 * two samples are significantly different. Note that the test is slightly distribution-dependent, so 
	 * prob_ks2d2s is only an estimate.
	 * @param x1 x1
	 * @param y1 y1
	 * @param n1 n1
	 * @param x2 x2
	 * @param y2 y2
	 * @param n2 n2
	 */
	public static void ks2d2s(double[] x1, double[] y1, int n1, double[] x2,
			double[] y2, int n2)
	// , float *d, float *prob)
	// Two-dimensional Kolmogorov-Smirnov test on two samples. Given the x and y
	// coordinates of
	// the first sample as n1 values in arrays x1[1..n1] and y1[1..n1], and
	// likewise for the second
	// sample, n2 values in arrays x2 and y2, this routine returns the
	// two-dimensional, two-sample
	// K-S statistic as d, and its significance level as prob. Small values of
	// prob show that the
	// two samples are significantly different. Note that the test is slightly
	// distribution-dependent, so
	// prob is only an estimate.
	{
		// void pearsn(float x[], float y[], unsigned long n, float *r, float
		// *prob, float *z);
		// float probks(float alam);
		// void quadct(float x, float y, float xx[], float yy[], unsigned long
		// nn,
		// float *fa, float *fb, float *fc, float *fd);

		int j = 0;
		double d1 = 0.0;
		double d2 = 0.0;
		// double dum = 0.0;
		// double dumm = 0.0;
		double fa = 0.0;
		double fb = 0.0;
		double fc = 0.0;
		double fd = 0.0;
		double ga = 0.0;
		double gb = 0.0;
		double gc = 0.0;
		double gd = 0.0;
		double r1 = 0.0;
		double r2 = 0.0;
		double rr = 0.0;
		double sqen = 0.0;

		d1 = 0.0;
		for (j = 1; j <= n1; j++) {// First, use points in the first sample as
									// origins.
			quadct(x1[j], y1[j], x1, y1, n1);// ,&fa,&fb,&fc,&fd);
			fa = fa_quadct;
			fb = fb_quadct;
			fc = fc_quadct;
			fd = fd_quadct;
			quadct(x1[j], y1[j], x2, y2, n2);// ,&ga,&gb,&gc,&gd);
			ga = fa_quadct;
			gb = fb_quadct;
			gc = fc_quadct;
			gd = fd_quadct;

			d1 = Math.max(d1, Math.abs(fa - ga));
			d1 = Math.max(d1, Math.abs(fb - gb));
			d1 = Math.max(d1, Math.abs(fc - gc));
			d1 = Math.max(d1, Math.abs(fd - gd));
		}
		d2 = 0.0;
		for (j = 1; j <= n2; j++) {// Then, use points in the second sample as
									// origins.
			quadct(x2[j], y2[j], x1, y1, n1);// ,&fa,&fb,&fc,&fd);
			fa = fa_quadct;
			fb = fb_quadct;
			fc = fc_quadct;
			fd = fd_quadct;
			quadct(x2[j], y2[j], x2, y2, n2);// ,&ga,&gb,&gc,&gd);
			ga = fa_quadct;
			gb = fb_quadct;
			gc = fc_quadct;
			gd = fd_quadct;

			d2 = Math.max(d2, Math.abs(fa - ga));
			d2 = Math.max(d2, Math.abs(fb - gb));
			d2 = Math.max(d2, Math.abs(fc - gc));
			d2 = Math.max(d2, Math.abs(fd - gd));
		}
		d_ks2d2s = 0.5 * (d1 + d2);// Average the K-S statistics.
		sqen = Math.sqrt(n1 * n2 / (double) (n1 + n2));
		pearsn(x1, y1, n1);// ,&r1,&dum,&dumm); //Get the linear correlation
							// coefficient for each sample.
		r1 = r_pearsn;
		// dum = prob_pearsn;
		// dumm = z_pearsn;
		pearsn(x2, y2, n2);// ,&r2,&dum,&dumm);
		r2 = r_pearsn;
		// dum = prob_pearsn;
		// dumm = z_pearsn;
		rr = Math.sqrt(1.0 - 0.5 * (r1 * r1 + r2 * r2));
		// Estimate the probability using the K-S probability function probks.
		prob_ks2d2s = probks(d_ks2d2s * sqen
				/ (1.0 + rr * (0.25 - 0.75 / sqen)));
	}

	/*
	 * Savitzky-Golay Smoothing Filters
	 * 
	 * In §13.5 we learned something about the construction and application of
	 * digital filters, but little guidance was given on which particular filter
	 * to use. That, of course, depends on what you want to accomplish by
	 * filtering. One obvious use for low-pass filters is to smooth noisy data.
	 * The premise of data smoothing is that one is measuring a variable that is
	 * both slowly varying and also corrupted by random noise. Then it can
	 * sometimes be useful to replace each data point by some kind of local
	 * average of surrounding data points. Since nearby points measure very
	 * nearly the same underlying value, averaging can reduce the level of noise
	 * without (much) biasing the value obtained. We must comment editorially
	 * that the smoothing of data lies in a murky area, beyond the fringe of
	 * some better posed, and therefore more highly recommended, techniques that
	 * are discussed elsewhere in this book. If you are fitting data to a
	 * parametric model, for example (see Chapter 15), it is almost always
	 * better to use raw data than to use data that has been pre-processed by a
	 * smoothing procedure. Another alternative to blind smoothing is so-called
	 * “optimal” or Wiener filtering, as discussed in §13.3 and more generally
	 * in §13.6. Data smoothing is probably most justified when it is used
	 * simply as a graphical technique, to guide the eye through a forest of
	 * data points all with large error bars; or as a means of making initial
	 * rough estimates of simple parameters from a graph. In this section we
	 * discuss a particular type of low-pass filter, well-adapted for data
	 * smoothing, and termed variously Savitzky-Golay [1], least-squares [2], or
	 * DISPO (Digital Smoothing Polynomial) [3] filters. Rather than having
	 * their properties defined in the Fourier domain, and then translated to
	 * the time domain, Savitzky-Golay filters derive directly from a particular
	 * formulation of the data smoothing problem in the time domain, as we will
	 * now see. Savitzky-Golay filters were initially (and are still often) used
	 * to render visible the relative widths and heights of spectral lines in
	 * noisy spectrometric data.
	 */
	/**
	 * Returns in c[1..np], in wrap-around order consistent with the argument respns in routine convlv, a set of 
	 * Savitzky-Golay filter coefficients. nl is the number of leftward (past) 
	 * data points used, while nr is the number of rightward (future) data points, making the total 
	 * number of data points used nl+nr+1. ld is the order of the derivative desired (e.g., ld = 0 
	 * for smoothed function). m is the order of the smoothing polynomial, also equal to the highest 
	 * conserved moment; usual values are m = 2 or m = 4.
	 * @param c c
	 * @param np np
	 * @param nl nl
	 * @param nr nr
	 * @param ld ld
	 * @param m m
	 */
	public static void savgol(double[] c, int np, int nl, int nr, int ld, int m)
	// Returns in c[1..np], in wrap-around order (N.B.!) consistent with the
	// argument respns in
	// routine convlv, a set of Savitzky-Golay filter coefficients. nl is the
	// number of leftward (past)
	// data points used, while nr is the number of rightward (future) data
	// points, making the total
	// number of data points used nl+nr+1. ld is the order of the derivative
	// desired (e.g., ld = 0
	// for smoothed function). m is the order of the smoothing polynomial, also
	// equal to the highest
	// conserved moment; usual values are m = 2 or m = 4.
	{
		// void lubksb(float **a, int n, int *indx, float b[]);
		// void ludcmp(float **a, int n, int *indx, float *d);
		int imj = 0;
		int ipj = 0;
		int j = 0;
		int k = 0;
		int kk = 0;
		int mm = 0;
		int[] indx = new int[m + 1];
		// double d = 0.0;
		double fac = 0.0;
		double sum = 0.0;
		double[][] a = new double[m + 1][m + 1];
		double[] b = new double[m + 1];
		failB = false;
		if (np < nl + nr + 1 || nl < 0 || nr < 0 || ld > m || nl + nr < m) {
			failB = true;
			failS = "bad args in savgol";
			return;
			// nrerror("bad args in savgol");
		}
		// indx=ivector(1,m+1);
		// a=matrix(1,m+1,1,m+1);
		// b=vector(1,m+1);
		for (ipj = 0; ipj <= (m << 1); ipj++) {// Set up the normal equations of
												// the desired least-squares
												// fit.
			sum = (ipj != 0 ? 0.0 : 1.0);// sum=(ipj ? 0.0 : 1.0);
			for (k = 1; k <= nr; k++)
				sum += Math.pow((double) k, (double) ipj);
			for (k = 1; k <= nl; k++)
				sum += Math.pow((double) -k, (double) ipj);
			mm = Math.min(ipj, 2 * m - ipj);
			for (imj = -mm; imj <= mm; imj += 2)
				// a[1+(ipj+imj)/2][1+(ipj-imj)/2]=sum;
				a[(ipj + imj) / 2][(ipj - imj) / 2] = sum;
		}
		// ludcmp(a,m+1,indx,&d); Solve them: LU decomposition.
		LinAEq.ludcmp(a, m + 1, indx);// ,&d);
		for (j = 1; j <= m + 1; j++)
			b[j - 1] = 0.0;// b[j]=0.0;
		b[ld] = 1.0;// b[ld+1]=1.0;
		// Right-hand side vector is unit vector, depending on which derivative
		// we want.
		LinAEq.lubksb(a, m + 1, indx, b);// Get one row of the inverse matrix.
		for (kk = 1; kk <= np; kk++)
			c[kk - 1] = 0.0;// c[kk]=0.0; Zero the output array (it may be
							// bigger than
							// number of coefficients).
		for (k = -nl; k <= nr; k++) {
			sum = b[0];// b[1]; Each Savitzky-Golay coefficient is the dot
						// product of powers of an integer with the inverse
						// matrix row.
			fac = 1.0;
			for (mm = 1; mm <= m; mm++)
				sum += b[mm] * (fac *= k);// sum += b[mm+1]*(fac *= k);
			kk = ((np - k) % np) + 1;// Store in wrap-around order.
			c[kk - 1] = sum;// c[kk]=sum;
		}
		// free_vector(b,1,m+1);
		// free_matrix(a,1,m+1,1,m+1);
		// free_ivector(indx,1,m+1);
	}
}
