package danfulea.math;

import danfulea.math.numerical.SpecFunc;
import danfulea.math.numerical.Stats;

/**
 * Statistics utilities. <br>
 * 
 * @author Dan Fulea, 02 JUN 2011
 * 
 */
public class StatsUtil {

	private static final int MAXITER = 10000;
	public static boolean failB = false;
	public static double confidenceLevel = 0.95;

	// true if different statistically, false otherwise!
	public static boolean differentB = false;

	/**
	 * Evaluate degrees of freedom based on given mean and standard deviation of
	 * mean. Usually, stdevOfMean = stdev(x)/SQRT(N);N=number of measurements on
	 * variable X and stdev(x) is standard deviation (also known as individual
	 * stdev)!
	 * 
	 * @param stdevOfMean standard deviation of mean
	 * @param mean the mean
	 * @return degrees of freedom
	 */
	public static double evaluateDegreesOfFreedom(double stdevOfMean,
			double mean) {
		failB = false;
		double f = mean * mean;
		if (stdevOfMean != 0.0) {
			f = f / (2.0 * stdevOfMean * stdevOfMean);
		} else {
			failB = true;
		}
		return f;
	}

	/**
	 * Evaluate deviation of a deviation based on given degrees of freedom and
	 * standard deviation (Gauss or Poisson or stdevOfMean). It is used in
	 * physics, in alpha, beta analysis.
	 * 
	 * @param stdev standard deviation
	 * @param f
	 *             degrees of freedom
	 * @return deviation of a deviation
	 */
	public static double evaluateDeviationOfDeviation(double stdev, double f) {
		failB = false;
		double aa = 0.0;
		if (f >= 1.0) {
			aa = stdev / Math.sqrt(2.0 * f);
		} else {
			failB = true;
		}
		return aa;
	}
	
	

	/**
	 * Compute if two dataset (suposed they have the same pooled variance) have
	 * different means! Degrees of freedom are evaluated inside this method.
	 * 
	 * @param ave1
	 *            , average means of dataset 1
	 * @param ave2
	 *            , average means of dataset 2
	 * @param sdevOfMean1
	 *            , standard deviation of mean 1
	 * @param sdevOfMean2
	 *            , standard deviation of mean 2
	 * @return true if the input means are different, false otherwise!
	 */
	public static boolean ttest(double ave1, double ave2, double sdevOfMean1,
			double sdevOfMean2) {
		failB = false;
		double prob_ttest = 0.0;
		double t_ttest = 0.0;
		double svar = 0.0;
		double df = 0.0;

		double f1 = evaluateDegreesOfFreedom(sdevOfMean1, ave1);
		if (failB) {
			return false;
		}
		double f2 = evaluateDegreesOfFreedom(sdevOfMean2, ave2);
		if (failB) {
			return false;
		}

		// evaluate variance: stdevOfMean=stdev/SQRT(N).
		double var1 = sdevOfMean1 * sdevOfMean1 * (f1 + 1);
		double var2 = sdevOfMean2 * sdevOfMean2 * (f2 + 1);

		// POOLED VARIANCE 1 tailed TTEST for mean comparison!!!!!
		df = f1 + f2;// Degrees of freedom.
		if (df <= 0) {
			failB = true;
			return false;
		}

		svar = (f1 * var1 + f2 * var2) / df;// Pooled variance.
		if (svar <= 0) {
			failB = true;
			return false;
		}

		t_ttest = Math.abs(ave1 - ave2)
				/ Math.sqrt(svar * (1.0 / (f1 + 1) + 1.0 / (f2 + 1)));

		prob_ttest = SpecFunc.betai(0.5 * df, 0.5, df
				/ (df + (t_ttest) * (t_ttest))); // See equation (6.4.9).

		prob_ttest = 1.0 - prob_ttest / 2.0;// 1-alpha/2.0 ONE TAILED!!!
		//System.out.println("f1= " +f1+" f2= "+f2);
		if (prob_ttest > confidenceLevel) {
			differentB = true;
		} else {
			differentB = false;
		}

		return differentB;
	}
	
	/**
	 * Compute if two dataset (suposed they have the same pooled variance) have
	 * different means!
	 * 
	 * @param ave1
	 *            , average means of dataset 1
	 * @param ave2
	 *            , average means of dataset 2
	 * @param stdevOfMean1
	 *            , standard deviation of mean of dataset 1
	 * @param stdevOfMean2
	 *            , standard deviation of mean of dataset 2
	 * @param f1
	 *            , degrees of freedom related to variance of dataset 1
	 * @param f2
	 *            , degrees of freedom related to variance of dataset 2
	 * @return true if the input means are different, false otherwise!
	 */
	public static boolean ttest_default_unc(double ave1, double ave2, double stdevOfMean1,
			double stdevOfMean2, double f1, double f2) {
		failB = false;
		double prob_ttest = 0.0;
		double t_ttest = 0.0;
		double svar = 0.0;
		double df = 0.0;

		// evaluate variance: stdevOfMean=stdev/SQRT(N).
		double var1 = stdevOfMean1 * stdevOfMean1 * (f1 + 1);
		double var2 = stdevOfMean2 * stdevOfMean2 * (f2 + 1);
		
		// POOLED VARIANCE 1 tailed TTEST for mean comparison!!!!!
		df = f1 + f2;// Degrees of freedom.
		if (df <= 0) {
			failB = true;
			return false;
		}

		svar = (f1 * var1 + f2 * var2) / df;// Pooled variance.
		if (svar <= 0) {
			failB = true;
			return false;
		}

		t_ttest = Math.abs(ave1 - ave2)
				/ Math.sqrt(svar * (1.0 / (f1 + 1) + 1.0 / (f2 + 1)));

		prob_ttest = SpecFunc.betai(0.5 * df, 0.5, df
				/ (df + (t_ttest) * (t_ttest))); // See equation (6.4.9).

		prob_ttest = 1.0 - prob_ttest / 2.0;// 1-alpha/2.0 ONE TAILED!!!
		// System.out.println(prob_ttest);
		if (prob_ttest > confidenceLevel) {
			differentB = true;
		} else {
			differentB = false;
		}

		return differentB;
	}

	/**
	 * Compute if two dataset (suposed they have the same pooled variance) have
	 * different means!
	 * 
	 * @param ave1
	 *            , average means of dataset 1
	 * @param ave2
	 *            , average means of dataset 2
	 * @param var1
	 *            , variance (standard deviation = SQRT(variance)) of dataset 1
	 * @param var2
	 *            , variance (standard deviation = SQRT(variance)) of dataset 2
	 * @param f1
	 *            , degrees of freedom related to variance of dataset 1
	 * @param f2
	 *            , degrees of freedom related to variance of dataset 2
	 * @return true if the input means are different, false otherwise!
	 */
	public static boolean ttest_default(double ave1, double ave2, double var1,
			double var2, double f1, double f2) {
		failB = false;
		double prob_ttest = 0.0;
		double t_ttest = 0.0;
		double svar = 0.0;
		double df = 0.0;

		// POOLED VARIANCE 1 tailed TTEST for mean comparison!!!!!
		df = f1 + f2;// Degrees of freedom.
		if (df <= 0) {
			failB = true;
			return false;
		}

		svar = (f1 * var1 + f2 * var2) / df;// Pooled variance.
		if (svar <= 0) {
			failB = true;
			return false;
		}

		t_ttest = Math.abs(ave1 - ave2)
				/ Math.sqrt(svar * (1.0 / (f1 + 1) + 1.0 / (f2 + 1)));

		prob_ttest = SpecFunc.betai(0.5 * df, 0.5, df
				/ (df + (t_ttest) * (t_ttest))); // See equation (6.4.9).

		prob_ttest = 1.0 - prob_ttest / 2.0;// 1-alpha/2.0 ONE TAILED!!!
		// System.out.println(prob_ttest);
		if (prob_ttest > confidenceLevel) {
			differentB = true;
		} else {
			differentB = false;
		}

		return differentB;
	}

	/**
	 * Evaluate if two deviations (usually experimental standard deviation and
	 * theoretical standard deviation, in other words Gauss and Poisson
	 * respectively)are different. It is used in physics, in alpha, beta
	 * analysis.
	 * 
	 * @param dev1
	 *            ,deviation 1 related to dataset (usually stdev Gauss or
	 *            experimental)
	 * @param dev2
	 *            ,deviation 2 related to dataset (usually stdev Poisson or
	 *            theoretical)
	 * @param f1
	 *            ,degrees of freedom related to variance 1 of dataset
	 * @param f2
	 *            ,degrees of freedom related to variance 2 of dataset
	 * @return true if the input means are different, false otherwise!
	 */
	public static boolean ttest_deviation(double dev1, double dev2, double f1,
			double f2) {
		failB = false;
		double prob_ttest = 0.0;
		double t_ttest = 0.0;
		double svar = 0.0;
		double df = 0.0;

		double var1 = evaluateDeviationOfDeviation(dev1, f1);//always get stdevOfMean!!
		var1 = var1 * var1* (f1 + 1);// variance!
		if (failB) {
			return false;
		}
		double var2 = evaluateDeviationOfDeviation(dev2, f2);
		var2 = var2 * var2* (f2 + 1);// variance!
		if (failB) {
			return false;
		}

		// POOLED VARIANCE 1 tailed TTEST for mean comparison!!!!!
		df = f1 + f2;// Degrees of freedom.
		if (df <= 0) {
			failB = true;
			return false;
		}

		svar = (f1 * var1 + f2 * var2) / df;// Pooled variance.
		if (svar <= 0) {
			failB = true;
			return false;
		}

		t_ttest = Math.abs(dev1 - dev2)
				/ Math.sqrt(svar * (1.0 / (f1 + 1) + 1.0 / (f2 + 1)));

		prob_ttest = SpecFunc.betai(0.5 * df, 0.5, df
				/ (df + (t_ttest) * (t_ttest))); // See equation (6.4.9).

		prob_ttest = 1.0 - prob_ttest / 2.0;// 1-alpha/2.0 ONE TAILED!!!
		// System.out.println(prob_ttest);
		if (prob_ttest > confidenceLevel) {
			differentB = true;
		} else {
			differentB = false;
		}

		return differentB;
	}

	/**
	 * Fisher test used in physics. All required default data are computed
	 * inside method.
	 * 
	 * @param ave1
	 *            , mean 1
	 * @param ave2
	 *            , mean 2
	 * @param sdevOfMean1
	 *            , standard deviation of mean 1
	 * @param sdevOfMean2
	 *            , standard deviation of mean 2
	 * @return true if variances are different, false otherwise!
	 */
	public static boolean ftest(double ave1, double ave2, double sdevOfMean1,
			double sdevOfMean2) {
		failB = false;
		double f_ftest = 0.0;

		double f1 = evaluateDegreesOfFreedom(sdevOfMean1, ave1);
		if (failB) {
			return false;
		}
		double f2 = evaluateDegreesOfFreedom(sdevOfMean2, ave2);
		if (failB) {
			return false;
		}

		// evaluate variance: stdevOfMean=stdev/SQRT(N).
		double var1 = sdevOfMean1 * sdevOfMean1 * (f1 + 1);
		double var2 = sdevOfMean2 * sdevOfMean2 * (f2 + 1);

		// now default:
		double df1 = f1;
		double df2 = f2;
		if (var1 > var2) {
			// Make F the ratio of the larger variance to the smaller one.
			if (var2 == 0) {
				failB = true;
				return false;
			}
			f_ftest = var1 / var2;
			df1 = f1;
			df2 = f2;
		} else {
			if (var1 == 0) {
				failB = true;
				return false;
			}
			f_ftest = var2 / var1;
			df1 = f2;
			df2 = f1;
		}

		double prob_ftest = 2.0 * SpecFunc.betai(0.5 * df2, 0.5 * df1, df2
				/ (df2 + df1 * (f_ftest)));
		if (prob_ftest > 1.0)
			prob_ftest = 2.0 - prob_ftest;

		prob_ftest = 1.0 - prob_ftest / 2.0;
		if (prob_ftest > confidenceLevel) {
			differentB = true;
		} else {
			differentB = false;
		}

		return differentB;
	}

	/**
	 * Fisher test. It is used in physics, in alpha, beta analysis.
	 * 
	 * @param var1
	 *            , variance 1
	 * @param var2
	 *            , variance 2
	 * @param f1
	 *            , degrees of freedom related to variance 1
	 * @param f2
	 *            , degrees of freedom related to variance 2
	 * @return true if the input variances are different, false otherwise!
	 */
	public static boolean ftest_default(double var1, double var2, double f1,
			double f2) {
		failB = false;
		double f_ftest = 0.0;
		double df1 = f1;
		double df2 = f2;
		if (var1 > var2) {
			// Make F the ratio of the larger variance to the smaller one.
			if (var2 == 0) {
				failB = true;
				return false;
			}
			f_ftest = var1 / var2;
			df1 = f1;
			df2 = f2;
		} else {
			if (var1 == 0) {
				failB = true;
				return false;
			}
			f_ftest = var2 / var1;
			df1 = f2;
			df2 = f1;
		}

		double prob_ftest = 2.0 * SpecFunc.betai(0.5 * df2, 0.5 * df1, df2
				/ (df2 + df1 * (f_ftest)));
		if (prob_ftest > 1.0)
			prob_ftest = 2.0 - prob_ftest;

		prob_ftest = 1.0 - prob_ftest / 2.0;
		if (prob_ftest > confidenceLevel) {
			differentB = true;
		} else {
			differentB = false;
		}

		return differentB;
	}

	/**
	 * Evaluate if two deviations (usually experimental standard deviation and
	 * theoretical standard deviation of mean in other words Gauss and Poisson
	 * respectively) are different. It is used in physics, in alpha, beta
	 * analysis.
	 * 
	 * @param sdevOfMean1
	 *            ,standard deviation of mean 1 related to dataset (usually
	 *            Gauss or experimental)
	 * @param sdevOfMean2
	 *            ,standard deviation of mean 2 related to dataset (usually
	 *            Poisson or theoretical)
	 * @param f1
	 *            , degrees of freedom related to deviation 1
	 * @param f2
	 *            , degrees of freedom related to deviation 2
	 * @return true if the input deviations are different, false otherwise!
	 */
	public static boolean ftest_deviation(double sdevOfMean1,
			double sdevOfMean2, double f1, double f2) {
		failB = false;
		double f_ftest = 0.0;
		// evaluate variance: stdevOfMean=stdev/SQRT(N).
		double var1 = sdevOfMean1 * sdevOfMean1 * (f1 + 1);
		double var2 = sdevOfMean2 * sdevOfMean2 * (f2 + 1);

		// now default:
		double df1 = f1;
		double df2 = f2;
		if (var1 > var2) {
			// Make F the ratio of the larger variance to the smaller one.
			if (var2 == 0) {
				failB = true;
				return false;
			}
			f_ftest = var1 / var2;
			df1 = f1;
			df2 = f2;
		} else {
			if (var1 == 0) {
				failB = true;
				return false;
			}
			f_ftest = var2 / var1;
			df1 = f2;
			df2 = f1;
		}

		double prob_ftest = 2.0 * SpecFunc.betai(0.5 * df2, 0.5 * df1, df2
				/ (df2 + df1 * (f_ftest)));
		if (prob_ftest > 1.0)
			prob_ftest = 2.0 - prob_ftest;

		prob_ftest = 1.0 - prob_ftest / 2.0;
		if (prob_ftest > confidenceLevel) {
			differentB = true;
		} else {
			differentB = false;
		}

		return differentB;
	}

	// ------------------------
	/**
	 * Get the coverage factor for final result expressed as mean +/- t multiply
	 * by the standard deviation of mean. Usually, degrees of freedom are
	 * computed as n-1 where n is the number of measurements on x variable.
	 * Equation is: result(or simply X) = MEAN(x) +/- t * STDEV(x)/SQRT(n).
	 * STDEV(x)/SQRT(n) is also known as standard deviation of mean!
	 * 
	 * @param df
	 *            , degrees of freedom
	 * @return t, student factor
	 */
	public static double getStudentFactor(double df) {

		failB = false;

		int iter = 0;

		double alpha = 0.0;
		double probab = 1.0 - alpha;// two-tailed!!

		double tlo = 0.0;
		double thi = 1000;
		double t = thi;

		while (true) {
			iter = iter + 1;// begin iteration!
			probab = 1.0 - alpha;// two-tailed!!
			if (probab > confidenceLevel) {
				thi = t;
			} else {
				tlo = t;
			}

			if (probab < confidenceLevel * 1.00005
					&& probab > confidenceLevel * 0.99995)
				break;

			t = (tlo + thi) / 2.0;// 1st guess of t factor
			alpha = SpecFunc.betai(0.5 * df, 0.5, df / (df + (t) * (t)));
			// -----fail safe
			if (iter > MAXITER) {
				failB = true;
				break;
			}
		}
		// test:
		// double
		// tint=1.812;//1.372;//1.812;//2.634;//2.228;//2.228;//0.700;//3.16927;//1.812;
		// alpha=SpecFunc.betai(0.5*df,0.5,df/(df+(tint)*(tint)));

		// display
		// System.out.println("alpha= " + alpha + " twotailed p%= "
		// + (1.0 - alpha) + " onetailed %= " + (1.0 - alpha / 2.0)
		// + " t= " + t + " iter " + iter + " fail? " + failB);// two
		// tailed!!
		// here, we have two-tailed ttest which is 1.0-alpha!
		// when comparing two datasets, we will use
		// one-tailed ttest which is 1.0-alpha/2.0.

		return t;
	}

	/**
	 
	 * Return degrees of freedom related to the coverage factor 
	 * (student t factor, two-tailed) involving in 
	 * extended uncertainty value!!
	 * 
	 * @param tst
	 *             student factor
	 * @return f - degrees of freedom 
	 */
	public static double getDegreesOfFreedomFromStudentFactor(double tst) {

		failB = false;

		int iter = 0;

		double alpha = 0.0;
		double probab = 1.0 - alpha;// two-tailed!!

		double flo = 0.0;
		double fhi = 100000.0;
		double f = fhi;

		while (true) {
			iter = iter + 1;// begin iteration!
			probab = 1.0 - alpha;// two-tailed!!
			if (probab > confidenceLevel) {
				fhi = f;
			} else {
				flo = f;
			}

			if (probab < confidenceLevel * 1.00005
					&& probab > confidenceLevel * 0.99995)
				break;

			f = (flo + fhi) / 2.0;// 1st guess of t factor
			alpha = SpecFunc.betai(0.5 * f, 0.5, f / (f + (tst) * (tst)));
			// -----fail safe
			if (iter > MAXITER) {
				failB = true;
				break;
			}
		}
		
		return f;
	}
	// ===========================================================================================
	// global variance of n stdev! It is used in physics, in alpha, beta
	// measurements.
	/**
	 * Computes the global variance of n stdev! It is used in physics, in alpha, beta measurements.
	 * @param ab the standard deviation array
	 * @param f the associated array of degrees of freedom
	 * @return the global variance.
	 */
	public static double getGlobalVariance(double[] ab, double[] f) {
		failB = false;
		double gv = 0.0;

		if (ab.length == f.length) {
			double numarator = 0.0;
			double numitor = 0.0;
			if (ab.length >= 1) {
				for (int i = 0; i < ab.length; i++) {
					numarator = numarator + ab[i] * ab[i] * f[i];
					numitor = numitor + f[i];
				}
				if (numitor != 0.0) {
					gv = numarator / numitor;
				} else {
					failB = true;
				}
			}
		}

		return gv;
	}

	// Welch - Satterwaite formula
	// If f=f(x1...xn) and if s1,...sn are the related uncertainties for
	// variables x1,...xn then:
	// Scompus=SUMi=1,n(df/dxi)^2*si^2+2*SUMi=1,n-1SUMj=i+1,N((df/dxi)*(df/dxj)*si*sj*rij)
	// rij=correlation coef. It is used in physics, in alpha, beta measurements.
	// abcompus is the final uncertainty based on error propagation formula
	// described above,
	// ab array are individual uncertainties of variables which are contained in
	// the final formula.
	// e.g. A=Net/eff where Net and eff present individual uncertainties.
	// f array are degrees of freedom associated with ab array uncertainties.
	/**
	 * Welch - Satterwaite formula for computing the degrees of freedom.
	 * If f=f(x1...xn) and if s1,...sn are the related uncertainties for variables x1,...xn then: 
	 * abcompus^2=SUMi=1,n(df/dxi)^2*si^2+2*SUMi=1,n-1SUMj=i+1,N((df/dxi)*(df/dxj)*si*sj*rij); 
	 * rij=correlation coefficient. It is used in physics, in alpha, beta measurements.
	 * @param abcompus the final uncertainty based on error propagation formula described above
	 * @param ab array are individual uncertainties of variables which are contained in the final formula. 
	 * e.g. A=Net/eff where Net and eff present individual uncertainties.
	 * @param f array are degrees of freedom associated with ab array uncertainties.
	 * @return the final overall degrees of freedom associated with final uncertainty abcompus.
	 */
	public static double getEffectiveDegreesOfFreedom(double abcompus,
			double[] ab, double[] f) {
		failB = false;
		boolean b = false;
		double num = 0.0;
		double gv = Math.pow(abcompus, 4);
		if (ab.length == f.length && ab.length >= 1) {
			for (int i = 0; i < f.length; i++)
				if (f[i] == 0.0) {
					b = false;
					failB = true;
					break;
				} else {
					num = num + Math.pow(ab[i], 4) / f[i];
					b = true;
				}

			if (b) {
				if (num != 0.0) {
					gv = gv / num;
				} else {
					failB = true;
				}
			}
		}

		return gv;
	}

	// Correlation coeff for 2 datasets.
	// It is used in physics, in alpha, beta measurements.
	/**
	 * Computes the correlation coefficient for 2 datasets.
	 * @param x the 1st dataset
	 * @param y the 2nd dataset
	 * @return the correlation coefficient
	 */
	public static double getDataCorrelation(double[] x, double[] y) {
		failB = false;
		double cor = 0.0;
		if (x.length >= 2 && y.length >= 2) {
			Stats.avevar(x, x.length);
			double xm = Stats.ave_avevar;
			double sx = Stats.var_avevar;
			sx = Math.sqrt(sx) / Math.sqrt(x.length);

			Stats.avevar(y, y.length);
			double ym = Stats.ave_avevar;
			double sy = Stats.var_avevar;
			sy = Math.sqrt(sy) / Math.sqrt(y.length);

			int n = Math.min(x.length, y.length);
			double num = sx
					* sy
					* Math.sqrt(x.length * y.length * (x.length - 1)
							* (y.length - 1));
			double numarator = 0.0;
			for (int i = 0; i < n; i++) {
				numarator = numarator + (x[i] - xm) * (y[i] - ym);
			}
			if (num != 0.0) {
				cor = numarator / num;
			} else {
				failB = true;
			}
		}

		return cor;
	}

	// Weughted mean based on default weight w.
	// It is used in physics, in alpha, beta measurements,
	// when we have several results of a variable (activity) each result having
	// its own uncertainty=>Overall final result is computed by weighted mean!!
	/**
	 * Weighted mean based on default weight w. It is used in physics, in alpha, beta measurements, 
	 * when we have several results of a variable (activity) each result having 
	 * its own uncertainty. The overall final result is computed by weighted mean!!
	 * @param x the array of means
	 * @param sx the array of their uncertainties
	 * @return the weighted mean.
	 */
	public static double getWeightedMean(double[] x, double[] sx) {
		failB = false;
		boolean b = false;
		double m = 0.0;

		if (x.length >= 1 && sx.length >= 1 && x.length == sx.length) {
			double[] w = new double[sx.length];
			double numarator = 0.0;
			double numitor = 0.0;
			for (int i = 0; i < sx.length; i++) {
				if (sx[i] != 0) {
					w[i] = 1.0 / (sx[i] * sx[i]);// default
					numarator = numarator + x[i] * w[i];
					numitor = numitor + w[i];
					b = true;

				} else {
					b = false;
					break;
				}
			}

			if (b) {
				m = numarator / numitor;
			} else {
				failB = true;
			}
		}

		return m;
	}

	// Weighted uncertainty based on default weight w.
	// It is used in physics, in alpha, beta measurements.
	/**
	 * Weighted uncertainty based on default weight w. It is used in physics, in alpha, beta measurements.
	 * @param sx the array of uncertainties
	 * @return the weighted uncertainty
	 */
	public static double getWeightedDeviation(double[] sx) {
		failB = false;
		boolean b = false;
		double m = 0.0;

		if (sx.length >= 1) {
			double[] w = new double[sx.length];
			double numarator = 1.0;
			double numitor = 0.0;
			for (int i = 0; i < sx.length; i++) {
				if (sx[i] != 0) {
					w[i] = 1.0 / (sx[i] * sx[i]);// default
					numitor = numitor + w[i];
					b = true;

				} else {
					b = false;
					break;
				}
			}

			if (b) {
				m = Math.sqrt(numarator / numitor);
			} else {
				failB = true;
			}
		}

		return m;
	}

	// ===========================================================================================
	// ===========================================================================================
	/**
	 * Get the Student factor for mean comparison of two datasets.
	 * 
	 * @param df
	 *            , degrees of freedom
	 * @return t, student factor
	 */
	public static double getStudentFactorForMeanComparison(double df) {

		failB = false;

		int iter = 0;

		double alpha = 0.0;
		double probab = 1.0 - alpha / 2.0;// one-tailed!!

		double tlo = 0.0;
		double thi = 1000;
		double t = thi;

		while (true) {
			iter = iter + 1;// begin iteration!
			probab = 1.0 - alpha / 2.0;// one-tailed!!
			if (probab > confidenceLevel) {
				thi = t;
			} else {
				tlo = t;
			}

			if (probab < confidenceLevel * 1.00005
					&& probab > confidenceLevel * 0.99995)
				break;

			t = (tlo + thi) / 2.0;// 1st guess of t factor
			alpha = SpecFunc.betai(0.5 * df, 0.5, df / (df + (t) * (t)));
			// -----fail safe
			if (iter > MAXITER) {
				failB = true;
				break;
			}
		}

		return t;
	}

	/**
	 * Get the Fisher factor according to the degrees of freedom and the
	 * required confidence level.
	 * 
	 * @param df1
	 *            , degrees of freedom for dataset 1
	 * @param df2
	 *            , degrees of freedom for dataset 2
	 * @return t, Fisher factor
	 */
	public static double getFisherFactor(double df1, double df2) {

		failB = false;

		int iter = 0;

		double alpha = 0.0;
		double probab = 1.0 - alpha;

		double tlo = 0.0;
		double thi = 1000;
		double t = thi;

		while (true) {
			iter = iter + 1;// begin iteration!
			probab = 1.0 - alpha;
			if (probab > confidenceLevel) {
				thi = t;
			} else {
				tlo = t;
			}
			// t=0-1=>probab decrease from 1 to 0.5
			// else if t>1 => probab increase from 0.5 to 1
			// but always t>1!!!

			if (probab < confidenceLevel * 1.00005
					&& probab > confidenceLevel * 0.99995)
				break;

			t = (tlo + thi) / 2.0;// 1st guess of f factor
			alpha = 2.0 * SpecFunc.betai(0.5 * df2, 0.5 * df1, df2
					/ (df2 + df1 * (t)));
			if (alpha > 1.0)
				alpha = 2.0 - alpha;
			alpha = alpha / 2.0;
			// -----fail safe
			if (iter > MAXITER) {
				failB = true;
				break;
			}
		}

		return t;
	}

}
