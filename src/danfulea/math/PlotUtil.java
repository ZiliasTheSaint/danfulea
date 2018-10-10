package danfulea.math;

import danfulea.math.numerical.Interpolator;
import danfulea.phys.HvlUtil;

/**
 * Utility class for plotting HVL function. 
 * @author Dan Fulea, 2006
 *
 */
@Deprecated
public class PlotUtil {
	// metoda spline directa
	// return perechile XY pentru trasarea graficului
	public static double[][] getSplineXY(double[] x, double[] y) {
		int ni = x.length;

		double[] a = new double[ni];
		double[] b = new double[ni];
		double[] c = new double[ni];
		double[] d = new double[ni];

		double[] xi = new double[x.length];
		double[] yi = new double[y.length];
		for (int i = 0; i < x.length; i++) {
			xi[i] = x[i];
			yi[i] = y[i];
		}
		Sort.qSort2(xi, yi);

		// /////////evaluare pas--pentru grafic
		double[] pas = new double[ni - 1];
		for (int i = 0; i <= ni - 2; i++)
			pas[i] = (xi[i + 1] - xi[i]) / 5;
		// pasul cel mai mic:
		double pasul = Sort.findValue(pas, ni - 1);
		// necesita sortare crescatoare->de cate ori incape acest pas
		double dbl = (xi[xi.length - 1] - xi[0]) / pasul;
		long j = Math.round(dbl);
		int n = 0;
		// partea intreaga a unui numar zecimal!!!
		if (j > dbl)
			n = (int) j - 1;
		else
			n = (int) j;// exemplu=10<->10*pasul acopera tot intervalul

		double[] xs = new double[n + 1];
		double[] ys = new double[n + 1];

		double[][] result = new double[xs.length][2];// 0-x;1-y

		xs[0] = xi[0];
		for (int i = 1; i <= n; i++)
			xs[i] = xi[0] + i * pasul;
		// gata evaluare tablouri pentru grafic
		// spline cu derivatele la capete deduse prin metoda Lagrange de
		// interpolare
		double dm = (yi[1] - yi[0]) / (xi[1] - xi[0])
				- ((yi[2] - yi[1]) / (xi[2] - xi[1])) + (yi[2] - yi[0])
				/ (xi[2] - xi[0]);
		double dp = -((yi[ni - 2] - yi[ni - 3]) / (xi[ni - 2] - xi[ni - 3]))
				+ (yi[ni - 1] - yi[ni - 2]) / (xi[ni - 1] - xi[ni - 2])
				+ (yi[ni - 1] - yi[ni - 3]) / (xi[ni - 1] - xi[ni - 3]);

		Interpolator.spline(xi, yi, ni, xs, ys, n + 1, a, b, c, d, dm, dp);

		for (int i = 0; i < xs.length; i++) {
			result[i][0] = xs[i];
			result[i][1] = ys[i];
		}

		return result;
	}

	public static double[][] getPolynomialXY(double[] x, double[] y, int n) {
		int ni = x.length;
		double[][] a = new double[n + 1][1];
		// sortez sirurile dupa x---just in case
		double[] xi = new double[x.length];
		double[] yi = new double[y.length];
		for (int i = 0; i < x.length; i++) {
			xi[i] = x[i];
			yi[i] = y[i];
		}
		Sort.qSort2(xi, yi);
		// ------------------------------------
		Interpolator.polynomial(xi, yi, n, a);
		double[] d = HvlUtil.convertMatrix(a, 0);

		// /////////evaluare pas--pentru grafic
		double[] pas = new double[ni - 1];
		for (int i = 0; i <= ni - 2; i++)
			pas[i] = (xi[i + 1] - xi[i]) / 5;
		// pasul cel mai mic:
		double pasul = Sort.findValue(pas, ni - 1);

		// necesita sortare crescatoare->de cate ori incape acest pas
		double dbl = (xi[xi.length - 1] - xi[0]) / pasul;
		long j = Math.round(dbl);
		int nn = 0;
		// partea intreaga a unui numar zecimal!!!
		if (j > dbl)
			nn = (int) j - 1;
		else
			nn = (int) j;// exemplu=10<->10*pasul acopera tot intervalul

		double[] xs = new double[nn + 1];
		double[] ys = new double[nn + 1];

		double[][] result = new double[xs.length][2];// 0-x;1-y

		xs[0] = xi[0];
		ys[0] = yi[0];
		result[0][0] = xs[0];
		result[0][1] = ys[0];

		for (int i = 1; i < xs.length; i++) {
			xs[i] = xi[0] + i * pasul;
			for (int k = 0; k < d.length; k++) {
				ys[i] = ys[i] + d[k] * Math.pow(xs[i], k);
			}
			result[i][0] = xs[i];
			result[i][1] = ys[i];
		}
		// gata evaluare tablouri pentru grafic

		return result;
	}

	// metoda spline revert
	// return perechile XY pentru trasarea graficului
	public static double[][] getSplineRevertXY(double[] x, double[] y) {
		int ni = x.length;

		double[] a = new double[ni];
		double[] b = new double[ni];
		double[] c = new double[ni];
		double[] d = new double[ni];

		double[] xi = new double[x.length];
		double[] yi = new double[y.length];
		for (int i = 0; i < x.length; i++) {
			xi[i] = x[i];
			yi[i] = y[i];
		}
		Sort.qSort2(yi, xi);

		// /////////evaluare pas--pentru grafic
		double[] pas = new double[ni - 1];
		for (int i = 0; i <= ni - 2; i++)
			pas[i] = (yi[i + 1] - yi[i]) / 5;
		// pasul cel mai mic:
		double pasul = Sort.findValue(pas, ni - 1);
		// necesita sortare crescatoare->de cate ori incape acest pas
		double dbl = (yi[yi.length - 1] - yi[0]) / pasul;
		long j = Math.round(dbl);
		int n = 0;
		// partea intreaga a unui numar zecimal!!!
		if (j > dbl)
			n = (int) j - 1;
		else
			n = (int) j;// exemplu=10<->10*pasul acopera tot intervalul

		double[] xs = new double[n + 1];
		double[] ys = new double[n + 1];

		double[][] result = new double[xs.length][2];// 0-x;1-y

		xs[0] = yi[0];
		for (int i = 1; i <= n; i++)
			xs[i] = yi[0] + i * pasul;
		// gata evaluare tablouri pentru grafic
		// spline cu derivatele la capete deduse prin metoda Lagrange de
		// interpolare
		double dm = (xi[1] - xi[0]) / (yi[1] - yi[0])
				- ((xi[2] - xi[1]) / (yi[2] - yi[1])) + (xi[2] - xi[0])
				/ (yi[2] - yi[0]);
		double dp = -((xi[ni - 2] - xi[ni - 3]) / (yi[ni - 2] - yi[ni - 3]))
				+ (xi[ni - 1] - xi[ni - 2]) / (yi[ni - 1] - yi[ni - 2])
				+ (xi[ni - 1] - xi[ni - 3]) / (yi[ni - 1] - yi[ni - 3]);

		Interpolator.spline(yi, xi, ni, xs, ys, n + 1, a, b, c, d, dm, dp);

		for (int i = 0; i < xs.length; i++) {
			result[i][0] = xs[i];
			result[i][1] = ys[i];
		}

		return result;
	}

	public static double[][] getPolynomialRevertXY(double[] x, double[] y, int n) {
		int ni = x.length;
		double[][] a = new double[n + 1][1];
		// sortez sirurile dupa x---just in case
		double[] xi = new double[x.length];
		double[] yi = new double[y.length];
		for (int i = 0; i < x.length; i++) {
			xi[i] = x[i];
			yi[i] = y[i];
		}
		Sort.qSort2(yi, xi);
		// ------------------------------------
		Interpolator.polynomial(yi, xi, n, a);
		double[] d = HvlUtil.convertMatrix(a, 0);
		// /////////evaluare pas--pentru grafic
		double[] pas = new double[ni - 1];
		for (int i = 0; i <= ni - 2; i++)
			pas[i] = (yi[i + 1] - yi[i]) / 5;
		// pasul cel mai mic:
		double pasul = Sort.findValue(pas, ni - 1);

		// necesita sortare crescatoare->de cate ori incape acest pas
		double dbl = (yi[yi.length - 1] - yi[0]) / pasul;
		long j = Math.round(dbl);
		int nn = 0;
		// partea intreaga a unui numar zecimal!!!
		if (j > dbl)
			nn = (int) j - 1;
		else
			nn = (int) j;// exemplu=10<->10*pasul acopera tot intervalul

		double[] xs = new double[nn + 1];
		double[] ys = new double[nn + 1];

		double[][] result = new double[xs.length][2];// 0-x;1-y

		xs[0] = yi[0];
		ys[0] = xi[0];
		result[0][0] = xs[0];
		result[0][1] = ys[0];

		for (int i = 1; i < xs.length; i++) {
			xs[i] = yi[0] + i * pasul;
			for (int k = 0; k < d.length; k++) {
				ys[i] = ys[i] + d[k] * Math.pow(xs[i], k);
			}
			result[i][0] = xs[i];
			result[i][1] = ys[i];
		}
		// gata evaluare tablouri pentru grafic

		return result;
	}

	public static double[][] getSplineDev(double[] x, double[] y) {
		int ni = x.length;

		double[] a = new double[ni];
		double[] b = new double[ni];
		double[] c = new double[ni];
		double[] d = new double[ni];

		double[] xi = new double[x.length];
		double[] yi = new double[y.length];
		for (int i = 0; i < x.length; i++) {
			xi[i] = x[i];
			yi[i] = y[i];
		}
		Sort.qSort2(xi, yi);

		double[] xs = new double[ni];
		double[] ys = new double[ni];
		for (int i = 0; i < x.length; i++)
			xs[i] = xi[i];

		double dm = (yi[1] - yi[0]) / (xi[1] - xi[0])
				- ((yi[2] - yi[1]) / (xi[2] - xi[1])) + (yi[2] - yi[0])
				/ (xi[2] - xi[0]);
		double dp = -((yi[ni - 2] - yi[ni - 3]) / (xi[ni - 2] - xi[ni - 3]))
				+ (yi[ni - 1] - yi[ni - 2]) / (xi[ni - 1] - xi[ni - 2])
				+ (yi[ni - 1] - yi[ni - 3]) / (xi[ni - 1] - xi[ni - 3]);

		Interpolator.spline(xi, yi, ni, xs, ys, ni, a, b, c, d, dm, dp);

		double[][] result = new double[xs.length][3];// 0-yi;1-yinterp;2-deviation=yinterp-yi

		for (int i = 0; i < xs.length; i++) {
			result[i][0] = yi[i];
			result[i][1] = ys[i];
			if (yi[i] != 0)
				result[i][2] = (ys[i] - yi[i]) * 100 / yi[i];// %!!!
			else
				result[i][2] = (ys[i] - yi[i]);
		}

		return result;
	}

	public static double[][] getPolynomialDev(double[] x, double[] y, int n) {
		int ni = x.length;
		double[][] a = new double[n + 1][1];
		// sortez sirurile dupa x---just in case
		double[] xi = new double[x.length];
		double[] yi = new double[y.length];
		for (int i = 0; i < x.length; i++) {
			xi[i] = x[i];
			yi[i] = y[i];
		}
		Sort.qSort2(xi, yi);
		// ------------------------------------
		Interpolator.polynomial(xi, yi, n, a);
		double[] d = HvlUtil.convertMatrix(a, 0);

		double[] xs = new double[ni];
		double[] ys = new double[ni];
		for (int i = 0; i < x.length; i++)
			xs[i] = xi[i];

		double[][] result = new double[xs.length][3];// 0-yi;1-yinterp;2-dev

		for (int i = 0; i < xs.length; i++) {
			for (int k = 0; k < d.length; k++) {
				ys[i] = ys[i] + d[k] * Math.pow(xs[i], k);// automat e
															// initializat cu
															// zero!!
			}
			result[i][0] = yi[i];
			result[i][1] = ys[i];
			if (yi[i] != 0)
				result[i][2] = (ys[i] - yi[i]) * 100 / yi[i];// %!!!
			else
				result[i][2] = (ys[i] - yi[i]);
		}
		// gata evaluare tablouri pentru grafic

		return result;
	}

	public static double[][] getSplineRevertDev(double[] x, double[] y) {
		int ni = x.length;

		double[] a = new double[ni];
		double[] b = new double[ni];
		double[] c = new double[ni];
		double[] d = new double[ni];

		double[] xi = new double[x.length];
		double[] yi = new double[y.length];
		for (int i = 0; i < x.length; i++) {
			xi[i] = x[i];
			yi[i] = y[i];
		}
		Sort.qSort2(yi, xi);

		double[] xs = new double[ni];
		double[] ys = new double[ni];
		for (int i = 0; i < x.length; i++)
			xs[i] = yi[i];

		double dm = (xi[1] - xi[0]) / (yi[1] - yi[0])
				- ((xi[2] - xi[1]) / (yi[2] - yi[1])) + (xi[2] - xi[0])
				/ (yi[2] - yi[0]);
		double dp = -((xi[ni - 2] - xi[ni - 3]) / (yi[ni - 2] - yi[ni - 3]))
				+ (xi[ni - 1] - xi[ni - 2]) / (yi[ni - 1] - yi[ni - 2])
				+ (xi[ni - 1] - xi[ni - 3]) / (yi[ni - 1] - yi[ni - 3]);

		Interpolator.spline(yi, xi, ni, xs, ys, ni, a, b, c, d, dm, dp);

		double[][] result = new double[xs.length][3];// 0-xi;1-xinterp;2-deviation=xinterp-xi

		for (int i = 0; i < xs.length; i++) {
			result[i][0] = xi[i];
			result[i][1] = ys[i];
			if (xi[i] != 0)
				result[i][2] = (ys[i] - xi[i]) * 100 / xi[i];
			else
				result[i][2] = (ys[i] - xi[i]);
		}

		return result;
	}

	public static double[][] getPolynomialRevertDev(double[] x, double[] y,
			int n) {
		int ni = x.length;
		double[][] a = new double[n + 1][1];
		// sortez sirurile dupa x---just in case
		double[] xi = new double[x.length];
		double[] yi = new double[y.length];
		for (int i = 0; i < x.length; i++) {
			xi[i] = x[i];
			yi[i] = y[i];
		}
		Sort.qSort2(yi, xi);
		// ------------------------------------
		Interpolator.polynomial(yi, xi, n, a);
		double[] d = HvlUtil.convertMatrix(a, 0);

		double[] xs = new double[ni];
		double[] ys = new double[ni];
		for (int i = 0; i < x.length; i++)
			xs[i] = yi[i];

		double[][] result = new double[xs.length][3];

		for (int i = 0; i < xs.length; i++) {
			for (int k = 0; k < d.length; k++) {
				ys[i] = ys[i] + d[k] * Math.pow(xs[i], k);
			}
			result[i][0] = xi[i];
			result[i][1] = ys[i];
			if (xi[i] != 0)
				result[i][2] = (ys[i] - xi[i]) * 100 / xi[i];
			else
				result[i][2] = (ys[i] - xi[i]);
		}

		return result;
	}

	public static double[] getPolyCoeff(double[] x, double[] y, int n) {
		//int ni = x.length;
		double[][] a = new double[n + 1][1];
		// sortez sirurile dupa x---just in case
		double[] xi = new double[x.length];
		double[] yi = new double[y.length];
		for (int i = 0; i < x.length; i++) {
			xi[i] = x[i];
			yi[i] = y[i];
		}
		Sort.qSort2(xi, yi);
		// ------------------------------------
		Interpolator.polynomial(xi, yi, n, a);
		double[] d = HvlUtil.convertMatrix(a, 0);

		return d;
	}

	public static double[][] getSplineCoeff(double[] x, double[] y) {
		int ni = x.length;

		double[][] result = new double[ni][4];

		double[] a = new double[ni];// 0
		double[] b = new double[ni];// 1
		double[] c = new double[ni];// 2
		double[] d = new double[ni];// 3

		double[] xi = new double[x.length];
		double[] yi = new double[y.length];
		for (int i = 0; i < x.length; i++) {
			xi[i] = x[i];
			yi[i] = y[i];
		}
		Sort.qSort2(xi, yi);

		double[] xs = new double[ni];
		double[] ys = new double[ni];
		for (int i = 0; i < x.length; i++)
			xs[i] = xi[i];

		double dm = (yi[1] - yi[0]) / (xi[1] - xi[0])
				- ((yi[2] - yi[1]) / (xi[2] - xi[1])) + (yi[2] - yi[0])
				/ (xi[2] - xi[0]);
		double dp = -((yi[ni - 2] - yi[ni - 3]) / (xi[ni - 2] - xi[ni - 3]))
				+ (yi[ni - 1] - yi[ni - 2]) / (xi[ni - 1] - xi[ni - 2])
				+ (yi[ni - 1] - yi[ni - 3]) / (xi[ni - 1] - xi[ni - 3]);

		Interpolator.spline(xi, yi, ni, xs, ys, ni, a, b, c, d, dm, dp);

		for (int i = 0; i < ni; i++) {
			result[i][0] = a[i];
			result[i][1] = b[i];
			result[i][2] = c[i];
			result[i][3] = d[i];
		}

		return result;
	}
}
