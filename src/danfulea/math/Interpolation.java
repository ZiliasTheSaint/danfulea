package danfulea.math;

/**
 * 
 * Handles widely used interpolation methods
 * @author Dan Fulea, 15 APR. 2011
 *
 */
public class Interpolation {

	/**
	 * Linear interpolation
	 * @param x1 first point x-value
	 * @param y1 first point y-value
	 * @param x2 second point x-value
	 * @param y2 second point y-value
	 * @param x desire point x-value
	 * @return desire point y-value
	 */
	public static double linInt(double x1, double y1, double x2, double y2,
			double x) {
		double result = -1.0;
		double[] mn = new double[2];
		// insucces
		mn[0] = -1.0;// m
		mn[1] = -1.0;// n
		double num = x1 - x2;
		if (num != 0.0) {
			mn[0] = (y1 - y2) / num;
			mn[1] = (x1 * y2 - y1 * x2) / num;
			result = mn[0] * x + mn[1];
		}
		return result;
	}

	/**
	 * Linear interpolation based on the natural logarithms of input values. If fail, then direct linear interpolation is used.
	 * @param x1 first point x-value
	 * @param y1 first point y-value
	 * @param x2 second point x-value
	 * @param y2 second point y-value
	 * @param x desire point x-value
	 * @return desire point y-value
	 */
	public static double linLogInt(double x1, double y1, double x2, double y2,
			double x) {
		double result = -1.0;
		double[] mn = new double[2];
		double x11 = Math.log(x1);// Ln (x1)...nuclides.net method of
									// interpolation!!!Not good at 0=>NaN
		double x22 = Math.log(x2);

		double yy1 = y1;
		double yy2 = y2;
		if (y1 == 0) {
			yy1 = -1.0E+30;
		} else {
			yy1 = Math.log(y1);
		}// -inf
		if (y2 == 0) {
			yy2 = -1.0E+30;
		} else {
			yy2 = Math.log(y2);
		}// -inf
		double y11 = yy1;// Math.log(y1);
		double y22 = yy2;// Math.log(y2);
		double xx = Math.log(x);
		// insucces
		mn[0] = -1.0;// m
		mn[1] = -1.0;// n
		double num = x11 - x22;
		if (num != 0.0) {
			mn[0] = (y11 - y22) / num;
			mn[1] = (x11 * y22 - y11 * x22) / num;
			result = mn[0] * xx + mn[1];
			result = Math.exp(result);
		} else {
			result = linInt(x1, y1, x2, y2, x);
		}
		return result;
	}

	
}
