package danfulea.math.numerical;

/**
 * Function interface
 * 
 * @author Dan Fulea <br>
 *         Created on 05.10.2006.
 */
public interface Function {
	void printSequence(String s);// for printing

	double F(double x);// y=f(x)

	double[] FD(double x);// 0=->y=f(x) and 1-> dy=f'(x)

	double MF(double[] x);// y=f(x1,x2,...)

	double[] DMF(double[] x);// y'1=df(x1,x2,...)/dx1;...the vector gradient
								// df[1..n] evaluated at the input point x
	// ==============3D func====================

	double F3D(double x, double y, double z);

	double yy1(double x);

	double yy2(double x);

	double z1(double x, double y);

	double z2(double x, double y);

	// end 3d====================================
	// non linear equation systems: root finding
	double[] vecfunc(int n, double[] x);

	// ==========================================
	// The user supplies a routine funcs(x,afunc,ma) that
	// returns the ma basis functions evaluated at x = x in the array
	// afunc[1..ma].
	double[] aF(double x, int ma);// poli fit

	double fdf(double x, double[] a, double[] dyda, int na);// gauss fit etc.
															// nonliniar
	// ===========================================================

	double[] derivF(double x, double[] y);

	// ============2point
	void load(double x1, double[] v, double[] y);

	void load1(double x1, double[] v, double[] y);

	void load2(double x2, double[] v, int nn2, double[] y);

	// double[] score(double x2, double[] y);
	void score(double x2, double[] y, double[] f);

	void difeq(int k, int k1, int k2, int jsf, int is1, int isf, int indexv[],
			int ne, double[][] s, double[][] y);

	// ================================================
	double g(double t);// g(t)=FREDHOLM

	double ak(double t, double s);// KERNEL

	double g(int k, double t);// voltera

	double ak(int k, int l, double t, double s);// voltera

	void kermom(double[] w, double y, int m);
}
