package danfulea.math.numerical;

/**
 * Class for solving integral equations.
 * Based on literature, including Numerical Recipes in C (Cambridge Univ.).
 * @author Dan Fulea, 26 OCT. 2006
 */
public class IntegralEq {
	public static Function func;
	public static boolean failB = false;
	public static String failS = "";

	public static double PI = 3.14159265;

	/*
	 * There is a close correspondence between linear integral equations, which
	 * specify linear, integral relations among functions in an
	 * infinite-dimensional function space, and plain old linear equations,
	 * which specify analogous relations among vectors in a finite-dimensional
	 * vector space. Because this correspondence lies at the heart of most
	 * computational algorithms, it is worth making it explicit as we recall how
	 * integral equations are classified. Fredholm equations involve definite
	 * integrals with fixed upper and lower limits. An inhomogeneous Fredholm
	 * equation of the first kind has the form
	 * 
	 * g(t) = Integr a->b of K(t, s)f(s) ds (18.0.1)
	 * 
	 * Here f(t) is the unknown function to be solved for, while g(t) is a known
	 * “right-hand side.” (In integral equations, for some odd reason, the
	 * familiar “right-hand side” is conventionally written on the left!) The
	 * function of two variables, K(t, s) is called the kernel. Equation
	 * (18.0.1) is analogous to the matrix equation K · f = g (18.0.2) whose
	 * solution is f = K-1·g, whereK-1 is the matrix inverse. Like equation
	 * (18.0.2), equation (18.0.1) has a unique solution whenever g is nonzero
	 * (the homogeneous case with g = 0 is almost never useful) and K is
	 * invertible. However, as we shall see, this latter condition is as often
	 * the exception as the rule. The analog of the finite-dimensional
	 * eigenvalue problem (K - sigma 1) · f = g is called a Fredholm equation of
	 * the second kind, usually written
	 * 
	 * f(t) = integral a->b K(t, s)f(s) ds + g(t) (18.0.4)
	 * 
	 * Again, the notational conventions do not exactly correspond: ? in
	 * equation (18.0.4) is 1/? in (18.0.3), while g is -g/?. If g (or g) is
	 * zero, then the equation is said to be homogeneous. If the kernel K(t, s)
	 * is bounded, then, like equation (18.0.3), equation (18.0.4) has the
	 * property that its homogeneous form has solutions for at most a
	 * denumerably infinite set ? = ?n, n = 1, 2, . . . , the eigenvalues. The
	 * corresponding solutions fn(t) are the eigenfunctions. The eigenvalues are
	 * real if the kernel is symmetric. In the inhomogeneous case of nonzero g
	 * (or g), equations (18.0.3) and (18.0.4) are soluble except when ? (or ?)
	 * is an eigenvalue — because the integral operator (or matrix) is singular
	 * then. In integral equations this dichotomy is called the Fredholm
	 * alternative. Fredholm equations of the first kind are often extremely
	 * ill-conditioned. Applying the kernel to a function is generally a
	 * smoothing operation, so the solution, which requires inverting the
	 * operator, will be extremely sensitive to small changes or errors in the
	 * input. Smoothing often actually loses information, and there is no way to
	 * get it back in an inverse operation. Specialized methods have been
	 * developed for such equations, which are often called inverse problems. In
	 * general, a method must augment the information given with some prior
	 * knowledge of the nature of the solution. This prior knowledge is then
	 * used, in one way or another, to restore lost information. We will
	 * introduce such techniques in §18.4. The Volterra equation of the first
	 * kind g(t) =  t a K(t, s)f(s) ds (18.0.6) has as its analog the matrix
	 * equation (now written out in components) k j =1 Kkjfj = gk (18.0.7)
	 * Comparing with equation (18.0.2), we see that the Volterra equation
	 * corresponds to a matrix K that is lower (i.e., left) triangular, with
	 * zero entries above the diagonal. As we know from Chapter 2, such matrix
	 * equations are trivially soluble by forward substitution. Techniques for
	 * solving Volterra equations are similarly straightforward. When
	 * experimental measurement noise does not dominate, Volterra equations of
	 * the first kind tend not to be ill-conditioned; the upper limit to the
	 * integral introduces a sharp step that conveniently spoils any smoothing
	 * properties of the kernel. The Volterra equation of the second kind is
	 * written f(t) =  t a K(t, s)f(s) ds + g(t) (18.0.8) whose matrix analog
	 * is the equation (K - 1) · f = g (18.0.9) with K lower triangular. The
	 * reason there is no ? in these equations is that (i) in the inhomogeneous
	 * case (nonzero g) it can be absorbed into K, while (ii) in the homogeneous
	 * case (g = 0), it is a theorem that Volterra equations of the second kind
	 * with bounded kernels have no eigenvalues with square-integrable
	 * eigenfunctions. We have specialized our definitions to the case of linear
	 * integral equations. The integrand in a nonlinear version of equation
	 * (18.0.1) or (18.0.6) would be K(t, s, f(s)) instead of K(t, s)f(s); a
	 * nonlinear version of equation (18.0.4) or (18.0.8) would have an
	 * integrandK(t, s, f(t), f(s)). Nonlinear Fredholm equations are
	 * considerably more complicated than their linear counterparts.
	 * Fortunately, they do not occur as frequently in practice and we shall by
	 * and large ignore them in this chapter. By contrast, solving nonlinear
	 * Volterra equations usually involves only a slight modification of the
	 * algorithm for linear equations, as we shall see. Almost all methods for
	 * solving integral equations numerically make use of quadrature rules,
	 * frequently Gaussian quadratures. This would be a good time for you to go
	 * back and review §4.5, especially the advanced material towards the end of
	 * that section. In the sections that follow,we first discuss Fredholm
	 * equations of the second kind with smooth kernels (§18.1). Nontrivial
	 * quadrature rules come into the discussion, but we will be dealing with
	 * well-conditioned systems of equations. We then return to Volterra
	 * equations (§18.2), and find that simple and straightforward methods are
	 * generally satisfactory for these equations
	 */
	/**
	 * Solves a linear Fredholm equation of the second kind. On input, a and b are the limits of 
	 * integration, and n is the number of points to use in the Gaussian quadrature. g and ak are 
	 * user-supplied external functions that respectively return g(t) and LAMBDAK(t,s). The routine returns 
	 * arrays t[1..n] and f[1..n] containing the abscissas ti of the Gaussian quadrature and the 
	 * solution f at these abscissas. Also returned is the array w[1..n] of Gaussian weights for use 
	 * with the Nystrom interpolation routine fredin.
	 * @param n n
	 * @param a a
	 * @param b b
	 * @param t t
	 * @param f f
	 * @param w w
	 */
	public static void fred2(int n, double a, double b, double[] t, double[] f,
			double[] w)// ,
	// float (*g)(float), float (*ak)(float, float))
	// Solves a linear Fredholm equation of the second kind. On input, a and b
	// are the limits of
	// integration, and n is the number of points to use in the Gaussian
	// quadrature. g and ak are
	// user-supplied external functions that respectively return g(t) and ?K(t,
	// s). The routine returns
	// arrays t[1..n] and f[1..n] containing the abscissas ti of the Gaussian
	// quadrature and the
	// solution f at these abscissas. Also returned is the array w[1..n] of
	// Gaussian weights for use
	// with the Nystrom interpolation routine fredin.
	{
		// void gauleg(float x1, float x2, float x[], float w[], int n);
		// void lubksb(float **a, int n, int *indx, float b[]);
		// void ludcmp(float **a, int n, int *indx, float *d);
		int i = 0;
		int j = 0;
		int[] indx = new int[n];
		// double d=0.0;
		double[][] omk = new double[n][n];
		// indx=ivector(1,n);
		// omk=matrix(1,n,1,n);
		Integrator.gauleg(a, b, t, w, n); // Replace gauleg with another routine
											// if not using
											// Gauss-Legendre quadrature.
		for (i = 1; i <= n; i++) {
			for (j = 1; j <= n; j++) // Form 1 - ?K.
			{
				// omk[i][j]=(float)(i == j)-(*ak)(t[i],t[j])*w[j];
				if (i == j)// true=1, false==0
					omk[i - 1][j - 1] = 1.0 - func.ak(t[i - 1], t[j - 1])
							* w[j - 1];
				else
					omk[i - 1][j - 1] = -func.ak(t[i - 1], t[j - 1]) * w[j - 1];
			}
			f[i - 1] = func.g(t[i - 1]);// f[i]=(*g)(t[i]);
		}
		LinAEq.ludcmp(omk, n, indx);// ,&d); Solve linear equations.
		LinAEq.lubksb(omk, n, indx, f);
		// free_matrix(omk,1,n,1,n);
		// free_ivector(indx,1,n);
	}

	/**
	 * Given arrays t[1..n] and w[1..n] containing the abscissas and weights of the Gaussian 
	 * quadrature, and given the solution array f[1..n] from fred2, this function returns the value of 
	 * f at x using the Nystrom interpolation formula. On input, a and b are the limits of integration, 
	 * and n is the number of points used in the Gaussian quadrature. g and ak are user-supplied 
	 * external functions that respectively return g(t) and LAMBDAK(t, s).
	 * @param x x
	 * @param n n
	 * @param a a
	 * @param b b
	 * @param t t
	 * @param f f
	 * @param w w
	 * @return the result
	 */
	public static double fredin(double x, int n, double a, double b,
			double[] t, double[] f, double[] w)// , float (*g)(float), float
												// (*ak)(float, float))
	// Given arrays t[1..n] and w[1..n] containing the abscissas and weights of
	// the Gaussian
	// quadrature, and given the solution array f[1..n] from fred2, this
	// function returns the value of
	// f at x using the Nystrom interpolation formula. On input, a and b are the
	// limits of integration,
	// and n is the number of points used in the Gaussian quadrature. g and ak
	// are user-supplied
	// external functions that respectively return g(t) and ?K(t, s).
	{
		int i;
		double sum = 0.0;
		for (i = 1; i <= n; i++)
			// sum += (*ak)(x,t[i])*w[i]*f[i];
			sum += func.ak(x, t[i - 1]) * w[i - 1] * f[i - 1];
		return func.g(x) + sum;
	}

	/*
	 * Volterra Equations
	 */
	/**
	 * Solves a set of m linear Volterra equations of the second kind using the extended trapezoidal rule. 
	 * On input, t0 is the starting point of the integration and n-1 is the number of steps of size h to 
	 * be taken. g(k,t) is a user-supplied external function that returns gk(t), while ak(k,l,t,s) 
	 * is another user-supplied external function that returns the (k, l) element of the matrix K(t, s). 
	 * The solution is returned in f[1..m][1..n], with the corresponding abscissas in t[1..n].
	 * @param n n
	 * @param m m
	 * @param t0 t0
	 * @param h h
	 * @param t t
	 * @param f f
	 */
	public static void voltra(int n, int m, double t0, double h, double[] t,
			double[][] f)// ,
	// float (*g)(int, float), float (*ak)(int, int, float, float))
	// Solves a set of m linear Volterra equations of the second kind using the
	// extended trapezoidal rule.
	// On input, t0 is the starting point of the integration and n-1 is the
	// number of steps of size h to
	// be taken. g(k,t) is a user-supplied external function that returns gk(t),
	// while ak(k,l,t,s)
	// is another user-supplied external function that returns the (k, l)
	// element of the matrix K(t, s).
	// The solution is returned in f[1..m][1..n], with the corresponding
	// abscissas in t[1..n].
	{
		// void lubksb(float **a, int n, int *indx, float b[]);
		// void ludcmp(float **a, int n, int *indx, float *d);
		int i = 0;
		int j = 0;
		int k = 0;
		int l = 0;
		int[] indx = new int[m];
		// double d=0.0;
		double sum = 0.0;
		double[][] a = new double[m][m];
		double[] b = new double[m];
		// indx=ivector(1,m);
		// a=matrix(1,m,1,m);
		// b=vector(1,m);
		t[0] = t0;// t[1]=t0;
		for (k = 1; k <= m; k++)
			f[k - 1][0] = func.g(k, t[0]);// f[k][1]=(*g)(k,t[1]); Initialize.
		for (i = 2; i <= n; i++) {// Take a step h.
			t[i - 1] = t[i - 2] + h;// t[i]=t[i-1]+h;
			for (k = 1; k <= m; k++) {
				// sum=(*g)(k,t[i]); Accumulate right-hand side of linear
				// equations in sum.
				sum = func.g(k, t[i - 1]);
				for (l = 1; l <= m; l++) {
					// sum += 0.5*h*(*ak)(k,l,t[i],t[1])*f[l][1];
					sum += 0.5 * h * func.ak(k, l, t[i - 1], t[0])
							* f[l - 1][0];
					for (j = 2; j < i; j++)
						// sum += h*(*ak)(k,l,t[i],t[j])*f[l][j];
						sum += h * func.ak(k, l, t[i - 1], t[j - 1])
								* f[l - 1][j - 1];
					// a[k][l]=(k == l)-0.5*h*(*ak)(k,l,t[i],t[i]); Left-hand
					// side goesin matrix a.
					if (k == l)
						a[k - 1][l - 1] = 1.0 - 0.5 * h
								* func.ak(k, l, t[i - 1], t[i - 1]);
					else
						a[k - 1][l - 1] = -0.5 * h
								* func.ak(k, l, t[i - 1], t[i - 1]);
				}
				b[k - 1] = sum;// b[k]=sum;
			}
			LinAEq.ludcmp(a, m, indx);// ,&d); Solve linear equations.
			LinAEq.lubksb(a, m, indx, b);
			for (k = 1; k <= m; k++)
				f[k - 1][i - 1] = b[k - 1];// f[k][i]=b[k];
		}
		// free_vector(b,1,m);
		// free_matrix(a,1,m,1,m);
		// free_ivector(indx,1,m);
	}

	/*
	 * Integral Equations with Singular Kernels
	 */
	/**
	 * Constructs in wghts[1..n] weights for the n-point equal-interval quadrature from 0 to (n-1)h 
	 * of a function f(x) times an arbitrary (possibly singular) weight function w(x) whose indefinite integral 
	 * moments Fn(y) are provided by the user-supplied routine kermom.
	 * @param wghts wghts
	 * @param n n
	 * @param h h
	 */ 
	public static void wwghts(double[] wghts, int n, float h)// ,
	// void (*kermom)(double [], double ,int))
	// Constructs in wghts[1..n] weights for the n-point equal-interval
	// quadrature from 0 to (n-1)h
	// of a function f(x) times an arbitrary (possibly singular) weight function
	// w(x) whose indefiniteintegral
	// moments Fn(y) are provided by the user-supplied routine kermom.
	{
		int j = 0;
		int k = 0;
		// double wold[5],wnew[5],w[5],hh,hi,c,fac,a,b;
		double hh = 0.0;
		double hi = 0.0;
		double c = 0.0;
		double fac = 0.0;
		double a = 0.0;
		double b = 0.0;
		double[] wold = new double[5];// [0] not defined!!
		double[] wnew = new double[5];
		double[] w = new double[5];
		// Double precision on internal calculations even though the interface
		// is in single precision.
		hh = h;
		hi = 1.0 / hh;
		for (j = 1; j <= n; j++)
			wghts[j - 1] = 0.0;// wghts[j]=0.0;
		// Zero all the weights so we can sum into them.
		func.kermom(wold, 0.0, 4); // Evaluate indefinite integrals at lower
									// end.
		if (n >= 4) {// Use highest available order.
			b = 0.0; // For another problem, you might change this lower limit.
			for (j = 1; j <= n - 3; j++) {
				c = j - 1;// This is called k in equation (18.3.5).
				a = b; // Set upper and lower limits for this step.
				b = a + hh;
				if (j == n - 3)
					b = (n - 1) * hh; // Last interval: go all the way to end.
				func.kermom(wnew, b, 4);
				for (fac = 1.0, k = 1; k <= 4; k++, fac *= hi)
					// Equation (18.3.4).
					w[k] = (wnew[k] - wold[k]) * fac;
				// wghts[j] += //( Equation (18.3.5).
				wghts[j - 1] += (((c + 1.0) * (c + 2.0) * (c + 3.0) * w[1]
						- (11.0 + c * (12.0 + c * 3.0)) * w[2] + 3.0
						* (c + 2.0) * w[3] - w[4]) / 6.0);
				wghts[j] += // wghts[j+1] +=
				((-c * (c + 2.0) * (c + 3.0) * w[1]
						+ (6.0 + c * (10.0 + c * 3.0)) * w[2] - (3.0 * c + 5.0)
						* w[3] + w[4]) * 0.5);
				wghts[j + 1] += // wghts[j+2] +=
				((c * (c + 1.0) * (c + 3.0) * w[1]
						- (3.0 + c * (8.0 + c * 3.0)) * w[2] + (3.0 * c + 4.0)
						* w[3] - w[4]) * 0.5);
				wghts[j + 2] += // wghts[j+3] +=
				((-c * (c + 1.0) * (c + 2.0) * w[1]
						+ (2.0 + c * (6.0 + c * 3.0)) * w[2] - 3.0 * (c + 1.0)
						* w[3] + w[4]) / 6.0);
				for (k = 1; k <= 4; k++)
					wold[k] = wnew[k];// Reset lower limits for moments.
			}
		} else if (n == 3) {// Lower-order cases; not recommended.
			func.kermom(wnew, hh + hh, 3);
			w[1] = wnew[1] - wold[1];
			w[2] = hi * (wnew[2] - wold[2]);
			w[3] = hi * hi * (wnew[3] - wold[3]);

			wghts[0] = w[1] - 1.5 * w[2] + 0.5 * w[3];// wghts[1]=w[1]-1.5*w[2]+0.5*w[3];
			wghts[1] = 2.0 * w[2] - w[3];// wghts[2]=2.0*w[2]-w[3];
			wghts[2] = 0.5 * (w[3] - w[2]);// wghts[3]=0.5*(w[3]-w[2]);
		} else if (n == 2) {
			func.kermom(wnew, hh, 2);
			// wghts[1]=wnew[1]-wold[1]-(wghts[2]=hi*(wnew[2]-wold[2]));
			wghts[0] = wnew[1] - wold[1]
					- (wghts[1] = hi * (wnew[2] - wold[2]));
		}
	}

	/*
	 * Worked Example: A Diagonally Singular Kernel As a particular example,
	 * consider the integral equation f(x) +  ? 0 K(x, y)f(y)dy = sinx
	 * (18.3.13) with the (arbitrarily chosen) nasty kernel K(x, y) = cos x cos
	 * y × -ln(x - y) y < x ?y -x y? x (18.3.14) which has a logarithmic
	 * singularity on the left of the diagonal, combined with a square-root
	 * discontinuity on the right. The first step is to do (analytically, in
	 * this case) the required moment integrals over the singular part of the
	 * kernel, equation (18.3.12). Since these integrals are done at a fixed
	 * value of x, we can use x as the lower limit. For any specified value of
	 * y, the required indefinite integral is then either
	 * 
	 * extern double x; Defined in quadmx. void kermom(double w[], double y, int
	 * m) Returns in w[1..m] the first m indefinite-integral moments of one row
	 * of the singular part of the kernel. (For this example, m is hard-wired to
	 * be 4.) The input variable y labels the column, while the global variable
	 * x is the row. We can take x as the lower limit of integration. Thus, we
	 * return the moment integrals either purely to the left or purely to the
	 * right of the diagonal. { double d,df,clog,x2,x3,x4,y2; if (y >= x) {
	 * d=y-x; df=2.0*sqrt(d)*d; w[1]=df/3.0; w[2]=df*(x/3.0+d/5.0);
	 * w[3]=df*((x/3.0 + 0.4*d)*x + d*d/7.0); w[4]=df*(((x/3.0 + 0.6*d)*x +
	 * 3.0*d*d/7.0)*x+d*d*d/9.0); } else { x3=(x2=x*x)*x; x4=x2*x2; y2=y*y;
	 * d=x-y; w[1]=d*((clog=log(d))-1.0); w[2] =
	 * -0.25*(3.0*x+y-2.0*clog*(x+y))*d;
	 * w[3]=(-11.0*x3+y*(6.0*x2+y*(3.0*x+2.0*y)) +6.0*clog*(x3-y*y2))/18.0;
	 * w[4]=(-25.0*x4+y*(12.0*x3+y*(6.0*x2+y*
	 * (4.0*x+3.0*y)))+12.0*clog*(x4-(y2*y2)))/48.0; } }
	 * 
	 * #define PI 3.14159265 double x; Communicates with kermom. void
	 * quadmx(float **a, int n) Constructs in a[1..n][1..n] the quadrature
	 * matrix for an example Fredholm equation of the second kind. The
	 * nonsingular part of the kernel is computed within this routine, while the
	 * quadrature weights which integrate the singular part of the kernel are
	 * obtained via calls to wwghts. An external routine kermom, which supplies
	 * indefinite-integral moments of the singular part of the kernel, is passed
	 * to wwghts. { void kermom(double w[], double y, int m); void wwghts(float
	 * wghts[], int n, float h, void (*kermom)(double [], double ,int)); int
	 * j,k; float h,*wt,xx,cx; wt=vector(1,n); h=PI/(n-1); for (j=1;j<=n;j++) {
	 * x=xx=(j-1)*h; Put x in global variable for use by kermom.
	 * wwghts(wt,n,h,kermom); cx=cos(xx); Part of nonsingular kernel. for
	 * (k=1;k<=n;k++) a[j][k]=wt[k]*cx*cos((k-1)*h); Put together all the pieces
	 * of the kernel. ++a[j][j]; Since equation of the second kind, there is
	 * diagonal piece independent of h. } free_vector(wt,1,n); } Finally, we
	 * solve the linear system for any particular right-hand side, here sin x.
	 * #include <stdio.h> #include <math.h> #include "nrutil.h" #define PI
	 * 3.14159265 #define N 40 Here the size of the grid is specified. int
	 * main(void) // Program fredex This sample program shows how to solve a
	 * Fredholm equation of the second kind using the product Nystrom method and
	 * a quadrature rule especially constructed for a particular, singular,
	 * kernel. { void lubksb(float **a, int n, int *indx, float b[]); void
	 * ludcmp(float **a, int n, int *indx, float *d); void quadmx(float **a, int
	 * n); float **a,d,*g,x; int *indx,j; indx=ivector(1,N); a=matrix(1,N,1,N);
	 * g=vector(1,N); quadmx(a,N); Make the quadrature matrix; all the action is
	 * here. ludcmp(a,N,indx,&d); Decompose the matrix. for (j=1;j<=N;j++)
	 * g[j]=sin((j-1)*PI/(N-1)); Construct the right hand side, here sin x.
	 * lubksb(a,N,indx,g); Backsubstitute. for (j=1;j<=N;j++) { Write out the
	 * solution. x=(j-1)*PI/(N-1); printf("%6.2d %12.6f %12.6f\n",j,x,g[j]); }
	 * free_vector(g,1,N); free_matrix(a,1,N,1,N); free_ivector(indx,1,N);
	 * return 0; }
	 */
}
