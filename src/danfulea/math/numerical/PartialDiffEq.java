package danfulea.math.numerical;

/**
 * Partial diff eq.
 * Based on literature, including Numerical Recipes in C (Cambridge Univ.).
 * @author Dan Fulea, 28 OCT. 2006
 */
public class PartialDiffEq {
	// public static Function func;
	public static boolean failB = false;
	public static String failS = "";

	public static int MAXITS = 1000;
	public static double EPS = 1.0e-5;

	/*
	 * Relaxation Methods for Boundary Value Problems
	 * 
	 * As we mentioned in §19.0, relaxation methods involve splitting the sparse
	 * matrix that arises from finite differencing and then iterating until a
	 * solution is found. There is another way of thinking about relaxation
	 * methods that is somewhat more physical. Suppose we wish to solve the
	 * elliptic equation Lu = rho where L represents some elliptic operator and
	 * ? is the source term. Rewrite the equation as a diffusion equation,
	 * 
	 * dpart u/dpart t = Lu - rho
	 * 
	 * An initial distribution u relaxes to an equilibrium solution as t › ?.
	 * This equilibrium has all time derivatives vanishing. Therefore it is the
	 * solution of the original elliptic problem (19.5.1). We see that all the
	 * machinery of §19.2, on diffusive initial value equations, can be brought
	 * to bear on the solution of boundary value problems by relaxation methods.
	 * 
	 * Let us apply this idea to our model problem (19.0.3). The diffusion
	 * equation is dpart u/dpart t = dpart2u/dpart x2 + dpart2u/dparty2 - rho
	 * ......................... un+1 j,l = 1 4 unj +1,l + un+1 j-1,l + un
	 * j,l+1 + un+1 j,l-1- ?2 4 ?j,l
	 * 
	 * A · x = b A = L + D + U The Gauss-Seidel method, equation (19.5.6),
	 * corresponds to the matrix decomposition (L + D) · x(r) = -U · x(r-1) + b
	 * 
	 * Successive Overrelaxation (SOR) Consider a general second-order elliptic
	 * equation in x and y, finite differenced on a square as for our model
	 * equation. Corresponding to each row of the matrix A is an equation of the
	 * form
	 * 
	 * aj,luj+1,l + bj,luj-1,l + cj,luj,l+1 + dj,luj,l-1 + ej,luj,l = fj,l
	 * (19.5.25)
	 * 
	 * Here we give a routine for SOR with Chebyshev acceleration. Our advice is
	 * to use SOR for trivial problems (e.g., 20 × 20), or for solving a larger
	 * problem once only, where ease of programming outweighs expense of
	 * computer time. Occasionally, the sparse matrix methods of §2.7 are useful
	 * for solving a set of difference equations directly. For production
	 * solution of large elliptic problems, however, multigrid is now almost
	 * always the method of choice.
	 */

	/**
	 * Successive overrelaxation solution for general 2nd order elliptic equation in x and y with Chebyshev acceleration. a, b, c, d, e, and f are input as the coefficients of the equation, each 
	 * dimensioned to the grid [1..jmax][1..jmax]. u is input as the initial guess to the solution, usually zero, and returns 
	 * with the final value. rjac is input as the spectral radius of the Jacobi iteration, or an estimate of it.
	 * @param a a
	 * @param b b
	 * @param c c
	 * @param d d
	 * @param e e
	 * @param f f
	 * @param u u
	 * @param jmax jmax
	 * @param rjac rjac
	 */
	public static void sor(double[][] a, double[][] b, double[][] c,
			double[][] d, double[][] e, double[][] f, double[][] u, int jmax,
			double rjac)
	// Successive overrelaxation solution of equation (19.5.25) with Chebyshev
	// acceleration. a, b, c,
	// d, e, and f are input as the coefficients of the equation, each
	// dimensioned to the grid size
	// [1..jmax][1..jmax]. u is input as the initial guess to the solution,
	// usually zero, and returns
	// with the final value. rjac is input as the spectral radius of the Jacobi
	// iteration, or an estimate
	// of it.
	{
		// void nrerror(char error_text[]);
		failB = false;
		int ipass = 0;
		int j = 0;
		int jsw = 0;
		int l = 0;
		int lsw = 0;
		int n = 0;
		double anorm = 0.0;
		double anormf = 0.0;
		double omega = 1.0;
		double resid = 0.0;
		// Double precision is a good idea for jmax bigger than about 25.
		for (j = 2; j < jmax; j++)
			// Compute initial norm of residual and terminate iteration when
			// norm has been reduced by
			// a factor EPS.
			for (l = 2; l < jmax; l++)
				// anormf += fabs(f[j][l]); Assumes initial u is zero.
				anormf += Math.abs(f[j - 1][l - 1]);
		for (n = 1; n <= MAXITS; n++) {
			anorm = 0.0;
			jsw = 1;
			for (ipass = 1; ipass <= 2; ipass++) {// Odd-even ordering.
				lsw = jsw;
				for (j = 2; j < jmax; j++) {
					for (l = lsw + 1; l < jmax; l += 2) {
						resid = a[j - 1][l - 1] * u[j][l - 1]// resid=a[j][l]*u[j+1][l]
								+ b[j - 1][l - 1] * u[j - 2][l - 1]// +b[j][l]*u[j-1][l]
								+ c[j - 1][l - 1] * u[j - 1][l]// +c[j][l]*u[j][l+1]
								+ d[j - 1][l - 1] * u[j - 1][l - 2]// +d[j][l]*u[j][l-1]
								+ e[j - 1][l - 1] * u[j - 1][l - 1]// +e[j][l]*u[j][l]
								- f[j - 1][l - 1];// -f[j][l];
						anorm += Math.abs(resid);
						// u[j][l] -= omega*resid/e[j][l];
						u[j - 1][l - 1] -= omega * resid / e[j - 1][l - 1];
					}
					lsw = 3 - lsw;
				}
				jsw = 3 - jsw;
				omega = (n == 1 && ipass == 1 ? 1.0 / (1.0 - 0.5 * rjac * rjac)
						: 1.0 / (1.0 - 0.25 * rjac * rjac * omega));
			}
			if (anorm < EPS * anormf)
				return;
		}
		// nrerror("MAXITS exceeded");
		failB = true;
		failS = "MAXITS exceeded";
		return;
	}

	/*
	 * Multigrid Methods for Boundary Value problems Practical multigrid methods
	 * were first introduced in the 1970s by Brandt. These methods can solve
	 * elliptic PDEs discretized on N grid points in O(N) operations. The
	 * “rapid” direct elliptic solvers discussed in §19.4 solve special kinds of
	 * elliptic equations in O(N logN) operations. The numerical coefficients in
	 * these estimates are such that multigrid methods are comparable to the
	 * rapid methods in execution speed. Unlike the rapid methods, however, the
	 * multigrid methods can solve general elliptic equations with nonconstant
	 * coefficients with hardly any loss in efficiency. Even nonlinear equations
	 * can be solved with comparable speed. Unfortunately there is not a single
	 * multigrid algorithm that solves all elliptic problems. Rather there is a
	 * multigrid technique that provides the framework for solving these
	 * problems. You have to adjust the various components of the algorithm
	 * within this framework to solve your specific problem. We can only give a
	 * brief introduction to the subject here. In particular, we will give two
	 * sample multigrid routines, one linear and one nonlinear. By following
	 * these prototypes and by perusing the references [1-4], you should be able
	 * to develop routines to solve your own problems. There are two related,
	 * but distinct, approaches to the use of multigrid techniques. The first,
	 * termed “the multigrid method,” is a means for speeding up the convergence
	 * of a traditional relaxation method, as defined by you on a grid of
	 * pre-specified fineness. In this case, you need define your problem (e.g.,
	 * evaluate its source terms) only on this grid. Other, coarser, grids
	 * defined by the method can be viewed as temporary computational adjuncts.
	 * The second approach, termed (perhaps confusingly) “the full multigrid
	 * (FMG) method,” requires you to be able to define your problem on grids of
	 * various sizes (generally by discretizing the same underlying PDE into
	 * different-sized sets of finitedifference equations). In this approach,
	 * the method obtains successive solutions on finer and finer grids. You can
	 * stop the solution either at a pre-specified fineness, or you can monitor
	 * the truncation error due to the discretization, quitting only when it is
	 * tolerably small. In this section we will first discuss the “multigrid
	 * method,” then use the concepts developed to introduce the FMG method. The
	 * latter algorithm is the one that we implement in the accompanying
	 * programs. ............ Full Multigrid Algorithm Smoothing, Restriction,
	 * and Prolongation Operators The most popular smoothing method, and the one
	 * you should try first, is Gauss-Seidel, since it usually leads to a good
	 * convergence rate. If we order the mesh points from 1 to N, then the
	 * Gauss-Seidel scheme is ui = - N  j =1 j=i Lijuj - fi 1 Lii i = 1, . .
	 * .,N
	 * 
	 * 
	 * #define NPRE 1 Number of relaxation sweeps before . . . #define NPOST 1 .
	 * . . and after the coarse-grid correction is computed. #define NGMAX 15
	 * 
	 * void mglin(double **u, int n, int ncycle) Full Multigrid Algorithm for
	 * solution of linear elliptic equation, @@@@@@@@@ here the model problem
	 * (19.0.6).@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ On input u[1..n][1..n]
	 * contains the right-hand side rho, while on output it returns the
	 * solution. The dimension n must be of the form 2j + 1 for some integer j.
	 * (j is actually the number of grid levels used in the solution, called ng
	 * below.) ncycle is the number of V-cycles to be used at each level. { void
	 * addint(double **uf, double **uc, double **res, int nf); void copy(double
	 * **aout, double **ain, int n); void fill0(double **u, int n); void
	 * interp(double **uf, double **uc, int nf); void relax(double **u, double
	 * **rhs, int n); void resid(double **res, double **u, double **rhs, int n);
	 * void rstrct(double **uc, double **uf, int nc); void slvsml(double **u,
	 * double **rhs); unsigned int j,jcycle,jj,jpost,jpre,nf,ng=0,ngrid,nn;
	 * double **ires[NGMAX+1],**irho[NGMAX+1],**irhs[NGMAX+1],**iu[NGMAX+1];
	 * nn=n; while (nn >>= 1) ng++; if (n != 1+(1L << ng))
	 * nrerror("n-1 must be a power of 2 in mglin."); if (ng > NGMAX)
	 * nrerror("increase NGMAX in mglin."); nn=n/2+1; ngrid=ng-1;
	 * irho[ngrid]=dmatrix(1,nn,1,nn); Allocate storage for r.h.s. on grid ng -
	 * 1, rstrct(irho[ngrid],u,nn); and fill it by restricting from the fine
	 * grid. while (nn > 3) { Similarly allocate storage and fill r.h.s. on all
	 * coarse grids. nn=nn/2+1; irho[--ngrid]=dmatrix(1,nn,1,nn);
	 * rstrct(irho[ngrid],irho[ngrid+1],nn); } nn=3; iu[1]=dmatrix(1,nn,1,nn);
	 * irhs[1]=dmatrix(1,nn,1,nn); slvsml(iu[1],irho[1]); Initial solution on
	 * coarsest grid. free_dmatrix(irho[1],1,nn,1,nn); ngrid=ng; for
	 * (j=2;j<=ngrid;j++) { Nested iteration loop. nn=2*nn-1;
	 * iu[j]=dmatrix(1,nn,1,nn); irhs[j]=dmatrix(1,nn,1,nn);
	 * ires[j]=dmatrix(1,nn,1,nn); interp(iu[j],iu[j-1],nn); Interpolate from
	 * coarse grid to next finer grid. copy(irhs[j],(j != ngrid ? irho[j] :
	 * u),nn); Set up r.h.s. for (jcycle=1;jcycle<=ncycle;jcycle++) { V-cycle
	 * loop. nf=nn; for (jj=j;jj>=2;jj--) { Downward stoke of the V. for
	 * (jpre=1;jpre<=NPRE;jpre++) Pre-smoothing. relax(iu[jj],irhs[jj],nf);
	 * resid(ires[jj],iu[jj],irhs[jj],nf); nf=nf/2+1;
	 * rstrct(irhs[jj-1],ires[jj],nf); Restriction of the residual is the next
	 * r.h.s. fill0(iu[jj-1],nf); Zero for initial guess in next relaxation. }
	 * slvsml(iu[1],irhs[1]); Bottom of V: solve on coarsest grid. nf=3; for
	 * (jj=2;jj<=j;jj++) { Upward stroke of V. nf=2*nf-1;
	 * addint(iu[jj],iu[jj-1],ires[jj],nf); Use res for temporary storage inside
	 * addint. for (jpost=1;jpost<=NPOST;jpost++) Post-smoothing.
	 * relax(iu[jj],irhs[jj],nf); } } } copy(u,iu[ngrid],n); Return solution in
	 * u. for (nn=n,j=ng;j>=2;j--,nn=nn/2+1) { free_dmatrix(ires[j],1,nn,1,nn);
	 * free_dmatrix(irhs[j],1,nn,1,nn); free_dmatrix(iu[j],1,nn,1,nn); if (j !=
	 * ng) free_dmatrix(irho[j],1,nn,1,nn); } free_dmatrix(irhs[1],1,3,1,3);
	 * free_dmatrix(iu[1],1,3,1,3); }
	 * 
	 * 
	 * void rstrct(double **uc, double **uf, int nc) Half-weighting restriction.
	 * nc is the coarse-grid dimension. The fine-grid solution is input in
	 * uf[1..2*nc-1][1..2*nc-1], the coarse-grid solution is returned in
	 * uc[1..nc][1..nc]. { int ic,iif,jc,jf,ncc=2*nc-1; for
	 * (jf=3,jc=2;jc<nc;jc++,jf+=2) { Interior points. for
	 * (iif=3,ic=2;ic<nc;ic++,iif+=2) {
	 * uc[ic][jc]=0.5*uf[iif][jf]+0.125*(uf[iif+1][jf]+uf[iif-1][jf]
	 * +uf[iif][jf+1]+uf[iif][jf-1]); } } for (jc=1,ic=1;ic<=nc;ic++,jc+=2) {
	 * Boundary points. uc[ic][1]=uf[jc][1]; uc[ic][nc]=uf[jc][ncc]; } for
	 * (jc=1,ic=1;ic<=nc;ic++,jc+=2) { uc[1][ic]=uf[1][jc];
	 * uc[nc][ic]=uf[ncc][jc]; } }
	 * 
	 * void interp(double **uf, double **uc, int nf) Coarse-to-fine prolongation
	 * by bilinear interpolation. nf is the fine-grid dimension. The coarsegrid
	 * solution is input as uc[1..nc][1..nc], where nc = nf/2 + 1. The fine-grid
	 * solution is returned in uf[1..nf][1..nf]. { int ic,iif,jc,jf,nc;
	 * nc=nf/2+1; for (jc=1,jf=1;jc<=nc;jc++,jf+=2) Do elements that are copies.
	 * for (ic=1;ic<=nc;ic++) uf[2*ic-1][jf]=uc[ic][jc]; for (jf=1;jf<=nf;jf+=2)
	 * Do odd-numbered columns, interpolating vertically. for
	 * (iif=2;iif<nf;iif+=2) uf[iif][jf]=0.5*(uf[iif+1][jf]+uf[iif-1][jf]); for
	 * (jf=2;jf<nf;jf+=2) Do even-numbered columns, interpolating horizontally.
	 * for (iif=1;iif <= nf;iif++)
	 * uf[iif][jf]=0.5*(uf[iif][jf+1]+uf[iif][jf-1]); }
	 * 
	 * 
	 * void addint(double **uf, double **uc, double **res, int nf) Does
	 * coarse-to-fine interpolation and adds result to uf. nf is the fine-grid
	 * dimension. The coarse-grid solution is input as uc[1..nc][1..nc], where
	 * nc = nf/2+1. The fine-grid solution is returned in uf[1..nf][1..nf].
	 * res[1..nf][1..nf] is used for temporary storage. { void interp(double
	 * **uf, double **uc, int nf); int i,j; interp(res,uc,nf); for
	 * (j=1;j<=nf;j++) for (i=1;i<=nf;i++) uf[i][j] += res[i][j]; }
	 * 
	 * void slvsml(double **u, double **rhs) Solution of the model problem on
	 * the coarsest grid, where h = 1 2 . The right-hand side is input in
	 * rhs[1..3][1..3] and the solution is returned in u[1..3][1..3]. { void
	 * fill0(double **u, int n); double h=0.5; fill0(u,3); u[2][2] =
	 * -h*h*rhs[2][2]/4.0; }
	 * 
	 * void relax(double **u, double **rhs, int n) Red-black Gauss-Seidel
	 * relaxation for model problem. Updates the current value of the solution
	 * u[1..n][1..n], using the right-hand side function rhs[1..n][1..n]. { int
	 * i,ipass,isw,j,jsw=1; double h,h2; h=1.0/(n-1); h2=h*h; for
	 * (ipass=1;ipass<=2;ipass++,jsw=3-jsw) { Red and black sweeps. isw=jsw; for
	 * (j=2;j<n;j++,isw=3-isw) for (i=isw+1;i<n;i+=2) Gauss-Seidel formula.
	 * u[i][j]=0.25*(u[i+1][j]+u[i-1][j]+u[i][j+1] +u[i][j-1]-h2*rhs[i][j]); } }
	 * 
	 * void resid(double **res, double **u, double **rhs, int n) Returns minus
	 * the residual for the model problem. Input quantities are u[1..n][1..n]
	 * and rhs[1..n][1..n], while res[1..n][1..n] is returned. { int i,j; double
	 * h,h2i; h=1.0/(n-1); h2i=1.0/(h*h); for (j=2;j<n;j++) Interior points. for
	 * (i=2;i<n;i++) res[i][j] = -h2i*(u[i+1][j]+u[i-1][j]+u[i][j+1]+u[i][j-1]-
	 * 4.0*u[i][j])+rhs[i][j]; for (i=1;i<=n;i++) Boundary points.
	 * res[i][1]=res[i][n]=res[1][i]=res[n][i]=0.0; }
	 * 
	 * void copy(double **aout, double **ain, int n) Copies ain[1..n][1..n] to
	 * aout[1..n][1..n]. { int i,j; for (i=1;i<=n;i++) for (j=1;j<=n;j++)
	 * aout[j][i]=ain[j][i]; }
	 * 
	 * void fill0(double **u, int n) Fills u[1..n][1..n] with zeros. { int i,j;
	 * for (j=1;j<=n;j++) for (i=1;i<=n;i++) u[i][j]=0.0; }
	 * 
	 * 
	 * 
	 * //================== Nonlinear Multigrid: The FAS Algorithm
	 * 
	 * 
	 * #define NPRE 1 Number of relaxation sweeps before . . . #define NPOST 1 .
	 * . . and after the coarse-grid correction is computed. #define ALPHA 0.33
	 * Relates the estimated truncation error to the norm of the residual.
	 * #define NGMAX 15 void mgfas(double **u, int n, int maxcyc) Full Multigrid
	 * Algorithm for FAS solution of nonlinear elliptic equation, here equation
	 * (19.6.44). On input u[1..n][1..n] contains the right-hand side ?, while
	 * on output it returns the solution. The dimension n must be of the form 2j
	 * + 1 for some integer j. (j is actually the number of grid levels used in
	 * the solution, called ng below.) maxcyc is the maximum number of V-cycles
	 * to be used at each level. { double anorm2(double **a, int n); void
	 * copy(double **aout, double **ain, int n); void interp(double **uf, double
	 * **uc, int nf); void lop(double **out, double **u, int n); void
	 * matadd(double **a, double **b, double **c, int n); void matsub(double
	 * **a, double **b, double **c, int n); void relax2(double **u, double
	 * **rhs, int n); void rstrct(double **uc, double **uf, int nc); void
	 * slvsm2(double **u, double **rhs); unsigned int
	 * j,jcycle,jj,jm1,jpost,jpre,nf,ng=0,ngrid,nn; double
	 * **irho[NGMAX+1],**irhs[NGMAX+1],**itau[NGMAX+1],
	 * *itemp[NGMAX+1],**iu[NGMAX+1]; double res,trerr; nn=n; while (nn >>= 1)
	 * ng++; if (n != 1+(1L << ng))
	 * nrerror("n-1 must be a power of 2 in mgfas."); if (ng > NGMAX)
	 * nrerror("increase NGMAX in mglin."); nn=n/2+1; ngrid=ng-1;
	 * irho[ngrid]=dmatrix(1,nn,1,nn); Allocate storage for r.h.s. on grid ng -
	 * 1, rstrct(irho[ngrid],u,nn); and fill it by restricting from the fine
	 * grid. while (nn > 3) { Similarly allocate storage and fill r.h.s. on all
	 * coarse grids. nn=nn/2+1; irho[--ngrid]=dmatrix(1,nn,1,nn);
	 * rstrct(irho[ngrid],irho[ngrid+1],nn); } nn=3; iu[1]=dmatrix(1,nn,1,nn);
	 * irhs[1]=dmatrix(1,nn,1,nn); itau[1]=dmatrix(1,nn,1,nn);
	 * itemp[1]=dmatrix(1,nn,1,nn); slvsm2(iu[1],irho[1]); Initial solution on
	 * coarsest grid. free_dmatrix(irho[1],1,nn,1,nn); ngrid=ng; for
	 * (j=2;j<=ngrid;j++) { Nested iteration loop. nn=2*nn-1;
	 * iu[j]=dmatrix(1,nn,1,nn); irhs[j]=dmatrix(1,nn,1,nn);
	 * itau[j]=dmatrix(1,nn,1,nn); itemp[j]=dmatrix(1,nn,1,nn);
	 * interp(iu[j],iu[j-1],nn); Interpolate from coarse grid to next finer
	 * grid. copy(irhs[j],(j != ngrid ? irho[j] : u),nn); Set up r.h.s. for
	 * (jcycle=1;jcycle<=maxcyc;jcycle++) { V-cycle loop. nf=nn; for
	 * (jj=j;jj>=2;jj--) { Downward stoke of the V. for
	 * (jpre=1;jpre<=NPRE;jpre++) Pre-smoothing. relax2(iu[jj],irhs[jj],nf);
	 * lop(itemp[jj],iu[jj],nf); Lh(uh). nf=nf/2+1; jm1=jj-1;
	 * rstrct(itemp[jm1],itemp[jj],nf); RLh(uh). rstrct(iu[jm1],iu[jj],nf);
	 * Ruh. lop(itau[jm1],iu[jm1],nf); LH(Ruh) stored temporarily in ?h.
	 * matsub(itau[jm1],itemp[jm1],itau[jm1],nf); Form ?h. if (jj == j)
	 * trerr=ALPHA*anorm2(itau[jm1],nf); Estimate truncation error ?.
	 * rstrct(irhs[jm1],irhs[jj],nf); fH.
	 * matadd(irhs[jm1],itau[jm1],irhs[jm1],nf); fH + ?h. }
	 * slvsm2(iu[1],irhs[1]); Bottom of V: Solve on coarsest grid. nf=3; for
	 * (jj=2;jj<=j;jj++) { Upward stroke of V. jm1=jj-1;
	 * rstrct(itemp[jm1],iu[jj],nf); Ruh.
	 * matsub(iu[jm1],itemp[jm1],itemp[jm1],nf); uH -Ruh. nf=2*nf-1;
	 * interp(itau[jj],itemp[jm1],nf); P(uH-Ruh) stored in ?h.
	 * matadd(iu[jj],itau[jj],iu[jj],nf); Form unew h . for
	 * (jpost=1;jpost<=NPOST;jpost++) Post-smoothing.
	 * relax2(iu[jj],irhs[jj],nf); } lop(itemp[j],iu[j],nf); Form residual dh.
	 * matsub(itemp[j],irhs[j],itemp[j],nf); res=anorm2(itemp[j],nf); if (res <
	 * trerr) break; No more V-cycles needed if residual small enough. } }
	 * copy(u,iu[ngrid],n); Return solution in u. for
	 * (nn=n,j=ng;j>=1;j--,nn=nn/2+1) { free_dmatrix(itemp[j],1,nn,1,nn);
	 * free_dmatrix(itau[j],1,nn,1,nn); free_dmatrix(irhs[j],1,nn,1,nn);
	 * free_dmatrix(iu[j],1,nn,1,nn); if (j != ng && j != 1)
	 * free_dmatrix(irho[j],1,nn,1,nn); } }
	 * 
	 * void relax2(double **u, double **rhs, int n) Red-black Gauss-Seidel
	 * relaxation for equation (19.6.44). The current value of the solution
	 * u[1..n][1..n] is updated, using the right-hand side function
	 * rhs[1..n][1..n]. { int i,ipass,isw,j,jsw=1; double foh2,h,h2i,res;
	 * h=1.0/(n-1); h2i=1.0/(h*h); foh2 = -4.0*h2i; for
	 * (ipass=1;ipass<=2;ipass++,jsw=3-jsw) { Red and black sweeps. isw=jsw; for
	 * (j=2;j<n;j++,isw=3-isw) { for (i=isw+1;i<n;i+=2) {
	 * res=h2i*(u[i+1][j]+u[i-1][j]+u[i][j+1]+u[i][j-1]-
	 * 4.0*u[i][j])+u[i][j]*u[i][j]-rhs[i][j]; u[i][j] -=
	 * res/(foh2+2.0*u[i][j]); Newton Gauss-Seidel formula. } } } }
	 * 
	 * 
	 * void slvsm2(double **u, double **rhs) Solution of equation (19.6.44) on
	 * the coarsest grid, where h = 1 2 . The right-hand side is input in
	 * rhs[1..3][1..3] and the solution is returned in u[1..3][1..3]. { void
	 * fill0(double **u, int n); double disc,fact,h=0.5; fill0(u,3);
	 * fact=2.0/(h*h); disc=sqrt(fact*fact+rhs[2][2]); u[2][2] =
	 * -rhs[2][2]/(fact+disc); }
	 * 
	 * void lop(double **out, double **u, int n) Given u[1..n][1..n], returns
	 * Lh(uh) for equation (19.6.44) in out[1..n][1..n]. { int i,j; double
	 * h,h2i; h=1.0/(n-1); h2i=1.0/(h*h); for (j=2;j<n;j++) Interior points. for
	 * (i=2;i<n;i++) out[i][j]=h2i*(u[i+1][j]+u[i-1][j]+u[i][j+1]+u[i][j-1]-
	 * 4.0*u[i][j])+u[i][j]*u[i][j]; for (i=1;i<=n;i++) Boundary points.
	 * out[i][1]=out[i][n]=out[1][i]=out[n][i]=0.0; }
	 * 
	 * void matadd(double **a, double **b, double **c, int n) Adds a[1..n][1..n]
	 * to b[1..n][1..n] and returns result in c[1..n][1..n]. { int i,j; for
	 * (j=1;j<=n;j++) for (i=1;i<=n;i++) c[i][j]=a[i][j]+b[i][j]; }
	 * 
	 * void matsub(double **a, double **b, double **c, int n) Subtracts
	 * b[1..n][1..n] from a[1..n][1..n] and returns result in c[1..n][1..n]. {
	 * int i,j; for (j=1;j<=n;j++) for (i=1;i<=n;i++) c[i][j]=a[i][j]-b[i][j]; }
	 * 
	 * #include <math.h> double anorm2(double **a, int n) Returns the Euclidean
	 * norm of the matrix a[1..n][1..n]. { int i,j; double sum=0.0; for
	 * (j=1;j<=n;j++) for (i=1;i<=n;i++) sum += a[i][j]*a[i][j]; return
	 * sqrt(sum)/n; }
	 */
}
