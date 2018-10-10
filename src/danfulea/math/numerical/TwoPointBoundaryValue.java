package danfulea.math.numerical;

/**
 * TwoPointBoundaryValue Problem.
 * Based on literature, including Numerical Recipes in C (Cambridge Univ.).
 * @author Dan Fulea, 25 OCT. 2006
 */
public class TwoPointBoundaryValue {
	public static Function func;
	public static boolean failB = false;
	public static String failS = "";

	public static double EPS = 1.0e-6;
	public static int nvar = 0;// Variables that you must define and set in your
								// main program.
	public static double x1 = 0.0;
	public static double x2 = 0.0;
	// int kmax,kount; Communicates with odeint.
	// float *xp,**yp,dxsav;
	/*
	 * public static int kmax_odeint=0; public static int kount_odeint=0; public
	 * static double[] xp_odeint; public static double[][] yp_odeint; public
	 * static double dxsav_odeint=0.0;
	 */
	// #define EPS 1.0e-6
	// Variables that you must define and set in your main program.
	public static int nn2_shootf = 0;
	public static int nvar_shootf = 0;
	public static double x1_shootf = 0.0;
	public static double x2_shootf = 0.0;
	public static double xf_shootf = 0.0;

	// int kmax,kount; Communicates with odeint.
	// float *xp,**yp,dxsav;

	/**
	 * Change the sign of 'a' based on value 'b'. If 'b' is greater than 0 then 'a' remains unchanged.
	 * @param a a
	 * @param b b 
	 * @return the result
	 */
	public static double SIGN(double a, double b) {
		return b >= 0.0 ? a : -a;
	}

	/*
	 * When ordinary differential equations are required to satisfy boundary
	 * conditions at more than one value of the independent variable, the
	 * resulting problem is called a two point boundary value problem. As the
	 * terminology indicates, the most common case by far is where boundary
	 * conditions are supposed to be satisfied at two points— usually the
	 * starting and ending values of the integration. However, the phrase “two
	 * point boundary value problem” is also used loosely to include more
	 * complicated cases, e.g., where some conditions are specified at
	 * endpoints, others at interior (usually singular) points. The crucial
	 * distinction between initial value problems (Chapter 16) and two point
	 * boundary value problems (this chapter) is that in the former case we are
	 * able to start an acceptable solution at its beginning (initial values)
	 * and just march it along by numerical integration to its end (final
	 * values); while in the present case, the boundary conditions at the
	 * starting point do not determine a unique solution to start with — and a
	 * “random” choice among the solutions that satisfy these (incomplete)
	 * starting boundary conditions is almost certain not to satisfy the
	 * boundary conditions at the other specified point(s). The “standard” two
	 * point boundary value problem has the following form: We desire the
	 * solution to a set of N coupled first-order ordinary differential
	 * equations, satisfying n1 boundary conditions at the starting point x1,
	 * and a remaining set of n2 = N - n1 boundary conditions at the final point
	 * x2. (Recall that all differential equations of order higher than first
	 * can be written as coupled sets of first-order equations, cf. §16.0.) The
	 * differential equations are dyi(x)/dx= gi(x, y1, y2, . . . , yN) i = 1, 2,
	 * . . .,N At x1, the solution is supposed to satisfy B1j(x1, y1, y2, . . .
	 * , yN) = 0 j = 1, . . . , n1 (17.0.2) while at x2, it is supposed to
	 * satisfy B2k(x2, y1, y2, . . . , yN) = 0 k = 1, . . . , n2 (17.0.3)
	 * 
	 * The Shooting Method In this section we discuss “pure” shooting, where the
	 * integration proceeds from x1 to x2, and we try to match boundary
	 * conditions at the end of the integration. In the next section, we
	 * describe shooting to an intermediate fitting point, where the solution to
	 * the equations and boundary conditions is found by launching “shots” from
	 * both sides of the interval and trying to match continuity conditions at
	 * some intermediate point.
	 * 
	 * Our implementation of the shooting method exactly implements
	 * multidimensional, globally convergent Newton-Raphson (§9.7). It seeks to
	 * zero n 2 functions of n2 variables. The functions are obtained by
	 * integrating N differential equations from x1 to x2. Let us see how this
	 * works: At the starting point x1 there are N starting values yi to be
	 * specified, but subject to n1 conditions. Therefore there are n2 = N -n1
	 * freely specifiable starting values. Let us imagine that these freely
	 * specifiable values are the components of a vector V that lives in a
	 * vector space of dimension n2. Then you, the user, knowing the functional
	 * form of the boundary conditions (17.0.2), can write a function that
	 * generates a complete set of N starting values y, satisfying the boundary
	 * conditions at x1, from an arbitrary vector value of V in which there are
	 * no restrictions on the n 2 component values. In other words, (17.0.2)
	 * converts to a prescription yi(x1) = yi(x1; V1, . . . , Vn2) i = 1, . .
	 * .,N (17.1.1) Below, the function that implements (17.1.1) will be called
	 * load. Notice that the components of V might be exactly the values of
	 * certain “free” components of y, with the other components of y determined
	 * by the boundary conditions. Alternatively, the components of V might
	 * parametrize the solutions that satisfy the starting boundary conditions
	 * in some other convenient way. Boundary conditions often impose algebraic
	 * relations among the y i, rather than specific values for each of them.
	 * Using some auxiliary set of parameters often makes it easier to “solve”
	 * the boundary relations for a consistent set of yi’s. It makes no
	 * difference which way you go, as long as your vector space of V’s
	 * generates (through 17.1.1) all allowed starting vectors y. Given a
	 * particular V, a particular y(x1) is thus generated. It can then be turned
	 * into a y(x2) by integrating the ODEs to x2 as an initial value problem
	 * (e.g., using Chapter 16’s odeint). Now, at x2, let us define a
	 * discrepancy vector F, also of dimension n2, whose components measure how
	 * far we are from satisfying the n2 boundary conditions at x2 (17.0.3).
	 * Simplest of all is just to use the right-hand sides of (17.0.3), Fk =
	 * B2k(x2, y) k = 1, . . . , n2 (17.1.2) As in the case of V, however, you
	 * can use any other convenient parametrization, as long as your space of
	 * F’s spans the space of possible discrepancies from the desired boundary
	 * conditions, with all components of F equal to zero if and only if the
	 * boundary conditions at x2 are satisfied. Below, you will be asked to
	 * supply a user-written function score which uses (17.0.3) to convert an
	 * N-vector of ending values y(x2) into an n2-vector of discrepancies F.
	 * Now, as far as Newton-Raphson is concerned, we are nearly in business. We
	 * want to find a vector value of V that zeros the vector value of F. We do
	 * this by invoking the globally convergent Newton’s method implemented in
	 * the routine newt of §9.7. Recall that the heart of Newton’s method
	 * involves solving the set of n2 linear equations J · delta lower case V =
	 * -F (17.1.3) and then adding the correction back, Vnew = Vold + delta
	 * lower case V (17.1.4) In (17.1.3), the Jacobian matrix J has components
	 * given by Jij = dpartialFi/dpartialVj (17.1.5) It is not feasible to
	 * compute these partial derivatives analytically. Rather, each requires a
	 * separate integration of the N ODEs, followed by the evaluation of
	 * dpartialFi/dpartialVj ~Fi(V1, . . . , Vj +?Vj, . . .) - Fi(V1, . . . ,
	 * Vj, . . .)/DELTA Vj (17.1.6) This is done automatically for you in the
	 * routine fdjac that comes with newt. The only input to newt that you have
	 * to provide is the routine vecfunc that calculates F by integrating the
	 * ODEs. Here is the appropriate routine, called shoot, that is to be passed
	 * as the actual argument in newt: //SEE RootFind line 1491!!
	 */
	/**
	 * Routine for use with newt to solve a two point boundary value problem for nvar coupled ODEs 
	 * by shooting from x1 to x2. Initial values for the nvar ODEs at x1 are generated from the n2 
	 * input coefficients v[1..n2], using the user-supplied routine load. The routine integrates the 
	 * ODEs to x2 using the Runge-Kutta method with tolerance EPS, initial stepsize h1, and minimum 
	 * stepsize hmin. At x2 it calls the user-supplied routine score to evaluate the n2 functions 
	 * f[1..n2] that ought to be zero to satisfy the boundary conditions at x2. The functions f 
	 * are returned on output. newt uses a globally convergent Newton’s method to adjust the values 
	 * of v until the functions f are zero. The user-supplied routine derivs(x,y,dydx) supplies derivative information to the ODE integrator.
	 * The first set of global variables above receives its values from the main program so that shoot can have 
	 * the syntax required for it to be the argument vecfunc of newt.
	 * @param n n
	 * @param v v
	 * @param f f
	 */
	public static void shoot(int n, double[] v, double[] f)
	// Routine for use with newt to solve a two point boundary value problem for
	// nvar coupled ODEs
	// by shooting from x1 to x2. Initial values for the nvar ODEs at x1 are
	// generated from the n2
	// input coefficients v[1..n2], using the user-supplied routine load. The
	// routine integrates the
	// ODEs to x2 using the Runge-Kutta method with tolerance EPS, initial
	// stepsize h1, and minimum
	// stepsize hmin. At x2 it calls the user-supplied routine score to evaluate
	// the n2 functions
	// f[1..n2] that ought to be zero to satisfy the boundary conditions at x2.
	// The functions f
	// are returned on output. newt uses a globally convergent Newton’s method
	// to adjust the values
	// of v until the functions f are zero. The user-supplied routine
	// derivs(x,y,dydx) supplies
	// derivative information to the ODE integrator (see Chapter 16). The first
	// set of global variables
	// above receives its values from the main program so that shoot can have
	// the syntax required
	// for it to be the argument vecfunc of newt.
	{
		// void derivs(float x, float y[], float dydx[]);
		// void load(float x1, float v[], float y[]);
		// void odeint(float ystart[], int nvar, float x1, float x2,
		// float eps, float h1, float hmin, int *nok, int *nbad,
		// void (*derivs)(float, float [], float []),
		// void (*rkqs)(float [], float [], int, float *, float, float,
		// float [], float *, float *, void (*)(float, float [], float [])));
		// void rkqs(float y[], float dydx[], int n, float *x,
		// float htry, float eps, float yscal[], float *hdid, float *hnext,
		// void (*derivs)(float, float [], float []));
		// void score(float xf, float y[], float f[]);

		// int nbad=0;
		// int nok=0;
		double h1 = 0.0;
		double hmin = 0.0;
		double[] y = new double[nvar];
		double[] f1 = new double[f.length];
		// y=vector(1,nvar);
		// kmax=0;//defaulrt and dxsav_odeint is already zero!!!!!!!!!!
		h1 = (x2 - x1) / 100.0;
		// load(x1,v,y);

		func.load(x1, v, y);// y=func.load(x1,v);
		OrdinaryDiffEq.func = func;// @@@@@@@@@@@@@@@@@
		OrdinaryDiffEq.odeint(y, nvar, x1, x2, EPS, h1, hmin);// ,&nok,&nbad,derivs,rkqs);
		// score(x2,y,f);

		// f1=func.score(x2,y);for(int i=1;i<=f.length;i++)f[i-1]=f1[i-1];
		func.score(x2, y, f1);
		for (int i = 1; i <= f.length; i++)
			f[i - 1] = f1[i - 1];

		// free_vector(y,1,nvar);
	}

	/*
	 * Shooting to a Fitting Point
	 * 
	 * The shooting method described in §17.1 tacitly assumed that the “shots”
	 * would be able to traverse the entire domain of integration, even at the
	 * early stages of convergence to a correct solution. In some problems it
	 * can happen that, for very wrong starting conditions, an initial solution
	 * can’t even get from x1 to x2 without encountering some incalculable, or
	 * catastrophic, result. For example, the argument of a square root might go
	 * negative, causing the numerical code to crash. Simple shooting would be
	 * stymied. A different, but related, case is where the endpoints are both
	 * singular points of the set of ODEs. One frequently needs to use special
	 * methods to integrate near the singular points, analytic asymptotic
	 * expansions, for example. In such cases it is feasible to integrate in the
	 * direction away from a singular point, using the special method to get
	 * through the first little bit and then reading off “initial” values for
	 * further numerical integration. However it is usually not feasible to
	 * integrate into a singular point, if only because one has not usually
	 * expended the same analytic effort to obtain expansions of “wrong”
	 * solutions near the singular point (those not satisfying the desired
	 * boundary condition). The solution to the above mentioned difficulties is
	 * shooting to a fitting point. Instead of integrating from x1 to x2, we
	 * integrate first from x1 to some point xf that is between x1 and x2; and
	 * second from x2 (in the opposite direction) to xf If (as before) the
	 * number of boundary conditions imposed at x 1 is n1, and the number
	 * imposed at x2 is n2, then there are n2 freely specifiable starting values
	 * at x1 and n1 freely specifiable starting values at x2. (If you are
	 * confused by this, go back to §17.1.) We can therefore define an n2-vector
	 * V(1) of starting parameters at x1, and a prescription load1(x1,v1,y) for
	 * mapping V(1) into a y that satisfies the boundary conditions at x1,
	 * yi(x1) = yi(x1; V(1)1, . . . , V(1)n2 ) i = 1, . . .,N (17.2.1) Likewise
	 * we can define an n1-vector V(2) of starting parameters at x2, and a
	 * prescription load2(x2,v2,y) for mapping V(2) into a y that satisfies the
	 * boundary conditions at x2, yi(x2) = yi(x2; V(2)1, . . . , V(2)n1 ) i = 1,
	 * . . .,N (17.2.2) We thus have a total of N freely adjustable parameters
	 * in the combination of V(1) and V(2). The N conditions that must be
	 * satisfied are that there be agreement in N components of y at xf between
	 * the values obtained integrating from one side and from the other, yi(xf
	 * ;V(1)) = yi(xf ;V(2)) i = 1, . . .,N (17.2.3) In some problems, the N
	 * matching conditions can be better described (physically, mathematically,
	 * or numerically) by usingN different functions Fi, i = 1. . .N, each
	 * possibly depending on the N components yi. In those cases, (17.2.3) is
	 * replaced by Fi[y(xf ;V(1))] = Fi[y(xf ;V(2))] i = 1, . . .,N
	 * 
	 * In the program below, the user-supplied function score(xf,y,f) is
	 * supposed to map an input N-vector y into an output N-vector F. In most
	 * cases, you can dummy this function as the identity mapping. Shooting to a
	 * fitting point uses globally convergent Newton-Raphson exactly as in
	 * §17.1. Comparing closely with the routine shoot of the previous section,
	 * you should have no difficulty in understanding the following routine
	 * shootf. The main differences in use are that you have to supply both
	 * load1 and load2. Also, in the calling program you must supply initial
	 * guesses for v1[1..n2] and v2[1..n1]. Once again a sample program
	 * illustrating shooting to a fitting point is given in §17.4.
	 */

	/**
	 * Routine for use with newt to solve a two point boundary value problem for nvar coupled 
	 * ODEs by shooting from x1 and x2 to a fitting point xf. Initial values for the nvar ODEs at 
	 * x1 (x2) are generated from the n2 (n1) coefficients v1 (v2), using the user-supplied routine 
	 * load1 (load2). The coefficients v1 and v2 should be stored in a single array v[1..n1+n2] 
	 * in the main program by statements of the form v1=v; and v2 = v[n2];. The input parameter 
	 * n = n1 + n2 = nvar. The routine integrates the ODEs to xf using the Runge-Kutta 
	 * method with tolerance EPS, initial stepsize h1, and minimum stepsize hmin. At xf it calls the 
	 * user-supplied routine score to evaluate the nvar functions f1 and f2 that ought to match 
	 * at xf. The differences f are returned on output. newt uses a globally convergent Newton’s 
	 * method to adjust the values of v until the functions f are zero. The user-supplied routine 
	 * derivs(x,y,dydx) supplies derivative information to the ODE integrator. 
	 * The first set of global variables above receives its values from the main program so that shoot 
	 * can have the syntax required for it to be the argument vecfunc of newt. Set nn2 = n2 in the main program.
	 * @param n n
	 * @param v v
	 * @param f f
	 */
	public static void shootf(int n, double[] v, double[] f)
	// Routine for use with newt to solve a two point boundary value problem for
	// nvar coupled
	// ODEs by shooting from x1 and x2 to a fitting point xf. Initial values for
	// the nvar ODEs at
	// x1 (x2) are generated from the n2 (n1) coefficients v1 (v2), using the
	// user-supplied routine
	// load1 (load2). The coefficients v1 and v2 should be stored in a single
	// array v[1..n1+n2]
	// in the main program by statements of the form v1=v; and v2 = &v[n2];. The
	// input parameter
	// n = n1 + n2 = nvar. The routine integrates the ODEs to xf using the
	// Runge-Kutta
	// method with tolerance EPS, initial stepsize h1, and minimum stepsize
	// hmin. At xf it calls the
	// user-supplied routine score to evaluate the nvar functions f1 and f2 that
	// ought to match
	// at xf. The differences f are returned on output. newt uses a globally
	// convergent Newton’s
	// method to adjust the values of v until the functions f are zero. The
	// user-supplied routine
	// derivs(x,y,dydx) supplies derivative information to the ODE integrator
	// (see Chapter 16).
	// The first set of global variables above receives its values from the main
	// program so that shoot
	// can have the syntax required for it to be the argument vecfunc of newt.
	// Set nn2 = n2 in
	// the main program.
	{
		// void derivs(float x, float y[], float dydx[]);
		// void load1(float x1, float v1[], float y[]);
		// void load2(float x2, float v2[], float y[]);
		// void odeint(float ystart[], int nvar, float x1, float x2,
		// float eps, float h1, float hmin, int *nok, int *nbad,
		// void (*derivs)(float, float [], float []),
		// void (*rkqs)(float [], float [], int, float *, float, float,
		// float [], float *, float *, void (*)(float, float [], float [])));
		// void rkqs(float y[], float dydx[], int n, float *x,
		// float htry, float eps, float yscal[], float *hdid, float *hnext,
		// void (*derivs)(float, float [], float []));
		// void score(float xf, float y[], float f[]);
		int i = 0;// int nbad=0;int nok=0;
		double h1 = 0.0;
		double hmin = 0.0;
		double[] f1 = new double[nvar_shootf];
		double[] f2 = new double[nvar_shootf];
		double[] y = new double[nvar_shootf];

		OrdinaryDiffEq.func = func;
		// f1=vector(1,nvar);
		// f2=vector(1,nvar);
		// y=vector(1,nvar);
		// kmax=0;//default and dxsav_odeint is already zero!!!!!!!!!!
		h1 = (x2_shootf - x1_shootf) / 100.0;
		// load1(x1,v,y); //Path from x1 to xf with best trial values v1.
		func.load1(x1_shootf, v, y);// y=func.load1(x1_shootf,v);
		OrdinaryDiffEq.odeint(y, nvar_shootf, x1_shootf, xf_shootf, EPS, h1,
				hmin);// ,&nok,&nbad,derivs,rkqs);
		// f1=func.score(x2_shootf,y);//score(xf,y,f1);
		func.score(x2_shootf, y, f1);
		// load2(x2,&v[nn2],y); Path from x2 to xf with best trial values v2.
		// y=func.load2(x2_shootf,v,nn2_shootf);//v de la elementul nn2 (-1) la
		// final
		func.load2(x2_shootf, v, nn2_shootf, y);
		OrdinaryDiffEq.odeint(y, nvar_shootf, x2_shootf, xf_shootf, EPS, h1,
				hmin);// ,&nok,&nbad,derivs,rkqs);
		// f2=func.score(xf_shootf,y);//score(xf,y,f2);
		func.score(xf_shootf, y, f2);
		for (i = 1; i <= n; i++)
			f[i - 1] = f1[i - 1] - f2[i - 1];// f[i]=f1[i]-f2[i];
		// free_vector(y,1,nvar);
		// free_vector(f2,1,nvar);
		// free_vector(f1,1,nvar);
	}

	/*
	 * In relaxation methods we replace ODEs by approximate finite-difference
	 * equations (FDEs) on a grid or mesh of points that spans the domain of
	 * interest. As a typical example, we could replace a general first-order
	 * differential equation dy/dx= g(x, y) (17.3.1) with an algebraic equation
	 * relating function values at two points k, k - 1: yk - yk-1 - (xk - xk-1)
	 * g [1/2 (xk + xk-1), 1/2 (yk + yk-1)] = 0 (17.3.2) The form of the FDE in
	 * (17.3.2) illustrates the idea, but not uniquely: There are many ways to
	 * turn the ODE into an FDE. When the problem involves N coupled first-order
	 * ODEs represented by FDEs on a mesh of M points, a solution consists of
	 * values for N dependent functions given at each of the M mesh points, or N
	 * × M variables in all. The relaxation method determines the solution by
	 * starting with a guess and improving it, iteratively. As the iterations
	 * improve the solution, the result is said to relax to the true solution.
	 * .................... The call statement also supplies solvde with the
	 * array y[1..nyj][1..nyk] containing the initial trial solution, and
	 * workspace arrays c[1..ne][1..ne-nb+1][1..m+1], s[1..ne][1..2*ne+1]. The
	 * array c is the blockbuster: It stores the unreduced elements of the
	 * matrix built up for the backsubstitution step. If there are m mesh
	 * points, then there will be m+1 blocks, each requiring ne rows and ne-nb+1
	 * columns. Although large, this is small compared with (ne×m)2 elements
	 * required for the whole matrix if we did not break it into blocks. We now
	 * describe the workings of the user-supplied function difeq. The synopsis
	 * of the function is void difeq(int k, int k1, int k2, int jsf, int is1,
	 * int isf, int indexv[], int ne, float **s, float **y); The only
	 * information passed from difeq to solvde is the matrix of derivatives
	 * s[1..ne][1..2*ne+1]; all other arguments are input to difeq and should
	 * not be altered. k indicates the current mesh point, or block number.
	 * k1,k2 label the first and last point in the mesh. If k=k1 or k>k2, the
	 * block involves the boundary conditions at the first or final points;
	 * otherwise the block acts on FDEs coupling variables at points k-1, k. The
	 * convention on storing information into the array s[i][j] follows that
	 * used in equations (17.3.8), (17.3.10), and (17.3.12): Rows i label
	 * equations, columns j refer to derivatives with respect to dependent
	 * variables in the solution. Recall that each equation will depend on the
	 * ne dependent variables at either one or two points. Thus, j runs from 1
	 * to either ne or 2*ne. The column ordering for dependent variables at each
	 * point must agree with the list supplied in indexv[j]. Thus, for a block
	 * not at a boundary, the first column multiplies?Y (l=indexv[1],k-1), and
	 * the column ne+1 multiplies?Y (l=indexv[1],k). is1,isf give the numbers of
	 * the starting and final rows that need to be filled in the s matrix for
	 * this block. jsf labels the column in which the difference equations Ej,k
	 * of equations (17.3.3)–(17.3.5) are stored. Thus, -s[i][jsf] is the vector
	 * on the right-hand side of the matrix. The reason for the minus sign is
	 * that difeq supplies the actual difference equation, Ej,k, not its
	 * negative. Note that solvde supplies a value for jsf such that the
	 * difference equation is put in the column just after all derivatives in
	 * the s matrix. Thus, difeq expects to find values entered into s[i][j] for
	 * rows is1 ? i ? isf and 1 ? j ? jsf. Finally, s[1..nsi][1..nsj] and
	 * y[1..nyj][1..nyk] supply difeq with storage for s and the solution
	 * variables y for this iteration. An example of how to use this routine is
	 * given in the next section
	 */

	/**
	 * Driver routine for solution of two point boundary value problems by relaxation. itmax is the 
	 * maximum number of iterations. conv is the convergence criterion. slowc controls the fraction of corrections actually used after each iteration. 
	 * scalv[1..ne] contains typical sizes for each dependent variable, used to weight errors. indexv[1..ne] lists the column 
	 * ordering of variables used to construct the matrix s[1..ne][1..2*ne+1] of derivatives. (The nb boundary conditions at the first mesh point must contain some 
	 * dependence on the first nb variables listed in indexv.) The problem involves ne equations for ne adjustable dependent 
	 * variables at each point. At the first mesh point there are nb boundary conditions. There are a 
	 * total of m mesh points. y[1..ne][1..m] is the two-dimensional array that contains the initial guess for all the dependent variables at each mesh point. On each 
	 * iteration, it is updated by the calculated correction. The arrays c[1..ne][1..ne-nb+1][1..m+1] and s supply dummy storage used by the relaxation code.
	 * @param itmax itmax
	 * @param conv conv
	 * @param slowc slowc
	 * @param scalv scalv
	 * @param indexv indexv
	 * @param ne ne
	 * @param nb nb
	 * @param m m
	 * @param y y
	 * @param c c
	 * @param s s
	 */
	public static void solvde(int itmax, double conv, double slowc,
			double[] scalv, int[] indexv, int ne, int nb, int m, double[][] y,
			double[][][] c, double[][] s)
	// Driver routine for solution of two point boundary value problems by
	// relaxation. itmax is the
	// maximum number of iterations. conv is the convergence criterion (see
	// text). slowc controls
	// the fraction of corrections actually used after each iteration.
	// scalv[1..ne] contains typical
	// sizes for each dependent variable, used to weight errors. indexv[1..ne]
	// lists the column
	// ordering of variables used to construct the matrix s[1..ne][1..2*ne+1] of
	// derivatives. (The
	// nb boundary conditions at the first mesh point must contain some
	// dependence on the first nb
	// variables listed in indexv.) The problem involves ne equations for ne
	// adjustable dependent
	// variables at each point. At the first mesh point there are nb boundary
	// conditions. There are a
	// total of m mesh points. y[1..ne][1..m] is the two-dimensional array that
	// contains the initial
	// guess for all the dependent variables at each mesh point. On each
	// iteration, it is updated by the
	// calculated correction. The arrays c[1..ne][1..ne-nb+1][1..m+1] and s
	// supply dummy
	// storage used by the relaxation code.
	{
		// void bksub(int ne, int nb, int jf, int k1, int k2, float ***c);
		// void difeq(int k, int k1, int k2, int jsf, int is1, int isf,
		// int indexv[], int ne, float **s, float **y);
		// void pinvs(int ie1, int ie2, int je1, int jsf, int jc1, int k,
		// float ***c, float **s);
		// void red(int iz1, int iz2, int jz1, int jz2, int jm1, int jm2, int
		// jmf,
		// int ic1, int jc1, int jcf, int kc, float ***c, float **s);
		// s[1..ne][1..2*ne+1];
		failB = false;
		double[][] ss = new double[ne][2 * ne + 1];
		for (int ii = 1; ii <= ne; ii++)
			for (int jj = 1; jj <= 2 * ne + 1; jj++)
				ss[ii - 1][jj - 1] = s[ii - 1][jj - 1];

		int ic1 = 0;
		int ic2 = 0;
		int ic3 = 0;
		int ic4 = 0;
		int it = 0;
		int j = 0;
		int j1 = 0;
		int j2 = 0;
		int j3 = 0;
		int j4 = 0;
		int j5 = 0;
		int j6 = 0;
		int j7 = 0;
		int j8 = 0;
		int j9 = 0;
		int jc1 = 0;
		int jcf = 0;
		int jv = 0;
		int k = 0;
		int k1 = 0;
		int k2 = 0;
		int km = 0;
		int kp = 0;
		int nvars = 0;
		int[] kmax = new int[ne];
		double err = 0.0;
		double errj = 0.0;
		double fac = 0.0;
		double vmax = 0.0;
		double vz = 0.0;
		double[] ermax = new double[ne];
		// kmax=ivector(1,ne);
		// ermax=vector(1,ne);
		k1 = 1; // Set up row and column markers.
		k2 = m;
		nvars = ne * m;
		j1 = 1;
		j2 = nb;
		j3 = nb + 1;
		j4 = ne;
		j5 = j4 + j1;
		j6 = j4 + j2;
		j7 = j4 + j3;
		j8 = j4 + j4;
		j9 = j8 + j1;
		ic1 = 1;
		ic2 = ne - nb;
		ic3 = ic2 + 1;
		ic4 = ne;
		jc1 = 1;
		jcf = ic3;
		for (it = 1; it <= itmax; it++) {// Primary iteration loop.
			k = k1; // Boundary conditions at first point.
			// func.difeq(k,k1,k2,j9,ic3,ic4,indexv,ne,s,y);
			func.difeq(k, k1, k2, j9, ic3, ic4, indexv, ne, ss, y);
			// pinvs(ic3,ic4,j5,j9,jc1,k1,c,s);
			pinvs(ic3, ic4, j5, j9, jc1, k1, c, ss);
			if (failB)
				return;
			for (k = k1 + 1; k <= k2; k++) {// Finite difference equations at
											// all point pairs.
				kp = k - 1;
				// difeq(k,k1,k2,j9,ic1,ic4,indexv,ne,s,y);
				func.difeq(k, k1, k2, j9, ic1, ic4, indexv, ne, ss, y);
				// red(ic1,ic4,j1,j2,j3,j4,j9,ic3,jc1,jcf,kp,c,s);
				red(ic1, ic4, j1, j2, j3, j4, j9, ic3, jc1, jcf, kp, c, ss);
				// pinvs(ic1,ic4,j3,j9,jc1,k,c,s);
				pinvs(ic1, ic4, j3, j9, jc1, k, c, ss);
				if (failB)
					return;
			}
			k = k2 + 1; // Final boundary conditions.
			// difeq(k,k1,k2,j9,ic1,ic2,indexv,ne,s,y);
			func.difeq(k, k1, k2, j9, ic1, ic2, indexv, ne, ss, y);
			// red(ic1,ic2,j5,j6,j7,j8,j9,ic3,jc1,jcf,k2,c,s);
			red(ic1, ic2, j5, j6, j7, j8, j9, ic3, jc1, jcf, k2, c, ss);
			// pinvs(ic1,ic2,j7,j9,jcf,k2+1,c,s);
			pinvs(ic1, ic2, j7, j9, jcf, k2 + 1, c, ss);
			if (failB)
				return;
			bksub(ne, nb, jcf, k1, k2, c); // Backsubstitution.
			err = 0.0;
			for (j = 1; j <= ne; j++) {// Convergence check, accumulate average
										// error.
				jv = indexv[j - 1];// jv=indexv[j];
				errj = vmax = 0.0;
				km = 0;
				for (k = k1; k <= k2; k++) {// Find point with largest error,
											// for each dependent variable.
					vz = Math.abs(c[jv - 1][0][k - 1]);// vz=Math.abs(c[jv][1][k]);
					if (vz > vmax) {
						vmax = vz;
						km = k;
					}
					errj += vz;
				}
				// err += errj/scalv[j]; Note weighting for each dependent
				// variable.
				err += errj / scalv[j - 1];
				// ermax[j]=c[jv][1][km]/scalv[j];
				ermax[j - 1] = c[jv - 1][0][km - 1] / scalv[j - 1];
				kmax[j - 1] = km;// kmax[j]=km;
			}
			err /= nvars;
			fac = (err > slowc ? slowc / err : 1.0);
			// Reduce correction applied when error is large.
			for (j = 1; j <= ne; j++) {// Apply corrections.
				jv = indexv[j - 1];// jv=indexv[j];
				for (k = k1; k <= k2; k++)
					// y[j][k] -= fac*c[jv][1][k];
					y[j - 1][k - 1] -= fac * c[jv - 1][0][k - 1];
			}
			// printf("\n%8s %9s %9s\n","Iter.","Error","FAC");
			// Summary of corrections for this step.
			// printf("%6d %12.6f %11.6f\n",it,err,fac);
			String str = "Iter: " + it + "; Error: " + err + "; FAC: " + fac;
			func.printSequence(str);
			if (err < conv) {// Point with largest error for each variable can
								// be monitored by writing out kmax and ermax.
								// free_vector(ermax,1,ne);
								// free_ivector(kmax,1,ne);
				return;
			}
		}
		// nrerror("Too many iterations in solvde"); Convergence failed.
		failB = true;
		failS = "Too many iterations in solvde";
		return;
	}

	/**
	 * Backsubstitution, used internally by solvde.
	 * @param ne ne
	 * @param nb nb
	 * @param jf jf
	 * @param k1 k1
	 * @param k2 k2
	 * @param c c
	 */
	public static void bksub(int ne, int nb, int jf, int k1, int k2,
			double[][][] c)
	// Backsubstitution, used internally by solvde.
	{
		int nbf = 0;
		int im = 0;
		int kp = 0;
		int k = 0;
		int j = 0;
		int i = 0;
		double xx = 0.0;

		nbf = ne - nb;
		im = 1;
		for (k = k2; k >= k1; k--) {// Use recurrence relations to eliminate
									// remaining dependences.
			if (k == k1)
				im = nbf + 1;
			kp = k + 1;
			for (j = 1; j <= nbf; j++) {
				xx = c[j - 1][jf - 1][kp - 1];// xx=c[j][jf][kp];
				for (i = im; i <= ne; i++)
					// c[i][jf][k] -= c[i][j][k]*xx;
					c[i - 1][jf - 1][k - 1] -= c[i - 1][j - 1][k - 1] * xx;
			}
		}
		for (k = k1; k <= k2; k++) {// Reorder corrections to be in column 1.
			kp = k + 1;
			for (i = 1; i <= nb; i++)
				// c[i][1][k]=c[i+nbf][jf][k];
				c[i - 1][0][k - 1] = c[i + nbf - 1][jf - 1][k - 1];
			for (i = 1; i <= nbf; i++)
				// c[i+nb][1][k]=c[i][jf][kp];
				c[i + nb - 1][0][k - 1] = c[i - 1][jf - 1][kp - 1];
		}
	}

	/**
	 * Diagonalize the square subsection of the s matrix, and store the recursion coefficients in c; used internally by solvde.
	 * @param ie1 ie1
	 * @param ie2 ie2
	 * @param je1 je1
	 * @param jsf jsf
	 * @param jc1 jc1
	 * @param k k
	 * @param c c
	 * @param s s
	 */
	public static void pinvs(int ie1, int ie2, int je1, int jsf, int jc1,
			int k, double[][][] c, double[][] s)
	// Diagonalize the square subsection of the s matrix, and store the
	// recursion coefficients in c;
	// used internally by solvde.
	{
		int js1 = 0;
		int jpiv = 0;
		int jp = 0;
		int je2 = 0;
		int jcoff = 0;
		int j = 0;
		int irow = 0;
		int ipiv = 0;
		int id = 0;
		int icoff = 0;
		int i = 0;
		int[] indxr = new int[ie2];// new int[ie2-ie1+1];
		double pivinv = 0.0;
		double piv = 0.0;
		double dum = 0.0;
		double big = 0.0;
		double[] pscl = new double[ie2];// new double[ie2-ie1+1];
		// indxr=ivector(ie1,ie2);
		// pscl=vector(ie1,ie2);
		failB = false;
		je2 = je1 + ie2 - ie1;
		js1 = je2 + 1;
		for (i = ie1; i <= ie2; i++) {// Implicit pivoting, as in §2.1.
			big = 0.0;
			for (j = je1; j <= je2; j++)
				// if (Math.abs(s[i][j]) > big) big=Math.abs(s[i][j]);
				if (Math.abs(s[i - 1][j - 1]) > big)
					big = Math.abs(s[i - 1][j - 1]);
			if (big == 0.0) {
				failB = true;
				failS = "Singular matrix - row all 0, in pinvs";
				return;
				// nrerror("Singular matrix - row all 0, in pinvs");
			}
			pscl[i - 1] = 1.0 / big;// pscl[i]=1.0/big;
			indxr[i - 1] = 0;// indxr[i]=0;
		}
		for (id = ie1; id <= ie2; id++) {
			piv = 0.0;
			for (i = ie1; i <= ie2; i++) {// Find pivot element.
				if (indxr[i - 1] == 0) // if (indxr[i] == 0)
				{
					big = 0.0;
					for (j = je1; j <= je2; j++) {
						// if (fabs(s[i][j]) > big)
						if (Math.abs(s[i - 1][j - 1]) > big) {
							jp = j;
							// big=fabs(s[i][j]);
							big = Math.abs(s[i - 1][j - 1]);
						}
					}
					// if (big*pscl[i] > piv)
					if (big * pscl[i - 1] > piv) {
						ipiv = i;
						jpiv = jp;
						// piv=big*pscl[i];
						piv = big * pscl[i - 1];
					}
				}
			}
			// if (s[ipiv][jpiv] == 0.0)
			if (s[ipiv - 1][jpiv - 1] == 0.0) {
				failB = true;
				failS = "Singular matrix in routine pinvs";
				return;
				// nrerror("Singular matrix in routine pinvs");
			}
			indxr[ipiv - 1] = jpiv;// indxr[ipiv]=jpiv; In place reduction. Save
									// column ordering.
			pivinv = 1.0 / s[ipiv - 1][jpiv - 1];// pivinv=1.0/s[ipiv][jpiv];
			// for (j=je1;j<=jsf;j++) s[ipiv][j] *= pivinv; Normalize pivot row.
			for (j = je1; j <= jsf; j++)
				s[ipiv - 1][j - 1] *= pivinv;
			s[ipiv - 1][jpiv - 1] = 1.0;// s[ipiv][jpiv]=1.0;
			for (i = ie1; i <= ie2; i++) {// Reduce nonpivot elements in column.
											// if (indxr[i] != jpiv)
				if (indxr[i - 1] != jpiv) {
					if (s[i - 1][jpiv - 1] != 0.0) // if (s[i][jpiv])
					{
						dum = s[i - 1][jpiv - 1];// dum=s[i][jpiv];
						for (j = je1; j <= jsf; j++)
							// s[i][j] -= dum*s[ipiv][j];
							s[i - 1][j - 1] -= dum * s[ipiv - 1][j - 1];
						s[i - 1][jpiv - 1] = 0.0;// s[i][jpiv]=0.0;
					}
				}
			}
		}
		jcoff = jc1 - js1;// Sort and store unreduced coefficients.
		icoff = ie1 - je1;
		for (i = ie1; i <= ie2; i++) {
			irow = indxr[i - 1] + icoff;// irow=indxr[i]+icoff;
			for (j = js1; j <= jsf; j++)
				// c[irow][j+jcoff][k]=s[i][j];
				c[irow - 1][j + jcoff - 1][k - 1] = s[i - 1][j - 1];
		}
		// free_vector(pscl,ie1,ie2);
		// free_ivector(indxr,ie1,ie2);
	}

	/**
	 * Reduce columns jz1-jz2 of the s matrix, using previous results as stored in the c matrix. Only 
	 * columns jm1-jm2,jmf are affected by the prior results. red is used internally by solvde.
	 * @param iz1 iz1
	 * @param iz2 iz2
	 * @param jz1 jz1
	 * @param jz2 jz2
	 * @param jm1 jm1
	 * @param jm2 jm2
	 * @param jmf jmf
	 * @param ic1 ic1
	 * @param jc1 jc1
	 * @param jcf jcf
	 * @param kc kc
	 * @param c c
	 * @param s s
	 */
	public static void red(int iz1, int iz2, int jz1, int jz2, int jm1,
			int jm2, int jmf, int ic1, int jc1, int jcf, int kc,
			double[][][] c, double[][] s)
	// Reduce columns jz1-jz2 of the s matrix, using previous results as stored
	// in the c matrix. Only
	// columns jm1-jm2,jmf are affected by the prior results. red is used
	// internally by solvde.
	{
		int loff = 0;
		int l = 0;
		int j = 0;
		int ic = 0;
		int i = 0;
		double vx = 0.0;

		loff = jc1 - jm1;
		ic = ic1;
		for (j = jz1; j <= jz2; j++) {// Loop over columns to be zeroed.
			for (l = jm1; l <= jm2; l++) {// Loop over columns altered.
											// vx=c[ic][l+loff][kc];
				vx = c[ic - 1][l + loff - 1][kc - 1];
				for (i = iz1; i <= iz2; i++)
					// s[i][l] -= s[i][j]*vx; Loop over rows.
					s[i - 1][l - 1] -= s[i - 1][j - 1] * vx;
			}
			// vx=c[ic][jcf][kc];
			vx = c[ic - 1][jcf - 1][kc - 1];
			for (i = iz1; i <= iz2; i++)
				// s[i][jmf] -= s[i][j]*vx; Plus final element.
				s[i - 1][jmf - 1] -= s[i - 1][j - 1] * vx;
			ic += 1;
		}
	}

	/*
	 * The best way to understand the algorithms of the previous sections is to
	 * see them employed to solve an actual problem. As a sample problem, we
	 * have selected the computation of spheroidal harmonics. (The more common
	 * name is spheroidal angle functions, but we prefer the explicit reminder
	 * of the kinship with spherical harmonics.) We will show how to find
	 * spheroidal harmonics, first by the method of relaxation (§17.3), and then
	 * by the methods of shooting (§17.1) and shooting to a fitting point
	 * (§17.2). Spheroidal harmonics typically arise when certain partial
	 * differential equations are solved by separation of variables in
	 * spheroidal coordinates. They satisfy the following differential equation
	 * on the interval -1 <= x <= 1: d/dx [(1 - x2)dS/dx ]+ (lambda - c^2x^2
	 * -m2/(1 - x2))S = 0
	 * 
	 * Here m is an integer, c is the “oblateness parameter,” and ? is the
	 * eigenvalue. Despite the notation, c2 can be positive or negative. For c2
	 * > 0 the functions are called “prolate,” while if c2 < 0 they are called
	 * “oblate.” The equation has singular points at x = ±1 and is to be solved
	 * subject to the boundary conditions that the solution be regular at x =
	 * ±1. Only for certain values of ?, the eigenvalues, will this be possible.
	 * If we consider first the spherical case, where c = 0, we recognize the
	 * differential equation for Legendre functions P m n (x). In this case the
	 * eigenvalues are ?mn = n(n + 1), n = m,m + 1, . . . . The integer n labels
	 * successive eigenvalues for fixed m: When n = m we have the lowest
	 * eigenvalue, and the corresponding eigenfunction has no nodes in the
	 * interval -1 < x < 1; when n = m + 1 we have the next eigenvalue, and the
	 * eigenfunction has one node inside (-1, 1); and so on. A similar situation
	 * holds for the general case c2 = 0. We write the eigenvalues of (17.4.1)
	 * as ?mn(c) and the eigenfunctions as Smn(x; c). For fixed m, n = m,m + 1,
	 * . . . labels the successive eigenvalues. The computation of ?mn(c) and
	 * Smn(x; c) traditionally has been quite difficult. Complicated recurrence
	 * relations, power series expansions, etc., can be found in references
	 * [1-3]. Cheap computing makes evaluation by direct solution of the
	 * differential equation quite feasible. The first step is to investigate
	 * the behavior of the solution near the singular points x = ±1.
	 * Substituting a power series expansion of the form S = (1 ± x)^alpha SUM
	 * k=0 to +INF ak(1 ± x)^k (17.4.2) in equation (17.4.1), we find that the
	 * regular solution has alpha = m/2. (Without loss of generality we can take
	 * m >= 0 since m -› -m is a symmetry of the equation.) We get an equation
	 * that is numerically more tractable if we factor out this behavior.
	 * Accordingly we set S = (1 - x^2)^m/2 * y (17.4.3) We then find from
	 * (17.4.1) that y satisfies the equation (1 - x2) d2y/dx2 - 2(m+ 1)x dy/dx
	 * + (µ - c^2x^2)y = 0 (17.4.4) where µ = lambda - m(m + 1) (17.4.5) Both
	 * equations (17.4.1) and (17.4.4) are invariant under the replacement x-›
	 * -x. Thus the functions S and y must also be invariant, except possibly
	 * for an overall scale factor. (Since the equations are linear, a constant
	 * multiple of a solution is also a solution.) Because the solutions will be
	 * normalized, the scale factor can only be ±1. If n-m is odd, there are an
	 * odd number of zeros in the interval (-1, 1). Thus we must choose the
	 * antisymmetric solution y(-x) = -y(x) which has a zero at x = 0.
	 * Conversely, if n - m is even we must have the symmetric solution. Thus
	 * ymn(-x) = (-1)^(n-m) ymn(x) (17.4.6) and similarly for Smn.
	 * 
	 * The boundary conditions on (17.4.4) require that y be regular at x = ±1.
	 * In other words, near the endpoints the solution takes the form y = a0 +
	 * a1(1 - x^2) + a2(1 - x^2)^2 + . . . (17.4.7) Substituting this expansion
	 * in equation (17.4.4) and letting x -› 1, we find that a1 = - [µ - c2]/4(m
	 * + 1) a0 (17.4.8) Equivalently, y'(1) =[ µ - c2]/2(m+ 1) y(1) (17.4.9) A
	 * similar equation holds at x = -1 with a minus sign on the right-hand
	 * side. The irregular solution has a different relation between function
	 * and derivative at the endpoints. Instead of integrating the equation from
	 * -1 to 1, we can exploit the symmetry (17.4.6) to integrate from 0 to 1.
	 * The boundary condition at x = 0 is y(0) = 0, n- m odd y'(0) = 0, n- m
	 * even (17.4.10) A third boundary condition comes from the fact that any
	 * constant multiple of a solution y is a solution. We can thus normalize
	 * the solution. We adopt the normalization that the function Smn has the
	 * same limiting behavior as P m n at x = 1: lim ;x›1 of (1 -
	 * x^2)^(-m/2)Smn(x; c) = lim; x›1 of (1 - x^2)^(-m/2)Pmn (x) (17.4.11)
	 * Various normalization conventions in the literature are tabulated by
	 * Flammer [1]. Imposing three boundary conditions for the second-order
	 * equation (17.4.4) turns it into an eigenvalue problem for ? or
	 * equivalently for µ. We write it in the standard form by setting y1 = y
	 * (17.4.12) y2 = y' (17.4.13) y3 = µ (17.4.14) Then y'1 = y2 (17.4.15) y'2
	 * =1/(1 - x2) [2x(m + 1)y2 - (y3 - c^2x^2)y1] (17.4.16) y'3 = 0
	 * 
	 * The boundary condition at x = 0 in this notation is y1 = 0, n- m odd y2 =
	 * 0, n- m even (17.4.18) At x = 1 we have two conditions: y2 = [y3 -
	 * c2]/2(m+ 1) y1 (17.4.19) y1 = lim; x›1 of (1 - x2)^(-m/2)Pmn (x) =
	 * (-1)m(n + m)!/[2^mm!(n - m)!] =gamma (17.4.20) We are now ready to
	 * illustrate the use of the methods of previous sections on this problem.
	 */
	public static int check_ln = 0;
	public static double f_ln = 0.0;
	public static double EPS2 = 1.0e-4;
	public static double ALF = 1.0e-4;
	public static int MAXITS = 200;
	public static double TOLF = 1.0e-4;
	public static double TOLMIN = 1.0e-6;
	public static double STPMX = 100.0;
	public static double TOLX = 1.0e-7;
	public static int nn_newt = 0;
	public static double[] fvec_newt;// =1.0e-4;
	// public void newt(float x[], int n, int *check,void (*vecfunc)(int, float
	// [], float []))

	/**
	 * Used internally. Given an initial guess x[1..n] for a root in n dimensions, find the root 
	 * by a globally convergent Newton’s method. The vector of functions to be zeroed, called fvec_newt
	 * in the routine below, is returned by the user-supplied routine. The 
	 * output quantity check is false (0) on a normal return and true (1) if the routine has 
	 * converged to a local minimum of the function fmin defined below. In this case try restarting 
	 * from a different initial guess.
	 * @param x x
	 * @param n n
	 * @param nameS nameS
	 */
	public static void newt(double[] x, int n, String nameS)
	// Given an initial guess x[1..n] for a root in n dimensions, find the root
	// by a globally convergent
	// Newton’s method. The vector of functions to be zeroed, called fvec[1..n]
	// in the routine
	// below, is returned by the user-supplied routine vecfunc(n,x,fvec). The
	// output quantity
	// check is false (0) on a normal return and true (1) if the routine has
	// converged to a local
	// minimum of the function fmin defined below. In this case try restarting
	// from a different initial
	// guess.
	{
		// void fdjac(int n, float x[], float fvec[], float **df,void
		// (*vecfunc)(int, float [], float []));
		// float fmin(float x[]);
		// void lnsrch(int n, float xold[], float fold, float g[], float p[],
		// float x[],float *f, float stpmax, int *check, float (*func)(float
		// []));
		// void lubksb(float **a, int n, int *indx, float b[]);
		// void ludcmp(float **a, int n, int *indx, float *d);

		int i = 0;
		int its = 0;
		int j = 0;
		int[] indx = new int[n];
		// double d=0.0;
		double den = 0.0;
		double f = 0.0;
		double fold = 0.0;
		double stpmax = 0.0;
		double sum = 0.0;
		double temp = 0.0;
		double test = 0.0;
		double[][] fjac = new double[n][n];
		double[] g = new double[n];
		double[] p = new double[n];
		double[] xold = new double[n];

		failB = false;
		// indx=ivector(1,n);
		// fjac=matrix(1,n,1,n);
		// g=vector(1,n);
		// p=vector(1,n);
		// xold=vector(1,n);
		// fvec=vector(1,n); Define global variables.

		nn_newt = n;
		fvec_newt = new double[nn_newt];
		// nrfuncv=vecfunc;
		f = fmin(x, nameS); // fvec is also computed by this call.
		test = 0.0; // Test for initial guess being a root. Use more stringent
					// test than simply TOLF.
		for (i = 1; i <= n; i++)
			// if (Math.abs(fvec[i]) > test) test=fabs(fvec[i]);
			if (Math.abs(fvec_newt[i - 1]) > test)
				test = Math.abs(fvec_newt[i - 1]);
		if (test < 0.01 * TOLF) {
			check_ln = 0;
			return;// FREERETURN
		}
		for (sum = 0.0, i = 1; i <= n; i++)
			// sum += SQR(x[i]); Calculate stpmax for line searches.
			sum += x[i - 1] * x[i - 1];
		stpmax = STPMX * Math.max(Math.sqrt(sum), (double) n);
		for (its = 1; its <= MAXITS; its++) {// Start of iteration loop.
			fjac = fdjac(n, x, fvec_newt, nameS);// fdjac(n,x,fvec,fjac,vecfunc);
			// If analytic Jacobian is available, you can replace the routine
			// fdjac below with your
			// own routine.
			for (i = 1; i <= n; i++) {// Compute gradientf for the line search.
				for (sum = 0.0, j = 1; j <= n; j++)
					sum += fjac[j - 1][i - 1] * fvec_newt[j - 1];// sum +=
																	// fjac[j][i]*fvec[j];
				g[i - 1] = sum;// g[i]=sum;
			}
			for (i = 1; i <= n; i++)
				xold[i - 1] = x[i - 1];// xold[i]=x[i]; Store x,
			fold = f;// and f.
			for (i = 1; i <= n; i++)
				p[i - 1] = -fvec_newt[i - 1];// p[i] = -fvec[i]; Right-hand side
												// for linear equations.
			LinAEq.ludcmp(fjac, n, indx);// ,&d);// Solve linear equations by LU
											// decomposition.
			LinAEq.lubksb(fjac, n, indx, p);
			// lnsrch(n,xold,fold,g,p,x,&f,stpmax,check,fmin);
			lnsrch(n, xold, fold, g, p, x, stpmax, "fmin", nameS);
			f = f_ln;
			// lnsrch returns new x and f. It also calculates fvec at the new x
			// when it calls fmin.
			test = 0.0; // Test for convergence on function values.
			for (i = 1; i <= n; i++)
				// if (fabs(fvec[i]) > test) test=fabs(fvec[i]);
				if (Math.abs(fvec_newt[i - 1]) > test)
					test = Math.abs(fvec_newt[i - 1]);
			if (test < TOLF) {
				check_ln = 0;
				return;// FREERETURN
			}
			if (check_ln != 0)// (*check)
			{// Check for gradient of f zero, i.e., spurious convergence.
				test = 0.0;
				den = Math.max(f, 0.5 * n);
				for (i = 1; i <= n; i++) {
					// temp=fabs(g[i])*FMAX(fabs(x[i]),1.0)/den;
					temp = Math.abs(g[i - 1])
							* Math.max(Math.abs(x[i - 1]), 1.0) / den;
					if (temp > test)
						test = temp;
				}
				check_ln = (test < TOLMIN ? 1 : 0);
				return;// FREERETURN
			}
			test = 0.0;// Test for convergence on dx.
			for (i = 1; i <= n; i++) {
				// temp=(fabs(x[i]-xold[i]))/FMAX(fabs(x[i]),1.0);
				temp = (Math.abs(x[i - 1] - xold[i - 1]))
						/ Math.max(Math.abs(x[i - 1]), 1.0);
				if (temp > test)
					test = temp;
			}
			if (test < TOLX)
				return;// FREERETURN
		}
		// nrerror("MAXITS exceeded in newt");
		failB = true;
		failS = "MAXITS exceeded in newt";
	}

	/**
	 * Used internally. Given an n-dimensional point xold[1..n], the value of the function and gradient there, fold 
	 * and g[1..n], and a direction p[1..n], finds a new point x[1..n] along the direction p from 
	 * xold where the function func has decreased “sufficiently.” The new function value is returned 
	 * in f. stpmax is an input quantity that limits the length of the steps so that you do not try to 
	 * evaluate the function in regions where it is undefined or subject to overflow. p is usually the 
	 * Newton direction. The output quantity check is false (0) on a normal exit. It is true (1) when 
	 * x is too close to xold. In a minimization algorithm, this usually signals convergence and can 
	 * be ignored. However, in a zero-finding algorithm the calling program should check whether the 
	 * convergence is spurious.
	 * @param n n
	 * @param xold xold
	 * @param fold fold
	 * @param g g
	 * @param p p
	 * @param x x
	 * @param stpmax stpmax
	 * @param namS namS
	 * @param nameS nameS
	 */
	public static void lnsrch(int n, double[] xold, double fold, double[] g,
			double[] p, double[] x, double stpmax, String namS, String nameS)// ,
																				// int
																				// *check,
																				// float
																				// (*func)(float
																				// []))
	// Given an n-dimensional point xold[1..n], the value of the function and
	// gradient there, fold
	// and g[1..n], and a direction p[1..n], finds a new point x[1..n] along the
	// direction p from
	// xold where the function func has decreased “sufficiently.” The new
	// function value is returned
	// in f. stpmax is an input quantity that limits the length of the steps so
	// that you do not try to
	// evaluate the function in regions where it is undefined or subject to
	// overflow. p is usually the
	// Newton direction. The output quantity check is false (0) on a normal
	// exit. It is true (1) when
	// x is too close to xold. In a minimization algorithm, this usually signals
	// convergence and can
	// be ignored. However, in a zero-finding algorithm the calling program
	// should check whether the
	// convergence is spurious. Some “difficult” problems may require double
	// precision in this routine.
	{
		int i = 0;
		double a = 0.0;
		double alam = 0.0;
		double alam2 = 0.0;
		double alamin = 0.0;
		double b = 0.0;
		double disc = 0.0;
		double f2 = 0.0;
		double rhs1 = 0.0;
		double rhs2 = 0.0;
		double slope = 0.0;
		double sum = 0.0;
		double temp = 0.0;
		double test = 0.0;
		double tmplam = 0.0;

		failB = false;

		check_ln = 0;
		for (sum = 0.0, i = 1; i <= n; i++)
			sum += p[i - 1] * p[i - 1];// sum += p[i]*p[i];
		sum = Math.sqrt(sum);
		if (sum > stpmax)
			for (i = 1; i <= n; i++)
				p[i - 1] *= stpmax / sum;// p[i] *= stpmax/sum; Scale if
											// attempted step is too big.
		for (slope = 0.0, i = 1; i <= n; i++)
			slope += g[i - 1] * p[i - 1];// slope += g[i]*p[i];
		if (slope >= 0.0) {
			// nrerror("Roundoff problem in lnsrch.");
			failB = true;
			failS = "Roundoff problem in lnsrch.";
			return;
		}
		test = 0.0; // Compute ?min.
		for (i = 1; i <= n; i++) {
			temp = Math.abs(p[i - 1]) / Math.max(Math.abs(xold[i - 1]), 1.0);// fabs(p[i])/FMAX(fabs(xold[i]),1.0);
			if (temp > test)
				test = temp;
		}
		alamin = TOLX / test;
		alam = 1.0;// Always try full Newton step first.
		for (;;) {// Start of iteration loop.
			for (i = 1; i <= n; i++)
				x[i - 1] = xold[i - 1] + alam * p[i - 1];// x[i]=xold[i]+alam*p[i];
			f_ln = localFunc(x, namS, nameS);// func.vecfunc(n,x);//*f=(*func)(x);
			if (alam < alamin) {// Convergence on ?x. For zero finding,the
								// calling program should verify the
								// convergence.
				for (i = 1; i <= n; i++)
					x[i - 1] = xold[i - 1];// x[i]=xold[i];
				check_ln = 1;
				return;
			} else if (f_ln <= fold + ALF * alam * slope)
				return; // Sufficient function decrease.
			else {// Backtrack.
				if (alam == 1.0)
					tmplam = -slope / (2.0 * (f_ln - fold - slope)); // First
																		// time.
				else {// Subsequent backtracks.
					rhs1 = f_ln - fold - alam * slope;
					rhs2 = f2 - fold - alam2 * slope;
					a = (rhs1 / (alam * alam) - rhs2 / (alam2 * alam2))
							/ (alam - alam2);
					b = (-alam2 * rhs1 / (alam * alam) + alam * rhs2
							/ (alam2 * alam2))
							/ (alam - alam2);
					if (a == 0.0)
						tmplam = -slope / (2.0 * b);
					else {
						disc = b * b - 3.0 * a * slope;
						if (disc < 0.0)
							tmplam = 0.5 * alam;
						else if (b <= 0.0)
							tmplam = (-b + Math.sqrt(disc)) / (3.0 * a);
						else
							tmplam = -slope / (b + Math.sqrt(disc));
					}
					if (tmplam > 0.5 * alam)
						tmplam = 0.5 * alam; // ? = 0.5?1.
				}
			}
			alam2 = alam;
			f2 = f_ln;
			alam = Math.max(tmplam, 0.1 * alam); // ? = 0.1?1.
		}// Try again.
	}

	/**
	 * Used internally. Computes forward-difference approximation to Jacobian. On input, x[1..n] is the point at 
	 * which the Jacobian is to be evaluated, fvec[1..n] is the vector of function values at the 
	 * point, and vecfunc(n,x,f) is a user-supplied routine that returns the vector of functions at 
	 * x. On output, df[1..n][1..n] is the Jacobian array.
	 * @param n n
	 * @param x x
	 * @param fvec fvec
	 * @param nameS nameS
	 * @return the result
	 */
	public static double[][] fdjac(int n, double[] x, double[] fvec,
			String nameS)// , double[][] df)
	// Computes forward-difference approximation to Jacobian. On input, x[1..n]
	// is the point at
	// which the Jacobian is to be evaluated, fvec[1..n] is the vector of
	// function values at the
	// point, and vecfunc(n,x,f) is a user-supplied routine that returns the
	// vector of functions at
	// x. On output, df[1..n][1..n] is the Jacobian array.
	{
		double[][] df = new double[n][n];

		int i = 0;
		int j = 0;
		double h = 0.0;
		double temp = 0.0;
		double[] f = new double[n];
		// f=vector(1,n);
		for (j = 1; j <= n; j++) {
			temp = x[j - 1];// temp=x[j];
			h = EPS2 * Math.abs(temp);
			if (h == 0.0)
				h = EPS2;
			x[j - 1] = temp + h;// x[j]=temp+h; Trick to reduce finite precision
								// error.
			h = x[j - 1] - temp;// h=x[j]-temp;
			// f=func.vecfunc(n,x);//(*vecfunc)(n,x,f);
			if (nameS.compareTo("shoot") == 0)
				shoot(n, x, f);
			else if (nameS.compareTo("shootf") == 0)
				shootf(n, x, f);
			x[j - 1] = temp;// x[j]=temp;
			for (i = 1; i <= n; i++)
				// df[i][j]=(f[i]-fvec[i])/h; Forward difference formula.
				df[i - 1][j - 1] = (f[i - 1] - fvec[i - 1]) / h;
		}
		// free_vector(f,1,n);
		return df;
	}

	/**
	 * Used internally. Returns f = 1/2 F · F at x. The vecfunc is a routine that returns the 
	 * vector of functions at x. It is set to point to a user-supplied routine in the calling program. 
	 * Global variables also communicate the function values back to the calling program.
	 * @param x x
	 * @param nameS nameS
	 * @return the result
	 */
	public static double fmin(double x[], String nameS)
	// Returns f = 1/2 F · F at x. The global pointer *nrfuncv points to a
	// routine that returns the
	// vector of functions at x. It is set to point to a user-supplied routine
	// in the calling program.
	// Global variables also communicate the function values back to the calling
	// program.
	{
		int i = 0;
		double sum = 0.0;
		// (*nrfuncv)(nn,x,fvec);public static void shoot(int n, double[] v,
		// double[] f)
		// fvec_newt=func.vecfunc(nn_newt,x);
		if (nameS.compareTo("shoot") == 0)
			shoot(nn_newt, x, fvec_newt);
		else if (nameS.compareTo("shootf") == 0)
			shootf(nn_newt, x, fvec_newt);
		for (sum = 0.0, i = 1; i <= nn_newt; i++)
			// sum += (fvec_newt[i])*(fvec_newt[i]);
			sum += (fvec_newt[i - 1]) * (fvec_newt[i - 1]);
		return 0.5 * sum;
	}

	/**
	 * Used internally. Local function.
	 * @param x x
	 * @param s s
	 * @param nameS nameS
	 * @return the result
	 */
	public static double localFunc(double[] x, String s, String nameS) {
		double res = 0.0;
		if (s.compareTo("fmin") == 0)
			res = fmin(x, nameS);
		return res;
	}
}
