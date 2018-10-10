package danfulea.math.numerical;

/**
 * Linear Algebric Equations class<br>
 * Based on Numerical Recipes in C (Cambridge Univ.).
 * 
 * @author Dan Fulea, 26 SEP. 2006.
 */
public class LinAEq {
	public static boolean failB = false;
	public static String failS = "";

	public static final double TINY = 1.0E-20;// a small number
	public static double dlu = 0.0;// used in LU decomposition

	// a11x1 + a12x2 + a13x3 + · · · + a1NxN = b1
	// a21x1 + a22x2 + a23x3 + · · · + a2NxN = b2
	// a31x1 + a32x2 + a33x3 + · · · + a3NxN = b3
	// · · · · · ·
	// aM1x1 + aM2x2 + aM3x3 + · · · + aMNxN = bM

	// @@Gauss Jordan Elimination method
	/*
	 * its principal weaknesses are (i) that it requires all the right-hand
	 * sides to be stored and manipulated at the same time, and (ii) that when
	 * the inverse matrix is not desired, Gauss-Jordan is three times slower
	 * than the best alternative technique for solving a single linear set . The
	 * method’s principal strength is that it is as stable as any other direct
	 * method, perhaps even a bit more stable when full pivoting is used. If you
	 * come along later with an additional right-hand side vector, you can
	 * multiply it by the inverse matrix, of course. This does give an answer,
	 * but one that is quite susceptible to roundoff error, not nearly as good
	 * as if the new vector had been included with the set of right-hand side
	 * vectors in the first instance. For these reasons, Gauss-Jordan
	 * elimination should usually not be your method of first choice, either for
	 * solving linear equations or for matrix inversion. The decomposition
	 * methods are better. Why do we give you Gauss-Jordan at all? Because it is
	 * straightforward, understandable, solid as a rock, and an exceptionally
	 * good 'psychological' backup for those times that something is going wrong
	 * and you think it might be your linear-equation solver.
	 */
	// #define SWAP(a,b) {temp=(a);(a)=(b);(b)=temp;}
	// public static void gaussj(float **a, int n, float **b, int m)
	/**
	 * Linear equation solution by Gauss-Jordan elimination. a[1..n][1..n] 
	 * is the input matrix. b[1..n][1..m] is input containing the m right-hand 
	 * side vectors. On output, a is replaced by its matrix inverse, and b is 
	 * replaced by the corresponding set of solution vectors. AX=B it is solved 
	 * simultaneously m linear equation systems. If m=1 we have one system of equations.
	 * @param a a
	 * @param n n
	 * @param b b
	 * @param m m
	 */
	public static void gaussj(double[][] a, int n, double[][] b, int m) {
		// Linear equation solution by Gauss-Jordan elimination. a[1..n][1..n]
		// is the input matrix. b[1..n][1..m] is input containing the m
		// right-hand
		// side vectors. On output, a is replaced by its matrix inverse, and b
		// is
		// replaced by the corresponding set of solution vectors.

		int j = 0;
		int i = 0;
		int k = 0;
		int icol = 0;
		int irow = 0;
		int l = 0;
		int ll = 0;
		int[] ipiv = new int[n];
		int[] indxc = new int[n];
		int[] indxr = new int[n];
		double big = 0.0;
		double temp = 0.0;
		double dum = 0.0;
		double pivinv = 0.0;

		failB = false;
		for (j = 1; j <= n; j++)
			ipiv[j - 1] = 0;// ipiv[j]=0;
		for (i = 1; i <= n; i++) {
			// This is the main loop over the columns to be reduced.
			big = 0.0;
			for (j = 1; j <= n; j++)
				// This is the outer loop of the search for a pivot element.
				// if (ipiv[j] != 1)
				if (ipiv[j - 1] != 1)
					for (k = 1; k <= n; k++) {
						if (ipiv[k - 1] == 0)// if (ipiv[k] == 0)
						{
							// if (Math.abs(a[j][k]) >= big)
							if (Math.abs(a[j - 1][k - 1]) >= big) {
								big = Math.abs(a[j - 1][k - 1]);// big=fabs(a[j][k]);
								irow = j;
								icol = k;
							}
						}
					}
			++(ipiv[icol - 1]);// ++(ipiv[icol]);
			// We now have the pivot element, so we interchange rows, if needed,
			// to put the pivot
			// element on the diagonal. The columns are not physically
			// interchanged, only relabeled:
			// indxc[i], the column of the ith pivot element, is the ith column
			// that is reduced, while
			// indxr[i] is the row in which that pivot element was originally
			// located. If indxr[i] =
			// indxc[i] there is an implied column interchange. With this form
			// of bookkeeping, the
			// solution b’s will end up in the correct order, and the inverse
			// matrix will be scrambled
			// by columns.
			if (irow != icol) {
				for (l = 1; l <= n; l++) // SWAP(a[irow][l],a[icol][l])
				{
					temp = a[irow - 1][l - 1];
					a[irow - 1][l - 1] = a[icol - 1][l - 1];
					a[icol - 1][l - 1] = temp;
				}
				for (l = 1; l <= m; l++) // SWAP(b[irow][l],b[icol][l])
				{
					temp = b[irow - 1][l - 1];
					b[irow - 1][l - 1] = b[icol - 1][l - 1];
					b[icol - 1][l - 1] = temp;
				}
			}
			indxr[i - 1] = irow;// indxr[i]=irow;
			// We are now ready to divide the pivot row by the
			// pivot element, located at irow and icol.
			indxc[i - 1] = icol;// indxc[i]=icol;
			// if (a[icol][icol] == 0.0) nrerror("gaussj: Singular Matrix");
			if (a[icol - 1][icol - 1] == 0.0) {
				// nrerror("gaussj: Singular Matrix");
				// System.out.println("gaussj: Singular Matrix");
				failS = "gaussj: Singular Matrix";
				failB = true;

				return;
			}
			pivinv = 1.0 / a[icol - 1][icol - 1];// pivinv=1.0/a[icol][icol];
			a[icol - 1][icol - 1] = 1.0;// a[icol][icol]=1.0;
			for (l = 1; l <= n; l++)
				a[icol - 1][l - 1] *= pivinv;// a[icol][l] *= pivinv;
			for (l = 1; l <= m; l++)
				b[icol - 1][l - 1] *= pivinv;// b[icol][l] *= pivinv;
			for (ll = 1; ll <= n; ll++)
				// Next, we reduce the rows...
				if (ll != icol) { // ...except for the pivot one, of course.
					dum = a[ll - 1][icol - 1];// dum=a[ll][icol];
					a[ll - 1][icol - 1] = 0.0;// a[ll][icol]=0.0;
					for (l = 1; l <= n; l++)
						a[ll - 1][l - 1] -= a[icol - 1][l - 1] * dum;// a[ll][l]
																		// -=
																		// a[icol][l]*dum;
					for (l = 1; l <= m; l++)
						b[ll - 1][l - 1] -= b[icol - 1][l - 1] * dum;// b[ll][l]
																		// -=
																		// b[icol][l]*dum;
				}
		}
		// This is the end of the main loop over columns of the reduction. It
		// only remains to unscramble
		// the solution in view of the column interchanges. We do this by
		// interchanging pairs of
		// columns in the reverse order that the permutation was built up.
		for (l = n; l >= 1; l--) {
			if (indxr[l - 1] != indxc[l - 1])// if (indxr[l] != indxc[l])
				for (k = 1; k <= n; k++) // SWAP(a[k][indxr[l]],a[k][indxc[l]]);
				{
					temp = a[k - 1][indxr[l - 1] - 1];
					a[k - 1][indxr[l - 1] - 1] = a[k - 1][indxc[l - 1] - 1];
					a[k - 1][indxc[l - 1] - 1] = temp;
				}
		} // And we are done.

		// free_ivector(ipiv,1,n);
		// free_ivector(indxr,1,n);
		// free_ivector(indxc,1,n);
	}

	// @@Gaussian elimination method
	/*
	 * The usefulness of Gaussian elimination with backsubstitution is primarily
	 * pedagogical. It stands between full elimination schemes such as
	 * Gauss-Jordan, and triangular decomposition schemes such as will be
	 * discussed in the next section. Gaussian elimination reduces a matrix not
	 * all the way to the identity matrix, but only halfway, to a matrix whose
	 * components on the diagonal and above (say) remain nontrivial. The
	 * advantage of Gaussian elimination and backsubstitution over Gauss-Jordan
	 * elimination is simply that the former is faster in raw operations count:
	 * Both Gaussian elimination and Gauss-Jordan elimination share the
	 * disadvantage that all right-hand sides must be known in advance. The LU
	 * decomposition method in the next section does not share that deficiency,
	 * and also has an equally small operations count, both for solution with
	 * any number of right-hand sides, and for matrix inversion. For this reason
	 * we will not implement the method of Gaussian elimination as a routine. =>
	 * TAKEN FROM T.BEU:
	 */
	/**
	 * Gaussian elimination with backsubstitution. The input data are the same as in gaussj.
	 * @param a a
	 * @param n n
	 * @param b b
	 * @param m m
	 * @return the determinant of a.
	 */
	public static double gauss(double[][] a, int n, double[][] b, int m) {
		int l = 0;
		double t = 0.0;
		double det = 1.0;

		// looking for not null pivot element!
		for (int k = 0; k < n; k++) {
			t = 0.0;
			for (int i = k; i < n; i++)
				if (t < Math.abs(a[i][k])) {
					t = Math.abs(a[i][k]);
					l = i;
				}
			// if required, swap l and k lines and
			// put pivot element on diagonal.
			if (l != k) {
				det = -det;// change det sign!
				for (int j = k; j < n; j++) {
					t = a[k][j];
					a[k][j] = a[l][j];
					a[l][j] = t;
				}
				for (int j = 0; j < m; j++) {
					t = b[k][j];
					b[k][j] = b[l][j];
					b[l][j] = t;
				}
			}
			// end pivot manipulation
			// divide pivot line
			t = 1 / a[k][k];
			det = det * a[k][k];
			for (int j = k + 1; j < n; j++)
				a[k][j] = a[k][j] * t;
			for (int j = 0; j < m; j++)
				b[k][j] = b[k][j] * t;
			// subtract non-pivot lines
			for (int i = k + 1; i < n; i++) {
				t = a[i][k];
				for (int j = k + 1; j < n; j++)
					a[i][j] = a[i][j] - a[k][j] * t;
				for (int j = 0; j < m; j++)
					b[i][j] = b[i][j] - b[k][j] * t;
			}
		}
		// build result array
		for (int k = n - 2; k >= 0; k--)
			// or->(int k=n-1; --k>=0; )//
			for (int j = 0; j < m; j++) {
				t = 0;
				for (int i = k + 1; i < n; i++)
					t = t + a[k][i] * b[i][j];
				b[k][j] = b[k][j] - t;
			}

		return det;
	}

	// @@ LU decomposition
	/*
	 * Suppose we are able to write the matrix A as a product of two matrices, L
	 * · U = A (2.3.1) where L is lower triangular (has elements only on the
	 * diagonal and below) and U is upper triangular (has elements only on the
	 * diagonal and above). For the case of a 4 × 4 matrix A, for example,
	 * equation (2.3.1) would look like this: (a11 0 0 0 a21 a22 0 0 a31 a32 a33
	 * 0 a41 a42 a43 a44) · ( ß11 ß12 ß13 ß14 0 ß22 ß23 ß24 0 0 ß33 ß34 0 0 0
	 * ß44 ) = ( a11 a12 a13 a14 a21 a22 a23 a24 a31 a32 a33 a34 a41 a42 a43 a44
	 * ) (2.3.2) We can use a decomposition such as (2.3.1) to solve the linear
	 * set A · x = (L · U) · x = L · (U · x) = b (2.3.3) by first solving for
	 * the vector y such that L · y = b (2.3.4) and then solving U · x = y
	 * (2.3.5) What is the advantage of breaking up one linear set into two
	 * successive ones? The advantage is that the solution of a triangular set
	 * of equations is quite trivial Thus, equation (2.3.4) can be solved by
	 * forward substitution as follows, y1 = b1/a11 yi =[1/aii][bi -sumj=1;j=i-1
	 * aijyj]; i = 2, 3, . . .,N (2.3.6) while (2.3.5) can then be solved by
	 * backsubstitution exactly, xN = yN/ßNN xi =[1/ßii] [yi -sumj=i+1;j=N
	 * ßijxj] i = N - 1,N - 2, . . . , 1 (2.3.7)
	 */
	/**
	 * Given a matrix a[1..n][1..n], this routine replaces it by the LU decomposition of a rowwise 
	 * permutation of itself [LU=A]. a and n are input. a is output, arranged as : 
	 * ß11 ß12 ß13 ß14;a21 ß22 ß23 ß24;a31 a32 ß33 ß34;a41 a42 a43 ß44. 
	 * indx[1..n] is an output vector that records the row permutation effected by the partial 
	 * pivoting; d is output as ±1 depending on whether the number of row interchanges was even 
	 * or odd, respectively. This routine is used in combination with lubksb to solve linear equations 
	 * or invert a matrix.
	 * @param a a
	 * @param n n
	 * @param indx indx
	 */
	public static void ludcmp(double[][] a, int n, int[] indx)// , float *d)
	{
		// Given a matrix a[1..n][1..n], this routine replaces it by the LU
		// decomposition of a rowwise
		// permutation of itself. a and n are input. a is output, arranged as :
		// ß11 ß12 ß13 ß14
		// a21 ß22 ß23 ß24
		// a31 a32 ß33 ß34
		// a41 a42 a43 ß44
		// indx[1..n] is an output vector that records the row permutation
		// effected by the partial
		// pivoting; d is output as ±1 depending on whether the number of row
		// interchanges was even
		// or odd, respectively. This routine is used in combination with lubksb
		// to solve linear equations
		// or invert a matrix.
		int i = 0;
		int imax = 0;
		int j = 0;
		int k = 0;
		double big = 0.0;
		double dum = 0.0;
		double sum = 0.0;
		double temp = 0.0;
		double[] vv = new double[n];// vv stores the implicit scaling of each
									// row.
		// vv=vector(1,n);
		failB = false;
		dlu = 1.0;// *d=1.0; No row interchanges yet.
		for (i = 1; i <= n; i++) { // Loop over rows to get the implicit scaling
									// information.
			big = 0.0;
			for (j = 1; j <= n; j++)
				// if ((temp=Math.abs(a[i][j])) > big) big=temp;
				if ((temp = Math.abs(a[i - 1][j - 1])) > big)
					big = temp;
			if (big == 0.0) {
				// nrerror("Singular matrix in routine ludcmp");
				// System.out.println("Singular matrix in routine ludcmp");
				failS = "Singular matrix in routine ludcmp";
				failB = true;

				return;
			}
			// No nonzero largest element.
			vv[i - 1] = 1.0 / big;// vv[i]=1.0/big; //Save the scaling.
		}
		for (j = 1; j <= n; j++) { // This is the loop over columns of Crout’s
									// method.
			for (i = 1; i < j; i++) {// This is equation (2.3.12) except for i =
										// j.
				sum = a[i - 1][j - 1];// sum=a[i][j];
				for (k = 1; k < i; k++)
					sum -= a[i - 1][k - 1] * a[k - 1][j - 1];// sum -=
																// a[i][k]*a[k][j];
				a[i - 1][j - 1] = sum;// a[i][j]=sum;
			}
			big = 0.0; // Initialize for the search for largest pivot element.
			for (i = j; i <= n; i++) { // This is i = j of equation (2.3.12) and
										// i = j+1. . .N
										// of equation (2.3.13).
				sum = a[i - 1][j - 1];// sum=a[i][j];
				for (k = 1; k < j; k++)
					sum -= a[i - 1][k - 1] * a[k - 1][j - 1];// sum -=
																// a[i][k]*a[k][j];
				a[i - 1][j - 1] = sum;// a[i][j]=sum;
				// if ( (dum=vv[i]*fabs(sum)) >= big) {
				if ((dum = vv[i - 1] * Math.abs(sum)) >= big) {
					// Is the figure of merit for the pivot better than the best
					// so far?
					big = dum;
					imax = i;
				}
			}
			if (j != imax) {// Do we need to interchange rows?
				for (k = 1; k <= n; k++) {// Yes, do so...
					dum = a[imax - 1][k - 1];// dum=a[imax][k];
					a[imax - 1][k - 1] = a[j - 1][k - 1];// a[imax][k]=a[j][k];
					a[j - 1][k - 1] = dum;// a[j][k]=dum;
				}
				dlu = -(dlu);// *d = -(*d); ...and change the parity of d.
				vv[imax - 1] = vv[j - 1];// vv[imax]=vv[j]; Also interchange the
											// scale factor.
			}
			indx[j - 1] = imax;// indx[j]=imax;
			if (a[j - 1][j - 1] == 0.0)
				a[j - 1][j - 1] = TINY;// if (a[j][j] == 0.0) a[j][j]=TINY;
			// If the pivot element is zero the matrix is singular (at least to
			// the precision of the
			// algorithm). For some applications on singular matrices, it is
			// desirable to substitute
			// TINY for zero.
			if (j != n) {// Now, finally, divide by the pivot element.
				dum = 1.0 / (a[j - 1][j - 1]);// dum=1.0/(a[j][j]);
				for (i = j + 1; i <= n; i++)
					a[i - 1][j - 1] *= dum;// for (i=j+1;i<=n;i++) a[i][j] *=
											// dum;
			}
		} // Go back for the next column in the reduction.
			// free_vector(vv,1,n);
	}

	// Here is the routine for forward substitution and backsubstitution,
	// implementing
	// equations (2.3.6) and (2.3.7).
	/**
	 * Solves the set of n linear equations A·X = B. Here a[1..n][1..n] is input, not as the matrix 
	 * A but rather as its LU decomposition, determined by the routine ludcmp. indx[1..n] is input 
	 * as the permutation vector returned by ludcmp. b[1..n] is input as the 
	 * right-hand side vector B, and returns with the solution vector X. a, n, and indx are not 
	 * modified by this routine and can be left in place for successive calls with different 
	 * right-hand sides b. This routine takes into account the possibility that b will begin with many zero 
	 * elements, so it is efficient for use in matrix inversion.
	 * @param a a
	 * @param n n
	 * @param indx indx
	 * @param b b
	 */
	public static void lubksb(double[][] a, int n, int[] indx, double[] b) {
		// Solves the set of n linear equations A·X = B. Here a[1..n][1..n] is
		// input, not as the matrix
		// A but rather as its LU decomposition, determined by the routine
		// ludcmp. indx[1..n] is input
		// as the permutation vector returned by ludcmp. b[1..n] is input as the
		// right-hand side vector
		// B, and returns with the solution vector X. a, n, and indx are not
		// modified by this routine
		// and can be left in place for successive calls with different
		// right-hand sides b. This routine takes
		// into account the possibility that b will begin with many zero
		// elements, so it is efficient for use
		// in matrix inversion.
		int i = 0;
		int ii = 0;
		int ip = 0;
		int j = 0;
		double sum = 0.0;
		for (i = 1; i <= n; i++) { // When ii is set to a positive value, it
									// will become the
									// index of the first nonvanishing element
									// of b. Wenow
									// do the forward substitution, equation
									// (2.3.6). The
									// only new wrinkle is to unscramble the
									// permutation as we go.
			ip = indx[i - 1];// ip=indx[i];
			sum = b[ip - 1];// sum=b[ip];
			b[ip - 1] = b[i - 1];// b[ip]=b[i];
			if (ii > 0)// if (ii)
				for (j = ii; j <= i - 1; j++)
					sum -= a[i - 1][j - 1] * b[j - 1];// for (j=ii;j<=i-1;j++)
														// sum -= a[i][j]*b[j];
			else if (sum != 0.0)
				ii = i;// if (sum) ii=i;
			// A nonzero element was encountered, so from now on we
			// will have to do the sums in the loop above.
			b[i - 1] = sum;// b[i]=sum;
		}
		for (i = n; i >= 1; i--) { // Now we do the backsubstitution, equation
									// (2.3.7).
			sum = b[i - 1];// sum=b[i];
			for (j = i + 1; j <= n; j++)
				sum -= a[i - 1][j - 1] * b[j - 1];// sum -= a[i][j]*b[j];
			b[i - 1] = sum / a[i - 1][i - 1];// b[i]=sum/a[i][i]; Store a
												// component of the solution
												// vector X.
		} // All done!
	}

	/*
	 * Using the above LU decomposition and backsubstitution routines, it is
	 * completely straightforward to find the inverse of a matrix column by
	 * column. The matrix y will now contain the inverse of the original matrix
	 * a, which will have been destroyed. Alternatively, there is nothing wrong
	 * with using a Gauss-Jordan routine like gaussj (§2.1) to invert a matrix
	 * in place, again destroying the original. Both methods have practically
	 * the same operations count.
	 */
	/**
	 * Using the above LU decomposition and backsubstitution routines, it is 
	 * completely straightforward to find the inverse of a matrix column by 
	 * column. The matrix y will now contain the inverse of the original matrix 
	 * a[1..N][1..N], which will have been destroyed.
	 * @param a a
	 * @param N N
	 * @param indx empty indx [1..N] to be used by LU routines.
	 * @return the result, y the inverse matrix of a.
	 */
	public static double[][] inverseMatrix(double[][] a, int N, int[] indx) {
		double[] col = new double[N];
		double[][] y = new double[N][N];
		ludcmp(a, N, indx);// ,&d); Decompose the matrix just once.
		if (failB)
			return y;// @@@@@@@@
		for (int j = 1; j <= N; j++) { // Find inverse by columns.
			for (int i = 1; i <= N; i++)
				col[i - 1] = 0.0;// col[i]=0.0;
			col[j - 1] = 1.0;// col[j]=1.0;
			lubksb(a, N, indx, col);
			for (int i = 1; i <= N; i++)
				y[i - 1][j - 1] = col[i - 1];// y[i][j]=col[i];
		}

		return y;
	}

	/*
	 * The determinant of an LU decomposed matrix is just the product of the
	 * diagonal elements,
	 */
	/**
	 * The determinant of an LU decomposed matrix is just the product of the 
	 * diagonal elements. Input data are the same as in inverseMatrix.
	 * @param a a
	 * @param N N
	 * @param indx indx
	 * @return the determinant
	 */
	public static double determinantMatrix(double[][] a, int N, int[] indx) {
		ludcmp(a, N, indx);// ,&d); //This returns d as ±1.
		if (failB)
			return -1;// Anyway the real output test is if failS<>""->naspa
		for (int j = 1; j <= N; j++)
			dlu *= a[j - 1][j - 1];// d *= a[j][j];

		return dlu;
	}

	/*
	 * @@TRIDIAG The special case of a system of linear equations that is
	 * tridiagonal, that is, has nonzero elements only on the diagonal plus or
	 * minus one column, is one that occurs frequently. Also common are systems
	 * that are band diagonal, with nonzero elements only along a few diagonal
	 * lines adjacent to the main diagonal (above and below). The resulting
	 * routine tridag is one that we will use in later chapters.
	 */
	/**
	 * Solves for a vector u[1..n] the tridiagonal linear set given by equation: 
	 * TRImatrix x uvecor = rvector. b is the main diagonal array. 
	 * TRImatrix have a2,b1,c1,...aN,bN,cN-1 tridiagonal elements (the rest are 0). 
	 * a[1..n], b[1..n], c[1..n], and r[1..n] are input vectors and are not modified. Solution 
	 * is stored in u.
	 * @param a a
	 * @param b b
	 * @param c c
	 * @param r r
	 * @param u u
	 * @param n n
	 */
	public static void tridag(double a[], double b[], double c[], double r[],
			double u[], int n) {
		/*
		 * ( b1 c1 0 · · · a2 b2 c2 · · · · · · · · · aN-1 bN-1 cN-1 · · · 0 aN
		 * bN) ·
		 * 
		 * (u1 u2 · · · uN-1 uN)
		 * 
		 * =
		 * 
		 * (r1 r2 · · · rN-1 rN )
		 */// (2.4.1)
			// Solves for a vector u[1..n] the tridiagonal linear set given by
			// equation (2.4.1). a[1..n],
			// b[1..n], c[1..n], and r[1..n] are input vectors and are not
			// modified.
		int j = 0;// unsigned long j;
		double bet = 0.0;
		double[] gam = new double[n];// gam=vector(1,n); One vector of
										// workspace, gam is needed.
		// if (b[1] == 0.0) nrerror("Error 1 in tridag");
		failB = false;

		if (b[0] == 0.0) {
			// nrerror("Error 1 in tridag");
			// System.out.println("Error 1 in tridag");
			failS = "Error 1 in tridag";
			failB = true;

			return;
		}
		// If this happens then you should rewrite your equations as a set of
		// order N - 1, with u2
		// trivially eliminated.
		u[0] = r[0] / (bet = b[0]);// u[1]=r[1]/(bet=b[1]);
		for (j = 2; j <= n; j++) {// Decomposition and forward substitution.
			gam[j - 1] = c[j - 2] / bet;// gam[j]=c[j-1]/bet;
			bet = b[j - 1] - a[j - 1] * gam[j - 1];// bet=b[j]-a[j]*gam[j];
			if (bet == 0.0) {// nrerror("Error 2 in tridag"); Algorithm fails;
								// see below.
								// System.out.println("Error 2 in tridag");
				failS = "Error 2 in tridag";
				failB = true;

				return;
			}
			u[j - 1] = (r[j - 1] - a[j - 1] * u[j - 2]) / bet;// u[j]=(r[j]-a[j]*u[j-1])/bet;
		}
		for (j = (n - 1); j >= 1; j--)
			u[j - 1] -= gam[j] * u[j];// u[j] -= gam[j+1]*u[j+1];
										// //Backsubstitution.
		// free_vector(gam,1,n);
	}

	// Alternative=>BEU
	/**
	 * Alternative approach for tridiagonal matrices. Solution is stored in d.
	 * @param a a
	 * @param b b
	 * @param c c
	 * @param d d
	 * @param n n
	 */
	public static void triDiag(double[] a, double[] b, double[] c, double[] d,
			int n) {
		double t = 0.0;
		c[0] = c[0] / b[0];
		d[0] = d[0] / b[0];
		for (int i = 1; i <= n - 2; i++) {
			t = b[i] - a[i] * c[i - 1];
			c[i] = c[i] / t;
			d[i] = (d[i] - a[i] * d[i - 1]) / t;
		}
		d[n - 1] = (d[n - 1] - a[n - 1] * d[n - 2])
				/ (b[n - 1] - a[n - 1] * c[n - 2]);
		for (int i = n - 2; i >= 0; i--)
			d[i] = d[i] - c[i] * d[i + 1];

	}

	// ################################BAND_NOT IMPLEMENTED
	/*
	 * @@@@LU IMPROVE Obviously it is not easy to obtain greater precision for
	 * the solution of a linear set than the precision of your computer’s
	 * floating-point word. Unfortunately, for large sets of linear equations,
	 * it is not always easy to obtain precision equal to, or even comparable
	 * to, the computer’s limit. In direct methods of solution, roundoff errors
	 * accumulate, and they are magnified to the extent that your matrix is
	 * close to singular. You can easily lose two or three significant figures
	 * for matrices which (you thought) were far from singular. If this happens
	 * to you, there is a neat trick to restore the full machine precision,
	 * called iterative improvement of the solution. The theory is very
	 * straightforward : Suppose that a vector x is the exact solution of the
	 * linear set A · x = b (2.5.1) You don’t, however, know x. You only know
	 * some slightly wrong solution x + dx, where dx is the unknownerror. When
	 * multiplied by the matrixA, your slightly wrong solution gives a product
	 * slightly discrepant from the desired right-hand side b, namely A · (x +
	 * dx) = b + db (2.5.2) Subtracting (2.5.1) from (2.5.2) gives A · dx = db
	 * But (2.5.2) can also be solved, trivially, for db. Substituting this into
	 * (2.5.3) gives A · dx = A · (x + dx) - b (2.5.4) In this equation, the
	 * whole right-hand side is known, since x + dx is the wrong solution that
	 * you want to improve. It is essential to calculate the right-hand side in
	 * double precision, since there will be a lot of cancellation in the
	 * subtraction of b. Then, we need only solve (2.5.4) for the error dx, then
	 * subtract this from the wrong solution to get an improved solution. An
	 * important extra benefit occurs if we obtained the original solution by LU
	 * decomposition. In this case we already have the LU decomposed form of A,
	 * and all we need do to solve (2.5.4) is compute the right-hand side and
	 * backsubstitute!
	 */
	/**
	 * Improves a solution vector x[1..n] of the linear set of equations AX = B. The matrix 
	 * a[1..n][1..n], and the vectors b[1..n] and x[1..n] are input, as is the dimension n. 
	 * Also input is alud[1..n][1..n], the LU decomposition of a as returned by ludcmp, and 
	 * the vector indx[1..n] also returned by that routine. On output, only x[1..n] is modified, 
	 * to an improved set of values.
	 * @param a a
	 * @param alud alud
	 * @param n n
	 * @param indx indx
	 * @param b b
	 * @param x x
	 */
	public static void mprove(double[][] a, double[][] alud, int n, int[] indx,
			double[] b, double[] x) {
		// Improves a solution vector x[1..n] of the linear set of equations A ·
		// X = B. The matrix
		// a[1..n][1..n], and the vectors b[1..n] and x[1..n] are input, as is
		// the dimension n.
		// Also input is alud[1..n][1..n], the LU decomposition of a as returned
		// by ludcmp, and
		// the vector indx[1..n] also returned by that routine. On output, only
		// x[1..n] is modified,
		// to an improved set of values.
		// void lubksb(float **a, int n, int *indx, float b[]);
		int j = 0;
		int i = 0;
		double sdp = 0.0;
		double[] r = new double[n];// r=vector(1,n);
		for (i = 1; i <= n; i++) {
			// Calculate the right-hand side, accumulating
			// the residual in double precision.
			sdp = -b[i - 1];// sdp = -b[i];
			for (j = 1; j <= n; j++)
				sdp += a[i - 1][j - 1] * x[j - 1];// sdp += a[i][j]*x[j];
			r[i - 1] = sdp;// r[i]=sdp;
		}
		lubksb(alud, n, indx, r); // Solve for the error term,
		for (i = 1; i <= n; i++)
			x[i - 1] -= r[i - 1];// x[i] -= r[i]; //and subtract it from the old
									// solution.
		// free_vector(r,1,n);
		// You should note that the routine ludcmp in §2.3 destroys the input
		// matrix as
		// it LU decomposes it. Since iterative improvement requires both the
		// original matrix
		// and its LU decomposition, you will need to copy A before calling
		// ludcmp. Likewise
		// lubksb destroys b in obtaining x, so make a copy of b also.
	}

	/**
	 * Improved LU method. The input data is the matrix a[1..n][1..n], an empty 
	 * indx[1..n] vector and right-hand vector b[1..n].
	 * @param a a
	 * @param n n
	 * @param indx indx
	 * @param b b
	 * @return the improved set of values.
	 */
	public static double[] improveLU(double[][] a, int n, int[] indx, double[] b) {
		double[][] alud = new double[n][n];
		double[][] acopy = new double[n][n];
		double[] bcopy = new double[n];
		double[] x = new double[n];
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				acopy[i][j] = a[i][j];
			}
			bcopy[i] = b[i];
		}
		ludcmp(a, n, indx);
		if (failB)
			return x;// !!!!!!!!!!!!!!!!!!!!!
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				alud[i][j] = a[i][j];
			}
			// bcopy[i]=b[i];
		}
		lubksb(a, n, indx, b);
		for (int i = 0; i < n; i++) {
			x[i] = b[i];// System.out.println(" i= "+i+" value= "+x[i]);//OK!!
		}

		mprove(acopy, alud, n, indx, bcopy, x);

		return x;// PERFECT!
	}

	/*
	 * @@Singular Value Decomposition There exists a very powerful set of
	 * techniques for dealing with sets of equations or matrices that are either
	 * singular or else numerically very close to singular. In many cases where
	 * Gaussian elimination and LU decomposition fail to give satisfactory
	 * results, this set of techniques, known as singular value decomposition,
	 * or SVD, will diagnose for you precisely what the problem is. In some
	 * cases, SVD will not only diagnose the problem, it will also solve it, in
	 * the sense of giving you a useful numerical answer, although, as we shall
	 * see, not necessarily “the” answer that you thought you should get. SVDis
	 * also the method of choice for solving most linear least-squares problems.
	 * We will outline the relevant theory in this section, but defer detailed
	 * discussion of the use of SVD in this application to Chapter 15, whose
	 * subject is the parametric modeling of data. SVD methods are based on the
	 * following theorem of linear algebra, whose proof is beyond our scope:
	 * AnyM ×N matrix A whose number of rowsM is greater than or equal to its
	 * number of columns N, can be written as the product of an M × N
	 * column-orthogonal matrix U, an N × N diagonal matrix W with positive or
	 * zero elements (the singular values), and the transpose of an N ×N
	 * orthogonal matrix V.
	 * 
	 * If the matrix A is square, N × N say, then U, V, andWare all square
	 * matrices of the same size. Their inverses are also trivial to compute: U
	 * and V are orthogonal, so their inverses are equal to their transposes; W
	 * is diagonal, so its inverse is the diagonal matrix whose elements are the
	 * reciprocals of the elements wj . The only thing that can go wrong with
	 * this construction is for one of the wj ’s to be zero, or (numerically)
	 * for it to be so small that its value is dominated by roundoff error and
	 * therefore unknowable. If more than one of the w j ’s have this problem,
	 * then the matrix is even more singular. So, first of all, SVD gives you a
	 * clear diagnosis of the situation.
	 * 
	 * In that case, the direct solution methods of LU decomposition or Gaussian
	 * elimination may actually give a formal solution to the set of equations
	 * (that is, a zero pivot may not be encountered); but the solution vector
	 * may have wildly large components whose algebraic cancellation, when
	 * multiplying by the matrix A, may give a very poor approximation to the
	 * right-hand vector b. In such cases, the solution vector x obtained by
	 * zeroing the small wj ’s and then using equation (2.6.7) is very often
	 * better (in the sense of the residual |A · x - b| being smaller) than both
	 * the direct-method solution and the SVD solution where the small wj ’s are
	 * left nonzero.
	 * 
	 * It may seem paradoxical that this can be so, since zeroing a singular
	 * value corresponds to throwing away one linear combination of the set of
	 * equations that we are trying to solve. The resolution of the paradox is
	 * that we are throwing away precisely a combination of equations that is so
	 * corrupted by roundoff error as to be at best useless; usually it is worse
	 * than useless since it “pulls” the solution vector way off towards
	 * infinity along some direction that is almost a nullspace vector. In doing
	 * this, it compounds the roundoff problem and makes the residual |A · x -
	 * b| larger. SVD cannot be applied blindly, then. You have to exercise some
	 * discretion in deciding at what threshold to zero the small wj ’s, and/or
	 * you have to have some idea what size of computed residual |A · x - b| is
	 * acceptable. As an example, here is a “backsubstitution” routine svbksb
	 * for evaluating equation (2.6.7) and obtaining a solution vector x from a
	 * right-hand side b, given that the SVD of a matrix A has already been
	 * calculated by a call to svdcmp. Note that this routine presumes that you
	 * have already zeroed the small wj ’s. It does not do this for you. If you
	 * haven’t zeroed the small wj ’s, then this routine is just as
	 * ill-conditioned as any direct method, and you are misusing SVD.
	 */

	/**
	 * Solves A·X = B for a vector X, where A is specified by the arrays 
	 * u[1..m][1..n], w[1..n], v[1..n][1..n] as returned by svdcmp. m and n are the dimensions of a, and 
	 * will be equal for square matrices. b[1..m] is the input right-hand side. x[1..n] is the 
	 * output solution vector. No input quantities are destroyed, so the routine may be called 
	 * sequentially with different b’s.
	 * @param u u
	 * @param w w
	 * @param v v
	 * @param m m
	 * @param n n
	 * @param b b
	 * @param x x
	 */
	public static void svbksb(double[][] u, double[] w, double[][] v, int m,
			int n, double[] b, double[] x)
	// Solves A·X = B for a vector X, where A is specified by the arrays
	// u[1..m][1..n], w[1..n],
	// v[1..n][1..n] as returned by svdcmp. m and n are the dimensions of a, and
	// will be equal for
	// square matrices. b[1..m] is the input right-hand side. x[1..n] is the
	// output solution vector.
	// No input quantities are destroyed, so the routine may be called
	// sequentially with different b’s.
	{
		int jj = 0;
		int j = 0;
		int i = 0;
		double s = 0.0;
		double[] tmp = new double[n];
		// tmp=vector(1,n);
		for (j = 1; j <= n; j++) { // Calculate U Tranpus x B.
			s = 0.0;
			if (w[j - 1] != 0.0)// if (w[j])
			{// Nonzero result only if wj is nonzero.
				for (i = 1; i <= m; i++)
					s += u[i - 1][j - 1] * b[i - 1];// s += u[i][j]*b[i];
				s /= w[j - 1];// s /= w[j]; This is the divide by wj .
			}
			tmp[j - 1] = s;// tmp[j]=s;
		}
		for (j = 1; j <= n; j++) {// Matrix multiply by V to get answer.
			s = 0.0;
			for (jj = 1; jj <= n; jj++)
				s += v[j - 1][jj - 1] * tmp[jj - 1];// s += v[j][jj]*tmp[jj];
			x[j - 1] = s;// x[j]=s;
		}
		// free_vector(tmp,1,n);
	}

	/**
	 * Given a matrix a[1..m][1..n], this routine computes its singular 
	 * value decomposition, A = U·W·V Transpus. The matrix U replaces a on output. The diagonal matrix 
	 * of singular values W is output as a vector w[1..n]. The matrix V (not the transpose V T ) is output 
	 * as v[1..n][1..n].
 	 * @param a a
	 * @param m m
	 * @param n n
	 * @param w w
	 * @param v v
	 */
	public static void svdcmp(double[][] a, int m, int n, double[] w,
			double[][] v) {
		// Given a matrix a[1..m][1..n], this routine computes its singular
		// value decomposition, A =
		// U·W·V Transpus. Thematrix U replaces a on output. The diagonal matrix
		// of singular values W is output
		// as a vector w[1..n]. Thematrix V (not the transpose V T ) is output
		// as v[1..n][1..n].

		// float pythag(float a, float b);
		int flag = 0;
		int i = 0;
		int its = 0;
		int j = 0;
		int jj = 0;
		int k = 0;
		int l = 0;
		int nm = 0;
		double anorm = 0.0;
		double c = 0.0;
		double f = 0.0;
		double g = 0.0;
		double h = 0.0;
		double s = 0.0;
		double scale = 0.0;
		double x = 0.0;
		double y = 0.0;
		double z = 0.0;// *rv1;
		double[] rv1 = new double[n];

		double gf = 0.0;

		failB = false;
		// Householder reduction to bidiagonal form.
		g = scale = anorm = 0.0;
		for (i = 1; i <= n; i++) {
			l = i + 1;
			rv1[i - 1] = scale * g;// rv1[i]=scale*g;
			g = s = scale = 0.0;
			if (i <= m) {
				for (k = i; k <= m; k++)
					scale += Math.abs(a[k - 1][i - 1]);// scale +=
														// fabs(a[k][i]);
				if (scale != 0.0)// if (scale)
				{
					for (k = i; k <= m; k++) {
						a[k - 1][i - 1] /= scale;// a[k][i] /= scale;
						s += a[k - 1][i - 1] * a[k - 1][i - 1];// s +=
																// a[k][i]*a[k][i];
					}
					f = a[i - 1][i - 1];// f=a[i][i];
					// g = -SIGN(sqrt(s),f);//cap11: g=(f >= 0.0 ? -sqrt(h) :
					// sqrt(h));
					g = (f >= 0.0 ? -Math.sqrt(s) : Math.sqrt(s));
					// DAR=>SIGN(a,b) Magnitude of a times sign of b.=>OK!!!
					h = f * g - s;
					a[i - 1][i - 1] = f - g;// a[i][i]=f-g;
					for (j = l; j <= n; j++) {
						for (s = 0.0, k = i; k <= m; k++)
							s += a[k - 1][i - 1] * a[k - 1][j - 1];// s +=
																	// a[k][i]*a[k][j];
						f = s / h;
						for (k = i; k <= m; k++)
							a[k - 1][j - 1] += f * a[k - 1][i - 1];// a[k][j] +=
																	// f*a[k][i];
					}
					for (k = i; k <= m; k++)
						a[k - 1][i - 1] *= scale;// a[k][i] *= scale;
				}
			}
			w[i - 1] = scale * g;// w[i]=scale *g;
			g = s = scale = 0.0;
			if (i <= m && i != n) {
				for (k = l; k <= n; k++)
					scale += Math.abs(a[i - 1][k - 1]);// scale +=
														// fabs(a[i][k]);
				if (scale != 0.0)// if (scale)
				{
					for (k = l; k <= n; k++) {
						a[i - 1][k - 1] /= scale;// a[i][k] /= scale;
						s += a[i - 1][k - 1] * a[i - 1][k - 1];// s +=
																// a[i][k]*a[i][k];
					}
					f = a[i - 1][l - 1];// f=a[i][l];
					// g = -SIGN(sqrt(s),f);
					g = (f >= 0.0 ? -Math.sqrt(s) : Math.sqrt(s));
					h = f * g - s;
					a[i - 1][l - 1] = f - g;// a[i][l]=f-g;
					for (k = l; k <= n; k++)
						rv1[k - 1] = a[i - 1][k - 1] / h;// rv1[k]=a[i][k]/h;
					for (j = l; j <= m; j++) {
						for (s = 0.0, k = l; k <= n; k++)
							s += a[j - 1][k - 1] * a[i - 1][k - 1];// s +=
																	// a[j][k]*a[i][k];
						for (k = l; k <= n; k++)
							a[j - 1][k - 1] += s * rv1[k - 1];// a[j][k] +=
																// s*rv1[k];
					}
					for (k = l; k <= n; k++)
						a[i - 1][k - 1] *= scale;// a[i][k] *= scale;
				}
			}
			// anorm=Math.max(anorm,(fabs(w[i])+fabs(rv1[i])));
			anorm = Math
					.max(anorm, (Math.abs(w[i - 1]) + Math.abs(rv1[i - 1])));
		}

		for (i = n; i >= 1; i--) {// Accumulation of right-hand transformations.
			if (i < n) {
				if (g != 0.0)// if (g)
				{
					for (j = l; j <= n; j++)
						// Double division to avoid possible underflow.
						v[j - 1][i - 1] = (a[i - 1][j - 1] / a[i - 1][l - 1])
								/ g;// v[j][i]=(a[i][j]/a[i][l])/g;
					for (j = l; j <= n; j++) {
						for (s = 0.0, k = l; k <= n; k++)
							s += a[i - 1][k - 1] * v[k - 1][j - 1];// s +=
																	// a[i][k]*v[k][j];
						for (k = l; k <= n; k++)
							v[k - 1][j - 1] += s * v[k - 1][i - 1];// v[k][j] +=
																	// s*v[k][i];
					}
				}
				for (j = l; j <= n; j++)
					v[i - 1][j - 1] = v[j - 1][i - 1] = 0.0;// v[i][j]=v[j][i]=0.0;
			}
			v[i - 1][i - 1] = 1.0;// v[i][i]=1.0;
			g = rv1[i - 1];// g=rv1[i];
			l = i;
		}

		for (i = Math.min(m, n); i >= 1; i--) { // Accumulation of left-hand
												// transformations.
			l = i + 1;
			g = w[i - 1];// g=w[i];
			for (j = l; j <= n; j++)
				a[i - 1][j - 1] = 0.0;// a[i][j]=0.0;
			if (g != 0.0) // if (g)
			{
				g = 1.0 / g;
				for (j = l; j <= n; j++) {
					for (s = 0.0, k = l; k <= m; k++)
						s += a[k - 1][i - 1] * a[k - 1][j - 1];// s +=
																// a[k][i]*a[k][j];
					f = (s / a[i - 1][i - 1]) * g;// f=(s/a[i][i])*g;
					for (k = i; k <= m; k++)
						a[k - 1][j - 1] += f * a[k - 1][i - 1];// a[k][j] +=
																// f*a[k][i];
				}
				for (j = i; j <= m; j++)
					a[j - 1][i - 1] *= g;// a[j][i] *= g;
			} else
				for (j = i; j <= m; j++)
					a[j - 1][i - 1] = 0.0;// a[j][i]=0.0;
			++a[i - 1][i - 1];// ++a[i][i];
		}

		for (k = n; k >= 1; k--) {
			// Diagonalization of the bidiagonal form: Loop over
			// singular values, and over allowed iterations.
			for (its = 1; its <= 30; its++) {
				flag = 1;
				for (l = k; l >= 1; l--) {// Test for splitting.
					nm = l - 1; // Note that rv1[1] is always zero.=>OK for *
								// (is breakeed before)!
					// if ((float)(fabs(rv1[l])+anorm) == anorm) {
					if ((double) (Math.abs(rv1[l - 1]) + anorm) == anorm) {
						flag = 0;
						break;
					}
					// if ((float)(fabs(w[nm])+anorm) == anorm) break;
					if ((double) (Math.abs(w[nm - 1]) + anorm) == anorm)
						break;// *
				}
				if (flag != 0)// if (flag)
				{
					c = 0.0; // Cancellation of rv1[l], if l > 1.
					s = 1.0;
					for (i = l; i <= k; i++) {
						f = s * rv1[i - 1];// f=s*rv1[i];
						rv1[i - 1] = c * rv1[i - 1];// rv1[i]=c*rv1[i];
						// if ((float)(fabs(f)+anorm) == anorm) break;
						if ((double) (Math.abs(f) + anorm) == anorm)
							break;
						g = w[i - 1];// g=w[i];
						h = pythag(f, g);
						w[i - 1] = h;// w[i]=h;
						h = 1.0 / h;
						c = g * h;
						s = -f * h;
						for (j = 1; j <= m; j++) {
							y = a[j - 1][nm - 1];// y=a[j][nm];
							z = a[j - 1][i - 1];// z=a[j][i];
							a[j - 1][nm - 1] = y * c + z * s;// a[j][nm]=y*c+z*s;
							a[j - 1][i - 1] = z * c - y * s;// a[j][i]=z*c-y*s;
						}
					}
				}
				z = w[k - 1];// z=w[k];
				if (l == k) {// Convergence.
					if (z < 0.0) {// Singular value is made nonnegative.
						w[k - 1] = -z;// w[k] = -z;
						for (j = 1; j <= n; j++)
							v[j - 1][k - 1] = -v[j - 1][k - 1];// v[j][k] =
																// -v[j][k];
					}
					break;
				}
				if (its == 30) {
					// nrerror("no convergence in 30 svdcmp iterations");
					// System.out.println("no convergence in 30 svdcmp iterations");
					failS = "no convergence in 30 svdcmp iterations";
					failB = true;

					return;
				}
				x = w[l - 1];// x=w[l]; //Shift from bottom 2-by-2 minor.
				nm = k - 1;
				y = w[nm - 1];// y=w[nm];
				g = rv1[nm - 1];// g=rv1[nm];
				h = rv1[k - 1];// h=rv1[k];
				f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y);
				g = pythag(f, 1.0);
				// @@@@@@@@@@@@@@@DED
				gf = (f >= 0.0 ? Math.sqrt(g) : -Math.sqrt(g));
				// @@@@@@@@@@@@@
				// f=((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x;
				f = ((x - z) * (x + z) + h * ((y / (f + gf)) - h)) / x;
				c = s = 1.0; // Next QR transformation:
				for (j = l; j <= nm; j++) {
					i = j + 1;
					g = rv1[i - 1];// g=rv1[i];
					y = w[i - 1];// y=w[i];
					h = s * g;
					g = c * g;
					z = pythag(f, h);
					rv1[j - 1] = z;// rv1[j]=z;
					c = f / z;
					s = h / z;
					f = x * c + g * s;
					g = g * c - x * s;
					h = y * s;
					y *= c;
					for (jj = 1; jj <= n; jj++) {
						x = v[jj - 1][j - 1];// x=v[jj][j];
						z = v[jj - 1][i - 1];// z=v[jj][i];
						v[jj - 1][j - 1] = x * c + z * s;// v[jj][j]=x*c+z*s;
						v[jj - 1][i - 1] = z * c - x * s;// v[jj][i]=z*c-x*s;
					}
					z = pythag(f, h);
					w[j - 1] = z;// w[j]=z; //Rotation can be arbitrary if z =
									// 0.
					if (z != 0.0)// if (z)
					{
						z = 1.0 / z;
						c = f * z;
						s = h * z;
					}
					f = c * g + s * y;
					x = c * y - s * g;

					for (jj = 1; jj <= m; jj++) {
						y = a[jj - 1][j - 1];// y=a[jj][j];
						z = a[jj - 1][i - 1];// z=a[jj][i];
						a[jj - 1][j - 1] = y * c + z * s;// a[jj][j]=y*c+z*s;
						a[jj - 1][i - 1] = z * c - y * s;// a[jj][i]=z*c-y*s;
					}
				}

				rv1[l - 1] = 0.0;// rv1[l]=0.0;
				rv1[k - 1] = f;// rv1[k]=f;
				w[k - 1] = x;// w[k]=x;
			}
		}
		// free_vector(rv1,1,n);
	}

	/**
	 * Computes (a^2 + b^2)^1/2 without destructive underflow or overflow.
	 * @param a a
	 * @param b b
	 * @return the result
	 */
	public static double pythag(double a, double b)
	// Computes (a2 + b2)1/2 without destructive underflow or overflow.
	{
		double absa = 0.0;
		double absb = 0.0;
		absa = Math.abs(a);
		absb = Math.abs(b);
		if (absa > absb)
			return absa * Math.sqrt(1.0 + (absb / absa) * (absb / absa));
		// else return (absb == 0.0 ? 0.0 : absb*sqrt(1.0+SQR(absa/absb)));
		else
			return (absb == 0.0 ? 0.0 : absb
					* Math.sqrt(1.0 + (absa / absb) * (absa / absb)));
	}

	// finaly how to use:
	/*
	 * #define N ... float wmax,wmin,**a,**u,*w,**v,*b,*x; int i,j; ...
	 * for(i=1;i<=N;i++) Copy a into u if you don’t want it to be destroyed. for
	 * j=1;j<=N;j++) u[i][j]=a[i][j]; svdcmp(u,N,N,w,v); SVD the square matrix
	 * a. wmax=0.0; Will be the maximum singular value obtained.
	 * for(j=1;j<=N;j++) if (w[j] > wmax) wmax=w[j]; This is where we set the
	 * threshold for singular values allowed to be nonzero. The constant is
	 * typical, but not universal. You have to experiment with your own
	 * application. wmin=wmax*1.0e-6; for(j=1;j<=N;j++) if (w[j] < wmin)
	 * w[j]=0.0; svbksb(u,w,v,N,N,b,x); Now we can backsubstitute.
	 */
	/**
	 * A more convenient use of SVD method. Solves AX=B.
	 * @param a input matrix a [1..N][1..N]
	 * @param N N
	 * @param b input right-hand vector b[1,,N]
	 * @param x the solution vector X.
	 */
	public static void svd(double[][] a, int N, double[] b, double[] x) {
		double[][] u = new double[N][N];
		double[] w = new double[N];// of n; a = matrix of m x n
		double[][] v = new double[N][N];// of n x n!!

		for (int i = 1; i <= N; i++)
			// Copy a into u if you don’t want it to be destroyed.
			for (int j = 1; j <= N; j++)
				u[i - 1][j - 1] = a[i - 1][j - 1];
		svdcmp(u, N, N, w, v); // SVD the square matrix a.
		if (failB)
			return;
		double wmax = 0.0; // Will be the maximum singular value obtained.
		for (int j = 1; j <= N; j++)
			if (w[j - 1] > wmax)
				wmax = w[j - 1];
		// This is where we set the threshold for singular values allowed to be
		// nonzero. The constant
		// is typical, but not universal. You have to experiment with your own
		// application.
		double wmin = wmax * 1.0e-6;
		for (int j = 1; j <= N; j++)
			if (w[j - 1] < wmin)
				w[j - 1] = 0.0;
		svbksb(u, w, v, N, N, b, x); // Now we can backsubstitute.
	}

	// ################################CYCLIC_NOT IMPLEMENTED
	// ################################INDEXED STORAGE_SPARSE SYSTEM_NOT
	// IMPLEMENTED
	// ################################CONJUGATED_GRADIENT_NOT IMPLEMENTED

	// @@VANDERMOND MATRICES
	/*
	 * A Vandermonde matrix of size N × N is completely determined by N
	 * arbitrary numbers x1, x2, . . . , xN, in terms of which its N^2
	 * components are the integer powers xi^j-1 , i,j = 1, . . . , N. Evidently
	 * there are two possible such forms, depending on whether we view the i’s
	 * as rows, j’s as columns, or vice versa. In the former case, we get a
	 * linear system of equations that looks like this, ( 1 x1 x1^2 · · · x1^N-1
	 * 1 x2 x2^2· · · x2^N-1 ............ 1 xN xN^2 · · · xN^N-1) · (c1 c2 ...
	 * cN) = (y1 y2 ... yN )
	 */
	/**
	 * Solves the Vandermonde linear system sum from i=1 to N xi^(k-1)wi = qk (k = 1,...,N). 
	 * Input consists of the vectors x[1..n] and q[1..n]; the vector w[1..n] is output.
	 * @param x x
	 * @param w w
	 * @param q q
	 * @param n n
	 */
	public static void vander(double[] x, double[] w, double[] q, int n)
	// Solves the Vandermonde linear system sum from i=1 to N xi^(k-1)wi = qk (k
	// = 1,...,N).
	// Input consists of the vectors x[1..n] and q[1..n]; the vector w[1..n] is
	// output.
	{
		int i = 0;
		int j = 0;
		int k = 0;
		double b = 0.0;
		double s = 0.0;
		double t = 0.0;
		double xx = 0.0;
		double[] c = new double[n];
		// c=dvector(1,n);
		if (n == 1)
			w[0] = q[0];// w[1]=q[1];
		else {
			for (i = 1; i <= n; i++)
				c[i - 1] = 0.0;// c[i]=0.0; //Initialize array.
			c[n - 1] = -x[0];// c[n] = -x[1]; //Coefficients of the master
								// polynomial are found by recursion.
			for (i = 2; i <= n; i++) {
				xx = -x[i - 1];// xx = -x[i];
				for (j = (n + 1 - i); j <= (n - 1); j++)
					c[j - 1] += xx * c[j];// c[j] += xx*c[j+1];
				c[n - 1] += xx;// c[n] += xx;
			}
			for (i = 1; i <= n; i++) {// Each subfactor in turn
				xx = x[i - 1];// xx=x[i];
				t = b = 1.0;
				s = q[n - 1];// s=q[n];
				for (k = n; k >= 2; k--) {// is synthetically divided,
					b = c[k - 1] + xx * b;// b=c[k]+xx*b;
					s += q[k - 2] * b;// s += q[k-1]*b;// matrix-multiplied by
										// the right-hand side,
					t = xx * t + b;
				}
				w[i - 1] = s / t;// w[i]=s/t; //and supplied with a denominator.
			}
		}
		// free_dvector(c,1,n);
	}

	// @@Toeplitz matrices
	/*
	 * An N × N Toeplitz matrix is specified by giving 2N - 1 numbers Rk, k = -N
	 * + 1, . . . ,-1, 0, 1, . . . , N - 1. Those numbers are then emplaced as
	 * matrix elements constant along the (upper-left to lower-right) diagonals
	 * of the matrix: ( R0 R-1 R-2 · · · R-(N-2) R-(N-1) R1 R0 R-1 · · · R-(N-3)
	 * R-(N-2) R2 R1 R0 · · · R-(N-4) R-(N-3) · · · · · · RN-2 RN-3 RN-4 · · ·
	 * R0 R-1 RN-1 RN-2 RN-3 · · · R1 R0 ) (2.8.8) The linear Toeplitz problem
	 * can thus be written as N sumj=1[R de index (i-j) ori xj] = yi (i = 1, . .
	 * . ,N) (2.8.9) where the xj ’s, j = 1, . . . , N, are the unknowns to be
	 * solved for. The Toeplitz matrix is symmetric if Rk = R-k for all k.
	 */
	// #define FREERETURN {free_vector(h,1,n);free_vector(g,1,n);return;}
	/**
	 * Solves the Toeplitz system sum j=1 to N of R(N+i-j)xj = yi (i = 1, . . . ,N). The Toeplitz matrix need 
	 * not be symmetric. y[1..n] and r[1..2*n-1] are input arrays; x[1..n] is the output array. 
	 * @param r r
	 * @param x x
	 * @param y y
	 * @param n n
	 */
	public static void toeplz(double[] r, double[] x, double[] y, int n)
	// Solves the Toeplitz system sum j=1 to N of
	// R(N+i-j)xj = yi (i = 1, . . . ,N). The Toeplitz matrix need
	// not be symmetric. y[1..n] and r[1..2*n-1] are input arrays;
	// x[1..n] is the output array.
	{
		int j = 0;
		int k = 0;
		int m = 0;
		int m1 = 0;
		int m2 = 0;
		double pp = 0.0;
		double pt1 = 0.0;
		double pt2 = 0.0;
		double qq = 0.0;
		double qt1 = 0.0;
		double qt2 = 0.0;
		double sd = 0.0;
		double sgd = 0.0;
		double sgn = 0.0;
		double shn = 0.0;
		double sxn = 0.0;
		double[] g = new double[n];
		double[] h = new double[n];
		// if (r[n] == 0.0) nrerror("toeplz-1 singular principal minor");

		failB = false;

		if (r[n - 1] == 0.0) {
			// nrerror("toeplz-1 singular principal minor");
			// System.out.println("toeplz-1 singular principal minor");
			failS = "toeplz-1 singular principal minor";
			failB = true;

			return;
		}
		// g=vector(1,n);
		// h=vector(1,n);
		x[0] = y[0] / r[n - 1];// x[1]=y[1]/r[n]; Initialize for the recursion.
		if (n == 1)
			return;// FREERETURN
		g[0] = r[n - 2] / r[n - 1];// g[1]=r[n-1]/r[n];
		h[0] = r[n] / r[n - 1];// h[1]=r[n+1]/r[n];
		for (m = 1; m <= n; m++) {// Main loop over the recursion.
			m1 = m + 1;
			sxn = -y[m1 - 1];// sxn = -y[m1]; Compute numerator and denominator
								// for x,
			sd = -r[n - 1];// sd = -r[n];
			for (j = 1; j <= m; j++) {
				sxn += r[n + m1 - j - 1] * x[j - 1];// sxn += r[n+m1-j]*x[j];
				sd += r[n + m1 - j - 1] * g[m - j];// sd += r[n+m1-j]*g[m-j+1];
			}
			if (sd == 0.0) {
				// nrerror("toeplz-2 singular principal minor");
				// System.out.println("toeplz-2 singular principal minor");
				failS = "toeplz-2 singular principal minor";
				failB = true;

				return;
			}
			x[m1 - 1] = sxn / sd;// x[m1]=sxn/sd; whence x.
			for (j = 1; j <= m; j++)
				x[j - 1] -= x[m1 - 1] * g[m - j];// x[j] -= x[m1]*g[m-j+1];
			if (m1 == n)
				return;// FREERETURN->OK@@@@@@@@@@@@!!!!!!!!!!!!!!!!!!!11
			sgn = -r[n - m1 - 1];// sgn = -r[n-m1]; //Compute numerator and
									// denominator for G and H,
			shn = -r[n + m1 - 1];// shn = -r[n+m1];
			sgd = -r[n - 1];// sgd = -r[n];
			for (j = 1; j <= m; j++) {
				sgn += r[n + j - m1 - 1] * g[j - 1];// sgn += r[n+j-m1]*g[j];
				shn += r[n + m1 - j - 1] * h[j - 1];// shn += r[n+m1-j]*h[j];
				sgd += r[n + j - m1 - 1] * h[m - j];// sgd +=
													// r[n+j-m1]*h[m-j+1];
			}
			if (sgd == 0.0) {
				// nrerror("toeplz-3 singular principal minor");
				// System.out.println("toeplz-3 singular principal minor");
				failS = "toeplz-3 singular principal minor";
				failB = true;

				return;
			}
			g[m1 - 1] = sgn / sgd;// g[m1]=sgn/sgd; whence G and H.
			h[m1 - 1] = shn / sd;// h[m1]=shn/sd;
			k = m - 1;// k=m;@@@@@@@@@@@@@@@@@@@@@@
			m2 = (m + 1) >> 1;
			pp = g[m1 - 1];// pp=g[m1];
			qq = h[m1 - 1];// qq=h[m1];
			for (j = 1; j <= m2; j++) {
				pt1 = g[j - 1];// pt1=g[j];
				pt2 = g[k];// pt2=g[k];@@@@@@@@@@@@@@@@@@@
				qt1 = h[j - 1];// qt1=h[j];
				qt2 = h[k];// qt2=h[k];@@@@@@@@@@@@@@@@@@@@@@@@@
				g[j - 1] = pt1 - pp * qt2;// g[j]=pt1-pp*qt2;
				g[k] = pt2 - pp * qt1;// g[k]=pt2-pp*qt1;@@@@@@@@@@@@@@@@@@@@@
				h[j - 1] = qt1 - qq * pt2;// h[j]=qt1-qq*pt2;
				h[k--] = qt2 - qq * pt1;// h[k--]=qt2-qq*pt1;@@@@@@@@@@@@@@@@@@
			}
		}// Back for another recurrence.

		// nrerror("toeplz - should not arrive here!");
		// System.out.println("toeplz - should not arrive here!");
		failS = "toeplz - should not arrive here!";
		failB = true;

		return;
	}

	// @@@@@@@Cholesky Decomposition
	/*
	 * If a square matrix A happens to be symmetric and positive definite, then
	 * it has a special, more efficient, triangular decomposition. Symmetric
	 * means that aij = aji for i, j = 1, . . . , N, while positive definite
	 * means that v · A · v > 0 for all vectors v
	 */
	/**
	 * Given a positive-definite symmetric matrix a[1..n][1..n], this routine 
	 * constructs its Cholesky decomposition, A = L · LT . On input, only the upper triangle of a need 
	 * be given; it is not modified. The Cholesky factor L is returned in the lower triangle of a, 
	 * except for its diagonal elements which are returned in p[1..n].
	 * @param a a
	 * @param n n
	 * @param p p
	 */
	public static void choldc(double[][] a, int n, double[] p)
	// Given a positive-definite symmetric matrix a[1..n][1..n], this routine
	// constructs its Cholesky
	// decomposition, A = L · LT . On input, only the upper triangle of a need
	// be given; it is not
	// modified. The Cholesky factor L is returned in the lower triangle of a,
	// except for its diagonal
	// elements which are returned in p[1..n].
	{
		// void nrerror(char error_text[]);
		int i = 0;
		int j = 0;
		int k = 0;
		double sum = 0.0;

		failB = false;
		for (i = 1; i <= n; i++) {
			for (j = i; j <= n; j++) {
				// for (sum=a[i][j],k=i-1;k>=1;k--) sum -= a[i][k]*a[j][k];
				for (sum = a[i - 1][j - 1], k = i - 1; k >= 1; k--)
					sum -= a[i - 1][k - 1] * a[j - 1][k - 1];
				if (i == j) {
					if (sum <= 0.0) // a, with rounding errors, is not positive
									// definite.
					{
						// nrerror("choldc failed");
						// System.out.println("choldc failed");
						failS = "choldc failed";
						failB = true;

						return;
					}
					p[i - 1] = Math.sqrt(sum);// p[i]=sqrt(sum);
				} else
					a[j - 1][i - 1] = sum / p[i - 1];// a[j][i]=sum/p[i];
			}
		}
	}

	/*
	 * You might at this point wonder about pivoting. The pleasant answer is
	 * that Cholesky decomposition is extremely stable numerically, without any
	 * pivoting at all. Failure of choldc simply indicates that the matrix A
	 * (or, with roundoff error, another very nearby matrix) is not positive
	 * definite. In fact, choldc is an efficient way to test whether a symmetric
	 * matrix is positive definite.
	 */
	/**
	 * Solves the set of n linear equations A · x = b, where a is a positive-definite symmetric matrix. 
	 * a[1..n][1..n] and p[1..n] are input as the output of the routine choldc. Only the lower 
	 * subdiagonal portion of a is accessed. b[1..n] is input as the right-hand side vector. The 
	 * solution vector is returned in x[1..n]. a, n, and p are not modified and can be left in place 
	 * for successive calls with different right-hand sides b. b is not modified unless you identify b and 
	 * x in the calling sequence, which is allowed.
	 * @param a a
	 * @param n n
	 * @param p p
	 * @param b b
	 * @param x x
	 */
	public static void cholsl(double[][] a, int n, double[] p, double[] b,
			double[] x)
	// Solves the set of n linear equations A · x = b, where a is a
	// positive-definite symmetric matrix.
	// a[1..n][1..n] and p[1..n] are input as the output of the routine choldc.
	// Only the lower
	// subdiagonal portion of a is accessed. b[1..n] is input as the right-hand
	// side vector. The
	// solution vector is returned in x[1..n]. a, n, and p are not modified and
	// can be left in place
	// for successive calls with different right-hand sides b. b is not modified
	// unless you identify b and
	// x in the calling sequence, which is allowed.
	{
		int i = 0;
		int k = 0;
		double sum = 0.0;
		for (i = 1; i <= n; i++) { // Solve L · y = b, storing y in x.
									// for (sum=b[i],k=i-1;k>=1;k--) sum -=
									// a[i][k]*x[k];
			for (sum = b[i - 1], k = i - 1; k >= 1; k--)
				sum -= a[i - 1][k - 1] * x[k - 1];
			x[i - 1] = sum / p[i - 1];// x[i]=sum/p[i];
		}
		for (i = n; i >= 1; i--) {// Solve LT · x = y.
									// for (sum=x[i],k=i+1;k<=n;k++) sum -=
									// a[k][i]*x[k];
			for (sum = x[i - 1], k = i + 1; k <= n; k++)
				sum -= a[k - 1][i - 1] * x[k - 1];
			x[i - 1] = sum / p[i - 1];// x[i]=sum/p[i];
		}
	}

	/*
	 * A typical use of choldc and cholsl is in the inversion of covariance
	 * matrices describing the fit of data to a model; see, e.g., §15.6. In
	 * this, and many other applications, one often needs L-1. The lower
	 * triangle of this matrix can be efficiently found from the output of
	 * choldc: for (i=1;i<=n;i++) { a[i][i]=1.0/p[i]; for (j=i+1;j<=n;j++) {
	 * sum=0.0; for (k=i;k<j;k++) sum -= a[j][k]*a[k][i]; a[j][i]=sum/p[j]; } }
	 */

	// @@@@QR decomposition,
	/*
	 * There is another matrix factorization that is sometimes very useful, the
	 * so-called QR decomposition, A = Q · R (2.10.1) Here R is upper
	 * triangular, while Q is orthogonal, that is, QT · Q = 1 (2.10.2) where QT
	 * is the transpose matrix of Q. Although the decomposition exists for a
	 * general rectangular matrix, we shall restrict our treatment to the case
	 * when all the matrices are square, with dimensions N × N. Like the other
	 * matrix factorizations we have met (LU, SVD, Cholesky), QR decomposition
	 * can be used to solve systems of linear equations. To solve A · x = b
	 * (2.10.3) first form QT · b and then solve R · x = QT · b (2.10.4) by
	 * backsubstitution. Since QR decomposition involves about twice as many
	 * operations as LU decomposition, it is not used for typical systems of
	 * linear equations. However, we will meet special cases where QR is the
	 * method of choice.
	 */
	public static int sing = 0;

	/**
	 * Constructs the QR decomposition of a[1..n][1..n]. The upper triangular matrix R is returned 
	 * in the upper triangle of a, except for the diagonal elements of R which are returned in 
	 * d[1..n]. The orthogonal matrix Q is represented as a product of n- 1 Householder matrices 
	 * Q1 . . .Qn-1, where Qj = 1-uj x uj/cj. The ith component of uj is zero for i = 1, . . . , j -1 
	 * while the nonzero components are returned in a[i][j] for i = j, . . . , n. sing returns as 
	 * true (1) if singularity is encountered during the decomposition, but the decomposition is still 
	 * completed in this case; otherwise it returns false (0).
	 * @param a a
	 * @param n n
	 * @param c c
	 * @param d d
	 */
	public static void qrdcmp(double[][] a, int n, double[] c, double[] d)// ,
																			// int
																			// *sing)
	// Constructs the QR decomposition of a[1..n][1..n]. The upper triangular
	// matrix R is returned
	// in the upper triangle of a, except for the diagonal elements of R which
	// are returned in
	// d[1..n]. The orthogonal matrix Q is represented as a product of n- 1
	// Householder matrices
	// Q1 . . .Qn-1, where Qj = 1-uj x uj/cj. The ith component of uj is zero
	// for i = 1, . . . , j -1
	// while the nonzero components are returned in a[i][j] for i = j, . . . ,
	// n. sing returns as
	// true (1) if singularity is encountered during the decomposition, but the
	// decomposition is still
	// completed in this case; otherwise it returns false (0).
	{
		int i = 0;
		int j = 0;
		int k = 0;
		double scale = 0.0;
		double sigma = 0.0;
		double sum = 0.0;
		double tau = 0.0;
		sing = 0;// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
		for (k = 1; k < n; k++) {
			scale = 0.0;
			for (i = k; i <= n; i++)
				scale = Math.max(scale, Math.abs(a[i - 1][k - 1]));// scale=FMAX(scale,fabs(a[i][k]));
			if (scale == 0.0) {// Singular case.
				sing = 1;
				c[k - 1] = d[k - 1] = 0.0;// c[k]=d[k]=0.0;
			} else { // Form Qk and Qk · A.
				for (i = k; i <= n; i++)
					a[i - 1][k - 1] /= scale;// a[i][k] /= scale;
				for (sum = 0.0, i = k; i <= n; i++)
					sum += (a[i - 1][k - 1]) * (a[i - 1][k - 1]);// sum +=
																	// SQR(a[i][k]);
				// sigma=SIGN(sqrt(sum),a[k][k]);
				// g = -SIGN(sqrt(s),f);//cap11: g=(f >= 0.0 ? -sqrt(h) :
				// sqrt(h));
				// g=(f >= 0.0 ? -Math.sqrt(s) : Math.sqrt(s));
				sigma = (a[k - 1][k - 1] >= 0.0 ? Math.sqrt(sum) : -Math
						.sqrt(sum));
				a[k - 1][k - 1] += sigma;// a[k][k] += sigma;
				c[k - 1] = sigma * a[k - 1][k - 1];// c[k]=sigma*a[k][k];
				d[k - 1] = -scale * sigma;// d[k] = -scale*sigma;
				for (j = k + 1; j <= n; j++) {
					for (sum = 0.0, i = k; i <= n; i++)
						sum += a[i - 1][k - 1] * a[i - 1][j - 1];// sum +=
																	// a[i][k]*a[i][j];
					tau = sum / c[k - 1];// tau=sum/c[k];
					for (i = k; i <= n; i++)
						a[i - 1][j - 1] -= tau * a[i - 1][k - 1];// a[i][j] -=
																	// tau*a[i][k];
				}
			}
		}
		d[n - 1] = a[n - 1][n - 1];// d[n]=a[n][n];
		if (d[n - 1] == 0.0)
			sing = 1;// if (d[n] == 0.0) *sing=1;
	}

	// The next routine, qrsolv, is used to solve linear systems.
	/**
	 * Solves the set of n linear equations A · x = b. a[1..n][1..n], c[1..n], and d[1..n] are 
	 * input as the output of the routine qrdcmp and are not modified. b[1..n] is input as the 
	 * right-hand side vector, and is overwritten with the solution vector on output.
	 * @param a a
	 * @param n n
	 * @param c c
	 * @param d d
	 * @param b b
	 */
	public static void qrsolv(double[][] a, int n, double[] c, double[] d,
			double[] b)
	// Solves the set of n linear equations A · x = b. a[1..n][1..n], c[1..n],
	// and d[1..n] are
	// input as the output of the routine qrdcmp and are not modified. b[1..n]
	// is input as the
	// right-hand side vector, and is overwritten with the solution vector on
	// output.
	{
		// void rsolv(float **a, int n, float d[], float b[]);
		int i = 0;
		int j = 0;
		double sum = 0.0;
		double tau = 0.0;
		for (j = 1; j < n; j++) { // Form QT · b.
			for (sum = 0.0, i = j; i <= n; i++)
				sum += a[i - 1][j - 1] * b[i - 1];// sum += a[i][j]*b[i];
			tau = sum / c[j - 1];// tau=sum/c[j];
			for (i = j; i <= n; i++)
				b[i - 1] -= tau * a[i - 1][j - 1];// b[i] -= tau*a[i][j];
		}
		rsolv(a, n, d, b); // Solve R · x = QT · b.
	}

	/**
	 * Solves the set of n linear equations R · x = b, where R is an upper triangular matrix stored in 
	 * a and d. a[1..n][1..n] and d[1..n] are input as the output of the routine qrdcmp and 
	 * are not modified. b[1..n] is input as the right-hand side vector, and is overwritten with the 
	 * solution vector on output.
	 * @param a a
	 * @param n n
	 * @param d d
	 * @param b b
	 */
	public static void rsolv(double[][] a, int n, double[] d, double[] b)
	// Solves the set of n linear equations R · x = b, where R is an upper
	// triangular matrix stored in
	// a and d. a[1..n][1..n] and d[1..n] are input as the output of the routine
	// qrdcmp and
	// are not modified. b[1..n] is input as the right-hand side vector, and is
	// overwritten with the
	// solution vector on output.
	{
		int i = 0;
		int j = 0;
		double sum = 0.0;
		b[n - 1] /= d[n - 1];// b[n] /= d[n];
		for (i = n - 1; i >= 1; i--) {
			for (sum = 0.0, j = i + 1; j <= n; j++)
				sum += a[i - 1][j - 1] * b[j - 1];// sum += a[i][j]*b[j];
			b[i - 1] = (b[i - 1] - sum) / d[i - 1];// b[i]=(b[i]-sum)/d[i];
		}
	}

	// @@@@@@Updating a QR decomposition
	/*
	 * Some numerical algorithms involve solving a succession of linear systems
	 * each of which differs only slightly from its predecessor. Instead of
	 * doing O(N3) operations each time to solve the equations from scratch, one
	 * can often update a matrix factorization in O(N2) operations and use the
	 * new factorization to solve the next set of linear equations. The LU
	 * decomposition is complicated to update because of pivoting. However, QR
	 * turns out to be quite simple for a very common kind of update, A › A + s
	 * x t
	 */
	/**
	 * Given matrices r[1..n][1..n] and qt[1..n][1..n], carry out a Jacobi rotation on rows 
	 * i and i + 1 of each matrix. a and b are the parameters of the rotation: 
	 * cos theta = a/sqrt(a2 + b2), sin theta = b/sqrt(a2 + b2).
	 * @param r r
	 * @param qt qt
	 * @param n n
	 * @param i i
	 * @param a a
	 * @param b b
	 */
	public static void rotate(double[][] r, double[][] qt, int n, int i,
			double a, double b)
	// Given matrices r[1..n][1..n] and qt[1..n][1..n], carry out a Jacobi
	// rotation on rows
	// i and i + 1 of each matrix. a and b are the parameters of the rotation:
	// cos theta = a/sqrt(a2 + b2),
	// sin theta = b/sqrt(a2 + b2).
	{
		int j = 0;
		double c = 0.0;
		double fact = 0.0;
		double s = 0.0;
		double w = 0.0;
		double y = 0.0;
		if (a == 0.0) {// Avoid unnecessary overflow or underflow.
			c = 0.0;
			s = (b >= 0.0 ? 1.0 : -1.0);
		} else if (Math.abs(a) > Math.abs(b))// (fabs(a) > fabs(b))
		{
			fact = b / a;
			// c=SIGN(1.0/sqrt(1.0+(fact*fact)),a);
			// g=(f >= 0.0 ? -Math.sqrt(s) : Math.sqrt(s));
			c = (a >= 0.0 ? 1.0 / Math.sqrt(1.0 + (fact * fact)) : -1.0
					/ Math.sqrt(1.0 + (fact * fact)));
			s = fact * c;
		} else {
			fact = a / b;
			// s=SIGN(1.0/sqrt(1.0+(fact*fact)),b);
			s = (b >= 0.0 ? 1.0 / Math.sqrt(1.0 + (fact * fact)) : -1.0
					/ Math.sqrt(1.0 + (fact * fact)));
			c = fact * s;
		}
		for (j = i; j <= n; j++) { // Premultiply r by Jacobi rotation.
			y = r[i - 1][j - 1];// y=r[i][j];
			w = r[i][j - 1];// w=r[i+1][j];
			r[i - 1][j - 1] = c * y - s * w;// r[i][j]=c*y-s*w;
			r[i][j - 1] = s * y + c * w;// r[i+1][j]=s*y+c*w;
		}
		for (j = 1; j <= n; j++) {// Premultiply qt by Jacobi rotation.
			y = qt[i - 1][j - 1];// y=qt[i][j];
			w = qt[i][j - 1];// w=qt[i+1][j];
			qt[i - 1][j - 1] = c * y - s * w;// qt[i][j]=c*y-s*w;
			qt[i][j - 1] = s * y + c * w;// qt[i+1][j]=s*y+c*w;
		}
	}

	/**
	 * Given the QR decomposition of some n × n matrix, calculates the QR decomposition of the 
	 * matrix Q·(R+u x v). The quantities are dimensioned as r[1..n][1..n], qt[1..n][1..n], 
	 * u[1..n], and v[1..n]. Note that QT is input and returned in qt.
	 * @param r r
	 * @param qt qt
	 * @param n n
	 * @param u u
	 * @param v v
	 */
	public static void qrupdt(double[][] r, double[][] qt, int n, double[] u,
			double[] v)
	// Given the QR decomposition of some n × n matrix, calculates the QR
	// decomposition of the
	// matrix Q·(R+u x v). The quantities are dimensioned as r[1..n][1..n],
	// qt[1..n][1..n],
	// u[1..n], and v[1..n]. Note that QT is input and returned in qt.
	{
		// void rotate(float **r, float **qt, int n, int i, float a, float b);
		int i = 0;
		int j = 0;
		int k = 0;
		for (k = n; k >= 1; k--) {// Find largest k such that u[k] != 0.
			if (u[k - 1] != 0.0)
				break;// if (u[k]) break;<=>if (w[j-1]!=0.0)//if (w[j])
		}
		if (k < 1)
			k = 1;
		for (i = k - 1; i >= 1; i--) {// Transform R + u x v to upper
										// Hessenberg.
			rotate(r, qt, n, i, u[i - 1], -u[i]);// rotate(r,qt,n,i,u[i],-u[i+1]);
			if (u[i - 1] == 0.0)
				u[i - 1] = Math.abs(u[i]);// if (u[i] == 0.0) u[i]=fabs(u[i+1]);
			else if (Math.abs(u[i - 1]) > Math.abs(u[i]))// (fabs(u[i]) >
															// fabs(u[i+1]))
				u[i - 1] = Math.abs(u[i - 1])
						* Math.sqrt(1.0 + (u[i] / u[i - 1]) * (u[i] / u[i - 1]));// u[i]=fabs(u[i])*sqrt(1.0+SQR(u[i+1]/u[i]));
			else
				u[i - 1] = Math.abs(u[i])
						* Math.sqrt(1.0 + (u[i - 1] / u[i]) * (u[i - 1] / u[i]));// u[i]=fabs(u[i+1])*sqrt(1.0+SQR(u[i]/u[i+1]));
		}
		for (j = 1; j <= n; j++)
			r[0][j - 1] += u[0] * v[j - 1];// r[1][j] += u[1]*v[j];
		for (i = 1; i < k; i++)
			// Transform upper Hessenberg matrix to upper triangular.
			rotate(r, qt, n, i, r[i - 1][i - 1], -r[i][i - 1]);// rotate(r,qt,n,i,r[i][i],-r[i+1][i]);
	}
	//==================
	//--rezolvarea sistemelor de ecuatii liniare prin metoda substitutiei Gauss
	//n=numar linii, m=numar coloane --pentru 1 sistemec =1!!!
	//a=matricea coeficientilor sistemului liniar b= matricea termenilor liberi
	/**
	 * Gaussian elimination with backsubstitution. The input data are the same as in gaussj. 
	 * This is gauss above...with different name. Must be a left-over.
	 * @param a a
	 * @param n n
	 * @param b b
	 * @param m m
	 * @return the determinant of a.
	 */
	  	public static double sysEqGauss(double[][] a, double[][] b, int n, int m)
	  	{
	        int l=0;
	        double t = 0.0;
	        double det = 1.0;

	        //cautare element pivot nenul
	        for (int k = 0; k<n; k++)
	        {
				t = 0.0;
				for (int i = k; i<n; i++)
				   if (t<Math.abs(a[i][k]))
				   {
					   t=Math.abs(a[i][k]);
					   l=i;
				   }
				//interschimb liniile l si k daca e cazul<->pun pivotul pe diagonala
				if(l!=k)
				{
					det = -det;//determinantul isi schimba semnul
	                for (int j = k; j<n; j++)
	                {
						t=a[k][j];
						a[k][j]=a[l][j];
	                    a[l][j]=t;
					}
					for (int j =0; j<m; j++)
					{
					    t=b[k][j];
	                    b[k][j]=b[l][j];
	                    b[l][j]=t;
					}
				}
				///gata cu pivotul
				////impart linia pivot
				t=1/a[k][k];
				det=det*a[k][k];
			    for (int j=k+1; j<n; j++)
				   a[k][j]=a[k][j]*t;
				for (int j=0; j<m; j++)
				   b[k][j]=b[k][j]*t;
			   //reduc liniile nepivot
			    for (int i=k+1; i<n; i++)
				{
			        t=a[i][k];
			        for (int j=k+1; j<n; j++)
			          a[i][j]=a[i][j]-a[k][j]*t;
			        for (int j=0; j<m; j++)
			           b[i][j]=b[i][j]-b[k][j]*t;
			    }
			}
			///construiesc tabloul cu rezultatul de interes
			for (int k=n-2; k>=0; k--)//mai subtil merge si asa->(int k=n-1; --k>=0; )//
			    for (int j=0; j<m; j++)
			    {
			        t=0;
			        for (int i=k+1; i<n; i++)
			          t=t+a[k][i]*b[i][j];
			        b[k][j]=b[k][j]-t;
			    }

	        return det;
		}

		//rezolvarea sistemelor de ecuatii ce prezinta matrice tridiagonala
		//a,b,c sunt vectorii celor 3 "diagonale" ale matricii
		//n=numar de linii
	  	/**
		 * Alternative approach for tridiagonal matrices. Solution is stored in d. 
		 * This is the same as triDiag....must be some left-over.
		 * @param a a
		 * @param b b
		 * @param c c
		 * @param d d
		 * @param n n
		 */
		public static void sysEqTriDiag(double[] a, double[] b, double[] c, double[] d, int n)
	  	{
			double t=0.0;
			c[0]=c[0]/b[0];
			d[0]=d[0]/b[0];
			for (int i=1; i<=n-2;i++)
			{
			    t=b[i]-a[i]*c[i-1];
			    c[i]=c[i]/t;
			    d[i]=(d[i]-a[i]*d[i-1])/t;
		    }
		    d[n-1]=(d[n-1]-a[n-1]*d[n-2])/(b[n-1]-a[n-1]*c[n-2]);
			for(int i=n-2; i>= 0; i--)
			  d[i]=d[i]-c[i]*d[i+1];

		}
}
