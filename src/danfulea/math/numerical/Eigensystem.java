package danfulea.math.numerical;

/**
 * Eigensystem class. 
 * Based on literature, including Numerical Recipes in C (Cambridge Univ.).
 * @author Dan Fulea, 13 OCT. 2006
 */
public class Eigensystem {
	public static boolean failB = false;
	public static String failS = "";

	public static int nrot_jacobi = 0;
	public static double RADIX = 2.0;

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
	 * An N × N matrix A is said to have an eigenvector x and corresponding
	 * eigenvalue lambda if A · x = lambda x=>vectori si valori proprii
	 * 
	 * Obviously any multiple of an eigenvector x will also be an eigenvector,
	 * but we won’t consider such multiples as being distinct eigenvectors. (The
	 * zero vector is not considered to be an eigenvector at all.) Evidently
	 * (11.0.1) can hold only if det |A - lambda1| = 0
	 * 
	 * 
	 * Jacobi Transformations of a Symmetric Matrix
	 * 
	 * In the following routine the n×n symmetric matrix a is stored as a[1..n]
	 * [1..n]. On output, the superdiagonal elements of a are destroyed, but the
	 * diagonal and subdiagonal are unchanged and give full information on the
	 * original symmetric matrix a. The vector d[1..n] returns the eigenvalues
	 * of a. During the computation, it contains the current diagonal of a. The
	 * matrix v[1..n][1..n] outputs the normalized eigenvector belonging to d[k]
	 * in its kth column. The parameter nrot is the number of Jacobi rotations
	 * that were needed to achieve convergence. #define ROTATE(a,i,j,k,l)
	 * g=a[i][j];h=a[k][l];a[i][j]=g-s*(h+g*tau);\ a[k][l]=h+s*(g-h*tau);
	 */

	/**
	 * Jacobi Transformations of a Symmetric Matrix
	 * @param a n×n symmetric matrix
	 * @param n dimension of matrix a
	 * @param d returns the eigenvalues of a.
	 * @param v outputs the normalized eigenvector belonging to d[k] in its kth column
	 */
	public static void jacobi(double[][] a, int n, double[] d, double[][] v)// ,
																			// int
																			// *nrot)
	// Computes all eigenvalues and eigenvectors of a real symmetric matrix
	// a[1..n][1..n]. On
	// output, elements of a above the diagonal are destroyed. d[1..n] returns
	// the eigenvalues of a.
	// v[1..n][1..n] is a matrix whose columns contain, on output, the
	// normalized eigenvectors of
	// a. nrot returns the number of Jacobi rotations that were required.
	{
		int j = 0;
		int iq = 0;
		int ip = 0;
		int i = 0;
		double tresh = 0.0;
		double theta = 0.0;
		double tau = 0.0;
		double t = 0.0;
		double sm = 0.0;
		double s = 0.0;
		double h = 0.0;
		double g = 0.0;
		double c = 0.0;
		double[] b = new double[n];
		double[] z = new double[n];

		failB = false;
		// b=vector(1,n);
		// z=vector(1,n);
		for (ip = 1; ip <= n; ip++) {// Initialize to the identity matrix.
			for (iq = 1; iq <= n; iq++)
				v[ip - 1][iq - 1] = 0.0;// v[ip][iq]=0.0;
			v[ip - 1][ip - 1] = 1.0;// v[ip][ip]=1.0;
		}
		for (ip = 1; ip <= n; ip++) {// Initialize b and d to the diagonalof a.
			b[ip - 1] = d[ip - 1] = a[ip - 1][ip - 1];// b[ip]=d[ip]=a[ip][ip];
			z[ip - 1] = 0.0;// z[ip]=0.0; This vector will accumulate terms of
							// the form tapq as in equation(11.1.14).
		}
		nrot_jacobi = 0;
		for (i = 1; i <= 50; i++) {
			sm = 0.0;
			for (ip = 1; ip <= n - 1; ip++) {// Sum off-diagonal elements.
				for (iq = ip + 1; iq <= n; iq++)
					sm += Math.abs(a[ip - 1][iq - 1]);// sm += fabs(a[ip][iq]);
			}
			if (sm == 0.0) {// The normal return, which relies on quadratic
							// convergence to machine underflow.
							// free_vector(z,1,n);
							// free_vector(b,1,n);
				return;
			}
			if (i < 4)
				tresh = 0.2 * sm / (n * n);// ...on the first three sweeps.
			else
				tresh = 0.0;// ...thereafter.
			for (ip = 1; ip <= n - 1; ip++) {
				for (iq = ip + 1; iq <= n; iq++) {
					g = 100.0 * Math.abs(a[ip - 1][iq - 1]);// g=100.0*fabs(a[ip][iq]);
					// After four sweeps, skip the rotation if the off-diagonal
					// element is small.
					// if (i > 4 && (float)(fabs(d[ip])+g) == (float)fabs(d[ip])
					if (i > 4
							&& (double) (Math.abs(d[ip - 1]) + g) == (double) Math
									.abs(d[ip - 1])
							// && (float)(fabs(d[iq])+g) == (float)fabs(d[iq]))
							&& (double) (Math.abs(d[iq - 1]) + g) == (double) Math
									.abs(d[iq - 1]))
						a[ip - 1][iq - 1] = 0.0;// a[ip][iq]=0.0;
					else if (Math.abs(a[ip - 1][iq - 1]) > tresh) // (fabs(a[ip][iq])
																	// > tresh)
					{
						h = d[iq - 1] - d[ip - 1];// h=d[iq]-d[ip];
						if ((double) (Math.abs(h) + g) == (double) Math.abs(h))
							t = (a[ip - 1][iq - 1]) / h;// t=(a[ip][iq])/h; t =
														// 1/(2?)
						else {
							theta = 0.5 * h / (a[ip - 1][iq - 1]);// theta=0.5*h/(a[ip][iq]);
																	// Equation
																	// (11.1.10).
							t = 1.0 / (Math.abs(theta) + Math.sqrt(1.0 + theta
									* theta));
							if (theta < 0.0)
								t = -t;
						}
						c = 1.0 / Math.sqrt(1 + t * t);
						s = t * c;
						tau = s / (1.0 + c);
						h = t * a[ip - 1][iq - 1];// h=t*a[ip][iq];
						z[ip - 1] -= h;// z[ip] -= h;
						z[iq - 1] += h;// z[iq] += h;
						d[ip - 1] -= h;// d[ip] -= h;
						d[iq - 1] += h;// d[iq] += h;
						a[ip - 1][iq - 1] = 0.0;// a[ip][iq]=0.0;
						for (j = 1; j <= ip - 1; j++) {// Case of rotations 1 =
														// j <
														// p.ROTATE(a,i,j,k,l)
														// ROTATE(a,j,ip,j,iq)
														// g=a[i][j];h=a[k][l];a[i][j]=g-s*(h+g*tau);
														// a[k][l]=h+s*(g-h*tau);
							g = a[j - 1][ip - 1];
							h = a[j - 1][iq - 1];
							a[j - 1][ip - 1] = g - s * (h + g * tau);
							a[j - 1][iq - 1] = h + s * (g - h * tau);
						}
						for (j = ip + 1; j <= iq - 1; j++) {// Case of rotations
															// p < j < q.
															// ROTATE(a,ip,j,j,iq)
							g = a[ip - 1][j - 1];
							h = a[j - 1][iq - 1];
							a[ip - 1][j - 1] = g - s * (h + g * tau);
							a[j - 1][iq - 1] = h + s * (g - h * tau);
						}
						for (j = iq + 1; j <= n; j++) {// Case of rotations q <
														// j = n.
														// ROTATE(a,ip,j,iq,j)
							g = a[ip - 1][j - 1];
							h = a[iq - 1][j - 1];
							a[ip - 1][j - 1] = g - s * (h + g * tau);
							a[iq - 1][j - 1] = h + s * (g - h * tau);
						}
						for (j = 1; j <= n; j++) {
							// ROTATE(v,j,ip,j,iq)
							g = v[j - 1][ip - 1];
							h = v[j - 1][iq - 1];
							v[j - 1][ip - 1] = g - s * (h + g * tau);
							v[j - 1][iq - 1] = h + s * (g - h * tau);
						}
						++(nrot_jacobi);
					}
				}
			}
			for (ip = 1; ip <= n; ip++) {
				b[ip - 1] += z[ip - 1];// b[ip] += z[ip];
				d[ip - 1] = b[ip - 1];// d[ip]=b[ip]; //Update d with the sum of
										// tapq,
				z[ip - 1] = 0.0;// z[ip]=0.0; //and reinitialize z.
			}
		}

		// nrerror("Too many iterations in routine jacobi");
		failB = true;
		failS = "Too many iterations in routine jacobi";
		return;
	}

	/*
	 * Note that the above routine assumes that underflows are set to zero. On
	 * machines where this is not true, the program must be modified. The
	 * eigenvalues are not ordered on output. If sorting is desired, the
	 * following routine can be invoked to reorder the output of jacobi or of
	 * later routines in this chapter. (The method, straight insertion, is N 2
	 * rather than N logN; but since you have just done an N3 procedure to get
	 * the eigenvalues, you can afford yourself this little indulgence.)
	 */

	/**
	 * Given the eigenvalues d[1..n] and eigenvectors v[1..n][1..n] as output from e.g. jacobi, this routine sorts the eigenvalues into
	 * descending order, and rearranges the columns of v correspondingly.
	 * @param d eigenvalues computed by e.g. jacobi
	 * @param v their corresponding eigenvectors
	 * @param n the dimension
	 */
	public static void eigsrt(double[] d, double[][] v, int n)
	// Given the eigenvalues d[1..n] and eigenvectors v[1..n][1..n] as output
	// from jacobi
	// (§11.1) or tqli (§11.3), this routine sorts the eigenvalues into
	// descending order, and rearranges
	// the columns of v correspondingly. The method is straight insertion.
	{
		int k = 0;
		int j = 0;
		int i = 0;
		double p = 0.0;

		for (i = 1; i < n; i++) {
			k = i;
			p = d[k - 1];// p=d[k=i];

			for (j = i + 1; j <= n; j++) {
				if (d[j - 1] >= p) {
					k = j;
					p = d[k - 1];// if (d[j] >= p) p=d[k=j];
				}
			}
			if (k != i) {
				d[k - 1] = d[i - 1];// d[k]=d[i];
				d[i - 1] = p;// d[i]=p;
				for (j = 1; j <= n; j++) {
					p = v[j - 1][i - 1];// p=v[j][i];
					v[j - 1][i - 1] = v[j - 1][k - 1];// v[j][i]=v[j][k];
					v[j - 1][k - 1] = p;// v[j][k]=p;
				}
			}
		}
	}

	/*
	 * Reduction of a Symmetric Matrix to Tridiagonal Form: Givens and
	 * Householder Reductions
	 * 
	 * Note that when dealing with a matrix whose elements vary over many orders
	 * of magnitude, it is important that the matrix be permuted, insofar as
	 * possible, so that the smaller elements are in the top left-hand corner.
	 * This is because the reduction is performed starting from the bottom
	 * right-hand corner, and a mixture of small and large elements there can
	 * lead to considerable rounding errors. The routine tred2 is designed for
	 * use with the routine tqli of the next section. tqli finds the eigenvalues
	 * and eigenvectors of a symmetric, tridiagonal matrix. The combination of
	 * tred2 and tqli is the most efficient known technique for finding all the
	 * eigenvalues and eigenvectors (or just all the eigenvalues) of a real,
	 * symmetric matrix. In the listing below, the statements indicated by
	 * comments are required only for subsequent computation of eigenvectors. If
	 * only eigenvalues are required, omission of the commented statements
	 * speeds up the execution time of tred2 by a factor of 2 for large n. In
	 */
	/**
	 * Householder reduction of a real, symmetric matrix a[1..n][1..n].
	 * On output, a is replaced by the orthogonal matrix Q effecting the transformation. d[1..n] returns 
	 * the diagonal elements of the tridiagonal matrix, and e[1..n] the off-diagonal elements, with 
	 * e[1]=0.
	 * @param a a
	 * @param n dimension
	 * @param d d
	 * @param e e
	 */
	public static void tred2(double[][] a, int n, double[] d, double[] e)
	// Householder reduction of a real, symmetric matrix a[1..n][1..n]. On
	// output, a is replaced
	// by the orthogonal matrix Q effecting the transformation. d[1..n] returns
	// the diagonal elements
	// of the tridiagonal matrix, and e[1..n] the off-diagonal elements, with
	// e[1]=0. Several
	// statements, as noted in comments, can be omitted if only eigenvalues are
	// to be found, in which
	// case a contains no useful information on output. Otherwise they are to be
	// included.
	{
		int l = 0;
		int k = 0;
		int j = 0;
		int i = 0;
		double scale = 0.0;
		double hh = 0.0;
		double h = 0.0;
		double g = 0.0;
		double f = 0.0;

		for (i = n; i >= 2; i--) {
			l = i - 1;
			h = scale = 0.0;
			if (l > 1) {
				for (k = 1; k <= l; k++)
					scale += Math.abs(a[i - 1][k - 1]);// scale +=
														// fabs(a[i][k]);
				if (scale == 0.0) // Skip transformation.
					e[i - 1] = a[i - 1][l - 1];// e[i]=a[i][l];
				else {
					for (k = 1; k <= l; k++) {
						a[i - 1][k - 1] /= scale;// a[i][k] /= scale; Use scaled
													// a’s for transformation.
						h += a[i - 1][k - 1] * a[i - 1][k - 1];// h +=
																// a[i][k]*a[i][k];
																// Form s in h.
					}
					f = a[i - 1][l - 1];// f=a[i][l];
					g = (f >= 0.0 ? -Math.sqrt(h) : Math.sqrt(h));
					e[i - 1] = scale * g;// e[i]=scale*g;
					h -= f * g;// Now h is equation (11.2.4).
					a[i - 1][l - 1] = f - g;// a[i][l]=f-g; Store u in the ith
											// row of a.
					f = 0.0;
					for (j = 1; j <= l; j++) {
						/*
						 * Next statement can be omitted if eigenvectors not
						 * wanted
						 */
						a[j - 1][i - 1] = a[i - 1][j - 1] / h;// a[j][i]=a[i][j]/h;
																// Store u/H in
																// ith column of
																// a.
						g = 0.0;// Form an element of A · u in g.
						for (k = 1; k <= j; k++)
							g += a[j - 1][k - 1] * a[i - 1][k - 1];// g +=
																	// a[j][k]*a[i][k];
						for (k = j + 1; k <= l; k++)
							g += a[k - 1][j - 1] * a[i - 1][k - 1];// g +=
																	// a[k][j]*a[i][k];
						e[j - 1] = g / h;// e[j]=g/h;
						f += e[j - 1] * a[i - 1][j - 1];// f += e[j]*a[i][j];
					}
					hh = f / (h + h);// Form K, equation (11.2.11).
					for (j = 1; j <= l; j++) {// Form q and store in e
												// overwriting p.
						f = a[i - 1][j - 1];// f=a[i][j];
						e[j - 1] = g = e[j - 1] - hh * f;// e[j]=g=e[j]-hh*f;
						for (k = 1; k <= j; k++)
							// Reduce a, equation (11.2.13).
							a[j - 1][k - 1] -= (f * e[k - 1] + g
									* a[i - 1][k - 1]);// a[j][k] -=
														// (f*e[k]+g*a[i][k]);
					}
				}
			} else
				e[i - 1] = a[i - 1][l - 1];// e[i]=a[i][l];
			d[i - 1] = h;// d[i]=h;
		}
		/* Next statement can be omitted if eigenvectors not wanted */
		d[0] = 0.0;// d[1]=0.0;
		e[0] = 0.0;// e[1]=0.0;
		/*
		 * Contents of this loop can be omitted if eigenvectors not wanted
		 * except for statement d[i]=a[i][i];
		 */
		for (i = 1; i <= n; i++) {// Begin accumulation of transformation
									// matrices.
			l = i - 1;
			if (d[i - 1] != 0.0)// (d[i])
			{// This block skipped when i=1.
				for (j = 1; j <= l; j++) {
					g = 0.0;
					for (k = 1; k <= l; k++)
						// Use u and u/H stored in a to form P·Q.
						g += a[i - 1][k - 1] * a[k - 1][j - 1];// g +=
																// a[i][k]*a[k][j];
					for (k = 1; k <= l; k++)
						a[k - 1][j - 1] -= g * a[k - 1][i - 1];// a[k][j] -=
																// g*a[k][i];
				}
			}
			d[i - 1] = a[i - 1][i - 1];// d[i]=a[i][i]; This statement remains.
			a[i - 1][i - 1] = 1.0;// a[i][i]=1.0; Reset row and column of a to
									// identity
			// matrix for next iteration.
			for (j = 1; j <= l; j++)
				a[j - 1][i - 1] = a[i - 1][j - 1] = 0.0;// a[j][i]=a[i][j]=0.0;
		}
	}

	/*
	 * Eigenvalues and Eigenvectors of a Tridiagonal Matrix
	 */
	/**
	 * QL algorithm with implicit shifts, to determine the eigenvalues and eigenvectors of a real, symmetric, 
	 * tridiagonal matrix, or of a real, symmetric matrix previously reduced by tred2. 
	 * On input, d[1..n] contains the diagonal elements of the tridiagonal matrix. On output, it returns 
	 * the eigenvalues. The vector e[1..n] inputs the subdiagonal elements of the tridiagonal matrix, 
	 * with e[1] arbitrary. On output e is destroyed. If the eigenvectors of a tridiagonal matrix are desired, 
	 * the matrix z[1..n][1..n] is input as the identity matrix. If the eigenvectors of a matrix 
	 * that has been reduced by tred2 are required, then z is input as the matrix output by tred2. 
	 * In either case, the kth column of z returns the normalized eigenvector corresponding to d[k].
	 * @param d d
	 * @param e e
	 * @param n n
	 * @param z z
	 */

	public static void tqli(double[] d, double[] e, int n, double[][] z)
	// QL algorithm with implicit shifts, to determine the eigenvalues and
	// eigenvectors of a real, symmetric,
	// tridiagonal matrix, or of a real, symmetric matrix previously reduced by
	// tred2 §11.2. On
	// input, d[1..n] contains the diagonal elements of the tridiagonal matrix.
	// On output, it returns
	// the eigenvalues. The vector e[1..n] inputs the subdiagonal elements of
	// the tridiagonal matrix,
	// with e[1] arbitrary. On output e is destroyed. When finding only the
	// eigenvalues, several lines
	// may be omitted, as noted in the comments. If the eigenvectors of a
	// tridiagonal matrix are desired,
	// the matrix z[1..n][1..n] is input as the identity matrix. If the
	// eigenvectors of a matrix
	// that has been reduced by tred2 are required, then z is input as the
	// matrix output by tred2.
	// In either case, the kth column of z returns the normalized eigenvector
	// corresponding to d[k].
	{
		// float pythag(float a, float b);=>LinAEq.pythag(double a, double b)
		int m = 0;
		int l = 0;
		int iter = 0;
		int i = 0;
		int k = 0;
		double s = 0.0;
		double r = 0.0;
		double p = 0.0;
		double g = 0.0;
		double f = 0.0;
		double dd = 0.0;
		double c = 0.0;
		double b = 0.0;

		failB = false;

		for (i = 2; i <= n; i++)
			e[i - 2] = e[i - 1];// e[i-1]=e[i]; //Convenient to renumber the
								// elements of e.
		e[n - 1] = 0.0;// e[n]=0.0;
		for (l = 1; l <= n; l++) {
			iter = 0;
			do {
				for (m = l; m <= n - 1; m++) {// Look for a single small
												// subdiagonal
												// element to split the matrix.
					dd = Math.abs(d[m - 1]) + Math.abs(d[m]);// dd=fabs(d[m])+fabs(d[m+1]);
					if (Math.abs(e[m - 1]) + dd == dd)
						break;// if ((float)(fabs(e[m])+dd) == dd) break;
				}
				if (m != l) {
					if (iter++ == 30) {
						// nrerror("Too many iterations in tqli");
						failS = "too many iterations in gaujac";
						failB = true;
						return;
					}
					g = (d[l] - d[l - 1]) / (2.0 * e[l - 1]);// g=(d[l+1]-d[l])/(2.0*e[l]);
																// Form shift.
					r = LinAEq.pythag(g, 1.0);
					// g=d[m]-d[l]+e[l]/(g+SIGN(r,g)); //This is dm - ks.

					double dbl = (g >= 0.0 ? r : -r);
					g = d[m - 1] - d[l - 1] + e[l - 1] / (g + dbl);
					s = c = 1.0;
					p = 0.0;
					for (i = m - 1; i >= l; i--) {// A plane rotation as in the
													// original
													// QL, followed by Givens
													// rotations to restore
													// tridiagonal form.
						f = s * e[i - 1];// f=s*e[i];
						b = c * e[i - 1];// b=c*e[i];
						e[i] = (r = LinAEq.pythag(f, g));// e[i+1]=(r=pythag(f,g));
						if (r == 0.0) {// Recover from underflow.
							d[i] -= p;// d[i+1] -= p;
							e[m - 1] = 0.0;// e[m]=0.0;
							break;
						}
						s = f / r;
						c = g / r;
						g = d[i] - p;// g=d[i+1]-p;
						r = (d[i - 1] - g) * s + 2.0 * c * b;// r=(d[i]-g)*s+2.0*c*b;
						d[i] = g + (p = s * r);// d[i+1]=g+(p=s*r);
						g = c * r - b;
						// Next loop can be omitted if eigenvectors not wanted
						for (k = 1; k <= n; k++) {// Form eigenvectors.
							f = z[k - 1][i];// f=z[k][i+1];
							z[k - 1][i] = s * z[k - 1][i - 1] + c * f;// z[k][i+1]=s*z[k][i]+c*f;
							z[k - 1][i - 1] = c * z[k - 1][i - 1] - s * f;// z[k][i]=c*z[k][i]-s*f;
						}
					}
					if (r == 0.0 && i >= l)
						continue;
					d[l - 1] -= p;// d[l] -= p;
					e[l - 1] = g;// e[l]=g;
					e[m - 1] = 0.0;// e[m]=0.0;
				}
			} while (m != l);
		}
	}

	/*
	 * Hermitian Matrices The complex analog of a real, symmetric matrix is a
	 * Hermitian matrix, satisfying equation (11.0.4). Jacobi transformations
	 * can be used to find eigenvalues and eigenvectors, as also can Householder
	 * reduction to tridiagonal form followed by QL iteration. Complex versions
	 * of the previous routines jacobi, tred2, and tqli are quite analogous to
	 * their real counterparts.
	 * 
	 * Reduction of a General Matrix to Hessenberg Form The algorithms for
	 * symmetric matrices, given in the preceding sections, are highly
	 * satisfactory in practice. By contrast, it is impossible to design equally
	 * satisfactory algorithms for the nonsymmetric case. There are two reasons
	 * for this. First, the eigenvalues of a nonsymmetric matrix can be very
	 * sensitive to small changes in the matrix elements. Second, the matrix
	 * itself can be defective, so that there is no complete set of
	 * eigenvectors. We emphasize that these difficulties are intrinsic
	 * properties of certain nonsymmetric matrices, and no numerical procedure
	 * can “cure” them. The best we can hope for are procedures that don’t
	 * exacerbate such problems. The presence of rounding error can only make
	 * the situation worse. With finiteprecision arithmetic, one cannot even
	 * design a foolproof algorithm to determine whether a given matrix is
	 * defective or not. Thus current algorithms generally try to find a
	 * complete set of eigenvectors, and rely on the user to inspect the
	 * results. If any eigenvectors are almost parallel, the matrix is probably
	 * defective.
	 * 
	 * The sensitivity of eigenvalues to rounding errors during the execution of
	 * some algorithms can be reduced by the procedure of balancing. The errors
	 * in the eigensystem found by a numerical procedure are generally
	 * proportional to the Euclidean norm of the matrix, that is, to the square
	 * root of the sum of the squares of the elements. The idea of balancing is
	 * to use similarity transformations to make corresponding rows and columns
	 * of the matrix have comparable norms, thus reducing the overall norm of
	 * the matrix while leaving the eigenvalues unchanged.@@@@@@@@@@ A symmetric
	 * matrix is already balanced. Balancing is a procedure with of order N2
	 * operations. Thus, the time taken by the procedure balanc, given below,
	 * should never be more than a few percent of the total time required to
	 * find the eigenvalues. It is therefore recommended that you always balance
	 * nonsymmetric matrices. It never hurts, and it can substantially improve
	 * the accuracy of the eigenvalues computed for a badly balanced matrix.
	 */
	/**
	 * Given a matrix a[1..n][1..n], this routine replaces it by a balanced 
	 * matrix with identical eigenvalues. A symmetric matrix is already balanced and is unaffected by 
	 * this procedure. Class member RADIX should be the machine’s floating-point radix.
	 * @param a a
	 * @param n n
	 */
	public static void balanc(double[][] a, int n)
	// Given a matrix a[1..n][1..n], this routine replaces it by a balanced
	// matrix with identical
	// eigenvalues. A symmetric matrix is already balanced and is unaffected by
	// this procedure. The
	// parameter RADIX should be the machine’s floating-point radix.
	{
		int last = 0;
		int j = 0;
		int i = 0;
		double s = 0.0;
		double r = 0.0;
		double g = 0.0;
		double f = 0.0;
		double c = 0.0;
		double sqrdx = RADIX * RADIX;

		last = 0;
		while (last == 0) {
			last = 1;
			for (i = 1; i <= n; i++) {// Calculate row and column norms.
				r = c = 0.0;
				for (j = 1; j <= n; j++)
					if (j != i) {
						c += Math.abs(a[j - 1][i - 1]);// Math.abs(a[j][i]);
						r += Math.abs(a[i - 1][j - 1]);// Math.abs(a[i][j]);
					}
				if (c != 0.0 && r != 0.0)// if (c && r)
				{// If both are nonzero,
					g = r / RADIX;
					f = 1.0;
					s = c + r;
					while (c < g) {// find the integer power of the machine
									// radix that
									// comes closest to balancing the matrix.
						f *= RADIX;
						c *= sqrdx;
					}
					g = r * RADIX;
					while (c > g) {
						f /= RADIX;
						c /= sqrdx;
					}
					if ((c + r) / f < 0.95 * s) {
						last = 0;
						g = 1.0 / f;
						for (j = 1; j <= n; j++)
							a[i - 1][j - 1] *= g;// a[i][j] *= g; Apply
													// similarity
													// transformation.
						for (j = 1; j <= n; j++)
							a[j - 1][i - 1] *= f;// a[j][i] *= f;
					}
				}
			}
		}
	}

	/**
	 * Reduction to Hessenberg form by the elimination method. The real, nonsymmetric matrix 
	 * a[1..n][1..n] is replaced by an upper Hessenberg matrix with identical eigenvalues. Recommended, 
	 * but not required, is that this routine be preceded by balanc. On output the 
	 * Hessenberg matrix is in elements a[i][j] with i = j+1. Elements with i greater than j+1 are to be 
	 * thought of as zero, but are returned with random values.
	 * @param a a
	 * @param n n
	 */
	public static void elmhes(double[][] a, int n)
	// Reduction to Hessenberg form by the elimination method. The real,
	// nonsymmetric matrix
	// a[1..n][1..n] is replaced by an upper Hessenberg matrix with identical
	// eigenvalues. Recommended,
	// but not required, is that this routine be preceded by balanc. On output,
	// the
	// Hessenberg matrix is in elements a[i][j] with i = j+1. Elements with i >
	// j+1 are to be
	// thought of as zero, but are returned with random values.
	{
		int m = 0;
		int j = 0;
		int i = 0;
		double y = 0.0;
		double x = 0.0;

		for (m = 2; m < n; m++) {// m is called r + 1 in the text.
			x = 0.0;
			i = m;
			for (j = m; j <= n; j++) {// Find the pivot.
				if (Math.abs(a[j - 1][m - 2]) > Math.abs(x)) // if
																// (fabs(a[j][m-1])
																// > fabs(x))
				{
					x = a[j - 1][m - 2];// x=a[j][m-1];
					i = j;
				}
			}
			if (i != m) {// Interchange rows and columns.
				for (j = m - 1; j <= n; j++) {
					// SWAP(a[i][j],a[m][j])
					y = a[i - 1][j - 1];
					a[i - 1][j - 1] = a[m - 1][j - 1];
					a[m - 1][j - 1] = y;
				}
				for (j = 1; j <= n; j++) {
					// SWAP(a[j][i],a[j][m])
					y = a[j - 1][i - 1];
					a[j - 1][i - 1] = a[j - 1][m - 1];
					a[j - 1][m - 1] = y;
				}
			}
			if (x != 0.0)// if (x)
			{// Carry out the elimination.
				for (i = m + 1; i <= n; i++) {
					if ((y = a[i - 1][m - 2]) != 0.0) // if ((y=a[i][m-1]) !=
														// 0.0)
					{
						y /= x;
						a[i - 1][m - 2] = y;// a[i][m-1]=y;
						for (j = m; j <= n; j++)
							a[i - 1][j - 1] -= y * a[m - 1][j - 1];// a[i][j] -=
																	// y*a[m][j];
						for (j = 1; j <= n; j++)
							a[j - 1][m - 1] += y * a[j - 1][i - 1];// a[j][m] +=
																	// y*a[j][i];
					}
				}
			}
		}
	}

	/*
	 * The QR Algorithm for Real Hessenberg Matrices
	 */

	/**
	 * Finds all eigenvalues of an upper Hessenberg matrix a[1..n][1..n]. On 
	 * input a can be exactly as output from elmhes. On output it is destroyed. The real 
	 * and imaginary parts of the eigenvalues are returned in wr[1..n] and wi[1..n], respectively.
	 * @param a a
	 * @param n n
	 * @param wr wr
	 * @param wi wi
	 */
	public void hqr(double[][] a, int n, double wr[], double wi[])
	// Finds all eigenvalues of an upper Hessenberg matrix a[1..n][1..n]. On
	// input a can be
	// exactly as output from elmhes §11.5; on output it is destroyed. The real
	// and imaginary parts
	// of the eigenvalues are returned in wr[1..n] and wi[1..n], respectively.
	{
		int nn = 0;
		int m = 0;
		int l = 0;
		int k = 0;
		int j = 0;
		int its = 0;
		int i = 0;
		int mmin = 0;
		double z = 0.0;
		double y = 0.0;
		double x = 0.0;
		double w = 0.0;
		double v = 0.0;
		double u = 0.0;
		double t = 0.0;
		double s = 0.0;
		double r = 0.0;
		double q = 0.0;
		double p = 0.0;
		double anorm = 0.0;

		failB = false;

		anorm = 0.0;// Compute matrix norm for possible use in locating
					// single small subdiagonal element.
		for (i = 1; i <= n; i++)
			for (j = Math.max(i - 1, 1); j <= n; j++)
				anorm += Math.abs(a[i - 1][j - 1]);// anorm += fabs(a[i][j]);
		nn = n;
		t = 0.0; // Gets changed only by an exceptional shift.
		while (nn >= 1) {// Begin search for next eigenvalue.
			its = 0;
			do {
				for (l = nn; l >= 2; l--) {// Begin iteration: look for single
											// small subdiagonal element.
					s = Math.abs(a[l - 2][l - 2]) + Math.abs(a[l - 1][l - 1]);// s=fabs(a[l-1][l-1])+fabs(a[l][l]);
					if (s == 0.0)
						s = anorm;
					if ((double) (Math.abs(a[l - 1][l - 2]) + s) == s)// ((double)(Math.abs(a[l][l-1])
																		// + s)
																		// == s)
					{
						a[l - 1][l - 2] = 0.0;// a[l][l-1]=0.0;
						break;
					}
				}
				x = a[nn - 1][nn - 1];// x=a[nn][nn];
				if (l == nn) {// One root found.
					wr[nn - 1] = x + t;// wr[nn]=x+t;
					wi[nn-- - 1] = 0.0;// wi[nn--]=0.0;
				} else {
					y = a[nn - 2][nn - 2];// y=a[nn-1][nn-1];
					w = a[nn - 1][nn - 2] * a[nn - 2][nn - 1];// w=a[nn][nn-1]*a[nn-1][nn];
					if (l == (nn - 1)) {// Two roots found...
						p = 0.5 * (y - x);
						q = p * p + w;
						z = Math.sqrt(Math.abs(q));
						x += t;
						if (q >= 0.0) {// ...a real pair.
							z = p + SIGN(z, p);
							wr[nn - 2] = wr[nn - 1] = x + z;// wr[nn-1]=wr[nn]=x+z;
							if (z != 0.0)
								wr[nn - 1] = x - w / z;// wr[nn]=x-w/z;
							wi[nn - 2] = wi[nn - 1] = 0.0;// wi[nn-1]=wi[nn]=0.0;
						} else {// ...a complex pair.
							wr[nn - 2] = wr[nn - 1] = x + p;// wr[nn-1]=wr[nn]=x+p;
							wi[nn - 2] = -(wi[nn - 1] = z);// wi[nn-1]=
															// -(wi[nn]=z);
						}
						nn -= 2;
					} else {// No roots found. Continue iteration.
						if (its == 30) {
							// nrerror("Too many iterations in hqr");
							failB = true;
							failS = "Too many iterations in hqr";

							return;
						}
						if (its == 10 || its == 20) {// Form exceptional shift.
							t += x;
							for (i = 1; i <= nn; i++)
								a[i - 1][i - 1] -= x;// a[i][i] -= x;
							// s=fabs(a[nn][nn-1])+fabs(a[nn-1][nn-2]);
							s = Math.abs(a[nn - 1][nn - 2])
									+ Math.abs(a[nn - 2][nn - 3]);
							y = x = 0.75 * s;
							w = -0.4375 * s * s;
						}
						++its;
						for (m = (nn - 2); m >= l; m--) {// Form shift and then
															// look for 2
															// consecutive small
															// subdiagonal
															// elements.
							z = a[m - 1][m - 1];// z=a[m][m];
							r = x - z;
							s = y - z;
							// p=(r*s-w)/a[m+1][m]+a[m][m+1]; Equation
							// (11.6.23).
							p = (r * s - w) / a[m][m - 1] + a[m - 1][m];
							q = a[m][m] - z - r - s;// q=a[m+1][m+1]-z-r-s;
							r = a[m + 1][m];// r=a[m+2][m+1];
							s = Math.abs(p) + Math.abs(q) + Math.abs(r); // Scale
																			// to
																			// prevent
																			// overflow
																			// or
																			// underflow.
							p /= s;
							q /= s;
							r /= s;
							if (m == l)
								break;
							// u=fabs(a[m][m-1])*(fabs(q)+fabs(r));
							u = Math.abs(a[m - 1][m - 2])
									* (Math.abs(q) + Math.abs(r));
							// v=fabs(p)*(fabs(a[m-1][m-1])+fabs(z)+fabs(a[m+1][m+1]));
							v = Math.abs(p)
									* (Math.abs(a[m - 2][m - 2]) + Math.abs(z) + Math
											.abs(a[m][m]));
							if ((double) (u + v) == v)
								break; // Equation (11.6.26).
						}
						for (i = m + 2; i <= nn; i++) {
							a[i - 1][i - 3] = 0.0;// a[i][i-2]=0.0;
							if (i != (m + 2))
								a[i - 1][i - 4] = 0.0;// a[i][i-3]=0.0;
						}
						for (k = m; k <= nn - 1; k++) {
							// Double QR step on rows l to nn and columns m to
							// nn.
							if (k != m) {
								p = a[k - 1][k - 2];// p=a[k][k-1]; Begin setup
													// of Householder vector.
								q = a[k][k - 2];// q=a[k+1][k-1];
								r = 0.0;
								if (k != (nn - 1))
									r = a[k + 1][k - 2];// r=a[k+2][k-1];
								if ((x = Math.abs(p) + Math.abs(q)
										+ Math.abs(r)) != 0.0) {
									p /= x; // Scale to prevent overflow or
											// underflow.
									q /= x;
									r /= x;
								}
							}
							if ((s = SIGN(Math.sqrt(p * p + q * q + r * r), p)) != 0.0) {
								if (k == m) {
									if (l != m)
										a[k - 1][k - 2] = -a[k - 1][k - 2];// a[k][k-1]
																			// =
																			// -a[k][k-1];
								} else
									a[k - 1][k - 2] = -s * x;// a[k][k-1] =
																// -s*x;
								p += s;// Equations (11.6.24).
								x = p / s;
								y = q / s;
								z = r / s;
								q /= p;
								r /= p;
								for (j = k; j <= nn; j++) {// Row modification.
									p = a[k - 1][j - 1] + q * a[k][j - 1];// p=a[k][j]+q*a[k+1][j];
									if (k != (nn - 1)) {
										p += r * a[k + 1][j - 1];// p +=
																	// r*a[k+2][j];
										a[k + 1][j] -= p * z;// a[k+2][j] -=
																// p*z;
									}
									a[k][j - 1] -= p * y;// a[k+1][j] -= p*y;
									a[k - 1][j - 1] -= p * x;// a[k][j] -= p*x;
								}
								mmin = nn < k + 3 ? nn : k + 3;
								for (i = l; i <= mmin; i++) {// Column
																// modification.
									p = x * a[i - 1][k - 1] + y * a[i - 1][k];// p=x*a[i][k]+y*a[i][k+1];
									if (k != (nn - 1)) {
										p += z * a[i - 1][k + 1];// p +=
																	// z*a[i][k+2];
										a[i - 1][k + 1] -= p * r;// a[i][k+2] -=
																	// p*r;
									}
									a[i - 1][k] -= p * q;// a[i][k+1] -= p*q;
									a[i - 1][k - 1] -= p;// a[i][k] -= p;
								}
							}
						}
					}
				}
			} while (l < nn - 1);
		}
	}
}
