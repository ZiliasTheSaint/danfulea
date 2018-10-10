package danfulea.math.numerical;

/**
 * Sorting class<br>
 * Based on literature, including Numerical Recipes in C (Cambridge Univ.).
 * @author Dan Fulea, 10 OCT. 2006.
 */
public class Sorting {
	public static boolean failB = false;
	public static String failS = "";

	public static double a_swap = 0.0;
	public static double b_swap = 0.0;
	public static int ia_swap = 0;
	public static int ib_swap = 0;
	public static int M_quick = 7;
	public static int NSTACK_quick = 50;

	public static int M_selip = 64;
	public static double BIG = 1.0e30;

	/**
	 * Sorts an array arr[1..n] into ascending numerical order, by straight insertion. n is input; arr 
	 * is replaced on output by its sorted rearrangement.
	 * @param n n
	 * @param arr arr
	 */
	public static void piksrt(int n, double arr[])
	// Sorts an array arr[1..n] into ascending numerical order, by straight
	// insertion. n is input; arr
	// is replaced on output by its sorted rearrangement.
	{
		int i = 0;
		int j = 0;
		double a = 0.0;

		for (j = 2; j <= n; j++) {// Pick out each element in turn.
			a = arr[j - 1];// arr[j];
			i = j - 1;
			while (i > 0 && arr[i - 1] > a)// while (i > 0 && arr[i] > a)
			{// Look for the place to insert it.
				arr[i] = arr[i - 1];// arr[i+1]=arr[i];
				i--;
			}
			arr[i] = a;// arr[i+1]=a; //Insert it.
		}
	}

	/**
	 * Sorts an array arr[1..n] into ascending numerical order, by straight insertion, while making 
	 * the corresponding rearrangement of the array brr[1..n].
	 * @param n n
	 * @param arr arr
	 * @param brr brr
	 */
	public static void piksr2(int n, double arr[], double brr[])
	// Sorts an array arr[1..n] into ascending numerical order, by straight
	// insertion, while making
	// the corresponding rearrangement of the array brr[1..n].
	{
		int i = 0;
		int j = 0;
		double a = 0.0;
		double b = 0.0;

		for (j = 2; j <= n; j++) {// Pick out each element in turn.
			a = arr[j - 1];// a=arr[j];
			b = brr[j - 1];// b=brr[j];
			i = j - 1;
			while (i > 0 && arr[i - 1] > a)// while (i > 0 && arr[i] > a)
			{// Look for the place to insert it.
				arr[i] = arr[i - 1];// arr[i+1]=arr[i];
				brr[i] = brr[i - 1];// brr[i+1]=brr[i];
				i--;
			}
			arr[i] = a;// arr[i+1]=a; Insert it.
			brr[i] = b;// brr[i+1]=b;
		}
	}

	/**
	 * Sorts an array a[] into ascending numerical order by Shell’s method (diminishing increment 
	 * sort). a is replaced on output by its sorted rearrangement. Normally, the argument n should 
	 * be set to the size of array a, but if n is smaller than this, then only the first n elements of a 
	 * are sorted. This feature is used in selip.
	 * @param n n
	 * @param a a
	 */
	public static void shell(int n, double a[])// (unsigned long n, float a[])
	// Sorts an array a[] into ascending numerical order by Shell’s method
	// (diminishing increment
	// sort). a is replaced on output by its sorted rearrangement. Normally, the
	// argument n should
	// be set to the size of array a, but if n is smaller than this, then only
	// the first n elements of a
	// are sorted. This feature is used in selip.
	{
		// /unsigned long i,j,inc;
		int i = 0;
		int j = 0;
		int inc = 0;
		double v = 0.0;
		inc = 1;// Determine the starting increment.
		do {
			inc *= 3;
			inc++;
		} while (inc <= n);
		do {// Loop over the partial sorts.
			inc /= 3;
			for (i = inc + 1; i <= n; i++) {// Outer loop of straight insertion.
				v = a[i - 1];// a[i];
				j = i;
				while (a[j - inc - 1] > v)// while (a[j-inc] > v)
				{// Inner loop of straight insertion.
					a[j - 1] = a[j - inc - 1];// a[j]=a[j-inc];
					j -= inc;
					if (j <= inc)
						break;
				}
				a[j - 1] = v;// a[j]=v;
			}
		} while (inc > 1);
	}

	// ===============================================================
	/**
	 * SWAP two numbers of double type
	 * @param a a
	 * @param b b
	 */
	public static void SWAP(double a, double b) {
		a_swap = a;
		b_swap = b;
		double temp = a_swap;
		a_swap = b_swap;
		b_swap = temp;
	}

	/**
	 * SWAP two numbers of int type
	 * @param a a
	 * @param b b
	 */
	public static void SWAP(int a, int b) {
		ia_swap = a;
		ib_swap = b;
		int temp = ia_swap;
		ia_swap = ib_swap;
		ib_swap = temp;
	}

	// ===============================================================
	/**
	 * Sorts an array arr[1..n] into ascending numerical order using the Quicksort algorithm. n is 
	 * input; arr is replaced on output by its sorted rearrangement.
	 * @param n n
	 * @param arr arr
	 */
	public static void sort(int n, double arr[])// (unsigned long n, float
												// arr[])
	// Sorts an array arr[1..n] into ascending numerical order using the
	// Quicksort algorithm. n is
	// input; arr is replaced on output by its sorted rearrangement.
	{
		failB = false;
		// unsigned long i,ir=n,j,k,l=1,*istack;
		int i = 0;
		int ir = n;
		int j = 0;
		int k = 0;
		int l = 1;
		int[] istack = new int[NSTACK_quick];
		int jstack = 0;
		double a = 0.0;
		// double temp = 0.0;
		// istack=lvector(1,NSTACK);
		for (;;) {// Insertion sort when subarray small enough.
			if (ir - l < M_quick) {
				for (j = l + 1; j <= ir; j++) {
					a = arr[j - 1];// a=arr[j];
					for (i = j - 1; i >= l; i--) {
						if (arr[i - 1] <= a)
							break;// if (arr[i] <= a) break;
						arr[i] = arr[i - 1];// arr[i+1]=arr[i];
					}
					arr[i] = a;// arr[i+1]=a;
				}
				if (jstack == 0)
					break;
				ir = istack[jstack-- - 1];// ir=istack[jstack--]; //Pop stack
											// and begin a new round of
											// partitioning.
				l = istack[jstack-- - 1];// l=istack[jstack--];
			} else {
				k = (l + ir) >> 1;// Choose median of left, center, and right
									// elements
				// as partitioning element a. Also rearrange so that a[l] ?
				// a[l+1] ? a[ir].
				// SWAP(arr[k],arr[l+1]);
				SWAP(arr[k - 1], arr[l]);
				arr[k - 1] = a_swap;
				arr[l] = b_swap;
				if (arr[l - 1] > arr[ir - 1])// if (arr[l] > arr[ir])
				{
					// SWAP(arr[l],arr[ir])
					SWAP(arr[l - 1], arr[ir - 1]);
					arr[l - 1] = a_swap;
					arr[ir - 1] = b_swap;
				}
				if (arr[l] > arr[ir - 1])// //if (arr[l+1] > arr[ir])
				{
					// SWAP(arr[l+1],arr[ir])
					SWAP(arr[l], arr[ir - 1]);
					arr[l] = a_swap;
					arr[ir - 1] = b_swap;
				}
				if (arr[l - 1] > arr[l])// if (arr[l] > arr[l+1])
				{
					// SWAP(arr[l],arr[l+1])
					SWAP(arr[l - 1], arr[l]);
					arr[l - 1] = a_swap;
					arr[l] = b_swap;
				}
				i = l + 1;// Initialize pointers for partitioning.
				j = ir;
				a = arr[l];// a=arr[l+1]; Partitioning element.
				for (;;) {// Beginning of innermost loop.
					do
						i++;
					while (arr[i - 1] < a);// while (arr[i] < a); Scan up to
											// find element > a.
					do
						j--;
					while (arr[j - 1] > a);// while (arr[j] > a); Scan down to
											// find element < a.
					if (j < i)
						break;// Pointers crossed. Partitioning complete.
					// SWAP(arr[i],arr[j]); Exchange elements.
					SWAP(arr[i - 1], arr[j - 1]);
					arr[i - 1] = a_swap;
					arr[j - 1] = b_swap;
				}// End of innermost loop.
				arr[l] = arr[j - 1];// arr[l+1]=arr[j]; Insert partitioning
									// element.
				arr[j - 1] = a;// arr[j]=a;
				jstack += 2;
				// Push pointers to larger subarray on stack, process smaller
				// subarray immediately.
				if (jstack > NSTACK_quick) {
					// nrerror("NSTACK too small in sort.");
					failB = true;
					failS = "NSTACK too small in sort.";
					return;
				}
				if (ir - i + 1 >= j - l) {
					istack[jstack - 1] = ir;// istack[jstack]=ir;
					istack[jstack - 2] = i;// istack[jstack-1]=i;
					ir = j - 1;
				} else {
					istack[jstack - 1] = j - 1;// istack[jstack]=j-1;
					istack[jstack - 2] = l;// istack[jstack-1]=l;
					l = i;
				}
			}
		}
		// free_lvector(istack,1,NSTACK);
	}

	// coordinate systems!!!
	/**
	 * Sorts an array arr[1..n] into ascending order using Quicksort, while making the corresponding 
	 * rearrangement of the array brr[1..n]. Useful for sampling a coordinate system (we sort x-value 
	 * and want y to be changed accordingly as they are pairs).
	 * @param n n
	 * @param arr arr
	 * @param brr brr
	 */
	public static void sort2(int n, double[] arr, double[] brr)
	// Sorts an array arr[1..n] into ascending order using Quicksort, while
	// making the corresponding
	// rearrangement of the array brr[1..n].
	{
		failB = false;
		int i = 0;
		int ir = n;
		int j = 0;
		int k = 0;
		int l = 1;
		int[] istack = new int[NSTACK_quick];
		int jstack = 0;
		double a = 0.0;
		// double temp = 0.0;
		double b = 0.0;

		for (;;) {// Insertion sort when subarray small enough.
			if (ir - l < M_quick) {
				for (j = l + 1; j <= ir; j++) {
					a = arr[j - 1];// a=arr[j];
					b = brr[j - 1];// b=brr[j];
					for (i = j - 1; i >= l; i--) {
						if (arr[i - 1] <= a)
							break;// if (arr[i] <= a) break;
						arr[i] = arr[i - 1];// arr[i+1]=arr[i];
						brr[i] = brr[i - 1];// brr[i+1]=brr[i];
					}
					arr[i] = a;// arr[i+1]=a;
					brr[i] = b;// brr[i+1]=b;
				}
				if (jstack == 0)// if (!jstack)//if (jstack == 0) break;
				{
					// free_lvector(istack,1,NSTACK);
					return;
				}
				ir = istack[jstack - 1];// ir=istack[jstack]; //Pop stack and
										// begin a new round of partitioning.
				l = istack[jstack - 2];// l=istack[jstack-1];
				jstack -= 2;
			} else {
				k = (l + ir) >> 1;// Choose median of left, center and right
									// elements
				// as partitioning element a. Also rearrange so that a[l] ?
				// a[l+1] ? a[ir].
				// SWAP(arr[k],arr[l+1])
				SWAP(arr[k - 1], arr[l]);
				arr[k - 1] = a_swap;
				arr[l] = b_swap;
				// SWAP(brr[k],brr[l+1])
				SWAP(brr[k - 1], brr[l]);
				brr[k - 1] = a_swap;
				brr[l] = b_swap;
				if (arr[l - 1] > arr[ir - 1])// if (arr[l] > arr[ir])
				{
					// SWAP(arr[l],arr[ir])
					SWAP(arr[l - 1], arr[ir - 1]);
					arr[l - 1] = a_swap;
					arr[ir - 1] = b_swap;
					// SWAP(brr[l],brr[ir])
					SWAP(brr[l - 1], brr[ir - 1]);
					brr[l - 1] = a_swap;
					brr[ir - 1] = b_swap;
				}
				if (arr[l] > arr[ir - 1]) // if (arr[l+1] > arr[ir])
				{
					// SWAP(arr[l+1],arr[ir])
					SWAP(arr[l], arr[ir - 1]);
					arr[l] = a_swap;
					arr[ir - 1] = b_swap;
					// SWAP(brr[l+1],brr[ir])
					SWAP(brr[l], brr[ir - 1]);
					brr[l] = a_swap;
					brr[ir - 1] = b_swap;
				}
				if (arr[l - 1] > arr[l])// if (arr[l] > arr[l+1])
				{
					// SWAP(arr[l],arr[l+1])
					SWAP(arr[l - 1], arr[l]);
					arr[l - 1] = a_swap;
					arr[l] = b_swap;
					// SWAP(brr[l],brr[l+1])
					SWAP(brr[l - 1], brr[l]);
					brr[l - 1] = a_swap;
					brr[l] = b_swap;
				}
				i = l + 1;// Initialize pointers for partitioning.
				j = ir;
				a = arr[l];// a=arr[l+1]; Partitioning element.
				b = brr[l];// b=brr[l+1];
				for (;;) {// Beginning of innermost loop.
					do
						i++;
					while (arr[i - 1] < a);// while (arr[i] < a); Scan up to
											// find element > a.
					do
						j--;
					while (arr[j - 1] > a);// while (arr[j] > a); Scan down to
											// find element < a.
					if (j < i)
						break;// Pointers crossed. Partitioning complete.
					// SWAP(arr[i],arr[j]) Exchange elements of both arrays.
					SWAP(arr[i - 1], arr[j - 1]);
					arr[i - 1] = a_swap;
					arr[j - 1] = b_swap;
					// SWAP(brr[i],brr[j])
					SWAP(brr[i - 1], brr[j - 1]);
					brr[i - 1] = a_swap;
					brr[j - 1] = b_swap;
				}// End of innermost loop.
				arr[l] = arr[j - 1];// arr[l+1]=arr[j]; Insert partitioning
									// element in both arrays.
				arr[j - 1] = a;// arr[j]=a;
				brr[l] = brr[j - 1];// brr[l+1]=brr[j];
				brr[j - 1] = b;// brr[j]=b;
				jstack += 2;
				// Push pointers to larger subarray on stack, process smaller
				// subarray immediately.
				if (jstack > NSTACK_quick) {
					// nrerror("NSTACK too small in sort2.");
					failB = true;
					failS = "NSTACK too small in sort2.";
					return;
				}
				if (ir - i + 1 >= j - l) {
					istack[jstack - 1] = ir;// istack[jstack]=ir;
					istack[jstack - 2] = i;// istack[jstack-1]=i;
					ir = j - 1;
				} else {
					istack[jstack - 1] = j - 1;// istack[jstack]=j-1;
					istack[jstack - 2] = l;// istack[jstack-1]=l;
					l = i;
				}
			}// else
		}// for
	}

	/**
	 * Sorts an array ra[1..n] into ascending numerical order using the Heapsort algorithm. n is 
	 * input; ra is replaced on output by its sorted rearrangement.
	 * @param n n
	 * @param ra ra
	 */
	public static void hpsort(int n, double[] ra)// unsigned long n, float ra[])
	// Sorts an array ra[1..n] into ascending numerical order using the Heapsort
	// algorithm. n is
	// input; ra is replaced on output by its sorted rearrangement.
	{
		// unsigned long i,ir,j,l;
		int i = 0;
		int ir = 0;
		int j = 0;
		int l = 0;
		double rra = 0.0;
		if (n < 2)
			return;

		l = (n >> 1) + 1;
		ir = n;
		// The index l will be decremented from its initial value down to 1
		// during the “hiring” (heap
		// creation) phase. Once it reaches 1, the index ir will be decremented
		// from its initial value
		// down to 1 during the “retirement-and-promotion” (heap selection)
		// phase.
		for (;;) {
			if (l > 1) {// Still in hiring phase.
				rra = ra[--l - 1];// rra=ra[--l];
			} else {// In retirement-and-promotion phase.
				rra = ra[ir - 1];// rra=ra[ir];// Clear a space at end of array.
				ra[ir - 1] = ra[0];// ra[ir]=ra[1]; Retire the top of the heap
									// into it.
				if (--ir == 1) {// Done with the last promotion.
					ra[0] = rra;// ra[1]=rra; The least competent worker of all!
					break;
				}
			}
			i = l; // Whether in the hiring phase or promotion phase, we
					// here set up to sift down element rra to its proper level.
			j = l + l;
			while (j <= ir) {
				// if (j < ir && ra[j] < ra[j+1]) j++; Compare to the better
				// underling.
				if (j < ir && ra[j - 1] < ra[j])
					j++;
				if (rra < ra[j - 1]) // if (rra < ra[j])
				{// Demote rra.
					ra[i - 1] = ra[j - 1];// ra[i]=ra[j];
					i = j;
					j <<= 1;
				} else
					break;// Found rra’s level. Terminate the sift-down.
			}
			ra[i - 1] = rra;// ra[i]=rra; Put rra into its slot.
		}
	}

	// void indexx(unsigned long n, float arr[], unsigned long indx[])
	/**
	 * Indexes an array arr[1..n], i.e., outputs the array indx[1..n] such that arr[indx[j]] is 
	 * in ascending order for j = 1, 2, . . . ,N. The input quantities n and arr are not changed.
	 * @param n n
	 * @param arr arr
	 * @param indx indx
	 */
	public static void indexx(int n, double[] arr, int[] indx)
	// Indexes an array arr[1..n], i.e., outputs the array indx[1..n] such that
	// arr[indx[j]] is
	// in ascending order for j = 1, 2, . . . ,N. The input quantities n and arr
	// are not changed.
	{
		// unsigned long i,indxt,ir=n,itemp,j,k,l=1;
		int i = 0;
		int indxt = 0;
		int ir = n;
		// int itemp = 0;
		int j = 0;
		int k = 0;
		int l = 1;
		int jstack = 0;
		int[] istack = new int[NSTACK_quick];
		double a = 0.0;
		// istack=ivector(1,NSTACK);
		failB = false;

		for (j = 1; j <= n; j++)
			indx[j - 1] = j;// indx[j]=j;
		for (;;) {
			if (ir - l < M_quick) {
				for (j = l + 1; j <= ir; j++) {
					indxt = indx[j - 1];// indxt=indx[j];
					a = arr[indxt - 1];// a=arr[indxt];
					for (i = j - 1; i >= l; i--) {
						// if (arr[indx[i]] <= a) break;
						if (arr[indx[i - 1] - 1] <= a)
							break;
						indx[i] = indx[i - 1];// indx[i+1]=indx[i];
					}
					indx[i] = indxt;// indx[i+1]=indxt;
				}
				if (jstack == 0)
					break;
				ir = istack[jstack-- - 1];// ir=istack[jstack--];
				l = istack[jstack-- - 1];// l=istack[jstack--];
			} else {
				k = (l + ir) >> 1;
				// SWAP(indx[k],indx[l+1]);
				SWAP(indx[k - 1], indx[l]);
				indx[k - 1] = ia_swap;
				indx[l] = ib_swap;
				if (arr[indx[l - 1] - 1] > arr[indx[ir - 1] - 1]) // if
																	// (arr[indx[l]]
																	// >
																	// arr[indx[ir]])
				{
					// SWAP(indx[l],indx[ir])
					SWAP(indx[l - 1], indx[ir - 1]);
					indx[l - 1] = ia_swap;
					indx[ir - 1] = ib_swap;
				}
				if (arr[indx[l] - 1] > arr[indx[ir - 1] - 1])// if
																// (arr[indx[l+1]]
																// >
																// arr[indx[ir]])
				{
					// SWAP(indx[l+1],indx[ir])
					SWAP(indx[l], indx[ir - 1]);
					indx[l] = ia_swap;
					indx[ir - 1] = ib_swap;
				}
				if (arr[indx[l - 1] - 1] > arr[indx[l] - 1])// if (arr[indx[l]]
															// > arr[indx[l+1]])
				{
					// SWAP(indx[l],indx[l+1])
					SWAP(indx[l - 1], indx[l]);
					indx[l - 1] = ia_swap;
					indx[l] = ib_swap;
				}
				i = l + 1;
				j = ir;
				indxt = indx[l];// indxt=indx[l+1];
				a = arr[indxt - 1];// a=arr[indxt];
				for (;;) {
					do
						i++;
					while (arr[indx[i - 1] - 1] < a);// while (arr[indx[i]] <
														// a);
					do
						j--;
					while (arr[indx[j - 1] - 1] > a);// while (arr[indx[j]] >
														// a);
					if (j < i)
						break;
					// SWAP(indx[i],indx[j])
					SWAP(indx[i - 1], indx[j - 1]);
					indx[i - 1] = ia_swap;
					indx[j - 1] = ib_swap;
				}
				indx[l] = indx[j - 1];// indx[l+1]=indx[j];
				indx[j - 1] = indxt;// indx[j]=indxt;
				jstack += 2;
				if (jstack > NSTACK_quick) {
					// nrerror("NSTACK too small in indexx.");
					failB = true;
					failS = "NSTACK too small in indexx.";
					return;
				}
				if (ir - i + 1 >= j - l) {
					istack[jstack - 1] = ir;// istack[jstack]=ir;
					istack[jstack - 2] = i;// istack[jstack-1]=i;
					ir = j - 1;
				} else {
					istack[jstack - 1] = j - 1;// istack[jstack]=j-1;
					istack[jstack - 2] = l;// istack[jstack-1]=l;
					l = i;
				}
			}
		}
		// free_ivector(istack,1,NSTACK);
	}

	/*
	 * A rank table is different from an index table. A rank table’s jth entry
	 * gives the rank of the jth element of the original array of keys, ranging
	 * from 1 (if that element was the smallest) to N (if that element was the
	 * largest). One can easily construct a rank table from an index table,
	 * however:
	 */
	// public static void rank(unsigned long n, unsigned long indx[], unsigned
	// long irank[])
	/**
	 * Given indx[1..n] as output from the routine indexx, returns an array irank[1..n], the 
	 * corresponding table of ranks.
	 * @param n n
	 * @param indx indx
	 * @param irank irank
	 */
	public static void rank(int n, int[] indx, int[] irank)
	// Given indx[1..n] as output from the routine indexx, returns an array
	// irank[1..n], the
	// corresponding table of ranks.
	{
		int j = 0;// unsigned long j;
		for (j = 1; j <= n; j++)
			irank[indx[j - 1] - 1] = j;// irank[indx[j]]=j;
	}

	// @@Selecting the Mth Largest
	// public static double select(unsigned long k, unsigned long n, float
	// arr[])
	// k=1,2,...,N
	/**
	 * Returns the kth smallest value in the array arr[1..n]. The input array will be rearranged 
	 * to have this value in location arr[k], with all smaller elements moved to arr[1..k-1] (in 
	 * arbitrary order) and all larger elements in arr[k+1..n] (also in arbitrary order).
	 * @param k k
	 * @param n n
	 * @param arr arr
	 * @return the result
	 */
	public static double select(int k, int n, double[] arr)
	// Returns the kth smallest value in the array arr[1..n]. The input array
	// will be rearranged
	// to have this value in location arr[k], with all smaller elements moved to
	// arr[1..k-1] (in
	// arbitrary order) and all larger elements in arr[k+1..n] (also in
	// arbitrary order).
	{
		// unsigned long i,ir,j,l,mid;
		int i = 0;
		int ir = 0;
		int j = 0;
		int l = 0;
		int mid = 0;
		double a = 0.0;
		// double temp = 0.0;

		l = 1;
		ir = n;
		for (;;) {
			if (ir <= l + 1) {// Active partition contains 1 or 2 elements.
				if (ir == l + 1 && arr[ir - 1] < arr[l - 1])// if (ir == l+1 &&
															// arr[ir] < arr[l])
				{// Case of 2el ements.
					// SWAP(arr[l],arr[ir])
					SWAP(arr[l - 1], arr[ir - 1]);
					arr[l - 1] = a_swap;
					arr[ir - 1] = b_swap;
				}
				return arr[k - 1];// return arr[k];
			} else {
				mid = (l + ir) >> 1; // Choose median of left, center, and right
										// elements
				// as partitioning element a. Also rearrange so that arr[l] <=
				// arr[l+1], arr[ir] >= arr[l+1].
				// SWAP(arr[mid],arr[l+1])
				SWAP(arr[mid - 1], arr[l]);
				arr[mid - 1] = a_swap;
				arr[l] = b_swap;
				if (arr[l - 1] > arr[ir - 1]) // if (arr[l] > arr[ir])
				{
					// SWAP(arr[l],arr[ir])
					SWAP(arr[l - 1], arr[ir - 1]);
					arr[l - 1] = a_swap;
					arr[ir - 1] = b_swap;
				}
				if (arr[l] > arr[ir - 1]) // if (arr[l+1] > arr[ir])
				{
					// SWAP(arr[l+1],arr[ir])
					SWAP(arr[l], arr[ir - 1]);
					arr[l] = a_swap;
					arr[ir - 1] = b_swap;
				}
				if (arr[l - 1] > arr[l]) // if (arr[l] > arr[l+1])
				{
					// SWAP(arr[l],arr[l+1])
					SWAP(arr[l - 1], arr[l]);
					arr[l - 1] = a_swap;
					arr[l] = b_swap;
				}
				i = l + 1;// Initialize pointers for partitioning.
				j = ir;
				a = arr[l];// a=arr[l+1]; Partitioning element.
				for (;;) {// Beginning of innermost loop.
					do
						i++;
					while (arr[i - 1] < a);// while (arr[i] < a); Scan up to
											// find element > a.
					do
						j--;
					while (arr[j - 1] > a);// while (arr[j] > a); Scan down to
											// find element < a.
					if (j < i)
						break; // Pointers crossed. Partitioning complete.
					// SWAP(arr[i],arr[j])
					SWAP(arr[i - 1], arr[j - 1]);
					arr[i - 1] = a_swap;
					arr[j - 1] = b_swap;
				}// End of innermost loop.
				arr[l] = arr[j - 1];// arr[l+1]=arr[j]; Insert partitioning
									// element.
				arr[j - 1] = a;// arr[j]=a;
				if (j >= k)
					ir = j - 1;// Keep active the partition that contains the
								// kth element.
				if (j <= k)
					l = i;
			}
		}
	}

	// public static double selip(unsigned long k, unsigned long n, float arr[])
	/**
	 * Returns the kth smallest value in the array arr[1..n]. The input array is not altered.
	 * @param k k
	 * @param n n
	 * @param arr arr
	 * @return the result
	 */
	public static double selip(int k, int n, double[] arr)
	// Returns the kth smallest value in the array arr[1..n]. The input array is
	// not altered.
	{
		// void shell(unsigned long n, float a[]);
		// unsigned long i,j,jl,jm,ju,kk,mm,nlo,nxtmm,*isel;
		int i = 0;
		int j = 0;
		int jl = 0;
		int jm = 0;
		int ju = 0;
		int kk = 0;
		int mm = 0;
		int nlo = 0;
		int nxtmm = 0;
		int[] isel = new int[M_selip + 2];
		double ahi = 0.0;
		double alo = 0.0;
		double sum = 0.0;
		double[] sel = new double[M_selip + 2];

		failB = false;

		if (k < 1 || k > n || n <= 0) {
			// nrerror("bad input to selip");
			failB = true;
			failS = "bad input to selip";
			return 1;
		}
		// isel=lvector(1,M+2);
		// sel=vector(1,M+2);
		kk = k;
		ahi = BIG;
		alo = -BIG;
		for (;;) {// Main iteration loop, until desired element is isolated.
			mm = nlo = 0;
			sum = 0.0;
			nxtmm = M_selip + 1;
			for (i = 1; i <= n; i++) {// Make a pass through the whole array.
				if (arr[i - 1] >= alo && arr[i - 1] <= ahi) // if (arr[i] >= alo
															// && arr[i] <= ahi)
				{
					// Consider only elements in the current brackets.
					mm++;
					if (arr[i - 1] == alo)
						nlo++;// if (arr[i] == alo) nlo++;// In case of ties for
								// low bracket.
					// Now use statistical procedure for selecting m in-range
					// elements with equal
					// probability, even without knowing in advance how many
					// there are!
					if (mm <= M_selip)
						sel[mm - 1] = arr[i - 1];// sel[mm]=arr[i];
					else if (mm == nxtmm) {
						nxtmm = mm + mm / M_selip;
						sel[((i + mm + kk) % M_selip)] = arr[i - 1];// sel[1 +
																	// ((i+mm+kk)
																	// %
																	// M)]=arr[i];
						// The % operation provides a somewhat random number.
					}
					sum += arr[i - 1];// sum += arr[i];
				}
			}
			if (kk <= nlo) {// Desired element is tied for lower bound;return
							// it.
							// FREEALL
				return alo;
			} else if (mm <= M_selip) {// All in-range elements were kept. So
										// return answer by direct method.
				shell(mm, sel);
				ahi = sel[kk - 1];// ahi = sel[kk];
				// FREEALL
				return ahi;
			}
			sel[M_selip] = sum / mm;// sel[M+1]=sum/mm; Augment selected set by
									// mean value (fixes degeneracies), and sort
									// it.
			shell(M_selip + 1, sel);
			sel[M_selip + 1] = ahi;// sel[M_selip+2]=ahi;
			for (j = 1; j <= M_selip + 2; j++)
				isel[j - 1] = 0;// isel[j]=0; Zero the count array.
			for (i = 1; i <= n; i++) {// Make another pass through the array.
				if (arr[i - 1] >= alo && arr[i - 1] <= ahi)// if (arr[i] >= alo
															// && arr[i] <= ahi)
				{// For each in-range element..
					jl = 0;
					ju = M_selip + 2;
					while (ju - jl > 1) {// ...find its position among the
											// select by bisection...
						jm = (ju + jl) / 2;
						if (arr[i - 1] >= sel[jm - 1])
							jl = jm;// if (arr[i] >= sel[jm]) jl=jm;
						else
							ju = jm;
					}
					isel[ju - 1]++;// isel[ju]++; //...and increment the
									// counter.
				}
			}
			j = 1;// Now we can narrow the bounds to just one bin, that is, by a
					// factor of order m.
			while (kk > isel[j - 1]) // while (kk > isel[j])
			{
				alo = sel[j - 1];// alo=sel[j];
				kk -= isel[j++ - 1];// kk -= isel[j++];
			}
			ahi = sel[j - 1];// ahi=sel[j];
		}
	}

	/*
	 * Of course neither of the above routines should be used for the trivial
	 * cases of finding the largest, or smallest, element in an array. Those
	 * cases, you code by hand as simple for loops. There are also good ways to
	 * code the case where k is modest in comparison to N, so that extra memory
	 * of order k is not burdensome. An example is to use the method of Heapsort
	 * (§8.3) to make a single pass through an array of length N while saving
	 * the m largest elements.
	 */
	// public static void hpsel(unsigned long m, unsigned long n, float arr[],
	// float heap[])
	/**
	 * Returns in heap[1..m] the largest m elements of the array arr[1..n], with 
	 * heap[1] guaranteed to be the the mth largest element. The array arr is not altered. For 
	 * efficiency, this routine should be used only when m less or equal than n/2.
	 * @param m m
	 * @param n n
	 * @param arr arr
	 * @param heap heap
	 */
	public static void hpsel(int m, int n, double[] arr, double[] heap)
	// Returns in heap[1..m] the largest m elements of the array arr[1..n], with
	// heap[1] guaranteed
	// to be the the mth largest element. The array arr is not altered. For
	// efficiency, this routine
	// should be used only when m <= n/2.
	{
		// void sort(unsigned long n, float arr[]);
		// void nrerror(char error_text[]);
		failB = false;
		// unsigned long i,j,k;
		int i = 0;
		int j = 0;
		int k = 0;
		double swap = 0.0;

		// if (m > n/2 || m < 1)
		// {
		// nrerror("probable misuse of hpsel");
		// }

		for (i = 1; i <= m; i++)
			heap[i - 1] = arr[i - 1];// heap[i]=arr[i];
		sort(m, heap); // Create initial heap by overkill! We assume m <= n.
		for (i = m + 1; i <= n; i++) {// For each remaining element...
			if (arr[i - 1] > heap[0])// if (arr[i] > heap[1])
			{// Put it on the heap?
				heap[0] = arr[i - 1];// heap[1]=arr[i];
				for (j = 1;;) {// Sift down.
					k = j << 1;
					if (k > m)
						break;
					// if (k != m && heap[k] > heap[k+1]) k++;
					if (k != m && heap[k - 1] > heap[k])
						k++;
					if (heap[j - 1] <= heap[k - 1])
						break;// if (heap[j] <= heap[k]) break;
					swap = heap[k - 1];// swap=heap[k];
					heap[k - 1] = heap[j - 1];// heap[k]=heap[j];
					heap[j - 1] = swap;// heap[j]=swap;
					j = k;
				}
			}
		}
	}

	// @@@@@@@@@Determination of Equivalence Classes
	/*
	 * void eclass(int nf[], int n, int lista[], int listb[], int m) Given m
	 * equivalences between pairs of n individual elements in the form of the
	 * input arrays lista[1..m] and listb[1..m], this routine returns in
	 * nf[1..n] the number of the equivalence class of each of the n elements,
	 * integers between 1 and n (not all such integers used). { int l,k,j; for
	 * (k=1;k<=n;k++) nf[k]=k; Initialize each element its own class. for
	 * (l=1;l<=m;l++) { For each piece of input information... j=lista[l]; while
	 * (nf[j] != j) j=nf[j]; Track first element up to its ancestor. k=listb[l];
	 * while (nf[k] != k) k=nf[k]; Track second element up to its ancestor. if
	 * (j != k) nf[j]=k; If they are not already related, make them so. } for
	 * (j=1;j<=n;j++) Final sweep up to highest ancestors. while (nf[j] !=
	 * nf[nf[j]]) nf[j]=nf[nf[j]]; }
	 * 
	 * 
	 * void eclazz(int nf[], int n, int (*equiv)(int, int)) Given a
	 * user-supplied boolean function equiv which tells whether a pair of
	 * elements, each in the range 1...n, are related, return in nf[1..n]
	 * equivalence class numbers for each element. { int kk,jj; nf[1]=1; for
	 * (jj=2;jj<=n;jj++) { Loop over first element of all pairs. nf[jj]=jj; for
	 * (kk=1;kk<=(jj-1);kk++) { Loop over second element of all pairs.
	 * nf[kk]=nf[nf[kk]]; Sweep it up this much. if ((*equiv)(jj,kk))
	 * nf[nf[nf[kk]]]=jj; Good exercise for the reader to figure out why this
	 * much ancestry is necessary! } } for (jj=1;jj<=n;jj++) nf[jj]=nf[nf[jj]];
	 * Only this much sweeping is needed finally. }
	 */
}
