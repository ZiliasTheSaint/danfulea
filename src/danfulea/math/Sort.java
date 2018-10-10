package danfulea.math;

import java.util.Vector;

/**
 * Sorting utilities.
 * 
 * @author Dan Fulea, 14 MAY 2011
 * 
 */
public class Sort {
	private static int nearestposition = 0;

	/**
	 * MATRIX sorting
	 * 
	 * @param a
	 *            the matrix vector
	 * @param ncols
	 *            matrix number of "columns"
	 * @param sortcolindex
	 *            "column" index for sorting
	 * @param lo0
	 *            low row index
	 * @param hi0
	 *            high row index
	 */
	private static void qSort(Vector<Object> a, int ncols, int sortcolindex,
			int lo0, int hi0) {

		int tmp = 0;
		int lo = lo0;
		int hi = hi0;

		if (lo >= hi) {
			return;
		}

		int midInt = (lo + hi) / 2;
		String midD = getValueAt(a, midInt, sortcolindex).toString();
		double mid = Convertor.stringToDouble(midD);

		while (lo < hi) {
			String aloS = getValueAt(a, lo, sortcolindex).toString();
			double alo = Convertor.stringToDouble(aloS);

			String ahiS = getValueAt(a, hi, sortcolindex).toString();
			double ahi = Convertor.stringToDouble(ahiS);

			while (lo < hi && alo < mid) {
				lo++;

				aloS = getValueAt(a, lo, sortcolindex).toString();
				alo = Convertor.stringToDouble(aloS);
			}

			while (lo < hi && ahi > mid) {
				hi--;

				ahiS = getValueAt(a, hi, sortcolindex).toString();
				ahi = Convertor.stringToDouble(ahiS);
			}
			if (lo < hi) {

				for (int j = 0; j < ncols; j++) {
					String tmpS = getValueAt(a, lo, j).toString();
					String setS = getValueAt(a, hi, j).toString();
					setValueAt(setS, a, lo, j);
					setValueAt(tmpS, a, hi, j);
				}

				lo++;
				hi--;
			}
		}
		if (hi < lo) {
			tmp = hi;
			hi = lo;
			lo = tmp;
		}
		// recursive sort
		qSort(a, ncols, sortcolindex, lo0, lo);
		qSort(a, ncols, sortcolindex, lo == lo0 ? lo + 1 : lo, hi0);
	}

	/**
	 * Matrix get value method
	 * 
	 * @param a
	 *            the matrix vector
	 * @param row
	 *            row index
	 * @param col
	 *            column index
	 * @return the value Object at specified row and col indexes.
	 */
	private static Object getValueAt(Vector<Object> a, int row, int col) {
		@SuppressWarnings("unchecked")
		Vector<Object> v = (Vector<Object>) a.elementAt(row);
		return v.elementAt(col);
	}

	/**
	 * MATRIX set value method
	 * 
	 * @param value
	 *            the value to be set
	 * @param a
	 *            the matrix vector
	 * @param row
	 *            row index
	 * @param col
	 *            col index
	 */
	private static void setValueAt(Object value, Vector<Object> a, int row,
			int col) {
		@SuppressWarnings("unchecked")
		Vector<Object> v = (Vector<Object>) a.elementAt(row);
		v.setElementAt(value, col);
	}

	/**
	 * public quick sort method for matrix vectors
	 * 
	 * @param a
	 *            the MATRIX vector
	 * @param ncols
	 *            number of columns
	 * @param sortcolindex
	 *            column index for sorting!
	 */
	public static void qSort(Vector<Object> a, int ncols, int sortcolindex) {
		// NOTES: Object parameters such as VECTORS, ARRAYS will be modified by
		// this
		// methods. Therefore the MATRIX will be affected...will be sorted!!!!
		qSort(a, ncols, sortcolindex, 0, a.size() - 1);
	}

	// gaseste cea mai apropiata valoare dintr-un sir fata de una
	// data.!!notsorted array!!
	// flag=true->valoarea din sir sa fie mai mica--LOWER THAN VALUE
	// flag=false->valoarea din sir sa fie mai mare--NOT LOWER THAN VALUE
	/**
	 * Find the nearest value from an array relative to the given value.
	 *  If flag is true then the array value should be lower than the given value. 
	 *  If flag is set to false then the array value should be greater than the given value. 
	 *  The input array is supposed to be already sorted.
	 * @param a the input array
	 * @param value the given value
	 * @param flag the flag
	 * @return the value from array
	 */
	public static double findNearestValue(double[] a, double value, boolean flag) {
		boolean b = true;
		int ip = 0;

		if (a.length > 1) {
			// double[] a1 = newQSort(a);
			while (b) {
				if (flag) {
					if ((a[ip] <= value) && (a[ip + 1] > value)) {
						break;
					}
				} else {
					if (ip > 0)
						if ((a[ip] >= value) && (a[ip - 1] < value)) {
							break;
						}
				}

				ip++;
				if (ip == a.length - 1) {
					b = false;
					break;
				}
			}
			nearestposition = ip;// ----------------
			return a[ip];
		} else {
			nearestposition = 0;// ---------------
			return a[0];
		}
	}

	 
	/**
	 * Same as findNearestValue but the input array is sorted in ascending order.
	  * @param a the input array
	 * @param value the given value
	 * @param flag the flag
	 * @return the value from array
	 */
	public static double findNearestValueSorted(double[] a, double value, boolean flag)
	{
		boolean b=true;
		int ip=0;

        if (a.length > 1)
        {
           double[] a1 = newQSort(a);
           while(b)
           {
			   if(flag)
			   {
                  if((a1[ip]<=value) && (a1[ip+1]>value))
                  {
					  break;
				  }
			   }
			   else
			   {
				  if (ip>0)
                  	if((a1[ip]>=value) && (a1[ip-1]<value))
                  	{
					  break;
				  	}
			   }

			   ip++;
			   if(ip==a1.length-1)
			   {
				  b=false;
				  break;
			   }
		   }
		   nearestposition=ip;//----------------
		   return a1[ip];
		}
		else
		{
		   nearestposition=0;//----------------
		   return a[0];
	    }
	}
	
	/**
	 * 
	 * @return the nearest position index corresponding to the nearest value.
	 */
	public static int getNearestPosition() {
		return nearestposition;
	}

	// ================The following are internal============================
	/**
	 * Quick sort algorithm
	 * @param a the array of doubles
	 * @param lo0 input lower index
	 * @param hi0 input upper index
	 */
	private static void qSort(double[] a, int lo0, int hi0)
	{
		double temp = 0.0;
		int tmp = 0;
  		int lo = lo0;
		int hi = hi0;

		if (lo >= hi)
		{
		    return;//parasirea metodei
	    }

		double mid = a[(lo + hi) / 2];
		while (lo < hi)
		{
	    	while (lo<hi && a[lo] < mid)
	    	{
				lo++;
	    	}
	    	while (lo<hi && a[hi] > mid)
	    	{
				hi--;
	    	}
	    	if (lo < hi)
	    	{
				temp = a[lo];
				a[lo] = a[hi];
				a[hi] = temp;
                lo++;
                hi--;
		    }
	    }
	    if (hi < lo)
	    {
	    	tmp = hi;
	    	hi = lo;
	    	lo = tmp;
		}
		qSort(a, lo0, lo);
		qSort(a, lo == lo0 ? lo+1 : lo, hi0);

    }
	
	/**
	 * Quick sort algorithm
	 * @param a the array of ints
	 * @param lo0 input lower index
	 * @param hi0 input upper index
	 */
	private static void qSortInt(int[] a, int lo0, int hi0)
	{
		int temp = 0;
		int tmp = 0;
  		int lo = lo0;
		int hi = hi0;

		if (lo >= hi)
		{
		    return;//parasirea metodei
	    }

		double mid = a[(lo + hi) / 2];
		while (lo < hi)
		{
	    	while (lo<hi && a[lo] < mid)
	    	{
				lo++;
	    	}
	    	while (lo<hi && a[hi] > mid)
	    	{
				hi--;
	    	}
	    	if (lo < hi)
	    	{
				temp = a[lo];
				a[lo] = a[hi];
				a[hi] = temp;
                lo++;
                hi--;
		    }
	    }
	    if (hi < lo)
	    {
	    	tmp = hi;
	    	hi = lo;
	    	lo = tmp;
		}
		qSortInt(a, lo0, lo);
		qSortInt(a, lo == lo0 ? lo+1 : lo, hi0);

    }
	// sorteaza 2 siruri dupa sortarea crescatoare a primului sir
	/**
	 * Quick sorts two arrays based on ascending sorting of the first array. 
	 * Changes in first array are reflected in changes in second array. Strictly speaking the second array is NOT sorted but 
	 * the element positions are changed as in the first array. 
	 * @param a first array
	 * @param b second array
	 * @param lo0 input lower index
	 * @param hi0 input upper index
	 */
	private static void qSort2(double[] a, double[] b, int lo0, int hi0) {
		double temp = 0.0;
		int tmp = 0;
		int lo = lo0;
		int hi = hi0;

		if (lo >= hi) {
			return;// parasirea metodei
		}

		double mid = a[(lo + hi) / 2];
		while (lo < hi) {
			while (lo < hi && a[lo] < mid) {
				lo++;
			}
			while (lo < hi && a[hi] > mid) {
				hi--;
			}
			if (lo < hi) {
				temp = a[lo];
				a[lo] = a[hi];
				a[hi] = temp;

				temp = b[lo];
				b[lo] = b[hi];
				b[hi] = temp;

				lo++;
				hi--;
			}
		}
		if (hi < lo) {
			tmp = hi;
			hi = lo;
			lo = tmp;
		}
		qSort2(a, b, lo0, lo);
		qSort2(a, b, lo == lo0 ? lo + 1 : lo, hi0);

	}

	//modifica sirul input--crescator
	/**
	 * Sort the array in ascending order. The input array is altered.
	 * @param a the array
	 * @return the sorted array (redundant)
	 */
    public static double[] qSort(double[] a)
    {
	    qSort(a, 0, a.length-1);
	    return a;
	}
    
	// sorteaza 2 siruri dupa sortarea crescatoare a primului sir
	// modifica sirurile input--crescator
    /**
     * Sort 2 arrays at the same time, based on the values stored in the 1st array which is 
     * sorted in ascending order. Useful when we have a series 
     * (x,y) both x, y arrays with y = f(x) and we want to sort x in ascending order and also change  
     * the corresponding y values. The input arrays are altered during the process.
     * @param a 1st array
     * @param b 2bd array
     */
	public static void qSort2(double[] a, double[] b) {
		qSort2(a, b, 0, a.length - 1);
	}
	
	// nu modifica sirul input--crescator
	/**
	* Sort the array in ascending order. The input array is unaltered.
	 * @param a the array
	 * @return the sorted array
	 */
    public static double[] newQSort(double[] a)
    {
		double[] a1 = new double[a.length];
		for (int i=0; i<a.length; i++)
		   a1[i] = a[i];

	    qSort(a1, 0, a1.length-1);
	    return a1;
	}

    /**
     * Sort the array in ascending order. The input array is unaltered.
	 * @param a the array of integers
	 * @return the sorted array
     */
    public static int[] newQSortInt(int[] a)
    {
		int[] a1 = new int[a.length];
		for (int i=0; i<a.length; i++)
		   a1[i] = a[i];

	    qSortInt(a1, 0, a1.length-1);
	    return a1;
	}
    
    /**
     * Find the ranked value from an array. First a new array is computed and sorted. Rank 1 corresponds 
     * to the maximum value. So, findValue(myArray, 1) returns the maximum value.
     * @param a the array
     * @param n the rank
     * @return the ranked value
     */
	public static double findValue(double[] a, int n) {
		if (a.length > 1) {
			double[] a1 = newQSort(a);

			if (n <= a.length && n > 0)
				return a1[a.length - n];
			else
				return a1[0];

		} else
			return a[0];
	}
	
	//============
	private static boolean zero=false;
	@Deprecated
	public static boolean validare()
	{
		return zero;
	}
	@Deprecated
	public static int findHvlPosition(double[] a, int hvl, double zeroVal)
	{
		zero=true;
		double y=0.0;
		if(hvl==1)
		   y=zeroVal/2;
		else
		   y=zeroVal/4;

		int ip=0;
        while(zero)
        {
           if (a[ip]<=y && a[ip+1]>=y)
             break;

           if (a[ip]>=y && a[ip+1]<=y)
             break;

           ip++;
           if (ip==a.length-1)
           {
			   zero=false;
			   break;
		   }
	    }

	    return ip;
	}
	
	 /**
     * Find the ranked index from an array. First a new array is computed and sorted. Rank 1 corresponds 
     * to the maximum value. So, findPosition(myArray, 1) returns the index associated to the maximum value.
     * @param a the array
     * @param n the rank
     * @return the ranked position (index)
     */
	public static int findPosition(double[] a, int n)
	{

        double[] pos = new double[a.length];
        double[] atemp = new double[a.length];

		if (a.length > 1)
		{
			for (int i=0; i<a.length; i++)
			{
				atemp[i]=a[i];
			    pos[i]=i;
			}


		    qSort2(atemp,pos);

		    if (n<=a.length && n>0)
			    return (int)pos[a.length-n];//.intValue();
			else
			    return (int)pos[0];//.intValue();

		}
		else
		    return (int)pos[0];
	}
}
