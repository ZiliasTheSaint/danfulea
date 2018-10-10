package danfulea.math;

import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.io.UnsupportedEncodingException;
import java.util.Vector;

import danfulea.utils.MessageRetriever;

/**
 * Computing primes and factorial. Not super efficient (no super fancy
 * algorithms) but gets the job done! A class implementing MessageRetriever 
 * interface can be used to display results (e.g. inside a JTextArea).
 * 
 * @author Dan Fulea, 12 AUG. 2012
 *
 */
public class PrimesAndFactorial {

	//public static Vector<Double> pV;
	public static MessageRetriever mr;

	/**
	 * Compute factorial in sequence. This is the core function.
	 * 
	 * @param start
	 *            from start
	 * @param n
	 *            to n
	 * @return the factorial
	 */
	private static long factorial(long start, long n) {
		long i;
		if (n <= 16) {
			long r = start;
			for (i = start + 1; i < start + n; i++)
				r *= i;
			return r;
		}
		i = n / 2;
		return factorial(start, i) * factorial(start + i, n - i);
	}

	/**
	 * Compute factorial of given n.
	 * 
	 * @param n
	 *            the long n
	 * @return the factorial
	 */
	public static long factorial(long n) {
		return factorial(1, n);
	}

	/**
	 * Checks if a long number is prime or not.
	 * 
	 * @param n
	 *            the number
	 * @return true if prime
	 */
	public static boolean isPrime(long n) {
		// we don’t have to try out all integers from 2 to n.
		// It is sufficient up to n/2 because n cannot be
		// divided by a number greater than its half!

		// Also, we only have to go up to the square root of n. Because:
		// if n isn’t prime it can be represented as p*q = n
		// if p > sqrt(n), that will mean that q < sqrt(n).

		if (n == 2)
			return true;
		if (n % 2 == 0)
			return false;// 2 is the only even prime number
		for (long i = 3; i * i <= n; i += 2) {
			if (n % i == 0)
				return false;
		}
		return true;
	}

	/**
	 * Retrieves prime numbers up to nmax. All primes are stored in pV Vector.
	 * 
	 * @param nmax
	 *            the end point for searching primes
	 * @param print
	 *            true if want to print results in primes.txt file in
	 *            application folder.
	 * @return the pV vector.
	 * @throws FileNotFoundException
	 *             can throw this exception
	 * @throws UnsupportedEncodingException
	 *             can throw this exception
	 */
	public static Vector<Double> getPrimes(long nmax, boolean print) throws FileNotFoundException, UnsupportedEncodingException {

		PrintWriter writer = new PrintWriter("primes.txt", "UTF-8");
		String str = "";
		double pp = 2.0;// first prime =2
		long index = 1;
		Vector<Double> pV = new Vector<Double>();
		pV.addElement(new Double(pp));

		str = " #" + index + "; " + pp;
		if (mr != null)
			mr.printSequence(str);
		// System.out.println(str);

		if (print)
			writer.println(str);

		pp = pp + 1.0;// 2nd prime = 3
		pV.addElement(new Double(pp));
		index++;

		str = " #" + index + "; " + pp;
		if (mr != null)
			mr.printSequence(str);
		// System.out.println(str);

		if (print)
			writer.println(str);

		boolean test = false;
		while (pp < nmax) {
			pp = pp + 2.0;// 5,7,...
			test = true;
			double sqrtpp = Math.sqrt(pp);
			for (int i = 0; i < pV.size(); i++) {// we need the vector here.
				double a = pV.elementAt(i).doubleValue();// for checkings
				if (a > sqrtpp)
					break;
				if (pp % a == 0) {
					test = false;
					break;
				}
			}
			if (test) {
				pV.addElement(new Double(pp));// so might as well make it public
												// static
				index++;

				str = " #" + index + "; " + pp;
				if (mr != null)
					mr.printSequence(str);
				// System.out.println(str);

				if (print)
					writer.println(str);
			}
		}

		writer.close();
		return pV;
	}
	
	// ========TESTING===============
	/*public static void main(String[] args) {
		// long f = PrimesAndFactorial.factorial(17);
		// System.out.println("fact = "+f);
		long startTime = System.currentTimeMillis();
		try {
			Vector<Double> pv = PrimesAndFactorial.getPrimes(1000, false);
			for (int i = 0; i < pv.size(); i++) {
				System.out.println(i + "; " + pv.elementAt(i));
			}
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (UnsupportedEncodingException e) { // TODO Auto-generated catch
													// block
			e.printStackTrace();
		}
		
		System.out.println("31? "+isPrime(31));
		long endTime = System.currentTimeMillis();
		String deltaT = TimeUtilities.timeElapsed(startTime, endTime);
		System.out.println(deltaT);
	}
	*/
	// ========END====================

}
