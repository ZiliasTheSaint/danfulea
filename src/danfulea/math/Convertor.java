package danfulea.math;

import java.math.BigDecimal;
import java.sql.Date;
import java.sql.Timestamp;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.text.NumberFormat;
import java.text.ParseException;
import java.text.SimpleDateFormat;
//import java.util.Date;
import java.util.Locale;

/**
 * The Convertor class is used for handling several number conversions.
 * 
 * @author Dan Fulea, 14 APR. 2011
 * 
 */

public class Convertor {

	/**
	 * Converts a string to int value.
	 * 
	 * @param value
	 *            the string
	 * @return the int value of input string
	 * @throws NumberFormatException
	 * can throws this exception
	 */
	public static int stringToInt(String value) throws NumberFormatException {
		value = value.trim();// Returns a copy of the string, with leading and
								// trailing whitespace omitted.
		if (value.length() == 0) {
			return 0;
		} else {
			return Integer.parseInt(value);
		}
	}

	/**
	 * Converts an int number to string.
	 * 
	 * @param i
	 *            the int value
	 * @return the string representation
	 */
	public static String intToString(int i) {
		return Integer.toString(i);
	}

	/**
	 * Converts an long number to string.
	 * 
	 * @param i
	 *            the long value
	 * @return the string representation
	 */
	public static String longToString(long i) {
		return Long.toString(i);
	}

	/**
	 * Converts a string to long value.
	 * 
	 * @param value
	 *            the string
	 * @return the long value of input string
	 * @throws NumberFormatException
	 * can throws this exception
	 */
	public static long stringToLong(String value) throws NumberFormatException {
		value = value.trim();
		if (value.length() == 0) {
			return 0;
		} else {
			return Long.parseLong(value);
		}
	}

	/**
	 * Converts a string to SQL Timestamp format.
	 * 
	 * @param value
	 *            the string in format: yyyy-[m]m-[d]d hh:mm:ss[.f...]. The
	 *            fractional seconds may be omitted. The leading zero for mm and
	 *            dd may also be omitted.
	 * 
	 * @return the timestamp
	 * @throws IllegalArgumentException
	 *             can throw this exception
	 */
	public static Timestamp stringToTimestamp(String value) throws IllegalArgumentException {
		return Timestamp.valueOf(value);
	}

	/**
	 * Converts a Timestamp to string. The timestamp can be created using new
	 * Timestamp(long millis) where millis is given by current time:
	 * System.currentTimeMillis()
	 * 
	 * @param ts
	 *            the Timestamp
	 * @return the String in format: yyyy-mm-dd hh:mm:ss.fffffffff, where
	 *         ffffffffff indicates nanoseconds.
	 */
	public static String timestampToString(Timestamp ts) {
		return ts.toString();
	}

	/**
	 * Converts a string to SQL Date format.
	 * 
	 * @param value
	 *            the string in format: "yyyy-[m]m-[d]d". The leading zero for
	 *            mm and dd may also be omitted.
	 * @return the Date
	 * @throws IllegalArgumentException
	 *             can throw this exception if the date given is not in the JDBC
	 *             date escape format (yyyy-[m]m-[d]d)
	 */
	public static Date stringToDate(String value) throws IllegalArgumentException {
		return Date.valueOf(value);
	}

	/**
	 * Converts a Date to string. The date can be created using new Date(long
	 * millis) where millis is given by current time: System.currentTimeMillis()
	 * 
	 * @param d
	 *            the date
	 * @return the string in format: yyyy-mm-dd format
	 */
	public static String dateToString(Date d) {
		return d.toString();
	}

	/**
	 * Gets the time in millis from a string representing the date and time in
	 * format: yyyy-MM-dd HH:mm:ss. This should be compared to
	 * System.currentTimeMillis() to calculate time passed. Also it can be used
	 * to get millis in order to build a java SQL Timestamp or Date object.
	 * 
	 * @param timeFormattedStr
	 *            the string representing the date to be processed.
	 * @return the time in millis since referrence date: January 1, 1970,
	 *         00:00:00 GMT
	 * @throws ParseException
	 *             can throw this exception
	 */
	public static long getMillisFromDateTime(String timeFormattedStr) throws ParseException {
		SimpleDateFormat formatter = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss");
		formatter.setLenient(false);

		String oldTime = timeFormattedStr;// e.g: "2012-07-11 10:55:21";
		java.util.Date oldDate = formatter.parse(oldTime);
		long oldMillis = oldDate.getTime();

		return oldMillis;
	}

	/**
	 * Converts a string to boolean (true or false)
	 * 
	 * @param value
	 *            the string
	 * @return the boolean value of input string as true only if the string is
	 *         "true" ignoring case. Otherwise returns false.
	 */
	public static boolean stringToBoolean(String value) {
		return Boolean.parseBoolean(value);// return true only if input is
											// "true" (ignoring case) and false
											// otherwise
	}

	/**
	 * Converts a boolean number to string.
	 * 
	 * @param b
	 *            the boolean value
	 * @return the string representation
	 */
	public static String booleanToString(boolean b) {
		return Boolean.toString(b);
	}

	/**
	 * Converts a string to short (small int) value.
	 * 
	 * @param value
	 *            the string
	 * @return the short value of input string
	 * @throws NumberFormatException
	 * can throws this exception
	 */
	public static short stringToShort(String value) throws NumberFormatException {
		value = value.trim();
		if (value.length() == 0) {
			return 0;
		} else {
			return Short.parseShort(value);
		}
	}

	/**
	 * Converts an short number to string.
	 * 
	 * @param i
	 *            the short value
	 * @return the string representation
	 */
	public static String shortToString(short i) {
		return Short.toString(i);
	}

	/**
	 * Converts a string to float value.
	 * 
	 * @param value
	 *            the string
	 * @return the float value of input string
	 * @throws NumberFormatException
	 * can throws this exception
	 */
	public static float stringToFloat(String value) throws NumberFormatException {
		value = value.trim();
		if (value.length() == 0) {
			return 0;
		} else {
			return Float.parseFloat(value);
		}
	}

	/**
	 * Converts a float number to string.
	 * 
	 * @param f
	 *            the float value
	 * @return the string representation
	 */
	public static String floatToString(float f) {
		return Float.toString(f);
	}

	/**
	 * Converts a string to Big Decimal
	 * @param value
	 * the string
	 * @return the Big Decimal value of String
	 * @throws NumberFormatException
	 * can throws this exception
	 */
	public static BigDecimal stringToBigDecimal(String value) throws NumberFormatException{
		value = value.trim();
		if (value.length() == 0) {
			return BigDecimal.ZERO;
		} else {
			// or Locale.getDefault() but we dont want JVM we want what
			// SQL uses i.e. want System locale:
			// Locale locale = new Locale(System.getProperty("user.language"),
			// System.getProperty("user.country"));
			
			//int point = value.indexOf(".");//Returns the index within this string of the first occurrence of the specified substring.			
			//int point = value.lastIndexOf(".");//last index, if no such character occurs in this string, then -1 is returned
			//int comma = value.lastIndexOf(",");//last index
			//if (point==-1 && comma==-1){//plain number, e.g. 675456
			//	return new BigDecimal(value);
			//} WE HAVE TO USE SEPARATORS SINCE WE CAN HAVE 868.989.989 but no decimal here!!!!!
			//System.out.println("point?="+point+"; comma="+comma);
			
			DecimalFormatSymbols dfs = new DecimalFormatSymbols(new Locale(System.getProperty("user.language"),
					 System.getProperty("user.country")));//Locale.US);
			char groupSeparatorChar = dfs.getGroupingSeparator();
			char decimalSeparatorChar = dfs.getDecimalSeparator();
			String groupSeparator;
			String decimalSeparator;
			
			if(groupSeparatorChar == '.')
	        {
	            groupSeparator = "\\" + groupSeparatorChar;
	        }
	        else
	        {
	            groupSeparator = Character.toString(groupSeparatorChar);	            
	        }

	        if(decimalSeparatorChar == '.')
	        {
	            decimalSeparator = "\\" + decimalSeparatorChar;
	        }
	        else
	        {
	            decimalSeparator = Character.toString(decimalSeparatorChar);
	        }
			
			String fixedString = value.replaceAll(groupSeparator , "");
			fixedString = fixedString.replaceAll(decimalSeparator , ".");
			return new BigDecimal(fixedString);
		}
	}

	/**
	 * Converts a big decimal to string
	 * 
	 * @param bd
	 *            the big decimal
	 * @return the string representation
	 */
	public static String bigDecimalToString(BigDecimal bd) {
		return bd.toString();
	}

	/**
	 * Converts a string to double value.
	 * 
	 * @param value
	 *            the string
	 * @return the double value of input string
	 * @throws NumberFormatException
	 * can throws this exception
	 */
	public static double stringToDouble(String value) throws NumberFormatException {
		value = value.trim();
		if (value.length() == 0) {
			return 0;
		} else {
			return Double.parseDouble(value);
		}
	}

	/**
	 * Converts a double number to string.
	 * 
	 * @param d
	 *            the double value
	 * @return the string representation
	 */
	public static String doubleToString(double d) {
		return Double.toString(d);
	}

	/**
	 * Converts ASCII int value to a String.
	 * 
	 * @param i
	 *            the ASCII integer
	 * @return the string representation
	 */
	public static String asciiToStr(int i) {
		char a[] = new char[1];
		a[0] = (char) i;
		return (new String(a)); // char to string
	}

	/**
	 * Converts a number to a string with significant digits.
	 * 
	 * @param number
	 *            the number to be formatted
	 * @param digits
	 *            the number of digits required!
	 * @return the string representation
	 */
	public static String formatNumber(double number, int digits) {
		NumberFormat nf = NumberFormat.getInstance(
				new Locale(System.getProperty("user.language"),
						 System.getProperty("user.country")));//Locale.US);
		nf.setMinimumFractionDigits(digits);
		nf.setMaximumFractionDigits(digits);
		nf.setGroupingUsed(false);// e.g. no 4,568.02 but 4568.02

		return nf.format(number);
	}

	/**
	 * Converts a number to a string using the scientific format.
	 * 
	 * @param number
	 *            the number
	 * @return the string representation
	 */
	public static String formatNumberScientific(double number) {
		String pattern = "0.###E0";
		DecimalFormatSymbols dfs = new DecimalFormatSymbols(
				new Locale(System.getProperty("user.language"),
						 System.getProperty("user.country")));//Locale.US);
		DecimalFormat nff = new DecimalFormat(pattern, dfs);

		return nff.format(number);
	}
}
