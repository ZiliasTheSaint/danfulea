package danfulea.utils;

import java.util.Calendar;

import danfulea.math.Convertor;

/**
 * Time utilities class.
 * 
 * @author Dan Fulea, 28 APR. 2011
 * 
 */
public class TimeUtilities {

	/**
	 * global member representing the day
	 */
	private int iday = 0;

	/**
	 * global member representing the month
	 */
	private int imonth = 0;

	/**
	 * global member representing the year
	 */
	private int iyear = 0;
	
	/**
	 * global member representing the day as String
	 */
	private String idayS = "";

	/**
	 * global member representing the month as String
	 */
	private String imonthS = "";

	/**
	 * global member representing the year as String
	 */
	private String iyearS = "";
	
	/**
	 * Constructor of a TimeUtilities object set as today.
	 */
	public TimeUtilities(){
		today();
	}
	
	/**
	 * Constructor of a TimeUtilities object set at a given date.
	 * @param day
	 * the day
	 * @param month
	 * the month
	 * @param year
	 * the year
	 */
	public TimeUtilities(int day, int month, int year){
		setDate(day, month, year);
	}
	
	/**
	 * Constructor based on formatted date yyyy-mm-dd.
	 * @param s
	 * the string encoding the date
	 */
	public TimeUtilities(String s){
		unformatDate(s);
	}
	
	/**
	 * 
	 * @return the day
	 */
	public int getDay(){
		return iday;
	}

	/**
	 * 
	 * @return the month
	 */
	public int getMonth(){
		return imonth;
	}
	
	/**
	 * 
	 * @return the year
	 */
	public int getYear(){
		return iyear;
	}
	
	//-------------------
	/**
	 * 
	 * @return the day as String
	 */
	public String getDayS(){
		return idayS;
	}

	/**
	 * 
	 * @return the month as String
	 */
	public String getMonthS(){
		return imonthS;
	}
	
	/**
	 * 
	 * @return the year as String
	 */
	public String getYearS(){
		return iyearS;
	}
	
	/**
	 * Sets the date.
	 * 
	 * @param day
	 *            the day
	 * @param month
	 *            the month
	 * @param year
	 *            the year
	 */
	public void setDate(int day, int month, int year) {
		iday = day;
		imonth = month;
		iyear = year;
		
		idayS = Convertor.intToString(iday);
		if (iday < 10)
			idayS = "0" + idayS;
		imonthS = Convertor.intToString(imonth);
		if (imonth < 10)
			imonthS = "0" + imonthS;
		iyearS = Convertor.intToString(iyear);
	}

	/**
	 * Sets date as today.
	 */
	public void today() {
		Calendar cal = Calendar.getInstance();
		iday = cal.get(Calendar.DAY_OF_MONTH);
		imonth = cal.get(Calendar.MONTH) + 1;// 0 index; January=0
		iyear = cal.get(Calendar.YEAR);
		
		idayS = Convertor.intToString(iday);
		if (iday < 10)
			idayS = "0" + idayS;
		imonthS = Convertor.intToString(imonth);
		if (imonth < 10)
			imonthS = "0" + imonthS;
		iyearS = Convertor.intToString(iyear);
	}

	/**
	 * Encodes the date using yyyy-mm-dd format, which is optimal for sorting!
	 * 
	 * @return the String representation of date
	 */
	public String formatDate() {
		idayS = Convertor.intToString(iday);
		if (iday < 10)
			idayS = "0" + idayS;

		imonthS = Convertor.intToString(imonth);
		if (imonth < 10)
			imonthS = "0" + imonthS;

		iyearS = Convertor.intToString(iyear);

		return iyearS + "-" + imonthS + "-" + idayS;
	}

	/**
	 * Decodes the string representation of date and sets the global members
	 * accordingly.
	 * 
	 * @param s
	 *            the string containing the encoded date.
	 */
	public void unformatDate(String s) {
		if (s.equals("")) {
			return;// no formated data
		}
		String[] result = s.split("-");

		idayS = result[2];
		imonthS = result[1];
		iyearS = result[0];

		iyear = Convertor.stringToInt(iyearS);
		imonth = Convertor.stringToInt(imonthS);
		iday = Convertor.stringToInt(idayS);
	}

	/**
	 * Compute time elapsed from start time to end time.
	 * 
	 * @param st
	 *            start time in millis taken from System.currentTimeMillis() when job starts
	 * @param et
	 *            end time in millis taken from System.currentTimeMillis() when job ends
	 * @return time elapsed
	 */
	public static String timeElapsed(long st, long et) {
		int delta = (new Long(et - st)).intValue();// ms
		String times = "";

		int sec = delta / 1000;// impartire intreaga->catul!!!
		int milis = delta % 1000;// restul impartirii intregi!!
		if (sec > 60) {
			int min = sec / 60;
			sec = sec % 60;
			if (min > 60) {
				int h = min / 60;
				min = min % 60;
				if (h > 24) {
					int z = h / 24;
					h = h % 24;
					times = z + " days " + h + " h, " + min + " min, " + sec + " sec, " + milis + " milis";
				} else {
					times = h + " h, " + min + " min, " + sec + " sec, " + milis + " milis";
				}
			} else {
				times = min + " min, " + sec + " sec, " + milis + " milis";
			}
		} else {
			times = sec + " sec, " + milis + " milis";
		}

		String seqStr = "******************************************" + " \n" +
		"Time elapsed: " + times;// +
																										// "
																										// \n";
		//System.out.println(seqStr);
		return seqStr;
	}

}
