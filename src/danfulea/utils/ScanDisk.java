package danfulea.utils;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DateFormat;
import java.util.Calendar;
import java.util.Date;
import java.util.Vector;

/**
 * 
 * ScanDisk utilities. It is designed to look for a particular file/folder
 * inside a computer<br>
 * or scan entire computer and generate a text file with file/folder structure!
 * <br>
 * 
 * @author Dan Fulea, 01 DEC. 2003
 * 
 */
public class ScanDisk {// implements MessageRetriever{

	/**
	 * Stores the number of folders
	 */
	private static int contorDirs;//an int (primitive type) has default value of 0 even if not explicitly set!!!

	/**
	 * Store the number of files and folders
	 */
	private static int contorFilesAndDirs;

	/**
	 * The interface for displaying the results
	 */
	public static MessageRetriever mr;

	/**
	 * Aborts scan execution!
	 */
	public static boolean stopB = false;
	//this class is not a thread, does not implement runnable so stopB doesn't have to be 
	//volatile as in good Java practice when communication between threads are of concern.
	//this stopB is like a volatile exitThreadRun flag in the sense that it is required
	//in methods which are time-consuming due to intensive calculations (see recursive methods
	//nelow). Again this is not a thread just a time-consuming collection of methods.
	
	
	/**
	 * Get the machine roots (drives), i.e. C:\, D:\...
	 * 
	 * @return the array of roots
	 */	
	public static String[] getDrives() {
		String[] s = { "" };
		File[] f = File.listRoots();// List the available filesystem roots.

		if (f != null) {
			s = new String[f.length];
			for (int i = 0; i < f.length; i++) {
				s[i] = f[i].toString();				
			}
		}

		return s;
	}

	/**
	 * Get files and folders from a root folder (e.g. C:\ or C:\Windows).
	 * 
	 * @param root
	 *            the root folder
	 * @return the files and folders (sub-folders not included).
	 */
	public static Vector<String> getFilesAndDirs(File root) {
		
		Vector<String> v = new Vector<String>();
		
		File[] f = root.listFiles();// list all files and folders

		if (f != null) {			
			for (int i = 0; i < f.length; i++) {
				if (f[i]!=null)
					v.addElement((String) f[i].toString());			
			}
		}

		return v;
	}
	
	/**
	 * Returns the folders from a root folder. It also counts folders and
	 * folders plus files.
	 * 
	 * @param root
	 *            the root folder
	 * @return the folders (sub-folders not included).
	 */
	public static Vector<String> getDirs(File root) {
		
		Vector<String> v = new Vector<String>();

		File[] f = root.listFiles();// list all files and folders
		
		if (f != null) {
			for (int i = 0; i < f.length; i++) {
				if (f[i]!=null && f[i].isDirectory()){ //the list method in File class may issue a null!!!
					v.addElement((String) f[i].toString());
					contorDirs++;
				}
				contorFilesAndDirs++;
			}
		}
		
		return v;		
	}
	
	/**
	 * Gets files inside a root.
	 * 
	 * @param root
	 * 			the root folder
	 * @return the files contained in root folder
	 */
	public static Vector<String> getFiles(File root) {
		
		Vector<String> v = new Vector<String>();

		File[] f = root.listFiles();// list all files and folders
		//can contain null files=>the reason of using Vectors!
				
		if (f != null) {
			for (int i = 0; i < f.length; i++) {
				if (f[i]!=null && f[i].isFile()) {
					v.addElement((String) f[i].toString());
				}
			}
		}
		
		return v;
	}

	/**
	 * Recursively scan all folders from a root, including sub-folders. While
	 * scanning, the absolute path of folders is displayed within any class
	 * objects which implements {@link danfulea.utils.MessageRetriever
	 * MessageRetriever} interface.
	 * <p>
	 * For saving the scan result in a text file (located in application
	 * folder), it also uses the FileWriter to write the absolute path of
	 * folders.
	 * </p>
	 * 
	 * @see danfulea.utils.ScanDisk#scanDirs(File, String, Vector, boolean)
	 * 
	 * @param root
	 *            the root
	 * @param fw
	 *            the FileWriter
	 * @throws IOException
	 *             can throw this exception
	 */
	private static void scanDirs(File root, FileWriter fw) throws IOException {
		
		Vector<String> s = getDirs(root);
		for (int i = 0; i < s.size(); i++) {			
			if (mr != null)
				mr.printSequence(s.elementAt(i));

			if (stopB)
				return;// force exit!!

			fw.write(s.elementAt(i) + "\n");
			File f = new File(s.elementAt(i));
			scanDirs(f, fw);
		}
	}

	/**
	 * Recursively scan all folders from a root, including sub-folders, in order
	 * to find a file or a folder.
	 * 
	 * @param root
	 *            the root
	 * @param target
	 *            the file we are looking for
	 * @param stock
	 *            the stock vector to store the results
	 * @param case_sensitiv
	 *            if search is case-sensitive or not
	 * @see danfulea.utils.ScanDisk#scanDirs(File, FileWriter)
	 */
	private static void scanDirs(File root, String target, Vector<String> stock, boolean case_sensitiv) {
		
		Vector<String> s = getDirs(root);
		findFile(root, target, stock, case_sensitiv);
		for (int i = 0; i < s.size(); i++) {
			File f = new File(s.elementAt(i));

			if (stopB)
				return;// force exit!!

			scanDirs(f, target, stock, case_sensitiv);
		}
	}	

	/**
	 * Finds if a specific target file or folder is located inside the root
	 * folder. Usage of wildcard characters (which is *) is implemented!
	 * 
	 * @param root
	 *            the root folder
	 * @param target
	 *            the file or folder we are searching for
	 * @param stock
	 *            the stock vector to store the results
	 * @param case_sensitiv
	 *            if search is case-sensitive or not
	 */
	private static void findFile(File root, String target, Vector<String> stock, boolean case_sensitiv) {
		//some tests first
		String[] s = root.list();// Returns an array of strings naming the files
									// and directories
		File[] f = root.listFiles();//list and listFiles do exactly the same=>redundant

		if (s == null)
			return;
		if (f == null)
			return;
		//test passed
		
		boolean findB = false;

		String queryStr = target;
		queryStr = queryStr.replaceAll("\\*", ".*");// replace one regexp with
													// another.
		
		//Reserved characters, aka meta characters are command characters that have 
		//special meaning in regexes must be quoted when you mean them literally, as 
		//just characters. This does not mean you must enclose them in quotation marks,
		//but rather you must specially mark them as meant literally by preceding them 
		//with a \. e.g. \- \+ \?. If you are unsure, quote. It won’t hurt to quote 
		//punctuation that does not need it. However, Don’t quote : in since \:… has
		//special meaning.
		//Unfortunately, the regex people used the same quoting character \ as the
		//designers of Java did for String literals. In a non-regex Java String literal,
		//every literal \ must be doubled. In a regex every literal \ must be doubled. 
		//So when you express a regex as a 	Java String literal, every literal \ must
		//be quadrupled! and written as \\\\.
		
		//So, we want to replace * with .* using String.replaceAll(String regexp, 
		//String replacement)
		//.* is simple=>it is a regexp given by the STring ".*" meaning what we want: 
		//find all.
		//* however is already a regexp with a different meaning of what we want (we want
		//find all). So, we must mark it as meant literally *, so we have to quote it
		//as \*. But in order to make it a regexp (to be used with replace all) we must 
		//double \ so we finally have "\\*".
		
		// if queryStr do not contain wildcard * then normal match is performed!!
		
		for (int i = 0; i < s.length; i++) {
			if (s[i]!=null){//because list method in File class may issue a null!!!
				if (case_sensitiv) {
		//matches() = tells whether or not this string matches the given regular expression
		//contains()=Returns true if and only if this string contains the specified sequence of char values.
					if (s[i].matches(queryStr)){// || s[i].contains(queryStr)) {
						String ss = f[i].toString();
						stock.addElement(ss);
						findB = true;
					}
				} else {
				
					if (s[i].matches("(?i:" + queryStr + ")")) //|| 
						//s[i].contains(queryStr))//CONTAINS IS CASE SENSITIVE!!!!
		//no need for contains...using *ming* is the same and we already allow wildcards!!!!
					{
						String ss = f[i].toString();
						stock.addElement(ss);
						findB = true;
					}
				}
			}

		}
		
		if (findB) {
			if (mr != null)
				mr.printSequence("found: " + stock.elementAt(stock.size() - 1));// .toString());
		}
		//if contains ignore case is needed, we can use Region Matches:
		//str is the string containing searchStr.:
		//final int len = searchStr.length();
	    //final int max = str.length() - len;
	    //for (int i = 0; i <= max; i++) {
	    //    if (str.regionMatches(true, i, searchStr, 0, len)) {//TRUE here means IGNORE CASE
	     //       return true;
	     //   }
	    //}
	    //return false;
	}

	/**
	 * 
	 * Begins the scan of root folder. Creates scan result file, ScanResult.txt,
	 * in application folder. In addition, the results are passed to any class
	 * object which implements {@link danfulea.utils.MessageRetriever
	 * MessageRetriever} interface.
	 * 
	 * @param root
	 *            the root folder
	 * @throws IOException
	 *             can throw this exception
	 */
	public static void createScanResultFileInCurrentDir(File root) throws IOException{//throws Exception {
		
		resetContors();
		
		File f = new File(System.getProperty("user.dir") + 
				System.getProperty("file.separator") +
				"ScanResult.txt");
		f.createNewFile();
		FileWriter fw = new FileWriter(f);

		fw.write("Folders structure: " + "\n");
		if (mr != null)
			mr.printSequence("Folders structure: ");

		Calendar CALENDAR = Calendar.getInstance();
		Date date = CALENDAR.getTime();
		String myString = DateFormat.getDateTimeInstance().format(date);

		fw.write("Start scanning at " + myString + "\n");
		if (mr != null)
			mr.printSequence("Start scanning at " + myString);

		if (isAvailable(root))
			scanDirs(root, fw);
		fw.write("Number of files and folders = " + contorFilesAndDirs + "\n");
		if (mr != null)
			mr.printSequence("Number of files and folders = " + contorFilesAndDirs);

		fw.write("Number of folders = " + contorDirs + "\n");
		if (mr != null)
			mr.printSequence("Number of folders = " + contorDirs);

		fw.write("Number of files = " + getFilesNum() + "\n");
		if (mr != null)
			mr.printSequence("Number of files = " + getFilesNum());

		Calendar CAL = Calendar.getInstance();
		Date d = CAL.getTime();
		String s = DateFormat.getDateTimeInstance().format(d);

		fw.write("Stop scanning at " + s + "\n");
		if (mr != null)
			mr.printSequence("Stop scanning at " + s);

		fw.close();
	}

	/**
	 * Checks if a root folder contains files or not
	 * 
	 * @param root
	 *            the root folder
	 * @return true if available, false otherwise
	 */

	public static boolean isAvailable(File root) {
		boolean b = true;
		try {
			File[] f = root.listFiles();//this can also contains null files!!!
			@SuppressWarnings("unused")
			int available = f.length;// required for forcing the exception if any!
		} catch (NullPointerException e) {
			b = false;
		}
		return b;
	}

	/**
	 * 
	 * Begins the computer scan. Creates scan result file, ScanResult.txt, in
	 * application folder. In addition, the results are passed to any class
	 * object which implements {@link danfulea.utils.MessageRetriever
	 * MessageRetriever} interface.
	 * 
	 * @throws IOException
	 *             can throw this exception
	 */
	public static void createComputerScanResultFileInCurrentDir() throws IOException{//throws Exception {
		//first thing first: reset counts
		resetContors();
		
		File f = new File(System.getProperty("user.dir") + 
				System.getProperty("file.separator") + 
				"ScanResult.txt");
		f.createNewFile();
		FileWriter fw = new FileWriter(f);

		fw.write("Folders structure: " + "\n");
		if (mr != null)
			mr.printSequence("Folders structure: ");

		File[] ff = File.listRoots();
		Calendar CALENDAR = Calendar.getInstance();
		Date date = CALENDAR.getTime();
		String myString = DateFormat.getDateTimeInstance().format(date);

		fw.write("Start scanning at " + myString + "\n");
		if (mr != null)
			mr.printSequence("Start scanning at " + myString);

		for (int i = 0; i < ff.length; i++) {
			//System.out.println("" + ff[i].toString());
			if (isAvailable(ff[i]))
				scanDirs(ff[i], fw);
		}
		fw.write("Number of files and folders = " + contorFilesAndDirs + "\n");
		if (mr != null)
			mr.printSequence("Number of files and folders = " + contorFilesAndDirs);

		fw.write("Number of folders = " + contorDirs + "\n");
		if (mr != null)
			mr.printSequence("Number of folders = " + contorDirs);

		fw.write("Number of files = " + getFilesNum() + "\n");// contorFiles +
																// "\n");
		if (mr != null)
			mr.printSequence("Number of files = " + getFilesNum());

		Calendar CAL = Calendar.getInstance();
		Date d = CAL.getTime();
		String s = DateFormat.getDateTimeInstance().format(d);

		fw.write("Stop scanning at " + s + "\n");
		if (mr != null)
			mr.printSequence("Stop scanning at " + s);

		fw.close();
	}

	/**
	 * 
	 * Begins the scan of root folder in order to find a file or a folder. The
	 * results are passed to any class object which implements
	 * {@link danfulea.utils.MessageRetriever MessageRetriever} interface.
	 * 
	 * @param root
	 *            The root folder to begin the search
	 * @param target
	 *            the file or folder we are looking for
	 * @param case_sensitiv
	 *            if the search is case sensitive or not.
	 */
	public static void scanForFile(File root, String target, boolean case_sensitiv) {
		Calendar CALENDAR = Calendar.getInstance();
		Date date = CALENDAR.getTime();
		String myString = DateFormat.getDateTimeInstance().format(date);
		if (mr != null)
			mr.printSequence("Start scanning at " + myString);

		Vector<String> v = new Vector<String>();
		if (isAvailable(root))
			scanDirs(root, target, v, case_sensitiv);

		Calendar CAL = Calendar.getInstance();
		Date d = CAL.getTime();
		String ss = DateFormat.getDateTimeInstance().format(d);
		if (mr != null)
			mr.printSequence("Stop scanning at " + ss);

		// return s;
	}

	/**
	 * 
	 * Begins the whole computer scan in order to find a file or a folder. The
	 * results are passed to any class object which implements
	 * {@link danfulea.utils.MessageRetriever MessageRetriever} interface.
	 * 
	 * @param target
	 *            the file we are looking for
	 * @param case_sensitiv
	 *            if search is case sensitive or not
	 */
	public static void scanComputerForFile(String target, boolean case_sensitiv) {
		Calendar CALENDAR = Calendar.getInstance();
		Date date = CALENDAR.getTime();
		String myString = DateFormat.getDateTimeInstance().format(date);
		if (mr != null)
			mr.printSequence("Start scanning at " + myString);
		
		Vector<String> v = new Vector<String>();
		File[] ff = File.listRoots();
		for (int i = 0; i < ff.length; i++) {
			if (isAvailable(ff[i]))
				scanDirs(ff[i], target, v, case_sensitiv);
		}

		Calendar CAL = Calendar.getInstance();
		Date d = CAL.getTime();
		String ss = DateFormat.getDateTimeInstance().format(d);
		if (mr != null)
			mr.printSequence("Stop scanning at " + ss);
		// return s;
	}

	/**
	 * 
	 * @return the folder counts
	 */
	public static int getDirsNum() {
		return contorDirs;
	}

	/**
	 * 
	 * @return the files counts
	 */
	public static int getFilesNum() {
		// return contorFiles;
		return contorFilesAndDirs - contorDirs;
	}

	/**
	 * 
	 * @return the folder plus file counts
	 */
	public static int getFilesAndDirsNum() {
		return contorFilesAndDirs;
	}

	/**
	 * reset counting
	 */
	public static void resetContors() {
		contorDirs = 0;		
		contorFilesAndDirs = 0;
	}

	// =================================================TESTING
	/*public static void main(String[] args) {
		new ScanDisk();
		String[] array = getDrives();

		// for (int i = 0; i<array.length; i++)
		// System.out.println(array[i]);

		File root = new File("C:\\");		
		Vector<String> v = getFilesAndDirs(root);//getFiles(root);//getDirs(root);
		for(int i = 0; i<v.size(); i++)
			System.out.println(v.elementAt(i));
		
		//scanForFile(root, "*mingw*", false);
	}
	
	public ScanDisk(){
		mr = this;
	}
	
	@Override
	public void printSequence(String s) {
		// TODO Auto-generated method stub
		System.out.println(s);
	}*/
	// ================================================END TESTING

	
}
