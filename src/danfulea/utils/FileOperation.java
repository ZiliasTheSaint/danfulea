package danfulea.utils;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.util.Date;
import java.util.ResourceBundle;
import java.util.Vector;

/**
 * Class used for handling copy operations! It is thread-safe because the only static member that matters, i.e. the buffer,
 * runs in synchronized block! So, this methods can run in multiple threads without causing a mess. 
 * It is based on some tutorials/examples freely available on the web.
 * 
 * @author Dan Fulea, 15 APR 2011
 * 
 */
public class FileOperation {
	private static final int BUFF_SIZE = 8 * 1024 * 1000;// 8 MB RAM used...should be
	// enough! 1024 B= 1 KB!
	private static final byte[] buffer = new byte[BUFF_SIZE];

	private static final String BASE_RESOURCE_CLASS = "danfulea.resources.DanfuleaResources";
	private static ResourceBundle resources = ResourceBundle.getBundle(BASE_RESOURCE_CLASS);
	public static MessageRetriever mr;

	/**
	 * The best way to copy files. Any NIO method (e.g. file channels) corrupts
	 * data. Tested!<br>
	 * Success copy for 3.5 Gb file, 400Mb file!! The copied files were opened
	 * using specific<br>
	 * softwares and we found data are NOT corrupted!!!! <br>
	 * 
	 * @param from
	 *            a string representing the source filename
	 * @param to
	 *            a string representing the destination filename
	 * @throws IOException
	 * 			  can throw this exception
	 */
	public static void copyFile(String from, String to) throws IOException {
		/*
		 * Assign null to variables that you are no longer using. This will, if
		 * there are no other references, allow this memory to be recycled
		 * (garbage collected). Because local variables are deallocated when a
		 * method returns, they are less of a problem.
		 * 
		 * Variables with the longest lifetime are static variables. Static
		 * variables are destroyed when the class is unloaded. This can occur
		 * when it's class loader is garbage collected.
		 * 
		 * ....................................... So any local variables such
		 * as source_file = new File(from) which
		 * are presented here (this method) are auto-marked for garbage
		 * collector (and subject to be freed from computer RAM) as soon as this
		 * method has return!!! The real process of garbage collecting can occur
		 * any time and is taken care auto (smart) from within java (JVM). Users
		 * don't need to worry about and take extra-caution measures! Anyway, we
		 * can optimize the code because we'll never know when memory
		 * reallocation has taken place in order to speed up the program.
		 * Therefore we make variable buffer static because each
		 * FileOperation.copy...call will use the same memory...and variable
		 * buffer does the hard work here!
		 */
		String str = "";
		File source_file = new File(from);
		File destination_file = new File(to);
		// making sure the source exists!
		if (!source_file.exists() || !source_file.isFile()) {
			str = resources.getString("copy.sourceExists");
			throw new IOException(str + from);
		}
		if (!source_file.canRead()) {
			str = resources.getString("copy.sourceRead");
			throw new IOException(str + from);
		}
		// making sure destination exists and the clone file can be written
		// there!
		File parentdir = parent(destination_file);// System.out.println(parentdir.getAbsolutePath());
		if (!parentdir.exists()) {
			str = resources.getString("copy.destParentExists");
			throw new IOException(str + parentdir.getAbsolutePath());
		}
		if (!parentdir.canWrite()) {
			str = resources.getString("copy.destParentWritable");
			throw new IOException(str + parentdir.getAbsolutePath());
		}
		// if (getFileExtension(to).compareTo("")==0){
		// not used because we can copy files without extension!
		// let the system (java) error to be caught!
		// str = resources.getString("copy.destFile");
		// throw new IOException(str + to);
		// }
		// now we can copy!
		double dCopied = 0.0;
		double dfileSize = (long) new File(from).length() / 1.0;// (1024*1000);//MB
		int ibatch = 0;
		int inc = 10;
		// double
		long timeIn = new Date().getTime();
		InputStream in = null;
		OutputStream out = null;
		try {
			in = new FileInputStream(from);
			out = new FileOutputStream(to);
			while (true) {
				synchronized (buffer) {//synchronized means can be used one at the time! Force all other threads to wait for this to finish!
					int amountRead = in.read(buffer);
					if (amountRead == -1) {
						break;
					}
					out.write(buffer, 0, amountRead);
					dCopied = dCopied + amountRead;// /(1024*1000);//in MB!
					if (dCopied >= dfileSize) {
						// str = "Copyed: " + "100" + " %";
						str = resources.getString("copy.copyed.100");
						if (mr != null)
							mr.printSequence(str);
					} else if (dCopied >= dfileSize * (ibatch + inc) / 100) {
						ibatch = ibatch + inc;
						// str = "Copyed: " + ibatch + " %";
						str = resources.getString("copy.copyed");
						str = str + ibatch + resources.getString("copy.proc");
						if (mr != null)
							mr.printSequence(str);
					}

				}//end synchronized!!!
			}
		} finally {
			if (in != null) {
				in.close();
			}
			if (out != null) {
				out.close();
			}
		}
		long timeOut = new Date().getTime();
		str = resources.getString("copy.time");
		// str = str + (long) new File(from).length() + " => "
		str = str + dfileSize + " => " + (timeOut - timeIn);
		// System.out.println("copy Time [ms] for a file of size [bytes] :"
		// System.out.println(str
		// + (long) new File(from).length() + " => " + (timeOut - timeIn));
		// System.out.println(str);
		if (mr != null)
			mr.printSequence(str);
	}

	/**
	 * Copies the content of a folder (recursively) into another folder
	 * 
	 * @param source_name
	 *            the source folder
	 * @param dest_name
	 *            the destination folder
	 * @throws IOException
	 *  			  can throw this exception
	 */
	public static void copyDir(String source_name, String dest_name) throws IOException {
		String str = "";
		File source_file = new File(source_name);
		File destination_file = new File(dest_name);

		// making sure the source exists and is readable!
		if (!source_file.exists() || !source_file.isDirectory()) {
			str = resources.getString("copy.sourceDirExists");
			throw new IOException(str + source_name);
		}
		if (!source_file.canRead()) {
			str = resources.getString("copy.sourceDirRead");
			throw new IOException(str + source_name);
		}
		// making sure destination exists and it is a folder
		// Always overwrites!
		if (!destination_file.exists() || !destination_file.isDirectory()) {
			str = resources.getString("copy.destExists");
			throw new IOException(str + dest_name);
		}

		// String[] s = getFiles(source_file);
		Vector<String> s = ScanDisk.getFiles(source_file);
		for (int i = 0; i < s.size(); i++) {// length; i++) {
			File dest = new File(s.elementAt(i));// s[i]);
			String new_dest_name = dest_name + System.getProperty("file.separator") + dest.getName();
			copyFile(s.elementAt(i), new_dest_name);// s[i], new_dest_name);
		}

		// String[] sdir = getDirs(source_file);
		Vector<String> sdir = ScanDisk.getDirs(source_file);
		for (int i = 0; i < sdir.size(); i++) {// length; i++) {
			File sdirfile = new File(sdir.elementAt(i));// [i]);
			String new_dest_name2 = dest_name + System.getProperty("file.separator") + sdirfile.getName();
			File newDir = new File(new_dest_name2);
			newDir.mkdir();

			source_name = sdir.elementAt(i);// [i];

			copyDir(source_name, new_dest_name2);// recursive
		}
	}

	/**
	 * Copy the whole folder into a parent folder creating a new folder having
	 * the same name as the source folder.
	 * 
	 * @param source_name
	 *            the source folder
	 * @param dest_name
	 *            the destination folder (parent)
	 * @throws IOException
	 * 			  can throw this exception
	 */
	public static void copyWholeDir(String source_name, String dest_name) throws IOException {
		String str = "";
		File source_file = new File(source_name);
		File destination_file = new File(dest_name);

		// making sure the source exists and is readable!
		if (!source_file.exists() || !source_file.isDirectory()) {
			str = resources.getString("copy.sourceDirExists");
			throw new IOException(str + source_name);
		}
		if (!source_file.canRead()) {
			str = resources.getString("copy.sourceDirRead");
			throw new IOException(str + source_name);
		}
		// making sure destination exists and it is a folder
		// Always overwrites!
		if (!destination_file.exists() || !destination_file.isDirectory()) {
			str = resources.getString("copy.destExists");
			throw new IOException(str + dest_name);
		}

		String new_dest_name = dest_name + System.getProperty("file.separator") + source_file.getName();
		File newDir = new File(new_dest_name);
		newDir.mkdir();

		copyDir(source_name, new_dest_name);
	}

	/**
	 * Retrieving the parent of the specified file.
	 * 
	 * @param f
	 *            the file
	 * @return the file directory (folder)
	 */
	public static File parent(File f) {
		String dirname = f.getParent();
		if (dirname == null) {
			if (f.isAbsolute())
				return new File(File.separator);
			else
				return new File(System.getProperty("user.dir"));

		}
		return new File(dirname);
	}

	/**
	 * Retrieving the file extension.
	 * 
	 * @param file
	 *            a String representation of file
	 * @return the file extension
	 */
	public static String getFileExtension(String file) {
		int idx = file.lastIndexOf(".");
		// handles unix hidden files and files without an extension
		if (idx < 1) {
			return "";
		}
		return file.substring(idx+1);//without point
	}
	
	/**
	 * Retrieving the file extension length including '.' character.
	 * 
	 * @param file
	 *            a String representation of file
	 * @return the extension length
	 */
	public static int getFileExtensionLength(String file) {
		int idx = file.lastIndexOf(".");
		// handles unix hidden files and files without an extension
		if (idx < 1) {
			return 0;
		}
		return file.length()-idx;// +1 is included in length which is not last index;
	}
}
