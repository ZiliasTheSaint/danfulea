package danfulea.crypto;

import java.io.File;
import java.io.FileInputStream;
import java.io.InputStream;
import java.security.MessageDigest;

/**
 * Utility class for checksum. It is thread-safe.
 * <p>
 * Based on various code posted on Internet. For example:<br>
 * Posted by Vlad at Wednesday, August 01, 2007 at
* <a href="http://vyshemirsky.blogspot.com/2007/08/computing-md5-digest-checksum-in-java.html">
* http://vyshemirsky.blogspot.com/2007/08/computing-md5-digest-checksum-in-java.html</a>
* <p>
 * @author Dan Fulea, 13 JUL. 2010
 *
 */
public class ChecksumUtilities {

	
	private boolean checkSuccess = true;//using static identifier =>not always thread safe
	//not always because it depends on how the static class (have static members and 
	//methods) is called and how threads are build.
	//In general, each thread has its own stack, and so its own copy of variables it 
	//can access. When the thread is created, it copies the value of 
    //all accessible variables in its own memory (including static one) so it is thread 
	//safe as long as you don't explicitly set threads to talk each other (e.g. using a 
	//"common" volatile variable) or using/updating same resources.
	//anyway, use static only when really needed for efficiency. In some cases using 
	//static is more efficient since it avoids object creation.
	
	private String checkError = "";
	
	/**
	 * 
	 * @return true if checkSuccess.
	 */
	public boolean getCheckSuccess(){
		return checkSuccess;
	}
	
	/**
	 * 
	 * @return the errors if any.
	 */
	public String getCheckError(){
		return checkError;
	}
	
	/**
	 * Returns the cheksum of a file using MD5 Mesage Digest
	 * 
	 * @param file
	 *            the file
	 * @return its checksum
	 */
	public String checkSum(File file) {
		
		checkSuccess = true;
		checkError = "";
		try {
			InputStream fin = new FileInputStream(file);
			MessageDigest md5 = MessageDigest.getInstance("MD5");
			byte[] buffer = new byte[1024];
			int read;

			do {
				read = fin.read(buffer);
				if (read > 0)
					md5.update(buffer, 0, read);
			} while (read != -1);

			fin.close();

			byte[] digest = md5.digest();
			if (digest == null)
				return null;

			String strDigest = "";// "0x";
			for (int i = 0; i < digest.length; i++) {
				strDigest += Integer.toString((digest[i] & 0xff) + 0x100, 16).substring(1);// .toUpperCase();
				// the parameter is an int, so digest[i] is promoted to an int
				// by leading left-most bit which can be 0 or 1
				// digest[i] & 0xff is necesary for not having negatives (if
				// leading bit is 1)
				// 0xff is the int literal 00 00 00 ff -8 x 4 =32 bits. Hitting
				// with my int converted byte, which for negative
				// is something like ff ff ff fe (if my byte is fe) then it
				// becomes 00 00 00 fe as it should be.
				// then adds 0x100 which is 1 0000 0000 (binary) which is 256 in
				// base 10.
				// 0x100 is added to ensure that the hex string is 3 chars long;
				// (3 bit-paterns)
				// this is needed because we want a string with 2 hex digits for
				// each byte. A byte
				// is stored on 2 bit-pattern, and a byte with leading bit
				// pattern 0, i.e. between 0 and 15(f in hex) in last
				// bit-pattern, would produce 1 char only.
				// Finally the leftmost (third) digit (which is obvious 1) is
				// discarded with substring(1) since it does not
				// belong to the byte!
				// 16 is the base 16!
			}
			return strDigest;
		} catch (Exception e) {
			checkSuccess = false;
			checkError = e.toString();
			e.printStackTrace();
			return null;
		}
	}
}
