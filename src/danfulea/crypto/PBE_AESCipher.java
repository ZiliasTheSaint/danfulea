package danfulea.crypto;

import java.io.InputStream;
import java.io.OutputStream;
import java.security.Key;
import java.security.MessageDigest;
import java.security.SecureRandom;

import javax.crypto.Cipher;
import javax.crypto.CipherInputStream;
import javax.crypto.CipherOutputStream;
import javax.crypto.spec.IvParameterSpec;
import javax.crypto.spec.SecretKeySpec;

/**
 * PBE_AESCipher class for handling the password based encryption.<br>
 * <br>
 * User password is transformed into a valid key for the encryption/decryption
 * algorithm. The application encryption/decryption algorithm is
 * "AES/CBC/PKCS5Padding" and it is considered the best available! It uses a
 * default 256 bits key, which is strong enough for any purpose. This 256 bits
 * key is generated using a message digest SHA-256 algorithm and it is based on
 * user input password! <br>
 * This application uses the advanced SecureRandom class for random number
 * generator. <br>
 * Moreover, (which can be considered redundant), it is used a 8 bytes salt and
 * it is performed 1000 iterations in order to derive the "fingerprints" of
 * initial user password using the message digest algorithm. However, initial
 * "fingerprints" (without using the salt and iterations) should be enough
 * because:<br>
 * <br>
 * THE ENCRYPTION KEY IS NEVER STORED BY THIS APPLICATIONS! IT EXISTS ONLY IN
 * THE MIND OF THE USER AND IT IS DERIVED AT RUNTIME FROM THE USER PASSWORD!
 * THEREFORE ANY ATTEMPT OF REVERSE ENGINEERING THE ENCRYPTION KEY WILL FAIL!
 * <br>
 * <br>
 * Message diggest algorithm MD5 generates a fixed 128 bits hash value!<br>
 * - It can be use for standard 128 bits AES key long!<br>
 * Message diggest algorithm SHA or SHA1 generates a fixed 160 bits hash value!
 * <br>
 * - No use, AES key must be of 128, 192 or 256 bits long!<br>
 * Message diggest algorithm SHA-256 generates a fixed 256 bits hash value!<br>
 * - It can be use for standard 256 bits AES key long!<br>
 * MD5 is considered weak but it is not important since the encryption key is
 * never stored or transmitted! Using MD5 affects only the strength of AES
 * encryption which will be based on 128 bits key instead of 256 bits key.
 * Although 256 bits AES encryption is the best, both encryptions are strong
 * enough for any purpose! <br>
 * The JRE provided by Sun Java(tm) does not work by default with 256 bits AES
 * encryption key, so if user plans to use the default JRE, he can only use the
 * MD5 message digest algorithms which generates an 128 bits AES encryption key.
 * <br>
 * In order to run the SHA-256 hash algorithm, which generates a 256 bits AES
 * encryption key, the user must download unlimited encryption policy jars and
 * override the default jars in JAVA_HOME folder as follows:<br>
 * - go to
 * <a href="http://java.sun.com/javase/downloads/index.jsp" target="_blank">
 * http://java.sun.com/javase/downloads/index.jsp</a> <br>
 * - scroll down and choose "download" for "Java Cryptography Extension (JCE)
 * Unlimited Strength Jurisdiction Policy "<br>
 * - Download "jce_policy-6.zip" (or later) and save it! Extract and copy the
 * two jars (local_policy.jar and US_export_policy.jar) on
 * JAVA_HOME/lib/security (e.g."C:\Program Files\Java\jre6\lib\security")
 * overwriting the existed jars!<br>
 * Now you have "unlimited strength cryptography" and you can use 256 bits AES
 * key!!<br>
 * <br>
 * A weak password can easily be broken!<br>
 * Best scenario: use a strong password (10 characters minimum) and 256 bits AES
 * encryption key (or 128 bits if standard JRE is installed) which means using
 * SHA-256 algorithm (or MD5 algorithm if standard JRE is installed) for message
 * digest.<br>
 * <br>
 * <br>
 * 
 * @author Dan Fulea, 28 MAR 2010
 */

public class PBE_AESCipher {
	
	//public constants
	public static final String MESSAGE_DIGGEST_ALGORITHM_SHA256 = "SHA-256";
	public static final String MESSAGE_DIGGEST_ALGORITHM_MD5 = "MD5";
	
	//private constants; If it cannot change, there is no point having one copy per 
	//instance so that's why they are also static.
	private static final String ENCRYPTION_ALGORITHM = "AES/CBC/PKCS5Padding";
	private static final String ENCRYPTION_KEY_ALGORITHM = "AES";	
	private static final int ENCRYPTION_IV_LENGTH = 16;
	private static final int SALT_LENGTH = 8;
	private static final int ITERATIONS = 1000;
	
	//instance variables
	public String MESSAGE_DIGGEST_ALGORITHM = "SHA-256";// default	
	private byte[] encryptionIV;	
	private byte[] salt;	
	private SecureRandom random;
	private String password;
	private Key key;

	/**
	 * Constructor, a password input is required!
	 * 
	 * @param password
	 *            the password
	 * @throws Exception
	 *             can throw this exception
	 */
	public PBE_AESCipher(String password) throws Exception {
		this.password = password;
	}
	
	/**
	 * Sets the message digest algorithm
	 * @param MESSAGE_DIGGEST_ALGORITHM
	 * the message digest algorithm
	 */
	public void setMessageDigestAlgorithm(String MESSAGE_DIGGEST_ALGORITHM){
		this.MESSAGE_DIGGEST_ALGORITHM=MESSAGE_DIGGEST_ALGORITHM;
	}

	/**
	 * Generates random bytes of fixed length.
	 * 
	 * @param length
	 *            of byte array
	 * @return the byte array
	 */
	private byte[] randomBytes(int length) {
		if (random == null) {
			random = new SecureRandom();
		}

		byte[] bytes = new byte[length];
		random.nextBytes(bytes);

		return bytes;
	}

	/**
	 * Reads bytes from an input stream.
	 * 
	 * @param length
	 *            of bytes array
	 * @param in
	 *            the input stream
	 * @return the bytes array
	 * @throws Exception
	 *             can throw this exception
	 */
	private byte[] readBytes(int length, InputStream in) throws Exception {
		byte[] bytes = new byte[length];
		int read = in.read(bytes);

		if (read != length) {
			throw new Exception("expected length != actual length");
		}

		return bytes;
	}

	/**
	 * Gets the secure fingerprints from user password.
	 * 
	 * @param password
	 *            a user password
	 * @param salt
	 *            the salt used, see class description
	 * @return the bytes array (fingerprints)
	 * @throws Exception
	 *             can throw this exception
	 */
	private byte[] getHash(String password, byte[] salt) throws Exception {
		byte[] input = null;

		MessageDigest digest = MessageDigest.getInstance(MESSAGE_DIGGEST_ALGORITHM);
		digest.reset();
		digest.update(salt);

		input = digest.digest(password.getBytes("UTF-8"));
		for (int i = 0; i < ITERATIONS; i++) {
			digest.reset();
			input = digest.digest(input);
		}

		// now we have a secure fingerprints of user password
		return input;
	}

	/**
	 * Encrypts the output stream.
	 * 
	 * @param outputStream
	 *            user output stream.
	 * @return the encrypted output stream.
	 * @throws Exception
	 *             can throw this exception
	 */
	public OutputStream encrypt(OutputStream outputStream) throws Exception {
		if (key == null) {
			// create the salt for the password hash algorithm
			salt = randomBytes(SALT_LENGTH); // default 8 byte salt

			// generate the IV for AES encryption
			encryptionIV = randomBytes(ENCRYPTION_IV_LENGTH);// AES blocksize =
																// 16 bytes!

			// turn the password into an AES key using its hash (fingerprints)!!
			byte[] keyBytes = getHash(password, salt);// a valid 256 bits long!
			key = new SecretKeySpec(keyBytes, ENCRYPTION_KEY_ALGORITHM);
		}

		// write the header to the output stream. this has the format
		// (without "|" delimeters):
		// PBE IV|ENCRYPTION IV|
		outputStream.write(salt);
		outputStream.write(encryptionIV);

		// now create the encryption cipher
		Cipher cipher = Cipher.getInstance(ENCRYPTION_ALGORITHM);
		cipher.init(Cipher.ENCRYPT_MODE, key, new IvParameterSpec(encryptionIV));

		return new CipherOutputStream(outputStream, cipher);
	}

	/**
	 * Decrypts the input stream.
	 * 
	 * @param in
	 *            the encrypted input stream
	 * @return the decrypted input stream
	 * @throws Exception
	 *             can throw this exception
	 */
	public InputStream decrypt(InputStream in) throws Exception {
		// Read the header of the encrypted file.
		byte[] salt = readBytes(SALT_LENGTH, in);
		byte[] encryptionIV = readBytes(ENCRYPTION_IV_LENGTH, in);

		// turn the password into an AES key using its hash (fingerprints)!!
		byte[] keyBytes = getHash(password, salt);
		Key key = new SecretKeySpec(keyBytes, ENCRYPTION_KEY_ALGORITHM);

		// If we haven't yet generated a key, just use this one
		if (this.key == null) {
			this.salt = salt;
			this.encryptionIV = encryptionIV;
			this.key = key;
		}

		// now create the decrypt cipher
		Cipher cipher = Cipher.getInstance(ENCRYPTION_ALGORITHM);
		cipher.init(Cipher.DECRYPT_MODE, key, new IvParameterSpec(encryptionIV));

		return new CipherInputStream(in, cipher);
	}

	// ---------------------------------Testing--------------------------
	/*
	public static void main(String args[]) {
		String myPass = "fhgisryu";// "qwertysrjhm,adnga,gbadf@344;alfjf54b))(9dswefhwer";
		try {

			PBE_AESCipher aesc = new PBE_AESCipher(myPass);

			File outfile = new File("test_encrypt"); ///
			FileOutputStream fos = new FileOutputStream(outfile);
			OutputStream os = aesc.encrypt(fos); //

			String ss = "I trust in Java(tm) Sun!";
			byte[] bs = ss.getBytes("UTF-8"); // now use os to write the data:
			os.write(bs);

			os.close();
			fos.close();
			
			//--------------------------
			aesc = new PBE_AESCipher(myPass);
			InputStream fin1 = new FileInputStream(outfile);
			InputStream fin = aesc.decrypt(fin1);

			StringBuffer sbuf = new StringBuffer();
			File f = new File("test_decrypted");
			OutputStream out = new FileOutputStream(f);

			byte buf[] = new byte[1024];
			int len;
			while ((len = fin.read(buf)) > 0) {
				sbuf.append(new String(buf, 0, len));
				out.write(buf, 0, len);
			}
			String sdec = sbuf.toString();
			System.out.println(sdec);

			out.close();
			fin.close();
			fin1.close(); // Works
		} catch (Exception e) {
			e.printStackTrace();
		}

	} // end main
	*/
}// end class