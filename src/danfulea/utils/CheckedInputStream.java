package danfulea.utils;

import java.io.FilterInputStream;
import java.io.InputStream;
import java.io.IOException;

import java.util.zip.Checksum;

/**
 * Utility class for streaming CRC32 checsum thanks to java sun tutorial;
 * 
 * @author Sun developer, modified by Dan Fulea 23 FEB.2015
 * Now, it is part of jdk...java.util.zip package
 */
@Deprecated
public class CheckedInputStream extends FilterInputStream {
	private Checksum cksum;

	/**
	 * Constructor
	 * 
	 * @param in
	 *            the input stream
	 * @param cksum
	 *            the checksum object
	 */
	public CheckedInputStream(InputStream in, Checksum cksum) {
		super(in);
		this.cksum = cksum;
	}

	/**
	 * Read and update the checksum
	 */
	public int read() throws IOException {
		int b = in.read();
		if (b != -1) {
			cksum.update(b);//
		}
		return b;
	}

	/**
	 * Read and update the checksum
	 */
	public int read(byte[] b) throws IOException {
		int len;
		len = in.read(b, 0, b.length);
		if (len != -1) {
			cksum.update(b, 0, len);
		}
		return len;
	}

	/**
	 * Read and update the checksum
	 */
	public int read(byte[] b, int off, int len) throws IOException {
		len = in.read(b, off, len);
		if (len != -1) {
			cksum.update(b, off, len);
		}
		return len;
	}

	/**
	 * 
	 * @return the checksum
	 */
	public Checksum getChecksum() {
		return cksum;
	}
}