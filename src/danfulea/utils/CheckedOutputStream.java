package danfulea.utils;

import java.io.FilterOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.util.zip.Checksum;

/**
 * Utility class for streaming CRC32 checsum thanks to java sun tutorial;
 * 
 * @author Sun developer, modified by Dan Fulea 23 FEB.2015
 * Now, it is part of jdk...java.util.zip package
 */
@Deprecated
public class CheckedOutputStream extends FilterOutputStream {
	private Checksum cksum;

	/**
	 * Constructor
	 * 
	 * @param out
	 *            the Output Stream
	 * @param cksum
	 *            the checksum object
	 */
	public CheckedOutputStream(OutputStream out, Checksum cksum) {
		super(out);
		this.cksum = cksum;//
	}

	/**
	 * Writes the stream and updates the checksum
	 */
	public void write(int b) throws IOException {
		out.write(b);
		cksum.update(b);
	}

	/**
	 * Writes the stream and updates the checksum
	 */
	public void write(byte[] b) throws IOException {
		out.write(b, 0, b.length);
		cksum.update(b, 0, b.length);
	}

	/**
	 * Writes the stream and updates the checksum
	 */
	public void write(byte[] b, int off, int len) throws IOException {
		out.write(b, off, len);
		cksum.update(b, off, len);
	}

	/**
	 * 
	 * @return the checksum
	 */
	public Checksum getChecksum() {
		return cksum;
	}
}