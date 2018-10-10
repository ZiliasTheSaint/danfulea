package danfulea.misc.graphics3d;

import javax.media.j3d.Appearance;
import javax.media.j3d.Shape3D;
import javax.media.j3d.TriangleStripArray;
import javax.vecmath.Color4f;
import javax.vecmath.Point3d;

/**
 * Class for visualize a 3D Moebius. It is based on old java3D tutorials freely available on the web (at that time).
 * 
 * @author Dan Fulea, 12 MAY. 2005
 *
 */

public class Moebius extends Shape3D {
	// stripLength is the number of points plotted in each segment of the strip.
	// The number of points stored in each strip is twice stripLength since each
	// strip stores points that zig-zag back and forth between two segments.
	private int stripLength;//, strips;
	private TriangleStripArray tsa;

	public Moebius(int strips, int stripLength, Color4f c1, Color4f c2,// )
			Appearance appearance)// )
	{
		// -----------THIS PART IS COMMON--------------------------------------
		tsa = new TriangleStripArray(2 * stripLength * strips,
				TriangleStripArray.COORDINATES | TriangleStripArray.COLOR_4,
				indexArray(strips, stripLength));
		;
		this.stripLength = stripLength;
		//this.strips = strips;
		for (int strip = 0; strip < strips; strip++) {
			for (int stripIndex = 0; stripIndex < stripLength; stripIndex++) {
				setCoordinate(strip, stripIndex, 0, c1, c2);
				setCoordinate(strip, stripIndex, 1, c1, c2);
			}
		}
		// set the geometry and appearance
		this.setGeometry(tsa);
		this.setAppearance(appearance);
		// -----------END COMMON--------------------------------------

	}

	// -----------THIS PART IS COMMON--------------------------------------
	private int[] indexArray(int strips, int stripLength) {
		int[] indexArray = new int[strips];
		for (int i = 0; i < strips; i++)
			indexArray[i] = 2 * stripLength;
		return indexArray;
	}

	// Divide by stripLength-1 to include both 0 and 360.
	private double radians(int index) {
		return index * 2 * Math.PI / (stripLength - 1);
	}

	private void setCoordinate(int strip, int stripIndex, int segment,
			Color4f c1, Color4f c2) {
		double x = x(strip, stripIndex, segment);
		double y = y(strip, stripIndex, segment);
		double z = z(strip, stripIndex, segment);
		Point3d p = new Point3d(x, y, z);
		tsa.setCoordinate(strip * 2 * stripLength + 2 * stripIndex + segment, p);
		boolean even = strip % 2 == 0;
		Color4f c = even ? (segment == 0 ? c1 : c2) : (segment == 0 ? c2 : c1);
		tsa.setColor(strip * 2 * stripLength + 2 * stripIndex + segment, c);
	}

	// -----------END COMMON--------------------------------------
	private double x(int strip, int index, int segment) {
		double a = radians(index * 2);
		return 0.7 * Math.sin(a) + 0.2 * Math.cos(a) * (segment * 2 - 1);
	}

	private double y(int strip, int index, int segment) {
		double a = radians(index * 2);
		return 0.3 * Math.sin(a) * (segment * 2 - 1);
	}

	private double z(int strip, int index, int segment) {
		double a = radians(index * 2);
		return 0.7 * Math.cos(a) + 0.2 * Math.cos(a) * (segment * 2 - 1);
	}

}