package danfulea.misc.graphics3d;

import javax.media.j3d.Appearance;
import javax.media.j3d.Shape3D;
import javax.media.j3d.TriangleStripArray;
import javax.vecmath.Color4f;
import javax.vecmath.Point3d;

/**
 * Class for visualize a 3D RightEllipsoid. It is based on old java3D tutorials freely available on the web (at that time).
 * 
 * @author Dan Fulea, 12 MAY. 2005
 *
 */

public class RightEllipsoid extends Shape3D {
	// stripLength is the number of points plotted in each segment of the strip.
	// The number of points stored in each strip is twice stripLength since each
	// strip stores points that zig-zag back and forth between two segments.
	private int stripLength, strips;
	private TriangleStripArray tsa;
	// particular constants
	private double a = 0.0;// great semiaxe x
	private double b = 0.0;// great semiaxe y
	private double c = 0.0;// great semiaxe z

	// general equiation of ellipsoid: x^2/a^2+y^2/b^2+z^2/c^2=1.0;

	public RightEllipsoid(int strips, int stripLength, Color4f c1, Color4f c2,// )
			Appearance appearance, double a, double b, double c)

	{
		this.a = a;
		this.b = b;
		this.c = c;
		// -----------THIS PART IS COMMON--------------------------------------
		tsa = new TriangleStripArray(2 * stripLength * strips,
				TriangleStripArray.COORDINATES | TriangleStripArray.COLOR_4,
				indexArray(strips, stripLength));
		;
		this.stripLength = stripLength;
		this.strips = strips;
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
	//private double distance(int strip, int index, int segment) {
	//	double x = x(strip, index, segment);
	//	double z = z(strip, index, segment);
	//	return Math.sqrt(x * x + z * z);
	//}

	private double radiusx(int strip, int segment) {
		
		// return (strip + 1 - segment)/(double)strips;
		// ----------------------------------------------------------------------
		double radius = 0.0;// a=0.4/2!!!
		if (strip >= (double) strips / 2)
			radius = a;
		else
			// radius=2*a*Math.abs((strip + 1 -
			// segment-(double)strips/2))/(double)strips;
			radius = 2 * a * Math.abs((strip + 1 - segment)) / (double) strips;
		return radius;
		// ----------------------------------------------------------------------
		
		// return (strip + segment)/(double)strips;
	}

	private double radiusz(int strip, int segment) {
		
		// return (strip + 1 - segment)/(double)strips;
		// ----------------------------------------------------------------------
		double radius = 0.0;// c=0.2/2!!!
		if (strip >= (double) strips / 2)
			radius = c;
		else
			// radius=2*c*Math.abs((strip + 1 -
			// segment-(double)strips/2))/(double)strips;
			radius = 2 * c * Math.abs((strip + 1 - segment)) / (double) strips;
		return radius;
		// ----------------------------------------------------------------------
		
		// return (strip + segment)/(double)strips;
	}

	private double x(int strip, int index, int segment) {
		return radiusx(strip, segment) * Math.cos(radians(index));
	}

	private double y(int strip, int index, int segment) {
		double yy = 0.0;
		if (strip >= (double) strips / 2)
			yy = 2 * b * (strip + 1 - segment - (double) strips / 2)
					/ (double) strips;
		else
			yy = 0.0;
		return yy;
	}

	private double z(int strip, int index, int segment) {
		return radiusz(strip, segment) * Math.sin(radians(index));
	}

}