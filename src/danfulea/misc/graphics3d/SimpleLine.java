package danfulea.misc.graphics3d;

import javax.media.j3d.Shape3D;
import javax.media.j3d.LineArray;
import javax.vecmath.Color4f;
import javax.vecmath.Point3d;

//import javax.vecmath.Point3f;

/**
 * Class for visualize a 3D SimpleLine. It is based on old java3D tutorials freely available on the web (at that time).
 * 
 * @author Dan Fulea, 12 MAY. 2005
 *
 */

public class SimpleLine extends Shape3D {
	public SimpleLine(String name, Point3d start, Point3d end, Color4f color1,
			Color4f color2)

	{
		LineArray la = new LineArray(2, LineArray.COORDINATES
				| LineArray.COLOR_4);
		la.setCoordinate(0, start);
		la.setCoordinate(1, end);
		la.setColor(0, color1);
		la.setColor(1, color2);
		this.setGeometry(la);
		this.setUserData(name);
		this.setCapability(Shape3D.ALLOW_APPEARANCE_READ);
		this.setCapability(Shape3D.ALLOW_GEOMETRY_READ);
	}

}