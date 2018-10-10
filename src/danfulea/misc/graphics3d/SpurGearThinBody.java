package danfulea.misc.graphics3d;

import javax.media.j3d.Appearance;

/**
 * Class for visualize a 3D SpurGearThinBody. It is based on old java3D tutorials freely available on the web (at that time).
 * 
 * @author Dan Fulea, 12 MAY. 2005
 *
 */

public class SpurGearThinBody extends SpurGear {

	/**
	 * Construct a SpurGearThinBody;
	 * 
	
	 * @param toothCount
	 *            number of teeth
	 * @param pitchCircleRadius
	 *            radius at center of teeth
	 * @param shaftRadius
	 *            radius of hole at center
	 * @param addendum
	 *            distance from pitch circle to top of teeth
	 * @param dedendum
	 *            distance from pitch circle to root of teeth
	 * @param gearThickness
	 *            thickness of the gear
	 */
	public SpurGearThinBody(int toothCount, float pitchCircleRadius,
			float shaftRadius, float addendum, float dedendum,
			float gearThickness) {
		this(toothCount, pitchCircleRadius, shaftRadius, addendum, dedendum,
				gearThickness, gearThickness, 0.25f, null);
	}

	/**
	 * Construct a SpurGearThinBody;
	 * 

	 * @param toothCount
	 *            number of teeth
	 * @param pitchCircleRadius
	 *            radius at center of teeth
	 * @param shaftRadius
	 *            radius of hole at center
	 * @param addendum
	 *            distance from pitch circle to top of teeth
	 * @param dedendum
	 *            distance from pitch circle to root of teeth
	 * @param gearThickness
	 *            thickness of the gear
	 * @param look
	 *            the gear's appearance
	 */
	public SpurGearThinBody(int toothCount, float pitchCircleRadius,
			float shaftRadius, float addendum, float dedendum,
			float gearThickness, Appearance look) {
		this(toothCount, pitchCircleRadius, shaftRadius, addendum, dedendum,
				gearThickness, gearThickness, 0.25f, look);
	}

	/**
	 * Construct a SpurGearThinBody;
	 * 

	 * @param toothCount
	 *            number of teeth
	 * @param pitchCircleRadius
	 *            radius at center of teeth
	 * @param shaftRadius
	 *            radius of hole at center
	 * @param addendum
	 *            distance from pitch circle to top of teeth
	 * @param dedendum
	 *            distance from pitch circle to root of teeth
	 * @param gearThickness
	 *            thickness of the gear
	 * @param toothTipThickness
	 *            thickness of the tip of the tooth
	 * @param look
	 *            the gear's appearance
	 */
	public SpurGearThinBody(int toothCount, float pitchCircleRadius,
			float shaftRadius, float addendum, float dedendum,
			float gearThickness, float toothTipThickness, Appearance look) {
		this(toothCount, pitchCircleRadius, shaftRadius, addendum, dedendum,
				gearThickness, toothTipThickness, 0.25f, look);
	}

	/**
	 * Construct a SpurGearThinBody;
	 * 

	 * @param toothCount
	 *            number of teeth
	 * @param pitchCircleRadius
	 *            radius at center of teeth
	 * @param shaftRadius
	 *            radius of hole at center
	 * @param addendum
	 *            distance from pitch circle to top of teeth
	 * @param dedendum
	 *            distance from pitch circle to root of teeth
	 * @param gearThickness
	 *            thickness of the gear
	 * @param toothTipThickness
	 *            thickness of the tip of the tooth
	 * @param toothToValleyAngleRatio
	 *            ratio of tooth valley to circular pitch (must be less or equal to .25)
	 * @param look
	 *            the gear's appearance object
	 */
	public SpurGearThinBody(int toothCount, float pitchCircleRadius,
			float shaftRadius, float addendum, float dedendum,
			float gearThickness, float toothTipThickness,
			float toothToValleyAngleRatio, Appearance look) {

		this(toothCount, pitchCircleRadius, shaftRadius, addendum, dedendum,
				gearThickness, toothTipThickness, 0.25f, look,
				0.6f * gearThickness, 0.75f * (pitchCircleRadius - shaftRadius));
	}

	/**
	 * Construct a SpurGearThinBody;
	 * 
	 * @param toothCount
	 *            number of teeth
	 * @param pitchCircleRadius
	 *            radius at center of teeth
	 * @param shaftRadius
	 *            radius of hole at center
	 * @param addendum
	 *            distance from pitch circle to top of teeth
	 * @param dedendum
	 *            distance from pitch circle to root of teeth
	 * @param gearThickness
	 *            thickness of the gear
	 * @param toothTipThickness
	 *            thickness of the tip of the tooth
	 * @param toothToValleyAngleRatio
	 *            ratio of tooth valley to circular pitch (must be less or equal to .25)
	 * @param look
	 *            the gear's appearance object
	 * @param bodyThickness
	 *            the thickness of the gear body
	 * @param crossSectionWidth
	 *            the width of the depressed portion of the gear's body
	 */
	public SpurGearThinBody(int toothCount, float pitchCircleRadius,
			float shaftRadius, float addendum, float dedendum,
			float gearThickness, float toothTipThickness,
			float toothToValleyAngleRatio, Appearance look,
			float bodyThickness, float crossSectionWidth) {

		super(toothCount, pitchCircleRadius, addendum, dedendum,
				toothToValleyAngleRatio);

		float diskCrossSectionWidth = (rootRadius - shaftRadius - crossSectionWidth) / 2.0f;
		float outerShaftRadius = shaftRadius + diskCrossSectionWidth;
		float innerToothRadius = rootRadius - diskCrossSectionWidth;

		// Generate the gear's body disks, first by the shaft, then in
		// the body and, lastly, by the teeth
		addBodyDisks(shaftRadius, outerShaftRadius, gearThickness, look);
		addBodyDisks(innerToothRadius, rootRadius, gearThickness, look);
		addBodyDisks(outerShaftRadius, innerToothRadius, bodyThickness, look);

		// Generate the gear's "shaft" equivalents the two at the teeth
		// and the two at the shaft
		addCylinderSkins(innerToothRadius, gearThickness, InwardNormals, look);
		addCylinderSkins(outerShaftRadius, gearThickness, OutwardNormals, look);

		// Generate the gear's interior shaft
		addCylinderSkins(shaftRadius, gearThickness, InwardNormals, look);

		// Generate the gear's teeth
		addTeeth(pitchCircleRadius, rootRadius, outsideRadius, gearThickness,
				toothTipThickness, toothToValleyAngleRatio, look);
	}

}
