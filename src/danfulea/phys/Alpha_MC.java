package danfulea.phys;

import danfulea.math.RandomCollection;

/**
 * The main class for Monte-Carlo calculation of alpha detector efficiency. <br>
 * Results are accurate ONLY for vacuum based alpha detector systems!! 
 * 
 * 
 * @author Dan Fulea, 05 APR. 2011
 * 
 */
public class Alpha_MC {

	private double DISTZ = 1.0;// distance from disk source to disk detector
	private double adet = 1.0;// detector radius
	private double asource = 1.0;// source radius
	private double WEIGHT = 1.0;// initial particle weight
	private int NCASE = 100;// number of histories
	private double nscore = 0.0;
	private double nscore2 = 0.0;
	private double solidAngle = 0.0;
	private double eff = 0.0;
	private double effp = 0.0;
	private double erreff = 0.0;// 2 sigma
	public static boolean STOPSIMULATION = false;

	/**
	 * Constructor.
	 * 
	 * @param DISTZ
	 *            source to detector distance in cm.
	 * @param asource
	 *            source diameter in cm
	 * @param adet
	 *            detector diameter in cm
	 * @param NCASE
	 *            number of simulations (histories)
	 */
	public Alpha_MC(double DISTZ, double asource, double adet, int NCASE) {
		RandomCollection.reset();
		STOPSIMULATION = false;// mark as run!
		this.setDISTZ(DISTZ);// DISTZ = DISTZ;// cm
		this.setAdet(adet / 2.0);// adet = adet / 2.0;// cm
		this.setAsource(asource / 2.0);// asource = asource / 2.0;// cm
		this.setNCASE(NCASE);
		nscore = 0.0;// initialize
		nscore2 = 0.0;

		for (int i = 1; i <= NCASE; i++) {
			if (STOPSIMULATION) {
				return;// force exit!
			}
			individualHistory();
		}

		nscore = nscore / NCASE;
		// stat
		nscore2 = nscore2 / NCASE;
		nscore2 = (nscore2 - nscore * nscore) / (NCASE - 1.0);
		if (nscore2 >= 0.)
			nscore2 = Math.sqrt(nscore2);
		if (nscore != 0.) {
			nscore2 = Math.min(nscore2 / nscore * 100., 99.9);
		} else {
			nscore2 = 99.9;
		}// procentual!

		// finaly:
		setEff(nscore);
		setEffp(eff * 100.0);
		setErreff(2.0 * nscore2);// 2 sigma
		setSolidAngle(eff * 4.0 * Math.PI);

		// System.out.println("alpha eff= " + effp + " %"
		// + " ;uncertainty [2sigma] in %: " + erreff);
		// System.out.println("solid angle= " + eff * 4.0 * Math.PI + " sr");

	}

	/**
	 * 
	 * @return the source to detector distance.
	 */
	public double getDISTZ() {
		return DISTZ;
	}

	/**
	 * 
	 * Set source to detector distance.
	 * @param dISTZ dISTZ
	 */
	public void setDISTZ(double dISTZ) {
		DISTZ = dISTZ;
	}

	/**
	 * 
	 *@return the detector radius.
	 */
	public double getAdet() {
		return adet;
	}

	/**
	 * 
	 * Set detector radius.
	 * @param adet adet
	 */
	public void setAdet(double adet) {
		this.adet = adet;
	}

	/**
	 * 
	 * @return the source radius.
	 */
	public double getAsource() {
		return asource;
	}

	/**
	 * 
	 * Set source radius.
	 * @param asource asource
	 */
	public void setAsource(double asource) {
		this.asource = asource;
	}

	/**
	 * 
	 * Set number of histories.
	 * @param nCASE nCASE
	 */
	public void setNCASE(int nCASE) {
		NCASE = nCASE;
	}

	/**
	 * 
	 * @return the number of histories.
	 */
	public int getNCASE() {
		return NCASE;
	}

	/**
	 * 
	 * Set solid angle. Should never be used. This is computed internally!
	 * @param solidAngle solidAngle
	 */
	private void setSolidAngle(double solidAngle) {
		this.solidAngle = solidAngle;
	}

	/**
	 * 
	 * @return the solid angle.
	 */
	public double getSolidAngle() {
		return solidAngle;
	}

	/**
	 * 
	 * @return the alpha efficiency [number less than 1.0].
	 */
	public double getEff() {
		return eff;
	}

	/**
	 * 
	 * Set alpha efficiency [number less than 1.0]. Should never be used. This is computed
	 * internally!
	 * @param eff eff
	 */
	private void setEff(double eff) {
		this.eff = eff;
	}

	/**
	 * 
	 * @return the alpha efficiency [%].
	 */
	public double getEffp() {
		return effp;
	}

	/**
	 * 
	 * Set alpha efficiency [%]. Should never be used. This is computed
	 * internally!
	 * 
	 * @param effp effp
	 */
	private void setEffp(double effp) {
		this.effp = effp;
	}

	/**
	 * 
	 * @return the alpha efficiency uncertainty [%].
	 */
	public double getErreff() {
		return erreff;
	}

	/**
	 * 
	 * Set alpha efficiency uncertainty [%]. Should never be used. This is
	 * computed internally!
	 * @param erreff erreff
	 */
	private void setErreff(double erreff) {
		this.erreff = erreff;
	}

	/**
	 * Computing a single, random history.
	 */
	private void individualHistory() {
		// evaluation of initial coordinates inside the source disk
		double r = RandomCollection.random01();
		double ro = asource * Math.sqrt(r);// dist to z axis; ro evaluation
		r = RandomCollection.random01();
		double phi1 = 2 * Math.PI * r;
		double x = ro * Math.cos(phi1);
		double y = ro * Math.sin(phi1);
		// double z = -DISTZ;// for consistency...not used!

		// evaluation of weight for variance reduction=>we asure that the solid
		// angle from the
		// radiation emission point encompass the detector!!
		double diam = 0.0;
		if (asource <= adet) {
			diam = 2.0 * adet;
		} else {
			diam = 2.0 * asource;
		}
		// cosines of theta max.
		double costmax = DISTZ / Math.sqrt(DISTZ * DISTZ + diam * diam);
		// maximum solid angle
		double us = 2.0 * Math.PI * (1.0 - costmax);
		WEIGHT = us / (4 * Math.PI);
		// this is the actual particle weight according to the above variance
		// reduction!!

		// actual costet evaluation (for polar angle evaluation) and the
		// traversed distance
		double dom = (1 - costmax) / 2;// >0 and <1/2, costmax <1 and >0,
										// thetamax<90!
		r = RandomCollection.random01();
		r = r * dom;
		double costet = 1 - 2 * r;// >0 always--- positive z axis!!
		double sintet = Math.sqrt(1 - costet * costet);
		double traversed_distance = DISTZ / costet;

		// actual azimutal angle evaluation
		r = RandomCollection.random01();
		double phi2 = 2 * Math.PI * r;// azimutal angle
		// @REDUNDANT
		// double tgtet=sintet/costet;.....REDUNDANT_CHECKED
		// double teta=Math.abs(Math.atan(tgtet));//-pi/2,pi/2
		// initial we have formal directional cosines u,v,w = 0,0,-1 and then we
		// have teta and phi2:
		// double u=Math.sin(teta)*Math.cos(phi2);//@@@@@@@@@@@@@@@@@@@@@@@@@@@
		// double v=Math.sin(teta)*Math.sin(phi2);
		// @REDUNDANT!!

		// initial directional cosines
		double u = sintet * Math.cos(phi2);
		double v = sintet * Math.sin(phi2);
		// double w = costet;// >0...for consistency, not used!

		// transport particle

		x = x + traversed_distance * u;
		y = y + traversed_distance * v;
		// z = 0.0;// for consistency, not used!!

		// test if hit the detector, if hit, the particle is scored!!
		if (x * x + y * y <= adet * adet) {
			// score
			nscore = nscore + WEIGHT;
			nscore2 = nscore2 + WEIGHT * WEIGHT;
		}

	}

}
