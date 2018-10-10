package danfulea.phys;

import java.util.ResourceBundle;

import danfulea.math.RandomCollection;
import danfulea.math.Sort;
import danfulea.math.StatsUtil;
import danfulea.math.numerical.Stats;

/**
 * Compute gamma detector efficiency using KERMA approximation. Since the geometry is cylindrical, 
 * i.e. RZ geometry, a much better approach is using EGSnrc (or GEANT4) toolkit. Therefore, 
 * most methods of this class are deprecated. 
 * @author Dan Fulea, 2005
 *
 */

public class GammaDetector {
	private static final String BASE_RESOURCE_CLASS = "danfulea.phys.resources.GammaDetectorResources";
	private static ResourceBundle resources = ResourceBundle
			.getBundle(BASE_RESOURCE_CLASS);
	private static String simulationTimeElapsed = "";// store simulation time
														// elapsed
	private static int n_total = 0;// stores the total number of simulated
									// photons for each eff eval
	private static int n_pair = 0;// store the total number of pair production
	private static int n_total_all = 0;// stores the total number of simulated
										// photons
	private static int n_photoefect = 0;// stores the total number of
										// photoelectric absorbed photons
	private static int n_fwhm = 0;// stores the total number of photons from
									// resolution consideration!
	private static double[] efficiency;// efficiency NOT in % plus its error!!!
	private static double pondt = 0.0;// total weight of photon
	private static double pond = 0.0;// reprezents the geometric efficiency.
										// This will be added to
	// total weight only if the photon loses its entire energy (fwhm discussion)
	private static double source_parcurs = 0.0;
	private static int isValidInteraction = 0;// false, 1=true
	private static double adet = 0.0;// detector radius
	private static double asource = 0.0;// source radius
	private static double hdet = 0.0;// detector height
	private static double hsource = 0.0;// source height
	private static double hsourceup = 0.0;// upper source height---MARINELLI
	private static double bsource = 0.0;// inner radius of beaker---MARINELLI
	private static boolean isM = false;

	private static double[] effind;

	private static boolean stopB = false;
	private static boolean geomB = false;// &&&&&&&&&&&
	private static double[] wrho;// &&&&&&&&
	private static double[] watt;// &&&&&&&&&
	// -------------------------------------------------------
	private static int knAlgor = 0;// Wielopolski

	/**
	 * Set Klein Nishina sampling algorithm
	 * @param kn kn
	 */
	@Deprecated
	public static void setKNAlgor(int kn) {
		knAlgor = kn;
	}

	// -------------------------------------------------------
	@Deprecated
	public static void setStop() {
		stopB = true;
	}

	@Deprecated
	public static void setStart() {
		stopB = false;
	}

	@Deprecated
	public static double[] getEfficiency() {
		return efficiency;
	}

	@Deprecated
	public static double[] getEfficiencySet() {
		return effind;
	}

	@Deprecated
	public static int getTotalPhotons() {
		return n_total_all;
	}

	@Deprecated
	public static int getResPhotons() {
		return n_fwhm;
	}

	@Deprecated
	public static int getPhotoAbsorbedPhotons() {
		return n_photoefect;
	}

	@Deprecated
	public static int getPairProduction() {
		return n_pair;
	}

	@Deprecated
	public static String getSimulationTimeElapsed() {
		return simulationTimeElapsed;
	}

	// ---------------------------------------------------------------------------------------------
	@Deprecated
	public static double comptonSim(double k) {
		double theta = 0.0;
		if (knAlgor == 0)
			theta = comptonSimW(k);// Wielopolski algor--the best!
		else if (knAlgor == 1)
			theta = comptonSimK(k);// Kahn algor--not quite good
		else if (knAlgor == 2)
			theta = comptonSimC(k);// Clasic algor--not quite good
		else
			theta = comptonSimE(k);// EGS4+EGSnrc code
		return theta;
	}

	/**
	 * Return polar angle using Klein-Nishina sampling based on EGS4 algorithm.  
	 * @param k energy in MeV
	 * @return the result
	 */
	public static double comptonSimE(double k)// EGS4
	{
		// int n_acc=0; //number of accepted Compton events
		double theta = 0.0;// the output value!!
		double m = 0.511;// rest mass of electron in MeV->mc2
		double ko = k / m;// Gamma energy in units of electron rest energy
		double broi = 1 + 2 * ko;// Needed for scattering angle sampling
		double br = 0.0;
		double temp = 0.0;
		double sinthe = 0.0;
		double costhe = 0.0;
		double rejf3 = 0.0;// rejection function
		// store the 3 random generated numbers used in simulation
		double r1 = 0.0;
		double r2 = 0.0;
		double r3 = 0.0;

		while (true)// added
		{// added
		// :RESAMPLE:
			if (ko > 2) // At high energies the original EGS4 method is most
						// efficient
			{
				double broi2 = broi * broi;
				double alph1 = Math.log(broi);
				double alph2 = ko * (broi + 1) / broi2;
				double alpha = alph1 / (alph1 + alph2);

				while (true)// n_acc==0)
				{

					r1 = RandomCollection.random01();
					r2 = RandomCollection.random01();
					if (r1 < alpha) { // "Use 1/br part"
						br = Math.exp(alph1 * r2) / broi;
					} else { // "Use the br part."
						br = Math.sqrt(r2 + (1 - r2) / broi2);
					}
					temp = (1 - br) / (ko * br);
					sinthe = Math.max(0., temp * (2 - temp));
					rejf3 = 1 - br * sinthe / (1 + br * br);
					r3 = RandomCollection.random01();
					// ------------------------------------------
					if (r3 < rejf3) {
						break;
						// if((br < 1./broi) || (br > 1))
						// {
						// n_acc=0;
						// }
						// else
						// {
						// n_acc++;
						// costhe = 1 - temp;
						// theta=Math.acos(costhe);
						// }
					}
					// otherwise, loop is repeated until it finds an appropriate
					// value for theta
				}
			} else// At low energies it is faster to sample br uniformely
			{
				double bro = 1. / broi;
				double bro1 = 1 - bro;
				double rejmax = broi + bro;

				while (true)// n_acc==0)
				{

					r1 = RandomCollection.random01();
					r2 = RandomCollection.random01();
					br = bro + bro1 * r1;
					temp = (1 - br) / (ko * br);
					sinthe = Math.max(0., temp * (2 - temp));
					rejf3 = (br + 1. / br - sinthe) / rejmax;
					// ------------------------------------------
					if (r2 < rejf3) {
						break;
						// if((br < 1./broi) || (br > 1))
						// {
						// n_acc=0;
						// }
						// else
						// {
						// n_acc++;
						// costhe = 1 - temp;
						// theta=Math.acos(costhe);
						// }
					}
					// otherwise, loop is repeated until it finds an appropriate
					// value for theta
				}
			}
			// if((br < 1./broi) || (br > 1))
			if ((br >= 1. / broi) && (br <= 1)) {
				// if( (br < 0.99999/broi) || (br > 1.00001 ))
				// {
				// write(6,*) ' sampled br outside of allowed range!
				// ',ko,1./broi,br;
				// }
				// goto :RESAMPLE: ;
				break;
			}
		}// main loop added
		costhe = 1 - temp;
		theta = Math.acos(costhe);

		return theta;
	}

	// -------------------------------------------------------------------------------------
	// function to simulate Compton scattering.
	// input:energy of photon
	// output: scattering angle of photon in radians
	// based on code written by Lesley Buckley Dec17,2001
	/**
	 * Return polar angle using Klein-Nishina sampling based on classic algorithm.  
	 * @param k energy in MeV
	 * @return the result
	 */
	public static double comptonSimC(double k)// @@@@@@@@@@@@@@@@@@@@@@@CLASIC
	{
		int n_acc = 0; // number of accepted Compton events
		double theta = 0.0;// the output value!!
		double theta_temp = 0.0;// temporary value
		double u = 0.0;// 1-u=cos(theta)
		double kprime = 0.0;// energy of scattered photon
		double m = 0.511;// rest mass of electron in MeV->mc2
		double r = 0.0;// store the random generated number used in simulation
		// -------------------------------------------------------------------
		// want to be sure that Compton events would follow the Klein-Nishina
		// distribution!!
		double w = 0.0;// the Klein-Nishina based weighted factor!
		double wmax = 2.0;// the maximum Klein-Nishina based weighted factor!
		// -------------------------------------------------------------------
		while (n_acc == 0) {
			r = RandomCollection.random01();
			u = (m / k) * (-1 + Math.pow(1 + 2 * k / m, r));
			kprime = k / (1 + u * k / m);
			theta_temp = Math.acos(1 - u);// radians
			w = 1 + (kprime / k)
					* ((kprime / k) - Math.pow(Math.sin(theta_temp), 2));
			if (w > wmax * r) {
				n_acc++;
				theta = theta_temp;
			}
			// otherwise, loop is repeated until it finds an appropriate value
			// for theta
		}

		return theta;
	}

	// based on Wielopolski algorithm (1987)--O.Sima book
	/**
	 * Return polar angle using Klein-Nishina sampling based on Wielopolski algorithm.  
	 * @param k energy in MeV
	 * @return the result
	 */
	public static double comptonSimW(double k) {
		double theta = 0.0;// the output value!!
		double en = k * 1000;// in keV here
		if (en == 0.0 || en == 1.0)
			return theta;
		en = Math.log(en);// natural log
		double[][] coef1 = (double[][]) resources
				.getObject("kn.wielopolski.r1");
		double[][] coef2 = (double[][]) resources
				.getObject("kn.wielopolski.r2");
		double[][] coef3 = (double[][]) resources
				.getObject("kn.wielopolski.r3");
		double[] a = new double[4];
		double[] b = new double[7];
		// double a0,a1,a2,a3=0.0;
		double r = RandomCollection.random01();
		if (r <= 0.39) {
			for (int j = 0; j < 4; j++) {
				for (int i = 0; i < 7; i++) {
					b[i] = coef1[j][i];
				}

				a[j] = b[0] + b[1] * r + b[2] * r * r + b[3] * r * r * r + b[4]
						* r * r * r * r + b[5] * Math.sqrt(r);
			}

		} else if (r <= 0.7) {
			for (int j = 0; j < 4; j++) {
				for (int i = 0; i < 7; i++) {
					b[i] = coef2[j][i];
				}

				a[j] = b[0] + b[1] * r + b[2] * r * r + b[3] * r * r * r + b[4]
						* Math.exp(-b[5] * (r - b[6]) * (r - b[6]));
			}

		} else {
			r = 1 - r;
			for (int j = 0; j < 4; j++) {
				for (int i = 0; i < 7; i++) {
					b[i] = coef3[j][i];
				}

				a[j] = b[0] + b[1] * r + b[2] * r * r + b[3] * r * r * r + b[4]
						* r * r * r * r + b[5] * Math.sqrt(r);
			}

		}
		theta = a[0] + a[1] / en + a[2] * en + a[3] * en * en;
		if (theta > 2 * Math.PI)
			theta = theta - 2 * Math.PI;
		if (theta > Math.PI)
			theta = Math.PI - (theta - Math.PI);

		return theta;
	}

	// based on Kahn algorithm (old)--O.Sima book
	// is good for <1MeV
	/**
	 * Return polar angle using Klein-Nishina sampling based on Kahn algorithm.  
	 * @param k energy in MeV
	 * @return the result
	 */
	public static double comptonSimK(double k) {
		int n_acc = 0; // number of accepted Compton events
		double theta = 0.0;// the output value!!
		double miu = 0.0;
		double m = 0.511;// rest mass of electron in MeV->mc2
		double r = RandomCollection.random01();
		k = k / m;
		// -------------------------------------------------------------------
		double d = (2 * k + 1) / (2 * k + 9);
		while (n_acc == 0) {
			if (r <= d) {
				double y = 1 + 2 * k * r;
				r = RandomCollection.random01();
				double dd = 4 * (1 / y - 1 / (y * y));
				if (r <= dd) {
					miu = 1 - (y - 1) / k;
					n_acc++;
				}
			} else {
				double y = (1 + 2 * k) / (1 + 2 * k * r);
				miu = 1 - (y - 1) / k;
				r = RandomCollection.random01();
				double dd = 0.5 * (miu * miu + 1 / y);
				if (r <= dd) {
					n_acc++;
				}
			}
		}
		theta = Math.acos(miu);// radians
		return theta;
	}

	// ------------------
	
	@Deprecated
	public static double rayleighSim(double k) {
		int n_acc = 0; // number of accepted Compton events
		double theta = 0.0;// the output value!!
		double theta_temp = 0.0;// temporary value
		double u = 0.0;// 1-u=cos(theta)
		double r = 0.0;// store the random generated number used in simulation
		// -------------------------------------------------------------------
		// want to be sure that Compton events would follow the RAYLEIGH
		// distribution!!
		double w = 0.0;// the RAYLEIGH based weighted factor!
		double wmax = 2.0;// the maximum RAYLEIGH based weighted factor!
		// -------------------------------------------------------------------
		while (n_acc == 0) {
			r = RandomCollection.random01();
			u = wmax * r;
			theta_temp = Math.acos(1 - u);// radians
			w = 1 + (Math.pow(Math.cos(theta_temp), 2));
			if (w > wmax * r) {
				n_acc++;
				theta = theta_temp;
			}
			// otherwise, loop is repeated until it finds an appropriate value
			// for theta
		}

		return theta;
	}

	// -----------------------

	// input System.currentTimeMillis() at start of event->e.g. simulation
	// output : time elapsed in String format
	@Deprecated
	public static String timeElapsed(long startTime) {
		String times = "";
		long currentTime = System.currentTimeMillis();
		// tipul int spre deosebire de long--->aprox 8-9 zile suporta ca numar
		// de sec!!
		int delta = (new Long(currentTime - startTime)).intValue();// time
																	// elapsed
																	// in
																	// milliseconds
		int sec = delta / 1000;// impartire intreaga->catul!!!
		int milis = delta % 1000;// restul impartirii intregi!!
		if (sec > 60) {
			int min = sec / 60;
			sec = sec % 60;
			if (min > 60) {
				int h = min / 60;
				min = min % 60;
				if (h > 24) {
					int z = h / 24;
					h = h % 24;
					times = z + " days " + h + " h, " + min + " min, " + sec
							+ " sec, " + milis + " milis";
				} else {
					times = h + " h, " + min + " min, " + sec + " sec, "
							+ milis + " milis";
				}
			} else {
				times = min + " min, " + sec + " sec, " + milis + " milis";
			}
		} else {
			times = sec + " sec, " + milis + " milis";
		}
		return times;
	}

	// calculeaza m si n din cele 2 puncte p1(x1,y1) si p2(x2,y2)
	// si returneaza y(x).
	// pentru actualul mod nu se poate da eroare
	@Deprecated
	public static double linInt(double x1, double y1, double x2, double y2,
			double x) {
		double result = -1.0;
		double[] mn = new double[2];
		// insucces
		mn[0] = -1.0;// m
		mn[1] = -1.0;// n
		double num = x1 - x2;
		if (num != 0.0) {
			mn[0] = (y1 - y2) / num;
			mn[1] = (x1 * y2 - y1 * x2) / num;
			result = mn[0] * x + mn[1];
		}
		return result;
	}

	
	// --GENERAL COMPUTATION OF NEW COORDINATE @@!!!!!!!!!!!!!!!!
	@Deprecated//Use EGS class instead!!
	public static double[] getUpdatedCoordinate(double theta0, double phi0,
			double dist, double[] oldcoord) {
		double[] result = new double[6];
		// theta: 0 - pi so sin(theta) is always positive
		// phi: 0 - 2pi
		// for 0-pi sin(phi) is +
		// for pi-2pi sin(phi) is -
		// ---------computing new directional cosines--------
		double cost = Math.cos(theta0);
		double sint = Math.sqrt(1.0 - cost * cost);// sqrt is faster than
													// sin!!!!
		double cosp = Math.cos(phi0);
		double sinp = 0.0;
		double ux0 = oldcoord[3];// miux
		double uy0 = oldcoord[4];// miuy
		double uz0 = oldcoord[5];// miuz
		double x = oldcoord[0];
		double y = oldcoord[1];
		double z = oldcoord[2];
		double ux = 0.0;// miux
		double uy = 0.0;// miuy
		double uz = 0.0;// miuz
		if (phi0 < Math.PI)
			sinp = Math.sqrt(1.0 - cosp * cosp);
		else
			sinp = -Math.sqrt(1.0 - cosp * cosp);

		if (Math.abs(uz0) > 0.99999)// miuz->normal incident
		{
			ux = sint * cosp;// new miux
			uy = sint * sinp;// new miuy
			if (uz0 < 0)
				uz = -cost;// new miuz
			else
				uz = cost;// new miuz
		} else {
			double temp = Math.sqrt(1.0 - uz0 * uz0);
			ux = sint * (ux0 * uz0 * cosp - uy0 * sinp) / temp + ux0 * cost;
			uy = sint * (uy0 * uz0 * cosp + ux0 * sinp) / temp + uy0 * cost;
			uz = -sint * cosp * temp + uz0 * cost;
		}
		// --------computing the new
		// ccordinate------------------------------------------
		x = x + dist * ux;
		y = y + dist * uy;
		z = z + dist * uz;
		result[0] = x;
		result[1] = y;
		result[2] = z;
		result[3] = ux;
		result[4] = uy;
		result[5] = uz;

		return result;
	}

	
	// this will simulate the volumic source (source sampling)!!
	@Deprecated//Use EGS class instead!!
	private static double[] getCylinderRandom() {
		double[] result = new double[6];
		double r = RandomCollection.random01();
		// center of SC is in the center of cylinder (source shape), det is
		// placed at -hsource/2!!
		double z1 = hsource * r - hsource / 2;// initial random z coord in
												// cylinder coordinates
		r = RandomCollection.random01();
		double ro1 = asource * Math.sqrt(r);// dist to z axis ro1:
		r = RandomCollection.random01();
		double phi1 = 2 * Math.PI * r;// random angle for x0,y0 eval. from phi1
										// and ro1!!
		double d = hsource / 2 + z1;// distance point 0 to detector----is >0
									// always
		// for MC variance reduction we impose that all photons emerged from
		// point 0 will go into
		// an solid angle us1 having a distance d and a radius=2*detradius or
		// 2*sourceradius in order
		// to make sure that this allways contain the detector!!
		double diam = 0.0;
		if (asource <= adet) {
			diam = 2 * adet;
		} else {
			diam = 2 * asource;
		}
		double us1 = 2 * Math.PI * (1 - d / Math.sqrt(d * d + diam * diam));
		double costmax = d / Math.sqrt(d * d + diam * diam);// cosines of theta
															// max.
		// costet evaluation-->polar angle
		double dom = (1 - costmax) / 2;// >0 and <1/2, costmax <1 and >0,
										// thetamax<90!
		r = RandomCollection.random01();
		r = r * dom;
		double costet = 2 * r - 1;// <0 always--- negativ z axis!!
		double sintet = Math.sqrt(1 - costet * costet);
		double tgtet = sintet / costet;
		r = RandomCollection.random01();
		double phi2 = 2 * Math.PI * r;// azimutal angle
		pond = us1 / (4 * Math.PI);// initial weight associated with this
									// variance reduction
		double teta = Math.abs(Math.atan(tgtet));// -pi/2,pi/2
		// initial we have formal directional cosines u,v,w = 0,0,-1 and then we
		// have teta and phi2:
		double u = Math.sin(teta) * Math.cos(phi2);
		double v = Math.sin(teta) * Math.sin(phi2);
		double w = costet;// <0
		source_parcurs = d / costet;
		source_parcurs = Math.abs(source_parcurs);// >0
		double x0 = ro1 * Math.cos(phi1);
		double y0 = ro1 * Math.sin(phi1);
		double x = x0 + u * source_parcurs;
		double y = y0 + v * source_parcurs;
		double z = -hsource / 2;
		// attenuation in source:
		// double l1=Math.abs((asource-ro1)/sintet);//strictly in source volume
		// if (l1<source_parcurs) //it also fly in air (neglected)
		// source_parcurs=l1;
		// -----------------------------------
		double ux = u;
		double uy = v;// double uz=w;
		// strictly in source volume
		double l1 = 0.0;
		double zer = (ux * x0 + uy * y0) * (ux * x0 + uy * y0)
				+ (ux * ux + uy * uy) * (asource * asource - x0 * x0 - y0 * y0);
		if (zer > 0.) {
			double s1 = (-(ux * x0 + uy * y0) - Math.sqrt((ux * x0 + uy * y0)
					* (ux * x0 + uy * y0) + (ux * ux + uy * uy)
					* (asource * asource - x0 * x0 - y0 * y0)))
					/ (ux * ux + uy * uy);
			double s2 = (-(ux * x0 + uy * y0) + Math.sqrt((ux * x0 + uy * y0)
					* (ux * x0 + uy * y0) + (ux * ux + uy * uy)
					* (asource * asource - x0 * x0 - y0 * y0)))
					/ (ux * ux + uy * uy);
			double s = 0.;
			if ((s1 < 0.) && (s2 > 0)) {
				s = s2;
			} else {
				s = Math.min(Math.abs(s1), Math.abs(s2));
			}
			l1 = s;// System.out.println(l1);
		} else
			l1 = source_parcurs;// something wrong is happened!!
		if (l1 < source_parcurs) // it also fly in air (neglected)
			source_parcurs = l1;
		// ----------------------------------
		// new pond based on attenuation in source volume will be based on this
		// source_parcurs!!
		result[0] = x;
		result[1] = y;
		result[2] = z;
		result[3] = u;
		result[4] = v;
		result[5] = w;

		return result;
	}

	
	// ------------------------------------------------------------------------------------------
	// hwin=total thickness of absorbers
	@Deprecated//Use EGS class instead!!
	private static double[] getCylinderRandom2(double delzup, String[] wtype,
			double[] wthick, double e_inci) {
		double hwin = delzup;// 0.5=(7.3-6.3)/2
		double[] result = new double[6];
		double r = RandomCollection.random01();
		// center of SC is in the center of cylinder (source shape), det is
		// placed at -hsource/2!!
		double z1 = hsource * r - hsource / 2;// initial random z coord in
												// cylinder coordinates
		r = RandomCollection.random01();
		double ro1 = asource * Math.sqrt(r);// dist to z axis ro1:
		r = RandomCollection.random01();
		double phi1 = 2 * Math.PI * r;// random angle for x0,y0 eval. from phi1
										// and ro1!!
		double d = hsource / 2 + z1 + hwin;// distance point 0 to detector----is
											// >0 always
		double dso = hsource / 2 + z1;// distance point to out of source
		// for MC variance reduction we impose that all photons emerged from
		// point 0 will go into
		// an solid angle us1 having a distance d and a radius=2*detradius or
		// 2*sourceradius in order
		// to make sure that this allways contain the detector!!
		double diam = 0.0;
		if (asource <= adet) {
			diam = 2 * adet;
		} else {
			diam = 2 * asource;
		}
		double us1 = 2 * Math.PI * (1 - d / Math.sqrt(d * d + diam * diam));
		double costmax = d / Math.sqrt(d * d + diam * diam);// cosines of theta
															// max.
		// costet evaluation-->polar angle
		double dom = (1 - costmax) / 2;// >0 and <1/2, costmax <1 and >0,
										// thetamax<90!
		r = RandomCollection.random01();
		r = r * dom;
		double costet = 2 * r - 1;// <0 always--- negativ z axis!!
		double sintet = Math.sqrt(1 - costet * costet);
		double tgtet = sintet / costet;
		r = RandomCollection.random01();
		double phi2 = 2 * Math.PI * r;// azimutal angle
		pond = us1 / (4 * Math.PI);// initial weight associated with this
									// variance reduction
		double teta = Math.abs(Math.atan(tgtet));// -pi/2,pi/2
		// initial we have formal directional cosines u,v,w = 0,0,-1 and then we
		// have teta and phi2:
		double u = Math.sin(teta) * Math.cos(phi2);
		double v = Math.sin(teta) * Math.sin(phi2);
		double w = costet;// <0
		source_parcurs = d / costet;
		source_parcurs = Math.abs(source_parcurs);// >0
		double x0 = ro1 * Math.cos(phi1);
		double y0 = ro1 * Math.sin(phi1);
		double x = x0 + u * source_parcurs;
		double y = y0 + v * source_parcurs;
		// double z=-hsource/2;
		double z = -hsource / 2 - hwin;// @@@@@@@@@@@@@@@@@@@@@@@@@@
		double source_parcurs2 = dso / costet;
		source_parcurs2 = Math.abs(source_parcurs2);// >0
		// attenuation in source:
		// double l1=Math.abs((asource-ro1)/sintet);//strictly in source volume
		// if (l1<source_parcurs) //it also fly in air (neglected)
		// source_parcurs=l1;
		// -----------------------------------
		double ux = u;
		double uy = v;// double uz=w;
		// strictly in source volume
		double l1 = 0.0;
		double zer = (ux * x0 + uy * y0) * (ux * x0 + uy * y0)
				+ (ux * ux + uy * uy) * (asource * asource - x0 * x0 - y0 * y0);
		if (zer > 0.) {
			double s1 = (-(ux * x0 + uy * y0) - Math.sqrt((ux * x0 + uy * y0)
					* (ux * x0 + uy * y0) + (ux * ux + uy * uy)
					* (asource * asource - x0 * x0 - y0 * y0)))
					/ (ux * ux + uy * uy);
			double s2 = (-(ux * x0 + uy * y0) + Math.sqrt((ux * x0 + uy * y0)
					* (ux * x0 + uy * y0) + (ux * ux + uy * uy)
					* (asource * asource - x0 * x0 - y0 * y0)))
					/ (ux * ux + uy * uy);
			double s = 0.;
			if ((s1 < 0.) && (s2 > 0)) {
				s = s2;
			} else {
				s = Math.min(Math.abs(s1), Math.abs(s2));
			}
			l1 = s;// System.out.println(l1);
		} else
			l1 = source_parcurs2;// something wrong is happened!!
		if (l1 < source_parcurs2) // it also fly in air (neglected)
			source_parcurs2 = l1;
		// detector attenuation---------------------------------------------
		// double ldetatt=source_parcurs-source_parcurs2;//no use only in Al
		// double windens=0.;
		// double winatt=0.;
		if (!geomB) {
			wrho = solvewinDens(wtype);
			watt = solvewinAtt(wtype, e_inci);
		}

		for (int i = 0; i < wthick.length; i++) {
			// double wrho=solvewinDens(wtype[i]);
			// double watt=solvewinAtt(wtype[i],e_inci);
			// pond=pond*Math.exp(-wrho*watt*wthick[i]);
			pond = pond
					* Math.exp(-wrho[i] * watt[i] * wthick[i]
							/ Math.abs(costet));
			// pondwin
		}
		// -------------------------------------------------------------------
		source_parcurs = source_parcurs2;// reset to in source attenuation
		// ----------------------------------
		// new pond based on attenuation in source volume will be based on this
		// source_parcurs!!
		result[0] = x;
		result[1] = y;
		result[2] = z;
		result[3] = u;
		result[4] = v;
		result[5] = w;
		// test=test+pond;
		return result;
	}

	
	@Deprecated//Use EGS class instead!!
	public static double[] solvewinDens(String[] s) {
		double[] result = new double[s.length];
		double[] windensity = (double[]) resources
				.getObject("window.densities");
		String[] wintype = (String[]) resources.getObject("window.type");
		for (int j = 0; j < s.length; j++) {
			for (int i = 0; i < wintype.length; i++) {
				if (s[j].compareTo(wintype[i]) == 0) {
					result[j] = windensity[i];
					break;
				}
			}
		}
		geomB = true;
		return result;
	}

	
	@Deprecated//Use EGS class instead!!
	public static double[] solvewinAtt(String[] s, double e_inc) {
		double[] result = new double[s.length];
		int index_low = 0;
		int index_high = 0;
		double[][] wintable = (double[][]) resources
				.getObject("detector.window.attenuationCoef");
		double[] energy = new double[wintable.length];
		// double[] total_coeff= new double[wintable.length];
		String[] wintype = (String[]) resources.getObject("window.type");
		for (int jj = 0; jj < s.length; jj++) {
			double[] total_coeff = new double[wintable.length];
			for (int i = 0; i < energy.length; i++) {
				energy[i] = wintable[i][0];
				for (int j = 0; j < wintype.length; j++) {
					if (s[jj].compareTo(wintype[j]) == 0) {
						total_coeff[i] = wintable[i][j + 1];
						break;
					}
				}
			}
			Sort.findNearestValue(energy, e_inc, true);
			index_low = Sort.getNearestPosition();
			if (index_low < energy.length - 1)
				index_high = index_low + 1;
			else
				index_high = index_low;
			result[jj] = linInt(energy[index_high], total_coeff[index_high],
					energy[index_low], total_coeff[index_low], e_inc);
		}
		return result;
	}

	
	@Deprecated//Use EGS class instead!!
	public static double solvewinDens(String s) {
		double result = 0.;
		double[] windensity = (double[]) resources
				.getObject("window.densities");
		String[] wintype = (String[]) resources.getObject("window.type");
		for (int i = 0; i < wintype.length; i++) {
			if (s.compareTo(wintype[i]) == 0) {
				result = windensity[i];
				break;
			}
		}

		return result;
	}

	
	@Deprecated//Use EGS class instead!!
	public static double solvewinAtt(String s, double e_inc) {
		double result = 0.;
		int index_low = 0;
		int index_high = 0;
		double[][] wintable = (double[][]) resources
				.getObject("detector.window.attenuationCoef");
		double[] energy = new double[wintable.length];
		double[] total_coeff = new double[wintable.length];
		String[] wintype = (String[]) resources.getObject("window.type");
		for (int i = 0; i < energy.length; i++) {
			energy[i] = wintable[i][0];
			for (int j = 0; j < wintype.length; j++) {
				if (s.compareTo(wintype[j]) == 0) {
					total_coeff[i] = wintable[i][j + 1];
					break;
				}
			}

			// if (s.compareTo("Al")==0)
			// {
			// total_coeff[i]=wintable[i][1];
			// }
		}
		Sort.findNearestValue(energy, e_inc, true);
		index_low = Sort.getNearestPosition();
		if (index_low < energy.length - 1)
			index_high = index_low + 1;
		else
			index_high = index_low;
		result = linInt(energy[index_high], total_coeff[index_high],
				energy[index_low], total_coeff[index_low], e_inc);
		return result;
	}

	
	// double[] wthick2==> pe orizontala
	@Deprecated//Use EGS class instead!!
	private static double[] getMarrinelliRandom2(double deltazup,
			double deltazdown, double deltaxy, String[] wtype, double[] wthick,
			double e_inci) {
		double[] result = new double[6];
		double r = RandomCollection.random01();
		double z1 = hsource * r - hsource / 2;// initial random z coord in
												// cylinder coordinates
		boolean infcyl = false;// test if inf cylinder exists
		// if(hsource>hsourceup+hdet)//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
		if (hsource > hsourceup + hdet + deltazup + deltazdown)
			infcyl = true;
		// if(z1>=hsource/2-hsourceup)//sup cylinder
		// case//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
		if (z1 >= hsource / 2 - hsourceup)// -deltazup)
		{
			result = getMCylSup2(z1, deltazup, wtype, wthick, e_inci);
		} else if (infcyl) {
			// if (z1<=hsource/2-hsourceup-hdet)//inf cylinder
			// case//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
			if (z1 <= hsource / 2 - hsourceup - hdet - deltazup - deltazdown) {
				result = getMCylInf2(z1, deltazup, deltazdown);
			} else// middle region
			{
				result = getMMidInf2(z1, deltaxy);
			}
		} else// middle region
		{
			result = getMMidInf2(z1, deltaxy);
		}

		return result;
	}
	
	
	@Deprecated//Use EGS class instead!!
	private static double[] getMCylSup2(double zi, double delzup,
			String[] wtype, double[] wthick, double e_inci) {
		double hwin = delzup;
		double[] result = new double[6];
		double z1 = zi;// initial random z coord in cylinder coordinates
		double r = RandomCollection.random01();
		double ro1 = asource * Math.sqrt(r);// dist to z axis ro1:
		r = RandomCollection.random01();
		double phi1 = 2 * Math.PI * r;// random angle for x0,y0 eval. from phi1
										// and ro1!!
		double d = z1 - (hsource / 2 - hsourceup - hwin);// distance point 0 to
															// detector>0
		double dso = z1 - (hsource / 2 - hsourceup);
		;// distance point to out of source
		// for MC variance reduction we impose that all photons emerged from
		// point 0 will go into
		// an solid angle us1 having a distance d and a radius=2*sourceradius in
		// order
		// to make sure that this allways contain the detector!!
		double diam = 2 * asource;
		double us1 = 2 * Math.PI * (1 - d / Math.sqrt(d * d + diam * diam));
		double costmax = d / Math.sqrt(d * d + diam * diam);// cosines of theta
															// max.
		// costet evaluation-->polar angle
		double dom = (1 - costmax) / 2;// >0 and <1/2, costmax <1 and >0,
										// thetamax<90!
		r = RandomCollection.random01();
		r = r * dom;
		double costet = 2 * r - 1;// <0 always--- negativ z axis!!
		double sintet = Math.sqrt(1 - costet * costet);
		double tgtet = sintet / costet;
		r = RandomCollection.random01();
		double phi2 = 2 * Math.PI * r;// azimutal angle
		pond = us1 / (4 * Math.PI);// initial weight associated with this
									// variance reduction
		double teta = Math.abs(Math.atan(tgtet));// -pi/2,pi/2
		// initial we have formal directional cosines u,v,w = 0,0,-1 and then we
		// have teta and phi2:
		double u = Math.sin(teta) * Math.cos(phi2);
		double v = Math.sin(teta) * Math.sin(phi2);
		double w = costet;// <0
		source_parcurs = d / costet;
		source_parcurs = Math.abs(source_parcurs);// >0
		double x0 = ro1 * Math.cos(phi1);
		double y0 = ro1 * Math.sin(phi1);
		double x = x0 + u * source_parcurs;
		double y = y0 + v * source_parcurs;
		// double z=hsource/2-hsourceup;-hwin;
		double z = hsource / 2 - hsourceup - hwin;// @@@@@@@@@@@@@@@@@@@@@@@@@@
		double source_parcurs2 = dso / costet;
		source_parcurs2 = Math.abs(source_parcurs2);// >0
		// attenuation is always in source:
		// detector attenuation---------------------------------------------
		// double ldetatt=source_parcurs-source_parcurs2;//no use only in Al
		// double windens=0.;
		// double winatt=0.;
		if (!geomB) {
			wrho = solvewinDens(wtype);
			watt = solvewinAtt(wtype, e_inci);
		}

		for (int i = 0; i < wthick.length; i++) {
			// double wrho=solvewinDens(wtype[i]);
			// double watt=solvewinAtt(wtype[i],e_inci);
			// pond=pond*Math.exp(-wrho*watt*wthick[i]);
			pond = pond
					* Math.exp(-wrho[i] * watt[i] * wthick[i]
							/ Math.abs(costet));
			// pondwin
		}
		// -------------------------------------------------------------------
		source_parcurs = source_parcurs2;// reset to in source attenuation
		// ----------------------------------
		result[0] = x;
		result[1] = y;
		result[2] = z;
		result[3] = u;
		result[4] = v;
		result[5] = w;

		return result;
	}
	
	@Deprecated//Use EGS class instead!!
	private static double[] getMCylInf2(double zi, double delup, double deldown) {
		// from behind=>no change!!
		double[] result = new double[6];
		double z1 = zi;// initial random z coord in cylinder coordinates
		double r = RandomCollection.random01();
		double ro1 = Math.sqrt((asource * asource - bsource * bsource) * r
				+ bsource * bsource);// dist to z axis ro1:
		r = RandomCollection.random01();
		double phi1 = 2 * Math.PI * r;// random angle for x0,y0 eval. from phi1
										// and ro1!!
		// double d=Math.abs(z1+(hsource/2-hsourceup-hdet));//distance point 0
		// to detector>0
		// double d=Math.abs(z1+(hsource/2-hsourceup-hdet-delup));
		// double dso=Math.abs(z1+(hsource/2-hsourceup-hdet-delup-deldown));
		double d = Math.abs(z1 - (hsource / 2 - hsourceup - hdet - delup));
		double dso = Math.abs(z1
				- (hsource / 2 - hsourceup - hdet - delup - deldown));
		// for MC variance reduction we impose that all photons emerged from
		// point 0 will go into
		// an solid angle us1 having a distance d and a radius=2*sourceradius in
		// order
		// to make sure that this allways contain the detector!!
		double diam = 2 * asource;
		double us1 = 2 * Math.PI * (1 - d / Math.sqrt(d * d + diam * diam));
		double costmax = d / Math.sqrt(d * d + diam * diam);// cosines of theta
															// max.
		// costet evaluation-->polar angle
		double dom = (1 - costmax) / 2;// >0 and <1/2, costmax <1 and >0,
										// thetamax<90!
		r = RandomCollection.random01();
		r = r * dom;
		double costet = -2 * r + 1;// >0 always--- positiv z axis!!!!!!!!!!!!!!
		double sintet = Math.sqrt(1 - costet * costet);
		double tgtet = sintet / costet;
		r = RandomCollection.random01();
		double phi2 = 2 * Math.PI * r;// azimutal angle
		pond = us1 / (4 * Math.PI);// initial weight associated with this
									// variance reduction
		double teta = Math.abs(Math.atan(tgtet));// -pi/2,pi/2
		// initial we have formal directional cosines u,v,w = 0,0,-1 and then we
		// have teta and phi2:
		double u = Math.sin(teta) * Math.cos(phi2);
		double v = Math.sin(teta) * Math.sin(phi2);
		double w = costet;// >0
		source_parcurs = d / costet;
		source_parcurs = Math.abs(source_parcurs);// >0
		double x0 = ro1 * Math.cos(phi1);
		double y0 = ro1 * Math.sin(phi1);
		double x = x0 + u * source_parcurs;
		double y = y0 + v * source_parcurs;
		// double z=hsource/2-hsourceup-hdet;
		double z = hsource / 2 - hsourceup - hdet - delup;
		double source_parcurs2 = dso / costet;
		source_parcurs2 = Math.abs(source_parcurs2);// >0
		// attenuation !is always in source:
		// double l1=Math.abs((-bsource+ro1)/sintet);//strictly in source volume
		// if (l1<source_parcurs) //it also fly in air (neglected)
		// source_parcurs=l1;
		// -----------------------------------
		double ux = u;
		double uy = v;// double uz=w;
		// strictly in source volume
		double l1 = 0.0;
		double zer = (ux * x0 + uy * y0) * (ux * x0 + uy * y0)
				+ (ux * ux + uy * uy) * (bsource * bsource - x0 * x0 - y0 * y0);
		if (zer > 0.) {
			double s1 = (-(ux * x0 + uy * y0) - Math.sqrt((ux * x0 + uy * y0)
					* (ux * x0 + uy * y0) + (ux * ux + uy * uy)
					* (bsource * bsource - x0 * x0 - y0 * y0)))
					/ (ux * ux + uy * uy);
			double s2 = (-(ux * x0 + uy * y0) + Math.sqrt((ux * x0 + uy * y0)
					* (ux * x0 + uy * y0) + (ux * ux + uy * uy)
					* (bsource * bsource - x0 * x0 - y0 * y0)))
					/ (ux * ux + uy * uy);
			double s = 0.;
			if ((s1 < 0.) && (s2 > 0)) {
				s = s2;
			} else {
				s = Math.min(Math.abs(s1), Math.abs(s2));
			}
			l1 = s;// System.out.println(l1);
		} else
			l1 = source_parcurs2;// something wrong is happened!!
		if (l1 < source_parcurs2) // it also fly in air (neglected)
			source_parcurs2 = l1;
		// ----------------------------------
		source_parcurs = source_parcurs2;// reset to in source attenuation
		// ----------------------------------
		result[0] = x;
		result[1] = y;
		result[2] = z;
		result[3] = u;
		result[4] = v;
		result[5] = w;

		return result;
	}
	
	@Deprecated//Use EGS class instead!!
	private static double[] getMMidInf2(double zi, double delxy) {
		double[] result = new double[6];
		double z1 = zi;// initial random z coord in cylinder coordinates
		double r = RandomCollection.random01();
		double ro1 = Math.sqrt((asource * asource - bsource * bsource) * r
				+ bsource * bsource);// dist to z axis ro1:
		r = RandomCollection.random01();
		double phi1 = 2 * Math.PI * r;// ///////for x0 and y0 eval
										// ???????????????????/
		double rdet = Math.sqrt(adet * adet + hdet * hdet / 4);// for us eval!!
		double d = ro1;// distance for us eval!!
		// to make sure that this allways contain the detector!!
		double diam = 2 * rdet;// for us eval!!
		double us1 = 2 * Math.PI * (1 - d / Math.sqrt(d * d + diam * diam));// OK!
		double costmax = d / Math.sqrt(d * d + diam * diam);// cosines of theta
															// max.
		// costet evaluation-->polar angle
		double dom = (1 - costmax) / 2;// >0 and <1/2, costmax <1 and >0,
										// thetamax<90!
		r = RandomCollection.random01();
		r = r * dom;
		double costet = -2 * r + 1;// >0 for instance
		double sintet = Math.sqrt(1 - costet * costet);
		double tgtet = sintet / costet;
		r = RandomCollection.random01();
		// for simplicity and due to simmetry of problem we choose that ***
		double phi2 = 2 * Math.PI * r;// phi1;//2*Math.PI*r;//azimutal angle
		pond = us1 / (4 * Math.PI);// initial weight associated with this
									// variance reduction
		double teta = Math.abs(Math.atan(tgtet));// -pi/2,pi/2
		// initial we have formal directional cosines u,v,w and then we have
		// teta and phi2:
		double xsi = Math.PI / 2;// shift with Pi/2 relative to z axis!!!
		double u = Math.sin(xsi) * Math.cos(phi1 + Math.PI);// -changing
															// direction
		double v = Math.sin(xsi) * Math.sin(phi1 + Math.PI);// -changing
															// direction
		double w = Math.cos(xsi);// 0!
		// "new " direction
		double ux = 0.0;
		double uy = 0.0;
		double uz = 0.0;
		if (Math.abs(w) > 0.99999)// miuz->normal incident
		{
			ux = Math.sin(teta) * Math.cos(phi2);// new miux
			uy = Math.sin(teta) * Math.sin(phi2);// new miuy
			if (w < 0)
				uz = -Math.cos(teta);// new miuz--never here---random
			else
				uz = Math.cos(teta);// new miuz
		} else {
			double temp = Math.sqrt(1.0 - w * w);
			ux = Math.sin(teta) * (u * w * Math.cos(phi2) - v * Math.sin(phi2))
					/ temp + u * Math.cos(teta);
			uy = Math.sin(teta) * (v * w * Math.cos(phi2) + u * Math.sin(phi2))
					/ temp + v * Math.cos(teta);
			uz = -Math.sin(teta) * Math.cos(phi2) * temp + w * Math.cos(teta);
		}
		// --------------------------------------
		// source_parcurs=(d-adet)/costet;//minimum for phi2 particular see ***
		// source_parcurs=Math.abs(source_parcurs)+0.1;//>0 and be sure to hit
		// if it is tend to!!
		double x0 = ro1 * Math.cos(phi1);
		double y0 = ro1 * Math.sin(phi1);
		// -----------------------------------
		double zer = (ux * x0 + uy * y0) * (ux * x0 + uy * y0)
				+ (ux * ux + uy * uy) * (adet * adet - x0 * x0 - y0 * y0);
		if (zer > 0.) {
			// source_parcurs=-(ux*x0+uy*y0);//System.out.println(
			// source_parcurs);>0
			double s1 = (-(ux * x0 + uy * y0) - Math.sqrt((ux * x0 + uy * y0)
					* (ux * x0 + uy * y0) + (ux * ux + uy * uy)
					* (adet * adet - x0 * x0 - y0 * y0)))
					/ (ux * ux + uy * uy);
			double s2 = (-(ux * x0 + uy * y0) + Math.sqrt((ux * x0 + uy * y0)
					* (ux * x0 + uy * y0) + (ux * ux + uy * uy)
					* (adet * adet - x0 * x0 - y0 * y0)))
					/ (ux * ux + uy * uy);
			double s = 0.;
			if ((s1 < 0.) && (s2 > 0)) {
				s = s2;
			} else {
				s = Math.min(Math.abs(s1), Math.abs(s2));
			}
			source_parcurs = s + 0.1;// System.out.println(source_parcurs);
			// aici - radical ca vreu prima atingere, aia mica, Nu a doua in
			// spate cilindru!!!!
			// source_parcurs=Math.abs(source_parcurs/(ux*ux+uy*uy))+0.025;//System.out.println(source_parcurs);
		} else
			source_parcurs = 0.;// not hit
			// ----------------------------------
		double x = x0 + ux * source_parcurs;
		double y = y0 + uy * source_parcurs;
		double z = z1 + uz * source_parcurs;
		// attenuation is always in source:
		// starting point
		result[0] = x;
		result[1] = y;
		result[2] = z;
		result[3] = ux;
		result[4] = uy;
		result[5] = uz;

		return result;
	}
	
	// ---------------------------------------------------------------------------------------
	@Deprecated//Use EGS class instead!!
	private static double[] getMarrinelliRandom() {
		double[] result = new double[6];
		double r = RandomCollection.random01();
		double z1 = hsource * r - hsource / 2;// initial random z coord in
												// cylinder coordinates
		boolean infcyl = false;// test if inf cylinder exists
		if (hsource > hsourceup + hdet)
			infcyl = true;
		if (z1 >= hsource / 2 - hsourceup)// sup cylinder case
		{
			result = getMCylSup(z1);
		} else if (infcyl) {
			if (z1 <= hsource / 2 - hsourceup - hdet)// inf cylinder case
			{
				result = getMCylInf(z1);
			} else// middle region
			{
				result = getMMidInf(z1);
			}
		} else// middle region
		{
			result = getMMidInf(z1);
		}

		return result;
	}
	
	@Deprecated//Use EGS class instead!!
	private static double[] getMCylSup(double zi) {
		double[] result = new double[6];
		double z1 = zi;// initial random z coord in cylinder coordinates
		double r = RandomCollection.random01();
		double ro1 = asource * Math.sqrt(r);// dist to z axis ro1:
		r = RandomCollection.random01();
		double phi1 = 2 * Math.PI * r;// random angle for x0,y0 eval. from phi1
										// and ro1!!
		double d = z1 - (hsource / 2 - hsourceup);// distance point 0 to
													// detector>0
		// for MC variance reduction we impose that all photons emerged from
		// point 0 will go into
		// an solid angle us1 having a distance d and a radius=2*sourceradius in
		// order
		// to make sure that this allways contain the detector!!
		double diam = 2 * asource;
		double us1 = 2 * Math.PI * (1 - d / Math.sqrt(d * d + diam * diam));
		double costmax = d / Math.sqrt(d * d + diam * diam);// cosines of theta
															// max.
		// costet evaluation-->polar angle
		double dom = (1 - costmax) / 2;// >0 and <1/2, costmax <1 and >0,
										// thetamax<90!
		r = RandomCollection.random01();
		r = r * dom;
		double costet = 2 * r - 1;// <0 always--- negativ z axis!!
		double sintet = Math.sqrt(1 - costet * costet);
		double tgtet = sintet / costet;
		r = RandomCollection.random01();
		double phi2 = 2 * Math.PI * r;// azimutal angle
		pond = us1 / (4 * Math.PI);// initial weight associated with this
									// variance reduction
		double teta = Math.abs(Math.atan(tgtet));// -pi/2,pi/2
		// initial we have formal directional cosines u,v,w = 0,0,-1 and then we
		// have teta and phi2:
		double u = Math.sin(teta) * Math.cos(phi2);
		double v = Math.sin(teta) * Math.sin(phi2);
		double w = costet;// <0
		source_parcurs = d / costet;
		source_parcurs = Math.abs(source_parcurs);// >0
		double x0 = ro1 * Math.cos(phi1);
		double y0 = ro1 * Math.sin(phi1);
		double x = x0 + u * source_parcurs;
		double y = y0 + v * source_parcurs;
		double z = hsource / 2 - hsourceup;
		// attenuation is always in source:
		result[0] = x;
		result[1] = y;
		result[2] = z;
		result[3] = u;
		result[4] = v;
		result[5] = w;

		return result;
	}
	
	@Deprecated//Use EGS class instead!!
	private static double[] getMCylInf(double zi) {
		double[] result = new double[6];
		double z1 = zi;// initial random z coord in cylinder coordinates
		double r = RandomCollection.random01();
		double ro1 = Math.sqrt((asource * asource - bsource * bsource) * r
				+ bsource * bsource);// dist to z axis ro1:
		r = RandomCollection.random01();
		double phi1 = 2 * Math.PI * r;// random angle for x0,y0 eval. from phi1
										// and ro1!!
		// double d=Math.abs(z1+(hsource/2-hsourceup-hdet));//distance point 0
		// to detector>0
		double d = Math.abs(z1 - (hsource / 2 - hsourceup - hdet));
		// for MC variance reduction we impose that all photons emerged from
		// point 0 will go into
		// an solid angle us1 having a distance d and a radius=2*sourceradius in
		// order
		// to make sure that this allways contain the detector!!
		double diam = 2 * asource;
		double us1 = 2 * Math.PI * (1 - d / Math.sqrt(d * d + diam * diam));
		double costmax = d / Math.sqrt(d * d + diam * diam);// cosines of theta
															// max.
		// costet evaluation-->polar angle
		double dom = (1 - costmax) / 2;// >0 and <1/2, costmax <1 and >0,
										// thetamax<90!
		r = RandomCollection.random01();
		r = r * dom;
		double costet = -2 * r + 1;// >0 always--- positiv z axis!!!!!!!!!!!!!!
		double sintet = Math.sqrt(1 - costet * costet);
		double tgtet = sintet / costet;
		r = RandomCollection.random01();
		double phi2 = 2 * Math.PI * r;// azimutal angle
		pond = us1 / (4 * Math.PI);// initial weight associated with this
									// variance reduction
		double teta = Math.abs(Math.atan(tgtet));// -pi/2,pi/2
		// initial we have formal directional cosines u,v,w = 0,0,-1 and then we
		// have teta and phi2:
		double u = Math.sin(teta) * Math.cos(phi2);
		double v = Math.sin(teta) * Math.sin(phi2);
		double w = costet;// >0
		source_parcurs = d / costet;
		source_parcurs = Math.abs(source_parcurs);// >0
		double x0 = ro1 * Math.cos(phi1);
		double y0 = ro1 * Math.sin(phi1);
		double x = x0 + u * source_parcurs;
		double y = y0 + v * source_parcurs;
		double z = hsource / 2 - hsourceup - hdet;
		// attenuation is always in source:
		// double l1=Math.abs((-bsource+ro1)/sintet);//strictly in source volume
		// if (l1<source_parcurs) //it also fly in air (neglected)
		// source_parcurs=l1;
		// -----------------------------------
		double ux = u;
		double uy = v;// double uz=w;
		// strictly in source volume
		double l1 = 0.0;
		double zer = (ux * x0 + uy * y0) * (ux * x0 + uy * y0)
				+ (ux * ux + uy * uy) * (bsource * bsource - x0 * x0 - y0 * y0);
		if (zer > 0.) {
			double s1 = (-(ux * x0 + uy * y0) - Math.sqrt((ux * x0 + uy * y0)
					* (ux * x0 + uy * y0) + (ux * ux + uy * uy)
					* (bsource * bsource - x0 * x0 - y0 * y0)))
					/ (ux * ux + uy * uy);
			double s2 = (-(ux * x0 + uy * y0) + Math.sqrt((ux * x0 + uy * y0)
					* (ux * x0 + uy * y0) + (ux * ux + uy * uy)
					* (bsource * bsource - x0 * x0 - y0 * y0)))
					/ (ux * ux + uy * uy);
			double s = 0.;
			if ((s1 < 0.) && (s2 > 0)) {
				s = s2;
			} else {
				s = Math.min(Math.abs(s1), Math.abs(s2));
			}
			l1 = s;// System.out.println(l1);
		} else
			l1 = source_parcurs;// something wrong is happened!!
		if (l1 < source_parcurs) // it also fly in air (neglected)
			source_parcurs = l1;
		// ----------------------------------
		result[0] = x;
		result[1] = y;
		result[2] = z;
		result[3] = u;
		result[4] = v;
		result[5] = w;

		return result;
	}
	
	@Deprecated//Use EGS class instead!!
	private static double[] getMMidInf(double zi) {
		double[] result = new double[6];
		double z1 = zi;// initial random z coord in cylinder coordinates
		double r = RandomCollection.random01();
		double ro1 = Math.sqrt((asource * asource - bsource * bsource) * r
				+ bsource * bsource);// dist to z axis ro1:
		r = RandomCollection.random01();
		double phi1 = 2 * Math.PI * r;// ///////for x0 and y0 eval
										// ???????????????????/
		double rdet = Math.sqrt(adet * adet + hdet * hdet / 4);// for us eval!!
		double d = ro1;// distance for us eval!!
		// to make sure that this allways contain the detector!!
		double diam = 2 * rdet;// for us eval!!
		double us1 = 2 * Math.PI * (1 - d / Math.sqrt(d * d + diam * diam));// OK!
		double costmax = d / Math.sqrt(d * d + diam * diam);// cosines of theta
															// max.
		// costet evaluation-->polar angle
		double dom = (1 - costmax) / 2;// >0 and <1/2, costmax <1 and >0,
										// thetamax<90!
		r = RandomCollection.random01();
		r = r * dom;
		double costet = -2 * r + 1;// >0 for instance
		double sintet = Math.sqrt(1 - costet * costet);
		double tgtet = sintet / costet;
		r = RandomCollection.random01();
		// for simplicity and due to simmetry of problem we choose that ***
		double phi2 = 2 * Math.PI * r;// phi1;//2*Math.PI*r;//azimutal angle
		pond = us1 / (4 * Math.PI);// initial weight associated with this
									// variance reduction
		double teta = Math.abs(Math.atan(tgtet));// -pi/2,pi/2
		// initial we have formal directional cosines u,v,w and then we have
		// teta and phi2:
		double xsi = Math.PI / 2;// shift with Pi/2 relative to z axis!!!
		double u = Math.sin(xsi) * Math.cos(phi1 + Math.PI);// -changing
															// direction
		double v = Math.sin(xsi) * Math.sin(phi1 + Math.PI);// -changing
															// direction
		double w = Math.cos(xsi);// 0!
		// "new " direction
		double ux = 0.0;
		double uy = 0.0;
		double uz = 0.0;
		if (Math.abs(w) > 0.99999)// miuz->normal incident
		{
			ux = Math.sin(teta) * Math.cos(phi2);// new miux
			uy = Math.sin(teta) * Math.sin(phi2);// new miuy
			if (w < 0)
				uz = -Math.cos(teta);// new miuz--never here---random
			else
				uz = Math.cos(teta);// new miuz
		} else {
			double temp = Math.sqrt(1.0 - w * w);
			ux = Math.sin(teta) * (u * w * Math.cos(phi2) - v * Math.sin(phi2))
					/ temp + u * Math.cos(teta);
			uy = Math.sin(teta) * (v * w * Math.cos(phi2) + u * Math.sin(phi2))
					/ temp + v * Math.cos(teta);
			uz = -Math.sin(teta) * Math.cos(phi2) * temp + w * Math.cos(teta);
		}
		// --------------------------------------
		// source_parcurs=(d-adet)/costet;//minimum for phi2 particular see ***
		// source_parcurs=Math.abs(source_parcurs)+0.1;//>0 and be sure to hit
		// if it is tend to!!
		double x0 = ro1 * Math.cos(phi1);
		double y0 = ro1 * Math.sin(phi1);
		// -----------------------------------
		double zer = (ux * x0 + uy * y0) * (ux * x0 + uy * y0)
				+ (ux * ux + uy * uy) * (adet * adet - x0 * x0 - y0 * y0);
		if (zer > 0.) {
			// source_parcurs=-(ux*x0+uy*y0);//System.out.println(
			// source_parcurs);>0
			double s1 = (-(ux * x0 + uy * y0) - Math.sqrt((ux * x0 + uy * y0)
					* (ux * x0 + uy * y0) + (ux * ux + uy * uy)
					* (adet * adet - x0 * x0 - y0 * y0)))
					/ (ux * ux + uy * uy);
			double s2 = (-(ux * x0 + uy * y0) + Math.sqrt((ux * x0 + uy * y0)
					* (ux * x0 + uy * y0) + (ux * ux + uy * uy)
					* (adet * adet - x0 * x0 - y0 * y0)))
					/ (ux * ux + uy * uy);
			double s = 0.;
			if ((s1 < 0.) && (s2 > 0)) {
				s = s2;
			} else {
				s = Math.min(Math.abs(s1), Math.abs(s2));
			}
			source_parcurs = s + 0.1;// System.out.println(source_parcurs);
			// aici - radical ca vreu prima atingere, aia mica, Nu a doua in
			// spate cilindru!!!!
			// source_parcurs=Math.abs(source_parcurs/(ux*ux+uy*uy))+0.025;//System.out.println(source_parcurs);
		} else
			source_parcurs = 0.;// not hit
			// ----------------------------------
		double x = x0 + ux * source_parcurs;
		double y = y0 + uy * source_parcurs;
		double z = z1 + uz * source_parcurs;
		// attenuation is always in source:
		// starting point
		result[0] = x;
		result[1] = y;
		result[2] = z;
		result[3] = ux;
		result[4] = uy;
		result[5] = uz;

		return result;
	}
	
	@Deprecated//Use EGS class instead!!
	private static boolean atSurfaceDet(double[] coord) {
		boolean b = true;
		if (isM) {
			double d2 = coord[0] * coord[0] + coord[1] * coord[1];
			if (d2 > adet * adet)
				b = false;
			if ((coord[2] < hsource / 2 - hsourceup - hdet)
					|| (coord[2] > hsource / 2 - hsourceup))
				b = false;
		} else {
			double d2 = coord[0] * coord[0] + coord[1] * coord[1];
			if (d2 > adet * adet)
				b = false;
		}
		return b;
	}
	
	@Deprecated
	private static int insideDet(double[] coord) {
		int b = 1;// true
		if (isM) {
			double d2 = coord[0] * coord[0] + coord[1] * coord[1];
			if ((coord[2] < hsource / 2 - hsourceup - hdet)
					|| (coord[2] > hsource / 2 - hsourceup))
				b = 0;// false
			if (d2 > adet * adet)
				b = 0;// false
		} else {
			double d2 = coord[0] * coord[0] + coord[1] * coord[1];
			if ((coord[2] < -hsource / 2 - hdet) || (coord[2] > -hsource / 2))
				b = 0;// false
			if (d2 > adet * adet)
				b = 0;// false
		}
		return b;
	}
	
	@Deprecated//Use EGS class instead!!
	// energy in MeV!!!!
	public static void gamaSim(int photon_number, int sourcecode,
			int sourcegeometrycode, int detectorcode, double gammaenergy,
			double det_radius, double source_radius, double det_height,
			double source_height, double source_inner_radius,
			double source_height_up, boolean sourceatt)

	{
		adet = det_radius;// detector radius
		asource = source_radius;// source radius
		hdet = det_height;// detector height
		hsource = source_height;// source height
		// ---------Marinelli-----------------
		hsourceup = source_height_up;// upper source height---MARINELLI
		bsource = source_inner_radius;// inner radius of beaker---MARINELLI
		// --------------------------------------------
		if (sourcegeometrycode == 1)
			isM = true;
		else
			isM = false;

		int npairprod = 0;// score the pair effect
		long startSimulationTime = System.currentTimeMillis();
		simulationTimeElapsed = "";
		int index_low = 0;
		int index_high = 0;
		double ph_interp = 0.0;// interpolation
		double incoh_interp = 0.0;// Compton=incoherent scattering!!
		double incoh_probab = 0.0;// probability
		double coh_probab = 0.0;// probability
		double coh_interp = 0.0;// RAYLEIGH=coherent scattering!!
		double total_interp = 0.0;// interpolation
		double ph_probab = 0.0;// probability
		double r = 0.0;// [0,1] random number
		double r0 = 0.0;// [0,1] random number
		double dist_to_int = 0.0;// store distance to interaction
		double phi = 0.0;// store the azimutal angle
		double phi2 = 0.0;// store the azimutal angle for 2nd annihlilation
							// photon
		double[] newcoord;// new coordinates
		double e_incident = 0.0;// init
		double e_scatt = 0.0;// store the energy of the Compton scattered
								// photons
		double endep = 0.0;// temporary energy deposition of each photon
		double[] oldcrd = new double[6];// old coordinates
		double[] newcoord2 = new double[6];
		;// new coordinates for 2nd annihilation photon
		int nrepeat = 10;// repeat 10 times for unc. eval
		effind = new double[nrepeat];
		double effmed = 0.0;
		double efferr = 0.0;
		double fwhm = 0.0;// error in peak energy due to "electronics"
		// here is the left side of peak at 95% confidence (2*fwhm/2)
		efficiency = new double[2];
		efficiency[0] = 0.0;
		efficiency[1] = 0.0;
		double theta = 0.0;// store the Compton simulation's scattering angle
		double theta2 = 0.0;// store the polar angle for 2nd annihlilation
							// photon
		int n_photons = 0;// number of photons being transported
		// ------read attenuation coefficient and density----------
		// source coefficients
		double sattc_interp = 0.0;
		int scod = sourcecode + 1;// handle further reading from the following
									// general table ([i][0]=first= energy)
		double[][] sourcecoeftable = (double[][]) resources
				.getObject("source.attenuationCoef");
		double[] sourcedensity = (double[]) resources
				.getObject("source.densities");// scod-1 will go
		// the source total linear attenuation coefficient
		double smiu = 0.;
		double smius = sourcedensity[scod - 1];// *mass total
												// ....@@@@@@@@@@@@@@@@@@@@@
		// detector coeficients
		double[][] coeftable = new double[0][0];
		double[] density = (double[]) resources.getObject("detector.densities");
		// detectorcode=0->NaI and 1->Ge
		if (detectorcode == 0) {
			coeftable = (double[][]) resources
					.getObject("detector.nai.attenuationCoef");
			fwhm = (53.0 * Math.sqrt(gammaenergy * 1000 / 661.66)) / 1000;// in
																			// MeV!!
			// System.out.println(""+fwhm);
		} else if (detectorcode == 1) {
			coeftable = (double[][]) resources
					.getObject("detector.ge.attenuationCoef");
			fwhm = (0.6097 + 0.0011 * gammaenergy * 1000) / 1000;// in MeV!!
		} else// anyother option->NaI
		{
			coeftable = (double[][]) resources
					.getObject("detector.nai.attenuationCoef");
			fwhm = (53.0 * Math.sqrt(gammaenergy * 1000 / 661.66)) / 1000;// in
																			// MeV!!
		}
		double dmiu = 0.;
		double dmius = density[detectorcode];// density*mass total
												// ....@@@@@@@@@@@@@@@@@@@@@@
		// coeftable[i][0])-->energy in MeV;
		// coeftable[i][2])-->Compton att. coef in cm2/g;
		// coeftable[i][3])-->Photoelectric abs. coef in cm2/g;etc.
		// ------read data and build arrays-------------------------------------
		double[] energy = new double[coeftable.length];
		double[] coh_coeff = new double[coeftable.length];
		double[] incoh_coeff = new double[coeftable.length];
		double[] ph_coeff = new double[coeftable.length];
		double[] pairn_coeff = new double[coeftable.length];
		double[] paire_coeff = new double[coeftable.length];
		double[] total_coeff = new double[coeftable.length];
		double[] satt_coeff = new double[sourcecoeftable.length];

		for (int i = 0; i < energy.length; i++) {
			energy[i] = coeftable[i][0];
			coh_coeff[i] = coeftable[i][1];
			incoh_coeff[i] = coeftable[i][2];
			ph_coeff[i] = coeftable[i][3];
			pairn_coeff[i] = coeftable[i][4];
			paire_coeff[i] = coeftable[i][5];
			total_coeff[i] = coeftable[i][6];
			satt_coeff[i] = sourcecoeftable[i][scod];
		}
		// ---------------END--------------------------------------------------
		boolean indet = false;// in detector?
		double rcoh = 0.0;// another ALEATOR NUMBER
		double rpair = 0.0;// another ALEATOR NUMBER
		n_photoefect = 0;
		n_pair = 0;
		n_total_all = 0;
		n_fwhm = 0;
		// -------MONTE CARLO
		// SIMULATION---------------------------------------------------
		for (int j = 0; j < nrepeat; j++)// run 10 times for uncertainty
											// evaluation
		{
			// temporary zero init
			n_total = 0;// for each eff eval
			pondt = 0.0;
			// ----------MC BASED LOOP---------------------------
			for (int i = 0; i < photon_number; i++) {
				npairprod = 0;// init
				endep = 0.0;// energy deposited by current photon
				n_total++;
				// n_total_all++;
				e_incident = gammaenergy;
				// -----------getting initial attenuation coeficient
				Sort.findNearestValue(energy, e_incident, true);
				index_low = Sort.getNearestPosition();
				if (index_low < energy.length - 1)
					index_high = index_low + 1;
				else
					index_high = index_low;
				sattc_interp = linInt(energy[index_high],
						satt_coeff[index_high], energy[index_low],
						satt_coeff[index_low], e_incident);
				smiu = smius * sattc_interp;// smiu=smiu*sattc_interp;
				// ---------END------------------------------------
				// ---------getting initial evaluation parameters------
				if (sourcegeometrycode == 0)// cylinder
				{
					oldcrd = getCylinderRandom();
					indet = atSurfaceDet(oldcrd);
				} else if (sourcegeometrycode == 1)// Marrinelli
				{
					oldcrd = getMarrinelliRandom();
					indet = atSurfaceDet(oldcrd);
				} else// default=cylinder
				{
					oldcrd = getCylinderRandom();
					indet = atSurfaceDet(oldcrd);
				}
				if (sourceatt)
					pond = pond * Math.exp(-smiu * source_parcurs);// decreasing
																	// weight
																	// due to
																	// attenuation
																	// in sorce
				// volume. attenuation in air, in source and detector walls etc.
				// are neglected
				// ------------END-----------------------------------
				// test if we are in detector or not
				if (indet) {
					n_total_all++;
					isValidInteraction = 1;// photon reaches the detector
											// surface
				} else
					isValidInteraction = 0;// exit or miss the detector

				n_photons = 1;// 0 only at photoel and exit!!
				while (n_photons == 1) {
					if (stopB)
						return;

					// evaluate ph_probab,ph_interp,compton_interp,total_interp
					// for given energy
					Sort.findNearestValue(energy, e_incident, true);
					index_low = Sort.getNearestPosition();
					if (index_low < energy.length - 1)
						index_high = index_low + 1;
					else
						index_high = index_low;
					// --------------choose interaction coefficients and
					// interaction probabilities----
					ph_interp = linInt(energy[index_high],
							ph_coeff[index_high], energy[index_low],
							ph_coeff[index_low], e_incident);
					incoh_interp = linInt(energy[index_high],
							incoh_coeff[index_high], energy[index_low],
							incoh_coeff[index_low], e_incident);
					total_interp = linInt(energy[index_high],
							total_coeff[index_high], energy[index_low],
							total_coeff[index_low], e_incident);
					coh_interp = linInt(energy[index_high],
							coh_coeff[index_high], energy[index_low],
							coh_coeff[index_low], e_incident);
					incoh_probab = (incoh_interp)
							/ (total_interp - ph_interp - coh_interp);// incoh
																		// from
																		// remaining
					coh_probab = (coh_interp) / (total_interp - ph_interp);// coh
																			// from
																			// REMAINING
					ph_probab = ph_interp / total_interp;// absortion from total
					// ---------------END-----------------------------------------------------------
					dmiu = dmius * total_interp;// dmiu=dmiu*total_interp;

					r = RandomCollection.random01();
					double rdist = RandomCollection.random01();
					while (rdist == 0.0)
						rdist = RandomCollection.random01();
					dist_to_int = -Math.log(rdist) / (dmiu);// [cm]!!@distance
															// to interaction

					if (isValidInteraction != 0)// ---------!!!!!!
					{
						r = RandomCollection.random01();
						if (r <= ph_probab) // photoelectric interaction
						{
							pondt = pondt + pond;// photons are detected by
													// summing in "photopeak "
							n_photoefect++;
							endep = endep + e_incident;
							// So=>absorbed photons->end hystory!!
							if (npairprod == 0)
								n_photons = 0;
							else // is 1 and continue de 2nd tracing
							{
								e_incident = 0.511;
								npairprod--;// score it!
								isValidInteraction = insideDet(newcoord2);
								// permute newcoordonate with old coordonate
								oldcrd[0] = newcoord2[0];
								oldcrd[1] = newcoord2[1];
								oldcrd[2] = newcoord2[2];
								oldcrd[3] = newcoord2[3];
								oldcrd[4] = newcoord2[4];
								oldcrd[5] = newcoord2[5];

							}
						} else // Compton scattering
						{
							// --------RAYLEIGH--elastic COHERENT----------
							rcoh = RandomCollection.random01();
							if (rcoh <= coh_probab) {
								theta = rayleighSim(e_incident);
								r0 = RandomCollection.random01();// [0-1]
								phi = 2 * r0 * Math.PI;// ->[0.2PI];
								newcoord = getUpdatedCoordinate(theta, phi,
										dist_to_int, oldcrd);
								isValidInteraction = insideDet(newcoord);
								// permute newcoordonate with old coordonate
								oldcrd[0] = newcoord[0];
								oldcrd[1] = newcoord[1];
								oldcrd[2] = newcoord[2];
								oldcrd[3] = newcoord[3];
								oldcrd[4] = newcoord[4];
								oldcrd[5] = newcoord[5];
							} else// COMPTON OR PAIR
							{
								if (e_incident < 1.022)// MeV prag!! here is
														// always Compton
								{
									theta = comptonSim(e_incident);// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
									e_scatt = e_incident
											/ (1 + (e_incident / 0.511)
													* (1 - Math.cos(theta)));
									endep = endep + e_incident - e_scatt;
									// continue tracing photon hystory
									e_incident = e_scatt;// new loop in while
									// azimutal angle
									r0 = RandomCollection.random01();// [0-1]
									phi = 2 * r0 * Math.PI;// ->[0.2PI];
									newcoord = getUpdatedCoordinate(theta, phi,
											dist_to_int, oldcrd);
									isValidInteraction = insideDet(newcoord);
									// permute newcoordonate with old coordonate
									oldcrd[0] = newcoord[0];
									oldcrd[1] = newcoord[1];
									oldcrd[2] = newcoord[2];
									oldcrd[3] = newcoord[3];
									oldcrd[4] = newcoord[4];
									oldcrd[5] = newcoord[5];
								} else // pair production could be
								{
									rpair = RandomCollection.random01();
									if (rpair <= incoh_probab) // Compton
									{
										theta = comptonSim(e_incident);// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
										e_scatt = e_incident
												/ (1 + (e_incident / 0.511)
														* (1 - Math.cos(theta)));
										endep = endep + e_incident - e_scatt;
										// continue tracing photon hystory
										e_incident = e_scatt;// new loop in
																// while
										// azimutal angle
										r0 = RandomCollection.random01();// [0-1]
										phi = 2 * r0 * Math.PI;// ->[0.2PI];
										newcoord = getUpdatedCoordinate(theta,
												phi, dist_to_int, oldcrd);
										isValidInteraction = insideDet(newcoord);
										// permute newcoordonate with old
										// coordonate
										oldcrd[0] = newcoord[0];
										oldcrd[1] = newcoord[1];
										oldcrd[2] = newcoord[2];
										oldcrd[3] = newcoord[3];
										oldcrd[4] = newcoord[4];
										oldcrd[5] = newcoord[5];
									} else// pair production
									{
										n_pair++;
										npairprod = 2;// score it
										endep = endep + e_incident - 1.022;
										r0 = RandomCollection.random01();// [0-1]
										theta = Math.PI * r0;
										theta2 = Math.PI - theta;
										r0 = RandomCollection.random01();// [0-1]
										phi = 2 * r0 * Math.PI;// ->[0.2PI];
										phi2 = Math.PI + phi;
										e_incident = 0.511;
										newcoord = getUpdatedCoordinate(theta,
												phi, dist_to_int, oldcrd);
										newcoord2 = getUpdatedCoordinate(
												theta2, phi2, dist_to_int,
												oldcrd);
										// first photon
										npairprod--;// score it!
										isValidInteraction = insideDet(newcoord);
										// permute newcoordonate with old
										// coordonate
										oldcrd[0] = newcoord[0];
										oldcrd[1] = newcoord[1];
										oldcrd[2] = newcoord[2];
										oldcrd[3] = newcoord[3];
										oldcrd[4] = newcoord[4];
										oldcrd[5] = newcoord[5];
									}

								}
							}
						}
					} else {
						if (npairprod == 0) {
							n_photons = 0;
							if (endep >= gammaenergy - fwhm) {
								n_fwhm++;
								pondt = pondt + pond;
							}
						} else // is 1 and continue de 2nd tracing
						{
							e_incident = 0.511;
							npairprod--;// score it!
							isValidInteraction = insideDet(newcoord2);
							// permute newcoordonate with old coordonate
							oldcrd[0] = newcoord2[0];
							oldcrd[1] = newcoord2[1];
							oldcrd[2] = newcoord2[2];
							oldcrd[3] = newcoord2[3];
							oldcrd[4] = newcoord2[4];
							oldcrd[5] = newcoord2[5];
						}
						// pathlength -->exit the phantom
						// n_photons=0;//exit region of interest, so we kill the
						// photon
						// taking into account the theoretical resolution
					}
				}// end while
			}// end for i=0 to photon_number
			effind[j] = pondt / n_total;
		}// end energy loop
		// calculate mean and variance
		Stats.avevar(effind, effind.length);
		effmed = Stats.ave_avevar;// Stat.mean(effind);
		efferr = Stats.var_avevar;// Stat.stDevOfMean(effind);
		efferr = Math.sqrt(efferr / effind.length);
		int glib = nrepeat - 1;// degree of freedom
		double t = StatsUtil.getStudentFactor(glib);// Stat.tTestCoefForFinal(glib,
													// 95);//95% confidence
													// level
		efferr = efferr * t;// unc for 95% confidence level
		efficiency[0] = effmed;
		efficiency[1] = efferr;

		simulationTimeElapsed = timeElapsed(startSimulationTime);
	}
	
	@Deprecated//Use EGS class instead!!
	// ---------------------------------------------------------------------
	public static void gamaSim(int photon_number, int sourcecode,
			int sourcegeometrycode, int detectorcode, double gammaenergy,
			double det_radius, double source_radius, double det_height,
			double source_height, double source_inner_radius,
			double source_height_up,
			// ////////////////////////////////////
			double deltazup, double deltazdown, double deltaxy,
			String[] winDetType, double[] winThickness,
			// ////////////////////////////////////
			boolean sourceatt)

	{
		geomB = false;
		adet = det_radius;// detector radius
		asource = source_radius;// source radius
		hdet = det_height;// detector height
		hsource = source_height;// source height
		// ---------Marinelli-----------------
		hsourceup = source_height_up;// upper source height---MARINELLI
		bsource = source_inner_radius;// inner radius of beaker---MARINELLI
		// --------------------------------------------
		if (sourcegeometrycode == 1)
			isM = true;
		else
			isM = false;

		int npairprod = 0;// score the pair effect
		long startSimulationTime = System.currentTimeMillis();
		simulationTimeElapsed = "";
		int index_low = 0;
		int index_high = 0;
		double ph_interp = 0.0;// interpolation
		double incoh_interp = 0.0;// Compton=incoherent scattering!!
		double incoh_probab = 0.0;// probability
		double coh_probab = 0.0;// probability
		double coh_interp = 0.0;// RAYLEIGH=coherent scattering!!
		double total_interp = 0.0;// interpolation
		double ph_probab = 0.0;// probability
		double r = 0.0;// [0,1] random number
		double r0 = 0.0;// [0,1] random number
		double dist_to_int = 0.0;// store distance to interaction
		double phi = 0.0;// store the azimutal angle
		double phi2 = 0.0;// store the azimutal angle for 2nd annihlilation
							// photon
		double[] newcoord;// new coordinates
		double e_incident = 0.0;// init
		double e_scatt = 0.0;// store the energy of the Compton scattered
								// photons
		double endep = 0.0;// temporary energy deposition of each photon
		double[] oldcrd = new double[6];// old coordinates
		double[] newcoord2 = new double[6];
		;// new coordinates for 2nd annihilation photon
		int nrepeat = 10;// repeat 10 times for unc. eval
		effind = new double[nrepeat];
		double effmed = 0.0;
		double efferr = 0.0;
		double fwhm = 0.0;// error in peak energy due to "electronics"
		// here is the left side of peak at 95% confidence (2*fwhm/2)
		efficiency = new double[2];
		efficiency[0] = 0.0;
		efficiency[1] = 0.0;
		double theta = 0.0;// store the Compton simulation's scattering angle
		double theta2 = 0.0;// store the polar angle for 2nd annihlilation
							// photon
		int n_photons = 0;// number of photons being transported
		// ------read attenuation coefficient and density----------
		// source coefficients
		double sattc_interp = 0.0;
		int scod = sourcecode + 1;// handle further reading from the following
									// general table ([i][0]=first= energy)
		double[][] sourcecoeftable = (double[][]) resources
				.getObject("source.attenuationCoef");
		double[] sourcedensity = (double[]) resources
				.getObject("source.densities");// scod-1 will go
		// the source total linear attenuation coefficient
		double smiu = 0.;
		double smius = sourcedensity[scod - 1];// *mass total
												// ....@@@@@@@@@@@@@@@@@@@@@
		// detector coeficients
		double[][] coeftable = new double[0][0];
		double[] density = (double[]) resources.getObject("detector.densities");
		// detectorcode=0->NaI and 1->Ge
		if (detectorcode == 0) {
			coeftable = (double[][]) resources
					.getObject("detector.nai.attenuationCoef");
			fwhm = (53.0 * Math.sqrt(gammaenergy * 1000 / 661.66)) / 1000;// in
																			// MeV!!
			// System.out.println(""+fwhm);
		} else if (detectorcode == 1) {
			coeftable = (double[][]) resources
					.getObject("detector.ge.attenuationCoef");
			fwhm = (0.6097 + 0.0011 * gammaenergy * 1000) / 1000;// in MeV!!
		} else// anyother option->NaI
		{
			coeftable = (double[][]) resources
					.getObject("detector.nai.attenuationCoef");
			fwhm = (53.0 * Math.sqrt(gammaenergy * 1000 / 661.66)) / 1000;// in
																			// MeV!!
		}
		double dmiu = 0.;
		double dmius = density[detectorcode];// density*mass total
												// ....@@@@@@@@@@@@@@@@@@@@@@
		// coeftable[i][0])-->energy in MeV;
		// coeftable[i][2])-->Compton att. coef in cm2/g;
		// coeftable[i][3])-->Photoelectric abs. coef in cm2/g;etc.
		// ------read data and build arrays-------------------------------------
		double[] energy = new double[coeftable.length];
		double[] coh_coeff = new double[coeftable.length];
		double[] incoh_coeff = new double[coeftable.length];
		double[] ph_coeff = new double[coeftable.length];
		double[] pairn_coeff = new double[coeftable.length];
		double[] paire_coeff = new double[coeftable.length];
		double[] total_coeff = new double[coeftable.length];
		double[] satt_coeff = new double[sourcecoeftable.length];

		for (int i = 0; i < energy.length; i++) {
			energy[i] = coeftable[i][0];
			coh_coeff[i] = coeftable[i][1];
			incoh_coeff[i] = coeftable[i][2];
			ph_coeff[i] = coeftable[i][3];
			pairn_coeff[i] = coeftable[i][4];
			paire_coeff[i] = coeftable[i][5];
			total_coeff[i] = coeftable[i][6];
			satt_coeff[i] = sourcecoeftable[i][scod];
		}
		// ---------------END--------------------------------------------------
		boolean indet = false;// in detector?
		double rcoh = 0.0;// another ALEATOR NUMBER
		double rpair = 0.0;// another ALEATOR NUMBER
		n_photoefect = 0;
		n_pair = 0;
		n_total_all = 0;
		n_fwhm = 0;
		// -------MONTE CARLO
		// SIMULATION---------------------------------------------------
		for (int j = 0; j < nrepeat; j++)// run 10 times for uncertainty
											// evaluation
		{
			// temporary zero init
			n_total = 0;// for each eff eval
			pondt = 0.0;
			// ----------MC BASED LOOP---------------------------
			for (int i = 0; i < photon_number; i++) {
				npairprod = 0;// init
				endep = 0.0;// energy deposited by current photon
				n_total++;
				// n_total_all++;
				e_incident = gammaenergy;
				// -----------getting initial attenuation coeficient
				Sort.findNearestValue(energy, e_incident, true);
				index_low = Sort.getNearestPosition();
				if (index_low < energy.length - 1)
					index_high = index_low + 1;
				else
					index_high = index_low;
				sattc_interp = linInt(energy[index_high],
						satt_coeff[index_high], energy[index_low],
						satt_coeff[index_low], e_incident);
				smiu = smius * sattc_interp;// smiu=smiu*sattc_interp;
				// ---------END------------------------------------
				// ---------getting initial evaluation parameters------
				if (sourcegeometrycode == 0)// cylinder
				{
					oldcrd = getCylinderRandom2(deltazup, winDetType,
							winThickness, e_incident);
					indet = atSurfaceDet(oldcrd);
				} else if (sourcegeometrycode == 1)// Marrinelli
				{
					oldcrd = getMarrinelliRandom2(deltazup, deltazdown,
							deltaxy, winDetType, winThickness, e_incident);
					indet = atSurfaceDet(oldcrd);
				} else// default=cylinder
				{
					oldcrd = getCylinderRandom2(deltazup, winDetType,
							winThickness, e_incident);
					indet = atSurfaceDet(oldcrd);
				}
				if (sourceatt)
					pond = pond * Math.exp(-smiu * source_parcurs);// decreasing
																	// weight
																	// due to
																	// attenuation
																	// in sorce
				// volume. attenuation in air, in source and detector walls etc.
				// are neglected
				// ------------END-----------------------------------
				// test if we are in detector or not
				if (indet) {
					n_total_all++;
					isValidInteraction = 1;// photon reaches the detector
											// surface
				} else
					isValidInteraction = 0;// exit or miss the detector

				n_photons = 1;// 0 only at photoel and exit!!
				while (n_photons == 1) {
					if (stopB)
						return;

					// evaluate ph_probab,ph_interp,compton_interp,total_interp
					// for given energy
					Sort.findNearestValue(energy, e_incident, true);
					index_low = Sort.getNearestPosition();
					if (index_low < energy.length - 1)
						index_high = index_low + 1;
					else
						index_high = index_low;
					// --------------choose interaction coefficients and
					// interaction probabilities----
					ph_interp = linInt(energy[index_high],
							ph_coeff[index_high], energy[index_low],
							ph_coeff[index_low], e_incident);
					incoh_interp = linInt(energy[index_high],
							incoh_coeff[index_high], energy[index_low],
							incoh_coeff[index_low], e_incident);
					total_interp = linInt(energy[index_high],
							total_coeff[index_high], energy[index_low],
							total_coeff[index_low], e_incident);
					coh_interp = linInt(energy[index_high],
							coh_coeff[index_high], energy[index_low],
							coh_coeff[index_low], e_incident);
					incoh_probab = (incoh_interp)
							/ (total_interp - ph_interp - coh_interp);// incoh
																		// from
																		// remaining
					coh_probab = (coh_interp) / (total_interp - ph_interp);// coh
																			// from
																			// REMAINING
					ph_probab = ph_interp / total_interp;// absortion from total
					// ---------------END-----------------------------------------------------------
					dmiu = dmius * total_interp;// dmiu=dmiu*total_interp;

					r = RandomCollection.random01();
					double rdist = RandomCollection.random01();
					while (rdist == 0.0)
						rdist = RandomCollection.random01();
					dist_to_int = -Math.log(rdist) / (dmiu);// [cm]!!@distance
															// to interaction

					if (isValidInteraction != 0)// ---------!!!!!!
					{
						r = RandomCollection.random01();
						if (r <= ph_probab) // photoelectric interaction
						{
							pondt = pondt + pond;// photons are detected by
													// summing in "photopeak "
							n_photoefect++;
							endep = endep + e_incident;
							// So=>absorbed photons->end hystory!!
							if (npairprod == 0)
								n_photons = 0;
							else // is 1 and continue de 2nd tracing
							{
								e_incident = 0.511;
								npairprod--;// score it!
								isValidInteraction = insideDet(newcoord2);
								// permute newcoordonate with old coordonate
								oldcrd[0] = newcoord2[0];
								oldcrd[1] = newcoord2[1];
								oldcrd[2] = newcoord2[2];
								oldcrd[3] = newcoord2[3];
								oldcrd[4] = newcoord2[4];
								oldcrd[5] = newcoord2[5];

							}
						} else // Compton scattering
						{
							// --------RAYLEIGH--elastic COHERENT----------
							rcoh = RandomCollection.random01();
							if (rcoh <= coh_probab) {
								theta = rayleighSim(e_incident);
								r0 = RandomCollection.random01();// [0-1]
								phi = 2 * r0 * Math.PI;// ->[0.2PI];
								newcoord = getUpdatedCoordinate(theta, phi,
										dist_to_int, oldcrd);
								isValidInteraction = insideDet(newcoord);
								// permute newcoordonate with old coordonate
								oldcrd[0] = newcoord[0];
								oldcrd[1] = newcoord[1];
								oldcrd[2] = newcoord[2];
								oldcrd[3] = newcoord[3];
								oldcrd[4] = newcoord[4];
								oldcrd[5] = newcoord[5];
							} else// COMPTON OR PAIR
							{
								if (e_incident < 1.022)// MeV prag!! here is
														// always Compton
								{
									theta = comptonSim(e_incident);// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
									e_scatt = e_incident
											/ (1 + (e_incident / 0.511)
													* (1 - Math.cos(theta)));
									endep = endep + e_incident - e_scatt;
									// continue tracing photon hystory
									e_incident = e_scatt;// new loop in while
									// azimutal angle
									r0 = RandomCollection.random01();// [0-1]
									phi = 2 * r0 * Math.PI;// ->[0.2PI];
									newcoord = getUpdatedCoordinate(theta, phi,
											dist_to_int, oldcrd);
									isValidInteraction = insideDet(newcoord);
									// permute newcoordonate with old coordonate
									oldcrd[0] = newcoord[0];
									oldcrd[1] = newcoord[1];
									oldcrd[2] = newcoord[2];
									oldcrd[3] = newcoord[3];
									oldcrd[4] = newcoord[4];
									oldcrd[5] = newcoord[5];
								} else // pair production could be
								{
									rpair = RandomCollection.random01();
									if (rpair <= incoh_probab) // Compton
									{
										theta = comptonSim(e_incident);// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
										e_scatt = e_incident
												/ (1 + (e_incident / 0.511)
														* (1 - Math.cos(theta)));
										endep = endep + e_incident - e_scatt;
										// continue tracing photon hystory
										e_incident = e_scatt;// new loop in
																// while
										// azimutal angle
										r0 = RandomCollection.random01();// [0-1]
										phi = 2 * r0 * Math.PI;// ->[0.2PI];
										newcoord = getUpdatedCoordinate(theta,
												phi, dist_to_int, oldcrd);
										isValidInteraction = insideDet(newcoord);
										// permute newcoordonate with old
										// coordonate
										oldcrd[0] = newcoord[0];
										oldcrd[1] = newcoord[1];
										oldcrd[2] = newcoord[2];
										oldcrd[3] = newcoord[3];
										oldcrd[4] = newcoord[4];
										oldcrd[5] = newcoord[5];
									} else// pair production
									{
										n_pair++;
										npairprod = 2;// score it
										endep = endep + e_incident - 1.022;
										r0 = RandomCollection.random01();// [0-1]
										theta = Math.PI * r0;
										theta2 = Math.PI - theta;
										r0 = RandomCollection.random01();// [0-1]
										phi = 2 * r0 * Math.PI;// ->[0.2PI];
										phi2 = Math.PI + phi;
										e_incident = 0.511;
										newcoord = getUpdatedCoordinate(theta,
												phi, dist_to_int, oldcrd);
										newcoord2 = getUpdatedCoordinate(
												theta2, phi2, dist_to_int,
												oldcrd);
										// first photon
										npairprod--;// score it!
										isValidInteraction = insideDet(newcoord);
										// permute newcoordonate with old
										// coordonate
										oldcrd[0] = newcoord[0];
										oldcrd[1] = newcoord[1];
										oldcrd[2] = newcoord[2];
										oldcrd[3] = newcoord[3];
										oldcrd[4] = newcoord[4];
										oldcrd[5] = newcoord[5];
									}

								}
							}
						}
					} else {
						if (npairprod == 0) {
							n_photons = 0;
							if (endep >= gammaenergy - fwhm) {
								n_fwhm++;
								pondt = pondt + pond;
							}
						} else // is 1 and continue de 2nd tracing
						{
							e_incident = 0.511;
							npairprod--;// score it!
							isValidInteraction = insideDet(newcoord2);
							// permute newcoordonate with old coordonate
							oldcrd[0] = newcoord2[0];
							oldcrd[1] = newcoord2[1];
							oldcrd[2] = newcoord2[2];
							oldcrd[3] = newcoord2[3];
							oldcrd[4] = newcoord2[4];
							oldcrd[5] = newcoord2[5];
						}
						// pathlength -->exit the phantom
						// n_photons=0;//exit region of interest, so we kill the
						// photon
						// taking into account the theoretical resolution
					}
				}// end while
			}// end for i=0 to photon_number
			effind[j] = pondt / n_total;
		}// end energy loop
		// calculate mean and variance
		Stats.avevar(effind, effind.length);
		effmed = Stats.ave_avevar;// Stat.mean(effind);
		efferr = Stats.var_avevar;// Stat.stDevOfMean(effind);
		efferr = Math.sqrt(efferr / effind.length);
		int glib = nrepeat - 1;// degree of freedom
		double t = StatsUtil.getStudentFactor(glib);// Stat.tTestCoefForFinal(glib,
													// 95);//95% confidence
													// level
		efferr = efferr * t;// unc for 95% confidence level
		efficiency[0] = effmed;
		efficiency[1] = efferr;
		System.out.println("resph : " + n_fwhm);
		simulationTimeElapsed = timeElapsed(startSimulationTime);
	}

	// ---------------------------------------------------------------------

}
