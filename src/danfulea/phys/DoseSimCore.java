package danfulea.phys;

import java.util.ResourceBundle;

import danfulea.math.RandomCollection;
import danfulea.math.Sort;

/**
 * Class for Monte Carlo simulation of radiation transport which is used in conjunction with MIRD5 Phantom 
 * class and XRay class to assess organ doses and effective dose for patients undergoing 
 * radiological examinations (radiodiagnostic - mammography, radiography and CT). This simulation 
 * uses KERMA approximation, i.e. only the photons are tracked and the dose is considered deposited wherever 
 * the photons are absorbed. This approximation is good enough for human body due to the fact that 
 * organs are large compared to the electron path-length. However, at boundaries this approach is not 
 * good and therefore a much more accurate method is needed. Such methods must take into account 
 * photons and electrons. EGSnrc and GEANT4 are both very good MC engines suitable for this task. 
 * As far as I know the EGSnrc is very good for simple geometries such as cylindrical (RZ) or spherical. 
 * Since MIRD5 phantom is anything but simple, EGSnrc cannot be used unless we develop complex geometry routines which is 
 * a pain. The CERN's GEANT4 on the other hand 
 * has build-in very good geometry routines which can easily handle all kind of situations and is based on the latest 
 * top notch data. GEANT4 as well as EGSnrc are in fact toolkits, so users must write their own routines 
 * to couple with those very good MC engines. GEANT4 requires C++ programming language, 
 * EGSnrc requires MORTRAN (a variant of fortran) and recently it was ported to C++.<br> 
 * Bottom line: Use this class only for gross (rough) estimations of doses and/or 
 * in educational purpose. For more accuracy use GEANT4. Personally I used this class inside a computer program  
 * called IradMed to perform various quality control tests and dose/risk assessments in my PhD thesis. Nowadays, 
 * I always use GEANT4 for such tasks. 
 * 
 * 
 * @author Dan Fulea, 10 APR. 2005
 */
public class DoseSimCore {
	private static final String BASE_RESOURCE_CLASS = "danfulea.phys.resources.DoseSimCoreResources";
	private static ResourceBundle resources = ResourceBundle
			.getBundle(BASE_RESOURCE_CLASS);
	private static String simulationTimeElapsed = "";// store simulation time
														// elapsed
	private static double total_energy_deposit = 0.0;// stores the total energy
														// deposition (MeV)
	private static double[][] energy_deposit;// stores the average energy
												// deposition per organ(MeV)
	// and its associated error (Poisson 2 sigma in %)
	private static int n_detected = 0;// stores the number of detected photons
										// (exit the breast)
	private static int n_total = 0;// stores the total number of simulated
									// photons
	private static int n_film = 0;// stores the number of photons which exit
									// towards film-side
	private static int n_film_direct = 0;// stores the number of photons which
											// exit towards film-side with no
											// interactions!
	private static int n_side = 0;// stores the number of photons which exit
									// around the phantom!
	private static int n_photoefect = 0;// stores the total number of
										// photoelectric absorbed photons
	private static double[][] organAbsorbedDose;// organ Doses in mGy!!!
	private static double[] effectiveDose;// effective Doses in mSv!!!
	private static String[] organs;// organe de interes!
	// ---------------------------------------------------CT------------------------
	private static double[][] orgAbsorbedDose;// organ Doses in mGy!!!
	private static int scannedslices;
	// ------------------------------------------------------------------------------
	private static boolean stopB = false;

	// ------------------------------------------------------
	// -------------------------------------------------------
	private static int knAlgor = 0;// Wielopolski

	/**
	 * Setup Klein-Nishina algorithm for Compton simulation routine.
	 * @param kn kn, the algoritm index
	 */
	public static void setKNAlgor(int kn) {
		knAlgor = kn;
	}

	// -------------------------------------------------------

	// ------------------------------------------------------

	/**
	 * Sets flag to stop simulation
	 */
	public static void setStop() {
		stopB = true;
	}

	/**
	 * Sets flag to start simulation
	 */
	public static void setStart() {
		stopB = false;
	}

	/**
	 * Get energy deposition matrix. Columns are energy and its uncertainty. Rows are 
	 * related to each organ location.
	 * @return the result
	 */
	public static double[][] getEnergyDeposition() {
		return energy_deposit;
	}

	// ---------------------CT------------------------------------
	/**
	 * Get organ absorbed dose for CT. Columns are dose and its uncertainty. Rows are 
	 * related to each organ location.
	 * @return the result
	 */
	public static double[][] getctOrganAbsorbedDose() {
		return orgAbsorbedDose;
	}

	/**
	 * Get the number of slices used in CT scan.
	 * @return the result
	 */
	public static int getctScannedSlices() {
		return scannedslices;
	}

	// ----------------------------------------------------------
	/**
	 * Get organ absorbed dose. Columns are dose and its uncertainty. Rows are 
	 * related to each organ location.
	 * @return the result
	 */
	public static double[][] getOrganAbsorbedDose() {
		return organAbsorbedDose;
	}

	/**
	 * Get the array of organs.
	 * @return the result
	 */
	public static String[] getOrgans() {
		return organs;
	}

	/**
	 * Get the effective dose [0] and its uncertainty [1].
	 * @return the result array
	 */
	public static double[] getEffectiveDose() {
		return effectiveDose;
	}

	/**
	 * Get total photons.
	 * @return the result.
	 */
	public static int getTotalPhotons() {
		return n_total;
	}

	/**
	 * Return the number of photons which hit the "detector" (phantom).
	 * @return the result
	 */
	public static int getDetectedPhotons() {
		return n_detected;
	}

	/**
	 * Returns the number of photons which reached the film (exit through the back of the phantom)  
	 * @return the result
	 */
	public static int getFilmPhotons() {
		return n_film;
	}

	/**
	 * Returns the number of photons which undergo photo-electric effect.
	 * @return the result
	 */
	public static int getPhotoAbsorbedPhotons() {
		return n_photoefect;
	}

	/**
	 * 
	 * @return the time taken to complete this simulation.
	 */
	public static String getSimulationTimeElapsed() {
		return simulationTimeElapsed;
	}

	/**
	 * Percentage of photons escaped through any side of phantom.
	 * @return the result
	 */
	public static double getSideEscapedProc() {
		if (n_total == 0)
			return 0.0;

		return 100. * n_side / n_total;
	}

	/**
	 * Percentage of photons escaped through phantom back.
	 * @return the result
	 */
	public static double getFilmEscapedProc() {
		if (n_total == 0)
			return 0.0;

		return 100. * n_film / n_total;
	}

	/**
	 * Percentage of photons which exit through phantom back but with no interactions 
	 * along their way through phantom relative to the photons which exit through phanotm back.
	 * @return the result
	 */
	public static double getFilmScatteredVsTotalFilmProc() {
		if (n_film == 0)
			return 0.0;

		return 100. * (n_film - n_film_direct) / n_film;
	}

	// ---------------------------------------------------------------------------------------------
	/**
	 * Dispatch the Compton simulation process based on the chosen Klein-Nishina algorithm.
	 * @param k the photon energy in MeV
	 * @return the polar angle of Compton scattering.
	 */
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
	 * Compton simulation process based on EGSnrc algorithm.
	 * @param k the photon energy in MeV
	 * @return the polar angle of Compton scattering.
	 */
	public static double comptonSimE(double k)// EGS4
	{
		//int n_acc = 0; // number of accepted Compton events
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

	// based on Wielopolski algorithm (1987)--O.Sima book
	/**
	 * Compton simulation process based on Wielopolski algorithm.
	 * @param k the photon energy in MeV
	 * @return the polar angle of Compton scattering.
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
	 * Compton simulation process based on Khan algorithm.
	 * @param k the photon energy in MeV
	 * @return the polar angle of Compton scattering.
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
	// -------------------------------------------------------------------------------------
	// function to simulate Compton scattering.
	// input:energy of photon
	// output: scattering angle of photon in radians
	// based on code written by Lesley Buckley Dec17,2001
	/**
	 * Compton simulation process based on "classic" algorithm. 
	 * Based on code written by Lesley Buckley Dec 17, 2001 in his/her PhD thesis.
	 * @param k the photon energy in MeV
	 * @return the polar angle of Compton scattering.
	 */
	public static double comptonSimC(double k) {
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

	// ------------------
	/**
	 * Rayleigh simulation process (coherent scattering) crude approach. 
	 * 
	 * @param k the photon energy in MeV, not used in this approach but in later development (if any)
	 * @return the polar angle of Rayleigh scattering.
	 */
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
	/**
	 * Format the time elapsed from a start time. 
	 * @param startTime the start time
	 * @return the formatted String representation of time elapsed.
	 */
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
	/**
	 * Perform linear interpolation between point (x1,y1) and (x2,y2) and return y for a given x.
	 * @param x1 x1
	 * @param y1 y1
	 * @param x2 x2
	 * @param y2 y2
	 * @param x x
	 * @return the result
	 */
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
	/**
	 * Being the core for transporting particle, this routine updates particle coordinates based on the old coordinates, the distance 
	 * to travel, the polar and azimuth angles. The term coordinate refers to a 6 dimensional array 
	 * containing the spatial position relative to the phantom coordinate system x,y,z and the direction 
	 * to travel described by direction cosines ux, uy and uz (or u, v, w). <br>The scenario is as follows: particle is located 
	 * at point A having coordinates cA. At this point, the Monte Carlo code figured out what kind of interaction the particle has, 
	 * e.g Compton scattering, and figured out the polar and azimuth angle related to that scattering process. 
	 * From interaction sampling we figured out what is the distance to the next interaction site (point B). Finally, transport particle to B  
	 * (i.e. update coordinates to cB). This process is then repeated (determine what kind of scattering (if any) happen at B and so on) 
	 * until particle escapes the phantom or it is absorbed by photo-electric effect.  
	 * @param theta0 polar angle of scattering
	 * @param phi0 azimuth angle of scattering
	 * @param dist distance to the next interaction location (radiation interaction that is).
	 * @param oldcoord current (and old) coordinates
	 * @return the new coordinates where particle is transported.
	 */
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

	// general spectrum--xrs, focus-skin distance
	// exposure =free in air without backscatter im miliGray
	/**
	 * Perform Monte Carlo simulation for mammography.
	 * @param photon_number the number of histories for each photon energy.
	 * @param xrs the link to XRaySpectrum class.
	 * @param airExposure the dose free in air at phantom (breast) entrance in mGy.
	 * @param fsd focus-breast entrance distance in cm (focus-top cylinder base, breast is simulated 
	 * as an right cylinder due to the compression device).
	 */
	public static void mamoSim(int photon_number, XRaySpectrum xrs,
			double airExposure, double fsd) {
		// -store the index of energy array corresponding to the lower than
		// e_incident
		// and higher than e_incident
		long startSimulationTime = System.currentTimeMillis();
		simulationTimeElapsed = "";
		int index_low = 0;
		int index_high = 0;
		//double e_low = 0.0;
		//double e_high = 0.0;
		double ph_interp = 0.0;
		//double incoh_interp = 0.0;// Compton=incoherent scattering!!
		double total_interp = 0.0;
		double ph_probab = 0.0;
		double r = 0.0;// [0,1] random number
		double dist_to_int = 0.0;// store distance to interaction
		Phantom phantom = new Phantom(Phantom.MAMO_INDEX);
		double breastThickness = phantom.getOrganThickness();// organism
		// -------------SC-------------------------------------------------------------
		double radius = phantom.getMaximumOrganDimension();// organism
		double theta0_max = 0.0;// @@@@@@@@@@@@@@@@@@@@@@@@@
		// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
		double costmax = fsd / Math.sqrt(fsd * fsd + radius * radius);
		// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@
		if (fsd > 0.0)
			theta0_max = Math.abs(Math.atan(radius / fsd));// idem
															// asin:-pi/2-pi/2,
															// acos:0,pi
		else
			theta0_max = 0.0;
		double r0 = 0.0;
		double phi0 = 0.0;
		double theta0 = 0.0;
		//double ux0 = 0.0;
		//double uy0 = 0.0;
		//double uz0 = 0.0;
		double[] oldcoord = new double[6];
		double phi = 0.0;
		double[] newcoord;
		int isValidInteraction = Phantom.REMAINDER;// initial enter in organism
		// ----------------END
		// SC-----------------------------------------------------
		double breastTransversalArea = phantom.getEntranceTransversalArea() * 10 * 10;// -for
																						// Area
																						// in
																						// mm2
		// -----------------------//organism
		double e_incident = 0.0;// init
		double e_scatt = 0.0;// store the energy of the Compton scattered
								// photons
		// -----------INIT!!---------------------------
		total_energy_deposit = 0.0;
		energy_deposit = new double[1][2];// only breast here
		energy_deposit[0][0] = 0.0;
		energy_deposit[0][1] = 0.0;
		organAbsorbedDose = new double[1][2];// only breast here
		organAbsorbedDose[0][0] = 0.0;
		organAbsorbedDose[0][1] = 0.0;
		// ---------------------------------
		//double enmedtemp = 0.0;// handle error eval
		double enindivtemp = 0.0;// handle error eval

		@SuppressWarnings("unused")
		int n_interaction = 0;
		double n_interaction_total = 0.0;
		//double dummy = 0.0;
		// ----------------------------------
		effectiveDose = new double[2];
		effectiveDose[0] = 0.0;
		effectiveDose[1] = 0.0;
		n_detected = 0;
		n_total = 0;
		n_film = 0;
		n_film_direct = 0;
		n_side = 0;
		n_photoefect = 0;
		double theta = 0.0;// store the Compton simulation's scattering angle
		int n_photons = 0;// number of photons being transported
		// ------read attenuation coefficient and density----------
		double[][] coeftable = (double[][]) resources
				.getObject("breast.diagnostic.attenuationCoef");
		// coeftable[i][0])-->energy in MeV;
		// coeftable[i][2])-->Compton att. coef in cm2/g;
		// coeftable[i][3])-->Photoelectric abs. coef in cm2/g;
		// coeftable[i][5])-->Total (without Rayleigh scattering) att. coef in
		// cm2/g;
		// read data:
		double[] energy = new double[coeftable.length];
		double[] ph_coeff = new double[coeftable.length];
		double[] incoh_coeff = new double[coeftable.length];
		double[] total_coeff = new double[coeftable.length];
		// ----------------------------------------------------------------------
		double[] coh_coeff = new double[coeftable.length];
		double[] total_coeff_r = new double[coeftable.length];// @@@@@@@RAYLEIGH

		for (int i = 0; i < energy.length; i++) {
			energy[i] = coeftable[i][0];
			incoh_coeff[i] = coeftable[i][2];
			ph_coeff[i] = coeftable[i][3];
			total_coeff[i] = coeftable[i][5];
			// --------------------------------------------------------@@@@@@@RAYLEIGH
			coh_coeff[i] = coeftable[i][1];
			total_coeff_r[i] = coeftable[i][4];
		}
		double coh_probab = 0.0;
		double coh_interp = 0.0;// RAYLEIGH=coherent scattering!!
		double total_interp_r = 0.0;
		double rcoh = 0.0;// ALEATOR NUMBER
		// ------------------------------------@@@@@@@@@@@@@@@@@@@@@@@@RAYLEIGH
		double density = ((Double) resources
				.getObject("breast.diagnostic.density")).doubleValue();
		double breastMass = breastTransversalArea * breastThickness * density
				/ 100;
		organs = (String[]) resources.getObject("mamo.organname");

		double endep = 0.0;
		double totalendep = 0.0;
		double endeperror = 0.0;
		// -------MONTE CARLO
		// SIMULATION---------------------------------------------------
		boolean interactB = false;
		for (int j = 0; j < (xrs.getXRayEnergies()).length; j++) {
			// for each energy
			endep = 0.0;
			totalendep = 0.0;
			endeperror = 0.0;
			n_interaction = 0;
			// ------------------
			for (int i = 0; i < photon_number; i++) {
				interactB = false;
				n_total++;
				e_incident = (xrs.getXRayEnergies())[j] / 1000;// incident_energy
																// in MeV;
				// ----------------------------------------------SC--------------------------------
				// generating an aleator theta0 and a polar angle phi!!!
				isValidInteraction = Phantom.REMAINDER;// incident flux-->always
														// true//----------------------
				r0 = RandomCollection.random01();// [0-1]
				phi0 = 2 * r0 * Math.PI;
				// for theta0->[0,theta0_max]
				theta0 = theta0_max * RandomCollection.random01();
				// @@@@@@@@@@@@@@
				// costet evaluation-->polar angle
				double dom = (1.0 - costmax) / 2.0;// >0 and <1/2, costmax <1
													// and >0, thetamax<90!
				r = RandomCollection.random01();// EGS4.random01();
				r = r * dom;
				double costet = 1.0 - 2.0 * r;// 2*r-1;//<0 always--- negativ z
												// axis!!
				double sintet = Math.sqrt(1.0 - costet * costet);
				double tgtet = sintet / costet;
				// r=RandomCollection.random01();//EGS4.random01();
				theta0 = Math.abs(Math.atan(tgtet));
				// @@@@@@@@@@@@@@
				oldcoord[2] = 0.0;// organism or breast in this case->z
				oldcoord[0] = fsd * Math.tan(theta0) * Math.cos(phi0);// ->x
				oldcoord[1] = fsd * Math.tan(theta0) * Math.sin(phi0);// ->y
				oldcoord[3] = Math.sin(theta0) * Math.cos(phi0);// miux->directional
																// cosines
				oldcoord[4] = Math.sin(theta0) * Math.sin(phi0);// miuy->directional
																// cosines
				// oldcoord[5]=Math.cos(theta0);//miuz->directional cosines
				oldcoord[5] = costet;
				// -------------------------------------------END
				// SC-------------------------------
				n_photons = 1;
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
						index_high = index_low;// even if e>0.150MeV!!!
					ph_interp = linInt(energy[index_high],
							ph_coeff[index_high], energy[index_low],
							ph_coeff[index_low], e_incident);
					//incoh_interp = linInt(energy[index_high],
					//		incoh_coeff[index_high], energy[index_low],
					//		incoh_coeff[index_low], e_incident);
					total_interp = linInt(energy[index_high],
							total_coeff[index_high], energy[index_low],
							total_coeff[index_low], e_incident);
					// @@@@@@@@@@@@@@@@@@@@@@@@@RAYLEIGH=coherent scattering!!
					coh_interp = linInt(energy[index_high],
							coh_coeff[index_high], energy[index_low],
							coh_coeff[index_low], e_incident);
					total_interp_r = linInt(energy[index_high],
							total_coeff_r[index_high], energy[index_low],
							total_coeff_r[index_low], e_incident);
					// -------------------------------------------
					coh_probab = (coh_interp) / (total_interp_r - ph_interp);// coh
																				// from
																				// REMAINING
					ph_probab = ph_interp / total_interp_r;// absortion from
															// total
					// --------------------------------------------------------------------------
					r = RandomCollection.random01();
					double rdist = RandomCollection.random01();
					while (rdist == 0.0)
						rdist = RandomCollection.random01();
					dist_to_int = -Math.log(rdist) / (density * total_interp);// [cm]
					if (isValidInteraction != Phantom.NO_ORGAN)// ---------!!!!!!
					{
						interactB = true;// allways force entrance
											// interecation->so pointless
											// n_film_direct!!!
						r = RandomCollection.random01();
						if (r <= ph_probab) // photoelectric interaction
						{
							n_interaction++;
							n_interaction_total = n_interaction_total + 1.0;
							n_photoefect++;
							totalendep = total_energy_deposit + e_incident;
							endep = endep + e_incident;
							//enmedtemp = endep / n_interaction;// per photon,
																// endep is sum
							enindivtemp = e_incident;
							// endeperror=endeperror+Math.pow(enindivtemp-enmedtemp,2);
							endeperror = endeperror + enindivtemp * enindivtemp;
							n_photons = 0;// absorbed photons->end hystory!!
						} else // Compton scattering
						{
							// --------RAYLEIGH--elastic COHERENT----------
							rcoh = RandomCollection.random01();
							if (rcoh <= coh_probab) {
								theta = rayleighSim(e_incident);
								r0 = RandomCollection.random01();// [0-1]
								phi = 2 * r0 * Math.PI;// [0.2PI];
								newcoord = getUpdatedCoordinate(theta, phi,
										dist_to_int, oldcoord);
								isValidInteraction = phantom
										.inWhatOrgan(newcoord);// HERE->organism
								// permute newcoordonate with old coordonate
								oldcoord[0] = newcoord[0];
								oldcoord[1] = newcoord[1];
								oldcoord[2] = newcoord[2];
								oldcoord[3] = newcoord[3];
								oldcoord[4] = newcoord[4];
								oldcoord[5] = newcoord[5];
							} else// COMPTON
							{
								n_interaction++;
								n_interaction_total = n_interaction_total + 1.0;
								theta = comptonSim(e_incident);
								e_scatt = e_incident
										/ (1 + (e_incident / 0.511)
												* (1 - Math.cos(theta)));
								totalendep = total_energy_deposit + e_incident
										- e_scatt;
								endep = endep + e_incident - e_scatt;
								//enmedtemp = endep / n_interaction;// per photon,
																	// endep is
																	// sum
								enindivtemp = e_incident - e_scatt;
								// endeperror=endeperror+Math.pow(enindivtemp-enmedtemp,2);
								endeperror = endeperror + enindivtemp
										* enindivtemp;
								// continue tracing photon hystory
								e_incident = e_scatt;// new loop in while
								// azimutal angle
								r0 = RandomCollection.random01();// [0-1]
								phi = 2 * r0 * Math.PI;// [0.2PI];
								newcoord = getUpdatedCoordinate(theta, phi,
										dist_to_int, oldcoord);
								isValidInteraction = phantom
										.inWhatOrgan(newcoord);// HERE->organism
								// permute newcoordonate with old coordonate
								oldcoord[0] = newcoord[0];
								oldcoord[1] = newcoord[1];
								oldcoord[2] = newcoord[2];
								oldcoord[3] = newcoord[3];
								oldcoord[4] = newcoord[4];
								oldcoord[5] = newcoord[5];
							}
						}
					} else {
						if (!interactB)// no interaction in phantom
						{
							if (oldcoord[2] > phantom.getOrganThickness())// ONLY
																			// HERE,
																			// MAMO
																			// PHANTOM!!
								n_film_direct++;
						}

						// if(Math.sqrt(oldcoord[0]*oldcoord[0]+oldcoord[1]*oldcoord[1])>2.*radius)
						// {
						// n_side++;
						// }
						// pathlength -->exit the phantom
						n_photons = 0;// exit region of interest, so we kill the
										// photon
						n_detected++;

						if (oldcoord[2] > phantom.getOrganThickness())// ONLY
																		// HERE,
																		// MAMO
																		// PHANTOM!!
							n_film++;

						n_side = n_detected - n_film;
					}
				}// end while
			}// end for i=0 to photon_number
				// standard deviation of mean
				// endeperror=Math.sqrt((endeperror/n_interaction)/(n_interaction-1));
			// scale spectrum
			// REAL NORM----v RAD!!
			// endeperror=endeperror/photon_number;
			// -----------------------------------
			// endeperror=endeperror*(xrs.getXRayIntensities())[j]/xrs.getNormalizedValue();
			// normalized per energy-------------------------------------
			endep = endep / photon_number;
			// @@@@@@@@@
			endeperror = endeperror / photon_number;
			endeperror = (endeperror - endep * endep) / (photon_number - 1.0);
			if (endeperror >= 0.0)
				endeperror = Math.sqrt(endeperror);
			// if(endep!=0.)
			// {
			// endeperror= Math.min(endeperror/endep*100.,99.9);
			// }else{endeperror=99.9;}
			endeperror = endeperror * (xrs.getXRayIntensities())[j]
					/ xrs.getNormalizedValue();
			// @@@@@@@
			// ----------------------------------------------
			total_energy_deposit = total_energy_deposit + totalendep;
			// ----------------------------------------------SCALE
			// SPECTRUM----//--------------ORGAN
			energy_deposit[0][0] = energy_deposit[0][0] + endep
					* (xrs.getXRayIntensities())[j] / xrs.getNormalizedValue();
			energy_deposit[0][1] = Math.sqrt(energy_deposit[0][1]
					* energy_deposit[0][1] + endeperror * endeperror);
			// ---------------
		}// end energy loop
		// //% and 2 here for infinit grdlib (Nmare) is 95% conf level
		energy_deposit[0][1] = Math.min(2 * 100 * energy_deposit[0][1]
				/ energy_deposit[0][0], 99.9);// %
		//
		// @@@@@@@@@@@@@@
		if (n_interaction_total > 0.0)
			energy_deposit[0][1] = Math.min(
					2.0 * 100 * Math.sqrt(n_interaction_total)
							/ n_interaction_total, 99.9);// %
		else
			energy_deposit[0][1] = 99.9;// no dose so high error
			// @@@@@@@@@@@@@@@
		double exposure = airExposure * 1000;// -mGy->microGy
		double mevtojoule = 1.60218E-13;// 1MeV=1,60218 � 10-13 J
		double factor = breastTransversalArea * exposure * mevtojoule
				* xrs.getPhotonFlux() / xrs.getAirKerma();
		energy_deposit[0][0] = energy_deposit[0][0] * factor;// Jouli
		organAbsorbedDose[0][0] = 1000 * 1000 * energy_deposit[0][0]
				/ breastMass;// 1g=10-03kk,Gy->mGy
		organAbsorbedDose[0][1] = energy_deposit[0][1];// %
		// --------------------------------
		double wt = ((Double) resources.getObject("breast.wt")).doubleValue();
		effectiveDose[0] = wt * organAbsorbedDose[0][0];
		effectiveDose[1] = organAbsorbedDose[0][1];// %

		simulationTimeElapsed = timeElapsed(startSimulationTime);
	}

	// ------------------------------------------------------------------------------------------------
	// Adus la piele!!!center[0]->x,1-y,2-z; entrance 0->horiz(X),1-vert(Y)
	// entranceXField--la plan median. centrul campului X is real!!
	/**
	 * When simulation starts, this routine computes initial coordinates of particle at patient (MIRD 5 phantom) entrance. Notice 
	 * that centerXField and entranceXField (improperly named "entrance") are related to phantom MID-PLANE. The reason for this is the 
	 * direct correlation with the MIRD 5 phantom equations.
	 * @param centerXField the Xray center position (x,y,z) at phantom mid-plane relative to phantom coordinate system
	 * @param entranceXField the Xray dimensions at phantom mid-plane (width and height)
	 * @param projection the projection index corresponding to AP (0), PA (2), LLAT (1) or RLAT (3)
	 * @param phantom the link to geometry, i.e. the MIRD5 phantom
	 * @param fsd the focus-skin (i.e. focus-entrance) distance in cm.
	 * @return the initial coordinate array.
	 */
	private static double[] getInitialCoord(double[] centerXField,
			double[] entranceXField, int projection, Phantom phantom, double fsd) {
		double[] result = new double[8];// /
		// ------convert all midplane dimension (X,Y) at surface entrance and
		// the vertex
		// first the vertex
		double cxf_z = 0.0;
		double vertexZ = phantom.getVertexCoordonate();
		cxf_z = vertexZ - centerXField[2];
		// ---------------
		double exf_x = 0.0;
		double exf_y = 0.0;
		double thickness = phantom.getOrganProjThickness(projection, cxf_z);// head
																			// or
																			// trunK
		if (fsd > 0) {//real phantom entrance
			exf_x = entranceXField[0] * fsd / (fsd + thickness / 2);
			exf_y = entranceXField[1] * fsd / (fsd + thickness / 2);
		} else// =0!
		{
			exf_x = entranceXField[0];
			exf_y = entranceXField[1];
		}
		// ---------------------------------------------------------
		double x0 = 0.0;
		double y0 = 0.0;
		double z0 = 0.0;
		// --@@@@@@@@@@@@@@@@@@@@----------------------------------
		double theta0 = 0.0;
		double phi0 = 0.0;
		// --@@@@@@@@@@@@@@@@@@@----------------------------------
		double a = 0.0;
		double b = 0.0;
		if (projection == 0 || projection == 2)// AP or PA
		{
			a = centerXField[0] - exf_x / 2;
			b = centerXField[0] + exf_x / 2;
			x0 = RandomCollection.uniform(a, b);

			a = cxf_z - exf_y / 2;
			b = cxf_z + exf_y / 2;
			z0 = RandomCollection.uniform(a, b);
			// --------------------------------
			y0 = phantom.getFittedCoord(projection, z0, x0);
			// ----------------------------------
			double rap = 0.0;
			if (fsd > 0) {
				rap = Math.abs((z0 - cxf_z) / (fsd));
				theta0 = Math.atan(rap);// 0-PI/2
			} else {
				theta0 = 0.0;
			}
			if (x0 - centerXField[0] != 0.0) {
				rap = Math.abs((z0 - cxf_z) / (x0 - centerXField[0]));
				// cadran 1-->0-Pi/2
				if (z0 - cxf_z >= 0 && x0 - centerXField[0] > 0) {
					phi0 = Math.atan(rap);
				}
				// cadran 2-->PI/2-Pi
				if (z0 - cxf_z >= 0 && x0 - centerXField[0] < 0) {
					phi0 = Math.PI - Math.atan(rap);
				}
				// cadran 3-->PI-3Pi/2
				if (z0 - cxf_z <= 0 && x0 - centerXField[0] < 0) {
					phi0 = Math.PI + Math.atan(rap);
				}
				// cadran 4-->PI-2Pi
				if (z0 - cxf_z <= 0 && x0 - centerXField[0] > 0) {
					phi0 = 2 * Math.PI - Math.atan(rap);
				}
			} else {
				phi0 = 0.0;
			}
			// --@@@@@@@@@@@@@@@@@@@@@@@@@@
		}
		if (projection == 1 || projection == 3)// LL or RL
		{
			a = centerXField[1] - exf_x / 2;
			b = centerXField[1] + exf_x / 2;
			y0 = RandomCollection.uniform(a, b);

			a = cxf_z - exf_y / 2;
			b = cxf_z + exf_y / 2;
			z0 = RandomCollection.uniform(a, b);
			// --------------------------------
			x0 = phantom.getFittedCoord(projection, z0, y0);
			// --@@@@@@@@@@@@@@@@@@@@@@@@@
			double rap = 0.0;
			if (fsd > 0)
				rap = Math.abs((z0 - cxf_z) / (fsd));
			else
				rap = 0.0;
			theta0 = Math.atan(rap);// 0-PI/2
			if (z0 - cxf_z != 0.0) {
				rap = Math.abs((y0 - centerXField[1]) / (z0 - cxf_z));
				// cadran 1-->0-Pi/2
				if (y0 - centerXField[1] >= 0 && z0 - cxf_z > 0) {
					phi0 = Math.atan(rap);
				}
				// cadran 2-->PI/2-Pi
				if (y0 - centerXField[1] >= 0 && z0 - cxf_z < 0) {
					phi0 = Math.PI - Math.atan(rap);
				}
				// cadran 3-->PI-3Pi/2
				if (y0 - centerXField[1] <= 0 && z0 - cxf_z < 0) {
					phi0 = Math.PI + Math.atan(rap);
				}
				// cadran 4-->PI-2Pi
				if (y0 - centerXField[1] <= 0 && z0 - cxf_z > 0) {
					phi0 = 2 * Math.PI - Math.atan(rap);
				}
			} else {
				phi0 = 0.0;
			}
			// --@@@@@@@@@@@@@@@@@@@@@@@@@@
		}
		// -----------------directional cosines
		double xsi = 0.0;// angle around x axis. if 0 and lambda =0->from bootom
							// to top
		// if-PI and lambda =0-> from top to bottom=caniocaudal
		double lambda = 0.0;// angle around z axis.
		// if 0 and xsi=Pi/2->LL
		// if PI/2 and xsi=Pi/2->AP
		// if PI and xsi=Pi/2->RL
		// if 3*PI/2 and xsi=Pi/2->PA
		if (projection == 0)// AP
		{
			xsi = Math.PI / 2;
			lambda = Math.PI / 2;
		} else if (projection == 1)// LL
		{
			xsi = Math.PI / 2;
			lambda = 0.0;
		} else if (projection == 2)// PA
		{
			xsi = Math.PI / 2;
			lambda = 3 * Math.PI / 2;
		} else if (projection == 3)// RL
		{
			xsi = Math.PI / 2;
			lambda = Math.PI;
		}
		double cost = Math.cos(theta0);
		double sint = Math.sqrt(1.0 - cost * cost);
		double cosp = Math.cos(phi0);
		double sinp = 0.0;// Math.sqrt(1.0-cosp*cosp);
		if (phi0 < Math.PI)
			sinp = Math.sqrt(1.0 - cosp * cosp);
		else
			sinp = -Math.sqrt(1.0 - cosp * cosp);

		double uz0 = Math.cos(xsi);
		double ux0 = Math.sin(xsi) * Math.cos(lambda);
		double uy0 = Math.sin(xsi) * Math.sin(lambda);
		double ux = 0.0;
		double uy = 0.0;
		double uz = 0.0;
		// rotation with theta and phi and compute the initial real directional
		// cosines:
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

		// -------------------------------------
		result[0] = x0;
		result[1] = y0;
		result[2] = z0;
		result[3] = theta0;
		result[4] = phi0;
		result[5] = ux;
		result[6] = uy;
		result[7] = uz;
		return result;
	}

	/**
	 * Computes the entrance field area in cm^2. Notice 
	 * that centerXField and entranceXField (improperly named "entrance") are related to phantom MID-PLANE. The reason for this is the 
	 * direct correlation with the MIRD 5 phantom equations.
	 * @param fsd focus-skin distance (i.e. focus-entrance distance) in cm.
	 * @param entranceXField the Xray dimensions at phantom mid-plane (width and height)
	 * @param centerXField the Xray center position (x,y,z) at phantom mid-plane relative to phantom coordinate system
	 * @param projection the projection index corresponding to AP (0), PA (2), LLAT (1) or RLAT (3)
	 * @param phantom the link to geometry, i.e. the MIRD5 phantom
	 * @return the entrance field area in cm^2
	 */
	private static double getEntranceField(double fsd, double[] entranceXField,
			double[] centerXField, int projection, Phantom phantom) {
		double result = 0.0;
		// ------convert all midplane dimension (X,Y) at surface entrance and
		// the vertex
		// first the vertex
		double cxf_z = 0.0;
		double vertexZ = phantom.getVertexCoordonate();
		cxf_z = vertexZ - centerXField[2];
		// ---------------
		double exf_x = 0.0;
		double exf_y = 0.0;
		double thickness = phantom.getOrganProjThickness(projection, cxf_z);// head
																			// or
																			// trunK
		if (fsd > 0) {
			exf_x = entranceXField[0] * fsd / (fsd + thickness / 2);
			exf_y = entranceXField[1] * fsd / (fsd + thickness / 2);
		} else// =0!
		{
			exf_x = entranceXField[0];
			exf_y = entranceXField[1];
		}
		// ---------------------------------------------------------
		result = exf_x * exf_y;
		return result;
	}

	// general spectrum
	/**
	 * Perform default Monte Carlo simulation for radiography. Notice 
	 * that centerXField and entranceXField (improperly named "entrance") are related to phantom MID-PLANE. The reason for this is the 
	 * direct correlation with the MIRD 5 phantom equations.
	 * @param photon_number the number of histories for each photon energy.
	 * @param xrs the link to XRaySpectrum class.
	 * @param airExposure the dose free in air at phantom entrance in mGy.
	 * @param fsd focus-skin distance (i.e. focus-entrance distance) in cm.
	 * @param phantom the link to geometry, i.e. the MIRD5 phantom
	 * @param projection the projection index corresponding to AP (0), PA (2), LLAT (1) or RLAT (3)
	 * @param centerXField the Xray center position (x,y,z) at phantom mid-plane relative to phantom coordinate system
	 * @param entranceXField the Xray dimensions at phantom mid-plane (width and height)
	 */
	public static void radSim(int photon_number, XRaySpectrum xrs,
			double airExposure, double fsd, Phantom phantom, int projection,
			double[] centerXField, double[] entranceXField) {
		long startSimulationTime = System.currentTimeMillis();
		simulationTimeElapsed = "";
		int index_low = 0;
		int index_high = 0;
		double ph_interp = 0.0;
		//double incoh_interp = 0.0;// Compton=incoherent scattering!!
		//double total_interp = 0.0;
		double ph_probab = 0.0;
		double r = 0.0;// [0,1] random number
		double dist_to_int = 0.0;// store distance to interaction
		// -------------SC-------------------------------------------------------------
		double r0 = 0.0;
		double[] oldcoord = new double[8];
		// -----------------------------
		double[] oldcrd = new double[6];
		//double theta0 = 0.0;
		//double phi0 = 0.0;
		// ---------------------
		double phi = 0.0;
		double[] newcoord;
		int isValidInteraction =Phantom.NO_ORGAN;// INIT ANYWAY
		// ----------------END
		// SC-----------------------------------------------------
		// Area in mm2
		double transversalArea = getEntranceField(fsd, entranceXField,
				centerXField, projection, phantom) * 100;
		// -----------
		double e_incident = 0.0;// init
		double e_scatt = 0.0;// store the energy of the Compton scattered
								// photons
		// -----------INIT ORGANS BASED VARIABLES!!---------------------------
		organs = (String[]) resources.getObject("rad.organname");
		energy_deposit = new double[organs.length][2];// all organs
		organAbsorbedDose = new double[organs.length][2];// all organs
		int[] contorInteraction = new int[organs.length];
		double[] contorInteraction_total = new double[organs.length];
		//double dummy = 0.0;
		double[] endep = new double[organs.length];// temporary
		double[] endeperror = new double[organs.length];// temporary
		for (int i = 0; i < organs.length; i++)// just in case
		{
			energy_deposit[i][0] = 0.0;
			energy_deposit[i][1] = 0.0;
			organAbsorbedDose[i][0] = 0.0;
			organAbsorbedDose[i][1] = 0.0;
			// ------------------------
			endep[i] = 0.0;
			endeperror[i] = 0.0;
			contorInteraction[i] = 0;
			contorInteraction_total[i] = 0.0;
		}
		// ---------------------------------
		//double enmedtemp = 0.0;// handle error eval
		double enindivtemp = 0.0;// handle error eval
		// ----------------------------------
		effectiveDose = new double[2];
		effectiveDose[0] = 0.0;
		effectiveDose[1] = 0.0;
		n_detected = 0;
		n_total = 0;
		n_film = 0;
		n_film_direct = 0;
		n_side = 0;
		n_photoefect = 0;
		double theta = 0.0;// store the Compton simulation's scattering angle
		int n_photons = 0;// number of photons being transported
		// ------read attenuation coefficient and density----------
		double[][] coeftable = new double[0][0];
		double[][] coeftableSk = new double[0][0];
		double[][] coeftableLu = new double[0][0];
		double[][] coeftableBr = (double[][]) resources
				.getObject("breast.diagnostic.attenuationCoef");
		double density = 0.0;
		double densitySk = 0.0;
		double densityLu = 0.0;
		double densityBr = ((Double) resources
				.getObject("breast.diagnostic.density")).doubleValue();

		if (phantom.getIndex() != Phantom.NEWBORN_INDEX) {
			coeftable = (double[][]) resources
					.getObject("commonTissue.diagnostic.attenuationCoef");
			coeftableSk = (double[][]) resources
					.getObject("skeletal_boneWithMarrow.diagnostic.attenuationCoef");
			coeftableLu = (double[][]) resources
					.getObject("lung.diagnostic.attenuationCoef");
			density = ((Double) resources
					.getObject("commonTissue.diagnostic.density"))
					.doubleValue();
			densitySk = ((Double) resources
					.getObject("skeletal_boneWithMarrow.diagnostic.density"))
					.doubleValue();
			densityLu = ((Double) resources
					.getObject("lung.diagnostic.density")).doubleValue();

		} else // MAMO WILL NEVER HAPPEN HERE!!
		{
			coeftable = (double[][]) resources
					.getObject("commonTissue.newborn.diagnostic.attenuationCoef");
			coeftableSk = (double[][]) resources
					.getObject("skeletal_boneWithMarrow.newborn.diagnostic.attenuationCoef");
			coeftableLu = (double[][]) resources
					.getObject("lung.newborn.diagnostic.attenuationCoef");
			density = ((Double) resources
					.getObject("commonTissue.newborn.diagnostic.density"))
					.doubleValue();
			densitySk = ((Double) resources
					.getObject("skeletal_boneWithMarrow.newborn.diagnostic.density"))
					.doubleValue();
			densityLu = ((Double) resources
					.getObject("lung.newborn.diagnostic.density"))
					.doubleValue();
		}
		double dens = 0.0;// handle densities
		// coeftable[i][0])-->energy in MeV;
		// coeftable[i][2])-->Compton att. coef in cm2/g;
		// coeftable[i][3])-->Photoelectric abs. coef in cm2/g;
		// coeftable[i][5])-->Total (without Rayleigh scattering) att. coef in
		// cm2/g;
		// read data:
		double[] energy = new double[coeftable.length];
		double[] ph_coeff = new double[coeftable.length];
		double[] incoh_coeff = new double[coeftable.length];
		double[] total_coeff = new double[coeftable.length];
		double[] ph_coeffSk = new double[coeftableSk.length];// same energy!!
		double[] incoh_coeffSk = new double[coeftableSk.length];
		double[] total_coeffSk = new double[coeftableSk.length];
		double[] ph_coeffLu = new double[coeftableLu.length];
		double[] incoh_coeffLu = new double[coeftableLu.length];
		double[] total_coeffLu = new double[coeftableLu.length];
		double[] ph_coeffBr = new double[coeftableBr.length];
		double[] incoh_coeffBr = new double[coeftableBr.length];
		double[] total_coeffBr = new double[coeftableBr.length];
		// ----------------------------------------------------------------------
		double[] coh_coeff = new double[coeftable.length];
		double[] total_coeff_r = new double[coeftable.length];// @@@@@@@RAYLEIGH
		double[] coh_coeffSk = new double[coeftableSk.length];
		double[] total_coeff_rSk = new double[coeftableSk.length];// @@@@@@@RAYLEIGH
		double[] coh_coeffLu = new double[coeftableLu.length];
		double[] total_coeff_rLu = new double[coeftableLu.length];// @@@@@@@RAYLEIGH
		double[] coh_coeffBr = new double[coeftableBr.length];
		double[] total_coeff_rBr = new double[coeftableBr.length];// @@@@@@@RAYLEIGH

		// separated ABM data
		double[][] abmraptable = (double[][]) resources
				.getObject("abm.sk.fraction");
		double[] abmrap = new double[abmraptable.length];
		double[] abmwtab = (double[]) resources.getObject("abm.weight");
		// NewBorn-index 1 and in tab is index 0!!
		double abmw = abmwtab[phantom.getIndex() - 1];

		for (int i = 0; i < energy.length; i++) {
			energy[i] = coeftable[i][0];
			incoh_coeff[i] = coeftable[i][2];
			ph_coeff[i] = coeftable[i][3];
			total_coeff[i] = coeftable[i][5];
			incoh_coeffSk[i] = coeftableSk[i][2];
			ph_coeffSk[i] = coeftableSk[i][3];
			total_coeffSk[i] = coeftableSk[i][5];
			incoh_coeffLu[i] = coeftableLu[i][2];
			ph_coeffLu[i] = coeftableLu[i][3];
			total_coeffLu[i] = coeftableLu[i][5];
			incoh_coeffBr[i] = coeftableBr[i][2];
			ph_coeffBr[i] = coeftableBr[i][3];
			total_coeffBr[i] = coeftableBr[i][5];

			abmrap[i] = abmraptable[i][1];
			// --------------------------------------------------------@@@@@@@RAYLEIGH
			coh_coeff[i] = coeftable[i][1];
			total_coeff_r[i] = coeftable[i][4];
			coh_coeffSk[i] = coeftableSk[i][1];
			total_coeff_rSk[i] = coeftableSk[i][4];
			coh_coeffLu[i] = coeftableLu[i][1];
			total_coeff_rLu[i] = coeftableLu[i][4];
			coh_coeffLu[i] = coeftableLu[i][1];
			total_coeff_rLu[i] = coeftableLu[i][4];
		}
		double abmrapinterp = 0.0;
		double abmtempen = 0.0;
		double abmtemperror = 0.0;
		// -------------------------------------------
		double coh_probab = 0.0;
		double coh_interp = 0.0;// RAYLEIGH=coherent scattering!!
		double total_interp_r = 0.0;
		double rcoh = 0.0;// ALEATOR NUMBER
		// ------------------------------------@@@@@@@@@@@@@@@@@@@@@@@@RAYLEIGH
		// -------MONTE CARLO
		// SIMULATION---------------------------------------------------
		boolean interactB = false;
		for (int j = 0; j < (xrs.getXRayEnergies()).length; j++) {
			// temporary zero init
			for (int k = 0; k < organs.length; k++) {
				endep[k] = 0.0;
				endeperror[k] = 0.0;
				contorInteraction[k] = 0;
			}
			// abm--------------------
			abmtempen = 0.0;
			abmtemperror = 0.0;
			// ----------MC BASED LOOP---------------------------
			for (int i = 0; i < photon_number; i++) {
				interactB = false;
				n_total++;
				e_incident = (xrs.getXRayEnergies())[j] / 1000;// incident_energy
																// in MeV;
				// ----------------------------------------------SC--------------------------------
				oldcoord = getInitialCoord(centerXField, entranceXField,
						projection, phantom, fsd);
				// prelucration----------@@@@@@@@@@----------------
				oldcrd[0] = oldcoord[0];// x
				oldcrd[1] = oldcoord[1];// y
				oldcrd[2] = oldcoord[2];// z
				//theta0 = oldcoord[3];
				//phi0 = oldcoord[4];
				oldcrd[3] = oldcoord[5];// ux
				oldcrd[4] = oldcoord[6];// uy
				oldcrd[5] = oldcoord[7];// uz
				// -------@@@@@@@@@@@@@@@----------------
				// test if we are in body or not
				boolean inbody = phantom.inBody(oldcoord);
				if (inbody)
					isValidInteraction = Phantom.REMAINDER;// at skin entrance
															// SKIN=REMAINDER->REMOVE
				else
					isValidInteraction = Phantom.NO_ORGAN;// exit or miss the
															// phantom

				// -------------------------------------------END
				// SC-------------------------------
				n_photons = 1;
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
						index_high = index_low;// even if e>0.150MeV!!!
					// --------------choose interaction
					// coefficients----------------------------------
					if (isValidInteraction == Phantom.SKELETON) {
						ph_interp = linInt(energy[index_high],
								ph_coeffSk[index_high], energy[index_low],
								ph_coeffSk[index_low], e_incident);
						//incoh_interp = linInt(energy[index_high],
						//		incoh_coeffSk[index_high], energy[index_low],
						//		incoh_coeffSk[index_low], e_incident);
						//total_interp = linInt(energy[index_high],
						//		total_coeffSk[index_high], energy[index_low],
						//		total_coeffSk[index_low], e_incident);
						dens = densitySk;
						// ABM
						abmrapinterp = linInt(energy[index_high],
								abmrap[index_high], energy[index_low],
								abmrap[index_low], e_incident);
						// @@@@@@@@@@@@@@@@@@@@@@@@@RAYLEIGH=coherent
						// scattering!!
						coh_interp = linInt(energy[index_high],
								coh_coeffSk[index_high], energy[index_low],
								coh_coeffSk[index_low], e_incident);
						total_interp_r = linInt(energy[index_high],
								total_coeff_rSk[index_high], energy[index_low],
								total_coeff_rSk[index_low], e_incident);

					} else if (isValidInteraction == Phantom.LUNGS) {
						ph_interp = linInt(energy[index_high],
								ph_coeffLu[index_high], energy[index_low],
								ph_coeffLu[index_low], e_incident);
						//incoh_interp = linInt(energy[index_high],
						//		incoh_coeffLu[index_high], energy[index_low],
						//		incoh_coeffLu[index_low], e_incident);
						//total_interp = linInt(energy[index_high],
						//		total_coeffLu[index_high], energy[index_low],
						//		total_coeffLu[index_low], e_incident);
						dens = densityLu;
						// @@@@@@@@@@@@@@@@@@@@@@@@@RAYLEIGH=coherent
						// scattering!!
						coh_interp = linInt(energy[index_high],
								coh_coeffLu[index_high], energy[index_low],
								coh_coeffLu[index_low], e_incident);
						total_interp_r = linInt(energy[index_high],
								total_coeff_rLu[index_high], energy[index_low],
								total_coeff_rLu[index_low], e_incident);

					} else if (isValidInteraction == Phantom.BREASTS) {
						ph_interp = linInt(energy[index_high],
								ph_coeffBr[index_high], energy[index_low],
								ph_coeffBr[index_low], e_incident);
						//incoh_interp = linInt(energy[index_high],
						//		incoh_coeffBr[index_high], energy[index_low],
						//		incoh_coeffBr[index_low], e_incident);
						//total_interp = linInt(energy[index_high],
						//		total_coeffBr[index_high], energy[index_low],
						//		total_coeffBr[index_low], e_incident);
						dens = densityBr;
						// @@@@@@@@@@@@@@@@@@@@@@@@@RAYLEIGH=coherent
						// scattering!!
						coh_interp = linInt(energy[index_high],
								coh_coeffBr[index_high], energy[index_low],
								coh_coeffBr[index_low], e_incident);
						total_interp_r = linInt(energy[index_high],
								total_coeff_rBr[index_high], energy[index_low],
								total_coeff_rBr[index_low], e_incident);

					} else// soft tissue-common or noorgan but this is later
							// resolved
					{
						ph_interp = linInt(energy[index_high],
								ph_coeff[index_high], energy[index_low],
								ph_coeff[index_low], e_incident);
						//incoh_interp = linInt(energy[index_high],
						//		incoh_coeff[index_high], energy[index_low],
						//		incoh_coeff[index_low], e_incident);
						//total_interp = linInt(energy[index_high],
						//		total_coeff[index_high], energy[index_low],
						//		total_coeff[index_low], e_incident);
						dens = density;
						// @@@@@@@@@@@@@@@@@@@@@@@@@RAYLEIGH=coherent
						// scattering!!
						coh_interp = linInt(energy[index_high],
								coh_coeff[index_high], energy[index_low],
								coh_coeff[index_low], e_incident);
						total_interp_r = linInt(energy[index_high],
								total_coeff_r[index_high], energy[index_low],
								total_coeff_r[index_low], e_incident);
					}
					// ph_probab=ph_interp/total_interp;
					// -------------------------------------------
					coh_probab = (coh_interp) / (total_interp_r - ph_interp);// coh
																				// from
																				// REMAINING
					ph_probab = ph_interp / total_interp_r;// absortion from
															// total
					// --------------------------------------------------------------------------
					// --END
					// COEFF--------------------------------------------------------------------
					r = RandomCollection.random01();
					double rdist = RandomCollection.random01();
					while (rdist == 0.0)
						rdist = RandomCollection.random01();
					dist_to_int = -Math.log(rdist) / (dens * total_interp_r);// [cm]!!@

					if (isValidInteraction != Phantom.NO_ORGAN)// ---------!!!!!!
					{
						interactB = true;
						r = RandomCollection.random01();
						if (r <= ph_probab) // photoelectric interaction
						{
							contorInteraction[isValidInteraction - 1]++;// interaction
																		// colission
							contorInteraction_total[isValidInteraction - 1] = contorInteraction_total[isValidInteraction - 1] + 1.0;
							n_photoefect++;
							endep[isValidInteraction - 1] = endep[isValidInteraction - 1]
									+ e_incident;// --!!!!!!!
							// from phantom if
							// isValidInteraction=phantom.REMAINDER=25->
							// in organs coresponds 24!!
							//enmedtemp = endep[isValidInteraction - 1]
							//		/ contorInteraction[isValidInteraction - 1];
							// per photon, endep is sum
							enindivtemp = e_incident;
							// endeperror[isValidInteraction-1]=endeperror[isValidInteraction-1]+
							// Math.pow(enindivtemp-enmedtemp,2);//error at
							// power 2
							endeperror[isValidInteraction - 1] = endeperror[isValidInteraction - 1]
									+ enindivtemp * enindivtemp;
							// separated ABM-----for skeleton
							// eval---------------------------------
							if (isValidInteraction == Phantom.SKELETON) {
								abmtempen = abmtempen + abmrapinterp * abmw
										* enindivtemp;
								// abmtemperror=abmtemperror+
								// Math.pow(abmrapinterp*abmw*enindivtemp-abmtempen/contorInteraction[isValidInteraction-1],2);
								abmtemperror = abmtemperror + abmrapinterp
										* abmw * enindivtemp * abmrapinterp
										* abmw * enindivtemp;
							}
							// So=>absorbed photons->end hystory!!
							n_photons = 0;
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
								isValidInteraction = phantom
										.inWhatOrgan(newcoord);// HERE->organism
								// permute newcoordonate with old coordonate
								oldcrd[0] = newcoord[0];
								oldcrd[1] = newcoord[1];
								oldcrd[2] = newcoord[2];
								oldcrd[3] = newcoord[3];
								oldcrd[4] = newcoord[4];
								oldcrd[5] = newcoord[5];
							} else// COMPTON
							{
								contorInteraction[isValidInteraction - 1]++;// interaction
																			// colission
								contorInteraction_total[isValidInteraction - 1] = contorInteraction_total[isValidInteraction - 1] + 1.0;
								theta = comptonSim(e_incident);
								e_scatt = e_incident
										/ (1 + (e_incident / 0.511)
												* (1 - Math.cos(theta)));
								endep[isValidInteraction - 1] = endep[isValidInteraction - 1]
										+ e_incident - e_scatt;

								//enmedtemp = endep[isValidInteraction - 1]
								//		/ contorInteraction[isValidInteraction - 1];
								enindivtemp = e_incident - e_scatt;
								// endeperror[isValidInteraction-1]=endeperror[isValidInteraction-1]+
								// Math.pow(enindivtemp-enmedtemp,2);
								endeperror[isValidInteraction - 1] = endeperror[isValidInteraction - 1]
										+ enindivtemp * enindivtemp;
								// separated ABM-----for skeleton
								// eval--------------------------------
								if (isValidInteraction == Phantom.SKELETON) {
									abmtempen = abmtempen + abmrapinterp * abmw
											* enindivtemp;
									// abmtemperror=abmtemperror+
									// Math.pow(abmrapinterp*abmw*enindivtemp-abmtempen/contorInteraction[isValidInteraction-1],2);
									abmtemperror = abmtemperror + abmrapinterp
											* abmw * enindivtemp * abmrapinterp
											* abmw * enindivtemp;
								}

								// continue tracing photon hystory
								e_incident = e_scatt;// new loop in while
								// azimutal angle
								r0 = RandomCollection.random01();// [0-1]
								phi = 2 * r0 * Math.PI;// ->[0.2PI];
								// ----------------------------------------RAD---------------
								newcoord = getUpdatedCoordinate(theta, phi,
										dist_to_int, oldcrd);
								// -------------------!!!!!!!!!_____________________________
								isValidInteraction = phantom
										.inWhatOrgan(newcoord);// retrieve the
																// organ
								// HERE->organisminclude inBody!
								// permute newcoordonate with old coordonate
								oldcrd[0] = newcoord[0];
								oldcrd[1] = newcoord[1];
								oldcrd[2] = newcoord[2];
								oldcrd[3] = newcoord[3];
								oldcrd[4] = newcoord[4];
								oldcrd[5] = newcoord[5];
							}
						}
					} else {
						// pathlength -->exit the phantom
						n_photons = 0;// exit region of interest, so we kill the
										// photon
						n_detected++;
						if (projection == 0)// corect are loc initial
											// mulajul----legs=virtual=cu trunk
						{
							if (oldcrd[1] > 0.0)// phantom.getOrganProjThickness(projection,oldcrd[2])/2)
							{
								if (!interactB)
									n_film_direct++;
								n_film++;
							}
						}
						if (projection == 2) {
							if (oldcrd[1] < 0.0)// -phantom.getOrganProjThickness(projection,oldcrd[2])/2)
							{
								if (!interactB)
									n_film_direct++;
								n_film++;
							}
						}
						if (projection == 1) {
							if (oldcrd[0] > 0.0)// phantom.getOrganProjThickness(projection,oldcrd[2])/2)
							{
								if (!interactB)
									n_film_direct++;
								n_film++;
							}
						}
						if (projection == 3) {
							if (oldcrd[0] < 0.0)// -phantom.getOrganProjThickness(projection,oldcrd[2])/2)
							{
								if (!interactB)
									n_film_direct++;
								n_film++;
							}
						}

						n_side = n_detected - n_film;
					}
				}// end while
			}// end for i=0 to photon_number
				// NOW WE HAVE DISTRIBUTION OF ENERGY IN ALL ORGANS FOR e.g
				// 10000 phgoton hystories.
				// error and we try the follow code.MAYBE DUE TO integer
				// multiplication !!!!!!!!!!!@@@@@@@@!!
			for (int k = 0; k < organs.length; k++) {
				/*
				 * if (contorInteraction[k]>0) { if (contorInteraction[k]>1) {
				 * //standard deviation of mean endeperror[k]=
				 * Math.sqrt((endeperror
				 * [k]/contorInteraction[k])/(contorInteraction[k]-1)); }//else
				 * is unchanged!! }//else unchanged
				 */
				// ---REAL @@ normalized per
				// energy-------------------------------------
				// endeperror[k]=endeperror[k]/photon_number;
				// scale spectrum
				// endeperror[k]=endeperror[k]*(xrs.getXRayIntensities())[j]/xrs.getNormalizedValue();
				// REAL @@ normalized per
				// energy-------------------------------------
				endep[k] = endep[k] / photon_number;
				// @@@@@@@@@@@@
				endeperror[k] = endeperror[k] / photon_number;
				endeperror[k] = (endeperror[k] - endep[k] * endep[k])
						/ (photon_number - 1.0);
				if (endeperror[k] >= 0.0)
					endeperror[k] = Math.sqrt(endeperror[k]);
				endeperror[k] = endeperror[k] * (xrs.getXRayIntensities())[j]
						/ xrs.getNormalizedValue();
				// @@@@@@@@@@@@@
				// ----------------------------------------------SCALE
				// SPECTRUM---------------
				endep[k] = endep[k] * (xrs.getXRayIntensities())[j]
						/ xrs.getNormalizedValue();
				// ----------------NOW REAL
				// VARIABLES------------------------------------------
				energy_deposit[k][0] = energy_deposit[k][0] + endep[k];
				energy_deposit[k][1] = Math.sqrt(energy_deposit[k][1]
						* energy_deposit[k][1] + endeperror[k] * endeperror[k]);
			}
			// separated ABM--from skeleton eval
			/*
			 * if (contorInteraction[phantom.SKELETON-1]>0) { if
			 * (contorInteraction[phantom.SKELETON-1]>1) { //standard deviation
			 * of mean abmtemperror=
			 * Math.sqrt((abmtemperror/contorInteraction[phantom
			 * .SKELETON-1])/(contorInteraction[phantom.SKELETON-1]-1)); }//else
			 * is unchanged!! }//else unchanged
			 */
			// ---REAL @@ normalized per
			// energy-------------------------------------
			// abmtemperror=abmtemperror/photon_number;
			// scale spectrum
			// abmtemperror=abmtemperror*(xrs.getXRayIntensities())[j]/xrs.getNormalizedValue();
			// ---REAL @@ normalized per
			// energy-------------------------------------
			abmtempen = abmtempen / photon_number;
			abmtemperror = abmtemperror / photon_number;
			abmtemperror = (abmtemperror - abmtempen * abmtempen)
					/ (photon_number - 1.0);
			if (abmtemperror >= 0.0)
				abmtemperror = Math.sqrt(abmtemperror);
			abmtemperror = abmtemperror * (xrs.getXRayIntensities())[j]
					/ xrs.getNormalizedValue();
			// ----------------------------------------------SCALE
			// SPECTRUM---------------
			abmtempen = abmtempen * (xrs.getXRayIntensities())[j]
					/ xrs.getNormalizedValue();
			// ----------------NOW REAL VARIABLES IN ABM NOT IN SKELETON-DIFF OF
			// 1!!!------------
			energy_deposit[Phantom.ACTVEBONEMARROW - 1][0] = energy_deposit[Phantom.ACTVEBONEMARROW - 1][0]
					+ abmtempen;
			energy_deposit[Phantom.ACTVEBONEMARROW - 1][1] = Math
					.sqrt(energy_deposit[Phantom.ACTVEBONEMARROW - 1][1]
							* energy_deposit[Phantom.ACTVEBONEMARROW - 1][1]
							+ abmtemperror * abmtemperror);
			// ---------------------
		}// end energy loop
		double exposure = airExposure * 1000;// -mGy->microGy
		double mevtojoule = 1.60218E-13;// 1MeV=1,60218 � 10-13 J
		double factor = transversalArea * exposure * mevtojoule
				* xrs.getPhotonFlux() / xrs.getAirKerma();
		// separated ABM
		if (energy_deposit[Phantom.ACTVEBONEMARROW - 1][0] != 0.0)
			energy_deposit[Phantom.ACTVEBONEMARROW - 1][1] = Math.min(2 * 100
					* energy_deposit[Phantom.ACTVEBONEMARROW - 1][1]
					/ energy_deposit[Phantom.ACTVEBONEMARROW - 1][0], 99.9);// %
		else
			energy_deposit[Phantom.ACTVEBONEMARROW - 1][1] = 99.9;// 0.0;
			// normalized at real XRay Spectrum:
		energy_deposit[Phantom.ACTVEBONEMARROW - 1][0] = energy_deposit[Phantom.ACTVEBONEMARROW - 1][0]
				* factor;// Jouli
		// organ dose:
		organAbsorbedDose[Phantom.ACTVEBONEMARROW - 1][0] = 1000 * 1000
				* energy_deposit[Phantom.ACTVEBONEMARROW - 1][0]
				/ phantom.getOrganMass(Phantom.ACTVEBONEMARROW - 1);
		// 1g=10-03kk,Gy->mGy
		organAbsorbedDose[Phantom.ACTVEBONEMARROW - 1][1] = energy_deposit[Phantom.ACTVEBONEMARROW - 1][1];// %
		// -------------------------------
		// 33333333333333333333333333333333333333333333333
		for (int k = 0; k < organs.length; k++) {
			if (contorInteraction_total[k] > 0.0)
				organAbsorbedDose[k][1] = Math.min(
						2.0 * 100.0 * Math.sqrt(contorInteraction_total[k])
								/ contorInteraction_total[k], 99.9);
			else
				organAbsorbedDose[k][1] = 99.9;
		}
		organAbsorbedDose[Phantom.ACTVEBONEMARROW - 1][1] = organAbsorbedDose[Phantom.SKELETON - 1][1];
		// 33333333333333333333333333333333333333333333333
		double[] wt = ((double[]) resources.getObject("wt.organname"));
		effectiveDose[0] = 0.0;
		double seff2 = 0.0;
		double errmin = organAbsorbedDose[Phantom.ACTVEBONEMARROW - 1][1];// initialisation
		for (int k = 0; k < organs.length; k++) {
			if (k != Phantom.ACTVEBONEMARROW - 1) {
				// //% and 2 here for infinit grdlib (Nmare) is 95% conf level
				if (energy_deposit[k][0] != 0.0)
					energy_deposit[k][1] = Math
							.min(2 * 100 * energy_deposit[k][1]
									/ energy_deposit[k][0], 99.9);// ;//%
				else
					energy_deposit[k][1] = 99.9;// 0.0;
				// normalized at real XRay Spectrum:
				energy_deposit[k][0] = energy_deposit[k][0] * factor;// Jouli
				// organ dose:
				organAbsorbedDose[k][0] = 1000 * 1000 * energy_deposit[k][0]
						/ phantom.getOrganMass(k);
				// 1g=10-03kk,Gy->mGy
				// organAbsorbedDose[k][1]=energy_deposit[k][1];//%
				errmin = Math.min(errmin, organAbsorbedDose[k][1]);
			}
			effectiveDose[0] = effectiveDose[0] + wt[k]
					* organAbsorbedDose[k][0];
			// --error->proc readus la doza
			seff2 = seff2
					+ Math.pow(wt[k] * organAbsorbedDose[k][0]
							* organAbsorbedDose[k][1] / 100, 2);
		}
		if (effectiveDose[0] != 0.0)
			effectiveDose[1] = errmin;// Math.min((100*Math.sqrt(seff2))/effectiveDose[0],99.9);//%
		else
			effectiveDose[1] = 99.9;// 0.0;

		simulationTimeElapsed = timeElapsed(startSimulationTime);
	}

	// -------------------@@@@@@@@@@NEW
	// ALGOR@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	/**
	 * Perform slightly modified Monte Carlo simulation for radiography using Rosenstein algorithm for evaluating dose. Notice 
	 * that centerXField and entranceXField (improperly named "entrance") are related to phantom MID-PLANE. The reason for this is the 
	 * direct correlation with the MIRD 5 phantom equations.
	 * @param photon_number the number of histories for each photon energy.
	 * @param xrs the link to XRaySpectrum class.
	 * @param airExposure the dose free in air at phantom entrance in mGy.
	 * @param fsd focus-skin distance (i.e. focus-entrance distance) in cm.
	 * @param phantom the link to geometry, i.e. the MIRD5 phantom
	 * @param projection the projection index corresponding to AP (0), PA (2), LLAT (1) or RLAT (3)
	 * @param centerXField the Xray center position (x,y,z) at phantom mid-plane relative to phantom coordinate system
	 * @param entranceXField the Xray dimensions at phantom mid-plane (width and height)
	 */
	public static void radWeightedSim(int photon_number, XRaySpectrum xrs,
			double airExposure, double fsd, Phantom phantom, int projection,
			double[] centerXField, double[] entranceXField) {
		long startSimulationTime = System.currentTimeMillis();
		simulationTimeElapsed = "";
		int index_low = 0;
		int index_high = 0;
		double ph_interp = 0.0;
		//double incoh_interp = 0.0;// Compton=incoherent scattering!!
		//double total_interp = 0.0;
		double ph_probab = 0.0;
		//double r = 0.0;// [0,1] random number
		double dist_to_int = 0.0;// store distance to interaction
		// -------------SC-------------------------------------------------------------
		double r0 = 0.0;
		double[] oldcoord = new double[8];
		// -----------------------------
		double[] oldcrd = new double[6];
		//double theta0 = 0.0;
		//double phi0 = 0.0;
		// ---------------------
		double phi = 0.0;
		double[] newcoord;
		int isValidInteraction = Phantom.NO_ORGAN;// INIT ANYWAY
		// ----------------END
		// SC-----------------------------------------------------
		// Area in mm2
		double transversalArea = getEntranceField(fsd, entranceXField,
				centerXField, projection, phantom) * 100;
		// -----------
		double e_incident = 0.0;// init
		double e_scatt = 0.0;// store the energy of the Compton scattered
								// photons
		// -----------INIT ORGANS BASED VARIABLES!!---------------------------
		organs = (String[]) resources.getObject("rad.organname");
		energy_deposit = new double[organs.length][2];// all organs
		organAbsorbedDose = new double[organs.length][2];// all organs
		int[] contorInteraction = new int[organs.length];
		double[] contorInteraction_total = new double[organs.length];
		double[] endep = new double[organs.length];// temporary
		double[] endeperror = new double[organs.length];// temporary
		for (int i = 0; i < organs.length; i++)// just in case
		{
			energy_deposit[i][0] = 0.0;
			energy_deposit[i][1] = 0.0;
			organAbsorbedDose[i][0] = 0.0;
			organAbsorbedDose[i][1] = 0.0;
			// ------------------------
			endep[i] = 0.0;
			endeperror[i] = 0.0;
			contorInteraction[i] = 0;
			contorInteraction_total[i] = 0.0;
		}
		// ---------------------------------
		//double enmedtemp = 0.0;// handle error eval
		double enindivtemp = 0.0;// handle error eval
		// ----------------------------------
		effectiveDose = new double[2];
		effectiveDose[0] = 0.0;
		effectiveDose[1] = 0.0;
		n_detected = 0;
		n_total = 0;
		n_film = 0;
		n_film_direct = 0;
		n_side = 0;
		n_photoefect = 0;
		double theta = 0.0;// store the Compton simulation's scattering angle
		int n_photons = 0;// number of photons being transported
		// ------read attenuation coefficient and density----------
		double[][] coeftable = new double[0][0];
		double[][] coeftableSk = new double[0][0];
		double[][] coeftableLu = new double[0][0];
		double[][] coeftableBr = (double[][]) resources
				.getObject("breast.diagnostic.attenuationCoef");
		double density = 0.0;
		double densitySk = 0.0;
		double densityLu = 0.0;
		double densityBr = ((Double) resources
				.getObject("breast.diagnostic.density")).doubleValue();

		if (phantom.getIndex() != Phantom.NEWBORN_INDEX) {
			coeftable = (double[][]) resources
					.getObject("commonTissue.diagnostic.attenuationCoef");
			coeftableSk = (double[][]) resources
					.getObject("skeletal_boneWithMarrow.diagnostic.attenuationCoef");
			coeftableLu = (double[][]) resources
					.getObject("lung.diagnostic.attenuationCoef");
			density = ((Double) resources
					.getObject("commonTissue.diagnostic.density"))
					.doubleValue();
			densitySk = ((Double) resources
					.getObject("skeletal_boneWithMarrow.diagnostic.density"))
					.doubleValue();
			densityLu = ((Double) resources
					.getObject("lung.diagnostic.density")).doubleValue();

		} else // MAMO WILL NEVER HAPPEN HERE!!
		{
			coeftable = (double[][]) resources
					.getObject("commonTissue.newborn.diagnostic.attenuationCoef");
			coeftableSk = (double[][]) resources
					.getObject("skeletal_boneWithMarrow.newborn.diagnostic.attenuationCoef");
			coeftableLu = (double[][]) resources
					.getObject("lung.newborn.diagnostic.attenuationCoef");
			density = ((Double) resources
					.getObject("commonTissue.newborn.diagnostic.density"))
					.doubleValue();
			densitySk = ((Double) resources
					.getObject("skeletal_boneWithMarrow.newborn.diagnostic.density"))
					.doubleValue();
			densityLu = ((Double) resources
					.getObject("lung.newborn.diagnostic.density"))
					.doubleValue();
		}
		double dens = 0.0;// handle densities
		// coeftable[i][0])-->energy in MeV;
		// coeftable[i][2])-->Compton att. coef in cm2/g;
		// coeftable[i][3])-->Photoelectric abs. coef in cm2/g;
		// coeftable[i][5])-->Total (without Rayleigh scattering) att. coef in
		// cm2/g;
		// read data:
		double[] energy = new double[coeftable.length];
		double[] ph_coeff = new double[coeftable.length];
		double[] incoh_coeff = new double[coeftable.length];
		double[] total_coeff = new double[coeftable.length];
		double[] ph_coeffSk = new double[coeftableSk.length];// same energy!!
		double[] incoh_coeffSk = new double[coeftableSk.length];
		double[] total_coeffSk = new double[coeftableSk.length];
		double[] ph_coeffLu = new double[coeftableLu.length];
		double[] incoh_coeffLu = new double[coeftableLu.length];
		double[] total_coeffLu = new double[coeftableLu.length];
		double[] ph_coeffBr = new double[coeftableBr.length];
		double[] incoh_coeffBr = new double[coeftableBr.length];
		double[] total_coeffBr = new double[coeftableBr.length];
		// ----------------------------------------------------------------------
		double[] coh_coeff = new double[coeftable.length];
		double[] total_coeff_r = new double[coeftable.length];// @@@@@@@RAYLEIGH
		double[] coh_coeffSk = new double[coeftableSk.length];
		double[] total_coeff_rSk = new double[coeftableSk.length];// @@@@@@@RAYLEIGH
		double[] coh_coeffLu = new double[coeftableLu.length];
		double[] total_coeff_rLu = new double[coeftableLu.length];// @@@@@@@RAYLEIGH
		double[] coh_coeffBr = new double[coeftableBr.length];
		double[] total_coeff_rBr = new double[coeftableBr.length];// @@@@@@@RAYLEIGH

		// separated ABM data
		double[][] abmraptable = (double[][]) resources
				.getObject("abm.sk.fraction");
		double[] abmrap = new double[abmraptable.length];
		double[] abmwtab = (double[]) resources.getObject("abm.weight");
		// NewBorn-index 1 and in tab is index 0!!
		double abmw = abmwtab[phantom.getIndex() - 1];

		double wold = 1.00;// weighted
		double wnew = 1.00;// weighted

		for (int i = 0; i < energy.length; i++) {
			energy[i] = coeftable[i][0];
			incoh_coeff[i] = coeftable[i][2];
			ph_coeff[i] = coeftable[i][3];
			total_coeff[i] = coeftable[i][5];
			incoh_coeffSk[i] = coeftableSk[i][2];
			ph_coeffSk[i] = coeftableSk[i][3];
			total_coeffSk[i] = coeftableSk[i][5];
			incoh_coeffLu[i] = coeftableLu[i][2];
			ph_coeffLu[i] = coeftableLu[i][3];
			total_coeffLu[i] = coeftableLu[i][5];
			incoh_coeffBr[i] = coeftableBr[i][2];
			ph_coeffBr[i] = coeftableBr[i][3];
			total_coeffBr[i] = coeftableBr[i][5];

			abmrap[i] = abmraptable[i][1];
			// --------------------------------------------------------@@@@@@@RAYLEIGH
			coh_coeff[i] = coeftable[i][1];
			total_coeff_r[i] = coeftable[i][4];
			coh_coeffSk[i] = coeftableSk[i][1];
			total_coeff_rSk[i] = coeftableSk[i][4];
			coh_coeffLu[i] = coeftableLu[i][1];
			total_coeff_rLu[i] = coeftableLu[i][4];
			coh_coeffLu[i] = coeftableLu[i][1];
			total_coeff_rLu[i] = coeftableLu[i][4];
		}
		double abmrapinterp = 0.0;
		double abmtempen = 0.0;
		double abmtemperror = 0.0;
		// -------------------------------------------
		double coh_probab = 0.0;
		double coh_interp = 0.0;// RAYLEIGH=coherent scattering!!
		double total_interp_r = 0.0;
		double rcoh = 0.0;// ALEATOR NUMBER
		// ------------------------------------@@@@@@@@@@@@@@@@@@@@@@@@RAYLEIGH
		// -------MONTE CARLO
		// SIMULATION---------------------------------------------------
		boolean interactB = false;
		for (int j = 0; j < (xrs.getXRayEnergies()).length; j++) {
			// temporary zero init
			for (int k = 0; k < organs.length; k++) {
				endep[k] = 0.0;
				endeperror[k] = 0.0;
				contorInteraction[k] = 0;
			}
			// abm--------------------
			abmtempen = 0.0;
			abmtemperror = 0.0;

			// ----------MC BASED LOOP---------------------------
			for (int i = 0; i < photon_number; i++) {
				interactB = false;
				n_total++;
				e_incident = (xrs.getXRayEnergies())[j] / 1000;// incident_energy
																// in MeV;
				// ----------------------------------------------SC--------------------------------
				oldcoord = getInitialCoord(centerXField, entranceXField,
						projection, phantom, fsd);
				// prelucration----------@@@@@@@@@@----------------
				oldcrd[0] = oldcoord[0];// x
				oldcrd[1] = oldcoord[1];// y
				oldcrd[2] = oldcoord[2];// z
				//theta0 = oldcoord[3];
				//phi0 = oldcoord[4];
				oldcrd[3] = oldcoord[5];
				oldcrd[4] = oldcoord[6];
				oldcrd[5] = oldcoord[7];
				// -------@@@@@@@@@@@@@@@----------------
				// test if we are in body or not
				boolean inbody = phantom.inBody(oldcoord);
				if (inbody)
					isValidInteraction = Phantom.REMAINDER;// at skin entrance
															// SKIN=REMAINDER->REMOVE
				else
					isValidInteraction = Phantom.NO_ORGAN;// exit or miss the
															// phantom
				// -------------------------------------------END
				// SC-------------------------------
				n_photons = 1;

				wold = 1.00;// new hystory
				wnew = 1.00;// new hystory

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
						index_high = index_low;// even if e>0.150MeV!!!
					// --------------choose interaction
					// coefficients----------------------------------
					if (isValidInteraction == Phantom.SKELETON) {
						ph_interp = linInt(energy[index_high],
								ph_coeffSk[index_high], energy[index_low],
								ph_coeffSk[index_low], e_incident);
						//incoh_interp = linInt(energy[index_high],
						//		incoh_coeffSk[index_high], energy[index_low],
						//		incoh_coeffSk[index_low], e_incident);
						//total_interp = linInt(energy[index_high],
						//		total_coeffSk[index_high], energy[index_low],
						//		total_coeffSk[index_low], e_incident);
						dens = densitySk;
						// ABM
						abmrapinterp = linInt(energy[index_high],
								abmrap[index_high], energy[index_low],
								abmrap[index_low], e_incident);
						// @@@@@@@@@@@@@@@@@@@@@@@@@RAYLEIGH=coherent
						// scattering!!
						coh_interp = linInt(energy[index_high],
								coh_coeffSk[index_high], energy[index_low],
								coh_coeffSk[index_low], e_incident);
						total_interp_r = linInt(energy[index_high],
								total_coeff_rSk[index_high], energy[index_low],
								total_coeff_rSk[index_low], e_incident);

					} else if (isValidInteraction == Phantom.LUNGS) {
						ph_interp = linInt(energy[index_high],
								ph_coeffLu[index_high], energy[index_low],
								ph_coeffLu[index_low], e_incident);
						//incoh_interp = linInt(energy[index_high],
						//		incoh_coeffLu[index_high], energy[index_low],
						//		incoh_coeffLu[index_low], e_incident);
						//total_interp = linInt(energy[index_high],
						//		total_coeffLu[index_high], energy[index_low],
						//		total_coeffLu[index_low], e_incident);
						dens = densityLu;
						// @@@@@@@@@@@@@@@@@@@@@@@@@RAYLEIGH=coherent
						// scattering!!
						coh_interp = linInt(energy[index_high],
								coh_coeffLu[index_high], energy[index_low],
								coh_coeffLu[index_low], e_incident);
						total_interp_r = linInt(energy[index_high],
								total_coeff_rLu[index_high], energy[index_low],
								total_coeff_rLu[index_low], e_incident);

					} else if (isValidInteraction == Phantom.BREASTS) {
						ph_interp = linInt(energy[index_high],
								ph_coeffBr[index_high], energy[index_low],
								ph_coeffBr[index_low], e_incident);
						//incoh_interp = linInt(energy[index_high],
						//		incoh_coeffBr[index_high], energy[index_low],
						//		incoh_coeffBr[index_low], e_incident);
						//total_interp = linInt(energy[index_high],
						//		total_coeffBr[index_high], energy[index_low],
						//		total_coeffBr[index_low], e_incident);
						dens = densityBr;
						// @@@@@@@@@@@@@@@@@@@@@@@@@RAYLEIGH=coherent
						// scattering!!
						coh_interp = linInt(energy[index_high],
								coh_coeffBr[index_high], energy[index_low],
								coh_coeffBr[index_low], e_incident);
						total_interp_r = linInt(energy[index_high],
								total_coeff_rBr[index_high], energy[index_low],
								total_coeff_rBr[index_low], e_incident);

					} else// soft tissue-common or noorgan but this is later
							// resolved
					{
						ph_interp = linInt(energy[index_high],
								ph_coeff[index_high], energy[index_low],
								ph_coeff[index_low], e_incident);
						//incoh_interp = linInt(energy[index_high],
						//		incoh_coeff[index_high], energy[index_low],
						//		incoh_coeff[index_low], e_incident);
						//total_interp = linInt(energy[index_high],
						//		total_coeff[index_high], energy[index_low],
						//		total_coeff[index_low], e_incident);
						dens = density;
						// @@@@@@@@@@@@@@@@@@@@@@@@@RAYLEIGH=coherent
						// scattering!!
						coh_interp = linInt(energy[index_high],
								coh_coeff[index_high], energy[index_low],
								coh_coeff[index_low], e_incident);
						total_interp_r = linInt(energy[index_high],
								total_coeff_r[index_high], energy[index_low],
								total_coeff_r[index_low], e_incident);
					}
					// -------------------------------------------
					coh_probab = (coh_interp) / (total_interp_r - ph_interp);// coh
																				// from
																				// REMAINING
					ph_probab = ph_interp / total_interp_r;// absortion from
															// total
					// --------------------------------------------------------------------------
					// --END
					// COEFF--------------------------------------------------------------------
					//r = RandomCollection.random01();
					double rdist = RandomCollection.random01();
					while (rdist == 0.0)
						rdist = RandomCollection.random01();
					dist_to_int = -Math.log(rdist) / (dens * total_interp_r);// [cm]!!@
					if (isValidInteraction != Phantom.NO_ORGAN)// ---------!!!!!!
					{
						interactB = true;
						wnew = wold * (1 - ph_probab);// probab of survive from
														// photoel absortion
						if (wnew < 0.003 || e_incident < 2.000E-03)// final
																	// absortion
						{
							n_photoefect++;
							contorInteraction[isValidInteraction - 1]++;// interaction
																		// colission
							contorInteraction_total[isValidInteraction - 1] = contorInteraction_total[isValidInteraction - 1] + 1.0;
							endep[isValidInteraction - 1] = endep[isValidInteraction - 1]
									+ e_incident;
							//enmedtemp = endep[isValidInteraction - 1]
							//		/ contorInteraction[isValidInteraction - 1];
							enindivtemp = e_incident;
							// endeperror[isValidInteraction-1]=endeperror[isValidInteraction-1]+
							// Math.pow(enindivtemp-enmedtemp,2);
							endeperror[isValidInteraction - 1] = endeperror[isValidInteraction - 1]
									+ enindivtemp * enindivtemp;
							// separated ABM-----for skeleton
							// eval--------------------------------
							if (isValidInteraction == Phantom.SKELETON) {
								abmtempen = abmtempen + abmrapinterp * abmw
										* e_incident;
								// abmtemperror=abmtemperror+
								// Math.pow(abmrapinterp*abmw*enindivtemp-abmtempen/contorInteraction[isValidInteraction-1],2);
								abmtemperror = abmtemperror + abmrapinterp
										* abmw * enindivtemp * abmrapinterp
										* abmw * enindivtemp;
							}
							n_photons = 0;
						} else// scatter
						{
							rcoh = RandomCollection.random01();
							if (rcoh <= coh_probab) {
								theta = rayleighSim(e_incident);
								r0 = RandomCollection.random01();// [0-1]
								phi = 2 * r0 * Math.PI;// [0.2PI];
								newcoord = getUpdatedCoordinate(theta, phi,
										dist_to_int, oldcrd);
								isValidInteraction = phantom
										.inWhatOrgan(newcoord);// HERE->organism
								// permute newcoordonate with old coordonate
								oldcrd[0] = newcoord[0];
								oldcrd[1] = newcoord[1];
								oldcrd[2] = newcoord[2];
								oldcrd[3] = newcoord[3];
								oldcrd[4] = newcoord[4];
								oldcrd[5] = newcoord[5];
							} else {
								contorInteraction[isValidInteraction - 1]++;// interaction
																			// colission
								contorInteraction_total[isValidInteraction - 1] = contorInteraction_total[isValidInteraction - 1] + 1.0;
								theta = comptonSim(e_incident);
								e_scatt = e_incident
										/ (1 + (e_incident / 0.511)
												* (1 - Math.cos(theta)));
								endep[isValidInteraction - 1] = endep[isValidInteraction - 1]
										+ wnew
										* ph_probab
										* e_incident
										+ wnew
										* (1 - ph_probab)
										* (e_incident - e_scatt);
								//enmedtemp = endep[isValidInteraction - 1]
								//		/ contorInteraction[isValidInteraction - 1];
								enindivtemp = wnew * ph_probab * e_incident
										+ wnew * (1 - ph_probab)
										* (e_incident - e_scatt);
								// endeperror[isValidInteraction-1]=endeperror[isValidInteraction-1]+
								// Math.pow(enindivtemp-enmedtemp,2);
								endeperror[isValidInteraction - 1] = endeperror[isValidInteraction - 1]
										+ enindivtemp * enindivtemp;
								// separated ABM-----for skeleton
								// eval--------------------------------
								if (isValidInteraction == Phantom.SKELETON) {
									abmtempen = abmtempen + abmrapinterp * abmw
											* (enindivtemp);
									// abmtemperror=abmtemperror+
									// Math.pow(enindivtemp*abmrapinterp*abmw-abmtempen/contorInteraction[isValidInteraction-1],2);
									abmtemperror = abmtemperror + abmrapinterp
											* abmw * enindivtemp * abmrapinterp
											* abmw * enindivtemp;
								}
								// continue tracing photon hystory
								e_incident = e_scatt;// new loop in while
								// azimutal angle
								r0 = RandomCollection.random01();// [0-1]
								phi = 2 * r0 * Math.PI;// [0.2PI];
								// ----------------------------------------RAD---------------
								newcoord = getUpdatedCoordinate(theta, phi,
										dist_to_int, oldcrd);
								// -------------------!!!!!!!!!_____________________________
								isValidInteraction = phantom
										.inWhatOrgan(newcoord);// retrieve the
																// organ
								// HERE->organisminclude inBody!
								// permute newcoordonate with old coordonate
								oldcrd[0] = newcoord[0];
								oldcrd[1] = newcoord[1];
								oldcrd[2] = newcoord[2];
								oldcrd[3] = newcoord[3];
								oldcrd[4] = newcoord[4];
								oldcrd[5] = newcoord[5];

								wold = wnew;
							}
						}// scatter
					} else// NoOrgan
					{
						// pathlength -->exit the phantom
						n_photons = 0;// exit region of interest, so we kill the
										// photon
						n_detected++;
						if (projection == 0)// corect are loc initial
											// mulajul----legs=virtual=cu trunk
						{
							if (oldcrd[1] > 0.0)// phantom.getOrganProjThickness(projection,oldcrd[2])/2)
							{
								if (!interactB)
									n_film_direct++;

								n_film++;
							}
						}
						if (projection == 2) {
							if (oldcrd[1] < 0.0)// -phantom.getOrganProjThickness(projection,oldcrd[2])/2)
							{
								if (!interactB)
									n_film_direct++;

								n_film++;
							}
						}
						if (projection == 1) {
							if (oldcrd[0] > 0.0)// phantom.getOrganProjThickness(projection,oldcrd[2])/2)
							{
								if (!interactB)
									n_film_direct++;

								n_film++;
							}
						}
						if (projection == 3) {
							if (oldcrd[0] < 0.0)// -phantom.getOrganProjThickness(projection,oldcrd[2])/2)
							{
								if (!interactB)
									n_film_direct++;

								n_film++;
							}
						}

						n_side = n_detected - n_film;
					}
				}// end while
			}// end for i=0 to photon_number
				// NOW WE HAVE DISTRIBUTION OF ENERGY IN ALL ORGANS FOR e.g
				// 10000 phgoton hystories.
				// error and we try the follow code.MAYBE DUE TO integer
				// multiplication !!!!!!!!!!!@@@@@@@@!!
			for (int k = 0; k < organs.length; k++) {
				/*
				 * if (contorInteraction[k]>0) { if (contorInteraction[k]>1) {
				 * //standard deviation of mean endeperror[k]=
				 * Math.sqrt((endeperror
				 * [k]/contorInteraction[k])/(contorInteraction[k]-1)); }//else
				 * is unchanged!! }//else unchanged
				 */
				// ---REAL @@ normalized per
				// energy-------------------------------------
				// endeperror[k]=endeperror[k]/photon_number;
				// scale spectrum
				// endeperror[k]=endeperror[k]*(xrs.getXRayIntensities())[j]/xrs.getNormalizedValue();
				// REAL @@ normalized per
				// energy-------------------------------------
				endep[k] = endep[k] / photon_number;
				// @@@@@@@@@@@@
				endeperror[k] = endeperror[k] / photon_number;
				endeperror[k] = (endeperror[k] - endep[k] * endep[k])
						/ (photon_number - 1.0);
				if (endeperror[k] >= 0.0)
					endeperror[k] = Math.sqrt(endeperror[k]);
				endeperror[k] = endeperror[k] * (xrs.getXRayIntensities())[j]
						/ xrs.getNormalizedValue();
				// @@@@@@@@@@@@@
				// ----------------------------------------------SCALE
				// SPECTRUM---------------
				endep[k] = endep[k] * (xrs.getXRayIntensities())[j]
						/ xrs.getNormalizedValue();
				// ----------------NOW REAL
				// VARIABLES------------------------------------------
				energy_deposit[k][0] = energy_deposit[k][0] + endep[k];
				energy_deposit[k][1] = Math.sqrt(energy_deposit[k][1]
						* energy_deposit[k][1] + endeperror[k] * endeperror[k]);
			}
			// separated ABM--from skeleton eval
			/*
			 * if (contorInteraction[phantom.SKELETON-1]>0) { if
			 * (contorInteraction[phantom.SKELETON-1]>1) { //standard deviation
			 * of mean abmtemperror=
			 * Math.sqrt((abmtemperror/contorInteraction[phantom
			 * .SKELETON-1])/(contorInteraction[phantom.SKELETON-1]-1)); }//else
			 * is unchanged!! }//else unchanged
			 */
			// ---REAL @@ normalized per
			// energy-------------------------------------
			// abmtemperror=abmtemperror/photon_number;
			// scale spectrum
			// abmtemperror=abmtemperror*(xrs.getXRayIntensities())[j]/xrs.getNormalizedValue();
			// ---REAL @@ normalized per
			// energy-------------------------------------
			abmtempen = abmtempen / photon_number;
			abmtemperror = abmtemperror / photon_number;
			abmtemperror = (abmtemperror - abmtempen * abmtempen)
					/ (photon_number - 1.0);
			if (abmtemperror >= 0.0)
				abmtemperror = Math.sqrt(abmtemperror);
			abmtemperror = abmtemperror * (xrs.getXRayIntensities())[j]
					/ xrs.getNormalizedValue();

			// ----------------------------------------------SCALE
			// SPECTRUM---------------
			abmtempen = abmtempen * (xrs.getXRayIntensities())[j]
					/ xrs.getNormalizedValue();
			// ----------------NOW REAL VARIABLES IN ABM NOT IN SKELETON-DIFF OF
			// 1!!!------------
			energy_deposit[Phantom.ACTVEBONEMARROW - 1][0] = energy_deposit[Phantom.ACTVEBONEMARROW - 1][0]
					+ abmtempen;
			energy_deposit[Phantom.ACTVEBONEMARROW - 1][1] = Math
					.sqrt(energy_deposit[Phantom.ACTVEBONEMARROW - 1][1]
							* energy_deposit[Phantom.ACTVEBONEMARROW - 1][1]
							+ abmtemperror * abmtemperror);
			// ---------------------
		}// end energy loop

		double exposure = airExposure * 1000;// -mGy->microGy
		double mevtojoule = 1.60218E-13;// 1MeV=1,60218 � 10-13 J
		double factor = transversalArea * exposure * mevtojoule
				* xrs.getPhotonFlux() / xrs.getAirKerma();
		// separated ABM
		if (energy_deposit[Phantom.ACTVEBONEMARROW - 1][0] != 0.0)
			energy_deposit[Phantom.ACTVEBONEMARROW - 1][1] = Math.min(2 * 100
					* energy_deposit[Phantom.ACTVEBONEMARROW - 1][1]
					/ energy_deposit[Phantom.ACTVEBONEMARROW - 1][0], 99.9);// %
		else
			energy_deposit[Phantom.ACTVEBONEMARROW - 1][1] = 99.9;// 0.0;
			// normalized at real XRay Spectrum:
		energy_deposit[Phantom.ACTVEBONEMARROW - 1][0] = energy_deposit[Phantom.ACTVEBONEMARROW - 1][0]
				* factor;// Jouli
		// organ dose:
		organAbsorbedDose[Phantom.ACTVEBONEMARROW - 1][0] = 1000 * 1000
				* energy_deposit[Phantom.ACTVEBONEMARROW - 1][0]
				/ phantom.getOrganMass(Phantom.ACTVEBONEMARROW - 1);
		// 1g=10-03kk,Gy->mGy
		organAbsorbedDose[Phantom.ACTVEBONEMARROW - 1][1] = energy_deposit[Phantom.ACTVEBONEMARROW - 1][1];// %
		// -------------------------------
		// 33333333333333333333333333333333333333333333333
		for (int k = 0; k < organs.length; k++) {
			if (contorInteraction_total[k] > 0.0)
				organAbsorbedDose[k][1] = Math.min(
						2.0 * 100.0 * Math.sqrt(contorInteraction_total[k])
								/ contorInteraction_total[k], 99.9);
			else
				organAbsorbedDose[k][1] = 99.9;
		}
		organAbsorbedDose[Phantom.ACTVEBONEMARROW - 1][1] = organAbsorbedDose[Phantom.SKELETON - 1][1];
		// 33333333333333333333333333333333333333333333333

		double[] wt = ((double[]) resources.getObject("wt.organname"));
		effectiveDose[0] = 0.0;
		double seff2 = 0.0;
		double errmin = organAbsorbedDose[Phantom.ACTVEBONEMARROW - 1][1];// initialisation
		for (int k = 0; k < organs.length; k++) {
			if (k != Phantom.ACTVEBONEMARROW - 1) {
				// //% and 2 here for infinit grdlib (Nmare) is 95% conf level
				if (energy_deposit[k][0] != 0.0)
					energy_deposit[k][1] = Math
							.min(2 * 100 * energy_deposit[k][1]
									/ energy_deposit[k][0], 99.9);// %
				else
					energy_deposit[k][1] = 99.9;// 0.0;
				// normalized at real XRay Spectrum:
				energy_deposit[k][0] = energy_deposit[k][0] * factor;// Jouli
				// organ dose:
				organAbsorbedDose[k][0] = 1000 * 1000 * energy_deposit[k][0]
						/ phantom.getOrganMass(k);
				// 1g=10-03kk,Gy->mGy
				// organAbsorbedDose[k][1]=energy_deposit[k][1];//%
				errmin = Math.min(errmin, organAbsorbedDose[k][1]);
			}
			effectiveDose[0] = effectiveDose[0] + wt[k]
					* organAbsorbedDose[k][0];
			seff2 = seff2
					+ Math.pow(wt[k] * organAbsorbedDose[k][0]
							* organAbsorbedDose[k][1] / 100, 2);
		}
		if (effectiveDose[0] != 0.0)
			effectiveDose[1] = errmin;// Math.min((100*Math.sqrt(seff2))/effectiveDose[0],99.9);//%
		else
			effectiveDose[1] = 99.9;// 0.0;

		simulationTimeElapsed = timeElapsed(startSimulationTime);
	}

	// --------------------CT-----------------------------------------------------------------------
	// general spectrum
	// here airExposure=CTDI(mGy),fsd=focus axe distance(cm)
	// startSlice z1(cm), endSlice z2(cm), nominal slice thickness(cm)
	// nomSliceThick should be always equal with slice thickness for CTDI
	// slicePerStep SHOULB BE ALWAYS 1, otherwise->dose is overestimated due to
	// high exposure!!!
	/**
	 * Perform Monte Carlo simulation for CT. Rotation in CT is roughly simulated by successive "radiographic" runs for each 4 
	 * projection types (AP, PA, LLAT, RLAT).
	 * 
	 * @param photon_number the number of histories for each photon energy and each projection type.
	 * @param xrs the link to XRaySpectrum class.
	 * @param airExposure the dose free in air per rotation measured at central axis (CTDI) in mGy.
	 * @param fsd focus-central axis distance (i.e. focus-where CTDI i measured distance) in cm.
	 * @param phantom the link to geometry, i.e. the MIRD5 phantom
	 * @param startSlice the start position in cm (z1) relative to phantom coordinate system
	 * @param endSlice the end position in cm (z2) relative to phantom coordinate system
	 * @param nomSliceThick the slice thickness in cm
	 * @param slicePerStep this SHOULD NORMALLY BE 1, it is related to how many slice per rotation are covered (kind of pitch factor)!
	 */
	public static void ctSim(int photon_number, XRaySpectrum xrs,
			double airExposure, double fsd,
			Phantom phantom, // int projection,
			// double[] centerXField, double[] entranceXField,
			double startSlice, double endSlice, double nomSliceThick,
			int slicePerStep) {
		long startSimulationTime = System.currentTimeMillis();
		simulationTimeElapsed = "";
		int index_low = 0;
		int index_high = 0;
		double ph_interp = 0.0;
		//double incoh_interp = 0.0;// Compton=incoherent scattering!!
		//double total_interp = 0.0;
		double ph_probab = 0.0;
		double r = 0.0;// [0,1] random number
		double dist_to_int = 0.0;// store distance to interaction
		// -------------SC-------------------------------------------------------------
		double r0 = 0.0;
		double[] oldcoord = new double[8];
		// -----------------------------
		double[] oldcrd = new double[6];
		//double theta0 = 0.0;
		//double phi0 = 0.0;
		// ---------------------
		double phi = 0.0;
		double[] newcoord;
		int isValidInteraction = Phantom.NO_ORGAN;// INIT ANYWAY
		// ----------------END
		// SC-----------------------------------------------------
		// Area in mm2
		double transversalArea = 0.0;// getEntranceField(fsd,entranceXField,centerXField,projection,phantom)*100;
		// -----------
		double e_incident = 0.0;// init
		double e_scatt = 0.0;// store the energy of the Compton scattered
								// photons
		// -----------INIT ORGANS BASED VARIABLES!!---------------------------
		organs = (String[]) resources.getObject("rad.organname");
		energy_deposit = new double[organs.length][2];// all organs
		organAbsorbedDose = new double[organs.length][2];// all organs
		orgAbsorbedDose = new double[organs.length][2];// all
														// organs---------PROJ
		int[] contorInteraction = new int[organs.length];
		double[] contorInteraction_total = new double[organs.length];
		double[] endep = new double[organs.length];// temporary
		double[] endeperror = new double[organs.length];// temporary
		for (int i = 0; i < organs.length; i++)// just in case
		{
			energy_deposit[i][0] = 0.0;
			energy_deposit[i][1] = 0.0;
			organAbsorbedDose[i][0] = 0.0;
			organAbsorbedDose[i][1] = 0.0;
			// ------------------------
			endep[i] = 0.0;
			endeperror[i] = 0.0;
			contorInteraction[i] = 0;
			contorInteraction_total[i] = 0.0;
			// ---------------------------------------------------------PROJ
			orgAbsorbedDose[i][0] = 0.0;
			orgAbsorbedDose[i][1] = 0.0;
		}
		// ---------------------------------
		//double enmedtemp = 0.0;// handle error eval
		double enindivtemp = 0.0;// handle error eval
		// ----------------------------------
		effectiveDose = new double[2];
		effectiveDose[0] = 0.0;
		effectiveDose[1] = 0.0;
		n_detected = 0;
		n_total = 0;
		n_film = 0;
		n_film_direct = 0;
		n_side = 0;
		n_photoefect = 0;
		double theta = 0.0;// store the Compton simulation's scattering angle
		int n_photons = 0;// number of photons being transported
		// ------read attenuation coefficient and density----------
		double[][] coeftable = new double[0][0];
		double[][] coeftableSk = new double[0][0];
		double[][] coeftableLu = new double[0][0];
		double[][] coeftableBr = (double[][]) resources
				.getObject("breast.diagnostic.attenuationCoef");
		double density = 0.0;
		double densitySk = 0.0;
		double densityLu = 0.0;
		double densityBr = ((Double) resources
				.getObject("breast.diagnostic.density")).doubleValue();

		if (phantom.getIndex() != Phantom.NEWBORN_INDEX) {
			coeftable = (double[][]) resources
					.getObject("commonTissue.diagnostic.attenuationCoef");
			coeftableSk = (double[][]) resources
					.getObject("skeletal_boneWithMarrow.diagnostic.attenuationCoef");
			coeftableLu = (double[][]) resources
					.getObject("lung.diagnostic.attenuationCoef");
			density = ((Double) resources
					.getObject("commonTissue.diagnostic.density"))
					.doubleValue();
			densitySk = ((Double) resources
					.getObject("skeletal_boneWithMarrow.diagnostic.density"))
					.doubleValue();
			densityLu = ((Double) resources
					.getObject("lung.diagnostic.density")).doubleValue();

		} else // MAMO WILL NEVER HAPPEN HERE!!
		{
			coeftable = (double[][]) resources
					.getObject("commonTissue.newborn.diagnostic.attenuationCoef");
			coeftableSk = (double[][]) resources
					.getObject("skeletal_boneWithMarrow.newborn.diagnostic.attenuationCoef");
			coeftableLu = (double[][]) resources
					.getObject("lung.newborn.diagnostic.attenuationCoef");
			density = ((Double) resources
					.getObject("commonTissue.newborn.diagnostic.density"))
					.doubleValue();
			densitySk = ((Double) resources
					.getObject("skeletal_boneWithMarrow.newborn.diagnostic.density"))
					.doubleValue();
			densityLu = ((Double) resources
					.getObject("lung.newborn.diagnostic.density"))
					.doubleValue();
		}
		double dens = 0.0;// handle densities
		// coeftable[i][0])-->energy in MeV;
		// coeftable[i][2])-->Compton att. coef in cm2/g;
		// coeftable[i][3])-->Photoelectric abs. coef in cm2/g;
		// coeftable[i][5])-->Total (without Rayleigh scattering) att. coef in
		// cm2/g;
		// read data:
		double[] energy = new double[coeftable.length];
		double[] ph_coeff = new double[coeftable.length];
		double[] incoh_coeff = new double[coeftable.length];
		double[] total_coeff = new double[coeftable.length];
		double[] ph_coeffSk = new double[coeftableSk.length];// same energy!!
		double[] incoh_coeffSk = new double[coeftableSk.length];
		double[] total_coeffSk = new double[coeftableSk.length];
		double[] ph_coeffLu = new double[coeftableLu.length];
		double[] incoh_coeffLu = new double[coeftableLu.length];
		double[] total_coeffLu = new double[coeftableLu.length];
		double[] ph_coeffBr = new double[coeftableBr.length];
		double[] incoh_coeffBr = new double[coeftableBr.length];
		double[] total_coeffBr = new double[coeftableBr.length];
		// ----------------------------------------------------------------------
		double[] coh_coeff = new double[coeftable.length];
		double[] total_coeff_r = new double[coeftable.length];// @@@@@@@RAYLEIGH
		double[] coh_coeffSk = new double[coeftableSk.length];
		double[] total_coeff_rSk = new double[coeftableSk.length];// @@@@@@@RAYLEIGH
		double[] coh_coeffLu = new double[coeftableLu.length];
		double[] total_coeff_rLu = new double[coeftableLu.length];// @@@@@@@RAYLEIGH
		double[] coh_coeffBr = new double[coeftableBr.length];
		double[] total_coeff_rBr = new double[coeftableBr.length];// @@@@@@@RAYLEIGH

		// separated ABM data
		double[][] abmraptable = (double[][]) resources
				.getObject("abm.sk.fraction");
		double[] abmrap = new double[abmraptable.length];
		double[] abmwtab = (double[]) resources.getObject("abm.weight");
		// NewBorn-index 1 and in tab is index 0!!
		double abmw = abmwtab[phantom.getIndex() - 1];

		for (int i = 0; i < energy.length; i++) {
			energy[i] = coeftable[i][0];
			incoh_coeff[i] = coeftable[i][2];
			ph_coeff[i] = coeftable[i][3];
			total_coeff[i] = coeftable[i][5];
			incoh_coeffSk[i] = coeftableSk[i][2];
			ph_coeffSk[i] = coeftableSk[i][3];
			total_coeffSk[i] = coeftableSk[i][5];
			incoh_coeffLu[i] = coeftableLu[i][2];
			ph_coeffLu[i] = coeftableLu[i][3];
			total_coeffLu[i] = coeftableLu[i][5];
			incoh_coeffBr[i] = coeftableBr[i][2];
			ph_coeffBr[i] = coeftableBr[i][3];
			total_coeffBr[i] = coeftableBr[i][5];

			abmrap[i] = abmraptable[i][1];
			// --------------------------------------------------------@@@@@@@RAYLEIGH
			coh_coeff[i] = coeftable[i][1];
			total_coeff_r[i] = coeftable[i][4];
			coh_coeffSk[i] = coeftableSk[i][1];
			total_coeff_rSk[i] = coeftableSk[i][4];
			coh_coeffLu[i] = coeftableLu[i][1];
			total_coeff_rLu[i] = coeftableLu[i][4];
			coh_coeffLu[i] = coeftableLu[i][1];
			total_coeff_rLu[i] = coeftableLu[i][4];
		}
		double abmrapinterp = 0.0;
		double abmtempen = 0.0;
		double abmtemperror = 0.0;
		// -------------------------------------------
		double coh_probab = 0.0;
		double coh_interp = 0.0;// RAYLEIGH=coherent scattering!!
		double total_interp_r = 0.0;
		double rcoh = 0.0;// ALEATOR NUMBER
		// ------------------------------------@@@@@@@@@@@@@@@@@@@@@@@@RAYLEIGH
		// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@CT
		double exposure = 0.0;// airExposure*1000;//-mGy->microGy
		double mevtojoule = 1.60218E-13;// 1MeV=1,60218 � 10-13 J
		double factor = 0.0;
		double fsdd = 0.0;// real fsd
		double[] entXField = new double[2];
		double[] cenXField = new double[3];
		double vrtx = 0.0;
		double zet = 0.0;
		double startS = 0.0;// startSlice;
		double endS = 0.0;// endSlice;
		// -------MONTE CARLO
		// SIMULATION---------------------------------------------------
		boolean interactB = false;
		// int projection=0;//0-AP,2-PA,1-LL(RL),3-RL(LL)
		int currentscannedslice = 0;
		scannedslices = new Double(Math.abs(startSlice - endSlice)
				/ nomSliceThick).intValue();
		startS = startSlice;// initial
		while (currentscannedslice < scannedslices) {
			for (int projection = 0; projection < 4; projection++) {
				// temporary zero init
				for (int k = 0; k < organs.length; k++) {
					energy_deposit[k][0] = 0.0;
					energy_deposit[k][1] = 0.0;
					organAbsorbedDose[k][0] = 0.0;
					organAbsorbedDose[k][1] = 0.0;
				}
				// construct the new center and XYField from start &end slice
				// ventrices!!
				cenXField[0] = 0.0;// xcenter
				cenXField[1] = 0.0;// ycenter
				if (slicePerStep != 0)
					vrtx = startS - slicePerStep * nomSliceThick;// @@@zet is
																	// !!!
				else
					vrtx = endSlice;
				endS = vrtx;// zet
				cenXField[2] = phantom.getVertexCoordonate() - (startS + endS)
						/ 2;// zcenter@@@@@@@
				// direct luat
				// zet=phantom.getVertexCoordonate()-cenXField[2];
				zet = (startS + endS) / 2;// zetcenter
				entXField[0] = phantom.getOrganProjWidth(projection, zet);
				if (slicePerStep != 0)
					entXField[1] = slicePerStep * nomSliceThick;
				else
					entXField[1] = scannedslices * nomSliceThick;// entranceXField[1];
					// here fsd means focus-axe distance
				fsdd = fsd - phantom.getOrganProjThickness(projection, zet) / 2;
				// transversalArea=getEntranceField(fsd,entranceXField,centerXField,projection,phantom)*100;
				transversalArea = getEntranceField(fsdd, entXField, cenXField,
						projection, phantom) * 100;
				for (int j = 0; j < (xrs.getXRayEnergies()).length; j++) {
					// temporary zero init
					for (int k = 0; k < organs.length; k++) {
						endep[k] = 0.0;
						endeperror[k] = 0.0;
						contorInteraction[k] = 0;
					}
					// abm--------------------
					abmtempen = 0.0;
					abmtemperror = 0.0;
					// ----------MC BASED LOOP---------------------------
					for (int i = 0; i < photon_number; i++) {
						interactB = false;
						n_total++;
						e_incident = (xrs.getXRayEnergies())[j] / 1000;// incident_energy
																		// in
																		// MeV;
						// ----------------------------------------------SC--------------------------------
						// oldcoord=getInitialCoord(centerXField,entranceXField,projection,phantom,fsd);
						oldcoord = getInitialCoord(cenXField, entXField,
								projection, phantom, fsdd);
						// prelucration----------@@@@@@@@@@----------------
						oldcrd[0] = oldcoord[0];// x
						oldcrd[1] = oldcoord[1];// y
						oldcrd[2] = oldcoord[2];// z
						//theta0 = oldcoord[3];
						//phi0 = oldcoord[4];
						oldcrd[3] = oldcoord[5];// ux
						oldcrd[4] = oldcoord[6];// uy
						oldcrd[5] = oldcoord[7];// uz
						// -------@@@@@@@@@@@@@@@----------------
						// test if we are in body or not
						boolean inbody = phantom.inBody(oldcoord);
						if (inbody)
							isValidInteraction = Phantom.REMAINDER;// at skin
																	// entrance
																	// SKIN=REMAINDER->REMOVE
						else
							isValidInteraction = Phantom.NO_ORGAN;// exit or
																	// miss the
																	// phantom

						// -------------------------------------------END
						// SC-------------------------------
						n_photons = 1;
						while (n_photons == 1) {
							if (stopB)
								return;

							// evaluate
							// ph_probab,ph_interp,compton_interp,total_interp
							// for given energy
							Sort.findNearestValue(energy, e_incident, true);
							index_low = Sort.getNearestPosition();
							if (index_low < energy.length - 1)
								index_high = index_low + 1;
							else
								index_high = index_low;// even if e>0.150MeV!!!
							// --------------choose interaction
							// coefficients----------------------------------
							if (isValidInteraction == Phantom.SKELETON) {
								ph_interp = linInt(energy[index_high],
										ph_coeffSk[index_high],
										energy[index_low],
										ph_coeffSk[index_low], e_incident);
								//incoh_interp = linInt(energy[index_high],
								//		incoh_coeffSk[index_high],
								//		energy[index_low],
								//		incoh_coeffSk[index_low], e_incident);
								//total_interp = linInt(energy[index_high],
								//		total_coeffSk[index_high],
								//		energy[index_low],
								//		total_coeffSk[index_low], e_incident);
								dens = densitySk;
								// ABM
								abmrapinterp = linInt(energy[index_high],
										abmrap[index_high], energy[index_low],
										abmrap[index_low], e_incident);
								// @@@@@@@@@@@@@@@@@@@@@@@@@RAYLEIGH=coherent
								// scattering!!
								coh_interp = linInt(energy[index_high],
										coh_coeffSk[index_high],
										energy[index_low],
										coh_coeffSk[index_low], e_incident);
								total_interp_r = linInt(energy[index_high],
										total_coeff_rSk[index_high],
										energy[index_low],
										total_coeff_rSk[index_low], e_incident);

							} else if (isValidInteraction == Phantom.LUNGS) {
								ph_interp = linInt(energy[index_high],
										ph_coeffLu[index_high],
										energy[index_low],
										ph_coeffLu[index_low], e_incident);
								//incoh_interp = linInt(energy[index_high],
								//		incoh_coeffLu[index_high],
								//		energy[index_low],
								//		incoh_coeffLu[index_low], e_incident);
								//total_interp = linInt(energy[index_high],
								//		total_coeffLu[index_high],
								//		energy[index_low],
									//	total_coeffLu[index_low], e_incident);
								dens = densityLu;
								// @@@@@@@@@@@@@@@@@@@@@@@@@RAYLEIGH=coherent
								// scattering!!
								coh_interp = linInt(energy[index_high],
										coh_coeffLu[index_high],
										energy[index_low],
										coh_coeffLu[index_low], e_incident);
								total_interp_r = linInt(energy[index_high],
										total_coeff_rLu[index_high],
										energy[index_low],
										total_coeff_rLu[index_low], e_incident);

							} else if (isValidInteraction == Phantom.BREASTS) {
								ph_interp = linInt(energy[index_high],
										ph_coeffBr[index_high],
										energy[index_low],
										ph_coeffBr[index_low], e_incident);
								//incoh_interp = linInt(energy[index_high],
								//		incoh_coeffBr[index_high],
								//		energy[index_low],
								//		incoh_coeffBr[index_low], e_incident);
								//total_interp = linInt(energy[index_high],
								//		total_coeffBr[index_high],
								//		energy[index_low],
								//		total_coeffBr[index_low], e_incident);
								dens = densityBr;
								// @@@@@@@@@@@@@@@@@@@@@@@@@RAYLEIGH=coherent
								// scattering!!
								coh_interp = linInt(energy[index_high],
										coh_coeffBr[index_high],
										energy[index_low],
										coh_coeffBr[index_low], e_incident);
								total_interp_r = linInt(energy[index_high],
										total_coeff_rBr[index_high],
										energy[index_low],
										total_coeff_rBr[index_low], e_incident);

							} else// soft tissue-common or noorgan but this is
									// later resolved
							{
								ph_interp = linInt(energy[index_high],
										ph_coeff[index_high],
										energy[index_low], ph_coeff[index_low],
										e_incident);
								//incoh_interp = linInt(energy[index_high],
								//		incoh_coeff[index_high],
								//		energy[index_low],
								//		incoh_coeff[index_low], e_incident);
								//total_interp = linInt(energy[index_high],
								//		total_coeff[index_high],
								//		energy[index_low],
								//		total_coeff[index_low], e_incident);
								dens = density;
								// @@@@@@@@@@@@@@@@@@@@@@@@@RAYLEIGH=coherent
								// scattering!!
								coh_interp = linInt(energy[index_high],
										coh_coeff[index_high],
										energy[index_low],
										coh_coeff[index_low], e_incident);
								total_interp_r = linInt(energy[index_high],
										total_coeff_r[index_high],
										energy[index_low],
										total_coeff_r[index_low], e_incident);
							}
							// ph_probab=ph_interp/total_interp;
							// -------------------------------------------
							coh_probab = (coh_interp)
									/ (total_interp_r - ph_interp);// coh from
																	// REMAINING
							ph_probab = ph_interp / total_interp_r;// absortion
																	// from
																	// total
							// --------------------------------------------------------------------------
							// --END
							// COEFF--------------------------------------------------------------------
							r = RandomCollection.random01();
							double rdist = RandomCollection.random01();
							while (rdist == 0.0)
								rdist = RandomCollection.random01();
							dist_to_int = -Math.log(rdist)
									/ (dens * total_interp_r);// [cm]!!@

							if (isValidInteraction != Phantom.NO_ORGAN)// ---------!!!!!!
							{
								interactB = true;
								r = RandomCollection.random01();
								if (r <= ph_probab) // photoelectric interaction
								{
									contorInteraction[isValidInteraction - 1]++;// interaction
																				// colission
									contorInteraction_total[isValidInteraction - 1] = contorInteraction_total[isValidInteraction - 1] + 1.0;
									n_photoefect++;
									endep[isValidInteraction - 1] = endep[isValidInteraction - 1]
											+ e_incident;// --!!!!!!!
									// from phantom if
									// isValidInteraction=phantom.REMAINDER=25->
									// in organs coresponds 24!!
									//enmedtemp = endep[isValidInteraction - 1]
									//		/ contorInteraction[isValidInteraction - 1];
									// per photon, endep is sum
									enindivtemp = e_incident;
									// endeperror[isValidInteraction-1]=endeperror[isValidInteraction-1]+
									// Math.pow(enindivtemp-enmedtemp,2);//error
									// at power 2
									endeperror[isValidInteraction - 1] = endeperror[isValidInteraction - 1]
											+ enindivtemp * enindivtemp;
									// separated ABM-----for skeleton
									// eval---------------------------------
									if (isValidInteraction == Phantom.SKELETON) {
										abmtempen = abmtempen + abmrapinterp
												* abmw * enindivtemp;
										// abmtemperror=abmtemperror+
										// Math.pow(abmrapinterp*abmw*enindivtemp-abmtempen/contorInteraction[isValidInteraction-1],2);
										abmtemperror = abmtemperror
												+ abmrapinterp * abmw
												* enindivtemp * abmrapinterp
												* abmw * enindivtemp;
									}
									// So=>absorbed photons->end hystory!!
									n_photons = 0;
								} else // Compton scattering
								{
									// --------RAYLEIGH--elastic
									// COHERENT----------
									rcoh = RandomCollection.random01();
									if (rcoh <= coh_probab) {
										theta = rayleighSim(e_incident);
										r0 = RandomCollection.random01();// [0-1]
										phi = 2 * r0 * Math.PI;// ->[0.2PI];
										newcoord = getUpdatedCoordinate(theta,
												phi, dist_to_int, oldcrd);
										isValidInteraction = phantom
												.inWhatOrgan(newcoord);// HERE->organism
										// permute newcoordonate with old
										// coordonate
										oldcrd[0] = newcoord[0];
										oldcrd[1] = newcoord[1];
										oldcrd[2] = newcoord[2];
										oldcrd[3] = newcoord[3];
										oldcrd[4] = newcoord[4];
										oldcrd[5] = newcoord[5];
									} else// COMPTON
									{
										contorInteraction[isValidInteraction - 1]++;// interaction
																					// colission
										contorInteraction_total[isValidInteraction - 1] = contorInteraction_total[isValidInteraction - 1] + 1.0;
										theta = comptonSim(e_incident);
										e_scatt = e_incident
												/ (1 + (e_incident / 0.511)
														* (1 - Math.cos(theta)));
										endep[isValidInteraction - 1] = endep[isValidInteraction - 1]
												+ e_incident - e_scatt;

										//enmedtemp = endep[isValidInteraction - 1]
										//		/ contorInteraction[isValidInteraction - 1];
										enindivtemp = e_incident - e_scatt;
										// endeperror[isValidInteraction-1]=endeperror[isValidInteraction-1]+
										// Math.pow(enindivtemp-enmedtemp,2);
										endeperror[isValidInteraction - 1] = endeperror[isValidInteraction - 1]
												+ enindivtemp * enindivtemp;
										// separated ABM-----for skeleton
										// eval--------------------------------
										if (isValidInteraction == Phantom.SKELETON) {
											abmtempen = abmtempen
													+ abmrapinterp * abmw
													* enindivtemp;
											// abmtemperror=abmtemperror+
											// Math.pow(abmrapinterp*abmw*enindivtemp-abmtempen/contorInteraction[isValidInteraction-1],2);
											abmtemperror = abmtemperror
													+ abmrapinterp * abmw
													* enindivtemp
													* abmrapinterp * abmw
													* enindivtemp;
										}

										// continue tracing photon hystory
										e_incident = e_scatt;// new loop in
																// while
										// azimutal angle
										r0 = RandomCollection.random01();// [0-1]
										phi = 2 * r0 * Math.PI;// ->[0.2PI];
										// ----------------------------------------RAD---------------
										newcoord = getUpdatedCoordinate(theta,
												phi, dist_to_int, oldcrd);
										// -------------------!!!!!!!!!_____________________________
										isValidInteraction = phantom
												.inWhatOrgan(newcoord);// retrieve
																		// the
																		// organ
										// HERE->organisminclude inBody!
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
							} else {
								// pathlength -->exit the phantom
								n_photons = 0;// exit region of interest, so we
												// kill the photon
								n_detected++;
								if (projection == 0)// corect are loc initial
													// mulajul----legs=virtual=cu
													// trunk
								{
									if (oldcrd[1] > 0.0)// phantom.getOrganProjThickness(projection,oldcrd[2])/2)
									{
										if (!interactB)
											n_film_direct++;
										n_film++;
									}
								}
								if (projection == 2) {
									if (oldcrd[1] < 0.0)// -phantom.getOrganProjThickness(projection,oldcrd[2])/2)
									{
										if (!interactB)
											n_film_direct++;

										n_film++;
									}
								}
								if (projection == 1) {
									if (oldcrd[0] > 0.0)// phantom.getOrganProjThickness(projection,oldcrd[2])/2)
									{
										if (!interactB)
											n_film_direct++;

										n_film++;
									}
								}
								if (projection == 3) {
									if (oldcrd[0] < 0.0)// -phantom.getOrganProjThickness(projection,oldcrd[2])/2)
									{
										if (!interactB)
											n_film_direct++;

										n_film++;
									}
								}

								n_side = n_detected - n_film;
							}
						}// end while
					}// end for i=0 to photon_number
						// NOW WE HAVE DISTRIBUTION OF ENERGY IN ALL ORGANS FOR
						// e.g 10000 phgoton hystories.
						// error and we try the follow code.MAYBE DUE TO integer
						// multiplication !!!!!!!!!!!@@@@@@@@!!
					for (int k = 0; k < organs.length; k++) {
						/*
						 * if (contorInteraction[k]>0) { if
						 * (contorInteraction[k]>1) { //standard deviation of
						 * mean endeperror[k]=
						 * Math.sqrt((endeperror[k]/contorInteraction
						 * [k])/(contorInteraction[k]-1)); }//else is
						 * unchanged!! }//else unchanged
						 */
						// ---REAL @@ normalized per
						// energy-------------------------------------
						// endeperror[k]=endeperror[k]/photon_number;
						// scale spectrum
						// endeperror[k]=endeperror[k]*(xrs.getXRayIntensities())[j]/xrs.getNormalizedValue();
						// REAL @@ normalized per
						// energy-------------------------------------
						endep[k] = endep[k] / photon_number;
						// @@@@@@@@@@@@
						endeperror[k] = endeperror[k] / photon_number;
						endeperror[k] = (endeperror[k] - endep[k] * endep[k])
								/ (photon_number - 1.0);
						if (endeperror[k] >= 0.0)
							endeperror[k] = Math.sqrt(endeperror[k]);
						endeperror[k] = endeperror[k]
								* (xrs.getXRayIntensities())[j]
								/ xrs.getNormalizedValue();
						// @@@@@@@@@@@@@

						// ----------------------------------------------SCALE
						// SPECTRUM---------------
						endep[k] = endep[k] * (xrs.getXRayIntensities())[j]
								/ xrs.getNormalizedValue();
						// ----------------NOW REAL
						// VARIABLES------------------------------------------
						energy_deposit[k][0] = energy_deposit[k][0] + endep[k];
						energy_deposit[k][1] = Math.sqrt(energy_deposit[k][1]
								* energy_deposit[k][1] + endeperror[k]
								* endeperror[k]);
					}
					// separated ABM--from skeleton eval
					/*
					 * if (contorInteraction[phantom.SKELETON-1]>0) { if
					 * (contorInteraction[phantom.SKELETON-1]>1) { //standard
					 * deviation of mean abmtemperror=
					 * Math.sqrt((abmtemperror/contorInteraction
					 * [phantom.SKELETON
					 * -1])/(contorInteraction[phantom.SKELETON-1]-1)); }//else
					 * is unchanged!! }//else unchanged
					 */
					// ---REAL @@ normalized per
					// energy-------------------------------------
					// abmtemperror=abmtemperror/photon_number;
					// scale spectrum
					// abmtemperror=abmtemperror*(xrs.getXRayIntensities())[j]/xrs.getNormalizedValue();
					// ---REAL @@ normalized per
					// energy-------------------------------------
					abmtempen = abmtempen / photon_number;
					abmtemperror = abmtemperror / photon_number;
					abmtemperror = (abmtemperror - abmtempen * abmtempen)
							/ (photon_number - 1.0);
					if (abmtemperror >= 0.0)
						abmtemperror = Math.sqrt(abmtemperror);
					abmtemperror = abmtemperror * (xrs.getXRayIntensities())[j]
							/ xrs.getNormalizedValue();

					// ----------------------------------------------SCALE
					// SPECTRUM---------------
					abmtempen = abmtempen * (xrs.getXRayIntensities())[j]
							/ xrs.getNormalizedValue();
					// ----------------NOW REAL VARIABLES IN ABM NOT IN
					// SKELETON-DIFF OF 1!!!------------
					energy_deposit[Phantom.ACTVEBONEMARROW - 1][0] = energy_deposit[Phantom.ACTVEBONEMARROW - 1][0]
							+ abmtempen;
					energy_deposit[Phantom.ACTVEBONEMARROW - 1][1] = Math
							.sqrt(energy_deposit[Phantom.ACTVEBONEMARROW - 1][1]
									* energy_deposit[Phantom.ACTVEBONEMARROW - 1][1]
									+ abmtemperror * abmtemperror);
					// ---------------------
				}// end energy loop
				if (slicePerStep != 0)
					exposure = airExposure * 1000 * slicePerStep / 4;// -mGy->microGy
																		// and
																		// for
																		// single
																		// projection
																		// !!!
				else
					exposure = airExposure * 1000 * scannedslices / 4;// -mGy->microGy
																		// and
																		// for
																		// single
																		// projection
																		// !!!
				mevtojoule = 1.60218E-13;// 1MeV=1,60218 � 10-13 J
				factor = transversalArea * exposure * mevtojoule
						* xrs.getPhotonFlux() / xrs.getAirKerma();
				// separated ABM
				if (energy_deposit[Phantom.ACTVEBONEMARROW - 1][0] != 0.0)
					energy_deposit[Phantom.ACTVEBONEMARROW - 1][1] = Math.min(2
							* 100
							* energy_deposit[Phantom.ACTVEBONEMARROW - 1][1]
							/ energy_deposit[Phantom.ACTVEBONEMARROW - 1][0],
							99.9);// %
				else
					energy_deposit[Phantom.ACTVEBONEMARROW - 1][1] = 99.9;// 0.0;
					// normalized at real XRay Spectrum:
				energy_deposit[Phantom.ACTVEBONEMARROW - 1][0] = energy_deposit[Phantom.ACTVEBONEMARROW - 1][0]
						* factor;// Jouli
				// organ dose:
				// ////////////////////////////////
				organAbsorbedDose[Phantom.ACTVEBONEMARROW - 1][0] = 1000 * 1000
						* energy_deposit[Phantom.ACTVEBONEMARROW - 1][0]
						/ phantom.getOrganMass(Phantom.ACTVEBONEMARROW - 1);
				// 1g=10-03kk,Gy->mGy
				organAbsorbedDose[Phantom.ACTVEBONEMARROW - 1][1] = energy_deposit[Phantom.ACTVEBONEMARROW - 1][1];// %
				// -------------@@@@@@@@@@@
				orgAbsorbedDose[Phantom.ACTVEBONEMARROW - 1][0] = orgAbsorbedDose[Phantom.ACTVEBONEMARROW - 1][0]
						+ organAbsorbedDose[Phantom.ACTVEBONEMARROW - 1][0];
				// 1000*1000*energy_deposit[phantom.ACTVEBONEMARROW-1][0]/phantom.getOrganMass(phantom.ACTVEBONEMARROW-1);
				// 1g=10-03kk,Gy->mGy
				// orgAbsorbedDose[phantom.ACTVEBONEMARROW-1][1]=Math.sqrt(orgAbsorbedDose[phantom.ACTVEBONEMARROW-1][1]*orgAbsorbedDose[phantom.ACTVEBONEMARROW-1][1]+
				// organAbsorbedDose[phantom.ACTVEBONEMARROW-1][1]*organAbsorbedDose[phantom.ACTVEBONEMARROW-1][1]*organAbsorbedDose[phantom.ACTVEBONEMARROW-1][0]*organAbsorbedDose[phantom.ACTVEBONEMARROW-1][0]/10000);
				// energy_deposit[phantom.ACTVEBONEMARROW-1][1]*energy_deposit[phantom.ACTVEBONEMARROW-1][1]*energy_deposit[phantom.ACTVEBONEMARROW-1][0]*energy_deposit[phantom.ACTVEBONEMARROW-1][0]/10000);
				orgAbsorbedDose[Phantom.ACTVEBONEMARROW - 1][1] = Math.max(
						orgAbsorbedDose[Phantom.ACTVEBONEMARROW - 1][1],
						energy_deposit[Phantom.ACTVEBONEMARROW - 1][1]);
				// -------------------------------
				for (int k = 0; k < organs.length; k++) {
					if (k != Phantom.ACTVEBONEMARROW - 1) {
						// //% and 2 here for infinit grdlib (Nmare) is 95% conf
						// level
						if (energy_deposit[k][0] != 0.0)
							energy_deposit[k][1] = Math.min(2 * 100
									* energy_deposit[k][1]
									/ energy_deposit[k][0], 99.9);// %
						else
							energy_deposit[k][1] = 99.9;// 0.0;
						// normalized at real XRay Spectrum:
						energy_deposit[k][0] = energy_deposit[k][0] * factor;// Jouli
						// organ dose:
						// /////////////////////
						organAbsorbedDose[k][0] = 1000 * 1000
								* energy_deposit[k][0]
								/ phantom.getOrganMass(k);
						// 1g=10-03kk,Gy->mGy
						organAbsorbedDose[k][1] = energy_deposit[k][1];// %
						// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
						orgAbsorbedDose[k][0] = orgAbsorbedDose[k][0]
								+ organAbsorbedDose[k][0];
						// 1000*1000*energy_deposit[k][0]/phantom.getOrganMass(k);
						// 1g=10-03kk,Gy->mGy
						// orgAbsorbedDose[k][1]=Math.sqrt(orgAbsorbedDose[k][1]*orgAbsorbedDose[k][1]+
						// organAbsorbedDose[k][1]*organAbsorbedDose[k][1]*organAbsorbedDose[k][0]*organAbsorbedDose[k][0]/10000);
						// energy_deposit[k][1]*energy_deposit[k][1]*energy_deposit[k][0]*energy_deposit[k][0]/10000);
						// ----------------------------------------------------------------
						// aprox:
						orgAbsorbedDose[k][1] = Math.max(orgAbsorbedDose[k][1],
								energy_deposit[k][1]);
					}
				}
			}// end projection loop
			// orgAbsorbedDose[k][1]=orgAbsorbedDose[k][1]*orgAbsorbedDose[k][0]/100.0;//number
			if (slicePerStep != 0)
				currentscannedslice = currentscannedslice + slicePerStep;
			else
				currentscannedslice = currentscannedslice
						+ new Double(Math.ceil(scannedslices)).intValue();
			startS = endS;
		}// end while loop

		// 33333333333333333333333333333333333333333333333
		for (int k = 0; k < organs.length; k++) {
			if (contorInteraction_total[k] > 0.0)
				orgAbsorbedDose[k][1] = Math.min(
						2.0 * 100.0 * Math.sqrt(contorInteraction_total[k])
								/ contorInteraction_total[k], 99.9);
			else
				orgAbsorbedDose[k][1] = 99.9;
		}
		orgAbsorbedDose[Phantom.ACTVEBONEMARROW - 1][1] = orgAbsorbedDose[Phantom.SKELETON - 1][1];
		// 33333333333333333333333333333333333333333333333

		double[] wt = ((double[]) resources.getObject("wt.organname"));
		effectiveDose[0] = 0.0;
		double seff2 = 0.0;
		double errmin = 99.9;// initialisation
		for (int k = 0; k < organs.length; k++) {
			effectiveDose[0] = effectiveDose[0] + wt[k] * orgAbsorbedDose[k][0];
			// seff2=seff2+Math.pow(wt[k]*orgAbsorbedDose[k][1],2);
			seff2 = seff2
					+ Math.pow(wt[k] * orgAbsorbedDose[k][1]
							* orgAbsorbedDose[k][0] / 100.0, 2);
			errmin = Math.min(errmin, orgAbsorbedDose[k][1]);
		}
		if (effectiveDose[0] != 0.0)
			effectiveDose[1] = errmin;// Math.min((100*Math.sqrt(seff2))/effectiveDose[0],99.9);//%
		else
			effectiveDose[1] = 99.9;// 0.0;

		simulationTimeElapsed = timeElapsed(startSimulationTime);
	}
}
