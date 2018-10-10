package danfulea.phys;

import java.util.Calendar;
import java.util.Date;

import danfulea.math.numerical.LinAEq;

/**
 * Several physics utilities including BATEMAN decay law. <br>
 * 
 * @author Dan Fulea, 04 JUL 2006.
 * 
 */
public class PhysUtilities {
	public static double[][] nt;// number of nuclides at time t in a series at a
								// time step
	public static double[][] at;// corresponding activities;
	public static double t;// time for evaluation
	public static int steps;// step interval, e.g 10
	public static double t12m = 0.0;// maximum activity for 2nd nuclide at this
									// time

	public static final int IAIR = 0;
	public static final int ISOFTTISSUE = 1;
	public static final int IALL = 2;
	public static final int IWATTER = 3;

	public static double[][] ntt;// related to att!!!
	public static double[][] att;// corresponding time-integrated activities;

	public static boolean failB = false;

	/**
	 * Return the final nuclide activity based on its initial activity a0 computed at day 
	 * d0 (1-31), month m0 (1-12), year y0 and nuclide half-life in seconds t12_sec. 
	 * The final activity is computed at day d, month m and year y. 
	 * @param a0 a0
	 * @param t12_sec t12_sec
	 * @param d0 d0
	 * @param m0 m0
	 * @param y0 y0
	 * @param d d
	 * @param m m
	 * @param y y
	 * @return the result
	 */
	public static double decayLaw(double a0, double t12_sec, int d0, int m0,
			int y0, int d, int m, int y) {
		double a = a0;
		failB = false;
		if (t12_sec == 0.0)
			failB = true;

		double delta = 0;
		int mm0 = m0 - 1;
		int mm = m - 1;
		Calendar cal0 = Calendar.getInstance();
		cal0.set(y0, mm0, d0);
		Date dt0 = cal0.getTime();

		Calendar cal = Calendar.getInstance();
		cal.set(y, mm, d);
		Date dt = cal.getTime();
		long delt = dt.getTime();
		delt = delt - dt0.getTime();
		delta = delt;
		delta = delta / 1000.0;// milisec->sec
		// System.out.println("delta seconds= "+delta);
		if (!failB) {
			double b = -Math.log(2.0) * delta / t12_sec;
			a = a * Math.exp(b);
		}

		return a;
	}

	/**
	 * Prerequisite to Bateman. It computes the decay constant for a given (array of) 
	 * half-life. LAMBDA = Ln(2)/HalfLife
	 * @param t nuclide half life array
	 * @return the result
	 */
	public static double[] computeDecayConstant(double[] t) {
		double[] result = new double[t.length];
		for (int i = 1; i <= t.length; i++) {
			result[i - 1] = Math.log(2.0) / t[i - 1];
			// System.out.println("t1/2[s]= "+t[i-1]+" decayConstant= "+result[i-1]);
		}
		return result;
	}

	/**
	 * Prerequisite to Bateman. It computes the partial decay constant.
	 * @param t nuclide half life array
	 * @param b nuclide branching ratio (in chain) array
	 * @return the result
	 */
	public static double[] computePartialDecayConstant(double[] t, double[] b) {
		double[] result = new double[t.length];
		for (int i = 1; i <= t.length; i++) {
			result[i - 1] = b[i - 1] * Math.log(2.0) / t[i - 1];
			// System.out.println("t1/2[s]= "+t[i-1]+" pdecayConstant= "+result[i-1]);
		}
		return result;
	}

	/**
	 * Prerequisite to Bateman. It computes the initial number of atoms for each nuclides. 
	 * @param m the nuclide mass array
	 * @param a the nuclide atomic mass array
	 * @return the result
	 */
	public static double[] computeInitialNuclideNumbers(double[] m, double[] a) {
		double navogadro = 6.02214199E26;// atoms/kmol!!
		double[] result = new double[m.length];
		for (int i = 1; i <= m.length; i++) {
			result[i - 1] = m[i - 1] * navogadro / a[i - 1];
			// System.out.println("m= "+m[i-1]+" no= "+result[i-1]);
		}

		return result;
	}

	/**
	 * Used internally by Bateman.
	 * @param i i
	 * @param n n
	 * @param pdc pdc
	 * @return the result
	 */
	private static double computePartialDecaySum(int i, int n, double[] pdc) {
		double result = 1.0;
		for (int j = i; j <= n; j++) {
			result = result * pdc[j - 1];
		}
		return result;
	}

	/**
	 * Used internally by Bateman.
	 * @param i i
	 * @param j j
	 * @param n n
	 * @param time time
	 * @param q0 q0
	 * @param dc dc
	 * @return the result
	 */
	private static double compute_1st_term(int i, int j, int n, double time,
			double[] q0, double[] dc) {
		double result = q0[i - 1] * Math.exp(-dc[j - 1] * time);
		double l = 1.0;
		for (int p = i; p <= n; p++) {
			if (p != j)
				l = l * (dc[p - 1] - dc[j - 1]);
		}
		return result / l;
	}

	/**
	 * Used internally by Bateman.
	 * @param i i
	 * @param j j
	 * @param n n
	 * @param time time
	 * @param s0 s0
	 * @param dc dc
	 * @return the result
	 */
	private static double compute_2nd_term(int i, int j, int n, double time,
			double[] s0, double[] dc) {
		double result = s0[i - 1] * (1.0 - Math.exp(-dc[j - 1] * time));
		double l = 1.0;
		for (int p = i; p <= n; p++) {
			if (p != j)
				l = l * (dc[p - 1] - dc[j - 1]);
		}
		return result / (dc[j - 1] * l);
	}

	/**
	 * Used internally by Bateman.
	 * @param i i
	 * @param j j
	 * @param n n
	 * @param time time
	 * @param q0 q0
	 * @param dc dc
	 * @return the result
	 */
	private static double compute_1st_term_t(int i, int j, int n, double time,
			double[] q0, double[] dc) {
		double result = q0[i - 1] * (1.0 - Math.exp(-dc[j - 1] * time))
				/ dc[j - 1];
		double l = 1.0;
		for (int p = i; p <= n; p++) {
			if (p != j)
				l = l * (dc[p - 1] - dc[j - 1]);
		}
		return result / l;
	}

	/**
	 * Used internally by Bateman.
	 * @param i i
	 * @param j j
	 * @param n n
	 * @param time time
	 * @param s0 s0
	 * @param dc dc
	 * @return the result
	 */
	private static double compute_2nd_term_t(int i, int j, int n, double time,
			double[] s0, double[] dc) {
		// in usual cases s=0!!!!!!=>this term does not matter!
		double result = s0[i - 1]
				* (time - (1.0 - Math.exp(-dc[j - 1] * time)) / dc[j - 1]);
		double l = 1.0;
		for (int p = i; p <= n; p++) {
			if (p != j)
				l = l * (dc[p - 1] - dc[j - 1]);
		}
		return result / (dc[j - 1] * l);
	}

	/**
	 * Compute chain activities using Bateman decay law.
	 * @param n, numbers of nuclides in chain.
	 * @param dcst, decay constant array for those nuclides
	 * @param pdcst, partial decay constant array for those nuclides
	 * @param n0, initial number of atoms for those nuclides. Common scenario, 
	 * only parent nuclide is present at initial time.  
	 * @param s, additional sources term which contributes by injecting nuclides in chain. 
	 * Common scenario, no external sources (zero array) just natural decay.  
	 */
	public static void bateman(int n, double[] dcst, double[] pdcst,
			double[] n0, double[] s) {
		double deltat = t / steps;
		double time = 0.0;
		nt = new double[n][steps + 1];
		at = new double[n][steps + 1];

		ntt = new double[n][steps + 1];
		att = new double[n][steps + 1];

		double d = 0.0, d1 = 0.0, d2 = 0.0, d12 = 0.0, dd = 0.0;
		double d1t = 0.0, d2t = 0.0, d12t = 0.0, ddt = 0.0;
		// ==========================================
		if (n == 2 && n0[1] == 0.0 && s[0] == 0.0 && s[1] == 0.0
				&& dcst[1] != dcst[0]) {
			t12m = (Math.log(dcst[1] / dcst[0])) * 1.0 / (dcst[1] - dcst[0]);
		}
		// ==========================================
		for (int is = 0; is <= steps; is++) {
			time = is * deltat;// compute time
			for (int m = 1; m <= n; m++) {
				nt[m - 1][is] = 0.0;
				ntt[m - 1][is] = 0.0;
				for (int i = 1; i <= m; i++) {
					d = computePartialDecaySum(i, m - 1, pdcst);

					d12 = 0.0;
					d12t = 0.0;
					for (int j = i; j <= m; j++) {
						d1 = compute_1st_term(i, j, m, time, n0, dcst);
						d2 = compute_2nd_term(i, j, m, time, s, dcst);
						d12 = d12 + d1 + d2;

						d1t = compute_1st_term_t(i, j, m, time, n0, dcst);
						d2t = compute_2nd_term_t(i, j, m, time, s, dcst);
						d12t = d12t + d1t + d2t;
					}// j=i->m
					dd = d12 * d;
					ddt = d12t * d;
					nt[m - 1][is] = nt[m - 1][is] + dd;
					ntt[m - 1][is] = ntt[m - 1][is] + ddt;
				}// i=1->m
				at[m - 1][is] = nt[m - 1][is] * dcst[m - 1];
				att[m - 1][is] = ntt[m - 1][is] * dcst[m - 1];
			}// m=1->n, all nuclides
		}// steps, all steps
	}

	// energy->MeV, Yield<1, activity->Bq, distance->cm, coeff->cm2/g
	// no attenuation in medium->AIR!!!
	/**
	 * Given photon energy in MeV, its yield (probability of emission) as a number less or equal to 1, the source activity in Bq 
	 * which emits those photons, the mass-energy absorbtion coefficient in cm^2/g taken 
	 * from literature (or by other means) and the distance from the radioactive source in cm, this routine returns 
	 * the dose debit due to absorbtion of gamma radiation at that distance in microGy/hour. 
	 * @param photonEnergy photonEnergy
	 * @param yield yield
	 * @param sourceActivity sourceActivity
	 * @param calculationPointDistance calculationPointDistance
	 * @param massEnergyAbsCoeff massEnergyAbsCoeff
	 * @return the result
	 */
	public static double gammaAbsorbedDoseDebit(double photonEnergy,
			double yield, double sourceActivity,
			double calculationPointDistance, double massEnergyAbsCoeff) {
		// for a given energy a massEnergyAbsCoeff interpolation coeff can be
		// computed by interpolation
		// fluenta=Y*Activ/(4 PI x^2), in fotoni/(m^2 sec) .Nr de fotoni la
		// distanta x de sursa per unitatea de timp si de suprafata

		double mevToJ = 1.60218E-13;
		double result = yield
				* sourceActivity
				* photonEnergy
				* massEnergyAbsCoeff
				/ (4.0 * Math.PI * calculationPointDistance * calculationPointDistance);
		// in mev/g/sec
		// now:
		result = result * mevToJ * 1000 * 1000;// in mGy/sec
		result = result * 1000 * 3600;// microGy/h
		return result;
	}

	// exposure time in sec.
	/**
	 * Given photon energy in MeV, its yield (probability of emission) as a number less or equal to 1, the source activity in Bq 
	 * which emits those photons, the mass-energy absorbtion coefficient in cm^2/g taken 
	 * from literature (or by other means), the distance from the radioactive source in cm and the exposure time in seconds at that location, this routine returns 
	 * the dose due to absorbtion of gamma radiation at that distance in microGy.  
	 * @param photonEnergy photonEnergy
	 * @param yield yield
	 * @param sourceActivity sourceActivity
	 * @param calculationPointDistance calculationPointDistance
	 * @param massEnergyAbsCoeff massEnergyAbsCoeff
	 * @param exposureTime exposureTime
	 * @return the result
	 */
	public static double gammaAbsorbedDose(double photonEnergy, double yield,
			double sourceActivity, double calculationPointDistance,
			double massEnergyAbsCoeff, double exposureTime) {
		double result = gammaAbsorbedDoseDebit(photonEnergy, yield,
				sourceActivity, calculationPointDistance, massEnergyAbsCoeff);
		result = result * exposureTime / 3600.0;// in microGy
		return result;
	}

	// compute maximum massic path length of beta electrons; semi empirical
	// formula
	// maximum beta Energy->in MeV; output in g/cm2
	/**
	 * Compute the maximum mass-path length of electrons (beta radiation) in g/cm^2 for given 
	 * energy (MeV). Semi-empirical formulae is used.
	 * @param maxBetaEnergy maxBetaEnergy
	 * @return the result
	 */
	public static double betaMaxMassPathLength(double maxBetaEnergy) {
		double result = 1.0;
		if (maxBetaEnergy <= 0.15) {
			result = 0.00667 * Math.pow(maxBetaEnergy, 5.0 / 3.0);
		} else if (maxBetaEnergy <= 0.8) {
			result = 0.407 * Math.pow(maxBetaEnergy, 1.38);
		} else {
			result = 0.542 * maxBetaEnergy - 0.133;
		}

		return result;
	}

	// density in g/cm3; output in cm
	/**
	 * Compute the maximum path length of electrons (beta radiation) in cm for given 
	 * energy (MeV) and medium density (g/cm^3). Semi-empirical formulae is used. 
	 * @param maxBetaEnergy maxBetaEnergy
	 * @param mediumDensity mediumDensity
	 * @return the result
	 */
	public static double betaMaxPathLength(double maxBetaEnergy,
			double mediumDensity) {
		double result = 1.0;
		result = (betaMaxMassPathLength(maxBetaEnergy)) / mediumDensity;
		return result;
	}

	// output in cm2/g
	/**
	 * Computes the beta mass attenuation coefficient (cm^2/g) for a given energy (MeV) and several 
	 * media (AIR, SOFTTISSUE or WATTER, Aluminum). Semi-empirical formulae is used. 
	 * @param maxBetaEnergy maxBetaEnergy
	 * @param mediumIndex mediumIndex
	 * @return the result
	 */
	public static double betaMassAttenuationCoeff(double maxBetaEnergy,
			int mediumIndex) {
		double result = 1.0;
		if (mediumIndex == IAIR) {
			result = 16.0 / Math.pow(maxBetaEnergy - 0.036, 1.37);
		} else if (mediumIndex == ISOFTTISSUE || mediumIndex == IWATTER) {
			result = 18.2 / Math.pow(maxBetaEnergy - 0.036, 1.37);
		} else if (mediumIndex == IALL) {
			result = 22.0 / Math.pow(maxBetaEnergy, 1.33);
		}

		return result;
	}

	// output in cm-1
	/**
	 * Computes the beta linear attenuation coefficient (cm^-1) for a given energy (MeV) and several 
	 * media (AIR, SOFTTISSUE or WATTER, Aluminum) with medium density given in g/cm^3. Semi-empirical formulae is used.
	 * @param maxBetaEnergy maxBetaEnergy
	 * @param mediumDensity mediumDensity
	 * @param mediumIndex mediumIndex
	 * @return the result
	 */
	public static double betaLinAttenuationCoeff(double maxBetaEnergy,
			double mediumDensity, int mediumIndex) {
		double result = 1.0;
		result = (betaMassAttenuationCoeff(maxBetaEnergy, mediumIndex))
				* mediumDensity;
		return result;
	}

	// Nuclides inside a sphere having radius=MaxPathLength,
	// generates beta rad. reaching at every point. Outside this range the beta
	// electrons are stopped!!
	// So, dose will be given by solving the integral
	// on volume; differential dose is given by::
	// dD=[ActivVol/(4*PI*r*r)]*dV*t*EffBetaEnergy*[1/MaxMassPathLength]*EXP(-LinAttMedium*r)
	// where dV=r*r*sin(theta)*dr*dtheta*dphi; r->0,MaxPathLength;
	// theta->0,PI/2; phi->0,2*PI
	// so if r<=MaxPathLength we have:
	// yield<1;energy->MeV;volActiv->Bq/cm3, time->sec
	// Note:Absorbed dose=Equivalent dose for photons,miuons and electrons
	// (factor=1)
	//Target is surrounded by radioactive cloud. All rad from a sphere of maximum radius 
	//equal to maxPathLength reach the target and contribute to dose!!
	
	//Valid scenarion is a softTissue surrounded by air or water. LinAttCoeff acts as a weight 
	//for radiation reaching the target hence: ActivVol/(4*PI*r*r)*dV*t*EXP(-LinAttMedium*r) 
	//is the weighted fluence, i.e. photons reaching the target due to absorbtion in surrounding 
	//medium. This must be multiplied with a mass-energy absorbtion LIKE coefficient for 
	//actual target (soft tissue) to compute dose. Since we lack this data, this method is deprecated. 
	@Deprecated
	public static double betaRadMaxDose(double yield, double betaMaxEnergy,
			double volActivity, double exposureTime, int mediumIndex,
			double mediumDensity) {
		double result = 0.0;
		double mevToJ = 1.60218E-13;
		double effectiveBetaEnergy = 0.4 * betaMaxEnergy;
		double linAttCoeff = betaLinAttenuationCoeff(betaMaxEnergy,
				mediumDensity, mediumIndex);// cm-1
		double maxPathLength = betaMaxPathLength(betaMaxEnergy, mediumDensity);// cm
		double maxMassPathLength = betaMaxMassPathLength(betaMaxEnergy);// g/cm2

		result = yield * volActivity * exposureTime * effectiveBetaEnergy
				* (1.0 - Math.exp(-linAttCoeff * maxPathLength))
				/ (linAttCoeff * maxMassPathLength);
		result = result * 1000. * 1000. * mevToJ;// mGy
		return result;
	}

	/**
	 * Compute Z equivalent for a mixture of elements. Used internally.
	 * @param z z array, z number for an element
	 * @param pz pz array, atomic fraction of an element in a compound.
	 * @return the result
	 */
	public static double zEquivalent(double[] z, double[] pz) {
		double result = 0.0;
		double sumpzz = 0.0;
		double sumpz = 0.0;
		for (int i = 1; i <= z.length; i++) {
			sumpz = sumpz + pz[i - 1];
			sumpzz = sumpzz + z[i - 1] * pz[i - 1];
		}
		result = sumpzz / sumpz;

		return result;
	}

	//Although good if and only if gammaMass reffers o tissue and gammaLin reffers to 
	//water as surrounding media, it is based on betaMaxEnergy and must be used in 
	//conjuction with betaRadMaxDose method which is deprecated due to lack of data.
	@Deprecated
	public static double betaBremsRadDose(double yield, double betaMaxEnergy,
			int mediumIndex, double calcPointDistance,
			double gammaMassEnergyAbsCoeff, double gammaLinAttCoeff,
			double volActivity, double exposureTime) {
		double result = 0.0;
		double zech = 1.0;
		if (mediumIndex == IAIR) {
			// -------
		} else if (mediumIndex == ISOFTTISSUE) {
			// ----------
		} else if (mediumIndex == IALL) {
			zech = 13.;
		} else if (mediumIndex == IWATTER) {
			double[] z = { 1., 8. };
			double[] pz = { 0.111009256, 0.0555087065764966 };
			zech = zEquivalent(z, pz);
		}
		double f = zech * betaMaxEnergy / 800.;// fraction of max beta energy
												// converted to brems!!
		// having of course betaMaxEnergy (calc point far away outside the
		// Rmax!! )!!

		double mevToJ = 1.60218E-13;

		result = yield * volActivity * exposureTime * f * betaMaxEnergy
				* gammaMassEnergyAbsCoeff
				* (1.0 - Math.exp(-gammaLinAttCoeff * calcPointDistance))
				/ (gammaLinAttCoeff);
		result = result * 1000. * 1000. * mevToJ;// mGy

		return result;
	}

	// point dose surrounding by radioactive medium!!
	// gammaMassEnergyAbsCoeff->in TISSUE human dose; cm2/g
	// gammaLinAttCoeff in medium, e.g watter,/cm-1
	/**
	 * Compute dose at a target point P surrounded by a radioactive medium. Common ussage is 
	 * to consider a soft tissue medium (simulates the human body) at P and surrounded by 
	 * air or water (air or water submersion). The radiation fluence proportional with 
	 * the number of photons reaching the target coming from volume dV at r distance from target is: 
	 * dPHI = yield*[volActivity/(4*PI*r*r)]*dV*exposureTime*EXP(-LinAttMedium*r). Dose at target is 
	 * dD = dPHI * gammaEnergy * gammaMassEnergyAbsCoeff. The infinitesimal volume 
	 * dV=r*r*sin(theta)*dr*dtheta*dphi; r in range 0,MaxValue; theta in range 0,PI; phi in range 0,2*PI 
	 * Integration is performed analytically here considering all radiations coming from within a sphere 
	 * of radius MaxValue reach the target.
	 * @param yield yield
	 * @param gammaEnergy gammaEnergy in MeV
	 * @param calcPointDistance is MaxValue considered in integration in cm
	 * @param gammaMassEnergyAbsCoeff gammaMassEnergyAbsCoeff in cm^2/g for target (e.g. soft tissue)
	 * @param gammaLinAttCoeff gammaLinAttCoeff for surrounding media in cm^-1 (e.g. air)
	 * @param volActivity volActivity in Bq/cm^3
	 * @param exposureTime exposureTime in seconds
	 * @return the absorbed dose in target in mGy
	 */
	public static double gammaRadDose(double yield, double gammaEnergy,
			double calcPointDistance, double gammaMassEnergyAbsCoeff,
			double gammaLinAttCoeff, double volActivity, double exposureTime) {
		double result = 0.0;

		double mevToJ = 1.60218E-13;

		result = yield * volActivity * exposureTime * gammaEnergy
				* gammaMassEnergyAbsCoeff
				* (1.0 - Math.exp(-gammaLinAttCoeff * calcPointDistance))
				/ (gammaLinAttCoeff);
		result = result * 1000. * 1000. * mevToJ;// mGy

		return result;
	}

	// U in kV, filtration in mmAL, absorbanmt distance is usually 8 cm from
	// focus!!
	// for normalised results time and ma must be 1
	/**
	 * Given an XRay tube voltage in kV, current in mA, Z for anode, exposure time in seconds, 
	 * mass-energy absorbtion coefficient in air (cm^2/g) at average photon energy in spectrum, 
	 * the tube total filtration in mmAl and the distance in cm, this routine returns the 
	 * kerma in air at that distance in microGy. Average photon enery ~ (2/3)voltage/1000. 
	 *   
	 * @param voltage voltage
	 * @param zAnode zAnode
	 * @param mA mA
	 * @param exposureTime exposureTime
	 * @param gammaMassAbsAirCoeff gammaMassAbsAirCoeff
	 * @param filtration filtration
	 * @param absorbantDist absorbantDist
	 * @return the result
	 */
	public static double xTubeAirKermaAtExit(double voltage, double zAnode,
			double mA, double exposureTime, double gammaMassAbsAirCoeff,
			double filtration, double absorbantDist) {
		double result = 0.0;
		double mevToJ = 1.60218E-13;
		double electronMaxEnergy = voltage / 1000.0;// MeV
		double f = zAnode * electronMaxEnergy / 800.;// fraction transferred to
														// brems
		double eCharge = 1.60217646E-19;// Coulomb
		double eqActiv = f * mA * 0.001 / eCharge;// f*I/e; mA->A!!
		// distrib function of X spectra:f= C*k^2(kmax-k)
		// medium energy from distrib function, is 2eU/3; Photon max
		// energy=electronMaxEnergy
		double averagePhotonEnergy = 2.0 * electronMaxEnergy / 3.0;
		// previously one must compute gammaMassAbsAirCoeff at
		// averagePhotonEnergy in air!!
		// dose is A*t/(4*PI*(filtr+x)^2)*kmean*(miu/rho)absAIR
		result = 1000.
				* 1000.
				* mevToJ
				* eqActiv
				* exposureTime
				* averagePhotonEnergy
				* gammaMassAbsAirCoeff
				/ (4. * Math.PI * (absorbantDist + filtration / 10.0) * (absorbantDist + filtration / 10.0));
		result = result * 1000.0;// ->microGy /mAs!!
		return result;
	}

	/**
	 * Test routine
	 */
	public static void test() {
		// number of nuclides
		int nuclideNumber = 2;
		// partial decay constant for n1->n2->...
		// ki=BRi x k; k=total decay constant, BRi is branching ratio for i
		// decay
		// e.g Sr90->100%->Y90->~100%->Zr90(stable)
		// br and hf MUST HAVE THE SAME LENGTH!!
		double[] br = { 1.0, 1.0 };
		// total decay constant for n1->etc->n2->etc->...taken from halflifetime
		// double[] hf={29.1*365.25*24*3600,2.67*3600};
		double[] hf = { 1600 * 365.25 * 24 * 3600, 3.8235 * 3600 };
		// -------------------------------------------------------
		double[] dc = PhysUtilities.computeDecayConstant(hf);
		double[] pdc = PhysUtilities.computePartialDecayConstant(hf, br);
		// ------------------------------------------------------
		// the number of initial nuclides:
		// 1 grams of n1 and 0 grams of n2
		double[] mass = { 1.0E-03, 0.0 };
		// double[] atomicMass= {89.907,90.0};
		double[] atomicMass = { 226.0, 222.0 };
		// -----------------------------------------
		double[] initialNuclideNumbers = PhysUtilities
				.computeInitialNuclideNumbers(mass, atomicMass);
		// ----------------------------------------
		// the number of nuclides added by additional sources (flux)!!
		double[] sources = { 0.0, 0.0 };
		PhysUtilities.steps = 100;
		PhysUtilities.t = 1 * 365.25 * 24 * 3600 / 12.0;// 5 years

		PhysUtilities.bateman(nuclideNumber, dc, pdc, initialNuclideNumbers,
				sources);

		System.out.println("theor tmax eq ideal " + PhysUtilities.t12m);
		for (int j = 0; j <= steps; j++) {
			for (int i = 1; i <= nuclideNumber; i++) {
				// System.out.println("nuc: "+i+" n(t)= "+nt[i-1][j]+" at step: "+j+" time[s]= "+j*t/steps);

				System.out.println("nuc: " + i + " a(t)= "
						+ PhysUtilities.at[i - 1][j] + " at step: " + j
						+ " time[s]= " + j * PhysUtilities.t
						/ PhysUtilities.steps + "; ore= " + j * PhysUtilities.t
						/ (3600. * PhysUtilities.steps) + "; zile= " + j
						* PhysUtilities.t / (24 * 3600. * PhysUtilities.steps));
			}

			if (Math.abs(at[0][j] - at[1][j]) / at[0][j] < 0.000001)
				break;
		}

		double en = 0.6617;// MeV;
		double y = 0.8521;
		double activ = 30778.57772;// Bq
		double x = 2.0;// cm
		double miuPerRho = 0.03229937;// cm2/g
		double time = 21600;// sec

		double debit = PhysUtilities.gammaAbsorbedDoseDebit(en, y, activ, x,
				miuPerRho);
		double dose = PhysUtilities.gammaAbsorbedDose(en, y, activ, x,
				miuPerRho, time);

		System.out.println("dose= " + dose + " microGy" + "; debit= " + debit
				+ " microGy/h ");

		double ben = 1.325;// mev
		double gmassAbs = 0.029;// g/cm2
		double glinAtt = 0.06215;// cm-1
		double y2 = 0.89;
		double dens = 1;// g/cm3, water
		System.out.println("maxpathlength= "
				+ PhysUtilities.betaMaxPathLength(ben, dens) + " cm ");
		// betaLinAttenuationCoeff(double maxBetaEnergy,double mediumDensity,int
		// mediumIndex)
		System.out.println("betaLinAttCoeff= "
				+ PhysUtilities.betaLinAttenuationCoeff(ben, dens,
						PhysUtilities.IWATTER) + " cm-1 ");
		double volac = 1.e6;// Bq/cm3
		double tm = 300;// sec
		double d1 = PhysUtilities.betaRadMaxDose(y2, ben, volac, tm,
				PhysUtilities.IWATTER, dens);
		System.out.println(" dose hemisphere= " + d1 / 2.0 + " mGy");
		x = 300.0;// cm
		double d2 = PhysUtilities.betaBremsRadDose(y2, ben,
				PhysUtilities.IWATTER, x, gmassAbs, glinAtt, volac, tm);
		System.out.println(" hemis brems dose at x=300 cm= " + d2 / 2.0
				+ " mGy");

		double gen = 1.46;// mev
		double y3 = 0.11;
		double gmassAbs1 = 0.02828;// g/cm2
		double glinAtt1 = 0.05859;// cm-1
		double d3 = PhysUtilities.gammaRadDose(y3, gen, x, gmassAbs1, glinAtt1,
				volac, tm);
		System.out.println(" hemis gamma dose at x=300 cm= " + d3 / 2.0
				+ " mGy");

		System.out.println(" total from 40K= " + (d1 + d2 + d3) / 2.0 + " mGy");

		double kv = 80;// kv
		double zanod = 74;// wolfram
		double itube = 1;// mA
		double exptime = 1;// sec
		double mmal = 2.5;// filtrare mmAl;
		double absdist = 8.0;// cm
		double gamaairabs = 4.10E-2;// cm2/g
		d2 = xTubeAirKermaAtExit(kv, zanod, itube, exptime, gamaairabs, mmal,
				absdist);
		System.out.println(" X tube airkerma= " + d2 + " microGy/mAs");
	}

	// =================
	private static boolean valid = false;

	/**
	 * Check the success of some methods.
	 * @return true or false
	 */
	public static boolean validate() {
		return valid;
	}
	
	/**
	 * Compute dead time based on two sources having known activities quotient!
	 * 
	 * @param r12 the fraction A1/A2 (or real counts rate quotient R1/R2) 
	 * @param Ro1 the observed (measured) counts rate for source 1 (counts per second)
	 * @param Ro2 the observed (measured) counts rate for source 2 (counts per second)
	 * @return the detection dead time in seconds
	 */
	public static double computeDeadTime(double r12, double Ro1, double Ro2){
		/*
		 * Real counts rate versus observed (measured) counts rate is:
		 * R1=Ro1/(1-tRo1);R2=Ro2/(1-tRo2); Let r12=R1/R2 (=A1/A2) then:
		 * r12=(Ro1/Ro2)[(1-tRo2)/(1-tRo1)]
		 * r12(Ro2/Ro1)(1-tRo1)= (1-tRo2)
		 * r12(Ro2/Ro1)-r12Ro2t=1-Ro2t
		 * r12(Ro2/Ro1)-1=Ro2t(r12-1) 
		 * t=[r12(Ro2/Ro1)-1]/[Ro2(r12-1)] or multiplying all with Ro1:
		 * t=[r12Ro2-Ro1]/[Ro1Ro2(r12-1)]
		 */
		valid = false;
		double result=0.0;
		double num=Ro1*Ro2*(r12-1.0);
		if (num!=0.0){
			result=(r12*Ro2-Ro1)/num;
			valid = true;
		}
		if (result<=0.0)
			valid=false;
		return result;
	}

	/**
	 * Compute dead time based on two sources having unknown activities quotient!
	 * Requires a sensitive measurements geometry: 
	 * First measure source 1 and record its observed counts rate r1.
	 * Then insert a second source and record counts rate r12.
	 * Finally remove source 1 and record counts rate r2.
	 * Note: the measurement geometry must be the same!!
	 * Alternatives : Construct 3 sources in same geometry !!! First two
	 * sources as above and the third being equivalent with the sum of first
	 * two sources!!  
	 * @param r12 the observed (measured) counts rate for source1+source2 (counts per second) 
	 * @param r1 the observed (measured) counts rate for source 1 (counts per second) 
	 * @param r2 the observed (measured) counts rate for source 2 (counts per second) 
	 * @return the detection dead time in seconds
	 */
	public static double computeDeadTime2(double r12, double r1, double r2){
		/*
		 * Real counts rate versus observed (measured) counts rate is:
		 * R1=r1/(1-tr1);R2=r2/(1-tr2); R12=r12/(1-tr12).
		 * R12=R1+R2 (by construction) and neglecting terms in t*t we have:
		 * r12(1-r1t)(1-r2t)=r1(1-r12t)(1-r2t)+r2(1-r12t)(1-r1t)
		 * t=~(r1+r2-r12)/2r1r2
		 */
		valid = false;
		double result=0.0;
		double num=2.0*r1*r2;
		if (num!=0.0){
			result=(r1+r2-r12)/num;
			valid = true;
		}
		if (result<=0.0)
			valid=false;
		return result;
	}
		
	// calculeaza timpul mort al aparatelor de masurare a radioactiv alfa sau
	// beta
	// intra q1=viteza de numarare a sursei 1;
	// q2= viteza de numarare a sursei 2;
	// q12=viteza de numarare a celor 2 surse masurate simultan;
	// f=viteza de numarare a fondului;
	@Deprecated
	public static double tmort(double q1, double q2, double q12, double f) {
		valid = false;
		double tmort = 0.0;
		double num = q12 * q12 - q1 * q1 - q2 * q2;
		if (num != 0) {
			tmort = (q1 + q2 - q12 - f) / num;
			valid = true;
		}
		return tmort;
	}

	// calculeaza corectia de timp mort
	// x->rata (imp/timp)
	/**
	 * Compute the counts rate corrected by deadTime in counts per seconds.
	 * @param x the displayed count rate (counts per seconds).
	 * @param tmort the dead time in seconds.
	 * @return the result
	 */
	public static double tmortCorOfRate(double x, double tmort) {
		valid = false;
		double tmCor = 0.0;
		double num = 1 - x * tmort;
		if (num != 0) {
			tmCor = x / num;
			valid = true;
		}
		return tmCor;
	}

	// calculeaza corectia de timp mort
	// x->impulsurile masurate in timpul t;
	/**
	 * Compute the counts corrected by dead time.  
	 * @param x the counts
	 * @param t time taken to measure the counts in seconds
	 * @param tmort dead time in seconds
	 * @return the result
	 */
	public static double tmortCorOfImp(double x, double t, double tmort) {
		valid = false;
		double tmCor = 0.0;
		if (t != 0.0) {
			double num = 1 - x * tmort / t;
			if (num != 0) {
				tmCor = x / num;
				valid = true;
			}
		}
		return tmCor;
	}

	// calculeaza rata minima semnificativa (decelabila de fond) pentru setul de
	// masuratori alfa-beta
	// bazat pe F=rata medie a fondului utilizat, Nf=numarul de masuratori de
	// fond, fiecare in
	// timpul tf, Nq= numarul de masuratori a probei in prezenta fondului
	// ,fiecare in
	// timpul tq
	/**
	 * Compute the minimum detectable count rate in counts per second. 
	 * @param F average count rate for background in counts per second 
	 * @param Nf number of successive measurements for background
	 * @param tf time for each background measurement in seconds
	 * @param Nq number of successive measurements for sample
	 * @param tq time for each sample measurement in seconds
	 * @return the result
	 */
	public static double minRateOfDec(double F, int Nf, double tf, int Nq,
			double tq) {
		valid = false;
		/*
		 * double mrd=3.29; if(Nq*tq!=0.0 && Nf*tf!=0.0 && F>0) {
		 * mrd=mrd*(1.645/
		 * (Nq*tq)+Math.sqrt(2.706025/(Nq*tq*Nq*tq)+F*(1/(Nq*tq)+1/(Nf*tf))));
		 * valid=true; }
		 */

		double mrd = 0.0;

		if (tq > 0.0 && F > 0)
			mrd = (2.706025 / tq) + 3.29 * Math.sqrt(F / tq);

		return mrd;
	}

	/**
	 * The formal uncertainty associated with the minimum detectable count rate in counts per second.
	 * @param F average count rate for background in counts per second
	 * @param Nf number of successive measurements for background
	 * @param tf time for each background measurement in seconds
	 * @param Nq number of successive measurements for sample
	 * @param tq time for each sample measurement in seconds
	 * @return the result
	 */
	public static double minRateOfDecError(double F, int Nf, double tf, int Nq,
			double tq) {
		valid = false;
		double mrd = 1.645;
		/*
		 * if(Nq*tq!=0.0 && Nf*tf!=0.0 && F>0) {
		 * mrd=mrd*(1/(Nq*tq)+1/(Nf*tf))*Math.sqrt(F/(Nf*tf));
		 * mrd=mrd/(Math.sqrt(2.706025/(Nq*tq*Nq*tq)+F*(1/(Nq*tq)+1/(Nf*tf))));
		 * valid=true; }
		 */
		if (tq > 0.0 && tf > 0.0)
			mrd = mrd / (Math.sqrt(tf * tq));
		if (Nf > 1)
			mrd = mrd / (Math.sqrt(Nf));

		return mrd;
	}

	// legea dezintegrarilor radioactive
	// activitatea initiala a0 la data data de zi,luna,an
	// Obs an>=1900 si <=9999!!!
	// timpul de injumatatire e dat de t12 in ani
	/*
	 * public static double decayLaw(double a0, double t12, int zi, int luna,
	 * int an) { double a = a0; double t12zile = t12 * 365.25; valid = false;
	 * int deltazile = 0; Calendar cal = Calendar.getInstance(); int d =
	 * cal.get(Calendar.DATE); int m = cal.get(Calendar.MONTH) + 1; int y =
	 * cal.get(Calendar.YEAR);
	 * 
	 * try { SerialDate start = SerialDate.createInstance(zi, luna, an);
	 * SerialDate end = SerialDate.createInstance(d, m, y); deltazile =
	 * SerialDateUtilities.dayCountActual(start, end);
	 * 
	 * if (t12 != 0.0) valid = true; } catch (Exception e) { valid = false; }
	 * 
	 * if (valid) { double b = -Math.log(2.0) * deltazile / t12zile; a = a *
	 * Math.exp(b); }
	 * 
	 * return a; }
	 */

	// intra
	// n1,n2 ordinele polinoamelor de fitare->eff=f(ln(E))
	// ec=energia de crossover
	// x=valorile energiei,y=valorile eff masurate in % (exASSAYER sau
	// GAMMA2000)
	// tl=tabloul termenilor liberi care va stoca rezultatul final!!!
	// t0-eff la ec, t1,..,tn1+1 coef prim polinom, tn1+1,,,,tn1+n2+2 coef al
	// 2-lea polinom
	
	@Deprecated
	public static void polyGammaEff(double[] x, double[] y, int n1, int n2,
			double ec, double[][] tl) {
		/*
		*/
		// transformare energie-ln(energie) pe fiecare subinterval dat de ec
		int k = 0;
		for (int i = 0; i < x.length; i++)
			if (x[i] < ec)
				k++;
			else
				break;
		double[] xe1 = new double[k];
		double[] xe2 = new double[x.length - k];
		double[] ye1 = new double[k];
		double[] ye2 = new double[x.length - k];

		for (int i = 0; i < k; i++) {
			if (x[i] == 0.0)
				xe1[i] = 0.0;
			else
				xe1[i] = Math.log(x[i]);
			if (y[i] == 0.0)
				ye1[i] = 0.0;
			else
				ye1[i] = Math.log(y[i] / 100);
			// System.out.println(xe1[i]);
		}
		// System.out.println("-----------------------------");
		for (int i = k; i < x.length; i++) {
			if (x[i] == 0.0)
				xe2[i - k] = 0.0;
			else
				xe2[i - k] = Math.log(x[i]);
			if (y[i] == 0.0)
				ye2[i - k] = 0.0;
			else
				ye2[i - k] = Math.log(y[i] / 100);
			// System.out.println(xe2[i-k]);
		}
		// System.out.println("-----------------------------");
		double xc = 0.0;
		if (ec > 0)
			xc = Math.log(ec);
		// System.out.println(xc);
		// end
		if (xe1.length != 0 && xe2.length != 0) {
			// construiesc sirul termenilor liberi ai sistemului de ecuatii
			// obtinut prin metoda celor mai mici patrate!
			double tls, cf = 0.0;
			for (int j = 0; j <= n1; j++) {
				tls = 0.0;
				for (int i = 0; i < xe1.length; i++)
					tls = tls + ye1[i] * Math.pow(xe1[i], j);
				tl[j][0] = tls;
			}
			for (int j = 0; j <= n2; j++) {
				tls = 0.0;
				for (int i = 0; i < xe2.length; i++)
					tls = tls + ye2[i] * Math.pow(xe2[i], j);// y[k+i]<->ye2[i]!!
				tl[n1 + 1 + j][0] = tls;
			}
			tl[n1 + n2 + 3 - 1][0] = 0.0;
			// System.out.println("-----------tl----------------");
			// for (int j=0; j<tl.length; j++)
			// System.out.println(tl[j][0]);
			// System.out.println("-----------coef----------------");
			// construiesc sirul elementelor matricii sistemului de ecuatii
			double[][] coef = new double[n1 + n2 + 3][n1 + n2 + 3];

			// prima jumatate
			for (int i = 0; i <= n1; i++) {
				// y300 termen
				coef[i][0] = -Math.pow(xc, i);
				for (int j = 0; j <= n1; j++) {
					cf = 0.0;
					for (int kk = 0; kk < xe1.length; kk++)
						cf = cf + Math.pow(xe1[kk], i + j);
					coef[i][j + 1] = cf + Math.pow(xc, i + j);
				}

				for (int j = 0; j <= n2; j++) {
					coef[i][n1 + 1 + j + 1] = 0.0;
				}
			}
			// a doua jumatate

			for (int i = 0; i <= n2; i++) {
				// y300 termen
				coef[n1 + 1 + i][0] = -Math.pow(xc, i);
				for (int j = 0; j <= n1; j++) {
					coef[n1 + 1 + i][j + 1] = 0.0;
				}

				for (int j = 0; j <= n2; j++) {
					cf = 0.0;
					for (int kk = 0; kk < xe2.length; kk++)
						cf = cf + Math.pow(xe2[kk], i + j);
					coef[n1 + i + 1][n1 + 1 + j + 1] = cf + Math.pow(xc, i + j);

				}
			}

			// ultima linie
			coef[n1 + n2 + 2][0] = 0.0;
			for (int j = 0; j <= n1; j++) {
				cf = Math.pow(xc, j);
				coef[n1 + n2 + 2][j + 1] = cf;
			}
			for (int j = 0; j <= n2; j++) {
				cf = -Math.pow(xc, j);
				coef[n1 + n2 + 2][j + n1 + 1 + 1] = cf;
			}

			// for (int i=0; i<coef.length; i++)
			// {
			// System.out.println("");
			// for (int j=0; j<n1+n2+3; j++)
			// System.out.print(coef[i][j]+" ");
			// }
			// System.out.println("");
			// apelez rezolvarea sistemului de ecuatii
			LinAEq.sysEqGauss(coef, tl, n1 + n2 + 3, 1);
		}

		if (xe1.length == 0 || xe2.length == 0) {
			int n = 0;
			if (xe1.length == 0)
				n = n2;
			else
				n = n1;

			double[] yy = new double[x.length];
			double[] xx = new double[x.length];
			for (int i = 0; i < x.length; i++) {
				if (x[i] == 0.0)
					xx[i] = 0.0;
				else
					xx[i] = Math.log(x[i]);
				if (y[i] == 0.0)
					yy[i] = 0.0;
				else
					yy[i] = Math.log(y[i] / 100);
			}
			// setlength(tl,n+1,1);
			double tls = 0.0;
			double cf = 0.0;
			for (int j = 0; j <= n; j++) {
				tls = 0.0;
				for (int i = 0; i < x.length; i++)
					tls = tls + yy[i] * Math.pow(xx[i], j);
				tl[j][0] = tls;
			}

			double[][] coef = new double[n + 1][n + 1];

			for (int i = 0; i <= n; i++) {
				for (int j = 0; j <= n; j++) {
					cf = 0;
					for (int kk = 0; kk < x.length; kk++)
						cf = cf + Math.pow(xx[kk], i + j);
					coef[i][j] = cf;
				}
			}

			LinAEq.sysEqGauss(coef, tl, n + 1, 1);
		}

		/*
		 * System.out.println("-----------tlfinal----------------"); for (int
		 * j=0; j<tl.length; j++) System.out.println(tl[j][0]);
		 */
		// returneaza pentru ecuatia de interpolare in aceasta ordine:
		// a0 + a1x + a2x^2 + a3x^3 + ... + anx^n
		// returneaza solutia in tabloul tl!!
	}
}
