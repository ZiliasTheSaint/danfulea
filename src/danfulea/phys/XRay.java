package danfulea.phys;

import java.io.FileInputStream;

/**
 * Core class for X-Ray spectrum generator, HVL and filtration calculations!
 * 
 * @author Dan Fulea, 28 MAR. 2007
 */
public class XRay {
	/**
	 * Set this to 1 for using SRS78 database. 0 means using old approach. Default is 1.
	 */
	public static int ICALC = 1;// default with new (ripple); 0-means OLD!!!
	public static String filename = "";
	public static final String dataS = "Data";
	public static final String data_filterS = "FILTERS";
	public static final String data_kermaS = "KERMA";
	public static final String data_wS = "SPECW";
	public static final String data_mS = "SPECMO";
	public static final String data_rS = "SPECRH";
	public static final String data_constantS = "CONSTANT";
	public static final String data_rippleS = "RIPPLE";
	public static String file_sep = System.getProperty("file.separator");
	// ==================================================================
	public static final String defaultSpectrumExt = ".SPC";
	public static final String filterKermaExt = ".CSV";
	public static final String rippleExt_R0 = ".R0";
	public static final String rippleExt_R = ".R";

	public static int ianod = 0;// 0=W,1=Mo,2=Rh
	public static int iripple = 0;// 0,5,...,30
	public static int NEIN = 0;
	public static int MAXENE = 1000;
	public static double[] EIN = new double[MAXENE];// kev
	public static double[] XVAL = new double[MAXENE];
	public static int NEINK = 0;
	public static double[] EINK = new double[MAXENE];
	public static double[] KVAL = new double[MAXENE];// uGy/photons/mm2
	public static int NEINATT = 0;
	public static double[] EINATT = new double[MAXENE];
	public static double[] ATTVAL = new double[MAXENE];
	public static double[][] EATT = new double[MAXENE][20];// maxim 20 media!!!
	public static double[][] ATT = new double[MAXENE][20];// mm-1
	public static double DE = 0.0;// delta energy, compute at spectrum load
	public static double KVP = 0.0;
	public static int NOATT = 0;
	public static int NEFF = 0;// effectiv no of energies
	public static int NEFF_kv = 0;
	public static double[] YS = new double[MAXENE];// temp energy yield in
													// spectra AFTER attenuation
	public static double[] YF = new double[MAXENE];// final energy yield in
													// spectra AFTER attenuation
	public static double[] TMM = new double[20];// atenuator:e.g. 20 mm Al
	public static double[] YFF = new double[MAXENE];
	public static double[] YF_temp = new double[MAXENE];
	public static double[] YF_kv = new double[MAXENE];
	public static double HVL1 = 0.0;
	public static double HVL2 = 0.0;
	public static double HVL1_TISS = 0.0;
	public static double HVL2_TISS = 0.0;
	public static double HOMMFAC = 0.0;
	public static double MEANSPECTRUMENERGY = 0.0;
	public static double KERMAPERMASAT750MM = 0.0;
	public static double THI_HVL1 = 120.0;
	public static double THI_HVL2 = 160.0;
	public static double p1_HVL = 0.0;
	public static double p2_HVL = 0.0;
	public static double eqFiltr_HVL1 = 0.0;
	public static double eqFiltr_HVL2 = 0.0;
	public static double eqFiltr = 0.0;
	public static double HVL1_E = 0.0;
	public static double HVL2_E = 0.0;
	public static double HVL1_TISS_E = 0.0;
	public static double HVL2_TISS_E = 0.0;
	public static double PLANE = 20.0;// default patient plane!!!
	public static double HVL1_kv = 0.0;
	public static double HVL2_kv = 0.0;
	public static double HVL1_TISS_kv = 0.0;
	public static double HVL2_TISS_kv = 0.0;
	public static int ikv_file = 0;
	public static int ianod_file = 0;
	public static int NEIN_kv = 0;
	public static double[] EIN_kv = new double[MAXENE];// kev
	public static double[] XVAL_kv = new double[MAXENE];
	public static double DE_kv = 0.0;
	/**
	 * Total photons in spectrum, given in photons/mAs/mm^2 at 75 cm distance from the tube.
	 */
	public static double TOT_PHOTONS = 0.0;// photons/mas/mm2 at 750 mm
	public static boolean ISOK = true;
	public static int MAXITER = 10000;

	// ========================================================================================

	/**
	 * Computes HVL1 and HVL2 for a given voltage. The anode angle and equivalent filtration must be known! As always, the anode 
	 * material is considered to be aluminum because equivalent filtration is given in terms of mm Al.
	 * @param ikvp ikvp
	 */
	public static void computeHVL12Kvp(int ikvp) {
		ikv_file = ikvp;
		String uanodS = "";
		if (ianod_file >= 10)
			uanodS = ianod_file + "0";
		else
			uanodS = "0" + ianod_file + "0";
		String kvS = "";
		if (ikvp < 100)// >10
			kvS = "0" + ikv_file;
		else
			kvS = "" + ikv_file;
		String filename = kvS + uanodS;
		XRay.readSpectrum_kv(filename);
		buildSpectra_kv(ikvp, eqFiltr);// eqFiltr must be known!!!
		computeHVL1_kv("AL", false);// =>HVL1_kv
		computeHVL2_kv("AL", false);// =>HVL2_kv
	}

	/**
	 * Computes equivalent filtration (to be considered as total tube filtration) as weighted mean from filtrations derived from HVL1 and HVL2. 
	 * In theory, filtrations computed from HVL1 and HVL2 should match so this is not quite necessary. If they don't match it will raise the 
	 * question which one is more appropriate to be used as total tube filtration and/or how to establish the weights for a mean value? One 
	 * empirical way is to consider radiation absorbed in soft tissue and based of tissue thickness make a comparison  with 
	 * HVL (radiation intensity drop at 1/2) and QVL (radiation intensity drop at 1/4). HVL1 is simply HVL and HVL2 is QVL-HVL so they 
	 * are related. So, if filtrations don't match this routine performs a fine tunning based on assumption the radiation will reach the human body! As stated above, 
	 * the filtrations should match regardless if they are computed from HVL, QVL or TVL (radiation drop at 1/10) so be sure your experimental setup 
	 * for HVL and QVL measurements is appropriate, there are no voltage variations etc. If however you still get different results for 
	 * filtrations, use this method to estimate a total tube filtration for further calculations such as Monte Carlo simulations for dose assessment but 
	 * be advise, there is something wrong with the tested XRay tube.    
	 * @return the result
	 */
	public static double computeFiltration() {
		double result = 0.0;
		double hvl1e = 0.0;
		// double hvl2e=0.0;

		buildSpectra(eqFiltr_HVL1);// based on unatenuated spectrum=>YF_TEMP
		computeHVL1_temp("TISS", true);// mm TISS equivalent on eqFiltr_HVL1 or
										// on hvl1e
		hvl1e = HVL1_TISS_E;
		buildSpectra(eqFiltr_HVL2);// based on unatenuated spectrum=>YF_TEMP
		computeHVL1_temp("TISS", true);
		computeHVL2_temp("TISS", true);// mm TISS equivalent on eqFiltr_HVL2 or
										// on hvl2e
		// hvl2e=HVL2_TISS_E;

		p1_HVL = 0.5;
		p2_HVL = 0.5;
		if (hvl1e >= PLANE)
			p1_HVL = 1.0;
		else
			p1_HVL = 1 - (PLANE - hvl1e) / PLANE;
		p2_HVL = 1 - p1_HVL;

		result = (eqFiltr_HVL1 * p1_HVL + eqFiltr_HVL2 * p2_HVL)
				/ (p1_HVL + p2_HVL);
		// report as mmAl:
		eqFiltr = result;
		System.out.println("========>eq filtr in mmAl = " + eqFiltr);
		return result;
	}

	/**
	 * Based on experimental values for HVL1 and HVL2, computes tube filtration based on HVL2.
	 * @param hvl1e hvl1e
	 * @param hvl2e hvl2e
	 * @return the result
	 */
	public static double computeFiltrationFromHVL2(double hvl1e, double hvl2e) {
		double df = 0.0;
		double t1 = 0.0;
		double t2 = 0.0;// double t3=0.0;
		int kx = 0;
		double tt = 0.0;
		double tlo = 0.0;
		double thi = 0.0;// MAXIMUM SEAK!!!
		String filename = "AL";
		readAttCoef(filename);// =>ATTVAL
		double K0 = 0.0;
		double qvl = hvl1e + hvl2e;

		// compute first XRAY VALUE without Attenuation
		for (int i = 1; i <= NEFF; i++) {
			YS[i - 1] = XVAL[i - 1];
			// and after HVL1 attenuation
			if (ATTVAL[i - 1] * qvl > 1440) {
				YFF[i - 1] = 0.0;
			}
			YFF[i - 1] = YS[i - 1] * Math.exp(-ATTVAL[i - 1] * qvl);
		}

		ISOK = true;
		while (true) {
			t1 = 0.0;
			t2 = 0.0;
			kx = kx + 1;
			for (int l = 1; l <= NEFF; l++) {
				// update YS with HVL1 mmAl
				t1 = t1 + YS[l - 1] * KVAL[l - 1];// sum yiKermai which must be
													// K0!!
				t2 = t2 + YFF[l - 1] * KVAL[l - 1];// sum yiKermai which must be
													// K0/2.0!!
			}
			K0 = t1;
			if (kx == 1) {
				thi = THI_HVL2;// maximum equivalent filtration
			} else {
				if (t2 > K0 / 4.0) {
					thi = tt;
				} else {
					tlo = tt;
				}
				if (t2 < (K0 / 4.0) * 1.00005 && t2 > (K0 / 4.0) * 0.99995)
					break;
			}
			df = (tlo + thi) / 2.0;// 1st guess of total equivalent filtration
			tt = df;
			// System.out.println("f= "+df+"  t1  "+t1+"  t2  "+t2);
			for (int i = 1; i <= NEFF; i++) {
				if (ATTVAL[i - 1] * tt > 1440) {
					YS[i - 1] = 0.0;
				} else {
					YS[i - 1] = XVAL[i - 1] * Math.exp(-ATTVAL[i - 1] * tt);
				}

				if (ATTVAL[i - 1] * qvl > 1440) {
					YFF[i - 1] = 0.0;
				}
				YFF[i - 1] = YS[i - 1] * Math.exp(-ATTVAL[i - 1] * qvl);
			}
			if (kx > MAXITER) {
				ISOK = false;
				System.out.println("ERR % = " + 100 * (t2 - K0 / 4.0) / K0
						/ 4.0);
				break;
			}
		}
		eqFiltr_HVL2 = df;
		// System.out.println("uGy/mAs at 750 mm in air= "+KERMAPERMASAT750MM);
		// System.out.println("Mean energy in keV = "+MEANSPECTRUMENERGY);
		System.out.println("========>eq filtr2 in mmAl = " + df);

		return df;
	}

	/**
	 * Based on experimental values for HVL1, computes tube filtration based on HVL1.
	 * @param hvl1e hvl1e
	 * @return the result
	 */
	public static double computeFiltrationFromHVL1(double hvl1e) {
		double df = 0.0;
		double t1 = 0.0;
		double t2 = 0.0;// double t3=0.0;
		// MEANSPECTRUMENERGY=0.0;KERMAPERMASAT750MM=0.0;
		int kx = 0;
		double tt = 0.0;
		double tlo = 0.0;
		double thi = 0.0;// MAXIMUM SEAK!!!
		String filename = "AL";
		readAttCoef(filename);// =>ATTVAL
		double K0 = 0.0;

		ISOK = true;
		// compute first XRAY VALUE without Attenuation
		for (int i = 1; i <= NEFF; i++) {
			YS[i - 1] = XVAL[i - 1];
			// and after HVL1 attenuation
			if (ATTVAL[i - 1] * hvl1e > 1440) {
				YFF[i - 1] = 0.0;
			}
			YFF[i - 1] = YS[i - 1] * Math.exp(-ATTVAL[i - 1] * hvl1e);
		}

		while (true) {
			t1 = 0.0;
			t2 = 0.0;// t3=0.0;
			kx = kx + 1;
			for (int l = 1; l <= NEFF; l++) {
				// update YS with HVL1 mmAl
				t1 = t1 + YS[l - 1] * KVAL[l - 1];// sum yiKermai which must be
													// K0!!
				// t2=t2+YS[l-1]*DE*l;//sum yiEi
				// t3=t3+YS[l-1];//sum yi
				t2 = t2 + YFF[l - 1] * KVAL[l - 1];// sum yiKermai which must be
													// K0/2.0!!
			}
			K0 = t1;
			if (kx == 1) {
				// K0=t1;
				// KERMAPERMASAT750MM=t1;
				// MEANSPECTRUMENERGY=t2/t3;
				thi = THI_HVL1;// maximum equivalent filtration
			} else {
				// if(t2>K0/2.0){tlo=tt;}else{thi=tt;}
				if (t2 > K0 / 2.0) {
					thi = tt;
				} else {
					tlo = tt;
				}
				if (t2 < (K0 / 2.0) * 1.00005 && t2 > (K0 / 2.0) * 0.99995)
					break;
			}
			df = (tlo + thi) / 2.0;// 1st guess of total equivalent filtration
			tt = df;
			// System.out.println("f= "+df+"  t1  "+t1+"  t2  "+t2);
			for (int i = 1; i <= NEFF; i++) {
				if (ATTVAL[i - 1] * tt > 1440) {
					YS[i - 1] = 0.0;
				} else {
					YS[i - 1] = XVAL[i - 1] * Math.exp(-ATTVAL[i - 1] * tt);
				}

				if (ATTVAL[i - 1] * hvl1e > 1440) {
					YFF[i - 1] = 0.0;
				}
				YFF[i - 1] = YS[i - 1] * Math.exp(-ATTVAL[i - 1] * hvl1e);
			}
			if (kx > MAXITER) {
				ISOK = false;
				System.out.println("ERR % = " + 100 * (t2 - K0 / 2.0) / K0
						/ 2.0);
				break;
			}
		}
		eqFiltr_HVL1 = df;
		// System.out.println("uGy/mAs at 750 mm in air= "+KERMAPERMASAT750MM);
		// System.out.println("Mean energy in keV = "+MEANSPECTRUMENERGY);
		System.out.println("========>eq filtr1 in mmAl = " + df);

		return df;
	}

	/**
	 * Computes theoretical HVL2. The XRay theory provides data of unattenuated XRay spectrum based on voltage, anode angle, anode material and voltage waveform ripple. 
	 * The attenuator (absorber) is usually considered  to be equivalent aluminum.
	 * @param filename the filename containing attenuation coefficients of the absorber (usually Al)
	 * @param tissB if this is true than HVL2 in tissue is computed and for consistency the filename must be properly assigned (TISS).
	 */
	public static void computeHVL2(String filename, boolean tissB) {
		double t1 = 0.0;
		double t2 = 0.0;
		double t3 = 0.0;
		MEANSPECTRUMENERGY = 0.0;
		KERMAPERMASAT750MM = 0.0;
		int kx = 0;
		double tt = 0.0;
		double tlo = 0.0;
		double thi = 0.0;// MAXIMUM SEAK!!!
		// String filename="AL";
		readAttCoef(filename);// =>ATTVAL

		resetYS();

		while (true) {
			t1 = 0.0;
			t2 = 0.0;
			t3 = 0.0;
			kx = kx + 1;
			for (int l = 1; l <= NEFF; l++) {
				t1 = t1 + YS[l - 1] * KVAL[l - 1];// sum yiKermai
				t2 = t2 + YS[l - 1] * DE * l;// sum yiEi
				t3 = t3 + YS[l - 1];// sum yi
			}
			if (kx == 1) {
				KERMAPERMASAT750MM = t1;
				MEANSPECTRUMENERGY = t2 / t3;
				thi = THI_HVL2;
			} else {
				if (t1 > KERMAPERMASAT750MM / 4.0) {
					tlo = tt;
				} else {
					thi = tt;
				}
				if (t1 < (KERMAPERMASAT750MM / 4.0) * 1.00005
						&& t1 > (KERMAPERMASAT750MM / 4.0) * 0.99995)
					break;
			}
			if (tissB) {
				HVL2_TISS = (tlo + thi) / 2.0;// 1st guess
				tt = HVL2_TISS;
			} else {
				HVL2 = (tlo + thi) / 2.0;// 1st guess
				tt = HVL2;
			}
			// tt=HVL2;
			for (int i = 1; i <= NEFF; i++) {
				if (ATTVAL[i - 1] * tt > 1440) {
					YS[i - 1] = 0.0;
				} else {
					YS[i - 1] = YF[i - 1] * Math.exp(-ATTVAL[i - 1] * tt);
				}
			}

			if (kx > MAXITER) {
				break;
			}

		}
		if (tissB)
			HVL2_TISS = HVL2_TISS - HVL1_TISS;// QVL-HVL1!!!
		else
			HVL2 = HVL2 - HVL1;// QVL-HVL1!!!

		if (!tissB)
			HOMMFAC = HVL1 / HVL2;

		System.out.println("uGy/mAs at 750 mm in air= " + KERMAPERMASAT750MM);
		System.out.println("Mean energy in keV = " + MEANSPECTRUMENERGY);
		if (!tissB) {
			System.out.println("HVL2 in mm" + filename + " = " + HVL2);
			System.out.println("HOMMFAC = " + HOMMFAC);
		} else {
			System.out.println("HVL2 in mm" + filename + " = " + HVL2_TISS);
		}
	}

	/**
	 * Basically same as above and called by computeFiltrationn(). The reason for this routine is the usage of different global variables.
	 * @param filename the filename containing attenuation coefficients of the absorber (usually Al)
	 * @param tissB if this is true than HVL2 in tissue is computed and for consistency the filename must be properly assigned (TISS).
	 */
	public static void computeHVL2_temp(String filename, boolean tissB) {
		double t1 = 0.0;
		double t2 = 0.0;
		double t3 = 0.0;
		// double MEANSPECTRUMENERGY_temp=0.0;
		double KERMAPERMASAT750MM_temp = 0.0;
		int kx = 0;
		double tt = 0.0;
		double tlo = 0.0;
		double thi = 0.0;// MAXIMUM SEAK!!!
		// String filename="AL";
		readAttCoef(filename);// =>ATTVAL

		resetYS_temp();

		while (true) {
			t1 = 0.0;
			t2 = 0.0;
			t3 = 0.0;
			kx = kx + 1;
			for (int l = 1; l <= NEFF; l++) {
				t1 = t1 + YS[l - 1] * KVAL[l - 1];// sum yiKermai
				t2 = t2 + YS[l - 1] * DE * l;// sum yiEi
				t3 = t3 + YS[l - 1];// sum yi
			}
			if (kx == 1) {
				KERMAPERMASAT750MM_temp = t1;
				// MEANSPECTRUMENERGY_temp=t2/t3;
				thi = THI_HVL2;
			} else {
				if (t1 > KERMAPERMASAT750MM_temp / 4.0) {
					tlo = tt;
				} else {
					thi = tt;
				}
				if (t1 < (KERMAPERMASAT750MM_temp / 4.0) * 1.00005
						&& t1 > (KERMAPERMASAT750MM_temp / 4.0) * 0.99995)
					break;
			}
			if (tissB) {
				HVL2_TISS_E = (tlo + thi) / 2.0;// 1st guess
				tt = HVL2_TISS_E;
			} else {
				HVL2_E = (tlo + thi) / 2.0;// 1st guess
				tt = HVL2_E;
			}
			// tt=HVL2;
			for (int i = 1; i <= NEFF; i++) {
				if (ATTVAL[i - 1] * tt > 1440) {
					YS[i - 1] = 0.0;
				} else {
					YS[i - 1] = YF_temp[i - 1] * Math.exp(-ATTVAL[i - 1] * tt);
				}
			}
			if (kx > MAXITER) {
				break;
			}
		}
		if (tissB)
			HVL2_TISS_E = HVL2_TISS_E - HVL1_TISS_E;// QVL-HVL1!!!
		else
			HVL2_E = HVL2_E - HVL1_E;// QVL-HVL1!!!

		// if(!tissB)
		// HOMMFAC=HVL1/HVL2;

		// System.out.println("uGy/mAs at 750 mm in air= "+KERMAPERMASAT750MM);
		// System.out.println("Mean energy in keV = "+MEANSPECTRUMENERGY);
		if (!tissB) {
			System.out.println("HVL2_E in mm" + filename + " = " + HVL2_E);
			// System.out.println("HOMMFAC = "+HOMMFAC);
		} else {
			System.out.println("HVL2_E in mm" + filename + " = " + HVL2_TISS_E);
		}
	}

	/**
	 * Basically same as above and called by computeHVL12Kvp routine. The reason for this routine is the usage of different global variables.
	 * @param filename the filename containing attenuation coefficients of the absorber (usually Al)
	 * @param tissB if this is true than HVL2 in tissue is computed and for consistency the filename must be properly assigned (TISS).
	 
	 */
	public static void computeHVL2_kv(String filename, boolean tissB) {
		double t1 = 0.0;
		double t2 = 0.0;
		double t3 = 0.0;
		// double MEANSPECTRUMENERGY_temp=0.0;
		double KERMAPERMASAT750MM_temp = 0.0;
		int kx = 0;
		double tt = 0.0;
		double tlo = 0.0;
		double thi = 0.0;// MAXIMUM SEAK!!!
		// String filename="AL";
		readAttCoef(filename);// =>ATTVAL

		resetYS_kv();

		while (true) {
			t1 = 0.0;
			t2 = 0.0;
			t3 = 0.0;
			kx = kx + 1;
			for (int l = 1; l <= NEFF_kv; l++) {
				t1 = t1 + YS[l - 1] * KVAL[l - 1];// sum yiKermai
				t2 = t2 + YS[l - 1] * DE_kv * l;// sum yiEi
				t3 = t3 + YS[l - 1];// sum yi
			}
			if (kx == 1) {
				KERMAPERMASAT750MM_temp = t1;
				// MEANSPECTRUMENERGY_temp=t2/t3;
				thi = THI_HVL2;
			} else {
				if (t1 > KERMAPERMASAT750MM_temp / 4.0) {
					tlo = tt;
				} else {
					thi = tt;
				}
				if (t1 < (KERMAPERMASAT750MM_temp / 4.0) * 1.00005
						&& t1 > (KERMAPERMASAT750MM_temp / 4.0) * 0.99995)
					break;
			}
			if (tissB) {
				HVL2_TISS_kv = (tlo + thi) / 2.0;// 1st guess
				tt = HVL2_TISS_kv;
			} else {
				HVL2_kv = (tlo + thi) / 2.0;// 1st guess
				tt = HVL2_kv;
			}
			// tt=HVL2;
			for (int i = 1; i <= NEFF_kv; i++) {
				if (ATTVAL[i - 1] * tt > 1440) {
					YS[i - 1] = 0.0;
				} else {
					YS[i - 1] = YF_kv[i - 1] * Math.exp(-ATTVAL[i - 1] * tt);
				}
			}
			if (kx > MAXITER) {
				break;
			}
		}
		if (tissB)
			HVL2_TISS_kv = HVL2_TISS_kv - HVL1_TISS_kv;// QVL-HVL1!!!
		else
			HVL2_kv = HVL2_kv - HVL1_kv;// QVL-HVL1!!!

		// if(!tissB)
		// HOMMFAC=HVL1/HVL2;

		// System.out.println("uGy/mAs at 750 mm in air= "+KERMAPERMASAT750MM);
		// System.out.println("Mean energy in keV = "+MEANSPECTRUMENERGY);
		if (!tissB) {
			System.out.println("HVL2_kv in mm" + filename + " = " + HVL2_kv);
			// System.out.println("HOMMFAC = "+HOMMFAC);
		} else {
			System.out.println("HVL2_kv in mm" + filename + " = "
					+ HVL2_TISS_kv);
		}
	}

	/**
	 * Computes theoretical HVL1. The XRay theory provides data of unattenuated XRay spectrum based on voltage, anode angle, anode material and voltage waveform ripple. 
	 * The attenuator (absorber) is usually considered  to be equivalent aluminum.
	 * @param filename the filename containing attenuation coefficients of the absorber (usually Al)
	 * @param tissB if this is true than HVL2 in tissue is computed and for consistency the filename must be properly assigned (TISS).
	 */
	public static void computeHVL1(String filename, boolean tissB) {
		double t1 = 0.0;
		double t2 = 0.0;
		double t3 = 0.0;
		MEANSPECTRUMENERGY = 0.0;
		KERMAPERMASAT750MM = 0.0;
		int kx = 0;
		double tt = 0.0;
		double tlo = 0.0;
		double thi = 20.0;// MAXIMUM SEAK!!!
		// String filename="AL";
		readAttCoef(filename);// =>ATTVAL

		resetYS();

		while (true) {
			t1 = 0.0;
			t2 = 0.0;
			t3 = 0.0;
			kx = kx + 1;
			for (int l = 1; l <= NEFF; l++) {
				t1 = t1 + YS[l - 1] * KVAL[l - 1];// sum yiKermai
				t2 = t2 + YS[l - 1] * DE * l;// sum yiEi
				t3 = t3 + YS[l - 1];// sum yi
			}
			if (kx == 1) {
				KERMAPERMASAT750MM = t1;
				MEANSPECTRUMENERGY = t2 / t3;
				thi = THI_HVL1;
			} else {
				if (t1 > KERMAPERMASAT750MM / 2.0) {
					tlo = tt;
				} else {
					thi = tt;
				}
				if (t1 < (KERMAPERMASAT750MM / 2.0) * 1.00005
						&& t1 > (KERMAPERMASAT750MM / 2.0) * 0.99995)
					break;
			}
			if (tissB) {
				HVL1_TISS = (tlo + thi) / 2.0;// 1st guess
				tt = HVL1_TISS;
			} else {
				HVL1 = (tlo + thi) / 2.0;// 1st guess
				tt = HVL1;
			}
			// tt=HVL1;
			for (int i = 1; i <= NEFF; i++) {
				if (ATTVAL[i - 1] * tt > 1440) {
					YS[i - 1] = 0.0;
				} else {
					YS[i - 1] = YF[i - 1] * Math.exp(-ATTVAL[i - 1] * tt);
				}
			}
			if (kx > MAXITER) {
				break;
			}
		}
		System.out.println("uGy/mAs at 750 mm in air= " + KERMAPERMASAT750MM);
		System.out.println("Mean energy in keV = " + MEANSPECTRUMENERGY);
		if (!tissB)
			System.out.println("HVL1 in mm" + filename + " = " + HVL1);
		else
			System.out.println("HVL1 in mm" + filename + " = " + HVL1_TISS);
	}

	/**
	 * Basically same as above and called by computeFiltrationn(). The reason for this routine is the usage of different global variables.
	 * @param filename the filename containing attenuation coefficients of the absorber (usually Al)
	 * @param tissB if this is true than HVL2 in tissue is computed and for consistency the filename must be properly assigned (TISS).
	 */
	public static void computeHVL1_temp(String filename, boolean tissB) {
		double t1 = 0.0;
		double t2 = 0.0;
		double t3 = 0.0;
		// double MEANSPECTRUMENERGY_temp=0.0;
		double KERMAPERMASAT750MM_temp = 0.0;
		int kx = 0;
		double tt = 0.0;
		double tlo = 0.0;
		double thi = 20.0;// MAXIMUM SEAK!!!
		// String filename="AL";
		readAttCoef(filename);// =>ATTVAL

		resetYS_temp();

		while (true) {
			t1 = 0.0;
			t2 = 0.0;
			t3 = 0.0;
			kx = kx + 1;
			for (int l = 1; l <= NEFF; l++) {
				t1 = t1 + YS[l - 1] * KVAL[l - 1];// sum yiKermai
				t2 = t2 + YS[l - 1] * DE * l;// sum yiEi
				t3 = t3 + YS[l - 1];// sum yi
			}
			if (kx == 1) {
				KERMAPERMASAT750MM_temp = t1;
				// MEANSPECTRUMENERGY_temp=t2/t3;
				thi = THI_HVL1;
			} else {
				if (t1 > KERMAPERMASAT750MM_temp / 2.0) {
					tlo = tt;
				} else {
					thi = tt;
				}
				if (t1 < (KERMAPERMASAT750MM_temp / 2.0) * 1.00005
						&& t1 > (KERMAPERMASAT750MM_temp / 2.0) * 0.99995)
					break;
			}
			if (tissB) {
				HVL1_TISS_E = (tlo + thi) / 2.0;// 1st guess
				tt = HVL1_TISS_E;
			} else {
				HVL1_E = (tlo + thi) / 2.0;// 1st guess
				tt = HVL1_E;
			}
			// tt=HVL1;
			for (int i = 1; i <= NEFF; i++) {
				if (ATTVAL[i - 1] * tt > 1440) {
					YS[i - 1] = 0.0;
				} else {
					YS[i - 1] = YF_temp[i - 1] * Math.exp(-ATTVAL[i - 1] * tt);
				}
			}
			if (kx > MAXITER) {
				break;
			}
		}

		if (!tissB)
			System.out.println("HVL1_E in mm" + filename + " = " + HVL1_E);
		else
			System.out.println("HVL1_E in mm" + filename + " = " + HVL1_TISS_E);
	}

	/**
	 * Basically same as above and called by computeHVL12Kvp routine. The reason for this routine is the usage of different global variables.
	 * @param filename the filename containing attenuation coefficients of the absorber (usually Al)
	 * @param tissB if this is true than HVL2 in tissue is computed and for consistency the filename must be properly assigned (TISS).
	 
	 */
	public static void computeHVL1_kv(String filename, boolean tissB) {
		double t1 = 0.0;
		double t2 = 0.0;
		double t3 = 0.0;
		// double MEANSPECTRUMENERGY_temp=0.0;
		double KERMAPERMASAT750MM_temp = 0.0;
		int kx = 0;
		double tt = 0.0;
		double tlo = 0.0;
		double thi = 20.0;// MAXIMUM SEAK!!!
		// String filename="AL";
		readAttCoef(filename);// =>ATTVAL

		resetYS_kv();

		while (true) {
			t1 = 0.0;
			t2 = 0.0;
			t3 = 0.0;
			kx = kx + 1;
			for (int l = 1; l <= NEFF_kv; l++) {
				t1 = t1 + YS[l - 1] * KVAL[l - 1];// sum yiKermai
				t2 = t2 + YS[l - 1] * DE_kv * l;// sum yiEi
				t3 = t3 + YS[l - 1];// sum yi
			}
			if (kx == 1) {
				KERMAPERMASAT750MM_temp = t1;
				// MEANSPECTRUMENERGY_temp=t2/t3;
				thi = THI_HVL1;
			} else {
				if (t1 > KERMAPERMASAT750MM_temp / 2.0) {
					tlo = tt;
				} else {
					thi = tt;
				}
				if (t1 < (KERMAPERMASAT750MM_temp / 2.0) * 1.00005
						&& t1 > (KERMAPERMASAT750MM_temp / 2.0) * 0.99995)
					break;
			}
			if (tissB) {
				HVL1_TISS_kv = (tlo + thi) / 2.0;// 1st guess
				tt = HVL1_TISS_kv;
			} else {
				HVL1_kv = (tlo + thi) / 2.0;// 1st guess
				tt = HVL1_kv;
			}
			// tt=HVL1;
			for (int i = 1; i <= NEFF_kv; i++) {
				if (ATTVAL[i - 1] * tt > 1440) {
					YS[i - 1] = 0.0;
				} else {
					YS[i - 1] = YF_kv[i - 1] * Math.exp(-ATTVAL[i - 1] * tt);
				}
			}
			if (kx > MAXITER) {
				break;
			}
		}

		if (!tissB)
			System.out.println("HVL1_kv in mm" + filename + " = " + HVL1_kv);
		else
			System.out.println("HVL1_kv in mm" + filename + " = "
					+ HVL1_TISS_kv);
	}

	/**
	 * Internally used.
	 */
	private static void resetYS() {
		for (int i = 1; i <= NEFF; i++) {
			YS[i - 1] = YF[i - 1];
		}
	}

	/**
	 * Internally used.
	 */
	private static void resetYS_temp() {
		for (int i = 1; i <= NEFF; i++) {
			YS[i - 1] = YF_temp[i - 1];
		}
	}

	/**
	 * Internally used.
	 */
	private static void resetYS_kv() {
		for (int i = 1; i <= NEFF_kv; i++) {
			YS[i - 1] = YF_kv[i - 1];
		}
	}

	/**
	 * Build the X-Ray spectrum by initializing global variables.
	 */
	public static void buildSpectra() {
		// initial YS==at spectrum Read
		Double d = new Double(KVP / DE);// number of energies!!
		int N = d.intValue();
		NEFF = N;
		for (int i = 1; i <= N; i++) {
			YF[i - 1] = 1.0;
		}// reset
		for (int j = 1; j <= NOATT; j++) {
			for (int i = 1; i <= N; i++) {
				if (ATT[i - 1][j - 1] * TMM[j - 1] > 1440) {
					YS[i - 1] = 0.0;
				}
				// else{YS[i-1]=YS[i-1]*Math.exp(-ATT[i-1][j-1]*TMM[j-1]);}
				else {
					YS[i - 1] = XVAL[i - 1]
							* Math.exp(-ATT[i - 1][j - 1] * TMM[j - 1]);
				}

				// YF[i-1]=YS[i-1];
				if (j == 1) {
					YF[i - 1] = YS[i - 1];
				} else {
					if (XVAL[i - 1] != 0.0)
						YF[i - 1] = YF[i - 1] * YS[i - 1] / XVAL[i - 1];// @@@@@@@@@@@@@@@@@@@@@@
					else
						YF[i - 1] = 0.0;
				}
			}
		}

		TOT_PHOTONS = 0.0;
		for (int i = 1; i <= N; i++) {
			TOT_PHOTONS = TOT_PHOTONS + YF[i - 1];
		}
	}

	/**
	 * Build X-Ray spectrum with a given filtration (total filtration in mm Al)
	 * @param filtr filtr
	 */
	public static void buildSpectra(double filtr) {
		// initial YS==at spectrum Read
		Double d = new Double(KVP / DE);// number of energies!!
		int N = d.intValue();
		// NEFF=N;

		String filename = "AL";
		readAttCoef(filename);// =>ATTVAL

		// for (int j=1;j<=NOATT;j++)
		// {
		for (int i = 1; i <= N; i++) {
			// if(ATT[i-1][j-1]*filtr>1440){YS[i-1]=0.0;}
			// else{YS[i-1]=XVAL[i-1]*Math.exp(-ATT[i-1][j-1]*filtr);}
			if (ATTVAL[i - 1] * filtr > 1440) {
				YS[i - 1] = 0.0;
			} else {
				YS[i - 1] = XVAL[i - 1] * Math.exp(-ATTVAL[i - 1] * filtr);
			}

			YF_temp[i - 1] = YS[i - 1];
		}
		// }

		// for (int
		// i=1;i<=N;i++){System.out.println("e: "+EIN[i-1]+" yf: "+YF[i-1]);}
	}

	/**
	 * Build X-Ray spectrum at desired voltage with a given filtration (total filtration in mm Al). Used by computeHVL12Kvp routine.
	 * @param ikvp ikvp
	 * @param filtr filtr
	 */
	public static void buildSpectra_kv(int ikvp, double filtr) {
		// initial YS==at spectrum Read
		Double d = new Double(ikvp / DE_kv);// number of energies!!
		int N = d.intValue();
		NEFF_kv = N;
		// System.out.println("neffkv="+NEFF_kv);
		String filename = "AL";
		readAttCoef(filename);// =>ATTVAL

		// for (int j=1;j<=NOATT;j++)
		// {
		for (int i = 1; i <= N; i++) {
			// if(ATT[i-1][j-1]*filtr>1440){YS[i-1]=0.0;}
			// else{YS[i-1]=XVAL[i-1]*Math.exp(-ATT[i-1][j-1]*filtr);}
			if (ATTVAL[i - 1] * filtr > 1440) {
				YS[i - 1] = 0.0;
			} else {
				YS[i - 1] = XVAL_kv[i - 1] * Math.exp(-ATTVAL[i - 1] * filtr);
			}

			YF_kv[i - 1] = YS[i - 1];
		}
		// }

		// for (int
		// i=1;i<=N;i++){System.out.println("e: "+EIN[i-1]+" yf: "+YF[i-1]);}
	}

	// =========================================================================================
	/**
	 * Reset global variables for re-use.
	 */
	public static void reset() {
		NEIN = 0;
		EIN = new double[MAXENE];
		XVAL = new double[MAXENE];
		NEINK = 0;
		EINK = new double[MAXENE];
		KVAL = new double[MAXENE];
		NEINATT = 0;
		EINATT = new double[MAXENE];
		ATTVAL = new double[MAXENE];
		EATT = new double[MAXENE][20];
		ATT = new double[MAXENE][20];
		DE = 0.0;
		KVP = 0.0;
		NOATT = 0;
		YS = new double[MAXENE];
		YF = new double[MAXENE];
		YFF = new double[MAXENE];
		HVL1 = 0.0;
		HVL2 = 0.0;
		HOMMFAC = 0.0;
		HVL1_TISS = 0.0;
		HVL2_TISS = 0.0;
		MEANSPECTRUMENERGY = 0.0;
		KERMAPERMASAT750MM = 0.0;
		NEFF = 0;
		p1_HVL = 0.0;
		p2_HVL = 0.0;
		eqFiltr_HVL1 = 0.0;
		eqFiltr_HVL2 = 0.0;
		YF_temp = new double[MAXENE];
		HVL1_E = 0.0;
		HVL2_E = 0.0;
		HVL1_TISS_E = 0.0;
		HVL2_TISS_E = 0.0;
		eqFiltr = 0.0;
		YF_kv = new double[MAXENE];
		NEFF_kv = 0;
		HVL1_kv = 0.0;
		HVL2_kv = 0.0;
		HVL1_TISS_kv = 0.0;
		HVL2_TISS_kv = 0.0;
		NEIN_kv = 0;
		EIN_kv = new double[MAXENE];
		XVAL_kv = new double[MAXENE];
		DE_kv = 0.0;
	}

	/**
	 * Convert a string to an int.
	 * @param value value
	 * @return the result
	 * @throws NumberFormatException may throw this exception
	 */
	public static int stringToInt(String value) throws NumberFormatException {
		value = value.trim();
		if (value.length() == 0) {
			return 0;
		} else {
			return Integer.parseInt(value);
		}
	}

	/**
	 * Convert an int to a String.
	 * @param i i
	 * @return the result
	 */
	public static String intToString(int i) {
		return Integer.toString(i);
	}

	/**
	 * Convert a string to a float.
	 * @param value value
	 * @return the result
	 * @throws NumberFormatException may throw this exception
	 */
	public static float stringToFloat(String value)
			throws NumberFormatException {
		value = value.trim();
		if (value.length() == 0) {
			return 0;
		} else {
			return Float.parseFloat(value);
		}
	}

	/**
	 * Convert a float to a String.
	 * @param f f
	 * @return the result
	 */
	public static String floatToString(float f) {
		return Float.toString(f);
	}

	/**
	 * Convert a string to a double.
	 * @param value value
	 * @return the result
	 * @throws NumberFormatException may throw this exception
	 */
	public static double stringToDouble(String value)
			throws NumberFormatException {
		value = value.trim();
		if (value.length() == 0) {
			return 0;
		} else {
			return Double.parseDouble(value);
		}
	}

	/**
	 * Convert a double to a String.
	 * @param d d
	 * @return the result
	 */
	public static String doubleToString(double d) {
		return Double.toString(d);
	}

	/**
	 * Read the spectrum from a file and initialize some global variables.
	 * @param filename filename
	 */
	public static void readSpectrum(String filename) {
		String filenam = "";// datas+file_sep+dataspec+file_sep+enerfile+defaultext;
		String ext = "";
		String subdir = "";
		String dir = "";

		if (iripple == 0) {
			ext = defaultSpectrumExt;
			subdir = data_constantS;
		} else if (iripple < 10) {
			ext = rippleExt_R0 + iripple;
			subdir = data_rippleS;
		} else {
			ext = rippleExt_R + iripple;
			subdir = data_rippleS;
		}
		if (ianod == 0) {
			dir = data_wS;
		} else if (ianod == 1) {
			dir = data_mS;
		} else if (ianod == 2) {
			dir = data_rS;
		}

		filenam = dataS + file_sep + dir + file_sep + subdir + file_sep
				+ filename + ext;
		// System.out.println(filenam);
		int iread = 0;
		int lnrr = 0;// line number
		int indx = 1;
		NEIN = 0;
		StringBuffer desc = new StringBuffer();
		boolean haveData = false;
		// String equals="=";
		char comma = ',';
		char lineSep = '\n';// System.getProperty("line.separator").charAt(0);

		boolean enB = false;
		boolean srB = false;

		try {
			FileInputStream in = new FileInputStream(filenam);

			while ((iread = in.read()) != -1) {
				if (!Character.isWhitespace((char) iread)
						&& ((char) iread != comma)) {
					desc.append((char) iread);
					haveData = true;
				} else {
					if (haveData)// we have data
					{
						haveData = false;// reset
						if (lnrr >= 1) {
							String s = desc.toString();
							if (!enB) {
								EIN[indx - 1] = stringToDouble(s);
								enB = true;
							} else if (!srB) {
								XVAL[indx - 1] = stringToDouble(s);
								YS[indx - 1] = XVAL[indx - 1];
								// System.out.println(EIN[indx-1]+"    "+XVAL[indx-1]);
								indx++;
								NEIN++;
								srB = true;
								enB = false;
								srB = false;
							}
						}

					}// have data

					if ((char) iread == lineSep)
						lnrr++;
					desc = new StringBuffer();
				}// else
			}// main while
			in.close();
			// System.out.println("  nein  "+NEIN);
		}// try
		catch (Exception exc) {
		}

		DE = EIN[1] - EIN[0];
	}

	/**
	 * Same as above but populates different global variables. Called by computeHVL12Kvp routine.
	 * @param filename filename
	 */
	public static void readSpectrum_kv(String filename) {
		String filenam = "";// datas+file_sep+dataspec+file_sep+enerfile+defaultext;
		String ext = "";
		String subdir = "";
		String dir = "";

		if (iripple == 0) {
			ext = defaultSpectrumExt;
			subdir = data_constantS;
		} else if (iripple < 10) {
			ext = rippleExt_R0 + iripple;
			subdir = data_rippleS;
		} else {
			ext = rippleExt_R + iripple;
			subdir = data_rippleS;
		}
		if (ianod == 0) {
			dir = data_wS;
		} else if (ianod == 1) {
			dir = data_mS;
		} else if (ianod == 2) {
			dir = data_rS;
		}

		filenam = dataS + file_sep + dir + file_sep + subdir + file_sep
				+ filename + ext;

		int iread = 0;
		int lnrr = 0;// line number
		int indx = 1;
		NEIN_kv = 0;
		StringBuffer desc = new StringBuffer();
		boolean haveData = false;
		// String equals="=";
		char comma = ',';
		char lineSep = '\n';// System.getProperty("line.separator").charAt(0);

		boolean enB = false;
		boolean srB = false;

		try {
			FileInputStream in = new FileInputStream(filenam);

			while ((iread = in.read()) != -1) {
				if (!Character.isWhitespace((char) iread)
						&& ((char) iread != comma)) {
					desc.append((char) iread);
					haveData = true;
				} else {
					if (haveData)// we have data
					{
						haveData = false;// reset
						if (lnrr >= 1) {
							String s = desc.toString();
							if (!enB) {
								EIN_kv[indx - 1] = stringToDouble(s);
								enB = true;
							} else if (!srB) {
								XVAL_kv[indx - 1] = stringToDouble(s);
								YS[indx - 1] = XVAL_kv[indx - 1];
								// System.out.println(EIN[indx-1]+"    "+XVAL[indx-1]);
								indx++;
								NEIN_kv++;
								srB = true;
								enB = false;
								srB = false;
							}
						}

					}// have data

					if ((char) iread == lineSep)
						lnrr++;
					desc = new StringBuffer();
				}// else
			}// main while
			in.close();
			// System.out.println("  nein  "+NEIN);
		}// try
		catch (Exception exc) {
		}

		DE_kv = EIN_kv[1] - EIN_kv[0];
	}

	/**
	 * Read kerma constants from a file.
	 * @param filename filename
	 */
	public static void readKerma(String filename) {
		String filenam = "";// datas+file_sep+dataspec+file_sep+enerfile+defaultext;
		// String ext="";
		// String subdir="";
		// String dir="";

		filenam = dataS + file_sep + data_kermaS + file_sep + filename
				+ filterKermaExt;

		int iread = 0;
		@SuppressWarnings("unused")
		int lnrr = 0;// line number
		int indx = 1;
		NEINK = 0;
		StringBuffer desc = new StringBuffer();
		boolean haveData = false;
		// String equals="=";
		char comma = ',';
		char lineSep = '\n';// System.getProperty("line.separator").charAt(0);

		boolean enB = false;
		boolean srB = false;

		try {
			FileInputStream in = new FileInputStream(filenam);

			while ((iread = in.read()) != -1) {
				if (!Character.isWhitespace((char) iread)
						&& ((char) iread != comma)) {
					desc.append((char) iread);
					haveData = true;
				} else {

					if (haveData)// we have data
					{
						haveData = false;// reset
						// if(lnrr>1)
						// {
						String s = desc.toString();
						if (!enB) {
							EINK[indx - 1] = stringToDouble(s);
							enB = true;
						} else if (!srB) {
							KVAL[indx - 1] = stringToDouble(s);
							// System.out.println(EINK[indx-1]+"    "+KVAL[indx-1]);
							indx++;
							NEINK++;
							srB = true;
							enB = false;
							srB = false;
						}
						// }

					}// have data

					if ((char) iread == lineSep)
						lnrr++;
					desc = new StringBuffer();
				}// else
			}// main while
			in.close();
			// System.out.println("  nein  "+NEINK);
		}// try
		catch (Exception exc) {
		}
	}

	/**
	 * Read attenuation coefficients from a file
	 * @param filename filename
	 */
	public static void readAttCoef(String filename) {
		String filenam = "";// datas+file_sep+dataspec+file_sep+enerfile+defaultext;
		// String ext="";
		// String subdir="";
		// String dir="";

		filenam = dataS + file_sep + data_filterS + file_sep + filename
				+ filterKermaExt;

		int iread = 0;
		@SuppressWarnings("unused")
		int lnrr = 0;// line number
		int indx = 1;
		NEINATT = 0;
		StringBuffer desc = new StringBuffer();
		boolean haveData = false;
		// String equals="=";
		char comma = ',';
		char lineSep = '\n';// System.getProperty("line.separator").charAt(0);

		boolean enB = false;
		boolean srB = false;

		try {
			FileInputStream in = new FileInputStream(filenam);

			while ((iread = in.read()) != -1) {
				if (!Character.isWhitespace((char) iread)
						&& ((char) iread != comma)) {
					desc.append((char) iread);
					haveData = true;
				} else {

					if (haveData)// we have data
					{
						haveData = false;// reset
						// if(lnrr>1)
						// {
						String s = desc.toString();
						if (!enB) {
							EINATT[indx - 1] = stringToDouble(s);
							enB = true;
						} else if (!srB) {
							ATTVAL[indx - 1] = stringToDouble(s);
							// System.out.println(EINATT[indx-1]+"    "+ATTVAL[indx-1]);
							indx++;
							NEINATT++;
							srB = true;
							enB = false;
							srB = false;
						}
						// }

					}// have data

					if ((char) iread == lineSep)
						lnrr++;
					desc = new StringBuffer();
				}// else
			}// main while
			in.close();
			// System.out.println("  nein  "+NEINATT);
		}// try
		catch (Exception exc) {
		}
	}

	/**
	 * Read attenuation coefficients from a file and set ATT matrix values based on j index. Useful when we have multiple attenuators.
	 * @param filename filename
	 * @param j j
	 */
	public static void readAttCoef(String filename, int j) {
		NOATT++;// each time attcoef are read!

		String filenam = "";// datas+file_sep+dataspec+file_sep+enerfile+defaultext;
		// String ext="";
		// String subdir="";
		// String dir="";

		filenam = dataS + file_sep + data_filterS + file_sep + filename
				+ filterKermaExt;

		int iread = 0;
		@SuppressWarnings("unused")
		int lnrr = 0;// line number
		int indx = 1;
		NEINATT = 0;
		StringBuffer desc = new StringBuffer();
		boolean haveData = false;
		// String equals="=";
		char comma = ',';
		char lineSep = '\n';// System.getProperty("line.separator").charAt(0);

		boolean enB = false;
		boolean srB = false;

		try {
			FileInputStream in = new FileInputStream(filenam);

			while ((iread = in.read()) != -1) {
				if (!Character.isWhitespace((char) iread)
						&& ((char) iread != comma)) {
					desc.append((char) iread);
					haveData = true;
				} else {

					if (haveData)// we have data
					{
						haveData = false;// reset
						// if(lnrr>1)
						// {
						String s = desc.toString();
						if (!enB) {
							EATT[indx - 1][j - 1] = stringToDouble(s);
							enB = true;
						} else if (!srB) {
							ATT[indx - 1][j - 1] = stringToDouble(s);
							// System.out.println(EATT[indx-1][j-1]+"    "+ATT[indx-1][j-1]);
							indx++;
							NEINATT++;
							srB = true;
							enB = false;
							srB = false;
						}
						// }

					}// have data

					if ((char) iread == lineSep)
						lnrr++;
					desc = new StringBuffer();
				}// else
			}// main while
			in.close();
			// System.out.println("  nein  "+NEINATT);
		}// try
		catch (Exception exc) {
			NOATT--;
		}
	}
}
