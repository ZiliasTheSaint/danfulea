package danfulea.phys;

import java.util.Vector;

import danfulea.math.StatsUtil;
import danfulea.math.numerical.Function;
import danfulea.math.numerical.ModelingData;

/**
 * Handling gamma ROI (gamma spectrometry). <br>
 * 
 * 
 * @author Dan Fulea, 23 APR. 2011.
 * 
 */
public class GammaRoi implements Function {
	private static double en_a3 = 0.0;
	private static double en_a2 = 0.0;
	private static double en_a1 = 0.0;
	private static double en_a0 = 0.0;
	//private static double ch_a3 = 0.0;
	//private static double ch_a2 = 0.0;
	//private static double ch_a1 = 0.0;
	//private static double ch_a0 = 0.0;
	private static double fwhm_a3 = 0.0;
	private static double fwhm_a2 = 0.0;
	private static double fwhm_a1 = 0.0;
	private static double fwhm_a0 = 0.0;
	private static double fwhm_overallProcentualError = 0.0;
	private static double eff_p1_a4 = 0.0;
	private static double eff_p1_a3 = 0.0;
	private static double eff_p1_a2 = 0.0;
	private static double eff_p1_a1 = 0.0;
	private static double eff_p1_a0 = 0.0;
	private static double eff_p2_a4 = 0.0;
	private static double eff_p2_a3 = 0.0;
	private static double eff_p2_a2 = 0.0;
	private static double eff_p2_a1 = 0.0;
	private static double eff_p2_a0 = 0.0;
	private static double eff_crossoverEnergy = 0.0;
	private static double eff_overallProcentualError = 0.0;

	private Vector<Double> channelV, pulsesV, bkgpulsesV = null;
	private double time = 1.0;
	// start, end channel is related to a valid gamma ROI
	// by default ROI edges, which is the lower and the upper
	// limit of ROI (the continuum background or the Compton background),
	// are set to those values.
	private double startChannel = 0.0;// in fact it is an integer
	private double startEnergy = 0.0;
	private double endChannel = 0.0;
	private double endEnergy = 0.0;
	private double centerChannel = 0.0;
	private double centerEnergy = 0.0;

	private double startEdgeChannel = 0.0;// in fact it is an integer
	private double startEdgeEnergy = 0.0;
	private double endEdgeChannel = 0.0;
	private double endEdgePulses = 0.0;
	private double startEdgePulses = 0.0;

	private double centroidChannel = 0.0;
	private double centroidEnergy = 0.0;

	private double peakChannel = 0.0;
	private double peakEnergy = 0.0;
	private double peakPulses = 0.0;
	private double endEdgeEnergy = 0.0;

	private double fwhmChannel = 0.0;
	private double fwhmChannelError = 0.0;
	private double fwhmEnergy = 0.0;
	private double fwhmEnergyError = 0.0;
	private double fwhmEnergyCalib = 0.0;
	private double fwhmEnergyCalibError = 0.0;
	private double resolution = 0.0;
	private double resolutionCalib = 0.0;
	private String significance = "No";
	private double bkgCounts = 0.0;
	private double bkgCountsError = 0.0;
	private double bkgCountsRate = 0.0;
	private double bkgCountsRateError = 0.0;
	private double grossCounts = 0.0;
	private double grossCountsError = 0.0;
	private double grossCountsRate = 0.0;
	private double grossCountsRateError = 0.0;
	private double comptonCounts = 0.0;
	private double comptonCountsError = 0.0;
	private double comptonCountsRate = 0.0;
	private double comptonCountsRateError = 0.0;
	private double netCounts = 0.0;
	private double netCountsError = 0.0;
	private double netCountsRate = 0.0;
	private double netCountsRateError = 0.0;
	private double confidenceLevel = 0.0;
	private String nuclide = "NoName";// NoName nuclide!
	private double yield = 0.0;
	private double efficiencyProcentual = 0.0;
	private double efficiencyProcentualError = 0.0;
	private double activity_Bq = 0.0;
	private double activity_BqError = 0.0;
	private double mda_Bq = 0.0;
	private double mda_BqError = 0.0;
	private String difference = "No";
	private double atomicMass = 0.0;
	private double halfLife = 0.0;
	private String halfLifeUnits = "s";// seconds
	private String netCalculationMethod = "Default";// Ge, GaussianFit
	private String mdaCalculationMethod = "Default";// Pasternack or Curie

	private int iMDA = 2;
	public static final int MDA_CALCULATION_PASTERNACK = 0;
	public static final int MDA_CALCULATION_CURIE = 1;
	public static final int MDA_CALCULATION_DEFAULT = 2;

	private int iNet = 0;
	public static final int NET_CALCULATION_NAI = 0;
	public static final int NET_CALCULATION_GE = 1;
	public static final int NET_CALCULATION_GAUSSIAN = 2;

	private double[] channel;
	private double[] pulses;

	private static final int MAXITER_mrqmin = 10000;
	public static AmbientalBkgRetriever abg;
	private boolean roiSet=false;

	/**
	 * Constructor. Scenario: The calling program makes use of insertRoi method and therefore 
	 * a new ROI (Region of Interest) is build (setup ROI around a peak). The calling program can also invoke updateRoiEdge 
	 * to set ROI edges at convenient places in order to build the continuum (Compton) background. 
	 * This is useful when we have multiple peaks cluttering together and want deconvolution (peak folding) by 
	 * selecting a ROI having few points (around peak) and set ROI edges far away and then perform Gaussian fit. 
	 * The calling program may also set ROI to desired nuclide and read energies and yields from a library. Also it 
	 * must transfer calibrations here (in energy, FWHM as well as efficiency). With all this, the 
	 * ROI is fully described and all ROI specific calculations can be easily done. 
	 *  
	 * 
	 * @param startChannel start channel
	 * @param endChannel end channel
	 */
	public GammaRoi(double startChannel, double endChannel) {
		channelV = new Vector<Double>();
		pulsesV = new Vector<Double>();
		bkgpulsesV = new Vector<Double>();
		setStartChannel(startChannel);
		setEndChannel(endChannel);
		setCenterChannel();
		// default
		setStartEdgeChannel(startChannel);
		setEndEdgeChannel(endChannel);
		// printCalib();//OK!
	}
	
	/**
	 * Reset the ambient background pulses 
	 */
	public void resetBkgData(){
		bkgpulsesV = new Vector<Double>();
	}

	// ====================
	/**
	 * 
	 * @return the pulses stored in startEgde channel
	 */
	public double getStartEdgePulses() {
		return this.startEdgePulses;
	}

	/**
	 * Set pulses at startEdge channel
	 * @param d d
	 */
	public void setStartEdgePulses(double d) {
		this.startEdgePulses = d;
	}

	// ---------------
	/**
	 * 
	 * @return the pulses stored in endEgde channel
	 */
	public double getEndEdgePulses() {
		return this.endEdgePulses;
	}

	/**
	 * Set pulses at endEdge channel
	 * @param d d
	 */
	public void setEndEdgePulses(double d) {
		this.endEdgePulses = d;
	}

	// =======================
	/**
	 * Set spectrum live time (take into account the detector dead-time, i.e. the time the detector does 
	 * not count any pulses; it is usually the acquisition time set by user - the gamma device take care 
	 * automatically of dead time by running a relatively longer real time. Here we set the LIVE time NOT 
	 * the REAL time).
	 * @param d d
	 */
	public void setLiveTime(double d) {
		time = d;
	}

	/**
	 * 
	 * @return the spectrum live time
	 */
	public double getLiveTime() {
		return time;
	}

	/**
	 * Add a channel number to channel vector.
	 * @param d the channel - an integer but treated as a double
	 */
	public void addChannel(double d) {
		channelV.addElement(new Double(d));
	}

	/**
	 * Add pulses corresponding to a channel to the pulses vector.
	 * @param d the pulses
	 */
	public void addPulses(double d) {
		pulsesV.addElement(new Double(d));
	}

	/**
	 * Add ambient background pulses corresponding to a channel to the background pulses vector.
	 * @param d the pulses (corrected to sample live time)
	 */
	public void addBkgPulses(double d) {
		bkgpulsesV.addElement(new Double(d));
	}

	/**
	 * Print to the console the in-use calibration coefficients.
	 */
	public void printCalib() {
		System.out.println("En " + en_a3 + " , " + en_a2 + " , " + en_a1
				+ " , " + en_a0);
		//System.out.println("Ch " + ch_a3 + " , " + ch_a2 + " , " + ch_a1
		//		+ " , " + ch_a0);
		System.out.println("FWHM " + fwhm_a3 + " , " + fwhm_a2 + " , "
				+ fwhm_a1 + " , " + fwhm_a0 + " , "
				+ fwhm_overallProcentualError);
		System.out.println("Eff1 " + eff_p1_a4 + " , " + eff_p1_a3 + " , "
				+ eff_p1_a2 + " , " + eff_p1_a1 + " , " + eff_p1_a0);
		System.out.println("Eff2 " + eff_p2_a4 + " , " + eff_p2_a3 + " , "
				+ eff_p2_a2 + " , " + eff_p2_a1 + " , " + eff_p2_a0);
		System.out.println(eff_crossoverEnergy + " ;err= "
				+ eff_overallProcentualError + "   " + channelV.size());

		for (int i = 0; i < channelV.size(); i++) {
			System.out.println("ch = " + channelV.elementAt(i) + " ;pulses = "
					+ pulsesV.elementAt(i));// +
			// " ;bkgpulses = "+bkgpulsesV.elementAt(i));
		}

		System.out.println("st " + startChannel + "  end " + endChannel);

		System.out.println("st edge " + startEdgeChannel + " ;endEdgeChannel "
				+ endEdgeChannel);
		System.out.println("stP edge " + startEdgePulses + " ;endEdgePulses "
				+ endEdgePulses);
	}
	
	/**
	 * Set energy calibration coefficients, energy [keV] = f(channel)
	 * @param e_a3 e_a3
	 * @param e_a2 e_a2
	 * @param e_a1 e_a1
	 * @param e_a0 e_a0
	 */
	public static void setEnergyCalibration(double e_a3, double e_a2, double e_a1,
			double e_a0){
		en_a3 = e_a3;
		en_a2 = e_a2;
		en_a1 = e_a1;
		en_a0 = e_a0;
	}

	/**
	 * Set FWHM calibration coefficients (FWHM = Full Width at Half Maximum), FWHM = f(centroid channel)
	 * @param f_a3 f_a3
	 * @param f_a2 f_a2
	 * @param f_a1 f_a1
	 * @param f_a0 f_a0
	 * @param f_overallProcentualError f_overallProcentualError
	 */
	public static void setFWHMCalibration(double f_a3, double f_a2, double f_a1, double f_a0,
			double f_overallProcentualError){
		fwhm_a3 = f_a3;
		fwhm_a2 = f_a2;
		fwhm_a1 = f_a1;
		fwhm_a0 = f_a0;
		fwhm_overallProcentualError = f_overallProcentualError;
	}
	
	/**
	 * Set efficiency calibration coefficients. The calibration uses two polynomials fit intersecting at a 
	 * crossover energy although in most cases a single polynomial and a crossover energy at 0 is good enough. 
	 * Ln(Eff) = f(Ln(centroid energy)). Energy is in keV, efficiency is normalized to 100 (a number less than 100, i.e. in %).
	 * @param e_p1_a4 e_p1_a4
	 * @param e_p1_a3 e_p1_a3
	 * @param e_p1_a2 e_p1_a2
	 * @param e_p1_a1 e_p1_a1
	 * @param e_p1_a0 e_p1_a0
	 * @param e_p2_a4 e_p2_a4
	 * @param e_p2_a3 e_p2_a3
	 * @param e_p2_a2 e_p2_a2
	 * @param e_p2_a1 e_p2_a1
	 * @param e_p2_a0 e_p2_a0
	 * @param e_crossoverEnergy e_crossoverEnergy
	 * @param e_overallProcentualError e_overallProcentualError
	 */
	public static void setEfficiencyCalibration(double e_p1_a4, double e_p1_a3,
			double e_p1_a2, double e_p1_a1, double e_p1_a0, double e_p2_a4,
			double e_p2_a3, double e_p2_a2, double e_p2_a1, double e_p2_a0,
			double e_crossoverEnergy, double e_overallProcentualError){
		eff_p1_a4 = e_p1_a4;
		eff_p1_a3 = e_p1_a3;
		eff_p1_a2 = e_p1_a2;
		eff_p1_a1 = e_p1_a1;
		eff_p1_a0 = e_p1_a0;
		eff_p2_a4 = e_p2_a4;
		eff_p2_a3 = e_p2_a3;
		eff_p2_a2 = e_p2_a2;
		eff_p2_a1 = e_p2_a1;
		eff_p2_a0 = e_p2_a0;
		eff_crossoverEnergy = e_crossoverEnergy;
		eff_overallProcentualError = e_overallProcentualError;
	}
	
	/**
	 * Set all calibrations at once. 
	 * @param e_a3 e_a3
	 * @param e_a2 e_a2
	 * @param e_a1 e_a1
	 * @param e_a0 e_a0
	 * @param f_a3 f_a3
	 * @param f_a2 f_a2
	 * @param f_a1 f_a1
	 * @param f_a0 f_a0
	 * @param f_overallProcentualError f_overallProcentualError
	 * @param e_p1_a4 e_p1_a4
	 * @param e_p1_a3 e_p1_a3
	 * @param e_p1_a2 e_p1_a2
	 * @param e_p1_a1 e_p1_a1
	 * @param e_p1_a0 e_p1_a0
	 * @param e_p2_a4 e_p2_a4
	 * @param e_p2_a3 e_p2_a3
	 * @param e_p2_a2 e_p2_a2
	 * @param e_p2_a1 e_p2_a1
	 * @param e_p2_a0 e_p2_a0
	 * @param e_crossoverEnergy e_crossoverEnergy
	 * @param e_overallProcentualError e_overallProcentualError
	 */
	public static void setCalibrations(double e_a3, double e_a2, double e_a1,
			double e_a0, //double c_a3, double c_a2, double c_a1, double c_a0,
			double f_a3, double f_a2, double f_a1, double f_a0,
			double f_overallProcentualError, double e_p1_a4, double e_p1_a3,
			double e_p1_a2, double e_p1_a1, double e_p1_a0, double e_p2_a4,
			double e_p2_a3, double e_p2_a2, double e_p2_a1, double e_p2_a0,
			double e_crossoverEnergy, double e_overallProcentualError) {
		en_a3 = e_a3;
		en_a2 = e_a2;
		en_a1 = e_a1;
		en_a0 = e_a0;
		//ch_a3 = c_a3;
		//ch_a2 = c_a2;
		//ch_a1 = c_a1;
		//ch_a0 = c_a0;
		fwhm_a3 = f_a3;
		fwhm_a2 = f_a2;
		fwhm_a1 = f_a1;
		fwhm_a0 = f_a0;
		fwhm_overallProcentualError = f_overallProcentualError;
		eff_p1_a4 = e_p1_a4;
		eff_p1_a3 = e_p1_a3;
		eff_p1_a2 = e_p1_a2;
		eff_p1_a1 = e_p1_a1;
		eff_p1_a0 = e_p1_a0;
		eff_p2_a4 = e_p2_a4;
		eff_p2_a3 = e_p2_a3;
		eff_p2_a2 = e_p2_a2;
		eff_p2_a1 = e_p2_a1;
		eff_p2_a0 = e_p2_a0;
		eff_crossoverEnergy = e_crossoverEnergy;
		eff_overallProcentualError = e_overallProcentualError;
	}

	/**
	 * Return FWHM at channel x (from calibration at ROI centroid)
	 * @param x x
	 * @return the result
	 */
	private double getFWHMFromChannel(double x) {
		double y = fwhm_a3 * Math.pow(x, 3) + fwhm_a2 * Math.pow(x, 2)
				+ fwhm_a1 * Math.pow(x, 1) + fwhm_a0;
		return y;
	}

	/**
	 * Return the energy in keV from calibration related to the channel x.
	 * @param x x
	 * @return the result
	 */
	private double getKevFromChannel(double x) {
		double y = en_a3 * Math.pow(x, 3) + en_a2 * Math.pow(x, 2) + en_a1
				* Math.pow(x, 1) + en_a0;
		return y;
	}

	/**
	 * Return the efficiency (%, number less than 100) computed from calibration at ROI centroid energy x (keV).
	 * @param x x
	 * @return the result
	 */
	public static double getEffFromEnergy(double x) {
		double y = 0.0;
		double e =0.0;
		if (x!=0.0){
			e=Math.log(x);
		}
		if (e<0.0)
			e=0.0;
		
		if (x < eff_crossoverEnergy) {
			y =  Math.exp(eff_p1_a4 * Math.pow(e, 4)
					+ eff_p1_a3 * Math.pow(e, 3) + eff_p1_a2
					* Math.pow(e, 2) + eff_p1_a1
					* Math.pow(e, 1) + eff_p1_a0);
		} else {
			y =  Math.exp(eff_p2_a4 * Math.pow(e, 4)
					+ eff_p2_a3 * Math.pow(e,3) + eff_p2_a2
					* Math.pow(e, 2) + eff_p2_a1
					* Math.pow(e, 1) + eff_p2_a0);
		}
		return y;
	}
	// -------------------------------------------------------------
	/**
	 * 
	 * @return the start channel
	 */
	public double getStartChannel() {
		return this.startChannel;
	}

	/**
	 * Set the start channel, an integer treated as a double
	 * @param startChannel startChannel
	 */
	public void setStartChannel(double startChannel) {
		this.startChannel = startChannel;
	}

	// ---------------
	/**
	 * Set the start energy in keV corresponding to the start channel
	 * @param d d
	 */
	public void setStartEnergy(double d) {
		this.startEnergy = d;
	}

	/**
	 * 
	 * @return the start energy in keV
	 */
	public double getStartEnergy() {
		return startEnergy;
	}

	// ---------------
	/**
	 * Return the start channel at edge
	 * @return the result
	 */
	public double getStartEdgeChannel() {
		return this.startEdgeChannel;
	}

	/**
	 * Set the start channel at edge
	 * @param d d
	 */
	public void setStartEdgeChannel(double d) {
		this.startEdgeChannel = d;
	}

	// ---------------
	/**
	 * Set the start energy at edge in keV
	 * @param d d
	 */
	public void setStartEdgeEnergy(double d) {
		this.startEdgeEnergy = d;
	}

	/**
	 * Return the start energy at edge in keV
	 * @return the result
	 */
	public double getStartEdgeEnergy() {
		return startEdgeEnergy;
	}

	// ---------------
	/**
	 * 
	 * @return the center channel, the median
	 */
	public double getCenterChannel() {
		return this.centerChannel;
	}

	/**
	 * Set the center channel (median channel)
	 */
	private void setCenterChannel() {
		this.centerChannel = (this.startChannel + this.endChannel) / 2.0;
	}

	// ---------------
	/**
	 * 
	 * @return the energy in keV associated to the central channel
	 */
	public double getCenterEnergy() {
		return this.centerEnergy;
	}

	/**
	 * Set the energy in keV associated to the central channel
	 * @param d d
	 */
	public void setCenterEnergy(double d) {
		this.centerEnergy = d;
	}

	// ---------------
	/**
	 * Set the centroid channel (very important quantity, which is the center channel weighted by pulses)
	 * @param d d
	 */
	public void setCentroidChannel(double d) {
		centroidChannel = d;
	}

	/**
	 * 
	 * @return the centroid channel
	 */
	public double getCentroidChannel() {
		return centroidChannel;
	}

	// ---------------
	/**
	 * 
	 * @return the energy in keV associated to the centroid channel
	 */
	public double getCentroidEnergy() {
		return this.centroidEnergy;
	}

	/**
	 * Set the energy in keV associated to the centroid channel
	 * @param d d
	 */
	public void setCentroidEnergy(double d) {
		this.centroidEnergy = d;
	}

	// ---------------
	/**
	 * Return the peak channel, which contain the maximum number of pulses in ROI.
	 * @return the result
	 */
	public double getPeakChannel() {
		return this.peakChannel;
	}

	/**
	 * Set the peak channel
	 * @param channel channel
	 */
	public void setPeakChannel(double channel) {
		this.peakChannel = channel;
	}

	// ---------------
	/**
	 * 
	 * @return the energy in keV associated to the peak channel
	 */
	public double getPeakEnergy() {
		return this.peakEnergy;
	}

	/**
	 * Set the energy in keV associated to the peak channel
	 * @param d d
	 */
	public void setPeakEnergy(double d) {
		this.peakEnergy = d;
	}

	/**
	 * 
	 * @return the pulses stored in peak channel
	 */
	public double getPeakPulses() {
		return this.peakPulses;
	}

	/**
	 * Set pulses associated to the peak channel
	 * @param pulses pulses
	 */
	public void setPeakPulses(double pulses) {
		this.peakPulses = pulses;
	}

	// ---------------
	/**
	 * 
	 * @return the end channel
	 */
	public double getEndChannel() {
		return this.endChannel;
	}

	/**
	 * Set the end channel
	 * @param endChannel endChannel
	 */
	public void setEndChannel(double endChannel) {
		this.endChannel = endChannel;
	}

	// .................
	/**
	 * Set the energy in keV associated to the end channel
	 * @param d d
	 */
	public void setEndEnergy(double d) {
		this.endEnergy = d;
	}

	/**
	 * 
	 * @return the energy in keV associated to the end channel
	 */
	public double getEndEnergy() {
		return endEnergy;
	}

	// ---------------
	/**
	 * 
	 * @return the end channel at edge
	 */
	public double getEndEdgeChannel() {
		return this.endEdgeChannel;
	}

	/**
	 * Set the end channel at edge
	 * @param d d
	 */
	public void setEndEdgeChannel(double d) {
		this.endEdgeChannel = d;
	}

	// ---------------
	/**
	 * Set the energy in keV associated to the end channel at edge
	 * @param d d
	 */
	public void setEndEdgeEnergy(double d) {
		this.endEdgeEnergy = d;
	}

	/**
	 * 
	 * @return the energy in keV associated to the end channel at edge
	 */
	public double getEndEdgeEnergy() {
		return endEdgeEnergy;
	}

	// ///////////////////////////////
	/**
	 * 
	 * @return the FWHM in channels associated to this ROI
	 */
	public double getFWHMChannel() {
		return fwhmChannel;
	}

	/**
	 * Set the FWHM in channels associated to this ROI
	 * @param d d
	 */
	public void setFWHMChannel(double d) {
		this.fwhmChannel = d;
	}

	// --------------
	/**
	 * 
	 * @return the FWHM error in channels
	 */
	public double getFWHMChannelError() {
		return fwhmChannelError;
	}

	/**
	 * Set the FWHM error in channels
	 * @param d d
	 */
	public void setFWHMChannelError(double d) {
		this.fwhmChannelError = d;
	}

	// --------------
	/**
	 * Set FWHM in keV associated to this ROI
	 * @param fwhmEnergy fwhmEnergy
	 */
	public void setFwhmEnergy(double fwhmEnergy) {
		this.fwhmEnergy = fwhmEnergy;
	}

	/**
	 * Return the FWHM in keV associated to this ROI
	 * @return the result
	 */
	public double getFwhmEnergy() {
		return fwhmEnergy;
	}

	/**
	 * Set the FWHM uncertainty in keV
	 * @param fwhmEnergyError fwhmEnergyError
	 */
	public void setFwhmEnergyError(double fwhmEnergyError) {
		this.fwhmEnergyError = fwhmEnergyError;
	}

	/**
	 * 
	 * @return the FWHM uncertainty in keV
	 */
	public double getFwhmEnergyError() {
		return fwhmEnergyError;
	}

	/**
	 * Set the FWHM energy based on FWHM calibration
	 * @param fwhmEnergyCalib fwhmEnergyCalib
	 */
	public void setFwhmEnergyCalib(double fwhmEnergyCalib) {
		this.fwhmEnergyCalib = fwhmEnergyCalib;
	}

	/**
	 * 
	 * @return the FWHM energy based on FWHM calibration
	 */
	public double getFwhmEnergyCalib() {
		return fwhmEnergyCalib;
	}

	/**
	 * Set the FWHM energy uncertainty based on FWHM calibration
	 * @param fwhmEnergyCalibError fwhmEnergyCalibError
	 */
	public void setFwhmEnergyCalibError(double fwhmEnergyCalibError) {
		this.fwhmEnergyCalibError = fwhmEnergyCalibError;
	}

	/**
	 * 
	 * @return the FWHM energy uncertainty based on FWHM calibration
	 */
	public double getFwhmEnergyCalibError() {
		return fwhmEnergyCalibError;
	}

	/**
	 * Set the ROI (peak) resolution being proportional to FWHM/ROIcentroid  
	 * @param resolution resolution
	 */
	public void setResolution(double resolution) {
		this.resolution = resolution;
	}

	/**
	 * 
	 * @return the resolution
	 */
	public double getResolution() {
		return resolution;
	}

	/**
	 * Set the ROI resolution based on FWHM calibration
	 * @param resolutionCalib resolutionCalib
	 */
	public void setResolutionCalib(double resolutionCalib) {
		this.resolutionCalib = resolutionCalib;
	}

	/**
	 * 
	 * @return the resolution computed using the FWHM calibration
	 */
	public double getResolutionCalib() {
		return resolutionCalib;
	}

	/**
	 * Significance is Yes or No and is related to the question: is FWHM computed from ROI statistically different than 
	 * the FWHM computed from calibration? A positive answer usually means there are multiple energy lines 
	 * in ROI, i.e. we have multiplets and not a clean ROI.    
	 * @param significance significance
	 */
	public void setSignificance(String significance) {
		this.significance = significance;
	}

	/**
	 * 
	 * @return the significance, i.e. yes or no to the statistic test of difference between FWHM computed and FWHM 
	 * taken from calibration
	 */
	public String getSignificance() {
		return significance;
	}

	/**
	 * Set ambient background counts (gross area)
	 * @param bkgCounts bkgCounts
	 */
	public void setBkgCounts(double bkgCounts) {
		this.bkgCounts = bkgCounts;
	}

	/**
	 * 
	 * @return the ambient background counts
	 */
	public double getBkgCounts() {
		return bkgCounts;
	}

	/**
	 * Set ambient background counts uncertainty
	 * @param bkgCountsError bkgCountsError
	 */
	public void setBkgCountsError(double bkgCountsError) {
		this.bkgCountsError = bkgCountsError;
	}

	/**
	 * 
	 * @return ambient background counts uncertainty
	 */
	public double getBkgCountsError() {
		return bkgCountsError;
	}

	/**
	 * Set ambient background counts rate
	 * @param bkgCountsRate bkgCountsRate
	 */
	public void setBkgCountsRate(double bkgCountsRate) {
		this.bkgCountsRate = bkgCountsRate;
	}

	/**
	 * 
	 * @return ambient background counts rate
	 */
	public double getBkgCountsRate() {
		return bkgCountsRate;
	}

	/**
	 * Set ambient background counts rate uncertainty
	 * @param bkgCountsRateError bkgCountsRateError
	 */
	public void setBkgCountsRateError(double bkgCountsRateError) {
		this.bkgCountsRateError = bkgCountsRateError;
	}

	/**
	 * 
	 * @return ambient background counts rate uncertainty
	 */
	public double getBkgCountsRateError() {
		return bkgCountsRateError;
	}

	/**
	 * Set ROI gross counts
	 * @param grossCounts grossCounts
	 */
	public void setGrossCounts(double grossCounts) {
		this.grossCounts = grossCounts;
	}

	/**
	 * 
	 * @return ROI gross counts
	 */
	public double getGrossCounts() {
		return grossCounts;
	}

	/**
	 * Set ROI gross counts uncertainty
	 * @param grossCountsError grossCountsError
	 */
	public void setGrossCountsError(double grossCountsError) {
		this.grossCountsError = grossCountsError;
	}

	/**
	 * 
	 * @return ROI gross counts uncertainty
	 */
	public double getGrossCountsError() {
		return grossCountsError;
	}

	/**
	 * Set ROI gross counts rate
	 * @param grossCountsRate grossCountsRate
	 */
	public void setGrossCountsRate(double grossCountsRate) {
		this.grossCountsRate = grossCountsRate;
	}

	/**
	 * 
	 * @return ROI gross counts rate
	 */
	public double getGrossCountsRate() {
		return grossCountsRate;
	}

	/**
	 * Set ROI gross counts rate uncertainty
	 * @param grossCountsRateError grossCountsRateError
	 */
	public void setGrossCountsRateError(double grossCountsRateError) {
		this.grossCountsRateError = grossCountsRateError;
	}

	/**
	 * 
	 * @return ROI gross counts rate uncertainty
	 */
	public double getGrossCountsRateError() {
		return grossCountsRateError;
	}

	/**
	 * Set continuum (Compton) background counts
	 * @param comptonCounts comptonCounts
	 */
	public void setComptonCounts(double comptonCounts) {
		this.comptonCounts = comptonCounts;
	}

	/**
	 * 
	 * @return continuum (Compton) background counts
	 */
	public double getComptonCounts() {
		return comptonCounts;
	}

	/**
	 * Set continuum (Compton) background counts uncertainty
	 * @param comptonCountsError comptonCountsError
	 */
	public void setComptonCountsError(double comptonCountsError) {
		this.comptonCountsError = comptonCountsError;
	}

	/**
	 * 
	 * @return continuum (Compton) background counts uncertainty
	 */
	public double getComptonCountsError() {
		return comptonCountsError;
	}

	/**
	 * Set continuum (Compton) background counts rate
	 * @param comptonCountsRate comptonCountsRate
	 */
	public void setComptonCountsRate(double comptonCountsRate) {
		this.comptonCountsRate = comptonCountsRate;
	}

	/**
	 * 
	 * @return continuum (Compton) background counts rate
	 */
	public double getComptonCountsRate() {
		return comptonCountsRate;
	}

	/**
	 * Set continuum (Compton) background counts rate uncertainty
	 * @param comptonCountsRateError comptonCountsRateError
	 */
	public void setComptonCountsRateError(double comptonCountsRateError) {
		this.comptonCountsRateError = comptonCountsRateError;
	}

	/**
	 * 
	 * @return continuum (Compton) background counts rate uncertainty
	 */
	public double getComptonCountsRateError() {
		return comptonCountsRateError;
	}

	/**
	 * Set ROI net counts (Net Area)
	 * @param netCounts netCounts
	 */
	public void setNetCounts(double netCounts) {
		this.netCounts = netCounts;
	}

	/**
	 * 
	 * @return ROI net counts
	 */
	public double getNetCounts() {
		return netCounts;
	}

	/**
	 * Set ROI net counts uncertainty
	 * @param netCountsError netCountsError
	 */
	public void setNetCountsError(double netCountsError) {
		this.netCountsError = netCountsError;
	}

	/**
	 * 
	 * @return ROI net counts uncertainty
	 */
	public double getNetCountsError() {
		return netCountsError;
	}

	/**
	 * Set ROI net counts rate
	 * @param netCountsRate netCountsRate
	 */
	public void setNetCountsRate(double netCountsRate) {
		this.netCountsRate = netCountsRate;
	}

	/**
	 * 
	 * @return ROI net counts rate
	 */
	public double getNetCountsRate() {
		return netCountsRate;
	}

	/**
	 * Set ROI net counts rate uncertainty
	 * @param netCountsRateError netCountsRateError
	 */
	public void setNetCountsRateError(double netCountsRateError) {
		this.netCountsRateError = netCountsRateError;
	}

	/**
	 * 
	 * @return ROI net counts rate uncertainty
	 */
	public double getNetCountsRateError() {
		return netCountsRateError;
	}

	/**
	 * Set confidence level for accuracy of statistical tests (usually 95%)
	 * @param confidenceLevel confidenceLevel
	 */
	public void setConfidenceLevel(double confidenceLevel) {
		this.confidenceLevel = confidenceLevel;
	}

	/**
	 * 
	 * @return confidence level
	 */
	public double getConfidenceLevel() {
		return confidenceLevel;
	}

	/**
	 * Assign ROI to a specific nuclide.
	 * @param nuclide nuclide
	 */
	public void setNuclide(String nuclide) {
		this.nuclide = nuclide;
	}

	/**
	 * 
	 * @return the nuclide
	 */
	public String getNuclide() {
		return nuclide;
	}

	/**
	 * Set probability of emission (yield) associated with the radiation energy
	 * @param yield yield
	 */
	public void setYield(double yield) {
		this.yield = yield;
	}

	/**
	 * 
	 * @return radiation yield
	 */
	public double getYield() {
		return yield;
	}

	/**
	 * Set detection efficiency associated to the ROI (%, number less than 100)
	 * @param efficiencyProcentual efficiencyProcentual
	 */
	public void setEfficiencyProcentual(double efficiencyProcentual) {
		this.efficiencyProcentual = efficiencyProcentual;
	}

	/**
	 * 
	 * @return efficiency associated to the ROI
	 */
	public double getEfficiencyProcentual() {
		return efficiencyProcentual;
	}

	/**
	 * Set efficiency uncertainty
	 * @param efficiencyProcentualError efficiencyProcentualError
	 */
	public void setEfficiencyProcentualError(double efficiencyProcentualError) {
		this.efficiencyProcentualError = efficiencyProcentualError;
	}

	/**
	 * 
	 * @return efficiency uncertainty
	 */
	public double getEfficiencyProcentualError() {
		return efficiencyProcentualError;
	}

	/**
	 * Set activity in Bq which gives this ROI. This is the purpose of gamma spectrommetry, compute sample 
	 * activity from ROI/ROIs.
	 * @param activity_Bq activity_Bq
	 */
	public void setActivity_Bq(double activity_Bq) {
		this.activity_Bq = activity_Bq;
	}

	/**
	 * 
	 * @return activity in Bq
	 */
	public double getActivity_Bq() {
		return activity_Bq;
	}

	/**
	 * Set activity uncertainty
	 * @param activity_BqError activity_BqError
	 */
	public void setActivity_BqError(double activity_BqError) {
		this.activity_BqError = activity_BqError;
	}

	/**
	 * 
	 * @return activity uncertainty
	 */
	public double getActivity_BqError() {
		return activity_BqError;
	}

	/**
	 * Set MDA (minimum detectable activity) in Bq associated with this ROI and the whole measurement process.
	 * @param mda_Bq mda_Bq
	 */
	public void setMda_Bq(double mda_Bq) {
		this.mda_Bq = mda_Bq;
	}

	/**
	 * 
	 * @return MDA (minimum detectable activity) in Bq
	 */
	public double getMda_Bq() {
		return mda_Bq;
	}

	/**
	 * Formally, set MDA uncertainty for statistic test. It is a formal quantity because this uncertainty is not 
	 * real in the sense that MDA has no "physical" error, but it make sense only for statistics.
	 * @param mda_BqError mda_BqError
	 */
	public void setMda_BqError(double mda_BqError) {
		this.mda_BqError = mda_BqError;
	}

	/**
	 * 
	 * @return MDA uncertainty
	 */
	public double getMda_BqError() {
		return mda_BqError;
	}

	/**
	 * Set difference to Yes or No as it answer the question: is the activity statistically different (greater) 
	 * than the MDA?
	 * @param difference difference
	 */
	public void setDifference(String difference) {
		this.difference = difference;
	}

	/**
	 * 
	 * @return Yes or No based on statistically difference between activity and MDA.
	 */
	public String getDifference() {
		return difference;
	}

	/**
	 * Set atomic mass for ROI associated nuclide.
	 * @param atomicMass atomicMass
	 */
	public void setAtomicMass(double atomicMass) {
		this.atomicMass = atomicMass;
	}

	/**
	 * 
	 * @return atomic mass for ROI associated nuclide.
	 */
	public double getAtomicMass() {
		return atomicMass;
	}

	/**
	 * Set half life for ROI associated nuclide.
	 * @param halfLife halfLife
	 */
	public void setHalfLife(double halfLife) {
		this.halfLife = halfLife;
	}

	/**
	 * 
	 * @return half life for ROI associated nuclide
	 */
	public double getHalfLife() {
		return halfLife;
	}

	/**
	 * Set half life measurement units for ROI associated nuclide.
	 * @param halfLifeUnits halfLifeUnits
	 */
	public void setHalfLifeUnits(String halfLifeUnits) {
		this.halfLifeUnits = halfLifeUnits;
	}

	/**
	 * 
	 * @return half life units
	 */
	public String getHalfLifeUnits() {
		return halfLifeUnits;
	}

	/**
	 * Set method for Net area calculation (NaI/Ge specific or Gaussian fit). Based on this selection 
	 * the FWHM is also computed differently (NaI specific, Ge specific or Gaussian fit). This is used for reports.
	 * @param netCalculationMethod netCalculationMethod
	 */
	public void setNetCalculationMethod(String netCalculationMethod) {
		this.netCalculationMethod = netCalculationMethod;
	}

	/**
	 * 
	 * @return method for Net area/FWHM calculation. This is used for reports.
	 */
	public String getNetCalculationMethod() {
		return netCalculationMethod;
	}

	/**
	 * Set method for MDA calculation (Pasternack, Curie, default). We strongly recommend the default method. 
	 * This is used for reports. 
	 * @param mdaCalculationMethod mdaCalculationMethod
	 */
	public void setMdaCalculationMethod(String mdaCalculationMethod) {
		this.mdaCalculationMethod = mdaCalculationMethod;
	}

	/**
	 * 
	 * @return the method for MDA calculation. This is used for reports.
	 */
	public String getMdaCalculationMethod() {
		return mdaCalculationMethod;
	}

	/**
	 * Set method for Net area calculation based on index. Used internally.
	 * @param imode imode
	 */
	public void setNetCalculationMethod(int imode) {
		this.iNet = imode;
	}

	/**
	 * Set method for MDA calculation based on index. Used internally.
	 * @param imode imode
	 */
	public void setMdaCalculationMethod(int imode) {
		this.iMDA = imode;
	}

	/**
	 * 
	 * @return method for Net area calculation based on index. Used internally.
	 */
	public int getNetCalculationMethod_internal() {// for fit!!
		return iNet;
	}

	/**
	 * 
	 * @return method for MDA calculation based on index. Used internally.
	 */
	private int getMdaCalculationMethod_internal() {
		return iMDA;
	}

	private double fwhmMinCh = 0.0;
	private double fwhmMaxCh = 0.0;

	/**
	 * Default (NaI specific) FWHM calculation.
	 * @return FWHM in channels
	 */
	private double defaultFWHMcomputation() {
		double c1 = channel[0];// start channel
		double c2 = channel[channel.length - 1];// end channel
		for (int i = 0; i < channel.length; i++) {
			// looking for start and end points for fwhm
			if ((pulses[i] < 0.5 * getPeakPulses())
					&& (channel[i] < getPeakChannel())) {
				c1 = channel[i];
			}

			if ((pulses[i] > 0.5 * getPeakPulses())
					&& (channel[i] > getPeakChannel())) {
				c2 = channel[i];
			}

			if ((pulses[i] < 0.5 * getPeakPulses())
					&& (channel[i] > getPeakChannel())) {
				break;
			}
		}

		double fwhmCh = c2 - c1;// not null
		fwhmMinCh = c1;
		fwhmMaxCh = c2;
		return fwhmCh;
	}

	/**
	 * FWHM calculation specific for Ge/HpGe detectors.
	 * @return FWHM in channels
	 */
	private double geFWHMcomputation() {
		double c1 = channel[0];// start channel
		double nc1 = pulses[0];// start channel pulses corrected for Compton
		for (int i = 0; i < channel.length; i++) {
			if ((pulses[i] > 0) && (channel[i] > 0)) {
				if (channel[i] < getPeakChannel()) {
					nc1 = pulses[i];
					c1 = channel[i];
				} else {
					c1 = channel[0];
					nc1 = 1.0;// 1 puls!!
				}
				break;
			}
		}

		if (nc1 == 0)
			nc1 = 1.0;// force!! just in case...never happen!

		// always peak>1.0 so if channel[i]>==peak channel...force nc1=1.0
		// so always getPeakPulses()/nc1>1!!! but just in case:
		if (getPeakPulses() / nc1 <= 1) {
			c1 = channel[0];
			nc1 = 1.0;// 1 puls!!
			// wrong peak set and also wrong results!!
		}

		// OK...see MAB_teor.pdf
		double fwhmCh = 2.355 * Math.sqrt(Math.abs((getPeakChannel() - c1)
				* (getPeakChannel() - c1)
				/ (2.0 * Math.log(getPeakPulses() / nc1))));

		fwhmMinCh = getPeakChannel() - fwhmCh / 2.0;
		fwhmMaxCh = getPeakChannel() + fwhmCh / 2.0; // solution

		return fwhmCh;
	}

	/**
	 * 
	 * @return the channel array
	 */
	public double[] getChannelData() {
		return channel;
	}

	/**
	 * 
	 * @return the corrected pulses array
	 */
	public double[] getComptonCorrectedPulses() {
		return pulses;
	}

	/**
	 * Some basic ROI computations (without nuclide assignment).
	 * Ambient BKG is channel by channel subtracted from sample spectrum.
	 */
	public void performBasicComputation() {
		int netMode = getNetCalculationMethod_internal();
		if (netMode == NET_CALCULATION_NAI) {
			setNetCalculationMethod("Default");
		} else if (netMode == NET_CALCULATION_GE) {
			setNetCalculationMethod("Ge_FWHM");
		} else if (netMode == NET_CALCULATION_GAUSSIAN) {
			setNetCalculationMethod("Gaussian_Fit");
		}
		// netMode is holding the Net area/FWHM calculation mode!

		// confidence level ...%->number!
		StatsUtil.confidenceLevel = getConfidenceLevel() / 100.0;

		// start energy
		double d = getStartChannel();
		d = getKevFromChannel(d);
		setStartEnergy(d);

		// end energy
		d = getEndChannel();
		d = getKevFromChannel(d);
		setEndEnergy(d);

		// center energy
		d = getCenterChannel();
		d = getKevFromChannel(d);
		setCenterEnergy(d);

		// start edge energy
		d = getStartEdgeChannel();
		d = getKevFromChannel(d);
		setStartEdgeEnergy(d);

		// end edge energy
		d = getEndEdgeChannel();
		d = getKevFromChannel(d);
		setEndEdgeEnergy(d);

		// spectrum live time
		double time = getLiveTime();// >0 always!

		// continuum Compton BKG COMPUTED BETWEEN
		// "REAL" START-END ROI CHANNELS!!
		double lowCh = getStartEdgeChannel();// edge
		double highCh = getEndEdgeChannel();// edge
		double lowPulses = getStartEdgePulses();
		double highPulses = getEndEdgePulses();
		double realLowCh = getStartChannel();// "real" ROI start channel
		double realHighCh = getEndChannel();// "real" end channel
		double realLowPulses, realHighPulses = 0.0;
		//@@AMBIENTAL BKG CORRECTION.......
		double alowPulses=abg.getAmbientalBKGPulsesAtChannel(lowCh);
		double ahighPulses=abg.getAmbientalBKGPulsesAtChannel(highCh);
		lowPulses=lowPulses-alowPulses;if(lowPulses<0.0) lowPulses=0.0;
		highPulses=highPulses-ahighPulses;if(highPulses<0.0) highPulses=0.0;
		//double lowPulsesErr=Math.sqrt(getStartEdgePulses()+alowPulses);
		//double highPulsesErr=Math.sqrt(getEndEdgePulses()+ahighPulses);
		//@@.........................
//System.out.println("low pulses= "+lowPulses+ " +- "+lowPulsesErr);
//System.out.println("high pulses= "+highPulses+ " +- "+highPulsesErr);
		// always highCh>lowCh
		realLowPulses = (lowPulses * highCh - highPulses * lowCh)
				/ (highCh - lowCh) + realLowCh * (highPulses - lowPulses)
				/ (highCh - lowCh);
		realHighPulses = (lowPulses * highCh - highPulses * lowCh)
				/ (highCh - lowCh) + realHighCh * (highPulses - lowPulses)
				/ (highCh - lowCh);
		//@@errors:
		//double realLowPulsesErr=Math.sqrt(
			//	lowPulsesErr*lowPulsesErr*(highCh-realLowCh)*(highCh-realLowCh)
				///((highCh - lowCh)*(highCh - lowCh))+
				//highPulsesErr*highPulsesErr*(realLowCh-lowCh)*(realLowCh-lowCh)
				///((highCh - lowCh)*(highCh - lowCh))
				//);
		//double realHighPulsesErr=Math.sqrt(
			//	lowPulsesErr*lowPulsesErr*(highCh-realHighCh)*(highCh-realHighCh)
				///((highCh - lowCh)*(highCh - lowCh))+
				//highPulsesErr*highPulsesErr*(realHighCh-lowCh)*(realHighCh-lowCh)
				///((highCh - lowCh)*(highCh - lowCh))
				//);
//System.out.println("real low pulses= "+realLowPulses+ " +- "+realLowPulsesErr);
//System.out.println("real high pulses= "+realHighPulses+ " +- "+realHighPulsesErr);

		//@@..........
		// we have linear (trapez) equations y1=mx1+n,....		
		// trapezoid area:		
		double comptonCnts = (realHighCh - realLowCh + 1.0)
				* (realLowPulses + realHighPulses) / 2.0;
		//@@modified:
		//double comptonCntsError = (realHighCh - realLowCh + 1.0)
			//	  * Math.sqrt(realLowPulsesErr*realLowPulsesErr + 
				//		  realHighPulsesErr*realHighPulsesErr) / 2.0;
		//the above may seem correct but is overestimated!!!!=>we select
		//start and end channel and build compton trapez!!=>but its unc
		//is physically related to trapez itself not to start and end 
		// pulses errors!
		double comptonCntsError=Math.sqrt(comptonCnts);
		//@@............
		// Compton counts computed between "real" start-end channels:
		// and with ambiental BKG correction!!!!!!!!!!!!!!!!!!!!!!!
		setComptonCounts(comptonCnts);
		setComptonCountsError(comptonCntsError);
//System.out.println("Compt= "+comptonCnts+" +- "+comptonCntsError);		
		setComptonCountsRate(comptonCnts / time);
		setComptonCountsRateError(comptonCntsError / time);

		// centroid, gross counts and ambiental bkg counts
		// also perform Compton subtraction and peak computation!
		double sumpc = 0.0;// sum of pulses x channel		
		double errp2 = 0.0;
		double centroid = 0.0;
		double gross = 0.0;// sum of pulses
		double bkgGross = 0.0;
		double p = 0.0;
		channel = new double[channelV.size()];
		pulses = new double[channelV.size()];
		double maxCh = 0.0;
		double maxP = 0.0;
		double maxPErr = 0.0;
		int ndat = channelV.size();// for Gaussian fit
		double[] ssig = new double[ndat];// for Gaussian fit
		double[] xx = new double[ndat];// for Gaussian fit
		double[] yy = new double[ndat];// for Gaussian fit
		double netArea0=0.0;
		double netArea0Err=0.0;
		for (int i = 0; i < channelV.size(); i++) {
			sumpc = sumpc + pulsesV.elementAt(i) * channelV.elementAt(i);
			gross = gross + pulsesV.elementAt(i);
			// bkg pulses are already scaled to sample spectrum.
			bkgGross = bkgGross + bkgpulsesV.elementAt(i);
			// construct channel data
			channel[i] = channelV.elementAt(i);
			
			// Compton pulses on channel i with ambiental BKG correction!!:
			// index based..0==x1; i==x; ...x1=x1, x2=x1+i(n); x=x1+i
			p = realLowPulses
					+
					// (i+realLowCh-realLowCh)*(realHighPulses-realLowPulses)/
					(i) * (realHighPulses - realLowPulses)
					/ (realHighCh - realLowCh);
			if (p<0.0)p=0.0;
			if (pulsesV.elementAt(i) - bkgpulsesV.elementAt(i)> p){
				pulses[i] = pulsesV.elementAt(i) - bkgpulsesV.elementAt(i) - p;
				// subtract ALL@@@@!
				// errp2 is variance of p (Compton pulses on channel i)
				//@@modified
				//errp2 = (1.0 - (i) / (realHighCh - realLowCh))
					//	* (1.0 - (i) / (realHighCh - realLowCh))
						//*realLowPulsesErr*realLowPulsesErr;
				//errp2 = errp2 + ((i) / (realHighCh - realLowCh))
					//	* ((i) / (realHighCh - realLowCh))
						//*realHighPulsesErr*realHighPulsesErr;
				errp2=p;//see above derivations!!!
				ssig[i] = Math.sqrt(pulsesV.elementAt(i) + bkgpulsesV.elementAt(i) + errp2);
				//@@
				//netArea0Err=netArea0Err+pulsesV.elementAt(i) + bkgpulsesV.elementAt(i) + errp2;
				
			} else{
				pulses[i] = 0.0;
				errp2=0.0;
				ssig[i]=1.0;//1 pulses in channel!!
				//netArea0Err=netArea0Err+0.0;
			}
			//netArea0=netArea0+pulses[i];//overestimated
			//NET AREA all correction performed!!
			// used in Gaussian fit routine:------------------
			xx[i] = channel[i];
			yy[i] = pulses[i];// ALL corrected pulses
			// data individual standard deviation:
			// -------------------------------------------------
			// peak corrected by ALL corrections=>valid peak!
			if (maxP < pulses[i]) {
				maxP = pulses[i];
				maxCh = channel[i];
				//@@added =ssig[i]
				//maxPErr=Math.sqrt(pulsesV.elementAt(i) + bkgpulsesV.elementAt(i) + errp2);
				maxPErr=ssig[i];
			}
		}
		netArea0=gross-bkgGross-comptonCnts;
		netArea0Err=Math.sqrt(gross+bkgGross+comptonCntsError*comptonCntsError);// sqrt!
		if (netArea0<0.0){
			netArea0=0.0;
			netArea0Err=0.0;
		}
		// Set-up peak data:
		setPeakPulses(maxP);
		setPeakChannel(maxCh);
		d = getKevFromChannel(maxCh);
		setPeakEnergy(d);

		// Set-up gross data and centroid:
		
		if (gross > 0) {
			centroid = sumpc / gross;
		} else {
			centroid = 0.0;
		}
		setCentroidChannel(centroid);
		d = getKevFromChannel(centroid);
		setCentroidEnergy(d);

		setGrossCounts(gross);
		setGrossCountsError(Math.sqrt(gross));
		setBkgCounts(bkgGross);
		setBkgCountsError(Math.sqrt(bkgGross));

		setGrossCountsRate(gross / time);
		setGrossCountsRateError(Math.sqrt(gross) / time);
		setBkgCountsRate(bkgGross / time);
		setBkgCountsRateError(Math.sqrt(bkgGross) / time);
//System.out.println("Compt= "+comptonCnts+" +- "+comptonCntsError);
//System.out.println("Gross= "+gross+" +- "+Math.sqrt(gross));
//System.out.println("AmbBKG= "+bkgGross+" +- "+Math.sqrt(bkgGross));
//System.out.println("Net= "+netArea0+" +- "+netArea0Err);
		// netArea and FWHM computation:
		double netArea = 0.0;
		double netAreaError = 0.0;
		if (iNet == NET_CALCULATION_NAI || iNet == NET_CALCULATION_GE) {
			// default net area calculation			
			netArea=netArea0;
			netAreaError=netArea0Err;

			if (iNet == NET_CALCULATION_NAI) {
				double fwhmCh = defaultFWHMcomputation();
				setFWHMChannel(fwhmCh);
				// setFWHMChannelError(fwhmCh*5.0/100.0);//default 5% error
				setFWHMChannelError(2.0);// 2 channels for c1, c2 evaluation!!!

				d = getKevFromChannel(fwhmMaxCh);
				d = d - getKevFromChannel(fwhmMinCh);
				setFwhmEnergy(d);

				d = getKevFromChannel(getPeakChannel());
				d = d - getKevFromChannel(getPeakChannel() - 2.0);
				// evaluate taking 2 channels arround the peak!
				setFwhmEnergyError(d);
			} else {
				double fwhmCh = geFWHMcomputation();
				setFWHMChannel(fwhmCh);
				// roi.errfwhmCh:= roi.fwhmCh*5.0/100.0;//default 5%
				setFWHMChannelError(1.0);// 1 channel reside from c1 guess!!

				d = getKevFromChannel(fwhmMaxCh);
				d = d - getKevFromChannel(fwhmMinCh);
				setFwhmEnergy(d);

				d = getKevFromChannel(getPeakChannel());
				d = d - getKevFromChannel(getPeakChannel() - 1.0);

				setFwhmEnergyError(d);
			}
		} else if (iNet == NET_CALCULATION_GAUSSIAN) {
			// Gaussian fit for both FWHM and Net Area calculation!
			ModelingData.func = this;// passes the required functions such as
										// fdf!!
			int mma = 3;// 3 coefficients!
			double[] acof = new double[mma];
			int[] iia = new int[mma];
			acof[0] = getPeakPulses();
			acof[1] = getPeakChannel();
			// initial guess:
			double fwhmCh = defaultFWHMcomputation();
			acof[2] = fwhmCh;// see c2-c1 which is not null!
			acof[2] = acof[2] / 2.355;// =>initial guess for sigma!
			// end initial guess
			iia[0] = 1;
			iia[1] = 1;
			iia[2] = 1;// ALL FIT!!
			double[][] ccovar = new double[mma][mma];
			double[][] alph = new double[mma][mma];

			ModelingData.mrqmin(xx, yy, ssig, ndat, acof, iia, mma, ccovar,
					alph, -1.0);
			if (ModelingData.failB)
				return;// fail!
			int iter = 1;
			for (int i = 1; i <= MAXITER_mrqmin; i++) {
				ModelingData.mrqmin(xx, yy, ssig, ndat, acof, iia, mma, ccovar,
						alph, ModelingData.alamda_mrqmin);
				if (ModelingData.failB)
					return;// fail!..never
				if (ModelingData.convB)
					break;// converge!
				iter++;
			}
			// due to the gamma spectrum nature, it can happened that mrgmin
			// will
			// not converge, but it is not a tragedy since the fit (which can be
			// viewed via view roi command) is OK in most cases !!!
			// of course, this is caused by error in pulses on all channels
			// SQRT(pulses) and user ROI set. Thats why in this case we will
			// compute netAreaError similar to grossError
			printSequence("gauss fit convergence = " + ModelingData.convB
					+ " after step no: " + iter + " chi2= "
					+ ModelingData.chisq_mrqmin);			
			ModelingData.mrqmin(xx, yy, ssig, ndat, acof, iia, mma, ccovar,
					alph, 0.0);
			if (ModelingData.failB)
				return;// fail!

			fwhmCh = 2.355 * acof[2];
			fwhmMinCh = getPeakChannel() - fwhmCh / 2.0;
			fwhmMaxCh = getPeakChannel() + fwhmCh / 2.0; // solution
			Nmax = acof[0];
			miu = acof[1];
			sgma = acof[2];
			netArea = Nmax * sgma * Math.sqrt(2.0 * Math.PI);// see teor!!
			// NET AREA COMPUTED HERE IS ALREADY ALL CORRECTED!!!!!
			if (ModelingData.convB)
				//@@modified.. in respect with peak pulses, 1 sigma error!!
				netAreaError = maxPErr* sgma//Math.sqrt(Nmax) * sgma
						* Math.sqrt(2.0 * Math.PI);
			else
				//@@modified.. in respect with all pulses affected by errors!!
				netAreaError = netArea0Err;//Math.sqrt(netArea);
			
			setFWHMChannel(fwhmCh);
			// roi.errfwhmCh:= roi.fwhmCh*5.0/100.0;//default 5%
			setFWHMChannelError(2.0);// from c1,c2 initial guess!

			d = getKevFromChannel(fwhmMaxCh);
			d = d - getKevFromChannel(fwhmMinCh);
			setFwhmEnergy(d);

			d = getKevFromChannel(getPeakChannel());//evaluate!
			d = d - getKevFromChannel(getPeakChannel() - 2.0);
			setFwhmEnergyError(d);			
		}

		double fwhmCalib = getFWHMFromChannel(getCentroidChannel());// in
																	// channels
		fwhmMinCh = getPeakChannel() - fwhmCalib / 2.0;
		fwhmMaxCh = getPeakChannel() + fwhmCalib / 2.0;
		d = getKevFromChannel(fwhmMaxCh);
		d = d - getKevFromChannel(fwhmMinCh);
		setFwhmEnergyCalib(d);
		setFwhmEnergyCalibError(d * fwhm_overallProcentualError / 100.0);

		if (getCentroidEnergy() != 0.0) {
			setResolution(100.0 * getFwhmEnergy() / getCentroidEnergy());
			setResolutionCalib(100.0 * getFwhmEnergyCalib()
					/ getCentroidEnergy());
		} else {
			setResolution(0.0);
			setResolutionCalib(0.0);
		}

		double fwhm_df = 0.0;// degrees of freedom!
		double fwhmCalib_df = StatsUtil.evaluateDegreesOfFreedom(
				getFwhmEnergyCalibError(), getFwhmEnergyCalib());
		if (StatsUtil.failB) fwhmCalib_df=10000.0;//stdev=0!!																	
		if (getFwhmEnergyError() != 0.0) {//No fail!!!
			fwhm_df = StatsUtil.evaluateDegreesOfFreedom(getFwhmEnergyError(),
					getFwhmEnergy());
		} else {
			fwhm_df = 10000;// infinity..means, accurate!
		}
		boolean diffB = StatsUtil.ttest_default_unc(getFwhmEnergy(),
				getFwhmEnergyCalib(), getFwhmEnergyError(),
				getFwhmEnergyCalibError(), fwhm_df, fwhmCalib_df);
		// System.out.println(diffB+"  "+StatsUtil.failB);//two
		// 0-es=>errors!!!but not real case!!
		if (diffB) {
			setSignificance("Yes");
		} else {
			setSignificance("No");
		}

		setNetCounts(netArea);
		setNetCountsError(netAreaError);

		setNetCountsRate(netArea / time);
		setNetCountsRateError(netAreaError / time);

		setRoiSet(true);//ROi is set not loaded from DB!!
		
		// double decayCorrection=0.5;System.out.println("n0 "+getNetCounts());
		// setNetCounts(getNetCounts()/decayCorrection);
		// setNetCountsError(getNetCountsError()/decayCorrection);
		// setNetCountsRate(getNetCountsRate()/decayCorrection);
		// setNetCountsRateError(getNetCountsRateError()/decayCorrection);
		// System.out.println("n "+getNetCounts());//WORKS

		// end..efficiency and other data are computed via updateRoi!!
	}

	/**
	 * Some basic operations (without nuclide assignment).
	 * It is based on net area subtraction from sample and ambient BKG, <br>
	 * which is not quite accurate!
	 * 
	 */
	@Deprecated
	public void performBasicComputation_old() {
		int netMode = getNetCalculationMethod_internal();
		if (netMode == NET_CALCULATION_NAI) {
			setNetCalculationMethod("Default");
		} else if (netMode == NET_CALCULATION_GE) {
			setNetCalculationMethod("Ge_FWHM");
		} else if (netMode == NET_CALCULATION_GAUSSIAN) {
			setNetCalculationMethod("Gaussian_Fit");
		}
		// netMode is holding the Net area/FWHM calculation mode!

		// confidence level ...%->number!
		StatsUtil.confidenceLevel = getConfidenceLevel() / 100.0;

		// start energy
		double d = getStartChannel();
		d = getKevFromChannel(d);
		setStartEnergy(d);

		// end energy
		d = getEndChannel();
		d = getKevFromChannel(d);
		setEndEnergy(d);

		// center energy
		d = getCenterChannel();
		d = getKevFromChannel(d);
		setCenterEnergy(d);

		// start edge energy
		d = getStartEdgeChannel();
		d = getKevFromChannel(d);
		setStartEdgeEnergy(d);

		// end edge energy
		d = getEndEdgeChannel();
		d = getKevFromChannel(d);
		setEndEdgeEnergy(d);

		// spectrum live time
		double time = getLiveTime();// >0 always!

		// continuum Compton BKG COMPUTED BETWEEN
		// "REAL" START-END ROI CHANNELS!!
		double lowCh = getStartEdgeChannel();// edge
		double highCh = getEndEdgeChannel();// edge
		double lowPulses = getStartEdgePulses();
		double highPulses = getEndEdgePulses();
		double realLowCh = getStartChannel();// "real" ROI start channel
		double realHighCh = getEndChannel();// "real" end channel
		double realLowPulses, realHighPulses = 0.0;
		// always highCh>lowCh
		realLowPulses = (lowPulses * highCh - highPulses * lowCh)
				/ (highCh - lowCh) + realLowCh * (highPulses - lowPulses)
				/ (highCh - lowCh);
		realHighPulses = (lowPulses * highCh - highPulses * lowCh)
				/ (highCh - lowCh) + realHighCh * (highPulses - lowPulses)
				/ (highCh - lowCh);
		// we have linear (trapez) equations y1=mx1+n,....
		// roundoff errors exist..but are small!!
		// trapezoid area:
		double comptonCnts = (realHighCh - realLowCh + 1.0)
				* (realLowPulses + realHighPulses) / 2.0;
		double comptonCntsError = (realHighCh - realLowCh + 1.0)
				* Math.sqrt(realLowPulses + realHighPulses) / 2.0;
		// Compton counts computed between "real" start-end channels:
		setComptonCounts(comptonCnts);
		setComptonCountsError(comptonCntsError);
		setComptonCountsRate(comptonCnts / time);
		setComptonCountsRateError(comptonCntsError / time);

		// centroid, gross counts and ambiental bkg counts
		// also perform Compton subtraction and peak computation!
		double sumpc = 0.0;// sum of pulses x channel
		double netBkg = 0.0;// ambiental net Area...if any
		double errnetBkg = 0.0;// its error
		double errp2 = 0.0;
		double centroid = 0.0;
		double gross = 0.0;// sum of pulses
		double bkgGross = 0.0;
		double p = 0.0;
		channel = new double[channelV.size()];
		pulses = new double[channelV.size()];
		double maxCh = 0.0;
		double maxP = 0.0;
		int ndat = channelV.size();// for Gaussian fit
		double[] ssig = new double[ndat];// for Gaussian fit
		double[] xx = new double[ndat];// for Gaussian fit
		double[] yy = new double[ndat];// for Gaussian fit
		for (int i = 0; i < channelV.size(); i++) {
			sumpc = sumpc + pulsesV.elementAt(i) * channelV.elementAt(i);
			gross = gross + pulsesV.elementAt(i);
			// bkg pulses are already scaled to sample spectrum.
			bkgGross = bkgGross + bkgpulsesV.elementAt(i);
			// construct channel data
			channel[i] = channelV.elementAt(i);
			// Compton pulses on channel i:
			// index based..0==x1; i==x; ...x1=x1, x2=x1+i(n); x=x1+i
			p = realLowPulses
					+
					// (i+realLowCh-realLowCh)*(realHighPulses-realLowPulses)/
					(i) * (realHighPulses - realLowPulses)
					/ (realHighCh - realLowCh);
			if (pulsesV.elementAt(i) > p)
				pulses[i] = pulsesV.elementAt(i) - p;// subtract continuum!
			else
				pulses[i] = 0.0;
			// used in Gaussian fit routine:------------------
			xx[i] = channel[i];
			yy[i] = pulses[i];// Compton corrected pulses
			// data individual standard deviation:
			// errp2 is variance of p (Compton pulses on channel i)
			errp2 = (1.0 - (i) / (realHighCh - realLowCh))
					* (1.0 - (i) / (realHighCh - realLowCh))
					* pulsesV.elementAt(0);
			errp2 = errp2 + ((i) / (realHighCh - realLowCh))
					* ((i) / (realHighCh - realLowCh))
					* pulsesV.elementAt(pulsesV.size() - 1);
			ssig[i] = Math.sqrt(pulsesV.elementAt(i) + errp2);// 1.0;
			// -------------------------------------------------
			// now the ambiental BKG computed between "real" start-end channels:
			p = bkgpulsesV.elementAt(0)
					+
					// (i+realLowCh-realLowCh)*(realHighPulses-realLowPulses)/
					(i)
					* (bkgpulsesV.elementAt(bkgpulsesV.size() - 1) - bkgpulsesV
							.elementAt(0)) / (realHighCh - realLowCh);

			if (bkgpulsesV.elementAt(i) > p) {
				// we have ambiental net area pulses=>store it!!
				netBkg = netBkg + bkgpulsesV.elementAt(i) - p;
				errp2 = (1.0 - (i) / (realHighCh - realLowCh))
						* (1.0 - (i) / (realHighCh - realLowCh))
						* bkgpulsesV.elementAt(0);
				errp2 = errp2 + ((i) / (realHighCh - realLowCh))
						* ((i) / (realHighCh - realLowCh))
						* bkgpulsesV.elementAt(bkgpulsesV.size() - 1);
				errnetBkg = errnetBkg + bkgpulsesV.elementAt(i) + errp2;
				// c=c1+c2=>errc2=errc12+errc22!! OK!
			}
			// peak corrected by Compton=>valid peak!
			if (maxP < pulses[i]) {
				maxP = pulses[i];
				maxCh = channel[i];
			}
		}
		// Adjust ambiental BKG error:
		errnetBkg = Math.sqrt(errnetBkg);// sqrt!

		// Set-up peak data:
		setPeakPulses(maxP);
		setPeakChannel(maxCh);
		d = getKevFromChannel(maxCh);
		setPeakEnergy(d);

		// Set-up gross data and centroid:
		// Note: ALWAYS GROSS COUNTS and BKG GROSS COUNTS ARE COMPUTED
		// for "real" START-END ROI channels!
		if (gross > 0) {
			centroid = sumpc / gross;
		} else {
			centroid = 0.0;
		}
		setCentroidChannel(centroid);
		d = getKevFromChannel(centroid);
		setCentroidEnergy(d);

		setGrossCounts(gross);
		setGrossCountsError(Math.sqrt(gross));
		setBkgCounts(bkgGross);
		setBkgCountsError(Math.sqrt(bkgGross));

		setGrossCountsRate(gross / time);
		setGrossCountsRateError(Math.sqrt(gross) / time);
		setBkgCountsRate(bkgGross / time);
		setBkgCountsRateError(Math.sqrt(bkgGross) / time);

		// netArea and FWHM computation:
		double netArea = 0.0;
		double netAreaError = 0.0;
		if (iNet == NET_CALCULATION_NAI || iNet == NET_CALCULATION_GE) {
			// default net area calculation
			netArea = gross - comptonCnts;// Compton Corection
			// errGross=SQRT(gross)=>errGross2=gross!
			netAreaError = Math.sqrt(gross + comptonCntsError
					* comptonCntsError);
			netArea = netArea - netBkg;// ambiental Bkg correction (if any)
			netAreaError = Math.sqrt(netAreaError * netAreaError + errnetBkg
					* errnetBkg);

			if (iNet == NET_CALCULATION_NAI) {
				double fwhmCh = defaultFWHMcomputation();
				setFWHMChannel(fwhmCh);
				// setFWHMChannelError(fwhmCh*5.0/100.0);//default 5% error
				setFWHMChannelError(2.0);// 2 channels for c1, c2 evaluation!!!

				d = getKevFromChannel(fwhmMaxCh);
				d = d - getKevFromChannel(fwhmMinCh);
				setFwhmEnergy(d);

				d = getKevFromChannel(getPeakChannel());
				d = d - getKevFromChannel(getPeakChannel() - 2.0);
				// evaluate taking 2 channels arround the peak!
				setFwhmEnergyError(d);
			} else {
				double fwhmCh = geFWHMcomputation();
				setFWHMChannel(fwhmCh);
				// roi.errfwhmCh:= roi.fwhmCh*5.0/100.0;//default 5%
				setFWHMChannelError(1.0);// 1 channel reside from c1 guess!!

				d = getKevFromChannel(fwhmMaxCh);
				d = d - getKevFromChannel(fwhmMinCh);
				setFwhmEnergy(d);

				d = getKevFromChannel(getPeakChannel());
				d = d - getKevFromChannel(getPeakChannel() - 1.0);

				setFwhmEnergyError(d);
			}
		} else if (iNet == NET_CALCULATION_GAUSSIAN) {
			// Gaussian fit for both FWHM and Net Area calculation!
			ModelingData.func = this;// passes the required functions such as
										// fdf!!
			int mma = 3;// 3 coefficients!
			double[] acof = new double[mma];
			int[] iia = new int[mma];
			acof[0] = getPeakPulses();
			acof[1] = getPeakChannel();
			// initial guess:
			double fwhmCh = defaultFWHMcomputation();
			acof[2] = fwhmCh;// see c2-c1 which is not null!
			acof[2] = acof[2] / 2.355;// =>initial guess for sigma!
			// end initial guess
			iia[0] = 1;
			iia[1] = 1;
			iia[2] = 1;// ALL FIT!!
			double[][] ccovar = new double[mma][mma];
			double[][] alph = new double[mma][mma];

			ModelingData.mrqmin(xx, yy, ssig, ndat, acof, iia, mma, ccovar,
					alph, -1.0);
			if (ModelingData.failB)
				return;// fail!
			int iter = 1;
			for (int i = 1; i <= MAXITER_mrqmin; i++) {
				ModelingData.mrqmin(xx, yy, ssig, ndat, acof, iia, mma, ccovar,
						alph, ModelingData.alamda_mrqmin);
				if (ModelingData.failB)
					return;// fail!..never
				if (ModelingData.convB)
					break;// converge!
				iter++;
			}
			// due to the gamma spectrum nature, it can happened that mrgmin
			// will
			// not converge, but it is not a tragedy since the fit (which can be
			// viewed via view roi command) is OK in most cases !!!
			// of course, this is caused by error in pulses on all channels
			// SQRT(pulses) and user ROI set. Thats why in this case we will
			// compute netAreaError similar to grossError
			printSequence("gauss fit convergence = " + ModelingData.convB
					+ " after step no: " + iter + " chi2= "
					+ ModelingData.chisq_mrqmin);
			// printSequence("Results are normal despite the appereance!!");
			ModelingData.mrqmin(xx, yy, ssig, ndat, acof, iia, mma, ccovar,
					alph, 0.0);
			if (ModelingData.failB)
				return;// fail!

			fwhmCh = 2.355 * acof[2];
			fwhmMinCh = getPeakChannel() - fwhmCh / 2.0;
			fwhmMaxCh = getPeakChannel() + fwhmCh / 2.0; // solution
			Nmax = acof[0];
			miu = acof[1];
			sgma = acof[2];
			netArea = Nmax * sgma * Math.sqrt(2.0 * Math.PI);// see teor!!
			// NET AREA COMPUTED HERE IS ALREADY COMPTON CORRECTED!!!!!
			if (ModelingData.convB)
				// in respect with peak pulses, 1 sigma error!!
				netAreaError = Math.sqrt(Nmax) * sgma
						* Math.sqrt(2.0 * Math.PI);
			else
				// in respect with all pulses affected by errors!!
				netAreaError = Math.sqrt(netArea);
			// Now we have CORRECTION only for Compton BKG!!
			// THE NEW AMBIENTAL BKG MUST BE SET...BETWEEN MIU-2SIGMA,
			// MIU+2SIGMA!!
			double stChannel = miu - 2.0 * sgma;// 95.44% coverage area!
			double endChannel = miu + 2.0 * sgma;// OK!
			// System.out.println("s , e "+stChannel+"  "+endChannel);
			// called the GammaAnalysisFrame where bkg is known!!!
			double[] ambBkgNetUnc = abg.getAmbientalNetAreaAndUnc(stChannel,
					endChannel);
			netBkg = ambBkgNetUnc[0];
			errnetBkg = ambBkgNetUnc[1];
			netArea = netArea - netBkg;// ambiental Bkg correction (if any)
			netAreaError = Math.sqrt(netAreaError * netAreaError + errnetBkg
					* errnetBkg);

			setFWHMChannel(fwhmCh);
			// roi.errfwhmCh:= roi.fwhmCh*5.0/100.0;//default 5%
			setFWHMChannelError(2.0);// from c1,c2 initial guess!

			d = getKevFromChannel(fwhmMaxCh);
			d = d - getKevFromChannel(fwhmMinCh);
			setFwhmEnergy(d);

			d = getKevFromChannel(getPeakChannel());//evaluate!
			d = d - getKevFromChannel(getPeakChannel() - 2.0);
			setFwhmEnergyError(d);			
		}

		double fwhmCalib = getFWHMFromChannel(getCentroidChannel());// in
																	// channels
		fwhmMinCh = getPeakChannel() - fwhmCalib / 2.0;
		fwhmMaxCh = getPeakChannel() + fwhmCalib / 2.0;
		d = getKevFromChannel(fwhmMaxCh);
		d = d - getKevFromChannel(fwhmMinCh);
		setFwhmEnergyCalib(d);
		setFwhmEnergyCalibError(d * fwhm_overallProcentualError / 100.0);

		if (getCentroidEnergy() != 0.0) {
			setResolution(100.0 * getFwhmEnergy() / getCentroidEnergy());
			setResolutionCalib(100.0 * getFwhmEnergyCalib()
					/ getCentroidEnergy());
		} else {
			setResolution(0.0);
			setResolutionCalib(0.0);
		}

		double fwhm_df = 0.0;// degrees of freedom!
		double fwhmCalib_df = StatsUtil.evaluateDegreesOfFreedom(
				getFwhmEnergyCalibError(), getFwhmEnergyCalib());// degrees of
																	// freedom!
		if (getFwhmEnergyError() != 0.0) {
			fwhm_df = StatsUtil.evaluateDegreesOfFreedom(getFwhmEnergyError(),
					getFwhmEnergy());
		} else {
			fwhm_df = 10000;// infinity..means, accurate!
		}
		boolean diffB = StatsUtil.ttest_default_unc(getFwhmEnergy(),
				getFwhmEnergyCalib(), getFwhmEnergyError(),
				getFwhmEnergyCalibError(), fwhm_df, fwhmCalib_df);
		// System.out.println(diffB+"  "+StatsUtil.failB);//two
		// 0-es=>errors!!!but not real case!!
		if (diffB) {
			setSignificance("Yes");
		} else {
			setSignificance("No");
		}

		setNetCounts(netArea);
		setNetCountsError(netAreaError);

		setNetCountsRate(netArea / time);
		setNetCountsRateError(netAreaError / time);

		setRoiSet(true);//ROi is set not loaded from DB!!
		
		// double decayCorrection=0.5;System.out.println("n0 "+getNetCounts());
		// setNetCounts(getNetCounts()/decayCorrection);
		// setNetCountsError(getNetCountsError()/decayCorrection);
		// setNetCountsRate(getNetCountsRate()/decayCorrection);
		// setNetCountsRateError(getNetCountsRateError()/decayCorrection);
		// System.out.println("n "+getNetCounts());//WORKS

		// end..efficiency and other data are computed via updateRoi!!
	}
	// =========IMPLEMENTATION OF FUNCTION INTERFACE===================
	/**
	 * Implementation of interface method. This is where to print text.
	 * @param s s
	 */
	public void printSequence(String s) {// for printing
		System.out.println(s);
	}

	// -------------------------------------------------------
	private double Nmax = 1.0;
	private double miu = 1.0;
	private double sgma = 1.0;

	/**
	 * The Gaussian fit function.
	 */
	public double F(double x) {// y=f(x)
		return Nmax * Math.exp(-(x - miu) * (x - miu) / (2.0 * sgma * sgma));
		// return 0.0;
	}

	// --------------------------------------------------------
	public double[] FD(double x) {// 0=->y=f(x) and 1-> dy=f'(x)
		return null;
	}

	public double MF(double[] x) {// y=f(x1,x2,...)
		return 0.0;
	}

	public double[] DMF(double[] x) {// y'1=df(x1,x2,...)/dx1;...the vector
										// gradient
		// df[1..n] evaluated at the input point x
		return null;
	}

	// ==============3D func====================
	public double F3D(double x, double y, double z) {
		return 0.0;
	}

	public double yy1(double x) {
		return 0.0;
	}

	public double yy2(double x) {
		return 0.0;
	}

	public double z1(double x, double y) {
		return 0.0;
	}

	public double z2(double x, double y) {
		return 0.0;
	}

	// end 3d====================================
	// non linear equation systems: root finding
	public double[] vecfunc(int n, double[] x) {
		return null;
	}

	// ==========================================
	public double[] aF(double x, int ma) {// poli fit
		// double[] y=new double[3];
		// y[0]=1.0;
		// y[1]=x;
		// y[2]=x*x;
		// return y;
		return null;
	}

	// fdf used in mrqmin-mrqcof routine!!!!!!!!!!!!!!!!!!
	/**
	 * Used in Gaussian fit (mrqmin-mrqcof routine). These are derivatives df/da.
	 */
	public double fdf(double x, double[] a, double[] dyda, int na) {
		// gauss fit etc. nonliniar
		double y = 0.0;
		// expt=Math.exp(-(x-a[1])*(x-a[1])/(2.0*a[2]*a[2]));
		double expt = Math.exp(-(x - a[1]) * (x - a[1]) / (2.0 * a[2] * a[2]));
		y = a[0] * expt;
		dyda[0] = expt;
		dyda[1] = y * ((x - a[1]) / (a[2] * a[2]));
		dyda[2] = y * ((x - a[1]) * (x - a[1]) / (a[2] * a[2] * a[2]));
		return y;
	}

	// ===========================================================
	public double[] derivF(double x, double[] y) {
		return null;
	}

	// ============2point
	public void load(double x1, double[] v, double[] y) {
		return;
	}

	public void load1(double x1, double[] v, double[] y) {
		return;
	}

	public void load2(double x2, double[] v, int nn2, double[] y) {
		return;
	}

	public void score(double x2, double[] y, double[] f) {
		return;
	}

	public void difeq(int k, int k1, int k2, int jsf, int is1, int isf,
			int indexv[], int ne, double[][] s, double[][] y) {
		return;
	}

	// ================================================
	public double g(double t) {// g(t)=FREDHOLM
		return 0.0;
	}

	public double ak(double t, double s) {// KERNEL
		return 0.0;
	}

	public double g(int k, double t) {// voltera
		return 0.0;
	}

	public double ak(int k, int l, double t, double s) {// voltera
		return 0.0;
	}

	public void kermom(double[] w, double y, int m) {
		return;
	}

	// =========END IMPLEMENTATION OF FUNCTION INTERFACE===================
	//**
	// * ROI nuclide assignment is performed here! <br>
	// * Requires good energy and efficiency calibration!
	// */
	/**
	 * Update ROI after nuclide assignment. The calling program must look in nuclide library and take all 
	 * energies and their yields which fall within ROI. Also take into account all corrections for accurate 
	 * ROI update.
	 * @param energies energies
	 * @param yields yields
	 * @param decayCorr decay correction
	 * @param energiesCorr energies used in coincidence correction
	 * @param yieldsCorr yields used in coincidence correction
	 * @param coinCorr coin correction
	 */
	public void updateRoi (double[] energies, double[] yields, double[] decayCorr, 
			double[] energiesCorr, double[] yieldsCorr,	double[] coinCorr) {
		performBasicComputation();//!!!!!!!!!!SET!!RESET!!!!
		int mdaMode = getMdaCalculationMethod_internal();
		if (mdaMode == MDA_CALCULATION_PASTERNACK) {
			setMdaCalculationMethod("Pasternack");
		} else if (mdaMode == MDA_CALCULATION_CURIE) {
			setMdaCalculationMethod("Curie");
		} else if (mdaMode == MDA_CALCULATION_DEFAULT) {
			setMdaCalculationMethod("Default");
		}
		// mdaMode is holding the MDA calculation mode!		
	
		// efficiency, yields and decayCorrection using weighted mean by yields!!
		Vector<Double> yldsV = new Vector<Double>();
		Vector<Double> effV = new Vector<Double>();
		Vector<Double> seffV = new Vector<Double>();
		Vector<Double> decayCorrV = new Vector<Double>();
		int lc = energies.length;
		double sumyields = 0.0;
		for (int i = 0; i < lc; i++) {
			if ((getStartEnergy() <= energies[i])
					&& (getEndEnergy() >= energies[i])) {
				// inside ROI!!
				sumyields = sumyields + yields[i];
				yldsV.addElement(yields[i]);
				double d = getEffFromEnergy(energies[i]);// in %!!
				effV.addElement(d);
				seffV.addElement(d * eff_overallProcentualError / 100.0);
				decayCorrV.addElement(decayCorr[i]);
			}
		}
		double[] weight = new double[effV.size()];
		double sumweight = 0.0;
		double sumweighteff = 0.0;
		double sumweightseff2 = 0.0;
		double sumweightdecayCorr = 0.0;
		for (int i = 0; i < effV.size(); i++) {
			if (sumyields != 0.0) {// always, for a good energy calibration!
				weight[i] = yldsV.elementAt(i) / sumyields;
			} else {
				weight[i] = 0.0;
			}

			sumweighteff = sumweighteff + weight[i] * effV.elementAt(i);
			sumweightseff2 = sumweightseff2 + weight[i] * weight[i]
					* seffV.elementAt(i) * seffV.elementAt(i);
			sumweight = sumweight + weight[i];
			sumweightdecayCorr = sumweightdecayCorr + weight[i]
					* decayCorrV.elementAt(i);
		}
		double decayCorrection = 1.0;
		if (sumweight != 0.0) {
			setEfficiencyProcentual(sumweighteff / sumweight);
			setEfficiencyProcentualError(Math.sqrt(sumweightseff2
					/ (sumweight * sumweight)));
			decayCorrection = sumweightdecayCorr / sumweight;
		} else {
			setEfficiencyProcentual(0.0);
			setEfficiencyProcentualError(0.0);
			decayCorrection = 1.0;
		}
		setYield(sumyields);

		setNetCounts(getNetCounts() / decayCorrection);
		setNetCountsError(getNetCountsError() / decayCorrection);
		setNetCountsRate(getNetCountsRate() / decayCorrection);
		setNetCountsRateError(getNetCountsRateError() / decayCorrection);

		// coincidence correction in the same manner as above:
		sumyields = 0.0;
		yldsV = new Vector<Double>();
		Vector<Double> coinCorrV = new Vector<Double>();
		int lcCorr = energiesCorr.length;
		for (int i = 0; i < lcCorr; i++) {
			if ((getStartEnergy() <= energiesCorr[i])
					&& (getEndEnergy() >= energiesCorr[i])) {
				// inside ROI!!
				sumyields = sumyields + yieldsCorr[i];
				yldsV.addElement(yieldsCorr[i]);

				coinCorrV.addElement(coinCorr[i]);
			}
		}
		weight = new double[yldsV.size()];
		sumweight = 0.0;
		double sumweightcoinCorr = 0.0;
		for (int i = 0; i < yldsV.size(); i++) {
			if (sumyields != 0.0) {// always, for a good energy calibration!
				weight[i] = yldsV.elementAt(i) / sumyields;
			} else {
				weight[i] = 0.0;
			}

			sumweight = sumweight + weight[i];
			sumweightcoinCorr = sumweightcoinCorr + weight[i]
					* coinCorrV.elementAt(i);
		}
		double coinCorrection = 1.0;
		if (sumweight != 0.0) {
			coinCorrection = sumweightcoinCorr / sumweight;
		} else {
			coinCorrection = 1.0;
		}
		//System.out.println("before "+getNetCounts());//WORKS
		setNetCounts(getNetCounts() / coinCorrection);
		setNetCountsError(getNetCountsError() / coinCorrection);
		setNetCountsRate(getNetCountsRate() / coinCorrection);
		setNetCountsRateError(getNetCountsRateError() / coinCorrection);
		//System.out.println("after "+getNetCounts());//WORKS
		double time = getLiveTime();// >0 always!
		double activity = 0.0;
		double activityError = 0.0;
		double net_df = 0.0;
		double eff_df = 0.0;
		double activity_df = 10000.0;
		double mda = 0.0;
		double mdaError = 0.0;
		double mda_df = 0.0;
		String diffS = "No";
		// activity and MDA:
		if (getEfficiencyProcentual() > 0.0 && (getYield() > 0.0)
				){//&& (getNetCountsRate() > 0.0)) {
			if (getNetCountsRate() > 0.0){
			activity = 100.0 * getNetCountsRate()
					/ (getYield() * getEfficiencyProcentual());
			activityError = activity
					* Math.sqrt(Math.pow(getNetCountsRateError()
							/ getNetCountsRate(), 2)
							+ Math.pow(getEfficiencyProcentualError()
									/ getEfficiencyProcentual(), 2));

			net_df = StatsUtil.evaluateDegreesOfFreedom(
					getNetCountsRateError(), getNetCountsRate());
			if (StatsUtil.failB) {
				net_df = 10000;// error=0=>infinit
			}
			eff_df = StatsUtil.evaluateDegreesOfFreedom(
					getEfficiencyProcentualError(), getEfficiencyProcentual());
			if (StatsUtil.failB) {
				eff_df = 10000;// error=0=>infinit
			}
			// construct terms for Welch-Satterthwaite formuls (linear
			// combination!)
			double abcompus0 = 0.0;
			if (activity != 0.0)
				abcompus0 = activityError / activity;
			double[] ab0 = new double[2];
			double[] f0 = new double[2];
			f0[0] = eff_df;
			f0[1] = net_df;
			if (getEfficiencyProcentual() != 0.0)
				ab0[0] = getEfficiencyProcentualError()
						/ getEfficiencyProcentual();
			else
				ab0[0] = 0.0;
			ab0[1] = getNetCountsRateError() / getNetCountsRate();

			activity_df = StatsUtil.getEffectiveDegreesOfFreedom(abcompus0,
					ab0, f0);

			// double tfactorA=StatsUtil.getStudentFactor(activity_df);
			// activityError=activityError*tfactorA;
			// Here we let 1 sigma error always!! in sample calculation
			// the final activities and uncertainties are set!!
			}
			// Detection limit:
			double ld = 0.0;
			double ldError = 0.0;
			double ld_df = 10000.0;
			//=========================================
			// From theory, comparison is made between SAMPLE and BKG
			// Here, we have 2 contribution to BKG: ambientBkg + continuumBkg
			// Therefore we have to add them.
			//=======================================
			double realBKG = getComptonCounts() + getBkgCounts();
			double realBKGerror = Math.sqrt(getComptonCountsError()*
					getComptonCountsError()+
					getBkgCountsError()*getBkgCountsError());
			// ==========================================
			if (mdaMode == MDA_CALCULATION_PASTERNACK) {
				ld = 3.29
						* (1.645 + Math
								.sqrt(2.706025 + 2.0 * realBKG))
						/ time;
				ldError = realBKGerror
						* 2.0
						* 1.645
						/ ((Math.sqrt(2.706025 + 2.0 * realBKG)) * time);
				ld_df = StatsUtil.evaluateDegreesOfFreedom(ldError, ld);
				if (StatsUtil.failB) ld_df=10000.0;//stdev=0!!
			} else if (mdaMode == MDA_CALCULATION_CURIE) {
				ld = (2.71 + 4.65 * Math.sqrt(realBKG)) / time;
				if(realBKG!=0.0)
					ldError = realBKGerror * 4.65
						* Math.sqrt(1.0 / (4.0 * realBKG)) / time;
				else
					ldError=0.0;
				ld_df = StatsUtil.evaluateDegreesOfFreedom(ldError, ld);
				if (StatsUtil.failB) ld_df=10000.0;//stdev=0!!
			} else if (mdaMode == MDA_CALCULATION_DEFAULT) {
				ld = (2.706025 + 3.29 * Math.sqrt(realBKG)) / time;
				if(realBKG!=0.0)
					ldError = realBKGerror * 1.645
						/ (Math.sqrt(realBKG) * time);
				else
					ldError=0.0;
				ld_df = StatsUtil.evaluateDegreesOfFreedom(ldError, ld);
				if (StatsUtil.failB) ld_df=10000.0;//stdev=0!!
			}

			mda = 100.0* ld / (getYield()*getEfficiencyProcentual());
			mdaError = mda
					* Math.sqrt(Math.pow(ldError / ld, 2)
							+ Math.pow(getEfficiencyProcentualError()
									/ getEfficiencyProcentual(), 2));

			// construct terms for Welch-Satterthwaite formuls (linear
			// combination!)
			double abcompus = 0.0;
			if (mda > 0.0)
				abcompus = mdaError / mda;
			double[] ab = new double[2];
			double[] f = new double[2];

			f[0] = eff_df;//
			f[1] = ld_df;
			ab[0] = getEfficiencyProcentualError() / getEfficiencyProcentual();
			ab[1] = ldError / ld;
			mda_df = StatsUtil.getEffectiveDegreesOfFreedom(abcompus, ab, f);

			// double tfactor=StatsUtil.getStudentFactor(mda_df);
			// mdaError=mdaError*tfactor;
			// Here we let 1 sigma error always!! in sample calculation
			// the final activities and uncertainties are set!!

			// Now the comparison!!!
			boolean diffB = StatsUtil.ttest_default_unc(activity, mda,
					activityError, mdaError, activity_df, mda_df);

			if (diffB && (activity > mda)) {
				diffS = "Yes";
			} else {
				diffS = "No";
			}

		}// if (getEfficiencyProcentual() > 0.0 ..............

		// now set activities:
		setActivity_Bq(activity);
		setActivity_BqError(activityError);

		setMda_Bq(mda);
		setMda_BqError(mdaError);

		setDifference(diffS);
	}

	/**
	 * Set true if this ROI is properly set in the sense that it is not loaded from DB.
	 * @param roiSet roiSet
	 */
	public void setRoiSet(boolean roiSet) {
		this.roiSet = roiSet;
	}

	/**
	 * 
	 * @return true if this ROI is properly set in the sense that it is not loaded from DB
	 */
	public boolean isRoiSet() {
		return roiSet;
	}

	// -------------------------end	
}
