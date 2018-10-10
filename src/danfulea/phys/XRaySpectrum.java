package danfulea.phys;

import java.awt.Color;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.util.ResourceBundle;

import javax.swing.BoxLayout;
import javax.swing.JPanel;

import danfulea.math.Convertor;
import danfulea.math.Sort;

import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.labels.StandardXYToolTipGenerator;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.StandardXYItemRenderer;
import org.jfree.chart.renderer.xy.XYItemRenderer;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

/**
 * Build X-Ray spectrum.
 * 
 * @author Dan Fulea, 11 APR. 2005
 */
public class XRaySpectrum {
	private static final String BASE_RESOURCE_CLASS = "danfulea.phys.resources.XRaySpectrum4Resources";//"danfulea.phys.resources.XRaySpectrum4Resources";
	private static ResourceBundle resources = ResourceBundle
			.getBundle(BASE_RESOURCE_CLASS);
	private double kv;// kilovoltage at console
	private double filtration;// X-Ray tube total equivalent filtration in mmAl
	private double anodAngle;// the anod angle
	// ----------------------------------------
	private double photonFlux = 0.0;// in photons/(mA*s*mm2)
	private double photonFlux_monoMamo = 0.0;// monoenergetic most common for
												// mamography!!
	private double airKerma = 0.0;// in microGy/(mA*s)
	private double airKerma_monoMamo = 0.0;// monoenergetic most common for
											// mamography!!
	private double normalizedValue = 0.0;
	private double[] xRayEnergies;// store energy in keV (x)
	private double[] xRayIntensities;// intensity in number of photons
	private double xRayEnergies_monoMamo = 0.0;// monoenergetic most common for
												// mamography!!
	private double xRayIntensities_monoMamo = 0.0;// monoenergetic most common
													// for mamography!!
	private int mode = 0;// handle which constructor is used!!
	// at corresponding energies normalized at normalizedValue (e.g. 100000!!)
	// (y)
	// private static final Dimension PREFERRED_SIZE = new Dimension(550, 550);
	private ChartPanel XRayChartPanel;

	// default spectrum monoenergetic for Mammography
	/**
	 * Constructor. Mono-energetic spectrum used in Mammography.
	 */
	public XRaySpectrum() {
		mode = 0;
		normalizedValue = ((Double) resources.getObject("normalizedSpectrum"))
				.doubleValue();
		photonFlux_monoMamo = ((Double) resources.getObject("mamo.common.flux"))
				.doubleValue();
		airKerma_monoMamo = ((Double) resources.getObject("mamo.common.kerma"))
				.doubleValue();
		xRayEnergies_monoMamo = ((Double) resources
				.getObject("mamo.common.energy")).doubleValue();
		xRayIntensities_monoMamo = ((Double) resources
				.getObject("mamo.common.intensity")).doubleValue();
	}

	/**
	 * Constructor. General X-Ray spectrum.
	 * @param kv tube voltage
	 * @param filtration tube filtration (total or equivalent in mmAl)
	 * @param anodAngle the anode angle. If not known, a good default value is 17 degrees.
	 */
	public XRaySpectrum(double kv, double filtration, double anodAngle) {
		mode = 1;
		this.kv = kv;
		this.filtration = filtration;
		this.anodAngle = anodAngle;
		generateSpectrum();
	}

	/**
	 * Get the number used for spectrum normalization (i.e. the intensities must be divided with this value). Used in Monte Carlo simulations. 
	 * If ICALC =1 this is given in photons/mAs/mm^2 at 75 cm distance from the tube.
	 * @return the result
	 */
	public double getNormalizedValue() {
		if (XRay.ICALC == 1)
			return XRay.TOT_PHOTONS;// idem with total photon flux at 75 cm!!
		return normalizedValue;
	}

	/**
	 * Return the photon flux for Monte Carlo simulations. If SRS78 is used (ICALC=1), this is given in photons/mAs/mm^2 at 75 cm distance from the tube.
	 * @return the result
	 */
	public double getPhotonFlux() {
		if (mode == 1) {
			if (XRay.ICALC == 1)
				return XRay.TOT_PHOTONS;
			return photonFlux;
		} else
			return photonFlux_monoMamo;
	}

	/**
	 * Return the air kerma for Monte Carlo simulations. If SRS78 is used (ICALC=1), this is given in uGy/mAs computed at 
	 * 75 cm distance from tube.
	 * @return the result
	 */
	public double getAirKerma() {
		if (mode == 1) {
			if (XRay.ICALC == 1)
				return XRay.KERMAPERMASAT750MM;
			return airKerma;
		} else
			return airKerma_monoMamo;
	}

	/**
	 * Return the array of X-Ray energies in keV.
	 * @return the result
	 */
	public double[] getXRayEnergies() {
		double[] result = new double[1];
		if (mode == 1) {
			if (XRay.ICALC == 1) {
				double[] xeney = new double[XRay.NEFF];
				for (int i = 1; i <= XRay.NEFF; i++) {
					xeney[i - 1] = XRay.EIN[i - 1];
				}
				return xeney;
			}
			return xRayEnergies;
		} else {
			result[0] = xRayEnergies_monoMamo;
			return result;
		}
	}

	/**
	 * Return the X-ray intensities. If SRS78 is used (ICALC=1), this is given in photons/mAs/mm^2 at 75 cm distance from the tube.
	 * @return the result
	 */
	public double[] getXRayIntensities() {
		double[] result = new double[1];
		if (mode == 1) {
			if (XRay.ICALC == 1) {
				double[] xeney = new double[XRay.NEFF];
				for (int i = 1; i <= XRay.NEFF; i++) {
					xeney[i - 1] = XRay.YF[i - 1];
				}
				return xeney;
			}

			return xRayIntensities;
		} else {

			result[0] = xRayIntensities_monoMamo;
			return result;
		}
	}

	/**
	 * Returns the graphic panel containing the spectrum plot.
	 * @return the result
	 */
	public JPanel getXRayPlotPanel() {
		JPanel jp;
		if (mode == 1)
			jp = createMainGraphPanel();
		else
			jp = new JPanel();// blank
		return jp;
	}

	// calculeaza m si n din cele 2 puncte p1(x1,y1) si p2(x2,y2)
	// si returneaza y(x).
	// pentru actualul mod nu se poate da eroare
	/**
	 * Linear interpolation.
	 * @param x1 first point x-value
	 * @param y1 first point y-value
	 * @param x2 second point x-value
	 * @param y2 second point y-value
	 * @param x desire point x-value
	 * @return desire point y-value
	 */
	private double linInt(double x1, double y1, double x2, double y2, double x) {
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

	// the most important method for this class
	/**
	 * Generates the unattenuated X-Ray spectrum. Called by constructor.
	 */
	private void generateSpectrum() {
		if (XRay.ICALC == 1) {
			// this.kv=kv;
			// this.filtration=filtration;
			// this.anodAngle=anodAngle;
			int ikv = (new Double(kv)).intValue();
			int ianod_file = (new Double(anodAngle)).intValue();
			String uanodS = "";
			if (anodAngle > 10) {
				if (XRay.ianod == 0)
					uanodS = ianod_file + "0";
				else if (XRay.ianod == 1)
					uanodS = ianod_file + "1";
				else if (XRay.ianod == 2)
					uanodS = ianod_file + "2";
			} else {
				if (XRay.ianod == 0)
					uanodS = "0" + ianod_file + "0";
				else if (XRay.ianod == 1)
					uanodS = "0" + ianod_file + "1";
				else if (XRay.ianod == 2)
					uanodS = "0" + ianod_file + "2";
			}
			String kvS = "";
			if (kv < 100)
				kvS = "0" + ikv;
			else
				kvS = "" + ikv;
			// XRay.ianod=0;===============DEFAULT
			// W==========================>EXTERNAL
			// XRay.iripple=0;=============DEFAULT 0,
			// CP===================>EXTERNAL
			// XRay.reset();==================>EXTERN
			String filename = kvS + uanodS;// "080170";
			XRay.ianod_file = ianod_file;

			XRay.KVP = kv;

			XRay.readSpectrum(filename);
			filename = "KERMAIR";// allways
			XRay.readKerma(filename);
			// filename="AL";//allways================>EXTERNAL???
			// XRay.readAttCoef(filename,1);//allways================>EXTERNAL???
			// XRay.TMM[0]=filtration;//allways================>EXTERNAL???
			// ======================================================
			XRay.buildSpectra();

			// in order to compute hermain air=>hvlcal
			XRay.computeHVL1("AL", false);

			return;
		}
		normalizedValue = ((Double) resources.getObject("normalizedSpectrum"))
				.doubleValue();
		double[] aArray = (double[]) resources.getObject("anodAngle");
		double[] fArray = (double[]) resources.getObject("filtration");
		double[] kArray = (double[]) resources.getObject("kv");
		// look for anodAngleIndexes
		Sort.findNearestValue(aArray, anodAngle, true);
		int indexAlow = Sort.getNearestPosition();
		int indexAhigh = indexAlow + 1;
		if (indexAlow == aArray.length - 1)// out of range
		{
			if (anodAngle < aArray[0]) {
				indexAlow = 0;
				indexAhigh = 1;
			} else// in the right outside!!
			{
				indexAlow = aArray.length - 2;
				indexAhigh = aArray.length - 1;
			}
		}
		double anodAngleMin = aArray[indexAlow];// -------------------------------------!!
		double anodAngleMax = aArray[indexAhigh];// -------------------------------------!!
		String anodAngleMins = Convertor.doubleToString(anodAngleMin);
		anodAngleMins = anodAngleMins.substring(0, anodAngleMins.length() - 2);
		String anodAngleMaxs = Convertor.doubleToString(anodAngleMax);
		anodAngleMaxs = anodAngleMaxs.substring(0, anodAngleMaxs.length() - 2);

		// look for kvIndexes
		Sort.findNearestValue(kArray, kv, true);
		int indexKlow = Sort.getNearestPosition();
		int indexKhigh = indexKlow + 1;
		if (indexKlow == kArray.length - 1)// out of range
		{
			if (kv < kArray[0]) {
				indexKlow = 0;
				indexKhigh = 1;
			} else// in the right outside!!
			{
				indexKlow = kArray.length - 2;
				indexKhigh = kArray.length - 1;
			}
		}
		double kvMin = kArray[indexKlow];// ----------------------------------------------!!
		double kvMax = kArray[indexKhigh];// ---------------------------------------------!!
		String kvMins = Convertor.doubleToString(kvMin);
		kvMins = kvMins.substring(0, kvMins.length() - 2);
		String kvMaxs = Convertor.doubleToString(kvMax);
		kvMaxs = kvMaxs.substring(0, kvMaxs.length() - 2);

		// look for filtrationIndexes
		Sort.findNearestValue(fArray, filtration, true);
		int indexFlow = Sort.getNearestPosition();
		int indexFhigh = indexFlow + 1;
		if (indexFlow == fArray.length - 1)// out of range
		{
			if (filtration < fArray[0]) {
				indexFlow = 0;
				indexFhigh = 1;
			} else// in the right outside!!
			{
				indexFlow = fArray.length - 2;
				indexFhigh = fArray.length - 1;
			}
		}
		double filtrationMin = fArray[indexFlow];// ------------------------------!!
		double filtrationMax = fArray[indexFhigh];// -----------------------------!!
		// --retrieving the resources reference
		String energyMinResources = "";
		String spectrumAminKminResources = "";// XRay-spectra
		String spectrumAmaxKminResources = "";// XRay-spectra
		if ((kvMins.compareTo("40") == 0) && (kvMaxs.compareTo("50") == 0)) {
			energyMinResources = "kv" + kvMins + ".energy.alter";
			spectrumAminKminResources = anodAngleMins + "kv" + kvMins
					+ ".spectrum.alter";// XRay-spectra
			spectrumAmaxKminResources = anodAngleMaxs + "kv" + kvMins
					+ ".spectrum.alter";// XRay-spectra
		} else {
			energyMinResources = "kv" + kvMins + ".energy";
			spectrumAminKminResources = anodAngleMins + "kv" + kvMins
					+ ".spectrum";// XRay-spectra
			spectrumAmaxKminResources = anodAngleMaxs + "kv" + kvMins
					+ ".spectrum";// XRay-spectra
		}
		String energyMaxResources = "kv" + kvMaxs + ".energy";
		String fixedAminKminResources = anodAngleMins + "kv" + kvMins
				+ ".fixed";// flux and kerma
		String fixedAminKmaxResources = anodAngleMins + "kv" + kvMaxs
				+ ".fixed";// flux and kerma
		String fixedAmaxKminResources = anodAngleMaxs + "kv" + kvMins
				+ ".fixed";// flux and kerma
		String fixedAmaxKmaxResources = anodAngleMaxs + "kv" + kvMaxs
				+ ".fixed";// flux and kerma
		String spectrumAminKmaxResources = anodAngleMins + "kv" + kvMaxs
				+ ".spectrum";// XRay-spectra
		String spectrumAmaxKmaxResources = anodAngleMaxs + "kv" + kvMaxs
				+ ".spectrum";// XRay-spectra
		// -----------INTERPOLATED SPECTRUM--the maximum energy
		// range!!-----------------------------
		double xriAminKminFmin = 0.0;
		double xriAminKminFmax = 0.0;
		double xriAminKminF = 0.0;// interp

		double xriAminKmaxFmin = 0.0;
		double xriAminKmaxFmax = 0.0;
		double xriAminKmaxF = 0.0;// interp

		double xriAmaxKminFmin = 0.0;
		double xriAmaxKminFmax = 0.0;
		double xriAmaxKminF = 0.0;// interp

		double xriAmaxKmaxFmin = 0.0;
		double xriAmaxKmaxFmax = 0.0;
		double xriAmaxKmaxF = 0.0;// interp

		double xriAminK = 0.0;
		double xriAmaxK = 0.0;

		double xriA = 0.0;

		xRayEnergies = (double[]) resources.getObject(energyMaxResources);
		xRayIntensities = new double[xRayEnergies.length];
		double[] xminRayEnergies = (double[]) resources
				.getObject(energyMinResources);
		double[] xmaxRayEnergies = (double[]) resources
				.getObject(energyMaxResources);
		double[][] xAminKminTable = (double[][]) resources
				.getObject(spectrumAminKminResources);
		double[][] xAmaxKminTable = (double[][]) resources
				.getObject(spectrumAmaxKminResources);
		double[][] xAminKmaxTable = (double[][]) resources
				.getObject(spectrumAminKmaxResources);
		double[][] xAmaxKmaxTable = (double[][]) resources
				.getObject(spectrumAmaxKmaxResources);

		for (int i = 0; i < xmaxRayEnergies.length; i++) {

			if (i < xminRayEnergies.length)// first index is energy reserved->+1
			{
				xriAminKminFmin = xAminKminTable[i][indexFlow + 1];
				xriAminKminFmax = xAminKminTable[i][indexFhigh + 1];
				xriAminKminF = linInt(filtrationMin, xriAminKminFmin,
						filtrationMax, xriAminKminFmax, filtration);
				xriAmaxKminFmin = xAmaxKminTable[i][indexFlow + 1];
				xriAmaxKminFmax = xAmaxKminTable[i][indexFhigh + 1];
				xriAmaxKminF = linInt(filtrationMin, xriAmaxKminFmin,
						filtrationMax, xriAmaxKminFmax, filtration);
			} else {
				xriAminKminFmin = 0.0;
				xriAminKminFmax = 0.0;
				xriAminKminF = 0.0;
				xriAmaxKminFmin = 0.0;
				xriAmaxKminFmax = 0.0;
				xriAmaxKminF = 0.0;
			}

			xriAminKmaxFmin = xAminKmaxTable[i][indexFlow + 1];
			xriAminKmaxFmax = xAminKmaxTable[i][indexFhigh + 1];
			xriAminKmaxF = linInt(filtrationMin, xriAminKmaxFmin,
					filtrationMax, xriAminKmaxFmax, filtration);
			xriAmaxKmaxFmin = xAmaxKmaxTable[i][indexFlow + 1];
			xriAmaxKmaxFmax = xAmaxKmaxTable[i][indexFhigh + 1];
			xriAmaxKmaxF = linInt(filtrationMin, xriAmaxKmaxFmin,
					filtrationMax, xriAmaxKmaxFmax, filtration);

			xriAminK = linInt(kvMin, xriAminKminF, kvMax, xriAminKmaxF, kv);
			xriAmaxK = linInt(kvMin, xriAmaxKminF, kvMax, xriAmaxKmaxF, kv);

			xriA = linInt(anodAngleMin, xriAminK, anodAngleMax, xriAmaxK,
					anodAngle);
			if (xmaxRayEnergies[i] <= kv)
				xRayIntensities[i] = xriA;
			else
				xRayIntensities[i] = 0.0;

			// System.out.println(" en[i]= "+xRayEnergies[i]+" min= "+xriAmaxKminF+
			// " max=" +xriAmaxKmaxF+" interpmin= "+xriAmaxK);//+
			// " interpmax= "+xriAmaxKmaxF);
		}
		// --------------------------------END SPECTRUM----------------------
		double[][] pfakAminKmin = (double[][]) resources
				.getObject(fixedAminKminResources);
		double[][] pfakAminKmax = (double[][]) resources
				.getObject(fixedAminKmaxResources);
		double[][] pfakAmaxKmin = (double[][]) resources
				.getObject(fixedAmaxKminResources);
		double[][] pfakAmaxKmax = (double[][]) resources
				.getObject(fixedAmaxKmaxResources);
		double pfAminKminFmin = pfakAminKmin[0][indexFlow];// photon flux
		double pfAminKminFmax = pfakAminKmin[0][indexFhigh];
		double pfAminKmaxFmin = pfakAminKmax[0][indexFlow];
		double pfAminKmaxFmax = pfakAminKmax[0][indexFhigh];
		double pfAmaxKminFmin = pfakAmaxKmin[0][indexFlow];
		double pfAmaxKminFmax = pfakAmaxKmin[0][indexFhigh];
		double pfAmaxKmaxFmin = pfakAmaxKmax[0][indexFlow];
		double pfAmaxKmaxFmax = pfakAmaxKmax[0][indexFhigh];
		double akAminKminFmin = pfakAminKmin[1][indexFlow];// airKerma
		double akAminKminFmax = pfakAminKmin[1][indexFhigh];
		double akAminKmaxFmin = pfakAminKmax[1][indexFlow];
		double akAminKmaxFmax = pfakAminKmax[1][indexFhigh];
		double akAmaxKminFmin = pfakAmaxKmin[1][indexFlow];
		double akAmaxKminFmax = pfakAmaxKmin[1][indexFhigh];
		double akAmaxKmaxFmin = pfakAmaxKmax[1][indexFlow];
		double akAmaxKmaxFmax = pfakAmaxKmax[1][indexFhigh];

		double pfAminKminF = linInt(filtrationMin, pfAminKminFmin,
				filtrationMax, pfAminKminFmax, filtration);
		double akAminKminF = linInt(filtrationMin, akAminKminFmin,
				filtrationMax, akAminKminFmax, filtration);
		double pfAmaxKminF = linInt(filtrationMin, pfAmaxKminFmin,
				filtrationMax, pfAmaxKminFmax, filtration);
		double akAmaxKminF = linInt(filtrationMin, akAmaxKminFmin,
				filtrationMax, akAmaxKminFmax, filtration);
		double pfAminKmaxF = linInt(filtrationMin, pfAminKmaxFmin,
				filtrationMax, pfAminKmaxFmax, filtration);
		double akAminKmaxF = linInt(filtrationMin, akAminKmaxFmin,
				filtrationMax, akAminKmaxFmax, filtration);
		double pfAmaxKmaxF = linInt(filtrationMin, pfAmaxKmaxFmin,
				filtrationMax, pfAmaxKmaxFmax, filtration);
		double akAmaxKmaxF = linInt(filtrationMin, akAmaxKmaxFmin,
				filtrationMax, akAmaxKmaxFmax, filtration);

		double pfAminK = linInt(kvMin, pfAminKminF, kvMax, pfAminKmaxF, kv);
		double akAminK = linInt(kvMin, akAminKminF, kvMax, akAminKmaxF, kv);
		double pfAmaxK = linInt(kvMin, pfAmaxKminF, kvMax, pfAmaxKmaxF, kv);
		double akAmaxK = linInt(kvMin, akAmaxKminF, kvMax, akAmaxKmaxF, kv);

		double pfA = linInt(anodAngleMin, pfAminK, anodAngleMax, pfAmaxK,
				anodAngle);
		double akA = linInt(anodAngleMin, akAminK, anodAngleMax, akAmaxK,
				anodAngle);

		photonFlux = pfA;
		airKerma = akA;
		// System.out.println("am: "+anodAngleMin+" amax: "+anodAngleMax);
		// System.out.println(" emin= "+energyMinResources+" emax= "+energyMaxResources);
		// System.out.println(" fxdmin= "+fixedMinResources+" fxdmax= "+fixedMaxResources);
		// System.out.println(" xmin= "+spectrumMinResources+" xmax= "+spectrumMaxResources);
	}

	// Auxiliary display

	/**
	 * Return the X-Ray chart.
	 * @return the result
	 */
	private JFreeChart getXRayChart() {
		XYSeries series = new XYSeries(resources.getString("graphics.XYSeries"));
		// XYSeries series1 = new
		// XYSeries(sa.resources.getString("graphics.XYSeries1"));
		XYSeries dotseries = new XYSeries(
				resources.getString("graphics.Points"));
		if (XRay.ICALC == 1) {
			for (int i = 0; i < XRay.NEFF; i++) {
				series.add(XRay.EIN[i], XRay.YF[i]);
				dotseries.add(XRay.EIN[i], XRay.YF[i]);
			}
		} else {
			for (int i = 0; i < xRayEnergies.length; i++) {
				series.add(xRayEnergies[i], xRayIntensities[i]);
				// series1.add(sa.loged[i],sa.lwg);
				dotseries.add(xRayEnergies[i], xRayIntensities[i]);
			}
		}

		XYSeriesCollection data = new XYSeriesCollection(series);
		// XYSeriesCollection dotdata = new XYSeriesCollection(dotseries);
		// data.addSeries(series1);

		NumberAxis xAxis = new NumberAxis(
				resources.getString("graphics.axes.1"));
		xAxis.setAutoRangeIncludesZero(false);
		NumberAxis yAxis = new NumberAxis(
				resources.getString("graphics.axes.2"));
		// XYPlot plot = new XYPlot(data, xAxis, yAxis, renderer);
		// plot.setOrientation(PlotOrientation.VERTICAL);

		// XYSeries dotseries = new
		// XYSeries(resources.getString("graphics.Points"));
		// for (int i=0;i<xRayEnergies.length;i++)
		// {
		// dotseries.add(xRayEnergies[i],xRayIntensities[i]);
		// }
		XYPlot plot = new XYPlot();
		plot.setOrientation(PlotOrientation.VERTICAL);
		plot.setBackgroundPaint(Color.lightGray);
		plot.setDomainGridlinePaint(Color.white);
		plot.setRangeGridlinePaint(Color.white);
		// allow chart movement by pressing CTRL and drag with mouse!
		plot.setDomainPannable(true);
		plot.setRangePannable(true);
		// 1st axis
		plot.setDomainAxis(0, xAxis);// the axis index;axis
		plot.setRangeAxis(0, yAxis);
		// DATASET AND RENDERER:
		int idataset = 0;
		plot.setDataset(idataset, data);// idataset=0!

		XYItemRenderer renderer = new StandardXYItemRenderer(
				StandardXYItemRenderer.LINES);
		renderer.setSeriesPaint(0, Color.RED);
		// renderer.setToolTipGenerator(new StandardXYToolTipGenerator());
		renderer.setSeriesPaint(0, Color.RED);
		renderer.setBaseToolTipGenerator(new StandardXYToolTipGenerator());
		plot.setRenderer(idataset, renderer);

		XYSeriesCollection dotdata = new XYSeriesCollection(dotseries);
		idataset = 1;// channel,pulses BKG
		plot.setDataset(idataset, dotdata);

		XYItemRenderer renderer2 = new StandardXYItemRenderer(
				StandardXYItemRenderer.SHAPES);
		renderer2.setSeriesPaint(0, Color.BLACK);
		plot.setRenderer(idataset, renderer2);
		// plot.setSecondaryDataset(0, dotdata);
		// plot.setSecondaryRenderer(0, renderer2);
		// @@@@@@@@@@@@@@@@@@@
		String sst = "";
		if (XRay.ICALC == 1)
			sst = resources.getString("graphics.Title2");
		else
			sst = resources.getString("graphics.Title");

		JFreeChart chart = new JFreeChart(sst, JFreeChart.DEFAULT_TITLE_FONT,
				plot, true);
		// chart.setBackgroundPaint(new GradientPaint(0, 0, Color.white, 1000,
		// 0, Color.green));
		return chart;
	}

	/**
	 * Create the graphic panel containing the X-Ray plot.
	 * @return the result
	 */
	private JPanel createMainGraphPanel() {
		XRayChartPanel = new ChartPanel(getXRayChart());
		XRayChartPanel.setPreferredSize(new java.awt.Dimension(600, 550));
		XRayChartPanel.setMinimumSize(new java.awt.Dimension(600, 550));
		JPanel pd = new JPanel();
		BoxLayout bl = new BoxLayout(pd, BoxLayout.X_AXIS);
		pd.setLayout(bl);
		pd.add(XRayChartPanel);

		GridBagLayout gbl1 = new GridBagLayout();
		GridBagConstraints cts1 = new GridBagConstraints();
		cts1.anchor = GridBagConstraints.WEST;// default
		cts1.fill = GridBagConstraints.BOTH;
		Insets ins1 = new Insets(5, 5, 5, 5);
		cts1.insets = ins1;

		JPanel mainP = new JPanel(gbl1);
		cts1.weightx = 1.0; // ponderi egale in maximizare toate pe orizontala
		cts1.gridheight = 1;
		cts1.weighty = 1.0;// maximizare pe verticala--atrage dupa sine tot
							// restul!!
		cts1.gridwidth = GridBagConstraints.REMAINDER;
		gbl1.setConstraints(pd, cts1);
		mainP.add(pd);
		mainP.setBackground(new Color(122, 181, 140));

		return mainP;
	}
}
