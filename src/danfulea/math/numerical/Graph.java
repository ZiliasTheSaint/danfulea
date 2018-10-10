package danfulea.math.numerical;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;

import javax.swing.BoxLayout;
import javax.swing.JFrame;
import javax.swing.JPanel;

import danfulea.utils.FrameUtilities;

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
 * Class for plotting user-supplied functions.
 * 
 * @author Dan Fulea, 15 OCT. 2006
 */
@SuppressWarnings("serial")
public class Graph extends JFrame {
	private static final Dimension PREFERRED_SIZE = new Dimension(800, 550);
	private ChartPanel polDirChartPanel;

	public static int interv = 200;
	public static double[] xarray;
	public static double[] yarray;

	private Function func;

	/**
	 * Constructor. Default point interval is 200.
	 * @param low lower x-value
	 * @param high higher x-value
	 * @param func the user function passed by a class implementing Function interface
	 */
	public Graph(double low, double high, Function func) {
		this.func = func;
		initPlotGraph(low, high);

		createGUI();

		// RefineryUtilities.centerFrameOnScreen(this);
		FrameUtilities.centerFrameOnScreen(this);
		setDefaultLookAndFeelDecorated(true);
		// createImageIcon(sa.resources.getString("form.icon.url"));

		setVisible(true);
		addWindowListener(new WindowAdapter() {
			public void windowClosing(WindowEvent e) {
				dispose();
				System.exit(0);// @@@@@@@@@@@@@@@@@
			}
		});
	}

	/**
	 * Sets up the window size.
	 */
	public Dimension getPreferredSize() {
		return PREFERRED_SIZE;
	}

	/**
	 * Evaluate the function and computes the array x, y to be plotted.
	 * @param low the lower x-value
	 * @param high the higher bound (x-value)
	 */
	public void initPlotGraph(double low, double high) {
		if (low >= high)
			return;
		double w = (high - low) / interv;
		xarray = new double[interv];
		yarray = new double[interv];
		for (int i = 0; i < interv; i++) {
			xarray[i] = low + i * w;// do not care if exceed UP
			double da = getY(xarray[i]);
			yarray[i] = da;
		}
	}

	/**
	 * Evaluate function at x
	 * @param x x
	 * @return the f(x)
	 */
	public double getY(double x) {
		return func.F(x);// 4.0-x*x;
	}

	/**
	 * GUI creation
	 */
	private void createGUI() {
		JPanel content = new JPanel(new BorderLayout());
		JPanel graphP = createMainGraphPanel();
		content.add(graphP, BorderLayout.CENTER);
		setContentPane(content);
		content.setOpaque(true); // content panes must be opaque
		pack();
	}

	/**
	 * 
	 * @return the chart
	 */
	private JFreeChart getPolyDirChartPanel() {
		XYSeries series = new XYSeries("series");

		for (int i = 0; i < interv; i++) {
			series.add(xarray[i], yarray[i]);
		}

		XYSeriesCollection data = new XYSeriesCollection(series);

		NumberAxis xAxis = new NumberAxis("x-axis");
		xAxis.setAutoRangeIncludesZero(false);
		NumberAxis yAxis = new NumberAxis("y-axis");
		XYItemRenderer renderer = new StandardXYItemRenderer(
				StandardXYItemRenderer.LINES);
		// ------------------------
		renderer.setSeriesPaint(0, Color.RED);
		// ---------------------------------
		// renderer.setToolTipGenerator(new StandardXYToolTipGenerator());
		renderer.setBaseToolTipGenerator(new StandardXYToolTipGenerator());
		XYPlot plot = new XYPlot(data, xAxis, yAxis, renderer);
		plot.setOrientation(PlotOrientation.VERTICAL);

		// XYSeries dotseries = new XYSeries("Lead Points");
		// for (int i=0;i<PEGS4A.energ.length;i++)
		// {
		// dotseries.add(PEGS4A.energ[i],PEGS4A.totg[i]);
		// }

		// XYSeriesCollection dotdata = new XYSeriesCollection(dotseries);
		// XYItemRenderer renderer2 = new
		// StandardXYItemRenderer(StandardXYItemRenderer.SHAPES);
		// renderer2.setSeriesPaint(0,Color.BLACK);
		// plot.setSecondaryDataset(0, dotdata);
		// plot.setSecondaryRenderer(0, renderer2);

		JFreeChart chart = new JFreeChart("y=f(x)",
				JFreeChart.DEFAULT_TITLE_FONT, plot, true);
		// chart.setBackgroundPaint(new GradientPaint(0, 0, Color.white, 1000,
		// 0, Color.YELLOW));//green));
		// chart.setBackgroundPaint(new GradientPaint(0, 0, Color.white, 1000,
		// 0, new Color(255,228,181)));
		return chart;
	}

	/**
	 * 
	 * @return the main GUI Panel
	 */
	private JPanel createMainGraphPanel() {
		polDirChartPanel = new ChartPanel(getPolyDirChartPanel());
		polDirChartPanel.setPreferredSize(new java.awt.Dimension(300, 250));
		polDirChartPanel.setMinimumSize(new java.awt.Dimension(300, 250));

		JPanel pd = new JPanel();
		BoxLayout bl = new BoxLayout(pd, BoxLayout.X_AXIS);
		pd.setLayout(bl);
		pd.add(polDirChartPanel);

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
		// mainP.setBackground(new Color(122,181,140));
		mainP.setBackground(Color.WHITE);

		return mainP;
	}
}
