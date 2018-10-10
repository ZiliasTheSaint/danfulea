package danfulea.phys.egs;
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
import javax.swing.JScrollPane;

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
import org.jfree.ui.RefineryUtilities;

/**
 * Plot cross sections computed in PEGS4A/PEGS4.
 * 
 * @author Dan Fulea, 22 AUG. 2005

 */
@SuppressWarnings("serial")
public class PEGS4AGraph extends JFrame{
	private static final Dimension PREFERRED_SIZE = new Dimension(800, 550);
	private ChartPanel polDirChartPanel;
	private ChartPanel elChartPanel;

	/**
	 * Constructor.
	 * @param low low energy for photon cross section
	 * @param high high energy for photon cross section
	 * @param lowe low energy for electron-positron cross section 
	 * @param highe high energy for electron-positron cross section
	 */
    public PEGS4AGraph(double low, double high, double lowe, double highe)
    {
		 PEGS4A.plotGamma(low,high);
		 PEGS4A.plotEP(lowe,highe);
		 createGUI();

		 RefineryUtilities.centerFrameOnScreen(this);
		 setDefaultLookAndFeelDecorated(true);
		 //createImageIcon(sa.resources.getString("form.icon.url"));

	     setVisible(true);
         addWindowListener(new WindowAdapter()
		 {
		     public void windowClosing(WindowEvent e)
		     {
		          dispose();
		          //System.exit(0);//@@@@@@@@@@@@@@@@@
		     }
         });
	}

    /**
	 * Sets up the window size.
	 */
    public Dimension getPreferredSize()
	{
	    return PREFERRED_SIZE;
    }


    /**
	 * GUI creation.
	 */
	private void createGUI()
	{
         JPanel content = new JPanel(new BorderLayout());
         JPanel graphP=createMainGraphPanel();
         content.add(graphP, BorderLayout.CENTER);
         setContentPane(new JScrollPane(content));
		 content.setOpaque(true); //content panes must be opaque
         pack();
	}

	/**
	 * Return the chart for photon cross sections
	 * @return the result
	 */
    private JFreeChart getPolyDirChartPanel()
    {
		XYSeries series = new XYSeries("photo");
		XYSeries series1 = new XYSeries("compton");
		XYSeries series2 = new XYSeries("pair");
		XYSeries series3 = new XYSeries("coherent");
		XYSeries series4 = new XYSeries("total");
		for (int i=0;i<PEGS4A.energ.length;i++)
		{
		   series.add(PEGS4A.energ[i],PEGS4A.photo[i]);
		   series1.add(PEGS4A.energ[i],PEGS4A.compt[i]);
		   series2.add(PEGS4A.energ[i],PEGS4A.pair[i]);
		   series3.add(PEGS4A.energ[i],PEGS4A.cohe[i]);
		   series4.add(PEGS4A.energ[i],PEGS4A.totg[i]);
        }

		XYSeriesCollection data = new XYSeriesCollection(series);
        data.addSeries(series1);
        data.addSeries(series2);
        data.addSeries(series3);
        data.addSeries(series4);

        NumberAxis xAxis = new NumberAxis("Energy [MeV]");
        xAxis.setAutoRangeIncludesZero(false);
        NumberAxis yAxis = new NumberAxis("Cross Sections");
        XYItemRenderer renderer = new StandardXYItemRenderer(StandardXYItemRenderer.LINES);
        //------------------------
        renderer.setSeriesPaint(0,Color.RED);
        renderer.setSeriesPaint(1,Color.BLUE);
        renderer.setSeriesPaint(2,Color.GREEN);
        renderer.setSeriesPaint(3,Color.ORANGE);
        renderer.setSeriesPaint(4,Color.BLACK);
        //---------------------------------
        //renderer.setToolTipGenerator(new StandardXYToolTipGenerator());
        renderer.setBaseToolTipGenerator(new StandardXYToolTipGenerator());
        XYPlot plot = new XYPlot(data, xAxis, yAxis, renderer);
        plot.setOrientation(PlotOrientation.VERTICAL);

        //XYSeries dotseries = new XYSeries("Lead Points");
        //for (int i=0;i<PEGS4A.energ.length;i++)
        //{
        //	dotseries.add(PEGS4A.energ[i],PEGS4A.totg[i]);
        //}

        //XYSeriesCollection dotdata = new XYSeriesCollection(dotseries);
        //XYItemRenderer renderer2 = new StandardXYItemRenderer(StandardXYItemRenderer.SHAPES);
        //renderer2.setSeriesPaint(0,Color.BLACK);
        //plot.setSecondaryDataset(0, dotdata);
        //plot.setSecondaryRenderer(0, renderer2);

        JFreeChart chart=new JFreeChart("Cross Sections (cm2/g)", JFreeChart.DEFAULT_TITLE_FONT, plot, true);
	    //chart.setBackgroundPaint(new GradientPaint(0, 0, Color.white, 1000, 0, Color.green));
        return chart;
 	}

    /**
	 * Return the chart for electron-positron cross sections
	 * @return the result
	 */
    private JFreeChart getElDirChartPanel()
    {
		XYSeries series = new XYSeries("brems");
		XYSeries series1 = new XYSeries("moller");
		XYSeries series2 = new XYSeries("bhabha");
		XYSeries series3 = new XYSeries("annih");
		XYSeries series4 = new XYSeries("total_E");
		XYSeries series5 = new XYSeries("total_P");
		for (int i=0;i<PEGS4A.energe.length;i++)
		{
		   series.add(PEGS4A.energe[i],PEGS4A.brem[i]);
		   series1.add(PEGS4A.energe[i],PEGS4A.moll[i]);
		   series2.add(PEGS4A.energe[i],PEGS4A.bhab[i]);
		   series3.add(PEGS4A.energe[i],PEGS4A.anni[i]);
		   series4.add(PEGS4A.energe[i],PEGS4A.tote[i]);
		   series5.add(PEGS4A.energe[i],PEGS4A.totp[i]);
        }

		XYSeriesCollection data = new XYSeriesCollection(series);
        data.addSeries(series1);
        data.addSeries(series2);
        data.addSeries(series3);
        data.addSeries(series4);
        data.addSeries(series5);

        NumberAxis xAxis = new NumberAxis("Energy [MeV]");
        xAxis.setAutoRangeIncludesZero(false);
        NumberAxis yAxis = new NumberAxis("Cross Sections");
        XYItemRenderer renderer = new StandardXYItemRenderer(StandardXYItemRenderer.LINES);
        //------------------------
        renderer.setSeriesPaint(0,Color.RED);
        renderer.setSeriesPaint(1,Color.BLUE);
        renderer.setSeriesPaint(2,Color.GREEN);
        renderer.setSeriesPaint(3,Color.ORANGE);
        renderer.setSeriesPaint(4,Color.BLACK);
        renderer.setSeriesPaint(5,Color.GRAY);
        //---------------------------------
        //renderer.setToolTipGenerator(new StandardXYToolTipGenerator());
        renderer.setBaseToolTipGenerator(new StandardXYToolTipGenerator());
        XYPlot plot = new XYPlot(data, xAxis, yAxis, renderer);
        plot.setOrientation(PlotOrientation.VERTICAL);

        //XYSeries dotseries = new XYSeries("Lead Points");
        //for (int i=0;i<PEGS4A.energe.length;i++)
        //{
        //	dotseries.add(PEGS4A.energe[i],PEGS4A.tote[i]);
        //}

        //XYSeriesCollection dotdata = new XYSeriesCollection(dotseries);
        //XYItemRenderer renderer2 = new StandardXYItemRenderer(StandardXYItemRenderer.SHAPES);
        //renderer2.setSeriesPaint(0,Color.BLACK);
        //plot.setSecondaryDataset(0, dotdata);
        //plot.setSecondaryRenderer(0, renderer2);

        JFreeChart chart=new JFreeChart("Cross Sections (cm2/g)", JFreeChart.DEFAULT_TITLE_FONT, plot, true);
	    //chart.setBackgroundPaint(new GradientPaint(0, 0, Color.white, 1000, 0, Color.green));
        return chart;
 	}

    /**
     * Create the plot panel.
     * @return the result
     */
	private JPanel createMainGraphPanel()
	{
		polDirChartPanel=new ChartPanel(getPolyDirChartPanel());
		polDirChartPanel.setPreferredSize(new java.awt.Dimension(300, 250));
		polDirChartPanel.setMinimumSize(new java.awt.Dimension(300, 250));
		elChartPanel=new ChartPanel(getElDirChartPanel());
		elChartPanel.setPreferredSize(new java.awt.Dimension(300, 250));
		elChartPanel.setMinimumSize(new java.awt.Dimension(300, 250));
		JPanel pd=new JPanel();
		BoxLayout bl = new BoxLayout(pd,BoxLayout.X_AXIS);
		pd.setLayout(bl);
		pd.add(polDirChartPanel);
		pd.add(elChartPanel);

        GridBagLayout gbl1 = new GridBagLayout();
        GridBagConstraints cts1 = new GridBagConstraints();
        cts1.anchor = GridBagConstraints.WEST;//default
        cts1.fill = GridBagConstraints.BOTH;
        Insets ins1 = new Insets(5,5,5,5);
        cts1.insets = ins1;

		JPanel mainP=new JPanel(gbl1);
		cts1.weightx = 1.0; //ponderi egale in maximizare toate pe orizontala
		cts1.gridheight = 1;
        cts1.weighty = 1.0;//maximizare pe verticala--atrage dupa sine tot restul!!
		cts1.gridwidth =GridBagConstraints.REMAINDER;
		gbl1.setConstraints(pd,cts1);
		mainP.add(pd);
		mainP.setBackground(new Color(122,181,140));

        return mainP;
	}

}
