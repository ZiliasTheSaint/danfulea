package danfulea.utils.gadgets;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.FontMetrics;
import java.awt.Graphics;
import java.awt.Polygon;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.Locale;

import javax.swing.JPanel;

/**
 * 
 * The panel containing the elegant clock. Based on tutorial/example freely available on the web.
 * 
 * @author Dan Fulea, 06 AUG. 2016
 * 
 */
public class ElegantClockPanel extends JPanel{

	private static final long serialVersionUID = 1L;
	
	private SimpleDateFormat formatter; // Formats the date displayed
	private Polygon lastsecondpolygon;//
	private Polygon lastminutepolygon;//
	private Polygon lasthourpolygon;//

	private int center_x;//
	private int center_y;//
	private Color hourcolor;//
	private Color minutecolor;//
	private Color secondcolor;//

	private int hourhand;//
	private int minutehand;//
	private int secondhand;//

	private int none = 0;
	private int slide = 1;
	private int on = 1;

	private Color fifteenminutecolor;//
	private Color fiveminutecolor;//
	private Color oneminutecolor;//

	private int fifteenminute;//
	private int fiveminute;//
	private int oneminute;//

	private int diamonds = 1;
	private int lines = 1;
	private int numbers = 2;

	private Color dialcolor;//
	private Color dialedgecolor;//

	private int dial;//
	private int withedge = 1;
	private int withoutedge = 2;
	
	/**
	 * Constructor and sets global members to default values.
	 */
	public ElegantClockPanel(){
		setDefaults();
	}
	
	/**
	 * Sets global members to default values
	 */
	private void setDefaults(){
		formatter = new SimpleDateFormat("MMMM dd yyyy '-' EEE '-'  HH:mm:ss",
				Locale.getDefault());
		secondcolor = Color.red;
		hourcolor = Color.black;
		minutecolor = Color.black;

		oneminutecolor = Color.darkGray;
		fiveminutecolor = new Color(184, 134, 11);
		fifteenminutecolor = new Color(0, 100, 0);

		secondhand = 1;
		// show minute hand: 0=Off;1=dupa sec;2=independent
		minutehand = 1;
		hourhand = 1;
		// oneminute to be lines if 1 or number if 2 or none if 0!
		oneminute = 1;
		// fiveminute to be diamonds if 1 or number if 2 or none if 0!
		fiveminute = 1;
		fifteenminute = 1;

		dial = 1;
		dialcolor = new Color(255, 250, 240);
		dialedgecolor = Color.darkGray;
	}

	/**
	 * Overrides default paintComponent method.
	 */
	public void paintComponent(Graphics g) {

		super.paintComponent(g);

		int n;

		Dimension dimension = getSize();

		center_x = dimension.width / 2;
		center_y = dimension.height / 2;

		lastsecondpolygon = new Polygon();
		lastminutepolygon = new Polygon();
		lasthourpolygon = new Polygon();

		if (dial == withedge) {
			g.setColor(dialedgecolor);
			g.fillOval(0, 0, dimension.width, dimension.height);

			g.setColor(dialcolor);
			g.fillOval((int) (dimension.width * .02),
					(int) (dimension.height * .02),
					(int) (dimension.width * .96),
					(int) (dimension.height * .96));
		} else if (dial == withoutedge) {
			g.setColor(dialcolor);
			g.fillOval(0, 0, dimension.width, dimension.height);
		}

		for (n = 0; n < 60; n++) {
			if (n % 15 == 0 && fifteenminute == diamonds) {
				g.setColor(fifteenminutecolor);
				g.fillPolygon(getDiamond(n, .4, 1, .8, .9));
			} else if (n % 15 == 0 && fifteenminute == numbers) {
				g.setColor(fifteenminutecolor);
				drawNumber(g, n, n == 0 ? 12 : n / 5, .17, .9);
			} else if (n % 5 == 0 && fiveminute == diamonds) {
				g.setColor(fiveminutecolor);
				g.fillPolygon(getDiamond(n, .2, .95, .85, .9));
			} else if (n % 5 == 0 && fiveminute == numbers) {
				g.setColor(fiveminutecolor);
				drawNumber(g, n, n == 0 ? 12 : n / 5, .1, .9);
			} else if (oneminute == lines) {
				//2PI = 60 seconds or minute=>1 sec(min) = 30/PI radians
				double sin = Math.sin(n * Math.PI / 30);
				double cos = Math.cos(n * Math.PI / 30);

				g.setColor(oneminutecolor);
				g.drawLine(center_x + (int) (center_x * .88 * sin),
						center_y - (int) (center_y * .88 * cos), center_x
								+ (int) (center_x * .92 * sin), center_y
								- (int) (center_y * .92 * cos));
			} else if (oneminute == numbers) {
				g.setColor(oneminutecolor);
				drawNumber(g, n, n, .05, .9);
			}
		}

		Polygon hourpolygon = new Polygon();
		Polygon minutepolygon = new Polygon();
		Polygon secondpolygon = new Polygon();

		double hours = 0;
		double minutes;
		double seconds;

		//get the time here----------------------
		Date date = new Date();

		formatter.applyPattern("s");
		try {
			seconds = Double.parseDouble(formatter.format(date));// 0-60
		} catch (NumberFormatException en) {
			seconds = 0;
		}
		formatter.applyPattern("m");
		try {
			minutes = Double.parseDouble(formatter.format(date));// 0-60
		} catch (NumberFormatException en) {
			minutes = 10;
		}
		formatter.applyPattern("h");
		try {
			hours = Double.parseDouble(formatter.format(date));// 1-12
		} catch (NumberFormatException en) {
			hours = 10;
		}
		//-----------------------------------------
		if (minutehand == slide)
			minutes = minutes + seconds / 60;
		if (hourhand == slide)
			hours = hours + minutes / 60;

		if (hourhand != none)
			hourpolygon = getHand(hours * 5, .4, .1, .05);
		if (minutehand != none)
			minutepolygon = getHand(minutes, .7, .2, .04);
		if (secondhand == on)
			secondpolygon = getHand(seconds, .8, .3, .03);

		if (dial != none)
			g.setColor(dialcolor);
		else
			g.setColor(getBackground());

		if (hourhand != none)
			drawHand(g, lasthourpolygon);
		if (minutehand != none)
			drawHand(g, lastminutepolygon);
		if (secondhand == on)
			drawHand(g, lastsecondpolygon);

		if (hourhand != none) {
			g.setColor(hourcolor);
			drawHand(g, hourpolygon);
		}

		if (minutehand != none) {
			g.setColor(minutecolor);
			drawHand(g, minutepolygon);
		}

		if (secondhand == on) {
			g.setColor(secondcolor);
			drawHand(g, secondpolygon);
		}

		lasthourpolygon = hourpolygon;
		lastminutepolygon = minutepolygon;
		lastsecondpolygon = secondpolygon;

		// Get the date to print at the bottom
		formatter.applyPattern("MMMM dd yyyy '-' EEE '-'  HH:mm:ss");
		String today = formatter.format(date);
		g.setColor(Color.BLACK);
		g.drawString(today, center_x - getSize().width / 4, center_y
				+ getSize().height / 4);
	}

	/**
	 * Diamond shape based on position around the clock (for major thicks - hours, minutes - as an option) used for internal computation.
	 * @param seconds seconds
	 * @param angle angle
	 * @param forwards forwards param
	 * @param backwards backward param
	 * @param sideways sideways param
	 * @return the diamond shape polygon object
	 */
	Polygon getDiamond(int seconds, double angle, double forwards,
			double backwards, double sideways) {
		Polygon polygon = new Polygon();

		double leftsin = Math.sin((seconds - angle) * Math.PI / 30);
		double leftcos = Math.cos((seconds - angle) * Math.PI / 30);
		double middlesin = Math.sin(seconds * Math.PI / 30);
		double middlecos = Math.cos(seconds * Math.PI / 30);
		double rightsin = Math.sin((seconds + angle) * Math.PI / 30);
		double rightcos = Math.cos((seconds + angle) * Math.PI / 30);

		polygon.addPoint(
				center_x + (int) (center_x * forwards * middlesin),
				center_y - (int) (center_y * forwards * middlecos));
		polygon.addPoint(center_x + (int) (center_x * sideways * leftsin),
				center_y - (int) (center_y * sideways * leftcos));
		polygon.addPoint(center_x
				+ (int) (center_x * backwards * middlesin), center_y
				- (int) (center_y * backwards * middlecos));
		polygon.addPoint(center_x + (int) (center_x * sideways * rightsin),
				center_y - (int) (center_y * sideways * rightcos));

		return polygon;
	}

	/**
	 * Draw number around the clock (as an option) used for internal computation.
	 * @param graphics the graphics
	 * @param seconds seconds
	 * @param number number
	 * @param size size
	 * @param radius radius
	 */
	void drawNumber(Graphics graphics, int seconds, int number,
			double size, double radius) {
		Font font = new Font("TimesRoman", Font.PLAIN, (int) (Math.min(
				center_x, center_y) * size));
		FontMetrics fontmetrics = getFontMetrics(font);

		String string = (new Integer(number)).toString();

		int x = center_x
				+ (int) (center_x * radius * Math.sin(seconds * Math.PI
						/ 30));
		int y = center_y
				- (int) (center_y * radius * Math.cos(seconds * Math.PI
						/ 30));

		x -= fontmetrics.stringWidth(string) / 2;
		y += fontmetrics.getAscent() / 2;

		graphics.setFont(font);
		graphics.drawString(string, x, y);
	}

	/**
	 * 
	 * @param seconds seconds
	 * @param forwards forwards
	 * @param backwards backwards
	 * @param sideways sideways
	 * @return the polygon shape depicted clock ticks (default second ticks - shapes around the clock where hands point to)
	 */
	Polygon getHand(double seconds, double forwards, double backwards,
			double sideways) {
		Polygon polygon = new Polygon();

		double sin = Math.sin(seconds * Math.PI / 30);
		double cos = Math.cos(seconds * Math.PI / 30);

		polygon.addPoint(center_x + (int) (center_x * forwards * sin),
				center_y - (int) (center_y * forwards * cos));
		polygon.addPoint(center_x + (int) (center_x * sideways * cos),
				center_y + (int) (center_y * sideways * sin));
		polygon.addPoint(center_x - (int) (center_x * backwards * sin),
				center_y + (int) (center_y * backwards * cos));
		polygon.addPoint(center_x - (int) (center_x * sideways * cos),
				center_y - (int) (center_y * sideways * sin));

		return polygon;
	}

	/**
	 * Draw clock hand
	 * @param graphics the graphics
	 * @param polygon the shape
	 */
	void drawHand(Graphics graphics, Polygon polygon) {
		// graphics.fillPolygon(polygon);
		graphics.drawLine(polygon.xpoints[0], polygon.ypoints[0],
				polygon.xpoints[2], polygon.ypoints[2]);
	}
}
