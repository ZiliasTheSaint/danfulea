package danfulea.utils.table;

import java.awt.Color;
import java.awt.Component;
import java.awt.Graphics;

import javax.swing.Icon;
import javax.swing.UIManager;

/**
 * The arrow icon.
 * Based on example code written by Joris Van den Bogaert see esus.com
 * @author Dan Fulea, 22 Mar. 2005
 */

public class BevelArrowIcon implements Icon {
	public static final int UP = 0; // direction
	public static final int DOWN = 1;

	private static final int DEFAULT_SIZE = 11;

	private Color edge1;
	private Color edge2;
	private Color fill;
	private int size;
	private int direction;

	/**
	 * Constructor
	 * @param direction direction
	 * @param isRaisedView isRaisedView
	 * @param isPressedView isPressedView
	 */
	public BevelArrowIcon(int direction, boolean isRaisedView,
			boolean isPressedView) {
		if (isRaisedView) {
			if (isPressedView) {
				init(UIManager.getColor("controlLtHighlight"),
						UIManager.getColor("controlDkShadow"),
						UIManager.getColor("controlShadow"), DEFAULT_SIZE,
						direction);
			} else {
				init(UIManager.getColor("controlHighlight"),
						UIManager.getColor("controlShadow"),
						UIManager.getColor("control"), DEFAULT_SIZE, direction);
			}
		} else {
			if (isPressedView) {
				init(UIManager.getColor("controlDkShadow"),
						UIManager.getColor("controlLtHighlight"),
						UIManager.getColor("controlShadow"), DEFAULT_SIZE,
						direction);
			} else {
				init(UIManager.getColor("controlShadow"),
						UIManager.getColor("controlHighlight"),
						UIManager.getColor("control"), DEFAULT_SIZE, direction);
			}
		}
	}

	/**
	 * Constructor
	 * @param edge1 edge1 color
	 * @param edge2 edge2 color
	 * @param fill fill color
	 * @param size size
	 * @param direction direction
	 */
	public BevelArrowIcon(Color edge1, Color edge2, Color fill, int size,
			int direction) {
		init(edge1, edge2, fill, size, direction);
	}

	/**
	 * Paint specific method.
	 */
	public void paintIcon(Component c, Graphics g, int x, int y) {
		switch (direction) {
		case DOWN:
			drawDownArrow(g, x, y);
			break;
		case UP:
			drawUpArrow(g, x, y);
			break;
		}
	}

	/**
	 * Return the icon width
	 * @return the result
	 */
	public int getIconWidth() {
		return size;
	}

	/**
	 * Return the icon height
	 * @return the result
	 */
	public int getIconHeight() {
		return size;
	}

	/**
	 * Internally used by 2nd constructor
	 * @param edge1 edge1 color
	 * @param edge2 edge2 color
	 * @param fill fill color
	 * @param size size
	 * @param direction direction
	 */
	private void init(Color edge1, Color edge2, Color fill, int size,
			int direction) {
		this.edge1 = edge1;
		this.edge2 = edge2;
		this.fill = fill;
		this.size = size;
		this.direction = direction;
	}

	/**
	 * Internally used by paintIcon.
	 * @param g g
	 * @param xo xo
	 * @param yo yo
	 */
	private void drawDownArrow(Graphics g, int xo, int yo) {
		g.setColor(edge1);
		g.drawLine(xo, yo, xo + size - 1, yo);
		g.drawLine(xo, yo + 1, xo + size - 3, yo + 1);
		g.setColor(edge2);
		g.drawLine(xo + size - 2, yo + 1, xo + size - 1, yo + 1);
		int x = xo + 1;
		int y = yo + 2;
		int dx = size - 6;
		while (y + 1 < yo + size) {
			g.setColor(edge1);
			g.drawLine(x, y, x + 1, y);
			g.drawLine(x, y + 1, x + 1, y + 1);
			if (0 < dx) {
				g.setColor(fill);
				g.drawLine(x + 2, y, x + 1 + dx, y);
				g.drawLine(x + 2, y + 1, x + 1 + dx, y + 1);
			}
			g.setColor(edge2);
			g.drawLine(x + dx + 2, y, x + dx + 3, y);
			g.drawLine(x + dx + 2, y + 1, x + dx + 3, y + 1);
			x += 1;
			y += 2;
			dx -= 2;
		}
		g.setColor(edge1);
		g.drawLine(xo + (size / 2), yo + size - 1, xo + (size / 2), yo + size
				- 1);
	}

	/**
	 * Internally used by paintIcon.
	 * @param g g
	 * @param xo xo
	 * @param yo yo
	 */
	private void drawUpArrow(Graphics g, int xo, int yo) {
		g.setColor(edge1);
		int x = xo + (size / 2);
		g.drawLine(x, yo, x, yo);
		x--;
		int y = yo + 1;
		int dx = 0;
		while (y + 3 < yo + size) {
			g.setColor(edge1);
			g.drawLine(x, y, x + 1, y);
			g.drawLine(x, y + 1, x + 1, y + 1);
			if (0 < dx) {
				g.setColor(fill);
				g.drawLine(x + 2, y, x + 1 + dx, y);
				g.drawLine(x + 2, y + 1, x + 1 + dx, y + 1);
			}
			g.setColor(edge2);
			g.drawLine(x + dx + 2, y, x + dx + 3, y);
			g.drawLine(x + dx + 2, y + 1, x + dx + 3, y + 1);
			x -= 1;
			y += 2;
			dx += 2;
		}
		g.setColor(edge1);
		g.drawLine(xo, yo + size - 3, xo + 1, yo + size - 3);
		g.setColor(edge2);
		g.drawLine(xo + 2, yo + size - 2, xo + size - 1, yo + size - 2);
		g.drawLine(xo, yo + size - 1, xo + size, yo + size - 1);
	}
}
