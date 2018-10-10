package danfulea.utils.table;

import java.awt.Color;
import java.awt.Component;
import java.awt.Graphics;

import javax.swing.Icon;

/**
 * The blank icon.
 * Based on example code written by Joris Van den Bogaert see esus.com
 * @author Dan Fulea, 22 Mar. 2005
 */
public class BlankIcon implements Icon {
	private Color fillColor;
	private int size;

	/**
	 * Constructor
	 */
	public BlankIcon() {
		this(null, 11);
	}

	/**
	 * Constructor
	 * @param color the color
	 * @param size the icon size
	 */
	public BlankIcon(Color color, int size) {
		fillColor = color;
		this.size = size;
	}

	/**
	 * Paint specific method.
	 */
	public void paintIcon(Component c, Graphics g, int x, int y) {
		if (fillColor != null) {
			g.setColor(fillColor);
			g.drawRect(x, y, size - 1, size - 1);
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
}
