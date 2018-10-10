package danfulea.utils.gadgets;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.FontMetrics;
import java.awt.Graphics;

import javax.swing.JPanel;

/**
 * 
 * The panel containing color adjustment elements
 * 
 * @author Dan Fulea, 06 AUG. 2016
 *
 */
public class ChangeColorPanel extends JPanel {
	private static final long serialVersionUID = 1L;

	private Dimension sizePan = new Dimension(350, 100);//default
	private int rectW, rectH;
	private int R = 0, G = 0, B = 0, A = 0;
	private boolean showText = true;
	private String text = "Java(tm) Power";

	/**
	 * Set the text to be displayed
	 * 
	 * @param text
	 *            the text
	 */
	public void setText(String text) {
		this.text = text;
	}

	/**
	 * Sets if text should be displayed or not.
	 * 
	 * @param showText
	 *            true or false
	 */
	public void setShowText(boolean showText) {
		this.showText = showText;
	}

	/**
	 * Sets R,G,B,A for color.
	 * 
	 * @param R
	 *            R component
	 * @param G
	 *            G component
	 * @param B
	 *            B component
	 * @param A
	 *            A component
	 */
	public void setRGBA(int R, int G, int B, int A) {
		this.R = R;
		this.G = G;
		this.B = B;
		this.A = A;
	}

	/**
	 * 
	 * @return panel dimension
	 */
	public Dimension getSizePan() {
		return sizePan;
	}

	/**
	 * 
	 * @param sizePan
	 *            the dimension
	 */
	public void setSizePan(Dimension sizePan) {
		this.sizePan = sizePan;
	}

	/**
	 * Overrides default paintComponent method.
	 */
	public void paintComponent(Graphics g) {

		super.paintComponent(g);// paint bkg
		// without it, artefacts appear and alpha becomes inadequate!!

		rectW = sizePan.width;
		rectH = sizePan.height;

		Color c = new Color(R, G, B, A);

		Font largefont;
		FontMetrics fmlarge;
		String str = "";
		int iascent = 0;
		int idescent = 0;

		if (showText) {

			g.setColor(Color.BLACK);
			str = text;//"Java Power";
			int ifont = 8;
			int ih = 0;
			int iw = 0;

			do {
				ifont++;
				largefont = new Font("Monospaced", Font.BOLD, ifont);
				g.setFont(largefont);
				fmlarge = g.getFontMetrics();
				ih = fmlarge.getHeight();
				iw = fmlarge.stringWidth(str);
				iascent = fmlarge.getAscent();
				idescent = fmlarge.getDescent();
			} while ((ih < rectH / 3) && (iw < rectW));

			g.drawString(str, (rectW - iw) / 2, rectH / 3 + rectH / 6 + (iascent + idescent) / 4);

			return;
		}

		g.drawRect(0, 0, rectW - 1, rectH - 1);
		g.setColor(c);
		g.fillRect(0 + 1, 0 + 1, rectW - 2, rectH - 2);

		largefont = new Font("Monospaced", Font.BOLD, 8);// 16);
		fmlarge = g.getFontMetrics();
		g.setFont(largefont);
		g.setColor(Color.BLACK);

		str = text;//"Java power";
		iascent = fmlarge.getAscent();
		idescent = fmlarge.getDescent();
		int ileading = fmlarge.getLeading();

		g.drawString(str, 0, rectH - iascent - idescent - ileading);
	}
}
