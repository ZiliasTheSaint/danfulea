package danfulea.utils;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.Toolkit;
import java.awt.Window;
import java.awt.event.ActionListener;
import java.io.InputStream;

import javax.swing.BorderFactory;
import javax.swing.ImageIcon;
import javax.swing.JButton;
import javax.swing.JFrame;
import javax.swing.border.Border;
import javax.swing.border.TitledBorder;

/**
 * JFrame utilities class!
 * 
 * @author Dan Fulea, 14 APR. 2011
 * 
 */
public class FrameUtilities {

	/**
	 * The size of buffer of 1 Mb which should be enough for icon images or
	 * regular images.
	 */
	private static final int BUFF_SIZE = 1024 * 1000;
	//final on primitives make them constants!

	/**
	 * the array of bytes buffer of BUFF_SIZE size.
	 */
	private static final byte[] tmp = new byte[BUFF_SIZE];
	//final in Java affects the variable, it has nothing to do with the object you are assigning to it.
	//we can modify the content of final tmp object (the array) just fine. The final variable of type array 
	//will always refer to the same array (kind of NO COPY, always ONE tmp array).
	//Here is also static in order to be accessed from static context!
	//there are cases where final variables are passed from constructor, so it is an instance final variable!

	//You are always allowed to initialize a final variable. The compiler makes sure that
	//you can do it only once. final is only about the reference itself, and not about the contents of the referenced object.
	//You can assign value to the final variable only one time. So compiler decides constructor is good choice!
	
	//Rule is that if you have initialized object to final reference variable then you cannot change it to 
	//refer different ,say, ArrayList object. (in this case ArrayList)
     //final class cannot be subclassed
     //final methods cannot be overridden.
	
	/**
	 * Setting the window image icon. Requires the URL of image file and the
	 * actual frame. Since we use a static buffer, the synchronized method is necessary.
	 * 
	 * @param URLstr
	 *            the URL of image icon
	 * @param frame
	 *            the frame to apply the icon
	 * 
	 */
	public synchronized static void createImageIcon(String URLstr, JFrame frame) {
		// tmp buffer should be global variable. Thus it does not have to store
		// buffer in RAM memory every
		// time this method is called!=>optimization!
		int size = 0;
		try {
			InputStream is = frame.getClass().getResourceAsStream(URLstr);

			// fill tmp buffer with data starting from offset size=0.
			while (is.available() > 0) {
				is.read(tmp, size, 1);// Reads up to 1 bytes of data from the
										// input stream into an array of bytes
				size++;// increment the offset
			}
			is.close();

			// copy data from buffer into actual bytes array.
			byte[] data = new byte[size];
			System.arraycopy(tmp, 0, data, 0, size);

			// build icon image
			ImageIcon icon = new ImageIcon(data);
			frame.setIconImage(icon.getImage());
		} catch (Exception exc) {
			exc.printStackTrace();
		}
	}

	/**
	 * Getting an image icon. Requires the URL of image file and the actual
	 * frame. Since we use a static buffer, the synchronized method is necessary.
	 * 
	 * @param URLstr
	 *            the URL of image icon
	 * @param frame
	 *            the frame to apply the icon
	 * @return the image icon
	 */
	public synchronized static ImageIcon getImageIcon(String URLstr, JFrame frame) {
		ImageIcon icon = null;

		int size = 0;
		try {
			InputStream is = frame.getClass().getResourceAsStream(URLstr);

			while (is.available() > 0) {
				is.read(tmp, size, 1);
				size++;
			}
			is.close();

			byte[] data = new byte[size];
			System.arraycopy(tmp, 0, data, 0, size);

			icon = new ImageIcon(data);
		} catch (Exception exc) {

		}

		return icon;
	}

	/**
	 * Centering the frame on screen.
	 * 
	 * @param window
	 *            the window
	 */
	public static void centerFrameOnScreen(Window window) {
		positionFrameOnScreen(window, 0.5D, 0.5D);
	}

	/**
	 * Setting up the frame position on screen.
	 * 
	 * @param window
	 *            the window
	 * @param d
	 *            anchor width fraction
	 * @param d1
	 *            anchor height fraction
	 */
	public static void positionFrameOnScreen(Window window, double d, double d1) {
		Dimension dimension = Toolkit.getDefaultToolkit().getScreenSize();
		Dimension dimension1 = window.getSize();
		int i = Math.max(dimension.width - dimension1.width, 0);
		int j = Math.max(dimension.height - dimension1.height, 0);
		int k = (int) (d * (double) i);
		int l = (int) (d1 * (double) j);
		window.setBounds(k, l, dimension1.width, dimension1.height);
	}

	/**
	 * Customizing buttons with image, action, tooltip and text.
	 * 
	 * @param imageName
	 *            imageName URL string
	 * @param actionCommand
	 *            action command
	 * @param toolTipText
	 *            tooltip text
	 * @param altText
	 *            button caption text
	 * @param frame
	 *            frame containing this button component
	 * @param al
	 *            action listener
	 * @return JButton the button
	 */
	public static JButton makeButton(String imageName, String actionCommand, String toolTipText, String altText,
			JFrame frame, ActionListener al) {

		ImageIcon img = FrameUtilities.getImageIcon(imageName, frame);
		JButton button = new JButton();
		button.setActionCommand(actionCommand);
		button.setToolTipText(toolTipText);
		button.addActionListener(al);
		button.setIcon(img);
		button.setText(altText);
		// button.setOpaque(true);//this should be disabled...
		// otherwise, Nimbus LF shows wrong button edges!!!
		return button;
	}

	/**
	 * Setting a standard panel border with a title.
	 * 
	 * @param title
	 *            title border text
	 * 
	 * @return the TitledBorder
	 */
	public static TitledBorder getGroupBoxBorder(String title) {
		return getGroupBoxBorder(title, Color.BLACK);
	}

	/**
	 * Setting a standard panel border with a title and custom color for the
	 * title.
	 * 
	 * @param title
	 *            title border text
	 * 
	 * @param textColor
	 *            title text color
	 * 
	 * @return the TitledBorder
	 */
	public static TitledBorder getGroupBoxBorder(String title, Color textColor) {
		Color c = Color.BLACK;
		Border LINE_BORDER = BorderFactory.createLineBorder(c, 2);
		TitledBorder tb = BorderFactory.createTitledBorder(LINE_BORDER, title);
		tb.setTitleColor(textColor);//
		Font fnt = tb.getTitleFont();
		Font f = fnt.deriveFont(Font.BOLD);
		tb.setTitleFont(f);
		return tb;
	}

}
