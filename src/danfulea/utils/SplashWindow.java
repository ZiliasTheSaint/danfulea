package danfulea.utils;

import java.awt.BorderLayout;
import java.awt.Image;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.lang.reflect.InvocationTargetException;

import javax.swing.ImageIcon;
import javax.swing.JLabel;
import javax.swing.JWindow;
import javax.swing.SwingUtilities;

/**
 * Splash window that displays an image for a given time or until
 * the user clicks on it with their mouse. Runs in a new thread therefore the calling application is initialized in parallel.
 * 
 * @author Dan Fulea 27. JUN. 2011 based on code written by 
 * Tony Colston Nov 17, 2000 http://www.javaworld.com/article/2077467/core-java/java-tip-104--make-a-splash-with-swing.html,
 * 
 */
@SuppressWarnings("serial")
public class SplashWindow extends JWindow {
	/** Contains the splash image */
	private JLabel splashLabel;

	/**
	 * Creates a new Splash window and displays it for the specified period.
	 * 
	 * @param splashImg
	 *            The splash image
	 * @param iDisplayMs
	 *            Time in milli-seconds to display splash window
	 */
	public SplashWindow(Image splashImg, int iDisplayMs) {
		initComponents(splashImg, iDisplayMs);
	}

	/**
	 * Initializes the window's GUI components and display the splash window for
	 * the specified period of time.
	 * 
	 * @param splashImg
	 *            The splash image
	 * @param iDisplayMs
	 *            Time in milli-seconds to display splash window
	 */
	private void initComponents(Image splashImg, int iDisplayMs) {
		getContentPane().setLayout(new BorderLayout(0, 0));
		splashLabel = new JLabel(new ImageIcon(splashImg));
		getContentPane().add(splashLabel, BorderLayout.CENTER);

		pack();

		setLocationRelativeTo(null);

		addMouseListener(new MouseAdapter() {
			public void mousePressed(MouseEvent e) {
				setVisible(false);
				dispose();
			}
		});

		final int iPauseMs = iDisplayMs;

		final Runnable closerRunner = new Runnable() {
			public void run() {
				setVisible(false);
				dispose();
			}
		};

		Runnable waitRunner = new Runnable() {
			public void run() {
				try {
					Thread.sleep(iPauseMs);
					SwingUtilities.invokeAndWait(closerRunner);
				} catch (InterruptedException e) { //
				} catch (InvocationTargetException e) { //
				}
			}
		};
		setVisible(true);
		toFront();//put this window on top
		Thread splashThread = new Thread(waitRunner, "SplashThread");
		splashThread.start();
	}
}
