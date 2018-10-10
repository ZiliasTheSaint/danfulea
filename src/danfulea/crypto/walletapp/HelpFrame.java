package danfulea.crypto.walletapp;

import java.awt.BorderLayout;
import java.awt.Dimension;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.util.ResourceBundle;

import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTextArea;

import danfulea.utils.FrameUtilities;

/**
 * The Help window displays some informations.
 * 
 * @author Dan Fulea, 29 MAR. 2010
 */

public class HelpFrame extends JFrame {

	private static final long serialVersionUID = 1L;
	protected MainForm mf;
	private static final String BASE_RESOURCE_CLASS = "danfulea.crypto.walletapp.resources.JWalletResources";
	protected ResourceBundle resources;
	private static final Dimension PREFERRED_SIZE = new Dimension(600, 400);

	private JTextArea textArea = new JTextArea();

	/**
	 * Constructor. Help window is connected to main window. <br>
	 * 
	 * @param mf
	 *            Main Form (frame)
	 */
	public HelpFrame(MainForm mf) {
		super("Help");
		this.mf = mf;
		resources = ResourceBundle.getBundle(BASE_RESOURCE_CLASS);

		createGUI();

		setDefaultLookAndFeelDecorated(true);
		FrameUtilities.createImageIcon(this.resources.getString("form.icon.url"), this);

		// centerScreen();
		FrameUtilities.centerFrameOnScreen(this);

		setVisible(true);

		mf.setVisible(false);
		final MainForm mff = mf;
		addWindowListener(new WindowAdapter() {
			public void windowClosing(WindowEvent e) {
				mff.setVisible(true);
				dispose();
			}
		});

	}

	/**
	 * Setting up the frame size.
	 */
	public Dimension getPreferredSize() {
		return PREFERRED_SIZE;
	}

	/**
	 * This method is called from within the constructor to initialize the form.
	 */
	private void createGUI() {
		textArea.setCaretPosition(0);
		textArea.setEditable(false);
		textArea.setLineWrap(true);
		textArea.setText(this.resources.getString("Infos"));
		textArea.setWrapStyleWord(true);

		JPanel resultP = new JPanel(new BorderLayout());
		resultP.add(new JScrollPane(textArea), BorderLayout.CENTER);
		resultP.setBackground(mf.bkgColor);

		getContentPane().add(resultP, BorderLayout.CENTER);
		pack();
	}
}
