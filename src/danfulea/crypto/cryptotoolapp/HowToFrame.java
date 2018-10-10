package danfulea.crypto.cryptotoolapp;

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
 * The HowTo window displays some tips about the application.
 * 
 * 
 * @author Dan Fulea, 19 MAY. 2011
 */

public class HowToFrame extends JFrame {

	private static final long serialVersionUID = 1L;
	private static final Dimension PREFERRED_SIZE = new Dimension(700, 700);	
	private static final String BASE_RESOURCE_CLASS = "danfulea.crypto.cryptotoolapp.resources.JCrypto2Resources";
	
	private JCrypto2 mf;
	private ResourceBundle resources;
	private JScrollPane jScrollPane1 = new JScrollPane();
	private JTextArea textArea = new JTextArea();

	/**
	 * Constructor. It is connected to the main window.
	 * 
	 * @param mf
	 *            the main frame
	 */
	public HowToFrame(JCrypto2 mf) {
		resources = ResourceBundle.getBundle(BASE_RESOURCE_CLASS);
		this.setTitle(resources.getString("HowTo.NAME"));
		// this.setResizable(false);

		this.mf = mf;

		textArea.setBackground(JCrypto2.textAreaBkgColor);
		textArea.setForeground(JCrypto2.textAreaForeColor);
		createGUI();

		setDefaultLookAndFeelDecorated(true);
		FrameUtilities.createImageIcon(this.resources.getString("form.icon.url"), this);
		// createImageIcon(this.resources.getString("form.icon.url"));

		FrameUtilities.centerFrameOnScreen(this);
		// centerScreen();

		setVisible(true);
		mf.setEnabled(false);

		setDefaultCloseOperation(JFrame.DO_NOTHING_ON_CLOSE);// not necessary,
																// exit normal!
		addWindowListener(new WindowAdapter() {
			public void windowClosing(WindowEvent e) {
				attemptExit();
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
	 * Exit method
	 */

	private void attemptExit() {
		mf.setEnabled(true);
		dispose();
	}

	/**
	 * This method is called from within the constructor to initialize the form.
	 */
	private void createGUI() {

		jScrollPane1.setBorder(new javax.swing.border.TitledBorder(
				new javax.swing.border.LineBorder(new java.awt.Color(0, 51, 255), 1, true),
				this.resources.getString("HowTo.title"), javax.swing.border.TitledBorder.CENTER,
				javax.swing.border.TitledBorder.TOP));
		jScrollPane1.setHorizontalScrollBarPolicy(javax.swing.ScrollPaneConstants.HORIZONTAL_SCROLLBAR_ALWAYS);
		jScrollPane1.setVerticalScrollBarPolicy(javax.swing.ScrollPaneConstants.VERTICAL_SCROLLBAR_ALWAYS);
		jScrollPane1.setAutoscrolls(true);
		textArea.setColumns(1);
		textArea.setEditable(false);

		textArea.setLineWrap(true);
		textArea.setRows(10);
		textArea.setText(this.resources.getString("HowTo"));
		textArea.setWrapStyleWord(true);
		jScrollPane1.setViewportView(textArea);

		JPanel jPanel2 = new JPanel();
		jPanel2.setLayout(new java.awt.BorderLayout());
		jPanel2.add(jScrollPane1, java.awt.BorderLayout.CENTER);
		jPanel2.setBackground(JCrypto2.fundal);

		getContentPane().add(jPanel2, java.awt.BorderLayout.CENTER);
		pack();
	}
}
