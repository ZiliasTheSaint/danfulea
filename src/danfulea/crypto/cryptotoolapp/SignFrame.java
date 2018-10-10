package danfulea.crypto.cryptotoolapp;

import java.awt.BorderLayout;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.KeyEvent;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.io.BufferedInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
//import java.security.InvalidAlgorithmParameterException;
import java.security.Signature;
import java.util.ResourceBundle;

import javax.swing.BoxLayout;
import javax.swing.JButton;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTextField;
import javax.swing.JToolBar;

import danfulea.utils.FrameUtilities;

/**
 * The sign file frame.
 * 
 * @author Dan Fulea, 24 APR. 2015
 */
@SuppressWarnings("serial")
public class SignFrame extends JFrame implements ActionListener, Runnable {

	private Thread computationTh = null;// computation thread!
	private Thread statusTh = null;// status display thread!
	private int delay = 100;
	private int frameNumber = -1;
	private String statusRunS = "";

	private static final Dimension PREFERRED_SIZE = new Dimension(700, 200);
	private JCrypto2 mf;
	private static final String BASE_RESOURCE_CLASS = "danfulea.crypto.cryptotoolapp.resources.JCrypto2Resources";
	private ResourceBundle resources;
	private JLabel statusL = new JLabel("Waiting...");
	private JTextField inputTf = new JTextField(25);
	private JButton browseB = new JButton("Browse...");
	private JButton signB = new JButton("Sign");
	private File infile = null;

	private String command = "";
	public static final String BROWSE_COMMAND = "BROWSE";
	public static final String SIGN_COMMAND = "SIGN";

	/**
	 * Constructor.
	 * 
	 * @param mf
	 *            the main frame.
	 */
	public SignFrame(JCrypto2 mf) {
		resources = ResourceBundle.getBundle(BASE_RESOURCE_CLASS);
		this.setTitle(resources.getString("SignFrame.NAME"));

		this.mf = mf;

		createGUI();

		setDefaultLookAndFeelDecorated(true);

		FrameUtilities.createImageIcon(this.resources.getString("form.icon.url"), this);
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
	 * Add buttons (in fact a panel and a label) to status bar.
	 * 
	 * @param toolBar
	 *            the status bar.
	 */
	private void addButtons(JToolBar toolBar) {
		JPanel toolP = new JPanel();
		toolP.setLayout(new FlowLayout(FlowLayout.LEFT, 5, 1));

		toolP.add(statusL);
		toolBar.add(toolP);
		statusL.setText("Waiting for your action!");
	}

	/**
	 * This method is called from within the constructor to initialize the form.
	 */
	private void createGUI() {
		inputTf.setEditable(false);

		this.browseB.setMnemonic(KeyEvent.VK_B);
		this.browseB.setActionCommand(BROWSE_COMMAND);
		this.browseB.addActionListener(this);

		this.signB.setMnemonic(KeyEvent.VK_S);
		this.signB.setActionCommand(SIGN_COMMAND);
		this.signB.addActionListener(this);

		JLabel label = new JLabel("Input file");
		JPanel pan = new JPanel();
		pan.setLayout(new BoxLayout(pan, BoxLayout.PAGE_AXIS));

		JPanel p0 = new JPanel();
		p0.setLayout(new FlowLayout(FlowLayout.CENTER, 20, 2));
		p0.add(label);
		p0.add(this.inputTf);
		p0.add(this.browseB);
		p0.add(this.signB);
		p0.setBackground(JCrypto2.fundal);
		pan.add(p0);

		JScrollPane scroller = new JScrollPane(pan);
		scroller.setOpaque(true);

		JPanel content = new JPanel(new BorderLayout());
		content.add(scroller);
		// Create the toolbar.
		JToolBar toolBar = new JToolBar();
		toolBar.setFloatable(false);
		addButtons(toolBar);
		content.add(toolBar, BorderLayout.PAGE_END);

		setContentPane(content);
		content.setOpaque(true); // content panes must be opaque
		pack();
	}

	/**
	 * Basic actions are performed here.
	 */
	public void actionPerformed(ActionEvent evt) {
		command = evt.getActionCommand();
		if (command.equals(BROWSE_COMMAND)) {
			browse();
		} else if (command.equals(SIGN_COMMAND)) {
			// sign();
			statusRunS = resources.getString("status.computing");
			startThread();
		}
	}

	/**
	 * Start the threads.
	 */
	private void startThread() {
		if (computationTh == null) {

			computationTh = new Thread(this);
			computationTh.start();// Allow one simulation at time!
		}

		if (statusTh == null) {
			statusTh = new Thread(this);
			statusTh.start();
		}
	}

	/**
	 * Stop threads.
	 */
	private void stopThread() {
		statusTh = null;
		frameNumber = 0;

		if (computationTh == null) {
			return;
		}
		computationTh = null;
	}

	/**
	 * Run threads.
	 */
	public void run() {
		Thread.currentThread().setPriority(Thread.NORM_PRIORITY);

		long startTime = System.currentTimeMillis();
		Thread currentThread = Thread.currentThread();

		while (currentThread == statusTh) {// if thread is status display
											// Thread!!
			frameNumber++;
			if (frameNumber % 2 == 0)
				statusL.setText(statusRunS + ".....");
			else
				statusL.setText(statusRunS);

			// Delay
			try {
				startTime += delay;
				Thread.sleep(Math.max(0, startTime - System.currentTimeMillis()));
			} catch (InterruptedException e) {
				break;
			}
		}

		if (currentThread == computationTh)
			sign();

		stopThread();
	}

	/**
	 * Gets the filename, without extension.
	 * 
	 * @param file
	 *            the file
	 * @return the file without extension
	 */
	private String getFileWithoutExtension(File file) {
		String name = file.getName();
		int lastIndexOf = name.lastIndexOf(".");
		if (lastIndexOf == -1) {
			return ""; // empty extension
		}
		return name.substring(0, lastIndexOf);
	}

	/**
	 * Sign the file, i.e. it generates its signature file.
	 */
	private void sign() {
		// System.out.println("Fire!");
		if (infile == null) {
			JOptionPane.showMessageDialog(this, "Please load a valid file", "No data to sign",
					JOptionPane.ERROR_MESSAGE);

			stopThread();
			statusL.setText("Waiting for your action!");
			return;
		}

		try {
			// String filename=infile.getName();
			// String extension = getFileExtension(infile);
			String fname = getFileWithoutExtension(infile);
			// System.out.println(filename+"; "+extension+"; "+fname);

			// Create a Signature object and initialize it with the private key
			Signature dsa = Signature.getInstance("SHA256withRSA");// , "SUN");
			dsa.initSign(mf.privatek);

			// Update and sign the data
			// String datafile=inputTf.getText();
			FileInputStream fis = new FileInputStream(infile);// (datafile);
			BufferedInputStream bufin = new BufferedInputStream(fis);
			byte[] buffer = new byte[1024];
			int len;
			while (bufin.available() != 0) {
				len = bufin.read(buffer);
				dsa.update(buffer, 0, len);
			}
			;
			bufin.close();

			// Now that all the data to be signed has been read in, generate a
			// signature for it
			byte[] realSig = dsa.sign();

			// save the signature in a file
			String sigfilename = fname + "_sig.sig";
			FileOutputStream sigfos = new FileOutputStream(sigfilename);
			sigfos.write(realSig);
			sigfos.close();

			statusL.setText("check application directory for the requested output file!");

		} catch (Exception e) {
			e.printStackTrace();
			stopThread();
			statusL.setText("Waiting for your action!");
		}
	}

	/**
	 * Open the file to be signed.
	 */
	private void browse() {
		String currentDir = System.getProperty("user.dir");
		JFileChooser chooser = new JFileChooser(currentDir);
		chooser.setFileSelectionMode(JFileChooser.FILES_ONLY);

		int returnVal = chooser.showOpenDialog(this);// parent=this frame
		if (returnVal == JFileChooser.APPROVE_OPTION) {
			String infilename = chooser.getSelectedFile().toString();

			statusL.setText("open: " + infilename);
			infile = chooser.getSelectedFile();
			inputTf.setText(infilename);
		}
	}
}
