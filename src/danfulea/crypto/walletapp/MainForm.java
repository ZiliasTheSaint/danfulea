package danfulea.crypto.walletapp;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.InputStream;
import java.io.OutputStream;
import java.util.ResourceBundle;

import javax.swing.BoxLayout;
import javax.swing.ButtonGroup;
import javax.swing.JButton;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JRadioButton;
import javax.swing.JScrollPane;
import javax.swing.JTextArea;
import javax.swing.JToolBar;

import danfulea.crypto.PBE_AESCipher;
import danfulea.utils.AboutFrame;
import danfulea.utils.ExampleFileFilter;
import danfulea.utils.FrameUtilities;
import danfulea.utils.LookAndFeel;

/**
 * GUI for running PBE_AESCipher. To create a jar file, choose: 
* this java file (plus all java files from this package), its resources (java + image file), PBE_AESCipher from danfulea.crypto, 
* Convertor.java from danfulea.math package, DanfuleaResources form danfulea.resources (images not needed), 
* FrameUtilities, AboutFrame, LookAndFeel and ExampleFileFilter from danfulea.utils. As manifest choose the file from 
* META-inf-MyWallet.
 * 
 * @author Dan Fulea, 29 MAR. 2010
 */

public class MainForm extends JFrame implements ActionListener, Runnable {

	private static final long serialVersionUID = 1L;
	protected Color bkgColor = new Color(230, 255, 210, 255);//new Color(180, 220, 150, 255);//new Color(230, 255, 210, 255);
	private JLabel statusL = new JLabel("Waiting...");
	//private static final Border STANDARD_BORDER = BorderFactory.createEmptyBorder(5, 5, 5, 5);
	//private static final Color c = Color.BLACK;
	//private static final Border LINE_BORDER = BorderFactory.createLineBorder(c, 2);
	private static final Dimension PREFERRED_SIZE = new Dimension(680, 550);
	private static final String FILESEPARATOR = System.getProperty("file.separator");
	private static final String currentDir = System.getProperty("user.dir");
	private static final String dataDir = "Data";

	private boolean encB = false;// decrypt mode, true=encrypt mode
	private Thread cryptTh;
	private File infile = null;
	protected String password = null;
	private String infilename = null;

	private static final String BASE_RESOURCE_CLASS = "danfulea.crypto.walletapp.resources.JWalletResources";
	protected ResourceBundle resources;
	private static final String EXIT_COMMAND = "EXIT";
	private static final String ABOUT_COMMAND = "ABOUT";
	private static final String OPEN_COMMAND = "OPEN";
	private static final String SAVE_COMMAND = "SAVE";
	private static final String NEW_COMMAND = "NEW";
	private static final String HELP_COMMAND = "HELP";

	private JRadioButton key128, key256;

	private JTextArea textArea = new JTextArea();

	/**
	 * Constructor of main frame!
	 */
	public MainForm() {
		super("JWallet");
		resources = ResourceBundle.getBundle(BASE_RESOURCE_CLASS);

		addWindowListener(new WindowAdapter() {
			public void windowClosing(WindowEvent e) {
				dispose();
				System.exit(0);
			}
		});

		JMenuBar menuBar = createMenuBar(resources);
		setJMenuBar(menuBar);

		createGUI();
		setDefaultLookAndFeelDecorated(true);
		FrameUtilities.createImageIcon(this.resources.getString("form.icon.url"), this);
		FrameUtilities.centerFrameOnScreen(this);
		setVisible(true);
	}

	/**
	 * Setting up the frame size.
	 */
	public Dimension getPreferredSize() {
		return PREFERRED_SIZE;
	}

	/**
	 * GUI creation.
	 */
	private void createGUI() {
		JPanel content = new JPanel(new BorderLayout());
		JPanel mainPanel = createMainPanel();
		content.add(mainPanel);
		// Create the toolBar
		JToolBar toolBar = new JToolBar();
		toolBar.setFloatable(false);
		initToolBar(toolBar);
		content.add(toolBar, BorderLayout.PAGE_START);
		// Create the statusbar.
		JToolBar statusBar = new JToolBar();
		statusBar.setFloatable(false);
		initStatusBar(statusBar);
		content.add(statusBar, BorderLayout.PAGE_END);

		setContentPane(content);
		content.setOpaque(true); // content panes must be opaque
		pack();
		textArea.requestFocusInWindow();
	}

	/**
	 * Initiates the tool bar!
	 * 
	 * @param toolBar
	 *            the toolbar
	 */
	private void initToolBar(JToolBar toolBar) {
		JButton button = null;
		JPanel toolP = new JPanel();
		toolP.setLayout(new FlowLayout(FlowLayout.LEFT, 5, 1));

		String buttonName = resources.getString("menu.file.new");
		String buttonToolTip = resources.getString("menu.file.new.toolTip");
		String buttonIconName = resources.getString("new.icon");
		button = FrameUtilities.makeButton(buttonIconName, NEW_COMMAND, buttonToolTip, buttonName, this, this);
		Character mnemonic = (Character) resources.getObject("menu.file.new.mnemonic");
		button.setMnemonic(mnemonic.charValue());
		toolP.add(button);

		buttonName = resources.getString("menu.file.open");
		buttonToolTip = resources.getString("menu.file.open.toolTip");
		buttonIconName = resources.getString("open.icon");
		button = FrameUtilities.makeButton(buttonIconName, OPEN_COMMAND, buttonToolTip, buttonName, this, this);
		mnemonic = (Character) resources.getObject("menu.file.open.mnemonic");
		button.setMnemonic(mnemonic.charValue());
		toolP.add(button);

		buttonName = resources.getString("menu.file.save");
		buttonToolTip = resources.getString("menu.file.save.toolTip");
		buttonIconName = resources.getString("save.icon");
		button = FrameUtilities.makeButton(buttonIconName, SAVE_COMMAND, buttonToolTip, buttonName, this, this);
		mnemonic = (Character) resources.getObject("menu.file.save.mnemonic");
		button.setMnemonic(mnemonic.charValue());
		toolP.add(button);

		JPanel toolP1 = new JPanel();
		toolP1.setLayout(new FlowLayout(FlowLayout.RIGHT, 5, 1));

		buttonName = resources.getString("menu.help.about");
		buttonToolTip = resources.getString("menu.help.about.toolTip");
		buttonIconName = resources.getString("about.icon");
		button = FrameUtilities.makeButton(buttonIconName, ABOUT_COMMAND, buttonToolTip, buttonName, this, this);
		mnemonic = (Character) resources.getObject("menu.help.about.mnemonic");
		button.setMnemonic(mnemonic.charValue());
		toolP1.add(button);

		buttonName = resources.getString("menu.help.help");
		buttonToolTip = resources.getString("menu.help.help.toolTip");
		buttonIconName = resources.getString("help.icon");
		button = FrameUtilities.makeButton(buttonIconName, HELP_COMMAND, buttonToolTip, buttonName, this, this);
		mnemonic = (Character) resources.getObject("menu.help.help.mnemonic");
		button.setMnemonic(mnemonic.charValue());
		toolP1.add(button);

		buttonName = resources.getString("menu.file.exit");
		buttonToolTip = resources.getString("menu.file.exit.toolTip");
		buttonIconName = resources.getString("exit.icon");
		button = FrameUtilities.makeButton(buttonIconName, EXIT_COMMAND, buttonToolTip, buttonName, this, this);
		mnemonic = (Character) resources.getObject("menu.file.exit.mnemonic");
		button.setMnemonic(mnemonic.charValue());
		toolP1.add(button);

		JPanel toolP2 = new JPanel();
		toolP2.setLayout(new BorderLayout());
		toolP2.add(toolP, BorderLayout.WEST);
		toolP2.add(toolP1, BorderLayout.EAST);

		toolP2.setBorder(FrameUtilities.getGroupBoxBorder(""));
		toolBar.add(toolP2);
	}

	/**
	 * Initiates the status bar!
	 * 
	 * @param toolBar
	 *            the status bar
	 */
	private void initStatusBar(JToolBar toolBar) {
		JPanel toolP = new JPanel();
		toolP.setLayout(new FlowLayout(FlowLayout.LEFT, 5, 1));

		toolP.add(statusL);
		toolBar.add(toolP);
		statusL.setText(resources.getString("status.wait"));
	}

	/**
	 * Creating the main panel!
	 * 
	 * @return the main panel
	 */
	private JPanel createMainPanel() {
		key128 = new JRadioButton(resources.getString("key.128"));
		key256 = new JRadioButton(resources.getString("key.256"));
		key128.setOpaque(true);
		key256.setOpaque(true);

		ButtonGroup group = new ButtonGroup();
		group.add(key128);
		group.add(key256);
		// Put the radio buttons in a row in a panel.
		JPanel buttP = new JPanel();
		buttP.setLayout(new FlowLayout(FlowLayout.CENTER, 20, 2));
		JLabel lbl = new JLabel(resources.getString("key.label"));
		buttP.add(lbl);
		buttP.add(key128);
		buttP.add(key256);
		key128.setBackground(bkgColor);
		key256.setBackground(bkgColor);
		buttP.setBackground(bkgColor);
		key256.setSelected(true);
		//PBE_AESCipher.MESSAGE_DIGGEST_ALGORITHM = PBE_AESCipher.MESSAGE_DIGGEST_ALGORITHM_SHA256;

		textArea.setCaretPosition(0);
		textArea.setEditable(true);
		textArea.setText("");
		textArea.setLineWrap(true);
		textArea.setWrapStyleWord(true);
		textArea.requestFocusInWindow();

		JPanel resultP = new JPanel(new BorderLayout());
		resultP.add(new JScrollPane(textArea), BorderLayout.CENTER);
		resultP.setBackground(bkgColor);

		JPanel p1 = new JPanel();
		p1.setLayout(new FlowLayout());
		JLabel label = new JLabel();
		label.setText(resources.getString("mainPanel.textArea.label"));
		p1.add(label);
		p1.setBackground(bkgColor);

		JPanel p2 = new JPanel();
		BoxLayout bld = new BoxLayout(p2, BoxLayout.Y_AXIS);
		p2.setLayout(bld);
		p2.add(buttP, null);
		p2.add(p1, null);
		p2.setBackground(bkgColor);

		JPanel mainP = new JPanel(new BorderLayout());
		mainP.add(p2, BorderLayout.NORTH);
		mainP.add(resultP, BorderLayout.CENTER);
		mainP.setBackground(bkgColor);
		return mainP;
	}

	/**
	 * Setting up actions!
	 */
	public void actionPerformed(ActionEvent evt) {
		String command = evt.getActionCommand();

		if (command.equals(ABOUT_COMMAND)) {
			about();
		} else if (command.equals(EXIT_COMMAND)) {
			attemptExit();
		} else if (command.equals(OPEN_COMMAND)) {
			openFile();
		} else if (command.equals(SAVE_COMMAND)) {
			saveFile();
		} else if (command.equals(NEW_COMMAND)) {
			clear();
		} else if (command.equals(HELP_COMMAND)) {
			help();
		}
	}

	/**
	 * Shows the help window!
	 */
	private void help() {
		new HelpFrame(this);
		textArea.requestFocusInWindow();
	}

	/**
	 * Clear all previously results!
	 */
	private void clear() {
		textArea.selectAll();
		textArea.replaceSelection("");
		textArea.requestFocusInWindow();
		statusL.setText(resources.getString("status.new"));
	}

	/**
	 * Open file!
	 */
	private void openFile() {
		// first clear the notes
		textArea.selectAll();
		textArea.replaceSelection("");

		String ext = resources.getString("file.extension");
		String pct = ".";
		String description = resources.getString("file.description");
		ExampleFileFilter eff = new ExampleFileFilter(ext, description);

		String myDir = currentDir + FILESEPARATOR + dataDir;
		// File select
		JFileChooser chooser = new JFileChooser(myDir);
		chooser.addChoosableFileFilter(eff);
		chooser.setFileSelectionMode(JFileChooser.FILES_ONLY);

		int returnVal = chooser.showOpenDialog(this);// parent=this frame
		if (returnVal == JFileChooser.APPROVE_OPTION) {
			infile = chooser.getSelectedFile();
			infilename = chooser.getSelectedFile().toString();

			int fl = infilename.length();
			String test = infilename.substring(fl - 4);// exstension lookup!!
			String ctest = pct + ext;
			if (test.compareTo(ctest) != 0)
				infilename = chooser.getSelectedFile().toString() + pct + ext;

			infile = new File(infilename);
			statusL.setText(resources.getString("status.open") + infilename);
		} else {
			textArea.requestFocusInWindow();
			return;
		}
		// end File select

		// Password prompt
		new PasswordFrame(PasswordFrame.OPEN_MODE, this);
		// end Password prompt
		textArea.requestFocusInWindow();
	}

	/**
	 * Attempt file opening from the password window!
	 */
	protected void attemptOpenFile()// running from PasswordFrame
	{
		encB = false;// decrypt the important notes in file!
		startThread();
	}

	/**
	 * Save file!
	 */
	private void saveFile() {
		String ext = resources.getString("file.extension");
		String pct = ".";
		String description = resources.getString("file.description");
		ExampleFileFilter eff = new ExampleFileFilter(ext, description);

		String myDir = currentDir + FILESEPARATOR + dataDir;
		// File select
		JFileChooser chooser = new JFileChooser(myDir);
		chooser.addChoosableFileFilter(eff);
		chooser.setFileSelectionMode(JFileChooser.FILES_ONLY);

		int returnVal = chooser.showSaveDialog(this);// parent=this frame
		if (returnVal == JFileChooser.APPROVE_OPTION) {
			infile = chooser.getSelectedFile();
			infilename = chooser.getSelectedFile().toString();

			int fl = infilename.length();
			String test = infilename.substring(fl - 4);// exstension lookup!!
			String ctest = pct + ext;
			if (test.compareTo(ctest) != 0)
				infilename = chooser.getSelectedFile().toString() + pct + ext;

			if (infile.exists()) {
				String title = resources.getString("dialog.overwrite.title");
				String message = resources.getString("dialog.overwrite.message");

				Object[] options = (Object[]) resources.getObject("dialog.overwrite.buttons");
				int result = JOptionPane.showOptionDialog(this, message, title, JOptionPane.YES_NO_OPTION,
						JOptionPane.QUESTION_MESSAGE, null, options, options[0]);
				if (result != JOptionPane.YES_OPTION) {
					textArea.requestFocusInWindow();
					return;
				}

			}

			infile = new File(infilename);
			statusL.setText(resources.getString("status.save") + infilename);
		} else {
			textArea.requestFocusInWindow();
			return;
		}
		// end File select

		// Password prompt
		new PasswordFrame(PasswordFrame.SAVE_MODE, this);
		// end Password prompt
		textArea.requestFocusInWindow();
	}

	/**
	 * Attempt file saving from the password window!
	 */
	protected void attemptSaveFile()// //running from PasswordFrame
	{
		encB = true;// encrypt the important notes in file!
		startThread();
	}

	/**
	 * Starting main thread!
	 */
	private void startThread() {
		if (cryptTh == null) {//prevents multiple thread calls for the same thread
			cryptTh = new Thread(this);
			
			cryptTh.start();
		}
		
	}

	/**
	 * Stop the thread!
	 */
	private void stopThread() {
		cryptTh = null;
	}

	/**
	 * Thread main method! Run!
	 */
	public void run() {
		Thread.currentThread().setPriority(Thread.NORM_PRIORITY);
		if (encB)// encrypt file
		{
			try {
				this.encryptFile();
			} catch (Exception e) {
				JOptionPane.showMessageDialog(this, e.getMessage(), resources.getString("error.encrypt"),
						JOptionPane.ERROR_MESSAGE);
			}
		} else// decrypt file
		{
			try {
				this.decryptFile();
			} catch (Exception e) {
				JOptionPane.showMessageDialog(this, e.getMessage(), resources.getString("error.decrypt"),
						JOptionPane.ERROR_MESSAGE);
			}

		}
		stopThread();
	}

	/**
	 * Application close!
	 */
	private void attemptExit() {

		String title = resources.getString("dialog.exit.title");
		String message = resources.getString("dialog.exit.message");

		Object[] options = (Object[]) resources.getObject("dialog.exit.buttons");
		int result = JOptionPane.showOptionDialog(this, message, title, JOptionPane.YES_NO_OPTION,
				JOptionPane.QUESTION_MESSAGE, null, options, options[0]);
		if (result == JOptionPane.YES_OPTION) {
			dispose();
			System.exit(0);
		}
	}

	/**
	 * Shows the about window!
	 */
	private void about() {
		new AboutFrame(resources);// this);
	}

	/**
	 * Decrypt file!
	 * 
	 * @throws Exception
	 *             may throw this exception
	 */
	private void decryptFile() throws Exception {
		
		PBE_AESCipher aesc = new PBE_AESCipher(password);
		
		if (key256.isSelected()) {
			//PBE_AESCipher.MESSAGE_DIGGEST_ALGORITHM = PBE_AESCipher.MESSAGE_DIGGEST_ALGORITHM_SHA256;
			aesc.setMessageDigestAlgorithm(PBE_AESCipher.MESSAGE_DIGGEST_ALGORITHM_SHA256);
		} else if (key128.isSelected()) {
			//PBE_AESCipher.MESSAGE_DIGGEST_ALGORITHM = PBE_AESCipher.MESSAGE_DIGGEST_ALGORITHM_MD5;
			aesc.setMessageDigestAlgorithm(PBE_AESCipher.MESSAGE_DIGGEST_ALGORITHM_MD5);
		}

		InputStream fin1 = new FileInputStream(infile);
		InputStream fin = aesc.decrypt(fin1);

		StringBuffer sbuf = new StringBuffer();

		byte buf[] = new byte[1024];
		int len;
		while ((len = fin.read(buf)) > 0) {
			sbuf.append(new String(buf, 0, len));
		}
		String sdec = sbuf.toString();

		textArea.selectAll();
		textArea.replaceSelection("");
		textArea.append(sdec);
		textArea.requestFocusInWindow();

		fin.close();
		fin1.close();

		aesc = null;
		password = null;

		statusL.setText(resources.getString("status.decrypt") + infilename);
	}

	/**
	 * Encrypt file!
	 * 
	 * @throws Exception
	 *             may throw this exception
	 */
	private void encryptFile() throws Exception {
		
		PBE_AESCipher aesc = new PBE_AESCipher(password);

		if (key256.isSelected()) {
			//PBE_AESCipher.MESSAGE_DIGGEST_ALGORITHM = PBE_AESCipher.MESSAGE_DIGGEST_ALGORITHM_SHA256;
			aesc.setMessageDigestAlgorithm(PBE_AESCipher.MESSAGE_DIGGEST_ALGORITHM_SHA256);
		} else if (key128.isSelected()) {
			//PBE_AESCipher.MESSAGE_DIGGEST_ALGORITHM = PBE_AESCipher.MESSAGE_DIGGEST_ALGORITHM_MD5;
			aesc.setMessageDigestAlgorithm(PBE_AESCipher.MESSAGE_DIGGEST_ALGORITHM_MD5);
		}
		
		FileOutputStream fos = new FileOutputStream(infile);
		OutputStream os = aesc.encrypt(fos);

		// String to encrypt
		String ss = textArea.getText();
		byte[] bs = ss.getBytes("UTF-8");
		// now use OutputStream (os) to write the data:
		os.write(bs);

		os.close();
		fos.close();

		aesc = null;
		password = null;

		statusL.setText(resources.getString("status.encrypt") + infilename);
	}

	/**
	 * Creating menu items!
	 * 
	 * @param resources
	 *            the resources
	 * @return the menu bar
	 */
	private JMenuBar createMenuBar(ResourceBundle resources) {

		// create the menus
		JMenuBar menuBar = new JMenuBar();

		String label;
		Character mnemonic;

		// first the file menu
		label = resources.getString("menu.file");
		mnemonic = (Character) resources.getObject("menu.file.mnemonic");
		JMenu fileMenu = new JMenu(label, true);
		fileMenu.setMnemonic(mnemonic.charValue());

		label = resources.getString("menu.file.new");
		mnemonic = (Character) resources.getObject("menu.file.new.mnemonic");
		JMenuItem newItem = new JMenuItem(label, mnemonic.charValue());
		newItem.setActionCommand(NEW_COMMAND);
		newItem.addActionListener(this);
		fileMenu.add(newItem);

		label = resources.getString("menu.file.open");
		mnemonic = (Character) resources.getObject("menu.file.open.mnemonic");
		JMenuItem openItem = new JMenuItem(label, mnemonic.charValue());
		openItem.setActionCommand(OPEN_COMMAND);
		openItem.addActionListener(this);
		fileMenu.add(openItem);

		label = resources.getString("menu.file.save");
		mnemonic = (Character) resources.getObject("menu.file.save.mnemonic");
		JMenuItem saveItem = new JMenuItem(label, mnemonic.charValue());
		saveItem.setActionCommand(SAVE_COMMAND);
		saveItem.addActionListener(this);
		fileMenu.add(saveItem);

		fileMenu.addSeparator();

		label = resources.getString("menu.file.exit");
		mnemonic = (Character) resources.getObject("menu.file.exit.mnemonic");
		JMenuItem exitItem = new JMenuItem(label, mnemonic.charValue());
		exitItem.setActionCommand(EXIT_COMMAND);
		exitItem.addActionListener(this);
		fileMenu.add(exitItem);

		// then the help menu
		label = resources.getString("menu.help");
		mnemonic = (Character) resources.getObject("menu.help.mnemonic");
		JMenu helpMenu = new JMenu(label);
		helpMenu.setMnemonic(mnemonic.charValue());

		label = resources.getString("menu.help.about");
		mnemonic = (Character) resources.getObject("menu.help.about.mnemonic");
		JMenuItem aboutItem = new JMenuItem(label, mnemonic.charValue());
		aboutItem.setActionCommand(ABOUT_COMMAND);
		aboutItem.addActionListener(this);
		helpMenu.add(aboutItem);
		helpMenu.addSeparator();

		label = resources.getString("menu.help.help");
		mnemonic = (Character) resources.getObject("menu.help.help.mnemonic");
		JMenuItem helpItem = new JMenuItem(label, mnemonic.charValue());
		helpItem.setActionCommand(HELP_COMMAND);
		helpItem.addActionListener(this);
		helpMenu.add(helpItem);

		// finally, glue together the menu and return it
		menuBar.add(fileMenu);
		menuBar.add(helpMenu);

		return menuBar;
	}

	public static void main(String[] args) {
		LookAndFeel.loadLookAndFeel();
		new MainForm();// MainForm mf = new MainForm();
	}
}
