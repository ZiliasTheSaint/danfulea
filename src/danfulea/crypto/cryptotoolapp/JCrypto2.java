package danfulea.crypto.cryptotoolapp;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.KeyEvent;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.security.GeneralSecurityException;
import java.security.InvalidAlgorithmParameterException;
import java.security.KeyFactory;
import java.security.KeyPair;
import java.security.KeyPairGenerator;
import java.security.NoSuchAlgorithmException;
import java.security.PrivateKey;
import java.security.PublicKey;
import java.security.SecureRandom;
import java.security.spec.InvalidKeySpecException;
import java.security.spec.PKCS8EncodedKeySpec;
import java.security.spec.X509EncodedKeySpec;
import java.util.ResourceBundle;

import javax.crypto.Cipher;
import javax.crypto.CipherInputStream;
import javax.crypto.CipherOutputStream;
import javax.crypto.KeyGenerator;
import javax.crypto.SecretKey;
import javax.crypto.spec.IvParameterSpec;
import javax.crypto.spec.SecretKeySpec;
import javax.swing.Box;
import javax.swing.BoxLayout;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTextArea;
import javax.swing.JTextField;
import javax.swing.JToolBar;

import danfulea.crypto.Base64Coder;
import danfulea.utils.AboutFrame;
import danfulea.utils.FrameUtilities;
import danfulea.utils.LookAndFeel;

/**
 * 
 * It is generally not advisable to use a public key encryption algorithm such
 * as RSA to directly encrypt files, since (i) public key encryption is slow,
 * and (ii) it will only let you encrypt small things.
 * <p>
 * The alternative, and commonly used approach, is to use a shared key algorithm
 * to encrypt/decrypt the files, and then use a public key algorithm to
 * encrypt/decrypt the (randomly generated) key used by the shared key
 * algorithm. This has the benefit of fast file encryption/decryption whilst
 * still requiring a non-shared private key to get access to the key needed to
 * decrypt the files.
 * </p>
 * <p>
 * Scenario: 1. Generate a key pair - 2 files namely: public.key and private.key
 * in current folder. 2. Generate (make) an AES key and save/encrypt-it using
 * above public.key into a file, e.q aeskey. 3. Use AES key to encrypt a
 * document and save its encrypted version. 4. Send the public key file, aes key
 * file and encrypted document to JohnDoe. At this stage JohnDoe can also use
 * public key to generate a new AES key and encrypt its document and send AES
 * key file and encrypted document to us. We can decrypt the John's AES key
 * using private key associated to public key and decrypt John's encrypted file
 * using John's decrypted AES key.
 * 
 * 5. Copy private key file to JohnDoe computer in private, in order to give him
 * the ability to decrypt our AES key and use it to decrypt our file. OR better:
 * John Doe creates a new pub-privkey and give us the public key and we use it
 * to encrypt message to John so it can read them.
 * 
 * 6. Now we can exchange messages/files.
 * </p>
 *To create a jar file, choose: 
* this java file (plus all java files from this package), its resources (java + image file), Base64Coder from danfulea.crypto, 
* Convertor.java from danfulea.math package, DanfuleaResources form danfulea.resources (images not needed), 
* FrameUtilities, AboutFrame, and LookAndFeel from danfulea.utils. As manifest choose the file from 
* META-inf-crypto.
 * @author Dan Fulea, 22 APR. 2015
 */

public class JCrypto2 extends JFrame implements ActionListener, Runnable {

	private static final long serialVersionUID = 1L;
	private static final String BASE_RESOURCE_CLASS = "danfulea.crypto.cryptotoolapp.resources.JCrypto2Resources";
	protected ResourceBundle resources;

	public static Color textAreaBkgColor = Color.white;// Color.black;can be changed before object initialization
	public static Color textAreaForeColor = Color.black;// Color.yellow;

	private static final String ENCRYPTION_ALGORITHM = "AES/CBC/PKCS5Padding";// "AES";//"AES/CBC/PKCS5Padding";//
	// for padding:
	private static final int ENCRYPTION_IV_LENGTH = 16;
	private static final String ENCRYPTIONKEY_ALGORITHM = "RSA";
	private static final String ENCRYPTION_KEY = "AES";
	private static final int AES_Key_Size = 256;// 128;// 256;
	
	private byte[] encryptionIV;
	private SecureRandom random;
	// --------------------------	
	private Cipher pkCipher, aesCipher;
	private byte[] aesKey;
	private SecretKeySpec aeskeySpec;

	private String publicKey = "";
	private String privateKey = "";
	private String encryptionKey = "";
	private String publicKeyFilePath = "";
	private String privateKeyFilePath = "";
	private String encryptionKeyFilePath = "";
	private String fileToEncryptPath = "";
	private String encryptedFilePath = "";
	private String decryptedFilePath = "";
	// ==============================Now gui
	private int runIndex = 0;

	private static final Dimension PREFERRED_SIZE = new Dimension(900, 750);
	private JButton encryptFileButton;
	private JButton decryptFileButton;
	private static final String ENCRYPT_FILE = "Encrypt file";
	private static final String DECRYPT_FILE = "Decrypt file";
	
	private JTextField infileTf = new JTextField(25);
	private JTextField outfileTf = new JTextField(15);
	private JButton browseB;
	private static final String BROWSE_FILE = "Browse file";
	private String infilename = "";
	private File infile = null;

	private static final String EXIT_COMMAND = "EXIT";
	private static final String ABOUT_COMMAND = "ABOUT";
	private Thread computationTh = null;// computation thread!
	private Thread statusTh = null;// status display thread!
	private int delay = 100;
	private int frameNumber = -1;
	private String statusRunS = "";
	private String command = "";
	private static final String GENERATE_PAIRKEY = "Generate public-private pair keys into the application folder";
	private static final String PLAINTEXT = "PlainText message:";
	private static final String CIPHERTEXT = "CipherText message:";
	private static final String ENCRYPT = "Encrypt";
	private static final String DECRYPT = "Decrypt";
	private static final String HOWTO_COMMAND = "HOWTO";
	private static final String SIGN_COMMAND = "SIGN";
	private static final String VERIFY_COMMAND = "VERIFY";
	// ====================
	protected static final Color fundal = new Color(230, 255, 210, 255);//new Color(180, 220, 150, 255);//new Color(230, 255, 210, 255);//background color shared among package classes
	private JLabel statusL = new JLabel("Waiting...");
	private JButton decryptButton;
	private JButton encryptButton;

	private JLabel cipherTextLabel = new JLabel(CIPHERTEXT);
	private JTextArea cipherTextArea;
	private JScrollPane cipherTextScroller;
	private JLabel plainTextLabel = new JLabel(PLAINTEXT);
	private JTextArea plainTextArea;
	private JScrollPane plainTextScroller;

	private static final String IMPLEMENTATION = "Encrypt/Decrypt algorithm using (AES/CBC/PKCS5PADDING) algorithm mode and 256 bit AES encryption key";
	private JLabel implementationLabel = new JLabel(IMPLEMENTATION);
	private JButton generateKeyButton;
	private static final int TEXT_ROWS = 5;//for text area. It relates to the size.
	private static final int TEXT_COLS = 20;

	private JLabel pubkeynameLabel = new JLabel("Public-private key pair name: ");

	private JTextField pubkeyNameTf = new JTextField(15);

	private JLabel pubkeyLabel = new JLabel("Use public key for encoding: ");
	private JLabel privkeyLabel = new JLabel("Use private key for decoding: ");
	private JTextField pubkeyTf = new JTextField(15);
	private JTextField privkeyTf = new JTextField(15);

	private static final String LOADPUBKEY = "Load public key";
	private static final String LOADPRIVKEY = "Load private key";
	private JButton loadPubKeyButton = new JButton(LOADPUBKEY);
	private JButton loadPrivKeyButton = new JButton(LOADPRIVKEY);

	private JLabel encryptionkeyLabel = new JLabel("Encryption key for decoding: ");
	private JTextField encryptionkeyTf = new JTextField(15);
	private static final String LOADENCRYPTIONKEY = "Load encryption key";
	private JButton loadEncryptionKeyButton = new JButton(LOADENCRYPTIONKEY);

	private JCheckBox useKeyCh;
	protected PublicKey publick;//use for verify signed file
	protected PrivateKey privatek;//use for sign file

	/**
	 * Constructor
	 */
	public JCrypto2() {

		resources = ResourceBundle.getBundle(BASE_RESOURCE_CLASS);
		this.setTitle(this.resources.getString("Application.NAME"));
		// -------------------------------------------
		addWindowListener(new WindowAdapter() {
			public void windowClosing(WindowEvent e) {
				dispose();
				System.exit(0);
			}
		});

		JMenuBar menuBar = createMenuBar(resources);
		setJMenuBar(menuBar);

		creareGUI();// main GUI
		// initEvent();//main event
		setDefaultLookAndFeelDecorated(true); // cu icoana pe forma
		FrameUtilities.createImageIcon(this.resources.getString("form.icon.url"), this);// "/Jad/images/personal.jpg");
		FrameUtilities.centerFrameOnScreen(this);
		setVisible(true);

		// ===============
		encryptionKey = "encryptionKey.key";
		// =========================
		initializeCipher();// init ciphers
		// ==========================
	}

	/**
	 * Basic actions are set here.
	 */
	public void actionPerformed(ActionEvent evt) {
		command = evt.getActionCommand();
		if (command.equals(GENERATE_PAIRKEY)) {
			statusRunS = resources.getString("status.computing");
			runIndex = 2;
			startThread();
			// this.generateKey();
		} else if (command.equals(ENCRYPT)) {
			this.doEncrypt();
		} else if (command.equals(DECRYPT)) {
			this.doDecrypt();
		} else if (command.equals(ENCRYPT_FILE)) {
			statusRunS = resources.getString("status.computing");
			runIndex = 0;// encB = true;
			startThread();
		} // ENCRYPT_FILE
		else if (command.equals(DECRYPT_FILE)) {
			statusRunS = resources.getString("status.computing");
			runIndex = 1;// encB = false;
			startThread();
		} // DECRYPT_FILE
		else if (command.equals(ABOUT_COMMAND)) {
			about();
		} else if (command.equals(EXIT_COMMAND)) {
			attemptExit();
		} else if (command.equals(BROWSE_FILE)) {
			loadFile();
		} else if (command.equals(LOADPUBKEY)) {
			loadPubKey();
		} else if (command.equals(LOADPRIVKEY)) {
			loadPrivKey();
		} else if (command.equals(HOWTO_COMMAND)) {
			howTo();
		} else if (command.equals(LOADENCRYPTIONKEY)) {
			loadEncryptionKey();
		} else if (command.equals(SIGN_COMMAND)) {
			sign();
		} else if (command.equals(VERIFY_COMMAND)) {
			verify();
		}
	}

	/**
	 * Verify the signed file!
	 */
	private void verify() {
		if ((publicKeyFilePath == null) || (publicKeyFilePath.length() < 1)) {
			JOptionPane.showMessageDialog(this, "Please load the public key of the provider", "No public key",
					JOptionPane.ERROR_MESSAGE);
			return;
		}

		try {
			// read public key
			File pubkeyfile = new File(publicKeyFilePath);
			byte[] encodedKey = new byte[(int) pubkeyfile.length()];
			// new FileInputStream(pubkeyfile).read(encodedKey);
			FileInputStream fis = new FileInputStream(pubkeyfile);
			fis.read(encodedKey);
			// create public key
			X509EncodedKeySpec publicKeySpec = new X509EncodedKeySpec(encodedKey);
			KeyFactory kf = KeyFactory.getInstance(ENCRYPTIONKEY_ALGORITHM);// "RSA");
			publick = kf.generatePublic(publicKeySpec);
			fis.close();
		} catch (Exception e) {
			e.printStackTrace();
			return;
		}

		// System.out.println("Succes!");
		new VerifySignedFile(this);
	}

	/**
	 * Sign the file.
	 */
	private void sign() {
		if ((privateKeyFilePath == null) || (privateKeyFilePath.length() < 1)) {
			JOptionPane.showMessageDialog(this,
					"Please load your private key coresponding to the public key used for encryption", "No private key",
					JOptionPane.ERROR_MESSAGE);
			return;
		}

		try {

			// read private key file
			File privateKeyFile = new File(privateKeyFilePath);
			byte[] encodedKey = new byte[(int) privateKeyFile.length()];
			// new FileInputStream(privateKeyFile).read(encodedKey);
			FileInputStream fis = new FileInputStream(privateKeyFile);
			fis.read(encodedKey);

			// create private key
			PKCS8EncodedKeySpec privateKeySpec = new PKCS8EncodedKeySpec(encodedKey);
			KeyFactory kf = KeyFactory.getInstance(ENCRYPTIONKEY_ALGORITHM);// "RSA");
			privatek = kf.generatePrivate(privateKeySpec);
			fis.close();
		} catch (Exception e) {
			e.printStackTrace();
			return;
		}

		// System.out.println("Succes!");
		new SignFrame(this);
	}

	/**
	 * Load the encryption key.
	 */
	private void loadEncryptionKey() {

		String currentDir = System.getProperty("user.dir");
		JFileChooser chooser = new JFileChooser(currentDir);
		chooser.setFileSelectionMode(JFileChooser.FILES_ONLY);

		int returnVal = chooser.showOpenDialog(this);// parent=this frame
		if (returnVal == JFileChooser.APPROVE_OPTION) {
			String infilename = chooser.getSelectedFile().toString();
			encryptionkeyTf.setText(infilename);
			statusL.setText("open: " + infilename);
			encryptionKeyFilePath = infilename;
		}
	}

	/**
	 * Display How to window.
	 */
	private void howTo() {
		new HowToFrame(this);
	}

	/**
	 * Load the public key.
	 */
	private void loadPubKey() {
		String currentDir = System.getProperty("user.dir");
		JFileChooser chooser = new JFileChooser(currentDir);
		chooser.setFileSelectionMode(JFileChooser.FILES_ONLY);

		int returnVal = chooser.showOpenDialog(this);// parent=this frame
		if (returnVal == JFileChooser.APPROVE_OPTION) {
			String infilename = chooser.getSelectedFile().toString();
			statusL.setText("open: " + infilename);

			publicKeyFilePath = infilename;
			pubkeyTf.setText(infilename);
		}
	}

	/**
	 * Load the private key.
	 */
	private void loadPrivKey() {
		String currentDir = System.getProperty("user.dir");
		JFileChooser chooser = new JFileChooser(currentDir);
		chooser.setFileSelectionMode(JFileChooser.FILES_ONLY);

		int returnVal = chooser.showOpenDialog(this);// parent=this frame
		if (returnVal == JFileChooser.APPROVE_OPTION) {
			String infilename = chooser.getSelectedFile().toString();
			statusL.setText("open: " + infilename);

			privateKeyFilePath = infilename;
			privkeyTf.setText(infilename);
		}
	}

	/**
	 * Load input file.
	 */
	private void loadFile() {
		String currentDir = System.getProperty("user.dir");
		JFileChooser chooser = new JFileChooser(currentDir);
		chooser.setFileSelectionMode(JFileChooser.FILES_ONLY);

		int returnVal = chooser.showOpenDialog(this);// parent=this frame
		if (returnVal == JFileChooser.APPROVE_OPTION) {
			infile = chooser.getSelectedFile();
			infilename = chooser.getSelectedFile().toString();
			infileTf.setText(infilename);
			statusL.setText("open: " + infilename);
		}

	}

	/**
	 * Starting threads.
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
	 * Stopping the threads.
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
			if (runIndex == 0)// enc file
			{
				try {
					this.encryptFile();
				} catch (IOException e) {
					JOptionPane.showMessageDialog(this, e.getMessage(), "IO error: AES.encryptFile()",
							JOptionPane.ERROR_MESSAGE);
				} catch (InvalidAlgorithmParameterException e) {
					JOptionPane.showMessageDialog(this, e.getMessage(), "Invalid parameter error: AES.encryptFile()",
							JOptionPane.ERROR_MESSAGE);
				} catch (Exception e) {
					e.printStackTrace();
					JOptionPane.showMessageDialog(this, e.getMessage(), "Error: AES.encryptFile()",
							JOptionPane.ERROR_MESSAGE);
				}
			} else if (runIndex == 1)// dec file
			{
				try {
					this.decryptFile();
				} catch (IOException e) {
					JOptionPane.showMessageDialog(this, e.getMessage(), "IO error: AES.decryptFile()",
							JOptionPane.ERROR_MESSAGE);
				} catch (InvalidAlgorithmParameterException e) {
					JOptionPane.showMessageDialog(this, e.getMessage(), "Invalid parameter error: AES.decryptFile()",
							JOptionPane.ERROR_MESSAGE);
				} catch (Exception e) {
					e.printStackTrace();
					JOptionPane.showMessageDialog(this, e.getMessage(), "Error: AES.decryptFile()",
							JOptionPane.ERROR_MESSAGE);
				}

			} else if (runIndex == 2) {// keypair generation
				generateKey();
			}
		stopThread();
	}

	/**
	 * Decrypt file.
	 * 
	 * @throws Exception
	 *             can throw this exception
	 */
	private void decryptFile() throws Exception {
		String nulls = "";
		if (outfileTf.getText().compareTo(nulls) == 0) {
			JOptionPane.showMessageDialog(this, "Set a valid name for output file", "No output file",
					JOptionPane.ERROR_MESSAGE);

			statusL.setText("Waiting for your action!");
			return;
		}

		if (infilename.compareTo(nulls) == 0) {
			JOptionPane.showMessageDialog(this, "Set a valid name for input file", "No input file",
					JOptionPane.ERROR_MESSAGE);

			statusL.setText("Waiting for your action!");
			return;
		}

		File original = infile;
		String extension = getFileExtension(original);

		if ((privateKeyFilePath == null) || (privateKeyFilePath.length() < 1)) {
			JOptionPane.showMessageDialog(this, "Please load a valid private key", "No private key",
					JOptionPane.ERROR_MESSAGE);
			return;
		}
		if ((encryptionKeyFilePath == null) || (encryptionKeyFilePath.length() < 1)) {
			JOptionPane.showMessageDialog(this, "Please load a valid encryption key", "No encryption key",
					JOptionPane.ERROR_MESSAGE);
			return;
		}

		// first load encryption key
		loadKey(new File(encryptionKeyFilePath), new File(privateKeyFilePath));

		String currentDir = System.getProperty("user.dir");
		String file_sep = System.getProperty("file.separator");
		String path = currentDir + file_sep;
		decryptedFilePath = path + outfileTf.getText() + extension;
		encryptedFilePath = infilename;

		// decrypt and save
		decrypt(new File(encryptedFilePath), new File(decryptedFilePath));

		JOptionPane.showMessageDialog(this, "Decrypt file finished", "File: " + original.getCanonicalFile().toString(),
				JOptionPane.INFORMATION_MESSAGE);
		statusL.setText("check application directory for the requested output file!");
	}

	/**
	 * Encrypt file.
	 * 
	 * @throws Exception
	 *             can throw this exception
	 */
	private void encryptFile() throws Exception {
		String nulls = "";
		if (outfileTf.getText().compareTo(nulls) == 0) {
			JOptionPane.showMessageDialog(this, "Set a valid name for output file", "No output file",
					JOptionPane.ERROR_MESSAGE);

			statusL.setText("Waiting for your action!");
			return;
		}

		if (infilename.compareTo(nulls) == 0) {
			JOptionPane.showMessageDialog(this, "Set a valid name for input file", "No input file",
					JOptionPane.ERROR_MESSAGE);

			statusL.setText("Waiting for your action!");
			return;
		}

		File original = infile;
		String extension = getFileExtension(original);

		if ((publicKeyFilePath == null) || (publicKeyFilePath.length() < 1)) {
			JOptionPane.showMessageDialog(this, "Please load a valid public key", "No public key",
					JOptionPane.ERROR_MESSAGE);
			statusL.setText("Waiting for your action!");
			return;
		}

		// first generate a new encryption key!
		if (useKeyCh.isSelected())
			makeKey();
		else {
			if ((encryptionKeyFilePath == null) || (encryptionKeyFilePath.length() < 1)) {
				JOptionPane.showMessageDialog(this, "Please load a valid encryption key", "No encryption key",
						JOptionPane.ERROR_MESSAGE);
				statusL.setText("Waiting for your action!");
				return;
			}

			if ((privateKeyFilePath == null) || (privateKeyFilePath.length() < 1)) {
				JOptionPane.showMessageDialog(this, "Please load a valid private key", "No private key",
						JOptionPane.ERROR_MESSAGE);
				statusL.setText("Waiting for your action!");
				return;
			}

			// load key
			loadKey(new File(encryptionKeyFilePath), new File(privateKeyFilePath));
		}

		String currentDir = System.getProperty("user.dir");
		String file_sep = System.getProperty("file.separator");
		String path = currentDir + file_sep;
		fileToEncryptPath = infilename;
		encryptedFilePath = path + outfileTf.getText() + extension;

		// encrypt and save the file
		encrypt(new File(fileToEncryptPath), new File(encryptedFilePath));

		// now encrypt and save the encryption key:
		String encryptedKeyFileS = path + encryptionKey;// "fixed name";
		saveKey(new File(encryptedKeyFileS), new File(publicKeyFilePath));

		JOptionPane.showMessageDialog(this, "Encrypt file finished", "File: " + original.getCanonicalFile().toString(),
				JOptionPane.INFORMATION_MESSAGE);
		statusL.setText("check application directory for the requested output file!");
	}

	/**
	 * Decrypt message.
	 */
	private void doDecrypt() {
		String cipherS = cipherTextArea.getText();
		if ((cipherS == null) || (cipherS.length() < 1)) {
			JOptionPane.showMessageDialog(this, "Please enter some ciphertext to encrypt", "No ciphertext to encrypt",
					JOptionPane.ERROR_MESSAGE);
			return;
		}
		if ((privateKeyFilePath == null) || (privateKeyFilePath.length() < 1)) {
			JOptionPane.showMessageDialog(this, "Please load a valid private key", "No private key",
					JOptionPane.ERROR_MESSAGE);
			return;
		}
		if ((encryptionKeyFilePath == null) || (encryptionKeyFilePath.length() < 1)) {
			JOptionPane.showMessageDialog(this, "Please load a valid encryption key", "No encryption key",
					JOptionPane.ERROR_MESSAGE);
			return;
		}

		try {
			// IV:
			String encIVS = cipherS.substring(0, 24);
			String cipherSS = cipherS.substring(24);

			this.encryptionIV = getBytesFromBase64(encIVS);

			// best raw data-string conversion, otherwise we may obtain bad
			// padding exception
			byte[] cipher2 = getBytesFromBase64(cipherSS);// cipherS);//
															// cipherS.getBytes();

			// load key
			loadKey(new File(encryptionKeyFilePath), new File(privateKeyFilePath));

			// decrypt
			String decrypted = decryptBytes(cipher2);

			plainTextArea.selectAll();
			plainTextArea.replaceSelection("");
			plainTextArea.setText(decrypted);
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	/**
	 * Encrypt message.
	 */
	private void doEncrypt() {
		String message = plainTextArea.getText();
		if ((message == null) || (message.length() < 1)) {
			JOptionPane.showMessageDialog(this, "Please enter some plaintext to encrypt", "No plaintext to encrypt",
					JOptionPane.ERROR_MESSAGE);
			return;
		}

		if ((publicKeyFilePath == null) || (publicKeyFilePath.length() < 1)) {
			JOptionPane.showMessageDialog(this, "Please load a valid public key", "No public key",
					JOptionPane.ERROR_MESSAGE);
			return;
		}
		try {
			// first generate a new encryption key!
			if (useKeyCh.isSelected())
				makeKey();
			else {
				if ((encryptionKeyFilePath == null) || (encryptionKeyFilePath.length() < 1)) {
					JOptionPane.showMessageDialog(this, "Please load a valid encryption key", "No encryption key",
							JOptionPane.ERROR_MESSAGE);
					return;
				}

				if ((privateKeyFilePath == null) || (privateKeyFilePath.length() < 1)) {
					JOptionPane.showMessageDialog(this, "Please load a valid private key", "No private key",
							JOptionPane.ERROR_MESSAGE);
					return;
				}

				// load key
				loadKey(new File(encryptionKeyFilePath), new File(privateKeyFilePath));
			}

			// encrypt
			byte[] cipher = encryptString(message);

			// encrtption IV
			String str = getBase64Text(encryptionIV);
			// encode for best raw data-string conversion, otherwise we may
			// obtain bad padding exception
			String cipherS = getBase64Text(cipher);// new String(cipher);
			String encodedMessage = str + cipherS;

			cipherTextArea.selectAll();
			cipherTextArea.replaceSelection("");
			cipherTextArea.setText(encodedMessage);// (cipherS);

			// cipherTextArea.setText(encodedMessage+"\n"+
			// encodedMessage.substring(0,24)+"\n"+encodedMessage.substring(24));//(cipherS);

			// now encrypt and save the encryption key:
			String currentDir = System.getProperty("user.dir");
			String file_sep = System.getProperty("file.separator");
			String path = currentDir + file_sep;
			String encryptedKeyFileS = path + encryptionKey;// "fixed name";

			saveKey(new File(encryptedKeyFileS), new File(publicKeyFilePath));
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	/**
	 * Preparing to generates the public-private key pair.
	 */
	private void generateKey() {
		publicKey = "publicKey.key";
		privateKey = "privateKey.key";

		String s = pubkeyNameTf.getText();
		if (!s.equals("")) {
			publicKey = s + "_public.key";
			privateKey = s + "_private.key";
		}

		generateKeyPairIntoFiles();

		plainTextArea.selectAll();
		plainTextArea.replaceSelection("");
		plainTextArea.setText("Public key filename: " + publicKey + "; Private key filename: " + privateKey);
		statusL.setText("check application directory for the requested output file!");
	}

	/**
	 * Creates the menu bar.
	 * 
	 * @param resources
	 *            the resources.
	 * @return the menu bar.
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

		label = resources.getString("menu.file.sign");
		mnemonic = (Character) resources.getObject("menu.file.sign.mnemonic");
		JMenuItem signItem = new JMenuItem(label, mnemonic.charValue());
		signItem.setActionCommand(SIGN_COMMAND);
		signItem.addActionListener(this);
		signItem.setToolTipText(resources.getString("menu.file.sign.toolTip"));
		fileMenu.add(signItem);

		label = resources.getString("menu.file.verify");
		mnemonic = (Character) resources.getObject("menu.file.verify.mnemonic");
		JMenuItem verifyItem = new JMenuItem(label, mnemonic.charValue());
		verifyItem.setActionCommand(VERIFY_COMMAND);
		verifyItem.addActionListener(this);
		verifyItem.setToolTipText(resources.getString("menu.file.verify.toolTip"));
		fileMenu.add(verifyItem);
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

		label = resources.getString("menu.help.howTo");
		mnemonic = (Character) resources.getObject("menu.help.howTo.mnemonic");
		JMenuItem howToItem = new JMenuItem(label, mnemonic.charValue());
		howToItem.setActionCommand(HOWTO_COMMAND);
		howToItem.addActionListener(this);
		// howToItem.setIcon(img);
		howToItem.setToolTipText(resources.getString("menu.help.howTo.toolTip"));
		helpMenu.add(howToItem);

		label = resources.getString("menu.help.about");
		mnemonic = (Character) resources.getObject("menu.help.about.mnemonic");
		JMenuItem aboutItem = new JMenuItem(label, mnemonic.charValue());
		aboutItem.setActionCommand(ABOUT_COMMAND);
		aboutItem.addActionListener(this);
		helpMenu.add(aboutItem);
		helpMenu.addSeparator();

		// finally, glue together the menu and return it
		menuBar.add(fileMenu);
		menuBar.add(helpMenu);

		return menuBar;
	}

	/**
	 * Exit the application.
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
	 * Display the about frame.
	 */
	private void about() {
		new AboutFrame(resources);// this);
	}

	/**
	 * Sets up the window size.
	 */
	public Dimension getPreferredSize() {
		return PREFERRED_SIZE;
	}

	/**
	 * GUI creation.
	 */
	private void creareGUI() {
		JPanel content = new JPanel(new BorderLayout());
		JPanel tabs = createTabs();
		// ======
		JScrollPane scroller = new JScrollPane(tabs);
		scroller.setOpaque(true);
		// =======
		content.add(scroller);// (tabs);
		// Create the toolbar.
		JToolBar toolBar = new JToolBar();
		toolBar.setFloatable(false);
		addButtons(toolBar);
		content.add(toolBar, BorderLayout.PAGE_END);

		setContentPane(new JScrollPane(content));// (content);
		content.setOpaque(true); // content panes must be opaque
		pack();

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
	 * Creates main panel.
	 * 
	 * @return the main panel.
	 */
	private JPanel createTabs() {
		// ---------
		useKeyCh = new JCheckBox("Always generate new encryption key for decoding", true);
		useKeyCh.setBackground(fundal);
		// -------------

		JPanel pan = new JPanel();
		pan.setLayout(new BoxLayout(pan, BoxLayout.PAGE_AXIS));

		// Add radio buttons to the panel.
		JPanel p0 = new JPanel();
		p0.setLayout(new FlowLayout(FlowLayout.CENTER));
		p0.add(this.implementationLabel);
		// p0.add(this.useKeyCh);
		p0.setBackground(fundal);
		pan.add(p0);

		JPanel p00 = new JPanel();
		p00.setLayout(new FlowLayout(FlowLayout.CENTER));
		String res00 = "Encryption option: ";
		JLabel label00 = new JLabel(res00);
		p00.add(label00);
		p00.add(this.useKeyCh);
		p00.setBackground(fundal);
		pan.add(p00);

		// Add blank space between components: 5 pixels high.
		Dimension spacer = new Dimension(0, 5);
		pan.add(Box.createRigidArea(spacer));

		JPanel p01 = new JPanel();
		p01.setLayout(new FlowLayout(FlowLayout.CENTER));
		p01.add(this.pubkeynameLabel);
		p01.add(this.pubkeyNameTf);
		p01.setBackground(fundal);
		pan.add(p01);

		pan.add(Box.createRigidArea(spacer));
		// Add button for generating key.
		this.generateKeyButton = new JButton(GENERATE_PAIRKEY);
		this.generateKeyButton.setMnemonic(KeyEvent.VK_G);
		this.generateKeyButton.setActionCommand(GENERATE_PAIRKEY);
		this.generateKeyButton.addActionListener(this);

		this.loadPubKeyButton.setMnemonic(KeyEvent.VK_L);
		this.loadPubKeyButton.setActionCommand(LOADPUBKEY);
		this.loadPubKeyButton.addActionListener(this);

		this.loadPrivKeyButton.setMnemonic(KeyEvent.VK_O);
		this.loadPrivKeyButton.setActionCommand(LOADPRIVKEY);
		this.loadPrivKeyButton.addActionListener(this);

		// encryptionkeyLabel
		this.loadEncryptionKeyButton.setMnemonic(KeyEvent.VK_E);
		this.loadEncryptionKeyButton.setActionCommand(LOADENCRYPTIONKEY);
		this.loadEncryptionKeyButton.addActionListener(this);

		JPanel p04 = new JPanel();
		p04.setLayout(new FlowLayout(FlowLayout.CENTER));
		p04.add(this.generateKeyButton);
		p04.setBackground(fundal);
		pan.add(p04);// this.generateKeyButton);

		// Add spacing:
		pan.add(Box.createRigidArea(spacer));

		JPanel p02 = new JPanel();
		p02.setLayout(new FlowLayout(FlowLayout.CENTER));
		p02.add(this.pubkeyLabel);
		p02.add(this.pubkeyTf);
		p02.add(this.loadPubKeyButton);
		p02.setBackground(fundal);
		pan.add(p02);
		pan.add(Box.createRigidArea(spacer));

		JPanel p03 = new JPanel();
		p03.setLayout(new FlowLayout(FlowLayout.CENTER));
		p03.add(this.privkeyLabel);
		p03.add(this.privkeyTf);
		p03.add(this.loadPrivKeyButton);
		p03.setBackground(fundal);
		pan.add(p03);
		pan.add(Box.createRigidArea(spacer));
		// pubkeyLabel

		// Add (labeled) plaintext area.
		JPanel p065 = new JPanel();
		p065.setLayout(new FlowLayout(FlowLayout.CENTER));
		p065.add(this.plainTextLabel);
		p065.setBackground(fundal);
		pan.add(p065);// this.plainTextLabel);
		this.plainTextArea = new JTextArea(TEXT_ROWS, TEXT_COLS);
		this.plainTextArea.setBackground(JCrypto2.textAreaBkgColor);
		this.plainTextArea.setForeground(JCrypto2.textAreaForeColor);
		this.plainTextScroller = new JScrollPane(this.plainTextArea);
		pan.add(this.plainTextScroller);

		// Add spacing:
		pan.add(Box.createRigidArea(spacer));

		// Add (labeled) ciphertext area.
		JPanel p07 = new JPanel();
		p07.setLayout(new FlowLayout(FlowLayout.CENTER));
		p07.add(this.cipherTextLabel);
		p07.setBackground(fundal);
		pan.add(p07);// this.cipherTextLabel);
		this.cipherTextArea = new JTextArea(TEXT_ROWS, TEXT_COLS);
		this.cipherTextArea.setBackground(JCrypto2.textAreaBkgColor);
		this.cipherTextArea.setForeground(JCrypto2.textAreaForeColor);
		this.cipherTextScroller = new JScrollPane(this.cipherTextArea);
		pan.add(this.cipherTextScroller);

		// Add spacing:
		pan.add(Box.createRigidArea(spacer));

		JPanel p022 = new JPanel();
		p022.setLayout(new FlowLayout(FlowLayout.CENTER));
		p022.add(this.encryptionkeyLabel);
		p022.add(this.encryptionkeyTf);
		p022.add(this.loadEncryptionKeyButton);
		// p022.add(this.useKeyCh);
		p022.setBackground(fundal);
		pan.add(p022);
		pan.add(Box.createRigidArea(spacer));

		// Add panel for holding encrypt, decrypt, and test buttons.
		JPanel processPanel = new JPanel();
		processPanel
				// .setLayout(new BoxLayout(processPanel, BoxLayout.LINE_AXIS));
				.setLayout(new FlowLayout(FlowLayout.CENTER));

		// Add encrypt button.
		this.encryptButton = new JButton(ENCRYPT);
		this.encryptButton.setMnemonic(KeyEvent.VK_N);
		this.encryptButton.setActionCommand(ENCRYPT);
		this.encryptButton.addActionListener(this);
		processPanel.add(this.encryptButton);

		new Dimension(5, 0);

		// Add decrypt button.
		this.decryptButton = new JButton(DECRYPT);
		this.decryptButton.setMnemonic(KeyEvent.VK_D);
		this.decryptButton.setActionCommand(DECRYPT);
		this.decryptButton.addActionListener(this);
		processPanel.add(this.decryptButton);
		// ----------------------
		this.browseB = new JButton(BROWSE_FILE);
		this.browseB.setMnemonic(KeyEvent.VK_O);
		this.browseB.setActionCommand(BROWSE_FILE);
		this.browseB.addActionListener(this);
		infileTf.setEditable(false);

		JPanel p111 = new JPanel();
		p111.setLayout(new FlowLayout(FlowLayout.CENTER));
		String res111 = "Input file: ";
		JLabel label111 = new JLabel(res111);
		p111.add(label111);
		p111.add(infileTf);
		p111.add(browseB);
		res111 = "Output file name: ";
		JLabel label112 = new JLabel(res111);
		p111.add(label112);
		p111.add(outfileTf);
		p111.setBackground(fundal);
		// pan.add(p111);
		// -------------------------
		// // Add encrypt file button.
		this.encryptFileButton = new JButton(ENCRYPT_FILE);
		this.encryptFileButton.setMnemonic(KeyEvent.VK_F);
		this.encryptFileButton.setActionCommand(ENCRYPT_FILE);
		this.encryptFileButton.addActionListener(this);
		this.decryptFileButton = new JButton(DECRYPT_FILE);
		this.decryptFileButton.setMnemonic(KeyEvent.VK_Y);
		this.decryptFileButton.setActionCommand(DECRYPT_FILE);
		this.decryptFileButton.addActionListener(this);
		processPanel.setBackground(fundal);
		// Add the process panel.
		pan.add(processPanel);
		pan.add(p111);
		JPanel p2 = new JPanel();
		p2.setLayout(new FlowLayout(FlowLayout.CENTER));
		p2.add(encryptFileButton);
		p2.add(decryptFileButton);
		p2.setBackground(fundal);
		pan.add(Box.createRigidArea(spacer));
		pan.add(p2);
		pan.setBackground(fundal);
		return pan;
	}

	// =================================
	/**
	 * Initialize the cipher object.
	 */
	private void initializeCipher() {
		try {
			// generate the IV for AES encryption with
			// padding==========================
			encryptionIV = randomBytes(ENCRYPTION_IV_LENGTH);// AES blocksize =
																// 16 bytes!
			// further:
			// write the header to the output stream.
			// ===========================================================================
			// create RSA public key cipher
			pkCipher = Cipher.getInstance(ENCRYPTIONKEY_ALGORITHM);// "RSA");
			// create AES shared key cipher
			aesCipher = Cipher.getInstance(ENCRYPTION_ALGORITHM);// "AES");
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	/**
	 * Generates random bytes of fixed length.
	 * 
	 * @param length
	 *            of byte array
	 * @return the byte array
	 */
	private byte[] randomBytes(int length) {
		if (random == null) {
			random = new SecureRandom();
		}

		byte[] bytes = new byte[length];
		random.nextBytes(bytes);

		return bytes;
	}

	/**
	 * Read bytes of an input stream.
	 * 
	 * @param length
	 *            of bytes array
	 * @param in
	 *            the input stream
	 * @return the bytes array
	 * @throws Exception
	 *             can throw this exception
	 */
	private byte[] readBytes(int length, InputStream in) throws Exception {
		byte[] bytes = new byte[length];
		int read = in.read(bytes);

		if (read != length) {
			throw new Exception("expected length != actual length");
		}

		return bytes;
	}

	/**
	 * Creates a new AES key
	 * 
	 * @throws NoSuchAlgorithmException
	 *             can throw this exception
	 */
	public void makeKey() throws NoSuchAlgorithmException {
		KeyGenerator kgen = KeyGenerator.getInstance(ENCRYPTION_KEY);// "AES");
		kgen.init(AES_Key_Size);
		SecretKey key = kgen.generateKey();
		aesKey = key.getEncoded();
		aeskeySpec = new SecretKeySpec(aesKey, ENCRYPTION_KEY);// "AES");
	}

	/**
	 * Decrypts an AES key from a file using an RSA private key
	 * 
	 * @param in
	 *            the input file
	 * @param privateKeyFile
	 *            the private key file
	 * @throws GeneralSecurityException
	 *             can throw this exception
	 * @throws IOException
	 *             can throw this exception
	 */
	public void loadKey(File in, File privateKeyFile) throws GeneralSecurityException, IOException {
		// read private key to be used to decrypt the AES key
		byte[] encodedKey = new byte[(int) privateKeyFile.length()];
		// new FileInputStream(privateKeyFile).read(encodedKey);
		FileInputStream fis = new FileInputStream(privateKeyFile);
		fis.read(encodedKey);

		// create private key
		PKCS8EncodedKeySpec privateKeySpec = new PKCS8EncodedKeySpec(encodedKey);
		KeyFactory kf = KeyFactory.getInstance(ENCRYPTIONKEY_ALGORITHM);// "RSA");
		PrivateKey pk = kf.generatePrivate(privateKeySpec);

		// read AES key
		pkCipher.init(Cipher.DECRYPT_MODE, pk);
		aesKey = new byte[AES_Key_Size / 8];
		CipherInputStream is = new CipherInputStream(new FileInputStream(in), pkCipher);
		is.read(aesKey);
		aeskeySpec = new SecretKeySpec(aesKey, ENCRYPTION_KEY);// "AES");

		fis.close();
		is.close();
	}

	/**
	 * Encrypts the AES key to a file using an RSA public key
	 * 
	 * @param out
	 *            the output file
	 * @param publicKeyFile
	 *            the public key file
	 * @throws IOException
	 *             can throw this exception
	 * @throws GeneralSecurityException
	 *             can throw this exception
	 */
	public void saveKey(File out, File publicKeyFile) throws IOException, GeneralSecurityException {
		// read public key to be used to encrypt the AES key
		byte[] encodedKey = new byte[(int) publicKeyFile.length()];
		// new FileInputStream(publicKeyFile).read(encodedKey);
		FileInputStream fis = new FileInputStream(publicKeyFile);
		fis.read(encodedKey);

		// create public key
		X509EncodedKeySpec publicKeySpec = new X509EncodedKeySpec(encodedKey);
		KeyFactory kf = KeyFactory.getInstance(ENCRYPTIONKEY_ALGORITHM);// "RSA");
		PublicKey pk = kf.generatePublic(publicKeySpec);

		// write AES key
		pkCipher.init(Cipher.ENCRYPT_MODE, pk);
		CipherOutputStream os = new CipherOutputStream(new FileOutputStream(out), pkCipher);
		os.write(aesKey);
		os.close();
		fis.close();
	}

	/**
	 * Encrypts and then copies the contents of a given file.
	 * 
	 * @param in
	 *            the input file
	 * @param out
	 *            the output file
	 * @throws Exception
	 *             can throw this exception
	 */
	public void encrypt(File in, File out) throws Exception {
		
		aesCipher.init(Cipher.ENCRYPT_MODE, aeskeySpec, new IvParameterSpec(encryptionIV));

		FileOutputStream fos = new FileOutputStream(out);
		fos.write(encryptionIV);

		FileInputStream is = new FileInputStream(in);
		CipherOutputStream os = new CipherOutputStream(fos, aesCipher);

		copy(is, os);

		os.close();
	}

	// ===========
	/**
	 * Decodes base64 text.
	 * 
	 * @param base64Text
	 *            the text
	 * @return the decoded byte array.
	 * @throws IOException
	 *             can throw this exception
	 */
	private byte[] getBytesFromBase64(String base64Text) throws IOException {
		return Base64Coder.decode(base64Text);
	}

	/**
	 * Encodes the array of bytes into base 64 text.
	 * 
	 * @param bytes
	 *            the bytes
	 * @return the encoded text.
	 */
	private String getBase64Text(byte[] bytes) {
		return new String(Base64Coder.encode(bytes));
	}

	// ============

	/**
	 * Encrypts the string into an array of bytes.
	 * 
	 * @param plainText
	 *            the input text.
	 * @return the encrypted array of bytes.
	 * @throws Exception
	 *             can throw this exception
	 */
	public byte[] encryptString(String plainText) throws Exception {
		
		aesCipher.init(Cipher.ENCRYPT_MODE, aeskeySpec, new IvParameterSpec(encryptionIV));
		
		byte[] result = aesCipher.doFinal(plainText.getBytes("UTF-8"));
		return result;
	}

	/**
	 * Decrypts and then copies the contents of a given file.
	 * 
	 * @param in
	 *            the input file
	 * @param out
	 *            the output file
	 * @throws Exception
	 *             can throw this exception
	 */
	public void decrypt(File in, File out) throws Exception {
		
		FileInputStream fis = new FileInputStream(in);
		byte[] encryptionIV = readBytes(ENCRYPTION_IV_LENGTH, fis);
		this.encryptionIV = encryptionIV;
		// ===========
		aesCipher.init(Cipher.DECRYPT_MODE, aeskeySpec, new IvParameterSpec(encryptionIV));

		CipherInputStream is = new CipherInputStream(fis, aesCipher);
		FileOutputStream os = new FileOutputStream(out);

		copy(is, os);

		is.close();
		os.close();
	}

	/**
	 * Decrypts the array of bytes
	 * 
	 * @param cipherText
	 *            the input array of bytes.
	 * @return the decrypted string
	 * @throws Exception
	 *             can throw this exception
	 */
	public String decryptBytes(byte[] cipherText) throws Exception {
		
		aesCipher.init(Cipher.DECRYPT_MODE, aeskeySpec, new IvParameterSpec(encryptionIV));

		byte[] original = aesCipher.doFinal(cipherText);
		return new String(original);// raw no "UTF-8");
		
	}

	/**
	 * Copies a stream.
	 * 
	 * @param is
	 *            the input stream
	 * @param os
	 *            the output stream
	 * @throws IOException
	 *             can throw this exception
	 */
	private void copy(InputStream is, OutputStream os) throws IOException {
		int i;
		byte[] b = new byte[1024];
		while ((i = is.read(b)) != -1) {
			os.write(b, 0, i);
		}
	}

	/**
	 * Get the file extension.
	 * 
	 * @param file
	 *            the file
	 * @return the file extension.
	 */
	private String getFileExtension(File file) {
		String name = file.getName();
		int lastIndexOf = name.lastIndexOf(".");
		if (lastIndexOf == -1) {
			return ""; // empty extension
		}
		return name.substring(lastIndexOf);
	}

	// ================================================PubPrivKey_forAESKEY
	/**
	 * Generates the public-private key pair.
	 */
	private void generateKeyPairIntoFiles() {
		String currentDir = System.getProperty("user.dir");
		String file_sep = System.getProperty("file.separator");
		String path = currentDir + file_sep;

		System.out.println("Path: " + path);

		try {
			KeyPairGenerator keyGen = KeyPairGenerator.getInstance(ENCRYPTIONKEY_ALGORITHM);// "RSA");
			keyGen.initialize(2048);// 2048 bit key length
			KeyPair generatedKeyPair = keyGen.genKeyPair();

			System.out.println("Generated Key Pair");
			dumpKeyPair(generatedKeyPair);

			saveKeyPair(path, generatedKeyPair);// save the key pair to files

			KeyPair loadedKeyPair = loadKeyPair(path, ENCRYPTIONKEY_ALGORITHM);// "RSA");
			System.out.println("Loaded Key Pair");
			dumpKeyPair(loadedKeyPair);// ok, it is verified is the same as
										// generated one!

		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			return;
		}
	}

	/**
	 * Get the hex representation for an input array of bytes.
	 * 
	 * @param b
	 *            the array of bytes
	 * @return the hex representation
	 */
	private String getHexString(byte[] b) {
		String result = "";
		for (int i = 0; i < b.length; i++) {
			result += Integer.toString((b[i] & 0xff) + 0x100, 16).substring(1);
		}
		return result;
	}

	/**
	 * Display the public-private key pairs in hex.
	 * 
	 * @param keyPair
	 *            the key pair
	 */
	private void dumpKeyPair(KeyPair keyPair) {
		PublicKey pub = keyPair.getPublic();
		System.out.println("Public Key: " + getHexString(pub.getEncoded()));

		PrivateKey priv = keyPair.getPrivate();
		System.out.println("Private Key: " + getHexString(priv.getEncoded()));
	}

	/**
	 * Save the key pair.
	 * 
	 * @param path
	 *            the path, i.e. where to save.
	 * @param keyPair
	 *            the key pair.
	 * @throws IOException
	 *             can throw this exception.
	 */
	public void saveKeyPair(String path, KeyPair keyPair) throws IOException {
		PrivateKey privateKey = keyPair.getPrivate();
		PublicKey publicKey = keyPair.getPublic();

		// Store Public Key.
		X509EncodedKeySpec x509EncodedKeySpec = new X509EncodedKeySpec(publicKey.getEncoded());
		FileOutputStream fos = new FileOutputStream(path + this.publicKey);// "public.key");
		fos.write(x509EncodedKeySpec.getEncoded());
		fos.close();

		// Store Private Key.
		PKCS8EncodedKeySpec pkcs8EncodedKeySpec = new PKCS8EncodedKeySpec(privateKey.getEncoded());
		fos = new FileOutputStream(path + this.privateKey);// "private.key");
		fos.write(pkcs8EncodedKeySpec.getEncoded());
		fos.close();
	}

	/**
	 * Load the key pair.
	 * 
	 * @param path
	 *            the path
	 * @param algorithm
	 *            the algorithm used for key pair generator.
	 * @return the key pair.
	 * @throws IOException
	 *             can throw this exception.
	 * @throws NoSuchAlgorithmException
	 *             can throw this exception.
	 * @throws InvalidKeySpecException
	 *             can throw this exception.
	 */
	public KeyPair loadKeyPair(String path, String algorithm)
			throws IOException, NoSuchAlgorithmException, InvalidKeySpecException {
		// Read Public Key.
		File filePublicKey = new File(path + this.publicKey);// "public.key");
		FileInputStream fis = new FileInputStream(path + this.publicKey);// "public.key");
		byte[] encodedPublicKey = new byte[(int) filePublicKey.length()];
		fis.read(encodedPublicKey);
		fis.close();

		// Read Private Key.
		File filePrivateKey = new File(path + this.privateKey);// "private.key");
		fis = new FileInputStream(path + this.privateKey);// "private.key");
		byte[] encodedPrivateKey = new byte[(int) filePrivateKey.length()];
		fis.read(encodedPrivateKey);
		fis.close();

		// Generate KeyPair.
		KeyFactory keyFactory = KeyFactory.getInstance(algorithm);
		X509EncodedKeySpec publicKeySpec = new X509EncodedKeySpec(encodedPublicKey);
		PublicKey publicKey = keyFactory.generatePublic(publicKeySpec);

		PKCS8EncodedKeySpec privateKeySpec = new PKCS8EncodedKeySpec(encodedPrivateKey);
		PrivateKey privateKey = keyFactory.generatePrivate(privateKeySpec);

		return new KeyPair(publicKey, privateKey);
	}

	// =========================================================PubPrivKey_forAES
	public static void main(String[] args) {
		LookAndFeel.loadLookAndFeel();
		//JCrypto2.textAreaBkgColor=Color.BLACK;
		//JCrypto2.textAreaForeColor=Color.YELLOW;
		new JCrypto2();
	}
}
