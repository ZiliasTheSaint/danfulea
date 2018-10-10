package danfulea.crypto.walletapp;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.util.ResourceBundle;

import javax.swing.BorderFactory;
import javax.swing.BoxLayout;
import javax.swing.JButton;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JPasswordField;
import javax.swing.border.Border;

import danfulea.utils.FrameUtilities;

/**
 * Password input class
 * 
 * @author Dan Fulea, 29 MAR. 2010
 */

public class PasswordFrame extends JFrame implements ActionListener {

	private static final long serialVersionUID = 1L;
	protected MainForm mf;
	private static final String BASE_RESOURCE_CLASS = "danfulea.crypto.walletapp.resources.JWalletResources";
	protected ResourceBundle resources;
	private static final Dimension PREFERRED_SIZE = new Dimension(400, 160);
	public static final int SAVE_MODE = 0;// used for saving=>a password check
											// is required
	public static final int OPEN_MODE = 1;// used for opening=>a password input
											// is required
	private JPasswordField pswPf = new JPasswordField();
	private JPasswordField psw_checkPf = new JPasswordField();
	private int mode = PasswordFrame.SAVE_MODE;
	private static final String OK_COMMAND = "OK";
	private static final String CANCEL_COMMAND = "CANCEL";
	private static final String ENTER_COMMAND = "ENTER";

	/**
	 * Constructor. Password window is connected to main window. <br>
	 * mode=0 - used for saving file (a password check is required). <br>
	 * mode=1 - used for opening a saved file (a password input is required).
	 * 
	 * @param mode
	 *            either to save ot to open
	 * @param mf
	 *            link to MainForm
	 */
	public PasswordFrame(int mode, MainForm mf) {
		this.mode = mode;
		this.mf = mf;
		resources = ResourceBundle.getBundle(BASE_RESOURCE_CLASS);

		createGUI();
		setTitle(resources.getString("passwordFrame.title"));
		setDefaultLookAndFeelDecorated(true);
		FrameUtilities.createImageIcon(this.resources.getString("form.icon.url"), this);
		FrameUtilities.centerFrameOnScreen(this);
		setVisible(true);

		mf.setEnabled(false);
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
	 * GUI creation.
	 */
	private void createGUI() {
		Dimension dim = new Dimension(150, 20);
		pswPf.setPreferredSize(dim);
		psw_checkPf.setPreferredSize(dim);
		// button panel
		JButton button = null;
		JPanel buttonPane = new JPanel();
		buttonPane.setLayout(new FlowLayout(FlowLayout.CENTER, 5, 1));

		String buttonName = resources.getString("passwordFrame.okB");
		String buttonToolTip = "";
		String buttonIconName = resources.getString("ok.icon");
		button = FrameUtilities.makeButton(buttonIconName, OK_COMMAND, buttonToolTip, buttonName, this, this);
		Character mnemonic = (Character) resources.getObject("passwordFrame.okB.mnemonic");
		button.setMnemonic(mnemonic.charValue());
		buttonPane.add(button);

		buttonName = resources.getString("passwordFrame.cancelB");
		buttonToolTip = "";
		buttonIconName = resources.getString("cancel.icon");
		button = FrameUtilities.makeButton(buttonIconName, CANCEL_COMMAND, buttonToolTip, buttonName, this, this);
		mnemonic = (Character) resources.getObject("passwordFrame.cancelB.mnemonic");
		button.setMnemonic(mnemonic.charValue());
		buttonPane.add(button);

		Border blackline;
		blackline = BorderFactory.createLineBorder(Color.black);
		buttonPane.setBorder(blackline);
		// end button panel

		// text panel
		JPanel textPane = new JPanel();
		textPane.setLayout(new BoxLayout(textPane, BoxLayout.Y_AXIS));

		if (mode == PasswordFrame.SAVE_MODE) {
			pswPf.setActionCommand(ENTER_COMMAND);
			pswPf.addActionListener(this);

			psw_checkPf.setActionCommand(OK_COMMAND);
			psw_checkPf.addActionListener(this);

			JLabel passL = new JLabel(resources.getString("passwordFrame.pass.new"));
			passL.setLabelFor(pswPf);

			JPanel pass1P = new JPanel();
			pass1P.setLayout(new FlowLayout(FlowLayout.CENTER));
			pass1P.add(passL);
			pass1P.add(pswPf);

			JLabel pass2L = new JLabel(resources.getString("passwordFrame.pass.check"));
			passL.setLabelFor(pswPf);

			JPanel pass2P = new JPanel();
			pass2P.setLayout(new FlowLayout(FlowLayout.CENTER));
			pass2P.add(pass2L);
			pass2P.add(psw_checkPf);

			textPane.add(pass1P);
			textPane.add(pass2P);
		} else if (mode == PasswordFrame.OPEN_MODE) {
			pswPf.setActionCommand(OK_COMMAND);
			pswPf.addActionListener(this);

			JLabel passL = new JLabel(resources.getString("passwordFrame.pass"));
			passL.setLabelFor(pswPf);

			JPanel pass1P = new JPanel();
			pass1P.setLayout(new FlowLayout(FlowLayout.CENTER));
			pass1P.add(passL);
			pass1P.add(pswPf);

			textPane.add(pass1P);
		}

		// end text panel
		JPanel content = new JPanel();
		content.setLayout(new BorderLayout());
		content.add(textPane, BorderLayout.CENTER);
		content.add(buttonPane, BorderLayout.SOUTH);
		content.setOpaque(true);

		getContentPane().add(content);
		pack();
	}

	/**
	 * Setting up actions!
	 */
	public void actionPerformed(ActionEvent evt) {
		String command = evt.getActionCommand();

		if (command.equals(OK_COMMAND)) {
			ok();
		} else if (command.equals(CANCEL_COMMAND)) {
			attemptExit();
		} else if (command.equals(ENTER_COMMAND)) {
			focus();
		}

	}

	/**
	 * Focus text field!
	 */
	private void focus() {
		psw_checkPf.requestFocus();
	}

	/**
	 * OK button!
	 */
	private void ok() {
		if (mode == PasswordFrame.OPEN_MODE) {
			mf.password = new String(pswPf.getPassword());
			mf.attemptOpenFile();
			attemptExit();
		} else if (mode == PasswordFrame.SAVE_MODE) {
			String pass1 = new String(pswPf.getPassword());
			String pass2 = new String(psw_checkPf.getPassword());
			int len = pass1.length();
			if (len < 6) {
				JOptionPane.showMessageDialog(null, resources.getString("dialog.inputPass.wrong"));
			} else {
				if (pass2.equals(pass1)) {
					mf.password = pass1;
					mf.attemptSaveFile();
					attemptExit();
				} else {
					JOptionPane.showMessageDialog(null, resources.getString("dialog.inputPass.wrong.check"));
				}

			}

		}
	}

	/**
	 * Window close!
	 */
	private void attemptExit() {
		final MainForm mff = mf;
		mff.setEnabled(true);
		dispose();
	}

}