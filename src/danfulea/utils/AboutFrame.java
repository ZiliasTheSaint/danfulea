package danfulea.utils;

import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.util.ResourceBundle;

import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTextArea;

import danfulea.resources.DanfuleaResources;

/**
 * This is the standard About Frame template. Just pass-in the resources
 * containing the required information to match different projects. The required
 * resource information are:
 * <p>
 * "About.NAME" - the title for About Frame, usually About
 * </p>
 * <p>
 * "form.icon.url" - the image icon URL, e.g: /danfulea/resources/duke.png
 * </p>
 * <p>
 * "logo.icon.url" - the logo image, e.g: /danfulea/resources/globe.gif
 * </p>
 * <p>
 * "Application.NAME" - the title for the application project e.g. MyApplication
 * </p>
 * <p>
 * "Author" - this should be Author:
 * </p>
 * <p>
 * "Author.name" - this is going to be you, the developer, e.g.: John Doe,
 * johndoe@johndoe.org
 * </p>
 * <p>
 * "Version" - this should be Version:
 * </p>
 * <p>
 * "Version.name" - the version of the project, e.g.: MyApplication 1.0
 * </p>
 * <p>
 * "License.name" - whatever the license type is, e.g: BSD License
 * </p>
 * <p>
 * "License" - whatever license you provide, e.g. copy and paste the BSD license
 * or Apache License 2.0 etc...just google for it!
 * </p>
 * 
 * @author Dan Fulea, 05 AUG. 2016
 *
 */
public class AboutFrame extends JFrame {
	
	private static final long serialVersionUID = 1L;

	/**
	 * The resources to be used
	 */
	private ResourceBundle resources;

	/**
	 * The AboutFrame constructor
	 * 
	 * @param resources
	 *            the resource bundle
	 */
	public AboutFrame(ResourceBundle resources) {
		this.resources = resources;

		this.setTitle(resources.getString("About.NAME"));
		this.setResizable(false);

		createGUI();

		setDefaultLookAndFeelDecorated(true);
		FrameUtilities.createImageIcon(this.resources.getString("form.icon.url"), this);

		FrameUtilities.centerFrameOnScreen(this);

		setVisible(true);

		setDefaultCloseOperation(JFrame.DO_NOTHING_ON_CLOSE);

		addWindowListener(new WindowAdapter() {
			public void windowClosing(WindowEvent e) {
				attemptExit();
			}
		});
	}

	/**
	 * Exit method
	 */
	private void attemptExit() {
		dispose();
	}

	/**
	 * This method is called from within the constructor to initialize the form.
	 */
	private void createGUI() {
		JLabel jLabel1 = new JLabel();
		JLabel jLabel2 = new JLabel();
		JLabel lbAuthor = new JLabel();
		JLabel lbVersion = new JLabel();
		JLabel jLabel7 = new JLabel();
		JPanel jPanel1 = new JPanel();
		JPanel jPanel2 = new JPanel();
		JPanel jPanel3 = new JPanel();
		JScrollPane jScrollPane1 = new JScrollPane();
		JTextArea textLicense = new JTextArea();

		jLabel1.setForeground(DanfuleaResources.foreColor);
		jLabel2.setForeground(DanfuleaResources.foreColor);
		lbAuthor.setForeground(DanfuleaResources.foreColor);
		lbVersion.setForeground(DanfuleaResources.foreColor);
		jLabel7.setForeground(DanfuleaResources.foreColor);
		textLicense.setBackground(DanfuleaResources.textAreaBkgColor);
		textLicense.setForeground(DanfuleaResources.textAreaForeColor);

		jPanel1.setLayout(new java.awt.BorderLayout());

		jLabel1.setHorizontalAlignment(javax.swing.SwingConstants.LEFT);
		jLabel1.setIcon(FrameUtilities.getImageIcon(this.resources.getString("logo.icon.url"), this));
		jLabel1.setText(this.resources.getString("Application.NAME"));
		jLabel1.setVerticalAlignment(javax.swing.SwingConstants.TOP);
		jPanel1.add(jLabel1, java.awt.BorderLayout.NORTH);
		jPanel1.setBackground(DanfuleaResources.bkgColor);

		jPanel3.setLayout(new java.awt.GridLayout(3, 2, 0, 4));

		jLabel2.setText(this.resources.getString("Author"));
		jLabel2.setHorizontalAlignment(javax.swing.SwingConstants.LEFT);

		jPanel3.add(jLabel2);

		lbAuthor.setText(this.resources.getString("Author.name"));

		jPanel3.add(lbAuthor);

		jLabel7.setText(this.resources.getString("Version"));
		jLabel7.setHorizontalAlignment(javax.swing.SwingConstants.LEFT);

		jPanel3.add(jLabel7);

		lbVersion.setText(this.resources.getString("Version.name"));

		jPanel3.add(lbVersion);
		jPanel3.setBackground(DanfuleaResources.bkgColor);

		jPanel1.add(jPanel3, java.awt.BorderLayout.SOUTH);

		getContentPane().add(jPanel1, java.awt.BorderLayout.NORTH);

		jPanel2.setLayout(new java.awt.BorderLayout());

		jScrollPane1.setBorder(new javax.swing.border.TitledBorder(
				new javax.swing.border.LineBorder(new java.awt.Color(0, 51, 255), 1, true), resources.getString("License.name"),//"BSD Licence",
				javax.swing.border.TitledBorder.CENTER, javax.swing.border.TitledBorder.TOP));
		jScrollPane1.setHorizontalScrollBarPolicy(javax.swing.ScrollPaneConstants.HORIZONTAL_SCROLLBAR_ALWAYS);
		jScrollPane1.setVerticalScrollBarPolicy(javax.swing.ScrollPaneConstants.VERTICAL_SCROLLBAR_ALWAYS);
		jScrollPane1.setAutoscrolls(true);
		textLicense.setColumns(1);
		textLicense.setEditable(false);

		textLicense.setLineWrap(true);
		textLicense.setRows(10);
		textLicense.setText(this.resources.getString("License"));
		textLicense.setWrapStyleWord(true);
		jScrollPane1.setViewportView(textLicense);

		jPanel2.add(jScrollPane1, java.awt.BorderLayout.CENTER);
		jPanel2.setBackground(DanfuleaResources.bkgColor);

		getContentPane().add(jPanel2, java.awt.BorderLayout.CENTER);
		pack();
	}

}
