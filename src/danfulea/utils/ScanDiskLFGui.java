package danfulea.utils;


import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.Window;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.io.File;
import java.io.FileWriter;
import java.util.ResourceBundle;

import javax.swing.BorderFactory;
import javax.swing.BoxLayout;
import javax.swing.ButtonGroup;
import javax.swing.JButton;
import javax.swing.JCheckBox;
import javax.swing.JComboBox;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JRadioButton;
import javax.swing.JScrollPane;
import javax.swing.JSlider;
import javax.swing.JTextArea;
import javax.swing.JTextField;
import javax.swing.JToolBar;
import javax.swing.SwingUtilities;
import javax.swing.UIManager;
import javax.swing.UIManager.LookAndFeelInfo;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import danfulea.utils.gadgets.ChangeColorPanel;
import danfulea.utils.gadgets.ElegantClockPanel;

/**
 * Scan disk and look and feel GUI!<br>
 * Also display some gadgets such as an elegant clock and color chooser! 
 * 
 * @author Dan Fulea, 04 JUL. 2011
 * 
 */
public class ScanDiskLFGui extends JFrame implements ActionListener, Runnable,
ItemListener, ChangeListener, MessageRetriever{
	
	private static final long serialVersionUID = 1L;
	private static final String BASE_RESOURCE_CLASS = "danfulea.resources.DanfuleaResources";
	private ResourceBundle resources;

	private JLabel statusL = new JLabel();

	private Thread statusTh = null;// status display thread!
	//the volatile variable will have only one main copy which will be updated by 
	//different threads and update made by one thread to the volatile variable will 
	//immediately reflect to the other Thread. The volatile keyword is used to ensure
	//prompt communication between threads. Since we did not use thread communications
	//but only starting different independent threads, no need to make it volatile.
	//anyway, if commmunication is desire as in Producer-Consumer template (see test package), it is best
	//to use synchronization then volatile. Use volatile only when 1 threads update the
	//volatile variable and all other threads just read it!! So it is useful for small things.
	private int delay = 100;
	private int frameNumber = -1;
	private String statusRunS = "";
	private String command = "";
	private Thread scanTh = null;// scan thread
	private Thread timer = null; // The thread that displays clock
	private boolean stopAnim=true;
	
	private int textAreaRowCount=0;
	private static final int maxLines=1000;
	private String initS="";
	
	private static final String RUN_COMMAND = "RUN";
	private static final String KILL_COMMAND = "KILL";
	private static final String SETLF_COMMAND = "LF";
	private static final String EXIT_COMMAND = "EXIT";

	public static final Dimension PREFERRED_SIZE = new Dimension(990, 650);
	private static final Dimension textAreaDimension = new Dimension(300, 330);
	private static final Dimension clockPanelSize = new Dimension(350, 180);
	private static final Dimension sizePan = new Dimension(350, 100);
	private JRadioButton scanRb, lookRb;
	private JCheckBox caseCb;
	private JTextField targetTf = new JTextField(5);
	private JTextArea textArea;

	//--------------------------------
	private	ElegantClockPanel advCp;
	// --------------------------------
	// --------ChangeColor related
	int R = 0, G = 0, B = 0, A = 0;
	Color c;
	JPanel pv;
	
	JSlider sR, sG, sB, sA;
	@SuppressWarnings("rawtypes")
	JComboBox selColorCb, lfCb;
	JLabel lR, lG, lB, lA;
	private static final Dimension sizeTxtCb = new Dimension(150, 21);
	private JLabel tf_hc = new JLabel("");
	private ChangeColorPanel cc;
	public static final int iBrightGold = 0xD9D919;
	public static final int iCopper = 0xB87333;
	public static final int iCoral = 0xFF7F00;
	public static final int iDustyRose = 0x856363;
	public static final int iForestGreen = 0x238E23;
	public static final int iKhaki = 0x9F9F5F;
	public static final int iMidnightBlue = 0x2F2F4F;
	public static final int iNeonPink = 0xFF6EC7;
	public static final int iSalmon = 0x6F4242;
	public static final int iTan = 0xDB9370;
	
	public static final String[] colors = { "AliceBlue", "AntiqueWhite",
			"Aquamarine", "Azure", "Beige", "Bisque", "Black",
			"BlanchedAlmond", "Blue", "BlueViolet", "BrightGold", "Brown",
			"Burlywood", "CadetBlue", "Chartreuse", "Chocolate", "Copper",
			"Coral", "CornflowerBlue", "Cornsilk", "Cyan", "DarkGoldenrod",
			"DarkGreen", "DarkGrey", "DarkKhaki", "DarkOliveGreen",
			"DarkOrange", "DarkOrchid", "DarkSalmon", "DarkSeaGreen",
			"DarkSlateBlue", "DarkSlateGray", "DarkTurquoise", "DarkViolet",
			"DeepPink", "DeepSkyBlue", "DimGray", "DodgerBlue", "DustyRose",
			"Firebrick", "FloralWhite", "ForestGreen", "Gainsboro",
			"GhostWhite", "Gold", "Goldenrod", "Gray", "Green", "GreenYellow",
			"Honeydew", "HotPink", "IndianRed", "Ivory", "Khaki", "Lavender",
			"LavenderBlush", "LawnGreen", "LemonChiffon", "LightBlue",
			"LightCoral", "LightCyan", "LightGoldenrod",
			"LightGoldenrodYellow", "LightGray", "LightPink", "LightSalmon",
			"LightSeaGreen", "LightSkyBlue", "LightSlateBlue",
			"LightSlateGray", "LightSteelBlue", "LightYellow", "LimeGreen",
			"Linen", "Magenta", "Maroon", "MediumAquamarine", "MediumBlue",
			"MediumOrchid", "MediumPurple", "MediumSeaGreen",
			"MediumSlateBlue", "MediumSpringGreen", "MediumTurquoise",
			"MediumVioletRed", "MidnightBlue", "MintCream", "MistyRose",
			"Moccasin", "NavajoWhite", "Navy", "NeonPink", "OldLace",
			"OliveDrab", "Orange", "OrangeRed", "Orchid", "PaleGoldenrod",
			"PaleGreen", "PaleTurquoise", "PaleVioletRed", "PapayaWhip",
			"PeachPuff", "Peru", "Pink", "Plum", "PowderBlue", "PreDawn",
			"Purple", "Red", "RosyBrown", "RoyalBlue", "SaddleBrown", "Salmon",
			"SandyBrown", "SeaGreen", "Seashell", "Sienna", "SkyBlue",
			"SlateBlue", "SlateGray", "Snow", "SpringGreen", "SteelBlue",
			"Tan", "Thistle", "Tomato", "Turquoise", "Violet", "Wheat",
			"White", "WhiteSmoke", "Yellow", "YellowGreen" };
	// color Panel
	private String filename = "JavaLookAndFeelLoader.laf";
	static String quaquaLF = "Quaqua";
	static String infoLF = "Info";

	//static String jtatooLF = "JTatoo";//this is licensed under GPL=>Not good
	//static String jtatooClassName = "com.jtattoo.plaf.smart.SmartLookAndFeel";

	static String seaGlass = "Sea Glass";
	static String seaGlassClassName = "com.seaglasslookandfeel.SeaGlassLookAndFeel";

	static String kunststoff = "Kunststoff";
	static String kunststoffClassName = "com.incors.plaf.kunststoff.KunststoffLookAndFeel";
	// ----------------------------------------------------------------------------------
	static String currentLF = UIManager.getLookAndFeel().getClass().getName();
	// -----------------------------------------------------------------------------------
	static String tonic = "Tonic";
	static String tonicClassName = "com.digitprop.tonic.TonicLookAndFeel";
	static String pgs = "Pgs";
	static String pgsClassName = "com.pagosoft.plaf.PgsLookAndFeel";
	static String office = "Office";
	static String officeClassName = "org.fife.plaf.Office2003.Office2003LookAndFeel";
	static String next = "Next";
	static String nextClassName = "nextlf.plaf.NextLookAndFeel";
	static String napkin = "Napkin";
	static String napkinClassName = "napkin.NapkinLookAndFeel";
	static String liquid = "Liquid";
	static String liquidClassName = "com.birosoft.liquid.LiquidLookAndFeel";

	Window parent = null;

	public static final String[] lfs = { 
			"CDE/Motif",
			//"GTK+",//only under Linux but similar to system
			//"Info",//useless
			"Java(tm)", //"JTatoo",//GPL License
			"Kunststoff", "Liquid",
			//"Metal",//only under Linux but similar to Java
			// "Motif",
			"Napkin",
			// "Next",
			"Nimbus", 
			//"Office",//not supported on Linux 
			"Pgs", 
			//"Quaqua",//lame and work only with mac
			//"Sea Glass",//glitches on Linux!!
			"System", 
			"Tonic",
			// "Windows_LF",
			//"Windows",//only under windows..use System here!
			//"Windows Classic" //only under Windows
			};
	
	/**
	 * Constructor, here the parent Window is set!
	 * @param frame
	 * 				the parent Window
	 */
	public ScanDiskLFGui(Window frame) {
		this();
		this.parent = frame;
		// final Window parent = frame;
		addWindowListener(new WindowAdapter() {
			public void windowClosing(WindowEvent e) {
				attemptExit();
			}
		});
	}

	/**
	 * Constructor with no Window, parent. Standalone test application
	 */
	public ScanDiskLFGui() {
		super("ScanDisk");//sets the title

		resources = ResourceBundle.getBundle(BASE_RESOURCE_CLASS);
		ScanDisk.mr = this;

		scanRb = new JRadioButton("Scan computer and print folder structure");
		scanRb.setToolTipText("The results are exported into ScanResult.txt file inside the application folder");

		lookRb = new JRadioButton(
				"Look for a file or folder into the computer -wildcard (*) usage enabled");
		caseCb = new JCheckBox("case sensitive");
		textArea = new JTextArea();
		textArea.setCaretPosition(0);
		textArea.setEditable(false);
		textArea.setText("");
		textArea.setWrapStyleWord(true);

		createGUI();
		setDefaultLookAndFeelDecorated(true);
		FrameUtilities.createImageIcon(
				this.resources.getString("form.icon.url"), this);
		FrameUtilities.centerFrameOnScreen(this);

		startClockThread();// /start the thread
		setVisible(true);
	}

	/**
	 * dispose this frame and exit.
	 */
	private void attemptExit() {
		parent.setVisible(true);
		stopThread();
		stopClockThread();
		try {
			UIManager.setLookAndFeel(currentLF);
		} catch (Exception exc) {
			// never
		}
		dispose();
	}

	/**
	 * GUI creation
	 */
	@SuppressWarnings({ "rawtypes", "unchecked" })
	private void createGUI() {
		Character mnemonic = null;
		JButton button = null;
		JLabel label = null;
		String buttonName = "";
		String buttonToolTip = "";
		String buttonIconName = "";

		ButtonGroup groupC = new ButtonGroup();
		groupC.add(scanRb);
		groupC.add(lookRb);
		scanRb.setSelected(true);

		JPanel textP = new JPanel();
		textP.setLayout(new FlowLayout(FlowLayout.CENTER, 10, 2));
		label = new JLabel("Looking for: ");
		textP.add(label);
		textP.add(targetTf);
		textP.add(caseCb);

		JPanel butP = new JPanel();
		butP.setLayout(new FlowLayout(FlowLayout.CENTER, 10, 2));
		buttonName = "Run";
		buttonToolTip = "";
		buttonIconName = null;//
		button = FrameUtilities.makeButton(buttonIconName, RUN_COMMAND,
				buttonToolTip, buttonName, this, this);
		mnemonic = new Character('R');
		button.setMnemonic(mnemonic.charValue());
		butP.add(button);
		buttonName = "Kill";
		buttonToolTip = "";
		buttonIconName = null;// "";
		button = FrameUtilities.makeButton(buttonIconName, KILL_COMMAND,
				buttonToolTip, buttonName, this, this);
		mnemonic = new Character('K');
		button.setMnemonic(mnemonic.charValue());
		butP.add(button);

		JPanel scanP = new JPanel();
		scanP.setLayout(new FlowLayout(FlowLayout.LEFT, 5, 1));
		scanP.add(scanRb, null);
		JPanel lookP = new JPanel();
		lookP.setLayout(new FlowLayout(FlowLayout.LEFT, 5, 1));
		lookP.add(lookRb, null);
		
		JPanel p1 = new JPanel();
		BoxLayout bl = new BoxLayout(p1, BoxLayout.Y_AXIS);
		p1.setLayout(bl);
		p1.add(scanP);
		p1.add(lookP);
		p1.add(textP, null);
		p1.add(butP, null);

		JPanel resultP = new JPanel(new BorderLayout());
		resultP.setPreferredSize(textAreaDimension);
		resultP.add(new JScrollPane(textArea), BorderLayout.CENTER);
		
		// =============			
		advCp = new ElegantClockPanel();
		advCp.setOpaque(true);

		JPanel westP = new JPanel();
		westP.setLayout(new BorderLayout());
		westP.add(advCp, BorderLayout.CENTER);
		advCp.setBackground(Color.WHITE);
		westP.setBackground(Color.WHITE);
		westP.setPreferredSize(clockPanelSize);
		// =============
		
		JPanel ps1, ps2, ps3, ps4, ps_n;
		selColorCb = new JComboBox(colors);
		String s = colors[6];// BLACK
		selColorCb.setSelectedItem((Object) s);
		selColorCb.setMaximumRowCount(10);
		selColorCb.setPreferredSize(sizeTxtCb);
		selColorCb.setOpaque(true);
		selColorCb.addItemListener(this);

		A = 255;
		sR = new JSlider(0, 255, R);
		sG = new JSlider(0, 255, G);
		sB = new JSlider(0, 255, B);
		sA = new JSlider(0, 255, A);
		sR.setOpaque(true);
		sG.setOpaque(true);
		sB.setOpaque(true);
		sA.setOpaque(true);
		sR.addChangeListener(this);
		sG.addChangeListener(this);
		sB.addChangeListener(this);
		sA.addChangeListener(this);

		lR = new JLabel("R= " + String.valueOf(R), JLabel.LEFT);// lR.setForeground(Color.white);
		lG = new JLabel("G= " + String.valueOf(G), JLabel.LEFT);// lG.setForeground(Color.white);
		lB = new JLabel("B= " + String.valueOf(B), JLabel.LEFT);// lB.setForeground(Color.white);
		lA = new JLabel("A= " + String.valueOf(A), JLabel.LEFT);// lA.setForeground(Color.white);

		ps1 = new JPanel();
		ps2 = new JPanel();
		ps3 = new JPanel();
		ps4 = new JPanel();
		ps_n = new JPanel();
		pv = new JPanel();
		JPanel pCb = new JPanel();

		ps1.setLayout(new FlowLayout(FlowLayout.CENTER));
		ps1.add("West", sR);
		sR.setBackground(Color.RED);
		ps1.add("Center", lR);
		ps1.setBackground(Color.RED);
		ps2.setLayout(new FlowLayout(FlowLayout.CENTER));
		ps2.add("West", sG);
		sG.setBackground(new Color(0, 100, 0));
		ps2.add("Center", lG);
		ps2.setBackground(new Color(0, 100, 0));
		ps3.setLayout(new FlowLayout(FlowLayout.CENTER));
		ps3.add("West", sB);
		sB.setBackground(Color.BLUE);
		ps3.add("Center", lB);
		ps3.setBackground(Color.BLUE);
		ps4.setLayout(new FlowLayout(FlowLayout.CENTER));
		ps4.add("West", sA);
		ps4.add("Center", lA);

		pCb.setLayout(new FlowLayout(FlowLayout.CENTER, 10, 2));
		JLabel lbl = new JLabel("Preset color: ");
		pCb.add(lbl);
		pCb.add(selColorCb);

		JPanel hexP = new JPanel();
		hexP.setLayout(new FlowLayout(FlowLayout.CENTER, 10, 2));
		lbl = new JLabel("Color HexCode: ");
		hexP.add(lbl);
		hexP.add(tf_hc);

		BoxLayout bl1 = new BoxLayout(ps_n, BoxLayout.Y_AXIS);
		ps_n.setLayout(bl1);
		ps_n.add(ps1);
		ps_n.add(ps2);
		ps_n.add(ps3);
		ps_n.add(ps4);
		ps_n.add(pCb);
		ps_n.add(hexP);
		ps_n.setBorder(FrameUtilities.getGroupBoxBorder("Color picker"));

		c = new Color(R, G, B, A);
		pv.setLayout(new BorderLayout());
		pv.setBorder(BorderFactory.createLineBorder(Color.black));
		pv.setBackground(c);
		pv.setOpaque(true);
		updateText(c);

		cc = new ChangeColorPanel();
		cc.setLayout(new BorderLayout());
		cc.setOpaque(true);
		pv.add(cc);
		pv.setOpaque(true);
		pv.setPreferredSize(sizePan);

		JPanel westAllP = new JPanel();
		BoxLayout bl11 = new BoxLayout(westAllP, BoxLayout.Y_AXIS);
		westAllP.setLayout(bl11);
		westAllP.add(westP);
		westAllP.add(ps_n);
		westAllP.add(pv);
		// ===========
		JPanel eastP = new JPanel(new BorderLayout());
		eastP.add(p1, BorderLayout.NORTH);
		eastP.add(resultP, BorderLayout.CENTER);

		JPanel mainP = new JPanel();
		mainP.setLayout(new FlowLayout(FlowLayout.CENTER, 10, 2));
		mainP.add(eastP);
		mainP.add(westAllP);

		JToolBar statusBar = new JToolBar();
		statusBar.setFloatable(false);
		initStatusBar(statusBar);

		// ===========
		lfCb = new JComboBox(lfs);
		String sss = lfs[1];// default
		lfCb.setSelectedItem((Object) sss);
		lfCb.setMaximumRowCount(10);
		lfCb.setPreferredSize(sizeTxtCb);
		lfCb.setOpaque(true);

		lbl = new JLabel("Choose look and feel: ");

		buttonName = "Set look and feel";
		buttonToolTip = "";
		buttonIconName = null;//
		button = FrameUtilities.makeButton(buttonIconName, SETLF_COMMAND,
				buttonToolTip, buttonName, this, this);
		mnemonic = new Character('S');
		button.setMnemonic(mnemonic.charValue());

		JPanel lfP = new JPanel();
		lfP.setLayout(new FlowLayout(FlowLayout.CENTER, 10, 2));
		lfP.add(lbl);
		lfP.add(lfCb);
		lfP.add(button);

		buttonName = "Exit";
		buttonToolTip = "";
		buttonIconName = null;//
		button = FrameUtilities.makeButton(buttonIconName, EXIT_COMMAND,
				buttonToolTip, buttonName, this, this);
		mnemonic = new Character('x');
		button.setMnemonic(mnemonic.charValue());
		lfP.add(button);
		// =========

		JPanel main2P = new JPanel();
		BoxLayout bl112 = new BoxLayout(main2P, BoxLayout.Y_AXIS);
		main2P.setLayout(bl112);
		main2P.add(mainP);
		main2P.add(lfP);

		JPanel content = new JPanel(new BorderLayout());
		content.add(statusBar, BorderLayout.PAGE_END);
		content.add(main2P, BorderLayout.CENTER);
		setContentPane(new JScrollPane(content));
		content.setOpaque(true); // content panes must be opaque
		pack();
	}

	/**
	 * Status bar initialization
	 * @param toolBar
	 * 			the status bar
	 */
	private void initStatusBar(JToolBar toolBar) {
		JPanel toolP = new JPanel();
		toolP.setLayout(new FlowLayout(FlowLayout.LEFT, 5, 1));

		toolP.add(statusL);
		toolBar.add(toolP);
		statusL.setText("Waiting for your action");
	}

	/**
	 * Overrides JFrame method in order to set the frame preferred dimension!
	 */
	public Dimension getPreferredSize() {
		return PREFERRED_SIZE;
	}

	/**
	 * Displays the chosen color in hex format
	 * @param color
	 * 			the color
	 */
	public final void updateText(Color color) {
		String s = Integer.toHexString(color.getRGB() & 0xffffff).toUpperCase();
		s = "000000" + s;
		s = s.substring(s.length() - 6);
		s = "#" + s;
		tf_hc.setText(s);
	}

	/**
	 * Standard actions are set here
	 */
	public void actionPerformed(ActionEvent arg0) {
		
		command = arg0.getActionCommand();
		if (command.equals(RUN_COMMAND)) {
			statusRunS = "Computing";
			startThread();
		} else if (command.equals(KILL_COMMAND)) {
			kill();
		} else if (command.equals(SETLF_COMMAND)) {
			setLf();
		} else if (command.equals(EXIT_COMMAND)) {
			attemptExit();
		}

	}

	/**
	 * Combobox actions are set here
	 */
	public void itemStateChanged(ItemEvent ie) {
		if (ie.getSource() == selColorCb) {
			changeColor();
		}
	}

	/**
	 * Slider actions are set here
	 */
	public void stateChanged(ChangeEvent e) {
		if (e.getSource() == sR) {
			R = sR.getValue();
		}
		if (e.getSource() == sG) {
			G = sG.getValue();
		}
		if (e.getSource() == sB) {
			B = sB.getValue();
		}
		if (e.getSource() == sA) {
			A = sA.getValue();
		}
		
		// invoke repaint for paintComponent method which is overriden				
		cc.setShowText(false);
		cc.setRGBA(R, G, B, A);//pass the RGBA
		cc.repaint();
		
		c = new Color(R,G,B,A);
		updateText(c);

		lR.setText("R= " + String.valueOf(R));
		lG.setText("G= " + String.valueOf(G));
		lB.setText("B= " + String.valueOf(B));
		lA.setText("A= " + String.valueOf(A));
	}

	/**
	 * Write the chosen look and feel into a file
	 * @param lnfName
	 * 				the file
	 */
	private void scheduleChanges(String lnfName) {
		// ---write file
		String fileSeparator = System.getProperty("file.separator");
		String curentDir = System.getProperty("user.dir");
		String filename = curentDir + fileSeparator + this.filename;
		try {
			File f = new File(filename);
			f.createNewFile();// Exceptie if nu poate crea(nu exista fisier<->Nu
								// sa determinat bine WinDir-ul!!
			FileWriter fw = new FileWriter(f);
			fw.write(lnfName);
			fw.close();
		} catch (Exception e) {

		}
		// ---
		// ----display scheduled changes
		
		textArea.selectAll();
		textArea.replaceSelection("");
		String text = "Restart the application in order to validate changes! ";
		textArea.append(text + " \n");
		textArea.append(lnfName);

	}

	/**
	 * Sets the look and feel
	 */
	private void setLf() {
		boolean foundB = false;
		try {
			String lnfName = (String) lfCb.getSelectedItem();

			JFrame.setDefaultLookAndFeelDecorated(true);

			for (LookAndFeelInfo info : UIManager.getInstalledLookAndFeels()) {
				if (lnfName.equals(info.getName())) {
					UIManager.setLookAndFeel(info.getClassName());
					foundB = true;
					break;
				}
			}

			if (!foundB) {
				if (lnfName.equals("System")) {
					foundB = true;
					UIManager.setLookAndFeel(UIManager
							.getSystemLookAndFeelClassName());
				} else if (lnfName.equals("Java(tm)")) {
					foundB = true;
					UIManager.setLookAndFeel(UIManager
							.getCrossPlatformLookAndFeelClassName());
				}
			}

			if (!foundB) {
				//if (lnfName.equals(jtatooLF)) {
					//lnfName = jtatooClassName;foundB = true;//licensed under GPL==not GOOD
				//} else if (lnfName.equals(seaGlass)) {
					//lnfName = seaGlassClassName;foundB = true;
				// else 
				if (lnfName.equals(kunststoff)) {
					lnfName = kunststoffClassName;foundB = true;
				} else if (lnfName.equals(tonic)) {
					lnfName = tonicClassName;foundB = true;
				} else if (lnfName.equals(pgs)) {
					lnfName = pgsClassName;foundB = true;
				//} else if (lnfName.equals(office)) {
					//lnfName = officeClassName;foundB = true;
				//} else if (lnfName.equals(next)) {
					//lnfName = nextClassName;foundB = true;
				} else if (lnfName.equals(napkin)) {
					lnfName = napkinClassName;foundB = true;
				} else if (lnfName.equals(liquid)) {
					lnfName = liquidClassName;foundB = true;
				}

				UIManager.setLookAndFeel(lnfName);
			}
			
			if(foundB){
				SwingUtilities.updateComponentTreeUI(this);// frame);

				//----!------------------
				scheduleChanges(lnfName);
				// --!!--------------------
				repaint();
			}//do nothing else!
		} catch (Exception exc) {
			exc.printStackTrace();
			try {
				UIManager.setLookAndFeel(currentLF);
			} catch (Exception e) {
				// never
			}
		}
	}

	/**
	 * The main method for scan/look for file
	 */
	private void executeRun() {
		//==============
		boolean isWinOS=false;
		boolean isLinux=false;
		String OS = System.getProperty("os.name").toLowerCase();
		//System.out.println(OS);//
		if ((OS.indexOf("windows") > -1) ||
			(OS.indexOf("nt") > -1)) {
			//windows
			isWinOS=true;//System.out.println("enter");
		}
		if (OS.indexOf("linux") > -1){
			isLinux=true;
		}
		//getting current root
		String root="";
		String parentF=System.getProperty("user.dir");
		//System.out.println(parentF+" ");//userDir!
		File f=new File(parentF);		
		while(true){
			parentF=f.getParent();			
			if (parentF==null){
				break;
			} else {
				f=new File(parentF);
				root=parentF;
			}
		}
		File rootFile=new File(root);
		//System.out.println(root+" "+OS);//!!!
		
		if (isLinux){
			File superF =null;
			//parent is not / but /home!!!
			parentF=System.getProperty("user.dir");
			//System.out.println(parentF+" ");//userDir!
			f=new File(parentF);
			//initialize!
			String superRoot=parentF;
			root=parentF;
			while(true){				
				parentF=f.getParent();
				//check superparent
				if (parentF!=null){//always
					superF=new File(parentF);
					superRoot=superF.getParent();
					if (superRoot==null){
						break;//here we exit the loop!
					} else{
						f=new File(superRoot);
						root=superRoot;						
					}
				} else {
					System.out.println("WE ARE IN / !!!");
					root="/home";//force..otherwise scan will fail anyway!
					break;
					//never enter here....
					//except we are in root /...this is not the case for application!
				}
				
			}
			rootFile=new File(root);
			//System.out.println(root+" "+OS);//!!!
		}
		//================
		textArea.selectAll();
		textArea.replaceSelection("");

		if (scanRb.isSelected()) {
			try {
				if (isWinOS){//System.out.println("enter");
					ScanDisk.createComputerScanResultFileInCurrentDir();
				}
				else//scan only current root drive not all computer!!!..gllitches in LINUX
					ScanDisk.createScanResultFileInCurrentDir(rootFile);
			} catch (Exception e) {
				statusL.setText("Unexpected error!");
				stopThread();
				e.printStackTrace();
				return;
			}
		} else if (lookRb.isSelected()) {
			boolean case_sensitiv = caseCb.isSelected();
			String target = targetTf.getText();
			if (targetTf.getText().compareTo("*") == 0)
				// desi asta da bine -aparent paradox:(cauza=comparare
				// lexicografica-v API)
				target = "*";
			try {
				if (isWinOS)
					ScanDisk.scanComputerForFile(target, case_sensitiv);
				else
					ScanDisk.scanForFile(rootFile, target, case_sensitiv);
			} catch (Exception e) {
				statusL.setText("Unexpected error!");
				stopThread();
				e.printStackTrace();
				return;
			}
		}

		stopThread();
		statusL.setText("Done!");
	}

	/**
	 * Stops the scan/look for files
	 */
	private void kill() {
		stopThread();// kill all threads
		statusL.setText("Process killed!");
	}

	/**
	 * Start the clock thread
	 */
	private void startClockThread() {
		if (timer == null) {
			timer = new Thread(this);
		}

		timer.start();
	}
	
	/**
	 * Stop the clock thread
	 */
	private void stopClockThread() {
		timer = null;
	}

	/**
	 * Start the scan thread
	 */
	private void startThread() {
		stopAnim=false;
		if (scanTh == null) {
			ScanDisk.stopB = false;
			scanTh = new Thread(this);
			scanTh.start();
		}

		if (statusTh == null) {
			statusTh = new Thread(this);
			statusTh.start();
		}

	}

	/**
	 * Stop the scan thread
	 */
	public void stopThread() {
		stopAnim=true;
		ScanDisk.stopB = true;
		statusTh = null;
		frameNumber = 0;
		// Stop thread
		scanTh = null;
	}
	
	/**
	 * Run method, required when working with threads
	 */
	public void run() {
		Thread.currentThread().setPriority(Thread.NORM_PRIORITY);

		long startTime = System.currentTimeMillis();
		Thread currentThread = Thread.currentThread();
		while (!stopAnim && currentThread == statusTh) {// if thread is status display
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

		if (currentThread == scanTh) {// if thread is the main
										// computation Thread!!
			executeRun();
		}

		while (timer == currentThread) {
			try {
				//Thread.currentThread();//just put it to sleep for 1000 ms Clocks shoud be updated every second!
				Thread.sleep(1000);
			} catch (InterruptedException e) {
			}

			advCp.repaint();
		}

	}
	
	/**
	 * Required method for MessageRetriever interface
	 */
	public void printSequence(String s) {
		if (textAreaRowCount==0)
			initS=1+"     "+ s;//save
		if (textAreaRowCount==1)
			initS=initS+"\n"+2+"     "+s;//save
		if (textAreaRowCount>maxLines){
			textArea.selectAll();
			textArea.replaceSelection("");
			textArea.append(initS+" \n");
			textArea.append(3+" ... "+"a maximum of about "+(maxLines+1)+" lines are displayed!"+" \n\n");
			textAreaRowCount=4;
		}
		textArea.append((textAreaRowCount+1)+"     "+s + " \n");
		textAreaRowCount++;
	}
	
	/**
	 * Changes the color based on pre-defined color selection!
	 */
	private void changeColor() {
		int ir = 0;
		int ig = 0;
		int ib = 0;
		int ia = 0;
		String s = (String) selColorCb.getSelectedItem();
		
		if (s.compareTo(colors[0]) == 0)// AliceBlue
		{
			c = new Color(240, 248, 255);
			ir = c.getRed();
			ig = c.getGreen();
			ib = c.getBlue();
			ia = c.getAlpha();
		} else if (s.compareTo(colors[1]) == 0)// Antique White
		{
			c = new Color(250, 235, 215);
			ir = c.getRed();
			ig = c.getGreen();
			ib = c.getBlue();
			ia = c.getAlpha();
		} else if (s.compareTo(colors[2]) == 0)// AquaMarine
		{
			c = new Color(127, 255, 212);
			ir = c.getRed();
			ig = c.getGreen();
			ib = c.getBlue();
			ia = c.getAlpha();
		} else if (s.compareTo(colors[3]) == 0)// Azure
		{
			c = new Color(240, 255, 255);
			ir = c.getRed();
			ig = c.getGreen();
			ib = c.getBlue();
			ia = c.getAlpha();
		} else if (s.compareTo(colors[4]) == 0)// Beige
		{
			c = new Color(245, 245, 220);
			ir = c.getRed();
			ig = c.getGreen();
			ib = c.getBlue();
			ia = c.getAlpha();
		} else if (s.compareTo(colors[5]) == 0)// Bisque
		{
			c = new Color(255, 228, 196);
			ir = c.getRed();
			ig = c.getGreen();
			ib = c.getBlue();
			ia = c.getAlpha();
		} else if (s.compareTo(colors[6]) == 0)// Black
		{
			c = Color.BLACK;
			ir = c.getRed();
			ig = c.getGreen();
			ib = c.getBlue();
			ia = c.getAlpha();
		} else if (s.compareTo(colors[7]) == 0)// BlanchedAlmond
		{
			c = new Color(255, 235, 205);
			ir = c.getRed();
			ig = c.getGreen();
			ib = c.getBlue();
			ia = c.getAlpha();
		} else if (s.compareTo(colors[8]) == 0)// Blue
		{
			c = Color.BLUE;
			ir = c.getRed();
			ig = c.getGreen();
			ib = c.getBlue();
			ia = c.getAlpha();
		} else if (s.compareTo(colors[9]) == 0)// BlueViolet
		{
			c = new Color(138, 43, 226);
			ir = c.getRed();
			ig = c.getGreen();
			ib = c.getBlue();
			ia = c.getAlpha();
		} else if (s.compareTo(colors[10]) == 0)// BrightGold
		{
			c = new Color(iBrightGold);
			ir = c.getRed();
			ig = c.getGreen();
			ib = c.getBlue();
			ia = c.getAlpha();
		} else if (s.compareTo(colors[11]) == 0)// Brown
		{
			c = new Color(165, 42, 42);
			ir = c.getRed();
			ig = c.getGreen();
			ib = c.getBlue();
			ia = c.getAlpha();
		} else if (s.compareTo(colors[12]) == 0)// BurlyWood
		{
			c = new Color(222, 184, 135);
			ir = c.getRed();
			ig = c.getGreen();
			ib = c.getBlue();
			ia = c.getAlpha();
		} else if (s.compareTo(colors[13]) == 0)// CadetBlue
		{
			c = new Color(95, 158, 160);
			ir = c.getRed();
			ig = c.getGreen();
			ib = c.getBlue();
			ia = c.getAlpha();
		} else if (s.compareTo(colors[14]) == 0)// Chartereuse
		{
			c = new Color(127, 255, 0);
			ir = c.getRed();
			ig = c.getGreen();
			ib = c.getBlue();
			ia = c.getAlpha();
		} else if (s.compareTo(colors[15]) == 0)// Chocolate
		{
			c = new Color(210, 105, 30);
			ir = c.getRed();
			ig = c.getGreen();
			ib = c.getBlue();
			ia = c.getAlpha();
		} else if (s.compareTo(colors[16]) == 0)// Copper
		{
			c = new Color(iCopper);
			ir = c.getRed();
			ig = c.getGreen();
			ib = c.getBlue();
			ia = c.getAlpha();
		} else if (s.compareTo(colors[17]) == 0)// Coral
		{
			c = new Color(iCoral);
			ir = c.getRed();
			ig = c.getGreen();
			ib = c.getBlue();
			ia = c.getAlpha();
		} else if (s.compareTo(colors[18]) == 0)// CornflowerBlue
		{
			c = new Color(100, 149, 237);
			ir = c.getRed();
			ig = c.getGreen();
			ib = c.getBlue();
			ia = c.getAlpha();
		} else if (s.compareTo(colors[19]) == 0)// Cornsilk
		{
			c = new Color(255, 248, 220);
			ir = c.getRed();
			ig = c.getGreen();
			ib = c.getBlue();
			ia = c.getAlpha();
		} else if (s.compareTo(colors[20]) == 0)// Cyan
		{
			c = Color.CYAN;
			ir = c.getRed();
			ig = c.getGreen();
			ib = c.getBlue();
			ia = c.getAlpha();
		} else if (s.compareTo(colors[21]) == 0)// DarkGoldenrod
		{
			c = new Color(184, 134, 11);
			ir = c.getRed();
			ig = c.getGreen();
			ib = c.getBlue();
			ia = c.getAlpha();
		} else if (s.compareTo(colors[22]) == 0)// DarkGreen
		{
			c = new Color(0, 100, 0);
			ir = c.getRed();
			ig = c.getGreen();
			ib = c.getBlue();
			ia = c.getAlpha();
		} else if (s.compareTo(colors[23]) == 0)// DarkGrey
		{
			c = Color.darkGray;
			ir = c.getRed();
			ig = c.getGreen();
			ib = c.getBlue();
			ia = c.getAlpha();
		} else if (s.compareTo(colors[24]) == 0)// DarkKhaki
		{
			c = new Color(189, 183, 107);
			ir = c.getRed();
			ig = c.getGreen();
			ib = c.getBlue();
			ia = c.getAlpha();
		} else if (s.compareTo(colors[25]) == 0)// DarkOliveGreen
		{
			c = new Color(85, 107, 47);
			ir = c.getRed();
			ig = c.getGreen();
			ib = c.getBlue();
			ia = c.getAlpha();
		} else if (s.compareTo(colors[26]) == 0)// DarkOrange
		{
			c = new Color(255, 140, 0);
			ir = c.getRed();
			ig = c.getGreen();
			ib = c.getBlue();
			ia = c.getAlpha();
		} else if (s.compareTo(colors[27]) == 0)// DarkOrchid
		{
			c = new Color(153, 50, 204);
			ir = c.getRed();
			ig = c.getGreen();
			ib = c.getBlue();
			ia = c.getAlpha();
		} else if (s.compareTo(colors[28]) == 0)// DarkSalmon
		{
			c = new Color(233, 150, 122);
			ir = c.getRed();
			ig = c.getGreen();
			ib = c.getBlue();
			ia = c.getAlpha();
		} else if (s.compareTo(colors[29]) == 0)// DarkSeaGreen
		{
			c = new Color(143, 188, 143);
			ir = c.getRed();
			ig = c.getGreen();
			ib = c.getBlue();
			ia = c.getAlpha();
		} else if (s.compareTo(colors[30]) == 0)// DarkSlateBlue
		{
			c = new Color(72, 61, 139);
			ir = c.getRed();
			ig = c.getGreen();
			ib = c.getBlue();
			ia = c.getAlpha();
		} else if (s.compareTo(colors[31]) == 0)// DarkSlateGray
		{
			c = new Color(47, 79, 79);
			ir = c.getRed();
			ig = c.getGreen();
			ib = c.getBlue();
			ia = c.getAlpha();
		} else if (s.compareTo(colors[32]) == 0)// DarkTurquoise
		{
			c = new Color(0, 206, 209);
			ir = c.getRed();
			ig = c.getGreen();
			ib = c.getBlue();
			ia = c.getAlpha();
		} else if (s.compareTo(colors[33]) == 0)// DarkViolet
		{
			c = new Color(148, 0, 211);
			ir = c.getRed();
			ig = c.getGreen();
			ib = c.getBlue();
			ia = c.getAlpha();
		} else if (s.compareTo(colors[34]) == 0)// DeepPink
		{
			c = new Color(255, 20, 147);
			ir = c.getRed();
			ig = c.getGreen();
			ib = c.getBlue();
			ia = c.getAlpha();
		} else if (s.compareTo(colors[35]) == 0)// DeepSkyBlue
		{
			c = new Color(0, 191, 255);
			ir = c.getRed();
			ig = c.getGreen();
			ib = c.getBlue();
			ia = c.getAlpha();
		} else if (s.compareTo(colors[36]) == 0)// DimGray
		{
			c = new Color(105, 105, 105);
			ir = c.getRed();
			ig = c.getGreen();
			ib = c.getBlue();
			ia = c.getAlpha();
		} else if (s.compareTo(colors[37]) == 0)// DodgerBlue
		{
			c = new Color(30, 144, 255);
			ir = c.getRed();
			ig = c.getGreen();
			ib = c.getBlue();
			ia = c.getAlpha();
		} else if (s.compareTo(colors[38]) == 0)// DustyRose
		{
			c = new Color(iDustyRose);
			ir = c.getRed();
			ig = c.getGreen();
			ib = c.getBlue();
			ia = c.getAlpha();
		} else if (s.compareTo(colors[39]) == 0)// firebrick
		{
			c = new Color(178, 34, 34);
			ir = c.getRed();
			ig = c.getGreen();
			ib = c.getBlue();
			ia = c.getAlpha();
		} else if (s.compareTo(colors[40]) == 0)// floralWhite
		{
			c = new Color(255, 250, 240);
			ir = c.getRed();
			ig = c.getGreen();
			ib = c.getBlue();
			ia = c.getAlpha();
		} else if (s.compareTo(colors[41]) == 0)// ForestGreen
		{
			c = new Color(iForestGreen);
			ir = c.getRed();
			ig = c.getGreen();
			ib = c.getBlue();
			ia = c.getAlpha();
		} else if (s.compareTo(colors[42]) == 0)// gainsboro
		{
			c = new Color(220, 220, 220);
			ir = c.getRed();
			ig = c.getGreen();
			ib = c.getBlue();
			ia = c.getAlpha();
		} else if (s.compareTo(colors[43]) == 0)// ghostwhite
		{
			c = new Color(248, 248, 255);
			ir = c.getRed();
			ig = c.getGreen();
			ib = c.getBlue();
			ia = c.getAlpha();
		} else if (s.compareTo(colors[44]) == 0)// gold
		{
			c = new Color(255, 215, 0);
			ir = c.getRed();
			ig = c.getGreen();
			ib = c.getBlue();
			ia = c.getAlpha();
		} else if (s.compareTo(colors[45]) == 0)// goldenrod
		{
			c = new Color(218, 165, 32);
			ir = c.getRed();
			ig = c.getGreen();
			ib = c.getBlue();
			ia = c.getAlpha();
		} else if (s.compareTo(colors[46]) == 0)// Gray
		{
			c = Color.GRAY;
			ir = c.getRed();
			ig = c.getGreen();
			ib = c.getBlue();
			ia = c.getAlpha();
		} else if (s.compareTo(colors[47]) == 0)// Green
		{
			c = Color.GREEN;
			ir = c.getRed();
			ig = c.getGreen();
			ib = c.getBlue();
			ia = c.getAlpha();
		} else if (s.compareTo(colors[48]) == 0)// GreenYellow
		{
			c = new Color(173, 255, 47);
			ir = c.getRed();
			ig = c.getGreen();
			ib = c.getBlue();
			ia = c.getAlpha();
		} else if (s.compareTo(colors[49]) == 0)// Honeydew
		{
			c = new Color(240, 255, 240);
			ir = c.getRed();
			ig = c.getGreen();
			ib = c.getBlue();
			ia = c.getAlpha();
		} else if (s.compareTo(colors[50]) == 0)// HotPink
		{
			c = new Color(255, 105, 180);
			ir = c.getRed();
			ig = c.getGreen();
			ib = c.getBlue();
			ia = c.getAlpha();
		} else if (s.compareTo(colors[51]) == 0)// IndianRed
		{
			c = new Color(205, 92, 92);
			ir = c.getRed();
			ig = c.getGreen();
			ib = c.getBlue();
			ia = c.getAlpha();
		} else if (s.compareTo(colors[52]) == 0)// Ivory
		{
			c = new Color(255, 255, 240);
			ir = c.getRed();
			ig = c.getGreen();
			ib = c.getBlue();
			ia = c.getAlpha();
		} else if (s.compareTo(colors[53]) == 0)// Khaki
		{
			c = new Color(iKhaki);
			ir = c.getRed();
			ig = c.getGreen();
			ib = c.getBlue();
			ia = c.getAlpha();
		} else if (s.compareTo(colors[54]) == 0)// Lavender
		{
			c = new Color(230, 230, 250);
			ir = c.getRed();
			ig = c.getGreen();
			ib = c.getBlue();
			ia = c.getAlpha();
		} else if (s.compareTo(colors[55]) == 0)// LavenderBlush
		{
			c = new Color(255, 240, 245);
			ir = c.getRed();
			ig = c.getGreen();
			ib = c.getBlue();
			ia = c.getAlpha();
		} else if (s.compareTo(colors[56]) == 0)// LawnGreen
		{
			c = new Color(124, 252, 0);
			ir = c.getRed();
			ig = c.getGreen();
			ib = c.getBlue();
			ia = c.getAlpha();
		} else if (s.compareTo(colors[57]) == 0)// LemonChiffon
		{
			c = new Color(255, 250, 205);
			ir = c.getRed();
			ig = c.getGreen();
			ib = c.getBlue();
			ia = c.getAlpha();
		} else if (s.compareTo(colors[58]) == 0)// LightBlue
		{
			c = new Color(173, 216, 230);
			ir = c.getRed();
			ig = c.getGreen();
			ib = c.getBlue();
			ia = c.getAlpha();
		} else if (s.compareTo(colors[59]) == 0)// LightCoral
		{
			c = new Color(240, 128, 128);
			ir = c.getRed();
			ig = c.getGreen();
			ib = c.getBlue();
			ia = c.getAlpha();
		} else if (s.compareTo(colors[60]) == 0)// LightCyan
		{
			c = new Color(224, 255, 255);
			ir = c.getRed();
			ig = c.getGreen();
			ib = c.getBlue();
			ia = c.getAlpha();
		} else if (s.compareTo(colors[61]) == 0)// LightGoldenrod
		{
			c = new Color(238, 221, 130);
			ir = c.getRed();
			ig = c.getGreen();
			ib = c.getBlue();
			ia = c.getAlpha();
		} else if (s.compareTo(colors[62]) == 0)// LightGoldenrodYellow
		{
			c = new Color(250, 250, 210);
			ir = c.getRed();
			ig = c.getGreen();
			ib = c.getBlue();
			ia = c.getAlpha();
		} else if (s.compareTo(colors[63]) == 0)// LightGray
		{
			c = Color.lightGray;
			ir = c.getRed();
			ig = c.getGreen();
			ib = c.getBlue();
			ia = c.getAlpha();
		} else if (s.compareTo(colors[64]) == 0)// LightPink
		{
			c = new Color(255, 182, 193);
			ir = c.getRed();
			ig = c.getGreen();
			ib = c.getBlue();
			ia = c.getAlpha();
		} else if (s.compareTo(colors[65]) == 0)// LightSalmon
		{
			c = new Color(255, 160, 122);
			ir = c.getRed();
			ig = c.getGreen();
			ib = c.getBlue();
			ia = c.getAlpha();
		} else if (s.compareTo(colors[66]) == 0)// LightSeaGreen
		{
			c = new Color(32, 178, 170);
			ir = c.getRed();
			ig = c.getGreen();
			ib = c.getBlue();
			ia = c.getAlpha();
		} else if (s.compareTo(colors[67]) == 0)// LightSkyBlue
		{
			c = new Color(135, 206, 250);
			ir = c.getRed();
			ig = c.getGreen();
			ib = c.getBlue();
			ia = c.getAlpha();
		} else if (s.compareTo(colors[68]) == 0)// LightSlateBlue
		{
			c = new Color(132, 112, 255);
			ir = c.getRed();
			ig = c.getGreen();
			ib = c.getBlue();
			ia = c.getAlpha();
		} else if (s.compareTo(colors[69]) == 0)// LightSlateGray
		{
			c = new Color(119, 136, 153);
			ir = c.getRed();
			ig = c.getGreen();
			ib = c.getBlue();
			ia = c.getAlpha();
		} else if (s.compareTo(colors[70]) == 0)// LightSteelBlue
		{
			c = new Color(176, 196, 222);
			ir = c.getRed();
			ig = c.getGreen();
			ib = c.getBlue();
			ia = c.getAlpha();
		} else if (s.compareTo(colors[71]) == 0)// LightYellow
		{
			c = new Color(255, 255, 224);
			ir = c.getRed();
			ig = c.getGreen();
			ib = c.getBlue();
			ia = c.getAlpha();
		} else if (s.compareTo(colors[72]) == 0)// LimeGreen
		{
			c = new Color(50, 205, 50);
			ir = c.getRed();
			ig = c.getGreen();
			ib = c.getBlue();
			ia = c.getAlpha();
		} else if (s.compareTo(colors[73]) == 0)// Linen
		{
			c = new Color(250, 240, 230);
			ir = c.getRed();
			ig = c.getGreen();
			ib = c.getBlue();
			ia = c.getAlpha();
		} else if (s.compareTo(colors[74]) == 0)// Magenta
		{
			c = Color.MAGENTA;
			ir = c.getRed();
			ig = c.getGreen();
			ib = c.getBlue();
			ia = c.getAlpha();
		} else if (s.compareTo(colors[75]) == 0)// Maroon
		{
			c = new Color(176, 48, 96);
			ir = c.getRed();
			ig = c.getGreen();
			ib = c.getBlue();
			ia = c.getAlpha();
		} else if (s.compareTo(colors[76]) == 0)// MediumAquamarine
		{
			c = new Color(102, 205, 170);
			ir = c.getRed();
			ig = c.getGreen();
			ib = c.getBlue();
			ia = c.getAlpha();
		} else if (s.compareTo(colors[77]) == 0)// MediumBlue
		{
			c = new Color(0, 0, 205);
			ir = c.getRed();
			ig = c.getGreen();
			ib = c.getBlue();
			ia = c.getAlpha();
		} else if (s.compareTo(colors[78]) == 0)// MediumOrchid
		{
			c = new Color(186, 85, 211);
			ir = c.getRed();
			ig = c.getGreen();
			ib = c.getBlue();
			ia = c.getAlpha();
		} else if (s.compareTo(colors[79]) == 0)// MediumPurple
		{
			c = new Color(147, 112, 219);
			ir = c.getRed();
			ig = c.getGreen();
			ib = c.getBlue();
			ia = c.getAlpha();
		} else if (s.compareTo(colors[80]) == 0)// MediumSeaGreen
		{
			c = new Color(60, 179, 113);
			ir = c.getRed();
			ig = c.getGreen();
			ib = c.getBlue();
			ia = c.getAlpha();
		} else if (s.compareTo(colors[81]) == 0)// MediumSlateBlue
		{
			c = new Color(123, 104, 238);
			ir = c.getRed();
			ig = c.getGreen();
			ib = c.getBlue();
			ia = c.getAlpha();
		} else if (s.compareTo(colors[82]) == 0)// MediumSpringGreen
		{
			c = new Color(0, 250, 154);
			ir = c.getRed();
			ig = c.getGreen();
			ib = c.getBlue();
			ia = c.getAlpha();
		} else if (s.compareTo(colors[83]) == 0)// MediumTurquoise
		{
			c = new Color(72, 209, 204);
			ir = c.getRed();
			ig = c.getGreen();
			ib = c.getBlue();
			ia = c.getAlpha();
		} else if (s.compareTo(colors[84]) == 0)// MediumVioletRed
		{
			c = new Color(199, 21, 133);
			ir = c.getRed();
			ig = c.getGreen();
			ib = c.getBlue();
			ia = c.getAlpha();
		} else if (s.compareTo(colors[85]) == 0)// MidnightBlue
		{
			c = new Color(iMidnightBlue);
			ir = c.getRed();
			ig = c.getGreen();
			ib = c.getBlue();
			ia = c.getAlpha();
		} else if (s.compareTo(colors[86]) == 0)// MintCream
		{
			c = new Color(245, 255, 250);
			ir = c.getRed();
			ig = c.getGreen();
			ib = c.getBlue();
			ia = c.getAlpha();
		} else if (s.compareTo(colors[87]) == 0)// MistyRose
		{
			c = new Color(255, 228, 225);
			ir = c.getRed();
			ig = c.getGreen();
			ib = c.getBlue();
			ia = c.getAlpha();
		} else if (s.compareTo(colors[88]) == 0)// Moccasin
		{
			c = new Color(255, 228, 181);
			ir = c.getRed();
			ig = c.getGreen();
			ib = c.getBlue();
			ia = c.getAlpha();
		} else if (s.compareTo(colors[89]) == 0)// NavajoWhite
		{
			c = new Color(255, 222, 173);
			ir = c.getRed();
			ig = c.getGreen();
			ib = c.getBlue();
			ia = c.getAlpha();
		} else if (s.compareTo(colors[90]) == 0)// Navy
		{
			c = new Color(0, 0, 128);
			ir = c.getRed();
			ig = c.getGreen();
			ib = c.getBlue();
			ia = c.getAlpha();
		} else if (s.compareTo(colors[91]) == 0)// NeonPink
		{
			c = new Color(iNeonPink);
			ir = c.getRed();
			ig = c.getGreen();
			ib = c.getBlue();
			ia = c.getAlpha();
		} else if (s.compareTo(colors[92]) == 0)// OldLace
		{
			c = new Color(253, 245, 230);
			ir = c.getRed();
			ig = c.getGreen();
			ib = c.getBlue();
			ia = c.getAlpha();
		} else if (s.compareTo(colors[93]) == 0)// OliveDrab
		{
			c = new Color(107, 142, 35);
			ir = c.getRed();
			ig = c.getGreen();
			ib = c.getBlue();
			ia = c.getAlpha();
		} else if (s.compareTo(colors[94]) == 0)// orange
		{
			c = Color.ORANGE;
			ir = c.getRed();
			ig = c.getGreen();
			ib = c.getBlue();
			ia = c.getAlpha();
		} else if (s.compareTo(colors[95]) == 0)// OrangeRed
		{
			c = new Color(255, 69, 0);
			ir = c.getRed();
			ig = c.getGreen();
			ib = c.getBlue();
			ia = c.getAlpha();
		} else if (s.compareTo(colors[96]) == 0)// Orchid
		{
			c = new Color(218, 112, 214);
			ir = c.getRed();
			ig = c.getGreen();
			ib = c.getBlue();
			ia = c.getAlpha();
		} else if (s.compareTo(colors[97]) == 0)// PaleGoldenrod
		{
			c = new Color(238, 232, 170);
			ir = c.getRed();
			ig = c.getGreen();
			ib = c.getBlue();
			ia = c.getAlpha();
		} else if (s.compareTo(colors[98]) == 0)// PaleGreen
		{
			c = new Color(152, 251, 152);
			ir = c.getRed();
			ig = c.getGreen();
			ib = c.getBlue();
			ia = c.getAlpha();
		} else if (s.compareTo(colors[99]) == 0)// PaleTurquoise
		{
			c = new Color(175, 238, 238);
			ir = c.getRed();
			ig = c.getGreen();
			ib = c.getBlue();
			ia = c.getAlpha();
		} else if (s.compareTo(colors[100]) == 0)// PaleVioletRed
		{
			c = new Color(219, 112, 147);
			ir = c.getRed();
			ig = c.getGreen();
			ib = c.getBlue();
			ia = c.getAlpha();
		} else if (s.compareTo(colors[101]) == 0)// PapayaWhip
		{
			c = new Color(255, 239, 213);
			ir = c.getRed();
			ig = c.getGreen();
			ib = c.getBlue();
			ia = c.getAlpha();
		} else if (s.compareTo(colors[102]) == 0)// PeachPuff
		{
			c = new Color(255, 218, 185);
			ir = c.getRed();
			ig = c.getGreen();
			ib = c.getBlue();
			ia = c.getAlpha();
		} else if (s.compareTo(colors[103]) == 0)// Peru
		{
			c = new Color(205, 133, 63);
			ir = c.getRed();
			ig = c.getGreen();
			ib = c.getBlue();
			ia = c.getAlpha();
		} else if (s.compareTo(colors[104]) == 0)// Pink
		{
			c = Color.PINK;
			ir = c.getRed();
			ig = c.getGreen();
			ib = c.getBlue();
			ia = c.getAlpha();
		}

		else if (s.compareTo(colors[105]) == 0)// Plum
		{
			c = new Color(0.6f, 0.3f, 0.3f);
			ir = c.getRed();
			ig = c.getGreen();
			ib = c.getBlue();
			ia = c.getAlpha();
		} else if (s.compareTo(colors[106]) == 0)// PowderBlue
		{
			c = new Color(176, 224, 230);
			ir = c.getRed();
			ig = c.getGreen();
			ib = c.getBlue();
			ia = c.getAlpha();
		} else if (s.compareTo(colors[107]) == 0)// PreDawn
		{
			c = new Color(0.4f, 0.5f, 0.5f);
			ir = c.getRed();
			ig = c.getGreen();
			ib = c.getBlue();
			ia = c.getAlpha();
		} else if (s.compareTo(colors[108]) == 0)// Purple
		{
			c = new Color(160, 32, 240);
			ir = c.getRed();
			ig = c.getGreen();
			ib = c.getBlue();
			ia = c.getAlpha();
		} else if (s.compareTo(colors[109]) == 0)// Red
		{
			c = Color.RED;
			ir = c.getRed();
			ig = c.getGreen();
			ib = c.getBlue();
			ia = c.getAlpha();
		} else if (s.compareTo(colors[110]) == 0)// RosyBrown
		{
			c = new Color(188, 143, 143);
			ir = c.getRed();
			ig = c.getGreen();
			ib = c.getBlue();
			ia = c.getAlpha();
		} else if (s.compareTo(colors[111]) == 0)// RoyalBlue
		{
			c = new Color(65, 105, 225);
			ir = c.getRed();
			ig = c.getGreen();
			ib = c.getBlue();
			ia = c.getAlpha();
		} else if (s.compareTo(colors[112]) == 0)// SaddleBrown
		{
			c = new Color(139, 69, 19);
			ir = c.getRed();
			ig = c.getGreen();
			ib = c.getBlue();
			ia = c.getAlpha();
		} else if (s.compareTo(colors[113]) == 0)// Salmon
		{
			c = new Color(iSalmon);
			ir = c.getRed();
			ig = c.getGreen();
			ib = c.getBlue();
			ia = c.getAlpha();
		} else if (s.compareTo(colors[114]) == 0)// SandyBrown
		{
			c = new Color(244, 164, 96);
			ir = c.getRed();
			ig = c.getGreen();
			ib = c.getBlue();
			ia = c.getAlpha();
		} else if (s.compareTo(colors[115]) == 0)// SeaGreen
		{
			c = new Color(46, 139, 87);
			ir = c.getRed();
			ig = c.getGreen();
			ib = c.getBlue();
			ia = c.getAlpha();
		} else if (s.compareTo(colors[116]) == 0)// Seashell
		{
			c = new Color(255, 245, 238);
			ir = c.getRed();
			ig = c.getGreen();
			ib = c.getBlue();
			ia = c.getAlpha();
		} else if (s.compareTo(colors[117]) == 0)// Sienna
		{
			c = new Color(160, 82, 45);
			ir = c.getRed();
			ig = c.getGreen();
			ib = c.getBlue();
			ia = c.getAlpha();
		} else if (s.compareTo(colors[118]) == 0)// SkyBlue
		{
			c = new Color(135, 206, 235);
			ir = c.getRed();
			ig = c.getGreen();
			ib = c.getBlue();
			ia = c.getAlpha();
		} else if (s.compareTo(colors[119]) == 0)// SlateBlue
		{
			c = new Color(106, 90, 205);
			ir = c.getRed();
			ig = c.getGreen();
			ib = c.getBlue();
			ia = c.getAlpha();
		} else if (s.compareTo(colors[120]) == 0)// SlateGray
		{
			c = new Color(112, 128, 144);
			ir = c.getRed();
			ig = c.getGreen();
			ib = c.getBlue();
			ia = c.getAlpha();
		} else if (s.compareTo(colors[121]) == 0)// Snow
		{
			c = new Color(255, 250, 250);
			ir = c.getRed();
			ig = c.getGreen();
			ib = c.getBlue();
			ia = c.getAlpha();
		} else if (s.compareTo(colors[122]) == 0)// SpringGreen
		{
			c = new Color(0, 255, 127);
			ir = c.getRed();
			ig = c.getGreen();
			ib = c.getBlue();
			ia = c.getAlpha();
		} else if (s.compareTo(colors[123]) == 0)// SteelBlue
		{
			c = new Color(70, 130, 180);
			ir = c.getRed();
			ig = c.getGreen();
			ib = c.getBlue();
			ia = c.getAlpha();
		} else if (s.compareTo(colors[124]) == 0)// Tan
		{
			c = new Color(iTan);
			ir = c.getRed();
			ig = c.getGreen();
			ib = c.getBlue();
			ia = c.getAlpha();
		} else if (s.compareTo(colors[125]) == 0)// Thistle
		{
			c = new Color(216, 191, 216);
			ir = c.getRed();
			ig = c.getGreen();
			ib = c.getBlue();
			ia = c.getAlpha();
		} else if (s.compareTo(colors[126]) == 0)// Tomato
		{
			c = new Color(255, 99, 71);
			ir = c.getRed();
			ig = c.getGreen();
			ib = c.getBlue();
			ia = c.getAlpha();
		} else if (s.compareTo(colors[127]) == 0)// Turquoise
		{
			c = new Color(64, 224, 208);
			ir = c.getRed();
			ig = c.getGreen();
			ib = c.getBlue();
			ia = c.getAlpha();
		} else if (s.compareTo(colors[128]) == 0)// Violet
		{
			c = new Color(238, 130, 238);
			ir = c.getRed();
			ig = c.getGreen();
			ib = c.getBlue();
			ia = c.getAlpha();
		} else if (s.compareTo(colors[129]) == 0)// Wheat
		{
			c = new Color(245, 222, 179);
			ir = c.getRed();
			ig = c.getGreen();
			ib = c.getBlue();
			ia = c.getAlpha();
		} else if (s.compareTo(colors[130]) == 0)// White
		{
			c = Color.WHITE;
			ir = c.getRed();
			ig = c.getGreen();
			ib = c.getBlue();
			ia = c.getAlpha();
		} else if (s.compareTo(colors[131]) == 0)// WhiteSmoke
		{
			c = new Color(245, 245, 245);
			ir = c.getRed();
			ig = c.getGreen();
			ib = c.getBlue();
			ia = c.getAlpha();
		} else if (s.compareTo(colors[132]) == 0)// Yellow
		{
			c = Color.YELLOW;
			ir = c.getRed();
			ig = c.getGreen();
			ib = c.getBlue();
			ia = c.getAlpha();
		} else if (s.compareTo(colors[133]) == 0)// YellowGreen
		{
			c = new Color(154, 205, 50);
			ir = c.getRed();
			ig = c.getGreen();
			ib = c.getBlue();
			ia = c.getAlpha();
		}

		sR.setValue(ir);
		sG.setValue(ig);
		sB.setValue(ib);
		sA.setValue(ia);// here a ChangeListener is auto-fire-ed!!!
		updateText(c);
	}

	
}
