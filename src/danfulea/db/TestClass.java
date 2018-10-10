package danfulea.db;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.FlowLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.math.BigDecimal;
import java.sql.Connection;
import java.sql.Date;
import java.sql.Timestamp;
import java.text.ParseException;
import java.util.ResourceBundle;

import javax.swing.BoxLayout;
import javax.swing.ImageIcon;
import javax.swing.JButton;
import javax.swing.JComboBox;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTable;
import javax.swing.JTextField;
import javax.swing.ListSelectionModel;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;

import danfulea.math.Convertor;
import danfulea.utils.AboutFrame;
import danfulea.utils.FileOperation;
import danfulea.utils.FrameUtilities;
import danfulea.utils.LookAndFeel;
import danfulea.utils.ScanDiskLFGui;

/**
 * Test Class for database operation. It also uses the generic About frame and
 * look and feel.
 * 
 * @author Dan Fulea, 04 AUG. 2016
 *
 */

@SuppressWarnings("unused")
public class TestClass extends JFrame implements ActionListener, ItemListener, ListSelectionListener {

	// general
	private static final long serialVersionUID = 1L;

	/**
	 * Frame preferred dimension
	 */
	private final Dimension PREFERRED_SIZE = new Dimension(1000, 500);

	/**
	 * Table preferred dimension
	 */
	private final Dimension tableDimension = new Dimension(400, 200);

	/**
	 * Combobox preferred dimension
	 */
	private final Dimension sizeCb = new Dimension(150, 27);

	/**
	 * Background color
	 */
	public static Color bkgColor = new Color(230, 255, 210, 255);// Linux mint
																	// green
																	// alike
	/**
	 * Foreground color
	 */
	public static Color foreColor = Color.black;// Color.white;

	// the resources
	/**
	 * The resource location
	 */
	private static final String BASE_RESOURCE_CLASS = "danfulea.resources.DanfuleaResources";

	/**
	 * the resource object
	 */
	private ResourceBundle resources;

	// Commands==========================================================
	/**
	 * the command string
	 */
	private String command = null;

	private static final String EXIT_COMMAND = "EXIT";
	private static final String ABOUT_COMMAND = "ABOUT";
	private static final String LOOKANDFEEL_COMMAND = "LOOKANDFEEL";
	private static final String ADD_COMMAND = "ADD";
	private static final String DELETE_COMMAND = "DELETE";
	private static final String DELETEALL_COMMAND = "DELETEALL";
	private static final String UPDATE_COMMAND = "UPDATE";
	private static final String CLEAR_COMMAND = "CLEAR";

	private static final String NESTED_ADD_COMMAND = "NESTED_ADD";
	private static final String NESTED_DELETE_COMMAND = "NESTED_DELETE";
	private static final String NESTED_DELETEALL_COMMAND = "NESTED_DELETEALL";
	private static final String NESTED_UPDATE_COMMAND = "NESTED_UPDATE";
	private static final String NESTED_CLEAR_COMMAND = "NESTED_CLEAR";
	// =====================================================================================

	// GUI elements
	private JTextField userNameTf = new JTextField(20);
	private JTextField passwordTf = new JTextField(20);
	private JComboBox<String> orderbyCb;
	private JTextField emailTf = new JTextField(20);
	private JTextField addressTf = new JTextField(20);
	private JComboBox<String> nestedorderbyCb;

	// DB data
	/**
	 * The connection
	 */
	private Connection dbcon = null;
	
	/**
	 * The database agent associated to main table
	 */
	private DatabaseAgentSupport dbagent;
	
	/**
	 * The database agent associated to nested table
	 */
	private DatabaseAgentSupport nesteddbagent;
	
	/**
	 * Database name
	 */
	private String maindbname;
	
	/**
	 * Database main table name
	 */
	private String maindbTable;
	
	/**
	 * Database nested table name
	 */
	private String nesteddbTable;
	
	///**
	// * Index of main table primary key column (first column having index of 0)
	// */
	//private int mainPrimaryKeyColumnIndex = 0;
	
	//**
	// * Index of nested table primary key column (first column having index of 0)
	// */
	//private int nestedPrimaryKeyColumnIndex = 0;

	/**
	 * Main table primary key column name
	 */
	private String mainTablePrimaryKey = "ID";
	
	/**
	 * Nested table primary key column name
	 */
	private String nestedTablePrimaryKey = "NRCRT";
	
	/**
	 * Shared column name for main table and nested table
	 */
	private String IDlink = "ID";
	
	/**
	 * The JTable component associated to main table
	 */
	private JTable mainTable;
	
	/**
	 * The column used for sorting data in main table (ORDER BY SQL syntax)
	 */
	private String orderbyS = "ID";
	
	/**
	 * The JTable component associated to nested table
	 */
	private JTable nestedTable;
	
	/**
	 * The column used for sorting data in nested table (ORDER BY SQL syntax)
	 */
	private String nestedorderbyS = "NRCRT";
	

	/**
	 * Constructor. Some initializations are set here.
	 */
	public TestClass() {
		resources = ResourceBundle.getBundle(BASE_RESOURCE_CLASS);

		maindbname = "testDb";
		
		maindbTable = "testTable3";//"testTable";//"testTable3";//"testTable";
		mainTablePrimaryKey = "ID";// "NRCRT";//"ID";//
		
		nesteddbTable = "testTable2";
		nestedTablePrimaryKey = "NRCRT";
		
		IDlink = "ID";

		// DERBY
		DatabaseAgent.ID_CONNECTION = DatabaseAgent.DERBY_CONNECTION;
		startDerbyConnection();// with mainDB

		// MySQL
		// DatabaseAgent.ID_CONNECTION = DatabaseAgent.MYSQL_CONNECTION;
		// dbcon = DatabaseAgent.getConnection(maindbname, "root", "");
		// dbagent = new DatabaseAgentSupport(dbcon, mainTablePrimaryKey,maindbTable);
		// nesteddbagent = new DatabaseAgentSupport(dbcon, nestedTablePrimaryKey,nesteddbTable);

		// PostgreSQL
		 //DatabaseAgent.ID_CONNECTION = DatabaseAgent.POSTGRESQL_CONNECTION;
		// dbcon = DatabaseAgent.getConnection(maindbname, "postgres", "admin");
		 //dbagent = new DatabaseAgentSupport(dbcon, mainTablePrimaryKey,maindbTable);//maindbname, "postgres", "admin");
		 //nesteddbagent = new DatabaseAgentSupport(dbcon, nestedTablePrimaryKey,nesteddbTable);
		
		dbagent.setHasValidAIColumn(false);//(true);//(false);//(true);
		nesteddbagent.setHasValidAIColumn(true);

		this.setTitle(resources.getString("Application.NAME"));// this.setTitle("Fast
																// DB");
		setDefaultCloseOperation(JFrame.DO_NOTHING_ON_CLOSE);
		addWindowListener(new WindowAdapter() {
			public void windowClosing(WindowEvent e) {
				attemptExit();
			}
		});

		performQueryDb();
		createGUI();
		setDefaultLookAndFeelDecorated(true);
		FrameUtilities.createImageIcon// ("/danfulea/resources/duke.png", this);
		(this.resources.getString("form.icon.url"), this);

		FrameUtilities.centerFrameOnScreen(this);

		setVisible(true);
	}

	// @SuppressWarnings("unused")
	/**
	 * Starts derby connection and initializes the agents.
	 */
	private void startDerbyConnection() {

		String datas = this.resources.getString("data.load");// "Data";
		String currentDir = System.getProperty("user.dir");
		String file_sep = System.getProperty("file.separator");
		String opens = currentDir + file_sep + datas;

		opens = opens + file_sep + maindbname;
		//-------------------------------------------------------------
		dbcon = DatabaseAgent.getConnection(opens, "", "");
		
		dbagent = new DatabaseAgentSupport(dbcon, mainTablePrimaryKey, 
				//mainPrimaryKeyColumnIndex,
				maindbTable);
		nesteddbagent = new DatabaseAgentSupport(dbcon, nestedTablePrimaryKey,
				//nestedPrimaryKeyColumnIndex, 
				nesteddbTable);		
	}

	/**
	 * At start-up, performs a query in database to populate JTables.
	 * 
	 */
	private void performQueryDb() {
		// init select for retrieving columns name
		// and also for populating initial table.

		dbagent.init();

		orderbyS = mainTablePrimaryKey;// when start-up...ID is default!!
		//the init method already perform select!!!
		/*String command = "select * from " + maindbTable + " ORDER BY " + orderbyS + " ASC";
		try {
			dbagent.select(command);
		} catch (Exception ex) {
			ex.printStackTrace();
		}*/

		mainTable = dbagent.getMainTable();
		// allow single selection of rows...not multiple rows!
		ListSelectionModel rowSM = mainTable.getSelectionModel();
		rowSM.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
		rowSM.addListSelectionListener(this);// listener!

		// =========================================================
		// not necessarily needs specific link here...display all data!
		nesteddbagent.init();

		nestedorderbyS = nestedTablePrimaryKey;// when start-up...ID is
												// default!!
		/*command = "select * from " + nesteddbTable + " ORDER BY " + nestedorderbyS + " ASC";
		try {
			nesteddbagent.select(command);
		} catch (Exception ex) {
			ex.printStackTrace();
		}*/

		nestedTable = nesteddbagent.getMainTable();
		ListSelectionModel rowSM2 = nestedTable.getSelectionModel();
		rowSM2.setSelectionMode(ListSelectionModel.MULTIPLE_INTERVAL_SELECTION);//ListSelectionModel.SINGLE_SELECTION);
		rowSM2.addListSelectionListener(this);
	}

	/**
	 * The frame preferred size
	 */
	public Dimension getPreferredSize() {
		return PREFERRED_SIZE;
	}

	/**
	 * Exit application method
	 */
	public void attemptExit() {
		// shut down connection
		//dbagent.shutdown();
		//nesteddbagent.shutdown();
		try{
			if (dbcon != null)
				dbcon.close();
		}catch (Exception e){
			e.printStackTrace();
		}
		
		DatabaseAgent.shutdownDerby();
		
		dispose();
		System.exit(0);
	}

	/**
	 * GUI creation
	 */
	private void createGUI() {
		JPanel content = new JPanel(new BorderLayout());

		JMenuBar menuBar = createMenuBar(resources);
		setJMenuBar(menuBar);

		JPanel mainPanel = createMainPanel();
		content.add(mainPanel, BorderLayout.CENTER);

		setContentPane(new JScrollPane(content));
		content.setOpaque(true); // content panes must be opaque
		pack();
	}

	/**
	 * Main panel GUI creation
	 * @return the main panel
	 */
	private JPanel createMainPanel() {
		Character mnemonic = null;
		JButton button = null;
		JLabel label = null;
		String buttonName = "";
		String buttonToolTip = "";
		String buttonIconName = "";

		orderbyCb = dbagent.getOrderByComboBox();
		orderbyCb.setMaximumRowCount(5);
		orderbyCb.setPreferredSize(sizeCb);
		orderbyCb.addItemListener(this);

		nestedorderbyCb = nesteddbagent.getOrderByComboBox();
		nestedorderbyCb.setMaximumRowCount(5);
		nestedorderbyCb.setPreferredSize(sizeCb);
		nestedorderbyCb.addItemListener(this);

		// =========================MAIN PANEL===========================
		JPanel p1P = new JPanel();
		p1P.setLayout(new FlowLayout(FlowLayout.CENTER, 20, 2));
		label = new JLabel("Username: ");// resources.getString("textfield.userName"));
		label.setForeground(TestClass.foreColor);
		p1P.add(label);
		p1P.add(userNameTf);
		p1P.setBackground(TestClass.bkgColor);

		JPanel p2P = new JPanel();
		p2P.setLayout(new FlowLayout(FlowLayout.CENTER, 20, 2));
		label = new JLabel("Password: ");// resources.getString("textfield.password"));
		label.setForeground(TestClass.foreColor);
		p2P.add(label);
		p2P.add(passwordTf);
		p2P.setBackground(TestClass.bkgColor);

		JPanel p3P = new JPanel();
		p3P.setLayout(new FlowLayout(FlowLayout.CENTER, 20, 2));
		buttonName = "Add";// resources.getString("button.add");
		buttonToolTip = "Add in database";// resources.getString("button.add.toolTip");
		buttonIconName = "/danfulea/resources/add.png";// resources.getString("img.insert");
		button = FrameUtilities.makeButton(buttonIconName, ADD_COMMAND, buttonToolTip, buttonName, this, this);
		mnemonic = new Character('A');// (Character)
										// resources.getObject("button.add.mnemonic");
		button.setMnemonic(mnemonic.charValue());
		p3P.add(button);
		buttonName = "Update";// resources.getString("button.remove");
		buttonToolTip = "Update database";// resources.getString("button.remove.toolTip");
		buttonIconName = "/danfulea/resources/cog.png";// resources.getString("img.close");
		button = FrameUtilities.makeButton(buttonIconName, UPDATE_COMMAND, buttonToolTip, buttonName, this, this);
		mnemonic = new Character('U');// (Character)
										// resources.getObject("button.remove.mnemonic");
		button.setMnemonic(mnemonic.charValue());
		p3P.add(button);
		buttonName = "Clear";// resources.getString("button.remove");
		buttonToolTip = "Clear fields";// resources.getString("button.remove.toolTip");
		buttonIconName = null;// "/jdf/resources/cog.png";//resources.getString("img.close");
		button = FrameUtilities.makeButton(buttonIconName, CLEAR_COMMAND, buttonToolTip, buttonName, this, this);
		mnemonic = new Character('C');// (Character)
										// resources.getObject("button.remove.mnemonic");
		button.setMnemonic(mnemonic.charValue());
		p3P.add(button);
		p3P.setBackground(TestClass.bkgColor);

		JPanel p4P = new JPanel();
		p4P.setLayout(new FlowLayout(FlowLayout.CENTER, 20, 2));
		buttonName = "Remove";// resources.getString("button.remove");
		buttonToolTip = "Remove from database";// resources.getString("button.remove.toolTip");
		buttonIconName = "/danfulea/resources/cross.png";// resources.getString("img.close");
		button = FrameUtilities.makeButton(buttonIconName, DELETE_COMMAND, buttonToolTip, buttonName, this, this);
		mnemonic = new Character('R');// (Character)
										// resources.getObject("button.remove.mnemonic");
		button.setMnemonic(mnemonic.charValue());
		p4P.add(button);
		buttonName = "Remove all";// resources.getString("button.removeAll");
		buttonToolTip = "Remove all from database";// resources.getString("button.removeAll.toolTip");
		buttonIconName = "/danfulea/resources/bin_empty.png";// resources.getString("img.delete.all");
		button = FrameUtilities.makeButton(buttonIconName, DELETEALL_COMMAND, buttonToolTip, buttonName, this, this);
		mnemonic = new Character('e');// (Character)
										// resources.getObject("button.removeAll.mnemonic");
		button.setMnemonic(mnemonic.charValue());
		p4P.add(button);
		label = new JLabel("Records count: ");// resources.getString("textfield.userName"));
		label.setForeground(TestClass.foreColor);
		p4P.add(label);
		// recordsCount.setText(""+rowCount);
		p4P.add(dbagent.getRecordsLabel());// recordsCount);
		p4P.setBackground(TestClass.bkgColor);
		// ===========
		JPanel orderP = new JPanel();
		orderP.setLayout(new FlowLayout(FlowLayout.CENTER, 20, 2));
		label = new JLabel("Sort by: ");
		label.setForeground(TestClass.foreColor);
		orderP.add(label);
		orderP.add(orderbyCb);
		orderP.setBackground(TestClass.bkgColor);

		JPanel suportP = new JPanel(new BorderLayout());
		suportP.setPreferredSize(tableDimension);
		JScrollPane scrollPane = new JScrollPane(mainTable);
		mainTable.setFillsViewportHeight(true);// fill the viewport, never
												// smaller then viewport

		suportP.add(orderP, BorderLayout.NORTH);
		suportP.add(scrollPane, BorderLayout.CENTER);
		suportP.setBackground(TestClass.bkgColor);

		JPanel p5P = new JPanel();
		BoxLayout blp5P = new BoxLayout(p5P, BoxLayout.Y_AXIS);
		p5P.setLayout(blp5P);
		p5P.setBorder(FrameUtilities.getGroupBoxBorder("Records: ", // resources.getString("records.border"),
				TestClass.foreColor));
		p5P.add(suportP);// (scrollPane);//suportSp);
		p5P.add(p4P);
		p5P.setBackground(TestClass.bkgColor);

		JPanel infoBoxP = new JPanel();
		BoxLayout bl03 = new BoxLayout(infoBoxP, BoxLayout.Y_AXIS);
		infoBoxP.setLayout(bl03);
		infoBoxP.add(p1P);
		infoBoxP.add(p2P);
		infoBoxP.add(p3P);
		infoBoxP.add(p5P);
		infoBoxP.setBackground(TestClass.bkgColor);
		infoBoxP.setBorder(FrameUtilities.getGroupBoxBorder("MainTable: ", // resources.getString("records.border"),
				TestClass.foreColor));
		// =========================END MAIN PANEL===========================
		// ====================NESTED PANEL===================================
		JPanel p1Pn = new JPanel();
		p1Pn.setLayout(new FlowLayout(FlowLayout.CENTER, 20, 2));
		label = new JLabel("Email: ");// resources.getString("textfield.userName"));
		label.setForeground(TestClass.foreColor);
		p1Pn.add(label);
		p1Pn.add(emailTf);
		p1Pn.setBackground(TestClass.bkgColor);

		JPanel p2Pn = new JPanel();
		p2Pn.setLayout(new FlowLayout(FlowLayout.CENTER, 20, 2));
		label = new JLabel("Address: ");// resources.getString("textfield.password"));
		label.setForeground(TestClass.foreColor);
		p2Pn.add(label);
		p2Pn.add(addressTf);
		p2Pn.setBackground(TestClass.bkgColor);

		JPanel p3Pn = new JPanel();
		p3Pn.setLayout(new FlowLayout(FlowLayout.CENTER, 20, 2));
		buttonName = "Add";// resources.getString("button.add");
		buttonToolTip = "Add in database";// resources.getString("button.add.toolTip");
		buttonIconName = "/danfulea/resources/add.png";// resources.getString("img.insert");
		button = FrameUtilities.makeButton(buttonIconName, NESTED_ADD_COMMAND, buttonToolTip, buttonName, this, this);
		mnemonic = new Character('A');// (Character)
										// resources.getObject("button.add.mnemonic");
		button.setMnemonic(mnemonic.charValue());
		p3Pn.add(button);
		buttonName = "Update";// resources.getString("button.remove");
		buttonToolTip = "Update database";// resources.getString("button.remove.toolTip");
		buttonIconName = "/danfulea/resources/cog.png";// resources.getString("img.close");
		button = FrameUtilities.makeButton(buttonIconName, NESTED_UPDATE_COMMAND, buttonToolTip, buttonName, this,
				this);
		mnemonic = new Character('U');// (Character)
										// resources.getObject("button.remove.mnemonic");
		button.setMnemonic(mnemonic.charValue());
		p3Pn.add(button);
		buttonName = "Clear";// resources.getString("button.remove");
		buttonToolTip = "Clear fields";// resources.getString("button.remove.toolTip");
		buttonIconName = null;// "/jdf/resources/cog.png";//resources.getString("img.close");
		button = FrameUtilities.makeButton(buttonIconName, NESTED_CLEAR_COMMAND, buttonToolTip, buttonName, this, this);
		mnemonic = new Character('C');// (Character)
										// resources.getObject("button.remove.mnemonic");
		button.setMnemonic(mnemonic.charValue());
		p3Pn.add(button);
		p3Pn.setBackground(TestClass.bkgColor);

		JPanel p4Pn = new JPanel();
		p4P.setLayout(new FlowLayout(FlowLayout.CENTER, 20, 2));
		buttonName = "Remove";// resources.getString("button.remove");
		buttonToolTip = "Remove from database";// resources.getString("button.remove.toolTip");
		buttonIconName = "/danfulea/resources/cross.png";// resources.getString("img.close");
		button = FrameUtilities.makeButton(buttonIconName, NESTED_DELETE_COMMAND, buttonToolTip, buttonName, this,
				this);
		mnemonic = new Character('R');// (Character)
										// resources.getObject("button.remove.mnemonic");
		button.setMnemonic(mnemonic.charValue());
		p4Pn.add(button);
		buttonName = "Remove all";// resources.getString("button.removeAll");
		buttonToolTip = "Remove all from database";// resources.getString("button.removeAll.toolTip");
		buttonIconName = "/danfulea/resources/bin_empty.png";// resources.getString("img.delete.all");
		button = FrameUtilities.makeButton(buttonIconName, NESTED_DELETEALL_COMMAND, buttonToolTip, buttonName, this,
				this);
		mnemonic = new Character('e');// (Character)
										// resources.getObject("button.removeAll.mnemonic");
		button.setMnemonic(mnemonic.charValue());
		p4Pn.add(button);
		label = new JLabel("Records count: ");// resources.getString("textfield.userName"));
		label.setForeground(TestClass.foreColor);
		p4Pn.add(label);
		// recordsCount.setText(""+rowCount);
		p4Pn.add(nesteddbagent.getRecordsLabel());// recordsCount);
		p4Pn.setBackground(TestClass.bkgColor);
		// ===========
		JPanel orderPn = new JPanel();
		orderPn.setLayout(new FlowLayout(FlowLayout.CENTER, 20, 2));
		label = new JLabel("Sort by: ");
		label.setForeground(TestClass.foreColor);
		orderPn.add(label);
		orderPn.add(nestedorderbyCb);
		orderPn.setBackground(TestClass.bkgColor);

		JPanel suportPn = new JPanel(new BorderLayout());
		suportPn.setPreferredSize(tableDimension);
		JScrollPane scrollPanen = new JScrollPane(nestedTable);
		mainTable.setFillsViewportHeight(true);// fill the viewport, never
												// smaller then viewport

		suportPn.add(orderPn, BorderLayout.NORTH);
		suportPn.add(scrollPanen, BorderLayout.CENTER);
		suportPn.setBackground(TestClass.bkgColor);

		JPanel p5Pn = new JPanel();
		BoxLayout blp5Pn = new BoxLayout(p5Pn, BoxLayout.Y_AXIS);
		p5Pn.setLayout(blp5Pn);
		p5Pn.setBorder(FrameUtilities.getGroupBoxBorder("Records: ", // resources.getString("records.border"),
				TestClass.foreColor));
		p5Pn.add(suportPn);// (scrollPane);//suportSp);
		p5Pn.add(p4Pn);
		p5Pn.setBackground(TestClass.bkgColor);

		JPanel infoBoxPn = new JPanel();
		BoxLayout bl03n = new BoxLayout(infoBoxPn, BoxLayout.Y_AXIS);
		infoBoxPn.setLayout(bl03n);
		infoBoxPn.add(p1Pn);
		infoBoxPn.add(p2Pn);
		infoBoxPn.add(p3Pn);
		infoBoxPn.add(p5Pn);
		infoBoxPn.setBackground(TestClass.bkgColor);
		infoBoxPn.setBorder(FrameUtilities.getGroupBoxBorder("Nested Table: ", // resources.getString("records.border"),
				TestClass.foreColor));
		// ===================END NESTED PANEL===============================
		JPanel mainP = new JPanel(new FlowLayout(FlowLayout.CENTER, 20, 2));// new
																			// BorderLayout());
		mainP.add(infoBoxP);// , BorderLayout.CENTER);// main dimension !!
		mainP.add(infoBoxPn);
		mainP.setBackground(TestClass.bkgColor);

		// =-----
		userNameTf.addActionListener(this);
		passwordTf.addActionListener(this);
		// ---
		emailTf.addActionListener(this);
		addressTf.addActionListener(this);

		return mainP;
	}

	/**
	 * Setting up the menu bar.
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
		ImageIcon img;
		String imageName = "";

		// the file menu
		label = resources.getString("menu.file");
		mnemonic = (Character) resources.getObject("menu.file.mnemonic");
		JMenu fileMenu = new JMenu(label, true);
		fileMenu.setMnemonic(mnemonic.charValue());

		imageName = resources.getString("img.close");
		img = FrameUtilities.getImageIcon(imageName, this);
		label = resources.getString("menu.file.exit");
		mnemonic = (Character) resources.getObject("menu.file.exit.mnemonic");
		JMenuItem exitItem = new JMenuItem(label, mnemonic.charValue());
		exitItem.setActionCommand(EXIT_COMMAND);
		exitItem.addActionListener(this);
		exitItem.setIcon(img);
		exitItem.setToolTipText(resources.getString("menu.file.exit.toolTip"));
		fileMenu.add(exitItem);

		// the help menu
		label = resources.getString("menu.help");
		mnemonic = (Character) resources.getObject("menu.help.mnemonic");
		JMenu helpMenu = new JMenu(label);
		helpMenu.setMnemonic(mnemonic.charValue());

		imageName = resources.getString("img.about");
		img = FrameUtilities.getImageIcon(imageName, this);
		label = resources.getString("menu.help.about");
		mnemonic = (Character) resources.getObject("menu.help.about.mnemonic");
		JMenuItem aboutItem = new JMenuItem(label, mnemonic.charValue());
		aboutItem.setActionCommand(ABOUT_COMMAND);
		aboutItem.addActionListener(this);
		aboutItem.setIcon(img);
		aboutItem.setToolTipText(resources.getString("menu.help.about.toolTip"));

		label = resources.getString("menu.help.LF");
		mnemonic = (Character) resources.getObject("menu.help.LF.mnemonic");
		JMenuItem lfItem = new JMenuItem(label, mnemonic.charValue());
		lfItem.setActionCommand(LOOKANDFEEL_COMMAND);
		lfItem.addActionListener(this);
		lfItem.setToolTipText(resources.getString("menu.help.LF.toolTip"));

		helpMenu.add(lfItem);
		helpMenu.addSeparator();
		helpMenu.add(aboutItem);

		// finally, glue together the menu and return it
		menuBar.add(fileMenu);
		menuBar.add(helpMenu);

		return menuBar;
	}

	/**
	 * Sets the links (column name link and its corresponding value) for nested table
	 */
	private void setLinks() {
		int selID = 0;
		int selRow = mainTable.getSelectedRow();
		if (selRow != -1) {
			selID = (Integer) mainTable.getValueAt(selRow, dbagent.getPrimaryKeyColumnIndex());//mainPrimaryKeyColumnIndex);// first
																						// column,0,
																						// is
																						// ID
		} else {
			return;
		}
		nesteddbagent.setLinks(IDlink, Convertor.intToString(selID));
	}

	/**
	 * Deletes a row from main table and its corresponding row from nested table
	 */
	private void delete() {
		int selID = 0;
		int selRow = mainTable.getSelectedRow();
		if (selRow != -1) {
			selID = (Integer) mainTable.getValueAt(selRow, dbagent.getPrimaryKeyColumnIndex());//mainPrimaryKeyColumnIndex);// first
																						// column,0,
																						// is
																						// ID
		} else {
			return;
		}

		// deleteRecord
		dbagent.delete(Convertor.intToString(selID));//, orderbyS);
		
		// now delete from nested
		nesteddbagent.setLinks(IDlink, Convertor.intToString(selID));
		nesteddbagent.delete(nesteddbTable, IDlink, Convertor.intToString(selID));
		
		dbagent.performSelection(orderbyS);
	}

	/**
	 * Deletes all data from both main and nested tables
	 */
	private void deleteAll() {
		// deleteAllRecords
		dbagent.deleteAll();//orderbyS);
		dbagent.performSelection(orderbyS);
		// delete from nested
		nesteddbagent.deleteAll();//nestedorderbyS);
		nesteddbagent.performSelection(nestedorderbyS);
	}

	/**
	 * Sorts data from main table
	 */
	private void sort() {
		orderbyS = (String) orderbyCb.getSelectedItem();
		// performSelection();
		dbagent.performSelection(orderbyS);
	}

	/**
	 * Sorts data from nested table
	 */
	private void sort2() {
		nestedorderbyS = (String) nestedorderbyCb.getSelectedItem();
		// link to main table
		setLinks();
		// end link to main table
		nesteddbagent.performSelection(nestedorderbyS);
	}

	/**
	 * Saves data in main table
	 */
	private void add() {
		String[] data = new String[2];
		data[0] = userNameTf.getText();// String username =
										// userNameTf.getText();//"johndoe2";
		data[1] = passwordTf.getText();// String password =
										// passwordTf.getText();//"1232";
		// test for maximum varchar!
		if (data[0].length() > 49 || data[1].length() > 49) {
			JOptionPane.showMessageDialog(this, "Too long text. Maximum 50 characters!", "Error",
					JOptionPane.ERROR_MESSAGE);
			return;
		}
		// ---------------------
		// insert new record and display all data order by a specific column
		dbagent.insert(data);//, orderbyS);
		dbagent.performSelection(orderbyS);
	}

	/**
	 * Updates data in main table
	 */
	private void update() {
		String[] data = new String[2];
		data[0] = userNameTf.getText();
		data[1] = passwordTf.getText();
		// --some testing
		if (data[0].length() > 49 || data[1].length() > 49) {
			JOptionPane.showMessageDialog(this, "Too long text. Maximum 50 characters!", "Error",
					JOptionPane.ERROR_MESSAGE);
			return;
		}
		// --------------------------
		int selID = 0;
		int selRow = mainTable.getSelectedRow();
		if (selRow != -1) {
			selID = (Integer) mainTable.getValueAt(selRow, dbagent.getPrimaryKeyColumnIndex());//mainPrimaryKeyColumnIndex);// first
																						// column,0,
																						// is
																						// ID
		} else {
			return;
		}

		dbagent.setSelectedRow(selRow);
		dbagent.update(data, Convertor.intToString(selID));//, orderbyS);
		dbagent.performSelection(orderbyS);
	}

	/**
	 * Clear some GUI fields
	 */
	private void clear() {
		userNameTf.setText("");
		passwordTf.setText("");
		userNameTf.requestFocusInWindow();

		// String command = "select * from "+maindbTable+" where username like
		// 'as'";
		// dbagent.select(command);//works!
	}
	// ================

	/**
	 * Saves data in nested table
	 */
	private void add2() {
		int selID = 0;
		int selRow = mainTable.getSelectedRow();
		if (selRow != -1) {
			selID = (Integer) mainTable.getValueAt(selRow, dbagent.getPrimaryKeyColumnIndex());//mainPrimaryKeyColumnIndex);// first
																						// column,0,
																						// is
																						// ID
		} else {
			return;
		}

		String[] data = new String[3];
		data[0] = Convertor.intToString(selID);
		data[1] = emailTf.getText();// String username =
									// userNameTf.getText();//"johndoe2";
		data[2] = addressTf.getText();// String password =
										// passwordTf.getText();//"1232";
		// test for maximum varchar!
		if (data[1].length() > 49 || data[2].length() > 49) {
			JOptionPane.showMessageDialog(this, "Too long text. Maximum 50 characters!", "Error",
					JOptionPane.ERROR_MESSAGE);
			return;
		}
		// ---------------------
		/// setLinks();
		nesteddbagent.setLinks(IDlink, Convertor.intToString(selID));
		// insert new record and display all data order by a specific column
		nesteddbagent.insert(data);//, nestedorderbyS);
		nesteddbagent.performSelection(nestedorderbyS);

	}

	/**
	 * Updates data in nested table
	 */
	private void update2() {
		int selID = 0;
		int selRow = mainTable.getSelectedRow();
		if (selRow != -1) {
			selID = (Integer) mainTable.getValueAt(selRow, dbagent.getPrimaryKeyColumnIndex());//mainPrimaryKeyColumnIndex);// first
																						// column,0,
																						// is
																						// ID
		} else {
			return;
		}

		String[] data = new String[3];
		data[0] = Convertor.intToString(selID);
		data[1] = emailTf.getText();
		data[2] = addressTf.getText();
		// --some testing
		if (data[0].length() > 49 || data[1].length() > 49) {
			JOptionPane.showMessageDialog(this, "Too long text. Maximum 50 characters!", "Error",
					JOptionPane.ERROR_MESSAGE);
			return;
		}
		// --------------------------
		int selID2 = 0;
		int selRow2 = nestedTable.getSelectedRow();
		if (selRow2 != -1) {
			selID2 = (Integer) nestedTable.getValueAt(selRow2, nesteddbagent.getPrimaryKeyColumnIndex());//nestedPrimaryKeyColumnIndex);// first
																							// column,0,
																							// is
																							// ID
		} else {
			return;
		}

		// setLinks();
		nesteddbagent.setLinks(IDlink, Convertor.intToString(selID));
		nesteddbagent.setSelectedRow(selRow2);
		nesteddbagent.update(data, Convertor.intToString(selID2));//, nestedorderbyS);
		nesteddbagent.performSelection(nestedorderbyS);
	}

	/**
	 * Clear some GUI fields
	 */
	private void clear2() {
		emailTf.setText("");
		addressTf.setText("");
		emailTf.requestFocusInWindow();
	}

	/**
	 * Delete a row from nested table
	 */
	private void delete2() {
		/*int selID = 0;
		int selRow = nestedTable.getSelectedRow();
		if (selRow != -1) {
			selID = (Integer) nestedTable.getValueAt(selRow, nestedPrimaryKeyColumnIndex);// first
																							// column,0,
																							// is
																							// ID
		} else {
			return;
		}*/

		
		int[] selRow = nestedTable.getSelectedRows();
		int n = selRow.length;
		String[] selID = new String[n];
		
		for (int i = 0; i<n; i++){
			if (selRow[i] != -1) {
				int ii = (Integer) nestedTable.getValueAt(selRow[i], nesteddbagent.getPrimaryKeyColumnIndex());//nestedPrimaryKeyColumnIndex);
				selID[i] = Convertor.intToString(ii);
			}
		}
		
		// deleteRecord
		setLinks();
		nesteddbagent.delete(selID);//, nestedorderbyS);
		nesteddbagent.performSelection(nestedorderbyS);
	}

	/**
	 * Delete all data from nested table. This is only shown here for DEMO purpose. It should never be used in real application since delete all must only came from main table.
	 */
	private void deleteAll2() {
		setLinks();
		nesteddbagent.deleteAll();//nestedorderbyS);
		nesteddbagent.performSelection(nestedorderbyS);
	}
	// ===============

	/**
	 * About frame
	 */
	private void about() {
		new AboutFrame(resources);
	}

	/**
	 * Sets the look and feel
	 */
	private void lookAndFeel() {
		 setVisible(false);// setEnabled(false);
		 new ScanDiskLFGui(this);
	}
	
	/**
	 * Main actions are set here
	 */
	public void actionPerformed(ActionEvent e) {
		// TODO Auto-generated method stub
		command = e.getActionCommand();
		if (command.equals(ABOUT_COMMAND)) {
			about();
		} else if (command.equals(EXIT_COMMAND)) {
			attemptExit();
		} else if (command.equals(LOOKANDFEEL_COMMAND)) {
			lookAndFeel();
		} else if (e.getSource() == passwordTf || command.equals(ADD_COMMAND)) {
			add();
			userNameTf.setText("");
			passwordTf.setText("");
			userNameTf.requestFocusInWindow();
		} else if (command.equals(UPDATE_COMMAND)) {
			update();
		} else if (command.equals(CLEAR_COMMAND)) {
			clear();
		} else if (command.equals(DELETE_COMMAND)) {
			delete();
		} else if (command.equals(DELETEALL_COMMAND)) {
			deleteAll();
		} else if (e.getSource() == userNameTf) {// press enter
			passwordTf.setText("");
			passwordTf.requestFocusInWindow();
		}
		// ==============
		else if (e.getSource() == addressTf || command.equals(NESTED_ADD_COMMAND)) {
			add2();
			emailTf.setText("");
			addressTf.setText("");
			emailTf.requestFocusInWindow();
		} else if (command.equals(NESTED_UPDATE_COMMAND)) {
			update2();
		} else if (command.equals(NESTED_CLEAR_COMMAND)) {
			clear2();
		} else if (command.equals(NESTED_DELETE_COMMAND)) {
			delete2();
		} else if (command.equals(NESTED_DELETEALL_COMMAND)) {
			deleteAll2();
		} else if (e.getSource() == emailTf) {// press enter
			addressTf.setText("");
			addressTf.requestFocusInWindow();
		}
	}
	
	/**
	 * JCombobox actions are set here
	 */
	public void itemStateChanged(ItemEvent e) {
		// TODO Auto-generated method stub
		if (e.getSource() == orderbyCb) {
			sort();
		} else if (e.getSource() == nestedorderbyCb) {
			sort2();
		}
	}
	
	/**
	 * JTable related actions are set here
	 */
	public void valueChanged(ListSelectionEvent e) {

		if (e.getSource() == mainTable.getSelectionModel()) {
			// firts reset fields
			userNameTf.setText("");
			passwordTf.setText("");
			emailTf.setText("");
			addressTf.setText("");

			int selRow = mainTable.getSelectedRow();
			if (selRow != -1) {
				int selID = (Integer) mainTable.getValueAt(selRow, dbagent.getPrimaryKeyColumnIndex());//mainPrimaryKeyColumnIndex);// first
																								// column,0,
																								// is
																								// ID
				userNameTf.setText((String) mainTable.getValueAt(selRow, 1));
				passwordTf.setText((String) mainTable.getValueAt(selRow, 2));

				// ===update nested===
				nesteddbagent.setLinks(IDlink, Convertor.intToString(selID));
				nesteddbagent.performSelection(nestedorderbyS);
			} else {
				return;
			}
		} else if (e.getSource() == nestedTable.getSelectionModel()) {
			emailTf.setText("");
			addressTf.setText("");
			int selRow = nestedTable.getSelectedRow();
			if (selRow != -1) {
				emailTf.setText((String) nestedTable.getValueAt(selRow, 2));
				addressTf.setText((String) nestedTable.getValueAt(selRow, 3));
			} else {
				return;
			}
		}
		// ---------------------------------------------------
	}

	public static void main(String[] args) {
		LookAndFeel.loadLookAndFeel();//load LF if available!
		new TestClass();
		//weve just created testTable3
		
		//String dbTable="testTable";
		//String dbName="testDb";	
		//======Database creation==
		//DatabaseAgent.dataFolderUsedByDerbyEmbedded = "Data";
		//DatabaseAgent.createDerbyDatabaseInDataFolder("TESTARE","","");//works
		//======END================
				
		//======Database delete============
		//DatabaseAgent.ID_CONNECTION = DatabaseAgent.POSTGRESQL_CONNECTION;
		//DatabaseAgent.deleteDatabase("testare", "postgres", "admin");
		//works only if in PostgreSQL we created db with lower-case: eg testare and not TESTARE
		//for derby simply delete folder manually!!!
				
		//DatabaseAgent.ID_CONNECTION = DatabaseAgent.MYSQL_CONNECTION;
		//DatabaseAgent.deleteDatabase("TEST", "root", "");//works
		//DatabaseAgent.deleteDatabase("TEST1", "root", "");//works regardless of lower-case
		//in fact in MySQL -phpAdmin, database is created always with lower case!!
		//==========END===================
				
		//========Table creation===============
		//String table = "testtable2";
		//String [] columns= {"nrcrt", "id", "email", "address"};
		//String [] types= {"integer", "integer", "VARCHAR(59)", "VARCHAR(59)"};
		//String [] indexes= {"Y", "N", "N", "N"};
				
		//String table = "testtable";
		//String [] columns= {"id", "username", "password"};
		//String [] types= {"integer", "VARCHAR(59)", "VARCHAR(59)"};
		//String [] indexes= {"Y", "N", "N"};
				
		//no indexes=>we will mimic them
		//String table = "testtable3";
		//String [] columns= {"id", "Name", "Surname"};
		//String [] types= {"integer", "VARCHAR(59)", "VARCHAR(59)"};//TEXT is Postqres CLOB
		//String [] indexes= {"N", "N", "N", "N"};
				
		//DatabaseAgent.ID_CONNECTION = DatabaseAgent.DERBY_CONNECTION;
		//DatabaseAgent.ID_CONNECTION = DatabaseAgent.MYSQL_CONNECTION;
		//DatabaseAgent.ID_CONNECTION = DatabaseAgent.POSTGRESQL_CONNECTION;
		//DatabaseAgent.createTable(dbName, table, columns, types, indexes, "", "");//works
		//DatabaseAgent.createTable(dbName, table, columns, types, indexes, "root", "");//works
		//DatabaseAgent.createTable(dbName, table, columns, types, indexes, "postgres", "admin");//works
				
		//========end==========================
				
		//==========Table drop====================
		//DatabaseAgent.deleteTable(dbName, table, "", "");//works
		//DatabaseAgent.deleteTable(dbName, table, "root", "");//works
		//DatabaseAgent.deleteTable(dbName, table, "postgres", "admin");//works
		//==========END=============================
		
		//============RESTART PK===================
		//DatabaseAgent.restartPKincrementation(dbName, dbTable, "",
		//		"", "ID", 118);//WORKS!!!
		
		//java reflection
		//Method
		/*Field[] flds = TestClass.class.getDeclaredFields();
		for (Field f:flds){
			System.out.println("field = "+f.getName());//works
		}*/
		
		/*try {
			String mydate = "2012-07-11 10:55:21.7878";
			long milis = Convertor.getMillisFromDateTime(mydate);
			Timestamp t = new Timestamp(milis);
			System.out.println("Milis: "+milis);
			//System.out.println(System.currentTimeMillis());
			System.out.println(Convertor.timestampToString(t));//ok!!!!
			
			Timestamp tt = Convertor.stringToTimestamp(mydate);
			System.out.println(Convertor.timestampToString(tt));
			
			Timestamp ttt = Convertor.stringToTimestamp(Convertor.timestampToString(tt));
			System.out.println(Convertor.timestampToString(ttt));
			
			mydate="2012-07-11";// 00:00:00";
			milis = Convertor.getMillisFromDateTime(mydate+" 00:00:00");//time is not important here
			System.out.println("Milis: "+milis);
			Date d = new Date(milis);
			System.out.println(Convertor.dateToString(d));
			
			Date dd = Convertor.stringToDate(mydate);
			System.out.println(Convertor.dateToString(dd));
			
			
		} catch (ParseException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		String bdS = "8.765.478,980";///fail
		bdS = "7865324.5";//OK
		//bdS = "7865324";//OK
		//BigDecimal bbb = new BigDecimal(7865324.5);
		BigDecimal bbb = Convertor.stringToBigDecimal(bdS);
		System.out.println(Convertor.bigDecimalToString(bbb));
		*/
		
		//int bb = 0xfe;//254 int literal
		/*int my23 = 0x23;//0010 0011 = 2^5+2+1=35
		System.out.println("my= "+my23);
		byte bb = new Integer(0x23).byteValue();//-2 byte
		byte xff = new Integer(0xff).byteValue();//-1 byte
		System.out.println("Byte bb= "+bb+"; "+xff);
		int promo = bb;//-2
		System.out.println("int promoted bb= "+promo);
		promo = bb & 0xff;//real int promotion 0xff is int literal
		System.out.println("int bb&0xff= "+promo);//254 again!!*/
		
		//so int promotion of byte fe = 1111 1110 in binary or 0xfe in hexadecimal
		//makes it: ff ff ff fe (8x4=32 bits) this is a negative number because
		//the left-most bit is 1 so, the algorithm for negative numbers requires
		//to reveres all bits and adding 1=>00 00 00 01 and plus 1 gives 00 00 00 10 which is 2,
		//thus -2...which is consistent to original byte! OK!!!!!!!!!!!!
		//if we apply the AND with 0xff int which is 00 00 00 ff we obtain:
		//00 00 00 fe which is 254! OK!
		//promotion are done in the following way:
		//read left-most bite in byte and fill the bytes in int from the left!!!
		//if byte is fe=>put 1=> ff ff ff fe
		//if byte is 23=>put 0=> 00 00 00 23
		//so the operation bb & 0xff is necesary to handle NEGATIVE bytes
		//byte bb = new Integer(0x23).byteValue();//-2 byte==35 in int
		//String sst = Integer.toString((bb & 0xff) + 0x100, 16).substring(1);
		//System.out.println("0x23&0xff++ 0x10016base substring1= "+sst);
		//sst = Integer.toString((bb & 0xff) + 0x100, 16);
		//System.out.println("0x23&0xff++ 0x10016base= "+sst);
		//String ss= "mama.mia";
		//System.out.println(FileOperation.getFileExtension(ss)+" "+
		//FileOperation.getFileExtensionLength(ss));
	}

}
