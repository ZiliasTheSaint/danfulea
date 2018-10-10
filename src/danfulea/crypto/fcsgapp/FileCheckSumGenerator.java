package danfulea.crypto.fcsgapp;

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
import java.io.FileWriter;
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
import javax.swing.JTextField;
import javax.swing.JToolBar;

import danfulea.crypto.ChecksumUtilities;
import danfulea.utils.AboutFrame;
import danfulea.utils.ExampleFileFilter;
import danfulea.utils.FrameUtilities;
import danfulea.utils.LookAndFeel;

/**
* File MD5SUM generator and Files comparison. To create a jar file, choose: 
* this java file, its resources (java + image file), ChecksumUtilities from danfulea.crypto, 
* Convertor.java from danfulea.math package, DanfuleaResources form danfulea.resources (images not needed), 
* FrameUtilities, AboutFrame, LookAndFeel and ExampleFileFilter from danfulea.utils. As manifest choose the file from 
* META-inf-fcsg.
* 
* @author Dan Fulea, 13 JUL. 2010.
*

*/

public class FileCheckSumGenerator extends JFrame implements ActionListener, Runnable
{	
	private static final long serialVersionUID = 1L;
	
	private Thread computationTh = null;// computation thread!
	private Thread statusTh = null;// status display thread!
	private int delay = 100;
	private int frameNumber = -1;
	private String statusRunS = "";
	
    private static final String EXIT_COMMAND = "EXIT";
    private static final String ABOUT_COMMAND = "ABOUT";
    private static final String OPEN_COMMAND = "OPEN";
    private static final String OPEN2_COMMAND = "OPEN2";
    private static final String GENERATE_COMMAND = "GENERATE";
    private static final String COMPARE_COMMAND = "COMPARE";
    private String command="";

	private Color bkgColor = new Color(230, 255, 210, 255);//new Color(180, 220, 150, 255);//new Color(230, 255, 210, 255);
	private static final String BASE_RESOURCE_CLASS="danfulea.crypto.fcsgapp.resources.FileCheckSumGeneratorResources";
	private ResourceBundle resources;
	private JLabel statusL= new JLabel("Waiting...");
		
    private static final Dimension PREFERRED_SIZE = new Dimension(900, 500);

    private JRadioButton anotherFileRb,md5sumFileRb;
    private JButton openB,open2B,generateB,compareB;
    private JTextArea textArea= new JTextArea();
    private JTextField openTf=new JTextField(35);
    private JTextField open2Tf=new JTextField(35);

    private boolean checkSuccess=true;
    private String checkError="";
        
	//private boolean stopAppend =false;
	//private boolean stopAnim=true;

   /**
    * Constructor
    */
	public FileCheckSumGenerator()
	{
        super("FileCheckSumGenerator");
        resources = ResourceBundle.getBundle(BASE_RESOURCE_CLASS);

		//the key to force attemptExit() method on close!!
        setDefaultCloseOperation(JFrame.DO_NOTHING_ON_CLOSE);
        //----------------
        addWindowListener(new WindowAdapter()
        {
            public void windowClosing(WindowEvent e)
            {
                //attemptExit();
                dispose();
                System.exit(0);
            }
        });

        JMenuBar menuBar = createMenuBar(resources);
        setJMenuBar(menuBar);

        createGUI();
        setDefaultLookAndFeelDecorated(true);
		FrameUtilities.createImageIcon(this.resources.getString("form.icon.url"),this);
		FrameUtilities.centerFrameOnScreen(this);
	    setVisible(true);


	}

	/**
	 * Setting up window size
	 */
    public Dimension getPreferredSize()
    {
        return PREFERRED_SIZE;
    }

    /**
     * GUI creation
     */
    private void createGUI()
    {
		JPanel content = new JPanel(new BorderLayout());
        JPanel mainPanel = createMainPanel();
        content.add(mainPanel);
        //Create the statusbar.
        JToolBar statusBar = new JToolBar();
        statusBar.setFloatable(false);
        initStatusBar(statusBar);
        content.add(statusBar, BorderLayout.PAGE_END);

        setContentPane(content);
        content.setOpaque(true); //content panes must be opaque
        pack();
	}

    /**
     * Main panel creation
     * @return main panel
     */
    private JPanel createMainPanel()
    {
		anotherFileRb=new JRadioButton(resources.getString("another.file"));
		md5sumFileRb=new JRadioButton(resources.getString("md5.file"));
		anotherFileRb.setOpaque(true);
		md5sumFileRb.setOpaque(true);

		openB=new JButton(resources.getString("menu.file.open"));
		Character mnemonic=(Character) resources.getObject("menu.file.open.mnemonic");
		openB.setMnemonic(mnemonic.charValue());
		openB.setActionCommand(OPEN_COMMAND);
		openB.addActionListener(this);
		openB.setToolTipText(resources.getString("menu.file.open.toolTip"));
		openB.setIcon(FrameUtilities.getImageIcon(resources.getString("open.icon"),this));

		compareB=new JButton(resources.getString("menu.file.compare"));
		mnemonic = (Character) resources.getObject("menu.file.compare.mnemonic");
		compareB.setMnemonic(mnemonic.charValue());
		compareB.setActionCommand(COMPARE_COMMAND);
		compareB.addActionListener(this);
		compareB.setToolTipText(resources.getString("menu.file.compare.toolTip"));
		compareB.setIcon(FrameUtilities.getImageIcon(resources.getString("compare.icon"),this));

		generateB=new JButton(resources.getString("menu.file.generate"));
		mnemonic = (Character) resources.getObject("menu.file.generate.mnemonic");
		generateB.setMnemonic(mnemonic.charValue());
		generateB.setActionCommand(GENERATE_COMMAND);
		generateB.addActionListener(this);
		generateB.setToolTipText(resources.getString("menu.file.generate.toolTip"));
		generateB.setIcon(FrameUtilities.getImageIcon(resources.getString("generate.icon"),this));

		open2B=new JButton(resources.getString("menu.file.open2"));
		mnemonic = (Character) resources.getObject("menu.file.open2.mnemonic");
		open2B.setMnemonic(mnemonic.charValue());
		open2B.setActionCommand(OPEN2_COMMAND);
		open2B.addActionListener(this);
		open2B.setToolTipText(resources.getString("menu.file.open2.toolTip"));
		open2B.setIcon(FrameUtilities.getImageIcon(resources.getString("open.icon"),this));

		open2Tf.setEnabled(false);
		openTf.setEnabled(false);

		ButtonGroup group = new ButtonGroup();
		group.add(anotherFileRb);
		group.add(md5sumFileRb);
		//Put the radio buttons in a row in a panel.
	    JPanel buttP = new JPanel();
		buttP.setLayout(new FlowLayout(FlowLayout.CENTER, 20,2));
		buttP.add(anotherFileRb);
		buttP.add(md5sumFileRb);
	    anotherFileRb.setBackground(bkgColor);
	    md5sumFileRb.setBackground(bkgColor);
		buttP.setBackground(bkgColor);
		md5sumFileRb.setSelected(true);

		JPanel openP = new JPanel();
		openP.setLayout(new FlowLayout(FlowLayout.CENTER, 20,2));
		openP.add(openTf);
		openP.add(openB);
		openP.add(generateB);
		openP.setBackground(bkgColor);
		openP.setBorder(FrameUtilities.getGroupBoxBorder(resources.getString("titleBorder.open")));

		JPanel open2P = new JPanel();
		open2P.setLayout(new FlowLayout(FlowLayout.CENTER, 20,2));
		open2P.add(open2Tf);
		open2P.add(open2B);
		open2P.add(compareB);
		open2P.setBackground(bkgColor);

	    JPanel op2P=new JPanel();
	    BoxLayout bld2 = new BoxLayout(op2P,BoxLayout.Y_AXIS);
    	op2P.setLayout(bld2);
    	op2P.add(buttP);
    	op2P.add(open2P);
    	op2P.setBackground(bkgColor);
    	op2P.setBorder(FrameUtilities.getGroupBoxBorder(resources.getString("titleBorder.open2")));

        textArea.setCaretPosition(0);
        textArea.setEditable(false);
        textArea.setText("");
        textArea.setLineWrap(true);
        textArea.setWrapStyleWord(true);

	    JPanel resultP=new JPanel(new BorderLayout());
	    resultP.add(new JScrollPane(textArea),BorderLayout.CENTER);
	    resultP.setBackground(bkgColor);

	    JPanel p1=new JPanel();
	    p1.setLayout(new FlowLayout());
	    JLabel label=new JLabel();
	    label.setText(resources.getString("mainPanel.textArea.label"));
	    p1.add(label);
	    p1.setBackground(bkgColor);

	    JPanel p2=new JPanel();
	    BoxLayout bld = new BoxLayout(p2,BoxLayout.Y_AXIS);
    	p2.setLayout(bld);
    	p2.add(openP, null);
    	p2.add(op2P, null);
    	p2.add(p1, null);
	    p2.setBackground(bkgColor);

        JPanel mainP=new JPanel(new BorderLayout());
        mainP.add(p2,BorderLayout.NORTH);
        mainP.add(resultP,BorderLayout.CENTER);
	    mainP.setBackground(bkgColor);
        return mainP;
    }

    /**
     * Status bar initialization
     * @param toolBar
     * the status bar
     */
    private void initStatusBar(JToolBar toolBar)
    {
        JPanel toolP = new JPanel();
        toolP.setLayout(new FlowLayout(FlowLayout.LEFT,5,1));

        toolP.add(statusL);
        toolBar.add(toolP) ;
        statusL.setText(resources.getString("status.wait"));
    }

    /**
     * Basic actions are performed here
     */
    public void actionPerformed(ActionEvent evt)
    {
        command = evt.getActionCommand();

        if (command.equals(ABOUT_COMMAND))
        {
			about();
		}
        else if (command.equals(EXIT_COMMAND))
        {
		    attemptExit();
		}
		else if (command.equals(OPEN_COMMAND))
		{
			openFile();
		}
		else if (command.equals(OPEN2_COMMAND))
		{
			openFile2();
		}
		else if (command.equals(GENERATE_COMMAND))
		{
			
			statusRunS = resources.getString("status.computing");
			startThread();
		}
		else if (command.equals(COMPARE_COMMAND))
		{
			
			statusRunS = resources.getString("status.computing");
			startThread();

		}
    }

    /**
     * Display about window
     */
    private void about()
    {
		new AboutFrame(resources);//this);
	}

    /**
     * Exit application
     */
    private void attemptExit()
    {
        String title = resources.getString("dialog.exit.title");
        String message = resources.getString("dialog.exit.message");

        Object[] options = (Object[])resources.getObject("dialog.exit.buttons");
		int result = JOptionPane.showOptionDialog(this,message, title,
		JOptionPane.YES_NO_OPTION, JOptionPane.QUESTION_MESSAGE,
		null, options, options[0]);
        if (result == JOptionPane.YES_OPTION)
        {
            dispose();
            System.exit(0);
        }
	}

    /**
     * Open main file
     */
    private void openFile()
    {
		String currentDir=System.getProperty("user.dir");
        JFileChooser chooser = new JFileChooser(currentDir);
        chooser.setFileSelectionMode(JFileChooser.FILES_ONLY);

		int returnVal = chooser.showOpenDialog(this);//parent=this frame
		if(returnVal == JFileChooser.APPROVE_OPTION)
		{
			chooser.getSelectedFile();
			String infilename= chooser.getSelectedFile().toString();
			openTf.setText(infilename);
			statusL.setText(resources.getString("status.open")+infilename);
		}

	}

    /**
     * Open second file
     */
    private void openFile2()
    {
		String currentDir=System.getProperty("user.dir");
        JFileChooser chooser = new JFileChooser(currentDir);
        chooser.setFileSelectionMode(JFileChooser.FILES_ONLY);

		int returnVal = chooser.showOpenDialog(this);//parent=this frame
		if(returnVal == JFileChooser.APPROVE_OPTION)
		{
			chooser.getSelectedFile();
			String infilename= chooser.getSelectedFile().toString();
			open2Tf.setText(infilename);
			statusL.setText(resources.getString("status.open2")+infilename);
		}
	}

    /**
     * Try to generate the MD5 checksum
     */
    private void generate()
    {
		String infileS=openTf.getText();
		if (infileS.equals(""))
		{
	        textArea.selectAll();
	        textArea.replaceSelection("");
	        String str=resources.getString("open.fail");
			textArea.append(str);
			statusL.setText(resources.getString("status.open.fail"));
	        //pw.stopAnimation();
	        //pw=null;
			stopThread();

			return;
		}

		String ext=resources.getString("file.extension");
		String pct=".";
		String description=resources.getString("file.description");
		ExampleFileFilter jpgFilter =new ExampleFileFilter(ext, description);
	    String filename="";
	    String currentDir=System.getProperty("user.dir");
	    System.getProperty("file.separator");
	    String opens=currentDir;
		JFileChooser chooser = new JFileChooser(opens);
		chooser.addChoosableFileFilter(jpgFilter);
		chooser.setFileSelectionMode(JFileChooser.FILES_ONLY);

		int returnVal = chooser.showSaveDialog(this);//parent=this frame
		if(returnVal == JFileChooser.APPROVE_OPTION)
		{
			File outfile=chooser.getSelectedFile();
			filename= chooser.getSelectedFile().toString();
			int fl=filename.length();
			String test=filename.substring(fl-4);//exstension lookup!!
			String ctest=pct+ext;
			if (test.compareTo(ctest)!=0)
				filename=chooser.getSelectedFile().toString()+pct+ext;

			if (outfile.exists())
			{
				String title = resources.getString("dialog.overwrite.title");
        		String message = resources.getString("dialog.overwrite.message");

        		Object[] options = (Object[])resources.getObject("dialog.overwrite.buttons");
				int result = JOptionPane.showOptionDialog(this,message, title,
				JOptionPane.YES_NO_OPTION, JOptionPane.QUESTION_MESSAGE,
				null, options, options[0]);
        		if (result != JOptionPane.YES_OPTION)
        		{
			        
        			stopThread();
        			statusL.setText(resources.getString("status.done"));
					return;
				}

			}

			File infile=new File(infileS);
			String s="";
			s=checkSum(infile);//key method
			if(!checkSuccess)
			{
	       		textArea.selectAll();
	       		textArea.replaceSelection("");
	       		String str=resources.getString("error")+checkError+"\n";
				textArea.append(str);
				statusL.setText(resources.getString("status.error"));
	       		//pw.stopAnimation();
	       		//pw=null;
				stopThread();
				
				return;
			}

	        textArea.selectAll();
	        textArea.replaceSelection("");
	        String str=resources.getString("md5sum.main.display")+"\n";
			textArea.append(str);
			textArea.append(s);

			try
			{
				FileWriter sigfos = new FileWriter(filename);
				sigfos.write(s+" ");//final string whitespace!
				sigfos.close();
				statusL.setText(resources.getString("status.save")+filename);
			}
			catch (Exception ex)
			{
	        	textArea.selectAll();
	        	textArea.replaceSelection("");
	        	str=resources.getString("error")+ex.toString()+"\n";
				textArea.append(str);
				statusL.setText(resources.getString("status.error"));
	        	
				stopThread();

				return;
			}
		}

		
		stopThread();
	}

    /**
     * Try comparing files
     */
    private void compare()
    {
		String checksum1="";
		String checksum2="";

		String str="";
		String infileS=openTf.getText();
		if (infileS.equals(""))
		{
	        textArea.selectAll();
	        textArea.replaceSelection("");
	        str=resources.getString("open.fail")+"\n";
			textArea.append(str);
			statusL.setText(resources.getString("status.open.fail"));
	       
			stopThread();

			return;
		}

		String infile2S=open2Tf.getText();
		if (infile2S.equals(""))
		{
	        textArea.selectAll();
	        textArea.replaceSelection("");
	        str=str+resources.getString("open2.fail")+"\n";
			textArea.append(str);
			statusL.setText(resources.getString("status.open.fail"));
	      
			stopThread();

			return;
		}

		if (md5sumFileRb.isSelected())
		{
			String ext=resources.getString("file.extension");
			String pct=".";
	    	String filename="";

			int i =0;
			//@SuppressWarnings("unused")
			//int lnr =0;//line number
			StringBuffer desc=new StringBuffer();
			boolean haveData=false;
			//--------------
			filename= infile2S;
			int fl=filename.length();
			String test=filename.substring(fl-4);//exstension lookup!!
			String ctest=pct+ext;
			if (test.compareTo(ctest)!=0)
				filename=infile2S+pct+ext;

			try
			{
				FileInputStream in = new FileInputStream(filename);
        	   	while ((i = in.read()) != -1)
        	   	{
					if (!Character.isWhitespace((char)i))
					{
			   			desc.append((char)i);
			   			haveData=true;
					}
					else
					{
						if (haveData)
						{
							haveData=false;//reset
							//lnr++;
							String s=desc.toString();
							checksum2=s;
						}//have data
						desc=new StringBuffer();
					}
        	   	}
				in.close();
			}
			catch (Exception e)
			{
	        	textArea.selectAll();
	        	textArea.replaceSelection("");
	        	str=resources.getString("error")+e.toString()+"\n";
				textArea.append(str);
				statusL.setText(resources.getString("status.error"));
	      
				stopThread();

				return;
			}

			File infile=new File(infileS);
			checksum1=checkSum(infile);
			if(!checkSuccess)
			{
	       		textArea.selectAll();
	       		textArea.replaceSelection("");
	       		String ss=resources.getString("error")+checkError+"\n";
				textArea.append(ss);
				statusL.setText(resources.getString("status.error"));
	       	
				stopThread();

				return;
			}

        	textArea.selectAll();
        	textArea.replaceSelection("");
        	str=resources.getString("md5sum.main.display")+"\n"+checksum1+"\n";
        	str=str+resources.getString("md5sum.display")+"\n"+checksum2+"\n"+"\n";

			if (checksum1.equals(checksum2))
			{
				str=str+resources.getString("md5sum.compare.success")+"\n";
			}
			else
			{
				str=str+resources.getString("md5sum.compare.fail")+"\n";
			}

			textArea.append(str);
			statusL.setText(resources.getString("status.compare"));

		}
		else if (anotherFileRb.isSelected())
		{
			File infile=new File(infileS);
			checksum1=checkSum(infile);
			if(!checkSuccess)
			{
	       		textArea.selectAll();
	       		textArea.replaceSelection("");
	       		String ss=resources.getString("error")+checkError+"\n";
				textArea.append(ss);
				statusL.setText(resources.getString("status.error"));
	       		
				stopThread();

				return;
			}

			File infile2=new File(infile2S);
			checksum2=checkSum(infile2);
			if(!checkSuccess)
			{
	       		textArea.selectAll();
	       		textArea.replaceSelection("");
	       		String ss=resources.getString("error")+checkError+"\n";
				textArea.append(ss);
				statusL.setText(resources.getString("status.error"));
	     
				stopThread();

				return;
			}

        	textArea.selectAll();
        	textArea.replaceSelection("");
        	str=resources.getString("md5sum.main.display")+"\n"+checksum1+"\n";
        	str=str+resources.getString("md5sum.second.display")+"\n"+checksum2+"\n"+"\n";

			if (checksum1.equals(checksum2))
			{
				str=str+resources.getString("md5sum.compare.success")+"\n";
			}
			else
			{
				str=str+resources.getString("md5sum.compare.fail")+"\n";
			}

			textArea.append(str);
			statusL.setText(resources.getString("status.compare"));

		}

		
		stopThread();

	}

    /**
     * The key method for generating the checksum
     * @param file
     * the file 
     * @return the file checksum
     */
	private String checkSum(File file)
	{
		ChecksumUtilities ch = new ChecksumUtilities();
		String str = ch.checkSum(file);//ChecksumUtilities.checkSum(file);
		checkSuccess = ch.getCheckSuccess();//ChecksumUtilities.checkSuccess;
		checkError = ch.getCheckError();//ChecksumUtilities.checkError;
		return str;		
	}

	/**
	 * Creating menu bar
	 * @param resources
	 * the resources
	 * @return the menu bar
	 */
    private JMenuBar createMenuBar(ResourceBundle resources)
    {

        // create the menus
        JMenuBar menuBar = new JMenuBar();

        String label;
        Character mnemonic;

        // first the file menu
        label = resources.getString("menu.file");
        mnemonic = (Character) resources.getObject("menu.file.mnemonic");
        JMenu fileMenu = new JMenu(label, true);
        fileMenu.setMnemonic(mnemonic.charValue());

        label = resources.getString("menu.file.open");
        mnemonic = (Character) resources.getObject("menu.file.open.mnemonic");
        JMenuItem openItem = new JMenuItem(label, mnemonic.charValue());
        openItem.setActionCommand(OPEN_COMMAND);
        openItem.addActionListener(this);
        fileMenu.add(openItem);

        label = resources.getString("menu.file.open2");
        mnemonic = (Character) resources.getObject("menu.file.open2.mnemonic");
        JMenuItem open2Item = new JMenuItem(label, mnemonic.charValue());
        open2Item.setActionCommand(OPEN2_COMMAND);
        open2Item.addActionListener(this);
        fileMenu.add(open2Item);

        label = resources.getString("menu.file.generate");
        mnemonic = (Character) resources.getObject("menu.file.generate.mnemonic");
        JMenuItem generateItem = new JMenuItem(label, mnemonic.charValue());
        generateItem.setActionCommand(GENERATE_COMMAND);
        generateItem.addActionListener(this);
        fileMenu.add(generateItem);

        label = resources.getString("menu.file.compare");
        mnemonic = (Character) resources.getObject("menu.file.compare.mnemonic");
        JMenuItem compareItem = new JMenuItem(label, mnemonic.charValue());
        compareItem.setActionCommand(COMPARE_COMMAND);
        compareItem.addActionListener(this);
        fileMenu.add(compareItem);

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

        // finally, glue together the menu and return it
        menuBar.add(fileMenu);
        menuBar.add(helpMenu);

        return menuBar;
    }

    /**
	 * Start the computation thread.
	 */
	private void startThread() {
		//stopAnim=false;
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
	 * Stop the computation thread.
	 */
	private void stopThread() {
		statusTh = null;
		frameNumber = 0;
		//stopAnim=true;//also redundant
		
		if (computationTh == null) {
			//stopAppend = false;// press kill button but simulation never
								// started!
			return;
		}
		
		computationTh = null;
		
		//if (stopAppend) {// kill button was pressed!
			//not used here..stopAppend redundant
		//	textArea.append(resources.getString("text.simulation.stop") + "\n");			
		//	stopAppend = false;
		//	String label = resources.getString("status.done");
		//	statusL.setText(label);
		//}
			
	}
	
	
	/**
	 * Run method ..current thread is initiated!
	 */
	public void run() {
		Thread.currentThread().setPriority(Thread.NORM_PRIORITY);// both thread
																	// same
																	// priority

		long startTime = System.currentTimeMillis();
		Thread currentThread = Thread.currentThread();
		//while (!stopAnim && currentThread == statusTh) {// if thread is status display
											// Thread!!
		while (currentThread == statusTh) {//works, stopAnim=redundant
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

		if (currentThread == computationTh) {// if thread is the main
												// computation Thread!!
			if (command.equals(GENERATE_COMMAND)) {
				generate();
			} else if (command.equals(COMPARE_COMMAND)) {
				compare();
			} 
		}
	}
	
	public static void main(String[] args)
    {
		LookAndFeel.loadLookAndFeel();
	    
	    new FileCheckSumGenerator();
    }
}