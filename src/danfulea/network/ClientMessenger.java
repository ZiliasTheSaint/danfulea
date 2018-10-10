package danfulea.network;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.io.BufferedReader;
import java.io.EOFException;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.net.InetAddress;
import java.net.Socket;
import java.net.UnknownHostException;

import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTextArea;
import javax.swing.JTextField;
import javax.swing.SwingUtilities;

import danfulea.math.Convertor;
import danfulea.utils.FrameUtilities;

/**
 * The messenger client associated with messenger server.
 * It is very simple and very efficient to use in an internal network environment, such as 
 * internal network of a company. The company users can use this messenger because is light, fast and 
 * the default setting is to accept connections only from inside internal network (e.g. server IP = "192.168.0.100").
 * Therefore, it is a better alternative to Skype, Yahoo messenger etc.
 * 
 * <p>
 * Don't forget to change the default serverIp text file to a real internal server IP!  
 * </p>
 * 
 * <p>
 * To create a jar file: In eclipse, select ClientMessenger.java (from danfulea.network package) and FrameUtilities 
 * (from danfulea.utils package) as well as Convertor.java (from danfulea.math package) using CTRL key and left click. Then right click on this selection and choose export.
 * Select Java/jar file and then browse to a desired location where you want to create the jar, give a name (e.g. ClientMsgr.jar). 
 * Next, next and choose the manifest: /danfulea/META-INF-client/MANIFEST.MF. That's it.
 * </p>
 * 
 * @author Dan Fulea, 11 AUG. 2016
 *
 */

public class ClientMessenger extends JFrame implements ActionListener{
	
	private static final long serialVersionUID = 1L;
	private static final Dimension PREFERRED_SIZE = new Dimension(600, 400);

	public static Color bkgColor = new Color(230, 255, 210, 255);// Linux mint
	
	private static final String EXIT_COMMAND = "EXIT";
	// area where you type the message to be send
	private JTextField userText;
	// display the conversion:
	private JTextArea chatWindow;
	
	private String serverIP;

	private Socket connection;
	private PrintWriter output;
	private BufferedReader input;
	private String clientName="CLIENT";
	
	private static String filename = "serverIp.txt";//store the IP
	
	/**
	 * The constructor. Default server IP is "localhost" or "127.0.0.1" which is the same thing.
	 * This is default setting for testing. In a real environment, use the internal server IP address, for
	 * example "192.168.0.100" (this means the server jar file must be place in a valid server service location, for example inside htdoc if using
	 * XAMPP server or inside www folder if using WampServer etc.).
	 * ClientMessenger class read the serverIp from a file located in same folder, thus please modify serverIp.txt accordingly!
	 */
	public ClientMessenger(){
		this.setTitle("Client-Messenger");

		// =============CHANGE THIS TO REAL IP WHEN DEPLOY ServerMessenger ON A SERVER!!
		serverIP = "localhost";//"127.0.0.1";// "192.168.0.100";//"127.0.0.1";// "192.168.0.100";
		serverIP = retrieveServerIp();
		String hostname = clientName;//just in case
		try
		{
		    InetAddress addr;
		    addr = InetAddress.getLocalHost();
		    hostname = addr.getHostName();		    
		}
		catch (UnknownHostException ex)
		{
		    System.out.println("Hostname can not be resolved");
		}
		clientName = hostname;
		// =================
		// the key to force attemptExit() method on close!!
		setDefaultCloseOperation(JFrame.DO_NOTHING_ON_CLOSE);
		addWindowListener(new WindowAdapter() {
			public void windowClosing(WindowEvent e) {
				attemptExit();
			}
		});

		createGUI();
		setDefaultLookAndFeelDecorated(true);
		
		FrameUtilities.centerFrameOnScreen(this);
		setVisible(true);
	}
	
	/**
	 * Read the server IP from serverIp.txt file. If fail or file not exist then localhost is used instead.
	 * @return the string representation of server IP.
	 */
	private String retrieveServerIp(){
		String ip = "";
		//----------------
		String fileSeparator = System.getProperty("file.separator");
		String curentDir = System.getProperty("user.dir");
		String filename1 = curentDir + fileSeparator + filename;//serverIp.txt
		
		File f = new File(filename1);
		int i = 0;
		
		//boolean foundB = false;
		if (f.exists()) {
			try {
				FileReader fr = new FileReader(f);
				while ((i = fr.read()) != -1) {
					String s1 = new String();
					s1 = Convertor.asciiToStr(i);
					ip = ip + s1;
				}
				fr.close();
			}catch (Exception e) {
				ip = "localhost";
			}
		} else {
			ip = "localhost";
		}
		//--------------
		//System.out.println(ip);
		return ip;
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
		
		setContentPane(new JScrollPane(content));
		content.setOpaque(true); // content panes must be opaque
		pack();
	}
	
	/**
	 * Main panel creation
	 * @return the main panel
	 */
	private JPanel createMainPanel() {
		JPanel panel = new JPanel(new BorderLayout());

		userText = new JTextField();
		userText.addActionListener(this);

		chatWindow = new JTextArea();
		chatWindow.setCaretPosition(0);		
		chatWindow.setEditable(false);
		chatWindow.setText("");
		chatWindow.setLineWrap(true);
		chatWindow.setWrapStyleWord(true);
		
		panel.add(userText, BorderLayout.NORTH);
		panel.add(new JScrollPane(chatWindow), BorderLayout.CENTER);
		panel.setBackground(bkgColor);

		return panel;

	}

	/**
	 * Actions are set up here
	 */
	public void actionPerformed(ActionEvent arg0) {
		
		String command = arg0.getActionCommand();
		if (command.equals(EXIT_COMMAND)) {
			attemptExit();
		} else if (arg0.getSource() == userText)// enterul
		{
			// do smthing
			sendMessage(userText.getText());
			userText.setText("");
			userText.requestFocusInWindow();
		}
		
		//Receive text from server AUTOMATICALLY, JUST BY PASSING FOLLOWING CODE HERE!!
		//not good though=>whileChatting	
	}
	
	/**
	 * Read and show message
	 * @throws IOException
	 * can throw this exception
	 */
	private void whileChatting() throws IOException{			
		
		do{
			//have a conversation
			String line = input.readLine();
			showMessage("\n" + line);			
		}while(true);
	}

	/**
	 * Default application close!
	 */
	private void attemptExit() {//closing window, not typing END

		try { // avoid memory leak-used when there is not need to
			if (output != null){
				sendMessage("END");//make sure to close connections on serveside
				output.close();// close stream to them
			}
			if (input != null)
				input.close();// close stream from them
			if (connection != null)
				connection.close();// close socket

		} catch (IOException ioe) {
			ioe.printStackTrace();
		}

		dispose();
		System.exit(0);
	}
	
	/**
	 * Exit method when client explicitly send END command.
	 */
	private void attemptExit2() {

		try { // avoid memory leak-used when there is not need to
			if (output != null){				
				output.close();// close stream to them
			}
			if (input != null)
				input.close();// close stream from them
			if (connection != null)
				connection.close();// close socket

		} catch (IOException ioe) {
			ioe.printStackTrace();
		}

		dispose();
		System.exit(0);
	}

	/**
	 * Main run method
	 */			
	public void startRunning() {		
		try {			
			connectToServer();
			setupStreams();// to actually have a conversation			
			whileChatting();// allow to send data -conversation-back and forth
		} catch (EOFException ex) {// Not really an error...just conversation ends!
			showMessage("\nClient ended the connection! ");
		} catch (IOException e) {// (Exception ex){			
			showMessage("\nCommunication failed! RESTART!");
		} finally {
			// nothing here..........closeCrap();
		}		
	}

	/**
	 * Sending message
	 * @param text
	 * the text to be send
	 */
	private void sendMessage(String text) {		
		try{			
			output.println(clientName+" - "+text);//"CLIENT - " +text);	
			if (text.equals("END")){
				attemptExit2();
			}
		}catch(Exception ex){
			showMessage("\nMessage send failed! RESTART!");
		}
	}
	
	/**
	 * Connect to the server
	 * @throws IOException
	 * can throw this exception
	 */
	private void connectToServer() throws IOException {
		showMessage("Attempting to connect...\n");
		// to a server and to a specific port given by the server admin!
		connection = new Socket(InetAddress.getByName(serverIP), 6789);

		showMessage("Connected to " + connection.getInetAddress().getHostName());
		showMessage("\nType END and press ENTER to terminate conversation (or close the window)!");
		// getInetAddress = returns your address from which you are now
		// connected to the server.
		// getHostNAme is actually your IP address converted to a string.
		// showMessage=accept strings
	}

	/**
	 * Set up communication streams
	 * @throws IOException
	 * can throw this exception
	 */
	private void setupStreams() throws IOException {
		output = new PrintWriter(connection.getOutputStream(), true);// auto-flush
		//flushes the stream..send all left-over bytes from buffer. 
		//Think of flushes the toilet :))))
		//push the rest of the crap (bytes of information) through.

		input = new BufferedReader(new InputStreamReader(connection.getInputStream()));		
		showMessage("\nStreams are now setup! \n");
	}

	/**
	 * Show message	coming from the server.		
	 * @param text
	 * the text 
	 */
	private void showMessage(final String text) {
		this.toFront();//highlight if minimized to discreetly alert client a message has been send to chat! 
		// update parts of the GUI not all, namely only update chatWindow:		
		SwingUtilities.invokeLater(				
				new Runnable() {
					public void run() {
						//int lll = chatWindow.getText().length();//testing
						//chatWindow.append(lll+" "+text);//testing						
						chatWindow.append(text);
						//to scroll down automatically even if it loose focus.
						int lll = chatWindow.getText().length();
							//////////for efficiency
						if (lll>1000000){//2147483647 = max int value
							lll=1;
							chatWindow.selectAll();
							chatWindow.replaceSelection("");
						}
							///////////////
						lll=lll-1;
						chatWindow.setCaretPosition(lll);	
						//END SCROLL AUTO
					}
				});// invoke later

	}

	/**
	 * Java(tm) entry point, main class
	 * @param args
	 * default args
	 */
	public static void main(String[] args) {
		ClientMessenger client = new ClientMessenger();
		client.startRunning();
	}
}
