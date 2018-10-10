package danfulea.network;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.io.IOException;
import java.net.ServerSocket;
import java.net.Socket;
import java.util.HashMap;

import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTextArea;
import javax.swing.SwingUtilities;

import danfulea.utils.FrameUtilities;

/**
 * A simple messenger server, a skype-like messenger
 * It is very simple and very efficient to use in an internal network environment, such as 
 * internal network of a company. The company users can use this messenger because is light, fast and 
 * the default setting is to accept connections only from inside internal network (e.g. server IP = "192.168.0.100").
 * Therefore, it is a better alternative to Skype, Yahoo messenger etc.
 * <p>
 * By default, it listens on port 6789 and accepts simultaneously 100 connections (change this if you need in 
 * {@link danfulea.network.ServerMessenger#startRunning()} method).
 * </p>
 * <p>
 * In a real environment, the server jar file must be place in a valid server service location, for example inside htdocs if using
 * XAMPP server or inside www folder if using WampServer etc.
 * </p>
 * <p>
 * To create a jar file: In eclipse, select ServerMessenger.java and ClientWorker.java (from danfulea.network package) and FrameUtilities 
 * (from danfulea.utils package) using CTRL key and left click. Then right click on this selection and choose export.
 * Select Java/jar file and then browse to a desired location where you want to create the jar, give a name (e.g. ServerMsgr.jar). 
 * Next, next and choose the manifest: /danfulea/META-INF-server/MANIFEST.MF. That's it.
 * </p>
 * 
 * @author Dan Fulea, 11 AUG. 2016
 *
 */
public class ServerMessenger extends JFrame implements ActionListener{

	private static final long serialVersionUID = 1L;
	private static final Dimension PREFERRED_SIZE = new Dimension(400, 200);

	public static Color bkgColor = new Color(230, 255, 210, 255);// Linux mint
	public static Color foreColor = Color.black;// Color.white;
	public static Color textAreaBkgColor = Color.white;// Color.black;
	public static Color textAreaForeColor = Color.black;// Color.yellow;
	
	private JTextArea chatWindow;

	private static final String EXIT_COMMAND = "EXIT";

	private ServerSocket server;
	private HashMap<Integer, Socket> clients;
	private int icount;
	private String clientName="";
	
	/**
	 * Constructor
	 */
	public ServerMessenger(){
		this.setTitle("Server-Messenger");
		clients = new HashMap<Integer, Socket>();//ArrayList<Socket>();
		icount=0;

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
	 * 
	 * @return the HashMap containing the clients (connections to this server)
	 */
	public HashMap<Integer, Socket> getClients(){
		return clients;
	}
	
	/**
	 * Delete a particular client from the map. The connection was lost from client messenger application (attemptExit method took care of that).
	 * @param index
	 * the index in the map
	 */
	public void deleteClient(int index){
		clients.remove(index);
	}
	
	/**
	 * 
	 * @return the client name
	 */
	public String getClientName(){
		return clientName;
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
	 * Creates main panel
	 * @return the main panel
	 */
	private JPanel createMainPanel() {
		JPanel panel = new JPanel(new BorderLayout());
		
		chatWindow = new JTextArea();
		chatWindow.setCaretPosition(0);
		chatWindow.setEditable(false);
		chatWindow.setText("");
		chatWindow.setLineWrap(true);
		chatWindow.setWrapStyleWord(true);

		panel.add(new JScrollPane(chatWindow), BorderLayout.CENTER);
		panel.setBackground(bkgColor);

		return panel;

	}

	/**
	 * Actions are set here
	 */
	public void actionPerformed(ActionEvent arg0) {
		// TODO Auto-generated method stub
		String command = arg0.getActionCommand();
		if (command.equals(EXIT_COMMAND)) {
			attemptExit();
		}
	}

	/**
	 * Application close!
	 */
	private void attemptExit() {
		
		dispose();
		System.exit(0);
	}

	/**
	 * show message on chat area
	 * @param text
	 * the text to be shown
	 */
	private void showMessage(final String text) {
		// update parts of the GUI not all, namely only update chatWindow:
		//this.
		SwingUtilities.invokeLater(
				// updates part of the GUI in a thread
				new Runnable() {
					public void run() {
						
						chatWindow.append(text);
						//to scroll down automatically even if it loose focus.
						int lll = chatWindow.getText().length();
						//for efficiency, clear window after 1 Mil rows======
						if (lll>1000000){//2147483647 = max int value
							lll=1;
							chatWindow.selectAll();
							chatWindow.replaceSelection("");
						}
						//===================================================
						lll=lll-1;
						chatWindow.setCaretPosition(lll);	
						//END SCROLL AUTO
					}
				});// invoke later

	}

	/**
	 * Server runs and listens for connections.
	 * By default it listen on port 6789 and accepts simultaneously 100 connections (change this if you need).
	 */
	public void startRunning() {
		try {
			server = new ServerSocket(6789, 100);
			while (true) {
				// run over and over again
				// server listen for connections every second!!
				
				ClientWorker w;
			    try{
			      //server.accept returns a client connection
			      Socket client = server.accept();
			      
			      clientName = client.getInetAddress().getHostName();//getting its name
			      showMessage("\nNow connected to " + clientName);
			      //add to hashmap======
			      clients.put(new Integer(icount),client);			      
			      //===================		      
			      w = new ClientWorker(client, chatWindow,icount, this);//put client in this thread
			      icount=icount+1;
			      /////////////////
			      Thread t = new Thread(w);
			      t.start();
			      
			      
			    } catch (IOException e) {
			      System.out.println("Accept failed: 4444");
			      System.exit(-1);
			    }
			}
		} catch (IOException e) {
			showMessage("\nCannot listen on default port! ");
			System.exit(-1);
		}
	}

	/**
	 * Java(tm) specific entry point, main class
	 * @param args default args
	 */
	public static void main(String[] args) {
		ServerMessenger server = new ServerMessenger();
		server.startRunning();
	}
}
