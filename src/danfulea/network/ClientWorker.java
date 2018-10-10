package danfulea.network;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.net.Socket;
import java.util.HashMap;

import javax.swing.JTextArea;

/**
 * The thread for messenger client which is connected to messenger server.
 * Each clients (connections to server) run in separate threads.
 * 
 * @author Dan Fulea, 11 AUG. 2016
 *
 */
public class ClientWorker implements Runnable {
	//The Java Virtual Machine continues to execute threads until either of the 
	//following occurs:
	//The exit method of class Runtime has been called and the security manager has 
	//permitted the exit operation to take place.
	//All threads that are not daemon threads have died, either by returning from 
	//the call to the run method or by throwing an exception that propagates beyond 
	//the run method.

	//SO BY REACHING THE END OF RUN method=>auto -stop the thread by JVM!!
	
	private Socket client;
	private JTextArea textArea;
	private int ihash;
	private ServerMessenger server;

	/**
	 * Constructor
	 * @param client
	 * the client connection, i.e. client socket
	 * @param textArea
	 * the GUI component where conversation is displayed
	 * @param ihash
	 * the client ID, unique for each client
	 * @param server
	 * the messenger server each client is connected to.
	 */
	ClientWorker(Socket client, JTextArea textArea, int ihash, ServerMessenger server) {
		this.client = client;
		this.textArea = textArea;
		this.ihash=ihash;
		this.server=server;
	}

	/**
	 * Thread run method
	 */
	public void run() {
		String line="";
		String endS="";
		BufferedReader in = null;
		
		try {
			in = new BufferedReader(new InputStreamReader(client.getInputStream()));			
		} catch (IOException e) {
			System.out.println("in or out failed");
			System.exit(-1);
		}

		do{
			try {
				line = in.readLine();
				////////////////////////////////////////////
				int len = line.length();
				endS = line.substring(len-3, len);//looking for END because when client send END message the connection is terminated
				////////////				
				if(endS.equals("END")){
					server.deleteClient(ihash);
				}
				//////////////
				
				// Send data back to client....all clients, the key of being a skype-like messenger!
				HashMap<Integer, Socket> clnts = server.getClients();				
				for (Integer itg : clnts.keySet()){
					int ii = itg.intValue();					
					Socket socket = clnts.get(ii);
					PrintWriter outp = new PrintWriter(socket.getOutputStream(), true);
					if(!endS.equals("END"))//not want to send END command to all!
						outp.println(line);
				}
				
				// Append data to text area
				textArea.append("\n"+line);
				//some housekeeping------------------------------------
				int lll = textArea.getText().length();
				//for efficiency, clear text when exceeds 1 mil. rows.
				if (lll>1000000){//2147483647 = max int value
					lll=1;
					textArea.selectAll();
					textArea.replaceSelection("");
				}
				
				lll=lll-1;
				textArea.setCaretPosition(lll);	
				//END HOUSEKEEPING--------------------------------------
			} catch (IOException e) {
				System.out.println("Read failed");
				System.exit(-1);
			}
		} while (!endS.equals("END"));
	}
}
