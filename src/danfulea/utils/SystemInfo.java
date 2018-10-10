package danfulea.utils;

import java.awt.DisplayMode;
import java.awt.GraphicsDevice;
import java.awt.GraphicsEnvironment;
import java.net.InetAddress;
import java.net.NetworkInterface;
import java.util.Enumeration;
import java.util.Properties;

/**
 * Utility class for displaying some system information.
 * 
 * @author Dan Fulea, 11 AUG. 2016
 *
 */
public class SystemInfo {

	private String machineName = "";
	private String machineIPAddress = "";
	private String machineMacAddress = "";

	/**
	 * The interface for displaying the results
	 */
	public static MessageRetriever mr;

	/**
	 * Available only if
	 * {@link danfulea.utils.SystemInfo#tryFetchingLocalNetwokInfo()} succeeds.
	 * 
	 * @return the basic machine network MAC address
	 */
	public String getMachineMacAddress() {// base network card
		return machineMacAddress;
	}

	/**
	 * Available only if
	 * {@link danfulea.utils.SystemInfo#tryFetchingLocalNetwokInfo()} succeeds.
	 * 
	 * @return the machine network name
	 */
	public String getMachineName() {
		return machineName;
	}

	/**
	 * Available only if
	 * {@link danfulea.utils.SystemInfo#tryFetchingLocalNetwokInfo()} succeeds.
	 * 
	 * @return the machine network local IP address
	 */
	public String getMachineIPAddress() {
		return machineIPAddress;
	}

	/**
	 * Print some system properties relevant to Java(tm).
	 */
	public static void printJavaProperties() {
		Properties p = System.getProperties();
		Enumeration<Object> keys = p.keys();
		while (keys.hasMoreElements()) {
			String key = (String) keys.nextElement();
			String value = (String) p.get(key);

			if (mr != null)
				mr.printSequence(key + ": " + value);

			// System.out.println(key + ": " + value);
		}
	}

	/**
	 * Print the current and supported display modes associated to video card.
	 */
	public static void printSupportedGraphicDisplays() {
		GraphicsEnvironment ge = GraphicsEnvironment.getLocalGraphicsEnvironment();
		GraphicsDevice de = ge.getDefaultScreenDevice();
		DisplayMode dm = de.getDisplayMode();
		int w = dm.getWidth();
		int h = dm.getHeight();
		int bd = dm.getBitDepth();
		int rate = dm.getRefreshRate();

		if (mr != null) {
			mr.printSequence("Current display mode:");
			mr.printSequence(w + " x " + h + "; color: " + bd + "; refresh rate [Hz]: " + rate);
		}
		// System.out.println("Current display mode:");
		// System.out.println(w+" x "+h+"; color: "+bd+"; refresh rate [Hz]:
		// "+rate);

		DisplayMode[] dmodes = de.getDisplayModes();
		int ww;
		int hh;
		int bdepth;
		int rrate;

		if (mr != null)
			mr.printSequence("Available display modes:");
		// System.out.println("Available display modes:");
		for (int i = 0; i < dmodes.length; i++) {
			ww = dmodes[i].getWidth();
			hh = dmodes[i].getHeight();
			bdepth = dmodes[i].getBitDepth();
			rrate = dmodes[i].getRefreshRate();
			if (mr != null)
				mr.printSequence(ww + " x " + hh + "; color: " + bdepth + "; refresh rate [Hz]: " + rrate);
			// System.out.println(ww+" x "+hh+"; color: "+bdepth+"; refresh rate
			// [Hz]: "+rrate);
		}
	}

	/**
	 * 
	 * @return true if Windows OS
	 */
	public static boolean isWindows() {

		String os = System.getProperty("os.name").toLowerCase();
		// windows
		return (os.indexOf("win") >= 0);

	}

	/**
	 * 
	 * @return true if Mac OS
	 */
	public static boolean isMac() {

		String os = System.getProperty("os.name").toLowerCase();
		// Mac
		return (os.indexOf("mac") >= 0);

	}

	/**
	 * 
	 * @return true if Linux OS
	 */
	public static boolean isLinux() {

		String os = System.getProperty("os.name").toLowerCase();
		// linux
		return (os.indexOf("nux") >= 0);

	}

	/**
	 * 
	 * @return true if Unix OS
	 */
	public static boolean isUnix() {

		String os = System.getProperty("os.name").toLowerCase();
		// unix
		return (os.indexOf("nix") >= 0);

	}

	/**
	 * Try to retrieve machine network name, relevant local IP address and its
	 * associated MAC address.
	 */
	public void tryFetchingLocalNetwokInfo() {
		byte[] macAddress;
		try {
			Enumeration<NetworkInterface> ifaces = NetworkInterface.getNetworkInterfaces();

			boolean found = false;
			while (ifaces.hasMoreElements()) {
				NetworkInterface iface = (NetworkInterface) ifaces.nextElement();
				// --------------
				Enumeration<NetworkInterface> virtualIfaces = iface.getSubInterfaces();

				while (virtualIfaces.hasMoreElements()) {
					NetworkInterface viface = (NetworkInterface) virtualIfaces.nextElement();

					if (mr != null)
						mr.printSequence(iface.getDisplayName() + " VIRT " + viface.getDisplayName());
					// System.out.println(iface.getDisplayName() + " VIRT " +
					// viface.getDisplayName());

					Enumeration<InetAddress> vaddrs = viface.getInetAddresses();
					while (vaddrs.hasMoreElements()) {
						InetAddress vaddr = (InetAddress) vaddrs.nextElement();

						if (mr != null)
							mr.printSequence("\t" + "VIRT IP= " + vaddr.toString());
						// System.out.println("\t" + "VIRT IP= " +
						// vaddr.toString());
					}
				}
				// ---------------------------
				// Any address in the range 127.xxx.xxx.xxx is a "loopback"
				// address.
				// /Any address in the range 192.168.xxx.xxx is a private (aka
				// site local) IP address
				// same for: 10.xxx.xxx.xxx; 72.16.xxx.xxx through
				// 172.31.xxx.xxx.
				if (mr != null)
					mr.printSequence("Real iface addresses: " + iface.getDisplayName());
				// System.out.println("Real iface addresses: " +
				// iface.getDisplayName());
				String ifacename = iface.getDisplayName();

				Enumeration<InetAddress> raddrs = iface.getInetAddresses();
				while (raddrs.hasMoreElements()) {
					InetAddress raddr = (InetAddress) raddrs.nextElement();

					if (mr != null)
						mr.printSequence("\t" + "IP= " + raddr.toString());
					// System.out.println("\t" + "IP= " + raddr.toString());
					String addS = raddr.toString();
					if ((addS.startsWith("/192.168.") || addS.startsWith("/10.") || addS.startsWith("/172.")) && // !ifacename.startsWith("Virtual")){
							ifacename.indexOf("Virt") < 0) {// not find virtual!
						machineIPAddress = addS.substring(1);
						machineName = raddr.getHostName();
						found = true;
					}
				}
				// ==========
				macAddress = iface.getHardwareAddress();
				StringBuilder mac = new StringBuilder();
				if (macAddress != null) {
					for (byte b : macAddress) {
						mac.append(String.format("%1$02X ", b));
					}
				}
				String macString = "" + mac;

				if (mr != null)
					mr.printSequence("MAC address (if any): " + macString + "\n");
				// System.out.println("MAC address (if any): " + macString +
				// "\n");
				if (found) {
					machineMacAddress = macString;
					found = false;
				}

				// ==========
			}
		} catch (Exception ex) {
			ex.printStackTrace();
		}
		// System.out.println("Found local IP Address: " + machineIPAddress);
		// System.out.println("MachineName:" + getMachineName() + " ;IP: " +
		// getMachineIPAddress() + "; MAC= " + getMachineMacAddress());
	}

	// ===========testing---------------------
	/* public static void main(String[] args) {
	 printJavaProperties();

	 System.out.println(isWindows()+" "+isUnix());
	 SystemInfo so = new SystemInfo();
	 so.tryFetchingLocalNetwokInfo();
	 System.out.println("MachineName:" + so.getMachineName() + " ;IP: " +
	 so.getMachineIPAddress() + "; MAC= " + so.getMachineMacAddress());

	 printSupportedGraphicDisplays();
	 }*/
	// ============END TESTING
}
