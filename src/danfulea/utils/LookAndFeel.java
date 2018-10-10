package danfulea.utils;

import java.io.File;
import java.io.FileReader;

import javax.swing.JFrame;
import javax.swing.UIManager;
import javax.swing.UIManager.LookAndFeelInfo;

import danfulea.math.Convertor;

/**
 * 
 * Utility class for loading look and feel (LF) in order to have nice-looking
 * GUI. Usage: just put "LookAndFeel.loadLookAndFeel();" into application main
 * class before instantiating the GUI.
 * 
 * 
 * @author Dan Fulea, 14 APR. 2011
 *
 */

public class LookAndFeel {
	/**
	 * The file located in application folder which stores the name of look and
	 * feel to be use!
	 */
	private static String filename = "JavaLookAndFeelLoader.laf";

	/**
	 * Setting a default system look and feel!
	 */
	public static void defaultLookAndFeel() {

		JFrame.setDefaultLookAndFeelDecorated(true);// the key to set Java Look
													// and Feel from the start!
		try {
			UIManager.setLookAndFeel(UIManager.getCrossPlatformLookAndFeelClassName());
		} catch (final Exception exc1) {
			System.err.println("Error defaultLookAndFell " + exc1);
		}
		// }
	}

	/**
	 * Try loading the look and feel!
	 */
	public static void loadLookAndFeel() {
		String fileSeparator = System.getProperty("file.separator");
		String curentDir = System.getProperty("user.dir");
		String filename1 = curentDir + fileSeparator + filename;

		File f = new File(filename1);
		int i = 0;
		String desiredLF = "";
		boolean foundB = false;
		if (f.exists()) {
			try {
				FileReader fr = new FileReader(f);
				while ((i = fr.read()) != -1) {
					String s1 = new String();
					s1 = Convertor.asciiToStr(i);
					desiredLF = desiredLF + s1;
				}
				fr.close();

				JFrame.setDefaultLookAndFeelDecorated(true);

				for (LookAndFeelInfo info : UIManager.getInstalledLookAndFeels()) {
					if (desiredLF.equals(info.getName())) {
						UIManager.setLookAndFeel(info.getClassName());
						foundB = true;
						break;
					}
				}

				if (!foundB) {
					if (desiredLF.equals("System")) {
						foundB = true;
						UIManager.setLookAndFeel(UIManager.getSystemLookAndFeelClassName());
					} else if (desiredLF.equals("Java")) {
						foundB = true;
						UIManager.setLookAndFeel(UIManager.getCrossPlatformLookAndFeelClassName());
					}
				}

				if (!foundB)
					UIManager.setLookAndFeel(desiredLF);

				// empty
				if (desiredLF.equals("")) {
					defaultLookAndFeel();
				}

			} catch (Exception e) {
				defaultLookAndFeel();
			}
		} else {
			defaultLookAndFeel();
		}
	}
}
