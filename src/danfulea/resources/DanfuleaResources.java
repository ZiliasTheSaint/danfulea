package danfulea.resources;

import java.awt.Color;
import java.util.ListResourceBundle;

/**
 * Stores resources for danfulea package!
 * 
 * @author Dan Fulea, 15.apr.2011
 * 
 */
public class DanfuleaResources extends ListResourceBundle {

	/**
	 * Default background color
	 */
	public static Color bkgColor = new Color(230, 255, 210, 255);//new Color(180, 220, 150, 255);//new Color(230, 255, 210, 255);// Linux mint
																	// green
																	// alike
	/**
	 * Default foreground color
	 */
	public static Color foreColor = Color.black;// Color.white;
	
	/**
	 * Default text area background color
	 */
	public static Color textAreaBkgColor = Color.white;// Color.black;
	
	/**
	 * Default text area foreground color
	 */
	public static Color textAreaForeColor = Color.black;// Color.yellow;
	
	/**
	 * Override method to return resources for danfulea package.
	 * 
	 * @return the resources.
	 */
	@Override
	protected Object[][] getContents() {
		// TODO Auto-generated method stub
		return CONTENTS;
	}

	/** The resources to be localized. */
	private static final Object[][] CONTENTS = {

			// Images=====================================================
			{ "form.icon.url", "/danfulea/resources/duke.png" },//AboutFrame requires this
			{ "logo.icon.url", "/danfulea/resources/globe.gif" },//AboutFrame requires this

			{ "img.zoom.in", "/danfulea/resources/zoom_in.png" },
			{ "img.zoom.out", "/danfulea/resources/zoom_out.png" },
			{ "img.pan.left", "/danfulea/resources/arrow_left.png" },
			{ "img.pan.up", "/danfulea/resources/arrow_up.png" },
			{ "img.pan.down", "/danfulea/resources/arrow_down.png" },
			{ "img.pan.right", "/danfulea/resources/arrow_right.png" },
			{ "img.pan.refresh", "/danfulea/resources/arrow_refresh.png" },

			{ "img.accept", "/danfulea/resources/accept.png" }, 
			{ "img.insert", "/danfulea/resources/add.png" },
			{ "img.delete", "/danfulea/resources/delete.png" },
			{ "img.delete.all", "/danfulea/resources/bin_empty.png" },
			{ "img.view", "/danfulea/resources/eye.png" },
			{ "img.set", "/danfulea/resources/cog.png" }, 
			{ "img.report", "/danfulea/resources/document_prepare.png" },
			{ "img.today", "/danfulea/resources/time.png" },
			{ "img.open.file", "/danfulea/resources/open_folder.png" },
			{ "img.open.database", "/danfulea/resources/database_connect.png" },
			{ "img.save.database", "/danfulea/resources/database_save.png" },
			{ "img.substract.bkg", "/danfulea/resources/database_go.png" },
			{ "img.close", "/danfulea/resources/cross.png" },
			{ "img.about", "/danfulea/resources/information.png" },
			{ "img.printer", "/danfulea/resources/printer.png" },
			// ================================================================

			// Database Test Class
			// resources======================================================
			{ "Application.NAME", "Fast DB - Test Class" }, //AboutFrame requires this
			{ "About.NAME", "About" },//AboutFrame requires this

			{ "Author", "Author:" },//AboutFrame requires this
			{ "Author.name", "Dan Fulea , fulea.dan@gmail.com" },//AboutFrame requires this

			{ "Version", "Version:" }, //AboutFrame requires this
			{ "Version.name", "Database TestClass 1.0" },//AboutFrame requires this

			{ "License.name", "BSD License"},
			{ "License",//AboutFrame requires this
					"Copyright (c) 2016, Dan Fulea, fulea.dan@gmail.com \nAll rights reserved.\n\nRedistribution and use in source and binary forms, "
					+ "with or without modification, are permitted provided that the following conditions are met:\n\n"
					+ "1. Redistributions of source code must retain the above copyright notice, this list of conditions and the "
					+ "following disclaimer.\n\n"
					+ "2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and "
					+ "the following disclaimer in the documentation and/or other materials provided with the distribution.\n\n"
					+ "3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or "
					+ "promote products derived from this software without specific prior written permission.\n\n"
					+ "THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS 'AS IS' AND ANY EXPRESS OR IMPLIED WARRANTIES,"
					+ " INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE"
					+ " ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, "
					+ "INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE "
					+ "GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY "
					+ "OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY "
					+ "WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.\n" },

			{ "menu.file", "File" }, 
			{ "menu.file.mnemonic", new Character('F') },
			
			{ "menu.file.exit", "Close" },			
			{ "menu.file.exit.mnemonic", new Character('C') },
			{ "menu.file.exit.toolTip", "Close the application" },
			
			{ "menu.help", "Help" },
			{ "menu.help.mnemonic", new Character('H') },

			{ "menu.help.about", "About..." },
			{ "menu.help.about.mnemonic", new Character('A') },
			{ "menu.help.about.toolTip",
					"Several informations about this application" },

			{ "menu.help.LF", "Look and feel..." },
			{ "menu.help.LF.mnemonic", new Character('L') },
			{ "menu.help.LF.toolTip", "Change application look and feel" },
			
			{ "data.load", "Data" },
			
			{ "add.button", "Add" },
			{ "add.button.mnemonic", new Character('A') },
			{ "add.button.toolTip", "Add data!" },

			// =====================================================================================			
			{ "dialog.exit.title", "Confirm..." },
			{ "dialog.exit.message", "Are you sure?" },
			{ "dialog.exit.buttons", new Object[] { "Yes", "No" } },
			{"dialog.exit.yes", "Yes"},
			{"dialog.exit.no", "No"},
			//====================================================================================
			
			{ "project.name", "danfulea" }, { "project.version", "3.0" }, 
			{ "project.info", "fulea.dan@gmail.com" },
			{ "project.copyright", "(C)opyright 2002-2016" + " Contributors" },
			{ "copy.time", "copy Time [ms] for a file of size [B] : " },
			{ "copy.sourceExists", "copyFile: Source file does not exist or is not a file: " },
			{ "copy.sourceRead", "copyFile: Source file cannot be read: " },
			{ "copy.destParentExists", "copyFile: Destination parent folder does not exist: " },
			{ "copy.destParentWritable", "copyFile: Destination parent folder is not writable: " },
			{ "copy.destFile", "copyFile: Destination is not a valid file: " },
			{ "copy.sourceDirExists", "copyFile: Source folder does not exist: " },
			{ "copy.sourceDirRead", "copyFile: Source folder cannot be read: " },
			{ "copy.destExists", "copyFile: Destination folder does not exist: " },
			{ "copy.copyed.100", "Copyed: 100 % " }, 
			{ "copy.copyed", "Copied: " },
			{ "copy.proc", " % " },

	};

}
