package danfulea.crypto.fcsgapp.resources;

import java.util.ListResourceBundle;

/**
 * resources for FileCheckSumGeneratorResources
 *
 *
 * @author Dan Fulea, 13 JUL. 2010
 *
 */

public class FileCheckSumGeneratorResources extends ListResourceBundle
{

    /**
     * Returns the array of strings in the resource bundle.
     *
     * @return the resources.
     */
    public Object[][] getContents()
    {
        return CONTENTS;
    }

    /** The resources to be localised. */
    private static final Object[][] CONTENTS =
    {

    	//{"form.icon.url", "/danfulea/crypto/fcsgapp/resources/Dan.ico"},//ico files not supported       
        {"form.icon.url", "/danfulea/crypto/fcsgapp/resources/personal.png"},//"/danfulea/crypto/fcsgapp/resources/personal.jpg"},
        {"logo.icon.url", "/danfulea/crypto/fcsgapp/resources/personal.png"},//"/danfulea/crypto/fcsgapp/resources/Kerrigan.jpg"},
       
        {"about.icon", "/danfulea/crypto/fcsgapp/resources/about.gif"},
		{"exit.icon", "/danfulea/crypto/fcsgapp/resources/exit.gif"},
		{"compare.icon", "/danfulea/crypto/fcsgapp/resources/libra.gif"},
		{"open.icon", "/danfulea/crypto/fcsgapp/resources/open.gif"},
		{"generate.icon", "/danfulea/crypto/fcsgapp/resources/save.gif"},
		
        {"plswait.title",      "Please wait!"},
        {"plswait.label",      "Work in progres..."},
        
        {"Application.NAME", "FileCheckSumGenerator: MD5SUM check and compare"},
        {"About.NAME", "About" },//AboutFrame requires this

        {"Author", "Author:"},
        {"Author.name", "Dan Fulea , fulea.dan@gmail.com"},        

        {"Version", "Version:"},
        {"Version.name", "FileCheckSumGenerator 0.5"},

        { "License.name", "BSD License"},
		{ "License",
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

		{"pleaseWait.label", "Work in progress!"},

        {"mainPanel.textArea.label", "Results:"},
        
        { "status.computing", "Computing..." },
        { "status.done", "Done! " },
        
        {"menu.file", "File"},
        {"menu.file.mnemonic", new Character('F') },

        {"menu.file.open", "Open the main file..."},
        {"menu.file.open.toolTip", "Open the main file in order to generate its MD5SUM file"},
        {"menu.file.open.mnemonic", new Character('O') },

        {"menu.file.open2", "Open the compare with file..."},
        {"menu.file.open2.toolTip", "Open the second file in order to make the comparison (direct main-second file comparison or main-MD5SUM file comparison)"},
        {"menu.file.open2.mnemonic", new Character('p') },

        {"menu.file.generate", "Generate MD5SUM file..."},
        {"menu.file.generate.mnemonic", new Character('G') },
        {"menu.file.generate.toolTip", "Generate the MD5SUM file related to main file"},

        {"menu.file.compare", "Compare"},
        {"menu.file.compare.mnemonic", new Character('C') },
        {"menu.file.compare.toolTip", "Compare the main file with the second file"},

        {"menu.file.exit", "Close"},
        {"menu.file.exit.mnemonic", new Character('C') },
        {"menu.file.exit.toolTip", "Close the application"},

        {"menu.help", "Help"},
        {"menu.help.mnemonic", new Character('H')},

        {"menu.help.about", "About..."},
        {"menu.help.about.mnemonic", new Character('A')},
        {"menu.help.about.toolTip", "Several informations about this application"},

        {"titleBorder.open", "Main file"},
        {"titleBorder.open2", "Compare with"},

        {"another.file", "another file"},
        {"md5.file", "MD5SUM file"},

		

		{"file.extension","md5"},
		{"file.description","md5sum file"},

		{"status.wait", "Waiting for your action!"},
		{"status.open", "Open: "},
		{"status.open2", "Compare with: "},
		{"status.open.fail", "Warning: Open file failure!"},
		{"status.save", "Save: "},
		{"status.compare", "Comparison done! "},
		{"status.error", "Unexpected error!"},

        {"dialog.exit.title", "Confirm..."},
        {"dialog.exit.message", "Are you sure?"},
        {"dialog.exit.buttons", new Object[]{"Yes","No"}},

        {"open.fail", "Please open the main file first!"},
        {"open2.fail", "Please open the compare with file!"},
        {"md5sum.main.display", "The main file checksum (MD5 hashes, hex): "},
        {"md5sum.second.display", "The second file checksum (MD5 hashes, hex): "},
        {"md5sum.display", "The pre-computed file checksum (MD5 hashes, hex): "},

        {"md5sum.compare.success", "Comparison status: Good file!"},
        {"md5sum.compare.fail", "Comparison status: Bad file!"},
        {"error", "Unexpected error: "},

        {"dialog.overwrite.title", "Overwriting..."},
        {"dialog.overwrite.message", "Are you sure?"},
        {"dialog.overwrite.buttons", new Object[]{"Yes","No"}},

    };

}
