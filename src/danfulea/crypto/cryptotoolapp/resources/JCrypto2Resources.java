package danfulea.crypto.cryptotoolapp.resources;

import java.util.ListResourceBundle;

/**
 * The resource collector.
 * 
 * @author Dan Fulea, 22 APR. 2015.
 */
public class JCrypto2Resources extends ListResourceBundle{
	/**
     * Returns the array of strings in the resource bundle.
     *
     * @return the resources.
     */
    public Object[][] getContents() {
        return CONTENTS;
    }

    /** The resources to be localised. */
    private static final Object[][] CONTENTS = {

         // menu labels...
        {"form.icon.url", "/danfulea/crypto/cryptotoolapp/resources/personal.png"},//"/danfulea/crypto/cryptotoolapp/resources/duke.png"},
        {"logo.icon.url", "/danfulea/crypto/cryptotoolapp/resources/personal.png"},//"/danfulea/crypto/cryptotoolapp/resources/globe.gif"},
        
        {"Application.NAME", "JCryptoTool -Sun Java(tm) implementation"},
        {"About.NAME", "About" },//AboutFrame requires this
        
        { "HowTo.NAME", "How to..." },
        { "SignFrame.NAME", "Sign data file to be send" },
        { "VerifySignedFile.NAME", "Verify signed files" },

        { "status.computing", "Computing..." },
        { "status.done", "Done! " },
        {"ds.PleaseWait.label", "Work in progress!"},

        {"Author", "Author:"},
        {"Author.name", "Dan Fulea,  fulea.dan@gmail.com"},
        {"Algorithms", "Algorithms:"},
        //{"Algorithms.name", "AES-256 bits for main encryption, RSA for AES-key encryption using public-private keys"},
        {"Algorithms.name", "AES-256 bits / RSA 2048bits for AES key"},
        {"Version", "Version:"},
        {"Version.name", "JCryptoTool 0.5"},
        
        { "License.name", "BSD License"},
        {"License", //"This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License (version 2) as published by the Free Software Foundation. \n\nThis program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details. \n"},
        "Copyright (c) 2014, Dan Fulea \nAll rights reserved.\n\nRedistribution and "
        + "use in source and binary forms, with or without modification, are permitted"
        + " provided that the following conditions are met:\n\n1. Redistributions of"
        + " source code must retain the above copyright notice, this list of conditions "
        + "and the following disclaimer.\n\n2. Redistributions in binary form must "
        + "reproduce the above copyright notice, this list of conditions and the "
        + "following disclaimer in the documentation and/or other materials provided with"
        + " the distribution.\n\n3. Neither the name of the copyright holder nor"
        + " the names of its contributors may be used to endorse or promote products "
        + "derived from this software without specific prior written permission.\n\nTHIS"
        + " SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS 'AS IS' AND"
        + " ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED "
        + "WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE "
        + "DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE "
        + "LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY,"
        + " OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF"
        + " SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS "
        + "INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN "
        + "CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)"
        + " ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF "
        + "THE POSSIBILITY OF SUCH DAMAGE.\n" },
        

        {"menu.file", "File"},
        {"menu.file.mnemonic", new Character('F') },
        {"menu.file.exit", "Close"},
        {"menu.file.exit.mnemonic", new Character('C') },
        {"menu.help", "Help"},
        {"menu.help.mnemonic", new Character('H')},
        {"menu.help.about", "About..."},
        {"menu.help.about.mnemonic", new Character('A')},
        
        { "menu.help.howTo", "How to..." },
		{ "menu.help.howTo.mnemonic", new Character('H') },
		{ "menu.help.howTo.toolTip",
				"Several hints and tips about this application" },

        { "menu.file.sign", "Sign files..." },
		{ "menu.file.sign.mnemonic", new Character('S') },
		{ "menu.file.sign.toolTip",
						"Sign files in order to be checked by receiver using our private key" },

		{ "menu.file.verify", "Verify signed files..." },
		{ "menu.file.verify.mnemonic", new Character('V') },
		{ "menu.file.verify.toolTip",
						"Verify the received signed files using the provider public key" },
//===========================
        {"dialog.exit.title", "Confirm..."},
        {"dialog.exit.message", "Are you sure?"},
        {"dialog.exit.buttons", new Object[]{"Yes","No"}},
        
        { "HowTo.title", "Hints and tips" },
		{
				"HowTo",
				"It is generally not advisable to use a public key encryption algorithm such as RSA to directly encrypt files, since\n"
				+"(a) public key encryption is slow, and (b) it will only let you encrypt small things.\n"
				+"The alternative, and commonly used approach, is to use a shared key algorithm to encrypt/decrypt the files," 
 +"and then use a public key algorithm to encrypt/decrypt the (randomly generated) key used by the shared key algorithm." 
 +"This has the benefit of fast file encryption/decryption whilst still requiring a non-shared private key to get" 
 +"access to the key needed to decrypt the files.\n\n"
				+"Common scenario:\n"
						+ "1. Generate a key pair=>2 files: public.key file and private.key file in the application current folder.\n"
						+ "2. By default, during encryption stage, it is auto-generated an AES encryption key and it is saved (and encrypted using the above public.key) into a file located in the application current folder\n"
						+ "3. During message/file encryption, the AES encryption key is used to encrypt a document and save its encrypted version into the aplication current folder. \n"
						+ "4. Send the public key file (optional for 5.1 below, send also: aes encryption key file and encrypted document) to JohnDoe. At this stage, JohnDoe can use public key to generate a new AES encryption key and "
 +"encrypt his document and send his encrypted AES key file and his encrypted document to us. We can decrypt the John's AES key using our private key associated to the public key and " 
 +"decrypt John's encrypted file using John's decrypted AES key.\n"
						+ "5.1 Copy our private key file to JohnDoe computer in private, in order to give him the ability to decrypt our AES key and use it to decrypt our file.\n"
						+ "OR BETTER:\n"
						+ "5.2 John Doe creates a new public-private keys and give us his own public key and we use it to encrypt message to John so he can decrypt and read them. \n"
						+ "6. Now we can exchange messages/files.\n"
						+"\n"
						+"7. Additional security option: In order to check if send/received files are un-altered, one can use File/sign and File/verify option and send/receive signature file along with encrypted data. \n"
						+"This will allow to sign/verify data using the same private-public key used in encryption. A medium computer can perform this operation on 600 Mb data in about 20 seconds. \n"
						+"Still, if time is the issue, one can download any checksum software, such as FCSG, generate checksum file corresponding to the data file and then sign/verify the checksum file instead of actual data file.\n"
						+"\n"
						+"Best case scenario: Use a checksum sofware such as FCSG and generate cheksum for your public key. Publish the cheksum result and upload the public key on an internet trusted server. JohnDoe can download the public key and he can verify (using for instance FCSG) its checksum. Then proceed with 1-6 or 1-7 steps."
 },


    };

}
