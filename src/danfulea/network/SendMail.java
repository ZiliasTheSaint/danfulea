package danfulea.network;

import java.util.Properties;

import javax.activation.DataHandler;
import javax.activation.DataSource;
import javax.activation.FileDataSource;
import javax.mail.BodyPart;
import javax.mail.Multipart;
import javax.mail.Session;
import javax.mail.Transport;
import javax.mail.internet.InternetAddress;
import javax.mail.internet.MimeBodyPart;
import javax.mail.internet.MimeMessage;
import javax.mail.internet.MimeMultipart;
import javax.mail.internet.MimeMessage.RecipientType;

import danfulea.utils.MessageRetriever;

/**
 * Sometimes it is useful to send mail directly from an application. For example,
 * if the application is designed to generate some work report, it is very useful
 * to send that report as attachment from within the application.
 * <p>
 * Sometimes, in order for the receiver to actually get the mail he must change some settings:
 * If he has an yahoo account, then he can go to https://login.yahoo.com/account/security#other-apps
 * and turn on "Allow apps that use less secure sign in". That's because external application
 * trying to send mail to yahoo is regarded as less secure. Same thing for gmail accounts:
 * go to https://www.google.com/settings/security/lesssecureapps and turn on the access for
 * less secure apps. Usually, if the receiver has a personal email-server these adjustments are not required.
 * 
 * @author Dan Fulea, 11 AUG. 2016
 *
 */
public class SendMail {

	private String from;
	private String password;
	private String to;
	private String subject;
	private String text;
	private String[] attachments;

	private String host;
	private String port;
	
	/**
	 * The interface for displaying the results
	 */
	public static MessageRetriever mr;
	
	/**
	 * Constructor of SendMail object
	 * @param from
	 * from e.g. "john.doe@gmail.com" or john.doe@johnDoeSite.org
	 * @param to
	 * to e.g. "john.doe@yahoo.com" or another.doe@anotherDoeSite.org
	 * @param subject
	 * the subject, e.g. Report
	 * @param text
	 * the actual text (message body)	 * 
	 * @param password
	 * the password associated to "from" recipient, e.g. johnDoePasswordForGmail
	 * @param attachments
	 * the attachment, e.g. {report1.xls, report1.doc};
	 */
	public SendMail(String from, String to, String subject, String text,
			String password, String[] attachments) {
		this.from = from;
		this.to = to;
		this.subject = subject;
		this.text = text;
		this.password = password;
		this.attachments = attachments;		
	}
	
	/**
	 * 
	 * @param host
	 * the smtp server used for sending email, e.g. personalServer.org
	 * @param port
	 * the port for email service, e.g. 25
	 */
	public void setHostAndPort(String host, String port) {
		this.host = host;
		this.port = port;
	}

	/**
	 * 
	 * @return the host
	 */
	public String getHost() {
		return this.host;
	}

	/**
	 * 
	 * @return the port
	 */
	public String getPort() {
		return this.port;
	}
	
	/**
	 * Try to send the email.	
	 * @param props
	 * Properties object holding the host (the smtp server used for sending email, e.g. personalServer.org) and the corresponding port used for sending email
	 * @throws Exception
	 * can throw this exception
	 */
	private void trySendingEMail(Properties props)//(String host, Properties props)
			throws Exception {
		//--------------
		String host = props.getProperty("mail.smtp.host");
		//----------------		
		Session session = Session.getDefaultInstance(props, null);
		MimeMessage message = new MimeMessage(session);

		InternetAddress fromAddress = new InternetAddress(from);
		InternetAddress toAddress = new InternetAddress(to);

		message.setFrom(fromAddress);
		message.setRecipient(RecipientType.TO, toAddress);
		message.setSubject(subject);
		
		// Create a message part to represent the body text
		BodyPart messageBodyPart = new MimeBodyPart();
		messageBodyPart.setText(text);

		// use a MimeMultipart as we need to handle the file attachments
		Multipart multipart = new MimeMultipart();

		// add the message body to the mime message
		multipart.addBodyPart(messageBodyPart);

		// add any file attachments to the message
		addAtachments(attachments, multipart);

		// Put all message parts in the message
		message.setContent(multipart);

		Transport transport = session.getTransport("smtp");
		transport.connect(host, from, password);
		transport.sendMessage(message, message.getAllRecipients());
		transport.close();
	}

	/**
	 * Add attachments
	 * 
	 * @param attachments
	 * the attachment String array
	 * @param multipart
	 * the multipart object which handles these attachments
	 * @throws Exception
	 * may throw this exception
	 */
	private void addAtachments(String[] attachments, Multipart multipart)
			throws Exception// MessagingException, AddressException
	{
		if (attachments == null)
			return;

		for (int i = 0; i <= attachments.length - 1; i++) {
			String filename = attachments[i];
			MimeBodyPart attachmentBodyPart = new MimeBodyPart();

			// use a JAF FileDataSource as it does MIME type detection
			DataSource source = new FileDataSource(filename);
			attachmentBodyPart.setDataHandler(new DataHandler(source));

			// assume that the filename you want to send is the same as the
			// actual file name - could alter this to remove the file path
			attachmentBodyPart.setFileName(filename);

			// add the attachment
			multipart.addBodyPart(attachmentBodyPart);
		}
	}

	/**
	 * Send mail as if it is sent by your gmail account
	 */
	public void sendFromGmail() {
		String host = "smtp.gmail.com";
		String port = "587";

		Properties props = System.getProperties();
		props.put("mail.smtp.host", host);
		//with TSL
		props.put("mail.smtp.port", port);
		props.put("mail.smtp.starttls.enable", "true");		

		props.put("mail.smtp.auth", "true");
		
		try {

			trySendingEMail(props);//trySendingEMail(host, props);

			if (mr != null)
				mr.printSequence("Email sent successfully!");
			//System.out.println("succes");
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

	}

	/**
	 * Send mail as if it is sent by your yahoo account
	 */
	public void sendFromYahoo() {
		String host = "smtp.mail.yahoo.com";
		String port = "465";

		Properties props = System.getProperties();
		props.put("mail.smtp.host", host);
		props.put("mail.smtp.port", port);

		// props.put("mail.smtp.starttls.enable", "true");
		props.put("mail.smtp.socketFactory.class",
				"javax.net.ssl.SSLSocketFactory");
		props.put("mail.smtp.ssl", "true");

		try {

			trySendingEMail(props);//trySendingEMail(host, props);

			if (mr != null)
				mr.printSequence("Email sent successfully!");
			//System.out.println("succes");
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

	}
	

	/**
	 * Send mail as if it is sent by your default (personal email server) account
	 */
	public void sendDefaultEMail() {// String host, String port) {

		String host = getHost();
		String port = getPort();

		Properties props = System.getProperties();
		props.put("mail.smtp.host", host);
		props.put("mail.smtp.port", port);

		try {

			trySendingEMail(props);//trySendingEMail(host, props);

			if (mr != null)
				mr.printSequence("Email sent successfully!");
			//System.out.println("succes");
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

	}
	//========TESTING
	/*public static void main(String[] args) {
		SendMail sendMail = null;
		String from = "john.doe@gmail.com";/////
		String password = "yourPassword";//////
		String to = "john.doe@yahoo.com";/////
		String subject = "Java(tm) mail final!";
		String message = "A java(tm) mail test message";
		String[] attachments =  null;//no attachments
		//String[] attachments = { "test.txt" };//works

		//sendMail = new SendMail(from, to, subject, message, password,
		//		attachments);
		//sendMail.sendFromGmail();//works!

		from = "john.doe@yahoo.com";////
		password ="yourPassword";/////
		to = "john.doe@gmail.com";/////

		subject = "Testing java(tm) mail !!!";
		message = "A java(tm) mail test message";

		//sendMail = new SendMail(from, to, subject, message, password,
		//		attachments);
		//sendMail.sendFromYahoo();//works

		from = "john.doe@personalServer.org";/////
		password = "yourPassword";///////////
		to = "john.doe@gmail.com";////////////////

		subject = "Test java mail from my server!!!";
		message = "Java mail test message and attach";
		sendMail = new SendMail(from, to, subject, message, password,
				attachments);		
		//sendMail.setHostAndPort("personalServer.org", "25");
		//sendMail.sendDefaultEMail();//works

	}*/
	//end TESTING
}
