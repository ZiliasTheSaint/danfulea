package danfulea.network;

import java.util.ResourceBundle;

import danfulea.utils.AlertBox;
import javafx.application.Application;
import javafx.beans.value.ChangeListener;
import javafx.beans.value.ObservableValue;
import javafx.event.ActionEvent;
import javafx.event.EventHandler;
import javafx.scene.Scene;
import javafx.scene.control.TextField;
import javafx.scene.image.Image;
import javafx.scene.layout.VBox;
import javafx.scene.web.WebEngine;
import javafx.scene.web.WebView;
import javafx.stage.Stage;
import javafx.stage.WindowEvent;

/**
 * This class is a demo class for building a very simple web browser.
 * Since a modern web browser requires support for html5 and java script to be actually useful, then swing technology is not enough.
 * As an alternative, one can use new javafx technology. Let's try it!
 * 
 * @author Dan Fulea, 12 AUG. 2016
 *
 */
public class SimpleWebBrowserDemo extends Application{//javaFX specific	requires start method!!
	
	private static final String BASE_RESOURCE_CLASS = "danfulea.resources.DanfuleaResources";
	private ResourceBundle resources;
	
	Stage window;
	private WebEngine myWebEngine;		
	
	/**
	 * Initiates GUI and Events
	 */
	private void initGuiAndEvents(){
		resources = ResourceBundle.getBundle(BASE_RESOURCE_CLASS);
		
		window.setTitle("my  web browser");
		window.getIcons().add(new Image(this.resources.getString("form.icon.url")));//"/path/to/stackoverflow.jpg"));
		
		//action using lambda expression
		//window.setOnCloseRequest(e-> {
        //	e.consume();//this prevent for always closing the window if user choose NO in a dialog box
        //	closeProgram();
        //});
		//same thing without lambda expression
		window.setOnCloseRequest(new EventHandler<WindowEvent>(){
			public void handle(WindowEvent event) {
				event.consume();//this prevent for always closing the window if user choose NO in a dialog box
				closeProgram();				
			}			
		});
		
		
		final TextField addressBar = new TextField();//e.g. http://google.com
		// action
		addressBar.setOnAction(new EventHandler<ActionEvent>() {
			public void handle(ActionEvent event) {
				String url = basicURLcheck(addressBar.getText());
				myWebEngine.load(url);//addressBar.getText());
			}
		});

		//web browser component, javafx build-in
		WebView browser = new WebView();
		browser.setMinSize(1000, 700);//its size!!!
		myWebEngine = browser.getEngine();
		// webEngine.load("http://google.com");
		//handling errors
		myWebEngine.getLoadWorker().exceptionProperty().addListener(new ChangeListener<Throwable>() {
		
			public void changed(ObservableValue<? extends Throwable> observableValue, Throwable oldException,
					Throwable exception) {
				System.out.println("WebView encountered an exception loading a page: " + exception);
			}
		});

		//the layout
		VBox root = new VBox();//space out with 20 pixels:new VBox(20)
		root.getChildren().setAll(addressBar, browser);
		window.setScene(new Scene(root));
		//window.setScene(new Scene(root, 1000, 700));
		window.show();
	}
	
	/**
	 * Closing program
	 */
	private void closeProgram(){
		//System.out.println("Closing program!");
		boolean result = AlertBox.display(resources.getString("dialog.exit.title"), 
				resources.getString("dialog.exit.message"));
		if (result)
			window.close();
	}
	
	/**
	 * Checks if http suntax is present in request and append it if it is not.
	 * @param url
	 * the request url
	 * @return the completed request
	 */
	private String basicURLcheck(String url){
		String str = url;
		if (!url.startsWith("http://"))
			str = "http://" + url;
		
		return str;
	}

	@Override
	public void start(Stage stage) throws Exception {// called from main launch(args)
		// in java fx entire window is called Stage! The content inside window is called the Scene
		//notice, it does not allow standard constructor! This is the constructor!!!!
		window = stage;	
		initGuiAndEvents();
		
	}
	
	//testing
	public static void main(String[] args) {
		
		SimpleWebBrowserDemo.launch(args);// Application method used for launching a standalone app!
						// It does call the start method above!!!
	}
	//end testing
}
