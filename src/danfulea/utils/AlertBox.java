package danfulea.utils;

import java.util.ResourceBundle;

import javafx.event.ActionEvent;
import javafx.event.EventHandler;
import javafx.geometry.Insets;
import javafx.geometry.Pos;
import javafx.scene.Scene;
import javafx.scene.control.Button;
import javafx.scene.control.Label;
import javafx.scene.image.Image;
import javafx.scene.layout.Background;
import javafx.scene.layout.BackgroundFill;
import javafx.scene.layout.CornerRadii;
import javafx.scene.layout.HBox;
import javafx.scene.layout.VBox;
import javafx.scene.paint.Color;
import javafx.stage.Modality;
import javafx.stage.Stage;

/**
 * A simple dialog box with two yes/no options (it uses fx technology).
 * 
 * @author Dan Fulea, 12 AUG.2016
 *
 */
public class AlertBox {

	private static final String BASE_RESOURCE_CLASS = "danfulea.resources.DanfuleaResources";
	public static Color bkgColor = Color.web("#E6FFD2");// new Color(230, 255, 210, 255);// Linux mint
	private static ResourceBundle resources;
	private static boolean answer;
    
	/**
	 * Display the dialog box
	 * @param title
	 * its title
	 * @param message
	 * its message
	 * @return true or false depending on what user choose (yes or no)
	 */
	public static boolean display(String title, String message) {
		resources = ResourceBundle.getBundle(BASE_RESOURCE_CLASS);
        final Stage window = new Stage();//creates new stage

        //Block events to other windows
        window.initModality(Modality.APPLICATION_MODAL);
        window.setTitle(title);
        window.getIcons().add(new Image(resources.getString("form.icon.url")));
        window.setMinWidth(250);

        Label label = new Label();
        label.setText(message);          
      
        Button yesButton = new Button(resources.getString("dialog.exit.yes"));//"Yes");
        Button noButton = new Button(resources.getString("dialog.exit.no"));//"No");

        //Clicking will set answer and close window
        //yesButton.setOnAction(e -> {
         //   answer = true;
          //  window.close();
        //});
        //without lambda expression because we like keep things clear not necessarly simple
        yesButton.setOnAction(new EventHandler<ActionEvent>() {
			public void handle(ActionEvent event) {
				answer = true;
				window.close();
			}
		});
        //noButton.setOnAction(e -> {
         //   answer = false;
          //  window.close();
        //});
        noButton.setOnAction(new EventHandler<ActionEvent>() {
			public void handle(ActionEvent event) {
				answer = false;
				window.close();
			}
		});
       
        
        HBox hlay = new HBox(20);
        hlay.getChildren().addAll(yesButton, noButton);//closeButton);
        hlay.setAlignment(Pos.CENTER);
        hlay.setBackground(new Background(new BackgroundFill(bkgColor, CornerRadii.EMPTY, Insets.EMPTY)));


        
        VBox layout = new VBox(20);
        layout.getChildren().addAll(label, hlay);//closeButton);
        layout.setAlignment(Pos.CENTER);
        layout.setBackground(new Background(new BackgroundFill(bkgColor, CornerRadii.EMPTY, Insets.EMPTY)));
        
        //Display window and wait for it to be closed before returning
        Scene scene = new Scene(layout);
        window.setScene(scene);
        window.showAndWait();
        
        return answer;
    }
}
