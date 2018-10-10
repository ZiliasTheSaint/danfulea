package danfulea.utils.table;

import javax.swing.table.DefaultTableCellRenderer;

/**
 * Handle Jtable cell number display. <br>
 * 
 * @author Dan Fulea, 15 May 2011
 * 
 */
@SuppressWarnings("serial")
public class UnformattedCellRenderer extends DefaultTableCellRenderer{
	
	/**
	 * Constructor
	 */
	public UnformattedCellRenderer(){
		super();
	}
	
	//override
	public void setValue(Object value){
		//display the value as it is...not formatted!!!!
		setText(value.toString());
	}

}
