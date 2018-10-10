package danfulea.db;

import javax.swing.table.DefaultTableModel;

/**
 * Simple table model for the database viewer, which is a JTable object.
 * 
 * @author Dan Fulea, 04 AUG. 2016
 */
public class SimpleTableModel extends DefaultTableModel {

	private static final long serialVersionUID = 1L;

	/**
	 * Constructor. This is the table model used by the JTable object to display
	 * the data stored in database table.
	 */
	public SimpleTableModel() {
		super();
	}

	/**
	 * We work with data taken from database, hence we do not want table cells
	 * to be editable!
	 */
	public boolean isCellEditable(int row, int col) {
		return false;
	}
}
