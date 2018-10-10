package danfulea.utils.table;

import javax.swing.table.DefaultTableModel;

/**
 * The sortable table model.
 * Based on example code written by Joris Van den Bogaert see esus.com
 * @author Dan Fulea, 22 Mar. 2005
 */
@SuppressWarnings("serial")
public class SortableTableModel extends DefaultTableModel {
	int[] indexes;
	TableSorter sorter;

	/**
	 * Constructor
	 */
	public SortableTableModel() {
		super();
	}

	/**
	 * Get value stored at some row and column.
	 * @param row row
	 * @param col col
	 * @return the object stored at (row, col)
	 */
	public Object getValueAt(int row, int col) {
		int rowIndex = row;
		if (indexes != null) {
			rowIndex = indexes[row];
		}

		return super.getValueAt(rowIndex, col);		
	}

	/**
	 * Set value at some row and col
	 * @param value value
	 * @param row row
	 * @param col col
	 */
	public void setValueAt(Object value, int row, int col) {
		int rowIndex = row;
		if (indexes != null) {
			rowIndex = indexes[row];
		}

		super.setValueAt(value, rowIndex, col);
	}

	/**
	 * Sort by column
	 * @param column column
	 * @param isAscent isAscent
	 */
	public void sortByColumn(int column, boolean isAscent) {
		if (sorter == null) {
			sorter = new TableSorter(this);
		}

		sorter.sort(column, isAscent);
		fireTableDataChanged();
	}

	/**
	 * Rerieve row indexes  
	 * @return the result
	 */
	public int[] getIndexes() {
		int n = getRowCount();
		if (indexes != null) {
			if (indexes.length == n) {
				return indexes;
			}
		}

		indexes = new int[n];
		for (int i = 0; i < n; i++) {
			indexes[i] = i;
		}

		return indexes;
	}
}
