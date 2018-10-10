package danfulea.utils;

import javax.swing.DefaultListModel;
import javax.swing.JList;

/**
 * List utilities class for handling operations with Strings.
 * 
 * @author Dan Fulea, 01 FEB. 2005
 *
 */
public class ListUtilities {

	/**
	 * Adds a String to list
	 * 
	 * @param s
	 *            the String
	 * @param dlm
	 *            the default list model
	 */
	@SuppressWarnings({ "rawtypes", "unchecked" })
	public static void add(String s, DefaultListModel dlm) {
		dlm.addElement(s);
	}

	/**
	 * Removes all elements from list
	 * 
	 * @param dlm
	 *            the default list model
	 */
	@SuppressWarnings("rawtypes")
	public static void removeAll(DefaultListModel dlm) {
		dlm.clear();
	}

	/**
	 * Gets the list selected value as a String
	 * 
	 * @param list
	 *            the list
	 * @return the String representation of selected value
	 */
	@SuppressWarnings("rawtypes")
	public static String getSelectedValueAsString(JList list) {
		return (String) list.getSelectedValue();
	}

	// ------------------------
	/**
	 * Selects an entry in the list based on its index
	 * 
	 * @param index
	 *            the index
	 * @param list
	 *            the list
	 */
	@SuppressWarnings("rawtypes")
	public static void select(int index, JList list) {
		list.setSelectedIndex(index);
	}

	/**
	 * Gets the list selected index
	 * 
	 * @param list
	 *            the list
	 * @return the index
	 */
	@SuppressWarnings("rawtypes")
	public static int getSelectedIndex(JList list) {
		return list.getSelectedIndex();
	}

	/**
	 * Removes an entry based on its index
	 * 
	 * @param index
	 *            the index
	 * @param dlm
	 *            the default list model
	 */
	@SuppressWarnings("rawtypes")
	public static void remove(int index, DefaultListModel dlm) {
		dlm.remove(index);
	}

	/**
	 * Gets the selected indexes as an array of ints (multiple list selction)
	 * 
	 * @param list
	 *            the list
	 * @return the array of selected indexes
	 */
	@SuppressWarnings("rawtypes")
	public static int[] getSelectedIndeces(JList list) {
		return list.getSelectedIndices();
	}

}
