package danfulea.db;

import java.awt.Color;
import java.awt.Component;
import java.awt.event.MouseEvent;
import java.sql.Connection;
import java.sql.Types;

import javax.swing.JComboBox;
import javax.swing.JLabel;
import javax.swing.JTable;
import javax.swing.JViewport;
import javax.swing.table.JTableHeader;
import javax.swing.table.TableCellRenderer;
import javax.swing.table.TableColumn;

import danfulea.math.Convertor;

/**
 * <p>
 * This class extends the functionality of the core DatabaseAgent. It is
 * intended to be as fast as possible for database manipulations by keeping
 * these operations at Database Management System (DMS like derby, MySQL,etc.)
 * level. This includes:
 * </p>
 * <p>
 * No updating id column to be always in order 1,2,3...
 * </p>
 * <p>
 * No fancy sorting data using JTable header
 * </p>
 * <p>
 * Working with tables having an auto-increment primary key column or  with tables having reserved one column of integer type for sole purpose 
 * of mimicking the missing auto-increment primary key column. Usually this column is named ID or NRCRT or Unique_ID etc. There are very few exception 
 * for this rule: when no single column act as unique but a combination of columns as in the case of nested tables which have no nested tables themselves. 
 * Moreover, if nested tables are huge as in the case of HpGe spectrum details then makes a lot o sense to not have a valid primary key. In this case, insertAll 
 * method is more appropriate to be used. 
 * </p>
 * <p>
 * Other feature is GUI support for viewing or sorting data, i.e JTable,
 * JCombobox and JLabel.
 * </p>
 * <p>
 * This class can be extended to provide a new layer of functionality such as
 * working with java objects, i.e all database operations are directly related
 * to java objects making the proper conversions to the required arrays (ORM, Object Relational Mapping technique). A
 * common scenario is to use the java class members as database table column
 * names. Therefore, build the column names array by accessing class member via
 * java reflection API (e.g. TestClass.class.getDeclaredFields();) and save/retrieve data via setters and getters.
 * Spring or Hibernate are java frameworks which are designed to work with java reflection and uses ORM technique.
 * BUT, according to Oracle documentation, the java reflection technique should be avoided in 
 * sections of code which are called frequently in performance-sensitive 
 * applications. Hence, we also discourage such techniques. In addition, modern 
 * applications uses GUI so why on earth do we want to fetch data from GUI fields
 * (e.g. JTextField) create an object out of it and save in database (ORM technique)
 * when we can simply save raw data directly?
 * </p>
 * 
 * @author Dan Fulea, 04 AUG. 2016
 *
 */
public class DatabaseAgentSupport{

	/**
	 * Connection object to a database.
	 */
	private Connection con = null;
	
	private static final int insertIndex = 0;
	private static final int updateIndex = 1;
	private static final int deleteIndex = 2;
	private static final int deleteAllIndex = 3;

	/**
	 * stores the index reflecting what database operation was used when
	 * performing select method to display data.
	 */
	private int dbOperationIndex = -1;

	/**
	 * stores current row selection in JTable to be preserved (highlighted)
	 * after update method.
	 */
	private int selectedRow = 0;

	/**
	 * sorts data by different columns
	 */
	private JComboBox<String> orderbyCb;

	/**
	 * holds records count
	 */
	private JLabel recordsCount = new JLabel();

	/**
	 * table holding the data
	 */
	private JTable mainTable;

	/**
	 * table model for the JTable which holds the data
	 */
	private final SimpleTableModel tableModel = new SimpleTableModel();

	/**
	 * The array of column names WITHOUT the auto-increment (primary key) column
	 * name
	 */
	private String[] usefullColumnNameArray;

	/**
	 * The array of column types WITHOUT the auto-increment (primary key) column
	 * type
	 */
	private Integer[] usefullColumnTypeArray;
	
	/**
	 * The number of columns in table WITHOUT the auto-increment (primary key) column
	 */
	private int usefullColumnCount=0;
	
	/**
	 * The array of column names 
	 */
	private String[] allColumnNameArray;

	/**
	 * The array of column types
	 */
	private Integer[] allColumnTypeArray;
	
	/**
	 * The number of columns in table
	 */
	private int allColumnCount=0;
	
	/**
	 * Database table.
	 */
	private String tableName="";
	
	/**
	 * Stores the auto-incremented primary key column name. The default value is
	 * "ID".
	*/
	private String primaryKey = "ID";

	/**
	 * stores the index of primary key column to be used when building
	 * usefullColumnNameArray. Default value is 0, i.e. the first column in
	 * table.
	 */
	private int primaryKeyColumnIndex = 0;//-1;//=>TEST ok -> updating in init
	
	/**
	 * If table was created with auto-increment option this is true.
	 */
	private boolean hasValidAIColumn =true;
	
	/**
	 * For backward compatibility with tables having no auto-increment column (not created in such way), mimicAI
	 * simulates the value of auto-incremented like column to be taken into account when performing INSERT operation.
	 */
	private int mimicAI = 0;
	

	/**
	 * Link column name used in a WHERE clause for a specific select operation.
	 * For example: ID=1; ID=IDLink, 1=IDValue.
	 */
	private String IDLink = null;

	/**
	 * Link value used in a WHERE clause for a specific select operation. For
	 * example: ID=1; ID=IDLink, 1=IDValue.
	 */
	private String IDValue = null;
	
	/**
	 * The number of rows count
	 */
	private int rowsCount=0;

	/**
	 * Constructor. DatabaseAgentSupport is related to a specific table in database. The name of the table as well as its primary key (fake or not) 
	 * must therefore be known. The connection to database must also be known in order to perform SQL tasks.
	 * 
	 * @param con
	 *            The connection to the database
	 * @param primaryKey
	 *            the primary key column name	 
	 * @param tableName the table name from database.
	 */
	public DatabaseAgentSupport(Connection con, String primaryKey,	String tableName) {	
		this.con = con;
		this.primaryKey = primaryKey;
		//this.primaryKeyColumnIndex = primaryKeyColumnIndex;
		this.tableName=tableName;
	}
		
	/**
	 * Return the rows count.
	 * @return the result
	 */
	public int getRowsCount(){
		return rowsCount;
	}
	
	/**
	 * Return the database table name this agent is pointing to.
	 * @return the result
	 */
	public String getDatabaseTableName(){
		return this.tableName;
	}
	
	/**
	 * Return the unique id column name (usually the primaryKey) of database table this agent is pointing to.
	 * @return the result
	 */
	public String getPrimaryKeyColumnName(){
		return this.primaryKey;
	}
	
	/**
	 * Return the connection to the database used by this agent
	 * @return the result
	 */
	public Connection getConnection(){
		return this.con;
	}
	
	/**
	 * Return the number of columns in table WITHOUT the auto-increment (primary key) column. <br>
	 * Useful when constructing column_values for insert/update without having to remember the 
	 * column count. E.g. column_values = new String[getUsefullColumnCount()]. 
	 * @return the result
	 */
	public int getUsefullColumnCount(){
		return usefullColumnCount;
	}
	
	/**
	 * Return the number of columns in table.<br>
	 * Useful when constructing column_values for insert/update without having to remember the 
	 * column count. E.g. column_values = new String[getAllColumnCount()]. 
	 * @return the result
	 */
	public int getAllColumnCount(){
		return allColumnCount;
	}
	
	/**
	 * 
	 * @return true if table was created with auto-increment column.
	 */
	public boolean isHasValidAIColumn() {
		return hasValidAIColumn;
	}

	/**
	 * If table was created with auto-increment column, hasValidAIColumn must be set to true. If set to false, then the fake primaryKey will 
	 * act as if it would be auto-incremented.
	 * @param hasValidAIColumn hasValidAIColumn
	 */
	public void setHasValidAIColumn(boolean hasValidAIColumn) {
		this.hasValidAIColumn = hasValidAIColumn;
	}
			
	/**
	 * 
	 * @return the index of auto-incremented (usually primary key) column
	 */
	public int getPrimaryKeyColumnIndex(){
		return this.primaryKeyColumnIndex;
	}

	/**
	 * Sets links used in a WHERE clause for a specific select operation. For
	 * example: ID=1; ID=IDLink, 1=IDValue. If links are not null, then the
	 * select operation with WHERE clause is used!
	 * 
	 * @param IDLink
	 *            the column name link
	 * @param IDValue
	 *            the value corresponding to column name link
	 */
	public void setLinks(String IDLink, String IDValue) {
		this.IDLink = IDLink;
		this.IDValue = IDValue;
	}

	/**
	 * Reset links to NO links
	 */
	public void resetLinks() {
		IDLink = null;
		IDValue = null;
	}

	/**
	 * Sets row in JTable to be highlighted after update method.
	 * 
	 * @param selectedRow
	 *            the selected row
	 */
	public void setSelectedRow(int selectedRow) {
		this.selectedRow = selectedRow;
	}

	/**
	 * Call this method to perform the first selection in database table and to
	 * generate basic stuff such as: the sorting JCombobox and the main JTable
	 *	 
	 */
	@SuppressWarnings("serial")
	public void init(){
		//orderby is a must for creating mimicAI
		String command = "select * from " + tableName + " ORDER BY " + primaryKey + " ASC";
		if (IDLink != null && IDValue != null) {
			command = "select * from " + tableName + " WHERE " + IDLink + "=" + IDValue + " ORDER BY "
					+ primaryKey + " ASC";
		}

		try {
			
			DatabaseAgent.select(con,command);
			
		} catch (Exception ex) {
			ex.printStackTrace();
		}
		rowsCount = DatabaseAgent.getRowCount();
		//1 column is reserved for autoincremented column.
		usefullColumnNameArray = new String[DatabaseAgent.getColumnCount() - 1];																	
		usefullColumnTypeArray = new Integer[DatabaseAgent.getColumnCount() - 1];
		
		usefullColumnCount = DatabaseAgent.getColumnCount() - 1;
		
		allColumnNameArray = new String[DatabaseAgent.getColumnCount()];																	
		allColumnTypeArray = new Integer[DatabaseAgent.getColumnCount()];
		
		allColumnCount = DatabaseAgent.getColumnCount();
		
		String[] columnNames = new String[DatabaseAgent.getColumnCount()];

		// loop through columns building usefullColumnNameArray and autoincremneted primarykey on the fly
		int k = 0;
		for (int i = 0; i < DatabaseAgent.getColumnCount(); i++) {// without id
			columnNames[i] = (String) DatabaseAgent.getColumnNames().elementAt(i);
			allColumnNameArray[i] = columnNames[i];
			allColumnTypeArray[i] = (Integer) DatabaseAgent.getColumnType().elementAt(i);
			if (columnNames[i].equals(primaryKey)){
				primaryKeyColumnIndex = i;
			} else {
				usefullColumnNameArray[k] = columnNames[i];
				usefullColumnTypeArray[k] = (Integer) DatabaseAgent.getColumnType().elementAt(i);
				k++;
			}	
		}		
		
		// now construct the JComboBox
		orderbyCb = new JComboBox<String>(columnNames);
		orderbyCb.setSelectedItem(primaryKey);

		// build the main JTable
		mainTable = new JTable(tableModel) {
			// tootip for cells in case its text is out of wiew
			// Implement table cell tool tips.
			public String getToolTipText(MouseEvent e) {
				String tip = null;
				java.awt.Point p = e.getPoint();
				int rowIndex = rowAtPoint(p);
				int colIndex = columnAtPoint(p);

				try {
					tip = getValueAt(rowIndex, colIndex).toString();
				} catch (RuntimeException e1) {
					// catch null pointer exception if mouse is over an empty
					// line
				}

				return tip;
			}
			
			//Implement table header tool tips.
		    protected JTableHeader createDefaultTableHeader() {
		        return new JTableHeader(columnModel) {
		            public String getToolTipText(MouseEvent e) {
		                //String tip = null;
		                java.awt.Point p = e.getPoint();
		                int index = columnModel.getColumnIndexAtX(p.x);
		                int realIndex = 
		                        columnModel.getColumn(index).getModelIndex();
		                return allColumnNameArray[realIndex];//columnToolTips[realIndex];
		            }
		        };
		    }
		    
		    //try to auto-resize column width to fit data and also alternate row colors
		    @Override
		       public Component prepareRenderer(TableCellRenderer renderer, int row, int column) {
		           Component component = super.prepareRenderer(renderer, row, column);
		           //==============================================================================
		           Color alternateColor = new Color(173,216,230);//lightblue
		           Color whiteColor = Color.WHITE;//#######################
		           if (!component.getBackground().equals(getSelectionBackground())){
			           Color bg = (row % 2 == 0 ? alternateColor : whiteColor);
			           component .setBackground(bg);
			         bg = null;
			        }
		           //==============================================================================
		           int rendererWidth = component.getPreferredSize().width;
		           TableColumn tableColumn = getColumnModel().getColumn(column);
		           tableColumn.setPreferredWidth(Math.max(rendererWidth + getIntercellSpacing().width, tableColumn.getPreferredWidth()));
		           return component;
		        }		   
		    
		    //the following are needed for the case when a table has width larger than its parent (many columns) then EXCEPT the case 
		    //we use AUTO_RESIZE_OFF the horizontal scrollbar does not appear and all columns are cluttered together in a small space. 
		    //That's because they are considered in auto-mode and STRANGELY enough, data width is ignored.
		    //But using AUTO_RESIZE_OFF then if parent width is larger then the table the table does not automatically fill the available 
		    //space.
		    /*@Override
		    public Dimension getPreferredSize() {//this is not needed since we change AUTORESIZEOFF for much better performance!!!
		      if (getParent () instanceof JViewport) {
		        if (
		          ((JViewport) getParent()).getWidth() > super.getPreferredSize().width)
		         {
		          return getMinimumSize();
		        }
		      }
		      return super.getPreferredSize(); 
		    }*/
		    		    
		    @Override
		    public boolean getScrollableTracksViewportWidth () {
		    	if (autoResizeMode == AUTO_RESIZE_OFF) {//if (autoResizeMode != AUTO_RESIZE_OFF) {
		    		if (getParent() instanceof JViewport) {
		    			return (((JViewport) getParent()).getWidth() > getPreferredSize().width);
		    			//if true=>the width of the Viewport determine the width of the table
		    		}
		    		return true;
		      }
		      return true;//false;
		    }
		};

		// some customization for main JTable===============================
		// make sure vertical column lines are displayed
		mainTable.setShowVerticalLines(true);
		// do not move columns by dragging:
		mainTable.getTableHeader().setReorderingAllowed(false);
		// auto-resize=>Don't care about dynamic resizing since it will involve
		// more time-consuming DB operations!!!!
		//to work with my custom suto-resize
		mainTable.setAutoResizeMode(JTable.AUTO_RESIZE_OFF);//.setAutoResizeMode(JTable.AUTO_RESIZE_ALL_COLUMNS);//
		mainTable.setAutoCreateColumnsFromModel(true);
		
		// try to alterante row color THE BELOW DOESNT WORK WELL=>ON SELECT ROW ...NOT VISIBLE
		/*mainTable.setDefaultRenderer(Object.class, new DefaultTableCellRenderer() {

	        @Override
	        public Component getTableCellRendererComponent(JTable table, 
	                Object value, boolean isSelected, boolean hasFocus,
	                int row, int column) {
	            Component c = super.getTableCellRendererComponent(table, 
	                value, isSelected, hasFocus, row, column);
	            c.setBackground(row%2==0 ? Color.white : Color.lightGray);//Color.yellow);                        
	            return c;
	        };
	    });*/
		//end customization===========================================================

		// set table to data taken from database.
		// getData() and getColumnNames() are available because we performed the select method
		tableModel.setDataVector(DatabaseAgent.getData(), DatabaseAgent.getColumnNames());
		// setup the record count JLabel getRowCount() is available because we performed the select method.
		recordsCount.setText("" + rowsCount);
		
		//mimic autoincremented column. NOTICE, if deleteall then mimicAI is restarted with 0 (but only after restart app)!!!
		if (rowsCount>0){
			//mimicAI = (int)DatabaseAgent.getValueAt(rowsCount-1, primaryKeyColumnIndex);
			mimicAI = ((Integer)DatabaseAgent.getValueAt(rowsCount-1, primaryKeyColumnIndex)).intValue();
		} else {
			mimicAI=0;
		}	

		// TESTING=====================
		//System.out.println("mimicAI= "+mimicAI);
		// for (int i = 0; i< usefullColumnNameArray.length; i++)
		// System.out.print(usefullColumnNameArray[i]+", ");//WORKS
		// ==================================
	}

	/**
	 * 
	 * @return the sorting GUI element
	 */
	public JComboBox<String> getOrderByComboBox() {
		return orderbyCb;
	}

	/**
	 * 
	 * @return the records count GUI element
	 */
	public JLabel getRecordsLabel() {
		return recordsCount;
	}

	/**
	 * 
	 * @return the JTable object holding data from database
	 */
	public JTable getMainTable() {
		return mainTable;
	}

	/**
	 * Performs selection in database table and sorts data by desired column. If
	 * link members (IDLink and IDValue) are not null, then a SELECT query with
	 * WHERE clause is performed. Otherwise, the SELECT all from database is
	 * used. Of course, ORDER BY clause is presented in either cases.
	 * 
	 * @param orderbyS
	 *            the ORDER BY column name
	 */
	public void performSelection(String orderbyS) {
	//(String table, String orderbyS) {
		String str = "select * from " + tableName + " ORDER BY " + orderbyS + " ASC";
		if (IDLink != null && IDValue != null) {
			str = "select * from " + tableName + " WHERE " + IDLink + "=" + IDValue + " ORDER BY " + orderbyS + " ASC";
		}

		try {
			//super.select(str);
			DatabaseAgent.select(con, str);
		} catch (Exception ex) {
			ex.printStackTrace();
		}
		rowsCount = DatabaseAgent.getRowCount();
		//int rowCount = DatabaseAgent.getRowCount();
		recordsCount.setText("" + rowsCount);// records

		tableModel.setDataVector(DatabaseAgent.getData(), DatabaseAgent.getColumnNames());// table

		// when coming from update, no last row selection
		if (dbOperationIndex == updateIndex) {
			mainTable.setRowSelectionInterval(selectedRow, selectedRow);
			dbOperationIndex = -1;// reset to avoid error when no updating
									// dbindex in operation like sorting via
									// JComboBox
			return;
		}

		// By default, select last row!
		if (rowsCount > 0) {
			mainTable.setRowSelectionInterval(rowsCount - 1, rowsCount - 1);
		}
	}

	/**
	 * Generic SELECT query is performed here. It can be used for searching in
	 * database table.
	 * @param command the generic SQL select query
	 */
	public void select(String command) {
		try {
			DatabaseAgent.select(con, command);
		} catch (Exception ex) {
			ex.printStackTrace();
		}
		rowsCount = DatabaseAgent.getRowCount();
		//int rowCount = DatabaseAgent.getRowCount();
		recordsCount.setText("" + rowsCount);// records

		tableModel.setDataVector(DatabaseAgent.getData(), DatabaseAgent.getColumnNames());// table

		// By default, select last row!
		if (rowsCount > 0) {
			mainTable.setRowSelectionInterval(rowsCount - 1, rowsCount - 1);
		}
	}
	
	/**
	 * Inserts a row in database table using all columns 
	 * and all their values. The value passed here must be in order relative to 
	 * the columns of the table. For instance a table with columns: ID, UserName, Password 
	 * the value array to be inserted must be: {ID_value, userName_value, password_value}. <br>
	 * 
	 * SOMETIMES IS USEFUL TO NOT HAVE A VALID PRIMARY KEY COLUMN SINCE A COMBINMATION OF COLUMNS IS UNIQUE. 
	 * THIS IS THE CASE FOR NESTED TABLES WHICH DOES NOT HAVE ANY NESTED TABLES TTHEMSELVES. 
	 * MOREOVER, IT IS USEFUL WHEN NESTED TABLE HAS HUGE ROWS (E.G. HpGe SPECTRUM DETAIL OF CHANNELS AND 
	 * PULSES WHICH CAN HAVE 16000 ROWS OR MORE). IN THIS CASE IT IS POINTLESS TO HAVE A COLUMN 
	 * AS PRIMARY KEY SINCE ROWNUMBER COLUMN PLUS LINK COLUMN (LINK TO MAIN SPECTRUM TABLE) IS UNIQUE.  
	 * 
	 * @param value
	 *            the array of values to be inserted	 
	 */
	public void insertAll(String[] value){//, String orderbyS) {
		insert(tableName, allColumnNameArray, allColumnTypeArray, value);
		
		dbOperationIndex = insertIndex;
		//performSelection(orderbyS);//BEST LET CALLING PROGRAM DO THIS IF NECESSARY!!!
	}

	/**
	 * Inserts a row in database table using columns
	 * and their values other than primary key column name and value. The value passed here must be in order relative to 
	 * the columns of the table. For instance a table with columns: ID, UserName, Password having ID as auto-incremented
	 * primary key, the value array to be inserted must be: {userName_value, password_value}
	 * 
	 * @param value
	 *            the array of values to be inserted	 
	 */
	public void insert(String[] value){//, String orderbyS) {
		if (hasValidAIColumn)
			insert(tableName, usefullColumnNameArray, usefullColumnTypeArray, value);
		else{
			//some work to do
			mimicAI=mimicAI+1;//mimicAI is increased!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!and stay increased!!!!
			int size = usefullColumnNameArray.length;
			String[] newValues = new String[size+1];
			String[] columns = new String[size+1];
			Integer[] columnTypes = new Integer[size+1];
			for (int i=0;i<size;i++){
				newValues[i]=value[i];
				columns[i]=usefullColumnNameArray[i];
				columnTypes[i]=usefullColumnTypeArray[i];
			}
			columns[size]=primaryKey;
			columnTypes[size]=Types.INTEGER;
			newValues[size]=Convertor.intToString(mimicAI);
			
			insert(tableName, columns, columnTypes, newValues);
		}
		dbOperationIndex = insertIndex;
		//performSelection(orderbyS);//BEST LET CALLING PROGRAM DO THIS IF NECESSARY!!!
	}

	/**
	 * Inserts a row in database table.
	 * 
	 * @param tabel
	 *            the table where to perform this operation
	 * @param cvalue
	 *            the array of column names
	 * @param tvalue
	 *            the array of column types
	 * @param value
	 *            the array of values to be inserted	 
	 */
	private void insert(String tabel, String cvalue[], Integer tvalue[], String[] value) {
		try {
			DatabaseAgent.insert(con,tabel, cvalue, tvalue, value);
		} catch (Exception e) {
			e.printStackTrace();
		}		
	}
		
	/**
	 * Updates a row using columns and their values
	 * other then primary key column name and value where primary key column has
	 * "IDValue" value. The value passed here must be in order relative to 
	 * the columns of the table. For instance a table with columns: ID, UserName, Password having ID as auto-incremented
	 * primary key, the value array to be updated must be: {userName_value, password_value}
	 * 
	 * @param value
	 *            the array of values to be updated
	 * @param IDValue
	 *            the reference value for corresponding primary key column
	 */
	public void update(String[] value, String IDValue){//, String orderbyS) {
		update(tableName, usefullColumnNameArray, usefullColumnTypeArray, value, 
				primaryKey, IDValue);
		dbOperationIndex = updateIndex;
		//performSelection(orderbyS);
	}
	
	/**
	 * Updates a row where a specific IDLink column
	 * name has a specific value IDValue.
	 * 
	 * @param tabel
	 *            the table where to perform this operation
	 * @param cvalue
	 *            the array of column names where values will be updated
	 * @param tvalue
	 *            the array of column types where values will be updated
	 * @param value
	 *            the array of values to be updated
	 * @param IDLink
	 *            the reference column name used as link
	 * @param IDValue
	 *            the reference value for corresponding IDLink column
	 */
	public void update(String tabel, String[] cvalue, Integer tvalue[], String[] value,
			String IDLink, String IDValue) {
		try {
			DatabaseAgent.update(con,tabel, cvalue, tvalue, value, IDLink, IDValue);
		} catch (Exception e) {
			e.printStackTrace();
		}		
		dbOperationIndex = updateIndex;
	}

	/**
	 * Deletes all data and displays nothing.
	 * 
	 */
	public void deleteAll(){//String orderbyS) {
		try {
			DatabaseAgent.deleteAll(con,tableName);
		} catch (Exception e) {
			e.printStackTrace();
		}
		dbOperationIndex = deleteAllIndex;
		//performSelection(orderbyS);// needed to display empty
	}
	
	/**
	 * Deletes multiple rows where primary key column has
	 * "IDValue" value.
	 * 
	 * @param IDValue
	 *            the reference value array for corresponding primary key column	 
	 */
	public void delete(String[] IDValue){//, String orderbyS) {
		int n = IDValue.length;
		String[] IDLink = new String[n];
		for (int i=0; i<n;i++){
			IDLink[i] = primaryKey;//based on priKey
		}
		try {
			DatabaseAgent.delete(con, tableName, IDLink, IDValue);
		} catch (Exception e) {
			e.printStackTrace();
		}	
		dbOperationIndex = deleteIndex;
		//performSelection(orderbyS);
	}

	/**
	 * Deletes a row where primary key column has
	 * "IDValue" value.
	 * 
	 * @param IDValue
	 *            the reference value for corresponding primary key column	 
	 */
	public void delete(String IDValue){//, String orderbyS) {
		
		delete(tableName, primaryKey, IDValue);
		dbOperationIndex = deleteIndex;
		//performSelection(orderbyS);
	}

	/**
	 * Deletes a row where a specific IDLink column
	 * name has a specific value IDValue.
	 * 
	 * @param tabel
	 *            the table where to perform this operation
	 * @param IDLink
	 *            the reference column name used as link
	 * @param IDValue
	 *            the reference value for corresponding IDLink column	 
	 */
	public void delete(String tabel, String IDLink, String IDValue) {
		try {
			DatabaseAgent.delete(con, tabel, IDLink, IDValue);
		} catch (Exception e) {
			e.printStackTrace();
		}		
		dbOperationIndex = deleteIndex;
	}

	///**
	// * This is mainly used for delete rows from nested tables. If 2 or more
	// * tables are linked through the IDLink column name then deleting data from
	// * one table (main) leads to delete the corresponding data from other
	// * (nested) table. This is done using WHERE clause, e.g. WHERE ID=1. Common
	// * usage: inside a delete row from main table method, use deleteNested
	// * method to also delete rows in nested tables.
	 //* 
	// * @param tabel
	// *            the table where to perform this operation
	// * @param IDLink
	// *            the reference column name used as link
	// * @param IDValue
	// *            the reference value for corresponding IDLink column
	// */
	//public void deleteNested(String tabel, String IDLink, String IDValue) {// ,
																			// String
																			// orderbyS){
	//	try {
	//		DatabaseAgent.delete(con,tabel, IDLink, IDValue);
	//	} catch (Exception e) {
	//		e.printStackTrace();
	//	}
	//	dbOperationIndex = deleteIndex;		
	//}
	
	/**
	 * Return the current auto-incremented primary key value
	 * @return the result
	 */
	public int getAIPrimaryKeyValue(){
		int result = -1;
		String s = "select MAX("+primaryKey+") from "+tableName;
		try{
			DatabaseAgent.select(con, s);
			result =(Integer) DatabaseAgent.getValueAt(0, primaryKeyColumnIndex);//ID	
		}catch (Exception ex){
			ex.printStackTrace();
		}
		 return result;
	}

}
