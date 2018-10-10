package danfulea.db;

import java.io.ByteArrayInputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.io.StringReader;
import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.ResultSetMetaData;
import java.sql.SQLException;
import java.sql.Statement;
import java.sql.Types;
import java.util.Vector;

import danfulea.math.Convertor;

/**
 * <p>
 * The DatabaseAgent class handles basic database operations: connection,
 * creation or removal of database and tables, select to retrieve database data
 * in form of data and column names vectors, insert, delete and update data into
 * tables.
 * </p>
 *
 * <p>
 * Supports derby embedded, MySql and PostgreSQL connections. Supports most used
 * data types like: integer, double precision, varchar, clob, real, etc.
 * </p>
 * <p>
 * Be advise that MySQL connector is licensed under GNU-GPL copyleft license
 * which means that by using it, you can only create free software and not any
 * kind of software (proprietary software included)! Also, using derby embedded
 * means the derby engine cannot run in multiple instances of JVM, i.e. cannot
 * run two or more instances of the same application because the previous one
 * has already established a connection to embedded derby. To avoid this, use
 * MySQL DB or PostgreSQL DB or use Derby Network Server. However, using
 * embedded derby can be very useful when we do not want to install additional
 * software (such as MySQL DB, PostgreSQL DB etc) in order for our application
 * to run. Just use embedded derby.jar and you are good to go. The price is one
 * instance of application is allowed!
 * </p>
 * 
 * <p>
 * It is intended to be fast therefore the standard scenario is to use tables
 * having a column created with auto-increment and primary key option. All other
 * columns must have no indexes. Alternatively, you always have to have an integer type column which 
 * mimics the auto-incremented standard column.
 * Generally, using synchronized to control database access is not a good idea. 
 * The database has all kinds of internal mechanisms to control concurrent access, so it is best to use them instead.
 * </p>
 * 
 * @author Dan Fulea, 04 AUG. 2016
 * 
 */

public class DatabaseAgent {

	public static final int DERBY_CONNECTION = 0;
	public static final int MYSQL_CONNECTION = 1;
	public static final int POSTGRESQL_CONNECTION = 2;

	/**
	 * Connection ID pointing to various database management systems (DMS) such
	 * as derby, MySQL or PostgreSQL. It points to embedded derby by default.
	 */
	public static int ID_CONNECTION = DERBY_CONNECTION;// common connection

	/**
	 * Specifies where embedded derby database is located relative to
	 * application folder. The default value is Data.
	 */
	public static String dataFolderUsedByDerbyEmbedded = "Data";// common derby
																// "workspace"

	//
	// Connection object used for all SQL operations. It is common for all
	// DatabaseAgent objects.->wrong->see below
	//The connection is 1 per database indeed so it seems appropriate to share it
	//between objects but what if we want to connect to multiple databases?
	//private static Connection con = null;

	/**
	 * Stores table column names.
	 */
	private static Vector<Object> colNumeVec;

	/**
	 * Stores table data.
	 */
	private static Vector<Object> dataVec;

	/**
	 * Stores table column types as Java class
	 */
	private static Vector<String> columnClass;

	/**
	 * Stores table column types as SQL types. This is used for operations like
	 * insert or update where PreparedStatement methods like setString or setInt
	 * are used.
	 */
	private static Vector<Integer> columnType;

	/**
	 * internal flag to check if derby is closed.
	 */
	private static boolean derbyClosed = true;

	/**
	 * Stores row count after a SELECT operation.
	 */
	private static int rowCount = 0;

	/**
	 * Stores column count after a SELECT operation.
	 */
	private static int colCount = 0;

	//**
	// * Stores the auto-incremented primary key column name. The default value is
	// * "ID".
	// */
	//private static String primaryKey = "ID";

	// ===========SOME UTILITIES================================================

	/**
	 * Creates derby database in Data folder (specified by
	 * dataFolderUsedByDerbyEmbedded) relative to the application folder.
	 * 
	 * @param dbName
	 *            database name to be created
	 * @param user
	 *            user name, can be null
	 * @param pass
	 *            password, can be null
	 */
	public static void createDerbyDatabaseInDataFolder(String dbName, String user, String pass) {

		Connection conng = null;

		String datas = dataFolderUsedByDerbyEmbedded;
		String currentDir = System.getProperty("user.dir");
		String file_sep = System.getProperty("file.separator");
		String opens = currentDir + file_sep + datas;

		opens = opens + file_sep + dbName;// e.g. : Data/dbname

		String protocol = "jdbc:derby:";// derby protocol
		String driver = "org.apache.derby.jdbc.EmbeddedDriver";

		try {
			Class.forName(driver).newInstance();
			conng = DriverManager.getConnection(protocol + opens + ";create=true", user, pass);

			if (conng != null) {
				conng.close();
			}

		} catch (InstantiationException e) {
			e.printStackTrace();
			return;
		} catch (IllegalAccessException e) {
			e.printStackTrace();
			return;
		} catch (ClassNotFoundException e) {
			e.printStackTrace();
			return;
		} catch (SQLException e) {
			e.printStackTrace();
			return;
		}

		System.out.println("Database created succssefully!");
	}

	/**
	 * Deletes a database. It doesn't work for embedded derby; in this case, you
	 * can manually delete the folder which contains the data!
	 * 
	 * @param dbName
	 *            database name to be deleted
	 * @param user
	 *            user name, can be null
	 * @param pass
	 *            password, can be null
	 */
	public static void deleteDatabase(String dbName, String user, String pass) {

		Connection conng = null;
		Statement s = null;

		String dbURL = "jdbc:mysql://localhost:3306/";
		String driver = "org.apache.derby.jdbc.EmbeddedDriver";

		try {
			if (ID_CONNECTION == DERBY_CONNECTION) {

				System.out.println("Cannot delete database using embedded driver."
						+ " Go to database folder and delete it manually!");
				return;
			} else if (ID_CONNECTION == MYSQL_CONNECTION) {
				dbURL = "jdbc:mysql://localhost:3306/";
				driver = "com.mysql.jdbc.Driver";

				Class.forName(driver).newInstance();

				String url = dbURL;// path to mysql server not to any specific
									// DB;
				// example of user: "root"
				conng = DriverManager.getConnection(url + "?user=" + user + "&password=" + pass);
			} else if (ID_CONNECTION == POSTGRESQL_CONNECTION) {
				dbURL = "jdbc:postgresql://localhost:5432/";
				driver = "org.postgresql.Driver";

				Class.forName(driver).newInstance();

				String url = dbURL;// path to server
				// example of user and password: "postgres" and "admin"
				conng = DriverManager.getConnection(url + "?user=" + user + "&password=" + pass);
			}

			conng.setAutoCommit(true);// false);//allow multiple operations
			// here we set to false because postgresql complains about it!!!
			s = conng.createStatement();

			s.execute("DROP DATABASE " + dbName);

			// conng.commit();//execute all operations. If autocommit is true,
			// then this is not needed!

			// ====================comments=================
			// conn.setAutoCommit(false);
			// statement.executeQuery(query);
			// statement.commit();//can also be conn.commit() if multiple
			// statements

			// is the same as:

			// conn.setAutoCommit(true);
			// statement.executeQuery(query);

			// If a connection is in auto-commit mode, then all its SQL
			// statements will be executed and committed as individual
			// transactions. Otherwise, its SQL statements are grouped into
			// transactions that are terminated by a call to either
			// the method commit or the method rollback.
			// commit = Makes all changes made since the previous
			// commit/rollback permanent and releases any database locks
			// currently held by this Connection object. This method should be
			// used only when auto-commit mode has been disabled.
			// rollback=Undoes all changes made in the current transaction and
			// releases any database locks currently held by this
			// Connection object. This method should be used only when
			// auto-commit mode has been disabled.

			// By default, new connections are in auto-commit mode! OK!
			// ===================end comments==================

			if (s != null)
				s.close();

			if (conng != null) {
				conng.close();
			}

		} catch (InstantiationException e) {
			e.printStackTrace();
			return;
		} catch (IllegalAccessException e) {
			e.printStackTrace();
			return;
		} catch (ClassNotFoundException e) {
			e.printStackTrace();
			return;
		} catch (SQLException e) {
			e.printStackTrace();
			return;
		}

		System.out.println("Database deleted!");
	}

	/**
	 * Deletes a table from database.
	 * 
	 * @param dbName
	 *            database name
	 * @param table
	 *            table name to be deleted
	 * @param user
	 *            user name, can be null
	 * @param pass
	 *            password, can be null
	 */
	public static void deleteTable(String dbName, String table, String user, String pass) {

		Connection conng = null;
		Statement s = null;

		String datas = dataFolderUsedByDerbyEmbedded;// default: Data folder
		String currentDir = System.getProperty("user.dir");
		String file_sep = System.getProperty("file.separator");
		String opens = currentDir + file_sep + datas;

		opens = opens + file_sep + dbName;

		String dbURL = "jdbc:mysql://localhost:3306/";// initialization
		String driver = "org.apache.derby.jdbc.EmbeddedDriver";// initialization

		try {
			if (ID_CONNECTION == DERBY_CONNECTION) {
				String protocol = "jdbc:derby:";
				driver = "org.apache.derby.jdbc.EmbeddedDriver";

				Class.forName(driver).newInstance();
				conng = DriverManager.getConnection(protocol + opens + ";create=false", user, pass);

			} else if (ID_CONNECTION == MYSQL_CONNECTION) {
				dbURL = "jdbc:mysql://localhost:3306/";
				driver = "com.mysql.jdbc.Driver";

				Class.forName(driver).newInstance();
				String url = dbURL + dbName;
				conng = DriverManager.getConnection(url + "?user=" + user + "&password=" + pass);

			} else if (ID_CONNECTION == POSTGRESQL_CONNECTION) {
				dbURL = "jdbc:postgresql://localhost:5432/";
				driver = "org.postgresql.Driver";

				Class.forName(driver).newInstance();
				String url = dbURL + dbName;
				conng = DriverManager.getConnection(url + "?user=" + user + "&password=" + pass);

			}

			conng.setAutoCommit(false);
			s = conng.createStatement();

			s.execute("drop table  " + table);

			conng.commit();

			if (s != null)
				s.close();

			if (conng != null) {
				conng.close();
			}

		} catch (InstantiationException e) {
			e.printStackTrace();
			return;
		} catch (IllegalAccessException e) {
			e.printStackTrace();
			return;
		} catch (ClassNotFoundException e) {
			e.printStackTrace();
			return;
		} catch (SQLException e) {
			e.printStackTrace();
			return;
		}

		System.out.println("Table deleted!");
	}

	/**
	 * Creates a table inside database.
	 * 
	 * @param dbName
	 *            database name
	 * @param table
	 *            table name to be created
	 * @param cname
	 *            table column names stored in an array of strings.
	 * @param ctype
	 *            table column types stored in an array of strings. Example:
	 *            {"integer", "DOUBLE PRECISION", "VARCHAR(59)", "CLOB"}. Note
	 *            that the corresponding CLOB type in PostgreSQL is TEXT.
	 * @param cindexes
	 *            table column indexes stored in an array of strings. To mark
	 *            the primary key column use "Y". Example: {"Y", "N", "N", "N"}.
	 * @param user
	 *            user name, can be null
	 * @param pass
	 *            password, can be null
	 */
	public static void createTable(String dbName, String table, String[] cname, String[] ctype, String[] cindexes,
			String user, String pass) {

		Connection conng = null;
		Statement s = null;

		String datas = dataFolderUsedByDerbyEmbedded;
		String currentDir = System.getProperty("user.dir");
		String file_sep = System.getProperty("file.separator");
		String opens = currentDir + file_sep + datas;

		opens = opens + file_sep + dbName;

		String dbURL = "jdbc:mysql://localhost:3306/";
		String driver = "org.apache.derby.jdbc.EmbeddedDriver";

		String str = "";

		try {
			if (ID_CONNECTION == DERBY_CONNECTION) {
				String protocol = "jdbc:derby:";
				driver = "org.apache.derby.jdbc.EmbeddedDriver";

				Class.forName(driver).newInstance();
				conng = DriverManager.getConnection(protocol + opens + ";create=false", user, pass);

				str = "create table " + table + " ( ";

				for (int i = 0; i < cname.length; i++) {
					str = str + cname[i] + " " + ctype[i];
					if (cindexes[i].equals("Y")) {
						str = str
								+ " not null primary key GENERATED BY DEFAULT AS IDENTITY (START WITH 1, INCREMENT BY 1)";
					}

					if (i < cname.length - 1) {
						str = str + ", ";
					} else {
						str = str + ")";
					}
				}

			} else if (ID_CONNECTION == MYSQL_CONNECTION) {
				dbURL = "jdbc:mysql://localhost:3306/";
				driver = "com.mysql.jdbc.Driver";

				Class.forName(driver).newInstance();

				String url = dbURL + dbName;
				conng = DriverManager.getConnection(url + "?user=" + user + "&password=" + pass);

				str = "create table " + table + " ( ";

				for (int i = 0; i < cname.length; i++) {
					str = str + cname[i] + " " + ctype[i];
					if (cindexes[i].equals("Y")) {
						str = str + " AUTO_INCREMENT PRIMARY KEY";
					}

					if (i < cname.length - 1) {
						str = str + ", ";
					} else {
						str = str + ")";
					}
				}

			} else if (ID_CONNECTION == POSTGRESQL_CONNECTION) {
				dbURL = "jdbc:postgresql://localhost:5432/";
				driver = "org.postgresql.Driver";

				Class.forName(driver).newInstance();

				String url = dbURL + dbName;

				conng = DriverManager.getConnection(url + "?user=" + user + "&password=" + pass);

				str = "create table " + table + " ( ";

				for (int i = 0; i < cname.length; i++) {

					if (cindexes[i].equals("Y")) {
						str = str + cname[i] + " SERIAL PRIMARY KEY";
					} else {
						str = str + cname[i] + " " + ctype[i];
					}

					if (i < cname.length - 1) {
						str = str + ", ";
					} else {
						str = str + ")";
					}
				}
			}

			conng.setAutoCommit(false);
			s = conng.createStatement();

			s.execute(str);

			conng.commit();

			if (s != null)
				s.close();

			if (conng != null) {
				conng.close();
			}

		} catch (InstantiationException e) {
			e.printStackTrace();
			return;
		} catch (IllegalAccessException e) {
			e.printStackTrace();
			return;
		} catch (ClassNotFoundException e) {
			e.printStackTrace();
			return;
		} catch (SQLException e) {
			e.printStackTrace();
			return;
		}

		System.out.println("Table created!");
	}

	/**
	 * Restart the primary key with a new int value.
	 * 
	 * @param dbName
	 *            database name
	 * @param table
	 *            table name
	 * @param user
	 *            user name, can be null
	 * @param pass
	 *            password, can be null
	 * @param pkColumn
	 *            primary key column name
	 * @param newInt
	 *            new starting integer
	 */
	public static void restartPK(String dbName, String table, String user, String pass, String pkColumn,
			int newInt) {

		Connection conng = null;
		Statement s = null;

		String datas = dataFolderUsedByDerbyEmbedded;// default: Data folder
		String currentDir = System.getProperty("user.dir");
		String file_sep = System.getProperty("file.separator");
		String opens = currentDir + file_sep + datas;

		opens = opens + file_sep + dbName;

		String dbURL = "jdbc:mysql://localhost:3306/";// initialization
		String driver = "org.apache.derby.jdbc.EmbeddedDriver";// initialization

		try {
			if (ID_CONNECTION == DERBY_CONNECTION) {
				String protocol = "jdbc:derby:";
				driver = "org.apache.derby.jdbc.EmbeddedDriver";

				Class.forName(driver).newInstance();
				conng = DriverManager.getConnection(protocol + opens + ";create=false", user, pass);

			} else if (ID_CONNECTION == MYSQL_CONNECTION) {
				dbURL = "jdbc:mysql://localhost:3306/";
				driver = "com.mysql.jdbc.Driver";

				Class.forName(driver).newInstance();
				String url = dbURL + dbName;
				conng = DriverManager.getConnection(url + "?user=" + user + "&password=" + pass);

			} else if (ID_CONNECTION == POSTGRESQL_CONNECTION) {
				dbURL = "jdbc:postgresql://localhost:5432/";
				driver = "org.postgresql.Driver";

				Class.forName(driver).newInstance();
				String url = dbURL + dbName;
				conng = DriverManager.getConnection(url + "?user=" + user + "&password=" + pass);

			}

			conng.setAutoCommit(false);
			s = conng.createStatement();

			String str = "alter table " + table + " alter " + pkColumn + " RESTART WITH " + newInt;
			s.execute(str);// "drop table " + table);

			conng.commit();

			if (s != null)
				s.close();

			if (conng != null) {
				conng.close();
			}

		} catch (InstantiationException e) {
			e.printStackTrace();
			return;
		} catch (IllegalAccessException e) {
			e.printStackTrace();
			return;
		} catch (ClassNotFoundException e) {
			e.printStackTrace();
			return;
		} catch (SQLException e) {
			e.printStackTrace();
			return;
		}

		System.out.println("PK restarted!");
	}
	/// ===========END UTILITIES============

	/**
	 * Setting up connection. For derby embedded the driver is
	 * org.apache.derby.jdbc.EmbeddedDriver and the protocol: jdbc:derby:. For
	 * MySQL, the driver is com.mysql.jdbc.Driver and the protocol:
	 * jdbc:mysql://localhost:3306/. For PostgreSQL, the driver is
	 * org.postgresql.Driver and the protocol:
	 * jdbc:postgresql://localhost:5432/. 
	 * 
	 * @param dbURL
	 *            database name to connect. For derby embedded, dbURL must be
	 *            the full path. For MySQL or PostgreSQL this is simply the
	 *            name, for example: testdb
	 * @param user
	 *            user name, can be null
	 * @param password
	 *            password, can be null
	 * @return the connection to database
	 */
	public static Connection getConnection(String dbURL, String user, String password) {
		// When an application accesses a Derby database using the Embedded
		// Derby JDBC driver, the Derby engine does not run in a
		// separate process, and there are no separate database processes to
		// start up and shut down. Instead, the Derby database engine
		// runs inside the same Java Virtual Machine (JVM) as the application.
		// So, Derby becomes part of the application just like any
		// other jar file that the application uses.
		// For derby embedded: only one JVM may boot ("open") that database, so
		// multiple applications running in different JVMs
		// cannot access the same database. There's one JVM per Java
		// application. There shouldn't be any connection between them
		// unless you establish one, e.g. with networking.=>run testClass class
		// fill fail because the previous tescClass is already
		// connected to derby!! So for allowing multiple instances=>use Derby
		// Network Server or MySQL or PostgreSQL DB!!!

		//if (con != null) {// is already connected; con is a static private
							// member!
		//	return;// prevent generating multiple connections to the same DB!!!
			//The reason: con is very expensive (resource wise) object. It is not efficient to have multiple con all over the place!
		//}
		Connection con = null;

		try {
			if (ID_CONNECTION == DERBY_CONNECTION)
				con = initDerby(dbURL, user, password);
			if (ID_CONNECTION == MYSQL_CONNECTION)
				con = initMysql(dbURL, user, password);
			if (ID_CONNECTION == POSTGRESQL_CONNECTION)
				con = initPostgresql(dbURL, user, password);
		} catch (InstantiationException e) {
			derbyClosed = true;
			e.printStackTrace();
		} catch (IllegalAccessException e) {
			derbyClosed = true;
			e.printStackTrace();
		} catch (ClassNotFoundException e) {
			derbyClosed = true;
			e.printStackTrace();
		} catch (SQLException e) {
			derbyClosed = true;
			e.printStackTrace();
		}
		
		return con;
	}

	/**
	 * Set-up the connection for embedded derby. Derby connector jar (derby.jar)
	 * must be placed in the application class-path.
	 * 
	 * @param dbURL
	 *            passed by constructor
	 * @param user
	 *            passed by constructor
	 * @param password
	 *            passed by constructor
	 * @throws InstantiationException
	 *             can throw this java exception
	 * @throws IllegalAccessException
	 *             can throw this java exception
	 * @throws ClassNotFoundException
	 *             can throw this java exception
	 * @throws SQLException
	 *             can throw this java exception
	 * @return the connection to database
	 */
	private static Connection initDerby(String dbURL, String user, String password)
			throws InstantiationException, IllegalAccessException, ClassNotFoundException, SQLException {
		if (derbyClosed) {
			String driver = "org.apache.derby.jdbc.EmbeddedDriver";
			// disable log file!
			System.setProperty("derby.stream.error.method", "danfulea.db.DatabaseAgent.disableDerbyLogFile");

			Class.forName(driver).newInstance();

			derbyClosed = false;
		}

		String protocol = "jdbc:derby:";
		// dbURL path must be set relative to the application folder. e.g.
		// mainApp/Data/database1
		Connection con = 
				DriverManager.getConnection(protocol + dbURL + ";create=false",
						user, password);
		return con;
	}

	/**
	 * Do not want to create derby log file. This must be public and static. If,
	 * for some reason, you want derby logs then comment-out this method as well
	 * as the following line inside initDerby method:
	 * System.setProperty("derby.stream.error.method",
	 * "danfulea.db.DatabaseAgent.disableDerbyLogFile");
	 * 
	 * @return a null OutputStream
	 */
	public static OutputStream disableDerbyLogFile() {
		return new OutputStream() {
			public void write(int b) throws IOException {
				// Ignore all log messages
			}
		};
	}

	/**
	 * Shuts down derby.
	 */
	public static void shutdownDerby() {

		try {
			// ------------------------this will close the connection first!
			//con.close();
			//con = null;// sets connection to null!!
			//connections are closed elsewhere!!!!!!!
			// ------------------------

			// the shutdown=true attribute shuts down Derby server
			DriverManager.getConnection("jdbc:derby:;shutdown=true");
			derbyClosed = true;

			// To shut down a specific database only, but keep the
			// engine running (for example for connecting to other
			// databases), specify a database in the connection URL:
			// DriverManager.getConnection("jdbc:derby:" + dbName +
			// ";shutdown=true");
		} catch (SQLException se) {
			if (((se.getErrorCode() == 50000) && ("XJ015".equals(se.getSQLState())))) {
				// we got the expected exception
				derbyClosed = true;// just in case
				// System.out.println("Derby shut down normally");
				// Note that for single database shutdown, the expected
				// SQL state is "08006", and the error code is 45000.
			} else {
				// if the error code or SQLState is different, we have
				// an unexpected exception (shutdown failed)
				derbyClosed = false;
				System.err.println("Derby did not shut down normally");
				se.printStackTrace();
			}
		}
	}

	//**
	// * Shuts down the connection. 
	// */
	//Connections are closed elsewhere!!!!!
	/*public void shutdown() {
		if (con == null) {
			// System.out.println("Connection already closed!");
			return;// connection already closed
		}
		if (ID_CONNECTION == DERBY_CONNECTION) {
			shutdownDerby();// this works because it is embedded
		} else {
			try {
				con.close();
				con = null;// mysql server must be stopped manually if needed!
			} catch (SQLException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
	}*/

	/**
	 * Set-up the connection for MySQL. MySQL connector jar must be placed in
	 * the application class-path.
	 * 
	 * @param dbName
	 *            passed by constructor
	 * @param user
	 *            passed by constructor
	 * @param password
	 *            passed by constructor
	 * @throws InstantiationException
	 *             can throw this java exception
	 * @throws IllegalAccessException
	 *             can throw this java exception
	 * @throws ClassNotFoundException
	 *             can throw this java exception
	 * @throws SQLException
	 *             can throw this java exception
	 * @return the connection to database
	 */
	private static Connection initMysql(String dbName, String user, String password)
			throws InstantiationException, IllegalAccessException, ClassNotFoundException, SQLException {
		// base url used mysql connector! e.g.
		// mysql-connector-java-3.1.0-alpha-bin.jar
		// jar file must be placed in classpath, e.g. mainApp/lib folder

		// String urlbase = "jdbc:mysql://localhost/";//for default single-port
		// open (which is mysql server port 3306) environment
		String urlbase = "jdbc:mysql://localhost:3306/";// we can have multiple
														// ports open so must be
														// explicit!
		String url = urlbase + dbName;

		Class.forName("com.mysql.jdbc.Driver").newInstance();

		Connection con = DriverManager.getConnection(url + "?user=" + user
				+ "&password=" + password);
		return con;
	}

	/**
	 * Set-up the connection for PostgreSQL. PostgreSQL connector jar must be
	 * placed in the application class-path.
	 * 
	 * @param dbName
	 *            passed by constructor
	 * @param user
	 *            passed by constructor
	 * @param password
	 *            passed by constructor
	 * @throws InstantiationException
	 *             can throw this java exception
	 * @throws IllegalAccessException
	 *             can throw this java exception
	 * @throws ClassNotFoundException
	 *             can throw this java exception
	 * @throws SQLException
	 *             can throw this java exception
	 * @return the connection to database    
	 *             
	 */
	private static Connection initPostgresql(String dbName, String user, String password)
			throws InstantiationException, IllegalAccessException, ClassNotFoundException, SQLException {
		// base url used postgresql connector! e.g.
		// postgresql-9.4.1208.jar
		// jar file must be placed in classpath, e.g. mainApp/lib folder

		String urlbase = "jdbc:postgresql://localhost:5432/";
		String url = urlbase + dbName;

		Class.forName("org.postgresql.Driver").newInstance();

		Connection con = DriverManager.getConnection(url + "?user=" +
		user + "&password=" + password);

		return con;
	}

	//**
	// * 
	// * @return the connection object
	// */
	/*public static Connection getConnection() {
		return con;
	}*/

	//**
	// * 
	// * @return the auto-incremented, primary key column name
	// */
	//NOT HERE, THIS IS THE BASE CLASS!!
	/*public String getPrimaryKey() {
		return primaryKey;
	}*/

	//**
	// * Set-up the auto-incremented, primary key column name
	// * 
	// * @param primaryKey
	// *            the auto-incremented, primary key column name
	// */
	//NOT HERE, THIS IS THE BASE CLASS!!
	/*public void setPrimaryKey(String primaryKey) {
		this.primaryKey = primaryKey;
	}*/

	/**
	 * This is the KEY method for all SELECT queries. It populates the required
	 * members, such as colNumeVec or dataVec, with data for displaying the
	 * result.
	 * @param con, the connection to database
	 * @param command
	 *            the command containing the SELECT query. Example: select *
	 *            from my_table.
	 * @throws SQLException
	 *             can throw this exception
	 */
	public static void select(Connection con, String command) throws SQLException {

		Vector<Object> randVec;
		colCount = 0;// init
		rowCount = 0;// init
		colNumeVec = new Vector<Object>();
		dataVec = new Vector<Object>();
		columnClass = new Vector<String>();
		columnType = new Vector<Integer>();

		PreparedStatement pstmt = con.prepareStatement(command, ResultSet.TYPE_SCROLL_SENSITIVE,
				ResultSet.CONCUR_UPDATABLE);
		// scrolable and updatable
		// Scroll sensitive specifies that a resultset is scrollable in either
		// direction (same as Insensitive)
		// and is affected by changes committed by other transactions or
		// statements within the same transaction.
		// CONCUR_UPDATABLE =>The ResultSet object can be updated using the
		// ResultSet interface. Eg. for result set to be able
		// to delete a row, rs.deleteRow() this must be uddatable!

		ResultSet rs = pstmt.executeQuery();
		ResultSetMetaData rsmd = rs.getMetaData();

		colCount = rsmd.getColumnCount();
		for (int i = 1; i <= colCount; i++) {
			colNumeVec.addElement(rsmd.getColumnName(i));
			columnClass.addElement(rsmd.getColumnClassName(i));
			columnType.addElement(rsmd.getColumnType(i));
			// System.out.println("table "+ rsmd.getTableName(i)+" ; "
			// +rsmd.getColumnClassName(i)+" ; "+rsmd.getColumnType(i));
		}

		while (rs.next()) {
			randVec = new Vector<Object>();
			for (int i = 1; i <= colCount; i++) {
				//the following OBJ is required for not having sql exception when dealing with 
				//LOB like columns....their object CAN BE RECEIVED ONLY ONCE!!!!! 
				Object obj=rs.getObject(i);
				if (obj == null)//if (rs.getObject(i) == null)
					randVec.addElement("");// reference must exists!
				else
					randVec.addElement(obj);//randVec.addElement(rs.getObject(i));
			} // end for

			dataVec.addElement(randVec);
			rowCount++;
		} // end while

		if (pstmt != null)
			pstmt.close();
		if (rs != null)
			rs.close();// rsmd as well!
	}

	/**
	 * 
	 * @return the number of rows after the SELECT query
	 */
	public static int getRowCount() {
		return rowCount;// rows;
	}

	/**
	 * 
	 * @return the number of columns after the SELECT query
	 */
	public static int getColumnCount() {
		return colCount;// cols;
	}

	/**
	 * 
	 * @param row
	 *            the row index
	 * @param col
	 *            the column index
	 * @return a specific value object stored at some row and some column
	 */
	public static Object getValueAt(int row, int col) {
		@SuppressWarnings("unchecked")
		Vector<Object> v = (Vector<Object>) dataVec.elementAt(row);
		return v.elementAt(col);
	}

	/**
	 * 
	 * @return the column names vector
	 */
	public static Vector<Object> getColumnNames() {
		return colNumeVec;
	}

	/**
	 * 
	 * @return the data vector
	 */
	public static Vector<Object> getData() {
		return dataVec;
	}

	/**
	 * 
	 * @return the column types as Java class vector
	 */
	public static Vector<String> getColumnClass() {
		return columnClass;
	}

	/**
	 * 
	 * @return the column types as SQL type vector
	 */
	public static Vector<Integer> getColumnType() {
		return columnType;
	}

	/**
	 * Set the PreparedStatement for insert/update in complete agreement with
	 * column type. For example, it knows when to use pstmt.setString or
	 * pstmt.setDouble.
	 * 
	 * @param pstmt
	 *            The PreparedStatement to be properly set for insert/update
	 *            operation.
	 * @param ival
	 *            the question mark placeholder index representing the pair columnName-columnType (starting with 1)
	 * @param val
	 *            the value to be set at that column
	 * @param tvalue
	 *            the column type. This is of java.sql.Types type.
	 * @throws SQLException
	 *             can throw this exception
	 */
	private static void setPreparedStatement(PreparedStatement pstmt, int ival, String val, int tvalue) throws SQLException {

		switch (tvalue) {
		case Types.INTEGER:
			pstmt.setInt(ival, Convertor.stringToInt(val));
			break;
		case Types.CLOB:
			StringReader reader = new StringReader(val);
			int length = val.length();
			pstmt.setCharacterStream(ival, reader, length);
			break;
		case Types.DOUBLE:
			pstmt.setDouble(ival, Convertor.stringToDouble(val));
			break;
		case Types.FLOAT:
			// pstmt.setFloat(ival, Convertor.stringToFloat(val));
			pstmt.setDouble(ival, Convertor.stringToDouble(val));
			break;
		case Types.VARCHAR:
			pstmt.setString(ival, val);
			break;
		// -----------------------------------------
		case Types.BIGINT:
			pstmt.setLong(ival, Convertor.stringToLong(val));
			break;
		case Types.TIMESTAMP:
			pstmt.setTimestamp(ival, Convertor.stringToTimestamp(val));
			break;
		case Types.DATE:
			pstmt.setDate(ival, Convertor.stringToDate(val));
			break;
		case Types.BLOB:
			ByteArrayInputStream barray = new ByteArrayInputStream(val.getBytes());// plattform
																					// default
																					// encoding
			length = val.length();
			pstmt.setBinaryStream(ival, barray, length);
			break;
		case Types.BINARY:
			pstmt.setBytes(ival, val.getBytes());
			break;
		case Types.BOOLEAN:
			pstmt.setBoolean(ival, Convertor.stringToBoolean(val));
			break;
		case Types.BIT:
			pstmt.setBoolean(ival, Convertor.stringToBoolean(val));
			break;
		case Types.REAL:// as float
			pstmt.setFloat(ival, Convertor.stringToFloat(val));
			break;
		case Types.DECIMAL:
			pstmt.setBigDecimal(ival, Convertor.stringToBigDecimal(val));
			break;
		case Types.NUMERIC:
			pstmt.setBigDecimal(ival, Convertor.stringToBigDecimal(val));
			break;
		case Types.TINYINT:
			pstmt.setShort(ival, Convertor.stringToShort(val));
			break;
		case Types.SMALLINT:
			pstmt.setShort(ival, Convertor.stringToShort(val));
			break;
		default:
			pstmt.setString(ival, val);
		}

	}

	/**
	 * Inserts a row in database table.
	 * 
	 * @param con, the connection to database
	 * @param tabel
	 *            the table where to perform this operation
	 * @param cvalue
	 *            the array of column names
	 * @param tvalue
	 *            the array of column types
	 * @param value
	 *            the array of values to be inserted
	 * @throws SQLException
	 *             can throw this exception
	 */
	public static void insert(Connection con, String tabel, String[] cvalue, Integer tvalue[], String[] value) throws SQLException {
		PreparedStatement pstmt = null;

		String s = "insert into " + tabel + " (";
		for (int i = 0; i < cvalue.length; i++) {
			s = s + cvalue[i];
			if (i < cvalue.length - 1) {
				s = s + ",";
			} else {
				s = s + ")";
			}
		}

		s = s + " values (";

		for (int i = 0; i < value.length; i++) {
			s = s + "?";// '" + value[i] + "'";
			// this syntax:"'" + value[i] + "'"- works only for strings, not for
			// integers
			if (i < value.length - 1) {
				s = s + ",";
			} else {
				s = s + ")";
			}
		}

		// for integer type, the sintax must be for example:
		// s="insert into testTable2 (ID,EMAIL,ADDRESS) values
		// (77,'desr@ss','afasd 6')";
		// instead of what is generated here (look at "77" data above and
		// below):
		// s="insert into testTable2 (ID,EMAIL,ADDRESS) values
		// ('77','desr@ss','afasd 6')";
		// Solution: use ?,?,...and after prepareStatement, fill with strings!
		pstmt = con.prepareStatement(s);

		// now fill with specific data (by types)
		for (int i = 1; i <= value.length; i++) {
			// pstmt.setString(i, value[i-1]);//works with derby and MySQL but
			// fails in PostrgeSQL. SOLUTION: check column types!
			setPreparedStatement(pstmt, i, value[i - 1], tvalue[i - 1]);
		}

		pstmt.executeUpdate();

		if (pstmt != null)
			pstmt.close();
	}

	//**
	// * Updates a row in database table. This is based on primary key column name
	// * as reference where update takes place. The reference is passed in WHERE
	// * clause, e.g. ID = 1 where ID is the reference column name link and 1 is
	// * its corresponding reference value (IDValue).
	// * 
	// * @param tabel
	// *            the table where to perform this operation
	// * @param cvalue
	// *            the array of column names where values will be updated
	// * @param tvalue
	// *            the array of column types where values will be updated
	// * @param value
	// *            the array of values to be updated
	// * @param IDValue
	// *            the reference value for corresponding primary key link
	 //* @throws SQLException
	// *             can throw this exception
	// */
	//NOT HERE, THIS IS THE BASE CLASS!!
	/*public void update(String tabel, String[] cvalue, Integer tvalue[], String[] value, String IDValue)
			throws SQLException {
		update(tabel, cvalue, tvalue, value, primaryKey, IDValue);
	}*/

	/**
	 * Updates a row in database table. This is based on a column name as
	 * reference where update takes place. The reference is passed in WHERE
	 * clause, e.g. ID = 1 where ID is the the reference column name link
	 * (IDLink) and 1 is its corresponding reference value (IDValue).
	 * 
	 * @param con, the connection to database
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
	 *            the reference value for corresponding link
	 * @throws SQLException
	 *             can throw this exception
	 */
	public static void update(Connection con, String tabel, String[] cvalue, Integer tvalue[], String[] value, String IDLink, String IDValue)
			throws SQLException {
		PreparedStatement pstmt = null;

		String s = "update " + tabel + " set ";
		for (int i = 0; i < cvalue.length; i++) {
			s = s + cvalue[i] + "=" + "?";// '" + value[i] + "'";
			// The syntax: "'" + value[i] + "'"; works only for strings, not for
			// integers
			if (i < cvalue.length - 1) {
				s = s + ",";
			}
		}

		s = s + " where " + IDLink + "=" + IDValue;

		pstmt = con.prepareStatement(s);

		for (int i = 1; i <= value.length; i++) {
			// pstmt.setString(i, value[i-1]);//works with derby and MySQL but
			// fails in PostrgeSQL
			setPreparedStatement(pstmt, i, value[i - 1], tvalue[i - 1]);
		}

		pstmt.executeUpdate();

		if (pstmt != null)
			pstmt.close();
	}

	//**
	// * Deletes a row from database table. This is based on primary key column
	// * name as reference where delete takes place. The reference is passed in
	// * WHERE clause, e.g. ID = 1 where ID is the reference column name link and
	// * 1 is its corresponding reference value (IDValue).
	// * 
	// * @param tabel
	// *            the table where to perform this operation
	// * @param IDValue
	// *            the reference value for corresponding primary key link
	// * @throws SQLException
	// *             can throw this exception
	// */
	//NOT HERE, THIS IS THE BASE CLASS!!
	/*
	public void delete(String tabel, String IDValue) throws SQLException {
		delete(tabel, primaryKey, IDValue);
	}
	*/
	/**
	 * Deletes a row from database table. This is based on a column name as
	 * reference where delete takes place. The reference is passed in WHERE
	 * clause, e.g. ID = 1 where ID is the reference column name link (IDLink)
	 * and 1 is its corresponding reference value (IDValue).
	 * 
	 * @param con, the connection to database
	 * @param tabel
	 *            the table where to perform this operation
	 * @param IDLink
	 *            the reference column name used as link
	 * @param IDValue
	 *            the reference value for corresponding link
	 * @throws SQLException
	 *             can throw this exception
	 */
	public static void delete(Connection con, String tabel, String IDLink, String IDValue) throws SQLException {
		String str = "DELETE FROM " + tabel + " WHERE " + IDLink + "=" + IDValue;
		Statement s = con.createStatement(ResultSet.TYPE_FORWARD_ONLY, ResultSet.CONCUR_UPDATABLE);
		s.executeUpdate(str);

		if (s != null)
			s.close();

		// another delete method but more time-consuming is:
		// ResultSet res = s.executeQuery("SELECT * FROM " + maindbTable+ "
		// WHERE ID="+ID);
		// while (res.next()) {
		// res.deleteRow();//only 1 record if ID=primary key
		// }
		// if (res != null)
		// res.close();
	}
	
	/**
	 * Deletes multiple rows from database table. This is based on a column name as
	 * reference where delete takes place. The reference is passed in WHERE
	 * clause, e.g. ID = 1 where ID is the reference column name link (IDLink)
	 * and 1 is its corresponding reference value (IDValue).
	 * 
	 * @param con, the connection to database
	 * @param tabel
	 *            the table where to perform this operation
	 * @param IDLink
	 *            the reference column name array used as link
	 * @param IDValue
	 *            the reference value array for corresponding link
	 * @throws SQLException
	 *             can throw this exception
	 */
	public static void delete(Connection con, String tabel, String[] IDLink, String[] IDValue) throws SQLException {
		int n = IDLink.length;
		Statement s = null;
		
		for (int i = 0; i<n; i++){
			String str = "DELETE FROM " + tabel + " WHERE " + IDLink[i] + "=" + IDValue[i];
			s = con.createStatement(ResultSet.TYPE_FORWARD_ONLY, ResultSet.CONCUR_UPDATABLE);
			s.executeUpdate(str);
		}
		

		if (s != null)
			s.close();
	}

	/**
	 * Deletes all rows from database table
	 * 
	 * @param con, the connection to database
	 * @param tabel
	 *            the table where to perform this operation
	 * @throws SQLException
	 *             can throw this exception
	 */
	public static void deleteAll(Connection con, String tabel) throws SQLException {
		String str = "DELETE FROM " + tabel;
		Statement s = con.createStatement(ResultSet.TYPE_FORWARD_ONLY, ResultSet.CONCUR_UPDATABLE);
		s.executeUpdate(str);

		if (s != null)
			s.close();
	}
}
