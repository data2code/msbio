#!/usr/bin/env python
from __future__ import print_function
from __future__ import absolute_import
import pandas as pd
import pandas.io.sql
import re
import util
import db
import os
import six
from six.moves import range
import setting
import time
import json
#str(db.__class__).find("Oracle")>=0 # Oracle connection
#str(db.__class__).find('MySQL')>=0

def no_warning():
    # https://stackoverflow.com/questions/17626694/remove-python-userwarning
    import warnings
    warnings.simplefilter("ignore")

def db_type(con):
    if str(con.__class__).find('MySQL')>=0:
        return 'MYSQL'
    elif str(con.__class__).find('Oracle')>=0:
        return 'ORACLE'
    elif str(con.__class__).find('pgdb')>=0:
        return 'POSTGRES'
    elif str(con.__class__).find('sqlite')>=0:
        return 'SQLITE'

def from_sql(con, sql, params=(), index_col=None, coerce_float=True, l_select=None, l_commit=True):
    """
    Same as Pandas.io.sql.read_frame, except it takes params for binding

    Returns a DataFrame corresponding to the result set of the query string.

    Optionally provide an index_col parameter to use one of the
    columns as the index. Otherwise will be 0 to len(results) - 1.

    Parameters
    ----------
    sql: string
        SQL query to be executed
    con: DB connection object, optional
    index_col: string, optional column name to use for the returned DataFrame object.

    By default, the method returns None if SQL is not a SELECT statement, i.e., the
    keyword SELECT is not the first word in the SQL.
    However, this logic maybe wrong for complex SQL, such as
        with mytable (select * from t) select * from mytable
    Set l_select=True to return the dataframe in such case
    Set l_select=False if you explicitly do not want a dataframe

    Example:
    import MySQLdb as mysql
    con=mysql.connect("LDDB","usrname","password","LDDB")
    df=from_sql("select * from Assay where id = %s and AssayName like %s", con, (3214, 'Lung%'))
    df=from_sql("select * from Assay where id = %(id)s and AssayName like %(name)s", con, {'id':3214, 'name':'Lung%'})

    import cx_Oracle as oracle
    con=oracle.connect("usrname/password@imcprod")
    df=from_sql("select * from LDDB_WAREHOUSE.Assay where Assay_id = :1 and Assay_Name like :2", con, (3214, 'Lung%'))
    df=from_sql("select * from LDDB_WAREHOUSEAssay where Assay_id = :id and Assay_Name like :name", con, {'id':3214, 'name':'Lung%'})
    """
    cur = con.cursor()
    # we assum sql statement is always using '?' style binding
    # need to convert that into Oracle style if it's Oracle
    pattern=re.compile(r'\?')
    db_src=db_type(con)
    if type(params) not in (list, tuple):
        # if user pass in (1) instead of (1,)
        params=[params]
    if db_src=='MYSQL':
        sql=re.sub(r'%', '%%', sql)
        # in case one use like 'E%', % needs to be escape, but one should really use like ?
        sql=pattern.sub('%s', sql)
    elif db_src=='ORACLE':
        cur.arraysize=1024*16
        i=1
        while pattern.search(sql):
            sql=pattern.sub(':%d' % i, sql, count=1)
            i+=1
        # reserve space for clob parameters if there are any, otherwise error will occur.
        import cx_Oracle
        if util.is_python3():
            inputsizes = [cx_Oracle.CLOB if (type(p) in (str, bytes)) and len(p) > 4000 else None for p in params]
        else:
            inputsizes = [cx_Oracle.CLOB if isinstance(p, basestring) and len(p) > 4000 else None for p in params]
        cur.setinputsizes(*inputsizes)
    elif db_src=='POSTGRES':
        sql=pattern.sub('%s', sql)
    cur.execute(sql, params)
    if l_select is None:
        l_select = re.search(r'^\s*select\s', sql, flags=re.I) is not None
        if not l_select and db_src=='MYSQL' and re.search(r'^\s*show\s', sql, flags=re.I) is not None:
            l_select=True
    if l_select:
        rows = cur.fetchall() #pandas.io.sql._safe_fetch(cur)
        rows = [x for x in rows ] #convert tuple to list
        columns = [col_desc[0] for col_desc in cur.description]
    cur.close()
    if l_commit: con.commit()

    if not l_select:
        return None

    result = pandas.DataFrame.from_records(rows, columns=columns, coerce_float=coerce_float)
    if index_col is not None:
        result = result.set_index(index_col)

    return result

def table_exists(con, s_table, s_db=""):
    db_src=db_type(con)
    if db_src=='MYSQL':
        #s_sql="describe "+(f"{s_db}.{s_table}" if s_db else s_table)
        #print(s_sql)
        #t=from_sql(con, s_sql, l_select=True)
        #return len(t)>0
        # abandon below, as user may not have permission to read information_schema
        if "." in s_table:
            s_db, s_table=s_table.split(".")
        s_sql="select count(*) from information_schema.tables where table_name=?"
        param=[s_table]
        if s_db:
            s_sql+=" and table_schema=?"
            param.append(s_db)
        t=from_sql(con, s_sql, param)
    elif db_src=='ORACLE':
        t=from_sql(con, "select count(*) from user_tables where table_name=upper(?)", [s_table])
    elif db_src=='SQLITE':
        t=from_sql(con, "select count(*) from sqlite_master where type='table' and name=?", [s_table])
    else:
        util.error_msg('Unsupported database engine!')
    return t.iloc[0,0]>0

def get_con(s_name, auth_path=None, db=None):
    """For MySQL, one can specify a default database to connect to"""
    con=None
    one=get_con_info(s_name, auth_path=auth_path)
    s_db=db or one['DB']
    if one['TYPE']=='MYSQL':
        import MySQLdb as mysql
        #print(one['HOST'], one['USR'], one['PWD'], s_db, one['PORT'])
        if not one['PORT'] or pd.isnull(one['PORT']):
            con=mysql.connect(one['HOST'], one['USR'], one['PWD'], s_db, port=3306, charset='utf8')
        else:
            con=mysql.connect(one['HOST'], one['USR'], one['PWD'], s_db, port=int(one['PORT']), charset='utf8')
    elif one['TYPE']=='POSTGRES':
        import pgdb
        # make sure you do:
        #module load postgresql/12.2
        #export LD_LIBRARY_PATH=.:/apps/GNU/postgresql/12.2/lib
        #print(one['CONNECT'])
        #con=pgdb.connect(one['CONNECT'])
        con=pgdb.Connection(one['CONNECT'])
    elif one['TYPE']=='ORACLE':
        import cx_Oracle as oracle
        con=oracle.connect(one['CONNECT'])
    else:
        util.error_msg('Unsupported database engine: %s' % one['TYPE'])
    return con

def open_sqlite(s_db):
    import sqlite3
    con=sqlite3.connect(s_db)
    return con

def sql_in(s_sql_left, s_sql_right, S_id, num=1000, con=None, params_before=None, params_after=None, l_commit=True):
    if con is None:
        util.error_msg('Database connection is not provided!')
    S_id=list(S_id)
    n=len(S_id)
    # in case it's multi-line SQL statement
    s_sql_left=re.sub('[\r\n]', ' ', s_sql_left)
    s_sql_right=re.sub('[\r\n]', ' ', s_sql_right)
    s_sql_right=re.sub(r'^\s*\)', '', s_sql_right)
    pat=re.compile(r'\s+(?P<col>[\w.]+)(?P<not>\s+NOT)?\s+IN\s*\(\s*$', re.I)
    m=re.search(pat, s_sql_left)
    if m is None:
        util.error_msg('Left SQL does not ends with IN statement: %s' % s_sql_left)
    s_sql_left=s_sql_left[:m.start()]+" "
    s_col=m.group('col')
    l_NOT=m.group('not') is not None
    # somethings SQL contains its own parameters, we need to provide parameter before/after the
    # S_id parameters
    params_before = params_before or []
    params_after=params_after or []
    # If oracle, either use ToTable() or run multiple SQL queries
    # If multiple SQL are run, results do not support GROUP BY, ORDER BY, DISTINCT
    # as they are applied to individual SQL runs
    # For other SQL servers, results will be exact, via multiple OR statements
    if db_type(con)=='ORACLE': # and len(S_id)>num:
        # Custom function used where TOTABLE is defined under MY_SCHEMA
        MY_SCHEMA=setting.db['MY_SCHEMA']

        ### The definition of user function
        #CREATE OR REPLACE FUNCTION MY_SCHEMA.ToTable (
        #   P_STR     IN CLOB,
        #   P_DELIM   IN VARCHAR2 DEFAULT ',' ,
        #   P_LIKE    IN INT DEFAULT 0
        #)
        #   RETURN t_NTtype
        #AS
        #  L_DATA t_NTtype := t_NTtype ();
        #  L_STR  CLOB := P_STR || P_DELIM;
        #  L_SUBSTR VARCHAR2(4000);
        #  L_STEP PLS_INTEGER := 0;
        #  L_THIS INT := 1;
        #  L_PREV INT := 0;
        #  L_END CHAR := CASE P_LIKE WHEN 0 THEN NULL ELSE '%' END;
        #BEGIN
        #  LOOP
        #    L_STEP := L_STEP + 1;
        #    L_THIS := DBMS_LOB.INSTR(L_STR, P_DELIM, L_PREV + 1, 1);
        #    EXIT WHEN L_THIS = 0;
        #    L_SUBSTR :=
        #    TRIM(
        #      DBMS_LOB.SUBSTR(
        #        L_STR,
        #        L_THIS - L_PREV - 1,
        #        L_PREV + 1
        #      )
        #    );
        #    L_PREV := L_THIS;
        #    L_DATA.EXTEND();
        #    L_DATA(L_STEP) := L_SUBSTR || L_END;
        #  END LOOP;
        #  RETURN L_DATA;
        #END;
        ###

        t=from_sql(con, "select * from SYS.User_Objects where object_name='TOTABLE'")
        if len(t):
            S_id = [str(id) for id in S_id]
            S_id_filtered = [id for id in S_id if len(id) > 20]  # 20 is entry length limit within ToTable()
            S_id = [','.join([id for id in S_id if len(id) <= 20])]
            if len(S_id_filtered) > 0:
                util.warn_msg('The following parameters were filtered out from the sql due to the length constraint: ' + str(S_id_filtered))
            if l_NOT:
                s_id = ('(NOT exists (select /*+ CARDINALITY(tmp_table, %d) */ 1 from Table('+MY_SCHEMA+'.ToTable(?)) tmp_table where tmp_table.COLUMN_VALUE = ' + s_col + '))') % len(S_id)
            else:
                s_id = ('(exists (select /*+ CARDINALITY(tmp_table, %d) */ 1 from Table('+MY_SCHEMA+'.ToTable(?)) tmp_table where tmp_table.COLUMN_VALUE = ' + s_col + '))') % len(S_id)
            return from_sql(con, s_sql_left+s_id+s_sql_right, params=params_before+S_id+params_after, l_commit=l_commit)
        #else:
        #    t=[]
        #    for i in xrange(0, n, num):
        #        j=min(n,i+num)
        #        s_id=",".join(["?"]*(j-i))
        #        t2=from_sql(con, s_sql_left+" "+s_col+" IN ("+s_id+") "+s_sql_right, params=params_before+S_id[i:j]+params_after)
        #        t.append(t2)
        #    return pd.concat(t, axis=0, ignore_index=True)
    S=[]
    for i in range(0, n, num):
        j=min(n,i+num)
        S.append(",".join(["?"]*(j-i)))
    if l_NOT:
        s_id="("+s_col+" NOT IN ("+(") AND "+s_col+" NOT IN (").join(S)+"))"
    else:
        s_id="("+s_col+" IN ("+(") OR "+s_col+" IN (").join(S)+"))"
    #if db_type(con) =='MYSQL' and not S:
    if len(S)==0:
        if l_NOT:
            s_id = "("+s_col+" NOT IN (''))"
        else:
            s_id = "("+s_col+" IN (''))"
    t=from_sql(con, s_sql_left+s_id+s_sql_right, params=params_before+S_id+params_after, l_commit=l_commit)
    return t

def call_proc(con, s_proc, params=()):
    """http://www.oracle.com/technetwork/articles/prez-stored-proc-084100.html"""
    cur = con.cursor()
    cur.callproc(s_proc, params)

def sql_in_old2(s_sql_left, s_sql_right, S_id, num=1000, con=None, params_before=None, params_after=None, l_commit=True):
    if con is None:
        util.error_msg('Database connection is not provided!')
    S_id=list(S_id)
    n=len(S_id)
    t=[]
    # in case it's multi-line SQL statement
    s_sql_left=re.sub('[\r\n]', ' ', s_sql_left)
    s_sql_right=re.sub('[\r\n]', ' ', s_sql_right)
    s_sql_right=re.sub(r'^\s*\)', '', s_sql_right)
    pat=re.compile(r'\s+([\w.]+)\s+IN\s*\(\s*$', re.I)
    m=re.search(pat, s_sql_left)
    if m is None:
        util.error_msg('Left SQL does not ends with IN statement: %s' % s_sql_left)
    s_sql_left=s_sql_left[:m.start()]+" "

    s_col=m.groups()[0]
    # somethings SQL contains its own parameters, we need to provide parameter before/after the
    # S_id parameters
    params_before = params_before or []
    params_after=params_after or []
    S=[]
    for i in range(0, n, num):
        j=min(n,i+num)
        S.append(",".join(["?"]*(j-i)))
    s_id="("+s_col+" IN ("+(") OR "+s_col+" IN (").join(S)+"))"
    #s_id = "(" + " OR ".join(["("+s_col + " IN (" + x +"))" for x in S]) + ")"
    #print s_sql_left+s_id+s_sql_right
    t=from_sql(con, s_sql_left+s_id+s_sql_right, params=params_before+S_id+params_after, l_commit=l_commit)
    return t

def sql_in_old1(s_sql_left, s_sql_right, S_id, num=1000, con=None, params_before=None, params_after=None):
    if con is None:
        util.error_msg('Database connection is not provided!')
    S_id=list(S_id)
    n=len(S_id)
    t=[]
    # somethings SQL contains its own parameters, we need to provide parameter before/after the
    # S_id parameters
    params_before = params_before or []
    params_after=params_after or []
    for i in range(0, n, num):
        j=min(n,i+num)
        s_id=",".join(["?"]*(j-i))
        t2=from_sql(con, s_sql_left+s_id+s_sql_right, params=params_before+S_id[i:j]+params_after)
        #if t2 is not None and len(t2)>0:
        t.append(t2)
    t=pd.concat(t, axis=0, ignore_index=True)
    return t

def create_table(con, s_table, tbl, S_type=None, c_rename = None):
    """Creates an empty table in sqlite database. The number of columns in a new table, their types and names
    will match the given pandas dataframe object.
    tbl     - pandas dataframe
    s_table - name of the table to be created"""
    cur = con.cursor()
    from string import join
    if S_type is None:
        S_type=util.col_types(tbl)
    fields = []
    for i,c_name in enumerate(tbl.columns):
        if S_type[i] =="i":
            f_type= "integer"
        elif S_type[i] =="r":
            f_type = "real"
        else:
            f_type = "text"
        if c_rename is not None:
            if c_name in c_rename: c_name = c_rename[c_name]
        fields.append(join([c_name,f_type]))
    query = "create table " + s_table + " (" + join(fields,",") + ")"
    cur.execute(query)
    con.commit()

def load_table(con, s_table, tbl, S_type=None):
    """Adds records from the existing pandas dataframe into a table in sqlite database.
    tbl     - pandas dataframe
    s_table - name of the table where new records will be added"""
    if S_type is None:
        S_type=util.col_types(tbl)
    ptype = {'i':int, 'r':float, 's':object}
    num_param = ','.join(['?'] * len(S_type))
    for i in range(len(S_type)):
        tbl.iloc[:,i] = tbl.iloc[:,i].astype(ptype[S_type[i]])
    param = [tuple(row) for row in tbl.values]
    cur = con.cursor()
    cur.executemany('insert into '+s_table+' values (' + num_param + ')', param)
    con.commit()

def get_con_info(s_name, auth_path=None):
    constr = os.environ.get(f'GNF_DB_{s_name.upper()}', None)
    if constr is not None:
        return json.loads(constr)
    auth_path = os.path.dirname(os.path.abspath(__file__)) + '/db.csv' if auth_path is None else auth_path
    csv_path=util.format_path(auth_path)
    if os.path.exists(csv_path):
        t_db=pd.read_csv(csv_path, dtype={'ID': object})
    elif os.path.exists(csv_path+".crypt"):
        from rsacipher import RSACipher
        s_dir=os.environ['HOME']+"/.ssh/"
        cipher=RSACipher(s_dir+"cheminfo_rsa.pem")
        data=cipher.read_encrypt_file(csv_path+".crypt")
        from io import BytesIO
        t_db=pd.read_csv(BytesIO(data))
    else:
        util.error_msg(f"Missing {auth_path}")
    t_db.fillna('', inplace=True)
    t_db=t_db[t_db['ID']==s_name]
    if len(t_db)==0:
        util.error_msg('Database %s is not defined!' % s_name)
        return None
    return t_db.iloc[0]

def get_curr_db(con):
    s=get_type(con)
    if s=='MYSQL':
        t=from_sql(con, 'select database()')
    elif s=='ORACLE':
        t=from_sql(con, "select sys_context('userenv','current_schema') from dual")
    else:
        return ""
    return t.iloc[0,0]

def get_databases(con):
    return from_sql(con, 'show databases', l_select=True).iloc[:,0].tolist()

def get_tables(con, s_db=None):
    if s_db is None:
        s_db=get_curr_db(con)
    s=get_type(con)
    if s=='MYSQL':
        return from_sql(con, 'show tables from %s' % s_db, l_select=True).iloc[:,0].tolist()
    elif s=='ORACLE':
        return from_sql(con, "SELECT DISTINCT OBJECT_NAME FROM ALL_OBJECTS WHERE OBJECT_TYPE = 'TABLE' AND OWNER = ?", params=(s_db))
    return None

def use_database(con, s_db):
    from_sql(con, f"use {s_db}")

def database_exists(con, s_db):
    t=from_sql(con, "SELECT SCHEMA_NAME FROM INFORMATION_SCHEMA.SCHEMATA WHERE SCHEMA_NAME = ?", params=[s_db])
    return len(t)>0

class DB:
    """This is a OO wrapper class"""

    def __init__(self, name="", auth_path=None, sql_lite=None, db=None):
        """If sql_lite is not None, it should contain the connection string, db and auth_path will be ignored"""
        self._name=name
        self._con=None
        self._sql_lite = None
        self._auth_path = None
        if sql_lite is not None:
            self._sql_lite = sql_lite
            one={"ID":name, "DB":db, "COMMENT": auth_path}
        else:
            one=get_con_info(name, auth_path=self._auth_path)
            self._auth_path= os.path.dirname(os.path.abspath(__file__)) + '/db.csv' if auth_path is None else auth_path
        self._db=db or one['DB']
        self._auto_commit=True

    def con(self):
        if self._con is None:
            no_warning()
            if self._sql_lite is None:
                self._con=get_con(self._name, auth_path=self._auth_path, db=self._db)
            else:
                import sqlite3
                self._con=sqlite3.connect(self._sql_lite)
        return self._con

    def check_con(self):
        if self._con is not None:
            try:
                self.from_sql('select 1 from dual')
            except Exception as e:
                self._con=None
        return self._con is not None

    def __enter__(self):
        return self

    def __exit__(self,exc_type, exc_value, traceback):
        if self._con is not None:
            self._con.close()

    def from_sql(self, sql, params=(), index_col=None, coerce_float=True, l_select=None):
        db_src=db_type(self.con())
        t=None
        if db_src=='MYSQL':
            n_try=3
            i_try=1
            while True:
                try:
                    return from_sql(self.con(), sql, params=params, index_col=index_col, coerce_float=coerce_float, l_select=l_select, l_commit=self._auto_commit)
                except Exception as e:
                    if ('connect' in e.args[1].lower() or int(e.args[0]) in (2006,)) and i_try<=n_try:
                        #2006, 'MySQL server has gone away'
                         # connection glitch, retry
                        print(e)
                        time.sleep(30)
                        i_try+=1
                        self._con=None
                    else:
                        raise e
                        break
        else:
            return from_sql(self.con(), sql, params=params, index_col=index_col, coerce_float=coerce_float, l_select=l_select, l_commit=self._auto_commit)
        return None

    def sql_in(self, s_sql_left, s_sql_right, S_id, num=1000, params_before=None, params_after=None):
        return sql_in(s_sql_left, s_sql_right, S_id, num=num, con=self.con(), params_before=params_before, params_after=params_after, l_commit=self._auto_commit)

    def load_table(self, s_table, tbl, S_type=None):
        return load_table(self.con(), s_table, tbl, S_type=S_type)

    def get_con_info(self):
        return get_con_info(self._name, auth_path=self._auth_path)

    def create_table(self, s_table, tbl, S_type=None, c_rename = None):
        return create_table(self.con(), s_table, tbl, S_type=S_typle, c_rename = c_rename)

    def call_proc(self, s_proc, params=()):
        return call_proc(self.con(), s_proc, params=params)

    def table_exists(self, s_table, s_db=None):
        if s_db is None: s_db=self._db
        return table_exists(self.con(), s_table, s_db=self._db)

    def database_exists(self, s_db):
        return database_exists(self.con(), s_db)

    def db_type(self):
        return db_type(self.con())

    def dump_databases(self, s_file, S_db, s_host=None, l_include_db_create=True, no_data=False):
        one=self.get_con_info()
        if s_host is not None: one['HOST']=s_host # sometimes we want to provide os.environ['HOSTNAME'] in order to work with local MySQL
        s_create="--databases" if l_include_db_create else "--no-create-db"
        if no_data: s_create+=" --no-data" # schema only
        if type(S_db) in (list, tuple):
            s_db=" ".join(S_db)
        else:
            s_db=str(S_db)
        if s_file.endswith('.gz'):
            s_cmd="""mysqldump -h %s -u %s -p'%s' %s %s |gzip > %s""" % (one['HOST'], one['USR'], one['PWD'], s_create, s_db, s_file)
        else:
            s_cmd="""mysqldump -h %s -u %s -p'%s' --result-file=%s %s %s""" % (one['HOST'], one['USR'], one['PWD'], s_file, s_create, s_db)
        #print(s_cmd)
        util.unix(s_cmd, l_print=True, l_error=True)

    def dump_tables(self, s_file, s_db, S_table, s_host=None, no_data=False):
        one=self.get_con_info()
        if s_host is not None: one['HOST']=s_host
        if type(S_table) in (list, tuple):
            s_table=" ".join(S_table)
        else:
            s_table=str(S_table)
        s_no_data=" --no-data" if no_data else ""
        if s_file.endswith('.gz'):
            s_cmd="""mysqldump -h %s -u %s -p'%s' %s %s %s %s | gzip > %s""" % (one['HOST'], one['USR'], one['PWD'], s_no_data, "", s_db, s_table, s_file)
        else:
            s_cmd="""mysqldump -h %s -u %s -p'%s' %s --result-file=%s %s %s %s""" % (one['HOST'], one['USR'], one['PWD'], s_no_data, s_file, "", s_db, s_table)
        #print(s_cmd)
        util.unix(s_cmd, l_print=True, l_error=True)

    def load_tables_file(self, s_file, s_target_db, s_host=None):
        one=self.get_con_info()
        if s_host is not None: one['HOST']=s_host
        if s_file.endswith('.gz'):
            s_cmd="""gunzip < %s | /usr/bin/mysql -h %s -u %s -p'%s' %s""" % (s_file, one['HOST'], one['USR'], one['PWD'], s_target_db)
        else:
            s_cmd="""/usr/bin/mysql -h %s -u %s -p'%s' %s < %s""" % (one['HOST'], one['USR'], one['PWD'], s_target_db, s_file)
        #print(s_cmd)
        util.unix(s_cmd, l_print=True, l_error=True)

    def get_curr_db(self):
        return get_curr_db(self.con())

    def get_databases(self):
        return get_databases(self.con())

    def get_tables(self, s_db=None):
        return get_tables(self.con(), s_db)

    def use_database(self, s_db):
        use_database(self.con(), s_db)
        self._db=s_db

if __name__ == '__main__':
    con=get_con('TEST_ORACLE')
    t= from_sql(con, "select dept_id from co.depts")
    # convert np.int64 to int, otherwise, oracle driver does not know how to handle
    S_id=[int(x) for x in t.dept_id[:10]]
    t=sql_in("select * from co.dept where dept_id in (", ")", S_id, con=con)
    #print t
    #exit()
    df= from_sql(con, "select * from co.dept where dept_id = ? and dept_Name like ?", params=[ 3214, 'Sales%' ])
    print("Testing Oracle...")
    print(df[:3])
    #exit()
    ##t=from_sql(con, 'select * from test.test_table where id=?', params=[1])
    #from_sql(con, "insert into test.test_table (id,val,Created,LastUpdated) values (null, ?, now(), now())", ['insert'])
    #from_sql(con, "delete from test.test_table where id>2")
    #t=from_sql(con, 'select * from test.test_table')
    print("testing MySQL ...")
    #print t[:3]
    con=get_con('TEST_POSTGRES')
    print (con)
    print ("Testing Postgres")
    t=from_sql(con, 'select * from co.depts where dept_id>? limit 5 offset 0', params=[100])
    print (t[:3])

