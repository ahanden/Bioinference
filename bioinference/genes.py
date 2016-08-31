#!/bin/python

import MySQLdb

"""
These methods are meant for the querying of the genes database.
"""

gene_conn = None # Static database connection (to save typing)

def connect(*args):
    """Stores a MySQLdb connection statically within the class. This is
    a shortcut for users so they don't have to constantly pass connection/
    cursor objects to methods over and over.

    Paramters
    ---------
    conn : a MySQLdb connection

    OR

    host : the host location of the database

    db : the name of the database
    
    user : the username for database access
    
    passwd : the password associated with the use
    """

    global gene_conn

    if len(args) == 1:
        if type(args[0]) is MySQLdb.conn:
            gene_conn = args[0]
        else:
            raise TypeError("connect() takes a MySQLdb.conn object as its first argument")

    elif len(args) == 4:
        gene_conn = MySQLdb.connect(
            host=args[0],
            db=args[1],
            user=args[2],
            passwd=args[3])

    else:
        raise TypeError("connect() takes 1 or 4 arguments (%d given)" % len(args))


def checkEID(eid, conn=None):
    """Returns valid Entrez IDs given an uncertain Entrez ID.

    Entrez IDs may, from time to time, be discontinued for many
    reasons. This method interfaces with the gene database to 
    check that a given Entrez ID is still valid, or to find what
    the current valid Entrez ID is for an old, discontinued ID.

    Paramters
    ---------
    eid : An Entrez ID
    
    conn : A MySQL connection (optional if you already called connect())

    Returns
    -------
    valid_eids : A set of valid Entrez IDs. The list will be
                 empty if no valid IDs were found.
    """

    if conn is None:
        if gene_conn is None:
            raise TypeError("checkEID() takes 2 arguments if connect() has not been called")
        conn = gene_conn
    cursor = conn.cursor()

    valid_eids = set()
    cursor.execute("SELECT EXISTS(SELECT * FROM genes WHERE entrez_id = %(eid)s)", {'eid': eid})
    if cursor.fetchone()[0] == 1:
        valid_eids.add(int(eid))
    else:
        cursor.execute("SELECT entrez_id FROM discontinued_genes WHERE discontinued_id = %(eid)s", {'eid': eid})
        valid_eids.update([int(row[0]) for row in cursor.fetchall()])
   
    return valid_eids

def getEID(symbol, conn=None):
    """Returns a valid Entrez ID given a gene symbol.

    Paramters
    ---------
    symbol : a gene symbol to convert
    
    conn : A MySQL connection (optional if you already called connect())

    Returns
    -------
    eids : a set of Entrez IDs - which will be empty if no valid IDs were found.
    """

    if conn is None:
        if gene_conn is None:
            raise TypeError("getEID() takes 2 arguments if connect() has not been called")
        conn = gene_conn
    cursor = conn.cursor()

    # Start with a direct query
    args = {"symbol": symbol}
    cursor.execute("SELECT entrez_id FROM genes WHERE symbol = %(symbol)s", args)
    eids = [int(row[0]) for row in cursor.fetchall()]

    # If that didn't work, search for a discontinued symbol
    if not eids:
        cursor.execute("SELECT entrez_id FROM discontinued_genes WHERE discontinued_symbol = %(symbol)s", args)
        eids = [int(row[0]) for row in cursor.fetchall()]

        # if THAT didn't work, search for synonyms
        if not eids:
            cursor.execute("SELECT entrez_id FROM gene_synonyms WHERE symbol = %(symbol)s", args)
            eids = [int(row[0]) for row in cursor.fetchall()]

    # Return whatever we found (if anything)
    return set(eids)

def getSymbol(eid, conn=None):
    """Returns a valid gene symbol given an Entrez ID.

    Paramters
    ---------
    eid : An Entrez ID
    
    conn : A MySQL connection (optional if you already called connect())

    Returns
    -------
    symbols : a set of symbols - which will be empty if no valid symbols were found.
    """

    if conn is None:
        if gene_conn is None:
            raise TypeError("getSymbol() takes 2 arguments if connect() has not been called")
        conn = gene_conn
    cursor = conn.cursor()

    valid_eids = checkEID(eid, conn)
    symbols = set()
    for gene in valid_eids:
        cursor.execute("SELECT symbol FROM genes WHERE entrez_id = %(eid)s AND symbol IS NOT NULL", {'eid': gene})
        symbols.update([row[0] for row in cursor.fetchall()])

    return symbols

def crossQuery(db, id, conn=None):
    """Returns valid Entrez IDs given a foreign database and identifier.

    Paramters
    ---------
    cursor : a MySQL cursor to the gene database

    id : the foregin identifier
    
    conn : A MySQL connection (optional if you already called connect())

    Returns
    -------
    eids : a set of matching Entrez IDs - which will be empty if no matching ids were found.
    """
    
    if conn is None:
        if gene_conn is None:
            raise TypeError("crossQuery() takes 2 arguments if connect() has not been called")
        conn = gene_conn
    cursor = conn.cursor()

    cursor.execute("SELECT entrez_id FROM gene_xrefs WHERE Xref_db = %(db)s AND Xref_id = %(id)s", {'db': db, 'id': id})
    return set([row[0] for row in cursor.fetchall()])
