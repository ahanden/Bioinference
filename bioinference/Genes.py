#!/bin/python

import MySQLdb
import warnings
from MySQLdb import connections

class Gene:

    __output_key = "eid"

    def __init__(self, eid, db):
        try:
            self.eid = int(eid)
        except ValueError:
            raise ValueError("Entrez ID must be an integer")

        if not isinstance(db, GeneDB):
            raise TypeError("db must be a GeneDB object, not %s" % type(db))

        self.db  = db
        self.cross = {}

    @staticmethod
    def set_output_source(source):
        """Defines what should be displayed when converting the gene to a
        string.

        Parameters
        ----------
        source : The key to use for the Gene's output. Same as __get__
        """

        Gene.__output_key = source

    def __repr__(self):
        return str(self[Gene.__output_key])
    __str__ = __repr__

    def __getitem__(self, key):
        if key == "eid":
            return self.eid
        elif key == "symbol":
            try:
                return self.symbol
            except AttributeError:
                if self.db is not None:
                    self.symbol = self.db.get_symbol(self.eid)
                    return self.symbol
                else:
                    raise KeyError("Gene does not have a database " + \
                                   "connection to look up its symbol")
        else:
            if key not in self.cross:
                if self.db is not None:
                    self.cross[key] = self.db.get_cross_id(self.eid, key)
                else:
                    raise KeyError("Gene does not have a database " + \
                                   "connection to look up its symbol")
            return self.cross[key]
               
    def __eq__(self, other):
        return self["eid"] == other["eid"]

    def __ne__(self, other):
        return self["eid"] != other["eid"]

    def __hash__(self):
        return hash((type(self), self["eid"]))

class GeneDB:
    def __init__(self, *args, **kwargs):
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

        if len(args) == 1 or "conn" in kwargs:
            conn = args[0] if len(args) == 1 else kwargs["conn"]

            if type(args[0]) is connections.Connection:
                self.conn = conn
            else:
                raise TypeError("connect() requires a MySQLdb.conn object")

        elif len(args) == 4:
            self.conn = MySQLdb.connect(
                host=args[0],
                db=args[1],
                user=args[2],
                passwd=args[3])
        else:
            try:
                self.conn = MySQLdb.connect(
                    host=kwargs["host"],
                    db=kwargs["db"],
                    user=kwargs["user"],
                    passwd=kwargs["passwd"])
            except KeyError:
                raise TypeError("You must pass host, db, user, and passwd " + \
                                "to connect()")

        self.cursor = self.conn.cursor()
        self.gene_cache = {}
        self.source_cache = {}


    def get_gene(self, gene_id, source="eid"):
        """Returns a Gene object from the given ID and type.

        Parameters
        ----------
        gene_id : A gene identifier

        source : The source of the identifier. Valid options are eid,
                 symbol, or any cross database name in the gene_xrefs table.

        Returns
        -------
        A Gene object
        """

        gene_id = str(gene_id)

        
        try:
            valid_gene = self.gene_cache[self.source_cache[source][gene_id]]
            return valid_gene
        except KeyError:
            pass
        
        valid_eid = None

        if source == "eid":
            try:
                eid = int(gene_id)
            except ValueError:
                raise ValueError("gene_id must be an integer if source " + \
                                 "is \"Entrez\"")

            self.cursor.execute("""
                SELECT EXISTS(
                    SELECT * 
                    FROM genes 
                    WHERE entrez_id = %(eid)s
                )""", {'eid': eid})
            if self.cursor.fetchone()[0] == 1:
                valid_eid = eid

            else:

                self.cursor.execute("""
                    SELECT entrez_id
                    FROM discontinued_genes
                    WHERE discontinued_id = %(eid)s""", {'eid': eid})
                row = self.cursor.fetchone()
                if row is not None:
                    valid_eid = row[0]
                else:
                    raise KeyError("Entrez ID %d was not found in the database" % eid)

        elif source == "symbol":

            args = {"symbol": gene_id}
            self.cursor.execute("""
                SELECT entrez_id
                FROM genes
                WHERE symbol = %(symbol)s""", args)
            row = self.cursor.fetchone()
            if row is not None:
                valid_eid = row[0]
            else:
                self.cursor.execute("""
                    SELECT entrez_id
                    FROM discontinued_genes
                    WHERE discontinued_symbol = %(symbol)s""", args)
                row = self.cursor.fetchone()
                if row is not None:
                    valid_eid = row[0]
                else:
                    self.cursor.execute("""
                        SELECT entrez_id
                        FROM gene_synonyms
                        WHERE symbol = %(symbol)s""", args)
                    row = self.cursor.fetchone()
                    if row is not None:
                        valid_eid = row[0]
                    else:
                        raise KeyError("Symbol %s not found in the database" % gene_id)
        else:
            self.cursor.execute("""
                SELECT entrez_id
                FROM gene_xrefs
                WHERE Xref_db = %(db)s
                AND Xref_id = %(id)s""", {'db': source, 'id': gene_id})
            row = self.cursor.fetchone()
            if row is not None:
                valid_eid = row[0]
            else:
                raise KeyError(("Gene ID %s from source %s was not found " + \
                                "in the database") % (gene_id, source))

        if valid_eid is None:
            raise KeyError("Unable to find a valid Entrez ID for %s from %s" % (gene_id, source))

        valid_eid = int(valid_eid)
        if source not in self.source_cache:
            self.source_cache[source] = {}
        self.source_cache[source][gene_id] = valid_eid
        self.gene_cache[valid_eid] = Gene(valid_eid, self)

        return self.gene_cache[valid_eid]

    def get_cross_id(self, entrez_id, xref_db):
        """Finds a gene's identifier from an outside database based on its
        Entrez ID.

        Paramters
        ---------
        entrez_id : An Entrez ID for a gene

        xref_db : The name of the external database to query from

        Returns
        -------
        A string containing the foreign identifier
        """
        
        try:
            entrez_id = int(entrez_id)
        except ValueError:
            raise ValueError("entrez_id must be an integer")

        self.cursor.execute("""
            SELECT entrez_id
            FROM gene_xrefs
            WHERE Xref_db = %(db)s
            AND entrez_id = %(eid)s""", {'db': xref_db, 'eid': entrez_id})
        row = self.cursor.fetchone()
        if row is not None:
            return row[0]
        
        raise KeyError("Unable to find an external identifer for database " + \
                       "%s using Entrez ID %d" % (xref_db, entrez_id))

    def get_symbol(self, entrez_id):
        """Finds a gene's official symbol (if it has one) given its Entrez ID.

        Parameters
        ----------
        entrez_id : A gene's Entrez ID

        Returns
        -------
        The gene's symbol (if it has one)
        """

        try:
            entrez_id = int(entrez_id)
        except ValueError:
            raise ValueError("entrez_id must be an integer")


        self.cursor.execute("""
            SELECT symbol
            FROM genes
            WHERE entrez_id = %(eid)s""", {'eid': entrez_id})
        row = self.cursor.fetchone()
        if row is not None:
            return row[0]
        raise KeyError("Entrez ID %d was not found in the database" % entrez_id)
