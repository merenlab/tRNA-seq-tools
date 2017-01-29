# -*- coding: utf-8
# pylint: disable=line-too-long

import sqlite3

import tRNASeqTools.filesnpaths as filesnpaths
from tRNASeqTools.errors import ConfigError

__author__ = "Steven Cui"
__copyright__ = "Copyright 2017, Meren Lab"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = 0.1
__maintainer__ = "Steven Cui"
__email__ = "stevencui729@gmail.com"


class DB:
    """This class handles basic database functions."""

    def __init__(self, db_path, client_version, new_database=False, ignore_version=False):
        """Initializes path and connects database."""
        self.db_path = db_path
        self.version = None

        if new_database:
            filesnpaths.is_output_file_writable(db_path)
        else:
            filesnpaths.is_file_exists(db_path)

        self.conn = sqlite3.connect(self.db_path)
        self.conn.text_factory = str

        self.cursor = self.conn.cursor()

        if new_database:
            self.create_self()
            self.set_version(client_version)
        else:

            self.version = self.get_version()
            if str(self.version) != str(client_version) and not ignore_version:
                raise ConfigError("It seems the database '%s' was generated when your client was at version %s,\
                                    however, your client now is at version %s. Which means this database file\
                                    cannot be used with this client anymore and needs to be upgraded to the\
                                    version %s :/"\
                                            % (self.db_path, self.version, client_version, client_version, self.version))


    def set_version(self, version):
        self.set_meta_value('version', version)
        self.commit()


    def get_version(self):
        try:
            return self.get_meta_value('version')
        except:
            raise ConfigError("%s does not seem to be a valid tRNA-seq database :/" % self.db_path)


    def set_meta_value(self, key, value):
        self._exec('''INSERT INTO self VALUES(?,?)''', (key, value,))
        self.commit()


    def get_meta_value(self, key):
        response = self._exec("""SELECT value FROM self WHERE key='%s'""" % key)
        rows = response.fetchall()

        if not rows:
            raise ConfigError("A value for '%s' does not seem to be set in table 'self'." % key)

        val = rows[0][0]

        if isinstance(val, type(None)):
            return None

        try:
            val = int(val)
        except ValueError:
            pass

        return val


    def _exec_many(self, sql_query, values):
        return self.cursor.executemany(sql_query, values)


    def _exec(self, sql_query, value=None):
        """Overwrite function used in place of execute() to execute sql
        queries.
        """
        if value:
            return_val = self.cursor.execute(sql_query, value)
        else:
            return_val = self.cursor.execute(sql_query)

        self.commit()
        return return_val


    def create_self(self):
        """Creates an empty default table."""
        self._exec("""CREATE TABLE self (key text, value text)""")
        self.commit()


    def create_table(self, table_name, fields, types):
        """Creates a table with the arguments given."""
        if len(fields) != len(types):
            print("error: fields and types different sizes")

        db_fields = ", ".join(["%s %s" % (t[0], t[1]) for t in zip(fields,
            types)])

        self._exec("""CREATE TABLE IF NOT EXISTS %s (%s)""" % (table_name,
            db_fields))
        self.commit()


    def get_all_rows_from_table(self, table):
        """Get all the rows from a table as a list."""
        result = self._exec("""SELECT * FROM %s""" % table)
        return result.fetchall()


    def get_table_structure(self, table):
        """Get the headers of a table as a list."""
        result = self._exec("""SELECT * FROM %s""" % table)
        return [t[0] for t in result.description]


    def get_table_as_dict(self, table, table_structure=None,
        string_the_key=False, keys_of_interest=None,
        omit_parent_column=False):
        """Returns a table's entire contents as a dict."""

        results_dict = {}
        rows = self.get_all_rows_from_table(table)

        if not table_structure:
            table_structure = self.get_table_structure(table)

        columns_to_return = list(range(0, len(table_structure)))

        if omit_parent_column:
            if "__parent__" in table_structure:
                columns_to_return.remove(table_structure.index("__parent__"))

        if len(columns_to_return) == 1:
            print("nothing left to return")

        for row in rows:
            entry = {}

            for i in columns_to_return[1:]:
                entry[table_structure[i]] = row[i]

            if string_the_key:
                results_dict[str(row[0])] = entry
            else:
                results_dict[row[0]] = entry

        return results_dict


    def get_some_rows_from_table_as_dict(self, table, where_clause,
        string_the_key=False):
        """Returns some of a table's contents as a dict. The contents are
        determined by the where-clause arguement, which will be used in an SQL
        query.
        """

        results_dict = {}
        rows = self._exec("""SELECT * FROM %s WHERE %s""" % (table,
            where_clause)).fetchall()
        table_structure = self.get_table_structure(table)
        columns_to_return = list(range(0, len(table_structure)))

        for row in rows:
            entry = {}

            for i in columns_to_return[1:]:
                entry[table_structure[i]] = row[i]

            if string_the_key:
                results_dict[str(row[0])] = entry
            else:
                results_dict[row[0]] = entry

        return results_dict


    def commit(self):
        """Commits the changes made to the database."""
        self.conn.commit()


    def disconnect(self):
        """Disconnects from the database."""
        self.conn.commit()
        self.conn.close()

