# -*- coding: utf-8
# pylint: disable=line-too-long


import os
import sqlite3

__author__ = "Steven Cui"
__copyright__ = "Copyright 2016, The University of Chicago"
__credits__ = []
__license__ = "GPL 3.0"
__version__ = 0.1
__maintainer__ = "Steven Cui"
__email__ = "stevencui729@gmail.com"


class DB:
    def __init__(self, db_path):
        self.db_path = db_path

        self.conn = sqlite3.connect(self.db_path)
        print "successfully connected"

        self.cursor = self.conn.cursor()


    def _exec(self, sql_query, value=None):
        if value:
            return_val = self.cursor.execute(sql_query, value)
        else:
            return_val = self.cursor.execute(sql_query)

        self.commit()
        return return_val


    def create_self(self):
        self._exec("""CREATE TABLE self (key text, value text)""")

    def create_table(self, table_name, fields, types):
        if len(fields) != len(types):
            print "error: fields and types different sizes"

        db_fields = ", ".join(["%s %s" % (t[0], t[1]) for t in zip(fields,
            types)])


        self._exec("""CREATE TABLE IF NOT EXISTS %s (%s)""" % (table_name,
            db_fields))
        self.commit()


    def get_all_rows_from_table(self, table):
        result = self._exec("""SELECT * FROM %s""" % table)
        return result.fetchall()

    def get_table_structure(self, table):
        result = self._exec("""SELECT * FROM %s""" % table)
        return [t[0] for t in result.description]

    def get_table_as_dict(self, table, table_structure=None,
        string_the_key=False, keys_of_interest=None,
        omit_parent_column=False):
       
        results_dict = {}
        rows = self.get_all_rows_from_table(table)

        if not table_structure:
            table_structure = self.get_table_structure(table)

        columns_to_return = range(0, len(table_structure))

        if omit_parent_column:
            if "__parent__" in table_structure:
                columns_to_return.remove(table_structure.index("__parent__"))

        if len(columns_to_return) == 1:
            print "nothing left to return"

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

        results_dict = {}
        rows = self._exec("""SELECT * FROM %s WHERE %s""" % (table,
            where_clause)).fetchall()
        table_structure = self.get_table_structure(table)
        columns_to_return = range(0, len(table_structure))

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
        self.conn.commit()

    def disconnect(self):
        self.conn.commit()
        self.conn.close()

