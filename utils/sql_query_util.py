"""Simplify operations against Jan 2022 genegraph.db SQLite file """
import os.path
import sqlite3
from typing import List

DEFAULT_DATABASE_FILE = "/GeneGraphDB/data/genegraph.db"


class SimplifiedSqlDb:
    def __init__(self, database_file=DEFAULT_DATABASE_FILE):
        if not os.path.exists(database_file):
            raise FileNotFoundError(f"Can't find database file '{database_file}'")
        self.conn = sqlite3.connect(database_file)
        self.conn.row_factory = sqlite3.Row

    def get_tables(self) -> List[str]:
        """List available database tables"""
        tables = []

        cur = self.conn.cursor()
        cur.execute("SELECT name FROM sqlite_master WHERE type='table';")
        for table_tup in cur.fetchall():
            tables.append(table_tup[0])
        cur.close()

        return tables

    def get_p100_windowed_neighbors(self, p100_hash) -> List[str]:
        """Get all proteins within 5k bp of specified p100"""
        neighbor_hashes = set()
        cur = self.conn.cursor()
        cur.execute(f"SELECT * FROM prot2protwindow WHERE p1hash is '{p100_hash}' OR p2hash is '{p100_hash}'")
        for row in cur.fetchall():
            for key in row.keys():
                neighbor_hashes.add(row[key])
        cur.close()
        neighbor_hashes.remove(p100_hash)
        return list(sorted(neighbor_hashes))

    def get_p30_cluster_for_p100(self, p100_hash) -> sqlite3.Row:
        cur = self.conn.cursor()
        cur.execute(f"SELECT * FROM clusters WHERE p100 is '{p100_hash}'")
        all_rows = cur.fetchall()
        assert len(all_rows) == 1, f"Found {len(all_rows)} clusters for p100 {p100_hash}"
        cur.close()
        row = all_rows[0]
        return row

    def get_p30_cluster_members(self, p30_hash) -> List[sqlite3.Row]:
        cur = self.conn.cursor()
        cur.execute(f"SELECT * FROM clusters WHERE p30 is '{p30_hash}'")
        rows = list(cur.fetchall())
        assert len(rows) > 0, f"No clusters found for p30 {p30_hash}"
        cur.close()
        return rows
