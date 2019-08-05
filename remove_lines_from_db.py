#!/usr/bin/env python3

import os
import sqlite3 as lite
try:
    DB_FILE = os.environ['SEARCH_DB']
except KeyError:
    print("environmental variable CMD_BS_DB_DEF_FILE must be defined")

con = lite.connect(DB_FILE)
con.execute("VACUUM")
with con:
    cur = con.cursor()
    """
    for rownum in range(899,6000):
        print(rownum)
        cur.execute("DELETE from PulsarSearch where Rownum=?", (rownum,))
        cur.execute("DELETE from Beamform where BSID=?", (rownum,))
        cur.execute("DELETE from Prepdata where BSID=?", (rownum,))
        cur.execute("DELETE from FFT where BSID=?", (rownum,))
        cur.execute("DELETE from Accel where BSID=?", (rownum,))
        cur.execute("DELETE from Fold where BSID=?", (rownum,))
        cur.execute("DELETE from Candidates where BSID=?", (rownum,))
        #if rownum%100 == 0:
        #    con.execute("VACUUM")
    """



