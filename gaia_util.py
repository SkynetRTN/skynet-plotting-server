from os import error
import sqlite3

from numpy import mat
from kdtree import Star, Star_tree
import math


def gaia_get_data(range):
    try:
        dec_min = float(range['mindec'])
        dec_max = float(range['maxdec'])
        ra_min = float(range['minra'])
        ra_max = float(range['maxra'])
    except:
        return {"error": "Input invalid type"}

    where = 'ra >= ? and ra <= ? and dec >= ? and dec <= ?'
    args = (ra_min, ra_max, dec_min, dec_max)
    sqlite_filename = 'gaia_clusters.sqlite'
    conn = sqlite3.connect(sqlite_filename)
    try:
        # Query RA/Dec region(s) using constraints defined above (should be
        # fast thanks to the indexes) in a subquery; the outer query returns
        # only sources within the circle of radius r using the haversine
        # formula, which is more accurate for small distances
        cur = conn.cursor()
        sources = cur.execute(
            'select * from (select * from gedr3dis where ' + where + ')', args).fetchall()
    except error as e:
        print(e)
    finally:
        conn.close()
    # Output sources in CSV
    # for source in sources:
        # source_id, ra, dec, r, pmra, pmdec
        # print(','.join(str(x) for x in source))
    return sources


# print(gaia_get_data({'minra': 90, 'maxra': 120, 'mindec': 45, 'maxdec': 50}))

def gaia_match(photometry, range):
    gaia_data = gaia_get_data(range)
    kd_tree = Star_tree(Star(gaia_data[0][1:]))
    for entry in gaia_data[1:]:
        kd_tree.insert(Star(entry[1:]))
    result = []
    for entry in photometry:
        target = Star([entry['ra'], entry['dec']]+[0, 0, 0])
        match = kd_tree.nn(target)
        if match[1] < 1.45444*10**(-5):
            dist = match[1]**0.5*2/math.pi*180
            result.insert({'id': entry['id'], 'range': dist, 'pm': {
                          'ra': match[0].pmra, 'dec': match[0].pmdec}})
            kd_tree.delete(match[0])
    return result


# print(gaia_match(
#     {'id': "src10", 'ra': 9.6, 'dec': -72}, {'minra': 0, 'maxra': 20, 'mindec': -90, 'maxdec': -50}))
