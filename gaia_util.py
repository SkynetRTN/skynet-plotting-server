from os import error
import sqlite3


from numpy import arcsin, cos, deg2rad, rad2deg, sin
from kdtree import Star, Star_tree
import math


def gaia_get_data(range):
    try:
        dec = float(range['dec'])
        ra = float(range['ra'])
        r = float(range['r'])
    except:
        return {"error": "Input invalid type"}
    if not 0 <= ra < 360:
        return {'error': 'Expected RA in the range [0, 360)'}
    if not -90 <= dec < 90:
        return {'error': 'Expected Dec in the range [-90, +90]'}
    if not 0 < r < 90:
        return {'error': 'Expected query radius in the range (0, 90)'}

    # Compute the RA/Dec query ranges; handle poles and RA=0/360 wrap
    dec_min, dec_max = dec - r, dec + r
    if dec_min < -90:
        # South Pole in FOV, use the whole RA range
        where = 'dec <= ?'
        args = (dec_max,)
    elif dec_max > 90:
        # North Pole in FOV, use the whole RA range
        where = 'dec >= ?'
        args = (dec_min,)
    else:
        # See http://janmatuschek.de/LatitudeLongitudeBoundingCoordinates
        dra = rad2deg(arcsin(sin(deg2rad(r))/cos(deg2rad(dec))))
        ra_min, ra_max = ra - dra, ra + dra
        if ra_max >= ra_min + 360:
            # RA spans the whole 360 range
            where = 'dec >= ? and dec <= ?'
            args = (dec_min, dec_max)
        elif ra_min < 0:
            # RA range encloses RA=0 => two separate RA ranges:
            # ra_min + 360 <= ra <= 360 and 0 <= ra <= ra_max
            where = '(ra >= ? or ra <= ?) and dec >= ? and dec <= ?'
            args = (ra_min + 360, ra_max, dec_min, dec_max)
        elif ra_max > 360:
            # RA range encloses RA=360 => two separate RA ranges:
            # ra_min <= ra <= 360 and 0 <= ra <= ra_max - 360
            where = '(ra >= ? or ra <= ?) and dec >= ? and dec <= ?'
            args = (ra_min, ra_max - 360, dec_min, dec_max)
        else:
            # RA range fully within [0, 360)
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
            'select * from (select * from gedr3dis where ' + where +
            ') where asin(sqrt(pow(sin(radians(dec - ?)/2), 2) + '
            'pow(sin(radians(ra - ?)/2), 2)*cos(radians(dec))*?)) <= ?',
            args + (dec, ra, cos(deg2rad(dec)), deg2rad(r)/2)
        ).fetchall()
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
    print(len(gaia_data))
    kd_tree = Star_tree(Star(gaia_data[0][1:]))
    for entry in gaia_data[1:]:
        kd_tree.insert(Star(entry[1:]))
    result = []
    count = 0
    print(len(photometry))
    for entry in photometry:
        target = Star([entry['ra'], entry['dec']]+[0, 0, 0])
        match = kd_tree.nn(target)
        # print(match[1])
        if match[1] < 1.45444*10**(-5):
            # if match[1] < 5.45444*10**(-5):
            dist = match[1]**0.5*2/math.pi*180
            result.append({'id': entry['id'], 'range': dist, 'pm': {
                          'ra': match[0].pmra, 'dec': match[0].pmdec}})
            kd_tree.delete(match[0])
        count += 1
    return result


# print(gaia_match(
    # [{'id': "src10", 'ra': 9.6, 'dec': -72}], {'minra': 0, 'maxra': 20, 'mindec': -90, 'maxdec': -50}))

print(
    len(gaia_get_data({'minra': 0, 'maxra': 290, 'mindec': 3, 'maxdec': 9})))
