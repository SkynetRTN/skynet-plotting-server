from os import error
import os
from platform import node
import sqlite3


from numpy import arcsin, cos, deg2rad, rad2deg, sin
import math
import numpy as np

from pro_tree import tree_matching_sci


def gaia_get_data(range):
    try:
        dec = float(range['dec'])
        ra = float(range['ra'])
        r = float(range['r'])
    except:
        raise error({"error": "Input invalid type"})
    if not 0 <= ra < 360:
        raise error({'error': 'Expected RA in the range [0, 360)'})
    if not -90 <= dec < 90:
        raise error({'error': 'Expected Dec in the range [-90, +90]'})
    if not 0 < r < 90:
        raise error({'error': 'Expected query radius in the range (0, 90)'})

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
    sqlite_filename = os.path.join(
        os.path.dirname(__file__), 'gaia_clusters.sqlite')
    conn = sqlite3.connect(sqlite_filename)
    try:
        # Query RA/Dec region(s) using constraints defined above (should be
        # fast thanks to the indexes) in a subquery; the outer query returns
        # only sources within the circle of radius r using the haversine
        # formula, which is more accurate for small distances
        cur = conn.cursor()
        conn.create_function('asin', 1, math.asin)
        conn.create_function('sqrt', 1, math.sqrt)
        conn.create_function('sin', 1, math.sin)
        conn.create_function('cos', 1, math.cos)
        conn.create_function('radians', 1, math.radians)
        conn.create_function('pow', 2, math.pow)
        sources = cur.execute(
            'select * from (select * from gedr3dis where ' + where +
            ') where asin(sqrt(pow(sin(radians(dec - ?)/2), 2) + '
            'pow(sin(radians(ra - ?)/2), 2)*cos(radians(dec))*?)) <= ?',
            args + (dec, ra, cos(deg2rad(dec)), deg2rad(r)/2)
        ).fetchall()
    except:
        raise error({'error': "Cannot Pull GAIA Database"})
    finally:
        conn.close()
    # Output sources in CSV
    # for source in sources:
        # source_id, ra, dec, r, pmra, pmdec
        # print(','.join(str(x) for x in source))
    return sources


# print(gaia_get_data({'minra': 90, 'maxra': 120, 'mindec': 45, 'maxdec': 50}))
def convert_gaia(data, range):
    ra = math.radians(data[1])
    dec = math.radians(data[2])
    return convert(ra, dec, range)


def convert_usr(data, range):
    ra = math.radians(data['ra'])
    dec = math.radians(data['dec'])
    return convert(ra, dec, range)


def convert(ra, dec, range):
    if range['wrap']:
        ra = ra - 2 * math.pi
    ra_c = math.radians(range['ra'])
    dec_c = math.radians(range['dec'])
    x = math.cos(dec) * math.sin(ra-ra_c)
    y = -math.sin(dec) * math.cos(ra-ra_c) * math.sin(dec_c) + \
        math.cos(dec-dec_c) * math.sin(dec_c)
    return [x, y]


def haversine(dec1, dec2, ra1, ra2):
    # approx of (theta/2)^2
    dec1 = math.radians(dec1)
    dec2 = math.radians(dec2)
    ra1 = math.radians(ra1)
    ra2 = math.radians(ra2)
    return math.sin((dec1 - dec2)/2)**2 + math.cos(dec1) * math.cos(dec2) * (math.sin((ra1 - ra2)/2))**2


def gaia_match(photometry, star_range):
    try:
        gaia_data: list[dict] = gaia_get_data(star_range)
    except error as e:
        raise error(e)
    try:
        nodes = []
        for data in gaia_data:
            nodes.append(convert_gaia(data, star_range))
        entrys = []
        for entry in photometry:
            entrys.append(convert_usr(entry, star_range))
    except:
        raise error({'error': 'Cannot Convert User Data for GAIA Matching '})
    try:
        nodes = np.array(nodes)
        entrys = np.array(entrys)
        matched = tree_matching_sci(nodes, entrys)
    except:
        raise error({'error': 'KD-Tree Failure'})
    result = []
    for i in range(len(photometry)):
        query = photometry[i]
        gaia = gaia_data[matched[i]]
        if haversine(query['dec'], gaia[2], query['ra'], gaia[1]) < 1.45444*10**(-5):
            result.append({'id': query['id'], 'range': gaia[3], 'pm': {
                          'ra': gaia[4], 'dec': gaia[5]}})
    return result


# print(gaia_match(
    # [{'id': "src10", 'ra': 9.6, 'dec': -72}], {'minra': 0, 'maxra': 20, 'mindec': -90, 'maxdec': -50}))

# print(
#     len(gaia_get_data({'minra': 0, 'maxra': 290, 'mindec': 3, 'maxdec': 9})))
