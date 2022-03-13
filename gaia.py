#!/usr/bin/env python3

"""
Example code for querying the local sqlite database of Gaia stars around
galactic clusters

Usage: gaia_clusters_query.py RA(h|hh:mm:ss) Dec(d|+dd:mm:ss) [r(arcmins)]
"""

from copy import error
import numbers
import sqlite3
import sys

from numpy import arcsin, cos, deg2rad, rad2deg, sin


def gaia_args_verify(args: list):
    error_message = ""
    if not 0 <= args[0] <= 360:
        error_message += "Expected RA in the range [0, 360] "
    if not -90 <= args[1] <= 90:
        error_message += "Expected Dec in the range [-90, +90] "
    if not 0 < args[2] < 90:
        error_message += "Expected query radius in the range (0, 90)"
    return error_message

# noinspection SqlResolve,SqlDerivedTableAlias


def main():
    # Parse RA/Dec, both will be in degrees
    try:
        ra = float(sys.argv[1])
    except IndexError:
        print('Missing RA', file=sys.stderr)
        sys.exit(1)
    except ValueError:
        h, m, s = sys.argv[1].split(':')
        ra = int(h) + int(m)/60 + float(s)/3600
    if not 0 <= ra < 24:
        print('Expected RA in the range [0, 24]', file=sys.stderr)
        sys.exit(1)
    ra *= 15
    try:
        dec = float(sys.argv[2])
    except IndexError:
        print('Missing Dec', file=sys.stderr)
        sys.exit(2)
    except ValueError:
        d, m, s = sys.argv[2].split(':')
        dec = (1 - 2*d.startswith('-'))*(
            abs(int(d)) + int(m)/60 + float(s)/3600)
    if not -90 <= dec < 90:
        print('Expected Dec in the range [-90, +90]', file=sys.stderr)
        sys.exit(2)
    try:
        r = float(sys.argv[3])/60
    except IndexError:
        r = 1/6
    if not 0 < r < 90:
        print('Expected query radius in the range (0, 90)', file=sys.stderr)
        sys.exit(2)

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

    # Connect to the local Gaia subset database
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
    for source in sources:
        # source_id, ra, dec, r, pmra, pmdec
        print(','.join(str(x) for x in source))


if __name__ == '__main__':
    main()
