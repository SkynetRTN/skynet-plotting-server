from astroquery.simbad import Simbad
from astroquery.vizier import Vizier
from astropy.coordinates import Angle
import astropy
import astropy.units as u
import astropy.coordinates as coord
import sqlite3


def scraper_query_object(query: str):
    simbad = Simbad()
    simbad.add_votable_fields('bibcodelist(2003-2013)')
    raw = simbad.query_object(query)
    return {'RA': hms2d(raw['RA'][0]), 'DEC': dms2d(raw['DEC'][0]), 'Range': 0}


def scraper_query_object_local(query: str):
    simbad_result = scraper_query_object(query)
    simbad_ra = simbad_result['RA']
    simbad_dec = simbad_result['DEC']
    print(simbad_result)
    sqlite_filename = 'MWSC.sqlite'
    conn = sqlite3.connect(sqlite_filename)
    delta = 0.05
    try:
        cur = conn.cursor()
        cluster = cur.execute('SELECT * FROM MWSC WHERE ra >= ? AND ra <=? AND dec >= ? AND dec <= ?',
                              [simbad_ra - delta, simbad_ra + delta, simbad_dec - delta, simbad_dec + delta]).fetchall()
        print(cluster)
        if cluster:
            cluster = cluster[0]
            result = {'RA': cluster[2], "DEC": cluster[3], "Range": cluster[5]}
        else:
            result = simbad_result
    except Exception as e:
        print(e)
        result = simbad_result
    return result


# print(scraper_query_vizier(115.44,-14.80,0.177))
def scraper_query_vizier(ra: float, dec: float, r: float) -> astropy.table:
    query_coords = coord.SkyCoord(ra=ra, dec=dec, unit=(u.deg, u.deg), frame='icrs')
    result = Vizier.query_region(query_coords, radius=Angle(r * u.deg), catalog='I/355')
    return result[0]


def scraper_query_vizier_gaia(ra: float, dec: float, r: float):
    t = scraper_query_vizier(ra, dec, r)
    result = []
    for row in t:
        result.append(dict(id=str(row['Source']),
                           isValid=True,
                           GMag=float(row['Gmag']),
                           Gerr=float(row['e_FG']),
                           Gra=float(row['RA_ICRS']),
                           Gdec=float(row['DE_ICRS']),
                           GBPMag=float(row['BPmag']),
                           GBPerr=float(row['e_FBP']),
                           GBPra=float(row['RA_ICRS']),
                           GBPdec=float(row['DE_ICRS']),
                           GRPMag=float(row['RPmag']),
                           GRPerr=float(row['e_FRP']),
                           GRPra=float(row['RA_ICRS']),
                           GRPdec=float(row['DE_ICRS']),
                           # B=float(row['Gmag']),
                           # Berr=float(row['e_FG']),
                           # Bra=float(row['RA_ICRS']),
                           # Bdec=float(row['DE_ICRS']),
                           # V=float(row['BPmag']),
                           # Verr=float(row['e_FBP']),
                           # Vra=float(row['RA_ICRS']),
                           # Vdec=float(row['DE_ICRS']),
                           # R=float(row['RPmag']),
                           # Rerr=float(row['e_FRP']),
                           # Rra=float(row['RA_ICRS']),
                           # Rdec=float(row['DE_ICRS']),
                           ))
    return {'data': result, 'filters': ['G', 'GRP', 'GBP']}


# scraper_query_vizier_gaia(115.44, -14.80, 0.177)

def hms2d(hms: str):
    hms = hms.split(' ')
    result = float(hms[0]) * 15 + float(hms[1]) * 0.25
    if len(hms) == 3:
        result += float(hms[2]) * 0.25 / 60
    return result


def dms2d(dms: str):
    dms = dms.split(' ')
    result = float(dms[0])
    is_positive = True
    if result < 0:
        result = -result
        is_positive = False
    result += float(dms[1]) / 60
    if len(dms) == 3:
        result += float(dms[2]) / 3600
    if is_positive:
        return result
    else:
        return -result


# t.show_in_browser()

# with open('scraper_sample.txt', 'w') as f:
#     f.write(str(result))
# f.close()

# scraper_query_object_local('M2')
def table_to_python(table):
    """Convert Astropy Table to Python dict.

    Numpy arrays are converted to lists, so that
    the output is JSON serialisable.

    Can work with multi-dimensional array columns,
    by representing them as list of list.
    """
    total_data = {}
    for name in table.colnames:
        data = table[name].tolist()
        total_data[name] = data
    return total_data
