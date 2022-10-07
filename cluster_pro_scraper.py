from astroquery.simbad import Simbad
from astroquery.vizier import Vizier
from astropy.coordinates import Angle
import astropy
import astropy.units as u
import astropy.coordinates as coord
import sqlite3

from gaia_util import gaia_match


def scraper_query_object(query: str):
    simbad = Simbad()
    simbad.add_votable_fields('bibcodelist(2003-2013)')
    raw = simbad.query_object(query)
    return {'RA': hms2d(raw['RA'][0]), 'DEC': dms2d(raw['DEC'][0]), 'Range': 0}


def scraper_query_object_local(query: str):
    simbad_result = scraper_query_object(query)
    simbad_ra = simbad_result['RA']
    simbad_dec = simbad_result['DEC']
    sqlite_filename = 'MWSC.sqlite'
    conn = sqlite3.connect(sqlite_filename)
    delta = 0.05
    try:
        cur = conn.cursor()
        cluster = cur.execute('SELECT * FROM MWSC WHERE ra >= ? AND ra <=? AND dec >= ? AND dec <= ?',
                              [simbad_ra - delta, simbad_ra + delta, simbad_dec - delta, simbad_dec + delta]).fetchall()
        if cluster:
            cluster = cluster[0]
            result = {'RA': cluster[2], "DEC": cluster[3], "Range": cluster[5]}
        else:
            result = simbad_result
    except Exception as e:
        result = simbad_result
    return result


# print(scraper_query_vizier(115.44,-14.80,0.177))
def scraper_query_vizier(ra: float, dec: float, r: float, catalog) -> astropy.table:
    query_coords = coord.SkyCoord(ra=ra, dec=dec, unit=(u.deg, u.deg), frame='icrs')
    # Skymapper by default does not report the error columns
    columns = []
    catalog_vizier = []
    filters = []
    if 'gaia' in catalog:
        filters += ['G', 'BP', 'RP']
        columns += ['Source', 'RA_ICRS', 'DE_ICRS', 'pmRA', 'e_pmRA', 'pmDE', 'e_pmDE', 'Dist']
        columns += ['Gmag', 'e_Gmag', 'BPmag', 'e_BPmag', 'RPmag', 'e_RPmag']
        catalog_vizier.append('I/355/gaiadr3')
    if '2mass' in catalog:
        filters += ['J', 'H', 'K']
        columns += ['RAJ2000', 'DEJ2000']
        columns += ['Jmag', 'e_Jmag', 'Hmag', 'e_Hmag', 'Kmag', 'e_Kmag']
        catalog_vizier.append('II/246/out')

    if catalog:
        vquery = Vizier(columns=columns, row_limit=30000)
        query = vquery.query_region(query_coords, radius=Angle(r * u.deg), catalog=catalog_vizier)
        query = query[0]
        result = []
        for row in query:
            result_row = dict(id=str(row['Source']), isValid=True)
            for filt in filters:
                result_row[filt + 'Mag'] = float(row[filt + 'mag'])
                result_row[filt + 'err'] = float(row['e_' + filt + 'mag'])
                result_row[filt + 'ra'] = float(row['RA_ICRS'])
                result_row[filt + 'dec'] = float(row['DE_ICRS'])
                result_row[filt + 'pmra'] = float(row['pmRA'])
                result_row[filt + 'pmdec'] = float(row['pmDE'])
                result_row[filt + 'dist'] = float(row['Dist'])
            result.append(result_row)
        return {'data': result, 'filters': filters}
    else:
        return []

def datatable_to_gaiatable(data, filters, star_range):
    query = []
    for row in data:
        for filter_name in filters:
            if row[filter_name + 'ra'] and row[filter_name + 'dec']:
                query.append({'id': row['id'], 'ra': row[filter_name + 'ra'], 'dec': row[filter_name + 'dec']})
                break
    gaia_result = gaia_match(query, star_range)
    table_index = 0
    gaia_index = 0

    result = data
    while table_index < len(result) and gaia_index < len(gaia_result):
        if result[table_index]['id'] == gaia_result[gaia_index]['id']:
            for filter_name in filters:
                result[table_index][filter_name + 'pmra'] = gaia_result[gaia_index]['pm']['ra']
                result[table_index][filter_name + 'pmdec'] = gaia_result[gaia_index]['pm']['dec']
                result[table_index][filter_name + 'dist'] = gaia_result[gaia_index]['range']
            table_index += 1
            gaia_index += 1
        else:
            result[table_index][filter_name + 'pmra'] = None
            result[table_index][filter_name + 'pmdec'] = None
            result[table_index][filter_name + 'dist'] = None
            table_index += 1

    return result


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


# data = scraper_query_vizier_gaia(ra=115.45, dec=-14.80, r=0.570)
#
# print(data)
