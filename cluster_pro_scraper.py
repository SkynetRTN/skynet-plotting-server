import os
import time

import astroquery.gaia
import numpy
import numpy as np
from astroquery.simbad import Simbad
from astroquery.vizier import Vizier
from astropy.coordinates import Angle
import astropy
import astropy.units as u
import astropy.coordinates as coord
import sqlite3
import grispy as gsp
import operator
from functools import reduce

from gaia_util import gaia_match

max_row_limit = 50000


def scraper_query_object(query: str):
    simbad = Simbad()
    simbad.add_votable_fields('bibcodelist(2003-2013)')
    raw = simbad.query_object(query)
    return {'RA': hms2d(raw['RA'][0]), 'DEC': dms2d(raw['DEC'][0]), 'Range': 0}


def scraper_query_object_local(query: str):
    simbad_result = scraper_query_object(query)
    simbad_ra = simbad_result['RA']
    simbad_dec = simbad_result['DEC']
    sqlite_filename = os.path.join(
        os.path.dirname(__file__), 'MWSC.sqlite')
    conn = sqlite3.connect(sqlite_filename)
    delta = 0.05
    cur = conn.cursor()
    cluster = cur.execute('SELECT * FROM MWSC WHERE ra >= ? AND ra <=? AND dec >= ? AND dec <= ?',
                          [simbad_ra - delta, simbad_ra + delta, simbad_dec - delta, simbad_dec + delta]).fetchall()
    if cluster:
        cluster = cluster[0]
        result = {'RA': cluster[2], "DEC": cluster[3], "Range": cluster[5]}
    else:
        result = simbad_result
    return result


def scraper_query(coordinates, constrain, catalog, file_keys, file_data):
    time_temp = time.time_ns()

    gaia = scraper_query_gaia_esac(coordinates, constrain)

    # print("GAIA Fetch Finished:")
    # print((time.time_ns() - time_temp)/(10**(-9)))
    # time_temp = time.time_ns()

    result_table = gaia['gaia_table']
    grid = gsp.GriSPy(gaia['np_coordinates'], N_cells=128)
    output_filters = []

    if file_keys:
        file_keys = [key.replace(" ", "") for key in file_keys]
        file_filters = [key[0:-3] for key in file_keys if 'err' in key]
        output_filters += file_filters
        # print(file_keys)

        if file_keys[0] == "id":
            # print(file_data[0])
            # print(file_keys)
            file_keys = ["source", "isValid"] + file_keys[2:]
            key_dtype = tuple(['U', 'U'] + ['<f8' for _ in range(0, len(file_keys) - 2)])
        else:
            # print(file_data[0])
            # print(file_keys)
            file_keys = file_keys[:-2] + ["source", "isValid"]
            key_dtype = tuple(['<f8' for _ in range(0, len(file_keys) - 2)] + ['U', 'U'])

        # file_keys = [key + "mag" if key in file_filters else key for key in file_keys]
        file_keys = [key.replace("Mag", "mag") if "Mag" in key else key for key in file_keys]
        file_keys = ["e_" + key[0:-3] + "mag" if "err" in key else key for key in file_keys]
        # print(file_keys)
        file_table = astropy.table.Table(rows=file_data, names=tuple(file_keys), dtype=key_dtype, masked=True)

        ra_keys = tuple([filt + "ra" for filt in file_filters])
        dec_keys = tuple([filt + "dec" for filt in file_filters])

        coord_dtype = tuple(['<f8' for _ in range(0, len(file_filters))])
        file_ras = astropy.table.Table(file_table[ra_keys], dtype=coord_dtype)
        file_decs = astropy.table.Table(file_table[dec_keys], dtype=coord_dtype)

        np_file_ras = np.array(file_ras[ra_keys[0]])
        np_file_decs = np.array(file_decs[dec_keys[0]])

        for i in range(0, len(ra_keys)):
            ra_key = ra_keys[i]
            dec_key = dec_keys[i]
            if i > 0:
                if i == 1:
                    np_file_ras = np.dstack((np_file_ras, np.array(file_ras[ra_key])))[0]
                    np_file_decs = np.dstack((np_file_decs, np.array(file_decs[dec_key])))[0]
                else:
                    np_file_ras = np.concatenate((np_file_ras, np.array([file_ras[ra_key]]).T), axis=1)
                    np_file_decs = np.concatenate((np_file_decs, np.array([file_decs[dec_key]]).T), axis=1)

        masked_np_file_ras = np.ma.masked_array(np_file_ras, np.isnan(np_file_ras))
        masked_np_file_decs = np.ma.masked_array(np_file_decs, np.isnan(np_file_decs))

        ras = numpy.average(masked_np_file_ras, axis=1)
        decs = numpy.average(masked_np_file_decs, axis=1)

        ra_col = astropy.table.column.Column(data=ras, name='RAJ2000')
        dec_col = astropy.table.column.Column(data=decs, name='DEJ2000')

        # file_table.show_in_browser()

        file_table.add_columns([ra_col, dec_col])

        file_table = file_table[tuple(['RAJ2000', 'DEJ2000'] +
                                      [filt + "mag" for filt in file_filters] +
                                      ["e_" + filt + "mag" for filt in file_filters])]

        mask = np.logical_and(abs(file_table['RAJ2000']) > 0.000001, abs(file_table['DEJ2000']) > 0.000001)

        file_table = file_table[mask]

        # file_table.show_in_browser()

        result_table = gaia_table_matching(grid, result_table, file_table)

    if 'gaia' in catalog:
        gaia_filters = ['G', 'BP', 'RP']
        output_filters += gaia_filters
        gaia_columns = filters_to_columns(gaia_filters, ['RA_ICRS', 'DE_ICRS'])
        gaia_table = scraper_query_vizier(coordinates, gaia_columns, 'I/355/gaiadr3')
        gaia_table.rename_column('RA_ICRS', 'RAJ2000')
        gaia_table.rename_column('DE_ICRS', 'DEJ2000')
        result_table = gaia_table_matching(grid, result_table, gaia_table)

    if 'apass' in catalog:
        apass_filters = ['V', 'B', 'g\'', 'r\'', 'i\'']
        output_filters += apass_filters
        apass_columns = filters_to_columns(apass_filters, ['RAJ2000', 'DEJ2000'])
        apass_table = scraper_query_vizier(coordinates, apass_columns, 'II/336')
        for filt in ['g', 'r', 'i']:
            apass_table.rename_column(filt + '_mag', filt + '\'mag')
            apass_table.rename_column('e_' + filt + '_mag', 'e_' + filt + '\'mag')
        result_table = gaia_table_matching(grid, result_table, apass_table)

    if 'twomass' in catalog:
        two_mass_filters = ['J', 'H', 'K']
        output_filters += two_mass_filters
        two_mass_columns = filters_to_columns(two_mass_filters, ['RAJ2000', 'DEJ2000'])
        two_mass_table = scraper_query_vizier(coordinates, two_mass_columns, 'II/246/out')
        result_table = gaia_table_matching(grid, result_table, two_mass_table)

    if 'wise' in catalog:
        wise_filters = ['W1', 'W2', 'W3', 'W4']
        output_filters += wise_filters
        wise_columns = filters_to_columns(wise_filters, ['RAJ2000', 'DEJ2000'])
        wise_table = scraper_query_vizier(coordinates, wise_columns, 'II/328/allwise')

        result_table = gaia_table_matching(grid, result_table, wise_table)

    output_filters = list(dict.fromkeys(output_filters))  # remove duplicates

    return astropy_table_to_result(result_table, output_filters)


def filters_to_columns(filters, ra_dec):
    result = ra_dec
    for filt in filters:
        result.append(filt + 'mag')
        result.append('e_' + filt + 'mag')
    return result


def gaia_table_matching(grid, table, target_query, join_type='left'):
    target_table = target_query
    target_cord = np.dstack((np.array(target_table['RAJ2000']), np.array(target_table['DEJ2000'])))[0]
    nn_dist, nn_indices = grid.nearest_neighbors(target_cord, n=1)
    nn_dist = np.concatenate(nn_dist)
    nn_indices = np.concatenate(nn_indices)
    nn_indices_filtered = [nn_indices[i] + 1 if nn_dist[i] < 0.000833 else 0 for i in range(0, len(nn_indices))]

    duplicate_cols = list(numpy.intersect1d(np.array(target_table.colnames),
                                            np.intersect1d(np.array(target_table.colnames),
                                                           np.array(table.colnames))))
    for col in duplicate_cols:
        table.rename_column(col, col + '1')
        target_table.rename_column(col, col + '2')

    nn_indices_col = astropy.table.column.Column(data=np.array(nn_indices_filtered), name='id')
    target_table.add_column(nn_indices_col, 0)

    joined_table = astropy.table.join(left=table, right=target_table, keys='id', join_type=join_type)

    joined_table = astropy.table.Table(joined_table, masked=True)

    for col in duplicate_cols:
        if col == 'RAJ2000' or col == 'DEJ2000':
            average = np.add(joined_table[col + '1'], joined_table[col + '2']) / 2
            new_col = astropy.table.column.Column(data=np.array(average), name=col)
            del joined_table[col + '1']
            del joined_table[col + '2']
            joined_table.add_column(new_col)
        else:
            if 'e_' not in col:
                e_col = 'e_' + col

                joined_table[e_col + '2'].filled(np.nan)

                compare_1 = joined_table[e_col + '1'].data
                compare_2 = joined_table[e_col + '2'].data
                compare_1.mask = np.ma.nomask
                compare_2 = compare_2.filled(np.nan)
                compare_small = compare_1 < compare_2
                b_nan = np.isnan(compare_2)
                table_1_bool = np.array(np.logical_xor(compare_small, b_nan))
                table_1_rows = np.where(table_1_bool)[0]
                table_1 = joined_table[table_1_rows]

                table_1.rename_column(col + '1', col)
                table_1.rename_column(e_col + '1', e_col)
                del table_1[col + '2']
                del table_1[e_col + '2']

                table_2_rows = np.where(np.invert(table_1_bool))[0]
                table_2 = joined_table[table_2_rows]
                table_2.rename_column(col + '2', col)
                table_2.rename_column(e_col + '2', e_col)
                del table_2[col + '1']
                del table_2[e_col + '1']

                joined_table = astropy.table.vstack([table_1, table_2], join_type='exact')

    return joined_table


def astropy_table_to_result(table, filters):
    result = []
    filter_mags = [filt + 'mag' for filt in filters]
    mag_columns = []
    for col in table.columns:
        if col in filter_mags:
            mag_columns.append(col)
    table = table[reduce(operator.or_, [~table[col].mask for col in mag_columns])]
    # print(table.info)
    for row in table:
        result_row = dict(id=str(row['id']), isValid=True)
        for filt in filters:
            result_row[filt + 'Mag'] = float(row[filt + 'mag'])
            result_row[filt + 'err'] = float(row['e_' + filt + 'mag'])
            result_row[filt + 'ra'] = float(row['ra'])
            result_row[filt + 'dec'] = float(row['dec'])
            result_row[filt + 'pmra'] = float(row['pmra'])
            result_row[filt + 'pmdec'] = float(row['pmdec'])
            result_row[filt + 'dist'] = float(row['Dist'])
        result.append(result_row)
    return {'data': result, 'filters': filters}


def coordinates_to_dist(ra: float, dec: float, r: float):
    return {'ra': ra, 'dec': dec, 'r': r}


def scraper_query_gaia_esac(coordinates, constrain):
    columns = ['ra', 'dec', 'pmra', 'pmra_error', 'pmdec', 'pmdec_error', 'parallax', 'parallax_error']

    ra = coordinates['ra']
    dec = coordinates['dec']
    r = u.Quantity(coordinates['r'], u.deg)
    query_coords = coord.SkyCoord(ra=ra, dec=dec, unit=(u.deg, u.deg), frame='icrs')

    astroquery.gaia.Gaia.MAIN_GAIA_TABLE = "gaiadr3.gaia_source"
    astroquery.gaia.Gaia.ROW_LIMIT = max_row_limit
    gaia_table = astroquery.gaia.Gaia.cone_search_async(coordinate=query_coords, radius=r, columns=columns)
    gaia_table = gaia_table.get_data()

    if len(gaia_table) == max_row_limit:
        raise Exception('Radius too big')

    mask_prl = np.logical_and(gaia_table['parallax'] > 0,
                              gaia_table['parallax_error'] < 1)
    mask_pmra = gaia_table['pmra_error'] / gaia_table['pmra'] < 0.000833
    mask_pmdec = gaia_table['pmdec_error'] / gaia_table['pmdec'] < 0.000833
    mask = np.logical_and(mask_prl, mask_pmra, mask_pmdec)
    gaia_table = gaia_table[mask]

    dist_col = astropy.table.column.Column(data=np.array(1000 / gaia_table['parallax']), name='Dist')
    gaia_table.add_column(dist_col)

    if constrain['distance']['min']:
        mask_dist = np.logical_and(constrain['distance']['min'] < gaia_table['Dist'],
                                   gaia_table['Dist'] < constrain['distance']['max'])
        gaia_table = gaia_table[mask_dist]
    if constrain['pmra']['min']:
        mask_pmra = np.logical_and(constrain['pmra']['min'] < gaia_table['pmra'],
                                   gaia_table['pmra'] < constrain['pmra']['max'])
        gaia_table = gaia_table[mask_pmra]
    if constrain['pmdec']['min']:
        mask_pmdec = np.logical_and(constrain['pmdec']['min'] < gaia_table['pmdec'],
                                    gaia_table['pmdec'] < constrain['pmdec']['max'])
        gaia_table = gaia_table[mask_pmdec]

    row_len = len(gaia_table['ra'])
    index_col = astropy.table.column.Column(data=np.array(range(1, row_len + 1)), name='id')
    gaia_table.add_column(index_col, 0)

    gaia_coordinates = np.dstack((np.array(gaia_table['ra']), np.array(gaia_table['dec'])))[0]

    return {'gaia_table': gaia_table, 'np_coordinates': gaia_coordinates}


def scraper_query_vizier(coordinates, columns, catalog_vizier):
    ra = coordinates['ra']
    dec = coordinates['dec']
    r = coordinates['r']
    query_coords = coord.SkyCoord(ra=ra, dec=dec, unit=(u.deg, u.deg), frame='icrs')
    vquery = Vizier(columns=columns, row_limit=max_row_limit)
    query = vquery.query_region(query_coords,
                                radius=Angle(r * u.deg),
                                catalog=catalog_vizier,
                                )[0]
    # print(len(query))
    if len(query) == max_row_limit:
        raise Exception('Radius too big')
    # print(query.info)
    return query


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
