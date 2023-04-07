import astropy.table
import astroquery.gaia
import numpy as np


def scraper_query_gaia_esac(coordinates, constrain):
    catalog = "gaiadr3.gaia_source"
    min_parallax = 1000 / float(constrain['distance']['max'])
    max_parallax = 1000 / float(constrain['distance']['min'])
    print(min_parallax, max_parallax)
    min_pmra = float(constrain['pmra']['min'])
    max_pmra = float(constrain['pmra']['max'])
    min_pmdec = float(constrain['pmdec']['min'])
    max_pmdec = float(constrain['pmdec']['max'])
    columns = "ra, dec, pmra, pmra_error, pmdec, pmdec_error, parallax, parallax_error, " \
              "phot_g_mean_mag, phot_bp_mean_mag, phot_rp_mean_mag, " \
              "phot_g_mean_flux_over_error, phot_bp_mean_flux_over_error, phot_rp_mean_flux_over_error"

    query = "SELECT {columns} FROM {catalog} WHERE 1=CONTAINS(POINT({ra}, {dec}),CIRCLE(ra, dec, {radius})) " \
            "AND pmra BETWEEN {min_pmra} AND {max_pmra} " \
            "AND pmdec BETWEEN {min_pmdec} AND {max_pmdec} " \
            "AND parallax BETWEEN {min_parallax} AND {max_parallax} " \
            "AND parallax_error < 1".format(
                columns=columns, catalog=catalog, ra=coordinates['ra'], dec=coordinates['dec'], radius=coordinates['r'],
                min_pmra=min_pmra, max_pmra=max_pmra, min_pmdec=min_pmdec, max_pmdec=max_pmdec,
                min_parallax=min_parallax, max_parallax=max_parallax)
    result = astroquery.gaia.Gaia.launch_job(query).results
    # create index column
    row_len = len(result['ra'])
    index_col = astropy.table.column.Column(data=np.array(range(1, row_len + 1)), name='id')
    result.add_column(index_col, 0)
    # create distance column
    dist_col = astropy.table.column.Column(data=np.array(1000 / result['parallax']), name='Dist')
    del result['parallax']
    del result['parallax_error']
    result.add_column(dist_col)
    # calculate errors
    g_err_col = astropy.table.column.Column(
        data=np.array(result['phot_g_mean_mag']/result['phot_g_mean_flux_over_error']), name='e_Gmag')
    del result['phot_g_mean_flux_over_error']
    result.add_column(g_err_col)
    bp_err_col = astropy.table.column.Column(
        data=np.array(result['phot_bp_mean_mag']/result['phot_bp_mean_flux_over_error']), name='e_BPmag')
    del result['phot_bp_mean_flux_over_error']
    result.add_column(bp_err_col)
    rp_err_col = astropy.table.column.Column(
        data=np.array(result['phot_rp_mean_mag']/result['phot_rp_mean_flux_over_error']), name='e_RPmag')
    del result['phot_rp_mean_flux_over_error']
    result.add_column(rp_err_col)
    # rename phot columns
    result.rename_column('phot_g_mean_mag', 'Gmag')
    result.rename_column('phot_bp_mean_mag', 'BPmag')
    result.rename_column('phot_rp_mean_mag', 'RPmag')
    return result


# {'distance': {'min': -42050, 'max': 44950}, 'pmra': {'min': -4.5809999999999995,'max': -3.315},'pmdec': {'min':
# -0.124,'max': 0.866}} ra = '115.45' dec = '-14.80' r = '3'

result = scraper_query_gaia_esac(
    {'ra': '115.45', 'dec': '-14.80', 'r': '0.1'},
    {'distance': {'min': 1410, 'max': 1480},
     'pmra': {'min': -4.5809999999999995, 'max': -3.315},
     'pmdec': {'min': -0.124, 'max': 0.866}}
)

print(result.info)
result.show_in_browser()
