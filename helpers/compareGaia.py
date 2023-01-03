import astroquery.gaia
import numpy as np
import astropy
import astropy.units as u
import astropy.coordinates as coord
import grispy as gsp

from cluster_pro_scraper import scraper_query_vizier
from gaia_util import gaia_get_data
from astropy.io import ascii


def scraper_query_gaia(coordinates, constrain):
    columns = ['RA_ICRS', 'DE_ICRS', 'pmRA', 'e_pmRA', 'pmDE', 'e_pmDE', 'Dist', 'Plx']

    gaia_table = scraper_query_vizier(coordinates, columns, 'I/355/gaiadr3', constrain)

    gaia_coordinates = np.dstack((np.array(gaia_table['RA_ICRS']), np.array(gaia_table['DE_ICRS'])))[0]

    return {'gaia_table': gaia_table, 'np_coordinates': gaia_coordinates}


def scraper_query_gaia_esac(coordinates, constrain):
    columns = ['ra', 'dec', 'pmra', 'pmra_error', 'pmdec', 'pmdec_error', 'parallax', 'parallax_error']

    ra = coordinates['ra']
    dec = coordinates['dec']
    r = u.Quantity(coordinates['r'], u.deg)
    query_coords = coord.SkyCoord(ra=ra, dec=dec, unit=(u.deg, u.deg), frame='icrs')

    astroquery.gaia.Gaia.MAIN_GAIA_TABLE = "gaiadr3.gaia_source"
    astroquery.gaia.Gaia.ROW_LIMIT = 50000
    gaia_table = astroquery.gaia.Gaia.cone_search_async(coordinate=query_coords, radius=r, columns=columns)
    gaia_table = gaia_table.get_data()

    # row_len = len(gaia_table['ra'])
    # index_col = astropy.table.column.Column(data=np.array(range(1, row_len + 1)), name='id')
    # gaia_table.add_column(index_col, 0)

    gaia_coordinates = np.dstack((np.array(gaia_table['ra']), np.array(gaia_table['dec'])))[0]

    return {'gaia_table': gaia_table, 'np_coordinates': gaia_coordinates}


ngc_2437_coord = {"ra": 115.44, "dec": -14.80, "r": 0.078}

esac = scraper_query_gaia_esac(ngc_2437_coord, None)
esac_table = esac['gaia_table']
esac_coord = esac['np_coordinates']

vizier = scraper_query_gaia(ngc_2437_coord, None)
vizier_table = vizier['gaia_table']
vizier_coord = vizier['np_coordinates']

local = gaia_get_data(ngc_2437_coord)

local = np.array(local)
local = astropy.table.Table(local, names=('source_id', 'ra_loc', 'dec_loc', 'dist_loc', 'pmra_loc', 'pmdec_loc'))
del local['source_id']
row_len = len(local['ra_loc'])
index_col = astropy.table.column.Column(data=np.array(range(1, row_len + 1)), name='id')
local.add_column(index_col, 0)

local_cord = np.dstack((np.array(local['ra_loc']), np.array(local['dec_loc'])))[0]

grid = gsp.GriSPy(local_cord, N_cells=128)

nn_dist, nn_indices = grid.nearest_neighbors(esac_coord, n=1)
nn_dist = np.concatenate(nn_dist)
nn_indices = np.concatenate(nn_indices)
nn_indices_filtered = [nn_indices[i]+1 if nn_dist[i] < 0.000833 else 0 for i in range(0, len(nn_indices))]

nn_indices_col = astropy.table.column.Column(data=np.array(nn_indices_filtered), name='id')
esac_table.add_column(nn_indices_col, 0)

joined_table = astropy.table.join(left=local, right=esac_table, keys='id', join_type='outer')
parl = joined_table['parallax']
parl_dist = 1000/parl
parl_dist = astropy.table.column.Column(data=np.array(parl_dist), name='dist_esac')
joined_table.add_column(parl_dist)

nn_dist, nn_indices = grid.nearest_neighbors(vizier_coord, n=1)
nn_dist = np.concatenate(nn_dist)
nn_indices = np.concatenate(nn_indices)
nn_indices_filtered = [nn_indices[i]+1 if nn_dist[i] < 0.000833 else 0 for i in range(0, len(nn_indices))]

nn_indices_col = astropy.table.column.Column(data=np.array(nn_indices_filtered), name='id')
vizier_table.add_column(nn_indices_col, 0)

joined_table = astropy.table.join(left=joined_table, right=vizier_table, keys='id', join_type='outer')

parl = joined_table['Plx']
parl_dist = 1000/parl
parl_dist = astropy.table.column.Column(data=np.array(parl_dist), name='dist_vizier')
joined_table.add_column(parl_dist)




# joined_table.show_in_browser()
# ascii.write(
#     joined_table, "compareGaia.csv", format="csv", fast_writer=False, overwrite=True
# )
