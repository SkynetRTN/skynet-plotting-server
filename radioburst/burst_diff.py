import astropy.table as pytable
from astropy.io import ascii
import numpy as np
import grispy as gsp
import math

from cluster_pro_scraper import scraper_query_gaia

photometry = pytable.Table.read('photometry.csv')
photometry.rename_column('dec_degs', 'dec')
photometry.add_column(pytable.column.Column(data=np.array(photometry['ra_hours']) * 15, name='ra'))
del photometry['id']

print(photometry[0])
print(photometry.info)


def query_range(ras, decs):
    max_ra = np.max(ras)
    min_ra = np.min(ras)
    max_dec = np.max(decs)
    min_dec = np.min(decs)

    delta_ra = max_ra - min_ra
    wrap = delta_ra > 180

    if wrap and min_ra - max_dec < 180:
        return [0, 90, 90 - abs(min_dec)]

    center_ra = (min_ra + max_ra) / 2
    # if min_ra - max_dec > 180:
    #     if center_ra > 180:
    #         center_ra -= 180
    #     else:
    #         center_ra += 180
    center_dec = (max_dec + min_dec) / 2

    return {'ra': center_ra, 'dec': center_dec, 'r': 0.4, 'wrap': wrap}



coord = query_range(photometry['ra'], photometry['dec'])

gaia = scraper_query_gaia(coord, None)
grid = gsp.GriSPy(gaia['np_coordinates'], N_cells=128)
result_table = gaia['gaia_table']

target_cord = np.dstack((np.array(photometry['ra']), np.array(photometry['dec'])))[0]

nn_dist, nn_indices = grid.nearest_neighbors(target_cord, n=1)

print(gaia['gaia_table'].info)

nn_dist = np.concatenate(nn_dist)
nn_indices = np.concatenate(nn_indices)
nn_indices_filtered = [nn_indices[i] if nn_dist[i] < 0.000833 else 0 for i in range(0, len(nn_indices))]

nn_indices_col = pytable.column.Column(data=np.array(nn_indices_filtered), name='id')
photometry.add_column(nn_indices_col, 0)

joined_table = pytable.join(left=result_table, right=photometry, keys='id', join_type='right')

# joined_table.show_in_browser()

ascii.write(joined_table, 'right_join.csv', format='csv', overwrite=True, fast_writer=False)
