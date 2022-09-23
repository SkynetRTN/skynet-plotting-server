from astroquery.simbad import Simbad
from astroquery.vizier import Vizier
from astropy.coordinates import Angle
import astropy
import astropy.units as u
import astropy.coordinates as coord


def scraper_query_object(query: str):
    simbad = Simbad()
    simbad.add_votable_fields('bibcodelist(2003-2013)')
    raw = simbad.query_object(query)
    # ['MAIN_ID', 'RA', 'DEC', 'RA_PREC', 'DEC_PREC',
    # 'COO_ERR_MAJA', 'COO_ERR_MINA', 'COO_ERR_ANGLE', 'COO_QUAL', 'COO_WAVELENGTH', 'COO_BIBCODE',
    # 'BIBLIST_2003_2013', 'SCRIPT_NUMBER_ID']
    # raw.pprint()
    return {'RA': hms2d(raw['RA'][0]), 'DEC': dms2d(raw['DEC'][0]), 'Range': 0}


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
                           G=float(row['Gmag']),
                           Gerr=float(row['e_FG']),
                           Gra=float(row['RA_ICRS']),
                           Gdec=float(row['DE_ICRS']),
                           GBP=float(row['BPmag']),
                           GBPerr=float(row['e_FBP']),
                           GBPra=float(row['RA_ICRS']),
                           GBPdec=float(row['DE_ICRS']),
                           GRP=float(row['RPmag']),
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
