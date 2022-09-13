from astroquery.simbad import Simbad


def scraper_query_object(query: str):
    simbad = Simbad()
    simbad.add_votable_fields('bibcodelist(2003-2013)')
    raw = simbad.query_object(query)
    # ['MAIN_ID', 'RA', 'DEC', 'RA_PREC', 'DEC_PREC',
    # 'COO_ERR_MAJA', 'COO_ERR_MINA', 'COO_ERR_ANGLE', 'COO_QUAL', 'COO_WAVELENGTH', 'COO_BIBCODE',
    # 'BIBLIST_2003_2013', 'SCRIPT_NUMBER_ID']
    # raw.pprint()
    return {'RA': hms2d(raw['RA'][0]), 'DEC': dms2d(raw['DEC'][0]), 'Range': 0}


def hms2d(hms: str):
    hms = hms.split(' ')
    result = float(hms[0])*15+float(hms[1])*0.25
    if len(hms) == 3:
        result += float(hms[2])*0.25/60
    return result

def dms2d(dms: str):
    dms = dms.split(' ')
    result = float(dms[0])
    is_positive = True
    if result < 0:
        result = -result
        is_positive = False
    result += float(dms[1])/60
    if len(dms) == 3:
        result += float(dms[2])/3600
    if is_positive:
        return result
    else:
        return -result
