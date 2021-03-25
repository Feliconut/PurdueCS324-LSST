"""
Local Searching Functionalities.
"""
from os import listdir
from os.path import join
from .io import DATA_PATH, fetch_lightcurve, fetch_locus


def all_loci():
    """
    Iterates Through all local locus.
    
    @return Generator[Locus, Any, Any]
    """
    for fname in listdir(join(DATA_PATH, 'loci')):
        try:
            yield fetch_locus(fname)
        except Exception as e:
            print(f'{repr(e)} fetching {fname}')


def all_lightcurves():
    """
    Iterates Through all local lightcurves.
    
    @return Generator[pd.Dataframe, Any, Any]
    """
    for fname in listdir(join(DATA_PATH, 'lightcurves')):
        try:
            locus_id, fmt = fname.split('.')
            assert fmt == 'lc'
            yield fetch_lightcurve(locus_id)
        except Exception as e:
            print(f'{repr(e)} fetching {fname}')


def make_antares_ztf_dict():
    result = dict()
    for locus in all_loci():
        result.setdefault(locus.locus_id, locus.properties['ztf_object_id'])
    return result