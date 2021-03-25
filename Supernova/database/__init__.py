__all__ = ['add_locus', 'fetch_locus', 'io', 'search']

from os import mkdir as __mkdir

print('Supernovae Database Initialization Start')


def __new_folder(path):
    try:
        __mkdir(path)
        print('created path ' + path)
        return True
    except FileExistsError:
        print(f'data folder {path} already exists')
        return False


__new_folder('DATA')
__new_folder('DATA/alerts')
__new_folder('DATA/lightcurves')
__new_folder('DATA/loci')

with open('DATA/.gitignore', 'w+') as f:
    f.write('*')

print('Supernovae Database Initialization Successful')

from .io import add_locus, fetch_locus, fetch_lightcurve, fetch_alert, reindex as __reindex
from .search import all_lightcurves, all_loci
__reindex()
