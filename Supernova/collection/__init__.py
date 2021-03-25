from os import mkdir as __mkdir

print('Supernovae Collections Initialization Start')


def __new_folder(path):
    try:
        __mkdir(path)
        print('created path ' + path)
        return True
    except FileExistsError:
        print(f'data folder {path} already exists')
        return False


__new_folder('src')
__new_folder('src/collections')
print('Supernovae Database Initialization Successful')

from .io import Collection, fetch_collection, remove_collection