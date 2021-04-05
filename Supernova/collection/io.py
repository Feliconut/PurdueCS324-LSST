from Supernova.database.io import fetch_locus
from typing import Iterable
from antares_client._api.models import Locus
from os.path import join
from os import listdir, remove

COLLECTIONS_PATH = join('src', 'collections')


class Collection():
    def __init__(self, name, locus_ids=None, doc='') -> None:
        self.name = name
        self.locus_ids = set(locus_ids) if locus_ids else set()
        self.doc = doc

    def write(self):
        fpath = join(COLLECTIONS_PATH, self.name)
        with open(fpath, 'w') as f:
            f.write(self.__repr__())

    def delete(self):
        remove_collection(self.name)

    def add(self, locus):
        if isinstance(locus, Locus):
            self.locus_ids.add(locus.locus_id)
        elif isinstance(locus, str):
            self.locus_ids.add(locus)
        else:
            raise ValueError(locus)

    def remove(self, locus):
        if isinstance(locus, Locus):
            self.locus_ids.remove(locus.locus_id)
        elif isinstance(locus, str):
            self.locus_ids.remove(locus)
        else:
            raise ValueError(locus)

    def extend(self, stream: Iterable[Locus]):
        if isinstance(stream, Collection):
            self.locus_ids.update(stream.locus_ids)
        else:
            for locus in stream:
                self.add(locus)
                self.write()

    def __enter__(self):
        # self._load()
        return self

    def __exit__(self, exc_type=None, exc_val=None, exc_tb=None):
        # self.flush()
        self.write()

    def __len__(self):
        return len(self.locus_ids)

    def __iter__(self):
        from ..database.io import fetch_locus
        for name in self.locus_ids:
            try:
                yield fetch_locus(name)
            except KeyError:
                print(f'Locus {name} skipped due to fail of fetching.')

    def __repr__(self):
        return (
        f'Collection('+\
        f'name={repr(self.name)},'+\
        f'doc={repr(self.doc)},'+\
        f'locus_ids={repr(self.locus_ids)},'+\
        ')'.replace('\n', ''))

    def pack(self):
        return pack(self)

    def __getitem__(self, index):
        for i in list(self.locus_ids)[index]:
            yield fetch_locus(i)

    def __contains__(self, locus: Locus):
        return locus.locus_id in self.locus_ids


def fetch_collection(name):
    try:
        with open(join(COLLECTIONS_PATH, name), 'r') as f:
            return eval(f.read())
    except FileNotFoundError:
        return Collection(name)
    except SyntaxError:
        print(f'local storage invalid. remove and retry.')
        remove_collection(name)


def save_collection(collection: Collection):
    collection.write()


def remove_collection(name):
    try:
        fpath = join(COLLECTIONS_PATH, name)
        remove(fpath)
    except FileNotFoundError:
        print(name, ' is alread removed, or does not exist.')


def pack(collection: Collection):
    """Pack up all locus data in the collection in a .npz file"""
    from ..database.pack import pack
    pack(collection.name, collection.locus_ids)