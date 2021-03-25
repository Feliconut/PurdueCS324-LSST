# %%
from Supernova import database as db
from Supernova.collection import Collection

from os import listdir
from os.path import join

import pandas as pd
# %%
TNS_SEARCH = 'src/tns_search'

# %% sync solution
db.io.alerts_off()


def add_locus_to_collection(c, locus_id):
    try:
        assert isinstance(locus_id, str)
        assert locus_id.startswith('ZTF')
        locus = db.fetch_locus(locus_id, debug=True)
        c.add(locus)
    except (KeyError, AssertionError) as e:
        print('failed ', locus_id, str(e))


def download_one(type_name):
    type_path = join(TNS_SEARCH, type_name)
    print('work on ', type_name)
    with Collection(type_name) as c:
        for file_name in listdir(type_path):
            fpath = join(type_path, file_name)
            print('work on', fpath)
            df = pd.read_csv(fpath)
            df.columns = list(map(lambda x: str(x).strip(), df.columns))
            for locus_id in df['Disc. Internal Name'].apply(
                    lambda x: str(x).strip()):
                add_locus_to_collection(c, locus_id)


def download_all():
    for type_name in listdir(TNS_SEARCH):
        download_one(type_name)


# %%

if __name__ == '__main__':
    download_all()

# %%