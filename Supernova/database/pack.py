from shutil import copy, make_archive, rmtree
from os import mkdir, remove
from os.path import join, exists
from .io import DATA_PATH, fetch_locus


def pack(name, locus_ids, include_alerts=False):
    if exists(name + '.zip'):
        raise FileExistsError(name + '.zip')
    DST_PATH = name + '_temp'
    print(f'Creating temp folder ./{name}_temp ...')
    mkdir(DST_PATH)
    mkdir(join(DST_PATH, 'loci'))
    mkdir(join(DST_PATH, 'lightcurves'))
    mkdir(join(DST_PATH, 'alerts'))
    print(f'Copying necessary files ...')
    for locus_id in locus_ids:
        copy(join(DATA_PATH, 'loci', locus_id), join(DST_PATH, 'loci',
                                                     locus_id))
        copy(join(DATA_PATH, 'lightcurves', locus_id + '.lc'),
             join(DST_PATH, 'lightcurves', locus_id + '.lc'))
        if include_alerts:
            for alert in fetch_locus(locus_id).alerts:
                alert_id = alert.alert_id
                copy(join(DATA_PATH, 'alerts', alert_id),
                     join(DST_PATH, 'alerts', alert_id))
    print(f'Making {name}.zip ...')
    make_archive(name, 'zip', DST_PATH)
    print(f'Complete. Clearing temp files')
    rmtree(DST_PATH)
    print(f'Complete.')
