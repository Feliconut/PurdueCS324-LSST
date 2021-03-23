"""
Local Storage of Locus, Lightcurve and Alert.

@Author: Xiaoyu Liu
"""

from os.path import join
from antares_client._api.models import Locus, Alert
from antares_client import search
from pandas.io.feather_format import read_feather

DATA_PATH = './DATA'
keep_alerts = True


def alerts_on():
    global keep_alerts
    keep_alerts = True


def alerts_off():
    global keep_alerts
    keep_alerts = False


def encode_alert(alert: Alert):
    return (
        f'Alert('+\
        f'alert_id={repr(alert.alert_id)},'+\
        f'mjd={repr(alert.mjd)},'+\
        f'properties={repr(alert.properties)}'+\
        ')'.replace('\n', ''))


def encode_locus(locus: Locus):
    if keep_alerts:
        for alert in locus.alerts:
            add_alert(alert)
    add_lightcurve(locus)
    return (
        f'Locus('+\
            (f'alerts=[fetch_alert(alert_id) for alert_id in {repr([alert.alert_id for alert in locus.alerts])}],' \
                if keep_alerts else 'alerts=[],')+\
            f'catalog_objects={repr(locus.catalog_objects)},'+\
            f'catalogs={repr(locus.catalogs)},'+\
            f'dec={repr(locus.dec)},'+\
            f'lightcurve=fetch_lightcurve({repr(locus.locus_id)}),'+\
            f'locus_id={repr(locus.locus_id)},'+\
            f'properties={repr(locus.properties)},'+\
            f'ra={repr(locus.ra)},'+\
            f'tags={repr(locus.tags)},'+\
            # f'timeseries={repr(locus.timeseries)},'+\
            f'watch_list_ids={repr(locus.watch_list_ids)},'+\
            f'watch_object_ids={repr(locus.watch_object_ids)},'+\
        ')'.replace('\n', ''))


def add_locus(locus: Locus):
    with open(join(DATA_PATH, 'loci', locus.locus_id), 'w+') as f:
        f.write(encode_locus(locus))


def add_alert(alert: Alert):
    with open(join(DATA_PATH, 'alerts', alert.alert_id.replace(':', '')),
              'w+') as f:
        f.write(encode_alert(alert))


def add_lightcurve(locus: Locus):
    locus.lightcurve.to_feather(
        join(DATA_PATH, 'lightcurves', locus.locus_id + '.lc'))


def fetch_locus(locus_id):
    try:
        with open(join(DATA_PATH, 'loci', locus_id), 'r') as f:
            return eval(f.read())
    except FileNotFoundError:
        print(f'Locus {locus_id} not found in local. Trying remote.')
        try:
            res = search.get_by_id(locus_id)
            if res:
                add_locus(res)
                return res
            res = search.get_by_ztf_object_id(locus_id)
            if res:
                add_locus(res)
                return res
        except Exception:
            pass
    raise KeyError(f'Locus {locus_id} not found in remote.')


def fetch_alert(alert_id):
    try:
        with open(join(DATA_PATH, 'alerts', alert_id.replace(':', '')),
                  'r') as f:
            return eval(f.read())
    except FileNotFoundError:
        pass
    raise KeyError(f'Alert {alert_id} not found in local.')


def fetch_lightcurve(locus_id):
    return read_feather(join(DATA_PATH, 'lightcurves', locus_id + '.lc'))