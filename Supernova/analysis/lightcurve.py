# %%
from enum import Enum, Flag, auto

from pandas.core.frame import DataFrame


class Band(Enum):
    R = 'R'
    G = 'g'


class AlertType(Flag):
    cand = auto()
    ulim = auto()
    llim = auto()


# %%
# extract a band.
def extract_info(lightcurve: DataFrame, band: Band, alertType: AlertType):
    '''
    @return A series of brightness value in the given band.
    Include jld? y
    '''
    lc = lightcurve

    lc =  (
        lc[(\
            lc['alert_id'].apply(lambda x: x.startswith('ztf_candidate'))&\
            lc['ant_passband'].apply(lambda x: x == band.value))][
                [
                    'ant_mjd', 'ant_mag', 'ant_magerr'
                ]
            ])

    lc['ant_mjd'] -= lc['ant_mjd'].min()

    lc = lc.reset_index(drop=True)
    lc = lc.rename(columns={
        'ant_mjd': 'mjd',
        'ant_mag': 'mag',
        'ant_magerr': 'err'
    })
    return lc


# %%
def get_date_range(lc):
    _julian_dates = lc['ant_mjd']
    duration = _julian_dates.max() - _julian_dates.min()
    return duration