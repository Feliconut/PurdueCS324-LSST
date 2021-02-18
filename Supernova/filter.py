# %%
from antares_client._api.models import Locus

# %%
def range(locus: Locus) -> bool:
    lc = locus.lightcurve
    _julian_dates = lc['ant_mjd']
    date_range = _julian_dates.max() - _julian_dates.min()

    return date_range < 200