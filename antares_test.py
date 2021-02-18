# %%
from antares_client.search import get_by_id, get_by_ztf_object_id

# Lookup by ANTARES ID
locus = get_by_id("ANT2020j7wo4")

# Lookup by ZTF Object ID
locus = get_by_ztf_object_id("ZTF20aafqubg")
# %%
from antares_client.search import cone_search
from astropy.coordinates import Angle, SkyCoord

center = SkyCoord("20h48m25.1805s 29d45m4.8361s")
radius = Angle("1s")

for locus in cone_search(center, radius):
    pass
# %%
# query
query = {
    "query": {
        "bool": {
            "filter": [
                {
                    "range": {
                        "properties.num_mag_values": {
                            "gte": 20,
                            # "lte": 100,
                        }
                    }
                },
                {
                    "term": {
                        # signal-noise ratio cut.
                        "tags": "refitt_newsources_snrcut"
                    }
                }
            ]
        }
    }
}
# %%
from antares_client.search import search


# Simple Filter
def filter_search(query, *criteria):
    i = 0
    sr = search(query)
    for locus in sr:
        i += 1
        for crit in criteria:
            if crit(locus):
                print(f'passed {i} loci')
                yield locus
                i = 0


def filter_1(locus):
    lc = locus.lightcurve
    _julian_dates = lc['ant_mjd']
    date_range = _julian_dates.max() - _julian_dates.min()

    return date_range < 200


# %%
filtered_search = filter_search(query, filter_1)
# %%
# treating a new filtered locus
locus = next(filtered_search)

# %%
lc = locus.lightcurve

# two types of alert_id: upper_limit and candidate.
_ulim_and_candidate = lc.__len__()
_candidate = lc[lc['alert_id'].apply(
    lambda x: x.startswith('ztf_candidate'))].__len__()
lc['alert_type'] = lc['alert_id'].apply(lambda x: x[:x.index(':')])

_julian_dates = lc['ant_mjd']
date_range = _julian_dates.max() - _julian_dates.min()

print(f'We have {_ulim_and_candidate - _candidate} limit alerts' +
      f'and {_candidate} candidate alerts.')
print(f'This event ranges {date_range:.2f} days.')
# %%
# Visualize a lightcurve
import matplotlib.pyplot as plt

fig = plt.figure(num=None, figsize=(12, 5), dpi=100)
for (pb, at), df in lc.groupby(['ant_passband', 'alert_type']):
    is_candidate = at == 'ztf_candidate'
    plt.errorbar(
        x=df['ant_mjd'],
        y=df['ant_mag'] if is_candidate else df['ant_maglim'],
        yerr=df['ant_magerr'],
        #  uplims=(at!='ztf_candidate'),
        label=pb + '  ' + at[4:],
        color=pb,
        fmt=('o' if is_candidate else '^') + '--' + pb.lower(),
        alpha=1 if is_candidate else 0.3)
plt.title(locus.properties['ztf_object_id'])

plt.legend()
plt.show()
# %%
