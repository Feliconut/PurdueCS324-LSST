"""
In this sub-project, I get a collection of bogus events used for testing.
"""
# %%
from Supernova.visualize import plt_lightcurve
from Supernova.database.io import fetch_lightcurve, fetch_locus
from antares_client._api.models import Locus
from Supernova.collection.io import fetch_collection
from Supernova.collection import Collection
from Supernova.antares_search import default_search, search
from Supernova.database import add_locus
from Supernova import filter

query = {"query": {"bool": {"filter": []}}}


# %% Bogus from Non-SNRcut stream. This might be biased.
def non_snrcut(locus: Locus) -> bool:
    return 'refitt_newsources_snrcut' not in locus.tags


non_snrcut_feed = filter.apply(
    search(query),
    non_snrcut,
    filter.lightcurve_datapoints(),
)

with fetch_collection('bogus') as c:
    while len(c) < 100:
        locus = next(non_snrcut_feed)
        c.add(locus)
        add_locus(locus)
        fetch_locus(locus.locus_id, debug=True)

# %% Manually extract bogus from antares stream.
feed = filter.apply(
    default_search(),
    filter.lightcurve_datapoints(min=10),
)

with fetch_collection('bogus_manual') as c:
    while len(c) < 40:
        locus = next(feed)
        if locus in c:
            continue
        plt_lightcurve(locus)
        # press N if I beleive this is SN, otherwise just press Enter.
        if input() != 'N':
            c.add(locus)
            add_locus(locus)
            fetch_locus(locus.locus_id, debug=True)

# %%
