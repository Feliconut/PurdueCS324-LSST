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