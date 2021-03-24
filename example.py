"""
This file shows most functions of the package Supernova.
"""

# %%
from Supernova import database as db
from Supernova import filter
from Supernova.antares_search import default_search
from Supernova.visualize import info, plt_lightcurve

# %%
db.io.alerts_off()
filtered_search = filter.apply(
    default_search(), 
    filter.date_range(200),
    filter.lightcurve_datapoints(50)
    )
# %%
# Saving Locus
for i in range(10):
    locus = next(filtered_search)
    db.add_locus(locus)

# %% Fetching Locus
locus_id = locus.locus_id
loaded_locus = db.fetch_locus(locus_id)
info(locus)
# %% Iterating through all local locus
for locus in db.all_loci():
    info(locus)
    plt_lightcurve(locus)
# %%
