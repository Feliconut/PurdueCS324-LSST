# %%
from Supernova.antares_search import query, filter_search
from Supernova import filter
from Supernova.visualize import info, plt_lightcurve
# %%
filtered_search = filter_search(
    query,
    filter.date_range, 
    debug=True)
# %%
# treating a new filtered locus
locus = next(filtered_search)
# %%
info(locus)
# %%
plt_lightcurve(locus)
# %%
from Supernova import database as db
locus_id = locus.locus_id
db.add_locus(locus)

loaded_locus = db.fetch_locus(locus)
# %%
