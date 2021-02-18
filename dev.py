# %%
from Supernova.antares_search import query, filter_search
from Supernova import filter
from Supernova.visualize import info, plt_lightcurve
# %%
filtered_search = filter_search(query, filter.range)
# %%
# treating a new filtered locus
locus = next(filtered_search)
# %%
info(locus)
# %%
plt_lightcurve(locus)
# %%
