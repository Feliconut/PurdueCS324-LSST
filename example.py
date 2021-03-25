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
# %% Use Collection to define a batch of loci
from Supernova.collection import Collection,fetch_collection

c = Collection('test', doc='it is just a test')
c.add(locus) # Add locus to the collection
c.add(locus.locus_id) # This works also, you can directly add antares id

c.add("ZTF19aauhwqv")
c.remove("ZTF19aauhwqv") # removes a locus

c.write() # Write to local disk

# %% loading from local disk
c_load = fetch_collection('test')

d = Collection('test2')
d.extend(c) # Merging collections

# %% iterating a collection
for locus in c:
    info(locus)

# %% directing a locus stream to a collection

e = Collection('temp')
feed = filter.apply(iter(c), filter.date_range(100))
e.extend(feed)

# %%
c_load.delete() # delete from local
c.delete()
d.delete()
e.delete()
# %%
