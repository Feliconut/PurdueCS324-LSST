# %%
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
def filter_search(query, *criteria, debug=False):
    i = 0
    sr = search(query)
    for locus in sr:
        i += 1
        for crit in criteria:
            if crit(locus):
                if debug: print(f'passed {i} loci')
                yield locus
                i = 0