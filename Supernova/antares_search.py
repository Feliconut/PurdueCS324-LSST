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
from typing import Generator
from antares_client._api.models import Locus
from antares_client.search import search
from marshmallow.fields import Boolean, Function, List

def default_search():
    return search(query)

def apply_filters(stream:Generator[Locus, None, None], *criteria):
    for locus in stream:
        for crit in criteria:
            if not crit(locus): continue
        yield locus