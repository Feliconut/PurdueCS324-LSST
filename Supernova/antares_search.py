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


def apply_filters(stream: Generator[Locus, None, None], *criteria):
    def validate(locus, *criteria):
        i = 0
        for crit in criteria:
            if not crit(locus):
                i += 1
                print(
                    f'locus {locus.locus_id} filtered out by {crit.__name__}')
                return False
        return True

    for locus in stream:
        if validate(locus):
            yield locus


# backward compatibility
def filter_search(query, *criteria):
    sr = search(query)
    return apply_filters(sr, *criteria)
