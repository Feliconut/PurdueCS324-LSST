# Supernova

## Locus and Lightcurve

This paragraph assumes knowledge of what is a locus and lightcurve.

## Locus Querying

Access locus data from a generator stream.

### Antares Query Stream

See annotation for `antares_search.py`

Sample Code
```python
from Supernova.antares_search import default_search
default_search() # this returns a Generator[Locus]

# Or you can specify the query
from Supernova.antares_search import search

query = {
    "query": {
        "bool": {
            "filter": [
                {
                    "range": {
                        "properties.num_mag_values": {
                            "gte": 20,
                            "lte": 100,
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
search(query) # this returns a Generator[Locus]
```

### TNS Query Stream

Not Implemented Yet, because we don't have API key now.

### Filtering the Query

see `filter.md`

## Local Locus Database

see `database.md`

## Locus Analysis

TBD