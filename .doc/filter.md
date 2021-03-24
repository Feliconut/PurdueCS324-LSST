# Filter System

## What is a filter

A **filter** is essentially a **function** that takes in a `Locus` and returns a `Boolean`.
`True` means pass the filter, `False` means filtered out by the filter.

### What does it act on

We **apply filters** to a **stream** of Loci, i.e. a `Generator[Locus]`.

### What does it produce

Applying a filter to a Loci stream gives us a `Generator[Locus]` which only yields filtered Loci.

## How to use a filter
```python
from Supernova.antares_search import default_search
from Supernova import filter

locus = next(default_search())

filter.date_range(100)(locus) # This returns a boolean.
```

Note that the function `filter.date_range` is a **meta-function** that produces a **filter**(function) with specified parameters. In this case, the input is the maximum numbers of julian days of this locus event.

### Where to Find Filters

All standard filters are kept in `filter.py`. These filters should explain themselves clearly in their docstring so we don't provide doc here.

### How to apply filters to a search

Use `filter.apply` method.

```python
from Supernova.antares_search import default_search
from Supernova import filter

filtered_search = filter.apply(
    default_search(), # The input Loci Stream
    filter.date_range(200), # All Filters to be applied
    filter.lightcurve_datapoints(50)
    )
```

## Creating a Filter

Just mimick the existing filters in order to create yours. Doc for this part can be added if you request.
