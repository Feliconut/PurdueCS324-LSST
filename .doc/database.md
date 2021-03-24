# Locus Database

**What can it do?** 

It keeps a whole bunch of `Locus` objects together with their lightcures in the local storage. You can save and load `Locus` and `lightcurve`, and search the database.

**Why it's important?**

We store locus data in local machine to avoid time and resource limitation of Antares.​
Querying each `Locus` online will be unworkable if we deal with a large amount of lightcurve data.

## I/O

Implemented in `database/io.py`

### API

We provide the following methods. Their usage should be apparant from naming. If need clarification, check out the source code.
1. `fetch_alert(alert_id) -> alert: Alert`
1. `fetch_locus(locus_id) -> locus: Locus`
1. `fetch_lightcurve(lightcurve_id) -> lightcurve: pd.DataFrame`
1. `add_alert(alert: Alert)`
1. `add_locus(locus: Locus)`
1. `add_lightcurve(locus: Locus)`

Since `Alert`s are numerous and usually not useful, not saving them can speed up the process. Use this command to supress alert saving:

`database.io.alert_off()`

Similarly, use this to turn on alert saving:

`database.io.alert_on()`

### The Saving and Loading Scheme

Each Locus object is serialized to a python expression that, upon evaluated, will produce the exactly same object.

Each alert is serialized similarly.

Each lightcurve is stored using `Dataframe.to_feather()` and `pd.read_feather()`. This is an implementation of Apache Arrow.

## Iteration and Querying

implemented in `database/search.py`.

`all_loci()` and `all_lightcurves()` are generators. In particular, `all_loci()` can be passed on to a filter. See `filter.md` for more information.
