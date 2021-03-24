'''
Filters that determine whether a Locus meets a particular criteria.
'''
# %%
from antares_client._api.models import Locus
from typing import Generator


def apply(stream: Generator[Locus, None, None], *filters, debug=False):
    '''
    Apply a series of filters to an iterable sequence of loci.
    
    @param stream Input sequence of loci.
    
    @param filters A sequence of filters you want to apply.
    If no filters is provided, then this function has no effect.

    @param debug If set True, then it will print out when each Locus gets filtered out.
    '''
    def validate(locus, *filters):
        for crit in filters:
            if not crit(locus):
                if debug:
                    print(
                        f'locus {locus.locus_id} filtered out by {crit.__name__}'
                    )
                return False
        return True

    for locus in stream:
        if validate(locus):
            yield locus


# %%
# Standard Filters.


def date_range(max=200, min=0):
    def filter(locus: Locus) -> bool:
        lc = locus.lightcurve
        _julian_dates = lc['ant_mjd']
        duration = _julian_dates.max() - _julian_dates.min()
        return min <= duration <= max

    return filter


def lightcurve_datapoints(min=20, max=9e9):
    def filter(locus: Locus) -> bool:
        lc = locus.lightcurve
        return min <= len(lc) <= max

    return filter
