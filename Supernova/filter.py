'''
Filters that determine whether a Locus meets a particular criteria.
'''
# %%
from antares_client._api.models import Locus
from typing import Generator


def apply(
    stream: Generator[Locus, None, None],
    *filters,
    debug=False,
):
    '''
    Apply a series of filters to an iterable sequence of loci.
    
    @param stream Input sequence of loci.
    
    @param filters A sequence of filters you want to apply.
    If no filters is provided, then this function has no effect.

    @param debug If set True, then it will print out when each Locus gets filtered out.
    '''
    if not debug:

        def validate(locus):
            for crit in filters:
                if not crit(locus):
                    # if debug:
                    #     print(
                    #         f'locus {locus.locus_id} filtered out by {crit.__name__}'
                    #     )
                    return False
            return True

        for locus in stream:
            if validate(locus):
                yield locus
    else:

        def validate(locus):
            for f in filters:
                if not f(locus):
                    counters[f] += 1
                    # if debug:
                    #     print(
                    #         f'locus {locus.locus_id} filtered out by {crit.__name__}'
                    #     )
                    return False
            counters['pass'] += 1
            return True

        counters = dict()
        for f in filters:
            counters.setdefault(f, 0)
            counters.setdefault('pass', 0)
        for locus in stream:
            if validate(locus):
                yield locus

        print(('->'.join((f'{f.__name__}({counters[f]})' for f in filters))),
              '->', f'pass({counters["pass"]})')

        return [counters[f] for f in filters] + [counters['pass']]


# %% Meta-filters


def logic_and(*fs):
    def composite_filter(locus: Locus) -> bool:
        for f in fs:
            if not f(locus):
                return False
        return True

    return composite_filter


def logic_or(*fs):
    def composite_filter(locus: Locus) -> bool:
        for f in fs:
            if f(locus):
                return True
        return False

    return composite_filter


def logic_not(f):
    def composite_filter(locus: Locus) -> bool:
        return not f(locus)

    return composite_filter


# %%
# Standard Filters.


def date_range(max=200, min=0):
    from .analysis.lightcurve import get_date_range

    def date_range_filter(locus: Locus) -> bool:
        duration = get_date_range(locus.lightcurve)
        return min <= duration <= max

    return date_range_filter


def lightcurve_datapoints(min=20, max=9e9):
    def lc_datapoints_filter(locus: Locus) -> bool:
        lc = locus.lightcurve
        return min <= len(lc) <= max

    return lc_datapoints_filter


def snr_slope_cut(snr_max=0.3, slope_max=0.0025):
    "True if likely non-bogus. "
    from .analysis.simple_classification import snr_and_linear

    def snr_slope_filter(locus: Locus) -> bool:
        snr, slope = snr_and_linear(locus)
        return snr > snr_max or slope > slope_max

    return snr_slope_filter
