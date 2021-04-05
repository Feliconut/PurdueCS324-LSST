"""
In this project, I try to find the best differentiator between SN and NON-SN.

Then I turn this into a filter that produces high-quality stream of unclassified events.
"""

# %%
from antares_client._api.models import Locus
from matplotlib import pyplot as plt
import numpy as np
from functools import partial, reduce

from pandas.core.frame import DataFrame

from Supernova import database as db
from Supernova import filter
from Supernova import collection
from Supernova.antares_search import default_search
from Supernova.collection import Collection, fetch_collection
from Supernova.visualize import plt_lightcurve


def compose(*funcs):
    return lambda x: reduce(lambda acc, f: f(acc), funcs, x)


# %%
# Goal: Find out the best differentiator between SN and NON-SN.
def feature_func(f):
    """
    Evaluate a quantitative measure of a locus.

    @param f The feature function that takes in a single locus and outputs a value.
    @return a function that takes in a collection name and outputs evaluated array of measure.
    """
    return lambda x: np.array(list(map(f, fetch_collection(x)[:40])))


# %% SNR_Comparison
from Supernova.analysis.lightcurve import Band, extract_info


def signal_noise_ratio(a, axis=0, ddof=0):
    a = np.asanyarray(a)
    if len(a) == 0: return 0
    a = np.diff(a)
    m = a.mean(axis)
    sd = a.std(axis=axis, ddof=ddof)
    return np.where(sd == 0, 0, m / sd)


def locus_snr(locus):
    band_mag = lambda band: extract_info(locus.lightcurve, band, None)['mag']

    f = lambda band: signal_noise_ratio(band_mag(band))

    return (f(Band.R)**2 + f(Band.G)**2)**.5


def feature_hist(feature,
                 locusType,
                 color='green',
                 bins=30,
                 alpha=0.3,
                 density=True,
                 clip_range=[0, 1]):
    plt.hist(
        feature_func(feature)(locusType).clip(*clip_range),
        color=color,
        bins=bins,
        alpha=alpha,
        density=density,
        label=locusType,
    )


# %%
SNR_hist = partial(feature_hist, locus_snr, clip_range=[0, 0.5])

def plt_all_types(plot_func):
    plot_func('Ia', color='green')
    plot_func('Ib', color='yellow')
    plot_func('II', color='blue')
    plot_func('IIb', color='grey')
    plot_func('IIn', color='pink')
    plot_func('IIp', color='orange')
    plot_func('bogus_manual', color='red',alpha=0.5)

plt_all_types(SNR_hist)

plt.legend()
plt.title('SNR (diff lc) for each type')
plt.show()


# %%
def linear_fitting(locus: Locus):
    from scipy.stats import linregress

    def surpress_error(f, val=0):
        try:
            return f()
        except Exception:
            return val

    band_mag = lambda band: extract_info(locus.lightcurve, band, None
                                         )[['mjd', 'mag']]
    slope = lambda data: surpress_error(lambda: linregress(data).slope, 0)

    f = lambda band: slope(band_mag(band))

    return (f(Band.R)**2 + f(Band.G)**2)**.5

    # Debugging Code

    # def f(band):
    #     df = band_mag(band)
    #     plt.plot(df['mjd'], df['mag'])
    #     slope = surpress_error(lambda: linregress(df).slope, 0)
    #     intercept = surpress_error(lambda: linregress(df).intercept, 0)

    #     x = np.linspace(0, max(df['mjd']))
    #     y = slope * x + intercept

    #     plt.plot(x, y, color='red')
    #     plt.show()
    #     print(slope)
    # return slope

    # res = (f(Band.R)**2 + f(Band.G)**2)**.5
    # if input() == 'C': raise KeyboardInterrupt()
    # return res


linear_hist = partial(feature_hist, linear_fitting, clip_range=[0, 0.02])
plt_all_types(linear_hist)
plt.legend()
plt.title('linear slope')

plt.show()

# %%
