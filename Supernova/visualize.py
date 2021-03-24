# %%
import matplotlib.pyplot as plt
from antares_client._api.models import Locus


# %%
def info(locus: Locus):
    'Prints a summary sentence for this locus.'
    lc = locus.lightcurve

    # two types of alert_id: upper_limit and candidate.
    _ulim_and_candidate = lc.__len__()
    _candidate = lc[lc['alert_id'].apply(
        lambda x: x.startswith('ztf_candidate'))].__len__()
    lc['alert_type'] = lc['alert_id'].apply(lambda x: x[:x.index(':')])

    _julian_dates = lc['ant_mjd']
    date_range = _julian_dates.max() - _julian_dates.min()

    print(
        f'Locus {locus.locus_id} have {_ulim_and_candidate - _candidate} limit alerts '
        + f'and {_candidate} candidate alerts.',
        end='')
    print(f'This event ranges {date_range:.2f} days.')


# %%
# Visualize a lightcurve
def plt_lightcurve(locus: Locus):
    plt.figure(num=None, figsize=(12, 5), dpi=100)
    lc = locus.lightcurve
    lc['alert_type'] = lc['alert_id'].apply(lambda x: x[:x.index(':')])
    for (pb, at), df in lc.groupby(['ant_passband', 'alert_type']):
        is_candidate = at == 'ztf_candidate'
        plt.errorbar(
            x=df['ant_mjd'],
            y=df['ant_mag'] if is_candidate else df['ant_maglim'],
            yerr=df['ant_magerr'],
            #  uplims=(at!='ztf_candidate'),
            label=pb + '  ' + at[4:],
            color=pb,
            fmt=('o' if is_candidate else '^') + '--' + pb.lower(),
            alpha=1 if is_candidate else 0.3)
    plt.title(locus.properties['ztf_object_id'])
    plt.legend()
    plt.show()