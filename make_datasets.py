"""
In this project, I evaluate and develop a filter scheme that
filters out most bogus and keep almost all SN.

Then we use this filter to generate a non-labelled dataset.
"""
# %%
import enum
from Supernova import filter
from Supernova.antares_search import default_search
from Supernova.visualize import plt_lightcurve
from Supernova.collection import fetch_collection
from Supernova import database as db

locus_Types = ['Ia', 'Ib', 'II', 'IIb', 'IIn', 'IIp', 'bogus', 'bogus_manual']

# %%

data = []

for locus_type in locus_Types:
    feed = filter.apply(
        fetch_collection(locus_type),
        filter.date_range(600, 10),
        filter.lightcurve_datapoints(0, 200),
        filter.snr_slope_cut(snr_max=0.25, slope_max=0.0020),
        debug=True,
    )
    try:
        while True:
            next(feed)
    except StopIteration as e:
        counter = e.args[0]
        data.append(counter)


# %%
def generate_plot():
    import numpy as np
    import matplotlib.pyplot as plt

    columns = ('date_range', 'lc_points', 'snr_slope', 'pass')
    rows = locus_Types

    # Get some pastel shades for the colors
    colors = plt.cm.Spectral(np.linspace(0, 0.5, len(rows)))

    n_rows = len(data)

    index = np.arange(len(columns)) + 0.3
    bar_width = 0.4

    # Initialize the vertical-offset for the stacked bar chart.
    y_offset = np.zeros(len(columns))

    data_perc = data.copy()
    for i, row in enumerate(data_perc):
        data_perc[i] = list(map(lambda x: x / sum(row), row))
    print(data_perc)
    # Plot bars and create text labels for the table
    # cell_text = []
    for row in range(n_rows):
        plt.bar(index,
                data_perc[row],
                bar_width,
                bottom=y_offset,
                color=colors[row],
                alpha=1)
        y_offset = y_offset + data_perc[row]
        # cell_text.append(['%d' % (x) for x in data[row]])

    # Add a table at the bottom of the axes
    the_table = plt.table(cellText=data,
                          rowLabels=rows,
                          rowColours=colors,
                          colLabels=columns,
                          loc='bottom')

    # Adjust layout to make room for the table:
    plt.subplots_adjust(left=0.2, bottom=0.2)

    plt.ylabel("Percentage")
    plt.xticks([])
    plt.title('Yield Rate Analysis of Filter')

    plt.show()


# %%
generate_plot()

# %%
for type, type_data in zip(locus_Types, data):
    pass_count = type_data[-1]
    total_count = sum(type_data)
    print(f'Type {type} Summary'.center(40, '-'))
    print(f'True Positive: {pass_count/total_count*100:.2f}%')
    print(f'False Negative: {(1-pass_count/total_count)*100:.2f}%')
# %%

expected_feed_size = 100

feed = filter.apply(
    default_search(),
    filter.date_range(600, 10),
    filter.lightcurve_datapoints(0, 200),
    filter.snr_slope_cut(snr_max=0.25, slope_max=0.0020),
    debug=True,
)

with fetch_collection('dataset_1') as c:
    while len(c) < expected_feed_size:
        locus = next(feed)
        c.add(locus)
        db.add(locus)
        print('add ', locus.locus_id)
