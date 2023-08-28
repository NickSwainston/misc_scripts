import pandas as pd
import matplotlib.pyplot as plt

def make_histogram_plots(all_data, hist_range, label, titles, plotname, xlabel):
    # Make histogram plots
    n_bins = 20
    colours = [
        'blue',
        'green',
        'orange',
        'purple',
    ]
    characters = [
        'a)',
        'b)',
        'c)',
        'd)',
        'e)'
    ]
    n_data = len(all_data) + 1

    fig, axes = plt.subplots(nrows=n_data, figsize=(5, 3*n_data))

    axes[0].hist(all_data, n_bins, density=True, histtype='bar', stacked=True, label=label, color=colours[:n_data-1])
    axes[0].text(0.1, 0.8, characters[0], transform=axes[0].transAxes, size=20)
    axes[0].set_title(titles[0])
    axes[0].legend(prop={'size': 10})
    axes[0].set_ylabel(f"Probability Density")
    axes[0].set_xlabel(f"${xlabel}$")

    for ai, df_col, colour, title in zip(range(1, n_data), all_data, colours, titles[1:]):
        #print(ai, n_data)
        axes[ai].hist(df_col, n_bins, histtype='bar', color=colour, range=hist_range)
        axes[ai].text(0.1, 0.8, characters[ai], transform=axes[ai].transAxes, size=20)
        axes[ai].set_title(title)
        axes[ai].set_ylabel(f"#")
        axes[ai].set_xlabel(f"${xlabel}$")

    fig.tight_layout()
    fig.savefig(plotname)
    plt.close(fig)

# Swainston data
df = pd.read_csv('all_pulsar_fits.csv')
spl_df   = df[df["Model"] == "simple_power_law"]

# Sett data
df = pd.read_csv('sp_index.csv')
# Alpha histogram
all_indexs = [
    spl_df ["a"],
    df["specindex"],
]
hist_range = (spl_df["a"].min(), spl_df["a"].max())
titles = [
    "Comparison",
    'Swainston et al. 2013',
    'This work',
]
make_histogram_plots(all_indexs, hist_range, label=["Swainston et al. 2013", "This work"], titles=titles, plotname="spectral_index_histogram.png", xlabel="\\alpha")
