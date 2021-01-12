import hail as hl
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from hail import Table


def plt_hist(expression: hl.Expression, bins: int = 50, range: list = None, threshold: float = None,
             title: str = None, x_label: str = None, y_label: str = None, log: bool = False):
    exprs = expression.collect()
    df = pd.DataFrame({'col': exprs})

    title = f'{title}' if title else ''
    if log is True:
        df = df[df['col'] != 0]
        plt.hist(np.log10(df.col), edgecolor='black', density=False, bins=bins, color='tab:blue')
    else:
        plt.hist(df['col'], edgecolor='black', density=False, bins=bins, color='tab:blue')
    if threshold: plt.axvline(x=threshold, color='red', linestyle='--')
    if range is not None: plt.xlim(xmin=range[0], xmax=range[1])
    plt.title(title)
    plt.ylabel(y_label if y_label else 'Frequency')
    plt.xlabel(x_label if x_label else '')

    # return type? Zan please assist

    # want to return a plot that can be used in the following ways:
    # crplt = plt_hist(mt_cases.variant_qc)  # this will return the plot
    # option 1
    # cr_base64 = plt_to_base64(crplt) this will convert the plot to base64 so it can be used in HTML report
    # option 2
    # crplt.show()  # show the plot
    # option 3
    # crplt.savefig('crplot.png') save the plot to a file


def qqplot(pvals, title: str = None):
    source = pvals._indices.source
    if isinstance(source, Table):
        ht = source.select(p_value=pvals)
    else:
        ht = source.select_rows(p_value=pvals).rows()

    ht = ht.key_by().select('p_value').key_by('p_value').persist()
    n = ht.count()
    ht = ht.annotate(
        observed_p=-hl.log10(ht['p_value']),
        expected_p=-hl.log10((hl.scan.count() + 1) / n)
    ).persist()

    p_val_pd = ht.to_pandas()
    maxi = max(p_val_pd['expected_p'].max(), p_val_pd['observed_p'].max())

    title = f'{title}' if title else 'QQ Plot'

    plt.scatter(p_val_pd['expected_p'], p_val_pd['observed_p'], c='black', s=0.5)
    plt.plot((0, maxi + 0.4), (0, maxi + 0.4), 'red')
    plt.xlim([0, maxi + 0.5])
    plt.ylim([0, maxi + 0.5])
    plt.title(title)
    plt.ylabel('Observed -log10(' + r'$p$' + ')')
    plt.xlabel('Expected -log10(' + r'$p$' + ')')


def fstat_plot(imputed_sex_ht: hl.Table, f_stat_x: float = 0.4, f_stat_y: float = 0.8):
    fstat_df = imputed_sex_ht.to_pandas()

    fstat_df['is_female'] = fstat_df['is_female'].astype(str)
    fstat_df['is_female'] = fstat_df['is_female'].replace(['True', 'False', 'None'], ['female', 'male', 'unspecified'])

    plt.hist(fstat_df['f_stat'],
             bins=40,
             histtype='bar',
             alpha=0.8,
             fill=True,
             color='tab:blue',
             edgecolor="k")
    plt.axvline(x=f_stat_y, color='red', linestyle='--')
    plt.axvline(x=f_stat_x, color='red', linestyle='--')


def manhattan_plot(pvals, significance_threshold: float = -np.log10(5E-08)):
    source = pvals._indices.source

    if isinstance(source, Table):
        ht = source.select(p_value=pvals)
    else:
        ht = source.select_rows(p_value=pvals).rows()

    data = ht.to_pandas()
    data = data.drop('alleles', 1)  # remove the 'allele' column
    data.columns = ['chromosome', 'position', 'p']  # rename columns
    data['chromosome'].replace({"X": 23, "Y": 24, "MT": 25}, inplace=True)
    data.dropna(subset=['p'], inplace=True)  # drop NAs as log10(val) won't work

    data['-log10(p_value)'] = -np.log10(data['p'])  # compute log10(pvals)
    data['chromosome'] = data['chromosome'].astype('category')
    data['ind'] = range(len(data))
    data_grouped = data.groupby('chromosome')

    fig = plt.figure(figsize=(15, 11))
    ax = fig.add_subplot()

    x_labels = []
    x_labels_pos = []

    colors = ['#E24E42', '#008F95']

    for num, (name, group) in enumerate(data_grouped):
        group.plot(kind='scatter', x='ind', y='-log10(p_value)', marker='o', color=colors[int(name) % len(colors)],
                   ax=ax, s=1000000/len(data))
        x_labels.append(name)
        x_labels_pos.append((group['ind'].iloc[-1] - (group['ind'].iloc[-1] - group['ind'].iloc[0]) / 2))

    ax.set_xticks(x_labels_pos)
    ax.set_xticklabels(x_labels)
    ax.set_xlim([0, len(data)])
    ax.set_ylim([0, data['-log10(p_value)'].max() + 1])
    ax.set_xlabel('Chromosome')
    plt.axhline(y=significance_threshold, color='red', linestyle='--', linewidth=0.5)
    plt.xticks(fontsize=9, rotation=90)
    plt.yticks(fontsize=7)

    # Usage: manhattan_plot(gwas.p_value)