import hail as hl
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from hail import Table
import base64
import io


def matplotlib_to_base64(plot):
    # converts a matplotlib fig to a base64 image (good for HTML reports)
    buffer = io.BytesIO()
    plot.savefig(buffer, format='PNG')
    buffer.seek(0)

    plt_base64 = base64.b64encode(buffer.read()).decode('ascii')
    return '<img src="data:image/png;base64,{}">'.format(plt_base64)


# TESTED AND WORKING
def plt_hist(expression: hl.Expression, bins: int = 50, pltrange: list = None, threshold: float = None,
             title: str = None, x_label: str = None, y_label: str = None, log: bool = False, figsize: tuple = (12, 8)):
    exprs = expression.collect()
    df = pd.DataFrame({'col': exprs})

    title = f'{title}' if title else ''
    fig = plt.figure(figsize=figsize)
    if log is True:
        df = df[df['col'] != 0]
        plt.hist(np.log10(df.col), edgecolor='black', density=False, bins=bins, color='tab:blue')
    else:
        plt.hist(df['col'], edgecolor='black', density=False, bins=bins, color='tab:blue')
    if threshold: plt.axvline(x=threshold, color='red', linestyle='--')
    if pltrange is not None: plt.xlim(xmin=pltrange[0], xmax=pltrange[1])
    plt.title(title)
    plt.ylabel(y_label if y_label else 'Frequency')
    plt.xlabel(x_label if x_label else '')
    plt.close()

    return fig


# TESTED AND WORKING
def qqplot(pvals, title: str = None, figsize: tuple = (10, 10)):
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

    fig = plt.figure(figsize=figsize)
    plt.scatter(p_val_pd['expected_p'], p_val_pd['observed_p'], c='black', s=0.5)
    plt.plot((0, maxi + 0.4), (0, maxi + 0.4), 'red')
    plt.xlim([0, maxi + 0.5])
    plt.ylim([0, maxi + 0.5])
    plt.title(title)
    plt.ylabel('Observed -log10(' + r'$p$' + ')')
    plt.xlabel('Expected -log10(' + r'$p$' + ')')
    plt.close()

    return fig


# TESTED AND WORKING
def fstat_plot(imputed_sex_ht: hl.Table, f_stat_x: float = 0.4, f_stat_y: float = 0.8, figsize: tuple = (12, 8)):
    fstat_df = imputed_sex_ht.to_pandas()

    fstat_df['is_female'] = fstat_df['is_female'].astype(str)
    fstat_df['is_female'] = fstat_df['is_female'].replace(['True', 'False', 'None'], ['female', 'male', 'unspecified'])

    fig = plt.figure(figsize=figsize)
    plt.hist(fstat_df['f_stat'],
             bins=40,
             histtype='bar',
             alpha=0.8,
             fill=True,
             color='tab:blue',
             edgecolor="k")
    plt.axvline(x=f_stat_y, color='red', linestyle='--')
    plt.axvline(x=f_stat_x, color='red', linestyle='--')
    plt.close()

    return fig


# TESTED AND WORKING
def manhattan_plot(pvals, significance_threshold: float = -np.log10(5E-08), title: str = None,
                   figsize: tuple = (17, 11), annotate_sig: bool = False):
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

    title = f'{title}' if title else 'Manhattan Plot'

    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot()

    x_labels = []
    x_labels_pos = []

    colors = ['#E24E42', '#008F95']

    for num, (name, group) in enumerate(data_grouped):
        group.plot(kind='scatter', x='ind', y='-log10(p_value)', marker='o', color=colors[int(name) % len(colors)],
                   ax=ax, s=1000000 / len(data))
        x_labels.append(name)
        x_labels_pos.append((group['ind'].iloc[-1] - (group['ind'].iloc[-1] - group['ind'].iloc[0]) / 2))

    ax.set_xticks(x_labels_pos)
    ax.set_xticklabels(x_labels)
    ax.set_xlim([0, len(data)])
    ax.set_ylim([0, data['-log10(p_value)'].max() + 1])
    ax.set_xlabel('Chromosome')
    plt.title(title)
    plt.axhline(y=significance_threshold, color='red', linestyle='--', linewidth=0.5)
    plt.xticks(fontsize=9, rotation=90)
    plt.yticks(fontsize=7)
    if annotate_sig is True:
        for index, row in data.iterrows():
            if row['-log10(p_value)'] >= significance_threshold:
                ax.annotate('{}:{}'.format(row['chromosome'], row['position']),
                            xy=(index, row['-log10(p_value)'] + 0.1),
                            bbox=dict(boxstyle="round", fc="0.8"), fontsize=8)
                # FIX: the annotation box keeps moving further away from the point (SNP) as you move from left to right
    plt.close()

    # Usage: manhattan_plot(gwas.p_value)
    # saving the plot to a file: man_plt.savefig('test.png', facecolor='w', dpi=300)

    return fig
