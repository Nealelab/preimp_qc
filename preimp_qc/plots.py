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

    title = f': {title}' if title else 'QQ Plot'

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
