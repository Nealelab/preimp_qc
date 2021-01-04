from plotnine import *
import hail as hl
import matplotlib.pyplot as plt
from qqman import qqman
import base64
import io


def plt_to_base64(plt):
    buffer = io.BytesIO()
    plt.save(buffer, format='PNG', verbose=False)
    buffer.seek(0)

    plt_base64 = base64.b64encode(buffer.read()).decode('ascii')
    return '<img src="data:image/png;base64,{}">'.format(plt_base64)


def plt_cr(df, threshold, title):

    plt_cr = ggplot(df, aes(x='call_rate')) + \
             geom_histogram(bins=40, color="black", fill="blue") + \
             geom_vline(xintercept=1 - threshold, linetype="dashed", color="red") + \
             labs(title=title, y="Frequency") + \
             theme_bw()

    return plt_cr


def cr_var_plts(mt, geno):

    from .test_qc import stats_split_mt

    mt_cases, mt_controls = stats_split_mt(mt)
    cas_cr_ht_row = mt_cases.rows()
    con_cr_ht_row = mt_controls.rows()

    cas_cr_var = cas_cr_ht_row.select(cas_cr_ht_row.variant_qc.call_rate)
    con_cr_var = con_cr_ht_row.select(con_cr_ht_row.variant_qc.call_rate)

    # convert the ht to pd (this conversion is the one taking most time, FIND EFFICIENT WAY)
    cas_var_df = cas_cr_var.to_pandas()
    con_var_df = con_cr_var.to_pandas()

    cas_var_plt = plt_cr(cas_var_df, geno, "Cases variant call rate")
    cas_var_plt64 = plt_to_base64(cas_var_plt)
    con_var_plt = plt_cr(con_var_df, geno, "Controls variant call rate")
    con_var_plt64 = plt_to_base64(con_var_plt)

    return cas_var_plt64, con_var_plt64


def cr_id_plts(mt, mind):
    from .test_qc import stats_split_mt
    mt_cases, mt_controls = stats_split_mt(mt)
    cas_cr_ht_col = mt_cases.cols()
    con_cr_ht_col = mt_controls.cols()

    # select the call_rate col for both samples and snps
    cas_cr_id = cas_cr_ht_col.select(cas_cr_ht_col.sample_qc.call_rate)
    con_cr_id = con_cr_ht_col.select(con_cr_ht_col.sample_qc.call_rate)

    # convert the ht to pd (this conversion is the one taking most time, FIND EFFICIENT WAY)
    cas_id_df = cas_cr_id.to_pandas()
    con_id_df = con_cr_id.to_pandas()

    cas_id_plt = plt_cr(cas_id_df, mind, "Cases sample call rate")
    cas_id_plt64 = plt_to_base64(cas_id_plt)
    con_id_plt = plt_cr(con_id_df, mind, "Controls sample call rate")
    con_id_plt64 = plt_to_base64(con_id_plt)

    return cas_id_plt64, con_id_plt64


def fstat_plt(imputed_sex_ht, female_thresh, male_thresh):
    fstat_df = imputed_sex_ht.to_pandas()
    fstat_df['is_female'] = fstat_df['is_female'].astype(str)
    fstat_df['is_female'] = fstat_df['is_female'].replace(['True', 'False', 'None'], ['female', 'male', 'unspecified'])

    sex_colors = {"male": "blue", "female": "purple", "unspecified": "red"}
    f_stat_plot = ggplot(fstat_df, aes(x='f_stat', fill='is_female')) + \
                  geom_histogram(bins=30, color="black") + \
                  geom_vline(xintercept=male_thresh, linetype="dashed", color="blue") + \
                  geom_vline(xintercept=female_thresh, linetype="dashed", color="purple") + \
                  scale_fill_manual(name="Sex", values=sex_colors) + \
                  labs(title="F-statistic distribution",
                       x="F-statistic", y="Frequency") + \
                  theme_bw()

    f_stat_plot64 = plt_to_base64(f_stat_plot)

    return f_stat_plot64


def man_qq_plts(mt):

    gwas_ht = hl.linear_regression_rows(y=mt.is_case,
                                        x=mt.GT.n_alt_alleles(),
                                        covariates=[1.0])

    pvals = gwas_ht.select(gwas_ht.p_value)
    man_df = pvals.to_pandas()

    man_df_pruned = man_df[['locus.contig', 'locus.position', 'p_value']]
    man_df_pruned.columns = ['CHR', 'BP', 'P']
    man_df_pruned = man_df_pruned.dropna()

    man_df_pruned = man_df_pruned.replace(to_replace=["X", "Y", "MT"],
                                          value=[23, 24, 25])

    buffer = io.BytesIO()
    figure, axes = plt.subplots(nrows=1, ncols=2, figsize=(25, 10))
    qqman.manhattan(man_df_pruned, ax=axes[0], xrotation=90.0, title="Manhattan plot")
    qqman.qqplot(man_df_pruned, ax=axes[1], title="QQ plot")

    figure.tight_layout()
    plt.savefig(buffer, format='PNG')
    plt.clf()
    plt.close()
    buffer.seek(0)

    plt_base64 = base64.b64encode(buffer.read()).decode('ascii')
    return '<img src="data:image/png;base64,{}">'.format(plt_base64)
