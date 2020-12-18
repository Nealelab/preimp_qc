from plotnine import *
import hail as hl
import matplotlib.pyplot as plt
from qqman import qqman


def plt_cr(df, threshold, title, prepost):

    plt_cr = ggplot(df, aes(x='call_rate')) + \
             geom_histogram(bins=40, color="black", fill="blue") + \
             geom_vline(xintercept=1 - threshold, linetype="dashed", color="red") + \
             labs(title="{} {}-QC".format(title, prepost), y="Frequency") + \
             theme_bw()

    return plt_cr


def cr_var_plts(mt, geno, prepost):

    from .test_qc import stats_split_mt

    mt_cases, mt_controls = stats_split_mt(mt)
    cas_cr_ht_row = mt_cases.rows()
    con_cr_ht_row = mt_controls.rows()

    cas_cr_var = cas_cr_ht_row.select(cas_cr_ht_row.variant_qc.call_rate)
    con_cr_var = con_cr_ht_row.select(con_cr_ht_row.variant_qc.call_rate)

    # convert the ht to pd (this conversion is the one taking most time, FIND EFFICIENT WAY)
    cas_var_df = cas_cr_var.to_pandas()
    con_var_df = con_cr_var.to_pandas()

    cas_var_plt = plt_cr(cas_var_df, geno, "Cases SNP Call Rate", prepost)
    con_var_plt = plt_cr(con_var_df, geno, "Controls SNP Call Rate", prepost)

    return cas_var_plt, con_var_plt


def cr_plts(mt_cases, mt_controls, mind, geno, prepost):
    cas_cr_ht_col = mt_cases.cols()
    con_cr_ht_col = mt_controls.cols()

    # select the call_rate col for both samples and snps
    cas_cr_id = cas_cr_ht_col.select(cas_cr_ht_col.sample_qc.call_rate)
    con_cr_id = con_cr_ht_col.select(con_cr_ht_col.sample_qc.call_rate)

    # convert the ht to pd (this conversion is the one taking most time, FIND EFFICIENT WAY)
    cas_id_df = cas_cr_id.to_pandas()
    con_id_df = con_cr_id.to_pandas()

    cas_id_plt = plt_cr(cas_id_df, mind, "Cases Sample Call Rate", prepost)
    con_id_plt = plt_cr(con_id_df, mind, "Controls Sample Call Rate", prepost)

    return cas_id_plt, con_id_plt


def fstat_plt(imputed_sex_ht, female_thresh, male_thresh, prepost):
    fstat_df = imputed_sex_ht.to_pandas()
    fstat_df['is_female'] = fstat_df['is_female'].astype(str)
    fstat_df['is_female'] = fstat_df['is_female'].replace(['True', 'False', 'None'], ['female', 'male', 'unspecified'])

    sex_colors = {"male": "blue", "female": "purple", "unspecified": "red"}
    f_stat_plot = ggplot(fstat_df, aes(x='f_stat', fill='is_female')) + \
                  geom_histogram(bins=30, color="black") + \
                  geom_vline(xintercept=male_thresh, linetype="dashed", color="blue") + \
                  geom_vline(xintercept=female_thresh, linetype="dashed", color="purple") + \
                  scale_fill_manual(name="Sex", values=sex_colors) + \
                  labs(title="F-statistic distribution {}-QC".format(prepost),
                       x="F-statistic", y="Frequency") + \
                  theme_bw()

    return f_stat_plot


def man_qq_plts(mt, prepost):

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

    figure, axes = plt.subplots(nrows=1, ncols=2, figsize=(25, 10))
    qqman.manhattan(man_df_pruned, ax=axes[0], xrotation=90.0, title="Manhattan plot {}-QC".format(prepost))
    qqman.qqplot(man_df_pruned, ax=axes[1], title="QQ plot {}-QC".format(prepost))

    return axes
