from plotnine import *
import hail as hl
import matplotlib.pyplot as plt
import base64
import io
from qqman import qqman


def cr_plts(mt_cases, mt_controls, mind, geno, prepost):
    import base64
    import io
    # takes bout 6 mins to generate all 4 plots, bokeh in Hail takes ~16 minutes (~4 mins each plot)
    # run using: cr_plts(mt_cases, mt_controls, 0.98, 0.98, 'pre', '/Users/lindokuhle/Desktop/preimpQC')
    # convert mts to tables
    cas_cr_ht_row = mt_cases.rows()
    cas_cr_ht_col = mt_cases.cols()
    con_cr_ht_row = mt_controls.rows()
    con_cr_ht_col = mt_controls.cols()

    # select the call_rate col for both samples and snps
    cas_cr_var = cas_cr_ht_row.select(cas_cr_ht_row.variant_qc.call_rate)
    cas_cr_id = cas_cr_ht_col.select(cas_cr_ht_col.sample_qc.call_rate)
    con_cr_var = con_cr_ht_row.select(con_cr_ht_row.variant_qc.call_rate)
    con_cr_id = con_cr_ht_col.select(con_cr_ht_col.sample_qc.call_rate)

    # convert the ht to pd (this conversion is the one taking most time, FIND EFFICIENT WAY)
    cas_var_df = cas_cr_var.to_pandas()
    cas_id_df = cas_cr_id.to_pandas()
    con_var_df = con_cr_var.to_pandas()
    con_id_df = con_cr_id.to_pandas()

    def plt_cr(df, threshold, title, prepost):
        plt = ggplot(df, aes(x='call_rate')) + \
              geom_histogram(bins=40, color="black", fill="blue") + \
              geom_vline(xintercept=1-threshold, linetype="dashed", color="red") + \
              labs(title="{} {}-QC".format(title, prepost), y="Frequency") + \
              theme_bw()
        return plt

    cas_var_plt = plt_cr(cas_var_df, geno, "Cases SNP Call Rate", prepost)
    cas_id_plt = plt_cr(cas_id_df, mind, "Cases Sample Call Rate", prepost)
    con_var_plt = plt_cr(con_var_df, geno, "Controls SNP Call Rate", prepost)
    con_id_plt = plt_cr(con_id_df, mind, "Controls Sample Call Rate", prepost)

    buffer1 = io.BytesIO()
    cas_var_plt.save(buffer1, format='PNG', verbose=False)
    buffer1.seek(0)
    cas_var_plt_base64 = base64.b64encode(buffer1.read()).decode('ascii')
    cas_var_base64 = '<img src="data:image/png;base64,{}">'.format(cas_var_plt_base64)

    buffer2 = io.BytesIO()
    cas_id_plt.save(buffer2, format='PNG', verbose=False)
    buffer2.seek(0)
    cas_id_plt_base64 = base64.b64encode(buffer2.read()).decode('ascii')
    cas_id_base64 = '<img src="data:image/png;base64,{}">'.format(cas_id_plt_base64)

    buffer3 = io.BytesIO()
    con_var_plt.save(buffer3, format='PNG', verbose=False)
    buffer3.seek(0)
    con_var_plt_base64 = base64.b64encode(buffer3.read()).decode('ascii')
    con_var_base64 = '<img src="data:image/png;base64,{}">'.format(con_var_plt_base64)

    buffer4 = io.BytesIO()
    con_id_plt.save(buffer4, format='PNG', verbose=False)
    buffer4.seek(0)
    con_id_plt_base64 = base64.b64encode(buffer4.read()).decode('ascii')
    con_id_base64 = '<img src="data:image/png;base64,{}">'.format(con_id_plt_base64)

    # cas_var_plt.save(filename="{}casSNPcr{}QC.png".format(outdir, prepost), height=5, width=5, units='in', dpi=1000)
    # cas_id_plt.save(filename="{}casIDcr{}QC.png".format(outdir, prepost), height=5, width=5, units='in', dpi=1000)
    # con_var_plt.save(filename="{}conSNPcr{}QC.png".format(outdir, prepost), height=5, width=5, units='in', dpi=1000)
    # con_id_plt.save(filename="{}conIDPcr{}QC.png".format(outdir, prepost), height=5, width=5, units='in', dpi=1000)

    return cas_var_base64, cas_id_base64, con_var_base64, con_id_base64


def fstat_plt(imputed_sex_ht, female_thresh, male_thresh, prepost):
    fstat_df = imputed_sex_ht.to_pandas()  # convert the imputed_sex table to a pandas df
    fstat_df['is_female'] = fstat_df['is_female'].astype(str)  # convert bools to str
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

    buffer = io.BytesIO()
    f_stat_plot.save(buffer, format='PNG')
    buffer.seek(0)

    plt_base64 = base64.b64encode(buffer.read()).decode('ascii')
    # outpltdir = outdir + "FstatDistribution" + prepost + "QC.png"
    # f_stat_plot.save(filename=outpltdir, height=5, width=5, units='in', dpi=1000)
    return '<img src="data:image/png;base64,{}">'.format(plt_base64)


def man_qq_plts(mt, prepost):
    # run using: man_qq_plts(mt, '/Users/lindokuhle/Desktop/preimpQC/', 'pre')
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
    qqman.manhattan(man_df_pruned, ax=axes[0], xrotation=90.0, title="Manhattan plot {}-QC".format(prepost))
    qqman.qqplot(man_df_pruned, ax=axes[1], title="QQ plot {}-QC".format(prepost))

    figure.tight_layout()
    plt.savefig(buffer, format='PNG')
    plt.clf()
    plt.close()
    buffer.seek(0)

    plt_base64 = base64.b64encode(buffer.read()).decode('ascii')
    # plt.savefig("./manhattan.png", format="png")
    # qqman.manhattan(man_df_pruned, xrotation=90.0, out="{}Man{}QC.png".format(outdir, prepost),
    # title="Manhattan plot {}-QC".format(prepost))
    # qqman.qqplot(man_df_pruned, out="{}QQplot{}QC.png".format(outdir, prepost), title="QQ plot {}-QC".format(prepost))
    return '<img src="data:image/png;base64,{}">'.format(plt_base64)
