import hail as hl
from typing import List, Tuple

import preimp_qc.test_qc as qc
import preimp_qc.test_plots as plt


def compute_qc_metrics(mt: hl.MatrixTable) -> hl.MatrixTable:
    """
    Compute per-sample metrics and common variant statistics useful for quality control
    :param mt: Hail MatrixTable
    :return: Hail MatrixTable with variant and sample qc metrics
    """
    mt = hl.variant_qc(mt)
    mt = hl.sample_qc(mt)

    return mt


# TESTED: WORKING
def stats_split_mt(mt: hl.MatrixTable) -> Tuple[List[int], hl.MatrixTable, hl.MatrixTable]:
    """
    Collect basic stat counts and split the MatrixTable by phenotype status (for plots)
    :param mt: Hail MatrixTable
    :return: basic stat counts, cases MatrixTable, controls MatrixTable
    """
    # 1. Sex
    n_females: List[str] = mt.filter_cols(mt.is_female == True).s.collect()
    n_males: List[str] = mt.filter_cols(mt.is_female == False).s.collect()
    n_sex_missing: List[str] = mt.filter_cols(hl.is_missing(mt.is_female)).s.collect()

    # 2. Phenotype status
    mt_cases: hl.MatrixTable = mt.filter_cols(mt.is_case == True)
    n_cases: int = mt_cases.count_cols()
    mt_controls: hl.MatrixTable = mt.filter_cols(mt.is_case == False)
    n_controls: int = mt_controls.count_cols()
    n_unknown_pheno: List[str] = mt.filter_cols(hl.is_missing(mt.is_case)).s.collect()

    # 3. Number of SNPs
    n_snps = mt.count_rows()

    counts: List[int] = [len(n_males), len(n_females), len(n_sex_missing), n_cases,
                         n_controls, len(n_unknown_pheno), n_snps]

    return counts, mt_cases, mt_controls


def run_qc(mt: hl.MatrixTable, dirname: str, basename: str, input_type: str, pre_geno: float, mind: float, fhet_y: int,
           fhet_x: int, geno: float, midi: float, maf: float, hwe_th_co: float, hwe_th_ca: float, qc_round: int,
           withpna: int = 0) -> hl.MatrixTable:
    """
    :param mt: Hail MatrixTable
    :param dirname:
    :param basename:
    :param input_type:
    :param pre_geno:
    :param mind:
    :param fhet_y:
    :param fhet_x:
    :param geno:
    :param midi:
    :param maf:
    :param hwe_th_co:
    :param hwe_th_ca:
    :param qc_round:
    :param withpna:
    :return:
    """

    # compute qc metrics
    mt = qc.compute_qc_metrics(mt)

    # Pre-qc counts
    pre_qc_counts = qc.collect_counts(mt)

    # pre-qc plots
    print("Generating pre-QC plots")
    pre_cas_var_base64, pre_con_var_base64 = plt.cr_var_plts(mt, geno)
    pre_cas_id_base64, pre_con_id_base64 = plt.cr_id_plts(mt, mind)

    pre_man_qq_base64 = plt.man_qq_plts(mt)

    # 1. SNP QC: call rate ≥ 0.95
    print("1. SNP QC: call rate ≥ 0.95")
    mt, var_pre_filter = qc.filter_var_cr(mt, pre_geno)
    print("Pre QC call rate < 0.95: {}".format(var_pre_filter['geno_removed']))
    print("Samples: {}".format(mt.count_cols()))

    # 2. Sample QC: call rate in cases or controls ≥ 0.98
    print("2. Sample QC: call rate in cases or controls ≥ 0.98")
    mt, id_cr_filter = qc.filter_sample_cr(mt, mind)
    print("Sample QC < 0.98: {}".format(id_cr_filter['sample_miss_cases'] + id_cr_filter['sample_miss_controls']))
    print("Samples: {}".format(mt.count_cols()))

    # 3. Sample QC: F_stats
    print("3. Sample QC: F_stats")
    
    mt, f_stat_results = qc.filter_sex_check(mt, fhet_y, fhet_x)
    print("Sex check filtered: {}".format(f_stat_results['sex_check_removed']))
    print("Samples: {}".format(mt.count_cols()))

    # 4. Sample QC: Sex violations (excluded) - genetic sex does not match pedigree sex
    print("4. Sample QC: Sex violations (excluded) - genetic sex does not match pedigree sex")
    mt, sex_violations = qc.sex_violations(mt, input_type)
    print("Sex violations: {}".format(sex_violations['sex_excluded']))
    print("Samples: {}".format(mt.count_cols()))

    # 5. Sample QC: Sex warnings (not excluded) - undefined phenotype / ambiguous genotypes
    print("# 5. Sample QC: Sex warnings (not excluded) - undefined phenotype / ambiguous genotypes")
    sex_warnings_count = qc.sex_warnings(mt, input_type)
    print("Sex warning: {}".format(sex_warnings_count))
    print("Samples: {}".format(mt.count_cols()))

    # 6. SNP QC: call rate ≥ 0.98
    print("# 6. SNP QC: call rate ≥ 0.98")
    mt, var_filter = qc.filter_var_cr(mt, geno)
    print("SNP QC call rate < 0.98: {}".format(var_filter['geno_removed']))
    print("Samples: {}".format(mt.count_cols()))

    # 7. SNP QC: missing difference > 0.02
    print("# 7. SNP QC: missing difference > 0.02")

    # 8. SNP QC: SNPs with no valid association p value are excluded (i.e., invariant SNP)
    print("# 8. SNP QC: SNPs with no valid association p value are excluded (i.e., invariant SNP)")
    if withpna == 0:
        mt, invariant_snps = qc.filter_invariant_snps(mt)
        print("Monormorphic SNPs: {}".format(invariant_snps['monomorphic_snps']))
        print("Samples: {}".format(mt.count_cols()))

    # 9. SNP QC: with MAF ≥ 0.01
    print("# 9. SNP QC: with MAF ≥ 0.01")
    mt, maf_results = qc.filter_maf(mt, maf)
    print("MAF: {}".format(maf_results['maf_removed']))
    print("Samples: {}".format(mt.count_cols()))

    # 10. SNP QC: Hardy-Weinberg equilibrium (HWE) in controls p value ≥ 1e-06
    print("# 10. SNP QC: Hardy-Weinberg equilibrium (HWE) in controls p value ≥ 1e-06")
    mt, hwe_con_results = qc.filter_hwe(mt, 'Control', hwe_th_co)
    print("HWE Controls: {}".format(hwe_con_results['maf_removed']))
    print("Samples: {}".format(mt.count_cols()))

    # 11. SNP QC: Hardy-Weinberg equilibrium (HWE) in cases p value ≥ 1e-10
    print("# 11. SNP QC: Hardy-Weinberg equilibrium (HWE) in cases p value ≥ 1e-10")
    mt, hwe_cas_results = qc.filter_hwe(mt, 'Case', hwe_th_ca)
    print("HWE Cases: {}".format(hwe_cas_results['maf_removed']))
    print("Samples: {}".format(mt.count_cols()))

    # Post-qc counts
    post_qc_counts = qc.collect_counts(mt)

    # Post-QC plots
    print("Generating post-QC plots")
    print("Generating variant call rate plots")
    pos_cas_var_base64, pos_con_var_base64 = plt.cr_var_plts(mt, geno)
    print("Generating sample call rate plots")
    pos_cas_id_base64, pos_con_id_base64 = plt.cr_id_plts(mt, mind)
    print("Generating Manhattand & QQ plots")
    pos_man_qq_base64 = plt.man_qq_plts(mt)

    # # pre_cas_var_base64, pre_cas_id_base64, pre_con_var_base64, pre_con_id_base64
    qc_plots_list = [pre_man_qq_base64, pos_man_qq_base64, pre_con_id_base64, pre_cas_id_base64, pos_con_id_base64, pos_cas_id_base64,
                     f_stat_results['sex_check_plot'], pre_con_var_base64, pre_cas_var_base64, pos_con_var_base64, pos_cas_var_base64]

    # Tables
    filter_counts_list = [var_pre_filter['geno_removed'], id_cr_filter['sample_miss_cases'] + id_cr_filter['sample_miss_controls'],
                          f_stat_results['sex_check_removed'], sex_violations['sex_excluded'], sex_warnings_count,
                          var_filter['geno_removed'], invariant_snps['monomorphic_snps'],
                          hwe_con_results['maf_removed'], hwe_cas_results['maf_removed']]
    size_of_sample_html, exlusion_overview_html = generate_tables(pre_qc_counts, post_qc_counts, filter_counts_list)
    qc_tables_list = [size_of_sample_html, exlusion_overview_html]

    outplink = dirname + basename + '_qc{}'.format(qc_round)
    hl.export_plink(mt, outplink)

    return qc_tables_list, qc_plots_list


def generate_tables(pre_qc_counts, post_qc_counts, filter_counts):
    import pandas as pd
    # this function is for generating the tables in the Generak Info section
    # HAVE TO WORK ON FLAGS TABLE

    # Size of Sample table
    pre_qc_pheno_counts = [pre_qc_counts[3], pre_qc_counts[4], pre_qc_counts[5]]
    post_qc_pheno_counts = [post_qc_counts[3], post_qc_counts[4], post_qc_counts[5]]
    ex_pheno = [x1 - x2 for x1, x2 in zip(pre_qc_pheno_counts, post_qc_pheno_counts)]
    pre_qc_sex_counts = [pre_qc_counts[0], pre_qc_counts[1], pre_qc_counts[2]]
    post_qc_sex_counts = [post_qc_counts[0], post_qc_counts[1], post_qc_counts[2]]
    ex_sex = [x1 - x2 for x1, x2 in zip(pre_qc_sex_counts, post_qc_sex_counts)]
    n_snps_pre = pre_qc_counts[6]
    n_snps_post = post_qc_counts[6]
    ex_snps = n_snps_pre - n_snps_post
    size_of_sample = [['Cases,Controls,Missing', pre_qc_pheno_counts, post_qc_pheno_counts, ex_pheno],
                      ['Males,Females,Unspec', pre_qc_sex_counts, post_qc_sex_counts, ex_sex],
                      ['SNPs', n_snps_pre, n_snps_post, ex_snps]]
    size_of_sampledf = pd.DataFrame(size_of_sample, columns=['Test', 'pre QC', 'post QC', 'exclusion-N'])
    size_of_sample_html = size_of_sampledf.to_html()

    # Exlusion overview table
    snps_cr_95 = ['SNPs: call rate < 0.950 (pre-filter)', filter_counts[0]]
    ids_cr = ['IDS: call rate (cases/controls) < 0.980', filter_counts[1]]
    ids_fhet = ['IDs: FHET outside +- 0.20 (cases/controls)', filter_counts[2]]
    ids_sex_violations = ['IDs: Sex violations -excluded- (N-tested)', filter_counts[3]]
    ids_sex_warnings = ['IDs: Sex warnings (undefined genotype/ambigous genotype)', filter_counts[4]]
    snps_cr_98 = ['SNPs: call rate < 0.980', filter_counts[5]]
    snps_monomorphic = ['SNPs: without valid association p-value (invariant)', filter_counts[6]]
    snps_hwe_con = ['SNPs: HWE-controls < -6', filter_counts[7]]
    snps_hwe_cas = ['SNPs: HWE-cases < -10', filter_counts[8]]
    exlusion_overview = [snps_cr_95, ids_cr, ids_fhet, ids_sex_violations, ids_sex_warnings, snps_cr_98,
                         snps_monomorphic, snps_hwe_con, snps_hwe_cas]
    exlusion_overviewdf = pd.DataFrame(exlusion_overview, columns=['Filter', 'N'])
    exlusion_overview_html = exlusion_overviewdf.to_html()

    return size_of_sample_html, exlusion_overview_html