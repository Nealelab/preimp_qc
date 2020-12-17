import hail as hl
from typing import List, Tuple
from preimp_qc.plots import cr_plts, man_qq_plts


def plink_to_mt(dirname: str, basename: str, reference: str = 'GRCh38') -> hl.MatrixTable:
    hl.init(default_reference=reference)
    mt: hl.MatrixTable = hl.import_plink(bed=dirname + basename + '.bed',
                                         bim=dirname + basename + '.bim',
                                         fam=dirname + basename + '.fam')
    return mt


# work out a way to deal with .vcf and .vcf.gz/bgz
def vcf_to_mt(dirname: str, vcf: str, annotations: str) -> hl.MatrixTable:
    hl.import_vcf(vcf).write('{}preimpQC.mt'.format(dirname), overwrite=True)
    mt = hl.read_matrix_table('{}preimpQC.mt'.format(dirname))
    ann = hl.import_table(annotations, impute=True).key_by('Sample')
    mt = mt.annotate_cols(annotations=ann[mt.s])
    # need to change reported sex to True/False, can update how this is done later, ideally don't want to hardcode
    # this will not work for unreported sex but will work for missing values
    mt = mt.annotate_cols(annotations=mt.annotate(Sex=hl.if_else((mt.annotations.Sex == 'F' | mt.annotations.Sex == 2 |
                                                                  mt.annotations.Sex == 'Female'), True, False)))
    # add a check to make sure file is formatted as we'd expect else quit and throw error
    return mt


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

    # Pre-qc counts
    pre_qc_counts, pre_mt_cas, pre_mt_con = stats_split_mt(mt)

    # pre-qc plots
    print("Generating pre-QC plots")
    pre_cas_var_base64, pre_cas_id_base64, pre_con_var_base64, pre_con_id_base64 = cr_plts(pre_mt_cas, pre_mt_con,
                                                                                           mind,
                                                                                           geno, "pre")
    pre_man_qq_base64 = man_qq_plts(mt, "pre")

    # 1. SNP QC: call rate ≥ 0.95
    print("1. SNP QC: call rate ≥ 0.95")
    pre_geno_cr = mt.filter_rows(mt.variant_qc.call_rate < (1 - pre_geno)).rsid.collect()
    if len(pre_geno_cr) > 0:
        mt = mt.filter_rows(hl.literal(pre_geno_cr).contains(mt['rsid']), keep=False)

    # 2. Sample QC: call rate in cases or controls ≥ 0.98
    print("2. Sample QC: call rate in cases or controls ≥ 0.98")
    sample_miss = mt.filter_cols(mt.sample_qc.call_rate < (1 - mind))
    sample_miss_cases = sample_miss.filter_cols(sample_miss.is_case == True).s.collect()
    sample_miss_controls = sample_miss.filter_cols(sample_miss.is_case == False).s.collect()
    n_sample_miss = len(sample_miss_cases) + len(sample_miss_controls)
    samples_miss = sample_miss_cases + sample_miss_controls
    if n_sample_miss > 0:
        mt = mt.filter_cols(hl.literal(samples_miss).contains(mt['s']), keep=False)

    # 3. Sample QC: F_stats
    print("3. Sample QC: F_stats")
    imputed_sex = hl.impute_sex(mt.GT)
    from preimp_qc.plots import fstat_plt
    f_stat_plot = fstat_plt(imputed_sex, fhet_y, fhet_x, "pre")
    f_stat_out = mt.filter_cols(((imputed_sex[mt.s].f_stat < fhet_x) & (mt.is_female == False) |
                                 (imputed_sex[mt.s].f_stat > fhet_y) & (mt.is_female == True))).s.collect()
    if len(f_stat_out) > 0:
        mt = mt.filter_cols(hl.literal(f_stat_out).contains(mt['s']), keep=False)

    # 4. Sample QC: Sex violations (excluded) - genetic sex does not match pedigree sex
    print("4. Sample QC: Sex violations (excluded) - genetic sex does not match pedigree sex")
    if input_type == "plink":
        # Verify that when sex info is missing value is set to None
        mt = mt.filter_cols((mt.is_female != imputed_sex[mt.s]) & (mt.is_female != None), keep=False)
    elif input_type == "vcf":
        # Verify that when meta file is read in, column formatting is kept
        mt = mt.filter_cols((mt.annotations.Sex != imputed_sex[mt.s]) & (mt.annotations.Sex != None), keep=False)

    # 5. Sample QC: Sex warnings (not excluded) - undefined phenotype / ambiguous genotypes
    print("# 5. Sample QC: Sex warnings (not excluded) - undefined phenotype / ambiguous genotypes")
    if input_type == "plink":
        undef_count = mt.aggregate_cols(hl.agg.counter(mt.is_female == None))
        print("Warning: {} individuals have undefined phenotype/ambiguous genotypes".format(undef_count))
    elif input_type == "vcf":
        undef_count = mt.aggregate_cols(hl.agg.counter(mt.annotations.Sex == None))
        print("Warning: {} individuals have undefined phenotype/ambiguous genotypes".format(undef_count))
    # Add an output of the individual IDs who have ambiguous sex so they can be removed?
    # Alternatively add option at this step to remove indivs with ambiguous sex?

    # 6. SNP QC: call rate ≥ 0.98
    print("# 6. SNP QC: call rate ≥ 0.98")
    geno_cr = mt.filter_rows(mt.variant_qc.call_rate < (1 - geno)).rsid.collect()
    if len(geno_cr) > 0:
        mt = mt.filter_rows(hl.literal(geno_cr).contains(mt['rsid']), keep=False)

    # 7. SNP QC: missing difference > 0.02
    print("# 7. SNP QC: missing difference > 0.02")


    # 8. SNP QC: SNPs with no valid association p value are excluded (i.e., invariant SNP)
    print("# 8. SNP QC: SNPs with no valid association p value are excluded (i.e., invariant SNP)")
    if withpna == 0:
        mt = mt.annotate_rows(MAC=hl.min(mt.variant_qc.AC))
        monomorphic_snps = mt.filter_rows(mt.MAC == 0).rsid.collect()
        if len(monomorphic_snps) > 0:
            mt = mt.filter_rows(hl.literal(monomorphic_snps).contains(mt['rsid']), keep=False)

    # 9. SNP QC: with MAF ≥ 0.01
    print("# 9. SNP QC: with MAF ≥ 0.01")
    mt = mt.annotate_rows(maf=hl.min(mt.variant_qc.AF))
    snps_maf = mt.filter_rows(mt.maf < maf).rsid.collect()
    if len(snps_maf) > 0:
        mt = mt.filter_rows(hl.literal(snps_maf).contains(mt['rsid']), keep=False)

    # 10. SNP QC: Hardy-Weinberg equilibrium (HWE) in controls p value ≥ 1e-06
    print("# 10. SNP QC: Hardy-Weinberg equilibrium (HWE) in controls p value ≥ 1e-06")
    hwe_controls = mt.filter_cols(mt.is_case == False)
    snps_hwe_controls = hwe_controls.filter_rows(hwe_controls.variant_qc.p_value_hwe < hwe_th_co).rsid.collect()
    if len(snps_hwe_controls) > 0:
        mt = mt.filter_rows(hl.literal(snps_hwe_controls).contains(mt['rsid']), keep=False)

    # 11. SNP QC: Hardy-Weinberg equilibrium (HWE) in cases p value ≥ 1e-10
    print("# 11. SNP QC: Hardy-Weinberg equilibrium (HWE) in cases p value ≥ 1e-10")
    hwe_cases = mt.filter_cols(mt.is_case == True)
    snps_hwe_cases = hwe_cases.filter_rows(hwe_cases.variant_qc.p_value_hwe < hwe_th_ca).rsid.collect()
    if len(snps_hwe_cases) > 0:
        mt = mt.filter_rows(hl.literal(snps_hwe_cases).contains(mt['rsid']), keep=False)

    # Post-qc counts
    post_qc_counts, post_mt_cas, post_mt_con = stats_split_mt(mt)

    # Post-QC plots
    print("Generating post-QC plots")
    pos_cas_var_base64, pos_cas_id_base64, pos_con_var_base64, pos_con_id_base64 = cr_plts(post_mt_cas, post_mt_con,
                                                                                           mind,
                                                                                           geno, "post")
    pos_man_qq_base64 = man_qq_plts(mt, "post")

    # # pre_cas_var_base64, pre_cas_id_base64, pre_con_var_base64, pre_con_id_base64
    qc_plots_list = [pre_man_qq_base64, pos_man_qq_base64, pre_con_id_base64, pre_cas_id_base64, pos_con_id_base64, pos_cas_id_base64,
                     f_stat_plot, pre_con_var_base64, pre_cas_var_base64, pos_con_var_base64, pos_cas_var_base64]

    # Tables
    filter_counts_list = [len(pre_geno_cr), n_sample_miss, len(f_stat_out), len(geno_cr), len(monomorphic_snps),
                          len(snps_hwe_controls), len(snps_hwe_cases)]
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
    #ids_sex_violations = ['IDs: Sex violations -excluded- (N-tested)', len(sex_violations)]
    snps_cr_98 = ['SNPs: call rate < 0.980', filter_counts[3]]
    snps_monomorphic = ['SNPs: without valid association p-value (invariant)', filter_counts[4]]
    snps_hwe_con = ['SNPs: HWE-controls < -6', filter_counts[5]]
    snps_hwe_cas = ['SNPs: HWE-cases < -10', filter_counts[6]]
    exlusion_overview = [snps_cr_95, ids_cr, ids_fhet, snps_cr_98, snps_monomorphic, snps_hwe_con,
                         snps_hwe_cas]
    exlusion_overviewdf = pd.DataFrame(exlusion_overview, columns=['Filter', 'N'])
    exlusion_overview_html = exlusion_overviewdf.to_html()

    return size_of_sample_html, exlusion_overview_html