import hail as hl
from typing import List, Tuple, Dict


def compute_qc_metrics(mt: hl.MatrixTable) -> hl.MatrixTable:
    """
    Compute per-sample metrics and common variant statistics useful for quality control
    :param mt: Hail MatrixTable
    :return: Hail MatrixTable with variant and sample qc metrics
    """
    mt = hl.variant_qc(mt)
    mt = hl.sample_qc(mt)

    return mt


def stats_split_mt(mt: hl.MatrixTable) -> Tuple[hl.MatrixTable, hl.MatrixTable]:
    mt_cases: hl.MatrixTable = mt.filter_cols(mt.is_case == True)
    mt_controls: hl.MatrixTable = mt.filter_cols(mt.is_case == False)

    return mt_cases, mt_controls


def collect_counts(mt: hl.MatrixTable) -> List[int]:
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
    n_cases: List[str] = mt.filter_cols(mt.is_case == True).s.collect()
    n_controls: List[str] = mt.filter_cols(mt.is_case == False).s.collect()
    n_unknown_pheno: List[str] = mt.filter_cols(hl.is_missing(mt.is_case)).s.collect()

    # 3. Number of SNPs
    n_snps = mt.count_rows()

    counts: List[int] = [len(n_males), len(n_females), len(n_sex_missing), len(n_cases),
                         len(n_controls), len(n_unknown_pheno), n_snps]

    return counts


def filter_var_cr(mt: hl.MatrixTable, geno: float) -> Tuple[hl.MatrixTable, Dict[str, int]]:
    # steps 1 and 6
    mt = compute_qc_metrics(mt)
    geno_cr_remove = mt.filter_rows(mt.variant_qc.call_rate < (1 - geno)).rsid.collect()
    if len(geno_cr_remove) > 0:
        mt = mt.filter_rows(hl.literal(geno_cr_remove).contains(mt['rsid']), keep=False)

    results = {
        'geno_removed': len(geno_cr_remove)
    }

    return mt, results


def filter_sample_cr(mt: hl.MatrixTable, mind: float) -> Tuple[hl.MatrixTable, Dict[[str, int], [str, int]]]:
    # step 2
    mt = compute_qc_metrics(mt)

    sample_miss = mt.filter_cols(mt.sample_qc.call_rate < (1 - mind))
    sample_miss_cases = sample_miss.filter_cols(sample_miss.is_case == True).s.collect()
    sample_miss_controls = sample_miss.filter_cols(sample_miss.is_case == False).s.collect()
    n_sample_miss = len(sample_miss_cases) + len(sample_miss_controls)
    samples_miss = sample_miss_cases + sample_miss_controls
    if n_sample_miss > 0:
        mt = mt.filter_cols(hl.literal(samples_miss).contains(mt['s']), keep=False)

    results = {
        'sample_miss_cases': len(sample_miss_cases),
        'sample_miss_controls': len(sample_miss_controls)
    }

    return mt, results


def filter_sex_check(mt, fhet_y, fhet_x):
    # step 3
    imputed_sex = hl.impute_sex(mt.GT)
    f_stat_out = mt.filter_cols(((imputed_sex[mt.s].f_stat < fhet_x) & (mt.is_female == False) |
                                 (imputed_sex[mt.s].f_stat > fhet_y) & (mt.is_female == True))).s.collect()
    if len(f_stat_out) > 0:
        mt = mt.filter_cols(hl.literal(f_stat_out).contains(mt['s']), keep=False)

    from .test_plots import fstat_plt
    import pandas as pd
    sex_check_plot = fstat_plt(imputed_sex, fhet_y, fhet_x)
    sex_check_table = pd.DataFrame(f_stat_out, columns=['SampleID'])

    results = {
        'sex_check_plot': sex_check_plot,
        'sex_check_table': sex_check_table
    }

    return mt, results


class FilterSex:
    def __init__(self, mt, input_type):
        self.mt = mt
        self.input_type = input_type

    def sex_violations(self):
        # step 4
        imputed_sex = hl.impute_sex(self.mt.GT)
        if self.input_type == "plink":
            # Verify that when sex info is missing value is set to None
            sex_exclude = self.mt.filter_cols(
                (self.mt.is_female != imputed_sex[self.mt.s]) & (self.mt.is_female != None)).s.collect()
        elif self.input_type == "vcf":
            # Verify that when meta file is read in, column formatting is kept
            sex_exclude = self.mt.filter_cols(
                (self.mt.annotations.Sex != imputed_sex[self.mt.s]) & (self.mt.annotations.Sex != None)).s.collect()
        if len(sex_exclude) > 0:
            mt = self.mt.filter_cols(hl.literal(sex_exclude).contains(self.mt['s']), keep=False)

        results = {
            'sex_excluded': len(sex_exclude)
        }

        return mt, results

    def sex_warnings(self):
        # step 5
        if self.input_type == "plink":
            undef_count = self.mt.aggregate_cols(hl.agg.counter(self.mt.is_female == None))
        elif self.input_type == "vcf":
            undef_count = self.mt.aggregate_cols(hl.agg.counter(self.mt.annotations.Sex == None))
            # print("Warning: {} individuals have undefined phenotype/ambiguous genotypes".format(undef_count))

        return undef_count


def filter_invariant_snps(mt: hl.MatrixTable) -> Tuple[hl.MatrixTable, dict[str, int]]:
    # step 8
    mt = compute_qc_metrics(mt)
    mt = mt.annotate_rows(MAC=hl.min(mt.variant_qc.AC))
    monomorphic_snps = mt.filter_rows(mt.MAC == 0).rsid.collect()
    if len(monomorphic_snps) > 0:
        mt = mt.filter_rows(hl.literal(monomorphic_snps).contains(mt['rsid']), keep=False)

    results = {
        'monomorphic_snps': len(monomorphic_snps)
    }

    return mt, results


def filter_maf(mt: hl.MatrixTable, maf: float) -> Tuple[hl.MatrixTable, dict[str, int]]:
    # step 9
    mt = compute_qc_metrics(mt)
    mt = mt.annotate_rows(maf=hl.min(mt.variant_qc.AF))
    maf_removed = mt.filter_rows(mt.maf < maf).rsid.collect()
    if len(maf_removed) > 0:
        mt = mt.filter_rows(hl.literal(maf_removed).contains(mt['rsid']), keep=False)

    results = {
        'maf_removed': len(maf_removed)
    }

    return mt, results


def filter_hwe(mt: hl.MatrixTable, pheno: str = None, hwe_threshold: float = None) -> Tuple[hl.MatrixTable, int]:
    # steps 10 and 11
    mt = compute_qc_metrics(mt)
    if pheno == 'Case':
        if hwe_threshold:
            hwe_thresh = hwe_threshold
        else:
            hwe_thresh = 1e-10
        mt = mt.filter_cols(mt.is_case == True)
    else:
        if hwe_threshold:
            hwe_thresh = hwe_threshold
        else:
            hwe_thresh = 1e-06
        mt = mt.filter_cols(mt.is_case == False)

    snps = mt.filter_rows(mt.variant_qc.p_value_hwe < hwe_thresh).rsid.collect()
    hwe_snps_removed = len(snps)
    if hwe_snps_removed > 0:
        mt = mt.filter_rows(hl.literal(snps).contains(mt['rsid']), keep=False)

    return mt, hwe_snps_removed
