import pandas as pd
import hail as hl
from typing import Tuple, Any, Dict, Optional

from .utils import gen_uid


#Tested and working
# Usage example (split mt by phenotype): split_mts = split_mt(mt.is_case), then to get the cases mt, cases_mt = split_mts[True]
def split_mt(category: hl.Expression) -> Dict[str, hl.MatrixTable]:
    mt = category._indices.source
    categories = mt.aggregate_cols(hl.struct(categories=hl.agg.collect_as_set(category)))

    result = {}
    for cat in categories['categories']:
        if cat is not None:
            result[cat] = mt.filter_cols(category == cat)
        else:
            result[cat] = mt.filter_cols(hl.is_missing(category))
    return result


# Tested and working
def summary_stats(mt: hl.MatrixTable) -> Tuple[hl.MatrixTable, Dict[str, Any]]:
    results = {}

    counts = mt.aggregate_cols(hl.struct(is_case=hl.agg.counter(mt.is_case),
                                         is_female=hl.agg.counter(mt.is_female)))

    mapping = {
        True: 'case',
        False: 'control',
        None: 'unknown'
    }
    is_case_counts = {mapping[pheno]: count for pheno, count in counts['is_case'].items()}
    results['is_case_counts'] = is_case_counts

    mapping = {
        True: 'female',
        False: 'male',
        None: 'unknown'
    }
    is_female_counts = {mapping[pheno]: count for pheno, count in counts['is_female'].items()}
    results['is_female_counts'] = is_female_counts

    n_variants = mt.count_rows()
    n_samples = mt.count_cols()

    results['n_variants'] = n_variants
    results['n_samples'] = n_samples

    return mt, results


def compute_p_values(mt, covariates=None, title: Optional[str] = None) -> Tuple[hl.MatrixTable, Dict[str, Any]]:
    if covariates is None:
        covariates = [1.0]

    results = {}

    gwas = hl.linear_regression_rows(y=mt.is_case,
                                     x=mt.GT.n_alt_alleles(),
                                     covariates=covariates)

    n_sig_variants = gwas.filter(gwas.p_value < 5E-8).count()
    results['n_sig_variants'] = n_sig_variants

    title = f': {title}' if title else ''
    qq_plot = hl.plot.qq(gwas.p_value, title=f'Q-Q Plot{title}')
    manhattan_plot = hl.plot.manhattan(gwas.p_value, title=f'Manhattan Plot{title}')

    results['qq_plot'] = qq_plot
    results['manhattan_plot'] = manhattan_plot

    return mt, results


def filter_variants_by_call_rate(mt: hl.MatrixTable, call_rate: float, *,
                                 title: Optional[str] = None) -> Tuple[hl.MatrixTable, Dict[str, Any]]:
    qc_dest = f'variant-qc-{gen_uid()}'
    mt = hl.variant_qc(mt, name=qc_dest)

    plot = hl.plot.histogram(mt[qc_dest].call_rate, title=title, legend='Call Rate')
    # FIXME: add vertical line for threshold

    excluded_variants = mt.filter_rows(mt[qc_dest].call_rate < call_rate).locus
    n_excluded_variants = hl.eval(excluded_variants.length())

    mt = mt.filter_rows(mt[qc_dest].call_rate >= call_rate, keep=True)

    mt.drop(qc_dest)

    results = {
        'excluded_variants': excluded_variants,
        'n_excluded_variants': n_excluded_variants,
        'call_rate': call_rate,
        'call_rate_plot': plot
    }

    return mt, results


def filter_monomorphic_variants(mt: hl.MatrixTable) -> Tuple[hl.MatrixTable, Dict[str, Any]]:
    qc_dest = f'variant-qc-{gen_uid()}'
    mt = hl.variant_qc(mt, name=qc_dest)

    excluded_variants = mt.filter_rows(hl.any(mt[qc_dest].AF == 1), keep=True).locus.collect()
    if len(excluded_variants) > 0:
        mt = mt.filter_rows(hl.literal(excluded_variants).contains(mt['locus']), keep=False)

    mt.drop(qc_dest)

    results = {
        'excluded_variants': excluded_variants,
        'n_excluded_variants': len(excluded_variants)
    }

    return mt, results


def filter_variants_by_maf(mt: hl.MatrixTable, maf: float,
                           title: Optional[str] = None) -> Tuple[hl.MatrixTable, Dict[str, Any]]:
    qc_dest = f'variant-qc-{gen_uid()}'
    mt = hl.variant_qc(mt, name=qc_dest)

    plot = hl.plot.histogram(hl.min(mt[qc_dest].AF), legend='Minor Allele Frequency',
                             log=True, title=title)
    # FIXME: add vertical line for threshold

    excluded_variants = mt.filter_rows(hl.min(mt[qc_dest].AF) < maf).locus.collect()
    if len(excluded_variants) > 0:
        mt = mt.filter_rows(hl.literal(excluded_variants).contains(mt['locus']), keep=False)

    mt.drop(qc_dest)

    results = {
        'excluded_variants': excluded_variants,
        'n_excluded_variants': len(excluded_variants),
        'maf': maf,
        'plot': plot
    }

    return mt, results


def filter_variants_by_hwe(mt: hl.MatrixTable, hwe: float,
                           title: Optional[str] = None) -> Tuple[hl.MatrixTable, Dict[str, Any]]:
    qc_dest = f'variant-qc-{gen_uid()}'
    mt = hl.variant_qc(mt, name=qc_dest)

    plot = hl.plot.histogram(mt[qc_dest].p_value_hwe, legend='HWE P-value (log10)',
                             log=True, title=title)
    # FIXME: add vertical line for threshold

    excluded_variants = mt.filter_rows(mt[qc_dest].p_value_hwe < hwe).locus.collect()
    if len(excluded_variants) > 0:
        mt = mt.filter_rows(hl.literal(excluded_variants).contains(mt['locus']), keep=False)

    mt.drop(qc_dest)

    results = {
        'excluded_variants': excluded_variants,
        'n_excluded_variants': len(excluded_variants),
        'hwe': hwe,
        'hwe_plot': plot
    }

    return mt, results


def filter_variants_by_call_rate_diff(mt: hl.MatrixTable, call_rate_diff: float) -> Tuple[hl.MatrixTable, Dict[str, Any]]:
    raise NotImplementedError


def filter_samples_by_call_rate(mt: hl.MatrixTable, call_rate: float, *,
                                title: Optional[str] = None) -> Tuple[hl.MatrixTable, Dict[str, Any]]:
    qc_dest = f'sample-qc-{gen_uid()}'
    mt = hl.sample_qc(mt, name=qc_dest)

    plot = hl.plot.histogram(mt[qc_dest].call_rate, title=title)
    # FIXME: add vertical line for threshold

    excluded_samples = mt.filter_cols(mt[qc_dest].call_rate < call_rate).s
    n_excluded_samples = hl.eval(excluded_samples.length())

    mt = mt.filter_cols(mt[qc_dest].call_rate >= call_rate, keep=True)

    mt.drop(qc_dest)

    results = {
        'excluded_samples': excluded_samples,
        'n_excluded_samples': n_excluded_samples,
        'call_rate': call_rate,
        'call_rate_plot': plot
    }

    return mt, results


def filter_samples_by_sex_violations(gt: hl.CallExpression, is_female: hl.BooleanExpression,
                                     f_stat_y: float, f_stat_x: float) -> Tuple[hl.MatrixTable, Dict[str, Any]]:
    mt_gt = gt._indices.source
    mt_is_female = is_female._indices.source
    assert mt_gt == mt_is_female
    mt = mt_gt

    imputed_sex_ht = hl.impute_sex(gt)
    f_stat_out = mt.filter_cols(((imputed_sex_ht[mt.s].f_stat < f_stat_x) & (is_female == False) |
                                 (imputed_sex_ht[mt.s].f_stat > f_stat_y) & (is_female == True))).s.collect()
    if len(f_stat_out) > 0:
        mt = mt.filter_cols(hl.literal(f_stat_out).contains(mt['s']), keep=False)
        
    excluded_samples = f_stat_out
    n_excluded_samples = len(f_stat_out)

    # FIXME: I wasn't sure what you wanted to use this for
    # wanted to have an option to write out these samples to a file, but we can use the excluded_samples variable for that
    # sex_check_table = pd.DataFrame(f_stat_out, columns=['SampleID'])

    results = {
        'excluded_samples': excluded_samples,
        'n_excluded_samples': n_excluded_samples,
        #'sex_check_table': sex_check_table,
        'f_stat_x': f_stat_x,
        'f_stat_y': f_stat_y
    }

    return mt, results
