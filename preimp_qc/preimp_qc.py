#!/usr/bin/env python

from bokeh.io import export_png
from bokeh.layouts import row, column

import argparse
from .functions import *
from .io import vcf_to_mt, plink_to_mt
from .utils import TimeLogger
from .report import Report


def initial_summary_stats(mt: hl.MatrixTable, artifacts_dir: str,
                          report: Report = None, verbose: bool = True) -> hl.MatrixTable:
    with TimeLogger('initial summary stats', verbose=verbose):
        mt, init_stats = summary_stats(mt)
        mt, p_values = compute_p_values(mt, title='Pre-QC')

        is_case_counts = init_stats['is_case_counts']
        is_female_counts = init_stats['is_female_counts']

        manqq_path = f'{artifacts_dir}/preqc/manqq.png'
        fig = row(p_values['qq_plot'],
                  p_values['manhattan_plot'])
        export_png(fig, filename=manqq_path)

        if report:
            with report.add_section() as s:
                s.add_header('Initial Summary Statistics')
                s.add_text('No QC has been done.')

                with s.add_section() as ss:
                    ss.add_header('Counts')
                    with ss.add_section() as sss:
                        sss.add_header('Phenotype Counts')
                        with sss.add_table() as t:
                            t.add_header_row(['Phenotype', 'Count'])
                            t.add_row(['Case', init_stats['case']])
                            t.add_row(['Control', is_case_counts['control']])
                            t.add_row(['Unknown', is_case_counts['unknown']])

                    with ss.add_section() as sss:
                        sss.add_header('Sex Counts')
                        with sss.add_table() as t:
                            t.add_header_row(['Sex', 'Count'])
                            t.add_row(['Female', is_female_counts['female']])
                            t.add_row(['Male', is_female_counts['male']])
                            t.add_row(['Unknown', is_female_counts['unknown']])

                    with sss.add_section() as sss:
                        sss.add_header('Data Counts')
                        with sss.add_table() as t:
                            t.add_header_row(['Variable', 'Count'])
                            t.add_row(['# of Samples', init_stats['n_samples']])
                            t.add_row(['# of Variants', init_stats['n_variants']])
                            t.add_row(['# of GWS Variants', p_values['n_sig_variants']])

                with s.add_section() as ss:
                    ss.add_header('Manhattan QQ-Plot')
                    ss.add_image(manqq_path)

    return mt


def variant_qc_1(mt: hl.MatrixTable, call_rate: float, artifacts_dir: str, *,
                 report: Report = None, verbose: bool = True) -> hl.MatrixTable:
    with TimeLogger(f'first pass variant qc (call_rate={call_rate})', verbose=verbose):
        mt, results = filter_variants_by_call_rate(mt, call_rate)

        # Not sure if we want to do this, but an example in case
        excluded_variants = results['excluded_variants']
        excluded_variants.export(f'{artifacts_dir}/varqc-pass-1/excluded-variants.tsv')

        export_png(results['call_rate_plot'], filename=f'{artifacts_dir}/varqc-pass-1/call-rate.png')

        if report:
            pass

    return mt


# TODO: I didn't spend a lot of time on this. You should assume all of this is wrong and
# it should be broken up into functions that go into functions.py
def sample_qc(mt: hl.MatrixTable, call_rate: float, artifacts_dir: str, *,
              report: Report = None, verbose: bool = True) -> hl.MatrixTable:
    with TimeLogger(f'sample qc (call_rate={call_rate})', verbose=verbose):
        mt, results = filter_samples_by_call_rate(mt, call_rate)
        sample_miss = mt.filter_cols(mt.sample_qc.call_rate < call_rate)
        sample_miss_cases = sample_miss.filter_cols(sample_miss.is_case == True).s.collect()
        sample_miss_controls = sample_miss.filter_cols(sample_miss.is_case == False).s.collect()
        n_sample_miss = len(sample_miss_cases) + len(sample_miss_controls)
        samples_miss = sample_miss_cases + sample_miss_controls
        if n_sample_miss > 0:
            mt = mt.filter_cols(hl.literal(samples_miss).contains(mt['s']), keep=False)

        if report:
            pass

    return mt


# TODO: I didn't spend a lot of time on this. You should assume all of this is wrong and
# it should be broken up into functions that go into functions.py
def variant_qc_2(mt: hl.MatrixTable, call_rate: float, hwe_case: float, hwe_control: float,
                 maf: float, call_rate_diff: float, artifacts_dir: str, *,
                 report: Report = None, verbose: bool = True) -> hl.MatrixTable:
    with TimeLogger(f'second pass variant qc (call_rate={call_rate})', verbose=verbose):
        mt, call_rate_results = filter_variants_by_call_rate(mt, call_rate)
        mt, monomorphic_results = filter_monomorphic_variants(mt)

        mts_by_pheno = split_mt(mt.is_case)
        mt, hwe_case_results = filter_variants_by_hwe(mt['case'], hwe_case)
        mt, hwe_control_results = filter_variants_by_hwe(mt['control'], hwe_control)

        mt, maf_results = filter_variants_by_maf(mt, maf)
        mt, call_rate_by_pheno_results = filter_variants_by_call_rate_diff(mt, call_rate_diff)

        # Not sure if we want to do this, but I gave an example in case
        excluded_variants = call_rate_results['excluded_variants']
        excluded_variants.export(f'{artifacts_dir}/varqc-pass-1/call-rate-excluded-variants.tsv')

        export_png(call_rate_results['call_rate_plot'], filename=f'{artifacts_dir}/varqc-pass-2/call-rate.png')
        export_png(hwe_case_results['hwe_plot'], filename=f'{artifacts_dir}/varqc-pass-2/hwe-case.png')
        export_png(hwe_control_results['hwe_plot'], filename=f'{artifacts_dir}/varqc-pass-2/hwe-control.png')
        export_png(maf_results['maf_plot'], filename=f'{artifacts_dir}/varqc-pass-2/maf.png')

        # TODO: more operations

        if report:
            pass

    return mt


def preimp_qc(mt, artifacts_dir, var_call_rate_1=0.95, sample_call_rate=0.98,
              var_call_rate_2=0.98, hwe_case=1E-6, hwe_control=1E-6, maf=0.01, report=None):
    mt = variant_qc_1(mt, var_call_rate_1, artifacts_dir, report=report)

    mt = sample_qc(mt, sample_call_rate, artifacts_dir, report=report)

    mt = variant_qc_2(mt, call_rate=var_call_rate_2, hwe_case=hwe_case,
                      hwe_control=hwe_control, maf=maf, call_rate_diff=0.02,
                      artifacts_dir=artifacts_dir, report=report)

    return mt


if __name__ == '__main__':
    # QUESTION: Do you want to keep the argument names here or go for a clearer
    # naming scheme than Ricopili had? Use `-` instead of `_` regardless as that's the standard.

    parser = argparse.ArgumentParser(description='preimp_qc V1.0')
    parser.add_argument('--artifacts', type=str, required=True)
    parser.add_argument('--dest', type=str, required=True)
    parser.add_argument('--input-type', type=str, required=True, choices=['vcf', 'plink'])
    parser.add_argument('--annotations', type=str)
    parser.add_argument('--reference', type=str, default='GRCh38')
    parser.add_argument('--qc-round', type=str, required=True)

    # required for QC
    parser.add_argument('--pre_geno', type=float, default=0.05,
                        help="include only SNPs with missing-rate < NUM (before "
                             "ID filter), important for post merge of multiple "
                             "platforms")
    parser.add_argument('--mind', type=float, default=0.02, help="include only IDs with missing-rate < NUM")
    parser.add_argument('--fhet-y', type=float, default=0.4, help="include only female IDs with fhet < NUM")
    parser.add_argument('--fhet-x', type=float, default=0.8, help="include only male IDs with fhet > NUM")
    parser.add_argument('--geno', type=float, default=0.02, help="include only SNPs with missing-rate < NUM")
    parser.add_argument('--midi', type=float, default=0.02, help="include only SNPs with missing-rate-difference ("
                                                                 "case/control) < NUM")
    parser.add_argument('--with-pna', type=int, default=0, choices=[0, 1], help="include monomorphic (invariant) SNPs")
    parser.add_argument('--maf', type=float, default=0.01, help="include only SNPs with MAF >= NUM")
    parser.add_argument('--hwe-th-con', type=float, default=1e-6, help="HWE_controls < NUM")
    parser.add_argument('--hwe-th-cas', type=float, default=1e-6, help="HWE_cases < NUM")

    args = parser.parse_args()

    # read input
    if args.inputType == 'plink':
        input_mt = plink_to_mt(args.dirname, args.basename, args.reference)
    else:
        assert args.input_type == 'vcf'
        input_mt = vcf_to_mt(args.dirname, args.vcf, args.annotations)

    # Note: All functions in preimp_qc assume there exists a locus, s, is_case, and is_female field
    # we might want to assert this is the case before calling preimp_qc so it's clear

    report = Report(f'{args.artifacts}/report.md')

    # TODO: fill in parameters
    mt = preimp_qc(input_mt, artifacts_dir=args.artifacts, report=report)

    mt.write(args.dest)
    report.write()
