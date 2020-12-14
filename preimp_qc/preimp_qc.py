#!/usr/bin/env python

import argparse
from preimp_qc.functions import compute_qc_metrics, plink_to_mt, vcf_to_mt, run_qc
from preimp_qc.report import write_html_report


def main():
    parser = argparse.ArgumentParser()
    parser = argparse.ArgumentParser(description='preimp_qc V1.0')
    parser.add_argument('--dirname', type=str, required=True)
    parser.add_argument('--basename', type=str, required=True)
    parser.add_argument('--inputType', type=str, required=True, choices=['vcf', 'plink'])
    parser.add_argument('--annotations', type=str)
    parser.add_argument('--reference', type=str, default='GRCh38')
    parser.add_argument('--qc_round', type=str, required=True)

    # required for QC
    parser.add_argument('--pre_geno', type=float, default=0.05,
                        help="include only SNPs with missing-rate < NUM (before "
                             "ID filter), important for post merge of multiple "
                             "platforms")
    parser.add_argument('--mind', type=float, default=0.02, help="include only IDs with missing-rate < NUM")
    parser.add_argument('--fhet_y', type=float, default=0.4, help="include only female IDs with fhet < NUM")
    parser.add_argument('--fhet_x', type=float, default=0.8, help="include only male IDs with fhet > NUM")
    parser.add_argument('--geno', type=float, default=0.02, help="include only SNPs with missing-rate < NUM")
    parser.add_argument('--midi', type=float, default=0.02, help="include only SNPs with missing-rate-difference ("
                                                                 "case/control) < NUM")
    parser.add_argument('--withpna', type=int, default=0, choices=[0, 1], help="include monomorphic (invariant) SNPs")
    parser.add_argument('--maf', type=float, default=0.01, help="include only SNPs with MAF >= NUM")
    parser.add_argument('--hwe_th_con', type=float, default=1e-6, help="HWE_controls < NUM")
    parser.add_argument('--hwe_th_cas', type=float, default=1e-6, help="HWE_cases < NUM")

    arg = parser.parse_args()

    # read input
    if arg.inputType == 'plink':
        input_mt = plink_to_mt(arg.dirname, arg.basename, arg.reference)

    if arg.inputType == 'vcf':
        input_mt = vcf_to_mt(arg.dirname, arg.vcf, arg.annotations)

    # compute sample and variant qc metrics
    print("Computing QC metrics")
    mt = compute_qc_metrics(input_mt)

    print("Running QC")
    qc_tables, qc_plots = run_qc(mt, arg.dirname, arg.basename, arg.pre_geno, arg.mind, arg.fhet_y, arg.fhet_x,
                                 arg.geno, arg.midi, arg.maf, arg.hwe_th_con, arg.hwe_th_cas, arg.qc_round, arg.withpna)

    print("Generating report")
    write_html_report(arg.dirname, arg.basename, qc_tables, qc_plots)

    print("\nDone running QC!")


if __name__ == '__main__':
    main()
