import hail as hl


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
    return mt
