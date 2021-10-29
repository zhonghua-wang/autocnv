import os

DEFAULT_SCORE = {
    'del': {
        '1A': 0, '1B': -0.6,
        '2A': 1, '2B': 0,
        '2C': -1, '2C-1': 0.9, '2C-2': 0,
        '2D': -1, '2D-1': 0, '2D-2': 0.9, '2D-3': 0.3, '2D-4': 0.9,
        '2E': 0, '2F': -1, '2G': 0, '2H': 0.15, '2I': 0, '2J': 0, '2K': 0.45,
        '3A': 0, '3B': 0.45, '3C': 0.9,
        '4O': -1,
        'PVS1': 0.9, 'PVS1_S': 0.45, 'PVS1_M': 0.3, 'PVS1_P': 0.15, 'PVS1_U': 0
    },
    'dup': {
        '1A': 0, '1B': -0.6,
        '2A': 1, '2B': 0,
        '2C': -1, '2C-1': 0.9, '2C-2': 0,
        '2D': -1, '2D-1': 0, '2D-2': 0.9, '2D-3': 0.3, '2D-4': 0.9,
        '2E': 0, '2F': -1, '2G': 0, '2H': 0, '2I': 0, '2J': 0, '2K': 0.45, '2L': 0,
        '3A': 0, '3B': 0.45, '3C': 0.9,
        '4O': -1,
        'PVS1': 0.9, 'PVS1_S': 0.45, 'PVS1_M': 0.3, 'PVS1_P': 0.15, 'PVS1_U': 0
    }
}



BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

CYTO_BAND_FILE = os.path.join(BASE_DIR, 'data', 'cyto-band.bed.gz')

GENE_EXON_DATABASE = os.path.join(BASE_DIR, 'data', 'exon.sorted.bed.gz')

GENE_DATABASE = os.path.join(BASE_DIR, 'data', 'gene.sorted.bed.gz')

OMIM_GENE_DATABASE = os.path.join(BASE_DIR, 'data', 'omim-gene.sorted.bed.gz')

FUNC_REGION_DATABASE = os.path.join(BASE_DIR, 'data', 'func-region.sorted.gz')

HI_GENE_DATABASE = os.path.join(BASE_DIR, 'data', 'hi-gene.sorted.bed.gz')

HI_EXON_DATABASE = os.path.join(BASE_DIR, 'data', 'hi-exon.sorted.bed.gz')

HI_CDS_DATABASE = os.path.join(BASE_DIR, 'data', 'hi-cds.sorted.bed.gz')

CLINVAR_PATHOGENIC_DATABASE = os.path.join(
    BASE_DIR, 'data', 'clinvar-pathogenic.sorted.vcf.gz'
)

UHI_GENE_DATABASE = os.path.join(BASE_DIR, 'data', 'uhi-gene.sorted.bed.gz')

HI_REGION_DATABASE = os.path.join(BASE_DIR, 'data', 'hi-region.sorted.bed.gz')

UHI_REGION_DATABASE = os.path.join(BASE_DIR, 'data', 'uhi-region.sorted.bed.gz')

DECIPHER_GENE_DATABASE = os.path.join(BASE_DIR, 'data', 'decipher-gene.sorted.bed.gz')

TS_GENE_DATABASE = os.path.join(BASE_DIR, 'data', 'ts-gene.sorted.bed.gz')

TS_REGION_DATABASE = os.path.join(BASE_DIR, 'data', 'ts-region.sorted.bed.gz')

UTS_GENE_DATABASE = os.path.join(BASE_DIR, 'data', 'uts-gene.sorted.bed.gz')

UTS_REGION_DATABASE = os.path.join(BASE_DIR, 'data', 'uts-region.sorted.bed.gz')

DGV_GAIN_DATABASE = os.path.join(BASE_DIR, 'data', 'dgv-gain.sorted.bed.gz')

DGV_LOSS_DATABASE = os.path.join(BASE_DIR, 'data', 'dgv-loss.sorted.bed.gz')

GNOMAD_DEL_DATABASE = os.path.join(BASE_DIR, 'data', 'gnomad-del.sorted.bed.gz')

GNOMAD_DUP_DATABASE = os.path.join(BASE_DIR, 'data', 'gnomad-dup.sorted.bed.gz')

CNV_SYNDROME_DEL_DATABASE = os.path.join(BASE_DIR, 'data', 'cnv-syndrome-del.bed.gz')
CNV_SYNDROME_DUP_DATABASE = os.path.join(BASE_DIR, 'data', 'cnv-syndrome-dup.bed.gz')


try:
    from autocnv.local_settings import *
except ImportError:
    pass
