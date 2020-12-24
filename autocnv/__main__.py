import argparse
from autocnv.annotate import AnnotateHelper

CNVKIT_COL_MAP = {
    'Chr': 'chr', 'Start': 'start', 'End': 'end', 'Detect': 'type'
}
CNV_MAP = {
    'Del': 'del', 'Dup': 'dup'
}

if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument(
        '--in', dest='input', required=True,
        help='input file (TSV format), columns: chr, start, end, type, error are required'
    )
    ap.add_argument(
        '--out', dest='output', required=True,
        help='annotated result file path'
    )
    ap.add_argument('--cnvkit', dest='cnvkit', help='use CNVkit result as input', action='store_true', default=False)
    args = ap.parse_args()
    anno = AnnotateHelper()
    if args.cnvkit:
        anno.annotation_file(args.input, args.output, col_map=CNVKIT_COL_MAP, cnv_map=CNV_MAP)
    else:
        anno.annotation_file(args.input, args.output)