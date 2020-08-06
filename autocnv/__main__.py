import argparse
from autocnv.annotate import AnnotateHelper


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
    args = ap.parse_args()
    anno = AnnotateHelper()
    anno.annotation_file(args.input, args.output)