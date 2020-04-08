from acit.database import DataBase
from pysam import VariantFile
from acit import settings
from autopvs1.cnv import CNVRecord, PVS1CNV
from autopvs1.utils import get_transcript
from autopvs1.read_data import transcripts
from autopvs1.strength import Strength
from collections import defaultdict
from itertools import chain
import operator
import pandas as pd
from acit.utils import ACITEncoder
import json
from collections import OrderedDict


SEP = '\n'
DEFAULT_EMPTY_VALUE = '-'

PVS1 = {
    Strength.VeryStrong: 'PVS1', Strength.Strong: 'PVS1_S', Strength.Moderate: 'PVS1_M',
    Strength.Supporting: 'PVS1_P', Strength.Unmet: 'PVS1_U'
}

SCORE_GROUP = {
    'del': {
        rule: 'G1' for rule in ('2A', '2B', '2C-1', '2C-2', '2D-1', '2D-2', '2D-3', '2D-4', '2E')
    },
    'dup': {}
}

PATHOGENICITY_LEVELS = [
    (operator.ge, 0.99, 'P'), (operator.ge, 0.9, 'LP'), (operator.gt, -0.9, 'VUS'),
    (operator.gt, -0.99, 'LB'), (operator.le, -0.99, 'B')
]


class AnnotateHelper:
    def __init__(self):
        self._gene_database = DataBase(settings.GENE_DATABASE)
        self._omim_gene_database = DataBase(settings.OMIM_GENE_DATABASE)
        self._func_region_database = DataBase(settings.FUNC_REGION_DATABASE)
        self._hi_gene_database = DataBase(settings.HI_GENE_DATABASE)
        self._hi_exon_database = DataBase(settings.HI_EXON_DATABASE)
        self._hi_cds_database = DataBase(settings.HI_CDS_DATABASE)
        self._clinvar_pathogenic_database = VariantFile(settings.CLINVAR_PATHOGENIC_DATABASE)
        self._uhi_gene_database = DataBase(settings.UHI_GENE_DATABASE)
        self._hi_region_database = DataBase(settings.HI_REGION_DATABASE)
        self._uhi_region_database = DataBase(settings.UHI_REGION_DATABASE)
        self._decipher_gene_database = DataBase(settings.DECIPHER_GENE_DATABASE)
        self._ts_gene_database = DataBase(settings.TS_GENE_DATABASE)
        self._ts_region_database = DataBase(settings.TS_REGION_DATABASE)
        self._uts_gene_database = DataBase(settings.UTS_GENE_DATABASE)
        self._uts_region_database = DataBase(settings.UTS_REGION_DATABASE)
        self._dgv_gain_database = DataBase(settings.DGV_GAIN_DATABASE)
        self._dgv_loss_database = DataBase(settings.DGV_LOSS_DATABASE)
        self._gnomad_del_database = DataBase(settings.GNOMAD_DEL_DATABASE)
        self._gnomad_dup_database = DataBase(settings.GNOMAD_DUP_DATABASE)

    @staticmethod
    def _norm_chrom(ch):
        """
        normalize chromosome name, eg. 2 -> chr2, 23 -> chrX
        :param ch: input chromosome name
        :return: normalized name
        >>> norm_chrom(2)
        'chr2'
        >>> norm_chrom('chr23')
        'chrX'
        """
        ch = str(ch).replace('chr', '')
        if ch == '23':
            return 'chrX'
        if ch == '24':
            return 'chrY'
        return f'chr{ch}'

    @staticmethod
    def _annotate_loss(**annotation):
        loss = dict()

        if len(annotation['overlap_genes']) + len(annotation['overlap_func_regions']) > 0:
            loss['1A'] = True
        else:
            loss['1B'] = True

        for gene, overlap, coverage in annotation['overlap_hi_genes']:
            # 是否覆盖整个基因
            if coverage == 1:
                loss['2A'] = True
            elif overlap < 1:  # 是否位于基因内部
                if any(
                        exon.last_exon == 'True'
                        for exon, *_ in annotation['overlap_hi_exons'][gene.gene_id]
                ):  # 是否覆盖末位外显子
                    if len(annotation['overlap_hi_exons'][gene.gene_id]) >= 2:
                        # 覆盖超过两个外显子
                        loss['2D-4'] = True
                    elif gene.gene_id in annotation['overlap_hi_cds'] \
                            and len(annotation['overlap_hi_cds'][gene.gene_id]) > 0:  # 是否覆盖CDS
                        if len(annotation['variants']) > 0:  # 末位外显子是否有致病变异
                            loss['2D-2'] = True
                        else:
                            loss['2D-3'] = True
                    else:
                        # 不覆盖CDS区
                        loss['2D-1'] = True
                # 未覆盖末位外显子
                elif gene.gene_id in annotation['overlap_hi_cds'] \
                        and len(annotation['overlap_hi_cds'][gene.gene_id]) > 0:  # 是否覆盖5'端CDS
                    loss['2C-1'] = True
                else:
                    loss['2C-2'] = True
            # 位于基因内部
            else:
                cnv = CNVRecord(
                    annotation['chromosome'], annotation['outer_start'],
                    annotation['outer_end'], annotation['func']
                )
                tx = get_transcript(gene.transcript, transcripts)
                pvs1 = PVS1CNV(cnv, None, tx)
                loss['2E'] = True
                loss[PVS1[pvs1.verify_DEL()[0]]] = True

        # 完全覆盖hi区域
        for region, overlap, coverage in annotation['overlap_hi_regions']:
            if coverage == 1:
                loss['2A'] = True
            elif set(region.omim_genes.split(',')) == \
                    set(gene.symbol for gene, *_ in annotation['overlap_omim_genes']):
                # 是否覆盖区域内所有OMIM致病基因
                loss['2A'] = True
            else:
                loss['2B'] = True

        # 包含预测HI基因
        if len(annotation['overlap_hi_genes']) + len(annotation['overlap_hi_regions']) == 0 \
                and len(annotation['overlap_decipher_genes']) > 0:
            loss['2H'] = True

        # 落入uhi基因
        for gene, overlap, coverage in annotation['overlap_uhi_genes']:
            if overlap == 1:
                loss['2F'] = True

        # 落入uhi区域
        for region, overlap, coverage in annotation['overlap_uhi_regions']:
            genes = set(gene.symbol for gene, *_ in annotation['overlap_genes'])
            if overlap == 1:  # 完全落入区域
                loss['2F'] = True
            elif set(region.genes.split(',')) == genes:  # 覆盖相同的基因
                loss['2F'] = True

        # 覆盖基因个数
        gene_count = len(annotation['overlap_genes'])
        if gene_count >= 35:
            loss['3C'] = True
        elif gene_count >= 25:
            loss['3B'] = True
        elif gene_count >= 0:
            loss['3A'] = True

        # DGV金标
        genes = set(gene.symbol for gene, *_ in annotation['overlap_genes'])
        for record, *_ in chain(annotation['dgv_loss_records'], annotation['gnomad_del_records']):
            if len(genes - set(record.genes.split(','))) == 0:
                loss['4O'] = True

        annotation['rules'] = loss
        return annotation

    @staticmethod
    def _annotate_gain(**annotation):
        gain = dict()

        if len(annotation['overlap_genes']) + len(annotation['overlap_func_regions']) > 0:
            gain['1A'] = True
        else:
            gain['1B'] = True

        for gene, overlap, coverage in annotation['overlap_ts_genes']:
            # 覆盖整个基因
            if coverage == 1:
                gain['2A'] = True

        # 完全覆盖ts区域
        for region, overlap, coverage in annotation['overlap_ts_regions']:
            if coverage == 1:  # 是否覆盖整改区域
                gain['2A'] = True
            elif set(region.omim_genes.split(',')) == \
                    set(gene.symbol for gene, *_ in annotation['overlap_omim_genes']):
                # 是否覆盖区域内所有OMIM致病基因
                gain['2A'] = True
            else:
                gain['2B'] = True

        # 落入uts基因
        for gene, overlap, coverage in annotation['overlap_uts_genes']:
            if overlap == 1:
                gain['2D'] = True

        # 落入uts区域
        for region, overlap, coverage in annotation['overlap_uts_regions']:
            genes = set(gene.symbol for gene, *_ in annotation['overlap_genes'])
            region_genes = set(region.genes.split(','))
            if overlap == coverage == 1:
                gain['2C'] = True
            elif len(genes - region_genes) > 0:  # 多
                gain['2G'] = True
            elif any(c < 1 for *_, c in annotation['overlap_genes']):  # 破坏蛋白编码基因
                gain['2E'] = True
            elif overlap == 1:
                gain['2D'] = True
            else:
                gain['2F'] = True

        hi_genes = set()
        for gene, overlap, coverage in annotation['overlap_hi_genes']:
            hi_genes.add(gene.symbol)
            if coverage == 1:
                gain['2H'] = True
            elif overlap == 1:
                cnv = CNVRecord(
                    annotation['chromosome'], annotation['outer_start'],
                    annotation['outer_end'], annotation['func']
                )
                tx = get_transcript(gene.transcript, transcripts)
                pvs1 = PVS1CNV(cnv, None, tx)
                gain['2I'] = True
                gain[PVS1[pvs1.verify_DUP()[0]]] = True

        for gene, overlap, coverage in annotation['overlap_genes']:
            if gene.symbol not in hi_genes and coverage != 1:
                gain['2L'] = True
                annotation['break_point_genes'].append(gene.symbol)

        # 覆盖基因个数
        gene_count = len(annotation['overlap_genes'])
        if gene_count >= 50:
            gain['3C'] = True
        elif gene_count >= 35:
            gain['3B'] = True
        elif gene_count >= 0:
            gain['3A'] = True

        # DGV金标
        genes = set(gene.symbol for gene, *_ in annotation['overlap_genes'])
        for record, *_ in chain(annotation['dgv_gain_records'], annotation['gnomad_dup_records']):
            if len(genes - set(record.genes.split(','))) == 0:
                gain['4O'] = True

        annotation['rules'] = gain
        return annotation

    @staticmethod
    def merge_score(func, **rules):
        groups = defaultdict(list)
        for rule, score in rules.items():
            try:
                groups[SCORE_GROUP[func][rule]].append(score)
            except KeyError:
                yield score
        for _, scores in groups.items():
            yield max(scores)

    @staticmethod
    def judge(func, **rules):
        rules = {
            rule: settings.DEFAULT_SCORE[func][rule] for rule, check in rules.items() if check
        }
        score = sum(AnnotateHelper.merge_score(func, **rules))
        for op, cutoff, level in PATHOGENICITY_LEVELS[:-1]:
            if op(score, cutoff):
                pathogenicity = level
                break
        else:
            pathogenicity = PATHOGENICITY_LEVELS[-1][2]
        return rules, score, pathogenicity

    def annotate(self, chromosome, start, end, func, error=0):
        annotation = dict(
            chromosome=chromosome, start=start, end=end,
            length=end - start, error=error,
            outer_start=start - error, outer_end=end + error,
            inner_start=start + error, inner_end=end - error,
            func=func, break_point_genes=list()
        )

        annotation['overlap_genes'] = list(self._gene_database.overlap(
            chromosome, annotation['outer_start'], annotation['outer_end']
        ))

        annotation['overlap_omim_genes'] = list(self._omim_gene_database.overlap(
            chromosome, annotation['outer_start'], annotation['outer_end']
        ))

        annotation['overlap_func_regions'] = list(self._func_region_database.overlap(
            chromosome, annotation['outer_start'], annotation['outer_end']
        ))

        annotation['overlap_hi_genes'] = list(self._hi_gene_database.overlap(
            chromosome, annotation['inner_start'], annotation['inner_end']
        ))

        annotation['overlap_hi_exons'] = self._hi_exon_database.overlap_groups(
            chromosome, annotation['inner_start'], annotation['inner_end'],
            lambda record: record[0].gene_id
        )

        annotation['overlap_hi_cds'] = self._hi_cds_database.overlap_groups(
            chromosome, annotation['inner_start'], annotation['inner_end'],
            lambda record: record[0].gene_id
        )

        try:
            annotation['variants'] = list(self._clinvar_pathogenic_database.fetch(
                chromosome, annotation['inner_start'], annotation['inner_end'])
            )
        except ValueError:
            annotation['variants'] = []

        annotation['overlap_hi_regions'] = list(self._hi_region_database.overlap(
            chromosome, annotation['inner_start'], annotation['inner_end']
        ))

        annotation['overlap_decipher_genes'] = list(self._decipher_gene_database.overlap(
            chromosome, annotation['inner_start'], annotation['inner_end']
        ))

        annotation['overlap_uhi_genes'] = list(self._uhi_gene_database.overlap(
            chromosome, annotation['outer_start'], annotation['outer_end']
        ))

        annotation['overlap_uhi_regions'] = list(self._uhi_region_database.overlap(
            chromosome, annotation['outer_start'], annotation['outer_end']
        ))

        annotation['overlap_ts_genes'] = list(self._ts_gene_database.overlap(
            chromosome, annotation['inner_start'], annotation['inner_end']
        ))

        annotation['overlap_ts_regions'] = list(self._ts_region_database.overlap(
            chromosome, annotation['inner_start'], annotation['inner_end']
        ))

        annotation['overlap_uts_genes'] = list(self._uts_gene_database.overlap(
            chromosome, annotation['outer_start'], annotation['outer_end']
        ))

        annotation['overlap_uts_regions'] = list(self._uts_region_database.overlap(
            chromosome, annotation['outer_start'], annotation['outer_end']
        ))

        annotation['dgv_gain_records'] = list(self._dgv_gain_database.overlap(
            chromosome, annotation['outer_start'], annotation['outer_end']
        ))

        annotation['dgv_loss_records'] = list(self._dgv_loss_database.overlap(
            chromosome, annotation['outer_start'], annotation['outer_end']
        ))

        annotation['gnomad_del_records'] = list(self._gnomad_del_database.overlap(
            chromosome, annotation['outer_start'], annotation['outer_end']
        ))

        annotation['gnomad_dup_records'] = list(self._gnomad_dup_database.overlap(
            chromosome, annotation['outer_start'], annotation['outer_end']
        ))

        if func == 'del':
            annotation = self._annotate_loss(**annotation)
        elif func == 'dup':
            annotation = self._annotate_gain(**annotation)
        else:
            raise ValueError('Unknown func `{}`'.format(func))

        annotation['rules'], annotation['score'], annotation['pathogenicity'] = self.judge(
            func, **annotation['rules']
        )

        return annotation

    def _serializer(self, anno_result):
        seri = {}
        seri['gene'] = ','.join(x[0].symbol for x in anno_result['overlap_genes'])
        seri['omim_gene'] = ','.join(x[0].symbol for x in anno_result['overlap_omim_genes'])
        seri['HI_gene'] = ','.join(f'{x[0].symbol}({x[1]:.2%};{x[2]:.2%})' for x in anno_result['overlap_hi_genes'])
        seri['HI_region'] = SEP.join(f'{x[0].name}({x[1]:.2%};{x[2]:.2%})' for x in anno_result['overlap_hi_regions'])
        seri['TS_gene'] = ','.join(f'{x[0].symbol}({x[1]:.2%};{x[2]:.2%})' for x in anno_result['overlap_ts_genes'])
        seri['TS_region'] = ','.join(f'{x[0].name}({x[1]:.2%};{x[2]:.2%})' for x in anno_result['overlap_ts_regions'])
        seri['Pred_HI_gene'] = ','.join(f'{x[0].symbol}({x[1]:.2%};{x[2]:.2%})' for x in anno_result['overlap_decipher_genes'])
        seri['auto_evidence'] = ','.join(sorted(anno_result['rules']))
        seri['auto_evidence_score'] = ','.join(f'{k}:{anno_result["rules"][k]}' for k in sorted(anno_result['rules']))
        seri['auto_score'] = anno_result['score']
        seri['auto_pathogenicity'] = anno_result['pathogenicity']
        return seri

    def _seri_anno(self, seri: pd.Series) -> pd.Series:
        anno_result = self.annotate(seri['chr'], seri['start'], seri['end'], seri['type'], seri['error'])
        return seri.append(pd.Series(self._serializer(anno_result)).replace('', '-').fillna(DEFAULT_EMPTY_VALUE))

    def annotation_file(self, file_path, result_path):
        """
        annotate specified file, required columns: chr, start, end, type, error
        :param file_path: input file (TSV)
        :param result_path: result file path (TSV)
        :return: -
        """

        if file_path.endswith('xlsx'):
            input_df = pd.read_excel(file_path)
        else:
            input_df = pd.read_csv(file_path, sep='\t')
        input_df['chr'] = input_df['chr'].map(self._norm_chrom)
        try:
            from tqdm import tqdm
            tqdm.pandas()
            input_df = input_df.progress_apply(self._seri_anno, axis=1)
        except ImportError:
            input_df = input_df.apply(self._seri_anno, axis=1)
        input_df.to_csv(result_path, sep='\t', index=False)
