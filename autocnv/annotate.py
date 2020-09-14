from autocnv.database import DataBase
from pysam import VariantFile
from autocnv import settings
from autopvs1.cnv import CNVRecord, PVS1CNV
from autopvs1.utils import get_transcript
from autopvs1.read_data import transcripts
from autopvs1.strength import Strength
from collections import defaultdict
from itertools import chain
import operator
import pandas as pd

SEP = '\n'
DEFAULT_EMPTY_VALUE = '-'

PVS1 = {
    Strength.VeryStrong: 'PVS1', Strength.Strong: 'PVS1_S', Strength.Moderate: 'PVS1_M',
    Strength.Supporting: 'PVS1_P', Strength.Unmet: 'PVS1_U'
}

# 计分分组配置，同一组证据仅计算最大分值
SCORE_GROUP = {
    'del': {
        rule: 'G1' for rule in ('2A', '2B', '2C-1', '2C-2', '2D-1', '2D-2', '2D-3', '2D-4', '2E')
    },
    'dup': {}
}

# 致病性判断分级配置
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
        self._clinvar_pathogenic_database = VariantFile(
            settings.CLINVAR_PATHOGENIC_DATABASE)
        self._uhi_gene_database = DataBase(settings.UHI_GENE_DATABASE)
        self._hi_region_database = DataBase(settings.HI_REGION_DATABASE)
        self._uhi_region_database = DataBase(settings.UHI_REGION_DATABASE)
        self._decipher_gene_database = DataBase(
            settings.DECIPHER_GENE_DATABASE)
        self._ts_gene_database = DataBase(settings.TS_GENE_DATABASE)
        self._ts_region_database = DataBase(settings.TS_REGION_DATABASE)
        self._uts_gene_database = DataBase(settings.UTS_GENE_DATABASE)
        self._uts_region_database = DataBase(settings.UTS_REGION_DATABASE)
        self._dgv_gain_database = DataBase(settings.DGV_GAIN_DATABASE)
        self._dgv_loss_database = DataBase(settings.DGV_LOSS_DATABASE)
        self._gnomad_del_database = DataBase(settings.GNOMAD_DEL_DATABASE)
        self._gnomad_dup_database = DataBase(settings.GNOMAD_DUP_DATABASE)
        self._cnv_syndrome_del_database = DataBase(settings.CNV_SYNDROME_DEL_DATABASE)
        self._cnv_syndrome_dup_database = DataBase(settings.CNV_SYNDROME_DUP_DATABASE)

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
        """
        计算拷贝数减少的CNV的证据项
        :param annotation: 已注释的CNV
        :return: 注释后的CNV
        """
        loss = dict()

        # Section 1

        if len(annotation['outer_overlap_genes']) + len(annotation['overlap_func_regions']) > 0:
            loss['1A'] = True
        else:
            loss['1B'] = True

        # Section 2

        # hi区域
        for region, overlap, coverage in annotation['overlap_hi_regions']:
            if coverage == 1:  # 完全覆盖区域
                loss['2A'] = True
            elif len(set(gene.symbol for gene, *_ in annotation['overlap_hi_genes'])) == 0:
                # 未覆盖hi基因
                loss['2B'] = True

        # hi基因
        for gene, overlap, coverage in annotation['overlap_hi_genes']:
            if coverage == 1:  # 完全覆盖基因
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
                        else:  # 末尾外显子无致病变异
                            loss['2D-3'] = True
                    else:
                        # 不覆盖CDS区
                        loss['2D-1'] = True
                # 未覆盖末位外显子
                elif gene.gene_id in annotation['overlap_hi_cds'] \
                        and len(annotation['overlap_hi_cds'][gene.gene_id]) > 0:  # 是否覆盖5'端CDS
                    loss['2C-1'] = True
                else:  # 未覆盖5'端CDS
                    loss['2C-2'] = True
            # 位于基因内部
            else:
                cnv = CNVRecord(
                    annotation['chromosome'], annotation['inner_start'],
                    annotation['inner_end'], annotation['func']
                )
                tx = get_transcript(gene.transcript, transcripts)
                pvs1 = PVS1CNV(cnv, None, tx)
                loss['2E'] = True
                #loss[PVS1[pvs1.verify_DEL()[0]]] = True
                loss['pvs1'] = PVS1[pvs1.verify_DEL()[0]]

        # 包含预测HI基因
        if len(annotation['overlap_hi_genes']) + len(annotation['overlap_hi_regions']) == 0 \
                and len(annotation['overlap_decipher_genes']) > 0:
            loss['2H'] = True

        # 落入uhi基因
        for gene, overlap, coverage in annotation['overlap_uhi_genes']:
            if overlap == 1:
                loss['2F'] = True

        # 落入uhi区域
        genes = set(gene.symbol for gene, *
                    _ in annotation['outer_overlap_genes'])
        for region, overlap, coverage in annotation['overlap_uhi_regions']:
            if len(genes - set(region.genes.split(','))) > 0:
                loss['2G'] = True
            else:
                loss['2F'] = True

        # Section 3

        # 覆盖基因个数
        gene_count = len(annotation['outer_overlap_genes'])
        if gene_count >= 35:
            loss['3C'] = True
        elif gene_count >= 25:
            loss['3B'] = True
        elif gene_count >= 0:
            loss['3A'] = True

        # Section 4

        # DGV金标和Gnomad
        genes = set(gene.symbol for gene, *
                    _ in annotation['outer_overlap_genes'])
        l, m = 0, 0
        for record, overlap, coverage in chain(
                annotation['dgv_loss_records'], annotation['gnomad_del_records']
        ):
            if overlap == 1 and any(
                float(v) >= 0.01 for f, v in record._asdict().items() if f.startswith('af')
            ):  # 完全覆盖待解读CNV且频率大于1%
                loss['4O'] = True
                break
            elif overlap >= 0.5 and len(genes - set(record.genes.split(','))) == 0:
                # 与待解读CNV重叠超过50%且覆盖全部蛋白编码基因
                if any(float(v) < 0.01 for f, v in record._asdict().items() if f.startswith('af')):
                    # 频率小于1%
                    m += 1
                else:
                    # 频率大于1%
                    l += 1
        else:
            if l > 0 and m == 0:  # 存在频率大于1%且不存在小于1%的CNV
                loss['4O'] = True

        annotation['rules'] = loss
        return annotation

    @staticmethod
    def _annotate_gain(**annotation):
        """
        计算拷贝数减少的CNV的证据项
        :param annotation: 已注释的CNV
        :return: 注释后的CNV
        """
        gain = dict()

        # Section 1

        if len(annotation['outer_overlap_genes']) + len(annotation['overlap_func_regions']) > 0:
            gain['1A'] = True
        else:
            gain['1B'] = True

        # Section 2

        # 完全覆盖ts区域
        for region, overlap, coverage in annotation['overlap_ts_regions']:
            if coverage == 1:  # 是否覆盖整改区域
                gain['2A'] = True
            elif len(set(gene.symbol for gene, *_ in annotation['overlap_ts_genes'])) == 0:
                # 未覆盖ts基因
                gain['2B'] = True

        for gene, overlap, coverage in annotation['overlap_ts_genes']:
            # 覆盖整个基因
            if coverage == 1:
                gain['2A'] = True

        # 落入uts基因
        for gene, overlap, coverage in annotation['overlap_uts_genes']:
            if overlap == 1:
                gain['2D'] = True

        # 落入uts区域
        for region, overlap, coverage in annotation['overlap_uts_regions']:
            genes = set(gene.symbol for gene, *
                        _ in annotation['inner_overlap_genes'])
            region_genes = set(region.genes.split(','))
            if overlap == coverage == 1:  # 与良性区域完全一致
                gain['2C'] = True
            elif len(genes - region_genes) > 0:  # 编码蛋白基因比良性区域多
                gain['2G'] = True
            # 破坏蛋白编码基因
            elif any(c < 1 for *_, c in annotation['inner_overlap_genes']):
                gain['2E'] = True
            elif overlap == 1:  # 被良性区域完全覆盖
                gain['2D'] = True
            else:
                gain['2F'] = True

        # hi基因
        hi_genes = set()
        for gene, overlap, coverage in annotation['overlap_hi_genes']:
            hi_genes.add(gene.symbol)
            if coverage == 1:  # 完全覆盖
                gain['2H'] = True
            elif overlap == 1:  # 两端均位于基因内
                cnv = CNVRecord(
                    annotation['chromosome'], annotation['inner_start'],
                    annotation['inner_end'], annotation['func']
                )
                tx = get_transcript(gene.transcript, transcripts)
                pvs1 = PVS1CNV(cnv, None, tx)
                gain['2I'] = True
                # gain[PVS1[pvs1.verify_DUP()[0]]] = True
                gain['pvs1'] = PVS1[pvs1.verify_DUP()[0]]

        # 非hi基因
        for gene, overlap, coverage in annotation['inner_overlap_genes']:
            if gene.symbol not in hi_genes and coverage != 1:
                gain['2L'] = True
                annotation['break_point_genes'].append(gene.symbol)

        # Section 3

        # 覆盖基因个数
        gene_count = len(annotation['inner_overlap_genes'])
        if gene_count >= 50:
            gain['3C'] = True
        elif gene_count >= 35:
            gain['3B'] = True
        elif gene_count >= 0:
            gain['3A'] = True

        # Section 4

        # DGV金标和Gnomad
        genes = set(gene.symbol for gene, *
                    _ in annotation['outer_overlap_genes'])
        l, m = 0, 0
        for record, overlap, coverage in chain(
                annotation['dgv_gain_records'], annotation['gnomad_dup_records']
        ):
            if overlap == 1 and any(
                    float(v) >= 0.01 for f, v in record._asdict().items() if f.startswith('af')
            ):  # 完全覆盖待解读CNV且频率大于1%
                gain['4O'] = True
                break
            elif overlap >= 0.5 and len(genes - set(record.genes.split(','))) == 0:
                # 与待解读CNV重叠超过50%且覆盖全部蛋白编码基因
                if any(float(v) < 0.01 for f, v in record._asdict().items() if
                       f.startswith('af')):
                    # 频率小于1%
                    m += 1
                else:
                    # 频率大于1%
                    l += 1
        else:
            if l > 0 and m == 0:  # 存在频率大于1%且不存在小于1%的CNV
                gain['4O'] = True

        annotation['rules'] = gain
        return annotation

    @staticmethod
    def merge_score(func, **rules):
        """
        整合所有证据项得分
        :param func: 变异类型
        :param rules: 证据项
        :return: 生成各证据项得分
        """
        groups = defaultdict(list)
        for rule, score in rules.items():
            try:  # 需要分组计分的证据项先收集起来
                groups[SCORE_GROUP[func][rule]].append(score)
            except KeyError:  # 无需分组计分的证据项直接计分
                yield score
        for _, scores in groups.items():  # 分组计分的证据项只计算最大分值
            yield max(scores)

    @staticmethod
    def judge(func, **rules):
        """
        判断给定的证据项组合最终的致病性
        :param func: 变异类型
        :param rules: 勾选的证据项
        :return: 证据项、得分和致病性
        """
        # 获取所有证据项得分
        # rules = {
        #     rule: settings.DEFAULT_SCORE[func][rule] for rule, check in rules.items() if check
        # }
        rules_value = {}
        for rule, check in rules.items():
            if check in PVS1.values():
                rules_value['pvs1'] = settings.DEFAULT_SCORE[func][check]
            elif check:
                rules_value[rule] = settings.DEFAULT_SCORE[func][rule]
        # 整合所有证据项得分
        score = sum(AnnotateHelper.merge_score(func, **rules_value))
        # 判断致病性
        for op, cutoff, level in PATHOGENICITY_LEVELS[:-1]:
            if op(score, cutoff):
                pathogenicity = level
                break
        else:
            pathogenicity = PATHOGENICITY_LEVELS[-1][2]
        return rules_value, score, pathogenicity

    def annotate(self, chromosome, start, end, func, error=0):
        """
        对给定CNV进行注释
        :param chromosome: 染色体编号
        :param start: 起始位置
        :param end: 终止位置
        :param func: 变异类型
        :param error: 误差值
        :return: 注释结果
        """
        annotation = dict(
            chromosome=chromosome, start=start, end=end,
            length=end - start, error=error,
            outer_start=start - error, outer_end=end + error,
            inner_start=start + error, inner_end=end - error,
            func=func, break_point_genes=list()
        )

        annotation['inner_overlap_genes'] = list(self._gene_database.overlap(
            chromosome, annotation['inner_start'], annotation['inner_end'],
        ))

        annotation['outer_overlap_genes'] = list(self._gene_database.overlap(
            chromosome, annotation['outer_start'], annotation['outer_end'],
        ))

        annotation['overlap_omim_genes'] = list(self._omim_gene_database.overlap(
            chromosome, annotation['inner_start'], annotation['inner_end']
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

        annotation['cnv_syndrome_loss'] = list(self._cnv_syndrome_del_database.overlap(
            chromosome, annotation['outer_start'], annotation['outer_end']
        ))
        annotation['cnv_syndrome_gain'] = list(self._cnv_syndrome_dup_database.overlap(
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
        # PVS1
        if func == 'del' and '2E' in annotation['rules'].keys():
            annotation['rules']['2E'] = annotation['rules'].get('pvs1')
        elif func == 'dup' and '2I' in annotation['rules'].keys():
            annotation['rules']['2I'] = annotation['rules'].get('pvs1')
        annotation['pvs1'] = annotation['rules'].pop('pvs1', None)

        return annotation

    def _serializer(self, anno_result):
        seri = {}
        seri['inner_gene'] = ','.join(
            x[0].symbol for x in anno_result['inner_overlap_genes'])
        seri['inner_omim_gene'] = ','.join(
            x[0].symbol for x in anno_result['overlap_omim_genes'])
        seri['HI_gene'] = ','.join(
            f'{x[0].symbol}({x[1]:.2%};{x[2]:.2%})' for x in anno_result['overlap_hi_genes'])
        seri['HI_region'] = SEP.join(
            f'{x[0].name}({x[1]:.2%};{x[2]:.2%})' for x in anno_result['overlap_hi_regions'])
        seri['TS_gene'] = ','.join(
            f'{x[0].symbol}({x[1]:.2%};{x[2]:.2%})' for x in anno_result['overlap_ts_genes'])
        seri['TS_region'] = ','.join(
            f'{x[0].name}({x[1]:.2%};{x[2]:.2%})' for x in anno_result['overlap_ts_regions'])
        seri['Pred_HI_gene'] = ','.join(
            f'{x[0].symbol}({x[1]:.2%};{x[2]:.2%})' for x in anno_result['overlap_decipher_genes'])
        seri['auto_evidence'] = ','.join(sorted(anno_result['rules']))
        seri['auto_evidence_score'] = ','.join(
            f'{k}:{anno_result["rules"][k]}' for k in sorted(anno_result['rules']))
        seri['benign_hi_gene'] = ','.join(
            f'{x[0].symbol}({x[1]:.2%};{x[2]:.2%})' for x in anno_result['overlap_uhi_genes'])
        seri['benign_hi_region'] = ','.join(
            f'{x[0].name}({x[1]:.2%};{x[2]:.2%})' for x in anno_result['overlap_uhi_regions'])
        seri['benign_ts_gene'] = ','.join(
            f'{x[0].symbol}({x[1]:.2%};{x[2]:.2%})' for x in anno_result['overlap_uts_genes'])
        seri['benign_ts_region'] = ','.join(
            f'{x[0].name}({x[1]:.2%};{x[2]:.2%})' for x in anno_result['overlap_uts_regions'])
        seri['dgv_loss_records'] = ','.join(
            f'{x[0].id}(af: {x[0].af})({x[1]:.2%};{x[2]:.2%})' for x in anno_result['dgv_loss_records']
        )
        seri['dgv_gain_records'] = ','.join(
            f'{x[0].id}(af: {x[0].af})({x[1]:.2%};{x[2]:.2%})' for x in anno_result['dgv_gain_records']
        )
        seri['gnomad_loss_records'] = ','.join(
            f'{x[0].name}(af: {x[0].af})({x[1]:.2%};{x[2]:.2%})' for x in anno_result['gnomad_del_records']
        )
        seri['gnomad_gain_records'] = ','.join(
            f'{x[0].name}(af: {x[0].af})({x[1]:.2%};{x[2]:.2%})' for x in anno_result['gnomad_dup_records']
        )
        seri['cnv_syndrome_gain'] = ','.join(
            f'{x[0].disease_name}({x[1]:.2%};{x[2]:.2%})' for x in anno_result['cnv_syndrome_gain']
        )
        seri['cnv_syndrome_loss'] = ','.join(
            f'{x[0].disease_name}({x[1]:.2%};{x[2]:.2%})' for x in anno_result['cnv_syndrome_loss']
        )
        seri['auto_score'] = anno_result['score']
        seri['auto_pathogenicity'] = anno_result['pathogenicity']
        seri['pvs1'] = anno_result['pvs1']
        return seri

    def _seri_anno(self, seri: pd.Series) -> pd.Series:
        anno_result = self.annotate(seri['chr'], seri['start'], seri['end'], seri['type'],
                                    seri['error'])
        return seri.append(
            pd.Series(self._serializer(anno_result)).replace('', '-').fillna(DEFAULT_EMPTY_VALUE))

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
        if result_path.endswith('xlsx'):
            input_df.to_excel(result_path, index=False)
        else:
            input_df.to_csv(result_path, sep='\t', index=False)
