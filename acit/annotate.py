from collections import defaultdict
from acit.database import DataBase
from pysam import VariantFile
from acit import settings
import re


class AnnotateHelper:
    def __init__(self):
        self._gene_database = DataBase(settings.GENE_DATABASE)
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

    @staticmethod
    def _annotate_loss(**annotation):
        loss = defaultdict(float)

        for gene, overlap, coverage in annotation['overlap_hi_genes']:
            # 覆盖整个基因
            if coverage == 1:
                # loss['2A'] = max(loss['2A'], 1)
                loss['2A'] = True
            # 不位于基因内部
            elif overlap < 1:
                # 覆盖末位外显子
                if any(
                        exon.last_exon == 'True'
                        for exon, *_ in annotation['overlap_hi_exons'][gene.gene_id]
                ):
                    # 覆盖超过两个外显子
                    if len(annotation['overlap_hi_exons'][gene.gene_id]) > 1:
                        # loss['2D-4'] = max((loss['2D-4'], 0.9))
                        loss['2D-4'] = True
                    # 仅覆盖末位外显子
                    elif len(annotation['variants']) > 0:
                        # loss['2D-2'] = max((loss['2D-2'], 0.9))
                        loss['2D-2'] = True
                    else:
                        # loss['2D-3'] = max((loss['2D-3'], 0.3))
                        loss['2D-3'] = True
                # 覆盖5'端CDS区
                elif gene.gene_id in annotation['overlap_hi_cds'] and \
                        len(annotation['overlap_hi_cds'][gene.gene_id]) > 0:
                    # loss['2C-1'] = max((loss['2C-1'], 0.9))
                    loss['2C-1'] = True

        # 完全覆盖hi区域
        for region, overlap, coverage in annotation['overlap_hi_regions']:
            if coverage == 1:
                # loss['2A'] = max((loss['2A'], 1))
                loss['2A'] = True

        # 落入uhi基因
        for gene, overlap, coverage in annotation['overlap_uhi_genes']:
            if overlap == 1:
                # loss['2F'] = -1
                loss['2F'] = True

        # 落入uhi区域
        for region, overlap, coverage in annotation['overlap_uhi_regions']:
            genes = set(gene.symbol for gene, *_ in annotation['overlap_genes'])
            if overlap == 1:  # 完全落入区域
                # loss['2F'] = -1
                loss['2F'] = True
            elif set(region.genes.split(',')) == genes:  # 覆盖相同的基因
                # loss['2F'] = -1
                loss['2F'] = True

        # 包含decipher基因
        if '2F' not in loss and len(annotation['overlap_decipher_genes']) > 0:
            # loss['2H'] = 0.15
            loss['2H'] = True

        # 覆盖基因个数
        gene_count = len(annotation['overlap_genes'])
        if gene_count >= 25:
            if gene_count < 35:
                # loss['3B'] = 0.45
                loss['3B'] = True
            else:
                # loss['3C'] = 0.9
                loss['3C'] = True

        annotation['rules'] = loss
        return annotation

    @staticmethod
    def _annotate_gain(**annotation):
        gain = defaultdict(float)

        for gene, overlap, coverage in annotation['overlap_ts_genes']:
            # 覆盖整个基因
            if coverage == 1:
                # gain['2A'] = max(gain['2A'], 1)
                gain['2A'] = True

        # 完全覆盖ts区域
        for region, overlap, coverage in annotation['overlap_ts_regions']:
            if coverage == 1:
                # gain['2A'] = max((gain['2A'], 1))
                gain['2A'] = True

        # 落入uts基因
        for gene, overlap, coverage in annotation['overlap_uts_genes']:
            if overlap == 1:
                # gain['2D'] = -1
                gain['2D'] = True

        # 落入uts区域
        for region, overlap, coverage in annotation['overlap_uts_regions']:
            genes = set(gene.symbol for gene, *_ in annotation['overlap_genes'])
            if overlap == 1:  # 完全落入区域
                # gain['2D'] = -1
                gain['2D'] = True
            elif set(region.genes.split(',')) == genes:  # 覆盖相同的基因
                # gain['2F'] = -1
                gain['2F'] = True

        # todo: 2.2.4

        # 覆盖基因个数
        gene_count = len(annotation['overlap_genes'])
        if gene_count >= 35:
            if gene_count < 50:
                # gain['3B'] = 0.45
                gain['3B'] = True
            else:
                # gain['3C'] = 0.9
                gain['3C'] = True

        annotation['rules'] = gain
        return annotation

    def annotate(self, chromosome, start, end, func, error=0):
        annotation = dict(
            chromosome=chromosome, start=start, end=end,
            length=end - start, error=error,
            outer_start=start - error, outer_end=end + error,
            inner_start=start + error, inner_end=end - error,
            func=func
        )

        annotation['overlap_genes'] = list(self._gene_database.overlap(
            chromosome, annotation['outer_start'], annotation['outer_end']
        ))
        # annotation['genes'] = '\n'.join(
        #     '{}({})/{:.2%}'.format(
        #         gene, tx, fraction
        #     ) for (*_, gene, tx), _, fraction in annotation['overlap_genes']
        # )

        annotation['overlap_hi_genes'] = list(self._hi_gene_database.overlap(
            chromosome, annotation['inner_start'], annotation['inner_end']
        ))
        # annotation['hi_genes'] = '\n'.join(
        #     '{}({})/{:.2%}'.format(
        #         gene, tx, fraction
        #     ) for (*_, gene, tx), _, fraction in annotation['overlap_hi_genes']
        # )

        annotation['overlap_hi_exons'] = self._hi_exon_database.overlap_groups(
            chromosome, annotation['inner_start'], annotation['inner_end'],
            lambda record: record[0].gene_id
        )
        # annotation['hi_exons'] = '\n'.join(
        #     '{}({}):{}'.format(
        #         gene, tx, re.sub(r'-.+-', '-', '-'.join(map(
        #             lambda x: 'EX{}'.format(x), sorted(int(record[0][6]) for record in records)
        #         )))
        #     ) for (gene, tx), records in annotation['overlap_hi_exons'].items()
        # )

        annotation['overlap_hi_cds'] = self._hi_cds_database.overlap_groups(
            chromosome, annotation['inner_start'], annotation['inner_end'],
            lambda record: record[0].gene_id
        )

        annotation['variants'] = list(self._clinvar_pathogenic_database.fetch(
            chromosome, annotation['inner_start'], annotation['inner_end']
        ))

        annotation['overlap_hi_regions'] = list(self._hi_region_database.overlap(
            chromosome, annotation['inner_start'], annotation['inner_end']
        ))

        annotation['overlap_uhi_genes'] = list(self._uhi_gene_database.overlap(
            chromosome, annotation['outer_start'], annotation['outer_end']
        ))

        annotation['overlap_uhi_regions'] = list(self._uhi_region_database.overlap(
            chromosome, annotation['inner_start'], annotation['inner_end']
        ))

        annotation['overlap_decipher_genes'] = list(self._decipher_gene_database.overlap(
            chromosome, annotation['inner_start'], annotation['inner_end']
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
            chromosome, annotation['inner_start'], annotation['inner_end']
        ))

        if func == 'del':
            annotation = self._annotate_loss(**annotation)
        elif func == 'dup':
            annotation = self._annotate_gain(**annotation)
        else:
            raise ValueError('Unknown func `{}`'.format(func))

        return annotation
