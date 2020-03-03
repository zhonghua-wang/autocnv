from acit.database import DataBase
from pysam import VariantFile
from acit import settings
from autopvs1.cnv import CNVRecord, PVS1CNV
from autopvs1.utils import get_transcript
from autopvs1.read_data import transcripts


class AnnotateHelper:
    def __init__(self):
        self._gene_database = DataBase(settings.GENE_DATABASE)
        self._omim_gene_database = DataBase(settings.OMIM_GENE_DATABASE)
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
        self._dgv_database = DataBase(settings.DGV_DATABASE)

    @staticmethod
    def _annotate_loss(**annotation):
        loss = dict()

        if len(annotation['overlap_genes']) + len(annotation['overlap_hi_regions']) > 0:
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
                annotation['autoPVS1'] = pvs1.verify_DEL()

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

        # 包含预测HI基因
        if len(annotation['overlap_hi_genes']) + len(annotation['overlap_hi_regions']) == 0 \
                and len(annotation['overlap_decipher_genes']) > 0:
            loss['2H'] = True

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
        for record, *_ in annotation['dgv_records']:
            if genes - set(record.genes.split(',')):
                loss['4O'] = True

        annotation['rules'] = loss
        return annotation

    @staticmethod
    def _annotate_gain(**annotation):
        gain = dict()

        if len(annotation['overlap_genes']) + len(annotation['overlap_ts_regions']) > 0:
            gain['1A'] = True
        else:
            gain['1B'] = True

        for gene, overlap, coverage in annotation['overlap_ts_genes']:
            # 覆盖整个基因
            if coverage == 1:
                gain['2A'] = True
            else:
                gain['2B'] = True

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
            if len(genes - region_genes) > 0:  # 多
                gain['2G'] = True
            elif len(region_genes - genes) > 0:  # 少
                if any(coverage != 1 for *_, coverage in annotation['overlap_genes']):
                    gain['2E'] = True
                else:
                    gain['2D'] = True
            else:
                gain['2C'] = True

        for gene, overlap, coverage in annotation['overlap_hi_genes']:
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
                annotation['autoPVS1'] = pvs1.verify_DUP()

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
        for record, *_ in annotation['dgv_records']:
            if genes - set(record.genes.split(',')):
                gain['4O'] = True

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

        annotation['overlap_omim_genes'] = list(self._omim_gene_database.overlap(
            chromosome, annotation['outer_start'], annotation['outer_end']
        ))

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

        try:
            annotation['variants'] = list(self._clinvar_pathogenic_database.fetch(
                chromosome, annotation['inner_start'], annotation['inner_end'])
            )
        except ValueError:
            annotation['variants'] = []

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

        annotation['dgv_records'] = list(self._dgv_database.overlap(
            chromosome, annotation['outer_start'], annotation['outer_end']
        ))

        if func == 'del':
            annotation = self._annotate_loss(**annotation)
        elif func == 'dup':
            annotation = self._annotate_gain(**annotation)
        else:
            raise ValueError('Unknown func `{}`'.format(func))

        return annotation
