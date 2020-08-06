from collections import defaultdict
from collections import namedtuple
import pysam


class DataBase:
    def __init__(self, path):
        self._path = path
        self._tbx = pysam.Tabixfile(path)
        self._fields = namedtuple('Record', self._tbx.header[-1].strip('#').split('\t'))

    def fetch(self, chrom, start, end):
        """
        查找并生成给定基因组位置的记录
        :param chrom: 染色体编号
        :param start: 起始位置
        :param end: 终止位置
        :return: 记录生成器
        """
        try:
            for record in self._tbx.fetch(chrom, start, end):
                chrom, start, end, *fields = record.split('\t')
                start, end = int(start), int(end)
                # 所得记录按照表头组装为namedtuple方便使用
                yield self._fields(chrom, start, end, *fields)
        except ValueError:
            yield from ()

    def overlap(self, chrom, start, end):
        """
        查找并生成给定基因组位置的记录，同时计算两者之间的重叠程度
        :param chrom: 染色体编号
        :param start: 起始位置
        :param end: 终止位置
        :return: 记录与重叠程度生成器
        """
        length = end - start
        for record in self.fetch(chrom, start, end):
            _, overlap_start, overlap_end, _ = sorted((start, end, record[1], record[2]))
            overlap = overlap_end - overlap_start
            yield record, overlap / length, overlap / (record[2] - record[1])

    def overlap_groups(self, chrom, start, end, key=None):
        """
        查找并生成给定基因组位置的记录，同时计算两者之间的重叠程度，返回按照key方法分组的字典
        :param chrom: 染色体编号
        :param start: 起始位置
        :param end: 终止位置
        :param key: 分组依据
        :return: 分组后的记录字典
        """
        groups = defaultdict(list)
        for record in self.overlap(chrom, start, end):
            groups[record if key is None else key(record)].append(record)
        return dict(groups)
