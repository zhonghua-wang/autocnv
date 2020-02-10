from collections import defaultdict
from collections import namedtuple
import pysam


class DataBase:
    def __init__(self, path):
        self._path = path
        self._tbx = pysam.Tabixfile(path)
        self._fields = namedtuple('Record', self._tbx.header[-1].strip('#').split('\t'))

    def fetch(self, chrom, start, end):
        try:
            for record in self._tbx.fetch(chrom, start, end):
                chrom, start, end, *fields = record.split('\t')
                start, end = int(start), int(end)
                yield self._fields(chrom, start, end, *fields)
        except ValueError:
            yield from ()

    def overlap(self, chrom, start, end):
        length = end - start
        for record in self.fetch(chrom, start, end):
            _, overlap_start, overlap_end, _ = sorted((start, end, record[1], record[2]))
            overlap = overlap_end - overlap_start
            yield record, overlap / length, overlap / (record[2] - record[1])

    def overlap_groups(self, chrom, start, end, key=None):
        groups = defaultdict(list)
        for record in self.overlap(chrom, start, end):
            groups[record if key is None else key(record)].append(record)
        return dict(groups)
