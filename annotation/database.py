import pysam


class DataBase:
    def __init__(self, path):
        self._path = path
        self._tbx = pysam.Tabixfile(path)

    def fetch(self, chrom, start, end):
        for record in self._tbx.fetch(chrom, start, end):
            chrom, start, end, *fields = record.split('\t')
            start, end = int(start), int(end)
            yield tuple([chrom, start, end] + fields)
    
    def overlap(self, chrom, start, end):
        length = end - start
        for record in self.fetch(chrom, start, end):
            _, overlap_start, overlap_end, _ = sorted((start, end, record[1], record[2]))
            overlap = overlap_end - overlap_start
            yield record, overlap / length, overlap / (record[2] - record[1])
