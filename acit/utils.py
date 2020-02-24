# -*- coding:utf-8 -*-
from pysam import VariantRecord
from json import JSONEncoder


class ACITEncoder(JSONEncoder):
    def default(self, o):
        if isinstance(o, VariantRecord):
            return [o.contig, o.pos, o.ref, o.alts]
        return JSONEncoder.default(self, o)
