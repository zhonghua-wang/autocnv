# -*- coding:utf-8 -*-
from pysam import VariantRecord
from json import JSONEncoder
from autopvs1.strength import Strength
from autocnv import settings

class ACITEncoder(JSONEncoder):
    def default(self, o):
        if isinstance(o, VariantRecord):
            return [o.contig, o.pos, o.ref, o.alts]
        if isinstance(o, Strength):
            return o.name
        return JSONEncoder.default(self, o)
