#%%
from autocnv.annotate import AnnotateHelper
import pandas as pd

test_file = '/Users/zhonghua/test/test.xlsx'
anno = AnnotateHelper()

#%%
# df = pd.read_excel(test_file)
# df = df[df['chromosome'] != 24]
# df['chromosome'] = df['chromosome'].map(anno._norm_chrom)
#
# df.apply(anno._seri_anno, axis=1)

annotation = anno.annotate('chrX', 6420555, 8153336, 'del')
d = anno.serializer(annotation)
print(d)

