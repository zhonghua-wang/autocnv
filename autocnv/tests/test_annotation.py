from autocnv.annotate import AnnotateHelper

annotate = AnnotateHelper()


def test_annotation1():
    # 包含TBX5基因全部（基因 zoom out 3x）
    annotation = annotate.annotate('chr12', 114737220, 114900761, 'del')
    assert '2A' in annotation['rules']


def test_annotation2():
    # TBX5基因5’端，覆盖CDS序列
    annotation = annotate.annotate('chr12', 114800000, 114900761, 'del')
    assert '2C-1' in annotation['rules']


def test_annotation3():
    # TBX5基因5’端，不覆盖CDS序列
    annotation = annotate.annotate('chr12', 114842000, 114900761, 'del')
    assert '2C-2' in annotation['rules']


def test_annotation4():
    # TBX5基因3'端，仅覆盖UTR区
    annotation = annotate.annotate('chr12', 114737220, 114793000, 'del')
    assert '2D-1' in annotation['rules']


def test_annotation5():
    # TBX5基因3'端，仅覆盖末位外显子（该外显子含低频致病变异）
    annotation = annotate.annotate('chr12', 114737220, 114800000, 'del')
    assert '2D-2' in annotation['rules']


def test_annotation6():
    # TBX5基因3'端，覆盖多个exon
    annotation = annotate.annotate('chr12', 114737220, 114835000, 'del')
    assert '2D-4' in annotation['rules']


def test_annotation7():
    # 包含TBX5中间3个exon
    annotation = annotate.annotate('chr12', 114800000, 114835000, 'del')
    assert '2E' in annotation['rules']


def test_annotation8():
    # TBX5基因5’端，不覆盖CDS序列；TBX3基因3'端，仅覆盖末位外显子（该外显子含低频致病变异）
    annotation = annotate.annotate('chr12', 114842000, 115111000, 'del')
    assert '2C-2' in annotation['rules']
    assert '2D-2' in annotation['rules']
    assert '4O' not in annotation['rules']


def test_annotation9():
    # TBX5基因5’端，不覆盖CDS序列；TBX3基因3'端，覆盖多个外显子
    annotation = annotate.annotate('chr12', 114842000, 115114000, 'del')
    assert '2C-2' in annotation['rules']
    assert '2D-4' in annotation['rules']


def test_annotation10():
    # TBX5基因5’端，覆盖CDS序列；TBX3基因3'端，仅覆盖末位外显子（该外显子含低频致病变异）
    annotation = annotate.annotate('chr12', 114800000, 115111000, 'del')
    assert '2C-1' in annotation['rules']
    assert '2D-2' in annotation['rules']


def test_annotation11():
    # TBX5基因5’端，覆盖CDS序列；TBX3基因3'端，覆盖多个外显子
    annotation = annotate.annotate('chr12', 114800000, 115114000, 'del')
    assert '2C-1' in annotation['rules']
    assert '2D-4' in annotation['rules']


def test_annotation12():
    # 包含MECP2基因全部（基因范围 zoom out 1.5x）
    annotation = annotate.annotate('chrX', 153268282, 153382170, 'dup')
    assert '2A' in annotation['rules']
    assert '2H' in annotation['rules']


def test_annotation13():
    # MECP2基因内CNV
    annotation = annotate.annotate('chrX', 153300000, 153363100, 'dup')
    assert '2I' in annotation['rules']


def test_4O():
    annotation = annotate.annotate('chr1', 196757278, 196796716, 'del')
    assert '4O' in annotation['rules']
    annotation = annotate.annotate('chr15', 22750305, 23226254, 'del')
    assert '4O' not in annotation['rules']
    annotation = annotate.annotate('chr1', 25584597, 25767647, 'del')
    assert '4O' not in annotation['rules']
    annotation = annotate.annotate('chr1', 148974342, 149441884, 'dup')
    assert '4O' in annotation['rules']
    annotation = annotate.annotate('chr15', 22750305, 23226254, 'dup')
    assert '4O' not in annotation['rules']
    annotation = annotate.annotate('chr1', 425, 69091, 'dup')
    assert '4O' not in annotation['rules']

def test_syndrome():
    annotation = annotate.annotate('chr23', 6420555, 8153336, 'del')


def test_random():
    annotation = annotate.annotate('chr11', 45904399, 46480747, 'del')
    assert '2H' in annotation['rules']
