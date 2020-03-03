from acit.annotate import AnnotateHelper

annotate = AnnotateHelper()


def test_annotation1():
    # 包含TBX5基因全部（基因 zoom out 3x）
    annotation = annotate.annotate('chr12', 114737220, 114900761, 'del')
    assert annotation['rules']['2A']


def test_annotation2():
    # TBX5基因5’端，覆盖CDS序列
    annotation = annotate.annotate('chr12', 114800000, 114900761, 'del')
    assert annotation['rules']['2C-1']


def test_annotation3():
    # TBX5基因5’端，不覆盖CDS序列
    annotation = annotate.annotate('chr12', 114842000, 114900761, 'del')
    assert annotation['rules']['2C-2']


def test_annotation4():
    # TBX5基因3'端，仅覆盖UTR区
    annotation = annotate.annotate('chr12', 114737220, 114793000, 'del')
    assert annotation['rules']['2D-1']


def test_annotation5():
    # TBX5基因3'端，仅覆盖末位外显子（该外显子含低频致病变异）
    annotation = annotate.annotate('chr12', 114737220, 114800000, 'del')
    assert annotation['rules']['2D-2']


def test_annotation6():
    # TBX5基因3'端，覆盖多个exon
    annotation = annotate.annotate('chr12', 114737220, 114835000, 'del')
    assert annotation['rules']['2D-4']


def test_annotation7():
    # 包含TBX5中间3个exon
    annotation = annotate.annotate('chr12', 114800000, 114835000, 'del')
    # todo: autoPVS1
    assert annotation['rules']['2E']


def test_annotation8():
    # TBX5基因5’端，不覆盖CDS序列；TBX3基因3'端，仅覆盖末位外显子（该外显子含低频致病变异）
    annotation = annotate.annotate('chr12', 114842000, 115111000, 'del')
    assert annotation['rules']['2C-2']
    assert annotation['rules']['2D-2']


def test_annotation9():
    # TBX5基因5’端，不覆盖CDS序列；TBX3基因3'端，覆盖多个外显子
    annotation = annotate.annotate('chr12', 114842000, 115114000, 'del')
    assert annotation['rules']['2C-2']
    assert annotation['rules']['2D-4']


def test_annotation10():
    # TBX5基因5’端，覆盖CDS序列；TBX3基因3'端，仅覆盖末位外显子（该外显子含低频致病变异）
    annotation = annotate.annotate('chr12', 114800000, 115111000, 'del')
    assert annotation['rules']['2C-1']
    assert annotation['rules']['2D-2']


def test_annotation11():
    # TBX5基因5’端，覆盖CDS序列；TBX3基因3'端，覆盖多个外显子
    annotation = annotate.annotate('chr12', 114800000, 115114000, 'del')
    assert annotation['rules']['2C-1']
    assert annotation['rules']['2D-4']


def test_annotation12():
    # 包含MECP2基因全部（基因范围 zoom out 1.5x）
    annotation = annotate.annotate('chrX', 153268282, 153382170, 'dup')
    assert annotation['rules']['2A']
    assert annotation['rules']['2H']


def test_annotation13():
    # MECP2基因内CNV
    annotation = annotate.annotate('chrX', 153287263, 153363188, 'dup')
    # todo: autoPVS1
    assert annotation['rules']['2I']
