# AutoCNV: an **Auto**matic **C**opy **N**umber **V**ariant Interpretation Tool

## Installation

### Pre-requirements

- `pysam`: install `pysam` via conda:

  ```shell
  conda config --add channels r
  conda config --add channels bioconda
  conda install pysam
  ```
  For detail install instruction of `pysam`, please refer to [the official document](https://pysam.readthedocs.io/en/latest/installation.html).

- install [`autopvs1`](https://github.com/JiguangPeng/autopvs1)
- install requirments `pip install -r requirements.txt`

### Usage
please refer to the test cases.

### Citation

```
@article{Fan2021,
    author = {Fan, Chunna and Wang, Zhonghua and Sun, Yan and Sun, Jun and Liu, Xi and Kang, Licheng and Xu, Yingshuo and Yang, Manqiu and Dai, Wentao and Song, Lijie and Wei, Xiaoming and Xiang, Jiale and Huang, Hui and Zhou, Meizhen and Zeng, Fanwei and Huang, Lin and Xu, Zhengfeng and Peng, Zhiyu},
    doi = {10.1186/s12864-021-08011-4},
    isbn = {1286402108011},
    issn = {14712164},
    journal = {BMC Genomics},
    keywords = {AutoCNV,CNV classification,CNV interpretation,Scoring},
    number = {1},
    pages = {1--12},
    pmid = {34615484},
    publisher = {BMC Genomics},
    title = {{AutoCNV: a semiautomatic CNV interpretation system based on the 2019 ACMG/ClinGen Technical Standards for CNVs}},
    volume = {22},
    year = {2021}
}
```