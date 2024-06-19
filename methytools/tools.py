import argparse
import gzip
from importlib import resources
import pandas as pd
import numpy as np
import time
from pathlib import Path


def load_center(element='tss'):
    if element == 'tss':
        center_file = resources.path('methytools.ref.hg38', 'center.tss.bed.gz')
    elif element == 'cgi':
        center_file = resources.path('methytools.ref.hg38', 'center.cgi.bed.gz')
    else:
        raise Exception('Unknown genomic element {}'.format(element))
    out = []
    with gzip.open(center_file.__str__(), 'rt') as f:
        for line in f:
            row = line.strip()
            item = row.split('\t')
            out.append([item[0], int(item[1]), int(item[2])])
    return pd.DataFrame(out)


def extract_signal(ratio_txt, element, out_matrix, up=10000, down=10000):
    center_file = load_center(element)
    ratio = pd.read_csv(ratio_txt, sep='\t')
    ratio_col = ratio.columns.to_list()
    ratio_col[0] = 'chrom'
    ratio.columns = ratio_col
    sample = ratio.columns[3:].to_list()
    SAMPLE = len(sample)

    top_count = np.zeros((up + down + 1, SAMPLE))
    top_ratio = np.zeros((up + down + 1, SAMPLE))

    ratio_map = dict()
    for chrom, data in ratio.groupby("chrom"):
        ratio_map[chrom] = data

    total = len(center_file)
    start = time.time()
    idx = 0
    for i in range(len(center_file)):
        if idx % 10000 == 0:
            print("index {}/{}, pass {} min".format(idx, total, round((time.time() - start) / 60), 4))
        idx += 1
        chrom = center_file.iloc[i, 0]
        tss = center_file.iloc[i, 1] - 1
        sign = center_file.iloc[i, 2]
        ratio2 = ratio_map[chrom]
        interval = [tss - up, tss + down]
        if sign == '-':
            interval = [tss - down, tss + up]
        data = ratio2.loc[(ratio2['start'] >= interval[0]) & (ratio2['start'] <= interval[1])]
        left = pd.DataFrame({
            'chrom': chrom,
            'start': list(range(interval[0], interval[1] + 1))
        })
        a = pd.merge(left, data, how='left', left_on=['chrom', 'start'], right_on=['chrom', 'start'])
        if sign == '-':
            a = a.iloc[::-1]
        a[a == -1] = np.nan
        mm = a.iloc[:, 3:].to_numpy()
        mm_count = mm.copy()
        mm_count[~np.isnan(mm)] = 1
        mm_count[np.isnan(mm)] = 0
        top_count += mm_count
        mm[np.isnan(mm)] = 0
        top_ratio += mm
    np.seterr(invalid='ignore')
    ratio_matrix = top_ratio / top_count
    matrix = pd.DataFrame(ratio_matrix, columns=ratio.columns[3:])
    matrix['tss'] = list(range(-up, down + 1))
    matrix.set_index('tss', inplace=True)
    with open(out_matrix.__str__(), "w") as f:
        matrix.to_csv(f, sep="\t")
    end = time.time()
    print()
    print("use {} min !!".format((end - start) / 60))
    return matrix


def analysis_genomic_element(ratio_txt, element, out, up=10000, down=10000):
    matrix = extract_signal(ratio_txt, element, Path(out), up, down)


if __name__ == '__main__':
    pass
