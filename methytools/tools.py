import pysam
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


def analysis_genomic_element_from_gz(ratio_txt, element, out_matrix, up=50, down=50):
    ratio_txt_csi = ratio_txt + '.csi'
    if not Path(ratio_txt_csi).exists():
        raise FileNotFoundError(f"{ratio_txt_csi} not exist!")
    center_file = load_center(element)
    tabix_file = pysam.TabixFile(ratio_txt, index=ratio_txt_csi)
    header = tabix_file.header[0]
    cols = header.split('\t')
    sample_count = len(cols) - 3
    with open(out_matrix, 'w') as fo:
        finals = []
        total = len(center_file)
        start = time.time()
        idx = 0
        for v in center_file.values:
            if idx % 10000 == 0:
                print(f"index {idx}/{total}, pass {round(time.time() - start)} s")
            idx += 1
            region = f'{v[0]}:{int(v[1]) + up}-{int(v[2]) + down}'
            _iterator = tabix_file.fetch(region=region)
            lines = []
            for line in _iterator:
                line = line.strip().split('\t')
                lines.append([float(x) for x in line[3:]])
            if len(lines):
                m = np.array(lines)
                m[m == -1] = np.nan
                m = np.nanmean(m, axis=0)
            else:
                m = np.array([np.nan] * sample_count)
            finals.append(m)
        finals = np.array(finals)
        finals = np.nanmean(finals, axis=0)
        fo.write('\t'.join(cols[3:]) + '\n')
        fo.write('\t'.join([str(round(x, 4)) for x in finals]) + '\n')


if __name__ == '__main__':
    pass
