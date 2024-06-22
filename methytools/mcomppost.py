import sys
import argparse
import numpy as np
import pandas as pd
import time
import pysam
from pathlib import Path
from .utils import get_param_from_doc


def mcomp_dmc_filter(dmc_file, matrix_bed, target_samples, dmc_file_filtered=None, percentile=0.8,
                     minimum_sample_keep_ratio=0.8, verbose=True):
    """
    This tool is employed to further refine and filter the DMCs identified by MOABS:mcomp.

    :param dmc_file: the dmc files generate by mcomp.
    :param matrix_bed: the bed format methylation matrix. This file should bed compressed by bgzip and index in csi format.
    :param target_samples: the names of samples in target group. This names should appear in the header of *matrix_bed*.
    :param dmc_file_filtered: the filtered dmc file. If not set, it will print to stdout.
    :param percentile: The p-th percentile of the low methylation group will be compared with the (1−p)-th percentile
                    of the high methylation group. If the former is less than the latter, the condition is satisfied.
    :param minimum_sample_keep_ratio: The minimum proportion of each group (target group) that must be retained.
                                    A sample will be removed if its value is -1.
    :param verbose: whether to print log to stderr
    :return:
    """
    matrix_bed = Path(matrix_bed)
    matrix_bed_csi = matrix_bed.with_suffix('.gz.csi')
    dmc = pd.read_csv(dmc_file, sep='\t')
    if dmc_file_filtered:
        fo = open(dmc_file_filtered, 'w')
    else:
        fo = sys.stdout
    tabix_file = pysam.TabixFile(str(matrix_bed), index=matrix_bed_csi)
    header = tabix_file.header[0]
    header_items = header.strip().split('\t')
    samples = pd.Series(header_items[3:])
    fo.write(header + "\n")
    i = 0  # the number of DMCs processed
    j = 0  # the number of DMCs passed the filter
    total = dmc.shape[0]  # the total number of DMCs
    start = time.time()
    for item in dmc.values:
        if i % 1000 == 0:
            passed_seconds = round(time.time() - start)
            processed_r = round(i / total, 4)
            keep_r = round(j / (i + 0.000001), 4)
            if verbose:
                sys.stderr.write(
                    f'Processed {i}/{total}={processed_r}\tkeep {j}/{i}={keep_r}\tpassed {passed_seconds}s\n')
        foldchange = item[10]
        iterator = tabix_file.fetch(region=f'{item[0]}:{item[1]}-{item[2]}')
        try:
            line = next(iterator)
        except StopIteration:
            continue
        line = line.strip().split('\t')
        full = pd.Series([float(x) for x in line[3:]])
        full[full == -1] = np.nan
        target = full[samples.isin(target_samples)]
        background = full[~samples.isin(target_samples)]
        target_keep = target[~target.isna()]
        background_keep = background[~background.isna()]
        if len(target) == 0 or len(background) == 0:
            raise Exception('No samples matched')
        if len(target_keep) / len(target) < minimum_sample_keep_ratio or len(background_keep) / len(
                background) < minimum_sample_keep_ratio:
            continue
        if foldchange < 0:
            left = target_keep.quantile(percentile)
            right = background_keep.quantile(1 - percentile)
        else:
            left = background_keep.quantile(percentile)
            right = target_keep.quantile(1 - percentile)
        if left < right:
            fo.write('\t'.join([str(x) for x in item]) + '\n')
            j = j + 1
        i += 1
    fo.close()


def get_args():
    parser = argparse.ArgumentParser(
        description='This tool is employed to further refine and filter the DMCs identified by MOABS:mcomp')
    parser.add_argument('-i', '--input', required=True, help=get_param_from_doc('dmc_file', mcomp_dmc_filter))
    parser.add_argument('-o', '--output', help=get_param_from_doc('dmc_file_filtered', mcomp_dmc_filter))
    parser.add_argument('-m', '--matrix-bed', required=True,
                        help=get_param_from_doc('matrix_bed', mcomp_dmc_filter))
    parser.add_argument('-t', '--target-samples', required=True,
                        help=get_param_from_doc('target_samples',
                                                mcomp_dmc_filter) + " The samples should split by ',', such as: sample1,sample2,sample3")
    parser.add_argument('-p', '--percentile', default=0.8, type=float,
                        help=get_param_from_doc('percentile', mcomp_dmc_filter))
    parser.add_argument('-k', '--minimum-keep', default=0.8, type=float,
                        help=get_param_from_doc('minimum_sample_keep_ratio', mcomp_dmc_filter))
    parser.add_argument('-v', '--verbose', default=True, help=get_param_from_doc('verbose', mcomp_dmc_filter))
    return parser.parse_args()


def main():
    args = get_args()
    target_samples = args.target_samples.split(',')
    mcomp_dmc_filter(args.input, args.matrix_bed, target_samples, args.output, args.percentile, args.minimum_keep,
                     args.verbose)


if __name__ == '__main__':
    main()
