import argparse

import numpy as np
import pandas as pd
from Bio import SeqIO


def get_lengths(filename):
    # 初始化序列长度列表
    lengths = []

    # 遍历fastq文件中的所有序列，并将它们的长度添加到列表中
    for record in SeqIO.parse(filename, "fastq"):
        lengths.append(len(record.seq))

    # 计算序列长度的统计信息
    mean_length = np.mean(lengths)
    median_length = np.median(lengths)
    std_length = np.std(lengths)

    # 将序列长度列表转换为Pandas Series对象
    lengths_series = pd.Series(lengths)

    # 计算每个长度的数量和占比
    counts = lengths_series.value_counts(normalize=True)

    # 获取占比最多的序列长度及其所占比例
    most_common_length = counts.index[0]
    most_common_length_percentage = counts.iloc[0]

    result_str = "Mean length: " + str(mean_length) + "\n" + \
                 "Median length: " + str(median_length) + "\n" + \
                 "Standard deviation: " + str(std_length) + "\n" + \
                 "Most common length: " + str(most_common_length) + "\n" + \
                 "Most common length percentage: " + str(most_common_length_percentage) + "\n"

    return result_str


def main():
    parser = argparse.ArgumentParser(description="Duplicate Base Kmer Tool")
    parser.add_argument("--file", type=str, help="Input fq1 file")
    parser.add_argument("--out", type=str, help="Input fq1 file")
    args = parser.parse_args()
    result_str=get_lengths(args.file)
    with open(args.out,"w")as file:
        print(result_str,file=file)


if __name__ == '__main__':
    main()
