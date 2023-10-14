import Levenshtein
import gzip
import math
from functools import partial
import os
from collections import namedtuple
import multiprocessing
from Bio.Align.substitution_matrices import Array
from Bio.Align import PairwiseAligner
import numpy as np
import argparse

# rule_1=fill mC   rule_2=fill C
# 这里额外保留了GG的情况，考虑到了BS的假阴性的情况，GG认为是G
class Alignment:
    alphabet = "ACGTN"
    rule_matrix_1 = np.array([
        [1.0, 0.0, 0.0, 0.0, 0.0],
        [0.0, 1.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0, 1.0, 0.0],
        [0.0, 0.0, 0.0, 0.0, 0.0]
    ])
    rule_matrix_2 = np.array([
        [1.0, 0.0, 0.0, 0.0, 0.0],
        [0.0, 1.0, 0.0, 0.0, 0.0],
        [1.0, 0.0, 1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0, 1.0, 0.0],
        [0.0, 0.0, 0.0, 0.0, 0.0]
    ])

    def __init__(self, sub_matrix, mode='global', open_gap_score=-3, extend_gap_score=-2):
        aligner = PairwiseAligner()
        aligner.mode = mode
        aligner.open_gap_score = open_gap_score
        aligner.extend_gap_score = extend_gap_score
        if sub_matrix == 'rule_matrix_1':
            aligner.substitution_matrix = Array(self.alphabet, 2, self.rule_matrix_1)
        elif sub_matrix == 'rule_matrix_2':
            aligner.substitution_matrix = Array(self.alphabet, 2, self.rule_matrix_2)
        else:
            raise ValueError(f"Invalid sub_matrix argument: {sub_matrix}")
        self.aligner = aligner
        self.type = sub_matrix


    def align(self, seq1, seq2):
        alignments = self.aligner.align(seq1[::-1],seq2[::-1])
        align_read1, align_read2 = alignments[0][0][::-1], alignments[0][1][::-1]
        if self.type == 'rule_matrix_1':
            score_rules = {
                ('A', 'A'): 1, ('A', 'C'): -10, ('A', 'G'): -10, ('A', 'T'): -10, ('A', 'N'): -10, ('A', '-'): -10,
                ('C', 'A'): -10, ('C', 'C'): 1, ('C', 'G'): -10, ('C', 'T'): -10, ('C', 'N'): -10, ('C', '-'): -10,
                ('G', 'A'): -20, ('G', 'C'): -10, ('G', 'G'): 1, ('G', 'T'): -10, ('G', 'N'): -10, ('G', '-'): -10,
                ('T', 'A'): -10, ('T', 'C'): 1, ('T', 'G'): -10, ('T', 'T'): 1, ('T', 'N'): -10, ('T', '-'): -10,
                ('N', 'A'): -10, ('N', 'C'): -10, ('N', 'G'): -10, ('N', 'T'): -10, ('N', 'N'): -10, ('N', '-'): -10,
                ('-', 'A'): -10, ('-', 'C'): -10, ('-', 'G'): -10, ('-', 'T'): -10, ('-', 'N'): -10, ('-', '-'): -10
            }

            scores = [score_rules.get((c1, c2)) for c1, c2 in zip(align_read1, align_read2)]
            max_index = max(range(len(scores)), key=lambda i: sum(scores[:i + 1])) + 1
            align_read1 = align_read1[:max_index]
            align_read2 = align_read2[:max_index]

            result = [c2 if ((c1 == c2) or (c1 == 'T' and c2 == 'C')) else 'N' for c1, c2 in zip(align_read1, align_read2)]

        elif self.type == 'rule_matrix_2':
            score_rules = {
                ('A', 'A'): 1, ('A', 'C'): -10, ('A', 'G'): -10, ('A', 'T'): -10, ('A', 'N'): -10, ('A', '-'): -10,
                ('C', 'A'): -10, ('C', 'C'): 1, ('C', 'G'): -10, ('C', 'T'): -10, ('C', 'N'): -10, ('C', '-'): -10,
                ('G', 'A'): 1, ('G', 'C'): -10, ('G', 'G'): 1, ('G', 'T'): -10, ('G', 'N'): -10, ('G', '-'): -10,
                ('T', 'A'): -10, ('T', 'C'): 1, ('T', 'G'): -10, ('T', 'T'): 1, ('T', 'N'): -10, ('T', '-'): -10,
                ('N', 'A'): -10, ('N', 'C'): -10, ('N', 'G'): -10, ('N', 'T'): -10, ('N', 'N'): -10, ('N', '-'): -10,
                ('-', 'A'): -10, ('-', 'C'): -10, ('-', 'G'): -10, ('-', 'T'): -10, ('-', 'N'): -10,  ('-', '-'): -10
            }
            scores = [score_rules.get((c1, c2)) for c1, c2 in zip(align_read1, align_read2)]
            max_index = max(range(len(scores)), key=lambda i: sum(scores[:i + 1])) + 1
            align_read1 = align_read1[:max_index]
            align_read2 = align_read2[:max_index]

            result = [c2 if ((c1 == c2) or (c1 == 'T' and c2 == 'C') ) else c1 if (c1 == 'G' and c2 == 'A') else 'N' for
                      c1, c2 in zip(align_read1, align_read2)]
        else:
            raise ValueError(f"Invalid sub_matrix argument: {self.type}")

        return ''.join(result), align_read1, align_read2



Record = namedtuple('Record', ['name', 'seq', 'qual', 'seq1', 'seq2'])
def complement_dna(dna):
    complement = str.maketrans('ATCG', 'TAGC')
    reverse = dna[::-1].translate(complement)
    return reverse

def check(str1, str2):
    length1 = len(str1)
    length = min(length1, len(str2))
    k = max(range(length + 1),
            key=lambda i: i if Levenshtein.hamming(str1[length1 - i:], str2[:i]) < i * 0.1 else False)
    return k

def trim_overlap(read1, read2):
    a1 = check(read1, complement_dna(read2))
    if a1 > 5:
        return read1[:-int(a1 / 2)], read2[:-int(a1 / 2)]
    else:
        a2 = check(complement_dna(read2), read1)
        if a2 > min(len(read2), len(read1)) * 0.8 and a2 > 5:
            return read1[:-int(a2 / 2)], read2[:-int(a2 / 2)]
    return read1, read2

def open_input_file(input_file):
    if input_file.endswith('.gz'):
        in_handle = gzip.open(input_file, "rt")
    else:
        in_handle = open(input_file)
    return in_handle

def correct_seq(seq, seq1, seq2, qual):
    q1=qual.split()[0]
    q2=qual.split()[1]

    temp_result = []
    temp = 0
    for char in seq1:
        if char == "-":
            temp_result.append(",")
        else:
            temp_result.append(q1[temp])
            temp += 1
    str_q1 = "".join(temp_result)

    temp_result = []
    temp = 0
    for char in seq2:
        if char == "-":
            temp_result.append(",")
        else:
            temp_result.append(q2[temp])
            temp += 1
    str_q2 = "".join(temp_result)

    # 可选部分，通过质量比较纠正seq
    c_qq=[]
    result_seq=list(seq)
    for index,char in enumerate(seq):
        if char=="N":
            c_qq.append(",")
            if str_q1[index]>str_q2[index] and seq1!="-" and seq1!="T":
                result_seq[index] = seq1[index]
            elif str_q2[index]>str_q1[index] and seq2!="-" and seq2!="A":
                result_seq[index] = seq2[index]
        else:
            c_qq.append(max(str_q1[index],str_q2[index]))
    c_qq = "".join(c_qq)
    return "".join(result_seq), c_qq+" "+str_q1+" "+str_q2


def process_chunk(chunks, alignment, min_lenth=50):
    result = []
    for record in chunks:
        read1, read2 = trim_overlap(record.seq.split()[0], record.seq.split()[1])
        seq, seq1, seq2 = alignment.align(read1, read2)
        if len(seq) < min_lenth:
            continue
        c_seq,c_q=correct_seq(seq, seq1, seq2, record.qual)

        result.append(Record(name=record.name, seq=c_seq.replace('-', 'N'), qual=c_q, seq1=seq1.replace('-', 'N'), seq2=seq2.replace('-', 'N')))
    return result
def process_records(buffer, alignment, min_lenth=50, processes=24):
    pool = multiprocessing.Pool(processes=processes)
    chunk_size = math.ceil(len(buffer) / processes)
    chunks = [buffer[i:i + chunk_size] for i in range(0, len(buffer), chunk_size)]
    partial_process_chunk = partial(process_chunk, alignment=alignment, min_lenth=min_lenth)
    results = pool.map(partial_process_chunk, chunks)
    pool.close()
    pool.join()
    return [record for result in results for record in result]

def cache_and_process(input_file1, input_file2, output_file, alignment, min_lenth=50, processes=24, chunk_size=1000000):
    # output_dir = os.path.dirname(output_file)
    # cut_f1_path = os.path.join(output_dir, output_file[:-2]+'cut_f1.fq')
    # cut_f2_path = os.path.join(output_dir, output_file[:-2]+'cut_f2.fq')
    cut_f1_path = output_file +'_cut_f1.fq'
    cut_f2_path = output_file +'_cut_f2.fq'



    report_n1 = 0

    in_handle1 = open_input_file(input_file1)
    in_handle2 = open_input_file(input_file2)

    with in_handle1, in_handle2, open(output_file+".fq", "w") as out_handle, open(cut_f1_path, "w") as f1, open(cut_f2_path, "w") as f2:


        buffer = []
        for line_num, (line1, line2) in enumerate(zip(in_handle1, in_handle2)):
            if line_num % 4 == 0:
                record_name1 = line1.strip().split()[0]
            elif line_num % 4 == 1:
                record_seq = line1.strip() + " " + line2.strip()
            elif line_num % 4 == 3:
                record_qual1 = line1.strip()
                record_qual2 = line2.strip()
                buffer.append(Record(name=record_name1, seq=record_seq, qual=record_qual1+" "+record_qual2, seq1="nothing",
                                     seq2="nothing"))
            if len(buffer) == chunk_size:
                print("process reads : " + str(int((line_num + 1) / 4)))
                temp = process_records(buffer, alignment, min_lenth, processes)
                report_n1 = report_n1 + len(temp)
                for seq in temp:
                    out_handle.write(f"{seq.name}\n")
                    out_handle.write(f"{seq.seq}\n")
                    out_handle.write(f"+\n")
                    out_handle.write(f"{seq.qual.split()[0]}\n")
                    f1.write(f"{seq.name}\n")
                    f1.write(f"{seq.seq1}\n")
                    f1.write(f"+\n")
                    f1.write(f"{seq.qual.split()[1]}\n")
                    f2.write(f"{seq.name}\n")
                    f2.write(f"{seq.seq2}\n")
                    f2.write(f"+\n")
                    f2.write(f"{seq.qual.split()[2]}\n")
                buffer = []
        if buffer:
            temp = process_records(buffer, alignment, min_lenth, processes)
            report_n1 = report_n1 + len(temp)
            for seq in temp:
                out_handle.write(f"{seq.name}\n")
                out_handle.write(f"{seq.seq}\n")
                out_handle.write(f"+\n")
                out_handle.write(f"{seq.qual.split()[0]}\n")
                f1.write(f"{seq.name}\n")
                f1.write(f"{seq.seq1}\n")
                f1.write(f"+\n")
                f1.write(f"{seq.qual.split()[1]}\n")
                f2.write(f"{seq.name}\n")
                f2.write(f"{seq.seq2}\n")
                f2.write(f"+\n")
                f2.write(f"{seq.qual.split()[2]}\n")



    message = "fq1: " + input_file1 + "\n" + "fq2: " + input_file2 + "\n" + "min_read_length: " + str(
        min_lenth) + "\n" + "input_reads: " + str(int((line_num + 1) / 4)) + "\n" + "resolved_reads: " + str(
        report_n1) + "\n" + "resolved ratio: " + str("{:.2f}%".format(report_n1 / int((line_num + 1) / 4) * 100))
    print(message)
    return message




# 定义mC版本的函数
def correct_seq_mC(seq, seq1, seq2, qual):
    q1=qual.split()[0]
    q2=qual.split()[1]

    temp_result = []
    temp = 0
    for char in seq1:
        if char == "-":
            temp_result.append(",")
        else:
            temp_result.append(q1[temp])
            temp += 1
    str_q1 = "".join(temp_result)

    temp_result = []
    temp = 0
    for char in seq2:
        if char == "-":
            temp_result.append(",")
        else:
            temp_result.append(q2[temp])
            temp += 1
    str_q2 = "".join(temp_result)

    # 可选部分，通过质量比较纠正seq
    c_qq=[]
    result_seq=list(seq)
    for index,char in enumerate(seq):
        if char=="N":
            c_qq.append(",")
            if str_q1[index]>str_q2[index] and seq1!="-" and seq1!="T":
                result_seq[index] = seq1[index]
            elif str_q2[index]>str_q1[index] and seq2!="-" :
                result_seq[index] = seq2[index]
        else:
            c_qq.append(max(str_q1[index],str_q2[index]))
    c_qq = "".join(c_qq)
    return "".join(result_seq), c_qq+" "+str_q1+" "+str_q2
def process_chunk_mC(chunks, alignment, min_lenth=50):
    result = []
    for record in chunks:
        read1, read2 = trim_overlap(record.seq.split()[0], record.seq.split()[1])
        seq, seq1, seq2 = alignment.align(read1, read2)
        if len(seq) < min_lenth:
            continue
        c_seq,c_q=correct_seq_mC(seq, seq1, seq2, record.qual)

        result.append(Record(name=record.name, seq=c_seq.replace('-', 'N'), qual=c_q, seq1=seq1.replace('-', 'N'), seq2=seq2.replace('-', 'N')))
    return result
def process_records_mC(buffer, alignment, min_lenth=50, processes=24):
    pool = multiprocessing.Pool(processes=processes)
    chunk_size = math.ceil(len(buffer) / processes)
    chunks = [buffer[i:i + chunk_size] for i in range(0, len(buffer), chunk_size)]
    partial_process_chunk = partial(process_chunk_mC, alignment=alignment, min_lenth=min_lenth)
    results = pool.map(partial_process_chunk, chunks)
    pool.close()
    pool.join()
    return [record for result in results for record in result]

def cache_and_process_mC(input_file1, input_file2, output_file, alignment, min_lenth=50, processes=24, chunk_size=1000000):
    # output_dir = os.path.dirname(output_file)
    # cut_f1_path = os.path.join(output_dir, output_file[:-2]+'cut_f1.fq')
    # cut_f2_path = os.path.join(output_dir, output_file[:-2]+'cut_f2.fq')
    cut_f1_path = output_file +'_cut_f1.fq'
    cut_f2_path = output_file +'_cut_f2.fq'



    report_n1 = 0

    in_handle1 = open_input_file(input_file1)
    in_handle2 = open_input_file(input_file2)

    with in_handle1, in_handle2, open(output_file+".fq", "w") as out_handle, open(cut_f1_path, "w") as f1, open(cut_f2_path, "w") as f2:


        buffer = []
        for line_num, (line1, line2) in enumerate(zip(in_handle1, in_handle2)):
            if line_num % 4 == 0:
                record_name1 = line1.strip().split()[0]
            elif line_num % 4 == 1:
                record_seq = line1.strip() + " " + line2.strip()
            elif line_num % 4 == 3:
                record_qual1 = line1.strip()
                record_qual2 = line2.strip()
                buffer.append(Record(name=record_name1, seq=record_seq, qual=record_qual1+" "+record_qual2, seq1="nothing",
                                     seq2="nothing"))
            if len(buffer) == chunk_size:
                print("process reads : " + str(int((line_num + 1) / 4)))
                temp = process_records_mC(buffer, alignment, min_lenth, processes)
                report_n1 = report_n1 + len(temp)
                for seq in temp:
                    out_handle.write(f"{seq.name}\n")
                    out_handle.write(f"{seq.seq}\n")
                    out_handle.write(f"+\n")
                    out_handle.write(f"{seq.qual.split()[0]}\n")
                    f1.write(f"{seq.name}\n")
                    f1.write(f"{seq.seq1}\n")
                    f1.write(f"+\n")
                    f1.write(f"{seq.qual.split()[1]}\n")
                    f2.write(f"{seq.name}\n")
                    f2.write(f"{seq.seq2}\n")
                    f2.write(f"+\n")
                    f2.write(f"{seq.qual.split()[2]}\n")
                buffer = []
        if buffer:
            temp = process_records_mC(buffer, alignment, min_lenth, processes)
            report_n1 = report_n1 + len(temp)
            for seq in temp:
                out_handle.write(f"{seq.name}\n")
                out_handle.write(f"{seq.seq}\n")
                out_handle.write(f"+\n")
                out_handle.write(f"{seq.qual.split()[0]}\n")
                f1.write(f"{seq.name}\n")
                f1.write(f"{seq.seq1}\n")
                f1.write(f"+\n")
                f1.write(f"{seq.qual.split()[1]}\n")
                f2.write(f"{seq.name}\n")
                f2.write(f"{seq.seq2}\n")
                f2.write(f"+\n")
                f2.write(f"{seq.qual.split()[2]}\n")



    message = "fq1: " + input_file1 + "\n" + "fq2: " + input_file2 + "\n" + "min_read_length: " + str(
        min_lenth) + "\n" + "input_reads: " + str(int((line_num + 1) / 4)) + "\n" + "resolved_reads: " + str(
        report_n1) + "\n" + "resolved ratio: " + str("{:.2f}%".format(report_n1 / int((line_num + 1) / 4) * 100))
    print(message)
    return message












def main():
    parser = argparse.ArgumentParser(description="Duplicate Base Kmer Tool")

    parser.add_argument("--fq1", type=str, help="Input fq1 file")
    parser.add_argument("--fq2", type=str, help="Input fq2 file")
    parser.add_argument("--outfile", type=str,default="restored_seq", help="Output w file")
    parser.add_argument("--rule", type=int, default=2, help="default: 2")
    parser.add_argument("--min", type=int, default=50, help="Position parameter")
    parser.add_argument("--parallel", type=int, default=24, help="default: 24")
    parser.add_argument("--chunk_size", type=int, default=1000000, help="default: 1000000")

    args = parser.parse_args()

    if args.rule == 2:
        alignment = Alignment("rule_matrix_2")
        mylog = cache_and_process(args.fq1, args.fq2, args.outfile, alignment, args.min, args.parallel, args.chunk_size)
    else:
        alignment = Alignment("rule_matrix_1")
        mylog = cache_and_process_mC(args.fq1, args.fq2, args.outfile, alignment, args.min, args.parallel, args.chunk_size)


    print(mylog)


if __name__ == '__main__':
    main()



