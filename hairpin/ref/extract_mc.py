import multiprocessing
import argparse
from functools import partial
from multiprocessing import Process
import pysam
import collections
import subprocess
from Bio import SeqIO



def split_sam_file(chr,sam_file):
    split_command = f"samtools view -b -h -o {sam_file}.mc.{chr} {sam_file} {chr}"
    index_command = f"samtools index {sam_file}.mc.{chr}"
    subprocess.run(split_command, shell=True, check=True)
    subprocess.run(index_command, shell=True, check=True)

def process_file(bam_file,dd,ss):
    sam = pysam.AlignmentFile(bam_file, "r")
    hmc_bed = open(bam_file + ".bed", "w")
    for read in sam.fetch():
        if read.is_unmapped:
            continue
        cut_seq = ss[read.query_name]
        seq = read.seq
        aligned_pairs = read.get_aligned_pairs(with_seq=True, matches_only=False)
        if read.is_reverse:
            cut_seq=reverse_complement(cut_seq)
            for query_pos, ref_pos, ref_base in aligned_pairs:
                if ref_pos is None or query_pos is None:
                    continue
                elif seq[query_pos] == "G" and cut_seq[query_pos] == "G" and dd[read.reference_name][ref_pos] == "G":
                    get_pos = read.reference_name + " " + str(ref_pos + 1) + " " + str(ref_pos + 1) + " "
                    if ref_pos == 0:
                        print(get_pos + str(1) + " " + "CN", file=hmc_bed)
                    elif ref_pos == 1 and dd[read.reference_name][ref_pos - 1] != "C":
                        print(get_pos + str(1) + " " + "CN", file=hmc_bed)
                    elif dd[read.reference_name][ref_pos - 1] == "C":
                        print(get_pos + str(1) + " " + "CpG", file=hmc_bed)
                    elif dd[read.reference_name][ref_pos - 1] == "N" :
                        print(get_pos + str(1) + " " + "CN", file=hmc_bed)
                    elif dd[read.reference_name][ref_pos - 2] == "C":
                        print(get_pos + str(1) + " " + "CHG", file=hmc_bed)
                    elif dd[read.reference_name][ref_pos - 2] == "N" :
                        print(get_pos + str(1) + " " + "CN", file=hmc_bed)
                    else:
                        print(get_pos + str(1) + " " + "CHH", file=hmc_bed)
                elif (seq[query_pos] == "G" and ref_base == "G") or (seq[query_pos] == "A" and ref_base == "g"):
                    get_pos = read.reference_name + " " + str(ref_pos + 1) + " " + str(ref_pos + 1) + " "
                    if ref_pos == 0:
                        print(get_pos + str(0) + " " + "CN", file=hmc_bed)
                    elif ref_pos == 1 and dd[read.reference_name][ref_pos - 1] != "C":
                        print(get_pos + str(0) + " " + "CN", file=hmc_bed)
                    elif dd[read.reference_name][ref_pos - 1] == "C":
                        print(get_pos + str(0) + " " + "CpG", file=hmc_bed)
                    elif dd[read.reference_name][ref_pos - 1] == "N" :
                        print(get_pos + str(0) + " " + "CN", file=hmc_bed)
                    elif dd[read.reference_name][ref_pos - 2] == "C":
                        print(get_pos + str(0) + " " + "CHG", file=hmc_bed)
                    elif dd[read.reference_name][ref_pos - 2] == "N" :
                        print(get_pos + str(0) + " " + "CN", file=hmc_bed)
                    else:
                        print(get_pos + str(0) + " " + "CHH", file=hmc_bed)
        else:
            for query_pos, ref_pos, ref_base in aligned_pairs:
                if ref_pos is None or query_pos is None:
                    continue
                elif seq[query_pos] == "C" and cut_seq[query_pos] == "C" and dd[read.reference_name][ref_pos] == "C":
                    get_pos = read.reference_name + " " + str(ref_pos + 1) + " " + str(ref_pos + 1) + " "
                    if ref_pos == len(dd[read.reference_name])-1:
                        print(get_pos + str(1) + " " + "CN", file=hmc_bed)
                    elif ref_pos == len(dd[read.reference_name])-2 and dd[read.reference_name][ref_pos+1]!="G":
                        print(get_pos + str(1) + " " + "CN", file=hmc_bed)
                    elif dd[read.reference_name][ref_pos +1] == "G":
                        print(get_pos + str(1) + " " + "CpG", file=hmc_bed)
                    elif dd[read.reference_name][ref_pos +1] == "N" :
                        print(get_pos + str(1) + " " + "CN", file=hmc_bed)
                    elif dd[read.reference_name][ref_pos +2] == "G":
                        print(get_pos + str(1) + " " + "CHG", file=hmc_bed)
                    elif dd[read.reference_name][ref_pos +2] == "N" :
                        print(get_pos + str(1) + " " + "CN", file=hmc_bed)
                    else:
                        print(get_pos + str(1) + " " + "CHH", file=hmc_bed)
                elif (seq[query_pos] == "C" and ref_base == "C") or (seq[query_pos] == "T" and ref_base == "c"):
                    get_pos = read.reference_name + " " + str(ref_pos + 1) + " " + str(ref_pos + 1) + " "
                    if ref_pos == len(dd[read.reference_name])-1:
                        print(get_pos + str(0) + " " + "CN", file=hmc_bed)
                    elif ref_pos == len(dd[read.reference_name])-2 and dd[read.reference_name][ref_pos+1]!="G":
                        print(get_pos + str(0) + " " + "CN", file=hmc_bed)
                    elif dd[read.reference_name][ref_pos +1] == "G":
                        print(get_pos + str(0) + " " + "CpG", file=hmc_bed)
                    elif dd[read.reference_name][ref_pos +1] == "N" :
                        print(get_pos + str(0) + " " + "CN", file=hmc_bed)
                    elif dd[read.reference_name][ref_pos +2] == "G":
                        print(get_pos + str(0) + " " + "CHG", file=hmc_bed)
                    elif dd[read.reference_name][ref_pos +2] == "N" :
                        print(get_pos + str(0) + " " + "CN", file=hmc_bed)
                    else:
                        print(get_pos + str(0) + " " + "CHH", file=hmc_bed)
    sam.close()
    hmc_bed.close()

def final_bed(hmc_bed):
    my_bed = collections.defaultdict(lambda: (0, 0))
    with open(hmc_bed) as fin, open(hmc_bed+".merge", "w") as fout:
        for line in fin:
            line = line.strip().split()
            chrom, start, end, hmc, hmc_type = line[0], int(line[1]), int(line[2]), int(line[3]), line[4]
            key = (chrom, start, hmc_type)
            if hmc == 1:
                my_bed[key] = (my_bed[key][0] + 1, my_bed[key][1])
            elif hmc == 0:
                my_bed[key] = (my_bed[key][0], my_bed[key][1] + 1)

        # 将处理结果写入输出文件
        for key, counts in my_bed.items():
            chrom, start, hmc_type = key
            meth_ratio = counts[0] / (counts[0] + counts[1]) * 100
            meth_counts = counts[0]
            unmeth_counts = counts[1]
            fout.write(f"{chrom}\t{start}\t{start}\t{meth_ratio:.2f}\t{meth_counts}\t{unmeth_counts}\t{hmc_type}\n")


def reverse_complement(dna):
    complement = str.maketrans('ATCG', 'TAGC')
    reverse = dna[::-1].translate(complement)
    return reverse
def readgenome(name):
    dd = {}
    with open(name, 'r') as input_fasta:
        for record in SeqIO.parse(input_fasta, 'fasta'):
            dd[record.id] = str(record.seq)
    return dd


def ppp(partial_process_chunk_chr, chunks_chr, args):
    processes = []  # 存储进程对象的列表
    for file in chunks_chr:
        p = Process(target=partial_process_chunk_chr, args=(file,))
        processes.append(p)
    # 按照最大进程数限制逐个启动进程
    for i in range(0, len(processes), args.parallel):
        batch = processes[i:i + args.parallel]
        for p in batch:
            p.start()
        for p in batch:
            p.join()

def get_cut1_seq(fastq_file):
    seq_dict = {}
    with open(fastq_file) as f:
        while True:
            line = f.readline().strip()
            if not line:
                break  # 到达文件末尾，跳出循环
            seq_id = line.split()[0][1:]  # 去掉'@'字符
            seq = f.readline().strip()
            f.readline()  # 跳过不需要的行
            f.readline()  # 跳过不需要的行
            seq_dict[seq_id] = seq
    return seq_dict



def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--parallel",
        help="Number of parallel threads to use for processing. Default: 20",
        default=20,
        type=int,
    )

    parser.add_argument(
        "--ref",
        help="ref fasta",
        default="/data/wangzc/cgs/mc/new_mm9_lam_mc.fa",
        type=str,
    )

    parser.add_argument(
        "--sam",
        help="sam file",
        type=str,
    )

    parser.add_argument(
        "--output",
        help="output file",
        type=str,
    )

    parser.add_argument(
        "--cutfq1",
        help="output file",
        type=str,
    )

    args = parser.parse_args()
    if args.output is None:
        args.output = args.sam+".all.bed.mc"



    dd = readgenome(args.ref)
    ss = get_cut1_seq(args.cutfq1)
    chr_list = list(dd.keys())
    index_part_function = partial(split_sam_file, sam_file=args.sam)
    with multiprocessing.Pool(processes=args.parallel) as pool:
        pool.map(index_part_function, chr_list)

    partial_process_file = partial(process_file, dd=dd,ss=ss)
    ppp(partial_process_file, [args.sam+".mc" + "." + _ for _ in chr_list ], args)



    with multiprocessing.Pool(processes=args.parallel) as pool:
        pool.map(final_bed, [args.sam+".mc" + "." + _ + ".bed" for _ in chr_list ])



    command = f"cat {args.sam}.mc.*.bed.merge > {args.output}"
    subprocess.run(command, shell=True)

    rm_list=[args.sam+".mc" + "." + _ for _ in chr_list ] + [args.sam +".mc"+ "." + _ + ".bed" for _ in chr_list ] +\
            [args.sam +".mc"+ "." + _ + ".bed.merge" for _ in chr_list ] + [args.sam +".mc"+ "." + _ + ".bai" for _ in chr_list ]
    rm_str=" ".join(rm_list)
    delete_command = f"rm -rf {rm_str}"
    subprocess.run(delete_command, shell=True)




if __name__ == '__main__':
    main()












