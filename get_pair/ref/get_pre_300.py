import pysam
import gzip
import argparse


def complement_dna(dna):
    complement = str.maketrans('ATCG', 'TAGC')
    reverse = dna[::-1].translate(complement)
    return reverse


def get_mapped_read_names(input_bam):
    bam = pysam.AlignmentFile(input_bam)
    mapped_read_names = [read.query_name for read in bam.fetch() if not read.is_unmapped]
    bam.close()
    return mapped_read_names


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="get pairs")

    parser.add_argument("file1", help="第一个文件路径")
    parser.add_argument("file2", help="第二个文件路径")
    parser.add_argument("sam1", type=str, help="sam1")
    parser.add_argument("sam2", type=str, help="sam2")
    parser.add_argument("outfile", help="outfile")

    args = parser.parse_args()



    list1 = get_mapped_read_names(args.sam1)
    list2 = get_mapped_read_names(args.sam2)
    list_merge = set(list1 + list2)
    
    
    with gzip.open(args.file1, 'rb') as file1, gzip.open(args.file2, 'rb') as file2, open(args.outfile, "w") as file3:
        for i, line in enumerate(zip(file1, file2)):
            if i % 4 == 0:
                name = (line[0].split()[0]).decode('utf-8')
            if i % 4 == 1:
                read1 = (line[0].strip()).decode('utf-8')
                read2 = (line[1].strip()).decode('utf-8')
                read = read1 + complement_dna(read2)
            if i % 4 == 3:
                q1 = (line[0].strip()).decode('utf-8')
                q2 = (line[1].strip()).decode('utf-8')
                qq = q1 + q2[::-1]
                if name[1:] not in list_merge:
                    print(name, file=file3)
                    print(read, file=file3)
                    print("+", file=file3)
                    print(qq, file=file3)






