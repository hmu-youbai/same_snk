import argparse
import os

from Bio import SeqIO

def readgenome(name):
    dd = {}
    with open(name, 'r') as input_fasta:
        for record in SeqIO.parse(input_fasta, 'fasta'):
            dd[record.id] = str(record.seq)
    return dd







def get_list(bed_file,num_start):
    list_for_yz=[]
    with open(bed_file) as file:
        temp = 0
        for line in file:
            line1=line.split()
            if line.startswith("#") or (not line1[num_start-1].isdigit()):
                continue
            temp = temp + 1
            list_for_yz.append(line.strip())
            if temp==200000:
                break

    if len(list_for_yz) != 200000:
        new_len = len(list_for_yz) - (len(list_for_yz) % 2)
    else:
        new_len = len(list_for_yz)

    return list_for_yz[:new_len], new_len





def start_yz(bed_file,dd,list_for_yz,new_len,num_chr,num_start,out_file):


    zuo1 = 0
    zuo0 = 0
    for i,line in enumerate(list_for_yz):
        line=line.split()
        if i < int(new_len/2):
            if dd[line[num_chr-1]][int(line[num_start-1]) - 1] in ["C", "G", "g", "c"]:
                zuo1 = zuo1 + 1
        elif i >= int(new_len/2):
            if dd[line[num_chr-1]][int(line[num_start-1])] in ["C", "G", "g", "c"]:
                zuo0 = zuo0 + 1



    if zuo1 == int(new_len/2) and zuo0 != int(new_len/2):
        print("参考基因组正确，bed文件start对应位点，start最小值从1开始，这是少见的情况，已纠正")
        with open(bed_file) as file, open(out_file, "w") as file2:
            for line in file:
                line1 = line.split()
                if line.startswith("#") or (not line1[num_start-1].isdigit()):
                    continue
                line = line.strip()
                line = line.split()
                new_line = '\t'.join(str(int(val) - 1) if idx == num_start - 1 else str(val) for idx, val in enumerate(line))

                if dd[line[num_chr - 1]][int(line[num_start - 1]) - 1] in ["C", "c"]:
                    new_line=new_line+"\t"+"+"
                    print(new_line, file=file2)
                if dd[line[num_chr - 1]][int(line[num_start - 1]) - 1] in ["G", "g"]:
                    new_line = new_line + "\t" + "-"
                    print(new_line, file=file2)



    elif zuo0 == int(new_len/2) and zuo1 != int(new_len/2):
        print("参考基因组正确，bed文件start对应位点，start最小值从0开始，这是最常见的情况")

        with open(bed_file) as file, open(out_file, "w") as file2:
            for line in file:
                line = line.strip()
                line1 = line.split()
                row_line=line
                if line.startswith("#") or (not line1[num_start-1].isdigit()):
                    continue
                line=line.split()
                if dd[line[num_chr-1]][int(line[num_start - 1]) ] in ["C", "c"]:
                    new_line=row_line+"\t"+"+"
                    print(new_line, file=file2)

                if dd[line[num_chr-1]][int(line[num_start - 1]) ] in ["G", "g"]:
                    new_line=row_line+"\t"+"-"
                    print(new_line, file=file2)


    elif zuo0 != int(new_len/2) and zuo1 != int(new_len/2):
        print("参考基因组错误")
    elif zuo0 == int(new_len/2) and zuo1 == int(new_len/2):
        print("错误，参考基因组正确，但是起点模式待定，不排除全正链的可能")
    else:
        print("错误，啊？亚托克斯音效！")







if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Description of your program')

    # 添加参数
    parser.add_argument('bed_file', type=str, help='Description of bed_file')
    parser.add_argument('ref', type=str, help='Description of dd')


    parser.add_argument('num_chr', type=int, help='Description of num_chr')
    parser.add_argument('num_start', type=int, help='Description of num_start')


    parser.add_argument('out_file', type=str, help='Description of dd')
    # 解析命令行参数
    args = parser.parse_args()


    dd = readgenome(args.ref)
    list_for_yz,new_len=get_list(args.bed_file,args.num_start)

    start_yz(args.bed_file, dd, list_for_yz, new_len, args.num_chr, args.num_start,args.out_file)


