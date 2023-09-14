import argparse





if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="get pairs")

    parser.add_argument("merge", help="第一个文件路径")
    parser.add_argument("mc_file", help="第二个文件路径")
    parser.add_argument("c_file",  help="sam1")

    args = parser.parse_args()


    substring1="TCTCTTATACAC"
    substring2="TTTTTTATATAT"
    
    mode=""
    with open(args.mc_file,"w")as file1,open(args.c_file,"w")as file2,open(args.merge)as file3:
        for i,line in enumerate(file3):
            if i %4==0:
                name=line.strip()
            if i%4==1:
                seq=line.strip()
            if i%4==3:
                q=line.strip()
                if seq[260:].find(substring1) != -1:
                    mode="mc"
                elif seq[260:].find(substring2) != -1:
                    mode="c"
    
                else:         # 计算一下，300*0.03=9
                    if seq.count("C") > len(seq)*0.03:
                        mode = "mc"
                    else:
                        mode = "c"
    
    
                if mode=="mc":
                    print(name,file=file1)
                    print(seq, file=file1)
                    print("+", file=file1)
                    print(q, file=file1)
    
                else:
                    print(name,file=file2)
                    print(seq, file=file2)
                    print("+", file=file2)
                    print(q, file=file2)

