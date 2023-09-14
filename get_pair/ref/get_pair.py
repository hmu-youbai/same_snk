import multiprocessing
from functools import partial
import argparse
import time
from collections import deque
from collections import defaultdict
from collections import Counter
import multiprocessing



import math
import Levenshtein
import numpy as np
import math



def find_dup_base_distance(my_list,len_mc,dis=4,num_for_diff=10):
    temp_list,duplicate_indices,len1=merge_lists(my_list[:len_mc],my_list[len_mc:])

    result_list=[]
    for _ in duplicate_indices:
        result_list.append((temp_list[_], temp_list[_]))

    if len(temp_list)<30:
        all_sci_dir = defaultdict(list)
        for i in range(len(temp_list)):
            for j in range(i+1,len(temp_list)):
                if Levenshtein.distance(temp_list[i], temp_list[j]) < dis:
                    all_sci_dir[str(i)].append(j)

        if len(all_sci_dir) != 0:
            for temp_key in all_sci_dir.keys():
                temp_read1 = temp_list[int(temp_key)]
                for temp_index in all_sci_dir[temp_key]:
                    temp_read2 = temp_list[temp_index]

                    if int(temp_key) < len1 and temp_index < len1:
                        if int(temp_key) in duplicate_indices:
                            result_list.append((temp_read2, temp_read1))
                        elif temp_index in duplicate_indices:
                            result_list.append((temp_read1, temp_read2))

                    elif int(temp_key) < len1 and temp_index >= len1:
                        result_list.append((temp_read1, temp_read2))
        return result_list


    all_del_pos_list = np.triu(np.ones((len(temp_list), len(temp_list)), dtype=bool), k=1)
    all_sci_dir = {}

    for read_num in range(num_for_diff):
        temp_line = np.array([Levenshtein.distance(temp_list[read_num], s2) for s2 in temp_list], dtype=np.int16)
        temp_condition = np.abs(np.subtract.outer(temp_line, temp_line)) < dis
        all_del_pos_list &= temp_condition
        temp_aa = np.where(temp_line < dis)[0]
        if len(temp_aa) != 0:
            all_sci_dir[str(read_num)] = temp_aa

    for _ in range(num_for_diff):
        temp_arr=all_sci_dir[str(_)].copy()
        for x, y in enumerate(temp_arr):
            if y > _:
                all_sci_dir[str(_)]=temp_arr[x:]
                break

    for read_num in range(num_for_diff, len(temp_list)):
        temp_index = np.where(all_del_pos_list[read_num] == 1)[0].astype(np.uint32)
        temp_line = np.zeros(len(temp_index), dtype=np.int16)
        for i, my_index in enumerate(temp_index):
            temp_line[i] = Levenshtein.distance(temp_list[read_num], temp_list[my_index])

        temp_aa = np.where(temp_line < dis)[0]
        if len(temp_aa) != 0:
            all_sci_dir[str(read_num)] = temp_index[temp_aa]

    if len(all_sci_dir) != 0:
        for temp_key in all_sci_dir.keys():
            temp_read1 = temp_list[int(temp_key)]
            for temp_index in all_sci_dir[temp_key]:
                temp_read2 = temp_list[temp_index]

                if int(temp_key) < len1 and temp_index < len1:
                    if int(temp_key) in duplicate_indices:
                        result_list.append((temp_read2, temp_read1))
                    elif temp_index in duplicate_indices:
                        result_list.append((temp_read1, temp_read2))

                elif int(temp_key) < len1 and temp_index >= len1:
                    result_list.append((temp_read1, temp_read2))




    return result_list



def merge_lists(list1, list2):
    index_dict = {item: i for i, item in enumerate(list1)}
    merged_list = list1.copy()
    duplicate_indices = []
    for item in list2:
        if item not in index_dict:
            merged_list.append(item)
        else:
            duplicate_indices.append(index_dict[item])
    return merged_list, set(duplicate_indices), len(list1)

def find_dup_base_distance_para(my_list,dis=4,num_for_diff=10):
    len1=int(my_list[-1])
    my_list=my_list[:-1]
    temp_list,duplicate_indices,my_len1=merge_lists(my_list[:len1],my_list[len1:])

    result_list=[]
    for _ in duplicate_indices:
        result_list.append((temp_list[_], temp_list[_]))

    if len(temp_list)<30:
        all_sci_dir = defaultdict(list)
        for i in range(len(temp_list)):
            for j in range(i+1,len(temp_list)):
                if Levenshtein.distance(temp_list[i], temp_list[j]) < dis:
                    all_sci_dir[str(i)].append(j)

        if len(all_sci_dir) != 0:
            for temp_key in all_sci_dir.keys():
                temp_read1 = temp_list[int(temp_key)]
                for temp_index in all_sci_dir[temp_key]:
                    temp_read2 = temp_list[temp_index]

                    if int(temp_key) < len1 and temp_index < len1:
                        if int(temp_key) in duplicate_indices:
                            result_list.append((temp_read2, temp_read1))
                        elif temp_index in duplicate_indices:
                            result_list.append((temp_read1, temp_read2))

                    elif int(temp_key) < len1 and temp_index >= len1:
                        result_list.append((temp_read1, temp_read2))
        return result_list


    all_del_pos_list = np.triu(np.ones((len(temp_list), len(temp_list)), dtype=bool), k=1)
    all_sci_dir = {}

    for read_num in range(num_for_diff):
        temp_line = np.array([Levenshtein.distance(temp_list[read_num], s2) for s2 in temp_list], dtype=np.int16)
        temp_condition = np.abs(np.subtract.outer(temp_line, temp_line)) < dis
        all_del_pos_list &= temp_condition
        temp_aa = np.where(temp_line < dis)[0]
        if len(temp_aa) != 0:
            all_sci_dir[str(read_num)] = temp_aa

    for _ in range(num_for_diff):
        temp_arr=all_sci_dir[str(_)].copy()
        for x, y in enumerate(temp_arr):
            if y > _:
                all_sci_dir[str(_)]=temp_arr[x:]
                break

    for read_num in range(num_for_diff, len(temp_list)):
        temp_index = np.where(all_del_pos_list[read_num] == 1)[0].astype(np.uint32)
        temp_line = np.zeros(len(temp_index), dtype=np.int16)
        for i, my_index in enumerate(temp_index):
            temp_line[i] = Levenshtein.distance(temp_list[read_num], temp_list[my_index])

        temp_aa = np.where(temp_line < dis)[0]
        if len(temp_aa) != 0:
            all_sci_dir[str(read_num)] = temp_index[temp_aa]

    if len(all_sci_dir) != 0:
        for temp_key in all_sci_dir.keys():
            temp_read1 = temp_list[int(temp_key)]
            for temp_index in all_sci_dir[temp_key]:
                temp_read2 = temp_list[temp_index]

                if int(temp_key) < len1 and temp_index < len1:
                    if int(temp_key) in duplicate_indices:
                        result_list.append((temp_read2, temp_read1))
                    elif temp_index in duplicate_indices:
                        result_list.append((temp_read1, temp_read2))

                elif int(temp_key) < len1 and temp_index >= len1:
                    result_list.append((temp_read1, temp_read2))




    return result_list



def process_records(all_read, len1, dis=4, num_for_diff=10, parallel=20, memory_limit=50):
    pool = multiprocessing.Pool(processes=parallel)

    max_matrix_size = min(20000, math.ceil(10000 * (math.sqrt((memory_limit / parallel) / 0.25))))
    row_size1 = min(max_matrix_size, math.ceil(len(all_read) / parallel))
    all_read1 = all_read[:len1]
    all_read2 = all_read[len1:]



    chunks1 = [all_read1[i:i + row_size1] for i in range(0, len(all_read1), row_size1)]
    chunks2 = [all_read2[i:i + row_size1] for i in range(0, len(all_read2), row_size1)]

    chunks = []
    for l1 in chunks1:
        for l2 in chunks2:
            ll=l1 + l2
            ll.append(str(len(l1)))
            chunks.append(ll)


    partial_process_chunk = partial(find_dup_base_distance_para,dis=dis, num_for_diff=num_for_diff)
    results = pool.map(partial_process_chunk, chunks)

    pool.close()
    pool.join()

    return [record for result in results for record in result]



# ---------------------------------------------------read——file---------------------------
def sliding_window_optimal(string, k):
    order = {'G': 0, 'A': 1, 'T': 2, 'N': 5}

    window = deque(i for i in string[:k])
    min_substring = window.copy()
    for i in range(k, len(string)):
        window.append(string[i])
        window.popleft()

        for char1, char2 in zip(min_substring, window):
            if order[char1] > order[char2]:
                min_substring = window.copy()
                break
            elif order[char1] < order[char2]:
                break
    return ''.join(min_substring)

def get_fastq(fq_file, replace_C = True):
    my_dir=defaultdict(list)
    with open(fq_file) as file1:
        for i,line in enumerate(file1):
            if i % 4==0:
                name=line.split()[0]
            if i % 4==1:
                if replace_C:
                    read = line.strip().replace("C", "T")
                else:
                    read=line.strip()
                my_dir[read].append(name)
    result_dir={}
    dup_dir={}
    for key, value in my_dir.items():
        if len(value) != 1:
            dup_dir[value[0]] = value[1:]
        result_dir[key] = value[0]
    return result_dir,dup_dir
def get_kmer_dir0(fq_dir,start=0,w=30 ,k=26):
    my_dir = defaultdict(list)
    for key,value in fq_dir.items():
        my_dir[sliding_window_optimal(key[start:start + w], k)].append(key)
    return my_dir


def get_kmer_dir_mc(fq_dir,start=0,w=30,indel=2):
    # indel=int(dis/2)

    my_dir = defaultdict(list)
    for key, value in fq_dir.items():
        my_dir[key[start+indel:start+w-indel]].append(key)
    return my_dir
def get_kmer_dir_c(fq_dir, start=0,w=30,indel=2):

    my_dir = defaultdict(list)
    for key, value in fq_dir.items():
        my_dir[key[start + indel:start + w - indel]].append(key)
        for i in range(1,indel+1):
            my_dir[key[start + indel+i:start + w - indel+i]].append(key)
            my_dir[key[start + indel-i:start + w - indel-i]].append(key)
    return my_dir

# 这里解释一下，仅支持两个错误，最多2处插入+2处缺失，暂时不支持修改
# 这里提供了可调整参数start,w  受限于任务的复杂情况，手动设置这些参数






# ---------------------------------------------main--------------------------------------------
# dis参数有效
def start_get_pair_in_group_dis(mc_dir,c_dir,mc_fq,c_fq,dis):
    print("start")
    result_list = []
    pp=0
    for key in mc_dir.keys():
        pp=pp+1
        if pp%1000==0:
            print(pp)
        if key in c_dir:
            list1 = mc_dir[key]
            list2 = c_dir[key]

            if len(list1) == len(list2) == 1:
                if Levenshtein.distance(list1[0], list2[0]) < dis:
                    result_list.append((mc_fq[list1[0]], c_fq[list2[0]]))

            elif len(list1) + len(list2) > 1000:
                my_list=list1+list2
                temp_result = process_records(my_list, len(list1),dis)
                for i in temp_result:
                    result_list.append((mc_fq[i[0]],c_fq[i[1]]))


            else:
                my_list = list1 + list2
                temp_result=find_dup_base_distance(my_list,len(list1),dis)
                for i in temp_result:
                    result_list.append((mc_fq[i[0]],c_fq[i[1]]))


    return  result_list








if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Duplicate Base Kmer Tool")

    parser.add_argument("file1", help="第一个mC文件路径")
    parser.add_argument("file2", nargs='?', help="第二个C文件路径（可选）")
    parser.add_argument("w_file", type=str, help="Output w file")
    parser.add_argument("pos", type=int, help="Position parameter")
    parser.add_argument("--dis", type=int, default=4, help="Distance parameter (default: 4)")
    parser.add_argument("--windows", type=int, default=30, help="Distance parameter (default: 30)")

    parser.add_argument("--replace_C", type=int, default=1, help=" (default: 1)")





    args = parser.parse_args()



    if args.file2:
        mc_fq, mc_dup = get_fastq(args.file1, args.replace_C)
        mc_dir = get_kmer_dir_mc(mc_fq, args.pos,args.windows)
        c_fq, c_dup = get_fastq(args.file2, args.replace_C)
        c_dir = get_kmer_dir_c(c_fq, args.pos,args.windows)

        my_list = start_get_pair_in_group_dis(mc_dir, c_dir, mc_fq, c_fq, args.dis)
        with open(args.w_file, "w") as file:
            for i in my_list:
                print(i[0] + " " + i[1], file=file)

    else:
        pass


        # mc_dir = get_kmer_dir(mc_fq, args.pos)
        # c_dir = get_kmer_dir(c_fq, args.pos)
        # my_list = get_pair_with_barcode(mc_dir, c_dir, mc_fq, c_fq, args.dis)
        # with open(args.w_file, "w") as file:
        #     for i in my_list:
        #         print(i[0] + " " + i[1], file=file)







# if __name__ == '__main__':
#     mc=[]
#     c=[]
#     with open("mc_32000")as file,open("c_18000")as file2:
#         for line in file:
#             line=line.strip()
#             mc.append(line)
#         for line in file2:
#             line=line.strip()
#             c.append(line)
#     my_list=mc+c
#     my_list.append('32000')
#
#     aa=process_records(my_list,32000)
#     with open("raw1","w")as file,open("raw2","w")as file2:
#         for i in aa:
#             print(i[0],file=file)
#             print(i[1],file=file2)



