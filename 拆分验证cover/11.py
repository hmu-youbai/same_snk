import os

# 定义输入BAM文件和采样比例列表
input_bam = "p0.bam"
sampling_ratios = [i/10 for i in range(1,10, 1)]  

# 定义文件路径变量
out_dir = "/home/wangzc/bam_cgs/out_cgs"

# 定义中间变量，存储qualimap的输出文件路径
qualimap_outputs = [out_dir + f"/qualimap_{sample}" for sample in sampling_ratios]

# 定义随机采样、排序和索引规则
rule all:
    input: ["/home/wangzc/bam_cgs" + f"/{sample}/genome_results.txt" for sample in sampling_ratios]

rule random_sampling:
    output: out_dir + "/sampled_{sample}.bam"
    params:
        sample_ratio="{sample}"
    shell:
        """
        samtools view -bs {params.sample_ratio} {input_bam} > {output}
        """

rule sort_and_index:
    input: out_dir + "/sampled_{sample}.bam"
    output: out_dir + "/sorted_{sample}.bam", out_dir + "/sorted_{sample}.bam.bai"
    shell:
        """
        samtools sort {input} -o {output[0]}
        samtools index {output[0]}
        """

# 定义Qualimap规则
rule qualimap_bamqc:
    input: out_dir + "/sorted_{sample}.bam"
    output: "/home/wangzc/bam_cgs" + "/{sample}/genome_results.txt"
    params:
        threads=12,
        java_mem_size="20G"
    shell:
        """
        qualimap bamqc -bam {input} -outdir {wildcards.sample} -outformat PDF:HTML -nt {params.threads} --java-mem-size={params.java_mem_size}
        """
