# Snakefile

# 定义要处理的所有样本文件列表
samples = [
    "C11_1.fastq.gz", "C11_2.fastq.gz", "C12_1.fastq.gz", "C12_2.fastq.gz",
    "C15_1.fastq.gz", "C15_2.fastq.gz", "C18_1.fastq.gz", "C18_2.fastq.gz",
    "C20_1.fastq.gz", "C20_2.fastq.gz", "C22_1.fastq", "C22_2.fastq",
    "C25_1.fastq", "C25_2.fastq", "C27_1.fastq.gz", "C27_2.fastq.gz",
    "C28_1.fastq.gz", "C28_2.fastq.gz", "C29_1.fastq.gz", "C29_2.fastq.gz",
    "C35_1.fastq.gz", "C35_2.fastq.gz", "C43_1.fastq.gz", "C43_2.fastq.gz",
    "C48_1.fastq.gz", "C48_2.fastq.gz", "C5_1.fastq", "C5_2.fastq"
]

# 定义文件路径变量
data_path = "/home/milab"
qc_path= "/home/milab/qc"
trimmed_path = "/home/milab/trimmed"

# 定义fastqc程序的路径


# 定义fastqc规则
rule fastqc:
    input:
        expand(data_path+"/{sample}",sample=samples )
    output:
        qc_html=expand(qc_path + "/{sample}_fastqc.html", sample=samples),
        qc_zip=expand(qc_path + "/{sample}_fastqc.zip", sample=samples)
    params:
        extra_options=""  # 可以添加其他需要的fastqc选项
    shell:
        """
        fastqc {input} -o {qc_path} {params.extra_options}
        """

# 定义all规则来运行所有的fastqc规则
rule all:
    input:
        expand(qc_path + "/{sample}_fastqc.html", sample=samples),
        expand(qc_path + "/{sample}_fastqc.zip", sample=samples)

# 如果你有其他的规则，可以在这里继续定义
# ...

# 如果有需要预处理样本的步骤（例如trimming等），你可以在这里定义更多规则
# ...

