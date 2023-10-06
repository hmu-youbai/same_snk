import pandas as pd
# Snakemake RNA-seq pipeline

# Define the input files

SAMPLES = [ "A104", "A56","A60","A65","A67","A70","A71","A77","A79","A80",
            "A83", "A84", "A86","A84-techrep", "A90","A91","A96","A99"]

# Define the rule to run fastp on each sample
rule all:
    input:
        "output/expression_matrix.txt"


rule run_fastp:
    input:
        r1="data/{sample}_1.fastq.gz",
        r2="data/{sample}_2.fastq.gz"
    output:
        r1_trimmed="output/{sample}_1.trimmed.fastq.gz",
        r2_trimmed="output/{sample}_2.trimmed.fastq.gz"
    params:
        report="output/{sample}_fastp_report.html",
        json="output/{sample}_fastp_report.json"
    threads: 4
    shell:
        """
        fastp --in1 {input.r1} --in2 {input.r2} --out1 {output.r1_trimmed} --out2 {output.r2_trimmed} \
              --html {params.report} --json {params.json} --thread {threads}
        """

# Define the rule to run STAR alignment on each sample
rule run_star:
    input:
        r1="output/{sample}_1.trimmed.fastq.gz",
        r2="output/{sample}_2.trimmed.fastq.gz"
    output:
        bam="output/{sample}.star.Aligned.sortedByCoord.out.bam"
    params:
        genome_index="/home/index/hg38_star"  
    threads: 8
    shell:
        """
        STAR --genomeDir {params.genome_index} --readFilesIn {input.r1} {input.r2} \
             --outFileNamePrefix output/{wildcards.sample}.star. --runThreadN {threads} \
             --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 10000000000 --readFilesCommand zcat
        """



rule generate_expression_matrix:
    input:
        bam="output/{sample}.star.Aligned.sortedByCoord.out.bam"
    output:
        expression_matrix="output/{sample}_expression_matrix.txt" 
    params:
        gtf_file="/home/index/gencode.v44.chr_patch_hapl_scaff.annotation.gtf"  
    shell:
        """
        featureCounts -T 8 -p -t exon -g gene_id -a {params.gtf_file} -o {output.expression_matrix} {input.bam}
        """



rule link_expression_matrix:
    input:
        expand("output/{sample}_expression_matrix.txt", sample=SAMPLES)
    output:
        "output/expression_matrix.txt"
    run:
        with open(output[0], "w") as out, open(input[0]) as f:
            out.write(f.read())
