import pandas as pd


configfile: "config.yaml"




# Define the rule to run fastp on each sample
rule all:
    input:
        expand("data/bam/{sample}.star.Aligned.sortedByCoord.out.bam",
            sample=config["SAMPLES"]),
        expand("data/count/{sample}_expression_matrix.txt",
            sample=config["SAMPLES"]),

rule run_fastp:
    input:
        r1="data/fastq/{sample}_1.fastq.gz",
        r2="data/fastq/{sample}_2.fastq.gz"
    output:
        r1_trimmed="data/trim/{sample}_1.trimmed.fastq.gz",
        r2_trimmed="data/trim/{sample}_2.trimmed.fastq.gz"
    params:
        report="data/trim/html/{sample}_fastp_report.html",
        json="data/trim/json/{sample}_fastp_report.json"
    conda:
        'py39'
    threads: 4
    log:
        "logs/trim/{sample}.log"
    shell:
        """
        (fastp --in1 {input.r1} --in2 {input.r2} --out1 {output.r1_trimmed} --out2 {output.r2_trimmed} \
              --html {params.report} --json {params.json} --thread {threads}) 2> {log}
        """

# Define the rule to run STAR alignment on each sample
rule run_star:
    input:
        r1="data/trim/{sample}_1.trimmed.fastq.gz",
        r2="data/trim/{sample}_2.trimmed.fastq.gz"
    output:
        bam="data/bam/{sample}.star.Aligned.sortedByCoord.out.bam"
    params:
        genome_index=config["star_index"]
    conda:
        'py39'
    threads: 8
    log:
        "logs/star/{sample}.log"
    shell:
        """
        (STAR --genomeDir {params.genome_index} --readFilesIn {input.r1} {input.r2} \
             --outFileNamePrefix data/bam/{wildcards.sample}.star. --runThreadN {threads} \
             --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 10000000000 --readFilesCommand zcat )  2> {log} 
        """



rule generate_expression_matrix:
    input:
        bam="data/bam/{sample}.star.Aligned.sortedByCoord.out.bam"
    output:
        expression_matrix="data/count/{sample}_expression_matrix.txt" 
    params:
        gtf_file=config["star_gtf"]
    conda:
        'py39'
    log:
        "log/featureCount/{sample}.log"
    shell:
        """
        (featureCounts -T 8 -p -t exon -g gene_id -a {params.gtf_file} -o {output.expression_matrix} {input.bam}) 2> {log}
        """


