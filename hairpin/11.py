configfile: "config.yaml"

from os.path import join
from os.path import exists

spikein=config["spikein"]



rule all: 
    input:
        expand("{trim_out_path}/{sample}.fq",
            sample=config["sampleList"],
            trim_out_path=config["trim_out_path"]),
        expand("{trim_out_path}/{sample}_cut_f1.fq",
            sample=config["sampleList"],
            trim_out_path=config["trim_out_path"]),
        expand("{trim_out_path}/{sample}_cut_f2.fq",
            sample=config["sampleList"],
            trim_out_path=config["trim_out_path"]),
        
        expand("data/bam/{sample}_sorted.bam",
            sample=config["sampleList"]),
        expand("data/bam/{sample}_sorted.bam.bai",
            sample=config["sampleList"]),

        expand("data/dedup/{sample}_dedup.bam",
            sample=config["sampleList"]),
        expand("data/hmc/{sample}.all.bed.hmc",
            sample=config["sampleList"]),
        expand("data/mc/{sample}.all.bed.mc",
            sample=config["sampleList"]),

        expand("spikein/hmc/{sample}.spikein.all.bed.hmc",
            sample=config["sampleList"]),
        expand("hairpin_results/hmc/{sample}.ref_nochrm.all.bed.hmc",
            sample=config["sampleList"]),

rule hairpin_cut: 
    input:
        fq1=expand("{fq_in_path}/{{sample}}_1.fq.gz",fq_in_path=config["fq_in_path"]),
        fq2=expand("{fq_in_path}/{{sample}}_2.fq.gz",fq_in_path=config["fq_in_path"])
    output:
        expand("{trim_out_path}/{{sample}}.fq",
            trim_out_path=config["trim_out_path"]),
        expand("{trim_out_path}/{{sample}}_cut_f1.fq",
            trim_out_path=config["trim_out_path"]),
        expand("{trim_out_path}/{{sample}}_cut_f2.fq",
            trim_out_path=config["trim_out_path"]),
                      
    log:
        "logs/hairpin_cut/{sample}.log"
    params:
        rule=config['rule'],
        out_file= expand("{trim_out_path}/{{sample}}",
            trim_out_path=config["trim_out_path"]),
    threads: 20
    shell:
        """
        (python3 ref/hairpin_cut.py --fq1 {input.fq1} --fq2 {input.fq2}  \
        --outfile {params.out_file} --rule {params.rule} --parallel {threads}) > {log} 2>&1
        """

rule bwa_map: 
    input:
        config["ref_genome"], # reference genome
        expand("{trim_out_path}/{{sample}}.fq",
            trim_out_path=config["trim_out_path"],)
    output:
        bam="data/bam/{sample}_sorted.bam",
        bai="data/bam/{sample}_sorted.bam.bai",
    log:
        "logs/bwa_mem/{sample}.log"
    threads: 20
    params:
        "-c 500 -U 17"   #bwa mem 默认参数,可省略
    shell:
        """
        (bwa mem -t {threads} {params} {input} | \
        samtools sort -@ {threads} -T data/bam/{wildcards.sample} - > {output.bam}) 2> {log}
        samtools index {output.bam}
        """

rule mark_duplicate: 
    input:
        "data/bam/{sample}_sorted.bam"
    output:
        "data/dedup/{sample}_dedup.bam"
    log:
        "logs/dedup/{sample}.log"
    threads: 10
    params:
        fq1=expand("{fq_in_path}/{{sample}}_1.fq.gz",fq_in_path=config["fq_in_path"]),
        fq2=expand("{fq_in_path}/{{sample}}_2.fq.gz",fq_in_path=config["fq_in_path"])
    shell:
        """
        (python3 ref/remove_dup.py --input {input} --output {output} --fq1 {params.fq1} --fq2 {params.fq2} ) 2> {log}
        samtools index {output}
        """

rule extract_hmc:
    input:
        "data/dedup/{sample}_dedup.bam"
    output:
        "data/hmc/{sample}.all.bed.hmc"
    log:
        "logs/extract_hmc/{sample}.log"
    threads: 20
    params:
        ref=config["ref_fa"]
    shell:
        """
        (python3 ref/extract_hmc.py --sam {input} --output {output} --ref {params.ref} --parallel {threads} ) 2> {log}
        """


rule extract_mc:
    input:
        bam="data/dedup/{sample}_dedup.bam",
        fq1=expand("{trim_out_path}/{{sample}}_cut_f1.fq",
            trim_out_path=config["trim_out_path"])
    output:
        "data/mc/{sample}.all.bed.mc"
    log:
        "logs/extract_mc/{sample}.log"
    threads: 20
    params:
        ref=config["ref_fa"]
    shell:
        """
        (python3 ref/extract_mc.py --sam {input.bam} --cutfq1 {input.fq1} --output {output} --ref {params.ref} --parallel {threads} ) 2> {log}
        """

rule call_spikein:
    input:
        hmc="data/hmc/{sample}.all.bed.hmc",
        mc="data/mc/{sample}.all.bed.mc"
    output:
        all_hmc="spikein/hmc/{sample}.spikein.all.bed.hmc",
        all_mc="spikein/mc/{sample}.spikein.all.bed.mc",
        cpg_hmc="spikein/hmc/{sample}.spikein.CpG.bed.hmc",
        cpg_mc="spikein/mc/{sample}.spikein.CpG.bed.mc"
    threads: 10
    params:
        awk_condition=" || ".join(f'$1 == "{value}"' for value in spikein),
        left="{",
        right="}", 
        other="print"
    shell:       
        """
        awk '{params.awk_condition} '  {input.hmc} > {output.all_hmc}
        awk '{params.awk_condition} && $7=="CpG" ' {input.hmc} > {output.cpg_hmc}
        awk '{params.awk_condition} ' {input.mc} > {output.all_mc}
        awk '{params.awk_condition} && $7=="CpG" ' {input.mc} > {output.cpg_mc}
        """ 



rule call_ref_no_chrm:
    input:
        hmc="data/hmc/{sample}.all.bed.hmc",
        mc="data/mc/{sample}.all.bed.mc"
    output:
        all_hmc="hairpin_results/hmc/{sample}.ref_nochrm.all.bed.hmc",
        all_mc="hairpin_results/mc/{sample}.ref_nochrm.all.bed.mc",
        cpg_hmc="hairpin_results/hmc/{sample}.ref_nochrm.CpG.bed.hmc",
        cpg_mc="hairpin_results/mc/{sample}.ref_nochrm.CpG.bed.mc"
    threads: 10
    params:
        awk_condition=" && ".join(f'$1 != "{value}"' for value in spikein),
    shell:
        """
        awk '{params.awk_condition} && $1 != "chrM" ' {input.hmc} > {output.all_hmc}
        awk '{params.awk_condition} && $1 != "chrM" && $7 == "CpG" ' {input.hmc} > {output.cpg_hmc}
        awk '{params.awk_condition} && $1 != "chrM" ' {input.mc} > {output.all_mc}
        awk '{params.awk_condition} && $1 != "chrM" && $7 == "CpG" ' {input.mc} > {output.cpg_mc}
        """

