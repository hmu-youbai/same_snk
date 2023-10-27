configfile: "config.yaml"

from os.path import join
from os.path import exists

rule all: 
    input:
        expand("data/trand_bed/{sample}.trand",
            sample=config["sampleList"]),
        expand("data/intersect/{sample}_{region_bed}",
            region_bed=config["region"],
            sample=config["sampleList"]),
        expand("data/merged/{sample}_merged.bed",
            sample=config["sampleList"]),







rule verify_ref_and_trands:
    input:
        "data/bed/{sample}"
    output:
        "data/trand_bed/{sample}.trand"

    params:
        input_chr=config["input_chr"],
        input_start=config["input_start"],
        input_ref=config["ref_fa"]
    shell:
        """
        python3 ref/verify_trand.py {input} {params.input_ref} {params.input_chr} {params.input_start} {output}
        """

rule intersect: 
    input:
        input_bed="data/trand_bed/{sample}.trand",
        input_region="region/{region_bed}",

    output:
        "data/intersect/{sample}_{region_bed}",
            
    shell:
        """
        bedtools intersect -a {input.input_bed} -b {input.input_region}  -wo  > {output}

        """

rule merge_bed:
    input:
        lambda wildcards: expand("data/intersect/{sample}_{region_bed}",
                      sample=wildcards.sample,
                      region_bed=config["region"])

    output:
        "data/merged/{sample}_merged.bed"

    shell:
        """
        cat {input} | awk '{{if ($(NF-1) == $(NF-7)) print $0}}' > {output}
        """

