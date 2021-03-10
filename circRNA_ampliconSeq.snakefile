from os.path import join

include: join("rules/init.smk")
include: join("rules/trim.smk")
include: join("rules/align.smk")
include: join("rules/quant.smk")


localrules: all

rule all:
    input:
        #trim
        expand(join(WORKDIR,"results","{sample}","trim","{sample}.R1.trim.fastq.gz"),sample=SAMPLES),
        expand(join(WORKDIR,"results","{sample}","trim","{sample}.R2.trim.fastq.gz"),sample=SAMPLES),
        #align
        expand(join(WORKDIR,"results","{sample}","bam","{sample}.bam"),sample=SAMPLES),
        #quant
        mergedquant=join(WORKDIR,"results","mergedquant.tsv")
        


