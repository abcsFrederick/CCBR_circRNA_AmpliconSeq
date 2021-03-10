# if config['bsjfa']:
#     BSJFA=config['bsjfa']
# else:
#     BSJFA=config['bsjfadefault']
# check_readaccess(BSJFA)
PRIMERINFO=config['primertsv']
check_readaccess(PRIMERINFO)
REFFA=config['reffa']
check_readaccess(REFFA)

rule create_primerbed:
    input:
        primertsv=PRIMERINFO
    output:
        primerbed=join(WORKDIR,"bowtie_index","primers.bed")
    params:
        script=join(SCRIPTSDIR,"build_primer_bed.py")
    envmodules: TOOLS["python"]["version"]
    shell:"""
python {script} {input.primertsv} > {output.primerbed}
"""

rule create_bsjfa:
    input:
        bed=rules.create_primerbed.output.primerbed,
        reffa=REFFA
    output:
        bsjfa=join(WORKDIR,"bowtie_index","bsj.fa")
    params:
        script=join(SCRIPTSDIR,"generate_bsj_fasta_from_primer_bed.py")
    container: "docker://nciccbr/ccbr_htseq_0.13.5:latest"
    shell:"""
python {params.script} \
--bed {input.bed} \
--reffa {input.reffa} \
--outfa {output.bsjfa}
"""

rule build_index:
    input:
        bsjfa=rules.create_bsjfa.output.bsjfa,
    output:
        join(WORKDIR,"bowtie_index","bsj.1.bt2")
    params:
        outdir=join(WORKDIR,"bowtie_index")
    envmodules: TOOLS["bowtie2"]["version"]
    threads: 2
    shell: """
# cp {input.bsjfa} {params.outdir}
cd {params.outdir}
bsjfa=$(basename {input.bsjfa})
bowtie2-build $bsjfa bsj
gzip -n {input.bsjfa}
"""

rule align:
    input:
        r1=rules.cutadapt.output.of1,
        r2=rules.cutadapt.output.of2,
        bt2_index=rules.build_index.output,
    output:
        bam=join(WORKDIR,"results","{sample}","bam","{sample}.bam")
    params:
        sample="{sample}",
        workdir=WORKDIR,
        outdir=join(WORKDIR,"results","{sample}","bowtie2"),
        peorse=get_peorse,
        bowtie2_alignment_parameters=config["bowtie2_alignment_parameters"],
        mem=MEMORY
    threads: 56
    envmodules: TOOLS["bowtie2"]["version"], TOOLS["sambamba"]["version"]
    shell:"""
bt2_index=$(echo "{input.bt2_index}"|awk -F".1.bt2" '{{print $1}}')
if [ "{params.peorse}" == "PE" ];then
    bowtie2 \
    --threads {threads} \
    {params.bowtie2_alignment_parameters} \
    -x $bt2_index \
    -U {input.r1},{input.r2}
else
    bowtie2 \
    --threads {threads} \
    {params.bowtie2_alignment_parameters} \
    -x $bt2_index \
    -U {input.r1}
fi \
| awk -F"\t" '{{if ($6 ~ /60M/ || $1 ~ /^@/){{print}}}}' \
| sambamba view --nthreads={threads} -S --format=bam -o=/dev/stdout /dev/stdin \
| sambamba sort --memory-limit={params.mem}G --tmpdir=/dev/shm --nthreads={threads} --out={output.bam} /dev/stdin
"""
