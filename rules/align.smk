if config['bsjfa']:
    BSJFA=config['bsjfa']
else:
    BSJFA=config['bsjfadefault']
check_readaccess(BSJFA)

rule build_index:
    input:
        reffa=BSJFA,
    output:
        join(WORKDIR,"bowtie_index","bsj.1.bt2")
    params:
        outdir=join(WORKDIR,"bowtie_index")
    envmodules: TOOLS["bowtie2"]["version"]
    threads: 2
    shell: """
cp {input.reffa} {params.outdir}
cd {params.outdir}
reffa=$(basename {input.reffa})
bowtie2-build $reffa bsj
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
