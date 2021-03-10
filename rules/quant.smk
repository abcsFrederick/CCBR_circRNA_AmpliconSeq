BSJFA=rules.create_bsjfa.output.bsjfa

rule salmon:
    input:
        bam=rules.align.output.bam,
        reffa=BSJFA,
    output:
        sf=join(WORKDIR,"results","{sample}","salmon","quant.sf")
    params:
        salmon_parameters=config["salmon_parameters"]
    envmodules: TOOLS["salmon"]["version"]
    shell:"""
outdir=$(dirname {output.sf})
salmon quant \
--threads {threads} \
{params.salmon_parameters} \
-t {input.reffa} \
-a {input.bam} \
--output $outdir
"""

rule aggregate_quant:
    input:
        expand(join(WORKDIR,"results","{sample}","salmon","quant.sf"),sample=SAMPLES)
    output:
        mergedquant=join(WORKDIR,"results","mergedquant.tsv"),
        # filteredmergedquant=join(WORKDIR,"results","mergedquant.filtered.tsv")
    envmodules: TOOLS["salmon"]["version"]
    shell:"""
names=""
quants=""
for i in {input};do
    outdir=$(dirname $i)
    quants="$quants $outdir"
    samplename=$(echo $outdir|awk -F"/" '{{print $(NF-1)}}')
    names="$names $samplename"
done
salmon quantmerge --quants $quants --names $names --column numreads -o {output.mergedquant}
"""


    