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
        filteredmergedquant=join(WORKDIR,"results","mergedquant.filtered.tsv")
    params:
        rowsumfilter=config['aggregate_quant_rowsum_filter']
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
head -n1 {output.mergedquant} > {output.filteredmergedquant}
tail -n +2 {output.mergedquant} |\
awk -F"\\t" -v r={params.rowsumfilter} -v OFS="\\t" '{{for(i=2;i<=NF;i++){{sum[NR]=sum[NR]+$i}};if (sum[NR] >= r) {{print sum[NR],$0}}}}' |\
sort -k1,1gr |\
cut -f2- >> {output.filteredmergedquant}
"""


    