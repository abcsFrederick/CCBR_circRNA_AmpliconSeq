# CCBR circRNA AmpliconSeq Snakemake Workflow
This is a snakemake workflow to process circRNA AmpliconSeq datasets generated to study circRNAs in KSHV and human hosts using divergent primers. Some basic usage instructions are as follows:

```
USAGE:
  bash ./run_circRNA_AmpliconSeq -m/--runmode=<RUNMODE> -w/--workdir=<WORKDIR>
Required Arguments:
1.  RUNMODE: [Type: String] Valid options:
    *) init : initialize workdir
    *) run : run with slurm
    *) reset : DELETE workdir dir and re-init it
    *) dryrun : dry run snakemake to generate DAG
    *) unlock : unlock workdir if locked by snakemake
    *) runlocal : run without submitting to sbatch
2.  WORKDIR: [Type: String]: Absolute or relative path to the output folder with write permissions.
```