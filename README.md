# CCBR circRNA AmpliconSeq Snakemake Workflow

Process circRNA AmpliconSeq datasets generated to study circRNAs in KSHV and human hosts using divergent primers

![GitHub issues](https://img.shields.io/github/issues/CCBR/CCBR_circRNA_AmpliconSeq)![forks](https://img.shields.io/github/forks/CCBR/CCBR_circRNA_AmpliconSeq)![stars](https://img.shields.io/github/stars/CCBR/CCBR_circRNA_AmpliconSeq)![LICENSE](https://img.shields.io/github/license/CCBR/CCBR_circRNA_AmpliconSeq)
[![release](https://img.shields.io/github/v/release/CCBR/CCBR_circRNA_AmpliconSeq?color=blue&label=latest%20release)](https://github.com/CCBR/CCBR_circRNA_AmpliconSeq/releases/latest)

## Usage

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
