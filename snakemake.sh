#! /bin/bash
# This file is is used to submit the snakemake job
# Usage: sbatch --cpus-per-task=32 --mem=64g snakemake.sh

module load snakemake samtools peakranger macs gem music PePr bedtools sicer R || exit 1

sbcmd="sbatch --cpus-per-task={threads} --mem={cluster.mem}"
sbcmd+=" --time={cluster.time} --job-name={cluster.jobname} --partition={cluster.partition}"
sbcmd+=" --out={cluster.out} {cluster.extra}"

snakemake -pr --keep-going --rerun-incomplete --local-cores $SLURM_CPUS_PER_TASK \
    --jobs --cluster-config cluster.json --cluster "$sbcmd" \
    --latency-wait 180 all
