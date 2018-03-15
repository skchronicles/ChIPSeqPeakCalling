##########################################################################################################
# Benchmarking of ChIP-Seq peak Callers
# Author: Skyler Kuhn (NIH/NCI) [C]
# CCR Collaborative Bioinformatics Resource
# Version 1.0.0
# See README.MD for more information
# USAGE:
#   sbatch --cpus-per-task=8 --mem=16g snakemake.sh
##########################################################################################################

import sys

shell.prefix("set -eo pipefail; ")

configfile: "config.yaml"
localrules: all

controls = config["controls"]
if controls is None:
    sys.exit("Controls are needed")

samples_narrow = config["samples_narrow"]
if samples_narrow is None:
    samples_narrow = []

samples_broad = config["samples_broad"]
if samples_broad is None:
    samples_broad = []

ALL_SAMPLES = samples_narrow + samples_broad + controls

ALL_BAM = expand("bam/{sample}.sorted.Q5DD.bam", sample = ALL_SAMPLES)
ALL_BAM.extend(expand("bam/{sample}.sorted.Q5DD.bam.bai", sample = ALL_SAMPLES))
ALL_BAM.extend(expand("bam/{sample}.sorted.Q5DD.bam.flagstat", sample = ALL_SAMPLES))
ALL_BAM.extend(("bam/control.sorted.Q5DD.bam", "bam/control.sorted.Q5DD.bam.bai", "bam/control.sorted.Q5DD.bam.flagstat"))


ALL_PEAKS = expand("peaks/mac2/narrow/{sample}_peaks.narrowPeak", sample = samples_narrow) + \
expand("peaks/mac2/broad/{sample}_peaks.broadPeak", sample = samples_broad)

rule all:
    input: ALL_PEAKS + ALL_BAM

rule merge_controls:
    input:   bam = expand("bam/{sample}.sorted.Q5DD.bam", sample = controls),
             bai = expand("bam/{sample}.sorted.Q5DD.bam.bai", sample = controls),
             flagstat = expand("bam/{sample}.sorted.Q5DD.bam.flagstat", sample = controls) 
    output:  bam = "bam/control.sorted.Q5DD.bam",
             bai = "bam/control.sorted.Q5DD.bam.bai",
             flagstat = "bam/control.sorted.Q5DD.bam.flagstat"
    log:     "log/merge_controls"
    threads: 2
    shell:
        '''
        inbam=( {input.bam} )
        if [[ ${{#inbam[@]}} -eq 1 ]]; then
            ln -s $(cd $(dirname {input.bam}) && pwd)/$(basename {input.bam}) {output.bam}
            ln -s $(cd $(dirname {input.bai}) && pwd)/$(basename {input.bai}) {output.bai}
            ln -s $(cd $(dirname {input.flagstat}) && pwd)/$(basename {input.flagstat}) {output.flagstat}
        else
            module load samtools/1.2
            samtools merge -r -@{threads} {output.bam} {input.bam}
        fi
        '''
	
rule MAC2_narrowPeaks:
    input:  "bam/{sample}.sorted.Q5DD.bam",
            "bam/control.sorted.Q5DD.bam"
    output: "peaks/mac2/narrow/{sample}_model.r", "peaks/mac2/narrow/{sample}_peaks.narrowPeak",
            "peaks/mac2/narrow/{sample}_peaks.xls", "peaks/mac2/narrow/{sample}_summits.bed",
            "peaks/mac2/narrow/{sample}_model.pdf"
    log:    "log/mac2/{sample}.find_narrow_peaks"
    threads: 2
    shell:
        '''
        module load macs/2.1.0.20150420 R
        macs2 callpeak -t {input[0]} \
        -c {input[1]} -f BAM -g {config[macs_g]} \
        --outdir peaks/mac2/narrow -n {wildcards.sample} -q 0.001 &> {log}
	cd peaks/mac2/narrow && Rscript {wildcards.sample}_model.r
        '''

rule MAC2_broadPeaks:
    input:  "bam/{sample}.sorted.Q5DD.bam", "bam/control.sorted.Q5DD.bam"
    output: "peaks/mac2/broad/{sample}_peaks.xls", "peaks/mac2/broad/{sample}_peaks.broadPeak"
    log:    "log/mac2/{sample}.find_broad_peaks"
    threads: 2
    shell:
        '''
        module load macs/2.1.0.20150420
        macs2 callpeak -t {input[0]} \
        -c {input[1]} -f BAM -g {config[macs_g]} \
        --broad --broad-cutoff 0.1 --nomodel --extsize 150 \
        --outdir peaks/mac2/broad -n {wildcards.sample} -q 0.001 &> {log}
        ''' 