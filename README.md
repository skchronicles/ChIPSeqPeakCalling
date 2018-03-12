# ChIP-Seq Peakcalling Benchmarking
   

## 1. PeakRanger (BCP)   
##### Description:  
PeakRanger is a multi-purporse software suite for analyzing next-generation sequencing (NGS) data. 
It contains the BCP tool for broad peak calling. BCP supports HTML-based annotation reports. 
BCP is installed on [Biowulf.](https://hpc.nih.gov/apps/peakranger.html) 

##### Loading PeakRanger on Biowulf:  

    module load peakranger
 
##### Running BCP:  

    peakranger bcp \
    --data {expt1.bam} \
    --control {control.bam} \
    --format bam \
    --output bcp_results
    
*We can also do nr, lc, ranger. ccat may be useful for broad peaks. nr and lc can get use some QC stats. Also check /data/CCBR_Pipeliner/db/PipeDB/bin/BCP_v1.1 ... refer /data/CCBR_Pipeliner/db/PipeDB/bin/BCP_v1.1/manual.pdf to run BCP_TF. This may be an older version of the BCP tool, which we should look into it and rule it out!*


## 2. MAC2 
##### Description:
MACS empirically models the length of the sequenced ChIP fragments, which tends to be shorter than sonication or library 
construction size estimates, and uses it to improve the spatial resolution of predicted binding sites. MACS also uses a 
dynamic Poisson distribution to effectively capture local biases in the genome sequence, allowing for more sensitive and 
robust prediction. MACS compares favorably to existing ChIP-Seq peak-finding algorithms and can be used for ChIP-Seq with 
or without control samples.MAC2 is installed on [Biowulf.](https://hpc.nih.gov/apps/macs.html)  

##### Loading MACs on Biowulf:

    module load macs


##### Running MAC2:  

    macs2 callpeak \
    -t {input1.bam} \
    -c {ctrl.bam} \
    --call-summits \
    -f BAM \
    -g {params.gsize} \
    -n {wildcards.name} \
    --outdir {macsn_dir}/{wildcards.group} \
    -B -q 0.01

*need to use the --nomodel --extsize $extsize. get extsize by running phantompeakqualtools ...something like this--> Rscript /data/CCBR_Pipeliner/3.0/Pipeliner/Results-template/Scripts/phantompeakqualtools/run_spp_nodups.R -c=$bamfile -p=${NTHREADS} -savp=${CC_PLOT} -out=${CC_SCORES}... the CC_Scores files will have the peaks, pick the first one for extsize*


## 3. SISSRS
##### Description: 
SISSRs is a software application for precise identification of genome-wide transcription factor binding
sites from ChIP-Seq data. SISSCRS is currently not installed Biowulf, but more information-- including installation details-- 
can be found on it's [homepage.](https://dir.nhlbi.nih.gov/papers/lmi/epigenomes/sissrs/SISSRs-Manual.pdf)

##### Running SISSRS:  

     sissrs.pl \
     -i {input-file.bed} \
     -o {output-file} \
     -s {genome-size}
   * {genome-size} is the effective genome size (or length): e.g. 3080436051 for hg18
  
*effective or mappable genome size is actually the genome size minus the Ns.. you can write a script to get this for hg18,hg38,mm9,mm10. This project seems abandoned to me... may not be worthwhile digging too much into SISSRs*

## 4. GEM 
##### Description: 
GEM is a high-resolution peak calling and motif discovery tool for ChIP-seq and ChIP-exo data. 
GEM is installed on [Biowulf.](https://hpc.nih.gov/apps/gem.html)

##### Loading GEM on Biowulf:

    module load gem

##### Running GEM:  

    java -Xmx10g -jar $GEMJAR --t 24 \
    --d ./Read_Distribution_default.txt \
    --g ./mm10.chrom.sizes 
    --genome /fdb/igenomes/Mus_musculus/UCSC/mm10/Sequence/Chromosomes/ \
    --s 2000000000 
    --expt SRX000540_mES_CTCF.bed \
    --ctrl SRX000543_mES_GFP.bed \
    --f BED \
    --out mouseCTCF --k_min 6 --k_max 13
    
* looks good ... does it take bams? *


## 5. MUSIC
##### Description: 
MUSIC is a tool for identification of enriched regions at multiple scales in the read depth signals from ChIP-Seq experiments. 
MUSIC is installed on [Biowulf.](https://hpc.nih.gov/apps/music.html)

##### Loading Music on Biowulf:

    module load samtools  # needed to convert to sam format
    module load music

##### Running Music:
    mkdir chip; mkdir input
    samtools view chip.bam | MUSIC -preprocess SAM stdin chip/ 
    samtools view input.bam | MUSIC -preprocess SAM stdin input/
    samtools view /directory/to/chip.bam | MUSIC -preprocess SAM stdin chip/ 
    samtools view /directory/to/input.bam | MUSIC -preprocess SAM stdin input/
    mkdir chip/sorted;mkdir chip/dedup;mkdir input/sorted;mkdir input/dedup
    MUSIC -sort_reads chip chip/sorted 
    MUSIC -sort_reads input input/sorted 
    MUSIC -remove_duplicates chip/sorted 2 chip/dedup 
    MUSIC -remove_duplicates input/sorted 2 input/dedup
      
    MUSIC -get_multiscale_broad_ERs \
    -chip chip/dedup \
    -control input/dedup \
    -mapp Mappability_36bp \
    -l_mapp 36 \
    -begin_l 1000 \
    -end_l 16000 \
    -step 1.5
   * This code tells MUSIC to identify the enriched regions starting from 1kb smoothing window length upto 16kb with 
   multiplicative factor of 1.5 using the default parameters for the remaining parameters. The ERs for each scale are dumped.
   
* You will have to generate appropriate mappability files (https://github.com/gersteinlab/MUSIC#multi-mappability-profile-generation) Our read lengths are 100 or 125, so generate for those. Also, make sure you have bowtie2 indices for the genomes. Look at https://github.com/gersteinlab/MUSIC#running-music-with-default-parameters-and-automatic-selection-of-l_p-parameter- to automate parameter selection.*
    

## 6. PePr
##### Description:
PePr is a ChIP-Seq Peak-calling and Prioritization pipeline that uses a sliding window approach and models read counts 
across replicates and between groups with a negative binomial distribution. PePr empirically estimates the optimal 
shift/fragment size and sliding window width, and estimates dispersion from the local genomic area. Regions with less 
variability across replicates are ranked more favorably than regions with greater variability. Optional post-processing 
steps are also made available to filter out peaks not exhibiting the expected shift size and/or to narrow the width of peaks. 
PePr is installed on [Biowulf.](https://hpc.nih.gov/apps/PePr.html)

##### Loading PePr on Biowulf:  
    module load PePr

##### Running Pepr:
    PePr -c chip_rep1.bam,chip_rep2.bam \
    -i input_rep1.bam,input_rep2.bam \
    -f bam \
    -n {expname}
    
 * --shiftsize Half the fragment size.. again comes for ppqt... this should be half ext size that we used for macs
 -f needs to be bampe for PE data...not to worry about this now*
