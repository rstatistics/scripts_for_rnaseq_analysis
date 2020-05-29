# scripts_for_rnaseq_analysis
This repo contains scripts to analyze rnaseq data

## Quantification with salmon
### install salmon
```
conda install -c bioconda salmon=1.2.1
```
### fetch transcripts with gffread
```
gffread gene.gtf -g genome.fa -w transcript.fa
```
### build salmon index
```
salmon index -t transcript.fa -i salmon-index
```
### quantification
```
ls *_1.fastq.gz | sed 's/_1.fastq.gz//' | \
        xargs -i salmon quant -i salmon-index -l A -1 {}_1.fastq.gz -2 {}_2.fastq.gz -o {}
```
