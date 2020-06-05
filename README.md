# scripts_for_rnaseq_analysis
This repo contains scripts to analyze rnaseq data

### remove large files with rsync
Asume targets is our directory to delete, make a empty directory and run the following command:
```bash
rsync --delete-before --force -d empty/ targets/
```
### vlookup function in unix
Table_1.txt
```
1GR_P1:001PI
:040VG_L1
:001PO_L3
1JPI_P1:001PO_L1
1JPI_P1:001PO_L2
```
Table_2.txt
```
1JPI_P1:001PO_L1    1401UC
1JPI_P1:001PO_L2    1401UC
1HIK_P2:001ER       1402UC
1GR_P1:001PI        1402UC
```
You can use awk to realize vlookup function
vlookup.awk
```bash
FNR==NR{
  a[$1]=$2
  next
}
{ if ($1 in a) {print $1"\t"a[$1]} else {print $1"\t""NA"} }
```
Then use vlookup.awk function
```bash
awk -f vlookup.awk Table2.txt Table1.txt
```
### one line for vlookup function
```bash
awk 'NR==FNR{a[$1]=$2; next}; {print $1,$1 in a ? a[$1] : "NA"}' Table_2.txt Table_1.txt
```
or
```bash
awk 'NR==FNR{a[$1]=$2; next}; { if($1 in a){print $1"\t"a[$1]} else {print $1"\t""NA"} }' Table_2.txt Table_1.txt
```



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
### bedToGene
```bash
bedtools intersect -wa -wb -a bed -b gene.gtf|perl -ne 'if(/gene_id "(\S+?)"/){print $1,"\n"}'|sort -u > GeneList
```
