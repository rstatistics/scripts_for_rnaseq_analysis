extract_splice_sites.py gene.gtf > gene.ss
extract_exons.py gene.gtf > gene.exon
hisat2-build --ss gene.ss --exon gene.exon genome.fa star-index/genome
