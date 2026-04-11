# Haplotype 0
gtime -a ../target/release/pgr align alnmap \
    --reference-fasta-path GCA_000001405.15_GRCh38_no_alt_analysis_set.PanSN.fa.gz \
    --assembly-contig-path hg002v1.1.mat_MT.PanSN.fa.gz \
    --output-prefix hg002_hap0 \
    --preset default

# Haplotype 1
gtime -a ../target/release/pgr align alnmap \
    --reference-fasta-path GCA_000001405.15_GRCh38_no_alt_analysis_set.PanSN.fa.gz \
    --assembly-contig-path hg002v1.1.pat.PanSN.fa.gz \
    --output-prefix hg002_hap1 \
    --preset default
