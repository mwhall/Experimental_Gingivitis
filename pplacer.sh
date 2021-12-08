#!/bin/bash

cd data/pplacer_files
#Thanks for the walkthrough! https://www.polarmicrobes.org/phylogenetic-placement-re-re-visited/

## Degap
seqmagick mogrify --ungap HOMD_16S_rRNA_RefSeq_V15.1.aligned.fasta
## Align
cmalign --dna -o homd.sto --outformat Pfam 16S_bacteria.cm HOMD_16S_rRNA_RefSeq_V15.1.aligned.fasta
## Convert to fasta format
seqmagick convert homd.sto HOMD_16S_rRNA_RefSeq_V15.1.aligned.fasta
seqmagick mogrify --deduplicate-sequences HOMD_16S_rRNA_RefSeq_V15.1.aligned.fasta
raxmlHPC-SSE3 -T 8 -m GTRGAMMA -s HOMD_16S_rRNA_RefSeq_V15.1.aligned.fasta -n homd.tre -f d -p 12345
raxmlHPC-SSE3 -T 2 -m GTRGAMMA -f I -t RAxML_bestTree.homd.tre -n root.homd.tre
raxmlHPC-SSE3 -T 8 -m GTRGAMMA -f J -p 12345 -t RAxML_rootedTree.root.homd.tre -n conf.root.homd.tre -s HOMD_16S_rRNA_RefSeq_V15.1.aligned.fasta
taxit create -l 16S_rRNA -P homd.refpkg --aln-fasta HOMD_16S_rRNA_RefSeq_V15.1.aligned.fasta --tree-stats RAxML_info.homd.tre --tree-file RAxML_fastTreeSH_Support.conf.root.homd.tre

#Get the significant queries from the GLM results and the rep seq file
for tax in `awk 'BEGIN{FS="\t"}{if (($12 <= 0.05) || ($13 <= 0.05) || ($14 <= 0.05) || ($15 <=0.05)) {print;}}' ../GLM_results.tsv | cut -f 1 | tail -n +2`; do grep -A 1 ${tax} ../representative_sequences.fasta; done > queries_from_GLM.fasta

cat HOMD_16S_rRNA_RefSeq_V15.1.aligned.fasta queries_from_GLM.fasta > query.homd.fasta
seqmagick mogrify --ungap query.homd.fasta
cmalign --dna -o query.homd.sto --outformat Pfam 16S_bacteria.cm query.homd.fasta
seqmagick convert query.homd.sto query.homd.fasta
pplacer -o query.homd.jplace -p --keep-at-most 20 -c homd.refpkg query.homd.fasta
guppy to_csv --point-mass --pp -o query.homd.csv query.homd.jplace
guppy tog --node-numbers --xml --pp -o query.homd.tog.phyloxml query.homd.jplace
