mafft --thread 8 ps_16S_ASVs.fasta > ps_16S_ASVs_aln.fasta
trimal -in ps_16S_ASVs_aln.fasta -out ps_16S_ASVs_aln_nogaps.fasta -gappyout
fasttree -nt -gtr -gamma ps_16S_ASVs_aln_nogaps.fasta > ps_16S_ASVs.tre
