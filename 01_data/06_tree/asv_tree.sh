mafft --thread 8 ps_all_ASVs.fasta > ps_all_ASVs_aln.fasta
trimal -in ps_all_ASVs_aln.fasta -out ps_all_ASVs_aln_nogaps.fasta -gappyout
fasttree -nt -gtr -gamma ps_all_ASVs_aln_nogaps.fasta > ps_all_ASVs.tre
