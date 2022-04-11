# 16/11/2015, 
# SASSI EDITS
# install.packages("seqinr")
library(seqinr)
genome=read.fasta(file="/Users/Sassi/Documents/Postdoc/Miseq/Genomic_Analysis/USA300/SNPs_genomes/Genome1.fasta")

pos_tabl=read.table("/Users/Sassi/Documents/Postdoc/Miseq/Genomic_Analysis/USA300/SNPs_genomes/tab.txt", header=T)
pos_tabl = as.data.frame(pos_tabl)

all_genomes=list()
for(i in 2:ncol(pos_tabl)){
  mut_genome = genome$Genome1
  mut_genome[pos_tabl$Position] = paste(pos_tabl[,i])
  all_genomes[[i-1]] = mut_genome
}
names(all_genomes) = colnames(pos_tabl[,2:ncol(pos_tabl)])
summary(all_genomes) # to check the content
# loop to write all genomes in fasta format
# default, files are saved in current directory, getwd()
for(i in 1:length(all_genomes)){
  write.fasta(all_genomes[[i]], names=names(all_genomes)[i], 
            file.out=paste(names(all_genomes)[i],".fasta", sep=""))
}

############################################################################
# Function for more general use
mut_genomes = function(genome_fasta,pos_table){
  all_genomes=list()
  for(i in 2:ncol(pos_table)){
    mut_genome = genome$Genome1
    mut_genome[pos_table$Position] = paste(pos_table[,i])
    all_genomes[[i-1]] = mut_genome
  return(all_genomes)  
  }
#############################################################################"


