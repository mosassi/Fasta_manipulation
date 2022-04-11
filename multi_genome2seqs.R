##################################################################################
# EDITS 17/11/2015
##################################################################################
install.packages("seqinr")
library(seqinr)
coor_seqs=read.table("extraction.txt", header=F)
colnames(coor_seqs)=c("name", "start", "end")

all_fasta_files = list.files(path = ".", pattern="fasta") # path = path to fasta directory
# tu peux remplacer "fasta" par "USA300", ca va prendre tout les fichiers contenant le nom USA300 (ex. USA300_1.fasta, USA300_2.fasta...)
all_fasta = list()
# read all fasta 
for(i in 1:length(all_fasta_files)){ all_fasta[i] = read.fasta(file=paste(all_fasta_files[i]))}
names(all_fasta) = all_fasta_files

#je declare une liste qui va contenir les seqs a extraire 

mySeqs = list()

# exctract and save 
for(j in 1:length(all_fasta)){
    for(i in 1:nrow(coor_seqs)){
      #mySeqs[[i]]=getFrag(genome[[1]], coor_seqs$start[i], coor_seqs$end[i], name=coor_seqs$name[i])
      mySeqs[[i]]=getFrag(all_fasta[[j]], coor_seqs$start[i], coor_seqs$end[i], name=coor_seqs$name[i])
      setTxtProgressBar(txtProgressBar(min=0, max=nrow(coor_seqs), style=3), i)
    }
    #write.fasta(mySeqs,names=names(all_fasta[j]), file.out=paste(names(all_fasta[j])))
    write.fasta(mySeqs, names=coor_seqs$name, file.out=paste("anzaWeh_",names(all_fasta[j]), sep=""))
}
# attention, tu utilise les memes noms de sotie que ceux d'entré, ceux la vont étres ecrasés, c'est déconseillé!!
# don't run : file.out=«  le nom du fasta pris.fasta



