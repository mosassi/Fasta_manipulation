install.packages("seqinr")
library(seqinr)
genome=read.fasta(file="/Users/Sassi/Documents/Postdoc/Miseq/Genomic_Analysis/USA300/sRNAs/X.fasta")
coor_seqs=read.table("Users/Sassi/Documents/Postdoc/Miseq/Genomic_Analysis/USA300/sRNAs/sRNA-pergenome/extraction.txt", header=F)
#je donne des noms aux 3 colonnes qu'on trouve dans ton fichier test_coor.xls
colnames(coor_seqs)=c("name", "start", "end")
getFrag(genome[[1]], 1, 10)

#je declare une liste qui va contenir les seqs a extraire 
mySeqs=vector("list", length=nrow(coor_seqs))
mySeqs = list(length=nrow(coor_seqs))


#j'utilise la fonction getFrag pour recuperer les fragments selon les coordonnées start & end que j'ai importé dans "coor_seqs"

for(i in 1:nrow(coor_seqs)){
    mySeqs[[i]]=getFrag(genome[[1]], coor_seqs$start[i], coor_seqs$end[i], name=coor_seqs$name[i])
    paste(mySeqs)
    #la commande suivante pour voir l'avancement de la boucle avec le %
    setTxtProgressBar(txtProgressBar(min=0, max=nrow(coor_seqs), style=3), i)
}

#je crée un fichier qui va contenir toutes les sequences extraites du génome
write.fasta(mySeqs,names=coor_seqs$name, file.out="sRNA_x.fasta")



