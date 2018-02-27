#!/usr/bin/env Rscript
#Sarah Kaddah
#01/02/2018
#Répertoire: results/espece/ setwd(dir="chemin")

#Librairies
library(dplyr)

#Lecture des données (.dat -> classification)
args = commandArgs(trailingOnly=TRUE)
dt = read.table(args[1])
nomFichier = args[1]
print(nomFichier)

#Je veux récupérer un fichier avec les noms des séquences + clusters triés + sans les petites familles pour travailler dessus.
#FamOrdre = dt[order(dt[,2], decreasing = T),]
#FamOrdre = FamOrdre %>% filter(dt[,2]>100)

#Afficher les dimensions du tableau
print (c("Dimension du tableau:",(dim(dt))))

#Nombre de familles
NbFamilles = max(dt[,2])
print(c("Nombre de familles:", NbFamilles))

#Bilan du nombre d'individu par famille
ResumeFamille = data.frame(table(dt[,2]))
colnames(ResumeFamille) <- c("NumCluster","NbSeq")
#print("Nombre d'individus par famille:")
#print(ResumeFamille)
#Ordonner la famille par ordre décroissant
ResFamOrdre = ResumeFamille[order(ResumeFamille[,2], decreasing = T),]
ResFamOrdre[,1] = as.numeric(ResFamOrdre[,1])

#Nombre de familles ayant une taille <100
NbPetitesFamilles = length(which(ResumeFamille[,2]<100))
print(c("Nombre de famille de taille < 100:", NbPetitesFamilles))

#############New
#Nombre de grandes familles
NbGrandesFamilles = NbFamilles - NbPetitesFamilles
print(c("Nombre de famille de taille > 100:", NbGrandesFamilles))

#Tableau contenant les petites familles
PetitesFamilles = ResumeFamille %>% filter(ResumeFamille[,2]<100)
print("Petites familles: cluster associé au nombre de séquences:")
print(PetitesFamilles[order(PetitesFamilles[,2]),])
#Nombre de séquences au total dans les petites familles
TotSeqPF = sum(PetitesFamilles[,2])
print(c("Nombre total de séquences dans les petites familles:", TotSeqPF))
PSeqPF = TotSeqPF*100/nrow(dt)
print(c("Pourcentage de séquences dans les petites familles:",PSeqPF)) 
############NEW
PSeqGF = 100 - PSeqPF
print(c("Pourcentage de séquences dans les grandes familles:",PSeqGF)) 

#Tableau contenant les grandes familles
Familles = ResumeFamille %>% filter(ResumeFamille[,2]>100)


dtOrd = dt[order(dt[,2], decreasing = T),]

#Récupérer un fichier ordonné par ordre décroissant de famille sans les petites familles
#fichierOut = paste(nomFichier,".order",sep=""))
#write.table(FamOrdre, file = paste(nomFichier,".order",sep=""),sep="\t")


########################REPRESENTATIONS GRAPHIQUES###########################

#Histogramme Taille en fonction du Cluster pour tout le jeu de données
png("img/hist_family_distribution_pm090.png")
hist(ResFamOrdre[,2], 
	main="Family size plotted against the cluster", 
	xlab="Cluster size", ylab="Number of sequence",
	border="darkblue", col="lightblue",	
	freq = T, breaks=seq(0,max(ResFamOrdre[,2])+100,100))
dev.off()

# Pour les grosses familles
#png("img/hist_family_distribution_pm090.png")
#hist(ResFamOrdre[,2], 
#	main="Family size plotted against the cluster", 
#	xlab="Cluster size", ylab="Number of sequence",
#	border="darkblue", col="lightblue",	
#	freq = T, breaks=seq(0,max(ResFamOrdre[,2])+2000,2000))
#dev.off()

# xlim=c(0,7000), ylim=c(0,100)
#-----------------------------------------------------------------------------

#Histogramme Taille en fonction du Cluster sans les petites familles
#svg("img/hist_sansPF.svg")
#hist(FamOrdre[,2], 
#	main="Taille en fonction du cluster\nfamilles > 100", 
#	xlab="Taille des clusters", ylab="Nombre de séquences",
#	#border="darkblue", 
#	col="lightblue",	
#	#col.axis = "darkblue", col.lab="darkblue", col.main="darkblue",
#	freq = T, breaks=seq(0,max(ResFamOrdre[,2])+100,100))
#dev.off()
#-----------------------------------------------------------------------------
#library(ggplot2)
#Representation de la distribution sous forme d'un barplot
#ggsave(file="img/barplot45.png",width=12, height=10)
#ggplot(ResumeFamille, aes(x=ResumeFamille[,1], y=ResumeFamille[,2])) + geom_bar(stat="identity", width=1, color="darkblue", fill="lightblue")+theme_minimal()+theme(axis.text.x = element_text(face="bold", color="black", size=8, angle=45), plot.title=element_text(face="bold"))+labs(title="Distribution des familles", x="Clusters", y="Nombre de séquences")
#dev.off()
#-----------------------------------------------------------------------------
########################BROUILLON###########################
#Nombre de séquences dans les familles <100
#ResumeFamilleOrdre = ResumeFamille[order(ResumeFamille[,2]),]
#reperer ceux qui appartiennent a un cluster particulier
#which(test[,2] == 34)
#reperer les indices pour avoir le nom de la sequence
#test[25752,1]
#Distribustion de la taille des famille
#Ordonner le tableau:  
#dtOrdre = dt[order(dt[,2]),]
#Récupérer les données dans un fichier
#write.table(df_out, file=args[2], rown.names=FALSE)
#write.table(ResumeFamille, sep="\t", file="comptage_famille.dat")
