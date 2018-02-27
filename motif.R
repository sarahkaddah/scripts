#!/usr/bin/env Rscript
#Sarah Kaddah
#08/02/2018
#Se placer dans le répertoire results/<espece>
#Script qui permet d'étudier les motifs

#=============================================================================
#================================VERSION 2====================================
#=============================================================================

library(dplyr)
args = commandArgs(trailingOnly=TRUE)

#Lecture des données
fichier.cenpb = args[1]
fichier.pja = args[2]
fichier.pkb = args [3]
dt.cenpb = read.table(fichier.cenpb, header = T)
dt.pja = read.table(fichier.pja, header = T)
dt.pkb = read.table(fichier.pkb, header = T)
print(c(fichier.cenpb, fichier.pja, fichier.pkb))

#Enlever les petites familles
dt.cenpb = dt.cenpb %>% filter(dt.cenpb[,2]>100)
dt.pja = dt.pja %>% filter (dt.pja[,2]>100)
dt.pkb = dt.pkb %>% filter (dt.pkb[,2]>100)

#Créer un fichier qui réunit les données
dt = cbind(dt.cenpb[,1], dt.cenpb[,2], dt.cenpb[,3], dt.pja[,3], dt.pkb[,3])
colnames(dt) = c("ID_family", "nb_seq", "CENP-B", "pJalpha", "pKbeta")
write.table("out/motifs.txt", sep = "\t")

#Enlever les familles qui n'ont pas le motif
dt.cenpb = dt.cenpb %>% filter (dt.cenpb[,3]>0)
dt.pja = dt.pja %>% filter (dt.pja[,3]>0)
dt.pkb = dt.pkb %>% filter (dt.pkb[,3]>0)

#Représentation des motifs 
dt1 = cbind(dt[,2], dt[,3])
dt2 = cbind(dt[,2], dt[,4])
dt3 = cbind(dt[,2], dt[,5])
dt1 = dt1[dt1[,2]>0,]
dt2 = dt2[dt2[,2]>0,]
dt3 = dt3[dt3[,2]>0,]
#dt.plot = rbind(dt1[dt1[,2]>0,], dt2[dt2[,2]>0,], dt3[dt3[,2]>0,])
dt.cenpb = dt.cenpb %>% filter (dt.cenpb[,3]>0)

png("img/motifs.png")
plot(dt1[dt1[,2]>0,], ylim = c(0, 1), xlim = c(0, max(dt[,2])), col = "red", xlab = "family size", ylab = "% motif")
points(dt2[dt2[,2]>0,], col = "blue")
points(dt3[dt3[,2]>0,], col = "green")
legend("topright", legend=c("CENP-B", "pJalpha", "pKbeta"),
       col=c("red", "blue", "green"), pch=1)
dev.off()

#==============================Similarity=======================================
#Lecture des données
fichier.similarity = args[4]
dt = read.table(fichier.similarity, header = T)

#Enlever les petites familles
dt = dt %>% filter (dt[,2]>100)

#Représentation de la diversité
png("img/similarity_090.png")
plot(x = dt[,"total_nb_seq"], y = dt[,"mean"], xlab="Taille des familles", ylab="Distance moyenne")
dev.off()


#plot(dt.cenpb[,2], dt.cenpb[,3], xlab = "family size", ylab = "% motif", ylim = c(0.1, 100))
#dev.off()
#png("img/pja_090.png")
#plot(dt.pja[,2], dt.pja[,3], xlab = "family size", ylab = "% motif", ylim = c(0.1, 100))
#dev.off()
#png("img/pkb_090.png")
#plot(dt.pkb[,2], dt.pkb[,3], xlab = "family size", ylab = "% motif", ylim = c(0.1, 100))
#dev.off()

#png("img/motifs_090.png")
#plot(dt[,2],dt[,3], col="red", ylim = c(0.01,1))
#points(dt[,2], dt[,4], col="blue")
#points(dt[,2], dt[,5], col="green")
#dev.off()

#plot(dt[,2],dt[,3], col="red", ylim = c(0,1))
#points(dt[,2], dt[,4], col="blue")
#points(dt[,2], dt[,5], col="green")

#dt23 = cbind(dt[,2], dt[,3])
#dt23[dt23[,2]> 0]
#dt23[dt23[,2]> 0,]
#dt[dt[,3]>0,2]
#dt[dt[,3]>0,3]

#=============================================================================
#================================VERSION 1====================================
#=============================================================================

#dt = read.table("Chlorocebus_Sabaeus/Chlorocebus_Sabaeus_5kmers_pm095_lda40000_fs0.0025_CENPB.motif", header = T)
#dt = read.table("Chlorocebus_Sabaeus/Chlorocebus_Sabaeus_5kmers_pm095_lda40000_fs0.0025_pJalpha.motif", header = T)

#Découpage du tableau en plusieurs parties:
#Valeur nulle pour le motif
#dt.null = dt %>% filter(dt[,3]==0)
#v.null = length(dt.null[,3])

#Supérieur à 0.5
#dt.sup = dt %>% filter((dt[,3]>0.5), (dt[,3]<1))
#v.sup = length(dt.sup[,3])

#Inférieur à 0.5
#dt.inf = dt %>% filter((dt[,3] <= 0.5), (dt[,3] > 0))
#v.inf = length(dt.inf[,3])

#Chez 100% des séquences
#dt.all = dt %>% filter(dt[,3] == 1)
#v.all = length(dt.all[,3])

#Proportion des familles qui ont le motif
#prop = round((v.all + v.inf + v.sup) / length(dt[,3]),2)
#prop.all = round(v.all/ length(dt[,3]),2)
#prop.sup = round(v.sup/ length(dt[,3]),2)
#prop.inf = round(v.inf/ length(dt[,3]),2)
#prop.null = round(v.null/ length(dt[,3]),2) 

#Affichage des résultats
#print("Pourcentage de séquences dans chaque famille ayant le motif")
#print(c("Proportion de familles ayant le motif:", prop))
#print(c("Nombre de familles 100%:", v.all, prop.all))
#print(c("Nombre de familles 0> et > 0.5%:", v.sup, prop.sup))
#print(c("Nombre de familles 1< et < 0.5%:", v.inf, prop.inf))
#print(c("Nombre de familles 0%:", v.null, prop.null))
