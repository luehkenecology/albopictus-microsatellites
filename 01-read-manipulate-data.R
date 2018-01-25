# working directory--------------------------------------------------------------------------------
setwd("C:/Users/RenkeLuehken/Google Drive/R/microsatellite")
#setwd("G:/NeuAll/R/microsatellite")

# load libraries--------------------------------------------------------------------------------
library(adegenet)
library(plyr)
library(tidyr)

# read data--------------------------------------------------------------------------------
dat <- read.table(file = "data/TEST.csv", sep=";", 
                  head = F, row.names = 1)

# delete rows with NAs
dat <- na.omit(dat)




#write.table(dat,"dat2.txt",row.names = T,col.names=FALSE, quote=FALSE)


# transform to structure format
for(i in 1:(ncol(dat)/2)){
  e <- print(dat[1,(1+(2*(i-1)))])
}

EE <- matrix(nrow = nrow(dat)*2, ncol = (ncol(dat)/2))
for(i in 1:(nrow(dat))){
  EE[(1+(2*(i-1))),] <- as.numeric(dat[i,])[seq(1, ncol(dat), 2)]
  EE[(2+(2*(i-1))),] <- as.numeric(dat[i,])[seq(2, ncol(dat), 2)]
}




population1 <- do.call(rbind, (strsplit(rownames(dat), "[.]")))
population <- apply(population1[,1:2], 1,
                    function(x) paste(x, collapse = "."))                 

sort(unique(population))
populationt <- revalue(population,
                       c("ALF.x" = "1",
                         "E.16" = "2",
                         "E.17" = "3",
                         "F.15" = "4",
                         "F.16" = "5",
                         "F.17" = "6",
                         "H.15" = "7",
                         "H.16" = "8",
                         "H.17" = "9",
                         "HB.17" = "10",
                         "HH.17" = "11",
                         "I.15" = "12",
                         "J.15" = "13",
                         "J.16" = "14",
                         "J.17" = "15",
                         "K.17" = "16",
                         "L.17" = "17",
                         "N.16" = "18",
                         "R.16" = "19",
                         "R.17" = "20",
                         "S.16" = "21",
                         "S.17" = "22"
                         ))
#
#population2 <- population[!(population %in% "ALF.x")]
#dat2 <- dat[!(population %in% "ALF.x"),]

EE2 <- matrix(nrow = nrow(dat)*2, ncol = 2+(ncol(dat)/2))
EE2[,1] <- rep(rownames(dat), each=2)
EE2[,2] <- rep(population, each=2)
EE2[,3:(2+(ncol(dat)/2))] <- EE

EL <- data.frame(EE2)

# remove Aedes albopictus lab colony
EL2 <- subset(EL, !(EL[,2] == 1))

write.table(EL2,"filename.txt",row.names = F,col.names=FALSE, quote=FALSE)

#
# plot (allele frequency)
#library(ggplot2)
#

colnames(EL)
# gather data
EL[,1]
colnames(data_gather)
data_gather <- gather(EL2[,-1], key = locus, 
                      value = allele,
                      X3,
                      X4,
                      X5,
                      X6,
                      X7,
                      X8,
                      X9,
                      X10,
                      X11,
                      X12,
                      X13,
                      X14,
                      X15,
                      X16,
                      X17,
                      X18,
                      -X2,
                      na.rm = T)

A <- ddply(data_gather, .(X2, locus, allele), summarise, freq = length(allele))
B <- ddply(A, .(X2, locus), summarise, sum = sum(freq))
C <- merge(A, B, by = c("X2"))
C$perc <- C$freq/C$sum*100

colnames(C)
C2 <- subset(C, locus.x == c("X3", "X4", "X9", "X18"))

#ggplot(C2, aes(x = as.factor(allele), y = perc, col = X2)) +
#  geom_point()+
#  facet_wrap(~locus.x, scales = "free") +
#  theme_bw()


## genid object
A <- matrix(nrow = nrow(dat), ncol = (ncol(dat)/2))
for(i in 1:(ncol(dat)/2)){
  
  A[,i] <- as.numeric(paste(dat[,(1+(2*(i-1)))], dat[,(2+(2*(i-1)))], sep = ""))
}

B <- data.frame(A)
colnames(B) <- c("CF8", "CF10", "CF12", "CF14",
                 "CF18", "CF22", "CF24", "CF26",
                 "CF28", "CF30","CF32",
                 "CF34",		"CF38","CF40"	,	"CF42",		"CF46")
rownames(B) <- rownames(dat)


#popE <- revalue(population, c("F.15" = "Freiburg",
#                              "F.16" = "Freiburg",
#                              "N.16" = "Novara",
#                              "H.16" = "Heidelberg",
#                              "H.15" = "Heidelberg",
#                              "S.16" = "Sinsheim",
#                              "R.16" = "Rieselfeld"))
library(adegenet)
obj <- df2genind(B, ploidy=2, ncode=3, pop = population2)

library("poppr")
library("ape") # To visualize the tree using the "nj" function
library("magrittr")

png(file = "figure/EE23.jpg",width = 6, height=5, units = 'in', res = 1000)
obj %>% 
  genind2genpop(pop = population2) %>%
  aboot(quiet = TRUE, sample = 1000, distance = nei.dist) 
dev.off()









popE <- revalue(population, c("F.15" = "Freiburg 2015",
                              "F.16" = "Freiburg 2016",
                              "N.16" = "ITA, Novara, 2016",
                              "I.15" = "ITA, Calabria, 2015",
                              "H.16" = "Heidelberg 2016",
                              "H.15" = "Heidelberg 2015",
                              "S.16" = "Sinsheim 2016",
                              "R.16" = "Freiburg, Rieselfeld, 2016",
                              "J.15" = "Jena 2015",
                              "J.16" = "Jena 2016"))
popE2 <- revalue(population, c("ALF.x" = "1",
                               "E.16" = "2",
                               "E.17" = "3",
                               "F.15" = "4",
                               "F.16" = "5",
                               "F.17" = "6",
                               "H.15" = "7",
                               "H.16" = "8",
                               "HB.17" = "9",
                               "HH.17" = "10",
                               "I.15" = "11",
                               "J.15" = "12",
                               "J.16" = "13",
                               "J.17" = "14",
                               "K.17" = "15",
                               "L.17" = "16",
                               "N.16" = "17",
                               "R.16" = "18",
                               "R.17" = "19",
                               "S.16" = "20",
                               "S.17" = "21"))

popE <- revalue(population, c("ALF.x" = "1",
                               "E.16" = "2",
                               "E.17" = "3",
                               "F.15" = "4",
                               "F.16" = "5",
                               "F.17" = "6",
                               "H.15" = "7",
                               "H.16" = "8",
                               "HB.17" = "9",
                               "HH.17" = "10",
                               "I.15" = "11",
                               "J.15" = "12",
                               "J.16" = "13",
                               "J.17" = "14",
                               "K.17" = "15",
                               "L.17" = "16",
                               "N.16" = "17",
                               "R.16" = "18",
                               "R.17" = "19",
                               "S.16" = "20",
                               "S.17" = "21"))

write.table(population,"pop.txt",
            row.names = F,
            col.names=FALSE,
            quote=FALSE)




# plot (allele frequency)
library(ggplot2)
library(tidyr)

colnames(EL)
# gather data
EL[,1]
colnames(data_gather)
data_gather <- gather(EL[,-1], key = locus, 
                      value = allele,
                      X3,
                      X4,
                      X5,
                      X6,
                      X7,
                      X8,
                      X9,
                      X10,
                      X11,
                      X12,
                      X13,
                      X14,
                      X15,
                      X16,
                      X17,
                      X18,
                      -X2,
                      na.rm = T)

A <- ddply(data_gather, .(X2, locus, allele), summarise, freq = length(allele))
B <- ddply(A, .(X2, locus), summarise, sum = sum(freq))
C <- merge(A, B, by = c("X2"))
C$perc <- C$freq/C$sum*100

colnames(C)
C2 <- subset(C, locus.x == c("X3", "X4", "X9", "X18"))

#ggplot(C2, aes(x = as.factor(allele), y = perc, col = X2)) +
#  geom_point()+
#  facet_wrap(~locus.x, scales = "free") +
#  theme_bw()


## genid object
A <- matrix(nrow = nrow(dat), ncol = (ncol(dat)/2))
for(i in 1:(ncol(dat)/2)){
   
  A[,i] <- as.numeric(paste(dat[,(1+(2*(i-1)))], dat[,(2+(2*(i-1)))], sep = ""))
}

B <- data.frame(A)
colnames(B) <- c("CF8", "CF10", "CF12", "CF14",
                 "CF18", "CF22", "CF24", "CF26",
                 "CF28", "CF30","CF32",
                 "CF34",		"CF38","CF40"	,	"CF42",		"CF46")
rownames(B) <- rownames(dat)


#popE <- revalue(population, c("F.15" = "Freiburg",
#                              "F.16" = "Freiburg",
#                              "N.16" = "Novara",
#                              "H.16" = "Heidelberg",
#                              "H.15" = "Heidelberg",
#                              "S.16" = "Sinsheim",
#                              "R.16" = "Rieselfeld"))
library(adegenet)
obj <- df2genind(B, ploidy=2, ncode=3, pop = popE)

rownames(B)
popE

l <- genind2genpop(obj)

#to <- adegenet::summary(obj)

##
##
mat <- as.matrix(obj)
xval <- xvalDapc(mat, popE)

grp <- find.clusters(obj, max.n.clust = 100) # 80 + 4
dapc1 <- dapc(obj, grp$grp) # 80 + 2

o <- data.frame(clust = 1:100, BIC = grp$Kstat)
png(file = "BIC.png",width = 4, height=3 ,
    units = 'in', res = 2000)
ggplot(o, aes(clust, BIC)) +
  geom_point(col = "blue") +
  geom_line(col = "blue") +
  xlab("Number of clusters") +
  geom_vline(xintercept=4, col = "red")+
  theme_classic()
dev.off()

library(grDevices)	
GetColorHexAndDecimal <- function(color)
{
  c <- col2rgb(color)
  sprintf("#%02X%02X%02X %3d %3d %3d", c[1],c[2],c[3], c[1], c[2], c[3])
}
GetColorHexAndDecimal("darkgreen")
GetColorHexAndDecimal("blue")
GetColorHexAndDecimal("gray")
GetColorHexAndDecimal("red")

COL <- revalue(grp$grp, c(
  "1" = "darkgreen",
  "2" = "blue",
  "3" = "gray",
  "4" = "red")
)

myCol <- c("blue","darkgreen","red","gray")

colouP <- c(rep("darkgreen", 58), #F.15+F16
            rep("gray", 33), # R.16
            rep("red", 41), # H.15+H.16
            rep("blue", 32), # S.16
            rep("orange", 32), # N.16
            rep("yellow", 10), # I.15
            rep("darkred", 19)) # J.15+J.16

png(file = "LL22.png",width = 6, height=4 ,
    units = 'in', res = 2000)
scatter(dapc1, scree.da=FALSE,posi.da="bottomright", bg="white",
        pch=15:18, cex = 2, cstar=0, 
        clab=0, cell = 0, col=myCol, 
        scree.pca=F,
        posi.pca = "bottomleft")
par(xpd=TRUE)
points(dapc1$ind.coord[,1], dapc1$ind.coord[,2], 
       pch=21, col = "black",
       bg=colouP)
dev.off()


# map

# Freiburg, Rieselfeld  47.989979, 7.832909
a <- cbind(c(7.843,
             7.843,
             7.792,
             8.652,
             8.652,
             8.875,
             11.598,
             11.598,
             16.837,
             8.634),
            c(48.016,
              48.016,
              47.997,
              49.411,
              49.411,
              49.246,
              50.898,
              50.898,
              39.473,
              45.453))

library(raster)
state.map.DEU <- getData('GADM' , country="DEU", level=0)
crs(state.map.DEU) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 
state.map.ITA <- getData('GADM' , country="ITA", level=0)
crs(state.map.ITA) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 
state.map.CHE <- getData('GADM' , country="CHE", level=0)
crs(state.map.CHE) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 
state.map.AUT <- getData('GADM' , country="AUT", level=0)
crs(state.map.AUT) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 
state.map.CZE <- getData('GADM' , country="CZE", level=0)
crs(state.map.CZE) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 
state.map.NLD <- getData('GADM' , country="NLD", level=0)
crs(state.map.NLD) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 
state.map.BEL <- getData('GADM' , country="BEL", level=0)
crs(state.map.BEL) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 
state.map.LUX <- getData('GADM' , country="LUX", level=0)
crs(state.map.LUX) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 
state.map.FRA <- getData('GADM' , country="FRA", level=0)
crs(state.map.FRA) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 
state.map.POL <- getData('GADM' , country="POL", level=0)
crs(state.map.POL) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 
state.map.SVK <- getData('GADM' , country="SVK", level=0)
crs(state.map.SVK) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 
state.map.HUN <- getData('GADM' , country="HUN", level=0)
crs(state.map.HUN) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 
state.map.HRV <- getData('GADM' , country="HRV", level=0)
crs(state.map.HRV) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 
state.map.SVN <- getData('GADM' , country="SVN", level=0)
crs(state.map.SVN) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 
state.map.BIH <- getData('GADM' , country="BIH", level=0)
crs(state.map.BIH) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 
state.map.MNE <- getData('GADM' , country="MNE", level=0)
crs(state.map.MNE) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 


ab <- rbind(state.map.DEU, state.map.ITA, state.map.CHE,
            state.map.AUT, state.map.CZE, state.map.NLD,
            state.map.LUX, state.map.BEL,state.map.FRA,
            state.map.POL, state.map.SVK, state.map.HUN,
            state.map.HRV,state.map.SVN, state.map.BIH,
            state.map.MNE)

plot(ab)

getwd()
png(file = "germany.png",
    width = 5, height=5, units = 'in', res =1000)
plot(ab)
plot()
points(a, pch = 19, col = "red", cex = 1)
dev.off()

library(raster)
write.table(pointDistance(a, longlat=TRUE)/1000, "du8.csv", sep = ";")



png(file = "LL_legend1.png",width = 6, height=4 ,
    units = 'in', res = 1000)
plot.new()
legend('bottomright',
       'groups',
       c("Cluster 1","Cluster 2","Cluster 3", "Cluster 4"), 
       pch=c(18, 16, 15, 17),
       #lty = c(1,2,3),
       col=c("gray","darkgreen", "blue", "red"))
       #bty ="n")

legend('bottomleft','groups',
       c("GER.FH.15/16", 
         "GER.FR.16", 
         "GER.H.15/16", 
         "GER.S.16", 
         "GER.J.15/16",
         "ITA.M.15", 
         "ITA.N.16"
         ),
       pch=c(rep(21, 4)), col = "black",
       pt.bg=c("darkgreen",
            "gray",
            "red",
            "blue",
            "darkred",
            "yellow",
            "orange"))
dev.off()


scatter(dapc1, pch=21,
        cstar=0,  clab=0,
        mstree=F, scree.da=F)

dapc1$ind.coord
dapc1$grp.coord


dataf3 <- data.frame(eins = dapc1$ind.coord[,1],
                     zwei = dapc1$ind.coord[,2],
                     cluster = dapc1$grp,
                     population = population)

nrow()
ggplot() +
  geom_point(data = dataf3, aes(eins, zwei, 
                                group = cluster, colour = cluster, size = 1.2), size = 1.2) +
  geom_point(data = dataf3, aes(eins, zwei, group = population, colour = population), size = 1.2) +
  geom_vline(xintercept = 0) +
  geom_hline(yintercept = 0) +
  xlab("") +
  ylab("") +
  theme_classic()

png(file = "EE7.jpg",width = 10, height=10, units = 'in', res = 1000)
scatter(dapc1, pch=20,
        cstar=0,  solid=1, cex=3, clab=0,
        mstree=F, scree.da=F)
#plot(type = "n")
points(dapc1$ind.coord[,1], dapc1$ind.coord[,2], 
       pch=21,
       cex=1, lwd=1, col="black",  
       bg = as.numeric(as.factor(popE)))
dev.off()

      

set.seed(4)

png(file = "contri1.jpg",width = 10, height=10, units = 'in', res = 1000)
contrib <- loadingplot(dapc1$var.contr, axis=1,
                       thres=.05, lab.jitter=1)
dev.off()

png(file = "contri2.jpg",width = 10, height=10, units = 'in', res = 1000)
contrib <- loadingplot(dapc1$var.contr, axis=2,
                       thres=.05, lab.jitter=1)
dev.off()

points(dapc1$grp.coord[,1], dapc1$grp.coord[,2], pch=4,
       cex=3, lwd=2, col=myCol)

scatter(dapc1,1,1, col=myCol, bg="white",
        scree.da=FALSE, legend=TRUE, solid=.4)
dapc1$posterior

mat.fst <- pairwise.fst(obj, res.type="matrix")
mat.fst

library(PopGenReport)
#popgenreport(obj, mk.null.all=TRUE, 
#             mk.locihz = T, mk.pdf=TRUE,
#                  path.pgr = getwd())



library("mmod")
data("nancycats")
nancycats


unique(popE2[1:381])

diversity_boot(mlg(obj, quiet = FALSE))

#str(e)
#e$PopHet
#library(PopGenReport)

# allelic richness
#popgenreport(obj, mk.null.all = T)

getwd()
###
###
library(hierfstat)
# observed heterozygosity (H0)
gpg <- basic.stats(genind2hierfstat(obj,pop=NULL))
gpg$overall
# Expected heterozygosity (HE)
Hs(genind2genpop(obj))

# Hardy Weinberg Equilibrium
library(pegas)
hw.test((obj))

# linkage disequilibrium
LD2(genind2loci(obj))


# null alleles
#null.all(obj)

x <- matrix(c("A", "A", "A", "a"), 2)
colnames(x) <- c("Loc1", NA)
y <- alleles2loci(x)
print(y, details = TRUE)

data(gtrunchier)
basic.stats(gtrunchier[,-1])


t <- genind2genpop(obj)
basic.stats(genind2genpop(obj))

basic.stats(t[,-1])


## converting a genind as data.frame 
genind2df(obj)
genind2df(obj, sep="/")
genind2df(obj, oneColPerAll=TRUE)



library(hierfstat)
if(require(adegenet)){
  data(nancycats)
  ## pairwise Fst
  mat.fst <- pairwise.fst(nancycats, res.type="matrix")
  mat.fst
}

?genind2structure 

data(gtrunchier)
?allelic.richness(genind2df(obj))

?data(gtrunchier)
?basic.stats(gtrunchier[,-1])


genind2structure <- function(obj, file="", pops=FALSE){
  if(!"genind" %in% class(obj)){
    warning("Function was designed for genind objects.")
  }
  
  # get the max ploidy of the dataset
  pl <- max(obj@ploidy)
  # get the number of individuals
  S <- nInd(obj)
  # column of individual names to write; set up data.frame
  tab <- data.frame(ind=rep(indNames(obj), each=pl))
  # column of pop ids to write
  if(pops){
    popnums <- 1:nPop(obj)
    names(popnums) <- as.character(unique(pop(obj)))
    popcol <- rep(popnums[as.character(pop(obj))], each=pl)
    tab <- cbind(tab, data.frame(pop=popcol))
  }
  loci <- locNames(obj) 
  # add columns for genotypes
  tab <- cbind(tab, matrix(-9, nrow=dim(tab)[1], ncol=nLoc(obj), dimnames=list(NULL,loci)))
  
  # begin going through loci
  for(L in loci){
    thesegen <- obj@tab[,grep(paste(L, ".", sep=""), dimnames(obj@tab)[[2]], fixed=TRUE)] # genotypes by locus
    al <- 1:dim(thesegen)[2] # numbered alleles
    for(s in 1:S){
      if(all(!is.na(thesegen[s,]))){
        tabrows <- (1:dim(tab)[1])[tab[[1]] == indNames(obj)[s]] # index of rows in output to write to
        tabrows <- tabrows[1:sum(thesegen[s,])] # subset if this is lower ploidy than max ploidy
        tab[tabrows,L] <- rep(al, times = thesegen[s,])
      }
    }
  }
  
  # export table
  write.table(tab, file=file, sep="\t", quote=FALSE, row.names=FALSE)
}







#####
####
#Dear Rudy, 
#
#there is currently no implementation for getting directly a p-value for each Fst. This can be done fairly easily, albeit it will take a bit of time to run. The idea is simple:
#  1) compute original Fst matrix
#2) permute groups randomly
#3) recompute Fst - this is one matrix under Ho: individuals are randomly distributed between the groups (panmixie)
#4) redo 2-3) xxx times to get a full reference distribution.
#5) compare observed values to the ref distribution to get a p-value


#Here is an example using a small dataset:

data(nancycats)
x <- nancycats[sample(1:nrow(nancycats@tab))] # keep 50 individuals
mat.obs <- pairwise.fst(x, res.type="matrix") # this takes about 17 seconds on my computer

NBPERM <- 100 # this is the number of permutations used for the p-values; for a publication at least 999 would be required.
mat.perm <- lapply(1:NBPERM, function(i) pairwise.fst(x, pop=sample(pop(x)), res.type="matrix"))
####

#mat.obs contains original Fst values, mat.perm is a list with NPERM matrices of permuted Fst values. To get e.g. right-tail p-values, we can just count the proportion of mat.obs >= mat.perm; e.g. for the first pair of populations:
  ####
mean(c(mat.obs[1,5] < na.omit(sapply(1:NBPERM, function(i) mat.perm[[i]][1,2])), TRUE))
[1] 0.2
####
In the above command, "na.omit(sapply(1:NBPERM, function(i) mat.perm[[i]][1,2]))" is a vector of permuted values for this pair of populations across all replicates; c(..., TRUE) is added because the observed value is always added to the permuted values (it is one of the possible permutations of groups).

In practice, it is easier to convert the results as objects of the class randtest (class for Monte Carlo test in ade4):
  ####
  library(ade4)
test12 <- as.randtest(na.omit(sapply(1:NBPERM, function(i) mat.perm[[i]][1,2])), mat.obs[1,2], alter="greater")
> test12
Monte-Carlo test
Call: as.randtest(sim = na.omit(sapply(1:NBPERM, function(i) mat.perm[[i]][1, 
                                                                           2])), obs = mat.obs[1, 2], alter = "greater")

Observation: 0.1548673 

Based on 9 replicates
Simulated p-value: 0.2 
Alternative hypothesis: greater 

Std.Obs Expectation    Variance 
1.91197702  0.08374942  0.00138354 
####

To have it done for all pairs of populations:
  ####
  allTests <- list()
for(i in 1:(nrow(mat.obs)-1)){
  for(j in 2:nrow(mat.obs)){
    allTests[[paste(rownames(mat.obs)[i],rownames(mat.obs)[j],sep="-")]] <- as.randtest(na.omit(sapply(1:NBPERM, function(k) mat.perm[[k]][i,j])), mat.obs[i,j], alter="greater")
  }
}
####

allTests is a list with one test of the class "randtest" for each pair of populations. The printing of the object gives the corresponding p-value; note that you can even plot these objects.

As a last note, getting the permuted matrices takes some time: NPERM x 17 seconds in my case would make it 4.7 hours for 1,000 permutations. This could be speeded up drastically by recoding everything in C, but I have no time for this at the moment. To reduce the time taken by the whole stuff, replace lapply by mclapply from the multicore package if you have a multiple core machine. On mine, it decreases the computational time fairly efficiently (from >170 seconds to 52seconds).

Best regards,

Thibaut

library(pvclust)
?pvclust
example(pvclust)

### example using Boston data in package MASS
data(Boston, package = "MASS")
## multiscale bootstrap resampling (non-parallel)
boston.pv <- pvclust(Boston, nboot=100, parallel=FALSE)
## CAUTION: nboot=100 may be too small for actual use.
## We suggest nboot=1000 or larger.
## plot/print functions will be useful for diagnostics.
## plot dendrogram with p-values
plot(boston.pv)

