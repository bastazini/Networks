####################################
#R code to reproduce models tested in Table 2 
#Bastazini, V.A.G., Ferreira, P.M., Azambuja, B.O., Casas, G., Debastiani, V.J., Guimarães, P.R. & Pillar, V.D. 2017. Untangling the Tangled Bank: A Novel Method for Partitioning the Effects of Phylogenies and Traits on Ecological Networks. 
#Evolutionary Biology  44(3): 312–324. DOI: 10.1007/s11692-017-9409-8
#Last updated: 2015-04-15
#Contacts: 
#bastazini.vinicius@gmail.com
#vanderleidebastiani@yahoo.com.br
####################################
###Packages
require(vegan)
require(SYNCSA)
####################################
### Function mantel.residuals
# it extracts residuals from a mantel matrix correlation (see Fig. 2)
####################################
mantel.residuals<-function(A,B){
  a<-as.vector(as.dist(A))
  b<-as.vector(as.dist(B))
  modelo<-lm(a~b)
  RESIDUOS<-residuals(modelo)
  N<-dim(as.matrix(A))[1]
  RES<-matrix(NA,N,N)
  OR<-0
  for(j in 1:N){
    for(i in 1:N){
      if(i==j){
        RES[i,j]=0
      }
      if(i>j){
        OR<-OR+1
        RES[i,j]=RESIDUOS[OR]
      }  
    }	
  }
  R<-sqrt(summary(modelo)$r.squared)
  Dis_RES<-as.dist(RES,diag=T)
  return(list(R=R,Residuals=Dis_RES))
}
####################################
###  Function MATRIX.P1 
#this is a modification of the function matrix.P from SYNCSA
#which allows the computation of matrices with some rows that sum to 0
####################################
matrix.p1<-function (comm, dist.spp) {
  matrix.w <- sweep(comm, 1, rowSums(comm), "/")
  for(i in 1:dim(comm)[1]){
    for(j in 1:dim(comm)[2]){
      if(is.nan(matrix.w[i,j])==TRUE){
        matrix.w[i,j]=as.numeric(0)
      }
    }
  }
  similar.phy <- 1 - (dist.spp/max(dist.spp))
  matrix.phy <- 1/colSums(similar.phy)
  matrix.q <- sweep(similar.phy, 1, matrix.phy, "*")
  matrix.P <- matrix.w %*% matrix.q
  return(list(matrix.w = matrix.w, matrix.q = matrix.q, matrix.P = matrix.P))
}
####################################
######### Data entry #########
####################################
birds_traits #Bird trait data with traits as rows and species as colums 
plants_traits #Plant trait data with traits as rows and species as colums
interactions #Interaction matrix, with plants as rows and birds as colums
birds_phylogeny #A species by species phylogenetic distance matrix
plants_phylogeny #A species by species phylogenetic distance matrix
temporal # A co-occurrence matrix, with plants as rows and birds as colums

############################################################################
######### REMOVING PHYLOGENETIC EFFECT FROM TRAITS #########
############################################################################
sa1<-vegdist(birds_traits,method="gower",na.rm=TRUE);sa1
pa1<-vegdist(birds_phylogeny,method="gower",na.rm=TRUE);pa1
sp1<-vegdist(plants_traits,method="gower",na.rm=TRUE);sp1
pp1<-vegdist(plants_phylogeny,method="gower",na.rm=TRUE);pp1

mant1<-mantel.residuals(sa1,pa1);mant1
mant2<-mantel.residuals(sp1,pp1);mant2
sa<-as.matrix(mant1$Residuals);sa # distance matriz between species described by their traits, removing the effect of the phylogeny 
sb<-as.matrix(mant2$Residuals);sb

##################################################################################################################################
#### PROBABILISTIC INTERACTION MATRICES WEIGHTED BY TRAITS AFTER REMOVING PHYLOGENETIC SIGNAL and accounting for temporal variation
##################################################################################################################################
ua<-matrix.p1(t(interactions),sa)
ya<-ua$matrix.P;ya
up<-matrix.p1(interactions,sb)
yp<-up$matrix.P;yp
ya<-t(ya)*temporal
yp<-yp*temporal
ya<-t(ya)
yp<-t(yp)

##################################################################################################################################
#### PROBABILISTIC INTERACTION MATRICES WEIGHTED BY PHYLOGENETIC RESEMBLANCE and accounting for temporal variation
##################################################################################################################################

dpa<-matrix.p1(t(interactions), birds_phylogeny)
pa<-dpa$matrix.P
dpp<-matrix.p1(interactions, plants_phylogeny)
pp<-dpp$matrix.P
pat<-t(pa)*temporal
pat<-t(pat)
ppt<-(pp)*temporal
ppt<-t(ppt)

##############################
#### MODELS TESTED IN TABLE 2
##############################
#Functional component, removing phylogenetic signal
ya.yp<cor.matrix(ua$matrix.w,ua$matrix.q,ya,yp,method="pearson",dist="euclidean",permutations=9999,norm=FALSE);ya.yp
#Phylogenetic signal in bird species
ya.pa<-cor.matrix(ua$matrix.w,ua$matrix.q,ya,pa,method="pearson",dist="euclidean",permutations=9999,norm=FALSE);ya.pa
#Phylogenetic signal in plant species
yp.pp<-cor.matrix(up$matrix.w,up$matrix.q,t(yp),pp,method="pearson",dist="euclidean",permutations=9999,norm=T);yp.pp
#Phylogenetic component
pa.pp<-cor.matrix(ua$matrix.w,ua$matrix.q,pat,ppt,method="pearson",dist="euclidean",permutations=9999,norm=FALSE);pa.pp
.

