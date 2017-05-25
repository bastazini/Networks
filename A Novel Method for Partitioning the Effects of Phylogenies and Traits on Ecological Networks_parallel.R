####################################
#R code to reproduce models tested in Table 2  in parallel
#Bastazini, V. A., Ferreira, P. M., Azambuja, B. O., Casas, G., Debastiani, V. J., Guimarães, P. R., & Pillar, V. D. (2017). Untangling the Tangled Bank: A Novel Method for Partitioning the Effects of Phylogenies and Traits on Ecological Networks. Evolutionary Biology, 1-13.
#DOI: 10.1007/s11692-017-9409-8
#Last updated: 2015-04-15
#Contacts: 
#bastazini.vinicius@gmail.com
#vanderleidebastiani@yahoo.com.br
####################################
#rm(list=ls())

require(geiger)
require(bipartite)
require(SYNCSA)
require(vegan)
require(mice)

####################################
###### EXTRACTING RESIDUALS #############
####################################
mantel.residuos<-function(A,B){
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
## MODIFIED MATRIX.P1 FUNCTION #####
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
## GENERATING SIMULATED DATA #####
####################################
gera_dados<-function(no_species1, no_species2, no_traits1, no_traits2, alpha1, alpha2){
	tree1<-sim.bdtree(b=0.10, d=0.0,stop="taxa",n=no_species1, extinct=FALSE)
	tree2<-sim.bdtree(b=0.10, d=0.0,stop="taxa",n=no_species2, extinct=FALSE)
	### gerar atributos##
	plantas_atributos<-matrix(NA,no_species1,no_traits1)
	for(i in 1:no_traits1){
		plantas_atributos[,i]<-(rTraitCont(tree1, model = "OU", sigma = 0.1, alpha = alpha1, theta = 0, ancestor = FALSE, root.value = 0))
	}
	aves_atributos<-matrix(NA,no_species2,no_traits2)
	for(i in 1:no_traits2){
		aves_atributos[,i]<-(rTraitCont(tree2, model = "OU", sigma = 0.1, alpha = alpha2, theta = 0, ancestor = FALSE, root.value = 0))
	}    
	###gerar matrix aleatória###
	interacao <- matrix(rbinom((no_species2*no_species1),1,0.5),no_species2,no_species1)
	RES<-list(tree1=tree1,tree2=tree2,plantas_atributos=plantas_atributos,aves_atributos=aves_atributos,interacao=interacao)
return(RES)
}

####################################
## ANALISYS #####
####################################
run<-function(i, no_species11, no_species22, no_traits11, no_traits22, alpha11_min, alpha11_max, alpha22_min, alpha22_max){
#	resultado<-matrix(NA, simulations, 12)
	resultado<-matrix(NA, 1, 12)
	i=1
#	for(i in 1:simulations){
		no_species1<-sample(no_species11,1)
		no_species2<-sample(no_species22,1)
		no_traits1<-sample(no_traits11,1)
		no_traits2<-sample(no_traits22,1)
		alpha1<-runif(1,alpha11_min, alpha11_max)
		alpha2<-runif(1,alpha22_min, alpha22_max)
##################################################################################################################################
		DATA<-gera_dados(no_species1=no_species1, no_species2=no_species2, no_traits1=no_traits1, no_traits2=no_traits2, alpha1=alpha1, alpha2=alpha2)
##################################################################################################################################
		sa1<-as.matrix(vegdist(DATA$aves_atributos,method="euclidean",na.rm=TRUE))
		pa1<-cophenetic(DATA$tree2)
		sp1<-as.matrix(vegdist(DATA$plantas_atributos,method="euclidean",na.rm=TRUE))
		pp1<-cophenetic(DATA$tree1)
		teste1<-mantel.residuos(as.matrix(sa1),pa1)
		teste2<-mantel.residuos(as.matrix(sp1),pp1)
		sa<-as.matrix(teste1$Residuals) # Distance matrix of species described by their traits, withoutphylogenetic signal
		sb<-as.matrix(teste2$Residuals)
##################################################################################################################################
######### RARP = Functional component removing phylogenetic signal ##########
##################################################################################################################################
		uares<-matrix.p1(t(DATA$interacao),sa)
		ra<-uares$matrix.P
		ubres<-matrix.p1(DATA$interacao,sb)
		rb<-t(ubres$matrix.P)
		ra.rb<-cor.matrix(uares$matrix.w,uares$matrix.q,ra,y=rb,method="pearson",dist="euclidean",permutations=9999,norm=FALSE)
##################################################################################################################################
######### PAPP = Phylogenetic component ##########
##################################################################################################################################		
		qa<-matrix.p1(t(DATA$interacao),pa1)
		dpa<-qa$matrix.P
		qb<-matrix.p1((DATA$interacao),pp1)
		dpb<-qb$matrix.P
		pa.pb<-cor.matrix(qa$matrix.w,qa$matrix.q,dpa,y=t(dpb),method="pearson",dist="euclidean",permutations=9999,norm=FALSE)
##################################################################################################################################
######### YAYP = MFunctional component without removing phylogenetic signal ##########
##################################################################################################################################		
		ua<-matrix.p1(t(DATA$interacao),sa1)
		dya<-ua$matrix.P
		ub<-matrix.p1((DATA$interacao),sp1)
		dyb<-ub$matrix.P
		ya.yb<-cor.matrix(ua$matrix.w,ua$matrix.q,dya,y=t(dyb),method="pearson",dist="euclidean",permutations=9999,norm=FALSE)
##################################################################################################################################
		resultado[i,1]<- no_species1
		resultado[i,2]<- no_species2
		resultado[i,3]<- no_traits1
		resultado[i,4]<- no_traits2
		resultado[i,5]<- alpha1
		resultado[i,6]<- alpha2
		resultado[i,7]<- ra.rb$Obs
		resultado[i,8]<- ra.rb$p
		resultado[i,9]<- pa.pb$Obs
		resultado[i,10]<- pa.pb$p
		resultado[i,11]<- ya.yb$Obs
		resultado[i,12]<- ya.yb$p
#	}
	colnames(resultado)=c("no_sp1","no_sp2","no_traits1","no_traits2","alpha1","alpha2","ra.rb_r","ra.rb_p","pa.pb_r","pa.pb_p","ya.yb_r","ya.yb_p")
	rownames(resultado)=rownames(resultado,do.NULL=FALSE,"Sim")
return(resultado)
}
####################################
####################################

run(no_species11=c(5:20), no_species22=c(5:20), no_traits11=c(1:8), no_traits22=c(1:8), alpha11_min=-0.1, alpha11_max=0.1, alpha22_min=-0.1, alpha22_max=0.1)

require(parallel)
parallel = 16
CL <- parallel::makeCluster(parallel, type = "PSOCK")

clusterExport(CL,"gera_dados")
clusterExport(CL,"sim.bdtree")
clusterExport(CL,"rTraitCont")
clusterExport(CL,"vegdist")
clusterExport(CL,"mantel.residuos")
clusterExport(CL,"matrix.p1")
clusterExport(CL,"cor.matrix")

simulations<-50000

RESULTADO<-parallel:: clusterApply(CL, vector("list",simulations), run, no_species11=c(5:20), no_species22=c(5:20), no_traits11=c(1:8), no_traits22=c(1:8), alpha11_min=-0.1, alpha11_max=0.1, alpha22_min=-0.1, alpha22_max=0.1)

parallel::stopCluster(CL)

RESULTADO[[1]]

RESULTADOS0<-t(sapply(seq_len(simulations), function(i) RESULTADO[[i]][1,]))
rownames(RESULTADOS0)<-rownames(RESULTADOS0,FALSE,"Sim.")
RESULTADOS0
head(RESULTADOS0)
tail(RESULTADOS0)

# Rodei ate aqui

# save.image("workspace_2017_01_25")
load("workspace_2017_01_25")

head(RESULTADOS0)
tail(RESULTADOS0)


#Running simulation
#RESULTADO1<-run(no_species11=c(5:20), no_species22=c(5:20), no_traits11=c(1:8), no_traits22=c(1:8), alpha11_min=-0.1, alpha11_max=0.1, alpha22_min=-0.1, alpha22_max=0.1, simulations = 50000)


##################################################################################################################################
######### Type I error ##########
##################################################################################################################################		

#Calculating 99,9%CI
simulations = 50000
#Bicaudal
bicaudal=3.291*(sqrt((0.05*(1-0.05))/simulations))
#Unicaudal
unicaudal=3.090*(sqrt((0.05*(1-0.05))/simulations))

# Example: Rejection rate for the phylogenetic component test
resu=as.data.frame(RESULTADO1)
pa.pb=resu$pa.pb_p<=0.05
rej=table(pa.pb) ["TRUE"]
rej/simulations
