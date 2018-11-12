###################################
#Loss of generalist plant species and functional diversity decreases the robustness of a seed dispersal network
#Bastazini V.A.G., V.J. Debastiani, B.O. Azambuja, P.R. Guimarães Jr., & V.D. Pillar (2018). 
#Environmental Conservation: 1-7. doi: 10.1017/S0376892918000334
#Last updated: 2017-06-04
#Contacts: 
#bastazini.vinicius@gmail.com
#vanderleidebastiani@yahoo.com.br
##################################
#packages
require(picante)
require(bipartite)


#### Estimating disctinct scenarios of species extinction and network robustness

## Estimating species extinction based on Functional Distinctiveness
functional=(hclust(dist(atributos)))
tree<-compute.brlen(as.phylo(functional))
plot(tree)

distinc<-evol.distinct(tree,type ="fair.proportion")
distinc

distinc.order<-distinc[order(distinc[,2], decreasing = TRUE), ]
distinc.order$Order<-1:length(tree$tip.label)
distinc.order

rownames(rede)==distinc.order[,1]carnivore
extinct.row<-distinc.order[match(rownames(rede),distinc.order[,1]),3]
extinct.row

rob.functional=second.extinct(rede, participant="lower", method="external", nrep=50,details=FALSE,ext.row=extinct.row)
robustness(rob.functional)
fit.hyperbolica(rob.functional)


##Estimating species extinction based on Phylogenetic Distinctiveness
tree<-rcoal(9) 
tree$tip.label=sample(rownames(Safariland))
filo=(hclust(filo))
tree<-compute.brlen(as.phylo(filo))
plot(tree)

distinc<-evol.distinct(tree,type ="fair.proportion")
distinc

distinc.order<-distinc[order(distinc[,2], decreasing = TRUE), ]
distinc.order$Order<-1:length(tree$tip.label)
distinc.order

rownames(rede)==distinc.order[,1]
extinct.row<-distinc.order[match(rownames(rede),distinc.order[,1]),3]
extinct.row

rob.filogenetica=second.extinct(rede, participant="lower", method="external", nrep=50,details=FALSE,ext.row=extinct.row)
robustness(rob.filogenetica)
fit.hyperbolica(rob.filogenetica)


##Random extinctions
rob.aleat=second.extinct(rede, participant="lower", method="random", nrep=1000,details=FALSE,ext.row=extinct.row)
robustness(rob.aleat)
fit.hyperbolica(rob.aleat)


##Species specialization
rob.grau=second.extinct(rede, participant="lower", method="degree", nrep=1000,details=FALSE,ext.row=extinct.row)
robustness(rob.grau)
fit.hyperbolica(rob.grau)
#
rob.abun=second.extinct(rede, participant="lower", method="abun", nrep=1000,details=FALSE,ext.row=extinct.row)
robustness(rob.abun)
fit.hyperbolica(rob.abun)



##############################
## Figure 2 Phylogenetic and Functional Distinctiveness
##############################

disti=read.table(file.choose(), h=T)
layout(matrix(c(1,2),1,2,byrow=T), widths=c(0.4,0.4))
par(mar=c(5.4,1,1.2,0))
plot(as.phylo(hclust(vegdist(filo,method="euclidean"))),cex=0.7)
par(mar=c(5,0,1,2))
barplot(t(disti), beside = TRUE, horiz=T,axisnames=F, xlim=c(0,0.6), xlab="Distinctiveness",
        legend = c("Phylogenetic","Functional"), args.legend = list(x="bottomright", bty = "n", ncol = 1,cex=0.8))

par(mar=c(5,3,1,2))
barplot(t(disti), beside = TRUE, horiz=T,axisnames=T, xlim=c(0,0.6), xlab="Distinctiveness",las=1)

rownames(filo)==rownames(disti)


##############################
## Figure 3 Attack Tolerance Curves for each scenario
##############################
par(mfrow=c(3,2))
fit.hyperbolica(rob.aleat)
legend("topright", legend = "i", cex=1.2, bty = "n")
fit.hyperbolica(rob.grau)
legend("topright", legend = "ii", cex=1.2,bty = "n")
fit.hyperbolica(rob.filogenetica)
legend("topright", legend = "iii",cex=1.2, bty = "n")
fit.hyperbolica(rob.funcional)
legend("topright", legend = "iv", cex=1.2,bty = "n")
fit.hyperbolica(rob.abun)
legend("topright", legend = "v", cex=1.2,bty = "n")


fit.hyperbolica=function (object, plot.it = TRUE, ...) 
{
  if (class(object) != "bipartite") 
    stop("This function cannot be meaningfully applied to objects of this class.")
  N <- colSums(object)
  if (all(object[-nrow(object), 2] == 1)) 
    y <- -object[, 3]
  else y <- -object[, 2]
  y <- (sum(y) - cumsum(y))/sum(y)
  x <- (object[, "no"]/max(object[, "no"]))
  fit <- try(nls(y ~ 1 - x^a, start = list(a = 1)))
  if (class(fit) == "try-error") 
    fit <- nls((y + rnorm(length(y), s = 0.01)) ~ 1 - x^a, 
               start = list(a = 1))
  if (plot.it) {
    par(mar = c(5, 5, 1, 1))
    plot(x, y, xlab = "Fraction of eliminated plants", 
         ylab = "Fraction of surviving birds", 
         axes = TRUE, type = "n", cex.lab = 1)
    
    
    points(x, y, ...)
    lines(seq(0, 1, 0.1), predict(fit, newdata = data.frame(x = seq(0, 
                                                                    1, 0.1))), col = "red", lwd = 2)
  }
  return(c(exponent = as.numeric(coef(fit)[1])))
}
