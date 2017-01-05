## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, cache=TRUE)
test=TRUE
CACHE=T

## ----eval=F--------------------------------------------------------------
#  install.packages("BIOMASS")

## ------------------------------------------------------------------------
library(BIOMASS)
require(knitr) # To build tables in this document

## ---- cache=CACHE--------------------------------------------------------
data(KarnatakaForest)
str(KarnatakaForest)
#
data(NouraguesHD)
str(NouraguesHD)

## ---- cache=CACHE--------------------------------------------------------
selecPlot<-KarnatakaForest$plotId%in%c("BSP2","BSP12","BSP14","BSP26","BSP28","BSP30","BSP34","BSP44","BSP63","BSP65")
KarnatakaForestsub<-droplevels(KarnatakaForest[selecPlot,])

## ----eval=test, cache=CACHE----------------------------------------------
Taxo<-correctTaxo(genus=KarnatakaForestsub$genus,species=KarnatakaForestsub$species)
KarnatakaForestsub$genusCorr<-Taxo$genusCorrected
KarnatakaForestsub$speciesCorr<-Taxo$speciesCorrected

## ----eval=test, cache=CACHE----------------------------------------------
APG<-getTaxonomy(KarnatakaForestsub$genusCorr, findOrder =T)
KarnatakaForestsub$familyAPG<-APG$family
KarnatakaForestsub$orderAPG<-APG$order

## ----eval=test, cache=CACHE----------------------------------------------
dataWD<-getWoodDensity(genus=KarnatakaForestsub$genusCorr,
             species=KarnatakaForestsub$speciesCorr,
             stand=KarnatakaForestsub$plotId)

## ----eval=test, cache=CACHE----------------------------------------------
LocalWoodDensity<-data.frame(genus=c("Ziziphus","Terminalia","Garcinia"),
                             species=c("oenopolia","bellirica","indica"),
                             wd=c(0.65,0.72,0.65))

dataWD<-getWoodDensity(genus=KarnatakaForestsub$genusCorr,
             species=KarnatakaForestsub$speciesCorr,
             family=KarnatakaForestsub$familyAPG,
             stand=KarnatakaForestsub$plotID,
             addWoodDensityData=LocalWoodDensity)

## ----eval=test, cache=CACHE----------------------------------------------
# At species level
sum(dataWD$levelWD=="species")
# At genus level
sum(dataWD$levelWD=="genus")
# At plot level
sum(!dataWD$levelWD%in%c("genus","species"))

## ----eval=F, cache=CACHE-------------------------------------------------
#  HDmodel <- modelHD(D=NouraguesHD$D,
#                     H =NouraguesHD$H,
#                     drawGraph=TRUE,
#                     useWeight=TRUE)

## ----echo=F, cache=CACHE-------------------------------------------------
modelHDplot <- function (D, H, method = NULL, useWeight = FALSE, drawGraph = FALSE){
  nbNonNA <- sum(!is.na(H))
  if (nbNonNA < 15) 
    stop(paste("The data has not enough height data (less than 15 non NA)"))
  Hdata <- data.frame(H, D)
  names(Hdata) <- c("H", "D")
  Hdata <- na.omit(Hdata)
  weight <- NULL
  D_Plot <- seq(from = floor(min(Hdata$D)), to = ceiling(max(Hdata$D)), 
                0.5)
  if (useWeight == TRUE) 
    weight <- (Hdata$D^2) * Hdata$H
  if (!is.null(method)) {
    RSElog <- NULL
    if (grepl("log", method)) {
      modSelected <- loglogFunction(Hdata, method)
      RSElog <- summary(modSelected)$sigma
      coeff <- summary(modSelected)$coefficients
      if (method == "log1") {
        Hpredict_plot <- exp(coeff[1] + 0.5 * RSElog^2 + 
                               coeff[2] * log(D_Plot))
        Hpredict <- exp(coeff[1] + 0.5 * RSElog^2 + 
                          coeff[2] * log(Hdata$D))
      }
      if (method == "log2") {
        Hpredict_plot <- exp(coeff[1] + 0.5 * RSElog^2 + 
                               coeff[2] * log(D_Plot) + coeff[3] * log(D_Plot)^2)
        Hpredict <- exp(coeff[1] + 0.5 * RSElog^2 + 
                          coeff[2] * log(Hdata$D) + coeff[3] * log(Hdata$D)^2)
      }
      if (method == "log3") {
        Hpredict_plot <- exp(coeff[1] + 0.5 * RSElog^2 + 
                               coeff[2] * log(D_Plot) + coeff[3] * log(D_Plot)^2 + 
                               coeff[4] * log(D_Plot)^3)
        Hpredict <- exp(coeff[1] + 0.5 * RSElog^2 + 
                          coeff[2] * log(Hdata$D) + coeff[3] * log(Hdata$D)^2 + 
                          coeff[4] * log(Hdata$D)^3)
      }
    }
    if (method == "weibull") {
      modSelected <- weibullFunction(Hdata, weight)
      coeff <- summary(modSelected)$coefficients
      a <- coeff[1]
      b <- coeff[2]
      c <- coeff[3]
      Hpredict_plot <- a * (1 - exp(-(D_Plot/b)^c))
      Hpredict <- a * (1 - exp(-(Hdata$D/b)^c))
    }
    if (method == "michaelis") {
      modSelected <- michaelisFunction(Hdata, weight)
      coeff <- summary(modSelected)$coefficients
      A <- coeff[1]
      B <- coeff[2]
      Hpredict_plot <- SSmicmen(D_Plot, A, B)
      Hpredict <- SSmicmen(Hdata$D, A, B)
    }
    if (drawGraph == TRUE) {
      par(mar = c(5, 5, 3, 3))
      plot(Hdata$D, Hdata$H, pch = 20, cex = 0.5, col = "grey50", 
           log = "xy", las = 1, xlab = "D (cm)", ylab = "H (m)", 
           cex.lab = 1.8, cex.axis = 1.5, main = paste("Selected model : ", 
                                                       method), cex.main = 2)
      lines(D_Plot, Hpredict_plot, lwd = 2, col = "blue")
      legend("bottomright", c("Data", "Model selected"), 
             lty = c(3, 1), lwd = c(3, 3), col = c("grey", 
                                                   "blue"), cex = 1.5)
    }
  }
  else {
    mod_log1 <- loglogFunction(Hdata, method = "log1")
    RSElog <- summary(mod_log1)$sigma
    coeff <- summary(mod_log1)$coefficients
    Hpredict_log1_plot <- exp(coeff[1] + 0.5 * RSElog^2 + 
                                coeff[2] * log(D_Plot))
    Hpredict_log1 <- exp(coeff[1] + 0.5 * RSElog^2 + coeff[2] * 
                           log(Hdata$D))
    mod_log2 <- loglogFunction(Hdata, method = "log2")
    RSElog <- summary(mod_log2)$sigma
    coeff <- summary(mod_log2)$coefficients
    Hpredict_log2_plot <- exp(coeff[1] + 0.5 * RSElog^2 + 
                                coeff[2] * log(D_Plot) + coeff[3] * log(D_Plot)^2)
    Hpredict_log2 <- exp(coeff[1] + 0.5 * RSElog^2 + coeff[2] * 
                           log(Hdata$D) + coeff[3] * log(Hdata$D)^2)
    mod_log3 <- loglogFunction(Hdata, method = "log3")
    RSElog <- summary(mod_log3)$sigma
    coeff <- summary(mod_log3)$coefficients
    Hpredict_log3_plot <- exp(coeff[1] + 0.5 * RSElog^2 + 
                                coeff[2] * log(D_Plot) + coeff[3] * log(D_Plot)^2 + 
                                coeff[4] * log(D_Plot)^3)
    Hpredict_log3 <- exp(coeff[1] + 0.5 * RSElog^2 + coeff[2] * 
                           log(Hdata$D) + coeff[3] * log(Hdata$D)^2 + coeff[4] * 
                           log(Hdata$D)^3)
    mod_wei <- weibullFunction(Hdata, weight)
    coeff <- summary(mod_wei)$coefficients
    a <- coeff[1]
    b <- coeff[2]
    c <- coeff[3]
    Hpredict_wei_plot <- a * (1 - exp(-(D_Plot/b)^c))
    Hpredict_wei <- a * (1 - exp(-(Hdata$D/b)^c))
    mod_mich <- michaelisFunction(Hdata, weight)
    coeff <- summary(mod_mich)$coefficients
    A <- coeff[1]
    B <- coeff[2]
    Hpredict_mich_plot <- SSmicmen(D_Plot, A, B)
    Hpredict_mich <- SSmicmen(Hdata$D, A, B)
    par(mar = c(5, 5, 3, 3))
    plot(Hdata$D, Hdata$H, pch = 20, cex = 0.5, col = "grey50", 
         log = "xy", las = 1, xlab = "D (cm)", ylab = "H (m)", 
         cex.lab = 1.8, cex.axis = 1.5, main = "Model comparison", 
         cex.main = 2)
    lines(D_Plot, Hpredict_log1_plot, lwd = 2, col = "blue")
    lines(D_Plot, Hpredict_log2_plot, lwd = 2, col = "green")
    lines(D_Plot, Hpredict_log3_plot, lwd = 2, col = "red")
    lines(D_Plot, Hpredict_wei_plot, lwd = 2, col = "orange")
    lines(D_Plot, Hpredict_mich_plot, lwd = 2, col = "purple")
    legend("bottomright", c("Log 1", "Log 2", "Log 3", "Weibull", 
                            "Michaelis"), lty = c(1, 1, 1, 1, 1), lwd = c(2, 
                                                                          2, 2, 2, 2), cex = 1.5, col = c("blue", "green", 
                                                                                                          "red", "orange", "purple"))
  }

}
HDmodel <- modelHDplot(D=NouraguesHD$D, 
                   H =NouraguesHD$H,
                   drawGraph=TRUE,
                   useWeight=TRUE)

## ---- cache=CACHE--------------------------------------------------------
HDmodel<-modelHD(D=NouraguesHD$D,
                 H=NouraguesHD$H,
                 method="log2",
                 useWeight =TRUE)

## ---- cache=CACHE--------------------------------------------------------
HDmodelPerPlot <- by(NouraguesHD,NouraguesHD$plotId,
                     function(x) modelHD(D=x$D,H=x$H, method="weibull",useWeight =T),
                     simplify=FALSE)                     
RSEmodels<-sapply(HDmodelPerPlot,function(x) x$RSE)
Coeffmodels<-lapply(HDmodelPerPlot,function(x) x$coefficients)
ResHD<-data.frame(Plot=names(unlist(RSEmodels)),
                  a=round(unlist(sapply(Coeffmodels,"[",1)),3),
                  b=round(unlist(sapply(Coeffmodels,"[",2)),3),
                  c=round(unlist(sapply(Coeffmodels,"[",3)),3),
                  RSE=round(unlist(RSEmodels),3))
kable(ResHD, row.names = F)

## ---- cache=CACHE--------------------------------------------------------
dataHlocal<-retrieveH(D=KarnatakaForestsub$D,
                      model =HDmodel)

## ---- cache=CACHE--------------------------------------------------------
dataHfeld<-retrieveH(D=KarnatakaForestsub$D,
                     region ="SEAsia")

## ---- eval=F, cache=CACHE------------------------------------------------
#  dataHchave<-retrieveH(D=KarnatakaForestsub$D,
#                        coord=cbind(KarnatakaForestsub$long,KarnatakaForestsub$lat))

## ---- cache=CACHE--------------------------------------------------------
KarnatakaForestsub$WD=dataWD$meanWD
KarnatakaForestsub$H=dataHlocal$H
KarnatakaForestsub$Hfeld=dataHfeld$H

## ----warning=F, cache=CACHE----------------------------------------------
AGBtree<-computeAGB(D=KarnatakaForestsub$D,
                    WD=KarnatakaForestsub$WD,
                    H =KarnatakaForestsub$H)

## ----warning=F, cache=CACHE----------------------------------------------
AGBPlotList<-by(KarnatakaForestsub, KarnatakaForestsub$plotId,
                function(x) computeAGB(D=x$D,WD=x$WD,H=x$H),
                simplify=F)
AGBplot<-sapply(AGBPlotList,sum) 

## ----warning=F, eval=F, cache=CACHE--------------------------------------
#  AGBPlotListChave<-by(KarnatakaForestsub, KarnatakaForestsub$plotId,
#                  function(x) computeAGB(D=x$D,WD=x$WD,coord =cbind(x$long, x$lat)),
#                  simplify=F)
#  AGBplotChave<-sapply(AGBPlotListChave,sum)

## ----warning=F, cache=CACHE----------------------------------------------
AGBPlotListFeld<-by(KarnatakaForestsub, KarnatakaForestsub$plotId,
                function(x) computeAGB(D=x$D,WD=x$WD,H=x$Hfeld),
                simplify=F)
AGBplotFeld<-sapply(AGBPlotListFeld,sum) 

## ---- cache=CACHE--------------------------------------------------------
KarnatakaForestsub$sdWD=dataWD$sdWD
KarnatakaForestsub$HfeldRSE=dataHfeld$RSE

## ---- cache=CACHE--------------------------------------------------------
resultMC<-AGBmonteCarlo(D=KarnatakaForestsub$D,WD=KarnatakaForestsub$WD,errWD = KarnatakaForestsub$sdWD,HDmodel=HDmodel,Dpropag ="chave2004")
meanAGBperplot<-by(resultMC$AGB_simu,KarnatakaForestsub$plotId,function(x) mean(apply(x, 2, sum))) 
credperplot<-by(resultMC$AGB_simu,KarnatakaForestsub$plotId,function(x) quantile(apply(x,2,sum, na.rm = T), probs = c(0.025, 0.975))) 
credinf<-sapply(credperplot,"[",1)
credsup<-sapply(credperplot,"[",2)
ord<-order(meanAGBperplot)
plot(meanAGBperplot[ord],pch=20,xlab="Plots",ylab="AGB (Mg/ha)",ylim=c(0,max(credsup)),las=1,cex.lab=1.3)
segments(1:length(ord),credinf[ord],1:length(ord),credsup[ord],col="red")

## ---- cache=CACHE--------------------------------------------------------
resultMC<-by(KarnatakaForestsub, KarnatakaForestsub$plotId,
             function(x) AGBmonteCarlo(D=x$D,WD=x$WD,H=x$H,errWD = x$sdWD,
                                       HDmodel=HDmodel,Dpropag ="chave2004"),
             simplify=F)
meanAGBperplot<-unlist(sapply(resultMC,"[",1))
credperplot<-sapply(resultMC,"[",4)
credinf<-sapply(credperplot,"[",1)
credsup<-sapply(credperplot,"[",2)
ord<-order(meanAGBperplot)
plot(meanAGBperplot[ord],pch=20,xlab="Plots",ylab="AGB (Mg/ha)",ylim=c(0,max(credsup)),las=1,cex.lab=1.3)
segments(1:length(ord),credinf[ord],1:length(ord),credsup[ord],col="red")

## ---- eval=F, cache=CACHE------------------------------------------------
#  resultMC<-by(KarnatakaForestsub, KarnatakaForestsub$plotId,
#               function(x) AGBmonteCarlo(D=x$D,WD=x$WD,errWD=x$sdWD, H=x$Hfeld,
#                                         errH=x$HfeldRSE, Dpropag="chave2004"),
#               simplify=F)
#  meanAGBperplot<-unlist(sapply(resultMC,"[",1))
#  credperplot<-sapply(resultMC,"[",4)
#  credinf<-sapply(credperplot,"[",1)
#  credsup<-sapply(credperplot,"[",2)
#  ord<-order(meanAGBperplot)
#  plot(meanAGBperplot[ord],pch=20,xlab="Plots",ylab="AGB (Mg/ha)",ylim=c(0,max(credsup)),las=1,cex.lab=1.3)
#  segments(1:length(ord),credinf[ord],1:length(ord),credsup[ord],col="red")

## ---- eval=F,cache=CACHE-------------------------------------------------
#  resultMC<-by(KarnatakaForestsub, KarnatakaForestsub$plotId,
#               function(x)AGBmonteCarlo(D=x$D,WD=x$WD,errWD=x$sdWD,
#                                        coord=cbind(x$long,x$lat),
#                                        Dpropag="chave2004"),
#               simplify=F)
#  meanAGBperplot<-unlist(sapply(resultMC,"[",1))
#  credperplot<-sapply(resultMC,"[",4)
#  credinf<-sapply(credperplot,"[",1)
#  credsup<-sapply(credperplot,"[",2)
#  ord<-order(meanAGBperplot)
#  plot(meanAGBperplot[ord],pch=20,xlab="Plots",ylab="AGB (Mg/ha)",ylim=c(0,max(credsup)),las=1,cex.lab=1.3)
#  segments(1:length(ord),credinf[ord],1:length(ord),credsup[ord],col="red")

## ---- cache=CACHE--------------------------------------------------------
NouraguesHD$Hmix<-NouraguesHD$H
NouraguesHD$RSEmix<-0.5
filt<-is.na(NouraguesHD$Hmix)
NouraguesHD$Hmix[filt]<- retrieveH(NouraguesHD$D,model = HDmodel)$H[filt]
NouraguesHD$RSEmix[filt]<-HDmodel$RSE

## ----eval=F, cache=CACHE-------------------------------------------------
#  resultMC<-by(NouraguesHD, NouraguesHD$plotId,
#               function(x)AGBmonteCarlo(D=x$D,WD=x$WD,errWD=x$sdWD,
#                                        H=NouraguesHD$Hmix,errH=NouraguesHD$RSEmix,
#                                        Dpropag="chave2004"),
#               simplify=F)
#  meanAGBperplot<-unlist(sapply(resultMC,"[",1))
#  credperplot<-sapply(resultMC,"[",4)
#  credinf<-sapply(credperplot,"[",1)
#  credsup<-sapply(credperplot,"[",2)
#  ord<-order(meanAGBperplot)
#  plot(meanAGBperplot[ord],pch=20,xlab="Plots",ylab="AGB (Mg/ha)",ylim=c(0,max(credsup)),las=1,cex.lab=1.3)
#  segments(1:length(ord),credinf[ord],1:length(ord),credsup[ord],col="red")

