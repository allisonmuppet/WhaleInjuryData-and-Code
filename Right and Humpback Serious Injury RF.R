# Model to predict outcomes for Right and Humpback Whale injuries based on known-outcome cases.
# Data for NEFSC Right Whales
# Script updated 12/31/2021, using Allison Henry's file dated 7/12/2020
# and loads package 'WhaleInjuryCovariates'
# and using Eric Archer's new GitHub version of rfPermute ()
# Added variable 'mother' 10/16/2020
# Separate covariate suites are used for entanglement and vessel strike models.

# make sure you have devtools installed
if (!require('devtools')) install.packages('devtools')

# install from GitHub
  devtools::install_github('EricArcher/rfPermute', force=TRUE)
  devtools::install_github('JimCarretta/WhaleInjuryCovariates', force=TRUE)

library(rfPermute)
library(tidyverse)
library(ggplot2)
library(ggforce)
library(devtools)
library(gridExtra)
library(WhaleInjuryCovariates)
library(grid)

 rm(list=ls())
 
 unlink(dir(pattern="output.dat"))
 
 size.RF = 1000       # how many RF trees to build
 model.region = "ALL"  # Select region-specific model? "ATL", "PAC", "ALL" ...
 
 setwd("c:/carretta/noaa/net_mort/serious injury/Serious Injury Proration RF")
 filenames = dir(pattern="2020 Proration Recalculation -")
 sp.prefix = substr(filenames, 32, 35)
 
 for (s in 1:length(filenames)) {  
   
   set.seed(1234)

 if (sp.prefix[s]=="HUWH") species="Humpback Whale"
 if (sp.prefix[s]=="RIWH") species="Right Whale"

  df = read.csv(filenames[s], header=T, na.strings=c("NA",""), stringsAsFactors = FALSE, blank.lines.skip = TRUE)
  
  df = WhaleInjuryCovariates(df)
  
  df$VessSz <- factor(df$VessSz)
  df$VessSpd <- factor(df$VessSpd)

 if (model.region=="ATL") data = data[data$Region=="ATL",]
 if (model.region=="PAC") data = data[data$Region=="PAC",]
  
# Disentanglement effort made? These cases are omitted from model-building as they are 'biased' due to intervention. 
# Such cases are part of novel data for which predictions are made.
  
  disentangle = grepl("disentangle|rescue|removed.*gear|gear.*removed|shorten|partial.*disentangl|telemetry", df$Narrative, ignore.case=T)
   disentangle = as.numeric(lapply(disentangle, as.numeric))
    intervention.cases = df[which(disentangle==1),]
    
  data.new = df[which(disentangle==0),]
   data.new$Health.status[grep("dead|decline", data.new$Health.status, ignore.case=TRUE)] <- "DEAD.DECLINE"
    data.new$Health.status[grep("recovered", data.new$Health.status, ignore.case=TRUE)] <- "RECOVERED"
     data.new$Health.status[grep("unknown", data.new$Health.status, ignore.case=TRUE)] <- "UNKNOWN"
  
  write.csv(intervention.cases, paste(species, ".Intervention.Cases.csv", sep=""), row.names = F)

 known = grep("DEAD.DECLINE|RECOVERED", data.new$Health.status, ignore.case=T)
 unknown = grep("UNKNOWN", data.new$Health.status, ignore.case=T)
 data.known = data.new[known,]
 data.unknown = data.new[unknown,]

 dead.decline.ind = grep("DEAD.DECLINE", data.known$Health.status, ignore.case=T)
 recovered.ind = grep("RECOVERED", data.known$Health.status, ignore.case=T)
 
 data.known$Health.status <- factor(data.known$Health.status)

 data.model = data.known
 data.test = data.unknown

 # output data used in model

 write.csv(data.model, paste(species, ".data.model.csv", sep=""), row.names=F)

# identify covariate columns used in randomForest function

 cols.covariates = which(names(data.model)%in%c("anchored", "calf.juv", "constricting", "decline", "extensive.severe", "fluke.peduncle", "gear.free", "head", "healing", "laceration.deep", 
                                               "laceration.shallow", "pectoral", "swim.dive", "trailing", "VessSpd", "VessSz", "wraps.multi", "wraps.no"))
 
# Use different covariate suites for entanglement + vessel strike models
 
 entangle.covariates = which(names(data.model)%in%c("anchored", "calf.juv", "constricting", "decline", "extensive.severe", "fluke.peduncle", "gear.free", "head", "healing", "laceration.deep", 
                                                    "laceration.shallow", "pectoral", "swim.dive", "trailing", "wraps.multi", "wraps.no"))
 
 vessel.covariates = which(names(data.model)%in%c("calf.juv", "decline", "extensive.severe", "fluke.peduncle", "gear.free", "head", "healing", "laceration.deep", 
                                                  "laceration.shallow", "pectoral", "swim.dive", "trailing", "VessSpd", "VessSz"))
 

# parse data into entanglement and vessel models

 data.entangle = data.model[which(data.model$CAUSE%in%c("EN","ET")),]
 write.csv(data.entangle, paste(species, ".data.entangle.csv", sep=""), row.names = F)

 data.vessel = data.model[which(data.model$CAUSE%in%c("VS")),]
 write.csv(data.vessel, paste(species, ".data.vessel.csv", sep=""), row.names=F)

 data.test.entangle = data.test[which(data.test$CAUSE%in%c("EN", "ET")),]
 data.test.vessel = data.test[which(data.test$CAUSE%in%c("VS")),]
 
 ################################## Create randomForest models ###############################################
 
 #### Entanglement Model
 
 sampsize = balancedSampsize(data.entangle$Health.status)

 sink(paste0(sp.prefix[s],".output.dat"))
 model.entangle = rfPermute(data.entangle$Health.status ~ ., data.entangle[,c(entangle.covariates)], sampsize=sampsize, ntree=size.RF, replace=FALSE, importance=TRUE, proximity=TRUE)
 model.entangle
 print(paste0(species, " Entanglement Model"))
 print(model.entangle)
 print("#")
 print("#")
 print(classPriors(model.entangle, sampsize=sampsize))
 print("#")
 print("#")

 pdf(file=paste(species, "all.plots.pdf"), width=8, height=6)
 plotImportance(model.entangle, imp.type="MeanDecreaseAccuracy", main=paste(species, "Entanglement Variable Importance"))
# plotImpPreds(model.entangle, data.entangle, class.col="Health.status", violin.alpha=1)
 plotInbag(model.entangle, sampsize=sampsize, replace=FALSE)
 plotTrace(model.entangle)
 plotPredictedProbs(model.entangle)

 data.entangle.output = data.frame(cbind(casePredictions(model.entangle)$predicted, casePredictions(model.entangle)$DEAD.DECLINE, casePredictions(model.entangle)$RECOVERED, data.entangle))
 names(data.entangle.output)[1:3]=c("Predicted", "DEAD.DECLINE", "RECOVERED")
 write.csv(data.entangle.output, paste(species, ".entangle.training.data.csv", sep=""), row.names = F)

### Use RF 'model.entangle' to predict outcomes for data.test.entangle cases

 entangled.unk.cases.prob = predict(model.entangle, data.test.entangle, "prob")
 entangled.unk.cases.response = predict(model.entangle, data.test.entangle, "response")

# append predictions to data

 data.test.entangle.new = cbind(entangled.unk.cases.response, entangled.unk.cases.prob, data.test.entangle)
 names(data.test.entangle.new)[1:3]=c("Predicted", "DEAD.DECLINE", "RECOVERED")
 write.csv(data.test.entangle.new, paste(species, ".Predicted.Unknown.Entanglement.Cases.csv", sep=""), row.names=F)
 
 ###### Vessel Strike Model
 
 sampsize = balancedSampsize(data.vessel$Health.status)
 
 sink(paste0(sp.prefix[s],".output.dat"), append=TRUE)
 model.vessel = rfPermute(data.vessel$Health.status ~ ., data.vessel[,c(vessel.covariates)], sampsize=sampsize, ntree=size.RF, replace=FALSE, importance=TRUE, proximity=TRUE)
 model.vessel
 print(paste0(species, " Vessel Model"))
 print(model.vessel)
 print("#")
 print("#")
 print(classPriors(model.vessel, sampsize=sampsize))
 sink()
 
 plotImportance(model.vessel, imp.type="MeanDecreaseAccuracy", main=paste(species, "Vessel Variable Importance"))
 plotInbag(model.vessel, sampsize=sampsize, replace=F)
 plotTrace(model.vessel)
 plotPredictedProbs(model.vessel)

 data.vessel.output = data.frame(cbind(casePredictions(model.vessel)$predicted, casePredictions(model.vessel)$DEAD.DECLINE, casePredictions(model.vessel)$RECOVERED, data.vessel))
 names(data.vessel.output)[1:3]=c("Predicted", "DEAD.DECLINE", "RECOVERED")
 write.csv(data.vessel.output, paste(species, ".vessel.training.data.csv", sep=""), row.names = F)  

### Use model.vessel to predict outcomes for data.test.vessel cases

 vessel.unk.cases.prob = predict(model.vessel, data.test.vessel, "prob")
 vessel.unk.cases.response = predict(model.vessel, data.test.vessel, "response")

# append predictions to data

 data.test.vessel.new = cbind(vessel.unk.cases.response, vessel.unk.cases.prob, data.test.vessel)
 names(data.test.vessel.new)[1:3]=c("Predicted", "DEAD.DECLINE", "RECOVERED")
 write.csv(data.test.vessel.new, paste(species, ".Predicted.Unknown.Vessel.Cases.csv", sep=""), row.names=F)
 
 dev.off() 
 
 ## Table 1 of manuscript
 
 EN.data <- data.new[grep("EN|ET", data.new$CAUSE, ignore.case=TRUE),]
 VS.data <- data.new[grep("VS", data.new$CAUSE, ignore.case=TRUE),]
 
 T1.line1 <- cbind.data.frame(
   species,
   "EN",
   nrow(EN.data),
   length(grep("DEAD.DECLINE|RECOVERED", EN.data$Health.status, ignore.case=TRUE)),
   length(grep("DEAD.DECLINE", EN.data$Health.status, ignore.case=TRUE)),
   length(grep("RECOVERED", EN.data$Health.status, ignore.case=TRUE)),
   length(grep("UNKNOWN", EN.data$Health.status, ignore.case=TRUE))
 )
 
 T1.line2 <- cbind.data.frame(
   species,
   "VS",
   nrow(VS.data),
   length(grep("DEAD.DECLINE|RECOVERED", VS.data$Health.status, ignore.case=TRUE)),
   length(grep("DEAD.DECLINE", VS.data$Health.status, ignore.case=TRUE)),
   length(grep("RECOVERED", VS.data$Health.status, ignore.case=TRUE)),
   length(grep("UNKNOWN", VS.data$Health.status, ignore.case=TRUE))
 )
 
 names(T1.line1) <- c("Species", "Injury Type", "Number Cases", "Known Outcome", "Dead.Decline", "Recovered", "Unknown Outcome")
 names(T1.line2) <- names(T1.line1)
 
 Table1 <- rbind.data.frame(T1.line1, T1.line2)
 write.csv(Table1, paste(species, "Table1.csv"), row.names=FALSE)
 
 Table1
 
 
 save.image(paste(species, "RF Serious Injury.RData"))  
 
 sink()   }
 
############## Plot Known + Unknown Outcome Entanglement PBR categories vs RF predictions 
############## of Death.Decline for Right + Humpback whales combined
 
 jpeg(file="RF Predictions vs PBR Assignments Entangle Allspp.jpg", width=2400, height=1200, units="px", res=300)
 
 
### function for use in ggplot to label means and sample sizes
       give.n <- function(x){
           return(c(y = mean(x), label = length(x)))
       }
       
  EG.unk.entangle = read.csv("Right Whale.Predicted.Unknown.Entanglement.Cases.csv")
  MN.unk.entangle = read.csv("Humpback Whale.Predicted.Unknown.Entanglement.cases.csv")

  Allspp.unk.entangle = rbind(EG.unk.entangle, MN.unk.entangle)
  Allspp.unk.entangle$PBR..Value = factor(Allspp.unk.entangle$PBR..Value)
  
  EG.known.entangle = read.csv("Right Whale.entangle.training.data.csv")
  MN.known.entangle = read.csv("Humpback Whale.entangle.training.data.csv")
  
  Allspp.known.entangle = rbind(EG.known.entangle, MN.known.entangle)
  Allspp.known.entangle$PBR..Value = factor(Allspp.known.entangle$PBR..Value)

  p1 <- ggplot(Allspp.known.entangle, aes(factor(PBR..Value), DEAD.DECLINE)) +
    geom_hline(yintercept = 0.50, color="gray") +
     geom_violin(fill="gray") +
      geom_sina()+
       ggtitle("(A) Entanglement: Known Outcomes") +
       labs (y= "RF Probability Dead.Decline", x = "PBR Value Assigned")
  
   p2 <- ggplot(Allspp.unk.entangle, aes(PBR..Value, DEAD.DECLINE)) +
      geom_hline(yintercept = 0.50, color="gray") +
       geom_violin(fill="gray") +
        geom_sina()+
         ggtitle("(B) Entanglement: Unknown Outcomes") +
          labs (y= "RF Probability Dead.Decline", x = "PBR Value Assigned")
   
   jpeg(file="RF Predictions vs PBR Assignments Entanglement Allspp.jpg", width=2400, height=1200, units="px", res=300)
   
   grid.arrange(p1, p2, ncol=2, nrow=1)
   
   dev.off()
  
### Vessel Strike comparisons
  
  EG.unk.vessel = read.csv("Right Whale.Predicted.Unknown.Vessel.Cases.csv")
  MN.unk.vessel = read.csv("Humpback Whale.Predicted.Unknown.Vessel.cases.csv")
  
  Allspp.unk.vessel = rbind(EG.unk.vessel, MN.unk.vessel)
  Allspp.unk.vessel$PBR..Value = factor(Allspp.unk.vessel$PBR..Value)
  
  EG.known.vessel = read.csv("Right Whale.vessel.training.data.csv")
  MN.known.vessel = read.csv("Humpback Whale.vessel.training.data.csv")
  
  Allspp.known.vessel = rbind(EG.known.vessel, MN.known.vessel)
  Allspp.known.vessel$PBR..Value = factor(Allspp.known.vessel$PBR..Value)
  
  p3 <- ggplot(Allspp.known.vessel, aes(factor(PBR..Value), DEAD.DECLINE)) +
    geom_hline(yintercept = 0.50, color="gray") +
     geom_violin(fill="gray") +
      geom_sina()+
    ggtitle("(A) Vessel Strike: Known Outcomes") +
    labs (y= "RF Probability Dead.Decline", x = "PBR Value Assigned")
  
  p4 <-ggplot(Allspp.unk.vessel, aes(PBR..Value, DEAD.DECLINE)) +
    geom_hline(yintercept = 0.50, color="gray") +
     geom_violin(fill="gray") +
      geom_sina()+
     ylim(0,1) +
    ggtitle("(B) Vessel Strike: Unknown Outcomes") +
    labs (y= "RF Probability Dead.Decline", x = "PBR Value Assigned")
  
  jpeg(file="RF Predictions vs PBR Assignments Vessel Allspp.jpg", width=2400, height=1200, units="px", res=300)
  
  grid.arrange(p3, p4, ncol=2, nrow=1)
  
  dev.off()
  
  ### Known entanglements
  
      onePBR.EN.known <- Allspp.known.entangle[Allspp.known.entangle$PBR..Value==1,]
      quantile(onePBR.EN.known$DEAD.DECLINE, c(0.025, 0.5, 0.975))
  
 ### Fraction of unknown entanglement and vessel strike outcomes assigned PBR = 0 with RF probability of DEAD.DECLINE > 0.5
     zeroPBR.EN <- Allspp.unk.entangle[Allspp.unk.entangle$PBR..Value==0,]
     Fr.Gr50.EN <- sum(zeroPBR.EN$DEAD.DECLINE > 0.5) / nrow(zeroPBR.EN)
     signif(Fr.Gr50.EN, 2)
     # sample size
     nrow(zeroPBR.EN)
     quantile(zeroPBR.EN$DEAD.DECLINE, c(0.025, 0.5, 0.975))
     
     zeroPBR.VS <- Allspp.unk.vessel[Allspp.unk.vessel$PBR..Value==0,]
     Fr.Gr50.VS <- sum(zeroPBR.VS$DEAD.DECLINE > 0.5) / nrow(zeroPBR.VS)
     signif(Fr.Gr50.VS, 2)
     # sample size
     nrow(zeroPBR.VS)
     quantile(zeroPBR.VS$DEAD.DECLINE, c(0.025, 0.5, 0.975))
     
 ### Fraction of unknown entanglement and vessel strike outcomes assigned PBR = 1 with RF probability of DEAD.DECLINE <= 0.5
     onePBR.EN <- Allspp.unk.entangle[Allspp.unk.entangle$PBR..Value==1,]
     Fr.LessEqual.50.EN <- sum(onePBR.EN$DEAD.DECLINE <= 0.5) / nrow(onePBR.EN)
     signif(Fr.LessEqual.50.EN, 3)
     quantile(onePBR.EN$DEAD.DECLINE, c(0.025, 0.5, 0.975))
     
     onePBR.VS <- Allspp.unk.vessel[Allspp.unk.vessel$PBR..Value==1,]
     Fr.LessEqual.50.VS <- sum(onePBR.VS$DEAD.DECLINE <= 0.5) / nrow(onePBR.VS)
     signif(Fr.LessEqual.50.VS, 3)
     nrow(onePBR.VS)
     quantile(onePBR.VS$DEAD.DECLINE, c(0.025, 0.5, 0.975))
     
# Fraction of unknown entanglements assigned PBR = 0.75
     
     three.qtrs.EN <- Allspp.unk.entangle[Allspp.unk.entangle$PBR..Value==0.75,]
     Fr.LessEqual.50.EN <- sum(three.qtrs.EN$DEAD.DECLINE <= 0.5) / nrow(three.qtrs.EN)
     signif(Fr.LessEqual.50.EN, 3)
     quantile(three.qtrs.EN$DEAD.DECLINE, c(0.025, 0.5, 0.975))
     

 table.L1 <- cbind.data.frame("Entangle Known",
                       sum(Allspp.known.entangle$PBR..Value==0),
                        sum(Allspp.known.entangle$PBR..Value==0 & Allspp.known.entangle$DEAD.DECLINE>0.5),
                         sum(Allspp.known.entangle$PBR..Value==0.75),
                          sum(Allspp.known.entangle$PBR..Value==0.75 & Allspp.known.entangle$DEAD.DECLINE<=0.5),
                           sum(Allspp.known.entangle$PBR..Value==1),
                            sum(Allspp.known.entangle$PBR..Value==1 & Allspp.known.entangle$DEAD.DECLINE<=0.5),
                              sum(as.numeric(as.character(Allspp.known.entangle$PBR..Value))),
                                sum(Allspp.known.entangle$Predicted=="DEAD.DECLINE"),
                                 sum(Allspp.known.entangle$DEAD.DECLINE)
                        )
 
 names(table.L1) <- c("Injury Type", "PBR=0", "RFprob>0.5", "PBR=0.75", "RFprob<=0.5", "PBR=1", "RFprob<=0.5", "Sum PBR", "Sum RF Majority DEAD.DECLINE", "Sum RF Prob DEAD.DECLINE")
 
  
 table.L2 <- cbind.data.frame("Entangle Unknown",
   sum(Allspp.unk.entangle$PBR..Value==0),
    sum(Allspp.unk.entangle$PBR..Value==0 & Allspp.unk.entangle$DEAD.DECLINE>0.5),
      sum(Allspp.unk.entangle$PBR..Value==0.75),
       sum(Allspp.unk.entangle$PBR..Value==0.75 & Allspp.unk.entangle$DEAD.DECLINE<=0.5),
        sum(Allspp.unk.entangle$PBR..Value==1),
         sum(Allspp.unk.entangle$PBR..Value==1 & Allspp.unk.entangle$DEAD.DECLINE<=0.5),
          sum(as.numeric(as.character(Allspp.unk.entangle$PBR..Value))),
            sum(Allspp.unk.entangle$Predicted=="DEAD.DECLINE"),
              sum(Allspp.unk.entangle$DEAD.DECLINE))
 
 table.L3 <- cbind.data.frame("Vessel Known",
                              sum(Allspp.known.vessel$PBR..Value==0),
                              sum(Allspp.known.vessel$PBR..Value==0 & Allspp.known.vessel$DEAD.DECLINE>0.5),
                              sum(Allspp.known.vessel$PBR..Value==0.75),
                              sum(Allspp.known.vessel$PBR..Value==0.75 & Allspp.known.vessel$DEAD.DECLINE<=0.5),
                              sum(Allspp.known.vessel$PBR..Value==1),
                              sum(Allspp.known.vessel$PBR..Value==1 & Allspp.known.vessel$DEAD.DECLINE<=0.5),
                              sum(as.numeric(as.character(Allspp.known.vessel$PBR..Value))),
                              sum(Allspp.known.vessel$Predicted=="DEAD.DECLINE"),
                              sum(Allspp.known.vessel$DEAD.DECLINE))
 
 table.L4 <- cbind.data.frame("Vessel Unknown",
                              sum(Allspp.unk.vessel$PBR..Value==0),
                              sum(Allspp.unk.vessel$PBR..Value==0 & Allspp.unk.vessel$DEAD.DECLINE>0.5),
                              sum(Allspp.unk.vessel$PBR..Value==0.75),
                              sum(Allspp.unk.vessel$PBR..Value==0.75 & Allspp.unk.vessel$DEAD.DECLINE<=0.5),
                              sum(Allspp.unk.vessel$PBR..Value==1),
                              sum(Allspp.unk.vessel$PBR..Value==1 & Allspp.unk.vessel$DEAD.DECLINE<=0.5),
                              sum(as.numeric(as.character(Allspp.unk.vessel$PBR..Value))),
                              sum(Allspp.unk.vessel$Predicted=="DEAD.DECLINE"),
                              sum(Allspp.unk.vessel$DEAD.DECLINE))
 
 names.Comps.Table <- c("Injury Type", "PBR=0", "RFprob>0.5", "PBR=0.75", "RFprob<=0.5", "PBR=1", "RFprob<=0.5", "Sum PBR", "Sum RF Majority DEAD.DECLINE", "Sum RF Prob DEAD.DECLINE")
 names(table.L1) <- names.Comps.Table
 names(table.L2) <- names.Comps.Table
 names(table.L3) <- names.Comps.Table
 names(table.L4) <- names.Comps.Table
 
 Comps.Table <- rbind.data.frame(table.L1, table.L2, table.L3, table.L4)
 
 write.csv(Comps.Table, "Table Compare RF vs PBR assignments.csv", row.names=FALSE)

 # Vessel Strike Unknown Outcome prorated PBR totals
 # include these proration factors
  prorated.VS.PBR <- c(0.14, 0.20, 0.36, 0.52, 0.56)
 # Number cases
  sum(Allspp.unk.vessel$PBR..Value %in% prorated.VS.PBR)
 # Indices of cases
   prorated.VS.PBR.ind <- which(Allspp.unk.vessel$PBR..Value %in% prorated.VS.PBR)
 # sum of PBR assignments
    sum(as.numeric(as.character(Allspp.unk.vessel$PBR..Value[prorated.VS.PBR.ind])))
 # sum of RF majority predictions
    table(Allspp.unk.vessel$Predicted[prorated.VS.PBR.ind])
  
  closeAllConnections()
  
  dev.off()
  
  source("Plots Publication.R")
  
  
    
     

