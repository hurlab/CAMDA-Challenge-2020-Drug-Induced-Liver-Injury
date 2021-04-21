

# a script that holds needed packages.


curr <- rstudioapi::getActiveDocumentContext()$path
setwd(dirname(curr))

#devtools::install_github("RomeroBarata/bimba")


library(VennDetail)
library(ggplot2)
library(caret)
library(caretEnsemble)
library(foreign)
library(readxl)
library(tidyverse)
library(parallel)
library(doParallel)
#library(GeneExpressionSignature)
library(glmnet)
library(mccr)
library(plyr)
library(gbm)
library(DMwR)
library(kernlab)
library(nnet)
library(future)
library(pROC)
library(naivebayes)
library(gridExtra)
library(mboost)
library(MASS)
library(rpart)
library(adabag)
library(mctest)
library(corpcor)
library(ROSE)
library(bimba)
library(mccr)
library(tidyr)



detachAllPackages <- function() {
    
    basic.packages <- c("package:stats","package:graphics","package:grDevices","package:utils","package:datasets","package:methods","package:base")
    
    package.list <- search()[ifelse(unlist(gregexpr("package:",search()))==1,TRUE,FALSE)]
    
    package.list <- setdiff(package.list,basic.packages)
    
    if (length(package.list)>0)  for (package in package.list) detach(package, character.only=TRUE)
    
}

# a function to calculate the matthew's correlation co-efficient ====
mcc <- function(x.model){
    cm <- confusionMatrix(x.model)[['table']]
    true.pos <- cm[1,1]
    false.pos <- cm[1,2]
    false.neg <- cm[2,1]
    true.neg <- cm[2,2]
    
    above <- true.pos*true.neg - false.pos*false.neg
    below <- as.double(true.pos+false.pos)*as.double(true.pos+false.neg)*as.double(true.neg+false.pos)*as.double(true.neg+false.neg)
    
    mcc <- above/sqrt(below)
    mcc
}

collate.metrics <- function(sum.list, x.list, dili.name='', dt='', write.file=F){
    
    # """ This function collates metrics for a list of models """
    
    out.roc <- list()
    out.sens <- list()
    out.spec <- list()
    out.mcc <- list()
    
    for(i in 1:length(sum.list)){
        roc <- as.data.frame(sum.list[[i]]$statistics$ROC[, 'Mean'])
        roc$dataset <- row.names(roc)
        row.names(roc) <- NULL
        out.roc[[i]] <- roc
        
        sens <- as.data.frame(sum.list[[i]]$statistics$Sens[, 'Mean'])
        sens$dataset <- row.names(sens)
        row.names(sens) <- NULL
        out.sens[[i]] <- sens
        
        spec <- as.data.frame(sum.list[[i]]$statistics$Spec[, 'Mean'])
        spec$dataset <- row.names(spec)
        row.names(spec) <- NULL
        out.spec[[i]] <- spec
    }
    
    for(i in 1:length(x.list)){
        mcc.i <- data.frame(dataset=names(x.list)[i], mcc=mcc(x.list[[i]]))
        out.mcc[[i]] <- mcc.i
    }
    
    roc <- do.call(rbind, out.roc)
    names(roc)[1] <- 'ROC'
    sens <- do.call(rbind, out.sens)
    names(sens)[1] <- 'Sens'
    spec <- do.call(rbind, out.spec)
    names(spec)[1] <- 'Spec'
    mcc <- do.call(rbind, out.mcc)
    
    first <- merge(roc, sens, by='dataset')
    mid <- merge(first, spec, by='dataset')
    last <- merge(mid, mcc, by='dataset')
    last <- last %>% arrange(desc(ROC)) 
    
    last$training_type <- ifelse(grepl('.*(\\brose\\b).*|.*(\\bup\\b).*|.*(\\bsmote\\b).*', last$dataset), 'resampled', 'non-resampled')
    names(last)[1] <- 'model'
    
    if(write.file==T){
        write.csv(last, file = paste('../output_files/', dt, '.', dili.name, '.csv', sep=''), row.names=F)
    } else {
        return(last)
    }
}

select.top.three <- function(model.metrics, x.list){
    
    tops <- model.metrics %>% 
        group_by(training_type) %>%
        dplyr::arrange(desc(ROC)) %>%
        filter(row_number() %in% c(1:3))
    
    tt <- x.list[names(x.list) %in% tops$model]
    
    return(list(top_three=tt, top_metrics=tops))
    
}

tt <- select.top.three(tst)

model.list.1[names(model.list.1) %in% tt$model]

# tst$training_type <- ifelse(grepl('.*(\\brose\\b).*|.*(\\bup\\b).*|.*(\\bsmote\\b).*', tst$dataset), 'resampled', 'non-resampled')
# 
# grepl('.*(\\brose\\b).*|.*(\\bup\\b).*|.*(\\bsmote\\b).*', tst$dataset)

training.set %>%
    filter(Training_Validation=='Training Set') %>% .$DILI5 %>%
    table()
