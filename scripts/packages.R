

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
library(RColorBrewer)
library(reshape2)
library(stringr)
library(gridExtra)
library(ggrepel)
library(ggpubr)
library(grid)
library(cowplot)
library(rcartocolor)



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


collate_metrics <- function(models_list){
    
    # get all metrics
    mets <- lapply(models_list, function(dt){
        lapply(dt, function(md){
            
            cc <- caret::confusionMatrix(table(md$pred[['pred']], md$pred[['obs']]))
            cc$byClass %>%
                as.data.frame() %>%
                t() %>%
                as.data.frame() %>%
                dplyr::mutate(Accuracy=cc$overall[['Accuracy']])
            
        })
    }) %>% unlist(., recursive = F) %>%
        do.call(rbind, .)
    
    # get ROC value
    class_sum <- lapply(models_list, function(dt){
        lapply(dt, function(md){
            twoClassSummary(md$pred, lev=levels(md$pred[['obs']]))
        })
    }) %>% unlist(., recursive = F) %>%
        do.call(rbind, .)
    
    # get mccs
    mccs <- lapply(models_list, function(dt){
        lapply(dt, function(md){
            data.frame(mcc=mcc(md))
        })
    }) %>%
        unlist(., recursive = F) %>%
        do.call(rbind, .)
    
    # merge all and select the columns needed
    metrics_df <- merge(mets, class_sum, by=0) %>%
        merge(., mccs, by.x='Row.names', by.y=0) %>%
        separate(col='Row.names', into=c('dili data', 'models'), sep='_') %>%
        mutate(models=sub('\\.', '_', .$models)) %>%
        separate(col='models', into=c('dili type', 'models'), sep='_') %>%
        dplyr::select(1:3, ROC, Sensitivity, Specificity, mcc, `Balanced Accuracy`, Accuracy, Precision, F1) 
    
    metrics_df$`training type` <- ifelse(grepl('.*(\\brose\\b).*|.*(\\bup\\b).*|.*(\\bsmote\\b).*', metrics_df$models), 
                                         'resampled', 'non-resampled')
    
    metrics_df <- relocate(metrics_df, `training type`, .before=ROC) %>%
        dplyr::arrange(desc(ROC))
    
    metrics_df
    
}

select.top.three <- function(model.metrics, x.list, n=3){
    
    tops <- model.metrics %>% 
        group_by(`training type`) %>%
        dplyr::arrange(desc(ROC)) %>%
        filter(row_number() %in% c(1:n))
    
    #x.list <- unlist(x.list, recursive = F)
    
    tt <- x.list[names(x.list) %in% tops$models]
    
    return(list(top_models=tt, top_metrics=tops))
    
}






