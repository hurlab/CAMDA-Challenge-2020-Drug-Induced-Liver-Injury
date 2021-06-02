

# a script that holds needed packages.
curr <- rstudioapi::getActiveDocumentContext()$path
setwd(dirname(curr))

#devtools::install_github("RomeroBarata/bimba")

# this script is used for loading necessary libraries and package sand installing them if not available
# source it at the beginning of your own main script
# it will also create three folders outside the scripts folders: figures, objects, & data

# download the packages if not available
if (!require("pacman")) {
    install.packages("pacman")
}

pacman::p_load(rstudioapi, doFuture, doSNOW, doParallel, ggplot2,
               plyr, caret, kernlab, tidyverse, parallel,
               gridExtra, pROC, future, mctest,
               tictoc, stringr, foreach, RColorBrewer,
               VennDetail, foreign, readxl, glmnet, mccr,
               gbm, DMwR, nnet, naivebayes, gridExtra, mboost, 
               MASS, rpart, adabag, mctest, corpcor, ROSE, bimba,
               rcartocolor, cowplot, grid, ggpubr, ggrepel, gridExtra, 
               reshape2, VennDetail, devtools)



if(!requireNamespace("BiocManager", quietly = TRUE)){
    install.packages("BiocManager")
}

if(!requireNamespace('VennDetail', quietly = T)){
    BiocManager::install("VennDetail", dependencies=T)
}

# if(!requireNamespace('PGSEA', quietly = T)){
#     BiocManager::install("PGSEA", force = T, dependencies=T, version = '3.12')
# }

if(!requireNamespace('GeneExpressionSignature', quietly = T)){
    stop('Install the package, GeneExpressionSignature.')
    #BiocManager::install("GeneExpressionSignature", force = T, dependencies=T)
}



# load the library(ies) and packages

require(VennDetail)
require(ggplot2)
require(caret)
require(foreign)
require(readxl)
require(tidyverse)
require(parallel)
require(doParallel)
require(GeneExpressionSignature)
require(glmnet)
require(mccr)
require(plyr)
require(gbm)
require(DMwR)
require(kernlab)
require(nnet)
require(future)
require(pROC)
require(naivebayes)
require(gridExtra)
require(mboost)
require(MASS)
require(rpart)
require(adabag)
require(mctest)
require(corpcor)
require(ROSE)
require(bimba)
require(tidyr)
require(RColorBrewer)
require(reshape2)
require(stringr)
require(gridExtra)
require(ggrepel)
require(ggpubr)
require(grid)
require(cowplot)
require(rcartocolor)



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






