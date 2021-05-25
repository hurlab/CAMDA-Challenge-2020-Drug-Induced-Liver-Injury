

# This script is used to carry out the soft voting benchmarking

soft_voting <- function(){
    
    curr <- rstudioapi::getActiveDocumentContext()$path
    setwd(dirname(curr))
    
    # prediction on test data
    mold_preds <- lapply(list_of_models$mold, function(dilitype){
        print(names(dilitype))
        lapply(dilitype, function(mod){
            tr_data <- mod$trainingData %>% tibble::rownames_to_column('CAM_ID')
            mold_test_temp <- mold_test[names(mold_test) %in% c(names(tr_data))]
            mold_test_temp <- mold_test_temp[complete.cases(mold_test_temp), ]
            p <- predict(mod, mold_test_temp[, !names(mold_test_temp) %in% 'CAM_ID'], type='prob')
            as.data.frame(cbind(drug=mold_test_temp$CAM_ID, p))
        })
    })
    
    tox_preds <- lapply(list_of_models$tox, function(dilitype){
        print(names(dilitype))
        y <- lapply(dilitype, function(mod){
            
            tr_data <- mod$trainingData %>% tibble::rownames_to_column('CAM_ID')
            tox_test_temp <- tox_test[names(tox_test) %in% c(names(tr_data))]
            tox_test_temp <- tox_test_temp[complete.cases(tox_test_temp), ]
            p <- predict(mod, tox_test_temp[, !names(tox_test_temp) %in% 'CAM_ID'], type='prob')
            as.data.frame(cbind(drug=tox_test_temp$CAM_ID, p))
        })
    })
    
    
    # set up validation ============================
    # target/label/training
    target.data <- read.csv('../data/targets-camda2020.csv',
                            header=T, stringsAsFactors = F, check.names = F)
    validation.set <- subset(target.data, Training_Validation=='') # 195 obs
    
    faers.data <- read.csv('../data/p2-faers-camda2020.csv',
                           header=T, stringsAsFactors=F, row.names = 'CAM_ID')
    
    faers.validate <- faers.data[rownames(faers.data) %in% validation.set$CAM_ID, ]
    
    faers.validate <- faers.validate %>% 
        mutate(ratio_dili_all=dili_gender_all/all_gender_all, 
               male_rate=(dili_gender_all*dili_gender_male_percentage)/(all_gender_male_percentage*all_gender_all),
               female_rate=(dili_gender_all*dili_gender_female_percentage)/(all_gender_female_percentage*all_gender_all)) 
    # ===============
    
    faers_preds <- lapply(list_of_models$faers, function(dilitype){
        print(names(dilitype))
        y <- lapply(dilitype, function(mod){
            
            tr_data <- mod$trainingData %>% tibble::rownames_to_column('CAM_ID')
            faers_test_temp <- faers.validate[names(faers.validate) %in% c(names(tr_data))]
            faers_test_temp <- faers_test_temp[complete.cases(faers_test_temp), ]
            
            p <- predict(mod, faers_test_temp[, !names(faers_test_temp) %in% 'CAM_ID'], type='prob') %>%
                tibble::rownames_to_column('drug') %>%
                as.data.frame()
            
            p
        })
    })
    
    
    faers_preds$`dili 3`$dili3.svmPoly3.up$drug <- faers_preds$`dili 3`$dili3.rf3.up$drug
    faers_preds$`dili 5`$dili5.svmPoly5$drug <- faers_preds$`dili 5`$dili5.rpart5.up$drug
    
    
    # =================
    for(ts in expr_test){
        
        expr_preds <- lapply(list_of_models$gene_expr, function(dilitype){
            lapply(dilitype, function(mod){
                
                names(ts)[names(ts) == 'normalized_name'] <- 'CAM_ID'
                
                tr_data <- mod$trainingData %>% tibble::rownames_to_column('CAM_ID')
                expr_test_temp <- ts[names(ts) %in% c(names(tr_data))]
                
                expr_test_temp <- expr_test_temp[complete.cases(expr_test_temp), ]
                p <- predict(mod, expr_test_temp[, !names(expr_test_temp) %in% 'CAM_ID'], type='prob')
                as.data.frame(cbind(drug=expr_test_temp$CAM_ID, p))
                
                # txx <- ts[names(ts) %in% names(mod$trainingData)]
                # predict(mod, txx, type='prob')
            })
        })
    }
    
    # fxn to merge the data in the list
    merge_data <- function(df_list, by=''){
        Reduce(function(x, y){
            merge(x=x, y=y, by=by, all=T)
        }, df_list)
    }
    
    #View(merge_data(df_list=mold_preds$`1`, by='drug'))
    
    merged_mold_preds <- lapply(mold_preds, function(x){
        as.data.frame(merge_data(df_list=x, by='drug'))
    })
    merged_tox_preds <- lapply(tox_preds, function(x){
        as.data.frame(merge_data(df_list = x, by='drug'))
    })
    
    merged_faers_preds <- lapply(faers_preds, function(x){
        as.data.frame(merge_data(df_list = x, by='drug'))
    })
    
    merged_expr_preds <- lapply(expr_preds, function(x){
        as.data.frame(merge_data(df_list = x, by='drug'))
    })
    
    # merge by dili type
    dili1_preds <- list(merged_mold_preds$`dili 1`, merged_tox_preds$`dili 1`, merged_faers_preds$`dili 1`, merged_expr_preds$`dili 1`)
    dili1_preds <- merge_data(df_list = dili1_preds, by='drug')
    
    dili3_preds <- list(merged_mold_preds$`dili 3`, merged_tox_preds$`dili 3`, merged_faers_preds$`dili 3`, merged_expr_preds$`dili 3`)
    dili3_preds <- merge_data(df_list = dili3_preds, by='drug')
    
    dili5_preds <- list(merged_mold_preds$`dili 5`, merged_tox_preds$`dili 5`, merged_faers_preds$`dili 5`, merged_expr_preds$`dili 5`)
    #dili5_preds[[1]] <- NULL
    dili5_preds <- merge_data(df_list = dili5_preds, by='drug')
    
    dili6_preds <- list(merged_mold_preds$`dili 6`, merged_tox_preds$`dili 6`, merged_faers_preds$`dili 6`, merged_expr_preds$`dili 6`)
    dili6_preds <- merge_data(df_list = dili6_preds, by='drug')
    
    names(dili1_preds) <- gsub('x0.*|Positive.*', 'Positive', names(dili1_preds))
    names(dili1_preds) <- gsub('x1.*|Negative.*', 'Negative', names(dili1_preds))
    
    names(dili3_preds) <- gsub('x0.*|Positive.*', 'Positive', names(dili3_preds))
    names(dili3_preds) <- gsub('x1.*|Negative.*', 'Negative', names(dili3_preds))
    
    names(dili5_preds) <- gsub('x0.*|Positive.*', 'Positive', names(dili5_preds))
    names(dili5_preds) <- gsub('x1.*|Negative.*', 'Negative', names(dili5_preds))
    
    names(dili6_preds) <- gsub('x0.*|Positive.*', 'Positive', names(dili6_preds))
    names(dili6_preds) <- gsub('x1.*|Negative.*', 'Negative', names(dili6_preds))
    
    n_1 <- list(mold_preds$`dili 1`, tox_preds$`dili 1`, faers_preds$`dili 1`, expr_preds$`dili 1`)
    n_1 <- lapply(n_1, function(dilitype){
        names(dilitype)
    })
    n_1 <- unlist(n_1)
    
    n_3 <- list(mold_preds$`dili 3`, tox_preds$`dili 3`, faers_preds$`dili 3`, expr_preds$`dili 3`)
    n_3 <- lapply(n_3, function(dilitype){
        names(dilitype)
    })
    n_3 <- unlist(n_3)
    
    n_5 <- list(mold_preds$`dili 5`, tox_preds$`dili 5`, faers_preds$`dili 5`, expr_preds$`dili 5`)
    n_5 <- lapply(n_5, function(dilitype){
        names(dilitype)
    })
    n_5 <- unlist(n_5)
    
    n_6 <- list(mold_preds$`dili 6`, tox_preds$`dili 6`, faers_preds$`dili 6`, expr_preds$`dili 6`)
    n_6 <- lapply(n_6, function(dilitype){
        names(dilitype)
    })
    n_6 <- unlist(n_6)
    
    
    n_1 <- do.call(c, lapply(n_1, function(x){
        paste0(x, c('_positive', '_negative'))
    }))
    n_3 <- do.call(c, lapply(n_3, function(x){
        paste0(x, c('_positive', '_negative'))
    }))
    n_5 <- do.call(c, lapply(n_5, function(x){
        paste0(x, c('_positive', '_negative'))
    }))
    n_6 <- do.call(c, lapply(n_6, function(x){
        paste0(x, c('_positive', '_negative'))
    }))
    
    
    names(dili1_preds) <- c('drug', n_1)
    names(dili3_preds) <- c('drug', n_3)
    names(dili5_preds) <- c('drug', n_5)
    names(dili6_preds) <- c('drug', n_6)
    
    preds <- list(dili1=dili1_preds, dili3=dili3_preds, dili5=dili5_preds, dili6=dili6_preds)
    
    # at this point, you can do soft voting and hard voting
    # soft voting, using probabilities
    soft_voting_predictions <- lapply(preds, function(x){
        
        drug <- x$drug
        pos <- x[, grepl('.*_positive$', names(x))] %>%
            rowMeans(., na.rm = T)
        neg <- x[, grepl('.*_negative$', names(x))] %>%
            rowMeans(., na.rm = T)
        
        bb <- cbind(drug, pos, neg) %>%
            as.data.frame() %>%
            mutate(prediction=ifelse(.$pos >= .$neg, 1, 0))
        
        bb$prediction
    }) %>%
        do.call(cbind, .) %>%
        as.data.frame() %>%
        rename_with(toupper)
    
    rownames(soft_voting_predictions) <- dili1_preds$drug
    soft_voting_predictions <- soft_voting_predictions %>%
        tibble::rownames_to_column('CAM_ID')
    
    
    
}

