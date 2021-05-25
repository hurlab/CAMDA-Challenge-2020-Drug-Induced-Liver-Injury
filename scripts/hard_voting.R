


# This script is used to carry out the hard voting benchmarking

hard_voting <- function(){
    
    curr <- rstudioapi::getActiveDocumentContext()$path
    setwd(dirname(curr))
    
    # prediction on test data
    mold_preds <- lapply(list_of_models$mold, function(dilitype){
        print(names(dilitype))
        lapply(dilitype, function(mod){
            tr_data <- mod$trainingData %>% tibble::rownames_to_column('CAM_ID')
            mold_test_temp <- mold_test[names(mold_test) %in% c(names(tr_data))]
            mold_test_temp <- mold_test_temp[complete.cases(mold_test_temp), ]
            p <- predict(mod, mold_test_temp[, !names(mold_test_temp) %in% 'CAM_ID'])
            as.data.frame(cbind(drug=mold_test_temp$CAM_ID, p))
        })
    })
    
    tox_preds <- lapply(list_of_models$tox, function(dilitype){
        print(names(dilitype))
        y <- lapply(dilitype, function(mod){
            
            tr_data <- mod$trainingData %>% tibble::rownames_to_column('CAM_ID')
            tox_test_temp <- tox_test[names(tox_test) %in% c(names(tr_data))]
            tox_test_temp <- tox_test_temp[complete.cases(tox_test_temp), ]
            p <- predict(mod, tox_test_temp[, !names(tox_test_temp) %in% 'CAM_ID'])
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
            
            #print(faers_test_temp[1:3, 1:3])
            
            p <- predict(mod, faers_test_temp[, !names(faers_test_temp) %in% 'CAM_ID']) 
            p <- as.data.frame(cbind(drug=rownames(faers_test_temp), p))
            
            # %>%
            #     tibble::rownames_to_column('drug') %>%
            #     as.data.frame()
            
            p
        })
    })
    
    # lapply(faers_preds, function(i){
    #     
    #     lapply(i, function(ii){
    #         rownames(ii) <- 
    #     })
    # })
    # 
    # 
    # faers_preds$`dili 3`$dili3.svmPoly3.up$drug <- faers_preds$`dili 3`$dili3.rf3.up$drug
    # faers_preds$`dili 5`$dili5.svmPoly5$drug <- faers_preds$`dili 5`$dili5.rpart5.up$drug
    
    
    # =================
    for(ts in expr_test){
        
        expr_preds <- lapply(list_of_models$gene_expr, function(dilitype){
            lapply(dilitype, function(mod){
                
                names(ts)[names(ts) == 'normalized_name'] <- 'CAM_ID'
                
                tr_data <- mod$trainingData %>% tibble::rownames_to_column('CAM_ID')
                expr_test_temp <- ts[names(ts) %in% c(names(tr_data))]
                
                expr_test_temp <- expr_test_temp[complete.cases(expr_test_temp), ]
                p <- predict(mod, expr_test_temp[, !names(expr_test_temp) %in% 'CAM_ID'])
                as.data.frame(cbind(drug=expr_test_temp$CAM_ID, p))
                
                # txx <- ts[names(ts) %in% names(mod$trainingData)]
                # predict(mod, txx, type='prob')
            })
        })
    }
    
    # grab the names by dili type
    tox_names <- lapply(tox_preds, function(nn){unlist(names(nn))})
    mold_names <- lapply(mold_preds, function(nn){unlist(names(nn))})
    faers_names <- lapply(faers_preds, function(nn){unlist(names(nn))})
    expr_names <- lapply(expr_preds, function(nn){unlist(names(nn))})
    
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
    
    # rename the columns
    merged_tox_preds <- lapply(names(tox_names), function(nn){
        names(merged_tox_preds[[nn]])[2:ncol(merged_tox_preds[[nn]])] <- tox_names[[nn]]
        merged_tox_preds[[nn]]
    })
    
    merged_faers_preds <- lapply(names(faers_names), function(nn){
        names(merged_faers_preds[[nn]])[2:ncol(merged_faers_preds[[nn]])] <- faers_names[[nn]]
        merged_faers_preds[[nn]]
    })
    
    merged_mold_preds <- lapply(names(mold_names), function(nn){
        names(merged_mold_preds[[nn]])[2:ncol(merged_mold_preds[[nn]])] <- mold_names[[nn]]
        merged_mold_preds[[nn]]
    })
    
    merged_expr_preds <- lapply(names(expr_names), function(nn){
        names(merged_expr_preds[[nn]])[2:ncol(merged_expr_preds[[nn]])] <- expr_names[[nn]]
        merged_expr_preds[[nn]]
    })
    
    nn <- paste0('dili ', c(1,3,5,6))
    names(merged_tox_preds) <- nn
    names(merged_faers_preds) <- nn
    names(merged_mold_preds) <- nn
    names(merged_expr_preds) <- nn
    
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
    
    
    preds <- list(dili1=dili1_preds, dili3=dili3_preds, dili5=dili5_preds, dili6=dili6_preds)
    
    hard_voting_predictions <- lapply(preds, function(pred){
        
        apply(pred, 1, function(e){
            table(e, useNA = 'no') 
        }) %>%
            t() %>%
            .[, 1:2] %>%
            as.data.frame() %>%
            dplyr::rename(positive=V1, negative=V2) %>%
            mutate(prediction=ifelse(positive >= negative, 1, 0)) %>%
            dplyr::select(prediction)
        
    }) %>%
        do.call(cbind, .)
    
    rownames(hard_voting_predictions) <- dili1_preds$drug
    hard_voting_predictions <- hard_voting_predictions %>%
        tibble::rownames_to_column('CAM_ID')
    names(hard_voting_predictions)[2:ncol(hard_voting_predictions)] <- paste0('DILI', c(1,3,5,6))
    
    
    return(hard_voting_predictions)
    
    
}





