


# create the directories needed


path <- '..'
directories <- c('models', 'predictions', 'evaluations', 'output_files', 'output_plots', 'data')

for(d in directories){
    if(!dir.exists(file.path(path, d))){
        dir.create(file.path(path, d), showWarnings = F)
    }else{
        print('Directory exists.')
    }
}

get_clusters <- function(){
    
    cluster_number <- readline(prompt='Kindly supply the number of cores to use to run the models (The number of cores used\
                               the slower the entire process will take to complete). This should be a digit: ')
    
    if(is.na(as.integer(cluster_number))){
        stop('You must supply a digit. For instance, "2", "45"')
    } else {
        return(as.integer(cluster_number))
    }
}

clusters <- get_clusters()

source('./mold_2.R')
source('./faers.R')
source('./tox_21.R')
source('./gene_expression.R')

print('Completed training models, and saved them in them in the models directory.')




