#' Convert SE to mlr task
#'
#' @param SE_obj SE class
#' @param param.Y.name Vector of dependent variable name
#' @param param.covariates Vector of coaraiate variable(s) name
#' @param param_positive_y_level if ClassifTask, value (character or numeric) to be considered as the positive factor outcome
#' @return mlr's \code{ClassifTask} or \code{RegrTask}
#' @details SE with multiple sub-assays, will include each sub-assay separately.
#'
#'Either ClassifTask or RegrTask will be returned, based on the type of the param.Y.name variable
#' @examples
#' data(airway, package="airway") # RangedSummarizedExperiment
#' airway
#' Fun_SE_to_taskFunc(airway, param.Y.name = 'dex', param.covariates = c('avgLength'), param_positive_y_level = 'trt')
#'
#' @export



Fun_SE_to_taskFunc<-function(SE_obj, param.Y.name, param.covariates, param_positive_y_level, task_return_format, task_type){
  # SE_obj<-SE; param.Y.name = 'vital_status'; param.covariates = c('gender','days_to_death')
  # SE_obj<-smallG_SE; param.Y.name = 'ALL.AML'; param.covariates = c('FAB'); param_positive_y_level = 'ALL'
  # task_return_format = 'functional'
  # task_type = 'classif'


  ## A. Extract data from SE: colData and assay(s)
  DF_ColData<-SE_obj %>% colData %>% data.frame# %>% dplyr::slice(1:10)# %>% select(-!!param.Y.name)
  DF_functionals<-DF_ColData[,c(param.covariates, param.Y.name), drop = FALSE]

  ## all non-numeric coariates must be factors!
  DF_functionals %<>% mutate_if(is.character, as.factor)
  # DF_functionals %>% str

  DF_exprsS<-assays(SE_obj) %>% as.list %>% map(t) ## assume single assay within each SE
  # str(DF_exprsS)
  # DF_exprsS %>% map(dim)
  # DF_exprsS %>% map(rownames) %>% gplots::venn()





  ## B. Returned task. two options: Standard (Non-Functionals), and Functional

  ## B.1. Standard / Non-Functionals:
  DF_non_functional<-cbind(DF_functionals, DF_exprsS$exprs %>% data.frame)

  ## make_task:
  ## Classif
  if( task_type == 'classif' ){
    DF_functionals[,param.Y.name]<-DF_functionals[,param.Y.name] %>% as.factor
    task_SE_NON_functionals<-makeClassifTask(data = DF_non_functional, target = param.Y.name, positive = param_positive_y_level) # with covariates!
  }

  ## Regr
  if( task_type == 'regr' ){
    task_SE_NON_functionals<-makeRegrTask(data = DF_non_functional, target = param.Y.name) # with covariates!
  }


  ## B.2. Functionals:

  ## Add assays to DF_functionals (currently it has only the covariates).
  ## not sure how to do it alternatively to 'for'
  for (i in 1: length(DF_exprsS) ){ # for each assay
    # i=1
    DF_functionals[,names(DF_exprsS)[[i]]]<-DF_exprsS[[i]]
  }
  # str(DF_functionals)

  ## Classif:
  if( task_type == 'classif'){
    DF_functionals[,param.Y.name]<-DF_functionals[,param.Y.name] %>% as.factor
    task_SE_functionals<-makeClassifTask(data = DF_functionals, target = param.Y.name, positive = param_positive_y_level) # with covariates!
  }
    ## Regr
  if( task_type == 'regr'){
     task_SE_functionals<-makeRegrTask(data = DF_functionals, target = param.Y.name) # with covariates!
  }


  if (task_return_format == 'functional') out<-task_SE_functionals
  if (task_return_format == 'dfcols')     out<-task_SE_NON_functionals

  return(out)
}
