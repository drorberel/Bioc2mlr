#' Functionals
#'
#' @param task task
#' @return \code{task}
#' @name functional_helpers

NULL


#' @details Convert functional task, to non-functional task
#' @examples
#' fd.features = list("UVVIS" = 3:136, "NIR" = 137:367)
#' fdf = makeFunctionalData(df, fd.features = fd.features)
#' tsk_functional = makeRegrTask("fuelsubset", data = fdf, target = "heatan")
#' tsk_functional
#' tsk_non_functional<-functional_to_NonFunctional_task_function(tsk_functional)



#' @export
#' @rdname functional_helpers

functional_to_NonFunctional_task_function<-function(task_functional){
  # task_functional = task_SE_Functional
  Matrix<-getTaskData(task_functional, functionals.as = "dfcols") # keep matrix

  if( task_functional$type =='classif' ) task_SE_NON<-makeClassifTask(data = Matrix, target = task_functional %>% getTaskTargetNames)
  ## Regr
  if(task_functional$type =='regr')      task_SE_NON<-makeRegrTask(   data = Matrix, target = task_functional %>% getTaskTargetNames)
  return(task_SE_NON)
}



#' @export
#' @rdname functional_helpers
next_F<-function(task){
  task
}







#' @export
#' @rdname functional_helpers

task_classif_to_Unsupervised<-function(task_non_functional){
  # task_non_functional = iris.task
  # task = task_SE_Functional
  Matrix<-getTaskData(task_non_functional)#, functionals.as = "dfcols") # keep matrix

  task_unsupervised<-makeClusterTask(data = Matrix)
  return(task_unsupervised)
}






#' @export
#' @rdname functional_helpers
Func_task_classif_to_regr<-function(task_class_any){
  # task_class_any = task
  if( (task_class_any %>% getTaskTargets %>% levels %>% length) != 2) 'must be only 2 levels'

  Features<-getTaskData(task_class_any, functionals.as = "matrix") # keep matrix
  Features[,task_class_any %>% getTaskTargetNames] %<>% {ifelse(. == task_class_any$task.desc$positive, 1, 0)}
  # Features$study_randomization %>% table

  ## Regr
  task_regr<-makeRegrTask(data = Features, target = task_class_any %>% getTaskTargetNames)
  task_regr
}




#' @export
#' @rdname functional_helpers
F_mlr_task_add_median_target<-function(task){
  # task = task_j
  DF<-task %>% getTaskData
  target<-task %>% getTaskTargetNames
  DF$y_median<-DF %>% select(!!target) %>% ifelse(.>median(.), 'above','below')

  # back to task
  if(task$type == cont)   task_out<-makeRegrTask(   id = task$id, data = DF, target = target)
  if(task$type == binary) task_out<-makeClassifTask(id = task$id, data = DF, target = target)

}

# task_j %>% F_mlr_task_add_median_target()

