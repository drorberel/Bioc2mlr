#' Helper functions
#' @param all_kind all kind...
#'
#' @name Wrapper_Filter_1_helpers
NULL

#' @return all kind...
#' @details TBA



#' @export
#' @rdname Wrapper_Filter_1_helpers

Func_update_args_univ_clusters<-function(lrn, args_vec, lrn.id){
  # lrn = lrn_PreProcess_glmnet; args_vec = list(50, 1, 'method.k', 30, 0.1, 'gene', 5)
  lrn$par.vals<-c(map2(lrn$par.vals[1:7], args_vec, function(.x, .y) .x<-.y), lrn$par.vals[8])
  lrn$id<-lrn.id
  lrn
}




#' @export
#' @rdname Wrapper_Filter_1_helpers

# getTaskSubjectId
getTaskSubjectId<-function(Task){
  Task %>% subsetTask(features=1) %>%  getTaskData() %>% rownames()
}


#' @export
#' @rdname Wrapper_Filter_1_helpers

Func_extract_feature_name<-function(DF_Assay_feature_names, Assay.Analyte.sep){
  # DF_Assay_feature_names=DF_Comb_by_Assay_tib$DF.i[[1]] %>% rownames
  DF_i_gene_loc<-DF_Assay_feature_names %>% str_locate(Assay.Analyte.sep) %>% data.frame %>% pull(end)
  DF_i_gene<-DF_Assay_feature_names %>% substr(DF_i_gene_loc+1, 10000)
  return(DF_i_gene)
}



#' @export
#' @rdname Wrapper_Filter_1_helpers

# Limma / (long) gene level only (with covariates)
Filter_M_Limma_eBayes_gene_REV_function<-function(DF.i, y.either, Covaraites_DF_const){
  # DF.i = x

  # DF.i = Task_Comb_by_Assay_tib_Gene$DF_i_fixed_rowname[[1]] #Task_Comb_by_Assay_tib$DF.i[[1]]
  ## d = getTaskData(task, target.extra = TRUE); DF.i=d$data[,5:10] %>% t %>% data.frame;
  ## DF.i=X.by.assay$DF.i[[1]]; y.either=d$target; Covaraites_DF_const; DF.i[1:3, 1:4]

  ## if samples (train/test) are all NA, remove before run
  observations_to_remove<-DF.i %>% map_lgl(~!all(is.na(.x)))
  DF.i<- DF.i[,observations_to_remove]
  y.either<-y.either[observations_to_remove]
  Covaraites_DF_const<-Covaraites_DF_const[,observations_to_remove]

  ncol(DF.i); length(y.either); ncol(Covaraites_DF_const)

  DF_i_T<-DF.i %>% t %>% data.frame
  my.fml.cont<-formula( '~y.either' )
  Design.cont<-model.matrix(my.fml.cont, data=DF_i_T) # make sure y in numeric, NOT factor!

  if( is.data.frame(Covaraites_DF_const) ){
    Covaraites_DF_const_T<-Covaraites_DF_const %>% t %>% data.frame
    cov.fml.only<-paste(rownames(Covaraites_DF_const), collapse='+')

    my.fml.cont<-formula(paste('~y.either', '+', cov.fml.only))
    Design.cont<-model.matrix(my.fml.cont, data=data.frame(DF_i_T, Covaraites_DF_const_T)) # same gene order. both were taken from DF.Design with no filtering/ordering
  }


  ## limma error messages are listed only at warnings print. require specific 'catch/tryCatch'.
  # Warning message:
  #    Partial NA coefficients for 23 probe(s)
  tryCatch.W.E <- function(expr){
    W <- NULL
    w.handler <- function(w){ # warning handler
      W <<- w
      invokeRestart("muffleWarning")
    }
    list(value = withCallingHandlers(tryCatch(expr, error = function(e) e),warning = w.handler), warning = W)
  }

  fit<-lmFit(DF.i, Design.cont, adjust.method ='fdr') # should always work, but will throw an error

  fit_try<-tryCatch.W.E( eBayes(fit) )
  if( (!is.null(fit_try$value$message)) ) { #$message %>% substr(1, 27))== 'Partial NA coefficients for') {
    fit<-NA
    out<-tibble() # empty tibble iwth 0 rows
  } else {
      fit<-fit_try$value

      if(nrow(DF.i)==1){ # Edge: for a single analyte, tidy(fit) fail to work!
        out<-topTable(fit, coef=2) %>%
          rownames_to_column('gene') %>%
          add_column(order.within.assay=1) %>%
          select(gene, t, P.Value, adj.P.Val, order.within.assay) %>%
          setNames(c('gene', 'estimate', 'p.value', 'p.value.adj', 'order.within.assay')) %>%
          as.tibble
        out$gene<-fit$Amean %>% names

      } else {
        out<-tidy(fit) %>% # filter(gene!=names(Y)) %>% # remove outcome from table, so later on could be matched into rowData
          filter(term=='y.either') %>%  # when covariates added, tidy will accumulate/multiple the table for each of the covariates. take only the ones with Y. It will be adjusted (equivalent to topTable(coef=2))
          mutate(p.value.adj = p.adjust(p.value, 'fdr')) %>%
          mutate(order.within.assay.before.rev = min_rank(p.value)) %>%
          mutate(order.within.assay = - order.within.assay.before.rev) %>% # reverse rank order
          select(gene, estimate, p.value, p.value.adj, order.within.assay) %>%
          arrange(p.value) #%>% head(n=param.Univ.filt.top.n.features) ### select all, order by p-value. top ranked will be selected outside by fw_abs
      }
  } # end tryCatch.W.E
  return(out)
}


#' @export
#' @rdname Wrapper_Filter_1_helpers

### (Short) Univ-glm assays
Filter_M_glm_gene_REV_function<-function(DF.i, y.either, Covaraites_DF_const){
  # DF.i = x
  # DF.i = Task_Comb_by_Assay_tib_Gene$DF_i_fixed_rowname[[1]] # Task_Comb_by_Assay_tib$DF.i[[1]]
  ## d = getTaskData(task, target.extra = TRUE); DF.i=d$data[,5:10] %>% t %>% data.frame;
  ## DF.i=X.by.assay$DF.i[[1]]; y.either=d$target; Covaraites_DF_const; DF.i[1:3, 1:4]

  ## if samples (train/test) are all NA, remove before run
  # DF.i<- DF.i[,DF.i %>% map_lgl(~!all(is.na(.x)))]
  # y.either<-y.either[DF.i %>% map_lgl(~!all(is.na(.x)))]
  ## see alternative fix below

  DF_i_T<-DF.i %>% t %>% data.frame



  My_family<-ifelse( length(table(y.either))==2 , 'binomial', 'gaussian') # & all(names(table(y.either))==c(0,1))

  Short_res_tib<-DF_i_T %>% map(~.x) %>% enframe('gene','X') %>%
    # mutate(frml = Analyte<-formula(paste(colnames(Y),'~',x))
    mutate(X_cov = X %>% map(~if(is.data.frame(Covaraites_DF_const)) data.frame(.x, Covaraites_DF_const %>% t))) %>%
    # Short_res_tib$X_cov[[1]] %>% str
    mutate(GLM_tidy = X_cov %>% map(~glm(y.either~., family=My_family, data = .x) %>% tidy %>% filter(term=='.x'))) %>%
    # .x = Short_res_tib$GLM_tidy[[14]]
    mutate(issue_null_glmfit_to_filter = GLM_tidy %>% map_chr(~ifelse(nrow(.x)==0, 'ISSUE','OK'))) %>%
    filter(issue_null_glmfit_to_filter=='OK') %>%
    mutate(estimate   =GLM_tidy %>% map_dbl(~ .x %>% pull(estimate))) %>%  ### This might throw an error. see how it was handle below for DF/screen approach
    mutate(p.value    =GLM_tidy %>% map_dbl(~ .x  %>% pull(p.value))) %>%
    mutate(p.value.adj = p.adjust(p.value, 'fdr'))

  Short_res_tib %<>% # add ranks, and then reverse order
    mutate(order.within.assay.before.rev = min_rank(p.value)) %>%
    mutate(order.within.assay            = - order.within.assay.before.rev) %>% # reverse rank order
    select(gene, estimate, p.value, p.value.adj, order.within.assay) %>%
    arrange(p.value) #%>% head(n=param.Univ.filt.top.n.features) ### select all, order by p-value. top ranked will be selected outside by fw_abs

  Short_res_tib
}





#' @export
#' @rdname Wrapper_Filter_1_helpers

# Module:
### load GMTs, and keep only relevant ones, and merge all (relevant) sets

Fun_Load_module_sets<-function(mod_all.sets.tib){
  # mydir = 'Z:/R_rhino/JDRF_Ensemble/saved_MAE/'
  # load(file=paste0(mydir, 'mod_all.sets.tib.Rdata')) # mod_all.sets.tib
  # mod_all.sets.tib$GMT.list %>% map(function(x) x[1] %>% map(function(x) x[1:4]))

  #param.Selected.Meta.sets<-c('BMT', 'c5-GO', 'c7-immunologic', 'h-Hallmark')
  modules_sets_tib<-mod_all.sets.tib #%>% filter( Gmt.name %in% param.Selected.Meta.sets) # Filter OUT meta-sets
  ## Before unlist, assure there is no overlap pathway names across all sets, if there is, make the name uniqe, e.g. for each set: unite('set_name', pathawy_name)
  #if( all(unlist(BMT$GMT.list ,recursive = FALSE) %>% names %>% table !=1)) '!!! some overlap pathway names across all sets'

  ## Combine all modules from ALL selected sets
  modules_sets_combined_l<-unlist(modules_sets_tib$GMT.list, recursive = FALSE) # assume list names at each of 4 sublists is unique, and no overlap across 4 lists
  # modules_sets_combined_l %>% head
  return(modules_sets_combined_l)
}
# modules_sets_combined_l<-Fun_Load_module_sets()




#' @export
#' @rdname Wrapper_Filter_1_helpers

##### Helper: ids2indices prep:
F_assay_i_module_row_index<-function(DF_i, modules_sets_combined_l_F){
  # DF_i = DF_by_Assay_tib$DF_i_fixed_rowname[[i]]; modules_sets_combined_l_F = DF_by_Assay_tib$Module_i_list[[i]]

  # ids2indices

  # gene:
  # DF_i=Task_Comb_by_Assay_tib$DF.i[[3]] ; modules_sets_combined_l_F=modules_sets_combined_l; modules_sets_combined_l_F %>% head

  # module
  # DF_i=DF_by_Assay_tib$DF_i[[1]] ; modules_sets_combined_l_F=DF_by_Assay_tib$Module_i_list[[1]]

  # DF_i %>% rownames %>% head


  ## Setup: extract assay name from feature name: can be done externally, though gene names are used only for matching with set names
  # DF_i_gene_loc<-rownames(DF_i) %>% str_locate(Assay.Analyte.sep) %>% data.frame %>% pull(end) ## old version, when it used to include Assay.Analyte.sep. pre-'fixed' gene name is now used
  DF_i_gene<-rownames(DF_i) # %>% substr(DF_i_gene_loc+1, 10000)

  ## Filter 1: pathway data base: only sets with at least 5 genes
  sets_indices_DF_rows<-limma::ids2indices(modules_sets_combined_l_F, DF_i_gene)
  sets_indices_DF_rows_passed<-sets_indices_DF_rows[sets_indices_DF_rows %>% map(length) >= 5]
  return(sets_indices_DF_rows_passed)
}



#' @export
#' @rdname Wrapper_Filter_1_helpers

####### Helper: camera:
Filter_M_Limma_camera_Module_REV_function<-function(DF_i, y.either, Covaraites_DF_const_F, Set_i, assay_i){
  # run camera:
  # DF_i=Task_Comb_by_Assay_tib_Module$DF.i[[1]]; y.either=d$target; Covaraites_DF_const_F<-Covaraites_DF_const; DF_i[1:3, 1:4]; Set_i=Task_Comb_by_Assay_tib_Module$set_i_F1_Symbol_to_row_ideces[[1]][1:5]; assay_i=Task_Comb_by_Assay_tib_Module$assay_i[[1]]

  ## remove subjects with 100% NAs
  crit1<-DF_i %>% map_lgl(~ !all(is.na(.x)))
  DF_i<-DF_i[,crit1]
  y.either<-y.either[crit1]

  DF_i_T<-DF_i %>% t %>% data.frame  # DF_i[1:3, 1:4]
  my.fml.cont<-formula( '~y.either' )
  Design.cont<-model.matrix(my.fml.cont, data=DF_i_T) # make sure y in numeric, NOT factor!

  if( is.data.frame(Covaraites_DF_const_F) ){
    Covaraites_DF_const_F<-Covaraites_DF_const_F[,crit1]
    Covaraites_DF_const_F_T<-Covaraites_DF_const_F %>% t %>% data.frame
    cov.fml.only<-paste(rownames(Covaraites_DF_const_F), collapse='+')

    my.fml.cont<-formula(paste('~y.either', '+', cov.fml.only))
    Design.cont<-model.matrix(my.fml.cont, data=data.frame(DF_i_T, Covaraites_DF_const_F_T)) # same gene order. both were taken from DF.Design with no filtering/ordering
  }

  res.camera<-camera(y=DF_i %>% as.matrix, index=Set_i, design=Design.cont) # rows=genes

  # head(res.camera)
  out<-res.camera %>% rownames_to_column('gene') %>%
    # mutate(gene2 = paste(assay_i %>% as.character(), Assay.Analyte.sep, gene, sep='.')) %>%
    mutate(estimate = NGenes) %>%
    mutate(p.value = PValue) %>%
    mutate(p.value.adj = FDR) %>%
    mutate(order.within.assay.before.rev = min_rank(p.value)) %>%
    mutate(order.within.assay = -order.within.assay.before.rev) %>% # reverse rank order
    #mutate(passed_camrea_top = order.within.assay %>% map_lgl(~.x >= -args$param.Univ.filt.top.n.features)) %>%
    #filter(passed_camrea_top) %>% ## !!!!!! This will cause trouble if later on wish to look back at all tested modules, e.g. for network/graph analysis
    select(gene, estimate, p.value, p.value.adj, order.within.assay)
  return(out)
}

# ff<-Filter_M_Limma_camera_Module_REV_function(DF_i, y.either, Covaraites_DF_const_F, Set_i, assay_i)







#' @export
#' @rdname Wrapper_Filter_1_helpers

# Helper1: PreProc: Train, Gene + module

Func_task_to_Passed_Features_DF_MaG<-function(task.multi, args){  # modules_sets_combined_l is also a parameter. for now will be a global 'floating' variable, but should be included in args!!
  # task.multi = task.int; args = args; #list(param.Univ.filt.top.n.features = param.Univ.filt.top.n.features, param.UnivClustRankTopN = param.UnivClustRankTopN, param.cluster_method_KH = param.cluster_method_KH, param.corrplot.n.clusters.k = param.corrplot.n.clusters.k, param.corrplot.n.clusters.h = param.corrplot.n.clusters.h, parame.gene.or.module = parame.gene.or.module)
  if( !exists('modules_sets_combined_l')) cat('modules sets are not loaded')

  d = getTaskData(task.multi, target.extra = TRUE) # str(d$data) %>% tail
  y.cont<-d$target %>% as.numeric
  y.bin <-d$target %>% factor %>% as.numeric -1
  # Covar<-d$data$Covariate_dem_sex # ! note: univglms does NOT allow covariates!
  if(task.multi$type=='classif') y.either<-y.bin
  if(task.multi$type=='regr')    y.either<-y.cont

  ## Extract covariates:
  Covaraites_location<-which(str_sub(colnames(d$data),1, 10)=='Covariate_')
  if(length(Covaraites_location)>0){
    colnames(d$data)[Covaraites_location]<-gsub('Covariate_', paste0('Covariate', Assay.Analyte.sep), colnames(d$data)[Covaraites_location] )
  }

  ## Step 0: setup: extract assay name from feature name:
  features_Gene_location<-colnames(d$data) %>% str_locate(Assay.Analyte.sep) %>% data.frame %>% pull(start)
  # table(features_Gene_location)
  Task_Comb_Assay_type_v<-colnames(d$data) %>% substr(1, features_Gene_location-1) ### !!! this include the covariates at the END!
  # table(Task_Comb_Assay_type_v)


  ## Step 1: within each assay, select features based on 'top' univariate score (by assay)
  Task_Comb_by_Assay_tib<-data.frame(assay_i = Task_Comb_Assay_type_v, d$data %>% t) %>% rownames_to_column('Analyte.name') %>%
    group_by(assay_i) %>% nest(.key = 'DF.and.Analyte') %>%
    mutate(DF.i = DF.and.Analyte %>% map(~.x %>% data.frame %>% column_to_rownames('Analyte.name') )) %>%
    # .x=Task_Comb_by_Assay_tib$DF.i[[1]]
    mutate(DF_i_fixed_rowname = DF.i %>% map(~.x %>% set_rownames(Func_extract_feature_name(.x %>% rownames, Assay.Analyte.sep) )))

  # Task_Comb_by_Assay_tib$DF.i[[1]][1:3, 1:4]
  # Task_Comb_by_Assay_tib$DF_i_fixed_rowname[[1]][1:3, 1:4]

  if(length(Covaraites_location)>0){
    Covaraites_DF_const<-Task_Comb_by_Assay_tib %>% filter(assay_i=='Covariate') %>% pull(DF_i_fixed_rowname) %>% .[[1]]
    Task_Comb_by_Assay_tib  %<>% filter(assay_i != 'Covariate') # used only to ADJUST the univariate analysis, but will be mandatory added into ML/multivariate step.
  } else Covaraites_DF_const<-NA # rm(Covaraites_DF_const)
  ## Covaraites_DF_const are actually NOT returned. will be added to ML regardless!


  ## add assay type (short/long from external input)

  # External assay table
  # param.assay.type.vec<-c('Short', 'Long', 'Short', 'Short', rep('Long', 11)) # to be automatically generated from MAE/SE_i@metadata
  Assays_used_task_oredered<-task.multi %>% getTaskFeatureNames %>% word(1, sep = Assay.Analyte.sep) %>% unique %>% .[!str_detect(., 'Covariate_')]
  if(length(args$param.assay.type.vec) != length(Assays_used_task_oredered)) print('param.assay.type.vec does NOT match data (task) assay names!')
  metadata.raw.table<-tibble(Assay = Assays_used_task_oredered, Assay.type = args$param.assay.type.vec)

  Task_Comb_by_Assay_tib %<>% # vlookup
    mutate(assay_i_type = assay_i %>% map_chr(~metadata.raw.table %>% filter(Assay==.x) %>% pull(Assay.type)))




  ## Step 3: limma filtering (either gene or module level)

  #### Gene level:

  ### 3.1.Gene:  filter by eBayes (with adjustment to covaraites)
  if(args$parame.gene.or.module=='gene'){
    Task_Comb_by_Assay_tib_Gene<-Task_Comb_by_Assay_tib %>%
      # Filter_M_Limma_eBayes_gene_REV_function(Task_Comb_by_Assay_tib$DF.i[[1]], y.either, Covaraites_DF_const)
      # ORIG, long only via limma:  mutate(Assay_i_stat_score = DF_i_fixed_rowname %>% map(~Filter_M_Limma_eBayes_gene_REV_function(.x, y.either, Covaraites_DF_const)) ) # param.Univ.filt.top.n.features
      mutate(Assay_i_stat_score = map2(DF_i_fixed_rowname, assay_i_type, function(x, y){
        # x = Task_Comb_by_Assay_tib_Gene$DF_i_fixed_rowname[[10]]
        if(y == 'Long')  out<-Filter_M_Limma_eBayes_gene_REV_function(x, y.either, Covaraites_DF_const)  # param.Univ.filt.top.n.features
        if(y == 'Short') out<-Filter_M_glm_gene_REV_function(         x, y.either, Covaraites_DF_const) # param.Univ.filt.top.n.features
        out
      }))
    # Task_Comb_by_Assay_tib_Gene$Assay_i_stat_score[[1]]


    ## Step 3.2.Gene: combine selected features accross all assays:
    Features_Master_all_F_Gene<-Task_Comb_by_Assay_tib_Gene %>% select(assay_i, Assay_i_stat_score) %>% unnest()
    ## 'order.within.assay' "internal" features will be level-up, later will be used for filtering

    ## Step 3.3.Gene: (actual) Univ filtering:
    Features_Master_all_F_Gene %<>%
      mutate(assay_i = assay_i %>% as.character) %>%
      unite(assay_i_gene , c(assay_i, gene), sep = Assay.Analyte.sep, remove = FALSE) %>%
      mutate(passed_univ_fw_abs = order.within.assay %>% map_lgl(~.x >= -args$param.Univ.filt.top.n.features))
    # table(Features_Master_all_F_Gene$assay_i, Features_Master_all_F_Gene$passed_univ_fw_abs)

    Features_Master_Passed_Univ_F_Gene<-Features_Master_all_F_Gene %>% filter(passed_univ_fw_abs) %>% pull(assay_i_gene)

    dataOUT_Gene<-task.multi %>% subsetTask(features = Features_Master_Passed_Univ_F_Gene) %>% getTaskData(target.extra = TRUE) %>% .$data # target will be added outside, at PreProc function
    # dataOUT_Gene[1:3, 1:4]
    # rownames(dataOUT_Gene)
    out<-list(
      Features_Master_all = Features_Master_all_F_Gene,
      DF_passed = dataOUT_Gene # outcome and covaraites will be added at PreProcess_...() function
    )
  } # end gene






  ### Module level:

  ### 3.1.Module: filter by camera (with adjustment to covaraites)
  if(args$parame.gene.or.module=='module'){

    ## prep: ids2indices
    Task_Comb_by_Assay_tib_Module<-Task_Comb_by_Assay_tib %>%
      mutate(set_i_F1_Symbol_to_row_ideces = DF_i_fixed_rowname %>% map(~F_assay_i_module_row_index(.x, modules_sets_combined_l)))
    # Task_Comb_by_Assay_tib_Module$set_i_F1_Symbol_to_row_ideces[[1]]

    Task_Comb_by_Assay_tib_Module %<>%
      mutate(Assay_i_stat_score = pmap(list(DF_i_fixed_rowname, set_i_F1_Symbol_to_row_ideces, assay_i),
                                       possibly( function(DF_i, set_i_F1_Symbol_to_row_ideces, assay_i)
                                         # DF_i=Task_Comb_by_Assay_tib_Module$DF_i_fixed_rowname[[1]]; set_i_F1_Symbol_to_row_ideces=Task_Comb_by_Assay_tib_Module$set_i_F1_Symbol_to_row_ideces[[1]]; assay_i=Task_Comb_by_Assay_tib_Module$assay_i[[1]]
                                         Filter_M_Limma_camera_Module_REV_function(DF_i, y.either, Covaraites_DF_const, set_i_F1_Symbol_to_row_ideces, assay_i)
                                         , otherwise = NA_real_)
      ))

    ## Step 3.3.Module: (actual) Univ filtering:
    # Task_Comb_by_Assay_tib_Module$Assay_i_stat_score[[1]] %>% head
    Task_Comb_by_Assay_tib_Module  %<>%
      mutate(Assay_i_passed_camera_names = Assay_i_stat_score %>% map(function(z) z %>%
                                                                        mutate(passed.univ.internal = order.within.assay %>% map_lgl(function(y) y >= -args$param.Univ.filt.top.n.features) ) %>%
                                                                        filter(passed.univ.internal) %>%
                                                                        pull(gene)
      ))
    # Task_Comb_by_Assay_tib_Module$Assay_i_passed_camera_names[[1]] %>% head

    ## Step 3.2.Module: combine selected features accross all assays:
    Features_Master_all_F_Module<-Task_Comb_by_Assay_tib_Module %>% select(assay_i, Assay_i_stat_score) %>% unnest()

    ## Step 3.3.Module: Univ filtering:
    Features_Master_all_F_Module %<>%
      unite(assay_i_gene , c(assay_i, gene), sep = Assay.Analyte.sep, remove = FALSE) %>% ## Add assay namae to for module feature names: was not done internally before, to allow raw data retreival

      mutate(passed_univ_fw_abs = order.within.assay %>% map_lgl(~.x >= -args$param.Univ.filt.top.n.features))
    # Features_Master_all_F_Module$passed_univ_fw_abs %>% table

    ## part 2:
    ### Step 3.4.Module: Raw data retreival, for selected genes/modules:
    Task_Comb_by_Assay_tib_Module %<>%
      mutate(set_i_F1_Symbol_to_row_ideces_selected = map2(set_i_F1_Symbol_to_row_ideces, Assay_i_passed_camera_names,
                                                           ~.x[ (.x %>% names) %in% .y ])) %>%
      # Task_Comb_by_Assay_tib_Module$set_i_F1_Symbol_to_row_ideces_selected %>% map(names) %>% venn

      ## Create DF_module: Collapse/average set rows
      mutate(DF_i_selected_feature = map2(DF.i,
                                          set_i_F1_Symbol_to_row_ideces_selected,
                                          function(DF, All_sets) All_sets %>% map(function(set_i) DF[set_i,] %>% summarize_all(mean) %>% t) %>% data.frame)
      ) %>%
      ## Add assay name to feature/module names
      mutate(DF_i_selected_feature_Assay_names_pre=map2(DF_i_selected_feature, assay_i,
                                                        # .x=Task_Comb_by_Assay_tib_Module$DF_i_selected_feature[[1]]; .y=Task_Comb_by_Assay_tib_Module$assay_i[[1]]
                                                        ~.x %>% set_colnames(paste0(.y, Assay.Analyte.sep, .x %>% names))
      ))

    dataOUT_Module<-Task_Comb_by_Assay_tib_Module$DF_i_selected_feature_Assay_names_pre %>% data.frame
    # rownames(dataOUT_Module)

    ## Features_and_DF_L has two elements:
    ##  $Features_Master_all is used for feature ranking within clustering, as reference list
    ##  $DF_passed is used for the actual filtering
    out<-list(
      Features_Master_all = Features_Master_all_F_Module,
      DF_passed = dataOUT_Module
    )

  } # end module


  return(out)
}

# TRY1<-Func_task_to_Passed_Features_DF_MaG(task, args)
# TRY1$Features_Master_all %>% filter(assay_i=='ADCC') %>% select(c(4:6))
# TRY1$Features_Master_all %>% filter(order.within.assay>=-2) %>% select(c(1,4:6))
# TRY1$DF_passed %>% tail







#' @export
#' @rdname Wrapper_Filter_1_helpers

# Helper2: module Predict (preProc)
Func_Collapse_Data_by_module_names<-function(DF_gene, Train_feature_names){  # modules_sets_combined_l is also a parameter. for now will be a global 'floating' variable, but should be included in args!!
  # DF_gene = dataOUT_Predict; Train_feature_names = control$feature_names

  # neither has covaraites. filtered out before.
  if( !exists('modules_sets_combined_l')) cat('modules sets are not loaded')
  # modules_sets_combined_l %>% head

  modules_sets_combined_l_FIXED<-modules_sets_combined_l %>% setNames(modules_sets_combined_l %>% names %>% make.names())
  # modules_sets_combined_l_FIXED %>% head


  ## Step 1: Break coimbined Train_feature_names by assay -> nest
  Module_by_assay<-Train_feature_names %>% data.frame %>% set_colnames('Module_merged') %>%
    separate(col=1, into=c('assay_i','Module'), sep = Assay.Analyte.sep) %>%
    group_by(assay_i) %>% nest(.key = 'Module_i_df') %>%
    mutate(Module_i_names = Module_i_df %>% map(~.x %>% pull(Module) ))
  # Module_by_assay$Module_i_names

  ## prep: filter Module Master
  Module_by_assay %<>%
    mutate(Module_i_list = Module_i_names %>% map(~modules_sets_combined_l_FIXED[names(modules_sets_combined_l_FIXED) %in% . ]))
  # Module_by_assay$Module_i_list[[1]]



  ## Step 2: Break combined DF_gene by assay -> nest
  features_Gene_location<-colnames(DF_gene) %>% str_locate(Assay.Analyte.sep) %>% data.frame %>% pull(start)
  # table(features_Gene_location)
  Task_Comb_Assay_type_v<-colnames(DF_gene) %>% substr(1, features_Gene_location-1) ### !!! this include the covariates at the END!
  # table(Task_Comb_Assay_type_v)


  ## Step 1: within each assay, select features based on 'top' univariate score (by assay)
  DF_by_Assay_tib<-data.frame(assay_i = Task_Comb_Assay_type_v, DF_gene %>% t) %>% rownames_to_column('Analyte.name') %>%
    group_by(assay_i) %>% nest(.key = 'DF.and.Analyte') %>%
    filter(assay_i %in% Module_by_assay$assay_i) %>%
    mutate(assay_i = assay_i %>% as.character) %>%
    mutate(DF_i = DF.and.Analyte %>% map(~.x %>% data.frame %>% column_to_rownames('Analyte.name') )) %>%
    # .x=DF_by_Assay_tib$DF_i[[1]]
    mutate(DF_i_fixed_rowname = DF_i %>% map(~.x %>% set_rownames(Func_extract_feature_name(.x %>% rownames, Assay.Analyte.sep) )))
  # DF_by_Assay_tib$DF_i[[1]][1:3, 1:4] %>% rownames
  # DF_by_Assay_tib$DF_i_fixed_rowname[[1]][1:3, 1:4] %>% rownames



  ## merge sets with df, all by assay
  ## filtered before to match the same assays. to make sure same order remain, do full_join
  DF_by_Assay_tib<-full_join(DF_by_Assay_tib, Module_by_assay, by='assay_i')
  # DF_by_Assay_tib %<>% add_column(Module_i_list = Module_by_assay$Module_i_list) # at Predict: Error: `.data` must have 5 rows, not 2

  # ids2indices
  DF_by_Assay_tib %<>%
    # rownames(DF_by_Assay_tib$DF_i_fixed_rowname[[i]])
    # names(DF_by_Assay_tib$Module_i_list[[i]])
    # i=1; F_assay_i_module_row_index(DF_by_Assay_tib$DF_i_fixed_rowname[[i]], DF_by_Assay_tib$Module_i_list[[i]])
    mutate(set_i_F1_Symbol_to_row_ideces = map2(DF_i_fixed_rowname, Module_i_list, ~F_assay_i_module_row_index(.x, .y)))
  # DF_by_Assay_tib$set_i_F1_Symbol_to_row_ideces[[1]]


  ## Create DF_module: Collapse/average set rows via summarize_all(mean)
  ## alternative: limma::averep()
  DF_by_Assay_tib %<>%
    mutate(DF_i_selected_feature = map2(DF_i,
                                        set_i_F1_Symbol_to_row_ideces,
                                        function(DF, All_sets) All_sets %>% map(function(set_i) DF[set_i,] %>% summarize_all(mean) %>% t) %>% data.frame)
    )
  # DF_by_Assay_tib$DF_i_selected_feature[[1]]

  ## Add assay name to feature/module names (internally)
  DF_by_Assay_tib %<>%
    mutate(DF_i_selected_feature_Assay_names_pre=map2(DF_i_selected_feature, assay_i,
                                                      # .x=DF_Comb_by_Assay_tib_Module$DF_i_selected_feature[[1]]; .y=DF_Comb_by_Assay_tib_Module$assay_i[[1]]
                                                      ~.x %>% set_colnames(paste0(.y, Assay.Analyte.sep, .x %>% names))
    ))
  dataOUT_Module<-DF_by_Assay_tib$DF_i_selected_feature_Assay_names_pre %>% data.frame
  # names(dataOUT_Module)

  return(dataOUT_Module)
}

# TRY3<-Func_Collapse_Data_by_module_names(DF_gene, Train_feature_names)
# TRY3 %>% names
# TRY3[1:3, 1:3]







#' @export
#' @rdname Wrapper_Filter_1_helpers

#' @details when cv.glmnet lambda.min return a null model (intercept only), alternatively return 'previous'
#' value of lambda that reults \code {param.LASSO.n.features.arbitrary} features
#'
# Alternative nearest lambda
# https://github.com/mlr-org/mlr/issues/1030
cvglmnet_alternative_lambda.min_function<-function(cv_fit, param.LASSO.n.features.arbitrary){
  # cv_fit = model.lasso.cv.tune$learner.model; param.LASSO.n.features.arbitrary = 5
  # cv_fit = fitted_models_lrn[[1]]
  # cv_fit = R_AB_model_blocking$models %>% map(getLearnerModel) %>% map(getLearnerModel) %>% .[[1]]
  # cv_fit = cvob1
  # cv_fit = train(makeLearner(cl= "regr.cvglmnet", par.vals = list(alpha=1) ), task_baked_UnivClust)$learner.model # this will work on baked, AFTER pre-processing was done

  # TargetName<-task_j %>% getTaskTargetNames;
  # cv_fit = cv.glmnet(baked_UnivClust$data %>% select(-TargetName) %>% as.matrix, task_j %>% getTaskTargets, alpha = 1, nfolds = task_j %>% getTaskSize)
  # crit.lambda.min.no.features.ALTERNATIVE<-(cv_fit %>% tidy %>% filter(lambda==cv_fit$lambda.min) %>% pull(nzero))==0
  max_nzero<-cv_fit %>% tidy %>% pull(nzero) %>% max

  lambda.min<-cv_fit$lambda.min
  lambda.1se<-cv_fit$lambda.1se

  ## if no match, pick the next (higher) value
  param.LASSO.n.features.arbitrary_next<-tidy(cv_fit)$nzero %>% .[which.max(.>param.LASSO.n.features.arbitrary)]
  param.LASSO.n.features.arbitrary<-ifelse(param.LASSO.n.features.arbitrary %in% tidy(cv_fit)$nzero, param.LASSO.n.features.arbitrary, param.LASSO.n.features.arbitrary_next)

  pushed_lambda<-cv_fit %>% tidy %>%
    filter(nzero==min(max_nzero, param.LASSO.n.features.arbitrary)) %>% # if arbitrary value is too high and not in the range at all
    top_n(n=1, wt=estimate) %>%
    pull(lambda)

  out<-list(lambda.min    = lambda.min,
            lambda.1se    = lambda.1se,
            pushed_lambda = pushed_lambda)
  # crit_lambda_min_no_features<-( (cv_fit %>% coef(s = "lambda.min") %>% tidy %>% nrow) == 1) # intercept only
  # out<-ifelse( !crit_lambda_min_no_features, default_lambda, pushed_lambda)
  return(out)
}

#' @examples
#' set.seed(1011)
#' cvob1=cv.glmnet(x = matrix(rnorm(100*20),100,20), y = rnorm(100))
#' cvob1 %>% tidy
#' cvob1 %>% plot
#' cvob1$lambda.min
#' cv_fit_example<-cvob1
#' alt_lambda_1<-cvglmnet_alternative_lambda.min_function(cvob1, param.LASSO.n.features.arbitrary = 6)
#' coef(cvob1, s = alt_lambda_1) %>% tidy


#' set.seed(3)
#' cvob2=cv.glmnet(x = matrix(rnorm(100*3),100, 5), y = rnorm(100, 5))
#' cvob2 %>% tidy
#' cvob2 %>% plot
#' cvob2$lambda.min
#' cv_fit_example<-cvob2
#' alt_lambda_2<-cvglmnet_alternative_lambda.min_function(cvob2, param.LASSO.n.features.arbitrary = 5)
#' coef(cvob2, s = alt_lambda_2) %>% tidy

