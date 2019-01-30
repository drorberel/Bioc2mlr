#' Customized preprocessing wrapper. wrappers around makePreprocWrapper: Univariate only; Univariate + Cluster
#'
#' @param learner A learner
#' @param train_F train_F
#' @param Predict_F Predict_F
#' @param param.Univ.filt.top.n.features param.Univ.filt.top.n.features
#' @param parame.gene.or.module parame.gene.or.module
#' @param param.LASSO.n.features.arbitrary param.LASSO.n.features.arbitrary
#'
#' @param param.UnivClustRankTopN param.UnivClustRankTopN
#' @param param.cluster_method_KH param.cluster_method_KH
#' @param param.corrplot.n.clusters.k param.corrplot.n.clusters.k
#' @param param.corrplot.n.clusters.h param.corrplot.n.clusters.h
#' @param parame.gene.or.module parame.gene.or.module
#' @param param.LASSO.n.features.arbitrary param.LASSO.n.features.arbitrary
#' @param param.assay.type.vec param.assay.type.vec
#'
#' @name Wrapper_Filter_3_makePrep_MaG
NULL



#' @return \code{learner} custom preprocessing wrapper with makePreprocWrapper() of class 'PreprocWrapper'
#' @details TBA
#' @references https://pat-s.github.io/mlr/articles/tutorial/devel/preproc.html#creating-the-preprocessing-wrapper
#' @examples
#' My_learner<-lrn
#' F_PreProc_1_Univ_Only_Train_MaG<-A
#' F_PreProc_3_UnivClust_Train_MaG<-B
#' Fun_lrn_univ_only_makePrep_MaG(        learner = My_learner, train_F = F_PreProc_1_Univ_Only_Train_MaG, Predict_F = F_PreProc13_BOTH_Predict_MaG, param.Univ.filt.top.n.features = 10, parame.gene.or.module = 'gene', param.LASSO.n.features.arbitrary = '5')
#' Fun_lrn_univ_Clusters_All_makePrep_MaG(learner = My_learner, train_F = F_PreProc_3_UnivClust_Train_MaG, Predict_F = F_PreProc13_BOTH_Predict_MaG, param.Univ.filt.top.n.features, param.UnivClustRankTopN, param.cluster_method_KH, param.corrplot.n.clusters.k, param.corrplot.n.clusters.h, parame.gene.or.module, param.LASSO.n.features.arbitrary, param.assay.type.vec)

#' @export
#' @rdname Wrapper_Filter_3_makePrep_MaG
# PreProc1, Univ Only
# Univ only, gene + module
Fun_lrn_univ_only_makePrep_MaG<-function(learner, train_F, Predict_F, param.Univ.filt.top.n.features, parame.gene.or.module, param.LASSO.n.features.arbitrary){
  # learner = lrn.glm.1.orig, train_F = F_PreProc_1_Univ_Only_Train_MaG; Predict_F = F_PreProc13_BOTH_Predict_MaG
  makePreprocWrapper(
    learner = learner,
    train   = train_F,
    predict = Predict_F,
    par.set = makeParamSet( # just declare them. Default values does not matter. will be over-ride at next step by par.vals
      makeDiscreteParam("param.Univ.filt.top.n.features",   values=c(1, 2, 3)),
      makeDiscreteParam("parame.gene.or.module",            values=c('gene', 'module')),
      makeDiscreteParam("param.LASSO.n.features.arbitrary", values=c(1,2,3))
    ),
    par.vals = list(
      param.Univ.filt.top.n.features   = param.Univ.filt.top.n.features,
      parame.gene.or.module            = parame.gene.or.module,
      param.LASSO.n.features.arbitrary = param.LASSO.n.features.arbitrary
    ) # once declared, set default values for them
  )
}



#' @export
#' @rdname Wrapper_Filter_3_makePrep_MaG

# PreProc3, Univ+Clust
Fun_lrn_univ_Clusters_All_makePrep_MaG<-function(learner, train_F, Predict_F, param.Univ.filt.top.n.features, param.UnivClustRankTopN, param.cluster_method_KH, param.corrplot.n.clusters.k, param.corrplot.n.clusters.h, parame.gene.or.module, param.LASSO.n.features.arbitrary, param.assay.type.vec) {
  # learner = lrn.glmnet.1.orig;   train_F = F_PreProc_3_UnivClust_Train_MaG, predict_F = F_PreProc13_BOTH_Predict_MaG

  makePreprocWrapper(
    learner = learner,
    train   = train_F,
    predict = Predict_F,

    ## these parameters are declared as ParamSet, so they could be iteratively benchmarked
    par.set = makeParamSet( # just declare them. Default values does not matter. will be over-ride at next step by par.vals
      makeDiscreteParam("param.Univ.filt.top.n.features", values=c(1, 2, 3)),
      makeDiscreteParam("param.UnivClustRankTopN",        values=c(1, 2, 3)),
      makeDiscreteParam("param.cluster_method_KH",        values=c('method.k', 'method.h')),
      makeDiscreteParam("param.corrplot.n.clusters.k",    values=c(10, 20, 30)),
      makeDiscreteParam("param.corrplot.n.clusters.h",    values=c(0.3, 0.4, 1)),
      makeDiscreteParam("parame.gene.or.module",          values=c('gene', 'module')),
      makeDiscreteParam("param.LASSO.n.features.arbitrary", values=c(1,2,3))
    ),

    par.vals = list(       # once declared, set default values for them, BUT, these will be over-written by tuneParams()+GridSearch()
      param.Univ.filt.top.n.features   = param.Univ.filt.top.n.features,
      param.UnivClustRankTopN          = param.UnivClustRankTopN,
      param.cluster_method_KH          = param.cluster_method_KH,
      param.corrplot.n.clusters.k      = param.corrplot.n.clusters.k,
      param.corrplot.n.clusters.h      = param.corrplot.n.clusters.h,
      parame.gene.or.module            = parame.gene.or.module,
      param.LASSO.n.features.arbitrary = param.LASSO.n.features.arbitrary
    )

  ) # end makePreprocWrapper
}






# x=matrix(rnorm(100*20),100,20); y=rnorm(100); cv.fit=glmnet(x,y, alpha=1)
# plot(cv.fit);
# cv.fit$lambda
# cv.fit %>% coef %>% tidy %>% head
# cv.fit %>%  tidy %>% head
# predict(cv.fit, x)[,1]

# param.LASSO.n.features.arbitrary<-4

# task.example<-makeRegrTask(id='example',data.frame(x,y), target='y')
# model<-train(makeLearner("regr.glmnet", par.vals = list(alpha=1, s=0.209)), task.example)
# model$learner.model %>% coef %>% tidy %>% head
# model$learner.model %>% plot
# predict(model, task.example)





