#' Run Partial Least Squares Regression on spectral data
#' 
#' @description Partial Least Squares Regression (plsr) is often used to analyze spectral data.
#' 
#' @param df Data frame containing metadata and spectral histogram data
#' @param resps Vector of response variables.
#' @param spectra Either one column name (in the case of long data) or a set of columns in the case of wide data.
#'  If a single character string is provided and it is not one of the column names then it is
#'  taken to be a pattern that will match some set of column names in the data to use (see examples).
#' @param train Proportion of data to use as training data. 
#' @param cv Number of cross validation iterations.
#' @param ... Further arguments passed to caret::train.
#' 
#' @details Note that columns that sum to 0 in the training or test data will be removed.
#' 
#' @import ggplot2
#' 
#' @keywords PLSR
#' @examples 
#' 
#' ## Not run:
#' 
#' file = "https://raw.githubusercontent.com/joshqsumner/pcvrTestData/main/pcvrTest1.csv"
#' df1<-read.pcv(file, "wide", T, multiValPattern = c("index_frequencies_index_ari",
#' "index_frequencies_index_ci_rededge", "npq_hist_NPQ", "yii_hist_Fq'/Fm'", "yii_hist_Fv/Fm"))
#' colnames(df1)<-sub("index_frequencies_index_ndvi.", "ndvi_", colnames(df1))
#' x<-pcv.plsr(df=df1, resps = "area.pixels", spectra = grepl("^ndvi_", colnames(df1)))
#' x<-pcv.plsr(df=df1, resps = "area.pixels", spectra = "^ndvi_")
#' 
#' ## End(Not run)
#' 
#' @export

pcv.plsr<-function(df, resps = NULL, spectra=NULL, train=0.8, cv=10, ...){
  if(is.character(spectra) && length(spectra)==1 && (!spectra %in% colnames(df))){
    spectra<-grepl(spectra, colnames(df))
  }
  if(is.numeric(resps)){resps<-colnames(df)[resps]}
  if(is.numeric(spectra)){resps<-colnames(df)[spectra]}
  outList<-lapply(resps, function(resp){
    sub <- cbind(data.frame(outcome=df[,resp]), df[,spectra])
    sub <- sub[complete.cases(sub), ]
    #sub <- sub[,colSums(sub)>0]
    #* Training Test Split and removing full 0 columns
    trainIndex<-sample(c(1:nrow(sub)), size=ceiling(train*nrow(sub)))
    train.sub<-sub[trainIndex,]
    train.non.0<-colSums(train.sub[,2:ncol(train.sub)])>0
    test.sub<-sub[-trainIndex,]
    test.non.0<-colSums(test.sub[,2:ncol(test.sub)])>0
    train.sub<-train.sub[, c(T, train.non.0 & test.non.0)]
    test.sub<-test.sub[, c(T, train.non.0 & test.non.0)]
    #* Fit model
    model <- caret::train(
      outcome ~., data = train.sub, method = "pls",
      scale = TRUE,
      trControl = caret::trainControl("cv", number = cv),
      tuneLength = 10)
    coef_df<-as.data.frame(coef(model$finalModel, ncomp =1:model$finalModel$ncomp))
    colnames(coef_df)<-paste0("coef_",1:model$finalModel$ncomp)
    coef_df$wavelength=rownames(coef_df)
    if(model$finalModel$ncomp==1){
      vip_df <- as.data.frame(vip(model$finalModel))
    } else {vip_df <- as.data.frame(t(vip(model$finalModel)))}
    colnames(vip_df)<-paste0("VIP_",1:model$finalModel$ncomp)
    vip_df<-cbind(vip_df,coef_df)
    rownames(vip_df)<-c(1:nrow(vip_df))
    # Make predictions
    predictions <- predict(model, test.sub)
    # Model performance metrics
    performance <- data.frame(
      RMSE = caret::RMSE(predictions, test.sub$outcome),
      Rsquare = caret::R2(predictions, test.sub$outcome),
      MAE = caret::MAE(predictions, test.sub$outcome),
      outcome = resp)
    plot <- ggplot2::ggplot(model)+labs(title=paste0(resp))
    list("model_performance" = performance,
                         "predictionTarget"=resp,
                         "model" = model,
                         "plot" = plot,
                         "number_of_components" = model$finalModel$ncomp,
                         "vip_df"=vip_df)
  })
  names(outList)<-resps
  return(outList)
}

vip <- function(object) {
  if (object$method != "oscorespls"){
    stop("Only implemented for orthogonal scores algorithm.  Refit with \"method = 'oscorespls'\"")}
  if (nrow(object$Yloadings) > 1){stop("VIP can only be calculated for a model with one response term")}
  SS <- c(object$Yloadings)^2 * colSums(object$scores^2)
  Wnorm2 <- colSums(object$loading.weights^2)
  SSW <- sweep(object$loading.weights^2, 2, SS / Wnorm2, "*")
  sqrt(nrow(SSW) * apply(SSW, 1, cumsum) / cumsum(SS))
}
