#' Multi Value Trait simulating function
#' 
#' @description growthSim can be used to help pick reasonable parameters for common
#'  growth models to use in prior distributions or to simulate data for example models/plots.
#' 
#' @param n_samples Number of samples per distribution to generate. Defaults to 10.
#' @param counts Number of counts per histogram, defaults to 1000.
#' @param max_bin The number of bins to return. Note that this is also the max value that will be
#' accepted in the distribution functions, with higher numbers being shrunk to this value. 
#' Defaults to 180.
#' @param dists A list of lists, with names corresponding to random deviate generating functions
#' and arguments to the function in the list values (see examples). Note that the n argument
#' does not need to be provided.
#' @param wide Boolean, should data be returned in wide format (the default)?
#' If FALSE then long data is returned.
#' @keywords multi-value
#' @return Returns a dataframe of example multi-value trait data simulated from specified distributions.
#' 
#' @importFrom graphics hist
#' @importFrom data.table melt as.data.table
#' 
#' @examples 
#' 
#' ## Not run:
#' library(extraDistr) # for rmixnorm
#' library(ggplot2)
#' n_samples = 10
#' counts = 1000
#' max_bin = 180
#' dists <- list(rmixnorm = list(mean=c(70, 150), sd = c(15,5), alpha = c(0.3, 0.7)),
#'               rnorm = list(mean = 90, sd = 3))
#' x <- mvSim(dists=dists, wide=FALSE)
#' dim(x)
#' x2 <- mvSim(dists=dists)
#' dim(x2)
#' 
#' ggplot(x, aes(x=as.numeric(sub("sim_", "", variable)),
#'               y = value, group = interaction(group,id), fill = group))+
#'   geom_col(position="identity", alpha=0.25)+
#'   pcv_theme()+
#'   labs(x="bin")
#' 
#' ## End(Not run)
#' @export

mvSim <- function(dists = list(rnorm = list(mean = 100, sd = 15)),
                  n_samples=10, counts=1000, max_bin = 180, wide = TRUE){
  vecs <- .makeVecs(dists, counts, n_samples)
  out <- .simFreqs(vecs, max_bin)
  if(!wide){
    out$id <- 1:nrow(out)
    out <- as.data.frame(data.table::melt(data.table::as.data.table(out),id.vars = c("group", "id")))
  }
  return(out)
}

#' internal vector making helper function
#' @keywords internal
#' @noRd

.makeVecs <- function(dists, counts, n_samples){
  funs <- names(dists)
  out <- lapply(funs, function(f){
    fun <- match.fun(f)
    dists[[f]]$n <- counts
    lapply(1:n_samples, function(i){do.call(fun, args = dists[[f]])})
  })
  names(out)<- names(dists)
  return(out)
}

#' internal histogram making helper function
#' @keywords internal
#' @noRd

.simFreqs<-function(vecs, max_bin){
  do.call(rbind, lapply(names(vecs), function(vecName){
    vec <- vecs[[vecName]]
    do.call(rbind,lapply(vec, function(v){
      v[v>max_bin] <- max_bin
      v[v<1] <- 1
      s1<-hist(v, breaks=seq(1,(max_bin+1),1), plot=FALSE)$counts
      s1d<-as.data.frame(cbind(data.frame(vecName), matrix(s1,nrow=1)))
      colnames(s1d)<-c("group", paste0("sim_",1:max_bin))
      s1d
    }))
  }))
}







