#' Read in lemnatech watering data from metadata.json files
#' 
#' @param file Path to a json file of lemnatech metadata.
#' @param envKey Character string representing the json key for environment data. By default this is set to "environment". Currently there are no situations where this makes sense to change.
#' @keywords watering, json
#' @import jsonlite
#' @return A data frame containing the bellwether watering data
#' @examples
#' wateringData<-bw.water("example.json") 
#' @export

bw.water<-function(file = NULL, envKey="environment"){
  meta <- jsonlite::fromJSON(txt = file)
  env<-as.data.frame(do.call(rbind, meta[[envKey]]))
  env$snapshot<-rownames(env)
  rownames(env)<-NULL
  env <- as.data.frame(apply(env, 2, as.character))
  env<-type.convert(env, as.is=T)
  if("timestamp" %in% colnames(env)){
    tryCatch({
      time <-as.POSIXct(env$timestamp)
      begin<-min(time,na.rm=T)
      DAS<-as.numeric((time-begin)/24/60/60)
      env$DAS<-DAS
    }, error = function(err){}, warning=function(warn){})
  }
  return(env)
}


