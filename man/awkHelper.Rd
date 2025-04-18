% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/awkHelper.R
\name{awkHelper}
\alias{awkHelper}
\title{subset helper function for use reading in large data, called in pcv.sub.read}
\usage{
awkHelper(inputFile, filters, awk = NULL)
}
\arguments{
\item{inputFile}{Path to csv file of plantCV output, should be provided internally in read.pcv}

\item{filters}{filtering conditions, see read.pcv for details. Format as
list("trait in area, perimeter", "other contains stringToMatch")}

\item{awk}{Optional awk command to use instead.}
}
\value{
Returns a character string representing a unix style awk statement
  which is typically passed to \code{pipe} or used as a connection in \code{data.table::fread}.
}
\description{
subset helper function for use reading in large data, called in pcv.sub.read
}
\details{
awkHelper attempts to make awk commands from human readable input.
Currently when filters are supplied the input file has quotes removed by `sed`
then is piped into awk, so an equivalent command line statement may be:
\code{sed 's/\"//g' pcvrTest2.csv | awk -F ','  '{ if (NR==1 || $18=="area") { print } }'}
}
\examples{
tryCatch(
  { # in case offline
    link1 <- "https://gist.githubusercontent.com/seankross/"
    link2 <- "a412dfbd88b3db70b74b/raw/5f23f993cd87c283ce766e7ac6b329ee7cc2e1d1/mtcars.csv"
    file <- paste0(link1, link2)
    awkHelper(file, list("gear in 4, 3"), awk = NULL)
    awkHelper(file, "gear contains 3", awk = NULL)
    # note that to be filtered the file has to exist on your local system, this example only shows
    # the output of awkHelper, which would then be executed by read.pcv on a unix system
    awkHelper(file, list("gear in 4, 3"), awk = "existing_command")
  },
  error = function(e) {
    message(e)
  }
)
}
\keyword{long}
\keyword{pcv}
\keyword{read.csv}
\keyword{wide}
