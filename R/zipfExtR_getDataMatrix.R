#' Convert a sample vector to a frequency matrix.
#'
#' Converts a sequence of values into a matrix of frequencies.
#'
#' @param values Vector of positive integer values.
#'
#' @return The matrix of frequencies associated to the vector \code{values}.
#'
#' @examples
#' data <- rmoezipf(100, 2.5, 1.3)
#' zipfExtR_getDataMatrix(data)
#'
#' @export
zipfExtR_getDataMatrix <- function(values){
  spectrumTable <- table(values)
  data_frame <- data.frame(spectrumTable)
  frequencies <- as.numeric(as.matrix(data_frame[,2]))
  values <- as.numeric(as.matrix(data_frame[,1]))
  data <- matrix(ncol = 2, nrow = length(values))
  data[,1] <- values
  data[,2] <- frequencies
  return(data)
}
