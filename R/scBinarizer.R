scBinarizer <- function(x, cutoff = 0) {
  if(any(x < 0)) {stop("Negative values in count matrix.")}

  if(!(cutoff == round(cutoff))) {stop("Cutoff is not a whole number.")}

  x[x > cutoff] <- 1
  x[x <= cutoff] <- 0

  return(x)
}
