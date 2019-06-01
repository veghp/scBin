\name{scBinarizer}
\alias{scBinarizer}
\title{Binarize count matrix using a cutoff}
\description{Replace counts greater than cutoff with 1,
  and counts less than or equal to cutoff with 0.}
\usage{scBinarizer(x, cutoff = 0)}
\arguments{
  \item{x}{count matrix. Rows: samples (cells). Columns: features (genes).}
  \item{cutoff}{numeric. The applied cutoff.}
  \item{logit}{boolean. Use False for UMI counts. If True, performs a library-size normalization, before applying the cutoff.}
}
\details{
}
\value{Binarized count matrix.
}
\references{
}
\author{Peter Vegh}
\note{
}

\seealso{
}
\examples{

%## The function is currently defined as
%function (x, cutoff = 0, logit = F) 
%{
    %if (any(x < 0)) {
        %stop("Negative values in count matrix.")
    %}
    %if(logit) {
      %x.row.sums <- rowSums(x)
      %x <- sweep(x, 1, x.row.sums, "/")
      %x <- x * 1000000 # cpm
    %}

    %x[x <= cutoff] <- 0
    %x[x > cutoff] <- 1 # order is very important

    %return(x)
  %}
}
