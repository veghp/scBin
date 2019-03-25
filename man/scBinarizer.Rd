\name{scBinarizer}
\alias{scBinarizer}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Binarize count matrix using a cutoff}
\description{Replace counts greater than cutoff with 1,
  and counts less than or equal to cutoff with 0.}
\usage{scBinarizer(x, cutoff = 0)}
\arguments{
  \item{x}{count matrix. Rows: samples (cells). Columns: features (genes).}
  \item{cutoff}{numeric. The applied cutoff.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{Binarized count matrix.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{Peter Vegh}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
%##---- Should be DIRECTLY executable !! ----
%##-- ==>  Define data, use random,
%##--	or do  help(data=index)  for the standard data sets.

%## The function is currently defined as
%function (x, cutoff = 0) 
%{
    %if (any(x < 0)) {
        %stop("Negative values in count matrix.")
    %}
    %if (!(cutoff == round(cutoff))) {
        %stop("Cutoff is not a whole number.")
    %}
    %x[x > cutoff] <- 1
    %x[x <= cutoff] <- 0
    %return(x)
  %}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
