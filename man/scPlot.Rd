\name{scPlot}
\alias{scPlot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot protein vs gene expression
%%  ~~function to do ... ~~
}
\description{
Plots protein expression for specified marker genes and calculates a count cutoff for expression, a confusion matrix, and generates a plot.
}
\usage{
scPlot(x, meta, markers, tnr = 0.95, jitter = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{count matrix. Rows: samples (cells). Columns: features (genes).}
  \item{meta}{metadata dataframe, where one of the columns is a binary vector of flow (protein) sorting bins. Count and metadata cells must be the same and in the same order.}
  \item{markers}{list of vectors of length 2. The first element of each vector is a marker gene column name in 'x', the second element is the matching marker protein name in 'meta'. This column of 'meta' must contain only 0 and 1.}
  \item{tnr}{desired true negative rate of nonzero counts.}
  \item{jitter}{jitters the dots on the plot if set TRUE.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{List of three:
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{cutoff }{the calculated cutoff for given specificity for nonzero counts}
  \item{confusion }{confusion matrix of all values after applying the cutoff}
  \item{plots }{the saved plots of protein expression vs gene count}
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
%function (x, meta, markers, tnr = 0.95, jitter = F) 
%{
    %if (any(x < 0)) {
        %stop("Negative values in count matrix.")
    %}
    %if (!identical(rownames(x), rownames(meta))) {
        %stop("Count matrix and metadata rownames are not identical.")
    %}
    %protein <- meta[, markers[[1]][2]] > 0
    %if (!all(protein == 0 | protein == 1)) {
        %stop("Protein metadata column must contain only 0 and 1.")
    %}
    %gene <- markers[[1]][1]
    %gene.counts <- x[, gene]
    %counts <- gene.counts[!protein]
    %counts.nonzero <- counts[counts > 0]
    %cutoff <- quantile(counts.nonzero, probs = tnr, type = 1)
    %plot(x = gene.counts, y = (if (jitter == T) {
        %jitter(meta[, markers[[1]][2]], amount = 0.1)
    %}
    %else {
        %meta[, markers[[1]][2]]
    %}), main = gene, xlab = "Gene count", ylab = "Protein", las = 1, 
        %yaxt = "n")
    %abline(v = cutoff, col = "red")
    %axis(2, at = 0:1, labels = 0:1, las = 1)
    %mtext(cutoff, side = 3, at = cutoff, col = "red")
    %saveplot <- recordPlot()
    %conf <- as.matrix(table(factor(as.integer(gene.counts > cutoff), 
        %levels = c(0, 1)), protein))
    %plotlist <- list(cutoff = cutoff, confusion = conf, plots = saveplot)
    %return(plotlist)
  %}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
