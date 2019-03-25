scPlot <- function(x, meta, markers, tnr = 0.95, jitter = F) {

  if(any(x < 0)) {stop("Negative values in count matrix.")}

  if(!identical(rownames(x), rownames(meta))) {
    stop("Count matrix and metadata rownames are not identical.")}

  # Will implement a loop on markers later.
  protein <- meta[, markers[[1]][2]] > 0
  if(!all(protein == 0 | protein == 1)) {
    stop("Protein metadata column must contain only 0 and 1.")}

  gene <- markers[[1]][1]
  gene.counts <- x[, gene]
  counts <- gene.counts[!protein] # need the protein-negative
  counts.nonzero <- counts[counts > 0] # for calculating cutoff

  # (TNR)th percentile
  cutoff <- quantile(counts.nonzero, probs = tnr, type = 1)
  
  plot(
    x = gene.counts,
    y = (if(jitter == T) {
      jitter(meta[, markers[[1]][2]], amount = 0.1)
      } else {meta[, markers[[1]][2]]}),
    main = gene,
    xlab = "Gene count",
    ylab = "Protein",
    las = 1,
    yaxt = "n")
  abline(v = cutoff, col = "red")
  axis(2, at = 0:1, labels = 0:1, las = 1)
  mtext(cutoff, side = 3, at = cutoff, col = "red")
  saveplot <- recordPlot()

  conf <- as.matrix(table(factor(as.integer(gene.counts > cutoff),
    levels=c(0, 1)), protein))


  plotlist <- list(cutoff = cutoff, confusion = conf, plots = saveplot)
  return(plotlist)
}
