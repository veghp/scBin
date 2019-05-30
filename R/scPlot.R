# Copyright 2019 Peter Vegh
# 
# This file is part of scBin.
# 
# scBin is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# scBin is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with scBin.  If not, see <https://www.gnu.org/licenses/>.


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
