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


scLogit <- function(x, meta, markers, prob = 0.95) {

  if(!identical(rownames(x), rownames(meta))) {
    stop("Count matrix and metadata rownames are not identical.")}

  # only one marker, [[1]], is used in this version:
  protein <- meta[, markers[[1]][2]] > 0
  if(!all(protein == 0 | protein == 1)) {
    stop("Protein metadata column must contain only 0 and 1.")}

  gene <- markers[[1]][1] # multiple marker testing will be added later

  x.row.sums <- rowSums(x)
  cpm <- sweep(x, 1, x.row.sums, "/")
  cpm <- cpm * 1000000

  refprotein <- meta[, markers[[1]][2]]
  refgene <- cpm[, gene]

  df <- data.frame(gene = refgene, protein = refprotein)

  # Apply glm logit model
  glm.fit <- glm(protein ~ gene, family = binomial(link = "logit"), data = df)

  # Create data for plotting the curve:
  newdat <- data.frame(gene = seq(min(df$gene), max(df$gene), len = 100))
  newdat$protein <- predict(glm.fit, newdata = newdat, type = "response")

  cutoff <- (log(prob / (1 - prob)) - coef(glm.fit)[1])  /  (coef(glm.fit)[2])

  plot(
    x = df$gene,
    y = df$protein,
    main = gene,
    xlab = "Gene count",
    ylab = "Protein",
    las = 1,
    yaxt = "n")
  abline(v = cutoff, col = "grey")
  axis(2, at = 0:1, labels = 0:1, las = 1)
  mtext(round(cutoff, digits = 1), side = 3, at = cutoff, col = "red")
  lines(protein ~ gene, newdat, col = "red", lwd = 2)

  saveplot <- recordPlot()


  plotlist <- list(cutoff = cutoff, confusion = NULL, plots = saveplot)
  return(plotlist)
}
