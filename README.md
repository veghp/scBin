# scBin

scBin (**s**ingle-**c**ell RNA-seq data **Bin**arizer) plots a chosen gene's counts and determines a cutoff that is used for binarizing the UMI count data.

Please note that this package is **under development.**

## Background

Investigations often require that we determine presence or absence of expression of a given or all genes in a sample, as opposed to differentially expressed genes in relation to another sample or condition. This context-free perspective of gene expression requires a tool that classifies each gene's expression in a single cell or a cell type as 'not expressed' (0) or 'expressed' (1). The gene used for determining a cutoff is one for which we have protein expression information (e.g. a flow cytometry marker). This solution therefore requires that we have sequenced at least one positive and one negative population for that gene, and also works under the somewhat strong assumption that mRNA numbers for all other genes will have similar effect as for the reference gene. A key factor in the efficiency of the method is the accuracy of the protein annotation, i.e. the positive and negative (flow sort) populations should be well-separated. Also important to note that the method does not demonstrate the _absence_ of gene expression.

## Method

We know true expression of the flow marker gene for each cell, thus we can build a confusion matrix (contingency table) and calculate a cutoff for given _true negative rate_ (TNR, specificity/SPC; 95% default). The cutoff is calculated for _nonzero_ samples. The confusion matrix is calculated from _all_ samples, after applying the cutoff. As only the nonzero samples were used for calculating specificity, the specificity in the confusion matrix is much higher than the set specificity value. This ensures that we have high confidence in our expression calls.

A _logit method_ will be added later to the package, for non-UMI data. This will build a logit model from the flow marker gene count numbers and protein expression, and uses that for classifying gene expression. However, the low count numbers and the presence of many dropouts in droplet methods make it hard to build a useful model for UMI data.

## Example

The example data has 21623 cells (and 23310 genes) isolated from skin dermis. CD49f (ITGA6) was used to separate the endothelial cells (I) from the fibroblasts (H), and this will be used to classify the rest of the genes.

![CD49f](flow.jpg)

Load data and metadata:
```R
library("scBin")
library("Matrix")
x <- readMM("counts.txt")
colnames(x) <- scan("counts.colnames.txt", what = character())
rownames(x) <- scan("counts.rownames.txt", what = character())
x <- t(x)
meta <- read.delim("metadata.txt", stringsAsFactors = F, header = T, row.names = 1)
```

The metadata dataframe (meta), must have a binary vector column of flow (protein) sorting bins, used for the classification.

```R
meta$itga6 <- meta$Flow == "CD49f_positive"
meta$itga6 <- as.integer(meta$itga6)
markers <- list(c("ITGA6", "itga6"))
```

Data rows: cells, columns: genes.

```R
scplot.x <- scPlot(x, meta, markers, tnr = 0.95, jitter = T)
scplot.x$confusion
#   protein
#    FALSE  TRUE
#  0  8352 12980
#  1     1   290
scplot.x$plots
```

![CD49f](itga6.png)

The calculated cutoff ('3') is shown in red. Binarize data:

```R
bin.matrix <- scBinarizer(x, cutoff = scplot.x$cutoff)
```

The binarized data can be used either for a new analysis of cells; or for augmenting cell annotations to decide which genes are expressed in each cell type.
