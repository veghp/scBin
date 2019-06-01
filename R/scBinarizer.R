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


scBinarizer <- function(x, cutoff = 0, logit = F) {
  if(any(x < 0)) {stop("Negative values in count matrix.")}

  if(logit) {
    x.row.sums <- rowSums(x)
    x <- sweep(x, 1, x.row.sums, "/")
    x <- x * 1000000 # cpm
  }

  x[x <= cutoff] <- 0
  x[x > cutoff] <- 1 # order is very important

  return(x)
}
