#!/usr/bin/env Rscript

.libPaths(c("C:/Users/ziemn/Documents/R/win-library/4.1",
            "C:/Program Files/R/R-4.1.2/library", .libPaths()))
print(.libPaths())
library(HiCcompare)
suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser()

# specify our desired options
# by default ArgumentParser will add an help option

parser$add_argument("--method", default="KRnorm",
                    help = "Normalization method [default \"%(default)s\"]")
parser$add_argument("--input",
                    help="Input file")
parser$add_argument("--output",
                    help="Output file")

# file_name <-  file.path("C:", "Users", "ziemn", "EwaTAD", "chr1_25kb.RAWobserved.txt" )
# file_name <-  file.path("C:", "Users", "ziemn", "EwaTAD", "test.txt" )
# out_file_name <- file.path("C:", "Users", "ziemn", "EwaTAD", "test1.csv" )
args <- parser$parse_args()
raw_data <- read.table(gzfile(args$input),head=FALSE)
names(raw_data) <- c('start1','start2','IF1')
mat1 = sparse2full(raw_data[, c('start1', 'start2', 'IF1')])


zeros1 = which(colSums(mat1) == 0)

if (length(zeros1) > 0) {
  cr.mat1 = mat1[-zeros1, -zeros1]
} else {
  cr.mat1 = mat1
}

sim1.kr = KRnorm(cr.mat1)

colnames(sim1.kr) = colnames(cr.mat1)

write.table(sim1.kr,file=args$output)

# sim1.kr = full2sparse(sim1.kr)




