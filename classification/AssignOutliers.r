#!/usr/bin/env Rscript
library(dbscan)
library(class)
library(optparse)
source("load_data.r")
reorder_labels_by_count <- function(x) {
  labels <- x
  class_by_count <- as.integer(names(sort(table(x), decreasing=T)))
  sub_list <- lapply(class_by_count, function(count, d=labels) {return(which(d==count))})
  for(i in 1:length(sub_list)) {
    labels <- replace(labels, sub_list[[i]], i)
  }
  return(labels)
}

refine_labels_hdbscan <- function(hdbcluster, k=16) {
  labels <- hdbcluster
  # the biggest cluster is almost always the noise class (0 in hdbscan)
  topk <- as.integer(names(sort(table(hdbcluster),decreasing=T)))[2:(k+1)]
  labels[which(!labels %in% topk)]<-NA
  sub_list <- lapply(topk, function(count, d=labels) {return(which(d==count))})
  for(i in 1:length(sub_list)) {
    labels <- replace(labels, sub_list[[i]], i)
  }
  labels[is.na(labels)] <- length(topk)+1
  return(labels)
}

get_knn <- function(data, k, labels) {
  UNSET <- max(labels)
  train <- data[which(labels < UNSET),]
  test <- data[which(labels == UNSET),]
  tl <- labels[which(labels<UNSET)]
  cl_knn <- knn(train, test, tl, k, prob=T)
  labels_new <- labels
  labels_new[which(labels_new == UNSET)] <- cl_knn
  return(labels_new)
}

########## Main function ########## 
option_list <- list(
  make_option(c("-i", "--input"), type="character", dest="input",
    help="Input file name for Fisher score of outliers."),
  make_option(c("-c", "--core"), type="character", dest="core",
    help="Input file name for Fisher score of core monomers."),
  make_option(c("-r", "--ref"), type="character", dest="ref",
    help="Input file name for kNN clustering of core monomers."),
  make_option(c("-o", "--output"), type="character", dest="output",
    help="Output file name."),
  make_option(c("-p", "--pca-dimension"), type="integer", dest="pca_dim",
    default=100, help="Number of PCs to be used for clustering."),
  make_option(c("-k", "--knn"), type="integer", dest="knn", default=3,
    help="K for k-NN class membership assignment.")
  )
parser <- OptionParser(option_list=option_list, add_help_option=T,
  description="An R script for assigning subtype to outlier monomers using k-NN.",
  usage="Usage: %prog [options] -i <input file1> -c <input file2> -r <input file3> -o <output file>")
opt <- parse_args(parser, args=commandArgs(trailingOnly=T))
if (is.null(opt$input) | is.null(opt$output) | is.null(opt$core)
    | is.null(opt$ref)) {
  print_help(parser)
  quit()
}

data <- read.table(opt$core, comment.char="", stringsAsFactors=F, row.names=1)
dup <- duplicated(data)
pca <- get_pca_fisher(data[!dup,])
class_knn <- read.table(opt$ref, comment.char="", stringsAsFactors=F, row.names=1)

dothers <- read.table(opt$input, comment.char="", stringsAsFactors=F, row.names=1)
fisher_others <- get_fisher(dothers)
pca_others <- predict(pca, fisher_others)
class_others <- knn(pca$x[, 1:opt$pca_dim], pca_others[, 1:opt$pca_dim]
  , class_knn[rownames(data)[!dup],], opt$knn, prob=T)

write.table(class_others, opt$output, row.names=rownames(dothers), sep="\t"
  , quote=F, col.names=F)
