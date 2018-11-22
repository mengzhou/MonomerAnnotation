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
    help="Input file name for Fisher score of core monomers."),
  make_option(c("-o", "--output"), type="character", dest="output",
    help="Output file name."),
  make_option(c("-p", "--pca-dimension"), type="integer", dest="pca_dim",
    default=100, help="Number of PCs to be used for clustering."),
  make_option(c("-m", "--min-points"), type="integer", dest="min_pts",
    default=20, help="MinPoints parameter for HDBSCAN."),
  make_option(c("-k", "--knn"), type="integer", dest="knn", default=3,
    help="K for k-NN class membership assignment.")
  )
parser <- OptionParser(option_list=option_list, add_help_option=T,
  description="An R script for monomer subtype classification using Fisher score vectors.",
  usage="Usage: %prog [options] -i <input file> -o <output file>")
opt <- parse_args(parser, args=commandArgs(trailingOnly=T))
if (is.null(opt$input) | is.null(opt$output)) {
  print_help(parser)
  quit()
}

# getting clusters by HDBSCAN and knn
data <- read.table(opt$input, comment.char="", stringsAsFactors=F, row.names=1)
dup <- duplicated(data)
pca <- get_pca_fisher(data[!dup,])
hdb <- hdbscan(pca$x[, 1:opt$pca_dim], minPts=opt$min_pts)
N <- max(hdb$cluster)
class_hdb <- refine_labels_hdbscan(hdb$cluster, N-1)
class_knn <- get_knn(pca$x[,1:opt$pca_dim], opt$knn, class_hdb)
# add back all duplicated data points
dup_span <- which(sapply(1:(length(dup)-1), function(x){return(xor(dup[x], dup[x+1]))}))
starts <- dup_span[seq(1, length(dup_span), 2)]
ends <- dup_span[seq(2, length(dup_span), 2)]
if(length(ends)<length(starts)) { ends <- c(ends, length(dup)) }
class_all <- rep(0, nrow(data))
names(class_all) <- rownames(data)
for(i in 1:N){class_all[rownames(pca$x)[which(class_knn == i)]] = i}
for(i in 1:length(starts)) {
  for(j in seq(starts[i]+1, ends[i])) {
    class_all[j] <- class_all[starts[i]]
  }
}
class_all <- reorder_labels_by_count(class_all)

write.table(class_all, opt$output, row.names=names(class_all), sep="\t", quote=F, col.names=F)
