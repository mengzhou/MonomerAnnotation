# Classification scripts
This directory includes scripts used for the analysis of monomer subtype
classification based on Fisher score vectors.

## Dependencies
The scripts are written in R. They are tested in R version 3.4.4.
Below is a list of required R packages:
  * dbscan (version 1.1-2)
  * class (version 7.3-14)
  * optparse (version 1.6.0)

## Instructions
Brief instructions for the scripts are given here.

### FisherClustering.r
This script does clustering analysis for the input Fisher score vectors.
The input file should be generated by the `fisher` program, which can
be found in the `RepeatProfileHMM` package included in the detection
directory. Note that in the monomer classification study only the "core"
monomers were used for this script.
Below is a list of parameters for this script:
  * `-i`: specify the input file path of the Fisher score vectors of the core
    monomers.
  * `-o`: specity the output file name.
  * `-p`: set the number of PCs (principal components) to be used. Default: 100.
  * `-m`: set the minimum points parameter for HDBSCAN. Default: 20.
  * `-k`: set the K for k-NN subtype assignment. Default: 3.

The output of this script is a two-column matrix, with the first column
being the names of input sequences, and the second column being cluster
(subtype) IDs. The IDs are reordered based on the counts from high to low.

### AssignOutliers.r
This script assigns cluster IDs to monomers which were not included in the
previous clustering analysis. These monomers are considered as outliers due
to their lengths. Therefore their subtypes will only be inferred using k-NN.
Below is a list of parameters for this script:
  * `-i`: specify the input file path for the Fisher score vectors of the
    outlier monomers.
  * `-c`: specify the Fisher score vectors of the core monomers. This should be
    the input file of `FisherClustering.r`.
  * `-r`: specify the clustering of core monomers. This should be the output
    of `FisherClustering.r`.
  * `-o`: specity the output file for the clustering result of outlier monomers.
  * `-p`: set the number of PCs used for PCA. This should match what was used
    by `FisherClustering.r`.
  * `-k`: set the number of k for k-NN subtype assignment. Default: 3.

The output of this script is similar to that of `FisherClustering.r`. Together
these two scripts can be used for subtype classification of all detected
monomers.
