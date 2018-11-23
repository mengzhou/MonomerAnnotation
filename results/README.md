# Results
Results of monomer detection and classification in the mm10 genome

# Monomer detection and classification results
The file `monomers.tar.gz` includes three BED files, corresponding to three
monomer types, A, G, and T, respectively. Below is an example of the results:
```
chr1    3112101 3112297 A_mm10|chr1:3111602-3112793|-|%2#1      18      -
chr1    3112297 3112411 A_mm10|chr1:3111602-3112793|-|%1#2      2       -
```
The 4th column indicates the type of monomer, as well as the coordinates of
the original region for posterior decoding. The last two numbers in this
column, mark the position order of each monomer in a promoter. For example,
`%2#1` means this monomer is positioned at the 2nd position from the 5'-end
of the promoter, and at the 1st position from the 3'-end of the promoter.
The 5th column indicates the subtype ID of each monomer.

# Files for profile-HMM parameters
The file `profileHMM.tar.gz` includes the parameter files for all three
types of monomers. These files are trained after the iterative process.
