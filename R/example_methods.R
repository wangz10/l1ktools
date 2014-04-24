# This script contains examples for reading .gctx files into R.
# In order for the script to work properly, make sure that the R working
# directory is the directory containing the script.

# source the io script
source("cmap/io.R")

# read the gctx file
ds = parse.gctx("../data/modzs_n272x978.gctx")

# inspect the matrix
print(ds@mat[1:10,1:10])

# inspect the row annotations
print(ds@rdesc[1:10,])

# inspect the column annotations
print(ds@cdesc[1:10,])
