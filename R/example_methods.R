# This script contains examples for reading .gctx files into R.
# In order for the script to work properly, make sure that the R working
# directory is the directory containing the script.

# source the io script
source("cmap/io.R")


# initialize a GCT object with a single signature
ds <- parse.gctx("/xchip/cogs/data/build/a2y13q1/modzs.gctx", cid="CPC005_A375_6H:BRD-A85280935-003-01-7:10")


# make a call to the CLUE API
require(RCurl)
require(rjson)

# count the number of sirolimus signatures in the MCF7 cell line
url <- 'https://api.clue.io/a2/siginfo?q={"pert_desc":"sirolimus","cell_id":"MCF7"}&user_key=lincsdemo&c=1'
response <- fromJSON(getURL(url, .opts = list(ssl.verifypeer = FALSE)))

# The '.opts = list(ssl.verifypeer = FALSE)' part disables the SSL certificate check.
