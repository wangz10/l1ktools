# load dependencies
library(rhdf5)


########################################
### GCT class and method definitions ###
########################################

setClass("GCT",
         representation(
             mat = "matrix",
             rid = "vector",
             cid = "vector",
             rdesc = "data.frame",
             cdesc = "data.frame",
             version = "character",
             src = "character"
         )
)


#### define some helper methods for parsing gctx files ###

# helper function to set all the row and column annotations to the correct data type
fix.datatypes <- function(meta) {
    # turn all warnings to errors so we can use the try statement to grab strings
    options(warn = 2)
    for (field.name in names(meta)) {
        # get the field
        field = meta[[field.name]]
        # check if it's numeric
        try({field = as.numeric(field)}, silent = TRUE)
        if (is.numeric(field)) {
            #check if it's an integer
            int.field = NULL
            try({int.field = as.integer(field)}, silent = TRUE)
            if ( ! is.null(int.field) && identical(int.field, field) )
                field = int.field
        }
        # insert back into the annotations
        meta[[field.name]] = field
    }
    options(warn = 0)
    return(meta)
}


# helper function for parsing row or column metadata
read.gctx.meta <- function(gctx_path, dimension="row", ids=NULL) {
  if (!(dimension %in% c("row", "col"))) {
    stop("dimension can be either row or col")
  }
  if (dimension == "row") {
    name <- "0/META/ROW"
  } else {
    name <- "0/META/COL"
  }
  raw_annots <- h5read(gctx_path, name=name) # returns a list
  fields <- names(raw_annots)
  # define an empty data frame of the correct dimensions
  annots <-  data.frame(matrix(nrow=length(raw_annots[[fields[1]]]), ncol=length(fields)))
  colnames(annots) <-  fields
  # loop through each field and fill the annots data.frame
  for (i in 1:length(fields)) {
    field <- fields[i]
    # remove any trailing spaces
    annots[,i] <- gsub("\\s*$", "", raw_annots[[field]], perl=T)
  } 
  annots <- fix.datatypes(annots)
  # use the id field to set the rownames
  rownames(annots) <- as.character(annots$id)
  # subset to the provided set of ids, if given
  if (is.null(ids)) {
    ids <- rownames(annots)
  } else {
    ids <- ids
  }
  annots <- droplevels(subset(annots, id %in% ids))
  return(annots)
}


# helper function for reading gctx row/col ids
read.gctx.ids <- function(gctx_path, dimension="row") {
  if (!(dimension %in% c("row", "col"))) {
    stop("dimension can be either row or col")
  }
  if (dimension == "row") {
    name <- "0/META/ROW/id"
  } else {
    name <- "0/META/COL/id"
  }
  # remove any spaces
  ids <- gsub("\\s*$", "", h5read(gctx_path, name=name), perl=T)
  return(ids)
}


# define the initialization method for the class
setMethod("initialize",
          signature = "GCT",
          definition = function(.Object, src, rid = NULL, cid = NULL) {
              # check to make sure it's either .gct or .gctx
              if (! (grepl(".gct$", src) || grepl(".gctx$", src) ))
                  stop("Either a .gct or .gctx file must be given")
              if (grepl(".gct$", src)) {
                  if ( ! is.null(rid) || !is.null(cid) )
                      stop("rid and cid values may only be given for .gctx files, not .gct files")
                  # parse the .gct
                  .Object@src = src
                  # get the .gct version by reading first line
                  .Object@version = scan(src, what = "", nlines = 1, sep = "\t", quiet = TRUE)[1]
                  # get matrix dimensions by reading second line
                  dimensions = scan(src, what = double(0), nlines = 1, skip = 1, sep = "\t", quiet = TRUE)
                  nrmat = dimensions[1]
                  ncmat = dimensions[2]
                  nrhd = dimensions[3]
                  nchd = dimensions[4]
                  # read in header line
                  header = scan(src, what = "", nlines = 1, skip = 2, sep = "\t", quote = NULL, quiet = TRUE)
                  # construct row header and column id's from the header line
                  if ( nrhd ) {
                      rhd = header[2:(nrhd+1)]
                      cid = header[-(nrhd+1):-1]
                  }
                  else {
                      rhd = NULL
                      cid = header[-1]
                  }
                  # read in the next set of headers (column annotations) and shape into a matrix
                  if ( nchd ) {
                      header = scan(src, what = "", nlines = nchd, skip = 3, sep = "\t", 
                                    quote = NULL, quiet = TRUE)		
                      header = matrix(header, nrow = nchd, 
                                      ncol = ncmat + nrhd + 1, byrow = TRUE)
                      # extract the column header and column descriptions
                      chd = header[,1]
                      cdesc = header[,-(nrhd+1):-1]
                      # need to transpose in the case where there's only one column annotation
                      if ( nchd == 1 )
                          cdesc = t(cdesc)
                  }
                  else {
                      chd = NULL
                      cdesc = data.frame()
                  }
                  # read in the data matrix and row descriptions, shape into a matrix
                  mat = scan(src, what = "", nlines = nrmat, 
                             skip = 3 + nchd, sep = "\t", quote = NULL, quiet = TRUE)
                  mat = matrix(mat, nrow = nrmat, ncol = ncmat + nrhd + 1, 
                               byrow = TRUE)
                  # Extract the row id's row descriptions, and the data matrix
                  rid = mat[,1]
                  if ( nrhd ) {
                      # need as.matrix for the case where there's only one row annotation
                      rdesc = as.matrix(mat[,2:(nrhd + 1)])
                      mat = matrix(as.numeric(mat[,-(nrhd + 1):-1]),
                                   nrow = nrmat, ncol = ncmat)
                  }
                  else {
                      rdesc = data.frame()
                      mat = matrix(as.numeric(mat[,-1]),
                                   nrow = nrmat, ncol = ncmat)
                  }
                  # assign names to the data matrix and the row and column descriptions
                  dimnames(mat) = list(rid, cid)
                  if ( nrhd ) {
                      dimnames(rdesc) = list(rid,rhd)
                      rdesc = as.data.frame(rdesc, stringsAsFactors = FALSE)
                  }
                  if ( nchd ) {
                      cdesc = t(cdesc)
                      dimnames(cdesc) = list(cid,chd)
                      cdesc = as.data.frame(cdesc, stringsAsFactors = FALSE)
                  }
                  # assign to the GCT slots
                  .Object@mat = mat
                  .Object@rid = rownames(mat)
                  .Object@cid = colnames(mat)
                  .Object@rdesc = fix.datatypes(rdesc)
                  .Object@cdesc = fix.datatypes(cdesc)
                  return(.Object)
              }
              else { 
                  # parse the .gctx
                  .Object@src = src
                  # if the rid's or column id's are .grp files, read them in
                  if ( length(rid) == 1 && grepl(".grp$", rid) )
                      rid <- parse.grp(rid)
                  if ( length(cid) == 1 && grepl(".grp$", cid) )
                      cid <- parse.grp(cid)
                  # get the row and column ids
                  all_rid <- read.gctx.ids(src, dimension="row")
                  all_cid <- read.gctx.ids(src, dimension="col")
                  # if rid or cid specified, read only those rows/columns
                  if (!is.null(rid)) {
                    ridx <- match(rid, all_rid)
                  } else {
                    ridx <- seq_along(all_rid)
                  }
                  if (!is.null(cid)) {
                    cidx <- match(cid, all_cid)
                  } else {
                    cidx <- seq_along(all_cid)
                  }
                  # read the data matrix
                  .Object@mat <- h5read(src, name="0/DATA/0/matrix", index=list(ridx, cidx))
                  # set the row and column ids
                  .Object@rid <- all_rid[ridx]
                  .Object@cid <- all_cid[cidx]
                  colnames(.Object@mat) <- all_cid[cidx]
                  rownames(.Object@mat) <- all_rid[ridx]
                  # get the meta data
                  .Object@rdesc <- read.gctx.meta(src, dimension="row", ids=rid)
                  .Object@cdesc <- read.gctx.meta(src, dimension="col", ids=cid)
                  # close any open handles and return the object
                  H5close()
                  return(.Object)
              }
          }
)


# function to parse a GCT(X)
# just instantiates a new GCT object
parse.gctx <- function(fname, rid = NULL, cid = NULL) {
    ds <- new("GCT", src = fname, rid = rid, cid = cid)
    return(ds)
}


### method to extract rows and cols from an existing GCT object ###
gct.extract <- function(ds, rid = NULL, cid = NULL) {
    if (! is.null(rid)) {
        # these will be ordered as rid (that's how the intersect works in R)
        rid.ds <- intersect(rid, ds@rid)
        if (! identical(rid.ds, rid)) {
            missings = setdiff(rid, rid.ds)
            warning("The following rid's were not found: ", 
                    paste(missings, collapse = ", "))
        }
        ds@mat = ds@mat[rid.ds,,drop = FALSE]
        ds@rid = rid.ds
        ds@rdesc = ds@rdesc[rid.ds,]
    }
    # extract the columns
    if (! is.null(cid)) {
        # will be ordered as cid
        cid.ds = intersect(cid, ds@cid)
        if (! identical(cid.ds, cid)) {
            missings = setdiff(cid, cid.ds)
            warning("The following cid's were not found: ",
                    paste(missings, collapse = ", "))
        }
        ds@mat = ds@mat[,cid.ds,drop = FALSE]
        ds@cid = cid.ds
        ds@cdesc = ds@cdesc[cid.ds,]
    }
    return(ds)
}


append.dim <- function(ofile, mat, extension="gct") {
  nc <- ncol(mat)
  nr <- nrow(mat)
  outFile <- basename(ofile)
  filename <- strsplit(outFile,'.',fixed=T)[[1]][1]
  ofile <- path.join(dirname(ofile),
                    sprintf('%s_n%dx%d.%s',filename,
                            nc, nr, extension))
  return(ofile)
}


# write a gct file to disk
write.gct <- function(ofile, ds, precision=4, appenddim=T, ver=3) {
  # gct must contain the following fields
  #          mat: Numeric data matrix [RxC]
  #          rid: Cell array of row ids
  #          rhd: Cell array of row annotation fieldnames
  #          rdesc: Cell array of row annotations
  #          cid: Cell array of column ids
  #          chd: Cell array of column annotation fieldnames
  #          cdesc: Cell array of column annotations
  #          version: GCT version string
  #          src: Source filename
  
  
  # append the dimensions of the data set, if desired
  if (appenddim) ofile <- append.dim(ofile, ds@mat, extension="gct")
  
  precision = floor(precision)
  cat(sprintf('Saving file to %s\n',ofile))
  cat(sprintf('Dimensions of matrix: [%dx%d]\n',nr,nc))
  cat(sprintf('Setting precision to %d\n',precision))
  
  # open file      
  if (ver==3) {
    nrdesc = dim(ds@rdesc)[2]
    ncdesc = dim(ds@cdesc)[2]
    colkeys = colnames(ds@cdesc)
    # append header
    cat(sprintf('#1.%d\n%d\t%d\t%d\t%d', ver, nr, nc, nrdesc, ncdesc),
        file=ofile,sep='\n')      
    # line 3: sample row desc keys and sample names
    cat(paste(c('id',colnames(ds@rdesc),ds@cid),collapse='\t'),
        file=ofile,sep='\n',append=T)
    # line 4 + ncdesc: sample desc
    filler = 'na'
    for (ii in 1:ncdesc) {
      if (is.numeric(ds@cdesc[,ii])) {
        cat(paste(c(colkeys[ii],rep(filler,nrdesc),
                    round(ds@cdesc[,ii],precision)),
                  collapse='\t'),
            file=ofile,sep='\n',append=T)  
      } else {
        cat(paste(c(colkeys[ii],rep(filler,nrdesc),
                    ds@cdesc[,ii]),
                  collapse='\t'),
            file=ofile,sep='\n',append=T)
      }
    }
  } else {
    # append header
    cat(sprintf('#1.%d\n%d\t%d\t%d\t%d', ver, nr, nc),
        file=ofile,sep='\n')      
    # line 3: sample row desc keys and sample names
    cat(paste(c('id','Description',ds@cid),collapse='\t'),
        file=ofile,sep='\n',append=T)
  }
  
  for (ii in 1:nr) {    
    # print rows
    cat(paste(c(ds@rid[ii],
                ds@rdesc[ii,],
                round(ds@mat[ii,],precision)),collapse='\t'),
        sep='\n',file=ofile,append=T)
  }
  cat(sprintf('Saved.\n'))  
}


# write a GCTX object
write.gctx <- function(ds, ofile, appenddim=T) {
  if (appenddim) ofile <- append.dim(ofile, ds@mat, extension="gctx")
  # check if the file already exists
  if (file.exists(ofile)) {
    message(paste(ofile, "exists, removing"))
    file.remove(ofile)
  }
  # start the file object
  h5createFile(ofile)
  # create all the necessary groups
  h5createGroup(ofile, "0")
  h5createGroup(ofile, "0/DATA")
  h5createGroup(ofile, "0/DATA/0")
  h5createGroup(ofile, "0/META")
  h5createGroup(ofile, "0/META/COL")
  h5createGroup(ofile, "0/META/ROW")
  
  # write data to each group
  h5write(ds@mat, ofile, "0/DATA/0/matrix")
  h5write(ds@rid, ofile, "0/META/ROW/id")
  h5write(ds@cid, ofile, "0/META/COL/id")
  write.meta(ofile, ds@cdesc, dimension="column")
  write.meta(ofile, ds@rdesc, dimension="row")
  h5write("GCTX1.0", ofile, "version")

  # close any open handles
  H5close()
}


# helper function to write a data.frame of meta data to gctx object
# makes an HDF5 entry for each column of the data.frame
write.meta <- function(ofile, df, dimension="row") {
  path <- if ((dimension=="row")) "0/META/ROW/" else "0/META/COL/"
  # loop through all columns
  fields <- names(df)
  if (length(fields) > 0) {
    for (i in 1:length(fields)) {
      field <- fields[i]
      v <- df[, i]
      # convert factors to character
      if(class(v) == "factor" || class(v) == "AsIs") {
        v <- as.character(v)
      }
      h5write(v, ofile, paste(path, field, sep=""))
    }
  }
}

###########################################
### functions for other CMap file types ###
###########################################

### function to read a .grp file and return a vector ###
parse.grp <- function(fname) {
    grp <- scan(fname, what = "", quote = NULL, quiet = TRUE)
    return(grp)
}


### function to read a .gmx file and return a list ###
parse.gmx <- function(fname) {
    tmp <- read.table(fname, sep = "\t", 
                     header = TRUE, stringsAsFactors = FALSE)
    # preallocate a list for the gmx
    L <- list()
    # loop over the first row of the .gmx
    for ( n in names(tmp) ) {
        # get all the values; remove empties at the end
        values <- tmp[[n]][-1]
        remove.idx <- values == ""
        values <- values[!remove.idx]
        # put in a list
        L[[n]] <- list(head = n,
                      desc = tmp[[n]][1], 
                      len = length(values), 
                      entry = values)
    }
    return(L)
}


### function to read a .gmt file and return a a list ###
parse.gmt <- function(fname) {
    gmt.lines <- scan(fname, what = "", sep = "\n",
                     quote = NULL, quiet = TRUE)
    tmp <- lapply(gmt.lines, function(x) unlist(strsplit(x, "\t")))
    mk.gmt.entry <- function(x) {
        L <- list()
        L[["head"]] <- x[1]
        L[["desc"]] <- x[2]
        l.entry <- x[-c(1:2)]
        idx <- l.entry != ""
        L[["entry"]] <- l.entry[idx]
        L[["len"]] <- length(L[["entry"]])
        return(L)
    }
    L <- lapply(tmp, function(x) mk.gmt.entry(x))
    names(L) <- unlist(lapply(L, function(x) x$head))
    return(L)
}


### function to write tab-delimited text files at a fixed numerical precision ###
mktbl <- function(tbl, ofile, precision = 4, col.names = TRUE, row.names = TRUE) {
    # format numeric floats to two decimal points; leave all others as is
    for (col in names(tbl))
        if (class(tbl[[col]]) == "numeric")
            if (! (all(round(tbl[[col]]) == tbl[[col]]))) {
                # the format string; tells us the precision
                fmt <- paste("%0.", as.character(precision), "f", sep = "")
                tbl[[col]] <- sprintf(fmt, tbl[[col]])
            }
    write.table(tbl, file = ofile, sep = "\t", quote = FALSE, 
                col.names = col.names, row.names = row.names)
}

########################################
### Other Misc. utility functions ######
########################################

### function to join a bunch of arguments into a file path (same as matlab's fullfile) ###
# function to join a bunch of arguments into a file path (same as matlab's fullfile)
path.join <- function(...) {
    args = c(...)
    args = sub("/$", "", args)
    path = paste(args, collapse = "/")
    return(path)
}
