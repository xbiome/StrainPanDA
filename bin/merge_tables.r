#!/usr/bin/env Rscript

library(getopt)
library(data.table)

####### a template for merging tables (2 columns) in one folder ##########
#######       merge according to 'Merge_Key'       #######################

spec <- matrix(c(
    'help','h',0,'logical','Show this help information',
    'colPattern','c',1,'character','Regex to match part of file names and use them as the column names [default: .+]',
    'folder', 'f', 1, 'character', 'Folder containing the inputs [default: cwd]',
    'fPattern', 'p', 1, 'character', 'File pattern for the input [default: *]',
    'withHeader', 'H', 1, 'logical', 'The input files contain a header line [default: FALSE]',
    'outFile', 'o', 1, 'character', 'Output file (will write in csv format) [default: stdout]',
    'nThreads', 't', 1, 'integer', 'Number of threads used by fread [default: 4]'
    ), byrow=T, ncol=5)

opt = getopt(spec);

usage_message <- 'This script combine a set of files with to columns into one using the first column as the merge key.\n'
if ( !is.null(opt$help) ) {
    cat(usage_message)
    cat(getopt(spec, usage=TRUE));
    q(status=1);
}
#set some reasonable defaults for the options that are needed,
#but were not specified.
if ( is.null(opt$colPattern ) ) {opt$colPattern = '[^_]+'}
if ( is.null(opt$folder ) ) { opt$folder = '.' }
if ( is.null(opt$fPattern ) ) { opt$fPattern = '*' }
if ( is.null(opt$withHeader ) ) { opt$withHeader = FALSE }
if ( is.null(opt$nThreads ) ) { opt$nThreads = 4 }


Merge_Key <- 'Gene'# merge according to this column
Col_Pattern <- opt$colPattern   # name the columns according to the file name
Folder <- opt$folder       # folder containing tables
File_Pattern <- opt$fPattern # Pattern for matching file names

# input files list
filenames <- list.files(Folder, pattern=File_Pattern, full.names=TRUE)

aux <- function(x) {
    ## closure to use "Merge_Key"
    names <- c(Merge_Key, regmatches(basename(x), regexpr(Col_Pattern, basename(x), perl=T)))
    out <- tryCatch({fread(x, header=opt$withHeader, sep='\t', check.name=F,
	            col.names=names,
		    nThread=opt$nThreads)
             },
             error = function(x) NA
	)
    out
}

df_list <- lapply(filenames, aux)
df_list <- Filter(function(x) any(!is.na(x)), df_list) # reomve empty

# merge DFs
setDTthreads(opt$nThreads)
df_merged <- Reduce(function(x,y)merge(x,y,all=T,by=Merge_Key), df_list)

# NA removal
df_merged[is.na(df_merged)] <- 0

## # write to table
if ( !is.null(opt$outFile ) ){
    fwrite(df_merged, opt$outFile, quote=F, nThread=opt$nThread)
}else{
    print(df_merged)
}

