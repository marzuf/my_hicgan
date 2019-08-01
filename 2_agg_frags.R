# Rscript agg_frags.R INPUT_FRAGS/KARPAS_DMSO_chr3_downsample16_merged.txt 10000 INPUT_AGG
# 0       3       60003   16      3       60070   1       NB500883:199:HKV3WBGX3:1:21210:20071:14953      30      ACGT
# 0       3       60009   16      3       75279296        1       NB500883:199:HKV3WBGX3:3:12506:3372:10845       30      ACGT
# 16      3       60015   0       3       75279018        1       NB500883:199:HKV3WBGX3:2:23306:21253:16242      23      ACGT
# 16      3       60015   0       3       75279018        1       NB500883:199:HKV3WBGX3:2:23306:21253:16242      23      ACGT

script_name <- "2_agg_frags.R"

startTime <- Sys.time()

cat("> START ", script_name, "\n")

require(data.table)
SSHFS <- FALSE
setDir <- ifelse(SSHFS, "/media/electron", "")

#### HARD-CODED
fileSuffixToRemove <- ".txt"
outSuffix <- "_agg.txt"
########################################################################################################
my_args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(my_args) == 3)
cat("... LOAD FRAG DATA \n")
frag_file <- my_args[1]
stopifnot(file.exists(frag_file))
binSize <- as.numeric(my_args[2])
stopifnot(!is.na(binSize))
out_dir <- my_args[3]

dir.create(out_dir, recursive = TRUE)
out_file <- file.path(out_dir, paste0(gsub(fileSuffixToRemove, "", basename(frag_file)), outSuffix))
########################################################################################################
frag_DT <- read.table(frag_file, header=FALSE, stringsAsFactors = FALSE)
frag_DT <- frag_DT[,1:6]
colnames(frag_DT) <- c("strandA", "chrA", "posA", "strandB", "chrB", "posB")

cat("... TEST CHROMO EQUALITY \n")
stopifnot(frag_DT$chrA == frag_DT$chrB)
frag_DT$strandA <- NULL
frag_DT$strandB <- NULL

cat("... TEST ORDERING \n")
stopifnot(frag_DT$posA <= frag_DT$posB)

cat("... TEST IF NUMERIC \n")
stopifnot(is.numeric(frag_DT$posA))
stopifnot(is.numeric(frag_DT$posB))

cat("... COMPUTE BIN \n")
frag_DT$binA <- frag_DT$posA%/%binSize
frag_DT$binB <- frag_DT$posB%/%binSize


frag_DT$posA <- frag_DT$posB <- NULL

cat("... AGGREGATE FRAGMENT \n")
agg_frag_DT <- aggregate(. ~ binA + binB, data = frag_DT, FUN=length)

cat("... TEST AGAIN CHROMO EQUALITY \n")
stopifnot(agg_frag_DT$chrA == agg_frag_DT$chrB)
agg_frag_DT$chrB <- NULL
colnames(agg_frag_DT)[colnames(agg_frag_DT) == "chrA"] <- "count"

cat("... TEST AGAIN ORDERING \n")
stopifnot(agg_frag_DT$binA <= agg_frag_DT$binB)

agg_frag_DT <- agg_frag_DT[order(agg_frag_DT$binA, agg_frag_DT$binB),]

cat("... write to file")
write.table(agg_frag_DT, file = out_file, col.names=FALSE, row.names=FALSE, quote=FALSE, append=FALSE, sep="\t")
cat("... written: ", out_file, "\n")

########################################################################################################
########################################################################################################
########################################################################################################
cat("*** DONE: ", script_name, "\n")
cat(paste0(startTime, "\n", Sys.time()))
