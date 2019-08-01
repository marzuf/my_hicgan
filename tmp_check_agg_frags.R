require(read.table)

SSHFS <- FALSE

setDir <- ifelse(SSHFS, "/media/electron", "")

binSize <- 10*10^3

cat("... LOAD CMAP DATA \n")
cmap_file <- file.path(setDir, "/mnt/pd2/shared/EZH2_Interactome/KARPAS/DMSO/Separated_filtered/Cmap_KARPAS_DMSO_chr1_library1_Rep1AB_mapq1_10kb.txt")
cmap_dt <- read.table(cmap_file, header=FALSE, stringsAsFactors = FALSE)
colnames(cmap_dt) <- c("binA", "binB", "count", "dist")
cat("... TEST IF NUMERIC \n")
stopifnot(is.numeric(cmap_dt$binA))
stopifnot(is.numeric(cmap_dt$binB))
stopifnot(is.numeric(cmap_dt$count))
stopifnot(is.numeric(cmap_dt$dist))
head(cmap_dt)

cat("... LOAD FRAG DATA \n")
frag_file <- file.path(setDir, "/mnt/pd2/shared/EZH2_Interactome/KARPAS/DMSO/Separated_filtered/KARPAS_DMSO_chr1_library1_Rep1AB.rmdup.threshold0_mapq1.txt")
frag_DT <- read.table(frag_file, header=FALSE, stringsAsFactors = FALSE)
frag_DT <- frag_DT[,1:6]
colnames(frag_DT) <- c("strandA", "chrA", "posA", "strandB", "chrB", "posB")
cat("... TEST CHROMO EQUALITY \n")
stopifnot(frag_DT$chrA == frag_DT$chrB)
frag_DT$strandA <- NULL
frag_DT$strandB <- NULL
cat("... TEST IF NUMERIC \n")
stopifnot(is.numeric(frag_DT$posA))
stopifnot(is.numeric(frag_DT$posB))
cat("... COMPUTE BIN \n")
frag_DT$binA <- frag_DT$posA%/%binSize
frag_DT$binB <- frag_DT$posB%/%binSize

head(frag_DT)

frag_DT$posA <- frag_DT$posB <- NULL
cat("... AGGREGATE FRAGMENT \n")
agg_frag_DT <- aggregate(. ~ binA + binB, data = frag_DT, FUN=length)
cat("... TEST CHROMO EQUALITY \n")
stopifnot(agg_frag_DT$chrA == agg_frag_DT$chrB)
agg_frag_DT$chrB <- NULL
colnames(agg_frag_DT)[colnames(agg_frag_DT) == "chrA"] <- "count"

head(agg_frag_DT)

cat(paste0("nrow(cmap_dt)\t=\t", nrow(cmap_dt) , "\n"))
cat(paste0("nrow(agg_frag_DT)\t=\t", nrow(agg_frag_DT) , "\n"))

cat("\n... retrieve some from Cmap\n")
cmap_dt[cmap_dt$binA %in% c(70:72) & cmap_dt$binB %in% c(70:72),]

cat("\n... retrieve some from agg_frag \n")
agg_frag_DT[agg_frag_DT$binA %in% c(70:72) & agg_frag_DT$binB %in% c(70:72),]

cmap_dt <- cmap_dt[order(cmap_dt$binA, cmap_dt$binB),]

agg_frag_DT <- agg_frag_DT[order(agg_frag_DT$binA, agg_frag_DT$binB),]

cat("\n... HEAD Cmap\n")
head(cmap_dt,10)
cat("\n... HEAD agg_frag\n")
head(agg_frag_DT,10)

rownames(cmap_dt) <- NULL
rownames(agg_frag_DT) <- NULL
all.equal(cmap_dt[,1:3], agg_frag_DT[,1:3])