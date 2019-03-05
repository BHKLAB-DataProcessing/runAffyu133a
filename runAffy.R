library(affy)
library(downloader)
library(PharmacoGxPrivate)
options(stringsAsFactors=FALSE)

load("/pfs/gdscCellInfo/celline.gdsc.RData")
celline <- celline.gdsc

file.paths <- file.path("/pfs/BrainArray/",c(list.files(pattern="hgu133ahsensg*", path="/pfs/BrainArray/"),
                list.files(pattern="pd.hgu133a.hs.ensg*", path="/pfs/BrainArray/")))

celfn <- list.celfiles("/pfs/gdscU133a", full.names=TRUE)
celfns <- list.celfiles("/pfs/gdscU133a", full.names=FALSE)
## experiments' names
names(celfn) <- names(celfns) <- gsub(".CEL.gz", "", celfns)
## chip type and date
chipt <- sapply(celfn, celfileChip)
chipd <- t(sapply(celfn, celfileDateHour))

message("Read sample information")
sampleinfo <- read.csv(file.path("/pfs/gdscU133a", "gdsc_ge_sampleinfo.txt"), sep="\t")
sampleinfo[sampleinfo == "" | sampleinfo == " "] <- NA
## curate cell line names
sampleinfo[sampleinfo[ , "Source.Name"] == "MZ2-MEL.", "Source.Name"] <- "MZ2-MEL"
sampleinfo[sampleinfo[ , "Source.Name"] == "KMS-12-BM", "Source.Name"] <- "KMS12-BM"
iix <- which(!duplicated(sampleinfo[ , "Source.Name"]) & !is.element(sampleinfo[ , "Source.Name"], celline[ , "Sample.name"]))
if(length(iix) > 0) {
  ## enrich the list of cell lines
  tt <- matrix(NA, nrow=length(iix), ncol=ncol(celline), dimnames=list(sampleinfo[iix, "Source.Name"], colnames(celline)))
  tt[ , "Sample.name"] <- sampleinfo[iix, "Source.Name"]
  celline <- rbind(celline, tt)
}
fn <- gsub(patter="[.]CEL", replacement="", x=sampleinfo[ , "Array.Data.File"])
if(any(!is.element(fn[!is.na(fn)], names(celfns)))) { stop("some CEL files are missing for the GDSC project") }
rownames(sampleinfo) <- fn
sampleinfo <- sampleinfo[names(celfn), , drop=FALSE]
sampleinfo <- data.frame("samplename"=names(celfns), "filename"=celfns, "chiptype"=chipt, "hybridization.date"=chipd[ , "day"], "hybridization.hour"=chipd[ , "hour"], "file.day"=celfile.timestamp[ , "file.day"], "file.hour"=celfile.timestamp[ , "file.hour"], "batchid"=NA, "cellid"=sampleinfo[ , "Source.Name"], sampleinfo)
  



## phenodata
load("/pfs/gdscU133a/celfile_timestamp.RData")

rownames(celfile.timestamp) <- basename(rownames(celfile.timestamp))

celfn <- list.files(pattern="*.cel.gz", path="/pfs/gdscU133a/", full.names=TRUE)

cgp.u133a <- just.rma(filenames=celfn, cdfname="hgu133ahsensgcdf")
save(cgp.u133a, compress=TRUE, file="GDSC_U219_ENSG_RAW.RData")
print(head(rownames(pData(cgp.u133a))))
pData(cgp.u133a) <- data.frame(pData(cgp.u133a), sampleinfo[match(colnames(exprs(cgp.u133a)), sampleinfo[ , "Array.Data.File"]), , drop=FALSE], celfile.timestamp[rownames(pData(cgp.u133a)), , drop=FALSE])
colnames(exprs(cgp.u133a)) <- rownames(pData(cgp.u133a)) <- colnames(exprs(cgp.u133a))
fData(cgp.u133a) <- data.frame("PROBE"=rownames(exprs(cgp.u133a)), "GENEID"=sapply(strsplit(rownames(exprs(cgp.u133a)), "_"), function (x) { return (x[[1]]) }), "BEST"=TRUE)
rownames(fData(cgp.u133a)) <- rownames(exprs(cgp.u133a))
cgp.u133a.ensg <- cgp.u133a
save(cgp.u133a.ensg, compress=TRUE, file="/pfs/out/GDSC_U133a_ENSG.RData")
