library(lassosum)
library(data.table)
library(methods)
library(magrittr)
library(parallel)
#threads
cl <- makeCluster(8)

sum.stat <- "/ssd/LachanceLab/GWASData/METAL/META_X_Auto_LDPred2_mod23.txt.gz"
bfile <- "/ssd/LachanceLab/UKBB/PLINK/UKB_BED_AFTER_QC_removed182IDs"
# Read in and process the covariates
covariate <-  fread("/ssd/LachanceLab/UKBB/PLINK/Chr22_Afr_PrCa_118/UKB_Afr_PrCa_Chr22.cov")
pcs <- fread("EUR.eigenvec") %>%
  setnames(., colnames(.), c("FID","IID", paste0("PC",1:10)))
# Need as.data.frame here as lassosum doesn't handle data.table 
# covariates very well
cov <- merge(covariate, pcs)

# We will need the EUR.hg19 file provided by lassosum 
# which are LD regions defined in Berisa and Pickrell (2015) for the European population and the hg19 genome.
ld.file <- "EUR.hg19"
# output prefix
prefix <- "EUR"
# Read in the target phenotype file
target.pheno <- fread("/ssd/LachanceLab/UKBB/PLINK/Chr22_Afr_PrCa_118/Afr_PrCa_118.pheno")[,c("FID", "IID", "PrCa")]
# Read in the summary statistics
ss <- fread(sum.stat)
# Remove P-value = 0, which causes problem in the transformation
ss <- ss[!P == 0]
# Transform the P-values into correlation
cor <- p2cor(p = ss$P,
             n = ss$N,
             sign = log(ss$OR)
)
fam <- fread(paste0(bfile, ".fam"))
fam[,ID:=do.call(paste, c(.SD, sep=":")),.SDcols=c(1:2)]


# Run the lassosum pipeline
# The cluster parameter is used for multi-threading
# You can ignore that if you do not wish to perform multi-threaded processing
out <- lassosum.pipeline(
  cor = cor,
  chr = ss$CHR,
  pos = ss$BP,
  A1 = ss$A1,
  A2 = ss$A2,
  ref.bfile = bfile,
  test.bfile = bfile,
  LDblocks = ld.file, 
  cluster=cl
)
# Store the R2 results
target.res <- validate(out, pheno = as.data.frame(target.pheno), covar=as.data.frame(cov))
# Get the maximum R2
r2 <- max(target.res$validation.table$value)^2