library(getopt)
##--------------------------------------------------------------
spec = matrix(c(
  'input','i', 1, "character","input file",
  'prefix','p', 2, "character","output prefix (default \"out\")",
  'help','h', 0, "logical","print help"
), byrow=TRUE, ncol=5)
opt = getopt(spec) 

if (!is.null(opt$help) || is.null(opt$input)) {
  cat(paste(getopt(spec, usage = T), "\n"))
  q()
}

input_file=opt$input
output_prefix="out"
if (!is.null(opt$prefix)) {
  output_prefix=opt$prefix
}

cat("----------------NDID is running----------------------------\n")
##----------------------------------------------------

library(fANCOVA)
library(qvalue)

##--------------------------------------------------------

yab <- read.delim(input_file) 

yab[yab$selfAvg_1 == 0, "selfAvg_1"] <- 0.5
yab[yab$selfAvg_2 == 0, "selfAvg_2"] <- 0.5

yab$D = 0.5*log2(yab$selfAvg_2) + 0.5*log2(yab$selfAvg_1)
outlier = boxplot(yab$D, plot=FALSE)$out

if(length(outlier) != 0) {
  yab = yab[-which(yab$D %in% outlier),]
} else {
  yab = yab
}

c1 = sum(yab$ipet_1)
c2 = sum(yab$ipet_2)

if(c1 > c2){
  yab$ipet_2 = yab$ipet_2 *(c1/c2)
} else{
  yab$ipet_1 = yab$ipet_1 *(c2/c1)
}

d1 = sum(yab$selfAvg_1)
d2 = sum(yab$selfAvg_2)

if(d1 > d2){
  yab$selfAvg_2 = yab$selfAvg_2 *(d1/d2)
} else{
  yab$selfAvg_1 = yab$selfAvg_1 *(d2/d1)
}

y1 = yab$ipet_1
y2 = yab$ipet_2

M = log2(y2/y1)
D = 0.5*log2(yab$selfAvg_2) + 0.5*log2(yab$selfAvg_1)

png(file="Before_Normalization.png", width=500, height=389)
model1 = loess.as(D,M, plot = TRUE,pch=20)
abline(h=0, lwd=2, lty=1, col="red")
dev.off()

log_y1_norm = log2(y1) + (model1$fitted/2)
log_y2_norm = log2(y2) - (model1$fitted/2)

y1_norm = 2^log_y1_norm
y2_norm = 2^log_y2_norm
M_norm = log2(y2_norm/y1_norm) 

png(file="After_Normalization.png", width=500, height=389)
model2 = loess.as(D, M_norm, plot = TRUE,pch=20)
dev.off

zi = (M_norm - mean(M_norm))/sd(M_norm)

p.val = 2*pnorm(abs(zi), lower.tail = FALSE)
fdr = p.adjust(p.val, method = "BH", n = length(p.val))

intensity = c()

for (i in 1:nrow(yab)) {
  if (M_norm[i] > 0) {
    intensity[i] = "1"
  } else 
    intensity[i] = "0"
} 

output = data.frame(yab$chrom1, yab$start1,yab$end1, 
                   yab$chrom2, yab$start2, yab$end2,
                   y1_norm, y2_norm,p.val, fdr,intensity)

colnames(output) <- c("chrom1","start1", "end1", "chrom2",
                      "start2", "end2", "Normalized_ipet_1","Normalized_ipet_2",
                      "P-value","p.adjust", "intensity type")

write.table(output, paste(output_prefix,"significant_interaction.txt",sep="_"), sep="\t", row.names = F, quote=FALSE)
