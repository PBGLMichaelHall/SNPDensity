

library(rMVP)
#Set the working directory where the Variant Calling File is located
setwd("/home/michael/Desktop/vcf18012021/")
getwd()


#Call the MVP Data function
main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  n <- args[1]
  for (i in 1:n) {
    dffile<-paste0("mvp.S",i,".vcf.geno.map")
    out<- paste0("mvp.S",i,".vcf")
MVP.Data(fileVCF=args[i],
           #filePhe="Phenotype.txt",
           fileKin=FALSE,
           filePC=FALSE,
           out=tempfile("outfile")
)
  
df <- read.table(file = dffile, header=TRUE)
MVP.Report.Density(df[,c(1:3)], bin.size = 10000000, col = c("blue", "yellow", "red"), memo = memo, file.type = "jpg", dpi=300)
  }
}

main()



