#Install rMVP package
install.packages("rMVP")

#Make the package available in R Library
library(rMVP)

#Set the working directory where the Variant Calling File is located
setwd("/home/michael/Desktop/vcf18012021/output/")
getwd()

#Call the MVP Data function

for (i in 1:14) {
sample<-paste0("S",i)
pathtosample<-paste0("../freebayes~bwa~IRGSP-1.0~",sample,"~HOM-VAR.vcf")
out<- paste0("mvp.",sample,".vcf")
memo<-paste0(sample)
dffile<-paste0("mvp.",sample,".vcf.geno.map")

MVP.Data(fileVCF=pathtosample,
         #filePhe="Phenotype.txt",
         fileKin=FALSE,
         filePC=FALSE,
         out=out
)

df <- read.table(file = dffile, header=TRUE)
MVP.Report.Density(df[,c(1:3)], bin.size = 10000000, col = c("blue", "yellow", "red"), memo = memo, file.type = "jpg", dpi=300)
}



