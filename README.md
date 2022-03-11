# SNPDensity[Total.pdf](https://github.com/PBGLMichaelHall/SNPDensity/files/8109884/Total.pdf)
# VCF File Header UpStream Filtering 
![Rice_UpstreamFiltering](https://user-images.githubusercontent.com/93121277/157855864-b4fa542c-3aae-4885-a098-e40038e0d592.png)
# VCF File Header with Samples Listed
![RiceHeader](https://user-images.githubusercontent.com/93121277/157856014-65cc5b6d-661b-4499-b82e-af6b3190f3f6.png)
# Separate Samples with Control Sample 15
![Rice1_15](https://user-images.githubusercontent.com/93121277/157856783-73a6f621-7fa9-4972-999b-66edeb903bfd.png)
# Run the RScript
```{r cars}
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

```
