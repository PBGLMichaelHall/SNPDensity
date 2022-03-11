# SNPDensity[Total.pdf](https://github.com/PBGLMichaelHall/SNPDensity/files/8109884/Total.pdf)
# VCF File Header UpStream Filtering 
![Rice_UpstreamFiltering](https://user-images.githubusercontent.com/93121277/157855864-b4fa542c-3aae-4885-a098-e40038e0d592.png)
# VCF File Header with Samples Listed
![RiceHeader](https://user-images.githubusercontent.com/93121277/157856014-65cc5b6d-661b-4499-b82e-af6b3190f3f6.png)
# Separate Samples with Control Sample 15
![Rice1_15](https://user-images.githubusercontent.com/93121277/157856783-73a6f621-7fa9-4972-999b-66edeb903bfd.png)
# Inspect the Data Structure Tree

![TREE](https://user-images.githubusercontent.com/93121277/157858372-a2b9b953-6913-4f2c-a2ae-9a9c6a6b7604.png)


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
# SNP Denisity Plots for Samples 1 - 14 against Control Sample 15
# Notice Chromosome 6 has magnitude of significance in TOTAL SNPs called
![SNP1](https://user-images.githubusercontent.com/93121277/157857578-1394e266-822e-4ae8-885c-79a80fdd5a50.png)
![SNP2](https://user-images.githubusercontent.com/93121277/157857581-2e4737b6-a004-404f-abec-a0c3beac40e0.png)
![SNP3](https://user-images.githubusercontent.com/93121277/157857586-b0ddf88f-dc03-4de4-a5b4-1e32db6585c9.png)
![SNP4](https://user-images.githubusercontent.com/93121277/157857589-d745d9d4-84bb-4add-a44c-e4e6da628365.png)
![SNP5](https://user-images.githubusercontent.com/93121277/157857592-ba020eeb-23e6-4a53-8be2-845efd9fc252.png)
![SNP6](https://user-images.githubusercontent.com/93121277/157857595-c0d9acc4-6a7e-4b05-aa2d-314e26ce1389.png)
![SNP7](https://user-images.githubusercontent.com/93121277/157857599-de6c6cf6-a664-46fa-9fd9-61a4b40abc8c.png)
![SNP8](https://user-images.githubusercontent.com/93121277/157857601-54026508-1cc6-4028-8467-e8647cd93b7b.png)

