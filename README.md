# Inspect the Data Structure Tree
![Screenshot from 2022-03-21 16-08-51](https://user-images.githubusercontent.com/93121277/159290869-30802323-21b5-401a-94e7-6c4a5f1c8c1d.png)


# The final output of the workflow
[Total.pdf](https://github.com/PBGLMichaelHall/SNPDensity/files/8316366/Total.pdf)



# VCF File Header UpStream Filtering 
![Rice_UpstreamFiltering](https://user-images.githubusercontent.com/93121277/157855864-b4fa542c-3aae-4885-a098-e40038e0d592.png)
# VCF File Header with Samples Listed
![RiceHeader](https://user-images.githubusercontent.com/93121277/157856014-65cc5b6d-661b-4499-b82e-af6b3190f3f6.png)
# Separate Samples with Control Sample 15
![Rice1_15](https://user-images.githubusercontent.com/93121277/157856783-73a6f621-7fa9-4972-999b-66edeb903bfd.png)

# RUN the R Script from Snakemake workflow
```r
library(magick)
library(vcfR)
library(rMVP)
#Set the working directory where the Variant Calling File is located
message("Setting Working Directory")
setwd("/home/michael/Desktop/SNPVCFPRACTICE")
message("Calling MVP Data function on all decompressed variant calls")

#Call the MVP Data function
message("Calling the MVP function")


sample<-"S1"
pathtosample <- "VCF/freebayes~bwa~IRGSP-1.0~S1~HOM-VAR.vcf"
out<- paste0("mvp.",sample,".vcf")
memo<-paste0(sample)
dffile<-paste0("mvp.",sample,".vcf.geno.map")

message("Making MVP data S1")
MVP.Data(fileVCF=pathtosample,
           #filePhe="Phenotype.txt",
           fileKin=FALSE,
           filePC=FALSE,
           out=out
)
message("Reading MVP Data S1")
df <- read.table(file = dffile, header=TRUE)
message("Making SNP Density Plots")
MVP.Report.Density(df[,c(1:3)], bin.size = 10000000, col = c("blue", "yellow", "red"), memo = memo, file.type = "jpg", dpi=300)


sample<-"S2"
pathtosample <- "VCF/freebayes~bwa~IRGSP-1.0~S2~HOM-VAR.vcf"
out<- paste0("mvp.",sample,".vcf")
memo<-paste0(sample)
dffile<-paste0("mvp.",sample,".vcf.geno.map")

message("Making MVP data S2")
MVP.Data(fileVCF=pathtosample,
           #filePhe="Phenotype.txt",
           fileKin=FALSE,
           filePC=FALSE,
           out=out
)
message("Reading MVP Data S2")
df <- read.table(file = dffile, header=TRUE)
message("Making SNP Density Plots")
MVP.Report.Density(df[,c(1:3)], bin.size = 10000000, col = c("blue", "yellow", "red"), memo = memo, file.type = "jpg", dpi=300)





sample<-"S3"
pathtosample <- "VCF/freebayes~bwa~IRGSP-1.0~S3~HOM-VAR.vcf"
out<- paste0("mvp.",sample,".vcf")
memo<-paste0(sample)
dffile<-paste0("mvp.",sample,".vcf.geno.map")

message("Making MVP data 3")
MVP.Data(fileVCF=pathtosample,
           #filePhe="Phenotype.txt",
           fileKin=FALSE,
           filePC=FALSE,
           out=out
)
message("Reading MVP Data 3")
df <- read.table(file = dffile, header=TRUE)
message("Making SNP Density Plots")
MVP.Report.Density(df[,c(1:3)], bin.size = 10000000, col = c("blue", "yellow", "red"), memo = memo, file.type = "jpg", dpi=300)


sample<-"S4"
pathtosample <- "VCF/freebayes~bwa~IRGSP-1.0~S4~HOM-VAR.vcf"
out<- paste0("mvp.",sample,".vcf")
memo<-paste0(sample)
dffile<-paste0("mvp.",sample,".vcf.geno.map")

message("Making MVP data 4")
MVP.Data(fileVCF=pathtosample,
           #filePhe="Phenotype.txt",
           fileKin=FALSE,
           filePC=FALSE,
           out=out
)
message("Reading MVP Data 4")
df <- read.table(file = dffile, header=TRUE)
message("Making SNP Density Plots")
MVP.Report.Density(df[,c(1:3)], bin.size = 10000000, col = c("blue", "yellow", "red"), memo = memo, file.type = "jpg", dpi=300)





sample<-"S5"
pathtosample <- "VCF/freebayes~bwa~IRGSP-1.0~S5~HOM-VAR.vcf"
out<- paste0("mvp.",sample,".vcf")
memo<-paste0(sample)
dffile<-paste0("mvp.",sample,".vcf.geno.map")

message("Making MVP data 5")
MVP.Data(fileVCF=pathtosample,
           #filePhe="Phenotype.txt",
           fileKin=FALSE,
           filePC=FALSE,
           out=out
)
message("Reading MVP Data 5")
df <- read.table(file = dffile, header=TRUE)
message("Making SNP Density Plots")
MVP.Report.Density(df[,c(1:3)], bin.size = 10000000, col = c("blue", "yellow", "red"), memo = memo, file.type = "jpg", dpi=300)


sample<-"S6"
pathtosample <- "VCF/freebayes~bwa~IRGSP-1.0~S6~HOM-VAR.vcf"
out<- paste0("mvp.",sample,".vcf")
memo<-paste0(sample)
dffile<-paste0("mvp.",sample,".vcf.geno.map")

message("Making MVP data 6")
MVP.Data(fileVCF=pathtosample,
           #filePhe="Phenotype.txt",
           fileKin=FALSE,
           filePC=FALSE,
           out=out
)
message("Reading MVP Data 6")
df <- read.table(file = dffile, header=TRUE)
message("Making SNP Density Plots")
MVP.Report.Density(df[,c(1:3)], bin.size = 10000000, col = c("blue", "yellow", "red"), memo = memo, file.type = "jpg", dpi=300)





sample<-"S7"
pathtosample <- "VCF/freebayes~bwa~IRGSP-1.0~S7~HOM-VAR.vcf"
out<- paste0("mvp.",sample,".vcf")
memo<-paste0(sample)
dffile<-paste0("mvp.",sample,".vcf.geno.map")

message("Making MVP data 7")
MVP.Data(fileVCF=pathtosample,
           #filePhe="Phenotype.txt",
           fileKin=FALSE,
           filePC=FALSE,
           out=out
)
message("Reading MVP Data 7")
df <- read.table(file = dffile, header=TRUE)
message("Making SNP Density Plots")
MVP.Report.Density(df[,c(1:3)], bin.size = 10000000, col = c("blue", "yellow", "red"), memo = memo, file.type = "jpg", dpi=300)


sample<-"S8"
pathtosample <- "VCF/freebayes~bwa~IRGSP-1.0~S8~HOM-VAR.vcf"
out<- paste0("mvp.",sample,".vcf")
memo<-paste0(sample)
dffile<-paste0("mvp.",sample,".vcf.geno.map")

message("Making MVP data 8")
MVP.Data(fileVCF=pathtosample,
           #filePhe="Phenotype.txt",
           fileKin=FALSE,
           filePC=FALSE,
           out=out
)
message("Reading MVP Data 8")
df <- read.table(file = dffile, header=TRUE)
message("Making SNP Density Plots")
MVP.Report.Density(df[,c(1:3)], bin.size = 10000000, col = c("blue", "yellow", "red"), memo = memo, file.type = "jpg", dpi=300)




sample<-"S9"
pathtosample <- "VCF/freebayes~bwa~IRGSP-1.0~S9~HOM-VAR.vcf"
out<- paste0("mvp.",sample,".vcf")
memo<-paste0(sample)
dffile<-paste0("mvp.",sample,".vcf.geno.map")

message("Making MVP data 9")
MVP.Data(fileVCF=pathtosample,
           #filePhe="Phenotype.txt",
           fileKin=FALSE,
           filePC=FALSE,
           out=out
)
message("Reading MVP Data 9")
df <- read.table(file = dffile, header=TRUE)
message("Making SNP Density Plots")
MVP.Report.Density(df[,c(1:3)], bin.size = 10000000, col = c("blue", "yellow", "red"), memo = memo, file.type = "jpg", dpi=300)





sample<-"S10"
pathtosample <- "VCF/freebayes~bwa~IRGSP-1.0~S10~HOM-VAR.vcf"
out<- paste0("mvp.",sample,".vcf")
memo<-paste0(sample)
dffile<-paste0("mvp.",sample,".vcf.geno.map")

message("Making MVP data 10")
MVP.Data(fileVCF=pathtosample,
           #filePhe="Phenotype.txt",
           fileKin=FALSE,
           filePC=FALSE,
           out=out
)
message("Reading MVP Data 10")
df <- read.table(file = dffile, header=TRUE)
message("Making SNP Density Plots")
MVP.Report.Density(df[,c(1:3)], bin.size = 10000000, col = c("blue", "yellow", "red"), memo = memo, file.type = "jpg", dpi=300)


sample<-"S11"
pathtosample <- "VCF/freebayes~bwa~IRGSP-1.0~S11~HOM-VAR.vcf"
out<- paste0("mvp.",sample,".vcf")
memo<-paste0(sample)
dffile<-paste0("mvp.",sample,".vcf.geno.map")

message("Making MVP data 11")
MVP.Data(fileVCF=pathtosample,
           #filePhe="Phenotype.txt",
           fileKin=FALSE,
           filePC=FALSE,
           out=out
)
message("Reading MVP Data 11")
df <- read.table(file = dffile, header=TRUE)
message("Making SNP Density Plots")
MVP.Report.Density(df[,c(1:3)], bin.size = 10000000, col = c("blue", "yellow", "red"), memo = memo, file.type = "jpg", dpi=300)





sample<-"S12"
pathtosample <- "VCF/freebayes~bwa~IRGSP-1.0~S12~HOM-VAR.vcf"
out<- paste0("mvp.",sample,".vcf")
memo<-paste0(sample)
dffile<-paste0("mvp.",sample,".vcf.geno.map")

message("Making MVP data 12")
MVP.Data(fileVCF=pathtosample,
           #filePhe="Phenotype.txt",
           fileKin=FALSE,
           filePC=FALSE,
           out=out
)
message("Reading MVP Data 12")
df <- read.table(file = dffile, header=TRUE)
message("Making SNP Density Plots")
MVP.Report.Density(df[,c(1:3)], bin.size = 10000000, col = c("blue", "yellow", "red"), memo = memo, file.type = "jpg", dpi=300)


sample<-"S13"
pathtosample <- "VCF/freebayes~bwa~IRGSP-1.0~S13~HOM-VAR.vcf"
out<- paste0("mvp.",sample,".vcf")
memo<-paste0(sample)
dffile<-paste0("mvp.",sample,".vcf.geno.map")

message("Making MVP data 13")
MVP.Data(fileVCF=pathtosample,
           #filePhe="Phenotype.txt",
           fileKin=FALSE,
           filePC=FALSE,
           out=out
)
message("Reading MVP Data 13")
df <- read.table(file = dffile, header=TRUE)
message("Making SNP Density Plots")
MVP.Report.Density(df[,c(1:3)], bin.size = 10000000, col = c("blue", "yellow", "red"), memo = memo, file.type = "jpg", dpi=300)


sample<-"S14"
pathtosample <- "VCF/freebayes~bwa~IRGSP-1.0~S14~HOM-VAR.vcf"
out<- paste0("mvp.",sample,".vcf")
memo<-paste0(sample)
dffile<-paste0("mvp.",sample,".vcf.geno.map")

message("Making MVP data 14")
MVP.Data(fileVCF=pathtosample,
           #filePhe="Phenotype.txt",
           fileKin=FALSE,
           filePC=FALSE,
           out=out
)
message("Reading MVP Data 14")
df <- read.table(file = dffile, header=TRUE)
message("Making SNP Density Plots")
MVP.Report.Density(df[,c(1:3)], bin.size = 10000000, col = c("blue", "yellow", "red"), memo = memo, file.type = "jpg", dpi=300)

```

# The Snakemake file


![Screenshot from 2022-03-21 16-13-01](https://user-images.githubusercontent.com/93121277/159291703-47df8d41-9a6b-4f10-a7c9-d7a44a6e18e3.png)


# In a SNAKEMAKE Workflow pipeline


![Screenshot from 2022-03-21 16-13-58](https://user-images.githubusercontent.com/93121277/159291848-757ea783-422c-4079-9551-5d695075ddbd.png)



# The Directed Acrylic Graph
![Screenshot from 2022-03-21 16-14-53](https://user-images.githubusercontent.com/93121277/159292012-ceddf55d-2502-48e6-8f4a-341b97d226ee.png)






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

