# Inspect the Data Structure Tree
![tree9](https://user-images.githubusercontent.com/93121277/158329889-1a9870e3-e656-4de3-94d7-4c3176136a43.png)
# SNPDensity[Total.pdf](https://github.com/PBGLMichaelHall/SNPDensity/files/8109884/Total.pdf)
# VCF File Header UpStream Filtering 
![Rice_UpstreamFiltering](https://user-images.githubusercontent.com/93121277/157855864-b4fa542c-3aae-4885-a098-e40038e0d592.png)
# VCF File Header with Samples Listed
![RiceHeader](https://user-images.githubusercontent.com/93121277/157856014-65cc5b6d-661b-4499-b82e-af6b3190f3f6.png)
# Separate Samples with Control Sample 15
![Rice1_15](https://user-images.githubusercontent.com/93121277/157856783-73a6f621-7fa9-4972-999b-66edeb903bfd.png)


# Decompress the Separated VCF Files

```r


library(R.utils)
message("SETTING WORKING DIRECTORY")
setwd("/home/michael/Desktop/SNPVCFPRACTICE")


message("CALLING THE DECOMPRESS FUNCTION IN R")
Decompress <- function(file){
vcf <- gunzip(file,remove=FALSE)
return(vcf)
}



Decompress(snakemake@input[[1]])
Decompress(snakemake@input[[2]])
Decompress(snakemake@input[[3]])
Decompress(snakemake@input[[4]])
Decompress(snakemake@input[[5]])
Decompress(snakemake@input[[6]])
Decompress(snakemake@input[[7]])
Decompress(snakemake@input[[8]])
Decompress(snakemake@input[[9]])
Decompress(snakemake@input[[10]])
Decompress(snakemake@input[[11]])
Decompress(snakemake@input[[12]])
Decompress(snakemake@input[[13]])


message("PREPARING TO REMOVE GUNZIPPED FILES FROM DIRECTORY")
#Decompress the VCF files
setwd("/home/michael/Desktop/SNPVCFPRACTICE/VCF/")
x <- 'rm *.gz'
system(x)

q()

```
# Run the RScript
```r
library(magick)
library(vcfR)
library(rMVP)
#Set the working directory where the Variant Calling File is located
message("Setting Working Directory")
setwd("/home/michael/Desktop/SNPVCFPRACTICE")
message("Calling MVP Data function on all decompressed variant calls")

#Call the MVP Data function
message("Calling the MVP for Loop function on all samples 1 - 13")

Loop <- function(file,sample){
sample<-sample
pathtosample <- file
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


Loop(snakemake@input[[1]],"S1")
Loop(snakemake@input[[2]],"S2")
Loop(snakemake@input[[3]],"S3")
Loop(snakemake@input[[4]],"S4")
Loop(snakemake@input[[5]],"S5")
Loop(snakemake@input[[6]],"S6")
Loop(snakemake@input[[7]],"S7")
Loop(snakemake@input[[8]],"S8")
Loop(snakemake@input[[9]],"S9")
Loop(snakemake@input[[10]],"S10")
Loop(snakemake@input[[11]],"S11")
Loop(snakemake@input[[12]],"S12")
Loop(snakemake@input[[13]],"S13")




message("CONVERTING JPG IMAGES TO PNG IMAGES")
#Convert .jpg image to png image

for (i in 1:13){
snp<-image_read(paste0("S",i,".SNP_Density.jpg"))
image <- image_convert(image = snp, "png")
image_write(image,paste0("S",i,".SNP_Density.png",format = "png"))
}

message("MAKING FINAL PDF REPORT")

y <- 'convert S1.SNP_Density.png S2.SNP_Density.png  S3.SNP_Density.png
 S4.SNP_Density.png  S5.SNP_Density.png  S6.SNP_Density.png  S7.SNP_Density.png
 S8.SNP_Density.png  S9.SNP_Density.png  S10.SNP_Density.png
 S11.SNP_Density.png  S12.SNP_Density.png  S13.SNP_Density.png Total.pdf'

system(y)

z <- 'rm *.map *.ind *.bin *.desc'
system(z)



```

# The Snakemake file

![snakemakefile2](https://user-images.githubusercontent.com/93121277/158329199-8de5b090-77ab-4bb8-ab2a-6cac8be26441.png)



# In a SNAKEMAKE Workflow pipeline


![JobStatus](https://user-images.githubusercontent.com/93121277/158328868-f1b2cc6d-3c10-4f08-a571-a3be494d0f6b.png)


# The Directed Acrylic Graph

![DAG](https://user-images.githubusercontent.com/93121277/158328677-4f7e59be-fb09-42b3-ab7e-29736a17772a.png)




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

