#A for loop for every sample in the VCF
for i in {1..14}
do
echo S$i;
#This command decompresses the vcf files
gunzip freebayes~bwa~IRGSP-1.0~S$i~HOM-VAR.vcf
#This command runs the R Script
Rscript Myscript.R
#If you do not have imagemagick download it
#sudo apt-get install imagemagick
#Convert the images into a pdf
convert output/S1.SNP_Density.jpg output/S2.SNP_Density.jpg output/S3.SNP_Density.jpg output/S4.SNP_Density.jpg output/S5.SNP_Density.jpg output/S6.SNP_Density.jpg output/S7.SNP_Density.jpg output/S8.SNP_Density.jpg output/S9.SNP_Density.jpg output/S10.SNP_Density.jpg output/S11.SNP_Density.jpg output/S12.SNP_Density.jpg output/S13.SNP_Density.jpg output/S14.SNP_Density.jpg / Total2.pdf
done



