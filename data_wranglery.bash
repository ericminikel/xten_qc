
. private-paths.bash

cat $fam_raw | grep -v $substring_to_grepv > sample.fam # remove some samples from the FAM for analysis
awk '{print $2}' sample.fam > sample.list # list of samples for the transmission analysis

## subset the desired WES data
# generate sample list. for the duplicated samples use sequencing run C1675 and not C1575
zcat $wes_raw_vcf | grep -m 1 ^#CHROM | tr '\t' '\n' | grep -w -f wgs_sample.list - | grep -v C1575 > wes_sample.list
# now subset to those individuals and Gencode CDS only
bsub -q week -P $RANDOM -o subs1.o -e subs1.e "java -Xmx2g -jar $gatkjar \
   -R $b37ref \
   -T SelectVariants \
   --variant $wes_raw_vcf \
   -o $wes_vcf \
   -L $gencode_cds \
   -env \
   -sf wes_sample.list"

# # to redo with -env if the above step was run without -env the first time
# cat wes_sample.list | sed 's/C1675:://g' > wes_sample.new.list
# bsub -q priority -P $RANDOM -o subsx.o -e subsx.e "java -Xmx2g -jar $gatkjar \
#    -R $b37ref \
#    -T SelectVariants \
#    --variant wes.backup.vcf.gz \
#    -o $wes_vcf \
#    -L $gencode_cds \
#    -env \
#    -sf wes_sample.new.list"

   

# now rename the C1675 samples to match the WGS VCF
zcat $wes_vcf | grep ^# | sed 's/C1675:://g' > new_wes_header.txt
tabix -r new_wes_header.txt $wes_vcf > reheadered_wes_vcf.gz
mv reheadered_wes_vcf.gz $wes_vcf
tabix $wes_vcf

# subset wgs vcf to gencode cds only
bsub -q week -P $RANDOM -o subs2.o -e subs2.e "java -Xmx2g -jar $gatkjar \
   -R $b37ref \
   -T SelectVariants \
   --variant $wgs_raw_vcf \
   -o $wgs_vcf \
   -L $gencode_cds \
   -env \
   -sf wgs_sample.list"

# split SNPs and INDELs in case we want to do separate GenotypeConcordance
bsub -q week -P $RANDOM -o subs3.o -e subs3.e "java -Xmx2g -jar $gatkjar \
   -R $b37ref \
   -T SelectVariants \
   --variant $wgs_vcf \
   -selectType INDEL \
   -o wgs.indels.vcf.gz"
bsub -q week -P $RANDOM -o subs4.o -e subs4.e "java -Xmx2g -jar $gatkjar \
   -R $b37ref \
   -T SelectVariants \
   --variant $wgs_vcf \
   -selectType SNP \
   -o wgs.snps.vcf.gz"
bsub -q week -P $RANDOM -o subs4.o -e subs4.e "java -Xmx2g -jar $gatkjar \
   -R $b37ref \
   -T SelectVariants \
   --variant $wes_vcf \
   -selectType INDEL \
   -o wes.indels.vcf.gz"
bsub -q week -P $RANDOM -o subs5.o -e subs5.e "java -Xmx2g -jar $gatkjar \
   -R $b37ref \
   -T SelectVariants \
   --variant $wes_vcf \
   -selectType SNP \
   -o wes.snps.vcf.gz"

# further subset the list of samples for the transmission analysis
bsub -q week -P $RANDOM -o subsA.o -e subsA.e "java -Xmx2g -jar $gatkjar \
   -R $b37ref \
   -T SelectVariants \
   --variant $wes_vcf \
   -o $wes_vcf_transmission \
   -sf sample.list"
bsub -q week -P $RANDOM -o subsB.o -e subsB.e "java -Xmx2g -jar $gatkjar \
   -R $b37ref \
   -T SelectVariants \
   --variant $wgs_vcf \
   -o $wgs_vcf_transmission \
   -sf sample.list"
