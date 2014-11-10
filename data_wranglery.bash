
. private-paths.bash

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
   -sf wes_sample.list"

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
   -sf wgs_sample.list"