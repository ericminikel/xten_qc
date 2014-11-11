
# GenotypeConcordance
bsub -q week -P $RANDOM -J gtconc -M 16000000 \
            -o gtconc.o \
            -e gtconc.e \
"java -Xmx15g -jar $gatkjar \
              -R $b37ref \
              -T GenotypeConcordance \
              -gfe 'GQ<30' \
              -gfe 'DP<10' \
              -gfc 'GQ<30' \
              -gfc 'DP<10' \
              -comp $wes_vcf \
              -eval $wgs_vcf \
              -moltenize \
              -o wgs.vs.wes.gq30dp10.molt"
# extract the one line with the overall summary of concordance for all samples
cat wgs.vs.wes.gq30dp10.molt | grep -A 1 ^#:GATKTable:GenotypeConcordance_EvalProportions > wgs.vs.wes.gq30dp10.molt.all.concordance.proportions
cat wgs.vs.wes.gq30dp10.molt | grep -A 400 ^#:GATKTable:GenotypeConcordance_EvalProportions | grep ^ALL >> wgs.vs.wes.gq30dp10.molt.all.concordance.proportions

## DepthOfCoverage
# WGS
bsub -q week -P $RANDOM -J doc -M 24000000 \
    -o docg.o \
    -e docg.e \
"java -Xmx23g -jar $gatkjar \
     -R $b37ref \
     -T DepthOfCoverage \
     -o wgs_doc_20_1 \
     -I $wgs_bams \
     -L $gencode_cds \
     --omitDepthOutputAtEachBase \
     --minBaseQuality 20 \
     --minMappingQuality 20 \
     --countType COUNT_FRAGMENTS"
# WES
bsub -q week -P $RANDOM -J doc -M 24000000 \
    -o doce.o \
    -e doce.e \
"java -Xmx23g -jar $gatkjar \
     -R $b37ref \
     -T DepthOfCoverage \
     -o wes_doc_20_1 \
     -I $wes_bams \
     -L $gencode_cds \
     --omitDepthOutputAtEachBase \
     --minBaseQuality 20 \
     --minMappingQuality 20 \
     --countType COUNT_FRAGMENTS"


mkdir jobtemp
mkdir wes_bybam
mkdir wgs_bybam
while read bam
do
    sname=`echo $bam | sed 's/.*\///'`
    bsub -q week -P $RANDOM -J doc_e -M 8000000 \
        -o jobtemp/job.wes_bybam.$sname.out \
        -e jobtemp/job.wes_bybam.$sname.err \
        "java -Xmx7g -jar $gatkjar \
             -R $b37ref \
             -T DepthOfCoverage \
             -o wes_bybam/cov_$sname \
             -I $wes_bams \
             -L $gencode_cds \
             --omitDepthOutputAtEachBase \
             --minBaseQuality 20 \
             --minMappingQuality 20 \
             --countType COUNT_FRAGMENTS"
done < $wes_bams
while read bam
do
    sname=`echo $bam | sed 's/.*\///'`
    bsub -q week -P $RANDOM -J doc_g -M 8000000 \
        -o jobtemp/job.wgs_bybam.$sname.out \
        -e jobtemp/job.wgs_bybam.$sname.err \
        "java -Xmx7g -jar $gatkjar \
             -R $b37ref \
             -T DepthOfCoverage \
             -o wgs_bybam/cov_$sname \
             -I $wgs_bams \
             -L $gencode_cds \
             --omitDepthOutputAtEachBase \
             --minBaseQuality 20 \
             --minMappingQuality 20 \
             --countType COUNT_FRAGMENTS"
done < $wgs_bams

# will need to redo these later once $wgs_vcf_transmission is finished
bsub -q week -P $RANDOM -J tx_wgs -M 8000000 -o tx_wgs.o -e tx_wgs.e "export PYTHONPATH=$PYTHONPATH; transmission.py --pedfile sample.fam --vcfpath $wgs_vcf --biallelic_only > wgs_tx.txt"
bsub -q week -P $RANDOM -J tx_wes -M 8000000 -o tx_wes.o -e tx_wes.e "export PYTHONPATH=$PYTHONPATH; transmission.py --pedfile sample.fam --vcfpath $wes_vcf --biallelic_only > wes_tx.txt"

