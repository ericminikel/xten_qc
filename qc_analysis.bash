
# next time try bsub -R rusage[mem=16G]

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
              -comp $wgs_vcf \
              -eval $wes_vcf \
              -moltenize \
              -o wes.vs.wgs.gq30dp10.molt"
cat wes.vs.wgs.gq30dp10.molt | grep -A 1 ^#:GATKTable:GenotypeConcordance_EvalProportions > wes.vs.wgs.gq30dp10.molt.all.concordance.proportions
cat wes.vs.wgs.gq30dp10.molt | grep -A 400 ^#:GATKTable:GenotypeConcordance_EvalProportions | grep ^ALL >> wes.vs.wgs.gq30dp10.molt.all.concordance.proportions


# GenotypeConcordance separately for SNPs and INDELs
bsub -q week -P $RANDOM -J gtc_snps -M 16000000 \
            -o gtconc_snps.o \
            -e gtconc_snps.e \
"java -Xmx15g -jar $gatkjar \
              -R $b37ref \
              -T GenotypeConcordance \
              -gfe 'GQ<30' \
              -gfe 'DP<30' \
              -gfc 'GQ<30' \
              -gfc 'DP<30' \
              -comp wes.snps.vcf.gz \
              -eval wgs.snps.vcf.gz \
              -moltenize \
              -o wgs.vs.wes.gq30dp30.snps.molt"
bsub -q week -P $RANDOM -J gtc_indels -M 16000000 \
            -o gtconc_indels.o \
            -e gtconc_indels.e \
"java -Xmx15g -jar $gatkjar \
              -R $b37ref \
              -T GenotypeConcordance \
              -gfe 'GQ<30' \
              -gfe 'DP<30' \
              -gfc 'GQ<30' \
              -gfc 'DP<30' \
              -comp wes.indels.vcf.gz \
              -eval wgs.indels.vcf.gz \
              -moltenize \
              -o wgs.vs.wes.gq30dp30.indels.molt"

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

## DiagnoseTargets instead of DepthOfCoverage
bsub -q week -P $RANDOM -J doc -M 24000000 \
    -o docg.o \
    -e docg.e \
"java -Xmx23g -jar $gatkjar \
     -R $b37ref \
     -T DiagnoseTargets \
     -BQ 20 \
     -MQ 20 \
     -min 10 \
     -o wgs_diagnosetargets.vcf \
     -I $wgs_bams \
     -L $gencode_cds"
bsub -q week -P $RANDOM -J doc -M 24000000 \
    -o doce.o \
    -e doce.e \
"java -Xmx23g -jar $gatkjar \
     -R $b37ref \
     -T DiagnoseTargets \
     -BQ 20 \
     -MQ 20 \
     -min 10 \
     -o wes_diagnosetargets.vcf \
     -I $wes_bams \
     -L $gencode_cds"


mkdir jobtemp
mkdir wes_bybam
mkdir wgs_bybam
while read bam
do
    sname=`echo $bam | sed 's/.*\///'`
    bsub -q week -P $RANDOM -J doc_e -M 16000000 \
        -o jobtemp/job.wes_bybam.$sname.out \
        -e jobtemp/job.wes_bybam.$sname.err \
        "java -Xmx15g -jar $gatkjar \
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
    bsub -q week -P $RANDOM -J doc_g -M 16000000 \
        -o jobtemp/job.wgs_bybam.$sname.out \
        -e jobtemp/job.wgs_bybam.$sname.err \
        "java -Xmx15g -jar $gatkjar \
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

bsub -q hour -P $RANDOM -J tx_wgs -M 8000000 -o tx_wgs.o -e tx_wgs.e "export PYTHONPATH=$PYTHONPATH; transmission.py --pedfile sample.fam --vcfpath $wgs_vcf_transmission --biallelic_only > wgs_tx.txt"
bsub -q hour -P $RANDOM -J tx_wes -M 8000000 -o tx_wes.o -e tx_wes.e "export PYTHONPATH=$PYTHONPATH; transmission.py --pedfile sample.fam --vcfpath $wes_vcf_transmission --biallelic_only > wes_tx.txt"

