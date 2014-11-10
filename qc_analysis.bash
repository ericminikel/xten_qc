
# GenotypeConcordance
bsub -q week -P $RANDOM -J gtconc -M 16000000 \
            -o gtconc.o \
            -e gtconc.e \
"java -Xmx15g -jar $gatkjar \
              -R $b37ref \
              -T GenotypeConcordance \
              -gfe 'GQ<30' \
              -gfe 'DP<10' \
              -comp $wes_vcf \
              -eval $wgs_vcf \
              -moltenize \
              -o wgs.vs.wes.gq30dp10.molt"

## DepthOfCoverage
# WGS
bsub -q week -P $RANDOM -J doc -M 16000000 \
    -o docg.o \
    -e docg.e \
"java -Xmx15g -jar $gatkjar \
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
bsub -q week -P $RANDOM -J doc -M 16000000 \
    -o doce.o \
    -e doce.e \
"java -Xmx15g -jar $gatkjar \
     -R $b37ref \
     -T DepthOfCoverage \
     -o wes_doc_20_1 \
     -I $wes_bams \
     -L $gencode_cds \
     --omitDepthOutputAtEachBase \
     --minBaseQuality 20 \
     --minMappingQuality 20 \
     --countType COUNT_FRAGMENTS"
