

bsub -q bweek -P $RANDOM -J gtconc -M 8000000 \
            -o gtconc.o \
            -e gtconc.e \
"java -Xmx8g -jar $gatkjar \
              -R $b37ref \
              -T GenotypeConcordance \
              -gfe 'GQ<30' \
              -gfe 'DP<10' \
              -comp $wes_vcf \
              -eval $wgs_vcf \
              -moltenize \
              -o wgs.vs.wes.gq30dp10.molt"

