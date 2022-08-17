#!/bin/bash

#Activate env
source ~/fastq2gvcf/bin/activate

##Set variables
REF=~/sympn29qe/ra57gic/ref/Bos_taurus.ARS-UCD1.2.dna.toplevel.fa
m=$1
n=${m%%_1.fastq.gz}
id=$(gunzip -c $HOME/data/${n}_1.fastq.gz  | head -n 1 | cut -f 3-4 -d ":" | sed 's/:/_/')

##Mapping
bwa-mem2 mem ${REF} -t 16 -M -R "@RG\tID:$id\tLB:${n}\tPL:ILLUMINA\tPM:HISIEQ\tSM:${n}" $HOME/data/${n}_1.fastq.gz $HOME/data/${n}_2.fastq.gz |samtools view -@16 -Shb - > $SCRATCH/SACP/work/${n}_unsorted.bam

##SortSam
picard -Xmx28g SortSam -I $SCRATCH/SACP/work/${n}_unsorted.bam -O $SCRATCH/SACP/work/${n}_sorted.bam -SORT_ORDER coordinate -CREATE_INDEX true -TMP_DIR $SCRATCH/SACP/tmp

##MarkDuplicates
picard -Xmx28g MarkDuplicates -I $SCRATCH/SACP/work/${n}_unsorted.bam -O $SCRATCH/SACP/work/${n}_sorted.bam -METRICS_FILE ${n}_sorted.bam.dedup.metrics.txt -VALIDATION_STRINGENCY LENIENT -OPTICAL_DUPLICATE_PIXEL_DISTANCE 100 -CREATE_INDEX true -TMP_DIR $SCRATCH/SACP/tmp

##GATK BaseRecalibrator
gatk --java-options "-Xms28G -Xmx28G -XX:+UseParallelGC -XX:ParallelGCThreads=28" BaseRecalibrator -R ${REF} -I $SCRATCH/SACP/work/${n}.dedup.bam --known-sites ~/SACP/scripts/sympn29qe/ra57gic/ref/new_known.vcf.gz --bqsr-baq-gap-open-penalty 45 -O $SCRATCH/SACP/work/${n}.recal_data.table --tmp-dir $SCRATCH/SACP/tmp

##Apply BQSR
gatk --java-options "-Xmx28G -XX:+UseParallelGC -XX:ParallelGCThreads=28" ApplyBQSR -I $SCRATCH/SACP/work/${n}.dedup.bam -bqsr $SCRATCH/SACP/work/${n}.recal_data.table -O $SCRATCH/SACP/work/${n}.recal.bam --tmp-dir $SCRATCH/SACP/tmp

#Variant calling:Create GVCF file
gatk --java-options "-Xmx28G -XX:+UseParallelGC -XX:ParallelGCThreads=28" HaplotypeCaller -R ${REF} -I $SCRATCH/SACP/work/${n}.recal.bam -O $SCRATCH/SACP/gvcf/${n}.raw.g.vcf --native-pair-hmm-threads 28 -pairHMM AVX_LOGLESS_CACHING --emit-ref-confidence GVCF --dont-use-soft-clipped-bases true --tmp-dir $SCRATCH/SACP/tmp

#Generate list of intervals 
bedtools makewindows -g Bos_taurus.ARS-UCD1.2.dna.toplevel.fa.fai -w 5000000 > Bos_taurus.ARS-UCD1.2_intervals.bed 

picard -Xmx28g BedToIntervalList I=Bos_taurus.ARS-UCD1.2_intervals.bed O=Bos_taurus.ARS-UCD1.2.interval_list SD=Bos_taurus.ARS-UCD1.2.dna.toplevel.dict 

gatk SplitIntervals -R $REF -L Bos_taurus.ARS-UCD1.2.interval_list --scatter-count 2000 -O $interval_files_folder

#Genotype GVCF - GenomicsDB
gatk --java-options "-Xmx20g -XX:+UseParallelGC -XX:ParallelGCThreads=20 -DGATK_STACKTRACE_ON_USER_EXCEPTION=true" GenomicsDBImport --genomicsdb-workspace-path $SCRATCH/SACP/${n}_my_database --batch-size 50 --tmp-dir $SCRATCH/SACP/tmp -L sympn29qe/ra57gic/ref2/intervals/${n}-scattered.interval_list --sample-name-map sample.map --reader-threads 5

#Genotype GenomicsDB
gatk --java-options "-Xmx20G -XX:+UseParallelGC -XX:ParallelGCThreads=20" GenotypeGVCFs -R ${REF} -V gendb://$SCRATCH/SACP/GenomicsDB/${n}_my_database -O $SCRATCH/SACP/gvcf/${n}.vcf.gz --tmp-dir $SCRATCH/SACP/tmp

#Gather resulting vcfs 
picard -Xmx28g GatherVcfs I=*.vcf.gz O=combined.vcf.gz #Where each vcf file is listed with I e.g. 0000.g.vcf.gz

#sortvcf
picard -Xmx28g SortVcf I=combined.vcf.gz O=combined_sorted.vcf.gz CREATE_INDEX=true TMP_DIR=../tmp

#create an index file
tabix -p vcf combined.vcf.gz

##Extract SNPs and Indels 
gatk --java-options "-Xmx28G -XX:+UseParallelGC -XX:ParallelGCThreads=28" SelectVariants -R ${REF} -V combined_sorted.vcf.gz --select-type-to-include SNP -O cattle_raw_snps.vcf.gz --tmp-dir tmp 

gatk --java-options "-Xmx28G -XX:+UseParallelGC -XX:ParallelGCThreads=28" SelectVariants -R ${REF} -V combined_sorted.vcf.gz --select-type-to-include INDEL -O cattle_raw_indels.vcf.gz --tmp-dir tmp 

##Plotting distribution of each parameter:Extract parameter values then plot in R
vcftools --gzvcf cattle_raw_snps.vcf.gz --out cattle_snps_MQ --get-INFO MQ 

vcftools --gzvcf cattle_raw_indels.vcf.gz --out cattle_indels_MQ --get-INFO MQ 

#Variant Filtration 
gatk --java-options "-Xmx28G -XX:+UseParallelGC -XX:ParallelGCThreads=28"  VariantFiltration \
-R ${REF} -V cattle_raw_snps.vcf.gz -O cattle_filtered_snps.vcf.gz --cluster-size 3 --cluster-window-size 10 \
--filter "QD < 2.0" --filter-name "QD2" \
--filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
--filter "QUAL < 30.0" --filter-name "QUAL30" \
--filter "FS > 60.0" --filter-name "FS60" \
--filter "MQ < 40.0" --filter-name "MQ40" \ 
--filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
--filter "SOR > 3.0" --filter-name "SOR3" \
--tmp-dir ../tmp 

#Filter out biallelic loci
vcftools --gzvcf cattle_filtered_snps.vcf.gz --min-alleles 2 --max-alleles 2 --recode --stdout | ./bgzip -c > cattle_biallelic_filtered.vcf.gz
