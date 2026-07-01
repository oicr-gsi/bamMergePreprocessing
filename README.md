# bamMergePreprocessing

Workflow to merge and preprocess lane level alignments.

## Overview

## Dependencies

* [samtools 1.15](http://www.htslib.org/)
* [gatk 4.1.6.0](https://gatk.broadinstitute.org)


## Usage

### Cromwell
```
java -jar cromwell.jar run bamMergePreprocessing.wdl --inputs inputs.json
```

### Inputs

#### Required workflow parameters:
Parameter|Value|Description
---|---|---
`inputBamFiles`|Array[bamFiles]|Array of objects describing sets of bams to merge together and the merged file name. These merged bams will be cocleaned together and output separately (by merged name).
`outputFileNamePrefix`|String|Prefix of output file name
`reference`|String|Path to reference file.
`referenceGenome`|String|The reference genome version for input sample


#### Optional workflow parameters:
Parameter|Value|Default|Description
---|---|---|---
`intervalsToParallelizeByString`|String|"chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY,chrM,NC,SPLIT,UNALIGNED"|Comma separated list of intervals to split by (e.g. chr1,chr2,chr3+chr4).
`doFilter`|Boolean|true|Enable/disable Samtools filtering.
`doMarkDuplicates`|Boolean|true|Enable/disable GATK4 MarkDuplicates.
`doBqsr`|Boolean|true|Enable/disable GATK baseQualityScoreRecalibration
`provisionBqsr`|Boolean|false|Enable/disable provision out bqsr report and table
`libType`|String|"dna"|Sequencing library type, e.g. 'dna' or 'rna'
`doBamMetrics`|Boolean|false|Enable/disable generation of bam metrics (samtools stats/flagstat/counts) at each processing stage.


#### Optional task parameters:
Parameter|Value|Default|Description
---|---|---|---
`prepareIntervals.lineSeparator`|String|","|Interval group separator - these are the intervals to split by.
`prepareIntervals.recordSeparator`|String|"+"|Interval interval group separator - this can be used to combine multiple intervals into one group.
`prepareIntervals.jobMemory`|Int|1|Memory allocated to job (in GB).
`prepareIntervals.cores`|Int|1|The number of cores to allocate to the job.
`prepareIntervals.timeout`|Int|1|Maximum amount of time (in hours) the task can run for.
`prepareIntervals.modules`|String|""|Environment module name and version to load (space separated) before command execution.
`inputBamMetrics.jobMemory`|Int|12|Memory allocated to job (in GB).
`inputBamMetrics.cores`|Int|1|The number of cores to allocate to the job.
`inputBamMetrics.timeout`|Int|6|Maximum amount of time (in hours) the task can run for.
`inputBamMetrics.modules`|String|"samtools/1.15"|Environment module name and version to load.
`subsetAndFilter.temporaryWorkingDir`|String|""|Where to write out intermediary bam files. Only the final preprocessed bam will be written to task working directory if this is set to local tmp.
`subsetAndFilter.filterSuffix`|String|".filtered"|Suffix to append to output file name when filtering is applied.
`subsetAndFilter.filterFlags`|Int|256|Samtools filter flags to apply.
`subsetAndFilter.minMapQuality`|Int?|None|Samtools minimum mapping quality filter to apply.
`subsetAndFilter.filterAdditionalParams`|String?|None|Additional parameters to pass to samtools.
`subsetAndFilter.oldStyle`|Boolean|false|Hidden option to revert to the previous dupMarking strategy, for assessment purposes.
`subsetAndFilter.jobMemory`|Int|36|Memory allocated to job (in GB).
`subsetAndFilter.minMemory`|Int|12|A minimum amount of memory allocated to the task, overrides the scaled RAM setting
`subsetAndFilter.overhead`|Int|8|Java overhead memory (in GB). jobMemory - overhead == java Xmx/heap memory.
`subsetAndFilter.cores`|Int|1|The number of cores to allocate to the job.
`subsetAndFilter.timeout`|Int|6|Maximum amount of time (in hours) the task can run for.
`subsetAndFilter.modules`|String|"samtools/1.15"|Environment module name and version to load (space separated) before command execution.
`baseQualityScoreRecalibration.intervals`|Array[String]|[]|One or more genomic intervals over which to operate.
`baseQualityScoreRecalibration.additionalParams`|String?|None|Additional parameters to pass to GATK BaseRecalibrator.
`baseQualityScoreRecalibration.jobMemory`|Int|24|Memory allocated to job (in GB).
`baseQualityScoreRecalibration.overhead`|Int|6|Java overhead memory (in GB). jobMemory - overhead == java Xmx/heap memory.
`baseQualityScoreRecalibration.cores`|Int|1|The number of cores to allocate to the job.
`baseQualityScoreRecalibration.timeout`|Int|6|Maximum amount of time (in hours) the task can run for.
`baseQualityScoreRecalibration.modules`|String|"gatk/4.1.6.0"|Environment module name and version to load (space separated) before command execution.
`markDuplicates.removeDuplicates`|Boolean|false|MarkDuplicates remove duplicates?
`markDuplicates.opticalDuplicatePixelDistance`|Int|100|MarkDuplicates optical distance.
`markDuplicates.markDuplicatesAdditionalParams`|String?|None|Additional parameters to pass to GATK MarkDuplicates.
`markDuplicates.jobMemory`|Int|36|Memory allocated to job (in GB).
`markDuplicates.minMemory`|Int|12|A minimum amount of memory allocated to the task, overrides the scaled RAM setting
`markDuplicates.overhead`|Int|6|Java overhead memory (in GB). jobMemory - overhead == java Xmx/heap memory.
`markDuplicates.cores`|Int|1|The number of cores to allocate to the job.
`markDuplicates.timeout`|Int|6|Maximum amount of time (in hours) the task can run for.
`markDuplicates.modules`|String|"gatk/4.1.6.0 samtools/1.15"|Environment module name and version to load (space separated) before command execution.
`mergeWithinInterval.additionalParams`|String?|None|Additional parameters to pass to GATK MergeSamFiles.
`mergeWithinInterval.comment`|String|""|Comment to add to the header of the merged bam file.
`mergeWithinInterval.jobMemory`|Int|24|Memory allocated to job (in GB).
`mergeWithinInterval.overhead`|Int|6|Java overhead memory (in GB). jobMemory - overhead == java Xmx/heap memory.
`mergeWithinInterval.cores`|Int|1|The number of cores to allocate to the job.
`mergeWithinInterval.timeout`|Int|6|Maximum amount of time (in hours) the task can run for.
`mergeWithinInterval.modules`|String|"gatk/4.1.6.0 samtools/1.15"|Environment module name and version to load (space separated) before command execution.
`splitNCigarString.refactorCigarString`|String|false|SplitNCigarReads refactor-cigar-string option.
`splitNCigarString.splitNCigarReadsAdditionalParams`|String?|None|Additional parameters to pass to GATK SplitNCigarReads.
`splitNCigarString.readFilters`|String?|None|Optional GATK read filters to apply.
`splitNCigarString.jobMemory`|Int|24|Memory allocated to job (in GB).
`splitNCigarString.overhead`|Int|6|Java overhead memory (in GB). jobMemory - overhead == java Xmx/heap memory.
`splitNCigarString.cores`|Int|1|The number of cores to allocate to the job.
`splitNCigarString.timeout`|Int|6|Maximum amount of time (in hours) the task can run for.
`splitNCigarString.modules`|String|"gatk/4.1.6.0 samtools/1.15"|Environment module name and version to load (space separated) before command execution.
`subsetAndFilterBamMetrics.jobMemory`|Int|12|Memory allocated to job (in GB).
`subsetAndFilterBamMetrics.cores`|Int|1|The number of cores to allocate to the job.
`subsetAndFilterBamMetrics.timeout`|Int|6|Maximum amount of time (in hours) the task can run for.
`subsetAndFilterBamMetrics.modules`|String|"samtools/1.15"|Environment module name and version to load.
`duplicateMarkedBamMetrics.jobMemory`|Int|12|Memory allocated to job (in GB).
`duplicateMarkedBamMetrics.cores`|Int|1|The number of cores to allocate to the job.
`duplicateMarkedBamMetrics.timeout`|Int|6|Maximum amount of time (in hours) the task can run for.
`duplicateMarkedBamMetrics.modules`|String|"samtools/1.15"|Environment module name and version to load.
`mergedWithinIntervalBamMetrics.jobMemory`|Int|12|Memory allocated to job (in GB).
`mergedWithinIntervalBamMetrics.cores`|Int|1|The number of cores to allocate to the job.
`mergedWithinIntervalBamMetrics.timeout`|Int|6|Maximum amount of time (in hours) the task can run for.
`mergedWithinIntervalBamMetrics.modules`|String|"samtools/1.15"|Environment module name and version to load.
`splitNCigarStringBamMetrics.jobMemory`|Int|12|Memory allocated to job (in GB).
`splitNCigarStringBamMetrics.cores`|Int|1|The number of cores to allocate to the job.
`splitNCigarStringBamMetrics.timeout`|Int|6|Maximum amount of time (in hours) the task can run for.
`splitNCigarStringBamMetrics.modules`|String|"samtools/1.15"|Environment module name and version to load.
`gatherBQSRReports.additionalParams`|String?|None|Additional parameters to pass to GATK GatherBQSRReports.
`gatherBQSRReports.jobMemory`|Int|24|Memory allocated to job (in GB).
`gatherBQSRReports.overhead`|Int|6|Java overhead memory (in GB). jobMemory - overhead == java Xmx/heap memory.
`gatherBQSRReports.cores`|Int|1|The number of cores to allocate to the job.
`gatherBQSRReports.timeout`|Int|6|Maximum amount of time (in hours) the task can run for.
`gatherBQSRReports.modules`|String|"gatk/4.1.6.0"|Environment module name and version to load (space separated) before command execution.
`analyzeCovariates.additionalParams`|String?|None|Additional parameters to pass to GATK AnalyzeCovariates
`analyzeCovariates.jobMemory`|Int|24|Memory allocated to job (in GB).
`analyzeCovariates.overhead`|Int|6|Java overhead memory (in GB). jobMemory - overhead == java Xmx/heap memory.
`analyzeCovariates.cores`|Int|1|The number of cores to allocate to the job.
`analyzeCovariates.timeout`|Int|6|Maximum amount of time (in hours) the task can run for.
`analyzeCovariates.modules`|String|"gatk/4.1.6.0"|Environment module name and version to load (space separated) before command execution.
`applyBaseQualityScoreRecalibration.additionalParams`|String?|None|Additional parameters to pass to GATK ApplyBQSR.
`applyBaseQualityScoreRecalibration.jobMemory`|Int|24|Memory allocated to job (in GB).
`applyBaseQualityScoreRecalibration.overhead`|Int|6|Java overhead memory (in GB). jobMemory - overhead == java Xmx/heap memory.
`applyBaseQualityScoreRecalibration.cores`|Int|1|The number of cores to allocate to the job.
`applyBaseQualityScoreRecalibration.timeout`|Int|6|Maximum amount of time (in hours) the task can run for.
`applyBaseQualityScoreRecalibration.modules`|String|"gatk/4.1.6.0 samtools/1.15"|Environment module name and version to load (space separated) before command execution.
`recalibratedBamMetrics.jobMemory`|Int|12|Memory allocated to job (in GB).
`recalibratedBamMetrics.cores`|Int|1|The number of cores to allocate to the job.
`recalibratedBamMetrics.timeout`|Int|6|Maximum amount of time (in hours) the task can run for.
`recalibratedBamMetrics.modules`|String|"samtools/1.15"|Environment module name and version to load.
`mergeAcrossIntervals.additionalParams`|String?|None|Additional parameters to pass to GATK MergeSamFiles.
`mergeAcrossIntervals.jobMemory`|Int|24|Memory allocated to job (in GB).
`mergeAcrossIntervals.overhead`|Int|6|Java overhead memory (in GB). jobMemory - overhead == java Xmx/heap memory.
`mergeAcrossIntervals.cores`|Int|1|The number of cores to allocate to the job.
`mergeAcrossIntervals.timeout`|Int|6|Maximum amount of time (in hours) the task can run for.
`mergeAcrossIntervals.modules`|String|"gatk/4.1.6.0 samtools/1.15"|Environment module name and version to load (space separated) before command execution.
`finalBamMetrics.jobMemory`|Int|12|Memory allocated to job (in GB).
`finalBamMetrics.cores`|Int|1|The number of cores to allocate to the job.
`finalBamMetrics.timeout`|Int|6|Maximum amount of time (in hours) the task can run for.
`finalBamMetrics.modules`|String|"samtools/1.15"|Environment module name and version to load.
`zippedBamMetrics.jobMemory`|Int|4|Memory allocated to job (in GB).
`zippedBamMetrics.cores`|Int|1|The number of cores to allocate to the job.
`zippedBamMetrics.timeout`|Int|1|Maximum amount of time (in hours) the task can run for.
`zippedDupMarkMetrics.jobMemory`|Int|4|Memory allocated to job (in GB).
`zippedDupMarkMetrics.cores`|Int|1|The number of cores to allocate to the job.
`zippedDupMarkMetrics.timeout`|Int|1|Maximum amount of time (in hours) the task can run for.


### Outputs

Output | Type | Description | Labels
---|---|---|---
`mergedBam`|File|the final merged bam.|
`mergedBamIndex`|File|the final merged bam index|
`recalibrationReport`|File?|Recalibration report pdf (if BQSR enabled).|
`recalibrationTable`|File?|Recalibration csv that was used by BQSR (if BQSR enabled).|
`markDuplicateMetricsZip`|File?|A tarball of markDuplicates metrics files across all intervals (if doMarkDuplicates enabled).|
`bamMetricsZip`|File?|A tarball of samtools stats/flagstat/mapcounts metrics files generated at each processing stage (if doBamMetrics enabled).|


./commands.txt found, printing out the content...
## Commands
 
 This section lists command(s) run by bamMergePreprocessing
 
 * Running bamMergePreprocessing workflow
 
 ### Running samstats on input bam files
 
 ```
     set -euo pipefail
     ### pipe the SAM file through flagstat (tee), stats (tee) and a counting operation
     samtools view -h ~{inputBam} |  tee >(samtools flagstat - > ~{prefix}.flagstats.txt) | tee >(samtools stats - > ~{prefix}.samstats.txt) | samtools view - | cut -f 2,3,7 | sort | uniq -c > ~{prefix}.mapcounts.txt
 ```
 
 ### Preparing intervals 
 
 ```
     set -euo pipefail
 
     ### intervals are separated by line or record separator, 
     echo "~{str}" | tr '~{lineSeparator}' '\n' | tr '~{recordSeparator}' '\t' > intervals
     
     ## this will generate a list of chrosomes or keywords in the intervals, removing any position information
     cat intervals | sed 's/\t/\n/g' | sed 's/:.*//' | sort -u > interval_contigs
     
     ### create a bed file from all contigs in the reference
     cat ~{refFai} | awk -v OFS="\t" '{ print $1, 1, $2 }' > contigs.bed
     
     ### create a file with the allowed keywords 
     echo -e "NC\nSPLIT\nUNALIGNED" > keywords
 
     ### are there any contigs in the supplied intervals that are NOT in the reference build.  if so, this should rais an concern
     #cat interval_contigs | grep -v -f <(cut -f 1 contigs.bed) | grep -v -f keywords > unknown_contigs
 
     ### nc.contigs.bed includes intervals NOT in the interval_contigs.
     ### this is returned by the task, and is used to subset the bam file with samtools view -L on the NC interval
     #cat contigs.bed | grep -vFw -f interval_contigs > nc.contigs.bed
     cat contigs.bed | grep "_" > nc.contigs.bed
     
 
     ### this is now being read from a file, instead of from stdout
     ####echo "~{str}" | tr '~{lineSeparator}' '\n' | tr '~{recordSeparator}' '\t'
 
     #### the python code block will read in the intervals and determine the size of each based on the contigs
     python3 <<CODE
     import re
     contigs={}
 
     total=0
     with open("contigs.bed","r") as contigbed:
         for line in contigbed:
             contig,start,end=line.strip().split("\t")
             contigs[contig]=int(end)-int(start)+1
             total=total + contigs[contig]
 
     cout=open("coefficients.txt","w")
   
     with open("intervals","r") as interval_set:
         for line in interval_set:
             intervals=line.strip().split(" ")
             interval_size=0
             for interval in intervals:
                 if ":" in interval:
                     contig,start,end=re.split(r'[:-]',interval)
                     size=int(end)-int(start)+1
                     interval_size=interval_size+size
                 elif interval == "NC":
                     with open("nc.contigs.bed","r") as ncbed:
                         for ncline in ncbed:
                             nc_contig,start,end=ncline.strip().split("\t")
                             interval_size=interval_size + int(end)-int(start)+1
                 else:
                     ## the interval should be a full contig
                     size=contigs.get(interval,0)
                     interval_size=interval_size+size
             coeff=interval_size/total
             cout.write(line.strip() + "\t" + str(coeff) + "\n")
     cout.close()
     CODE
 ```
 
 ### Subsetting and filtering bam files, given an interval
 
 ```
     set -euxo pipefail
     
     ### write to local file
     #### dev fix to get rid of canonical chromosome
 
     #cat ~{ncBed} > nc.bed
     #sleep 10
 
     samtools view -b ~{exprString} ~{filterString} ~{intervalsString} ~{inputBam} ~{samtoolsInterval} > ~{outputFileNamePrefix}.bam
     samtools index ~{outputFileNamePrefix}.bam ~{outputFileNamePrefix}.bai
 ```
 
 ### Marking Duplicates
 
 ```
     set -euo pipefail
     gatk --java-options "-Xmx~{jobMemory - overhead}G" MarkDuplicates \
     ~{sep=" " prefix("--INPUT=", inputBams)}  \
     --OUTPUT ~{outputFileNamePrefix}.bam \
     --METRICS_FILE="~{outputFileNamePrefix}.metrics" \
     --VALIDATION_STRINGENCY=SILENT \
     --REMOVE_DUPLICATES=~{removeDuplicates} \
     --OPTICAL_DUPLICATE_PIXEL_DISTANCE=~{opticalDuplicatePixelDistance} \
     --CREATE_INDEX=true \
     ~{markDuplicatesAdditionalParams}
 
 ```
 
 ### Splitting reads based on cigar string
 
 ```
     set -euo pipefail
 
     gatk --java-options "-Xmx~{jobMemory - overhead}G" SplitNCigarReads \
     --input ~{inputBam}  \
     --output ~{outputFileNamePrefix}.bam \
     --reference ~{reference} \
     --create-output-bam-index true \
     --refactor-cigar-string ~{refactorCigarString} \
     ~{splitNCigarReadsAdditionalParams}
 
 ```
 
 ### Creating a tarball
 
 ```
     set -euo pipefail
     mkdir ./files/
     files="~{sep="," inputFiles}"    
     IFS=',' read -ra f <<< "$files"
     for f in ${f[@]}
     do
       cp $f ./files/
     done
     tar czf  ~{zipName}.tar.gz ./files/*
 
 ```
 
 
 ### Merging bam files
 
 ```
     set -euo pipefail
 
     gatk --java-options "-Xmx~{jobMemory - overhead}G" MergeSamFiles \
     ~{sep=" " prefix("--INPUT=", bams)} \
     --OUTPUT="~{outputFileNamePrefix}.bam" \
     --CREATE_INDEX=true \
     --SORT_ORDER=coordinate \
     --ASSUME_SORTED=false \
     --USE_THREADING=true \
     --COMMENT="'~{comment}'" \
     --VALIDATION_STRINGENCY=SILENT \
     ~{additionalParams}
 ```
 
 ### Base Recalibration
 
 ```
     set -euo pipefail
 
     gatk --java-options "-Xmx~{jobMemory - overhead}G" BaseRecalibrator \
     --reference ~{reference} \
     ~{sep=" " prefixedIntervals} \
     ~{sep=" " prefix("--input=", bams)} \
     ~{sep=" " prefix("--known-sites ", knownSites)} \
     --output=gatk.recalibration.csv \
     ~{additionalParams}
 ```
 
 ### Gathering Base Quality Score Recalibration Reports
 
 ```
     set -euo pipefail
 
     gatk --java-options "-Xmx~{jobMemory - overhead}G" GatherBQSRReports \
         ~{sep=" " prefix("--input=", recalibrationTables)} \
         --output ~{outputFileNamePrefix}.gatk.recalibration.csv \
         ~{additionalParams}
 ```
 
 ### Analysis of Covariates
 
 ```
     set -euo pipefail
 
     gatk --java-options "-Xmx~{jobMemory - overhead}G" AnalyzeCovariates \
     --bqsr-recal-file=~{recalibrationTable} \
     --plots-report-file ~{outputFileNamePrefix}.gatk.recalibration.pdf \
     ~{additionalParams}
 ```
 
 ### Appliying Base Quality Score Recalibration
 
 ```
     set -euo pipefail
 
     gatk --java-options "-Xmx~{jobMemory - overhead}G" ApplyBQSR \
     --bqsr-recal-file=~{recalibrationTable} \
     ~{sep=" " prefix("--input=", [bam])} \
     --output ~{outputFileNamePrefix}.bam \
     ~{additionalParams}
 
 ```
 
 ## Support

For support, please file an issue on the [Github project](https://github.com/oicr-gsi) or send an email to gsi@oicr.on.ca .

_Generated with generate-markdown-readme (https://github.com/oicr-gsi/gsi-wdl-tools/)_