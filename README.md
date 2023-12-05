# bamMergePreprocessing



## Overview

## Dependencies

* [samtools 1.9](http://www.htslib.org/)
* [gatk 4.1.6.0](https://gatk.broadinstitute.org)
* [gatk 3.6-0](https://gatk.broadinstitute.org)


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
`intervalsToParallelizeByString`|String|Comma separated list of intervals to split by (e.g. chr1,chr2,chr3+chr4).
`reference`|String|Path to reference file.
`referenceGenome`|String|The reference genome version for input sample


#### Optional workflow parameters:
Parameter|Value|Default|Description
---|---|---|---
`doFilter`|Boolean|true|Enable/disable Samtools filtering.
`doMarkDuplicates`|Boolean|true|Enable/disable GATK4 MarkDuplicates.
`doBqsr`|Boolean|false|Enable/disable GATK baseQualityScoreRecalibration


#### Optional task parameters:
Parameter|Value|Default|Description
---|---|---|---
`splitStringToArray.lineSeparator`|String|","|Interval group separator - these are the intervals to split by.
`splitStringToArray.recordSeparator`|String|"+"|Interval interval group separator - this can be used to combine multiple intervals into one group.
`splitStringToArray.jobMemory`|Int|1|Memory allocated to job (in GB).
`splitStringToArray.cores`|Int|1|The number of cores to allocate to the job.
`splitStringToArray.timeout`|Int|1|Maximum amount of time (in hours) the task can run for.
`splitStringToArray.modules`|String|""|Environment module name and version to load (space separated) before command execution.
`preprocessBam.temporaryWorkingDir`|String|""|Where to write out intermediary bam files. Only the final preprocessed bam will be written to task working directory if this is set to local tmp.
`preprocessBam.filterSuffix`|String|".filtered"|Suffix to use for filtered bams.
`preprocessBam.filterFlags`|Int|260|Samtools filter flags to apply.
`preprocessBam.minMapQuality`|Int?|None|Samtools minimum mapping quality filter to apply.
`preprocessBam.filterAdditionalParams`|String?|None|Additional parameters to pass to samtools.
`preprocessBam.dedupSuffix`|String|".deduped"|Suffix to use for markDuplcated bams
`preprocessBam.removeDuplicates`|Boolean|false|MarkDuplicates remove duplicates?
`preprocessBam.opticalDuplicatePixelDistance`|Int|100|MarkDuplicates optical distance.
`preprocessBam.markDuplicatesAdditionalParams`|String?|None|Additional parameters to pass to GATK MarkDuplicates.
`preprocessBam.jobMemory`|Int|48|Memory allocated to job (in GB).
`preprocessBam.overhead`|Int|8|Java overhead memory (in GB). jobMemory - overhead == java Xmx/heap memory.
`preprocessBam.cores`|Int|1|The number of cores to allocate to the job.
`preprocessBam.timeout`|Int|6|Maximum amount of time (in hours) the task can run for.
`preprocessBam.modules`|String|"samtools/1.9 gatk/4.1.6.0"|Environment module name and version to load (space separated) before command execution.
`filterBam.temporaryWorkingDir`|String|""|Where to write out intermediary bam files. Only the final preprocessed bam will be written to task working directory if this is set to local tmp.
`filterBam.filterSuffix`|String|".filtered"|Suffix to use for filtered bams.
`filterBam.filterFlags`|Int|260|Samtools filter flags to apply.
`filterBam.minMapQuality`|Int?|None|Samtools minimum mapping quality filter to apply.
`filterBam.filterAdditionalParams`|String?|None|Additional parameters to pass to samtools.
`filterBam.jobMemory`|Int|48|Memory allocated to job (in GB).
`filterBam.overhead`|Int|8|Java overhead memory (in GB). jobMemory - overhead == java Xmx/heap memory.
`filterBam.cores`|Int|1|The number of cores to allocate to the job.
`filterBam.timeout`|Int|6|Maximum amount of time (in hours) the task can run for.
`filterBam.modules`|String|"samtools/1.9 gatk/4.1.6.0"|Environment module name and version to load (space separated) before command execution.
`markDuplicates.removeDuplicates`|Boolean|false|MarkDuplicates remove duplicates?
`markDuplicates.opticalDuplicatePixelDistance`|Int|100|MarkDuplicates optical distance.
`markDuplicates.markDuplicatesAdditionalParams`|String?|None|Additional parameters to pass to GATK MarkDuplicates.
`markDuplicates.jobMemory`|Int|24|Memory allocated to job (in GB).
`markDuplicates.overhead`|Int|6|Java overhead memory (in GB). jobMemory - overhead == java Xmx/heap memory.
`markDuplicates.cores`|Int|1|The number of cores to allocate to the job.
`markDuplicates.timeout`|Int|6|Maximum amount of time (in hours) the task can run for.
`markDuplicates.modules`|String|"gatk/4.1.6.0"|Environment module name and version to load (space separated) before command execution.
`markDuplicates.dedupSuffix`|String|".deduped"|Suffix to use for markDuplcated bams
`mergeMultipleBam.baseName`|String|basename(bams[0])|The base name for output files
`mergeMultipleBam.additionalParams`|String?|None|Additional parameters to pass to GATK MergeSamFiles.
`mergeMultipleBam.jobMemory`|Int|24|Memory allocated to job (in GB).
`mergeMultipleBam.overhead`|Int|6|Java overhead memory (in GB). jobMemory - overhead == java Xmx/heap memory.
`mergeMultipleBam.cores`|Int|1|The number of cores to allocate to the job.
`mergeMultipleBam.timeout`|Int|6|Maximum amount of time (in hours) the task can run for.
`mergeMultipleBam.modules`|String|"gatk/4.1.6.0"|Environment module name and version to load (space separated) before command execution.
`baseQualityScoreRecalibration.intervals`|Array[String]|[]|One or more genomic intervals over which to operate.
`baseQualityScoreRecalibration.additionalParams`|String?|None|Additional parameters to pass to GATK BaseRecalibrator.
`baseQualityScoreRecalibration.outputFileName`|String|"gatk.recalibration.csv"|Recalibration table file name.
`baseQualityScoreRecalibration.jobMemory`|Int|24|Memory allocated to job (in GB).
`baseQualityScoreRecalibration.overhead`|Int|6|Java overhead memory (in GB). jobMemory - overhead == java Xmx/heap memory.
`baseQualityScoreRecalibration.cores`|Int|1|The number of cores to allocate to the job.
`baseQualityScoreRecalibration.timeout`|Int|6|Maximum amount of time (in hours) the task can run for.
`baseQualityScoreRecalibration.modules`|String|"gatk/4.1.6.0"|Environment module name and version to load (space separated) before command execution.
`gatherBQSRReports.additionalParams`|String?|None|Additional parameters to pass to GATK GatherBQSRReports.
`gatherBQSRReports.outputFileName`|String|"gatk.recalibration.csv"|Recalibration table file name.
`gatherBQSRReports.jobMemory`|Int|24|Memory allocated to job (in GB).
`gatherBQSRReports.overhead`|Int|6|Java overhead memory (in GB). jobMemory - overhead == java Xmx/heap memory.
`gatherBQSRReports.cores`|Int|1|The number of cores to allocate to the job.
`gatherBQSRReports.timeout`|Int|6|Maximum amount of time (in hours) the task can run for.
`gatherBQSRReports.modules`|String|"gatk/4.1.6.0"|Environment module name and version to load (space separated) before command execution.
`analyzeCovariates.additionalParams`|String?|None|Additional parameters to pass to GATK AnalyzeCovariates
`analyzeCovariates.outputFileName`|String|"gatk.recalibration.pdf"|Recalibration report file name.
`analyzeCovariates.jobMemory`|Int|24|Memory allocated to job (in GB).
`analyzeCovariates.overhead`|Int|6|Java overhead memory (in GB). jobMemory - overhead == java Xmx/heap memory.
`analyzeCovariates.cores`|Int|1|The number of cores to allocate to the job.
`analyzeCovariates.timeout`|Int|6|Maximum amount of time (in hours) the task can run for.
`analyzeCovariates.modules`|String|"gatk/4.1.6.0"|Environment module name and version to load (space separated) before command execution.
`applyBaseQualityScoreRecalibration.outputFileName`|String|basename(bam,".bam")|Output files will be prefixed with this.
`applyBaseQualityScoreRecalibration.suffix`|String|".recalibrated"|Suffix to use for recalibrated bams.
`applyBaseQualityScoreRecalibration.additionalParams`|String?|None|Additional parameters to pass to GATK ApplyBQSR.
`applyBaseQualityScoreRecalibration.jobMemory`|Int|24|Memory allocated to job (in GB).
`applyBaseQualityScoreRecalibration.overhead`|Int|6|Java overhead memory (in GB). jobMemory - overhead == java Xmx/heap memory.
`applyBaseQualityScoreRecalibration.cores`|Int|1|The number of cores to allocate to the job.
`applyBaseQualityScoreRecalibration.timeout`|Int|6|Maximum amount of time (in hours) the task can run for.
`applyBaseQualityScoreRecalibration.modules`|String|"gatk/4.1.6.0"|Environment module name and version to load (space separated) before command execution.
`mergeBams.baseName`|String|basename(bams[0])|The base name for output files
`mergeBams.additionalParams`|String?|None|Additional parameters to pass to GATK MergeSamFiles.
`mergeBams.jobMemory`|Int|24|Memory allocated to job (in GB).
`mergeBams.overhead`|Int|6|Java overhead memory (in GB). jobMemory - overhead == java Xmx/heap memory.
`mergeBams.cores`|Int|1|The number of cores to allocate to the job.
`mergeBams.timeout`|Int|6|Maximum amount of time (in hours) the task can run for.
`mergeBams.modules`|String|"gatk/4.1.6.0"|Environment module name and version to load (space separated) before command execution.


### Outputs

Output | Type | Description
---|---|---
`mergedBam`|File|the final merged bam.
`mergedBamIndex`|File|the final merged bam index
`recalibrationReport`|File?|Recalibration report pdf (if BQSR enabled).
`recalibrationTable`|File?|Recalibration csv that was used by BQSR (if BQSR enabled).


## Commands
 This section lists command(s) run by bamMergePreprocessing workflow
 
 * Running bamMergePreprocessing
 
 
 <<<
     set -euo pipefail
 
     echo "~{str}" | tr '~{lineSeparator}' '\n' | tr '~{recordSeparator}' '\t'
   >>>
 <<<
     set -euxo pipefail
 
     # filter
     if [ "~{doFilter}" = true ]; then
       outputBam="~{workingDir}~{baseFileName}~{filterSuffix}.bam"
       outputBamIndex="~{workingDir}~{baseFileName}~{filterSuffix}.bai"
       samtools view -b \
       -F ~{filterFlags} \
       ~{"-q " + minMapQuality} \
       ~{filterAdditionalParams} \
       ~{inputBam} \
       ~{sep=" " intervals} > $outputBam
       samtools index $outputBam $outputBamIndex
 
       # set inputs for next step
       inputBam=$outputBam
       inputBamIndex=$outputBamIndex
     else
       outputBam="~{workingDir}~{baseFileName}.bam"
       outputBamIndex="~{workingDir}~{baseFileName}.bai"
       samtools view -b \
       ~{inputBam} \
       ~{sep=" " intervals} > $outputBam
       samtools index $outputBam $outputBamIndex
 
       # set inputs for next step
       inputBam=$outputBam
       inputBamIndex=$outputBamIndex
     fi
 
     # mark duplicates
     if [ "~{doMarkDuplicates}" = true ]; then
       gatk --java-options "-Xmx~{jobMemory - overhead}G" MarkDuplicates \
       --INPUT=$inputBam  \
       --OUTPUT="~{markDuplicatesFileName}.bam" \
       --METRICS_FILE="~{outputFileName}.metrics" \
       --VALIDATION_STRINGENCY=SILENT \
       --REMOVE_DUPLICATES=~{removeDuplicates} \
       --OPTICAL_DUPLICATE_PIXEL_DISTANCE=~{opticalDuplicatePixelDistance} \
       --CREATE_INDEX=true \
       ~{markDuplicatesAdditionalParams}
     fi
 
   >>>
 <<<
     set -euxo pipefail
 
     # filter
     if [ "~{doFilter}" = true ]; then
       outputBam="~{workingDir}~{baseFileName}~{filterSuffix}.bam"
       outputBamIndex="~{workingDir}~{baseFileName}~{filterSuffix}.bai"
       samtools view -b \
       -F ~{filterFlags} \
       ~{"-q " + minMapQuality} \
       ~{filterAdditionalParams} \
       ~{inputBam} \
       ~{sep=" " intervals} > $outputBam
       samtools index $outputBam $outputBamIndex
 
       # set inputs for next step
       inputBam=$outputBam
       inputBamIndex=$outputBamIndex
     else
       outputBam="~{workingDir}~{baseFileName}.bam"
       outputBamIndex="~{workingDir}~{baseFileName}.bai"
       samtools view -b \
       ~{inputBam} \
       ~{sep=" " intervals} > $outputBam
       samtools index $outputBam $outputBamIndex
     fi
   >>>
 <<<
     set -euo pipefail
       gatk --java-options "-Xmx~{jobMemory - overhead}G" MarkDuplicates \
       ~{sep=" " prefix("--INPUT=", inputBams)}  \
       --OUTPUT ~{outputFileName}~{dedupSuffix}.bam \
       --METRICS_FILE="~{outputFileName}~{dedupSuffix}.metrics" \
       --VALIDATION_STRINGENCY=SILENT \
       --REMOVE_DUPLICATES=~{removeDuplicates} \
       --OPTICAL_DUPLICATE_PIXEL_DISTANCE=~{opticalDuplicatePixelDistance} \
       --CREATE_INDEX=true \
       ~{markDuplicatesAdditionalParams}
     >>>
 <<<
     set -euo pipefail
     baseName=~{baseName}
     outputBamSuffix="${baseName#*.}"
     gatk --java-options "-Xmx~{jobMemory - overhead}G" MergeSamFiles \
     ~{sep=" " prefix("--INPUT=", bams)} \
     --OUTPUT="~{outputFileName}.$outputBamSuffix" \
     --CREATE_INDEX=true \
     --SORT_ORDER=coordinate \
     --ASSUME_SORTED=false \
     --USE_THREADING=true \
     --VALIDATION_STRINGENCY=SILENT \
     ~{additionalParams}
   >>>
 <<<
     set -euo pipefail
 
     gatk --java-options "-Xmx~{jobMemory - overhead}G" BaseRecalibrator \
     --reference ~{reference} \
     ~{sep=" " prefixedIntervals} \
     ~{sep=" " prefix("--input=", bams)} \
     ~{sep=" " prefix("--known-sites ", knownSites)} \
     --output=~{outputFileName} \
     ~{additionalParams}
   >>>
 <<<
     set -euo pipefail
 
     gatk --java-options "-Xmx~{jobMemory - overhead}G" GatherBQSRReports \
     --input ~{recalibrationTables} \
     --output ~{outputFileName} \
     ~{additionalParams}
   >>>
 <<<
     set -euo pipefail
 
     gatk --java-options "-Xmx~{jobMemory - overhead}G" AnalyzeCovariates \
     --bqsr-recal-file=~{recalibrationTable} \
     --plots-report-file ~{outputFileName} \
     ~{additionalParams}
   >>>
 <<<
     set -euo pipefail
 
     gatk --java-options "-Xmx~{jobMemory - overhead}G" ApplyBQSR \
     --bqsr-recal-file=~{recalibrationTable} \
     ~{sep=" " prefix("--input=", [bam])} \
     --output ~{outputFileName}~{suffix}.bam \
     ~{additionalParams}
   >>>
 ## Support

For support, please file an issue on the [Github project](https://github.com/oicr-gsi) or send an email to gsi@oicr.on.ca .

_Generated with generate-markdown-readme (https://github.com/oicr-gsi/gsi-wdl-tools/)_
