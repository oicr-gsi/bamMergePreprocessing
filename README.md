# bamMergePreprocessing



## Overview

## Dependencies

* [samtools 1.9](http://www.htslib.org/)
* [gatk 4.1.6.0](https://gatk.broadinstitute.org)
* [gatk 3.6-0](https://gatk.broadinstitute.org)
* [python 3.7](https://www.python.org)


## Usage

### Cromwell
```
java -jar cromwell.jar run bamMergePreprocessing.wdl --inputs inputs.json
```

### Inputs

#### Required workflow parameters:
Parameter|Value|Description
---|---|---
`inputGroups`|Array[InputGroup]|Array of objects describing sets of bams to merge together and the merged file name. These merged bams will be cocleaned together and output separately (by merged name).
`intervalsToParallelizeByString`|String|Comma separated list of intervals to split by (e.g. chr1,chr2,chr3+chr4).
`reference`|String|Path to reference file.
`realignerTargetCreator.knownIndels`|Array[String]|Array of input VCF files with known indels.
`indelRealign.knownAlleles`|Array[String]|Array of input VCF files with known indels.
`baseQualityScoreRecalibration.knownSites`|Array[String]|Array of VCF with known polymorphic sites used to exclude regions around known polymorphisms from analysis.


#### Optional workflow parameters:
Parameter|Value|Default|Description
---|---|---|---
`doFilter`|Boolean|true|Enable/disable Samtools filtering.
`doMarkDuplicates`|Boolean|true|Enable/disable GATK4 MarkDuplicates.
`doSplitNCigarReads`|Boolean|false|Enable/disable GATK4 SplitNCigarReads.
`doIndelRealignment`|Boolean|true|Enable/disable GATK3 RealignerTargetCreator + IndelRealigner.
`doBqsr`|Boolean|true|Enable/disable GATK4 BQSR.


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
`preprocessBam.filterSuffix`|String|".filter"|Suffix to use for filtered bams.
`preprocessBam.filterFlags`|Int|260|Samtools filter flags to apply.
`preprocessBam.minMapQuality`|Int?|None|Samtools minimum mapping quality filter to apply.
`preprocessBam.filterAdditionalParams`|String?|None|Additional parameters to pass to samtools.
`preprocessBam.markDuplicatesSuffix`|String|".deduped"|Suffix to use for duplicate marked bams.
`preprocessBam.removeDuplicates`|Boolean|false|MarkDuplicates remove duplicates?
`preprocessBam.opticalDuplicatePixelDistance`|Int|100|MarkDuplicates optical distance.
`preprocessBam.markDuplicatesAdditionalParams`|String?|None|Additional parameters to pass to GATK MarkDuplicates.
`preprocessBam.splitNCigarReadsSuffix`|String|".split"|Suffix to use for SplitNCigarReads bams.
`preprocessBam.refactorCigarString`|Boolean|false|SplitNCigarReads refactor cigar string?
`preprocessBam.readFilters`|Array[String]|[]|SplitNCigarReads read filters
`preprocessBam.splitNCigarReadsAdditionalParams`|String?|None|Additional parameters to pass to GATK SplitNCigarReads.
`preprocessBam.jobMemory`|Int|24|Memory allocated to job (in GB).
`preprocessBam.overhead`|Int|6|Java overhead memory (in GB). jobMemory - overhead == java Xmx/heap memory.
`preprocessBam.cores`|Int|1|The number of cores to allocate to the job.
`preprocessBam.timeout`|Int|6|Maximum amount of time (in hours) the task can run for.
`preprocessBam.modules`|String|"samtools/1.9 gatk/4.1.6.0"|Environment module name and version to load (space separated) before command execution.
`realignerTargetCreator.downsamplingType`|String?|None|Type of read downsampling to employ at a given locus (NONE|ALL_READS|BY_SAMPLE).
`realignerTargetCreator.additionalParams`|String?|None|Additional parameters to pass to GATK RealignerTargetCreator.
`realignerTargetCreator.jobMemory`|Int|24|Memory allocated to job (in GB).
`realignerTargetCreator.overhead`|Int|6|Java overhead memory (in GB). jobMemory - overhead == java Xmx/heap memory.
`realignerTargetCreator.cores`|Int|1|The number of cores to allocate to the job.
`realignerTargetCreator.timeout`|Int|6|Maximum amount of time (in hours) the task can run for.
`realignerTargetCreator.modules`|String|"gatk/3.6-0"|Environment module name and version to load (space separated) before command execution.
`realignerTargetCreator.gatkJar`|String|"$GATK_ROOT/GenomeAnalysisTK.jar"|Path to GATK jar.
`indelRealign.additionalParams`|String?|None|Additional parameters to pass to GATK IndelRealigner.
`indelRealign.jobMemory`|Int|24|Memory allocated to job (in GB).
`indelRealign.overhead`|Int|6|Java overhead memory (in GB). jobMemory - overhead == java Xmx/heap memory.
`indelRealign.cores`|Int|1|The number of cores to allocate to the job.
`indelRealign.timeout`|Int|6|Maximum amount of time (in hours) the task can run for.
`indelRealign.modules`|String|"python/3.7 gatk/3.6-0"|Environment module name and version to load (space separated) before command execution.
`indelRealign.gatkJar`|String|"$GATK_ROOT/GenomeAnalysisTK.jar"|Path to GATK jar.
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
`collectFilesBySample.jobMemory`|Int|1|Memory allocated to job (in GB).
`collectFilesBySample.cores`|Int|1|The number of cores to allocate to the job.
`collectFilesBySample.timeout`|Int|1|Maximum amount of time (in hours) the task can run for.
`collectFilesBySample.modules`|String|"python/3.7"|Environment module name and version to load (space separated) before command execution.
`mergeSplitByIntervalBams.additionalParams`|String?|None|Additional parameters to pass to GATK MergeSamFiles.
`mergeSplitByIntervalBams.jobMemory`|Int|24|Memory allocated to job (in GB).
`mergeSplitByIntervalBams.overhead`|Int|6|Java overhead memory (in GB). jobMemory - overhead == java Xmx/heap memory.
`mergeSplitByIntervalBams.cores`|Int|1|The number of cores to allocate to the job.
`mergeSplitByIntervalBams.timeout`|Int|6|Maximum amount of time (in hours) the task can run for.
`mergeSplitByIntervalBams.modules`|String|"gatk/4.1.6.0"|Environment module name and version to load (space separated) before command execution.


### Outputs

Output | Type | Description
---|---|---
`outputGroups`|Array[OutputGroup]|Array of objects with outputIdentifier (from inputGroups) and the final merged bam and bamIndex.
`recalibrationReport`|File?|Recalibration report pdf (if BQSR enabled).
`recalibrationTable`|File?|Recalibration csv that was used by BQSR (if BQSR enabled).


## Commands
 
 This section lists command(s) run by bamMergePreprocessing
 
 * Running bamMergePreprocessing workflow
 
### Parsing Records
 
 ```
     set -euo pipefail
 
     echo "~{str}" | tr '~{lineSeparator}' '\n' | tr '~{recordSeparator}' '\t'
 ```
 
 ### Filtering and marking Duplicates
 
 ```
     set -euxo pipefail
     inputBams="~{sep=" " bams}"
     inputBamIndexes="~{sep=" " bamIndexes}"
 
     # filter
     if [ "~{doFilter}" = true ]; then
       outputBams=()
       outputBamIndexes=()
       for inputBam in $inputBams; do
         filename="$(basename $inputBam ".bam")"
         outputBam="~{workingDir}${filename}.filtered.bam"
         outputBamIndex="~{workingDir}${filename}.filtered.bai"
         samtools view -b \
         -F ~{filterFlags} \
         ~{"-q " + minMapQuality} \
         ~{filterAdditionalParams} \
         $inputBam \
         ~{sep=" " intervals} > $outputBam
         samtools index $outputBam $outputBamIndex
         outputBams+=("$outputBam")
         outputBamIndexes+=("$outputBamIndex")
       done
       # set inputs for next step
       inputBams=("${outputBams[@]}")
       inputBamIndexes=("${outputBamIndexes[@]}")
     else
       outputBams=()
       outputBamIndexes=()
       for inputBam in $inputBams; do
         filename="$(basename $inputBam ".bam")"
         outputBam="~{workingDir}${filename}.bam"
         outputBamIndex="~{workingDir}${filename}.bai"
         samtools view -b \
         $inputBam \
         ~{sep=" " intervals} > $outputBam
         samtools index $outputBam $outputBamIndex
         outputBams+=("$outputBam")
         outputBamIndexes+=("$outputBamIndex")
       done
       # set inputs for next step
       inputBams=("${outputBams[@]}")
       inputBamIndexes=("${outputBamIndexes[@]}")
     fi
 
     # mark duplicates
     if [ "~{doMarkDuplicates}" = true ]; then
       outputBams=()
       outputBamIndexes=()
       gatk --java-options "-Xmx~{jobMemory - overhead}G" MarkDuplicates \
       ${inputBams[@]/#/--INPUT } \
       --OUTPUT="~{markDuplicatesFilePath}.bam" \
       --METRICS_FILE="~{outputFileName}.metrics" \
       --VALIDATION_STRINGENCY=SILENT \
       --REMOVE_DUPLICATES=~{removeDuplicates} \
       --OPTICAL_DUPLICATE_PIXEL_DISTANCE=~{opticalDuplicatePixelDistance} \
       --CREATE_INDEX=true \
       ~{markDuplicatesAdditionalParams}
       outputBams+=("~{markDuplicatesFilePath}.bam")
       outputBamIndexes+=("~{markDuplicatesFilePath}.bai")
       # set inputs for next step
       inputBams=("${outputBams[@]}")
       inputBamIndexes=("${outputBamIndexes[@]}")
     fi
 
     # split N cigar reads
     if [ "~{doSplitNCigarReads}" = true ]; then
       outputBams=()
       outputBamIndexes=()
       gatk --java-options "-Xmx~{jobMemory - overhead}G" SplitNCigarReads \
       ${inputBams[@]/#/--input=} \
       --output="~{splitNCigarReadsFilePath}.bam" \
       --reference ~{reference} \
       ~{sep=" " prefix("--intervals ", intervals)} \
       ~{sep=" " prefixedReadFilters} \
       --create-output-bam-index true \
       --refactor-cigar-string ~{refactorCigarString} \
       ~{splitNCigarReadsAdditionalParams}
       outputBams+=("~{splitNCigarReadsFilePath}.bam")
       outputBamIndexes+=("~{splitNCigarReadsFilePath}.bai")
       # set inputs for next step
       inputBams=("${outputBams[@]}")
       inputBamIndexes=("${outputBamIndexes[@]}")
     fi
 
     # catch all - need to merge filtered+split bams if MarkDuplicates or SplitNCigarReads isn't called
     if [ "~{doMarkDuplicates}" = false ] && [ "~{doSplitNCigarReads}" = false ]; then
       gatk --java-options "-Xmx~{jobMemory - overhead}G" MergeSamFiles \
       ${inputBams[@]/#/--INPUT=} \
       --OUTPUT="~{filteredFileName}.bam" \
       --CREATE_INDEX=true \
       --SORT_ORDER=coordinate \
       --ASSUME_SORTED=false \
       --USE_THREADING=true \
       --VALIDATION_STRINGENCY=SILENT
     fi
 ```
 
 ### Merging bam files
 
 ```
     set -euo pipefail
 
     gatk --java-options "-Xmx~{jobMemory - overhead}G" MergeSamFiles \
     ~{sep=" " prefix("--INPUT=", bams)} \
     --OUTPUT="~{outputFileName}~{suffix}.bam" \
     --CREATE_INDEX=true \
     --SORT_ORDER=coordinate \
     --ASSUME_SORTED=false \
     --USE_THREADING=true \
     --VALIDATION_STRINGENCY=SILENT \
     ~{additionalParams}
 ```
 
 ### realignerTargetCreator processing
 
 ```
     set -euo pipefail
 
     java -Xmx~{jobMemory - overhead}G -jar ~{gatkJar} --analysis_type RealignerTargetCreator \
     --reference_sequence ~{reference} \
     ~{sep=" " prefix("--intervals ", intervals)} \
     ~{sep=" " prefix("--input_file ", bams)} \
     ~{sep=" " prefix("--known ", knownIndels)} \
     --out realignerTargetCreator.intervals \
     ~{"--downsampling_type " + downsamplingType} \
     ~{additionalParams}
 ```
 
 ### 
 
 ```
     set -euo pipefail
 
     # generate gatk nWayOut file
     python3 <<CODE
     import os
     import csv
 
     with open('~{write_lines(bams)}') as f:
         bamFiles = f.read().splitlines()
 
     nWayOut = []
     for bam in bamFiles:
         fileName = os.path.basename(bam)
         realignedFileName = os.path.splitext(fileName)[0] + ".realigned.bam"
         nWayOut.append([fileName, realignedFileName])
 
     with open('input_output.map', 'w') as f:
         tsv_writer = csv.writer(f, delimiter='\t')
         tsv_writer.writerows(nWayOut)
     CODE
 
     java -Xmx~{jobMemory - overhead}G -jar ~{gatkJar} --analysis_type IndelRealigner \
     --reference_sequence ~{reference} \
     ~{sep=" " prefix("--intervals ", intervals)} \
     ~{sep=" " prefix("--input_file ", bams)} \
     --targetIntervals ~{targetIntervals} \
     ~{sep=" " prefix("--knownAlleles ", knownAlleles)} \
     --bam_compression 0 \
     --nWayOut input_output.map \
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
     --output=~{outputFileName} \
     ~{additionalParams}
 ```
 
 ### Gathering Base Quality Score Recalibration Reports
 
 ```
     set -euo pipefail
 
     gatk --java-options "-Xmx~{jobMemory - overhead}G" GatherBQSRReports \
     ~{sep=" " prefix("--input=", recalibrationTables)} \
     --output ~{outputFileName} \
     ~{additionalParams}
 ```
 
 ### Analysis of Covariates
 
 ```
     set -euo pipefail
 
     gatk --java-options "-Xmx~{jobMemory - overhead}G" AnalyzeCovariates \
     --bqsr-recal-file=~{recalibrationTable} \
     --plots-report-file ~{outputFileName} \
     ~{additionalParams}
 ```
 
 ### Appliying Base Quality Score Recalibration
 
 ```
     set -euo pipefail
 
     gatk --java-options "-Xmx~{jobMemory - overhead}G" ApplyBQSR \
     --bqsr-recal-file=~{recalibrationTable} \
     ~{sep=" " prefix("--input=", [bam])} \
     --output ~{outputFileName}~{suffix}.bam \
     ~{additionalParams}
 
 ```
 ### Assemble a list of files by Identifier
 
 ```
     set -euo pipefail
 
     python3 <<CODE
     import json
     import os
     import re
 
     with open('~{write_json(wrappedInputGroups)}') as f:
         inputGroups = json.load(f)
     with open('~{write_lines(bams)}') as f:
         bamFiles = f.read().splitlines()
     with open('~{write_lines(bamIndexes)}') as f:
         bamIndexFiles = f.read().splitlines()
 
     filesByOutputIdentifier = []
     for outputIdentifier in [inputGroup['outputIdentifier'] for inputGroup in inputGroups['inputGroups']]:
         # select bams and bamIndexes for outputIdentifier (preprocessBam prefixes the outputIdentifier, so include that too)
         bams = [bam for bam in bamFiles if re.match("^" + outputIdentifier + "\.", os.path.basename(bam))]
         bais = [bai for bai in bamIndexFiles if re.match("^" + outputIdentifier + "\.", os.path.basename(bai))]
 
         fileNames = list(set([os.path.splitext(os.path.basename(f))[0] for f in bams + bais]))
         if len(fileNames) != 1:
             raise Exception("Unable to determine unique fileName from fileNames = [" + ','.join(f for f in fileNames) + "]")
         else:
             fileName = fileNames[0]
 
         filesByOutputIdentifier.append({
             'outputIdentifier': outputIdentifier,
             'outputFileName': fileName,
             'bams': bams,
             'bamIndexes': bais})
 
     # wrap the array into collectionGroups object
     wrappedFilesByOutputIdentifier = {'collectionGroups': filesByOutputIdentifier}
 
     with open('filesByOutputIdentifier.json', 'w') as f:
         json.dump(wrappedFilesByOutputIdentifier, f, indent=4)
     CODE
 ```
 ## Support

For support, please file an issue on the [Github project](https://github.com/oicr-gsi) or send an email to gsi@oicr.on.ca .

_Generated with generate-markdown-readme (https://github.com/oicr-gsi/gsi-wdl-tools/)_
