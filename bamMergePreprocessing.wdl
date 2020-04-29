version 1.0

workflow bamMergePreprocessing {

  input {
    Array[InputGroup] inputGroups
    String intervalsToParallelizeByString
    Boolean doFilter = true
    Boolean doMarkDuplicates = true
    Boolean doSplitNCigarReads = false
    Boolean doIndelRealignment = true
    Boolean doBqsr = true
    String reference
  }

  call splitStringToArray {
    input:
      str = intervalsToParallelizeByString
  }
  Array[Array[String]] intervalsToParallelizeBy = splitStringToArray.out

  scatter (intervals in intervalsToParallelizeBy) {
    scatter (i in inputGroups) {
      scatter(bamAndBamIndexInput in i.bamAndBamIndexInputs) {
        File inputGroupBam = bamAndBamIndexInput.bam
        File inputGroupBamIndex = bamAndBamIndexInput.bamIndex
      }
      Array[File] inputGroupBams = inputGroupBam
      Array[File] inputGroupBamIndexes = inputGroupBamIndex
      call preprocessBam {
        input:
          bams = inputGroupBams,
          bamIndexes = inputGroupBamIndexes,
          intervals = intervals,
          outputFileName = i.outputIdentifier,
          reference = reference,
          doFilter = doFilter,
          doMarkDuplicates = doMarkDuplicates,
          doSplitNCigarReads = doSplitNCigarReads
      }
    }
    Array[File] preprocessedBams = preprocessBam.preprocessedBam
    Array[File] preprocessedBamIndexes = preprocessBam.preprocessedBamIndex

    # indel realignment combines samples (nWayOut) and is parallized by chromosome
    if(doIndelRealignment) {
      call realignerTargetCreator {
        input:
          bams = preprocessedBams,
          bamIndexes = preprocessedBamIndexes,
          intervals = intervals,
          reference = reference
      }

      call indelRealign {
        input:
          bams = preprocessedBams,
          bamIndexes = preprocessedBamIndexes,
          intervals = intervals,
          targetIntervals = realignerTargetCreator.targetIntervals,
          reference = reference
      }
      Array[File] indelRealignedBams = indelRealign.indelRealignedBams
      Array[File] indelRealignedBamIndexes = indelRealign.indelRealignedBamIndexes
    }

    if(doBqsr) {
      call baseQualityScoreRecalibration {
        input:
          bams = select_first([indelRealignedBams, preprocessedBams]),
          reference = reference
      }
    }
    Array[File] processedBamsByInterval = select_first([indelRealignedBams, preprocessedBams])
    Array[File] processedBamIndexesByInterval = select_first([indelRealignedBamIndexes, preprocessedBamIndexes])
    File? recalibrationTableByInterval = baseQualityScoreRecalibration.recalibrationTable
  }
  Array[File] processedBams = flatten(processedBamsByInterval)
  Array[File] processedBamIndexes = flatten(processedBamIndexesByInterval)

  if(doBqsr) {
    call gatherBQSRReports {
      input:
        recalibrationTables = select_all(recalibrationTableByInterval)
    }

    call analyzeCovariates {
      input:
        recalibrationTable = gatherBQSRReports.recalibrationTable
    }

    scatter(bam in processedBams) {
      call applyBaseQualityScoreRecalibration {
        input:
          recalibrationTable = gatherBQSRReports.recalibrationTable,
          bam = bam
      }
    }
    Array[File] recalibratedBams = applyBaseQualityScoreRecalibration.recalibratedBam
    Array[File] recalibratedBamIndexes = applyBaseQualityScoreRecalibration.recalibratedBamIndex
  }

  call collectFilesBySample {
    input:
      inputGroups = inputGroups,
      bams = select_first([recalibratedBams, processedBams]),
      bamIndexes = select_first([recalibratedBamIndexes, processedBamIndexes])
  }

  scatter(o in collectFilesBySample.filesByOutputIdentifier.collectionGroups) {
    if(length(o.bams) > 1) {
      call mergeBams as mergeSplitByIntervalBams {
        input:
          bams = o.bams,
          outputFileName = o.outputFileName,
          suffix = "" # collectFilesBySample task generates the file name
      }
    }
    OutputGroup outputGroup = { "outputIdentifier": o.outputIdentifier,
                                "bam": select_first([mergeSplitByIntervalBams.mergedBam, o.bams[0]]),
                                "bamIndex": select_first([mergeSplitByIntervalBams.mergedBamIndex, o.bamIndexes[0]])}
  }

  output {
    Array[OutputGroup] outputGroups = outputGroup
    File? recalibrationReport = analyzeCovariates.recalibrationReport
  }
}

task splitStringToArray {
  input {
    String str
    String lineSeparator = ","
    String recordSeparator = "+"

    Int jobMemory = 1
    Int cores = 1
    Int timeout = 1
    String modules = ""
  }

  command <<<
    set -euo pipefail

    echo "~{str}" | tr '~{lineSeparator}' '\n' | tr '~{recordSeparator}' '\t'
  >>>

  output {
    Array[Array[String]] out = read_tsv(stdout())
  }

  runtime {
    memory: "~{jobMemory} GB"
    cpu: "~{cores}"
    timeout: "~{timeout}"
    modules: "~{modules}"
  }

  parameter_meta {

  }

  meta {
    output_meta: {

    }
  }
}

task preprocessBam {
  input {
    Boolean doFilter = true
    Boolean doMarkDuplicates = true
    Boolean doSplitNCigarReads = false

    String outputFileName

    # by default write tmp files to the current working directory (cromwell task directory)
    # $TMPDIR is set by Cromwell
    # $TMP is set by Univa
    String temporaryWorkingDir = ""

    Array[File] bams
    Array[File] bamIndexes
    Array[String] intervals

    # filter parameters
    String filterSuffix = ".filter"
    Int filterFlags = 260
    Int? minMapQuality
    String? filterAdditionalParams

    # mark duplicates
    String markDuplicatesSuffix = ".deduped"
    Boolean removeDuplicates = false
    Int opticalDuplicatePixelDistance = 100
    String? markDuplicatesAdditionalParams

    # split N cigar reads
    String splitNCigarReadsSuffix = ".split"
    String reference
    Boolean refactorCigarString = false
    Array[String] readFilters = []
    String? splitNCigarReadsAdditionalParams

    Int jobMemory = 24
    Int overhead = 6
    Int cores = 1
    Int timeout = 6
    String modules = "samtools/1.9 gatk/4.1.6.0"
  }

  String workingDir = if temporaryWorkingDir == "" then "" else "~{temporaryWorkingDir}/"

  String baseFileName = "~{outputFileName}"

  String filteredFileName = if doFilter then
                            "~{baseFileName}.filter"
                           else
                            "~{baseFileName}"
  String filteredFilePath = if doMarkDuplicates || doSplitNCigarReads then
                            "~{workingDir}~{filteredFileName}"
                           else "~{filteredFileName}"

  String markDuplicatesFileName = if doMarkDuplicates then
                                  "~{filteredFileName}.deduped"
                                 else
                                  "~{filteredFileName}"
  String markDuplicatesFilePath = if doSplitNCigarReads then
                                  "~{workingDir}~{markDuplicatesFileName}"
                                 else
                                  "~{markDuplicatesFileName}"

  String splitNCigarReadsFileName = if doSplitNCigarReads then
                                    "~{markDuplicatesFileName}.split"
                                   else
                                    "~{markDuplicatesFileName}"
  String splitNCigarReadsFilePath = if false then # there are no downstream steps, so don't write to temp dir
                                    "~{workingDir}~{splitNCigarReadsFileName}"
                                   else
                                    "~{splitNCigarReadsFileName}"

  # workaround for this issue https://github.com/broadinstitute/cromwell/issues/5092
  # ~{sep = " " prefix("--read-filter ", readFilters)}
  Array[String] prefixedReadFilters = prefix("--read-filter ", readFilters)

  command <<<
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
  >>>

  output {
    File preprocessedBam = if doSplitNCigarReads then
                            "~{splitNCigarReadsFilePath}.bam"
                           else if doMarkDuplicates then
                            "~{markDuplicatesFilePath}.bam"
                           else if doFilter then
                            "~{filteredFilePath}.bam"
                           else "~{filteredFileName}.bam"
    File preprocessedBamIndex = if doSplitNCigarReads then
                                  "~{splitNCigarReadsFilePath}.bai"
                                else if doMarkDuplicates then
                                  "~{markDuplicatesFilePath}.bai"
                                else if doFilter then
                                  "~{filteredFilePath}.bai"
                                else "~{filteredFileName}.bai"
    File? markDuplicateMetrics = "~{outputFileName}.metrics"
  }

  runtime {
    memory: "~{jobMemory} GB"
    cpu: "~{cores}"
    timeout: "~{timeout}"
    modules: "~{modules}"
  }
}

task mergeBams {
  input {
    Array[File] bams
    String outputFileName
    String suffix = ".merge"
    String? additionalParams

    Int jobMemory = 24
    Int overhead = 6
    Int cores = 1
    Int timeout = 6
    String modules = "gatk/4.1.6.0"
  }

  command <<<
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
  >>>

  output {
    File mergedBam = "~{outputFileName}~{suffix}.bam"
    File mergedBamIndex = "~{outputFileName}~{suffix}.bai"
  }

  runtime {
    memory: "~{jobMemory} GB"
    cpu: "~{cores}"
    timeout: "~{timeout}"
    modules: "~{modules}"
  }

  parameter_meta {

  }

  meta {
    output_meta: {

    }
  }
}

task realignerTargetCreator {
  input {
    Array[File] bams
    Array[File] bamIndexes
    String reference
    Array[String] knownIndels
    Array[String] intervals
    String? downsamplingType
    String? additionalParams

    Int jobMemory = 24
    Int overhead = 6
    Int cores = 1
    Int timeout = 6

    # use gatk3 for now: https://github.com/broadinstitute/gatk/issues/3104
    String modules = "gatk/3.6-0"
    String gatkJar = "$GATK_ROOT/GenomeAnalysisTK.jar"
  }

  command <<<
    set -euo pipefail

    #--phone_home NO_ET --gatk_key /.mounts/labs/PDE/data/gatkAnnotationResources/GATK_public.key --logging_level INFO
    java -Xmx~{jobMemory - overhead}G -jar ~{gatkJar} --analysis_type RealignerTargetCreator \
    --reference_sequence ~{reference} \
    ~{sep=" " prefix("--intervals ", intervals)} \
    ~{sep=" " prefix("--input_file ", bams)} \
    ~{sep=" " prefix("--known ", knownIndels)} \
    --out realignerTargetCreator.intervals \
    ~{"--downsampling_type " + downsamplingType} \
    ~{additionalParams}
  >>>

  output {
    File targetIntervals = "realignerTargetCreator.intervals"
  }

  runtime {
    memory: "~{jobMemory} GB"
    cpu: "~{cores}"
    timeout: "~{timeout}"
    modules: "~{modules}"
  }

  parameter_meta {

  }

  meta {
    output_meta: {

    }
  }
}

task indelRealign {
  input {
    Array[File] bams
    Array[File] bamIndexes
    Array[String] intervals
    String reference
    Array[String] knownAlleles
    File targetIntervals
    String? additionalParams

    Int jobMemory = 24
    Int overhead = 6
    Int cores = 1
    Int timeout = 6

    # use gatk3 for now: https://github.com/broadinstitute/gatk/issues/3104
    String modules = "gatk/3.6-0"
    String gatkJar = "$GATK_ROOT/GenomeAnalysisTK.jar"
  }

  command <<<
    set -euo pipefail

    # generate gatk nWayOut file
    python <<CODE
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

    #--phone_home NO_ET  --gatk_key /.mounts/labs/PDE/data/gatkAnnotationResources/GATK_public.key --logging_level INFO
    java -Xmx~{jobMemory - overhead}G -jar ~{gatkJar} --analysis_type IndelRealigner \
    --reference_sequence ~{reference} \
    ~{sep=" " prefix("--intervals ", intervals)} \
    ~{sep=" " prefix("--input_file ", bams)} \
    --targetIntervals ~{targetIntervals} \
    ~{sep=" " prefix("--knownAlleles ", knownAlleles)} \
    --bam_compression 0 \
    --nWayOut input_output.map \
    ~{additionalParams}
  >>>

  output {
    Array[File] indelRealignedBams = glob("*.bam")
    Array[File] indelRealignedBamIndexes = glob("*.bai")
  }

  runtime {
    memory: "~{jobMemory} GB"
    cpu: "~{cores}"
    timeout: "~{timeout}"
    modules: "~{modules}"
  }

  parameter_meta {

  }

  meta {
    output_meta: {

    }
  }
}

task baseQualityScoreRecalibration {
  input {
    Array[File] bams
    String reference
    Array[String] intervals = []
    Array[String] knownSites
    String? additionalParams
    String outputFileName = "gatk.recalibration.csv"

    Int jobMemory = 24
    Int overhead = 6
    Int cores = 1
    Int timeout = 6
    String modules = "gatk/4.1.6.0"
  }

  # workaround for this issue https://github.com/broadinstitute/cromwell/issues/5092
  # ~{sep=" " prefix("--intervals ", intervals)}
  Array[String] prefixedIntervals = prefix("--intervals ", intervals)

  command <<<
    set -euo pipefail

    gatk --java-options "-Xmx~{jobMemory - overhead}G" BaseRecalibrator \
    --reference ~{reference} \
    ~{sep=" " prefixedIntervals} \
    ~{sep=" " prefix("--input=", bams)} \
    ~{sep=" " prefix("--known-sites ", knownSites)} \
    --output=~{outputFileName} \
    ~{additionalParams}
  >>>

  output {
    File recalibrationTable = outputFileName
  }

  runtime {
    memory: "~{jobMemory} GB"
    cpu: "~{cores}"
    timeout: "~{timeout}"
    modules: "~{modules}"
  }

  parameter_meta {

  }

  meta {
    output_meta: {

    }
  }
}

task gatherBQSRReports {
  input {
    Array[File] recalibrationTables
    String? additionalParams
    String outputFileName = "gatk.recalibration.csv"

    Int jobMemory = 24
    Int overhead = 6
    Int cores = 1
    Int timeout = 6
    String modules = "gatk/4.1.6.0"
  }

  command <<<
    set -euo pipefail

    gatk --java-options "-Xmx~{jobMemory - overhead}G" GatherBQSRReports \
    ~{sep=" " prefix("--input=", recalibrationTables)} \
    --output ~{outputFileName} \
    ~{additionalParams}
  >>>

  output {
    File recalibrationTable = outputFileName
  }

  runtime {
    memory: "~{jobMemory} GB"
    cpu: "~{cores}"
    timeout: "~{timeout}"
    modules: "~{modules}"
  }

  parameter_meta {

  }

  meta {
    output_meta: {

    }
  }
}

task analyzeCovariates {
  input {
    File recalibrationTable
    String? additionalParams
    String outputFileName = "gatk.recalibration.pdf"

    Int jobMemory = 24
    Int overhead = 6
    Int cores = 1
    Int timeout = 6
    String modules = "gatk/4.1.6.0"
  }

  command <<<
    set -euo pipefail

    gatk --java-options "-Xmx~{jobMemory - overhead}G" AnalyzeCovariates \
    --bqsr-recal-file=~{recalibrationTable} \
    --plots-report-file ~{outputFileName} \
    ~{additionalParams}
  >>>

  output {
    File recalibrationReport = outputFileName
  }

  runtime {
    memory: "~{jobMemory} GB"
    cpu: "~{cores}"
    timeout: "~{timeout}"
    modules: "~{modules}"
  }

  parameter_meta {

  }

  meta {
    output_meta: {

    }
  }
}

task applyBaseQualityScoreRecalibration {
  input {
    File recalibrationTable
    File bam
    String outputFileName = basename(bam, ".bam")
    String suffix = ".recalibrated"
    String? additionalParams

    Int jobMemory = 24
    Int overhead = 6
    Int cores = 1
    Int timeout = 6
    String modules = "gatk/4.1.6.0"
  }

  command <<<
    set -euo pipefail

    gatk --java-options "-Xmx~{jobMemory - overhead}G" ApplyBQSR \
    --bqsr-recal-file=~{recalibrationTable} \
    ~{sep=" " prefix("--input=", [bam])} \
    --output ~{outputFileName}~{suffix}.bam \
    ~{additionalParams}
  >>>

  output {
    File recalibratedBam = outputFileName + suffix + ".bam"
    File recalibratedBamIndex = outputFileName + suffix + ".bai"
  }

  runtime {
    memory: "~{jobMemory} GB"
    cpu: "~{cores}"
    timeout: "~{timeout}"
    modules: "~{modules}"
  }

  parameter_meta {

  }

  meta {
    output_meta: {

    }
  }
}

task collectFilesBySample {
  input {
    Array[InputGroup] inputGroups
    Array[File] bams
    Array[File] bamIndexes

    Int jobMemory = 1
    Int cores = 1
    Int timeout = 1
    String modules = "python/3.6"
  }

  InputGroups wrappedInputGroups = {"inputGroups": inputGroups}

  command <<<
    set -euo pipefail

    python <<CODE
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
  >>>

  output {
    CollectionGroups filesByOutputIdentifier = read_json("filesByOutputIdentifier.json")
  }

  runtime {
    memory: "~{jobMemory} GB"
    cpu: "~{cores}"
    timeout: "~{timeout}"
    modules: "~{modules}"
  }

  parameter_meta {

  }

  meta {
    output_meta: {

    }
  }
}

struct BamAndBamIndex {
  File bam
  File bamIndex
}

struct InputGroup {
  String outputIdentifier
  Array[BamAndBamIndex]+ bamAndBamIndexInputs
}

struct InputGroups {
  Array[InputGroup] inputGroups
}

struct CollectionGroup {
  String outputIdentifier
  String outputFileName
  Array[File] bams
  Array[File] bamIndexes
}

struct CollectionGroups {
  Array[CollectionGroup] collectionGroups
}

struct OutputGroup {
  String outputIdentifier
  File bam
  File bamIndex
}
