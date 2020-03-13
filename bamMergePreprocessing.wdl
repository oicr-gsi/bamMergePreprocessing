version 1.0

workflow bamMergePreprocessing {

  input {
    Array[InputGroup] inputGroups
    String intervalsToParallelizeByString
    Boolean doFilterBam = true
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

  scatter (i in inputGroups) {
    if(doFilterBam){
      scatter(f in i.inputFiles) {
        call filterBam {
          input:
            bam = f,
            outputFileName = if length(i.inputFiles) == 1 then i.outputIdentifier else basename(f, ".bam"),
            suffix = ".filter"
        }
      }
    }

    if(length(i.inputFiles) == 1) {
      if(defined(filterBam.filteredBam)) {
        File filteredSingleBam = select_first([filterBam.filteredBam])[0]
        File filteredSingleBamIndex = select_first([filterBam.filteredBamIndex])[0]
      }
      if(!defined(filterBam.filteredBam)) {
        call renameAndIndexBam {
          input:
            bam = i.inputFiles[0],
            outputFileName = i.outputIdentifier,
            suffix = ""
        }
      }
    }
    if(length(i.inputFiles) > 1) {
      call mergeBams {
        input:
          bams = select_first([filterBam.filteredBam, i.inputFiles]),
          outputFileName = i.outputIdentifier,
          suffix = if doFilterBam then ".filter.merge" else ".merge"
      }
    }
    File singleBam = select_first([mergeBams.mergedBam, filteredSingleBam, renameAndIndexBam.renamedBam])
    File singleBamIndex = select_first([mergeBams.mergedBamIndex, filteredSingleBamIndex, renameAndIndexBam.renamedBamIndex])

    if(doMarkDuplicates) {
      call markDuplicates {
        input:
          bam = singleBam
      }
    }

    if(doSplitNCigarReads) {
        call splitNCigarReads {
          input:
            bam = select_first([markDuplicates.dedupedBam, singleBam]),
            reference = reference
        }
    }

    File preprocessedBam = select_first([splitNCigarReads.splitBam, markDuplicates.dedupedBam, singleBam])
    File preprocessedBamIndex = select_first([splitNCigarReads.splitBamIndex, markDuplicates.dedupedBamIndex, singleBam])
  }
  Array[File] preprocessedBams = preprocessedBam
  Array[File] preprocessedBamIndexes = preprocessedBamIndex

  if(doIndelRealignment) {
    scatter (intervals in intervalsToParallelizeBy) {
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
    }
    Array[File] indelRealignedBams = flatten(indelRealign.indelRealignedBams)
    Array[File] indelRealignedBamIndexes = flatten(indelRealign.indelRealignedBamIndexes)
  }

  if(doBqsr) {
    call baseQualityScoreRecalibration {
      input:
        bams = select_first([indelRealignedBams, preprocessedBams]),
        reference = reference
    }

    call analyzeCovariates {
      input:
        recalibrationTables = [baseQualityScoreRecalibration.recalibrationTable]
    }

    scatter(f in select_first([indelRealignedBams, preprocessedBams])) {
      call applyBaseQualityScoreRecalibration {
        input:
          recalibrationTables = [baseQualityScoreRecalibration.recalibrationTable],
          bam = f
      }
    }
    Array[File] recalibratedBams = applyBaseQualityScoreRecalibration.recalibratedBam
    Array[File] recalibratedBamIndexes = applyBaseQualityScoreRecalibration.recalibratedBamIndex
  }

  call collectFilesBySample {
    input:
      inputGroups = inputGroups,
      bams = select_first([recalibratedBams, indelRealignedBams, preprocessedBams]),
      bamIndexes = select_first([recalibratedBamIndexes, indelRealignedBamIndexes, preprocessedBamIndexes])
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

task filterBam {
  input {
    File bam
    String outputFileName = basename(bam, ".bam")
    String suffix = ".filter"
    Int filterFlags = 260
    Int? minMapQuality
    String? additionalParams

    Int jobMemory = 24
    Int cores = 1
    Int timeout = 6
    String modules = "samtools/1.9"
  }

  command <<<
    set -euo pipefail

    samtools view -b \
    -F ~{filterFlags} \
    ~{"-q " + minMapQuality} \
    ~{additionalParams} \
    ~{bam} > ~{outputFileName}~{suffix}.bam

    samtools index ~{outputFileName}~{suffix}.bam ~{outputFileName}~{suffix}.bai
  >>>

  output {
    File filteredBam = outputFileName + suffix + ".bam"
    File filteredBamIndex = outputFileName + suffix + ".bai"
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

task renameAndIndexBam {
  input {
    File bam
    String outputFileName
    String suffix = ""

    Int jobMemory = 24
    Int cores = 1
    Int timeout = 6
    String modules = "samtools/1.9"
  }

  command <<<
    set -euo pipefail

    ln -s ~{bam} ~{outputFileName}~{suffix}.bam

    samtools index ~{outputFileName}~{suffix}.bam ~{outputFileName}~{suffix}.bai
  >>>

  output {
    File renamedBam = outputFileName + suffix + ".bam"
    File renamedBamIndex = outputFileName + suffix + ".bai"
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
    String modules = "gatk/4.1.5.0"
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

task markDuplicates {
  input {
    File bam
    String outputFileName = basename(bam, ".bam")
    String suffix = ".deduped"
    Boolean removeDuplicates = false
    Int opticalDuplicatePixelDistance = 100
    String? additionalParams

    Int jobMemory = 24
    Int overhead = 6
    Int cores = 1
    Int timeout = 6
    String modules = "gatk/4.1.5.0"
  }
  command <<<
    set -euo pipefail

    gatk --java-options "-Xmx~{jobMemory - overhead}G" MarkDuplicates \
    --INPUT="~{bam}" \
    --OUTPUT="~{outputFileName}~{suffix}.bam" \
    --METRICS_FILE="~{outputFileName}.metrics" \
    --VALIDATION_STRINGENCY=SILENT \
    --REMOVE_DUPLICATES=~{removeDuplicates} \
    --OPTICAL_DUPLICATE_PIXEL_DISTANCE=~{opticalDuplicatePixelDistance} \
    --CREATE_INDEX=true \
    ~{additionalParams}
  >>>

  output {
    File dedupedBam = "~{outputFileName}~{suffix}.bam"
    File dedupedBamIndex = "~{outputFileName}~{suffix}.bai"
    File metricsFile = "~{outputFileName}.metrics"
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

task splitNCigarReads {
  input {
    File bam
    String outputFileName = basename(bam, ".bam")
    String suffix = ".split"
    String reference
    Boolean refactorCigarString = false
    Array[String] readFilters = []
    String? additionalParams

    Int jobMemory = 24
    Int overhead = 6
    Int cores = 1
    Int timeout = 6
    String modules = "gatk/4.1.5.0"
  }

  # workaround for this issue https://github.com/broadinstitute/cromwell/issues/5092
  # ~{sep = " " prefix("--read-filter ", readFilters)}
  Array[String] prefixedReadFilters = prefix("--read-filter ", readFilters)

  command <<<
    set -euo pipefail

    gatk --java-options "-Xmx~{jobMemory - overhead}G" SplitNCigarReads \
    --input="~{bam}" \
    --output="~{outputFileName}~{suffix}.bam" \
    --reference ~{reference} \
    ~{sep = " " prefixedReadFilters} \
    --create-output-bam-index true \
    --refactor-cigar-string ~{refactorCigarString} \
    ~{additionalParams}
  >>>

  output {
    File splitBam = "~{outputFileName}~{suffix}.bam"
    File splitBamIndex = "~{outputFileName}~{suffix}.bai"
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
    String? additionalParams

    Int jobMemory = 24
    Int overhead = 6
    Int cores = 1
    Int timeout = 6
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
    --downsampling_type NONE \
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
    String modules = "gatk/4.1.5.0"
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

task analyzeCovariates {
  input {
    Array[File] recalibrationTables
    String? additionalParams
    String outputFileName = "gatk.recalibration.pdf"

    Int jobMemory = 24
    Int overhead = 6
    Int cores = 1
    Int timeout = 6
    String modules = "gatk/4.1.5.0"
  }

  command <<<
    set -euo pipefail

    gatk --java-options "-Xmx~{jobMemory - overhead}G" AnalyzeCovariates \
    ~{sep=" " prefix("--bqsr-recal-file=", recalibrationTables)} \
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
    Array[File] recalibrationTables
    File bam
    String outputFileName = basename(bam, ".bam")
    String suffix = ".recalibrated"
    String? additionalParams

    Int jobMemory = 24
    Int overhead = 6
    Int cores = 1
    Int timeout = 6
    String modules = "gatk/4.1.5.0"
  }

  command <<<
    set -euo pipefail

    gatk --java-options "-Xmx~{jobMemory - overhead}G" ApplyBQSR \
    ~{sep=" " prefix("--bqsr-recal-file=", recalibrationTables)} \
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

    with open('~{write_json(wrappedInputGroups)}') as f:
        inputGroups = json.load(f)
    with open('~{write_lines(bams)}') as f:
        bamFiles = f.read().splitlines()
    with open('~{write_lines(bamIndexes)}') as f:
        bamIndexFiles = f.read().splitlines()

    filesByOutputIdentifier = []
    for outputIdentifier in [inputGroup['outputIdentifier'] for inputGroup in inputGroups['inputGroups']]:
        bams = [bam for bam in bamFiles if outputIdentifier in bam]
        bais = [bai for bai in bamIndexFiles if outputIdentifier in bai]
        fileNames = list(set([os.path.splitext(os.path.basename(f))[0] for f in bams + bais]))
        if len(fileNames) != 1:
            raise Exception("Unable to determine fileName, fileNames = [" + ','.join(f for f in fileNames) + "]")
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

struct InputGroup {
  String outputIdentifier
  Array[File]+ inputFiles
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
