version 1.0

struct bamFiles {
  File bam
  File bamIndex
}

struct GenomeResources {
  Array[String] known_indels
  Array[String] known_alleles
  Array[String] known_sites
}

workflow bamMergePreprocessing {

  input {
    Array[bamFiles] inputBamFiles
    String outputFileNamePrefix
    String intervalsToParallelizeByString
    Boolean doFilter = true
    Boolean doMarkDuplicates = true
    Boolean doBqsr = false
    String reference
    String reference_genome
  }

  parameter_meta {
    inputBamFiles: "Array of objects describing sets of bams to merge together and the merged file name. These merged bams will be cocleaned together and output separately (by merged name)."
    intervalsToParallelizeByString: "Comma separated list of intervals to split by (e.g. chr1,chr2,chr3+chr4)."
    doFilter: "Enable/disable Samtools filtering."
    doMarkDuplicates: "Enable/disable GATK4 MarkDuplicates."
    doBqsr: "Enable/disable GATK baseQualityScoreRecalibration"
    reference: "Path to reference file."
    outputFileNamePrefix: "Prefix of output file name"
    reference_genome: "The reference genome version for input sample"
  }

  meta {
    author: "Michael Laszloffy"
    email: "michael.laszloffy@oicr.on.ca"
    description: ""
    dependencies: [
      {
        name: "samtools/1.9",
        url: "http://www.htslib.org/"
      },
      {
        name: "gatk/4.1.6.0",
        url: "https://gatk.broadinstitute.org"
      },
      {
        name: "gatk/3.6-0",
        url: "https://gatk.broadinstitute.org"
      },
      {
       name: "python/3.7",
       url: "https://www.python.org"
      }
    ]
    output_meta: {
      mergedBam: "the final merged bam.",
      mergedBamIndex: "the final merged bam index",
      recalibrationReport: "Recalibration report pdf (if BQSR enabled).",
      recalibrationTable: "Recalibration csv that was used by BQSR (if BQSR enabled)."
    }
  }

  Map[String,GenomeResources] resources = {
    "hg19": {
      "known_indels": ["/.mounts/labs/gsi/modulator/sw/data/hg19-dbsnp-leftaligned-138/dbsnp_138.hg19.leftAligned.vcf.gz"],
      "known_alleles": ["/.mounts/labs/gsi/modulator/sw/data/hg19-dbsnp-leftaligned-138/dbsnp_138.hg19.leftAligned.vcf.gz"],
      "known_sites": ["/.mounts/labs/gsi/modulator/sw/data/hg19-dbsnp-leftaligned-138/dbsnp_138.hg19.leftAligned.vcf.gz"]
    },
    "hg38": {
      "known_indels": ["/.mounts/labs/gsi/modulator/sw/data/hg38-dbsnp-138/dbsnp_138.hg38.vcf.gz"],
      "known_alleles": ["/.mounts/labs/gsi/modulator/sw/data/hg38-dbsnp-138/dbsnp_138.hg38.vcf.gz"],
      "known_sites": ["/.mounts/labs/gsi/modulator/sw/data/hg38-dbsnp-138/dbsnp_138.hg38.vcf.gz"]
    },
    "mm10": {
      "known_indels": ["/.mounts/labs/gsi/modulator/sw/data/mm10-dbsnp-150/dbsnp_150.mm10.vcf.gz"],
      "known_alleles": ["/.mounts/labs/gsi/modulator/sw/data/mm10-dbsnp-150/dbsnp_150.mm10.vcf.gz"],
      "known_sites": ["/.mounts/labs/gsi/modulator/sw/data/mm10-dbsnp-150/dbsnp_150.mm10.vcf.gz"]
    }
  }

  call splitStringToArray {
    input:
      str = intervalsToParallelizeByString
  }
  Array[Array[String]] intervalsToParallelizeBy = splitStringToArray.out

  scatter (intervals in intervalsToParallelizeBy) {
    scatter (i in inputBamFiles) {
      call preprocessBam {
        input:
          inputBam = i.bam,
          inputBamIndex = i.bamIndex,
          outputFileName = outputFileNamePrefix,
          intervals = intervals,
          reference = reference,
          doFilter = doFilter,
          doMarkDuplicates = doMarkDuplicates
      }
    }
    Array[File] preprocessedBams = preprocessBam.preprocessedBam
    Array[File] preprocessedBamIndexes = preprocessBam.preprocessedBamIndex

    if(doBqsr) {
      call baseQualityScoreRecalibration {
        input:
          bams = preprocessedBams,
          reference = reference,
          knownSites = resources[reference_genome].known_sites
      }
    }
    File? recalibrationTableByInterval = baseQualityScoreRecalibration.recalibrationTable
  }
  Array[File] processedBams = flatten(preprocessedBams)
  Array[File] processedBamIndexes = flatten(preprocessedBamIndexes)

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

  Array[File] bamsToMerge = select_first([recalibratedBams, processedBams])
  String outputFileName = basename(bamsToMerge[0], ".bam")
  call mergeBams {
    input:
      bams = bamsToMerge,
      outputFileName = outputFileName
  }

  output {
    File mergedBam = mergeBams.mergedBam
    File mergedBamIndex = mergeBams.mergedBamIndex
    File? recalibrationReport = analyzeCovariates.recalibrationReport
    File? recalibrationTable = gatherBQSRReports.recalibrationTable
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
    str: "Interval string to split (e.g. chr1,chr2,chr3+chr4)."
    lineSeparator: "Interval group separator - these are the intervals to split by."
    recordSeparator: "Interval interval group separator - this can be used to combine multiple intervals into one group."
    jobMemory: "Memory allocated to job (in GB)."
    cores: "The number of cores to allocate to the job."
    timeout: "Maximum amount of time (in hours) the task can run for."
    modules: "Environment module name and version to load (space separated) before command execution."
  }
}

task preprocessBam {
  input {
    Boolean doFilter = true
    Boolean doMarkDuplicates = true
    String outputFileName

    # by default write tmp files to the current working directory (cromwell task directory)
    # $TMPDIR is set by Cromwell
    # $TMP is set by Univa
    String temporaryWorkingDir = ""

    File inputBam
    File inputBamIndex
    Array[String] intervals

    # filter parameters
    String filterSuffix = ".filtered"
    Int filterFlags = 260
    Int? minMapQuality
    String? filterAdditionalParams 

    # markDuplicates parameters
    String dedupSuffix = ".deduped"
    Boolean removeDuplicates = false
    Int opticalDuplicatePixelDistance = 100
    String? markDuplicatesAdditionalParams

    String reference
    Int jobMemory = 48
    Int overhead = 8
    Int cores = 1
    Int timeout = 6
    String modules = "samtools/1.9 gatk/4.1.6.0"
  }

  String workingDir = if temporaryWorkingDir == "" then "" else "~{temporaryWorkingDir}/"

  String baseFileName = "~{outputFileName}"

  String filteredFileName = if doFilter then
                            "~{baseFileName}~{filterSuffix}"
                           else
                            "~{baseFileName}"
  String filteredFilePath = if doMarkDuplicates then
                            "~{workingDir}~{filteredFileName}"
                           else "~{filteredFileName}"

  String markDuplicatesFileName = if doMarkDuplicates then
                                  "~{filteredFileName}~{dedupSuffix}"
                                 else
                                  "~{filteredFileName}"

  command <<<
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

  output {
    File preprocessedBam = if doMarkDuplicates then
                            "~{markDuplicatesFileName}.bam"
                           else if doFilter then
                            "~{filteredFilePath}.bam"
                           else "~{filteredFileName}.bam"
    File preprocessedBamIndex = if doMarkDuplicates then
                                  "~{markDuplicatesFileName}.bai"
                                else if doFilter then
                                  "~{filteredFilePath}.bai"
                                else "~{filteredFileName}.bai"
    File? markDuplicateMetrics = "~{markDuplicatesFileName}.metrics"
  }

  runtime {
    memory: "~{jobMemory} GB"
    cpu: "~{cores}"
    timeout: "~{timeout}"
    modules: "~{modules}"
  }

  parameter_meta {
    doFilter: "Enable/disable Samtools filtering."
    outputFileName: "Output files will be prefixed with this."
    temporaryWorkingDir: "Where to write out intermediary bam files. Only the final preprocessed bam will be written to task working directory if this is set to local tmp."
    inputBam: "bam files to process."
    inputBamIndex: "index files for input bam."
    intervals: "One or more genomic intervals over which to operate."
    filterSuffix: "Suffix to use for filtered bams."
    filterFlags: "Samtools filter flags to apply."
    dedupSuffix: "Suffix to use for markDuplcated bams"
    minMapQuality: "Samtools minimum mapping quality filter to apply."
    filterAdditionalParams: "Additional parameters to pass to samtools."
    reference: "Path to reference file."
    jobMemory:  "Memory allocated to job (in GB)."
    overhead: "Java overhead memory (in GB). jobMemory - overhead == java Xmx/heap memory."
    cores: "The number of cores to allocate to the job."
    timeout: "Maximum amount of time (in hours) the task can run for."
    modules: "Environment module name and version to load (space separated) before command execution."
  }
}

task mergeBams {
  input {
    Array[File] bams
    String outputFileName
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
    --OUTPUT="~{outputFileName}.bam" \
    --CREATE_INDEX=true \
    --SORT_ORDER=coordinate \
    --ASSUME_SORTED=false \
    --USE_THREADING=true \
    --VALIDATION_STRINGENCY=SILENT \
    ~{additionalParams}
  >>>

  output {
    File mergedBam = "~{outputFileName}.bam"
    File mergedBamIndex = "~{outputFileName}.bai"
  }

  runtime {
    memory: "~{jobMemory} GB"
    cpu: "~{cores}"
    timeout: "~{timeout}"
    modules: "~{modules}"
  }

  parameter_meta {
    bams: "Array of bam files to merge together."
    outputFileName: "Output files will be prefixed with this."
    additionalParams: "Additional parameters to pass to GATK MergeSamFiles."
    jobMemory: "Memory allocated to job (in GB)."
    overhead: "Java overhead memory (in GB). jobMemory - overhead == java Xmx/heap memory."
    cores: "The number of cores to allocate to the job."
    timeout: "Maximum amount of time (in hours) the task can run for."
    modules: "Environment module name and version to load (space separated) before command execution."
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
    bams: "Array of bam files to produce a recalibration table for."
    reference: "Path to reference file."
    intervals: "One or more genomic intervals over which to operate."
    knownSites: "Array of VCF with known polymorphic sites used to exclude regions around known polymorphisms from analysis."
    additionalParams: "Additional parameters to pass to GATK BaseRecalibrator."
    outputFileName: "Recalibration table file name."
    jobMemory:  "Memory allocated to job (in GB)."
    overhead: "Java overhead memory (in GB). jobMemory - overhead == java Xmx/heap memory."
    cores: "The number of cores to allocate to the job."
    timeout: "Maximum amount of time (in hours) the task can run for."
    modules: "Environment module name and version to load (space separated) before command execution."
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
    recalibrationTables: "Array of recalibration tables to merge."
    additionalParams: "Additional parameters to pass to GATK GatherBQSRReports."
    outputFileName: "Recalibration table file name."
    jobMemory:  "Memory allocated to job (in GB)."
    overhead: "Java overhead memory (in GB). jobMemory - overhead == java Xmx/heap memory."
    cores: "The number of cores to allocate to the job."
    timeout: "Maximum amount of time (in hours) the task can run for."
    modules: "Environment module name and version to load (space separated) before command execution."
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
    recalibrationTable: "Recalibration table to produce report for."
    additionalParams: "Additional parameters to pass to GATK AnalyzeCovariates"
    outputFileName: "Recalibration report file name."
    jobMemory:  "Memory allocated to job (in GB)."
    overhead: "Java overhead memory (in GB). jobMemory - overhead == java Xmx/heap memory."
    cores: "The number of cores to allocate to the job."
    timeout: "Maximum amount of time (in hours) the task can run for."
    modules: "Environment module name and version to load (space separated) before command execution."
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
    recalibrationTable: "Recalibration table to apply to all input bams."
    bam: "Bam file to recalibrate."
    outputFileName: "Output files will be prefixed with this."
    suffix: "Suffix to use for recalibrated bams."
    additionalParams: "Additional parameters to pass to GATK ApplyBQSR."
    jobMemory:  "Memory allocated to job (in GB)."
    overhead: "Java overhead memory (in GB). jobMemory - overhead == java Xmx/heap memory."
    cores: "The number of cores to allocate to the job."
    timeout: "Maximum amount of time (in hours) the task can run for."
    modules: "Environment module name and version to load (space separated) before command execution."
  }
}

