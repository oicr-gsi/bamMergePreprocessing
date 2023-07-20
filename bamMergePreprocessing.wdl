version 1.0

struct bamFiles {
  File bam
  File bamIndex
}

struct outputGroup {
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
    Boolean doSplitNCigarReads = false
    Boolean doBqsr = false
    String reference
    String reference_genome
  }

  parameter_meta {
    inputBamFiles: "Array of objects describing sets of bams to merge together and the merged file name. These merged bams will be cocleaned together and output separately (by merged name)."
    intervalsToParallelizeByString: "Comma separated list of intervals to split by (e.g. chr1,chr2,chr3+chr4)."
    doFilter: "Enable/disable Samtools filtering."
    doMarkDuplicates: "Enable/disable GATK4 MarkDuplicates."
    doSplitNCigarReads: "Enable/disable GATK4 SplitNCigarReads."
    reference: "Path to reference file."
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
    ],
    output_meta: {
      outputBamFile: "the final merged bam and bamIndex.",
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
    "hgmm10": {
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

  scatter (i in inputBamFiles) {
    scatter (intervals in intervalsToParallelizeBy) {
      call preprocessBam {
        input:
          inputBam = i.bam,
          inputBamIndex = i.bamIndex,
          outputFileName = outputFileNamePrefix,
          intervals = intervals,
          reference = reference,
          doFilter = doFilter,
          doMarkDuplicates = doMarkDuplicates,
          doSplitNCigarReads = doSplitNCigarReads
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
    File? recalibrationTableByLane = baseQualityScoreRecalibration.recalibrationTable
  }
  Array[File] processedBams = flatten(preprocessedBams)
  Array[File] processedBamIndexes = flatten(preprocessedBamIndexes)

  if(doBqsr) {
    call gatherBQSRReports {
      input:
        recalibrationTables = select_all(recalibrationTableByLane)
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

  call mergeBams {
    input:
      bams = select_first([recalibratedBams, processedBams]),
      outputFileName = outputFileNamePrefix,
      suffix = "" # collectFilesBySample task generates the file name
  }

    outputGroup outputBamfile = { 
                              "bam": mergeBams.mergedBam,
                              "bamIndex": mergeBams.mergedBamIndex
                            }

  output {
    outputGroup outputBamFile = outputBamfile
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
    Boolean doSplitNCigarReads = false

    String outputFileName

    # by default write tmp files to the current working directory (cromwell task directory)
    # $TMPDIR is set by Cromwell
    # $TMP is set by Univa
    String temporaryWorkingDir = ""

    File inputBam
    File inputBamIndex
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
                            "~{baseFileName}.filtered"
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

    # filter
    if [ "~{doFilter}" = true ]; then
      outputBam="~{workingDir}~{baseFileName}.filtered.bam"
      outputBamIndex="~{workingDir}~{baseFileName}.filtered.bai"
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
      --INPUT $inputBam \
      --OUTPUT="~{markDuplicatesFilePath}.bam" \
      --METRICS_FILE="~{outputFileName}.metrics" \
      --VALIDATION_STRINGENCY=SILENT \
      --REMOVE_DUPLICATES=~{removeDuplicates} \
      --OPTICAL_DUPLICATE_PIXEL_DISTANCE=~{opticalDuplicatePixelDistance} \
      --CREATE_INDEX=true \
      ~{markDuplicatesAdditionalParams}

      # set inputs for next step
      inputBam=$outputBam
      inputBamIndex=$outputBamIndex
    fi

    # split N cigar reads
    if [ "~{doSplitNCigarReads}" = true ]; then
      gatk --java-options "-Xmx~{jobMemory - overhead}G" SplitNCigarReads \
      --input $inputBam \
      --output="~{splitNCigarReadsFilePath}.bam" \
      --reference ~{reference} \
      ~{sep=" " prefix("--intervals ", intervals)} \
      ~{sep=" " prefixedReadFilters} \
      --create-output-bam-index true \
      --refactor-cigar-string ~{refactorCigarString} \
      ~{splitNCigarReadsAdditionalParams}
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

  parameter_meta {
    doFilter: "Enable/disable Samtools filtering."
    doMarkDuplicates: "Enable/disable GATK4 MarkDuplicates."
    doSplitNCigarReads: "Enable/disable GATK4 SplitNCigarReads."
    outputFileName: "Output files will be prefixed with this."
    temporaryWorkingDir: "Where to write out intermediary bam files. Only the final preprocessed bam will be written to task working directory if this is set to local tmp."
    inputBam: "bam files to process."
    inputBamIndex: "index files for input bam."
    intervals: "One or more genomic intervals over which to operate."
    filterSuffix: "Suffix to use for filtered bams."
    filterFlags: "Samtools filter flags to apply."
    minMapQuality: "Samtools minimum mapping quality filter to apply."
    filterAdditionalParams: "Additional parameters to pass to samtools."
    markDuplicatesSuffix: "Suffix to use for duplicate marked bams."
    removeDuplicates: "MarkDuplicates remove duplicates?"
    opticalDuplicatePixelDistance: "MarkDuplicates optical distance."
    markDuplicatesAdditionalParams: "Additional parameters to pass to GATK MarkDuplicates."
    splitNCigarReadsSuffix: "Suffix to use for SplitNCigarReads bams."
    reference: "Path to reference file."
    refactorCigarString: "SplitNCigarReads refactor cigar string?"
    readFilters: "SplitNCigarReads read filters"
    splitNCigarReadsAdditionalParams: "Additional parameters to pass to GATK SplitNCigarReads."
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

