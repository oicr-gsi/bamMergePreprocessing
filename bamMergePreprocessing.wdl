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
    String referenceGenome
  }

  parameter_meta {
    inputBamFiles: "Array of objects describing sets of bams to merge together and the merged file name. These merged bams will be cocleaned together and output separately (by merged name)."
    intervalsToParallelizeByString: "Comma separated list of intervals to split by (e.g. chr1,chr2,chr3+chr4)."
    doFilter: "Enable/disable Samtools filtering."
    doMarkDuplicates: "Enable/disable GATK4 MarkDuplicates."
    doBqsr: "Enable/disable GATK baseQualityScoreRecalibration"
    reference: "Path to reference file."
    outputFileNamePrefix: "Prefix of output file name"
    referenceGenome: "The reference genome version for input sample"
  }

  meta {
    author: "Michael Laszloffy and Gavin Peng"
    email: "michael.laszloffy@oicr.on.ca and gpeng@oicr.on.ca"
    description: "Workflow to merge and preprocess lane level alignments."
    dependencies: [
      {
        name: "samtools/1.9",
        url: "http://www.htslib.org/"
      },
      {
        name: "gatk/4.1.6.0",
        url: "https://gatk.broadinstitute.org"
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

  scatter (interval in flatten(intervalsToParallelizeBy)) {
    
    if (length(inputBamFiles) == 1) {
      call getChrCoefficient as coeffForPreprocess {
        input: 
          chromosome = interval,
          bamFile = inputBamFiles[0].bam
      }

      call preprocessBam {
        input:
          inputBam = inputBamFiles[0].bam,
          inputBamIndex = inputBamFiles[0].bamIndex,
          outputFileName = outputFileNamePrefix,
          interval = interval,
          reference = reference,
          scaleCoefficient = coeffForPreprocess.coeff,
          doFilter = doFilter,
          doMarkDuplicates = doMarkDuplicates
      }
      File preprocessedBam = preprocessBam.preprocessedBam
      File preprocessedBamIndex = preprocessBam.preprocessedBamIndex
    }

    if (length(inputBamFiles) > 1) {
      
      scatter (i in inputBamFiles) {
        call getChrCoefficient as coeffForFilter {
          input:
            chromosome = interval,
            bamFile = i.bam
        }

        call filterBam {
          input:
            inputBam = i.bam,
            inputBamIndex = i.bamIndex,
            outputFileName = outputFileNamePrefix,
            interval = interval,
            scaleCoefficient = coeffForFilter.coeff,
            reference = reference,
            doFilter = doFilter
        }
      }
      Array[File] filteredBams = filterBam.filteredBam
      Array[File] filteredBamIndexes = filterBam.filteredBamIndex

      if (doMarkDuplicates) {
        call markDuplicates {
          input:
          inputBams = filteredBams,
          outputFileName = outputFileNamePrefix+".filtered"
        }
      }

      if (!doMarkDuplicates) {
        call mergeBams as mergeMultipleBam {
        input:
          bams = filteredBams,
          outputFileName = outputFileNamePrefix
        }
      }
      File filterDedupedBam = select_first([markDuplicates.dedupedBam, mergeMultipleBam.mergedBam])
      File filterDedupedBamIndex = select_first([markDuplicates.dedupedBamIndex, mergeMultipleBam.mergedBamIndex])
    }

    File processedBam = select_first([preprocessedBam, filterDedupedBam])
    File processedBamIndex = select_first([preprocessedBamIndex, filterDedupedBamIndex])
  }

  Array[File] processedBams = processedBam
  Array[File] processedBamIndexes = processedBamIndex

  if(doBqsr) {
      call baseQualityScoreRecalibration {
        input:
          bams = processedBams,
          reference = reference,
          knownSites = resources[referenceGenome].known_sites
      }

    File recalibrationTableByInterval = baseQualityScoreRecalibration.recalibrationTable

    call gatherBQSRReports {
      input:
        recalibrationTables = recalibrationTableByInterval
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

  call mergeBams {
    input:
      bams = bamsToMerge,
      outputFileName = outputFileNamePrefix
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

# ================================================================
#  Scaling coefficient - use to scale RAM allocation by chromosome
# ================================================================
task getChrCoefficient {
  input {
    Int memory = 2
    Int timeout = 1
    String chromosome
    String modules = "samtools/1.14"
    File bamFile
  }

  parameter_meta {
    bamFile: ".bam file to process, we just need the header"
    timeout: "Hours before task timeout"
    chromosome: "Chromosome to check"
    memory: "Memory allocated for this job"
    modules: "Names and versions of modules to load"
  }

  command <<<
    CHROM_LEN=$(samtools view -H ~{bamFile} | grep ^@SQ | cut -f 2,3 | grep -v _ | grep -w ~{chromosome} | cut -f 2 | sed 's/LN://')
    LARGEST=$(samtools view -H ~{bamFile} | grep ^@SQ | cut -f 2,3 | grep -v _ | cut -f 2 | sed 's/LN://' | sort -n | tail -n 1)
    echo | awk -v chrom_len=$CHROM_LEN -v largest=$LARGEST '{print int((chrom_len/largest + 0.1) * 10)/10}'
  >>>

  runtime {
    memory:  "~{memory} GB"
    modules: "~{modules}"
    timeout: "~{timeout}"
  }

  output {
    String coeff = read_string(stdout())
  }

  meta {
    output_meta: {
      coeff: "Length ratio as relative to the largest chromosome."
    }
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
    String interval

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
    Int jobMemory = 36
    Int overhead = 8
    Int cores = 1
    Int timeout = 6
    Float scaleCoefficient = 1.0
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
      ~{interval} > $outputBam
      samtools index $outputBam $outputBamIndex

      # set inputs for next step
      inputBam=$outputBam
      inputBamIndex=$outputBamIndex
    else
      outputBam="~{workingDir}~{baseFileName}.bam"
      outputBamIndex="~{workingDir}~{baseFileName}.bai"
      samtools view -b \
      ~{inputBam} \
      ~{interval} > $outputBam
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
    memory: "~{round(jobMemory * scaleCoefficient) + 8} GB"
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
    interval: "Genomic interval over which to operate."
    filterSuffix: "Suffix to use for filtered bams."
    dedupSuffix: "Suffix to use for markDuplcated bams"
    removeDuplicates: "MarkDuplicates remove duplicates?"
    opticalDuplicatePixelDistance: "MarkDuplicates optical distance."
    markDuplicatesAdditionalParams: "Additional parameters to pass to GATK MarkDuplicates."
    filterFlags: "Samtools filter flags to apply."
    minMapQuality: "Samtools minimum mapping quality filter to apply."
    filterAdditionalParams: "Additional parameters to pass to samtools."
    reference: "Path to reference file."
    jobMemory:  "Memory allocated to job (in GB)."
    overhead: "Java overhead memory (in GB). jobMemory - overhead == java Xmx/heap memory."
    cores: "The number of cores to allocate to the job."
    scaleCoefficient: "Chromosome-dependent RAM scaling coefficient"
    timeout: "Maximum amount of time (in hours) the task can run for."
    modules: "Environment module name and version to load (space separated) before command execution."
  }
}

task filterBam {
  input {
    Boolean doFilter = true
    String outputFileName

    # by default write tmp files to the current working directory (cromwell task directory)
    # $TMPDIR is set by Cromwell
    # $TMP is set by Univa
    String temporaryWorkingDir = ""

    File inputBam
    File inputBamIndex
    String interval

    # filter parameters
    String filterSuffix = ".filtered"
    Int filterFlags = 260
    Int? minMapQuality
    String? filterAdditionalParams

    String reference
    Int jobMemory = 48
    Int overhead = 8
    Int cores = 1
    Int timeout = 6
    Float scaleCoefficient = 1.0
    String modules = "samtools/1.9 gatk/4.1.6.0"
  }

  String workingDir = if temporaryWorkingDir == "" then "" else "~{temporaryWorkingDir}/"

  String baseFileName = "~{outputFileName}"

  String filteredFileName = if doFilter then
                            "~{baseFileName}~{filterSuffix}"
                           else
                            "~{baseFileName}"

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
      ~{interval} > $outputBam
      samtools index $outputBam $outputBamIndex

      # set inputs for next step
      inputBam=$outputBam
      inputBamIndex=$outputBamIndex
    else
      outputBam="~{workingDir}~{baseFileName}.bam"
      outputBamIndex="~{workingDir}~{baseFileName}.bai"
      samtools view -b \
      ~{inputBam} \
      ~{interval} > $outputBam
      samtools index $outputBam $outputBamIndex
    fi
  >>>

  output {
    File filteredBam = "~{filteredFileName}.bam"
    File filteredBamIndex = "~{filteredFileName}.bai"
  }

  runtime {
    memory: "~{round(jobMemory * scaleCoefficient)} GB"
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
    interval: "Genomic interval over which to operate."
    filterSuffix: "Suffix to use for filtered bams."
    filterFlags: "Samtools filter flags to apply."
    minMapQuality: "Samtools minimum mapping quality filter to apply."
    filterAdditionalParams: "Additional parameters to pass to samtools."
    reference: "Path to reference file."
    jobMemory:  "Memory allocated to job (in GB)."
    overhead: "Java overhead memory (in GB). jobMemory - overhead == java Xmx/heap memory."
    cores: "The number of cores to allocate to the job."
    scaleCoefficient: "Chromosome-dependent RAM scaling coefficient"
    timeout: "Maximum amount of time (in hours) the task can run for."
    modules: "Environment module name and version to load (space separated) before command execution."
  }
}

task markDuplicates {
  input {
    Array[File]inputBams
    String outputFileName
    Boolean removeDuplicates = false
    Int opticalDuplicatePixelDistance = 100
    String? markDuplicatesAdditionalParams
    Int jobMemory = 24
    Int overhead = 6
    Int cores = 1
    Int timeout = 6
    String modules = "gatk/4.1.6.0"
    String dedupSuffix = ".deduped"
  }

    command <<<
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

  output {
    File dedupedBam = outputFileName + dedupSuffix + ".bam"
    File dedupedBamIndex = outputFileName + dedupSuffix + ".bai"
  }

  runtime {
    memory: "~{jobMemory} GB"
    cpu: "~{cores}"
    timeout: "~{timeout}"
    modules: "~{modules}"
  }

  parameter_meta {
    inputBams: "Array of bam files to go through markDuplicates."
    dedupSuffix: "Suffix to use for markDuplcated bams"
    removeDuplicates: "MarkDuplicates remove duplicates?"
    opticalDuplicatePixelDistance: "MarkDuplicates optical distance."
    markDuplicatesAdditionalParams: "Additional parameters to pass to GATK MarkDuplicates."
    jobMemory: "Memory allocated to job (in GB)."
    overhead: "Java overhead memory (in GB). jobMemory - overhead == java Xmx/heap memory."
    cores: "The number of cores to allocate to the job."
    timeout: "Maximum amount of time (in hours) the task can run for."
    modules: "Environment module name and version to load (space separated) before command execution."
  }
}

task mergeBams {
  input {
    Array[File] bams
    String baseName = basename(bams[0])
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

  output {
    File mergedBam = glob("*.bam")[0]
    File mergedBamIndex = glob("*.bai")[0]
  }

  runtime {
    memory: "~{jobMemory} GB"
    cpu: "~{cores}"
    timeout: "~{timeout}"
    modules: "~{modules}"
  }

  parameter_meta {
    bams: "Array of bam files to merge together."
    baseName: "The base name for output files"
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
    File recalibrationTables
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
    --input ~{recalibrationTables} \
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
    recalibrationTables: "Recalibration tables to merge."
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
