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

workflow bamMergePreprocessingMOD {

    input {
        Array[bamFiles] inputBamFiles
        String outputFileNamePrefix
        String intervalsToParallelizeByString = "chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY,chrM,NC,SPLIT,UNALIGNED"
        Boolean doFilter = true
        Boolean doMarkDuplicates = true
        Boolean doBqsr = true
        Boolean provisionBqsr = false
        String libType = "dna"
        String reference
        String referenceGenome
        Boolean doBamMetrics = false
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
        provisionBqsr: "Enable/disable provision out bqsr report and table"
    }

    output {
        File mergedBam = finalBam
        File mergedBamIndex = finalBamIndex
        File? recalibrationReport = analyzeCovariates.recalibrationReport
        File? recalibrationTable = gatherBQSRReports.recalibrationTable
        File? markDuplicateMetricsZip = zippedDupMarkMetrics.zipFile
        File? bamMetricsZip = zippedBamMetrics.zipFile
    }

 

    meta {
        author: "Michael Laszloffy, Gavin Peng & Lawrence Heisler"
        email: "mlaszloffy@oicr.on.ca  gpeng@oicr.on.ca. lheisler@oicr.on.ca"
        description: "Workflow to merge and preprocess lane level alignments."
        dependencies: [
            {
                name: "samtools/1.15",
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
            recalibrationTable: "Recalibration csv that was used by BQSR (if BQSR enabled).",
            markDuplicateMetrics: "A tarball of markDuplicates metrics file of all chromsomes"
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

    # given a set of intervals for parallelization, split into an array that can be scattered across
    call prepareIntervals {
        input:
            str = intervalsToParallelizeByString,
            refFai = "~{reference + '.fai'}"
        }
    

    String suffixFilter = if doFilter then ".filtered" else ""
    String suffixDupMarked = if doMarkDuplicates then ".dupmarked" else ".dupunmarked"
    String suffixRecalibrated = if doBqsr then ".recalibrated" else ""	
    String suffixSplitNCigarString = if libType == "rna" then ".split" else ""	
    #Array[Array[String]] intervalsToParallelizeBy = prepareIntervals.intervals
    String WorkflowVersion = "v2.1"
    String header_comment = "CallReady BAM file generated from the bamMergePreprocessing ~{WorkflowVersion} Workflow. Filtering=~{doFilter},DuplicateMarking=~{doMarkDuplicates},BQSR=~{doBqsr},libType=~{libType}"

    if(doBamMetrics){
       scatter (inputBamFile in inputBamFiles){
           call bamMetrics as inputBamMetrics {
               input:
               inputBam = inputBamFile.bam,
           }
       }
    }

    ###################################################
    ### scatter TWICE, first across the set of intervalsm then the set of input bam files
    ### subsetAndFilter will process each input bam splitting to the intervals before filtering
    ### baseQualityRecalibration will process each filtered per interval bam file to generate a recalibration table
    ###
    ##################################################
    Array[String] intervals =  flatten(prepareIntervals.intervals)
    scatter (interval in intervals){
        scatter (i in inputBamFiles) {
            String bamId = basename(i.bam,".bam")
            String filterPrefix = "~{bamId + suffixFilter + "." + interval}"

            ### this is repetitive and should give the same results, do once before scattering with the reference information
            #call getChrCoefficient as coeffForPreprocess {
            #    input:
            #    chromosome = interval,
            #    bamFile = i.bam
            #}
  
            ###################################################
            ### subsetAndFilter will subset the input bam by interval
            ### intervals are a chromosome + coordinates : these extract reads where both map on the same chromosome
            ### there are several keywords that indicate special treatment
            ### NC : this will subset a set of non-canonical chromosomes together using -L bedfile
            ### SPLIT : this will subset pairs on different chromosomes
            ### UNALIGNED : this will subset when both reads are UNALIGNED.  Unaligned reads paired with a mapped read are extracted with the mapped read
            ### filtering by defaults is -F256, removing non-primary mappings.  this can be modified
            ###
            ##################################################
            call subsetAndFilter {
                input:
                    inputBam = i.bam,
                    inputBamIndex = i.bamIndex,
                    outputFileNamePrefix =  filterPrefix,
                    interval = interval,
                    ncBed = prepareIntervals.ncbed,
                    reference = reference,
                    scaleCoefficient = prepareIntervals.intervalCoefficients[interval],
                    doFilter = doFilter
            }

            #### bqsr can be done on these per-interval,per-input bams files to generate the table, and later collected across all
            #### is there any issue to generate these tables BEFORE duplicate marking?
            if(doBqsr) {
                ### bqsr takes a list of bam files, capture the single bam file into an array
                call baseQualityScoreRecalibration {
                    input:
                    bams = [subsetAndFilter.bam],
                    reference = reference,
                    knownSites = resources[referenceGenome].known_sites
                }
           }
          
       }     #### END by inputBam scatter. ####
       Array[File] subsetAndFilterBams_byInterval = subsetAndFilter.bam
       

        ###################################################
        ### MERGING AND DUPLICATE MARKING, within Intervals
        ### Merging is EITHER done with MarkDuplicates, or with MergeSamFiles
        ### both will collect the subsetAndFilter bam files within an interval
        ### markDuplicates will gather the bams from the same interval
        ### each interval will produce a merged bam
        ###
        ##################################################
        if(doMarkDuplicates) {
            String markdupPrefix = "~{outputFileNamePrefix + suffixFilter + suffixDupMarked + "." + interval}"
            ### within the interval markDuplicates across all bams.  this will also merge within the interval
            call markDuplicates {
                input:
                inputBams = subsetAndFilter.bam,
                outputFileNamePrefix = markdupPrefix,
                scaleCoefficient = prepareIntervals.intervalCoefficients[interval]
           }
        }
        if(!doMarkDuplicates) {
            String mergePrefix = "~{outputFileNamePrefix + suffixFilter + "." + interval}"
            call mergeBams as mergeWithinInterval {
                input:
                bams = subsetAndFilter.bam,
                outputFileNamePrefix = mergePrefix
            }
        }


        ###################################################
        ### rna libraries should be processed with splitNCigarString
        ### this will split reads into one or more new reads
        ### this is done AFTER duplicate marking or merging the data
        ### each interval will produce a merged bam
        ###
        ##################################################
        if (libType == "rna") {
            ### this processed the data through SplitNCigarString
            ## splitNCigarString will take a single bam file, it should not be used for merging

            String splitPrefix = "~{outputFileNamePrefix + suffixFilter + suffixDupMarked + suffixSplitNCigarString + "." + interval}"

            call splitNCigarString {
               input:
               inputBam = select_first([markDuplicates.bam, mergeWithinInterval.bam]),
               inputBamIndex = select_first([markDuplicates.bamIndex, mergeWithinInterval.bamIndex]),
               outputFileNamePrefix =  splitPrefix,
               reference = reference
           }
        }

        #### collect the bams for each intervals into these variables
        #### this is the bam that will go to the next step
        File mergedIntervalBam = select_first([splitNCigarString.bam, markDuplicates.bam, mergeWithinInterval.bam])
        File mergedIntervalBamIndex = select_first([splitNCigarString.bamIndex, markDuplicates.bamIndex, mergeWithinInterval.bamIndex])

        ### organize by interval to later process through ApplyBQSR
        Pair[String,File] IntervalToBam = (interval,mergedIntervalBam)        
        
        #### this array of bams are for bamStats
        Array[File] mergedPerIntervalBams = select_all([splitNCigarString.bamIndex, markDuplicates.bamIndex, mergeWithinInterval.bamIndex])

        
         

    } ### END by Interval scatter  ###


    ###################################################
    ### run bamMetrics on any bams generated in the scatter
    ### tthis will include, subsetAndFilter, markDuplicates, mergeWithinInterval and splitNCigarString
    ##################################################
    if(doBamMetrics){
        ### collect all bam output
        Array[File] subsetAndFilterBams = flatten(subsetAndFilterBams_byInterval)
        scatter(subsetAndFilterBam in subsetAndFilterBams){
            call bamMetrics as subsetAndFilterBamMetrics{
                input:
                    inputBam = subsetAndFilterBam
            }
        }
        if(doMarkDuplicates){
           #Array[File] duplicateMarkedBam = [select_first(markDuplicates.bam)]
           Array[File] duplicateMarkedBams = select_all(markDuplicates.bam)
           scatter(duplicateMarkedBam in duplicateMarkedBams){
               call bamMetrics as duplicateMarkedBamMetrics{
                   input:
                       inputBam = duplicateMarkedBam
               }
           }
        }
        if(!doMarkDuplicates){
            #Array[File] mergedWithinIntervalBams = [select_first(mergeWithinInterval.bam)]
            Array[File] mergedWithinIntervalBams = select_all(mergeWithinInterval.bam)
            scatter(mergeWithinIntervalBam in mergedWithinIntervalBams){
                call bamMetrics as mergedWithinIntervalBamMetrics{
                    input:
                        inputBam = mergeWithinIntervalBam
                }
            }
        }

        if(libType == "rna"){
            #Array[File] splitNCigarStringBams = [select_first(splitNCigarString.bam)]
            Array[File] splitNCigarStringBams = select_all(splitNCigarString.bam)
            scatter(splitNCigarStringBam in splitNCigarStringBams){
                call bamMetrics as splitNCigarStringBamMetrics{
                    input:
                        inputBam = splitNCigarStringBam
                }
            }
        }
    }


    ###################################################
    ### all per interval, per input processing done, producing a merged bam file
    ### recalibration still needs to be applied (if requested), and summary of the recalibration
    ##################################################
    if(doBqsr) {
        ####  each byInterval, byInputbam should have a recalibration tabl
        ##### these can now be gathered together
        ####. recalibration tables are optional, and are in an array of array of file.  flatten + select all

        #### check if this is correct, these are being generated in a double scatter, and i may not be gathering correctly
        Array[File] recalibrationTables = select_all(flatten(baseQualityScoreRecalibration.recalibrationTable))

        call gatherBQSRReports {
            input:
            recalibrationTables = recalibrationTables,
            outputFileNamePrefix = "~{outputFileNamePrefix + suffixFilter + suffixDupMarked + suffixSplitNCigarString + suffixRecalibrated}"
        }
        call analyzeCovariates {
            input:
            recalibrationTable = gatherBQSRReports.recalibrationTable,
            outputFileNamePrefix = "~{outputFileNamePrefix + suffixFilter + suffixDupMarked + suffixSplitNCigarString + suffixRecalibrated}"
        }
        ### applying to the final merge bam, scattering across the perInterval bams
        scatter(itb in IntervalToBam){
            String intervalName = itb.left
            File intervalBam = itb.right
            
            String applyBqsrPrefix = "~{outputFileNamePrefix + suffixFilter + suffixDupMarked + suffixSplitNCigarString + suffixRecalibrated + "." + intervalName}"
            call applyBaseQualityScoreRecalibration {
                input:
                recalibrationTable = gatherBQSRReports.recalibrationTable,
                bam = intervalBam,
                outputFileNamePrefix = applyBqsrPrefix
            }
        }

        if(doBamMetrics){
            Array[File] recalibratedBams = flatten([applyBaseQualityScoreRecalibration.recalibratedBam])
            scatter(recalibratedBam in recalibratedBams){
                call bamMetrics as recalibratedBamMetrics{
                    input:
                       inputBam = recalibratedBam
                }
            }
        }
    }


    ###################################################
    ### merge across intervals
    ### recalibration still needs to be applied (if requested), and summary of the recalibration
    ### this is selecting the output from either applyBQSR or the gathered mergeWitiinIntervalBam 
    ### (either markDups.ba or mergedWithinInterval.bam)
    ##################################################
    String finalPrefix = "~{outputFileNamePrefix + suffixFilter + suffixDupMarked + suffixSplitNCigarString + suffixRecalibrated + ".merged"}"
    call mergeBams as mergeAcrossIntervals {
        input:
            bams = select_first([applyBaseQualityScoreRecalibration.recalibratedBam,mergedIntervalBam]),
            outputFileNamePrefix = finalPrefix,
            ### i would like to add a url but GATK seems to fail when including : character, calling it a tagged argument.  Its not clear how to escape
            comment = "CallReady BAM file generated from the bamMergePreprocessingWorkflow. Filtering=~{doFilter},DuplicateMarking=~{doMarkDuplicates},BQSR=~{doBqsr},libType=~{libType}"
    }
    
    ### metrics on the final merge, and zip everything together
    if(doBamMetrics){
        ##### generated the metrics for the final output
        Array[File] finalBams = [mergeAcrossIntervals.bam]
        scatter (finalBam in finalBams){
            call bamMetrics as finalBamMetrics{
                input:
                    inputBam = finalBam
            }
        }

        ##### collect ALL metrics files here, and send to make zip
        Array[File] inputBamMetricFiles = flatten(select_all([inputBamMetrics.stats, inputBamMetrics.flagstats,inputBamMetrics.counts]))
        Array[File] subsetAndFilterBamMetricFiles = flatten(select_all([subsetAndFilterBamMetrics.stats, subsetAndFilterBamMetrics.flagstats,subsetAndFilterBamMetrics.counts]))
        Array[File] duplicateMarkedBamMetricsFiles = flatten(select_all([duplicateMarkedBamMetrics.stats, duplicateMarkedBamMetrics.flagstats,duplicateMarkedBamMetrics.counts]))
        Array[File] mergedWithinIntervalBamMetricsFiles = flatten(select_all([mergedWithinIntervalBamMetrics.stats, mergedWithinIntervalBamMetrics.flagstats,mergedWithinIntervalBamMetrics.counts]))
        Array[File] splitNCigarStringBamMetricsFiles = flatten(select_all([splitNCigarStringBamMetrics.stats, splitNCigarStringBamMetrics.flagstats,splitNCigarStringBamMetrics.counts]))
        Array[File] recalibratedBamMetricsFiles = flatten(select_all([recalibratedBamMetrics.stats, recalibratedBamMetrics.flagstats,recalibratedBamMetrics.counts]))
        Array[File] finalBamMetricsFiles = flatten([finalBamMetrics.stats, finalBamMetrics.flagstats,finalBamMetrics.counts])
        
        Array[File] allBamMetricFiles = flatten([inputBamMetricFiles,subsetAndFilterBamMetricFiles,duplicateMarkedBamMetricsFiles,mergedWithinIntervalBamMetricsFiles,splitNCigarStringBamMetricsFiles,recalibratedBamMetricsFiles,finalBamMetricsFiles]) 
        call makeZip as zippedBamMetrics{
            input:
               inputFiles = allBamMetricFiles,
               zipName = "~{outputFileNamePrefix + '.bamMetrics'}"
        }
    }



    File finalBam = mergeAcrossIntervals.bam
    File finalBamIndex = mergeAcrossIntervals.bamIndex

    ### collect the dedup metrics into a zip file, only if doMarkDuplicates
    if(doMarkDuplicates) {
        ### to manage optional files
        Array[File] metricFiles = select_all(markDuplicates.metrics)
        call makeZip as zippedDupMarkMetrics{
            input:
                inputFiles = metricFiles,
                zipName = "~{outputFileNamePrefix + '.markDuplicatesMetrics'}"
        }
    }
  
  
  ### for debugging
  #File finalBam = inputBamFiles[0].bam
  #File finalBamIndex = inputBamFiles[0].bamIndex
  
}

# =========================================================================
#  run samstats on a bam file
#  this is primarily used for the input bams
#  it possibly be repurposed for other bam files generated in the workflow
#     which are currently paired with the samtools stats command
# =========================================================================
task bamMetrics {
  input{
    File inputBam
    Int jobMemory = 12
    Int cores = 1
    Int timeout = 6
    String modules = "samtools/1.15"
  }
  String prefix = basename(inputBam,".bam")
  command <<<
    set -euo pipefail
    ### pipe the SAM file through flagstat (tee), stats (tee) and a counting operation
    samtools view -h ~{inputBam} |  tee >(samtools flagstat - > ~{prefix}.flagstats.txt) | tee >(samtools stats - > ~{prefix}.samstats.txt) | samtools view - | cut -f 2,3,7 | sort | uniq -c > ~{prefix}.mapcounts.txt

  >>>

  output {
    File stats = "~{prefix}.samstats.txt"
    File flagstats = "~{prefix}.flagstats.txt"
    File counts = "~{prefix}.mapcounts.txt"
  }

  runtime {
    memory: "~{jobMemory} GB"
    cpu: "~{cores}"
    timeout: "~{timeout}"
    modules: "~{modules}"
  }

  
}



# =============================================================================
#  split the string of intervals to an array
#  review intervals for keywords
#  prepare the NC list from the index.fai
#. generated coefficients from the interval size to use for memory management
# =============================================================================
task prepareIntervals {
  input {
    String str
    String lineSeparator = ","
    String recordSeparator = "+"
    String refFai
    Int jobMemory = 1
    Int cores = 1
    Int timeout = 1
    String modules = ""
  }

  command <<<
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
  >>>

  output {
    ### WDL 1.0 does not support keys function, otherwise the intervals could be extracted from intervalCoefficients
    Array[Array[String]] intervals = read_tsv("intervals")
    Map[String,String] intervalCoefficients = read_map("coefficients.txt")
    File ncbed = "nc.contigs.bed" 
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
    refFai: "The reference fai file, with a list of the contigs available in the build"
    jobMemory: "Memory allocated to job (in GB)."
    cores: "The number of cores to allocate to the job."
    timeout: "Maximum amount of time (in hours) the task can run for."
    modules: "Environment module name and version to load (space separated) before command execution."
  }
}




# ================================================================
#  given an interval, subset the bam files using samtools
#  if filtering is required, apply the filters here
# ================================================================
task subsetAndFilter  {
  input {
    Boolean doFilter = true
    String outputFileNamePrefix

    # by default write tmp files to the current working directory (cromwell task directory)
    # $TMPDIR is set by Cromwell
    # $TMP is set by Univa
    String temporaryWorkingDir = ""

    File inputBam
    File inputBamIndex
    String interval
    File ncBed

    # filter parameters
    String filterSuffix = ".filtered"
    Int filterFlags = 256
    Int? minMapQuality
    String? filterAdditionalParams

    String reference
	
	## this is a hidden option that we revert back to the previous dupMarking strategy.  For assessment.  can likely be removed
    Boolean oldStyle = false

    Int jobMemory = 36
    Int minMemory = 12
    Int overhead = 8
    Int cores = 1
    Int timeout = 6
    Float scaleCoefficient = 1.0
    String modules = "samtools/1.15"

  }

  ### determine the expr value
  ### define the samtool parameters for various interval

  
  String exprString0 = if interval == "SPLIT" then "--expr 'rname != rnext'" else "--expr 'rname == rnext'"
  String exprString = if oldStyle then "" else exprString0
  
  ### to handle optional values
  String flags = if defined(filterFlags) then '~{"-F " + filterFlags}' else ""
  String quality = if defined(minMapQuality) then '~{"-q" + minMapQuality}' else ""
  String addParams = if defined(filterAdditionalParams) then '~{"" + filterAdditionalParams}' else ""
  String unaligned = if interval == "UNALIGNED" then "-f 12" else ""
  String filterString = if doFilter then "~{flags} ~{quality} ~{addParams} ~{unaligned}" else ""

  ### use the interval unless one of these keywords   SPLIT UNALIGNED NC
  String samtoolsInterval = if (interval == "SPLIT" || interval == "UNALIGNED" || interval == "NC") then "" else interval
  #String intervalsString = if interval == "NC" then "-L nc.bed" else ""
  String intervalsString = if interval == "NC" then "-L ~{ncBed}" else ""

  Int allocatedMemory = if minMemory > round(jobMemory * scaleCoefficient) then minMemory else round(jobMemory * scaleCoefficient)

  command <<<
    set -euxo pipefail
    
    ### write to local file
    #### dev fix to get rid of canonical chromosome

    #cat ~{ncBed} > nc.bed
    #sleep 10

    samtools view -b ~{exprString} ~{filterString} ~{intervalsString} ~{inputBam} ~{samtoolsInterval} > ~{outputFileNamePrefix}.bam
    samtools index ~{outputFileNamePrefix}.bam ~{outputFileNamePrefix}.bai



  >>>

  output {
    File bam = "~{outputFileNamePrefix}.bam"
    File bamIndex = "~{outputFileNamePrefix}.bai"
  }
  runtime {
    memory: "~{allocatedMemory} GB"
    cpu: "~{cores}"
    timeout: "~{timeout}"
    modules: "~{modules}"
  }

  parameter_meta {
    doFilter: "Enable/disable Samtools filtering."
    outputFileNamePrefix: "Output files will be prefixed with this."
    temporaryWorkingDir: "Where to write out intermediary bam files. Only the final preprocessed bam will be written to task working directory if this is set to local tmp."
    inputBam: "bam files to process."
    inputBamIndex: "index files for input bam."
    interval: "Genomic interval over which to operate."
    filterFlags: "Samtools filter flags to apply."
    minMapQuality: "Samtools minimum mapping quality filter to apply."
    filterAdditionalParams: "Additional parameters to pass to samtools."
    reference: "Path to reference file."
    jobMemory:  "Memory allocated to job (in GB)."
    minMemory: "A minimum amount of memory allocated to the task, overrides the scaled RAM setting"
    overhead: "Java overhead memory (in GB). jobMemory - overhead == java Xmx/heap memory."
    cores: "The number of cores to allocate to the job."
    scaleCoefficient: "Chromosome-dependent RAM scaling coefficient"
    timeout: "Maximum amount of time (in hours) the task can run for."
    modules: "Environment module name and version to load (space separated) before command execution."
  }
}



# ================================================================
#  given an set of bam files, markDuplicates across all and merge
#  this is generally done within an interval, across all input bam files
# ================================================================

task markDuplicates {
  input {
    Array[File]inputBams
    String outputFileNamePrefix
    Boolean removeDuplicates = false
    Int opticalDuplicatePixelDistance = 100
    String? markDuplicatesAdditionalParams
    Int jobMemory = 36
    Int minMemory = 12
    Int overhead = 6
    Int cores = 1
    Int timeout = 6
    String modules = "gatk/4.1.6.0 samtools/1.15"
    Float scaleCoefficient = 1.0
  }

    Int allocatedMemory = if minMemory > round(jobMemory * scaleCoefficient) then minMemory else round(jobMemory * scaleCoefficient)

    command <<<
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

    
    >>>

  output {
    File bam = "~{outputFileNamePrefix}.bam"
    File bamIndex = "~{outputFileNamePrefix}.bai"
    File metrics = "~{outputFileNamePrefix}.metrics"
  }

  runtime {
    memory: "~{jobMemory} GB"
    cpu: "~{cores}"
    timeout: "~{timeout}"
    modules: "~{modules}"
  }

  parameter_meta {
    inputBams: "Array of bam files to go through markDuplicates."
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


# ================================================================
#  given an set of bam files, reds with N in the cigar string are
#  split into multiple reads
#  
# 
# ================================================================

task splitNCigarString {
  input {
    #Array[File]inputBams
    File inputBam
    File inputBamIndex
    String outputFileNamePrefix
    String refactorCigarString = false
    String? splitNCigarReadsAdditionalParams
    String? readFilters
    String reference
    Int jobMemory = 24
    Int overhead = 6
    Int cores = 1
    Int timeout = 6
    String modules = "gatk/4.1.6.0 samtools/1.15"
  }

    command <<<
    set -euo pipefail

    gatk --java-options "-Xmx~{jobMemory - overhead}G" SplitNCigarReads \
    --INPUT= ~{inputBam}  \
    --OUTPUT ~{outputFileNamePrefix}.bam \
    --REFERENCE ~{reference} \
    --CREATE_INDEX=true \
    --REFACTOR-CIGAR-STRING ~{refactorCigarString} \
    ~{splitNCigarReadsAdditionalParams}

    
    >>>

  output {
    File bam = "~{outputFileNamePrefix}.bam"
    File bamIndex = "~{outputFileNamePrefix}.bai"
  }

  runtime {
    memory: "~{jobMemory} GB"
    cpu: "~{cores}"
    timeout: "~{timeout}"
    modules: "~{modules}"
  }

  parameter_meta {
    jobMemory: "Memory allocated to job (in GB)."
    overhead: "Java overhead memory (in GB). jobMemory - overhead == java Xmx/heap memory."
    cores: "The number of cores to allocate to the job."
    timeout: "Maximum amount of time (in hours) the task can run for."
    modules: "Environment module name and version to load (space separated) before command execution."
  }
}



task makeZip {
  input {
    Array[File] inputFiles
    String zipName
    Int jobMemory = 4
    Int cores = 1
    Int timeout = 1
  }

  command <<<
    set -euo pipefail
    mkdir ./files/
    files="~{sep="," inputFiles}"    
    IFS=',' read -ra f <<< "$files"
    for f in ${f[@]}
    do
      cp $f ./files/
    done
    tar czf  ~{zipName}.tar.gz ./files/*
  >>>

  output {
    File zipFile = "./~{zipName}.tar.gz"
  }

  runtime {
    memory: "~{jobMemory} GB"
    cpu: "~{cores}"
    timeout: "~{timeout}"
  }

  parameter_meta {
    inputFiles: "Array of input files to zip together."
    jobMemory: "Memory allocated to job (in GB)."
    cores: "The number of cores to allocate to the job."
    timeout: "Maximum amount of time (in hours) the task can run for."
  }
}



# ================================================================
#  given an set of bam files, merge and Index
#  this is within an interval, and then across Intervals
# ================================================================

task mergeBams {
  input {
    Array[File] bams
    String outputFileNamePrefix
    String? additionalParams
    String comment = ""
    Int jobMemory = 24
    Int overhead = 6
    Int cores = 1
    Int timeout = 6
    String modules = "gatk/4.1.6.0 samtools/1.15"

  }

  command <<<
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


  >>>

  output {
    File bam = "~{outputFileNamePrefix}.bam"
    File bamIndex = "~{outputFileNamePrefix}.bai"

  }
  runtime {
    memory: "~{jobMemory} GB"
    cpu: "~{cores}"
    timeout: "~{timeout}"
    modules: "~{modules}"
  }

  parameter_meta {
    bams: "Array of bam files to merge together."
    outputFileNamePrefix: "Output files will be prefixed with this."
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
    --output=gatk.recalibration.csv \
    ~{additionalParams}
  >>>

  output {
    File recalibrationTable = "gatk.recalibration.csv"
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
    String outputFileNamePrefix

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
    --output ~{outputFileNamePrefix}.gatk.recalibration.csv \
    ~{additionalParams}
  >>>

  output {
    File recalibrationTable = "~{outputFileNamePrefix}.gatk.recalibration.csv"
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
    outputFileNamePrefix: "prefix for Recalibration table file name."
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
    String outputFileNamePrefix

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
    --plots-report-file ~{outputFileNamePrefix}.gatk.recalibration.pdf \
    ~{additionalParams}
  >>>

  output {
    File recalibrationReport = "~{outputFileNamePrefix}.gatk.recalibration.pdf"
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
    outputFileNamePrefix: "Recalibration report file name."
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
    String outputFileNamePrefix
    String? additionalParams

    Int jobMemory = 24
    Int overhead = 6
    Int cores = 1
    Int timeout = 6
    String modules = "gatk/4.1.6.0 samtools/1.15"

  }

  command <<<
    set -euo pipefail

    gatk --java-options "-Xmx~{jobMemory - overhead}G" ApplyBQSR \
    --bqsr-recal-file=~{recalibrationTable} \
    ~{sep=" " prefix("--input=", [bam])} \
    --output ~{outputFileNamePrefix}.bam \
    ~{additionalParams}

  >>>

  output {
    File recalibratedBam = "~{outputFileNamePrefix + '.bam'}"
    File recalibratedBamIndex = "~{outputFileNamePrefix + '.bai'}"
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
    outputFileNamePrefix: "Output files will be prefixed with this."
    additionalParams: "Additional parameters to pass to GATK ApplyBQSR."
    jobMemory:  "Memory allocated to job (in GB)."
    overhead: "Java overhead memory (in GB). jobMemory - overhead == java Xmx/heap memory."
    cores: "The number of cores to allocate to the job."
    timeout: "Maximum amount of time (in hours) the task can run for."
    modules: "Environment module name and version to load (space separated) before command execution."
  }
}
