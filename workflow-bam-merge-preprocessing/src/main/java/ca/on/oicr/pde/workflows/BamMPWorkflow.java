package ca.on.oicr.pde.workflows;

import ca.on.oicr.pde.tools.gatk3.AnalyzeCovariates;
import ca.on.oicr.pde.tools.gatk3.BaseRecalibrator;
import ca.on.oicr.pde.tools.gatk3.IndelRealigner;
import ca.on.oicr.pde.tools.gatk3.PrintReads;
import ca.on.oicr.pde.tools.gatk3.RealignerTargetCreator;
import ca.on.oicr.pde.utilities.workflows.SemanticWorkflow;
import com.google.common.collect.HashMultimap;
import com.google.common.collect.Iterables;
import com.google.common.collect.Multimap;
import com.google.common.collect.Sets;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Set;
import net.sourceforge.seqware.common.util.Log;
import net.sourceforge.seqware.pipeline.workflowV2.model.Command;
import net.sourceforge.seqware.pipeline.workflowV2.model.Job;
import net.sourceforge.seqware.pipeline.workflowV2.model.SqwFile;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.lang3.StringUtils;
import org.apache.commons.lang3.tuple.Pair;

/**
 *
 * Merge ➜ Sort ➜ Mark Duplicates (optional q-score filtering) ➜ Realign
 * Indels ➜ Base Q score recalibration
 * <p>
 */
public class BamMPWorkflow extends SemanticWorkflow {

    private String dataDir;
    private String tmpDir;
    private String rDir;
    private String java;
    private String picardJar;
    private Integer picardMarkDupMem;
    private Integer picardMergeMem;
    private Integer picardMemOverhead;
    private boolean doDedup = true;
    private boolean doRemoveDups = true;
    private boolean doFilter = true;
    private boolean doSplitNTrim = false; // flag for split and trim ; set to false unless changed in the workflow.ini
    private final boolean doIndelRealignment = true; //currently this is not an optional step
    private final boolean assumeSorted = false; // we should have assumeSorted set to false ALWAYS (due to re-structuring)
    private boolean picardMergeUseThreading = true;
    private boolean manualOutput = false;
    private String samtools;
    private Integer samtoolsMem;
    private Integer samtoolsFlag;
    private String samtoolsMinMapQuality;
    private String alignerName;
    private String filterOtherParams;
    private String picardDedupOtherParams;
    private String picardMergeOtherParams;
    private String picardSortOtherParams;
    private String picardSortOrder;
    private static final String REMOVE_DUPLICATES = "REMOVE_DUPLICATES=";
    private final static String BAM_METATYPE = "application/bam";
    private final static String BAI_METATYPE = "application/bam-index";
    private final static String PDF_METATYPE = "application/pdf";
    private final static String TXT_METATYPE = "text/plain";

    private List<String> chrSizesList;
    private Set<Set<String>> chrSizes;

    //References
    private String refFasta;
    private String dbsnpVcf;

    // SPLIT CIGAR
    private Integer splitCigarXmxg;
    private Integer splitCigarRMQF;
    private Integer splitCigarRMQT;
    private Boolean flagReassignOneMappingQuality;

    //GATK
    private Integer gatkRealignTargetCreatorXmx;
    private Integer gatkIndelRealignerXmx;
    private Integer gatkPrintReadsXmx;
    private Integer gatkBaseRecalibratorXmx;
    private Integer gatkBaseRecalibratorMem;
    private Integer gatkBaseRecalibratorNct;
    private Integer gatkBaseRecalibratorSmp;
    private Integer gatkOverhead;
    private Integer downsamplingCoverage;
    private Integer preserveQscoresLessThan;
    private Integer intervalPadding;

    private String queue;
    private String realignerTargetCreatorParams;
    private String indelRealignerParams;
    private String baseRecalibratorParams;
    private String printReadsParams;
    private String analyzeCovariatesParams;
    private String downsamplingType;
    private String binDir;
    private String gatk;
    private String gatkKey;
    private Boolean doBQSR;
    private List<String> intervalFilesList;
    private Set<String> intervalFiles;
    private Set<String> bqsrCovariates;

    //Ontology-related variables
    private static final String EDAM = "EDAM";
    private static final Map<String, Set<String>> cvTerms;

    private List<String> outputIdentifiers;

    private Map<String, Set<String>> inputFilesByGroup;
    private Map<String, Set<String>> outputSwidsByGroup;
    private final Map<String, Set<String>> provisionedFilesByGroup = new HashMap<>();

    static {
        cvTerms = new HashMap<>();
        cvTerms.put(EDAM, new HashSet<>(Arrays.asList("BAM", "BAI", "plain text format (unformatted)",
                "Read pre-processing", "Sequence alignment refinement",
                "Sequence alignment metadata")));
    }

    /**
     * Here we need to register all of our CV terms for attaching to result
     * files (Formats, Data, Processes)
     *
     * @return myTerms
     */
    @Override
    protected Map<String, Set<String>> getTerms() {
        Map<String, Set<String>> myTerms = new HashMap<>();
        myTerms.putAll(cvTerms);
        return myTerms;
    }

    /**
     * Launched from setupDirectory();
     */
    private void BamMPWorkflow() {
        java = getProperty("java");
        dataDir = "data";
        tmpDir = getProperty("tmp_dir");
        rDir = getProperty("r_dir");
        manualOutput = Boolean.parseBoolean(getProperty("manual_output"));

        //Flags for controlling different steps:
        doFilter = Boolean.parseBoolean(getProperty("do_sam_filter"));
        doRemoveDups = Boolean.parseBoolean(getOptionalProperty("do_remove_duplicates", "false"));
        doDedup = doRemoveDups ? true : Boolean.parseBoolean(getProperty("do_mark_duplicates"));
        doSplitNTrim = Boolean.parseBoolean(getProperty("do_split_trim_reassign_quality"));

        //Picard
        picardJar = getProperty("picard_jar");
        picardMemOverhead = Integer.parseInt(getProperty("picard_mem_overhead_mb"));

        picardMarkDupMem = Integer.parseInt(getProperty("picard_mark_duplicates_mem_mb"));
        picardDedupOtherParams = getOptionalProperty("picard_mark_duplicates_other_params", "");

        picardMergeMem = Integer.parseInt(getProperty("picard_merge_sam_mem_mb"));
        picardMergeUseThreading = Boolean.parseBoolean(getProperty("picard_merge_use_threading"));
        picardMergeOtherParams = getOptionalProperty("picard_merge_other_params", "");

        picardSortOrder = getProperty("sort_order");
        picardSortOtherParams = getOptionalProperty("picard_sort_other_params", "");

        filterOtherParams = getOptionalProperty("samtools_filter_other_params", "");

        analyzeCovariatesParams = getOptionalProperty("gatk_analyze_covariates_params", null);
        binDir = getWorkflowBaseDir() + "/bin/";
        doBQSR = Boolean.valueOf(getOptionalProperty("do_bqsr", "true"));
        //Samtools
        samtools = getProperty("samtools");
        samtoolsMem = Integer.parseInt(getProperty("samtools_mem_mb"));
        samtoolsFlag = Integer.parseInt(getProperty("samtools_filter_flag"));
        samtoolsMinMapQuality = getOptionalProperty("samtools_min_map_quality", null);

        alignerName = getOptionalProperty("aligner_name", "");

        //References
        refFasta = getProperty("ref_fasta");
        dbsnpVcf = getProperty("gatk_dbsnp_vcf");

        //split trim reassign mapping quality
        splitCigarXmxg = Integer.parseInt(getProperty("split_cigar_Xmxg"));
        splitCigarRMQF = Integer.parseInt(getProperty("split_cigar_RMQF"));
        splitCigarRMQT = Integer.parseInt(getProperty("split_cigar_RMQT"));
        flagReassignOneMappingQuality = Boolean.parseBoolean(getProperty("reassign_One_Mapping_Quality"));

        //GATK
        gatkRealignTargetCreatorXmx = Integer.parseInt(getProperty("gatk_realign_target_creator_xmx"));
        gatkIndelRealignerXmx = Integer.parseInt(getProperty("gatk_indel_realigner_xmx"));
        gatkPrintReadsXmx = Integer.parseInt(getProperty("gatk_print_reads_xmx"));
        gatkBaseRecalibratorXmx = Integer.parseInt(getProperty("gatk_base_recalibrator_xmx"));
        gatkBaseRecalibratorMem = Integer.parseInt(getProperty("gatk_base_recalibrator_mem"));
        gatkBaseRecalibratorNct = Integer.parseInt(getProperty("gatk_base_recalibrator_nct"));
        gatkBaseRecalibratorSmp = Integer.parseInt(getProperty("gatk_base_recalibrator_smp"));
        gatkOverhead = Integer.parseInt(getProperty("gatk_sched_overhead_mem"));
        downsamplingCoverage = hasPropertyAndNotNull("downsampling_coverage") ? Integer.parseInt(getProperty("downsampling_coverage")) : null;
        preserveQscoresLessThan = hasPropertyAndNotNull("preserve_qscores_less_than") ? Integer.parseInt(getProperty("preserve_qscores_less_than")) : null;
        intervalPadding = hasPropertyAndNotNull("interval_padding") ? Integer.parseInt(getProperty("interval_padding")) : null;

        queue = getOptionalProperty("queue", "");
        realignerTargetCreatorParams = getOptionalProperty("gatk_realigner_target_creator_params", null);
        indelRealignerParams = getOptionalProperty("gatk_indel_realigner_params", null);
        baseRecalibratorParams = getOptionalProperty("gatk_base_recalibrator_params", null);
        printReadsParams = getOptionalProperty("gatk_print_reads_params", null);
        downsamplingType = getOptionalProperty("downsampling_type", null);
        gatk = getOptionalProperty("gatk_jar", binDir);
        gatkKey = getProperty("gatk_key");

        intervalFilesList = Arrays.asList(StringUtils.split(getOptionalProperty("interval_files", ""), ","));
        intervalFiles = new HashSet<>(intervalFilesList);
        if (intervalFiles.size() != intervalFilesList.size()) {
            throw new RuntimeException("Duplicate interval_files detected");
        }

        bqsrCovariates = Sets.newHashSet(StringUtils.split(getProperty("bqsr_covariates"), ","));

        chrSizesList = Arrays.asList(StringUtils.split(getProperty("chr_sizes"), ","));
        chrSizes = new LinkedHashSet<>();
        for (String c : chrSizesList) {
            List<String> cc = Arrays.asList(StringUtils.split(c, "+"));
            Set<String> cc2 = new LinkedHashSet<>(cc);
            chrSizes.add(cc2);
        }

        if (chrSizes.size() != chrSizesList.size()) {
            throw new RuntimeException("Duplicate chr_sizes detected.");
        }

        //output identifiers
        outputIdentifiers = Arrays.asList(StringUtils.split(getProperty("output_identifiers"), ";"));
        if (hasDuplicates(outputIdentifiers)) {
            error("Duplicates detected in output_identifiers");
        }

        //input files
        inputFilesByGroup = getMapOfSets(outputIdentifiers, Arrays.asList(StringUtils.split(getProperty("input_files"), ";")), ",");

        //IUS-LimsKey SWIDs to link output files to
        String outputIusLimsKeysProperty = getOptionalProperty("output_ius_lims_keys", "");
        if (!outputIusLimsKeysProperty.isEmpty()) {
            outputSwidsByGroup = getMapOfSets(outputIdentifiers, Arrays.asList(StringUtils.split(outputIusLimsKeysProperty, ";")), ",");
        } else {
            outputSwidsByGroup = Collections.EMPTY_MAP;
        }

    }

    @Override
    public void setupDirectory() {
        try {
            BamMPWorkflow();
            addDirectory(dataDir);
            addDirectory(tmpDir);
            if (!dataDir.endsWith("/")) {
                dataDir += "/";
            }
            if (!tmpDir.endsWith("/")) {
                tmpDir += "/";
            }
        } catch (Exception e) {
            error("Error in setupDirectory:\n" + e);
        }
    }

    @Override
    public Map<String, SqwFile> setupFiles() {

        int id = 0;
        for (Entry<String, Set<String>> e : inputFilesByGroup.entrySet()) {
            String outputGroup = e.getKey();
            Set<String> provisionedPaths = new HashSet<>();
            for (String filePath : e.getValue()) {
                if (!"bam".equals(FilenameUtils.getExtension(filePath))) {
                    error("Unsupported input file: " + filePath);
                }
                SqwFile bam = createFile("file_in_" + id++);
                bam.setSourcePath(filePath);
                bam.setType(BAM_METATYPE);
                bam.setIsInput(true);
                provisionedPaths.add(bam.getProvisionedPath());
            }
            provisionedFilesByGroup.put(outputGroup, provisionedPaths);
        }

        return getFiles();
    }

    @Override
    public void buildWorkflow() {

        //setup operations performed on output group map
        Map<String, String> outputNameByOutputGroup = new HashMap<>();
        for (String outputIdentifier : outputIdentifiers) {
            outputNameByOutputGroup.put(outputIdentifier, outputIdentifier + ".");
        }

        //used to collect input files (and associated job) to the next step of the workflow - partitioned by the output group
        Multimap<String, Pair<String, Job>> inputFileAndJobToNextStepByOutputGroup = null;

        //for each output group: merge, sort, filter, dedup, index
        Multimap<String, Pair<String, Job>> mergedBamsByGroup = HashMultimap.create();
        for (String outputGroup : outputIdentifiers) {

            Job parentJob;
            String inputFilePath;

            if (inputFilesByGroup.get(outputGroup).size() > 1) {
                String outputName = outputNameByOutputGroup.get(outputGroup) + "merged.sorted.";
                outputNameByOutputGroup.put(outputGroup, outputName);
                String outputFilePath = dataDir + outputName + "bam";
                Job jobMergeSort = mergeSamFiles(java,
                        picardJar,
                        picardMergeMem,
                        tmpDir,
                        picardSortOrder,
                        assumeSorted,
                        picardMergeUseThreading, //use threading
                        picardMergeOtherParams,
                        outputFilePath,
                        provisionedFilesByGroup.get(outputGroup).toArray(new String[0]));
                jobMergeSort.setQueue(getOptionalProperty("queue", ""));
                parentJob = jobMergeSort;
                inputFilePath = outputFilePath;
            } else {
                String outputName = outputNameByOutputGroup.get(outputGroup) + "sorted.";
                outputNameByOutputGroup.put(outputGroup, outputName);
                String outputFilePath = dataDir + outputName + "bam";
                Job jobSort = sortSamFile(java,
                        picardJar,
                        picardMergeMem,
                        tmpDir,
                        picardSortOrder,
                        outputFilePath,
                        Iterables.getOnlyElement(provisionedFilesByGroup.get(outputGroup)),
                        picardSortOtherParams);
                jobSort.setQueue(getOptionalProperty("queue", ""));
                parentJob = jobSort;
                inputFilePath = outputFilePath;
            }

            // Filter each BAM file
            if (doFilter) {
                String outputName = outputNameByOutputGroup.get(outputGroup) + "filter.";
                outputNameByOutputGroup.put(outputGroup, outputName);
                String outputFilePath = dataDir + outputName + "bam";
                Job jobFilter = samtoolsFilterReads("SamtoolsFilter", inputFilePath, outputFilePath);
                jobFilter.addParent(parentJob);
                parentJob = jobFilter;
                inputFilePath = outputFilePath;
            }

            //deduplicate BAM file
            if (doDedup) {
                String outputName = outputNameByOutputGroup.get(outputGroup) + "deduped.";
                outputNameByOutputGroup.put(outputGroup, outputName);
                String outputFilePath = dataDir + outputName + "bam";
                Job jobDedup = markDuplicates(java,
                        picardJar,
                        picardMarkDupMem,
                        tmpDir,
                        inputFilePath,
                        outputFilePath,
                        dataDir + outputName + "metrics",
                        picardDedupOtherParams);
                jobDedup.getCommand().addArgument(REMOVE_DUPLICATES + Boolean.toString(doRemoveDups));
                jobDedup.setQueue(getOptionalProperty("queue", ""));
                jobDedup.addParent(parentJob);

                SqwFile metricFile = createOutputFile(dataDir + outputName + "metrics", TXT_METATYPE, manualOutput);
                attachCVterms(metricFile, EDAM, "plain text format (unformatted),Sequence alignment metadata");
                jobDedup.addFile(metricFile);

                //link this file to its LIMS metadata (IUS-LimsKeys)
                if (outputSwidsByGroup.containsKey(outputGroup)) {
                    metricFile.setParentAccessions(outputSwidsByGroup.get(outputGroup));
                }

                parentJob = jobDedup;
                inputFilePath = outputFilePath;
            }

            Job jobIdx = getIndexBamJob(inputFilePath);
            jobIdx.addParent(parentJob);
            parentJob = jobIdx;

            // SplitNTrim N Cigars and Reassign mapping quality 
            if (doSplitNTrim) {
                String outputName = outputNameByOutputGroup.get(outputGroup) + "split.";
                outputNameByOutputGroup.put(outputGroup, outputName);
                String outputFilePath = dataDir + outputName + "bam";
                Job jobSplitNCigar = splitNCigarReads(inputFilePath, outputFilePath);
                jobSplitNCigar.addParent(parentJob);
                parentJob = jobSplitNCigar;
                inputFilePath = outputFilePath;
            }
            mergedBamsByGroup.put(outputGroup, Pair.of(inputFilePath, parentJob));

            inputFileAndJobToNextStepByOutputGroup = mergedBamsByGroup;
        }

        if (doIndelRealignment) {
            for (String outputGroup : outputIdentifiers) {
                String outputName = outputNameByOutputGroup.get(outputGroup) + "realign.";
                outputNameByOutputGroup.put(outputGroup, outputName);
            }

            Multimap<String, Pair<String, Job>> realignedBamsByGroup = HashMultimap.create();
            for (Set<String> chrSize : chrSizes) {
                if (chrSize != null && chrSize.contains("unmapped") && chrSize.size() == 1) {
                    realignedBamsByGroup.putAll(getUnmappedBamJob(inputFileAndJobToNextStepByOutputGroup));
                } else {
                    realignedBamsByGroup.putAll(indelRealignJob(inputFileAndJobToNextStepByOutputGroup, chrSize));
                }
            }

            inputFileAndJobToNextStepByOutputGroup = realignedBamsByGroup;
        }

        if (doBQSR) {
            for (String outputGroup : outputIdentifiers) {
                String outputName = outputNameByOutputGroup.get(outputGroup) + "recal.";
                outputNameByOutputGroup.put(outputGroup, outputName);
            }

            Multimap<String, Pair<String, Job>> recalibratedBamsByGroup = baseQRecalibrateJob(inputFileAndJobToNextStepByOutputGroup);
            inputFileAndJobToNextStepByOutputGroup = recalibratedBamsByGroup;
        }

        //Merge into final output bam file (Realigned or Recalibrated reads)
        for (Entry<String, Collection<Pair<String, Job>>> e : inputFileAndJobToNextStepByOutputGroup.asMap().entrySet()) {

            String outputGroup = e.getKey();
            String outputFileName = outputNameByOutputGroup.get(outputGroup);
            String bamOutputFilePath = dataDir + outputFileName + "bam";
            String baiOutputFilePath = dataDir + outputFileName + "bai";

            String[] inputBams = getLeftCollection(inputFileAndJobToNextStepByOutputGroup.get(outputGroup)).toArray(new String[0]);

            Job jobMergeFinal = mergeSamFiles(java,
                    picardJar,
                    picardMergeMem,
                    tmpDir,
                    picardSortOrder,
                    assumeSorted,
                    picardMergeUseThreading, //use threading
                    picardMergeOtherParams,
                    bamOutputFilePath,
                    inputBams);
            jobMergeFinal.setQueue(getOptionalProperty("queue", ""));
            jobMergeFinal.getParents().addAll(getRightCollection(inputFileAndJobToNextStepByOutputGroup.get(outputGroup)));

            // Annotate and provision final bam and its index
            SqwFile finalBam = createOutputFile(bamOutputFilePath, BAM_METATYPE, manualOutput);
            if (!alignerName.isEmpty()) {
                finalBam.getAnnotations().put("aligner", alignerName);
            }
            attachCVterms(finalBam, EDAM, "BAM,Sequence alignment refinement");
            jobMergeFinal.addFile(finalBam);

            //link bam file to its LIMS metadata (IUS-LimsKeys)
            if (outputSwidsByGroup.containsKey(outputGroup)) {
                finalBam.setParentAccessions(outputSwidsByGroup.get(outputGroup));
            }

            SqwFile finalBai = createOutputFile(baiOutputFilePath, BAI_METATYPE, manualOutput);
            attachCVterms(finalBai, EDAM, "BAI");
            jobMergeFinal.addFile(finalBai);

            //link bai file to its LIMS metadata (IUS-LimsKeys)
            if (outputSwidsByGroup.containsKey(outputGroup)) {
                finalBai.setParentAccessions(outputSwidsByGroup.get(outputGroup));
            }

        }
    }

//=======================Jobs as functions===================
    /**
     * <p>
     * Filters out the reads according to the samtools flag.
     * http://picard.sourceforge.net/explain-flags.html</p>
     * <p>
     * <p>
     * Example command line</p>
     * <p>
     * {@code samtools view -b -F 260 > output.bam}
     *
     * @param jobName    the name of the samtools filter job
     * @param inputFile  the input bam file
     * @param outputFile the output bam file
     *
     * @return samtools filter reads job
     */
    protected Job samtoolsFilterReads(String jobName, String inputFile, String outputFile) {
        Job job = getWorkflow().createBashJob(jobName + samtoolsFlag);
        Command cmd = job.getCommand();
        cmd.addArgument(samtools).addArgument("view -b");
        if (samtoolsFlag != null) {
            cmd.addArgument("-F");
            cmd.addArgument(samtoolsFlag.toString());
        }
        if (samtoolsMinMapQuality != null) {
            cmd.addArgument("-q");
            cmd.addArgument(samtoolsMinMapQuality);
        }
        if (filterOtherParams != null) {
            job.getCommand().addArgument(filterOtherParams);
        }
        job.getCommand().addArgument(inputFile).addArgument(" > " + outputFile);
        job.setMaxMemory(String.valueOf(samtoolsMem));
        job.setQueue(getOptionalProperty("queue", ""));
        return job;
    }

    // splitNCigar
    protected Job splitNCigarReads(String inputFile, String outputFile) {
        /**
         * This method wraps the following command: java -Xmx16G -jar
         * $GATKROOT/GenomeAnalysisTK.jar \ -T SplitNCigarReads \ -R
         * /.mounts/labs/PDE/data/gatkAnnotationResources/hg19_random.fa \ -I
         * TGL01_0001_Pb_R_PE_466_EX.filter.dedupped.bam \ -o
         * TGL01_0001_Pb_R_PE_466_EX.filter.dedupped.split.bam \ -rf
         * ReassignOneMappingQuality \ -RMQF 255 \ -RMQT 60 \ -U
         * ALLOW_N_CIGAR_READS -fixNDN
         * <p>
         */
        Job jobSplitCigar = getWorkflow().createBashJob("split_n_trim_reassign");

        String cmd = java + " -Xmx" + (splitCigarXmxg).toString() + "G" + " -jar "
                + gatk + " -T SplitNCigarReads"
                + " -R " + refFasta
                + " -I " + inputFile
                + " -o " + outputFile;

        if (flagReassignOneMappingQuality) {
            cmd += " -rf " + "ReassignOneMappingQuality"
                    + " -RMQF " + splitCigarRMQF
                    + " -RMQT " + splitCigarRMQT;
        }

        cmd += " -U ALLOW_N_CIGAR_READS"
                + " -fixNDN";

        jobSplitCigar.setCommand(cmd);
        jobSplitCigar.setMaxMemory(Integer.toString((splitCigarXmxg + gatkOverhead) * 1024));
        jobSplitCigar.setQueue(getOptionalProperty("queue", ""));

        return jobSplitCigar;
    }

    protected Multimap<String, Pair<String, Job>> indelRealignJob(Multimap<String, Pair<String, Job>> inputFilesByGroup, Set<String> chrSizes) {

        if (chrSizes != null && !chrSizes.isEmpty() && chrSizes.contains("unmapped")) {
            error("indel realignment is not supported for \"unmapped\" - \"unmapped\" should not be grouped with other \"chr_sizes\" targets");
        }

        Multimap<String, Pair<String, Job>> realignedBamsByGroup = HashMultimap.create();

        RealignerTargetCreator.Builder realignerTargetCreatorCommandBuilder = new RealignerTargetCreator.Builder(java, gatkRealignTargetCreatorXmx + "g", tmpDir, gatk, gatkKey, dataDir)
                .setReferenceSequence(refFasta)
                .addInputBamFiles(getLeftCollection(inputFilesByGroup.values()))
                .setKnownIndels(dbsnpVcf)
                .addIntervalFiles(intervalFiles)
                .setIntervalPadding(intervalPadding)
                .setDownsamplingCoverageThreshold(downsamplingCoverage)
                .setDownsamplingType(downsamplingType)
                .setExtraParameters(realignerTargetCreatorParams);
        for (String chrSize : chrSizes) {
            realignerTargetCreatorCommandBuilder.addInterval(chrSize);
        }

        RealignerTargetCreator realignerTargetCreatorCommand = realignerTargetCreatorCommandBuilder.build();

        Job realignerTargetCreatorJob = getWorkflow().createBashJob("GATKRealignerTargetCreator")
                .setMaxMemory(Integer.toString((gatkRealignTargetCreatorXmx + gatkOverhead) * 1024))
                .setQueue(queue);
        realignerTargetCreatorJob.getCommand().setArguments(realignerTargetCreatorCommand.getCommand());
        realignerTargetCreatorJob.getParents().addAll(getRightCollection(inputFilesByGroup.values()));

        Multimap<String, String> inputBamFilesByGroup = HashMultimap.create();
        for (Entry<String, Collection<Pair<String, Job>>> e : inputFilesByGroup.asMap().entrySet()) {
            for (Pair<String, Job> p : e.getValue()) {
                inputBamFilesByGroup.put(e.getKey(), p.getLeft());
            }
        }

        IndelRealigner.Builder indelRealignerCommandBuilder = new IndelRealigner.Builder(java, gatkIndelRealignerXmx + "g", tmpDir, gatk, gatkKey, dataDir)
                .setReferenceSequence(refFasta)
                .addInputBamFiles(inputBamFilesByGroup)
                .addKnownIndelFile(dbsnpVcf)
                .addIntervalFiles(intervalFiles)
                .setIntervalPadding(intervalPadding)
                .setTargetIntervalFile(realignerTargetCreatorCommand.getOutputFile())
                .setExtraParameters(indelRealignerParams);
        for (String chrSize : chrSizes) {
            indelRealignerCommandBuilder.addInterval(chrSize);
        }

        IndelRealigner indelRealignerCommand = indelRealignerCommandBuilder.build();

        Job indelRealignerJob = getWorkflow().createBashJob("GATKIndelRealigner")
                .setMaxMemory(Integer.toString((gatkIndelRealignerXmx + gatkOverhead) * 1024))
                .setQueue(queue)
                .addParent(realignerTargetCreatorJob);
        indelRealignerJob.getCommand().setArguments(indelRealignerCommand.getCommand());

        for (Entry<String, Collection<String>> e : indelRealignerCommand.getOutputFilesByGroup().entrySet()) {
            String outputGroup = e.getKey();
            for (String outputFile : e.getValue()) {
                realignedBamsByGroup.put(outputGroup, Pair.of(outputFile, indelRealignerJob));
            }
        }

        return realignedBamsByGroup;
    }

    protected Multimap<String, Pair<String, Job>> baseQRecalibrateJob(Multimap<String, Pair<String, Job>> realignedBamsByGroup) {

        //GATK Base Recalibrator ( https://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_gatk_tools_walkers_bqsr_BaseRecalibrator.php )
        BaseRecalibrator baseRecalibratorCommand = new BaseRecalibrator.Builder(java, gatkBaseRecalibratorXmx + "m", tmpDir, gatk, gatkKey, dataDir)
                .setReferenceSequence(refFasta)
                .setCovariates(bqsrCovariates)
                .addKnownSite(dbsnpVcf)
                .addInputFiles(getLeftCollection(realignedBamsByGroup.values()))
                .addIntervalFiles(intervalFiles)
                .setIntervalPadding(intervalPadding)
                .setNumCpuThreadsPerDataThread(gatkBaseRecalibratorNct)
                .setExtraParameters(baseRecalibratorParams)
                .build();
        Job baseRecalibratorJob = getWorkflow().createBashJob("GATKBaseRecalibrator")
                .setMaxMemory(gatkBaseRecalibratorMem.toString())
                .setThreads(gatkBaseRecalibratorSmp)
                .setQueue(queue);
        baseRecalibratorJob.getParents().addAll(getRightCollection(realignedBamsByGroup.values()));
        baseRecalibratorJob.getCommand().setArguments(baseRecalibratorCommand.getCommand());

        SqwFile recalibrationData = createOutputFile(baseRecalibratorCommand.getOutputFile(), TXT_METATYPE, manualOutput);
        attachCVterms(recalibrationData, EDAM, "plain text format (unformatted),Read pre-processing,Sequence alignment refinement,Sequence alignment metadata");
        baseRecalibratorJob.addFile(recalibrationData);

        //GATK Analyze Covariates ( https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_bqsr_AnalyzeCovariates.php )
        AnalyzeCovariates analyzeCovariatesCommand = new AnalyzeCovariates.Builder(java, "4g", tmpDir, gatk, gatkKey, rDir, dataDir)
                .setReferenceSequence(refFasta)
                .setRecalibrationTable(baseRecalibratorCommand.getOutputFile())
                //setOutputFileName(outputName)
                .setExtraParameters(analyzeCovariatesParams)
                .build();
        Job analyzeCovariatesJob = getWorkflow().createBashJob("GATKAnalyzeCovariates")
                .setMaxMemory(Integer.toString((4 + gatkOverhead) * 1024))
                .setQueue(queue)
                .addParent(baseRecalibratorJob);
        analyzeCovariatesJob.getCommand().setArguments(analyzeCovariatesCommand.getCommand());

        SqwFile recalibrationReport = createOutputFile(analyzeCovariatesCommand.getPlotsReportFile(), PDF_METATYPE, manualOutput);
        attachCVterms(recalibrationReport, EDAM, "Read pre-processing,Sequence alignment refinement");
        analyzeCovariatesJob.addFile(recalibrationReport);

        Multimap<String, Pair<String, Job>> recalibratedBams = HashMultimap.create();
        for (Entry<String, Collection<Pair<String, Job>>> e : realignedBamsByGroup.asMap().entrySet()) {

            String outputGroup = e.getKey();
            Collection<String> inputFiles = getLeftCollection(e.getValue());

            for (String inputFile : inputFiles) {

                PrintReads printReadsCommand = new PrintReads.Builder(java, gatkPrintReadsXmx + "g", tmpDir, gatk, gatkKey, dataDir)
                        .setReferenceSequence(refFasta)
                        .setPreserveQscoresLessThan(preserveQscoresLessThan)
                        .setCovariatesTablesFile(baseRecalibratorCommand.getOutputFile())
                        .addInputFile(inputFile)
                        .setIntervalPadding(intervalPadding)
                        .setExtraParameters(printReadsParams).build();

                Job printReadsJob = getWorkflow().createBashJob("GATKTableRecalibration")
                        .setMaxMemory(Integer.toString((gatkPrintReadsXmx + gatkOverhead) * 1024))
                        .setQueue(queue);

                printReadsJob.addParent(baseRecalibratorJob);
                printReadsJob.getCommand().setArguments(printReadsCommand.getCommand());

                recalibratedBams.put(outputGroup, Pair.of(printReadsCommand.getOutputFile(), printReadsJob));
            }
        }

        return recalibratedBams;
    }

    protected Multimap<String, Pair<String, Job>> getUnmappedBamJob(Multimap<String, Pair<String, Job>> inputFilesByGroup) {
        Multimap<String, Pair<String, Job>> unmappedFilesByGroup = HashMultimap.create();
        for (Entry<String, Collection<Pair<String, Job>>> e : inputFilesByGroup.asMap().entrySet()) {
            String outputGroup = e.getKey();
            String inputFile = Iterables.getOnlyElement(e.getValue()).getLeft();
            Job parentJob = Iterables.getOnlyElement(e.getValue()).getRight();

            PrintReads printReadsCommand;
            printReadsCommand = new PrintReads.Builder(java, gatkPrintReadsXmx + "g", tmpDir, gatk, gatkKey, dataDir)
                    .setReferenceSequence(refFasta)
                    .addInputFile(inputFile)
                    .addInterval("unmapped")
                    .setExtraParameters(printReadsParams)
                    .build();

            Job printReadsJob = getWorkflow().createBashJob("GATKPrintReadsUnmapped")
                    .setMaxMemory(Integer.toString((gatkPrintReadsXmx + gatkOverhead) * 1024))
                    .setQueue(queue);
            printReadsJob.addParent(parentJob);
            printReadsJob.getCommand().setArguments(printReadsCommand.getCommand());

            unmappedFilesByGroup.put(outputGroup, Pair.of(printReadsCommand.getOutputFile(), printReadsJob));
        }

        return unmappedFilesByGroup;
    }

    /**
     * Index will be written to the same directory as the bam
     *
     * @param inputFilePath the input bam file
     *
     * @return index bam file job
     */
    protected Job getIndexBamJob(String inputFilePath) {
        String outputFilePath = FilenameUtils.getPath(inputFilePath) + FilenameUtils.getBaseName(inputFilePath) + ".bai";

        Job jobIndex = getWorkflow().createBashJob("index_bam");
        jobIndex.setCommand(java + " -Xmx3G -jar "
                + picardJar + " BuildBamIndex"
                + " VALIDATION_STRINGENCY=SILENT"
                + " I=" + inputFilePath
                + " O=" + outputFilePath);
        jobIndex.setMaxMemory(getProperty("index_mem"));
        jobIndex.setQueue(getOptionalProperty("queue", ""));

        return jobIndex;
    }

    private <T, S> Set<T> getLeftCollection(Collection<Pair<T, S>> pairs) {
        Set<T> ts = new HashSet<>();
        for (Pair<T, S> p : pairs) {
            ts.add(p.getLeft());
        }
        return ts;
    }

    private <S, T> Set<T> getRightCollection(Collection<Pair<S, T>> pairs) {
        Set<T> ts = new HashSet<>();
        for (Pair<S, T> p : pairs) {
            ts.add(p.getRight());
        }
        return ts;
    }

    private Map<String, Set<String>> getMapOfSets(List<String> keys, List<String> values, String valuesDelimiter) {
        if (keys.isEmpty() || values.isEmpty() || keys.size() != values.size()) {
            error("key size (" + keys.size() + ") != values size (" + values.size() + ")");
        }
        Map<String, Set<String>> map = new HashMap<>();
        for (int i = 0; i < keys.size(); i++) {
            List<String> valsList = Arrays.asList(values.get(i).split(valuesDelimiter));
            Set<String> vals = new LinkedHashSet<>(valsList);
            if (map.put(keys.get(i), vals) != null) {
                error("Duplicate key detected = [" + keys.get(i));
            }
        }
        return map;
    }

    private <T> boolean hasDuplicates(Collection<T> c) {
        Set<T> s = new HashSet<>(c);
        return c.size() != s.size();
    }

    private void error(String msg) {
        Log.error(msg);
        setWorkflowInvalid();
    }

    public Job markDuplicates(String java, String picardJar, int memoryMb, String tmpDir,
            String input, String output, String metricsFile, String otherParams) {
        String command = String.format("%s -Xmx%dM -jar %s MarkDuplicates INPUT=%s OUTPUT=%s "
                + "VALIDATION_STRINGENCY=SILENT TMP_DIR=%s METRICS_FILE=%s",
                java, memoryMb, picardJar, input, output, tmpDir, metricsFile, otherParams);
        Job job = getWorkflow().createBashJob("PicardMarkDuplicates");
        job.getCommand().addArgument(command);
        if (otherParams != null) {
            job.getCommand().addArgument(otherParams);
        }
        job.setMaxMemory(String.valueOf(memoryMb + picardMemOverhead));
        return job;
    }

    public Job mergeSamFiles(String java, String picardJar, int memoryMb, String tmpDir,
            String sortOrder, boolean assumeSorted, boolean useThreading, String otherParams,
            String output, String... input) {
        String command = String.format("java -Xmx%dM -jar %s MergeBamAlignment "
                + "OUTPUT=%s "
                + "VALIDATION_STRINGENCY=SILENT "
                + "TMP_DIR=%s "
                + "SORT_ORDER=%s "
                + "CREATE_INDEX=true",
                memoryMb, picardJar, output, tmpDir, sortOrder);
        Job job = getWorkflow().createBashJob("PicardMergeBam");
        job.getCommand().addArgument(command);
        for (String in : input) {
            job.getCommand().addArgument(String.format("INPUT=%s", in));
        }
        if (assumeSorted) {
            job.getCommand().addArgument(String.format("ASSUME_SORTED=true"));
        }
        if (useThreading) {
            job.getCommand().addArgument(String.format("USE_THREADING=true"));
        }
        if (otherParams != null) {
            job.getCommand().addArgument(otherParams);
        }
        job.setMaxMemory(String.valueOf(memoryMb + picardMemOverhead));
        return job;
    }

    public Job sortSamFile(String java, String picardJar, int memoryMb,
            String tmpDir, String sortOrder, String output, String input, String otherParams) {
        String command = String.format("java -Xmx%dM -jar %s SortSam "
                + "OUTPUT=%s "
                + "INPUT=%s "
                + "VALIDATION_STRINGENCY=SILENT "
                + "TMP_DIR=%s "
                + "SORT_ORDER=%s "
                + "CREATE_INDEX=true",
                memoryMb, picardJar, output, input, tmpDir, sortOrder);
        Job job = getWorkflow().createBashJob("PicardSortBam");
        job.getCommand().addArgument(command);
        if (otherParams != null) {
            job.getCommand().addArgument(otherParams);
        }
        job.setMaxMemory(String.valueOf(memoryMb + picardMemOverhead));
        return job;
    }

}
