package ca.on.oicr.pde.workflows;

import ca.on.oicr.pde.tools.gatk3.AnalyzeCovariates;
import ca.on.oicr.pde.tools.gatk3.BaseRecalibrator;
import ca.on.oicr.pde.tools.gatk3.IndelRealigner;
import ca.on.oicr.pde.tools.gatk3.PrintReads;
import ca.on.oicr.pde.tools.gatk3.PrintReadsUnmapped;
import ca.on.oicr.pde.tools.gatk3.RealignerTargetCreator;
import ca.on.oicr.pde.utilities.workflows.SemanticWorkflow;
import ca.on.oicr.pde.utilities.workflows.jobfactory.PicardTools;
import com.google.common.collect.HashMultimap;
import com.google.common.collect.Multimap;
import com.google.common.collect.Sets;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import net.sourceforge.seqware.common.module.ReturnValue;
import net.sourceforge.seqware.common.util.Log;
import net.sourceforge.seqware.pipeline.workflowV2.model.Command;
import net.sourceforge.seqware.pipeline.workflowV2.model.Job;
import net.sourceforge.seqware.pipeline.workflowV2.model.SqwFile;
import net.sourceforge.seqware.pipeline.workflowV2.model.Workflow;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.lang3.StringUtils;
import org.apache.commons.lang3.tuple.Pair;

/**
 *
 * Merge -> Sort -> Mark Duplicates (optional q-score filtering) -> Realign
 * Indels -> Base Q score recalibration
 *
 */
public class BamMPWorkflow extends SemanticWorkflow {

    private final List<String> inputBamFiles = new LinkedList<>();
    private Multimap<String, Pair<String, Job>> realignedBams;
    private Multimap<String, Pair<String, Job>> recalibratedBams;
    private String finalOutput;

    private String dataDir, tmpDir, rDir;
    private String java, markDuplicatesJar, mergeSamFilesJar, sortSamFilesJar;
    private Integer picardMarkDupMem, picardMergeMem;
    // we should have assumeSorted set to false ALWAYS (due to re-structuring)
    private boolean doDedup = true, doRemoveDups = true, doFilter = true, assumeSorted = false, useThreading = true, manualOutput = false;
    private String samtools;
    private String picard_dir;
    private Integer samtoolsMem;
    private Integer samtoolsFlag;
    private String samtoolsMinMapQuality;
    private String identifier;
    private String alignerName;
    private Workflow wf;
    private PicardTools picard;
    private String filterOtherParams, dedupOtherParams, mergeOtherParams, sortOtherParams;
    private String sortOrder;
    private static final String REMOVE_DUPLICATES = "REMOVE_DUPLICATES=";
    private final static String BAM_METATYPE = "application/bam";
    private final static String BAI_METATYPE = "application/bam-index";
    private final static String PDF_METATYPE = "application/pdf";
    private final static String TXT_METATYPE = "text/plain";

    private List<String> chrSizesList;
    private Set<String> chrSizes;

    //References
    private String refFasta;
    private String dbsnpVcf;

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

    static {
        cvTerms = new HashMap<String, Set<String>>();
        cvTerms.put(EDAM, new HashSet<String>(Arrays.asList("BAM", "BAI", "plain text format (unformatted)",
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
        Map<String, Set<String>> myTerms = new HashMap<String, Set<String>>();
        myTerms.putAll(cvTerms);
        return myTerms;
    }

    /**
     * Launched from setupDirectory();
     */
    private void BamMPWorkflow() {
        identifier = getProperty("identifier");
        java = getProperty("java");
        dataDir = "data";
        tmpDir = getProperty("tmp_dir");
        rDir = getProperty("r_dir");
        manualOutput = Boolean.parseBoolean(getProperty("manual_output"));

        //Flags for controlling different steps:
        doFilter = Boolean.parseBoolean(getProperty("do_sam_filter"));
        doRemoveDups = Boolean.parseBoolean(getOptionalProperty("do_remove_duplicates", "false"));
        doDedup = doRemoveDups ? true : Boolean.parseBoolean(getProperty("do_mark_duplicates"));

        //Picard
        markDuplicatesJar = getProperty("picard_mark_duplicates");
        picardMarkDupMem = Integer.parseInt(getProperty("picard_mark_duplicates_mem_mb"));
        mergeSamFilesJar = getProperty("picard_merge_sam");
        sortSamFilesJar = getProperty("picard_sort_sam");
        picardMergeMem = Integer.parseInt(getProperty("picard_merge_sam_mem_mb"));
        filterOtherParams = getOptionalProperty("samtools_filter_other_params", "");
        dedupOtherParams = getOptionalProperty("picard_mark_duplicates_other_params", "");
        mergeOtherParams = getOptionalProperty("picard_merge_other_params", "");
        analyzeCovariatesParams = getOptionalProperty("gatk_analyze_covariates_params", null);
        picard_dir = getProperty("picard_dir");
        wf = this.getWorkflow();
        picard = new PicardTools(wf);
        binDir = this.getWorkflowBaseDir() + "/bin/";
        doBQSR = Boolean.valueOf(getOptionalProperty("do_bqsr", "true"));
        //Samtools
        samtools = getProperty("samtools");
        samtoolsMem = Integer.parseInt(getProperty("samtools_mem_mb"));
        samtoolsFlag = Integer.parseInt(getProperty("samtools_filter_flag"));
        samtoolsMinMapQuality = getOptionalProperty("samtools_min_map_quality", null);

        alignerName = getOptionalProperty("aligner_name", "");
        sortOtherParams = getOptionalProperty("picard_sort_other_params", "");
        sortOrder = getProperty("sort_order");
        useThreading = Boolean.parseBoolean(getProperty("picard_merge_use_threading"));

        this.chrSizesList = Arrays.asList(StringUtils.split(getProperty("chr_sizes"), ","));
        this.chrSizes = new HashSet<>(chrSizesList);
        if (chrSizes.size() != chrSizesList.size()) {
            throw new RuntimeException("Duplicate chr_sizes detected.");
        }

    }

    @Override
    public Map<String, SqwFile> setupFiles() {

        //References
        this.refFasta = getProperty("ref_fasta");
        this.dbsnpVcf = getProperty("gatk_dbsnp_vcf");

        //GATK
        this.gatkRealignTargetCreatorXmx = Integer.parseInt(getProperty("gatk_realign_target_creator_xmx"));
        this.gatkIndelRealignerXmx = Integer.parseInt(getProperty("gatk_indel_realigner_xmx"));
        this.gatkPrintReadsXmx = Integer.parseInt(getProperty("gatk_print_reads_xmx"));
        this.gatkBaseRecalibratorXmx = Integer.parseInt(getProperty("gatk_base_recalibrator_xmx"));
        this.gatkBaseRecalibratorMem = Integer.parseInt(getProperty("gatk_base_recalibrator_mem"));
        this.gatkBaseRecalibratorNct = Integer.parseInt(getProperty("gatk_base_recalibrator_nct"));
        this.gatkBaseRecalibratorSmp = Integer.parseInt(getProperty("gatk_base_recalibrator_smp"));
        this.gatkOverhead = Integer.parseInt(getProperty("gatk_sched_overhead_mem"));
        this.downsamplingCoverage = hasPropertyAndNotNull("downsampling_coverage") ? Integer.parseInt(getProperty("downsampling_coverage")) : null;
        this.preserveQscoresLessThan = hasPropertyAndNotNull("preserve_qscores_less_than") ? Integer.parseInt(getProperty("preserve_qscores_less_than")) : null;
        this.intervalPadding = hasPropertyAndNotNull("interval_padding") ? Integer.parseInt(getProperty("interval_padding")) : null;

        this.queue = getOptionalProperty("queue", "");
        this.realignerTargetCreatorParams = getOptionalProperty("gatk_realigner_target_creator_params", null);
        this.indelRealignerParams = getOptionalProperty("gatk_indel_realigner_params", null);
        this.baseRecalibratorParams = getOptionalProperty("gatk_base_recalibrator_params", null);
        this.printReadsParams = getOptionalProperty("gatk_print_reads_params", null);
        this.downsamplingType = getOptionalProperty("downsampling_type", null);
        this.gatk = getOptionalProperty("gatk_jar", binDir);
        this.gatkKey = getProperty("gatk_key");

        this.intervalFilesList = Arrays.asList(StringUtils.split(getOptionalProperty("interval_files", ""), ","));
        this.intervalFiles = new HashSet<>(intervalFilesList);
        this.bqsrCovariates = Sets.newHashSet(StringUtils.split(getProperty("bqsr_covariates"), ","));

        List<String> inputFilesList = Arrays.asList(StringUtils.split(getProperty("input_files"), ","));
        Set<String> inputFilesSet = new HashSet<>(inputFilesList);

        if (inputFilesList.size() != inputFilesSet.size()) {
            throw new RuntimeException("Duplicate files detected in input_files");
        }

        Map<String, String> bams = new HashMap<>();

        for (String f : inputFilesSet) {
            String fileExtension = FilenameUtils.getExtension(f);
            String fileKey = FilenameUtils.removeExtension(f);
            if (null != fileExtension) {
                switch (fileExtension) {
                    case "bam":
                        bams.put(fileKey, f);
                        break;
                    default:
                        throw new RuntimeException("Unsupported input file type");
                }
            }
        }

        int id = 0;
        for (Map.Entry<String, String> e : bams.entrySet()) {
            String bamFilePath = e.getValue();

            SqwFile bam = this.createFile("file_in_" + id++);
            bam.setSourcePath(bamFilePath);
            bam.setType(BAM_METATYPE);
            bam.setIsInput(true);

            inputBamFiles.add(bam.getProvisionedPath());
        }

        return this.getFiles();
    }

    @Override
    public void setupDirectory() {
        try {
            BamMPWorkflow();
            this.addDirectory(dataDir);
            this.addDirectory(tmpDir);
            if (!dataDir.endsWith("/")) {
                dataDir += "/";
            }
            if (!tmpDir.endsWith("/")) {
                tmpDir += "/";
            }
        } catch (Exception e) {
            e.printStackTrace();
            Log.error("Error in setupDirectory", e);
            ret.setReturnValue(ReturnValue.INVALIDFILE);
        }
    }

    @Override
    public void buildWorkflow() {

        if (intervalFiles.size() != intervalFilesList.size()) {
            throw new RuntimeException("Duplicate interval_files detected");
        }
        // use these two strings to construct the file names for the intermediate
        // files. This way the file name describes the operations performed on it
        String operationsOnMergedFile = identifier + ".";
        String outputFile = "";
        String inputFile = outputFile;
        List<Job> upstreamJobs = new ArrayList<Job>();
        Job jobMergeSort;

        // Merge + sort BAM files - move into a function
        if (inputBamFiles.size() > 1) {
            operationsOnMergedFile += "merged.sorted.";
            outputFile = dataDir + "/" + operationsOnMergedFile;
            String[] bamProvisionedPath = new String[inputBamFiles.size()];
            for (int s = 0; s < inputBamFiles.size(); s++) {
                bamProvisionedPath[s] = inputBamFiles.get(s).toString();
            }

            jobMergeSort = picard.mergeSamFiles(this.java,
                    mergeSamFilesJar,
                    picardMergeMem,
                    tmpDir,
                    sortOrder,
                    assumeSorted,
                    useThreading, //use threading
                    mergeOtherParams,
                    outputFile + "bam",
                    bamProvisionedPath);
            jobMergeSort.setQueue(getOptionalProperty("queue", ""));
            inputFile = outputFile;
            upstreamJobs.add(jobMergeSort);
        } else if (inputBamFiles.size() == 0) {
            Log.fatal("There are no BAM files to launch!");
            ret.setExitStatus(ReturnValue.INVALIDPARAMETERS);
        } else if (inputBamFiles.size() == 1) {
            inputFile = inputBamFiles.get(0).toString();
            // Sort BAM file this is mandatory
            operationsOnMergedFile += "sorted.";
            outputFile = dataDir + "/" + operationsOnMergedFile;
            Job jobSort = picard.sortSamFile(this.java,
                    sortSamFilesJar,
                    picardMergeMem,
                    tmpDir,
                    sortOrder,
                    outputFile + "bam",
                    inputFile,
                    sortOtherParams);
            jobSort.setQueue(getOptionalProperty("queue", ""));
            inputFile = outputFile;
            upstreamJobs.add(jobSort);
        }

        // Clean merged/sorted BAM
        String dir = inputFile.substring(0, inputFile.lastIndexOf("/")) + "/";
        Job upstreamJob = upstreamJobs.isEmpty() ? null : (Job) upstreamJobs.get(0); // There suppsed to be one!

        // Filter each BAM file
        if (doFilter) {
            operationsOnMergedFile += "filter.";
            outputFile = dir + operationsOnMergedFile;
            Job jobFilter = samtoolsFilterReads("SamtoolsFilter", inputFile + "bam", outputFile + "bam");
            if (null != upstreamJob) {
                jobFilter.addParent(upstreamJob);
            }
            upstreamJob = jobFilter;
            inputFile = outputFile;
        }

        //deduplicate BAM file
        if (doDedup) {
            operationsOnMergedFile += "deduped.";
            outputFile = dir + operationsOnMergedFile;
            Job jobDedup = picard.markDuplicates(this.java,
                    markDuplicatesJar,
                    picardMarkDupMem,
                    tmpDir,
                    inputFile + "bam",
                    outputFile + "bam",
                    dir + operationsOnMergedFile + "metrics",
                    dedupOtherParams);
            jobDedup.getCommand().addArgument(REMOVE_DUPLICATES + Boolean.toString(this.doRemoveDups));
            jobDedup.setQueue(getOptionalProperty("queue", ""));
            SqwFile metricFile = this.createOutputFile(dir + operationsOnMergedFile + "metrics", TXT_METATYPE, this.manualOutput);

            this.attachCVterms(metricFile, EDAM, "plain text format (unformatted),Sequence alignment metadata");
            jobDedup.addFile(metricFile);
            if (null != upstreamJob) {
                jobDedup.addParent(upstreamJob);
            }
            upstreamJob = jobDedup;
            inputFile = outputFile;
        }

        Job jobIdx = this.indexBamJob(inputFile);

        if (null != upstreamJob) {
            jobIdx.addParent(upstreamJob);
        }

        operationsOnMergedFile += "realigned.";

        // Indel Realignment Job
        LinkedList<String> inputs = new LinkedList();
        inputs.add(inputFile + "bam");
        this.recalibratedBams = HashMultimap.create();
        this.realignedBams = HashMultimap.create();
        Multimap<String, Pair<String, Job>> inputBams;

        // We split by chromosome here
        for (String chrSize : chrSizes) {
            if ("unmapped".equals(chrSize)) {
                Pair<String, Job> unmappedBamJob = getUnmappedBamJob(jobIdx, inputs);
                realignedBams.put(chrSize, unmappedBamJob);
            } else {
                this.indelRealignJob(inputs, chrSize, jobIdx);
            }
        }

        if (doBQSR) {
            //Conditional Recalibration job
            operationsOnMergedFile += "recal.";
            this.baseQRecalibrateJob(operationsOnMergedFile);
            inputBams = this.recalibratedBams;

        } else {
            inputBams = this.realignedBams;
        }

        // Use picard for indexing
        this.finalOutput = this.dataDir + operationsOnMergedFile;
        Object[] finalInputs = this.getLeftCollection(inputBams.values()).toArray();
        String[] inputBamsArray = new String[finalInputs.length];
        for (int s = 0; s < finalInputs.length; s++) {
            inputBamsArray[s] = finalInputs[s].toString();
        }

        //Merge into final output bam file (Realigned or Recalibrated reads)
        Job jobMergeFinal = picard.mergeSamFiles(this.java,
                mergeSamFilesJar,
                picardMergeMem,
                tmpDir,
                sortOrder,
                assumeSorted,
                useThreading, //use threading
                mergeOtherParams,
                this.finalOutput + "bam",
                inputBamsArray);
        jobMergeFinal.setQueue(getOptionalProperty("queue", ""));
        jobMergeFinal.getParents().addAll(this.getRightCollection(inputBams.values()));

        // Annotate and provision final bam and its index
        SqwFile finalBam = this.createOutputFile(this.finalOutput + "bam", BAM_METATYPE, manualOutput);
        if (!this.alignerName.isEmpty()) {
            finalBam.getAnnotations().put("aligner", this.alignerName);
        }

        this.attachCVterms(finalBam, EDAM, "BAM,Sequence alignment refinement");
        jobMergeFinal.addFile(finalBam);

        Job jobIdx2 = this.indexBamJob(this.finalOutput);
        jobIdx2.addParent(jobMergeFinal);

        SqwFile finalBai = this.createOutputFile(this.finalOutput + "bai", BAI_METATYPE, manualOutput);
        this.attachCVterms(finalBai, EDAM, "BAI");
        jobIdx2.addFile(finalBai);
    }

    //=======================Jobs as functions===================
    /**
     * <p>
     * Filters out the reads according to the samtools flag.
     * http://picard.sourceforge.net/explain-flags.html</p>
     *
     * <p>
     * Example command line</p>
     *
     * <code>samtools view -b -F 260 > output.bam</code>
     *
     * @param jobName    the name of the samtools filter job
     * @param inputFile  the input bam file
     * @param outputFile the output bam file
     *
     * @return
     */
    protected Job samtoolsFilterReads(String jobName, String inputFile, String outputFile) {
        Job job = this.getWorkflow().createBashJob(jobName + samtoolsFlag);
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

    protected void indelRealignJob(LinkedList<String> inputFile, String chrSize, Job parentJob) {

        RealignerTargetCreator realignerTargetCreatorCommand = new RealignerTargetCreator.Builder(java, gatkRealignTargetCreatorXmx + "g", tmpDir, gatk, gatkKey, dataDir)
                .setReferenceSequence(refFasta)
                .addInputBamFiles(inputFile)
                .setKnownIndels(dbsnpVcf)
                .addInterval(chrSize)
                .addIntervalFiles(intervalFiles)
                .setIntervalPadding(intervalPadding)
                .setDownsamplingCoverageThreshold(downsamplingCoverage)
                .setDownsamplingType(downsamplingType)
                .setOutputFileName("gatk" + (chrSize != null ? "." + chrSize.replace(":", "-") : ""))
                .setExtraParameters(realignerTargetCreatorParams)
                .build();
        Job realignerTargetCreatorJob = getWorkflow().createBashJob("GATKRealignerTargetCreator")
                .setMaxMemory(Integer.toString((gatkRealignTargetCreatorXmx + gatkOverhead) * 1024))
                .setQueue(queue);
        realignerTargetCreatorJob.getCommand().setArguments(realignerTargetCreatorCommand.getCommand());
        realignerTargetCreatorJob.addParent(parentJob);

        IndelRealigner indelRealignerCommand = new IndelRealigner.Builder(java, gatkIndelRealignerXmx + "g", tmpDir, gatk, gatkKey, dataDir)
                .setReferenceSequence(refFasta)
                .addInputBamFiles(inputFile)
                .addKnownIndelFile(dbsnpVcf)
                .addInterval(chrSize)
                .addIntervalFiles(intervalFiles)
                .setIntervalPadding(intervalPadding)
                .setTargetIntervalFile(realignerTargetCreatorCommand.getOutputFile())
                .setExtraParameters(indelRealignerParams)
                .build();
        Job indelRealignerJob = getWorkflow().createBashJob("GATKIndelRealigner")
                .setMaxMemory(Integer.toString((gatkIndelRealignerXmx + gatkOverhead) * 1024))
                .setQueue(this.queue)
                .addParent(realignerTargetCreatorJob);
        indelRealignerJob.getCommand().setArguments(indelRealignerCommand.getCommand());

        if (realignedBams.containsKey(chrSize)) {
            throw new RuntimeException("Unexpected state: Duplicate interval key.");
        }
        for (String outputFile : indelRealignerCommand.getOutputFiles()) {
            realignedBams.put(chrSize, Pair.of(outputFile, indelRealignerJob));
        }

    }

    protected void baseQRecalibrateJob(String OperationOnFile) {

        //GATK Base Recalibrator ( https://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_gatk_tools_walkers_bqsr_BaseRecalibrator.php )
        BaseRecalibrator baseRecalibratorCommand = new BaseRecalibrator.Builder(java, gatkBaseRecalibratorXmx + "m", tmpDir, gatk, gatkKey, dataDir)
                .setReferenceSequence(refFasta)
                .setCovariates(bqsrCovariates)
                .addKnownSite(dbsnpVcf)
                .addInputFiles(getLeftCollection(realignedBams.values()))
                .addIntervalFiles(intervalFiles)
                .setIntervalPadding(intervalPadding)
                .setNumCpuThreadsPerDataThread(gatkBaseRecalibratorNct)
                .setExtraParameters(baseRecalibratorParams)
                .build();
        Job baseRecalibratorJob = getWorkflow().createBashJob("GATKBaseRecalibrator")
                .setMaxMemory(gatkBaseRecalibratorMem.toString())
                .setThreads(gatkBaseRecalibratorSmp)
                .setQueue(queue);
        baseRecalibratorJob.getParents().addAll(getRightCollection(realignedBams.values()));
        baseRecalibratorJob.getCommand().setArguments(baseRecalibratorCommand.getCommand());

        SqwFile recalibrationData = createOutputFile(baseRecalibratorCommand.getOutputFile(), TXT_METATYPE, manualOutput);
        this.attachCVterms(recalibrationData, EDAM, "plain text format (unformatted),Read pre-processing,Sequence alignment refinement,Sequence alignment metadata");
        baseRecalibratorJob.addFile(recalibrationData);

        //GATK Analyze Covariates ( https://www.broadinstitute.org/gatk/guide/tooldocs/org_broadinstitute_gatk_tools_walkers_bqsr_AnalyzeCovariates.php )
        String outputName = OperationOnFile.endsWith(".") ? OperationOnFile.substring(0, OperationOnFile.lastIndexOf(".")) : OperationOnFile;
        AnalyzeCovariates analyzeCovariatesCommand = new AnalyzeCovariates.Builder(java, "4g", tmpDir, gatk, gatkKey, rDir, dataDir)
                .setReferenceSequence(refFasta)
                .setRecalibrationTable(baseRecalibratorCommand.getOutputFile())
                .setOutputFileName(outputName)
                .setExtraParameters(analyzeCovariatesParams)
                .build();
        Job analyzeCovariatesJob = getWorkflow().createBashJob("GATKAnalyzeCovariates")
                .setMaxMemory(Integer.toString((4 + gatkOverhead) * 1024))
                .setQueue(queue)
                .addParent(baseRecalibratorJob);
        analyzeCovariatesJob.getCommand().setArguments(analyzeCovariatesCommand.getCommand());

        SqwFile recalibrationReport = createOutputFile(analyzeCovariatesCommand.getPlotsReportFile(), PDF_METATYPE, manualOutput);
        this.attachCVterms(recalibrationReport, EDAM, "Read pre-processing,Sequence alignment refinement");
        analyzeCovariatesJob.addFile(recalibrationReport);

        for (Map.Entry<String, Pair<String, Job>> e : realignedBams.entries()) {

            String chrSize = e.getKey();
            String inputBam = e.getValue().getLeft();

            PrintReads.Builder printReadsCommandBuilder;
            printReadsCommandBuilder = new PrintReads.Builder(java, gatkPrintReadsXmx + "g", tmpDir, gatk, gatkKey, dataDir)
                    .setReferenceSequence(refFasta)
                    .setPreserveQscoresLessThan(preserveQscoresLessThan)
                    .setCovariatesTablesFile(baseRecalibratorCommand.getOutputFile())
                    .addInputFile(inputBam)
                    .setIntervalPadding(intervalPadding)
                    .setExtraParameters(printReadsParams);

            //workaround for GATK error when -L unmapped and an empty input bam is provided
            if ("unmapped".equals(chrSize)) {
                // do not set -L/--intervals for "unmapped"
            } else {
                printReadsCommandBuilder.addInterval(chrSize);
            }

            PrintReads printReadsCommand = printReadsCommandBuilder.build();

            Job printReadsJob = getWorkflow().createBashJob("GATKTableRecalibration")
                    .setMaxMemory(Integer.toString((gatkPrintReadsXmx + gatkOverhead) * 1024))
                    .setQueue(queue);

            printReadsJob.addParent(analyzeCovariatesJob);
            printReadsJob.getCommand().setArguments(printReadsCommand.getCommand());
            recalibratedBams.put(chrSize, Pair.of(printReadsCommand.getOutputFile(), printReadsJob));
        }
    }

    protected Pair<String, Job> getUnmappedBamJob(Job parent, List<String> inputBams) {
        PrintReadsUnmapped printReadsCommand;
        printReadsCommand = new PrintReadsUnmapped.Builder(java, gatkPrintReadsXmx + "g", tmpDir, gatk, gatkKey, dataDir)
                .setReferenceSequence(refFasta)
                .addInputFiles(inputBams)
                .addInterval("unmapped")
                .setExtraParameters(printReadsParams)
                .build();

        Job printReadsJob = getWorkflow().createBashJob("GATKPrintReadsUnmapped")
                .setMaxMemory(Integer.toString((gatkPrintReadsXmx + gatkOverhead) * 1024))
                .setQueue(queue);

        printReadsJob.addParent(parent);
        printReadsJob.getCommand().setArguments(printReadsCommand.getCommand());
        return Pair.of(printReadsCommand.getOutputFile(), printReadsJob);
    }

    /**
     * Note that file path should not have extension
     *
     * @param inputFile
     *
     * @return
     */
    protected Job indexBamJob(String inputFile) {
        Job jobIndex = this.getWorkflow().createBashJob("index_bam");
        jobIndex.setCommand(this.java + " -Xmx3G -jar "
                + this.picard_dir + "BuildBamIndex.jar"
                + " I=" + inputFile + "bam"
                + " O=" + inputFile + "bai");
        jobIndex.setMaxMemory("5000");
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

}
