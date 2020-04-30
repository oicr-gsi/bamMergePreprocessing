package ca.on.oicr.pde.deciders;

import ca.on.oicr.pde.deciders.GroupableFileFactory.GroupableFile;
import com.google.common.base.Function;
import com.google.common.base.Joiner;
import com.google.common.collect.Iterables;
import com.google.common.collect.ListMultimap;
import com.google.common.collect.MultimapBuilder;
import com.google.common.collect.Multimaps;
import com.google.common.collect.SetMultimap;
import com.google.common.collect.Sets;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import net.sourceforge.seqware.common.hibernate.FindAllTheFiles.Header;
import net.sourceforge.seqware.common.module.ReturnValue;
import org.apache.logging.log4j.LogManager;

/**
 * This decider has several modifications to allow it to group files properly
 * and select the appropriate files for processing.
 * <p>
 * <ol> <li> separateFiles: captures information about the file necessary for
 * grouping and for filtering. It also starts filtering the files, choosing only
 * the most recent file from each sequencer run/lane/barcode/metatype for
 * subsequent processing</li> <li>handleGroupByAttribute: Makes file groups from
 * the same donor LIBRARY_TEMPLATE_TYPE and GROUP_ID.</li> <li>checkFileDetails:
 * restricts processing to only those files that have the correct tissue type
 * and library template type</li> <li>customizeRun: creates the INI file using
 * those files that passed the previous steps. </ol>
 *
 *
 * @author mtaschuk
 */
public class BamMPDecider extends MergingDecider {

    private String ltt = "";
    private List<String> tissueTypes = null;
    private Boolean doFilter = true, doDedup = true, doRemoveDups = true, doSplitNTrim = false;
    private String queue = null;

    private String chrSizes = null;
    private Integer intervalPadding = null;
    private double standCallConf = 30;
    private double standEmitConf = 1;
    private String downsamplingType = null;
    private String dbSNPfile = null;
    private String currentTemplate;
    private Boolean doBQSR = true;

    private Integer filterFlag = 260;
    private Integer minMapQuality = null;
    private final static String[] GATK_DT = {"NONE", "ALL_READS", "BY_SAMPLE"};
    private final static String BAM_METATYPE = "application/bam";
    private final static String TRANSCRIPTOME_SUFFIX = "Aligned.toTranscriptome.out";

    private final Function<List<String>, String> listOfStringsToStringFunction = new Function<List<String>, String>() {
        @Override
        public String apply(List<String> input) {
            String s = Joiner.on(",").skipNulls().join(input);
            if (s.isEmpty()) {
                return null;
            } else {
                return s;
            }
        }
    };

    public BamMPDecider() {
        super(LogManager.getLogger(BamMPDecider.class));
        files = new HashMap<>();

        defineArgument("library-template-type", "Restrict the processing to samples of a particular template type, e.g. WG, EX, TS", false);
        defineArgument("tissue-type", "Restrict the processing to samples of particular tissue types, e.g. P, R, X, C. Multiple values can be comma-separated. Tissue types are processed individually (all R's together, all C's together)", false);
        defineArgument("use-tissue-prep", "Use tissue prep metadata for grouping", false);
        defineArgument("use-tissue-region", "Use tissue region metadata for grouping", false);
        defineArgument("do-filter", "Whether to filter reads out of the BAM file using the value of sam-filter-flag. Default: true", false);
        defineArgument("sam-filter-flag", "The SAM flag to use to remove reads from the input BAM files. Default: 260 (unmapped reads; not primary alignment)", false);
        defineArgument("min-map-quality", "Remove reads with map quality less than input value, e.g. 30", false);
        defineArgument("do-mark-duplicates", "Whether to mark duplicates in the BAM file. Default: true", false);
        defineArgument("do-remove-duplicates", "Whether to remove duplicates in the BAM file. Default: true", false);
        defineArgument("group-by-aligner", "Flag to enable/disable grouping by aligner. Default: true", false);
        parser.accepts("queue", "Optional: Override the default queue setting (production) setting it to something else").withRequiredArg();
        defineArgument("chr-sizes", "Comma separated list of chromosome intervals used to parallelize indel realigning and variant calling. Default: By chromosome", false);
        defineArgument("interval-padding", "Amount of padding to add to each interval (chr-sizes and interval-file determined by decider) in bp. Default: 100", false);
        defineArgument("stand-emit-conf", "Emission confidence threshold to pass to GATK. Default 1", false);
        defineArgument("stand-call-conf", "Calling confidence threshold to pass to GATK. Default 30.", false);
        defineArgument("downsampling", "Set whether or not the variant caller should downsample the reads. Default: false for TS, true for everything else", false);
        defineArgument("dbsnp", "Specify the absolute path to the dbSNP vcf.", false);
        defineArgument("disable-bqsr", "Disable BQSR (BaseRecalibrator + PrintReads steps) and pass indel realigned BAMs directly to variant calling.", false);
        parser.accepts("do-split-and-trim", "Optional: Whether to run splitNcigar reads on the BAM file. Default: false").withRequiredArg();
    }

    @Override
    public ReturnValue init() {
        setMetaType(Arrays.asList("application/bam"));
        ltt = getArgument("library-template-type");
        String tt = getArgument("tissue-type");
        if (!tt.isEmpty()) {
            tissueTypes = Arrays.asList(tt.split(","));
        }
        if (options.has("do-filter")) {
            doFilter = Boolean.valueOf(getArgument("do-filter"));
        }

        if (options.has("do-remove-duplicates")) {
            doRemoveDups = Boolean.valueOf(getArgument("do-remove-duplicates"));
        }

        if (this.options.has("queue")) {
            this.queue = options.valueOf("queue").toString();
        }

        if (!doRemoveDups && options.has("do-mark-duplicates")) {
            doDedup = Boolean.valueOf(getArgument("do-mark-duplicates"));
        }

        if (options.has("sam-filter-flag")) {
            filterFlag = Integer.valueOf(getArgument("sam-filter-flag"));
        }
        if (options.has("min-map-quality")) {
            minMapQuality = Integer.parseInt(getArgument("min-map-quality"));
        }
        if (options.has("group-by-aligner")) {
            groupableFileFactory.setGroupByAligner(Boolean.valueOf(getArgument("group-by-aligner")));
        }

        // GP-597 ==== Tissue prep and Tissue region flags
        if (options.has("use-tissue-prep")) {
            groupableFileFactory.setGroupByTissuePrep(Boolean.valueOf(getArgument("use-tissue-prep")));
        }

        if (options.has("use-tissue-region")) {
            groupableFileFactory.setGroupByTissueRegion(Boolean.valueOf(getArgument("use-tissue-region")));
        }
        // GP-597 ends

        if (options.has("do-remove-duplicates")) {
            doRemoveDups = Boolean.valueOf(getArgument("do-remove-duplicates"));
        }

        if (options.has("chr-sizes")) {
            this.chrSizes = getArgument("chr-sizes");
        }

        if (options.has("interval-padding")) {
            this.intervalPadding = Integer.valueOf(getArgument("interval-padding"));
        } else {
            this.intervalPadding = 100;
        }

        if (options.has("stand-emit-conf")) {
            this.standEmitConf = Double.valueOf(getArgument("stand-emit-conf"));
        } else {
            this.standEmitConf = 1;
        }

        if (options.has("stand-call-conf")) {
            this.standCallConf = Double.valueOf(getArgument("stand-call-conf"));
        } else {
            this.standCallConf = 30;
        }

        if (options.has("disable-bqsr")) {
            this.doBQSR = !Boolean.parseBoolean(getArgument("disable-bqsr"));
        }

        if (options.has("dbsnp")) {
            this.dbSNPfile = getArgument("dbsnp");
        }

        if (options.has("downsampling")) {
            if (getArgument("downsampling").equalsIgnoreCase("false")) {
                this.downsamplingType = "NONE";
            } else if (getArgument("downsampling").equalsIgnoreCase("true")) {
                //do nothing, downsampling is performed by default
            } else {
                throw new RuntimeException("--downsampling parameter expects true/false.");
            }
        }

        if (options.has("do-split-and-trim")) {
            doSplitNTrim = Boolean.valueOf(getArgument("do-split-and-trim"));
        }

        return super.init();
    }

    @Override
    public boolean checkFilePassesFilterBeforeGrouping(FileAttributes fileAttributes) {
        boolean metatypeOK = false;
        boolean bamtypeOK = false;

        try {
            if (fileAttributes.getMetatype().equals(BAM_METATYPE)) {
                metatypeOK = true;
            }
            if (!fileAttributes.getPath().contains(TRANSCRIPTOME_SUFFIX)) {
                bamtypeOK = true;
            }
        } catch (Exception e) {
            log.error("Error checking a file");
            return false;
        }

        return metatypeOK && bamtypeOK;
    }

    @Override
    protected boolean checkFilePassesFilterAfterGrouping(FileAttributes attributes) {
        if (attributes.basename().contains(TRANSCRIPTOME_SUFFIX)) {
            log.debug("Found Transcriptome-aligned reads");
            return false;
        }

        String currentTemplateType = attributes.getOtherAttribute(Header.SAMPLE_TAG_PREFIX.getTitle() + "geo_library_source_template_type");
        if (tissueTypes == null || tissueTypes.contains(attributes.getLimsValue(Lims.TISSUE_TYPE))) {
            if (!ltt.isEmpty()) {
                if (!attributes.getLimsValue(Lims.LIBRARY_TEMPLATE_TYPE).equals(ltt)) {
                    return false;
                }
            }
        } else {
            return false;
        }

        this.currentTemplate = currentTemplateType;
        return true;
    }

    @Override
    public ReturnValue customizeWorkflowRun(WorkflowRun run) {
        return new ReturnValue();
    }

    @Override
    public ReturnValue customizeWorkflowRunAfterCreatingSwidsToLinkWorkflowRunTo(WorkflowRun run) {
        ReturnValue rv = new ReturnValue();

        //the following data structures needs keys sorted - use a tree map to sort keys
        ListMultimap<String, String> inputFilePathsByGroup = MultimapBuilder.treeKeys().linkedListValues().build();
        ListMultimap<String, String> iusLimsKeysByGroup = MultimapBuilder.treeKeys().linkedListValues().build();
        SetMultimap<String, String> outputIdentifierByGroup = MultimapBuilder.treeKeys().hashSetValues().build();

        SetMultimap<String, String> alignerByGroup = MultimapBuilder.treeKeys().hashSetValues().build();

        for (FileAttributes inputFile : run.getFiles()) {
            GroupableFile file = fileSwaToFile.get(inputFile.getOtherAttribute(Header.FILE_SWA.getTitle()));

            String fileGroup = file.getGroupByAttribute();

            inputFilePathsByGroup.put(fileGroup, inputFile.getPath());
            iusLimsKeysByGroup.put(fileGroup, run.getInputIusToOutputIus().get(inputFile.getOtherAttribute(Header.IUS_SWA.getTitle())));

            String idRaw = getCombinedFileName(new FileAttributes[]{inputFile});
            String[] tokens = idRaw.split("_");
            if (tokens.length > 5 && tokens[tokens.length - 1].contains("ius")) {
                idRaw = idRaw.substring(0, idRaw.lastIndexOf("_"));
            }
            outputIdentifierByGroup.put(fileGroup, idRaw);

            alignerByGroup.put(fileGroup, file.getParentWf());
        }

        String alignerName = Iterables.getFirst(Sets.newHashSet(alignerByGroup.values()), "");
        run.addProperty("aligner_name", alignerName);

        //Build the output_identifiers, input_files, and output_ius_lims_keys ini properties for each "output group" (GroupableFileFactory).
        //Each group is separated by the ";" delimter. The ordering of how the ini property is constructed is important, as these ini properties 
        //are split in the workflow and associated together - this is why keys of the following multimaps are ordered by group name.
        //The following example is how these following ini properties should be constructed:
        //input_files=group1_input1,group1_input2;group2_input1
        //output_identifiers=group1_output_name;group2_ouput_name
        //output_ius_lims_keys=group1_input1_ius,group1_input2_ius;group2_input1_ius
        String inputFilePathsStringByGroup = Joiner.on(";").join(Iterables.transform(Multimaps.asMap(inputFilePathsByGroup).values(), listOfStringsToStringFunction));
        run.addProperty("input_files", inputFilePathsStringByGroup);
        String outputIdentifierStringByGroup = Joiner.on(";").join(outputIdentifierByGroup.values());
        run.addProperty("output_identifiers", outputIdentifierStringByGroup);

        // add legacy "identifier" ini property value for BMPP 1.0.3
        run.addProperty("identifier", outputIdentifierStringByGroup);

        String iusLimsKeysStringByGroup = Joiner.on(";").skipNulls().join(Iterables.transform(Multimaps.asMap(iusLimsKeysByGroup).values(), listOfStringsToStringFunction));
        if (iusLimsKeysByGroup.keySet().size() > 1 && !iusLimsKeysStringByGroup.isEmpty()) {
            //when there is only one group, output files should use the IUS-LimsKeys associated with the workflow run
            //where there is more than one group, the workflow needs to assoicate the output files with the approriate IUS-LimsKeys
            //when in test/dry-run mode, the iusLimsKeysStringByGroup will be empty
            run.addProperty("output_ius_lims_keys", iusLimsKeysStringByGroup);
        }

        run.addProperty("do_mark_duplicates", doDedup.toString());
        run.addProperty("do_remove_duplicates", doRemoveDups.toString());
        run.addProperty("do_split_trim_reassign_quality", doSplitNTrim.toString());
        run.addProperty("do_sam_filter", doFilter.toString());
        run.addProperty("samtools_filter_flag", filterFlag.toString());
        run.addProperty("stand-emit-conf", String.valueOf(this.standEmitConf));
        run.addProperty("stand-call-conf", String.valueOf(this.standCallConf));
        run.addProperty("do_bqsr", String.valueOf(this.doBQSR));
        run.addProperty("interval-padding", String.valueOf(this.intervalPadding));

        if (this.chrSizes != null && !this.chrSizes.isEmpty()) {
            run.addProperty("chr_sizes", this.chrSizes);
        }

        if (this.dbSNPfile != null && !dbSNPfile.isEmpty()) {
            run.addProperty("gatk_dbsnp_vcf", this.dbSNPfile);
        }

        if (minMapQuality != null) {
            run.addProperty("samtools_min_map_quality", minMapQuality.toString());
        }

        if (this.currentTemplate != null && this.currentTemplate.equals("TS")) {
            this.downsamplingType = GATK_DT[0];
        }
        if (downsamplingType != null && !downsamplingType.isEmpty()) {
            run.addProperty("downsampling_type", downsamplingType);
        }

        if (this.queue != null) {
            run.addProperty("queue", this.queue);
        } else {
            run.addProperty("queue", " ");
        }

        return rv;
    }

    public static void main(String args[]) {
        List<String> params = new ArrayList<>();
        params.add("--plugin");
        params.add(BamMPDecider.class.getCanonicalName());
        params.add("--");
        params.addAll(Arrays.asList(args));
        System.out.println("Parameters: " + Arrays.deepToString(params.toArray()));
        net.sourceforge.seqware.pipeline.runner.PluginRunner.main(params.toArray(new String[params.size()]));

    }
}
