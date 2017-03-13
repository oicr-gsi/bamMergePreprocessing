package ca.on.oicr.pde.deciders;

import java.text.ParseException;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Date;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import net.sourceforge.seqware.common.hibernate.FindAllTheFiles.Header;
import net.sourceforge.seqware.common.module.ReturnValue;
import net.sourceforge.seqware.common.util.Log;

/**
 * This decider has several modifications to allow it to group files properly
 * and select the appropriate files for processing.
 *
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
public class BamMPDecider extends OicrDecider {

    private Map<String, BeSmall> fileSwaToSmall;
    private String ltt = "";
    private List<String> tissueTypes = null;
    private SimpleDateFormat format = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss.S");
    private Boolean doFilter = true, doDedup = true, doRemoveDups = true;
    private Boolean useTPrep = true;
    private Boolean useTRegion = true;
    private Boolean groupByAligner = true;
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
    private String outputPrefix;
    private String outputDir;

    public BamMPDecider() {
        super();
        files = new HashMap<String, FileAttributes>();
        fileSwaToSmall = new HashMap<String, BeSmall>();
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
        defineArgument("output-path", "Optional: the path where the files should be copied to "
                + "after analysis. Corresponds to output-prefix in INI file. Default: ./", false);
        defineArgument("output-folder", "Optional: the name of the folder to put the output into relative to "
                + "the output-path. Corresponds to output-dir in INI file. Default: seqware-results", false);
        parser.accepts("queue", "Optional: Override the default queue setting (production) setting it to something else").withRequiredArg();
        defineArgument("chr-sizes", "Comma separated list of chromosome intervals used to parallelize indel realigning and variant calling. Default: By chromosome", false);
        defineArgument("interval-padding", "Amount of padding to add to each interval (chr-sizes and interval-file determined by decider) in bp. Default: 100", false);
        defineArgument("stand-emit-conf", "Emission confidence threshold to pass to GATK. Default 1", false);
        defineArgument("stand-call-conf", "Calling confidence threshold to pass to GATK. Default 30.", false);
        defineArgument("downsampling", "Set whether or not the variant caller should downsample the reads. Default: false for TS, true for everything else", false);
        defineArgument("dbsnp", "Specify the absolute path to the dbSNP vcf.", false);
        defineArgument("disable-bqsr", "Disable BQSR (BaseRecalibrator + PrintReads steps) and pass indel realigned BAMs directly to variant calling.", false);
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
            groupByAligner = Boolean.valueOf(getArgument("group-by-aligner"));
        }

        // GP-597 ==== Tissue prep and Tissue region flags
        if (options.has("use-tissue-prep")) {
            this.useTPrep = Boolean.valueOf(getArgument("use-tissue-prep"));
        }

        if (options.has("use-tissue-region")) {
            this.useTRegion = Boolean.valueOf(getArgument("use-tissue-region"));
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

        if (this.options.has("output-path")) {
            this.outputPrefix = options.valueOf("output-path").toString();
            if (!this.outputPrefix.endsWith("/")) {
                this.outputPrefix += "/";
            }
        } else { this.outputPrefix = "./";}

        if (this.options.has("output-folder")) {
            this.outputDir = options.valueOf("output-folder").toString();
        } else { this.outputDir = "seqware-results"; }
        
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

        return super.init();
    }

    /**
     * Final check
     *
     * @param commaSeparatedFilePaths
     * @param commaSeparatedParentAccessions
     * @return
     */
    @Override
    protected ReturnValue doFinalCheck(String commaSeparatedFilePaths, String commaSeparatedParentAccessions) {

        if (this.currentTemplate != null && this.currentTemplate.equals("TS")) {
            this.downsamplingType = GATK_DT[0];
        }
        
        boolean haveRefAlignedBam = false; 
        String[] filePaths = commaSeparatedFilePaths.split(",");
        for (String p : filePaths) {
            for (BeSmall bs : fileSwaToSmall.values()) {
                if (!bs.getPath().equals(p)) {
                    continue;
                }
                if (!haveRefAlignedBam) {
                    haveRefAlignedBam = !p.contains(TRANSCRIPTOME_SUFFIX);}                   
                }

        }
        
        if (!haveRefAlignedBam) {
            Log.error("The Decider was not able to find Reference-aligned reads (bam file), WON'T RUN");
            return new ReturnValue(ReturnValue.INVALIDPARAMETERS);
        }


        return super.doFinalCheck(commaSeparatedFilePaths, commaSeparatedParentAccessions);
    }

    @Override
    protected boolean checkFileDetails(FileAttributes attributes) {
        boolean rv = super.checkFileDetails(attributes);

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
        return rv;
    }

    /**
     * This method is extended in the GATK decider so that only the most recent
     * file for each sequencer run, lane, barcode and filetype is kept.
     */
    @Override
    public Map<String, List<ReturnValue>> separateFiles(List<ReturnValue> vals, String groupBy) {
        Map<String, ReturnValue> iusDeetsToRV = new HashMap<String, ReturnValue>();

        //Iterate through the potential files
        for (ReturnValue currentRV : vals) {
            //set aside information needed for subsequent processing
            boolean metatypeOK = false;

            for (int f = 0; f < currentRV.getFiles().size(); f++) {
                try {
                    if (currentRV.getFiles().get(f).getMetaType().equals(BAM_METATYPE)) {
                        metatypeOK = true;
                    }
                } catch (Exception e) {
                    Log.stderr("Error checking a file");
                    continue;
                }
            }
            if (!metatypeOK) {
                continue; // Go to the next value
            }
            BeSmall currentSmall = new BeSmall(currentRV);
            fileSwaToSmall.put(currentRV.getAttribute(groupBy), currentSmall);

            //make sure you only have the most recent single file for each
            //sequencer run + lane + barcode + meta-type
            String fileDeets = currentSmall.getIusDetails();
            Date currentDate = currentSmall.getDate();

            //if there is no entry yet, add it
            if (iusDeetsToRV.get(fileDeets) == null) {
                Log.debug("Adding file " + fileDeets + " -> \n\t" + currentSmall.getPath());
                iusDeetsToRV.put(fileDeets, currentRV);
            } //if there is an entry, compare the current value to the 'old' one in
            //the map. if the current date is newer than the 'old' date, replace
            //it in the map
            else {
                ReturnValue oldRV = iusDeetsToRV.get(fileDeets);
                BeSmall oldSmall = fileSwaToSmall.get(oldRV.getAttribute(Header.FILE_SWA.getTitle()));
                Date oldDate = oldSmall.getDate();
                if (currentDate.after(oldDate)) {
                    Log.debug("Adding file " + fileDeets + " -> \n\t" + currentSmall.getDate()
                            + "\n\t instead of file "
                            + "\n\t" + oldSmall.getDate());
                    iusDeetsToRV.put(fileDeets, currentRV);
                } else {
                    Log.debug("Disregarding file " + fileDeets + " -> \n\t" + currentSmall.getDate()
                            + "\n\tas older than duplicate sequencer run/lane/barcode in favour of "
                            + "\n\t" + oldSmall.getDate());
                    Log.debug(currentDate + " is before " + oldDate);
                }
            }
        }
        //only use those files that entered into the iusDeetsToRV
        //since it's a map, only the most recent values
        List<ReturnValue> newValues = new ArrayList<ReturnValue>(iusDeetsToRV.values());
        Map<String, List<ReturnValue>> map = new HashMap<String, List<ReturnValue>>();
        //group files according to the designated header (e.g. sample SWID)
        for (ReturnValue r : newValues) {
            String currVal = fileSwaToSmall.get(r.getAttribute(Header.FILE_SWA.getTitle())).getGroupByAttribute();

            List<ReturnValue> vs = map.get(currVal);
            if (vs == null) {
                vs = new ArrayList<ReturnValue>();
            }
            vs.add(r);
            map.put(currVal, vs);
        }

        return map;
    }

    @Override
    public ReturnValue customizeRun(WorkflowRun run) {
        String inputFiles = "";
        for (FileAttributes atts : run.getFiles()) {
            if (!inputFiles.isEmpty()) {
                inputFiles += ",";
            }
            String filePath = atts.getPath();
            if (!filePath.contains(TRANSCRIPTOME_SUFFIX))
                inputFiles += filePath;
        }

        //Use aligner name, if available
        String[] filePaths = inputFiles.split(",");
        if (filePaths.length > 0) {
            for (BeSmall bs : fileSwaToSmall.values()) {
                if (!bs.getPath().equals(filePaths[0])) {
                    continue;
                }
                String alignerName = bs.getParentWf();
                if (!alignerName.isEmpty() && this.groupByAligner) {
                    run.addProperty("aligner_name", alignerName);
                }
                break;
            }
        }

        run.addProperty("input_files", inputFiles);

        // Get read of _ius addition when dealing with single file
        String idRaw = getCombinedFileName(run.getFiles());
        String[] tokens = idRaw.split("_");
        if (tokens.length > 5 && tokens[tokens.length - 1].contains("ius")) {
            idRaw = idRaw.substring(0, idRaw.lastIndexOf("_"));
        }

        run.addProperty("identifier", idRaw);
        run.addProperty("do_mark_duplicates", doDedup.toString());
        run.addProperty("do_remove_duplicates", doRemoveDups.toString());
        run.addProperty("do_sam_filter", doFilter.toString());
        run.addProperty("samtools_filter_flag", filterFlag.toString());
        run.addProperty("stand-emit-conf", String.valueOf(this.standEmitConf));
        run.addProperty("stand-call-conf", String.valueOf(this.standCallConf));
        run.addProperty("do_bqsr", String.valueOf(this.doBQSR));
        run.addProperty("interval-padding", String.valueOf(this.intervalPadding));

        if (this.chrSizes != null && !this.chrSizes.isEmpty()) {
            run.addProperty("chr-sizes", this.chrSizes);
        }     

        if (this.dbSNPfile != null && !dbSNPfile.isEmpty()) {
            run.addProperty("gatk_dbsnp_vcf", this.dbSNPfile);
        }

        if (minMapQuality != null) {
            run.addProperty("samtools_min_map_quality", minMapQuality.toString());
        }
        
        if (downsamplingType != null && !downsamplingType.isEmpty()) {
            run.addProperty("downsampling_type", downsamplingType);
        }
        
        if (this.queue != null) {
            run.addProperty("queue", this.queue);
        } else {
            run.addProperty("queue", " ");
        }

        return new ReturnValue();
    }

    private class BeSmall {

        private Date date = null;
        private String iusDetails = null;
        private String parentWf = "";
        private String groupByAttribute = null;
        private String path = null;

        public BeSmall(ReturnValue rv) {
            try {
                date = format.parse(rv.getAttribute(Header.PROCESSING_DATE.getTitle()));
            } catch (ParseException ex) {
                Log.error("Bad date!", ex);
                ex.printStackTrace();
            }
            FileAttributes fa = new FileAttributes(rv, rv.getFiles().get(0));
            iusDetails = fa.getLibrarySample() + fa.getSequencerRun() + fa.getLane() + fa.getBarcode();
            String wfName = rv.getAttribute(Header.WORKFLOW_NAME.getTitle());
            groupByAttribute = fa.getDonor() + ":" + fa.getLimsValue(Lims.TISSUE_ORIGIN) + ":" + fa.getLimsValue(Lims.LIBRARY_TEMPLATE_TYPE);

            if (null != wfName && !wfName.isEmpty() && groupByAligner) {
                this.parentWf = wfName;
                groupByAttribute.concat(":" + this.parentWf);
            }

            if (null != fa.getLimsValue(Lims.TISSUE_TYPE)) {
                groupByAttribute = groupByAttribute.concat(":" + fa.getLimsValue(Lims.TISSUE_TYPE));
            }

            if (null != fa.getLimsValue(Lims.TISSUE_PREP) && useTPrep) {
                groupByAttribute = groupByAttribute.concat(":" + fa.getLimsValue(Lims.TISSUE_PREP));
            }

            if (null != fa.getLimsValue(Lims.TISSUE_REGION) && useTRegion) {
                groupByAttribute = groupByAttribute.concat(":" + fa.getLimsValue(Lims.TISSUE_REGION));
            }

            if (null != fa.getLimsValue(Lims.GROUP_ID)) {
                groupByAttribute = groupByAttribute.concat(":" + fa.getLimsValue(Lims.GROUP_ID));
            }

            if (null != fa.getLimsValue(Lims.TARGETED_RESEQUENCING)) {
                groupByAttribute = groupByAttribute.concat(":" + fa.getLimsValue(Lims.TARGETED_RESEQUENCING));
            }

            //Grouping by workflow name (we don't care about version)
            path = rv.getFiles().get(0).getFilePath() + "";
        }

        public Date getDate() {
            return date;
        }

        public void setDate(Date date) {
            this.date = date;
        }

        public String getGroupByAttribute() {
            return groupByAttribute;
        }

        public void setGroupByAttribute(String groupByAttribute) {
            this.groupByAttribute = groupByAttribute;
        }

        public String getIusDetails() {
            return iusDetails;
        }

        public void setIusDetails(String iusDetails) {
            this.iusDetails = iusDetails;
        }

        public String getPath() {
            return path;
        }

        public void setPath(String path) {
            this.path = path;
        }

        public String getParentWf() {
            return parentWf;
        }
    }

    public static void main(String args[]) {
        List<String> params = new ArrayList<String>();
        params.add("--plugin");
        params.add(BamMPDecider.class.getCanonicalName());
        params.add("--");
        params.addAll(Arrays.asList(args));
        System.out.println("Parameters: " + Arrays.deepToString(params.toArray()));
        net.sourceforge.seqware.pipeline.runner.PluginRunner.main(params.toArray(new String[params.size()]));

    }
}
