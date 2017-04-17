package ca.on.oicr.pde.deciders;

import java.text.ParseException;
import java.text.SimpleDateFormat;
import java.util.Date;
import net.sourceforge.seqware.common.hibernate.FindAllTheFiles;
import net.sourceforge.seqware.common.module.ReturnValue;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

public class GroupableFileFactory {

    private final SimpleDateFormat format = new SimpleDateFormat("yyyy-MM-dd HH:mm:ss.S");
    private final Logger log = LogManager.getLogger(GroupableFileFactory.class);

    private boolean groupByAligner = true;
    private boolean groupByTissuePrep = true;
    private boolean groupByTissueRegion = true;

    public GroupableFileFactory() {
    }

    public void setGroupByAligner(boolean groupByAligner) {
        this.groupByAligner = groupByAligner;
    }

    public void setGroupByTissuePrep(boolean groupByTissuePrep) {
        this.groupByTissuePrep = groupByTissuePrep;
    }

    public void setGroupByTissueRegion(boolean groupByTissueRegion) {
        this.groupByTissueRegion = groupByTissueRegion;
    }

    public GroupableFile getGroupableFile(ReturnValue rv) {
        return new GroupableFile(rv);
    }

    public class GroupableFile {

        private Date date = null;
        private String iusDetails = null;
        private String parentWf = "";
        private String groupByAttribute = null;
        private String path = null;

        private GroupableFile(ReturnValue rv) {
            FileAttributes fa = new FileAttributes(rv, rv.getFiles().get(0));

            try {
                date = format.parse(rv.getAttribute(FindAllTheFiles.Header.PROCESSING_DATE.getTitle()));
            } catch (ParseException ex) {
                log.error("Bad date!", ex);
            }

            iusDetails = fa.getLibrarySample() + fa.getSequencerRun() + fa.getLane() + fa.getBarcode();

            groupByAttribute = fa.getDonor() + ":" + fa.getLimsValue(Lims.TISSUE_ORIGIN) + ":" + fa.getLimsValue(Lims.LIBRARY_TEMPLATE_TYPE);

            String wfName = rv.getAttribute(FindAllTheFiles.Header.WORKFLOW_NAME.getTitle());
            this.parentWf = wfName;
            if (groupByAligner && null != wfName && !wfName.isEmpty()) {
                //Grouping by workflow name (we don't care about version)

                groupByAttribute = groupByAttribute.concat(":" + this.parentWf);
            }

            if (null != fa.getLimsValue(Lims.TISSUE_TYPE)) {
                groupByAttribute = groupByAttribute.concat(":" + fa.getLimsValue(Lims.TISSUE_TYPE));
            }

            if (groupByTissuePrep && null != fa.getLimsValue(Lims.TISSUE_PREP)) {
                groupByAttribute = groupByAttribute.concat(":" + fa.getLimsValue(Lims.TISSUE_PREP));
            }

            if (groupByTissueRegion && null != fa.getLimsValue(Lims.TISSUE_REGION)) {
                groupByAttribute = groupByAttribute.concat(":" + fa.getLimsValue(Lims.TISSUE_REGION));
            }

            if (null != fa.getLimsValue(Lims.GROUP_ID)) {
                groupByAttribute = groupByAttribute.concat(":" + fa.getLimsValue(Lims.GROUP_ID));
            }

            if (null != fa.getLimsValue(Lims.TARGETED_RESEQUENCING)) {
                groupByAttribute = groupByAttribute.concat(":" + fa.getLimsValue(Lims.TARGETED_RESEQUENCING));
            }

            path = rv.getFiles().get(0).getFilePath() + "";
        }

        public Date getDate() {
            return date;
        }

        public String getGroupByAttribute() {
            return groupByAttribute;
        }

        public String getIusDetails() {
            return iusDetails;
        }

        public String getPath() {
            return path;
        }

        public String getParentWf() {
            return parentWf;
        }
    }

}
