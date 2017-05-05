package ca.on.oicr.pde.workflows;

import ca.on.oicr.pde.testing.workflow.DryRun;
import ca.on.oicr.pde.testing.workflow.TestDefinition;
import java.io.File;
import java.io.IOException;
import net.sourceforge.seqware.pipeline.workflowV2.AbstractWorkflowDataModel;
import org.apache.commons.io.FileUtils;
import org.testng.annotations.Test;

/**
 *
 * @author mlaszloffy
 */
public class BamMPWorkflowTest {

    public BamMPWorkflowTest() {
    }

    @Test
    public void validateRegressionTestDefinition() throws IllegalAccessException, InstantiationException, IOException, Exception {
        TestDefinition td = TestDefinition.buildFromJson(FileUtils.readFileToString(new File("src/test/resources/tests.json")));
        for (TestDefinition.Test t : td.getTests()) {
            DryRun d = new DryRun(System.getProperty("bundleDirectory"), t.getParameters(), BamMPWorkflow.class);
            AbstractWorkflowDataModel workflowModel = d.buildWorkflowModel();
            d.validateWorkflow();
        }
    }

}
