package ca.on.oicr.pde.tools.gatk3;

import ca.on.oicr.pde.tools.common.AbstractCommand;
import com.google.common.base.Joiner;
import java.util.Collection;
import java.util.LinkedList;
import java.util.List;
import org.apache.commons.io.FilenameUtils;
import org.apache.commons.lang3.RandomStringUtils;

/**
 *
 * @author mlaszloffy
 */
public class PrintReadsUnmapped extends AbstractCommand {

    private String outputFile;

    private PrintReadsUnmapped() {
    }

    public String getOutputFile() {
        return outputFile;
    }

    public static class Builder extends AbstractGatkBuilder<Builder> {

        private String covariatesTablesFile;
        private Integer preserveQscoresLessThan;
        List<String> inputFiles = new LinkedList<>();

        public Builder(String javaPath, String maxHeapSize, String tmpDir, String gatkJarPath, String gatkKey, String outputDir) {
            super(javaPath, maxHeapSize, tmpDir, gatkJarPath, gatkKey, outputDir);
        }

        public Builder addInputFile(String inputFile) {
            inputFiles.add(inputFile);
            return this;
        }

        public Builder addInputFiles(Collection<String> inputFiles) {
            this.inputFiles.addAll(inputFiles);
            return this;
        }

        public PrintReadsUnmapped build() {
            String intervalString;
            if (!intervals.isEmpty()) {
                intervalString = Joiner.on("_").join(intervals).replace(":", "-");
            } else {
                intervalString = RandomStringUtils.randomAlphanumeric(4);
            }

            //GP-604: When merging, use generic name
            String outputFilePath = "unmapped.bam";
            if (inputFiles.size() == 1) {
                outputFilePath = FilenameUtils.getBaseName(inputFiles.get(0)) + "_" + intervalString + ".bam";
            }

            List<String> c = build("PrintReads");

            if (covariatesTablesFile != null) {
                c.add("--BQSR");
                c.add(covariatesTablesFile);
            }

            //GP-604: Support multiple inputs
            for (String inFile : inputFiles) {
                c.add("-I");
                c.add(inFile);
            }

            if (preserveQscoresLessThan != null) {
                c.add("--preserve_qscores_less_than");
                c.add(preserveQscoresLessThan.toString());
            }

            c.add("--out");
            c.add(outputFilePath);

            PrintReadsUnmapped cmd = new PrintReadsUnmapped();
            cmd.command.addAll(c);
            cmd.outputFile = outputFilePath;
            return cmd;
        }

    }
}
