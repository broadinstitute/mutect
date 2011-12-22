package org.broadinstitute.cga.tools.gatk.walkers.cancer.mutect;

import net.sf.picard.util.FormatUtil;
import org.broadinstitute.sting.gatk.contexts.AlignmentContext;

import java.io.PrintStream;

/**
 * Handles writing of wiggle format coverage data to an underlying PrintStream from GATK
 */
public class CoverageWiggleFileWriter {
    protected PrintStream ps;

    private int lastContigIndex = -1;
    private long lastPosition = -1;

    private final FormatUtil fmt = new FormatUtil();


    public CoverageWiggleFileWriter(PrintStream ps) {
        this(ps, "SomaticCoverage");
    }

    public CoverageWiggleFileWriter(PrintStream ps, String trackName) {
        this.ps = ps;

        safePrint("track type=wiggle_0 name=" + trackName + "\n");
    }

    private boolean hasPrintStream() {
        return (ps != null);
    }

    private void safePrint(String s) {
        if (hasPrintStream()) {
            ps.print(s);
        }
    }

    public void writeCoverage(final AlignmentContext context, final boolean isBaseCovered) {
        writeCoverage(context, isBaseCovered?1:0);
    }

    public void possiblyGenerateHeader(StringBuilder sb, AlignmentContext context) {
        if (this.lastContigIndex != context.getLocation().getContigIndex() ||
            this.lastPosition + 1 != context.getPosition()) {
                this.lastContigIndex = context.getLocation().getContigIndex();
                sb.append("fixedStep").append(" ")
                  .append("chrom=").append(context.getContig()).append(" ")
                  .append("start=").append(context.getPosition()).append(" ")
                  .append("step=1")
                  .append("\n");
        }
        this.lastPosition = context.getPosition();
    }

    public void writeCoverage(final AlignmentContext context, final int depth) {
        if (hasPrintStream()) {
            final StringBuilder sb = new StringBuilder();
            possiblyGenerateHeader(sb,context);

            sb.append(depth).append("\n");

            safePrint(sb.toString());
        }
    }

    public void writeCoverage(final AlignmentContext context, final double power) {
        if (hasPrintStream()) {
            final StringBuilder sb = new StringBuilder();
            possiblyGenerateHeader(sb,context);

            sb.append(fmt.format(power)).append("\n");

            safePrint(sb.toString());
        }
    }

}
