package org.broadinstitute.cga.tools.gatk.walkers.cancer.mutect;

import org.broadinstitute.sting.utils.BaseUtils;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;

/**
 * Created by IntelliJ IDEA.
 * User: kcibul
 * Date: 3/29/12
 * Time: 4:47 PM
 * To change this template use File | Settings | File Templates.
 */
public class RecalibratedLocalQualityScores {
    protected int[] positiveStrandBaseCounts; 
    protected byte[] positiveStrandQualities;

    protected int[] negativeStrandBaseCounts;
    protected byte[] negativeStrandQualities;
    byte ref;

    public RecalibratedLocalQualityScores(byte ref, ReadBackedPileup pileup) {
        this.ref = ref;
        this.positiveStrandBaseCounts =  pileup.getPositiveStrandPileup().getBaseCounts();
        this.negativeStrandBaseCounts =  pileup.getNegativeStrandPileup().getBaseCounts();
        
        this.positiveStrandQualities = calculateEmpiricalQualities(ref, positiveStrandBaseCounts);
        this.negativeStrandQualities = calculateEmpiricalQualities(ref, negativeStrandBaseCounts);
    }

    static protected byte[] calculateEmpiricalQualities(byte ref, int[] counts) {
        byte[] newQualities = new byte[BaseUtils.BASES.length];
        int total = 0;
        for(int c : counts) { total += c; }

        for( byte base : BaseUtils.BASES) {
            int baseIndex = BaseUtils.simpleBaseToBaseIndex(base);
            int baseCount = (base == ref)?0:counts[baseIndex];

            // calculate a phred-like empirical quality score
            double qscore = -10 * Math.log10(((double)baseCount + 1d) / ((double) total + 4d));

            // TODO: should we carry this forward and allow fractional phred scores?
            newQualities[baseIndex] = (byte) Math.round(qscore);
        }

        return newQualities;
    }

    /**
     * Get the local empirical quality score for the supplied read
     * @param pe read
     * @return new phred-scale quality score
     */
    public byte getQual(PileupElement pe) {
        int baseIndex = BaseUtils.simpleBaseToBaseIndex(pe.getBase());
        return (pe.getRead().getReadNegativeStrandFlag())?negativeStrandQualities[baseIndex]:positiveStrandQualities[baseIndex];
    }

    @Override
    public String toString() {
        return "RecalibratedLocalQualityScores{" +
                "ref= " + (char) ref +
                " +counts=(" + arrayToString(positiveStrandBaseCounts) + ") " +
                " -counts=(" + arrayToString(negativeStrandBaseCounts) + ") " +
                " +quals=(" + arrayToString(positiveStrandQualities) + ") " +
                " -quals=(" + arrayToString(negativeStrandQualities) + ") " +
                '}';
    }

    private StringBuilder arrayToString(int[] data) {
        StringBuilder sb = new StringBuilder();
        if (data.length == 0) { return sb; }
        sb.append(data[0]);
        for(int i=1; i<data.length; i++) { sb.append(",").append(data[i]); }
        return sb;
    }

    private StringBuilder arrayToString(byte[] data) {
        StringBuilder sb = new StringBuilder();
        if (data.length == 0) { return sb; }
        sb.append(data[0]);
        for(int i=1; i<data.length; i++) { sb.append(",").append(data[i]); }
        return sb;
    }

}
