package org.broadinstitute.cga.tools.gatk.walkers.cancer.mutect;

import net.sf.samtools.CigarElement;
import net.sf.samtools.CigarOperator;
import net.sf.samtools.SAMRecord;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class SequenceUtils {

    public static String createSequenceContext(ReferenceContext ref, int size) {
        // create a context of 3 bases before, then 'x' then three bases after
        int offset = ref.getLocus().getStart() - ref.getWindow().getStart();
        StringBuilder sb = new StringBuilder(7);

        for(byte b : Arrays.copyOfRange(ref.getBases(), Math.max(0, offset - size), offset)) {
            sb.append(Character.toUpperCase((char)b));
        }
        sb.append('x');
        for(byte b : Arrays.copyOfRange(ref.getBases(), offset + 1, Math.min(ref.getBases().length,offset + 1 + size))) {
            sb.append(Character.toUpperCase((char)b));
        }
        return sb.toString();
    }

    public static int[] getStrandContingencyTable(LocusReadPile refPile, LocusReadPile mutantPile) {
        // Construct a 2x2 contingency table of
        //            pos     neg
        //      REF    a       b
        //      MUT    c       d
        //
        // and return an array of {a,b,c,d}

        int a = 0, b = 0, c = 0, d = 0;
        for (SAMRecord rec : refPile.finalPileupReads) {
            if (rec.getReadNegativeStrandFlag()) { b++;} else { a++; }
        }
        for (SAMRecord rec : mutantPile.finalPileupReads) {
            if (rec.getReadNegativeStrandFlag()) { d++;} else { c++; }
        }

        return new int[]{a,b,c,d};
    }


    public static List<Integer> getForwardOffsetsInRead(ReadBackedPileup p) {
        return getOffsetsInRead(p, true);
    }

    public static List<Integer> getReverseOffsetsInRead(ReadBackedPileup p) {
        return getOffsetsInRead(p, false);
    }

    public static List<Integer> getOffsetsInRead(ReadBackedPileup p, boolean useForwardOffsets) {
        List<Integer> positions = new ArrayList<Integer>();
        for(PileupElement pe : p) {

            positions.add(
                    Math.abs((int)(p.getLocation().getStart() - (useForwardOffsets?pe.getRead().getAlignmentStart():pe.getRead().getAlignmentEnd())))
            );
        }

        return positions;
    }

    // TODO: see if this is cheaper with a GATKSAMRecord...
    public static boolean isReadHeavilySoftClipped(SAMRecord rec, float threshold) {
        int total = 0;
        int clipped = 0;
        for(CigarElement ce : rec.getCigar().getCigarElements()) {
            total += ce.getLength();
            if (ce.getOperator() == CigarOperator.SOFT_CLIP) {
                clipped += ce.getLength();
            }
        }

        return ((float) clipped / (float)total >= threshold);
    }

}
