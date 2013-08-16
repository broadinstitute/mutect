package org.broadinstitute.cga.tools.gatk.walkers.cancer.mutect;

import net.sf.picard.reference.IndexedFastaSequenceFile;
import net.sf.samtools.CigarElement;
import net.sf.samtools.CigarOperator;
import net.sf.samtools.SAMFileHeader;
import org.broadinstitute.sting.BaseTest;
import org.broadinstitute.sting.gatk.contexts.ReferenceContext;
import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.GenomeLocParser;
import org.broadinstitute.sting.utils.Utils;
import org.broadinstitute.sting.utils.fasta.CachingIndexedFastaSequenceFile;
import org.broadinstitute.sting.utils.pileup.PileupElement;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileup;
import org.broadinstitute.sting.utils.pileup.ReadBackedPileupImpl;
import org.broadinstitute.sting.utils.sam.ArtificialSAMUtils;
import org.broadinstitute.sting.utils.sam.GATKSAMRecord;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.Test;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.*;

/**
 * Unit tests for MuTect
 */
public class MuTectWalkerUnitTest extends BaseTest {
    public byte ref = 'A';
    public byte alt = 'T';

    // example genome loc parser for this test, can be deleted if you don't use the reference
    private GenomeLocParser genomeLocParser;

    // example fasta index file, can be deleted if you don't use the reference
    private IndexedFastaSequenceFile seq;

    @BeforeClass
    public void setup() throws FileNotFoundException {
        // sequence
        seq = new CachingIndexedFastaSequenceFile(new File(hg19Reference));
        genomeLocParser = new GenomeLocParser(seq);
    }

    @Test
    public void basicPileupTest() {
        final SAMFileHeader header = ArtificialSAMUtils.createArtificialSamHeader(seq.getSequenceDictionary());
        final GenomeLoc myLocation = genomeLocParser.createGenomeLoc("1", 10);

        final int pileupSize = 100;
        final int readLength = 10;
        final List<GATKSAMRecord> reads = new LinkedList<GATKSAMRecord>();
        for ( int i = 0; i < pileupSize; i++ ) {
            final GATKSAMRecord read = ArtificialSAMUtils.createArtificialRead(header, "myRead" + i, 0, 1, readLength);
            final byte[] bases = Utils.dupBytes((byte) 'A', readLength);
            bases[0] = (byte)(i % 2 == 0 ? 'A' : 'C'); // every other base is a C

            // set the read's bases and quals
            read.setReadBases(bases);
            read.setBaseQualities(Utils.dupBytes((byte)30, readLength));
            reads.add(read);
        }

        // create a pileup with all reads having offset 0
        final ReadBackedPileup pileup = new ReadBackedPileupImpl(myLocation, reads, 0);
        // TODO -- add some tests here using pileup

        // this code ensures that the pileup example is correct.  Can be deleted
        Assert.assertEquals(pileup.getNumberOfElements(), pileupSize);
        int nA = 0, nC = 0;
        for ( final PileupElement p : pileup ) {
            if ( p.getBase() == 'A' ) nA++;
            if ( p.getBase() == 'C' ) nC++;
        }
        Assert.assertEquals(nA, pileupSize / 2);
        Assert.assertEquals(nC, pileupSize / 2);
    }

    @Test(enabled = true)
    public void simulate() {
        SAMFileHeader header = ArtificialSAMUtils.createArtificialSamHeader(1, 1, 1000);
        GenomeLocParser genomeLocParser;
        genomeLocParser = new GenomeLocParser(header.getSequenceDictionary());
        GenomeLoc loc = genomeLocParser.createGenomeLoc("chr1", 500);

        final int readLen = 100;
        final double qScoreMean = 35;

        final double qScoreStandardDeviation = 0.000001;

        final int trials = 100; // set to 10,000 for paper

        final Integer[] depths = new Integer[]{0, 5,10,15,20,25,30,35,40,45,50};
        final Double[] trueAlleleFractions = new Double[]{0.1, 0.2, 0.4, 0.5};

        ReferenceContext refContext = new ReferenceContext(genomeLocParser, loc, (byte)ref);
        Random r = new Random();
        System.out.println("f\tdepth\ttrials\tsuccess");
        for (int depth : depths) {
            for (double f : trueAlleleFractions) {
                int success = 0;
                for (int i=0; i< trials;i++) {
                    int alts = getBinomial(depth, f);

                    final byte[] bases = new byte[depth];
                    Arrays.fill(bases,(byte)ref);
                    Arrays.fill(bases,0,alts,(byte)alt);

                    final byte[] qscores = new byte[depth];
                    
                    // instead draw from a gaussian, and pick integer qualities
                    // such that they average to the non-integer quality
                    for(int j=0;j<depth;j++) {

                        double qScoreDouble = getGaussian(r, qScoreMean, qScoreStandardDeviation);
                        byte qScore = (byte) Math.floor(qScoreDouble);

                        double frac = qScore - Math.floor(qScore);
                        if (r.nextDouble() <= frac) { qScore++; }
                        
                        qscores[j] = qScore;
                        
                        // now that we have the quality score, check to see if
                        // we should flip the "true" base
                        double eps = Math.pow(10, -qScore/10);

                        // pick a random number, if it's smaller than eps/3
                        // we made a mistake (to the other base)
                        if (r.nextDouble() < eps/3) {
                            bases[j] = (byte)((bases[j] == ((byte)ref))?alt : ref);
                        }
                    }

                    ReadBackedPileup tumor  = createReadBackedPileup(header, loc, readLen, qscores, bases);
                    LocusReadPile tp = new LocusReadPile(tumor, (char)ref, 0,0,false, false, false);

                    // to use f or f-star?
                    double observedF = tp.estimateAlleleFraction((char)ref, (char)alt);

                    double newLod = tp.calculateAltVsRefLOD((byte) alt, observedF, 0);

                    VariableAllelicRatioGenotypeLikelihoods tumorFStarGl = tp.calculateLikelihoods(observedF, tp.finalPileup);
                    double lod = tp.getHetVsRef(tumorFStarGl, (char)ref, (char)alt);
                    
//                    System.out.println("Using Q="+qScore+" depth="+depth+" alts=" + alts + " trueF = " + f +
//                            " observedF=" + observedF +" and got GATK LOD:" + lod + " and new LOD:" + newLod +
//                            " with lodAlt = " + lodAlt + " and lodRef = " + lodRef);

                    Assert.assertEquals(newLod, lod, 0.0001, "GATK and local LOD scores do not agree");


                    if (lod >= 6.3) { success++; }
                }
        
                System.out.println("" + f + "\t" + depth + "\t" + trials + "\t" + success);
            }
        }
    }

    public static ReadBackedPileup createReadBackedPileup(final SAMFileHeader header,
                                                          final GenomeLoc loc,
                                                          final int readLen,
                                                          final byte[] pileQscores,
                                                          final byte[] pileBases) {
        
        boolean[] negativeStrandFlags = new boolean[pileBases.length];
        Arrays.fill(negativeStrandFlags, false);
        
        return createReadBackedPileup(header,loc,readLen, pileQscores, pileBases, negativeStrandFlags);
    }
    
    public static ReadBackedPileup createReadBackedPileup(final SAMFileHeader header,
                                                                final GenomeLoc loc, 
                                                                final int readLen, 
                                                                final byte[] pileQscores,
                                                                final byte[] pileBases,
                                                                final boolean[] pileNegStrandFlag) {

        final int pos = loc.getStart();
        
        final List<PileupElement> pileupElements = new ArrayList<PileupElement>();
        for ( int i = 0; i < pileBases.length; i++ ) {
            final String readName = "read" + i;
            int leftStart = pos - ((int)readLen/2);
            int offset = pos - leftStart;

            GATKSAMRecord read = ArtificialSAMUtils.createArtificialRead(header, readName, 0, leftStart, readLen);
            read.setMappingQuality(60);
            read.setReadNegativeStrandFlag(pileNegStrandFlag[i]);

            byte[] bases = read.getReadBases();
            bases[offset] = pileBases[i];
            read.setReadBases(bases);

            byte[] quals = read.getBaseQualities();
            quals[offset] = pileQscores[i];
            read.setBaseQualities(quals);

            PileupElement pe =
                    new PileupElement(read, pos - leftStart, new CigarElement(readLen, CigarOperator.M), 0, 0);

            pileupElements.add(pe);
        }

        Collections.sort(pileupElements);
        return new ReadBackedPileupImpl(loc, pileupElements);
    }

    @Test
    public void testPileupComparator() {
        byte A = (byte) 'A';
        byte C = (byte) 'C';
        byte T = (byte) 'T';
        byte Q30 = (byte) 30;
        byte Q60 = (byte) 60;

        MuTect.PileupComparatorByAltRefQual c = new MuTect.PileupComparatorByAltRefQual(A);

        // recall:
        //      A<B -> -1
        //      A=B ->  0
        //      A>B ->  1
        Assert.assertEquals(c.internalCompare(A, Q30, A, Q60), 1);
        Assert.assertEquals(c.internalCompare(A, Q60, A, Q30), -1);
        Assert.assertEquals(c.internalCompare(A, Q60, A, Q60), 0);

        Assert.assertEquals(c.internalCompare(C, Q30, C, Q60), 1);
        Assert.assertEquals(c.internalCompare(C, Q60, C, Q30), -1);
        Assert.assertEquals(c.internalCompare(C, Q60, C, Q60), 0);

        Assert.assertEquals(c.internalCompare(A, Q30, C, Q60), -1);
        Assert.assertEquals(c.internalCompare(A, Q60, C, Q30), -1);
        Assert.assertEquals(c.internalCompare(A, Q60, C, Q60), -1);

        Assert.assertEquals(c.internalCompare(C, Q30, A, Q60), 1);
        Assert.assertEquals(c.internalCompare(C, Q60, A, Q30), 1);
        Assert.assertEquals(c.internalCompare(C, Q60, A, Q60), 1);

        Assert.assertEquals(c.internalCompare(C, Q60, T, Q60), -1);
        Assert.assertEquals(c.internalCompare(T, Q60, C, Q60), 1);

    }

    public static int getBinomial(int n, double p) {
        int x = 0;
        for(int i = 0; i < n; i++) {
            if(Math.random() < p)
                x++;
        }
        return x;
    }

    private double getGaussian(Random r, double aMean, double aVariance){
        return aMean + r.nextGaussian() * aVariance;
    }

    private void dumpArray(int[] a) {
        System.out.println("---------------");
        for(int i : a) System.out.println(i);
        System.out.println("---------------");

    }
}
