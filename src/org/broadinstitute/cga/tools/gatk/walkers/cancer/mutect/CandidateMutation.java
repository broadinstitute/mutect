/*
 * By downloading the PROGRAM you agree to the following terms of use:
 *
 * BROAD INSTITUTE SOFTWARE LICENSE AGREEMENT
 * FOR ACADEMIC NON-COMMERCIAL RESEARCH PURPOSES ONLY
 *
 * This Agreement is made between the Broad Institute, Inc. with a principal address at 7 Cambridge Center, Cambridge, MA 02142 ("BROAD") and the LICENSEE and is effective at the date the downloading is completed ("EFFECTIVE DATE").
 * WHEREAS, LICENSEE desires to license the PROGRAM, as defined hereinafter, and BROAD wishes to have this PROGRAM utilized in the public interest, subject only to the royalty-free, nonexclusive, nontransferable license rights of the United States Government pursuant to 48 CFR 52.227-14; and
 * WHEREAS, LICENSEE desires to license the PROGRAM and BROAD desires to grant a license on the following terms and conditions.
 * NOW, THEREFORE, in consideration of the promises and covenants made herein, the parties hereto agree as follows:
 *
 * 1. DEFINITIONS
 * 1.1	"PROGRAM" shall mean copyright in the object code and source code known as MuTect and related documentation, if any, as they exist on the EFFECTIVE DATE and can be downloaded from http://www.broadinstitute.or/cancer/cga/mutect on the EFFECTIVE DATE.
 *
 * 2. LICENSE
 * 2.1   Grant. Subject to the terms of this Agreement, BROAD hereby grants to LICENSEE, solely for academic non-commercial research purposes, a non-exclusive, non-transferable license to: (a) download, execute and display the PROGRAM and (b) create bug fixes and modify the PROGRAM.
 * LICENSEE hereby automatically grants to BROAD a non-exclusive, royalty-free, irrevocable license to any LICENSEE bug fixes or modifications to the PROGRAM with unlimited rights to sublicense and/or distribute.  LICENSEE agrees to provide any such modifications and bug fixes to BROAD promptly upon their creation.
 * The LICENSEE may apply the PROGRAM in a pipeline to data owned by users other than the LICENSEE and provide these users the results of the PROGRAM provided LICENSEE does so for academic non-commercial purposes only.  For clarification purposes, academic sponsored research is not a commercial use under the terms of this Agreement.
 * 2.2  No Sublicensing or Additional Rights. LICENSEE shall not sublicense or distribute the PROGRAM, in whole or in part, without prior written permission from BROAD.  LICENSEE shall ensure that all of its users agree to the terms of this Agreement.  LICENSEE further agrees that it shall not put the PROGRAM on a network, server, or other similar technology that may be accessed by anyone other than the LICENSEE and its employees and users who have agreed to the terms of this agreement.
 * 2.3  License Limitations. Nothing in this Agreement shall be construed to confer any rights upon LICENSEE by implication, estoppel, or otherwise to any computer software, trademark, intellectual property, or patent rights of BROAD, or of any other entity, except as expressly granted herein. LICENSEE agrees that the PROGRAM, in whole or part, shall not be used for any commercial purpose, including without limitation, as the basis of a commercial software or hardware product or to provide services. LICENSEE further agrees that the PROGRAM shall not be copied or otherwise adapted in order to circumvent the need for obtaining a license for use of the PROGRAM.
 *
 * 3. OWNERSHIP OF INTELLECTUAL PROPERTY
 * LICENSEE acknowledges that title to the PROGRAM shall remain with BROAD. The PROGRAM is marked with the following BROAD copyright notice and notice of attribution to contributors. LICENSEE shall retain such notice on all copies.  LICENSEE agrees to include appropriate attribution if any results obtained from use of the PROGRAM are included in any publication.
 *
 * Copyright 2012 Broad Institute, Inc.
 * Notice of attribution:  The MuTect program was made available through the generosity of the Cancer Genome Analysis group at the Broad Institute, Inc. and is published at doi: 10.1038/nbt.2514.
 *
 * LICENSEE shall not use any trademark or trade name of BROAD, or any variation, adaptation, or abbreviation, of such marks or trade names, or any names of officers, faculty, students, employees, or agents of BROAD except as states above for attribution purposes.
 *
 * 4. INDEMNIFICATION
 * LICENSEE shall indemnify, defend, and hold harmless BROAD, and their respective officers, faculty, students, employees, associated investigators and agents, and their respective successors, heirs and assigns, ("Indemnitees"), against any liability, damage, loss, or expense (including reasonable attorneys fees and expenses) incurred by or imposed upon any of the Indemnitees in connection with any claims, suits, actions, demands or judgments arising out of any theory of liability (including, without limitation, actions in the form of tort, warranty, or strict liability and regardless of whether such action has any factual basis) pursuant to any right or license granted under this Agreement.
 *
 * 5. NO REPRESENTATIONS OR WARRANTIES
 * THE PROGRAM IS DELIVERED "AS IS."  BROAD MAKES NO REPRESENTATIONS OR WARRANTIES OF ANY KIND CONCERNING THE PROGRAM OR THE COPYRIGHT, EXPRESS OR IMPLIED, INCLUDING, WITHOUT LIMITATION, WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, NONINFRINGEMENT, OR THE ABSENCE OF LATENT OR OTHER DEFECTS, WHETHER OR NOT DISCOVERABLE. BROAD EXTENDS NO WARRANTIES OF ANY KIND AS TO PROGRAM CONFORMITY WITH WHATEVER USER MANUALS OR OTHER LITERATURE MAY BE ISSUED FROM TIME TO TIME.
 * IN NO EVENT SHALL BROAD OR ITS RESPECTIVE DIRECTORS, OFFICERS, EMPLOYEES, AFFILIATED INVESTIGATORS AND AFFILIATES BE LIABLE FOR INCIDENTAL OR CONSEQUENTIAL DAMAGES OF ANY KIND, INCLUDING, WITHOUT LIMITATION, ECONOMIC DAMAGES OR INJURY TO PROPERTY AND LOST PROFITS, REGARDLESS OF WHETHER BROAD SHALL BE ADVISED, SHALL HAVE OTHER REASON TO KNOW, OR IN FACT SHALL KNOW OF THE POSSIBILITY OF THE FOREGOING.
 *
 * 6. ASSIGNMENT
 * This Agreement is personal to LICENSEE and any rights or obligations assigned by LICENSEE without the prior written consent of BROAD shall be null and void.
 *
 * 7. MISCELLANEOUS
 * 7.1 Export Control. LICENSEE gives assurance that it will comply with all United States export control laws and regulations controlling the export of the PROGRAM, including, without limitation, all Export Administration Regulations of the United States Department of Commerce. Among other things, these laws and regulations prohibit, or require a license for, the export of certain types of software to specified countries.
 * 7.2 Termination. LICENSEE shall have the right to terminate this Agreement for any reason upon prior written notice to BROAD. If LICENSEE breaches any provision hereunder, and fails to cure such breach within thirty (30) days, BROAD may terminate this Agreement immediately. Upon termination, LICENSEE shall provide BROAD with written assurance that the original and all copies of the PROGRAM have been destroyed, except that, upon prior written authorization from BROAD, LICENSEE may retain a copy for archive purposes.
 * 7.3 Survival. The following provisions shall survive the expiration or termination of this Agreement: Articles 1, 3, 4, 5 and Sections 2.2, 2.3, 7.3, and 7.4.
 * 7.4 Notice. Any notices under this Agreement shall be in writing, shall specifically refer to this Agreement, and shall be sent by hand, recognized national overnight courier, confirmed facsimile transmission, confirmed electronic mail, or registered or certified mail, postage prepaid, return receipt requested.  All notices under this Agreement shall be deemed effective upon receipt.
 * 7.5 Amendment and Waiver; Entire Agreement. This Agreement may be amended, supplemented, or otherwise modified only by means of a written instrument signed by all parties. Any waiver of any rights or failure to act in a specific instance shall relate only to such instance and shall not be construed as an agreement to waive any rights or fail to act in any other instance, whether or not similar. This Agreement constitutes the entire agreement among the parties with respect to its subject matter and supersedes prior agreements or understandings between the parties relating to its subject matter.
 * 7.6 Binding Effect; Headings. This Agreement shall be binding upon and inure to the benefit of the parties and their respective permitted successors and assigns. All headings are for convenience only and shall not affect the meaning of any provision of this Agreement.
 * 7.7 Governing Law. This Agreement shall be construed, governed, interpreted and applied in accordance with the internal laws of the Commonwealth of Massachusetts, U.S.A., without regard to conflict of laws principles.
 */

package org.broadinstitute.cga.tools.gatk.walkers.cancer.mutect;

import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;
import org.broadinstitute.sting.gatk.walkers.genotyper.DiploidGenotype;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;

public class CandidateMutation {
    private GenomeLoc location;
    private String sequenceContext;
    private char refAllele;
    private boolean dbsnpSite = false;
    private boolean cosmicSite = false;
    private VariantContext panelOfNormalsVC;
    private VariantContext dbsnpVC;
    private boolean covered = false;

    private double power;
    private double tumorPower;
    private double normalPower;
    private double normalPowerWithSNPPrior;
    private double normalPowerNoSNPPrior;

    private char altAllele = 'N';
    private String tumorSampleName = "TUMOR";
    private String normalSampleName = "NORMAL";

    private double contaminationFraction;

    private double contaminantLod;
    private int score;

    private int tumorQ20Count;
    private int normalQ20Count;

    private int totalPairs;
    private int improperPairs;
    private int mapQ0Reads;
    private int initialTumorRefCounts;
    private int initialTumorAltCounts;
    private int initialTumorRefQualitySum;
    private int initialTumorAltQualitySum;
    private int initialTumorNonRefQualitySum;
    private int initialTumorReadDepth;
    private int initialNormalRefCounts;
    private int initialNormalAltCounts;
    private int initialNormalRefQualitySum;
    private int initialNormalAltQualitySum;
    private int tumorRefMaxMapQ;
    private int tumorAltMaxMapQ;
    private int initialNormalReadDepth;
    private DiploidGenotype initialNormalBestGenotype;

    private double initialTumorLod;
    private double initialNormalLod;

    private double tumorF;
    private double tumorLodFStar;
    private double tumorLodFStarForward;
    private double tumorLodFStarReverse;

    private double normalF;

    private double powerToDetectPositiveStrandArtifact;
    private double powerToDetectNegativeStrandArtifact;

    private int[] strandContingencyTable;

    private List<Integer> tumorAltForwardOffsetsInRead;
    private List<Integer> tumorAltReverseOffsetsInRead;

    private Double tumorForwardOffsetsInReadMedian;
    private Double tumorForwardOffsetsInReadMad;
    private Double tumorReverseOffsetsInReadMedian;
    private Double tumorReverseOffsetsInReadMad;

    private int tumorInsertionCount;
    private int tumorDeletionCount;

    private List<String> rejectionReasons = new ArrayList<String>();
    private boolean rejected = false; // summary judgement... keep or reject the site

    public CandidateMutation(GenomeLoc location, char refAllele) {
        this.location = location;
        this.refAllele = refAllele;
    }

    public int getScore() {
        return score;
    }

    public boolean isGermlineAtRisk() {
        return (dbsnpSite && !cosmicSite);
    }

    public int getCountOfNormalsObservedIn() {
        if (panelOfNormalsVC == null) { return 0; }

        // it's either an integer, or a collection of integers
        int count = 0;
        Object o = panelOfNormalsVC.getAttribute("AC");
        if (o == null) { return 0; }
        if (o instanceof String) { return Integer.parseInt((String) o); }
        if (o instanceof Collection) {
            for(String s : ((Collection<String>) o)) {
                count += Integer.parseInt(s);
            }
            return count;
        }
        throw new RuntimeException("Unexpected value processing panel of normals allele count: " + o);
    }

    public void setRejectionReasons(List<String> rejectionReasons) {
        if (rejectionReasons != null && rejectionReasons.size() > 0) {
            setRejected(true);
        }
        this.rejectionReasons = rejectionReasons;
    }

    public void addRejectionReason(String reason) {
        setRejected(true);
        getRejectionReasons().add(reason);
    }




    // -------------------------------------------------------------------------
    // GENERATED CODE BELOW THIS POINT
    // -------------------------------------------------------------------------

    public int[] getStrandContingencyTable() {
        return strandContingencyTable;
    }

    public void setStrandContingencyTable(int[] strandContingencyTable) {
        this.strandContingencyTable = strandContingencyTable;
    }

    public GenomeLoc getLocation() {
        return location;
    }

    public boolean isDbsnpSite() {
        return dbsnpSite;
    }

    public void setDbsnpSite(boolean dbsnpSite) {
        this.dbsnpSite = dbsnpSite;
    }

    public VariantContext getDbsnpVC() {
        return dbsnpVC;
    }

    public void setDbsnpVC(VariantContext dbsnpVC) {
        this.dbsnpVC = dbsnpVC;
    }

    public boolean isCovered() {
        return covered;
    }

    public void setCovered(boolean covered) {
        this.covered = covered;
    }

    public String getSequenceContext() {
        return sequenceContext;
    }

    public void setSequenceContext(String sequenceContext) {
        this.sequenceContext = sequenceContext;
    }

    public char getRefAllele() {
        return refAllele;
    }

    public char getAltAllele() {
        return altAllele;
    }

    public void setAltAllele(char altAllele) {
        this.altAllele = altAllele;
    }

    public boolean isRejected() {
        return rejected;
    }

    public void setRejected(boolean rejected) {
        this.rejected = rejected;
    }

    public double getInitialTumorLod() {
        return initialTumorLod;
    }

    public void setInitialTumorLod(double initialTumorLod) {
        this.initialTumorLod = initialTumorLod;
    }

    public double getInitialNormalLod() {
        return initialNormalLod;
    }

    public void setInitialNormalLod(double initialNormalLod) {
        this.initialNormalLod = initialNormalLod;
    }

    public double getTumorLodFStar() {
        return tumorLodFStar;
    }

    public void setTumorLodFStar(double tumorLodFStar) {
        this.tumorLodFStar = tumorLodFStar;
    }


    public double getTumorLodFStarForward() {
        return tumorLodFStarForward;
    }

    public void setTumorLodFStarForward(double tumorLodFStarForward) {
        this.tumorLodFStarForward = tumorLodFStarForward;
    }

    public double getTumorLodFStarReverse() {
        return tumorLodFStarReverse;
    }

    public void setTumorLodFStarReverse(double tumorLodFStarReverse) {
        this.tumorLodFStarReverse = tumorLodFStarReverse;
    }

    public double getTumorF() {
        return tumorF;
    }

    public void setTumorF(double tumorF) {
        this.tumorF = tumorF;
    }

    public double getNormalF() {
        return normalF;
    }

    public void setNormalF(double normalF) {
        this.normalF = normalF;
    }

    public int getInitialTumorRefQualitySum() {
        return initialTumorRefQualitySum;
    }

    public void setInitialTumorRefQualitySum(int initialTumorRefQualitySum) {
        this.initialTumorRefQualitySum = initialTumorRefQualitySum;
    }

    public int getInitialTumorAltQualitySum() {
        return initialTumorAltQualitySum;
    }

    public void setInitialTumorAltQualitySum(int initialTumorAltQualitySum) {
        this.initialTumorAltQualitySum = initialTumorAltQualitySum;
    }

    public int getInitialTumorNonRefQualitySum() {
        return initialTumorNonRefQualitySum;
    }

    public void setInitialTumorNonRefQualitySum(int initialTumorNonRefQualitySum) {
        this.initialTumorNonRefQualitySum = initialTumorNonRefQualitySum;
    }

    public int getInitialNormalRefQualitySum() {
        return initialNormalRefQualitySum;
    }

    public void setInitialNormalRefQualitySum(int initialNormalRefQualitySum) {
        this.initialNormalRefQualitySum = initialNormalRefQualitySum;
    }

    public int getInitialNormalAltQualitySum() {
        return initialNormalAltQualitySum;
    }

    public void setInitialNormalAltQualitySum(int initialNormalAltQualitySum) {
        this.initialNormalAltQualitySum = initialNormalAltQualitySum;
    }

    public DiploidGenotype getInitialNormalBestGenotype() {
        return initialNormalBestGenotype;
    }

    public void setInitialNormalBestGenotype(DiploidGenotype initialNormalBestGenotype) {
        this.initialNormalBestGenotype = initialNormalBestGenotype;
    }

    public int getInitialTumorReadDepth() {
        return initialTumorReadDepth;
    }

    public void setInitialTumorReadDepth(int initialTumorReadDepth) {
        this.initialTumorReadDepth = initialTumorReadDepth;
    }

    public int getInitialNormalReadDepth() {
        return initialNormalReadDepth;
    }

    public void setInitialNormalReadDepth(int initialNormalReadDepth) {
        this.initialNormalReadDepth = initialNormalReadDepth;
    }

    public String getTumorSampleName() {
        return tumorSampleName;
    }

    public void setTumorSampleName(String tumorSampleName) {
        this.tumorSampleName = tumorSampleName;
    }

    public String getNormalSampleName() {
        return normalSampleName;
    }

    public void setNormalSampleName(String normalSampleName) {
        this.normalSampleName = normalSampleName;
    }

    public List<String> getRejectionReasons() {
        return rejectionReasons;
    }

    public int getInitialTumorRefCounts() {
        return initialTumorRefCounts;
    }

    public void setInitialTumorRefCounts(int initialTumorRefCounts) {
        this.initialTumorRefCounts = initialTumorRefCounts;
    }

    public int getInitialTumorAltCounts() {
        return initialTumorAltCounts;
    }

    public void setInitialTumorAltCounts(int initialTumorAltCounts) {
        this.initialTumorAltCounts = initialTumorAltCounts;
    }

    public int getInitialNormalRefCounts() {
        return initialNormalRefCounts;
    }

    public void setInitialNormalRefCounts(int initialNormalRefCounts) {
        this.initialNormalRefCounts = initialNormalRefCounts;
    }

    public int getInitialNormalAltCounts() {
        return initialNormalAltCounts;
    }

    public void setInitialNormalAltCounts(int initialNormalAltCounts) {
        this.initialNormalAltCounts = initialNormalAltCounts;
    }

    public int getTumorQ20Count() {
        return tumorQ20Count;
    }

    public void setTumorQ20Count(int tumorQ20Count) {
        this.tumorQ20Count = tumorQ20Count;
    }

    public int getNormalQ20Count() {
        return normalQ20Count;
    }

    public void setNormalQ20Count(int normalQ20Count) {
        this.normalQ20Count = normalQ20Count;
    }

    public int getTumorInsertionCount() {
        return tumorInsertionCount;
    }

    public void setTumorInsertionCount(int tumorInsertionCount) {
        this.tumorInsertionCount = tumorInsertionCount;
    }

    public int getTumorDeletionCount() {
        return tumorDeletionCount;
    }

    public void setTumorDeletionCount(int tumorDeletionCount) {
        this.tumorDeletionCount = tumorDeletionCount;
    }

    public int getTotalPairs() {
        return totalPairs;
    }

    public void setTotalPairs(int totalPairs) {
        this.totalPairs = totalPairs;
    }

    public int getImproperPairs() {
        return improperPairs;
    }

    public void setImproperPairs(int improperPairs) {
        this.improperPairs = improperPairs;
    }

    public int getMapQ0Reads() {
        return mapQ0Reads;
    }

    public void setMapQ0Reads(int mapQ0Reads) {
        this.mapQ0Reads = mapQ0Reads;
    }

    public double getContaminationFraction() {
        return contaminationFraction;
    }

    public void setContaminationFraction(double contaminationFraction) {
        this.contaminationFraction = contaminationFraction;
    }

    public double getContaminantLod() {
        return contaminantLod;
    }

    public void setContaminantLod(double contaminantLod) {
        this.contaminantLod = contaminantLod;
    }

    public List<Integer> getTumorAltForwardOffsetsInRead() {
        return tumorAltForwardOffsetsInRead;
    }

    public void setTumorAltForwardOffsetsInRead(List<Integer> tumorAltForwardOffsetsInRead) {
        this.tumorAltForwardOffsetsInRead = tumorAltForwardOffsetsInRead;
    }

    public List<Integer> getTumorAltReverseOffsetsInRead() {
        return tumorAltReverseOffsetsInRead;
    }

    public void setTumorAltReverseOffsetsInRead(List<Integer> tumorAltReverseOffsetsInRead) {
        this.tumorAltReverseOffsetsInRead = tumorAltReverseOffsetsInRead;
    }

    public Double getTumorForwardOffsetsInReadMedian() {
        return tumorForwardOffsetsInReadMedian;
    }

    public void setTumorForwardOffsetsInReadMedian(Double tumorForwardOffsetsInReadMedian) {
        this.tumorForwardOffsetsInReadMedian = tumorForwardOffsetsInReadMedian;
    }

    public Double getTumorForwardOffsetsInReadMad() {
        return tumorForwardOffsetsInReadMad;
    }

    public void setTumorForwardOffsetsInReadMad(Double tumorForwardOffsetsInReadMad) {
        this.tumorForwardOffsetsInReadMad = tumorForwardOffsetsInReadMad;
    }

    public Double getTumorReverseOffsetsInReadMedian() {
        return tumorReverseOffsetsInReadMedian;
    }

    public void setTumorReverseOffsetsInReadMedian(Double tumorReverseOffsetsInReadMedian) {
        this.tumorReverseOffsetsInReadMedian = tumorReverseOffsetsInReadMedian;
    }

    public Double getTumorReverseOffsetsInReadMad() {
        return tumorReverseOffsetsInReadMad;
    }

    public void setTumorReverseOffsetsInReadMad(Double tumorReverseOffsetsInReadMad) {
        this.tumorReverseOffsetsInReadMad = tumorReverseOffsetsInReadMad;
    }

    public double getPower() {
        return power;
    }

    public void setPower(double power) {
        this.power = power;
    }

    public double getTumorPower() {
        return tumorPower;
    }

    public void setTumorPower(double tumorPower) {
        this.tumorPower = tumorPower;
    }

    public double getNormalPower() {
        return normalPower;
    }

    public void setNormalPower(double normalPower) {
        this.normalPower = normalPower;
    }

    public boolean isCosmicSite() {
        return cosmicSite;
    }

    public void setCosmicSite(boolean cosmicSite) {
        this.cosmicSite = cosmicSite;
    }

    public boolean isSeenInPanelOfNormals() {
        return (panelOfNormalsVC != null);
    }

    public VariantContext getPanelOfNormalsVC() {
        return panelOfNormalsVC;
    }

    public void setPanelOfNormalsVC(VariantContext panelOfNormalsVC) {
        this.panelOfNormalsVC = panelOfNormalsVC;
    }

    public double getPowerToDetectPositiveStrandArtifact() {
        return powerToDetectPositiveStrandArtifact;
    }

    public void setPowerToDetectPositiveStrandArtifact(double powerToDetectPositiveStrandArtifact) {
        this.powerToDetectPositiveStrandArtifact = powerToDetectPositiveStrandArtifact;
    }

    public double getPowerToDetectNegativeStrandArtifact() {
        return powerToDetectNegativeStrandArtifact;
    }

    public void setPowerToDetectNegativeStrandArtifact(double powerToDetectNegativeStrandArtifact) {
        this.powerToDetectNegativeStrandArtifact = powerToDetectNegativeStrandArtifact;
    }

    public double getNormalPowerWithSNPPrior() {
        return normalPowerWithSNPPrior;
    }

    public void setNormalPowerWithSNPPrior(double normalPowerWithSNPPrior) {
        this.normalPowerWithSNPPrior = normalPowerWithSNPPrior;
    }

    public double getNormalPowerNoSNPPrior() {
        return normalPowerNoSNPPrior;
    }

    public void setNormalPowerNoSNPPrior(double normalPowerNoSNPPrior) {
        this.normalPowerNoSNPPrior = normalPowerNoSNPPrior;
    }

    public int getTumorAltMaxMapQ() {
        return tumorAltMaxMapQ;
    }

    public void setTumorAltMaxMapQ(int tumorAltMaxMapQ) {
        this.tumorAltMaxMapQ = tumorAltMaxMapQ;
    }

    public int getTumorRefMaxMapQ() {
        return tumorRefMaxMapQ;
    }

    public void setTumorRefMaxMapQ(int tumorRefMaxMapQ) {
        this.tumorRefMaxMapQ = tumorRefMaxMapQ;
    }
}

