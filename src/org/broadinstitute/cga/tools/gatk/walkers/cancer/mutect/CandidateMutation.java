package org.broadinstitute.cga.tools.gatk.walkers.cancer.mutect;

import org.broadinstitute.sting.utils.GenomeLoc;
import org.broadinstitute.sting.utils.variantcontext.VariantContext;
import org.broadinstitute.sting.gatk.walkers.genotyper.DiploidGenotype;

import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;

public class CandidateMutation {
    private GenomeLoc location;
    private char refAllele;
    private boolean dbsnpSite = false;
    private boolean cosmicSite = false;
    private VariantContext panelOfNormalsVC;
    private boolean covered = false;

    private double power;
    private double tumorPower;
    private double normalPower;
    private double normalPowerWithSNPPrior;
    private double normalPowerNoSNPPrior;

    private char altAllele = 'N';
    private char priorBasePositiveDirection;
    private char priorBaseNegativeDirection;

    private boolean positiveDirectionPowered;
    private boolean negativeDirectionPowered;
    private String tumorSampleName;
    private String normalSampleName;

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
    private int initialNormalReadDepth;
    private DiploidGenotype initialNormalBestGenotype;

    private double initialTumorLod;
    private double initialNormalLod;

    private double tumorF;
    private double tumorLodFStar;
    private double tumorLodFStarForward;
    private double tumorLodFStarReverse;

    private double normalF;
    private double normalLodFStar;

    private double normalArtifactPower;
    private double normalArtifactLod;

    private RankSumTest.Result tumorQualityRankSumTest;
    private RankSumTest.Result tumorReadPositionRankSumTest;
    private Map<Integer, Double> classicSkewScoresAndOffsets;
    private Map<Integer, Double> fisherSkewScoresAndOffsets;

    private double powerToDetectPositiveStrandArtifact;
    private double powerToDetectNegativeStrandArtifact;

    private MuTectWalker.FisherData strandBias;
    private MuTectWalker.FisherData perfectStrandBias;

    private MuTectWalker.FisherData clippingBias;

    private List<Integer> tumorAltForwardOffsetsInRead;
    private List<Integer> tumorAltReverseOffsetsInRead;

    private Double tumorForwardOffsetsInReadMedian;
    private Double tumorForwardOffsetsInReadMad;
    private Double tumorReverseOffsetsInReadMedian;
    private Double tumorReverseOffsetsInReadMad;

    private int tumorInsertionCount;
    private int tumorDeletionCount;

    private List<String> poweredFilters = new ArrayList<String>();

    private List<String> rejectionReasons = new ArrayList<String>();
    private boolean rejected = false; // summary judgement... keep or reject the site

    public CandidateMutation(GenomeLoc location, char refAllele) {
        this.location = location;
        this.refAllele = refAllele;
    }

    public int getScore() {
        return score;
    }

    public boolean isPositiveDirectionAtRisk() {
        return getPriorBasePositiveDirection() == getAltAllele();
    }

    public boolean isNegativeDirectionAtRisk() {
        return getPriorBaseNegativeDirection() == getAltAllele();
    }
    
    public boolean isGermlineAtRisk() {
        return (dbsnpSite && !cosmicSite);
    }

    // -------------------------------------------------------------------------
    // GENERATED CODE BELOW THIS POINT
    // -------------------------------------------------------------------------


    public boolean isPositiveDirectionPowered() {
        return positiveDirectionPowered;
    }

    public void setPositiveDirectionPowered(boolean positiveDirectionPowered) {
        this.positiveDirectionPowered = positiveDirectionPowered;
    }

    public boolean isNegativeDirectionPowered() {
        return negativeDirectionPowered;
    }

    public void setNegativeDirectionPowered(boolean negativeDirectionPowered) {
        this.negativeDirectionPowered = negativeDirectionPowered;
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

    public boolean isCovered() {
        return covered;
    }

    public void setCovered(boolean covered) {
        this.covered = covered;
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

    public char getPriorBasePositiveDirection() {
        return priorBasePositiveDirection;
    }

    public void setPriorBasePositiveDirection(char priorBasePositiveDirection) {
        this.priorBasePositiveDirection = priorBasePositiveDirection;
    }

    public char getPriorBaseNegativeDirection() {
        return priorBaseNegativeDirection;
    }

    public void setPriorBaseNegativeDirection(char priorBaseNegativeDirection) {
        this.priorBaseNegativeDirection = priorBaseNegativeDirection;
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

    public double getNormalLodFStar() {
        return normalLodFStar;
    }

    public void setNormalLodFStar(double normalLodFStar) {
        this.normalLodFStar = normalLodFStar;
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

    public RankSumTest.Result getTumorQualityRankSumTest() {
        return tumorQualityRankSumTest;
    }

    public void setTumorQualityRankSumTest(RankSumTest.Result tumorQualityRankSumTest) {
        this.tumorQualityRankSumTest = tumorQualityRankSumTest;
    }

    public RankSumTest.Result getTumorReadPositionRankSumTest() {
        return tumorReadPositionRankSumTest;
    }

    public void setTumorReadPositionRankSumTest(RankSumTest.Result tumorReadPositionRankSumTest) {
        this.tumorReadPositionRankSumTest = tumorReadPositionRankSumTest;
    }

    public MuTectWalker.FisherData getStrandBias() {
        return strandBias;
    }

    public void setStrandBias(MuTectWalker.FisherData strandBias) {
        this.strandBias = strandBias;
    }

    public MuTectWalker.FisherData getPerfectStrandBias() {
        return perfectStrandBias;
    }

    public void setPerfectStrandBias(MuTectWalker.FisherData perfectStrandBias) {
        this.perfectStrandBias = perfectStrandBias;
    }

    public MuTectWalker.FisherData getClippingBias() {
        return clippingBias;
    }

    public void setClippingBias(MuTectWalker.FisherData clippingBias) {
        this.clippingBias = clippingBias;
    }

    public List<String> getRejectionReasons() {
        return rejectionReasons;
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


    public List<String> getPoweredFilters() {
        return poweredFilters;
    }

    public void setPoweredFilters(List<String> poweredFilters) {
        this.poweredFilters = poweredFilters;
    }

    public void addPoweredFilter(String filter) {
        getPoweredFilters().add(filter);        
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

    public Map<Integer, Double> getClassicSkewScoresAndOffsets() {
        return classicSkewScoresAndOffsets;
    }

    public void setClassicSkewScoresAndOffsets(Map<Integer, Double> classicSkewScoresAndOffsets) {
        this.classicSkewScoresAndOffsets = classicSkewScoresAndOffsets;
    }

    public Map<Integer, Double> getFisherSkewScoresAndOffsets() {
        return fisherSkewScoresAndOffsets;
    }

    public void setFisherSkewScoresAndOffsets(Map<Integer, Double> fisherSkewScoresAndOffsets) {
        this.fisherSkewScoresAndOffsets = fisherSkewScoresAndOffsets;
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

    public VariantContext getPanelOfNormalsVC() {
        return panelOfNormalsVC;
    }

    public void setPanelOfNormalsVC(VariantContext panelOfNormalsVC) {
        this.panelOfNormalsVC = panelOfNormalsVC;
    }

    public double getNormalArtifactPower() {
        return normalArtifactPower;
    }

    public void setNormalArtifactPower(double normalArtifactPower) {
        this.normalArtifactPower = normalArtifactPower;
    }

    public double getNormalArtifactLod() {
        return normalArtifactLod;
    }

    public void setNormalArtifactLod(double normalArtifactLod) {
        this.normalArtifactLod = normalArtifactLod;
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
}

