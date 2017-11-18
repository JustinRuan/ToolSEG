package edu.whut.info.analysis;

import edu.whut.info.dataset.Segment;
import edu.whut.info.util.BioToolbox;

import java.util.*;
import java.util.logging.Logger;

/**
 * Created by Justin on 2016/10/25.
 */
public class ResultAnalysis {
    private Logger m_log;

    private List<Set<Segment>> correctSegs;
    private List<Set<Segment>> results;
    private List<Integer> lengths;
    private List<List<Long>> breakPoints;
    private double[] P;
    private double[] N;
    private double[] TP;
    private double[] FP;

    private List<Double> goldStandard;
    private List<Double> testValues;

    public ResultAnalysis() {
        m_log = Logger.getLogger("segment");
        correctSegs = new LinkedList<>();
        results = new LinkedList<>();
        lengths = new LinkedList<>();
        breakPoints = new ArrayList<>();
        goldStandard = new LinkedList<>();
        testValues = new LinkedList<>();
    }

    public void clear() {
        correctSegs.clear();
        results.clear();
        lengths.clear();
        goldStandard.clear();
        testValues.clear();
    }

    public void addBreakPoints(List<Long> breakPoints) {
        this.breakPoints.add(breakPoints);
    }

    public void addStandard(short chrId, double[] values, List<Long> breakPoints) {
        Set<Segment> tmpSet = new TreeSet<>();

        int begin, end;
        begin = 0;
        end = 0;
        for (long t : breakPoints) {
            end = (int) t;
            Segment seg = new Segment();
            seg.setChr_id(chrId);
            seg.setRange(begin, end);
            BioToolbox.refreshSegment(seg, values);
            tmpSet.add(seg);

            begin = end;
        }

        begin = end;
        end = values.length;
        Segment seg = new Segment();
        seg.setChr_id(chrId);
        seg.setRange(begin, end);
        BioToolbox.refreshSegment(seg, values);
        tmpSet.add(seg);

        correctSegs.add(tmpSet);
        lengths.add(end);
    }

    public void addResult(Set<Segment> result) {
//        for (Segment seg : result){
//            if (seg.isLog2){
//                seg.CopyNumber = Math.pow(2,seg.CopyNumber);
//            }
//        }
        results.add(result);
    }

    public void prepareTest() {
        Iterator<Set<Segment>> itCorrSeg = correctSegs.iterator();
        Iterator<Set<Segment>> itResult = results.iterator();
        Iterator<Integer> itLength = lengths.iterator();
        for (int i = 0; i < lengths.size(); i++) {
            Set<Segment> corrSeg = itCorrSeg.next();
            Set<Segment> result = itResult.next();
            for (Segment seg : corrSeg) {
                int start = seg.Start();
                int end = seg.End();
                int step = 5;
                while (step <= 20) {
                    int pos = start + step;
                    double trueValue = getValue(pos, corrSeg);
                    //boolean tmpTrue = (Math.abs(trueValue - 2) < 0.02);
                    if (start != 0) {
                        goldStandard.add(trueValue);
                        testValues.add(getValue(pos, result));
                    }
                    pos = end - step;
                    trueValue = getValue(pos, corrSeg);
                    //boolean tmpTrue = (Math.abs(trueValue - 2) < 0.02);
                    if (end != lengths.get(i)) {
                        goldStandard.add(trueValue);
                        testValues.add(getValue(pos, result));
                    }
                    step += 5;
                }
            }
        }
//        Random rnd = new Random(100);
//
//        while (itLength.hasNext()){
//            int len = itLength.next();
//            Set<Segment> corrSeg = itCorrSeg.next();
//            Set<Segment> result = itResult.next();
//
//            int pos = rnd.nextInt(20);
//            while (pos < len){
//                double trueValue = getValue(pos, corrSeg);
//                //boolean tmpTrue = (Math.abs(trueValue - 2) < 0.02);
//                goldStandard.add(trueValue);
//
//                testValues.add(getValue(pos, result));
//
//                pos += 50 + rnd.nextInt(200);
//            }
    }

    public void prepareTest1() {
        Iterator<Set<Segment>> itCorrSeg = correctSegs.iterator();
        Iterator<Set<Segment>> itResult = results.iterator();
        Iterator<Integer> itLength = lengths.iterator();

        Random rnd = new Random(100);

        while (itLength.hasNext()) {
            int len = itLength.next();
            Set<Segment> corrSeg = itCorrSeg.next();
            Set<Segment> result = itResult.next();

            int pos = rnd.nextInt(20);
            while (pos < len) {
                double trueValue = getValue(pos, corrSeg);
                //boolean tmpTrue = (Math.abs(trueValue - 2) < 0.02);
                goldStandard.add(trueValue);

                testValues.add(getValue(pos, result));

                pos += 50 + rnd.nextInt(200);
            }
        }
    }

    private double getValue(int pos, Set<Segment> segments) {
        for (Segment seg : segments) {
            if (seg.isIncluded(pos)) {
                return seg.CopyNumber;
            }
        }
        System.out.println("Error! getValue(int pos, Set<Segment> segments)");
        return -1;
    }

    public void analysisResult() {
        //int sampleCount = goldStandard.size();
        final int step = 20;
        P = new double[step];
        N = new double[step];
        TP = new double[step];
        FP = new double[step];

        m_log.info(String.format("\n......Roc date ......\n"));

        for (int k = 0; k < step; k++) {
            double Tolerance = 0.0001 + 0.0001 * k * k * k;
//            double Tolerance = 0.001 + 0.001 * k * k * k;
            Iterator<Double> itTrue = goldStandard.iterator();
            Iterator<Double> itValue = testValues.iterator();
            while (itTrue.hasNext()) {
                double trueValue, testValue;
                trueValue = itTrue.next();
                testValue = itValue.next();

                boolean trueTag = (Math.abs(trueValue - 2) < 0.05);
                boolean testTag = (Math.abs(testValue - 2) < Tolerance);

                if (trueTag) {
                    P[k]++;
                    if (testTag) {
                        TP[k]++;
                    }
                } else {
                    N[k]++;
                    if (testTag) {
                        FP[k]++;
                    }
                }
            }

            double fpr = FP[k] / N[k];
            double tpr = TP[k] / P[k];
            m_log.info(String.format("\t% 4d\t%f\t%f", k, fpr, tpr));
        }
    }

    public void correctSegStatistics(){
        m_log.info(String.format("\n......correctSegs Statistics ......sample count = %d \n ",
                correctSegs.size()));

        int count = 0;
        StringBuilder result = new StringBuilder();
        for (Set<Segment> setSeg : correctSegs){
            count += setSeg.size();
            for (Segment seg :setSeg){
                result.append(String.format("\t%.4f\t%.4f\t%d\n",seg.CopyNumber,seg.stdCopyNumber,seg.length()));
            }
        }

        m_log.info(String.format("\n. Segment count = %d \n ",count));
        m_log.info("\n" + result.toString());

    }
}
