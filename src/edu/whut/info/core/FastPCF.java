package edu.whut.info.core;

import edu.whut.info.dataset.Result;
import edu.whut.info.dataset.Segment;
import edu.whut.info.util.BioToolbox;

import java.util.*;
import java.util.logging.Logger;

/**
 * Created by Liu on 2016/3/21.
 */
public class FastPCF implements SegmentCutter {
    private final boolean Show_Debug = false;
    private int winK;
    private double penaltyR;
    private Logger m_log;
    private int MIN_SEG_LENGTH = 256;
    private double percent;
    private String methodName = "FastPCF";
    private double[] chromosome;

    public FastPCF(int winK, double percent, double penaltyR) {
        this.winK = winK;
        this.percent = percent;
        this.penaltyR = penaltyR;
        m_log = Logger.getLogger("segment");
    }

    public void setPenaltyR(double penaltyR) {
        this.penaltyR = penaltyR;
    }

    @Override
    public void splitChromosome(double[] data, Set<Segment> result, short chrId) {
        prepareCopyNumberSegment(data);
        Segment input = new Segment();
        input.setChr_id(chrId);
        input.Seg_id = 0;
        input.setRange(0, chromosome.length);

        //m_log.info(String.format("segment by %s", methodName));
        ArrayList<Integer> candidates = highPassFiltering(chromosome, winK, percent);
        List<Long> minIndex = new ArrayList<>();
        ArrayList<Double> arrAk = new ArrayList<>();
        ArrayList<Double> arrCk = new ArrayList<>();
        ArrayList<Double> arrDk = new ArrayList<>();
        ArrayList<Double> arrEk = new ArrayList<>();
        ArrayList<Double> uk = new ArrayList<>(candidates.size());
        uk.add(0.0);

        int prior, next;
        prior = candidates.get(0);
        for (int i = 1; i < candidates.size(); i++) {
            next = candidates.get(i);
            uk.add(sum(chromosome, prior, next));
            prior = next;
        }

        arrEk.add(0.0);
        minIndex.add(0L);
        next = candidates.get(0);
        for (int i = 1; i < uk.size(); i++) {
            plus(arrAk, uk.get(i));
            prior = next;
            next = candidates.get(i);

            plus(arrCk, next - prior);
            multiply(arrAk, arrCk, arrDk);
            minDk(arrDk, arrEk, minIndex);
        }

        //output.add(input.getSubSegment(Candidates.get(start).intValue(), Length));
        candidates.add(chromosome.length);
        minIndex.add((long) (candidates.size() - 1));

        int start = minIndex.get(minIndex.size() - 2).intValue();
        for (int stop = minIndex.size() - 1; stop > 0; ) {

            int temp = stop - 1;
            while (candidates.get(stop) - candidates.get(start) <= MIN_SEG_LENGTH) {
                start = minIndex.get(temp).intValue();
                temp = (start < temp) ? start : temp - 1;

                if (candidates.get(start) < MIN_SEG_LENGTH) {
                    start = 0;
                    break;
                }
            }

            result.add(input.getSubSegment(candidates.get(start).intValue(), candidates.get(stop).intValue()));
            stop = start;
            start = minIndex.get(stop).intValue();

        }
        int i = 1;
        for (Segment seg : result) {
            BioToolbox.refreshSegment(seg,data);

            if (Show_Debug) m_log.info(seg.getCharacterString());
            // m_log.info(String.format("the %2d segment:\t start=%6d\t end=%6d\t mean=%.4f\t std=%.4f", i, seg.range.Start,seg.range.End,seg.HalfCopyNumber,seg.stdHalfCopyNumber));
            i++;
        }
    }

//    public void refreshSegment(Segment seg) {
//        if (seg.isDirty) {
//            double[] ms;
//            ms = BioToolbox.robustMean(chromosome, seg.Start(), seg.End(), 0);
//            seg.CopyNumber = ms[0];
//            seg.stdCopyNumber = ms[1];
//            seg.isDirty = false;
//        }
//    }

    @Override
    public String getMethodName() {
        return methodName;
    }

    @Override
    public void prepareCopyNumberSegment(double[] data) {
        chromosome = new double[data.length];
        for (int i = 0; i < data.length; i++) {
            chromosome[i] = BioToolbox.log2(data[i]);
        }


    }

    private double sum(double[] arr, int a, int b) {
        double sum = 0.;

        for (int i = a; i < b; i++) {
            sum += arr[i];
        }
        return sum;
    }

    public ArrayList<Long> highPassFiltering(Collection<Double> input, int width, double percent) {
        int length = input.size();
        double[] values = new double[length];

        int i = 0;
        for (double a : input) {
            values[i++] = a;
        }

        double[] diff = BioToolbox.difference(values, width);
        double[] sortDiff = new double[length];
        for (i = 0; i < length; i++) {
            values[i] = Math.abs(diff[i]);
            sortDiff[i] = values[i];
        }
        Arrays.parallelSort(sortDiff);

        double threshold = sortDiff[(int) Math.round(sortDiff.length * percent + 0.5)];

        ArrayList<Long> index = new ArrayList<>();//potential breakpoints
        index.add(0L);
        for (i = width >> 1; i < values.length - (width >> 1); i++) {
            if (values[i] > threshold)
                index.add((long) i);
        }
        index.add((long) length);
        return index;
    }

    private void plus(ArrayList<Double> vector, double a) {
        for (int i = 0; i < vector.size(); i++) {
            vector.set(i, vector.get(i) + a);
        }
        vector.add(a);
    }

    private void multiply(ArrayList<Double> arrA, ArrayList<Double> arrCk, ArrayList<Double> arrDk) {
        int length = arrA.size();
        arrDk.add(0.0);
        for (int i = 0; i < length; i++) {
            arrDk.set(i, -1 * arrA.get(i) * arrA.get(i) / arrCk.get(i));
        }
    }

    private void minDk(ArrayList<Double> arrDk, ArrayList<Double> arrE, List<Long> minIndex) {
        double min = Integer.MAX_VALUE;
        long index = 0;

        int length = arrDk.size() - 1;
        for (int i = 0; i < length; i++) {
            double temp = arrDk.get(i) + arrE.get(i) + penaltyR;
            if (temp < min) {
                min = temp;
                index = i;
            }
        }
        arrE.add(min);
        minIndex.add(index);
    }

    public ArrayList<Integer> highPassFiltering(double[] data, int width, double percent) {
        int length = data.length;
        double[] values = new double[length];
        for (int i = 0; i < length; i++)
            values[i] = data[i];
        double[] diff = BioToolbox.difference(values, width);
        double[] sortDiff = new double[length];
        for (int i = 0; i < length; i++) {
            values[i] = Math.abs(diff[i]);
            sortDiff[i] = values[i];
        }
        Arrays.parallelSort(sortDiff);

        double threshold = sortDiff[(int) Math.round(sortDiff.length * percent + 0.5)];

        ArrayList<Integer> index = new ArrayList<>();//potential breakpoints
        index.add(0);
        for (int i = width >> 1; i < values.length - (width >> 1); i++) {
            if (values[i] > threshold)
                index.add(i);
        }
        index.add(length);

        return index;
    }

    @Override
    public List<Result> getResult() {
        return null;
    }
}
