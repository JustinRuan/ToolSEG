package edu.whut.info.core;

import edu.whut.info.dataset.Result;
import edu.whut.info.dataset.Segment;
import edu.whut.info.util.BioToolbox;

import java.util.ArrayList;
import java.util.List;
import java.util.Set;
import java.util.logging.Logger;

/**
 * Created by Liu on 2016/3/21.
 */
public class PCFM implements SegmentCutter {
    private final boolean Show_Debug = false;
    private String methodName = "PCF";
    private double penaltyR;
    private int MIN_SEG_LENGTH = 256;
    private Logger m_log;
    private double[] chromosome;

    public PCFM(double penaltyR) {
        this.penaltyR = penaltyR;
        m_log = Logger.getLogger("segment");
    }

    @Override

    public void splitChromosome(double[] data, Set<Segment> result, short chrId) {
        //m_log.info(String.format("Segment by%s: ", methodName));
        Segment input = new Segment();
        input.setChr_id(chrId);
        input.Seg_id = 0;
        double[] arr = BioToolbox.mean_std(data);
        input.CopyNumber = arr[0];
        input.stdCopyNumber = arr[1];
        input.setRange(0, data.length);
        chromosome = data.clone();
        int IndexLength = data.length;
        ArrayList<Double> VectorA = new ArrayList<>(IndexLength);
        ArrayList<double[]> VectorE = new ArrayList<>(IndexLength + 1);
        ArrayList<Double> VectorDk = new ArrayList<>(IndexLength);
        VectorE.add(new double[]{0., 0.});

        List<Long> minIndex = new ArrayList<>(IndexLength);//VecortE ??????????????
        //int offset = MIN_SEG_LENGTH / MIN_SEG_STEP;

        for (double a : chromosome) {

            plus(VectorA, a);
            multiply(VectorA, VectorDk);
            minDk(VectorDk, VectorE, minIndex);
        }

        int start = minIndex.get(minIndex.size() - 1).intValue();
        for (int stop = minIndex.size(); stop > 0; ) {

            int temp = stop - 1;
            while (stop - start <= MIN_SEG_LENGTH) {
                start = minIndex.get(temp).intValue();
                temp = (start < temp) ? start : temp - 1;

                if (start < MIN_SEG_LENGTH) {
                    start = 0;
                    break;
                }
            }

            result.add(input.getSubSegment(start, stop));
            stop = start;
            start = minIndex.get(stop).intValue();

        }

        StringBuilder sb = new StringBuilder();
        for (double[] E : VectorE) {
            sb.append(String.format("%f\t%f\n", E[0], E[1]));

        }
        m_log.info(sb.toString());

        int i = 1;
        for (Segment seg : result) {
            BioToolbox.refreshSegment(seg, data);
            if (Show_Debug) m_log.info(seg.getCharacterString());
            // m_log.info(String.format("the %2d segment:\t start=%6d\t end=%6d\t mean=%.4f\t std=%.4f", i, seg.range.Start,seg.range.End,seg.HalfCopyNumber,seg.stdHalfCopyNumber));
            i++;
        }
    }

    private void plus(ArrayList<Double> vector, double a) {
        for (int i = 0; i < vector.size(); i++) {
            vector.set(i, vector.get(i) + a);
        }
        vector.add(a);
    }

    private void multiply(ArrayList<Double> arrA, ArrayList<Double> arrDk) {
        int length = arrA.size();
        arrDk.add(0.0);
        for (int i = 0; i < length; i++) {
            arrDk.set(i, -1 * arrA.get(i) * arrA.get(i) / (length - i));
        }
    }

    private void minDk(ArrayList<Double> arrDk, ArrayList<double[]> arrE, List<Long> minIndex) {
        double min = Integer.MAX_VALUE;
        long index = 0;
        double[] minE = new double[2];

        int length = arrDk.size() - 1;
        for (int i = 0; i < length; i++) {
            double[] E = new double[2];
            E[0] = arrDk.get(i) + arrE.get(i)[0];
            E[1] = arrE.get(i)[1] + penaltyR;
            //double temp = arrDk.get(i) + arrE.get(i) + penaltyR;
            double temp = E[0] + E[1];
            if (temp < min) {
                min = temp;
                index = i;
                minE = E;
            }
        }
        arrE.add(minE);
        minIndex.add(index);
    }

    @Override
    public String getMethodName() {
        return methodName;
    }

    @Override
    public void prepareCopyNumberSegment(double[] data) {

    }

    @Override
    public List<Result> getResult() {
        return null;
    }
}
