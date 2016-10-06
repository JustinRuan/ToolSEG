package edu.whut.info.core;

import edu.whut.info.dataset.Segment;
import edu.whut.info.util.BioToolbox;
import smile.regression.LASSO;

import java.util.*;
import java.util.logging.Logger;

/**
 * Created by Liu on 2016/3/21.
 */
public class LassoSegment implements SegmentCutter {
    private String methodName = "LassoSegment";
    private Logger m_log;
    private double lamda;
    private double tol;
    private int maxIter;
    private double[] chromosome;

    public LassoSegment(double lamda, double tol, int maxIter) {
        this.lamda = lamda;
        this.tol = tol;
        this.maxIter = maxIter;
        m_log = Logger.getLogger("segment");
    }

    @Override
    public void splitChromosome(double[] data, Set<Segment> result, short chrId) {
        m_log.info(String.format("segment by %s", methodName));
        //chromosome=data;
        prepareCopyNumberSegment(data);
        Segment input = new Segment();
        input.setChr_id(chrId);
        input.Seg_id = 0;
        double[] arr = BioToolbox.mean_std(data);
        input.CopyNumber = arr[0];
        input.stdCopyNumber = arr[1];
        input.setRange(0, data.length);
        Queue<Segment> splittingSegments = new LinkedList<Segment>();
        int totLength = data.length;//Ⱦɫ���ܳ�
        int K = totLength / 2000;
        int temp = 0;
        ArrayList<Long> potentials = new ArrayList<>();
        for (int k = 0; k < K; k++) {
            Segment seg = input.getSubSegment(temp, temp + 2000);
            splittingSegments.add(seg);
            temp += 2000;
        }
        if (temp != totLength)
            splittingSegments.add(input.getSubSegment(temp, totLength));
        while (!splittingSegments.isEmpty()) {
            Segment curSeg = splittingSegments.poll();
            ArrayList<Double> values = getsubdata(data, curSeg.Start(), curSeg.End());
            double[] y = new double[values.size()];
            int n = 0;
            for (double a : values) {
                y[n] = a;
                n++;
            }
            //get possible change-points by Lasso

            double[][] x = new double[y.length][y.length];
            for (int i = 0; i < y.length; i++) {
                for (int j = 0; j <= i; j++) {
                    x[i][j] = 1.;
                }
            }
            LASSO lasso = new LASSO(x, y, lamda, tol, maxIter);
            double[] coffienent = lasso.coefficients();
            ArrayList<Integer> Kmax = new ArrayList<>(50);//suppose the max number of change-point

            for (int i = 1; i < coffienent.length; i++) {
                if (Math.abs(coffienent[i]) > 0.01) {
                    Kmax.add(i);
                    //Kmax.add(curSeg.Index.get(i));
                    if (Kmax.size() == 50)
                        break;
                }
            }
            // for the case of number of change-point is less than 20
            for (int i = 0; i < Kmax.size(); i++) {
                if (Kmax.get(i) == 0) {
                    Kmax.remove(i);
                }
            }
            int start = 0;
            for (int stop : Kmax) {
                result.add(input.getSubSegment(curSeg.Start() + start, curSeg.Start() + stop));
                start = stop;
            }
            if (Kmax.size() == 0)
                result.add(input.getSubSegment(curSeg.Start(), curSeg.End()));
            else
                result.add(input.getSubSegment(curSeg.Start() + start, curSeg.End()));

        }

        //  mergeSegment(result,0.3);
        int i = 1;
        for (Segment seg : result) {
            refreshSegment(seg);
            m_log.info(seg.getCharacterString());
            //  m_log.info(String.format("the %2d segment:\t start=%6d\t end=%6d\t mean=%.4f\t std=%.4f", i, seg.range.Start,seg.range.End,seg.HalfCopyNumber,seg.stdHalfCopyNumber));
            i++;
        }
    }

    private void mergeSegment(Set<Segment> output, double threshold) {
        Set<Segment> temp = new TreeSet<Segment>();
        Iterator<Segment> it = output.iterator();
        Segment preSegment = it.next();
        while (it.hasNext()) {
            Segment nextSegment = it.next();
            refreshSegment(preSegment);
            refreshSegment(nextSegment);
            if (Math.abs(preSegment.CopyNumber - nextSegment.CopyNumber) <= threshold || preSegment.length() < 20 || nextSegment.length() < 20) {
                preSegment.setRange(preSegment.Start(), nextSegment.End());
            } else {
                temp.add(preSegment);
                preSegment = nextSegment;
            }
        }
        temp.add(preSegment);
        output.clear();
        output.addAll(temp);
    }

    public void refreshSegment(Segment seg) {
        if (seg.isDirty) {
            double[] ms;
            ms = BioToolbox.robustMean(chromosome, seg.Start(), seg.End(), 16);
            seg.CopyNumber = ms[0];
            seg.stdCopyNumber = ms[1];

            seg.isDirty = false;
        }
    }

    public ArrayList<Double> getsubdata(double[] data, int start, int stop) {
        ArrayList<Double> result = new ArrayList<>();
        for (int i = start; i < stop; i++)
            result.add(data[i]);
        return result;
    }

    @Override
    public String getMethodName() {
        return methodName;
    }

    @Override
    public void prepareCopyNumberSegment(double[] data) {
        chromosome = new double[data.length];
        for (int i = 0; i < data.length; i++) {
            chromosome[i] = data[i];
            data[i] = BioToolbox.log2(data[i]);
        }

    }
}
