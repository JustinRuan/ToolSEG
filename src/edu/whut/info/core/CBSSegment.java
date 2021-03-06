//Copyright 2017 Jun Ruan
//
//        Licensed under the Apache License, Version 2.0 (the "License");
//        you may not use this file except in compliance with the License.
//        You may obtain a copy of the License at
//
//        http://www.apache.org/licenses/LICENSE-2.0
//
//        Unless required by applicable law or agreed to in writing, software
//        distributed under the License is distributed on an "AS IS" BASIS,
//        WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//        See the License for the specific language governing permissions and
//        limitations under the License.
package edu.whut.info.core;

import edu.whut.info.dataset.Result;
import edu.whut.info.dataset.Segment;
import edu.whut.info.util.BioToolbox;

import java.util.*;
import java.util.logging.Logger;

/**
 * Created by Liu on 2016/3/21.
 */
public class CBSSegment implements SegmentCutter {
    private boolean Show_Debug = true;
    private int MIN_SEG_LENGTH;
    private int MIN_SEG_STEP = 2;
    private String name = "CBSegmentCutter";
    private Logger m_log;
    private double[] chromosome;
    private double[] IntegralCN;
    private List<Result> results;

    public CBSSegment(int MIN_SEG_LENGTH, int MIN_SEG_STEP) {
        this.MIN_SEG_LENGTH = MIN_SEG_LENGTH;
        this.MIN_SEG_STEP = MIN_SEG_STEP;
        m_log = Logger.getLogger("segment");
        results=new LinkedList<>();
    }

    @Override
    public void enableShowDebug(boolean value) {
        Show_Debug = value;
    }

    @Override
    public void splitChromosome(double[] data, Set<Segment> result, short chrId) {
        //m_log.info(String.format("segment by %s", name));
        //chromosome=data;
        prepareCopyNumberSegment(data);
        Segment input = new Segment();
        input.setChr_id(chrId);
        input.Seg_id = 0;
        double[] arr = BioToolbox.mean_std(chromosome);
        input.CopyNumber = arr[0];
        input.stdCopyNumber = arr[1];
        input.setRange(0, data.length);
        Queue<Segment> splittingSegments = new LinkedList<Segment>();
        //splittingSegments.add(input);
        int totLength = data.length;//
        final int L = 2000;//100000;//2000
        int K = totLength / L;
        int temp = 0;
        ArrayList<Long> potentials = new ArrayList<>();
        for (int k = 0; k < K; k++) {
            Segment seg = input.getSubSegment(temp, temp + L);
            splittingSegments.add(seg);
            temp += L;
        }
        if (temp != totLength)
            splittingSegments.add(input.getSubSegment(temp, totLength));
        List<Integer> Index = new ArrayList<Integer>();
        while (!splittingSegments.isEmpty()) {
            Segment curSeg = splittingSegments.poll();
            //refreshSegment(curSeg,chromosome);
            BioToolbox.refreshSegment(curSeg, chromosome);

            int segLength = curSeg.length();
            double maxZ = 0.;
            int Start = 0, Stop = 0;
            // getThresh(curSeg);
            if (curSeg.End() - curSeg.Start() < 2 * MIN_SEG_LENGTH) {
                result.add(curSeg);
            } else {
                maxZ = 0.0;
                Index.clear();
                for (int i = curSeg.Start(); i < curSeg.End(); i += MIN_SEG_STEP) {
                    Index.add(i);
                }
                Index.set(Index.size() - 1, curSeg.End());

                int indexLength = Index.size();
                double zij;
                int offset = MIN_SEG_LENGTH / MIN_SEG_STEP;
                for (int i = 0; i < indexLength - offset; i++) {
                    for (int j = i + offset; j < indexLength - offset; j++) {
                        //if ((i<MIN_SEG_LENGTH/MIN_SEG_STEP) || (j>indexLength - MIN_SEG_LENGTH/MIN_SEG_STEP))
                        if ((Index.get(j) - Index.get(i) + 1) >= segLength - MIN_SEG_LENGTH)
                            continue;
                        zij = calculateZ1(curSeg, Index.get(i), Index.get(j));
                        if (zij > maxZ) {
                            maxZ = zij;
                            Start = Index.get(i);
                            Stop = Index.get(j);
                        }
                    }
                }
                // getThresh(curSeg);
                double splittingThreshold = getPValue(curSeg, maxZ);
                // if (maxZ > curSeg.splittingThreshold)
                if (splittingThreshold <= 0.1) {//splitting
                    List<Segment> newSegs = new LinkedList<Segment>();
                    if ((Stop < Start) || (Stop == Start)) {//debug info
                        System.out.println("splitChromosome Error!");
                        System.exit(-1);
                    }
                    splitSegment(curSeg, newSegs, Start, Stop);
                    if (newSegs.size() > 0) {
                        splittingSegments.addAll(newSegs);
                    }
                } else {
                    //add to result
                    curSeg.isReady = true;
                    result.add(curSeg);
                }
            }
        }
        mergeSegment(result, 0.2);
        int i = 1;
        for (Segment seg : result) {
            seg.isDirty = true;
            BioToolbox.refreshSegment(seg,data);
            if (Show_Debug) m_log.info(seg.getCharacterString());
            // m_log.info(String.format("the %2d segment:\t start=%6d\t end=%6d\t mean=%.4f\t std=%.4f", i, seg.range.Start,seg.range.End,seg.HalfCopyNumber,seg.stdHalfCopyNumber));
            i++;
        }
    }

    private void mergeSegment(Set<Segment> output, double threshold) {
        Set<Segment> temp = new TreeSet<Segment>();
        Iterator<Segment> it = output.iterator();
        Segment preSegment = it.next();
        while (it.hasNext()) {
            Segment nextSegment = it.next();
//            refreshSegment(preSegment,chromosome);
//            refreshSegment(nextSegment,chromosome);
            BioToolbox.refreshSegment(preSegment, chromosome);
            BioToolbox.refreshSegment(preSegment, chromosome);
            if (Math.abs(preSegment.CopyNumber - nextSegment.CopyNumber) <= threshold) {
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

    private void splitSegment(Segment input, List<Segment> output, int start, int stop) {

        int s1 = (start - input.Start() > MIN_SEG_LENGTH) ? 1 : 0;
        int s2 = (input.End() - stop > MIN_SEG_LENGTH) ? 1 : 0;
        int s = (s1 << 1) | s2;

        switch (s) {
            case 0x00:
            case 0x03:
                output.add(input.getSubSegment(input.Start(), start));
                output.add(input.getSubSegment(start, stop));
                output.add(input.getSubSegment(stop, input.End()));
                break;
            case 0x02:
                //cut the first part to new segment, if the third part is too short
                output.add(input.getSubSegment(input.Start(), start));
                output.add(input.getSubSegment(start, input.End()));
                break;
            case 0x01:
                //cut the third part to new segment, if the first part is too short
                output.add(input.getSubSegment(input.Start(), stop));
                output.add(input.getSubSegment(stop, input.End()));
                break;
//                case 0x00:
//                    // if the first and third part all are too short.
//                    // do nothing
//                    break;
            default:

        }
    }

    public double calculateZ1(Segment seg, int start, int stop) {
        if (stop < start) {
            int temp = start;
            start = stop;
            stop = temp;
        }
        int len = stop - start;
        int total = seg.length();
        double sum1 = 0, sum2 = 0;
        sum1 = IntegralCN[stop] - IntegralCN[start];
        sum2 = IntegralCN[total] - sum1;
        double a = Math.sqrt(1. / len + 1. / (total - len));
        double b = sum1 / len - sum2 / (total - len);
        return 1 / Math.sqrt(1. / len + 1. / (total - len)) * Math.abs(sum1 / len - sum2 / (total - len));
    }

    public void refreshSegment(Segment seg,double[] data) {
        if (seg.isDirty) {
            double[] ms;
            ms = BioToolbox.robustMean(data, seg.Start(), seg.End(), 16);
            seg.CopyNumber = ms[0];
            seg.stdCopyNumber = ms[1];

            seg.isDirty = false;
        }
    }

    public double getPValue(Segment seg, double maxZ) {
        ArrayList<Double> values = getsubdata(chromosome, seg.Start(), seg.End());
        int count = 0;
        for (int n = 0; n < 100; n++) {
            Collections.shuffle(values, new Random(1));
            List<Integer> Index = new ArrayList<>();
            double max = 0;
            int segLength = seg.length();
            for (int i = 0; i < segLength; i += MIN_SEG_STEP) {
                Index.add(i);
            }
            Index.set(Index.size() - 1, segLength - 1);
            int indexLength = Index.size();
            double zij;
            List<Double> IntegralCN = BioToolbox.resetIntegralCN(values);
            int offset = MIN_SEG_LENGTH / MIN_SEG_STEP;
            for (int i = 0; i < indexLength - offset; i++) {
                for (int j = i + offset; j < indexLength - offset; j++) {
                    if ((Index.get(j) - Index.get(i) + 1) >= segLength - MIN_SEG_LENGTH)
                        continue;
                    zij = cal(IntegralCN, Index.get(i), Index.get(j));
                    if (zij > max) {
                        max = zij;
                    }
                }
            }
            if (max >= maxZ)
                count++;
        }
        double splittingThreshold = count * 1. / 100;
        return splittingThreshold;
    }

    public ArrayList<Double> getsubdata(double[] data, int start, int stop) {
        ArrayList<Double> result = new ArrayList<>();
        for (int i = start; i < stop; i++)
            result.add(data[i]);
        return result;
    }

    public double cal(List<Double> IntegralCN, int start, int stop) {
        if (stop < start) {
            int temp = start;
            start = stop;
            stop = temp;
        }
        int len = stop - start;
        int total = IntegralCN.size() - 1;
        double sum1 = 0, sum2 = 0;
        sum1 = IntegralCN.get(stop) - IntegralCN.get(start);
        sum2 = IntegralCN.get(total) - sum1;
        double a = Math.sqrt(1. / len + 1. / (total - len));
        double b = sum1 / len - sum2 / (total - len);
        return 1 / Math.sqrt(1. / len + 1. / (total - len)) * Math.abs(sum1 / len - sum2 / (total - len));
    }

    @Override
    public String getMethodName() {
        return name;
    }

    @Override
    public void prepareCopyNumberSegment(double[] data) {
        chromosome = new double[data.length];
        int length = data.length;
        IntegralCN = new double[length + 1];
        IntegralCN[0] = 0.0;
        double sum = 0;
        for (int j = 0; j < length; j++) {
            chromosome[j] = data[j];
            //chromosome[j] = BioToolbox.log2(data[j]);
            sum += chromosome[j];
            IntegralCN[j + 1] = sum;
        }
    }

    @Override
    public List<Result> getResult() {
        return null;
    }
}
