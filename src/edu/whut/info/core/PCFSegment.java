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

import java.util.ArrayList;
import java.util.List;
import java.util.Set;
import java.util.logging.Logger;

/**
 * Created by Liu on 2016/3/21.
 */
public class PCFSegment implements SegmentCutter {
    private boolean Show_Debug = false;
    private String methodName = "PCFSegment";
    private double penaltyR;
    private int MIN_SEG_LENGTH = 256;
    private Logger m_log;
    private double[] chromosome;

    public PCFSegment(double penaltyR) {
        this.penaltyR = penaltyR;
        m_log = Logger.getLogger("segment");
    }

    @Override
    public void enableShowDebug(boolean value) {
        Show_Debug = value;
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
        ArrayList<Double> VectorE = new ArrayList<>(IndexLength + 1);
        ArrayList<Double> VectorDk = new ArrayList<>(IndexLength);
        VectorE.add(0.0);

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

    private void minDk(ArrayList<Double> arrDk, ArrayList<Double> arrE, List<Long> minIndex) {
        double min = Integer.MAX_VALUE;
        long index = 0;
        //arrE.add(0.0);
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
