package edu.whut.info.core;


import edu.whut.info.util.BioToolbox;

import java.util.List;

/**
 * Created by Liu on 2015/12/30.
 */
public class SegmentProcessor {
    public static void smoothing(List<Double> seg, int w) {
        double[] cnValues = new double[seg.size()];
        int i = 0;
        for (double value : seg) {
            cnValues[i++] = value;
        }

        cnValues = BioToolbox.GaussianBlur(cnValues, w, 1.0);
        i = 0;
        for (; i < seg.size(); i++) {
            seg.set(i, cnValues[i]);
        }
    }
}
