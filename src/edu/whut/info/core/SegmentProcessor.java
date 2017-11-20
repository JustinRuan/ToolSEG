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
