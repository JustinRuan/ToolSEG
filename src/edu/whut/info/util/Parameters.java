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
package edu.whut.info.util;

import java.util.ArrayList;

/**
 * @author Justin
 * @version 1.0
 * @created 08-11 -2014 22:09:01
 */
public class Parameters {

    public ArrayList<GeneRange> roiArray;

    public int smoothingCopyNumber_Width;
    public int rou_WindowsSize;
    public double PVALUE_THRESH;
    public int MIN_SEG_LENGTH;
    public int MIN_SEG_STEP;
    public double COMBINE_FACTOR;
    public double ALPHA_POSITION;

    public boolean enable_rou_chart;
    public boolean enable_copynumber_histogram;
    public boolean enable_probeset_chart;

    public boolean enableCorrectCrossTalk;
    public boolean enableScanRou;
    public boolean OnlySNP;
    public int MinimumCountForBaseLine;

    public double widthFourthPart;
    public double widthNormalPart;
    public int limitAttenuationM;

    public Parameters() {
        roiArray = new ArrayList<GeneRange>();

        smoothingCopyNumber_Width = 128;
        rou_WindowsSize = 256;
        PVALUE_THRESH = 0.05;
        MIN_SEG_LENGTH = 256;
        MIN_SEG_STEP = 8;
        COMBINE_FACTOR = 0.2;


        OnlySNP = true;
        enable_copynumber_histogram = true;
        enableCorrectCrossTalk = false;
        enableScanRou = false;
        MinimumCountForBaseLine = 500;
        widthFourthPart = 0.1;
        widthNormalPart = 0.15;
        limitAttenuationM = 6;
    }


}