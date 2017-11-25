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

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.*;
import java.util.logging.Logger;

public class Import_Results implements SegmentCutter {
    private boolean Show_Debug = false;
    private String methodName = "Import_Results";
    private Logger m_log;
    private Map<Short, Set<Segment>> Results;

    public Import_Results() {
        m_log = Logger.getLogger("segment");
        try {
            importResultFromFile();
        } catch (Exception e) {
            System.out.println(e.getMessage());
        }

    }

    private void importResultFromFile() throws Exception {
        Results = new HashMap<>();

        List<String> filenames = new ArrayList<>();
//        filenames.add(".//data//LTestdata_20161026_151708_result_CBS.csv");
//        filenames.add(".//data//LTestdata_20161026_152635_result_CBS.csv");
//        filenames.add(".//data//LTestdata_20161026_153427_result_CBS.csv");
        filenames.add(".//data//N_Testdata_20161201_161019_result_CBS.csv");
        filenames.add(".//data//N_Testdata_20161201_164912_result_CBS.csv");
        filenames.add(".//data//R_Testdata_20161026_151627_result_CBS.csv");
        filenames.add(".//data//R_Testdata_20161026_152405_result_CBS.csv");
//        filenames.add(".//data//R_Testdata_20161026_153220_result_CBS.csv");
        filenames.add(".//data//R2_Testdata_20170219_122916_result_CBS.csv");
//        filenames.add(".//data//STestdata_20161026_152515_result_CBS.csv");
//        filenames.add(".//data//STestdata_20161026_153258_result_CBS.csv");
//        filenames.add(".//data//N_Testdata_20161201_161019_result_pcf.csv");
//        filenames.add(".//data//N_Testdata_20161201_164912_result_pcf.csv");
//        filenames.add(".//data//R_Testdata_20161026_151627_result_pcf.csv");
//        filenames.add(".//data//R_Testdata_20161026_152405_result_pcf.csv");
////        filenames.add(".//data//R_Testdata_20161026_153220_result_pcf.csv");
//        filenames.add(".//data//R2_Testdata_20170219_122916_result_pcf.csv");
        short chrId = 1;
        short global_chrId = 1;

        for (String filepath : filenames) {
            BufferedReader in = new BufferedReader(new FileReader(filepath));
            in.readLine();
            int lastPos = 0;
            String s = "";
            Set<Segment> oneChromo = new TreeSet<Segment>();
            while ((s = in.readLine()) != null) {
                String[] arr = s.split(",");
                if (arr.length == 4) {
                    short id = Short.valueOf(arr[0]);

                    Segment seg = new Segment();
                    seg.Chr_id = id;
                    int start = Integer.valueOf(arr[1]) - lastPos;
                    int end = Integer.valueOf(arr[2]) - lastPos + 1;
                    seg.setRange(start, end);
                    seg.CopyNumber = Double.valueOf(arr[3]);
                    seg.isDirty = false;

                    if ((chrId == id) || oneChromo.isEmpty()) {
                        oneChromo.add(seg);
                        chrId = id;
                    } else {
                        Results.put(global_chrId, oneChromo);
                        global_chrId++;
                        oneChromo = new TreeSet<Segment>();
                        chrId = id;
                        lastPos = Integer.valueOf(arr[1]) - 1;
                        start = Integer.valueOf(arr[1]) - lastPos;
                        end = Integer.valueOf(arr[2]) - lastPos + 1;
                        seg.setRange(start, end);

                        oneChromo.add(seg);
                    }
                } else {
                    throw new Exception("The format of import file is wrong!");
                }

            }
            Results.put(global_chrId, oneChromo);
//            chrId++;
            global_chrId++;
        }
    }


    @Override
    public void enableShowDebug(boolean value) {
        Show_Debug = value;
    }

    @Override
    public void splitChromosome(double[] data, Set<Segment> result, short chrId) {
        Set<Segment> temp = Results.get(chrId);
        if (temp.size() == 0) {
            System.out.println("");
        }
        result.addAll(Results.get(chrId));
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
