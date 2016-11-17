package edu.whut.info.analysis;

import edu.whut.info.dataset.Result;

import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;
import java.util.logging.Logger;

/**
 * Created by Justin on 2016/10/25.
 */
public class ResultAnalysis {
    private Logger m_log;
    private List<int[]> segRanges;
    private List<Boolean> includeBpt;
    private List<Double> values;

    public ResultAnalysis() {
        m_log = Logger.getLogger("segment");
        segRanges = new ArrayList<>();
        includeBpt = new ArrayList<>();
        values = new ArrayList<>();
    }

    private void clear(){
        segRanges.clear();
        includeBpt.clear();
        values.clear();
    }

    public void initialize(List<Long> breakPoints, int chrLength){
        final int Tolerance = 10;
        clear();

        int begin,end;
        begin = 0;
        end = 0;
        for (long t : breakPoints) {
            end = (int) t - 1 - Tolerance;
            segRanges.add(new int[]{begin, end});
            includeBpt.add(false);
            values.add(0.);

            begin = end;
            end = (int) t - 1 + Tolerance;
            segRanges.add(new int[]{begin,end});
            includeBpt.add(true);
            values.add(0.);

            begin = end;
        }

        begin = end;
        end = chrLength;
        segRanges.add(new int[]{begin, end});
        includeBpt.add(false);
        values.add(0.);

    }

    private final boolean isIncluded(int pos, int[] range){
        if (pos >= range[0] && pos < range[1]){
            return true;
        }else
            return false;
    }

    public void analysisResult(List<Result> results, int chrid) {
        int k = 0;
        int count = segRanges.size();
        for (Result r : results){
            for (int i = k; i<count; i++){
                if (isIncluded(r.pos,segRanges.get(i))){
                    double maxV = values.get(k);
                    maxV = r.value1 > maxV ? r.value1:maxV;
                    values.set(k,maxV);
                    break;
                }else{
                    k++;
                }
            }
        }

        for (int i = 0; i < segRanges.size(); i++){
            // chr id, range[S, E], isBreak, value
            int[] range = segRanges.get(i);
            m_log.info(String.format("\t% 4d\t% 4d\t%d\t%s\t%f",chrid,range[0],range[1],includeBpt.get(i).toString(),values.get(i)));
        }

    }
}
