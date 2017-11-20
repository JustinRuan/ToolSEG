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
package edu.whut.info.dataset;


import edu.whut.info.util.GeneRange;

import java.io.Serializable;

/**
 * @author Justin
 * @version 2.0
 * @created 12-03-2016 11:58:21
 */
public class Segment implements Comparable<Segment>, Serializable {
    private static final long serialVersionUID = 4905497107751450036L;
    private static int IDCounter = 0;
    public int Seg_id;
    public short Chr_id;

    public GeneRange range;
    public boolean isReady;
    public boolean isDirty;

    public int getStart(){
        return  range.Start;
    }

    public double CopyNumber;
    public double stdCopyNumber;

    //public double purityCorrCopyNumber;

    public Segment() {
        isReady = false;
        isDirty = true;
        range = new GeneRange();
        Seg_id = ++IDCounter ;
    }

    public static void resetID(){
        IDCounter = 0;
    }

    @Override
    public int compareTo(Segment o) {
        // TODO Auto-generated method stub
        int c = Short.compare(Chr_id, o.Chr_id);
        return (c == 0) ? Integer.compare(this.getStart(), o.getStart()) : c;
    }

    public void setChr_id(short id) {
        this.Chr_id = id;
        range.Chr_id = id;
    }

    public void setRange(int start, int end) {
        range.Start = start;
        range.End = end;
    }

    public int length() {
        return range.length();
    }

    public int Start() {
        return range.Start;
    }

    public int End() {
        return range.End;
    }

    public Segment getSubSegment(int start, int stop) {
        if (start == stop) {// for debug
            System.out.println("getSubSegment : Start == Stop!");
            System.exit(-1);
        }
        Segment newSeg = new Segment();
        //newSeg.Seg_id = start;
        newSeg.setChr_id(this.Chr_id);
        newSeg.setRange(start, stop);
        return newSeg;
    }

    public String getCharacterString() {
        return String.format("Seg_id = % 6d : %02d;[% 6d - % 6d]\t Length =% 6d;\t robustAvg = %.4f \t robustStd = %.4f",
                Seg_id, Chr_id, Start(), End(), length(), CopyNumber, stdCopyNumber);
    }

    public final boolean isIncluded(int pos) {
        return pos >= range.Start && pos < range.End;
    }

}


