package edu.whut.info.dataset;


import edu.whut.info.util.GeneRange;

import java.io.Serializable;

/**
 * @author Justin
 * @version 2.0
 * @created 12-03月-2016 11:58:21
 */
public class Segment implements Comparable<Segment>, Serializable {
    private static final long serialVersionUID = 4905497107751450036L;
    public static int IDCounter = 0;
    public int Seg_id;
    public short Chr_id;

    public GeneRange range;
    public boolean isReady;
    public boolean isDirty;
    public double HalfCopyNumber;
    public double stdHalfCopyNumber;
    public double purityCorrCopyNumber;

    public Segment() {
        isReady = false;
        isDirty = true;
        range = new GeneRange();
        Seg_id = ++IDCounter;
    }

    public int getStart() {
        return range.Start;
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
        if (start == stop) {//调试用
            System.out.println("getSubSegment : Start == Stop!");
            System.exit(-1);
        }
        Segment newSeg = new Segment();
        //newSeg.Seg_id = start;
        newSeg.setChr_id(this.Chr_id);
        newSeg.setRange(start, stop);
        return newSeg;
    }

    //    public Segment getSubSegment1(int start, int stop,double maxZ) {
//        if (start == stop) {//调试用
//            System.out.println("getSubSegment : Start == Stop!");
//            System.exit(-1);
//        }
//        Segment newSeg = new Segment();
//        newSeg.Seg_id = start;
//        newSeg.setChr_id(this.Chr_id);
//        newSeg.setRange(start, stop);
//        newSeg.maxZ=maxZ;
//        return newSeg;
//    }
    public String getCharacterString() {
        return String.format("Seg_id = % 6d : %02d;[% 6d - % 6d]\t Length =% 6d;\t robustAvg = %.4f \t robustStd = %.4f",
                Seg_id, Chr_id, Start(), End(), length(), HalfCopyNumber, stdHalfCopyNumber);
    }

}


