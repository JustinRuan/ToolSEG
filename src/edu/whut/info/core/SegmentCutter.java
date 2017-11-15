package edu.whut.info.core;

import edu.whut.info.dataset.Result;
import edu.whut.info.dataset.Segment;

import java.util.List;
import java.util.Set;

/**
 * @author Liu
 * @version 1.0
 * @created 11-ʮ��-2015 16:32:15
 */
public interface SegmentCutter {
    void enableShowDebug(boolean value);
    void splitChromosome(double[] data, Set<Segment> result, short chrId);
    String getMethodName();
    void prepareCopyNumberSegment(double[] data);
    List<Result> getResult();

}