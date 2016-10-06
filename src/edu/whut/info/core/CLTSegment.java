package edu.whut.info.core;

import edu.whut.info.dataset.Segment;
import edu.whut.info.util.BioToolbox;
import org.apache.commons.math3.distribution.NormalDistribution;

import java.util.*;
import java.util.logging.Logger;

/**
 * Created by Liu on 2016/3/19.
 */
public class CLTSegment implements SegmentCutter {
    private double PVALUE_THRESH = 0.05;
    private int MIN_SEG_LENGTH = 256;
    private int MIN_SEG_STEP = 8;
    // private double COMBINE_FACTOR = 0.2;
    private Logger m_log;
    private double[] IntegralCN;
    private double[] chromosome;
    private List<Double> robustSTD;
    private String methodname = "CLTSegment";

    public CLTSegment(double pvalueThre, int minSegLen, int minStep) {
        m_log = Logger.getLogger("segment");
        this.PVALUE_THRESH = pvalueThre;
        this.MIN_SEG_LENGTH = minSegLen;
        this.MIN_SEG_STEP = minStep;
        // this.COMBINE_FACTOR = combinefactor;
        this.robustSTD = new ArrayList<>();
    }

    @Override
    public void splitChromosome(double[] data, Set<Segment> result, short Chr_id) {
        m_log.info(String.format("segment by %s", methodname));
        //chromosome=data;
        prepareCopyNumberSegment(data);
        Segment input = new Segment();
        input.setChr_id(Chr_id);
        input.Seg_id = 0;
        double[] arr = BioToolbox.mean_std(data);
        input.CopyNumber = arr[0];
        input.stdCopyNumber = arr[1];
        input.setRange(0, data.length);
        splitChromosome(input, result);
        m_log.info(String.format("#### #### #### chr %02d First Step: Loci Count = %05d; \t Segments Count = %d",
                Chr_id, data.length, result.size()));
    }

    @Override
    public String getMethodName() {
        return methodname;
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
            data[j] = BioToolbox.log2(data[j]);
            sum += data[j];
            IntegralCN[j + 1] = sum;
        }
    }

    private void splitChromosome(Segment input, Set<Segment> output) {

        if (chromosome.length == 0) return;
        double[] arr = {input.CopyNumber, input.stdCopyNumber};
        double maxZ = 0;
        int Start = 0;
        int Stop = 0;

        Queue<Segment> splittingSegments = new LinkedList<Segment>();
        Set<Segment> scrapSegments = new LinkedHashSet<Segment>();

        splittingSegments.add(input);
        List<Integer> Index = new ArrayList<Integer>();
        int newSegCount = 0;

        while (!splittingSegments.isEmpty()) {

            Segment curSeg = splittingSegments.poll();
            int segLength = curSeg.length();

            refreshSegment(curSeg);

            //?���̫��
            if (segLength <= 2 * MIN_SEG_LENGTH) {
                scrapSegments.add(curSeg);
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
                //int offset = 1;
                for (int i = 0; i < indexLength - offset; i++) {
                    for (int j = i + offset; j < indexLength - offset; j++) {
                        //������ڽӽ�������
                        //if ((i<MIN_SEG_LENGTH/MIN_SEG_STEP) || (j>indexLength - MIN_SEG_LENGTH/MIN_SEG_STEP))
                        if ((Index.get(j) - Index.get(i) + 1) >= segLength - MIN_SEG_LENGTH)
                            continue;

                        zij = calculateZ(curSeg, Index.get(i), Index.get(j));
                        if (zij > maxZ) {
                            maxZ = zij;
                            Start = Index.get(i);
                            Stop = Index.get(j);
                        }
                    }
                }
                double splittingThreshold = getThreshold(curSeg, PVALUE_THRESH);
                if (maxZ > splittingThreshold) {
                    //splitting
                    List<Segment> newSegs = new LinkedList<Segment>();

                    if ((Stop < Start) || (Stop == Start)) {//debug info
                        System.out.println("splitChromosome Error!");
                        System.exit(-1);
                    }

                    splitSegment(curSeg, newSegs, Start, Stop);
                    splittingSegments.addAll(newSegs);

//                        //debug info
//                        for (Segment seg : splittingSegments){
//                            System.out.print(seg.range);
//                            System.out.print('\t');
//                        }
//                        System.out.print('\n');

                } else {
                    //add to result
                    curSeg.isReady = true;
                    output.add(curSeg);

                    newSegCount++;
                }
            }

            if (splittingSegments.isEmpty()) {
                if (newSegCount > 0) {
                    splittingSegments.addAll(mergeScrapSegments(scrapSegments, output));
                    scrapSegments.clear();
                    newSegCount = 0;
                } else {
                    output.addAll(scrapSegments);
                }
            }
        }
        mergeSegment(output, 0.3);
        int i = 1;
        for (Segment seg : output) {
            refreshSegment(seg);
            m_log.info(seg.getCharacterString());
            //  m_log.info(String.format("the %2d segment:\t start=%6d\t start=%6d\t mean=%.4f\t std=%.4f", i, seg.range.Start,seg.range.End,seg.HalfCopyNumber,seg.stdHalfCopyNumber));
            i++;
        }

    }

    private void mergeSegment(Set<Segment> output, double threshold) {
        Set<Segment> temp = new TreeSet<Segment>();
        Iterator<Segment> it = output.iterator();
        Segment preSegment = it.next();
        while (it.hasNext()) {
            Segment nextSegment = it.next();
            refreshSegment(preSegment);
            refreshSegment(nextSegment);
            if (Math.abs(preSegment.CopyNumber - nextSegment.CopyNumber) <= threshold || preSegment.length() < 20 || nextSegment.length() < 20) {
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

    public double calculateZ(Segment seg, int start, int stop) {

        int len = stop - start;

        double sum = 0;

        //��ͷ�Χ��Start��Stop - 1��������Stop
        sum = IntegralCN[stop] - IntegralCN[start];

        ///return Math.abs((sum - len*avgHalfCopyNumber)/(stdHalfCopyNumber*Math.sqrt(len)))
        double cn = seg.CopyNumber;
        double std = seg.stdCopyNumber;

        return Math.abs((sum - len * cn) / (std * Math.sqrt(len)));
    }

    public void refreshSegment(Segment seg) {
        if (seg.isDirty) {
            double[] ms;
            ms = BioToolbox.robustMean(chromosome, seg.Start(), seg.End(), 16);
            seg.CopyNumber = ms[0];
            seg.stdCopyNumber = ms[1];

            seg.isDirty = false;
        }
    }

    private double calculateZ1(Segment seg, int start, int stop, double[] arr) {

        int len = stop - start;

        double sum = 0;

        //��ͷ�Χ��Start��Stop - 1��������Stop
        sum = IntegralCN[stop] - IntegralCN[start];

        ///return Math.abs((sum - len*avgHalfCopyNumber)/(stdHalfCopyNumber*Math.sqrt(len)))
        double cn = arr[0];
        double std = arr[1];

        return Math.abs((sum - len * cn) / (std * Math.sqrt(len)));
    }

    private double getThreshold(Segment seg, double pValue) {
        NormalDistribution norm = new NormalDistribution(0, 1);

        // threshold of the z-score given the p-value cutoff
        // cn.length * (cn.length - 1) is used for correcting multiple
        // testing
        double cc = (double) seg.length() * (seg.length() - 1);
        //double cc = 1;
        return Math.abs(norm.inverseCumulativeProbability(pValue / cc));
    }

    private List<Segment> mergeScrapSegments(Set<Segment> segments, Set<Segment> readySegments) {

        List<Segment> tempList = new LinkedList<Segment>();
        List<Segment> result = new LinkedList<Segment>();

        tempList.addAll(segments);
        tempList.addAll(readySegments);
        Collections.sort(tempList);

        Segment mergeSeg, nextSeg;
        Iterator<Segment> its = tempList.listIterator();

        boolean isMerging;
        while (its.hasNext()) {
            mergeSeg = its.next();
            isMerging = false;

            if (!mergeSeg.isReady) {
                while (its.hasNext()) {
                    nextSeg = its.next();
                    if (!nextSeg.isReady) {
                        //mergeSeg.Loci.putAll(nextSeg.Loci);
                        mergeSeg.setRange(mergeSeg.Start(), nextSeg.End());
                        mergeSeg.isDirty = true;
                        isMerging = true;
                    } else {
                        if (isMerging) {
                            result.add(mergeSeg);
                        } else {
                            readySegments.add(mergeSeg);
                        }
                        isMerging = false;
                        mergeSeg = null;
                        break;
                    }
                }
                if (mergeSeg != null) {
                    if (isMerging) {
                        result.add(mergeSeg);
                    } else {
                        readySegments.add(mergeSeg);
                    }
                }
            }

        }

        return result;
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
}

