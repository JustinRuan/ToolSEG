package edu.whut.info.core;

import edu.whut.info.dataset.BinaryTreeNode;
import edu.whut.info.dataset.LinkedBinaryTreeNode;
import edu.whut.info.dataset.Segment;
import edu.whut.info.util.BioToolbox;
import edu.whut.info.util.Histogram;
import org.apache.commons.math3.distribution.NormalDistribution;

import java.io.*;
import java.util.*;
import java.util.function.BooleanSupplier;
import java.util.logging.Logger;

/**
 * Created by Justin on 2016/3/12.
 * Version : 23
 */
//Central limit theorem
public class BCLT implements SegmentCutter {
    private double PVALUE_THRESH = 0.05;
    private int MIN_SEG_LENGTH = 200;
    private String methodname = "BCLT";
    private double[] chromosome;
    private Logger m_log;
    private double[] IntegralCN;
    private double robustSTD;

    private NormalDistribution norm;
    private SegTree zTree;

    private double[] diff;

    public BCLT(double pvalueThre, int minSegLen) {
        m_log = Logger.getLogger("segment");

        this.PVALUE_THRESH = pvalueThre;
        this.MIN_SEG_LENGTH = minSegLen;
        norm = new NormalDistribution(0, 1);
        zTree = new SegTree();
    }

    @Override
    public void splitChromosome(double[] data, Set<Segment> result, short Chr_id) {
        m_log.info(String.format("segment by %s", methodname));
        chromosome = data;
        prepareCopyNumberSegment(data);

        Segment.resetID();
        Segment input = new Segment();
        input.setChr_id(Chr_id);
        input.Seg_id = 1;//根节点

        double[] arr = BioToolbox.mean_std(data);
        input.HalfCopyNumber = arr[0];
        input.stdHalfCopyNumber = arr[1];
        input.setRange(0, data.length);

        LinkedBinaryTreeNode<Integer> root = zTree.setRoot(input);
        splitChromosome(root, input, result, -1);

        printZHistogram();
        m_log.info(zTree.toString());

        m_log.info(String.format("#### #### #### chr %02d: Loci Count = %05d; \t Segments Count = %d",
                Chr_id, data.length, result.size()));

        zTree.clear();

        int i = 1;
        for (Segment seg : result) {
            refreshSegment(seg);
            m_log.info(seg.getCharacterString());
            // m_log.info(String.format("the %2d segment:\t start=%6d\t end=%6d\t mean=%.4f\t std=%.4f", i, seg.range.Start,seg.range.End,seg.HalfCopyNumber,seg.stdHalfCopyNumber));
            i++;
        }

    }


    @Override
    public String getMethodName() {
        return methodname;
    }

    @Override
    public void prepareCopyNumberSegment(double[] data) {
        int length = data.length;
        IntegralCN = new double[length + 1];
        IntegralCN[0] = 0.0;
        double sum = 0;
        for (int j = 0; j < length; j++) {
            sum += data[j];
            IntegralCN[j + 1] = sum;
        }
        robustSTD = getRobustStd();

    }

    private double getRobustStd() {
        int count = chromosome.length;
        if (count > 0) {

            diff = new double[count - 1];
            for (int k = 1; k < count; k++) {
                diff[k - 1] = chromosome[k] - chromosome[k - 1];
            }
            double[] df = diff.clone();

            Arrays.parallelSort(df);
            final double limit = 0.01;
            double[] ap = Arrays.copyOfRange(df, (int) (limit * count), (int) ((1 - limit) * count));
            double[] temp = BioToolbox.mean_std(ap);
            double std = temp[1] / Math.sqrt(2.0);
            m_log.info(String.format("Robust Std = %f", std));
            return std;

        }
        return -1;
    }

    private void splitChromosome(LinkedBinaryTreeNode<Integer> parent, Segment input, Set<Segment> output, double parentZ) {

        if (chromosome.length == 0) return;
        refreshSegment(input);

        double lamada = 0.2;
        double maxZ = 0;
        double maxZsearch = 0;
        int maxPos = 0;
        double z1, z2, z3;

        if (input.length() <= MIN_SEG_LENGTH) {
            //分段完成
            input.isReady = true;
            output.add(input);
            zTree.setZValue(parent, 0);
            zTree.setBreakPosition(parent, -1, 'C');
            return;
        }

        for (int j = input.Start() + MIN_SEG_LENGTH; j < input.End() - MIN_SEG_LENGTH; j++) {
            //如果过于接近整个段,就忽略

            z1 = calculateZ(input, input.Start(), j);
            z2 = calculateZ(input, j, input.End());

            z3 = Math.pow(Math.sqrt(Math.abs(z1)) + Math.sqrt(Math.abs(z2)), 2.0);

            if (z3 > maxZsearch) {
                maxZ = Math.max(Math.abs(z1), Math.abs(z2));
                maxZsearch = z3;
                maxPos = j;
            }
        }

        char windowModel = ' ';

        if (maxZ < robustSTD) {
            windowModel = 'W';

            //可以开启，辅助加窗二分
            List<Integer> winSet = new LinkedList<>();
            int w = input.length() >> 1;
            while (w > 10 * MIN_SEG_LENGTH) {
                w = w >> 1;
                winSet.add(w);
            }

            maxZsearch = 0;
            for (int width : winSet) {
                for (int j = input.Start() + MIN_SEG_LENGTH; j < input.End() - MIN_SEG_LENGTH; j++) {
                    //如果过于接近整个段,就忽略
                    int left = Math.max(input.Start(), j - width);
                    int right = Math.min(input.End(), j + width);
                    z1 = calculateZ(input, left, j);
                    z2 = calculateZ(input, j, right);

                    //z3 = Math.pow(Math.sqrt(Math.abs(z1)) + Math.sqrt(Math.abs(z2)), 2.0);
                    z3 = Math.max(Math.abs(z1), Math.abs(z2));
                    if (z3 > maxZsearch) {
                        maxZ = z3;
                        maxZsearch = z3;
                        maxPos = j;
                    }
                }
            }
        }
        zTree.setZValue(parent, maxZ);
        zTree.setBreakPosition(parent, maxPos, windowModel);

        if (maxZ >= robustSTD && maxPos > 0) {
            int newBreak = maxPos;
            //二分
            Segment front = input.getSubSegment(input.Start(), newBreak);
            LinkedBinaryTreeNode<Integer> left = zTree.setLeft(parent, front);
            splitChromosome(left, front, output, maxZ);

            Segment behind = input.getSubSegment(newBreak, input.End());
            LinkedBinaryTreeNode<Integer> right = zTree.setRight(parent, behind);
            splitChromosome(right, behind, output, maxZ);

        } else {
            //分段完成
            input.isReady = true;
            output.add(input);
        }
    }

    private void printZHistogram() {

        Histogram h = new Histogram(60, 0, 3);
        Map<Integer, Double> zPos = new TreeMap<>();
        for (Map.Entry<Integer, Double> kv : zTree.zMap.entrySet()) {
            int id = kv.getKey();
            double z = kv.getValue();
            if (z > 0) {
                h.addDataPoint(z);
                zPos.put(zTree.breakPositions.get(id), z/robustSTD);
            }
        }

        m_log.info(h.toString("z Value : robustSTD"));
    }

    private void save2File(double[][] data, int length) {
        BufferedWriter fw = null;

        try {
            File file1 = new File(".\\result\\zValueArray.txt");
            fw = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(file1, true), "UTF-8"));

            for (int i = 0; i < length; i++) {
                fw.append(String.format("%f\t%f", data[i][0], data[i][1]));
                fw.newLine();
            }

            fw.flush(); // 全部写入缓存中的内容
        } catch (Exception e) {
            e.printStackTrace();
        } finally {
            if (fw != null) {
                try {
                    fw.close();
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }
        }

    }

    public double calculateZ(Segment seg, int start, int stop) {

        int len = stop - start;

        double sum = 0;

        //求和范围是Start，Stop - 1；不包括Stop
        sum = IntegralCN[stop] - IntegralCN[start];

        ///return Math.abs((sum - len*avgHalfCopyNumber)/(stdHalfCopyNumber*Math.sqrt(len)))
        double cn = seg.HalfCopyNumber;

        //double cc = (double) seg.length() * (seg.length() - 1);
//        double cc = (double) len * (len - 1);
        double cc = (double) len;
//        double cc = 1;
        double thresh = Math.abs(norm.inverseCumulativeProbability(PVALUE_THRESH / cc));

        return ((sum - len * cn) / (thresh * Math.sqrt(len)));
    }

    public void refreshSegment(Segment seg) {
        if (seg.isDirty) {
            double[] ms;
            ms = BioToolbox.robustMean(chromosome, seg.Start(), seg.End(), 0);
//            if (ms == null){
//                System.out.println();
//            }
            seg.HalfCopyNumber = ms[0];
            seg.stdHalfCopyNumber = ms[1];

            seg.isDirty = false;
        }
    }


    public class SegTree {
        public Map<Integer, Double> zMap; // <id, zvalue of parent>
        public Map<Integer, Segment> segments; //
        public Map<Integer, Integer> breakPositions; //
        public Map<Integer, Integer> breakType;
        private LinkedBinaryTreeNode<Integer> mTreeRoot;

        public SegTree() {
            zMap = new TreeMap<>();
            segments = new HashMap<>();
            breakPositions = new TreeMap<>();
            breakType = new TreeMap<>();
            mTreeRoot = null;
        }

        public LinkedBinaryTreeNode<Integer> setRoot(Segment seg) {
            clear();

            mTreeRoot = new LinkedBinaryTreeNode<Integer>(1);
            segments.put(1, seg);
            return mTreeRoot;
        }

        public LinkedBinaryTreeNode<Integer> getRoot() {
            return mTreeRoot;
        }

        public void setZValue(LinkedBinaryTreeNode<Integer> id, double z) {
            zMap.put(id.getData(), z);
        }

        public void setBreakPosition(LinkedBinaryTreeNode<Integer> id, int pos, int type) {
            breakPositions.put(id.getData(), pos);
            breakType.put(id.getData(),type);
        }

        public void clear() {
            zMap.clear();
            segments.clear();
            breakPositions.clear();
            breakType.clear();
            mTreeRoot = null;
        }

        public LinkedBinaryTreeNode<Integer> setRight(LinkedBinaryTreeNode<Integer> parent, Segment seg) {

            LinkedBinaryTreeNode<Integer> child = new LinkedBinaryTreeNode<>(seg.Seg_id);
            parent.setRight(child);
            segments.put(seg.Seg_id, seg);
            return child;
        }

        public LinkedBinaryTreeNode<Integer> setLeft(LinkedBinaryTreeNode<Integer> parent, Segment seg) {

            LinkedBinaryTreeNode<Integer> child = new LinkedBinaryTreeNode<>(seg.Seg_id);
            parent.setLeft(child);
            segments.put(seg.Seg_id, seg);

            return child;
        }


        public String toString() {
            StringBuilder result = new StringBuilder();
            result.append('\n');

            LinkedList<Double> zArray = new LinkedList<Double>(zMap.values());
            Collections.sort(zArray);

            Map<Integer, String> zString = new LinkedHashMap<>();
            for (Map.Entry<Integer, Double> kv : zMap.entrySet()) {
                double z = kv.getValue();
                int id = kv.getKey();

                int len = zArray.indexOf(z);
                int tm = breakType.get(id);

                StringBuilder temp = new StringBuilder();
                for (int j = 0; j < len; j++) {
                    temp.append("-");
                }

                zString.put(id, String.format("\t z = %.4f [%c] <-%s [%02d]",
                        z, tm, temp.toString(), zArray.size() - len));
            }

            mTreeRoot.traverseInorder(new BinaryTreeNode.Visitor() {
                @Override
                public void visit(BinaryTreeNode node) {
                    int leftid, rightid, breakpos, startpos, endpos;

                    if (node.getRight() == null || node.getLeft() == null) {
                        startpos = segments.get(node.getData()).Start();
                        endpos = segments.get(node.getData()).End();
                        result.append(String.format("Seg: %4d = (XXXX:XXXX) <% 8d><--------><% 8d> %s",
                                node.getData(), startpos, endpos, zString.get(node.getData())));
                    } else {
                        leftid = (Integer) node.getLeft().getData();
                        rightid = (Integer) node.getRight().getData();
                        breakpos = segments.get(node.getRight().getData()).Start();
                        startpos = segments.get(node.getData()).Start();
                        endpos = segments.get(node.getData()).End();
                        result.append(String.format("Seg: %4d = (% 4d:% 4d) <% 8d><% 8d><% 8d> %s",
                                node.getData(), leftid, rightid, startpos, breakpos, endpos, zString.get(node.getData())));
                    }

                    result.append('\n');
                }
            });
            return result.toString();
        }
    }


}
