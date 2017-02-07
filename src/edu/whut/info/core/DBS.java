package edu.whut.info.core;

import edu.whut.info.dataset.BinaryTreeNode;
import edu.whut.info.dataset.LinkedBinaryTreeNode;
import edu.whut.info.dataset.Result;
import edu.whut.info.dataset.Segment;
import edu.whut.info.util.BioToolbox;
import edu.whut.info.util.Histogram;
import org.apache.commons.math3.distribution.NormalDistribution;

import java.io.*;
import java.util.*;
import java.util.logging.Logger;

/**
 * Created by Justin on 2016/3/12.
 * Version : 23
 */
//Central limit theorem
public class DBS implements SegmentCutter {
    private double PVALUE_THRESH = 0.05;
    private int MIN_SEG_LENGTH = 200;
    private double LAMBDA = 0.1;
    private String methodname = "DBS";
    private double[] chromosome;
    private Logger m_log;
    private double[] IntegralCN;
    private double robustSTD;

    private NormalDistribution norm;
    private SegTree zTree;

    private double[] diff;
    private final boolean Show_Debug = false;

    public DBS(double pvalueThre, int minSegLen, double lambda) {
        m_log = Logger.getLogger("segment");

        this.PVALUE_THRESH = pvalueThre;
        this.MIN_SEG_LENGTH = minSegLen;
        this.LAMBDA = lambda;
        norm = new NormalDistribution(0, 1);
        zTree = new SegTree();
    }

    @Override
    public void splitChromosome(double[] data, Set<Segment> result, short Chr_id) {
        zTree.clear();
        if (Show_Debug) m_log.info(String.format("segment by %s", methodname));
        chromosome = data;
        prepareCopyNumberSegment(data);

        Segment.resetID();
        Segment input = new Segment();
        input.setChr_id(Chr_id);
        input.Seg_id = 1;//根节点

        double[] arr = BioToolbox.mean_std(data);
        input.CopyNumber = arr[0];
        input.stdCopyNumber = arr[1];
        input.setRange(0, data.length);

        LinkedBinaryTreeNode<Integer> root = zTree.setRoot(input);
        splitChromosome(root, input, result, -1);

        double[] disLeaf = zTree.calculatDistance();
        double quality = disLeaf[0] - disLeaf[1];
        if (Show_Debug) {
            printZHistogram();
            m_log.info(zTree.toString());
            m_log.info(String.format("Quality = %.4f", quality));

            m_log.info(String.format("#### #### #### chr %02d: Loci Count = %05d; \t Segments Count = %d",
                    Chr_id, data.length, result.size()));
        }

        if (quality < 0){
            decideBreakPoints(result,Chr_id,disLeaf[1]);
        }


        if (Show_Debug) {
            int i = 1;
            for (Segment seg : result) {
                refreshSegment(seg);
                m_log.info(seg.getCharacterString());
                i++;
            }
        }
    }

    private void decideBreakPoints(Set<Segment> output, Short Chr_id, double thresh) {

        Set<Integer> breakPositions = new TreeSet<>();
        if (Show_Debug) {
            m_log.info(String.format("final thresh of Z = %.4f", thresh));
        }
        for (Map.Entry<Integer, Double> kv : zTree.zMap.entrySet()) {
            if (kv.getValue() > thresh + LAMBDA) {
                breakPositions.add(zTree.breakPositions.get(kv.getKey()));
            }
        }
        breakPositions.add(0);
        breakPositions.add(chromosome.length);

        Set<Segment> temp = new TreeSet<>();

        int start = 0;
        for (int end : breakPositions) {
            if (end > 0) {
                Segment newSeg = new Segment();
                newSeg.setChr_id(Chr_id);
                newSeg.setRange(start, end);

                refreshSegment(newSeg);
                temp.add(newSeg);
            }
            start = end;
        }
        output.clear();;
        output.addAll(temp);
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
        //robustSTD = Math.max(0.01,getRobustStd() + LAMBDA) ;
        double temp = getRobustStd();
        robustSTD = temp + 0.02;
        if (Show_Debug) m_log.info(String.format("Revised Robust Std = %f, robustSTD = %.3f", robustSTD, temp));
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
            if (Show_Debug) m_log.info(String.format("Robust Std = %f", std));
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

        if (input.length() < MIN_SEG_LENGTH) {
            //分段完成
            input.isReady = true;
            output.add(input);

            //并没有实际意义的计算，只是为了输出一个Z值
            for (int j = input.Start() + 1; j < input.End() - 1; j++) {
                //如果过于接近整个段,就忽略

                z1 = calculateZ(input, input.Start(), j);
                z2 = calculateZ(input, j, input.End());

                z3 = Math.pow(Math.sqrt(Math.abs(z1)) + Math.sqrt(Math.abs(z2)), 2.0);

                if (z3 > maxZsearch) {
                    maxZ = Math.max(Math.abs(z1), Math.abs(z2));
                    maxZsearch = z3;
                }
            }
            //没有意义的代码结束了！！！

            zTree.setZValue(parent, maxZ);
            zTree.setBreakPosition(parent, -1, 'C');
            return;
        }

        for (int j = input.Start() + 1; j < input.End() - 1; j++) {
        //for (int j = input.Start() + MIN_SEG_LENGTH; j < input.End() - MIN_SEG_LENGTH; j++) {
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
                //for (int j = input.Start() + 1; j < input.End() - 1; j++) {
                for (int j = input.Start() + MIN_SEG_LENGTH; j < input.End() - MIN_SEG_LENGTH; j++) {
                    //如果过于接近整个段,就忽略
                    int left = Math.max(input.Start(), j - width);
                    int right = Math.min(input.End(), j + width);
                    z1 = calculateZ(input, left, j);
                    z2 = calculateZ(input, j, right);

                    z3 = Math.pow(Math.sqrt(Math.abs(z1)) + Math.sqrt(Math.abs(z2)), 2.0);

                    if (z3 > maxZsearch) {
                        maxZ = Math.max(Math.abs(z1), Math.abs(z2));
                        maxZsearch = z3;
                        maxPos = j;
                    }
                }
            }
        }
//        if (maxZ > 0) {
//            zTree.setZValue(parent, maxZ);
//            zTree.setBreakPosition(parent, maxPos, windowModel);
//        } else {
//            zTree.setZValue(parent, input.stdCopyNumber);
//            zTree.setBreakPosition(parent, -1, 'C');
//        }
        if (maxZ <= 0) {
            zTree.setZValue(parent, 0);
            zTree.setBreakPosition(parent, -1, 'C');
        }

        if (maxZ > robustSTD && maxPos > 0 ) {
            zTree.setZValue(parent, maxZ);
            zTree.setBreakPosition(parent, maxPos, windowModel);


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

            zTree.setZValue(parent, maxZ);
            zTree.setBreakPosition(parent, maxPos, windowModel);
        }
    }

    private void printZHistogram() {
        Histogram h = new Histogram(60, 0, 3);
        for (Map.Entry<Integer, Double> kv : zTree.zMap.entrySet()) {
            int id = kv.getKey();
            double z = kv.getValue();
            if (z > 0) {
                h.addDataPoint(z);
            }
        }
        if (Show_Debug) m_log.info(h.toString("z Value"));
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
        double cn = seg.CopyNumber;

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
            seg.CopyNumber = ms[0];
            seg.stdCopyNumber = ms[1];

            seg.isDirty = false;
        }
    }

    @Override
    public List<Result> getResult() {

        return zTree.getResult();
    }


    public class SegTree {
        public Map<Integer, Double> zMap; // <id, zvalue of parent>
        public Map<Integer, Segment> segments; //
        public Map<Integer, Integer> breakPositions; //
        public Map<Integer, Integer> breakType;
        public LinkedBinaryTreeNode<Integer> mTreeRoot;

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
            breakType.put(id.getData(), type);
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

//        public double calculateRate(BinaryTreeNode node){
//            List<Double> zArray = new ArrayList<>();
//
//            List<Integer> path = getPath2Root(node);
//            for(int value : path){
//                double z = zMap.get(value);
//                zArray.add(z);
//            }
//            int count = zArray.size();
//
//            if (count == 1) {
//                return zArray.get(0);
//            }else if (count == 2){
//                return Math.abs(zArray.get(1) - zArray.get(0));
//            }else{
//                int step = (count >>1) + 1;
//                step = Math.min(step,3);
//
//                double sum1 = 0,sum2 = 0;
//                for (int i = 0; i < step; i++){
//                    sum1 += zArray.get(i);
//                    sum2 += zArray.get(i+1);
//                }
//                return Math.abs(sum2 - sum1)/step;
//            }
//        }

        public double calculateRate(BinaryTreeNode node) {
            BinaryTreeNode parent = node.getParent();

//            if (parent == null) {
//                return zMap.get(node.getData());
//            }else{
//                return Math.abs(zMap.get(parent.getData()) - zMap.get(node.getData()));
//            }
            if (parent == null) {
                return 1;
            } else {
                double r1 = segments.get(node.getData()).stdCopyNumber;
                double r2 = segments.get(parent.getData()).stdCopyNumber;
                return (1. - r1 / r2);
            }
        }

        public List<Integer> getPath2Root(BinaryTreeNode node) {
            List<Integer> temp = new ArrayList<>();

            BinaryTreeNode mySelf = node;
            while (mySelf != null) {
                temp.add((Integer) mySelf.getData());
                mySelf = mySelf.getParent();
            }
            return temp;
        }

        public double[] calculatDistance() {
            if (zTree.zMap.size() == 1) return new double[]{0,0};

            final double[] result = new double[2];
            result[0] = Double.MAX_VALUE;
            result[1] = Double.MIN_VALUE;

            mTreeRoot.traverseInorder(new BinaryTreeNode.Visitor() {
                @Override
                public void visit(BinaryTreeNode node) {
                    double z = zMap.get(node.getData());
                    int type = breakType.get(node.getData());
                    double std = segments.get(node.getData()).stdCopyNumber;

                    if (node.getLeft() != null && node.getRight() != null) {
                        result[0] = (z < result[0]) ? z : result[0];//非叶子节点
                    } else {//叶子节点
                        result[1] = (std > result[1]) ? std : result[1];
                    }
                }
            });
            return result;
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

                StringBuilder temp = new StringBuilder();
                for (int j = 0; j < len; j++) {
                    temp.append("-");
                }

//                zString.put(id, String.format("\t z = %.4f [%c] <-%s [%02d]",
//                        z, tm, temp.toString(), zArray.size() - len));
                zString.put(id, String.format("\t<-%s [%02d]",
                        temp.toString(), zArray.size() - len));
            }

            mTreeRoot.traverseInorder(new BinaryTreeNode.Visitor() {
                @Override
                public void visit(BinaryTreeNode node) {
                    int leftid, rightid, breakpos, startpos, endpos;

                    int id = (int) node.getData();
                    double z = zMap.get(id);
                    //double rate = calculateRate(node);
                    double std = segments.get(id).stdCopyNumber;
                    int tm = breakType.get(id);

                    if (node.getRight() == null || node.getLeft() == null) {
                        startpos = segments.get(node.getData()).Start();
                        endpos = segments.get(node.getData()).End();
                        result.append(String.format("Seg: %4d = (XXXX:XXXX) <% 8d><--------><% 8d> (id,z,r) = (% 4d,%.4f,%+1.4f)[%c]%s",
                                node.getData(), startpos, endpos, node.getData(), z, std, tm, zString.get(node.getData())));
                    } else {
                        leftid = (Integer) node.getLeft().getData();
                        rightid = (Integer) node.getRight().getData();
                        breakpos = segments.get(node.getRight().getData()).Start();
                        startpos = segments.get(node.getData()).Start();
                        endpos = segments.get(node.getData()).End();
                        result.append(String.format("Seg: %4d = (% 4d:% 4d) <% 8d><% 8d><% 8d> (id,z,r) = (% 4d,%.4f,%+1.4f)[%c]%s",
                                node.getData(), leftid, rightid, startpos, breakpos, endpos, node.getData(), z, std, tm, zString.get(node.getData())));
                    }

                    result.append('\n');
                }
            });
            return result.toString();
        }

        public List<Result> getResult() {

            List<Result> rList = new LinkedList<>();

            mTreeRoot.traverseInorder(new BinaryTreeNode.Visitor() {
                @Override
                public void visit(BinaryTreeNode node) {
                    int breakpos;

                    int id = (int) node.getData();
                    double z = zMap.get(id);
                    //int tm = breakType.get(id);

                    Result r = new Result();
                    r.id = (int) node.getData();
                    r.value1 = z / robustSTD;

                    if (node.getRight() == null || node.getLeft() == null) {
                        r.isBreakPoint = false;
                        r.pos = 0;
                    } else {
                        breakpos = segments.get(node.getRight().getData()).Start();
                        r.isBreakPoint = true;
                        r.pos = breakpos;
                    }

                    rList.add(r);
                }
            });

            int[] posArray = new int[rList.size() + 2];
            posArray[0] = 0;
            posArray[1] = 0;
            posArray[posArray.length - 1] =  chromosome.length-1;
            for (int i = 1; i < rList.size(); i++){
                posArray[i+1] = rList.get(i).pos;
            }

            for (int i = 0; i < rList.size(); i++){
                if (rList.get(i).pos == 0){
                    rList.get(i).range = new int[]{posArray[i], posArray[i+2]};
                }
            }

            return rList;
        }
    }


}
