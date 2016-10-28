package edu.whut.info.core;

import com.google.common.primitives.Doubles;
import edu.whut.info.dataset.Chromosome;
import edu.whut.info.dataset.Result;
import edu.whut.info.dataset.Segment;
import edu.whut.info.util.BioToolbox;
import org.apache.commons.math3.analysis.function.Max;

import javax.imageio.ImageIO;
import java.awt.*;
import java.awt.image.BufferedImage;
import java.io.File;
import java.text.SimpleDateFormat;
import java.util.*;
import java.util.List;
import java.util.logging.Logger;

/**
 * Created by Liu on 2016/3/19.
 */
public class CNSegment {
    public SegmentCutter cutter;
    private Logger m_log;
    private List<Set<Segment>> result;
    private List<double[]> cacheSample;

    public CNSegment() {
        result = new ArrayList<>();
        m_log = Logger.getLogger("segment");
    }

    public List<Set<Segment>> getResult() {
        return result;
    }

    public void setCutter(SegmentCutter cutter) {
        this.cutter = cutter;
    }

    public void preprocessing(ArrayList<Chromosome> chros, double ratio, int method) {
        cacheSample = new ArrayList<double[]>();

        for (int i = 0; i < chros.size(); i++) {
            if (chros.get(i).probes.size() == 0)
                continue;
            List<Double> values = chros.get(i).probes;

            double[] cnArray = Doubles.toArray(values);
            if (ratio > 0) {
                cnArray = BioToolbox.limitFiltering(cnArray, ratio);
            }

            //transform
            switch (method) {
                case 1:
                    for (int j = 0; j < cnArray.length; j++) {
                        cnArray[j] = BioToolbox.log2(cnArray[j]);
                    }
                    break;
                case 2:
                    for (int j = 0; j < cnArray.length; j++) {
                        cnArray[j] = Math.pow(2, cnArray[j]);
                    }
                    break;
                default:
                    //do nothing;
            }
            cacheSample.add(cnArray);
        }
    }

    public void analysisResult(List<Long> breakPoints, List<Result> results, int chrid) {
        List<Integer> bps = new LinkedList<>();
        for (long t : breakPoints) {
            bps.add((int) t - 1);
        }

        for (Result r : results) {
            if (r.pos > 0) {
                int Maxdist = Integer.MAX_VALUE;
                int index = -1;
                for (int i = 0; i < bps.size(); i++) {
                    int d = Math.abs(bps.get(i) - r.pos);
                    if (d < Maxdist){
                        Maxdist = d;
                        index = i;
                    }
                }
                if (Maxdist < 30 ){
                    r.nearestBreakPoint = bps.get(index);
                    bps.remove(index);
                }
            }
        }

        if (!bps.isEmpty()){
            for (int pos : bps){
                Result r = new Result();
                r.pos = 0;
                r.nearestBreakPoint = pos;

                Result temp = null;
                for (Result r2 : results) {
                    if (r2.pos > pos) {
                        break;
                    }else if (r2.pos == 0){
                        temp = r2;
                    }
                }
                if (temp != null){
                    r.value1 = temp.value1;
                }else{
                    r.value1 = -1;
                }
                results.add(r);
            }
        }

        for (Result r : results){
            m_log.info(String.format("\t% 4d\t% 4d\t%d\t%d\t%f",chrid,r.id,r.pos,r.nearestBreakPoint,r.value1));
        }
    }
    public void printOriginalSegment(Chromosome chro){
        Set<Segment> result=new TreeSet<>();
        List<Long> changepoints=chro.changepoints;
        Segment input = new Segment();
        input.setChr_id(chro.chrId);
        input.setRange(0, chro.getLength());
        long start=0,end=0;
        for(int i=0;i<changepoints.size();i++){
             end=changepoints.get(i);
             result.add(input.getSubSegment((int)start,(int)end));
            start=end;
        }
        result.add(input.getSubSegment((int)end,chro.getLength()));
        int i = 1;
        double[] data=new double[chro.getLength()];
        for (int j=0;j<chro.getLength();j++)
            data[j]=chro.probes.get(j);
        for (Segment seg : result) {
            BioToolbox.refreshSegment(seg,data);
            m_log.info(seg.getCharacterString());
            i++;
        }
    }
    public void splitChromosome(ArrayList<Chromosome> chros, double ratio, int method) {
        result.clear();
        //loading(chros);
        preprocessing(chros, ratio, method);

        int nums = cacheSample.size();
        for (int i = 0; i < nums; i++) {
            result.add(new TreeSet<Segment>());
        }

        ArrayList<Short> chrIds = new ArrayList<>();
        for (int i = 0; i < nums; i++)
            chrIds.add(chros.get(i).chrId);
        for (int i = 0; i < nums; i++) {
           // printOriginalSegment(chros.get(i));
            long t1 = System.currentTimeMillis();
            cutter.splitChromosome(cacheSample.get(i), result.get(i), chrIds.get(i));
            long t2 = System.currentTimeMillis();
            //m_log.info(String.format("Time = %d ms", t2 - t1));

            List<Long> breakPoints = chros.get(i).changepoints;
            List<Result> rList = cutter.getResult();

            analysisResult(breakPoints, rList, i);
        }

        //画出分段结果
        //drawProbeSets(chros, result, method);


        //只分一段
        //    cutter.splitChromosome(cacheSample.get(1),result.get(1),chrIds.get(1));
        // drawChromosome(chros.get(2),result.get(1),30);
    }

    public void splitChromosome1(Chromosome chro, Set<Segment> result) {
        List<Double> values = chro.probes;
        int length = values.size();
        double[] cnArray = new double[length];
        for (int j = 0; j < length; j++) {
            cnArray[j] = values.get(j);
        }
        cutter.splitChromosome(cnArray, result, chro.chrId);
    }

    public List<Set<Segment>> resultclone() {
        List<Set<Segment>> temp = new ArrayList<>();
        for (int i = 0; i < result.size(); i++) {
            Set<Segment> chros = new TreeSet<Segment>();
            for (Segment segm : result.get(i)) {
                Segment seg = new Segment();
                seg.range.Start = segm.Start();
                seg.range.End = segm.End();
                seg.stdCopyNumber = segm.stdCopyNumber;
                seg.CopyNumber = segm.CopyNumber;
                chros.add(segm);
            }
            temp.add(chros);
        }
        return temp;
    }

    public void loading(ArrayList<Chromosome> chros) {
        cacheSample = new ArrayList<double[]>();
        for (int i = 0; i < chros.size(); i++) {
            if (chros.get(i).probes.size() == 0)
                continue;
            List<Double> values = chros.get(i).probes;
            int length = values.size();
            double[] cnArray = new double[length];
            for (int j = 0; j < length; j++) {
                cnArray[j] = values.get(j);
            }
            cacheSample.add(cnArray);
        }
    }

    public Set<Segment> getSegment(int i) {
        return result.get(i);
    }

//    public String drawChromosome(Chromosome chro, Set<Segment> segments, int step) {
//        int maxLength = 0;
//        maxLength = chro.probes.size();
//        int estimatedLength = maxLength;
//        int width = (int) (1.15 * estimatedLength / step);
//        int height = 400;
//        BufferedImage image = new BufferedImage(width, height, BufferedImage.TYPE_INT_RGB);
//        Graphics2D g = (Graphics2D) image.getGraphics();
//        g.setColor(Color.WHITE);
//        g.fillRect(0, 0, width, height);
//        Stroke bs;
//        bs = new BasicStroke(2.0f);
//        g.setStroke(bs);
//        int x, sx;
//        Font f1 = new Font(null, Font.BOLD, 24);
//        Font f2 = new Font(null, Font.BOLD, 12);
//        x = 100;
//        sx = 100 * step;
//        g.setColor(Color.BLACK);
//        g.setFont(f1);
//        for (double cnValue : chro.probes) {
//            int y = calculateYPosition(0, cnValue);
//            g.setColor(Color.DARK_GRAY);
//            g.fillRect(x, y, 2, 2);
//            sx++;
//            x = sx / step;
//        }
//        x = 100;
//        sx = 100 * step;
//        if (segments.size() != 0)
//            for (Segment s : segments) {
//                double cnValue = s.CopyNumber * 2;
//                //int y = calculateYPosition(index,s.robustMeanHalfCopyNumber);
//                int y2 = calculateYPosition(0, BioToolbox.log2(cnValue));
//                g.setColor(Color.RED);
//                int w = (int) (s.length() / step + 0.5);
//                g.drawLine(x, y2, x + w, y2);
//                g.setFont(f2);
//                int textPos = (w > 30) ? x + w / 2 : x;
//                g.drawString(String.format("[%.2f]", cnValue), textPos, y2 - 6);
//                sx = sx + s.length();
//                x = (int) ((sx / step) + 0.5);
//            }
//        String filename = "";
//        try {
//            Date date = new Date();
//            SimpleDateFormat df = new SimpleDateFormat("yyyyMMdd_HHmmss");
//            filename = String.format(".%1$sResult%1$sChromosome_%2$s.png", File.separator, df.format(date));
//            ImageIO.write(image, "PNG", new File(filename));
//        } catch (Exception ex) {
//            ex.printStackTrace();
//        }
//        return filename;
//    }

    private int calculateYPosition(int row, double value) {
        // return (int) (value * 50) + 400*(row - 1);
        return (int) (340 - value * 50) + 400 * row;
    }

    public void drawProbeSets(ArrayList<Chromosome> chros, List<Set<Segment>> result, int method) {
        if (chros.size() == 0)
            return;

        int maxLength = 0;
        for (Chromosome chro : chros) {
            maxLength = Math.max(maxLength, chro.probes.size());
        }
        int estimatedLength = (int) (1.15 * maxLength);

        int step;
        if (estimatedLength < 2000) {
            step = 1;
        } else {
            step = estimatedLength / 2000;
        }

        int width = estimatedLength / step;
        int height = 400 * result.size();
        BufferedImage image = new BufferedImage(width, height, BufferedImage.TYPE_INT_RGB);
        Graphics2D g = (Graphics2D) image.getGraphics();
        g.setColor(Color.WHITE);
        g.fillRect(0, 0, width, height);

        Stroke bs;
        bs = new BasicStroke(2.0f);
//        bs2 = new BasicStroke(1, BasicStroke.CAP_BUTT,
//                BasicStroke.JOIN_BEVEL, 0,
//                new float[]{16, 4}, 0);
        g.setStroke(bs);

        int x, sx;

        Font f1 = new Font(null, Font.BOLD, 24);
        Font f2 = new Font(null, Font.BOLD, 12);
        for (int i = 0; i < chros.size(); i++) {
            Chromosome chro = chros.get(i);
            x = 100;
            sx = 100 * step;
            g.setColor(Color.BLACK);
            g.setFont(f1);
            for (double cnValue : chro.probes) {
                int y = calculateYPosition(i, cnValue);
                g.setColor(Color.DARK_GRAY);
                g.fillRect(x, y, 2, 2);
                sx++;
                x = sx / step;
            }

            x = 100;
            sx = 100 * step;
            for (Segment s : result.get(i)) {
                double cnValue;

                //transform
                switch (method) {
                    case 1:
                        cnValue = Math.pow(2, s.CopyNumber);
                        break;
                    case 2:
                        cnValue = BioToolbox.log2(s.CopyNumber);
                        break;
                    default:
                        cnValue = s.CopyNumber;
                }

                int y2 = calculateYPosition(i, cnValue);
                g.setColor(Color.RED);
                int w = (int) (s.length() / step + 0.5);
                g.drawLine(x, y2, x + w, y2);
                g.setFont(f2);
                int textPos = (w > 30) ? x + w / 2 : x;
                g.drawString(String.format("[%.2f]", cnValue), textPos, y2 - 6);
                sx = sx + s.length();
                x = (int) ((sx / step) + 0.5);
            }
        }
        try {
            Date date = new Date();
            SimpleDateFormat df = new SimpleDateFormat("yyyyMMdd_HHmmss");
            String filename = String.format(".%1$sPictures%1$s%2$s%3$s.png",
                    File.separator, cutter.getMethodName(), df.format(date));

            ImageIO.write(image, "PNG", new File(filename));
        } catch (Exception ex) {
            ex.printStackTrace();
        }

    }

}

