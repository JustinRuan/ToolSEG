package edu.whut.info.core;

import be.ac.ulg.montefiore.run.jahmm.Hmm;
import be.ac.ulg.montefiore.run.jahmm.ObservationReal;
import edu.whut.info.dataset.ObservationRealEx;
import edu.whut.info.dataset.Result;
import edu.whut.info.dataset.Segment;
import edu.whut.info.util.BioToolbox;
import edu.whut.info.util.ForwardBackwardCalculatorByBigDecimal;
import edu.whut.info.util.NewOpdfGaussian;
import net.sf.javaml.clustering.KMedoids;
import net.sf.javaml.core.Dataset;
import net.sf.javaml.core.DefaultDataset;
import net.sf.javaml.core.DenseInstance;
import net.sf.javaml.core.Instance;
import net.sf.javaml.distance.EuclideanDistance;
import org.apache.commons.math3.ml.clustering.CentroidCluster;

import java.math.BigDecimal;
import java.math.MathContext;
import java.util.*;
import java.util.logging.Logger;

/**
 * Created by Liu on 2016/3/20.
 */
public class HMMSegment implements SegmentCutter {
    private String methodName = "HMMSegment";
    private Logger m_log;
    private double centerProb;
    private double[] chromosome;

    public HMMSegment(double centerProb) {
        this.centerProb = centerProb;
        m_log = Logger.getLogger("segment");
    }

    private static double[] caculateMeanStd(Dataset data) {
        int len = data.size();
        double sumSquare = 0, sum = 0;
        for (Instance elem : data) {
            double value = elem.value(0);
            sum += value;
            sumSquare += value * value;
        }
        double mean = sum / len;
        return new double[]{mean, Math.sqrt(sumSquare / len - mean * mean)};
    }

    @Override
    public void splitChromosome(double[] data, Set<Segment> result, short chrId) {
        //chromosome=data;
        m_log.info(String.format("segment by %s", methodName));
        Segment input = new Segment();
        input.setChr_id(chrId);
        input.Seg_id = 0;
        double[] arr = BioToolbox.mean_std(data);
        input.CopyNumber = arr[0];
        input.stdCopyNumber = arr[1];
        input.setRange(0, data.length);
        chromosome = data;
        splitChromosome(input, result);
        m_log.info(String.format("#### #### #### chr %02d First Step: Loci Count = %05d; \t Segments Count = %d",
                chrId, data.length, result.size()));
    }

    @Override
    public void prepareCopyNumberSegment(double[] data) {

    }

    public ArrayList<Double> getsubdata(double[] data, int start, int stop) {
        ArrayList<Double> result = new ArrayList<>();
        for (int i = start; i < stop; i++)
            result.add(data[i]);
        return result;
    }

    public void splitChromosome1(Segment input, Set<Segment> output) {
        m_log.info(String.format(String.valueOf(centerProb)));
        ArrayList<Double> values = new ArrayList<Double>(getsubdata(chromosome, input.Start(), input.End()));
        SegmentProcessor.smoothing(values, 16);
        List<ObservationReal> Obser = new ArrayList<ObservationReal>();
        for (double d : values) {
            Obser.add(new ObservationReal(d));
        }
        int L = Obser.size();
        Hmm<ObservationReal> bestHmm;
        ArrayList<Double> arr = new ArrayList<Double>();
        List<Map.Entry<Hmm<ObservationReal>, Double>> Hmms = new ArrayList<>();
        int bestIndex = -1;
        double minPsi = Double.MAX_VALUE;
        for (int i = 1; i < 6; i++) {
            Hmm<ObservationReal> hmm = getObservationRealHmm(i, Obser);
            ForwardBackwardCalculatorByBigDecimal fbc = new ForwardBackwardCalculatorByBigDecimal(Obser, hmm);
            BigDecimal result2 = fbc.probabilityByBigDecimal();
            //����log
            int powers10 = result2.round(new MathContext(1)).scale() * -1;
            double result = result2.movePointLeft(powers10).doubleValue();
            double psi = 0;
//            if (result==0)
//                 psi=-1.0*powers10+2.0*i*i/Obser.size();
//            else
            psi = -1.0 * (Math.log10(result) + powers10) + 2.0 * i * i / Obser.size();
            Hmms.add(new HashMap.SimpleEntry<Hmm<ObservationReal>, Double>(hmm, psi));
            m_log.info(String.format("When K = %d, probability = %.3e, psi = %f", i, result2, psi));
            //get the min value
            if (psi < minPsi) {
                bestIndex = i;
                minPsi = psi;
            }
        }
        m_log.info("the state number of the bestHmm is " + bestIndex);
        bestHmm = Hmms.get(bestIndex - 1).getKey();
        int[] states = bestHmm.mostLikelyStateSequence(Obser);

        int start = 0;
        int startstate = states[0];
        for (int stop = 1; stop < states.length; stop++) {
            if (states[stop] != startstate) {
                output.add(input.getSubSegment(start, stop));
                start = stop;
                startstate = states[stop];
            }
        }
        output.add(input.getSubSegment(start, states.length));
        mergeSegment(output, 0.15);
        //merge the segment
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
            if (Math.abs(preSegment.CopyNumber - nextSegment.CopyNumber) <= threshold) {
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

    public void refreshSegment(Segment seg) {
        if (seg.isDirty) {
            double[] ms;
            ms = BioToolbox.robustMean(chromosome, seg.Start(), seg.End(), 16);
            seg.CopyNumber = ms[0];
            seg.stdCopyNumber = ms[1];

            seg.isDirty = false;
        }
    }

    public void splitChromosome(Segment input, Set<Segment> output) {
        m_log.info(String.format(String.valueOf(centerProb)));
        double[] temp = new double[chromosome.length];
        for (int i = 0; i < chromosome.length; i++)
            temp[i] = BioToolbox.log2(chromosome[i]);
        double[] data = BioToolbox.GaussianBlur(temp, 128, 1);
        List<Integer> candidates = highPassFiltering(data, 200, 0.8);
        List<Double> values = new ArrayList<Double>();
        for (int a : candidates)
            values.add(data[a]);

        List<ObservationReal> Obser = new ArrayList<ObservationReal>();
        for (double d : values) {
            Obser.add(new ObservationReal(d));
        }
        int L = Obser.size();
        Hmm<ObservationReal> bestHmm;
        ArrayList<Double> arr = new ArrayList<Double>();
        List<Map.Entry<Hmm<ObservationReal>, Double>> Hmms = new ArrayList<>();
        int bestIndex = -1;
        double minPsi = Double.MAX_VALUE;
        for (int i = 1; i < 6; i++) {
            Hmm<ObservationReal> hmm = getObservationRealHmm(i, Obser);
            ForwardBackwardCalculatorByBigDecimal fbc = new ForwardBackwardCalculatorByBigDecimal(Obser, hmm);
            BigDecimal result2 = fbc.probabilityByBigDecimal();
            //����log
            int powers10 = result2.round(new MathContext(1)).scale() * -1;
            double result = result2.movePointLeft(powers10).doubleValue();
            double psi = 0;
//            if (result==0)
//                 psi=-1.0*powers10+2.0*i*i/Obser.size();
//            else
            // psi = -1.0 *(Math.log10(result) + powers10)  + 2.0 * i *i/ Obser.size();
            psi = -1.0 * (Math.log10(result) + powers10);
            Hmms.add(new HashMap.SimpleEntry<Hmm<ObservationReal>, Double>(hmm, psi));
            m_log.info(String.format("When K = %d, probability = %.3e, psi = %f", i, result2, psi));
            //get the min value
            if (psi < minPsi) {
                bestIndex = i;
                minPsi = psi;
            }
        }
        //bestIndex=3;
        m_log.info("the state number of the bestHmm is " + bestIndex);
        bestHmm = Hmms.get(bestIndex - 1).getKey();
        int[] states = bestHmm.mostLikelyStateSequence(Obser);

        int start = 0;
        int startstate = states[0];
        for (int stop = 1; stop < states.length; stop++) {
            if (states[stop] != startstate) {
                output.add(input.getSubSegment(start, candidates.get(stop)));
                start = candidates.get(stop);
                startstate = states[stop];
            }
        }
        output.add(input.getSubSegment(start, chromosome.length));
        //   mergeSegment(output,0.3);
        //merge the segment
        int i = 1;
        for (Segment seg : output) {
            refreshSegment(seg);
            m_log.info(seg.getCharacterString());
            //  m_log.info(String.format("the %2d segment:\t start=%6d\t start=%6d\t mean=%.4f\t std=%.4f", i, seg.range.Start,seg.range.End,seg.HalfCopyNumber,seg.stdHalfCopyNumber));
            i++;
        }
    }

    public ArrayList<Integer> highPassFiltering(double[] data, int width, double percent) {
        int length = data.length;
        double[] values = new double[length];
        for (int i = 0; i < length; i++)
            values[i] = data[i];
        double[] diff = BioToolbox.difference(values, width);
        double[] sortDiff = new double[length];
        for (int i = 0; i < length; i++) {
            values[i] = Math.abs(diff[i]);
            sortDiff[i] = values[i];
        }
        Arrays.parallelSort(sortDiff);

        double threshold = sortDiff[(int) Math.round(sortDiff.length * percent + 0.5)];

        ArrayList<Integer> index = new ArrayList<>();//potential breakpoints
        index.add(0);
        for (int i = width >> 1; i < values.length - (width >> 1); i++) {
            if (values[i] > threshold)
                index.add(i);
        }
        // index.add(length);

        return index;
    }

    private Hmm<ObservationReal> getObservationRealHmm(List<ObservationReal> oseq) {
        double min = Double.MAX_VALUE;
        int k = 0;
        for (int K = 1; K < 6; K++) {
            KMedoids kms = new KMedoids(K, 1000, new EuclideanDistance());
            //KMeans kms = new KMeans(K,1000,new EuclideanDistance());
            Dataset data = new DefaultDataset();
            for (ObservationReal elem : oseq) {
                data.add(new DenseInstance(new double[]{elem.value}));
            }
            Dataset[] clustingResult = kms.cluster(data);
            double sum = 0;
            for (int i = 0; i < K; i++) {
                double mean = caculateMeanStd(clustingResult[i])[0];
                double std = caculateMeanStd(clustingResult[i])[1];
                m_log.info(String.format("center = %.4f; std = %.4f", mean, std));
                sum += std * std * clustingResult[i].size();
            }
            if (sum < min) {
                min = sum;
                k = K;
            }
        }
        m_log.info(String.format("the best cluster number is  %2d: ", k));
        return getObservationRealHmm(k, oseq);
    }

    private Hmm<ObservationReal> getObservationRealHmm(int K, List<ObservationReal> oseq) {
        KMedoids kms = new KMedoids(K, 1000, new EuclideanDistance());
        //KMeans kms = new KMeans(K,1000,new EuclideanDistance());
        Dataset data = new DefaultDataset();
        for (ObservationReal elem : oseq) {
            data.add(new DenseInstance(new double[]{elem.value}));
        }

        Dataset[] clustingResult = kms.cluster(data);
        double[][] mean_std = new double[K][2];
        for (int i = 0; i < K; i++) {
            mean_std[i] = caculateMeanStd(clustingResult[i]);
        }

        double[] pi = new double[K];
        double[][] A = new double[K][K];
        double m_pi = 1.0 / K;
        double diagRatio = centerProb; //�Խ����ϵ�Ԫ����ǶԽ�����Ԫ�ص�Ȩ�ر���8:1
        double mb = 1.0 / (diagRatio + K - 1);
        double ma = diagRatio * mb;
        for (int i = 0; i < K; i++) {
            pi[i] = m_pi;
            for (int j = 0; j < K; j++) {
                if (i == j) {
                    A[i][j] = ma;
                } else {
                    A[i][j] = mb;
                }
            }
        }

        List<NewOpdfGaussian> opdf = new ArrayList<>();
        for (int i = 0; i < K; i++) {
            double mean = mean_std[i][0];
            double std = mean_std[i][1] + Math.sqrt(0.4);
            opdf.add(new NewOpdfGaussian(mean, std));
            m_log.info(String.format("center = %.4f; std = %.4f", mean, std));
        }

        Hmm<ObservationReal> hmm = new Hmm<ObservationReal>(pi, A, opdf);
        return hmm;
    }

    private double getSquareLost(List<CentroidCluster<ObservationRealEx>> list) {
        double lost = 0;
        for (CentroidCluster<ObservationRealEx> li : list) {
            double mean = li.getCenter().getPoint()[0];
            List<ObservationRealEx> ll = li.getPoints();
            double sum = 0;
            for (ObservationRealEx obe : ll) {
                sum += (obe.value - mean) * (obe.value - mean);
            }
            lost += sum;
        }
        return lost;
    }

    private double std(CentroidCluster<ObservationRealEx> li) {
        double mean = li.getCenter().getPoint()[0];
        List<ObservationRealEx> ll = li.getPoints();
        int len = ll.size();
        double sum = 0;
        for (ObservationRealEx obe : ll) {
            sum += obe.value * obe.value;
        }
        return Math.sqrt(sum / len - mean * mean);
    }

    @Override
    public String getMethodName() {
        return methodName;
    }

    @Override
    public List<Result> getResult() {
        return null;
    }
}
