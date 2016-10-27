/**
 *
 */
package edu.whut.info.util;

import edu.whut.info.dataset.Segment;
import org.apache.commons.math3.analysis.function.Gaussian;
import org.apache.commons.math3.util.MathArrays;

import java.util.*;

/**
 * The type Bio toolbox.
 *
 * @author Justin
 */
public class BioToolbox {
    /**
     * The constant ChromosomeLength.:Approximate number of each chromosome SNPs
     */
    public final static int[] ChromosomeLength
            = new int[]{16,
            150000, 160000, 130000, 120000,
            120000, 120000, 100000, 100000,
            90000, 100000, 90000, 90000,
            90000, 70000, 60000, 60000,
            50000, 60000, 30000, 50000,
            30000, 30000,
            16,
            16, 16, 16};

    public final static double LOG_2 = Math.log(2);
    public final double RobustStdMultiplier = 2.5;

    public static double log2(double value) {
        return Math.log(value) / LOG_2;
    }

    /**
     * Two-dimensional array convert to a one-dimensional array
     *
     * @param Chr_id chromosome id [1~22]
     * @param Start  Position in the own chromosome
     * @return Global unique number of Loci
     */
    public static Long getGlobalLocation(short Chr_id, int Start) {
        return ((long) Chr_id << 480) | (long) Start;
    }

    /**
     * Get chr id and start position.
     *
     * @param pos Global unique number of Loci
     * @return [Chr_id, Start]
     */
    public static int[] getChrIdAndStartPosition(long pos) {
        int Chr_id = (int) (pos >> 480);
        //int Start = (int)(pos & 0xffffffff);
        int Start = (int) pos;

        return new int[]{Chr_id, Start};
    }

    /**
     * Mean double.
     *
     * @param <T>    the type parameter
     * @param values the values
     * @return the double
     */
    public static <T extends Number> double mean(List<T> values) {
        int len = values.size();
        if (len > 0) {
            double sum = 0;
            for (T d : values) {
                sum += d.doubleValue();
            }

            return sum / len;
        } else
            return 0;
    }

    public static double[] difference(double[] x, int width) {

        int len = width | 0x01; //make sure the width of window is odd.
        int halfLen = len >> 1;

        double[] h = new double[len];
        double f = 1.0 / len;

        Arrays.fill(h, 0, halfLen, f);
        h[halfLen] = 0.0;
        Arrays.fill(h, halfLen + 1, len, -f);

        double[] result = MathArrays.convolve(x, h);

        return Arrays.copyOfRange(result, halfLen, x.length + halfLen);
    }
//
//	public static double mean(List<Float> values) {
//
//		if (values.size() > 0) {
//			double sum = 0;
//			for (Float d : values) {
//				sum += d;
//			}
//
//			return sum/values.size();
//		}
//		else
//			return 0;
//	}

    /**
     * mean of Square values.
     *
     * @param <T>    the type parameter
     * @param values the values
     * @return the double
     */
    public static <T extends Number> double squareMean(List<T> values) {
        int len = values.size();
        if (len > 0) {
            double sum = 0;
            double v;
            for (T d : values) {
                v = d.doubleValue();
                sum += v * v;
            }

            return sum / len;
        } else
            return 0;
    }

    /**
     * Std double.
     *
     * @param <T>    the type parameter
     * @param values the values
     * @return the double
     */
    public static <T extends Number> double std(List<T> values) {
        //using this equation : Var(x) = E(x^2)-(E(x))^2
        int len = values.size();
        if (len > 0) {
            double m = 0;
            double s = 0;
            double v;
            for (T d : values) {
                v = d.doubleValue();
                m += v;
                s += v * v;
            }

            return Math.sqrt((s - m * m / len) / len);
        } else
            return 0;
    }

    /**
     * get Mean and std of the List values.
     *
     * @param <T>    the type parameter
     * @param values the values
     * @return the double [ ]
     */
    public static <T extends Number> double[] mean_std(List<T> values) {
        //using this equation : Var(x) = E(x^2)-(E(x))^2
        int len = values.size();
        if (len > 0) {
            double m = 0;
            double s = 0;
            double v;
            for (T d : values) {
                v = d.doubleValue();
                m += v;
                s += v * v;
            }
            double mean = m / len;

            return new double[]{mean, Math.sqrt(s / len - mean * mean)};
        } else
            return null;
    }

    public static double[] mean_std(double[] values) {
        //using this equation : Var(x) = E(x^2)-(E(x))^2
        int len = values.length;
        if (len > 0) {
            double m = 0;
            double s = 0;
            double v;
            for (double d : values) {
                m += d;
                s += d * d;
            }
            double mean = m / len;

            return new double[]{mean, Math.sqrt(s / len - mean * mean)};
        } else
            return null;
    }

    /**
     * Std double. Var(x) = E(x^2)-(E(x))^2
     *
     * @param smean the smean of Square values
     * @param mean  the mean
     * @return double
     */
    public static double std(double smean, double mean) {
        return Math.sqrt(smean - mean * mean);
    }

    /**
     * get minimum and maximum double [ ].
     *
     * @param values the values
     * @return the double [ ]
     */
    public static double[] minmax(List<Double> values) {

        if (values.size() > 0) {
            double[] result = new double[2];
            result[0] = values.get(0);
            result[1] = values.get(0);

            for (Double d : values) {
                if (d < result[0]) result[0] = d;
                if (d > result[1]) result[1] = d;
            }

            return result;
        } else
            return null;
    }

//	public static String toSQLString(String tablename, ProbeSet ps) {
//
//		StringBuilder str = new StringBuilder();
//
//		str.append(String.format("INSERT INTO %s_ProbeSets ", tablename));
//		str.append("(Loci, ProbeSet_ID, OtherProbeSetID, Chr_id, Start, Is_CNProbe, Is_SNPProbe) VALUES ");
//		str.append(String.format("(%d, '%s', '%s', %d, %d, %d, %d);\n",ps.Loci,ps.ProbeSetID,ps.OtherProbeSetID,
//															ps.Chr_id,ps.Start,ps.Is_CNProbe,ps.Is_SNPProbe));
//
//		str.append(String.format("NSERT INTO %s_Intensity", tablename));
//		str.append("(Loci, isNormalSample, mean_A, mean_B, var_A, var_B, mean_AA, mean_BB, mean_CN, var_CN) VALUES ");
//		str.append(String.format("(%d, %d, %f, %f, %f, %f, %f, %f, %f, %f);\n",ps.Loci,1,ps.data.Normal.mean_A,ps.data.Normal.mean_B,
//					ps.data.Normal.std_A,ps.data.Normal.std_B,ps.data.Normal.mean_AA,ps.data.Normal.mean_BB,ps.data.Normal.mean_CN,ps.data.Normal.std_CN));
//
//		str.append(String.format("NSERT INTO %s_Intensity", tablename));
//		str.append("(Loci, isNormalSample, mean_A, mean_B, var_A, var_B, mean_AA, mean_BB, mean_CN, var_CN) VALUES ");
//		str.append(String.format("(%d, %d, %f, %f, %f, %f, %f, %f, %f, %f);\n",ps.Loci,0,ps.data.Tumor.mean_A,ps.data.Tumor.mean_B,
//					ps.data.Tumor.std_A,ps.data.Tumor.std_B,ps.data.Tumor.mean_AA,ps.data.Tumor.mean_BB,ps.data.Tumor.mean_CN,ps.data.Tumor.std_CN));
//
//		return str.toString();
//	}

//	public static String createNewTableSQLString(String dbpath, String tablename) {
//
//		StringBuilder tableString = new StringBuilder();
//
//		tableString.append(String.format("DROP TABLE IF EXISTS %s_ProbeSets;\n", tablename));
//		tableString.append(String.format("CREATE TABLE %s_ProbeSets (", tablename));
//		tableString.append("Loci integer,");
//		tableString.append("ProbeSet_ID varchar(100),");
//		tableString.append("OtherProbeSetID varchar(100),");
//		tableString.append("Chr_id integer,");
//		tableString.append("Start integer,");
//		tableString.append("Is_CNProbe boolean,");
//		tableString.append("Is_SNPProbe boolean,");
//		tableString.append("PRIMARY KEY (Loci)); \n");
//
//		tableString.append(String.format("DROP TABLE IF EXISTS %s_Intensity;", tablename));
//		tableString.append(String.format("CREATE TABLE %s_Intensity (", tablename));
//		tableString.append("Loci integer,");
//		tableString.append("isNormalSample boolean,");
//		tableString.append("mean_A real,");
//		tableString.append("mean_B real,");
//		tableString.append("var_A real,");
//		tableString.append("var_B real,");
//		tableString.append("mean_AA real,");
//		tableString.append("mean_BB real,");
//		tableString.append("mean_CN real,");
//		tableString.append("var_CN real,");
//		tableString.append("PRIMARY KEY (Loci, isNormalSample));\n");
//
//		return tableString.toString();
//	}

    /**
     * Save Biocollection to a Bin file.
     *
     * @param bc       the bc
     * @param fileanme the fileanme and path
     */
//    public static void saveBioCollectionToFile(BioCollection bc, String fileanme) {
//        System.out.println("Saving ... ...");
//        try {
//            FileOutputStream fos = new FileOutputStream(fileanme);
//            ObjectOutputStream oos = new ObjectOutputStream(fos);
//
//            oos.writeObject(bc);
//
//            oos.close();
//
//        } catch (IOException ex) {
//            ex.printStackTrace();
//        }
//    }


    /**
     * Read Biocollection from a Bin file.
     *
     * @param filename the filename
     * @return the bio collection
     */
//    public static BioCollection readBioCollectionFromFile(String filename) {
//        System.out.println("Reading ... ...");
//        try {
//            FileInputStream fis = new FileInputStream(filename);
//
//            ObjectInputStream ois = new ObjectInputStream(fis);
//
//            return (BioCollection) ois.readObject();
//        } catch (Exception ex) {
//            ex.printStackTrace();
//
//        }
//        return null;
//    }


//
//	public static void genotypeCalling(ArrayList<Chromosome> data) {
//		HashMap<Pair<Short,Long>, Double[]> LogIntensity;
//		final double log2 = Math.log(2);
//
//		LogIntensity = new HashMap<Pair<Short,Long>, Double[]>();
//		double sumA = 0;
//		double sumB = 0;
//		double minA = 1000;
//		double minB = 1000;
//		int count = 0;
//
//		for (int i = 1; i< 27; i++) {
//			if (!data.get(i).probes.isEmpty()) {
//
//				for (ProbeSet ps:data.get(i).probes.values()) {
//					if (ps.data.Is_Valid) {
//						Double[] values = new Double[2];
//						values[0] = (ps.data.Normal.mean_A > 0)? Math.log(ps.data.Normal.mean_A)/log2:0;
//						values[1] = (ps.data.Normal.mean_B > 0)? Math.log(ps.data.Normal.mean_B)/log2:0;
//
//						LogIntensity.put(Pair.create(ps.Chr_id, ps.Loci), values);
//
//						sumA += values[0];
//						sumB += values[1];
//
//						if (values[0] > 0) minA = Math.min(minA, values[0]);
//						if (values[1] > 0) minB = Math.min(minB, values[1]);
//						count ++;
//					}
//				}
//			}
//		}
//
//		double meanA = sumA/count;
//		double meanB = sumB/count;
//		double k = meanB / meanA;
//
//		//y - y1 = k(x-x1)
//		//y = k(x - minA) + minB;
//		//kx - y + minB - k*minA = 0;
//		double C = minB - k*minA;
//
//		Map<Double,Pair<Short,Long>> distance = new TreeMap<Double,Pair<Short,Long>>();
//		for (Map.Entry<Pair<Short,Long>, Double[]> ap : LogIntensity.entrySet()) {
//			double d;
//			Double[] x = ap.getValue();
//			d = Math.abs(k*x[0] - x[1] + C);
//			distance.put(d, ap.getKey());
//		}
//
//		int ABcount = count / 5;
//		Iterator<Pair<Short,Long>> ABit = distance.values().iterator();
//		for (int i = 0; i < ABcount; i++) {
//			Pair<Short,Long> p = ABit.next();
//			data.get(p.getFirst()).probes.get(p.getSecond()).data.Normal.isGenotypeAB = true;
//		}
////		int A2B2count = distance.size() -  2*count / 3;
////		ListIterator<Pair<Short,Long>> A2B2it = distance.values().
////		for (int i = distance.size() - 1; i> A2B2count; i--){
////
////		}
//
//	}

    /**
     * Blur double [ ].
     *
     * @param x     the x
     * @param width the width of window
     * @return the double [ ]
     */
    public static double[] blur(double[] x, int width) {

        int len = width | 0x01; //make sure the width of window is odd.
        int halfLen = len >> 1;

        double[] h = new double[len];
        double f = 1.0 / len;
        for (int i = 0; i < len; i++) {
            h[i] = f;
        }

        double[] result = MathArrays.convolve(x, h);

        return Arrays.copyOfRange(result, halfLen, x.length + halfLen);

    }

    /**
     * Gaussian blur.
     *
     * @param <T>   the type parameter
     * @param x     the x
     * @param width the width of window
     * @param sigma the sigma
     * @return the double [ ]
     */
    public static <T extends Number> double[] GaussianBlur(T[] x, int width, double sigma) {
        int len = width | 0x01;
        int halfLen = len >> 1;

        double[] h = new double[len];

        Gaussian g = new Gaussian(0, sigma);
        double step = 3.0 / halfLen;
        double pos = step;

        h[halfLen] = g.value(0);
        for (int i = 1; i <= halfLen; i++, pos += step) {
            double v = g.value(pos);
            h[halfLen + i] = v;
            h[halfLen - i] = v;
        }

        h = MathArrays.normalizeArray(h, 1);
        //泛型的代价？
        double[] xx = new double[x.length];
        for (int i = 0; i < x.length; i++) {
            xx[i] = x[i].doubleValue();
        }
        double[] result = MathArrays.convolve(xx, h);

        return Arrays.copyOfRange(result, halfLen, x.length + halfLen);
    }

    /**
     * Gaussian blur.
     *
     * @param x     the x
     * @param width the width
     * @param sigma the sigma
     * @return the double [ ]
     */
    public static double[] GaussianBlur(double[] x, int width, double sigma) {

        int len = width | 0x01;
        int halflen = len >> 1;

        double[] h = new double[len];

        Gaussian g = new Gaussian(0, sigma);
        double step = 3.0 / halflen;
        double pos = step;

        h[halflen] = g.value(0);
        for (int i = 1; i <= halflen; i++, pos += step) {
            double v = g.value(pos);
            h[halflen + i] = v;
            h[halflen - i] = v;
        }

        h = MathArrays.normalizeArray(h, 1);

        double[] result = MathArrays.convolve(x, h);

        return Arrays.copyOfRange(result, halflen, x.length + halflen);

    }

    /**
     * Gaussian blur.
     *
     * @param x     the x
     * @param width the width
     * @param sigma the sigma
     * @return the double [ ]
     */
    public static double[] GaussianBlur(int[] x, int width, double sigma) {

        int len = width | 0x01;
        int halflen = len >> 1;

        double[] h = new double[len];

        Gaussian g = new Gaussian(0, sigma);
        double step = 3.0 / halflen;
        double pos = step;

        h[halflen] = g.value(0);
        for (int i = 1; i <= halflen; i++, pos += step) {
            double v = g.value(pos);
            h[halflen + i] = v;
            h[halflen - i] = v;
        }

        h = MathArrays.normalizeArray(h, 1);

        double[] xx = new double[x.length];
        for (int i = 0; i < x.length; i++) {
            xx[i] = x[i];
        }

        double[] result = MathArrays.convolve(xx, h);

        return Arrays.copyOfRange(result, halflen, x.length + halflen);

    }

//	public static double[] robustMean(List<Double> values){
//		List<Double> tempArray = new LinkedList<Double>();
//
//		double meanValue = BioToolbox.mean(values);
//		double stdValue = BioToolbox.std(values);
//
//		double down = 2.0 * stdValue;
//
//		Iterator<Double> iter = values.iterator();
//		while (iter.hasNext()){
//			double cn = iter.next();
//			if (Math.abs(cn - meanValue) < down){
//				tempArray.add(cn);
//			}
//		}
//
//		double m = BioToolbox.mean(values);
//		double s = BioToolbox.std(values);
//		int count = tempArray.size();
//		return new double[]{m/count, Math.sqrt((s-m*m/count)/count)};
//	}

    /**
     * Robust mean.
     *
     * @param values  the values
     * @param maxLoop the max loop
     * @return the double [ ]
     */
    public static double[] robustMean(List<Double> values, int maxLoop) {
        List<Double> tempArray = new LinkedList<Double>(values);
        int count = 0;
        double[] temp;
        for (int i = 0; i < maxLoop; i++) {
            temp = BioToolbox.mean_std(tempArray);
            double meanValue = temp[0];
            double stdValue = temp[1];

            double down = 2.0 * stdValue;
            count = tempArray.size();

            Iterator<Double> iter = tempArray.iterator();
            while (iter.hasNext()) {
                double cn = iter.next();
                if (Math.abs(cn - meanValue) > down) {
                    iter.remove();
                }
            }
            if (tempArray.size() == count) {
                break;
            }
        }

        return BioToolbox.mean_std(tempArray);
    }

    public static double[] getdataSub(double[] values, int start, int stop) {
        double[] result = new double[stop - start];
        for (int i = start, j = 0; i < stop; i++) {
            result[j++] = values[i];
        }
        return result;
    }

    public static double[] robustMean(double[] values, int start, int stop, int maxLoop) {
        List<Double> tempArray = new LinkedList<Double>();
        for (int i = start; i < stop; i++) {
            tempArray.add(values[i]);
        }

        int count;
        double[] temp;
        for (int i = 0; i < maxLoop; i++) {
            temp = BioToolbox.mean_std(tempArray);
            double meanValue = temp[0];
            double stdValue = temp[1];

            double down = 2.0 * stdValue;
            count = tempArray.size();

            Iterator<Double> iter = tempArray.iterator();
            while (iter.hasNext()) {
                double cn = iter.next();
                if (Math.abs(cn - meanValue) > down) {
                    iter.remove();
                }
            }
            if (tempArray.size() == count) {
                break;
            }
        }

        return BioToolbox.mean_std(tempArray);
    }

    public static List<Double> resetIntegralCN(List<Double> values) {
        List<Double> IntegralCN = new ArrayList<Double>();
        IntegralCN.add(0.0);
        double sum = 0;
        for (double v : values) {
            sum += v;
            IntegralCN.add(sum);
        }
        return IntegralCN;
    }

    /**
     * Divide double.
     *
     * @param a the a
     * @param b the b
     * @return the double
     */
    public static double divide(double a, double b) {
        if (b != 0) {
            return a / b;
        } else if ((a == 0) && (b == 0)) {
            return 1;
        }
        return Double.MAX_VALUE;

    }

    public static double meanEx(double[] value) {
        int count = 0;
        double sum = 0;
        for (double x : value) {
            if (x > 0) {
                count++;
                sum += x;
            }
        }
        if (count > 0) {
            return sum / count;
        } else {
            return 0;
        }
    }

    public static double[] medianFiltering(double[] input,int width){
        int size = width | 0x01;
        int length = input.length;
        double[] result = new double[length];
        int halflen = size >> 1;

        for (int i = 0; i < halflen; i++){
            result[i] = input[i];
            result[length - i - 1] = input[length - i - 1];
        }

        for (int i = halflen; i < length - halflen; i++){
            double[] temp =  Arrays.copyOfRange(input, i - halflen, i + halflen);
            Arrays.sort(temp);
            result[i] = temp[halflen];
        }

        return result;
    }
    public static double[] limitFiltering(double[] input,double ratio){
        int length = input.length;
        double[] result = new double[length];

        double temp[] = mean_std(input);
        double up = temp[0] + ratio * temp[1];
        double down = temp[0] - ratio * temp[1];

        for (int i = 0; i < length; i++){
            if (input[i] > up){
                result[i] = up;
            }else if (input[i] < down){
                result[i] = down;
            }else{
                result[i] = input[i];
            }
        }
        return result;
    }
    public static void refreshSegment(Segment seg,double[] chromosome) {
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



}

