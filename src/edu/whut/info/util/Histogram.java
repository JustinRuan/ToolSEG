package edu.whut.info.util;

import org.apache.commons.math3.util.MathArrays;
import org.jfree.data.xy.XYSeries;

import java.io.File;
import java.text.SimpleDateFormat;
import java.util.Arrays;
import java.util.Date;

/**
 * The type Histogram.
 */
public class Histogram {

    private int[] freq;   // freq[i] = # occurences of value i
    private double step;
    private double[] thresh;
    private int bin;

    private int count;
    private double down;
    private double up;

    /**
     * Instantiates a new Histogram.
     *
     * @param N the number of bins
     * @param Q the lower limit of the range
     * @param P the upper limit of the range
     */
// Create a new histogram.
    public Histogram(int N, double Q, double P) {

        bin = N;
        // the first is x<Q
        //the last is x>P
        freq = new int[N + 2]; // [x<Q,[Q:(P-Q)/N:P],x>P]

        thresh = new double[N + 1];

        step = (P - Q) / (double) N;
        down = Q;
        up = P;

        thresh[0] = down;
        for (int i = 1; i < N + 1; i++) {
            thresh[i] = thresh[i - 1] + step;
        }
        count = 0;

    }

    /**
     * Gets the number of bin.
     *
     * @return the bin
     */
    public int getBin() {
        return bin;
    }

    /**
     * Gets freq at the index.
     *
     * @param index the index
     * @return the freq
     */
    public int getFreq(int index) {
        return freq[index];
    }

    /**
     * Gets thresh.
     *
     * @param index the index
     * @return the thresh
     */
    public double getThresh(int index) {
        return thresh[index];
    }

    /**
     * Get freq array.
     *
     * @return the int [ ]
     */
    public int[] getFreqArray() {
        return Arrays.copyOfRange(freq, 1, bin + 1);
    }

    /**
     * Gets count.
     *
     * @return the count
     */
    public int getCount() {
        return count;
    }

    public void clear() {
        Arrays.fill(freq, 0);
    }

    public double percent(double position) {
        int sum = 0;
        for (int i = 1; i < thresh.length - 1; i++) {
            if (position > thresh[i]) {
                sum += freq[i];
            }
        }
        return ((double) sum) / count;
    }

    /**
     * Get peak range.
     *
     * @return the double [ ]
     */
    public double getPeakRange(int width) {

        double[] newFreq = BioToolbox.GaussianBlur(freq, width, 1);

        int maxIndex = 0;
        double peak = 0;

        for (int i = 1; i < newFreq.length - 1; i++) {
            if (peak < newFreq[i]) {
                peak = newFreq[i];
                maxIndex = i;
            }
        }

        return (thresh[maxIndex - 1] + thresh[maxIndex]) / 2;
    }

    /**
     * Get peak range by percent.
     *
     * @param peakPercent the ratio of the area under the freq curve at the peak point, the range is [0,1]
     * @return the double [ ]
     */
    public double[] getPeakRangeByPercent(double peakPercent) {

        double[] newFreq = BioToolbox.GaussianBlur(freq, 5, 1);

        int maxIndex = 0;
        double peak = 0;
        for (int i = 0; i < newFreq.length; i++) {
            if (peak < newFreq[i]) {
                peak = newFreq[i];
                maxIndex = i;
            }
        }

        double[] p1 = MathArrays.normalizeArray(newFreq, 1);

        int front = maxIndex - 1;
        int rear = maxIndex + 1;
        double sum = p1[maxIndex];

        while (sum < peakPercent) {
            if ((front == 0) && (rear == p1.length - 1)) break;

            if (front == 0) {
                sum += p1[rear];
                rear++;
            } else if (rear == p1.length - 1) {
                sum += p1[front];
                front--;
            } else if (p1[front] > p1[rear]) {
                sum += p1[front];
                front--;
            } else {
                sum += p1[rear];
                rear++;
            }
        }
        double[] th = new double[]{thresh[front - 1], thresh[rear]};

        return th;
    }

    /**
     * Add data array.
     *
     * @param values the values
     */
    public void addDataArray(Double[] values) {
        for (Double d : values) {
            addDataPoint(d);
        }
    }

    /**
     * Add data point.
     *
     * @param value the value
     */
    public void addDataPoint(double value) {
        count++;

        if (value < down) {
            freq[0]++;
            return;
        }

        int i = (int) ((value - down) / step);

        if (i < bin) {
            freq[i + 1]++;
        } else {
            freq[bin + 1]++;
        }

    }

    public void addDataPoint(double value, int c) {
        count += c;

        if (value < down) {
            freq[0] += c;
            return;
        }

        int i = (int) ((value - down) / step);

        if (i < bin) {
            freq[i + 1] += c;
        } else {
            freq[bin + 1] += c;
        }

    }

    /**
     * To string.
     *
     * @param name the name
     * @return the string
     */
    public String toString(String name) {
        int spacePosition = name.indexOf(" ");
        String rowText;
        if (spacePosition > 0) {
            rowText = name.substring(0, spacePosition) + "...";
        } else {
            rowText = name + "...";
        }


        StringBuilder s = new StringBuilder();
        s.append(String.format("##### %s ; count = %d #####\n", name, count));

        s.append(String.format("%s Histogram %4d:\t [%s : %.5f]\t freq = %d \n",
                rowText, 0, "  -Inf ", thresh[0], freq[0]));

        for (int i = 1; i < bin + 1; i++) {
            if (freq[i] > 0) {
                s.append(String.format("%s Histogram %4d:\t [%.5f : %.5f]\t freq = %d \n",
                        rowText, i, thresh[i - 1], thresh[i], freq[i]));
            }

        }
        s.append(String.format("%s Histogram %4d:\t [%.5f : %s]\t freq = %d \n",
                rowText, bin + 1, thresh[bin], "  +Inf ", freq[bin + 1]));
        return s.toString();
    }

    public void saveToPicture(String title, String sampleId, String result, String xStr, String yStr) {
        String temp = String.format("%s: %s", sampleId, result);
        XYSeries xyseries = new XYSeries("Rou");

        for (int i = 1; i < this.getBin(); i++) {
            xyseries.add(this.getThresh(i), this.getFreq(i));
        }
        ChartXYLine cxy = new ChartXYLine();
        cxy.addLine(xyseries);
        cxy.createChart(temp, xStr, yStr);

        Date date = new Date();
        SimpleDateFormat df = new SimpleDateFormat("yyyyMMdd_HHmmss");
        String filename = String.format("..%1$sData%1$s%2$s_%3$s_%4$s.png",
                File.separator, title, sampleId, df.format(date));

        cxy.saveAsFile(filename, 600, 400);
    }


    /**
     * Get thresh.
     *
     * @param r the r
     * @return the double [ ]
     */
    public double[] getThresh(double r) {
        double[] th = new double[2];

        //int border = (int) (r * count) - freq[bin+1] - freq[0];
        int border = (int) (r * count);

        if (border < 0) {
            System.out.println("The Thresh of Histogram is Error \n");
            th[0] = down;
            th[1] = 2 * up;
            return th;
        }


        int[] fq = new int[bin];
        System.arraycopy(freq, 1, fq, 0, bin);

        Arrays.sort(fq);

        int threshValue = 0;
        for (int v : fq) {
            border -= v;
            if (border < 0) {
                threshValue = v;
                break;
            }
        }

        int i;
        for (i = 1; i < bin; i++) {
            if (freq[i] > threshValue) {
                th[0] = thresh[i];
                break;
            }
        }
        for (i = bin - 1; i > 1; i--) {
            if (freq[i] > threshValue) {
                th[1] = thresh[i + 1];
                break;
            }
        }

        if (th[1] > th[0]) return th;
        else return null;
    }

}
//        for (int i = maxIndex - 1; i > 0; i--) { //��ҪҪ�Ľ������ܿ絽��һ������ȥ��
//            if (Math.abs(newFreq[i] / newFreq[maxIndex] - 1) < peakRatio) {
//                th[0] = thresh[i - 1];
//            }
//        }
//        for (int i = maxIndex + 1; i < bin; i++) {
//            if (Math.abs(newFreq[i] / newFreq[maxIndex] - 1) < peakRatio) {
//                th[1] = thresh[i];
//            }
//        }