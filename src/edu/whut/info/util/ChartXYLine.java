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
package edu.whut.info.util;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartUtilities;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.data.xy.XYDataset;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;

import java.awt.*;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;

/**
 * Created by Justin on 2015/1/21.
 */
public class ChartXYLine {
    public JFreeChart jfreechart;
    private XYDataset dataset;

    public ChartXYLine() {
        XYSeriesCollection xySeriesCollection = new XYSeriesCollection();
        //Histogram h, String name, String result
//        XYSeries xyseries = new XYSeries(name);
//
//        for (int i = 1; i < h.getBin(); i++) {
//            xyseries.add(h.getThresh(i), h.getFreq(i));
//        }
//
//        xySeriesCollection.addSeries(xyseries);

        dataset = xySeriesCollection;
//        createChart(name, result);
    }

    public void addLine(XYSeries xyline) {
        ((XYSeriesCollection) dataset).addSeries(xyline);
    }

    public void saveAsFile(String outputPath,
                           int weight, int height) {
        FileOutputStream out = null;
        try {
            File outFile = new File(outputPath);
            if (!outFile.getParentFile().exists()) {
                outFile.getParentFile().mkdirs();
            }
            out = new FileOutputStream(outputPath);

            // 保存为PNG
            ChartUtilities.writeChartAsPNG(out, jfreechart, weight, height);
            // 保存为JPEG
            // ChartUtilities.writeChartAsJPEG(out, chart, weight, height);
            out.flush();
        } catch (FileNotFoundException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();
        } finally {
            if (out != null) {
                try {
                    out.close();
                } catch (IOException e) {
                    // do nothing
                }
            }
        }
    }

    public void createChart(String result, String xStr, String yStr) {
        jfreechart = ChartFactory.createXYLineChart(result, xStr, yStr,
                dataset,
                PlotOrientation.VERTICAL,
                false,
                true,
                false);
        jfreechart.setBackgroundPaint(Color.white);

    }
}
