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

/**
 * Created by Liu on 2016/6/2.
 */

import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.plot.FastScatterPlot;
import org.jfree.ui.RefineryUtilities;

import javax.swing.*;
import java.awt.*;
import java.util.List;

/**
 * A demo of the fast scatter plot.
 */
public class Sandian extends JFrame {

    /**
     * A constant for the number of items in the sample dataset.
     */
    private static final int COUNT = 500000;

    /**
     * The data.
     */
    private float[][] data = new float[2][COUNT];
    private String name;
    private List<Double> probes;

    /**
     * Creates a new fast scatter plot demo.
     *
     * @param title the frame title.
     */
    public Sandian(final String title, List<Double> probes) {
        super(title);
        this.name = title;
        this.probes = probes;
        populateData(probes);
        final NumberAxis domainAxis = new NumberAxis("Loci");
        domainAxis.setAutoRangeIncludesZero(false);
        final NumberAxis rangeAxis = new NumberAxis("copunumber");
        rangeAxis.setAutoRangeIncludesZero(false);
        final FastScatterPlot plot = new FastScatterPlot(this.data, domainAxis, rangeAxis);
        plot.setPaint(Color.BLACK);
        final JFreeChart chart = new JFreeChart(name, plot);
//        chart.setLegend(null);

        // force aliasing of the rendered content..
        chart.getRenderingHints().put
                (RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);

        final ChartPanel panel = new ChartPanel(chart, true);
        panel.setPreferredSize(new Dimension(500, 270));
        //      panel.setHorizontalZoom(true);
        //    panel.setVerticalZoom(true);
        panel.setMinimumDrawHeight(10);
        panel.setMaximumDrawHeight(2000);
        panel.setMinimumDrawWidth(20);
        panel.setMaximumDrawWidth(2000);

        setContentPane(panel);
    }

    /**
     * Populates the data array with random values.
     */
    private void populateData(List<Double> probes) {
        int count = 0;
        for (int i = 0; i < probes.size(); i++) {
            double value = probes.get(i);
            if (value > 5.5)
                continue;
            final float x = (float) count;
            this.data[0][count] = x;
            this.data[1][count] = (float) value;
            count++;
        }
    }

    public void draw() {

        final Sandian demo = new Sandian(name, probes);
        demo.pack();
        RefineryUtilities.centerFrameOnScreen(demo);
        demo.setVisible(true);
    }
}

