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
 * Created by Liu on 2015/11/26.
 */

import be.ac.ulg.montefiore.run.jahmm.ObservationReal;
import be.ac.ulg.montefiore.run.jahmm.Opdf;

import java.text.NumberFormat;
import java.util.Arrays;
import java.util.Collection;
import java.util.Iterator;

public class NewOpdfGaussian implements Opdf<ObservationReal> {
    private static final long serialVersionUID = 1L;
    private NewGaussianDistribution distribution;

    public NewOpdfGaussian() {
        this.distribution = new NewGaussianDistribution();
    }

    public NewOpdfGaussian(double mean, double variance) {
        this.distribution = new NewGaussianDistribution(mean, variance);
    }

    public double mean() {
        return this.distribution.mean();
    }

    public double variance() {
        return this.distribution.variance();
    }

    public double probability(ObservationReal o) {
        return this.distribution.cumulativeProbability(o.value);
    }

    public ObservationReal generate() {
        return new ObservationReal(this.distribution.generate());
    }

    public void fit(ObservationReal... oa) {
        this.fit(Arrays.asList(oa));
    }

    public void fit(Collection<? extends ObservationReal> co) {
        double[] weights = new double[co.size()];
        Arrays.fill(weights, 1.0D / (double) co.size());
        this.fit(co, weights);
    }

    public void fit(ObservationReal[] o, double[] weights) {
        this.fit(Arrays.asList(o), weights);
    }

    public void fit(Collection<? extends ObservationReal> co, double[] weights) {
        if (!co.isEmpty() && co.size() == weights.length) {
            double mean = 0.0D;
            int i = 0;

            ObservationReal variance;
            for (Iterator var7 = co.iterator(); var7.hasNext(); mean += variance.value * weights[i++]) {
                variance = (ObservationReal) var7.next();
            }

            double var12 = 0.0D;
            i = 0;

            double d;
            for (Iterator var9 = co.iterator(); var9.hasNext(); var12 += d * d * weights[i++]) {
                ObservationReal o = (ObservationReal) var9.next();
                d = o.value - mean;
            }

            this.distribution = new NewGaussianDistribution(mean, var12);
        } else {
            throw new IllegalArgumentException();
        }
    }

    public NewOpdfGaussian clone() {
        try {
            return (NewOpdfGaussian) super.clone();
        } catch (CloneNotSupportedException var2) {
            throw new AssertionError(var2);
        }
    }

    public String toString() {
        return this.toString(NumberFormat.getInstance());
    }

    public String toString(NumberFormat numberFormat) {
        return "Gaussian distribution --- Mean: " + numberFormat.format(this.distribution.mean()) + " Variance " + numberFormat.format(this.distribution.variance());
    }
}
