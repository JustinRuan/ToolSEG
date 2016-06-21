package edu.whut.info.util;

/**
 * Created by Liu on 2015/11/26.
 */

import be.ac.ulg.montefiore.run.distributions.RandomDistribution;
import org.apache.commons.math3.special.Erf;
import org.apache.commons.math3.util.FastMath;

import java.util.Random;

public class NewGaussianDistribution implements RandomDistribution {
    private static final Random randomGenerator = new Random();
    private static final long serialVersionUID = 9127329839769283975L;
    private static final double SQRT2 = FastMath.sqrt(2.0D);
    private double mean;
    private double deviation;
    private double variance;

    public NewGaussianDistribution() {
        this(0.0D, 1.0D);
    }

    public NewGaussianDistribution(double mean, double variance) {
        if (variance <= 0.0D) {
            throw new IllegalArgumentException("Variance must be positive");
        } else {
            this.mean = mean;
            this.variance = variance;
            this.deviation = Math.sqrt(variance);
        }
    }

    public double mean() {
        return this.mean;
    }

    public double variance() {
        return this.variance;
    }

    public double generate() {
        return randomGenerator.nextGaussian() * this.deviation + this.mean;
    }

    public double probability(double n) {
        double expArg = -0.5D * (n - this.mean) * (n - this.mean) / this.variance;
        return Math.pow(6.283185307179586D * this.variance, -0.5D) * Math.exp(expArg);
    }

    public double cumulativeProbability(double x) {
        double dev = x - this.mean;
        double p = FastMath.abs(dev) > 40.0D * this.deviation ? (dev < 0.0D ? 0.0D : 1.0D) : 0.5D * (1.0D + Erf.erf(dev / (this.deviation * SQRT2)));
        if (p > 0.5)
            p = 1 - p;
        return p;

    }
}
