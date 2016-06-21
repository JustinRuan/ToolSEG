package edu.whut.info.dataset;

import be.ac.ulg.montefiore.run.jahmm.ObservationReal;

/**
 * Created by Liu on 2015/11/19.
 */


public class ObservationRealEx extends ObservationReal implements org.apache.commons.math3.ml.clustering.Clusterable {

    public ObservationRealEx(double value) {
        super(value);
    }

    @Override
    public double[] getPoint() {
        return new double[]{this.value};
    }
}
