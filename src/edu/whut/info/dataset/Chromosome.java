/**
 *
 */
package edu.whut.info.dataset;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;

/**
 * @author Justin
 */
public class Chromosome implements Serializable {

    private static final long serialVersionUID = -3488836842658361449L;
    public String sampleId;
    public short chrId;
    public List<Double> probes;
    public List<Long> changepoints;
    public Double baseline;

    public Chromosome() {
        probes = new ArrayList<>();
        changepoints = new ArrayList<>();
        baseline = -1.0;
    }

    public void setChrId(short chrId) {
        this.chrId = chrId;
    }

    public void setSampleId(String sampleId) {
        this.sampleId = sampleId;
    }

    public void clear() {
        probes.clear();
    }

    public int getLength() {
        return probes.size();
    }

    public void setProbes(List<Double> ps) {
        probes = ps;
    }

}
