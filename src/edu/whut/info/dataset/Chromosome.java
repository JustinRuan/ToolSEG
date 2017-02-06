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

    public void sampling(){
        List<Double> newProbes = new ArrayList<>();
        List<Long> newCpts = new ArrayList<>();
        for (int i = 0; i < probes.size();i++){
            if (i % 100 == 0){
                newProbes.add(probes.get(i));
            }
        }
        for (Long pos : changepoints){
            if (pos > 100){
                if (newCpts.size() > 0){
                    Long last = newCpts.get(newCpts.size()-1);
                    Long current = pos /100;

                    if (current - last > 3){
                        newCpts.add(current);
                    }
                }else{
                    newCpts.add(pos/100);
                }
            }

        }
        probes = newProbes;
        changepoints = newCpts;
    }


}
