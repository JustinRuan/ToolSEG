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
