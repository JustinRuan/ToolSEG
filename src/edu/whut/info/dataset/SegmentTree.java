package edu.whut.info.dataset;

import java.util.*;

/**
 * Created by Justin on 2016/4/17.
 */
public class SegmentTree {
    private Map<Integer, Double> zMap; // <id, zvalue>
    private Map<Integer, Integer> zPosition;
    private Map<Integer, Integer> parent; // < child, parent>, root is 1.
    public Map<Integer, Segment> segments; //

    public SegmentTree() {
        zMap = new TreeMap<>();
        parent = new HashMap<>();
        zPosition = new TreeMap<>();
        segments = new HashMap<>();
    }

    public void clear(){
        zMap.clear();
        parent.clear();
        zPosition.clear();
        segments.clear();
    }

    public void addZValue(int id, double value) {
        zMap.put(id, value);
    }

    public void addZPosition(int id, int value){
        zPosition.put(id, value);
    }

    public int getZPosition(int id){
        return zPosition.get(id);
    }

    public void addSegment(int id, Segment seg){
        segments.put(id, seg);
    }


    public double getZValueThreshold(){
        double sum = 0;
       for (double value :zMap.values()){
           sum += value;
       }
        return zMap.size() > 0 ? sum / zMap.size() : -1;
    }

    public Set<Map.Entry<Integer, Double>> getZMapEntry(){
        return  zMap.entrySet();
    }

    public void addChild(int childId,int parentId){
        parent.put(childId, parentId);
    }

    public double[] getZValue2Root(int id) {
        if (id == 0) return null;

        List<Integer> path = new LinkedList<>();
        path.add(id);
        getPath(id, path);

        double[] result = new double[path.size()];
        for (int i = 0; i < result.length; i++) {
            result[i] = zMap.get(path.get(i));
        }

        return result;
    }

    public boolean hasChild(int id){
        return parent.values().contains(id);
    }

    public List<Integer> getChilds(int id){
        List<Integer> temp = new LinkedList<>();

        for (Map.Entry<Integer,Integer> elem : parent.entrySet()){
            if (elem.getValue() == id){
                temp.add(elem.getKey());
                if (temp.size() == 2) break;
            }
        }

        return temp;
    }

    public List<Integer> getPath(int childid){
        List<Integer> temp = new LinkedList<>();

        getPath(childid, temp);

        return temp;
    }

    public void getPath(int childId, List<Integer> path) {
        int parentId = parent.get(childId);

        if (parentId > 1) {
            path.add(parentId);
            getPath(parentId, path);
        }else{
            path.add(parentId);
        }
    }

    public int getParent(int childId) {
        return parent.get(childId);
    }

    public double getParentZ(int childId){ // id of root is 1
        if (childId == 1) return Double.MAX_VALUE;

        int id = getParent(childId);
        if (id == 1) return Double.MAX_VALUE;

        return zMap.get(id);
    }

    public Map<Integer, Double> getZMap() {
        return zMap;
    }

    public String getZMapString() {

        StringBuilder result = new StringBuilder();
        result.append('\n');

        LinkedList<Double> zArray = new LinkedList<Double>(zMap.values());
        Collections.sort(zArray);
        for (Map.Entry<Integer, Double> kv : zMap.entrySet()) {
            double z = kv.getValue();
            int pos = kv.getKey();

            List<Integer> childs = new LinkedList<>();
            String child = "( X: X)";
            for (Map.Entry<Integer, Integer> element : parent.entrySet()){
                if (element.getValue() == pos){
                    childs.add(element.getKey());
                    if (childs.size() == 2){
                        child = String.format("(% 4d:% 4d)",childs.get(0),childs.get(1));
                        break;
                    }
                }
            }

            for (int i = 0; i < zArray.size(); i++) {
                if (z == zArray.get(i)) {
                    StringBuilder temp = new StringBuilder();
                    for (int j = 0; j < i; j++) {
                        temp.append("-");
                    }
                    //m_log.info(String.format("[% 8d] z = %f <-%s [%02d]", pos, z, temp.toString(), zArray.size() - i));
                    result.append(String.format("\t\t\t[Seg Id = % 8d >> %s] z = %f <-%s [%02d]", pos,child, z, temp.toString(), zArray.size() - i));
                    result.append('\n');
                    break;
                }
            }
        }
        return result.toString();
    }

//    public String getZMapString() {
//
//        StringBuilder result = new StringBuilder();
//        result.append('\n');
//
//        LinkedList<Double> zArray = new LinkedList<Double>(zMap.values());
//        Collections.sort(zArray);
//        for (Map.Entry<Integer, Double> kv : zMap.entrySet()) {
//            double z = kv.getValue();
//            int id = kv.getKey();
//
//            for (int i = 0; i < zArray.size(); i++) {
//                if (z == zArray.get(i)) {
//                    StringBuilder temp = new StringBuilder();
//                    for (int j = 0; j < i; j++) {
//                        temp.append("-");
//                    }
//                    int begin = segments.get(id).Start();
//                    int end = segments.get(id).End();
//                    int breakpos = zPosition.get(id);
//                    result.append(String.format("\t\t\tSeg [% 8d\t>>% 8d\t>>% 8d] z = %f <-%s [%02d]",
//                            begin, breakpos, end, z, temp.toString(), zArray.size() - i));
//                    result.append('\n');
//                    break;
//                }
//            }
//        }
//        return result.toString();
//    }


}
