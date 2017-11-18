package edu.whut.info.dataset;


import edu.whut.info.util.BioToolbox;

import java.io.File;
import java.sql.*;
import java.util.*;
import java.util.logging.Logger;

/**
 * @author Liu
 * @version 1.0
 * @created 10-10-2015 14:57:00
 */
public class ChromoGenerator {

    //  private Map<Long, ProbeSet> segmentTemplate;
    private String filepath;
    private List<Integer> lengthSet;
    private List<Double> cnSet;
    private List<Double> stdSet;
    private Logger m_log;

    public ChromoGenerator(String filepath) {
        m_log = Logger.getLogger("segment");
        this.filepath = filepath;
        // segmentTemplate = new TreeMap<>();
        lengthSet = new ArrayList<>();
        cnSet = new ArrayList<>();
        stdSet = new ArrayList<>();
    }

    public ChromoGenerator() {
        m_log = Logger.getLogger("segment");
        this.filepath = filepath;
        //segmentTemplate = new TreeMap<>();
        lengthSet = new ArrayList<>();
        cnSet = new ArrayList<>();
        stdSet = new ArrayList<>();
    }

    public void setLengthSet(List<Integer> lengthSet) {
        this.lengthSet = lengthSet;
    }

    public void setCnSet(List<Double> cnSet) {
        this.cnSet = cnSet;
    }

    public void setStdSet(List<Double> stdSet) {
        this.stdSet = stdSet;
    }

    public void setFilepath(String filepath) {
        this.filepath = filepath;
    }

    public Chromosome chromoGeneratorByRandom(double normalPercent) {
        m_log.info("generated segment: ");
        Chromosome chro = new Chromosome();
        long startIndex = 1;
        Random random = new Random(47);
        Random randMean = new Random(234);
        ArrayList<Double> values1 = new ArrayList<>();
        for (int i = 0; i < lengthSet.size(); i++) {
            double std = stdSet.get(i);
            ArrayList<Double> values = new ArrayList<>();
            int length = lengthSet.get(i);
            double cn = cnSet.get(i);
            cn = cn * (1 - normalPercent) + 2 * normalPercent;

            //double redund=0;
            double copynumber;
            double down = Math.max(0, cn - 3.5 * std);
            double up = cn + 3.5 * std;
            for (int j = 0; j < length; j++) {
                //   ProbeSet p1 = new ProbeSet();
                do {
                    copynumber = std * random.nextGaussian() + cn;
                } while ((copynumber > up) || (copynumber < down));

                // p1.halfCN = copynumber;
                //  p1.Loci=startIndex;
                chro.probes.add(copynumber);
                startIndex++;
                values.add(copynumber);
                values1.add(copynumber);
            }
            chro.changepoints.add(startIndex);
            m_log.info(String.format("the original %2d segment:\t length=%8d\t mean=%.4f\t std=%.4f\t", i, length, BioToolbox.mean(values), BioToolbox.std(values)));
        }
//        System.out.println("copy number sequence: ");
//        int i=0;
//        for (double d:values1){
//            System.out.print(String.format("%.3f",d)+"   ");
//            i++;
//            if(i%10==0)
//                System.out.println();
//        }
        chro.changepoints.remove(chro.changepoints.get(chro.changepoints.size() - 1));

        return chro;
    }

    public Chromosome chromoGeneratorByAllTemp(double normalPercent) {
        Chromosome chro = new Chromosome();
        Map<Integer, List<Double>> rawdata = new HashMap<>();
        for (int i = 1; i <= 4; i++)
            rawdata.put(i, new ArrayList<Double>());
        List<String> datas = new ArrayList<>();
        datas.add(String.format(".%1$sdata%1$ssourcedata.db", File.separator));
        datas.add(String.format(".%1$sdata%1$schro1.db", File.separator));
        for (String s : datas)
            merge(s, rawdata);
        int segNum = 0;
        long startIndex = 1;
        for (int i = 0; i < lengthSet.size(); i++) {
            int length = lengthSet.get(i);
            double cn = cnSet.get(i);
            int cn1 = (int) cn;
            List<Double> values = rawdata.get(cn1);
            int size = values.size();
            Random r = new Random(23);
            List<Double> temp=new ArrayList<>();
            for (int j = 0; j < length; j++) {
                int lg = j % size;
                double copyNumberHalf = values.get(lg) * (1 - normalPercent) + 2 * normalPercent;
                chro.probes.add(copyNumberHalf);
                temp.add(copyNumberHalf);
                startIndex++;
            }
            chro.changepoints.add(startIndex);
            segNum++;
            m_log.info(String.format("the original %2d segment:\t length=%8d\t mean=%.4f\t std=%.4f\t", segNum, length, BioToolbox.mean(temp), BioToolbox.std(temp)));
        }
        chro.changepoints.remove(chro.changepoints.get(chro.changepoints.size() - 1));
        return chro;
    }

    private void merge(String fileName, Map<Integer, List<Double>> rawdata) {
        try {
            Class.forName("org.sqlite.JDBC");
            String ConnectionString = String.format("jdbc:sqlite:%s", fileName);
            Connection c = null;
            c = DriverManager.getConnection(ConnectionString);
            Statement stmt = null;
            stmt = c.createStatement();
            for (int cn = 1; cn <= 4; cn++) {
                ResultSet rs;
                ArrayList<Double> values = new ArrayList<>();
                String sql = String.format("Select loci,Value,copynumber from probesets WHERE Value =%s ORDER by loci", cn * 1.0);
                //    String sql=String.format("Select ID_REF,VALUE,CHROMOSOME,POSITION,\"SMOOTH SIGNAL\" as" +
                //             " CopyNumber from \"GSM721892-5637\" WHERE VALUE =%d"+
                //                "  And CHROMOSOME = 15 And substr(ID_REF,0,4)=\"SNP\" ORDER by POSITION",cn);
                rs = stmt.executeQuery(sql);
                // rs=stmt.executeQuery("Select exp(\"LOG2 RATIO\") as CopyNumber from \"GSM721892-5637\" WHERE CHROMOSOME = 1 And substr(ID_REF,0,4)=\"SNP\" AND VALUE =4\n" +
                //         "ORDER by POSITION");
                while (rs.next()) {
                    double lRR = rs.getDouble("copynumber");
                    values.add(lRR);
                }
                rawdata.get(cn).addAll(values);
                //m_log.info(String.format("the original %2d sourcedata:\t length=%8d\t mean=%.4f\t std=%.4f\t", cn, values.size(), BioToolbox.mean(values), BioToolbox.std(values)));
                //rawdata.put(cn, values);
            }
        } catch (ClassNotFoundException e) {
            e.printStackTrace();
        } catch (SQLException e) {
            e.printStackTrace();
        }
    }

    public Chromosome chromoGeneratorByShift(double normalPercent) {
        List<Double> template = new ArrayList<>();
        readTemplate(template);
        Chromosome chro = new Chromosome();
        long startIndex = 1;
        //  ArrayList<Long> Index = new ArrayList<Long>(segmentTemplate.keySet());
        //Set<Long> segmentIndex=new TreeSet<>();
        for (int i = 0; i < lengthSet.size(); i++) {
            ArrayList<Double> values = new ArrayList<>();
            int length = lengthSet.get(i);
            double cn = cnSet.get(i);
            double redund = 0;
            if (cn < 2)
                redund = -1 * (2 - cn);
            else
                redund = cn - 2;
            for (int j = 0; j < length; j++) {
                int index = j % template.size();
                double p = template.get(index);
                p = (Math.abs(p * 2 + redund)) * (1 - normalPercent) + 2 * normalPercent;
                chro.probes.add(p);
                values.add(p);
                startIndex++;
            }
            m_log.info(String.format("the original %2d segment:\t length=%8d\t mean=%.4f\t std=%.4f\t", i, length, BioToolbox.mean(values), BioToolbox.std(values)));
            chro.changepoints.add(startIndex);
        }
        chro.changepoints.remove(chro.changepoints.get(chro.changepoints.size() - 1));
        return chro;
    }

    private void readTemplate(List<Double> template) {
        try {
            ArrayList<Double> values = new ArrayList<>();
            Class.forName("org.sqlite.JDBC");
            String filepath = String.format(".%1$sdata%1$sdatasource.db", File.separator);
            String ConnectionString = String.format("jdbc:sqlite:%s", filepath);
            Connection c = DriverManager.getConnection(ConnectionString);
            Statement stmt = c.createStatement();
            ResultSet rs;
            rs = stmt.executeQuery("SELECT loci,copynumber,Allele_A,Allele_B FROM chromosome");
            while (rs.next()) {
                Long loci = rs.getLong("loci");
                double cp = rs.getDouble("copynumber");
                double Allele_A = rs.getDouble("Allele_A");
                double Allele_B = rs.getDouble("Allele_B");
                template.add(cp);
            }
            rs.close();
            stmt.close();
            c.close();
            //  System.out.println("mean and std: " + BioToolbox.mean(values) + " " + BioToolbox.std(values));
        } catch (Exception e) {
            System.err.println(e.getClass().getName() + ": " + e.getMessage());
            System.exit(0);
        }
    }
}