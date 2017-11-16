package edu.whut.info.ui;

import edu.whut.info.core.*;
import edu.whut.info.dataset.Chromosome;
import edu.whut.info.dataset.Segment;
import edu.whut.info.util.BioLogger;
import edu.whut.info.util.BioToolbox;
import edu.whut.info.util.Controller;
import edu.whut.info.util.Sandian;

import javax.imageio.ImageIO;
import javax.swing.*;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.image.BufferedImage;
import java.io.*;
import java.text.SimpleDateFormat;
import java.util.*;
import java.util.List;
import java.util.logging.Logger;

/**
 * Created by Justin on 2016/3/18.
 */
public class Mainform {
    public int segmentMethod;
    public int nums = -1;
    public short chrId = 1;
    public int dataGenerate;
    public int workModel = 0;
    public ArrayList<Chromosome> chros;
    public SegmentCutter cutterAlgorithm;
    List<Integer> lengths;
    List<Double> cn;
    List<Double> variance;
    Controller myCore;
    double centerProb;
    double penaltyR;
    double gamap;
    String filepath = "";
    private JFrame m_frame;
    private JRadioButton byRandomRadioButton;
    private JRadioButton byNeutralTemplateRadioButton;
    private JRadioButton byMultiTemplateRadioButton;
    private JSpinner contaminationv;
    private JButton saveButton;
    private JButton showButton;
    private JButton generateButton;
    private JFormattedTextField segmentText;
    private JTabbedPane SegmentMethod;
    private JButton btnTestButton;
    private JButton segment;
    private JTextField inputfile;
    private JButton output;
    private JTextField txtOutputfile;
    private JButton input;
    private JCheckBox outliers;
    private JPanel potential;
    private JPanel MainPanel;
    private JPanel FastPCF;
    private JPanel CLT;
    private JPanel Lasso;
    private JPanel CBS;
    private JPanel PCF;
    private JTextField center;
    private JTextField gama;
    private JTextField cltlengthtxt;
    private JTextField lamdatxt;
    private JTextField maxItertxt;
    private JTextField cbsLengthtxt;
    private JTextField pgamatxt;
    private JTabbedPane tabbedPane1;
    private JTextField cltsteptxt;
    private JTextField cltpvaluetxt;
    private JButton btnSegment;
    private JTextField segmentlengthtxt;
    private JTextField cnvaluetxt;
    private JTextField variancetxt;
    private JTextField percenttxt;
    private JTextField tolerancetxt;
    private JButton btnTestAll;
    private JTextField cbsTimesText;
    private JTextField cbsStepText;
    private JRadioButton byTemplateShiftRadioButton;
    private JCheckBox chkIsTestdata;
    private JButton clearButton;
    private JTextField txtbctlamda;
    private JButton showReultButton;
    private JLabel cltLengtxt;
    private JPanel DBS;
    private JTextField BCLTLengtxt;
    private JTextField BCLTpvaluetxt;
    private JTextField outliertxt;
    private JRadioButton bypassRadioButton;
    private JRadioButton log2RadioButton;
    private JRadioButton power2RadioButton;
    private JTextField BCLTLambdatxt;
    private JCheckBox chkIsMultiFiles;
    private JButton btnClear;
    private JButton btnLoad;
    private JButton BtnCount;
    private JLabel LblInfo;
    private Logger m_log;

    private int transformMethod;

    public Mainform() {
        Date date = new Date();
        SimpleDateFormat df = new SimpleDateFormat("yyyyMMdd_HHmmss");
        //String filename = String.format(".%sResult_%s.log", File.separator, df.format(date));
        String filename = String.format("Result_%s.log", df.format(date));
        String path = String.format(".%sResult", File.separator);
        BioLogger log = new BioLogger(path, filename);
        m_log = Logger.getLogger("segment");
        m_frame = new JFrame("ToolSEG 1.2");
        cltsteptxt.setText("4");
        cltlengthtxt.setText("256");
        cbsLengthtxt.setText("100");
        cbsStepText.setText("8");
        cltpvaluetxt.setText("0.05");
        m_frame.setContentPane(MainPanel);
        m_frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        m_frame.pack();
        m_frame.setVisible(true);
        lengths = new ArrayList<>();
        cn = new ArrayList<>();
        variance = new ArrayList<>();
        myCore = new Controller();
        segmentMethod = 0;
        chros = new ArrayList<>();
        percenttxt.setText("0.5");
        pgamatxt.setText("20");
        lamdatxt.setText("20");
        tolerancetxt.setText("0.01");
        maxItertxt.setText("20");
        center.setText("8");
        gama.setText("20");

        //DBS params
        BCLTLengtxt.setText("50");
        BCLTpvaluetxt.setText("0.05");
        BCLTLambdatxt.setText("0.02");

        //preprocess
        outliertxt.setText("2.5");


//        segmentlengthtxt.setText("20000,4000,500,10000,5000");
//        cnvaluetxt.setText("1,2,1,2,3");
//        variancetxt.setText("0.2,0.2,0.3,0.3,0.2");
        cnvaluetxt.setText("1,3,2,4,1,3,2");
        variancetxt.setText("0.4,0.4,0.4,0.4,0.4,0.4,0.4");
        segmentlengthtxt.setText("3000,2000,3000,5000,4000,5000,4000");


        System.out.println();
        btnSegment.setEnabled(false);
        btnSegment.setToolTipText("before segment you need generate data or choose input file");
        byRandomRadioButton.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                dataGenerate = 0;
            }
        });
        byMultiTemplateRadioButton.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                dataGenerate = 1;
            }
        });
        byTemplateShiftRadioButton.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                dataGenerate = 2;
            }
        });
        generateButton.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                readLenCn();
                myCore.generator.setCnSet(cn);
                myCore.generator.setLengthSet(lengths);
                myCore.generator.setStdSet(variance);
                double normalpercent = Double.valueOf(percenttxt.getText());
                Chromosome chro = null;
                m_log.info(String.format("chro num %d",chrId));
                switch (dataGenerate) {
                    case 0:
                        chro = myCore.generator.chromoGeneratorByRandom(normalpercent);
                        break;
                    case 1:
                        chro = myCore.generator.chromoGeneratorByAllTemp(normalpercent);
                        break;
                    case 2:
                        chro = myCore.generator.chromoGeneratorByShift(normalpercent);
                }
                chro.setChrId(chrId);
                chros.add(chro);
                btnSegment.setEnabled(true);
                nums++;
                chrId++;
            }
        });
        btnSegment.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                switch (segmentMethod) {
                    case 4:
                        centerProb = Double.valueOf(center.getText());
                        cutterAlgorithm = new HMMSegment(centerProb);
                        break;
                    case 2:
                        penaltyR = Double.valueOf(gama.getText());
                        cutterAlgorithm = new FastPCF(200, 0.8, penaltyR);
                        break;
                    case 3:
                        int length = Integer.valueOf(cltlengthtxt.getText());
                        int step = Integer.valueOf(cltsteptxt.getText());
                        double p = Double.valueOf(cltpvaluetxt.getText());
                        cutterAlgorithm = new BACOMSegment(p, length, step);
                        break;
                    case 5:
                        double lamda = Double.valueOf(lamdatxt.getText());
                        double tol = Double.valueOf(tolerancetxt.getText());
                        int maxIter = Integer.valueOf(maxItertxt.getText());
                        cutterAlgorithm = new LassoSegment(lamda, tol, maxIter);
                        break;
                    case 0:
                        int cbslength = Integer.valueOf(cbsLengthtxt.getText());
                        int cbsstep = Integer.valueOf(cbsStepText.getText());
                        cutterAlgorithm = new CBSSegment(cbslength, cbsstep);
                        break;
                    case 1:
                        gamap = Double.valueOf(pgamatxt.getText());
                        cutterAlgorithm = new PCFSegment(gamap);
                        break;
                    case 6:
                        length = Integer.valueOf(BCLTLengtxt.getText());
                        p = Double.valueOf(BCLTpvaluetxt.getText());
                        double lambda = Double.valueOf(BCLTLambdatxt.getText());
                        cutterAlgorithm = new DBS(p, length, lambda);
                        break;
                    case 7:
                        gamap = Double.valueOf(pgamatxt.getText());
                        cutterAlgorithm = new PCFM(gamap);
                        break;
                }

                double ratio = 0;
                if (outliers.isSelected()) {
                    ratio = Double.valueOf(outliertxt.getText());

                }
                myCore.CnSegment.cutter = cutterAlgorithm;
                switch (workModel) {
                    case 0:
                        myCore.CnSegment.splitChromosome(chros, ratio, transformMethod, chkIsTestdata.isSelected());
                        break;
                    case 1:
//                        try {
//                            boolean flag = chkIsMultiFiles.isSelected();
//                            if (!flag) {
//                                ArrayList<Chromosome> tmpData = readCNtxt(filepath);
//                                chros.addAll(tmpData);
//                            } else {
//                                ArrayList<String> files = readFilePath(filepath);
//                                for (String path : files) {
//                                    ArrayList<Chromosome> tmpData = readCNtxt(path);
//                                    chros.addAll(tmpData);
//                                }
//                            }
//                        } catch (IOException e1) {
//                            e1.printStackTrace();
//                        }

                        myCore.CnSegment.splitChromosome(chros, ratio, transformMethod, chkIsTestdata.isSelected());
                }
            }
        });
        input.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                JFileChooser jfc = new JFileChooser(".");
                jfc.setFileSelectionMode(JFileChooser.FILES_AND_DIRECTORIES);
                //jfc.addChoosableFileFilter(new TxtFileFilter());
                jfc.showDialog(new JLabel(), "Select");
                File file = jfc.getSelectedFile();
                filepath = file.getAbsolutePath();

                inputfile.setText(file.getName());
                if (file.isDirectory()) {
                    System.out.println("Directory:" + file.getAbsolutePath());
                } else if (file.isFile()) {
                    System.out.println("File:" + file.getAbsolutePath());
                }
                btnSegment.setEnabled(true);
                System.out.println(jfc.getSelectedFile().getName());
            }
        });
        SegmentMethod.addChangeListener(new ChangeListener() {
            @Override
            public void stateChanged(ChangeEvent e) {
                segmentMethod = SegmentMethod.getSelectedIndex();
            }
        });
        btnTestAll.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                ArrayList<SegmentCutter> cutters = new ArrayList<SegmentCutter>();
                centerProb = Double.valueOf(center.getText());
                cutterAlgorithm = new HMMSegment(centerProb);
                cutters.add(cutterAlgorithm);
                penaltyR = Double.valueOf(gama.getText());
                cutterAlgorithm = new FastPCF(200, 0.8, penaltyR);
                cutters.add(cutterAlgorithm);
                int minlength = Integer.valueOf(cltlengthtxt.getText());
                int minStep = Integer.valueOf(cltsteptxt.getText());
                double pvalue = Double.valueOf(cltpvaluetxt.getText());
                // double bcltlamda=Double.valueOf(txtbctlamda.getText());
                cutterAlgorithm = new BACOMSegment(pvalue, minlength, minStep);
                cutters.add(cutterAlgorithm);
                double lamda = Double.valueOf(lamdatxt.getText());
                double tol = Double.valueOf(tolerancetxt.getText());
                int maxIter = Integer.valueOf(maxItertxt.getText());
                cutterAlgorithm = new LassoSegment(lamda, tol, maxIter);
                cutters.add(cutterAlgorithm);
                int cbslength = Integer.valueOf(cbsLengthtxt.getText());
                int cbsstep = Integer.valueOf(cbsStepText.getText());
                cutterAlgorithm = new CBSSegment(cbslength, cbsstep);
                cutters.add(cutterAlgorithm);
                gamap = Double.valueOf(pgamatxt.getText());
                cutterAlgorithm = new PCFSegment(gamap);
                ArrayList<Chromosome> data = null;
                switch (workModel) {
                    case 0:
                        data = chros;
                        break;
                    case 1:
                        try {
                            data = readCNtxt(filepath);
                        } catch (IOException e1) {
                            e1.printStackTrace();
                        }
                }
                Date date = new Date();
                SimpleDateFormat df = new SimpleDateFormat("yyyyMMdd_HHmmss");
                String filename = String.format("data%1$sResult_%2$s.txt",
                        File.separator, df.format(date));
                List<List<Set<Segment>>> result = new ArrayList<List<Set<Segment>>>();
                File output = new File(filename);
                try {
                    output.createNewFile(); //
                    BufferedWriter out = new BufferedWriter(new FileWriter(output));
                    out.write("Method\t\t numberOfSeg\t\t foundSeg\t\t Tolerance\t\t truePoint\t\t segTime\t\t \r\n");

                    int width = 0;
                    if (outliers.isSelected()) {
                        width = Integer.valueOf(outliertxt.getText());

                    }

                    for (SegmentCutter cutter : cutters) {
                        myCore.CnSegment.setCutter(cutter);
                        myCore.CnSegment.splitChromosome(data, width, transformMethod, chkIsTestdata.isSelected());
                        List<Set<Segment>> temp = myCore.CnSegment.resultclone();
                        // temp=myCore.CnSegment.getResult();
                        result.add(temp);
                        //myCore.CnSegment.getResult().clear();
                    }
                    out.flush();
                    out.close();
                } catch (IOException e1) {
                    e1.printStackTrace();
                }
                for (int i = 0; i < data.size(); i++) {
                    List<Set<Segment>> onechromBy = new ArrayList<Set<Segment>>();
                    onechromBy.add(result.get(0).get(i));
                    onechromBy.add(result.get(1).get(i));
                    onechromBy.add(result.get(2).get(i));
                    onechromBy.add(result.get(3).get(i));
                    onechromBy.add(result.get(4).get(i));
                    drawSegment(onechromBy, data.get(i), 10);
                }
            }

        });
        showButton.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
//                String filename = drawChromosome(chros.get(nums), 20);
//                new ImageFrame(filename);
                new Sandian("way3 50%", chros.get(nums).probes).draw();
            }
        });
        output.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                JFileChooser jfc = new JFileChooser(".");
                jfc.setFileSelectionMode(JFileChooser.FILES_AND_DIRECTORIES);
                //jfc.addChoosableFileFilter(new TxtFileFilter());
                jfc.showDialog(new JLabel(), "Select");
                File file = jfc.getSelectedFile();
                txtOutputfile.setText(file.getName());
                if (file.isDirectory()) {
                    System.out.println("Directory:" + file.getAbsolutePath());
                } else if (file.isFile()) {
                    System.out.println("File:" + file.getAbsolutePath());
                }
                System.out.println(jfc.getSelectedFile().getName());
            }
        });
        saveButton.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                try {
                    writeTotxt(chros);
                } catch (IOException e1) {
                    e1.printStackTrace();
                }
            }
        });

        tabbedPane1.addChangeListener(new ChangeListener() {
            @Override
            public void stateChanged(ChangeEvent e) {
                workModel = tabbedPane1.getSelectedIndex();
            }
        });
        clearButton.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                chros.clear();
                nums = -1;
            }
        });
        bypassRadioButton.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                transformMethod = 0;
            }
        });
        log2RadioButton.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                transformMethod = 1;
            }
        });
        power2RadioButton.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                transformMethod = 2;
            }
        });
        btnLoad.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                try {
                    boolean flag = chkIsMultiFiles.isSelected();
                    if (!flag) {
                        ArrayList<Chromosome> tmpData = readCNtxt(filepath);
//                        System.out.println("Loaded " + filepath);
                        chros.addAll(tmpData);
                    } else {
                        ArrayList<String> files = readFilePath(filepath);
                        for (String path : files) {
                            ArrayList<Chromosome> tmpData = readCNtxt(path);
                            chros.addAll(tmpData);
                        }
                    }
                } catch (IOException e1) {
                    e1.printStackTrace();
                }
                BtnCount.doClick();
            }
        });
        btnClear.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                chros.clear();
                BtnCount.doClick();
            }
        });
        BtnCount.addActionListener(new ActionListener() {
            @Override
            public void actionPerformed(ActionEvent e) {
                int count = chros.size();
                if (count > 1)
                    LblInfo.setText(String.format("[%03d] Samples are loaded.", count));
                else if (count == 1)
                    LblInfo.setText(String.format("[%03d] Samples is loaded.", 1));
                else
                    LblInfo.setText(String.format("None Sample is loaded."));
            }
        });
    }

    private void writeTotxt(List<Chromosome> sample) throws IOException {
        Date date = new Date();
        SimpleDateFormat df = new SimpleDateFormat("yyyyMMdd_HHmmss");
        String filename = String.format("data%1$sTestdata_%2$s.txt",
                File.separator, df.format(date));
        File output = new File(filename);
        output.createNewFile(); // Create a new file
        BufferedWriter out = new BufferedWriter(new FileWriter(output));
        out.write("chrId\t\t Loci\t\t CNValue\r\n");
        for (Chromosome chro : sample) {
            int index = 1;
            if (chro.probes.size() != 0) {
                for (double p : chro.probes) {
                    out.write(String.format("%2d\t\t %6d\t\t %.6f\n", chro.chrId, index, p));
                    index++;
                }
                StringBuffer sb = new StringBuffer();
                sb.append("changepoint: ");
                for (long lg : chro.changepoints)
                    sb.append(lg + "  ");
                out.write(sb.toString() + "\r\n");
            }
        }
        out.flush();
        out.close();
    }

    public String drawChromosome(Chromosome chro, int step) {
        int maxLength = 0;
        maxLength = chro.probes.size();
        int estimatedLength = maxLength;
        int width = (int) (1.15 * estimatedLength / step);
        int height = 400;
        BufferedImage image = new BufferedImage(width, height, BufferedImage.TYPE_INT_RGB);
        Graphics2D g = (Graphics2D) image.getGraphics();
        g.setColor(Color.WHITE);
        g.fillRect(0, 0, width, height);
        Stroke bs;
        bs = new BasicStroke(2.0f);
        g.setStroke(bs);
        int x, sx;
        Font f1 = new Font(null, Font.BOLD, 24);
        Font f2 = new Font(null, Font.BOLD, 12);
        x = 10;
        sx = 10 * step;
        g.setColor(Color.BLACK);
        g.setFont(f1);
        for (double cnValue : chro.probes) {
            int y = calculateYPosition(0, cnValue);
            g.setColor(Color.red);
            g.fillRect(x, y, 2, 2);
            sx++;
            x = sx / step;
        }
        x = 10;
        sx = 10 * step;
//        for (Segment seg:chro.segments){
//            if (seg.StartLoci==0){
//                sx = sx + seg.Loci.size();
//                x = sx / step;
//                continue;
//            }
//
//            if (chro.isTruePoint(30L,seg)){
//                g.setColor(Color.BLUE);
//                g.drawLine(x,100,x,400);
//            }
//            else{
//                g.setColor(Color.GREEN);
//                g.drawLine(x,100,x,400);
//            }
//            sx = sx + seg.Loci.size();
//            x = (int) ((sx / step) + 0.5);
//        }
        String filename = "";
        try {
            Date date = new Date();
            SimpleDateFormat df = new SimpleDateFormat("yyyyMMdd_HHmmss");
            filename = String.format(".%1$sResult%1$sChromosome_%2$s.png", File.separator, df.format(date));
            ImageIO.write(image, "PNG", new File(filename));
        } catch (Exception ex) {
            ex.printStackTrace();
        }
        return filename;
    }

    private void drawSegment(List<Set<Segment>> result, Chromosome chro, int step) {
        int maxLength = chro.getLength();
        int estimatedLength = maxLength;
        int width = (int) (1.5 * estimatedLength / step);
        int size = result.size();
        int height = 200 * (result.size() + 1);
        BufferedImage image = new BufferedImage(width, height, BufferedImage.TYPE_INT_RGB);
        Graphics2D g = (Graphics2D) image.getGraphics();
        g.setColor(Color.WHITE);
        g.fillRect(0, 0, width, height);
        Stroke bs;
        bs = new BasicStroke(2.0f);
        g.setStroke(bs);
        int x, sx;
        Font f1 = new Font(null, Font.BOLD, 24);
        Font f2 = new Font(null, Font.BOLD, 12);
        x = 100;
        sx = 100 * step;
        g.setColor(Color.BLACK);
        g.setFont(f1);
        double cnValue;
        int count = 0;
        for (int i = 0; i < result.size(); i++) {
            g.setColor(Color.BLACK);
            g.setFont(f1);
            switch (count % size) {
                case 0:
                    g.drawString(String.format("Chr_id = % 2d", i % size + 1), 10, 200 * (i % size + 1) - 150);
                    g.drawString("HMM", 20, 200 * (i % size + 1) - 130);
                    break;
                case 1:
                    g.drawString("FastPCF", 20, 200 * (i % size + 1) - 130);
                    break;
                case 2:
                    g.drawString("CLT", 20, 200 * (i % size + 1) - 130);
                    break;
                case 3:
                    g.drawString("Lasso", 20, 200 * (i % size + 1) - 100);
                    break;
                case 4:
                    g.drawString("CBS", 20, 200 * (i % size + 1) - 100);
            }
            x = 100;
            sx = 100 * step;
            for (double p : chro.probes) {
                cnValue = p;
                double lrr = BioToolbox.log2(cnValue);
                int y = calculateYPosition(i, lrr);
                g.setColor(Color.DARK_GRAY);
                g.fillRect(x, y, 2, 2);
                sx++;
                x = sx / step;
            }
            x = 100;
            sx = 100 * step;

            for (Segment s : result.get(i)) {
                cnValue = s.CopyNumber;
                int y2 = calculateYPosition(i, BioToolbox.log2(cnValue));
                g.setColor(Color.RED);
                int w = (int) (s.length() / step + 0.5);
                g.drawLine(x, y2, x + w, y2);
                g.setFont(f2);
                int textPos = (w > 30) ? x + w / 2 : x;
                g.drawString(String.format("[%.2f]", cnValue), textPos, y2 - 6);
                sx = sx + s.length();
                x = (int) ((sx / step) + 0.5);
            }
            count++;
        }
        try {
            Date date = new Date();
            SimpleDateFormat df = new SimpleDateFormat("yyyyMMdd_HHmmss");
            String filename = String.format(".%1$sResult%1$sSegmentChr%2$s%3$s.png", File.separator, chro.chrId, df.format(date));

            ImageIO.write(image, "PNG", new File(filename));
        } catch (Exception ex) {
            ex.printStackTrace();
        }
    }

    private int calculateYPosition(int row, double value) {
        // return (int) (value * 50) + 400*(row - 1);
        return (int) (150 - value * 40) + 200 * row;
    }

    public int[] totalTP(Chromosome chro, Set<Segment> segments, Long tolerance) {
        int[] result = new int[2];
        ArrayList<Long> points = new ArrayList<>(chro.changepoints);
        for (Segment seg : segments) {
            if (seg.Start() == 0)
                continue;
            int temp = seg.Start();
            Long down = Math.max(0, temp - tolerance);
            Long changepoint = 0L;
            Long up = temp + tolerance;
            boolean flag = false;
            for (Long lg : points) {
                if (lg >= down && lg <= up) {
                    changepoint = lg;
                    points.remove(lg);
                    flag = true;
                    break;
                }
            }
            if (flag == true) {
                result[1] += Math.abs(temp - changepoint);
                result[0]++;//number of true change point
            }
        }
        return result;
    }

    private ArrayList<String> readFilePath(String filepath) throws IOException {
        ArrayList<String> result = new ArrayList<>();
        BufferedReader in = new BufferedReader(new FileReader(filepath));
        String s = in.readLine();
        while (s != null && !s.isEmpty()) {
            result.add(s);
            s = in.readLine();
        }
        return result;
    }

    private ArrayList<Chromosome> readCNtxt(String filepath) throws IOException {
        BufferedReader in = new BufferedReader(new FileReader(filepath));
        in.readLine();
        String s = "";
        ArrayList<Chromosome> data = new ArrayList<>();
        for (int i = 0; i < 100; i++) {
            data.add(new Chromosome());
        }
        Chromosome chro = new Chromosome();
        int chrId = 1;
        chro.setChrId((short) chrId);
        ArrayList<Long> changePoints = new ArrayList<>();
        while ((s = in.readLine()) != null) {
            if (s.substring(0, 12).equals("changepoint:")) {
                s = s.substring(12);
                String[] arr = s.split("\\s+");
                List<String> arrayList = Arrays.asList(arr);
                //ArrayList<String> newarr=new ArrayList<>();
                for (String ss : arrayList) {
                    if (!(ss.matches("\\s+") || ss.isEmpty()))
                        changePoints.add(Long.valueOf(ss));
                }
            } else {
                String[] arr = s.split("\\s+");
                List<String> arrayList = Arrays.asList(arr);
                ArrayList<String> newarr = new ArrayList<>();
                for (String ss : arrayList) {
                    if (!(ss.matches("\\s+") || ss.isEmpty()))
                        newarr.add(ss);
                }
                if (chrId == Integer.valueOf(newarr.get(0))) {
//                    ProbeSet p=new ProbeSet();
//                    p.data.CopyNumberHalf=Double.valueOf(newarr.get(2));
//                    chro.probes.put(Long.valueOf(newarr.get(1)),p);
                    double cnhalf = Double.valueOf(newarr.get(2));
                    chro.probes.add(cnhalf);
                } else {
                    chrId = Integer.valueOf(newarr.get(0));
                    chro.changepoints = changePoints;
                    data.set(chro.chrId, chro);
                    chro = new Chromosome();
                    changePoints = new ArrayList<>();
                    chro.setChrId((short) chrId);
//                    ProbeSet p=new ProbeSet();
//                    p.data.CopyNumberHalf=Double.valueOf(newarr.get(2));
//                    chro.probes.put(Long.valueOf(newarr.get(1)), p);
                    double cnhalf = Double.valueOf(newarr.get(2));
                    chro.probes.add(cnhalf);
                }
            }
        }
        chro.changepoints = changePoints;
        data.set(chro.chrId, chro);
        ArrayList<Chromosome> result = new ArrayList<>();
        for (Chromosome chro1 : data) {
            if (chro1.probes.size() != 0)
                result.add(chro1);
        }
        return result;
    }

    private void readLenCn() {
        lengths.clear();
        cn.clear();
        variance.clear();
        String strlens = segmentlengthtxt.getText().replaceAll("\\s{1,}", "");
        String[] strlenss = strlens.split(",");
        String cns = cnvaluetxt.getText().replaceAll("\\s{1,}", "");
        String[] cnvalues = cns.split(",");
        String variances = variancetxt.getText().replaceAll("\\s{1,}", "");
        String[] varis = variances.split(",");
        for (String length : strlenss) {
            if (!length.isEmpty())
                lengths.add(Integer.valueOf(length));
        }

        for (String cnv : cnvalues) {
            if (!cnv.isEmpty())
                cn.add(Double.valueOf(cnv));
        }

        for (String var : varis) {
            if (!var.isEmpty())
                variance.add(Double.valueOf(Double.valueOf(var)));
        }

    }

}
