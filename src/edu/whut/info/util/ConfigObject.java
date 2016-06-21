package edu.whut.info.util;


import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.NodeList;
import org.xml.sax.SAXException;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.ParserConfigurationException;
import javax.xml.transform.Transformer;
import javax.xml.transform.TransformerFactory;
import javax.xml.transform.dom.DOMSource;
import javax.xml.transform.stream.StreamResult;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;
import java.util.regex.Matcher;
import java.util.regex.Pattern;


/**
 * The type Config object.
 *
 * @author Justin
 * @version 1.0
 * @created 08 -11月-2014 10:18:54
 */
public class ConfigObject {

    /**
     * The Params.
     */
    public Parameters params;
    /**
     * The Sample files path
     */
    public Map<String, String[]> SampleFiles; //first is Normal, second is tumor
    public Set<String> TestingSampleIDs;
    private int SampleSize;
    private byte AnalysisMode;

    private String AnnotationFile;
    private String CDFFile;
    private String DBPath;

    private String samplePath;

    /**
     * Instantiates a new Config object.
     */
    public ConfigObject() {

        params = new Parameters();
        SampleFiles = new TreeMap<String, String[]>();
        TestingSampleIDs = new TreeSet<String>();
    }

    /**
     * Gets annotation file.path
     *
     * @return the annotation file
     */
    public String getAnnotationFile() {
        return AnnotationFile;
    }

    public void setAnnotationFile(String annotationFile) {
        AnnotationFile = annotationFile;
    }

    /**
     * Gets cDF file.path
     *
     * @return the cDF file
     */
    public String getCDFFile() {
        return CDFFile;
    }

    public void setCDFFile(String CDFFile) {
        this.CDFFile = CDFFile;
    }

    /**
     * Gets dB file.path
     *
     * @return the dB file
     */
    public String getDBPath() {
        return DBPath;
    }

    public void setDBPath(String DBPath) {
        this.DBPath = DBPath;
    }

    /**
     * Gets sample size.
     *
     * @return the sample size
     */
    public int getSampleSize() {
        return SampleSize;
    }

    /**
     * Gets analysis mode.
     *
     * @return the analysis mode
     */
    public byte getAnalysisMode() {
        return AnalysisMode;
    }

    /**
     * Gets sample path.
     *
     * @return the sample path
     */
    public String getSamplePath() {
        return samplePath;
    }

    public void setSamplePath(String samplePath) {
        this.samplePath = samplePath;
    }

    /**
     * Read config.
     *
     * @param filename the config XML filename and path
     */
    public void readConfig(String filename) {

        try {
            DocumentBuilderFactory dbf = DocumentBuilderFactory.newInstance();
            DocumentBuilder db = dbf.newDocumentBuilder();
            Document document = db.parse(filename);

            //Params
            NodeList nl = document.getElementsByTagName("ROI");
            int n = nl.getLength();
            for (int i = 0; i < n; i++) {
                String Chr_id = ((Element) nl.item(i)).getAttribute("Chr_id").trim();
                //Disable Strat，End Param
//                String Start = ((Element) nl.item(i)).getAttribute("Start").trim();
//                String End = ((Element) nl.item(i)).getAttribute("End").trim();
//
//                params.roiArray.add(new GeneRange(Short.parseShort(Chr_id), Long.parseLong(Start),
//                        Long.parseLong(End)));
                params.roiArray.add(new GeneRange(Short.parseShort(Chr_id), 0,
                        0));
            }
//            params.roi.Chr_id = Short.parseShort(nl.item(0).getTextContent());
//            nl = document.getElementsByTagName("Start");
//            params.roi.Start = Long.parseLong(nl.item(0).getTextContent());
//            nl = document.getElementsByTagName("End");
//            params.roi.End = Long.parseLong(nl.item(0).getTextContent());

            nl = document.getElementsByTagName("parameter");
            n = nl.getLength();
            for (int i = 0; i < n; i++) {
                String paramName = ((Element) nl.item(i)).getAttribute("name").trim();
                String paramValue = nl.item(i).getTextContent().trim();
                if (paramValue.isEmpty()) continue;

                if (paramName.equalsIgnoreCase("smoothingCopyNumber_Width")) {
                    params.smoothingCopyNumber_Width = Integer.parseInt(paramValue);
                } else if (paramName.equalsIgnoreCase("rou_WindowsSize")) {
                    params.rou_WindowsSize = Integer.parseInt(paramValue);
                } else if (paramName.equalsIgnoreCase("PVALUE_THRESH")) {
                    params.PVALUE_THRESH = Double.parseDouble(paramValue);
                } else if (paramName.equalsIgnoreCase("MIN_SEG_LENGTH")) {
                    params.MIN_SEG_LENGTH = Integer.parseInt(paramValue);
                } else if (paramName.equalsIgnoreCase("MIN_SEG_STEP")) {
                    params.MIN_SEG_STEP = Integer.parseInt(paramValue);
                } else if (paramName.equalsIgnoreCase("COMBINE_FACTOR")) {
                    params.COMBINE_FACTOR = Double.parseDouble(paramValue);
                } else if (paramName.equalsIgnoreCase("ALPHA_POSITION")) {
                    params.ALPHA_POSITION = Double.parseDouble(paramValue);
                } else if (paramName.equalsIgnoreCase("rou_chart")) {
                    params.enable_rou_chart = (Integer.parseInt(paramValue) > 0);
                } else if (paramName.equalsIgnoreCase("probesets_chart")) {
                    params.enable_probeset_chart = (Integer.parseInt(paramValue) > 0);
                }
            }

            //testing sample ids
            nl = document.getElementsByTagName("testingSample");
            n = nl.getLength();
            for (int i = 0; i < n; i++) {
                String temp = ((Element) nl.item(i)).getAttribute("id").trim();
                if (temp.isEmpty()) {
                    String startStr = ((Element) nl.item(i)).getAttribute("start_id").trim();
                    String endStr = ((Element) nl.item(i)).getAttribute("end_id").trim();
                    if (startStr.length() == endStr.length()) {
                        String regEx = "\\d+";
                        Pattern p = Pattern.compile(regEx);

                        Matcher mStart = p.matcher(startStr);
                        Matcher mEnd = p.matcher(endStr);
                        if (mStart.find() && mEnd.find()) {
                            String firstPart = startStr.substring(0, mStart.start());
                            int start = Integer.valueOf(mStart.group());
                            int end = Integer.valueOf(mEnd.group());

                            int len = startStr.length() - firstPart.length();
                            String tempP = String.format("%%s%%0%dd", len);
                            for (int j = start; j <= end; j++) {
                                TestingSampleIDs.add(String.format(tempP, firstPart, j));
                            }
                        }
                    }
                } else {
                    TestingSampleIDs.add(temp);
                }

            }

            //Files
            nl = document.getElementsByTagName("Annotation");
            AnnotationFile = nl.item(0).getTextContent().trim();
            nl = document.getElementsByTagName("CDF");
            CDFFile = nl.item(0).getTextContent().trim();

            nl = document.getElementsByTagName("DB");
            DBPath = nl.item(0).getTextContent().trim();
            if (!DBPath.endsWith(System.getProperty("file.separator"))) {
                DBPath = DBPath + System.getProperty("file.separator");
            }

            //samples
            nl = document.getElementsByTagName("sample_path");
            samplePath = nl.item(0).getTextContent().trim();

            if (!samplePath.endsWith(System.getProperty("file.separator"))) {
                samplePath = samplePath + System.getProperty("file.separator");
            }

            nl = document.getElementsByTagName("sample_pairs");
            SampleSize = nl.getLength();

            for (int i = 0; i < SampleSize; i++) {
                String t = nl.item(i).getTextContent();
                String[] ts = t.split(";");
                String groupid = ((Element) nl.item(i)).getAttribute("id").trim();

                if (ts.length == 2) {
                    for (int j = 0; j < 2; j++) {
                        //ts[j] = samplePath + ts[j].trim();
                        ts[j] = ts[j].trim();
                    }
                    SampleFiles.put(groupid, ts);
                } else {
                    throw new Exception(groupid + " are Unpaired samples.");
                }

            }

        } catch (ParserConfigurationException e) {
            e.printStackTrace();
        } catch (SAXException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        } catch (IOException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        } catch (Exception e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }
    }

    public void writeConfig(String filename) {
        try {
            DocumentBuilderFactory dbf = DocumentBuilderFactory.newInstance();
            DocumentBuilder db = dbf.newDocumentBuilder();
            Document document = db.parse(filename);

            // samples part
            NodeList nl = document.getElementsByTagName("samples");
            document.getDocumentElement().removeChild(nl.item(0));
            nl = document.getElementsByTagName("selected");
            document.getDocumentElement().removeChild(nl.item(0));


            Element selectedItems = document.createElement("selected");
            for (String id : TestingSampleIDs) {
                Element item = document.createElement("testingSample");
                item.setAttribute("id", id);
                selectedItems.appendChild(item);
            }
            document.getDocumentElement().appendChild(selectedItems);

            Element samples = document.createElement("samples");
            for (Map.Entry<String, String[]> record : SampleFiles.entrySet()) {
                Element sample = document.createElement("sample_pairs");
                sample.setAttribute("id", record.getKey());
                sample.setTextContent(String.format("%s;%s", record.getValue()[0], record.getValue()[1]));
                samples.appendChild(sample);
            }

            document.getDocumentElement().appendChild(samples);

            TransformerFactory transFactory = TransformerFactory.newInstance();
            Transformer transFormer = transFactory.newTransformer();
            DOMSource domSource = new DOMSource(document);
            //File file = new File("modify_"+filename);
            File file = new File(filename);
            if (file.exists()) {
                file.delete();
            }
            file.createNewFile();
            FileOutputStream out = new FileOutputStream(file);
            StreamResult xmlResult = new StreamResult(out);
            transFormer.transform(domSource, xmlResult);

        } catch (ParserConfigurationException e) {
            e.printStackTrace();
        } catch (SAXException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        } catch (IOException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        } catch (Exception e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }
    }

    public void addNewSamplePair(String id, String normalFile, String tumorFile) {
        if (id.length() > 0 && normalFile.length() > 0 && tumorFile.length() > 0)
            SampleFiles.put(id, new String[]{normalFile, tumorFile});
    }

    public void deleteSamplePair(String id) {
        SampleFiles.remove(id);
        TestingSampleIDs.remove(id);
    }

}