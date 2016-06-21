package edu.whut.info.util;

import javax.swing.*;
import java.io.File;
import java.text.SimpleDateFormat;
import java.util.Date;

/**
 * Created by Liu on 2016/1/5.
 */
public class ImageFrame extends JFrame {
    public ImageFrame(String filename) {
        super("second frame");
        Date date = new Date();
        SimpleDateFormat df = new SimpleDateFormat("yyyyMMdd_HHmmss");
        String fileName1 = String.format(".%1$sResult%1$sChromosome1_%2$s.png", File.separator, df.format(date));
        ImageTool.compressImage(filename, fileName1, 640, 640);
        Icon icon = new ImageIcon(fileName1);
        JLabel label = new JLabel();
        // setSize(imagePanel.getWidth(),imagePanel.getHeight());
        this.add(label);
        label.setIcon(icon);
        pack();
        setVisible(true);
    }
}
