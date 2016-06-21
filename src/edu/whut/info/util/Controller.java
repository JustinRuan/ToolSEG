package edu.whut.info.util;

import edu.whut.info.core.CNSegment;
import edu.whut.info.dataset.ChromoGenerator;

/**
 * Created by Liu on 2016/3/18.
 */
public class Controller {
    public ChromoGenerator generator;
    public CNSegment CnSegment;

    public Controller() {
        generator = new ChromoGenerator();
        CnSegment = new CNSegment();
    }
}
