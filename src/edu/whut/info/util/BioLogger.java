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

package edu.whut.info.util;

import java.io.File;
import java.io.IOException;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.logging.*;

/**
 * Created by Justin on 2014/12/21.
 */
public class BioLogger {
    private Logger m_log;

    public BioLogger(String path, String filename) {
        m_log = Logger.getLogger("segment");
//        m_log = Logger.getGlobal();
        m_log.setLevel(Level.ALL);

        ConsoleHandler consoleHandler = new ConsoleHandler();
        consoleHandler.setLevel(Level.ALL);
        consoleHandler.setFormatter(new BioLogHander());
        m_log.addHandler(consoleHandler);

        try {
            String logfilename;
            if (!path.endsWith(File.separator)) {
                logfilename = String.format("%s%s%s", path, File.separator, filename);
            } else {
                logfilename = String.format("%s%s", path, filename);
            }
            //String logfilename = "/bai/ruanjun/Data/result.log";
            System.out.println(logfilename);
            FileHandler fileHandler = new FileHandler(logfilename, 500000, 10, true);

            fileHandler.setLevel(Level.ALL);

            fileHandler.setFormatter(new BioLogHander());
            m_log.addHandler(fileHandler);

            m_log.setUseParentHandlers(false);
        } catch (IOException e) {
            e.printStackTrace();
///////
            System.err.println("\n\n\n" + e.getClass().getName() + ": " + e.getMessage());
            System.exit(0);
        }

    }

    private class BioLogHander extends Formatter {
        private final DateFormat format = new SimpleDateFormat("hh:mm:ss:SSS");

        @Override
        public String format(LogRecord record) {
            return String.format("[%8d][%s] %s-> %s\n",
                    record.getSequenceNumber(), format.format(new Date()), record.getLevel(), record.getMessage());

        }
    }
}
