package uk.ac.ox.well.indiana.utils.io.table;

import com.google.common.base.Joiner;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class TableWriter {
    private PrintStream out;

    private List<String> header = null;

    public TableWriter(String fileToWrite) {
        initializeTable(new File(fileToWrite));
    }

    public TableWriter(File fileToWrite) {
        initializeTable(fileToWrite);
    }

    public TableWriter(PrintStream streamToWrite) {
        out = streamToWrite;
    }

    private void initializeTable(File fileToWrite) {
        try {
            out = new PrintStream(fileToWrite);
        } catch (IOException e) {
            throw new RuntimeException("Error while trying to write '" + fileToWrite.getAbsolutePath() + "'");
        }
    }

    public void addEntry(Map<String, String> entryMap) {
        if (header == null) {
            header = new ArrayList<String>(entryMap.keySet());

            out.println(Joiner.on("\t").join(header));
        }

        List<String> entryList = new ArrayList<String>();

        for (String field : header) {
            String entry = entryMap.get(field);
            entryList.add(entry);
        }

        out.println(Joiner.on("\t").join(entryList));
    }
}
