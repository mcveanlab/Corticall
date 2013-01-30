package uk.ac.ox.well.indiana.utils.io.utils;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.Iterator;

public class TableReader extends LineReader implements Iterator<HashMap<String, String>>, Iterable<HashMap<String, String>> {
    private String[] header;

    public TableReader(File fileToRead) {
        super(fileToRead);

        header = this.nextRecord.split("\t");
        try {
            this.nextRecord = br.readLine();
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

    public HashMap<String, String> next() {
        String line = this.getNextRecord();
        String[] fields = line.split("\t");

        HashMap<String, String> entry = new HashMap<String, String>();

        if (fields.length == header.length) {
            for (int i = 0; i < fields.length; i++) {
                entry.put(header[i], fields[i]);
            }
        } else if (fields.length == header.length + 1) {
            entry.put("", fields[0]);
            for (int i = 1; i < fields.length; i++) {
                entry.put(header[i-1], fields[i]);
            }
        } else {
            throw new RuntimeException("Unable to load table: header has " + header.length + " fields but data has " + fields.length + ".");
        }

        return entry;
    }

    public void remove() {
        throw new UnsupportedOperationException();
    }

    public Iterator<HashMap<String, String>> iterator() {
        return this;
    }
}
