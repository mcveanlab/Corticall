package uk.ac.ox.well.indiana.utils.io.utils;

import java.io.*;
import java.util.Iterator;

public class LineReader extends Reader {
    protected BufferedReader br = null;

    protected String nextRecord = null;

    public LineReader(File fileToRead) {
        try {
            br = new BufferedReader(new FileReader(fileToRead));

            nextRecord = br.readLine();
        } catch (FileNotFoundException e) {
            throw new RuntimeException(e);
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }

    @Override
    public int read(char[] chars, int i, int i1) throws IOException {
        return 0;
    }

    @Override
    public void close() throws IOException {
        br.close();
    }

    public boolean hasNext() {
        return nextRecord != null;
    }

    public String getNextRecord() {
        try {
            String currentRecord = nextRecord;

            nextRecord = br.readLine();
            if (nextRecord == null) {
                close();
            }

            return currentRecord;
        } catch (IOException e) {
            return null;
        }
    }
}
