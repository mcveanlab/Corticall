package uk.ac.ox.well.indiana.utils.io.table;

import ch.qos.logback.classic.Logger;
import it.unimi.dsi.io.ByteBufferInputStream;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.nio.channels.FileChannel;
import java.util.*;

/**
 * Streams very large text-based tables by first making a pass through the file to find all the line breaks
 * and once again to determine the positions of those breaks (and thus the record lengths).
 */
public class TableReader implements Iterable<Map<String, String>>, Iterator<Map<String, String>> {
    private ByteBufferInputStream mappedRecordBuffer;

    private List<Long> lineBreakPositions;
    private String[] header;
    private int nextRecordIndex;
    private int firstRecord = 1;

    private Logger log;

    public TableReader(String fileToRead) {
        loadTable(new File(fileToRead), null, null);
    }

    public TableReader(File fileToRead) {
        loadTable(fileToRead, null, null);
    }

    public TableReader(String fileToRead, String[] header) {
        loadTable(new File(fileToRead), header, null);
    }

    public TableReader(File fileToRead, String[] header) {
        loadTable(fileToRead, header, null);
    }

    public TableReader(String fileToRead, Logger log) {
        loadTable(new File(fileToRead), null, log);
    }

    public TableReader(File fileToRead, Logger log) {
        loadTable(fileToRead, null, log);
    }

    public TableReader(String fileToRead, String[] header, Logger log) {
        loadTable(new File(fileToRead), header, log);
    }

    public TableReader(File fileToRead, String[] header, Logger log) {
        loadTable(fileToRead, header, log);
    }

    private void loadTable(File fileToRead, String[] header, Logger log) {
        this.log = log;

        try {
            if (this.log != null) { this.log.info("Loading table records"); }

            FileInputStream fis = new FileInputStream(fileToRead);
            mappedRecordBuffer = ByteBufferInputStream.map(fis.getChannel(), FileChannel.MapMode.READ_ONLY);

            lineBreakPositions = new ArrayList<Long>();

            byte[] character = new byte[1];
            for (long position = 0; position < mappedRecordBuffer.length(); position++) {
                mappedRecordBuffer.read(character);

                if (character[0] == '\n') {
                    lineBreakPositions.add(position);
                }
            }

            mappedRecordBuffer.position(0);

            if (header == null) {
                byte[] headerBuffer = new byte[lineBreakPositions.get(0).intValue()];
                mappedRecordBuffer.read(headerBuffer);

                this.header = (new String(headerBuffer)).split("\t");

                firstRecord = 1;
            } else {
                this.header = header;
            }

            moveToBeginningOfRecords();
        } catch (FileNotFoundException e) {
            throw new RuntimeException("Could not find file '" + fileToRead.getAbsolutePath() + "': " + e);
        } catch (IOException e) {
            throw new RuntimeException("Error while reading '" + fileToRead.getAbsolutePath() + "': " + e);
        }
    }

    private void moveToBeginningOfRecords() {
        nextRecordIndex = firstRecord;

        if (nextRecordIndex > 0) {
            mappedRecordBuffer.position(lineBreakPositions.get(nextRecordIndex - 1).intValue() + 1);
        } else {
            mappedRecordBuffer.position(0);
        }
    }

    @Override
    public Iterator<Map<String, String>> iterator() {
        moveToBeginningOfRecords();

        return this;
    }

    @Override
    public boolean hasNext() {
        return nextRecordIndex < lineBreakPositions.size();
    }

    @Override
    public Map<String, String> next() {
        long start = lineBreakPositions.get(nextRecordIndex - 1) + 1;
        long end = lineBreakPositions.get(nextRecordIndex);
        int length = (int) (end - start);

        nextRecordIndex++;

        try {
            byte[] record = new byte[length];
            mappedRecordBuffer.read(record);
            mappedRecordBuffer.read(new byte[1]);

            String[] fields = (new String(record)).split("\t");
            //Map<String, String> entry = new HashMap<String, String>();
            Map<String, String> entry = new LinkedHashMap<String, String>();

            for (int i = 0; i < header.length; i++) {
                entry.put(header[i], fields[i]);
            }

            return entry;
        } catch (IOException e) {
            throw new RuntimeException("Unable to read record from table");
        }
    }

    @Override
    public void remove() {
        throw new UnsupportedOperationException();
    }

    public int size() {
        return lineBreakPositions.size();
    }
}
