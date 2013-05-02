package uk.ac.ox.well.indiana.utils.io.utils;

import it.unimi.dsi.io.ByteBufferInputStream;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.nio.channels.FileChannel;
import java.util.*;

public class TableReader2 extends ArrayList<Map<String, String>> {
    public TableReader2(String fileToRead) {
        loadTable(new File(fileToRead));
    }

    public TableReader2(File fileToRead) {
        loadTable(fileToRead);
    }

    private void loadTable(File fileToRead) {
        try {
            FileInputStream fis = new FileInputStream(fileToRead);
            ByteBufferInputStream mappedRecordBuffer = ByteBufferInputStream.map(fis.getChannel(), FileChannel.MapMode.READ_ONLY);

            if (mappedRecordBuffer.length() < Integer.MAX_VALUE) {
                byte[] buffer = new byte[(int) mappedRecordBuffer.length()];

                mappedRecordBuffer.read(buffer);

                String[] header = null;

                int start = 0;
                for (int end = 0; end < buffer.length; end++) {
                    if (buffer[end] == '\n') {
                        byte[] subbuffer = new byte[end - start];
                        System.arraycopy(buffer, start, subbuffer, 0, end - start);

                        String[] fields = (new String(subbuffer)).split("\t");

                        if (start == 0) { // we're processing the header
                            header = fields;
                        } else {
                            Map<String, String> record = new HashMap<String, String>();

                            for (int i = 0; i < header.length; i++) {
                                record.put(header[i], fields[i]);
                            }

                            this.add(record);
                        }

                        start = end + 1;
                    }
                }
            } else {
                throw new RuntimeException("Handling for tables greater than " + Integer.MAX_VALUE + " is not currently supported");
            }
        } catch (FileNotFoundException e) {
            throw new RuntimeException("Could not find file '" + fileToRead.getAbsolutePath() + "': " + e);
        } catch (IOException e) {
            throw new RuntimeException("Error while reading '" + fileToRead.getAbsolutePath() + "': " + e);
        }
    }
}
