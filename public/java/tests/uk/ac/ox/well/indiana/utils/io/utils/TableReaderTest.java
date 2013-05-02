package uk.ac.ox.well.indiana.utils.io.utils;

import org.testng.annotations.Test;

import java.util.Map;

public class TableReaderTest {
    @Test
    public void readFile() {
        TableReader2 tr = new TableReader2("testdata/smallReferencePanel.table");

        for (Map<String, String> te : tr) {
            System.out.println(te);
        }

        System.out.println("done");
    }
}
