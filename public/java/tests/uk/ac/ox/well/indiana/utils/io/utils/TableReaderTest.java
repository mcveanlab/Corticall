package uk.ac.ox.well.indiana.utils.io.utils;

import com.carrotsearch.sizeof.RamUsageEstimator;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.Test;
import uk.ac.ox.well.indiana.utils.io.table.TableReader;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;
import java.util.Random;

public class TableReaderTest {
    private File smallTable;
    private ArrayList<Map<String, String>> comparisonData;

    private String alphabet = "ACGT";
    private File largeTable;
    private Random randomGenerator;

    @BeforeClass
    public void loadSmallReferencePanelComparison() {
        smallTable = new File("testdata/smallReferencePanel.table");
        comparisonData = new ArrayList<Map<String, String>>();

        LineReader lr = new LineReader(smallTable);

        String[] header = null;

        while (lr.hasNext()) {
            String l = lr.getNextRecord();

            String[] fields = l.split("\t");

            if (header == null) {
                header = fields;
            } else {
                Map<String, String> te = new HashMap<String, String>();

                for (int i = 0; i < fields.length; i++) {
                    te.put(header[i], fields[i]);
                }

                comparisonData.add(te);
            }
        }
    }

    private void initializeRandomRecordGenerator() {
        randomGenerator = new Random(0);
    }

    private String generateRandomStringOfLengthN(int n) {
        StringBuilder sb = new StringBuilder();

        for (int i = 0; i < n; i++) {
            sb.append(alphabet.charAt(randomGenerator.nextInt(alphabet.length())));
        }

        return sb.toString();
    }

    private Map<String, String> generateRandomRecord() {
        Map<String, String> record = new HashMap<String, String>();

        record.put("kmer", generateRandomStringOfLengthN(30));
        record.put("gene", generateRandomStringOfLengthN(50));
        record.put("contig", generateRandomStringOfLengthN(70));

        return record;
    }

    @BeforeClass(enabled = false)
    public void writeLargeTable() {
        try {
            largeTable = new File("testdata/indianaLargeTableTest.table");

            if (!largeTable.exists()) {
                System.out.println("Generating data for large table reading test...");

                initializeRandomRecordGenerator();

                BufferedWriter bw = new BufferedWriter(new FileWriter(largeTable));

                bw.write("kmer\tgene\tcontig\n");

                int records = 0;
                while (largeTable.length() < Integer.MAX_VALUE) {
                    Map<String, String> te = generateRandomRecord();

                    bw.write(te.get("kmer") + "\t" + te.get("gene") + "\t" + te.get("contig") + "\n");

                    records++;
                }

                for (int i = 0; i < 5000; i++) {
                    Map<String, String> te = generateRandomRecord();

                    bw.write(te.get("kmer") + "\t" + te.get("gene") + "\t" + te.get("contig") + "\n");

                    records++;
                }

                bw.close();

                System.out.println("Wrote " + records + " records to " + largeTable.getAbsolutePath() + " (" + RamUsageEstimator.humanReadableUnits(largeTable.length()) + ")");
            } else {
                System.out.println("Using previously generated large table " + largeTable.getAbsolutePath());
            }
        } catch (IOException e) {
            throw new RuntimeException("Unable to create test table file for readLargeTable");
        }
    }

    @Test
    public void readSmallTable() {
        TableReader tr = new TableReader(smallTable);

        int index = 0;
        for (Map<String, String> te : tr) {
            Map<String, String> comp = comparisonData.get(index);

            Assert.assertEquals(comp, te);

            index++;
        }
    }

    @Test(enabled = false)
    public void readLargeTable() {
        TableReader tr = new TableReader(largeTable);

        LineReader lr = new LineReader(largeTable);
        String[] header = lr.getNextRecord().split("\t");

        int lineNum = 2;
        for (Map<String, String> te : tr) {
            String[] fields = lr.getNextRecord().split("\t");

            Map<String, String> comp = new HashMap<String, String>();
            for (int i = 0; i < header.length; i++) {
                comp.put(header[i], fields[i]);
            }

            Assert.assertEquals(comp, te, "lineNum = " + lineNum);

            lineNum++;
        }
    }
}
