package uk.ac.ox.well.cortexjdk.utils.io.fq;

import htsjdk.samtools.fastq.FastqRecord;
import org.jetbrains.annotations.NotNull;

import java.io.*;
import java.util.Iterator;
import java.util.Scanner;
import java.util.zip.GZIPInputStream;

public class FqReader implements Iterator<FastqRecord>, Iterable<FastqRecord> {
    private InputStream is;
    private Scanner sc;

    public FqReader(File fileToRead) {
        try {
            is = new GZIPInputStream(new FileInputStream(fileToRead), 4096*64);

            moveToBeginningOfRecords();
        } catch (FileNotFoundException e) {
            throw new RuntimeException("Could not find file '" + fileToRead.getAbsolutePath() + "': " + e);
        } catch (IOException e) {
            throw new RuntimeException("Error while reading '" + fileToRead.getAbsolutePath() + "': " + e);
        }
    }

    private void moveToBeginningOfRecords() {
        sc = new Scanner(is);
    }

    @NotNull
    @Override
    public Iterator<FastqRecord> iterator() {
        moveToBeginningOfRecords();

        return this;
    }

    @Override
    public boolean hasNext() {
        return sc.hasNext();
    }

    private String nextLine() {
        return sc.nextLine();
    }

    @Override
    public FastqRecord next() {
        String seqHeaderPrefix = nextLine();
        String seqLine = nextLine();
        String qualHeaderPrefix = nextLine();
        String qualLine = nextLine();

        return new FastqRecord(seqHeaderPrefix, seqLine, qualHeaderPrefix, qualLine);
    }
}
