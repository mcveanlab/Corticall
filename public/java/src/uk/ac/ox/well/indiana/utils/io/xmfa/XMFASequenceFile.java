package uk.ac.ox.well.indiana.utils.io.xmfa;

import com.google.common.base.Joiner;
import net.sf.picard.reference.ReferenceSequence;
import uk.ac.ox.well.indiana.utils.io.utils.LineReader;

import java.io.File;
import java.io.IOException;
import java.util.*;

public class XMFASequenceFile implements Iterable<XMFARecord>, Iterator<XMFARecord> {
    private File xmfaFile;
    private List<XMFARecord> records;

    private int index = 0;

    public XMFASequenceFile(String xmfaFilePath) {
        xmfaFile = new File(xmfaFilePath);
        loadXMFAFile(this.xmfaFile);
    }

    public XMFASequenceFile(File xmfaFile) {
        this.xmfaFile = xmfaFile;
        loadXMFAFile(this.xmfaFile);
    }

    private void loadXMFAFile(File xmfaFile) {
        LineReader lr = new LineReader(xmfaFile);

        records = new ArrayList<XMFARecord>();

        XMFARecord record = new XMFARecord();
        String fastaLine = null;
        List<String> contigPieces = new ArrayList<String>();

        String line;
        while ((line = lr.getNextRecord()) != null) {
            if (!line.startsWith("#")) {
                if (line.startsWith(">") || line.startsWith("=")) {
                    // Process previous contig
                    if (fastaLine != null) {
                        String seq = Joiner.on("").join(contigPieces);

                        String[] fastaLinePieces = fastaLine.split("\\s+");
                        String[] locus = fastaLinePieces[1].split("[:-]");
                        int contig = Integer.valueOf(locus[1]);
                        String name = new File(fastaLinePieces[3]).getName().replace(".fasta", "");

                        ReferenceSequence rseq = new ReferenceSequence(fastaLine, contig, seq.getBytes());

                        record.put(name, rseq);
                    }

                    if (line.startsWith(">")) {
                        // Set up for new contig
                        fastaLine = line;
                    } else {
                        // Sequence family is done; move on
                        records.add(record);

                        record = new XMFARecord();
                        fastaLine = null;
                    }

                    contigPieces.clear();
                } else {
                    // Add a line to the current contig
                    contigPieces.add(line);
                }
            }
        }

        records.add(record);
    }

    public int getNumRecords() { return records.size(); }

    @Override
    public Iterator<XMFARecord> iterator() {
        index = 0;
        return this;
    }

    @Override
    public boolean hasNext() {
        return index < records.size();
    }

    @Override
    public XMFARecord next() {
        if (this.hasNext()) {
            XMFARecord record = records.get(index);

            index++;

            return record;
        }

        return null;
    }

    @Override
    public void remove() {
        throw new UnsupportedOperationException();
    }
}
