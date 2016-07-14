package uk.ac.ox.well.indiana.commands.playground;

import htsjdk.samtools.*;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.exceptions.IndianaException;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexKmer;
import uk.ac.ox.well.indiana.utils.progress.ProgressMeter;
import uk.ac.ox.well.indiana.utils.progress.ProgressMeterFactory;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

public class index extends Module {
    @Argument(fullName="bam", shortName="b", doc="BAM")
    public File SAMFILE;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        SamReader sreader = SamReaderFactory.make()
                .setOption(SamReaderFactory.Option.INCLUDE_SOURCE_IN_RECORDS, true)
                .setOption(SamReaderFactory.Option.CACHE_FILE_BASED_INDEXES, true)
                .open(SAMFILE);

        CortexKmer ck = new CortexKmer("GCACAAGCCTTTTCCATATACTTCCTATATTTTATATCAACATCACA");
        int kmerSize = 47;

        List<Chunk> chunks = new ArrayList<>();

        ProgressMeter pm = new ProgressMeterFactory()
                .header("Processing reads...")
                .updateRecord(1000000)
                .message("reads processed")
                .make(log);

        Set<SAMRecord> readsContainingKmer = new HashSet<>();

        for (SAMRecord sr : sreader) {
            BAMFileSpan sfs = (BAMFileSpan) sr.getFileSource().getFilePointer();

            for (int i = 0; i <= sr.getReadString().length() - kmerSize; i++) {
                CortexKmer cki = new CortexKmer(sr.getReadString().substring(i, i + kmerSize));

                if (ck.equals(cki)) {
                    chunks.addAll(sfs.getChunks());

                    readsContainingKmer.add(sr);
                }
            }

            pm.update("reads processed, " + readsContainingKmer.size() + " reads contain kmer");
        }

        log.info("Found {} reads", readsContainingKmer.size());

        chunks = Chunk.optimizeChunkList(chunks, 0);

        try {
            sreader.close();
        } catch (IOException e) {
            throw new IndianaException("Couldn't close file", e);
        }

        sreader = SamReaderFactory.makeDefault()
                .open(SAMFILE);

        int numReads = 0;
        for (Chunk chunk : chunks) {
            log.info("chunk: {}", chunk);
            SAMFileSpan sfs = new BAMFileSpan(chunk);

            SAMRecordIterator sri = sreader.indexing().iterator(sfs);

            while (sri.hasNext()) {
                SAMRecord sr = sri.next();

                numReads++;
            }

            log.info("Recovered {} reads from chunks", numReads);
        }

        log.info("");
    }
}
