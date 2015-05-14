package uk.ac.ox.well.indiana.commands.gg;

import htsjdk.samtools.SAMFileReader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexKmer;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashSet;

public class ReadsContainingGivenKmers extends Module {
    @Argument(fullName="bam", shortName="b", doc="BAM file")
    public SAMFileReader BAM;

    @Argument(fullName="listOfKmers", shortName="l", doc="List of kmers")
    public HashSet<CortexKmer> KMER_LIST;

    @Output
    public File out;

    @Override
    public void execute() {
        int kmerSize = KMER_LIST.iterator().next().length();

        SAMFileWriter sfw = new SAMFileWriterFactory()
                .setCreateIndex(true)
                .makeBAMWriter(BAM.getFileHeader(), true, out);

        log.info("Looking for {} kmers in reads...", KMER_LIST.size());
        int numReads = 0, numReadsFound = 0;
        for (SAMRecord read : BAM) {
            if (numReads % 1000000 == 0) {
                log.info("  {} reads seen, {} reads with novel kmers", numReads, numReadsFound);
            }
            numReads++;

            String readSeq = read.getReadString();

            for (int i = 0; i <= readSeq.length() - kmerSize; i++) {
                CortexKmer ck = new CortexKmer(readSeq.substring(i, i + kmerSize));

                if (KMER_LIST.contains(ck)) {
                    sfw.addAlignment(read);

                    numReadsFound++;
                }
            }
        }

        log.info("  {} reads seen, {} reads with novel kmers", numReads, numReadsFound);

        sfw.close();
    }
}
