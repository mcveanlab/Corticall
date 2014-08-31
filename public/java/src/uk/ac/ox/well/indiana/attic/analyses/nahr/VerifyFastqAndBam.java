package uk.ac.ox.well.indiana.attic.analyses.nahr;

import htsjdk.samtools.SAMFileReader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.fastq.FastqReader;
import htsjdk.samtools.fastq.FastqRecord;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.exceptions.IndianaException;

import java.io.PrintStream;

public class VerifyFastqAndBam extends Module {
    @Argument(fullName="end1", shortName="e1", doc="Fastq (end1)")
    public FastqReader END1;

    @Argument(fullName="end2", shortName="e2", doc="Fastq (end2)")
    public FastqReader END2;

    @Argument(fullName="bam", shortName="b", doc="BAM file")
    public SAMFileReader BAM;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        int numFr1 = 0, numFr2 = 0, numBr = 0;

        log.info("Reading fastq end1...");
        for (FastqRecord fr : END1) { numFr1++; }

        log.info("Reading fastq end2...");
        for (FastqRecord fr : END2) { numFr2++; }

        log.info("Reading BAM...");
        for (SAMRecord sr : BAM)  {
            if (!sr.isSecondaryOrSupplementary()) {
                numBr++;
            }
        }

        log.info("");
        log.info("Results:");
        log.info("  numFr1: {}", numFr1);
        log.info("  numFr2: {}", numFr2);
        log.info("   numFr: {}", numFr1 + numFr2);
        log.info("   numBr: {}", numBr);

        if (numFr1 + numFr2 != numBr) {
            throw new IndianaException("Number of records in fastqs and BAMs do not match (" + numFr1 + " + " + numFr2 + " != " + numBr + ")");
        }

        log.info("");
        log.info("Number of records in Fastq and BAM files match.");

        out.println("valid");
    }
}