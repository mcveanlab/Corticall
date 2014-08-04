package uk.ac.ox.well.indiana.attic.analyses.nahr;

import htsjdk.samtools.reference.FastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;

import java.io.PrintStream;
import java.util.*;

public class IncorporateReaprContigs extends Module {
    @Argument(fullName="unbrokenContigs", shortName="uc", doc="Unbroken contigs")
    public FastaSequenceFile UNBROKEN_CONTIGS;

    @Argument(fullName="brokenContigs", shortName="bc", doc="Broken contigs")
    public FastaSequenceFile BROKEN_CONTIGS;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        Map<String, ReferenceSequence> finalContigs = new LinkedHashMap<String, ReferenceSequence>();

        ReferenceSequence rseq;
        while ((rseq = UNBROKEN_CONTIGS.nextSequence()) != null) {
            String name = rseq.getName().split(" ")[0];
            name = name.replaceAll("-", "_");

            finalContigs.put(name, rseq);
        }

        int numContigsBeforeRemoval = finalContigs.size();
        int numContigsAdded = 0;

        Set<String> contigsToRemove = new HashSet<String>();
        while ((rseq = BROKEN_CONTIGS.nextSequence()) != null) {
            String name = rseq.getName().split(" ")[0];
            name = name.replaceAll("-", "_");

            String seq = new String(rseq.getBases());
            String[] seqs = seq.split("N+");

            if (seqs.length > 1) {
                contigsToRemove.add(name);

                for (int i = 0; i < seqs.length; i++) {
                    String newName = name + "." + i;
                    ReferenceSequence rseqb = new ReferenceSequence(newName, i, seqs[i].getBytes());

                    finalContigs.put(newName, rseqb);

                    numContigsAdded++;
                }
            }
        }

        for (String contigToRemove : contigsToRemove) {
            finalContigs.remove(contigToRemove);
        }

        log.info("num contigs to remove: {}", contigsToRemove.size());
        log.info("num contigs to addition: {}", numContigsAdded);
        log.info("num contigs before processing: {}", numContigsBeforeRemoval);
        log.info("num contigs after processing: {}", finalContigs.size());

        for (String contigName : finalContigs.keySet()) {
            rseq = finalContigs.get(contigName);

            out.println(">" + contigName);
            out.println(new String(rseq.getBases()));
        }
    }
}
