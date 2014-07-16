package uk.ac.ox.well.indiana.attic.analyses.LongContigs;

import htsjdk.samtools.reference.FastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.table.TableWriter;
import uk.ac.ox.well.indiana.utils.sequence.SequenceUtils;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collection;
import java.util.LinkedHashMap;
import java.util.Map;

public class ContigStats2 extends Module {
    @Argument(fullName="contigs", shortName="c", doc="Contigs")
    public ArrayList<File> CONTIGS;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        TableWriter tw = new TableWriter(out);

        for (File f : CONTIGS) {
            log.info("id={}", f.getName());

            String[] id = f.getName().split("\\.");
            String idsample, idalg, idfam, idseed;

            if (f.getName().contains("supernode")) {
                idsample = id[0];
                idalg = "supernodes";
                idfam = "whole_genome";
                idseed = "none";
            } else {
                idsample = id[0];
                idalg = id[2];
                idfam = id[3];
                idseed = id[4];
            }

            Collection<String> contigs = new ArrayList<String>();

            FastaSequenceFile fasta = new FastaSequenceFile(f, true);
            ReferenceSequence rseq;
            while ((rseq = fasta.nextSequence()) != null) {
                String contig = new String(rseq.getBases());

                if (contig.length() > 0) {
                    contigs.add(contig);
                }
            }

            Map<String, String> te = new LinkedHashMap<String, String>();
            int minLength = SequenceUtils.minLength(contigs);
            int maxLength = SequenceUtils.maxLength(contigs);
            int n50Length = SequenceUtils.computeN50Length(contigs);
            int n50Value = SequenceUtils.computeN50Value(contigs);

            te.put("sample", idsample);
            te.put("algorithm", idalg);
            te.put("locus", idfam);
            te.put("seed", idseed);
            te.put("isuniqued", String.valueOf(f.getName().contains("unique")));
            te.put("numContigs", String.valueOf(contigs.size()));
            te.put("minLength", String.valueOf(minLength));
            te.put("maxLength", String.valueOf(maxLength));
            te.put("n50Length", String.valueOf(n50Length));
            te.put("n50Value", String.valueOf(n50Value));

            tw.addEntry(te);
        }
    }
}
