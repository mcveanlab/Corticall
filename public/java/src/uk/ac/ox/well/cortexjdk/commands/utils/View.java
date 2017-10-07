package uk.ac.ox.well.cortexjdk.commands.utils;

import com.google.common.base.Joiner;
import htsjdk.samtools.util.Interval;
import uk.ac.ox.well.cortexjdk.commands.Module;
import uk.ac.ox.well.cortexjdk.utils.alignment.kmer.KmerLookup;
import uk.ac.ox.well.cortexjdk.utils.arguments.Argument;
import uk.ac.ox.well.cortexjdk.utils.arguments.Output;
import uk.ac.ox.well.cortexjdk.utils.io.cortex.graph.CortexGraph;
import uk.ac.ox.well.cortexjdk.utils.io.cortex.graph.CortexKmer;
import uk.ac.ox.well.cortexjdk.utils.io.cortex.graph.CortexRecord;
import uk.ac.ox.well.cortexjdk.utils.sequence.SequenceUtils;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Set;
import java.util.TreeSet;

public class View extends Module {
    @Argument(fullName="graph", shortName="g", doc="Cortex graph")
    public CortexGraph GRAPH;

    @Argument(fullName="record", shortName="r", doc="Record to find and print", required=false)
    public ArrayList<String> RECORDS;

    @Argument(fullName="headerOnly", shortName="H", doc="Only print the file header", required=false)
    public Boolean HEADER_ONLY = false;

    @Argument(fullName="lookup", shortName="l", doc="Reference lookup", required=false)
    public KmerLookup LOOKUP;

    @Output
    public PrintStream out;

    private String recordToString(String sk, CortexRecord cr) {
        String kmer = cr.getKmerAsString();
        String cov = "";
        String ed = "";

        boolean fw = sk.equals(kmer);

        if (!fw) {
            kmer = SequenceUtils.reverseComplement(kmer);
        }

        for (int coverage : cr.getCoverages()) {
            cov += " " + coverage;
        }

        for (String edge : cr.getEdgeAsStrings()) {
            ed += " " + (fw ? edge : SequenceUtils.reverseComplement(edge));
        }

        Set<String> lss = new TreeSet<>();
        if (LOOKUP != null) {
            Set<Interval> loci = LOOKUP.find(kmer);

            if (loci != null && loci.size() > 0) {
                for (Interval locus : loci) {
                    String ls = locus.getContig() + ":" + locus.getStart() + "-" + locus.getEnd() + ":" + (locus.isPositiveStrand() ? "+" : "-");
                    lss.add(ls);
                }
            }
        }
        String lssCombined = Joiner.on(";").join(lss);

        return cr.getKmerAsString() + ": " + kmer + " " + cov + " " + ed + " " + lssCombined;
    }

    @Override
    public void execute() {
        if (HEADER_ONLY) {
            out.println(GRAPH);
        } else {
            if (RECORDS != null && !RECORDS.isEmpty()) {
                for (String seq : RECORDS) {
                    for (int i = 0; i <= seq.length() - GRAPH.getKmerSize(); i++) {
                        String kmer = seq.substring(i, i + GRAPH.getKmerSize());

                        CortexKmer ck = new CortexKmer(kmer);
                        CortexRecord cr = GRAPH.findRecord(ck);

                        if (cr != null) {
                            out.println(recordToString(kmer, cr));
                        } else {
                            out.println(kmer + ": missing");
                        }
                    }
                }
            } else {
                for (CortexRecord cr : GRAPH) {
                    out.println(cr);
                }
            }
        }
    }
}
