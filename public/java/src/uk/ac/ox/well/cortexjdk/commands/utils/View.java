package uk.ac.ox.well.cortexjdk.commands.utils;

import com.google.common.base.Joiner;
import htsjdk.samtools.util.Interval;
import uk.ac.ox.well.cortexjdk.commands.Module;
import uk.ac.ox.well.cortexjdk.utils.alignment.reference.IndexedReference;
import uk.ac.ox.well.cortexjdk.utils.arguments.Argument;
import uk.ac.ox.well.cortexjdk.utils.arguments.Output;
import uk.ac.ox.well.cortexjdk.utils.io.graph.cortex.CortexGraph;
import uk.ac.ox.well.cortexjdk.utils.io.graph.links.CortexLinks;
import uk.ac.ox.well.cortexjdk.utils.io.graph.links.CortexLinksRecord;
import uk.ac.ox.well.cortexjdk.utils.kmer.CanonicalKmer;
import uk.ac.ox.well.cortexjdk.utils.io.graph.cortex.CortexRecord;
import uk.ac.ox.well.cortexjdk.utils.sequence.SequenceUtils;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;

public class View extends Module {
    @Argument(fullName="graph", shortName="g", doc="Cortex graph")
    public CortexGraph GRAPH;

    @Argument(fullName="record", shortName="r", doc="Record to find and print", required=false)
    public ArrayList<String> RECORDS;

    @Argument(fullName="links", shortName="l", doc="Links", required=false)
    public ArrayList<CortexLinks> LINKS;

    @Argument(fullName="headerOnly", shortName="H", doc="Only print the file header", required=false)
    public Boolean HEADER_ONLY = false;

    @Argument(fullName="lookup", shortName="l", doc="Reference lookup", required=false)
    public IndexedReference LOOKUP;

    @Output
    public PrintStream out;

    private String recordToString(String sk, CortexRecord cr, List<CortexLinksRecord> cl) {
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

        ed += " (l:" + cl.size() + ")";

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

                        CanonicalKmer ck = new CanonicalKmer(kmer);
                        CortexRecord cr = GRAPH.findRecord(ck);

                        List<CortexLinksRecord> cl = new ArrayList<>();
                        if (LINKS != null && LINKS.size() > 0) {
                            for (CortexLinks links : LINKS) {
                                if (links.containsKey(ck)) {
                                    cl.add(links.get(ck));
                                }
                            }
                        }

                        if (cr != null) {
                            out.println(recordToString(kmer, cr, cl));
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
