package uk.ac.ox.well.indiana.analyses.kmerSharing;

import com.google.common.base.Joiner;
import net.sf.picard.reference.IndexedFastaSequenceFile;
import net.sf.picard.util.Interval;
import uk.ac.ox.well.indiana.tools.Tool;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.cortex.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.CortexRecord;
import uk.ac.ox.well.indiana.utils.io.gff.GFF3;
import uk.ac.ox.well.indiana.utils.io.gff.GFF3Record;
import uk.ac.ox.well.indiana.utils.sequence.SequenceUtils;

import java.io.PrintStream;
import java.util.HashMap;
import java.util.TreeSet;

public class CreateKmerReferencePanel extends Tool {
    @Argument(fullName="reference", shortName="R", doc="Reference genome")
    public IndexedFastaSequenceFile REFERENCE;

    @Argument(fullName="gff", shortName="gff", doc="GFF")
    public GFF3 GFF;

    @Argument(fullName="cortexGraph", shortName="cg", doc="Cortex graph")
    public CortexGraph CORTEX_GRAPH;

    @Argument(fullName="type", shortName="t", doc="Type of panel to create (UNIQUE or SHARED)")
    public String TYPE;

    @Output
    public PrintStream out;

    private class KmerInfo {
        public TreeSet<String> gene = new TreeSet<String>();
        public TreeSet<String> domains = new TreeSet<String>();
    }

    private HashMap<String, KmerInfo> getGeneKmersAndInfo() {
        HashMap<String, KmerInfo> kmers = new HashMap<String, KmerInfo>();

        for (GFF3Record r : GFF) {
            if ("gene".equalsIgnoreCase(r.getType())) {
                String seq = new String(REFERENCE.getSubsequenceAt(r.getSeqid(), r.getStart(), r.getEnd()).getBases());

                for (int i = 0; i < seq.length() - CORTEX_GRAPH.getKmerSize(); i++) {
                    String kmer = SequenceUtils.getAlphanumericallyLowestOrientation(seq.substring(i, i + CORTEX_GRAPH.getKmerSize()));

                    if (!kmers.containsKey(kmer)) {
                        kmers.put(kmer, new KmerInfo());
                    }
                    KmerInfo ki = kmers.get(kmer);

                    ki.gene.add(r.getAttributes().get("ID"));

                    Interval interval = new Interval(r.getSeqid(), r.getStart() + i, r.getStart() + i);

                    for (GFF3Record d : GFF.getType("domain", GFF.getOverlapping(interval))) {
                        ki.domains.add(d.getAttributes().get("DOMAIN_CLASS"));
                    }
                }
            }
        }

        return kmers;
    }

    private boolean isSharedKmer(CortexRecord cr) {
        int[] coverages = cr.getCoverages();

        boolean hasCoverage = false;
        int numColorsWithKmer = 0;
        boolean hasZeroOrUnitCoverageInColors = true;

        int totalCoverageInROI = 0;

        for (int color = 1; color < CORTEX_GRAPH.getNumColors(); color++) {
            hasCoverage |= (coverages[color] > 0);
            numColorsWithKmer += coverages[color] == 0 ? 0 : 1;
            hasZeroOrUnitCoverageInColors &= (coverages[color] <= 1);
            totalCoverageInROI += coverages[color];
        }

        boolean allCoverageIsInROI = (totalCoverageInROI == coverages[0]);

        return (hasCoverage && allCoverageIsInROI && numColorsWithKmer > 1 && hasZeroOrUnitCoverageInColors);
    }

    private boolean isUniqueKmer(CortexRecord cr) {
        int[] coverages = cr.getCoverages();

        boolean hasCoverage = false;
        int numColorsWithKmer = 0;
        int totalCoverageInROI = 0;
        boolean hasUnitCoverageInColor = false;

        for (int color = 1; color < CORTEX_GRAPH.getNumColors(); color++) {
            hasCoverage |= (coverages[color] > 0);
            numColorsWithKmer += coverages[color] == 0 ? 0 : 1;
            totalCoverageInROI += coverages[color];

            if (coverages[color] == 1) {
                hasUnitCoverageInColor = true;
            }
        }

        boolean allCoverageIsInROI = (totalCoverageInROI == coverages[0]);

        return (hasCoverage && allCoverageIsInROI && numColorsWithKmer == 1 && hasUnitCoverageInColor);
    }

    @Override
    public int execute() {
        if (!TYPE.equalsIgnoreCase("SHARED") && !TYPE.equalsIgnoreCase("UNIQUE")) {
            throw new RuntimeException("Argument for -type must be either SHARED or UNIQUE");
        }

        HashMap<String, KmerInfo> kmers = getGeneKmersAndInfo();

        out.println("kmer\tgenes\tdomains");

        int recordNum = 0;
        for (CortexRecord cr : CORTEX_GRAPH) {
            if (recordNum % (CORTEX_GRAPH.getNumRecords()/10) == 0) {
                log.info("processed {}/{} records", recordNum, CORTEX_GRAPH.getNumRecords());
            }
            recordNum++;

            if ( (TYPE.equalsIgnoreCase("SHARED") && isSharedKmer(cr) && kmers.containsKey(cr.getKmerString())) ||
                 (TYPE.equalsIgnoreCase("UNIQUE") && isUniqueKmer(cr) && kmers.containsKey(cr.getKmerString())) ) {
                String kmer = cr.getKmerString();
                KmerInfo ki = kmers.get(kmer);

                String genes = Joiner.on(",").join(ki.gene);
                String domains = ki.domains.isEmpty() ? "none" : Joiner.on(",").join(ki.domains);

                out.println(kmer + "\t" + genes + "\t" + domains);
            }
        }

        return 0;
    }
}
