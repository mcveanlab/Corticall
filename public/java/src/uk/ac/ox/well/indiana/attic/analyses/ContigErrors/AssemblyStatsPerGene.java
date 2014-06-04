package uk.ac.ox.well.indiana.attic.analyses.ContigErrors;

import net.sf.picard.reference.FastaSequenceFile;
import net.sf.picard.reference.IndexedFastaSequenceFile;
import net.sf.picard.util.Interval;
import net.sf.picard.util.IntervalTree;
import net.sf.picard.util.IntervalTreeMap;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.gff.GFF3;
import uk.ac.ox.well.indiana.utils.io.gff.GFF3Record;
import uk.ac.ox.well.indiana.utils.io.table.TableReader;
import uk.ac.ox.well.indiana.utils.io.table.TableWriter;
import uk.ac.ox.well.indiana.utils.sequence.SequenceUtils;

import java.io.File;
import java.io.PrintStream;
import java.util.*;

public class AssemblyStatsPerGene extends Module {
    @Argument(fullName="contigs", shortName="c", doc="Aligned contigs (BAM)")
    public TreeMap<String, SAMFileReader> CONTIGS;

    @Argument(fullName="ref", shortName="r", doc="Reference FASTA")
    public IndexedFastaSequenceFile REF;

    @Argument(fullName="gff", shortName="gff", doc="GFF file")
    public GFF3 GFF;

    @Argument(fullName="geneClusters", shortName="clusters", doc="Clusters file")
    public File GENE_CLUSTERS;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        // Load gene cluster IDs
        log.info("Loading gene cluster ID information...");

        Map<String, String> geneIDToClusterID = new HashMap<String, String>();
        TableReader tr = new TableReader(GENE_CLUSTERS);
        for (Map<String, String> te : tr) {
            String groupName = te.get("groupName");
            String[] groupMembers = te.get("groupMembers").split(",");

            for (String groupMember : groupMembers) {
                geneIDToClusterID.put(groupMember, groupName.replace(" ", "_"));

                //log.info("{} -> {}", groupMember, groupName);
            }
        }

        // Load BAM data
        log.info("Loading aligned contigs...");

        Map<String, IntervalTreeMap<SAMRecord>> contigs = new HashMap<String, IntervalTreeMap<SAMRecord>>();
        for (String id : CONTIGS.keySet()) {
            log.info("  {}", id);

            contigs.put(id, new IntervalTreeMap<SAMRecord>());

            for (SAMRecord contig : CONTIGS.get(id)) {
                Interval interval = new Interval(contig.getReferenceName(), contig.getAlignmentStart(), contig.getAlignmentEnd());
                contigs.get(id).put(interval, contig);
            }
        }

        // Writing table
        TableWriter tw = new TableWriter(out);

        for (GFF3Record gr : GFF) {
            if (gr.getType().equals("gene")) {
                Map<String, String> te = new LinkedHashMap<String, String>();

                String seq = new String(REF.getSubsequenceAt(gr.getSeqid(), gr.getStart(), gr.getEnd()).getBases());

                String geneName = gr.getAttribute("ID");

                te.put("gene", geneName);
                te.put("clusterid", geneIDToClusterID.containsKey(geneName) ? geneIDToClusterID.get(geneName) : geneName);
                te.put("length", String.valueOf(seq.length()));
                te.put("gcFrac", String.valueOf(SequenceUtils.fractionGC(seq)));

                Interval interval = new Interval(gr.getSeqid(), gr.getStart(), gr.getEnd());

                for (String id : contigs.keySet()) {
                    Collection<SAMRecord> contigCollection = contigs.get(id).getOverlapping(interval);

                    Map<Integer, Boolean> isCovered = new HashMap<Integer, Boolean>();
                    for (int i = gr.getStart(); i <= gr.getEnd(); i++) {
                        isCovered.put(i, false);
                    }

                    int numPerfect = 0;
                    for (SAMRecord contig : contigCollection) {
                        String md = contig.getStringAttribute("MD");
                        if (md.equals(String.valueOf(contig.getReadLength()))) {
                            numPerfect++;
                        }

                        for (int i = contig.getAlignmentStart(); i <= contig.getAlignmentEnd(); i++) {
                            isCovered.put(i, true);
                        }
                    }

                    int totalLength = 0;
                    int coveredLength = 0;
                    for (Integer pos : isCovered.keySet()) {
                        if (isCovered.containsKey(pos) && isCovered.get(pos)) {
                            coveredLength++;
                        }
                        totalLength++;
                    }

                    te.put(id + ".numContigs", String.valueOf(contigCollection.size()));
                    te.put(id + ".numPerfectContigs", String.valueOf(numPerfect));
                    te.put(id + ".fractionPerfectContigs", String.valueOf((float) numPerfect / (float) contigCollection.size()));

                    te.put(id + ".totalLength", String.valueOf(totalLength));
                    te.put(id + ".coveredLength", String.valueOf(coveredLength));
                    te.put(id + ".fractionCovered", String.valueOf((float) coveredLength / (float) totalLength));
                }

                tw.addEntry(te);
            }
        }
    }
}
