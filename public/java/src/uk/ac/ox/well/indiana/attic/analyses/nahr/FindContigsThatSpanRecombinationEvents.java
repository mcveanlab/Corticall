package uk.ac.ox.well.indiana.attic.analyses.nahr;

import net.sf.picard.util.Interval;
import net.sf.samtools.CigarElement;
import net.sf.samtools.CigarOperator;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.table.TableReader;
import uk.ac.ox.well.indiana.utils.io.table.TableWriter;

import java.io.File;
import java.io.PrintStream;
import java.util.*;

public class FindContigsThatSpanRecombinationEvents extends Module {
    @Argument(fullName="bam", shortName="b", doc="BAM file")
    public ArrayList<SAMFileReader> BAMS;

    @Argument(fullName="recombLoci", shortName="rl", doc="Recomb loci table")
    public File RECOMB_LOCI;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        Map<String, Set<Interval>> loci = new HashMap<String, Set<Interval>>();

        TableReader tr = new TableReader(RECOMB_LOCI);

        for (Map<String, String> te : tr) {
            String[] sampleFields = te.get("sample").split("/");
            String sampleName = sampleFields[1];

            Interval locus = new Interval(te.get("chrom"), Integer.valueOf(te.get("pos_max")), Integer.valueOf(te.get("pos_max")));

            if (!loci.containsKey(sampleName)) {
                loci.put(sampleName, new HashSet<Interval>());
            }

            loci.get(sampleName).add(locus);
        }

        TableWriter tw = new TableWriter(out);

        for (SAMFileReader BAM : BAMS) {
            String sampleName = BAM.getFileHeader().getReadGroups().iterator().next().getSample();

            Map<String, Integer> numAlignments = new HashMap<String, Integer>();
            List<SAMRecord> contigs = new ArrayList<SAMRecord>();
            for (SAMRecord contig : BAM) {
                if (!numAlignments.containsKey(contig.getReadName())) {
                    numAlignments.put(contig.getReadName(), 1);
                } else {
                    numAlignments.put(contig.getReadName(), numAlignments.get(contig.getReadName()) + 1);
                }

                contigs.add(contig);
            }

            //for (SAMRecord contig : BAM) {
            for (SAMRecord contig : contigs) {
                Interval contigLocus = new Interval(contig.getReferenceName(), contig.getAlignmentStart(), contig.getAlignmentEnd());

                if (loci.containsKey(sampleName)) {
                    for (Interval knownLocus : loci.get(sampleName)) {
                        if (knownLocus.intersects(contigLocus)) {
                            boolean isClipped = false;
                            for (CigarElement ce : contig.getCigar().getCigarElements()) {
                                if (ce.getOperator().equals(CigarOperator.S) || ce.getOperator().equals(CigarOperator.H)) {
                                    isClipped = true;
                                    break;
                                }
                            }

                            boolean perfectAlignment = false;
                            if (contig.getCigar().getCigarElements().size() == 1 && contig.getCigar().getCigarElement(0).getOperator().equals(CigarOperator.M)) {
                                String md = contig.getStringAttribute("MD");
                                try {
                                    int length = Integer.parseInt(md);
                                    if (length == contig.getReadLength()) {
                                        perfectAlignment = true;
                                    }
                                } catch (NumberFormatException e) {}
                            }

                            Map<String, String> te = new LinkedHashMap<String, String>();
                            te.put("sampleName", sampleName);
                            te.put("contigName", contig.getReadName());
                            te.put("contigLocus", contigLocus.toString());
                            te.put("knownLocus", knownLocus.toString());
                            te.put("contigLength", String.valueOf(contig.getReadLength()));
                            te.put("numAlignments", String.valueOf(numAlignments.get(contig.getReadName())));
                            te.put("isClipped", isClipped ? "1" : "0");
                            te.put("perfectAlignment", perfectAlignment ? "1" : "0");
                            te.put("cigar", contig.getCigarString());
                            te.put("md", contig.getStringAttribute("MD"));

                            tw.addEntry(te);
                        }
                    }
                }
            }
        }
    }
}
