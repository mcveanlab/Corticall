package uk.ac.ox.well.indiana.attic.analyses.ContigErrors;

import htsjdk.samtools.*;
import htsjdk.samtools.util.IntervalTree;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.table.TableWriter;

import java.io.PrintStream;
import java.util.*;
import java.util.regex.Pattern;

public class AlignedContigErrorStats extends Module {
    @Argument(fullName="bam", shortName="b", doc="Contigs (in BAM format)")
    public HashMap<String, SAMFileReader> BAMS;

    @Output
    public PrintStream out;

    private Pattern p = Pattern.compile("[0-9]+|[A-Z]|\\^[A-Z]+");

    private boolean isPerfectAlignment(SAMRecord contig) {
        String md = contig.getStringAttribute("MD");

        try {
            int length = Integer.parseInt(md);
            if (length == contig.getReadLength()) {
                return true;
            }
        } catch (NumberFormatException e) {}

        return false;
    }

    @Override
    public void execute() {
        TableWriter tw = new TableWriter(out);

        for (String id : BAMS.keySet()) {
            int alignedContigs = 0;
            int perfectContigs = 0;
            int unalignedContigs = 0;
            int totalGenomeSize = 0;

            SAMFileReader sfr = BAMS.get(id);

            Map<String, Integer> refNameSize = new HashMap<String, Integer>();
            for (SAMSequenceRecord sr : sfr.getFileHeader().getSequenceDictionary().getSequences()) {
                refNameSize.put(sr.getSequenceName(), sr.getSequenceLength());

                totalGenomeSize += sr.getSequenceLength();
            }

            Map<String, IntervalTree<String>> coveredIntervals = new TreeMap<String, IntervalTree<String>>();
            Map<String, Integer> numAlignedReads = new HashMap<String, Integer>();
            Map<String, Integer> numPerfectReads = new HashMap<String, Integer>();

            log.info("Processing {} reads...", id);
            for (SAMRecord contig : sfr) {
                if (contig.getAlignmentStart() > 0) {
                    alignedContigs++;

                    String refName = contig.getReferenceName();
                    int aStart = contig.getAlignmentStart();
                    int aEnd = contig.getAlignmentEnd();
                    if (aEnd == 0) {
                        aEnd = aStart;
                        for (CigarElement ce : contig.getCigar().getCigarElements()) {
                            CigarOperator co = ce.getOperator();

                            if (co.consumesReferenceBases()) {
                                aEnd += ce.getLength();
                            }
                        }
                    }

                    if (!numAlignedReads.containsKey(refName)) {
                        numAlignedReads.put(refName, 0);
                    }
                    numAlignedReads.put(refName, numAlignedReads.get(refName) + 1);

                    if (isPerfectAlignment(contig)) {
                        perfectContigs++;

                        if (!numPerfectReads.containsKey(refName)) {
                            numPerfectReads.put(refName, 0);
                        }

                        numPerfectReads.put(refName, numPerfectReads.get(refName) + 1);
                    }

                    if (!coveredIntervals.containsKey(refName)) {
                        coveredIntervals.put(refName, new IntervalTree<String>());
                    }

                    IntervalTree.Node<String> node = coveredIntervals.get(refName).minOverlapper(aStart, aEnd);
                    if (node == null) {
                        coveredIntervals.get(refName).put(aStart, aEnd, null);
                    } else {
                        int newStart = node.getStart() < aStart ? node.getStart() : aStart;
                        int newEnd = node.getEnd() > aEnd ? node.getEnd() : aEnd;

                        coveredIntervals.get(refName).remove(node.getStart(), node.getEnd());
                        coveredIntervals.get(refName).put(newStart, newEnd, null);
                    }
                } else {
                    unalignedContigs++;
                }
            }

            int totalCoveredLoci = 0;
            int totalPerfectContigs = 0;
            for (String refName : coveredIntervals.keySet()) {
                int coveredLoci = 0;
                for (IntervalTree.Node<String> node : coveredIntervals.get(refName)) {
                    coveredLoci += node.getLength();
                }
                totalCoveredLoci += coveredLoci;
                totalPerfectContigs += (numPerfectReads.containsKey(refName)) ? numPerfectReads.get(refName) : 0;

                float fractionCovered = (float) coveredLoci / (float) refNameSize.get(refName);
                float fractionPerfectContigs = (numPerfectReads.containsKey(refName)) ? (float) numPerfectReads.get(refName) / (float) numAlignedReads.get(refName) : 0;

                Map<String, String> entry = new LinkedHashMap<String, String>();
                entry.put("id", id);
                entry.put("refName", refName);
                entry.put("totalLength", String.valueOf(refNameSize.get(refName)));
                entry.put("coveredLength", String.valueOf(coveredLoci));
                entry.put("fractionCovered", String.valueOf(fractionCovered));
                entry.put("contigsUnaligned", "0");
                entry.put("contigsAligned", String.valueOf(numAlignedReads.get(refName)));
                entry.put("perfectContigs", String.valueOf(numPerfectReads.get(refName)));
                entry.put("fractionPerfectContigs", String.valueOf(fractionPerfectContigs));
                tw.addEntry(entry);
            }

            float totalFractionCovered = (float) totalCoveredLoci / (float) totalGenomeSize;
            float totalFractionPerfectContigs = (float) totalPerfectContigs / (float) alignedContigs;

            Map<String, String> entry = new LinkedHashMap<String, String>();
            entry.put("id", id);
            entry.put("refName", "total");
            entry.put("totalLength", String.valueOf(totalGenomeSize));
            entry.put("coveredLength", String.valueOf(totalCoveredLoci));
            entry.put("fractionCovered", String.valueOf(totalFractionCovered));
            entry.put("contigsUnaligned", String.valueOf(unalignedContigs));
            entry.put("contigsAligned", String.valueOf(alignedContigs));
            entry.put("perfectContigs", String.valueOf(perfectContigs));
            entry.put("fractionPerfectContigs", String.valueOf(totalFractionPerfectContigs));
            tw.addEntry(entry);
        }
    }
}
