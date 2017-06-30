package uk.ac.ox.well.indiana.commands.playground.assemblies;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.reference.FastaSequenceFile;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import org.jetbrains.annotations.NotNull;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.containers.ContainerUtils;
import uk.ac.ox.well.indiana.utils.io.gff.GFF3;
import uk.ac.ox.well.indiana.utils.io.gff.GFF3Record;

import java.io.PrintStream;
import java.util.*;

/**
 * Created by kiran on 30/06/2017.
 */
public class DetermineAssemblyLayout extends Module {
    @Argument(fullName="reference", shortName="r", doc="Reference")
    public FastaSequenceFile REFERENCE;

    @Argument(fullName="draft", shortName="d", doc="Reference")
    public FastaSequenceFile DRAFT;

    @Argument(fullName="gff", shortName="g", doc="GFF3 file")
    public GFF3 GFF;

    @Argument(fullName="exons", shortName="e", doc="Aligned exons")
    public SamReader EXONS;

    @Output
    public PrintStream out;

    @Output(fullName="agpout", shortName="ao", doc="AGP out")
    public PrintStream aout;

    @Override
    public void execute() {
        Map<String, String> refSeqs = new HashMap<>();
        ReferenceSequence rseq;
        while ((rseq = DRAFT.nextSequence()) != null) {
            refSeqs.put(rseq.getName().split("\\s+")[0], rseq.getBaseString());
        }

        Map<String, String> contigSeqs = new HashMap<>();
        while ((rseq = DRAFT.nextSequence()) != null) {
            contigSeqs.put(rseq.getName().split("\\s+")[0], rseq.getBaseString());
        }

        Map<String, GFF3Record> grrecs = new HashMap<>();
        Map<String, SAMRecord> srrecs = new HashMap<>();

        for (GFF3Record gr : GFF) {
            if (gr.getType().equals("exon")) {
                grrecs.put(gr.getAttribute("ID"), gr);
            }
        }

        Map<String, Map<String, Integer>> freq = new HashMap<>();
        Map<String, Map<String, Integer>> strands = new HashMap<>();
        Map<String, Map<String, Set<GFF3Record>>> exonPlacements = new HashMap<>();

        for (SAMRecord sr : EXONS) {
            String id = sr.getReadName();
            srrecs.put(id, sr);

            if (grrecs.containsKey(id) && isHighQualityPlacement(sr)) {
                if (!freq.containsKey(sr.getReferenceName())) {
                    freq.put(sr.getReferenceName(), new HashMap<>());
                    strands.put(sr.getReferenceName(), new HashMap<>());
                    exonPlacements.put(sr.getReferenceName(), new HashMap<>());
                }

                if (!freq.get(sr.getReferenceName()).containsKey(grrecs.get(id).getSeqid())) {
                    freq.get(sr.getReferenceName()).put(grrecs.get(id).getSeqid(), 0);
                    strands.get(sr.getReferenceName()).put(grrecs.get(id).getSeqid(), 0);
                    exonPlacements.get(sr.getReferenceName()).put(grrecs.get(id).getSeqid(), new HashSet<>());
                }

                freq.get(sr.getReferenceName()).put(grrecs.get(id).getSeqid(), freq.get(sr.getReferenceName()).get(grrecs.get(id).getSeqid()) + sr.getReadLength());

                int mult = sr.getReadNegativeStrandFlag() ? -1 : 1;
                strands.get(sr.getReferenceName()).put(grrecs.get(id).getSeqid(), strands.get(sr.getReferenceName()).get(grrecs.get(id).getSeqid()) + mult*sr.getReadLength());

                ContainerUtils.add(exonPlacements.get(sr.getReferenceName()), grrecs.get(id).getSeqid(), grrecs.get(id));
            }
        }

        Map<String, Set<ContigInfo>> layout = new TreeMap<>();
        for (String contigName : freq.keySet()) {
            Map.Entry<String, Integer> chrNameAndCount = entriesSortedByValues(freq.get(contigName)).iterator().next();

            int strandBases = strands.get(contigName).get(chrNameAndCount.getKey());
            int lowestPos = Integer.MAX_VALUE;
            int highestPos = 0;

            for (GFF3Record gr : exonPlacements.get(contigName).get(chrNameAndCount.getKey())) {
                if (gr.getStart() < lowestPos) { lowestPos = gr.getStart(); }
                if (gr.getEnd() > highestPos) { highestPos = gr.getEnd(); }
            }

            if (!layout.containsKey(chrNameAndCount.getKey())) {
                layout.put(chrNameAndCount.getKey(), new TreeSet<>());
            }

            ContainerUtils.add(layout, chrNameAndCount.getKey(), new ContigInfo(contigName, contigSeqs.get(contigName), lowestPos, highestPos, strandBases >= 0));
        }

        //aout.println("##agp-version\t2.0");

        for (String chr : layout.keySet()) {
            log.info("chr: {} {}", chr, refSeqs.get(chr).length());

            for (ContigInfo ci : layout.get(chr)) {
                log.info("     {}", ci);
            }
        }
    }

    private class ContigInfo implements Comparable<ContigInfo> {
        public String contigName;
        public String contigSeq;
        public int start;
        public int end;
        public boolean isForward;

        public ContigInfo(String contigName, String contigSeq, int start, int end, boolean isForward) {
            this.contigName = contigName;
            this.contigSeq = contigSeq;
            this.start = start;
            this.end = end;
            this.isForward = isForward;
        }

        @Override
        public int compareTo(@NotNull ContigInfo o) {
            if (start == o.start) return 0;
            return start < o.start ? -1 : 1;
        }

        @Override
        public String toString() {
            return "ContigInfo{" +
                    "contigName='" + contigName + '\'' +
                    ", length=" + (contigSeq == null ? 0 : contigSeq.length()) +
                    ", start=" + start +
                    ", end=" + end +
                    ", isForward=" + isForward +
                    '}';
        }
    }

    private static boolean isHighQualityPlacement(SAMRecord sr) {
        int numClips = 0;

        for (CigarElement ce : sr.getCigar().getCigarElements()) {
            if (ce.getOperator().isClipping()) {
                numClips++;
            }
        }

        return !sr.getReadUnmappedFlag() && numClips == 0;
    }

    private static <K,V extends Comparable<? super V>> SortedSet<Map.Entry<K,V>> entriesSortedByValues(Map<K,V> map) {
        SortedSet<Map.Entry<K,V>> sortedEntries = new TreeSet<>(
                new Comparator<Map.Entry<K,V>>() {
                    @Override public int compare(Map.Entry<K,V> e1, Map.Entry<K,V> e2) {
                        int res = e2.getValue().compareTo(e1.getValue());
                        return res != 0 ? res : 1;
                    }
                }
        );

        sortedEntries.addAll(map.entrySet());
        return sortedEntries;
    }
}
