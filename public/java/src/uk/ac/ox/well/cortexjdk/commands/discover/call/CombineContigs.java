package uk.ac.ox.well.cortexjdk.commands.discover.call;

import htsjdk.samtools.reference.FastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.util.StringUtil;
import uk.ac.ox.well.cortexjdk.commands.Module;
import uk.ac.ox.well.cortexjdk.utils.arguments.Argument;
import uk.ac.ox.well.cortexjdk.utils.arguments.Output;
import uk.ac.ox.well.cortexjdk.utils.io.graph.cortex.CortexGraph;
import uk.ac.ox.well.cortexjdk.utils.io.graph.cortex.CortexRecord;
import uk.ac.ox.well.cortexjdk.utils.io.graph.links.CortexLinks;
import uk.ac.ox.well.cortexjdk.utils.kmer.CanonicalKmer;
import uk.ac.ox.well.cortexjdk.utils.sequence.SequenceUtils;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

public class CombineContigs extends Module {
    @Argument(fullName = "contigs", shortName = "c", doc = "Contigs")
    public FastaSequenceFile CONTIGS;

    @Argument(fullName = "partitions", shortName = "p", doc = "Partitions")
    public File PARTITIONS_FILE;

    @Argument(fullName = "roi", shortName = "r", doc = "ROI")
    public CortexGraph ROIS;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        FastaSequenceFile partitions = new FastaSequenceFile(PARTITIONS_FILE, true);
        List<ReferenceSequence> qseqs = new ArrayList<>();
        ReferenceSequence qseq;
        while ((qseq = partitions.nextSequence()) != null) {
            qseqs.add(qseq);
        }

        Set<CanonicalKmer> rois = loadRois(ROIS);

        ReferenceSequence rseq;
        while ((rseq = CONTIGS.nextSequence()) != null) {
            log.info("{}", rseq.getName());

            Set<CanonicalKmer> nks = getNovelKmerList(rois, rseq.getBaseString());

            ReferenceSequence bestSeq = null;
            int bestOverlap = 0;

            for (int i = 0; i < qseqs.size(); i++) {
                qseq = qseqs.get(i);

                Set<CanonicalKmer> qks = getNovelKmerList(rois, qseq.getBaseString());

                int overlap = computeOverlap(nks, qks);

                if (overlap > bestOverlap) {
                    bestOverlap = overlap;
                    bestSeq = qseq;
                }
            }

            if (bestOverlap > 0) {
                boolean isFwd = getOrientation(rseq.getBaseString(), bestSeq.getBaseString());

                String best = isFwd ? bestSeq.getBaseString() : SequenceUtils.reverseComplement(bestSeq.getBaseString());

                String newContig;

                if (rseq.getBaseString().contains(best)) {
                    newContig = rseq.getBaseString();
                } else if (best.contains(rseq.getBaseString())) {
                    newContig = best;
                } else {
                    int ib = 0, ir = 0;

                    for (int km = 1; km < 5 && ib >= 0 && ir >= 0 && km*ROIS.getKmerSize() < rseq.getBaseString().length() && km*ROIS.getKmerSize() < best.length(); km++) {
                        int k = km*ROIS.getKmerSize();
                        String rsk = rseq.getBaseString().substring(0, k);
                        String bsk = best.substring(0, k);

                        ib = rseq.getBaseString().indexOf(bsk);
                        ir = best.indexOf(rsk);
                    }

                    String mergedContig = "";
                    if (ir >= 0 && ib == -1) {
                        mergedContig = best.substring(0, ir) + rseq.getBaseString();
                    } else if (ib >= 0 && ir == -1) {
                        mergedContig = rseq.getBaseString().substring(0, ib) + best;
                    } else {
                        mergedContig = rseq.length() > best.length() ? rseq.getBaseString() : best;
                    }

                    log.debug("{} {}", ib, ir);
                    log.debug("{} {} {}", rseq.length(), best.length(), mergedContig.length());
                    log.debug("{}{}", ir == -1 ? "" : StringUtil.repeatCharNTimes(' ', ir), rseq.getBaseString());
                    log.debug("{}{}", ib == -1 ? "" : StringUtil.repeatCharNTimes(' ', ib), best);
                    log.debug("{}", mergedContig);
                    log.debug("");

                    if (mergedContig.length() > rseq.length() && mergedContig.length() > best.length()) {
                        newContig = mergedContig;
                    } else {
                        newContig = rseq.length() > best.length() ? rseq.getBaseString() : best;
                    }
                }

                out.println(">" + rseq.getName().split(" ")[0] + " len=" + (newContig.length() - ROIS.getKmerSize() + 1));
                out.println(newContig);
            }
        }
    }

    private boolean getOrientation(String rseq, String bestSeq) {
        Set<String> rsk = new HashSet<>();
        for (int i = 0; i <= rseq.length() - ROIS.getKmerSize(); i++) {
            rsk.add(rseq.substring(i, i + ROIS.getKmerSize()));
        }

        Set<String> qfwd = new HashSet<>();
        for (int i = 0; i <= bestSeq.length() - ROIS.getKmerSize(); i++) {
            qfwd.add(bestSeq.substring(i, i + ROIS.getKmerSize()));
        }

        Set<String> qrev = new HashSet<>();
        for (int i = 0; i <= bestSeq.length() - ROIS.getKmerSize(); i++) {
            qrev.add(SequenceUtils.reverseComplement(bestSeq.substring(i, i + ROIS.getKmerSize())));
        }

        return computeStringOverlap(rsk, qfwd) > computeStringOverlap(rsk, qrev);
    }

    private Set<CanonicalKmer> getNovelKmerList(Set<CanonicalKmer> rois, String seq) {
        Set<CanonicalKmer> nks = new HashSet<>();

        for (int i = 0; i <= seq.length() - ROIS.getKmerSize(); i++) {
            CanonicalKmer ck = new CanonicalKmer(seq.substring(i, i + ROIS.getKmerSize()));
            if (rois.contains(ck)) {
                nks.add(ck);
            }
        }

        return nks;
    }

    private Set<CanonicalKmer> loadRois(CortexGraph rois) {
        Set<CanonicalKmer> nks = new HashSet<>();
        for (CortexRecord cr : rois) {
            nks.add(cr.getCanonicalKmer());
        }

        return nks;
    }

    private int computeOverlap(Set<CanonicalKmer> nks, Set<CanonicalKmer> qks) {
        int overlap = 0;

        for (CanonicalKmer ck : qks) {
            if (nks.contains(ck)) {
                overlap++;
            }
        }

        return overlap;
    }

    private int computeStringOverlap(Set<String> nks, Set<String> qks) {
        int overlap = 0;

        for (String ck : qks) {
            if (nks.contains(ck)) {
                overlap++;
            }
        }

        return overlap;
    }
}
