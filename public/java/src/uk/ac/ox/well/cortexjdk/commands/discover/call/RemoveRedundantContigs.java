package uk.ac.ox.well.cortexjdk.commands.discover.call;

import htsjdk.samtools.reference.FastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import uk.ac.ox.well.cortexjdk.commands.Module;
import uk.ac.ox.well.cortexjdk.utils.arguments.Argument;
import uk.ac.ox.well.cortexjdk.utils.arguments.Output;
import uk.ac.ox.well.cortexjdk.utils.io.graph.cortex.CortexGraph;
import uk.ac.ox.well.cortexjdk.utils.io.graph.cortex.CortexRecord;
import uk.ac.ox.well.cortexjdk.utils.kmer.CanonicalKmer;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

public class RemoveRedundantContigs extends Module {
    @Argument(fullName="contigs", shortName="c", doc="Contigs")
    public FastaSequenceFile CONTIGS;

    @Argument(fullName="roi", shortName="r", doc="ROIs")
    public CortexGraph ROI;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        Set<CanonicalKmer> rois = new HashSet<>();
        for (CortexRecord cr : ROI) {
            rois.add(cr.getCanonicalKmer());
        }

        List<ReferenceSequence> rseqs = new ArrayList<>();

        ReferenceSequence rseq;
        while ((rseq = CONTIGS.nextSequence()) != null) {
            rseqs.add(rseq);
        }

        Set<ReferenceSequence> toRemove = new HashSet<>();

        for (int i = 0; i < rseqs.size(); i++) {
            ReferenceSequence rseqi = rseqs.get(i);
            String seqi = rseqi.getBaseString();
            Set<CanonicalKmer> cksi = getUsedCanonicalKmers(seqi, rois);

            for (int j = i+1; j < rseqs.size(); j++) {
                ReferenceSequence rseqj = rseqs.get(j);
                String seqj = rseqj.getBaseString();
                Set<CanonicalKmer> cksj = getUsedCanonicalKmers(seqj, rois);

                double pctOverlap = pctOverlap(cksi, cksj);

                if (pctOverlap >= 0.1) {
                    log.info("{}={} {}={} {} {} {}", rseqi.getName().split(" ")[0], rseqi.length(), rseqj.getName().split(" ")[0], rseqj.length(), cksi.size(), cksj.size(), pctOverlap);

                    if (rseqi.length() < rseqj.length()) {
                        toRemove.add(rseqi);
                    } else {
                        toRemove.add(rseqj);
                    }
                }
            }
        }

        for (int i = 0; i < rseqs.size(); i++) {
            ReferenceSequence rseqi = rseqs.get(i);

            if (!toRemove.contains(rseqi)) {
                log.info("  retained: {}", rseqi.getName());
                out.println(">" + rseqi.getName());
                out.println(rseqi.getBaseString());
            } else {
                log.info(" * skipped: {}", rseqi.getName());
            }
        }
    }

    private Set<CanonicalKmer> getUsedCanonicalKmers(String seq, Set<CanonicalKmer> rois) {
        Set<CanonicalKmer> used = new HashSet<>();

        for (int i = 0; i <= seq.length() - ROI.getKmerSize(); i++) {
            CanonicalKmer ck = new CanonicalKmer(seq.substring(i, i + ROI.getKmerSize()));

            if (rois.contains(ck)) {
                used.add(ck);
            }
        }

        return used;
    }

    private double pctOverlap(Set<CanonicalKmer> x, Set<CanonicalKmer> y) {
        Set<CanonicalKmer> z = new HashSet<>();
        z.addAll(x);
        z.addAll(y);

        int overlap = 0;
        for (CanonicalKmer ck : z) {
            if (x.contains(ck) && y.contains(ck)) {
                overlap++;
            }
        }

        return 100.0 * (double) overlap / (double) z.size();
    }
}
