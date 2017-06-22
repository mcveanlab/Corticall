package uk.ac.ox.well.indiana.commands.caller.call;

import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.reference.FastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.util.Interval;
import org.jetbrains.annotations.NotNull;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.alignment.kmer.KmerLookup;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.containers.ContainerUtils;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexKmer;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexRecord;

import java.io.PrintStream;
import java.util.*;

/**
 * Created by kiran on 22/06/2017.
 */
public class IdentifyNahrEvents extends Module {
    @Argument(fullName="graph", shortName="g", doc="Graph")
    public CortexGraph GRAPH;

    @Argument(fullName="child", shortName="c", doc="Child")
    public String CHILD;

    @Argument(fullName="parents", shortName="p", doc="Parents")
    public ArrayList<String> PARENTS;

    @Argument(fullName="roi", shortName="r", doc="ROI")
    public CortexGraph ROI;

    @Argument(fullName="refs", shortName="R", doc="References")
    public HashMap<String, KmerLookup> LOOKUPS;

    @Argument(fullName="sequences", shortName="s", doc="Contigs")
    public FastaSequenceFile CONTIGS;

    @Argument(fullName="validated", shortName="v", doc="Validated NAHR event")
    public FastaSequenceFile NAHR;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        Map<CortexKmer, Set<String>> expectedNovelKmers = new HashMap<>();
        ReferenceSequence rseq;
        while ((rseq = NAHR.nextSequence()) != null) {
            String seq = rseq.getBaseString();
            for (int i = 0; i <= seq.length() - GRAPH.getKmerSize(); i++) {
                CortexKmer ck = new CortexKmer(seq.substring(i, i + GRAPH.getKmerSize()));
                CortexRecord rr = ROI.findRecord(ck);

                if (rr != null) {
                    expectedNovelKmers.put(rr.getCortexKmer(), new TreeSet<>());
                }
            }
        }

        Map<String, String> contigs = new HashMap<>();

        while ((rseq = CONTIGS.nextSequence()) != null) {
            String[] name = rseq.getName().split("\\s+");
            String seq = rseq.getBaseString();
            contigs.put(name[0], seq);

            for (int i = 0; i <= seq.length() - GRAPH.getKmerSize(); i++) {
                CortexKmer ck = new CortexKmer(seq.substring(i, i + GRAPH.getKmerSize()));

                if (expectedNovelKmers.containsKey(ck)) {
                    ContainerUtils.add(expectedNovelKmers, ck, name[0]);
                }
            }
        }

        Set<String> cnames = new HashSet<>();
        for (CortexKmer ck : expectedNovelKmers.keySet()) {
            log.info("ck={} rseqs={}", ck, expectedNovelKmers.get(ck));
            cnames.addAll(expectedNovelKmers.get(ck));
        }

        Map<String, String> enc = createContigEncoding();

        for (String cname : cnames) {
            String contig = contigs.get(cname);

            log.info("length {}", contig.length());
            log.info("contig     {}", contig);

            for (String background : LOOKUPS.keySet()) {
                String anntig = annotateContig(LOOKUPS.get(background), contig, enc);
                log.info("anntig {} {}", background, anntig);
            }
        }
    }

    private String annotateContig(KmerLookup kl, String contig, Map<String, String> contigEncoding) {
        StringBuilder sb = new StringBuilder();

        for (int i = 0; i <= contig.length() - GRAPH.getKmerSize(); i++) {
            String sk = contig.substring(i, i + GRAPH.getKmerSize());
            CortexKmer ck = new CortexKmer(sk);
            CortexRecord rr = ROI.findRecord(ck);
            Set<Interval> loci = kl.findKmer(sk);

            if (rr != null) {
                sb.append(".");
            } else if (loci.size() == 1) {
                Interval locus = loci.iterator().next();
                sb.append(contigEncoding.get(locus.getContig()));
            } else {
                sb.append("_");
            }
        }

        return sb.toString();
    }

    @NotNull
    private Map<String, String> createContigEncoding() {
        Map<String, String> contigEncoding = new HashMap<>();
        String alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789";
        Random r = new Random();
        for (String key : LOOKUPS.keySet()) {
            Set<String> usedCodes = new HashSet<>();

            for (SAMSequenceRecord ssr : LOOKUPS.get(key).getReferenceSequence().getSequenceDictionary().getSequences()) {
                String sn = ssr.getSequenceName();
                String c;
                do {
                    c = String.valueOf(alphabet.charAt(r.nextInt(alphabet.length())));
                } while(usedCodes.contains(c) && usedCodes.size() < alphabet.length());

                if (!usedCodes.contains(c)) {
                    contigEncoding.put(sn, c);
                    usedCodes.add(c);
                }
            }
        }
        return contigEncoding;
    }
}
