package uk.ac.ox.well.cortexjdk.commands.utils;

import com.google.common.base.Joiner;
import htsjdk.samtools.reference.FastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import org.jgrapht.DirectedGraph;
import org.jgrapht.graph.DefaultDirectedGraph;
import org.jgrapht.graph.DefaultEdge;
import uk.ac.ox.well.cortexjdk.commands.Module;
import uk.ac.ox.well.cortexjdk.utils.arguments.Argument;
import uk.ac.ox.well.cortexjdk.utils.arguments.Output;
import uk.ac.ox.well.cortexjdk.utils.io.graph.cortex.CortexGraph;
import uk.ac.ox.well.cortexjdk.utils.io.graph.cortex.CortexRecord;
import uk.ac.ox.well.cortexjdk.utils.kmer.CanonicalKmer;
import uk.ac.ox.well.cortexjdk.utils.kmer.CortexByteKmer;
import uk.ac.ox.well.cortexjdk.utils.sequence.SequenceUtils;
import uk.ac.ox.well.cortexjdk.utils.traversal.TraversalUtils;

import java.io.PrintStream;
import java.util.*;

public class ToGfa1 extends Module {
    @Argument(fullName="graph", shortName="g", doc="Cortex graph")
    public CortexGraph GRAPH;

    @Argument(fullName="fasta", shortName="f", doc="Fasta (unitigs, or kmers, etc.)")
    public FastaSequenceFile FASTA;

    @Argument(fullName="sampleName", shortName="s", doc="Sample name", required=false)
    public String SAMPLE_NAME;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        Map<CanonicalKmer, CortexRecord> crs = new HashMap<>();
        for (CortexRecord cr : GRAPH) {
            crs.put(cr.getCanonicalKmer(), cr);
        }

        DirectedGraph<String, DefaultEdge> g = new DefaultDirectedGraph<>(DefaultEdge.class);
        Map<String, String> beginningKmers = new HashMap<>();
        Map<String, String> endingKmers = new HashMap<>();

        Map<String, Integer> seqNames = new HashMap<>();
        Map<String, Boolean> positiveStrand = new HashMap<>();

        ReferenceSequence rseq;
        int index = 0;
        while ((rseq = FASTA.nextSequence()) != null) {
            for (String v : Arrays.asList(rseq.getBaseString(), SequenceUtils.reverseComplement(rseq.getBaseString()))) {
                g.addVertex(v);

                String beginningSk = v.substring(0, GRAPH.getKmerSize());
                String endingSk = v.substring(v.length() - GRAPH.getKmerSize(), v.length());

                beginningKmers.put(beginningSk, v);
                endingKmers.put(endingSk, v);

                seqNames.put(v, index);
            }

            positiveStrand.put(rseq.getBaseString(), true);
            positiveStrand.put(SequenceUtils.reverseComplement(rseq.getBaseString()), false);

            index++;
        }

        int sampleColor = SAMPLE_NAME == null ? 0 : GRAPH.getColorForSampleName(SAMPLE_NAME);

        Map<String, Integer> averageCoverages = new HashMap<>();

        for (String v : g.vertexSet()) {
            int cov = 0;

            for (int i = 0; i <= v.length() - GRAPH.getKmerSize(); i++) {
                String sk = v.substring(i, i + GRAPH.getKmerSize());
                CanonicalKmer ck = new CanonicalKmer(sk);

                if (crs.containsKey(ck)) {
                    CortexRecord cr = crs.get(ck);
                    cov += cr.getCoverage(sampleColor);
                }
            }

            int avCov = (int) ((float) cov / (float) (v.length() - GRAPH.getKmerSize() + 1));

            averageCoverages.put(v, avCov);

            String beginningSk = v.substring(0, GRAPH.getKmerSize());
            CortexRecord beginningCr = crs.getOrDefault(new CanonicalKmer(beginningSk), null);

            Set<CortexByteKmer> ins = TraversalUtils.getAllPrevKmers(beginningCr, !beginningCr.getCanonicalKmer().getKmerAsString().equals(beginningSk)).get(sampleColor);
            for (CortexByteKmer cbk : ins) {
                String sk = new String(cbk.getKmer());

                if (endingKmers.containsKey(sk)) {
                    String v0 = endingKmers.get(sk);

                    g.addEdge(v0, v);
                }
            }

            String endingSk = v.substring(v.length() - GRAPH.getKmerSize(), v.length());
            CortexRecord endingCr = crs.getOrDefault(new CanonicalKmer(endingSk), null);

            Set<CortexByteKmer> outs = TraversalUtils.getAllNextKmers(endingCr, !endingCr.getCanonicalKmer().getKmerAsString().equals(endingSk)).get(sampleColor);
            for (CortexByteKmer cbk : outs) {
                String sk = new String(cbk.getKmer());

                if (beginningKmers.containsKey(sk)) {
                    String v1 = beginningKmers.get(sk);

                    g.addEdge(v, v1);
                }
            }
        }

        out.println(Joiner.on("\t").join("H", "VN:Z:1.0"));

        Set<Integer> seen = new HashSet<>();
        for (String v : g.vertexSet()) {
            int vid = seqNames.get(v);
            if (!seen.contains(vid)) {
                int avCov = averageCoverages.get(v);
                int totCov = avCov * v.length();

                out.println(Joiner.on("\t").join("S", vid, v, String.format("RC:i:%d", totCov), String.format("AC:i:%d", avCov)));

                seen.add(vid);
            }
        }

        for (DefaultEdge e : g.edgeSet()) {
            String vs = g.getEdgeSource(e);
            int vsid = seqNames.get(vs);
            String vsStrand = positiveStrand.get(vs) ? ":" : "-";

            String vt = g.getEdgeTarget(e);
            int vtid = seqNames.get(vt);
            String vtStrand = positiveStrand.get(vt) ? ":" : "-";

            out.println(Joiner.on("\t").join("L", vsid, vsStrand, vtid, vtStrand, GRAPH.getKmerSize() + "M"));
        }
    }
}
