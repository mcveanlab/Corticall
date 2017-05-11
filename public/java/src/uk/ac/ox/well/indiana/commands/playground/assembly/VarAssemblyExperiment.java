package uk.ac.ox.well.indiana.commands.playground.assembly;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.util.StringUtil;
import org.jgrapht.DirectedGraph;
import org.jgrapht.graph.DefaultDirectedGraph;
import org.jgrapht.traverse.DepthFirstIterator;
import org.jgrapht.traverse.TopologicalOrderIterator;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexKmer;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexRecord;
import uk.ac.ox.well.indiana.utils.io.gff.GFF3;
import uk.ac.ox.well.indiana.utils.io.gff.GFF3Record;
import uk.ac.ox.well.indiana.utils.sequence.CortexUtils;
import uk.ac.ox.well.indiana.utils.stoppingconditions.NahrStopper;
import uk.ac.ox.well.indiana.utils.traversal.AnnotatedEdge;
import uk.ac.ox.well.indiana.utils.traversal.AnnotatedVertex;
import uk.ac.ox.well.indiana.utils.traversal.TraversalEngine;
import uk.ac.ox.well.indiana.utils.traversal.TraversalEngineFactory;

import java.util.*;

/**
 * Created by kiran on 08/05/2017.
 */
public class VarAssemblyExperiment extends Module {
    @Argument(fullName="graph", shortName="g", doc="Graph")
    public CortexGraph GRAPH;

    @Argument(fullName="child", shortName="c", doc="Child sample")
    public String CHILD;

    @Argument(fullName="parent", shortName="p", doc="Parent sample")
    public ArrayList<String> PARENTS;

    @Argument(fullName="roi", shortName="i", doc="ROI")
    public CortexGraph ROI;

    @Argument(fullName="reference", shortName="r", doc="Reference")
    public IndexedFastaSequenceFile REF;

    @Argument(fullName="refgraph", shortName="rg", doc="Reference graph")
    public CortexGraph REFGRAPH;

    @Argument(fullName="geneModel", shortName="m", doc="GFF file")
    public GFF3 GFF;

    @Override
    public void execute() {
        int childColor = GRAPH.getColorForSampleName(CHILD);
        Set<Integer> parentColors = new HashSet<>(GRAPH.getColorsForSampleNames(PARENTS));

        log.info("Colors:");
        log.info(" -   child: {}", GRAPH.getColorForSampleName(CHILD));
        log.info(" - parents: {}", GRAPH.getColorsForSampleNames(PARENTS));

        /*
        TraversalEngine e = new TraversalEngineFactory()
                .graph(GRAPH)
                .traversalSamples(CHILD)
                .recruitmentSamples(PARENTS)
                .make();
                */
        //DirectedGraph<AnnotatedVertex, AnnotatedEdge> vg = new DefaultDirectedGraph<>(AnnotatedEdge.class);
        Map<String, Set<String>> vars = new HashMap<>();

        for (GFF3Record gr : GFF) {
            if (gr.getType().equals("gene") &&
                    gr.getAttribute("description").contains("VAR") &&
                    !gr.getAttribute("description").contains("pseudogene") &&
                    !gr.getAttribute("description").contains("truncated") &&
                    !gr.getAttribute("description").contains("putative") &&
                    !gr.getAttribute("description").contains("like")) {

                log.info("{}", gr);

                String seq = REF.getSubsequenceAt(gr.getSeqid(), gr.getStart(), gr.getEnd()).getBaseString();
                vars.put(gr.getAttribute("ID"), new HashSet<>());

                StringBuilder sb = new StringBuilder();
                StringBuilder cb = new StringBuilder();
                int missing = 0;

                for (int i = 0; i <= seq.length() - GRAPH.getKmerSize(); i++) {
                    String sk = seq.substring(i, i + GRAPH.getKmerSize());
                    CortexKmer ck = new CortexKmer(sk);
                    CortexRecord cr = GRAPH.findRecord(ck);

                    if (cr != null) {
                        if (cr.getCoverage(childColor) > 0) {
                            cb.append("1");
                        } else {
                            cb.append(".");
                            missing++;
                        }
                    } else {
                        cb.append(".");
                        missing++;
                    }

                    if (cr == null) {
                        cr = REFGRAPH.findRecord(ck);
                    }

                    sb.append(cr == null ? "." : 1);
                }

                log.info("  {}", seq);
                log.info("  {}", sb.toString());
                log.info("  {}", cb.toString());
                log.info("  missing: {}", missing);

                for (int i = 0; i <= seq.length() - GRAPH.getKmerSize() - 1; i++) {
                    String sk0 = seq.substring(i, i + GRAPH.getKmerSize());
                    String sk1 = seq.substring(i + 1, i + 1 + GRAPH.getKmerSize());
                }
            }
        }

        int i = 0;
        for (GFF3Record gr : GFF) {
            if (gr.getType().equals("gene") &&
                gr.getAttribute("description").contains("VAR") &&
                !gr.getAttribute("description").contains("pseudogene") &&
                !gr.getAttribute("description").contains("truncated") &&
                !gr.getAttribute("description").contains("putative") &&
                !gr.getAttribute("description").contains("like")) {
                String seq = REF.getSubsequenceAt(gr.getSeqid(), gr.getStart(), gr.getEnd()).getBaseString();

                log.info("gr: {} {} {}", i, gr, seq);

                String sl = seq.substring(0, GRAPH.getKmerSize());
                for (int j = 0; j <= seq.length() - GRAPH.getKmerSize(); j++) {
                    String sk = seq.substring(j, j + GRAPH.getKmerSize());
                    CortexKmer ck = new CortexKmer(sk);

                    CortexRecord cr = GRAPH.findRecord(ck);

                    if (cr == null) {
                        log.info("{}:{} {} {} {} {}", gr.getAttribute("ID"), j, sk, null, null, null);
                    } else {
                        log.info("{}:{} {} {} {} {}", gr.getAttribute("ID"), j, sk, cr.getCoverage(GRAPH.getColorForSampleName(CHILD)), cr.getEdgesAsString(GRAPH.getColorForSampleName(CHILD)), cr);
                    }

                    if (cr != null && cr.getCoverage(childColor) == 0) {
                        Set<String> nextKmers = CortexUtils.getNextKmers(GRAPH, sl, childColor);

                        for (String nextKmer : nextKmers) {
                            DirectedGraph<AnnotatedVertex, AnnotatedEdge> g = CortexUtils.dfs(GRAPH, null, nextKmer, childColor, parentColors, null, NahrStopper.class, 0, true, new HashSet<>());

                            if (g != null) {
                                DepthFirstIterator<AnnotatedVertex, AnnotatedEdge> dfi = new DepthFirstIterator<>(g, new AnnotatedVertex(sl));

                                StringBuilder sb = new StringBuilder();
                                StringBuilder ab = new StringBuilder(StringUtil.repeatCharNTimes(' ', GRAPH.getKmerSize() - 1));

                                while (dfi.hasNext()) {
                                    AnnotatedVertex av = dfi.next();
                                    CortexRecord ar = GRAPH.findRecord(new CortexKmer(av.getKmer()));
                                    log.info("{} {} {} {} {} {} {} {} {}",
                                            av.getKmer(),
                                            ar.getCortexKmer(),
                                            ar.getCoverage(8), ar.getCoverage(0), ar.getCoverage(17),
                                            ar.getEdgesAsString(8), ar.getEdgesAsString(0), ar.getEdgesAsString(17),
                                            ROI.findRecord(ar.getCortexKmer()) == null ? "" : "*"
                                    );

                                    if (sb.length() == 0) {
                                        sb.append(av.getKmer());
                                    } else {
                                        sb.append(av.getKmer().substring(av.getKmer().length() - 1));
                                    }

                                    //if ()
                                }


                                /*
                                for (String s : vars.values()) {
                                    log.info(" - {}", s);
                                }
                                log.info(" - {}", sb.toString());
                                */

                                log.info("");
                            }
                        }
                    }

                    sl = sk;
                }

                i++;
            }
        }

        /*
        ReferenceSequence rseq;
        while ((rseq = REF.nextSequence()) != null) {
            String seq = rseq.getBaseString();

            for (int i = 0; i <= seq.length() - GRAPH.getKmerSize(); i++) {
                String sk = seq.substring(i, i + GRAPH.getKmerSize());
                CortexKmer ck = new CortexKmer(sk);

                CortexRecord cr = GRAPH.findRecord(ck);

                if (cr == null) {
                    log.info("{}:{} {} {} {} {}",
                            rseq.getName(),
                            i,
                            sk,
                            null,
                            null,
                            null
                    );
                } else {
                    log.info("{}:{} {} {} {} {}",
                            rseq.getName(),
                            i,
                            sk,
                            cr.getCoverage(GRAPH.getColorForSampleName(CHILD)),
                            cr.getEdgesAsString(GRAPH.getColorForSampleName(CHILD)),
                            cr
                    );
                }
            }
        }
        */
    }
}
