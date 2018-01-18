package uk.ac.ox.well.cortexjdk.commands.assemble.correct;

import com.google.common.collect.Sets;
import htsjdk.samtools.fastq.FastqRecord;
import org.apache.commons.math3.util.Pair;
import org.jgrapht.Graph;
import org.jgrapht.GraphPath;
import org.jgrapht.graph.DirectedWeightedPseudograph;
import uk.ac.ox.well.cortexjdk.commands.Module;
import uk.ac.ox.well.cortexjdk.utils.arguments.Argument;
import uk.ac.ox.well.cortexjdk.utils.arguments.Output;
import uk.ac.ox.well.cortexjdk.utils.io.graph.cortex.*;
import uk.ac.ox.well.cortexjdk.utils.io.reads.Reads;
import uk.ac.ox.well.cortexjdk.utils.kmer.CanonicalKmer;
import uk.ac.ox.well.cortexjdk.utils.kmer.CortexByteKmer;
import uk.ac.ox.well.cortexjdk.utils.progress.ProgressMeter;
import uk.ac.ox.well.cortexjdk.utils.progress.ProgressMeterFactory;
import uk.ac.ox.well.cortexjdk.utils.stoppingrules.GapClosingStopper;
import uk.ac.ox.well.cortexjdk.utils.traversal.*;

import java.io.File;
import java.util.*;

import static uk.ac.ox.well.cortexjdk.utils.traversal.TraversalEngineConfiguration.GraphCombinationOperator.OR;
import static uk.ac.ox.well.cortexjdk.utils.traversal.TraversalEngineConfiguration.TraversalDirection.BOTH;

public class CorrectGraph extends Module {
    @Argument(fullName="graph", shortName="g", doc="Graph")
    public CortexMap GRAPH;

    @Argument(fullName="clean", shortName="c", doc="Clean graph")
    public CortexGraph CLEAN;

    /*
    @Argument(fullName="dirty", shortName="d", doc="Dirty graph")
    public CortexGraph DIRTY;

    @Argument(fullName="ref", shortName="R", doc="Reference graphs", required = false)
    public ArrayList<CortexGraph> REFS;
    */

    @Argument(fullName="reads", shortName="r", doc="Reads")
    public ArrayList<File> READS;

    @Output
    public File out;

    @Override
    public void execute() {
        /*
        List<CortexGraph> graphs = new ArrayList<>();
        graphs.add(DIRTY);
        if (REFS != null) { graphs.addAll(REFS); }

        CortexCollection cc = new CortexCollection(graphs);
        */

        List<Integer> refColors = new ArrayList<>();
        //for (int c = 1; c < cc.getNumColors(); c++) {
        for (int c = 2; c < GRAPH.getNumColors(); c++) {
            refColors.add(c);
        }

        TraversalEngine e = new TraversalEngineFactory()
                //.traversalColor(0)
                .traversalColor(1)
                .recruitmentColors(refColors)
                .traversalDirection(BOTH)
                .combinationOperator(OR)
                .stoppingRule(GapClosingStopper.class)
                .graph(GRAPH)
                .make();

        Map<CanonicalKmer, CortexRecord> corrections = new TreeMap<>();

        for (File readsFile : READS) {
            Reads reads = new Reads(readsFile);

            ProgressMeter pm = new ProgressMeterFactory()
                    //.updateRecord(100000)
                    .updateRecord(GRAPH.getNumRecords() / 10)
                    .header("Processing '" + readsFile.getAbsolutePath() + "'")
                    .message("reads processed")
                    .make(log);

            int numGaps = 0, numGapsClosed = 0;

            for (Pair<FastqRecord, FastqRecord> p : reads) {
                for (FastqRecord fr : Arrays.asList(p.getFirst(), p.getSecond())) {
                    if (fr != null) {
                        String seq = fr.getReadString();
                        //CortexVertex[] cvs = new CortexVertex[seq.length() - CLEAN.getKmerSize() + 1];
                        CortexVertex[] cvs = new CortexVertex[seq.length() - GRAPH.getKmerSize() + 1];

                        //for (int i = 0; i <= seq.length() - CLEAN.getKmerSize(); i++) {
                        for (int i = 0; i <= seq.length() - GRAPH.getKmerSize(); i++) {
                            //String sk = seq.substring(i, i + CLEAN.getKmerSize());
                            String sk = seq.substring(i, i + GRAPH.getKmerSize());

                            if (!sk.contains("N")) {
                                CanonicalKmer ck = new CanonicalKmer(sk);

                                //CortexRecord cr = CLEAN.findRecord(ck);
                                CortexRecord cr = GRAPH.findRecord(ck);
                                //if (cr != null) {
                                if (cr != null && cr.getCoverage(0) > 0) {
                                    //CortexRecord dr = cc.findRecord(ck);
                                    //CortexVertex cv = new CortexVertex(new CortexByteKmer(sk), dr);
                                    CortexVertex cv = new CortexVertex(new CortexByteKmer(sk), cr);

                                    cvs[i] = cv;
                                }
                            }
                        }

                        Set<Pair<CortexVertex, CortexVertex>> lk = new LinkedHashSet<>();
                        Set<Pair<Integer, Integer>> li = new LinkedHashSet<>();

                        for (int i = 0; i < cvs.length; i++) {
                            if (cvs[i] == null) {
                                CortexVertex a = null;
                                int ai = -1;
                                for (int j = i - 1; j >= 0; j--) {
                                    if (cvs[j] != null) {
                                        a = cvs[j];
                                        ai = j;
                                        break;
                                    }
                                }

                                CortexVertex b = null;
                                int bi = -1;
                                for (int j = i + 1; j < cvs.length; j++) {
                                    if (cvs[j] != null) {
                                        b = cvs[j];
                                        bi = j;
                                        break;
                                    }
                                }

                                if (a != null && b != null) {
                                    lk.add(new Pair<>(a, b));
                                    li.add(new Pair<>(ai, bi));
                                }
                            }
                        }

                        for (Pair<CortexVertex, CortexVertex> l : lk) {
                            numGaps++;

                            CortexVertex cv0 = l.getFirst();
                            CortexVertex cv1 = l.getSecond();

                            DirectedWeightedPseudograph<CortexVertex, CortexEdge> pt = new DirectedWeightedPseudograph<>(CortexEdge.class);
                            pt.addVertex(cv1);

                            e.getConfiguration().setPreviousTraversal(pt);

                            DirectedWeightedPseudograph<CortexVertex, CortexEdge> g = e.dfs(l.getFirst().getKmerAsString());

                            if (g != null && !cv0.equals(cv1)) {
                                //PathFinder pf = new PathFinder(g, 0);
                                PathFinder pf = new PathFinder(g, 1);
                                GraphPath<CortexVertex, CortexEdge> gp = pf.getPathFinder(cv0, cv1);

                                if (gp != null) {
                                    List<CortexVertex> lv = gp.getVertexList();

                                    for (int i = 0; i < lv.size(); i++) {
                                        CortexVertex v0 = i == 0 ? null : lv.get(i - 1);
                                        CortexVertex v1 = lv.get(i);
                                        CortexVertex v2 = i == lv.size() - 1 ? null : lv.get(i + 1);

                                        Set<String> inEdges = new HashSet<>();
                                        if (v0 != null) {
                                            String s0 = v0.getKmerAsString();
                                            String inEdge = s0.substring(0, 1);
                                            inEdges.add(inEdge);
                                        }

                                        Set<String> outEdges = new HashSet<>();
                                        if (v2 != null) {
                                            String s2 = v2.getKmerAsString();
                                            String outEdge = s2.substring(s2.length() - 1, s2.length());
                                            outEdges.add(outEdge);
                                        }

                                        if (corrections.containsKey(v1.getCanonicalKmer())) {
                                            CortexRecord ocr = corrections.get(v1.getCanonicalKmer());

                                            if (v1.getKmerAsString().equals(ocr.getCanonicalKmer().getKmerAsString())) {
                                                inEdges.addAll(ocr.getInEdgesAsStrings(0, false));
                                                outEdges.addAll(ocr.getOutEdgesAsStrings(0, false));
                                            } else {
                                                inEdges.addAll(ocr.getOutEdgesAsStrings(0, true));
                                                outEdges.addAll(ocr.getInEdgesAsStrings(0, true));
                                            }
                                        }

                                        CortexRecord cr = new CortexRecord(
                                                v1.getKmerAsString(),
                                                Arrays.asList(v1.getCortexRecord().getCoverage(0)),
                                                Arrays.asList(inEdges),
                                                Arrays.asList(outEdges)
                                        );

                                        corrections.put(v1.getCanonicalKmer(), cr);
                                    }

                                    numGapsClosed++;
                                }
                            }
                        }

                        pm.update("reads processed, " + numGaps + " gaps, " + numGapsClosed + " closed");
                    }
                }

                //if (pm.pos() > 100000) { break; }
            }
        }

        CortexRecord[] ccrs = new CortexRecord[corrections.size()];
        int cpos = 0;
        for (CanonicalKmer ck : corrections.keySet()) {
            ccrs[cpos] = corrections.get(ck);
            cpos++;
        }

        CortexGraphWriter cgw = new CortexGraphWriter(out);
        cgw.setHeader(CLEAN.getHeader());

        cpos = 0;

        ProgressMeter pm2 = new ProgressMeterFactory()
                .header("Writing corrected graph")
                .message("records written")
                .maxRecord(CLEAN.getNumRecords())
                .make(log);

        int numAdded = 0;
        for (CortexRecord cr : CLEAN) {
            if (cpos >= ccrs.length || cr.getKmerAsString().compareTo(ccrs[cpos].getKmerAsString()) < 0) {
                cgw.addRecord(cr);
            } else {
                while (cpos < ccrs.length && cr.getKmerAsString().compareTo(ccrs[cpos].getKmerAsString()) >= 0) {
                    cgw.addRecord(ccrs[cpos]);
                    cpos++;

                    numAdded++;
                }
            }

            pm2.update("records written, " + numAdded + "/" + ccrs.length + " added");
        }

        for (int i = cpos; i < ccrs.length; i++) {
            cgw.addRecord(ccrs[i]);
        }

        cgw.close();
    }
}
