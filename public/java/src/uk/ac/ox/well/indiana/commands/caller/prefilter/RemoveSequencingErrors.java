package uk.ac.ox.well.indiana.commands.caller.prefilter;

import org.eclipse.collections.impl.factory.Sets;
import org.jgrapht.DirectedGraph;
import org.jgrapht.traverse.TopologicalOrderIterator;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexKmer;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexRecord;
import uk.ac.ox.well.indiana.utils.progress.ProgressMeter;
import uk.ac.ox.well.indiana.utils.progress.ProgressMeterFactory;
import uk.ac.ox.well.indiana.utils.sequence.CortexUtils;
import uk.ac.ox.well.indiana.utils.stoppingconditions.PathClosingStopper;
import uk.ac.ox.well.indiana.utils.stoppingconditions.PathOpeningStopper;
import uk.ac.ox.well.indiana.utils.traversal.AnnotatedEdge;
import uk.ac.ox.well.indiana.utils.traversal.AnnotatedVertex;

import java.util.*;

public class RemoveSequencingErrors extends Module {
    @Argument(fullName="graph", shortName="g", doc="Graph")
    public CortexGraph GRAPH;

    @Argument(fullName="parents", shortName="p", doc="Parents")
    public ArrayList<String> PARENTS;

    @Argument(fullName="child", shortName="c", doc="Child")
    public String CHILD;

    @Argument(fullName="roi", shortName="r", doc="ROIs")
    public CortexGraph ROI;

    @Override
    public void execute() {
        int childColor = GRAPH.getColorForSampleName(CHILD);
        Set<Integer> parentColors = new HashSet<>(GRAPH.getColorsForSampleNames(PARENTS));

        log.info("Colors:");
        log.info(" -   child: {}", GRAPH.getColorForSampleName(CHILD));
        log.info(" - parents: {}", GRAPH.getColorsForSampleNames(PARENTS));

        ProgressMeter pm = new ProgressMeterFactory()
                .header("Finding shared kmers")
                .message("records processed")
                .maxRecord(ROI.getNumRecords())
                .make(log);

        Set<CortexKmer> errorKmers = new HashSet<>();

        for (CortexRecord rr : ROI) {
            if (!errorKmers.contains(rr.getCortexKmer())) {
                DirectedGraph<AnnotatedVertex, AnnotatedEdge> g = CortexUtils.dfs_and(GRAPH, rr.getKmerAsString(), childColor, parentColors, PathOpeningStopper.class);

                for (AnnotatedVertex v : g.vertexSet()) {
                    errorKmers.add(new CortexKmer(v.getKmer()));
                }

                AnnotatedVertex ar = new AnnotatedVertex(rr.getKmerAsString());

                AnnotatedVertex af = null;
                while (g.inDegreeOf(ar) == 1) {
                    if (g.outDegreeOf(ar) > 1) {
                        af = new AnnotatedVertex(ar.getKmer());
                        break;
                    }

                    ar = g.getEdgeSource(g.incomingEdgesOf(ar).iterator().next());
                }

                ar = new AnnotatedVertex(rr.getKmerAsString());
                AnnotatedVertex al = null;
                while (g.outDegreeOf(ar) == 1) {
                    if (g.inDegreeOf(ar) > 1) {
                        al = new AnnotatedVertex(ar.getKmer());
                        break;
                    }

                    ar = g.getEdgeTarget(g.outgoingEdgesOf(ar).iterator().next());
                }

                log.info("{} {}", af, al);

                TopologicalOrderIterator<AnnotatedVertex, AnnotatedEdge> toi = new TopologicalOrderIterator<>(g);

                AnnotatedVertex avFirst = null;
                AnnotatedVertex avLast = null;

                log.info("{}", rr);
                while (toi.hasNext()) {
                    AnnotatedVertex v = toi.next();

                    if (avFirst == null) { avFirst = v; }
                    avLast = v;

                    CortexRecord cr = GRAPH.findRecord(new CortexKmer(v.getKmer()));

                    log.info(" - {} {} {} {} {} {} {} {}",
                            cr.getKmerAsString(),
                            cr.getCoverage(8), cr.getCoverage(0), cr.getCoverage(17),
                            cr.getEdgesAsString(8), cr.getEdgesAsString(0), cr.getEdgesAsString(17),
                            new CortexKmer(v.getKmer()).equals(rr.getCortexKmer())
                    );
                }
                log.info("");

                if (af != null && al != null) {
                    Set<String> nextKmers = CortexUtils.getNextKmers(GRAPH, af.getKmer(), childColor);

                    for (String nextKmer : nextKmers) {
                        DirectedGraph<AnnotatedVertex, AnnotatedEdge> gc = CortexUtils.dfs(GRAPH, null, nextKmer, childColor, parentColors, g, PathClosingStopper.class, 0, true, new HashSet<>());

                        if (gc != null) {
                            TopologicalOrderIterator<AnnotatedVertex, AnnotatedEdge> toi2 = new TopologicalOrderIterator<>(gc);

                            while (toi2.hasNext()) {
                                AnnotatedVertex toi2av = toi2.next();
                                CortexKmer toi2ck = new CortexKmer(toi2av.getKmer());
                                CortexRecord toi2cr = GRAPH.findRecord(toi2ck);
                                log.info("  {} {} {} {} {} {} {} {}", toi2av,
                                        toi2cr.getKmerAsString(),
                                        toi2cr.getCoverage(8), toi2cr.getCoverage(0), toi2cr.getCoverage(17),
                                        toi2cr.getEdgesAsString(8), toi2cr.getEdgesAsString(0), toi2cr.getEdgesAsString(17)
                                );
                            }
                        }
                    }
                }

                log.info("");

                /*
                CortexRecord crFirst = GRAPH.findRecord(new CortexKmer(avFirst.getKmer()));
                CortexRecord crLast  = GRAPH.findRecord(new CortexKmer(avLast.getKmer()));

                if (crFirst != null && crLast != null) {
                    log.info("{} {} {} {} {} {} {}",
                            crFirst.getKmerAsString(),
                            crFirst.getCoverage(8), crFirst.getCoverage(0), crFirst.getCoverage(17),
                            crFirst.getEdgesAsString(8), crFirst.getEdgesAsString(0), crFirst.getEdgesAsString(17)
                    );
                    log.info("{} {} {} {} {} {} {}",
                            crLast.getKmerAsString(),
                            crLast.getCoverage(8), crLast.getCoverage(0), crLast.getCoverage(17),
                            crLast.getEdgesAsString(8), crLast.getEdgesAsString(0), crLast.getEdgesAsString(17)
                    );
                    log.info("");
                }
                */
            }

            pm.update();
        }

        log.info("Found {} shared kmers", errorKmers.size());

        log.info("Writing...");

        /*
        CortexGraphWriter cgw = new CortexGraphWriter(out);
        cgw.setHeader(ROI.getHeader());

        CortexGraphWriter cgo = new CortexGraphWriter(shared_out);
        cgo.setHeader(ROI.getHeader());

        int numKept = 0, numExcluded = 0;
        for (CortexRecord rr : ROI) {
            if (!sharedKmers.contains(rr.getCortexKmer())) {
                cgw.addRecord(rr);
                numKept++;
            } else {
                cgo.addRecord(rr);
                numExcluded++;
            }
        }

        cgw.close();
        cgo.close();

        log.info("  {}/{} ({}%) kept, {}/{} ({}%) excluded",
                numKept,     ROI.getNumRecords(), 100.0f * (float) numKept / (float) ROI.getNumRecords(),
                numExcluded, ROI.getNumRecords(), 100.0f * (float) numExcluded / (float) ROI.getNumRecords()
        );
        */
    }
}
