package uk.ac.ox.well.indiana.commands.gg;

import org.jgrapht.DirectedGraph;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.exceptions.IndianaException;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexKmer;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexRecord;
import uk.ac.ox.well.indiana.utils.io.utils.LineReader;
import uk.ac.ox.well.indiana.utils.sequence.CortexUtils;

import java.io.File;
import java.io.PrintStream;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

public class FilterNovels extends Module {
    @Argument(fullName="clean", shortName="c", doc="Graph (clean)")
    public CortexGraph CLEAN;

    @Argument(fullName="dirty", shortName="d", doc="Graph (dirty)")
    public CortexGraph DIRTY;

    @Argument(fullName="novelKmers", shortName="n", doc="Novel kmers")
    public CortexGraph NOVEL_KMERS;

    @Argument(fullName="contaminants", shortName="x", doc="Contaminants")
    public CortexGraph REJECTED_KMERS;

    @Argument(fullName="lowerThreshold", shortName="l", doc="Lower coverage threshold")
    public String LOWER_THRESHOLD = "10";

    @Argument(fullName="upperThreshold", shortName="u", doc="Upper coverage threshold")
    public String UPPER_THRESHOLD = "200";

    @Argument(fullName="performAdjacencyPass", shortName="p", doc="Perform adjacency pass")
    public Boolean PERFORM_ADJACENCY_PASS = false;

    @Argument(fullName="excludeDirtyRecords", shortName="edr", doc="Exclude dirty records")
    public Boolean EXCLUDE_DIRTY_RECORDS = false;

    @Output
    public PrintStream out;

    @Output(fullName="rejectedOut", shortName="ro", doc="Rejected out")
    public PrintStream rout;

    private int loadThreshold(String threshold) {
        File f = new File(threshold);
        if (f.exists()) {
            LineReader lr = new LineReader(f);
            threshold = lr.getNextRecord();
        }

        return Integer.valueOf(threshold);
    }

    @Override
    public void execute() {
        int lowerThreshold = loadThreshold(LOWER_THRESHOLD);
        int upperThreshold = loadThreshold(UPPER_THRESHOLD);

        if (lowerThreshold > 10) { lowerThreshold = 10; }

//        if (REJECTED_KMERS.getNumRecords() > 100000) {
//            throw new IndianaException("Too many contaminants (" + REJECTED_KMERS.getNumRecords() + ")");
//        }

        log.info("Filtering novel kmers...");
        log.info("  {} kmers to start with", NOVEL_KMERS.getNumRecords());

        log.info("Examining coverage...");
        Map<CortexKmer, CortexRecord> remainingRecords = new HashMap<CortexKmer, CortexRecord>();
        int coverageOutliers = 0;
        for (CortexRecord cr : NOVEL_KMERS) {
            if (cr.getCoverage(0) < lowerThreshold || cr.getCoverage(0) > upperThreshold) {
                coverageOutliers++;
            } else {
                remainingRecords.put(cr.getCortexKmer(), cr);
            }
        }

        log.info("  {} kmers outside coverage limits {} and {}", coverageOutliers, lowerThreshold, upperThreshold);
        log.info("  {} kmers remaining", remainingRecords.size());

        Set<CortexKmer> contaminatingKmers = new HashSet<CortexKmer>();

        log.info("Exploring contaminants...");
        log.info("  {} to start with", REJECTED_KMERS.getNumRecords());
        for (CortexRecord cr : REJECTED_KMERS) {
            if (!contaminatingKmers.contains(cr.getCortexKmer())) {
                DirectedGraph<AnnotatedVertex, AnnotatedEdge> dfs = CortexUtils.dfs(CLEAN, DIRTY, cr.getKmerAsString(), 0, null, ContaminantStopper.class);

                for (AnnotatedVertex rv : dfs.vertexSet()) {
                    CortexRecord rr = CLEAN.findRecord(new CortexKmer(rv.getKmer()));
                    if (rr == null) {
                        rr = DIRTY.findRecord(new CortexKmer(rv.getKmer()));
                    }

                    CortexKmer rk = rr.getCortexKmer();

                    if (rr.getCoverage(0) > 0 && rr.getCoverage(1) == 0 && rr.getCoverage(2) == 0) {
                        if (remainingRecords.containsKey(rk) && !contaminatingKmers.contains(rk)) {
                            contaminatingKmers.add(rr.getCortexKmer());
                        }
                    }
                }

                contaminatingKmers.add(cr.getCortexKmer());
            }
        }
        log.info("  {} contaminants found", contaminatingKmers.size());

        Set<CortexKmer> orphanedKmers = new HashSet<CortexKmer>();

        log.info("Finding orphans...");
        for (CortexRecord cr : remainingRecords.values()) {
            if (remainingRecords.containsKey(cr.getCortexKmer()) && !contaminatingKmers.contains(cr.getCortexKmer()) && !orphanedKmers.contains(cr.getCortexKmer())) {
                String stretch = CortexUtils.getSeededStretch(CLEAN, DIRTY, cr.getKmerAsString(), 0, true);

                Set<CortexKmer> novelKmers = new HashSet<CortexKmer>();
                boolean isOrphaned = true;

                for (int i = 0; i <= stretch.length() - CLEAN.getKmerSize(); i++) {
                    CortexKmer ak = new CortexKmer(stretch.substring(i, i + CLEAN.getKmerSize()));
                    CortexRecord ar = CLEAN.findRecord(ak);
                    if (ar == null) {
                        ar = DIRTY.findRecord(ak);
                    }

                    if (ar != null) {
                        if (ar.getCoverage(0) > 0 && ar.getCoverage(1) == 0 && ar.getCoverage(2) == 0) {
                            if (remainingRecords.containsKey(ak) && !contaminatingKmers.contains(ak) && !orphanedKmers.contains(ak)) {
                                novelKmers.add(ak);
                            }
                        } else if (ar.getCoverage(1) > 0 || ar.getCoverage(2) > 0) {
                            isOrphaned = false;
                            break;
                        }
                    }
                }

                if (isOrphaned) {
                    orphanedKmers.addAll(novelKmers);
                }
            }
        }
        log.info("  {} orphaned kmers", orphanedKmers.size());

        log.info("Looking for overcleaning...");
        Set<CortexKmer> overcleanedKmers = new HashSet<CortexKmer>();
        for (CortexRecord cr : remainingRecords.values()) {
            if (remainingRecords.containsKey(cr.getCortexKmer()) && !contaminatingKmers.contains(cr.getCortexKmer()) && !orphanedKmers.contains(cr.getCortexKmer())) {
                Set<CortexKmer> kmers = new HashSet<CortexKmer>();
                boolean hasTaintedNovelKmers = false;

                String stretch = CortexUtils.getNovelStretch(CLEAN, cr.getKmerAsString(), 0, true);
                for (int i = 0; i <= stretch.length() - CLEAN.getKmerSize(); i++) {
                    String sk = stretch.substring(i , i + CLEAN.getKmerSize());
                    CortexKmer ck = new CortexKmer(sk);
                    CortexRecord cleanRecord = CLEAN.findRecord(ck);
                    CortexRecord dirtyRecord = DIRTY.findRecord(ck);

                    if (cleanRecord != null && dirtyRecord != null && CortexUtils.isNovelKmer(cleanRecord, 0) && !CortexUtils.isNovelKmer(dirtyRecord, 0)) {
                        hasTaintedNovelKmers = true;
                    }

                    if (remainingRecords.containsKey(ck)) {
                        kmers.add(ck);
                    }
                }

                if (hasTaintedNovelKmers) {
                    overcleanedKmers.addAll(kmers);
                }
            }
        }
        log.info("  {} overcleaned kmers", overcleanedKmers.size());

        log.info("Looking for adjacent kmers to discard...");
        Set<CortexKmer> adjacentToRejection = new HashSet<CortexKmer>();
        if (PERFORM_ADJACENCY_PASS) {
            boolean addedStuff;
            int passes = 0;
            do {
                addedStuff = false;
                for (CortexRecord cr : remainingRecords.values()) {
                    CortexKmer ck = cr.getCortexKmer();

                    if (!contaminatingKmers.contains(ck) && !orphanedKmers.contains(ck) && !adjacentToRejection.contains(ck) && !overcleanedKmers.contains(ck)) {
                        String sk = cr.getKmerAsString();

                        Set<String> adjKmers = new HashSet<String>();
                        adjKmers.addAll(CortexUtils.getPrevKmers(CLEAN, sk, 0));
                        adjKmers.addAll(CortexUtils.getNextKmers(CLEAN, sk, 0));

                        for (String adjKmer : adjKmers) {
                            CortexKmer cAdjKmer = new CortexKmer(adjKmer);

                            if (contaminatingKmers.contains(cAdjKmer) || orphanedKmers.contains(cAdjKmer) || adjacentToRejection.contains(cAdjKmer) || overcleanedKmers.contains(cAdjKmer)) {
                                adjacentToRejection.add(ck);

                                addedStuff = true;
                                break;
                            }
                        }
                    }
                }

                log.info("  pass {}, {} kmers", passes, adjacentToRejection.size());

                passes++;
            } while (addedStuff);
        }
        log.info("  {} adjacent kmers", adjacentToRejection.size());


        int contams = 0, orphans = 0, adj = 0, overcleaned = 0, dirty = 0, count = 0;
        for (CortexKmer ck : remainingRecords.keySet()) {
            if (contaminatingKmers.contains(ck)) {
                rout.println(">contams" + contams + "\n" + ck.getKmerAsString());
                contams++;
            } else if (orphanedKmers.contains(ck)) {
                rout.println(">orphans" + orphans + "\n" + ck.getKmerAsString());
                orphans++;
            } else if (adjacentToRejection.contains(ck)) {
                rout.println(">adj" + adj + "\n" + ck.getKmerAsString());
                adj++;
            } else if (overcleanedKmers.contains(ck)) {
                rout.println(">overcleaned" + overcleaned + "\n" + ck.getKmerAsString());
                overcleaned++;
            } else if (EXCLUDE_DIRTY_RECORDS && CLEAN.findRecord(ck) == null) {
                rout.println(">dirty" + dirty + "\n" + ck.getKmerAsString());
                dirty++;
            } else {
                out.println(">" + count + "\n" + ck.getKmerAsString());
                count++;
            }
        }

        log.info("Results");
        log.info("  {} before filtering", NOVEL_KMERS.getNumRecords());
        log.info("  {} rejected for coverage", coverageOutliers);
        log.info("  {} rejected for contamination", contams);
        log.info("  {} rejected for orphan status", orphans);
        log.info("  {} rejected for overcleaning", overcleaned);
        log.info("  {} rejected for adjacency", adj);
        log.info("  {} rejected for dirtiness", dirty);
        log.info("  {} remaining", count);
    }
}
