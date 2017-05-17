package uk.ac.ox.well.indiana.commands.playground.nahr;

import com.google.api.client.repackaged.com.google.common.base.Joiner;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.StringUtil;
import org.apache.commons.lang3.tuple.Pair;
import org.jgrapht.DirectedGraph;
import org.jgrapht.graph.DefaultDirectedGraph;
import org.jgrapht.graph.DefaultEdge;
import org.mapdb.BTreeMap;
import org.mapdb.DB;
import org.mapdb.DBMaker;
import org.mapdb.Serializer;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.commands.playground.index.KmerIndex;
import uk.ac.ox.well.indiana.utils.alignment.kmer.KmerLookup;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.exceptions.IndianaException;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexBinaryKmer;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexKmer;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexRecord;
import uk.ac.ox.well.indiana.utils.sequence.SequenceUtils;
import uk.ac.ox.well.indiana.utils.stoppingconditions.PartitionStopper;
import uk.ac.ox.well.indiana.utils.stoppingconditions.ShoreStopper;
import uk.ac.ox.well.indiana.utils.stoppingconditions.ShoreStopper2;
import uk.ac.ox.well.indiana.utils.traversal.*;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.*;

public class NahrExperiment extends Module {
    @Argument(fullName="reference", shortName="R", doc="Genome")
    public File REF;

    @Argument(fullName="graph", shortName="g", doc="Graph")
    public CortexGraph GRAPH;

    @Argument(fullName="parents", shortName="p", doc="Parents")
    public ArrayList<String> PARENTS;

    @Argument(fullName="child", shortName="c", doc="Child")
    public String CHILD;

    @Argument(fullName="roi", shortName="r", doc="ROI")
    public CortexGraph ROI;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        IndexedFastaSequenceFile ref;
        try {
            ref = new IndexedFastaSequenceFile(REF);
        } catch (FileNotFoundException e) {
            throw new IndianaException("Reference not found");
        }

        DB db = DBMaker
                .fileDB(REF.getAbsoluteFile() + ".kmerdb")
                .fileMmapEnable()
                .readOnly()
                .make();

        BTreeMap<long[], int[]> kmerIndex = db.treeMap("index")
                .keySerializer(Serializer.LONG_ARRAY)
                .valueSerializer(Serializer.INT_ARRAY)
                .counterEnable()
                .open();

        TraversalEngine e = new TraversalEngineFactory()
                .traversalColor(GRAPH.getColorForSampleName(CHILD))
                .joiningColors(GRAPH.getColorsForSampleNames(PARENTS))
                .recruitmentColors(GRAPH.getColorsForSampleNames(PARENTS))
                .traversalDirection(TraversalEngineConfiguration.TraversalDirection.BOTH)
                .combinationOperator(TraversalEngineConfiguration.GraphCombinationOperator.OR)
                .rois(ROI)
                .graph(GRAPH)
                .stopper(PartitionStopper.class)
                .make();

        String[] recombPartners = {
                ref.getSubsequenceAt("Pf3D7_01_v3", 29510, 38156).getBaseString(),
                SequenceUtils.reverseComplement(ref.getSubsequenceAt("Pf3D7_02_v3", 915560, 926711).getBaseString()),
        };

        DirectedGraph<String, DefaultEdge> d = new DefaultDirectedGraph<>(DefaultEdge.class);

        boolean keepGoing = true;
        int pi = 0;
        int pos[] = {0, 0};
        boolean recomb = false;

        String sk = recombPartners[pi].substring(pos[pi], pos[pi] + GRAPH.getKmerSize());

        StringBuilder haplotype = new StringBuilder(sk.substring(0, GRAPH.getKmerSize() - 1));
        StringBuilder copying = new StringBuilder(StringUtil.repeatCharNTimes('0', sk.length() - 1));

        do {
            int[] l = kmerIndex.get(new CortexBinaryKmer(sk.getBytes()).getBinaryKmer());

            List<String> loci = new ArrayList<>();

            if (l != null) {
                for (int i = 0; i < l.length - 1; i += 2) {
                    String locus = ref.getSequenceDictionary().getSequence(l[i]).getSequenceName() + ":" + l[i+1];
                    loci.add(locus);
                }
            } else {
                loci.add(".");
                recomb = true;
            }

            if (recomb) {
                int q0 = recombPartners[0].indexOf(sk, pos[0]);
                int q1 = recombPartners[1].indexOf(sk, pos[1]);

                if (q0 >= 0 && q1 == -1) {
                    pi = 0;
                    pos[pi] = q0;
                    recomb = false;
                } else if (q0 == -1 && q1 >= 0) {
                    pi = 1;
                    pos[pi] = q1;
                    recomb = false;
                } else {
                    int a = 0;
                }
            }

            log.info("{} {} {} {}", pi, pos, sk, Joiner.on(",").join(loci));
            haplotype.append(sk.substring(sk.length() - 1, sk.length()));
            if (loci.contains(".")) {
                copying.append(".");
            } else {
                copying.append(String.valueOf(pi));
            }

            pos[pi]++;

            String nk = recombPartners[pi].substring(pos[pi], pos[pi] + GRAPH.getKmerSize());

            Set<CortexVertex> nvs = e.getNextVertices(sk);
            Set<String> nvks = new HashSet<>();

            for (CortexVertex nv : nvs) { nvks.add(nv.getSk()); }

            if (nvks.contains(nk)) {
                sk = nk;
            } else {
                if (nvks.size() == 1) {
                    sk = nvks.iterator().next();
                } else {
                    sk = recombPartners[pi].substring(pos[pi] + 1, pos[pi] + 1 + GRAPH.getKmerSize());
                }
            }
        } while (haplotype.length() < 4000 && pos[0] + 1 + GRAPH.getKmerSize() < recombPartners[0].length() && pos[1] + 1 + GRAPH.getKmerSize() < recombPartners[1].length());

        log.info("hap: {}", haplotype);
        log.info("cop: {}", copying);

        db.close();

        /*
        Map<CortexKmer, CortexRecord> rois = new TreeMap<>();
        for (CortexRecord cr : ROI) {
            rois.put(cr.getCortexKmer(), cr);
        }

        Set<CortexKmer> seen = new HashSet<>();
        int group = 0;
        for (CortexKmer ck : rois.keySet()) {
            if (!seen.contains(ck)) {
                DirectedGraph<CortexVertex, CortexEdge> d = e.dfs(ck.getKmerAsString());

                String contig = e.getContig(d, ck.getKmerAsString(), GRAPH.getColorForSampleName(CHILD));
                StringBuilder qual = new StringBuilder(StringUtil.repeatCharNTimes('I', GRAPH.getKmerSize() - 1));

                for (int i = 0; i <= contig.length() - GRAPH.getKmerSize(); i++) {
                    String sk = contig.substring(i, i + GRAPH.getKmerSize());

                    int[] l = kmerIndex.get(new CortexBinaryKmer(sk.getBytes()).getBinaryKmer());

                    if (l != null) {
                        String seqid = ref.getSequenceDictionary().getSequence(l[0]).getSequenceName();
                        int pos = l[1] + 1;

                        qual.append("I");

                        log.info("{} {} {} {} {}",
                                group,
                                i,
                                sk,
                                seqid,
                                pos
                        );
                    } else {
                        qual.append("!");

                        log.info("{} {} {} {} {}",
                                group,
                                i,
                                sk,
                                ".",
                                "."
                        );
                    }
                }

                out.println("@HWI" + group);
                out.println(contig);
                out.println("+");
                out.println(qual);

                for (CortexVertex cv : d.vertexSet()) {
                    if (!seen.contains(cv.getCk())) {
                        seen.add(cv.getCk());
                    }
                }

                group++;
            }
        }
        */


        /*
        String start = "ATGGTGACGCAAAGTAGTGGTGGGGGTGCTGCTGGTAGTAGTGGTGA";
        CortexBinaryKmer cbk = new CortexBinaryKmer(start.getBytes());
        int[] l = kmerIndex.get(cbk.getBinaryKmer());
        String seqName = ref.getSequenceDictionary().getSequence(l[0]).getSequenceName();
        int pos = l[1];

        Map<Integer, Integer> donorIds = new HashMap<>();
        int maxDonorId = 0;
        int curDonorId = 0;
        donorIds.put(l[0], curDonorId);

        String sk = ref.getSubsequenceAt(seqName, pos + 1, pos + GRAPH.getKmerSize()).getBaseString();
        StringBuilder haplotype = new StringBuilder(sk);
        StringBuilder donor = new StringBuilder(StringUtil.repeatCharNTimes('0', sk.length()));

        boolean keepGoing = true;
        do {
            pos++;

            log.info("{} {} {} {}", seqName, pos + 1, sk, haplotype.length());

            String nextKmerFromRef = ref.getSubsequenceAt(seqName, pos + 1, pos + GRAPH.getKmerSize()).getBaseString();
            Set<CortexVertex> nextVerticesFromGraph = e.getNextVertices(sk);
            Set<String> nextKmersFromGraph = new HashSet<>();
            for (CortexVertex nextVertex : nextVerticesFromGraph) {
                nextKmersFromGraph.add(nextVertex.getSk());
            }

            if (nextKmersFromGraph.contains(nextKmerFromRef)) {
                donor.append(curDonorId);
                haplotype.append(nextKmerFromRef.substring(nextKmerFromRef.length() - 1, nextKmerFromRef.length()));

                sk = nextKmerFromRef;
            } else if (nextKmersFromGraph.size() == 1) {
                String nk = nextKmersFromGraph.iterator().next();

                DirectedGraph<CortexVertex, CortexEdge> d = e.dfs(nk);
                if (d != null) {
                    String contig = e.getContig(d, nk, GRAPH.getColorForSampleName(CHILD));

                    for (int i = 0; i <= contig.length() - GRAPH.getKmerSize(); i++) {
                        String ak = contig.substring(i, i + GRAPH.getKmerSize());

                        log.info("{} {} {} {}", ".", i, ak, haplotype.length());

                        int[] ki = kmerIndex.get(new CortexBinaryKmer(ak.getBytes()).getBinaryKmer());
                        if (ki != null) {
                            if (!donorIds.containsKey(ki[0])) {
                                maxDonorId++;
                                donorIds.put(ki[0], maxDonorId);
                            }
                            curDonorId = donorIds.get(ki[0]);
                            donor.append(curDonorId);

                            seqName = ref.getSequenceDictionary().getSequence(ki[0]).getSequenceName();
                            pos = ki[1];
                        } else {
                            donor.append(".");
                        }

                        haplotype.append(nextKmerFromRef.substring(nextKmerFromRef.length() - 1, nextKmerFromRef.length()));

                        sk = ak;
                    }
                } else {
                    log.info("{} {} {} {}", ".", 0, nk, haplotype.length());

                    int[] ki = kmerIndex.get(new CortexBinaryKmer(nk.getBytes()).getBinaryKmer());
                    if (ki != null) {
                        if (!donorIds.containsKey(ki[0])) {
                            maxDonorId++;
                            donorIds.put(ki[0], maxDonorId);
                        }
                        curDonorId = donorIds.get(ki[0]);
                        donor.append(curDonorId);

                        seqName = ref.getSequenceDictionary().getSequence(ki[0]).getSequenceName();
                        pos = ki[1];
                    } else {
                        donor.append(".");
                    }

                    haplotype.append(nk.substring(nk.length() - 1, nk.length()));

                    sk = nk;
                }
            } else {
                TraversalEngine f = new TraversalEngineFactory()
                        .traversalColor(GRAPH.getColorForSampleName(CHILD))
                        .joiningColors(GRAPH.getColorsForSampleNames(PARENTS))
                        .recruitmentColors(GRAPH.getColorsForSampleNames(PARENTS))
                        .traversalDirection(TraversalEngineConfiguration.TraversalDirection.FORWARD)
                        .combinationOperator(TraversalEngineConfiguration.GraphCombinationOperator.OR)
                        .rois(ROI)
                        .graph(GRAPH)
                        .stopper(ShoreStopper2.class)
                        .make();

                Map<String, DirectedGraph<CortexVertex, CortexEdge>> ds = new HashMap<>();

                for (String nk : nextKmersFromGraph) {
                    DirectedGraph<CortexVertex, CortexEdge> d = f.dfs(nk);

                    if (d != null) { ds.put(nk, d); }
                }

                if (ds.size() == 1) {
                    String nk = ds.keySet().iterator().next();
                    DirectedGraph<CortexVertex, CortexEdge> d = ds.get(nk);

                    String contig = e.getContig(d, nk, GRAPH.getColorForSampleName(CHILD));

                    for (int i = 0; i <= contig.length() - GRAPH.getKmerSize(); i++) {
                        String ak = contig.substring(i, i + GRAPH.getKmerSize());

                        int[] ki = kmerIndex.get(new CortexBinaryKmer(ak.getBytes()).getBinaryKmer());
                        if (ki != null) {
                            if (!donorIds.containsKey(ki[0])) {
                                maxDonorId++;
                                donorIds.put(ki[0], maxDonorId);
                            }
                            curDonorId = donorIds.get(ki[0]);
                            donor.append(curDonorId);

                            seqName = ref.getSequenceDictionary().getSequence(ki[0]).getSequenceName();
                            pos = ki[1];
                        } else {
                            donor.append(".");
                        }

                        log.info("{} {} {} {}", ".", seqName, pos, ak, haplotype.length());

                        haplotype.append(nextKmerFromRef.substring(nextKmerFromRef.length() - 1, nextKmerFromRef.length()));

                        sk = ak;
                    }
                } else {
                    log.info("{}", haplotype);
                    log.info("{}", donor);

                    log.info("");
                }
            }

            if (cvs.size() == 1) {

            }
        } while (keepGoing);
        */
    }
}
