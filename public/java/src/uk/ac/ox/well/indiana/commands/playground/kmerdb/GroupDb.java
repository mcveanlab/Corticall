package uk.ac.ox.well.indiana.commands.playground.kmerdb;

import org.jgrapht.DirectedGraph;
import org.jgrapht.graph.DefaultDirectedGraph;
import org.jgrapht.graph.DefaultEdge;
import org.mapdb.DB;
import org.mapdb.DBMaker;
import org.mapdb.HTreeMap;
import org.mapdb.Serializer;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexKmer;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexRecord;
import uk.ac.ox.well.indiana.utils.progress.ProgressMeter;
import uk.ac.ox.well.indiana.utils.progress.ProgressMeterFactory;
import uk.ac.ox.well.indiana.utils.sequence.CortexUtils;
import uk.ac.ox.well.indiana.utils.sequence.SequenceUtils;

import java.io.File;
import java.util.HashMap;
import java.util.Map;

public class GroupDb extends Module {
    @Argument(fullName="graph", shortName="g", doc="Novel kmer graph")
    public CortexGraph NOVEL_GRAPH;

    @Output
    public File out;

    @Override
    public void execute() {
        DirectedGraph<String, DefaultEdge> dg = new DefaultDirectedGraph<>(DefaultEdge.class);

        ProgressMeter pm = new ProgressMeterFactory()
                .maxRecord(NOVEL_GRAPH.getNumRecords())
                .header("Processing kmers...")
                .make(log);

        for (CortexRecord cr : NOVEL_GRAPH) {
            linkAdjacentKmers(dg, cr);

            pm.update();
        }

        pm = new ProgressMeterFactory()
                .maxRecord(NOVEL_GRAPH.getNumRecords())
                .header("Building contigs from kmers...")
                .make(log);

        Map<CortexKmer, Integer> kmerGroups = new HashMap<>();
        Map<CortexKmer, String> contigs = new HashMap<>();

        int id = 0;
        for (CortexRecord cr : NOVEL_GRAPH) {
            if (!kmerGroups.containsKey(cr.getCortexKmer())) {
                String contig = buildContig(dg, kmerGroups, cr);

                for (int i = 0; i <= contig.length() - NOVEL_GRAPH.getKmerSize(); i++) {
                    CortexKmer ck = new CortexKmer(contig.substring(i, i + NOVEL_GRAPH.getKmerSize()));

                    contigs.put(ck, contig);
                    kmerGroups.put(ck, id);
                }

                id++;
            }

            pm.update();
        }

        if (out.exists()) { out.delete(); }

        DB db = DBMaker
                .fileDB(out)
                .fileMmapEnable()
                .make();

        db.atomicInteger("kmerSize", NOVEL_GRAPH.getKmerSize()).create();
        db.atomicInteger("kmerBits", NOVEL_GRAPH.getKmerBits()).create();

        HTreeMap<long[], Map<String, Object>> nkdb = db.hashMap("novelKmers")
                .keySerializer(Serializer.LONG_ARRAY)
                .valueSerializer(Serializer.JAVA)
                .create();

        pm = new ProgressMeterFactory()
                .maxRecord(NOVEL_GRAPH.getNumRecords())
                .header("Writing novel kmer database...")
                .message("records written")
                .make(log);

        for (CortexRecord cr : NOVEL_GRAPH) {
            long[] bk = CortexUtils.encodeBinaryKmer(cr.getKmerAsBytes());

            Map<String, Object> h = new HashMap<>();
            h.put("group", kmerGroups.get(cr.getCortexKmer()));
            h.put("contig", contigs.get(cr.getCortexKmer()));

            nkdb.put(bk, h);

            pm.update();
        }

        db.commit();
        nkdb.close();
        db.close();
    }

    private String buildContig(DirectedGraph<String, DefaultEdge> dg, Map<CortexKmer, Integer> usedKmers, CortexRecord cr) {
        CortexKmer ck = cr.getCortexKmer();
        String fw = ck.getKmerAsString();

        StringBuilder contig = new StringBuilder(fw);

        while (dg.outDegreeOf(fw) == 1) {
            fw = dg.getEdgeTarget(dg.outgoingEdgesOf(fw).iterator().next());
            contig.append(fw.charAt(fw.length() - 1));

            usedKmers.put(new CortexKmer(fw), -1);
        }

        fw = ck.getKmerAsString();

        while (dg.inDegreeOf(fw) == 1) {
            fw = dg.getEdgeSource(dg.incomingEdgesOf(fw).iterator().next());
            contig.insert(0, fw.charAt(0));

            usedKmers.put(new CortexKmer(fw), -1);
        }

        return contig.toString();
    }

    private void linkAdjacentKmers(DirectedGraph<String, DefaultEdge> dg, CortexRecord cr) {
        String thisFw = cr.getKmerAsString();

        dg.addVertex(thisFw);

        for (String nextFw : CortexUtils.getNextKmers(cr, thisFw, 0)) {
            dg.addVertex(nextFw);
            dg.addEdge(thisFw, nextFw);
        }

        for (String prevFw : CortexUtils.getPrevKmers(cr, thisFw, 0)) {
            dg.addVertex(prevFw);
            dg.addEdge(prevFw, thisFw);
        }

        String thisRc = SequenceUtils.reverseComplement(thisFw);

        dg.addVertex(thisRc);

        for (String nextRc : CortexUtils.getNextKmers(cr, thisRc, 0)) {
            dg.addVertex(nextRc);
            dg.addEdge(thisRc, nextRc);
        }

        for (String prevRc : CortexUtils.getPrevKmers(cr, thisRc, 0)) {
            dg.addVertex(prevRc);
            dg.addEdge(prevRc, thisRc);
        }
    }
}
