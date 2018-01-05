package uk.ac.ox.well.cortexjdk.commands.simulate.generators;

import uk.ac.ox.well.cortexjdk.commands.Module;
import uk.ac.ox.well.cortexjdk.utils.arguments.Argument;
import uk.ac.ox.well.cortexjdk.utils.io.graph.cortex.CortexGraph;
import uk.ac.ox.well.cortexjdk.utils.io.graph.cortex.CortexRecord;
import uk.ac.ox.well.cortexjdk.utils.io.table.TableReader;
import uk.ac.ox.well.cortexjdk.utils.kmer.CanonicalKmer;

import java.io.File;
import java.util.Map;

public class EvaluateSim extends Module {
    @Argument(fullName="graph", shortName="g", doc="Graph")
    public CortexGraph GRAPH;

    @Argument(fullName="sim", shortName="s", doc="Simulation")
    public File SIM;

    @Argument(fullName="calls", shortName="c", doc="Calls")
    public File CALLS;

    @Argument(fullName="childName", shortName="cn", doc="Child name")
    public String CHILD_NAME;

    @Override
    public void execute() {
        int childColor = GRAPH.getColorForSampleName(CHILD_NAME);

        TableReader tr = new TableReader(SIM);

        for (Map<String, String> te : tr) {
            if (!te.get("type").equals("RECOMB")) {
                log.info("{}", te);

                String sleft = te.get("sleft");
                String sright = te.get("sright");
                String oldAllele = te.get("old");
                String newAllele = te.get("new");

                String oldHap = sleft + oldAllele + sright;
                String newHap = sleft + newAllele + sright;

                for (int i = 0; i <= oldHap.length() - GRAPH.getKmerSize(); i++) {
                    String sk = oldHap.substring(i, i + GRAPH.getKmerSize());
                    CanonicalKmer ck = new CanonicalKmer(sk);

                    CortexRecord cr = GRAPH.findRecord(ck);

                    log.info("  old {} {} {}", i, sk, cr.toString(0, 1, childColor));
                }

                for (int i = 0; i <= newHap.length() - GRAPH.getKmerSize(); i++) {
                    String sk = newHap.substring(i, i + GRAPH.getKmerSize());
                    CanonicalKmer ck = new CanonicalKmer(sk);

                    CortexRecord cr = GRAPH.findRecord(ck);

                    log.info("  new {} {} {}", i, sk, cr.toString(0, 1, childColor));
                }

                log.info("");
            }
        }
    }
}
