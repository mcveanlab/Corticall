package uk.ac.ox.well.cortexjdk.commands.simulate;

import htsjdk.samtools.SAMRecord;
import uk.ac.ox.well.cortexjdk.commands.Module;
import uk.ac.ox.well.cortexjdk.utils.alignment.reference.IndexedReference;
import uk.ac.ox.well.cortexjdk.utils.arguments.Argument;
import uk.ac.ox.well.cortexjdk.utils.io.graph.cortex.CortexGraph;
import uk.ac.ox.well.cortexjdk.utils.io.graph.cortex.CortexRecord;
import uk.ac.ox.well.cortexjdk.utils.io.table.TableReader;
import uk.ac.ox.well.cortexjdk.utils.kmer.CanonicalKmer;
import uk.ac.ox.well.cortexjdk.utils.sequence.SequenceUtils;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class SimToVCF extends Module {
    @Argument(fullName="sim", shortName="s", doc="Simulation")
    public File SIM;

    @Argument(fullName="background", shortName="b", doc="Background", required=false)
    public HashMap<Integer, IndexedReference> BACKGROUNDS;

    @Override
    public void execute() {
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

                if (te.get("type").equals("LARGE_INV")) {
                    oldHap = sleft + oldAllele + sright;
                }

                IndexedReference ref = BACKGROUNDS.get(Integer.valueOf(te.get("parent")));

                List<SAMRecord> srs = sort(threshold(ref.align(oldHap.toUpperCase()), 0));
                for (SAMRecord sr : srs) {
                    if (sr.getMappingQuality() > 0) {
                        log.info("  sr: {}", sr.getSAMString().trim());
                    }
                }

                //log.info("");
            }
        }
    }

    private List<SAMRecord> threshold(List<SAMRecord> srs, int mqThreshold) {
        List<SAMRecord> remaining = new ArrayList<>();

        srs.forEach(r -> { if (r.getMappingQuality() > mqThreshold) remaining.add(r); });

        return remaining;
    }

    private List<SAMRecord> sort(List<SAMRecord> srs) {
        srs.sort((o1, o2) -> {
            if (o1.getMappingQuality() != o2.getMappingQuality()) {
                return o1.getMappingQuality() > o2.getMappingQuality() ? -1 : 1;
            }

            if (!o1.getIntegerAttribute("NM").equals(o2.getIntegerAttribute("NM"))) {
                return o1.getIntegerAttribute("NM") < o2.getIntegerAttribute("NM") ? -1 : 1;
            }

            return 0;
        });

        return srs;
    }
}
