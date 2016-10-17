package uk.ac.ox.well.indiana.commands.playground.kmerdb;

import com.google.common.base.Joiner;
import org.mapdb.*;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.cortex.collection.CortexCollection;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.*;
import uk.ac.ox.well.indiana.utils.progress.ProgressMeter;
import uk.ac.ox.well.indiana.utils.progress.ProgressMeterFactory;
import uk.ac.ox.well.indiana.utils.sequence.CortexUtils;

import java.io.File;
import java.util.*;

public class BuildDB extends Module {
    @Argument(fullName="dirtyGraph", shortName="d", doc="Dirty graph")
    public CortexCollection DIRTY;

    @Argument(fullName="childColor", shortName="cc", doc="Child color")
    public Integer CHILD_COLOR = 0;

    @Argument(fullName="parentColor", shortName="pc", doc="Parent color")
    public HashSet<Integer> PARENT_COLORS = new HashSet<>();

    @Argument(fullName="coverageLowerLimit", shortName="l", doc="Coverage lower limit")
    public Integer COVERAGE_LOWER_LIMIT = 1;

    @Output
    public File out;

    @Override
    public void execute() {
        if (PARENT_COLORS.isEmpty()) {
            for (int c = 0; c < DIRTY.getNumColors(); c++) {
                if (c != CHILD_COLOR) {
                    PARENT_COLORS.add(c);
                }
            }
        }

        log.info("Colors: (child) {} (parents) {}", CHILD_COLOR, Joiner.on(",").join(PARENT_COLORS));

        CortexHeader ch = new CortexHeader();
        ch.setVersion(6);
        ch.setNumColors(DIRTY.getNumColors());
        ch.setKmerSize(DIRTY.getKmerSize());
        ch.setKmerBits(DIRTY.getKmerBits());

        for (int c = 0; c < DIRTY.getNumColors(); c++) {
            CortexColor cc = new CortexColor();
            cc.setCleanedAgainstGraph(false);
            cc.setCleanedAgainstGraphName("");
            cc.setErrorRate(0);
            cc.setLowCovgKmersRemoved(false);
            cc.setLowCovgSupernodesRemoved(false);
            cc.setSampleName("color" + c);
            cc.setTipClippingApplied(false);
            cc.setTotalSequence(0);

            ch.addColor(cc);
        }

        CortexGraphWriter cgw = new CortexGraphWriter(out);
        cgw.setHeader(ch);

        //List<long[]> binaryKmers = new ArrayList<>();

        ProgressMeter pm = new ProgressMeterFactory()
                .header("Processing graph...")
                .message("records processed")
                .updateRecord(DIRTY.getGraph(0).getNumRecords() / 10)
                .make(log);

        long numNovelRecords = 0;
        for (CortexRecord cr : DIRTY) {
            if (CortexUtils.isNovelKmer(cr, CHILD_COLOR, PARENT_COLORS) && cr.getCoverage(CHILD_COLOR) > COVERAGE_LOWER_LIMIT) {
                //binaryKmers.add(cr.getBinaryKmer());

                cgw.addRecord(cr);

                numNovelRecords++;
            }

            pm.update("records processed (" + numNovelRecords + " novel so far)");
        }

        log.info("  {} total novel records", numNovelRecords);

        cgw.close();

        /*
        if (out.exists()) {
            out.delete();
        }

        DB db = DBMaker
                .fileDB(out)
                .fileMmapEnable()
                .make();

        //Atomic.Integer aKmerSize = db.atomicInteger("kmerSize", DIRTY.getKmerSize()).create();
        //Atomic.Integer aKmerBits = db.atomicInteger("kmerBits", DIRTY.getKmerBits()).create();

        db.atomicInteger("kmerSize", DIRTY.getKmerSize()).create();
        db.atomicInteger("kmerBits", DIRTY.getKmerBits()).create();

        HTreeMap<long[], Map<String, Object>> nkdb = db.hashMap("novelKmers")
                .keySerializer(Serializer.LONG_ARRAY)
                .valueSerializer(Serializer.JAVA)
                .create();

        pm = new ProgressMeterFactory()
                .maxRecord(binaryKmers.size())
                .header("Writing records...")
                .message("records written")
                .make(log);

        for (long[] bk : binaryKmers) {
            nkdb.put(bk, new HashMap<>());

            pm.update();
        }

        db.commit();
        nkdb.close();
        db.close();
        */
    }
}
