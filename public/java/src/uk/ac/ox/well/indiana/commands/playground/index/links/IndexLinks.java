package uk.ac.ox.well.indiana.commands.playground.index.links;

import com.google.common.base.Joiner;
import htsjdk.samtools.reference.FastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import org.mapdb.BTreeMap;
import org.mapdb.DB;
import org.mapdb.DBMaker;
import org.mapdb.Serializer;
import org.mapdb.serializer.SerializerArrayTuple;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexBinaryKmer;
import uk.ac.ox.well.indiana.utils.io.cortex.links.CortexGraphLinks;
import uk.ac.ox.well.indiana.utils.io.cortex.links.CortexJunctionsRecord;
import uk.ac.ox.well.indiana.utils.io.cortex.links.CortexLinks;
import uk.ac.ox.well.indiana.utils.io.cortex.links.CortexLinksRecord;
import uk.ac.ox.well.indiana.utils.progress.ProgressMeter;
import uk.ac.ox.well.indiana.utils.progress.ProgressMeterFactory;
import uk.ac.ox.well.indiana.utils.sequence.SequenceUtils;

import java.io.File;
import java.util.Arrays;
import java.util.NavigableSet;
import java.util.Set;

public class IndexLinks extends Module {
    @Argument(fullName="links", shortName="l", doc="Links")
    public CortexGraphLinks LINKS;

    @Override
    public void execute() {
        File dbFile = new File(LINKS.getCortexLinksFile().getAbsolutePath() + ".linkdb");

        DB db = DBMaker
                .fileDB(dbFile)
                .fileMmapEnable()
                .make();

        NavigableSet<Object[]> linkIndex = db.treeSet("links")
                .serializer(new SerializerArrayTuple(Serializer.STRING, Serializer.STRING))
                .create();

        ProgressMeter pm = new ProgressMeterFactory()
                .header("Processing links")
                .message("links processed")
                .maxRecord(LINKS.getNumKmersWithLinks())
                .make(log);

        for (CortexLinksRecord clr : LINKS) {
            for (CortexJunctionsRecord cjr : clr.getJunctions()) {
                linkIndex.add(new Object[]{clr.getKmerAsString(), cjr.toString()});
            }
            db.commit();

            pm.update();
        }

        db.close();

        /*
        Set subset = linkIndex.subSet(
                new Object[]{example},        // lower interval bound
                new Object[]{example, null});

        for (Object o : subset) {
            Object[] oa = (Object[]) o;
            String kmer = (String) oa[0];
            String linkText = (String) oa[1];

            log.info("{} {} {}", example, kmer, linkText);
        }
        */
    }
}
