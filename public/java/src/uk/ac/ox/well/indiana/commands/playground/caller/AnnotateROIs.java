package uk.ac.ox.well.indiana.commands.playground.caller;

import com.google.common.base.Joiner;
import htsjdk.samtools.*;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.commands.playground.index.KmerIndex;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexKmer;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexRecord;
import uk.ac.ox.well.indiana.utils.progress.ProgressMeter;
import uk.ac.ox.well.indiana.utils.progress.ProgressMeterFactory;

import java.io.File;
import java.io.PrintStream;
import java.util.*;

public class AnnotateROIs extends Module {
    @Argument(fullName="combined", shortName="com", doc="Graph (combined)")
    public CortexGraph COMBINED;

    @Argument(fullName="graph", shortName="g", doc="Graph")
    public CortexGraph GRAPH;

    @Argument(fullName="roi", shortName="r", doc="ROI")
    public CortexGraph ROI;

    @Argument(fullName="contam", shortName="con", doc="Contam")
    public CortexGraph CONTAM;

    @Argument(fullName="bam", shortName="b", doc="BAM")
    public File BAM_FILE;

    @Argument(fullName="child", shortName="c", doc="Child")
    public String CHILD;

    @Argument(fullName="parents", shortName="p", doc="Parents")
    public ArrayList<String> PARENTS;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        int childColor = GRAPH.getColorForSampleName(CHILD);
        Set<Integer> parentColors = new HashSet<>();
        for (String parent : PARENTS) {
            parentColors.add(GRAPH.getColorForSampleName(parent));
        }

        Map<CortexKmer, Map<String, Object>> anns = new HashMap<>();

        KmerIndex ki = new KmerIndex(BAM_FILE, ROI.getKmerSize(), false);
        SamReader sr = SamReaderFactory.make()
                .setOption(SamReaderFactory.Option.INCLUDE_SOURCE_IN_RECORDS, true)
                .setOption(SamReaderFactory.Option.CACHE_FILE_BASED_INDEXES, true)
                .validationStringency(ValidationStringency.SILENT)
                .open(BAM_FILE);

        Map<String, CortexGraph> graphs = new HashMap<>();
        graphs.put("roi", ROI);
        graphs.put("contam", CONTAM);

        Map<String, Set<CortexRecord>> recset = new HashMap<>();
        for (String label : graphs.keySet()) {
            log.info("Processing {} records...", label);

            recset.put(label, new HashSet<>());

            for (CortexRecord cr : graphs.get(label)) {
                recset.get(label).add(cr);
            }
        }

        log.info("Processing dirty records...");
        recset.put("dirty", new HashSet<>());
        for (CortexRecord cr : COMBINED) {
            if (cr.getCoverage(0) > 0 && cr.getCoverage(1) == 0) {
                recset.get("dirty").add(cr);
            }
        }

        log.info("Processing clean records...");
        recset.put("clean", new HashSet<>());
        for (CortexRecord cr : GRAPH) {
            if (recset.get("clean").size() > 10000) {
                break;
            }

            if (cr.getInDegree(childColor) <= 1 && cr.getOutDegree(childColor) <= 1 && cr.getCoverage(childColor) > 0) {
                boolean isInAllSamples = true;

                for (int c = 0; c < cr.getNumColors(); c++) {
                    if (cr.getCoverage(c) == 0) {
                        isInAllSamples = false;
                        break;
                    }
                }

                if (isInAllSamples) {
                    recset.get("clean").add(cr);
                }
            }
        }

        for (String label : recset.keySet()) {
            ProgressMeter pm = new ProgressMeterFactory()
                    .header("Processing " + label + " records")
                    .maxRecord(recset.get(label).size())
                    .make(log);

            for (CortexRecord cr : recset.get(label)) {
                pm.update();
                //Map<String, Object> ann = new HashMap<>();
                //ann.put("")

                List<long[]> chunks = ki.find(cr.getKmerAsBytes());

                List<Byte> quals = new ArrayList<>();
                int qualMin = 0;

                for (long[] c : chunks) {
                    SAMFileSpan sfs = new BAMFileSpan(new Chunk(c[0], c[1]));

                    SAMRecordIterator recs = sr.indexing().iterator(sfs);

                    while (recs.hasNext()) {
                        SAMRecord rec = recs.next();

                        for (int i = 0; i <= rec.getReadLength() - ROI.getKmerSize(); i++) {
                            CortexKmer ck = new CortexKmer(rec.getReadString().substring(i, i + ROI.getKmerSize()));

                            if (cr.getCortexKmer().equals(ck)) {
                                byte minQual = 40;

                                for (int j = i; j < i + ROI.getKmerSize(); j++) {
                                    byte qual = rec.getBaseQualities()[j];

                                    if (qual < minQual) {
                                        minQual = qual;
                                    }
                                }

                                quals.add(minQual);
                                qualMin += minQual;
                            }
                        }
                    }

                    recs.close();
                }

                //log.info("{} {} {}", cr, qualMin / quals.size(), Joiner.on(", ").join(quals));

                out.println(Joiner.on("\t").join(label, cr.getKmerAsString(), cr.getCoverage(0), qualMin / quals.size()));
            }
        }
    }
}
