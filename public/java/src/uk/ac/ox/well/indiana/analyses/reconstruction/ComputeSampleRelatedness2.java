package uk.ac.ox.well.indiana.analyses.reconstruction;

import com.google.common.base.Joiner;
import uk.ac.ox.well.indiana.tools.Tool;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.containers.DataFrame;
import uk.ac.ox.well.indiana.utils.io.cortex.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.CortexKmer;
import uk.ac.ox.well.indiana.utils.io.cortex.CortexMap;
import uk.ac.ox.well.indiana.utils.io.cortex.CortexRecord;
import uk.ac.ox.well.indiana.utils.io.table.TableReader;
import uk.ac.ox.well.indiana.utils.io.table.TableWriter;

import java.io.*;
import java.util.*;

public class ComputeSampleRelatedness2 extends Tool {
    @Argument(fullName="cortexGraph", shortName="cg", doc="Cortex graph")
    public ArrayList<File> CORTEX_GRAPHS;

    @Argument(fullName="startIndex", shortName="si", doc="Starting index")
    public Integer START_INDEX = 0;

    @Argument(fullName="endIndex", shortName="ei", doc="Ending index (-1 to process entire list)")
    public Integer END_INDEX = -1;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        DataFrame<String, String, Float> relatedness = new DataFrame<String, String, Float>(0.0f);

        int startIndex = START_INDEX;
        int endIndex = (END_INDEX == -1) ? CORTEX_GRAPHS.size() : END_INDEX;

        for (int f1 = startIndex; f1 < endIndex; f1++) {
            log.info("Computing relatedness elements for {} ({})", f1, CORTEX_GRAPHS.get(f1));
            CortexMap cortexMap1 = new CortexMap(CORTEX_GRAPHS.get(f1));

            for (int f2 = f1; f2 < CORTEX_GRAPHS.size(); f2++) {
                log.info("\t vs {} ({})", f2, CORTEX_GRAPHS.get(f2));
                CortexMap cortexMap2 = new CortexMap(CORTEX_GRAPHS.get(f2));

                Set<CortexKmer> kmers = new HashSet<CortexKmer>(cortexMap1.keySet().size() + cortexMap2.keySet().size());
                kmers.addAll(cortexMap1.keySet());
                kmers.addAll(cortexMap2.keySet());

                for (CortexKmer kmer : kmers) {
                    if (cortexMap1.containsKey(kmer) && cortexMap2.containsKey(kmer)) {
                        CortexRecord cortexRecord1 = cortexMap1.get(kmer);
                        CortexRecord cortexRecord2 = cortexMap2.get(kmer);

                        for (int color1 = 0; color1 < cortexMap1.getGraph().getNumColors(); color1++) {
                            for (int color2 = 0; color2 < cortexMap2.getGraph().getNumColors(); color2++) {
                                if (cortexRecord1.getCoverage(color1) > 0 && cortexRecord2.getCoverage(color2) > 0) {
                                    String sample1 = cortexMap1.getGraph().getColor(color1).getSampleName();
                                    String sample2 = cortexMap2.getGraph().getColor(color2).getSampleName();

                                    float count = relatedness.get(sample1, sample2);
                                    relatedness.set(sample1, sample2, count + 1.0f);
                                    relatedness.set(sample2, sample1, count + 1.0f);
                                }
                            }
                        }
                    }
                }
            }
        }

        out.println(relatedness);
    }
}
