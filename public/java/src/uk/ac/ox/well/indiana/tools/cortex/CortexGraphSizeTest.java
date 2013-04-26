package uk.ac.ox.well.indiana.tools.cortex;

import com.carrotsearch.sizeof.RamUsageEstimator;
import uk.ac.ox.well.indiana.tools.Tool;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.io.cortex.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.CortexKmer;
import uk.ac.ox.well.indiana.utils.io.cortex.CortexMap;
import uk.ac.ox.well.indiana.utils.io.cortex.CortexRecord;
import uk.ac.ox.well.indiana.utils.performance.PerformanceUtils;
import uk.ac.ox.well.indiana.utils.sequence.SequenceUtils;

import java.util.*;

public class CortexGraphSizeTest extends Tool {
    @Argument(fullName="cortexGraph", shortName="cg", doc="Cortex graph")
    public CortexMap CORTEX_MAP;

    @Override
    public void execute() {
        String testkmer = SequenceUtils.reverseComplement("GGGAAAGAAAATGGTGAAAACCAAATTATAA");
        log.info("{}: {}", testkmer, CORTEX_MAP.get(testkmer));

        log.info("Sizeof map: {}", RamUsageEstimator.sizeOf(CORTEX_MAP));
    }
}
