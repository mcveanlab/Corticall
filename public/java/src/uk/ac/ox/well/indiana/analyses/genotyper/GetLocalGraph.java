package uk.ac.ox.well.indiana.analyses.genotyper;

import org.jgrapht.DirectedGraph;
import org.jgrapht.ext.DOTExporter;
import org.jgrapht.graph.DefaultEdge;
import uk.ac.ox.well.indiana.tools.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.assembly.CortexEdge;
import uk.ac.ox.well.indiana.utils.assembly.CortexGraphWalker;
import uk.ac.ox.well.indiana.utils.assembly.CortexKmerIDProvider;
import uk.ac.ox.well.indiana.utils.io.cortex.CortexKmer;
import uk.ac.ox.well.indiana.utils.io.cortex.CortexMap;

import java.io.PrintStream;
import java.io.PrintWriter;

public class GetLocalGraph extends Module {
    @Argument(fullName="cortexGraph", shortName="cg", doc="Cortex graph")
    public CortexMap CORTEX_MAP;

    @Argument(fullName="kmer", shortName="k", doc="Kmer to seed the graph reconstruction")
    public String KMER;

    @Argument(fullName="color", shortName="c", doc="Color to process")
    public Integer COLOR = 0;

    @Argument(fullName="maxForksLeft", shortName="mfl", doc="Maximum forks to follow to the left")
    public Integer MAX_FORKS_LEFT = 1;

    @Argument(fullName="maxForksRight", shortName="mfr", doc="Maximum forks to follow to the right")
    public Integer MAX_FORKS_RIGHT = 1;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        if (!CORTEX_MAP.getGraph().hasColor(COLOR)) {
            throw new RuntimeException("Graph '" + CORTEX_MAP.getGraph().getCortexFile().getAbsolutePath() + "' does not contain color '" + COLOR + "'");
        }

        CortexGraphWalker cgw = new CortexGraphWalker(CORTEX_MAP);

        CortexKmer ck = new CortexKmer(KMER);

        DirectedGraph<CortexKmer, CortexEdge> g = cgw.buildLocalGraph(COLOR, ck, MAX_FORKS_LEFT, MAX_FORKS_RIGHT);

        DOTExporter<CortexKmer, CortexEdge> exporter = new DOTExporter<CortexKmer, CortexEdge>(new CortexKmerIDProvider(), new CortexKmerIDProvider(), null);
        exporter.export(new PrintWriter(out), g);
    }
}
