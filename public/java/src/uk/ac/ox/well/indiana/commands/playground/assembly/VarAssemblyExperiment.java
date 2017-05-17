package uk.ac.ox.well.indiana.commands.playground.assembly;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.util.StringUtil;
import org.jgrapht.DirectedGraph;
import org.jgrapht.Graphs;
import org.jgrapht.graph.DefaultDirectedGraph;
import org.jgrapht.traverse.DepthFirstIterator;
import org.jgrapht.traverse.TopologicalOrderIterator;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexKmer;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexRecord;
import uk.ac.ox.well.indiana.utils.io.cortex.links.CortexLinksMap;
import uk.ac.ox.well.indiana.utils.io.gff.GFF3;
import uk.ac.ox.well.indiana.utils.io.gff.GFF3Record;
import uk.ac.ox.well.indiana.utils.sequence.CortexUtils;
import uk.ac.ox.well.indiana.utils.stoppingconditions.NahrStopper;
import uk.ac.ox.well.indiana.utils.stoppingconditions.ShoreStopper;
import uk.ac.ox.well.indiana.utils.traversal.*;

import java.io.File;
import java.util.*;

/**
 * Created by kiran on 08/05/2017.
 */
public class VarAssemblyExperiment extends Module {
    @Argument(fullName="graph", shortName="g", doc="Graph")
    public CortexGraph GRAPH;

    @Argument(fullName="parents", shortName="p", doc="Parents")
    public ArrayList<String> PARENTS;

    @Argument(fullName="child", shortName="c", doc="Child")
    public String CHILD;

    @Argument(fullName="roi", shortName="r", doc="ROI")
    public CortexGraph ROI;

    @Output
    public File out;

    @Override
    public void execute() {
        CortexLinksMap lm = new CortexLinksMap("PG0063-C.ERR019060.k47.3D7.ref.links.raw.ctp.gz");

        int childColor = GRAPH.getColorForSampleName(CHILD);
        Set<Integer> parentColors = new HashSet<>(GRAPH.getColorsForSampleNames(PARENTS));

        Map<String, String> varSeqs = new LinkedHashMap<>();
        varSeqs.put("PF3D7_0100100", "ATGGTGACGCAAAGTAGTGGTGGGGGTGCTGCTGGTAGTAGTGGTGAGGAAGATGCCAAACATGTATTGGATGAATTTGGGCAACAAGTGTACAATGAAAAAGTGGAAAAGTATGCTAATTCTAAAATATATAAAGAGGCGTTGAAAGGAGATTTGTCACAAGCATCAATTTTGAGCGAATTAGCTGGCACCTATAAACCATGTGCCCTTGAATATGAATATTATAAGCATACTAATGGCGGTGGTAAGGGTAAAAGGTATCCGTGTACAGAGTTAGGTGAAAAAGTAGAACCACGTTTTTCGGATACACTTGGTGGTCAGTGTACTAACAAAAAAATAGAAGGTAATAAATATATTAAAGGTAAGGATGTTGGTGCTTGTGCACCATACCGACGTCTACATCTATGTAGTCATAATTTGGAAAGTATACAAACAAATAATTATAATAGTGGTAATGCTAAACATAATTTATTGGTAGATGTGTGTATGGCAGCCAAATACGAAGGGGACTCAATAAAAAACTATTATCCAAAGTATCAAAGAACATATCCTGATACTAATTCTCAATTATGTACCGTGTTGGCACGAAGTTTTGCTGATATAGGAGATATCGTACGCGGTAAAGATCTGTATCTCGGTAATCCACAAGAAAGTACACAAAGAATAATATTAGAAAATAATTTGAAAGATATTTTCGCGAAAATACATAGTGACGTGATGTCAACGAGCGGGAGTAATGGGAGGGCGCTACAAAAACGCTACAAAGATACTGATAATTATTATGAATTGAGAGAAGATTGGTGGGCACTTAATAGAGACCAAGTATGGAAAGCTATCACATGCAATGCTGGGGGTGGTAATAGATATTTTCGACAAACATGTGGTTCAGGAGAATGGGCTAAAGACAAATGCCGGTGTAAGGACGACAAGGTCCCCACATATTTTGACTATGTGCCACAGTATCTTCGCTGGTTCGAGGAATGGGCCGAAGATTTTTGTAGATTAAGGAAACATAAATTAAAAGATGCTAAAAACAAATGTCGTGGAGATAGTGGTAACGATAGATATTGTGATCTTAATAGGTATGATTGCACACAAACTATTAGAGGAAATGAACATTTTGTTGAAAAGGATGATTGTAAAGGTTGTCAGTATTCGTGCGCTCATTTTGTGAACTGGATAGATAACCAAAAACTAGAATTTGAAAAACAAAAAGAAAAATATACAAAAGAAATTAAAAAAAAGCATCCAACAACCATAATAATAAAAACTGCAAATCGAAAAACAACTATTAATAACTTATATGTAAAAGAATTTTATAAAAAACTTCAAGAGAAATATGGAGATGTCGAAAATTTTTTACAAAAATTAAATGAAGAACAAATATGCAAAAATCAACCGTACAATGATGAAAGTAGTATTGATATTAATTTCAAAAGTATTAAAGATATTGACATATTTTCTCATACGGAATACTGTCAAGCATGTCCATGGTGCGGGGCCAAACGTAAAGGTAAAGGATGGGAACCTAAAGAAAAAACCTGTGGAAAAACAAAGACATACGATCCTAAGAAAACAACGAATATACCAATACTTACCCCTTATATATCACAGCAAAGTATACTAAAAAAATATAATAAATTTTGTAATGGTAATGGTGGAAATGGTGCACCTGCTACTGCAACTGGTGGTGGTCAAATTAAAAATTGGCAATGTCATTATGAAGGTGATAATAATGATAATTGTGTAGAAGGAGAATGGAAAGAGTTTAAAGAGGGTAAAAACGTTATGTCCTATAATGCTTTTTTTTGGAAGTGGGTTCATGATATGTTAATCGACTCTATGCAATGGAGAAATGAACATGGAAATTGTATAAATAAAGATAATGACAACACATGTAAAAATTCATGCAAAAGACCATGTGAATGTTTTAAAAGATGGGTAGATCAAAAAAAAAAAAACGAATGGGAGGCAATAAAAGACCATTTTAAAAAGCAAAATATTGCAGCTGAAACACAATGTGATCCTGGCGTAACTCTTCAATGGGTTTTGATATTAGACTTTTTGAAAGACGAATCCACAGAAGATAAGGAAAATAAGGTGAGTGCAGAGGAGGCAAAGGAAATAAAACACCTTCGCCAAATGTTGCAACAAGCAGGCGTTGATGATCCTGCTGCTTTTGCTCGTCCGTGTACTGAAGATGGTGTCGCTGAACAGGACACTATAATGGATAAATTGCTCAATCGCGAAGAAAACGATGCCACTGAATGCAAAAAATGCGACAAACCACCACCAGCACCCACTGCAGGAGATCGTGGCCCTGGAGCCCGCGCCGACCCCCACGACGTCCAACAGCCACGACCTCCTGGTAGTGGCCCCGGCACGGACGCCAACGACGAAGACGATGATGACGATGATGACGATGATGACGAAGAAGACGGTGAAGCCAAAGAAGAAGAAGAAGACGAGGAAAAACAAGAGGACGTCCACCAGGAGGAAAAGGCAAAGAAGGAAGAACCACAAAAAGAGGAGGTGGCACGAACACCAAAAGACGATGTAAATGTGTGCAATATAGTGAACAATGTGTTTACAGACGGCAGTAGTCTCCAAGCAGCGTGCTCTCTCAAATATGGCAAAAACGCACCCACAAGTTGGAAGTGTGTCACACCAAGTGGTAACACGAGTGACACCACTGTCAAAAGTGGTGACACCACCGGTGGTAGTATTTGTGTGCCACCCAGGAGACGACGATTATATGTCACACCACTAACGAGATTGACAGGTGGTGACAGTACCACACAGGCGTCACAGGCGAGTGAGGTACAGACACAAGCACGTGGTAGTAACACGGATAAGTCACCAGGTAGTAGTGAGGCAGCACAAGGTGACGGCGTGTCGAAAGACCCACAAAAGGCACTACTCAAAGCTTTTGTTGAGTCTGCAGCAGTTGAAACCTTCTTCCTATGGGATAGATATAAAAAAATAAAAGAGAAGGAGAAAAAGGAAAAAAAGAAAACATATGAACAAATATATGAATCAACCGACTATGACGATGAAGAAAAAGATCCACAAGAAGAATTAAAAAAAGGAATAATCCCTGATGAGTTTAAGCGTCAAATGTTTTATACGTTAGGTGACTACAAAGATATATTATACAGTGGTGATACGGTGAATGGTGGTAATGAGGACAAAATAAAAAAAGCTATAAATAACTATTTTCAAAAAATTCGTGAACAATCTTCTAGTGATAACAACCCATCTCCTCGTAGTGTCAAAACCCCTTCAACTAGTGACAAGGACCCTCAAACCTGGTGGAATGCACACGCCCCTTCCATCTGGAATGCTATGGTATGTGCTTTAACATATGATACAAACAGTGGCGGAGAGGGCAAAACCACAACTATTACGCAGGATCCTAATTTGAAAACTGCACTTTGGGACGAAAACGGCAAAAAACCCCTCAAAACCAAATACCAATATGATAGTGTCACAATTGGTGCTAGTGGTGCCAAACCCCAAACCAAAGCCAAACCCACTGGTGGTGACACCCCCCTCACCCAATTTGTGTTACGCCCCACCTACTTCCGATACCTTGAAGAATGGGGTCAAAATTTTTGTAAAAAACGAACAGAGATGTTGGAGAAAATAAAATATGAGTGTAAAGTAGGACAAGGTCGTGGTGGTCGTAAACAAAAAACCCCACAATGTAGTTGTTATGGGGAAAATTGTGACGATCAGCTTGACGACAATCCTAGTACTGATGCGGATTTAAAATGTCCTGGTTGTGGAAGAGAATGTAGAAAATATAAAAAATGGATAGAAAAAAAAAAAGAGGAATTTACTAAACAATCAAATGTATATGAAGAACAAAAAACAAAATGCCAAAAGGAAAGTAAGAGTGCTAAAGGTAATAACCACGGTAATGAATTTTGTGGAACACAAGGAACGTGCGATACAGCTGGAGACTTTTTAAATAGGTTAAAAAGTGGACCATGTAAAAAGGAGAATGGAAAGGATAATCAAGAGGATGAAATAAATTTTAAGGATGAAGATAAAACATTTGGACATGAAAATTATTGTGCTCCATGTCCTGTATTTAAAGATATATGTAAAAAAAAGGATTGCCGTAATGCTTCTAACAATATGTGCAATGGAAAAGATTTTATTACTGCAGAAGATATTAAAATAATGGACAGCAGTAGTGAAGAAGTTAATATGCTTGTGAGTGATAACGATACAAATAAATTTGATGGTGGTTTAGACGCTTGTAAAGATGCCCATATCTTTAAAGGTATTAAAGAAAATAAATGGTCATGTGGTAACGTATGTGGTTATAATGTGTGTAAACCGAAAAAAGTTAATGGGGAAAAAGGTAGTGGGGAAAACAATGATCAAATTATAACAATTAGAGGTTTGGTTACACATTGGGTACAAAATTTTTTAGACGATTATAATAAAATTAGAACAAAATTAAAGCCATGTAGGAATAATGGTGAGGTATCCAAATGTATAAAAGATTGTGTGAAAAAATGGGTAGAAAAAAAAACTGAAGAATGGCCAAAAATACGAGATCGTTACTTGGAACCATATAAAAGTGATGATGGCTATAACAAAAAATCTTTGGTTAGAAGTTTTATGGAGACCTTGATACCTCTAATGGATCTTACAAATGGTAAGGAAAAGATTCAAGAATTAAATAAGTTCCTTAGGTCATATGAATGTAATTGCGCTGATAACTCACAACAAAAAGGTGATACACCAAAAGACATCGTAGAATGTTTGCTTGAAAAGCTTGAAGATAAAGCAAACAAGTGTAAAACCCAAACTAGTGGTACCGACTGTCACCCCTCCACCCCCCTTGAAGATGACGATGAACCCCTTGAAGAAACAGAAGAAAATACTGTGGAACAACCGAACATTTGTCCAACAAAACAACCACAACCAGAGAAAGAAGACGGTTGTGAAGCAGCACCAACAACAGCAGAAGAAACGTCACCAACAGCAACTAGTGAAGGCACAGAGAACCAATCCCCTCCACCTCCTCCTCCAGCACCAGCACCAGCACCAGCACCGGCACCAGAAAAATCACAACCAAAAGAAGACAAAAAAGTGGAACCACAACCCAAACCACAACCAACAAACCCCCCCCCAAATTTGTTCAACAACCCCGCTGTTATACCCGCCCTCATGTCTTCTACCATCATGTGGAGTATTGGCATCGGTTTTGCTGCATTCACTTATTTTCTTCTAAAGgtattatatatatatatgtatatgtggggatgtgttttttttatatgtatttgtggggtgtgtttggatatatatatatgtatatgtgtttctgtatatgtgttttctgtatatgtatgtgttcgtatgtttggatatatatttgtgtatatgtatgtgttttatatatattttatatatatgtatttatattgataaagaaaaaaatgaaaaaaagaaaaaaaaaaatttattaaaataaaaaaaaaaaaaaaaaaaaggagaaaaatattttaaaaataataaaaattaaaataaaaatataaattttgataaaataaaaaatgaaaaatattatcaaaaagaaattaaaaaaaattttatatataaaaaaaatgattataaaaaaaatttattagaaataaaataaaaaaaaatttattaaataaaaacaaaaaaaaaaaaaaaggagaaaaatattttaaaaataataaaaattaaaataaaaatataaattttgataaaatgaaaaatgaaaaatattatcaaaaagaaattaaaaaaaattttatatataaaaaaaatgattataaaaaaaatttattagaaataaaataaaaacaaatttattaaataaaaaaaaaaaatgttaaaaaaaaatatatatatcataaaataaaaaaaaaagaaaaaaatatattaaaaaaaaaatatatatcataaaataaaaagaaaagggaaaaaatgtttaaaaaaaaaataaaaataaaaaaaaaaaaaaaaaaaaaaaaaaaattaaaaaaaaaaaataaaattaaaataaaaaaaaataaaaaaaatttaattaaataaaaaaaaaaaaatttaattaaataaaaaaaaaaaaattaaattaaatacatgcacatatacatatacgcatacatatacatatacacataaatatatatattatatatatatacccataactacattcacatatacacatacatatatatattatatatatatatacccataactacatacatatatacattaacaaacacatagatatacataaatacatatatacattaacaaacacatatatatacctaaatacatatatacatacacatatatgttcatttttttttttagAAAAAAACCAAATCATCTGTTGGAAATTTATTCCAAATACTGCAAATACCCCAAAACGATTATGGAATACCAACATTGAAATCCAAAAATAGGTACATACCATATAGAAGTGGTACATATAAAGGCAAAACATATATTTATATGGAAGGAGATAGCAGTGGAGATGAAAAATATGCATTTATGTCTGATACTACTGATGTAACTTCCTCAGAAAGTGAGTATGAAGAATTGGATATTAATGATATATATGTACCGCATGCTCCTAAATATAAAACATTAATTGAAGTAGTACTTGAACCTAGTGGTAACAACACAACAGCTAGTGGTAAAAACACACCTAGTGATACACAAAATGATATACCAACTAGTGATACACCACCACCCATTACTGATAATGAGTGGAATACATTGAAAGATGAATTTATATCACAATATCTACAAAGTGAACAACCAAAGGATGTACCAAATGATTATAAAAGTGGTGATATTCCATTGAATACACAACCGAATACTTTATATTTTAATAAACCTGAAGAAAAACCTTTTATTACTTCTATTCATGATAGGGATTTATATACTGGAGAACAAATTAGTTATAATATTCATATGAGTACTAATACTATGGATGATCCAAAATATGTATCAAATAATGTATATTCTGGTATAGATTTAATTAATGACGCACTAAATGGTGATTATGACATTTACGATGAAATATTGAAACGAAAAGAAAATGAATTATTTGGAACAAATCATGTGAAACAAACAAGTATACATAGTGTTGCCAAACTAACAAATAGTGACCCCATCCACAACCAACTGGAACTATTCCATAAATGGTTAGATAGACATAGAAATATGTGTGAAAAGTGGAAAAATGATAATGAGCGGTTAGCCAAATTAAAAGAAGAGTGGGAAAATGAGACACATAGTGGTAACACTCACCCTAGTGATAGTAACAAAACGTTAAATACTGATGTTTCTATACAGATAGATATGGATCATGAAAAACGAATGAAGGAATTTACTAATATGGATACTATCTTGGATGATTTGAAAACATATAATGAACCTTATTATGATGTGCAAGATGATATTTATTATGATGTAAATGATCATGATGCATCAACTGTGGATAGTAATAATATGGATGTTCCCAGTAGAGTACAAATTGAAATGGATGTAAATACGAAATTGGTGAAAGAGAAATATCCTATAGCCGATGTATGGGATATATAA".toUpperCase());
        varSeqs.put("PF3D7_0223500", "ATGGGGAGTGGTAAGGGCGGTGATCCGCAGGATGAAAGTGTCAAACATATGTTTGATAGGATAGGAGAAGATGTGTACGAGCAAGTGAAAAGTGAAACTGTAAATTATGTTAGTGAATTGGAAGGAAAGTTGTCACTAGCACCAATTTTGGGTGTGGAATCAGGTAGCACCAATGAAACATGCAACCTTGTACAGGATTATTATAATAAGCCTGTTTATGGTAACAGTAACAGGTATCCGTGCAAAAATTTAAAAGGAATTACAAATGAAGAACGTTTTTCGGATACACTTGGTGGCCAGTGTACTAACAAAAAAATAAAAGGTAATGAATATAGTACTAAAAGTGGTAAAGATTGTGGAGCATGTGCACCATACCGACGTCTACATTTATGTAGTCATAATTTGGAATCTATAGACACAACGTCGATGACGCATAAGTTGTTGTTAGAGGTGTGTATGGCAGCAAAATACGAAGGAAACTCAATAGATACACATTATCCACAACATCAACGAACTAATGAGGATTCTCCTTCTCAAATATGTACTATGTTGGCACGAAGTTTTGCAGATATAGGTGATATTGTAAGAGGAAAAGATTTATTTTATGGTAATAGCAAAGAAAAAGAAAAAAGAGATGAATTAGAAACCAATTTGAAAACAATTTTCGGGAAAATACATGAAAAATTGAAGGATAAGGAAGGAGCAGAAACTCGTTACGGAAGTGATACTACAAATTATTATCAATTACGAGAAGACTGGTGGTATGCGAATCGCGCCACAGTGTGGGAAGCTATCACGTGCGACGTTCATGGTTCTGACTATTTTCGACAAACATGTGGTGATAAAGAAACCACTGCAACTCGGGTTAAAGACAAATGCCGCTGTAAGGACGAAAACGGCAAAAAGCCCGGCTCAAATGCCGACCAAGTCCCCACATATTTTGACTACGTGCCGCAGTATCTTCGCTGGTTCGAGGAATGGGCAGAAGACTTTTGTAGGAAAAAAAAAAAGAAATTAGAAAAGTTGGAACAACAGTGTCGCGATTACAAACAAAATTTATATTGTAGTGGTAATGGCTACGATTGCACAAAAACTATATACAAAAAAGGTAAACTTGTTATAGGTGAACATTGTACAAACTGTTCTGTTTGGTGTCGTCTGTATGAATCTTGGATAGATAACCAAAAACTAGAATTTCTAAAACAAAAACAAAAATACGAAACAGAAATATCAAATAGCGGTAGTTGTGGTGGGAGTGGTGGTGTTAAGGGTAGGAATAGGAAAAAACGGGGTGCAGGTGTAGAAACTGCTACTAATTATGATGGGTATGAAAAAAAATTTTATAAAGAACTGAAAGAAAGTGAGTATGGAAAAGTCGATGATTTTTTAAAATTATTAAATAATGAAGATGTATGCAAAAAAATTAAGGATGAAAAAGAAAAAATTGATTTTACCAAACCTGCTGATAAAAATAGTAATAATGAAGGAACATTTTATCATTCGGAATATTGTAAACCGTGTCCCGACTGTGGGGTCAAACGTAAAGATAATCAATGGAAAGATAAATATGATGGCAAGTGCACACGTGGAAAACTTTATGAGCCTGCAAGTGGCGCACAAGGTACTCCTATTAAAATCCTTAAAAGTGGTGAAAAACAAAAAGAAATTGAAACAAAATTAAAAGCGTTTTGCGATCAAACAAATGGTGATACAACAAATAGTGTTGCTAGAGGCGGTGGCGCTGATGGTAGTGGTAGTAAGAGTAATAGTAAGGAACTGTATGAAGAATGGAAATGTTATAACGAGGTACAGAAAGTTAAAGATGATAAAAATGGAGAGGAAGAGGATGAAGACGAGGAAGATGTAGACAAGGTAAAAAAAGCAGGCGGATTATGTATATTGGAAAACAAAAAACATGAAAGTAGAAATAATTCTTCAAATGAACCTGAGCAATTCCAAAAGACATTCCATGATTTTTTTTACTTTTGGATAGGACGTTTTTTGAACGATTCTATGTATTGGAGAGGAAAAGTTAACAGTTGTATAAATAATCCTAAGCGAAAGAAATGTAGAAATGAATGTAAGGATGATTGTGGTTGTTTTAAAGAATGGATTGGAAAAAAGAAAGAAGAATGGGAAAATATAAAAAAACATTTTAAAACGCAAGAAGCTTTTAAGAATAAACGAGAAAATAGCGGAATTGACATGTTCAGCGGACTAATGGATTCTGCTGATGTTGTTCTTGAATTGGCTTTGGAATTAGAACAACTTTTCCAAGATATTAAAGATGGTTATGGGGATGTAAAGGAATTAAAAGGAATTAAAGAACTGTTGGATGAGGAAAAAAAAAAAAAACAAGCAGAAGAAGCAGTTGTTGTTGTTGTTGCCGACAATCAAAAGAAGACCACAATTGATAAATTACTACAACATGAAGGAGACGATGCCAATAACTGCCTAAAAACACACAAAGAAAAATGCGAAGAAACGCAACCAAAACCACCCGGCGCTGGAGGTCCTGGTGCCCCCTCCGAAACCGGAGAAACCACTACACTTGAGGACGAAGAAGAAGAAGAAGACGAAGAAGAAGACGCAGGCGACGAAGTCGAGGAGGGGGAGACGGTGGACACCACAGAAGGGGATGAGACAGAGACGGTGGAGCAGCCGGTGAAGGACACGGACAGGGAGGGGGAGGAGGAAGAGGCAAAGAAGGCAACAGATACGACTACATCACTAGACGTTTGCGACACAGTGAAAAACGCACTCACAAACAACGACAATCTCACTGATGCATGTAAACTAAAATACGGTCCAGGTGGAAAGGAAAGATTCCCCAATTGGAAATGTGTATCAAGTGGTGAAAAAAGTGTTGCCACTGCCGGTAGTAGTGGTGCCACTGGCAAAAGTGGTGATAAGGGTGCCATTTGTGTGCCACCCAGGAGGCGACGACTATACGTGGGTGGGTTAACCAAGTTGACAAGTGCTGGCACGTCTAGTGAGTCACCACAGGGGGGTAGTGAGTCATCACGGGCGAGTGATGTGTCACAAGGTAACGGCGGCGACGACATCACCACCACCGAGTCATTACGTAAGTGGTTTATAGAGACGGCAGCTATAGAGACTTTTTTCTTATGGCATAGATATAAAAAAGAGTGGGAGGCACAAAAGAAGGCGGAACTACAACGAAATGGATTACTACTCGGCACAGGTGCTAGCCTCAACCTTGGTGGTGATGACTCCAACCCCCAAACACAATTACAAAAAAGTGGTACCATACCCCTCGATTTCTTGAGATTAATGTTTTATACTTTAGGTGATTATAGAGATATTTTGGTACGAGGTGTTGCTGACGACAAAAACGGTGGCAACAACATAATACTTAATGCGAGTGGTAACAAGGATGAAAAACAGAAAATGGAGAAAATACAAGAGAAAATAGAACAAATTCTTCCAACTAGTGGTAACAAAGAAACTCGTGGCCCCCAAAATAGTGTCAATGACCGTCAATCCTTGTGGGATAGAATCGCCGAACATGTTTGGCATGGAATGGTTTGCGCATTAACATATAAAGATGACGACAATGGCCTCAAAGGCGTCGTAAAAAAACCACAAAAGATTGAAAATCCGGAGAAACTTTGGAACGAAACAACCAAAAAACCCAAAGACGAGAAATACCAATACCAAACTGCCAAACTCGAAGATGAAAGTGGCGAAAAACGACCAGACTCCTCAGCCAGTGGTACGAAATTAACCGACTTCATCAAACGCCCCCCTTATTTCCGTTACCTTGAAGAATGGGGTGAAAATTTTTGTAAAAAACGAACAGAGATGTTGGGGAAGATAAAGGAGGATTGCTACAAAAATGGTGGACGTTGTAGTGGTGATGGTTTGAAATGTAACGAAATAGTTATAGATAAGGAAAAAATTTTTGGCGATTTACTTTGTCCGACGTGTGCCAGACATTGTAGATTTTATAAAAAGTGGATAAACACAAAAAGGGACGAATTTAATAAACAATCAAATGCATATTCTGAACAAAAAAAAAAATACGAAGAGGAAAATGATAGTGCTCAAAAGAATAATGGAGTTTGCGGAACACTAAAAGATGACGCTGCAGAATTTTTAAATAGGTTAAAAAACGGACCATGTAAAAATGAGAGTGAAGAGAATAAAAAAGCAGAGGATGAAATAGATTTTAAGAAACCAGATGATACATTTAAAGATGCAGATAATTGTAAACCATGTTCTGAATTTAAAATTAAATGTGAAAATCATAATTGCAGCAGTGGTGGTAATACACAAGGGAAGTGCGATGGAAAAACGACTATTGCTGCAACAGAAATTGAAAATATAAAAACAAATACTAAAGAAGTTACTATGCTTGTGAGTGATGACAGTAAAAGTGCAACGGAATTTAAGGATGGTTTAAGCGAATGTAAAGATAAAGGTATATTTAAAGGTATTAGAAAAGATGAATGGGAATGTGGCAAAGTATGTGGTGTAGATATATGTAATCTGAAAAAAAAAGATAACATTGGGAAAGAAAGCGATAAAAAATATATCATAATGAAAGAATTGCTTAAACGATGGTTAGAATATTTTTTAGAAGATTATAATAAAATTAAACATAAAATTTCACATTGTACGAAAAATGGTAAAGGATCCAAATGTATAAAAGGTTGCGTAGATAAATGGGTACAACAGAAAAAGGAAGAATGGAAACAAATAAAAGAACGTTTCAATGAACAATATAAAAGTAAAACCTCAGATGAATATTTTAACGTTAAAAGTTTTTTGGAGACCTGGATACCTAAAATTGCTGTTGTAAATGATCAAGATAATGTTATAAAATTAAGTAAGTTCGGTAATTCTTGTGGATGTAGTGCCAGTGCGATCTCAACAAATGGTAATGAGGAGGATGCTATAGATTGTATGATTAAAAAGCTTGAAAAAAAAATTGACGAATGCAAAAGGAAACCTGGCGAAAATAGTGGTCAAACATGTAACGAAACACTAACACATCCCCTTGACGTTCAGGATGAAGATGAACCCCTTGAAGAAACAGAAGAAAACCCAGTGGGAAAACAACACCCATCATTTTGTCCGCCAGTGGAAGATAAAAAAAAAGAGGAAGAAGGAGAAACTTGTACACCGGCATCACCAGCACCAGCACCAGCACCAGCACCAGCATCTCCATCCCCGACACCGGCCCCTGCGGATGAACCGTTTGACCCAACTATACTACAAACAACCATTCCTTTAGGTATTGCGCTGGCATTAGGATCCATTGCTTTTTTATTTTTGAAGgtaatatatatatgtgtggtatatatgtatatatatatgtgtttctgtatatatatgtatgtgtgggtgtgtttggatatatatatatgtgtatgtataagtgtttgtgtatatgtatgtgatttatatatattttatatatatgtatttatattgaaaaagaaaaaaaaaaaaaaaaaaaaaaaaaaaaatttattaaaataaaaaaaaaaaaaaaaaaaaagagaaagattttaaaaataataaaaattataataaaaatataaattttgatagaataaaaaatgaaaaatattatcaaaaaaaaattaaaaaaaattttatatataaaaaaaatttattagaaataaaataaaaacaaaagaagaaaaaaaaaacattaaaaaaaaaaaaaatatatatatcataaaataaaaaaaaattaaaaaaatgttaaaaaaaaaatatatatcataaaataaaaaaaaaattaaaaaaatgttaaaaaaaaatatatatatcataaaataaaaaaaaaattaaaaaatttaattaaataaaaaaaaataataaataaaaaaatttaattaaataaaaaaaaaaaattaaaaaaaataaaataaaaaaaaaaaataaaaaaattaaaaaaaaaaaaaaaaaaaaaaatattttattcatacacatacatatacacatatatatatacatatattatatacatacacatatacctacatacatatacaaacctacttatacatacatacctcttttattattttagAAAAAAACTAAACACCCTGTCGACCTTTTCAGTGTTATTAATATCCCCAAAAGTGATTATGATATACCGACAAAACTTTCACCCAATAGATATATACCTTATACTAGTGGTAAATACAGAGGCAAACGGTACATTTACCTTGAAGGAGATAGTGGAACTGATAGTGGTTACACCGATCATTATAGTGATATTACTTCATCTTCCGAAAGTGAGTATGAAGAAATGGATATTAATGATATATATGTACCTGGTAGTCCTAAATATAAAACATTGATAGAAGTAGTACTTGAACCTAGTGGTAACAACACAACAGCTAGTGATACACAAAATGATATACAAAATGATGGTATACCTAGCAATAAATTTAGTGATAATGAATGGAATACATTGAAAGATGATTTTATATCTAATATGTTACAAAATCAACCAAAGGATGTACCAAATGATTATAAAAGTGGAGATATTCCATTCAATACACAACCGAATACTTTATATTTTGATAAACCTGAAGAAAAACCTTTTATTACTTCTATTCATGATAGAAATTTACTTAACGGAGAAGAATATAGTTATAATGTTAATATGAGTACTAATAGTATGGATGATCCAAAATATGTATCAAATAATGTATATTCTGGTATAGATTTAATTAATGATTCACTAAGTGGTAACAAACATATTGATATATATGATGAAGTTTTGAAACGAAAAGAAAATGAATTATTTGGAACAAATCATGTGAAACATACGAGTATACATAGTGTTGCAAAAAATACAAACAGTGATCCTATACTCAATCAAATAAATTTGTTCCATACATGGTTAGATAGACATAGAGATATGTGCGAAAAGTGGGAAAATCATCACGAACGATTAGCCAAATTGAAAGAAGAGTGGGAAAATGAGACACATAGTGGTAACACTCACCCTAGTGATAGTAACAAAACGTTAAATACTGATGTTTCTATACAAATACATATGGATAATCCTAAACCTATAAATCAATTTACTAATATGGATACTATCTTGGAGGATCTGGACAAACCATTTAATGAACCCTACTATTATGATATGTATGACGATGATATTTATTATGATGTAAATGATCATGATACATCAACTGTGGATACTAATGCTATGGATGTACCTAGTAAAGTACAAATTGAAATGGATGTAAATACCAAATTGGTGAAAGAGAAATATCCTATAGCAGATGTATGGGATATATAA".toUpperCase());

        String sk = varSeqs.get("PF3D7_0100100").substring(0, GRAPH.getKmerSize());

        TraversalEngine e = new TraversalEngineFactory()
                .traversalDirection(TraversalEngineConfiguration.TraversalDirection.FORWARD)
                .combinationOperator(TraversalEngineConfiguration.GraphCombinationOperator.OR)
                .traversalColor(childColor)
                .recruitmentColors(parentColors)
                .stopper(NahrStopper.class)
                .graph(GRAPH)
                .rois(ROI)
                .links(lm)
                .make();

        DirectedGraph<CortexVertex, CortexEdge> dfs1 = e.dfs(sk);
        String contig1 = e.getContig(dfs1, sk, childColor);

        String lk2 = contig1.substring(contig1.length() - GRAPH.getKmerSize(), contig1.length());
        String nk2 = CortexUtils.getNextKmer(GRAPH, lk2, 0, false);
        DirectedGraph<CortexVertex, CortexEdge> dfs2 = e.dfs(nk2);
        String contig2 = e.getContig(dfs2, nk2, childColor);

        String lk3 = contig1.substring(contig2.length() - GRAPH.getKmerSize(), contig2.length());
        String nk3 = CortexUtils.getNextKmer(GRAPH, lk3, 0, false);
        DirectedGraph<CortexVertex, CortexEdge> dfs3 = e.dfs(nk3);
        String contig3 = e.getContig(dfs3, nk3, childColor);

        StringBuilder sb = new StringBuilder(StringUtil.repeatCharNTimes(' ', GRAPH.getKmerSize() - 1));
        StringBuilder rb = new StringBuilder(contig1.substring(0, GRAPH.getKmerSize() - 1));
        StringBuilder cb = new StringBuilder(contig1.substring(0, GRAPH.getKmerSize() - 1));
        List<Integer> po = new ArrayList<>();

        int i0Last = 0;
        int i1Last = 0;

        for (String contig : Arrays.asList(contig1, contig2, contig3)) {
            for (int i = 0; i <= contig.length() - GRAPH.getKmerSize(); i++) {
                String sk0 = contig.substring(i, i + GRAPH.getKmerSize());
                CortexKmer ck0 = new CortexKmer(sk0);

                int i0 = varSeqs.get("PF3D7_0100100").indexOf(sk0, i0Last);
                int i1 = varSeqs.get("PF3D7_0223500").indexOf(sk0, i1Last);
                boolean isNovel = ROI.findRecord(ck0) != null;

                if (isNovel) {
                    sb.append(".");
                    rb.append(" ");
                    po.add(0);
                } else if (i0 >= 0 && i1 == -1) {
                    sb.append("0");
                    rb.append(varSeqs.get("PF3D7_0100100").substring(i0 + GRAPH.getKmerSize() - 1, i0 + GRAPH.getKmerSize()));
                    po.add(i0);
                    i0Last = i0;
                } else if (i1 >= 0 && i0 == -1) {
                    sb.append("1");
                    rb.append(varSeqs.get("PF3D7_0223500").substring(i1 + GRAPH.getKmerSize() - 1, i1 + GRAPH.getKmerSize()).toLowerCase());
                    po.add(i1);
                    i1Last = i1;
                } else {
                    sb.append(" ");
                    rb.append("_");
                    po.add(0);
                }

                cb.append(sk0.substring(sk0.length() - 1, sk0.length()));
            }

            sb.append("|");
            rb.append("|");
            cb.append("|");

            log.info("");
        }

        String llk = cb.subSequence(cb.length() - GRAPH.getKmerSize() - 1, cb.length() - 1).toString();

        log.info("{}", rb);
        log.info("{}", cb);
        log.info("{}", sb);
        log.info("{}", cb.length());
        log.info("{} {}", llk, GRAPH.findRecord(new CortexKmer(llk)));
        log.info("");

        int pos[] = { 0, 0 };
        int pi = 0;
        String[] ids = { "PF3D7_0100100", "PF3D7_0223500" };
        boolean keepGoing = true;
        while (keepGoing) {
            String qk = varSeqs.get(ids[pi]).substring(pos[pi], pos[pi] + GRAPH.getKmerSize());
            CortexRecord qr = GRAPH.findRecord(new CortexKmer(qk));

            Set<CortexVertex> vs = e.getNextVertices(qk);

            log.info("qr: {} {} {} {} {} {} {}", pos, qk, qr.getKmerAsString(), qr.getCoverage(0), qr.getCoverage(childColor), qr.getEdgesAsString(0), qr.getEdgesAsString(childColor));

            if (vs.size() == 1 && vs.iterator().next().getSk().equals(varSeqs.get(ids[pi]).substring(pos[pi] + 1, pos[pi] + 1 + GRAPH.getKmerSize()))) {
                pos[pi]++;
            } else {
                TraversalEngine f = new TraversalEngineFactory()
                        .traversalDirection(TraversalEngineConfiguration.TraversalDirection.FORWARD)
                        .combinationOperator(TraversalEngineConfiguration.GraphCombinationOperator.OR)
                        .traversalColor(childColor)
                        .joiningColors(parentColors)
                        .recruitmentColors(parentColors)
                        .stopper(ShoreStopper.class)
                        .graph(GRAPH)
                        .rois(ROI)
                        .links(lm)
                        .make();

                boolean found = false;
                for (CortexVertex cv : vs) {
                    DirectedGraph<CortexVertex, CortexEdge> d = f.dfs(cv.getSk());
                    if (d != null) {
                        String contig = f.getContig(d, cv.getSk(), childColor);

                        String rk = contig.substring(0, GRAPH.getKmerSize());
                        for (int j = 0; j <= contig.length() - GRAPH.getKmerSize(); j++) {
                            rk = contig.substring(j, j + GRAPH.getKmerSize());
                            CortexRecord rr = GRAPH.findRecord(new CortexKmer(rk));

                            log.info(" -- {} {} {} {} {} {} {}", j, rk, rr.getKmerAsString(), rr.getCoverage(0), rr.getCoverage(childColor), rr.getEdgesAsString(0), rr.getEdgesAsString(childColor));

                            int q0 = varSeqs.get("PF3D7_0100100").indexOf(rk, pos[0]);
                            int q1 = varSeqs.get("PF3D7_0223500").indexOf(rk, pos[1]);

                            if (q0 >= 0 && q1 == -1) {
                                pi = 0;
                                pos[pi] = q0;
                                found = true;
                                break;
                            } else if (q0 == -1 && q1 >= 0) {
                                pi = 1;
                                pos[pi] = q1;
                                found = true;
                                break;
                            } else if (q0 == -1 && q1 == -1) {
                                int a = 1;
                            } else if (q0 >= 0 && q1 >= 0) {
                                int a = 1;
                            }
                        }

                        if (!found) {
                            DirectedGraph<CortexVertex, CortexEdge> d2 = f.dfs(rk);
                            if (d2 != null) {
                                contig = f.getContig(d, cv.getSk(), childColor);

                                for (int j = 0; j <= contig.length() - GRAPH.getKmerSize(); j++) {
                                    rk = contig.substring(j, j + GRAPH.getKmerSize());
                                    CortexRecord rr = GRAPH.findRecord(new CortexKmer(rk));

                                    log.info(" -- {} {} {} {} {} {} {}", j, rk, rr.getKmerAsString(), rr.getCoverage(0), rr.getCoverage(childColor), rr.getEdgesAsString(0), rr.getEdgesAsString(childColor));

                                    int q0 = varSeqs.get("PF3D7_0100100").indexOf(rk, pos[0]);
                                    int q1 = varSeqs.get("PF3D7_0223500").indexOf(rk, pos[1]);

                                    if (q0 >= 0 && q1 == -1) {
                                        pi = 0;
                                        pos[pi] = q0;
                                        found = true;
                                        break;
                                    } else if (q0 == -1 && q1 >= 0) {
                                        pi = 1;
                                        pos[pi] = q1;
                                        found = true;
                                        break;
                                    } else if (q0 == -1 && q1 == -1) {
                                        int a = 1;
                                    } else if (q0 >= 0 && q1 >= 0) {
                                        int a = 1;
                                    }
                                }
                            }
                        }
                    }
                }

                if (!found) {
                    pos[pi]++;
                }
            }

            if (pos[pi] + GRAPH.getKmerSize() + 1 >= varSeqs.get(ids[pi]).length()) {
                keepGoing = false;
            }
        }

        log.info("");

        /*
        Map<String, List<String>> varSeqList = new HashMap<>();
        Map<String, Integer> varSeqIndex = new HashMap<>();
        Map<String, String> varSeqId = new HashMap<>();
        for (String id : varSeqs.keySet()) {
            String seq = varSeqs.get(id);
            varSeqList.put(id, new ArrayList<>());

            for (int i = 0; i <= seq.length() - GRAPH.getKmerSize(); i++) {
                String sk = seq.substring(i, i + GRAPH.getKmerSize());

                varSeqList.get(id).add(sk);

                if (!varSeqIndex.containsKey(sk)) {
                    varSeqIndex.put(sk, i);
                    varSeqId.put(sk, id);
                }
            }
        }

        String id = "PF3D7_0100100";
        int i = 0;
        String lk = varSeqList.get(id).get(0);

        do {
            String sk = varSeqList.get(id).get(i);

            CortexKmer ck = new CortexKmer(sk);
            CortexRecord cr = GRAPH.findRecord(ck);

            if (cr.getCoverage(childColor) > 0) {
                log.info("{} {} {} {} {} {} {} {} {} {} {}", id, i, sk, cr.getCortexKmer(), cr.getCoverage(childColor), cr.getCoverage(0), cr.getCoverage(17), cr.getEdgesAsString(childColor), cr.getEdgesAsString(0), cr.getEdgesAsString(17), id);
                lk = sk;
                i++;
            } else {
                TraversalEngine e = new TraversalEngineFactory()
                        .traversalDirection(TraversalEngineConfiguration.TraversalDirection.FORWARD)
                        .combinationOperator(TraversalEngineConfiguration.GraphCombinationOperator.OR)
                        .traversalColor(childColor)
                        .recruitmentColors(parentColors)
                        .stopper(NahrStopper.class)
                        .graph(GRAPH)
                        .rois(ROI)
                        .links(lm)
                        .make();

                DirectedGraph<CortexVertex, CortexEdge> dfs = e.dfs(lk);

                String contig = e.getContig(dfs, lk, childColor);

                for (int j = 0; j <= contig.length() - GRAPH.getKmerSize(); j++) {
                    String skn = contig.substring(j, j + GRAPH.getKmerSize());
                    CortexKmer ckn = new CortexKmer(skn);
                    CortexRecord crn = GRAPH.findRecord(ckn);

                    log.info("{} {} {} {} {} {} {} {} {} {} {}", "recon", j, skn, crn.getCortexKmer(), crn.getCoverage(childColor), crn.getCoverage(0), crn.getCoverage(17), crn.getEdgesAsString(childColor), crn.getEdgesAsString(0), crn.getEdgesAsString(17), varSeqId.get(skn));

                    lk = skn;
                }

                if (varSeqIndex.containsKey(lk)) {
                    i = varSeqIndex.get(lk);
                    id = varSeqId.get(lk);
                } else {
                    i = -1;
                }

                log.info("---");
            }
        } while (i >= 0 && i < varSeqs.get(id).length());

        log.info("");

        for (String id : varSeqs.keySet()) {
            String seq = varSeqs.get(id);

            Set<String> visited = new HashSet<>();

            String lk = seq.substring(0, GRAPH.getKmerSize());
            for (int i = 0; i <= seq.length() - GRAPH.getKmerSize(); i++) {
                String sk = seq.substring(i, i + GRAPH.getKmerSize());

                if (!visited.contains(sk)) {
                    CortexKmer ck = new CortexKmer(sk);
                    CortexRecord cr = GRAPH.findRecord(ck);

                    if (cr.getCoverage(childColor) > 0) {
                        log.info("{} {} {} {} {} {} {} {} {} {} {}", id, i, sk, cr.getCortexKmer(), cr.getCoverage(childColor), cr.getCoverage(0), cr.getCoverage(17), cr.getEdgesAsString(childColor), cr.getEdgesAsString(0), cr.getEdgesAsString(17), varSeqSet.get(sk));
                        lk = sk;
                    } else {
                        TraversalEngine e = new TraversalEngineFactory()
                                .traversalDirection(TraversalEngineConfiguration.TraversalDirection.FORWARD)
                                .combinationOperator(TraversalEngineConfiguration.GraphCombinationOperator.OR)
                                .traversalColor(childColor)
                                .recruitmentColors(parentColors)
                                .stopper(new NahrStopper(ROI))
                                .graph(GRAPH)
                                .links(lm)
                                .make();

                        DirectedGraph<CortexVertex, CortexEdge> dfs = e.dfs(lk);

                        String contig = e.getContig(dfs, lk, childColor);

                        for (int j = 0; j <= contig.length() - GRAPH.getKmerSize(); j++) {
                            String skn = contig.substring(j, j + GRAPH.getKmerSize());
                            CortexKmer ckn = new CortexKmer(skn);
                            CortexRecord crn = GRAPH.findRecord(ckn);

                            log.info("{} {} {} {} {} {} {} {} {} {} {}", "recon", j, skn, crn.getCortexKmer(), crn.getCoverage(childColor), crn.getCoverage(0), crn.getCoverage(17), crn.getEdgesAsString(childColor), crn.getEdgesAsString(0), crn.getEdgesAsString(17), varSeqSet.get(skn));

                            visited.add(skn);

                            lk = skn;
                        }

                        log.info("---");
                    }
                }

                visited.add(sk);
            }

            log.info("");
        }

        log.info("");
        */
    }
}
