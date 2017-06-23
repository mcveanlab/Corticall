package uk.ac.ox.well.indiana.commands.caller.call;

import htsjdk.samtools.reference.FastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.alignment.kmer.KmerLookup;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexKmer;
import uk.ac.ox.well.indiana.utils.sequence.SequenceUtils;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;

/**
 * Created by kiran on 23/06/2017.
 */
public class MergeContigs extends Module {
    @Argument(fullName="graph", shortName="g", doc="Graph")
    public CortexGraph GRAPH;

    @Argument(fullName="child", shortName="c", doc="Child")
    public String CHILD;

    @Argument(fullName="parents", shortName="p", doc="Parents")
    public ArrayList<String> PARENTS;

    @Argument(fullName="sequences", shortName="s", doc="Contigs")
    public FastaSequenceFile CONTIGS;

    @Argument(fullName="drafts", shortName="d", doc="Drafts")
    public HashMap<String, KmerLookup> LOOKUPS;

    @Override
    public void execute() {
        Map<CortexKmer, String> contigs = new HashMap<>();
        ReferenceSequence fw;
        while ((fw = CONTIGS.nextSequence()) != null) {
            ReferenceSequence rc = new ReferenceSequence(fw.getName(), fw.getContigIndex(), SequenceUtils.reverseComplement(fw.getBases()));

            String fwFirst = fw.getBaseString().substring(0, GRAPH.getKmerSize());
            String fwLast = fw.getBaseString().substring(fw.length() - GRAPH.getKmerSize(), fw.length());

            String rcFirst = rc.getBaseString().substring(0, GRAPH.getKmerSize());
            String rcLast = rc.getBaseString().substring(rc.length() - GRAPH.getKmerSize(), rc.length());

            log.info("{}", fwFirst);
            log.info("{}", fwLast);
            log.info("{}", rcFirst);
            log.info("{}", rcLast);
        }
    }
}
