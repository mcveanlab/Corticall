package uk.ac.ox.well.cortexjdk.commands.discover.eval;

import com.google.common.base.Joiner;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.reference.FastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.util.Interval;
import org.apache.commons.lang3.tuple.Triple;
import org.jetbrains.annotations.NotNull;
import org.jgrapht.GraphPath;
import org.jgrapht.alg.shortestpath.DijkstraShortestPath;
import org.jgrapht.graph.DirectedWeightedPseudograph;
import uk.ac.ox.well.cortexjdk.commands.Module;
import uk.ac.ox.well.cortexjdk.utils.alignment.reference.IndexedReference;
import uk.ac.ox.well.cortexjdk.utils.arguments.Argument;
import uk.ac.ox.well.cortexjdk.utils.arguments.Output;
import uk.ac.ox.well.cortexjdk.utils.io.graph.cortex.CortexGraph;
import uk.ac.ox.well.cortexjdk.utils.io.graph.cortex.CortexRecord;
import uk.ac.ox.well.cortexjdk.utils.io.graph.links.CortexLinks;
import uk.ac.ox.well.cortexjdk.utils.io.table.TableReader;
import uk.ac.ox.well.cortexjdk.utils.kmer.CanonicalKmer;
import uk.ac.ox.well.cortexjdk.utils.kmer.CortexByteKmer;
import uk.ac.ox.well.cortexjdk.utils.progress.ProgressMeter;
import uk.ac.ox.well.cortexjdk.utils.progress.ProgressMeterFactory;
import uk.ac.ox.well.cortexjdk.utils.stoppingrules.BubbleClosingStopper;
import uk.ac.ox.well.cortexjdk.utils.stoppingrules.NovelContinuationStopper;
import uk.ac.ox.well.cortexjdk.utils.traversal.*;

import java.io.File;
import java.io.PrintStream;
import java.util.*;

import static uk.ac.ox.well.cortexjdk.utils.traversal.TraversalEngineConfiguration.GraphCombinationOperator.OR;
import static uk.ac.ox.well.cortexjdk.utils.traversal.TraversalEngineConfiguration.TraversalDirection.FORWARD;

/**
 * Created by kiran on 30/08/2017.
 */
public class CompareContigLengths extends Module {
    @Argument(fullName="variants", shortName="v", doc="Variants")
    public File VARIANTS;

    @Argument(fullName="kmers", shortName="k", doc="Kmers")
    public File KMERS;

    @Argument(fullName="fasta", shortName="f", doc="Fasta")
    public LinkedHashMap<String, FastaSequenceFile> FASTAS;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        TableReader tr = new TableReader(VARIANTS);

        Map<String, String> ids = new HashMap<>();
        for (Map<String, String> te : tr) {
            if (!te.get("type").equals("RECOMB")) {
                ids.put(te.get("index"), te.get("type"));
            }
        }

        Map<CanonicalKmer, Triple<String, String, Map<String, Integer>>> m = new LinkedHashMap<>();

        TableReader tk = new TableReader(KMERS, "index", "numKmers", "kmerIndex", "kmer");
        for (Map<String, String> te : tk) {
            m.put(new CanonicalKmer(te.get("kmer")), Triple.of(te.get("numKmers"), ids.get(te.get("index")), new LinkedHashMap<>()));
        }

        for (String key : FASTAS.keySet()) {
            FastaSequenceFile fa = FASTAS.get(key);
            ReferenceSequence rseq;
            while ((rseq = fa.nextSequence()) != null) {
                String[] names = rseq.getName().split(" ");

                for (String name : names) {
                    if (name.startsWith("seed=")) {
                        String[] pieces = name.split("=");
                        String seed = pieces[1];
                        CanonicalKmer ck = new CanonicalKmer(seed);

                        if (m.containsKey(ck)) {
                            m.get(ck).getRight().put(key, rseq.length());
                        }
                    }
                }
            }
        }

        out.println(Joiner.on("\t").join("kmer", "id", "type", Joiner.on("\t").join(FASTAS.keySet())));
        for (CanonicalKmer ck : m.keySet()) {
            Triple<String, String, Map<String, Integer>> t = m.get(ck);

            List<String> fields = new ArrayList<>();
            fields.add(ck.getKmerAsString());
            fields.add(t.getLeft());
            fields.add(t.getMiddle());

            for (String key : FASTAS.keySet()) {
                fields.add(String.valueOf(t.getRight().get(key)));
            }

            out.println(Joiner.on("\t").join(fields));

            //log.info("{} {} {} {}", ck, t.getLeft(), t.getMiddle(), t.getRight());
        }
    }
}
