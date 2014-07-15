package uk.ac.ox.well.indiana.attic.analyses.geo;

import com.google.common.base.Joiner;
import htsjdk.samtools.reference.FastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.containers.DataFrame;

import java.io.PrintStream;
import java.util.*;

public class ComputeSequenceRelatedness extends Module {
    @Argument(fullName="sequences", shortName="s", doc="Sequences")
    public FastaSequenceFile FASTA;

    @Argument(fullName="kmerSize", shortName="ks", doc="Kmer size")
    public Integer KMER_SIZE = 11;

    @Argument(fullName="seqNames", shortName="sn", doc="Seq names", required=false)
    public HashSet<String> SEQ_NAMES;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        Map<String, Set<String>> kmerToSeqMap = new HashMap<String, Set<String>>();

        ReferenceSequence rseq;
        while ((rseq = FASTA.nextSequence()) != null) {
            String[] names = rseq.getName().split("\\s+");
            String name = names[0];

            if (SEQ_NAMES == null || SEQ_NAMES.size() == 0 || SEQ_NAMES.contains(name)) {
                String bases = new String(rseq.getBases());

                for (int i = 0; i <= bases.length() - KMER_SIZE; i++) {
                    String kmer = bases.substring(i, i + KMER_SIZE);

                    if (!kmerToSeqMap.containsKey(kmer)) {
                        kmerToSeqMap.put(kmer, new HashSet<String>());
                    }

                    kmerToSeqMap.get(kmer).add(name);
                }
            }
        }

        DataFrame<String, String, Float> rmat = new DataFrame<String, String, Float>(0.0f);

        Map<String, Float> totalKmers = new HashMap<String, Float>();

        int kmersProcessed = 0;
        for (String kmer : kmerToSeqMap.keySet()) {
            if (kmersProcessed % (kmerToSeqMap.size() / 10) == 0) {
                log.info("Processed {}/{} kmers", kmersProcessed, kmerToSeqMap.size());
            }
            kmersProcessed++;

            for (String name1 : kmerToSeqMap.get(kmer)) {
                for (String name2 : kmerToSeqMap.get(kmer)) {
                    rmat.set(name1, name2, rmat.get(name1, name2) + 1.0f);
                }

                if (!totalKmers.containsKey(name1)) {
                    totalKmers.put(name1, 0.0f);
                }
                totalKmers.put(name1, totalKmers.get(name1) + 1.0f);
            }
        }

        Collection<String> rowNames = rmat.getRowNames();
        Collection<String> colNames = rmat.getColNames();

        out.println("\t" + Joiner.on("\t").join(colNames));

        for (String name1 : rowNames) {
            Collection<String> fields = new ArrayList<String>();
            fields.add(name1);

            for (String name2 : colNames) {
                float intersection = rmat.get(name1, name2);
                float total = intersection / (totalKmers.get(name1) < totalKmers.get(name2) ? totalKmers.get(name1) : totalKmers.get(name2));
                float distance = 1.0f - total;

                //rmat.set(name1, name2, 1.0f - total);
                //rmat.set(name2, name1, 1.0f - total);

                fields.add(String.valueOf(distance));
            }

            out.println(Joiner.on("\t").join(fields));
        }

        //out.println(rmat);
    }
}
