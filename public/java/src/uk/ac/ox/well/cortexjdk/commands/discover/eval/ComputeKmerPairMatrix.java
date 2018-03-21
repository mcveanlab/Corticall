package uk.ac.ox.well.cortexjdk.commands.discover.eval;

import htsjdk.samtools.reference.FastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import uk.ac.ox.well.cortexjdk.commands.Module;
import uk.ac.ox.well.cortexjdk.utils.arguments.Argument;
import uk.ac.ox.well.cortexjdk.utils.arguments.Output;
import uk.ac.ox.well.cortexjdk.utils.io.table.TableReader;
import uk.ac.ox.well.cortexjdk.utils.kmer.CanonicalKmer;
import uk.ac.ox.well.cortexjdk.utils.progress.ProgressMeter;
import uk.ac.ox.well.cortexjdk.utils.progress.ProgressMeterFactory;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class ComputeKmerPairMatrix extends Module {
    @Argument(fullName="variantTable", shortName="v", doc="Variant table")
    public File VARIANT_TABLE;

    @Argument(fullName="contigs", shortName="c", doc="Contigs")
    public FastaSequenceFile CONTIGS;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        TableReader tr = new TableReader(VARIANT_TABLE, "id", "length", "index", "kmer");

        Map<CanonicalKmer, Integer> ids = new HashMap<>();
        Map<CanonicalKmer, Integer> matrixIndices = new HashMap<>();
        int matrixIndex = 0;

        for (Map<String, String> te : tr) {
            CanonicalKmer ck = new CanonicalKmer(te.get("kmer"));
            int id = Integer.valueOf(te.get("id"));

            ids.put(ck, id);
            matrixIndices.put(ck, matrixIndex);

            matrixIndex++;
        }

        int[][] m = new int[matrixIndex][matrixIndex];
        int kmerSize = ids.keySet().iterator().next().length();

        List<ReferenceSequence> rseqs = new ArrayList<>();
        ReferenceSequence aseq;
        while ((aseq = CONTIGS.nextSequence()) != null) {
            rseqs.add(aseq);
        }

//        ProgressMeter pm = new ProgressMeterFactory()
//                .maxRecord(rseqs.size())

        for (ReferenceSequence rseq : rseqs) {
            log.info("{}", rseq.getName());

            String seq = rseq.getBaseString();

            List<CanonicalKmer> cks = new ArrayList<>();
            for (int i = 0; i <= seq.length() - kmerSize; i++) {
                CanonicalKmer ck = new CanonicalKmer(seq.substring(i, i + kmerSize));

                if (ids.containsKey(ck)) {
                    cks.add(ck);
                }
            }

            for (int i = 0; i < cks.size(); i++) {
                for (int j = 0; j < cks.size(); j++) {
                    boolean expectedAdjacent = ids.get(cks.get(i)).equals(ids.get(cks.get(j)));

                    int idxi = matrixIndices.get(cks.get(i));
                    int idxj = matrixIndices.get(cks.get(j));

                    m[idxi][idxj] = expectedAdjacent ? 1 : -1;
                    m[idxj][idxi] = expectedAdjacent ? 1 : -1;
                }
            }
        }

        for (int i = 0; i < matrixIndex; i++) {
            for (int j = 0; j < matrixIndex; j++) {
                out.print(m[i][j]);
                if (j < matrixIndex - 1) {
                    out.print("\t");
                }
            }
            out.println();
        }
    }
}
