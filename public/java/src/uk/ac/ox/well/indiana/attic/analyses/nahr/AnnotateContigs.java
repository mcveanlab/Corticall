package uk.ac.ox.well.indiana.attic.analyses.nahr;

import com.google.common.base.Joiner;
import net.sf.picard.reference.FastaSequenceFile;
import net.sf.picard.reference.ReferenceSequence;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.cortex.CortexKmer;
import uk.ac.ox.well.indiana.utils.io.cortex.CortexMap;
import uk.ac.ox.well.indiana.utils.io.cortex.CortexRecord;
import uk.ac.ox.well.indiana.utils.sequence.SequenceUtils;

import java.io.PrintStream;
import java.util.*;

public class AnnotateContigs extends Module {
    @Argument(fullName="contigs", shortName="c", doc="Contigs (in FASTA format)")
    public FastaSequenceFile CONTIGS;

    @Argument(fullName="cortexGraph", shortName="cg", doc="Cortex graph for sample")
    public CortexMap SAMPLE;

    @Argument(fullName="parents", shortName="p", doc="Parents (in Cortex format)")
    public CortexMap PARENTS;

    //@Argument(fullName="parentsSupp", shortName="ps", doc="Supplemental parental information")
    //public CortexMap PARENTS_SUPP;

    @Output
    public PrintStream out;

    private boolean isContiguous(String prevStr, String curStr, int color) {
        CortexKmer prevKmer = new CortexKmer(prevStr);
        CortexKmer curKmer = new CortexKmer(curStr);

        CortexRecord prevCr = PARENTS.get(prevKmer);
        CortexRecord curCr = PARENTS.get(curKmer);

        /*
        if (prevCr == null) {
            prevCr = PARENTS_SUPP.get(prevKmer);
        }

        if (curCr == null) {
            curCr = PARENTS_SUPP.get(curKmer);
        }
        */

        if (prevCr != null && curCr != null) {
            Set<CortexKmer> computedCurOutKmers = new HashSet<CortexKmer>();
            if (!prevKmer.isFlipped()) {
                for (String outEdge : prevCr.getOutEdgesAsStrings(color)) {
                    CortexKmer outKmer = new CortexKmer(prevKmer.getKmerAsString().substring(1, prevKmer.length()) + outEdge);
                    computedCurOutKmers.add(outKmer);
                }
            } else {
                for (String outEdge : prevCr.getInEdgesAsStrings(color)) {
                    CortexKmer outKmer = new CortexKmer(outEdge + prevKmer.getKmerAsString().substring(0, prevKmer.length() - 1));
                    computedCurOutKmers.add(outKmer);
                }
            }

            Set<CortexKmer> computedPrevInKmers = new HashSet<CortexKmer>();
            if (!curKmer.isFlipped()) {
                for (String inEdge : curCr.getInEdgesAsStrings(color)) {
                    CortexKmer inKmer = new CortexKmer(inEdge + curKmer.getKmerAsString().substring(0, curKmer.length() - 1));
                    computedPrevInKmers.add(inKmer);
                }
            } else {
                for (String inEdge : curCr.getOutEdgesAsStrings(color)) {
                    CortexKmer inKmer = new CortexKmer(curKmer.getKmerAsString().substring(1, curKmer.length()) + inEdge);
                    computedPrevInKmers.add(inKmer);
                }
            }

            return (computedCurOutKmers.contains(curKmer) && computedPrevInKmers.contains(prevKmer));
        }

        return false;
    }

    @Override
    public void execute() {
        int kmerSize = PARENTS.keySet().iterator().next().length();

        out.println("contigName\tseq\tkmerOrigin\tkmerContiguity\tkmerCoverage");

        ReferenceSequence rseq;
        while ((rseq = CONTIGS.nextSequence()) != null) {
            if (rseq.length() > 0) {
                String seq = new String(rseq.getBases());

                StringBuilder annotation = new StringBuilder();
                List<String> coverages = new ArrayList<String>();

                for (int i = 0; i <= seq.length() - kmerSize; i++) {
                    CortexKmer kmer = new CortexKmer(seq.substring(i, i + kmerSize));

                    coverages.add(String.valueOf(SAMPLE.get(kmer).getCoverage(0)));

                    if (!PARENTS.containsKey(kmer)) {
                        annotation.append(".");
                    } else {
                        CortexRecord cr = PARENTS.get(kmer);

                        if (cr.getCoverage(0) == 0 && cr.getCoverage(1) > 0) {
                            annotation.append("1");
                        } else if (cr.getCoverage(0) > 0 && cr.getCoverage(1) == 0) {
                            annotation.append("0");
                        } else {
                            annotation.append("B");
                        }
                    }
                }

                StringBuilder contiguity = new StringBuilder();
                contiguity.append("0");
                for (int i = 1; i <= seq.length() - kmerSize; i++) {
                    String prevStr = seq.substring(i - 1, i - 1 + kmerSize);
                    String curStr  = seq.substring(i, i + kmerSize);

                    char colorChar = annotation.charAt(i);

                    boolean isContiguous;
                    switch (colorChar) {
                        case '0': isContiguous = isContiguous(prevStr, curStr, 0); break;
                        case '1': isContiguous = isContiguous(prevStr, curStr, 1); break;
                        case 'B': isContiguous = isContiguous(prevStr, curStr, 0) || isContiguous(prevStr, curStr, 1); break;
                        default:  isContiguous = false; break;
                    }

                    contiguity.append(isContiguous ? "1" : "0");
                }

                String coverage = Joiner.on(",").join(coverages);

                out.println(rseq.getName().split("\\s+")[0] + "\t" + seq + "\t" + annotation.toString() + "\t" + contiguity.toString() + "\t" + coverage);
            }
        }
    }
}
