package uk.ac.ox.well.indiana.commands.variants;

import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import org.mapdb.Atomic;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;

import java.io.PrintStream;
import java.util.*;

public class KmersFromVariants extends Module {
    @Argument(fullName="reference", shortName="r", doc="Reference")
    public IndexedFastaSequenceFile REFERENCE;

    @Argument(fullName="vcf", shortName="v", doc="VCF")
    public VCFFileReader VCF;

    @Argument(fullName="kmerSize", shortName="ks", doc="Kmer size")
    public Integer KMER_SIZE = 47;

    @Argument(fullName="windowSize", shortName="w", doc="Window size")
    public Integer WINDOW_SIZE = 50;

    @Output
    public PrintStream out;

    private List<List<Integer>> generateCombinatoricLists(List<Integer> indices) {
        List<List<Integer>> loos = new ArrayList<List<Integer>>();

        for (int i = 0; i < indices.size(); i++) {
            List<Integer> loo = new ArrayList<Integer>();

            for (int j = 0; j < indices.size(); j++) {
                if (i != j) {
                    loo.add(indices.get(j));
                }
            }

            if (loo.size() > 0) {
                loos.add(loo);
            }

            if (loo.size() > 1) {
                loos.addAll(generateCombinatoricLists(loo));
            }
        }

        return loos;
    }

    private Set<String> recursivelyGenerateCombinations(List<VariantContext> affectingVariants, List<String> alleles, int pos) {
        List<Integer> indices = new ArrayList<Integer>();

        for (int i = 0; i < affectingVariants.size(); i++) { indices.add(i); }

        Set<List<Integer>> loos = new HashSet<List<Integer>>();
        loos.add(new ArrayList<Integer>());
        loos.add(indices);
        loos.addAll(generateCombinatoricLists(indices));

        int start = pos - WINDOW_SIZE >= 0 ? pos - WINDOW_SIZE : 0;
        int stop = pos + WINDOW_SIZE < alleles.size() ? pos + WINDOW_SIZE : alleles.size();

        Set<String> haplotypes = new HashSet<String>();
        for (List<Integer> loo : loos) {
            Map<Integer, VariantContext> vcMap = new HashMap<Integer, VariantContext>();
            for (int index : loo) {
                VariantContext vc = affectingVariants.get(index);
                vcMap.put(vc.getStart() - 1, vc);
            }

            StringBuilder sb = new StringBuilder();
            for (int i = start; i <= stop; i++) {
                if (vcMap.containsKey(i)) {
                    VariantContext vc = vcMap.get(i);

                    if (vc.isSNP() || vc.isSimpleInsertion()) {
                        sb.append(vc.getGenotype(0).getGenotypeString());
                    } else if (vc.isSimpleDeletion()) {
                        sb.append(vc.getGenotype(0).getGenotypeString());
                        i += vc.getReference().length();
                    }
                } else {
                    sb.append(alleles.get(i));
                }
            }

            haplotypes.add(sb.toString());
        }

        return haplotypes;
    }

    @Override
    public void execute() {
        Map<String, Map<Integer, Set<VariantContext>>> allVcs = new HashMap<String, Map<Integer, Set<VariantContext>>>();

        for (VariantContext vc : VCF) {
            if (!vc.getGenotype(0).isHomRef() && !vc.getGenotype(0).getGenotypeString().equals(".")) {
                if (!allVcs.containsKey(vc.getChr())) {
                    allVcs.put(vc.getChr(), new HashMap<Integer, Set<VariantContext>>());
                }

                if (!allVcs.get(vc.getChr()).containsKey(vc.getStart() - 1)) {
                    allVcs.get(vc.getChr()).put(vc.getStart() - 1, new HashSet<VariantContext>());
                }

                allVcs.get(vc.getChr()).get(vc.getStart() - 1).add(vc);
            }
        }

        ReferenceSequence rseq;
        while ((rseq = REFERENCE.nextSequence()) != null) {
            String name = rseq.getName().split("\\s+")[0];
            String seq = new String(rseq.getBases());

            List<String> alleles = new ArrayList<String>();

            for (int i = 0; i < seq.length(); i++) {
                alleles.add(seq.substring(i, i + 1));
            }

            for (int i = 0; i <= seq.length() - KMER_SIZE; i++) {
                Set<VariantContext> vcs = allVcs.containsKey(name) && allVcs.get(name).containsKey(i) ? allVcs.get(name).get(i) : null;

                if (vcs != null) {
                    Set<VariantContext> affectingVariantsSet = new HashSet<VariantContext>(vcs);

                    for (int j = i; j >= 0 && j >= i - WINDOW_SIZE; j--) {
                        if (allVcs.containsKey(name) && allVcs.get(name).containsKey(j)) {
                            affectingVariantsSet.addAll(allVcs.get(name).get(j));
                        }
                    }

                    for (int j = i; j < seq.length() && j <= i + WINDOW_SIZE; j++) {
                        if (allVcs.containsKey(name) && allVcs.get(name).containsKey(j)) {
                            affectingVariantsSet.addAll(allVcs.get(name).get(j));
                        }
                    }

                    List<VariantContext> affectingVariants = new ArrayList<VariantContext>(affectingVariantsSet);

                    log.info("{} {} {} {}", name, i, alleles.get(i), affectingVariantsSet.size());
                    for (VariantContext vc : affectingVariantsSet) {
                        //log.info("  {}", vc);
                    }

                    Set<String> haplotypes = affectingVariants.size() < 11 ? recursivelyGenerateCombinations(affectingVariants, alleles, i) : new HashSet<String>();

                    Set<String> kmers = new HashSet<String>();
                    for (String haplotype : haplotypes) {
                        for (int j = 0; j <= haplotype.length() - KMER_SIZE; j++) {
                            String kmer = haplotype.substring(j, j + KMER_SIZE);
                            kmers.add(kmer);
                        }
                    }

                    int kindex = 0;
                    for (String kmer : kmers) {
                        out.println(">" + name + ":" + i + "." + kindex);
                        out.println(kmer);

                        kindex++;
                    }

                    int nextpos = i;
                    for (VariantContext vc : affectingVariants) {
                        if (vc.getStart() - 1 > nextpos) {
                            nextpos = vc.getStart() - 1;
                        }
                    }

                    i = nextpos + 1;
                }
            }
        }
    }
}
