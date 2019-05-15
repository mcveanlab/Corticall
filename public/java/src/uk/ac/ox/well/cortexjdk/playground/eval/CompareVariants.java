package uk.ac.ox.well.cortexjdk.playground.eval;

import com.google.common.base.Joiner;
import com.google.common.collect.Sets;
import htsjdk.samtools.reference.FastaSequenceFile;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import uk.ac.ox.well.cortexjdk.commands.Module;
import uk.ac.ox.well.cortexjdk.utils.arguments.Argument;
import uk.ac.ox.well.cortexjdk.utils.arguments.Output;
import uk.ac.ox.well.cortexjdk.utils.io.graph.cortex.CortexGraph;
import uk.ac.ox.well.cortexjdk.utils.io.graph.cortex.CortexRecord;
import uk.ac.ox.well.cortexjdk.utils.kmer.CanonicalKmer;
import uk.ac.ox.well.cortexjdk.utils.sequence.SequenceUtils;

import java.io.PrintStream;
import java.util.*;

public class CompareVariants extends Module {
    @Argument(fullName="graph", shortName="g", doc="Graph")
    public CortexGraph GRAPH;

    @Argument(fullName="reference", shortName="R", doc="Reference")
    public ArrayList<FastaSequenceFile> REFERENCE;

    @Argument(fullName="comp", shortName="c", doc="Comp VCF")
    public VCFFileReader COMP;

    @Argument(fullName="eval", shortName="e", doc="Eval VCF")
    public VCFFileReader EVAL;

    @Argument(fullName="kmerSize", shortName="ks", doc="Kmer size")
    public Integer KMER_SIZE = 47;

    @Argument(fullName="windowSize", shortName="ws", doc="Window size")
    public Integer WINDOW_SIZE = 50;

    @Output
    public PrintStream out;

    private List<List<Integer>> generateCombinatoricLists(List<Integer> indices) {
        List<List<Integer>> loos = new ArrayList<>();

        for (int i = 0; i < indices.size(); i++) {
            List<Integer> loo = new ArrayList<>();

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

    private Set<String> recursivelyGenerateCombinations(Set<VariantContext> variants, List<String> alleles, int pos, int seqLength, Map<String, List<String>> allAlleles) {
        List<VariantContext> affectingVariants = new ArrayList<>(variants);
        List<Integer> indices = new ArrayList<>();

        int startPos = pos;
        int endPos = pos;
        for (int i = 0; i < affectingVariants.size(); i++) {
            indices.add(i);

            startPos = Math.min(startPos, affectingVariants.get(i).getStart());
            endPos = Math.max(endPos, affectingVariants.get(i).getStart());
        }

        Set<List<Integer>> loos = new HashSet<>();
        loos.add(new ArrayList<>());
        loos.add(indices);
        loos.addAll(generateCombinatoricLists(indices));

        int start = startPos - WINDOW_SIZE >= 0 ? startPos - WINDOW_SIZE : 0;
        int stop = endPos + WINDOW_SIZE < seqLength ? endPos + WINDOW_SIZE : seqLength - 1;

        Set<String> haplotypes = new HashSet<>();
        for (List<Integer> loo : loos) {
            Map<Integer, VariantContext> vcMap = new HashMap<>();
            for (int index : loo) {
                VariantContext vc = affectingVariants.get(index);
                vcMap.put(vc.getStart() - 1, vc);
            }

            StringBuilder sb = new StringBuilder();
            for (int i = start; i <= stop; i++) {
                if (vcMap.containsKey(i)) {
                    VariantContext vc = vcMap.get(i);

                    if (vc.isSNP() || vc.isSimpleInsertion()) {
                        sb.append(vc.getAltAlleleWithHighestAlleleCount());
                    } else if (vc.isSimpleDeletion()) {
                        sb.append(vc.getAltAlleleWithHighestAlleleCount());
                        i += vc.getReference().length() - 1;
                    } else if (vc.isSymbolicOrSV()) {
                        String[] pieces = vc.getAltAlleleWithHighestAlleleCount().getDisplayString().split("[\\[\\]]");

                        boolean extendDirectionIsLeft = vc.getAltAlleleWithHighestAlleleCount().getDisplayString().contains("]");
                        boolean joinedAfter = pieces[0].matches("^\\w+");
                        boolean isRC = (extendDirectionIsLeft && joinedAfter) || (!extendDirectionIsLeft && !joinedAfter);

                        String t = joinedAfter ? pieces[0] : pieces[pieces.length - 1];

                        int bndStart = vc.getAttributeAsInt("END", 0);
                        int bndEnd = extendDirectionIsLeft ? bndStart - WINDOW_SIZE : bndStart + WINDOW_SIZE;
                        int bndDir = extendDirectionIsLeft ? -1 : 1;

                        StringBuilder bndBuild = new StringBuilder();

                        for (int j = bndStart; j != bndEnd && j > 0 && j < allAlleles.get(vc.getAttributeAsString("CHR2", vc.getContig())).size(); j += bndDir) {
                            bndBuild.append(allAlleles.get(vc.getAttributeAsString("CHR2", vc.getContig())).get(j - 1));
                        }

                        String bnd = bndBuild.toString();
                        if (isRC) {
                            bnd = SequenceUtils.reverseComplement(bnd);
                        }

                        sb.append(t);
                        sb.append(bnd);

                        /*
                        if (vc.hasAttribute("CONSENSUS")) {
                            SmithWaterman sw = new SmithWaterman();
                            String[] res = sw.getAlignment(sb.toString(), vc.getAttributeAsString("CONSENSUS", ""));

                            log.info("{}", res[0]);
                            log.info("{}", res[1]);
                            log.info("");
                        }
                        */
                    }

                    if (vc.hasAttribute("CONSENSUS")) {
                        haplotypes.add(vc.getAttributeAsString("CONSENSUS", ""));
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
        Set<VariantContext> compVcs = new LinkedHashSet<>();
        Set<VariantContext> evalVcs = new HashSet<>();

        log.info("Graph info:");
        int childColor = 0;
        for (int c = 0; c < GRAPH.getNumColors(); c++) {
            if (GRAPH.getSampleName(c).contains("child")) {
                childColor = c;
                break;
            }
        }

        List<Integer> parentColors = new ArrayList<>();
        parentColors.add(GRAPH.getColorForSampleName("PG0004-CW.ERR012788"));
        parentColors.add(GRAPH.getColorForSampleName("PG0008-CW.ERR012840"));

        log.info("  child:   {}", childColor);
        log.info("  parents: {}", Joiner.on(",").join(parentColors));

        log.info("Loading variants...");
        COMP.iterator().stream().forEach(compVcs::add);
        EVAL.iterator().stream().forEach(evalVcs::add);

        log.info("  {} {}", compVcs.size(), evalVcs.size());

        Map<String, List<String>> allAlleles = new LinkedHashMap<>();
        Map<String, String> allSeqs = new HashMap<>();

        log.info("Loading reference...");
        ReferenceSequence rseq;
        for (FastaSequenceFile ref : REFERENCE) {
            log.info("  {}", ref);

            while ((rseq = ref.nextSequence()) != null) {
                String name = rseq.getName().split("\\s+")[0];
                String seq = new String(rseq.getBases());

                List<String> alleles = new ArrayList<>();

                for (int i = 0; i < seq.length(); i++) {
                    alleles.add(seq.substring(i, i + 1));
                }

                allAlleles.put(name, alleles);
                allSeqs.put(name, seq);
            }
        }

        log.info("Processing variants...");
        for (VariantContext compVc : compVcs) {
            Set<String> haplotypes = recursivelyGenerateCombinations(Sets.newHashSet(compVc), allAlleles.get(compVc.getContig()), compVc.getStart(), allSeqs.get(compVc.getContig()).length(), allAlleles);

            log.info("{}", compVc);
            log.info("old: {}", compVc.getAttributeAsString("OLD_HAP", ""));
            log.info("     {}{}{}",
                    compVc.getAttributeAsString("SLEFT", "").toLowerCase(),
                    compVc.getAttributeAsString("SOLD", "").equalsIgnoreCase(".") ? "" : compVc.getAttributeAsString("SOLD", ""),
                    compVc.getAttributeAsString("SRIGHT", "").toLowerCase()
            );
            log.info("new: {}", compVc.getAttributeAsString("NEW_HAP", ""));
            log.info("     {}{}{}",
                    compVc.getAttributeAsString("SLEFT", "").toLowerCase(),
                    compVc.getAttributeAsString("SNEW", "").equalsIgnoreCase(".") ? "" : compVc.getAttributeAsString("SNEW", ""),
                    compVc.getAttributeAsString("SRIGHT", "").toLowerCase()
            );

            int oldKmers = 0;
            int kmersWithParentalCoverage = 0;

            for (int i = 0; i < compVc.getAttributeAsString("OLD_HAP", "").length() - GRAPH.getKmerSize(); i++) {
                String sk = compVc.getAttributeAsString("OLD_HAP", "").substring(i, i + GRAPH.getKmerSize());
                CanonicalKmer ck = new CanonicalKmer(sk);
                CortexRecord cr = GRAPH.findRecord(ck);

                if (cr != null) {
                    boolean coverageInParent = false;
                    for (int c : parentColors) {
                        coverageInParent |= cr.getCoverage(c) > 0;
                    }

                    if (coverageInParent) {
                        kmersWithParentalCoverage++;
                    }
                }

                oldKmers++;
            }

            int newKmers = 0;
            int kmersWithChildCoverage = 0;

            for (int i = 0; i < compVc.getAttributeAsString("NEW_HAP", "").length() - GRAPH.getKmerSize(); i++) {
                String sk = compVc.getAttributeAsString("NEW_HAP", "").substring(i, i + GRAPH.getKmerSize());
                CanonicalKmer ck = new CanonicalKmer(sk);
                CortexRecord cr = GRAPH.findRecord(ck);

                if (cr != null && cr.getCoverage(childColor) > 0) {
                    kmersWithChildCoverage++;
                }

                //log.info("{}", cr);

                newKmers++;
            }

            for (String haplotype : haplotypes) {
                log.info("{}", haplotype);
            }
            log.info("  old: {}/{}, new: {}/{}", kmersWithParentalCoverage, oldKmers, kmersWithChildCoverage, newKmers);
            log.info("");
        }
    }
}

