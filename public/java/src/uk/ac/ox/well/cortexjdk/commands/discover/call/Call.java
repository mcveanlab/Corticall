package uk.ac.ox.well.cortexjdk.commands.discover.call;

import com.google.common.base.Joiner;
import htsjdk.samtools.*;
import htsjdk.samtools.reference.FastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.StringUtil;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.SortingVariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.VCFHeader;
import org.apache.commons.lang3.tuple.Triple;
import org.apache.commons.math3.util.Pair;
import org.jetbrains.annotations.NotNull;
import org.jgrapht.DirectedGraph;
import org.jgrapht.Graph;
import org.jgrapht.Graphs;
import org.jgrapht.alg.ConnectivityInspector;
import org.jgrapht.graph.DefaultDirectedGraph;
import org.jgrapht.graph.DefaultEdge;
import org.jgrapht.graph.DirectedWeightedPseudograph;
import uk.ac.ox.well.cortexjdk.commands.Module;
import uk.ac.ox.well.cortexjdk.utils.alignment.mosaic.MosaicAligner;
import uk.ac.ox.well.cortexjdk.utils.alignment.reference.IndexedReference;
import uk.ac.ox.well.cortexjdk.utils.alignment.sw.SmithWaterman;
import uk.ac.ox.well.cortexjdk.utils.arguments.Argument;
import uk.ac.ox.well.cortexjdk.utils.arguments.Output;
import uk.ac.ox.well.cortexjdk.utils.containers.ContainerUtils;
import uk.ac.ox.well.cortexjdk.utils.io.graph.cortex.CortexGraph;
import uk.ac.ox.well.cortexjdk.utils.io.graph.cortex.CortexRecord;
import uk.ac.ox.well.cortexjdk.utils.io.graph.links.CortexLinks;
import uk.ac.ox.well.cortexjdk.utils.kmer.CanonicalKmer;
import uk.ac.ox.well.cortexjdk.utils.kmer.CortexByteKmer;
import uk.ac.ox.well.cortexjdk.utils.sequence.SequenceUtils;
import uk.ac.ox.well.cortexjdk.utils.stoppingrules.ContigStopper;
import uk.ac.ox.well.cortexjdk.utils.stoppingrules.DestinationStopper;
import uk.ac.ox.well.cortexjdk.utils.traversal.*;

import java.io.File;
import java.io.PrintStream;
import java.util.*;

import static htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder.OutputType.VCF;
import static uk.ac.ox.well.cortexjdk.utils.traversal.TraversalEngineConfiguration.GraphCombinationOperator.OR;
import static uk.ac.ox.well.cortexjdk.utils.traversal.TraversalEngineConfiguration.TraversalDirection.FORWARD;
import static uk.ac.ox.well.cortexjdk.utils.traversal.TraversalEngineConfiguration.TraversalDirection.REVERSE;

public class Call extends Module {
    @Argument(fullName = "graph", shortName = "g", doc = "Graph")
    public CortexGraph GRAPH;

    @Argument(fullName = "links", shortName = "l", doc = "Links", required=false)
    public ArrayList<CortexLinks> LINKS;

    @Argument(fullName = "rois", shortName = "r", doc = "Rois")
    public CortexGraph ROIS;

    @Argument(fullName="partitions", shortName="p", doc="Partitions")
    public FastaSequenceFile PARTITIONS;

    @Argument(fullName="mother", shortName="m", doc="Mother's sample name")
    public LinkedHashSet<String> MOTHER;

    @Argument(fullName="father", shortName="f", doc="Father's sample name")
    public LinkedHashSet<String> FATHER;

    @Argument(fullName="background", shortName="b", doc="Background", required=false)
    public HashMap<String, IndexedReference> BACKGROUNDS;

    @Argument(fullName="contigName", shortName="cn", doc="Contigs to process", required=false)
    public HashSet<String> CONTIG_NAMES;

    //@Argument(fullName="reference", shortName="R", doc="Reference", required=false)
    //public HashMap<String, IndexedReference> REFERENCE;

    @Output
    public File out;

    @Override
    public void execute() {
        Set<CanonicalKmer> rois = loadRois(ROIS);
        List<ReferenceSequence> rseqs = loadPartitions();

        List<SAMSequenceRecord> ssrs = new ArrayList<>();
        for (String id : BACKGROUNDS.keySet()) {
            IndexedReference ir = BACKGROUNDS.get(id);
            ssrs.addAll(ir.getReferenceSequence().getSequenceDictionary().getSequences());
            ssrs.add(new SAMSequenceRecord(id + "_unknown", rseqs.size()));
        }
        SAMSequenceDictionary ssd = new SAMSequenceDictionary(ssrs);

        Map<String, Integer> sid = new HashMap<>();
        for (int i = 0; i < ssrs.size(); i++) {
            sid.put(ssrs.get(i).getSequenceName(), i);
        }

        VariantContextWriter vcw = new VariantContextWriterBuilder()
                .setOutputFile(out)
                .setOutputFileType(VCF)
                .setOption(Options.DO_NOT_WRITE_GENOTYPES)
                .setOption(Options.ALLOW_MISSING_FIELDS_IN_HEADER)
                .unsetOption(Options.INDEX_ON_THE_FLY)
                .build();

        VCFHeader vcfHeader = new VCFHeader();
        vcfHeader.setSequenceDictionary(ssd);
        vcw.writeHeader(vcfHeader);

        MosaicAligner ma = new MosaicAligner(0.15, 0.99, 0.0001,  0.001);
        MosaicAligner mb = new MosaicAligner(0.35, 0.99, 1e-5, 0.001);

        Set<VariantContext> svcs = new TreeSet<>((v1, v2) -> {
            int sid0 = sid.get(v1.getContig());
            int sid1 = sid.get(v2.getContig());

            if (sid0 != sid1) { return sid0 < sid1 ? -1 : 1; }
            if (v1.getStart() != v2.getStart()) { return v1.getStart() < v2.getStart() ? -1 : 1; }

            return 0;
        });

        for (int rseqIndex = 0; rseqIndex < rseqs.size(); rseqIndex++) {
            ReferenceSequence rseq = rseqs.get(rseqIndex);

            List<CortexVertex> w = loadChildWalk(rseq, GRAPH);
            List<Triple<Integer, Integer, List<CortexVertex>>> sections = sectionContig(rois, w, 100, 500);

            log.info("Partition {}/{} (sections={}, fullname={})", (rseqIndex+1), rseqs.size(), sections.size(), rseq.getName());

            for (int sectionIndex = 0; sectionIndex < sections.size(); sectionIndex++) {
                Triple<Integer, Integer, List<CortexVertex>> section = sections.get(sectionIndex);

                List<CortexVertex> ws = section.getRight();
                List<Pair<Integer, Integer>> regions = getRegions(rois, ws);

                Map<String, String> targets = new HashMap<>();
                for (Set<String> parentName : Arrays.asList(MOTHER, FATHER)) {
                    Map<String, String> unlabeledTargets = assembleCandidateHaplotypes(rois, ws, regions, parentName);
                    Map<String, String> labeledTargets = labelTargets(unlabeledTargets);

                    targets.putAll(labeledTargets);
                }

                if (targets.size() > 0) {
                    String trimmedQuery = trimQuery(ws, targets, rois);

                    List<Triple<String, String, Pair<Integer, Integer>>> lps = mb.align(trimmedQuery, targets);
                    log.info("\n{}\n{}", makeNoveltyTrack(rois, trimmedQuery, lps), mb);

                    //List<Triple<String, String, Pair<Integer, Integer>>> lpsb = mb.align(trimmedQuery, targets);
                    //log.info("\n{}\n{}", makeNoveltyTrack(rois, trimmedQuery, lpsb), mb);

                    List<Pair<Integer, Integer>> nrs = getNoveltyRegions(rois, trimmedQuery, lps);
                    List<VariantContext> vcs = getDNMList(lps, nrs, rseqIndex, sectionIndex, rois);

                    for (VariantContext vc : vcs) {
                        Map<String, Object> attrs = new HashMap<>(vc.getAttributes());
                        attrs.remove("NOVELS");
                        log.info("{}:{}-{} type={} alleles=[{}] attr={{}}",
                                vc.getContig(),
                                vc.getStart(),
                                vc.getEnd(),
                                vc.getType(),
                                Joiner.on(", ").join(vc.getAlleles()),
                                Joiner.on(", ").withKeyValueSeparator('=').join(attrs)
                        );

                        svcs.add(vc);
                    }
                    log.info("");
                }
            }
        }

        for (VariantContext vc : svcs) {
            vcw.add(vc);
        }

        vcw.close();
    }

    @NotNull
    private List<ReferenceSequence> loadPartitions() {
        List<ReferenceSequence> rseqs = new ArrayList<>();
        ReferenceSequence rseq;
        while ((rseq = PARTITIONS.nextSequence()) != null) {
            String[] name = rseq.getName().split(" ");
            if (CONTIG_NAMES == null || CONTIG_NAMES.contains(name[0])) {
                rseqs.add(rseq);
            }
        }
        return rseqs;
    }

    private Map<String, String> labelTargets(Map<String, String> newTargets) {
        Map<String, String> targets = new HashMap<>();
        for (String label : newTargets.keySet()) {
            String background = label.split(":")[0];
            String seq = newTargets.get(label);

            String newLabel = label;

            if (BACKGROUNDS.containsKey(background)) {
                SAMRecord a = getBestAlignment(seq, background);

                if (a != null) {
                    newLabel = String.format("%s:%s:%d-%d:%s:%s", background, a.getReferenceName(), a.getAlignmentStart(), a.getAlignmentEnd(), a.getReadNegativeStrandFlag() ? "-" : "+", a.getCigarString());
                }
            }

            targets.put(newLabel, seq);
        }

        return targets;
    }

    private boolean hasAlignedBase(List<Triple<String, String, Pair<Integer, Integer>>> lps, int j) {
        for (int i = 1; i < lps.size(); i++) {
            if (j < lps.get(i).getMiddle().length() && lps.get(0).getMiddle().charAt(j) == lps.get(i).getMiddle().charAt(j)) {
                return true;
            }
        }

        return false;
    }

    private List<VariantContext> getDNMList(List<Triple<String, String, Pair<Integer, Integer>>> lps, List<Pair<Integer, Integer>> nrs, int partitionIndex, int sectionIndex, Set<CanonicalKmer> rois) {
        List<VariantContext> vcs = new ArrayList<>();

        Set<String> bs = new HashSet<>();
        for (String ba : BACKGROUNDS.keySet()) {
            for (SAMSequenceRecord ssr : BACKGROUNDS.get(ba).getReferenceSequence().getSequenceDictionary().getSequences()) {
                bs.add(ssr.getSequenceName());
            }
        }

        for (Pair<Integer, Integer> nr : nrs) {
            String background = "";
            int bi = -1;
            int alleleStart = -1;
            int alleleEnd = -1;
            StringBuilder childAllele = new StringBuilder();
            StringBuilder parentAllele = new StringBuilder();
            boolean isAtEnd = false;
            boolean isAtStart = false;
            boolean isExactRepeat = false;
            List<CanonicalKmer> cks = new ArrayList<>();


            int leftLimit = nr.getFirst();
            while (leftLimit > 0 && !hasAlignedBase(lps, leftLimit)) {
                leftLimit--;
            }

            int rightLimit = nr.getSecond();
            while (rightLimit < lps.get(0).getMiddle().length() - 1 && !hasAlignedBase(lps, rightLimit)) {
                rightLimit++;
            }

            for (int j = leftLimit; j <= rightLimit || (j < lps.get(0).getMiddle().length() && (childAllele.length() > 0 || parentAllele.length() > 0)); j++) {
                String novelRegion = lps.get(0).getMiddle().substring(leftLimit, rightLimit).replaceAll("[- ]", "");
                for (int i = 0; i <= novelRegion.length() - GRAPH.getKmerSize(); i++) {
                    CanonicalKmer ck = new CanonicalKmer(novelRegion.substring(i, i + GRAPH.getKmerSize()));
                    if (rois.contains(ck)) {
                        cks.add(ck);
                    }
                }

                for (int i = 1; i < lps.size(); i++) {
                    if (j < lps.get(i).getMiddle().length() && lps.get(i).getMiddle().charAt(j) != ' ') {
                        if (lps.get(0).getMiddle().charAt(j) != lps.get(i).getMiddle().charAt(j)) {
                            if (j == 0) { isAtStart = true; }
                            if (j == lps.get(i).getMiddle().length() - 1) { isAtEnd = true; }

                            if (alleleStart == -1) {
                                bi = i;
                                background = lps.get(bi).getLeft();
                                alleleStart = j;
                            }
                            alleleEnd = j;

                            if (lps.get(0).getMiddle().charAt(j) != '-') {
                                childAllele.append(lps.get(0).getMiddle().charAt(j));
                            }

                            if (lps.get(bi).getMiddle().charAt(j) != '-') {
                                parentAllele.append(lps.get(bi).getMiddle().charAt(j));
                            }
                        } else {
                            if (alleleStart >= 0) {
                                if (!childAllele.toString().equals(parentAllele.toString()) && childAllele.toString().equals(parentAllele.toString().toUpperCase())) {
                                    parentAllele = new StringBuilder();
                                    isExactRepeat = true;
                                }

                                boolean positiveStrand = BACKGROUNDS == null || background.contains(":+:");
                                int variantStart = positiveStrand ? alleleStart : alleleEnd;
                                if (parentAllele.length() != 1 || childAllele.length() != 1) {
                                    variantStart = positiveStrand ? alleleStart - 1 : alleleEnd + 1;
                                }

                                if (parentAllele.length() != 1 || childAllele.length() != 1) {
                                    if (positiveStrand) {
                                        if (variantStart >= 0 && variantStart < lps.get(bi).getMiddle().length()) {
                                            parentAllele.insert(0, lps.get(bi).getMiddle().charAt(variantStart));
                                            childAllele.insert(0, lps.get(0).getMiddle().charAt(variantStart));
                                            if (isAtEnd) {
                                                childAllele.append(".");
                                            }
                                        } else {
                                            variantStart = alleleEnd + 1;
                                            parentAllele.append(lps.get(bi).getMiddle().charAt(variantStart));
                                            childAllele.append(lps.get(0).getMiddle().charAt(variantStart));
                                            childAllele.insert(0, ".");
                                        }
                                    } else {
                                        if (variantStart >= 0 && variantStart < lps.get(bi).getMiddle().length()) {
                                            parentAllele.append(lps.get(bi).getMiddle().charAt(variantStart));
                                            childAllele.append(lps.get(0).getMiddle().charAt(variantStart));
                                            if (isAtStart) {
                                                childAllele.insert(0, ".");
                                            }
                                        } else {
                                            variantStart = alleleStart - 1;
                                            parentAllele.insert(0, lps.get(bi).getMiddle().charAt(variantStart));
                                            childAllele.insert(0, lps.get(0).getMiddle().charAt(variantStart));
                                            childAllele.append(".");
                                        }
                                    }
                                }

                                String[] bpieces = background.split(":");
                                String back = bpieces[0];
                                String chrName = bpieces[1];
                                int aStart = alleleStart;

                                if (BACKGROUNDS.containsKey(back) && bs.contains(bpieces[1])) {
                                    String ref = bpieces[1];
                                    int start = Integer.valueOf(bpieces[2].split("-")[0]);
                                    int stop = Integer.valueOf(bpieces[2].split("-")[1]);
                                    boolean isNegative = bpieces[3].equals("-");

                                    Interval it = new Interval(ref, start, stop, isNegative, null);

                                    List<SAMRecord> srs = BACKGROUNDS.get(back).align(lps.get(bi).getMiddle().replaceAll("[- ]", ""));
                                    SAMRecord sr = null;
                                    for (SAMRecord samRecord : srs) {
                                        Interval srit = new Interval(samRecord.getReferenceName(), samRecord.getAlignmentStart(), samRecord.getAlignmentEnd(), samRecord.getReadNegativeStrandFlag(), null);

                                        if (it.intersects(srit)) {
                                            sr = samRecord;
                                        }
                                    }

                                    if (sr != null) {
                                        Map<Integer, Pair<Character, Interval>> coordMapping = new HashMap<>();
                                        SmithWaterman sw = new SmithWaterman();
                                        String[] as = sw.getAlignment(lps.get(bi).getMiddle().replaceAll("[- ]", "").toUpperCase(), sr.getReadString());

                                        List<CigarElement> ces = sr.getCigar().getCigarElements();
                                        int ce = 0, cl = 0, q = 0, r = 0, ap = sr.getReadNegativeStrandFlag() ? sr.getAlignmentEnd() + 1: sr.getAlignmentStart() + 1;
                                        for (q = 0; q < as[0].length(); q++, cl++, r++) {
                                            if (cl == ces.get(ce).getLength()) {
                                                ce++;
                                                cl = 0;
                                            }

                                            while (r < lps.get(bi).getMiddle().length() && (lps.get(bi).getMiddle().charAt(r) == ' ' || lps.get(bi).getMiddle().charAt(r) == '-')) {
                                                r++;
                                            }

                                            if (as[1].charAt(q) != '-' || ces.get(ce).getOperator().equals(CigarOperator.M) || ces.get(ce).getOperator().equals(CigarOperator.X)) {
                                                coordMapping.put(r, Pair.create(as[1].charAt(q), new Interval(sr.getReferenceName(), ap, ap)));

                                                if (sr.getReadNegativeStrandFlag()) {
                                                    ap--;
                                                } else {
                                                    ap++;
                                                }
                                            }
                                        }

                                        //chrName = coordMapping.get(variantStart).getSecond().getContig();
                                        //aStart = coordMapping.get(variantStart).getSecond().getStart();
                                        chrName = sr.getReferenceName();
                                        if (coordMapping.getOrDefault(variantStart, null) != null) {
                                            aStart = coordMapping.get(variantStart).getSecond().getStart();
                                        } else {
                                            for (int w = 1; w <= variantStart; w++) {
                                                if (coordMapping.containsKey(variantStart - w)) {
                                                    aStart = coordMapping.get(variantStart - w).getSecond().getStart();
                                                    break;
                                                }

                                                if (coordMapping.containsKey(variantStart + w)) {
                                                    aStart = coordMapping.get(variantStart + w).getSecond().getStart();
                                                    break;
                                                }
                                            }
                                        }
                                    }
                                }

                                String refAllele = positiveStrand ? parentAllele.toString() : SequenceUtils.reverseComplement(parentAllele.toString());
                                String altAllele = positiveStrand ? childAllele.toString() : SequenceUtils.reverseComplement(childAllele.toString());

                                VariantContextBuilder vcb = new VariantContextBuilder()
                                        .chr(chrName)
                                        .start(aStart)
                                        .alleles(refAllele, altAllele)
                                        .attribute("PARTITION", partitionIndex)
                                        .attribute("SECTION", sectionIndex)
                                        .attribute("NOVELS", Joiner.on(",").join(cks))
                                        .noGenotypes();

                                if (isAtEnd || isAtStart) {
                                    vcb.attribute("SVTYPE", "BND");
                                    vcb.stop(aStart);
                                } else {
                                    vcb.computeEndFromAlleles(Arrays.asList(Allele.create(refAllele, true), Allele.create(altAllele)), aStart);
                                }

                                if (isExactRepeat) {
                                    vcb.attribute("IS_EXACT_REPEAT", isExactRepeat);
                                }

                                vcs.add(vcb.make());

                                background = "";
                                bi = -1;
                                alleleStart = -1;
                                alleleEnd = -1;
                                childAllele = new StringBuilder();
                                parentAllele = new StringBuilder();
                                isAtEnd = false;
                                isAtStart = false;
                                //isRepeat = true;
                                cks = new ArrayList<>();
                            }
                        }
                    }
                }
            }

            if (alleleStart >= 0) {
                boolean positiveStrand = BACKGROUNDS == null || background.contains(":+:");
                int variantStart = positiveStrand ? alleleStart : alleleEnd;
                if (parentAllele.length() != 1 || childAllele.length() != 1) {
                    variantStart = positiveStrand ? alleleStart - 1 : alleleEnd + 1;
                }

                if (parentAllele.length() != 1 || childAllele.length() != 1) {
                    if (positiveStrand) {
                        if (variantStart >= 0 && variantStart < lps.get(bi).getMiddle().length()) {
                            parentAllele.insert(0, lps.get(bi).getMiddle().charAt(variantStart));
                            childAllele.insert(0, lps.get(0).getMiddle().charAt(variantStart));
                            if (isAtEnd) {
                                childAllele.append(".");
                            }
                        } else {
                            variantStart = alleleEnd + 1;
                            parentAllele.append(lps.get(bi).getMiddle().charAt(variantStart));
                            childAllele.append(lps.get(0).getMiddle().charAt(variantStart));
                            childAllele.insert(0, ".");
                        }
                    } else {
                        if (variantStart >= 0 && variantStart < lps.get(bi).getMiddle().length()) {
                            parentAllele.append(lps.get(bi).getMiddle().charAt(variantStart));
                            childAllele.append(lps.get(0).getMiddle().charAt(variantStart));
                            if (isAtStart) {
                                childAllele.insert(0, ".");
                            }
                        } else {
                            variantStart = alleleStart - 1;
                            parentAllele.insert(0, lps.get(bi).getMiddle().charAt(variantStart));
                            childAllele.insert(0, lps.get(0).getMiddle().charAt(variantStart));
                            childAllele.append(".");
                        }
                    }
                }

                String[] bpieces = background.split(":");
                String back = bpieces[0];
                String chrName = bpieces[1];
                int aStart = alleleStart;

                if (BACKGROUNDS.containsKey(back) && bs.contains(bpieces[1])) {
                    String ref = bpieces[1];
                    int start = Integer.valueOf(bpieces[2].split("-")[0]);
                    int stop = Integer.valueOf(bpieces[2].split("-")[1]);
                    boolean isNegative = bpieces[3].equals("-");

                    Interval it = new Interval(ref, start, stop, isNegative, null);

                    List<SAMRecord> srs = BACKGROUNDS.get(back).align(lps.get(bi).getMiddle().replaceAll("[- ]", ""));
                    SAMRecord sr = null;
                    for (SAMRecord samRecord : srs) {
                        Interval srit = new Interval(samRecord.getReferenceName(), samRecord.getAlignmentStart(), samRecord.getAlignmentEnd(), samRecord.getReadNegativeStrandFlag(), null);

                        if (it.intersects(srit)) {
                            sr = samRecord;
                        }
                    }

                    if (sr != null) {
                        Map<Integer, Pair<Character, Interval>> coordMapping = new HashMap<>();
                        SmithWaterman sw = new SmithWaterman();
                        String[] as = sw.getAlignment(lps.get(bi).getMiddle().replaceAll("[- ]", "").toUpperCase(), sr.getReadString());

                        List<CigarElement> ces = sr.getCigar().getCigarElements();
                        int ce = 0, cl = 0, q = 0, r = 0, ap = sr.getReadNegativeStrandFlag() ? sr.getAlignmentEnd() + 1: sr.getAlignmentStart() + 1;
                        for (q = 0; q < as[0].length(); q++, cl++, r++) {
                            if (cl == ces.get(ce).getLength()) {
                                ce++;
                                cl = 0;
                            }

                            while (r < lps.get(bi).getMiddle().length() && (lps.get(bi).getMiddle().charAt(r) == ' ' || lps.get(bi).getMiddle().charAt(r) == '-')) {
                                r++;
                            }

                            if (as[1].charAt(q) != '-' || ces.get(ce).getOperator().equals(CigarOperator.M) || ces.get(ce).getOperator().equals(CigarOperator.X)) {
                                coordMapping.put(r, Pair.create(as[1].charAt(q), new Interval(sr.getReferenceName(), ap, ap)));

                                if (sr.getReadNegativeStrandFlag()) {
                                    ap--;
                                } else {
                                    ap++;
                                }
                            }
                        }

                        //chrName = coordMapping.get(variantStart).getSecond().getContig();
                        //aStart = coordMapping.get(variantStart).getSecond().getStart();
                        chrName = sr.getReferenceName();
                        if (coordMapping.getOrDefault(variantStart, null) != null) {
                            aStart = coordMapping.get(variantStart).getSecond().getStart();
                        } else {
                            for (int w = 1; w <= variantStart; w++) {
                                if (coordMapping.containsKey(variantStart - w)) {
                                    aStart = coordMapping.get(variantStart - w).getSecond().getStart();
                                    break;
                                }

                                if (coordMapping.containsKey(variantStart + w)) {
                                    aStart = coordMapping.get(variantStart + w).getSecond().getStart();
                                    break;
                                }
                            }
                        }
                    }
                }

                String refAllele = positiveStrand ? parentAllele.toString() : SequenceUtils.reverseComplement(parentAllele.toString());
                String altAllele = positiveStrand ? childAllele.toString() : SequenceUtils.reverseComplement(childAllele.toString());

                VariantContextBuilder vcb = new VariantContextBuilder()
                        .chr(chrName)
                        .start(aStart)
                        .alleles(refAllele, altAllele)
                        .attribute("PARTITION", partitionIndex)
                        .attribute("SECTION", sectionIndex)
                        .attribute("NOVELS", Joiner.on(",").join(cks))
                        .noGenotypes();

                if (isAtEnd || isAtStart) {
                    vcb.attribute("SVTYPE", "BND");
                    vcb.stop(aStart);
                } else {
                    vcb.computeEndFromAlleles(Arrays.asList(Allele.create(refAllele, true), Allele.create(altAllele)), aStart);
                }

                if (isExactRepeat) {
                    vcb.attribute("IS_EXACT_REPEAT", isExactRepeat);
                }

                vcs.add(vcb.make());
            }

            cks = new ArrayList<>();
            for (int j = leftLimit; j < rightLimit; j++) {
                String novelRegion = lps.get(0).getMiddle().substring(leftLimit, rightLimit).replaceAll("[- ]", "");
                for (int i = 0; i <= novelRegion.length() - GRAPH.getKmerSize(); i++) {
                    CanonicalKmer ck = new CanonicalKmer(novelRegion.substring(i, i + GRAPH.getKmerSize()));
                    if (rois.contains(ck)) {
                        cks.add(ck);
                    }
                }

                int stateBefore = 0;
                int posBefore = 0;
                StringBuilder befSequence = new StringBuilder();
                for (int i = 1; i < lps.size(); i++) {
                    if (j < lps.get(i).getMiddle().length() && lps.get(i).getMiddle().charAt(j) != ' ' && j + 1 >= lps.get(i).getMiddle().length()) {
                        stateBefore = i;
                        posBefore = j;

                        for (int k = j; k >= leftLimit; k--) {
                            befSequence.insert(0, lps.get(0).getMiddle().charAt(k));
                            posBefore = k;

                            if (Character.toUpperCase(lps.get(i).getMiddle().charAt(k)) == Character.toUpperCase(lps.get(0).getMiddle().charAt(k))) {
                                break;
                            }
                        }

                        break;
                    }
                }

                int stateAfter = 0;
                int posAfter = 0;
                StringBuilder aftSequence = new StringBuilder();
                for (int i = 1; i < lps.size(); i++) {
                    if (j + 1 < lps.get(i).getMiddle().length() && lps.get(i).getMiddle().charAt(j) == ' ' && lps.get(i).getMiddle().charAt(j + 1) != ' ') {
                        stateAfter = i;
                        posAfter = j + 1;

                        for (int k = j + 1; k < rightLimit; k--) {
                            aftSequence.append(lps.get(0).getMiddle().charAt(k));
                            posAfter = k;

                            if (Character.toUpperCase(lps.get(i).getMiddle().charAt(k)) == Character.toUpperCase(lps.get(0).getMiddle().charAt(k))) {
                                break;
                            }
                        }

                        break;
                    }
                }

                if (stateBefore != 0 && stateAfter != 0) {
                    String[] pieces0 = lps.get(stateBefore).getLeft().split("[:-]");
                    String[] pieces1 = lps.get(stateAfter).getLeft().split("[:-]");

                    String id0 = String.format("bnd_p%d_s%d_j%d", partitionIndex, sectionIndex, posBefore);
                    String id1 = String.format("bnd_p%d_s%d_j%d", partitionIndex, sectionIndex, posAfter);

                    List<Allele> alleles0;
                    if (pieces1.length > 3) {
                        alleles0 = Arrays.asList(Allele.create((byte) befSequence.charAt(0), true), Allele.create(befSequence.append("]").append(pieces1[1]).append(":").append(pieces1[2]).append("]").toString()));
                    } else {
                        alleles0 = Arrays.asList(Allele.create((byte) befSequence.charAt(0), true), Allele.create(befSequence.append("]").append(pieces1[0]).append(":").append(posBefore).append("]").toString()));
                    }

                    VariantContext vc0 = new VariantContextBuilder()
                            .chr(pieces0.length > 3 ? pieces0[1] : (pieces0[0] + "_unknown"))
                            .start(pieces0.length > 3 ? Integer.valueOf(pieces0[3]) : partitionIndex)
                            .stop(pieces0.length > 3 ? Integer.valueOf(pieces0[3]) : partitionIndex)
                            .alleles(alleles0)
                            .attribute("PARTITION", partitionIndex)
                            .attribute("SECTION", sectionIndex)
                            .attribute("NOVELS", Joiner.on(",").join(cks))
                            .attribute("SVTYPE", "BND")
                            .attribute("MATEID", id1)
                            .attribute("STATE", lps.get(stateBefore).getLeft())
                            .id(id0)
                            .noGenotypes()
                            .make();

                    List<Allele> alleles1;
                    if (pieces0.length > 3) {
                        alleles1 = Arrays.asList(Allele.create((byte) aftSequence.charAt(0), true), Allele.create(aftSequence.append("]").append(pieces0[1]).append(":").append(pieces0[3]).append("]").toString()));
                    } else {
                        alleles1 = Arrays.asList(Allele.create((byte) aftSequence.charAt(0), true), Allele.create(aftSequence.append("]").append(pieces0[0]).append(":").append(posAfter).append("]").toString()));
                    }

                    VariantContext vc1 = new VariantContextBuilder()
                            .chr(pieces1.length > 3 ? pieces1[1] : (pieces1[0] + "_unknown"))
                            .start(pieces1.length > 3 ? Integer.valueOf(pieces1[2]) : partitionIndex)
                            .stop(pieces1.length > 3 ? Integer.valueOf(pieces1[2]) : partitionIndex)
                            .alleles(alleles1)
                            .attribute("PARTITION", partitionIndex)
                            .attribute("SECTION", sectionIndex)
                            .attribute("NOVELS", Joiner.on(",").join(cks))
                            .attribute("SVTYPE", "BND")
                            .attribute("MATEID", id0)
                            .attribute("STATE", lps.get(stateAfter).getLeft())
                            .id(id1)
                            .noGenotypes()
                            .make();

                    vcs.add(vc0);
                    vcs.add(vc1);
                }
            }
        }

        return vcs;
    }

    private SAMRecord getBestAlignment(String query, String background) {
        List<SAMRecord> b = BACKGROUNDS.get(background).align(query.replaceAll("[- ]", ""));

        List<SAMRecord> a = new ArrayList<>();
        for (SAMRecord sr : b) {
            if (sr.getMappingQuality() > 0) {
                a.add(sr);
            }
        }

        a.sort((s1, s2) -> {
            int s1length = s1.getAlignmentEnd() - s1.getAlignmentStart();
            int nm1 = s1.getIntegerAttribute("NM");
            int mq1 = s1.getMappingQuality();

            int s2length = s2.getAlignmentEnd() - s2.getAlignmentStart();
            int nm2 = s2.getIntegerAttribute("NM");
            int mq2 = s1.getMappingQuality();

            if (s1length != s2length) {
                return s1length > s2length ? -1 : 1;
            }

            if (nm1 != nm2) {
                return nm1 < nm2 ? -1 : 1;
            }

            if (mq1 != mq2) {
                return mq1 > mq2 ? -1 : 1;
            }

            return 0;
        });

        return a.size() > 0 ? a.get(0) : null;
    }

    private String trimQuery(List<CortexVertex> ws, Map<String, String> targets, Set<CanonicalKmer> rois) {
        int firstIndex = Integer.MAX_VALUE, lastIndex = 0;
        int firstNovel = -1, lastNovel = -1;

        Map<CanonicalKmer, Integer> pos = new HashMap<>();
        for (int i = 0; i < ws.size(); i++) {
            pos.put(ws.get(i).getCanonicalKmer(), i);

            if (rois.contains(ws.get(i).getCanonicalKmer())) {
                if (firstNovel == -1) {
                    firstNovel = i;
                }

                lastNovel = i;
            }
        }

        for (String target : targets.values()) {
            for (int i = 0; i <= target.length() - GRAPH.getKmerSize(); i++) {
                String sk = target.substring(i, i + GRAPH.getKmerSize());
                CanonicalKmer ck = new CanonicalKmer(sk);

                if (pos.containsKey(ck)) {
                    int index = pos.get(ck);

                    if (index < firstIndex) { firstIndex = index; }
                    if (index > lastIndex) { lastIndex = index; }

                }
            }
        }

        if (firstNovel < firstIndex) { firstIndex = firstNovel; }
        if (lastNovel > lastIndex) { lastIndex = lastNovel; }

        return TraversalUtils.toContig(ws.subList(firstIndex, lastIndex));
    }

    private String makeNoveltyTrack(Set<CanonicalKmer> rois, String query, List<Triple<String, String, Pair<Integer, Integer>>> lps) {
        int maxLength = 0;
        for (Triple<String, String, Pair<Integer, Integer>> lp : lps) {
            String name = String.format("%s (%d-%d)", lp.getLeft(), lp.getRight().getFirst(), lp.getRight().getSecond());
            maxLength = Math.max(maxLength, name.length());
        }

        StringBuilder sb = new StringBuilder(StringUtil.repeatCharNTimes(' ', query.length() + 1));
        for (int i = 0; i <= query.length() - GRAPH.getKmerSize(); i++) {
            CanonicalKmer ck = new CanonicalKmer(query.substring(i, i + GRAPH.getKmerSize()));

            if (rois.contains(ck)) {
                for (int j = i; j <= i + GRAPH.getKmerSize(); j++) {
                    sb.setCharAt(j, '*');
                }
            }
        }

        for (int i = 0; i < lps.get(0).getMiddle().length(); i++) {
            if (lps.get(0).getMiddle().charAt(i) == '-') {
                sb.insert(i, sb.charAt(i) == '*' ? '*' : ' ');
            }
        }

        return String.format("%" + maxLength + "s %s", "novel", sb.toString());
    }

    private List<Pair<Integer, Integer>> getNoveltyRegions(Set<CanonicalKmer> rois, String query, List<Triple<String, String, Pair<Integer, Integer>>> lps) {
        String noveltyTrack = makeNoveltyTrack(rois, query, lps);
        noveltyTrack = noveltyTrack.replaceAll("^\\s+novel ", "");

        List<Pair<Integer, Integer>> regions = new ArrayList<>();
        int start = -1;
        int stop = noveltyTrack.length() - 1;
        for (int i = 0; i < noveltyTrack.length(); i++) {
            if (noveltyTrack.charAt(i) == '*') {
                if (start == -1) { start = i; }
                stop = i;
            } else {
                if (start >= 0) {
                    regions.add(Pair.create(start, stop));
                    start = -1;
                    stop = noveltyTrack.length() - 1;
                }
            }
        }

        if (start >= 0) {
            regions.add(Pair.create(start, stop));
        }

        return regions;
    }

    @NotNull
    private Map<String, String> assembleCandidateHaplotypes(Set<CanonicalKmer> rois, List<CortexVertex> ws, List<Pair<Integer, Integer>> regions, Set<String> parentName) {
        Map<String, String> targets = new HashMap<>();

        DirectedWeightedPseudograph<CortexVertex, CortexEdge> g = new DirectedWeightedPseudograph<>(CortexEdge.class);

        List<List<CortexVertex>> walks = new ArrayList<>();
        List<String> contigs = new ArrayList<>();

        // First build the non-novel stretches
        for (int i = 0; i < regions.size(); i++) {
            int lastEnd = i == 0 ? 0 : regions.get(i-1).getSecond();
            int nextStart = i == regions.size() - 1 ? ws.size() - 1 : regions.get(i+1).getFirst();

            DirectedWeightedPseudograph<CortexVertex, CortexEdge> gl = assembleLeft(rois, parentName, ws, regions.get(i), lastEnd);
            DirectedWeightedPseudograph<CortexVertex, CortexEdge> gr = assembleRight(rois, parentName, ws, regions.get(i), nextStart);

            Graphs.addGraph(g, gl);
            Graphs.addGraph(g, gr);

            List<CortexVertex> wl = gl.vertexSet().size() == 0 ? new ArrayList<>() : TraversalUtils.toWalk(gl, gl.vertexSet().iterator().next().getKmerAsString(), gl.edgeSet().iterator().next().getColor());
            List<CortexVertex> wr = gr.vertexSet().size() == 0 ? new ArrayList<>() : TraversalUtils.toWalk(gr, gr.vertexSet().iterator().next().getKmerAsString(), gr.edgeSet().iterator().next().getColor());

            if (wl.size() > 0 && !contigs.contains(TraversalUtils.toContig(wl))) {
                walks.add(wl);
                contigs.add(TraversalUtils.toContig(wl));
            }

            if (wr.size() > 0 && !contigs.contains(TraversalUtils.toContig(wr))) {
                walks.add(wr);
                contigs.add(TraversalUtils.toContig(wr));
            }
        }

        targets.putAll(assembleGapClosedCandidates(parentName, g, walks, ws));

        if (targets.size() == 0) {
            targets.putAll(assembleDovetailCandidates(parentName, g, walks));
        }

        if (targets.size() == 0) {
            for (int i = 0; i < walks.size(); i++) {
                List<CortexVertex> w = walks.get(i);

                String cc = TraversalUtils.toContig(w);
                String id = String.format("%s:%s_unknown:contig%d_%s", parentName.iterator().next(), parentName.iterator().next(), i, "gapped");

                targets.put(id, cc);
            }
        }

        return targets;
    }

    private Map<String, String> assembleDovetailCandidates(Set<String> parentName, DirectedWeightedPseudograph<CortexVertex, CortexEdge> g, List<List<CortexVertex>> walks) {
        // Now try dovetail overlaps with remaining gaps
        DirectedWeightedPseudograph<CortexVertex, CortexEdge> gm = new DirectedWeightedPseudograph<>(CortexEdge.class);
        Graphs.addGraph(gm, g);

        int numDovetailOverlaps = 0;

        for (int i = 0; i < walks.size() - 1; i++) {
            String cb = TraversalUtils.toContig(walks.get(i));
            String ca = TraversalUtils.toContig(walks.get(i+1));

            int longestMatch = 0;
            for (int j = 1; j < cb.length() && j < ca.length(); j++) {
                if (cb.substring(cb.length() - j, cb.length()).equals(ca.substring(0, j))) {
                    longestMatch = j;
                }
            }

            if (longestMatch > 11) {
                String cc = cb + ca.substring(longestMatch, ca.length());

                CortexVertex v0 = null;
                for (int j = 0; j <= cc.length() - GRAPH.getKmerSize(); j++) {
                    String sk = cc.substring(j, j + GRAPH.getKmerSize());
                    CortexRecord cr = GRAPH.findRecord(sk);

                    CortexVertex v = new CortexVertexFactory()
                            .bases(sk)
                            .record(cr)
                            .make();

                    gm.addVertex(v);

                    if (v0 != null) {
                        gm.addEdge(v0, v, new CortexEdge(v0, v, GRAPH.getColorsForSampleNames(parentName).iterator().next(), 1.0f));
                    }

                    v0 = v;
                }

                numDovetailOverlaps++;
            }
        }

        Map<String, String> targets = new HashMap<>();
        if (numDovetailOverlaps > 0) {
            targets.putAll(extractTargets(parentName, gm, "dovetail"));
        }

        return targets;
    }

    private Map<String, String> assembleGapClosedCandidates(Set<String> parentName, DirectedWeightedPseudograph<CortexVertex, CortexEdge> g, List<List<CortexVertex>> walks, List<CortexVertex> ws) {
        DirectedGraph<CortexByteKmer, DefaultEdge> gw = new DefaultDirectedGraph<>(DefaultEdge.class);

        gw.addVertex(ws.get(0).getKmerAsByteKmer());
        for (int i = 1; i < ws.size(); i++) {
            gw.addVertex(ws.get(i).getKmerAsByteKmer());
            gw.addEdge(ws.get(i-1).getKmerAsByteKmer(), ws.get(i).getKmerAsByteKmer());
        }

        // Now bridge any gaps that we can
        DirectedWeightedPseudograph<CortexVertex, CortexEdge> gm = new DirectedWeightedPseudograph<>(CortexEdge.class);
        Graphs.addGraph(gm, g);

        int numGapsFilled = 0;

        for (int i = 0; i < walks.size() - 1; i++) {
            DirectedWeightedPseudograph<CortexVertex, CortexEdge> gg = new DirectedWeightedPseudograph<>(CortexEdge.class);
            for (int j = i + 1; j < walks.size() && gg.vertexSet().size() == 0; j++) {
                gg = assembleMiddle(parentName, walks.get(i), walks.get(j), gw);
            }

            if (gg.vertexSet().size() > 0) {
                Graphs.addGraph(gm, gg);
                numGapsFilled++;
            } else {
                TraversalEngine ef = new TraversalEngineFactory()
                        .traversalColors(GRAPH.getColorsForSampleNames(parentName))
                        .traversalDirection(FORWARD)
                        .combinationOperator(OR)
                        .stoppingRule(ContigStopper.class)
                        .maxBranchLength(500)
                        .graph(GRAPH)
                        .links(LINKS)
                        .make();

                TraversalEngine er = new TraversalEngineFactory()
                        .traversalColors(GRAPH.getColorsForSampleNames(parentName))
                        .traversalDirection(REVERSE)
                        .combinationOperator(OR)
                        .stoppingRule(ContigStopper.class)
                        .maxBranchLength(500)
                        .graph(GRAPH)
                        .links(LINKS)
                        .make();

                if (i < walks.size() - 1) {
                    Graphs.addGraph(gm, ef.dfs(walks.get(i).get(walks.get(i).size() - 1).getKmerAsString()));
                }

                if (i > 0) {
                    Graphs.addGraph(gm, er.dfs(walks.get(i).get(walks.get(i).size() - 1).getKmerAsString()));
                }
            }
        }

        Map<String, String> targets = new HashMap<>();
        if (numGapsFilled > 0) {
            targets.putAll(extractTargets(parentName, gm, "gapfilled"));
        }

        return targets;
    }

    private Map<String, String> extractTargets(Set<String> parentName, DirectedWeightedPseudograph<CortexVertex, CortexEdge> g, String suffix) {
        Map<String, String> targets = new HashMap<>();
        ConnectivityInspector<CortexVertex, CortexEdge> ci = new ConnectivityInspector<>(g);
        int targetIndex = 0;
        for (Set<CortexVertex> cs : ci.connectedSets()) {
            List<CortexVertex> cl = TraversalUtils.toWalk(g, cs.iterator().next().getKmerAsString(), g.edgeSet().iterator().next().getColor());

            String cc = TraversalUtils.toContig(cl);
            String id = String.format("%s:%s_unknown:contig%d_%s", parentName.iterator().next(), parentName.iterator().next(), targetIndex, suffix);

            if (cc.length() > 0) {
                targets.put(id, cc);
                targetIndex++;
            }
        }

        return targets;
    }

    private DirectedWeightedPseudograph<CortexVertex, CortexEdge> assembleMiddle(Set<String> parentName, List<CortexVertex> wLeft, List<CortexVertex> wRight, DirectedGraph<CortexByteKmer, DefaultEdge> gw) {
        int childColor = GRAPH.getColorForSampleName(ROIS.getSampleName(0));

        for (int i = wLeft.size() - 1; i > 0; i--) {
            Map<Integer, Set<CortexByteKmer>> cbkOut = TraversalUtils.getAllNextKmers(wLeft.get(i).getCortexRecord(), !wLeft.get(i).getKmerAsString().equals(wLeft.get(i).getCanonicalKmer().getKmerAsString()));
            Set<CortexByteKmer> outgoingEdges = new HashSet<>();
            for (int c : GRAPH.getColorsForSampleNames(parentName)) {
                outgoingEdges.addAll(cbkOut.get(c));
            }

            if (gw.containsVertex(wLeft.get(i).getKmerAsByteKmer())) {
                outgoingEdges.removeAll(Graphs.successorListOf(gw, wLeft.get(i).getKmerAsByteKmer()));
            } else {
                outgoingEdges.removeAll(cbkOut.get(childColor));
            }

            if (outgoingEdges.size() > 0) {
                for (int j = 0; j < wRight.size() - 1; j++) {
                    Map<Integer, Set<CortexByteKmer>> cbkIn = TraversalUtils.getAllPrevKmers(wRight.get(j).getCortexRecord(), !wRight.get(j).getKmerAsString().equals(wRight.get(j).getCanonicalKmer().getKmerAsString()));
                    Set<CortexByteKmer> incomingEdges = new HashSet<>();
                    for (int c : GRAPH.getColorsForSampleNames(parentName)) {
                        incomingEdges.addAll(cbkIn.get(c));
                    }

                    if (gw.containsVertex(wRight.get(j).getKmerAsByteKmer())) {
                        incomingEdges.removeAll(Graphs.predecessorListOf(gw, wRight.get(j).getKmerAsByteKmer()));
                    } else {
                        incomingEdges.removeAll(cbkIn.get(childColor));
                    }

                    if (incomingEdges.size() > 0) {
                        TraversalEngine ef = new TraversalEngineFactory()
                                .traversalColors(GRAPH.getColorsForSampleNames(parentName))
                                .traversalDirection(FORWARD)
                                .combinationOperator(OR)
                                .stoppingRule(DestinationStopper.class)
                                .maxBranchLength(5000)
                                .graph(GRAPH)
                                .links(LINKS)
                                .make();

                        TraversalEngine er = new TraversalEngineFactory()
                                .traversalColors(GRAPH.getColorsForSampleNames(parentName))
                                .traversalDirection(REVERSE)
                                .combinationOperator(OR)
                                .stoppingRule(DestinationStopper.class)
                                .maxBranchLength(5000)
                                .graph(GRAPH)
                                .links(LINKS)
                                .make();

                        DirectedWeightedPseudograph<CortexVertex, CortexEdge> g = ef.dfs(wLeft.get(i-1).getKmerAsString(), wRight.get(j+1).getKmerAsString());
                        if (g == null || g.vertexSet().size() == 0) {
                            g = er.dfs(wRight.get(j+1).getKmerAsString(), wLeft.get(i-1).getKmerAsString());
                        }

                        if (g != null && g.vertexSet().size() > 0) {
                            return g;
                        }
                    }
                }
            }
        }

        return new DirectedWeightedPseudograph<>(CortexEdge.class);
    }

    private DirectedWeightedPseudograph<CortexVertex, CortexEdge> assembleLeft(Set<CanonicalKmer> rois, Set<String> parentName, List<CortexVertex> ws, Pair<Integer, Integer> region, int lastEnd) {
        DirectedWeightedPseudograph<CortexVertex, CortexEdge> gl = new DirectedWeightedPseudograph<>(CortexEdge.class);

        for (int j = region.getFirst() - 1; j > lastEnd && !rois.contains(ws.get(j).getCanonicalKmer()); j--) {
            boolean hasCoverage = false;
            for (int c : GRAPH.getColorsForSampleNames(parentName)) {
                if (ws.get(j).getCortexRecord().getCoverage(c) > 0) {
                    hasCoverage = true;
                    break;
                }
            }

            if (hasCoverage) {
                TraversalEngine e = new TraversalEngineFactory()
                        .traversalColors(GRAPH.getColorsForSampleNames(parentName))
                        .traversalDirection(REVERSE)
                        .combinationOperator(OR)
                        .stoppingRule(ContigStopper.class)
                        .maxBranchLength((region.getFirst() - 1) - (lastEnd + 1))
                        .graph(GRAPH)
                        .links(LINKS)
                        .make();

                DirectedWeightedPseudograph<CortexVertex, CortexEdge> g = e.dfs(ws.get(j).getKmerAsString());

                if (g != null) {
                    if (g.vertexSet().size() > gl.vertexSet().size()) {
                        gl = g;
                    } else if (gl.vertexSet().size() > 0 && g.vertexSet().size() <= gl.vertexSet().size()) {
                        break;
                    }
                }
            }
        }

        return gl;
    }

    private DirectedWeightedPseudograph<CortexVertex, CortexEdge> assembleRight(Set<CanonicalKmer> rois, Set<String> parentName, List<CortexVertex> ws, Pair<Integer, Integer> region, int nextStart) {
        DirectedWeightedPseudograph<CortexVertex, CortexEdge> gr = new DirectedWeightedPseudograph<>(CortexEdge.class);

        for (int j = region.getSecond() + 1; j < nextStart && !rois.contains(ws.get(j).getCanonicalKmer()); j++) {
            boolean hasCoverage = false;
            for (int c : GRAPH.getColorsForSampleNames(parentName)) {
                if (ws.get(j).getCortexRecord().getCoverage(c) > 0) {
                    hasCoverage = true;
                    break;
                }
            }

            if (hasCoverage) {
                TraversalEngine e = new TraversalEngineFactory()
                        .traversalColors(GRAPH.getColorsForSampleNames(parentName))
                        .traversalDirection(FORWARD)
                        .combinationOperator(OR)
                        .stoppingRule(ContigStopper.class)
                        .maxBranchLength((nextStart - 1) - (region.getSecond() + 1))
                        .graph(GRAPH)
                        .links(LINKS)
                        .make();

                DirectedWeightedPseudograph<CortexVertex, CortexEdge> g = e.dfs(ws.get(j).getKmerAsString());

                if (g != null) {
                    if (g.vertexSet().size() > gr.vertexSet().size()) {
                        gr = g;
                    } else if (gr.vertexSet().size() > 0 && g.vertexSet().size() <= gr.vertexSet().size()) {
                        break;
                    }
                }
            }
        }

        return gr;
    }

    private Set<CanonicalKmer> loadRois(CortexGraph rg) {
        Set<CanonicalKmer> rois = new HashSet<>();

        for (CortexRecord rr : rg) {
            rois.add(rr.getCanonicalKmer());
        }

        return rois;
    }

    private List<CortexVertex> loadChildWalk(ReferenceSequence seq, CortexGraph graph) {
        List<CortexVertex> w = new ArrayList<>();

        Map<String, Integer> seenCount = new HashMap<>();

        String contig = seq.getBaseString();
        for (int i = 0; i <= contig.length() - graph.getKmerSize(); i++) {
            String sk = contig.substring(i, i + graph.getKmerSize());

            if (!seenCount.containsKey(sk)) {
                seenCount.put(sk, 0);
            } else {
                ContainerUtils.increment(seenCount, sk);
            }

            w.add(new CortexVertexFactory()
                    .bases(sk)
                    .record(graph.findRecord(sk))
                    .copyIndex(seenCount.get(sk))
                    .make());
        }

        return w;
    }

    private List<Triple<Integer, Integer, List<CortexVertex>>> sectionContig(Set<CanonicalKmer> rois, List<CortexVertex> w, int window, int novelDistanceSplit) {
        List<Pair<Integer, Integer>> regions = getRegions(rois, w);

        int subcontigStart = regions.get(0).getFirst() - window;
        if (subcontigStart < 0) { subcontigStart = 0; }

        int subcontigStop  = regions.get(regions.size() - 1).getSecond() + window;
        if (subcontigStop >= w.size()) { subcontigStop = w.size() - 1; }

        List<Pair<Integer, Integer>> sections = new ArrayList<>();
        for (int i = 0; i < regions.size() - 1; i++) {
            if (regions.get(i+1).getFirst() - regions.get(i).getSecond() > novelDistanceSplit) {
                sections.add(Pair.create(subcontigStart, regions.get(i).getSecond() + window));

                subcontigStart = regions.get(i+1).getFirst() - window;
            }
        }

        sections.add(Pair.create(subcontigStart, subcontigStop));

        List<Triple<Integer, Integer, List<CortexVertex>>> wss = new ArrayList<>();
        for (Pair<Integer, Integer> section : sections) {
            List<CortexVertex> ws = new ArrayList<>();
            for (int i = section.getFirst(); i <= section.getSecond(); i++) {
                ws.add(w.get(i));
            }

            wss.add(Triple.of(section.getFirst(), section.getSecond(), ws));
        }

        return wss;
    }

    @NotNull
    private List<Pair<Integer, Integer>> getRegions(Set<CanonicalKmer> rois, List<CortexVertex> cvs) {
        List<Pair<Integer, Integer>> regions = new ArrayList<>();

        int regionStart = -1;
        int regionStop = 0;
        for (int i = 0; i < cvs.size(); i++) {
            CanonicalKmer currentKmer = cvs.get(i).getCanonicalKmer();

            if (rois.contains(currentKmer)) {
                if (regionStart == -1) { regionStart = i; }
                regionStop = i;
            } else {
                if (regionStart > -1) {
                    regions.add(new Pair<>(regionStart, regionStop));

                    regionStart = -1;
                    regionStop = 0;
                }
            }
        }

        if (regionStart > -1) {
            regions.add(new Pair<>(regionStart, regionStop));
        }

        return regions;
    }
}
