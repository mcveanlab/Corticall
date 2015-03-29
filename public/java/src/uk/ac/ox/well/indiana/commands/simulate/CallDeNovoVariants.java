package uk.ac.ox.well.indiana.commands.simulate;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileReader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;
import htsjdk.samtools.util.StringUtil;
import org.jgrapht.DirectedGraph;
import org.jgrapht.graph.DefaultEdge;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.table.TableReader;
import uk.ac.ox.well.indiana.utils.sequence.SequenceUtils;

import java.io.File;
import java.io.PrintStream;
import java.util.*;

public class CallDeNovoVariants extends Module {
    @Argument(fullName="metrics", shortName="m", doc="Annotated contigs")
    public File METRICS;

    @Argument(fullName="ref0", shortName="r0", doc="Reference for parent 0")
    public IndexedFastaSequenceFile REF0;

    @Argument(fullName="ref1", shortName="r1", doc="Reference for parent 1")
    public IndexedFastaSequenceFile REF1;

    @Argument(fullName="bam0", shortName="b0", doc="BAM for parent 0")
    public SAMFileReader BAM0;

    @Argument(fullName="bam1", shortName="b1", doc="BAM for parent 1")
    public SAMFileReader BAM1;

    @Output
    public File out;

    /*
    @Output(fullName="vcfOut0", shortName="vo0", doc="VCF out for variants vs parent 0")
    public File vo0;

    @Output(fullName="vcfOut1", shortName="vo1", doc="VCF out for variants vs parent 1")
    public File vo1;
    */

    private class Entry {
        private Map<String, String> attributes;
        private SAMRecord sr0;
        private SAMRecord sr1;
    }

    private Map<String, Entry> loadEntriesWithNovelKmers() {
        Map<String, Entry> entries = new LinkedHashMap<String, Entry>();

        TableReader tr = new TableReader(METRICS);

        for (Map<String, String> te : tr) {
            int refNovel = Integer.valueOf(te.get("refNovel"));

            if (refNovel > 0) {
                String contigName = te.get("contigName");

                Entry e = new Entry();
                e.attributes = te;

                entries.put(contigName, e);
            }
        }

        for (SAMRecord sr0 : BAM0) {
            if (entries.containsKey(sr0.getReadName())) {
                entries.get(sr0.getReadName()).sr0 = sr0;
            }
        }

        for (SAMRecord sr1 : BAM1) {
            if (entries.containsKey(sr1.getReadName())) {
                entries.get(sr1.getReadName()).sr1 = sr1;
            }
        }

        return entries;
    }

    private class AnnotatedKmer {
        private String sk;
        private boolean isNovel;

        public AnnotatedKmer(String sk, boolean isNovel) {
            this.sk = sk;
            this.isNovel = isNovel;
        }

        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;

            AnnotatedKmer kmer = (AnnotatedKmer) o;

            if (isNovel != kmer.isNovel) return false;
            if (sk != null ? !sk.equals(kmer.sk) : kmer.sk != null) return false;

            return true;
        }

        @Override
        public int hashCode() {
            int result = sk != null ? sk.hashCode() : 0;
            result = 31 * result + (isNovel ? 1 : 0);
            return result;
        }
    }

    private class AnnotatedEdge extends DefaultEdge {
        private boolean isInParent;
        private boolean isInChild;

        public AnnotatedEdge() {}
    }

    private void printGraph(DirectedGraph<AnnotatedKmer, AnnotatedEdge> g, PrintStream fout) {
        fout.println("digraph G {");

        for (AnnotatedKmer v : g.vertexSet()) {
            String extra = "";
            fout.println("    \"" + v.sk + "\"" + extra + ";");
        }

        for (AnnotatedEdge e : g.edgeSet()) {
            AnnotatedKmer vs = g.getEdgeSource(e);
            AnnotatedKmer vt = g.getEdgeTarget(e);

            //String extra = "";
            if (e.isInParent) {
                fout.println("    \"" + vs.sk + "\" -> \"" + vt.sk + "\" [color=blue];");
            }

            if (e.isInChild) {
                fout.println("    \"" + vs.sk + "\" -> \"" + vt.sk + "\" [color=red];");
            }
        }

        fout.println("}");
    }

    private enum PutativeVariantType { SNP, INSERTION, DELETION, INVERSION, STR_CON, STR_EXP, TD }

    private class PutativeVariant {
        private String refAllele = "";
        private int refPos = -1;
        private String refPrevBase = "";

        private String altAllele = "";
        private int altPos = -1;
        private String altPrevBase = "";

        private PutativeVariantType pt = PutativeVariantType.SNP;
        private int numContainedNovelKmers = 0;
        private int numFlankingNovelKmers = 0;

        private Set<String> plausibleAlleles = new HashSet<String>();

        private boolean isFilteredOut = false;

        @Override
        public String toString() {
            return "PutativeVariant{" +
                    "pt=" + pt +
                    ", refAllele='" + refAllele + '\'' +
                    ", refPos=" + refPos +
                    ", refPrevBase='" + refPrevBase + '\'' +
                    ", altAllele='" + altAllele + '\'' +
                    ", altPos=" + altPos +
                    ", altPrevBase='" + altPrevBase + '\'' +
                    ", numContainedNovelKmers=" + numContainedNovelKmers +
                    ", numFlankingNovelKmers=" + numFlankingNovelKmers +
                    '}';
        }

        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;

            PutativeVariant that = (PutativeVariant) o;

            if (altPos != that.altPos) return false;
            if (numContainedNovelKmers != that.numContainedNovelKmers) return false;
            if (numFlankingNovelKmers != that.numFlankingNovelKmers) return false;
            if (refPos != that.refPos) return false;
            if (altAllele != null ? !altAllele.equals(that.altAllele) : that.altAllele != null) return false;
            if (altPrevBase != null ? !altPrevBase.equals(that.altPrevBase) : that.altPrevBase != null) return false;
            if (plausibleAlleles != null ? !plausibleAlleles.equals(that.plausibleAlleles) : that.plausibleAlleles != null)
                return false;
            if (pt != that.pt) return false;
            if (refAllele != null ? !refAllele.equals(that.refAllele) : that.refAllele != null) return false;
            if (refPrevBase != null ? !refPrevBase.equals(that.refPrevBase) : that.refPrevBase != null) return false;

            return true;
        }

        @Override
        public int hashCode() {
            int result = refAllele != null ? refAllele.hashCode() : 0;
            result = 31 * result + refPos;
            result = 31 * result + (refPrevBase != null ? refPrevBase.hashCode() : 0);
            result = 31 * result + (altAllele != null ? altAllele.hashCode() : 0);
            result = 31 * result + altPos;
            result = 31 * result + (altPrevBase != null ? altPrevBase.hashCode() : 0);
            result = 31 * result + (pt != null ? pt.hashCode() : 0);
            result = 31 * result + numContainedNovelKmers;
            result = 31 * result + numFlankingNovelKmers;
            result = 31 * result + (plausibleAlleles != null ? plausibleAlleles.hashCode() : 0);
            return result;
        }
    }

    private boolean isComplementaryMutationType(PutativeVariant pvi, PutativeVariant pvj) {
        if (pvi.pt == PutativeVariantType.SNP && pvj.pt == PutativeVariantType.SNP) {
            return true;
        } else if (pvi.pt == PutativeVariantType.DELETION && pvj.pt == PutativeVariantType.INSERTION) {
            return true;
        } else if (pvi.pt == PutativeVariantType.INSERTION && pvj.pt == PutativeVariantType.DELETION) {
            return true;
        }

        return false;
    }

    private boolean isComplementaryAllele(PutativeVariant pvi, PutativeVariant pvj) {
        if (pvi.pt == PutativeVariantType.SNP && pvj.pt == PutativeVariantType.SNP) {
            return pvi.refAllele.equals(SequenceUtils.reverseComplement(pvj.altAllele)) && pvi.altAllele.equals(SequenceUtils.reverseComplement(pvj.refAllele));
        } else if (pvi.pt == PutativeVariantType.DELETION && pvj.pt == PutativeVariantType.INSERTION) {
            //return pvi.refAllele.equals(SequenceUtils.reverseComplement(pvj.altAllele));

            for (String delAllele : pvi.plausibleAlleles) {
                for (String insAllele : pvj.plausibleAlleles) {
                    String insAlleleRc = SequenceUtils.reverseComplement(insAllele);

                    if (delAllele.equals(insAlleleRc)) {
                        return true;
                    }
                }
            }

            return false;
        } else if (pvi.pt == PutativeVariantType.INSERTION && pvj.pt == PutativeVariantType.DELETION) {
            //return pvi.altAllele.equals(SequenceUtils.reverseComplement(pvj.refAllele));

            for (String insAllele : pvi.plausibleAlleles) {
                for (String delAllele : pvj.plausibleAlleles) {
                    String delAlleleRc = SequenceUtils.reverseComplement(delAllele);

                    if (insAllele.equals(delAlleleRc)) {
                        return true;
                    }
                }
            }

            return false;
        }

        return false;
    }

    private void showEntry(Entry e, int refIndex) {
        IndexedFastaSequenceFile ref = null;
        SAMRecord sr = null;

        if (!e.attributes.get("ref" + refIndex + "Locus").equals("NA") && !e.attributes.get("ref" + refIndex + "Locus").contains("*")) {
            ref = (refIndex == 0) ? REF0 : REF1;
            sr = (refIndex == 0) ? e.sr0 : e.sr1;
        }

        if (ref != null && sr != null && !sr.getReadUnmappedFlag()) {
            log.info("ref{}: {}", refIndex, sr.getSAMString());

            String refSequence = new String(ref.getSubsequenceAt(sr.getReferenceName(), sr.getAlignmentStart(), sr.getAlignmentEnd()).getBases());
            String contig = sr.getReadString();

            int kmerSize = e.attributes.get("seq").length() - e.attributes.get("kmerOrigin").length() + 1;
            String kmerOrigin = e.attributes.get("kmerOrigin");
            if (sr.getReadNegativeStrandFlag()) {
                kmerOrigin = SequenceUtils.reverse(kmerOrigin);
            }

            StringBuilder psb = new StringBuilder(StringUtil.repeatCharNTimes(' ', refSequence.length() + 10));
            StringBuilder rsb = new StringBuilder(refSequence);
            StringBuilder asb = new StringBuilder();
            StringBuilder csb = new StringBuilder(contig);
            StringBuilder ksb = new StringBuilder(kmerOrigin);
            ksb.append(StringUtil.repeatCharNTimes(' ', kmerSize));

            for (int i = 0; i < refSequence.length(); i += 20) {
                String pos = String.valueOf(i);
                psb.replace(i, i + pos.length(), pos);
            }

            int rpos = 0, cpos = 0;
            if (sr.getCigar().getCigarElement(0).getOperator().equals(CigarOperator.S)) {
                csb.delete(0, sr.getCigar().getCigarElement(0).getLength());
                ksb.delete(0, sr.getCigar().getCigarElement(0).getLength());

                cpos += sr.getCigar().getCigarElement(0).getLength();
            } else if (sr.getCigar().getCigarElement(0).getOperator().equals(CigarOperator.H)) {
                ksb.delete(0, sr.getCigar().getCigarElement(0).getLength());
            }

            List<PutativeVariant> pvs = new ArrayList<PutativeVariant>();
            int pos = 0;
            for (CigarElement ce : sr.getCigar().getCigarElements()) {
                if (!ce.getOperator().equals(CigarOperator.S) && !ce.getOperator().equals(CigarOperator.H)) {
                    if (ce.getOperator().equals(CigarOperator.M)) {
                        for (int i = 0; i < ce.getLength(); i++) {
                            if (rsb.charAt(pos + i) == csb.charAt(pos + i)) {
                                asb.append('|');
                            } else {
                                asb.append(' ');

                                String refAllele = String.valueOf(refSequence.charAt(rpos + i));
                                String altAllele = String.valueOf(contig.charAt(cpos + i));

                                PutativeVariant pv = new PutativeVariant();
                                pv.refAllele = refAllele;
                                pv.refPos = rpos + i;
                                pv.altAllele = altAllele;
                                pv.altPos = cpos + i;
                                pv.pt = PutativeVariantType.SNP;
                                pv.plausibleAlleles = new LinkedHashSet<String>();
                                pv.plausibleAlleles.add(altAllele);

                                pvs.add(pv);
                            }
                        }

                        rpos += ce.getLength();
                        cpos += ce.getLength();
                    } else if (ce.getOperator().equals(CigarOperator.I)) {
                        asb.append('v');
                        asb.append(StringUtil.repeatCharNTimes(' ', ce.getLength() - 1));

                        rsb.insert(pos, StringUtil.repeatCharNTimes('-', ce.getLength()));
                        psb.insert(pos, StringUtil.repeatCharNTimes(' ', ce.getLength()));

                        String altAllele = contig.substring(cpos, cpos + ce.getLength());

                        PutativeVariant pv = new PutativeVariant();
                        pv.refAllele = "";
                        pv.refPos = rpos;
                        pv.altAllele = altAllele;
                        pv.altPos = cpos;
                        pv.pt = PutativeVariantType.INSERTION;
                        pv.plausibleAlleles = new LinkedHashSet<String>();
                        pv.plausibleAlleles.add(altAllele);

                        pvs.add(pv);

                        cpos += ce.getLength();
                    } else if (ce.getOperator().equals(CigarOperator.D)) {
                        asb.append('^');
                        asb.append(StringUtil.repeatCharNTimes(' ', ce.getLength() - 1));

                        String refAllele = refSequence.substring(rpos, rpos + ce.getLength());

                        PutativeVariant pv = new PutativeVariant();
                        pv.refAllele = refAllele;
                        pv.refPos = rpos;
                        pv.altAllele = "";
                        pv.altPos = cpos;
                        pv.pt = PutativeVariantType.DELETION;
                        pv.plausibleAlleles = new LinkedHashSet<String>();
                        pv.plausibleAlleles.add(refAllele);

                        for (int i = 0; i < refAllele.length(); i++) {
                            if (rpos + i + 1 < refSequence.length() && rpos + i + 1 + ce.getLength() < refSequence.length() && cpos + i < contig.length() && refAllele.charAt(i) == contig.charAt(cpos + i)) {
                                String plausibleAllele = refSequence.substring(rpos + i + 1, rpos + i + 1 + ce.getLength());
                                pv.plausibleAlleles.add(plausibleAllele);
                            } else {
                                break;
                            }
                        }

                        pvs.add(pv);

                        csb.insert(pos, StringUtil.repeatCharNTimes('-', ce.getLength()));
                        ksb.insert(pos, StringUtil.repeatCharNTimes(' ', ce.getLength()));

                        rpos += ce.getLength();
                    }

                    pos += ce.getLength();
                }
            }

            log.info("  {}", psb.toString());
            log.info("  {}", rsb.toString());
            log.info("  {}", asb.toString());
            log.info("  {}", csb.toString());
            log.info("  {}", ksb.toString());

            IntervalTreeMap<Integer> inversions = new IntervalTreeMap<Integer>();
            for (int i = 0; i < pvs.size(); i++) {
                for (int j = pvs.size() - 1; j > i; j--) {
                    PutativeVariant pvi = pvs.get(i);
                    PutativeVariant pvj = pvs.get(j);

                    if (isComplementaryMutationType(pvi, pvj) && isComplementaryAllele(pvi, pvj)) {
                        String putativeInversion = SequenceUtils.reverseComplement(refSequence.substring(pvi.refPos + pvi.refAllele.length(), pvj.refPos));

                        if (contig.contains(putativeInversion)) {
                            Interval inversion = new Interval(e.attributes.get("contigName"), i, j);

                            if (!inversions.containsContained(inversion) && !inversions.containsOverlapping(inversion)) {
                                inversions.put(inversion, null);
                            }

                            break;
                        }
                    }
                }
            }

            if (inversions.size() > 0) {
                for (Interval inversion : inversions.keySet()) {
                    PutativeVariant pvi = pvs.get(inversion.getStart());
                    PutativeVariant pvj = pvs.get(inversion.getEnd());

                    // shrink
                    int refStart = pvi.refPos;
                    int refEnd = pvj.refPos + pvj.refAllele.length();
                    boolean found = false;
                    do {
                        String refAlleleR = SequenceUtils.reverseComplement(refSequence.substring(refStart + 1, refEnd));
                        String refAlleleL = SequenceUtils.reverseComplement(refSequence.substring(refStart, refEnd - 1));

                        if (contig.contains(refAlleleR)) {
                            refStart = refStart + 1;
                            found = true;
                        } else if (contig.contains(refAlleleL)) {
                            refEnd = refEnd - 1;
                            found = true;
                        } else {
                            refStart = refStart + 1;
                            refEnd = refEnd - 1;
                        }
                    } while (!found);

                    // grow right
                    int rightShift = 0;
                    boolean keepGoing = true;
                    do {
                        String refAlleleR = (refEnd + rightShift < refSequence.length()) ? SequenceUtils.reverseComplement(refSequence.substring(refStart, refEnd + rightShift)) : null;

                        if (refAlleleR != null && contig.contains(refAlleleR)) {
                            refEnd = refEnd + rightShift;
                            rightShift++;
                        } else {
                            keepGoing = false;
                        }
                    } while (keepGoing);

                    int leftShift = 0;
                    keepGoing = true;

                    do {
                        String refAlleleL = (refStart - leftShift >= 0) ? SequenceUtils.reverseComplement(refSequence.substring(refStart - leftShift, refEnd)) : null;

                        if (refAlleleL != null && contig.contains(refAlleleL)) {
                            refStart = refStart - leftShift;
                            leftShift++;
                        } else {
                            keepGoing = false;
                        }
                    } while (keepGoing);

                    String refAllele = refSequence.substring(refStart, refEnd);
                    String altAllele = SequenceUtils.reverseComplement(refAllele);

                    for (int i = inversion.getStart(); i <= inversion.getEnd(); i++) {
                        pvs.get(i).isFilteredOut = true;
                    }

                    PutativeVariant pv = new PutativeVariant();
                    pv.refAllele = refAllele;
                    pv.refPos = pvi.refPos;
                    pv.altAllele = altAllele;
                    pv.altPos = pvi.altPos;
                    pv.pt = PutativeVariantType.INVERSION;

                    pvs.add(pv);
                }
            }

            for (int i = 0; i < pvs.size(); i++) {
                if (!pvs.get(i).isFilteredOut && pvs.get(i).pt != PutativeVariantType.SNP && pvs.get(i).pt != PutativeVariantType.INVERSION) {
                    String allele = (pvs.get(i).pt == PutativeVariantType.DELETION) ? pvs.get(i).refAllele : pvs.get(i).altAllele;

                    String finalRepeatingUnit = "";
                    int finalRepeatUnitLength = allele.length();
                    //int finalNumRepeats = 0;

                    for (int repeatUnitLength = allele.length() - 1; repeatUnitLength >= 1; repeatUnitLength--) {
                        if (allele.length() % repeatUnitLength == 0) {
                            String repeatingUnit = allele.substring(0, repeatUnitLength);

                            int numRepeats = 0;
                            boolean repeatIsComplete = true;
                            for (int j = 0; j < allele.length(); j += repeatUnitLength) {
                                if (allele.substring(j, j + repeatUnitLength).equals(repeatingUnit)) {
                                    numRepeats++;
                                } else {
                                    repeatIsComplete = false;
                                    break;
                                }
                            }

                            if (repeatIsComplete) {
                                if (repeatUnitLength < finalRepeatUnitLength) {
                                    finalRepeatingUnit = repeatingUnit;
                                    //finalNumRepeats = numRepeats;
                                    finalRepeatUnitLength = repeatUnitLength;
                                }
                            }
                        }
                    }

                    int prevBasesStart1 = pvs.get(i).refPos - finalRepeatUnitLength;
                    int prevBasesEnd1 = pvs.get(i).refPos;
                    String prevBases1 = (prevBasesStart1 > 0) ? refSequence.substring(prevBasesStart1, prevBasesEnd1) : null;

                    int nextBasesStart1 = pvs.get(i).refPos + allele.length();
                    int nextBasesEnd1 = pvs.get(i).refPos + allele.length() + finalRepeatUnitLength;
                    String nextBases1 = (nextBasesEnd1 < refSequence.length()) ? refSequence.substring(nextBasesStart1, nextBasesEnd1) : null;

                    boolean matchesExistingStr = prevBases1 != null && nextBases1 != null && (finalRepeatingUnit.equals(prevBases1) || finalRepeatingUnit.equals(nextBases1));

                    if (!finalRepeatingUnit.isEmpty() && matchesExistingStr) {
                        pvs.get(i).pt = pvs.get(i).pt == PutativeVariantType.DELETION ? PutativeVariantType.STR_CON : PutativeVariantType.STR_EXP;
                    }
                }
            }

            for (int i = 0; i < pvs.size(); i++) {
                if (!pvs.get(i).isFilteredOut && pvs.get(i).pt == PutativeVariantType.INSERTION) {
                    String allele = pvs.get(i).altAllele;

                    int prevBasesStart1 = pvs.get(i).refPos - allele.length();
                    int prevBasesEnd1 = pvs.get(i).refPos;

                    int nextBasesStart1 = pvs.get(i).refPos;
                    int nextBasesEnd1 = pvs.get(i).refPos + allele.length();

                    String prevBases = prevBasesStart1 >= 0 ? refSequence.substring(prevBasesStart1, prevBasesEnd1) : null;
                    String nextBases = nextBasesEnd1 < refSequence.length() ? refSequence.substring(nextBasesStart1, nextBasesEnd1) : null;

                    log.info("TD: {} {} {}", allele, prevBases, nextBases);
                }
            }

            for (int i = 0; i < pvs.size(); i++) {
                if (!pvs.get(i).isFilteredOut) {
                    log.info("  variant: {}", pvs.get(i));
                }
            }

            log.info("");
        }
    }

    private void showEntry(Entry e) {
        log.info("{}:", e.attributes.get("contigName"));

        showEntry(e, 0);
        showEntry(e, 1);
    }

    @Override
    public void execute() {
        Map<String, Entry> entries = loadEntriesWithNovelKmers();

        for (String contigName : entries.keySet()) {
            showEntry(entries.get(contigName));
        }
    }
}
