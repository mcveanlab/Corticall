package uk.ac.ox.well.indiana.commands.simulate;

import com.google.common.base.Joiner;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileReader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;
import htsjdk.samtools.util.StringUtil;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.containers.DataTable;
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

    @Argument(fullName="show", shortName="s", doc="Show selected events")
    public Boolean SHOW = false;

    @Argument(fullName="knownTable", shortName="k", doc="Table of known events", required=false)
    public File KNOWN_TABLE;

    @Output
    public File out;

    private class Entry {
        private Map<String, String> attributes;
        private SAMRecord sr0;
        private SAMRecord sr1;
    }

    private Map<String, Entry> loadEntriesWithNovelKmers() {
        Map<String, Entry> entries = new LinkedHashMap<String, Entry>();

        TableReader tr = new TableReader(METRICS);

        for (Map<String, String> te : tr) {
            //int refNovel = Integer.valueOf(te.get("refNovel"));

            //if (refNovel > 0) {
                String contigName = te.get("contigName");

                Entry e = new Entry();
                e.attributes = te;

                entries.put(contigName, e);
            //}
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

    private enum PutativeVariantType { SNP, INS, DEL, INV, STR_CON, STR_EXP, TD, UNCLASSIFIED }

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

        private int distanceFromVariantLeft = Integer.MAX_VALUE;
        private int distanceFromVariantRight = Integer.MAX_VALUE;

        private SAMRecord sr;
        private StringBuilder psb;
        private StringBuilder rsb;
        private StringBuilder asb;
        private StringBuilder csb;
        private StringBuilder ksb;

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

        public void showAlignment(PrintStream out) {
            out.println(sr.getSAMString());
            out.println(psb);
            out.println(rsb);
            out.println(asb);
            out.println(csb);
            out.println(ksb);
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

    private Map<String, Set<Map<String, String>>> loadExistingEvents() {
        Map<String, Set<Map<String, String>>> knownVariants = new HashMap<String, Set<Map<String, String>>>();

        if (KNOWN_TABLE != null) {
            TableReader tr = new TableReader(KNOWN_TABLE);

            for (Map<String, String> te : tr) {
                String contig1 = te.get("contig1");

                if (!knownVariants.containsKey(contig1)) {
                    knownVariants.put(contig1, new HashSet<Map<String, String>>());
                }

                knownVariants.get(contig1).add(te);
            }
        }

        return knownVariants;
    }

    private boolean isComplementaryMutationType(PutativeVariant pvi, PutativeVariant pvj) {
        if (pvi.pt == PutativeVariantType.SNP && pvj.pt == PutativeVariantType.SNP) {
            return true;
        } else if (pvi.pt == PutativeVariantType.DEL && pvj.pt == PutativeVariantType.INS) {
            return true;
        } else if (pvi.pt == PutativeVariantType.INS && pvj.pt == PutativeVariantType.DEL) {
            return true;
        }

        return false;
    }

    private boolean isComplementaryAllele(PutativeVariant pvi, PutativeVariant pvj) {
        if (pvi.pt == PutativeVariantType.SNP && pvj.pt == PutativeVariantType.SNP) {
            return pvi.refAllele.equals(SequenceUtils.reverseComplement(pvj.altAllele)) && pvi.altAllele.equals(SequenceUtils.reverseComplement(pvj.refAllele));
        } else if (pvi.pt == PutativeVariantType.DEL && pvj.pt == PutativeVariantType.INS) {
            for (String delAllele : pvi.plausibleAlleles) {
                for (String insAllele : pvj.plausibleAlleles) {
                    String insAlleleRc = SequenceUtils.reverseComplement(insAllele);

                    if (delAllele.equals(insAlleleRc)) {
                        return true;
                    }
                }
            }

            return false;
        } else if (pvi.pt == PutativeVariantType.INS && pvj.pt == PutativeVariantType.DEL) {
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

    private Map<Integer, PutativeVariant> getVariants(Entry e, int refIndex) {
        Map<Integer, PutativeVariant> goodVariants = new TreeMap<Integer, PutativeVariant>();

        IndexedFastaSequenceFile ref = null;
        SAMRecord sr = null;

        log.debug("{}:", e.attributes.get("contigName"));

        if (!e.attributes.get("ref" + refIndex + "Locus").equals("NA") && !e.attributes.get("ref" + refIndex + "Locus").contains("*")) {
            ref = (refIndex == 0) ? REF0 : REF1;
            sr = (refIndex == 0) ? e.sr0 : e.sr1;
        }

        if (ref != null && sr != null && !sr.getReadUnmappedFlag()) {
            log.debug("  ref{}: {}", refIndex, sr.getSAMString().trim());

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

            if (sr.getCigar().getCigarElement(sr.getCigar().numCigarElements() - 1).getOperator().equals(CigarOperator.H)) {
                int hclength = sr.getCigar().getCigarElement(sr.getCigar().numCigarElements() - 1).getLength();
                ksb.delete(ksb.length() - hclength - 1, ksb.length());
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

                        //String prevBase = String.valueOf(contig.charAt(cpos - 1));
                        String altAllele = contig.substring(cpos, cpos + ce.getLength());

                        PutativeVariant pv = new PutativeVariant();
                        pv.refAllele = "";
                        pv.refPos = rpos;
                        pv.altAllele = altAllele;
                        pv.altPos = cpos;
                        pv.pt = PutativeVariantType.INS;
                        pv.plausibleAlleles = new LinkedHashSet<String>();
                        pv.plausibleAlleles.add(altAllele);

                        pvs.add(pv);

                        cpos += ce.getLength();
                    } else if (ce.getOperator().equals(CigarOperator.D)) {
                        asb.append('^');
                        asb.append(StringUtil.repeatCharNTimes(' ', ce.getLength() - 1));

                        //String prevBase = String.valueOf(refSequence.charAt(rpos - 1));
                        String refAllele = refSequence.substring(rpos, rpos + ce.getLength());

                        PutativeVariant pv = new PutativeVariant();
                        pv.refAllele = refAllele;
                        pv.refPos = rpos;
                        pv.altAllele = "";
                        pv.altPos = cpos;
                        pv.pt = PutativeVariantType.DEL;
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

            log.debug("  {}", psb.toString());
            log.debug("  {}", rsb.toString());
            log.debug("  {}", asb.toString());
            log.debug("  {}", csb.toString());
            log.debug("  {}", ksb.toString());

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
                    pv.pt = PutativeVariantType.INV;

                    pvs.add(pv);
                }
            }

            for (int i = 0; i < pvs.size(); i++) {
                if (!pvs.get(i).isFilteredOut && pvs.get(i).pt != PutativeVariantType.SNP && pvs.get(i).pt != PutativeVariantType.INV) {
                    String allele = (pvs.get(i).pt == PutativeVariantType.DEL) ? pvs.get(i).refAllele : pvs.get(i).altAllele;

                    String finalRepeatingUnit = allele;
                    int finalRepeatUnitLength = allele.length();

                    for (int repeatUnitLength = allele.length() - 1; repeatUnitLength >= 1; repeatUnitLength--) {
                        if (allele.length() % repeatUnitLength == 0) {
                            String repeatingUnit = allele.substring(0, repeatUnitLength);

                            boolean repeatIsComplete = true;
                            for (int j = 0; j < allele.length(); j += repeatUnitLength) {
                                if (!allele.substring(j, j + repeatUnitLength).equals(repeatingUnit)) {
                                    repeatIsComplete = false;
                                    break;
                                }
                            }

                            if (repeatIsComplete) {
                                if (repeatUnitLength < finalRepeatUnitLength) {
                                    finalRepeatingUnit = repeatingUnit;
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

                    if (!finalRepeatingUnit.isEmpty() && matchesExistingStr && finalRepeatUnitLength >= 2) {
                        pvs.get(i).pt = pvs.get(i).pt == PutativeVariantType.DEL ? PutativeVariantType.STR_CON : PutativeVariantType.STR_EXP;
                    }
                }
            }

            for (int i = 0; i < pvs.size(); i++) {
                if (!pvs.get(i).isFilteredOut && pvs.get(i).pt == PutativeVariantType.INS) {
                    String allele = pvs.get(i).altAllele;

                    int prevBasesStart1 = pvs.get(i).refPos - allele.length();
                    int prevBasesEnd1 = pvs.get(i).refPos;

                    int nextBasesStart1 = pvs.get(i).refPos;
                    int nextBasesEnd1 = pvs.get(i).refPos + allele.length();

                    String prevBases = prevBasesStart1 >= 0 ? refSequence.substring(prevBasesStart1, prevBasesEnd1) : null;
                    String nextBases = nextBasesEnd1 < refSequence.length() ? refSequence.substring(nextBasesStart1, nextBasesEnd1) : null;

                    if (allele.length() >= 10 && (allele.equals(prevBases) || allele.equals(nextBases))) {
                        pvs.get(i).pt = PutativeVariantType.TD;
                    }
                }
            }

            for (int i = 0; i < pvs.size(); i++) {
                for (int j = 0; j < pvs.get(i).altAllele.length(); j++) {
                    int kopos = pvs.get(i).altPos + j;

                    if (kopos < kmerOrigin.length() && kmerOrigin.charAt(kopos) == '.') {
                        pvs.get(i).numContainedNovelKmers++;
                    }
                }

                for (int leftpos = pvs.get(i).altPos - 1; leftpos >= 0; leftpos--) {
                    if (leftpos >= 0 && leftpos < kmerOrigin.length() && kmerOrigin.charAt(leftpos) == '.') {
                        pvs.get(i).numFlankingNovelKmers++;
                    } else {
                        break;
                    }
                }

                for (int rightpos = pvs.get(i).altPos + pvs.get(i).altAllele.length(); rightpos < kmerOrigin.length(); rightpos++) {
                    if (rightpos >= 0 && rightpos < kmerOrigin.length() && kmerOrigin.charAt(rightpos) == '.') {
                        pvs.get(i).numFlankingNovelKmers++;
                    } else {
                        break;
                    }
                }

                pvs.get(i).psb = psb;
                pvs.get(i).rsb = rsb;
                pvs.get(i).asb = asb;
                pvs.get(i).csb = csb;
                pvs.get(i).ksb = ksb;
                pvs.get(i).sr  = sr;
            }

            for (int i = 0; i < pvs.size(); i++) {
                if (!pvs.get(i).isFilteredOut) {
                    log.debug("  variant: {}", pvs.get(i));

                    goodVariants.put(pvs.get(i).altPos, pvs.get(i));
                }
            }

            List<Integer> variantPositions = new ArrayList<Integer>(goodVariants.keySet());
            Collections.sort(variantPositions);

            for (int i = 0; i < variantPositions.size(); i++) {
                int i0 = i - 1;
                int i1 = i + 1;

                int pos0 = (i0 >= 0) ? variantPositions.get(i0) : -1;
                int posT = variantPositions.get(i);
                int pos1 = (i1 < variantPositions.size()) ? variantPositions.get(i1) : -1;

                if (goodVariants.containsKey(pos0)) {
                    goodVariants.get(posT).distanceFromVariantLeft = goodVariants.get(posT).altPos - goodVariants.get(pos0).altPos;
                }

                if (goodVariants.containsKey(pos1)) {
                    goodVariants.get(posT).distanceFromVariantRight = goodVariants.get(pos1).altPos - goodVariants.get(posT).altPos;
                }
            }

            log.debug("");
        }

        return goodVariants;
    }

    private Map<Integer, PutativeVariant> callDeNovoVariants(Entry e) {
        Map<Integer, PutativeVariant> dn = new TreeMap<Integer, PutativeVariant>();

        Map<Integer, PutativeVariant> v0 = getVariants(e, 0);
        Map<Integer, PutativeVariant> v1 = getVariants(e, 1);

        Set<Integer> indices = new TreeSet<Integer>();
        indices.addAll(v0.keySet());
        indices.addAll(v1.keySet());

        for (int index : indices) {
            if (v0.containsKey(index) && v1.containsKey(index)) {
                if (v0.get(index).pt == PutativeVariantType.SNP) {
                    if (v0.get(index).numContainedNovelKmers > 0 &&
                        v0.get(index).altAllele.equals(v1.get(index).altAllele) &&
                        v0.get(index).distanceFromVariantLeft > 50 &&
                        v0.get(index).distanceFromVariantRight > 50
                    ) {
                        dn.put(index, v0.get(index));
                    } else {
                        v0.get(index).pt = PutativeVariantType.UNCLASSIFIED;
                        dn.put(index, v0.get(index));
                    }
                } else {
                    dn.put(index, v0.get(index));
                }
            } else if (v0.containsKey(index) && !v1.containsKey(index) && (v0.get(index).numContainedNovelKmers > 0 || v0.get(index).numFlankingNovelKmers > 0)) {
                if (v0.get(index).pt == PutativeVariantType.SNP) {
                    v0.get(index).pt = PutativeVariantType.UNCLASSIFIED;
                    dn.put(index, v0.get(index));
                } else {
                    dn.put(index, v0.get(index));
                }
            } else if (!v0.containsKey(index) && v1.containsKey(index) && (v1.get(index).numContainedNovelKmers > 0 || v1.get(index).numFlankingNovelKmers > 0)) {
                if (v1.get(index).pt == PutativeVariantType.SNP) {
                    v1.get(index).pt = PutativeVariantType.UNCLASSIFIED;
                    dn.put(index, v1.get(index));
                } else {
                    dn.put(index, v1.get(index));
                }
            }
        }

        return dn;
    }

    private boolean fuzzyEquals(PutativeVariant v, Map<String, String> ke) {
        String vRefAllele = v.refAllele;
        String vAltAllele = v.altAllele;
        String kRefAllele = ke.get("refAllele");
        String kAltAllele = ke.get("altAllele");

        boolean found = false;

        if (ke.get("event").equals("SNP") || ke.get("event").equals("INV")) {
            found = vRefAllele.equals(kRefAllele) && vAltAllele.equals(kAltAllele);
        } else if (ke.get("event").equals("INS") || ke.get("event").equals("STR_EXP") || ke.get("event").equals("TD")) {
            found = vAltAllele.equals(kAltAllele) || vAltAllele.equals(kAltAllele.substring(1, kAltAllele.length())) || vAltAllele.equals(kAltAllele.substring(0, kAltAllele.length() - 1));
        } else if (ke.get("event").equals("DEL") || ke.get("event").equals("STR_CON")) {
            found = vRefAllele.equals(kRefAllele) || vRefAllele.equals(kRefAllele.substring(1, kRefAllele.length())) || vRefAllele.equals(kRefAllele.substring(0, kRefAllele.length() - 1));
        }

        return found;
    }

    @Override
    public void execute() {
        Map<String, Entry> entries = loadEntriesWithNovelKmers();
        Map<String, Set<Map<String, String>>> knownEvents = loadExistingEvents();

        Map<PutativeVariantType, Integer> counts = new TreeMap<PutativeVariantType, Integer>();
        Map<String, Set<String>> found = new TreeMap<String, Set<String>>();
        Map<String, Set<String>> total = new TreeMap<String, Set<String>>();

        Map<String, String> foundClasses = new TreeMap<String, String>();
        Map<String, String> realClasses = new TreeMap<String, String>();

        for (String contigName : knownEvents.keySet()) {
            for (Map<String, String> ke : knownEvents.get(contigName)) {
                String id = ke.get("id");
                String event = ke.get("event");

                if (!found.containsKey(event)) { found.put(event, new HashSet<String>()); }
                if (!total.containsKey(event)) { total.put(event, new HashSet<String>()); }

                total.get(event).add(id);

                realClasses.put(id, event);
            }
        }

        for (String contigName : entries.keySet()) {
            Entry e = entries.get(contigName);

            Map<Integer, PutativeVariant> vs = callDeNovoVariants(e);
            Set<Map<String, String>> kes = (knownEvents.containsKey(contigName)) ? knownEvents.get(contigName) : null;

            boolean sawEvent = false;
            for (PutativeVariant v : vs.values()) {
                if (knownEvents.containsKey(contigName)) {
                    for (Map<String, String> ke : kes) {

                        if (fuzzyEquals(v, ke)) {
                            String id = ke.get("id");
                            String event = ke.get("event");

                            found.get(event).add(id);

                            foundClasses.put(id, v.pt.toString());
                        }
                    }
                }

                if (!counts.containsKey(v.pt)) {
                    counts.put(v.pt, 0);
                } else {
                    counts.put(v.pt, counts.get(v.pt) + 1);
                }

                if (v.pt == PutativeVariantType.UNCLASSIFIED) {
                    sawEvent = true;
                }
            }

            if (SHOW && sawEvent) {
                PutativeVariant v = vs.values().iterator().next();

                log.info("{}", v.sr.getSAMString().trim());
                log.info("{}", v.psb.toString());
                log.info("{}", v.rsb.toString());
                log.info("{}", v.asb.toString());
                log.info("{}", v.csb.toString());
                log.info("{}", v.ksb.toString());

                for (PutativeVariant v1 : vs.values()) {
                    log.info("{}", v1);
                }

                log.info("");
            }
        }

        log.info(" counts {}", Joiner.on("; ").withKeyValueSeparator("=").join(counts));

        log.info("  found");
        for (String event : found.keySet()) {
            log.info("    {} {}/{}", event, found.get(event).size(), total.get(event).size());
        }

        DataTable confusion = new DataTable("confusion", "confusion matrix");
        confusion.addColumns("real", "found", "count");
        for (String id : realClasses.keySet()) {
            if (realClasses.containsKey(id) && foundClasses.containsKey(id)) {
                String realEvent = realClasses.get(id);
                String foundEvent = foundClasses.get(id);

                String pk = realEvent + "." + foundEvent;
                confusion.set(pk, "real", realEvent);
                confusion.set(pk, "found", foundEvent);
                confusion.increment(pk, "count");
            }
        }

        log.info("  confusion\n{}", confusion);
    }
}
