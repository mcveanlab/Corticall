package uk.ac.ox.well.indiana.commands.simulate;

import com.google.common.base.Joiner;
import htsjdk.samtools.*;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.StringUtil;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import htsjdk.variant.variantcontext.writer.Options;
import htsjdk.variant.variantcontext.writer.VariantContextWriter;
import htsjdk.variant.variantcontext.writer.VariantContextWriterBuilder;
import htsjdk.variant.vcf.*;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.containers.DataTable;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexKmer;
import uk.ac.ox.well.indiana.utils.io.table.TableReader;
import uk.ac.ox.well.indiana.utils.io.utils.LineReader;
import uk.ac.ox.well.indiana.utils.sequence.SequenceUtils;

import java.io.File;
import java.util.*;

public class CallDeNovoVariants extends Module {
    @Argument(fullName="metrics", shortName="m", doc="Annotated contigs")
    public File METRICS;

    @Argument(fullName="knowns", shortName="k", doc="Table of known events", required=false)
    public VCFFileReader KNOWNS;

    @Argument(fullName="ref0", shortName="r0", doc="Reference for parent 0")
    public IndexedFastaSequenceFile REF0;

    @Argument(fullName="ref1", shortName="r1", doc="Reference for parent 1")
    public IndexedFastaSequenceFile REF1;

    @Argument(fullName="bam0", shortName="b0", doc="BAM for parent 0")
    public SAMFileReader BAM0;

    @Argument(fullName="bam1", shortName="b1", doc="BAM for parent 1")
    public SAMFileReader BAM1;

    @Argument(fullName="proximityThreshold", shortName="t", doc="Proximity threshold variant filter")
    public Integer PROXIMITY_THRESHOLD = 47;

    @Argument(fullName="novelKmerVariantMap", shortName="n", doc="Novel kmer variant map")
    public File NOVEL_KMER_VARIANT_MAP;

    @Output
    public File out;

    private class VariantInfo {
        private String variantId;
        private String vclass;
        private String vchr;
        private int vstart;
        private int vstop;

        @Override
        public boolean equals(Object o) {
            if (this == o) return true;
            if (o == null || getClass() != o.getClass()) return false;

            VariantInfo that = (VariantInfo) o;

            if (vstart != that.vstart) return false;
            if (vstop != that.vstop) return false;
            if (variantId != null ? !variantId.equals(that.variantId) : that.variantId != null) return false;
            if (vchr != null ? !vchr.equals(that.vchr) : that.vchr != null) return false;
            if (vclass != null ? !vclass.equals(that.vclass) : that.vclass != null) return false;

            return true;
        }

        @Override
        public int hashCode() {
            int result = variantId != null ? variantId.hashCode() : 0;
            result = 31 * result + (vclass != null ? vclass.hashCode() : 0);
            result = 31 * result + (vchr != null ? vchr.hashCode() : 0);
            result = 31 * result + vstart;
            result = 31 * result + vstop;
            return result;
        }
    }

    private Set<String> weird = new HashSet<String>();

    private Map<CortexKmer, Set<VariantInfo>> loadNovelKmerVariantMap() {
        TableReader tr = new TableReader(NOVEL_KMER_VARIANT_MAP);

        Map<CortexKmer, Set<VariantInfo>> vis = new HashMap<CortexKmer, Set<VariantInfo>>();
        for (Map<String, String> te : tr) {
            CortexKmer ck = new CortexKmer(te.get("kmer"));

            VariantInfo vi = new VariantInfo();
            vi.variantId = te.get("variantId");
            vi.vclass = te.get("vclass");
            vi.vchr = te.get("vchr");
            vi.vstart = Integer.valueOf(te.get("vstart"));
            vi.vstop = Integer.valueOf(te.get("vstop"));

            //vis.put(ck, vi);

            if (!vis.containsKey(ck)) {
                vis.put(ck, new HashSet<VariantInfo>());
            }
            vis.get(ck).add(vi);
        }

        return vis;
    }

    private Map<String, Set<VariantContext>> loadExistingEvents() {
        Map<String, Set<VariantContext>> knownVariants = new HashMap<String, Set<VariantContext>>();

        for (VariantContext vc : KNOWNS) {
            if (!knownVariants.containsKey(vc.getChr())) {
                knownVariants.put(vc.getChr(), new HashSet<VariantContext>());
            }

            knownVariants.get(vc.getChr()).add(vc);
        }

        return knownVariants;
    }

    private Map<String, Set<SAMRecord>> loadBamRecords(SAMFileReader bam) {
        Map<String, Set<SAMRecord>> records = new HashMap<String, Set<SAMRecord>>();

        for (SAMRecord sr : bam) {
            String contigName = sr.getReadName();

            if (!records.containsKey(contigName)) {
                records.put(contigName, new HashSet<SAMRecord>());
            }

            records.get(contigName).add(sr);
        }

        return records;
    }

    private List<CigarElement> getCigar(SAMRecord alignment) {
        List<CigarElement> ces = new ArrayList<CigarElement>();

        if (alignment.getReadNegativeStrandFlag()) {
            for (int i = alignment.getCigar().numCigarElements() - 1; i >= 0; i--) {
                ces.add(alignment.getCigar().getCigarElement(i));
            }
        } else {
            ces = alignment.getCigar().getCigarElements();
        }

        return ces;
    }

    private StringBuilder getQuery(ReferenceSequence qseq, SAMRecord alignment) {
        StringBuilder query = new StringBuilder(new String(qseq.getBases()));

        List<CigarElement> cigar = getCigar(alignment);

        CigarElement ceLast = cigar.get(cigar.size() - 1);
        if (ceLast.getOperator().equals(CigarOperator.H) || ceLast.getOperator().equals(CigarOperator.S)) {
            query.delete(query.length() - ceLast.getLength(), query.length());
        }

        CigarElement ceFirst = cigar.get(0);
        if (ceFirst.getOperator().equals(CigarOperator.H) || ceFirst.getOperator().equals(CigarOperator.S)) {
            query.delete(0, ceFirst.getLength());
        }

        return query;
    }

    private StringBuilder getQueryPositionTrack(ReferenceSequence qseq, SAMRecord alignment) {
        StringBuilder qpos = new StringBuilder(StringUtil.repeatCharNTimes(' ', qseq.length()));

        for (int i = 0; i < qseq.length(); i += 20) {
            String posString = String.valueOf(i);
            qpos.replace(i, i + posString.length(), posString);
        }

        List<CigarElement> cigar = getCigar(alignment);

        CigarElement ceLast = cigar.get(cigar.size() - 1);
        if (ceLast.getOperator().equals(CigarOperator.H) || ceLast.getOperator().equals(CigarOperator.S)) {
            qpos.delete(qpos.length() - ceLast.getLength(), qpos.length());
        }

        CigarElement ceFirst = cigar.get(0);
        if (ceFirst.getOperator().equals(CigarOperator.H) || ceFirst.getOperator().equals(CigarOperator.S)) {
            qpos.delete(0, ceFirst.getLength());
        }

        return qpos;
    }

    private StringBuilder getKmerOrigin(ReferenceSequence query, String kmerOrigin, SAMRecord alignment) {
        StringBuilder ko = new StringBuilder(kmerOrigin);
        ko.append(StringUtil.repeatCharNTimes(ko.charAt(ko.length() - 1), query.length() - ko.length()));

        List<CigarElement> cigar = getCigar(alignment);

        CigarElement ceLast = cigar.get(cigar.size() - 1);
        if (ceLast.getOperator().equals(CigarOperator.H) || ceLast.getOperator().equals(CigarOperator.S)) {
            ko.delete(ko.length() - ceLast.getLength(), ko.length());
        }

        CigarElement ceFirst = cigar.get(0);
        if (ceFirst.getOperator().equals(CigarOperator.H) || ceFirst.getOperator().equals(CigarOperator.S)) {
            ko.delete(0, ceFirst.getLength());
        }

        return ko;
    }

    private StringBuilder getTarget(SAMRecord alignment, IndexedFastaSequenceFile ref) {
        ReferenceSequence tseq = ref.getSubsequenceAt(alignment.getReferenceName(), alignment.getAlignmentStart(), alignment.getAlignmentEnd());
        StringBuilder target = new StringBuilder(new String(tseq.getBases()));

        if (alignment.getReadNegativeStrandFlag()) {
            target = new StringBuilder(SequenceUtils.reverseComplement(target.toString()));
        }

        return target;
    }

    private String getCigarString(SAMRecord alignment) {
        StringBuilder cigarString = new StringBuilder();

        List<CigarElement> cigar = getCigar(alignment);
        for (CigarElement ce : cigar) {
            cigarString.append(String.format("%d%s", ce.getLength(), ce.getOperator()));
        }

        return cigarString.toString();
    }

    private int[] align(ReferenceSequence qseq, String kmerOrigin, SAMRecord alignment, IndexedFastaSequenceFile ref) {
        StringBuilder target = getTarget(alignment, ref);
        StringBuilder query = getQuery(qseq, alignment);
        StringBuilder ko = getKmerOrigin(qseq, kmerOrigin, alignment);

        StringBuilder matches = new StringBuilder(StringUtil.repeatCharNTimes(' ', target.length()));

        StringBuilder qpos = getQueryPositionTrack(qseq, alignment);
        StringBuilder tpos = new StringBuilder(StringUtil.repeatCharNTimes(' ', target.length() > query.length() ? target.length() : query.length()));
        //StringBuilder qpos = new StringBuilder(StringUtil.repeatCharNTimes(' ', target.length() > query.length() ? target.length() : query.length()));

        String qss = new String(qseq.getBases());
        String tss = target.toString();

        int[] posarray = new int[qseq.length()];
        for (int i = 0; i < posarray.length; i++) {
            posarray[i] = -1;
        }

        for (int i = 0; i < target.length(); i += 20) {
            //String posString = String.valueOf(alignment.getAlignmentStart() + i);
            String posString = String.valueOf(i);
            tpos.replace(i, i + posString.length(), posString);
            //qpos.replace(i, i + posString.length(), posString);
        }

        List<CigarElement> ces = getCigar(alignment);
        int pos = 0, qp = 0, tp = 0;
        if (ces.get(0).getOperator().equals(CigarOperator.S) || ces.get(0).getOperator().equals(CigarOperator.H)) {
            qp = ces.get(0).getLength();
        }

        for (CigarElement ce : getCigar(alignment)) {
            if (ce.getOperator().equals(CigarOperator.M)) {
                for (int i = 0; i < ce.getLength(); i++) {
                    posarray[qp + i] = tp + i;

                    if (target.charAt(i + pos) == query.charAt(i + pos)) {
                        matches.setCharAt(i + pos, '|');
                    }
                }

                pos += ce.getLength();
                tp += ce.getLength();
                qp += ce.getLength();
            } else if (ce.getOperator().equals(CigarOperator.I)) {
                for (int i = 0; i < ce.getLength(); i++) {
                    posarray[qp + i] = tp;
                }

                target.insert(pos, StringUtil.repeatCharNTimes('-', ce.getLength()));
                tpos.insert(pos, StringUtil.repeatCharNTimes(' ', ce.getLength()));
                matches.insert(pos, StringUtil.repeatCharNTimes(' ', ce.getLength()));
                pos += ce.getLength();
                qp += ce.getLength();
            } else if (ce.getOperator().equals(CigarOperator.D)) {
                for (int i = 0; i < ce.getLength(); i++) {
                    posarray[qp] = tp + i;
                }

                query.insert(pos, StringUtil.repeatCharNTimes('-', ce.getLength()));
                qpos.insert(pos, StringUtil.repeatCharNTimes(' ', ce.getLength()));
                ko.insert(pos, StringUtil.repeatCharNTimes(' ', ce.getLength()));
                pos += ce.getLength();
                tp += ce.getLength();
            }
        }

        log.debug("{} {}:{}-{} {} {}", qseq.getName(), alignment.getReferenceName(), alignment.getAlignmentStart(), alignment.getAlignmentEnd(), alignment.getReadNegativeStrandFlag() ? "-" : "+", getCigarString(alignment));
        log.debug("{}", tpos);
        log.debug("{}", target);
        log.debug("{}", matches);
        log.debug("{}", query);
        log.debug("{}", ko);
        log.debug("{}", qpos);

        StringBuilder sa = new StringBuilder();
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < qss.length(); i++) {
            if (posarray[i] != -1) {
                sa.append(tss.charAt(posarray[i]));
                sb.append(qss.charAt(i));
            }
        }

        log.debug("{}", sa);
        log.debug("{}", sb);

        return posarray;
    }

    private Set<VariantContext> call(ReferenceSequence qseq, String kmerOrigin, SAMRecord alignment, IndexedFastaSequenceFile ref, int refIndex) {
        StringBuilder target = getTarget(alignment, ref);
        StringBuilder query = getQuery(qseq, alignment);
        StringBuilder ko = getKmerOrigin(qseq, kmerOrigin, alignment);

        Set<VariantContext> vcs = new HashSet<VariantContext>();

        List<CigarElement> ces = getCigar(alignment);
        int pos = 0, tpos = 0, qpos = 0;
        if (ces.get(0).getOperator().equals(CigarOperator.S) || ces.get(0).getOperator().equals(CigarOperator.H)) {
            qpos = ces.get(0).getLength();
        }

        for (CigarElement ce : ces) {
            if (ce.getOperator().equals(CigarOperator.M)) {
                for (int i = 0; i < ce.getLength(); i++) {
                    if (query.charAt(pos + i) != target.charAt(pos + i)) {
                        String targetAllele = String.valueOf(target.charAt(pos + i));
                        String queryAllele = String.valueOf(query.charAt(pos + i));

                        int novelKmersContained = ko.charAt(pos + i) == '.' ? 1 : 0;
                        int novelKmersLeft = 0;
                        int novelKmersRight = 0;

                        for (int q = pos + i - 1; q >= 0; q--) {
                            if (ko.charAt(q) == '.') { novelKmersLeft++; }
                            else { break; }
                        }

                        for (int q = pos + i + 1; q < query.length(); q++) {
                            if (ko.charAt(q) == '.') { novelKmersRight++; }
                            else { break; }
                        }

                        VariantContext vc = (new VariantContextBuilder())
                                .chr(qseq.getName())
                                .start(qpos + i)
                                .stop(qpos + i)
                                .alleles(targetAllele, queryAllele)
                                .attribute("novelKmersContained", novelKmersContained)
                                .attribute("novelKmersLeft", novelKmersLeft)
                                .attribute("novelKmersRight", novelKmersRight)
                                .attribute("chr" + refIndex, alignment.getReferenceName())
                                .attribute("start" + refIndex, alignment.getAlignmentStart() + tpos + i)
                                .attribute("end" + refIndex, alignment.getAlignmentStart() + tpos + i)
                                .make();

                        vcs.add(vc);
                    }
                }

                pos += ce.getLength();
                tpos += ce.getLength();
                qpos += ce.getLength();
            } else if (ce.getOperator().equals(CigarOperator.I)) {
                String prevBase = "N";
                String targetAllele, queryAllele;

                if (pos > 0) {
                    targetAllele = target.substring(pos - 1, pos);
                    queryAllele = query.substring(pos - 1, pos + ce.getLength());
                } else {
                    targetAllele = prevBase;
                    queryAllele = prevBase + query.substring(pos, pos + ce.getLength());
                }

                int novelKmersContained = 0;
                int novelKmersLeft = 0;
                int novelKmersRight = 0;

                for (int q = pos - 1; q < pos + ce.getLength(); q++) {
                    if (q >= 0 && ko.charAt(q) == '.') { novelKmersContained++; }
                }

                for (int q = pos - 2; q >= 0; q--) {
                    if (q >= 0 && ko.charAt(q) == '.') { novelKmersLeft++; }
                    else { break; }
                }

                for (int q = pos + ce.getLength(); q < query.length(); q++) {
                    if (q >= 0 && ko.charAt(q) == '.') { novelKmersRight++; }
                    else { break; }
                }

                VariantContext vc = (new VariantContextBuilder())
                        .chr(qseq.getName())
                        .start(qpos - 1)
                        .stop(qpos - 1)
                        .alleles(targetAllele, queryAllele)
                        //.stop(qpos - 1 + ce.getLength())
                        //.alleles(queryAllele, targetAllele)
                        .attribute("novelKmersContained", novelKmersContained)
                        .attribute("novelKmersLeft", novelKmersLeft)
                        .attribute("novelKmersRight", novelKmersRight)
                        .attribute("chr" + refIndex, alignment.getReferenceName())
                        .attribute("start" + refIndex, alignment.getAlignmentStart() + tpos)
                        .attribute("end" + refIndex, alignment.getAlignmentStart() + tpos)
                        .make();

                vcs.add(vc);

                target.insert(pos, StringUtil.repeatCharNTimes('-', ce.getLength()));
                pos += ce.getLength();
                qpos += ce.getLength();
            } else if (ce.getOperator().equals(CigarOperator.D)) {
                String prevBase = "N";

                String targetAllele, queryAllele;
                if (pos > 0) {
                    targetAllele = target.substring(pos - 1, pos + ce.getLength());
                    queryAllele = query.substring(pos - 1, pos);
                } else {
                    targetAllele = prevBase + target.substring(pos, pos + ce.getLength());
                    queryAllele = prevBase;
                }

                int novelKmersContained = (pos > 0 && ko.charAt(pos - 1) == '.') ? 1 : 0;
                int novelKmersLeft = 0;
                int novelKmersRight = 0;

                for (int q = pos - 2; q >= 0; q--) {
                    if (q >= 0 && ko.charAt(q) == '.') { novelKmersLeft++; }
                    else { break; }
                }

                for (int q = pos + ce.getLength() + 1; q < query.length(); q++) {
                    if (q >= 0 && ko.charAt(q) == '.') { novelKmersRight++; }
                    else { break; }
                }

                VariantContext vc = (new VariantContextBuilder())
                        .chr(qseq.getName())
                        .start(qpos - 1)
                        .stop(qpos - 1 + ce.getLength())
                        .alleles(targetAllele, queryAllele)
                        //.stop(qpos - 1)
                        //.alleles(queryAllele, targetAllele)
                        .attribute("novelKmersContained", novelKmersContained)
                        .attribute("novelKmersLeft", novelKmersLeft)
                        .attribute("novelKmersRight", novelKmersRight)
                        .attribute("chr" + refIndex, alignment.getReferenceName())
                        .attribute("start" + refIndex, alignment.getAlignmentStart() + tpos)
                        .attribute("end" + refIndex, alignment.getAlignmentStart() + tpos)
                        .make();

                vcs.add(vc);

                query.insert(pos, StringUtil.repeatCharNTimes('-', ce.getLength()));
                ko.insert(pos, StringUtil.repeatCharNTimes(' ', ce.getLength()));
                pos += ce.getLength();
                tpos += ce.getLength();
            }
        }

        return vcs;
    }

    private Set<VariantContext> callInversions(ReferenceSequence qseq, String kmerOrigin, SAMRecord alignment, IndexedFastaSequenceFile ref, int kmerSize, Set<VariantContext> vcs, int refIndex) {
        Map<Integer, VariantContext> vcMap = new TreeMap<Integer, VariantContext>();
        for (VariantContext vc : vcs) {
            vcMap.put(vc.getStart(), vc);
        }

        StringBuilder target = getTarget(alignment, ref);
        StringBuilder query = getQuery(qseq, alignment);
        StringBuilder ko = getKmerOrigin(qseq, kmerOrigin, alignment);

        Map<String, Set<Integer>> qKmers = new HashMap<String, Set<Integer>>();

        String qStr = new String(qseq.getBases());
        for (int i = qStr.length() - kmerSize; i >= 0; i--) {
            String kmer = qStr.substring(i, i + kmerSize);

            if (!qKmers.containsKey(kmer)) {
                qKmers.put(kmer, new TreeSet<Integer>(new Comparator<Integer>() {
                    @Override
                    public int compare(Integer o1, Integer o2) {
                        return o1 < o2 ? 1 : -1;
                    }
                }));
            }
            qKmers.get(kmer).add(i + kmerSize - 1);
        }

        for (int tStart = 0; tStart <= target.length() - kmerSize; tStart++) {
            String tKmer = SequenceUtils.reverseComplement(target.substring(tStart, tStart + kmerSize));

            if (qKmers.containsKey(tKmer)) {
                for (int qEnd : qKmers.get(tKmer)) {
                    int tEnd = tStart + 1;
                    int qStart = qEnd - 1;

                    while (tEnd < target.length() && qStart >= 0 && target.charAt(tEnd) == SequenceUtils.complement(qStr.charAt(qStart))) {
                        tEnd++;
                        qStart--;
                    }

                    if (tStart - kmerSize >= 0 && tEnd + kmerSize < target.length() && qStart - kmerSize >= 0 && qEnd + kmerSize < qStr.length()) {
                        String tLeftFlank = target.substring(tStart - kmerSize, tStart);
                        String tRightFlank = target.substring(tEnd, tEnd + kmerSize);

                        String qLeftFlank = qStr.substring(qStart - kmerSize + 1, qStart + 1);
                        String qRightFlank = qStr.substring(qEnd + 1, qEnd + kmerSize + 1);

                        String targetCandidateFw = target.substring(tStart, tEnd);
                        String targetCandidateRc = SequenceUtils.reverseComplement(target.substring(tStart, tEnd));
                        String queryCandidate = qStr.substring(qStart + 1, qEnd + 1);

                        if (tLeftFlank.equals(qLeftFlank) && tRightFlank.equals(qRightFlank) && queryCandidate.length() >= 10 && targetCandidateRc.equals(queryCandidate) && !targetCandidateFw.equals(queryCandidate)) {
                            Set<Integer> variantPositions = new TreeSet<Integer>();
                            for (int q = qStart; q <= qEnd; q++) {
                                if (vcMap.containsKey(q)) {
                                    variantPositions.add(q);
                                }
                            }

                            if (variantPositions.size() > 0) {
                                for (int variantPosition : variantPositions) {
                                    vcMap.remove(variantPosition);
                                }

                                int novelKmersContained = 0, novelKmersLeft = 0, novelKmersRight = 0;

                                for (int q = qStart; q < qEnd; q++) {
                                    if (q >= 0 && q < kmerOrigin.length() && kmerOrigin.charAt(q) == '.') { novelKmersContained++; }
                                }

                                for (int q = qStart - 1; q >= 0; q--) {
                                    if (q >= 0 && q < kmerOrigin.length() && kmerOrigin.charAt(q) == '.') { novelKmersLeft++; }
                                    else { break; }
                                }

                                for (int q = qEnd; q < query.length(); q++) {
                                    if (q >= 0 && q < kmerOrigin.length() && kmerOrigin.charAt(q) == '.') { novelKmersRight++; }
                                    else { break; }
                                }

                                VariantContext vc = (new VariantContextBuilder())
                                        .chr(alignment.getReadName())
                                        .start(qStart)
                                        .stop(qEnd - 1)
                                        .alleles(targetCandidateFw, queryCandidate)
                                        .attribute("novelKmersContained", novelKmersContained)
                                        .attribute("novelKmersLeft", novelKmersLeft)
                                        .attribute("novelKmersRight", novelKmersRight)
                                        .attribute("chr" + refIndex, alignment.getReferenceName())
                                        .attribute("start" + refIndex, alignment.getAlignmentStart() + tStart)
                                        .attribute("end" + refIndex, alignment.getAlignmentStart() + tEnd)
                                        .make();

                                vcMap.put(qStart, vc);
                            }
                        }
                    }
                }
            }
        }

        return new HashSet<VariantContext>(vcMap.values());
    }

    private Set<VariantContext> refine(ReferenceSequence qseq, String kmerOrigin, SAMRecord alignment, IndexedFastaSequenceFile ref, Set<VariantContext> vcs, int refIndex) {
        if (vcs.size() > 0) {
            vcs = callInversions(qseq, kmerOrigin, alignment, ref, 10, vcs, refIndex);
        }

        return vcs;
    }

    private Set<VariantContext> filter(ReferenceSequence qseq, String kmerOrigin, SAMRecord alignment, IndexedFastaSequenceFile ref, Set<VariantContext> vcs, int window, int minCount) {
        StringBuilder target = getTarget(alignment, ref);
        StringBuilder query = getQuery(qseq, alignment);

        Map<Integer, VariantContext> vcMap = new TreeMap<Integer, VariantContext>();
        for (VariantContext vc : vcs) {
            int pos = vc.getStart();
            vcMap.put(pos, vc);
        }

        List<CigarElement> ces = getCigar(alignment);
        CigarElement ceFirst = ces.get(0);
        int offsetLeft = (ceFirst.getOperator().equals(CigarOperator.S) || ceFirst.getOperator().equals(CigarOperator.H)) ? ceFirst.getLength() : 0;

        CigarElement ceLast = ces.get(ces.size() - 1);
        int offsetRight = (ceLast.getOperator().equals(CigarOperator.S) || ceLast.getOperator().equals(CigarOperator.H)) ? ceLast.getLength() : 0;

        Set<Integer> variantPositions = new TreeSet<Integer>();

        int firstVariantPos = query.length();
        int lastVariantPos = -1;

        for (int pos : vcMap.keySet()) {
            if (pos < firstVariantPos) { firstVariantPos = pos; }
            if (pos > lastVariantPos) { lastVariantPos = pos; }
        }

        if (firstVariantPos <= window + offsetLeft) {
            int endOfLastVariant = -1;
            Set<Integer> varPosLeft = new TreeSet<Integer>();
            for (int pos = 0; pos < query.length(); pos++) {
                if (vcMap.containsKey(pos)) {
                    int startOfThisVariant = vcMap.get(pos).getStart();
                    int endOfThisVariant = vcMap.get(pos).getEnd();

                    if (endOfLastVariant == -1 || startOfThisVariant - endOfLastVariant <= window) {
                        endOfLastVariant = endOfThisVariant;

                        varPosLeft.add(pos);
                    } else {
                        break;
                    }
                }
            }

            if (varPosLeft.size() > minCount) {
                variantPositions.addAll(varPosLeft);
            }
        }

        if (lastVariantPos >= target.length() - window - offsetRight) {
            int startOfLastVariant = -1;
            Set<Integer> varPosRight = new TreeSet<Integer>();
            for (int pos = query.length() - 1; pos > query.length(); pos++) {
                if (vcMap.containsKey(pos)) {
                    int startOfThisVariant = vcMap.get(pos).getStart();
                    int endOfThisVariant = vcMap.get(pos).getEnd();

                    if (startOfLastVariant == -1 || startOfLastVariant - endOfThisVariant <= window) {
                        startOfLastVariant = startOfThisVariant;

                        varPosRight.add(pos);
                    } else {
                        break;
                    }
                }
            }

            if (varPosRight.size() > minCount) {
                variantPositions.addAll(varPosRight);
            }
        }

        for (int variantPosition : variantPositions) {
            VariantContext vc = new VariantContextBuilder(vcMap.get(variantPosition))
                    .filter("PARTIAL")
                    .make();

            vcMap.put(variantPosition, vc);
        }

        return new HashSet<VariantContext>(vcMap.values());
    }

    private boolean canContainNovelKmers(VariantContext vc) {
        return vc.isSNP() || vc.isMNP() || vc.isSimpleInsertion();
    }

    private Set<VariantContext> markDenovos(Set<VariantContext> vcs0, Set<VariantContext> vcs1) {
        Set<VariantContext> vcsFinal = new HashSet<VariantContext>();

        Map<Interval, Set<VariantContext>> vcIntervalMap = new HashMap<Interval, Set<VariantContext>>();

        for (VariantContext vc : vcs0) {
            Interval interval = new Interval(vc.getChr(), vc.getStart(), vc.getEnd());
            if (!vcIntervalMap.containsKey(interval)) {
                vcIntervalMap.put(interval, new HashSet<VariantContext>());
            }

            vcIntervalMap.get(interval).add(vc);
        }

        for (VariantContext vc : vcs1) {
            Interval interval = new Interval(vc.getChr(), vc.getStart(), vc.getEnd());
            if (!vcIntervalMap.containsKey(interval)) {
                vcIntervalMap.put(interval, new HashSet<VariantContext>());
            }

            vcIntervalMap.get(interval).add(vc);
        }

        for (Interval interval : vcIntervalMap.keySet()) {
            Set<VariantContext> vcs = vcIntervalMap.get(interval);

            Map<String, Object> attrs = new HashMap<String, Object>();
            for (VariantContext vca : vcs) {
                attrs.putAll(vca.getAttributes());
            }

            VariantContext vc = vcs.iterator().next();

            int novelKmersContained = vc.getAttributeAsInt("novelKmersContained", 0);
            int novelKmersLeft = vc.getAttributeAsInt("novelKmersLeft", 0);
            int novelKmersRight = vc.getAttributeAsInt("novelKmersRight", 0);

            boolean isDeNovo = novelKmersContained > 0 && (novelKmersLeft > 0 || novelKmersRight > 0);

            VariantContext vcFinal = (new VariantContextBuilder(vc))
                    .attributes(attrs)
                    .attribute("DENOVO", isDeNovo)
                    .make();

            vcsFinal.add(vcFinal);
        }

        return vcsFinal;
    }

    private Set<VariantContext> annotateVariants(Set<VariantContext> vcsFinal) {
        Set<VariantContext> vcsFinalSorted = new TreeSet<VariantContext>(new Comparator<VariantContext>() {
            @Override
            public int compare(VariantContext o1, VariantContext o2) {
                return o1.getStart() < o2.getStart() ? -1 : 1;
            }
        });
        vcsFinalSorted.addAll(vcsFinal);

        List<VariantContext> vcsOrdered = new ArrayList<VariantContext>(vcsFinalSorted);

        Set<VariantContext> vcsFinalAnnotated = new HashSet<VariantContext>();
        for (int i = 0; i < vcsOrdered.size(); i++) {
            VariantContext vcLeft = i > 0 ? vcsOrdered.get(i-1) : null;
            VariantContext vcCenter = vcsOrdered.get(i);
            VariantContext vcRight = i < vcsOrdered.size() - 1 ? vcsOrdered.get(i+1) : null;

            int distanceToLeftVariant = Integer.MAX_VALUE;
            int distanceToRightVariant = Integer.MAX_VALUE;

            if (vcLeft != null) {
                distanceToLeftVariant = vcCenter.getStart() - vcLeft.getEnd();
            }

            if (vcRight != null) {
                distanceToRightVariant = vcRight.getStart() - vcCenter.getEnd();
            }

            VariantContext vcCenterMod = (new VariantContextBuilder(vcCenter))
                    .attribute("distanceToLeftVariant", distanceToLeftVariant)
                    .attribute("distanceToRightVariant", distanceToRightVariant)
                    .make();

            vcsFinalAnnotated.add(vcCenterMod);
        }

        return vcsFinalAnnotated;
    }

    private Set<VariantContext> filterVariants(Set<VariantContext> vcsFinal, int proximityThreshold) {
        Set<VariantContext> vcsFinalFiltered = new HashSet<VariantContext>();
        for (VariantContext vc : vcsFinal) {
            int distanceToLeftVariant = vc.getAttributeAsInt("distanceToLeftVariant", Integer.MAX_VALUE);
            int distanceToRightVariant = vc.getAttributeAsInt("distanceToRightVariant", Integer.MAX_VALUE);

            VariantContext vcFiltered = vc;
            if ((distanceToLeftVariant < proximityThreshold || distanceToRightVariant < proximityThreshold)) {
                vcFiltered = (new VariantContextBuilder(vc))
                        .filter("PROXIMITY")
                        .make();
            }

            vcsFinalFiltered.add(vcFiltered);
        }

        return vcsFinalFiltered;
    }

    private Set<VariantContext> classifyVariants(Set<VariantContext> vcsFinal, String seq) {
        Set<VariantContext> classifiedVCs = new HashSet<VariantContext>();

        for (VariantContext vc : vcsFinal) {
            if (vc.isIndel()) {
                String allele = vc.isSimpleInsertion() ? vc.getAlternateAllele(0).getBaseString() : vc.getReference().getBaseString();
                allele = allele.substring(1, allele.length());

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

                String prevBases = (vc.getStart() - finalRepeatUnitLength >= 0) ? seq.substring(vc.getStart() - finalRepeatUnitLength + 1, vc.getStart() + 1) : "";
                String nextBases = (vc.getEnd() + vc.getAlternateAllele(0).length() + finalRepeatUnitLength < seq.length()) ? seq.substring(vc.getStart() + allele.length(), vc.getStart() + allele.length() + finalRepeatUnitLength) : "";

                //int numRepeats = (finalRepeatUnitLength > 0) ? allele.length() / finalRepeatUnitLength : 0;

                String event = "unknown";
                if (finalRepeatingUnit.equals(prevBases) || finalRepeatingUnit.equals(nextBases)) {
                    if (vc.isSimpleInsertion()) {
                        if (allele.length() >= 10) {
                            event = "TD";
                        } else if (allele.length() >= 2 && allele.length() <= 5) {
                            event = "STR_EXP";
                        } else {
                            event = "INS";
                        }
                    } else if (allele.length() >= 2 && allele.length() <= 5) {
                        event = "STR_CON";
                    } else {
                        event = "DEL";
                    }

                    VariantContext classifiedVC = (new VariantContextBuilder(vc))
                            .attribute("event", event)
                            .attribute("repUnit", finalRepeatingUnit)
                            .make();

                    classifiedVCs.add(classifiedVC);
                } else {
                    if (vc.isSimpleInsertion()) {
                        event = "INS";
                    } else {
                        event = "DEL";
                    }

                    VariantContext classifiedVC = (new VariantContextBuilder(vc))
                            .attribute("event", event)
                            .make();

                    classifiedVCs.add(classifiedVC);
                }
            } else if (vc.isSNP()) {
                VariantContext classifiedVC = (new VariantContextBuilder(vc))
                        .attribute("event", "SNP")
                        .make();

                classifiedVCs.add(classifiedVC);
            } else if (vc.isMNP()) {
                VariantContext classifiedVC = (new VariantContextBuilder(vc))
                        .attribute("event", "INV")
                        .make();

                classifiedVCs.add(classifiedVC);
            }
        }

        return classifiedVCs;
    }

    private Set<VariantContext> verifyVariants(Set<VariantContext> vcsFinal, String seq, String contigName, SAMRecord a0, SAMRecord a1, int[] p0, int[] p1) {
        Set<VariantContext> classifiedVCs = new HashSet<VariantContext>();

        String t0 = a0 == null ? null : getTarget(a0, REF0).toString();
        String t1 = a1 == null ? null : getTarget(a1, REF1).toString();

        for (VariantContext vc : vcsFinal) {
            if (t0 != null && t1 != null && vc.isSNP() && vc.getAttributeAsBoolean("DENOVO", false) && ((vc.hasAttribute("chr0") && !vc.hasAttribute("chr1")) || (!vc.hasAttribute("chr0") && vc.hasAttribute("chr1")))) {
                int start = vc.getStart();

                int p0start = p0[start];
                String coAllele0 = "";
                if (p0start >= 0) {
                    //coAllele0 = String.valueOf(t0.charAt(p0start));
                    coAllele0 = t0.substring(p0start, p0start + vc.getReference().length());
                }

                int p1start = p1[start];
                String coAllele1 = "";
                if (p1start >= 0) {
                    //coAllele1 = String.valueOf(t1.charAt(p1start));
                    coAllele1 = t1.substring(p1start, p1start + vc.getReference().length());
                }

//              log.info("{} {} {} {} {}", start, vc.getReference(), vc.getAlternateAllele(0), coAllele0, coAllele1);
//              log.info("{}", coAllele0.equals(coAllele1));
//              log.info("--");

                if (!coAllele0.equals(coAllele1)) {
                    vc = (new VariantContextBuilder(vc))
                            .filter("INCONSISTENT")
                            .make();

                    classifiedVCs.add(vc);
                }
            } else {
                classifiedVCs.add(vc);
            }
        }

        return classifiedVCs;
    }

    private SAMSequenceDictionary buildSAMSequenceDictionary() {
        TableReader tr = new TableReader(METRICS);

        List<SAMSequenceRecord> ssrs = new ArrayList<SAMSequenceRecord>();
        for (Map<String, String> te : tr) {
            String contigName = te.get("contigName");
            String contig = te.get("seq");
            int contigLength = contig.length();

            SAMSequenceRecord ssr = new SAMSequenceRecord(contigName, contigLength);
            ssrs.add(ssr);

        }

        return new SAMSequenceDictionary(ssrs);
    }

    private Set<VariantInfo> getRelevantVariants(Map<CortexKmer, Set<VariantInfo>> vis, ReferenceSequence qseq, String kmerOrigin) {
        String seq = new String(qseq.getBases());

        int kmerSize = seq.length() - kmerOrigin.length() + 1;

        Set<VariantInfo> relevantVariants = new HashSet<VariantInfo>();
        for (int i = 0; i < kmerOrigin.length(); i++) {
            if (kmerOrigin.charAt(i) == '.') {
                CortexKmer ck = new CortexKmer(seq.substring(i, i + kmerSize));

                if (vis.containsKey(ck)) {
                    relevantVariants.addAll(vis.get(ck));
                }
            }
        }

        return relevantVariants;
    }

    @Override
    public void execute() {
        log.info("Loading known events...");
        Map<String, Set<VariantContext>> knownEvents = loadExistingEvents();
        int numEvents = 0;
        for (String contigName : knownEvents.keySet()) {
            numEvents += knownEvents.get(contigName).size();
        }
        log.info("  {} contigs with {} known events", knownEvents.size(), numEvents);

        log.info("Loading alignments...");
        Map<String, Set<SAMRecord>> align0 = loadBamRecords(BAM0);
        Map<String, Set<SAMRecord>> align1 = loadBamRecords(BAM1);
        log.info("  ref0: {}", align0.size());
        log.info("  ref1: {}", align1.size());

        log.info("Calling variants...");
        Map<String, VariantContext> variants = new HashMap<String, VariantContext>();

        log.info("Loading debugging info...");
        Map<CortexKmer, Set<VariantInfo>> vis = loadNovelKmerVariantMap();
        LineReader lr = new LineReader(new File("contigs.weird.txt"));
        while (lr.hasNext()) {
            String l = lr.getNextRecord();

            weird.add(l);
        }

        DataTable dt = new DataTable("recall", "recall");
        dt.addColumns("event", "known", "called", "filtered");

        DataTable dt2 = new DataTable("recall", "recall");
        dt2.addColumns("event", "known", "called", "filtered");

        VariantContextWriterBuilder vcwb = new VariantContextWriterBuilder();
        vcwb.setOutputFile(out);
        vcwb.unsetOption(Options.INDEX_ON_THE_FLY);
        vcwb.setOption(Options.ALLOW_MISSING_FIELDS_IN_HEADER);
        VariantContextWriter vcw = vcwb.build();

        Set<VCFHeaderLine> headerLines = new HashSet<VCFHeaderLine>();
        VCFHeader header = new VCFHeader(headerLines);
        header.setSequenceDictionary(buildSAMSequenceDictionary());

        vcw.writeHeader(header);

        TableReader tr = new TableReader(METRICS);
        int numContigsSeen = 0, numContigsProcessed = 0;
        for (Map<String, String> te : tr) {
            if (numContigsSeen % (tr.size()/10) == 0) {
                log.info("  {}/{} ({})", numContigsSeen, tr.size(), numContigsProcessed);
            }
            numContigsSeen++;

            //if (numContigsProcessed > 300) { break; }

            String contigName = te.get("contigName");
            int numNovelKmers = Integer.valueOf(te.get("refNovel"));
            String seq = te.get("seq");
            String kmerOrigin = te.get("kmerOrigin");
            ReferenceSequence qseq = new ReferenceSequence(contigName, 0, seq.getBytes());

            int numKnownDeNovos = 0, numKnownCompleteDenovos = 0;
            if (knownEvents.containsKey(contigName)) {
                Map<Integer, VariantContext> sortedKnowns = new TreeMap<Integer, VariantContext>();
                for (VariantContext ke : knownEvents.get(contigName)) {
                    int variantStart = ke.getStart();
                    sortedKnowns.put(variantStart, ke);

                    if (ke.isNotFiltered() && !ke.getAttributeAsString("denovo", "unknown").equals("unknown")) {
                        numKnownDeNovos++;

                        if (ke.getAttributeAsBoolean("isComplete", false)) {
                            numKnownCompleteDenovos++;
                        }
                    }
                }
                log.debug("known denovo variants: {}, complete: {}", numKnownDeNovos, numKnownCompleteDenovos);
                for (int variantStart : sortedKnowns.keySet()) {
                    VariantContext ke = sortedKnowns.get(variantStart);

                    if (ke.isNotFiltered() && !ke.getAttributeAsString("denovo", "unknown").equals("unknown") && ke.getAttributeAsBoolean("isComplete", false)) {
                        log.debug("  {}", ke);
                    }

                    dt.set(ke.getAttributeAsString("event", "unknown"), "event", ke.getAttributeAsString("event", "unknown"));
                    dt.increment(ke.getAttributeAsString("event", "unknown"), "known");

                    dt2.set(ke.getAttributeAsString("event", "unknown"), "event", ke.getAttributeAsString("event", "unknown"));
                    dt2.increment(ke.getAttributeAsString("event", "unknown"), "known");
                }
            }

            if (numNovelKmers > 0) {
                numContigsProcessed++;

                Set<VariantContext> vcs0 = new HashSet<VariantContext>();
                Set<VariantContext> vcs1 = new HashSet<VariantContext>();

                SAMRecord a0 = null, a1 = null;
                int[] p0 = null, p1 = null;
                if (align0.containsKey(contigName) && align0.get(contigName).size() == 1) {
                    a0 = align0.get(contigName).iterator().next();

                    if (!a0.getReadUnmappedFlag() && a0.getReadLength() >= 47) {
                        p0 = align(qseq, kmerOrigin, a0, REF0);
                        vcs0 = call(qseq, kmerOrigin, a0, REF0, 0);
                        vcs0 = refine(qseq, kmerOrigin, a0, REF0, vcs0, 0);
                        vcs0 = filter(qseq, kmerOrigin, a0, REF0, vcs0, 10, 5);
                        log.debug("");
                    }

                }

                if (align1.containsKey(contigName) && align1.get(contigName).size() == 1) {
                    a1 = align1.get(contigName).iterator().next();

                    if (!a1.getReadUnmappedFlag() && a1.getReadLength() >= 47) {
                        p1 = align(qseq, kmerOrigin, a1, REF1);
                        vcs1 = call(qseq, kmerOrigin, a1, REF1, 1);
                        vcs1 = refine(qseq, kmerOrigin, a1, REF1, vcs1, 1);
                        vcs1 = filter(qseq, kmerOrigin, a1, REF1, vcs1, 10, 5);
                        log.debug("");
                    }
                }

                Set<VariantContext> vcsMarked = markDenovos(vcs0, vcs1);
                Set<VariantContext> vcsAnnotated = annotateVariants(vcsMarked);
                Set<VariantContext> vcsFiltered = filterVariants(vcsAnnotated, PROXIMITY_THRESHOLD);
                Set<VariantContext> vcsClassified = classifyVariants(vcsFiltered, seq);
                Set<VariantContext> vcsFinal = verifyVariants(vcsClassified, seq, contigName, a0, a1, p0, p1);

                int numDeNovos = 0;
                for (VariantContext vc : vcsFinal) {
                    if (vc.getAttributeAsBoolean("DENOVO", false) && vc.isNotFiltered()) {
                        numDeNovos++;
                    }
                }

                log.debug("called denovo variants: {} out of {} total", numDeNovos, vcsFinal.size());
                Set<VariantContext> vcsFinalSorted = new TreeSet<VariantContext>(new Comparator<VariantContext>() {
                    @Override
                    public int compare(VariantContext o1, VariantContext o2) {
                        return o1.getStart() < o2.getStart() ? -1 : 1;
                    }
                });
                vcsFinalSorted.addAll(vcsFinal);

                for (VariantContext vc : vcsFinalSorted) {
                    vcw.add(vc);

                    if (vc.getAttributeAsBoolean("DENOVO", false) && vc.isNotFiltered()) {
                        String event = vc.getAttributeAsString("event", "unknown");
                        dt.set(event, "event", event);
                        dt.increment(event, "called");

                        log.debug("  id=noname event={} type={} length={} variantFound=TRUE matchedPos=irrelevant variantStart={} variantEnd={} vc={}",
                                vc.getAttributeAsString("event", "unknown"),
                                vc.getType(),
                                vc.getAlternateAllele(0).length(),
                                vc.getStart(),
                                vc.getEnd(),
                                vc
                        );
                    }
                }

                log.debug("filtered out variants: {} out of {} total", vcsFinal.size() - numDeNovos, vcsFinal.size());
                for (VariantContext vc : vcsFinalSorted) {
                    if (!vc.getAttributeAsBoolean("DENOVO", false) || !vc.isNotFiltered()) {
                        String event = vc.getAttributeAsString("event", "unknown");
                        dt.set(event, "event", event);
                        dt.increment(event, "filtered");

                        log.debug("  id={} event={} type={} length={} variantFound=TRUE matchedPos=irrelevant variantStart={} variantEnd={} vc={}",
                                Joiner.on(",").join(vc.getFilters()),
                                vc.getAttributeAsString("event", "unknown"),
                                vc.getType(),
                                vc.getAlternateAllele(0).length(),
                                vc.getStart(),
                                vc.getEnd(),
                                vc
                        );
                    }
                }

                if (numDeNovos > 0) {
                    log.debug("  {} {} {}", contigName, numDeNovos, numKnownDeNovos);
                }

                for (VariantContext vc : vcsFinalSorted) {
                    String id = String.format("%s:%d-%d;%s:%d-%d",
                            vc.getAttributeAsString("chr0", "none"), vc.getAttributeAsInt("start0", 0), vc.getAttributeAsInt("end0", 0),
                            vc.getAttributeAsString("chr1", "none"), vc.getAttributeAsInt("start1", 0), vc.getAttributeAsInt("end1", 0)
                    );

                    variants.put(id, vc);
                }

                Set<VariantInfo> relevantVariants = getRelevantVariants(vis, qseq, kmerOrigin);
                log.debug("variant sources of novel kmers: {}", relevantVariants.size());
                for (VariantInfo vi : relevantVariants) {
                    for (Set<VariantContext> vcs : knownEvents.values()) {
                        for (VariantContext vc : vcs) {
                            if (vc.getAttributeAsString("id", "unknown").equals(vi.variantId)) {
                                log.debug("  {} {}", Joiner.on(",").join(vc.getFilters()), vc);
                            }
                        }
                    }
                }

                log.debug("-----");
                getRelevantVariants(vis, qseq, kmerOrigin);
            }
        }

        for (String id : variants.keySet()) {
            VariantContext vc = variants.get(id);

            String event = vc.getAttributeAsString("event", "unknown");

            dt2.set(event, "event", event);
            if (vc.getAttributeAsBoolean("DENOVO", false) && vc.isNotFiltered()) {
                dt2.increment(event, "called");
            } else {
                dt2.increment(event, "filtered");
            }
        }

        log.info("\n{}", dt);
        log.info("");
        log.info("\n{}", dt2);

        /*
        for (String id : variants.keySet()) {
            VariantContext vc = variants.get(id);

            vcw.add(vc);
        }
        */

        vcw.close();
    }
}
