package uk.ac.ox.well.indiana.commands.simulate;

import com.google.common.base.Joiner;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileReader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.util.StringUtil;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.io.table.TableReader;
import uk.ac.ox.well.indiana.utils.sequence.SequenceUtils;

import java.io.File;
import java.util.*;

public class CallDeNovoVariants2 extends Module {
    @Argument(fullName="metrics", shortName="m", doc="Annotated contigs")
    public File METRICS;

    @Argument(fullName="knownTable", shortName="k", doc="Table of known events", required=false)
    public File KNOWN_TABLE;

    @Argument(fullName="ref0", shortName="r0", doc="Reference for parent 0")
    public IndexedFastaSequenceFile REF0;

    @Argument(fullName="ref1", shortName="r1", doc="Reference for parent 1")
    public IndexedFastaSequenceFile REF1;

    @Argument(fullName="bam0", shortName="b0", doc="BAM for parent 0")
    public SAMFileReader BAM0;

    @Argument(fullName="bam1", shortName="b1", doc="BAM for parent 1")
    public SAMFileReader BAM1;

    private Map<String, Set<Map<String, String>>> loadExistingEvents() {
        Map<String, Set<Map<String, String>>> knownVariants = new HashMap<String, Set<Map<String, String>>>();

        if (KNOWN_TABLE != null) {
            TableReader tr = new TableReader(KNOWN_TABLE);

            for (Map<String, String> te : tr) {
                if (te.get("variantFound").equals("TRUE")) {
                    String contig1 = te.get("contig1");

                    if (!knownVariants.containsKey(contig1)) {
                        knownVariants.put(contig1, new HashSet<Map<String, String>>());
                    }

                    knownVariants.get(contig1).add(te);
                }
            }
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

    private void align(ReferenceSequence qseq, String kmerOrigin, SAMRecord alignment, IndexedFastaSequenceFile ref) {
        StringBuilder target = getTarget(alignment, ref);
        StringBuilder query = getQuery(qseq, alignment);
        StringBuilder ko = getKmerOrigin(qseq, kmerOrigin, alignment);

        StringBuilder matches = new StringBuilder(StringUtil.repeatCharNTimes(' ', target.length()));
        StringBuilder tpos = new StringBuilder(StringUtil.repeatCharNTimes(' ', target.length() > query.length() ? target.length() : query.length()));
        StringBuilder qpos = new StringBuilder(StringUtil.repeatCharNTimes(' ', target.length() > query.length() ? target.length() : query.length()));

        for (int i = 0; i < target.length(); i += 10) {
            String posString = String.valueOf(i);
            tpos.replace(i, i + posString.length(), posString);
            qpos.replace(i, i + posString.length(), posString);
        }

        int pos = 0;
        for (CigarElement ce : getCigar(alignment)) {
            if (ce.getOperator().equals(CigarOperator.M)) {
                for (int i = 0; i < ce.getLength(); i++) {
                    if (target.charAt(i + pos) == query.charAt(i + pos)) {
                        matches.setCharAt(i + pos, '|');
                    }
                }

                pos += ce.getLength();
            } else if (ce.getOperator().equals(CigarOperator.I)) {
                target.insert(pos, StringUtil.repeatCharNTimes('-', ce.getLength()));
                tpos.insert(pos, StringUtil.repeatCharNTimes(' ', ce.getLength()));
                matches.insert(pos, StringUtil.repeatCharNTimes(' ', ce.getLength()));
                pos += ce.getLength();
            } else if (ce.getOperator().equals(CigarOperator.D)) {
                query.insert(pos, StringUtil.repeatCharNTimes('-', ce.getLength()));
                qpos.insert(pos, StringUtil.repeatCharNTimes(' ', ce.getLength()));
                ko.insert(pos, StringUtil.repeatCharNTimes(' ', ce.getLength()));
                pos += ce.getLength();
            }
        }

        log.info("{} {}:{}-{} {} {}", qseq.getName(), alignment.getReferenceName(), alignment.getAlignmentStart(), alignment.getAlignmentEnd(), alignment.getReadNegativeStrandFlag() ? "-" : "+", getCigarString(alignment));
        log.info("{}", tpos);
        log.info("{}", target);
        log.info("{}", matches);
        log.info("{}", query);
        log.info("{}", ko);
        log.info("{}", qpos);
    }

    private Set<VariantContext> call(ReferenceSequence qseq, String kmerOrigin, SAMRecord alignment, IndexedFastaSequenceFile ref) {
        StringBuilder target = getTarget(alignment, ref);
        StringBuilder query = getQuery(qseq, alignment);
        StringBuilder ko = getKmerOrigin(qseq, kmerOrigin, alignment);

        Set<VariantContext> vcs = new HashSet<VariantContext>();

        int pos = 0, tpos = 0, qpos = 0;
        for (CigarElement ce : getCigar(alignment)) {
            if (ce.getOperator().equals(CigarOperator.M)) {
                for (int i = 0; i < ce.getLength(); i++) {
                    if (query.charAt(pos + i) != target.charAt(pos + i)) {
                        String targetAllele = String.valueOf(target.charAt(pos + i));
                        String queryAllele = String.valueOf(query.charAt(pos + i));

                        VariantContext vc = (new VariantContextBuilder())
                                .chr(qseq.getName())
                                .start(qpos + i)
                                .stop(qpos + i)
                                .alleles(targetAllele, queryAllele)
                                .make();

                        vcs.add(vc);
                    }
                }

                pos += ce.getLength();
                tpos += ce.getLength();
                qpos += ce.getLength();
            } else if (ce.getOperator().equals(CigarOperator.I)) {
                String targetAllele = target.substring(pos - 1, pos);
                String queryAllele = query.substring(pos - 1, pos + ce.getLength());

                VariantContext vc = (new VariantContextBuilder())
                        .chr(qseq.getName())
                        .start(qpos)
                        .stop(qpos)
                        .alleles(targetAllele, queryAllele)
                        .make();

                vcs.add(vc);

                target.insert(pos, StringUtil.repeatCharNTimes('-', ce.getLength()));
                pos += ce.getLength();
                qpos += ce.getLength();
            } else if (ce.getOperator().equals(CigarOperator.D)) {
                String targetAllele = target.substring(pos - 1, pos + ce.getLength());
                String queryAllele = query.substring(pos - 1, pos);

                VariantContext vc = (new VariantContextBuilder())
                        .chr(qseq.getName())
                        .start(qpos)
                        .stop(qpos + ce.getLength())
                        .alleles(targetAllele, queryAllele)
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

    private Set<VariantContext> callInversions(StringBuilder target, StringBuilder query, int kmerSize, Set<VariantContext> vcs) {
        Map<Integer, VariantContext> vcMap = new TreeMap<Integer, VariantContext>();
        for (VariantContext vc : vcs) {
            vcMap.put(vc.getStart(), vc);
        }

        for (int i = 0; i <= target.length() - kmerSize; i++) {
            String tFw = target.substring(i, i + kmerSize);
            String tRc = SequenceUtils.reverseComplement(tFw);

            for (int j = query.length() - kmerSize; j > i; j--) {
                String kFw = query.substring(j, j + kmerSize);

                if (tRc.equals(kFw)) {
                    String targetCandidateFw = target.substring(i, j + kmerSize);
                    String targetCandidateRc = SequenceUtils.reverseComplement(target.substring(i, j + kmerSize));
                    String queryCandidate = query.substring(i, j + kmerSize);

                    Set<Integer> variantPositions = new TreeSet<Integer>();
                    for (int q = i; q < j + kmerSize; q++) {
                        if (vcMap.containsKey(q)) {
                            variantPositions.add(q);
                        }
                    }

                    int sRc = SequenceUtils.numSegregatingSites(targetCandidateRc, queryCandidate);
                    int sFw = SequenceUtils.numSegregatingSites(targetCandidateFw, queryCandidate);
                    if (sRc == 0 && sFw > 0 && variantPositions.size() > 0) {
                        for (int variantPosition : variantPositions) {
                            vcMap.remove(variantPosition);
                        }

                        String chr = vcs.iterator().next().getChr();

                        VariantContext vc = (new VariantContextBuilder())
                                .chr(chr)
                                .start(i)
                                .stop(j + kmerSize - 1)
                                .alleles(targetCandidateFw, queryCandidate)
                                .make();

                        vcMap.put(i, vc);
                    }
                }
            }
        }

        return new HashSet<VariantContext>(vcMap.values());
    }

    private Set<VariantContext> refine(ReferenceSequence qseq, String kmerOrigin, SAMRecord alignment, IndexedFastaSequenceFile ref, Set<VariantContext> vcs) {
        StringBuilder target = getTarget(alignment, ref);
        StringBuilder query = getQuery(qseq, alignment);

        if (vcs.size() > 0) {
            vcs = callInversions(target, query, 10, vcs);
        }

        return vcs;
    }

    private Set<VariantContext> filter(ReferenceSequence qseq, String kmerOrigin, SAMRecord alignment, IndexedFastaSequenceFile ref, Set<VariantContext> vcs) {
        StringBuilder target = getTarget(alignment, ref);
        StringBuilder query = getQuery(qseq, alignment);

        Map<Integer, VariantContext> vcMap = new TreeMap<Integer, VariantContext>();
        Set<Integer> variantPositions = new TreeSet<Integer>();
        for (VariantContext vc : vcs) {
            int pos = vc.getStart();
            vcMap.put(pos, vc);
        }

        int pos = 0;
        int window = 10;
        boolean foundVariants;
        do {
            foundVariants = false;
            int shift = 0;
            for (int i = pos; i < pos + window; i++) {
                if (vcMap.containsKey(i)) {
                    foundVariants = true;

                    shift += vcMap.get(i).getEnd() - vcMap.get(i).getStart();

                    variantPositions.add(i);
                }
            }

            pos += window + shift;
        } while (foundVariants && pos + window < target.length());

        pos = target.length() - 1;
        do {
            foundVariants = false;
            int shift = 0;
            for (int i = pos; i >= pos - window; i--) {
                if (vcMap.containsKey(i)) {
                    foundVariants = true;

                    shift += vcMap.get(i).getEnd() - vcMap.get(i).getStart();

                    variantPositions.add(i);
                }
            }

            pos -= window - shift;
        } while (foundVariants && pos - window >= 0);

        for (int variantPosition : variantPositions) {
            vcMap.remove(variantPosition);
        }

        return new HashSet<VariantContext>(vcMap.values());
    }

    @Override
    public void execute() {
        log.info("Loading known events...");
        Map<String, Set<Map<String, String>>> knownEvents = loadExistingEvents();
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
        int multipleAlignments0 = 0;
        int multipleAlignments1 = 0;
        TableReader tr = new TableReader(METRICS);
        for (Map<String, String> te : tr) {
            String contigName = te.get("contigName");
            String seq = te.get("seq");
            String kmerOrigin = te.get("kmerOrigin");
            ReferenceSequence qseq = new ReferenceSequence(contigName, 0, seq.getBytes());

            Set<VariantContext> vcs0 = new HashSet<VariantContext>();
            Set<VariantContext> vcs1 = new HashSet<VariantContext>();

            if (align0.containsKey(contigName) && align0.get(contigName).size() == 1) {
                SAMRecord alignment = align0.get(contigName).iterator().next();

                align(qseq, kmerOrigin, alignment, REF0);
                vcs0 = call(qseq, kmerOrigin, alignment, REF0);
                vcs0 = refine(qseq, kmerOrigin, alignment, REF0, vcs0);
                vcs0 = filter(qseq, kmerOrigin, alignment, REF0, vcs0);
                log.info("");
            }

            if (align1.containsKey(contigName) && align1.get(contigName).size() == 1) {
                SAMRecord alignment = align1.get(contigName).iterator().next();

                align(qseq, kmerOrigin, alignment, REF1);
                vcs1 = call(qseq, kmerOrigin, alignment, REF1);
                vcs1 = refine(qseq, kmerOrigin, alignment, REF1, vcs1);
                vcs1 = filter(qseq, kmerOrigin, alignment, REF1, vcs1);
                log.info("");
            }

            log.info("known variants: {}", knownEvents.containsKey(contigName) ? knownEvents.get(contigName).size() : 0);
            if (knownEvents.containsKey(contigName)) {
                for (Map<String, String> ke : knownEvents.get(contigName)) {
                    log.info("  {}", Joiner.on(" ").withKeyValueSeparator("=").join(ke));
                }
            }
            log.info("called variants: {} {}", vcs0.size(), vcs1.size());

            log.info("-----");
        }

        log.info("  {}/{} contigs with multiple alignments", multipleAlignments0, tr.size());
        log.info("  {}/{} contigs with multiple alignments", multipleAlignments1, tr.size());
    }
}
