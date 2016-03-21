package uk.ac.ox.well.indiana.commands.gg;

import com.google.common.base.Joiner;
import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.reference.FastaSequenceFile;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.alignment.kmer.KmerLookup;
import uk.ac.ox.well.indiana.utils.alignment.pairwise.ExternalAligner;
import uk.ac.ox.well.indiana.utils.alignment.pairwise.LastzAligner;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.containers.ContainerUtils;
import uk.ac.ox.well.indiana.utils.exceptions.IndianaException;
import uk.ac.ox.well.indiana.utils.sequence.SequenceUtils;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.*;

public class AlignVars extends Module {
    @Argument(fullName="vars", shortName="v", doc="Vars")
    public FastaSequenceFile VARS;

    @Argument(fullName="ref", shortName="r", doc="Reference")
    public File REF;

    @Argument(fullName="kmerSize", shortName="k", doc="Kmer size")
    public Integer KMER_SIZE = 47;

    @Override
    public void execute() {
        KmerLookup kl = new KmerLookup(REF, KMER_SIZE);

        IndexedFastaSequenceFile ifsf;
        try {
            ifsf = new IndexedFastaSequenceFile(REF);
        } catch (FileNotFoundException e) {
            throw new IndianaException("Could not open ref");
        }

        Map<String, Integer> lengths = new HashMap<String, Integer>();
        ReferenceSequence rs;
        while ((rs = ifsf.nextSequence()) != null) {
            String name = rs.getName().split("\\s+")[0];
            name = name.replaceAll("_v3", "")
                    .replaceAll("Pf3D7_", "")
                    .replaceAll("PfHB3_", "")
                    .replaceAll("PfDD2_", "")
                    .replaceAll("_T[TLR]", "");

            lengths.put(name, rs.length());
        }

        IntervalTreeMap<String> itm = new IntervalTreeMap<String>();

        itm.put(new Interval("01", 0, 92900), "5p");
        itm.put(new Interval("01", lengths.get("01") - 41000, lengths.get("01")), "3p");
        itm.put(new Interval("02", 0, 105800), "5p");
        itm.put(new Interval("02", lengths.get("02") - 63000, lengths.get("02")), "3p");
        itm.put(new Interval("03", 0, 70630), "5p");
        itm.put(new Interval("03", lengths.get("03") - 37000, lengths.get("03")), "3p");
        itm.put(new Interval("04", 0, 91420), "5p");
        itm.put(new Interval("04", lengths.get("04") - 65000, lengths.get("04")), "3p");
        itm.put(new Interval("05", 0, 37900), "5p");
        itm.put(new Interval("05", lengths.get("05") - 20000, lengths.get("05")), "3p");
        itm.put(new Interval("06", 0, 72350), "5p");
        itm.put(new Interval("06", lengths.get("06") - 72000, lengths.get("06")), "3p");
        itm.put(new Interval("07", 0, 77100), "5p");
        itm.put(new Interval("07", lengths.get("07") - 60000, lengths.get("07")), "3p");
        itm.put(new Interval("08", 0, 73560), "5p");
        itm.put(new Interval("08", lengths.get("08") - 55000, lengths.get("08")), "3p");
        itm.put(new Interval("09", 0, 79100), "5p");
        itm.put(new Interval("09", lengths.get("09") - 61000, lengths.get("09")), "3p");
        itm.put(new Interval("10", 0, 68970), "5p");
        itm.put(new Interval("10", lengths.get("10") - 43000, lengths.get("10")), "3p");
        itm.put(new Interval("11", 0, 110000), "5p");
        itm.put(new Interval("11", lengths.get("11") - 88000, lengths.get("11")), "3p");
        itm.put(new Interval("12", 0, 60300), "5p");
        itm.put(new Interval("12", lengths.get("12") - 88000, lengths.get("12")), "3p");
        itm.put(new Interval("13", 0, 74413), "5p");
        itm.put(new Interval("13", lengths.get("13") - 102000, lengths.get("13")), "3p");
        itm.put(new Interval("14", 0, 35774), "5p");
        itm.put(new Interval("14", lengths.get("14") - 36000, lengths.get("14")), "3p");

        ExternalAligner ea = new LastzAligner();

        ReferenceSequence rseq;
        while ((rseq = VARS.nextSequence()) != null) {
            String seq = new String(rseq.getBases());
            StringBuilder sb = new StringBuilder();

            List<Set<Interval>> alignments = new ArrayList<Set<Interval>>();

            Map<String, Integer> chrs = new HashMap<String, Integer>();

            for (int i = 0; i <= seq.length() - KMER_SIZE; i++) {
                String fw = seq.substring(i, i + KMER_SIZE);
                String rc = SequenceUtils.reverseComplement(fw);

                Set<Interval> intervals = new HashSet<Interval>();
                intervals.addAll(kl.findKmer(fw));
                intervals.addAll(kl.findKmer(rc));

                if (intervals.size() == 1) {
                    for (Interval interval : intervals) {
                        ContainerUtils.increment(chrs, interval.getSequence());
                    }
                }

                if (intervals.size() == 0) {
                    sb.append(".");
                } else if (intervals.size() == 1) {
                    sb.append("1");
                } else {
                    sb.append("+");
                }

                alignments.add(intervals);
            }

            String mostCommonChr = ContainerUtils.mostCommonKey(chrs);

            StringBuilder newseq = new StringBuilder(seq);
            for (int i = 0; i < sb.length(); i++) {
                if (sb.charAt(i) == '.') {
                    newseq.setCharAt(i, '.');
                }
            }

            for (int i = 0; i < sb.length() - 1; i++) {
                if (sb.charAt(i) != '.' && sb.charAt(i + 1) == '.') {
                    for (int j = i; j < i + KMER_SIZE; j++) {
                        newseq.setCharAt(j, seq.charAt(j));
                    }
                }
            }

            String[] pieces = newseq.toString().split("\\.+");

            String aseq = new String(ifsf.getSequence(mostCommonChr).getBases());
            int lowestPosition = Integer.MAX_VALUE;
            for (String piece : pieces) {
                List<SAMRecord> sas = ea.align(piece, aseq);
                SAMRecord sa = getBestAlignment(sas);

                if (sa != null) {
                    lowestPosition = sa.getAlignmentStart();
                }
            }

            String chr = mostCommonChr
                    .replaceAll("_v3", "")
                    .replaceAll("Pf3D7_", "")
                    .replaceAll("PfHB3_", "")
                    .replaceAll("PfDD2_", "");

            Interval it = new Interval(chr, lowestPosition, lowestPosition);

            String tag = "core";
            if (itm.containsOverlapping(it)) {
                tag = itm.getOverlapping(it).iterator().next();
            }

            log.info("{} {}:{} {} {}", rseq.getName(), mostCommonChr, lowestPosition, tag, aseq.length());
        }
    }

    private int getAlignmentLength(SAMRecord sr) {
        int length = 0;
        for (CigarElement ce : sr.getCigar().getCigarElements()) {
            if (ce.getOperator().equals(CigarOperator.M)) {
                length += ce.getLength();
            }
        }

        return length;
    }

    private SAMRecord getBestAlignment(List<SAMRecord> sas) {
        int bestAlignmentLength = 0;
        SAMRecord sr = null;

        for (SAMRecord sa : sas) {
            int al = getAlignmentLength(sa);

            if (al > bestAlignmentLength) {
                bestAlignmentLength = al;
                sr = sa;
            }
        }

        return sr;
    }
}
