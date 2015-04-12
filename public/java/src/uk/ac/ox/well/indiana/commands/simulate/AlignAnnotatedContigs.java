package uk.ac.ox.well.indiana.commands.simulate;

import htsjdk.samtools.SAMFileReader;
import htsjdk.samtools.SAMFileWriter;
import htsjdk.samtools.SAMFileWriterFactory;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.reference.ReferenceSequence;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.alignment.pairwise.LastzAligner;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.table.TableReader;

import java.io.File;
import java.io.PrintStream;
import java.util.*;

public class AlignAnnotatedContigs extends Module {
    @Argument(fullName="metrics", shortName="m", doc="Contig metrics")
    public File METRICS;

    @Argument(fullName="reference0", shortName="r0", doc="Reference sequence 0")
    public File REF0;

    @Argument(fullName="reference1", shortName="r1", doc="Reference sequence 1")
    public File REF1;

    @Argument(fullName="bam0", shortName="b0", doc="Bam 0")
    public SAMFileReader BAM0;

    @Argument(fullName="bam1", shortName="b1", doc="Bam 1")
    public SAMFileReader BAM1;

    @Output(fullName="bamOut0", shortName="bo0", doc="Bam out 0")
    public File BAMOUT0;

    @Output(fullName="bamOut1", shortName="bo1", doc="Bam out 1")
    public File BAMOUT1;

    private char toLowerCase(char a) {
        return String.valueOf(a).toLowerCase().charAt(0);
    }

    private String softMaskSequence(String seq, String kmerOrigin) {
        int kmerSize = seq.length() - kmerOrigin.length() + 1;

        StringBuilder softMaskedSeq = new StringBuilder(seq);

        for (int i = 0; i < kmerOrigin.length(); i++) {
            char code = kmerOrigin.charAt(i);
            if (code == 'R' || code == '.') {
                softMaskedSeq.setCharAt(i + kmerSize - 1, toLowerCase(softMaskedSeq.charAt(i + kmerSize - 1)));
            }
        }

        return softMaskedSeq.toString();
    }

    private boolean shouldRealign(Map<String, String> te, int refIndex) {
        String refLocus = te.get("ref" + refIndex + "Locus");
        boolean isUnaligned = refLocus.equals("*:0-0") || refLocus.equals("NA");
        boolean clippedInRef = te.get("clippedInRef" + refIndex).equals("1");

        //if (te.get("contigName").contains("contig-5999")) {
            //log.info("Hi!");
        //}

        return (isUnaligned || clippedInRef);
    }

    private void replaceAlignments(SAMFileReader bamIn, Map<String, String[]> replacements, File bamOut) {
        SAMFileWriter sfw = new SAMFileWriterFactory().setCreateIndex(true).makeBAMWriter(bamIn.getFileHeader(), false, bamOut);

        Set<String> isWritten = new HashSet<String>();
        for (SAMRecord orig : bamIn) {
            String contigName = orig.getReadName();

            if (!isWritten.contains(contigName)) {
                if (replacements.containsKey(contigName)) {
                    String[] fields = replacements.get(contigName);

                    SAMRecord sr = new SAMRecord(bamIn.getFileHeader());
                    sr.setReadName(fields[0]);
                    sr.setReadNegativeStrandFlag(fields[1].equals("16"));
                    sr.setReferenceName(fields[2]);
                    sr.setAlignmentStart(Integer.valueOf(fields[3]));
                    sr.setMappingQuality(Integer.valueOf(fields[4]));
                    sr.setCigarString(fields[5]);
                    sr.setReadString(fields[9]);
                    sr.setAttribute("XO", orig.getCigarString());

                    sfw.addAlignment(sr);

                    isWritten.add(contigName);
                } else if (!orig.isSecondaryOrSupplementary()) {
                    sfw.addAlignment(orig);

                    isWritten.add(contigName);
                }
            }
        }

        sfw.close();
    }

    @Override
    public void execute() {
        log.info("Examining contigs...");
        Set<ReferenceSequence> queries0 = new HashSet<ReferenceSequence>();
        Set<ReferenceSequence> queries1 = new HashSet<ReferenceSequence>();

        TableReader tr = new TableReader(METRICS);
        for (Map<String, String> te : tr) {
            String contigName = te.get("contigName");
            String seq = te.get("seq");
            String kmerOrigin = te.get("kmerOrigin");
            //String softMaskedSeq = softMaskSequence(seq, kmerOrigin);

            //ReferenceSequence query = new ReferenceSequence(contigName, 0, softMaskedSeq.getBytes());
            ReferenceSequence query = new ReferenceSequence(contigName, 0, seq.getBytes());

            //if (contigName.equals("contig-5999")) {
                if (shouldRealign(te, 0)) { queries0.add(query); }
                if (shouldRealign(te, 1)) { queries1.add(query); }
            //}

            //if (queries0.size() == 10) { break; }
        }
        log.info("  found {} realignment candidates in ref0", queries0.size());
        log.info("  found {} realignment candidates in ref1", queries1.size());

        LastzAligner la = new LastzAligner();
        log.info("Realigning {}/{} contigs to ref0...", queries0.size(), tr.size());
        Map<String, String[]> results0 = la.alignAll(queries0, REF0);
        log.info("  found {} replacements", results0.size());

        log.info("Realigning {}/{} contigs to ref1...", queries1.size(), tr.size());
        Map<String, String[]> results1 = la.alignAll(queries1, REF1);
        log.info("  found {} replacements", results1.size());

        log.info("Writing alignments...");
        replaceAlignments(BAM0, results0, BAMOUT0);
        replaceAlignments(BAM1, results1, BAMOUT1);
    }
}
