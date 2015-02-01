package uk.ac.ox.well.indiana.attic.analyses.nahr;

import htsjdk.samtools.*;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.util.StringUtil;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.table.TableReader;
import uk.ac.ox.well.indiana.utils.sequence.SequenceUtils;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.Map;

public class ExtractRefSequenceForContig extends Module {
    @Argument(fullName="metrics", shortName="m", doc="Metrics file")
    public File METRICS;

    @Argument(fullName="reference", shortName="r", doc="Reference file")
    public IndexedFastaSequenceFile REF;

    @Output
    public File out;

    private Cigar cigarStringToInvertedCigar(String cs) {
        List<CigarElement> ces = new ArrayList<CigarElement>();

        int start = 0;
        int stop = 1;

        while (stop < cs.length()) {
            char c = cs.charAt(stop);

            if (c == 'M' || c == 'S' || c == 'H' || c == 'I' || c == 'D') {
                int length = Integer.valueOf(cs.substring(start, stop));
                CigarOperator co = CigarOperator.characterToEnum(c);

                if (c == 'D') { co = CigarOperator.I; }
                if (c == 'I') { co = CigarOperator.D; }

                ces.add(new CigarElement(length, co));

                start = stop + 1;
            }

            stop++;
        }

        return new Cigar(ces);
    }

    private String formatString(String seq, Cigar cigar) {
        StringBuilder sb = new StringBuilder();

        int offset = 0;

        for (CigarElement ce : cigar.getCigarElements()) {
            String piece = "";

            switch (ce.getOperator()) {
                case M:
                case S:
                    piece = seq.substring(offset, offset + ce.getLength());
                    offset += ce.getLength();
                    break;
                case I:
                    piece = StringUtil.repeatCharNTimes('-', ce.getLength());
                    offset += ce.getLength();
                    break;
                case D:
                    break;
            }

            sb.append(piece);
        }

        return sb.toString();
    }

    @Override
    public void execute() {
        TableReader tr = new TableReader(METRICS);

        Map<String, String> firstRecord = tr.iterator().next();

        SAMReadGroupRecord srgr = new SAMReadGroupRecord(firstRecord.get("sampleName") + "." + firstRecord.get("accession"));
        srgr.setSample(firstRecord.get("sampleName"));

        SAMFileHeader sfh = new SAMFileHeader();
        sfh.addReadGroup(srgr);
        sfh.setSortOrder(SAMFileHeader.SortOrder.coordinate);

        List<SAMSequenceRecord> ssrs = new ArrayList<SAMSequenceRecord>();
        for (Map<String, String> te : tr) {
            String contigName = te.get("contigName").replaceAll("_", "-");
            int contigLength = te.get("seq").length();

            SAMSequenceRecord ssr = new SAMSequenceRecord(contigName, contigLength);
            ssrs.add(ssr);
        }

        SAMSequenceDictionary ssd = new SAMSequenceDictionary(ssrs);
        sfh.setSequenceDictionary(ssd);

        SAMFileWriterFactory sfwf = new SAMFileWriterFactory();
        sfwf.setCreateIndex(true);

        SAMFileWriter sfw = sfwf.makeBAMWriter(sfh, false, out);

        log.info("Processing contigs...");
        int numContigs = 0;
        for (Map<String, String> te : tr) {
            if (numContigs % (tr.size()/10) == 0) {
                log.info("  {}/{} contigs", numContigs, tr.size());
            }
            numContigs++;

            String canonicalLocus = te.get("canonicalLocus");
            String canonicalCigar = te.get("cigarCanonical");
            String canonicalIsRc  = te.get("isRcCanonical");

            if (!canonicalLocus.isEmpty() && !canonicalLocus.equals("*:0-0") && !canonicalCigar.isEmpty() && !canonicalIsRc.isEmpty()) {
                Cigar invertedCigar = cigarStringToInvertedCigar(canonicalCigar);

                String[] pieces = canonicalLocus.split("[:-]");
                int start = Integer.valueOf(pieces[1]);
                int stop = Integer.valueOf(pieces[2]);
                int prePad = 0;
                int postPad = 0;

                for (int i = 0; i < invertedCigar.getCigarElements().size(); i++) {
                    CigarElement ce = invertedCigar.getCigarElement(i);

                    if (ce.getOperator().equals(CigarOperator.S)) {
                        if (i == 0) {
                            prePad += ce.getLength();
                        } else {
                            postPad += ce.getLength();
                        }
                    }
                }

                int alignmentStart = prePad + 1;

                while (start > 1 && prePad > 0) {
                    start--;
                    prePad--;
                }

                while (stop < REF.getSequence(pieces[0]).length() && postPad > 0) {
                    stop++;
                    postPad--;
                }

                ReferenceSequence rseq = REF.getSubsequenceAt(pieces[0], start, stop);
                String ref = new String(rseq.getBases());
                if (canonicalIsRc.equals("1")) {
                    ref = SequenceUtils.reverseComplement(ref);
                }
                ref = StringUtil.repeatCharNTimes('N', prePad) + ref + StringUtil.repeatCharNTimes('N', postPad);

                SAMRecord sr = new SAMRecord(sfh);
                sr.setAttribute("RG", srgr.getReadGroupId());
                sr.setReadString(ref);
                sr.setBaseQualityString(StringUtil.repeatCharNTimes('I', ref.length()));
                sr.setReadName(te.get("contigName") + ".ref");
                sr.setReadPairedFlag(false);
                sr.setNotPrimaryAlignmentFlag(false);
                sr.setCigar(invertedCigar);
                sr.setReferenceName(te.get("contigName").replaceAll("_", "-"));
                sr.setAlignmentStart(alignmentStart);
                sr.setMappingQuality(60);

                sfw.addAlignment(sr);
            }
        }

        sfw.close();
    }
}
