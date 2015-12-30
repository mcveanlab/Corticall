package uk.ac.ox.well.indiana.commands.gg;

import com.google.common.base.Joiner;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.IntervalTreeMap;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.alignment.pairwise.BwaAligner;
import uk.ac.ox.well.indiana.utils.alignment.pairwise.ExternalAligner;
import uk.ac.ox.well.indiana.utils.alignment.pairwise.LastzAligner;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.containers.ContainerUtils;
import uk.ac.ox.well.indiana.utils.exceptions.IndianaException;
import uk.ac.ox.well.indiana.utils.io.gff.GFF3;
import uk.ac.ox.well.indiana.utils.io.gff.GFF3Record;
import uk.ac.ox.well.indiana.utils.io.table.TableReader;
import uk.ac.ox.well.indiana.utils.io.utils.LineReader;
import uk.ac.ox.well.indiana.utils.sequence.SequenceUtils;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.net.URLDecoder;
import java.net.URLEncoder;
import java.util.*;

public class TransferAnnotations extends Module {
    @Argument(fullName="reference", shortName="r", doc="Reference")
    public ArrayList<IndexedFastaSequenceFile> REFS;

    @Argument(fullName="gff", shortName="g", doc="Reference")
    public ArrayList<GFF3> GFFS;

    @Argument(fullName="descriptions", shortName="d", doc="Descriptions GFF", required=false)
    public GFF3 DESCRIPTIONS;

    @Argument(fullName="summaries", shortName="s", doc="Summaries txt", required=false)
    public File SUMMARIES;

    @Argument(fullName="target", shortName="t", doc="Target reference")
    public File TARGET;

    @Output
    public PrintStream out;

    private Map<String, String> loadDescriptions() {
        Map<String, String> desc = new HashMap<String, String>();

        if (DESCRIPTIONS != null) {
            for (GFF3Record gr : DESCRIPTIONS) {
                if (gr.getType().equals("gene")) {
                    desc.put(gr.getAttribute("ID"), gr.getAttribute("description"));
                }
            }
        }

        if (SUMMARIES != null) {
            LineReader lr = new LineReader(SUMMARIES);
            lr.getNextRecord();

            while(lr.hasNext()) {
                String[] pieces = lr.getNextRecord().split("\t");

                String geneName = pieces[0];
                String description = URLEncoder.encode(pieces[7]);

                if (!desc.containsKey(geneName)) {
                    desc.put(geneName, description);
                }
            }
        }

        return desc;
    }

    @Override
    public void execute() {
        Map<String, String> desc = loadDescriptions();
        List<ReferenceSequence> exons = new ArrayList<ReferenceSequence>();

        ExternalAligner la = new BwaAligner();
        //ExternalAligner la = new LastzAligner();

        log.info("Loading exons...");

        for (int i = 0; i < REFS.size(); i++) {
            IndexedFastaSequenceFile ref = REFS.get(i);
            GFF3 gff = GFFS.get(i);

            int index = 0;

            for (GFF3Record gr : gff) {
                if (gr.getType().equals("exon")) {
                    String id = gr.getAttribute("ID");
                    boolean fw = gr.getStrand().equals(GFF3Record.Strand.POSITIVE);
                    String chr = gr.getSeqid();
                    int start = gr.getStart();
                    int stop = gr.getEnd();

                    String seq = new String(ref.getSubsequenceAt(chr, start, stop).getBases());

                    if (!fw) {
                        seq = SequenceUtils.reverseComplement(seq);
                    }

                    ReferenceSequence exon = new ReferenceSequence(id, index, seq.getBytes());
                    exons.add(exon);

                    index++;
                }
            }

            log.info("  {}: {} exons", i, exons.size());
        }

        log.info("Aligning...");
        List<SAMRecord> recs = la.align(exons, TARGET);

        Collections.sort(recs, new Comparator<SAMRecord>() {
            @Override
            public int compare(SAMRecord o1, SAMRecord o2) {
                if (o1.getReferenceIndex().equals(o2.getReferenceIndex())) {
                    return o1.getAlignmentStart() < o2.getAlignmentStart() ? -1 : 1;
                }

                if (o1.getReferenceIndex() < o2.getReferenceIndex()) {
                    return -1;
                } else if (o1.getReferenceIndex() > o2.getReferenceIndex()) {
                    return 1;
                }

                if (o1.getReadName().contains("PF3D7")) {
                    return -1;
                } else if (o2.getReadName().contains("PF3D7")) {
                    return 1;
                }

                return 0;
            }
        });

        Map<String, Integer> numAlignments = new HashMap<String, Integer>();
        for (SAMRecord rec : recs) {
            ContainerUtils.increment(numAlignments, rec.getReadName());
        }

        out.println("##gff-version 3");
        out.println("# feature-ontology so.obo");
        out.println("# attribute-ontology gff3_attributes.obo");

        GFF3Record currentRecord = null;

        Map<String, GFF3Record> genes = new HashMap<String, GFF3Record>();

        for (SAMRecord rec : recs) {
            if (!rec.getReadUnmappedFlag() && numAlignments.get(rec.getReadName()) == 1) {
                String geneName = rec.getReadName();
                geneName = geneName.replaceAll("exon[_:]", "");
                geneName = geneName.replaceAll("T\\d+:\\d+", "");
                geneName = geneName.replaceAll("-\\d+", "");

                if (currentRecord == null) {
                    currentRecord = new GFF3Record();
                    currentRecord.setSeqid(rec.getReferenceName());
                    currentRecord.setSource("combined");
                    currentRecord.setType("exon");
                    currentRecord.setStart(rec.getAlignmentStart());
                    currentRecord.setEnd(rec.getAlignmentEnd());
                    currentRecord.setScore(".");
                    currentRecord.setStrand(rec.getReadNegativeStrandFlag() ? GFF3Record.Strand.NEGATIVE : GFF3Record.Strand.POSITIVE);
                    currentRecord.setPhase(".");
                    currentRecord.setAttribute("ID", rec.getReadName());
                    currentRecord.setAttribute("names", rec.getReadName());

                    GFF3Record geneRecord = new GFF3Record();
                    geneRecord.setSeqid(rec.getReferenceName());
                    geneRecord.setSource("combined");
                    geneRecord.setType("gene");
                    geneRecord.setStart(rec.getAlignmentStart());
                    geneRecord.setEnd(rec.getAlignmentEnd());
                    geneRecord.setScore(".");
                    geneRecord.setStrand(rec.getReadNegativeStrandFlag() ? GFF3Record.Strand.NEGATIVE : GFF3Record.Strand.POSITIVE);
                    geneRecord.setPhase(".");
                    geneRecord.setAttribute("ID", geneName);

                    genes.put(geneName, geneRecord);
                } else {
                    if (currentRecord.getStart() == rec.getAlignmentStart()) {
                        if (rec.getAlignmentEnd() > currentRecord.getEnd()) {
                            currentRecord.setEnd(rec.getAlignmentEnd());
                        }

                        if (genes.containsKey(geneName) && currentRecord.getSeqid().equals(genes.get(geneName).getSeqid())) {
                            if (currentRecord.getStart() < genes.get(geneName).getStart()) {
                                genes.get(geneName).setStart(currentRecord.getStart());
                            }

                            if (currentRecord.getEnd() > genes.get(geneName).getEnd()) {
                                genes.get(geneName).setEnd(currentRecord.getEnd());
                            }
                        }

                        currentRecord.setAttribute("names", currentRecord.getAttribute("names") + "," + rec.getReadName());
                    } else {
                        out.println(currentRecord);

                        currentRecord = new GFF3Record();
                        currentRecord.setSeqid(rec.getReferenceName());
                        currentRecord.setSource("combined");
                        currentRecord.setType("exon");
                        currentRecord.setStart(rec.getAlignmentStart());
                        currentRecord.setEnd(rec.getAlignmentEnd());
                        currentRecord.setScore(".");
                        currentRecord.setStrand(rec.getReadNegativeStrandFlag() ? GFF3Record.Strand.NEGATIVE : GFF3Record.Strand.POSITIVE);
                        currentRecord.setPhase(".");
                        currentRecord.setAttribute("ID", rec.getReadName());
                        currentRecord.setAttribute("names", rec.getReadName());

                        if (!genes.containsKey(geneName)) {
                            GFF3Record geneRecord = new GFF3Record();
                            geneRecord.setSeqid(rec.getReferenceName());
                            geneRecord.setSource("combined");
                            geneRecord.setType("gene");
                            geneRecord.setStart(rec.getAlignmentStart());
                            geneRecord.setEnd(rec.getAlignmentEnd());
                            geneRecord.setScore(".");
                            geneRecord.setStrand(rec.getReadNegativeStrandFlag() ? GFF3Record.Strand.NEGATIVE : GFF3Record.Strand.POSITIVE);
                            geneRecord.setPhase(".");
                            geneRecord.setAttribute("ID", geneName);

                            genes.put(geneName, geneRecord);
                        } else {
                            if (currentRecord.getSeqid().equals(genes.get(geneName).getSeqid())) {
                                if (currentRecord.getStart() < genes.get(geneName).getStart()) {
                                    genes.get(geneName).setStart(currentRecord.getStart());
                                }

                                if (currentRecord.getEnd() > genes.get(geneName).getEnd()) {
                                    genes.get(geneName).setEnd(currentRecord.getEnd());
                                }
                            }
                        }
                    }
                }
            }
        }

        out.println(currentRecord);

        for (String geneName : genes.keySet()) {
            if (desc.containsKey(geneName)) {
                genes.get(geneName).setAttribute("description", desc.get(geneName));
            } else {
                genes.get(geneName).setAttribute("description", "predicted");
            }

            if (genes.get(geneName).getInterval().length() < 60000) {
                out.println(genes.get(geneName));
            }
        }
    }
}
