package uk.ac.ox.well.indiana.commands.playground.assemblies;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.reference.FastaSequenceFile;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.io.gff.GFF3;
import uk.ac.ox.well.indiana.utils.io.gff.GFF3Record;

import java.util.HashMap;
import java.util.Map;

/**
 * Created by kiran on 30/06/2017.
 */
public class DetermineAssemblyLayout extends Module {
    @Argument(fullName="draft", shortName="d", doc="Reference")
    public FastaSequenceFile DRAFT;

    @Argument(fullName="gff", shortName="g", doc="GFF3 file")
    public GFF3 GFF;

    @Argument(fullName="exons", shortName="e", doc="Aligned exons")
    public SamReader EXONS;

    @Override
    public void execute() {
        Map<String, GFF3Record> grrecs = new HashMap<>();
        Map<String, SAMRecord> srrecs = new HashMap<>();

        for (GFF3Record gr : GFF) {
            if (gr.getType().equals("exon")) {
                grrecs.put(gr.getAttribute("ID"), gr);
            }
        }

        for (SAMRecord sr : EXONS) {
            srrecs.put(sr.getReadName(), sr);
        }

        for (String id : grrecs.keySet()) {
            if (srrecs.containsKey(id)) {
                log.info("gff: {}", grrecs.get(id));
                log.info("sam: {}", srrecs.get(id));
                log.info("--");
            }
        }
    }
}
