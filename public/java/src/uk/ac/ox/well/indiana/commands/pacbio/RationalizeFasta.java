package uk.ac.ox.well.indiana.commands.pacbio;

import com.google.common.base.Joiner;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.reference.FastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.alignment.pairwise.BwaAligner;
import uk.ac.ox.well.indiana.utils.alignment.pairwise.ExternalAligner;
import uk.ac.ox.well.indiana.utils.alignment.pairwise.LastzAligner;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.table.TableReader;
import uk.ac.ox.well.indiana.utils.sequence.SequenceUtils;

import java.io.File;
import java.io.PrintStream;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class RationalizeFasta extends Module {
    @Argument(fullName="query", shortName="q", doc="Query")
    public File QUERY;

    @Argument(fullName="target", shortName="t", doc="Target")
    public File TARGET;

    @Argument(fullName="Mapping", shortName="m", doc="Mapping")
    public File MAPPING;

    @Output
    public PrintStream out;

    private Map<String, ReferenceSequence> loadRef(File fasta) {
        Map<String, ReferenceSequence> contigs = new HashMap<String, ReferenceSequence>();

        FastaSequenceFile f = new FastaSequenceFile(fasta, true);
        ReferenceSequence rseq;

        while ((rseq = f.nextSequence()) != null) {
            contigs.put(rseq.getName(), rseq);
        }

        return contigs;
    }

    @Override
    public void execute() {
        Map<String, ReferenceSequence> qref = loadRef(QUERY);
        Map<String, ReferenceSequence> tref = loadRef(TARGET);

        //TableReader tr = new TableReader(MAPPING, "refchr", "altchrs", "newchr");
        TableReader tr = new TableReader(MAPPING);

        for (Map<String, String> te : tr) {
            int counter = 0;

            log.info("{} {}", te.get("refchr"), Joiner.on(", ").withKeyValueSeparator("=").join(te));

            for (String altchr : te.get("altchrs").split(",")) {
                boolean rc = altchr.contains("*");
                altchr = altchr.replaceAll("\\*", "");

                String query = new String(qref.get(altchr).getBases());

                if (rc) {
                    query = SequenceUtils.reverseComplement(query);
                }

                out.println(">" + te.get("newchr") + "_0" + counter);
                out.println(query);

                counter++;

                /*
                String refchr = te.get("refchr");
                String target = new String(tref.get(refchr).getBases());

                log.info("{} {} {}", refchr, altchr, rc);

                ExternalAligner la = new LastzAligner();
                List<SAMRecord> srs = la.align(query, target);

                for (SAMRecord sr : srs) {
                    log.info("  {}:{}-{}", sr.getReferenceName(), sr.getAlignmentStart(), sr.getAlignmentEnd());
                }
                */
            }
        }
    }
}
