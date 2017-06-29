package uk.ac.ox.well.indiana.commands.playground.assemblies;

import com.google.common.base.Joiner;
import htsjdk.samtools.reference.FastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.util.StringUtil;
import uk.ac.ox.well.indiana.Indiana;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.exceptions.IndianaException;
import uk.ac.ox.well.indiana.utils.io.table.TableReader;
import uk.ac.ox.well.indiana.utils.sequence.SequenceUtils;

import java.io.File;
import java.io.PrintStream;
import java.util.*;
import java.util.regex.Pattern;

public class Contiguate extends Module {
    @Argument(fullName="draft", shortName="d", doc="Draft")
    public FastaSequenceFile DRAFT;

    @Argument(fullName="alignments", shortName="a", doc="MUMmer alignment coordinates")
    public File ALIGNMENTS;

    @Argument(fullName="patternFind", shortName="pf", doc="Find pattern for contig names")
    public ArrayList<String> FIND_PATTERNS;

    @Argument(fullName="patternReplace", shortName="pr", doc="Replacement pattern for contig names")
    public ArrayList<String> REPLACEMENT_PATTERNS;

    @Output
    public PrintStream out;

    @Output(fullName="agpout", shortName="ao", doc="AGP out")
    public PrintStream aout;

    @Override
    public void execute() {
        log.info("Loading draft contigs...");

        Map<String, ReferenceSequence> seqsFwd = new HashMap<>();
        Map<String, ReferenceSequence> seqsRev = new HashMap<>();

        ReferenceSequence seqFwd;
        while ((seqFwd = DRAFT.nextSequence()) != null) {
            String name = seqFwd.getName().split("\\s+")[0];

            ReferenceSequence seqRev = new ReferenceSequence(name, seqFwd.getContigIndex(), SequenceUtils.reverseComplement(seqFwd.getBases()));

            seqsFwd.put(name, seqFwd);
            seqsRev.put(name, seqRev);
        }

        log.info("  {} contigs loaded", seqsFwd.size());

        log.info("Loading alignments...");
        Map<String, Set<ReferenceSequence>> alignments = new LinkedHashMap<>();
        Set<String> used = new HashSet<>();

        TableReader tr = new TableReader(ALIGNMENTS, "S1", "E1", "S2", "E2", "LEN1", "LEN2", "IDY", "LENR", "LENQ", "COVR", "COVQ", "STRANDR", "STRANDQ", "CONR", "CONQ", "REST");
        while (tr.hasNext()) {
            Map<String, String> m = tr.next();

            String rname = m.get("CONR");
            String nname = rname;
            for (int i = 0; i < FIND_PATTERNS.size(); i++) {
                nname = nname.replaceAll(FIND_PATTERNS.get(i), REPLACEMENT_PATTERNS.get(i));
            }

            if (rname.equals(nname) || nname.isEmpty()) {
                throw new IndianaException("Regex for chromosome name transformation didn't work (rname='" + rname + "' nname='" + nname + "' find='" + FIND_PATTERNS + "', replace='" + REPLACEMENT_PATTERNS + "')");
            }

            log.info("  transform {} -> {}", rname, nname);

            String qname = m.get("CONQ");
            boolean positiveStrand = m.get("STRANDQ").equals("1");
            ReferenceSequence seq = positiveStrand ? seqsFwd.get(qname) : seqsRev.get(qname);

            if (!used.contains(qname)) {
                if (!alignments.containsKey(nname)) {
                    alignments.put(nname, new LinkedHashSet<>());
                }
                alignments.get(nname).add(seq);

                used.add(qname);
            }
        }

        log.info("Writing FASTA and AGP files...");

        aout.println("##agp-version\t2.0");

        String curname = null;
        int curoffset = 1;
        int curindex = 1;
        for (String nname : alignments.keySet()) {
            if (curname == null) {
                curname = nname;
                curoffset = 1;
                curindex = 1;
            }

            StringBuilder sb = new StringBuilder();
            for (ReferenceSequence rseq : alignments.get(nname)) {
                if (curindex > 1) {
                    aout.println(Joiner.on("\t").join(nname, curoffset, curoffset + 100 - 1, curindex, "U", 100, "scaffold", "yes", "align_xgenus"));
                    sb.append(StringUtil.repeatCharNTimes('N', 100));

                    curoffset += 100;
                    curindex++;
                }

                aout.println(Joiner.on("\t").join(nname, curoffset, curoffset + rseq.length() - 1, curindex, "W", rseq.getName(), 1, rseq.length(), "+"));
                sb.append(rseq.getBaseString());

                curoffset += rseq.length();
                curindex++;
            }

            out.println(">" + nname);
            out.println(sb.toString());
        }
    }
}
