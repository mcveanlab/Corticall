package uk.ac.ox.well.indiana.attic.analyses.nahr;

import com.google.common.base.Joiner;
import net.sf.picard.reference.IndexedFastaSequenceFile;
import net.sf.picard.util.Interval;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.alignment.exact.ExactLookup;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.xmfa.XMFARecord;
import uk.ac.ox.well.indiana.utils.io.xmfa.XMFASequenceFile;
import uk.ac.ox.well.indiana.utils.sequence.SequenceUtils;

import java.io.PrintStream;
import java.util.*;

public class AlignXMFAToFasta extends Module {
    @Argument(fullName="xfma", shortName="xf", doc="XMFA file")
    public XMFASequenceFile XF;

    @Argument(fullName="fasta", shortName="f", doc="FASTA file")
    public HashMap<String, IndexedFastaSequenceFile> FASTA;

    @Argument(fullName="bam", shortName="b", doc="BAM file")
    public HashMap<String, SAMFileReader> BAM;

    @Argument(fullName="xmfaName", shortName="xn", doc="XMFA name")
    public String XMFA_NAME;

    @Argument(fullName="fastaName", shortName="fn", doc="FASTA name")
    public String FASTA_NAME;

    @Output
    public PrintStream out;

    private class Loc {
        public String contigName = "";
        public int start = 0;
        //public int stop;
        //public boolean isFw;

        @Override
        public String toString() {
            return contigName + ":" + start;
        }

        @Override
        public int hashCode() {
            int h1 = contigName.hashCode();
            int h2 = Integer.valueOf(start).hashCode();

            return h1 + h2;
        }

        @Override
        public boolean equals(Object obj) {
            Loc l = (Loc) obj;

            return (contigName.equals(l.contigName) && start == l.start);
        }
    }

    @Override
    public void execute() {
        Map<Loc, Loc> coordinateMap = new HashMap<Loc, Loc>();

        ExactLookup seqRef = new ExactLookup(FASTA.get(FASTA_NAME));
        ExactLookup seqAlt = new ExactLookup(FASTA.get(XMFA_NAME));

        log.info("Processing XMFA file...");

        int alt = 0;
        int noalt = 0;

        int index = 0;
        for (XMFARecord xr : XF) {
            if (xr.containsKey(FASTA_NAME) && xr.containsKey(XMFA_NAME)) {
                String fwref = new String(xr.get(FASTA_NAME).getBases());
                Interval lref = seqRef.find(fwref);
                if (lref == null) {
                    lref = seqRef.find(SequenceUtils.reverseComplement(fwref));
                }

                String fwalt = new String(xr.get(XMFA_NAME).getBases());
                Interval lalt = seqAlt.find(fwalt);
                if (lalt == null) {
                    lalt = seqRef.find(SequenceUtils.reverseComplement(fwalt));
                }

                //log.info("ref: {}, alt: {}", lref, lalt);

                if (lalt == null) { noalt++; }
                else { alt++; }

                if (lref != null && lalt != null) {
                    // Get approximate mapping between ref coordinates and alt coordinates
                    for (int i = 0; i < fwref.replaceAll("-", "").length(); i++) {
                        Loc locref = new Loc();
                        locref.contigName = lref.getSequence();
                        locref.start = lref.getStart() + i;

                        Loc localt = new Loc();
                        localt.contigName = lalt.getSequence();
                        localt.start = lalt.getStart() + i;

                        //log.info("ref: {} {}", lref.getName(), )

                        coordinateMap.put(localt, locref);
                    }
                }

                //out.println(lref.contigName + " " + lref.start + " " + lref.stop + " color=black");
            }

            if (index % (XF.getNumRecords()/5) == 0) {
                log.info("  processed {}/{} records", index, XF.getNumRecords());
            }
            index++;
        }

        log.info("Alt: {}, noalt: {}", alt, noalt);

        Map<String, Set<Loc>> contigRef = new HashMap<String, Set<Loc>>();
        Map<String, Set<Loc>> contigAlt = new HashMap<String, Set<Loc>>();

        for (String bamId : BAM.keySet()) {
            for (SAMRecord br : BAM.get(bamId)) {
                String contigName = br.getReadName();

                if (!contigRef.containsKey(contigName)) {
                    contigRef.put(contigName, new HashSet<Loc>());
                }

                if (!contigAlt.containsKey(contigName)) {
                    contigAlt.put(contigName, new HashSet<Loc>());
                }

                Loc loc = new Loc();
                loc.contigName = br.getReferenceName();
                loc.start = br.getAlignmentStart();

                if (bamId.equals("3D7")) {
                    contigRef.get(contigName).add(loc);
                } else if (coordinateMap.containsKey(loc)) {
                    contigAlt.get(contigName).add(coordinateMap.get(loc));
                } else {
                    contigAlt.get(contigName).add(loc);
                }
            }
        }

        for (String contigName : contigRef.keySet()) {
            out.println(contigName + "\t" + Joiner.on(";").join(contigRef.get(contigName)) + "\t" + Joiner.on(";").join(contigAlt.get(contigName)));
        }
    }
}
