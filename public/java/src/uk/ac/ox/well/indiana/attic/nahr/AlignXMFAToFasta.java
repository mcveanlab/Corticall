package uk.ac.ox.well.indiana.attic.nahr;

import net.sf.picard.reference.IndexedFastaSequenceFile;
import net.sf.picard.reference.ReferenceSequence;
import net.sf.picard.util.Interval;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.alignment.exact.ExactLookup;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.xmfa.XMFARecord;
import uk.ac.ox.well.indiana.utils.io.xmfa.XMFASequenceFile;
import uk.ac.ox.well.indiana.utils.sequence.SequenceUtils;

import java.io.PrintStream;
import java.util.HashMap;
import java.util.Map;

public class AlignXMFAToFasta extends Module {
    @Argument(fullName="xfma", shortName="xf", doc="XMFA file")
    public XMFASequenceFile XF;

    @Argument(fullName="fasta", shortName="f", doc="FASTA file")
    public HashMap<String, IndexedFastaSequenceFile> FASTA;

    @Argument(fullName="xmfaName", shortName="xn", doc="XMFA name")
    public String XMFA_NAME;

    @Argument(fullName="fastaName", shortName="fn", doc="FASTA name")
    public String FASTA_NAME;

    @Output
    public PrintStream out;

    private class Loc {
        public String contigName;
        public int start;
        public int stop;
        public boolean isFw;

        @Override
        public String toString() {
            return contigName + ":" + start + "-" + stop;
        }

        @Override
        public int hashCode() {
            return contigName.hashCode() + Integer.valueOf(start).hashCode() - Integer.valueOf(stop).hashCode() + Boolean.valueOf(isFw).hashCode();
        }

        @Override
        public boolean equals(Object obj) {
            Loc l = (Loc) obj;

            return (contigName.equals(l.contigName) && start == l.start && stop == l.stop && isFw == l.isFw);
        }
    }

    private Loc findContig(String fw, Map<String, String> seqs) {
        fw = fw.replaceAll("-", "");
        String rc = SequenceUtils.reverseComplement(fw);

        Loc l = new Loc();

        for (String contigName : seqs.keySet()) {
            l.contigName = contigName;

            String contig = seqs.get(l.contigName);
            l.isFw = true;

            if (contig.contains(fw)) {
                l.start = contig.indexOf(fw);
            } else if (contig.contains(rc)) {
                l.start = contig.indexOf(rc);
                l.isFw = false;
            }

            if (l.start > 0) {
                l.stop = l.start + fw.length();

                break;
            }
        }

        return l;
    }

    private Map<String, String> loadReference(IndexedFastaSequenceFile fasta) {
        Map<String, String> seqs = new HashMap<String, String>();
        ReferenceSequence rseq;
        while ((rseq = fasta.nextSequence()) != null) {
            String[] names = rseq.getName().split("\\s+");
            String seq = new String(rseq.getBases());

            seqs.put(names[0], seq);
        }

        return seqs;
    }

    @Override
    public void execute() {
        Map<Loc, Loc> coordinateMap = new HashMap<Loc, Loc>();

        log.info("Hashing {}...", FASTA_NAME);
        ExactLookup seqRef = new ExactLookup(FASTA.get(FASTA_NAME));

        log.info("Hashing {}...", XMFA_NAME);
        ExactLookup seqAlt = new ExactLookup(FASTA.get(XMFA_NAME));

        log.info("Processing XMFA file...");

        int index = 0;
        for (XMFARecord xr : XF) {
            if (xr.containsKey(FASTA_NAME) && xr.containsKey(XMFA_NAME)) {
                String fwref = new String(xr.get(FASTA_NAME).getBases());
                Interval lref = seqRef.find(fwref);
                //Loc lref = findContig(fwref, seqRef);

                String fwalt = new String(xr.get(XMFA_NAME).getBases());
                Interval lalt = seqAlt.find(fwalt);
                //Loc lalt = findContig(fwalt, seqAlt);

                log.info("ref: {}, alt: {}", lref, lalt);

                /*
                // Get approximate mapping between ref coordinates and alt coordinates
                for (int i = 0; i < fwref.replaceAll("-", "").length(); i++) {
                    Loc locref = new Loc();
                    locref.contigName = lref.contigName;
                    locref.start = lref.start + i;
                    locref.stop = locref.start + 1;
                    locref.isFw = true;

                    Loc localt = new Loc();
                    localt.contigName = lalt.contigName;
                    localt.start = lalt.start + i;
                    localt.stop = localt.start + 1;
                    localt.isFw = true;

                    coordinateMap.put(locref, localt);
                }

                out.println(lref.contigName + " " + lref.start + " " + lref.stop + " color=black");
                */
            }

            if (index % (XF.getNumRecords()/5) == 0) {
                log.info("  processed {}/{} records", index, XF.getNumRecords());
            }
            index++;
        }
    }
}
