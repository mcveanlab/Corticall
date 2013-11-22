package uk.ac.ox.well.indiana.analyses.roi;

import com.google.common.base.Joiner;
import net.sf.picard.reference.IndexedFastaSequenceFile;
import uk.ac.ox.well.indiana.tools.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.gff.GFF3;
import uk.ac.ox.well.indiana.utils.io.gff.GFF3Record;
import uk.ac.ox.well.indiana.utils.io.table.TableReader;
import uk.ac.ox.well.indiana.utils.sequence.SequenceUtils;

import java.io.File;
import java.io.PrintStream;
import java.util.*;

public class ExtractRegionOfInterest extends Module {
    @Argument(fullName="domainName", shortName="dn", doc="Domain name")
    public String DOMAIN_NAME;

    @Argument(fullName="reference", shortName="R", doc="Reference fasta")
    public IndexedFastaSequenceFile REFERENCE;

    @Argument(fullName="gff", shortName="gff", doc="GFF file")
    public GFF3 GFF;

    @Argument(fullName="domains", shortName="domains", doc="Domains file")
    public File DOMAINS;

    @Argument(fullName="translate", shortName="t", doc="Translate genomic sequence into amino acids")
    public Boolean TRANSLATE = false;

    @Output
    public PrintStream out;

    @Override
    public void execute() {
        Map<String, Map<String, String>> domains = new HashMap<String, Map<String, String>>();

        TableReader tr = new TableReader(DOMAINS);
        for (Map<String, String> te : tr) {
            if (te.get("domain").contains(DOMAIN_NAME)) {
                domains.put(te.get("newid"), te);
            }
        }

        for (GFF3Record gr : GFF) {
            if (gr.getType().equals("gene") && domains.containsKey(gr.getAttribute("ID"))) {
                List<GFF3Record> exons = new ArrayList<GFF3Record>(GFF3.getType("exon", GFF.getOverlapping(gr)));

                Collections.sort(exons, new Comparator<GFF3Record>() {
                    @Override
                    public int compare(GFF3Record o1, GFF3Record o2) {
                        if (o1.getStart() < o2.getStart()) {
                            return -1;
                        } else if (o1.getStart() > o2.getStart()) {
                            return 1;
                        }

                        return 0;
                    }
                });

                StringBuilder cdsb = new StringBuilder();
                for (GFF3Record exon : exons) {
                    String seq = new String(REFERENCE.getSubsequenceAt(exon.getSeqid(), exon.getStart(), exon.getEnd()).getBases());

                    cdsb.append(seq);
                }

                String cds = cdsb.toString();
                if (gr.getStrand() == GFF3Record.Strand.NEGATIVE) {
                    cds = SequenceUtils.reverseComplement(cds);
                }

                int dstart = Integer.valueOf(domains.get(gr.getAttribute("ID")).get("aa_start"));
                int dend = Integer.valueOf(domains.get(gr.getAttribute("ID")).get("aa_end"));

                StringBuilder domainSeq = new StringBuilder();

                if (TRANSLATE) {
                    StringBuilder aaSeqB = new StringBuilder();

                    for (int i = 0; i < cds.length() - 3; i += 3) {
                        String codon = cds.substring(i, i + 3);
                        String aa = SequenceUtils.codonToAminoAcid(codon);

                        aaSeqB.append(aa);
                    }

                    for (int i = dstart; i <= dend; i++) {
                        domainSeq.append(aaSeqB.charAt(i));
                    }
                } else {
                    short[] aa = new short[cds.length()];
                    short aaPos = 1;
                    for (int i = 0, cp = 0; i < cds.length(); i++, cp++) {
                        if (cp > 2) {
                            cp = 0;
                            aaPos++;
                        }

                        aa[i] = aaPos;
                    }

                    for (int i = 0; i < cds.length(); i++) {
                        if (aa[i] >= dstart && aa[i] <= dend) {
                            domainSeq.append(cds.charAt(i));
                        }
                    }
                }

                out.println(">" + gr.getAttribute("ID") + "_" + domains.get(gr.getAttribute("ID")).get("domain"));
                out.println(domainSeq);
            }
        }
    }
}
