package uk.ac.ox.well.indiana.attic.analyses.nahr;

import net.sf.picard.reference.FastaSequenceFile;
import net.sf.picard.reference.IndexedFastaSequenceFile;
import net.sf.picard.reference.ReferenceSequence;
import net.sf.picard.util.Interval;
import net.sf.picard.util.IntervalTreeMap;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.table.TableReader;
import uk.ac.ox.well.indiana.utils.io.table.TableWriter;
import uk.ac.ox.well.indiana.utils.io.utils.LineReader;
import uk.ac.ox.well.indiana.utils.sequence.SequenceUtils;

import java.io.File;
import java.io.PrintStream;
import java.util.*;

public class FindRecombinationEventsInVCF extends Module {
    @Argument(fullName="vcf", shortName="vcf", doc="VCF file")
    public File VCF;

    @Argument(fullName="ref", shortName="r", doc="Reference FASTA")
    public IndexedFastaSequenceFile REF;

    @Argument(fullName="recombLoci", shortName="rl", doc="Recomb loci table")
    public File RECOMB_LOCI;

    @Argument(fullName="contigs", shortName="c", doc="Contigs")
    public ArrayList<File> CONTIGS;

    @Output
    public PrintStream out;

    private IntervalTreeMap<Map<String, String>> loadVCF() {
        IntervalTreeMap<Map<String, String>> vcfRecords = new IntervalTreeMap<Map<String, String>>();

        LineReader lr = new LineReader(VCF);

        String[] header = null;
        String line;
        while ((line = lr.getNextRecord()) != null) {
            if (line.startsWith("##")) {
                // do nothing
            } else if (line.startsWith("#")) {
                // header
                header = line.split("\t");
                for (int i = 9; i < header.length; i++) {
                    String[] names = header[i].split("/");
                    header[i] = names[1];
                }
                header[0] = header[0].replaceAll("#", "");
            } else {
                if (header != null) {
                    String[] fields = line.split("\t");

                    for (int i = 9; i < fields.length; i++) {
                        String[] values = fields[i].split(":");
                        fields[i] = values[0];
                    }

                    Map<String, String> entry = new HashMap<String, String>();
                    for (int i = 0; i < fields.length; i++) {
                        entry.put(header[i], fields[i]);
                    }

                    if (entry.get("FILTER").equals("PASS")) {
                        Interval interval = new Interval(entry.get("CHROM"), Integer.valueOf(entry.get("POS")), Integer.valueOf(entry.get("POS")));

                        vcfRecords.put(interval, entry);
                    }
                }
            }
        }

        return vcfRecords;
    }

    private List<String> getSampleNames() {
        List<String> sampleNames = new ArrayList<String>();

        LineReader lr = new LineReader(VCF);
        String line;
        while ((line = lr.getNextRecord()) != null) {
            if (line.startsWith("##")) {
                // do nothing
            } else if (line.startsWith("#")) {
                String[] header = line.split("\t");
                for (int i = 9; i < header.length; i++) {
                    String[] names = header[i].split("/");
                    header[i] = names[1];
                }
                header[0] = header[0].replaceAll("#", "");

                for (int i = 9; i < header.length; i++) {
                    sampleNames.add(header[i]);
                }
            }
        }

        return sampleNames;
    }

    private Map<String, Set<Interval>> loadRecombTable() {
        TableReader tr = new TableReader(RECOMB_LOCI);

        Map<String, Set<Interval>> loci = new HashMap<String, Set<Interval>>();

        for (Map<String, String> te : tr) {
            String[] sampleFields = te.get("sample").split("/");
            String sampleName = sampleFields[1];

            if (!loci.containsKey(sampleName)) {
                loci.put(sampleName, new TreeSet<Interval>());
            }

            Interval interval = new Interval(te.get("chrom"), Integer.valueOf(te.get("pos_max")) - 20, Integer.valueOf(te.get("pos_max")) + 21);
            loci.get(sampleName).add(interval);
        }

        return loci;
    }

    @Override
    public void execute() {
        log.info("Loading VCF...");
        IntervalTreeMap<Map<String, String>> vcf = loadVCF();
        List<String> sampleNames = getSampleNames();

        log.info("Loading recombination events...");
        Map<String, Set<Interval>> recombTable = loadRecombTable();

        log.info("Extracting diagnostic kmers from HR sites...");
        Map<String, List<String>> recombKmers = new HashMap<String, List<String>>();
        Map<String, List<Interval>> recombSites = new HashMap<String, List<Interval>>();

        for (String sampleName : sampleNames) {
            log.debug("recombs: {} {}", sampleName, recombTable.get(sampleName));

            if (recombTable.containsKey(sampleName)) {
                recombKmers.put(sampleName, new ArrayList<String>());
                recombSites.put(sampleName, new ArrayList<Interval>());

                for (Interval interval : recombTable.get(sampleName)) {
                    Interval vinterval = interval;

                    boolean found = false;
                    for (int j = vinterval.getEnd() - 21; j > 0 && !found; j--) {
                        Interval jinterval = new Interval(interval.getSequence(), j, j);

                        Collection<Map<String, String>> records = vcf.getOverlapping(jinterval);
                        for (Map<String, String> record : records) {
                            String allele = record.get(sampleName);

                            log.debug("  {}", jinterval);

                            if (!allele.equals("0") && !allele.equals("./.")) {
                                found = true;

                                if (j != vinterval.getStart()) {
                                    vinterval = new Interval(jinterval.getSequence(), jinterval.getStart() - 20, jinterval.getStart() + 21);
                                }
                            }
                        }
                    }

                    String refhap = new String(REF.getSubsequenceAt(vinterval.getSequence(), vinterval.getStart(), vinterval.getEnd()).getBases());
                    StringBuilder modHapBuilder = new StringBuilder();

                    Map<Interval, String> modRefAlleles = new TreeMap<Interval, String>();
                    Map<Interval, String> modAltAlleles = new TreeMap<Interval, String>();

                    Collection<Map<String, String>> records = vcf.getOverlapping(vinterval);
                    for (Map<String, String> record : records) {
                        if (!record.get(sampleName).equals("./.") && !record.get(sampleName).equals("0")) {
                            Interval rinterval = new Interval(record.get("CHROM"), Integer.valueOf(record.get("POS")), Integer.valueOf(record.get("POS")));

                            String refAllele = record.get("REF");
                            String[] altAlleles = record.get("ALT").split(",");
                            int altIndex = Integer.valueOf(record.get(sampleName)) - 1;

                            String altAllele = altAlleles[altIndex];

                            modRefAlleles.put(rinterval, refAllele);
                            modAltAlleles.put(rinterval, altAllele);
                        }
                    }

                    for (int i = 0; i < refhap.length(); i++) {
                        Interval rinterval = new Interval(vinterval.getSequence(), vinterval.getStart() + i, vinterval.getStart() + i);

                        if (!modRefAlleles.containsKey(rinterval)) {
                            modHapBuilder.append(refhap.charAt(i));
                        } else {
                            String refAllele = modRefAlleles.get(rinterval);
                            String altAllele = modAltAlleles.get(rinterval);
                            modHapBuilder.append(altAllele);
                            i += refAllele.length() - 1;
                        }
                    }

                    log.debug("{} {} {} {}", interval, vinterval, modRefAlleles, modAltAlleles);
                    log.debug("ref: {}", refhap);
                    log.debug("alt: {}", modHapBuilder.toString());

                    recombKmers.get(sampleName).add(modHapBuilder.toString());
                    recombSites.get(sampleName).add(vinterval);
                }
            }
        }

        log.info("Finding contigs that contain diagnostic kmers...");
        TableWriter tw = new TableWriter(out);
        for (File contigsFile : CONTIGS) {
            String sampleName = contigsFile.getName().replace(".contigs.unique.fasta", "");
            log.info("  {}", sampleName);

            FastaSequenceFile contigs = new FastaSequenceFile(contigsFile, true);

            ReferenceSequence rseq;
            while ((rseq = contigs.nextSequence()) != null) {
                String seq = new String(rseq.getBases());

                List<String> altHaps = recombKmers.get(sampleName);
                List<Interval> altInts = recombSites.get(sampleName);

                if (recombKmers.containsKey(sampleName)) {
                    for (int i = 0; i < altHaps.size(); i++) {
                        Interval altInt = altInts.get(i);
                        String altHapFw = altHaps.get(i);
                        String altHapRc = SequenceUtils.reverseComplement(altHapFw);

                        if (seq.contains(altHapFw) || seq.contains(altHapRc)) {
                            Map<String, String> entry = new LinkedHashMap<String, String>();

                            entry.put("sampleName", sampleName);
                            entry.put("contigName", rseq.getName());
                            entry.put("locus", altInt.toString());
                            entry.put("diagnosticKmer", SequenceUtils.alphanumericallyLowestOrientation(altHapFw));

                            tw.addEntry(entry);
                        }
                    }
                }
            }
        }

    }
}
