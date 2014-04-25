package uk.ac.ox.well.indiana.attic.analyses.nahr;

import net.sf.picard.reference.FastaSequenceFile;
import net.sf.picard.reference.IndexedFastaSequenceFile;
import net.sf.picard.reference.ReferenceSequence;
import net.sf.picard.util.Interval;
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

    @Argument(fullName="parent1", shortName="p1", doc="Parent 1")
    public String PARENT1;

    @Argument(fullName="parent2", shortName="p2", doc="Parent 2")
    public String PARENT2;

    @Argument(fullName="ref", shortName="r", doc="Reference FASTA")
    public IndexedFastaSequenceFile REF;

    @Argument(fullName="metrics", shortName="m", doc="Metrics file")
    public File METRICS;

    @Output
    public PrintStream out;

    private List<Map<String, String>> loadVCF() {
        List<Map<String, String>> vcfRecords = new ArrayList<Map<String, String>>();

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
                        vcfRecords.add(entry);
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

    @Override
    public void execute() {
        log.info("Loading VCF...");
        List<Map<String, String>> vcfRecords = loadVCF();
        log.info("  loaded {} records", vcfRecords.size());

        log.info("Finding recombination events...");
        Map<String, Set<String>> recombSeqs = new HashMap<String, Set<String>>();
        Map<String, Map<String, Interval>> recombLoci = new HashMap<String, Map<String, Interval>>();

        List<String> sampleNames = getSampleNames();
        for (String sampleName : sampleNames) {
            recombSeqs.put(sampleName, new HashSet<String>());
            recombLoci.put(sampleName, new HashMap<String, Interval>());

            if (!sampleName.equals(PARENT1) && !sampleName.equals(PARENT2)) {
                Set<Integer> switchIndices = new HashSet<Integer>();

                for (int i = 2; i < vcfRecords.size(); i++) {
                    Map<String, String> pentry = vcfRecords.get(i-1);
                    Map<String, String> entry = vcfRecords.get(i);

                    String prevCopyingFrom = pentry.get(sampleName).equals(pentry.get(PARENT1)) ? PARENT1 : PARENT2;
                    String prevCopiedAllele = pentry.get(prevCopyingFrom);

                    String curCopyingFrom = entry.get(sampleName).equals(entry.get(PARENT1)) ? PARENT1 : PARENT2;
                    String curCopiedAllele = entry.get(curCopyingFrom);

                    boolean isSwitch = false;

                    //if (pentry.get("CHROM").equals(entry.get("CHROM")) && !pentry.get(sampleName).equals(entry.get(sampleName))) {
                    if (pentry.get("CHROM").equals(entry.get("CHROM")) && !prevCopiedAllele.equals(curCopiedAllele)) {
                        if (pentry.get(sampleName).equals("0")) { switchIndices.add(i); }
                        else { switchIndices.add(i-1); }

                        isSwitch = true;
                    }

                    log.debug("{} {}:{} ::::: {} {} {} {} {} {} {}", sampleName, entry.get("CHROM"), entry.get("POS"), entry.get(PARENT1), entry.get(PARENT2), entry.get(sampleName), isSwitch, i, entry.get("REF"), entry.get("ALT"));
                }

                for (int index : switchIndices) {
                    Map<String, String> pentry = vcfRecords.get(index - 1);
                    Map<String, String> entry = vcfRecords.get(index);
                    int bhapLength = 0;
                    int fhapLength = 0;
                    Interval bInterval = null;
                    Interval interval = new Interval(entry.get("CHROM"), Integer.valueOf(entry.get("POS")), Integer.valueOf(entry.get("POS")));
                    Interval fInterval = null;
                    int bIndex = 0;
                    int fIndex = 0;

                    for (int b = index - 2; b > 0; b--) {
                        Map<String, String> bentry = vcfRecords.get(b);

                        if (bentry.get("CHROM").equals(pentry.get("CHROM"))) {
                            if (bentry.get(sampleName).equals(pentry.get(sampleName))) {
                                bhapLength++;
                                bIndex = b;
                                bInterval = new Interval(bentry.get("CHROM"), Integer.valueOf(bentry.get("POS")), Integer.valueOf(bentry.get("POS")));

                                //log.debug("{} {} {} {} {} {}", index, b, bhapLength, bInterval, bentry.get(sampleName), pentry.get(sampleName));
                            } else {
                                break;
                            }
                        } else {
                            break;
                        }
                    }

                    for (int f = index + 1; f < vcfRecords.size(); f++) {
                        Map<String, String> fentry = vcfRecords.get(f);

                        if (entry.get("CHROM").equals(fentry.get("CHROM"))) {
                            if (fentry.get(sampleName).equals(entry.get(sampleName))) {
                                fhapLength++;
                                fIndex = f;
                                fInterval = new Interval(fentry.get("CHROM"), Integer.valueOf(fentry.get("POS")), Integer.valueOf(fentry.get("POS")));
                            } else {
                                break;
                            }
                        } else {
                            break;
                        }
                    }

                    if (bInterval != null && fInterval != null) {
                        int length = fInterval.getEnd() - bInterval.getStart();

                        if (length > 20000 && bhapLength > 100 && fhapLength > 100) {
                            log.debug("RECOMB ::::: {} {} {} {} {} {} {} {}", sampleName, index, bhapLength, fhapLength, interval, bInterval, fInterval, length);

                            Interval shortInterval = new Interval(interval.getSequence(), interval.getStart() - 20, interval.getStart() + 21);

                            String refHap = new String(REF.getSubsequenceAt(shortInterval.getSequence(), shortInterval.getStart(), shortInterval.getEnd()).getBases());

                            Map<Integer, String> mods = new HashMap<Integer, String>();
                            Map<Integer, String> modsRef = new HashMap<Integer, String>();
                            Map<Integer, String> modsAlt = new HashMap<Integer, String>();

                            for (int q = bIndex; q <= fIndex; q++) {
                                Map<String, String> qentry = vcfRecords.get(q);

                                Interval qInterval = new Interval(qentry.get("CHROM"), Integer.valueOf(qentry.get("POS")), Integer.valueOf(qentry.get("POS")));

                                if (qInterval.intersects(shortInterval) && !qentry.get(sampleName).equals("0") && !qentry.get(sampleName).equals("./.")) {
                                    String[] alts = qentry.get("ALT").split(",");
                                    int altIndex = Integer.valueOf(qentry.get(sampleName)) - 1;
                                    String alt = alts[altIndex];

                                    int offset = Integer.valueOf(qentry.get("POS")) - shortInterval.getStart();
                                    mods.put(offset, qentry.get("REF") + " (" + alt + ")");
                                    modsRef.put(offset, qentry.get("REF"));
                                    modsAlt.put(offset, alt);

                                    log.debug("mod: found {} at offset {} ({}:{})", qentry.get("REF"), offset, qentry.get("CHROM"), qentry.get("POS"));
                                }
                            }

                            StringBuilder modHapBuilder = new StringBuilder();
                            StringBuilder altHapBuilder = new StringBuilder();
                            for (int i = 0; i < refHap.length(); i++) {
                                if (!mods.containsKey(i)) {
                                    modHapBuilder.append(" ");
                                } else {
                                    modHapBuilder.append(mods.get(i));
                                }

                                if (!modsRef.containsKey(i)) {
                                    altHapBuilder.append(refHap.charAt(i));
                                } else {
                                    altHapBuilder.append(modsAlt.get(i));
                                    i += modsRef.get(i).length() - 1;
                                }
                            }

                            log.debug("refhap: {}", refHap);
                            log.debug("modhap: {}", modHapBuilder.toString());
                            log.debug("althap: {}", altHapBuilder.toString());

                            recombSeqs.get(sampleName).add(altHapBuilder.toString());
                            recombLoci.get(sampleName).put(altHapBuilder.toString(), shortInterval);
                        }
                    }
                }
            }

            log.info("  found {} events in {}", recombSeqs.get(sampleName).size(), sampleName);
        }

        log.info("Looking for contigs that span recombination events...");
        Map<String, Set<String>> recombsFound = new HashMap<String, Set<String>>();

        TableWriter tw = new TableWriter(out);
        TableReader tr = new TableReader(METRICS);
        for (Map<String, String> te : tr) {
            String sampleName = te.get("sampleName");
            String seq = te.get("seq");

            for (String recombFw : recombSeqs.get(sampleName)) {
                String recombRc = SequenceUtils.reverseComplement(recombFw);

                if (seq.contains(recombFw) || seq.contains(recombRc)) {
                    String contigName = te.get("contigName");

                    if (!recombsFound.containsKey(sampleName)) {
                        recombsFound.put(sampleName, new HashSet<String>());
                    }

                    recombsFound.get(sampleName).add(contigName);

                    Map<String, String> entry = new LinkedHashMap<String, String>();
                    entry.put("sampleName", sampleName);
                    entry.put("contigName", contigName);
                    entry.put("locus", recombLoci.get(sampleName).get(recombFw).toString());
                    entry.put("spansRecombinationEvent", "1");

                    tw.addEntry(entry);
                }
            }
        }

        for (String sampleName : recombsFound.keySet()) {
            log.info("  {}: {}", sampleName, recombsFound.get(sampleName).size());
        }
    }
}
