package uk.ac.ox.well.indiana.attic.analyses.nahr;

import net.sf.picard.util.Interval;
import net.sf.picard.util.IntervalTreeMap;
import net.sf.samtools.CigarElement;
import net.sf.samtools.CigarOperator;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.exceptions.IndianaException;
import uk.ac.ox.well.indiana.utils.io.table.TableReader;
import uk.ac.ox.well.indiana.utils.io.table.TableWriter;
import uk.ac.ox.well.indiana.utils.io.utils.LineReader;

import java.io.File;
import java.io.PrintStream;
import java.util.*;

public class FindContigsThatSpanRecombinationEvents extends Module {
    @Argument(fullName="bam", shortName="b", doc="BAM file")
    public ArrayList<SAMFileReader> BAMS;

    @Argument(fullName="recombLoci", shortName="rl", doc="Recomb loci table")
    public File RECOMB_LOCI;

    @Argument(fullName="vcf", shortName="vcf", doc="VCF file")
    public File VCF;

    @Argument(fullName="parent1", shortName="p1", doc="Parent 1")
    public String PARENT1;

    @Argument(fullName="parent2", shortName="p2", doc="Parent 2")
    public String PARENT2;

    @Output
    public PrintStream out;

    private IntervalTreeMap<Map<String, String>> loadVCF(File vcf) {
        IntervalTreeMap<Map<String, String>> vcfRecords = new IntervalTreeMap<Map<String, String>>();

        LineReader lr = new LineReader(vcf);

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

                    Map<String, String> entry = new LinkedHashMap<String, String>();
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

    private Set<String> getSampleNames(File vcf) {
        Set<String> sampleNames = new LinkedHashSet<String>();

        LineReader lr = new LineReader(vcf);
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
        Map<String, Set<Interval>> loci = new HashMap<String, Set<Interval>>();

        IntervalTreeMap<Map<String, String>> vcf = loadVCF(VCF);
        Set<String> sampleNames = getSampleNames(VCF);

        TableReader tr = new TableReader(RECOMB_LOCI);

        for (Map<String, String> te : tr) {
            String[] sampleFields = te.get("sample").split("/");
            String sampleName = sampleFields[1];

            if (sampleNames.contains(sampleName)) {
                Interval leftInterval = null;
                Interval centerInterval = new Interval(te.get("chrom"), Integer.valueOf(te.get("pos_max")), Integer.valueOf(te.get("pos_max")));
                Interval rightInterval = null;

                Map<String, String> leftVariant = null;
                Map<String, String> centerVariant = null;
                Map<String, String> rightVariant = null;

                String leftHapId = null;
                String centerHapId = null;
                String rightHapId = null;

                int i;

                i = 1;
                do {
                    leftInterval = new Interval(te.get("chrom"), Integer.valueOf(te.get("pos_max")) - i, Integer.valueOf(te.get("pos_max")) - i);

                    if (vcf.containsKey(leftInterval)) { leftVariant = vcf.get(leftInterval); }
                    i++;
                } while (leftVariant == null);

                centerVariant = vcf.get(centerInterval);

                i = 1;
                do {
                    rightInterval = new Interval(te.get("chrom"), Integer.valueOf(te.get("pos_max")) + i, Integer.valueOf(te.get("pos_max")) + i);

                    if (vcf.containsKey(rightInterval)) { rightVariant = vcf.get(rightInterval); }
                    i++;
                } while (rightVariant == null);

                leftHapId = leftVariant.get(sampleName).equals(leftVariant.get(PARENT1)) ? PARENT1 : PARENT2;
                centerHapId = centerVariant.get(sampleName).equals(centerVariant.get(PARENT1)) ? PARENT1 : PARENT2;
                rightHapId = rightVariant.get(sampleName).equals(rightVariant.get(PARENT1)) ? PARENT1 : PARENT2;

                //log.info("{} {}", sampleName, centerInterval);
                //log.info("  left: {}:{} p1={} p2={} c={} {}", leftVariant.get("CHROM"), leftVariant.get("POS"), leftVariant.get(PARENT1), leftVariant.get(PARENT2), leftVariant.get(sampleName), leftVariant.get(sampleName).equals(leftVariant.get(PARENT1)) ? PARENT1 : PARENT2);
                //log.info("center: {}:{} p1={} p2={} c={} {}", centerVariant.get("CHROM"), centerVariant.get("POS"), centerVariant.get(PARENT1), centerVariant.get(PARENT2), centerVariant.get(sampleName), centerVariant.get(sampleName).equals(centerVariant.get(PARENT1)) ? PARENT1 : PARENT2);
                //log.info(" right: {}:{} p1={} p2={} c={} {}", rightVariant.get("CHROM"), rightVariant.get("POS"), rightVariant.get(PARENT1), rightVariant.get(PARENT2), rightVariant.get(sampleName), rightVariant.get(sampleName).equals(rightVariant.get(PARENT1)) ? PARENT1 : PARENT2);
                //log.info("");

                if (!loci.containsKey(sampleName)) {
                    loci.put(sampleName, new HashSet<Interval>());
                }

                if (!leftHapId.equals(centerHapId)) {
                    loci.get(sampleName).add(new Interval(te.get("chrom"), leftInterval.getStart(), centerInterval.getEnd()));
                } else if (!centerHapId.equals(rightHapId)) {
                    //loci.get(sampleName).add(new Interval(te.get("chrom"), centerInterval.getStart(), rightInterval.getEnd()));

                    throw new IndianaException("Weird haplotype structure");
                }
            }
        }

        TableWriter tw = new TableWriter(out);

        for (SAMFileReader BAM : BAMS) {
            String sampleName = BAM.getFileHeader().getReadGroups().iterator().next().getSample();

            Map<String, Integer> numAlignments = new HashMap<String, Integer>();
            List<SAMRecord> contigs = new ArrayList<SAMRecord>();
            for (SAMRecord contig : BAM) {
                if (!numAlignments.containsKey(contig.getReadName())) {
                    numAlignments.put(contig.getReadName(), 1);
                } else {
                    numAlignments.put(contig.getReadName(), numAlignments.get(contig.getReadName()) + 1);
                }

                contigs.add(contig);
            }

            for (SAMRecord contig : contigs) {
                Interval contigLocus = new Interval(contig.getReferenceName(), contig.getAlignmentStart(), contig.getAlignmentEnd());

                if (loci.containsKey(sampleName)) {
                    for (Interval knownLocus : loci.get(sampleName)) {
                        if (knownLocus.intersects(contigLocus)) {
                            boolean isClipped = false;
                            for (CigarElement ce : contig.getCigar().getCigarElements()) {
                                if (ce.getOperator().equals(CigarOperator.S) || ce.getOperator().equals(CigarOperator.H)) {
                                    isClipped = true;
                                    break;
                                }
                            }

                            boolean perfectAlignment = false;
                            if (contig.getCigar().getCigarElements().size() == 1 && contig.getCigar().getCigarElement(0).getOperator().equals(CigarOperator.M)) {
                                String md = contig.getStringAttribute("MD");
                                try {
                                    int length = Integer.parseInt(md);
                                    if (length == contig.getReadLength()) {
                                        perfectAlignment = true;
                                    }
                                } catch (NumberFormatException e) {}
                            }

                            Map<String, String> te = new LinkedHashMap<String, String>();
                            te.put("sampleName", sampleName);
                            te.put("contigName", contig.getReadName());
                            te.put("contigLocus", contigLocus.toString());
                            te.put("knownLocus", knownLocus.toString());
                            te.put("contigLength", String.valueOf(contig.getReadLength()));
                            te.put("knownLocusLength", String.valueOf(knownLocus.length()));
                            te.put("numAlignments", String.valueOf(numAlignments.get(contig.getReadName())));
                            te.put("isClipped", isClipped ? "1" : "0");
                            te.put("pctKnownRegionSpanned", String.valueOf(((float) knownLocus.getIntersectionLength(contigLocus))/((float) knownLocus.length())));
                            te.put("perfectAlignment", perfectAlignment ? "1" : "0");
                            te.put("cigar", contig.getCigarString());
                            te.put("md", contig.getStringAttribute("MD"));

                            tw.addEntry(te);
                        }
                    }
                }
            }
        }
    }
}
