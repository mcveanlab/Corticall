package uk.ac.ox.well.indiana.attic.analyses.nahr;

import htsjdk.samtools.*;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.cortex.CortexKmer;
import uk.ac.ox.well.indiana.utils.io.table.TableReader;
import uk.ac.ox.well.indiana.utils.io.table.TableWriter;

import java.io.File;
import java.io.PrintStream;
import java.util.*;

public class ComputeParentalContributionTracks extends Module {
    @Argument(fullName="contigs", shortName="c", doc="Aligned contigs (BAM)")
    public SAMFileReader CONTIGS;

    @Argument(fullName="ann", shortName="ann", doc="Contig annotations")
    public File ANN;

    @Argument(fullName="kmerCodeName", shortName="kcn", doc="Kmer code names", required=false)
    public HashMap<String, String> KMER_CODE_NAMES;

    @Argument(fullName="contigNames", shortName="cn", doc="Contigs to select", required=false)
    public HashSet<String> CONTIG_NAMES;

    @Output
    public PrintStream out;

    @Output(fullName="bout", shortName="bo", doc="Kmer BAM out")
    public File bout;

    @Output(fullName="bedout", shortName="beo", doc="BED file out")
    public PrintStream beout;

    @Output(fullName="sout", shortName="so", doc="Stats out")
    public PrintStream sout;

    private class IGVEntry implements Comparable<IGVEntry> {
        public String chromosome;
        public int start;
        public String parentageName;
        public int parentage;

        @Override
        public int compareTo(IGVEntry o) {
            if (!chromosome.equals(o.chromosome)) { return chromosome.compareTo(o.chromosome); }

            return Integer.valueOf(start).compareTo(o.start);
        }

        @Override
        public boolean equals(Object o) {
            if (o instanceof IGVEntry) {
                IGVEntry oi = (IGVEntry) o;
                return (chromosome.equals(oi.chromosome) && start == oi.start && parentage == oi.parentage);
            }

            return false;
        }

        @Override
        public int hashCode() {
            return chromosome.hashCode() + Integer.valueOf(start).hashCode() - Integer.valueOf(parentage).hashCode();
        }
    }

    @Override
    public void execute() {
        log.info("Initializing kmer code map...");
        Map<Character, Integer> kmerCodeIndices = new HashMap<Character, Integer>();
        kmerCodeIndices.put('0', 1);
        kmerCodeIndices.put('A', 3);
        kmerCodeIndices.put('B', 4);
        kmerCodeIndices.put('C', 5);
        kmerCodeIndices.put('_', 6);
        kmerCodeIndices.put('.', 7);
        kmerCodeIndices.put('1', 9);

        Map<Character, String> kmerCodeNames = new LinkedHashMap<Character, String>();
        kmerCodeNames.put('0', "ref0");
        kmerCodeNames.put('A', "repetitive");
        kmerCodeNames.put('B', "both");
        kmerCodeNames.put('C', "lowcoverage");
        kmerCodeNames.put('_', "lowconfidence");
        kmerCodeNames.put('.', "novel");
        kmerCodeNames.put('1', "ref1");

        if (KMER_CODE_NAMES != null) {
            for (Character c : kmerCodeNames.keySet()) {
                String cStr = String.valueOf(c);
                if (KMER_CODE_NAMES.containsKey(cStr)) {
                    kmerCodeNames.put(c, KMER_CODE_NAMES.get(cStr));
                }
            }
        }

        for (Character c : kmerCodeNames.keySet()) {
            log.info("  {} {}: {}", c, kmerCodeIndices.get(c), kmerCodeNames.get(c));
        }

        log.info("Loading annotated contigs...");
        Map<String, Map<String, String>> annotatedContigs = new HashMap<String, Map<String, String>>();
        int kmerSize = 0;

        TableReader tr = new TableReader(ANN);
        for (Map<String, String> te : tr) {
            String contigName = te.get("contigName");

            if (kmerSize == 0) {
                kmerSize = te.get("seq").length() - te.get("kmerOrigin").length() + 1;
            }

            annotatedContigs.put(contigName, te);
        }

        log.info("    contigs: {}", annotatedContigs.size());
        log.info("  kmer size: {}", kmerSize);

        log.info("Computing kmer inheritance information...");

        SAMFileHeader sfh = CONTIGS.getFileHeader();
        for (Character c : kmerCodeNames.keySet()) {
            SAMReadGroupRecord rgr = new SAMReadGroupRecord(kmerCodeNames.get(c));
            rgr.setSample(kmerCodeNames.get(c));
            sfh.addReadGroup(rgr);
        }

        SAMFileWriterFactory sfwf = new SAMFileWriterFactory();
        sfwf.setCreateIndex(true);
        SAMFileWriter sfw = sfwf.makeBAMWriter(sfh, false, bout);

        TableWriter tw = new TableWriter(sout);

        Set<IGVEntry> igvEntries = new TreeSet<IGVEntry>();
        int numContigs = 0;
        for (SAMRecord contig : CONTIGS) {
            if (CONTIG_NAMES == null || CONTIG_NAMES.isEmpty() || CONTIG_NAMES.contains(contig.getReadName())) {
                Map<String, String> te = annotatedContigs.get(contig.getReadName());

                if (annotatedContigs.containsKey(contig.getReadName())) {
                    String seq = contig.getReadString();

                    //log.debug("  te: {}", te);

                    String annSeq = te.get("seq");
                    String kmerOrigin = te.get("kmerOrigin");

                    Map<CortexKmer, Character> kmerCodes = new HashMap<CortexKmer, Character>();
                    for (int i = 0; i < kmerOrigin.length(); i++) {
                        CortexKmer kmer = new CortexKmer(annSeq.substring(i, i + kmerSize));
                        Character code = kmerOrigin.charAt(i);

                        kmerCodes.put(kmer, code);
                    }

                    Map<Character, Integer> kmerStats = new HashMap<Character, Integer>();
                    for (Character c : kmerCodeNames.keySet()) {
                        kmerStats.put(c, 0);
                    }

                    boolean changed = false;

                    // We want to be able to examine soft-clipped regions as well.
                    List<CigarElement> ces = new ArrayList<CigarElement>();
                    for (CigarElement ce : contig.getCigar().getCigarElements()) {
                        if (ce.getOperator().equals(CigarOperator.S)) {
                            ces.add(new CigarElement(ce.getLength(), CigarOperator.M));
                            changed = true;
                        } else {
                            ces.add(ce);
                        }
                    }

                    if (changed) {
                        CigarElement firstCe = contig.getCigar().getCigarElements().get(0);

                        if (firstCe.getOperator().equals(CigarOperator.S)) {
                            contig.setAlignmentStart(contig.getAlignmentStart() - firstCe.getLength());
                        }

                        contig.setCigar(new Cigar(ces));
                    }

                    for (AlignmentBlock ab : contig.getAlignmentBlocks()) {
                        for (int i = ab.getReadStart() - 1; i < ab.getReadStart() + ab.getLength(); i++) {
                            if (i + kmerSize < seq.length()) {
                                CortexKmer kmer = new CortexKmer(seq.substring(i, i + kmerSize));

                                SAMRecord skmer = new SAMRecord(CONTIGS.getFileHeader());
                                skmer.setReadBases(seq.substring(i, i + kmerSize).getBytes());

                                List<CigarElement> cigarElements = new ArrayList<CigarElement>();
                                cigarElements.add(new CigarElement(kmerSize, CigarOperator.M));
                                Cigar cigar = new Cigar(cigarElements);

                                skmer.setReadName(contig.getReadName() + "." + kmer.getKmerAsString());
                                skmer.setReferenceName(contig.getReferenceName());
                                skmer.setCigar(cigar);
                                skmer.setReadPairedFlag(false);
                                skmer.setDuplicateReadFlag(false);
                                skmer.setMateNegativeStrandFlag(false);
                                skmer.setAlignmentStart(ab.getReferenceStart() - ab.getReadStart() + 1 + i);
                                skmer.setAttribute("RG", "none");
                                skmer.setMappingQuality(0);

                                Character c = kmerCodes.get(kmer);
                                String codeName = kmerCodeNames.get(c);

                                String parentReadGroupId = null;
                                String sampleReadGroupId = null;
                                for (SAMReadGroupRecord rgr : sfh.getReadGroups()) {
                                    if (rgr.getSample().equals(codeName)) {
                                        parentReadGroupId = rgr.getReadGroupId();
                                    }

                                    if (rgr.getSample().equals(contig.getReadGroup().getSample())) {
                                        sampleReadGroupId = rgr.getReadGroupId();
                                    }
                                }

                                skmer.setAttribute("RG", parentReadGroupId != null ? parentReadGroupId : sampleReadGroupId);
                                skmer.setMappingQuality(99);

                                sfw.addAlignment(skmer);

                                kmerStats.put(c, kmerStats.get(c) + 1);

                                IGVEntry igvEntry = new IGVEntry();
                                igvEntry.chromosome = contig.getReferenceName();
                                igvEntry.start = ab.getReferenceStart() - ab.getReadStart() + i;
                                igvEntry.parentageName = kmerCodeNames.get(c);
                                igvEntry.parentage = kmerCodeIndices.get(c);
                                igvEntries.add(igvEntry);
                            }
                        }
                    }

                    if (!contig.isSecondaryOrSupplementary()) {
                        beout.println(contig.getReferenceName() + "\t" + contig.getAlignmentStart() + "\t" + contig.getAlignmentEnd() + "\t" + contig.getReadName() + "." + contig.getReadGroup().getSample());

                        if (numContigs % (annotatedContigs.size() / 10) == 0) {
                            log.info("  processed {}/{} contigs", numContigs, annotatedContigs.size());
                        }
                        numContigs++;
                    }

                    Map<String, String> stats = new LinkedHashMap<String, String>();
                    stats.put("contigName", contig.getReadName());
                    stats.put("sampleName", contig.getReadGroup().getSample());
                    for (Character c : kmerCodeNames.keySet()) {
                        stats.put(kmerCodeNames.get(c), String.valueOf(kmerStats.get(c)));
                    }
                    tw.addEntry(stats);
                }
            }
        }

        log.info("Writing kmer inheritance information...");
        out.printf("%s\t%s\t%s\t%s\t%s\n", "Chromosome", "Start", "End", "Feature", "Parentage");
        for (IGVEntry igvEntry : igvEntries) {
            out.printf("%s\t%d\t%d\t%s\t%d\n", igvEntry.chromosome, igvEntry.start, igvEntry.start + 1, igvEntry.parentageName, igvEntry.parentage);
        }

        sfw.close();
    }
}
