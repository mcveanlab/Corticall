package uk.ac.ox.well.indiana.attic.analyses.nahr;

import net.sf.samtools.*;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.exceptions.IndianaException;
import uk.ac.ox.well.indiana.utils.io.cortex.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.CortexKmer;
import uk.ac.ox.well.indiana.utils.io.cortex.CortexMap;
import uk.ac.ox.well.indiana.utils.io.cortex.CortexRecord;
import uk.ac.ox.well.indiana.utils.io.table.TableWriter;

import java.io.File;
import java.io.PrintStream;
import java.util.*;

public class ComputeParentalContributionTracks extends Module {
    @Argument(fullName="contigs", shortName="c", doc="Aligned contigs (BAM)")
    public SAMFileReader CONTIGS;

    @Argument(fullName="parentGraph", shortName="pg", doc="Parent graph")
    public LinkedHashMap<String, CortexGraph> PARENT_GRAPH;

    @Argument(fullName="maskedKmers", shortName="m", doc="Masked kmers")
    public CortexMap MASKED_KMERS;

    //@Argument(fullName="contigNames", shortName="cn", doc="Contig names")
    //public HashSet<String> CONTIG_NAMES;

    @Output
    public PrintStream out;

    @Output(fullName="bout", shortName="bo", doc="Kmer BAM out")
    public File bout;

    @Output(fullName="bedout", shortName="beo", doc="BED file out")
    public PrintStream beout;

    @Output(fullName="sout", shortName="so", doc="Stats out")
    public PrintStream sout;

    @Override
    public void execute() {
        if (PARENT_GRAPH.size() != 2) {
            throw new IndianaException("Must supply two (and only two) parental graphs");
        }

        log.info("Initializing...");
        int kmerSize = PARENT_GRAPH.values().iterator().next().getKmerSize();

        Map<String, Byte> parentIndices = new TreeMap<String, Byte>();
        byte index = 3;
        for (String parentName : PARENT_GRAPH.keySet()) {
            parentIndices.put(parentName, index);

            index -= 2;
        }
        parentIndices.put("shared", (byte) 2);
        parentIndices.put("ambiguous", (byte) 6);

        log.info("  kmer size: {}", kmerSize);
        log.info("  parent indices: {}", parentIndices);

        log.info("Loading contigs...");

        Set<SAMRecord> contigs = new HashSet<SAMRecord>();
        Set<CortexKmer> contigKmers = new HashSet<CortexKmer>();

        for (SAMRecord contig : CONTIGS) {
            String seq = contig.getReadString();
            for (int i = 0; i <= seq.length() - kmerSize; i++) {
                CortexKmer kmer = new CortexKmer(seq.substring(i, i + kmerSize));

                contigKmers.add(kmer);
            }

            contigs.add(contig);
        }

        log.info("  loaded {} contigs with {} unique kmers", contigs.size(), contigKmers.size());

        log.info("Examining kmer parentage...");

        Map<CortexKmer, String> kmerParentage = new HashMap<CortexKmer, String>();
        for (String parentName : PARENT_GRAPH.keySet()) {
            CortexGraph parentGraph = PARENT_GRAPH.get(parentName);

            int recordsSeen = 0;
            for (CortexRecord cr : parentGraph) {
                CortexKmer ck = cr.getKmer();

                if (MASKED_KMERS.containsKey(ck)) {
                    kmerParentage.put(ck, "ambiguous");
                } else {
                    if (!kmerParentage.containsKey(ck)) {
                        kmerParentage.put(ck, parentName);
                    } else {
                        kmerParentage.put(ck, "shared");
                    }
                }

                if (recordsSeen % (parentGraph.getNumRecords() / 2) == 0) {
                    log.info("  {}: processed {}/{} (~{}%) records", parentName, recordsSeen, parentGraph.getNumRecords(), String.format("%.2f", 100.0f*((float) recordsSeen)/((float) parentGraph.getNumRecords())));
                }
                recordsSeen++;
            }
            //log.info("  {}: processed {}/{} (~{}%) records", parentName, recordsSeen, parentGraph.getNumRecords(), String.format("%.2f", 100.0f*((float) recordsSeen)/((float) parentGraph.getNumRecords())));
        }

        log.info("Constructing inheritance vectors...");

        SAMFileHeader sfh = CONTIGS.getFileHeader();
        SAMReadGroupRecord rgrNone = new SAMReadGroupRecord("none");
        rgrNone.setSample("none");
        SAMReadGroupRecord rgrAmb = new SAMReadGroupRecord("ambiguous");
        rgrAmb.setSample("ambiguous");
        sfh.addReadGroup(rgrNone);
        sfh.addReadGroup(rgrAmb);

        SAMFileWriterFactory sfwf = new SAMFileWriterFactory();
        sfwf.setCreateIndex(true);
        SAMFileWriter sfw = sfwf.makeBAMWriter(sfh, false, bout);

        TableWriter tw = new TableWriter(sout);

        out.printf("%s\t%s\t%s\t%s\t%s\n", "Chromosome", "Start", "End", "Feature", "Parentage");
        int withZero = 0;
        for (SAMRecord contig : contigs) {
            //if (CONTIG_NAMES.contains(contig.getReadName())) {
                String seq = contig.getReadString();

                beout.println(contig.getReferenceName() + "\t" + contig.getAlignmentStart() + "\t" + contig.getAlignmentEnd());

                int p1 = 0, p2 = 0, shared = 0, ambiguous = 0, none = 0;
                for (AlignmentBlock ab : contig.getAlignmentBlocks()) {
                    //log.info("{}: {} {} {}", contig.getReadName(), ab.getReadStart(), ab.getReferenceStart(), ab.getLength());

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

                            if (kmerParentage.containsKey(kmer)) {
                                String parentName = kmerParentage.get(kmer);
                                byte parentIndex = parentIndices.get(parentName);

                                if (parentIndex == 3) { p1++; }
                                if (parentIndex == 1) { p2++; }
                                if (parentIndex == 2) { shared++; }
                                if (parentIndex == -1) { ambiguous++; }

                                String parentReadGroupId = null;
                                String sampleReadGroupId = null;
                                for (SAMReadGroupRecord rgr : sfh.getReadGroups()) {
                                    if (rgr.getSample().equals(parentName)) {
                                        parentReadGroupId = rgr.getReadGroupId();
                                    }

                                    if (rgr.getSample().equals(contig.getReadGroup().getSample())) {
                                        sampleReadGroupId = rgr.getReadGroupId();
                                    }
                                }

                                out.printf("%s\t%d\t%d\t%s\t%d\n", contig.getReferenceName(), contig.getAlignmentStart() + i, contig.getAlignmentStart() + i + 1, "parentage", parentIndex);

                                skmer.setAttribute("RG", parentReadGroupId != null ? parentReadGroupId : sampleReadGroupId);
                                skmer.setMappingQuality(99);
                            } else {
                                none++;
                            }

                            sfw.addAlignment(skmer);
                        }
                    }
                }

                //log.info("{}:{} length={} p1={} p2={} shared={} none={}", contig.getReferenceName(), contig.getAlignmentStart(), contig.getReadLength(), p1, p2, shared, none);
                Map<String, String> stats = new LinkedHashMap<String, String>();
                stats.put("contigName", contig.getReadName());
                stats.put("sampleName", contig.getReadGroup().getSample());
                stats.put("p1", String.valueOf(p1));
                stats.put("p2", String.valueOf(p2));
                stats.put("shared", String.valueOf(shared));
                stats.put("ambiguous", String.valueOf(ambiguous));
                stats.put("none", String.valueOf(none));
                tw.addEntry(stats);

                if (p1 == 0 || p2 == 0) {
                    withZero++;
                }
            //}
        }

        sfw.close();
    }
}
