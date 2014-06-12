package uk.ac.ox.well.indiana.attic.analyses.nahr;

import net.sf.samtools.AlignmentBlock;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.exceptions.IndianaException;
import uk.ac.ox.well.indiana.utils.io.cortex.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.CortexKmer;
import uk.ac.ox.well.indiana.utils.io.cortex.CortexRecord;

import java.io.PrintStream;
import java.util.*;

public class ComputeParentalContributionTracks extends Module {
    @Argument(fullName="contigs", shortName="c", doc="Aligned contigs (BAM)")
    public SAMFileReader CONTIGS;

    @Argument(fullName="parentGraph", shortName="pg", doc="Parent graph")
    public LinkedHashMap<String, CortexGraph> PARENT_GRAPH;

    @Output
    public PrintStream out;

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

        log.info("  kmer size: {}", kmerSize);
        log.info("  parent indices: {}", parentIndices);

        log.info("Loading contigs...");

        Set<SAMRecord> contigs = new HashSet<SAMRecord>();
        Set<CortexKmer> contigKmers = new HashSet<CortexKmer>();

        for (SAMRecord contig : CONTIGS) {
            log.info("{}", contig.getSAMString());

            for (AlignmentBlock ab : contig.getAlignmentBlocks()) {
                log.info("  {} {} {}", ab.getReadStart(), ab.getReferenceStart(), ab.getLength());
            }

            log.info("");

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

                if (!kmerParentage.containsKey(ck)) {
                    kmerParentage.put(ck, parentName);
                } else {
                    kmerParentage.put(ck, "shared");
                }

                if (recordsSeen % (parentGraph.getNumRecords() / 2) == 0) {
                    log.info("  {}: processed {}/{} (~{}%) records", parentName, recordsSeen, parentGraph.getNumRecords(), 100*recordsSeen/parentGraph.getNumRecords());
                }
                recordsSeen++;
            }
            log.info("  {}: processed {}/{} (~{}%) records", parentName, recordsSeen, parentGraph.getNumRecords(), 100*recordsSeen/parentGraph.getNumRecords());
        }

        log.info("Constructing inheritance vectors...");


        out.printf("%s\t%s\t%s\t%s\t%s\n", "Chromosome", "Start", "End", "Feature", "Parentage");
        for (SAMRecord contig : contigs) {
            String seq = contig.getReadString();

            Map<String, Byte[]> wig = new HashMap<String, Byte[]>();
            for (String parentName : parentIndices.keySet()) {
                wig.put(parentName, new Byte[seq.length()]);
            }

            for (int i = 0; i <= seq.length() - kmerSize; i++) {
                CortexKmer kmer = new CortexKmer(seq.substring(i, i + kmerSize));

                if (kmerParentage.containsKey(kmer)) {
                    String parentName = kmerParentage.get(kmer);
                    byte parentIndex = parentIndices.get(parentName);

                    int j = i;
                    //for (int j = i; j < i + kmerSize; j++) {
                        wig.get(parentName)[j] = parentIndex;
                    //}
                }
            }

            log.info("{}: {}", contig.getReadName(), contig.getReadString());

            for (String parentName : wig.keySet()) {
                StringBuilder sb = new StringBuilder();

                for (int i = 0; i < wig.get(parentName).length; i++) {
                    Byte wigValue = wig.get(parentName)[i];

                    if (wigValue == null) {
                        sb.append(" ");
                    } else {
                        sb.append(wigValue);
                    }
                }

                log.info("{}: {}", contig.getReadName(), sb.toString());
            }

            log.info("");

            Byte[] finalWig = new Byte[seq.length()];
            for (int i = 0; i < finalWig.length; i++) {
                if (wig.get("shared")[i] != null) {
                    finalWig[i] = parentIndices.get("shared");
                } else {
                    for (String parentName : wig.keySet()) {
                        if (!parentName.equals("shared") && wig.get(parentName)[i] != null) {
                            finalWig[i] = parentIndices.get(parentName);
                        }
                    }
                }
            }

            StringBuilder sb = new StringBuilder();
            for (int i = 0; i < finalWig.length; i++) {
                if (finalWig[i] == null) {
                    sb.append(" ");
                } else {
                    sb.append(finalWig[i]);

                    out.printf("%s\t%d\t%d\t%s\t%d\n", contig.getReferenceName(), contig.getAlignmentStart() + i, contig.getAlignmentStart() + i + 1, "parentage", finalWig[i]);
                }
            }

            log.info("{}: {}", contig.getReadName(), sb.toString());

            log.info("");


        }
    }
}
