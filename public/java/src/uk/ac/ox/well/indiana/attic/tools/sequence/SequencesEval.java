package uk.ac.ox.well.indiana.attic.tools.sequence;

import htsjdk.samtools.reference.FastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.alignment.sw.NeedlemanWunsch;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.table.TableReader;
import uk.ac.ox.well.indiana.utils.performance.PerformanceUtils;
import uk.ac.ox.well.indiana.utils.sequence.SequenceUtils;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;

public class SequencesEval extends Module {
    @Argument(fullName="fasta", shortName="f", doc="Fasta file with sequences to compare")
    public FastaSequenceFile FASTA;

    @Argument(fullName="metaData", shortName="md", doc="Metadata on the sequences")
    public File METADATA;

    @Output
    public PrintStream out;

    private class SequenceRecord {
        public String strainId;
        public String accession;
        public String sequence;
    }

    private HashMap<Integer, ArrayList<SequenceRecord>> loadData(FastaSequenceFile fasta, File metaDataFile) {
        HashMap<String, ReferenceSequence> sequences = new HashMap<String, ReferenceSequence>();

        ReferenceSequence seq;
        while ((seq = fasta.nextSequence()) != null) {
            String[] name = seq.getName().split("\\|");

            sequences.put(name[3], seq);
        }

        TableReader metadata = new TableReader(metaDataFile);

        HashMap<Integer, ArrayList<SequenceRecord>> results = new HashMap<Integer, ArrayList<SequenceRecord>>();

        for (Map<String, String> entry : metadata) {
            Integer segment = Integer.valueOf(entry.get("segment_number"));

            SequenceRecord record = new SequenceRecord();
            record.strainId = entry.get("strain_id");
            record.accession = entry.get("gb_accession");
            record.sequence = new String(sequences.get(record.accession).getBases());

            if (!results.containsKey(segment)) {
                results.put(segment, new ArrayList<SequenceRecord>());
            }

            results.get(segment).add(record);
        }

        return results;
    }

    @Override
    public void execute() {
        log.info(PerformanceUtils.getMemoryUsageStats());

        HashMap<Integer, ArrayList<SequenceRecord>> data = loadData(FASTA, METADATA);

        for (int segmentNumber : data.keySet()) {
            log.info("Sequence: {}", segmentNumber);

            ArrayList<SequenceRecord> records = data.get(segmentNumber);

            HashMap<String, HashMap<String, Integer>> comparisonMatrix = new HashMap<String, HashMap<String, Integer>>();

            for (int i = 0; i < records.size(); i++) {
                //if (i % (records.size() / 10) == 0) {
                    log.info("Sequence {}/{}", i, records.size());
                //}

                for (int j = i + 1; j < records.size(); j++) {
                    log.info("\tSequence {}/{}", j, records.size());

                    SequenceRecord si = records.get(i);
                    SequenceRecord sj = records.get(j);

                    NeedlemanWunsch nw = new NeedlemanWunsch(si.sequence, sj.sequence);

                    int score = nw.getAlignmentScore();
                    String[] alignment = nw.getAlignment();

//                    log.info("{} {}:", si.strainId, sj.strainId);
//                    log.info("score: {}", score);
//                    log.info("si: {}", alignment[0]);
//                    log.info("sj: {}", alignment[1]);
//                    log.info("L: {}, S: {}", alignment[0].length(), SequenceUtils.numSegregatingSites(alignment[0], alignment[1]));

                    if (!comparisonMatrix.containsKey(si.strainId)) {
                        comparisonMatrix.put(si.strainId, new HashMap<String, Integer>());
                    }

                    comparisonMatrix.get(si.strainId).put(sj.strainId, SequenceUtils.numSegregatingSites(alignment[0], alignment[1]));
                }
            }

            for (String strainId1 : comparisonMatrix.keySet()) {
                for (String strainId2 : comparisonMatrix.get(strainId1).keySet()) {
                    out.printf("%d\t%s\t%s\t%d\n", segmentNumber, strainId1, strainId2, comparisonMatrix.get(strainId1).get(strainId2));
                }
            }
        }

//        log.info(PerformanceUtils.getMemoryUsageStats());
    }
}
