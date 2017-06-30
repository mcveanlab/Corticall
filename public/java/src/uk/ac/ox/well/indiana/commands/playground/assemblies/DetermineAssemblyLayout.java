package uk.ac.ox.well.indiana.commands.playground.assemblies;

import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.reference.FastaSequenceFile;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.containers.ContainerUtils;
import uk.ac.ox.well.indiana.utils.io.gff.GFF3;
import uk.ac.ox.well.indiana.utils.io.gff.GFF3Record;

import java.util.*;

/**
 * Created by kiran on 30/06/2017.
 */
public class DetermineAssemblyLayout extends Module {
    @Argument(fullName="draft", shortName="d", doc="Reference")
    public FastaSequenceFile DRAFT;

    @Argument(fullName="gff", shortName="g", doc="GFF3 file")
    public GFF3 GFF;

    @Argument(fullName="exons", shortName="e", doc="Aligned exons")
    public SamReader EXONS;

    @Override
    public void execute() {
        Map<String, GFF3Record> grrecs = new HashMap<>();
        Map<String, SAMRecord> srrecs = new HashMap<>();

        for (GFF3Record gr : GFF) {
            if (gr.getType().equals("exon")) {
                grrecs.put(gr.getAttribute("ID"), gr);
            }
        }

        Map<String, Map<String, Integer>> freq = new HashMap<>();

        for (SAMRecord sr : EXONS) {
            String id = sr.getReadName();
            srrecs.put(id, sr);

            if (grrecs.containsKey(id) && sr.getCigar().numCigarElements() == 1) {
                //log.info("gff: {}", grrecs.get(id));
                //log.info("sam: {}", srrecs.get(id).getSAMString());
                //log.info("--");

                if (!freq.containsKey(sr.getReferenceName())) {
                    freq.put(sr.getReferenceName(), new HashMap<>());
                }

                if (!freq.get(sr.getReferenceName()).containsKey(grrecs.get(id).getSeqid())) {
                    freq.get(sr.getReferenceName()).put(grrecs.get(id).getSeqid(), 0);
                }

                freq.get(sr.getReferenceName()).put(grrecs.get(id).getSeqid(),
                        freq.get(sr.getReferenceName()).get(grrecs.get(id).getSeqid()) + sr.getReadLength());
            }
        }

        for (String contigName : freq.keySet()) {
            log.info("{} {}", contigName, entriesSortedByValues(freq.get(contigName)));
        }
    }

    static <K,V extends Comparable<? super V>> SortedSet<Map.Entry<K,V>> entriesSortedByValues(Map<K,V> map) {
        SortedSet<Map.Entry<K,V>> sortedEntries = new TreeSet<>(
                new Comparator<Map.Entry<K,V>>() {
                    @Override public int compare(Map.Entry<K,V> e1, Map.Entry<K,V> e2) {
                        int res = e2.getValue().compareTo(e1.getValue());
                        return res != 0 ? res : 1;
                    }
                }
        );

        sortedEntries.addAll(map.entrySet());
        return sortedEntries;
    }
}
