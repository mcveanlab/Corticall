package uk.ac.ox.well.cortexjdk.playground;

import htsjdk.samtools.*;
import uk.ac.ox.well.cortexjdk.commands.Module;
import uk.ac.ox.well.cortexjdk.commands.index.alignedbam.KmerIndex;
import uk.ac.ox.well.cortexjdk.utils.arguments.Argument;
import uk.ac.ox.well.cortexjdk.utils.kmer.CanonicalKmer;
import uk.ac.ox.well.cortexjdk.utils.io.graph.links.CortexJunctionsRecord;
import uk.ac.ox.well.cortexjdk.utils.io.graph.links.CortexLinksIterable;
import uk.ac.ox.well.cortexjdk.utils.io.graph.links.CortexLinksRecord;

import java.io.File;
import java.util.*;

public class FindReadsFromLinks extends Module {
    @Argument(fullName="bam", shortName="b", doc="BAM")
    public File SAM_FILE;

    @Argument(fullName="links", shortName="l", doc="Links")
    public CortexLinksIterable LINKS;

    @Argument(fullName="kmerSize", shortName="k", doc="Kmer size")
    public Integer KMER_SIZE = 47;

    @Override
    public void execute() {
        SamReader sr = SamReaderFactory.make()
                .setOption(SamReaderFactory.Option.INCLUDE_SOURCE_IN_RECORDS, true)
                .setOption(SamReaderFactory.Option.CACHE_FILE_BASED_INDEXES, true)
                .open(SAM_FILE);

        KmerIndex ki = new KmerIndex(SAM_FILE, KMER_SIZE, false);

        for (CortexLinksRecord clr : LINKS) {
            for (CortexJunctionsRecord cjr : clr.getJunctions()) {
                String seq = null; // cjr.getSeq();

                Map<String, SAMRecord> srs = new HashMap<>();

                if (clr.getKmer().equals(new CanonicalKmer("ATAATTATCGCGATAATAATCGTAATAAAAATAATTATCGCGATAAG"))) {
                    log.info("{} {} {} {}", clr.getKmer(), seq.length(), clr.getJunctions().size(), seq);
                    log.info("  {}", cjr);

                    for (int i = 0; i <= seq.length() - KMER_SIZE; i++) {
                        String sk = seq.substring(i, i + KMER_SIZE);

                        List<long[]> ls = ki.find(sk);

                        for (long[] l : ls) {
                            SAMFileSpan sfs = new BAMFileSpan(new Chunk(l[0], l[1]));

                            SAMRecordIterator recs = sr.indexing().iterator(sfs);

                            while (recs.hasNext()) {
                                SAMRecord s = recs.next();
                                srs.put(s.getReadName() + ":" + s.getReadNegativeStrandFlag(), s);
                            }

                            recs.close();

                            //log.info("   {} {}", i, Joiner.on(",").join(readNames));
                        }
                    }
                }

                Map<String, SAMRecord> mrs = new HashMap<>();

                for (String sn : srs.keySet()) {
                    SAMRecord s = srs.get(sn);

                    SAMRecordIterator mates = sr.queryAlignmentStart(s.getMateReferenceName(), s.getMateAlignmentStart());

                    while (mates.hasNext()) {
                        SAMRecord m = mates.next();

                        if (sn.equals(m.getReadName() + ":" + (!m.getReadNegativeStrandFlag()))) {
                            mrs.put(m.getReadName() + ":" + m.getReadNegativeStrandFlag(), m);
                        }
                    }

                    mates.close();
                }

                Set<SAMRecord> ars = new HashSet<>();
                ars.addAll(srs.values());
                ars.addAll(mrs.values());

                int i = 0;
                for (SAMRecord s : ars) {
                    log.info("    {} {}", i, s.getSAMString().trim());
                    i++;
                }
            }
        }
    }
}
