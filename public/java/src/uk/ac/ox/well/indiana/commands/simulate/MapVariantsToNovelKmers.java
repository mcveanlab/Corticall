package uk.ac.ox.well.indiana.commands.simulate;

import htsjdk.samtools.reference.FastaSequenceFile;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexKmer;
import uk.ac.ox.well.indiana.utils.io.table.TableReader;
import uk.ac.ox.well.indiana.utils.io.table.TableWriter;
import uk.ac.ox.well.indiana.utils.io.utils.LineReader;

import java.io.File;
import java.io.PrintStream;
import java.util.*;

public class MapVariantsToNovelKmers extends Module {
    @Argument(fullName="reference", shortName="r", doc="Reference")
    public IndexedFastaSequenceFile REFERENCE;

    @Argument(fullName="novelKmers", shortName="n", doc="List of novel kmers")
    public File NOVEL_KMERS;

    @Argument(fullName="bed", shortName="b", doc="Bed file")
    public File BED;

    @Output
    public PrintStream out;

    private Map<String, String> parseInfoField(String info) {
        Map<String, String> infoMap = new HashMap<String, String>();
        String[] entries = info.split(";");

        for (String entry : entries) {
            String[] fields = entry.split("=");

            infoMap.put(fields[0], fields[1]);
        }

        return infoMap;
    }

    private Set<CortexKmer> loadNovelKmers() {
        LineReader lr = new LineReader(NOVEL_KMERS);

        Set<CortexKmer> kmers = new HashSet<CortexKmer>();
        while (lr.hasNext()) {
            String l = lr.getNextRecord();

            CortexKmer ck = new CortexKmer(l);
            kmers.add(ck);
        }

        return kmers;
    }

    @Override
    public void execute() {
        log.info("Loading novel kmers...");
        Set<CortexKmer> novelKmers = loadNovelKmers();
        log.info("  {} novel kmers seen", novelKmers.size());

        log.info("Assigning novel kmers to variants...");
        int kmerSize = novelKmers.iterator().next().getKmerAsBytes().length;

        TableWriter tw = new TableWriter(out);

        TableReader tr = new TableReader(BED, new String[] {"chrom", "start", "stop", "info"});
        for (Map<String, String> te : tr) {
            Map<String, String> infoMap = parseInfoField(te.get("info"));

            String chrom = te.get("chrom");
            int start = Integer.valueOf(te.get("start")) + 1 - (kmerSize-1);
            int stop  = Integer.valueOf(te.get("stop")) + 1 + (kmerSize-1);

            if (start < 0) { start = 0; }
            if (stop > REFERENCE.getSequence(chrom).length()) {
                stop = REFERENCE.getSequence(chrom).length();
            }

            String seq = new String(REFERENCE.getSubsequenceAt(chrom, start, stop).getBases());

            for (int i = 0; i <= seq.length() - kmerSize; i++) {
                CortexKmer ck = new CortexKmer(seq.substring(i, i + kmerSize));

                if (novelKmers.contains(ck)) {
                    String vclass = infoMap.get("denovo");
                    if (!infoMap.get("nahr").equals("unknown")) { vclass = "NAHR"; }
                    if (vclass.equals("unknown")) { vclass = "inherited_" + infoMap.get("type"); }

                    Map<String, String> twe = new LinkedHashMap<String, String>();
                    twe.put("kmer", ck.getKmerAsString());
                    twe.put("variantId", infoMap.get("id"));
                    twe.put("vclass", vclass);
                    twe.put("vchr", chrom);
                    twe.put("vstart", String.valueOf(start));
                    twe.put("vstop", String.valueOf(stop));
                    tw.addEntry(twe);

                    novelKmers.remove(ck);
                }
            }
        }

        for (CortexKmer ck : novelKmers) {
            Map<String, String> twe = new LinkedHashMap<String, String>();
            twe.put("kmer", ck.getKmerAsString());
            twe.put("variantId", "NA");
            twe.put("vclass", "NA");
            twe.put("vchr", "NA");
            twe.put("vstart", "NA");
            twe.put("vstop", "NA");
            tw.addEntry(twe);
        }

        log.info("  {} novel kmers unassigned to variants", novelKmers.size());
    }
}
