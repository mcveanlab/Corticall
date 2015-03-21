package uk.ac.ox.well.indiana.commands.simulate;

import com.google.common.base.Joiner;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.reference.FastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.tribble.readers.TabixIteratorLineReader;
import htsjdk.tribble.readers.TabixReader;
import htsjdk.tribble.util.TabixUtils;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.exceptions.IndianaException;
import uk.ac.ox.well.indiana.utils.io.table.TableReader;
import uk.ac.ox.well.indiana.utils.io.table.TableWriter;
import uk.ac.ox.well.indiana.utils.io.utils.LineReader;
import uk.ac.ox.well.indiana.utils.statistics.misc.StatisticsOnStream;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.*;

public class LiftoverFromRefToChild extends Module {
    @Argument(fullName="reference", shortName="r", doc="Reference")
    public FastaSequenceFile REF;

    @Argument(fullName="vcf", shortName="v", doc="VCF")
    public VCFFileReader VCF;

    @Argument(fullName="tableToLift", shortName="t", doc="Table to lift")
    public File TABLE;

    @Argument(fullName="seed", shortName="g", doc="Seed for RNG", required=false)
    public Long SEED;

    @Output
    public PrintStream out;

    private Random rng;

    private Map<String, Map<Integer, Set<VariantContext>>> loadVariants() {
        Map<String, Map<Integer, Set<VariantContext>>> variants = new HashMap<String, Map<Integer, Set<VariantContext>>>();
        for (VariantContext vc : VCF) {
            String chr = vc.getChr();
            int start = vc.getStart();

            if (!variants.containsKey(chr)) {
                variants.put(chr, new HashMap<Integer, Set<VariantContext>>());
            }

            if (!variants.get(chr).containsKey(start)) {
                variants.get(chr).put(start - 1, new HashSet<VariantContext>());
            }

            variants.get(chr).get(start - 1).add(vc);
        }

        return variants;
    }

    private class Entry {
        public int numReadsWithErrors;
        public int numReads;
        public int numFragmentsWithErrors;
        public int numFragments;
    }

    private Map<String, Map<Integer, Entry>> loadTable() {
        Map<String, Map<Integer, Entry>> table = new HashMap<String, Map<Integer, Entry>>();

        LineReader lr = new LineReader(TABLE);
        lr.getNextRecord();

        while (lr.hasNext()) {
            String line = lr.getNextRecord();
            String[] fields = line.split("\\s+");

            String chr = fields[0];
            int pos = Integer.valueOf(fields[1]);

            Entry e = new Entry();
            e.numReadsWithErrors = Integer.valueOf(fields[2]);
            e.numReads = Integer.valueOf(fields[3]);
            e.numFragmentsWithErrors = Integer.valueOf(fields[4]);
            e.numFragments = Integer.valueOf(fields[5]);

            if (!table.containsKey(chr)) {
                table.put(chr, new HashMap<Integer, Entry>());
            }

            if (!table.get(chr).containsKey(pos)) {
                table.get(chr).put(pos, e);
            }
        }

        return table;
    }

    private List<Entry> computeEntries(Map<Integer, Entry> t, int pos, int length) {
        StatisticsOnStream numReadsWithErrors = new StatisticsOnStream();
        StatisticsOnStream numReads = new StatisticsOnStream();
        StatisticsOnStream numFragmentsWithErrors = new StatisticsOnStream();
        StatisticsOnStream numFragments = new StatisticsOnStream();

        for (int i = 0; i < length; i++) {
            if (t.containsKey(pos - i)) {
                Entry e = t.get(pos - i);

                numReadsWithErrors.push(e.numReadsWithErrors);
                numReads.push(e.numReads);
                numFragmentsWithErrors.push(e.numFragmentsWithErrors);
                numFragments.push(e.numFragments);
            }
        }

        List<Entry> entries = new ArrayList<Entry>();
        for (int i = 0; i < length; i++) {
            Entry e = new Entry();

            e.numReadsWithErrors = (int) Math.ceil(rng.nextGaussian()*numReadsWithErrors.getStandardDeviation() + numReadsWithErrors.getMean());
            e.numReads = (int) Math.floor(rng.nextGaussian()*numReads.getStandardDeviation() + numReads.getMean());
            e.numFragmentsWithErrors = (int) Math.ceil(rng.nextGaussian()*numFragmentsWithErrors.getStandardDeviation() + numFragmentsWithErrors.getMean());
            e.numFragments = (int) Math.floor(rng.nextGaussian()*numFragments.getStandardDeviation() + numFragments.getMean());

            if (e.numReadsWithErrors < 0) { e.numReadsWithErrors = 0; }
            if (e.numReads < 0) { e.numReads = 0; }
            if (e.numFragmentsWithErrors < 0) { e.numFragmentsWithErrors = 0; }
            if (e.numFragments < 0) { e.numFragments = 0; }

            if (e.numReadsWithErrors > e.numReads) {
                int nr = e.numReads;
                e.numReads = e.numReadsWithErrors;
                e.numReadsWithErrors = nr;
            }

            if (e.numFragmentsWithErrors > e.numFragments) {
                int nf = e.numFragments;
                e.numFragments = e.numFragmentsWithErrors;
                e.numFragmentsWithErrors = nf;
            }

            entries.add(e);
        }

        return entries;
    }

    @Override
    public void execute() {
        if (SEED == null) { SEED = System.currentTimeMillis(); }
        rng = new Random(SEED);

        TableWriter tw = new TableWriter(out);

        log.info("Loading variants...");
        Map<String, Map<Integer, Set<VariantContext>>> variants = loadVariants();

        log.info("Loading table...");
        Map<String, Map<Integer, Entry>> table = loadTable();

        log.info("Processing bases...");
        ReferenceSequence rseq;
        while ((rseq = REF.nextSequence()) != null) {
            String[] names = rseq.getName().split("\\s+");
            String name = names[0];

            log.info("  {}", name);

            String seq = new String(rseq.getBases());
            List<String> alleles = new ArrayList<String>(seq.length());

            for (int i = 0; i < seq.length(); i++) {
                alleles.add(String.valueOf(seq.charAt(i)));
            }

            for (int i = 0; i < seq.length(); i++) {
                if (variants.containsKey(name) && variants.get(name).containsKey(i)) {
                    for (VariantContext vc : variants.get(name).get(i)) {
                        String alleleToRemove = vc.getReference().getBaseString();
                        String alleleToAdd = !vc.isVariant() ? vc.getReference().getBaseString() : vc.getAlternateAllele(0).getBaseString();

                        for (int j = 0; j < alleleToRemove.length(); j++) {
                            alleles.set(i + j, "");
                        }

                        alleles.set(i, alleleToAdd);
                    }
                }
            }

            int gpos = 0;
            for (int i = 0; i < seq.length(); i++) {
                Entry e = table.get(name).get(i);

                List<Entry> es = (alleles.get(i).length() > 1) ? computeEntries(table.get(name), i, alleles.get(i).length()) : null;

                for (int q = 0; q < alleles.get(i).length(); q++) {
                    Map<String, String> liftedTe = new LinkedHashMap<String, String>();

                    liftedTe.put("chr", name);
                    //liftedTe.put("oldpos", String.valueOf(i));
                    liftedTe.put("pos", String.valueOf(gpos + q));
                    //liftedTe.put("q", String.valueOf(q));

                    if (q == 0) {
                        liftedTe.put("numReadsWithErrors", String.valueOf(e.numReadsWithErrors));
                        liftedTe.put("numReads", String.valueOf(e.numReads));
                        liftedTe.put("numFragmentsWithErrors", String.valueOf(e.numFragmentsWithErrors));
                        liftedTe.put("numFragments", String.valueOf(e.numFragments));
                    } else {
                        liftedTe.put("numReadsWithErrors", String.valueOf(es.get(q).numReadsWithErrors));
                        liftedTe.put("numReads", String.valueOf(es.get(q).numReads));
                        liftedTe.put("numFragmentsWithErrors", String.valueOf(es.get(q).numFragmentsWithErrors));
                        liftedTe.put("numFragments", String.valueOf(es.get(q).numFragments));
                    }

                    tw.addEntry(liftedTe);
                }

                gpos += alleles.get(i).length();
            }
        }
    }
}
