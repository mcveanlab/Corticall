package uk.ac.ox.well.indiana.attic.analyses.nahr;

import com.google.common.base.Joiner;
import htsjdk.samtools.reference.FastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.cortex.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.CortexKmer;
import uk.ac.ox.well.indiana.utils.io.cortex.CortexMap;
import uk.ac.ox.well.indiana.utils.io.cortex.CortexRecord;
import uk.ac.ox.well.indiana.utils.io.table.TableReader;
import uk.ac.ox.well.indiana.utils.io.table.TableWriter;

import java.io.File;
import java.io.PrintStream;
import java.util.*;

public class AnnotateContigsQuickly extends Module {
    @Argument(fullName="contigs", shortName="c", doc="Contigs (in FASTA format)")
    public FastaSequenceFile CONTIGS;

    @Argument(fullName="cortexCoverage", shortName="cc", doc="Cortex graph with coverage for sample")
    public CortexGraph COVERAGE;

    @Argument(fullName="parents", shortName="p", doc="Parents (in Cortex format)")
    public CortexGraph PARENTS;

    @Argument(fullName="maskedKmers", shortName="m", doc="Masked kmers")
    public CortexGraph MASKED_KMERS;

    @Argument(fullName="covThresholds", shortName="ct", doc="Coverage thresholds table")
    public File COV_THRESHOLDS;

    @Output
    public PrintStream out;

    private class KmerInfo {
        public boolean isMasked = false;
        public int coverageInSample = 0;
        public int coverageInParent0 = 0;
        public int coverageInParent1 = 0;
    }

    @Override
    public void execute() {
        int kmerSize = COVERAGE.getKmerSize();

        // Section: load coverage thresholds
        log.info("Load coverage thresholds...");
        int sampleP0Threshold = 0;
        int sampleP1Threshold = 0;
        int p1Threshold = 0;
        int p2Threshold = 0;

        TableReader tr = new TableReader(COV_THRESHOLDS);
        for (Map<String, String> te : tr) {
            String sample = te.get("sample");

            if (COVERAGE.getColor(0).getSampleName().contains(sample)) {
                sampleP0Threshold = Integer.valueOf(te.get("p1"));
                sampleP1Threshold = Integer.valueOf(te.get("p2"));
            } else if (PARENTS.getColor(0).getSampleName().contains(sample)) {
                p1Threshold = Integer.valueOf(te.get("p1"));
            } else if (PARENTS.getColor(1).getSampleName().contains(sample)) {
                p2Threshold = Integer.valueOf(te.get("p2"));
            }
        }

        int sampleThreshold = sampleP0Threshold > sampleP1Threshold ? sampleP0Threshold : sampleP1Threshold;

        log.info("  sampleP0Threshold: {}", sampleP0Threshold);
        log.info("  sampleP1Threshold: {}", sampleP1Threshold);
        log.info("    sampleThreshold: {}", sampleThreshold);
        log.info("        p1Threshold: {}", p1Threshold);
        log.info("        p2Threshold: {}", p2Threshold);

        // Section: load contigs
        log.info("Loading contigs...");
        List<ReferenceSequence> contigs = new ArrayList<ReferenceSequence>();
        HashMap<CortexKmer, KmerInfo> contigKmers = new HashMap<CortexKmer, KmerInfo>();

        ReferenceSequence rseqa;
        while ((rseqa = CONTIGS.nextSequence()) != null) {
            contigs.add(rseqa);

            String seqa = new String(rseqa.getBases());
            for (int i = 0; i <= seqa.length() - kmerSize; i++) {
                CortexKmer kmera = new CortexKmer(seqa.substring(i, i + kmerSize));

                contigKmers.put(kmera, new KmerInfo());
            }
        }

        log.info("  {} contigs loaded", contigs.size());
        log.info("  {} unique contig kmers", contigKmers.size());

        // Section: fill in kmer attributes (mask status)
        log.info("Loading masked kmers...");
        int maskCount = 0;
        for (CortexRecord cr : MASKED_KMERS) {
            if (contigKmers.containsKey(cr.getKmer())) {
                contigKmers.get(cr.getKmer()).isMasked = true;

                maskCount++;
            }
        }

        log.info("  {} masked kmers loaded", maskCount);

        // Section: fill in kmer attributes (kmer coverage in sample)
        log.info("Loading kmer coverage information from sample...");
        int recIndex = 0;
        int recWithCoverageInfo = 0;
        for (CortexRecord cr : COVERAGE) {
            if (recIndex % (COVERAGE.getNumRecords() / 10) == 0) {
                log.info("  {}/{} records", recIndex, COVERAGE.getNumRecords());
            }
            recIndex++;

            if (contigKmers.containsKey(cr.getKmer())) {
                contigKmers.get(cr.getKmer()).coverageInSample = cr.getCoverage(0);

                recWithCoverageInfo++;
            }
        }

        log.info("  processed {} records, {} with relevant coverage information", recIndex, recWithCoverageInfo);

        // Section: fill in kmer attributes (kmer coverage in parents)
        log.info("Loading kmer coverage information from parents...");
        recIndex = 0;
        recWithCoverageInfo = 0;
        for (CortexRecord cr : PARENTS) {
            if (recIndex % (PARENTS.getNumRecords() / 10) == 0) {
                log.info("  {}/{} records", recIndex, PARENTS.getNumRecords());
            }
            recIndex++;

            if (contigKmers.containsKey(cr.getKmer())) {
                int cov0 = cr.getCoverage(0);
                int cov1 = cr.getCoverage(1);

                contigKmers.get(cr.getKmer()).coverageInParent0 = cov0;
                contigKmers.get(cr.getKmer()).coverageInParent1 = cov1;

                recWithCoverageInfo++;
            }
        }

        log.info("  processed {} records, {} with relevant coverage information", recIndex, recWithCoverageInfo);

        // Section: write annotation information for each contig
        log.info("Writing contig annotations...");

        out.println("contigName\tseq\tkmerOrigin\tkmerContiguity\tkmerCoverage");

        int contigsSeen = 0;
        for (ReferenceSequence rseq : contigs) {
            if (contigsSeen % (contigs.size() / 10) == 0) {
                log.info("  annotated {}/{} contigs", contigsSeen, contigs.size());
            }
            contigsSeen++;

            if (rseq.length() > 0) {
                String seq = new String(rseq.getBases());

                StringBuilder annotation = new StringBuilder();
                List<String> coverages = new ArrayList<String>();

                for (int i = 0; i <= seq.length() - kmerSize; i++) {
                    CortexKmer kmer = new CortexKmer(seq.substring(i, i + kmerSize));

                    boolean isMasked = contigKmers.get(kmer).isMasked;
                    int cov = contigKmers.get(kmer).coverageInSample;
                    int cov0 = contigKmers.get(kmer).coverageInParent0;
                    int cov1 = contigKmers.get(kmer).coverageInParent1;

                    coverages.add(String.valueOf(cov));

                    if (isMasked) {
                        annotation.append("A");
                    } else if (cov < sampleThreshold) {
                        annotation.append("C");
                    } else if (cov0 > 0 || cov1 > 0) {
                        if (cov0 == 0 && cov1 >= p2Threshold) {
                            annotation.append("1");
                        } else if (cov0 >= p1Threshold && cov1 == 0) {
                            annotation.append("0");
                        } else if (cov0 >= p1Threshold && cov1 >= p2Threshold) {
                            annotation.append("B");
                        } else {
                            annotation.append("_");
                        }
                    } else {
                        annotation.append(".");
                    }
                }

                StringBuilder contiguity = new StringBuilder();
                contiguity.append("0");
                for (int i = 1; i <= seq.length() - kmerSize; i++) {
                    /*
                    String prevStr = seq.substring(i - 1, i - 1 + kmerSize);
                    String curStr  = seq.substring(i, i + kmerSize);

                    char colorChar = annotation.charAt(i);

                    boolean isContiguous = true;

                    switch (colorChar) {
                        case '0': isContiguous = isContiguous(prevStr, curStr, 0); break;
                        case '1': isContiguous = isContiguous(prevStr, curStr, 1); break;
                        case 'B': isContiguous = isContiguous(prevStr, curStr, 0) || isContiguous(prevStr, curStr, 1); break;
                        default:  isContiguous = false; break;
                    }

                    contiguity.append(isContiguous ? "1" : "0");
                    */

                    contiguity.append("1");
                }

                String coverage = Joiner.on(",").join(coverages);

                if (annotation.length() == 0) { annotation.append("."); }
                if (coverage.length() == 0) { coverage = "0"; }

                //if (seq.length() >= kmerSize) {
                out.println(rseq.getName().split("\\s+")[0] + "\t" + seq + "\t" + annotation.toString() + "\t" + contiguity.toString() + "\t" + coverage);
                //}
            }
        }
    }
}
