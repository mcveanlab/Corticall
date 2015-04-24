package uk.ac.ox.well.indiana.attic.analyses.nahr;

import com.google.common.base.Joiner;
import htsjdk.samtools.reference.FastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.*;
import uk.ac.ox.well.indiana.utils.io.table.TableReader;

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

    @Argument(fullName="otherKs", shortName="ok", doc="Graphs at other kmer sizes", required=false)
    public HashSet<CortexGraph> OTHER_KS;

    @Output
    public PrintStream out;

    @Output(fullName="novelKmerList", shortName="ko", doc="List of kmers considered novel")
    public PrintStream kout;

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

            int p1 = Integer.valueOf(te.get("p1"));
            int p2 = Integer.valueOf(te.get("p2"));

            if (COVERAGE.getColor(0).getSampleName().contains(sample)) {
                sampleP0Threshold = p1;
                sampleP1Threshold = p2;
            } else if (p1 >  0 && p2 == 0) {
                p1Threshold = p1;
            } else if (p1 == 0 && p2 >  0) {
                p2Threshold = p2;
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
            if (contigKmers.containsKey(cr.getCortexKmer())) {
                contigKmers.get(cr.getCortexKmer()).isMasked = true;

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

            if (contigKmers.containsKey(cr.getCortexKmer())) {
                contigKmers.get(cr.getCortexKmer()).coverageInSample = cr.getCoverage(0);

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

            if (contigKmers.containsKey(cr.getCortexKmer())) {
                int cov0 = cr.getCoverage(0);
                int cov1 = cr.getCoverage(1);

                contigKmers.get(cr.getCortexKmer()).coverageInParent0 = cov0;
                contigKmers.get(cr.getCortexKmer()).coverageInParent1 = cov1;

                recWithCoverageInfo++;
            }
        }

        log.info("  processed {} records, {} with relevant coverage information", recIndex, recWithCoverageInfo);

        // Section: write annotation information for each contig
        log.info("Writing contig annotations...");

        out.println("contigName\tseq\tkmerOrigin\tkmerCoverage");

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
                    boolean missingAtHigherK = false;
                    boolean missingAtLowerK = false;

                    coverages.add(String.valueOf(cov));

                    if (OTHER_KS != null && !OTHER_KS.isEmpty()) {
                        for (CortexGraph og : OTHER_KS) {
                            if (og.getKmerSize() != kmerSize && i <= seq.length() - og.getKmerSize()) {
                                String sk = seq.substring(i, i + og.getKmerSize());
                                CortexKmer ck = new CortexKmer(sk);
                                CortexRecord ogr = og.findRecord(ck);

                                if (ogr == null) {
                                    if (og.getKmerSize() > kmerSize) {
                                        missingAtHigherK = true;
                                    } else if (og.getKmerSize() < kmerSize) {
                                        missingAtLowerK = true;
                                    }

                                    break;
                                }
                            }
                        }
                    }

                    if (isMasked) {
                        annotation.append("R");     // masked
                    } else if (cov < sampleThreshold) {
                        annotation.append("C");     // insufficient coverage in sample
                    } else if (missingAtHigherK) {
                        annotation.append("M");     // kmer starting here is missing at higher K
                    } else if (missingAtLowerK) {
                        annotation.append("m");     // kmer starting here is missing at lower K
                    } else if (cov0 > 0 || cov1 > 0) {
                        if (cov0 == 0 && cov1 >= p2Threshold) {
                            annotation.append("1"); // confidently parent 1
                        } else if (cov0 >= p1Threshold && cov1 == 0) {
                            annotation.append("0"); // confidently parent 0
                        } else if (cov0 >= p1Threshold && cov1 >= p2Threshold) {
                            annotation.append("B"); // confidently shared
                        } else {
                            annotation.append("_"); // probably parental but could not be confidently determined
                        }
                    } else {
                        annotation.append(".");     // novel
                        kout.println(kmer.getKmerAsString());
                    }
                }

                String coverage = Joiner.on(",").join(coverages);

                if (annotation.length() == 0) { annotation.append("."); }
                if (coverage.length() == 0) { coverage = "0"; }

                out.println(rseq.getName().split("\\s+")[0] + "\t" + seq + "\t" + annotation.toString() + "\t" + coverage);
            }
        }
    }
}
