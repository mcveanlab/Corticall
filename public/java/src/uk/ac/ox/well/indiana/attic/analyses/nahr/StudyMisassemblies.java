package uk.ac.ox.well.indiana.attic.analyses.nahr;

import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMFileReader;
import htsjdk.samtools.SAMRecord;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexKmer;
import uk.ac.ox.well.indiana.utils.io.cortex.graph.CortexRecord;
import uk.ac.ox.well.indiana.utils.io.table.TableReader;

import java.io.File;
import java.io.PrintStream;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

public class StudyMisassemblies extends Module {
    @Argument(fullName="contigs", shortName="c", doc="Contigs (in BAM format)")
    public SAMFileReader CONTIGS;

    @Argument(fullName="reapr", shortName="r", doc="Reaped contigs (in FASTA format)")
    public SAMFileReader REAPR;

    @Argument(fullName="cortexCoverage", shortName="cc", doc="Cortex graph with coverage for sample")
    public CortexGraph COVERAGE;

    @Argument(fullName="parents", shortName="p", doc="Parents (in Cortex format)")
    public CortexGraph PARENTS;

    @Argument(fullName="maskedKmers", shortName="m", doc="Masked kmers")
    public CortexGraph MASKED_KMERS;

    @Argument(fullName="covThresholds", shortName="ct", doc="Coverage thresholds table")
    public File COV_THRESHOLDS;

    @Argument(fullName="otherKs", shortName="ok", doc="Graphs at other kmer sizes")
    public HashSet<CortexGraph> OTHER_KS;

    @Output
    public PrintStream out;

    private int clipLength(SAMRecord read) {
        int cl = 0;

        for (CigarElement ce : read.getCigar().getCigarElements()) {
            if (ce.getOperator().equals(CigarOperator.SOFT_CLIP)) {
                int clnew = ce.getLength();

                if (clnew > cl) {
                    cl = clnew;
                }
            }
        }

        return cl;
    }

    private int numMismatches(SAMRecord read) {
        return (read.getIntegerAttribute("NM") != null) ? read.getIntegerAttribute("NM") : 0;
    }

    private boolean isInteresting(SAMRecord read) {
        return numMismatches(read) > 30 || clipLength(read) > 30;
    }

    private int count(String ann, char b) {
        int num = 0;

        for (int i = 0; i < ann.length(); i++) {
            if (ann.charAt(i) == b) { num++; }
        }

        return num;
    }

    @Override
    public void execute() {
        // Section: load coverage thresholds
        log.info("Load coverage thresholds...");
        int kmerSize = COVERAGE.getKmerSize();
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

        // Section: annotate contigs
        log.info("Annotating contigs...");

        Set<String> reapedContigs = new HashSet<String>();

        int contigsProcessed = 0;
        for (SAMRecord rseq : REAPR) {
            if (rseq.getReadName().matches(".+\\.\\d")) {
                String name = rseq.getReadName().replaceAll("\\.\\d$", "").replaceAll("_", "-");

                reapedContigs.add(name);
            }
        }

        int contigsAccepted = 0;
        for (SAMRecord rseq : CONTIGS) {
            if (reapedContigs.contains(rseq.getReadName()) || isInteresting(rseq)) {
                String seq = rseq.getReadString();
                StringBuilder annotation = new StringBuilder();

                contigsAccepted++;

                for (int i = 0; i <= seq.length() - kmerSize; i++) {
                    String sk = seq.substring(i, i + kmerSize);
                    CortexKmer ck = new CortexKmer(sk);

                    CortexRecord crMasked = MASKED_KMERS.findRecord(ck);
                    CortexRecord crSample = COVERAGE.findRecord(ck);
                    CortexRecord crParent = PARENTS.findRecord(ck);

                    boolean isMasked = crMasked != null;
                    int cov = crSample != null ? crSample.getCoverage(0) : 0;
                    int cov0 = crParent != null ? crParent.getCoverage(0) : 0;
                    int cov1 = crParent != null ? crParent.getCoverage(1) : 0;
                    boolean missingAtHigherK = false;
                    boolean missingAtLowerK = false;

                    for (CortexGraph cg : OTHER_KS) {
                        if (cg.getKmerSize() != kmerSize && i <= seq.length() - cg.getKmerSize()) {
                            String sko = seq.substring(i, i + cg.getKmerSize());
                            CortexKmer cko = new CortexKmer(sko);

                            CortexRecord crOther = cg.findRecord(cko);

                            if (crOther == null) {
                                if (cg.getKmerSize() > kmerSize) {
                                    missingAtHigherK = true;
                                } else {
                                    missingAtLowerK = true;
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
                    }
                }

                String ann = annotation.toString();

                /*
                log.info("name={}, loc={}:{}, nm={}, cl={}, R={}, C={}, M={}, m={}, 0={}, 1={}, B={}, _={}, .={}",
                        rseq.getReadName(),
                        rseq.getReferenceName(),
                        rseq.getAlignmentStart(),
                        numMismatches(rseq),
                        clipLength(rseq),
                        count(ann, 'R'),
                        count(ann, 'C'),
                        count(ann, 'M'),
                        count(ann, 'm'),
                        count(ann, '0'),
                        count(ann, '1'),
                        count(ann, 'B'),
                        count(ann, '_'),
                        count(ann, '.')
                );
                */

                log.info("name={}, loc={}:{}, nm={}, cl={}, R={}, M={}, m={}",
                        rseq.getReadName(),
                        rseq.getReferenceName(),
                        rseq.getAlignmentStart(),
                        numMismatches(rseq),
                        clipLength(rseq),
                        count(ann, 'R'),
                        count(ann, 'M'),
                        count(ann, 'm')
                );
            }
        }

        log.info("  {} contigs", contigsProcessed);
        log.info("  {} reaped", reapedContigs.size());
    }
}
