package uk.ac.ox.well.indiana.analyses.kmerSharing;

import com.google.common.base.Joiner;
import net.sf.picard.reference.IndexedFastaSequenceFile;
import net.sf.picard.util.Interval;
import uk.ac.ox.well.indiana.tools.Tool;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.io.cortex.CortexGraph;
import uk.ac.ox.well.indiana.utils.io.cortex.CortexRecord;
import uk.ac.ox.well.indiana.utils.io.gff.GFF3;
import uk.ac.ox.well.indiana.utils.io.gff.GFF3Record;
import uk.ac.ox.well.indiana.utils.sequence.SequenceUtils;

import java.io.PrintStream;
import java.util.*;

public class ComputeKmerSharingDistribution extends Tool {
    @Argument(fullName="cortexGraph", shortName="cg", doc="Cortex graph")
    public CortexGraph CORTEX_GRAPH;

    @Argument(fullName="colors", shortName="c", doc="Colors to process")
    public ArrayList<Integer> COLORS;

    @Argument(fullName="gff", shortName="gff", doc="GFF of regions of interest")
    public GFF3 GFF;

    @Argument(fullName="genes", shortName="genes", doc="IDs of genes to evaluate", required=false)
    public ArrayList<String> GENES;

    @Argument(fullName="reference", shortName="R", doc="Reference genome")
    public IndexedFastaSequenceFile FASTA;

    @Argument(fullName="o1", shortName="o1", doc="Sharing between genes")
    public PrintStream o1;

    @Argument(fullName="o2", shortName="o2", doc="Sharing between domains")
    public PrintStream o2;

    //@Output
    //public PrintStream out;

    public ArrayList<String> getGeneList(GFF3 gff, ArrayList<String> genes) {
        ArrayList<String> geneList = genes;

        if (geneList == null) {
            geneList = new ArrayList<String>();

            for (GFF3Record r : gff) {
                if ("gene".equalsIgnoreCase(r.getType())) {
                    geneList.add(r.getAttributes().get("ID"));
                }
            }
        }

        return geneList;
    }

    private class KmerInfo {
        public String gene;
        public String domain;

        public KmerInfo(String gene, String domain) {
            this.gene = gene;
            this.domain = domain;
        }

        public String toString() {
            return this.gene + ":" + this.domain;
        }
    }

    public void execute() {
        HashMap<String, ArrayList<KmerInfo>> kmerInfoMap = new HashMap<String, ArrayList<KmerInfo>>();

        for (String gene : getGeneList(GFF, GENES)) {
            GFF3Record geneRecord = GFF.getRecord(gene);

            log.info("geneRecord={}", geneRecord);

            Collection<GFF3Record> domainRecords = GFF.getType("domain", GFF.getOverlapping(geneRecord));

            for (GFF3Record domainRecord : domainRecords) {
                String seq = new String(FASTA.getSubsequenceAt(domainRecord.getSeqid(), domainRecord.getStart(), domainRecord.getEnd()).getBases());

                for (int i = 0; i < seq.length() - CORTEX_GRAPH.getKmerSize(); i++) {
                    Interval interval = new Interval(domainRecord.getSeqid(), domainRecord.getStart() + i, domainRecord.getStart() + i);

                    String kmer = SequenceUtils.alphanumericallyLowestOrientation(seq.substring(i, i + CORTEX_GRAPH.getKmerSize()));

                    if (!kmerInfoMap.containsKey(kmer)) {
                        kmerInfoMap.put(kmer, new ArrayList<KmerInfo>());
                    }

                    KmerInfo kmerInfo;

                    if (GFF.getType("exon", GFF.getOverlapping(interval)).size() > 0) {
                        kmerInfo = new KmerInfo(geneRecord.getAttributes().get("ID"), domainRecord.getAttributes().get("DOMAIN_CLASS"));
                    } else {
                        kmerInfo = new KmerInfo(geneRecord.getAttributes().get("ID"), "none");
                    }

                    kmerInfoMap.get(kmer).add(kmerInfo);
                }
            }
        }

        TreeMap<Integer, Integer> numGeneHistogram = new TreeMap<Integer, Integer>();
        TreeMap<String, Integer> domainHistogram = new TreeMap<String, Integer>();

        for (CortexRecord cr : CORTEX_GRAPH) {
            int[] coverages = cr.getCoverages();

            boolean hasCoverage = false;
            int numColorsWithKmer = 0;
            boolean hasZeroOrUnitCoverageInColors = true;

            int totalCoverageInROI = 0;

            for (int color : COLORS) {
                hasCoverage |= (coverages[color] > 0);
                numColorsWithKmer += coverages[color] == 0 ? 0 : 1;
                hasZeroOrUnitCoverageInColors &= (coverages[color] <= 1);
                totalCoverageInROI += coverages[color];
            }

            boolean allCoverageIsInROI = (totalCoverageInROI == coverages[0]);

            if (hasCoverage && allCoverageIsInROI && numColorsWithKmer > 1 && hasZeroOrUnitCoverageInColors) {
                String kmer = cr.getKmerAsString();

                if (kmerInfoMap.containsKey(kmer)) {
                    TreeSet<String> domains = new TreeSet<String>();

                    for (KmerInfo kmerInfo : kmerInfoMap.get(kmer)) {
                        domains.add(kmerInfo.domain);
                    }

                    //out.println(kmer + "\t" + numColorsWithKmer + "\t" + Joiner.on(",").join(domains));

                    String domainsString = Joiner.on(",").join(domains);

                    if (!numGeneHistogram.containsKey(numColorsWithKmer)) {
                        numGeneHistogram.put(numColorsWithKmer, 0);
                    }

                    if (!domainHistogram.containsKey(domainsString)) {
                        domainHistogram.put(domainsString, 0);
                    }

                    numGeneHistogram.put(numColorsWithKmer, numGeneHistogram.get(numColorsWithKmer) + 1);
                    domainHistogram.put(domainsString, domainHistogram.get(domainsString) + 1);
                }
            }
        }

        o1.println("num_genes\tcount");
        for (Integer numGene : numGeneHistogram.keySet()) {
            o1.println(numGene + "\t" + numGeneHistogram.get(numGene));
        }

        o1.println();

        o2.println("domain\tcount");
        for (String domain : domainHistogram.keySet()) {
            o2.println(domain + "\t" + domainHistogram.get(domain));
        }

        o2.println();
    }
}
