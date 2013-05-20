package uk.ac.ox.well.indiana.analyses.reconstruction;

import cern.colt.matrix.DoubleMatrix1D;
import com.google.common.base.Joiner;
import uk.ac.ox.well.indiana.tools.Tool;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;
import uk.ac.ox.well.indiana.utils.containers.DataFrame;
import uk.ac.ox.well.indiana.utils.io.utils.LineReader;
import uk.ac.ox.well.indiana.utils.io.utils.TableReader2;
import uk.ac.ox.well.indiana.utils.statistics.PCA;

import java.io.File;
import java.io.PrintStream;
import java.util.*;

public class ComputeSampleRelatedness extends Tool {
    @Argument(fullName="kmerSharingMatrix", shortName="ksm", doc="Kmer sharing matrix")
    public File KMER_SHARING_MATRIX;

    @Argument(fullName="geneLists", shortName="gl", doc="Gene lists")
    public ArrayList<File> GENE_LISTS;

    @Argument(fullName="geneClasses", shortName="gc", doc="Gene classes")
    public ArrayList<String> GENE_CLASSES;

    @Output
    public PrintStream out;

    private Map<String, String> loadGeneListsAndClasses() {
        Map<String, String> geneToClass = new HashMap<String, String>();

        for (int i = 0; i < GENE_LISTS.size(); i++) {
            LineReader lr = new LineReader(GENE_LISTS.get(i));

            String line;
            while ((line = lr.getNextRecord()) != null) {
                geneToClass.put(line, GENE_CLASSES.get(i));
            }
        }

        return geneToClass;
    }

    private String getGeneClass(String commaDelimitedGenes, Map<String, String> geneClasses) {
        String[] genes = commaDelimitedGenes.split(",");

        Set<String> classes = new TreeSet<String>();
        for (String gene : genes) {
            classes.add(geneClasses.get(gene));
        }

        return Joiner.on(",").join(classes);
    }

    @Override
    public void execute() {
        Map<String, String> geneToClass = loadGeneListsAndClasses();

        Map<String, DataFrame<String, String, Float>> d = new HashMap<String, DataFrame<String, String, Float>>();
        Set<String> samples = new HashSet<String>();
        int count = 0;

        TableReader2 tr = new TableReader2(KMER_SHARING_MATRIX);
        for (Map<String, String> te : tr) {
            if (samples.isEmpty()) {
                for (String field : te.keySet()) {
                    if (!field.equals("kmer") && !field.equals("isPanelKmer") && !field.equals("genes")) {
                        samples.add(field);
                    }
                }
            }

            String geneClass = getGeneClass(te.get("genes"), geneToClass);
            if (!d.containsKey(geneClass)) {
                d.put(geneClass, new DataFrame<String, String, Float>(0.0f));
            }

            for (String field1 : te.keySet()) {
                if (!field1.equals("kmer") && !field1.equals("isPanelKmer") && !field1.equals("genes")) {
                    for (String field2 : te.keySet()) {
                        if (!field2.equals("kmer") && !field2.equals("isPanelKmer") && !field2.equals("genes")) {
                            if (te.get(field1).equals("1") && te.get(field2).equals("1")) {
                                float value = d.get(geneClass).get(field1, field2);
                                d.get(geneClass).set(field1, field2, value + 1.0f);
                            }
                        }
                    }
                }
            }


            count++;
            if (count % (tr.size() / 5) == 0) {
                log.info("Processed {}/{} records", count, tr.size());
            }
        }

        out.printf("%s\t%s\t%s\t%s\t%s\t%s\n", "sample", "geneClass", "kmerCount", "PC1", "PC2", "PC3");

        for (String geneClass : d.keySet()) {
            log.info("Performing PCA on {}", geneClass);

            PCA<String, String> pca = new PCA<String, String>(d.get(geneClass));

            for (String sample : samples) {
                DoubleMatrix1D ev = pca.getEigenvector(sample);

                out.printf("%s\t%s\t%d\t%f\t%f\t%f\n", sample, geneClass, d.get(geneClass).getNumRows(), ev.get(0), ev.get(1), ev.get(2));
            }
        }
    }
}
