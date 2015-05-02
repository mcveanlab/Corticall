package uk.ac.ox.well.indiana.utils.io.fasta;

import htsjdk.samtools.reference.FastaSequenceIndex;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.util.Interval;
import uk.ac.ox.well.indiana.Indiana;
import uk.ac.ox.well.indiana.utils.exceptions.IndianaException;
import uk.ac.ox.well.indiana.utils.io.kmerindex.KmerIndex;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.HashMap;
import java.util.Map;

public class MultiIndexedFastaSequenceFile extends IndexedFastaSequenceFile {
    private File file;
    private Map<Integer, KmerIndex> kmerIndices = new HashMap<Integer, KmerIndex>();

    public MultiIndexedFastaSequenceFile(File file, FastaSequenceIndex index) {
        super(file, index);
        this.file = file;
    }

    public MultiIndexedFastaSequenceFile(File file) throws FileNotFoundException {
        super(file);
        this.file = file;
    }

    private void loadKmerIndex(int kmerSize) {
        if (!kmerIndices.containsKey(kmerSize)) {
            File indexFile = new File(file.getAbsolutePath() + ".k" + kmerSize + ".idx");

            if (!indexFile.exists()) {
                throw new IndianaException("Index for kmer size " + kmerSize + " not found.");
            }

            KmerIndex ki = new KmerIndex(indexFile);

            kmerIndices.put(kmerSize, ki);
        }
    }

    public Interval find(String kmer) {
        loadKmerIndex(kmer.length());

        KmerIndex ki = kmerIndices.get(kmer.length());

        return ki.findLocus(kmer);
    }
}
