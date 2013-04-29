package uk.ac.ox.well.indiana.utils.sequence;

import net.sf.picard.reference.FastaSequenceFile;
import net.sf.picard.reference.ReferenceSequence;

import java.util.*;

/**
 * A set of utilities for dealing with genomic sequences.
 */
public class SequenceUtils {
    /**
     * Private constructor - this class cannot be instantiated!
     */
    private SequenceUtils() {}

    /**
     * Get the complement of a single nucleotide
     *
     * @param nucleotide  the nucleotide that should be complemented
     * @return  the complement of the specified nucleotide
     */
    public static byte complement(byte nucleotide) {
        byte rcBase = 'N';

        switch(nucleotide) {
            case 'A': rcBase = 'T'; break;
            case 'a': rcBase = 't'; break;

            case 'C': rcBase = 'G'; break;
            case 'c': rcBase = 'g'; break;

            case 'G': rcBase = 'C'; break;
            case 'g': rcBase = 'c'; break;

            case 'T': rcBase = 'A'; break;
            case 't': rcBase = 'a'; break;

            case 'N': rcBase = 'N'; break;
            case 'n': rcBase = 'n'; break;

            case '.': rcBase = '.'; break;
        }

        return rcBase;
    }

    /**
     * Get the reverse complement of the sequence.
     *
     * @param sequence  the sequence that should be reverse complemented
     * @return  the reverse complement of the sequence
     */
    public static byte[] reverseComplement(byte[] sequence) {
        byte[] rc = new byte[sequence.length];

        for (int i = 0; i < sequence.length; i++) {
            rc[sequence.length - 1 - i] = complement(sequence[i]);
        }

        return rc;
    }

    /**
     * Get the reverse complement of the sequence.
     *
     * @param sequence  the sequence that should be reverse complemented
     * @return  the reverse complement of the sequence
     */
    public static String reverseComplement(String sequence) {
        return new String(reverseComplement(sequence.getBytes()));
    }

    /**
     * Get the alphanumerically lowest orientation of the specified sequence.
     *
     * @param sequence  the sequence to process
     * @return  the alphanumerically lowest orientation of the given sequence
     */
    public static byte[] alphanumericallyLowestOrientation(byte[] sequence) {
        byte[] rc = reverseComplement(sequence);

        String kmerStr = new String(sequence);
        String rcStr = new String(rc);

        return (kmerStr.compareTo(rcStr) < 0) ? sequence : rc;
    }

    /**
     * Get the alphanumerically lowest orientation of the specified sequence.
     *
     * @param sequence  the sequence to process
     * @return  the alphanumerically lowest orientation of the given sequence
     */
    public static String alphanumericallyLowestOrientation(String sequence) {
        return new String(alphanumericallyLowestOrientation(sequence.getBytes()));
    }

    /**
     * For given FASTA files, load all of the sequences as the hashcodes of the alphanumerically lowest kmers mapped to the contig name
     *
     * @param fastas  list of FASTA files to process
     * @param kmerSize  the kmer size to use
     * @return  a map of kmer hashcodes to sequence names
     */
    public static Map<Integer, String> loadSequenceCodesAsAlphanumericallyLowestKmers(List<FastaSequenceFile> fastas, int kmerSize) {
        Map<Integer, String> kmerHash = new HashMap<Integer, String>();

        for (FastaSequenceFile fasta : fastas) {
            kmerHash.putAll(loadSequenceCodesAsAlphanumericallyLowestKmers(fasta, kmerSize));
        }

        return kmerHash;
    }

    /**
     * For given FASTA files, load all of the sequences as the hashcodes of the alphanumerically lowest kmers mapped to the contig name
     *
     * @param fasta  FASTA file to process
     * @param kmerSize  the kmer size to use
     * @return  a map of kmer hashcodes to sequence names
     */
    public static Map<Integer, String> loadSequenceCodesAsAlphanumericallyLowestKmers(FastaSequenceFile fasta, int kmerSize) {
        Map<Integer, String> kmerHash = new HashMap<Integer, String>();

        ReferenceSequence seq;
        while ((seq = fasta.nextSequence()) != null) {
            kmerHash.putAll(loadSequenceCodesAsAlphanumericallyLowestKmers(seq, kmerSize));
        }

        return kmerHash;
    }

    /**
     * For given FASTA files, load all of the sequences as the hashcodes of the alphanumerically lowest kmers mapped to the contig name
     *
     * @param seq  the reference sequence from the FASTA file
     * @param kmerSize  the kmer size to use
     * @return  a map of kmer hashcodes to sequence names
     */
    public static Map<Integer, String> loadSequenceCodesAsAlphanumericallyLowestKmers(ReferenceSequence seq, int kmerSize) {
        Map<Integer, String> kmerHash = new HashMap<Integer, String>();

        byte[] contig = seq.getBases();
        String[] name = seq.getName().split("\\s+");

        for (int i = 0; i < contig.length - kmerSize; i++) {
            String kmer = new String(alphanumericallyLowestOrientation(Arrays.copyOfRange(contig, i, i + kmerSize)));

            kmerHash.put(kmer.hashCode(), name[0]);
        }

        return kmerHash;
    }

    /**
     * Compute the number of segregating sites between two sequences.
     * @param s1  sequence 1
     * @param s2  sequence 2
     * @return  the number of segregating sites
     */
    public static int numSegregatingSites(String s1, String s2) {
        if (s1.length() != s2.length()) {
            throw new RuntimeException("Cannot compute number of segreating sites between sequences of different lengths (" + s1.length() + " vs " + s2.length() + ")");
        }

        int S = 0;
        for (int i = 0; i < s1.length(); i++) {
            S += (s1.charAt(i) != s2.charAt(i)) ? 1 : 0;
        }

        return S;
    }

    /**
     * Private method to compute N50 information
     *
     * @param sequences  the sequences to process
     * @param lengthOrValue  if true, compute length; else, compute value
     * @return  the N50 length or value
     */
    private static int computeN50Metric(Collection<String> sequences, boolean lengthOrValue) {
        int totalLength = 0;

        List<Integer> lengths = new ArrayList<Integer>();
        for (String seq : sequences) {
            totalLength += seq.length();

            lengths.add(seq.length());
        }

        Collections.sort(lengths, new Comparator<Integer>() {
            @Override
            public int compare(Integer o1, Integer o2) {
                return o2.compareTo(o1);
            }
        });

        int n50Length = 0;
        int n50Value = 0;

        for (Integer length : lengths) {
            n50Length += length;
            n50Value = length;

            if (n50Length >= totalLength/2) {
                break;
            }
        }

        return lengthOrValue ? n50Length : n50Value;
    }

    /**
     * Compute N50 length
     * @param sequences  the sequences to process
     * @return  the N50 length
     */
    public static int computeN50Length(Collection<String> sequences) {
        return computeN50Metric(sequences, true);
    }

    /**
     * Compute N50 value
     * @param sequences  the sequences to process
     * @return  the N50 value
     */
    public static int computeN50Value(Collection<String> sequences) {
        return computeN50Metric(sequences, false);
    }
}
