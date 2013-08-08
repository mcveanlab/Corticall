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
     * Returns a String version of a nucleotide byte
     *
     * @param nucleotide  the byte version of the nucleotide
     * @return  the String version of the nucleotide
     */
    public static String nucleotideByteToString(byte nucleotide) {
        String base = "N";

        switch (nucleotide) {
            case 'A': base = "A"; break;
            case 'a': base = "a"; break;

            case 'C': base = "C"; break;
            case 'c': base = "c"; break;

            case 'G': base = "G"; break;
            case 'g': base = "g"; break;

            case 'T': base = "T"; break;
            case 't': base = "t"; break;

            case 'N': base = "N"; break;
            case 'n': base = "n"; break;

            case '.': base = "."; break;
        }

        return base;
    }
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
    private static int computeN50Metric(Collection<? extends CharSequence> sequences, boolean lengthOrValue) {
        int totalLength = 0;

        List<Integer> lengths = new ArrayList<Integer>();
        for (CharSequence seq : sequences) {
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
    public static int computeN50Length(Collection<? extends CharSequence> sequences) {
        return computeN50Metric(sequences, true);
    }

    /**
     * Compute N50 value (the longest length for which the collection of all contigs of that length or longer contains at least half of the total of the lengths of the contigs)
     * @param sequences  the sequences to process
     * @return  the N50 value
     */
    public static int computeN50Value(Collection<? extends CharSequence> sequences) {
        return computeN50Metric(sequences, false);
    }

    /**
     * Compute the maximum length in a collection of sequences
     *
     * @param sequences  the sequences to process
     * @return  the maximum sequence length observed
     */
    public static int maxLength(Collection<? extends CharSequence> sequences) {
        int length = 0;

        for (CharSequence seq : sequences) {
            if (seq.length() > length) {
                length = seq.length();
            }
        }

        return length;
    }

    private static final byte[] nucleotides = {'A', 'C', 'G', 'T'};
    private static Random randomGenerator = new Random(0);

    /**
     * Create a random nucleotide sequence of length N
     *
     * @param n  the length of the sequence to generate
     * @return  the generated sequence
     */
    public static byte[] generateRandomNucleotideSequenceOfLengthN(int n) {
        byte[] sb = new byte[n];

        for (int i = 0; i < n; i++) {
            sb[i] = nucleotides[randomGenerator.nextInt(nucleotides.length)];
        }

        return sb;
    }

    public static Collection<byte[]> generateSequencesWithEditDistance1(byte[] source) {
        Collection<byte[]> dest = new HashSet<byte[]>();

        for (int pos = 0; pos < source.length; pos++) {
            for (byte nucleotide : nucleotides) {
                if (source[pos] != nucleotide) {
                    byte[] newSeq = new byte[source.length];
                    System.arraycopy(source, 0, newSeq, 0, source.length);
                    newSeq[pos] = nucleotide;

                    dest.add(newSeq);
                }
            }
        }

        return dest;
    }

    public static Collection<byte[]> generateSequencesWithEditDistance2(byte[] source) {
        Collection<byte[]> dest = new LinkedHashSet<byte[]>();

        for (int pos1 = 0; pos1 < source.length; pos1++) {
            for (byte nucleotide1 : nucleotides) {
                if (source[pos1] != nucleotide1) {
                    for (int pos2 = 0; pos2 < source.length; pos2++) {
                        if (pos2 != pos1) {
                            for (byte nucleotide2 : nucleotides) {
                                if (source[pos2] != nucleotide2) {
                                    byte[] newSeq = new byte[source.length];
                                    System.arraycopy(source, 0, newSeq, 0, source.length);

                                    newSeq[pos1] = nucleotide1;
                                    newSeq[pos2] = nucleotide2;

                                    dest.add(newSeq);
                                }
                            }
                        }
                    }
                }
            }
        }

        return dest;
    }

    public static int editDistance(byte[] s1, byte[] s2) {
        int editDistance = 0;

        for (int i = 0; i < s1.length; i++) {
            if (s1[i] != s2[i]) {
                editDistance++;
            }
        }

        return editDistance;
    }
}
