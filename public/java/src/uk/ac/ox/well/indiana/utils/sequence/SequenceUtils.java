package uk.ac.ox.well.indiana.utils.sequence;

import htsjdk.samtools.reference.FastaSequenceFile;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import uk.ac.ox.well.indiana.utils.io.gff.GFF3Record;

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
     * Get the complement of a single nucleotide
     *
     * @param nucleotide  the nucleotide that should be complemented
     * @return  the complement of the specified nucleotide
     */
    public static char complement(char nucleotide) {
        char rcBase = 'N';

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
     * Get the reverse of the sequence.
     *
     * @param sequence  the sequence that should be reverse complemented
     * @return  the reverse complement of the sequence
     */
    public static byte[] reverse(byte[] sequence) {
        byte[] rc = new byte[sequence.length];

        for (int i = 0; i < sequence.length; i++) {
            rc[sequence.length - 1 - i] = sequence[i];
        }

        return rc;
    }

    /**
     * Get the reverse of the sequence.
     *
     * @param sequence  the sequence that should be reverse complemented
     * @return  the reverse complement of the sequence
     */
    public static String reverse(String sequence) {
        return new String(reverse(sequence.getBytes()));
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
     * Get the complement of the sequence.
     *
     * @param sequence  the sequence that should be complemented
     * @return  the complement of the sequence
     */
    public static byte[] complement(byte[] sequence) {
        byte[] c = new byte[sequence.length];

        for (int i = 0; i < sequence.length; i++) {
            c[i] = complement(sequence[i]);
        }

        return c;
    }

    /**
     * Get the complement of the sequence.
     *
     * @param sequence  the sequence that should be complemented
     * @return  the complement of the sequence
     */
    public static String complement(String sequence) {
        return new String(complement(sequence.getBytes()));
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
        float totalLength = 0;

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

        float n50Length = 0;
        float n50Value = 0;

        for (Integer length : lengths) {
            n50Length += length;
            n50Value = length;

            if (n50Length >= totalLength/2) {
                break;
            }
        }

        return lengthOrValue ? (int) n50Length : (int) n50Value;
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
     * Compute the minimum length in a collection of sequences
     *
     * @param sequences  the sequences to process
     * @return  the minimum sequence length observed
     */
    public static int minLength(Collection<? extends CharSequence> sequences) {
        int length = Integer.MAX_VALUE;

        for (CharSequence seq : sequences) {
            if (seq.length() < length) {
                length = seq.length();
            }
        }

        return length;
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

    /**
     * Compute the mean length in a collection of sequences
     *
     * @param sequences  the sequences to process
     * @return  the mean sequence length observed
     */
    public static float meanLength(Collection<? extends CharSequence> sequences) {
        float length = 0;
        float numSequences = 0;

        for (CharSequence seq : sequences) {
            length += seq.length();
            numSequences++;
        }

        return length / numSequences;
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

    private static final Map<String, String> codonToAminoAcidMap;
    static {
        Map<String, String> c2a = new HashMap<String, String>();
        c2a.put("GCT", "A");
        c2a.put("GCC", "A");
        c2a.put("GCA", "A");
        c2a.put("GCG", "A");

        c2a.put("CGT", "R");
        c2a.put("CGC", "R");
        c2a.put("CGA", "R");
        c2a.put("CGG", "R");
        c2a.put("AGA", "R");
        c2a.put("AGG", "R");

        c2a.put("AAT", "N");
        c2a.put("AAC", "N");

        c2a.put("GAT", "D");
        c2a.put("GAC", "D");

        c2a.put("TGT", "C");
        c2a.put("TGC", "C");

        c2a.put("CAA", "Q");
        c2a.put("CAG", "Q");

        c2a.put("GAA", "E");
        c2a.put("GAG", "E");

        c2a.put("GGT", "G");
        c2a.put("GGC", "G");
        c2a.put("GGA", "G");
        c2a.put("GGG", "G");

        c2a.put("CAT", "H");
        c2a.put("CAC", "H");

        c2a.put("ATT", "I");
        c2a.put("ATC", "I");
        c2a.put("ATA", "I");

        c2a.put("TTA", "L");
        c2a.put("TTG", "L");
        c2a.put("CTT", "L");
        c2a.put("CTC", "L");
        c2a.put("CTA", "L");
        c2a.put("CTG", "L");

        c2a.put("AAA", "K");
        c2a.put("AAG", "K");

        c2a.put("ATG", "M");

        c2a.put("TTT", "F");
        c2a.put("TTC", "F");

        c2a.put("CCT", "P");
        c2a.put("CCC", "P");
        c2a.put("CCA", "P");
        c2a.put("CCG", "P");

        c2a.put("TCT", "S");
        c2a.put("TCC", "S");
        c2a.put("TCA", "S");
        c2a.put("TCG", "S");
        c2a.put("AGT", "S");
        c2a.put("AGC", "S");

        c2a.put("ACT", "T");
        c2a.put("ACC", "T");
        c2a.put("ACA", "T");
        c2a.put("ACG", "T");

        c2a.put("TGG", "W");

        c2a.put("TAT", "Y");
        c2a.put("TAC", "Y");

        c2a.put("GTT", "V");
        c2a.put("GTC", "V");
        c2a.put("GTA", "V");
        c2a.put("GTG", "V");

        c2a.put("TAA", "*");
        c2a.put("TGA", "*");
        c2a.put("TAG", "*");

        codonToAminoAcidMap = Collections.unmodifiableMap(c2a);
    }

    public static String codonToAminoAcid(String codon) {
        if (codon.contains("N") && codon.length() == 3) {
            return "N";
        } else if (codonToAminoAcidMap.containsKey(codon)) {
            return codonToAminoAcidMap.get(codon);
        } else {
            throw new RuntimeException("Did not find corresponding amino acid for '" + codon + "'");
        }
    }

    public static String longestCommonSubstring(String S1, String S2) {
        // Taken from http://en.wikibooks.org/wiki/Algorithm_Implementation/Strings/Longest_common_substring#Java

        int Start = 0;
        int Max = 0;

        for (int i = 0; i < S1.length(); i++) {
            for (int j = 0; j < S2.length(); j++) {
                int x = 0;
                while (S1.charAt(i + x) == S2.charAt(j + x)) {
                    x++;
                    if (((i + x) >= S1.length()) || ((j + x) >= S2.length())) break;
                }

                if (x > Max) {
                    Max = x;
                    Start = i;
                }
            }
        }

        return S1.substring(Start, (Start + Max));
    }

    public static String extractGeneSequence(GFF3Record gene, IndexedFastaSequenceFile ref) {
        return new String(ref.getSubsequenceAt(gene.getSeqid(), gene.getStart(), gene.getEnd()).getBases());
    }

    public static String extractCodingSequence(Collection<GFF3Record> exons, IndexedFastaSequenceFile ref) {
        List<GFF3Record> exonList = new ArrayList<GFF3Record>(exons);
        Collections.sort(exonList, new Comparator<GFF3Record>() {
            public int compare(GFF3Record c1, GFF3Record c2) {
                return c1.getStart() - c2.getStart();
            }
        });

        StringBuilder cdsb = new StringBuilder();

        GFF3Record.Strand strand = GFF3Record.Strand.POSITIVE;
        for (GFF3Record exon : exonList) {
            String exonseq = new String(ref.getSubsequenceAt(exon.getSeqid(), exon.getStart(), exon.getEnd()).getBases());

            cdsb.append(exonseq);

            strand = exon.getStrand();
        }

        String cds = cdsb.toString();
        if (strand.equals(GFF3Record.Strand.NEGATIVE)) {
            cds = SequenceUtils.reverseComplement(cds);
        }

        return cds;
    }

    public static String translateCodingSequence(String cds) {
        StringBuilder tr = new StringBuilder();

        for (int pos = 0, aaPos = 0; pos < cds.length() - 3; pos += 3, aaPos++) {
            String codon = cds.substring(pos, pos + 3);
            String aa = codonToAminoAcid(codon);

            if (!"*".equals(aa) || pos < cds.length() - 3) {
                tr.append(aa);
            }
        }

        return tr.toString();
    }

    public static float fractionGC(String seq) {
        int gcBases = 0;

        for (int i = 0; i < seq.length(); i++) {
            if (seq.charAt(i) == 'G' || seq.charAt(i) == 'g' || seq.charAt(i) == 'C' || seq.charAt(i) == 'c') {
                gcBases++;
            }
        }

        return (float) gcBases / (float) seq.length();
    }

    public static Set<String> kmerizeSequence(String seq, int kmerSize) {
        Set<String> kmers = new HashSet<String>();
        for (int i = 0; i <= seq.length() - kmerSize; i++) {
            String kmer = seq.substring(i, i + kmerSize);

            kmers.add(kmer);
        }

        return kmers;
    }
}
