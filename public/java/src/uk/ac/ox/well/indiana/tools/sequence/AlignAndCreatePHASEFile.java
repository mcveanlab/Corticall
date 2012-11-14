package uk.ac.ox.well.indiana.tools.sequence;

import net.sf.picard.reference.FastaSequenceFile;
import net.sf.picard.reference.ReferenceSequence;
import uk.ac.ox.well.indiana.tools.Tool;
import uk.ac.ox.well.indiana.utils.alignment.NeedlemanWunsch;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.TreeMap;

public class AlignAndCreatePHASEFile extends Tool {
    @Argument(fullName="fasta", shortName="f", doc="Fasta file to align and turn into a PHASE file")
    public FastaSequenceFile FASTA;

    @Argument(fullName="withRespectTo", shortName="wrt", doc="Contig to which sequences should be aligned")
    public String WRT;

    @Argument(fullName="contigsToProcess", shortName="c", doc="Contigs to align (should all be from the same segment")
    public HashSet<String> CONTIGS;

    @Output
    public PrintStream out;

    private HashMap<String, String> loadSequences(FastaSequenceFile fasta, String wrt, HashSet<String> contigs) {
        HashMap<String, String> sequences = new HashMap<String, String>();

        ReferenceSequence seq;
        while ((seq = fasta.nextSequence()) != null) {
            String[] names = seq.getName().split("\\|");
            String name = names[3];

            if (wrt.equals(name) || contigs.contains(name)) {
                sequences.put(name, new String(seq.getBases()));
            }
        }

        return sequences;
    }

    private void addToMutationMap(TreeMap<Integer, HashMap<String, String>> mutationMap, int position, String sample, String mutation) {
        if (!mutationMap.containsKey(position)) {
            mutationMap.put(position, new HashMap<String, String>());
        }

        mutationMap.get(position).put(sample, mutation);
    }

    @Override
    public int execute() {
        HashMap<String, String> sequences = loadSequences(FASTA, WRT, CONTIGS);

        log.info("Sequences loaded: {}", sequences.size());

        TreeMap<Integer, HashMap<String, String>> mutationMap = new TreeMap<Integer, HashMap<String, String>>();

        String wrtSequence = sequences.get(WRT);
        sequences.remove(WRT);

        int seqIndex = 0;
        for (String name : sequences.keySet()) {
            if (seqIndex % (sequences.size() / 100) == 0) {
                log.info("Processing sequence {}/{}", seqIndex, sequences.size());
            }
            seqIndex++;

            String sequence = sequences.get(name);

            NeedlemanWunsch nw = new NeedlemanWunsch(wrtSequence, sequence);
            String[] alignments = nw.getAlignment();

            int pos = 1;
            int i = 0;
            while (i < alignments[0].length()) {
                if (alignments[1].charAt(i) == '-') { // deletion
                    int length = 1;
                    for (int j = i + 1; j < alignments[1].length(); j++) {
                        if (alignments[1].charAt(j) == '-') { length++; }
                        else { break; }
                    }

                    String deletedBases = alignments[0].substring(i, i + length);

                    //mutations.add(pos + "D" + length + deletedBases);
                    addToMutationMap(mutationMap, pos, name, pos + "D" + length + deletedBases);

                    i += length;
                    pos += length;
                } else if (alignments[0].charAt(i) == '-') { // insertion
                    int length = 1;
                    for (int j = i + 1; j < alignments[0].length(); j++) {
                        if (alignments[0].charAt(j) == '-') { length++; }
                        else { break; }
                    }

                    String insertedBases = alignments[1].substring(i, i + length);

                    //mutations.add((pos-1) + "I" + length + insertedBases);
                    addToMutationMap(mutationMap, pos, name, (pos-1) + "I" + length + insertedBases);

                    i += length;
                } else if (alignments[0].charAt(i) != alignments[1].charAt(i)) { // mismatch
                    //mutations.add(pos + "M" + alignments[0].charAt(i) + ">" + alignments[1].charAt(i));
                    addToMutationMap(mutationMap, pos, name, pos + "M" + alignments[0].charAt(i) + ">" + alignments[1].charAt(i));

                    i++;
                    pos++;
                } else { // match
                    i++;
                    pos++;
                }
            }
        }

        ArrayList<Integer> biallelicLoci = new ArrayList<Integer>();
        for (Integer pos : mutationMap.keySet()) {
            HashSet<String> alleles = new HashSet<String>(mutationMap.get(pos).values());

            if (alleles.size() == 1) { // biallelic mutation
                biallelicLoci.add(pos);
            }
        }

        log.info("Mutations: {}, Biallelic mutations: {}", mutationMap.size(), biallelicLoci.size());

        for (Integer pos : mutationMap.keySet()) {
            out.print("\tloc_" + pos);
        }

        out.println();

        for (String name : sequences.keySet()) {
            out.print(name);

            for (Integer pos : mutationMap.keySet()) {
                String allele = mutationMap.get(pos).get(name);

                if (allele == null) {
                    out.print("\t0");
                } else {
                    out.print("\t" + allele);
                }

            }

            out.println();
        }

        /*
        out.println("0");
        out.println(sequences.size());
        out.println(biallelicLoci.size());

        out.print("P");
        for (Integer pos : biallelicLoci) {
            out.print(" " + pos);
        }
        out.println();

        for (Integer pos : biallelicLoci) {
            out.print("S");
        }
        out.println();

        for (String name : sequences.keySet()) {
            for (Integer pos : biallelicLoci) {
                out.print(mutationMap.get(pos).containsKey(name) ? 1 : 0);
            }
            out.println();
        }
        */

        return 0;
    }
}
