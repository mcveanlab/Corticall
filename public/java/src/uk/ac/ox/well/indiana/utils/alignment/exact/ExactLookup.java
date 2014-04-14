package uk.ac.ox.well.indiana.utils.alignment.exact;

import com.javacodegeeks.stringsearch.BM;
import net.sf.picard.reference.ReferenceSequence;
import net.sf.picard.reference.ReferenceSequenceFile;
import net.sf.picard.util.Interval;

import java.util.*;

public class ExactLookup {
    private Map<String, String> ref = new HashMap<String, String>();

    public ExactLookup(ReferenceSequenceFile sequences) {
        ReferenceSequence rseq;
        while ((rseq = sequences.nextSequence()) != null) {
            String name = rseq.getName().split("\\s+")[0];
            String seq = new String(rseq.getBases());

            ref.put(name, seq);
        }
    }

    public Interval find(String s) {
        String fw = s.replaceAll("-", "");

        String finalName = null;
        int finalPos = -1;
        int numHomes = 0;

        for (String name : ref.keySet()) {
            String seq = ref.get(name);
            List<Integer> pos = BM.findAll(fw, seq);

            numHomes += pos.size();

            if (pos.size() == 1) {
                finalName = name;
                finalPos = pos.iterator().next();
            }
        }

        return (finalName != null && numHomes == 1) ? new Interval(finalName, finalPos, finalPos + fw.length()) : null;
    }
}
