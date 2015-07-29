package uk.ac.ox.well.indiana.utils.alignment.exact;

import com.javacodegeeks.stringsearch.BM;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.samtools.reference.ReferenceSequenceFile;
import htsjdk.samtools.util.Interval;
import uk.ac.ox.well.indiana.utils.sequence.SequenceUtils;

import java.util.HashMap;
import java.util.List;
import java.util.Map;

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
        String rc = SequenceUtils.reverseComplement(fw);

        String finalName = null;
        int finalPos = -1;
        int numHomes = 0;
        boolean isNegative = false;

        if (!fw.isEmpty() && !rc.isEmpty()) {
            for (String name : ref.keySet()) {
                String seq = ref.get(name);
                List<Integer> posFw = BM.findAll(fw, seq);
                List<Integer> posRc = BM.findAll(rc, seq);

                numHomes += posFw.size() + posRc.size();

                if (posFw.size() == 1) {
                    finalName = name;
                    finalPos = posFw.iterator().next();
                } else if (posRc.size() == 1) {
                    finalName = name;
                    finalPos = posRc.iterator().next();
                    isNegative = true;
                }

                if (numHomes > 1) {
                    break;
                }
            }
        }

        return (finalName != null && numHomes == 1) ? new Interval(finalName, finalPos, finalPos + fw.length(), isNegative, "unknown") : null;
    }
}
