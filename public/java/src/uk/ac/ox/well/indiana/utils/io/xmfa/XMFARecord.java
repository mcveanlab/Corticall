package uk.ac.ox.well.indiana.utils.io.xmfa;

import com.google.common.base.Joiner;
import htsjdk.samtools.reference.ReferenceSequence;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

public class XMFARecord extends HashMap<String, ReferenceSequence> {
    @Override
    public String toString() {
        List<String> pieces = new ArrayList<String>();

        for (String key : keySet()) {
            ReferenceSequence seq = get(key);

            String[] refnamePieces = seq.getName().split("\\s+");
            String locus = refnamePieces[1];
            String id = new File(seq.getName().split("\\s+")[3]).getName().replace(".fasta", "");

            pieces.add(key + "=XFMARecord{" + locus + " " + id + "}");
        }

        return Joiner.on("; ").join(pieces);
    }
}
