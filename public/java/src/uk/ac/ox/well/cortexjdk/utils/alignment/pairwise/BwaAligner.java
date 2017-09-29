package uk.ac.ox.well.cortexjdk.utils.alignment.pairwise;

import com.github.lindenb.jbwa.jni.AlnRgn;
import com.github.lindenb.jbwa.jni.BwaIndex;
import com.github.lindenb.jbwa.jni.BwaMem;
import com.github.lindenb.jbwa.jni.ShortRead;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.reference.FastaSequenceFile;
import uk.ac.ox.well.cortexjdk.utils.exceptions.CortexJDKException;
import uk.ac.ox.well.cortexjdk.utils.packageutils.InternalLibraryResource;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

public class BwaAligner {
    private static final InternalLibraryResource bwajni = new InternalLibraryResource("/libbwajni.jnilib");
    private static final InternalLibraryResource bwaso = new InternalLibraryResource("/libbwajni.so");

    private final BwaIndex index;
    private final BwaMem mem;
    private final SAMFileHeader header;

    public BwaAligner(String ref) {
        System.loadLibrary("bwajni");

        try {
            index = new BwaIndex(new File(ref));

            FastaSequenceFile fa = new FastaSequenceFile(new File(ref), true);
            header = new SAMFileHeader();
            header.setSequenceDictionary(fa.getSequenceDictionary());

            mem = new BwaMem(index);
        } catch (IOException e) {
            throw new CortexJDKException("Could not initialize bwajni library");
        }
    }

    public List<SAMRecord> align(String query) {
        ShortRead read = new ShortRead("unknown", query.getBytes(), new byte[0]);

        List<SAMRecord> alignments = new ArrayList<>();
        try {
            for (AlnRgn a : mem.align(read)) {
                SAMRecord rec = new SAMRecord(header);

                rec.setReadName(read.getName());
                rec.setReadBases(read.getBases());
                rec.setBaseQualities(read.getQualities());
                rec.setReadNegativeStrandFlag(a.getStrand() == '-');
                rec.setReferenceName(a.getChrom());
                rec.setAlignmentStart((int) a.getPos());
                rec.setMappingQuality(a.getMQual());
                rec.setCigarString(a.getCigar());
                rec.setAttribute("NM", a.getNm());
                rec.setSupplementaryAlignmentFlag(a.getSecondary() >= 0);

                alignments.add(rec);
            }
        } catch (IOException e) {
            throw new CortexJDKException("Failed when aligning '" + query + "'");
        }

        return alignments;
    }

    public List<SAMRecord> align(List<String> queries) {
        List<SAMRecord> alignments = new ArrayList<>();

        queries.forEach(q -> alignments.addAll(align(q)));

        return alignments;
    }

    public void close() {
        index.close();
        mem.dispose();
    }
}
