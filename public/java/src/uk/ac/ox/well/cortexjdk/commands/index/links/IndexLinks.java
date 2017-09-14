package uk.ac.ox.well.cortexjdk.commands.index.links;

import htsjdk.samtools.util.BlockCompressedOutputStream;
import org.apache.commons.math3.util.Pair;
import uk.ac.ox.well.cortexjdk.commands.Module;
import uk.ac.ox.well.cortexjdk.utils.arguments.Argument;
import uk.ac.ox.well.cortexjdk.utils.exceptions.CortexJDKException;
import uk.ac.ox.well.cortexjdk.utils.io.cortex.graph.CortexBinaryKmer;
import uk.ac.ox.well.cortexjdk.utils.io.cortex.graph.CortexByteKmer;
import uk.ac.ox.well.cortexjdk.utils.io.cortex.graph.CortexRecord;
import uk.ac.ox.well.cortexjdk.utils.io.cortex.links.CortexLinksIterable;
import uk.ac.ox.well.cortexjdk.utils.io.cortex.links.CortexLinksRecord;
import uk.ac.ox.well.cortexjdk.utils.progress.ProgressMeter;
import uk.ac.ox.well.cortexjdk.utils.progress.ProgressMeterFactory;

import java.io.IOException;
import java.nio.ByteBuffer;
import java.nio.channels.FileChannel;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.*;

import static java.nio.file.StandardOpenOption.CREATE_NEW;
import static java.nio.file.StandardOpenOption.READ;
import static java.nio.file.StandardOpenOption.WRITE;

/**
 * Created by kiran on 13/09/2017.
 */
public class IndexLinks extends Module {
    @Argument(fullName="links", shortName="l", doc="Links")
    public CortexLinksIterable LINKS;

    @Argument(fullName="source", shortName="s", doc="Link source")
    public String SOURCE;

    @Override
    public void execute() {
        Path bgzipPath = Paths.get(LINKS.getCortexLinksFile().getAbsolutePath().replace(".ctp.gz", ".ctp.bgz"));
        Path indexPath = Paths.get(LINKS.getCortexLinksFile().getAbsolutePath().replace(".ctp.gz", ".ctp.bgz.idx"));

        BlockCompressedOutputStream bc = new BlockCompressedOutputStream(bgzipPath.toFile().getAbsolutePath(), 9);

        log.info("Writing bgzipped links and link index to:");
        log.info("  - {}", bgzipPath);
        log.info("  - {}", indexPath);

        try (FileChannel fc = FileChannel.open(indexPath, CREATE_NEW, READ, WRITE)) {
            storeHeader(fc, LINKS, SOURCE);
            storeRecordsAndIndex(fc, bc, LINKS);

            fc.close();
            bc.close();
        } catch (IOException e) {
            throw new CortexJDKException("IOException", e);
        }

        /*
        log.info("Verifying index");

        try {
            BlockCompressedInputStream bi = new BlockCompressedInputStream(bgzipPath.toFile().getAbsoluteFile());

            BinaryFile bf = new BinaryFile(indexPath.toFile(), "r");

            byte[] magicWordStart = new byte[6];
            bf.read(magicWordStart);

            int numColors = bf.readInt();
            int kmerSize = bf.readInt();
            long numKmersInGraph = bf.readLong();
            long numKmersWithLinks = bf.readLong();
            long numLinkBytes = bf.readLong();

            byte[] source = new byte[bf.readInt()];
            bf.read(source);

            for (int c = 0; c < numColors; c++) {
                byte[] sn = new byte[bf.readInt()];
                bf.read(sn);
            }

            byte[] magicWordEnd = new byte[6];
            bf.read(magicWordEnd);

            log.info("{} {} {}", magicWordStart, magicWordEnd, Arrays.equals(magicWordStart, magicWordEnd));
            log.info("  kmerSize: {}", kmerSize);
            log.info("  numKmersInGraph: {}", numKmersInGraph);
            log.info("  numKmersWithLinks: {}", numKmersWithLinks);
            log.info("  numLinkBytes: {}", numLinkBytes);

            int numKmerBits = CortexRecord.getKmerBits(kmerSize);
            for (int i = 0; i < numKmersWithLinks; i++) {
                long[] kmer = new long[numKmerBits];
                for (int j = 0; j < numKmerBits; j++) {
                    kmer[j] = bf.readLong();
                }

                long pos = bf.readLong();
                int length = bf.readInt();

                String sk = new String(CortexRecord.decodeBinaryKmer(kmer, kmerSize, numKmerBits));

                bi.seek(pos);
                byte[] recbuf = new byte[length];
                bi.read(recbuf);
                CortexLinksRecord clr = new CortexLinksRecord(recbuf);

                log.info("{} {} {} {} {} {}", numKmersWithLinks, i, sk, pos, length, clr);
            }
        } catch (IOException e) {
            throw new CortexJDKException("IOException", e);
        }
        */
    }

    private void storeHeader(FileChannel fc, CortexLinksIterable links, String source) throws IOException {
        int capacity = 6 + 4 + 4 + 8 + 8 + 8;
        capacity += 4 + source.length();
        for (int c = 0; c < links.getNumColors(); c++) {
            capacity += 4 + links.getColor(c).getSampleName().length();
        }
        capacity += 6;

        ByteBuffer bb = ByteBuffer.allocateDirect(capacity);

        bb.put("LNKIDX".getBytes());
        bb.putInt(links.getNumColors());
        bb.putInt(links.getKmerSize());
        bb.putLong(links.getNumKmersInGraph());
        bb.putLong(links.getNumKmersWithLinks());
        bb.putLong(links.getLinkBytes());

        bb.putInt(source.length());
        bb.put(source.getBytes());

        for (int c = 0; c < links.getNumColors(); c++) {
            String sn = links.getColor(c).getSampleName();
            bb.putInt(sn.length());
            bb.put(sn.getBytes());
        }

        bb.put("LNKIDX".getBytes());

        bb.flip();
        fc.write(bb);
    }

    private void storeRecordsAndIndex(FileChannel fc, BlockCompressedOutputStream bc, CortexLinksIterable links) throws IOException {
        int numKmerBits = CortexRecord.getKmerBits(links.getKmerSize());
        int bufferSize = (int) (links.getNumKmersWithLinks() * ((8*numKmerBits) + 8 + 4));

        log.info("Allocating kmer table ({} records, {} bytes)", links.getNumKmersWithLinks(), bufferSize);

        bc.write(links.getJSONHeader().getBytes());
        bc.write("\n".getBytes());
        bc.write(links.getComments().getBytes());
        bc.write("\n".getBytes());

        ProgressMeter pm = new ProgressMeterFactory()
                .header("Processing links")
                .message("links")
                .maxRecord(links.getNumKmersWithLinks())
                .make(log);

        Map<CortexByteKmer, Pair<Long, Integer>> kmerTable = new TreeMap<>();

        for (CortexLinksRecord clr : links) {
            CortexByteKmer cbk = clr.getKmerAsByteKmer();
            kmerTable.put(cbk, new Pair<>(bc.getPosition(), clr.toString().length()));

            bc.write(clr.toString().getBytes());

            pm.update();
        }

        ByteBuffer bb = ByteBuffer.allocateDirect(bufferSize);
        for (CortexByteKmer cb : kmerTable.keySet()) {
            CortexBinaryKmer cbk = new CortexBinaryKmer(cb.getKmer());
            for (long l : cbk.getBinaryKmer()) {
                bb.putLong(l);
            }
            bb.putLong(kmerTable.get(cb).getFirst());
            bb.putInt(kmerTable.get(cb).getSecond());
        }

        bb.flip();
        fc.write(bb);
    }
}
