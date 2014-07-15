package uk.ac.ox.well.indiana.attic.analyses.ContigErrors;

import com.google.common.base.Joiner;
import htsjdk.samtools.fastq.FastqReader;
import htsjdk.samtools.fastq.FastqRecord;
import uk.ac.ox.well.indiana.commands.Module;
import uk.ac.ox.well.indiana.utils.arguments.Argument;
import uk.ac.ox.well.indiana.utils.arguments.Output;

import java.io.File;
import java.io.PrintStream;
import java.util.ArrayList;

public class FastqToFastaAndQual extends Module {
    @Argument(fullName="fastq", shortName="fq", doc="Fastq file")
    public ArrayList<File> FASTQS;

    @Output(fullName="fasta", shortName="fo", doc="Fasta output file")
    public PrintStream fo;

    @Output(fullName="qual", shortName="qo", doc="Qual output file")
    public PrintStream qo;

    @Override
    public void execute() {
        for (File fastq : FASTQS) {
            FastqReader fqr = new FastqReader(fastq);

            for (FastqRecord fr : fqr) {
                String[] header = fr.getReadHeader().split("\\s+");
                String fLine = fr.getReadString();
                String qLine = fr.getBaseQualityString();

                Integer[] quals = new Integer[qLine.length()];
                for (int i = 0; i < qLine.length(); i++) {
                    quals[i] = ((int) qLine.charAt(i)) - 33;
                }

                fo.println(">" + header[1]);
                fo.println(fLine);

                qo.println(">" + header[1]);
                qo.println(Joiner.on(" ").join(quals));
            }
        }
    }
}
