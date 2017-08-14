package uk.ac.ox.well.cortexjdk.utils.file;

import java.io.File;

public class FileAndPathUtils {
    private FileAndPathUtils() {}

    public enum FileType {
        CTX, FASTA, FASTQ, FASTQGZ, BAM, UNKNOWN;
    }

    public static FileType getFileType(File file) {
        return getFileType(file.getName());
    }

    public static FileType getFileType(String filename) {
        if (filename.endsWith(".ctx")) {
            return FileType.CTX;
        } else if (filename.endsWith(".fasta") || filename.endsWith(".fa")) {
            return FileType.FASTA;
        } else if (filename.endsWith(".fastq") || filename.endsWith(".fq")) {
            return FileType.FASTQ;
        } else if (filename.endsWith(".fastq.gz") || filename.endsWith(".fq.gz")) {
            return FileType.FASTQGZ;
        } else if (filename.endsWith(".bam")) {
            return FileType.BAM;
        }

        return FileType.UNKNOWN;
    }
}
