version 1.0

import "tasks/Find.wdl" as Find
import "tasks/Utils.wdl" as Utils
import "tasks/PreprocessReads.wdl" as PR
import "tasks/AlignReads.wdl" as AR
import "tasks/AssembleReads.wdl" as ASM
import "tasks/AlignedMetrics.wdl" as AM
import "tasks/AssemblyMetrics.wdl" as ASMM
import "tasks/Finalize.wdl" as FF

workflow AlignToParents {
    input {
        String gcs_input_dir
        Array[Object] refs
        String gcs_out_root_dir
    }

    String outdir = sub(gcs_out_root_dir, "/$", "")

    call Find.FindFastqs { input: gcs_input_dir = gcs_input_dir }

    scatter (i in range(length(FindFastqs.end1))) {
        String cross = FindFastqs.cross
        File end1 = FindFastqs.end1[i]
        File end2 = FindFastqs.end2[i]

        String ID = basename(end1, ".end1.fq.gz")
        String SM = basename(end1, ".end1.fq.gz")
        String PL = "ILLUMINA"
        String RG = "@RG\\tID:~{ID}\\tSM:~{SM}\\tPL:~{PL}"

        scatter (entry in refs) {
            String ref_name = entry['ref_name']

            call AR.BwaMem {
                input:
                    end1           = end1,
                    end2           = end2,
                    ref_dict       = entry['dict'],
                    ref_fasta      = entry['fa'],
                    ref_fasta_amb  = entry['amb'],
                    ref_fasta_ann  = entry['ann'],
                    ref_fasta_bwt  = entry['bwt'],
                    ref_fasta_fai  = entry['fai'],
                    ref_fasta_pac  = entry['pac'],
                    ref_fasta_sa   = entry['sa'],
                    RG             = "@RG\\tID:~{ref_name}\\tSM:~{SM}\\tPL:ILLUMINA",
                    prefix         = SM + "." + ref_name,
                    cpus           = 8
            }

            call FF.FinalizeToDir as FinalizeBam {
                input:
                    files  = [ BwaMem.aligned_bam, BwaMem.aligned_bai ],
                    outdir = outdir + "/aligned/"
            }
        }
    }
}
