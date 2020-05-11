version 1.0

import "tasks/Find.wdl" as Find
import "tasks/Utils.wdl" as Utils
import "tasks/PreprocessReads.wdl" as PR
import "tasks/AlignReads.wdl" as AR
import "tasks/AssembleReads.wdl" as ASM
import "tasks/AlignedMetrics.wdl" as AM
import "tasks/AssemblyMetrics.wdl" as ASMM
import "tasks/Finalize.wdl" as FF

workflow AlignValidationProgeny {
    input {
        String sample_name
        File end1
        File end2
        File asm
        Array[Object] refs
        String gcs_out_root_dir
    }

    String outdir = sub(gcs_out_root_dir, "/$", "")

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
                RG             = "@RG\\tID:~{ref_name}\\tSM:~{sample_name}\\tPL:ILLUMINA",
                prefix         = sample_name + "." + ref_name,
                cpus           = 8
        }

        call FF.FinalizeToDir as FinalizeBam {
            input:
                files  = [ BwaMem.aligned_bam, BwaMem.aligned_bai ],
                outdir = outdir + "/aligned/"
        }

#        call AR.Minimap2 {
#            input:
#                reads      = [ asm ],
#                ref_fasta  = entry['fa'],
#                RG         = "@RG\\tID:~{ref_name}\\tSM:~{sample_name}\\tPL:PACBIO",
#                map_preset = "map-pb",
#                prefix     = ref_name
#        }

        ##########
        # Finalize
        ##########

#        call FF.FinalizeToDir as FinalizeContigs {
#            input:
#                files  = [ Assemble.contigs_without_links, Assemble.contigs_with_se_links, Assemble.contigs_with_se_and_pe_links, Assemble.contigs_with_all_links ],
#                outdir = outdir + "/" + ID + "/wo_preprocessing/contigs"
#        }
#
#        call FF.FinalizeToDir as FinalizeModContigs {
#            input:
#                files  = [ AssembleMod.contigs_without_links, AssembleMod.contigs_with_se_links, AssembleMod.contigs_with_se_and_pe_links, AssembleMod.contigs_with_all_links ],
#                outdir = outdir + "/" + ID + "/w_preprocessing/contigs"
#        }
    }
}
