version 1.0

import "tasks/Find.wdl" as Find
import "tasks/Utils.wdl" as Utils
import "tasks/AssembleReads.wdl" as ASM
import "tasks/Finalize.wdl" as FF

workflow PfAsm {
    input {
        String gcs_input_dir

        String gcs_out_root_dir
    }

    String outdir = sub(gcs_out_root_dir, "/$", "")

    call Find.FindFastqs { input: gcs_input_dir = gcs_input_dir }

    scatter (i in range(length(FindFastqs.end1))) {
        String cross = FindFastqs.cross
        File end1 = FindFastqs.end1[i]
        File end2 = FindFastqs.end2[i]
        String sample_name = basename(end1, ".end1.fq.gz")

        call ASM.Assemble {
            input:
                cross       = cross,
                end1        = end1,
                end2        = end2,
                sample_name = sample_name
        }
    }

#    ##########
#    # Finalize
#    ##########
#
#    call FF.FinalizeToDir as FinalizeSVs {
#        input:
#            files = [ CallSVs.pbsv_vcf,       CallSVs.pbsv_tbi,
#                      CallSVs.sniffles_vcf,   CallSVs.sniffles_tbi,
#                      CallSVs.svim_vcf,       CallSVs.svim_tbi ],
#            outdir = outdir + "/" + DIR[0] + "/variants"
#    }
#
#    call FF.FinalizeToDir as FinalizeSmallVariants {
#        input:
#            files = [ CallSmallVariants.longshot_vcf, CallSmallVariants.longshot_tbi ],
#            outdir = outdir + "/" + DIR[0] + "/variants"
#    }
#
#    call FF.FinalizeToDir as FinalizeMergedRuns {
#        input:
#            files = [ MergeRuns.merged_bam, MergeRuns.merged_bai ],
#            outdir = outdir + "/" + DIR[0] + "/alignments"
#    }
}
