version 1.0

import "tasks/Find.wdl" as Find
import "tasks/Utils.wdl" as Utils
import "tasks/PreprocessReads.wdl" as PR
import "tasks/AlignReads.wdl" as AR
import "tasks/AssembleReads.wdl" as ASM
import "tasks/AlignedMetrics.wdl" as AM
import "tasks/Finalize.wdl" as FF

workflow ProcessPfCross {
    input {
        String gcs_input_dir

        Array[Map[String, File]] parents

        File ref_dict
        File ref_fasta
        File ref_fasta_amb
        File ref_fasta_ann
        File ref_fasta_bwt
        File ref_fasta_fai
        File ref_fasta_pac
        File ref_fasta_sa

        File tandem_repeat_bed
        File ref_flat
        File dbsnp_vcf
        File dbsnp_tbi

        String mt_chr_name
        File metrics_locus

        String gcs_out_root_dir
    }

    String outdir = sub(gcs_out_root_dir, "/$", "")

    scatter (entry in parents) {
        call ASM.BuildFromRef { input: ref = entry['ref'], sample_name = entry['sample_name'] }
    }

    call Find.FindFastqs { input: gcs_input_dir = gcs_input_dir }

    scatter (i in range(length(FindFastqs.end1))) {
        String cross = FindFastqs.cross
        File end1 = FindFastqs.end1[i]
        File end2 = FindFastqs.end2[i]

        String ID  = basename(end1, ".end1.fq.gz")
        String SM  = basename(end1, ".end1.fq.gz")
        String PL  = "ILLUMINA"
        String RG = "@RG\\tID:~{ID}\\tSM:~{SM}\\tPL:~{PL}"

        call Utils.CountFastqRecords as CountEnd1 { input: fastq = end1 }
        call Utils.CountFastqRecords as CountEnd2 { input: fastq = end2 }
        call Utils.ReadLength { input: fq = end1 }

        call AR.BwaMem {
            input:
                end1           = end1,
                end2           = end2,
                ref_dict       = ref_dict,
                ref_fasta      = ref_fasta,
                ref_fasta_amb  = ref_fasta_amb,
                ref_fasta_ann  = ref_fasta_ann,
                ref_fasta_bwt  = ref_fasta_bwt,
                ref_fasta_fai  = ref_fasta_fai,
                ref_fasta_pac  = ref_fasta_pac,
                ref_fasta_sa   = ref_fasta_sa,
                RG             = RG,
                prefix         = ID,
                cpus           = 8
        }

        call AM.AlignedMetrics as FullAlignedMetrics {
            input:
                aligned_bam    = BwaMem.aligned_bam,
                aligned_bai    = BwaMem.aligned_bai,
                ref_fasta      = ref_fasta,
                ref_dict       = ref_dict,
                per            = "sample",
                type           = "full",
                label          = ID,
                gcs_output_dir = outdir + "/" + ID
        }

        call ASM.Assemble as Assemble {
            input:
                end1        = end1,
                end2        = end2,
                sample_name = SM,
                parents     = parents,
                ref_ctxs    = BuildFromRef.ref_ctx
        }

        call PR.Bfc { input: end1 = end1, end2 = end2, prefix = "~{ID}.bfc" }

        call PR.Downsample {
            input:
                end1        = Bfc.bfc_end1,
                end2        = Bfc.bfc_end2,
                num_records = CountEnd1.num_records,
                read_length = ReadLength.read_length,
                target_cov  = 100,
                prefix      = "~{ID}.ds"
        }

        call PR.MergedPairedEndReads { input: end1 = Downsample.ds_end1, end2 = Downsample.ds_end2 }

        call PR.TrimFq as TrimEnd1   { input: fq = MergedPairedEndReads.notCombined_1, prefix = "~{ID}.trimmed.end1"          }
        call PR.TrimFq as TrimEnd2   { input: fq = MergedPairedEndReads.notCombined_2, prefix = "~{ID}.trimmed.end2"          }
        call PR.TrimFq as TrimMerged { input: fq = MergedPairedEndReads.extendedFrags, prefix = "~{ID}.trimmed.extendedFrags" }

        call Utils.CountFastqRecords as CountTrimEnd1   { input: fastq = TrimEnd1.trimmed_fq   }
        call Utils.CountFastqRecords as CountTrimEnd2   { input: fastq = TrimEnd2.trimmed_fq   }
        call Utils.CountFastqRecords as CountTrimMerged { input: fastq = TrimMerged.trimmed_fq }

        call AR.BwaMem as BwaMemTrimmed {
            input:
                end1           = TrimEnd1.trimmed_fq,
                end2           = TrimEnd2.trimmed_fq,
                ref_dict       = ref_dict,
                ref_fasta      = ref_fasta,
                ref_fasta_amb  = ref_fasta_amb,
                ref_fasta_ann  = ref_fasta_ann,
                ref_fasta_bwt  = ref_fasta_bwt,
                ref_fasta_fai  = ref_fasta_fai,
                ref_fasta_pac  = ref_fasta_pac,
                ref_fasta_sa   = ref_fasta_sa,
                RG             = RG,
                prefix         = ID,
                cpus           = 8
        }

        call AR.BwaMemSingle as BwaMemMerged {
            input:
                reads          = TrimMerged.trimmed_fq,
                ref_dict       = ref_dict,
                ref_fasta      = ref_fasta,
                ref_fasta_amb  = ref_fasta_amb,
                ref_fasta_ann  = ref_fasta_ann,
                ref_fasta_bwt  = ref_fasta_bwt,
                ref_fasta_fai  = ref_fasta_fai,
                ref_fasta_pac  = ref_fasta_pac,
                ref_fasta_sa   = ref_fasta_sa,
                RG             = RG,
                prefix         = ID,
                cpus           = 8
        }

        call AR.MergeBams { input: bams = [  BwaMemTrimmed.aligned_bam, BwaMemMerged.aligned_bam ] }

        call AM.AlignedMetrics as ModAlignedMetrics {
            input:
                aligned_bam    = MergeBams.merged_bam,
                aligned_bai    = MergeBams.merged_bai,
                ref_fasta      = ref_fasta,
                ref_dict       = ref_dict,
                per            = "sample",
                type           = "mod",
                label          = ID,
                gcs_output_dir = outdir + "/" + ID
        }

        call ASM.Assemble as AssembleMod {
            input:
                end1        = TrimEnd1.trimmed_fq,
                end2        = TrimEnd2.trimmed_fq,
                unpaired    = TrimMerged.trimmed_fq,
                sample_name = SM,
                parents     = parents,
                ref_ctxs    = BuildFromRef.ref_ctx
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
