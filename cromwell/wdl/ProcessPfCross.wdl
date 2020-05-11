version 1.0

import "tasks/Find.wdl" as Find
import "tasks/Utils.wdl" as Utils
import "tasks/PreprocessReads.wdl" as PR
import "tasks/AlignReads.wdl" as AR
import "tasks/AssembleReads.wdl" as ASM
import "tasks/AlignedMetrics.wdl" as AM
import "tasks/AssemblyMetrics.wdl" as ASMM
import "tasks/Finalize.wdl" as FF

workflow ProcessPfCross {
    input {
        String gcs_input_dir

        Array[Object] parents

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
        call ASM.BuildFromRef { input: ref = entry['ref'], ref_name = entry['ref_name'] }
    }

    call Find.FindFastqs { input: gcs_input_dir = gcs_input_dir }

    scatter (i in range(length(FindFastqs.end1))) {
        String cross = FindFastqs.cross
        File end1 = FindFastqs.end1[i]
        File end2 = FindFastqs.end2[i]

        String ID = basename(end1, ".end1.fq.gz")
        String SM = basename(end1, ".end1.fq.gz")
        String PL = "ILLUMINA"
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
        call PR.MergedPairedEndReads { input: end1 = Bfc.bfc_end1, end2 = Bfc.bfc_end2 }

        call PR.TrimFq as TrimEnd1   { input: fq = MergedPairedEndReads.notCombined_1, prefix = "~{ID}.trimmed.end1"          }
        call PR.TrimFq as TrimEnd2   { input: fq = MergedPairedEndReads.notCombined_2, prefix = "~{ID}.trimmed.end2"          }
        call PR.TrimFq as TrimMerged { input: fq = MergedPairedEndReads.extendedFrags, prefix = "~{ID}.trimmed.extendedFrags" }

        call ASM.Assemble as AssembleMod {
            input:
                end1        = TrimEnd1.trimmed_fq,
                end2        = TrimEnd2.trimmed_fq,
                unpaired    = TrimMerged.trimmed_fq,
                sample_name = SM,
                parents     = parents,
                ref_ctxs    = BuildFromRef.ref_ctx
        }

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

        call AR.MergeBams { input: bams = [ BwaMemTrimmed.aligned_bam, BwaMemMerged.aligned_bam ] }

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

        scatter (entry in parents) {
            call ASMM.AssemblyMetrics {
                input:
                    ref = entry['ref'],
                    asms = [ Assemble.contigs_without_links, Assemble.contigs_with_se_links, Assemble.contigs_with_se_and_pe_links, Assemble.contigs_with_all_links,
                             AssembleMod.contigs_without_links, AssembleMod.contigs_with_se_links, AssembleMod.contigs_with_se_and_pe_links, AssembleMod.contigs_with_all_links ]
            }

            String ref_name = entry['ref_name']
            call FF.FinalizeToDir as FinalizeQuast {
                input:
                    files = [ AssemblyMetrics.report_txt, AssemblyMetrics.report_html, AssemblyMetrics.report_pdf, AssemblyMetrics.transposed_report_txt ],
                    outdir = outdir + "/" + ID + "/contigs/quast/" + ref_name
            }
        }

        ##########
        # Finalize
        ##########

        call FF.FinalizeToDir as FinalizeContigs {
            input:
                files  = [ Assemble.contigs_without_links, Assemble.contigs_with_se_links, Assemble.contigs_with_se_and_pe_links, Assemble.contigs_with_all_links ],
                outdir = outdir + "/" + ID + "/wo_preprocessing/contigs"
        }

        call FF.FinalizeToDir as FinalizeModContigs {
            input:
                files  = [ AssembleMod.contigs_without_links, AssembleMod.contigs_with_se_links, AssembleMod.contigs_with_se_and_pe_links, AssembleMod.contigs_with_all_links ],
                outdir = outdir + "/" + ID + "/w_preprocessing/contigs"
        }

        call FF.FinalizeToDir as FinalizeGraph {
            input:
                files = flatten([ [Assemble.final_ctx, Assemble.se_ctp, Assemble.pe_ctp], Assemble.ref_ctps ]),
                outdir = outdir + "/" + ID + "/wo_preprocessing/graph"
        }

        call FF.FinalizeToDir as FinalizeModGraph {
            input:
                files = flatten([ [AssembleMod.final_ctx, AssembleMod.se_ctp, AssembleMod.pe_ctp], AssembleMod.ref_ctps ]),
                outdir = outdir + "/" + ID + "/wo_preprocessing/graph"
        }
    }
}
