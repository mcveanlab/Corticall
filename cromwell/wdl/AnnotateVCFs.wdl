version 1.0

import "tasks/Find.wdl" as Find
import "tasks/Utils.wdl" as Utils
import "tasks/PreprocessReads.wdl" as PR
import "tasks/AlignReads.wdl" as AR
import "tasks/AssembleReads.wdl" as ASM
import "tasks/AlignedMetrics.wdl" as AM
import "tasks/AssemblyMetrics.wdl" as ASMM
import "tasks/Finalize.wdl" as FF

workflow AnnotateVCFs {
    input {
        String gcs_input_dir
        String cross

        String gcs_out_root_dir
    }

    String indir = sub(gcs_input_dir, "/$", "")
    String outdir = sub(gcs_out_root_dir, "/$", "")

    call FindVCFs { input: indir = indir }

    scatter (vcf in read_lines(FindVCFs.vcf_list)) {
        call AnnotateVCF { input: vcf = vcf, cross = cross }

        String bn = basename(vcf, ".calls.ann.vcf")

        call FF.FinalizeToDir as FinalizeCalls {
            input:
                files  = [ AnnotateVCF.snpeff_vcf, AnnotateVCF.genes, AnnotateVCF.summary ],
                outdir = outdir + "/" + cross + "/" + bn
        }
    }
}

task FindVCFs {
    input {
        String indir

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 1

    command <<<
        set -euxo pipefail

        gsutil ls ~{indir}/**ann.vcf > vcf_list.txt
    >>>

    output {
        File vcf_list = "vcf_list.txt"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             2,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        0,
        docker:             "quay.io/corticall/utils:0.1.1"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task AnnotateVCF {
    input {
        File vcf
        String cross

        RuntimeAttr? runtime_attr_override
    }

    String ann = basename(vcf, "vcf") + "snpeff.vcf"

    Int disk_size = 1 + 4*ceil(size(vcf, "GB"))

    command <<<
        set -euxo pipefail

        java -jar /snpeff/snpEff.jar pf.~{cross} ~{vcf} > ~{ann}
    >>>

    output {
        File snpeff_vcf = ann
        File genes = "snpEff_genes.txt"
        File summary = "snpEff_summary.html"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             2,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  1,
        max_retries:        0,
        docker:             "quay.io/corticall/snpeff:0.1.0"
    }
    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}
