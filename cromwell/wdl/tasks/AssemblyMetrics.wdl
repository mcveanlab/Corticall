version 1.0

import "Structs.wdl"
import "Finalize.wdl" as FF

workflow AssemblyMetrics {
    input {
        File ref
        Array[File] asms
    }

    call Quast { input: ref = ref, asms = asms }

    output {
        File report_txt = Quast.report_txt
        File report_html = Quast.report_html
        File report_pdf = Quast.report_pdf
        File transposed_report_txt = Quast.transposed_report_txt
    }
}

task Quast {
    input {
        File ref
        Array[File] asms

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 10
    Int cpus = 4

    command <<<
        set -euxo pipefail

        quast -o quast_results -r ~{ref} -m 1 -t ~{cpus} ~{sep=' ' asms}
    >>>

    output {
        File report_txt = "quast_results/report.txt"
        File report_html = "quast_results/report.html"
        File report_pdf = "quast_results/report.pdf"
        File transposed_report_txt = "quast_results/transposed_report.txt"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          cpus,
        mem_gb:             4,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "quay.io/biocontainers/quast:5.0.2--1"
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
