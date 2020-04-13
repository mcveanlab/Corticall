version 1.0

import "Structs.wdl"

task MergedPairedEndReads {
    input {
        File end1
        File end2

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 10

    String prefix = basename(end1, ".end1.fq.gz")

    command <<<
        set -euxo pipefail

        flash ~{end1} ~{end2} -o ~{prefix}

        gzip ~{prefix}.extendedFrags.fastq
        gzip ~{prefix}.notCombined_1.fastq
        gzip ~{prefix}.notCombined_2.fastq
    >>>

    output {
        File extendedFrags = "~{prefix}.extendedFrags.fastq.gz"
        File notCombined_1 = "~{prefix}.notCombined_1.fastq.gz"
        File notCombined_2 = "~{prefix}.notCombined_2.fastq.gz"
        File hist = "~{prefix}.hist"
        File histogram = "~{prefix}.histogram"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             8,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "quay.io/corticall/utils:0.1.0"
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

task TrimFq {
    input {
        File fq
        String prefix

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 4*ceil(size(fq, "GB"))

    command <<<
        set -euxo pipefail

        seqtk trimfq ~{fq} | gzip > ~{prefix}.fq.gz
    >>>

    output {
        File trimmed_fq = "~{prefix}.fq.gz"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             4,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "quay.io/corticall/utils:0.1.0"
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

task Bfc {
    input {
        File end1
        File end2
        String prefix

        Int num_cpus = 4

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 4*ceil(size(end1, "GB") + size(end2, "GB"))

    command <<<
        set -euxo pipefail

        seqtk mergepe ~{end1} ~{end2} | gzip > merged.fq.gz
        bfc -s 23m -k47 -t~{num_cpus} merged.fq.gz | awk '{ print $1 }' | paste - - - - - - - - | tee >(cut -f 1-4 | tr '\t' '\n' | gzip > ~{prefix}.end1.fq.gz) | cut -f 5-8 | tr '\t' '\n' | gzip > ~{prefix}.end2.fq.gz
    >>>

    output {
        File bfc_end1 = "~{prefix}.end1.fq.gz"
        File bfc_end2 = "~{prefix}.end2.fq.gz"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          num_cpus,
        mem_gb:             8,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "quay.io/corticall/utils:0.1.0"
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

task Downsample {
    input {
        File end1
        File end2
        Int num_records
        Int read_length
        Int target_cov
        String prefix

        Int seed = 0

        RuntimeAttr? runtime_attr_override
    }

    Int est_records_to_sample = ceil(target_cov*23000000.0/(2.0*read_length))
    Int final_records_to_sample = if num_records <= est_records_to_sample then 0 else est_records_to_sample

    Int disk_size = 4*ceil(size(end1, "GB") + size(end2, "GB"))

    command <<<
        set -euxo pipefail

        if [ "~{final_records_to_sample}" == "0" ]
        then
            mv ~{end1} ~{prefix}.end1.fq.gz
            mv ~{end2} ~{prefix}.end2.fq.gz
        else
            seqtk sample -s ~{seed} ~{end1} ~{final_records_to_sample} | gzip > ~{prefix}.end1.fq.gz
            seqtk sample -s ~{seed} ~{end2} ~{final_records_to_sample} | gzip > ~{prefix}.end2.fq.gz
        fi
    >>>

    output {
        File ds_end1 = "~{prefix}.end1.fq.gz"
        File ds_end2 = "~{prefix}.end2.fq.gz"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             8,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "quay.io/corticall/utils:0.1.0"
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

