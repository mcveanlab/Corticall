version 1.0

import "Structs.wdl"

task BwaMem {
    input {
        File end1
        File end2

        File ref_dict
        File ref_fasta
        File ref_fasta_amb
        File ref_fasta_ann
        File ref_fasta_bwt
        File ref_fasta_fai
        File ref_fasta_pac
        File ref_fasta_sa

        String RG

        String prefix = "out"
        Int cpus = 4
        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 20*ceil(size(end1, "GB") + size(end2, "GB") + size(ref_dict, "GB") + size(ref_fasta, "GB") + size(ref_fasta_amb, "GB") + size(ref_fasta_ann, "GB") + size(ref_fasta_bwt, "GB") + size(ref_fasta_fai, "GB") + size(ref_fasta_pac, "GB") + size(ref_fasta_sa, "GB"))

    command <<<
        set -euxo pipefail

        bwa mem -R "~{RG}" -t ~{cpus} ~{ref_fasta} ~{end1} ~{end2} | samtools fixmate - - | samtools sort -o ~{prefix}.bam
        samtools index ~{prefix}.bam
    >>>

    output {
        File aligned_bam = "~{prefix}.bam"
        File aligned_bai = "~{prefix}.bam.bai"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          cpus,
        mem_gb:             30,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "quay.io/corticall/align:0.1.0"
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

task BwaMemSingle {
    input {
        File reads

        File ref_dict
        File ref_fasta
        File ref_fasta_amb
        File ref_fasta_ann
        File ref_fasta_bwt
        File ref_fasta_fai
        File ref_fasta_pac
        File ref_fasta_sa

        String RG

        String prefix = "out"
        Int cpus = 4
        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 20*ceil(size(reads, "GB") + size(ref_dict, "GB") + size(ref_fasta, "GB") + size(ref_fasta_amb, "GB") + size(ref_fasta_ann, "GB") + size(ref_fasta_bwt, "GB") + size(ref_fasta_fai, "GB") + size(ref_fasta_pac, "GB") + size(ref_fasta_sa, "GB"))

    command <<<
        set -euxo pipefail

        bwa mem -R "~{RG}" -t ~{cpus} ~{ref_fasta} ~{reads} | samtools sort -o ~{prefix}.bam
        samtools index ~{prefix}.bam
    >>>

    output {
        File aligned_bam = "~{prefix}.bam"
        File aligned_bai = "~{prefix}.bam.bai"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          cpus,
        mem_gb:             30,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "quay.io/corticall/align:0.1.0"
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

task MergeBams {
    input {
        Array[File] bams

        String prefix = "out"
        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 4*ceil(size(bams, "GB"))

    command <<<
        set -euxo pipefail

        samtools merge -p -c -@2 --no-PG ~{prefix}.bam ~{sep=" " bams}
        samtools index ~{prefix}.bam
    >>>

    output {
        File merged_bam = "~{prefix}.bam"
        File merged_bai = "~{prefix}.bam.bai"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             20,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-align:0.1.26"
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

# A simple task to covert SAM-format alignment to PAF format
task SAMtoPAF {
    input {
        File sam_formatted_file
        File? index

        RuntimeAttr? runtime_attr_override
    }

    # currently we only support bam (not sam or cram)
    String prefix = basename(sam_formatted_file, ".bam")

    Int disk_size = 2*ceil(size(sam_formatted_file, "GB"))

    command <<<
        set -eu

        filename=$(basename -- ~{sam_formatted_file})
        extension="${filename##*.}"
        if [[ "$extension" == "sam" ]]; then
            /minimap2-2.17_x64-linux/k8 \
                /minimap2-2.17_x64-linux/paftools.js \
                sam2paf \
                -L \
                ~{sam_formatted_file} \
                > ~{prefix}".paf"
        elif [[ "$extension" == "bam" ]]; then
            samtools view -h ~{sam_formatted_file} | \
                /minimap2-2.17_x64-linux/k8 \
                /minimap2-2.17_x64-linux/paftools.js \
                sam2paf \
                -L \
                - \
                > ~{prefix}".paf"
        else
            echo "Currently we only support SAM or BAM (not CRAM)." && exit 1;
        fi
    >>>

    output {
        File pat_formatted_file = "~{prefix}.paf"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             4,
        disk_gb:            "~{disk_size}",
        boot_disk_gb:       10,
        preemptible_tries:  1,
        max_retries:        0,
        docker:             "us.gcr.io/broad-dsp-lrma/lr-align:0.1.26"
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
