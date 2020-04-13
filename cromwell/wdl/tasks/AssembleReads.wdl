version 1.0

import "Structs.wdl"

workflow Assemble {
    input {
        File end1
        File end2
        File? unpaired
        String sample_name
        Array[Map[String, File]] parents
        Array[File] ref_ctxs
    }

    call Build {
        input:
            end1        = end1,
            end2        = end2,
            unpaired    = unpaired,
            sample_name = sample_name,
    }

    call Clean {
        input:
            raw_ctx     = Build.raw_ctx,
            sample_name = sample_name
    }

    call RecoverExcludedKmers {
        input:

    }

    call InferEdges {
        input:
            clean_ctx   = Clean.clean_ctx,
            sample_name = sample_name
    }

    call ThreadSe {
        input:
            infer_ctx   = InferEdges.infer_ctx,
            sample_name = sample_name,
            ends        = select_all([ end1, end2, unpaired ])
    }

    call ThreadPe {
        input:
            infer_ctx   = InferEdges.infer_ctx,
            sample_name = sample_name,
            end1        = end1,
            end2        = end2,
            unpaired    = unpaired,
            ctps        = [ ThreadSe.ctp ]
    }

    scatter (entry in parents) {
        call ThreadRef {
            input:
                infer_ctx   = InferEdges.infer_ctx,
                sample_name = sample_name,
                ref_name    = entry['sample'],
                ref         = entry['ref']
        }
    }

    call Contigs {
        input:
            infer_ctx = InferEdges.infer_ctx,
            sample_name = sample_name,
            ctps = flatten([[ThreadPe.ctp], ThreadRef.ctp ]),
    }

    output {
        File raw_ctx = Build.raw_ctx
    }
}

task Build {
    input {
        File end1
        File end2
        File? unpaired
        String sample_name

        Int mem_mult = 8

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 10*ceil(size(select_all([end1, end2, unpaired]), "GB"))

    Int mem_calc = mem_mult*ceil(size(end1, "GB") + size(end2, "GB"))
    Int mem_full = if mem_calc < 16 then 16 else mem_calc
    Int mem = mem_full - 2

    String f = select_first([unpaired, ""])

    command <<<
        set -euxo pipefail

        if [ -z "~{f}" ]
        then
            mccortex63 build -m ~{mem}G -k 47 -s ~{sample_name} -2 ~{end1}:~{end2} -S ~{sample_name}.raw.ctx
        else
            gsutil cp ~{f} se.fq.gz
            mccortex63 build -m ~{mem}G -k 47 -s ~{sample_name} -1 se.fq.gz -2 ~{end1}:~{end2} -S ~{sample_name}.raw.ctx
        fi
    >>>

    output {
        File raw_ctx = "~{sample_name}.raw.ctx"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             mem_full,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "quay.io/corticall/mccortex:0.2.0"
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

task Clean {
    input {
        File raw_ctx
        String sample_name

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 4*ceil(size(raw_ctx, "GB"))
    Int mem_calc = 4*ceil(size(raw_ctx, "GB"))
    Int mem_full = if mem_calc < 16 then 16 else mem_calc
    Int mem = mem_full - 2

    command <<<
        set -euxo pipefail

        mccortex63 clean \
            -m ~{mem}G \
            -S \
            --tips \
            --unitigs \
            --fallback 1 \
            --covg-before ~{sample_name}.covg_before.csv \
            --covg-after ~{sample_name}.covg_after.csv \
            --len-before ~{sample_name}.len_before.csv \
            --len-after ~{sample_name}.len_after.csv \
            -o ~{sample_name}.clean.ctx \
            ~{raw_ctx}
    >>>

    output {
        File clean_ctx   = "~{sample_name}.clean.ctx"
        File covg_before = "~{sample_name}.covg_before.csv"
        File covg_after  = "~{sample_name}.covg_after.csv"
        File len_before  = "~{sample_name}.len_before.csv"
        File len_after   = "~{sample_name}.len_after.csv"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             mem_full,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "quay.io/corticall/mccortex:0.2.0"
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

task RecoverExcludedKmers {
    input {
        Array[File] ref_ctxs
        File dirty_ctx
        String sample_name

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 4*ceil(size(ref_ctxs, "GB") + size(dirty_ctx, "GB"))
    Int mem = 8

    command <<<
        set -euxo pipefail

    >>>

    output {
        File infer_ctx = "~{sample_name}.inferedges.ctx"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             mem,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "quay.io/corticall/mccortex:0.2.0"
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

task InferEdges {
    input {
        File clean_ctx
        String sample_name

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 4*ceil(size(clean_ctx, "GB"))
    Int mem_full = 2*ceil(size(clean_ctx, "GB"))
    Int mem = mem_full - 1

    command <<<
        set -euxo pipefail

        mccortex63 inferedges -m ~{mem}G -o ~{sample_name}.inferedges.ctx ~{clean_ctx}
    >>>

    output {
        File infer_ctx = "~{sample_name}.inferedges.ctx"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             mem_full,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "quay.io/corticall/mccortex:0.2.0"
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

task ThreadSe {
    input {
        File infer_ctx
        String sample_name
        Array[File] ends

        RuntimeAttr? runtime_attr_override
    }

    Int num_cpus = 8
    Int disk_size = 4*ceil(size(infer_ctx, "GB"))
    Int mem_full = 8*ceil(size(infer_ctx, "GB"))
    Int mem = mem_full - 1

    command <<<
        set -euxo pipefail

        mccortex63 thread -m ~{mem}G -t ~{num_cpus} -W -E -1 ~{sep=' -1 ' ends} -o ~{sample_name}.se.ctp.gz ~{infer_ctx}
    >>>

    output {
        File ctp = "~{sample_name}.se.ctp.gz"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          num_cpus,
        mem_gb:             mem_full,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "quay.io/corticall/mccortex:0.2.0"
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

task ThreadPe {
    input {
        File infer_ctx
        String sample_name
        File end1
        File end2
        File? unpaired
        Array[File] ctps

        RuntimeAttr? runtime_attr_override
    }

    Int num_cpus = 16
    Int disk_size = 4*ceil(size(infer_ctx, "GB"))
    Int mem_full = 8*ceil(size(infer_ctx, "GB"))
    Int mem = mem_full - 1

    String f = select_first([unpaired, ""])

    command <<<
        set -euxo pipefail

        if [ -z "~{f}" ]
        then
            mccortex63 thread -m ~{mem}G -t ~{num_cpus} -W -E -2 ~{end1}:~{end2} -p ~{sep=' -p ' ctps} -o ~{sample_name}.pe.ctp.gz ~{infer_ctx}
        else
            gsutil cp ~{f} se.fq.gz
            mccortex63 thread -m ~{mem}G -t ~{num_cpus} -W -E -1 se.fq.gz -2 ~{end1}:~{end2} -p ~{sep=' -p ' ctps} -o ~{sample_name}.pe.ctp.gz ~{infer_ctx}
        fi
    >>>

    output {
        File ctp = "~{sample_name}.pe.ctp.gz"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          num_cpus,
        mem_gb:             mem_full,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "quay.io/corticall/mccortex:0.2.0"
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

task ThreadRef {
    input {
        File infer_ctx
        String sample_name
        String ref_name
        File ref

        RuntimeAttr? runtime_attr_override
    }

    Int num_cpus = 8
    Int disk_size = 4*ceil(size(infer_ctx, "GB"))
    Int mem_full = 8*ceil(size(infer_ctx, "GB"))
    Int mem = mem_full - 1

    command <<<
        set -euxo pipefail

        mccortex63 thread -m ~{mem}G -t ~{num_cpus} -W -E -1 ~{sep=' -1 ' ref} -o ~{sample_name}.~{ref_name}.ref.ctp.gz ~{infer_ctx}
    >>>

    output {
        File ctp = "~{sample_name}.~{ref_name}.ref.ctp.gz"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          num_cpus,
        mem_gb:             mem_full,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "quay.io/corticall/mccortex:0.2.0"
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

task Contigs {
    input {
        File infer_ctx
        String sample_name
        Array[File] ctps

        RuntimeAttr? runtime_attr_override
    }

    Int num_cpus = 4
    Int disk_size = 4*ceil(size(infer_ctx, "GB"))
    Int mem_full = 8*ceil(size(infer_ctx, "GB"))
    Int mem = mem_full - 1

    command <<<
        set -euxo pipefail

        mccortex63 contigs -m ~{mem}G -t ~{num_cpus} -p ~{sep=' -p ' ctps} -G 23000000 -M -o ~{sample_name}.contigs.fa ~{infer_ctx}
    >>>

    output {
        File contigs = "~{sample_name}.contigs.fa"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          num_cpus,
        mem_gb:             mem_full,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "quay.io/corticall/mccortex:0.2.0"
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

task BuildFromRef {
    input {
        File ref
        String sample_name

        Int mem_mult = 8

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 10*ceil(size(select_all(ref, "GB")))

    Int mem_calc = mem_mult*ceil(size(ref, "GB"))
    Int mem_full = if mem_calc < 16 then 16 else mem_calc
    Int mem = mem_full - 2

    command <<<
        set -euxo pipefail

        mccortex63 build -m ~{mem}G -k 47 -s ~{sample_name} -1 ~{ref} -S ~{sample_name}.ref.ctx
    >>>

    output {
        File ref_ctx = "~{sample_name}.ref.ctx"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             mem_full,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "quay.io/corticall/mccortex:0.2.0"
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
