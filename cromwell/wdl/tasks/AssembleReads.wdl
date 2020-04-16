version 1.0

import "Structs.wdl"

workflow Assemble {
    input {
        File end1
        File end2
        File? unpaired
        String sample_name
        Array[Object] parents
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
            ref_ctxs    = ref_ctxs,
            clean_ctx   = Clean.clean_ctx,
            dirty_ctx   = Build.raw_ctx,
            sample_name = sample_name
    }

    call InferEdges {
        input:
            clean_ctx   = RecoverExcludedKmers.recovered_ctx,
            sample_name = sample_name
    }

    call PopBubbles {
        input:
            infer_ctx   = InferEdges.infer_ctx,
            sample_name = sample_name
    }

    call Contigs as ContigsWithoutLinks {
        input:
            popped_ctx  = PopBubbles.popped_ctx,
            label       = "wo_links",
            sample_name = sample_name
    }

    call ThreadSe {
        input:
            popped_ctx  = PopBubbles.popped_ctx,
            sample_name = sample_name,
            ends        = select_all([ end1, end2, unpaired ])
    }

    call ThreadPe {
        input:
            popped_ctx  = PopBubbles.popped_ctx,
            sample_name = sample_name,
            end1        = end1,
            end2        = end2,
            ctps        = [ ThreadSe.ctp ]
    }

    scatter (entry in parents) {
        call ThreadRef {
            input:
                popped_ctx  = PopBubbles.popped_ctx,
                sample_name = sample_name,
                ref_name    = entry['ref_name'],
                ref         = entry['ref']
        }
    }

    call Contigs as ContigsWithSeLinks {
        input:
            popped_ctx  = PopBubbles.popped_ctx,
            sample_name = sample_name,
            label       = "w_se_links",
            ctps        = [ ThreadSe.ctp ]
    }

    call Contigs as ContigsWithSeAndPeLinks {
        input:
            popped_ctx  = PopBubbles.popped_ctx,
            sample_name = sample_name,
            label       = "w_se_and_pe_links",
            ctps        = [ ThreadSe.ctp, ThreadPe.ctp ]
    }

    call Contigs as ContigsWithAllLinks {
        input:
            popped_ctx  = PopBubbles.popped_ctx,
            sample_name = sample_name,
            label       = "w_all_links",
            ctps        = flatten([ [ThreadSe.ctp, ThreadPe.ctp], ThreadRef.ctp ]),
    }

    output {
        File final_ctx = PopBubbles.popped_ctx
        File se_ctp = ThreadSe.ctp
        File pe_ctp = ThreadPe.ctp
        Array[File] ref_ctps = ThreadRef.ctp
        File contigs_without_links = ContigsWithoutLinks.contigs
        File contigs_with_se_links = ContigsWithSeLinks.contigs
        File contigs_with_se_and_pe_links = ContigsWithSeAndPeLinks.contigs
        File contigs_with_all_links = ContigsWithAllLinks.contigs
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

    Array[File] fqs = select_all([end1, end2, unpaired])

    Int disk_size = 10*ceil(size(fqs, "GB"))
    Int mem_calc = mem_mult*ceil(size(end1, "GB") + size(end2, "GB"))
    Int mem_full = if mem_calc < 16 then 16 else mem_calc
    Int mem = mem_full - 2

    command <<<
        set -euxo pipefail

        mccortex63 build -m ~{mem}G -k 47 -s ~{sample_name} -1 ~{sep=' -1 ' fqs} -S ~{sample_name}.raw.ctx
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
        docker:             "quay.io/corticall/mccortex:0.3.0"
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
        docker:             "quay.io/corticall/mccortex:0.3.0"
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
        File clean_ctx
        File dirty_ctx
        String sample_name

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 4*ceil(size(ref_ctxs, "GB") + size(dirty_ctx, "GB"))
    Int mem_full = 8
    Int mem = mem_full - 1

    command <<<
        set -euxo pipefail

        java -Xmx~{mem}g -jar /CortexJDK/dist/cortexjdk.jar Join -g ~{clean_ctx} -g ~{sep=' -g ' ref_ctxs} -o pedigree.ctx
        java -Xmx~{mem}g -jar /CortexJDK/dist/cortexjdk.jar RecoverExcludedKmers -g pedigree.ctx -d ~{dirty_ctx} -o ~{sample_name}.recovered.ctx
    >>>

    output {
        File recovered_ctx = "~{sample_name}.recovered.ctx"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             mem_full,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "quay.io/corticall/corticall:0.1.0"
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
    Int mem_calc = 2*ceil(size(clean_ctx, "GB"))
    Int mem_full = if mem_calc < 4 then 4 else mem_calc
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
        docker:             "quay.io/corticall/mccortex:0.3.0"
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

task PopBubbles {
    input {
        File infer_ctx
        String sample_name

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 4*ceil(size(infer_ctx, "GB"))
    Int mem_full = 2*ceil(size(infer_ctx, "GB"))
    Int mem = mem_full - 1

    command <<<
        set -euxo pipefail

        mccortex63 popbubbles -m ~{mem}G -o ~{sample_name}.popped.ctx ~{infer_ctx}
    >>>

    output {
        File popped_ctx = "~{sample_name}.popped.ctx"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             mem_full,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "quay.io/corticall/mccortex:0.3.0"
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
        File popped_ctx
        String sample_name
        Array[File] ends

        RuntimeAttr? runtime_attr_override
    }

    Int num_cpus = 8
    Int disk_size = 4*ceil(size(popped_ctx, "GB") + size(ends, "GB"))
    Int mem_full = 8*ceil(size(popped_ctx, "GB"))
    Int mem = mem_full - 1

    command <<<
        set -euxo pipefail

        mccortex63 thread -m ~{mem}G -t ~{num_cpus} -W -E -1 ~{sep=' -1 ' ends} -o ~{sample_name}.se.ctp.gz ~{popped_ctx}
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
        docker:             "quay.io/corticall/mccortex:0.3.0"
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
        File popped_ctx
        String sample_name
        File end1
        File end2
        Array[File] ctps

        RuntimeAttr? runtime_attr_override
    }

    Int num_cpus = 16
    Int disk_size = 4*ceil(size(popped_ctx, "GB") + size(end1, "GB") + size(end2, "GB") + size(ctps, "GB"))
    Int mem_full = 8*ceil(size(popped_ctx, "GB"))
    Int mem = mem_full - 1

    command <<<
        set -euxo pipefail

        mccortex63 thread -m ~{mem}G -t ~{num_cpus} -W -E -2 ~{end1}:~{end2} -p ~{sep=' -p ' ctps} -o ~{sample_name}.pe.ctp.gz ~{popped_ctx}

        mccortex63 links -T link.stats.txt -L 1000 ~{sample_name}.pe.ctp.gz
        LINK_THRESH=$(grep 'suggested_cutoff=' link.stats.txt | grep -oE '[0-9,]+$')
        mccortex63 links --clean $LINK_THRESH --out ~{sample_name}.pe.clean.ctp.gz ~{sample_name}.pe.ctp.gz
    >>>

    output {
        File ctp = "~{sample_name}.pe.clean.ctp.gz"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          num_cpus,
        mem_gb:             mem_full,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "quay.io/corticall/mccortex:0.3.0"
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
        File popped_ctx
        String sample_name
        String ref_name
        File ref

        RuntimeAttr? runtime_attr_override
    }

    Int num_cpus = 8
    Int disk_size = 4*ceil(size(popped_ctx, "GB"))
    Int mem_full = 8*ceil(size(popped_ctx, "GB"))
    Int mem = mem_full - 1

    command <<<
        set -euxo pipefail

        mccortex63 thread -m ~{mem}G -t ~{num_cpus} -W -E -1 ~{sep=' -1 ' ref} -o ~{sample_name}.~{ref_name}.ref.ctp.gz ~{popped_ctx}
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
        docker:             "quay.io/corticall/mccortex:0.3.0"
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
        File popped_ctx
        Array[File] ctps = []
        String sample_name
        String label

        RuntimeAttr? runtime_attr_override
    }

    Int num_cpus = 16
    Int disk_size = 4*ceil(size(popped_ctx, "GB"))
    Int mem_full = 8*ceil(size(popped_ctx, "GB"))
    Int mem = mem_full - 1

    command <<<
        set -euxo pipefail

        if [ "~{length(ctps)}" -eq "0" ]
        then
            mccortex63 contigs -m ~{mem}G -t ~{num_cpus} -G 23332839 -M -o ~{sample_name}.contigs.fa ~{popped_ctx}
        else
            mccortex63 contigs -m ~{mem}G -t ~{num_cpus} -p ~{sep=' -p ' ctps} -G 23332839 -M -o ~{sample_name}.contigs.fa ~{popped_ctx}
        fi

        cd-hit-est -T 0 -M 4000 -c 0.95 -i ~{sample_name}.contigs.fa -o ~{sample_name}.~{label}.contigs.dedup.fa
    >>>

    output {
        File contigs = "~{sample_name}.~{label}.contigs.dedup.fa"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          num_cpus,
        mem_gb:             mem_full,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "quay.io/corticall/mccortex:0.3.0"
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
        String ref_name

        Int mem_mult = 8

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 10*ceil(size(ref, "GB"))

    Int mem_calc = mem_mult*ceil(size(ref, "GB"))
    Int mem_full = if mem_calc < 16 then 16 else mem_calc
    Int mem = mem_full - 2

    command <<<
        set -euxo pipefail

        mccortex63 build -m ~{mem}G -k 47 -s ~{ref_name} -1 ~{ref} -S ~{ref_name}.ref.ctx
    >>>

    output {
        File ref_ctx = "~{ref_name}.ref.ctx"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             mem_full,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  0,
        max_retries:        0,
        docker:             "quay.io/corticall/mccortex:0.3.0"
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
