version 1.0

import "Structs.wdl"
import "Finalize.wdl" as FF

workflow AlignedMetrics {
    input {
        File aligned_bam
        File aligned_bai

        File ref_fasta
        File ref_dict

        String per
        String type
        String label

        String? gcs_output_dir
    }

    call MakeChrIntervalList { input: ref_dict = ref_dict }

    scatter (chr_info in MakeChrIntervalList.chrs) {
        call MosDepth {
            input:
                bam = aligned_bam,
                bai = aligned_bai,
                chr = chr_info[0]
        }

        call SummarizeDepth { input: regions = MosDepth.regions }
    }

    call FlagStats as AlignedFlagStats { input: bam = aligned_bam }

    call ReadLength { input: bam = aligned_bam }

    if (defined(gcs_output_dir)) {
        String outdir = sub(gcs_output_dir + "", "/$", "") + "/metrics/~{per}_~{type}_~{label}"

        call FF.FinalizeToDir as FFFlagStats  { input: outdir = outdir + "/flag_stats/", files = [ AlignedFlagStats.flag_stats ] }
        call FF.FinalizeToDir as FFReadLength { input: outdir = outdir + "/read_length/", files = [ ReadLength.read_length ] }

        call FF.FinalizeToDir as FFCoverageFullDist { input: outdir = outdir + "/coverage/", files = MosDepth.full_dist }
        call FF.FinalizeToDir as FFCoverageGlobalDist { input: outdir = outdir + "/coverage/", files = MosDepth.global_dist }
        call FF.FinalizeToDir as FFCoverageRegionDist { input: outdir = outdir + "/coverage/", files = MosDepth.region_dist }
        call FF.FinalizeToDir as FFCoverageRegions { input: outdir = outdir + "/coverage/", files = MosDepth.regions }
        call FF.FinalizeToDir as FFCoverageRegionsCsi { input: outdir = outdir + "/coverage/", files = MosDepth.regions_csi }
        call FF.FinalizeToDir as FFCoverageQuantizedDist { input: outdir = outdir + "/coverage/", files = MosDepth.quantized_dist }
        call FF.FinalizeToDir as FFCoverageQuantized { input: outdir = outdir + "/coverage/", files = MosDepth.quantized }
        call FF.FinalizeToDir as FFCoverageQuantizedCsi { input: outdir = outdir + "/coverage/", files = MosDepth.quantized_csi }

        call FF.FinalizeToDir as FFDepthSummaries { input: outdir = outdir + "/coverage_summaries/", files = SummarizeDepth.cov_summary }
    }

    output {
        File aligned_flag_stats = AlignedFlagStats.flag_stats
        File read_length        = ReadLength.read_length

        Array[File] coverage_full_dist      = MosDepth.full_dist
        Array[File] coverage_global_dist    = MosDepth.global_dist
        Array[File] coverage_region_dist    = MosDepth.region_dist
        Array[File] coverage_regions        = MosDepth.regions
        Array[File] coverage_regions_csi    = MosDepth.regions_csi
        Array[File] coverage_quantized_dist = MosDepth.quantized_dist
        Array[File] coverage_quantized      = MosDepth.quantized
        Array[File] coverage_quantized_csi  = MosDepth.quantized_csi
    }
}

task MakeChrIntervalList {
    input {
        File ref_dict

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 10

    command <<<
        set -euxo pipefail

        grep '^@SQ' ~{ref_dict} | awk '{ print $2 "\t" 1 "\t" $3 }' | sed 's/[SL]N://g' | grep -v -e random -e chrUn -e decoy -e alt -e HLA -e EBV > chrs.txt
    >>>

    output {
        Array[Array[String]] chrs = read_tsv("chrs.txt")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             1,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "quay.io/corticall/metrics:0.1.0"
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

task MosDepth {
    input {
        File bam
        File bai
        String chr
        Int? window_size

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 2*ceil(size(bam, "GB") + size(bai, "GB"))
    Int ws = select_first([window_size, 500])
    String basename = basename(bam, ".bam")
    String prefix = "~{basename}.coverage.~{chr}"

    command <<<
        set -euxo pipefail

        mosdepth -t 4 -c ~{chr} -n -x -Q 1 ~{prefix}.full ~{bam}
        mosdepth -t 4 -c ~{chr} -n -x -Q 1 -b ~{ws} ~{prefix} ~{bam}

        export MOSDEPTH_Q0=NO_COVERAGE   # 0 -- defined by the arguments to --quantize
        export MOSDEPTH_Q1=LOW_COVERAGE  # 1..4
        export MOSDEPTH_Q2=CALLABLE      # 5..149
        export MOSDEPTH_Q3=HIGH_COVERAGE # 150 ...

        mosdepth -t 4 -c ~{chr} -n -x -Q 1 --quantize 0:1:5:150: ~{prefix}.quantized ~{bam}

        ls -lah
    >>>

    output {
        File full_dist      = "~{prefix}.full.mosdepth.global.dist.txt"
        File global_dist    = "~{prefix}.mosdepth.global.dist.txt"
        File region_dist    = "~{prefix}.mosdepth.region.dist.txt"
        File regions        = "~{prefix}.regions.bed.gz"
        File regions_csi    = "~{prefix}.regions.bed.gz.csi"
        File quantized_dist = "~{prefix}.quantized.mosdepth.global.dist.txt"
        File quantized      = "~{prefix}.quantized.quantized.bed.gz"
        File quantized_csi  = "~{prefix}.quantized.quantized.bed.gz.csi"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             8,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "quay.io/biocontainers/mosdepth:0.2.4--he527e40_0"
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

task SummarizeDepth {
    input {
        File regions

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 2*ceil(size(regions, "GB"))
    String chrName = sub(basename(regions, ".regions.bed.gz"), "out.coverage.", "")

    command <<<
        set -euxo pipefail

        ((echo 'chr start stop cov_mean cov_sd cov_q1 cov_median cov_q3 cov_iqr') && \
         (zcat ~{regions} | datamash first 1 first 2 last 3 mean 4 sstdev 4 q1 4 median 4 q3 4 iqr 4)) | \
         column -t > ~{chrName}.summary.txt
    >>>

    output {
        File cov_summary = "~{chrName}.summary.txt"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             2,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "quay.io/corticall/metrics:0.1.0"
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

task CoverageTrack {
    input {
        File bam
        File bai
        String chr
        String start
        String end

        RuntimeAttr? runtime_attr_override
    }

    String basename = basename(bam, ".bam")
    Int disk_size = 2*ceil(size(bam, "GB") + size(bai, "GB"))

    command <<<
        set -euxo pipefail

        samtools depth -a ~{bam} -r ~{chr}:~{start}-~{end} | bgzip > ~{basename}.coverage.~{chr}_~{start}_~{end}.txt.gz
        tabix -p bed ~{basename}.coverage.~{chr}_~{start}_~{end}.txt.gz
    >>>

    output {
        File coverage = "~{basename}.coverage.~{chr}_~{start}_~{end}.txt.gz"
        File coverage_tbi = "~{basename}.coverage.~{chr}_~{start}_~{end}.txt.gz.tbi"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             4,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "quay.io/corticall/metrics:0.1.0"
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

task FlagStats {
    input {
        File bam

        RuntimeAttr? runtime_attr_override
    }

    String basename = basename(bam, ".bam")
    Int disk_size = 2*ceil(size(bam, "GB"))

    command <<<
        set -euxo pipefail

        samtools flagstat ~{bam} > ~{basename}.flag_stats.txt
    >>>

    output {
        File flag_stats = "~{basename}.flag_stats.txt"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             4,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "quay.io/corticall/metrics:0.1.0"
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

task ReadLength {
    input {
        String bam

        RuntimeAttr? runtime_attr_override
    }

    String basename = basename(bam, ".bam")
    Int disk_size = 2*ceil(size(bam, "GB"))

    command <<<
        export GCS_OAUTH_TOKEN=$(gcloud auth application-default print-access-token)
        samtools view ~{bam} | head -1 | awk '{ print length($10) }' > ~{basename}.read_length.txt
    >>>

    output {
        File read_length = "~{basename}.read_length.txt"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             4,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
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

task BamToBed {
    input {
        File bam
        File bai

        RuntimeAttr? runtime_attr_override
    }

    String bed = basename(bam, ".bam") + ".bed"
    Int disk_size = 4*ceil(size(bam, "GB") + size(bai, "GB"))

    command <<<
        set -euxo pipefail

        bedtools bamtobed -i ~{bam} > ~{bed}
    >>>

    output {
        File bedfile = bed
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             8,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        1,
        docker:             "quay.io/corticall/metrics:0.1.0"
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

