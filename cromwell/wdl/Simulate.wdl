version 1.0

import "tasks/Find.wdl" as Find
import "tasks/Utils.wdl" as Utils
import "tasks/PreprocessReads.wdl" as PR
import "tasks/AlignReads.wdl" as AR
import "tasks/AssembleReads.wdl" as ASM
import "tasks/AlignedMetrics.wdl" as AM
import "tasks/AssemblyMetrics.wdl" as ASMM
import "tasks/Finalize.wdl" as FF

workflow Simulate {
    input {
        Int num_samples = 1000
        #Array[Int] covs = [ 120, 100, 80, 60, 40, 30, 20, 10, 5 ]
        Array[Int] covs = [ 120 ]
        Object refs
        String gcs_out_root_dir
    }

    String outdir = sub(gcs_out_root_dir, "/$", "")

    call PrepareParentResources as PrepareHB3 { input: dict = refs['HB3']['dict'], gff = refs['HB3']['gff'] }
    call PrepareParentResources as PrepareDD2 { input: dict = refs['DD2']['dict'], gff = refs['DD2']['gff'] }

    scatter (index in range(1)) {
        call SimulateHaploidChild {
            input:
                p = [ "HB3", "DD2" ],
                r = [ refs['HB3']['fasta'], refs['DD2']['fasta'] ],
                f = [ refs['HB3']['fai'],   refs['DD2']['fai'] ],
                g = [ refs['HB3']['gff'],   refs['DD2']['gff'] ],
                d = [ refs['HB3']['dict'],  refs['DD2']['dict'] ],
                c = [ PrepareHB3.chrs,      PrepareDD2.chrs ],
                m = [ 0.499806172257622, 0.807757414029959, 0.929297453545848, 1.06255217143551, 1.20641344962446,
                      1.2815130851906,   1.30862778925489,  1.33637900760334,  1.40569169260802, 1.55242281561059,
                      1.90505374014534,  2.13950216259681,  2.79687360203287,  3.16560944406413 ],
                i = index
        }

        call SimulatePerfectReads as ReadsParent0 { input: fasta = SimulateHaploidChild.p0_fasta }
        call AddQuals as AddQualsP0E1 { input: fa = ReadsParent0.end1 }
        call AddQuals as AddQualsP0E2 { input: fa = ReadsParent0.end2 }

        call SimulatePerfectReads as ReadsParent1 { input: fasta = SimulateHaploidChild.p1_fasta }
        call AddQuals as AddQualsP1E1 { input: fa = ReadsParent1.end1 }
        call AddQuals as AddQualsP1E2 { input: fa = ReadsParent1.end2 }

        call GetTruthROIs {
            input:
                child_fasta = SimulateHaploidChild.child_fasta,
                parent0_fasta = refs['HB3']['fasta'],
                parent1_fasta = refs['DD2']['fasta'],
        }

        call VerifyTruthROIs {
            input:
                rois_truth_ctx = GetTruthROIs.rois_truth_ctx,
                child_novelkmers = SimulateHaploidChild.child_novelkmers
        }

        scatter (cov in covs) {
            call SimulateReads { input: child_fasta = SimulateHaploidChild.child_fasta, cov = cov, index = index }
            call BuildGraph { input: end1 = SimulateReads.end1, end2 = SimulateReads.end2, sample_name = "child" + index }
            call ThreadReads { input: end1 = SimulateReads.end1, end2 = SimulateReads.end2, ctx = BuildGraph.ctx }
            call ThreadRef as ThreadRef0 { input: ref = SimulateHaploidChild.p0_fasta, ctx = BuildGraph.ctx, ref_name = "HB3" }
            call ThreadRef as ThreadRef1 { input: ref = SimulateHaploidChild.p1_fasta, ctx = BuildGraph.ctx, ref_name = "DD2" }
            call Join { input: HB3_ctx = refs['HB3']['ctx'], DD2_ctx = refs['DD2']['ctx'], child_ctx = BuildGraph.ctx }

            call FindROIs { input: trio_ctx = Join.trio_ctx, sample_name = "child" + index }
            call FindOrphans { input: trio_ctx = Join.trio_ctx, rois_ctx = FindROIs.rois_ctx }
            call FindTips { input: trio_ctx = Join.trio_ctx, rois_ctx = FindROIs.rois_ctx }
            call FindDust { input: trio_ctx = Join.trio_ctx, rois_ctx = FindROIs.rois_ctx }
            call FindLowCoverage { input: rois_ctx = FindROIs.rois_ctx }
            call FindLowComplexity { input: trio_ctx = Join.trio_ctx, rois_ctx = FindROIs.rois_ctx }

            call RemoveKmers {
                input:
                    rois_ctx = FindROIs.rois_ctx,
                    bad_ctx = [ FindOrphans.ctx, FindTips.ctx, FindDust.ctx, FindLowCoverage.ctx, FindLowComplexity.ctx ]
            }

            call Partition as PartitionWithoutLinks {
                input:
                    child_ctx = BuildGraph.ctx,
                    rois_filtered_ctx = RemoveKmers.rois_filtered_ctx
            }

            call Partition as PartitionWithLinks {
                input:
                    child_ctx = BuildGraph.ctx,
                    links = [ ThreadReads.links, ThreadRef0.links, ThreadRef1.links ],
                    rois_filtered_ctx = RemoveKmers.rois_filtered_ctx
            }

            call Corticall as CorticallWithoutLinks {
                input:
                    joined = Join.trio_ctx,
                    rois = RemoveKmers.rois_filtered_ctx,
                    partitions = PartitionWithoutLinks.partitions,
                    refs = refs
            }

            call Corticall as CorticallWithLinks {
                input:
                    joined = Join.trio_ctx,
                    rois = RemoveKmers.rois_filtered_ctx,
                    partitions = PartitionWithLinks.partitions,
                    links = [ ThreadReads.links, ThreadRef0.links, ThreadRef1.links ],
                    refs = refs
            }

            call JoinCtxs { input: ctxs = [ GetTruthROIs.rois_truth_ctx, RemoveKmers.rois_filtered_ctx ] }
            call ComputeVenn { input: joined_txt = JoinCtxs.joined_txt, sample = index, cov = cov }
            call EvaluateAccuracy {
                input:
                    child_ctx = BuildGraph.ctx,
                    child_novelkmers = SimulateHaploidChild.child_novelkmers,
                    rois_filtered_ctx = RemoveKmers.rois_filtered_ctx,
                    links = [ ThreadReads.links, ThreadRef0.links, ThreadRef1.links ],
                    sample = "child" + index,
                    cov = cov
            }

#            call FF.FinalizeToDir as FinalizeStats {
#                input:
#                    files  = [ EvaluateAccuracy.stats ],
#                    outdir = outdir + "/child" + index + "/" + cov
#            }

            call AddQuals as AddQualsEnd1 { input: fa = SimulateReads.end1 }
            call AddQuals as AddQualsEnd2 { input: fa = SimulateReads.end2 }

            String ID = "sim"
            String SM = "child" + index
            String PL = "ILLUMINA"
            String RG = "@RG\\tID:~{ID}\\tSM:~{SM}\\tPL:~{PL}"

            Array[File] end1s = [ AddQualsEnd1.fq, AddQualsP0E1.fq, AddQualsP1E1.fq ]
            Array[File] end2s = [ AddQualsEnd2.fq, AddQualsP0E2.fq, AddQualsP1E2.fq ]

            scatter (ref_name in [ 'HB3', 'DD2' ]) {
                scatter (p in zip(end1s, end2s)) {
                    call AR.BwaMem as BwaMem {
                        input:
                            end1 = p.left,
                            end2 = p.right,

                            ref_dict = refs[ref_name]['dict'],
                            ref_fasta = refs[ref_name]['fasta'],
                            ref_fasta_amb = refs[ref_name]['amb'],
                            ref_fasta_ann = refs[ref_name]['ann'],
                            ref_fasta_bwt = refs[ref_name]['bwt'],
                            ref_fasta_fai = refs[ref_name]['fai'],
                            ref_fasta_pac = refs[ref_name]['pac'],
                            ref_fasta_sa = refs[ref_name]['sa'],

                            RG = RG,
                    }

                    call Delly {
                        input:
                            bam = BwaMem.aligned_bam,
                            bai = BwaMem.aligned_bai,
                            ref = refs[ref_name]['fasta'],
                            label = ref_name
                    }

                    call Manta {
                        input:
                            bam = BwaMem.aligned_bam,
                            bai = BwaMem.aligned_bai,
                            ref = refs[ref_name]['fasta'],
                            fai = refs[ref_name]['fai'],
                            label = ref_name
                    }

#                    call Lumpy {
#                        input:
#                            bam = BwaMem.aligned_bam,
#                            bai = BwaMem.aligned_bai,
#                            ref = refs[ref_name]['fasta'],
#                            fai = refs[ref_name]['fai'],
#                    }

                    call HC {
                        input:
                            bam  = BwaMem.aligned_bam,
                            bai  = BwaMem.aligned_bai,
                            ref  = refs[ref_name]['fasta'],
                            dict = refs[ref_name]['dict'],
                            fai  = refs[ref_name]['fai'],
                            label = ref_name
                    }

                    call Gridss {
                        input:
                            bam = BwaMem.aligned_bam,
                            bai = BwaMem.aligned_bai,

                            ref_dict = refs[ref_name]['dict'],
                            ref_fasta = refs[ref_name]['fasta'],
                            ref_fasta_amb = refs[ref_name]['amb'],
                            ref_fasta_ann = refs[ref_name]['ann'],
                            ref_fasta_bwt = refs[ref_name]['bwt'],
                            ref_fasta_fai = refs[ref_name]['fai'],
                            ref_fasta_pac = refs[ref_name]['pac'],
                            ref_fasta_sa = refs[ref_name]['sa'],

                            label = ref_name
                    }
                }
            }
        }
    }
}

task SimulateHaploidChild {
    input {
        Array[String] p
        Array[File] r
        Array[File] f
        Array[File] g
        Array[File] d
        Array[File] c
        Array[Float] m
        Int i

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 1

    command <<<
        set -euxo pipefail

        java -Xmx8g -jar /usr/local/bin/corticall.jar SimulateHaploidChild \
            -p ~{p[0]} -p ~{p[1]} \
            -r1 ~{r[0]} -r2 ~{r[1]} \
            -g1 ~{g[0]} -g2 ~{g[1] }\
            -c1 ~{c[0]} -c2 ~{c[1]} \
            -m ~{sep=',' m} \
            -s ~{i} \
            -v 3 \
            -o child.fasta \
            -fo1 p0.fasta \
            -fo2 p1.fasta \
            -vo child.variants.txt \
            -to child.variants.vcf \
            -ko child.novelkmers.txt
    >>>

    output {
        File child_fasta = "child.fasta"
        File p0_fasta = "p0.fasta"
        File p1_fasta = "p1.fasta"
        File child_variants = "child.variants.txt"
        File child_vcf = "vhild.variants.vcf"
        File child_novelkmers = "child.novelkmers.txt"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             10,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        0,
        docker:             "quay.io/corticall/corticall:0.1.5"
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

task PrepareParentResources {
    input {
        File dict
        File gff

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 1

    command <<<
        set -euxo pipefail

        cat ~{dict} | grep -v -e '^\@HD' -e '_00' | awk '{ print $2 }' | sed 's/SN://' > chrs.list
        cat ~{gff} | grep -i pfemp1 | grep polypeptide | grep -v putative | awk '{ print $9 }' | sed 's/:/ /'g | awk '{ print $1 }' | sed 's/ID=//g' | sed 's/.1$//g' > ids.list
        ((grep '^##[a-z]' ~{gff}) && (grep -f ids.list ~{gff} | grep -v -e Parent -e polypeptide -e pseudogene)) > var.gff
    >>>

    output {
        File chrs = "chrs.list"
        File ids = "ids.list"
        File var_gff = "var.gff"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             10,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        0,
        docker:             "ubuntu:latest"
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

task GetTruthROIs {
    input {
        File child_fasta
        File parent0_fasta
        File parent1_fasta

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 1 + 4*ceil(size([child_fasta, parent0_fasta, parent1_fasta], "GB"))

    command <<<
        set -euxo pipefail

        mccortex63 build -m 10G -k47 -S \
            -s parent0 -1 ~{parent0_fasta} \
            -s parent1 -1 ~{parent1_fasta} \
            -s child -1 ~{child_fasta} \
            joined.ctx

        java -Xmx8g -jar /usr/local/bin/corticall.jar FindROIs \
            -g joined.ctx \
            -p parent0 \
            -p parent1 \
            -c child \
            -o rois_truth.ctx
    >>>

    output {
        File rois_truth_ctx = "rois_truth.ctx"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             12,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        0,
        docker:             "quay.io/corticall/corticall:0.1.5"
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

task VerifyTruthROIs {
    input {
        File rois_truth_ctx
        File child_novelkmers

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 1 + 4*ceil(size([rois_truth_ctx, child_novelkmers], "GB"))

    command <<<
        set -euxo pipefail

        awk '{ print $4 }' ~{child_novelkmers} > nk.txt

        mccortex63 build -m 2G -k47 -S -s nk -1 nk.txt nk.ctx
        mccortex63 join -m 10G -o joined.ctx ~{rois_truth_ctx} nk.ctx
        mccortex63 view -k joined.ctx > joined.txt
    >>>

    output {
        File joined_ctx = "joined.ctx"
        File joined_txt = "joined.txt"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             4,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        0,
        docker:             "quay.io/corticall/corticall:0.1.5"
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

task ListFiles {
    input {
        String gcs_input_dir

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 1

    command <<<
        set -euxo pipefail

        gsutil ls ~{gcs_input_dir}/**child.new.fasta > child_files.txt
        sed 's/child.new/parent0.new/' child_files.txt > parent0_files.txt
        sed 's/child.new/parent1.new/' child_files.txt > parent1_files.txt
        sed 's/fasta/variants.txt/' child_files.txt > variant_files.txt
        sed 's/fasta/novelkmers.txt/' child_files.txt > novelkmers_files.txt

    >>>

    output {
        Array[File] child = read_lines("child_files.txt")
        Array[File] parent0 = read_lines("parent0_files.txt")
        Array[File] parent1 = read_lines("parent1_files.txt")
        Array[File] variants = read_lines("variant_files.txt")
        Array[File] novelkmers = read_lines("novelkmers_files.txt")
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             1,
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

task SimulatePerfectReads {
    input {
        File fasta

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 1

    command <<<
        set -euxo pipefail

        readsim -r ~{fasta} -l 76 -i 250 -d 120 -g 0 -e 0 parent
    >>>

    output {
        File end1 = "parent.1.fa.gz"
        File end2 = "parent.2.fa.gz"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             4,
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

task SimulateReads {
    input {
        File child_fasta
        Int cov
        Int index

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 1

    command <<<
        set -euxo pipefail

        readsim -r ~{child_fasta} -l 76 -i 250 -d ~{cov} -g ~{index} -e 0.005 child
    >>>

    output {
        File end1 = "child.1.fa.gz"
        File end2 = "child.2.fa.gz"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             4,
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

task BuildGraph {
    input {
        File end1
        File end2
        String sample_name

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 10*ceil(size([end1, end2], "GB"))

    command <<<
        set -euxo pipefail

        mccortex63 build -m 10G -S -k 47 -s ~{sample_name} -2 ~{end1}:~{end2} ~{sample_name}.raw.ctx
        mccortex63 clean -m 10G -S -B 2 -o ~{sample_name}.clean.ctx ~{sample_name}.raw.ctx
        mccortex63 inferedges -m 10G -o ~{sample_name}.infer.ctx ~{sample_name}.clean.ctx
        mccortex63 sort -m 10G -o ~{sample_name}.sorted.ctx ~{sample_name}.infer.ctx
    >>>

    output {
        File ctx = "~{sample_name}.sorted.ctx"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             12,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
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

task ThreadReads {
    input {
        File end1
        File end2
        File ctx

        RuntimeAttr? runtime_attr_override
    }

    Int cpus = 4
    Int disk_size = 4*ceil(size([end1, end2, ctx], "GB"))

    command <<<
        set -euxo pipefail

        mccortex63 thread -m 10G -t ~{cpus} -1 ~{end1} -1 ~{end1} -o se.ctp.gz ~{ctx}
        mccortex63 thread -m 10G -t ~{cpus} -2 ~{end1}:~{end1} -p se.ctp.gz -o pe.ctp.gz ~{ctx}

        java -Xmx8g -jar /usr/local/bin/corticall.jar IndexLinks -s reads -l pe.ctp.gz
    >>>

    output {
        File links = "pe.ctp.bgz"
        File links_idx = "pe.ctp.bgz.idx"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          cpus,
        mem_gb:             12,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        0,
        docker:             "quay.io/corticall/corticall:0.1.5"
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
        File ref
        File ctx
        String ref_name

        RuntimeAttr? runtime_attr_override
    }

    Int cpus = 8
    Int disk_size = 10*ceil(size([ref, ctx], "GB"))

    command <<<
        set -euxo pipefail

        mccortex63 thread -m 10G -t ~{cpus} -1 ~{ref} -o ref.ctp.gz ~{ctx}
        java -Xmx8g -jar /usr/local/bin/corticall.jar IndexLinks -s ~{ref_name} -l ref.ctp.gz
    >>>

    output {
        File links = "ref.ctp.bgz"
        File links_idx = "ref.ctp.bgz.idx"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          cpus,
        mem_gb:             12,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        0,
        docker:             "quay.io/corticall/corticall:0.1.5"
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

task Join {
    input {
        File HB3_ctx
        File DD2_ctx
        File child_ctx

        RuntimeAttr? runtime_attr_override
    }

    Int cpus = 8
    Int disk_size = 4*ceil(size([HB3_ctx, DD2_ctx, child_ctx], "GB"))

    command <<<
        set -euxo pipefail

        java -Xmx8g -jar /usr/local/bin/corticall.jar Join -g ~{HB3_ctx} -g ~{DD2_ctx} -g ~{child_ctx} -o trio.ctx
    >>>

    output {
        File trio_ctx = "trio.ctx"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             10,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        0,
        docker:             "quay.io/corticall/corticall:0.1.5"
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

task FindROIs {
    input {
        File trio_ctx
        String sample_name

        RuntimeAttr? runtime_attr_override
    }

    Int cpus = 8
    Int disk_size = 4*ceil(size(trio_ctx, "GB"))

    command <<<
        set -euxo pipefail

        java -Xmx8g -jar /usr/local/bin/corticall.jar FindROIs -g ~{trio_ctx} -p HB3 -p DD2 -c ~{sample_name} -o rois.ctx
    >>>

    output {
        File rois_ctx = "rois.ctx"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             10,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        0,
        docker:             "quay.io/corticall/corticall:0.1.5"
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

task FindOrphans {
    input {
        File trio_ctx
        File rois_ctx

        RuntimeAttr? runtime_attr_override
    }

    Int cpus = 8
    Int disk_size = 4*ceil(size(trio_ctx, "GB"))

    command <<<
        set -euxo pipefail

        java -Xmx8g -jar /usr/local/bin/corticall.jar FindOrphans -g ~{trio_ctx} -p HB3 -p DD2 -r ~{rois_ctx} -o rois.orphans.ctx
    >>>

    output {
        File ctx = "rois.orphans.ctx"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             10,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        0,
        docker:             "quay.io/corticall/corticall:0.1.5"
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

task FindTips {
    input {
        File trio_ctx
        File rois_ctx

        RuntimeAttr? runtime_attr_override
    }

    Int cpus = 8
    Int disk_size = 4*ceil(size(trio_ctx, "GB"))

    command <<<
        set -euxo pipefail

        java -Xmx8g -jar /usr/local/bin/corticall.jar FindTips -g ~{trio_ctx} -p HB3 -p DD2 -r ~{rois_ctx} -o rois.tips.ctx
    >>>

    output {
        File ctx = "rois.tips.ctx"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             10,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        0,
        docker:             "quay.io/corticall/corticall:0.1.5"
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

task FindLowCoverage {
    input {
        File rois_ctx
        Int min = 10

        RuntimeAttr? runtime_attr_override
    }

    Int cpus = 8
    Int disk_size = 4*ceil(size(rois_ctx, "GB"))

    command <<<
        set -euxo pipefail

        java -Xmx8g -jar /usr/local/bin/corticall.jar FindLowCoverage -r ~{rois_ctx} -m ~{min} -o rois.low_coverage.ctx
    >>>

    output {
        File ctx = "rois.low_coverage.ctx"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             10,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        0,
        docker:             "quay.io/corticall/corticall:0.1.5"
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

task FindDust {
    input {
        File trio_ctx
        File rois_ctx

        RuntimeAttr? runtime_attr_override
    }

    Int cpus = 8
    Int disk_size = 4*ceil(size(trio_ctx, "GB"))

    command <<<
        set -euxo pipefail

        java -Xmx8g -jar /usr/local/bin/corticall.jar FindDust -g ~{trio_ctx} -p HB3 -p DD2 -r ~{rois_ctx} -o rois.dust.ctx
    >>>

    output {
        File ctx = "rois.dust.ctx"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             10,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        0,
        docker:             "quay.io/corticall/corticall:0.1.5"
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

task FindLowComplexity {
    input {
        File trio_ctx
        File rois_ctx

        Float threshold = 0.703

        RuntimeAttr? runtime_attr_override
    }

    Int cpus = 8
    Int disk_size = 4*ceil(size(trio_ctx, "GB"))

    command <<<
        set -euxo pipefail

        java -Xmx8g -jar /usr/local/bin/corticall.jar FindLowComplexity -g ~{trio_ctx} -p HB3 -p DD2 -r ~{rois_ctx} -t ~{threshold} -o rois.low_complexity.ctx
    >>>

    output {
        File ctx = "rois.low_complexity.ctx"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             10,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        0,
        docker:             "quay.io/corticall/corticall:0.1.5"
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

task RemoveKmers {
    input {
        File rois_ctx
        Array[File] bad_ctx

        RuntimeAttr? runtime_attr_override
    }

    Int cpus = 8
    Int disk_size = 4*ceil(size(rois_ctx, "GB") + size(bad_ctx, "GB"))

    command <<<
        set -euxo pipefail

        java -Xmx8g -jar /usr/local/bin/corticall.jar Remove -g ~{rois_ctx} -s ~{sep=' -s ' bad_ctx} -o rois.filtered.ctx
    >>>

    output {
        File rois_filtered_ctx = "rois.filtered.ctx"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             10,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        0,
        docker:             "quay.io/corticall/corticall:0.1.5"
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

task Partition {
    input {
        File child_ctx
        Array[File] links = []
        File rois_filtered_ctx

        RuntimeAttr? runtime_attr_override
    }

    Int cpus = 2
    Int disk_size = 1 + 10*ceil(size([child_ctx, rois_filtered_ctx], "GB"))
    Int mem_full = 10
    Int mem = mem_full - 1

    command <<<
        set -euxo pipefail

        mccortex63 view -k ~{rois_filtered_ctx} | awk '{ print $1 }' > novels.txt

        if [ ~{length(links)} == 0 ]
        then
            mccortex63 contigs -f -m ~{mem}G -M -s novels.txt -o contigs.fa ~{child_ctx}
        else
            mccortex63 contigs -f -m ~{mem}G -p ~{sep=' -p ' links} -M -s novels.txt -o contigs.fa ~{child_ctx}
        fi

        cd-hit-est -T 0 -M 4000 -c 0.95 -i contigs.fa -o contigs.dedup.fa
    >>>

    output {
        File partitions = "contigs.dedup.fa"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          cpus,
        mem_gb:             mem_full,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
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

task JoinCtxs {
    input {
        Array[File] ctxs

        RuntimeAttr? runtime_attr_override
    }

    Int cpus = 8
    Int disk_size = 4*ceil(size(ctxs, "GB"))

    command <<<
        set -euxo pipefail

        java -Xmx8g -jar /usr/local/bin/corticall.jar Join -g ~{sep=' -g ' ctxs} -o joined.ctx
        mccortex63 view -k joined.ctx > joined.txt
    >>>

    output {
        File joined_ctx = "joined.ctx"
        File joined_txt = "joined.txt"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          4,
        mem_gb:             10,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        0,
        docker:             "quay.io/corticall/corticall:0.1.5"
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

task ComputeVenn {
    input {
        File joined_txt
        String sample
        Int cov

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 1 + 4*ceil(size(joined_txt, "GB"))

    command <<<
        set -euxo pipefail

        UNIQUE_TO_TRUTH=$(awk '$2 > 0 && $3 == 0' ~{joined_txt} | wc -l)
        UNIQUE_TO_EVAL=$(awk '$2 == 0 && $3 > 0' ~{joined_txt} | wc -l)
        OVERLAP=$(awk '$2 > 0 && $3 > 0' ~{joined_txt} | wc -l)

        echo ~{sample} ~{cov} $UNIQUE_TO_TRUTH $UNIQUE_TO_EVAL $OVERLAP > venn.txt
    >>>

    output {
        File venn = "venn.txt"
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

task EvaluateAccuracy {
    input {
        File child_ctx
        File child_novelkmers
        File rois_filtered_ctx
        Array[File] links
        String sample
        Int cov

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 1 + 4*ceil(size([child_ctx, child_novelkmers], "GB"))
    Int mem_full = 8
    Int mem = mem_full - 1

    command <<<
        set -euxo pipefail

        TOTAL_NUM_ROIS=$(cat ~{child_novelkmers} | wc -l)

        for ID in $(awk '{ print $1 }' ~{child_novelkmers} | uniq)
        do
            TYPE=$(awk -v ID=$ID '{ if ($1 == ID) print $5 }' ~{child_novelkmers} | uniq)
            LENGTH=$(awk -v ID=$ID '{ if ($1 == ID) print length($9) }' ~{child_novelkmers} | uniq)

            awk -v ID=$ID '{ if ($1 == ID) print $4 }' ~{child_novelkmers} > exp_kmers.txt

            NUM_KMERS_EXPECTED=$(cat exp_kmers.txt | wc -l)

            mccortex63 build -f -m ~{mem}G -k47 -S -s truth -1 exp_kmers.txt truth.ctx
            #mccortex63 join -f -m ~{mem}G -o truth_and_filtered.ctx truth.ctx ~{rois_filtered_ctx}
            java -Xmx~{mem}g -jar /usr/local/bin/corticall.jar Join -o truth_and_filtered.ctx -g truth.ctx -g ~{rois_filtered_ctx}

            mccortex63 view -k truth_and_filtered.ctx | awk '{ if ($2 > 0 && $3 > 0) print $1 }' > novels.txt

            NUM_NOVELS=$(cat novels.txt | wc -l)

            mccortex63 contigs -f -m ~{mem}G -p ~{sep=' -p ' links} -M -s novels.txt -o contigs.fa ~{child_ctx}
            mccortex63 build -f -m ~{mem}G -k47 -S -s eval -1 contigs.fa eval.ctx
            #mccortex63 join -f -m ~{mem}G -o joint.ctx truth.ctx eval.ctx
            java -Xmx~{mem}g -jar /usr/local/bin/corticall.jar Join -o joint.ctx -g truth.ctx -g eval.ctx

            NUM_KMERS_FOUND=$(mccortex63 view -k joint.ctx | awk '$2 > 0 && $3 > 0' | wc -l)

            echo ~{sample} ~{cov} $TOTAL_NUM_ROIS $ID $TYPE $LENGTH $NUM_KMERS_EXPECTED $NUM_NOVELS $NUM_KMERS_FOUND >> stats.txt
        done
    >>>

    output {
        File stats = "stats.txt"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             mem_full,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        0,
        docker:             "quay.io/corticall/corticall:0.1.5"
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

task Corticall {
    input {
        File joined
        File rois
        File partitions
        Object refs
        Array[File] links = []

        RuntimeAttr? runtime_attr_override
    }

    File dict3D7 = refs["3D7"]["dict"]
    File ref3D7 = refs["3D7"]["fasta"]
    File amb3D7 = refs["3D7"]["amb"]
    File ann3D7 = refs["3D7"]["ann"]
    File bwt3D7 = refs["3D7"]["bwt"]
    File fai3D7 = refs["3D7"]["fai"]
    File pac3D7 = refs["3D7"]["pac"]
    File sa3D7 = refs["3D7"]["sa"]
    File ctx3D7 = refs["3D7"]["ctx"]
    File gff3D7 = refs["3D7"]["gff"]
    File so3D7 = refs["3D7"]["so"]

    File dictHB3 = refs["HB3"]["dict"]
    File refHB3 = refs["HB3"]["fasta"]
    File ambHB3 = refs["HB3"]["amb"]
    File annHB3 = refs["HB3"]["ann"]
    File bwtHB3 = refs["HB3"]["bwt"]
    File faiHB3 = refs["HB3"]["fai"]
    File pacHB3 = refs["HB3"]["pac"]
    File saHB3 = refs["HB3"]["sa"]
    File ctxHB3 = refs["HB3"]["ctx"]
    File gffHB3 = refs["HB3"]["gff"]
    File soHB3 = refs["HB3"]["so"]

    File dictDD2 = refs["DD2"]["dict"]
    File refDD2 = refs["DD2"]["fasta"]
    File ambDD2 = refs["DD2"]["amb"]
    File annDD2 = refs["DD2"]["ann"]
    File bwtDD2 = refs["DD2"]["bwt"]
    File faiDD2 = refs["DD2"]["fai"]
    File pacDD2 = refs["DD2"]["pac"]
    File saDD2 = refs["DD2"]["sa"]
    File ctxDD2 = refs["DD2"]["ctx"]
    File gffDD2 = refs["DD2"]["gff"]
    File soDD2 = refs["DD2"]["so"]

    Int disk_size = 1 + 10*ceil(size([joined, rois, partitions,
                                      dict3D7, ref3D7, amb3D7, ann3D7, bwt3D7, fai3D7, pac3D7, sa3D7, ctx3D7, gff3D7, so3D7,
                                      dictHB3, refHB3, ambHB3, annHB3, bwtHB3, faiHB3, pacHB3, saHB3, ctxHB3, gffHB3, soHB3,
                                      dictDD2, refDD2, ambDD2, annDD2, bwtDD2, faiDD2, pacDD2, saDD2, ctxDD2, gffDD2, soDD2], "GB"))
    Int mem_full = 8
    Int mem = mem_full - 1

    command <<<
        set -euxo pipefail

        mccortex63 view -k ~{rois} | wc -l
        grep -c '>' ~{partitions}

        java -Xmx~{mem}g -jar /usr/local/bin/corticall.jar Call \
            -g ~{joined} \
            -r ~{rois} \
            -p ~{partitions} \
            -b HB3 -b DD2 \
            -R 3D7:~{ref3D7} \
            -R HB3:~{refHB3} \
            -R DD2:~{refDD2} \
            -o corticall.vcf \
            -ao acct.txt
    >>>

    output {
        File vcf = "corticall.vcf"
        File acct = "acct.txt"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             mem_full,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        0,
        docker:             "quay.io/corticall/corticall:0.1.5"
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

task Delly {
    input {
        File bam
        File bai
        File ref
        String label

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 1 + 4*ceil(size([bam, ref], "GB"))

    command <<<
        set -euxo pipefail

        delly call -o sv.bcf -g ~{ref} ~{bam}
        bcftools view -O v -o delly.~{label}.vcf sv.bcf
    >>>

    output {
        File vcf = "delly.~{label}.vcf"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             4,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        0,
        docker:             "quay.io/corticall/sv:0.1.1"
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

task Manta {
    input {
        File bam
        File bai
        File ref
        File fai
        String label

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 1 + 4*ceil(size([bam, ref], "GB"))

    command <<<
        set -euxo pipefail

        configManta.py --bam ~{bam} --referenceFasta ~{ref} --runDir ./manta
        ./manta/runWorkflow.py

        zcat manta/results/variants/candidateSV.vcf.gz > manta.~{label}.vcf
    >>>

    output {
        File vcf = "manta.~{label}.vcf"
        File smallindels_vcf = "manta/results/variants/candidateSmallIndels.vcf.gz"
        File diploidsv_vcf = "manta/results/variants/diploidSV.vcf.gz"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          2,
        mem_gb:             16,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        0,
        docker:             "quay.io/corticall/sv:0.1.1"
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

task Lumpy {
    input {
        File bam
        File bai
        File ref
        File fai
        String label

        RuntimeAttr? runtime_attr_override
    }

    Int num_cpus = 2
    Int disk_size = 1 + 4*ceil(size([bam, ref], "GB"))

    command <<<
        set -euxo pipefail

        smoove call -x --name sim --fasta ~{ref} -p ~{num_cpus} --genotype ~{bam}
        bcftools view -O v -o lumpy.~{label}.vcf sim-smoove.genotyped.vcf.gz
    >>>

    output {
        File vcf = "lumpy.~{label}.vcf"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          num_cpus,
        mem_gb:             8,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        0,
        docker:             "brentp/smoove"
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

task HC {
    input {
        File bam
        File bai
        File ref
        File dict
        File fai
        String label

        RuntimeAttr? runtime_attr_override
    }

    Int num_cpus = 2
    Int disk_size = 1 + 4*ceil(size([bam, ref], "GB"))

    command <<<
        set -euxo pipefail

        gatk HaplotypeCaller -R ~{ref} -I ~{bam} -ploidy 1 -DF WellformedReadFilter -O hc.~{label}.vcf
    >>>

    output {
        File vcf = "hc.~{label}.vcf"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          num_cpus,
        mem_gb:             8,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        0,
        docker:             "quay.io/corticall/sv:0.1.1"
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

task Gridss {
    input {
        File bam
        File bai

        File ref_dict
        File ref_fasta
        File ref_fasta_amb
        File ref_fasta_ann
        File ref_fasta_bwt
        File ref_fasta_fai
        File ref_fasta_pac
        File ref_fasta_sa

        String label

        RuntimeAttr? runtime_attr_override
    }

    Int num_cpus = 2
    Int disk_size = 1 + 10*ceil(size([bam, bai, ref_fasta], "GB"))

    command <<<
        set -euxo pipefail

        /opt/gridss/gridss.sh --reference ~{ref_fasta} --output gridss.~{label}.vcf --assembly gridss.bam ~{bam}
    >>>

    output {
        File vcf = "gridss.~{label}.vcf"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          num_cpus,
        mem_gb:             8,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  2,
        max_retries:        0,
        docker:             "quay.io/corticall/gridss:0.1.0"
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

task AddQuals {
    input {
        File fa

        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 10*ceil(size(fa, "GB"))

    command <<<
        set -euxo pipefail

        FILE="~{fa}"

        if [[ "$FILE" =~ \.fasta$ ]] || [[ "$FILE" =~ \.fa$ ]]; then
            cat $FILE | python3 /usr/local/bin/cat_as_fastq.py | gzip > reads.fq.gz
        elif [[ "$FILE" =~ \.fasta.gz$ ]] || [[ "$FILE" =~ \.fa.gz$ ]]; then
            zcat $FILE | python3 /usr/local/bin/cat_as_fastq.py | gzip > reads.fq.gz
        fi
    >>>

    output {
        File fq = "reads.fq.gz"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             1,
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
