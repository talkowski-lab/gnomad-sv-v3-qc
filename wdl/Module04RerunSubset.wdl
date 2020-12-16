version 1.0

import "Module04.wdl" as module04
import "Tasks0506.wdl" as tasks0506

workflow Module04RerunSubset {
  input {

    String cohort
    String batch

    # Newline-deliminted list of variant IDs to regenotype with Module 04
    File vids_list

    # Cohort vcfs, which define variants to be genotyped with Module 04
    File cohort_pesr_vcf
    File cohort_depth_vcf

    # Original genotyped Module 04 outputs for this batch
    File genotyped_batch_pesr_vcf
    File genotyped_batch_depth_vcf

    String sv_base_mini_docker
  }

  # These two will call cache
  call SubsetVcfByVID as SubsetCohortDepth {
    input:
      vcf=cohort_depth_vcf,
      vcf_index=cohort_depth_vcf + ".tbi",
      vids_list=vids_list,
      output_prefix=cohort + ".merged_depth",
      sv_base_mini_docker=sv_base_mini_docker
  }

  call SubsetVcfByVID as SubsetCohortPesr {
    input:
      vcf=cohort_pesr_vcf,
      vcf_index=cohort_pesr_vcf + ".tbi",
      vids_list=vids_list,
      output_prefix=cohort + ".merged_pesr",
      sv_base_mini_docker=sv_base_mini_docker
  }

  call module04.Module04 {
    input:
      batch=batch,
      cohort_pesr_vcf=SubsetCohortPesr.filtered_vcf,
      cohort_depth_vcf=SubsetCohortDepth.filtered_vcf,
      sv_base_mini_docker=sv_base_mini_docker
  }

  call SubsetVcfByVID as SubsetGenotypedBatchDepth {
    input:
      vcf=genotyped_batch_depth_vcf,
      vcf_index=genotyped_batch_depth_vcf + ".tbi",
      vids_list=vids_list,
      output_prefix=batch + ".genotyped_batch_depth",
      sv_base_mini_docker=sv_base_mini_docker
  }

  call SubsetVcfByVID as SubsetGenotypedBatchPesr {
    input:
      vcf=genotyped_batch_pesr_vcf,
      vcf_index=genotyped_batch_pesr_vcf + ".tbi",
      vids_list=vids_list,
      output_prefix=batch + ".genotyped_batch_pesr",
      sv_base_mini_docker=sv_base_mini_docker
  }

  call tasks0506.ConcatVcfs as ConcatDepthVcfs {
    input:
      vcfs = [Module04.genotyped_depth_vcf, SubsetGenotypedBatchDepth.complement_vcf],
      vcfs_idx = [Module04.genotyped_depth_vcf + ".tbi", SubsetGenotypedBatchDepth.complement_vcf + ".tbi"],
      outfile_prefix=batch + ".module04_rerun_subset.depth",
      sv_base_mini_docker=sv_base_mini_docker
  }

  call tasks0506.ConcatVcfs as ConcatPesrVcfs {
    input:
      vcfs = [Module04.genotyped_pesr_vcf, SubsetGenotypedBatchPesr.complement_vcf],
      vcfs_idx = [Module04.genotyped_pesr_vcf + ".tbi", SubsetGenotypedBatchPesr.complement_vcf + ".tbi"],
      outfile_prefix=batch + ".module04_rerun_subset.pesr",
      sv_base_mini_docker=sv_base_mini_docker
  }

  output {
    File sr_bothside_pass = Module04.sr_bothside_pass
    File sr_background_fail = Module04.sr_background_fail
    File regeno_coverage_medians = Module04.regeno_coverage_medians

    File genotyped_depth_vcf = ConcatDepthVcfs.concat_vcf
    File genotyped_depth_vcf_index = ConcatDepthVcfs.concat_vcf_idx

    File genotyped_pesr_vcf = ConcatPesrVcfs.concat_vcf
    File genotyped_pesr_vcf_index = ConcatPesrVcfs.concat_vcf_idx
  }
}

task SubsetVcfByVID {
  input {
    File vcf
    File vcf_index
    File vids_list
    String output_prefix
    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_override
  }
  String output_filename_1 = output_prefix + ".filtered.vcf.gz"
  String output_filename_2 = output_prefix + ".complement.vcf.gz"

  Float input_size = size(vcf, "GB")
  RuntimeAttr runtime_default = object {
                                  mem_gb: 3.75,
                                  disk_gb: ceil(10.0 + input_size * 2.0),
                                  cpu_cores: 1,
                                  preemptible_tries: 3,
                                  max_retries: 1,
                                  boot_disk_gb: 10
                                }
  RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])
  runtime {
    memory: "~{select_first([runtime_override.mem_gb, runtime_default.mem_gb])} GB"
    disks: "local-disk ~{select_first([runtime_override.disk_gb, runtime_default.disk_gb])} HDD"
    cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
    preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
    maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
    docker: sv_base_mini_docker
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
  }

  command <<<
    set -euxo pipefail
    bcftools view -i 'ID==@~{vids_list}' ~{vcf} -O z -o ~{output_filename_1}
    tabix ~{output_filename_1}
    bcftools view -e 'ID==@~{vids_list}' ~{vcf} -O z -o ~{output_filename_2}
    tabix ~{output_filename_2}
  >>>

  output {
    File filtered_vcf = output_filename_1
    File filtered_vcf_index = output_filename_1 + ".tbi"
    File complement_vcf = output_filename_2
    File complement_vcf_index = output_filename_2 + ".tbi"
  }
}
