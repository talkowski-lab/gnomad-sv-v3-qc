version 1.0

# Author: Ryan Collins <rlcollins@g.harvard.edu>

import "Structs.wdl"
import "Tasks0506.wdl" as MiniTasks
import "ClusterSingleChromosome.wdl" as VcfClusterTasks

# Workflow to run parallelized vcf clustering for a single chromosome
workflow VcfClusterSingleChrom {
  input {
    Array[File] vcfs
    String prefix
    Int dist
    Float frac
    Float sample_overlap
    File? exclude_list
    Array[String] batches
    Int sv_size
    Array[String] sv_types
    String contig
    Int max_shards_per_chrom_svtype
    Int min_variants_per_shard_per_chrom_svtype
    Boolean subset_sr_lists
    File bothside_pass
    File background_fail
    File empty_file

    String sv_pipeline_docker
    String sv_base_mini_docker

    # overrides for local tasks
    RuntimeAttr? runtime_override_localize_vcfs
    RuntimeAttr? runtime_override_join_vcfs
    RuntimeAttr? runtime_override_fix_multiallelic
    RuntimeAttr? runtime_override_fix_ev_tags

    # overrides for MiniTasks
    RuntimeAttr? runtime_override_subset_bothside_pass
    RuntimeAttr? runtime_override_subset_background_fail

    # overrides for VcfClusterTasks
    RuntimeAttr? runtime_override_subset_sv_type
    RuntimeAttr? runtime_override_concat_sv_types
    RuntimeAttr? runtime_override_shard_vcf_precluster
    RuntimeAttr? runtime_override_svtk_vcf_cluster
    RuntimeAttr? runtime_override_get_vcf_header_with_members_info_line
    RuntimeAttr? runtime_override_concat_shards
  }

  scatter (i in range(length(vcfs))) {
    File vcf_indexes_ = vcfs[i] + ".tbi"
  }
  
  #Stream each vcf & join into a single vcf
  call LocalizeContigVcfs {
    input:
      vcfs=vcfs,
      vcf_indexes = vcf_indexes_,
      batches=batches,
      contig=contig,
      prefix=prefix,
      sv_pipeline_docker=sv_pipeline_docker,
      runtime_attr_override=runtime_override_localize_vcfs
  }
  call JoinVcfs {
    input:
      vcfs=LocalizeContigVcfs.out,
      contig=contig,
      prefix=prefix,
      sv_pipeline_docker=sv_pipeline_docker,
      runtime_attr_override=runtime_override_join_vcfs
  }
  call FixMultiallelicRecords {
    input:
      joined_vcf=JoinVcfs.out,
      batch_contig_vcfs=LocalizeContigVcfs.out,
      contig=contig,
      prefix=prefix,
      sv_pipeline_docker=sv_pipeline_docker,
      runtime_attr_override=runtime_override_fix_multiallelic
  }
  call FixEvidenceTags {
    input:
      vcf=FixMultiallelicRecords.out,
      contig=contig,
      prefix=prefix,
      sv_base_mini_docker=sv_base_mini_docker,
      runtime_attr_override=runtime_override_fix_ev_tags
  }

  #Run vcfcluster per chromosome
  call VcfClusterTasks.ClusterSingleChrom as ClusterSingleChrom {
    input:
      vcf=FixEvidenceTags.out,
      vcf_index=FixEvidenceTags.out_index,
      contig=contig,
      prefix=prefix,
      max_shards=max_shards_per_chrom_svtype,
      min_per_shard=min_variants_per_shard_per_chrom_svtype,
      dist=dist,
      frac=frac,
      sample_overlap=sample_overlap,
      exclude_list=exclude_list,
      sv_size=sv_size,
      sv_types=sv_types,
      sv_pipeline_docker=sv_pipeline_docker,
      sv_base_mini_docker=sv_base_mini_docker,
      runtime_override_subset_sv_type=runtime_override_subset_sv_type,
      runtime_override_concat_sv_types=runtime_override_concat_sv_types,
      runtime_override_shard_vcf_precluster=runtime_override_shard_vcf_precluster,
      runtime_override_svtk_vcf_cluster=runtime_override_svtk_vcf_cluster,
      runtime_override_get_vcf_header_with_members_info_line=runtime_override_get_vcf_header_with_members_info_line,
      runtime_override_concat_shards=runtime_override_concat_shards
  }

  String filtered_bothside_pass_name = prefix + "." + contig + ".pass.VIDs.list"
  String filtered_background_fail_name = prefix + "." + contig + ".fail.VIDs.list"
  if(subset_sr_lists) {
    #Subset bothside_pass & background_fail to chromosome of interest
    call MiniTasks.SubsetVariantList as SubsetBothsidePass {
      input:
        vid_list=bothside_pass,
        vcf=FixEvidenceTags.out,
        outfile_name=prefix + "." + contig + ".pass.VIDs.list",
        sv_base_mini_docker=sv_base_mini_docker,
        runtime_attr_override=runtime_override_subset_bothside_pass
    }
    call MiniTasks.SubsetVariantList as SubsetBackgroundFail {
      input:
        vid_list=background_fail,
        vcf=FixEvidenceTags.out,
        outfile_name=prefix + "." + contig + ".fail.VIDs.list",
        sv_base_mini_docker=sv_base_mini_docker,
        runtime_attr_override=runtime_override_subset_background_fail
    }
  }

  output {
    File clustered_vcf = ClusterSingleChrom.clustered_vcf
    File clustered_vcf_idx = ClusterSingleChrom.clustered_vcf_idx
    File filtered_bothside_pass = select_first([SubsetBothsidePass.filtered_vid_list, empty_file])
    File filtered_background_fail = select_first([SubsetBackgroundFail.filtered_vid_list, empty_file])
  }
}


# Shard batch VCFs, pulling down only this contig
task LocalizeContigVcfs {
  input {
    Array[File] vcfs
    Array[File] vcf_indexes
    Array[String] batches
    String contig
    String prefix
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  RuntimeAttr runtime_default = object {
    mem_gb: 3.75,
    disk_gb: ceil(10 + size(vcfs, "GiB") * 1.2),
    cpu_cores: 1,
    preemptible_tries: 0,
    max_retries: 1,
    boot_disk_gb: 10
  }
  RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])

  runtime {
    memory: select_first([runtime_override.mem_gb, runtime_default.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_override.disk_gb, runtime_default.disk_gb]) + " HDD"
    cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
    preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
    maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
    docker: sv_pipeline_docker
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
  }
  
  command <<<
    set -euxo pipefail

    # See Issue #52 "Use GATK to retrieve VCF records in JoinContigFromRemoteVcfs"
    # https://github.com/broadinstitute/gatk-sv/issues/52

    #Remote tabix all vcfs to chromosome of interest
    paste ~{write_lines(batches)} ~{write_lines(vcfs)} | while read BATCH VCF_PATH; do
      BATCH_VCF="$BATCH.~{contig}.subsetted.vcf.gz"
      tabix -h "$VCF_PATH" "~{contig}" \
        | sed "s/AN=[0-9]*;//g" \
        | sed "s/AC=[0-9]*;//g" \
        | bgzip \
        > $BATCH_VCF
    done
    VCFS_LIST="subsetted_vcfs.list"
    ls *.~{contig}.subsetted.vcf.gz > $VCFS_LIST

    #Sanity check to make sure all subsetted VCFs have same number of records
    # crazy ' || printf ""' statement to avoid pipefail if grep encounters no matching lines
    while read VCF; do
      zcat "$VCF" | (grep -Ev "^#" || printf "") | wc -l
    done < $VCFS_LIST \
    > records_per_vcf.txt

    if [ $( sort records_per_vcf.txt | uniq | wc -l ) -gt 1 ]; then
      1>&2 echo "ERROR: INCONSISTENT NUMBER OF RECORDS PER VCF DETECTED"
      cat records_per_vcf.txt
      exit 1
    fi
  >>>

  output {
    Array[File] out = glob("*.~{contig}.subsetted.vcf.gz")
  }
}

# Merge contig vcfs across batches
task JoinVcfs {
  input {
    Array[File] vcfs
    String contig
    String prefix
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  Float input_size = size(vcfs, "GiB")
  Float input_size_ratio = 3.0
  Float base_disk_gb = 10.0
  RuntimeAttr runtime_default = object {
                                  mem_gb: 1.0,
                                  disk_gb: ceil(base_disk_gb + input_size * input_size_ratio),
                                  cpu_cores: 1,
                                  preemptible_tries: 0,
                                  max_retries: 1,
                                  boot_disk_gb: 10
                                }
  RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])

  runtime {
    memory: select_first([runtime_override.mem_gb, runtime_default.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_override.disk_gb, runtime_default.disk_gb]) + " HDD"
    cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
    preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
    maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
    docker: sv_pipeline_docker
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
  }

  command <<<
    set -euo pipefail

    python3 <<CODE | bgzip > ~{prefix}.~{contig}.joined.vcf.gz
    import sys
    import gzip

    fl = open("~{write_lines(vcfs)}")
    files = [gzip.open(f.strip(), 'rb') for f in fl.readlines()]
    lines_zip = zip(*files)

    for linesb in lines_zip:
      lines = [l.decode('utf-8') for l in linesb]
      ex = lines[0]
      if ex.startswith('##'):
        sys.stdout.write(ex)
      else:
        sys.stdout.write(ex.strip())
        if len(lines) > 1:
          sys.stdout.write('\t')
          out_lines = [l.strip().split('\t', 9)[-1] for l in lines[1:]]
          sys.stdout.write("\t".join(out_lines))
        sys.stdout.write('\n')
    CODE
    tabix ~{prefix}.~{contig}.joined.vcf.gz
  >>>

  output {
    File out = "~{prefix}.~{contig}.joined.vcf.gz"
    File out_index = "~{prefix}.~{contig}.joined.vcf.gz.tbi"
  }
}

# Add in max CN state to multiallelics
task FixMultiallelicRecords {
  input {
    File joined_vcf
    Array[File] batch_contig_vcfs
    String contig
    String prefix
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  Float input_size = size(joined_vcf, "GiB") * 2 + size(batch_contig_vcfs, "GiB")
  Float input_size_fraction = 2.0
  Float base_disk_gb = 10.0
  RuntimeAttr runtime_default = object {
                                  mem_gb: 3.75,
                                  disk_gb: ceil(base_disk_gb + input_size * input_size_fraction),
                                  cpu_cores: 1,
                                  preemptible_tries: 3,
                                  max_retries: 1,
                                  boot_disk_gb: 10
                                }
  RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])

  runtime {
    memory: select_first([runtime_override.mem_gb, runtime_default.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_override.disk_gb, runtime_default.disk_gb]) + " HDD"
    cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
    preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
    maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
    docker: sv_pipeline_docker
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
  }

  command <<<
    set -euxo pipefail
    /opt/sv-pipeline/04_variant_resolution/scripts/make_concordant_multiallelic_alts.py \
      ~{joined_vcf} \
      ~{write_lines(batch_contig_vcfs)} \
      ~{prefix}.~{contig}.fixed_multiallelics.vcf.gz
    tabix ~{prefix}.~{contig}.fixed_multiallelics.vcf.gz
  >>>

  output {
    File out = "~{prefix}.~{contig}.fixed_multiallelics.vcf.gz"
    File out_index = "~{prefix}.~{contig}.fixed_multiallelics.vcf.gz.tbi"
  }
}

# Convert EV field from String to Integer
task FixEvidenceTags {
  input {
    File vcf
    String contig
    String prefix
    String sv_base_mini_docker
    RuntimeAttr? runtime_attr_override
  }

  Float input_size = size(vcf, "GiB")
  Float input_size_ratio = 2.0
  Float base_disk_gb = 10.0
  RuntimeAttr runtime_default = object {
                                  mem_gb: 1.0,
                                  disk_gb: ceil(base_disk_gb + input_size * input_size_ratio),
                                  cpu_cores: 1,
                                  preemptible_tries: 3,
                                  max_retries: 1,
                                  boot_disk_gb: 10
                                }
  RuntimeAttr runtime_override = select_first([runtime_attr_override, runtime_default])

  runtime {
    memory: select_first([runtime_override.mem_gb, runtime_default.mem_gb]) + " GiB"
    disks: "local-disk " + select_first([runtime_override.disk_gb, runtime_default.disk_gb]) + " HDD"
    cpu: select_first([runtime_override.cpu_cores, runtime_default.cpu_cores])
    preemptible: select_first([runtime_override.preemptible_tries, runtime_default.preemptible_tries])
    maxRetries: select_first([runtime_override.max_retries, runtime_default.max_retries])
    docker: sv_base_mini_docker
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
  }

  command <<<
    set -euxo pipefail
    zcat ~{vcf} \
      | sed -e 's/:RD,PE,SR/:7/g' \
      | sed -e 's/:PE,SR/:6/g' \
      | sed -e 's/:RD,SR/:5/g' \
      | sed -e 's/:RD,PE/:3/g' \
      | sed -e 's/:PE\t/:2\t/g' -e 's/:SR\t/:4\t/g' -e 's/:RD\t/:1\t/g' \
      | sed -e 's/ID=EV,Number=.,Type=String/ID=EV,Number=1,Type=Integer/g' \
      | bgzip \
      > ~{prefix}.~{contig}.unclustered.vcf.gz
    tabix ~{prefix}.~{contig}.unclustered.vcf.gz
  >>>

  output {
    File out = "~{prefix}.~{contig}.unclustered.vcf.gz"
    File out_index = "~{prefix}.~{contig}.unclustered.vcf.gz.tbi"
  }
}
