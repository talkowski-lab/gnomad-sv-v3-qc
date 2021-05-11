version 1.0

import "Structs.wdl"

workflow DropRedundantCNVs {
  input {
    File vcf
    String contig
    String sv_pipeline_docker
  }

  call DropRedundantCNVs_1 {
    input:
      vcf=vcf,
      sv_pipeline_docker=sv_pipeline_docker
  }

  call DropRedundantCNVs_2 {
    input:
      intervals_preclustered_bed=DropRedundantCNVs_1.intervals_preclustered_bed,
      sv_pipeline_docker=sv_pipeline_docker
  }

  call DropRedundantCNVs_3 {
    input:
      intervals_preclustered_subset_bed=DropRedundantCNVs_2.intervals_preclustered_subset_bed,
      sv_pipeline_docker=sv_pipeline_docker
  }

  call DropRedundantCNVs_4 {
    input:
      intervals_preclustered_subset_melted_bed=DropRedundantCNVs_3.intervals_preclustered_subset_melted_bed,
      sv_pipeline_docker=sv_pipeline_docker
  }

  call DropRedundantCNVs_5 {
    input:
      intervals_clustered_bed=DropRedundantCNVs_4.intervals_clustered_bed,
      intervals_preclustered_bed=DropRedundantCNVs_1.intervals_preclustered_bed,
      sv_pipeline_docker=sv_pipeline_docker
  }

  call DropRedundantCNVs_6 {
    input:
      intervals_preclustered_bed=DropRedundantCNVs_1.intervals_preclustered_bed,
      sv_pipeline_docker=sv_pipeline_docker
  }

  call DropRedundantCNVs_7 {
    input:
      step2_intervals_preclustered_subset_bed=DropRedundantCNVs_6.step2_intervals_preclustered_subset_bed,
      sv_pipeline_docker=sv_pipeline_docker
  }

  call DropRedundantCNVs_8 {
    input:
      vcf=vcf,
      vids_to_remove_list_1=DropRedundantCNVs_5.vids_to_remove_list_1,
      intervals_preclustered_bed=DropRedundantCNVs_1.intervals_preclustered_bed,
      step2_variants_to_resolve_list=DropRedundantCNVs_7.step2_variants_to_resolve_list,
      sv_pipeline_docker=sv_pipeline_docker
  }

  call DropRedundantCNVs_9 {
    input:
      vcf=vcf,
      records_to_add_vcf=DropRedundantCNVs_8.records_to_add_vcf,
      vids_to_remove_list_2=DropRedundantCNVs_8.vids_to_remove_list_2,
      contig=contig,
      sv_pipeline_docker=sv_pipeline_docker
  }

  output {
    File cleaned_vcf_shard = DropRedundantCNVs_9.cleaned_vcf_shard
  }
}

task DropRedundantCNVs_1 {
  input {
    File vcf
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  Float input_size = size(vcf, "GB")
  RuntimeAttr runtime_default = object {
                                  mem_gb: 7.5,
                                  disk_gb: ceil(10.0 + input_size * 2),
                                  cpu_cores: 1,
                                  preemptible_tries: 0,
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
    docker: sv_pipeline_docker
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
  }

  command <<<
    set -euxo pipefail

    ###PREP FILES
    #Convert full VCF to BED intervals
    #Ignore CPX events with UNRESOLVED filter status
    svtk vcf2bed --split-cpx --info SVTYPE \
      <(bcftools view -e 'INFO/SVTYPE == "CPX" && FILTER == "UNRESOLVED"' ~{vcf}) - \
      | grep -e '^#\|DEL\|DUP\|CNV\|CPX' \
      | awk -v OFS="\t" '{ if ($5=="CN0") print $1, $2, $3, $4, "DEL", $5"\n"$1, $2, $3, $4, "DUP", $5; \
        else if ($5=="DEL" || $5=="DUP") print $1, $2, $3, $4, $6, $5 }' \
      | sort -Vk1,1 -k2,2n -k3,3n -k4,4V \
      | bgzip -c \
      > intervals.preclustered.bed.gz
  >>>

  output {
    File intervals_preclustered_bed = "intervals.preclustered.bed.gz"
  }
}

task DropRedundantCNVs_2 {
  input {
    File intervals_preclustered_bed
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  Float input_size = size(intervals_preclustered_bed, "GB")
  RuntimeAttr runtime_default = object {
                                  mem_gb: 7.5,
                                  disk_gb: ceil(10.0 + input_size * 5.0),
                                  cpu_cores: 1,
                                  preemptible_tries: 0,
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
    docker: sv_pipeline_docker
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
  }

  command <<<
    set -euxo pipefail

    ###REMOVE CNVS REDUNDANT WITH COMPLEX EVENTS
    #Subset to only variants that share some overlap (at least 10% recip) with at least one CPX variant
    bedtools intersect -wa -r -f 0.1 \
      -a ~{intervals_preclustered_bed} \
      -b <( zcat ~{intervals_preclustered_bed} | fgrep "CPX" ) \
      | sort -Vk1,1 -k2,2n -k3,3n -k4,4V \
      | uniq \
      | bgzip -c \
      > intervals.preclustered.subset.bed.gz
  >>>

  output {
    File intervals_preclustered_subset_bed = "intervals.preclustered.subset.bed.gz"
  }
}


task DropRedundantCNVs_3 {
  input {
    File intervals_preclustered_subset_bed
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  Float input_size = size(intervals_preclustered_subset_bed, "GB")
  RuntimeAttr runtime_default = object {
                                  mem_gb: 7.5,
                                  disk_gb: ceil(10.0 + input_size * 2.0),
                                  cpu_cores: 1,
                                  preemptible_tries: 0,
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
    docker: sv_pipeline_docker
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
  }

  command <<<
    set -euxo pipefail

    #Melt subsetted variants
    while read chr start end VID samples CNV; do
      echo -e "${samples}" \
      | sed 's/,/\n/g' \
      | awk -v OFS="\t" -v chr=${chr} -v start=${start} -v end=${end} -v VID=${VID} -v CNV=${CNV} \
        '{ print chr, start, end, VID, $1, CNV }'
    done < <( zcat ~{intervals_preclustered_subset_bed} ) \
      | bgzip -c \
      > intervals.preclustered.subset.melted.bed.gz
  >>>

  output {
    File intervals_preclustered_subset_melted_bed = "intervals.preclustered.subset.melted.bed.gz"
  }
}


task DropRedundantCNVs_4 {
  input {
    File intervals_preclustered_subset_melted_bed
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  Float input_size = size(intervals_preclustered_subset_melted_bed, "GB")
  RuntimeAttr runtime_default = object {
                                  mem_gb: 7.5,
                                  disk_gb: ceil(10.0 + input_size * 2.0),
                                  cpu_cores: 1,
                                  preemptible_tries: 0,
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
    docker: sv_pipeline_docker
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
  }

  command <<<
    set -euxo pipefail

    #Cluster BED intervals (50% RO)
    svtk bedcluster -f 0.5 \
      ~{intervals_preclustered_subset_melted_bed} - \
      | bgzip -c \
      > intervals.clustered.bed.gz
  >>>

  output {
    File intervals_clustered_bed = "intervals.clustered.bed.gz"
  }
}


task DropRedundantCNVs_5 {
  input {
    File intervals_clustered_bed
    File intervals_preclustered_bed
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  Float input_size = size([intervals_clustered_bed, intervals_preclustered_bed], "GB")
  RuntimeAttr runtime_default = object {
                                  mem_gb: 7.5,
                                  disk_gb: ceil(10.0 + input_size * 2.0),
                                  cpu_cores: 1,
                                  preemptible_tries: 0,
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
    docker: sv_pipeline_docker
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
  }

  command <<<
    set -euo pipefail

    #Get list of all variants that cluster with a complex variant,
    # evaluate sample overlap from original intervals file,
    # and, if overlap >50%, write that ID to be stripped from the output VCF
    while read VIDs; do
      #Get nonredundant list of sample IDs involved in any clustered variant
      echo -e "${VIDs}" | sed 's/,/\n/g' \
        | fgrep -wf - <( zcat ~{intervals_preclustered_bed} ) \
        | cut -f5 | sort | uniq \
        > nonredundant_samples.list

      #Iterate over VIDs and print non-CPX VID if sample overlap >50%
      while read VID samples; do
        #Get list of samples in variant
        echo -e "${samples}" | sed 's/,/\n/g' \
          | sort | uniq > query_samples.list
          nsamp=$( cat query_samples.list | wc -l )

        #Compare
        frac=$( fgrep -wf query_samples.list \
          nonredundant_samples.list | wc -l \
          | awk -v nsamp=${nsamp} '{ print 100*($1/nsamp) }' \
          | cut -f1 -d\. )
        if [ ${frac} -ge 50 ]; then
          echo "${VID}"
        fi

        #Clean up
        rm query_samples.list

      done < <( echo -e "${VIDs}" | sed 's/,/\n/g' \
        | fgrep -wf - <( zcat ~{intervals_preclustered_bed} ) \
        | cut -f4,5 | sort | uniq | fgrep -v "CPX" )

      #Clean up
      rm nonredundant_samples.list

    done < <( zcat ~{intervals_clustered_bed} \
      | cut -f7 | fgrep "CPX" | grep -e "DEL\|DUP" ) \
      | sort -V | uniq \
      > VIDs_to_remove.list
  >>>

  output {
    File vids_to_remove_list_1= "VIDs_to_remove.list"
  }
}



task DropRedundantCNVs_6 {
  input {
    File intervals_preclustered_bed
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  Float input_size = size(intervals_preclustered_bed, "GB")
  RuntimeAttr runtime_default = object {
                                  mem_gb: 7.5,
                                  disk_gb: ceil(10.0 + input_size * 2.0),
                                  cpu_cores: 1,
                                  preemptible_tries: 0,
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
    docker: sv_pipeline_docker
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
  }

  command <<<
    set -euxo pipefail

    ###FIND REMAINING REDUNDANT CNVS WITH STRONG (80%) OVERLAP IN SAMPLES AND SIZE
    #Find CNV intervals that have 80% reciprocal overlap
    bedtools intersect -wa -wb -r -f 0.8 \
      -a ~{intervals_preclustered_bed} \
      -b ~{intervals_preclustered_bed} \
      | awk -v FS="\t" '{ if ($4!=$10 && $6==$12) print $0 }' \
      | awk -v OFS="\t" '$4 ~ /DEL|DUP/ { print $0 }' \
      | awk -v OFS="\t" '$10 ~ /DEL|DUP/ { print $0 }' \
      | bgzip -c \
      > step2.intervals.preclustered.subset.bed.gz
  >>>

  output {
    File step2_intervals_preclustered_subset_bed = "step2.intervals.preclustered.subset.bed.gz"
  }
}


task DropRedundantCNVs_7 {
  input {
    File step2_intervals_preclustered_subset_bed
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  Float input_size = size(step2_intervals_preclustered_subset_bed, "GB")
  RuntimeAttr runtime_default = object {
                                  mem_gb: 7.5,
                                  disk_gb: ceil(10.0 + input_size * 2.0),
                                  cpu_cores: 1,
                                  preemptible_tries: 0,
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
    docker: sv_pipeline_docker
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
  }

  command <<<
    set -euo pipefail

    #Determine which events share 80% sample overlap
    while read VIDa sa VIDb sb; do
      na=$( echo -e "${sa}" | sed 's/,/\n/g' | sort | uniq | wc -l )
      nb=$( echo -e "${sb}" | sed 's/,/\n/g' | sort | uniq | wc -l )
      denom=$( echo -e "${sa},${sb}" | sed 's/,/\n/g' | sort | uniq | wc -l )
      numer=$( echo -e "${sa}" | sed 's/,/\n/g' | { fgrep -wf - <( echo -e "${sb}" | sed 's/,/\n/g' ) || true; } \
        | sort | uniq | wc -l )

      if [ ${denom} -gt 0 ]; then
        ovr=$(( 100 * ${numer} / ${denom} ))
      fi
      if [ -z ${ovr} ]; then
        ovr=0
      fi
      if [ ${ovr} -ge 80 ]; then
        echo -e "${VIDa}\n${VIDb}" \
          | sort | uniq | paste -s -d,
      fi
    done < <( zcat ~{step2_intervals_preclustered_subset_bed} \
      | cut -f4,5,10,11 ) \
      | sort | uniq \
      > step2.variants_to_resolve.list
  >>>

  output {
    File step2_variants_to_resolve_list = "step2.variants_to_resolve.list"
  }
}


task DropRedundantCNVs_8 {
  input {
    File vcf
    File vids_to_remove_list_1
    File intervals_preclustered_bed
    File step2_variants_to_resolve_list
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  Float input_size = size([vcf, intervals_preclustered_bed, step2_variants_to_resolve_list], "GB")
  RuntimeAttr runtime_default = object {
                                  mem_gb: 7.5,
                                  disk_gb: ceil(10.0 + input_size * 5.0),
                                  cpu_cores: 1,
                                  preemptible_tries: 0,
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
    docker: sv_pipeline_docker
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
  }

  command <<<
    set -euxo pipefail

    python3 <<CODE
    import sys
    import pysam
    import gzip

    sys.stderr.write("Reading step2...\n")
    with open("~{step2_variants_to_resolve_list}") as f:
        vids_sets_to_resolve = [set(line.strip().split(',')) for line in f.readlines()]
        vids_list = sorted(list(set([x for y in vids_sets_to_resolve for x in y])))

    sys.stderr.write("Reading vids to remove...\n")
    with open("~{vids_to_remove_list_1}") as f:
        vids_to_remove = set([line.strip() for line in f.readlines()])

    sys.stderr.write("Reading preclustered intervals...\n")
    with gzip.open("~{intervals_preclustered_bed}") as f:
        intervals = {}
        for lineb in f:
            tokens = lineb.decode('utf-8').strip().split('\t')
            vid = tokens[3]
            intervals[vid] = tokens

    sys.stderr.write("Finding partners...\n")
    partners = {}
    all_partners = set([])
    for vid in vids_list:
        # get all other variants from clusters containing this variant
        partners[vid] = set([p for vset in vids_sets_to_resolve if vid in vset for p in vset])
        all_partners.update(partners[vid])

    vids_to_remove.update(all_partners)
    with open("vids_to_remove_2.list", 'w') as f:
        f.writelines(sorted([v+"\n" for v in vids_to_remove]))

    sys.stderr.write("Scanning vcf...\n")
    vcf = pysam.VariantFile("~{vcf}")
    records = {r.id: r for r in vcf if r.id in all_partners}

    def count_gts(record):
        result = [0, 0, 0]
        num_samples = len(record.samples)
        for g in [record.samples[i]['GT'] for i in range(num_samples)]:
            if g == (0, 0):
                result[1] += 1
            elif g == (None, None):
                result[2] += 1
            else:
                result[0] += 1
        return result

    def get_best_score_vid(scores):
        return sorted(scores.items(), key=lambda x: x[1])[-1][0]

    sys.stderr.write("Generating records...\n")
    with open("records_to_add.vcf", 'w') as f:
        processed_vids = set([])
        for vid in vids_list:
            if vid in processed_vids:
                continue
            vid_partners = partners[vid]
            processed_vids.update(vid_partners)
            partner_intervals = [intervals[p] for p in vid_partners]
            most_samples_vid = sorted(partner_intervals, key=lambda x : len(x[4].split(',')))[-1][3]
            x = sorted(partner_intervals, key=lambda x : len(x[4].split(',')))
            best_genotype_vid = None
            best_non_ref = -1
            best_ref = -1
            scores = {p: count_gts(records[p]) for p in vid_partners}
            scores_non_ref = {p: scores[p][0] for p in vid_partners if scores[p][0] > 0}
            scores_ref = {p: scores[p][1] for p in vid_partners if scores[p][1] > 0}
            scores_no_call = {p: scores[p][2] for p in vid_partners if scores[p][2] > 0}
            if len(scores_non_ref) > 0:
                best_genotype_vid = get_best_score_vid(scores_non_ref)
            elif len(scores_ref) > 0:
                best_genotype_vid = get_best_score_vid(scores_ref)
            else:
                best_genotype_vid = get_best_score_vid(scores_no_call)
            s1 = str(records[most_samples_vid]).split('\t')[0:9]
            s2 = str(records[best_genotype_vid]).split('\t', 9)
            f.write("\t".join(s1) + "\t" + s2[9])
    CODE

    cat ~{vids_to_remove_list_1} vids_to_remove_2.list | sort | uniq > vids_to_remove.list
  >>>

  output {
    File records_to_add_vcf = "records_to_add.vcf"
    File vids_to_remove_list_2 = "vids_to_remove.list"
  }
}


task DropRedundantCNVs_9 {
  input {
    File vcf
    File records_to_add_vcf
    File vids_to_remove_list_2
    String contig
    String sv_pipeline_docker
    RuntimeAttr? runtime_attr_override
  }

  String outfile_name = contig + ".shard.no_CNV_redundancies.vcf.gz"

  Float input_size = size([], "GB")
  RuntimeAttr runtime_default = object {
                                  mem_gb: 7.5,
                                  disk_gb: ceil(10.0 + input_size * 2.0),
                                  cpu_cores: 1,
                                  preemptible_tries: 0,
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
    docker: sv_pipeline_docker
    bootDiskSizeGb: select_first([runtime_override.boot_disk_gb, runtime_default.boot_disk_gb])
  }

  command <<<
    set -euxo pipefail

    ###CLEAN UP FINAL OUTPUT
    zcat ~{vcf} \
      | { fgrep -wvf ~{vids_to_remove_list_2} || true; } \
      | cat - ~{records_to_add_vcf} \
      | vcf-sort \
      | bgzip -c \
      > ~{outfile_name}
  >>>

  output {
    File cleaned_vcf_shard = outfile_name
  }
}

