version 1.0

import "GvsUtils.wdl" as Utils


workflow GvsCalculatePrecisionAndSensitivity {
  input {
    String? git_branch_or_tag
    File input_vcf_fofn
    String output_basename

    Array[String] chromosomes = ["chr20"]

    Array[String] sample_names
    Array[File] truth_vcfs
    Array[File] truth_vcf_indices
    Array[File] truth_beds

    File ref_fasta

    String? basic_docker
    String? variants_docker
    String? gatk_docker
    String? real_time_genomics_docker
    String? gotc_imputation_docker
  }

  parameter_meta {
    input_vcf_fofn: "A FOFN (file of file names) of input VCFs."
    output_basename: "The base name for the output files generated by the pipeline."
    chromosomes: "The chromosome(s) on which to analyze precision and sensitivity. The default value for this is `['chr20']`."
    sample_names: "A list of the sample names that are controls and that will be used for the analysis. For every element on the list of sample names there must be a corresponding element on the list of `truth_vcfs`, `truth_vcf_indices`, and `truth_beds`."
    truth_vcfs: "A list of the VCFs that contain the truth data used for analyzing the samples in `sample_names`."
    truth_vcf_indices: "A list of the VCF indices for the truth data VCFs supplied above."
    truth_beds: "A list of the bed files for the truth data used for analyzing the samples in `sample_names`."
    ref_fasta: "The cloud path for the reference fasta sequence."
  }

  # Always call `GetToolVersions` to get the git hash for this run as this is a top-level-only WDL (i.e. there are
  # no calling WDLs that might supply `git_hash`).
  call Utils.GetToolVersions {
    input:
      git_branch_or_tag = git_branch_or_tag,
  }

  String effective_basic_docker = select_first([basic_docker, GetToolVersions.basic_docker])
  String effective_variants_docker = select_first([variants_docker, GetToolVersions.variants_docker])
  String effective_gatk_docker = select_first([gatk_docker, GetToolVersions.gatk_docker])
  String effective_real_time_genomics_docker = select_first([real_time_genomics_docker, GetToolVersions.real_time_genomics_docker])
  String effective_gotc_imputation_docker = select_first([gotc_imputation_docker, GetToolVersions.gotc_imputation_docker])

  if ((length(sample_names) != length(truth_vcfs)) || (length(sample_names) != length(truth_vcf_indices)) || (length(sample_names) != length(truth_beds))) {
    call Utils.TerminateWorkflow {
      input:
        message = "The inputs 'sample_names', 'truth_vcfs', 'truth_vcf_indices', and 'truth_beds' must all contain the same number of elements",
        basic_docker = effective_basic_docker,
    }
  }

  call CountInputVcfs {
    input:
      input_vcf_fofn = input_vcf_fofn,
      basic_docker = effective_basic_docker,
  }

  scatter(i in range(CountInputVcfs.num_vcfs)) {
    call IsVcfOnChromosomes {
      input:
        input_vcf_fofn = input_vcf_fofn,
        index = i,
        chromosomes = chromosomes,
        variants_docker = effective_variants_docker,
    }
  }

  call GatherVcfs {
    input:
      input_vcfs = flatten(IsVcfOnChromosomes.output_vcf),
      output_basename = output_basename,
      gatk_docker = effective_gatk_docker,
  }

  scatter(i in range(length(sample_names))) {
    String sample_name = sample_names[i]
    String output_sample_basename = output_basename + "." + sample_name

    call SelectVariants {
      input:
        input_vcf = GatherVcfs.output_vcf,
        input_vcf_index = GatherVcfs.output_vcf_index,
        chromosomes = chromosomes,
        sample_name = sample_name,
        output_basename = output_sample_basename,
        gatk_docker = effective_gatk_docker,
    }

    call Add_AS_MAX_VQS_SCORE_ToVcf {
      input:
        input_vcf = SelectVariants.output_vcf,
        output_basename = output_sample_basename + ".maxas",
        variants_docker = effective_variants_docker,
    }

    call IsVQSRLite {
      input:
        input_vcf = Add_AS_MAX_VQS_SCORE_ToVcf.output_vcf,
        basic_docker = effective_basic_docker,
    }

    call BgzipAndTabix {
      input:
        input_vcf = Add_AS_MAX_VQS_SCORE_ToVcf.output_vcf,
        output_basename = output_sample_basename + ".maxas",
        gotc_imputation_docker = effective_gotc_imputation_docker,
    }

    call EvaluateVcf as EvaluateVcfFiltered {
      input:
        input_vcf = BgzipAndTabix.output_vcf,
        input_vcf_index = BgzipAndTabix.output_vcf_index,
        truth_vcf = truth_vcfs[i],
        truth_vcf_index = truth_vcf_indices[i],
        truth_bed = truth_beds[i],
        chromosomes = chromosomes,
        output_basename = sample_name + "-bq_roc_filtered",
        is_vqsr_lite = IsVQSRLite.is_vqsr_lite,
        ref_fasta = ref_fasta,
        real_time_genomics_docker = effective_real_time_genomics_docker,
    }

    call EvaluateVcf as EvaluateVcfAll {
      input:
        input_vcf = BgzipAndTabix.output_vcf,
        input_vcf_index = BgzipAndTabix.output_vcf_index,
        truth_vcf = truth_vcfs[i],
        truth_vcf_index = truth_vcf_indices[i],
        truth_bed = truth_beds[i],
        chromosomes = chromosomes,
        all_records = true,
        output_basename = sample_name + "-bq_all",
        is_vqsr_lite = IsVQSRLite.is_vqsr_lite,
        ref_fasta = ref_fasta,
        real_time_genomics_docker = effective_real_time_genomics_docker,
    }
  }

  call CollateReports {
    input:
      all_reports = EvaluateVcfAll.report,
      filtered_reports = EvaluateVcfFiltered.report,
      basic_docker = effective_basic_docker,
  }

  output {
    File report = CollateReports.report
    Array[Array[File]] filtered_eval_outputs = EvaluateVcfFiltered.outputs
    Array[Array[File]] all_eval_outputs = EvaluateVcfAll.outputs
    String recorded_git_hash = GetToolVersions.git_hash
  }
}

task IsVcfOnChromosomes {
  input {
    File input_vcf_fofn
    Int index
    Array[String] chromosomes
    String variants_docker
  }

  command <<<
    # Prepend date, time and pwd to xtrace log entries.
    PS4='\D{+%F %T} \w $ '
    set -o errexit -o nounset -o pipefail -o xtrace

    # Get a one-based line number from the zero-based input index.
    line_number=$((~{index} + 1))
    input_vcf_path=$(sed -n "${line_number}p" ~{input_vcf_fofn})

    input_vcf=$(basename ${input_vcf_path})
    gcloud storage cp ${input_vcf_path} ${input_vcf}
    cat ${input_vcf} | gunzip | grep -v '^#' | cut -f 1 | sort | uniq > chrom.txt
    NL=$(cat chrom.txt | wc -l)
    if [ $NL -ne 1 ]; then
      echo "${input_vcf_path} has either no records or records on multiple chromosomes"
      exit 1
    fi

    mkdir output
    output_vcf_name=${input_vcf}
    touch output/${output_vcf_name}
    VCF_CHR=$(cat chrom.txt)
    chromosomes=( ~{sep=' ' chromosomes} )

    for i in "${chromosomes[@]}"
    do
      if [ $VCF_CHR = $i ]; then
        cp ${input_vcf} output/${output_vcf_name}
        echo "Including ${input_vcf_path} as it is on $i."
        break
      else
        touch output/${output_vcf_name}
        echo "Skipping ${input_vcf_path} as it is not on $i."
      fi
    done
  >>>

  runtime {
    docker: variants_docker
    disks: "local-disk 100 HDD"
    memory: "2 GiB"
    preemptible: 3
  }
  output {
    Array[File] output_vcf = glob("output/*")
  }
}

task GatherVcfs {
  input {
    Array[File] input_vcfs
    String output_basename

    String gatk_docker
    Int cpu = 1
    Int memory_mb = 7500
    Int disk_size_gb = ceil(3*size(input_vcfs, "GiB")) + 500
  }
  Int command_mem = memory_mb - 1000
  Int max_heap = memory_mb - 500

  command <<<
    # Prepend date, time and pwd to xtrace log entries.
    PS4='\D{+%F %T} \w $ '
    set -o errexit -o nounset -o pipefail -o xtrace

    CHR_VCFS_ARG=""
    for file in ~{sep=' ' input_vcfs}
    do
      if [ -s $file ]; then
        CHR_VCFS_ARG+=" --INPUT $file "
      fi
    done
    echo $CHR_VCFS_ARG

    # --REORDER_INPUT_BY_FIRST_VARIANT means that the vcfs supplied here need not be ordered by location.
    gatk --java-options "-Xms~{command_mem}m -Xmx~{max_heap}m" \
      GatherVcfs \
        --REORDER_INPUT_BY_FIRST_VARIANT \
        $CHR_VCFS_ARG \
        --OUTPUT ~{output_basename}.vcf.gz

    tabix ~{output_basename}.vcf.gz
  >>>

  runtime {
    docker: gatk_docker
    cpu: cpu
    memory: "${memory_mb} MiB"
    disks: "local-disk ${disk_size_gb} HDD"
    bootDiskSizeGb: 15
    preemptible: 1
  }

  output {
    File output_vcf = "~{output_basename}.vcf.gz"
    File output_vcf_index = "~{output_basename}.vcf.gz.tbi"
  }
}

task SelectVariants {
  input {
    File input_vcf
    File input_vcf_index
    Array[String] chromosomes
    String sample_name

    String output_basename

    String gatk_docker
    Int cpu = 1
    Int memory_mb = 7500
    Int disk_size_gb = ceil(2*size(input_vcf, "GiB")) + 500
  }
  Int command_mem = memory_mb - 1000
  Int max_heap = memory_mb - 500

  command <<<
    gatk --java-options "-Xms~{command_mem}m -Xmx~{max_heap}m" \
      SelectVariants \
        -L ~{sep=' -L ' chromosomes} \
        -V ~{input_vcf} \
        --sample-name ~{sample_name} \
        --select-type-to-exclude NO_VARIATION \
        -O ~{output_basename}.vcf.gz
  >>>

  runtime {
    docker: gatk_docker
    cpu: cpu
    memory: "${memory_mb} MiB"
    disks: "local-disk ${disk_size_gb} HDD"
    bootDiskSizeGb: 15
    preemptible: 1
  }

  output {
    File output_vcf = "~{output_basename}.vcf.gz"
    File output_vcf_index = "~{output_basename}.vcf.gz.tbi"
  }
}

task Add_AS_MAX_VQS_SCORE_ToVcf {
  input {
    File input_vcf
    String output_basename

    String variants_docker
    Int cpu = 1
    Int memory_mb = 3500
    Int disk_size_gb = ceil(2*size(input_vcf, "GiB")) + 500
  }

  command <<<
    set -e

    python3 /app/add_max_as_vqs_score.py ~{input_vcf} > ~{output_basename}.vcf
  >>>
  runtime {
    docker: variants_docker
    cpu: cpu
    memory: "${memory_mb} MiB"
    disks: "local-disk ${disk_size_gb} HDD"
  }

  output {
    File output_vcf = "~{output_basename}.vcf"
  }
}

task IsVQSRLite {
  input {
    File input_vcf
    String basic_docker
  }

  String is_vqsr_lite_file = "is_vqsr_lite_file.txt"

  command {
    set +e

    # See if there are any non-header lines that contain the string 'AS_VQS_SENS'. If so, grep will return 0 else 1
    grep -v '^#' ~{input_vcf} | grep AS_VQS_SENS > /dev/null
    if [[ $? -eq 0 ]]; then
      echo "true" > ~{is_vqsr_lite_file}
    else
      echo "false" > ~{is_vqsr_lite_file}
    fi
  }

  runtime {
    docker: basic_docker
    disks: "local-disk 10 HDD"
    memory: "2 GiB"
    preemptible: 3
  }
  output {
    Boolean is_vqsr_lite = read_boolean(is_vqsr_lite_file)
  }
}

task BgzipAndTabix {
  input {
    File input_vcf
    String output_basename

    String gotc_imputation_docker
    Int cpu = 1
    Int memory_mb = 3500
    Int disk_size_gb = ceil(3 * size(input_vcf, "GiB")) + 500
  }

  command {
    # note that bgzip has an option (-i) to index the bgzipped output, but this file is not a tabix file
    # note also that we use '-c' so that bgzip doesn't create the bgzipped file in place, rather it's in a location
    # where it's easy to output from the task.
    bgzip -c ~{input_vcf} > ~{output_basename}.vcf.gz
    tabix ~{output_basename}.vcf.gz
  }
  runtime {
    docker: gotc_imputation_docker
    cpu: cpu
    memory: "${memory_mb} MiB"
    disks: "local-disk ${disk_size_gb} HDD"

  }
  output {
    File output_vcf = "~{output_basename}.vcf.gz"
    File output_vcf_index = "~{output_basename}.vcf.gz.tbi"
  }
}

task EvaluateVcf {
  input {
    File input_vcf
    File input_vcf_index
    File truth_vcf
    File truth_vcf_index
    File truth_bed

    Boolean all_records = false
    Array[String] chromosomes

    File ref_fasta

    String output_basename

    Boolean is_vqsr_lite

    String real_time_genomics_docker
    Int cpu = 1
    Int memory_mb = 3500
    Int disk_size_gb = ceil(2 * size(ref_fasta, "GiB")) + 500
  }

  String max_score_field_tag = if (is_vqsr_lite == true) then 'MAX_AS_VQS_SENS' else 'MAX_AS_VQSLOD'

  command <<<
    set -e -o pipefail

    chromosomes=( ~{sep=' ' chromosomes} )

    echo "Creating .bed file to control which chromosomes should be evaluated."
    for i in "${chromosomes[@]}"
    do
      echo "$i	0	300000000" >> chromosomes.to.eval.txt
    done

    rtg format --output human_REF_SDF ~{ref_fasta}

    rtg vcfeval \
      --bed-regions chromosomes.to.eval.txt \
      ~{if all_records then "--all-records" else ""} \
      --roc-subset snp,indel \
      --vcf-score-field=INFO.~{max_score_field_tag} \
      ~{if is_vqsr_lite then "--sort-order ascending" else "--sort-order descending"} \
      -t human_REF_SDF \
      -b ~{truth_vcf} \
      -e ~{truth_bed}\
      -c ~{input_vcf} \
      -o ~{output_basename}

    # Touch a file with the name of the sample in that directory, so that it's identifiable among the globbed outputs.
    touch ~{output_basename}/~{output_basename}

    touch report.txt
    for type in "snp" "indel"
      do
        d=$(cat ~{output_basename}/${type}_roc.tsv.gz | gunzip | tail -1 | cut -f3,5,6,7)
        echo -e "~{output_basename}\t$type\t$d" >> report.txt
    done
  >>>
  runtime {
    docker: real_time_genomics_docker
    cpu: cpu
    memory: "${memory_mb} MiB"
    disks: "local-disk ${disk_size_gb} HDD"

  }
  output {
    File coverage = "chromosomes.to.eval.txt"
    File report = "report.txt"
    Array[File] outputs = glob("~{output_basename}/*")
  }
}

task CollateReports {
  input {
    Array[File] all_reports
    Array[File] filtered_reports
    String basic_docker
  }

  command {
    set -e -o pipefail

    echo "sample  type  FPs FNs precision sensitivity"
    while read -r a;
    do
      cat $a
    done < ~{write_lines(all_reports)}

    while read -r a;
    do
      cat $a
    done < ~{write_lines(filtered_reports)}
  }

  runtime {
    docker: basic_docker
    disks: "local-disk 100 HDD"
    memory: "2 GiB"
    preemptible: 3
  }
  output {
    File report = stdout()
  }
}


task CountInputVcfs {
  input {
    File input_vcf_fofn
    String basic_docker
  }
  command <<<
    wc -l < ~{input_vcf_fofn} > num_vcfs.txt
  >>>
  output {
    Int num_vcfs = read_int("num_vcfs.txt")
  }
  runtime {
    docker: basic_docker
  }
}