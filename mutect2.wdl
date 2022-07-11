version 1.0

workflow Mutect2_base{

      
# Mutect2 inputs
    input{
      File? intervals
      File ref_fasta
      File ref_fai
      File ref_dict
      File tumor_reads
      File tumor_reads_index
      File? normal_reads
      File? normal_reads_index
      File pon
      File pon_idx
      Int scatter_count = 50
      File gnomad
      File gnomad_idx

      Boolean compress

      String gatk_docker
      File? gatk_override
    }

    String output_basename = basename(basename(tumor_reads, ".bam"),".cram")  #hacky way to strip either .bam or .cram
    String unfiltered_name = output_basename + "-unfiltered"
    String filtered_name = output_basename + "-filtered"


    call SplitIntervals {
        input:
            intervals = intervals,
            ref_fasta = ref_fasta,
            ref_fai = ref_fai,
            ref_dict = ref_dict,
            scatter_count = scatter_count,
            gatk_docker = gatk_docker
    }

    scatter (subintervals in SplitIntervals.interval_files ) {
        call M2 {
            input:
                intervals = subintervals,
                ref_fasta = ref_fasta,
                ref_fai = ref_fai,
                ref_dict = ref_dict,
                tumor_reads = tumor_reads,
                tumor_reads_index = tumor_reads_index,
                normal_reads = normal_reads,
                normal_reads_index = normal_reads_index,
                pon = pon,
                pon_idx = pon_idx,
                gnomad = gnomad,
                gnomad_idx = gnomad_idx,
                compress = compress,
                gatk_docker = gatk_docker
        }
    }


    call MergeVCFs {
        input:
            input_vcfs = M2.unfiltered_vcf,
            input_vcf_indices = M2.unfiltered_vcf_idx,
            output_name = unfiltered_name,
            compress = compress,
            gatk_docker = gatk_docker
    }

    call MergeStats { 
            input: stats = M2.stats,
                   gatk_docker = gatk_docker
                   }

    call Filter {
        input:
            ref_fasta = ref_fasta,
            ref_fai = ref_fai,
            ref_dict = ref_dict,
            intervals = intervals,
            unfiltered_vcf = MergeVCFs.merged_vcf,
            unfiltered_vcf_idx = MergeVCFs.merged_vcf_idx,
            output_name = filtered_name,
            compress = compress,
            mutect_stats = MergeStats.merged_stats,
            gatk_docker = gatk_docker
    }

    output {


        File mutect_stats = MergeStats.merged_stats
        File merged_vcf = MergeVCFs.merged_vcf
        File Filtered_vcf = Filter.filtered_vcf
    }
}

task SplitIntervals {
    input{
      File? intervals
      File ref_fasta
      File ref_fai
      File ref_dict
      Int scatter_count
      String gatk_docker


      File? gatk_override
    }

    command {
        set -e
        export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk_override}

        mkdir interval-files
        gatk --java-options "-Xmx8g" SplitIntervals -R ~{ref_fasta} ~{"-L " + intervals} -scatter ~{scatter_count} -O interval-files
        cp interval-files/*.interval_list .
    }

    runtime {
          docker: gatk_docker
          
    }

    output {
        Array[File] interval_files = glob("*.interval_list")
    }
}

task M2 {
    input{
   File? intervals
      File ref_fasta
      File ref_fai
      File ref_dict
      File tumor_reads
      File tumor_reads_index
      File? normal_reads
      File? normal_reads_index
      File? pon
      File? pon_idx
      File? gnomad
      File? gnomad_idx

      Boolean compress


      File? gatk_override

      String gatk_docker
    }

    String output_vcf = "output" + if compress then ".vcf.gz" else ".vcf"
    String output_vcf_idx = output_vcf + if compress then ".tbi" else ".idx"

    String output_stats = output_vcf + ".stats"


    # Mem is in units of GB but our command and memory runtime values are in MB

      command <<<
        set -e

        export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk_override}

        # We need to create these files regardless, even if they stay empty
        touch bamout.bam
        touch f1r2.tar.gz
        echo "" > normal_name.txt

        gatk --java-options "-Xmx8g" GetSampleName -R ~{ref_fasta} -I ~{tumor_reads} -O tumor_name.txt -encode 
        tumor_command_line="-I ~{tumor_reads} -tumor `cat tumor_name.txt`"

        if [[ ! -z "~{normal_reads}" ]]; then
            gatk --java-options "-Xmx8g" GetSampleName -R ~{ref_fasta} -I ~{normal_reads} -O normal_name.txt -encode 
            normal_command_line="-I ~{normal_reads} -normal `cat normal_name.txt`"
        fi

        gatk --java-options "-Xmx8g" Mutect2 -R ~{ref_fasta} $tumor_command_line $normal_command_line ~{"--germline-resource " + gnomad} ~{"-pon " + pon} ~{"-L " + intervals} -O "~{output_vcf}" --bam-output bamout.bam

        m2_exit_code=$?


        # the script only fails if Mutect2 itself fails
        exit $m2_exit_code
    >>>

    runtime {
          docker: gatk_docker
	  cpu: 1

          }

    output {
        File unfiltered_vcf = "~{output_vcf}"
        File unfiltered_vcf_idx = "~{output_vcf_idx}"
        File output_bamOut = "bamout.bam"
        String tumor_sample = read_string("tumor_name.txt")
        String normal_sample = read_string("normal_name.txt")
        File stats = "~{output_stats}"
        File f1r2_counts = "f1r2.tar.gz"
        Array[File] tumor_pileups = glob("*tumor-pileups.table")
        Array[File] normal_pileups = glob("*normal-pileups.table")
    }
}

task MergeVCFs {
    input{
      Array[File] input_vcfs
      Array[File] input_vcf_indices
      String output_name
      Boolean compress
      String gatk_docker
      File? gatk_override
    }    

    String output_vcf = output_name + if compress then ".vcf.gz" else ".vcf"
    String output_vcf_idx = output_vcf + if compress then ".tbi" else ".idx"

    # using MergeVcfs instead of GatherVcfs so we can create indices
    # WARNING 2015-10-28 15:01:48 GatherVcfs  Index creation not currently supported when gathering block compressed VCFs.

    command {
        set -e
        export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk_override}
        gatk --java-options "-Xmx8g" MergeVcfs -I ~{sep=' -I ' input_vcfs} -O ~{output_vcf}
    }

     runtime {
          docker: gatk_docker

          }

    output {
        File merged_vcf = "~{output_vcf}"
        File merged_vcf_idx = "~{output_vcf_idx}"
    }
}

task MergeStats {
    input{
      Array[File]+ stats
      String gatk_docker
      File? gatk_override
    }

    command {
        set -e
        export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk_override}


        gatk --java-options "-Xmx8g" MergeMutectStats -stats ~{sep=" -stats " stats} -O merged.stats
    }

       runtime {
          docker: gatk_docker

          }

    output {
        File merged_stats = "merged.stats"
    }
}

task Filter {
    input{
      File? intervals
      File ref_fasta
      File ref_fai
      File ref_dict
      File unfiltered_vcf
      File unfiltered_vcf_idx
      String output_name
      Boolean compress
      File? mutect_stats
      String gatk_docker
      File? gatk_override
    }


    String output_vcf = output_name + if compress then ".vcf.gz" else ".vcf"
    String output_vcf_idx = output_vcf + if compress then ".tbi" else ".idx"



    command {
        set -e

        export GATK_LOCAL_JAR=~{default="/root/gatk.jar" gatk_override}

        gatk --java-options "-Xmx8g" FilterMutectCalls -V ~{unfiltered_vcf} -R ~{ref_fasta} -O ~{output_vcf} ~{"-stats " + mutect_stats} --filtering-stats filtering.stats
    }
    runtime{
        docker: gatk_docker

          }

    output {
        File filtered_vcf = "~{output_vcf}"
        File filtered_vcf_idx = "~{output_vcf_idx}"
        File filtering_stats = "filtering.stats"
    }
}


