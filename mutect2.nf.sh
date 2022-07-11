workflow mutect2 {

    take:
        aligned_bams
        

    main:

        ch_wgs_scatter_intervals = Channel.from(file("${params.wgs_scatter_interval_dir}/*/scattered.interval_list")).map{ it -> [ file("${it}").getParent().getBaseName(), it ]  }
        
        ch_genome_fasta = Channel.value(file(params.genome_fasta))
        ch_genome_dict = Channel.value(file(params.genome_dict))
        ch_genome_index = Channel.value(file(params.genome_index))
        
        ch_tiny_vcf = Channel.value(file(params.tiny_vcf))
        ch_tiny_vcf_index = Channel.value(file(params.tiny_vcf_index))
        
        ch_gnomad_vcf = Channel.value(file(params.gnomad_vcf))
        ch_gnomad_vcf_index = Channel.value(file(params.gnomad_vcf_index))

        ch_hotspot_given_alleles_vcf = Channel.value(file(params.hotspot_given_alleles_vcf))
        ch_hotspot_given_alleles_vcf_index = Channel.value(file(params.hotspot_given_alleles_vcf_index))
        
        ch_vep_data_dir = Channel.value(file(params.vep_data_dir))
        
        ch_default_WXS_pon_vcf = Channel.value(file(params.default_WXS_pon_vcf))
        ch_default_WXS_pon_vcf_index = Channel.value(file(params.default_WXS_pon_vcf_index))
        
        ch_default_WGS_pon_vcf = Channel.value(file(params.default_WGS_pon_vcf))
        ch_default_WGS_pon_vcf_index = Channel.value(file(params.default_WGS_pon_vcf_index))


        sample_bams = make_samples ( aligned_bams, params.dna_legacy_file_sheet  )

        callPON( sample_bams.combine(ch_wgs_scatter_intervals)
                    .filter{        aliquot, bam, bai, patient, source, seqtype, gender, tmnm, interval_name, interval ->
                                    tmnm == 'normal' },
                ch_genome_fasta, ch_genome_dict, ch_genome_index )
        callPON.out
            .map{           aliquot, patient, source, seqtype, gender, tmnm, interval_name, vcf, vcfidx, cmds ->
                    return[ aliquot, patient, source, seqtype, gender, tmnm,                vcf, vcfidx ] }.groupTuple( by: [0, 1, 2, 3, 4, 5]) | \
        mergePON

        mergePON.out
            .map{           aliquot, patient, source, seqtype, gender, tmnm, vcf, vcfidx, cmds ->
                return[ source, seqtype, vcf,  vcfidx ] }.groupTuple( by: [0, 1]) | \
        createPON

        collectArtifacts(sample_bams, ch_genome_fasta, ch_genome_dict, ch_genome_index)

        sample_bams.map{ aliquot_barcode, bam_path, bai_path, patient_barcode, source_barcode, sequence_type, gender, tumor_or_normal -> return [aliquot_barcode, bam_path, bai_path ]}
            .combine(collectArtifacts.out
                    .map{       aliquot, patient, source, seqtype, gender, tmnm, alt_table, ref_metrics, alt_metrics, cmds ->
                        return[ aliquot, patient, source, seqtype, gender, tmnm, alt_table, ref_metrics, alt_metrics    ]}
                    , by:0) | orientationPriors
 
        pileupSummaryKnownSites(sample_bams, ch_tiny_vcf, ch_tiny_vcf_index)


        make_estimateContamination_input( make_pairs ( params.dna_legacy_file_sheet ), sample_bams, pileupSummaryKnownSites.out ) | \
        estimateContamination


        // -----------------------------------------
        //
        // Run callSingleSampleSNV with defaultPon
        //
        // -----------------------------------------
        callSingleSampleSNV_defaultPon( make_callSingleSampleSNV_input( sample_bams, params.dna_legacy_file_sheet, createPON.out, orientationPriors.out, ch_wgs_scatter_intervals ),
                                        ch_genome_fasta, ch_genome_dict, ch_genome_index,
                                        ch_gnomad_vcf, ch_gnomad_vcf_index,
                                        ch_hotspot_given_alleles_vcf, ch_hotspot_given_alleles_vcf_index,
                                        ch_default_WXS_pon_vcf, ch_default_WXS_pon_vcf_index,
                                        ch_default_WGS_pon_vcf, ch_default_WGS_pon_vcf_index )

        mergeSingleSampleSNV_defaultPon(make_mergeSingleSampleSNV_input(callSingleSampleSNV_defaultPon.out))

        filterSingleSampleSNV_defaultPon(make_filterSingleSampleSNV_input(estimateContamination.out, mergeSingleSampleSNV_defaultPon.out))

        vcfToMafssSNV_defaultPon(make_vcfToMafssSNV_input(filterSingleSampleSNV_defaultPon.out, make_pairs ( params.dna_legacy_file_sheet )),
                                ch_genome_fasta, ch_genome_dict, ch_genome_index,
                                ch_gnomad_vcf, ch_gnomad_vcf_index,
                                ch_vep_data_dir)



     
    emit:
        vcfToMafssSNV_defaultPon_out = vcfToMafssSNV_defaultPon.out
        filtered_vcf = filterSingleSampleSNV_defaultPon.out
}



process callPON {

    tag "${aliquot_barcode}-${interval_name}"

    publishDir "${params.scratch_dir}/results/${params.workflow_name}/callPON/${aliquot_barcode}/${interval_name}", pattern: "*{command,exitcode}*", mode: 'copy'

    label "callPON"
    
    input:
        tuple   val(aliquot_barcode), path(bam_path), path(bai_path), val(patient_barcode), val(source_barcode), val(sequence_type), val(gender), val(tumor_or_normal), val(interval_name), path(interval_list)
        file(genome_fasta)
        file(genome_dict)
        file(genome_index)
    output:
        tuple   val(aliquot_barcode),                                 val(patient_barcode), val(source_barcode), val(sequence_type), val(gender), val(tumor_or_normal), val(interval_name), path("${aliquot_barcode}.${interval_name}.pon.vcf"), path("${aliquot_barcode}.${interval_name}.pon.vcf.idx"), path("*{command,exitcode}*",hidden: true)
    
    script:
        """
        #!/bin/bash
        gatk --java-options ${params.callPON_java_arg} Mutect2 \
        -R ${genome_fasta} \
        -I ${bam_path} \
        -L ${interval_list} \
        --seconds-between-progress-updates ${params.seconds_between_progress_updates} \
        -O ${aliquot_barcode}.${interval_name}.pon.vcf
        """
        //        source activate mutect2

}

process mergePON {
        
    tag "${aliquot_barcode}"

    publishDir "${params.scratch_dir}/results/${params.workflow_name}/mergePON/${aliquot_barcode}", mode: 'copy'

    label "mergePON"
    
    input:
        tuple   val(aliquot_barcode), val(patient_barcode), val(source_barcode), val(sequence_type), val(gender), val(tumor_or_normal), path(vcf_files), path(vcf_idx_files)
    
    output:
        tuple   val(aliquot_barcode), val(patient_barcode), val(source_barcode), val(sequence_type), val(gender), val(tumor_or_normal), path("${aliquot_barcode}.pon.vcf"), path("${aliquot_barcode}.pon.vcf.idx"), path("*{command,exitcode}*",hidden: true)
    
    script:
        def input_files = vcf_files.collect{"--INPUT ${it}"}.join(' ')
        """
        #!/bin/bash
        source activate mutect2
        gatk --java-options ${params.mergePON_java_arg} MergeVcfs \
        ${input_files} \
        -O ${aliquot_barcode}.pon.vcf
        """

}

process createPON {
    
    tag "${source_barcode}-${sequence_type}"

    publishDir "${params.scratch_dir}/results/${params.workflow_name}/createPON/${source_barcode}/${sequence_type}", mode: 'copy'

    label "createPON"
    
    input:
        tuple   val(source_barcode), val(sequence_type), path(vcf_files), path(vcf_idx_files)

    output:
        tuple   val(source_barcode), val(sequence_type), path("${source_barcode}-${sequence_type}.pon.vcf"), path("${source_barcode}-${sequence_type}.pon.vcf.idx"), path("*{command,exitcode}*", hidden: true)
    
    when:
        source_barcode

    script:
        def input_vcfs = vcf_files.collect{"-vcfs ${it}"}.join(' ')
        """
        #!/bin/bash
        source activate mutect2
        ${params.nf_home}/opt/gatk/gatk-4.1.0.0/gatk --java-options ${params.createPON_java_arg} CreateSomaticPanelOfNormals \
        ${input_vcfs} \
        --duplicate-sample-strategy CHOOSE_FIRST \
        -O ${source_barcode}-${sequence_type}.pon.vcf \
        --QUIET true
        """
}

process collectArtifacts {
    
    tag "${aliquot_barcode}"

    publishDir "${params.scratch_dir}/results/${params.workflow_name}/collectArtifacts/${aliquot_barcode}", mode: 'copy'

    label "collectArtifacts"
    
    input:
        tuple   val(aliquot_barcode), path(bam_path), path(bai_path), val(patient_barcode), val(source_barcode), val(sequence_type), val(gender), val(tumor_or_normal)
        file(genome_fasta)
        file(genome_dict)
        file(genome_index)

    output:
        tuple   val(aliquot_barcode),                                 val(patient_barcode), val(source_barcode), val(sequence_type), val(gender), val(tumor_or_normal), path("${aliquot_barcode}.alt.tsv"), path("${aliquot_barcode}.ref.metrics"), path("${aliquot_barcode}.alt.metrics"), path("*{command,exitcode}*", hidden: true)
    
    script:

        """
        #!/bin/bash
        source activate mutect2
        ${params.nf_home}/opt/gatk/gatk-4.1.0.0/gatk --java-options ${params.collectArtifacts_java_arg} CollectF1R2Counts \
            -R ${genome_fasta} \
            -I ${bam_path} \
            -alt-table ${aliquot_barcode}.alt.tsv \
            -ref-hist ${aliquot_barcode}.ref.metrics \
            -alt-hist ${aliquot_barcode}.alt.metrics
        """


}


process orientationPriors {
    
    tag "${aliquot_barcode}"

    publishDir "${params.scratch_dir}/results/${params.workflow_name}/orientationPriors/${aliquot_barcode}", mode: 'copy'
    
    label "orientationPriors"
    
    input:
        tuple   val(aliquot_barcode), path(bam_path), path(bai_path), val(patient_barcode), val(source_barcode), val(sequence_type), val(gender), val(tumor_or_normal), path(alt_table), path(ref_metrics), path(alt_metrics)

    output:
        tuple   val(aliquot_barcode),                                 val(patient_barcode), val(source_barcode), val(sequence_type), val(gender), val(tumor_or_normal), path("${aliquot_barcode}.priors.tsv"), path("*{command,exitcode}*", hidden: true)
    
    script:
        //def gatk = "/projects/kimh/opt/gatk/gatk-4.1.0.0/gatk"
        """
        #!/bin/bash
        source activate mutect2
        ${params.nf_home}/opt/gatk/gatk-4.1.0.0/gatk --java-options ${params.orientationPriors_java_arg} LearnReadOrientationModel \
            -alt-table ${alt_table} \
            -ref-hist ${ref_metrics} \
            -alt-hist ${alt_metrics} \
            -O ${aliquot_barcode}.priors.tsv
        """


}

process pileupSummaryKnownSites {
    
    tag "${aliquot_barcode}"

    publishDir "${params.scratch_dir}/results/${params.workflow_name}/pileupSummaryKnownSites/${aliquot_barcode}", mode: 'copy'
    
    label "pileupSummaryKnownSites"
    
    input:
        tuple   val(aliquot_barcode), path(bam_path), path(bai_path), val(patient_barcode), val(source_barcode), val(sequence_type), val(gender), val(tumor_or_normal)
        file(tiny_vcf)
        file(tiny_vcf_index)
  
    output:
        tuple   val(aliquot_barcode),                                 val(patient_barcode), val(source_barcode), val(sequence_type), val(gender), val(tumor_or_normal), path("${aliquot_barcode}.pileupsummaries.txt"), path("*{command,exitcode}*", hidden: true)
    
    script:
        """
        #!/bin/bash
        source activate mutect2
        ${params.nf_home}/opt/gatk/gatk-4.1.0.0/gatk --java-options ${params.pileupSummaryKnownSites_java_arg} GetPileupSummaries \
        -I ${bam_path} \
        -V ${tiny_vcf} \
        -L ${tiny_vcf} \
        -O ${aliquot_barcode}.pileupsummaries.txt \
        --seconds-between-progress-updates ${params.seconds_between_progress_updates}
        """
}

process estimateContamination {
    
    tag "${pair_barcode}"

    publishDir "${params.scratch_dir}/results/${params.workflow_name}/estimateContamination/${pair_barcode}", mode: 'copy'

    label "estimateContamination"
    
    input:
        tuple   val(pair_barcode), val(patient_barcode), val(tumor_barcode), val(normal_barcode), path(tumor_pileupsum_file), path(normal_pileupsum_file)

    output:
        tuple   val(pair_barcode), val(patient_barcode), val(tumor_barcode), val(normal_barcode), path("${pair_barcode}.contamination.txt"), path("${pair_barcode}.segmentation.txt"), path("*{command,exitcode}*", hidden: true)
    
    when:
        pair_barcode

    script:
        """
        #!/bin/bash
        source activate mutect2
        ${params.nf_home}/opt/gatk/gatk-4.1.0.0/gatk --java-options ${params.estimateContamination_java_arg} CalculateContamination \
        -I ${tumor_pileupsum_file} \
        --matched-normal ${normal_pileupsum_file} \
        --output ${pair_barcode}.contamination.txt \
        --tumor-segmentation ${pair_barcode}.segmentation.txt
        """
}

process callSingleSampleSNV_defaultPon {
    
    tag "${pair_barcode}-${interval_name}"

    publishDir "${params.scratch_dir}/results/${params.workflow_name}/callSingleSampleSNV_defaultPon/${pair_barcode}/${interval_name}"

    label "callSingleSampleSNV_defaultPon"
    
    input:
        tuple   val(pair_barcode), val(tumor_barcode), val(tumor_rgsm), path(tumor_bam), path(tumor_bai), val(normal_barcode), val(normal_rgsm), path(normal_bam), path(normal_bai), val(patient_barcode), path(tumor_orientation_prior), val(pon_tmp), val(pon_tmpidx), val(source_barcode), val(sequence_type), val(interval_name), path(interval_list)
        file(genome_fasta)
        file(genome_dict)
        file(genome_index)
        file(gnomad_vcf)
        file(gnomad_vcf_index)
        file(given_alleles_vcf)
        file(given_alleles_vcf_index)
        
        file(default_WXS_pon_vcf)
        file(default_WXS_pon_vcf_index)

        file(default_WGS_pon_vcf)
        file(default_WGS_pon_vcf_index)

    output:
        tuple   val(pair_barcode), val(source_barcode), val(sequence_type), val(interval_name), path(interval_list), path("${pair_barcode}.${interval_name}.vcf"), path("${pair_barcode}.${interval_name}.vcf.idx"), path("*{command,exitcode}*", hidden: true)

    when:
        normal_bam

    script:
        def pon = selectMutect2DefaultPon ( sequence_type, default_WGS_pon_vcf, default_WXS_pon_vcf )
        """
        #!/bin/bash
        source activate mutect2
        ${params.nf_home}/opt/gatk/gatk-4.1.0.0/gatk --java-options ${params.callSingleSampleSNV_java_arg} Mutect2 \
            -R ${genome_fasta} \
            -I ${tumor_bam} \
            -I ${normal_bam} \
            --orientation-bias-artifact-priors ${tumor_orientation_prior} \
            -L ${interval_list} \
            --normal-sample ${normal_rgsm} \
            --panel-of-normals ${pon} \
            --germline-resource ${gnomad_vcf} \
            --genotyping-mode GENOTYPE_GIVEN_ALLELES \
            --genotype-filtered-alleles true \
            --alleles ${given_alleles_vcf} \
            --dont-use-soft-clipped-bases true \
            --standard-min-confidence-threshold-for-calling 20 \
            -O ${pair_barcode}.${interval_name}.vcf \
            --seconds-between-progress-updates ${params.seconds_between_progress_updates}
 
        """

}

process mergeSingleSampleSNV_defaultPon {
    
    tag "${pair_barcode}"

    publishDir "${params.scratch_dir}/results/${params.workflow_name}/mergeSingleSampleSNV_defaultPon/${pair_barcode}", mode: 'copy'

    label "mergeSingleSampleSNV"
    
    input:
        tuple   val(pair_barcode), val(source_barcode), val(sequence_type), path(interval_vcfs), path(interval_vcfidxs)

    output:
        tuple   val(pair_barcode), val(source_barcode), val(sequence_type), path("${pair_barcode}.vcf"), path("${pair_barcode}.vcf.idx"), path("*{command,exitcode}*",hidden: true)

    when:
        interval_vcfs

    script:
        def input_vcfs = interval_vcfs.collect{"-I ${it}"}.join(' ')
        """
        #!/bin/bash
        source activate mutect2
        ${params.nf_home}/opt/gatk/gatk-4.1.0.0/gatk --java-options ${params.mergeSingleSampleSNV_java_arg} MergeVcfs \
        ${input_vcfs} \
        -O ${pair_barcode}.vcf
        """
}


process filterSingleSampleSNV_defaultPon {
    
    tag "${pair_barcode}"

    publishDir "${params.scratch_dir}/results/${params.workflow_name}/filterSingleSampleSNV_defaultPon/${pair_barcode}", mode: 'copy'
    
    label "filterSingleSampleSNV"
    
    input:
        tuple   val(pair_barcode), val(patient), val(source_barcode), val(sequence_type), path(unfiltered_vcf), path(unfiltered_vcfidx), path(pair_contamination_file), path(pair_segmentation_file)

    output:
        tuple   val(pair_barcode), val(patient), val(source_barcode), val(sequence_type), path("${pair_barcode}.filterstats.tsv"), path("${pair_barcode}.filtered.vcf.gz"), path("${pair_barcode}.filtered.vcf.gz.tbi"), path("*{command,exitcode}*", hidden: true)
    
    when:
        pair_barcode

    script:
        """
        #!/bin/bash
        source activate mutect2
        ${params.nf_home}/opt/gatk/gatk-4.1.0.0/gatk --java-options ${params.filterSingleSampleSNV_java_arg} FilterMutectCalls \
            -V ${unfiltered_vcf} \
            --tumor-segmentation ${pair_segmentation_file} \
            --contamination-table ${pair_contamination_file} \
            --stats ${pair_barcode}.filterstats.tsv \
            -O ${pair_barcode}.filtered.vcf.gz \
            --seconds-between-progress-updates ${params.seconds_between_progress_updates}
        """
}


process vcfToMafssSNV_defaultPon {
    
    tag "${pair_barcode}"

    publishDir "${params.work_dir}/results/${params.workflow_name}/vcfToMafssSNV_defaultPon/${pair_barcode}", mode: 'copy'

    publishDir "${params.scratch_dir}/results/${params.workflow_name}/vcfToMafssSNV_defaultPon/${pair_barcode}", mode: 'copy'
    
    label "vcfToMafssSNV"
    
    input:
        tuple   val(pair_barcode), val(patient), val(source_barcode), val(sequence_type), path(filtered_vcf), path(filtered_vcfidx), val(tumor_sm), val(normal_sm)
        file(genome_fasta)
        file(genome_dict)
        file(genome_index)
        file(gnomad_vcf)
        file(gnomad_vcf_index)
        file(vep_data_dir)


    output:
        tuple   val(pair_barcode), val(patient), val(source_barcode), val(sequence_type), path("${pair_barcode}.filtered.maf"), path("*{command,exitcode}*", hidden: true)
    
    when:
        pair_barcode

    script:
        def FILTER_VCF = gnomad_vcf
        def SPECIES = "homo_sapiens"
        def NCBI_BUILD = params.genome
        def VEP_FORKS = 4
        def VEP_DATA = vep_data_dir
        def VEP_CACHE_VER = params.vep_cache_ver
        def TM_ID = tumor_sm
        def NM_ID = normal_sm
        def MAF_CENTER = "JAX"
        def REF_FASTA = genome_fasta
        """
        #!/bin/bash
        source activate ensembl_vep
        export VEP_PATH=\$(dirname `which vcf2maf.pl`)        
        gunzip -f ${filtered_vcf} && \
        vcf2maf.pl \
        --input-vcf ${pair_barcode}.filtered.vcf \
        --output-maf ${pair_barcode}.filtered.maf \
        --vep-path \${VEP_PATH} \
        --vep-data ${VEP_DATA} \
        --vep-forks ${VEP_FORKS} \
        --filter-vcf ${FILTER_VCF} \
        --tumor-id ${TM_ID} \
        --normal-id ${NM_ID} \
        --species ${SPECIES} \
        --ncbi-build ${NCBI_BUILD} \
        --maf-center ${MAF_CENTER} \
        --ref-fasta ${REF_FASTA} \
        --cache-version ${VEP_CACHE_VER}
        """
}


              // functions



def make_pairs ( prealign_sample_sheet  ) {
    
    tumors_ch = Channel.fromPath(prealign_sample_sheet)
        .splitCsv(header:true)
        .map{ row-> tuple(  row.aliquot_barcode,
                            row.sample_barcode,
                            row.patient_barcode,
                            row.tumor_or_normal,
                            row.source_barcode,
                            row.sequence_type,
                            row.gender ) }
        .filter {   output_barcode, sample, patient, tmnm, source, seqtype, gender -> 
                    tmnm == "tumor" }
        .map {      output_barcode, sample, patient, tmnm, source, seqtype, gender -> 
        return([    patient, output_barcode, sample, source, seqtype, gender ] ) }
        .unique()

    normals_ch = Channel.fromPath(prealign_sample_sheet)
        .splitCsv(header:true)
        .map{ row-> tuple(  row.aliquot_barcode,
                            row.sample_barcode,
                            row.patient_barcode,
                            row.tumor_or_normal ) }
        .filter {   output_barcode, sample, patient, tmnm -> 
                    tmnm == "normal" }
        .map {      output_barcode, sample, patient, tmnm -> 
        return([    patient, output_barcode, sample ] ) }
        .unique()

    pairs_ch = tumors_ch.combine( normals_ch, by: [0] )
        .map {      patient, tumor_aliquot, tumor_sample, source, seqtype, gender, normal_aliquot, normal_sample ->
                    pair = tumor_aliquot + "-" + normal_aliquot
        return [    pair, tumor_aliquot, tumor_sample, normal_aliquot, normal_sample, patient, source, seqtype ] }
        .unique()

    return pairs_ch

}


def make_samples ( aligned_bam_ch, prealign_sample_sheet  ) {
    
    samples_ch = Channel
        .fromPath(prealign_sample_sheet)
        .splitCsv(header:true)
        .map{ row-> tuple(  row.aliquot_barcode,
                            row.patient_barcode,
                            row.source_barcode,
                            row.sequence_type,
                            row.gender,
                            row.tumor_or_normal ) }
        .map {      aliquot, patient, source, seqtype, gender, tmnm -> 
        return([    aliquot, patient, source, seqtype, gender, tmnm ] ) }
        .distinct()

    aligned_bam_ch
        .map{           aliquot, bam, bai, md5 ->
                return [aliquot, bam, bai ]}
        .combine( samples_ch, by: 0 )
        .map{           aliquot, bam, bai, patient, source, seqtype, gender, tmnm ->
                return [aliquot, bam, bai, patient, source, seqtype, gender, tmnm ]}
        
}


def make_callSingleSampleSNV_input (sample_bams, samplesheet_path, createPON_out, orientationPriors_out, ch_wgs_scatter_intervals) {

    // ============================================================================
    //
    // Create input to callSingleSampleSNV
    // Run callSingleSampleSNV.
    //        
    // ============================================================================
    
    // Prepare input to callSingleSampleSNV =====================
    pairs = make_pairs ( samplesheet_path )
    // OUTPUT: pair, tm_aliquot, tm_sample, nm_aliquot, nm_sample, patient, source, seqtype

    // Add tumor_bam and tumor_bai
    pairs_wbam = pairs
        .map {      pair, tm_aliquot,       tm_sample, nm_aliquot, nm_sample, patient, source, seqtype ->
            return[       tm_aliquot, pair, tm_sample, nm_aliquot, nm_sample, patient, source, seqtype ] }
        .combine ( sample_bams
                    .map {      aliquot, bam, bai, patient, source, seqtype, gender, tmnm ->
                        return[ aliquot, bam, bai ] }, by: 0 )
        .map {                        tm_aliquot, pair, tm_sample, nm_aliquot, nm_sample, patient, source, seqtype, tm_bam, tm_bai ->
            return[ nm_aliquot, pair, tm_aliquot,       tm_sample,             nm_sample, patient, source, seqtype, tm_bam, tm_bai ] }
        .combine ( sample_bams
                    .map {      aliquot, bam, bai, patient, source, seqtype, gender, tmnm ->
                        return[ aliquot, bam, bai ] }, by: 0 )
        .map {      nm_aliquot, pair, tm_aliquot, tm_sample,                             nm_sample,                 patient, source, seqtype, tm_bam, tm_bai, nm_bam, nm_bai ->
            return[             pair, tm_aliquot, tm_sample, tm_bam, tm_bai, nm_aliquot, nm_sample, nm_bam, nm_bai, patient, source, seqtype ] }
    
    // Add tm_prior and pon to pairs_wbam
    pairs_wbam_wprior_wpon_tmp = pairs_wbam
        .map {      pair, tm_aliquot, tm_sample, tm_bam, tm_bai, nm_aliquot, nm_sample, nm_bam, nm_bai, patient, source, seqtype ->
            return[ tm_aliquot, pair, tm_sample, tm_bam, tm_bai, nm_aliquot, nm_sample, nm_bam, nm_bai, patient, source, seqtype ] }
        .combine ( orientationPriors_out
                    .map{ aliquot, patient, source, seqtype, gender, tmnm, prior, cmds ->
                        return [ aliquot, prior ] }, by: 0 )
        .map {      tm_aliquot, pair, tm_sample, tm_bam, tm_bai, nm_aliquot, nm_sample, nm_bam, nm_bai, patient, source, seqtype, prior ->
                return[ source, seqtype, pair, tm_aliquot, tm_sample, tm_bam, tm_bai, nm_aliquot, nm_sample, nm_bam, nm_bai, patient, prior ] }
        .combine ( createPON_out
                .map{           source, seqtype, pon, ponidx, cmds ->
                    return [    source, seqtype, pon, ponidx    ]}
                , by: [0,1] 
                )

    pairs_wbam_wprior_wpon_tmp
        .map { source, seqtype, pair, tm_aliquot, tm_sample, tm_bam, tm_bai, nm_aliquot, nm_sample, nm_bam, nm_bai, patient, prior, pon, ponidx ->
            return[pair, tm_aliquot, tm_sample, tm_bam, tm_bai, nm_aliquot, nm_sample, nm_bam, nm_bai, patient, prior, pon, ponidx, source, seqtype]}
        .combine(ch_wgs_scatter_intervals)

}


def make_mergeSingleSampleSNV_input (callSingleSampleSNV_out) {
    callSingleSampleSNV_out
        .map{       output_barcode, source, seqtype, interval_name, interval, interval_vcf, interval_vcfidx, cmds ->
            return[ output_barcode, source, seqtype, interval_vcf, interval_vcfidx]}
        .groupTuple(by:[0,1,2])     

}


def make_estimateContamination_input( pairs, sample_bams, pileupSummaryKnownSites_out ) {
    
    ch_bams = sample_bams.map{ aliquot_barcode, bam_path, bai_path, patient_barcode, source_barcode, sequence_type, gender, tumor_or_normal ->
                       return[ aliquot_barcode, bam_path, bai_path ]}

    pileupSummaryKnownSites_out_wbams = ch_bams
        .combine(pileupSummaryKnownSites_out
                    .map{   aliquot, patient, source, seqtype, gender, tmnm, pileup, cmds -> 
                    return[ aliquot, patient, source, seqtype, gender, tmnm, pileup]
                    },by:0)
    // pileupSummaryKnownSites_out:         aliquot,           patient, source, seqtype, gender, tmnm, pileup
    // pileupSummaryKnownSites_out_wbams:   aliquot, bam, bai, patient, source, seqtype, gender, tmnm, pileup

    pairs.map{ pair, tm_aliquot, tm_sample, nm_aliquot, nm_sample, patient, source, seqtype -> return[tm_aliquot, pair, patient, nm_aliquot] }
        .combine( pileupSummaryKnownSites_out_wbams.map{ aliquot, bam, bai, patient, source, seqtype, gender, tmnm, pileup -> return[aliquot, pileup]}, 
                    by:0)
        .map{ tm_aliquot, pair, patient, nm_aliquot, tm_pileup -> return[nm_aliquot, pair, patient, tm_aliquot, tm_pileup] }
        .combine( pileupSummaryKnownSites_out_wbams.map{ aliquot, bam, bai, patient, source, seqtype, gender, tmnm, pileup -> return[aliquot, pileup]}, 
                    by:0)
        .map{       nm_aliquot, pair, patient, tm_aliquot,              tm_pileup, nm_pileup -> 
            return[             pair, patient, tm_aliquot, nm_aliquot,  tm_pileup, nm_pileup]}

}


def make_filterSingleSampleSNV_input(estimateContamination_out, mergeSingleSampleSNV_out) {
        
    m2_ms_out2=mergeSingleSampleSNV_out.map{pair, source, seqtype, m2_vcf, m2_vcfidx, cmds -> return[pair, source, seqtype, m2_vcf, m2_vcfidx]}
    est_out2=estimateContamination_out.map{pair, patient, tm_aliquot, nm_aliquot, pair_contam, pair_seg, cmds -> return[pair, patient, pair_contam, pair_seg]}
    
    m2_ms_out2.combine(est_out2,by: 0)
        .map{   pair, source, seqtype, m2_vcf, m2_vcfidx, patient, pair_contam, pair_seg -> 
        return[ pair, patient, source, seqtype, m2_vcf, m2_vcfidx, pair_contam, pair_seg ]
    }
    //filterSingleSampleSNV_input.view()
    //pair, patient, source, seqtype, m2_vcf, m2_vcfidx, pair_contam, pair_seg

}


def make_vcfToMafssSNV_input(filterSingleSampleSNV_out, pairs) {

    filterSingleSampleSNV_out
            .map{       pair, patient, source, seqtype, pair_filterstats, pair_filtered_vcf_gz, pair_filtered_vcf_gz_tbi, cmds -> 
                return[ pair, patient, source, seqtype, pair_filtered_vcf_gz, pair_filtered_vcf_gz_tbi]}
            .combine(pairs.map{pair, tm_aliquot, tm_sm, nm_aliquot, nm_sm, patient, source, seqtype -> return[pair, tm_sm, nm_sm]}, by:0)

}

def hasExtension(it, extension) {

    it.toString().toLowerCase().endsWith(extension.toLowerCase())

}


// Return file if it exists
def returnFile(it) {
    if (!file(it).exists()) exit 1, "Missing file in TSV file: ${it}, see --help for more information"
    return file(it)
}


def selectMutect2DefaultPon ( seq_type, default_WGS_pon_vcf, default_WXS_pon_vcf ) {

    if ( seq_type == "WGS" ) {

        return default_WGS_pon_vcf

    } else if ( seq_type == "WXS" ) {

        return default_WXS_pon_vcf

    } else {
        throw new IllegalArgumentException("Unknown seq_type ${seq_type}")
    }
}

