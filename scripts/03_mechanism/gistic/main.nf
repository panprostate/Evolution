params.trajectories = ["ordering_1", "ordering_2", "ordering_3", "all"]
// params.trajectories = ["ordering_3"]
params.trajectories_dir = "${projectDir}/../../../outputs/02_trajectories/"
params.bb_copynumber_dir = "${projectDir}/../../../data/raw/Somatic_variants/SCNA/SCNA_with_Brass_Delly_SVs_20260212/Subclonal_SCNA"
params.ta_threshold = 0.2
params.td_threshold = 0.2
params.refgene = "${projectDir}/refgenefiles/hg19.mat"

process BB_TO_SEGFILE {
    publishDir "${projectDir}/gistic_input"
    conda "${projectDir}/envs/ppcg.yml"
    memory 32.Gb
    time 30.min
    cpus 10

    input: 
    file bb_copynumber_dir
    file trajectories_dir
    file src_dir
    val trajectory


    output: 
    file "*gistic_input.txt"

    script: 
    """
    Rscript ${projectDir}/scripts/bb_to_segfile.R \
        --ncores ${task.cpus} \
        --bb_copynumber_dir ${bb_copynumber_dir} \
        --trajectories_dir ${trajectories_dir} \
        --trajectory ${trajectory}
    """
}

process RUN_GISTIC {
    publishDir "${projectDir}/gistic_output"
    conda "${projectDir}/envs/gistic.yml"
    memory 60.Gb
    time 8.hours
    cpus 8

    input:
    tuple file(segfile), val(trajectory)
    file refgene_file
    val ta_threshold 
    val td_threshold

    output:
    path "${trajectory}/*"

    script:
    """
    unset LC_MESSAGES
    unset LC_PAPER
    unset LC_NAME
    unset LC_ADDRESS
    unset LC_TELEPHONE
    unset LC_MEASUREMENT
    unset LC_IDENTIFICATION

    export LC_ALL=C
    export LANG=C
    export LANGUAGE=C

    echo "After cleanup:"
    locale || true

    mkdir -p ${trajectory}
    gistic2 -b ${trajectory} \
            -seg ${segfile} \
            -refgene ${refgene_file} \
            -ta ${ta_threshold} -td ${td_threshold}
    """
}

workflow {
    segment_ch = BB_TO_SEGFILE(
        file(params.bb_copynumber_dir),
        file(params.trajectories_dir),
        file("${projectDir}/src/"),
        Channel.fromList(params.trajectories) 
    )

    segment_ch.view()
    segment_ch = segment_ch.map { file -> 
        def prefix = file.baseName.replaceAll('_gistic_input', '')
        [file, prefix]
    }

    RUN_GISTIC(
        segment_ch, 
        file(params.refgene),
        Channel.value(params.ta_threshold), Channel.value(params.td_threshold)
    )
}