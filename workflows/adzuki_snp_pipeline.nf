workflow ADZUKI_SNP_PIPELINE {
    take:
    samples_ch
    reference_ch

    main:
    validated_samples_ch = samples_ch.view {
        meta,
        fastq_1,
        fastq_2 ->
        "Validated sample=${meta.id}, " +
        "read_group=${meta.read_group_id}, " +
        "reads=${fastq_1.name},${fastq_2.name}"
    }

    validated_reference_ch = reference_ch.view {
        meta,
        fasta ->
        "Validated reference=${meta.id}, fasta=${fasta.name}"
    }

    emit:
    samples  = validated_samples_ch
    reference = validated_reference_ch
}
