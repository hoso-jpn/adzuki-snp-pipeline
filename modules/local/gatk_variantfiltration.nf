// Issue #17: GATK_VARIANTFILTRATION *is* the hard-filtering step of this
// pipeline (and, aliased as GATK_VARIANTFILTRATION_GS, of the GS panel
// lineage), so `filters` is a required contract rather than an optional
// hint: an empty list is a configuration error that would silently
// publish a cohort.<type>.filtered.vcf.gz in which every record is
// unfiltered PASS, i.e. an artifact whose name claims hard filtering was
// applied when none was. That is strictly worse than failing, so it
// fail-fasts here instead of being treated as a successful no-op.
//
// The checks below run at this process boundary, not only in
// workflows/adzuki_snp_pipeline.nf, because the workflow's literal
// snpHardFilters()/indelHardFilters() lists are only one of the callers
// this module has to be safe for -- the other being direct module use,
// as tests/modules/gatk_variantfiltration.nf.test does.
def validatedHardFilters(meta, filters) {
    def variant_type = meta instanceof Map
        ? (meta['variant_type'] ?: 'unknown')
        : 'unknown'

    // List, not Collection: the contract is an *ordered* list. The
    // position numbers in every diagnostic below ("hard filter 2 of 7"),
    // and the order the --filter-name/--filter-expression pairs reach
    // GATK in, are only meaningful and reproducible if the input has a
    // defined iteration order -- which an unordered Collection such as a
    // Set does not guarantee. snpHardFilters()/indelHardFilters() both
    // return Groovy list literals (ArrayList), so this rejects nothing
    // the pipeline actually passes.
    if (!(filters instanceof List)) {
        error(
            "hard filters for ${variant_type} must be a list, got " +
                "${filters == null ? 'null' : filters.getClass().getName()}"
        )
    }

    if (filters.isEmpty()) {
        error(
            "no hard filters were configured for ${variant_type}; " +
                "GATK_VARIANTFILTRATION is the hard-filtering step, so an " +
                "empty filter list is a configuration error, not a no-op"
        )
    }

    def seen_names = [:]

    return filters.withIndex().collect { filter, index ->
        def position = index + 1

        if (!(filter instanceof Map)) {
            error(
                "hard filter ${position} of ${filters.size()} for " +
                    "${variant_type} must be a map with 'name' and " +
                    "'expression' keys, got ${filter}"
            )
        }

        def name = requiredFilterField(
            filter,
            'name',
            position,
            filters.size(),
            variant_type,
        )
        def expression = requiredFilterField(
            filter,
            'expression',
            position,
            filters.size(),
            variant_type,
        )

        // A filter name becomes a VCF FILTER-column value and a
        // ##FILTER=<ID=...> header entry, where the format itself
        // forbids whitespace and ';' (the FILTER field's own
        // separator). Unlike an expression -- which is GATK JEXL and
        // legitimately needs '<', '>', '&&', '||', '(', ')', '.', '-',
        // spaces and even single-quoted string literals -- a name has
        // no reason to contain anything outside this set, so it takes
        // the stricter contract. Every filter name this pipeline
        // actually configures (SNP_QD_LOW ... INDEL_READPOSRANKSUM_LOW)
        // already satisfies it.
        if (!(name ==~ /[A-Za-z0-9._-]+/)) {
            error(
                "invalid hard filter name '${name}' for ${variant_type}: " +
                    "expected only letters, digits, '.', '_' or '-' " +
                    "(a VCF FILTER id cannot contain whitespace or ';')"
            )
        }

        // 'PASS' and '.' are the VCF-reserved values for "passed all
        // filters" and "no filters applied". Tagging a record with
        // either would make a *failing* record indistinguishable from a
        // passing one, and GATK_SELECTPASSVARIANTS's --exclude-filtered
        // would then keep records this step meant to drop.
        if (name in ['PASS', '.']) {
            error(
                "hard filter name '${name}' for ${variant_type} is reserved " +
                    "by the VCF FILTER field and cannot be used as a filter id"
            )
        }

        // Two filters sharing a FILTER id make the resulting VCF
        // unauditable: a record tagged 'QD_LOW' no longer identifies
        // which of the two thresholds it actually failed, and the
        // per-tag counts in bin/summarize_filter_qc.py silently merge
        // them.
        if (seen_names.containsKey(name)) {
            error(
                "duplicate hard filter name '${name}' for ${variant_type}: " +
                    "filters ${seen_names[name]} and ${position} share it, so " +
                    "the resulting VCF FILTER tag would not identify which " +
                    "expression a record failed"
            )
        }
        seen_names[name] = position

        return [name: name, expression: expression]
    }
}

def requiredFilterField(filter, field, position, total, variant_type) {
    def context = "hard filter ${position} of ${total} for ${variant_type}"

    if (!filter.containsKey(field)) {
        error("${context} is missing '${field}': ${filter}")
    }

    def value = filter[field]

    if (value == null) {
        error("${context} has a null '${field}': ${filter}")
    }

    // Not a bare Groovy truthiness check (`if (!name || !expression)`),
    // which was the previous contract: truthiness rejects null and ""
    // but *accepts* a whitespace-only "   " (a non-empty String is
    // truthy), and would conversely reject a legitimate value that
    // merely happens to be falsy in Groovy. It also silently accepts
    // non-String values -- a list or a number would be stringified
    // into the GATK argument. So the type and the blank check are both
    // explicit. CharSequence, not String, because the workflow builds
    // expressions by interpolation ("QD < ${params.snp_filter_qd_min}"),
    // which yields a GString.
    if (!(value instanceof CharSequence)) {
        error(
            "${context} has a non-string '${field}' " +
                "(${value.getClass().getName()}): ${filter}"
        )
    }

    def text = value.toString()

    if (text.isBlank()) {
        error("${context} has a blank '${field}': ${filter}")
    }

    return text
}

// Issue #17: POSIX single-quote quoting -- wrap the whole value in
// single quotes and rewrite each embedded single quote as '\'' (close,
// backslash-escaped quote, reopen). Inside single quotes the shell
// treats every other byte literally, so the JEXL operators the hard
// filters legitimately use ('<', '>', '&&', '||', '(', ')', '.', '-',
// spaces) and JEXL's own single-quoted string literals reach GATK
// unchanged -- no scientific expression has to be rewritten to be
// shell-safe, and none can terminate its own quoting and inject a
// command. The previous "--filter-name '${name}'" interpolation could:
// a single quote anywhere in a name or expression (reachable from the
// CLI, since the expressions embed --snp_filter_* params) ended the
// quoted string early and handed the rest to the shell.
def shellQuote(value) {
    return "'" + value.toString().replace("'", "'\\''") + "'"
}

process GATK_VARIANTFILTRATION {
    tag "${meta.id}:${meta.variant_type}"
    label 'process_medium'

    container params.containers.gatk

    input:
    tuple(
        val(meta),
        path(vcf),
        path(vcf_index),
        val(filters)
    )

    output:
    tuple(
        val(meta),
        path("${meta.id}.${meta.variant_type}.filtered.vcf.gz"),
        path("${meta.id}.${meta.variant_type}.filtered.vcf.gz.tbi"),
        emit: vcf
    )
    // Issue #52: this process is aliased as both GATK_VARIANTFILTRATION
    // (primary lineage) and GATK_VARIANTFILTRATION_GS (GS lineage) in
    // workflows/adzuki_snp_pipeline.nf, and each alias can be overridden
    // independently (withName selectors match by the included name, not
    // the module's own process name). Emitting task.container -- read from
    // whichever alias actually invoked this task -- rather than the
    // `container` directive's default lets a consumer distinguish the two
    // aliases' effective containers instead of assuming they always match.
    val(task.container), emit: container_id

    script:
    memory_gb = Math.max(
        1,
        task.memory.toGiga().intValue() - 1
    )

    filter_arguments = validatedHardFilters(meta, filters)
        .collect { filter ->
            "--filter-name ${shellQuote(filter.name)} " +
                "--filter-expression ${shellQuote(filter.expression)}"
        }
        .join(" \\\n        ")

    """
    gatk --java-options "-Xmx${memory_gb}g" VariantFiltration \
        --variant ${vcf} \
        ${filter_arguments} \
        --output ${meta.id}.${meta.variant_type}.filtered.vcf.gz \
        --create-output-variant-index true
    """
}
