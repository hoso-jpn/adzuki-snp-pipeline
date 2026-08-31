/**
 * Issue #52: reads the GS lineage's default container references straight
 * out of `conf/containers.config` -- the single source of truth those
 * defaults are declared in -- instead of copying the pinned digests into
 * test code.
 *
 * A test that hardcoded its own copy of a production digest could not
 * prove what Issue #52's acceptance criterion actually asks for ("change
 * only conf/containers.config and both the running container and the
 * manifest follow"): it would pass just as happily against a pipeline
 * that had gone back to hardcoding the same digest in every module, and
 * it would have to be edited by hand every time a pinned image was
 * legitimately bumped. Parsing the source of truth makes the assertion
 * "the manifest equals whatever conf/containers.config currently says",
 * which is the property under test.
 */
class GsContainerDefaults {

    /** Tool keys `conf/containers.config` is required to declare. */
    static final Set<String> TOOLS = ['bcftools', 'gatk', 'python'] as Set

    private static final java.util.regex.Pattern ENTRY =
        ~/(?m)^\s*(bcftools|gatk|python)\s*:\s*'([^']+)'\s*,?\s*$/

    /**
     * The `params.containers` map as currently declared in
     * `conf/containers.config`, keyed by tool name.
     *
     * Deliberately a text parse rather than a Nextflow config evaluation:
     * the point is to read the file a human edits when bumping an image,
     * with no chance of an override from a profile or another config file
     * silently substituting a different value into the "default" the test
     * is asserting against.
     */
    static Map<String, String> read(Object projectDir) {
        def text = new File("${projectDir}/conf/containers.config").text
        def defaults = [:]
        ENTRY.matcher(text).each { _match, tool, reference -> defaults[tool] = reference }

        assert defaults.keySet() == TOOLS:
            "conf/containers.config did not declare exactly ${TOOLS}, got ${defaults.keySet()}"
        defaults.each { tool, reference ->
            assert reference.contains('@sha256:'):
                "conf/containers.config's '${tool}' default is not digest-pinned: ${reference}"
        }
        return defaults
    }

    /**
     * The same image, referenced by digest alone (tag dropped).
     *
     * This is the transformation the Issue #52 test fixtures apply to
     * every default when they exercise "change the source of truth" --
     * `tests/pipeline/fixtures/gs_container_default_change.config` and the
     * per-module override fixtures each apply the identical
     * `replaceFirst(/:[^@]+@/, '@')` in Nextflow config, and the tests use
     * this method to say what they therefore expect to see recorded.
     *
     * The resulting reference is a distinct string from the default while
     * still naming the exact same image content by digest: the run needs
     * no additional image pull, and the test never depends on a mutable
     * tag resolving to the same bytes it did when the test was written.
     */
    static String digestOnly(String reference) {
        return reference.replaceFirst(/:[^@]+@/, '@')
    }
}
