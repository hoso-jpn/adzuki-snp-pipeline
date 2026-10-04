# Reference bundle validation (Issue #68)

The workflow now gates mapping and variant calling on the sequence itself,
not only agreement between the FAI and dictionary's contig names/lengths.
`VALIDATE_REFERENCE_BUNDLE` streams the uncompressed FASTA, reconstructs the
FAI's name, length, byte offset, bases per line and bytes per line, and checks
every dictionary M5 against the uppercase sequence. Missing M5 is a failure.
Duplicate/empty contigs, unsupported bases and non-indexable wrapping fail.

`BWA_MEM2_INDEX` writes an SHA-256 receipt for the FASTA and all five index
files in the same task that built them, plus the actual BWA version. Validation
checks those bytes and records the effective indexing container. It emits
`reference/<reference_id>.reference_bundle.json`, with a version, passed
result, sequence fingerprints, file hashes and a deterministic bundle
fingerprint. The final run manifest includes this file's checksum. Failure
emits no verified FAI/index channel, so mapping and HaplotypeCaller cannot run.

## Reuse

First run with no `bwa_index_prefix`: the controlled index build emits the
manifest automatically. Preserve the FASTA, FAI, dictionary, five BWA files
and this manifest as one immutable bundle. Reuse it by passing
`--bwa_index_prefix`, `--reference_fai`, `--reference_dict` and
`--reference_bundle_manifest`. The manifest must match *all* supplied bytes;
do not regenerate just the dictionary (its header bytes may change) while
retaining the previous manifest. A bundle with missing/unknown provenance
must be regenerated from FASTA; renaming its files does not validate it.

The manifest is an attestation from a trusted, controlled build, not a
cryptographic signature or proof obtained by decoding the BWA index. An
intentionally forged manifest is outside this accidental-substitution check.
Do not bless arbitrary customer indexes by writing new hashes into a JSON
file. Verify the build's provenance or regenerate the indexes.

The assembly accession/name is descriptive metadata. The FASTA and bundle
fingerprints identify the bytes; matching names or contig lengths do not
establish coordinate compatibility across assemblies.

## Verification

Python fixtures cover same-length sequence replacement, wrong/missing M5,
FAI offset change, index substitution, incomplete manifest, CRLF and receipt
roundtrip. Their index files are deliberately fake: these are byte-contract
tests, not mapping accuracy tests. The nf-test suite builds a real index with
the pinned BWA container and revalidates the resulting manifest, and tests
that a wrong M5 never reaches mapping/calling. Run it with:

```
nf-test test tests/modules/validate_reference_bundle.nf.test tests/pipeline/adzuki_snp_pipeline.nf.test --tag issue68_reference_bundle --profile test,docker
```

This change does not establish variant accuracy or a commercial SLA.
