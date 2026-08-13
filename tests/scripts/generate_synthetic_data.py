from __future__ import annotations

import gzip
import io
import random
from pathlib import Path


RANDOM_SEED = 20260812
CONTIG_LENGTH = 5000
READ_LENGTH = 100
FRAGMENT_LENGTH = 300


def alternate_base(reference_base: str) -> str:
    substitutions = {
        "A": "C",
        "C": "G",
        "G": "T",
        "T": "A",
    }
    return substitutions[reference_base]


def make_sequence(
    random_generator: random.Random,
    length: int,
) -> str:
    return "".join(
        random_generator.choice("ACGT")
        for _ in range(length)
    )


def reverse_complement(sequence: str) -> str:
    translation = str.maketrans("ACGT", "TGCA")
    return sequence.translate(translation)[::-1]


def write_fasta(
    path: Path,
    contigs: dict[str, str],
) -> None:
    lines: list[str] = []

    for name, sequence in contigs.items():
        lines.append(f">{name}")
        lines.extend(
            sequence[index:index + 80]
            for index in range(
                0,
                len(sequence),
                80,
            )
        )

    path.write_text(
        "\n".join(lines) + "\n",
        encoding="utf-8",
    )


def build_fastq(
    read_group_id: str,
    contig_segments: list[tuple[str, str, list[int]]],
    contig_variants: dict[str, dict[int, str]],
) -> tuple[str, str]:
    read1_records: list[str] = []
    read2_records: list[str] = []

    for contig_name, contig_sequence, fragment_starts in contig_segments:
        variants = contig_variants.get(
            contig_name,
            {},
        )

        for number, start in enumerate(
            fragment_starts,
            start=1,
        ):
            fragment = contig_sequence[
                start:start + FRAGMENT_LENGTH
            ]

            if len(fragment) != FRAGMENT_LENGTH:
                raise ValueError(
                    "Fragment exceeds contig: "
                    f"{read_group_id}, {contig_name}, {start}"
                )

            fragment_bases = list(fragment)

            for position, alternate in variants.items():
                relative_position = position - start

                if 0 <= relative_position < FRAGMENT_LENGTH:
                    fragment_bases[relative_position] = alternate

            fragment = "".join(fragment_bases)
            read1 = fragment[:READ_LENGTH]
            read2 = reverse_complement(
                fragment[-READ_LENGTH:]
            )
            read_name = (
                f"{read_group_id}_{contig_name}_{number:03d}"
            )
            quality = "I" * READ_LENGTH

            read1_records.extend(
                [
                    f"@{read_name}/1",
                    read1,
                    "+",
                    quality,
                ]
            )
            read2_records.extend(
                [
                    f"@{read_name}/2",
                    read2,
                    "+",
                    quality,
                ]
            )

    return (
        "\n".join(read1_records) + "\n",
        "\n".join(read2_records) + "\n",
    )


def write_deterministic_gzip(
    path: Path,
    content: str,
) -> None:
    with path.open("wb") as raw_file:
        with gzip.GzipFile(
            filename="",
            mode="wb",
            fileobj=raw_file,
            mtime=0,
        ) as gzip_file:
            with io.TextIOWrapper(
                gzip_file,
                encoding="utf-8",
                newline="\n",
            ) as text_file:
                text_file.write(content)


def main() -> None:
    project_dir = Path(__file__).resolve().parents[2]
    reference_dir = (
        project_dir / "tests/data/reference"
    )
    reads_dir = project_dir / "tests/data/reads"
    variants_dir = project_dir / "tests/data/variants"

    reference_dir.mkdir(
        parents=True,
        exist_ok=True,
    )
    reads_dir.mkdir(
        parents=True,
        exist_ok=True,
    )
    variants_dir.mkdir(
        parents=True,
        exist_ok=True,
    )

    random_generator = random.Random(
        RANDOM_SEED
    )
    contigs = {
        "chrSynthetic1": make_sequence(
            random_generator,
            CONTIG_LENGTH,
        ),
        "chrSynthetic2": make_sequence(
            random_generator,
            CONTIG_LENGTH,
        ),
    }

    write_fasta(
        reference_dir / "synthetic.fa",
        contigs,
    )

    # alt_dp/ref_dp are the exact per-sample depths this pinned GATK
    # container (broadinstitute/gatk:4.6.2.0, digest-pinned) produces
    # for this fixture, confirmed by actually running the pipeline
    # end to end (generated- and prebuilt-reference-index paths
    # produced bit-identical cohort VCFs): GATK's realized
    # HaplotypeCaller/GenotypeGVCFs depth at a site can land one read
    # below a naive count of fragments whose read1 spans the
    # position, because it is derived from GATK's internal
    # reference-confidence-band accounting rather than a plain
    # alignment pileup. These are exact contract values, not a
    # fragment-count estimate or a lower bound: expected_variants.tsv
    # is the canonical source that tests/pipeline/*.nf.test reads
    # directly, so a real depth regression (from, for example, a
    # duplicate-marking or mapping change) must fail the test rather
    # than silently pass a `>=` floor.
    variant_positions = {
        "sample_a": (
            "chrSynthetic1",
            1500,
            8,
            5,
        ),
        "sample_b": (
            "chrSynthetic2",
            1600,
            6,
            8,
        ),
    }

    sample_ids = sorted(variant_positions.keys())

    # sample_variants is contig-aware (sample_id -> contig_name ->
    # {position: alternate}) so that a sample's own ALT is never
    # applied while building its reference-supporting reads on the
    # other sample's variant contig.
    sample_variants: dict[str, dict[str, dict[int, str]]] = {}
    variant_sites: dict[str, dict[str, object]] = {}

    for sample_id, (
        contig_name,
        position,
        alt_dp,
        ref_dp,
    ) in variant_positions.items():
        reference_base = contigs[contig_name][position]
        alternate = alternate_base(reference_base)

        sample_variants[sample_id] = {
            contig_name: {
                position: alternate,
            },
        }
        variant_sites[sample_id] = {
            "contig": contig_name,
            "position": position,
            "ref": reference_base,
            "alt": alternate,
            "alt_dp": alt_dp,
            "ref_dp": ref_dp,
        }

    # Each read group lists (contig_name, fragment_starts) segments.
    # sample_a_L001/L002 share library_id=lib_a and intentionally
    # reuse fragment starts 400/900 on chrSynthetic1 across lanes to
    # exercise cross-lane duplicate marking; their other chrSynthetic1
    # starts cover sample_a's own ALT at position 1500. Both lanes
    # also carry reference-supporting fragments on chrSynthetic2 so
    # sample_a has a confident 0/0 call at sample_b's ALT site.
    # sample_b_L001 keeps its own ALT-covering fragments on
    # chrSynthetic2 and gains reference-supporting fragments on
    # chrSynthetic1 for sample_a's ALT site.
    read_groups = {
        "sample_a_L001": (
            "sample_a",
            [
                (
                    "chrSynthetic1",
                    [400, 900, 1430, 1450, 1470, 1490],
                ),
                (
                    "chrSynthetic2",
                    [1530, 1550, 1570, 1590],
                ),
            ],
        ),
        "sample_a_L002": (
            "sample_a",
            [
                (
                    "chrSynthetic1",
                    [400, 900, 1440, 1460, 1480, 1500],
                ),
                (
                    "chrSynthetic2",
                    [1540, 1560, 1580, 1600],
                ),
            ],
        ),
        "sample_b_L001": (
            "sample_b",
            [
                (
                    "chrSynthetic2",
                    [1510, 1525, 1540, 1555, 1570, 1585],
                ),
                (
                    "chrSynthetic1",
                    [1410, 1425, 1440, 1455, 1470, 1485],
                ),
            ],
        ),
    }

    expected_variant_rows = [
        "\t".join(
            [
                "contig",
                "position",
                "ref",
                "alt",
                "alt_sample_id",
                "alt_genotype",
                "alt_sample_dp",
                "ref_sample_id",
                "ref_genotype",
                "ref_sample_dp",
                "site_dp",
                "ac",
                "an",
                "af",
            ]
        )
    ]

    for alt_sample_id in sample_ids:
        site = variant_sites[alt_sample_id]
        contig_name = site["contig"]
        position = site["position"]
        alt_dp = site["alt_dp"]
        ref_dp = site["ref_dp"]
        ref_sample_id = next(
            sample_id
            for sample_id in sample_ids
            if sample_id != alt_sample_id
        )

        expected_variant_rows.append(
            "\t".join(
                [
                    contig_name,
                    str(position + 1),
                    site["ref"],
                    site["alt"],
                    alt_sample_id,
                    "1/1",
                    str(alt_dp),
                    ref_sample_id,
                    "0/0",
                    str(ref_dp),
                    str(alt_dp + ref_dp),
                    "2",
                    "4",
                    "0.5",
                ]
            )
        )

    (
        variants_dir / "expected_variants.tsv"
    ).write_text(
        "\n".join(expected_variant_rows) + "\n",
        encoding="utf-8",
    )

    for read_group_id, (
        sample_id,
        fragment_specs,
    ) in read_groups.items():
        contig_segments = [
            (
                contig_name,
                contigs[contig_name],
                fragment_starts,
            )
            for contig_name, fragment_starts in fragment_specs
        ]
        read1, read2 = build_fastq(
            read_group_id,
            contig_segments,
            sample_variants[sample_id],
        )

        write_deterministic_gzip(
            reads_dir
            / f"{read_group_id}_R1.fastq.gz",
            read1,
        )
        write_deterministic_gzip(
            reads_dir
            / f"{read_group_id}_R2.fastq.gz",
            read2,
        )

    samplesheet = (
        "sample_id,read_group_id,fastq_1,fastq_2,"
        "library_id,platform,platform_unit\n"
        "sample_a,sample_a_L001,"
        "tests/data/reads/sample_a_L001_R1.fastq.gz,"
        "tests/data/reads/sample_a_L001_R2.fastq.gz,"
        "lib_a,ILLUMINA,flowcell1.L001.ATCACG\n"
        "sample_a,sample_a_L002,"
        "tests/data/reads/sample_a_L002_R1.fastq.gz,"
        "tests/data/reads/sample_a_L002_R2.fastq.gz,"
        "lib_a,ILLUMINA,flowcell1.L002.ATCACG\n"
        "sample_b,sample_b_L001,"
        "tests/data/reads/sample_b_L001_R1.fastq.gz,"
        "tests/data/reads/sample_b_L001_R2.fastq.gz,"
        "lib_b,ILLUMINA,flowcell1.L001.CGATGT\n"
    )

    (
        project_dir / "tests/data/samplesheet.csv"
    ).write_text(
        samplesheet,
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
