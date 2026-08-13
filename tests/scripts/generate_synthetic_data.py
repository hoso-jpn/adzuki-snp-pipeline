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
    contig_sequence: str,
    fragment_starts: list[int],
    variants: dict[int, str],
) -> tuple[str, str]:
    read1_records: list[str] = []
    read2_records: list[str] = []

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
                f"{read_group_id}, {start}"
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
            f"{read_group_id}_{number:03d}"
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

    variant_positions = {
        "sample_a": (
            "chrSynthetic1",
            1500,
        ),
        "sample_b": (
            "chrSynthetic2",
            1600,
        ),
    }

    sample_variants: dict[str, dict[int, str]] = {}
    expected_variant_rows = [
        "sample_id\tcontig\tposition\tref\talt"
    ]

    for sample_id, (
        contig_name,
        position,
    ) in variant_positions.items():
        reference_base = contigs[contig_name][position]
        alternate = alternate_base(reference_base)

        sample_variants[sample_id] = {
            position: alternate,
        }
        expected_variant_rows.append(
            "\t".join(
                [
                    sample_id,
                    contig_name,
                    str(position + 1),
                    reference_base,
                    alternate,
                ]
            )
        )

    (
        variants_dir / "expected_variants.tsv"
    ).write_text(
        "\n".join(expected_variant_rows) + "\n",
        encoding="utf-8",
    )

    read_groups = {
        "sample_a_L001": (
            "sample_a",
            "chrSynthetic1",
            [400, 900, 1430, 1450, 1470, 1490],
        ),
        "sample_a_L002": (
            "sample_a",
            "chrSynthetic1",
            [400, 900, 1440, 1460, 1480, 1500],
        ),
        "sample_b_L001": (
            "sample_b",
            "chrSynthetic2",
            [1510, 1525, 1540, 1555, 1570, 1585],
        ),
    }

    for read_group_id, (
        sample_id,
        contig_name,
        fragment_starts,
    ) in read_groups.items():
        read1, read2 = build_fastq(
            read_group_id,
            contigs[contig_name],
            fragment_starts,
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
