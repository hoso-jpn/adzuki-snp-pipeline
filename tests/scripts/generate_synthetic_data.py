from __future__ import annotations

import gzip
import io
import random
from pathlib import Path


RANDOM_SEED = 20260812
CONTIG_LENGTH = 5000
READ_LENGTH = 100
FRAGMENT_LENGTH = 300


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

    reference_dir.mkdir(
        parents=True,
        exist_ok=True,
    )
    reads_dir.mkdir(
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

    read_groups = {
        "sample_a_L001": (
            "chrSynthetic1",
            [400, 900, 1400, 1900],
        ),
        "sample_a_L002": (
            "chrSynthetic1",
            [400, 900, 2500, 3000],
        ),
        "sample_b_L001": (
            "chrSynthetic2",
            [500, 1200, 2200, 3200],
        ),
    }

    for read_group_id, (
        contig_name,
        fragment_starts,
    ) in read_groups.items():
        read1, read2 = build_fastq(
            read_group_id,
            contigs[contig_name],
            fragment_starts,
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
