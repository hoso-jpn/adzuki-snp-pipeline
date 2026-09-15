"""Region sets and strata for the Issue #65 quality-evidence evaluation.

Coordinates follow BED throughout: 0-based, half-open `[start, end)`. VCF
positions are 1-based and are converted at the one place that reads them
(`RegionSet.contains_span`). Every region set is bound to one reference by its
ordered contig names and lengths, and any interval outside that reference is a
hard error rather than something silently clipped.

Two kinds of strata exist and are kept apart on purpose:

* a **partition** is a group of region sets that are pairwise disjoint and
  together cover a stated universe exactly (for example GC bins or depth bins
  over the callable universe). Each position belongs to exactly one member, so
  per-member counts add up to the universe's count.
* a **tag** is a single region set that may overlap any other tag (repeat,
  low-mappability, callable, ...). A variant can carry several tags at once,
  so tag counts are not additive and are never presented as if they were.

`definition_hash` identifies a region set by the reference it is bound to, its
merged intervals, and how it was made, so the same inputs always give the same
hash and any change to membership changes it.
"""

from __future__ import annotations

import bisect
import hashlib
import json
from dataclasses import dataclass, field
from pathlib import Path


class MalformedRegionError(ValueError):
    """A region input is not a valid set of intervals on the declared reference."""


def read_fai(path: Path) -> tuple[tuple[str, int], ...]:
    contigs: list[tuple[str, int]] = []
    seen: set[str] = set()
    for number, line in enumerate(Path(path).read_text(encoding="utf-8").splitlines(), start=1):
        if not line:
            continue
        fields = line.split("\t")
        if len(fields) < 2 or not fields[1].isdigit() or fields[0] in seen:
            raise MalformedRegionError(f"{path}: line {number} is not a valid FAI entry")
        seen.add(fields[0])
        contigs.append((fields[0], int(fields[1])))
    if not contigs:
        raise MalformedRegionError(f"{path}: no contigs")
    return tuple(contigs)


def _merge(intervals: list[tuple[int, int]]) -> list[tuple[int, int]]:
    merged: list[tuple[int, int]] = []
    for start, end in sorted(intervals):
        if merged and start <= merged[-1][1]:
            if end > merged[-1][1]:
                merged[-1] = (merged[-1][0], end)
        else:
            merged.append((start, end))
    return merged


@dataclass
class RegionSet:
    """Merged, sorted intervals per contig on one reference."""

    name: str
    contigs: tuple[tuple[str, int], ...]
    intervals: dict[str, list[tuple[int, int]]] = field(default_factory=dict)
    provenance: dict[str, object] = field(default_factory=dict)

    def __post_init__(self) -> None:
        lengths = dict(self.contigs)
        normalized: dict[str, list[tuple[int, int]]] = {}
        for contig, spans in self.intervals.items():
            if contig not in lengths:
                raise MalformedRegionError(
                    f"{self.name}: contig {contig!r} is not in the reference"
                )
            for start, end in spans:
                if not (0 <= start < end <= lengths[contig]):
                    raise MalformedRegionError(
                        f"{self.name}: interval {contig}:{start}-{end} is outside 0..{lengths[contig]}"
                    )
            merged = _merge(list(spans))
            if merged:
                normalized[contig] = merged
        self.intervals = normalized
        self._starts = {contig: [s for s, _ in spans] for contig, spans in normalized.items()}

    # -- construction ------------------------------------------------------

    @classmethod
    def from_bed(cls, path: Path, name: str, contigs: tuple[tuple[str, int], ...]) -> RegionSet:
        intervals: dict[str, list[tuple[int, int]]] = {}
        for number, line in enumerate(Path(path).read_text(encoding="utf-8").splitlines(), start=1):
            if not line or line.startswith(("#", "track", "browser")):
                continue
            fields = line.split("\t")
            if len(fields) < 3 or not fields[1].isdigit() or not fields[2].isdigit():
                raise MalformedRegionError(f"{path}: line {number} is not a BED interval")
            intervals.setdefault(fields[0], []).append((int(fields[1]), int(fields[2])))
        return cls(
            name=name, contigs=contigs, intervals=intervals, provenance={"bed": Path(path).name}
        )

    @classmethod
    def whole_reference(cls, name: str, contigs: tuple[tuple[str, int], ...]) -> RegionSet:
        return cls(name=name, contigs=contigs, intervals={c: [(0, n)] for c, n in contigs})

    # -- queries -----------------------------------------------------------

    def bases(self) -> int:
        return sum(end - start for spans in self.intervals.values() for start, end in spans)

    def _index_of(self, contig: str, position0: int) -> int | None:
        starts = self._starts.get(contig)
        if not starts:
            return None
        i = bisect.bisect_right(starts, position0) - 1
        if i >= 0 and position0 < self.intervals[contig][i][1]:
            return i
        return None

    def contains_span(self, contig: str, pos1: int, length: int) -> str:
        """Where a 1-based span of `length` bases sits: `inside`, `outside` or `boundary`.

        `boundary` means the span touches the region but is not wholly inside
        one merged interval; comparisons exclude such spans with their own
        reason rather than counting them on either side of the boundary.
        """
        start0, end0 = pos1 - 1, pos1 - 1 + max(length, 1)
        first = self._index_of(contig, start0)
        if first is not None and end0 <= self.intervals[contig][first][1]:
            return "inside"
        spans = self.intervals.get(contig, [])
        i = bisect.bisect_left(self._starts.get(contig, []), end0)
        for start, end in spans[max(0, i - 1) : i + 1]:
            if start < end0 and start0 < end:
                return "boundary"
        return "outside"

    # -- set algebra -------------------------------------------------------

    def _check_same_reference(self, other: RegionSet) -> None:
        if self.contigs != other.contigs:
            raise MalformedRegionError(f"{self.name} and {other.name} are on different references")

    def union(self, other: RegionSet, name: str) -> RegionSet:
        self._check_same_reference(other)
        out = {
            c: list(self.intervals.get(c, [])) + list(other.intervals.get(c, []))
            for c, _ in self.contigs
        }
        return RegionSet(name, self.contigs, out)

    def intersect(self, other: RegionSet, name: str) -> RegionSet:
        self._check_same_reference(other)
        out: dict[str, list[tuple[int, int]]] = {}
        for contig, _ in self.contigs:
            a, b = self.intervals.get(contig, []), other.intervals.get(contig, [])
            i = j = 0
            spans = []
            while i < len(a) and j < len(b):
                start, end = max(a[i][0], b[j][0]), min(a[i][1], b[j][1])
                if start < end:
                    spans.append((start, end))
                if a[i][1] < b[j][1]:
                    i += 1
                else:
                    j += 1
            out[contig] = spans
        return RegionSet(name, self.contigs, out)

    def complement(self, name: str) -> RegionSet:
        out: dict[str, list[tuple[int, int]]] = {}
        for contig, length in self.contigs:
            cursor, spans = 0, []
            for start, end in self.intervals.get(contig, []):
                if start > cursor:
                    spans.append((cursor, start))
                cursor = end
            if cursor < length:
                spans.append((cursor, length))
            out[contig] = spans
        return RegionSet(name, self.contigs, out)

    def subtract(self, other: RegionSet, name: str) -> RegionSet:
        return self.intersect(other.complement(f"not_{other.name}"), name)

    # -- identity ----------------------------------------------------------

    def definition(self, reference_sha256: str) -> dict[str, object]:
        return {
            "coordinate_system": "BED 0-based half-open",
            "reference_sha256": reference_sha256,
            "contigs": [[c, n] for c, n in self.contigs],
            "intervals": [[c, s, e] for c, _ in self.contigs for s, e in self.intervals.get(c, [])],
            "provenance": self.provenance,
        }

    def definition_hash(self, reference_sha256: str) -> str:
        text = json.dumps(self.definition(reference_sha256), sort_keys=True, separators=(",", ":"))
        return "sha256:" + hashlib.sha256(text.encode("utf-8")).hexdigest()

    def to_bed(self) -> str:
        return "".join(
            f"{c}\t{s}\t{e}\n" for c, _ in self.contigs for s, e in self.intervals.get(c, [])
        )


@dataclass(frozen=True)
class Stratum:
    name: str
    kind: str  # "partition" or "tag"
    group: str  # partition group name; for tags, the tag family
    region: RegionSet


def validate_partition(members: list[RegionSet], universe: RegionSet) -> None:
    """A partition's members must be pairwise disjoint and cover the universe exactly."""
    covered = 0
    for i, a in enumerate(members):
        a._check_same_reference(universe)
        for b in members[i + 1 :]:
            if a.intersect(b, "overlap").bases():
                raise MalformedRegionError(f"partition members {a.name} and {b.name} overlap")
        if a.subtract(universe, "outside").bases():
            raise MalformedRegionError(f"partition member {a.name} extends outside {universe.name}")
        covered += a.bases()
    if covered != universe.bases():
        raise MalformedRegionError(
            f"partition covers {covered} of {universe.bases()} bases of {universe.name}"
        )
