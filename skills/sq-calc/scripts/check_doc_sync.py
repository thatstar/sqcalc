#!/usr/bin/env python3
"""Check the sq-calc skill documents against the sqcalc sources.

The skill is distributed with the code, so its references are expected to stay
in step with what they describe:

* every command line flag handled by ``src/sqc_options.f90`` appears in
  ``references/usage.md``, and every flag documented there exists in the parser
* every user facing CMake option in ``CMakeLists.txt`` appears in
  ``references/building.md``, and every option documented there exists in the
  build system

Usage: check_doc_sync.py [checkout-root]

The checkout root defaults to the repository this skill lives in.  Exits 0 when
the documents are in sync and 1 otherwise, printing what drifted.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

SKILL_DIR = Path(__file__).resolve().parent.parent

# Flags: short or long option spellings in the option parser.
FLAG_IN_SOURCE = re.compile(r"'(-{1,2}[a-z][a-z0-9-]*)'")
FLAG_IN_DOC = re.compile(r"(?<![\w-])(-{1,2}[a-z][a-z0-9-]*)")

# CMake: `option(NAME ...)` and `set(NAME ... CACHE <type> ...)`.
CMAKE_OPTION = re.compile(r"^\s*option\(\s*([A-Za-z_][A-Za-z0-9_]*)", re.M)
CMAKE_CACHE_VAR = re.compile(
    r"^\s*set\(\s*([A-Za-z_][A-Za-z0-9_]*)\b[^)]*CACHE\s+\w+", re.M
)
CMAKE_IN_DOC = re.compile(r"(?<![\w-])-D([A-Za-z_][A-Za-z0-9_]*)")

# Only the knobs a user is meant to pass, not CMake's own bookkeeping or the
# FINUFFT-internal switches the build sets itself.
USER_FACING = re.compile(r"^(SQC_|CPM_|FFTW_ROOT$|CMAKE_BUILD_TYPE$)")
# Options a build accepts without sqcalc's CMakeLists.txt defining them: CPM's
# <dependency>_SOURCE overrides are read by the vendored FINUFFT.
EXTRA_CMAKE_VARS = {"CPM_CCCL_SOURCE"}


def parser_flags(text: str) -> set[str]:
    """Command line flags accepted by the option parser."""
    flags: set[str] = set()
    for line in text.splitlines():
        if "case (" in line:
            flags |= set(FLAG_IN_SOURCE.findall(line))
    return flags


def cmake_options(text: str) -> set[str]:
    """User facing CMake options and cache variables."""
    names = set(CMAKE_OPTION.findall(text)) | set(CMAKE_CACHE_VAR.findall(text))
    return {name for name in names if USER_FACING.match(name)}


def documented_flags(text: str) -> set[str]:
    return set(FLAG_IN_DOC.findall(text))


def documented_cmake_vars(text: str) -> set[str]:
    return set(CMAKE_IN_DOC.findall(text))


def report(label: str, documented: set[str], actual: set[str], path: Path) -> bool:
    """Print the drift in one direction; return True when nothing drifted."""
    missing = sorted(actual - documented)
    extra = sorted(documented - actual)
    if not missing and not extra:
        return True
    print(f"{label}: {path} is out of sync")
    if missing:
        print(f"  missing from {path.name}: {', '.join(missing)}")
    if extra:
        print(f"  no longer in the sources: {', '.join(extra)}")
    return False


def main() -> int:
    root = Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else SKILL_DIR.parent.parent
    parser = root / "src" / "sqc_options.f90"
    cmake = root / "CMakeLists.txt"
    usage = SKILL_DIR / "references" / "usage.md"
    building = SKILL_DIR / "references" / "building.md"

    for path in (parser, cmake, usage, building):
        if not path.is_file():
            print(f"error: {path} not found", file=sys.stderr)
            return 2

    flags = parser_flags(parser.read_text())
    docs_flags = documented_flags(usage.read_text())
    options = cmake_options(cmake.read_text()) | EXTRA_CMAKE_VARS
    docs_cmake = documented_cmake_vars(building.read_text())

    ok = report("command line flags", docs_flags, flags, usage)
    ok &= report("CMake options", docs_cmake, options, building)

    if not ok:
        return 1
    print(
        f"in sync: {len(flags)} command line flags and {len(options)} CMake options"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
