#!/usr/bin/env python3
# sqcalc - structure factors from LAMMPS dump trajectories
# Copyright (C) 2026 Rui Su, Hangzhou Dianzi University
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""Check the sq-calc skill documents against the sqcalc sources.

This is a development tool for this repository, not part of the distributed
skill: run it after touching the option parser or the build system to keep the
skill references in step with what they describe.

* every command line flag handled by ``src/sqc_options.f90`` appears in
  ``skills/sq-calc/references/usage.md``, and every flag documented there
  exists in the parser
* every flag is documented under the subcommand that accepts it: the
  ``flag_scope`` function of the parser says whether an option is global or
  belongs to ``static`` or ``dyn``, and the option tables of ``usage.md`` and
  ``README.md`` have to agree with it
* the built-in ``--help`` lists exactly the options of each scope, so a flag
  cannot be forgotten in ``print_usage`` while the documents know about it
* every user facing CMake option in ``CMakeLists.txt`` appears in
  ``skills/sq-calc/references/building.md``, and every option documented there
  exists in the build system

Usage: check_doc_sync.py [checkout-root]

The checkout root defaults to the repository this script lives in.  Exits 0
when the documents are in sync and 1 otherwise, printing what drifted.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent

# Flags: short or long option spellings in the option parser.
FLAG_IN_SOURCE = re.compile(r"'(-{1,2}[a-z][a-z0-9-]*)'")
FLAG_IN_DOC = re.compile(r"(?<![\w-])(-{1,2}[a-z][a-z0-9-]*)")

# Which subcommand accepts an option.  `both` is the global set.
SCOPE_BOTH = "both"
SCOPE_STATIC = "static"
SCOPE_DYN = "dyn"

# The parser states the scope of every flag that is not global in flag_scope();
# a flag it does not mention is accepted by both subcommands.
FLAG_SCOPE_FUNC = re.compile(r"function\s+flag_scope\b(.*?)\n\s*end function", re.S)
FLAG_SCOPE_GROUP = re.compile(r"case\s*\((.*?)\)\s*scope\s*=\s*scope_(\w+)", re.S)
SCOPE_NAMES = {"static": SCOPE_STATIC, "dyn": SCOPE_DYN, "both": SCOPE_BOTH}

# The option tables of the two documents and the section heading or lead-in
# that introduces each of them.
USAGE_SECTIONS = (
    (re.compile(r"global options\s*$", re.I), SCOPE_BOTH),
    (re.compile(r"sqcalc static.*options\s*$", re.I), SCOPE_STATIC),
    (re.compile(r"sqcalc dyn.*options\s*$", re.I), SCOPE_DYN),
)
README_SECTIONS = (
    (re.compile(r"global to both subcommands", re.I), SCOPE_BOTH),
    (re.compile(r"belong to `sqcalc static`", re.I), SCOPE_STATIC),
    (re.compile(r"pair histogram", re.I), SCOPE_STATIC),
    (re.compile(r"belong to `sqcalc dyn`", re.I), SCOPE_DYN),
)

# The built-in --help is assembled by hand in print_usage(); each heading there
# introduces the options of one scope.
PRINT_USAGE = re.compile(r"subroutine\s+print_usage\b(.*?)\n\s*end subroutine", re.S | re.I)
HELP_HEADINGS = (
    (re.compile(r"global options:"), SCOPE_BOTH),
    (re.compile(r"static options:"), SCOPE_STATIC),
    (re.compile(r"dyn options:"), SCOPE_DYN),
)
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


def source_scopes(text: str) -> dict[str, str] | None:
    """Scope of every non-global flag, from flag_scope().

    Returns None when the function cannot be read, so that a rework of the
    parser turns the scope check off instead of breaking the build.
    """
    func = FLAG_SCOPE_FUNC.search(text)
    if not func:
        return None
    scopes: dict[str, str] = {}
    groups = FLAG_SCOPE_GROUP.findall(func.group(1))
    if not groups:
        return None
    for flags, scope in groups:
        name = SCOPE_NAMES.get(scope)
        if name is None:
            return None
        for flag in FLAG_IN_SOURCE.findall(flags):
            scopes[flag] = name
    return scopes


def heading_sections(text: str) -> list[tuple[str, str]]:
    """(heading, body) of every level two or deeper markdown section."""
    parts = re.split(r"^#{2,}\s+(.*)$", text, flags=re.M)
    return [(parts[i].strip(), parts[i + 1]) for i in range(1, len(parts) - 1, 2)]


def leadin_sections(text: str) -> list[tuple[str, str]]:
    """(lead-in paragraph, following table) of the tables of a document."""
    blocks = [block for block in re.split(r"\n\s*\n", text) if block.strip()]
    return [(blocks[i], blocks[i + 1])
            for i in range(len(blocks) - 1)
            if blocks[i + 1].lstrip().startswith("|")]


def documented_sections(text: str, patterns, headings: bool) -> list[tuple[str, str, str]]:
    """(key, table text, claimed scope) of the option tables of one document."""
    found: list[tuple[str, str, str]] = []
    for key, body in heading_sections(text) if headings else leadin_sections(text):
        for pattern, scope in patterns:
            if pattern.search(key):
                found.append((key.strip(), body, scope))
                break
    return found


def check_scope(scopes: dict[str, str], sections, path: Path) -> bool:
    """Do the documented option tables agree with the scope of the parser?"""
    problems: list[str] = []
    exact = {SCOPE_BOTH: set(), SCOPE_STATIC: set(), SCOPE_DYN: set()}
    for key, body, claimed in sections:
        for flag in sorted(set(FLAG_IN_DOC.findall(body))):
            actual = scopes.get(flag)
            if actual is None:
                continue  # unknown flags are the flag set check's business
            fits = actual == claimed if claimed == SCOPE_BOTH else \
                actual in (claimed, SCOPE_BOTH)
            if not fits:
                problems.append(
                    f'  {flag} belongs to {actual}, but "{key}" lists it as {claimed}'
                )
            if actual == claimed:
                exact[claimed].add(flag)
    for flag, actual in sorted(scopes.items()):
        if flag not in exact[actual]:
            problems.append(
                f"  {flag} belongs to {actual}, but no {actual} table documents it"
            )
    if not problems:
        return True
    print(f"option scopes: {path} is out of sync")
    for line in problems:
        print(line)
    return False


def help_blocks(text: str) -> dict[str, set[str]] | None:
    """Flags listed under each heading of the built-in --help.

    Returns None when print_usage cannot be read, so that a rework of the help
    turns this check off instead of breaking the build.
    """
    func = PRINT_USAGE.search(text)
    if not func:
        return None
    blocks: dict[str, set[str]] = {}
    current: str | None = None
    for line in func.group(1).splitlines():
        for pattern, scope in HELP_HEADINGS:
            if pattern.search(line):
                current = scope
                blocks.setdefault(scope, set())
                break
        if current is not None:
            blocks[current] |= set(FLAG_IN_DOC.findall(line))
    return blocks or None


def check_help(scopes: dict[str, str], blocks: dict[str, set[str]], path: Path) -> bool:
    """Does the built-in --help list exactly the options of each scope?"""
    problems: list[str] = []
    for scope in (SCOPE_BOTH, SCOPE_STATIC, SCOPE_DYN):
        shown = blocks.get(scope)
        if shown is None:
            problems.append(f'  print_usage has no "{scope} options:" block')
            continue
        expected = {flag for flag, actual in scopes.items() if actual == scope}
        missing = sorted(expected - shown)
        extra = sorted(shown - expected)
        if missing:
            problems.append(
                f"  {scope} options missing from the built-in help: {', '.join(missing)}"
            )
        if extra:
            problems.append(
                f"  listed as {scope} options in the built-in help but not accepted "
                f"there: {', '.join(extra)}"
            )
    if not problems:
        return True
    print(f"built-in help: {path} is out of sync")
    for line in problems:
        print(line)
    return False


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
    root = Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else ROOT
    skill = root / "skills" / "sq-calc"
    parser = root / "src" / "sqc_options.f90"
    cmake = root / "CMakeLists.txt"
    usage = skill / "references" / "usage.md"
    building = skill / "references" / "building.md"
    readme = root / "README.md"

    for path in (parser, cmake, usage, building, readme):
        if not path.is_file():
            print(f"error: {path} not found", file=sys.stderr)
            return 2

    parser_text = parser.read_text()
    flags = parser_flags(parser_text)
    docs_flags = documented_flags(usage.read_text())
    options = cmake_options(cmake.read_text()) | EXTRA_CMAKE_VARS
    docs_cmake = documented_cmake_vars(building.read_text())

    ok = report("command line flags", docs_flags, flags, usage)
    ok &= report("CMake options", docs_cmake, options, building)

    stated = source_scopes(parser_text)
    if stated is None:
        print("note: flag_scope() not found in the parser, skipping the scope check")
    else:
        scopes = {flag: stated.get(flag, SCOPE_BOTH) for flag in flags}
        for path, patterns, headings in ((usage, USAGE_SECTIONS, True),
                                         (readme, README_SECTIONS, False)):
            sections = documented_sections(path.read_text(), patterns, headings)
            if not sections:
                print(f"note: no option table found in {path.name}, "
                      "skipping its scope check")
                continue
            ok &= check_scope(scopes, sections, path)
        blocks = help_blocks(parser_text)
        if blocks is None:
            print("note: print_usage() not found in the parser, "
                  "skipping the built-in help check")
        else:
            ok &= check_help(scopes, blocks, parser)

    if not ok:
        return 1
    print(
        f"in sync: {len(flags)} command line flags and {len(options)} CMake options"
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
