#!/usr/bin/env python3
# sqcalc - structure factors from LAMMPS dump trajectories
# Copyright (C) 2026 Rui Su, Hangzhou Dianzi University
#
# SPDX-License-Identifier: GPL-3.0-or-later

r"""Check the markdown formulas for constructs that break MathJax rendering.

GitHub renders these documents with MathJax 3, which is not a LaTeX engine:

* the old font commands (``\\rm``, ``\\bf``, ``\\it``, ...) are gone, so a
  formula using them fails with "Undefined control sequence";
* ``\\[ ... \\]`` and ``\\( ... \\)`` are not delimiters GitHub recognises, so
  they are shown literally instead of being typeset;
* whatever sits inside a code span is literal text, so LaTeX there is shown
  verbatim rather than rendered.

It also flags a ``$$ ... $$`` block that holds more than one relation (two
equations separated by ``\qquad`` end up on a single line): split the block or
wrap the lines in ``\begin{aligned}``.

This is a development tool for this repository: run it after editing the
formulas in ``README.md`` or the skill references, or let the test suite do it.

Usage: check_readme_math.py [FILES...]

With no arguments the README and the skill references are checked.  Exits 0
when nothing suspicious is found and 1 otherwise, printing file and line.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
DEFAULT_FILES = [
    ROOT / "README.md",
    *sorted((ROOT / "skills").rglob("*.md")),
]

# Font commands that MathJax 3 dropped in favour of \mathrm, \mathbf, ...
OLD_FONT_COMMANDS = re.compile(r"\\(rm|bf|it|tt|sl|sf|cal|mit|mbf)\b")
# LaTeX display delimiters GitHub does not treat as math.
LATEX_DELIMITERS = re.compile(r"\\[\[\]()]")
# A backslash command, used to spot LaTeX that sits in a code span (where the
# markdown shows it as source instead of typesetting it).
LATEX_COMMAND = re.compile(r"\\[A-Za-z]+")
# Math environments MathJax 3 supports out of the box.
KNOWN_ENVIRONMENTS = {
    "align", "align*", "aligned", "array", "Bmatrix", "bmatrix", "cases",
    "gathered", "matrix", "pmatrix", "smallmatrix", "split", "subarray",
    "Vmatrix", "vmatrix",
}
ENVIRONMENT = re.compile(r"\\begin\{([^}]*)\}")
# A raw ``<`` or ``>`` directly in front of a letter is read as an HTML tag by
# the renderer (``\sum_{a<b}`` looks like ``<b>`` there): the math is split and
# MathJax reports "Extra open brace or missing close brace".  Write ``\lt`` /
# ``\gt``, or put spaces around the sign, so no tag can be recognised.
RAW_ANGLE = re.compile(r"[<>][A-Za-z]")


def angle_message(token: str) -> str:
    return ("%s looks like an HTML tag inside math; use \\lt / \\gt or spaces "
            "around the sign" % token)


def split_code_spans(line: str) -> list[tuple[str, bool]]:
    """Split a line into (text, is_code) segments.

    A table row contains several code spans, so the segments have to be taken
    in order rather than matched with a regular expression: the text between
    two spans is prose and does contain formulas.
    """
    segments = []
    position = 0
    while position < len(line):
        tick = line.find("`", position)
        if tick < 0:
            segments.append((line[position:], False))
            break
        close = line.find("`", tick + 1)
        if close < 0:                      # an unmatched backtick: treat as text
            segments.append((line[position:], False))
            break
        segments.append((line[position:tick], False))
        segments.append((line[tick + 1:close], True))
        position = close + 1
    return segments


def report(path: Path, line_number: int, message: str) -> None:
    try:
        name = path.resolve().relative_to(ROOT).as_posix()
    except ValueError:                      # a file outside the checkout
        name = str(path)
    print("%s:%d: %s" % (name, line_number, message))


def check_complaints(body: str) -> list[str]:
    """Complaints about the body of one display block (a $$ ... $$ pair)."""
    complaints = []
    if "\\begin{" not in body:
        relations = body.count("=")
        if relations > 1:
            complaints.append("display block holds %d relations; split it or use "
                              "\\begin{aligned}" % relations)
    complaints.extend("%s does not exist in MathJax 3, use \\mathrm{...} or "
                      "\\text{...}" % match.group(0)
                      for match in OLD_FONT_COMMANDS.finditer(body))
    complaints.extend("%s is not a GitHub math delimiter, use $$ or $" % match.group(0)
                      for match in LATEX_DELIMITERS.finditer(body))
    complaints.extend(angle_message(match.group(0)) for match in RAW_ANGLE.finditer(body))
    complaints.extend("unknown math environment %s" % match.group(0)
                      for match in ENVIRONMENT.finditer(body)
                      if match.group(1) not in KNOWN_ENVIRONMENTS)
    return complaints


def check(path: Path) -> int:
    """Return the number of complaints for one file."""
    complaints = 0
    dollars = 0
    in_fence = False
    in_display = False
    in_inline = False
    display_start = 0
    display_body: list[str] = []
    for number, line in enumerate(path.read_text().splitlines(), start=1):
        if line.lstrip().startswith("```"):
            in_fence = not in_fence
            continue
        if in_fence:
            # Shell examples legitimately contain $ (commands, substitution),
            # and markdown does not typeset anything inside a fence.
            continue

        stripped = line.strip()
        if not in_display and stripped.startswith("$$"):
            dollars += 2*stripped.count("$$")
            if stripped.count("$$") >= 2:          # $$ ... $$ on one line
                for message in check_complaints(stripped.replace("$$", "")):
                    report(path, number, message)
                    complaints += 1
                continue
            in_display, display_start, display_body = True, number, [stripped[2:]]
            continue
        if in_display:
            dollars += 2*stripped.count("$$")
            closing = stripped.endswith("$$")
            display_body.append(stripped[:-2] if closing else line)
            if closing:
                in_display = False
                for message in check_complaints("".join(display_body)):
                    report(path, display_start, message)
                    complaints += 1
            continue

        for text, is_code in split_code_spans(line):
            if is_code:
                if LATEX_COMMAND.search(text):
                    report(path, number, "LaTeX inside a code span is shown "
                           "verbatim: `%s`" % text)
                    complaints += 1
                continue
            dollars += text.count("$")
            # Inline math is scanned character by character, because that is
            # the only place a raw angle bracket turns into an HTML tag.
            position = 0
            while position < len(text):
                if text[position] == "$":
                    in_inline = not in_inline
                elif (in_inline and text[position] in "<>" and position + 1 < len(text)
                      and text[position + 1].isalpha()):
                    report(path, number, angle_message(text[position:position + 2]))
                    complaints += 1
                position += 1
            for match in OLD_FONT_COMMANDS.finditer(text):
                report(path, number, "%s does not exist in MathJax 3, use "
                       "\\mathrm{...} or \\text{...}" % match.group(0))
                complaints += 1
            for match in LATEX_DELIMITERS.finditer(text):
                report(path, number, "%s is not a GitHub math delimiter, use $$ "
                       "or $" % match.group(0))
                complaints += 1
            for match in ENVIRONMENT.finditer(text):
                if match.group(1) not in KNOWN_ENVIRONMENTS:
                    report(path, number, "unknown math environment %s" % match.group(0))
                    complaints += 1
    if dollars % 2:
        report(path, 0, "odd number of $ delimiters (%d): a formula is not closed"
               % dollars)
        complaints += 1
    if in_display:
        report(path, display_start, "display block is never closed with $$")
        complaints += 1
    return complaints


def main() -> int:
    files = [Path(argument) for argument in sys.argv[1:]] or DEFAULT_FILES
    complaints = 0
    checked = 0
    for path in files:
        if not path.is_file():
            print("error: %s not found" % path, file=sys.stderr)
            return 2
        complaints += check(path)
        checked += 1
    if complaints:
        print("%d suspicious formula(s) in %d file(s)" % (complaints, checked))
        return 1
    print("formulas render cleanly: %d file(s) checked" % checked)
    return 0


if __name__ == "__main__":
    sys.exit(main())
