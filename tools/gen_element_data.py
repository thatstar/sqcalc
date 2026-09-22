#!/usr/bin/env python3
# sqcalc - structure factors from LAMMPS dump trajectories
# Copyright (C) 2026 Rui Su, Hangzhou Dianzi University
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""Generate src/sqc_element_data.f90 from the debyer element tables.

The tabulated numbers (masses, IT92 X-ray coefficients, NN92 bound coherent
neutron scattering lengths) are physical data compiled by the debyer project;
see the references below.  The C source itself is not redistributed here, only
the derived Fortran data table that this script rewrites mechanically.

The element table carries one row per element (mass, neutron length).  The
X-ray form factors are stored per species instead: the neutral atom of every
element plus the ions and valence states of the IT92 table, so that
`-m 1:Si4+,2:O2-` picks a different row than `-m 1:Si,2:O`.  Source rows that
merely repeat another row under a different name ("H'" beside "D", "Siv"
beside "Si") become aliases of that row.  Every emitted species is checked
against the electron-count sum rule f(0) = Z - q before it is written.

Usage:
    tools/gen_element_data.py /path/to/debyer/debyer/atomtables.c src/sqc_element_data.f90

References
----------
* masses: CRC Handbook of Chemistry and Physics (63rd/70th ed.)
* X-ray:  International Tables for Crystallography, Vol. C (1992),
          Table 6.1.1.4 - analytical approximation to the scattering factors
          ("IT92" / Cromer-Mann a_i, b_i, c).
* neutron: Neutron News 3 (1992) 29-37, NIST table of bound coherent
          scattering lengths ("NN92").
"""

import hashlib
import re
import sys

SPECIES_LEN = 8  # "Sival" is the longest label
F0_TOL = 0.1  # electron-count check; the rounded Tl and Pb rows reach 0.06

# Element symbol plus an optional charge ("O2-", "Na1+") or valence suffix.
SPECIES_RE = re.compile(r"^([A-Z][a-z]?)(?:(\d+)([+-])|([+-])|(val))?$")


def strip_comments(text):
    return re.sub(r"/\*.*?\*/", "", text, flags=re.S)


def grab_block(text, name):
    start = text.index(name)
    start = text.index("{", start)
    depth = 0
    for idx in range(start, len(text)):
        if text[idx] == "{":
            depth += 1
        elif text[idx] == "}":
            depth -= 1
            if depth == 0:
                return text[start : idx + 1]
    raise ValueError("unterminated block " + name)


def parse_pse(block):
    rows = {}
    for line in block.splitlines():
        m = re.match(r'\s*\{\s*(\d+)\s*,\s*"([^"]*)"\s*,\s*"([^"]*)"\s*,\s*([0-9.]+)\s*\}', line)
        if m:
            rows[m.group(2)] = (int(m.group(1)), m.group(3), float(m.group(4)))
    return rows


def parse_it92(block):
    rows = {}
    # rows look like: { "H", { a,a,a,a }, { b,b,b,b }, c },
    pattern = re.compile(
        r'\{\s*"([^"]*)"\s*,\s*\{([^}]*)\}\s*,\s*\{([^}]*)\}\s*,\s*([-0-9.eE+]+)\s*\}'
    )
    for m in pattern.finditer(block):
        a = [float(x) for x in m.group(2).split(",")]
        b = [float(x) for x in m.group(3).split(",")]
        rows[m.group(1)] = (a, b, float(m.group(4)))
    return rows


def parse_nn92(block):
    rows = {}
    pattern = re.compile(
        r'\{\s*"([^"]*)"\s*,\s*([-0-9.eE+]+)\s*,\s*([-0-9.eE+]+)\s*,\s*([-0-9.eE+]+)\s*\}'
    )
    for m in pattern.finditer(block):
        rows[m.group(1)] = float(m.group(2))
    return rows


def fortran_real(value):
    text = repr(float(value))
    if "e" in text or "E" in text:
        mantissa, exponent = re.split("[eE]", text)
        return "%sd%s" % (mantissa, int(exponent))
    return text + "d0"


def split_species_label(label):
    """Split "O2-", "Si4+", "Cval" or "Na" into (element, charge, valence).

    Returns None for labels that are not an element plus suffix, such as the
    two duplicate rows of the source table ("H'" and "Siv").
    """
    m = SPECIES_RE.match(label)
    if not m:
        return None
    element, digits, sign, shorthand, valence = m.groups()
    if valence:
        return element, 0, True
    if sign:
        return element, int(digits) * (1 if sign == "+" else -1), False
    if shorthand:
        return element, 1 if shorthand == "+" else -1, False
    return element, 0, False


def species_name(element, charge, valence):
    """Canonical label of the species an element, charge and flag describe."""
    if valence:
        return element + "val"
    if charge == 0:
        return element
    return "%s%d%s" % (element, abs(charge), "+" if charge > 0 else "-")


def same_coefficients(first, second, tol=1e-12):
    a1, b1, c1 = first
    a2, b2, c2 = second
    return (
        all(abs(x - y) < tol for x, y in zip(a1, a2))
        and all(abs(x - y) < tol for x, y in zip(b1, b2))
        and abs(c1 - c2) < tol
    )


def build_species(symbols, it92):
    """Return (species, aliases) from the element list and the IT92 table.

    Species is a list of (label, element index, charge, valence, a, b, c) with
    the rows of one element adjacent: the neutral atom, then the valence
    states, then the anions and cations by ascending charge magnitude, which is
    the order the source table uses.  Aliases map a duplicate source label onto
    the canonical species that carries the same coefficients.
    """
    index_of = {symbol: i + 1 for i, symbol in enumerate(symbols)}
    species = []
    used = set()
    for symbol in symbols:
        if symbol in it92:
            a, b, c = it92[symbol]
            species.append((symbol, index_of[symbol], 0, False, a, b, c))
            used.add(symbol)
        rows = []
        for label, (a, b, c) in it92.items():
            if label in used:
                continue
            parsed = split_species_label(label)
            if parsed is None:
                continue
            element, charge, valence = parsed
            if element != symbol or (charge == 0 and not valence):
                continue
            rows.append((charge, label, valence, a, b, c))
        def order(row):
            charge, label, valence = row[0], row[1], row[2]
            if valence:
                kind = 1
            elif charge < 0:
                kind = 2
            else:
                kind = 3
            return (kind, abs(charge), label)

        rows.sort(key=order)
        for charge, label, valence, a, b, c in rows:
            # The Fortran lookup rebuilds the label from the element, charge
            # and valence flag, so a source row spelled any other way would be
            # emitted as a state that no name can select.
            canonical = species_name(symbol, charge, valence)
            if label != canonical:
                raise SystemExit(
                    "gen_element_data.py: source label %r is not the canonical %r"
                    % (label, canonical)
                )
            species.append((label, index_of[symbol], charge, valence, a, b, c))
            used.add(label)

    aliases = []
    for label, coefficients in it92.items():
        if label in used:
            continue
        twins = [row for row in species if same_coefficients(coefficients, row[4:7])]
        if not twins:
            raise SystemExit("gen_element_data.py: no canonical row for %r" % label)
        aliases.append((label, twins[0][0]))
    aliases.sort()
    return species, aliases


def audit_species(species, symbols, pse):
    """Check f(0) = Z - q for every row; return the largest deviation."""
    worst = 0.0
    for label, element, charge, valence, a, b, c in species:
        z = pse[symbols[element - 1]][0]
        expected = z if valence else z - charge
        deviation = abs(sum(a) + c - expected)
        worst = max(worst, deviation)
        if deviation > F0_TOL:
            raise SystemExit(
                "gen_element_data.py: %s has f(0) = %.4f, expected %.1f"
                % (label, sum(a) + c, expected)
            )
    return worst


def main():
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    source = open(sys.argv[1], "rb").read()
    digest = hashlib.sha256(source).hexdigest()
    text = strip_comments(source.decode())
    out_path = sys.argv[2]

    pse = parse_pse(grab_block(text, "pse_table[]"))
    it92 = parse_it92(grab_block(text, "it92_table[]"))
    nn92 = parse_nn92(grab_block(text, "nn92_table[]"))

    symbols = sorted(pse, key=lambda s: (pse[s][0], s))
    species, aliases = build_species(symbols, it92)
    worst = audit_species(species, symbols, pse)

    neutral_species = [0] * len(symbols)
    for index, row in enumerate(species, start=1):
        if row[2] == 0 and not row[3]:
            neutral_species[row[1] - 1] = index

    lines = []
    add = lines.append
    add("! sqcalc - structure factors from LAMMPS dump trajectories")
    add("! Copyright (C) 2026 Rui Su, Hangzhou Dianzi University")
    add("!")
    add("! SPDX-License-Identifier: GPL-3.0-or-later")
    add("!")
    add("! Generated by tools/gen_element_data.py -- do not edit by hand.")
    add("!")
    add("! Sources of the tabulated data (via the debyer project):")
    add("!   masses  : CRC Handbook of Chemistry and Physics, 63rd/70th ed.")
    add("!   X-ray   : International Tables for Crystallography C (1992), 6.1.1.4")
    add("!   neutron : Neutron News 3 (1992) 29-37 (NIST bound coherent lengths)")
    add("!")
    add("! atomtables.c sha256 %s" % digest)
    add("! electron-count audit: max |f(0) - (Z - q)| = %.4f" % worst)
    add("module sqc_element_data")
    add("   use sqc_kinds, only: rk")
    add("   implicit none")
    add("   private")
    add("")
    add("   integer, parameter, public :: n_elements = %d" % len(symbols))
    add("   integer, parameter, public :: n_species = %d" % len(species))
    add("   integer, parameter, public :: n_species_alias = %d" % len(aliases))
    add("")
    add("   ! Atomic number, symbol, name and atomic mass of each entry.")
    add("   integer, parameter, public :: element_z(n_elements) = [ &")
    add(_wrapped(["%d" % pse[s][0] for s in symbols], "      "))
    add("   character(len=2), parameter, public :: element_symbol(n_elements) = [ &")
    add(_wrapped(['"%-2s"' % s for s in symbols], "      "))
    add("   character(len=16), parameter, public :: element_name(n_elements) = [ &")
    add(_wrapped(['"%-16s"' % pse[s][1] for s in symbols], "      "))
    add("   real(rk), parameter, public :: element_mass(n_elements) = [ &")
    add(_wrapped([fortran_real(pse[s][2]) for s in symbols], "      "))
    add("")
    add("   ! Bound coherent neutron scattering length [fm]; 0 if not tabulated.")
    add("   real(rk), parameter, public :: element_neutron_b(n_elements) = [ &")
    add(_wrapped([fortran_real(nn92.get(s, 0.0)) for s in symbols], "      "))
    add("   logical, parameter, public :: element_has_neutron(n_elements) = [ &")
    add(_wrapped(["%s" % (".true." if s in nn92 else ".false.") for s in symbols], "      "))
    add("")
    add("   ! X-ray species: the neutral atom of every element plus the ions")
    add("   ! and valence states of the IT92 table.  species_element indexes the")
    add("   ! element table above, species_charge is the formal charge (0 for")
    add("   ! neutral and valence states) and species_valence marks the latter.")
    add("   character(len=%d), parameter, public :: species_label(n_species) = [ &" % SPECIES_LEN)
    add(_wrapped(['"%-*s"' % (SPECIES_LEN, row[0]) for row in species], "      "))
    add("   integer, parameter, public :: species_element(n_species) = [ &")
    add(_wrapped(["%d" % row[1] for row in species], "      "))
    add("   integer, parameter, public :: species_charge(n_species) = [ &")
    add(_wrapped(["%d" % row[2] for row in species], "      "))
    add("   logical, parameter, public :: species_valence(n_species) = [ &")
    add(_wrapped(["%s" % (".true." if row[3] else ".false.") for row in species], "      "))
    add("")
    add("   ! X-ray IT92 (Cromer-Mann) coefficients: f(q) = sum_i a_i exp(-b_i (q/4pi)^2) + c")
    add("   real(rk), parameter, public :: species_xray_a(4, n_species) = reshape([ &")
    add(_wrapped([fortran_real(x) for row in species for x in row[4]], "      ", closer=" &"))
    add("     ], [4, n_species])")
    add("   real(rk), parameter, public :: species_xray_b(4, n_species) = reshape([ &")
    add(_wrapped([fortran_real(x) for row in species for x in row[5]], "      ", closer=" &"))
    add("     ], [4, n_species])")
    add("   real(rk), parameter, public :: species_xray_c(n_species) = [ &")
    add(_wrapped([fortran_real(row[6]) for row in species], "      "))
    add("")
    add("   ! Index of the neutral species of each element, 0 when the element")
    add("   ! has no X-ray row.")
    add("   integer, parameter, public :: element_neutral_species(n_elements) = [ &")
    add(_wrapped(["%d" % index for index in neutral_species], "      "))
    add("")
    add("   ! Source labels that repeat another row under a different name.")
    add("   character(len=%d), parameter, public :: species_alias_from(n_species_alias) = [ &" % SPECIES_LEN)
    add(_wrapped(['"%-*s"' % (SPECIES_LEN, row[0]) for row in aliases], "      "))
    add("   character(len=%d), parameter, public :: species_alias_to(n_species_alias) = [ &" % SPECIES_LEN)
    add(_wrapped(['"%-*s"' % (SPECIES_LEN, row[1]) for row in aliases], "      "))
    add("")
    add("end module sqc_element_data")
    add("")

    with open(out_path, "w") as handle:
        handle.write("\n".join(lines))

    missing_xray = [s for s in symbols if s not in it92]
    missing_neutron = [s for s in symbols if s not in nn92]
    n_ions = sum(1 for row in species if row[2] != 0)
    n_valence = sum(1 for row in species if row[3])
    sys.stderr.write(
        "%s: %d elements (%d without X-ray data, %d without neutron data)\n"
        % (out_path, len(symbols), len(missing_xray), len(missing_neutron))
    )
    sys.stderr.write(
        "%s: %d species (%d neutral, %d ions, %d valence), %d aliases: %s\n"
        % (
            out_path,
            len(species),
            len(species) - n_ions - n_valence,
            n_ions,
            n_valence,
            len(aliases),
            ", ".join("%s->%s" % pair for pair in aliases),
        )
    )


def _wrapped(items, indent, per_line=4, closer=" ]"):
    """Format a Fortran array constructor body, a few items per line.

    The last line is terminated with `closer`: either the closing bracket, or a
    continuation marker when the caller writes the bracket on the next line.
    """
    out = []
    for start in range(0, len(items), per_line):
        chunk = items[start : start + per_line]
        last = start + per_line >= len(items)
        sep = "" if last else ","
        out.append(indent + ", ".join(chunk) + sep + (closer if last else " &"))
    return "\n".join(out)


if __name__ == "__main__":
    main()
