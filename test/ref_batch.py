#!/usr/bin/env python3
# sqcalc - structure factors from LAMMPS dump trajectories
# Copyright (C) 2026 Rui Su, Hangzhou Dianzi University
#
# SPDX-License-Identifier: GPL-3.0-or-later

"""Run several numpy reference cases in a single process.

The reference scripts are independent implementations: each of them parses and
unwraps the trajectory itself.  Running them one per case therefore pays for
the same trajectory again and again, which dominates the cost of the small
cases.  This driver takes a JSON list of

  [{"case": "sqw", "args": ["--input", "traj.dump", ...]}, ...]

and calls the matching reference's ``main()`` in this process, so the frame
cache in ref_sqw hands the parsed trajectory to every case that asks for it.
The arguments are the ones the standalone scripts take.
"""

import argparse
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ref_fsqt  # noqa: E402  (sibling modules)
import ref_ngp  # noqa: E402
import ref_s4  # noqa: E402
import ref_sqw  # noqa: E402

CASES = {"sqw": ref_sqw, "s4": ref_s4, "fsqt": ref_fsqt, "ngp": ref_ngp}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--jobs", required=True,
                        help="JSON file with the list of reference cases")
    args = parser.parse_args()
    with open(args.jobs) as handle:
        jobs = json.load(handle)

    for job in jobs:
        module = CASES[job["case"]]
        argv = [module.__name__] + [str(value) for value in job["args"]]
        saved = sys.argv
        sys.argv = argv
        try:
            module.main()
        finally:
            sys.argv = saved


if __name__ == "__main__":
    main()
