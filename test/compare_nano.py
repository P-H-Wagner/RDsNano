#!/usr/bin/env python3
"""
Compare two nano files produced by run_nano_rds_cfg.py (e.g. before and after a code change).

Usage:  python3 compare_nano.py old.root new.root [--tree Events] [--skip ds_custom_m_err]

Reports, per branch: missing branches, different number of entries, and the
number of events where the values differ (NaN == NaN counts as equal).
Needs uproot and numpy (e.g. from an LCG view or `pip install --user uproot`).
"""
import argparse
import numpy as np
import uproot
import awkward as ak

parser = argparse.ArgumentParser()
parser.add_argument("old")
parser.add_argument("new")
parser.add_argument("--tree", default="Events")
parser.add_argument("--skip", nargs="*", default=[], help="branches expected to differ")
args = parser.parse_args()

t_old = uproot.open(args.old)[args.tree]
t_new = uproot.open(args.new)[args.tree]

b_old, b_new = set(t_old.keys()), set(t_new.keys())
for b in sorted(b_old - b_new): print("only in old:", b)
for b in sorted(b_new - b_old): print("only in new:", b)

if t_old.num_entries != t_new.num_entries:
    print("different number of events: %d vs %d" % (t_old.num_entries, t_new.num_entries))

n_bad = 0
for b in sorted(b_old & b_new):
    if b in args.skip:
        continue
    a = ak.flatten(t_old[b].array(), axis=None)
    c = ak.flatten(t_new[b].array(), axis=None)
    if len(a) != len(c):
        print("%-40s different length: %d vs %d" % (b, len(a), len(c)))
        n_bad += 1
        continue
    a, c = np.asarray(a), np.asarray(c)
    if a.dtype.kind == "f":
        diff = ~((a == c) | (np.isnan(a) & np.isnan(c)))
    else:
        diff = a != c
    if diff.any():
        print("%-40s %d / %d values differ" % (b, diff.sum(), len(a)))
        n_bad += 1

print("\n%d branches differ" % n_bad if n_bad else "\nall compared branches identical")
