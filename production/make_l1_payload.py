#!/usr/bin/env python3
"""Build data/l1_seeds_bparking2018.json for the L1 seed table (plugins/L1SeedsTableProducer.cc) from the
inputs of the Bmmm ntuplizer, read with Bmmm's own loader (Bmmm/Analysis/python/L1Seeds2018.py), so that
both ntuples carry the same L1 decisions and prescales:

  prescales : {seed: {run: [[ls_first, ls_last, prescale], ...]}}   brilcalc tables merged over the ten
              parking HLT paths (union of the lumi sections, as Bmmm), run-length encoded
  run_menu  : {run: menu name} for 2018 (l1menus/goodRuns2013to2022ByYear.json)
  menus     : {menu name: {seed: algorithm bit}} for the 2018 menus, restricted to the seeds
  mc_menu   : the menu used for MC (L1Menu_Collisions2018_v1_0_0, as Bmmm)

Usage: make_l1_payload.py <Bmmm/Analysis/data> <out.json>"""
import sys, os, json, pickle, importlib.util

datadir, out = sys.argv[1:3]
spec = importlib.util.spec_from_file_location(
    "L1Seeds2018", os.path.join(datadir, "..", "python", "L1Seeds2018.py"))
l1 = importlib.util.module_from_spec(spec); spec.loader.exec_module(l1)

PATHS = ["HLT_Mu7_IP4", "HLT_Mu8_IP3", "HLT_Mu8_IP5", "HLT_Mu8_IP6", "HLT_Mu8p5_IP3p5", "HLT_Mu9_IP4",
         "HLT_Mu9_IP5", "HLT_Mu9_IP6", "HLT_Mu10p5_IP3p5", "HLT_Mu12_IP6"]
tables, src = l1.load_tables(PATHS, datadir)

prescales = {}
for seed, runs in sorted(tables.items()):
    dst = prescales.setdefault(seed, {})
    for run, lss in sorted(runs.items()):
        segs = []
        for ls in sorted(lss):
            ps = float(lss[ls])
            if segs and segs[-1][1] == ls - 1 and segs[-1][2] == ps:
                segs[-1][1] = ls
            else:
                segs.append([int(ls), int(ls), ps])
        dst[str(run)] = segs

with open(os.path.join(datadir, "l1menus", "goodRuns2013to2022ByYear.json")) as fin:
    good = json.load(fin)
run_menu = {str(r["run_number"]): r["l1_menu"] for r in good["2018"]}

menus = {}
for name in l1.MENUS_2018:
    with open(os.path.join(datadir, "l1menus", "%s.pickle" % name), "rb") as fin:
        m = pickle.load(fin)
    menus[name] = {s: int(m[s]) for s in prescales if s in m}

payload = {"description": "L1 seeds of the 2018 B-parking paths: brilcalc prescales (merged over the ten "
                          "HLT paths, run-length encoded [ls_first, ls_last, ps]), run -> L1 menu, seed -> bit. "
                          "Built by production/make_l1_payload.py from the Bmmm inputs (%s)." % src,
           "seeds": sorted(prescales), "mc_menu": l1.MC_MENU,
           "prescales": prescales, "run_menu": run_menu, "menus": menus}
with open(out, "w") as fout:
    json.dump(payload, fout, separators=(",", ":"))
nseg = sum(len(v) for r in prescales.values() for v in r.values())
print("seeds", len(prescales), "runs", len(next(iter(prescales.values()))), "segments", nseg,
      "menus", {k: len(v) for k, v in menus.items()}, "bytes", os.path.getsize(out))
