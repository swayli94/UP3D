# ★★★ THIS DATASET IS WITHDRAWN — it is a FIRST-ORDER solution

Every CSV and figure in this directory was produced with
`NITFO = 1000` and `NCYC = 1000`, which makes CFL3D run **first order in
space on the finest grid for the entire solve**, even though the deck
requests `RKAP0 = 0.3333` (κ = 1/3 MUSCL).

`cfl3d/libs/resid.F:141`

```fortran
if (icyc.ge.nitfo+1 .and. level.ge.lglobal) then  rkap = rkap0
else if (level.ge.lglobal) then                   rkap = -3.   ! first order
```

and `cfl3d/dist/mgblk.F:416` sets `nitfo = nitfo1(iseq)` immediately before
`do 7000 icyc=1,ncyc` — so `icyc` **restarts at 1 on every mesh-sequence
level** and is compared against that level's own `nitfo`. With
`nitfo = ncyc`, `icyc >= nitfo+1` is never satisfied.

## What this means for the numbers here

- Everything — cl, cd, cm, the shock positions, the Cp curves — is a
  first-order-accurate solution. The grid-convergence table's implied orders
  (cd `p = 0.68`, shocks 4.90 → 0.43) are the orders **of a first-order
  scheme**, which is why they look the way they do.
- ★★ The experiment bias is **confounded**: `experiment_bias.csv` reports the
  shock 0.025–0.072 c aft of the measurement and calls the direction
  "as predicted", on the viscous-displacement argument. First-order numerical
  dissipation displaces a captured shock the **same way**, so two mechanisms
  push in one direction and the single-sided direction test cannot separate
  them. The reading is not wrong, it is *unattributed*.
- ★ The suction-peak correction in `compare_m6_experiment.py` still stands —
  `cp_min` really is non-asymptotic on every station — but its stated reason
  ("a finite grid under-resolves a leading-edge peak") is incomplete: first-order
  dissipation is a stronger source pointing the same way.
- The **withdrawn y/b = 0.99 shock position** and the `detector_premise`
  machinery are unaffected by this: that defect is in the detector's premise,
  not in the flow solution, and the same check carries forward.

## What is NOT affected

The four **2-D** datasets (`euler_naca0012`, `euler_rae2822`,
`rans_naca0012`, `rans_rae2822`, i.e. D05/D06/D08/D09). Their recipes carry
`NITFO` 0 / 500 / 1500 against `NCYC` 1000 / 2000 / 2000 — every one switches
to the requested order. Verified by reading the generated decks, and now
asserted at import time in `cfl3d_runner.py`.

## Replacement

A four-rung second-order ladder (L1–L4, the new `L4` being 11.278 M points)
is being run. `wing3d_otip.py::write_inp` now **raises** on
`nitfo >= ncyc`, so this cannot recur silently.

Record: `docs/dev_phase_six/20260903-0200-cfl3d-first-order-defect.md`.
