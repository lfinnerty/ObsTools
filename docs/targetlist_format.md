# Target-list format

A target list is a CSV file with a header row and one row per **planet**. Columns are matched
**by name**, case-insensitively and in any order. Extra columns are kept: they are carried
through to the outputs and can be used for weighting, e.g. `--contrast-column`.

## Standard columns

| Column | Required | Units / values | Notes |
|---|---|---|---|
| `Name` | yes | text | Planet name, e.g. `WASP-121 b`. The host name (for starlists) is this without the trailing planet letter. |
| `RA` | yes | `hh:mm:ss.s`, or degrees | J2000 |
| `Dec` | yes | `±dd:mm:ss`, or degrees | J2000 |
| `Period (day)` | yes | days | |
| `a (AU)` | | AU | Informational |
| `T transit/inf conj (JD)` | yes | JD or BJD | Transit, or inferior conjunction for non-transiting planets. Phase 0. |
| `Obs status` | | text | Informational; `KPIC obs status` is also accepted |
| `Eccentric?` | | Y/Yes/N/No | Defaults to `e > 0` if missing. Eccentric orbits use `e` and `omega` in `nights` and `windows`. |
| `Kp max (km/s)` | yes | km/s | **Total** orbital velocity 2πa/P, i.e. Kp for sin i = 1 (see below) |
| `Transiting?` | | Y/N | Informational |
| `Kmag (<9)` | | mag | Used by `windows --weight-kmag` |
| `e` | | | Eccentricity |
| `omega` | | deg | The star's argument of periastron (the RV convention, as in the NASA Exoplanet Archive), for eccentric orbits |
| `inc` | | deg | Informational. The planners use `--inclination` for all targets. |

Minimal example:

```csv
Name,RA,Dec,Period (day),T transit/inf conj (JD),Kp max (km/s),Kmag (<9)
WASP-33 b,02:26:51.06,+37:33:01.6,1.21987,2459883.198788,213,7.468
```

`hrccs-plan targetlist archive` writes the full standard header. Legacy lists, where the
first 14 columns are in the order above, are read unchanged.

## Kp conventions

- **`Kp max` is the total orbital velocity** 2πa/P. The planners multiply it by sin i, from
  `--inclination`, which defaults to 90°.
- For transiting planets, leave `--inclination` at 90: sin i ≈ 1.
- For non-transiting planets, the inclination is unknown. Choose a value to plan with (e.g.
  `--inclination 30`), and see [tutorial 4](tutorials/04_nontransiting.md).
- If a list already contains the observed Kp sin i instead, as KPIC-Synthgen's
  `make_downstream_runfiles.py` writes, keep `--inclination 90`.

## Velocity and phase convention

- **Phase 0** is `T transit/inf conj` and **phase 0.5** is secondary eclipse (for circular
  orbits). For eccentric orbits `T transit/inf conj` is still the transit (inferior
  conjunction), and the geometry uses the orbital phase (0.5 at secondary eclipse); see
  docs/method.md.
- The planet's velocity relative to the star is v = Kp sin(2π phase). It is positive, i.e.
  receding, just after transit, and decreases through secondary eclipse.
- For eccentric orbits, the Keplerian curve is used: v = −Kp (cos(f + ω) + e cos ω)/√(1 − e²), with Kp = `Kp max` sin i, so `Kp max` stays 2πa/P.
