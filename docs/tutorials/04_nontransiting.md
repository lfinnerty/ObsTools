# Tutorial 4: non-transiting planets

Non-transiting planets are detected in HRCCS by their orbital motion alone, so planning has
to deal with three unknowns:

- **the inclination**, which sets Kp = Kp_max sin i and how much the phase curve varies;
- **how much time you need**, since a bright target may be detected in a few hours;
- **the ephemeris uncertainty**, because the conjunction time can drift by hours by the time
  you observe.

The example is τ Boo b, which is non-transiting and has a measured inclination of about 45°.

```bash
hrccs-plan targetlist archive "tau Boo" -o tauboo --circular-below 0.05
```

## Inclination

`--inclination 45` applies to both the planet velocity (Kp = Kp_max sin 45°) and the
Lambertian phase curve. At low inclination the phase curve varies less, since the phase angle
only spans 90° − i to 90° + i, so f stays well above zero at every phase.

## Capped windows

```bash
hrccs-plan windows tauboo --site keck --start 2027-03-01 --end 2027-06-30 --inclination 45 \
    --max-hours 4 --airmass-k 0.05 --seeing-exponent 0.6 --sort snr --top 6
```

```
6 of 25 windows with relative S/N >= 0.3
rank Date       target         moon   illum   hrs UT start  UT end   ph0   ph1     v0     v1    dv  f_pl    X Xmax   S/N  ecl
   1 2027-03-21 tau Boo b      bright  1.00  4.02    10.13   14.13  0.48  0.53     17    -19    36  0.75 1.05 1.15  1.00
   2 2027-03-11 tau Boo b      dark    0.17  4.02    11.02   15.02  0.47  0.52     22    -14    36  0.75 1.05 1.17  1.00
   3 2027-03-31 tau Boo b      dark    0.28  4.02     9.27   13.27  0.48  0.53     11    -25    36  0.75 1.05 1.18  1.00
   4 2027-03-01 tau Boo b      dark    0.29  4.02    11.45   15.45  0.46  0.51     32     -4    36  0.75 1.05 1.15  0.99
   5 2027-04-10 tau Boo b      dark    0.23  4.02     8.38   12.38  0.49  0.54      5    -30    36  0.75 1.05 1.23  0.99
   6 2027-06-02 tau Boo b      dark    0.03  4.02     6.35   10.35  0.47  0.52     23    -12    36  0.75 1.08 1.34  0.99
```

- **How the cap works.** `--max-hours 4` keeps, on each night, the 4-hour stretch that
  maximizes the S/N score. That score comes from the brightness-weighted time
  (phase curve × airmass weight) and the velocity change Δv.
- **Where the windows land.** They centre on phase 0.5, where the dayside faces us. Here they
  sit near the zenith (X ≈ 1.05) because of the airmass weighting.
- **Ties.** Many nights score almost the same, so you can schedule around the Moon and other
  programs.
- **Relative S/N.** Values are normalized within one run, so S/N from capped and uncapped runs
  can't be compared directly.

## Ephemeris uncertainty

```bash
hrccs-plan ephemeris tauboo --date 2027-04-20T10:00 --t0-err 0.004 --period-err 0.0000047
```

```
target                n    T_conj (JD) sigma_t [h]  phase sigma_phase
tau Boo b          2539  2461515.15582        0.30  0.230       0.004
```

- **What it computes.** σ_t = √(σ_T0² + (n σ_P)²) after n orbits, and σ_phase = σ_t/P.
- **τ Boo b** has a well-determined RV ephemeris, so σ_phase is tiny.
- **Photometric candidates** are different. For example, planets from a TESS phase-curve
  search with periods known to about 10⁻⁴ d can reach σ_phase ≈ 0.1–0.2 after a few years.
  That is as large as a 4-hour window's phase coverage, so the window may miss the phases it
  was planned for.
- **Per-target errors.** If the target list has error columns (default names `T err (d)` and
  `P err (d)`), leave out `--t0-err` and `--period-err` and the per-target values are used.

Ways to reduce the risk:

- refine the ephemeris first;
- use longer windows;
- fit the planet's phase or conjunction time as a free parameter in the cross-correlation
  analysis.
