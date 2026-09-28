# Tutorial 3: ranking a semester of emission windows

`hrccs-plan windows` finds, for each target and night, the time when the planet's dayside is
visible at night. It then ranks those windows by an estimate of the planet S/N. The method is
in [docs/method.md](../method.md).

## Rank 2027A (May–July) at Gemini North

This uses the `my_planets` list from [tutorial 2](02_targetlists.md):

```bash
hrccs-plan windows my_planets --site gemini-n --start 2027-05-01 --end 2027-07-31 \
    --duration 2 --weight-kmag --airmass-k 0.05 --seeing-exponent 0.6 --top 8 --csv 2027A_my_planets.csv
```

```
8 of 70 windows with relative S/N >= 0.3
rank Date       target         moon   illum   hrs UT start  UT end   ph0   ph1     v0     v1    dv  f_pl    X Xmax   S/N  ecl
   1 2027-07-25 HD 189733 b    bright  0.53  6.08     6.38   14.45  0.42  0.57     72    -68   140  0.95 1.29 2.00  1.00    Y
   2 2027-07-14 HD 189733 b    bright  0.91  5.95     6.55   14.48  0.47  0.62     32   -102   133  0.92 1.21 2.00  0.95    Y
   3 2027-07-23 HD 189733 b    bright  0.72  8.17     6.40   14.55  0.52  0.67    -20   -136   116  0.82 1.23 1.97  0.94
   4 2027-07-16 HD 189733 b    bright  0.99  6.12     6.45   12.55  0.37  0.48    114     18    96  0.89 1.19 1.97  0.80
   5 2027-06-24 HD 189733 b    bright  0.67  4.48     7.87   14.33  0.48  0.60     22    -89   111  0.93 1.10 1.99  0.78    Y
   6 2027-06-13 HD 189733 b    bright  0.80  5.73     8.58   14.30  0.53  0.64    -31   -118    87  0.86 1.20 2.00  0.70
   7 2027-06-15 HD 189733 b    bright  0.93  3.88     8.45   14.32  0.43  0.54     64    -39   103  0.97 1.29 2.00  0.70    Y
   8 2027-07-12 HD 189733 b    bright  0.76  7.80     6.68   14.47  0.57  0.71    -63   -149    86  0.69 1.21 2.00  0.66
Saved /home/you/hrccs_work/output/best_windows/2027A_my_planets.csv and /home/you/hrccs_work/output/best_windows/2027A_my_planets.txt
```

What the options do:

- **`--duration 2`** excludes a 2-hour secondary eclipse around phase 0.5. For transiting
  planets, use the real eclipse duration.
- **`--weight-kmag`** weights each target by its K-band photon flux. HD 189733 (K = 5.5) then
  outranks KELT-20 (K = 7.4), because the S/N is compared *between* targets.
- **`--airmass-k 0.05 --seeing-exponent 0.6`** down-weights time at high airmass, through
  extinction and slit losses.
- **The columns** are: the usable hours and UT range; the orbital phase (`ph0`, `ph1`) and
  planet velocity (`v0`, `v1`, `dv`) at the start and end; the mean phase-curve factor
  (`f_pl`); the mean and maximum airmass (`X`, `Xmax`); the relative S/N; and whether the
  window crosses the eclipse (`ecl`).
- **Sorting** puts bright time first, since IR HRCCS is insensitive to moonlight, then orders
  by S/N. Use `--sort snr` to rank by S/N alone.
- **Missing targets.** WASP-121 (Dec −39°) does not appear, because it is not observable at
  night from Maunakea in May–July.

## One schedule per target

To choose a night for *each* target, normalize the S/N per target and write one file per
target:

```bash
hrccs-plan windows my_planets --site gemini-n --start 2027-05-01 --end 2027-07-31 --duration 2 \
    --airmass-k 0.05 --seeing-exponent 0.6 --normalize-per-target --snr-min 0.9 --sort snr \
    --csv 2027A_per_target.csv --split-targets
```

```
7 of 70 windows with relative S/N >= 0.9
rank Date       target         moon   illum   hrs UT start  UT end   ph0   ph1     v0     v1    dv  f_pl    X Xmax   S/N  ecl
   1 2027-07-25 HD 189733 b    bright  0.53  6.08     6.38   14.45  0.42  0.57     72    -68   140  0.95 1.29 2.00  1.00    Y
   2 2027-07-27 KELT-20 b      dark    0.32  7.58     6.47   14.03  0.51  0.60    -13   -102    89  0.93 1.22 1.99  1.00
   3 2027-07-13 KELT-20 b      bright  0.84  6.02     6.47   14.47  0.48  0.58     19    -80    99  0.96 1.20 1.79  0.97    Y
   4 2027-07-14 HD 189733 b    bright  0.91  5.95     6.55   14.48  0.47  0.62     32   -102   133  0.92 1.21 2.00  0.95    Y
   5 2027-07-06 KELT-20 b      dark    0.16  5.82     6.62   14.42  0.47  0.56     33    -65    98  0.98 1.25 2.00  0.95    Y
   6 2027-07-23 HD 189733 b    bright  0.72  8.17     6.40   14.55  0.52  0.67    -20   -136   116  0.82 1.23 1.97  0.94
   7 2027-07-20 KELT-20 b      bright  0.93  6.80     7.72   14.50  0.51  0.59    -13    -94    81  0.94 1.21 2.00  0.92
Saved .../output/best_windows/2027A_per_target.csv and .../2027A_per_target.txt
Saved .../output/best_windows/2027A_per_target_HD189733b.csv and .../2027A_per_target_HD189733b.txt
Saved .../output/best_windows/2027A_per_target_KELT20b.csv and .../2027A_per_target_KELT20b.txt
```

Each `.txt` file starts with the command that made it, so every result can be reproduced.

## Weighting by planet contrast

If your target list has a column with the expected planet/star contrast, the ranking can also
weight by it: `--weight-contrast --contrast-column "K contrast (ppm)"`. That column name is the
default.

`hrccs_planner.contrast.dayside_contrast` computes a blackbody estimate from the Cowan & Agol
(2011) dayside temperature, as in [tutorial 6](06_python_library.md).

## Other useful options

| Option | Effect |
|---|---|
| `--alt-min 35` | stricter altitude limit (default 30°, airmass 2) |
| `--site keck2` | Maunakea with the Keck II Nasmyth-deck pointing limit |
| `--twilight -12` | allow nautical twilight |
| `--moon-sep-min 20` | larger Moon avoidance |
| `--delta-v-min 50` | require more planet-velocity change |
| `--max-hours 4` | cap each window at 4 hours ([tutorial 4](04_nontransiting.md)) |
| `--top N` | keep the best N windows |
