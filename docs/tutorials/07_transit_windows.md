# Tutorial 7: transit windows

For transmission spectroscopy you need the whole transit, plus an out-of-transit baseline to
normalize against. `hrccs-plan transits` finds every transit that fits in a night with its
baseline, with airmass < 2 and darkness throughout. It then ranks the windows by the in-transit
airmass.

The example is TOI-1408 b in 2027A at Gemini North. Its transit is grazing, with T14 = 1.46 h,
and 2 h of baseline is requested.

## Baseline split around the transit

```bash
hrccs-plan transits TOI-1408b --site gemini-n --start 2027-02-01 --end 2027-07-31 --duration 1.46 \
    --baseline 2 --baseline-mode split --airmass-k 0.05 --seeing-exponent 0.6 --sort snr \
    --csv transits_TOI-1408b_2027A_split.csv --gemini-tw TOI-1408b_2027A_transit_tw_split.txt
```

```
1 of 1 transit windows with relative S/N >= 0.3 (baseline 2 h, split)
rank Date       target         moon   illum UT mid   T14   pre  post UT start  UT end   hrs  X_tr  Xmax   S/N
   1 2027-07-29 TOI-1408 b     dark    0.13  10.03  1.46  1.00  1.00     8.30   11.76  3.46  1.67  1.82  1.00

TW for TOI-1408 b
--------------------------------------------
2027-07-30 08:15:00 3:30
```

- **The columns:**
  - `UT mid` is mid-transit, in UT hours relative to 24:00 UT on the local date.
  - `pre` and `post` are the baselines before and after the transit.
  - `X_tr` is the mean airmass during the transit, and `Xmax` the highest airmass in the whole
    window.
- **The timing window** is the whole observation, transit plus baseline, rounded out to
  15 minutes. Its date is in UT, the day after the HST night.

## Baseline on either side

With `--baseline-mode any`, the baseline can go on either side of the transit, so more
transits qualify:

```bash
hrccs-plan transits TOI-1408b --site gemini-n --start 2027-02-01 --end 2027-07-31 --duration 1.46 \
    --baseline 2 --baseline-mode any --airmass-k 0.05 --seeing-exponent 0.6 --sort snr \
    --csv transits_TOI-1408b_2027A_any.csv --gemini-tw TOI-1408b_2027A_transit_tw_any.txt
```

```
3 of 3 transit windows with relative S/N >= 0.3 (baseline 2 h, any)
rank Date       target         moon   illum UT mid   T14   pre  post UT start  UT end   hrs  X_tr  Xmax   S/N
   1 2027-07-29 TOI-1408 b     dark    0.13  10.03  1.46  1.00  1.00     8.30   11.76  3.46  1.67  1.82  1.00
   2 2027-06-28 TOI-1408 b     dark    0.28  10.68  1.46  0.57  1.43     9.38   12.84  3.46  1.78  2.00  0.98
   3 2027-07-20 TOI-1408 b     bright  0.93  13.64  1.46  1.83  0.17    11.08   14.54  3.46  1.81  1.96  0.97

TW for TOI-1408 b
--------------------------------------------
2027-06-29 09:15:00 3:45
2027-07-21 11:00:00 3:45
2027-07-30 08:15:00 3:30
```

On 06-28 the target only rises above airmass 2 shortly before ingress, so most of the baseline
is taken after egress. On 07-20, morning twilight follows the transit, so most of the baseline
comes first. When an even split fits, as on 07-29, it is used.

## Notes

- **Per-target durations:** `hrccs-plan targetlist archive` writes the archive's T14 as a
  `T14 (h)` column; pass `--duration-column "T14 (h)"`. You can also add the column by hand.
- **Comparing targets:** `--weight-kmag` (stellar brightness) and `--weight-contrast` with a
  column such as an expected transmission signal in ppm rank several targets against each other.
- **Ephemeris errors:** the windows use T0 and P only. Check `hrccs-plan ephemeris` if the
  ephemeris is old. If σ_t is comparable to the baseline, request a longer baseline to keep
  margin.
- **Non-transiting targets** (`Transiting? = N`) are skipped.
