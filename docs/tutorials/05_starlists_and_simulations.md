# Tutorial 5: from windows to the telescope

## Starlists

`hrccs-plan starlist` looks up the host stars in SIMBAD with one batched query, which is
cached, and writes a starlist with J2000 positions and proper motions.

The Keck format, for the hosts in a target list:

```bash
hrccs-plan starlist keck --list my_planets -o my_planets_keck.txt
```

```
Querying simbad ...
HD 189733       20 00 43.713    +22 42 39.07    2000 pmra=-0.0002    pmdec=-0.2503   Rmag=7.13
WASP-121        07 10 24.060    -39 05 50.57    2000 pmra=-0.0003    pmdec=0.0257    Vmag=10.51
KELT-20         19 38 38.735    +31 13 09.22    2000 pmra=0.0002     pmdec=-0.0063   Vmag=7.58
Saved /home/you/hrccs_work/starlists/my_planets_keck.txt
```

The Magellan (TCS) format, for star names given directly:

```bash
hrccs-plan starlist magellan --names WASP-121 "HD 189733" -o my_planets_magellan.txt
```

```
# RA Dec equinox RApm Decpm offset rot RA_probe1 Dec_probe1 equinox RA_probe2 Dec_probe2 equinox pm_epoch
###name hh:mm:ss.s sdd:mm:ss yyyy.0 s.ss s.ss angle mode hh:mm:ss.s sdd:mm:ss yyyy.0 hh:mm:ss.s sdd:mm:ss yyyy.0 yyyy.0
001 WASP-121 07:10:24.06 -39:05:50.57 2000.0 -0.0 0.026 0.0 GRV 0 0 2000.0 0 0 2000.0 2017.5 # V=10.5
002 HD_189733 20:00:43.71 +22:42:39.07 2000.0 -0.0 -0.25 0.0 GRV 0 0 2000.0 0 0 2000.0 2017.5 # V=7.6
```

- **Proper motion in RA** is written as the rate of change of the RA *coordinate* in seconds
  of time per year, μα*/(15000 cos δ). The proper motion in Dec is in arcsec/yr.
- **Check your observatory's format.** Starlist formats vary by observatory and instrument,
  so check the current documentation, and adjust the settings or edit the file if needed.

## Simulating the best windows

The ranked-window CSV (`output/best_windows/*.csv`) is meant to be read by other tools.

- **Columns:** `hst_date` (the local date of the night at any site), `target`, `ut_start`,
  `ut_end` (hours relative to date 24:00 UTC), phases, velocities, `rel_snr` and the airmass
  columns.
- **KPIC-Synthgen** (`make_blind_uhj_runfiles.py --windows ... --targetlist ...`) turns the
  best window per target into IGRINS runfiles, to simulate the observation end to end.
- **KPIC-Synthgen's `make_downstream_runfiles.py`** writes a target list for a simulated
  planet into the workspace. It also prints the `hrccs-plan windows` command to rank its
  nights.
