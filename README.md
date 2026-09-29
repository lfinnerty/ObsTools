# hrccs_planner

Observation planning for **high-resolution cross-correlation spectroscopy (HRCCS)** of
exoplanet atmospheres: when is a planet observable, how much of its orbit (and velocity
change) can you cover in a night, and which nights give the best signal for dayside
emission. Built for KPIC, IGRINS-2, WINERED and similar spectrographs, and usable at any
observatory.

- **`hrccs-plan nights`**: a plot per night of the airmass and the planet's line-of-sight
  velocity for every target observable on its dayside.
- **`hrccs-plan windows`**: ranks nights, and the best window in each night, by the expected
  planet S/N. It accounts for:
  - the phase curve (with inclination) and the planet-velocity change;
  - the stellar K magnitude and the planet contrast;
  - airmass (extinction and slit losses);
  - moon separation, twilight and secondary eclipse.

  Windows can be capped at a maximum length.
- **`hrccs-plan transits`**: ranks transit windows for transmission spectroscopy: the transit
  plus an out-of-transit baseline, split around the transit or on either side, all above
  airmass 2 and weighted by the in-transit airmass.
- **Gemini timing windows:** `windows` and `transits` can both write their results
  (`--gemini-tw`) in the format of the Gemini PIT proposal's Scheduling field.
- **`hrccs-plan ephemeris`**: the orbital phase and its uncertainty on a date, which matters
  for non-transiting planets.
- **`hrccs-plan phase-sensitivity`**: how much S/N the best dayside windows lose to that
  uncertainty, assuming the analysis refits the conjunction. Compare scenarios, e.g. the
  current ephemeris and one improved by new RVs.
- **`hrccs-plan targetlist` / `catalog` / `starlist`**: target lists from the NASA Exoplanet
  Archive, TIC look-ups, and Keck and Magellan starlists with SIMBAD proper motions. Catalog
  queries are cached.

![Example nightly plan](docs/images/nights_example.png)

## Install

```bash
git clone https://github.com/lfinnerty/ObsTools.git
cd ObsTools
pip install -e .
```

This needs Python ≥ 3.9, numpy, astropy, astroquery and matplotlib. See
[docs/install.md](docs/install.md) for conda and development installs.

## Quick start

```bash
hrccs-plan init ~/hrccs_work               # workspace with an example target list
export HRCCS_WORKSPACE=~/hrccs_work
hrccs-plan nights examples --site keck --start 2027-07-10 --end 2027-07-16
hrccs-plan windows examples --site gemini-n --start 2027-05-01 --end 2027-07-31 \
    --weight-kmag --airmass-k 0.05 --seeing-exponent 0.6 --csv 2027A.csv
```

Plots go to `~/hrccs_work/plots/`, and ranked windows to `~/hrccs_work/output/best_windows/`.
Run `hrccs-plan COMMAND --help` for every option.

## Documentation

- [Installation and the workspace](docs/install.md)
- [Target-list format](docs/targetlist_format.md)
- [Method: usable time, S/N score, weights](docs/method.md)
- Tutorials:
  1. [Your first night plan](docs/tutorials/01_first_night_plan.md)
  2. [Target lists from the NASA Exoplanet Archive](docs/tutorials/02_targetlists.md)
  3. [Ranking a semester of emission windows](docs/tutorials/03_ranking_windows.md)
  4. [Non-transiting planets: inclination, capped windows, ephemeris errors](docs/tutorials/04_nontransiting.md)
  5. [From windows to the telescope: starlists and simulations](docs/tutorials/05_starlists_and_simulations.md)
  6. [Using the Python library](docs/tutorials/06_python_library.md)
  7. [Transit windows](docs/tutorials/07_transit_windows.md)

## Roadmap

- Emission-window ranking for eccentric orbits.
- Instrument presets: slit width, overheads, saturation-aware exposure times.
- A PyPI release.

## History

This repository was ObsTools, a set of observing scripts. `best_window.py` and
`observing_planner.py` still work as thin wrappers around `hrccs-plan windows` and
`hrccs-plan nights`. See [CHANGELOG.md](CHANGELOG.md).
