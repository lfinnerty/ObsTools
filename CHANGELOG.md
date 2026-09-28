# Changelog

## Unreleased

- `windows --gemini-tw FILE`: write the ranked windows as Gemini PIT timing windows
  (`YYYY-MM-DD HH:MM:SS H:MM`, UT start and duration, one block per target) for the proposal's
  Scheduling field; `--tw-round` sets the rounding (default 15 minutes).

## 0.1.0 — ObsTools becomes `hrccs_planner`

The scripts are now an installable package (`pip install -e .`) with one command-line
tool, `hrccs-plan`, plus documentation, tutorials and tests.

### Where things went

| Before | Now |
|---|---|
| `observing_planner.py SITE LIST ...` | `hrccs-plan nights LIST --site SITE ...` (the old script is a wrapper) |
| `best_window.py SITE LIST ...` | `hrccs-plan windows LIST --site SITE ...` (the old script is a wrapper) |
| `planet_props_exoarchive.py NAMES -o NAME` | `hrccs-plan targetlist archive NAMES -o NAME` |
| `get_tic_data.py` (hard-coded IDs) | `hrccs-plan catalog tic ID ...` |
| `targetlist_maker_keck.py` (hard-coded names) | `hrccs-plan starlist keck --list LIST` or `--names ...` |
| `targetlist_maker_magellan.py` (hard-coded names) | `hrccs-plan starlist magellan --list LIST` or `--names ...` |
| `targetlist_maker_tess_blind.py` | moved to the TESS blind-search code (`tess_targetlist.py`) |
| `noise_estimation/` | removed (target-specific scripts); blackbody contrast is in `hrccs_planner.contrast` |
| `targetlists/`, `starlists/`, `other_data/`, `plots/`, `output/` in the repo | a separate workspace (`--workspace` / `$HRCCS_WORKSPACE`); see docs/install.md |

### Unchanged results

`hrccs-plan windows` reproduces the ranked-window CSVs of `best_window.py` exactly. Likewise,
`hrccs-plan nights` selects the same targets on the same nights, with the same curves, as
`observing_planner.py`. Both are checked by the tests.

### Changes in behaviour

- **Target lists** are read by column name instead of position. Legacy headers
  (`KPIC obs status`) and quoted names containing commas now work. Missing or invalid values
  give an error naming the file and line.
- **`nights`:**
  - With `--require-conjunction`, the airmass panel now shows only the targets drawn in the
    velocity panel. Before, it showed every target observable on its dayside.
  - The Moon is drawn as a dotted black line, so it isn't confused with a target.
  - The title names the site and gives the night's local date.
- **`windows`:**
  - The text table's date column is labelled `Date` (the local date of the night), not
    `HST date`. The CSV column is still `hst_date`, for compatibility.
  - Eccentric targets are reported once per run, not once per night.
- **Starlists:**
  - The Keck `pmra` is now μα*/(15000 cos δ), in seconds of time per year. SIMBAD's pmra
    already includes cos δ; the old script left out the 1/cos δ.
  - Magellan names have spaces replaced by `_`, which the space-delimited format needs, and
    the proper-motion conversion uses the exact declination.
  - Coordinates are rounded rather than truncated.
- **New options:**
  - `targetlist archive --circular-below E`: treat e < E as circular. The Archive's small
    eccentricities otherwise make `windows` skip planets such as WASP-121 b.
  - `windows`: `--twilight`, `--sort snr`, `--site keck2` (Keck II pointing limits).
  - `hrccs-plan ephemeris`: phase and phase uncertainty on a date.
  - `hrccs-plan init`: create a workspace.
- **Catalog queries** are batched and cached in `<workspace>/cache/`.
