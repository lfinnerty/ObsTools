# Installation and the workspace

## Install

With pip, in any Python ≥ 3.9 environment:

```bash
git clone https://github.com/lfinnerty/ObsTools.git
cd ObsTools
pip install -e .
```

With conda, first make an environment with the dependencies:

```bash
conda create -n hrccs python=3.11 numpy astropy astroquery matplotlib
conda activate hrccs
pip install -e .
```

To develop (tests and linting), install the dev extras and run the tests:

```bash
pip install -e ".[dev]"
```

```bash
pytest
```

The tests don't need network access. The window and night tests check the package against
reference outputs from the original ObsTools scripts.

Check the install:

```bash
hrccs-plan --version
```

```bash
hrccs-plan sites
```

## Network access

- **Planning** (`nights`, `windows`, `ephemeris`) runs offline. Site coordinates are built
  in, and astropy computes the Sun and Moon positions. Astropy may download Earth-orientation
  (IERS) tables the first time; it caches them.
- **Catalog commands** (`targetlist archive`, `catalog tic`, `starlist`) query the NASA
  Exoplanet Archive, MAST or SIMBAD. Each command sends one batched query. The result is cached
  in `<workspace>/cache/` and reused on later runs; pass `--refresh` to query again. Failed or
  empty queries are not cached.

## The workspace

Your target lists and all outputs live in a **workspace** directory, separate from the
package:

```
<workspace>/
    targetlists/<name>_targetlist.csv   your target lists
    starlists/                          telescope starlists
    output/best_windows/                ranked windows (CSV + text table)
    plots/                              nightly plots
    cache/                              cached catalog queries
```

The workspace is chosen in this order:

1. the `--workspace DIR` option;
2. the `HRCCS_WORKSPACE` environment variable;
3. the current directory.

Create one, with a copy of the example target list:

```bash
hrccs-plan init ~/hrccs_work
```

Then make it the default for your shell:

```bash
export HRCCS_WORKSPACE=~/hrccs_work
```

To keep it across sessions, add that `export` line to your `~/.bashrc`.

A target list is named by its file name without `_targetlist.csv`, so
`targetlists/my_planets_targetlist.csv` is `my_planets`. A path to any CSV file also works. The
bundled example list, `examples`, is found even if it's not in your workspace.

A workspace can be kept in its own (private) git repository, so that target lists and
observing plans are versioned separately from the code.
