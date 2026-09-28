# Tutorial 6: using the Python library

Everything the CLI does is available from Python.

```python
from hrccs_planner import get_site, read_targets
from hrccs_planner.visibility import dates_between
from hrccs_planner.windows import WindowOptions, rank_windows
from hrccs_planner.workspace import resolve_targetlist

targets = read_targets(resolve_targetlist('my_planets'))   # or read_targets('path/to/list.csv')
kelt20 = [t for t in targets if t.name == 'KELT-20 b']
opts = WindowOptions(duration=2.0, max_hours=4, airmass_k=0.05, seeing_exponent=0.6, sort='snr', top=3)
rows, n_all = rank_windows(kelt20, get_site('gemini-n'), dates_between('2027-06-01', '2027-07-31'), opts)
for r in rows:
    print(f"{r['date']}  UT {r['ut0']:5.2f}-{r['ut1']:5.2f}  phase {r['ph0']:.2f}-{r['ph1']:.2f}  "
          f"dv {r['dv']:3.0f} km/s  S/N {r['snr']:.2f}")
```

```
2027-07-13  UT  8.95-12.95  phase 0.51-0.56  dv  50 km/s  S/N 1.00
2027-07-20  UT  7.72-11.72  phase 0.51-0.56  dv  50 km/s  S/N 1.00
2027-06-01  UT 10.32-14.32  phase 0.44-0.49  dv  50 km/s  S/N 0.99
```

With the 2-hour eclipse excluded, each 4-hour window sits on one side of secondary eclipse.

## Modules

| Module | Contents |
|---|---|
| `targets` | `Target`, `read_targets`, `write_targets` (columns matched by name; extra columns kept in `Target.extra`) |
| `sites` | `get_site`, `SITES`, the built-in site registry and night time grids, `keck2_min_altitude` |
| `ephemeris` | `orbital_phase`, `planet_rv_circular`, `planet_rv_eccentric`, `solve_kepler`, `lambert_phase`, `propagate_conjunction` |
| `visibility` | `night_sky` (Sun and Moon over a night), `airmass`, `moon_too_close`, `dates_between` |
| `windows` | `WindowOptions`, `rank_windows`, `best_subwindow`, `write_windows`, `read_windows` |
| `nights` | `NightOptions`, `observable_targets`, `plot_night`, `plan_nights` |
| `contrast` | `dayside_temperature`, `dayside_contrast`, `band_planck`, `semi_major_axis`, `kp_max` |
| `catalogs` | `query_exoarchive`, `exoarchive_target`, `query_tic`, `query_simbad` (all cached) |
| `starlists` | `keck_lines`, `magellan_lines`, `write_starlist` |

## Example: add an expected-contrast column

```python
from hrccs_planner import read_targets, write_targets
from hrccs_planner.contrast import dayside_contrast, semi_major_axis, RSUN

targets = read_targets('my_list.csv')
for t in targets:
    teff, rstar, mstar = 6500., 1.5, 1.4        # from your own catalog
    a_rs = semi_major_axis(t.period, mstar) / (rstar * RSUN)
    c_k, tday = dayside_contrast(teff, rstar, a_rs, rp_rjup=1.5, band='K')
    t.extra['K contrast (ppm)'] = f'{c_k * 1e6:.0f}'
write_targets('my_list_with_contrast.csv', targets)
```

Then rank the windows with `hrccs-plan windows my_list_with_contrast.csv --weight-contrast ...`.
