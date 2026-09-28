# Tutorial 2: target lists from the NASA Exoplanet Archive

`hrccs-plan targetlist archive` builds a target list from the archive's *Planetary Systems
Composite Parameters* table (`pscomppars`). It makes one batched query for all the names you
give, which can be host names or planet names.

```bash
hrccs-plan targetlist archive WASP-121 KELT-20 "HD 189733" -o my_planets --circular-below 0.05
```

```
Querying exoarchive ...
Wrote 3 planet(s) to /home/you/hrccs_work/targetlists/my_planets_targetlist.csv
```

```
Name,RA,Dec,Period (day),a (AU),T transit/inf conj (JD),Obs status,Eccentric?,Kp max (km/s),Transiting?,Kmag (<9),e,omega,inc
HD 189733 b,20:00:43.71,+22:42:35.19,2.21857567,0.03126,2453955.5255511,None,No,153,Y,5.541,0.0,20.0,85.71
WASP-121 b,07:10:24.06,-39:05:50.17,1.27492504,0.02571,2460245.02038,None,No,219,Y,9.374,0.0,88.95,88.49
KELT-20 b,19:38:38.74,+31:13:09.12,3.4741085,0.0542,2457503.120049,None,No,170,Y,7.415,0.0,,86.12
```

What each column comes from:

- **`Kp max`** is 2πa/P, the total orbital velocity.
- **`T transit/inf conj`** is the archive's transit midpoint.
- **`Kmag`** is the system K magnitude.

## Near-circular orbits

The archive gives WASP-121 b e = 0.0085. Any e > 0 marks a planet as eccentric, and
`hrccs-plan windows` does not yet rank eccentric orbits, so it would skip WASP-121 b.
`--circular-below 0.05` writes orbits with e < 0.05 as circular.

You can also edit the list by hand: set `Eccentric?` to `No`.

## Caching

Running the same command again doesn't query the archive a second time. The result is read
from `<workspace>/cache/`:

```bash
hrccs-plan targetlist archive WASP-121 KELT-20 "HD 189733" -o my_planets --circular-below 0.05
```

```
Wrote 3 planet(s) to /home/you/hrccs_work/targetlists/my_planets_targetlist.csv
```

Use `--refresh` to query again, for example after the archive has been updated.

## Other sources

A target list is a plain CSV, so you can also:

- **write or edit one by hand.** Only `Name`, `RA`, `Dec`, `Period (day)`,
  `T transit/inf conj (JD)` and `Kp max (km/s)` are required; see
  [the format reference](../targetlist_format.md).
- **add your own columns,** such as notes, priorities or an expected contrast. They are kept,
  and a contrast column can weight the ranking (`windows --weight-contrast --contrast-column NAME`).
- **generate one from a survey.** For example, `tess_targetlist.py` in the TESS blind-search
  repository writes lists with ephemeris-error and contrast columns.

`hrccs-plan catalog tic ID ...` looks up TIC coordinates, magnitudes and stellar parameters
for TESS targets that aren't in the archive.
