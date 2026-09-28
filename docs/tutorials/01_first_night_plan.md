# Tutorial 1: your first night plan

In this tutorial you'll plot, for a week of nights at Maunakea, which planets in a target list
can be observed on their dayside, and when.

## 1. Create a workspace

```bash
hrccs-plan init ~/hrccs_work
export HRCCS_WORKSPACE=~/hrccs_work
```

```
Workspace: /home/you/hrccs_work
  copied example targetlists/examples_targetlist.csv
Use it with --workspace /home/you/hrccs_work, or: export HRCCS_WORKSPACE=/home/you/hrccs_work
```

`examples` is a list of nine well-known hot and ultra-hot Jupiters from the NASA Exoplanet
Archive: WASP-121 b, WASP-33 b, KELT-9 b, KELT-20 b, WASP-76 b, WASP-189 b, WASP-18 b,
HD 80606 b and τ Boo b. Open `~/hrccs_work/targetlists/examples_targetlist.csv` to see the
format, which is described in [the format reference](../targetlist_format.md).

## 2. Plot a week of nights

```bash
hrccs-plan nights examples --site keck --start 2027-07-10 --end 2027-07-16
```

```
Making plot for 2027-07-10 at MKO: WASP-76 b, KELT-20 b, WASP-189 b, WASP-33 b
Making plot for 2027-07-11 at MKO: KELT-9 b, WASP-33 b
Making plot for 2027-07-12 at MKO: WASP-76 b, WASP-189 b, tau Boo b, WASP-33 b
Making plot for 2027-07-13 at MKO: KELT-20 b
Making plot for 2027-07-14 at MKO: WASP-76 b, KELT-9 b
Making plot for 2027-07-15 at MKO: WASP-189 b, tau Boo b
Making plot for 2027-07-16 at MKO: KELT-20 b
7 plot(s) in /home/you/hrccs_work/plots
```

`--site` accepts a site, a telescope or an instrument name: `keck`, `gemini-n`, `magellan`,
`vlt`, `gemini-s`, and others. Run `hrccs-plan sites` for the full list. Any astropy site name
also works.

## 3. Read the plot

![2027-07-12 at Maunakea](../images/nights_example.png)

- **Top panel:** airmass versus UT. The Sun and Moon are drawn too. Light gray is between
  sunset and sunrise; dark gray is astronomical night (Sun below −18°).
- **Bottom panel:** each planet's velocity relative to its star, shown while the target is
  above 15°. A target appears only if:
  - it is above 30° for part of the night;
  - its velocity changes by more than 30 km/s while it is up;
  - its velocity is *decreasing*, which means the planet is on the far side of its orbit and
    its dayside faces us.

On this night, WASP-33 b passes through secondary eclipse (v = 0) just before sunrise, and
WASP-189 b and τ Boo b are observable early in the night.

## Options you'll use often

- `--date 2027-07-12`: a single night.
- `--require-conjunction`: only planets whose velocity crosses zero, i.e. secondary eclipse,
  during the night.
- `--inclination 30`: scales every Kp by sin i, for non-transiting planets. See
  [tutorial 4](04_nontransiting.md).
- `--outdir DIR`: write the plots somewhere other than `<workspace>/plots`.

Next, [tutorial 2](02_targetlists.md) shows how to build your own target list.
