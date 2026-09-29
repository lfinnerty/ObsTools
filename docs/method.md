# Method

## Nights and time grids

- **Night labels.** A night is labelled by the local date on which it starts, which is the
  HST date at Maunakea.
- **The grid.** Times are sampled over 24 h. `windows` uses 1-minute steps and `nights` uses
  100 points. The grid starts at `date 23:59:59 UTC + h0`, with h0 = −3 h at Maunakea and
  −6 h at the Chilean sites, so it covers the whole local night. For other astropy sites, h0
  starts the grid at local noon.
- **Sun and Moon.** Positions come from astropy (`get_sun`, `get_body`).
- **Airmass** is 1/sin(alt).

## `windows`: ranking nights for dayside emission

For each target and night, the **usable time** is the set of samples where all of the
following hold:

- the target is above `--alt-min` (default 30°, airmass 2), and above the site's pointing
  limit if it has one (`--site keck2` applies the Keck II Nasmyth-deck limit);
- the Sun is below `--twilight` (default −18°, astronomical twilight);
- the planet is on its dayside half, orbital phase 0.25–0.75;
- the planet is not in secondary eclipse, i.e. |phase − 0.5| is at least half of `--duration`.

**Eccentric orbits** (`Eccentric? = Yes`, with `e` and `omega`) use the orbital phase u/2π,
where u = f + ω − π/2 is the orbital angle from transit (true anomaly f, from Kepler's
equation): the dayside is u/2π = 0.25–0.75, the phase curve uses it too, and the phases in
the output are orbital phases. The eclipse is excluded within `--duration`/2 of its own time,
which is generally not half an orbit after transit. The velocity is the Keplerian curve,
v = −Kp (cos(f + ω) + e cos ω)/√(1 − e²) with Kp = Kp_max sin i (docs/targetlist_format.md). The phase curve is geometric: the dayside flux doesn't follow the changing
irradiation. Very eccentric planets are on their dayside only briefly (HD 80606 b: 1.8 d of
its 111 d orbit, around periastron).

Nights are dropped if either of these holds:

- the Moon comes within `--moon-sep-min` (10°) while the target is up at night;
- the planet velocity changes by less than `--delta-v-min` (30 km/s) over the usable time.

### Relative S/N

For a photon-noise-limited observation, the matched-filter S/N of the planet signal is
estimated as

S/N ∝ √(Σ f² w Δt · Δv) × W

where each term is:

- **f** is the Lambertian phase curve, normalized to 1 at full phase:
  f = (sin α + (π − α) cos α)/π, with the phase angle α given by
  cos α = −sin i cos(2π phase).
- **w** is the optional airmass weight on S/N², w = 10^(−0.4 k (X − 1)) · X^(−p):
  - `--airmass-k k` sets the extinction in mag/airmass (≈ 0.05 in H and K at Maunakea);
  - `--seeing-exponent p` sets the slit losses. When the seeing (FWHM ∝ X^0.6) is wider than
    the slit, the fraction of light through the slit scales as 1/FWHM, so p = 0.6.
- **Δv** is the range of planet velocity covered. A larger Δv moves the planet lines across
  more pixels relative to the stellar and telluric lines, so less of the planet signal is
  removed in the detrending.
- **W** is an optional per-target weight:
  - `--weight-kmag` gives 10^(−0.2 (K − 6)), the photon flux relative to K = 6;
  - `--weight-contrast` multiplies by a contrast column, e.g. the expected dayside
    planet/star flux ratio.

The S/N is normalized to the best window overall. With `--normalize-per-target`, it is
normalized to each target's own best window instead, and `--snr-min` then applies per target.

This is a relative score: use it to compare nights and, with the weights, targets. For
absolute S/N predictions, simulate the observation (e.g. with KPIC-Synthgen; see
[tutorial 5](tutorials/05_starlists_and_simulations.md)).

### Capping the window length

`--max-hours H` limits each night's window to at most H hours of clock time, for targets that
are detected well enough without the whole night.

- **How the stretch is chosen.** It is the contiguous stretch of length ≤ H, within the usable
  time, that maximizes Σ f² w Δt · Δv. Gaps inside it, such as an excluded eclipse, still
  count toward H.
- **Where it lands.** Because f peaks at phase 0.5, capped windows usually centre on
  secondary eclipse, trading some Δv for the brightest part of the orbit.

### Output

For each window, the output gives:

- the date, Moon phase and illumination;
- the usable hours and UT start/end;
- the phase and planet velocity at the start and end, and Δv;
- the mean f;
- the mean and maximum airmass;
- the relative S/N;
- whether the window crosses secondary eclipse.

With `--gemini-tw FILE`, the windows are also written as Gemini PIT timing windows (UT start and
duration, one block per target) for the proposal's Scheduling field.

Rows are sorted with bright time first (IR observations are insensitive to moonlight, so
bright nights are easier to get), then by S/N. `--sort snr` sorts by S/N only.

## `phase-sensitivity`: what ephemeris uncertainty costs a window

A window is planned from the predicted ephemeris and observed at fixed clock times. If the true
conjunction is off by δ in phase, the observation covers phases shifted by δ. When the HRCCS
analysis refits the conjunction time, δ changes only *what is covered*: the phase curve f and
the planet-velocity range Δv. For a planned window the S/N becomes

S/N(δ) ∝ √(Σ f(phase + δ)² w Δt · Δv(phase + δ)),

evaluated on the window's own usable samples, with the same weights as `windows`.

- **Which windows.** Each target's `--top-windows` best windows (default 5), ranked as
  `windows` ranks them with per-target normalization. For eccentric orbits δ shifts the time
  since transit, and the dayside, eclipse, phase curve and velocities follow the orbit as in
  `windows`.
- **Averaging.** S/N(δ)/S/N(0) is averaged over δ ~ N(0, σ_phase). The command reports the mean
  and the 10th percentile, i.e. the S/N of an unlucky ephemeris, both averaged over the windows.
- **Where σ comes from.** Either the T0 and period uncertainty columns (`T err (d)`,
  `P err (d)`), propagated to each window as in `ephemeris`, or any number of
  `--sigma-column` columns giving σ_t in hours at the observing epoch. Each column is a
  scenario, e.g. the current ephemeris and one improved by new RVs; the last column of the
  output is the change in mean S/N relative to the first scenario.
- **What to expect.** Windows centred on secondary eclipse sit where f and the velocity slope
  both peak, so a small δ costs S/N only at second order: roughly ½ (S″/S) σ_phase². Losses
  become noticeable once σ_t is a sizeable fraction of the window (σ_phase ≳ 0.03). Windows
  off the eclipse are first-order sensitive, which widens the spread (the 10th percentile)
  more than it lowers the mean.
- **Assumptions.** The S/N is photon-noise limited, and the analysis re-derives the conjunction
  (a fixed, wrong ephemeris would also misplace the planet's velocity track). The inclination
  (`--inclination`) sets how strongly f varies; lower inclinations make phase errors matter less.

## `transits`: ranking transit windows

A **transit window** is the transit, from first to fourth contact (`--duration`, T14 in hours),
centred on T_c = T0 + n P. It is extended by an out-of-transit baseline of `--baseline` hours
in total:

- **`--baseline-mode split`** (the default) puts half the baseline before ingress and half
  after egress.
- **`--baseline-mode any`** accepts any division, including all of the baseline on one side.
  The most even division that fits is used; among equally even ones, the one with the lowest
  mean airmass.

The whole window, transit plus baseline, must meet these conditions:

- it is at night (Sun below `--twilight`);
- the target is above `--alt-min` (30°, i.e. airmass < 2) and the site's pointing limit;
- the Moon stays at least `--moon-sep-min` from the target.

The relative S/N counts only the in-transit time, weighted by airmass:

S/N ∝ √(Σ_in-transit w(X) Δt) × W

- **w(X)** is the same airmass weight as for emission, set by `--airmass-k` and
  `--seeing-exponent`.
- **W** is the optional per-target weight: `--weight-kmag`, and `--weight-contrast` with a
  column such as an expected transmission signal.
- **No phase curve or velocity term** is used. A transit is always seen at the same phases, so
  the planet's illuminated fraction doesn't enter.

Transit times use only T0 and the period, so eccentric orbits work as long as T0 is a transit
time. Targets with `Transiting? = N` are skipped. Per-target durations can come from a column
(`--duration-column`), which overrides `--duration` where it is filled in.

`--gemini-tw` writes each window, transit plus baseline, as a Gemini timing window.

## `nights`: nightly plots

For each night, `nights` plots every target that meets all of the following:

- the Moon stays at least `--moon-sep-min` away while the target is up at night;
- the target is above `--alt-limit` (30°) for part of the night (Sun below 0°);
- its velocity changes by more than `--delta-v-min` while it is above that limit;
- its velocity decreases over the night, i.e. the planet is on its dayside, around secondary
  eclipse.

The top panel shows airmass. The bottom panel shows the planet velocity wherever the target
is above 15°. Twilight is shaded light gray (Sun below 0°) and dark gray (Sun below −18°).
With `--require-conjunction`, only targets whose velocity crosses zero during the night are
shown.

`--shade CSV` or `--window HH:MM-HH:MM` shades observing windows. Targets with a shaded window
are always shown, with the velocity and dayside criteria skipped. `--compact` draws a small
publication-style figure on a 2-minute grid (instead of 100 points), so the twilight and window
edges are exact.

## `ephemeris`: phase uncertainty

- **What it computes.** The nearest conjunction to a date and its uncertainty,
  σ_t = √(σ_T0² + (n σ_P)²) after n orbits, with σ_phase = σ_t/P.
- **Where the errors come from.** Use `--t0-err` and `--period-err`, or per-target columns
  (default names `T err (d)` and `P err (d)`).
- **When it matters.** For non-transiting planets with photometric ephemerides, n σ_P usually
  dominates by the time of the observation. If σ_phase is comparable to a capped window's
  phase coverage, the window may miss the phases it was planned for.
