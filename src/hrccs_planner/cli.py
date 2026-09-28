"""hrccs-plan: command-line interface.

    hrccs-plan init [DIR]                  create a workspace with example target lists
    hrccs-plan sites                       list the built-in sites
    hrccs-plan nights LIST --site S ...    nightly airmass / planet-velocity plots
    hrccs-plan windows LIST --site S ...   rank observing windows for dayside emission
    hrccs-plan transits LIST --site S ...  rank transit windows (transit + baseline)
    hrccs-plan ephemeris LIST --date D     orbital phase and its uncertainty on a date
    hrccs-plan targetlist archive NAME ... target list from the NASA Exoplanet Archive
    hrccs-plan catalog tic ID ...          TIC coordinates and magnitudes
    hrccs-plan starlist keck|magellan ...  telescope starlists (SIMBAD positions and proper motions)

Run `hrccs-plan COMMAND --help` for the options of each command.
"""
import argparse
import csv
import os
import sys

from . import __version__
from .sites import ALIASES, SITES, get_site
from .targets import read_targets, write_targets
from .visibility import dates_between
from .workspace import ENV_VAR, init_workspace, output_path, resolve_targetlist, workspace_dir


def _dates(args, parser):
    if getattr(args, 'date', None):
        if args.start_date or args.end_date:
            parser.error('--date cannot be combined with --start-date/--end-date')
        return [args.date]
    if not (args.start_date and args.end_date):
        parser.error('give --date, or both --start-date and --end-date')
    return dates_between(args.start_date, args.end_date)


def _add_dates(p, single=True):
    if single:
        p.add_argument('--date', help='a single night (local date, YYYY-MM-DD)')
    p.add_argument('--start-date', '--start', dest='start_date', help='first night (local date, YYYY-MM-DD)')
    p.add_argument('--end-date', '--end', dest='end_date', help='last night (local date, YYYY-MM-DD)')


def cmd_init(args, parser):
    root, copied = init_workspace(args.dir or workspace_dir(args.workspace))
    print(f'Workspace: {root}')
    for p in copied:
        print('  copied example', os.path.relpath(p, root))
    print(f'Use it with --workspace {root}, or: export {ENV_VAR}={root}')


def cmd_sites(args, parser):
    for key, site in SITES.items():
        print(f'{key:10s} {site.name:32s} lat {site.location.lat.deg:7.3f}  lon {site.location.lon.deg:8.3f}  '
              f'aliases: {", ".join(ALIASES[key]) or "-"}')
    print('Any astropy site name (EarthLocation.get_site_names()) also works.')


def cmd_windows(args, parser):
    from .windows import WindowOptions, rank_windows, split_path, table_lines, write_windows
    path = resolve_targetlist(args.targetlist, args.workspace)
    targets = read_targets(path)
    site = get_site(args.site)
    opts = WindowOptions(
        duration=args.duration, alt_min=args.alt_min, twilight=args.twilight, delta_v_min=args.delta_v_min,
        bright_threshold=args.bright_threshold, moon_sep_min=args.moon_sep_min, inclination=args.inclination,
        weight_kmag=args.weight_kmag, contrast_column=args.contrast_column if args.weight_contrast else None,
        max_hours=args.max_hours, airmass_k=args.airmass_k, seeing_exponent=args.seeing_exponent,
        normalize_per_target=args.normalize_per_target, snr_min=args.snr_min, sort=args.sort, top=args.top)
    try:
        rows, n_all = rank_windows(targets, site, _dates(args, parser), opts)
    except ValueError as err:
        parser.error(str(err))
    first = f'{len(rows)} of {n_all} windows with relative S/N >= {args.snr_min}'
    print('\n'.join(table_lines(rows, first)))
    if args.csv:
        out = output_path(args.csv, 'output/best_windows', args.workspace)
        print('Saved', ' and '.join(write_windows(out, rows, first)))
        if args.split_targets:
            norm = "this target's best window" if args.normalize_per_target else 'the best window of all targets'
            for name in dict.fromkeys(r['name'] for r in rows):
                sub = [r for r in rows if r['name'] == name]
                files = write_windows(split_path(out, name), sub,
                                      f'{name}: {len(sub)} windows (relative S/N normalized to {norm})')
                print('Saved', ' and '.join(files))
    _write_gemini_tw(args, rows, 'output/best_windows')


def cmd_transits(args, parser):
    from .transits import TransitOptions, rank_transits, table_lines, write_transits
    from .windows import split_path
    targets = read_targets(resolve_targetlist(args.targetlist, args.workspace))
    if args.duration is None and args.duration_column is None:
        parser.error('give the transit duration: --duration HOURS, or --duration-column COLUMN')
    opts = TransitOptions(
        duration=args.duration, duration_column=args.duration_column, baseline=args.baseline,
        baseline_mode=args.baseline_mode, alt_min=args.alt_min, twilight=args.twilight,
        bright_threshold=args.bright_threshold, moon_sep_min=args.moon_sep_min, weight_kmag=args.weight_kmag,
        contrast_column=args.contrast_column if args.weight_contrast else None, airmass_k=args.airmass_k,
        seeing_exponent=args.seeing_exponent, normalize_per_target=args.normalize_per_target,
        snr_min=args.snr_min, sort=args.sort, top=args.top)
    rows, n_all = rank_transits(targets, get_site(args.site), _dates(args, parser), opts)
    first = (f'{len(rows)} of {n_all} transit windows with relative S/N >= {args.snr_min} '
             f'(baseline {args.baseline:g} h, {args.baseline_mode})')
    print('\n'.join(table_lines(rows, first)))
    if args.csv:
        out = output_path(args.csv, 'output/transit_windows', args.workspace)
        print('Saved', ' and '.join(write_transits(out, rows, first)))
        if args.split_targets:
            norm = "this target's best window" if args.normalize_per_target else 'the best window of all targets'
            for name in dict.fromkeys(r['name'] for r in rows):
                sub = [r for r in rows if r['name'] == name]
                print('Saved', ' and '.join(write_transits(split_path(out, name), sub,
                                                           f'{name}: {len(sub)} transit windows (relative S/N normalized to {norm})')))
    _write_gemini_tw(args, rows, 'output/transit_windows')


def _write_gemini_tw(args, rows, kind):
    if not args.gemini_tw:
        return
    from .windows import gemini_timing_windows
    out = output_path(args.gemini_tw, kind, args.workspace)
    text = gemini_timing_windows(rows, args.tw_round)
    with open(out, 'w') as fh:
        fh.write(text)
    print(f'\nGemini PIT timing windows (paste into the Scheduling field), saved to {out}:\n')
    print(text, end='')


def cmd_nights(args, parser):
    import matplotlib

    from .nights import NightOptions, plan_nights
    matplotlib.use('Agg')
    targets = read_targets(resolve_targetlist(args.targetlist, args.workspace))
    opts = NightOptions(alt_limit=args.alt_limit, delta_v_min=args.delta_v_min,
                        require_conjunction=args.require_conjunction, moon_sep_min=args.moon_sep_min,
                        inclination=args.inclination)
    outdir = args.outdir or os.path.join(workspace_dir(args.workspace), 'plots')
    from .nights import parse_window, read_shading
    shading = {}
    if args.shade:
        paths = [p if os.path.dirname(p) or os.path.exists(p) else
                 next((c for c in (os.path.join(workspace_dir(args.workspace), 'output', d, p)
                                   for d in ('best_windows', 'transit_windows')) if os.path.exists(c)), p)
                 for p in args.shade]
        shading = read_shading(paths)
    if args.window:
        shading[None] = [parse_window(w) + (None, None) for w in args.window]
    saved = plan_nights(targets, get_site(args.site), _dates(args, parser), outdir, opts, shading=shading,
                        compact=args.compact, fmt=args.format, title=args.title, site_label=args.site_label)
    print(f'{len(saved)} plot(s) in {outdir}')


def cmd_ephemeris(args, parser):
    from astropy.time import Time

    from .ephemeris import propagate_conjunction
    targets = read_targets(resolve_targetlist(args.targetlist, args.workspace))
    t = Time(args.date if 'T' in args.date else args.date + 'T12:00:00', scale='utc').jd
    print(f"{'target':16s} {'n':>6s} {'T_conj (JD)':>14s} {'sigma_t [h]':>11s} {'phase':>6s} {'sigma_phase':>11s}")
    for tg in targets:
        t0_err = args.t0_err if args.t0_err is not None else float(tg.extra.get(args.t0_err_column) or 0)
        p_err = args.period_err if args.period_err is not None else float(tg.extra.get(args.period_err_column) or 0)
        n, tn, st, sp = propagate_conjunction(tg.t0, t0_err, tg.period, p_err, t)
        phase = ((t - tg.t0) / tg.period) % 1
        print(f'{tg.name[:16]:16s} {n:6.0f} {tn:14.5f} {st * 24:11.2f} {phase:6.3f} {sp:11.3f}')
    if args.t0_err is None and args.period_err is None:
        print(f"(errors from the '{args.t0_err_column}' and '{args.period_err_column}' columns; 0 if absent)")


def cmd_targetlist_archive(args, parser):
    from .catalogs import exoarchive_target, query_exoarchive
    rows = query_exoarchive(args.names, args.workspace, args.refresh)
    name = args.output or '_'.join(n.replace(' ', '') for n in args.names)
    name = name[:-4] if name.endswith('.csv') else name
    name = name[:-11] if name.endswith('_targetlist') else name
    out = output_path(name + '_targetlist.csv', 'targetlists', args.workspace)
    write_targets(out, [exoarchive_target(r, args.circular_below) for r in rows])
    print(f'Wrote {len(rows)} planet(s) to {out}')


def cmd_catalog_tic(args, parser):
    from .catalogs import query_tic
    rows = query_tic(args.tic, args.workspace, args.refresh)
    out = output_path(args.output, 'output', args.workspace)
    cols = ['ID', 'ra', 'dec', 'Tmag', 'Vmag', 'Kmag', 'Teff', 'rad', 'mass', 'd']
    with open(out, 'w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(['TIC_ID', 'RA_deg', 'Dec_deg', 'Tmag', 'Vmag', 'Kmag', 'Teff', 'Rstar', 'Mstar', 'dist_pc'])
        for tic in args.tic:
            r = rows.get(int(tic))
            w.writerow([tic] + (['' if r is None else r[c] for c in cols[1:]]))
            if r is None:
                print('Not found in the TIC:', tic)
    print(f'Wrote {len(args.tic)} row(s) to {out}')


def cmd_starlist(args, parser):
    from .starlists import write_starlist
    if args.names:
        names = args.names
    elif args.targetlist:
        names = list(dict.fromkeys(t.host for t in read_targets(resolve_targetlist(args.targetlist, args.workspace))))
    else:
        parser.error('give a target list (--list) or --names')
    out = output_path(args.output, 'starlists', args.workspace)
    print(write_starlist(out, names, args.format, args.workspace, args.refresh), end='')
    print('Saved', out)


def build_parser():
    parser = argparse.ArgumentParser(prog='hrccs-plan', description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--version', action='version', version=f'hrccs_planner {__version__}')
    common = argparse.ArgumentParser(add_help=False)
    common.add_argument('--workspace', help=f'workspace directory (default: ${ENV_VAR} or the current directory)')
    sub = parser.add_subparsers(dest='command', required=True, metavar='COMMAND')

    p = sub.add_parser('init', parents=[common], help='create a workspace with example target lists')
    p.add_argument('dir', nargs='?', help='workspace directory (default: --workspace, or the current directory)')
    p.set_defaults(func=cmd_init)

    p = sub.add_parser('sites', parents=[common], help='list the built-in sites')
    p.set_defaults(func=cmd_sites)

    p = sub.add_parser('windows', parents=[common], help='rank observing windows for dayside emission',
                       description='Rank nights (and the window within each night) for dayside emission '
                                   'HRCCS of circular-orbit planets. See docs/method.md.')
    p.add_argument('targetlist', help='target-list name (<workspace>/targetlists/<name>_targetlist.csv) or path')
    p.add_argument('--site', required=True, help='site, telescope or instrument (hrccs-plan sites)')
    _add_dates(p, single=False)
    p.add_argument('--duration', type=float, default=0., help='eclipse duration [h] excluded around phase 0.5 (default 0)')
    p.add_argument('--alt-min', type=float, default=30., help='minimum target altitude [deg] (default 30, airmass 2)')
    p.add_argument('--twilight', type=float, default=-18., help='Sun altitude defining night [deg] (default -18)')
    p.add_argument('--delta-v-min', type=float, default=30., help='minimum planet velocity change [km/s] (default 30)')
    p.add_argument('--bright-threshold', type=float, default=0.5, help='Moon illumination at mid-night above which a night is bright time (default 0.5)')
    p.add_argument('--moon-sep-min', type=float, default=10., help='minimum target-Moon separation [deg] while the target is up at night (default 10)')
    p.add_argument('--inclination', type=float, default=90., help='orbital inclination [deg] for all targets: Kp = Kp_max sin(i), and the phase curve (default 90)')
    p.add_argument('--max-hours', type=float, help='maximum wall-clock length of a window [h]: the best stretch of this length is used')
    p.add_argument('--weight-kmag', action='store_true', help='weight S/N by the stellar K-band photon flux, 10^(-0.2 (Kmag - 6)); targets without Kmag are skipped')
    p.add_argument('--weight-contrast', action='store_true', help='weight S/N by the planet contrast in --contrast-column')
    p.add_argument('--contrast-column', default='K contrast (ppm)', help="target-list column used by --weight-contrast (default 'K contrast (ppm)')")
    p.add_argument('--airmass-k', type=float, default=0., help='extinction [mag/airmass] in the airmass weight 10^(-0.4 k (X-1)) on S/N^2 (default 0; ~0.05 in H/K at Maunakea)')
    p.add_argument('--seeing-exponent', type=float, default=0., help='slit-loss weight X^(-p) on S/N^2, for seeing (FWHM ~ X^0.6) wider than the slit (default 0; 0.6 in that limit)')
    p.add_argument('--normalize-per-target', action='store_true', help="normalize S/N to each target's own best window (--snr-min then applies per target)")
    p.add_argument('--snr-min', type=float, default=0.3, help='drop windows with relative S/N below this (default 0.3; 0 keeps all)')
    p.add_argument('--sort', choices=['moon', 'snr'], default='moon', help='moon: bright time first, then S/N (default); snr: by S/N only')
    p.add_argument('--top', type=int, help='keep only the top N windows')
    p.add_argument('--csv', help='save the windows to this CSV (in <workspace>/output/best_windows/ unless a directory is given), and the table to the same name with .txt')
    p.add_argument('--split-targets', action='store_true', help='with --csv, also write one CSV/.txt per target (<csv>_<target>.csv)')
    p.add_argument('--gemini-tw', metavar='FILE', help='also write the windows as Gemini PIT timing windows (UT start and duration per line, one block per target) for the proposal Scheduling field (in <workspace>/output/best_windows/ unless a directory is given)')
    p.add_argument('--tw-round', type=int, default=15, metavar='MIN', help='round timing-window starts down and ends up to MIN minutes (default 15; 0 for exact minutes)')
    p.set_defaults(func=cmd_windows)

    p = sub.add_parser('transits', parents=[common], help='rank transit windows (transit + out-of-transit baseline)',
                       description='Rank transit windows for transmission HRCCS: the transit (T14) plus an '
                                   'out-of-transit baseline, all at night and above --alt-min. The S/N counts the '
                                   'in-transit time weighted by airmass. See docs/method.md.')
    p.add_argument('targetlist', help='target-list name (<workspace>/targetlists/<name>_targetlist.csv) or path')
    p.add_argument('--site', required=True, help='site, telescope or instrument (hrccs-plan sites)')
    _add_dates(p, single=False)
    p.add_argument('--duration', type=float, help='transit duration T14 [h], first to fourth contact, for all targets')
    p.add_argument('--duration-column', help="target-list column with T14 [h] (overrides --duration where filled; 'T14 (h)' is written by targetlist archive)")
    p.add_argument('--baseline', type=float, default=2., help='total out-of-transit baseline [h] (default 2)')
    p.add_argument('--baseline-mode', choices=['split', 'any'], default='split',
                   help='split: half the baseline before and half after the transit (default); '
                        'any: any division, including all on one side (the most even one that fits is used)')
    p.add_argument('--alt-min', type=float, default=30., help='minimum altitude [deg] throughout the window (default 30, airmass 2)')
    p.add_argument('--twilight', type=float, default=-18., help='Sun altitude defining night [deg] (default -18)')
    p.add_argument('--bright-threshold', type=float, default=0.5, help='Moon illumination above which a night is bright time (default 0.5)')
    p.add_argument('--moon-sep-min', type=float, default=10., help='minimum target-Moon separation [deg] (default 10)')
    p.add_argument('--weight-kmag', action='store_true', help='weight S/N by 10^(-0.2 (Kmag - 6)); targets without Kmag are skipped')
    p.add_argument('--weight-contrast', action='store_true', help='weight S/N by the value in --contrast-column (e.g. an expected transmission signal)')
    p.add_argument('--contrast-column', default='K contrast (ppm)', help="column used by --weight-contrast (default 'K contrast (ppm)')")
    p.add_argument('--airmass-k', type=float, default=0., help='extinction [mag/airmass] in the in-transit airmass weight on S/N^2 (default 0)')
    p.add_argument('--seeing-exponent', type=float, default=0., help='slit-loss weight X^(-p) on S/N^2 (default 0; 0.6 for seeing wider than the slit)')
    p.add_argument('--normalize-per-target', action='store_true', help="normalize S/N to each target's own best window")
    p.add_argument('--snr-min', type=float, default=0.3, help='drop windows with relative S/N below this (default 0.3)')
    p.add_argument('--sort', choices=['moon', 'snr'], default='moon', help='moon: bright time first, then S/N (default); snr: by S/N only')
    p.add_argument('--top', type=int, help='keep only the top N windows')
    p.add_argument('--csv', help='save to this CSV (in <workspace>/output/transit_windows/ unless a directory is given) and the table to .txt')
    p.add_argument('--split-targets', action='store_true', help='with --csv, also write one CSV/.txt per target')
    p.add_argument('--gemini-tw', metavar='FILE', help='also write the windows (transit + baseline) as Gemini PIT timing windows for the proposal Scheduling field')
    p.add_argument('--tw-round', type=int, default=15, metavar='MIN', help='round timing-window starts down and ends up to MIN minutes (default 15; 0 for exact)')
    p.set_defaults(func=cmd_transits)

    p = sub.add_parser('nights', parents=[common], help='nightly airmass and planet-velocity plots',
                       description='Plot, for each night, the airmass and planet velocity of every target '
                                   'observable on its dayside. See docs/method.md.')
    p.add_argument('targetlist', help='target-list name or path')
    p.add_argument('--site', required=True, help='site, telescope or instrument (hrccs-plan sites)')
    _add_dates(p)
    p.add_argument('--outdir', help='directory for the plots (default <workspace>/plots)')
    p.add_argument('--alt-limit', type=float, default=30., help='altitude a target must exceed at night [deg] (default 30)')
    p.add_argument('--delta-v-min', type=float, default=30., help='minimum planet velocity change [km/s] (default 30)')
    p.add_argument('--require-conjunction', action='store_true', help='only targets whose velocity crosses zero (conjunction) during the night')
    p.add_argument('--moon-sep-min', type=float, default=10., help='minimum target-Moon separation [deg] (default 10)')
    p.add_argument('--inclination', type=float, default=90., help='orbital inclination [deg] for all targets (default 90)')
    p.add_argument('--shade', nargs='+', metavar='CSV', help='shade the observing windows in these `windows`/`transits` CSVs (matched by night and target; bare names are looked up in <workspace>/output/)')
    p.add_argument('--window', nargs='+', metavar='HH:MM-HH:MM', help='shade this UT window on every plotted night (e.g. 07:49-14:33)')
    p.add_argument('--compact', action='store_true', help='small figure for proposals/papers (6.4 x 5.2 in, shared UT axis, dark hours only, airmass-2 line)')
    p.add_argument('--format', choices=['png', 'pdf', 'svg'], default='png', help='output format (default png)')
    p.add_argument('--title', help="figure title (default generated; '' for none)")
    p.add_argument('--site-label', help='site name for the compact title (default the site registry name, e.g. "Gemini North")')
    p.set_defaults(func=cmd_nights)

    p = sub.add_parser('ephemeris', parents=[common], help='orbital phase and its uncertainty on a date')
    p.add_argument('targetlist', help='target-list name or path')
    p.add_argument('--date', required=True, help='UTC date or time (YYYY-MM-DD[THH:MM:SS]); a date means 12:00 UTC')
    p.add_argument('--t0-err', type=float, help='T0 uncertainty [d] for all targets')
    p.add_argument('--period-err', type=float, help='period uncertainty [d] for all targets')
    p.add_argument('--t0-err-column', default='T err (d)', help="column with the T0 uncertainty (default 'T err (d)')")
    p.add_argument('--period-err-column', default='P err (d)', help="column with the period uncertainty (default 'P err (d)')")
    p.set_defaults(func=cmd_ephemeris)

    p = sub.add_parser('targetlist', help='make target lists from catalogs')
    tsub = p.add_subparsers(dest='source', required=True, metavar='SOURCE')
    q = tsub.add_parser('archive', parents=[common], help='from the NASA Exoplanet Archive (pscomppars)')
    q.add_argument('names', nargs='+', help='host or planet names, e.g. WASP-121 "HD 189733"')
    q.add_argument('-o', '--output', help='target-list name (default: the names joined by _)')
    q.add_argument('--circular-below', type=float, default=0., help='treat orbits with e below this as circular (default 0: any e > 0 is eccentric; `windows` skips eccentric orbits)')
    q.add_argument('--refresh', action='store_true', help='query again instead of using the cache')
    q.set_defaults(func=cmd_targetlist_archive)

    p = sub.add_parser('catalog', help='catalog look-ups')
    csub = p.add_subparsers(dest='catalog', required=True, metavar='CATALOG')
    q = csub.add_parser('tic', parents=[common], help='TIC-8 coordinates, magnitudes and stellar parameters')
    q.add_argument('tic', nargs='+', type=int, help='TIC IDs')
    q.add_argument('-o', '--output', default='tic_data.csv', help='output CSV (default <workspace>/output/tic_data.csv)')
    q.add_argument('--refresh', action='store_true', help='query again instead of using the cache')
    q.set_defaults(func=cmd_catalog_tic)

    p = sub.add_parser('starlist', parents=[common], help='telescope starlist from SIMBAD')
    p.add_argument('format', choices=['keck', 'magellan'])
    p.add_argument('--list', dest='targetlist', help='target list whose host stars to include')
    p.add_argument('--names', nargs='+', help='star names (instead of --list)')
    p.add_argument('-o', '--output', required=True, help='output file (in <workspace>/starlists/ unless a directory is given)')
    p.add_argument('--refresh', action='store_true', help='query again instead of using the cache')
    p.set_defaults(func=cmd_starlist)
    return parser


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        args.func(args, parser)
    except FileNotFoundError as err:
        parser.exit(2, f'hrccs-plan: error: {err}\n')


def legacy(command, argv):
    """Entry point for the old scripts: `SCRIPT SITE LIST [options]`."""
    old = {'windows': 'best_window.py', 'nights': 'observing_planner.py'}[command]
    print(f'Note: {old} is deprecated; use `hrccs-plan {command} LIST --site SITE ...`', file=sys.stderr)
    if len(argv) < 2 or argv[0].startswith('-'):
        main([command, '--help'])
    site, targetlist, rest = argv[0], argv[1], argv[2:]
    main([command, targetlist, '--site', site] + rest)
