"""The workspace: the directory holding your target lists and outputs.

    <workspace>/
        targetlists/<name>_targetlist.csv   target lists
        starlists/                          telescope starlists
        output/best_windows/                ranked observing windows
        plots/                              nightly planning plots
        cache/                              cached catalog queries

The workspace is, in order of precedence: the --workspace option, the
HRCCS_WORKSPACE environment variable, or the current directory.
"""
import glob
import os
import shutil

ENV_VAR = 'HRCCS_WORKSPACE'
SUBDIRS = ['targetlists', 'starlists', 'output/best_windows', 'plots', 'cache']
EXAMPLES = os.path.join(os.path.dirname(__file__), 'data', 'targetlists')


def workspace_dir(path=None):
    return os.path.abspath(os.path.expanduser(path or os.environ.get(ENV_VAR) or os.getcwd()))


def subdir(name, path=None, create=False):
    d = os.path.join(workspace_dir(path), name)
    if create:
        os.makedirs(d, exist_ok=True)
    return d


def resolve_targetlist(name, path=None):
    """Path of a target list given as a file path or as a name
    (<workspace>/targetlists/<name>_targetlist.csv, or a bundled example)."""
    if os.path.isfile(name):
        return name
    for d in (subdir('targetlists', path), EXAMPLES):
        for cand in (f'{name}_targetlist.csv', f'{name}.csv', name):
            p = os.path.join(d, cand)
            if os.path.isfile(p):
                return p
    have = sorted(os.path.basename(p)[:-len('_targetlist.csv')]
                  for p in glob.glob(os.path.join(subdir('targetlists', path), '*_targetlist.csv')))
    examples = sorted(os.path.basename(p)[:-len('_targetlist.csv')]
                      for p in glob.glob(os.path.join(EXAMPLES, '*_targetlist.csv')))
    raise FileNotFoundError(
        f"no target list '{name}' in {subdir('targetlists', path)} (workspace: {workspace_dir(path)}; "
        f"set --workspace or ${ENV_VAR}). Available: {', '.join(have) or 'none'}; "
        f"examples: {', '.join(examples)}")


def output_path(filename, kind, path=None):
    """`filename` unchanged if it includes a directory, else <workspace>/<kind>/filename."""
    if os.path.dirname(filename):
        out = filename
    else:
        out = os.path.join(subdir(kind, path), filename)
    os.makedirs(os.path.dirname(os.path.abspath(out)), exist_ok=True)
    return out


def init_workspace(path, examples=True):
    """Create the workspace directories and copy the example target lists."""
    root = workspace_dir(path)
    for d in SUBDIRS:
        os.makedirs(os.path.join(root, d), exist_ok=True)
    copied = []
    if examples:
        for src in glob.glob(os.path.join(EXAMPLES, '*.csv')):
            dst = os.path.join(root, 'targetlists', os.path.basename(src))
            if not os.path.exists(dst):
                shutil.copy(src, dst)
                copied.append(dst)
    return root, copied
