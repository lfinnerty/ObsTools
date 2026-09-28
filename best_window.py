"""Deprecated: use `hrccs-plan windows LIST --site SITE ...` (same options).

    python best_window.py SITE LIST [options]  ==  hrccs-plan windows LIST --site SITE [options]
"""
import sys

from hrccs_planner.cli import legacy

if __name__ == '__main__':
    legacy('windows', sys.argv[1:])
