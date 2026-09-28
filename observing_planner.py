"""Deprecated: use `hrccs-plan nights LIST --site SITE ...` (same options).

    python observing_planner.py SITE LIST [options]  ==  hrccs-plan nights LIST --site SITE [options]
"""
import sys

from hrccs_planner.cli import legacy

if __name__ == '__main__':
    legacy('nights', sys.argv[1:])
