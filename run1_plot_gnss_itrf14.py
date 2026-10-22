#!/usr/bin/env python3
"""Run 1: inspect GNSS station coverage and velocities in global ITRF2014.

Usage:
    python3 run1_plot_gnss_itrf14.py
    python3 run1_plot_gnss_itrf14.py --region 67 113 22 45
    python3 run1_plot_gnss_itrf14.py --show

Uses the author's original v6.4 _itrf14.dat files, including horizontal
components from both 2D and 3D GNSS records. The map extent, symbol scaling,
and styling match run1_plot_gnss.py for direct visual comparison.

Output: GNSS_analysis/output/run1_itrf14/run1_gnss_itrf14.png and .pdf.
The shared plotting implementation lives in run1_plot_gnss.py; keep both
scripts in the same folder. No interpolation or frame transformation is done.
"""

from pathlib import Path

from run1_plot_gnss import main


if __name__ == "__main__":
    main(frame="itrf14", script_name=Path(__file__).name)
