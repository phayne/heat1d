"""Export utilities for figures and data."""

import numpy as np
from PySide6.QtWidgets import QFileDialog


def export_figure(fig, parent=None):
    """Save a matplotlib figure via file dialog.

    Parameters
    ----------
    fig : matplotlib.figure.Figure
    parent : QWidget or None
    """
    path, filt = QFileDialog.getSaveFileName(
        parent,
        "Export Figure",
        "",
        "PNG Image (*.png);;PDF Document (*.pdf);;SVG Vector (*.svg)",
    )
    if not path:
        return
    fig.savefig(path, dpi=300, bbox_inches="tight")


def export_temperature_csv(run_record, parent=None):
    """Export temperature data for a single run as CSV.

    Parameters
    ----------
    run_record : RunRecord
    parent : QWidget or None
    """
    path, _ = QFileDialog.getSaveFileName(
        parent, "Export Data", f"{run_record.label}.csv", "CSV Files (*.csv)"
    )
    if not path:
        return
    model = run_record.model
    header = "local_time_hr," + ",".join(
        f"z={z:.4f}m" for z in model.profile.z
    )
    data = np.column_stack([model.lt, model.T])
    np.savetxt(path, data, delimiter=",", header=header, comments="")


def export_multi_csv(run_records, parent=None):
    """Export temperature data for multiple runs.

    Each run gets its own section with a comment header.
    """
    path, _ = QFileDialog.getSaveFileName(
        parent, "Export Data", "heat1d_comparison.csv", "CSV Files (*.csv)"
    )
    if not path:
        return

    with open(path, "w") as f:
        for rec in run_records:
            model = rec.model
            f.write(f"# Run: {rec.label}\n")
            f.write(f"# Planet: {rec.planet_name}, Lat: {rec.lat_deg} deg\n")
            header = "local_time_hr," + ",".join(
                f"z={z:.4f}m" for z in model.profile.z
            )
            f.write(header + "\n")
            data = np.column_stack([model.lt, model.T])
            for row in data:
                f.write(",".join(f"{v:.6g}" for v in row) + "\n")
            f.write("\n")
