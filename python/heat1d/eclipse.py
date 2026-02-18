"""Eclipse geometry for satellite bodies.

Computes the fraction of the solar disk obscured by a parent body
(e.g., Earth eclipsing the Sun as seen from the Moon) using
circle-circle overlap geometry.  All functions are pure math with
no network or API dependencies.

The eclipse fraction is applied as a multiplicative reduction to
the absorbed flux: ``flux *= (1 - eclipse_fraction)``.
"""

import numpy as np


# ---------------------------------------------------------------------------
# Satellite detection
# ---------------------------------------------------------------------------

def is_satellite(body_id):
    """Determine if a Horizons body ID represents a natural satellite.

    Horizons convention: satellite IDs are ``xNN`` where ``x`` is the
    planet number (1-9) and ``NN`` is 01-98.  ``x99`` is the planet
    barycenter, not a satellite.

    Parameters
    ----------
    body_id : str or int
        Horizons body ID.

    Returns
    -------
    bool
    """
    try:
        bid = int(body_id)
    except (ValueError, TypeError):
        return False
    if bid < 100 or bid >= 1000:
        return False
    remainder = bid % 100
    return 0 < remainder < 99


def get_parent_body_id(satellite_id):
    """Derive the parent body Horizons ID from a satellite ID.

    Examples: 301 (Moon) → 399 (Earth), 502 (Europa) → 599 (Jupiter),
    606 (Titan) → 699 (Saturn).

    Parameters
    ----------
    satellite_id : str or int
        Horizons satellite body ID.

    Returns
    -------
    str
        Parent body Horizons ID.
    """
    bid = int(satellite_id)
    return str((bid // 100) * 100 + 99)


# ---------------------------------------------------------------------------
# Circle-circle overlap geometry
# ---------------------------------------------------------------------------

def circle_overlap_area(r1, r2, d):
    """Area of intersection of two circles.

    Parameters
    ----------
    r1 : float or np.ndarray
        Radius of the first circle.
    r2 : float or np.ndarray
        Radius of the second circle.
    d : float or np.ndarray
        Distance between circle centers.

    Returns
    -------
    float or np.ndarray
        Overlap area.
    """
    r1 = np.asarray(r1, dtype=float)
    r2 = np.asarray(r2, dtype=float)
    d = np.asarray(d, dtype=float)

    # Broadcast all inputs to common shape
    r1, r2, d = np.broadcast_arrays(r1, r2, d)

    # Output array
    area = np.zeros_like(d)

    # Case 1: No overlap
    no_overlap = d >= r1 + r2

    # Case 2: One circle fully inside the other
    r_min = np.minimum(r1, r2)
    small_in_large = (~no_overlap) & (d + r_min <= np.maximum(r1, r2))
    area[small_in_large] = np.pi * r_min[small_in_large] ** 2

    # Case 3: Partial overlap (lens formula)
    partial = ~no_overlap & ~small_in_large
    if np.any(partial):
        dp = d[partial]
        r1p = r1[partial]
        r2p = r2[partial]

        # Clamp arguments to arccos to [-1, 1] for numerical safety
        arg1 = (dp ** 2 + r1p ** 2 - r2p ** 2) / (2.0 * dp * r1p)
        arg2 = (dp ** 2 + r2p ** 2 - r1p ** 2) / (2.0 * dp * r2p)
        arg1 = np.clip(arg1, -1.0, 1.0)
        arg2 = np.clip(arg2, -1.0, 1.0)

        # Lens area formula
        term1 = r1p ** 2 * np.arccos(arg1)
        term2 = r2p ** 2 * np.arccos(arg2)

        # Heron-like term under the square root
        s = (-dp + r1p + r2p) * (dp + r1p - r2p) * (dp - r1p + r2p) * (dp + r1p + r2p)
        s = np.maximum(s, 0.0)  # guard against tiny negatives from rounding
        term3 = 0.5 * np.sqrt(s)

        area[partial] = term1 + term2 - term3

    return area


def compute_eclipse_fraction(toi_deg, sun_ang_diam_arcsec,
                             ib_ang_diam_arcsec):
    """Compute fraction of the solar disk obscured by the interfering body.

    Uses the T-O-I (Target-Observer-Interfering body) angle from JPL
    Horizons to determine the angular separation between the Sun and
    parent body centers, then computes the circle-circle overlap
    fraction.

    Parameters
    ----------
    toi_deg : np.ndarray
        T-O-I angle [degrees] from Horizons (Quantity 25).
        Negative values indicate the Sun center is behind the IB.
    sun_ang_diam_arcsec : np.ndarray
        Sun angular diameter [arcseconds] (Quantity 13).
    ib_ang_diam_arcsec : np.ndarray
        Interfering body angular diameter [arcseconds].

    Returns
    -------
    np.ndarray
        Eclipse fraction in [0, 1] at each timestep.
    """
    toi_deg = np.asarray(toi_deg, dtype=float)
    sun_ang_diam_arcsec = np.asarray(sun_ang_diam_arcsec, dtype=float)
    ib_ang_diam_arcsec = np.asarray(ib_ang_diam_arcsec, dtype=float)

    # Angular separation in arcseconds (use absolute value of T-O-I)
    separation = np.abs(toi_deg) * 3600.0

    # Angular radii
    r_sun = sun_ang_diam_arcsec / 2.0
    r_ib = ib_ang_diam_arcsec / 2.0

    # Overlap area of Sun and IB disks
    overlap = circle_overlap_area(r_sun, r_ib, separation)

    # Fraction of Sun disk obscured
    sun_area = np.pi * r_sun ** 2
    # Guard against zero sun area (shouldn't happen, but be safe)
    frac = np.where(sun_area > 0, overlap / sun_area, 0.0)

    return np.clip(frac, 0.0, 1.0)
