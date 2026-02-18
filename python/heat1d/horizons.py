"""JPL Horizons API integration for SPICE-based illumination.

Queries the JPL Horizons system for the Sun's apparent position as seen
from a surface point on any solar system body.  The resulting solar
elevation and distance are converted to absorbed surface flux using the
same angle-dependent albedo model used elsewhere in heat1d.

The flux array produced by this module plugs directly into the existing
``Model(flux_series=..., flux_dt=...)`` interface — no changes to the
core thermal model are required.

Uses only the Python standard library (``urllib``, ``json``).
"""

import json
import math
import urllib.parse
import urllib.request
from datetime import datetime

import numpy as np

from .properties import albedoVar

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

HORIZONS_API_URL = "https://ssd.jpl.nasa.gov/api/horizons.api"

#: Map from ``planets`` package names to Horizons body IDs.
HORIZONS_BODY_IDS = {
    "Moon": "301",
    "Mercury": "199",
    "Venus": "299",
    "Earth": "399",
    "Mars": "499",
    "Jupiter": "599",
    "Saturn": "699",
    "Uranus": "799",
    "Neptune": "899",
    "Pluto": "999",
    "Titan": "606",
    "Europa": "502",
    "Ganymede": "503",
    "Triton": "801",
    "Bennu": "2101955",
}


# ---------------------------------------------------------------------------
# Exceptions
# ---------------------------------------------------------------------------

class HorizonsError(Exception):
    """Raised when a Horizons query fails or the response cannot be parsed."""
    pass


# ---------------------------------------------------------------------------
# Body ID resolution
# ---------------------------------------------------------------------------

def get_body_id(planet_name, body_id_override=None):
    """Resolve a Horizons body ID from a planet name or explicit override.

    Parameters
    ----------
    planet_name : str
        Name from the ``planets`` package (e.g. ``"Moon"``).
    body_id_override : str or None
        If provided, returned directly (takes precedence).

    Returns
    -------
    str
        Horizons-compatible body ID string.

    Raises
    ------
    ValueError
        If *planet_name* is not in the built-in mapping and no override
        is given.
    """
    if body_id_override is not None:
        return str(body_id_override)
    if planet_name in HORIZONS_BODY_IDS:
        return HORIZONS_BODY_IDS[planet_name]
    raise ValueError(
        f"No Horizons body ID for '{planet_name}'. "
        f"Known bodies: {', '.join(sorted(HORIZONS_BODY_IDS))}. "
        f"Use --body-id to specify an explicit Horizons ID."
    )


# ---------------------------------------------------------------------------
# Step-size helpers
# ---------------------------------------------------------------------------

def compute_step_size(output_interval_s=None, planet_day_s=None,
                      default_steps=480):
    """Determine an appropriate Horizons step-size string.

    Parameters
    ----------
    output_interval_s : float or None
        Desired output interval [seconds].
    planet_day_s : float or None
        Diurnal period [seconds], used when *output_interval_s* is None.
    default_steps : int
        Fallback number of steps per diurnal period.

    Returns
    -------
    step_str : str
        Horizons-format step size (e.g. ``"10m"``, ``"1h"``).
    dt_seconds : float
        Actual time spacing [seconds].
    """
    if output_interval_s is not None:
        dt = output_interval_s
    elif planet_day_s is not None:
        dt = planet_day_s / default_steps
    else:
        dt = 600.0  # 10-minute default

    minutes = max(int(round(dt / 60.0)), 1)
    return f"{minutes}m", minutes * 60.0


# ---------------------------------------------------------------------------
# Horizons API query
# ---------------------------------------------------------------------------

def query_horizons(body_id, lon_deg, lat_deg, start_time, stop_time,
                   step_size="10m", timeout=60, quantities="4,20",
                   command="10"):
    """Query JPL Horizons for a target's position from a body surface point.

    Parameters
    ----------
    body_id : str
        Horizons body ID for the observer body (e.g. ``"301"`` for Moon).
    lon_deg : float
        East longitude [degrees].
    lat_deg : float
        Latitude [degrees].
    start_time : str
        UTC start, e.g. ``"2024-06-01 12:00"``.
    stop_time : str
        UTC stop, e.g. ``"2024-07-01 00:00"``.
    step_size : str
        Horizons step size (e.g. ``"10m"``, ``"1h"``).
    timeout : int
        HTTP timeout [seconds].
    quantities : str
        Horizons QUANTITIES string.  ``"4,20"`` = AZ/EL + distance;
        ``"4,13,20,25"`` = also angular diameter and eclipse geometry;
        ``"13"`` = angular diameter only (for parent body queries).
    command : str
        Horizons COMMAND (target body).  ``"10"`` = Sun (default),
        or a body ID like ``"399"`` for Earth.

    Returns
    -------
    dict
        Parsed ephemeris data.  Keys depend on *quantities* requested.

    Raises
    ------
    HorizonsError
        If the API returns an error or the response cannot be parsed.
    """
    params = {
        "format": "json",
        "COMMAND": f"'{command}'",
        "OBJ_DATA": "'NO'",
        "MAKE_EPHEM": "'YES'",
        "EPHEM_TYPE": "'OBSERVER'",
        "CENTER": f"'coord@{body_id}'",
        "COORD_TYPE": "'GEODETIC'",
        "SITE_COORD": f"'{lon_deg},{lat_deg},0'",
        "START_TIME": f"'{start_time}'",
        "STOP_TIME": f"'{stop_time}'",
        "STEP_SIZE": f"'{step_size}'",
        "QUANTITIES": f"'{quantities}'",
        "CSV_FORMAT": "'YES'",
        "SUPPRESS_RANGE_RATE": "'YES'",
        "APPARENT": "'AIRLESS'",
    }

    url = HORIZONS_API_URL + "?" + urllib.parse.urlencode(params)

    try:
        req = urllib.request.Request(url)
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            raw = resp.read().decode("utf-8")
    except urllib.error.URLError as exc:
        raise HorizonsError(f"Cannot reach Horizons API: {exc}") from exc
    except Exception as exc:
        raise HorizonsError(f"HTTP request failed: {exc}") from exc

    try:
        data = json.loads(raw)
    except json.JSONDecodeError as exc:
        raise HorizonsError(f"Invalid JSON from Horizons: {exc}") from exc

    result_text = data.get("result", "")
    if not result_text:
        raise HorizonsError("Empty result from Horizons API")

    return _parse_horizons_response(result_text, quantities=quantities)


def _parse_horizons_response(result_text, quantities="4,20"):
    """Parse the text inside a Horizons JSON ``result`` field.

    Expects ``CSV_FORMAT='YES'`` and ``SUPPRESS_RANGE_RATE='YES'``.
    The column layout depends on the *quantities* requested.

    Supported layouts:

    - ``"4,20"``: Date, f1, f2, Azi(3), Elev(4), delta(5)
    - ``"4,13,20,25"``: Date, f1, f2, Azi(3), Elev(4), Ang-diam(5),
      delta(6), T-O-I(7), IB_Illu%(8)
    - ``"13"``: Date, f1, f2, Ang-diam(3)

    Parameters
    ----------
    result_text : str
        The ``result`` value from the Horizons JSON response.
    quantities : str
        Which QUANTITIES were requested (determines column layout).

    Returns
    -------
    dict
        Parsed ephemeris.  Always includes ``times_utc``, ``dt_seconds``,
        ``n_samples``.  Other keys depend on *quantities*.

    Raises
    ------
    HorizonsError
        On parse failure or if ``$$SOE``/``$$EOE`` markers are missing.
    """
    # Column index maps for supported quantity combinations
    # Columns 0-2 are always: Date, flag1, flag2
    COLUMN_MAPS = {
        "4,20": {
            "min_cols": 6,
            "azimuth": 3, "elevation": 4, "delta": 5,
        },
        "4,13,20,25": {
            "min_cols": 9,
            "azimuth": 3, "elevation": 4, "sun_ang_diam": 5,
            "delta": 6, "toi": 7, "ib_illu_pct": 8,
        },
        "13": {
            "min_cols": 4,
            "ang_diam": 3,
        },
    }

    col_map = COLUMN_MAPS.get(quantities)
    if col_map is None:
        raise HorizonsError(
            f"Unsupported quantities layout: '{quantities}'. "
            f"Supported: {', '.join(sorted(COLUMN_MAPS))}"
        )

    lines = result_text.split("\n")

    # Check for error messages (no $$SOE marker)
    soe_idx = None
    eoe_idx = None
    for i, line in enumerate(lines):
        if "$$SOE" in line:
            soe_idx = i
        if "$$EOE" in line:
            eoe_idx = i

    if soe_idx is None or eoe_idx is None:
        msg = result_text.strip()
        if len(msg) > 300:
            msg = msg[:300] + "..."
        raise HorizonsError(f"Horizons returned an error: {msg}")

    # Parse data rows between $$SOE and $$EOE
    times_utc = []
    columns = {k: [] for k in col_map if k != "min_cols"}
    min_cols = col_map["min_cols"]

    for i in range(soe_idx + 1, eoe_idx):
        line = lines[i].strip()
        if not line:
            continue

        parts = [p.strip() for p in line.split(",")]
        if len(parts) < min_cols:
            continue

        times_utc.append(parts[0])
        for key, idx in col_map.items():
            if key == "min_cols":
                continue
            columns[key].append(float(parts[idx]))

    if not times_utc:
        raise HorizonsError("No data rows found between $$SOE and $$EOE")

    n = len(times_utc)

    # Compute dt from timestamps
    if n >= 2:
        t0 = _parse_horizons_date(times_utc[0])
        t1 = _parse_horizons_date(times_utc[1])
        dt_seconds = (t1 - t0).total_seconds()
    else:
        dt_seconds = 0.0

    # Build result dict
    result = {
        "times_utc": times_utc,
        "dt_seconds": dt_seconds,
        "n_samples": n,
    }

    # Map internal column names to output key names
    key_rename = {
        "azimuth": "azimuth_deg",
        "elevation": "solar_elevation_deg",
        "delta": "observer_range_au",
        "sun_ang_diam": "sun_ang_diam_arcsec",
        "toi": "toi_deg",
        "ib_illu_pct": "ib_illu_pct",
        "ang_diam": "ang_diam_arcsec",
    }

    for key, values in columns.items():
        out_key = key_rename.get(key, key)
        result[out_key] = np.array(values)

    return result


def _parse_horizons_date(date_str):
    """Parse a Horizons date string like ``'2024-Jun-01 12:00'``.

    Returns
    -------
    datetime
    """
    return datetime.strptime(date_str.strip(), "%Y-%b-%d %H:%M")


# ---------------------------------------------------------------------------
# Flux computation
# ---------------------------------------------------------------------------

def horizons_to_flux(solar_elevation_deg, observer_range_au, planet):
    """Convert Horizons ephemeris to absorbed surface flux.

    Applies the angle-dependent albedo model from Keihm (1984) /
    Vasavada et al. (2012) / Hayne et al. (2017, Eq. A8), matching the
    physics in :meth:`Model.surfFlux` and
    :func:`generate_flux.compute_flux_array`.

    Parameters
    ----------
    solar_elevation_deg : np.ndarray
        Solar elevation angle [degrees].  Negative = below horizon.
    observer_range_au : np.ndarray
        Sun–observer distance [AU].
    planet : object
        Planet object (needs ``S``, ``albedo``, ``albedoCoef``, ``rAU``).

    Returns
    -------
    flux : np.ndarray
        Absorbed surface flux [W/m²].
    """
    elev = np.asarray(solar_elevation_deg, dtype=float)
    r = np.asarray(observer_range_au, dtype=float)

    # Solar zenith angle (= incidence angle for horizontal surface)
    zenith_rad = np.deg2rad(90.0 - elev)
    cos_z = np.cos(zenith_rad)

    # Sun below horizon → zero flux
    night = elev <= 0.0
    cos_z[night] = 0.0

    # Incidence angle for albedo model (clip to [0, π/2])
    inc = np.clip(zenith_rad, 0.0, np.pi / 2.0)

    # Angle-dependent albedo
    a_coef = planet.albedoCoef[0]
    b_coef = planet.albedoCoef[1]
    A_var = albedoVar(planet.albedo, a_coef, b_coef, inc)

    # Absorbed flux correction factor
    f = (1.0 - A_var) / (1.0 - planet.albedo)

    # Baseline absorbed flux
    Sabs = planet.S * (1.0 - planet.albedo)

    # Inverse-square distance scaling
    flux = f * Sabs * (planet.rAU / r) ** 2 * cos_z

    # Ensure no negative flux from numerical edge cases
    flux[night] = 0.0

    return flux


# ---------------------------------------------------------------------------
# Eclipse detection
# ---------------------------------------------------------------------------

def query_parent_body(parent_id, satellite_id, lon_deg, lat_deg,
                      start_time, stop_time, step_size="10m", timeout=60):
    """Query Horizons for the parent body's angular diameter.

    Uses the **same observer location** (``lon_deg``, ``lat_deg``) on the
    satellite surface as the main Sun query, so that the angular diameter
    is computed from the correct surface point.

    Parameters
    ----------
    parent_id : str
        Horizons body ID for the parent (e.g. ``"399"`` for Earth).
    satellite_id : str
        Horizons body ID for the satellite (e.g. ``"301"`` for Moon).
    lon_deg : float
        East longitude [degrees] on the satellite surface.
    lat_deg : float
        Latitude [degrees] on the satellite surface.
    start_time : str
        UTC start time.
    stop_time : str
        UTC stop time.
    step_size : str
        Horizons step size.
    timeout : int
        HTTP timeout [seconds].

    Returns
    -------
    dict
        Keys: ``ang_diam_arcsec`` (ndarray), ``times_utc`` (list),
        ``n_samples`` (int), ``dt_seconds`` (float).
    """
    return query_horizons(
        body_id=satellite_id,
        lon_deg=lon_deg,
        lat_deg=lat_deg,
        start_time=start_time,
        stop_time=stop_time,
        step_size=step_size,
        timeout=timeout,
        quantities="13",
        command=parent_id,
    )


def apply_horizons_eclipses(flux, sun_result, parent_result):
    """Apply eclipse flux reduction from Horizons geometry.

    Parameters
    ----------
    flux : np.ndarray
        Absorbed flux array [W/m²] (modified in-place).
    sun_result : dict
        Enhanced Sun query result with ``toi_deg`` and
        ``sun_ang_diam_arcsec``.
    parent_result : dict
        Parent body query result with ``ang_diam_arcsec``.

    Returns
    -------
    dict
        Eclipse metadata: ``n_eclipses``, ``total_eclipse_samples``,
        ``max_fraction``, ``eclipse_fraction`` (the full array).
    """
    from .eclipse import compute_eclipse_fraction

    toi_deg = sun_result["toi_deg"]
    sun_diam = sun_result["sun_ang_diam_arcsec"]
    ib_diam = parent_result["ang_diam_arcsec"]

    # Align array lengths (truncate to shorter if Horizons returned
    # different sample counts due to rounding)
    n = min(len(flux), len(toi_deg), len(ib_diam))
    toi_deg = toi_deg[:n]
    sun_diam = sun_diam[:n]
    ib_diam = ib_diam[:n]

    frac = compute_eclipse_fraction(toi_deg, sun_diam, ib_diam)

    # Apply reduction
    flux[:n] *= (1.0 - frac)

    # Count distinct eclipse events (contiguous runs of frac > 0)
    eclipsed = frac > 0.0
    n_eclipses = 0
    in_eclipse = False
    for e in eclipsed:
        if e and not in_eclipse:
            n_eclipses += 1
            in_eclipse = True
        elif not e:
            in_eclipse = False

    return {
        "n_eclipses": n_eclipses,
        "total_eclipse_samples": int(np.sum(eclipsed)),
        "max_fraction": float(frac.max()) if len(frac) > 0 else 0.0,
        "eclipse_fraction": frac,
    }


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

def fetch_solar_flux(planet_name, lon_deg, lat_deg, start_time, stop_time,
                     body_id=None, output_interval_s=None, planet_day_s=None,
                     planet=None, timeout=60, eclipses=True,
                     parent_body_id=None):
    """Query Horizons and return an absorbed-flux array ready for the Model.

    This is the primary entry point used by the CLI.

    When *eclipses* is ``True`` and the body is a natural satellite,
    automatically queries for eclipse geometry and reduces flux during
    eclipse events.  Both the Sun query and the parent-body query use
    the same observer latitude/longitude, so the eclipse fraction
    correctly depends on the surface position.

    Parameters
    ----------
    planet_name : str
        Planet name (e.g. ``"Moon"``).
    lon_deg : float
        East longitude [degrees].
    lat_deg : float
        Latitude [degrees].
    start_time : str
        UTC start time.
    stop_time : str
        UTC stop time.
    body_id : str or None
        Explicit Horizons body ID override.
    output_interval_s : float or None
        Desired output interval [seconds] (controls query step size).
    planet_day_s : float or None
        Diurnal period [seconds] (used to compute default step size).
    planet : object or None
        Planet object for albedo model.  If ``None``, loaded from
        ``planets`` package using *planet_name*.
    timeout : int
        HTTP timeout [seconds].
    eclipses : bool
        If ``True`` (default), detect and apply eclipses for satellites.
    parent_body_id : str or None
        Explicit parent body ID override for eclipse detection.

    Returns
    -------
    flux : np.ndarray
        Absorbed surface flux [W/m²].
    dt : float
        Time spacing [seconds].
    metadata : dict
        Query metadata including optional ``eclipse_info``.

    Raises
    ------
    HorizonsError
        If the query fails.
    """
    from .eclipse import get_parent_body_id as _get_parent, is_satellite

    # Resolve body ID
    bid = get_body_id(planet_name, body_id_override=body_id)

    # Resolve planet object
    if planet is None:
        import planets as planets_pkg
        planet = getattr(planets_pkg, planet_name, None)
        if planet is None:
            raise HorizonsError(
                f"Cannot load planet '{planet_name}' from planets package. "
                f"Provide a planet object via the 'planet' parameter."
            )

    # Compute step size
    step_str, dt = compute_step_size(
        output_interval_s=output_interval_s,
        planet_day_s=planet_day_s or getattr(planet, "day", None),
    )

    # Determine if eclipse detection applies
    detect_eclipses = eclipses and (
        is_satellite(bid) or parent_body_id is not None
    )
    if detect_eclipses:
        quantities = "4,13,20,25"
        pid = parent_body_id or _get_parent(bid)
    else:
        quantities = "4,20"

    # Query Horizons for Sun position
    result = query_horizons(
        body_id=bid,
        lon_deg=lon_deg,
        lat_deg=lat_deg,
        start_time=start_time,
        stop_time=stop_time,
        step_size=step_str,
        timeout=timeout,
        quantities=quantities,
    )

    # Use the actual dt from the response (in case of rounding)
    if result["dt_seconds"] > 0:
        dt = result["dt_seconds"]

    # Convert to absorbed flux
    flux = horizons_to_flux(
        result["solar_elevation_deg"],
        result["observer_range_au"],
        planet,
    )

    # Apply eclipse reductions if applicable
    eclipse_info = None
    if detect_eclipses:
        parent_result = query_parent_body(
            parent_id=pid,
            satellite_id=bid,
            lon_deg=lon_deg,
            lat_deg=lat_deg,
            start_time=start_time,
            stop_time=stop_time,
            step_size=step_str,
            timeout=timeout,
        )
        eclipse_info = apply_horizons_eclipses(flux, result, parent_result)

    metadata = {
        "body_id": bid,
        "n_samples": result["n_samples"],
        "start_time": start_time,
        "stop_time": stop_time,
        "step_size": step_str,
        "times_utc": result["times_utc"],
        "eclipse_info": eclipse_info,
    }

    return flux, dt, metadata
