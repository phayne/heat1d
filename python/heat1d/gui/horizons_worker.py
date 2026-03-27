"""Background worker for Horizons body name search."""

import json
import re
import urllib.parse
import urllib.request

from PySide6.QtCore import QThread, Signal

HORIZONS_API_URL = "https://ssd.jpl.nasa.gov/api/horizons.api"


class HorizonsSearchWorker(QThread):
    """Search for a body by name via the Horizons API.

    Emits ``found`` with a list of ``(body_id, display_name)`` tuples,
    or ``error`` with a message string.
    """

    found = Signal(list)    # list of (body_id_str, display_name)
    error = Signal(str)

    def __init__(self, name, timeout=15, parent=None):
        super().__init__(parent)
        self.name = name
        self.timeout = timeout

    def run(self):
        try:
            results = search_horizons_body(self.name, timeout=self.timeout)
            self.found.emit(results)
        except Exception as exc:
            self.error.emit(str(exc))


def search_horizons_body(name, timeout=15):
    """Query JPL Horizons for body name resolution.

    Parameters
    ----------
    name : str
        Body name to search (e.g. "Enceladus", "Moon", "Bennu").
    timeout : int
        HTTP timeout [seconds].

    Returns
    -------
    list of (str, str, dict)
        List of ``(body_id, display_name, obj_data)`` tuples.
        ``obj_data`` is a dict of physical/orbital parameters parsed
        from the Horizons OBJ_DATA response (empty for multi-match
        disambiguation results).  Empty list if no match found.
    """
    params = {
        "format": "json",
        "COMMAND": f"'{name}'",
        "OBJ_DATA": "'YES'",
        "MAKE_EPHEM": "'NO'",
    }

    url = HORIZONS_API_URL + "?" + urllib.parse.urlencode(params)

    req = urllib.request.Request(url)
    with urllib.request.urlopen(req, timeout=timeout) as resp:
        raw = resp.read().decode("utf-8")

    data = json.loads(raw)
    result_text = data.get("result", "")

    if not result_text:
        return []

    return _parse_search_result(result_text)


def _parse_search_result(result_text):
    """Parse Horizons search response into (body_id, name) pairs.

    Handles three cases:
    1. Unique match — extract body ID from the header
    2. Multiple matches — parse the disambiguation table
    3. No match — return empty list
    """
    lines = result_text.split("\n")

    # Case: multiple matches — look for disambiguation table
    # Handles both "Multiple major-bodies match" and "Multiple small-bodies match"
    is_small_body_list = "Multiple small-bodies match" in result_text
    if is_small_body_list \
       or "Multiple major-bodies match" in result_text \
       or "Multiple matches" in result_text:
        results = []
        in_table = False
        for line in lines:
            stripped = line.strip()
            # Table starts after a line with "ID#" / "Record #" or dashes
            if any(h in stripped for h in ("ID#", "Record #")):
                in_table = True
                continue
            if stripped.startswith("---") and not in_table:
                in_table = True
                continue
            if in_table:
                if not stripped or stripped.startswith("---"):
                    if results:
                        break
                    continue
                # Parse "  301  Moon  " or "  269  Justitia  A887 SA"
                parts = stripped.split(None, 1)
                if len(parts) >= 2:
                    try:
                        raw_id = int(parts[0])
                        # Small body catalog numbers → NAIF IDs to avoid
                        # collision with satellite ID range (100-999)
                        if is_small_body_list:
                            body_id = str(2000000 + raw_id)
                            display = parts[1].strip() + " [small body]"
                        else:
                            body_id = parts[0]
                            display = parts[1].strip()
                        results.append((body_id, display, {}))
                    except ValueError:
                        pass
        return results

    # Case: unique match — look for body ID in the header
    for line in lines:
        # Format 1: "Target body name: Moon (301)    {source: ...}"
        if "Target body name:" in line:
            after = line.split("Target body name:")[1]
            if "(" in after and ")" in after:
                start = after.index("(") + 1
                end = after.index(")")
                body_id = after[start:end].strip()
                name_part = after[:after.index("(")].strip()
                return [(body_id, name_part, _parse_obj_data(result_text))]

        # Format 2: " Revised: ...   Enceladus / (Saturn)   602"
        # The body ID is the last token on the "Revised:" line
        if "Revised:" in line:
            # Look for name + ID pattern like "Enceladus / (Saturn)   602"
            # or "Moon                    301"
            parts = line.split()
            if parts and parts[-1].isdigit():
                body_id = parts[-1]
                # Extract the name: between "Revised: ..." and the ID
                after_revised = line.split("Revised:")[1]
                # Remove date portion (e.g., "Jan 26, 2022")
                name_match = re.search(
                    r'\d{4}\s+(.*?)\s+' + re.escape(body_id), after_revised
                )
                if name_match:
                    name_part = name_match.group(1).strip().rstrip("/").strip()
                else:
                    name_part = body_id
                return [(body_id, name_part, _parse_obj_data(result_text))]

    # Case: single asteroid/small-body match
    # Header line: "JPL/HORIZONS               269 Justitia (A887 SA)     2026-..."
    # Or "Rec #:     269"
    for line in lines:
        if "Rec #:" in line:
            m = re.search(r'Rec #:\s*(-?\d+)', line)
            if m:
                catalog_num = int(m.group(1))
                # Convert catalog number to NAIF ID to avoid collision
                # with satellite ID range (100-999)
                body_id = str(2000000 + catalog_num)
                # Extract name from JPL/HORIZONS header line
                name_part = str(catalog_num)
                for hline in lines:
                    if "JPL/HORIZONS" in hline:
                        # "JPL/HORIZONS  269 Justitia (A887 SA)  2026-Feb-16 ..."
                        hm = re.search(
                            r'JPL/HORIZONS\s+.*?(\d+)\s+(.+?)\s+\d{4}-',
                            hline,
                        )
                        if hm:
                            name_part = hm.group(2).strip()
                        break
                # Tag as small body — SPICE surface ephemeris not available
                name_part += " [small body]"
                return [(body_id, name_part, _parse_obj_data(result_text))]

    # Case: no useful match found
    if "No matches found" in result_text or "Cannot find" in result_text:
        return []

    return []


def _parse_obj_data(result_text):
    """Extract physical and orbital parameters from Horizons OBJ_DATA.

    Parses the ``Asteroid physical parameters`` and ``IAU76/J2000``
    orbital elements sections that Horizons returns when
    ``OBJ_DATA='YES'``.

    Parameters
    ----------
    result_text : str
        Full Horizons result text.

    Returns
    -------
    dict
        Keys: ``radius_km``, ``rotation_period_h``, ``albedo``,
        ``semi_major_axis_au``, ``eccentricity``, ``orbital_period_yr``,
        ``spectral_type``.  Values are ``None`` when not available.
    """
    data = {
        "radius_km": None,
        "rotation_period_h": None,
        "albedo": None,
        "semi_major_axis_au": None,
        "eccentricity": None,
        "orbital_period_yr": None,
        "spectral_type": None,
    }

    # --- Physical parameters ---
    # Format: "RAD= 25.364", "ROTPER= 33.128", "ALBEDO= .061"
    m = re.search(r'RAD=\s*([\d.]+)', result_text)
    if m:
        data["radius_km"] = float(m.group(1))

    m = re.search(r'ROTPER=\s*([\d.]+)', result_text)
    if m:
        data["rotation_period_h"] = float(m.group(1))
    # If "ROTPER= n.a." → no match, stays None

    m = re.search(r'ALBEDO=\s*([.\d]+)', result_text)
    if m:
        data["albedo"] = float(m.group(1))

    m = re.search(r'STYP=\s*(\S+)', result_text)
    if m and m.group(1).lower() != "n.a.":
        data["spectral_type"] = m.group(1)

    # --- Orbital elements ---
    # These appear in a block after "IAU76/J2000 helio. ecliptic osc. elements"
    # Format: "   EC= .2131484921430992", "   A= 2.616505270376438"
    #         "   PER= 4.23244"
    m = re.search(r'\bEC=\s*([\d.]+)', result_text)
    if m:
        data["eccentricity"] = float(m.group(1))

    # Semi-major axis: match "A= 2.616..." but not "MA=" or "ADIST="
    m = re.search(r'(?<![A-Z])\bA=\s*([\d.]+)', result_text)
    if m:
        data["semi_major_axis_au"] = float(m.group(1))

    m = re.search(r'\bPER=\s*([\d.]+)', result_text)
    if m:
        data["orbital_period_yr"] = float(m.group(1))

    return data
