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
    list of (str, str)
        List of ``(body_id, display_name)`` tuples.  Empty list if
        no match found.
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
    if "Multiple major-bodies match" in result_text \
       or "Multiple small-bodies match" in result_text \
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
                        int(parts[0])
                        results.append((parts[0], parts[1].strip()))
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
                return [(body_id, name_part)]

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
                return [(body_id, name_part)]

    # Case: single asteroid/small-body match
    # Header line: "JPL/HORIZONS               269 Justitia (A887 SA)     2026-..."
    # Or "Rec #:     269"
    for line in lines:
        if "Rec #:" in line:
            m = re.search(r'Rec #:\s*(-?\d+)', line)
            if m:
                body_id = m.group(1)
                # Extract name from JPL/HORIZONS header line
                name_part = body_id
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
                return [(body_id, name_part)]

    # Case: no useful match found
    if "No matches found" in result_text or "Cannot find" in result_text:
        return []

    return []
