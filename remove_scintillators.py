#!/usr/bin/env python3
"""Remove only the two seven-sector trigger scintillator placements from develop."""
import hashlib
import json
import re
import shutil
import sys
from pathlib import Path


def remove_placements(text):
    pattern = r'      <loop for="j" from="0" to="6" step="1">\s*<physvol name="TRIG_plane([03])_\[j\+1\]">.*?</physvol>\s*</loop>'
    matches = list(re.finditer(pattern, text, re.S))
    if sorted(m.group(1) for m in matches) != ["0", "3"]:
        raise ValueError("Expected exactly the develop plane 0 and plane 3 scintillator placements")
    result = re.sub(pattern, "      <!-- Trigger scintillator placements removed for this campaign. -->", text, flags=re.S)
    if re.search(r'<physvol\s+name="TRIG_', result):
        raise ValueError("A trigger scintillator placement remains")
    # No GEM placement, position or material is rewritten.
    gem = r'<physvol name="GEM_plane.*?</physvol>'
    if re.findall(gem, text, re.S) != re.findall(gem, result, re.S):
        raise ValueError("GEM placements changed")
    return result


if __name__ == "__main__":
    path = Path("geometry/tracking/GEMTSDetectors.gdml")
    original = path.read_text()
    if sys.argv[1:2] == ["--verify"]:
        receipt = json.loads(Path("geometry-variant.json").read_text())
        if receipt["modified_sha256"] != hashlib.sha256(path.read_bytes()).hexdigest():
            raise SystemExit("Geometry hash differs from build receipt")
        if re.search(r'<physvol\s+name="TRIG_', original):
            raise SystemExit("Trigger scintillator placements present")
        mother = Path("geometry/mollerMother.gdml").read_text()
        if '<file name="tracking/GEMTSDetectors.gdml"/>' not in mother:
            raise SystemExit("Mother geometry does not load the verified tracking geometry")
        for name in ("geometry-variant.json", "geometry-variant.patch", "build-provenance.txt"):
            shutil.copyfile(name, Path(sys.argv[2]) / name)
        print("GEOMETRY_VERIFIED no trigger scintillator placements; develop GEMs unchanged")
        raise SystemExit(0)
    changed = remove_placements(original)
    path.write_text(changed)
    Path("geometry-variant.json").write_text(json.dumps(dict(
        variant="develop_no_trigger_scintillators", changed_file=str(path),
        removed_placements=14, removed_groups=["TRIG_plane0", "TRIG_plane3"],
        original_sha256=hashlib.sha256(original.encode()).hexdigest(),
        modified_sha256=hashlib.sha256(changed.encode()).hexdigest(),
        gem_placements_unchanged=True), indent=2) + "\n")
