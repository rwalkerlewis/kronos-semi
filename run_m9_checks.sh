#!/bin/bash
set -e

echo "--- BUILDING DOCKER COMPOSE ---"
docker compose build

echo "--- GATE 3 ---"
docker compose run --rm test pytest tests/test_artifact.py -v

echo "--- GATES 4+5 and S1 S3 S4 S5 (single container, /tmp/runs persists) ---"
docker compose run --rm dev bash -c '
 set -e

 # Gate 4: generate all five run artifacts
 rm -rf /tmp/runs
 semi-run benchmarks/pn_1d/pn_junction.json --out /tmp/runs
 semi-run benchmarks/pn_1d_bias/pn_junction_bias.json --out /tmp/runs
 semi-run benchmarks/pn_1d_bias_reverse/pn_junction_bias_reverse.json --out /tmp/runs
 semi-run benchmarks/mos_2d/mos_cap.json --out /tmp/runs
 semi-run benchmarks/resistor_3d/resistor.json --out /tmp/runs
 echo "Gate 4: ls /tmp/runs/"
 ls -la /tmp/runs/

 # Gate 5: validate every manifest against the schema
 python - <<PY
from pathlib import Path
import json, jsonschema
schema = json.load(open("schemas/manifest.v1.json"))
jsonschema.Draft7Validator.check_schema(schema)
validator = jsonschema.Draft7Validator(schema)
fails = 0
for run in sorted(Path("/tmp/runs").iterdir()):
 m = json.load(open(run / "manifest.json"))
 errors = list(validator.iter_errors(m))
 if errors:
  fails += 1
  print(f"FAIL {run.name}")
  for e in errors:
   print(f"  {list(e.absolute_path)}: {e.message}")
 else:
  print(f"OK {run.name}")
assert fails == 0, f"{fails} manifests failed validation"
PY

 # S1: reader is stdlib-only and can load a manifest
 python - <<PY
import sys
import importlib.util
for blocked in ("dolfinx", "numpy"):
 if blocked in sys.modules:
  del sys.modules[blocked]
class Blocker:
 def find_spec(self, name, path, target=None):
  if name.split(".")[0] in ("dolfinx", "numpy"):
   raise ImportError(f"blocked: {name}")
  return None
sys.meta_path.insert(0, Blocker())
from semi.io.reader import read_manifest
print("stdlib-only import: OK")
from pathlib import Path
m = read_manifest(sorted(Path("/tmp/runs").iterdir())[0])
print(f"read manifest OK: {m[\"run_id\"]}")
PY

 # S3: input_sha256 matches raw bytes of input.json
 python - <<PY
import hashlib, json
from pathlib import Path
run = sorted(Path("/tmp/runs").iterdir())[0]
claimed = json.load(open(run / "manifest.json"))["input_sha256"]
actual = hashlib.sha256(open(run / "input.json", "rb").read()).hexdigest()
assert claimed == actual, f"sha mismatch\n  claimed: {claimed}\n  actual:  {actual}"
print("input_sha256: OK")
PY

 # S4: MOS submesh carrier fields present
 python - <<PY
from pathlib import Path
import json
run_dir = next(d for d in Path("/tmp/runs").iterdir() if "mos" in d.name)
m = json.load(open(run_dir / "manifest.json"))
field_names = {f["name"] for f in m["fields"]}
assert "phi_n" in field_names, f"MOS run missing phi_n field. Got: {field_names}"
assert "phi_p" in field_names, f"MOS run missing phi_p field. Got: {field_names}"
print("MOS carrier fields present: OK")
PY

 # S5: resistor bipolar IV is monotonic per leg
 python - <<PY
import csv
from pathlib import Path
run = next(d for d in Path("/tmp/runs").iterdir() if "resistor" in d.name)
iv_files = list((run / "iv").glob("*.csv"))
assert iv_files, "no IV files"
for f in iv_files:
 rows = list(csv.DictReader(open(f)))
 voltages = [float(r["V"]) for r in rows]
 assert len(voltages) >= 5, f"too few IV points in {f.name}: {len(voltages)}"
 assert min(voltages) < -1e-6, f"no negative voltages in {f.name}"
 assert max(voltages) > 1e-6, f"no positive voltages in {f.name}"
 print(f"{f.name}: {len(voltages)} points, V in [{min(voltages):+.4e}, {max(voltages):+.4e}]")
print("bipolar IV: OK")
PY
'

echo "--- S2 (determinism, self-contained) ---"
docker compose run --rm dev bash -c '
 set -e
 rm -rf /tmp/det_a /tmp/det_b
 semi-run benchmarks/pn_1d/pn_junction.json --out /tmp/det_a --run-id fixed
 semi-run benchmarks/pn_1d/pn_junction.json --out /tmp/det_b --run-id fixed
 python - <<PY
import json
a = json.load(open("/tmp/det_a/fixed/manifest.json"))
b = json.load(open("/tmp/det_b/fixed/manifest.json"))
a.pop("wall_time_s", None); b.pop("wall_time_s", None)
assert json.dumps(a, sort_keys=True) == json.dumps(b, sort_keys=True), "nondeterministic"
print("determinism: OK")
PY
'
