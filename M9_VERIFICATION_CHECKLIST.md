# M9 Verification Checklist

The Claude Code agent completed M9 structurally but explicitly deferred
gates 3, 4, and 5 because they require dolfinx. Run this checklist
inside Docker before merging M9 or starting M10.

```bash
# From the repo root, on the M9 branch:
docker compose build # if image isn't current
```

## Gate 3: FEM artifact tests on all five benchmarks

```bash
docker compose run --rm test pytest tests/test_artifact.py -v
```

**Expected:** 3 pure-Python tests + 5 parametrized FEM tests, total 8
passed, 0 skipped, 0 failed. If any test skips inside Docker, the
`pytest.mark.skipif` condition is wrong and must be tightened to only
skip when dolfinx truly isn't importable.

## Gates 4+5 and sanity checks S1, S3, S4, S5 (single container)

Gates 4 and 5 and the sanity checks that read from `/tmp/runs` are
combined into one `docker compose run` invocation so that the data
generated in Gate 4 is still present when Gates 5 and S1/S3/S4/S5
validate it. Splitting them across separate `docker compose run` calls
would silently discard `/tmp/runs` when the Gate-4 container exits.

```bash
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
```

**Gate 4 expected:** Five subdirectories under `/tmp/runs/`, each
containing `manifest.json`, `input.json`, `mesh/`, `fields/`, `iv/`,
`convergence/`, `logs/`. Exit code 0 for every command.

**Gate 5 expected:** Five `OK` lines, zero `FAIL`, exit 0.

**S1 expected:** Two `OK` lines. If this fails, M10's subprocess design
is broken and the agent's claim about `reader.py` being stdlib-only is
wrong.

**S3 expected:** `input_sha256: OK`. Common bug is hashing the
parsed-and-re-serialized dict instead of the raw bytes of the original
file.

**S4 expected:** `MOS carrier fields present: OK`. If this fails the MOS
submesh case was not handled — `phi_n`/`phi_p` live on the silicon
submesh, not the parent mesh.

**S5 expected:** `bipolar IV: OK`. Confirms the bipolar sweep was
serialized correctly.

## S2: Determinism (self-contained, separate invocation)

S2 generates its own data, so it does not depend on `/tmp/runs` and can
run in its own container.

```bash
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
```

**Expected:** `determinism: OK`. If not, the manifest contains
non-deterministic content (timestamps, unsorted iteration order,
un-rounded floats) and downstream diff-based tooling will break.

---

## Reporting back

If all gates pass, M9 is truly done — merge and start M10.

If any gate fails, capture the failure output and feed it back to the
agent as the first instruction of a follow-up session. Do not paper
over failures by relaxing the acceptance tests; they were set
deliberately.
