import json
from pathlib import Path


def read_manifest(run_dir: Path) -> dict:
    """Load and JSON-schema-validate manifest.json from a run dir."""
    run_dir = Path(run_dir)
    with open(run_dir / "manifest.json") as f:
        data = json.load(f)
    try:
        import jsonschema
        schema_path = Path(__file__).parent.parent.parent / "schemas" / "manifest.v1.json"
        if schema_path.exists():
            with open(schema_path) as sf:
                schema = json.load(sf)
            jsonschema.Draft7Validator(schema).validate(data)
    except ImportError:  # pragma: no cover
        pass
    return data
