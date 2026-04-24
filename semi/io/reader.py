import json
import warnings
from pathlib import Path

# Major version of the manifest schema (schemas/manifest.v1.json) this engine
# build accepts when reading run artifacts back. Symmetric with
# `semi.schema.ENGINE_SUPPORTED_SCHEMA_MAJOR` on the input side.
ENGINE_SUPPORTED_MANIFEST_MAJOR = 1


def read_manifest(run_dir: Path) -> dict:
    """Load and JSON-schema-validate manifest.json from a run dir.

    Also enforces a major-version gate on the manifest's ``schema_version``
    field: mismatched majors raise ``ValueError`` so a future engine cannot
    silently misinterpret an old artifact tree. Missing ``schema_version``
    is tolerated (with a warning) to keep pre-M11 manifests readable.
    """
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

    if "schema_version" in data:
        try:
            major = int(str(data["schema_version"]).split(".")[0])
        except (ValueError, AttributeError) as exc:
            raise ValueError(
                f"manifest schema_version {data['schema_version']!r} is not a valid "
                f"semver string"
            ) from exc
        if major != ENGINE_SUPPORTED_MANIFEST_MAJOR:
            raise ValueError(
                f"manifest schema_version {data['schema_version']!r} has major "
                f"{major}; engine supports major {ENGINE_SUPPORTED_MANIFEST_MAJOR}"
            )
    else:
        warnings.warn(
            f"manifest at {run_dir} has no schema_version field; assuming compatible",
            stacklevel=2,
        )
    return data
