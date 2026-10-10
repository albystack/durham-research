"""Convert the verified Hamilton compact archive into the public data bundle.

Usage: python analysis/export_data.py path/to/results-compact.tar.gz
Only measurement tables and their original verification metadata are read.
The input archive is retained locally; full tiling arrays stay on Hamilton.
"""

import argparse
import gzip
import hashlib
import io
import json
import tarfile
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
ARCHIVE_SHA256 = "90e30c0f97e0e22132110dfb35a0a6e5e7007d5bbe4c8270be3a10c9977c0254"


def digest(blob):
    return hashlib.sha256(blob).hexdigest()


def export(archive_path):
    """Verify original table hashes and preserve batch-prefixed pair identities."""
    assert digest(archive_path.read_bytes()) == ARCHIVE_SHA256
    destination = ROOT / "data"
    destination.mkdir(exist_ok=True)
    frames, certificates = [], {}
    with tarfile.open(archive_path) as archive:
        names = {m.name for m in archive.getmembers() if m.isfile()}

        def read(name):
            return archive.extractfile(name).read()

        for name in sorted(names):
            if name.endswith("/complete.json"):
                receipt = json.loads(read(name))
                batch_path = str(Path(name).parent)
                attempt = batch_path + "/" + receipt["attempt"]
                for filename, expected in receipt["files"].items():
                    assert digest(read(attempt + "/" + filename)) == expected
                table = read(attempt + "/cells/cell-00/measurements.csv")
                frame = pd.read_csv(
                    io.BytesIO(table),
                    dtype={
                        key: str
                        for key in (
                            "environment_seed",
                            "replica1_seed",
                            "replica2_seed",
                            "coin_seed",
                        )
                    },
                )
                frame["environment_id"] = (
                    Path(batch_path).name + ":" + frame.environment_id
                )
                frame["batch"] = receipt["batch"]["index"]
                frames.append(frame)
                certificates[Path(batch_path).name] = {
                    "receipt": receipt,
                    "raw_verification": json.loads(
                        read(attempt + "/raw_verification.json")
                    ),
                    "manifest": json.loads(
                        read(attempt + "/cells/cell-00/manifest.json")
                    ),
                }
        frame = pd.concat(frames, ignore_index=True)
        assert len(frame) == 516096 and frame.environment_id.nunique() == 258048
        assert len(certificates) == 504
        assert not frame.duplicated(["environment_id", "probe"]).any()
        # mtime=0 makes the compressed export reproducible byte for byte.
        payload = frame.to_csv(index=False, lineterminator="\n").encode()
        (destination / "measurements.csv.gz").write_bytes(
            gzip.compress(payload, mtime=0)
        )
        provenance = {
            "original_archive_sha256": ARCHIVE_SHA256,
            "launch": json.loads(read("launch.json")),
            "source_manifest": json.loads(read("code/source_manifest.json")),
            "scheduler_accounting": read("completed-accounting.txt").decode(),
            "batches": certificates,
        }
        (destination / "provenance.json.gz").write_bytes(
            gzip.compress(
                json.dumps(provenance, sort_keys=True, indent=2).encode(), mtime=0
            )
        )
        for name in ("analysis.json", "intervals.csv"):
            (destination / name).write_bytes(read("production/analysis/" + name))
    print(f"Exported {len(frame):,} rows from {len(certificates)} verified batches.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("archive", type=Path)
    export(parser.parse_args().archive)
