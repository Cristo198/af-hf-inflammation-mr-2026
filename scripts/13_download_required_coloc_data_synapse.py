"""Download required UKB-PPP files for formal FGF5/LPA colocalization.

This script requires a Synapse Personal Access Token with View and Download
permissions. Set it before running:

PowerShell:
    $env:SYNAPSE_AUTH_TOKEN="your-token"

Then run:
    python scripts/13_download_required_coloc_data_synapse.py
"""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[1]
OUT_DIR = PROJECT_ROOT / "data" / "raw" / "ukbppp_pgwas"

FILES = {
    "syn51469935": "FGF5_P12034_OID20490_v1_Inflammation.tar",
    "syn52361776": "LPA_P08519_OID30747_v1_Inflammation_II.tar",
}
MAP_FILES = {
    "syn51397015": "olink_rsid_map_mac5_info03_b0_7_chr4_patched_v2.tsv.gz",
    "syn51397018": "olink_rsid_map_mac5_info03_b0_7_chr6_patched_v2.tsv.gz",
}


def run(command: list[str]) -> None:
    print(" ".join(command))
    subprocess.run(command, check=True)


def main() -> None:
    token = os.environ.get("SYNAPSE_AUTH_TOKEN")
    if not token:
        raise SystemExit("Missing SYNAPSE_AUTH_TOKEN. Create a Synapse token with View + Download and set it first.")

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    map_dir = PROJECT_ROOT / "data" / "raw" / "ukbppp_metadata"
    map_dir.mkdir(parents=True, exist_ok=True)

    try:
        import synapseclient  # noqa: F401
    except ImportError:
        run([sys.executable, "-m", "pip", "install", "synapseclient"])

    run(["synapse", "login", "-p", token])
    for syn_id, filename in FILES.items():
        if (OUT_DIR / filename).exists():
            print(f"Already present: {OUT_DIR / filename}")
            continue
        run(["synapse", "get", syn_id, "--downloadLocation", str(OUT_DIR)])

    for syn_id, filename in MAP_FILES.items():
        if (map_dir / filename).exists():
            print(f"Already present: {map_dir / filename}")
            continue
        run(["synapse", "get", syn_id, "--downloadLocation", str(map_dir)])

    print(f"Downloaded required files to {OUT_DIR}")
    print(f"Downloaded required rsID maps to {map_dir}")


if __name__ == "__main__":
    main()
