"""Extract UKB-PPP Olink Explore inflammation-panel proteins.

The script downloads the CRAN source package `pQTLdata`, extracts the bundled
Olink Explore 3072 assay list, parses the XLSX without external dependencies,
and writes the full panel metadata plus the Inflammation/Inflammation_II subset.
"""

from __future__ import annotations

import csv
import io
import tarfile
import urllib.request
import zipfile
from pathlib import Path
from xml.etree import ElementTree as ET


PROJECT_ROOT = Path(__file__).resolve().parents[1]
CRAN_URL = "https://cran.r-project.org/src/contrib/pQTLdata_0.6.tar.gz"
XLSX_PATH = "pQTLdata/inst/Olink/assay-list-olink-explore-3072.xlsx"

FULL_OUT = PROJECT_ROOT / "data" / "metadata" / "olink_explore_3072_panel_metadata.csv"
PANEL_OUT = PROJECT_ROOT / "data" / "exposure" / "ukbppp_olink_inflammation_panel_proteins.csv"
SUMMARY_OUT = PROJECT_ROOT / "data" / "exposure" / "ukbppp_olink_inflammation_panel_summary.md"


def download_xlsx() -> bytes:
    with urllib.request.urlopen(CRAN_URL, timeout=60) as response:
        package_bytes = response.read()

    with tarfile.open(fileobj=io.BytesIO(package_bytes), mode="r:gz") as tar:
        member = tar.extractfile(XLSX_PATH)
        if member is None:
            raise FileNotFoundError(f"{XLSX_PATH} not found in {CRAN_URL}")
        return member.read()


def shared_strings(zfile: zipfile.ZipFile) -> list[str]:
    ns = {"m": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
    root = ET.fromstring(zfile.read("xl/sharedStrings.xml"))
    strings: list[str] = []
    for si in root.findall("m:si", ns):
        strings.append("".join((t.text or "") for t in si.findall(".//m:t", ns)))
    return strings


def column_index(cell_ref: str) -> int:
    letters = "".join(ch for ch in cell_ref if ch.isalpha())
    index = 0
    for letter in letters:
        index = index * 26 + ord(letter.upper()) - 64
    return index - 1


def parse_xlsx(xlsx_bytes: bytes) -> list[dict[str, str]]:
    ns = {"m": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
    with zipfile.ZipFile(io.BytesIO(xlsx_bytes)) as zfile:
        strings = shared_strings(zfile)
        root = ET.fromstring(zfile.read("xl/worksheets/sheet1.xml"))

    rows: dict[int, list[str]] = {}
    for row in root.findall(".//m:row", ns):
        row_number = int(row.attrib["r"])
        values = [""] * 4
        for cell in row.findall("m:c", ns):
            col = column_index(cell.attrib.get("r", "A1"))
            if col >= 4:
                continue
            node = cell.find("m:v", ns)
            value = node.text if node is not None else ""
            if cell.attrib.get("t") == "s" and value:
                value = strings[int(value)]
            values[col] = value.strip()
        rows[row_number] = values

    header = rows[3]
    records: list[dict[str, str]] = []
    for row_number in sorted(rows):
        if row_number < 4:
            continue
        values = rows[row_number]
        if not any(values):
            continue
        records.append(dict(zip(header, values)))
    return records


def write_csv(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    records = parse_xlsx(download_xlsx())
    for record in records:
        record["source"] = "pQTLdata_0.6::Olink_Explore_3072"
        record["source_url"] = CRAN_URL

    fieldnames = ["UniProt ID", "Protein name", "Gene name", "Explore 384 panel", "source", "source_url"]
    write_csv(FULL_OUT, records, fieldnames)

    inflammation = [
        row
        for row in records
        if row["Explore 384 panel"] in {"Inflammation", "Inflammation_II"}
    ]
    write_csv(PANEL_OUT, inflammation, fieldnames)

    counts: dict[str, int] = {}
    for row in inflammation:
        counts[row["Explore 384 panel"]] = counts.get(row["Explore 384 panel"], 0) + 1

    SUMMARY_OUT.parent.mkdir(parents=True, exist_ok=True)
    SUMMARY_OUT.write_text(
        "\n".join(
            [
                "# UKB-PPP Olink Explore Inflammation Panel Summary",
                "",
                f"Source: `{CRAN_URL}`",
                "",
                f"Full Olink Explore metadata rows: {len(records)}",
                f"Main exposure panel rows: {len(inflammation)}",
                "",
                "Panel counts:",
                *(f"- {panel}: {count}" for panel, count in sorted(counts.items())),
                "",
                f"Full metadata: `{FULL_OUT}`",
                f"Main exposure list: `{PANEL_OUT}`",
            ]
        ),
        encoding="utf-8",
    )

    print(f"Wrote {len(records)} full metadata rows to {FULL_OUT}")
    print(f"Wrote {len(inflammation)} inflammation-panel rows to {PANEL_OUT}")


if __name__ == "__main__":
    main()

