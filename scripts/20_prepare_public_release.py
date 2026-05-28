from __future__ import annotations

import csv
import shutil
import zipfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
RELEASE_ROOT = ROOT / "public_release"
RELEASE_NAME = "cardio_no_lab_paper_project_public_v0_1_2026-05-28"
OUT = RELEASE_ROOT / RELEASE_NAME
ZIP_PATH = RELEASE_ROOT / f"{RELEASE_NAME}.zip"


ROOT_FILES = [
    ("README_PUBLIC.md", "README.md"),
    ("DATA_ACCESS.md", "DATA_ACCESS.md"),
    ("LICENSE", "LICENSE"),
    ("CITATION.cff", "CITATION.cff"),
    (".zenodo.json", ".zenodo.json"),
    (".gitignore", ".gitignore"),
    ("requirements.txt", "requirements.txt"),
    ("analysis_plan.md", "analysis_plan.md"),
    ("protocol.md", "protocol.md"),
    ("data_sources.md", "data_sources.md"),
    ("references_seed.md", "references_seed.md"),
    ("manuscript_draft_en_v0_1_2026-05-28.md", "manuscript_draft_en_v0_1_2026-05-28.md"),
    ("strobe_mr_checklist_2026-05-28.md", "strobe_mr_checklist_2026-05-28.md"),
    ("submission_statements_2026-05-28.md", "submission_statements_2026-05-28.md"),
    ("submission_formatting_notes_2026-05-28.md", "submission_formatting_notes_2026-05-28.md"),
]


INCLUDE_GLOBS = [
    ("scripts/*.py", "scripts"),
    ("tables/*.csv", "tables"),
    ("tables/*.md", "tables"),
    ("results/figures/*.svg", "results/figures"),
    ("results/coloc/*.csv", "results/coloc"),
    ("results/coloc/*.md", "results/coloc"),
    ("results/candidates/*", "results/candidates"),
    ("results/mediation/*summary.md", "results/mediation"),
    ("results/mediation/af_mediated_effects_fgf5_lpa.csv", "results/mediation"),
    ("results/mr/*summary*.md", "results/mr"),
    ("results/mr/fdr_significant_preliminary.csv", "results/mr"),
    ("results/mr/shared_candidate_preliminary.csv", "results/mr"),
    ("results/overlap/*", "results/overlap"),
    ("results/prioritization/*", "results/prioritization"),
    ("results/qc/*summary.md", "results/qc"),
    ("results/qc/candidate_instrument_qc.csv", "results/qc"),
    ("results/replication/*.md", "results/replication"),
    ("results/replication/*.csv", "results/replication"),
    ("results/replication/full_panel/*summary.md", "results/replication/full_panel"),
    ("results/replication/ukb_opengwas/*summary.md", "results/replication/ukb_opengwas"),
    ("results/reverse_mr/*summary.md", "results/reverse_mr"),
    ("results/reverse_mr/candidate_reverse_mr_results.csv", "results/reverse_mr"),
    ("results/sensitivity/*summary.md", "results/sensitivity"),
    ("results/sensitivity/candidate_sensitivity_summary.csv", "results/sensitivity"),
    ("results/text/*.md", "results/text"),
]


EXCLUDE_PARTS = {
    ".git",
    "__pycache__",
    "data",
    "raw",
    "public_release",
}

EXCLUDE_SUFFIXES = {
    ".gz",
    ".zip",
    ".tar",
    ".tgz",
    ".tbi",
    ".vcf",
    ".bcf",
    ".pyc",
}


def is_allowed(path: Path) -> bool:
    rel = path.relative_to(ROOT)
    if any(part in EXCLUDE_PARTS for part in rel.parts):
        return False
    if path.suffix.lower() in EXCLUDE_SUFFIXES:
        return False
    return True


def is_release_file_allowed(path: Path) -> bool:
    rel = path.relative_to(OUT)
    if any(part in {".git", "__pycache__"} for part in rel.parts):
        return False
    if path.suffix.lower() in EXCLUDE_SUFFIXES:
        return False
    return True


def copy_file(src: Path, dst_rel: Path, copied: list[tuple[str, int]]) -> None:
    if not src.exists() or not src.is_file():
        return
    if not is_allowed(src):
        return
    dst = OUT / dst_rel
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dst)
    copied.append((str(dst_rel).replace("\\", "/"), dst.stat().st_size))


def main() -> None:
    RELEASE_ROOT.mkdir(exist_ok=True)
    OUT.mkdir(parents=True, exist_ok=True)
    copied: list[tuple[str, int]] = []

    for src_rel, dst_rel in ROOT_FILES:
        copy_file(ROOT / src_rel, Path(dst_rel), copied)

    for pattern, dst_dir in INCLUDE_GLOBS:
        for src in sorted(ROOT.glob(pattern)):
            if src.is_file():
                copy_file(src, Path(dst_dir) / src.name, copied)

    manifest_path = OUT / "PUBLIC_RELEASE_MANIFEST.csv"
    with manifest_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["relative_path", "bytes"])
        for rel, size in sorted(set(copied)):
            writer.writerow([rel, size])
    copied.append(("PUBLIC_RELEASE_MANIFEST.csv", manifest_path.stat().st_size))

    if ZIP_PATH.exists():
        ZIP_PATH.unlink()
    with zipfile.ZipFile(ZIP_PATH, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        for file_path in sorted(OUT.rglob("*")):
            if file_path.is_file() and is_release_file_allowed(file_path):
                zf.write(file_path, file_path.relative_to(OUT.parent))

    print(f"Release folder: {OUT}")
    print(f"Release zip: {ZIP_PATH}")
    print(f"Files copied: {len(set(rel for rel, _ in copied))}")


if __name__ == "__main__":
    main()
