"""
Download raw 16S FASTQ files for ISS Veggie lettuce rhizosphere reanalysis.

Fetches root-zone samples (Root*, Rhizo*, Pillow*) from OSD-267 (VEG-01A)
and OSD-269 (VEG-03A). Skips ITS, leaves, seeds.

Run in background:
    python download_veggie_data.py
"""
from __future__ import annotations

import json
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from urllib.request import urlopen, Request

BASE = Path(__file__).parent.resolve()
OSDR_BASE = "https://osdr.nasa.gov"

# Root-zone sample keywords (case-insensitive substring match)
ROOT_ZONE_KEYS = ["root", "rhizo", "pillow"]

# Include amplicon "raw" files only. OSD-267 files carry explicit 16S suffix;
# OSD-269 uses generic "Amplicon" in the filename (study is 16S).
def is_wanted(name: str, osd: int) -> bool:
    low = name.lower()
    if not low.endswith((".fastq.gz", ".fq.gz")):
        return False
    # root-zone only
    if not any(k in low for k in ROOT_ZONE_KEYS):
        return False
    # exclude explicit ITS
    if "_its_" in low or "_its." in low or "/its" in low:
        return False
    # prefer raw over filtered; OSD-267 raw has no "_raw" token -> check "_R1_16S" / "_R2_16S"
    if "_filtered" in low:
        return False
    return True


def list_files(osd: int) -> list[dict]:
    url = f"{OSDR_BASE}/osdr/data/osd/files/{osd}?format=json"
    with urlopen(url, timeout=60) as r:
        d = json.load(r)
    study = next(iter(d["studies"].values()))
    return study.get("study_files", []) or study.get("files", [])


def download_one(remote_url: str, dest: Path) -> tuple[str, int, str]:
    url = OSDR_BASE + remote_url if remote_url.startswith("/") else remote_url
    try:
        req = Request(url, headers={"User-Agent": "pepper-aa/1.0"})
        with urlopen(req, timeout=600) as r, open(dest, "wb") as fh:
            while True:
                chunk = r.read(1 << 20)
                if not chunk:
                    break
                fh.write(chunk)
        return (dest.name, dest.stat().st_size, "ok")
    except Exception as e:  # noqa: BLE001
        if dest.exists():
            dest.unlink()
        return (dest.name, 0, f"error: {e}")


def main() -> int:
    total_planned_mb = 0.0
    queued: list[tuple[int, dict, Path]] = []
    for osd in (267, 269):
        files = list_files(osd)
        out_dir = BASE / f"OSD-{osd}" / "raw"
        out_dir.mkdir(parents=True, exist_ok=True)
        meta_zip_url = None
        for f in files:
            name = f["file_name"]
            size = (f.get("file_size") or 0) / 1e6
            # grab metadata ISA zip for both
            if name.endswith("ISA.zip") and f"OSD-{osd}" in name:
                meta_dest = BASE / f"OSD-{osd}" / name
                if not meta_dest.exists():
                    queued.append((osd, f, meta_dest))
                continue
            if not is_wanted(name, osd):
                continue
            dest = out_dir / name
            if dest.exists() and dest.stat().st_size > 0:
                continue  # already downloaded
            queued.append((osd, f, dest))
            total_planned_mb += size
    print(f"[download] {len(queued)} files queued, ~{total_planned_mb:.1f} MB", flush=True)
    n_ok = 0; n_fail = 0; bytes_ok = 0
    with ThreadPoolExecutor(max_workers=4) as pool:
        futs = {pool.submit(download_one, item[1]["remote_url"], item[2]): item for item in queued}
        for i, fut in enumerate(as_completed(futs), 1):
            name, sz, st = fut.result()
            if st == "ok":
                n_ok += 1; bytes_ok += sz
                print(f"[download] {i}/{len(queued)} ok: {name} ({sz/1e6:.1f} MB)", flush=True)
            else:
                n_fail += 1
                print(f"[download] {i}/{len(queued)} FAIL: {name} :: {st}", flush=True)
    print(f"[download] done: ok={n_ok} fail={n_fail} total_ok={bytes_ok/1e6:.1f} MB",
          flush=True)
    return 0 if n_fail == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
