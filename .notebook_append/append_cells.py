#!/usr/bin/env python3
"""Append new analysis cells to Uchannel_notebook.ipynb IN PLACE (idempotent).
Usage: python3 append_cells.py <notebook.ipynb> <new_cells.json>
Makes a timestamped .bak first; skips if the marker cell is already present."""
import json, sys, shutil, os

nb_path, cells_path = sys.argv[1], sys.argv[2]
MARKER = "New analyses (2026-07 session): polarity synthesis"

nb = json.load(open(nb_path))
new_cells = json.load(open(cells_path))

# normalize source to list-of-lines (nbformat-friendly) and ensure required keys
def norm(cell):
    src = cell.get("source", "")
    if isinstance(src, str):
        cell["source"] = src.splitlines(keepends=True)
    cell.setdefault("metadata", {})
    if cell["cell_type"] == "code":
        cell.setdefault("execution_count", None)
        cell.setdefault("outputs", [])
    return cell
new_cells = [norm(c) for c in new_cells]

already = any(MARKER in "".join(
    c.get("source","") if isinstance(c.get("source",""),str) else "".join(c.get("source",[])))
    for c in nb.get("cells", []))
if already:
    print("SKIP: marker already present; notebook not modified.")
    sys.exit(0)

bak = nb_path + ".bak_before_append"
if not os.path.exists(bak):
    shutil.copy2(nb_path, bak)
    print("backup ->", os.path.basename(bak))

nb["cells"].extend(new_cells)
tmp = nb_path + ".tmp"
with open(tmp, "w") as f:
    json.dump(nb, f, indent=1, ensure_ascii=False)
os.replace(tmp, nb_path)
print(f"appended {len(new_cells)} cells; notebook now has {len(nb['cells'])} cells.")
