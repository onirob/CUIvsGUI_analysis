from __future__ import annotations

import json
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PATTERNS = {
    "email": re.compile(r"(?i)\b[A-Z0-9._%+-]+@[A-Z0-9.-]+\.[A-Z]{2,}\b"),
    "24-character platform identifier": re.compile(r"(?i)\b[0-9a-f]{24}\b"),
    "UUID": re.compile(r"(?i)\b[0-9a-f]{8}-[0-9a-f]{4}-[1-5][0-9a-f]{3}-[89ab][0-9a-f]{3}-[0-9a-f]{12}\b"),
    "user-specific Windows path": re.compile(r"(?i)\b[A-Z]:(?:\\{1,2}|/)Users(?:\\{1,2}|/)"),
    "absolute Unix filesystem path": re.compile(r"(?<![A-Za-z0-9:/])/(?:home|Users|tmp|var|opt|mnt|workspace)/[^\s\"'<>]+"),
}
TEXT_SUFFIXES = {".md", ".txt", ".csv", ".py", ".r", ".tex", ".ipynb", ".cff", ""}
failures = []

for path in ROOT.rglob("*"):
    if not path.is_file() or path.suffix.lower() not in TEXT_SUFFIXES:
        continue
    relative = path.relative_to(ROOT)
    text = path.read_text(encoding="utf-8", errors="ignore")
    for label, pattern in PATTERNS.items():
        if pattern.search(text):
            failures.append(f"{relative}: {label}")
    if path.suffix.lower() == ".ipynb":
        notebook = json.loads(text)
        if not notebook.get("metadata", {}).get("public_release", {}).get("outputs_preserved"):
            failures.append(f"{relative}: missing public-release output-preservation marker")

if failures:
    print("PUBLIC RELEASE CHECK FAILED")
    for failure in failures:
        print("-", failure)
    sys.exit(1)

print("PUBLIC RELEASE CHECK PASSED")
