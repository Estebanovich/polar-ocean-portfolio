"""Minimal editor for MITgcm namelist files: set one key inside one group (&GROUP ... &)."""
from __future__ import annotations

import re


def set_param(text: str, group: str, key: str, value: str) -> str:
    """Replace ``key=...`` in group ``group`` (comment lines ignored), or add ``key=value,``
    as the last line of the group. Raises ValueError if the group is missing or the key is
    set more than once in it."""
    lines = text.splitlines(keepends=True)
    start = next((n for n, ln in enumerate(lines)
                  if re.match(rf"\s*&{group}\s*$", ln, re.I)), None)
    if start is None:
        raise ValueError(f"no namelist group &{group}")
    end = next(n for n in range(start + 1, len(lines))
               if re.match(r"\s*[&/]\s*$", lines[n]))
    hits = [n for n in range(start + 1, end)
            if not lines[n].lstrip().startswith("#")
            and re.match(rf"\s*{re.escape(key)}\s*=", lines[n], re.I)]
    if len(hits) > 1:
        raise ValueError(f"&{group} sets {key} {len(hits)} times")
    new = f" {key}={value},\n"
    if hits:
        lines[hits[0]] = new
    else:
        lines.insert(end, new)
    return "".join(lines)
