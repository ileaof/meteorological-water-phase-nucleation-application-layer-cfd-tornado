"""Local download cache (task requirement 9: avoid repeated downloads).

A tiny, dependency-free content store keyed by ``(source, key)``.  Sources check
:meth:`Cache.path` and skip the network when the file already exists; ``--offline`` refuses any
missing file rather than reaching out.
"""
from __future__ import annotations

import json
import os
import time


class Cache:
    def __init__(self, directory="data/cache", offline=False):
        self.directory = os.path.abspath(directory)
        self.offline = bool(offline)
        os.makedirs(self.directory, exist_ok=True)
        self._index_path = os.path.join(self.directory, "index.json")
        self.index = {}
        if os.path.exists(self._index_path):
            try:
                with open(self._index_path, "r", encoding="utf-8") as f:
                    self.index = json.load(f)
            except Exception:
                self.index = {}

    def path(self, source, key, ext=""):
        """Deterministic cache path for a ``(source, key)`` pair."""
        safe = "".join(c if c.isalnum() or c in "._-" else "_" for c in "%s__%s" % (source, key))
        sub = os.path.join(self.directory, source)
        os.makedirs(sub, exist_ok=True)
        return os.path.join(sub, safe + ext)

    def has(self, source, key, ext=""):
        return os.path.exists(self.path(source, key, ext))

    def record(self, source, key, path, meta=None):
        """Register a fetched file in the index (source, time, metadata)."""
        self.index[self.path(source, key)] = {"source": source, "key": key,
                                               "path": path, "fetched_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
                                               "meta": meta or {}}
        with open(self._index_path, "w", encoding="utf-8") as f:
            json.dump(self.index, f, indent=2)

    def require_offline_ok(self, source, key, ext=""):
        """In offline mode, raise if the file is not already cached/present."""
        if self.offline and not self.has(source, key, ext):
            raise FileNotFoundError(
                "offline mode: %s data for %r not found in cache %s -- place the file there "
                "or run without --offline" % (source, key, self.directory))
