"""The storm-watch monitor loop (realtime / historical / replay).

Polls alerts, deduplicates (SQLite), filters by event/severity, classifies the trigger level,
and advances a case per new alert — respecting the resource limits and the safe automation
defaults.  Resilient: an API failure is logged and retried with exponential backoff; it never
kills the monitor.
"""
from __future__ import annotations

import os
import time

from .config import StormWatchConfig
from .db import WatchDB
from . import alerts as alerts_mod
from . import casemachine
from .notify import notify


class StormWatchMonitor:
    def __init__(self, sw=None, base_case_config=None, db_path=None, offline=False,
                 logger=print, max_n=24):
        self.sw = sw or StormWatchConfig()
        self.base = base_case_config
        self.offline = bool(offline)
        self.logger = logger
        self.max_n = max_n
        os.makedirs(self.sw.workdir, exist_ok=True)
        self.db = WatchDB(db_path or os.path.join(self.sw.workdir, "storm_watch.sqlite"))
        self._fail_streak = 0
        self.last_poll = None

    # ---- one polling cycle ----
    def poll_once(self, alerts=None):
        """Process a batch of alerts (fetched, or supplied for replay/historical).  Returns a
        summary dict.  Deduplicates, filters, classifies, and advances cases within limits."""
        self.last_poll = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
        if alerts is None:
            try:
                alerts = alerts_mod.fetch_active_alerts(self.sw, offline=self.offline)
                self._fail_streak = 0
            except Exception as e:                        # API down -> log, back off, keep running
                self._fail_streak += 1
                self.logger("[monitor] alert fetch failed (%s); backoff" % e)
                return {"fetched": 0, "error": str(e), "fail_streak": self._fail_streak}
        new, advanced, skipped = 0, 0, 0
        for al in alerts:
            if not al.alert_id or self.db.alert_seen(al.alert_id):
                skipped += 1
                continue
            level = alerts_mod.classify_level(al, self.sw)
            if level is None:                             # filtered out (event/severity)
                self.db.record_alert(al.alert_id, al.event, al.severity, "filtered",
                                     alerts_mod.confirmation_stage(al), _payload(al), status="filtered")
                skipped += 1
                continue
            new += 1
            self.db.record_alert(al.alert_id, al.event, al.severity, level,
                                 alerts_mod.confirmation_stage(al), _payload(al))
            notify("Storm-Watch: %s" % al.event, "%s [%s]" % (al.headline or al.affected_area, level),
                   enabled=self.sw.notifications.desktop,
                   logfile=os.path.join(self.sw.workdir, "notifications.log"))
            if self.db.active_case_count() >= self.sw.resource_limits.maximum_active_cases:
                self.logger("[monitor] resource limit: %d active cases; queueing %s"
                            % (self.db.active_case_count(), al.alert_id))
                cid = casemachine.case_id_for(al)
                self.db.upsert_case(cid, al.alert_id, al.event, "DETECTED", al.severity, level,
                                    0.0, 0.0, self.sw.workdir)
                self.db.enqueue(cid, al.severity)
                continue
            try:
                casemachine.advance_case(al, level, self.sw, self.db, base=self.base,
                                         offline=self.offline, logger=self.logger, max_n=self.max_n)
                advanced += 1
            except Exception as e:                        # never let one case kill the loop
                self.logger("[monitor] case error for %s: %s" % (al.alert_id, e))
        return {"fetched": len(alerts), "new": new, "advanced": advanced, "skipped": skipped,
                "active_cases": self.db.active_case_count(), "queued": self.db.queue_size()}

    # ---- replay / historical ----
    def replay(self, path):
        """Feed a saved alert file through the pipeline (no network) — full-automation test."""
        self.logger("[monitor] REPLAY %s" % path)
        return self.poll_once(alerts_mod.load_alerts_file(path))

    # ---- continuous loop ----
    def run(self, max_iterations=None):
        """Realtime loop.  ``max_iterations`` bounds it for tests; otherwise runs until stopped.
        Exponential backoff (capped) after consecutive fetch failures."""
        it = 0
        while max_iterations is None or it < max_iterations:
            summary = self.poll_once()
            self.logger("[monitor] poll @%s: %s" % (self.last_poll, summary))
            it += 1
            if max_iterations is not None and it >= max_iterations:
                break
            delay = self.sw.alert_poll_seconds
            if self._fail_streak:
                delay = min(delay * (2 ** self._fail_streak), 3600)
            time.sleep(delay)
        return it

    def status(self):
        return {"enabled": self.sw.enabled, "last_poll": self.last_poll,
                "active_cases": self.db.active_case_count(), "queued": self.db.queue_size(),
                "fail_streak": self._fail_streak, "workdir": self.sw.workdir,
                "auto": {"download": self.sw.actions.auto_download,
                         "preprocess": self.sw.actions.auto_preprocess,
                         "simulate": self.sw.actions.auto_simulate}}

    def close(self):
        self.db.close()


def _payload(al):
    return {"event": al.event, "severity": al.severity, "certainty": al.certainty,
            "urgency": al.urgency, "effective": al.effective_time, "onset": al.onset_time,
            "expires": al.expiration_time, "headline": al.headline, "area": al.affected_area,
            "office": al.issuing_office, "has_polygon": bool(al.polygon)}
