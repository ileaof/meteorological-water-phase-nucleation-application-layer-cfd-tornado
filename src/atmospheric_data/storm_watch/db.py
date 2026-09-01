"""Local SQLite persistence for storm-watch (dedup, cases, state transitions, queue).

The same alert is never processed twice (its id is stored); cases persist across monitor
restarts; every state transition, failure and retry is logged (task: "registre horário,
transições, falhas e tentativas").
"""
from __future__ import annotations

import json
import os
import sqlite3
import time

_SEVERITY_RANK = {"extreme": 4, "severe": 3, "moderate": 2, "minor": 1, "unknown": 0}


def _now():
    return time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())


class WatchDB:
    def __init__(self, path):
        os.makedirs(os.path.dirname(os.path.abspath(path)) or ".", exist_ok=True)
        self.path = path
        self.conn = sqlite3.connect(path)
        self.conn.row_factory = sqlite3.Row
        self._init()

    def _init(self):
        c = self.conn
        c.executescript("""
        CREATE TABLE IF NOT EXISTS alerts(
            alert_id TEXT PRIMARY KEY, event TEXT, severity TEXT, level TEXT,
            confirmation TEXT, first_seen TEXT, updated TEXT, status TEXT, payload TEXT);
        CREATE TABLE IF NOT EXISTS cases(
            case_id TEXT PRIMARY KEY, alert_id TEXT, name TEXT, state TEXT, severity TEXT,
            level TEXT, center_lat REAL, center_lon REAL, workdir TEXT, created TEXT,
            updated TEXT, retries INTEGER DEFAULT 0);
        CREATE TABLE IF NOT EXISTS transitions(
            id INTEGER PRIMARY KEY AUTOINCREMENT, case_id TEXT, from_state TEXT, to_state TEXT,
            time TEXT, note TEXT);
        CREATE TABLE IF NOT EXISTS queue(
            case_id TEXT PRIMARY KEY, priority INTEGER, enqueued TEXT);
        """)
        c.commit()

    # ---- alerts (dedup) ----
    def alert_seen(self, alert_id):
        return self.conn.execute("SELECT 1 FROM alerts WHERE alert_id=?", (alert_id,)).fetchone() is not None

    def record_alert(self, alert_id, event, severity, level, confirmation, payload, status="processed"):
        n = _now()
        row = self.conn.execute("SELECT first_seen FROM alerts WHERE alert_id=?", (alert_id,)).fetchone()
        first = row["first_seen"] if row else n
        self.conn.execute(
            "INSERT OR REPLACE INTO alerts VALUES(?,?,?,?,?,?,?,?,?)",
            (alert_id, event, severity, level, confirmation, first, n, status, json.dumps(payload)))
        self.conn.commit()

    def list_alerts(self, limit=50):
        return [dict(r) for r in self.conn.execute(
            "SELECT alert_id,event,severity,level,confirmation,updated,status FROM alerts "
            "ORDER BY updated DESC LIMIT ?", (limit,)).fetchall()]

    # ---- cases + state machine ----
    def upsert_case(self, case_id, alert_id, name, state, severity, level, lat, lon, workdir):
        n = _now()
        exists = self.conn.execute("SELECT state FROM cases WHERE case_id=?", (case_id,)).fetchone()
        if exists is None:
            self.conn.execute("INSERT INTO cases(case_id,alert_id,name,state,severity,level,"
                              "center_lat,center_lon,workdir,created,updated) VALUES(?,?,?,?,?,?,?,?,?,?,?)",
                              (case_id, alert_id, name, state, severity, level, lat, lon, workdir, n, n))
            self.conn.execute("INSERT INTO transitions(case_id,from_state,to_state,time,note) "
                              "VALUES(?,?,?,?,?)", (case_id, "", state, n, "created"))
        self.conn.commit()

    def set_state(self, case_id, state, note=""):
        cur = self.conn.execute("SELECT state FROM cases WHERE case_id=?", (case_id,)).fetchone()
        prev = cur["state"] if cur else ""
        n = _now()
        self.conn.execute("UPDATE cases SET state=?, updated=? WHERE case_id=?", (state, n, case_id))
        self.conn.execute("INSERT INTO transitions(case_id,from_state,to_state,time,note) "
                          "VALUES(?,?,?,?,?)", (case_id, prev, state, n, note))
        self.conn.commit()

    def bump_retry(self, case_id):
        self.conn.execute("UPDATE cases SET retries = retries + 1 WHERE case_id=?", (case_id,))
        self.conn.commit()

    def get_case(self, case_id):
        r = self.conn.execute("SELECT * FROM cases WHERE case_id=?", (case_id,)).fetchone()
        return dict(r) if r else None

    def list_cases(self, limit=100):
        return [dict(r) for r in self.conn.execute(
            "SELECT * FROM cases ORDER BY updated DESC LIMIT ?", (limit,)).fetchall()]

    def active_case_count(self):
        active = ("DETECTED", "DOWNLOADING", "READY_FOR_PREPROCESSING", "PREPROCESSING",
                  "READY_FOR_SIMULATION", "SIMULATING", "VALIDATING")
        q = "SELECT COUNT(*) n FROM cases WHERE state IN (%s)" % ",".join("?" * len(active))
        return self.conn.execute(q, active).fetchone()["n"]

    def transitions(self, case_id):
        return [dict(r) for r in self.conn.execute(
            "SELECT from_state,to_state,time,note FROM transitions WHERE case_id=? ORDER BY id",
            (case_id,)).fetchall()]

    # ---- priority queue (by severity) ----
    def enqueue(self, case_id, severity):
        self.conn.execute("INSERT OR REPLACE INTO queue VALUES(?,?,?)",
                          (case_id, _SEVERITY_RANK.get(str(severity).lower(), 0), _now()))
        self.conn.commit()

    def dequeue(self):
        r = self.conn.execute("SELECT case_id FROM queue ORDER BY priority DESC, enqueued ASC "
                              "LIMIT 1").fetchone()
        if r is None:
            return None
        self.conn.execute("DELETE FROM queue WHERE case_id=?", (r["case_id"],))
        self.conn.commit()
        return r["case_id"]

    def queue_size(self):
        return self.conn.execute("SELECT COUNT(*) n FROM queue").fetchone()["n"]

    def close(self):
        self.conn.close()
