"""Local desktop notifications (storm-watch) — best-effort, no external service required."""
from __future__ import annotations

import os


def notify(title, message, enabled=True, logfile=None):
    """Emit a local notification.  Tries a native toast (plyer / win10toast) if installed;
    always logs to the console and, optionally, to ``logfile``.  Never raises."""
    line = "[notify] %s — %s" % (title, message)
    print(line, flush=True)
    if logfile:
        try:
            os.makedirs(os.path.dirname(os.path.abspath(logfile)) or ".", exist_ok=True)
            with open(logfile, "a", encoding="utf-8") as f:
                f.write(line + "\n")
        except Exception:
            pass
    if not enabled:
        return
    try:                                                  # optional native toast
        from plyer import notification as _n
        _n.notify(title=title, message=message, timeout=8)
        return
    except Exception:
        pass
    try:
        from win10toast import ToastNotifier
        ToastNotifier().show_toast(title, message, duration=8, threaded=True)
    except Exception:
        pass                                              # console log already emitted
