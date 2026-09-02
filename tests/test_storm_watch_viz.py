"""storm-watch situational map (storm_watch.viz): the alert + auto-domain figure."""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))

DATA = os.path.join(os.path.dirname(__file__), "data", "sample_tornado_alert.json")


def test_map_alert_file_writes_a_png(tmp_path):
    from atmospheric_data.storm_watch import viz
    from atmospheric_data.storm_watch.config import StormWatchConfig
    out = str(tmp_path / "map.png")
    p = viz.map_alert_file(DATA, cfg=StormWatchConfig(), out_path=out)
    assert p == out and os.path.exists(out)
    assert os.path.getsize(out) > 5000            # a real figure, not an empty file


def test_plot_alert_domain_uses_the_auto_domain(tmp_path):
    from atmospheric_data.storm_watch import viz, alerts as A, domain as D
    from atmospheric_data.storm_watch.config import StormWatchConfig
    cfg = StormWatchConfig()
    alert = next(a for a in A.load_alerts_file(DATA) if a.polygon)
    d = D.build_domain(alert, cfg)
    assert d["width_km"] > 0 and "center_lat" in d
    out = str(tmp_path / "m.png")
    viz.plot_alert_domain(alert, cfg, out, level="confirmed")
    assert os.path.exists(out)


def test_no_polygon_raises(tmp_path):
    from atmospheric_data.storm_watch import viz
    from atmospheric_data.storm_watch.alerts import Alert
    from atmospheric_data.storm_watch.config import StormWatchConfig
    import pytest
    with pytest.raises(ValueError):
        viz.plot_alert_domain(Alert(alert_id="x"), StormWatchConfig(), str(tmp_path / "n.png"))


if __name__ == "__main__":
    import tempfile
    d = tempfile.mkdtemp()

    class _P:
        def __truediv__(self, n): return os.path.join(d, n)
    test_map_alert_file_writes_a_png(_P()); print("ok map_alert_file")
    test_plot_alert_domain_uses_the_auto_domain(_P()); print("ok plot_alert_domain")
    print("ALL VIZ TESTS PASSED")
