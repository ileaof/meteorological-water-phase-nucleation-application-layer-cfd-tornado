"""TLS trust resolution for the download sources.

A TLS-inspecting middlebox (consumer antivirus such as Norton/Avast, or a corporate
proxy such as Zscaler) terminates HTTPS locally and re-signs every server certificate
with a *private* root.  That root is installed in the operating system's trust store --
so ``curl`` and the browser are happy -- but it is **not** in :mod:`certifi`, which is
the bundle :mod:`requests` uses by default.  The symptom is a download path that works
from the shell and fails from Python on the identical URL::

    SSLError: CERTIFICATE_VERIFY_FAILED - unable to get local issuer certificate

:func:`ca_bundle` resolves this by handing :mod:`requests` a bundle that is
``certifi`` **plus the roots the operating system already trusts**.  Verification
still happens, and only against trust decisions the machine has already made --
certificate checking is never disabled anywhere in this module.

Resolution order (first hit wins):

1. ``MET_H2O_CA_BUNDLE`` -- explicit PEM path, for CI or an unusual store.
2. a merged ``certifi`` + OS-root bundle, written once into the cache directory
   (plus any PEMs in ``MET_H2O_EXTRA_CA_DIR``).
3. ``certifi``'s own bundle, when the OS store yields nothing extra.
"""
from __future__ import annotations

import os
import ssl
import tempfile
import time

#: Windows system stores worth harvesting roots from.
_WINDOWS_STORES = ("ROOT", "CA")

#: Rebuild the merged bundle when the cached copy is older than this.
_MAX_AGE_S = 7 * 24 * 3600

_MEMO = {}


def _certifi_path():
    try:
        import certifi
    except ImportError:
        return None
    return certifi.where()


def _os_root_pems():
    """PEM text for every root the OS trusts, or ``[]`` where unavailable.

    ``ssl.enum_certificates`` is Windows-only; on Linux/macOS OpenSSL already reads
    the system store, so there is nothing to merge and ``[]`` is the right answer.
    """
    if not hasattr(ssl, "enum_certificates"):
        return []
    pems = []
    for store in _WINDOWS_STORES:
        try:
            entries = ssl.enum_certificates(store)
        except Exception:
            continue
        for der, _encoding, _trust in entries:
            try:
                pems.append(ssl.DER_cert_to_PEM_cert(der))
            except Exception:
                continue  # a malformed store entry must not break the whole bundle
    return pems


def _extra_dir_pems():
    """PEM text from every ``*.pem``/``*.crt`` in ``MET_H2O_EXTRA_CA_DIR``."""
    directory = os.environ.get("MET_H2O_EXTRA_CA_DIR")
    if not directory or not os.path.isdir(directory):
        return []
    pems = []
    for name in sorted(os.listdir(directory)):
        if not name.lower().endswith((".pem", ".crt", ".cer")):
            continue
        try:
            with open(os.path.join(directory, name), "r", encoding="utf-8") as f:
                text = f.read()
        except Exception:
            continue
        if "BEGIN CERTIFICATE" in text:
            pems.append(text)
    return pems


def _default_cache_dir():
    return os.path.join(tempfile.gettempdir(), "met_h2o_tls")


def ca_bundle(refresh=False, cache_dir=None):
    """Return a CA bundle path suitable for ``requests.get(..., verify=...)``.

    Falls back to :mod:`certifi`'s bundle (or ``True``, meaning library default) when
    the OS store adds nothing.  The merged file is cached, so the store is harvested
    once rather than per request.
    """
    explicit = os.environ.get("MET_H2O_CA_BUNDLE")
    if explicit and os.path.exists(explicit):
        return explicit

    key = (cache_dir or _default_cache_dir(), bool(refresh))
    if not refresh and key in _MEMO:
        return _MEMO[key]

    certifi_path = _certifi_path()
    directory = cache_dir or _default_cache_dir()
    merged = os.path.join(directory, "ca_bundle_merged.pem")

    if not refresh and os.path.exists(merged):
        if time.time() - os.path.getmtime(merged) < _MAX_AGE_S:
            _MEMO[key] = merged
            return merged

    extras = _os_root_pems() + _extra_dir_pems()
    if not extras:
        result = certifi_path or True
        _MEMO[key] = result
        return result

    chunks = []
    if certifi_path and os.path.exists(certifi_path):
        try:
            with open(certifi_path, "r", encoding="utf-8") as f:
                chunks.append(f.read())
        except Exception:
            pass
    chunks.extend(extras)

    try:
        os.makedirs(directory, exist_ok=True)
        # Write-then-rename so a concurrent reader never sees a half-written bundle.
        tmp = merged + ".%d.tmp" % os.getpid()
        with open(tmp, "w", encoding="utf-8") as f:
            f.write("\n".join(c.strip() + "\n" for c in chunks if c and c.strip()))
        os.replace(tmp, merged)
    except Exception:
        result = certifi_path or True
        _MEMO[key] = result
        return result

    _MEMO[key] = merged
    return merged


def describe():
    """Diagnostic summary of what trust resolution found -- for troubleshooting."""
    bundle = ca_bundle()
    return {
        "bundle": bundle if isinstance(bundle, str) else "<library default>",
        "explicit_env": os.environ.get("MET_H2O_CA_BUNDLE") or None,
        "extra_ca_dir": os.environ.get("MET_H2O_EXTRA_CA_DIR") or None,
        "certifi": _certifi_path(),
        "os_roots_found": len(_os_root_pems()),
        "enum_certificates": hasattr(ssl, "enum_certificates"),
    }


__all__ = ["ca_bundle", "describe"]
