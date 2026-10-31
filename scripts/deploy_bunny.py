"""Sync a built static site to a bunny.net storage zone and purge the pull zone.

Uploads every file in the build directory, deletes remote files that are no
longer in the build, then purges the pull zone cache. Uses only the standard
library.

Usage:
    python scripts/deploy_bunny.py docs/_build/html

Environment variables:
    BUNNY_STORAGE_ZONE      Storage zone name.
    BUNNY_STORAGE_PASSWORD  Storage zone password (FTP & API Access > Password).
    BUNNY_STORAGE_ENDPOINT  Storage API hostname for the zone's primary region.
                            Default: storage.bunnycdn.com (Falkenstein). Others
                            include ny.storage.bunnycdn.com, la.storage.bunnycdn.com,
                            uk.storage.bunnycdn.com, sg.storage.bunnycdn.com.
    BUNNY_PULL_ZONE_ID      Pull zone id to purge (optional; skipped if unset).
    BUNNY_API_KEY           Account API key, needed for the purge.
    BUNNY_API_BASE          Override API URLs (for testing). Optional.
"""

from __future__ import annotations

import hashlib
import json
import mimetypes
import os
import pathlib
import sys
import urllib.error
import urllib.parse
import urllib.request
from concurrent.futures import ThreadPoolExecutor

RETRIES = 3
WORKERS = 8


def _request(method: str, url: str, key: str, data: bytes | None = None, headers=None):
    headers = {"AccessKey": key, **(headers or {})}
    last_error = None
    for _ in range(RETRIES):
        req = urllib.request.Request(url, data=data, method=method, headers=headers)
        try:
            with urllib.request.urlopen(req, timeout=60) as resp:
                return resp.status, resp.read()
        except urllib.error.HTTPError as e:
            if e.code == 404 and method in ("GET", "DELETE"):
                return 404, b""
            last_error = f"{method} {url} -> HTTP {e.code}: {e.read()[:300]!r}"
            if e.code < 500 and e.code != 429:
                break
        except urllib.error.URLError as e:
            last_error = f"{method} {url} -> {e.reason}"
    raise RuntimeError(last_error)


class StorageZone:
    def __init__(self, base: str, zone: str, password: str):
        self.root = f"{base.rstrip('/')}/{urllib.parse.quote(zone)}/"
        self.password = password

    def url(self, path: str) -> str:
        return self.root + urllib.parse.quote(path)

    def list_files(self, prefix: str = "") -> list[str]:
        """Returns every file path in the zone below 'prefix', recursively."""
        status, body = _request("GET", self.url(prefix), self.password)
        if status == 404:
            return []
        files = []
        for obj in json.loads(body):
            name = obj["ObjectName"]
            if obj["IsDirectory"]:
                files += self.list_files(f"{prefix}{name}/")
            else:
                files.append(prefix + name)
        return files

    def upload(self, path: str, content: bytes) -> None:
        content_type = mimetypes.guess_type(path)[0] or "application/octet-stream"
        _request(
            "PUT",
            self.url(path),
            self.password,
            data=content,
            headers={
                "Content-Type": content_type,
                "Checksum": hashlib.sha256(content).hexdigest().upper(),
            },
        )

    def delete(self, path: str) -> None:
        _request("DELETE", self.url(path), self.password)


def main(build_dir: str) -> None:
    build = pathlib.Path(build_dir)
    if not (build / "index.html").is_file():
        sys.exit(f"{build} does not look like a built site (no index.html).")

    env = os.environ
    api_base = env.get("BUNNY_API_BASE")
    endpoint = env.get("BUNNY_STORAGE_ENDPOINT") or "storage.bunnycdn.com"
    storage = StorageZone(
        api_base or f"https://{endpoint}",
        env["BUNNY_STORAGE_ZONE"],
        env["BUNNY_STORAGE_PASSWORD"],
    )

    local = {
        p.relative_to(build).as_posix(): p
        for p in sorted(build.rglob("*"))
        if p.is_file()
    }
    remote = set(storage.list_files())

    def upload(item):
        rel, path = item
        storage.upload(rel, path.read_bytes())

    with ThreadPoolExecutor(WORKERS) as pool:
        list(pool.map(upload, local.items()))
    print(f"Uploaded {len(local)} files.")

    stale = sorted(remote - set(local))
    with ThreadPoolExecutor(WORKERS) as pool:
        list(pool.map(storage.delete, stale))
    print(f"Deleted {len(stale)} stale files.")

    pull_zone = env.get("BUNNY_PULL_ZONE_ID")
    if pull_zone:
        api = api_base or "https://api.bunny.net"
        _request(
            "POST",
            f"{api.rstrip('/')}/pullzone/{urllib.parse.quote(pull_zone)}/purgeCache",
            env["BUNNY_API_KEY"],
            data=b"",
        )
        print(f"Purged pull zone {pull_zone}.")
    else:
        print("BUNNY_PULL_ZONE_ID not set; skipped cache purge.")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    main(sys.argv[1])
