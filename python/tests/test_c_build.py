"""Tests for compiling the C backend (c_wrapper.build_c)."""

import platform
import subprocess
import sys
from pathlib import Path

import pytest

from heat1d.c_wrapper import C_DIR, _native_arch_prefix, build_c


def _hardware_arch():
    """Architecture of the machine, not of this interpreter."""
    if _under_rosetta():
        return "arm64"
    return platform.machine()


def _under_rosetta():
    """True when this interpreter is an x86_64 process on Apple Silicon."""
    if sys.platform != "darwin":
        return False
    r = subprocess.run(["sysctl", "-n", "sysctl.proc_translated"],
                       capture_output=True, text=True)
    return r.returncode == 0 and r.stdout.strip() == "1"


class TestNativeArchPrefix:
    """The build must target the hardware, not the interpreter.

    A Python under Rosetta reports x86_64, so an un-prefixed `make`
    compiles x86_64 objects that cannot link against the arm64 Homebrew
    libraries the Makefile points at.
    """

    def test_matches_translation_state(self):
        prefix = _native_arch_prefix()
        if _under_rosetta():
            assert prefix == ["arch", "-arm64"]
        else:
            assert prefix == []

    def test_empty_off_macos(self, monkeypatch):
        monkeypatch.setattr(sys, "platform", "linux")
        assert _native_arch_prefix() == []

    def test_empty_when_sysctl_missing(self, monkeypatch):
        """Intel Macs have no sysctl.proc_translated key."""
        monkeypatch.setattr(sys, "platform", "darwin")

        def _fail(*a, **kw):
            return subprocess.CompletedProcess(a, 1, stdout="", stderr="")

        monkeypatch.setattr(subprocess, "run", _fail)
        assert _native_arch_prefix() == []

    def test_survives_missing_sysctl_binary(self, monkeypatch):
        monkeypatch.setattr(sys, "platform", "darwin")

        def _raise(*a, **kw):
            raise OSError("no sysctl")

        monkeypatch.setattr(subprocess, "run", _raise)
        assert _native_arch_prefix() == []


class TestBuildC:

    def test_returns_existing_exe_without_building(self, monkeypatch, tmp_path):
        """The fast path must not shell out at all."""
        exe = tmp_path / "heat1d"
        exe.touch()

        def _boom(*a, **kw):
            raise AssertionError("build_c rebuilt an existing executable")

        monkeypatch.setattr(subprocess, "run", _boom)
        assert build_c(c_dir=tmp_path) == exe

    def test_failure_reports_both_attempts(self, tmp_path):
        """An empty dir has no Makefile, so both attempts fail."""
        with pytest.raises(RuntimeError) as exc:
            build_c(c_dir=tmp_path)
        msg = str(exc.value)
        assert "first attempt" in msg
        assert "after clean" in msg


@pytest.mark.slow
class TestBuildCIntegration:
    """Actually compiles; needs a toolchain plus fftw and libyaml."""

    def test_clean_build_produces_native_binaries(self):
        subprocess.run([*_native_arch_prefix(), "make", "-C", str(C_DIR),
                        "clean"], capture_output=True, text=True)
        exe = build_c()
        assert exe.exists()
        assert Path(C_DIR, "test_validate").exists()

        if sys.platform == "darwin":
            out = subprocess.run(["file", "-b", str(exe)],
                                 capture_output=True, text=True).stdout
            # Built for the hardware, which under Rosetta is NOT
            # platform.machine()
            assert _hardware_arch() in out, out

    def test_recovers_from_stale_objects(self):
        """Objects from a different arch are newer than their sources, so
        make reuses them and the link fails; build_c cleans and retries."""
        if not _under_rosetta():
            pytest.skip("needs Rosetta to plant x86_64 objects cheaply")
        subprocess.run([*_native_arch_prefix(), "make", "-C", str(C_DIR),
                        "clean"], capture_output=True, text=True)
        # Plant wrong-architecture objects (no arch prefix under Rosetta)
        subprocess.run(["make", "-C", str(C_DIR), "heat1d.o"],
                       capture_output=True, text=True)
        Path(C_DIR, "heat1d").unlink(missing_ok=True)

        exe = build_c()
        assert exe.exists()
